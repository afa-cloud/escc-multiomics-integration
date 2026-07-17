#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
full_sha <- "--full-sha" %in% args

project_root <- normalizePath(getwd(), mustWork = TRUE)
data_root <- Sys.getenv("RESEARCH_DATA_ROOT", unset = path.expand("~/ResearchDataHub"))
data_root <- normalizePath(data_root, mustWork = TRUE)

inventory_path <- file.path(project_root, "results", "escc_local_datahub_inventory.tsv")
catalog_path <- file.path(data_root, "CATALOG.tsv")
output_path <- file.path(project_root, "results", "escc_local_datahub_inventory_validation.tsv")
summary_output_path <- file.path(
  project_root, "results", "escc_local_datahub_inventory_validation_summary.tsv"
)

required_inventory_columns <- c(
  "dataset_key", "status", "independence_group", "local_path"
)
required_catalog_columns <- c(
  "dataset_key", "status", "local_path", "manifest_path"
)

inventory <- fread(inventory_path, sep = "\t", colClasses = "character", na.strings = NULL)
catalog <- fread(catalog_path, sep = "\t", colClasses = "character", na.strings = NULL)

missing_inventory_columns <- setdiff(required_inventory_columns, names(inventory))
missing_catalog_columns <- setdiff(required_catalog_columns, names(catalog))
if (length(missing_inventory_columns) > 0L) {
  stop("项目清单缺少字段：", paste(missing_inventory_columns, collapse = ", "))
}
if (length(missing_catalog_columns) > 0L) {
  stop("统一仓库目录缺少字段：", paste(missing_catalog_columns, collapse = ", "))
}
if (anyDuplicated(inventory$dataset_key)) {
  stop("项目清单存在重复 dataset_key。")
}
if (anyDuplicated(catalog$dataset_key)) {
  stop("统一仓库 CATALOG.tsv 存在重复 dataset_key。")
}

catalog_subset <- catalog[, .(
  dataset_key,
  catalog_status = status,
  catalog_local_path = local_path,
  manifest_path
)]
validation <- merge(
  inventory,
  catalog_subset,
  by = "dataset_key",
  all.x = TRUE,
  sort = FALSE
)

validate_manifest <- function(dataset_path, manifest_path, run_sha = FALSE) {
  result <- list(
    manifest_exists = FALSE,
    manifest_file_count = 0L,
    manifest_total_bytes = 0,
    files_present = FALSE,
    sizes_match = FALSE,
    sha_checked = run_sha,
    sha_match = if (run_sha) FALSE else NA,
    problem_files = ""
  )
  if (!nzchar(manifest_path) || !file.exists(manifest_path)) {
    return(result)
  }

  result$manifest_exists <- TRUE
  manifest <- fread(manifest_path, sep = "\t", colClasses = "character", na.strings = NULL)
  size_column <- if ("size_bytes" %in% names(manifest)) {
    "size_bytes"
  } else if ("file_size" %in% names(manifest)) {
    "file_size"
  } else {
    NA_character_
  }
  needed <- c("relative_path", "sha256")
  if (!all(needed %in% names(manifest)) || is.na(size_column)) {
    result$problem_files <- "manifest_missing_required_columns"
    return(result)
  }

  file_paths <- file.path(dataset_path, manifest$relative_path)
  present <- file.exists(file_paths)
  observed_sizes <- rep(NA_real_, length(file_paths))
  observed_sizes[present] <- file.info(file_paths[present])$size
  expected_sizes <- suppressWarnings(as.numeric(manifest[[size_column]]))
  size_ok <- present & !is.na(expected_sizes) & observed_sizes == expected_sizes

  result$manifest_file_count <- nrow(manifest)
  result$manifest_total_bytes <- sum(expected_sizes, na.rm = TRUE)
  result$files_present <- all(present)
  result$sizes_match <- all(size_ok)

  sha_ok <- rep(NA, length(file_paths))
  if (run_sha && all(present)) {
    observed_sha <- vapply(
      file_paths,
      digest,
      character(1),
      algo = "sha256",
      file = TRUE,
      serialize = FALSE
    )
    sha_ok <- tolower(observed_sha) == tolower(manifest$sha256)
    result$sha_match <- all(sha_ok)
  }

  problem <- manifest$relative_path[!present | !size_ok | (run_sha & !is.na(sha_ok) & !sha_ok)]
  result$problem_files <- paste(problem, collapse = ";")
  result
}

validation[, catalog_match := !is.na(catalog_status)]
validation[, status_match := catalog_match & status == catalog_status]
validation[, path_match := catalog_match & normalizePath(local_path, mustWork = FALSE) ==
  normalizePath(catalog_local_path, mustWork = FALSE)]
validation[, path_exists := dir.exists(local_path)]
validation[, documentation_exists := file.exists(file.path(local_path, "DATASET.md")) |
  file.exists(file.path(local_path, "QUARANTINE.md"))]

manifest_results <- lapply(seq_len(nrow(validation)), function(i) {
  validate_manifest(
    dataset_path = validation$local_path[i],
    manifest_path = validation$manifest_path[i],
    run_sha = full_sha
  )
})
manifest_table <- rbindlist(manifest_results, fill = TRUE)
validation <- cbind(validation, manifest_table)

validation[, independence_group_record_count := .N, by = independence_group]
validation[, hard_validation_pass :=
  catalog_match & status_match & path_match & path_exists & documentation_exists &
    manifest_exists & files_present & sizes_match & (!full_sha | sha_match)]
validation[, analysis_eligible := hard_validation_pass & status == "verified"]
validation[, validation_status := fifelse(
  !hard_validation_pass,
  "failed",
  fifelse(status == "verified", "passed", "blocked_by_catalog_status")
)]

output_columns <- c(
  "dataset_key", "status", "independence_group", "independence_group_record_count",
  "catalog_match", "status_match", "path_match", "path_exists", "documentation_exists",
  "manifest_exists", "manifest_file_count", "files_present", "sizes_match",
  "manifest_total_bytes",
  "sha_checked", "sha_match", "problem_files", "hard_validation_pass",
  "analysis_eligible", "validation_status", "local_path"
)
fwrite(validation[, ..output_columns], output_path, sep = "\t", na = "NA", quote = FALSE)

validation_summary <- data.table(
  metric = c(
    "catalog_records", "inventory_records", "independence_groups",
    "verified_analysis_eligible_records", "catalog_status_blocked_records",
    "hard_validation_failures", "manifest_files", "manifest_total_bytes",
    "manifest_total_gib"
  ),
  value = c(
    nrow(catalog), nrow(validation), uniqueN(validation$independence_group),
    sum(validation$analysis_eligible),
    sum(validation$validation_status == "blocked_by_catalog_status"),
    sum(!validation$hard_validation_pass),
    sum(validation$manifest_file_count),
    sum(validation$manifest_total_bytes),
    sum(validation$manifest_total_bytes) / 1024^3
  ),
  full_sha_checked = full_sha,
  generated_date = as.character(Sys.Date())
)
fwrite(validation_summary, summary_output_path, sep = "\t", na = "NA", quote = FALSE)

failed <- validation[hard_validation_pass == FALSE]
message("数据记录数：", nrow(validation))
message("独立性分组数：", uniqueN(validation$independence_group))
message("manifest 文件数：", sum(validation$manifest_file_count))
message("manifest 总字节：", format(sum(validation$manifest_total_bytes), scientific = FALSE))
message("可进入分析的 verified 记录数：", sum(validation$analysis_eligible))
message("被目录状态阻断的记录数：", sum(validation$validation_status == "blocked_by_catalog_status"))
message("硬校验失败记录数：", nrow(failed))
message("输出：", output_path)
message("汇总输出：", summary_output_path)

if (nrow(failed) > 0L) {
  stop("本机数据资产校验失败；请查看输出表中的 problem_files 和状态列。")
}
