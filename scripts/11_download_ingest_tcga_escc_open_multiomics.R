#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(fs)
  library(httr2)
  library(jsonlite)
})

# TCGA-ESCC GDC DR45 开放多组学下载与入库。
#
# 为避免误触发数 GB 下载，必须显式选择运行模式：
#   Rscript scripts/11_download_ingest_tcga_escc_open_multiomics.R --metadata-only
#   Rscript scripts/11_download_ingest_tcga_escc_open_multiomics.R --download
#
# --metadata-only 仅锁定 DR45、查询七层文件、展平实体并断言 95/94/76 交集；
# --download 在上述基础上使用 .part 断点续传，逐文件核对 GDC MD5
# 和本地 SHA256，通过全量复制校验后才登记 CATALOG.tsv 与 data/datasets.tsv。

args <- commandArgs(trailingOnly = TRUE)
allowed_args <- c("--metadata-only", "--download")
if (length(args) != 1L || !args %in% allowed_args) {
  stop(
    "必须且只能指定 --metadata-only 或 --download；",
    "不允许默认触发大规模下载。"
  )
}
download_enabled <- identical(args, "--download")
download_workers <- suppressWarnings(as.integer(Sys.getenv("GDC_DOWNLOAD_WORKERS", "6")))
if (is.na(download_workers) || download_workers < 1L || download_workers > 12L) {
  stop("GDC_DOWNLOAD_WORKERS 必须是 1–12 的整数")
}

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required_project_files <- c(
  file.path(project_root, "PROJECT_INDEX.md"),
  file.path(project_root, "results", "tcga_escc_patient_whitelist.tsv"),
  file.path(project_root, "data", "datasets.tsv")
)
if (!all(file_exists(required_project_files))) {
  stop("必须从项目根目录运行；缺失：", paste(required_project_files[!file_exists(required_project_files)], collapse = ";"))
}

data_root <- Sys.getenv("RESEARCH_DATA_ROOT", path.expand("~/ResearchDataHub"))
data_root <- normalizePath(data_root, winslash = "/", mustWork = TRUE)
catalog_path <- file.path(data_root, "CATALOG.tsv")
if (!file_exists(catalog_path)) stop("ResearchDataHub/CATALOG.tsv 不可读：", catalog_path)

retrieval_date <- "2026-07-11"
snapshot_version <- "gdc_dr45_retrieved_20260711"
dataset_key <- "GDC_TCGA_ESCA_open_multiomics_gdc_dr45_retrieved_20260711"
if (format(Sys.Date(), "%Y-%m-%d") != retrieval_date) {
  stop("快照日期与实际启动日期不一致；请先显式升级版本目录和 retrieval_date")
}
stage_root <- file.path(
  data_root, "_incoming", "GDC", "TCGA-ESCA", snapshot_version
)
canonical_root <- file.path(
  data_root, "datasets", "public", "GDC", "TCGA-ESCA", snapshot_version
)
source_root <- file.path(stage_root, "00_source")
metadata_dir <- file.path(stage_root, "10_metadata")
reusable_dir <- file.path(stage_root, "20_reusable")
manifest_dir <- file.path(stage_root, "90_manifests")
dir_create(c(source_root, metadata_dir, reusable_dir, manifest_dir), recurse = TRUE)
stage_lock <- paste0(stage_root, ".run.lock")
if (!dir.create(stage_lock, showWarnings = FALSE)) {
  stop("同一 TCGA DR45 暂存目录已有另一个脚本实例运行：", stage_lock)
}
on.exit(unlink(stage_lock, recursive = TRUE, force = TRUE), add = TRUE)

whitelist <- fread(
  required_project_files[2L],
  colClasses = "character",
  na.strings = NULL
)
if (!"patient" %in% names(whitelist)) stop("ESCC 白名单缺少 patient 字段")
whitelist_patients <- sort(unique(whitelist$patient[nzchar(whitelist$patient)]))
if (length(whitelist_patients) != 96L || anyDuplicated(whitelist_patients)) {
  stop("DR45 ESCC 白名单必须是 96 个唯一患者")
}
if (!all(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", whitelist_patients))) {
  stop("白名单含非法 TCGA 患者条形码")
}

status_url <- "https://api.gdc.cancer.gov/status"
files_url <- "https://api.gdc.cancer.gov/files"
data_url <- "https://api.gdc.cancer.gov/data"
user_agent <- "ESCC-multiomics-project/1.0"

file_fields <- c(
  "file_id", "file_name", "file_size", "md5sum", "access", "state",
  "created_datetime", "updated_datetime", "data_category", "data_type",
  "data_format", "experimental_strategy", "platform",
  "analysis.analysis_id", "analysis.submitter_id", "analysis.workflow_type",
  "analysis.workflow_version", "analysis.state", "analysis.created_datetime",
  "analysis.updated_datetime", "associated_entities.entity_id",
  "associated_entities.entity_submitter_id", "associated_entities.entity_type",
  "associated_entities.case_id", "cases.case_id", "cases.submitter_id",
  "cases.project.project_id", "cases.samples.sample_id",
  "cases.samples.submitter_id", "cases.samples.sample_type",
  "cases.samples.tissue_type", "cases.samples.tumor_descriptor",
  "cases.samples.portions.portion_id", "cases.samples.portions.submitter_id",
  "cases.samples.portions.analytes.analyte_id",
  "cases.samples.portions.analytes.submitter_id",
  "cases.samples.portions.analytes.analyte_type",
  "cases.samples.portions.analytes.aliquots.aliquot_id",
  "cases.samples.portions.analytes.aliquots.submitter_id",
  "cases.samples.portions.analytes.aliquots.analyte_type"
)
expand_fields <- c("analysis", "associated_entities", "cases.samples.portions.analytes.aliquots")

layers <- data.table(
  layer_id = c(
    "rna_star", "mirna_bcgsc", "methylation_hm450", "masked_maf",
    "cnv_gene_ascat2", "cnv_segment_ascat2", "rppa"
  ),
  data_category = c(
    "Transcriptome Profiling", "Transcriptome Profiling", "DNA Methylation",
    "Simple Nucleotide Variation", "Copy Number Variation",
    "Copy Number Variation", "Proteome Profiling"
  ),
  data_type = c(
    "Gene Expression Quantification", "miRNA Expression Quantification",
    "Methylation Beta Value", "Masked Somatic Mutation",
    "Gene Level Copy Number", "Allele-specific Copy Number Segment",
    "Protein Expression Quantification"
  ),
  workflow_type = c(
    "STAR - Counts", "BCGSC miRNA Profiling",
    "SeSAMe Methylation Beta Estimation",
    "Aliquot Ensemble Somatic Variant Merging and Masking",
    "ASCAT2", "ASCAT2", NA_character_
  ),
  platform = c(
    "Illumina", "Illumina", "Illumina Human Methylation 450", "Illumina",
    "Affymetrix SNP 6.0", "Affymetrix SNP 6.0", "RPPA"
  ),
  experimental_strategy = c(
    "RNA-Seq", "miRNA-Seq", "Methylation Array", "WXS",
    "Genotyping Array", "Genotyping Array", "Reverse Phase Protein Array"
  ),
  data_format = c("TSV", "TXT", "TXT", "MAF", "TSV", "TXT", "TSV"),
  expected_file_count = c(95L, 96L, 96L, 96L, 96L, 96L, 78L),
  expected_patient_count = c(95L, 95L, 96L, 96L, 96L, 96L, 78L),
  expected_association_count = c(95L, 96L, 96L, 192L, 192L, 192L, 78L),
  expected_association_type = c(
    "aliquot", "aliquot", "aliquot", "aliquot", "aliquot", "aliquot", "portion"
  )
)

filter_clause <- function(field, values) {
  list(op = "in", content = list(field = field, value = as.list(values)))
}

build_filter <- function(layer) {
  clauses <- list(
    filter_clause("cases.project.project_id", "TCGA-ESCA"),
    filter_clause("cases.submitter_id", whitelist_patients),
    filter_clause("cases.samples.sample_type", "Primary Tumor"),
    filter_clause("access", "open"),
    filter_clause("state", "released"),
    filter_clause("data_category", layer$data_category),
    filter_clause("data_type", layer$data_type),
    filter_clause("platform", layer$platform),
    filter_clause("experimental_strategy", layer$experimental_strategy),
    filter_clause("data_format", layer$data_format)
  )
  if (!is.na(layer$workflow_type)) {
    clauses <- append(clauses, list(filter_clause("analysis.workflow_type", layer$workflow_type)))
  }
  list(op = "and", content = clauses)
}

query_definitions <- lapply(seq_len(nrow(layers)), function(i) {
  layer <- layers[i]
  list(
    layer_id = layer$layer_id,
    endpoint = files_url,
    filters = build_filter(layer),
    fields = file_fields,
    expand = expand_fields,
    size = 2000L,
    sort = "file_id:asc",
    format = "JSON",
    expected_file_count = layer$expected_file_count,
    note = "GDC Data Release 由 /status 锁定，不得在 files filter 中加 data_release。"
  )
})
names(query_definitions) <- layers$layer_id

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

md5_file <- function(path) {
  digest(path, algo = "md5", file = TRUE, serialize = FALSE)
}

read_raw_file <- function(path) {
  readBin(path, what = "raw", n = as.numeric(file_info(path)$size))
}

write_raw_if_absent <- function(bytes, path) {
  if (file_exists(path)) return(invisible(FALSE))
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  writeBin(bytes, tmp)
  if (!file.rename(tmp, path)) stop("无法原子写入：", path)
  invisible(TRUE)
}

write_text_identical_or_stop <- function(text, path) {
  text <- paste0(text, ifelse(endsWith(text, "\n"), "", "\n"))
  if (file_exists(path)) {
    old <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    old <- paste0(old, ifelse(endsWith(old, "\n"), "", "\n"))
    if (!identical(old, text)) stop("已有查询定义与当前脚本不一致，拒绝覆盖：", path)
    return(invisible(FALSE))
  }
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  writeLines(text, tmp, useBytes = TRUE)
  if (!file.rename(tmp, path)) stop("无法原子写入：", path)
  invisible(TRUE)
}

atomic_fwrite <- function(x, path) {
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "")
  if (!file.rename(tmp, path)) stop("无法原子更新：", path)
}

scalar <- function(x, field = NULL) {
  if (!is.null(field)) {
    if (is.null(x) || !is.list(x)) return(NA_character_)
    x <- x[[field]]
  }
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  value <- x[[1L]]
  if (is.null(value) || length(value) == 0L || is.list(value)) return(NA_character_)
  as.character(value[[1L]])
}

as_object_list <- function(x) {
  if (is.null(x) || length(x) == 0L) return(list())
  if (is.list(x) && !is.null(names(x))) return(list(x))
  x
}

collapse_unique <- function(x) {
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))
  if (length(x) == 0L) "" else paste(x, collapse = ";")
}

status_response <- request(status_url) |>
  req_user_agent(user_agent) |>
  req_retry(max_tries = 5L) |>
  req_perform()
status_raw_live <- resp_body_raw(status_response)
status_live <- fromJSON(rawToChar(status_raw_live), simplifyVector = TRUE)
status_ok <- identical(tolower(as.character(status_live$status)), "ok") ||
  identical(tolower(as.character(status_live$status)), "healthy")
if (!status_ok ||
    as.integer(status_live$data_release_version$major) != 45L ||
    as.integer(status_live$data_release_version$minor) != 0L ||
    !identical(as.character(status_live$data_release_version$release_date), "2025-12-04")) {
  stop("GDC /status 不是预期的 OK + DR45.0 (2025-12-04)；拒绝继续。")
}

status_path <- file.path(metadata_dir, "gdc_status_dr45.json")
if (file_exists(status_path)) {
  frozen_status <- fromJSON(rawToChar(read_raw_file(status_path)), simplifyVector = TRUE)
  if (as.integer(frozen_status$data_release_version$major) != 45L ||
      as.integer(frozen_status$data_release_version$minor) != 0L ||
      !identical(as.character(frozen_status$data_release_version$release_date), "2025-12-04")) {
    stop("已有冻结 status 不属于 DR45，拒绝覆盖。")
  }
} else {
  write_raw_if_absent(status_raw_live, status_path)
}

query_definition_path <- file.path(metadata_dir, "gdc_open_multiomics_queries_dr45.json")
query_definition_json <- toJSON(
  list(
    dataset_key = dataset_key,
    retrieval_date = retrieval_date,
    data_release = list(major = 45L, minor = 0L, release_date = "2025-12-04"),
    whitelist_source = "results/tcga_escc_patient_whitelist.tsv",
    whitelist_patients = whitelist_patients,
    queries = query_definitions
  ),
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null"
)
write_text_identical_or_stop(query_definition_json, query_definition_path)

# data.table 表达式与函数参数同名时显式避免作用域歧义。
query_results <- setNames(vector("list", nrow(layers)), layers$layer_id)
for (layer_name in layers$layer_id) {
  definition <- query_definitions[[layer_name]]
  raw_path <- file.path(metadata_dir, paste0("gdc_files_", layer_name, "_dr45.json"))
  if (!file_exists(raw_path)) {
    response <- request(files_url) |>
      req_user_agent(user_agent) |>
      req_url_query(
        filters = toJSON(definition$filters, auto_unbox = TRUE),
        fields = paste(definition$fields, collapse = ","),
        expand = paste(definition$expand, collapse = ","),
        size = definition$size,
        sort = definition$sort,
        format = definition$format
      ) |>
      req_retry(max_tries = 5L) |>
      req_perform()
    write_raw_if_absent(resp_body_raw(response), raw_path)
  }
  parsed <- fromJSON(rawToChar(read_raw_file(raw_path)), simplifyVector = FALSE)
  observed_total <- as.integer(parsed$data$pagination$total)
  expected_total <- layers[layer_id == layer_name, expected_file_count]
  if (observed_total != expected_total || length(parsed$data$hits) != expected_total) {
    stop(layer_name, " 文件数不符合 DR45 预期：observed=", observed_total, ", expected=", expected_total)
  }
  query_results[[layer_name]] <- list(raw_path = raw_path, parsed = parsed, hits = parsed$data$hits)
}

flatten_hierarchy <- function(hit) {
  rows <- list()
  row_index <- 0L
  for (case in as_object_list(hit$cases)) {
    case_id <- scalar(case, "case_id")
    case_submitter_id <- scalar(case, "submitter_id")
    project_id <- scalar(case$project, "project_id")
    samples <- as_object_list(case$samples)
    if (length(samples) == 0L) samples <- list(list())
    for (sample in samples) {
      portions <- as_object_list(sample$portions)
      if (length(portions) == 0L) portions <- list(list())
      for (portion in portions) {
        analytes <- as_object_list(portion$analytes)
        if (length(analytes) == 0L) analytes <- list(list())
        for (analyte in analytes) {
          aliquots <- as_object_list(analyte$aliquots)
          if (length(aliquots) == 0L) aliquots <- list(list())
          for (aliquot in aliquots) {
            row_index <- row_index + 1L
            rows[[row_index]] <- data.table(
              case_id = case_id,
              case_submitter_id = case_submitter_id,
              project_id = project_id,
              sample_id = scalar(sample, "sample_id"),
              sample_submitter_id = scalar(sample, "submitter_id"),
              sample_type = scalar(sample, "sample_type"),
              tissue_type = scalar(sample, "tissue_type"),
              tumor_descriptor = scalar(sample, "tumor_descriptor"),
              portion_id = scalar(portion, "portion_id"),
              portion_submitter_id = scalar(portion, "submitter_id"),
              analyte_id = scalar(analyte, "analyte_id"),
              analyte_submitter_id = scalar(analyte, "submitter_id"),
              analyte_type = scalar(analyte, "analyte_type"),
              aliquot_id = scalar(aliquot, "aliquot_id"),
              aliquot_submitter_id = scalar(aliquot, "submitter_id"),
              aliquot_analyte_type = scalar(aliquot, "analyte_type")
            )
          }
        }
      }
    }
  }
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

match_entity_to_hierarchy <- function(entity, hierarchy) {
  entity_id <- scalar(entity, "entity_id")
  entity_submitter_id <- scalar(entity, "entity_submitter_id")
  entity_type <- tolower(scalar(entity, "entity_type"))
  entity_case_id <- scalar(entity, "case_id")
  id_column <- switch(
    entity_type,
    aliquot = "aliquot_id",
    analyte = "analyte_id",
    portion = "portion_id",
    sample = "sample_id",
    case = "case_id",
    NA_character_
  )
  submitter_column <- switch(
    entity_type,
    aliquot = "aliquot_submitter_id",
    analyte = "analyte_submitter_id",
    portion = "portion_submitter_id",
    sample = "sample_submitter_id",
    case = "case_submitter_id",
    NA_character_
  )
  matched <- hierarchy[0]
  if (!is.na(id_column) && !is.na(entity_id)) {
    matched <- hierarchy[get(id_column) == entity_id]
  }
  if (nrow(matched) == 0L && !is.na(submitter_column) && !is.na(entity_submitter_id)) {
    matched <- hierarchy[get(submitter_column) == entity_submitter_id]
  }
  if (nrow(matched) > 0L && !is.na(entity_case_id)) {
    restricted <- matched[case_id == entity_case_id]
    if (nrow(restricted) > 0L) matched <- restricted
  }
  if (nrow(matched) == 0L) {
    matched <- data.table(
      case_id = entity_case_id, case_submitter_id = NA_character_, project_id = NA_character_,
      sample_id = NA_character_, sample_submitter_id = NA_character_, sample_type = NA_character_,
      tissue_type = NA_character_, tumor_descriptor = NA_character_, portion_id = NA_character_,
      portion_submitter_id = NA_character_, analyte_id = NA_character_,
      analyte_submitter_id = NA_character_, analyte_type = NA_character_, aliquot_id = NA_character_,
      aliquot_submitter_id = NA_character_, aliquot_analyte_type = NA_character_
    )
    matched[, hierarchy_match := FALSE]
  } else {
    matched[, hierarchy_match := TRUE]
  }
  matched[, `:=`(
    associated_entity_id = entity_id,
    associated_entity_submitter_id = entity_submitter_id,
    associated_entity_type = entity_type,
    associated_entity_case_id = entity_case_id
  )]
  matched
}

file_rows <- list()
entity_rows <- list()
file_index <- 0L
entity_index <- 0L

for (layer_name in layers$layer_id) {
  hits <- query_results[[layer_name]]$hits
  for (hit in hits) {
    file_index <- file_index + 1L
    file_id <- scalar(hit, "file_id")
    file_name <- scalar(hit, "file_name")
    if (is.na(file_id) || is.na(file_name) || basename(file_name) != file_name) {
      stop(layer_name, " 含缺失或不安全的 file_id/file_name")
    }
    analysis <- hit$analysis
    hierarchy <- flatten_hierarchy(hit)
    case_submitters <- unique(hierarchy$case_submitter_id)
    case_submitters <- case_submitters[!is.na(case_submitters) & case_submitters %in% whitelist_patients]
    if (length(case_submitters) != 1L) {
      stop(layer_name, "/", file_id, " 不能唯一定位到 96 例白名单中的患者")
    }
    associated_entities <- as_object_list(hit$associated_entities)
    if (length(associated_entities) == 0L) stop(layer_name, "/", file_id, " 缺失 associated_entities")
    associated_ids <- vapply(associated_entities, scalar, character(1), field = "entity_id")
    if (anyNA(associated_ids) || anyDuplicated(associated_ids)) {
      stop(layer_name, "/", file_id, " associated_entities 缺 ID 或重复")
    }
    for (entity_number in seq_along(associated_entities)) {
      entity_index <- entity_index + 1L
      matched <- match_entity_to_hierarchy(associated_entities[[entity_number]], hierarchy)
      matched[, `:=`(
        layer_id = layer_name,
        file_id = file_id,
        file_name = file_name,
        associated_entity_order = entity_number,
        hierarchy_path_order = seq_len(.N)
      )]
      entity_rows[[entity_index]] <- matched
    }
    target_relative_path <- file.path("00_source", layer_name, file_id, file_name)
    file_rows[[file_index]] <- data.table(
      layer_id = layer_name,
      patient_id = case_submitters,
      file_id = file_id,
      file_name = file_name,
      file_size = as.numeric(scalar(hit, "file_size")),
      gdc_md5 = tolower(scalar(hit, "md5sum")),
      access = scalar(hit, "access"),
      state = scalar(hit, "state"),
      created_datetime = scalar(hit, "created_datetime"),
      updated_datetime = scalar(hit, "updated_datetime"),
      data_category = scalar(hit, "data_category"),
      data_type = scalar(hit, "data_type"),
      data_format = scalar(hit, "data_format"),
      experimental_strategy = scalar(hit, "experimental_strategy"),
      platform = scalar(hit, "platform"),
      analysis_id = scalar(analysis, "analysis_id"),
      analysis_submitter_id = scalar(analysis, "submitter_id"),
      workflow_type = scalar(analysis, "workflow_type"),
      workflow_version = scalar(analysis, "workflow_version"),
      analysis_state = scalar(analysis, "state"),
      analysis_created_datetime = scalar(analysis, "created_datetime"),
      analysis_updated_datetime = scalar(analysis, "updated_datetime"),
      associated_entity_count = length(associated_entities),
      target_relative_path = target_relative_path,
      source_url = paste0(data_url, "/", file_id),
      download_status = "not_downloaded",
      local_size_bytes = NA_real_,
      local_md5 = "",
      local_sha256 = ""
    )
  }
}

file_manifest <- rbindlist(file_rows, use.names = TRUE, fill = TRUE)
entity_map <- rbindlist(entity_rows, use.names = TRUE, fill = TRUE)
setcolorder(entity_map, c(
  "layer_id", "file_id", "file_name", "associated_entity_order",
  "associated_entity_id", "associated_entity_submitter_id", "associated_entity_type",
  "associated_entity_case_id", "hierarchy_match", "hierarchy_path_order",
  "case_id", "case_submitter_id", "project_id", "sample_id", "sample_submitter_id",
  "sample_type", "tissue_type", "tumor_descriptor", "portion_id", "portion_submitter_id",
  "analyte_id", "analyte_submitter_id", "analyte_type", "aliquot_id",
  "aliquot_submitter_id", "aliquot_analyte_type"
))

if (any(!entity_map$hierarchy_match)) {
  failed <- unique(entity_map[hierarchy_match == FALSE, paste(layer_id, file_id, associated_entity_id, sep = "/")])
  stop("associated_entities 无法完整映射到 case/sample/portion/analyte/aliquot：", paste(failed, collapse = ";"))
}

entity_aggregate <- entity_map[, .(
  associated_entity_ids = collapse_unique(associated_entity_id),
  associated_entity_submitter_ids = collapse_unique(associated_entity_submitter_id),
  associated_entity_types = collapse_unique(associated_entity_type),
  associated_entity_case_ids = collapse_unique(associated_entity_case_id),
  case_ids = collapse_unique(case_id),
  case_submitter_ids = collapse_unique(case_submitter_id),
  sample_ids = collapse_unique(sample_id),
  sample_submitter_ids = collapse_unique(sample_submitter_id),
  sample_types = collapse_unique(sample_type),
  portion_ids = collapse_unique(portion_id),
  portion_submitter_ids = collapse_unique(portion_submitter_id),
  analyte_ids = collapse_unique(analyte_id),
  analyte_submitter_ids = collapse_unique(analyte_submitter_id),
  analyte_types = collapse_unique(analyte_type),
  aliquot_ids = collapse_unique(aliquot_id),
  aliquot_submitter_ids = collapse_unique(aliquot_submitter_id)
), by = .(layer_id, file_id)]
file_manifest <- merge(file_manifest, entity_aggregate, by = c("layer_id", "file_id"), all.x = TRUE, sort = FALSE)
file_manifest[, layer_order__ := match(layer_id, layers$layer_id)]
setorder(file_manifest, layer_order__, patient_id, file_id)
file_manifest[, layer_order__ := NULL]
entity_map[, layer_order__ := match(layer_id, layers$layer_id)]
setorder(entity_map, layer_order__, file_id, associated_entity_order, hierarchy_path_order)
entity_map[, layer_order__ := NULL]

observed_files <- file_manifest[, .(
  file_count = uniqueN(file_id),
  patient_count = uniqueN(patient_id)
), by = layer_id]
observed_associations <- unique(entity_map[, .(
  layer_id, file_id, associated_entity_id, associated_entity_type, sample_type
)])[, .(
  association_count = .N,
  association_types = collapse_unique(associated_entity_type)
), by = layer_id]
layer_checks <- merge(layers, observed_files, by = "layer_id", all.x = TRUE, sort = FALSE)
layer_checks <- merge(layer_checks, observed_associations, by = "layer_id", all.x = TRUE, sort = FALSE)
if (any(layer_checks$file_count != layer_checks$expected_file_count) ||
    any(layer_checks$patient_count != layer_checks$expected_patient_count) ||
    any(layer_checks$association_count != layer_checks$expected_association_count) ||
    any(layer_checks$association_types != layer_checks$expected_association_type)) {
  stop("七层文件/患者/associated_entities 计数或类型不符合 DR45 冻结审计")
}

association_by_sample <- unique(entity_map[, .(
  layer_id, file_id, associated_entity_id, associated_entity_type, sample_type
)])[, .N, by = .(layer_id, sample_type)]
expected_sample_associations <- data.table(
  layer_id = c(
    "rna_star", "mirna_bcgsc", "methylation_hm450",
    rep(c("masked_maf", "cnv_gene_ascat2", "cnv_segment_ascat2"), each = 3L),
    "rppa"
  ),
  sample_type = c(
    "Primary Tumor", "Primary Tumor", "Primary Tumor",
    rep(c("Primary Tumor", "Blood Derived Normal", "Solid Tissue Normal"), 3L),
    "Primary Tumor"
  ),
  expected_N = c(95L, 96L, 96L, rep(c(96L, 88L, 8L), 3L), 78L)
)
association_check <- merge(
  expected_sample_associations,
  association_by_sample,
  by = c("layer_id", "sample_type"),
  all = TRUE
)
association_check[is.na(N), N := 0L]
association_check[is.na(expected_N), expected_N := 0L]
if (any(association_check$N != association_check$expected_N)) {
  stop("MAF/CNV 配对正常或其他层 associated_entities 样本类型计数不符合审计")
}
paired_input_layers <- c("masked_maf", "cnv_gene_ascat2", "cnv_segment_ascat2")
per_file_pairing_check <- unique(entity_map[
  layer_id %in% paired_input_layers,
  .(layer_id, file_id, associated_entity_id, sample_type)
])[, .(
  associated_entity_count = .N,
  primary_tumor_count = sum(sample_type == "Primary Tumor"),
  matched_normal_count = sum(sample_type %in% c(
    "Blood Derived Normal", "Solid Tissue Normal"
  ))
), by = .(layer_id, file_id)]
if (nrow(per_file_pairing_check) != 288L ||
    any(per_file_pairing_check$associated_entity_count != 2L) ||
    any(per_file_pairing_check$primary_tumor_count != 1L) ||
    any(per_file_pairing_check$matched_normal_count != 1L)) {
  stop("MAF/CNV 存在未满足每文件 1 tumor + 1 matched normal 的记录")
}
entity_map[, analysis_role := fifelse(
  layer_id %in% c("masked_maf", "cnv_gene_ascat2", "cnv_segment_ascat2") &
    sample_type %in% c("Blood Derived Normal", "Solid Tissue Normal"),
  "matched_normal_associated_input",
  fifelse(
    layer_id %in% c("masked_maf", "cnv_gene_ascat2", "cnv_segment_ascat2") & sample_type == "Primary Tumor",
    "tumor_associated_input",
    fifelse(sample_type == "Primary Tumor", "primary_tumor_measurement", "context_only_not_expression_normal")
  )
)]

availability <- data.table(patient_id = whitelist_patients)
for (layer_name in layers$layer_id) {
  counts <- file_manifest[layer_id == layer_name, .(file_count = uniqueN(file_id)), by = patient_id]
  count_name <- paste0(layer_name, "_file_count")
  flag_name <- paste0(layer_name, "_available")
  setnames(counts, "file_count", count_name)
  availability <- merge(availability, counts, by = "patient_id", all.x = TRUE, sort = FALSE)
  set(availability, which(is.na(availability[[count_name]])), count_name, 0L)
  availability[, (flag_name) := get(count_name) > 0L]
}
setorder(availability, patient_id)
availability[, driver_core := rna_star_available & masked_maf_available & cnv_gene_ascat2_available]
availability[, five_layer_core :=
  driver_core & mirna_bcgsc_available & methylation_hm450_available
]
availability[, protein_deep_subset := five_layer_core & rppa_available]

intersection_summary <- data.table(
  analysis_set = c("clinical_whitelist", "driver_core", "five_layer_core", "protein_deep_subset"),
  included_layers = c(
    "clinical", "rna_star;masked_maf;cnv_gene_ascat2",
    "rna_star;mirna_bcgsc;methylation_hm450;masked_maf;cnv_gene_ascat2",
    "rna_star;mirna_bcgsc;methylation_hm450;masked_maf;cnv_gene_ascat2;rppa"
  ),
  patient_count = c(
    nrow(availability), sum(availability$driver_core),
    sum(availability$five_layer_core), sum(availability$protein_deep_subset)
  ),
  expected_patient_count = c(96L, 95L, 94L, 76L),
  boundary = c(
    "各层按实际可用患者分析，不要求每例七层齐全",
    "ASCAT2 segment 是同一 CNV 语义的片段层，不另算独立组学",
    "miRNA 双 aliquot 患者不得把两个文件当两例",
    "RPPA 缺失只限制蛋白深描子集，不屏蔽其他层证据"
  )
)
if (any(intersection_summary$patient_count != intersection_summary$expected_patient_count)) {
  stop("DR45 患者交集不符合 96/95/94/76 审计")
}

file_manifest_path <- file.path(metadata_dir, "tcga_escc_open_multiomics_file_manifest.tsv")
entity_map_path <- file.path(metadata_dir, "tcga_escc_open_multiomics_entity_map.tsv")
layer_checks_path <- file.path(metadata_dir, "tcga_escc_open_multiomics_layer_checks.tsv")
association_checks_path <- file.path(metadata_dir, "tcga_escc_open_multiomics_association_checks.tsv")
availability_path <- file.path(reusable_dir, "tcga_escc_open_multiomics_patient_availability.tsv")
intersection_path <- file.path(reusable_dir, "tcga_escc_open_multiomics_intersections.tsv")
atomic_fwrite(file_manifest, file_manifest_path)
atomic_fwrite(entity_map, entity_map_path)
atomic_fwrite(layer_checks, layer_checks_path)
atomic_fwrite(association_check, association_checks_path)
atomic_fwrite(availability, availability_path)
atomic_fwrite(intersection_summary, intersection_path)

if (!download_enabled) {
  message(
    "已完成 DR45 元数据冻结和断言；未下载数据文件，未登记为 verified。\n",
    "进入完整下载需显式重跑：--download"
  )
  quit(save = "no", status = 0L)
}

curl_bin <- Sys.which("curl")
if (!nzchar(curl_bin)) stop("未找到 curl，无法执行可断点下载")

download_one <- function(row) {
  target <- file.path(stage_root, row$target_relative_path)
  part <- paste0(target, ".part")
  dir_create(dirname(target), recurse = TRUE)
  if (file_exists(target)) {
    observed_size <- as.numeric(file_info(target)$size)
    observed_md5 <- tolower(md5_file(target))
    if (observed_size != row$file_size || observed_md5 != row$gdc_md5) {
      stop("已有正式目标文件与 GDC 大小/MD5 不符，拒绝覆盖：", target)
    }
    return(list(
      status = "verified", size = observed_size,
      md5 = observed_md5, sha256 = sha256_file(target)
    ))
  }
  if (file_exists(part) && as.numeric(file_info(part)$size) > row$file_size) {
    stop(".part 大于 GDC 声明大小，请人工核查：", part)
  }
  if (file_exists(part) && as.numeric(file_info(part)$size) == row$file_size) {
    part_md5 <- tolower(md5_file(part))
    if (part_md5 != row$gdc_md5) {
      stop("完整大小的 .part 与 GDC MD5 不符，拒绝覆盖：", part)
    }
    part_sha256 <- sha256_file(part)
    if (!file.rename(part, target)) stop("无法原子提升已完整 .part：", target)
    return(list(
      status = "verified", size = row$file_size,
      md5 = part_md5, sha256 = part_sha256
    ))
  }
  status <- system2(
    curl_bin,
    c(
      "--fail", "--location", "--silent", "--show-error",
      "--retry", "5", "--retry-delay", "5",
      "--continue-at", "-", "--output", shQuote(part), shQuote(row$source_url)
    )
  )
  if (!identical(status, 0L)) stop("GDC 下载失败；.part 已保留可续传：", part)
  observed_size <- as.numeric(file_info(part)$size)
  observed_md5 <- tolower(md5_file(part))
  if (observed_size != row$file_size || observed_md5 != row$gdc_md5) {
    stop("GDC 文件大小或 MD5 校验失败；拒绝提升 .part：", part)
  }
  observed_sha256 <- sha256_file(part)
  if (!file.rename(part, target)) stop("无法原子提升 .part：", target)
  list(status = "verified", size = observed_size, md5 = observed_md5, sha256 = observed_sha256)
}

download_batches <- split(
  seq_len(nrow(file_manifest)),
  ceiling(seq_len(nrow(file_manifest)) / download_workers)
)
completed_count <- 0L
for (batch in download_batches) {
  safe_download <- function(i) {
    tryCatch(
      list(ok = TRUE, index = i, result = download_one(file_manifest[i])),
      error = function(e) list(ok = FALSE, index = i, error = conditionMessage(e))
    )
  }
  batch_results <- if (download_workers == 1L || length(batch) == 1L) {
    lapply(batch, safe_download)
  } else {
    parallel::mclapply(
      batch,
      safe_download,
      mc.cores = min(download_workers, length(batch)),
      mc.preschedule = FALSE
    )
  }

  for (item in batch_results) {
    if (isTRUE(item$ok)) {
      i <- item$index
      result <- item$result
      file_manifest[i, `:=`(
        download_status = result$status,
        local_size_bytes = result$size,
        local_md5 = result$md5,
        local_sha256 = result$sha256
      )]
    }
  }
  atomic_fwrite(file_manifest, file_manifest_path)
  failures <- Filter(function(item) !isTRUE(item$ok), batch_results)
  if (length(failures) > 0L) {
    stop(
      "GDC 并发批次存在失败；成功文件已留存，.part 可续传：",
      paste(vapply(
        failures,
        function(item) paste0("row=", item$index, " ", item$error),
        character(1)
      ), collapse = "; ")
    )
  }
  completed_count <- completed_count + length(batch)
  if (completed_count %% (download_workers * 5L) == 0L ||
      completed_count == nrow(file_manifest)) {
    message("已校验 ", completed_count, "/", nrow(file_manifest), " 个 GDC 文件")
  }
}

if (!all(file_manifest$download_status == "verified") ||
    any(file_manifest$local_size_bytes != file_manifest$file_size) ||
    any(tolower(file_manifest$local_md5) != tolower(file_manifest$gdc_md5)) ||
    any(!grepl("^[0-9a-f]{64}$", file_manifest$local_sha256))) {
  stop("全量文件校验未通过，不得入库")
}

provenance <- data.table(
  output_scope = c(
    "10_metadata/gdc_files_<layer>_dr45.json",
    "10_metadata/tcga_escc_open_multiomics_file_manifest.tsv",
    "10_metadata/tcga_escc_open_multiomics_entity_map.tsv",
    "20_reusable/tcga_escc_open_multiomics_patient_availability.tsv",
    "20_reusable/tcga_escc_open_multiomics_intersections.tsv",
    "00_source/<layer>/<file_id>/<file_name>"
  ),
  corresponding_source_file = c(
    "GDC /files API", "10_metadata/gdc_files_<layer>_dr45.json",
    "10_metadata/gdc_files_<layer>_dr45.json",
    "10_metadata/tcga_escc_open_multiomics_file_manifest.tsv",
    "20_reusable/tcga_escc_open_multiomics_patient_availability.tsv",
    "GDC /data/<file_id>"
  ),
  generating_script = "scripts/11_download_ingest_tcga_escc_open_multiomics.R",
  key_parameters = c(
    "DR45.0 release_date=2025-12-04; seven exact filters; no data_release files clause",
    "one row per GDC file; preserve analysis and collapsed entity hierarchy",
    "one row per associated entity and matched case/sample/portion/analyte/aliquot path",
    "96-patient ESCC whitelist; per-layer availability; no complete-case requirement",
    "driver=95; five-layer=94; protein-deep=76",
    ".part resume; exact GDC size+MD5; local SHA256"
  ),
  software = paste0(
    "R ", getRversion(), "; data.table ", packageVersion("data.table"),
    "; httr2 ", packageVersion("httr2"), "; jsonlite ", packageVersion("jsonlite"),
    "; curl ", system2(curl_bin, "--version", stdout = TRUE, stderr = FALSE)[1L]
  ),
  generated_date = retrieval_date,
  regenerable = "yes"
)
provenance_path <- file.path(reusable_dir, "PROVENANCE.tsv")
atomic_fwrite(provenance, provenance_path)

dataset_md <- c(
  "# TCGA-ESCC GDC DR45 开放多组学数据说明",
  "",
  "## 基本信息",
  "",
  paste0("- `dataset_key`：`", dataset_key, "`。"),
  "- 来源：GDC `TCGA-ESCA`；数据发布 `45.0`，发布日期 `2025-12-04`。",
  "- 病例母集：96 例主疾病诊断确认的食管鳞状细胞癌（ESCC）。",
  "- 开放文件：RNA 95、miRNA 96、HM450 96、masked MAF 96、ASCAT2 gene 96、ASCAT2 segment 96、RPPA 78。",
  "- 当前完整性状态：`verified`；所有数据文件均通过 GDC 声明大小、GDC MD5 和本地 SHA256。",
  "",
  "## 目录与主要文件",
  "",
  "- `00_source/<layer>/<file_id>/<file_name>`：GDC 开放数据原始字节。",
  "- `10_metadata/gdc_files_<layer>_dr45.json`：七层 GDC `/files` 原始 API JSON。",
  "- `10_metadata/gdc_open_multiomics_queries_dr45.json`：白名单、精确过滤条件、fields 和断言。",
  "- `10_metadata/tcga_escc_open_multiomics_file_manifest.tsv`：每个 GDC 文件的展平元数据、实体概要和校验。",
  "- `10_metadata/tcga_escc_open_multiomics_entity_map.tsv`：`associated_entities` 到 case/sample/portion/analyte/aliquot 的完整映射。",
  "- `20_reusable/tcga_escc_open_multiomics_patient_availability.tsv`：96 例患者逐层可用性。",
  "- `20_reusable/tcga_escc_open_multiomics_intersections.tsv`：95 例驱动核心、94 例五层核心和 76 例蛋白深描子集。",
  "- `90_manifests/MANIFEST.tsv`：相对路径、大小、SHA256、来源和完整性。",
  "",
  "## 推荐使用场景",
  "",
  "- 体细胞突变/CNV—RNA/miRNA/甲基化—RPPA 的患者级驱动与调控轴候选。",
  "- 先做各层最大可用样本分析，再使用 95/94/76 完整个案子集进行跨层整合。",
  "",
  "## 纳入、排除与证据边界",
  "",
  "- `/files` 查询仅限 96 例白名单、`Primary Tumor`、`open`、`released`；病理白名单来自独立 DR45 临床入库。",
  "- MAF 和 CNV 每个文件保留肿瘤与配对正常两个 aliquot 来源；96 个配对正常中 88 个为 Blood Derived Normal、8 个为 Solid Tissue Normal；它们是调用/估计输入，不是表达层正常对照。",
  "- miRNA 有 96 个文件但只对应 95 例患者；重复 aliquot 需按预先规则选主，不得伪增患者数。",
  "- ASCAT2 gene 和 segment 是同一 CNV 语义的两个表示，不得当成两个独立组学证据。",
  "- RPPA 是 78 例子集；缺 RPPA 不得否定其他组学层证据。",
  "- 公开多组学整合用于假设生成和校准，不单独证明因果机制或治疗靶点。",
  "",
  "## 来源、引用与许可",
  "",
  "- GDC：https://portal.gdc.cancer.gov/projects/TCGA-ESCA",
  "- GDC API：https://api.gdc.cancer.gov/",
  "- 使用时应引用 GDC、TCGA ESCA 研究及对应数据管线，并遵守 GDC 条款。",
  paste0("- 使用项目：`", project_root, "`。")
)
dataset_path <- file.path(stage_root, "DATASET.md")
writeLines(dataset_md, dataset_path, useBytes = TRUE)

raw_api_paths <- c(
  status_path,
  vapply(query_results, function(x) x$raw_path, character(1))
)
derived_paths <- c(
  query_definition_path, file_manifest_path, entity_map_path, layer_checks_path, association_checks_path,
  availability_path, intersection_path, provenance_path, dataset_path
)
source_files <- file.path(stage_root, file_manifest$target_relative_path)
manifest_files <- c(raw_api_paths, derived_paths, source_files)
relative_paths <- as.character(path_rel(manifest_files, start = stage_root))
manifest <- data.table(
  relative_path = relative_paths,
  file_level = fifelse(grepl("/", relative_paths), sub("/.*$", "", relative_paths), "root"),
  size_bytes = as.numeric(file_info(manifest_files)$size),
  sha256 = vapply(manifest_files, sha256_file, character(1)),
  source_url = "generated within the project",
  download_date = retrieval_date,
  file_status = "generated_verified",
  corresponding_source_file = "10_metadata/gdc_files_<layer>_dr45.json",
  generation_method = "scripts/11_download_ingest_tcga_escc_open_multiomics.R",
  notes = ""
)
manifest[relative_path %in% path_rel(raw_api_paths, start = stage_root), `:=`(
  source_url = "https://api.gdc.cancer.gov/status or /files",
  file_status = "verified",
  corresponding_source_file = "",
  generation_method = "GDC API raw response retained without modification",
  notes = "DR45.0 frozen metadata"
)]
source_lookup <- file_manifest[, .(
  relative_path = target_relative_path,
  source_url,
  gdc_md5,
  file_id
)]
manifest[source_lookup, on = "relative_path", `:=`(
  source_url = i.source_url,
  file_status = "verified",
  corresponding_source_file = "",
  generation_method = "GDC /data/<file_id>; curl resumable download",
  notes = paste0("file_id=", i.file_id, ";gdc_md5=", i.gdc_md5)
)]
setorder(manifest, relative_path)
manifest_path <- file.path(manifest_dir, "MANIFEST.tsv")
atomic_fwrite(manifest, manifest_path)

# 不复制 .part 或任何未在 MANIFEST.tsv 中的暂存文件。
canonical_copy_sources <- c(manifest_files, manifest_path)
canonical_copy_relatives <- as.character(path_rel(canonical_copy_sources, start = stage_root))
dir_create(canonical_root, recurse = TRUE)
for (i in seq_along(canonical_copy_sources)) {
  source <- canonical_copy_sources[i]
  destination <- file.path(canonical_root, canonical_copy_relatives[i])
  dir_create(dirname(destination), recurse = TRUE)
  source_size <- as.numeric(file_info(source)$size)
  source_sha256 <- sha256_file(source)
  if (file_exists(destination)) {
    destination_matches <-
      as.numeric(file_info(destination)$size) == source_size &&
      sha256_file(destination) == source_sha256
    if (destination_matches) next
  }
  copy_temp <- paste0(destination, ".copying-", Sys.getpid())
  if (file_exists(copy_temp)) file_delete(copy_temp)
  if (!file.copy(source, copy_temp, overwrite = FALSE, copy.mode = TRUE, copy.date = TRUE)) {
    stop("复制到规范临时路径失败：", copy_temp)
  }
  if (as.numeric(file_info(copy_temp)$size) != source_size ||
      sha256_file(copy_temp) != source_sha256) {
    stop("规范临时副本大小或 SHA256 不一致：", copy_temp)
  }
  if (file_exists(destination)) file_delete(destination)
  file_move(copy_temp, destination)
  if (!file_exists(destination)) stop("无法原子提升规范副本：", destination)
}

observed_canonical_files <- sort(as.character(path_rel(
  dir_ls(canonical_root, recurse = TRUE, type = "file", all = TRUE),
  start = canonical_root
)))
expected_canonical_files <- sort(canonical_copy_relatives)
if (!identical(observed_canonical_files, expected_canonical_files)) {
  stop("规范目录存在未登记或缺失文件，拒绝登记")
}

canonical_manifest_path <- file.path(canonical_root, "90_manifests", "MANIFEST.tsv")
canonical_manifest <- fread(canonical_manifest_path, colClasses = "character", na.strings = NULL)
canonical_files <- file.path(canonical_root, canonical_manifest$relative_path)
if (!all(file_exists(canonical_files)) ||
    any(as.character(as.numeric(file_info(canonical_files)$size)) != canonical_manifest$size_bytes) ||
    any(vapply(canonical_files, sha256_file, character(1)) != canonical_manifest$sha256)) {
  stop("规范路径全量大小/SHA256 校验失败，不得登记")
}

catalog_row <- data.table(
  dataset_key = dataset_key,
  record_type = "patient_level_open_multiomics_dataset",
  access_level = "public",
  source = "GDC",
  accession = "TCGA-ESCA",
  version = snapshot_version,
  species = "Homo sapiens",
  disease = "esophageal squamous cell carcinoma",
  tissue = "primary tumor; paired normal aliquots retained only as MAF/CNV associated inputs",
  assay = "RNA-Seq; miRNA-Seq; HM450 methylation array; WXS masked MAF; ASCAT2 CNV; RPPA",
  sample_summary = "96 ESCC whitelist; files RNA=95 miRNA=96 HM450=96 MAF=96 ASCAT2 gene=96 segment=96 RPPA=78; intersections 95/94/76",
  available_levels = "00_source;10_metadata;20_reusable",
  status = "verified",
  local_path = canonical_root,
  manifest_path = canonical_manifest_path,
  source_url = "https://portal.gdc.cancer.gov/projects/TCGA-ESCA;https://api.gdc.cancer.gov/",
  download_date = retrieval_date,
  last_verified = retrieval_date,
  license_or_access = "public GDC files; cite GDC/TCGA and follow GDC data use terms",
  projects_using = project_root,
  recommended_use = "patient-level driver/transcription/epigenetic/protein integration with maximum-available and 95/94/76 complete-case subsets",
  limitations = "TCGA-ESCA mixed histology was filtered by DR45 clinical whitelist; miRNA has duplicate aliquot for one patient; MAF/CNV paired normals are associated calling inputs, not expression controls; RPPA is a 78-patient subset",
  notes = "Exact DR45 filters, raw API JSON, full associated entity hierarchy, GDC MD5 and local SHA256 retained."
)

project_datasets_path <- file.path(project_root, "data", "datasets.tsv")
project_row <- data.table(
  logical_name = "tcga_escc_gdc_open_multiomics_dr45",
  dataset_key = dataset_key,
  accession = "TCGA-ESCA",
  version = snapshot_version,
  data_level = "00_source;10_metadata;20_reusable",
  central_path = canonical_root,
  project_purpose = "ESCC 患者级驱动事件—转录/miRNA/甲基化—蛋白跨层轴及异质性",
  inclusion_exclusion = "仅 96 例 DR45 ESCC 白名单的 Primary Tumor 开放 released 文件；分层最大样本与 95/94/76 交集并行",
  limitations = "miRNA 一例双 aliquot；MAF/CNV 配对正常是调用输入而非表达正常；ASCAT2 gene/segment 不独立计数；RPPA 仅 78 例",
  local_link = "",
  last_verified = retrieval_date
)

update_global_and_project_indexes <- function() {
  catalog_lock <- paste0(catalog_path, ".update.lock")
  if (!dir.create(catalog_lock, showWarnings = FALSE)) {
    stop("CATALOG.tsv 正在被其他进程更新：", catalog_lock)
  }
  on.exit(unlink(catalog_lock, recursive = TRUE, force = TRUE), add = TRUE)
  catalog <- fread(catalog_path, colClasses = "character", na.strings = NULL)
  dataset_key_value <- dataset_key
  existing_catalog <- catalog[dataset_key == dataset_key_value]
  if (nrow(existing_catalog) == 0L) {
    catalog <- rbindlist(list(catalog, catalog_row), use.names = TRUE, fill = FALSE)
    atomic_fwrite(catalog, catalog_path)
  } else if (
    nrow(existing_catalog) != 1L ||
      existing_catalog$status != "verified" ||
      existing_catalog$local_path != canonical_root ||
      existing_catalog$manifest_path != canonical_manifest_path
  ) {
    stop("CATALOG.tsv 已有同 dataset_key 但不是同一 verified 规范记录")
  }

  # 项目索引在同一全局锁内重读和更新，避免 CATALOG 与项目定位表分裂。
  project_datasets <- fread(project_datasets_path, colClasses = "character", na.strings = NULL)
  existing_project <- project_datasets[dataset_key == dataset_key_value]
  if (nrow(existing_project) == 0L) {
    project_datasets <- rbindlist(list(project_datasets, project_row), use.names = TRUE, fill = FALSE)
    atomic_fwrite(project_datasets, project_datasets_path)
  } else if (nrow(existing_project) != 1L || existing_project$central_path != canonical_root) {
    stop("data/datasets.tsv 已有同 dataset_key 但规范路径不一致")
  }
}
update_global_and_project_indexes()

message(
  "已完成 TCGA-ESCC GDC DR45 七层开放数据下载、GDC MD5 + 本地 SHA256、",
  "associated entity 层级映射、95/94/76 交集、规范复制和双索引登记。"
)
