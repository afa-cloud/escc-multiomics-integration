#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(limma)
})

options(stringsAsFactors = FALSE, scipen = 999)

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("无法唯一定位当前脚本。")
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE)
project_root <- dirname(dirname(script_path))
if (!file.exists(file.path(project_root, "PROJECT_INDEX.md")) ||
    !dir.exists(file.path(project_root, "results")) ||
    !dir.exists(file.path(project_root, "_work", "checks"))) {
  stop("项目根目录标记或规范输出目录缺失：", project_root)
}

data_root <- Sys.getenv("RESEARCH_DATA_ROOT", path.expand("~/ResearchDataHub"))
target_dataset_key <- "METABOLOMICS_WORKBENCH_PR001876_retrieved_20260711"
canonical_root <- file.path(
  data_root, "datasets", "public", "Metabolomics_Workbench", "PR001876",
  "retrieved_20260711"
)
catalog_path <- file.path(data_root, "CATALOG.tsv")
manifest_path <- file.path(canonical_root, "90_manifests", "MANIFEST.tsv")
results_dir <- file.path(project_root, "results")
qa_dir <- file.path(project_root, "_work", "checks")

required_paths <- c(canonical_root, catalog_path, manifest_path)
if (!all(file.exists(required_paths) | dir.exists(required_paths))) {
  stop("缺少 PR001876 规范路径、CATALOG 或 manifest。")
}

catalog <- fread(catalog_path, colClasses = "character", na.strings = NULL)
catalog_row <- catalog[dataset_key == target_dataset_key]
if (nrow(catalog_row) != 1L || catalog_row$status != "verified" ||
    normalizePath(catalog_row$local_path, winslash = "/", mustWork = TRUE) !=
      normalizePath(canonical_root, winslash = "/", mustWork = TRUE)) {
  stop("PR001876 的 CATALOG 记录不唯一、未验证或规范路径不一致。")
}

analysis_spec <- data.table(
  study_id = c("ST003025", "ST003027", "ST003027"),
  analysis_id = c("AN004960", "AN004962", "AN004963"),
  platform = c("targeted GC-MS", "targeted LC-MS", "targeted LC-MS"),
  ion_mode = c("GC-EI positive", "positive", "negative"),
  transform = c("log2", "asinh", "asinh"),
  evidence_family = c(
    "ST003025_targeted_GCMS",
    "ST003027_targeted_LCMS_combined_modes",
    "ST003027_targeted_LCMS_combined_modes"
  ),
  within_analysis_evidence_unit = c(
    "ST003025_targeted_tissue_32",
    "ST003027_targeted_tissue_shared_32",
    "ST003027_targeted_tissue_shared_32"
  ),
  gate_independence_group = "PR001876_early_stage_tissue_subject_overlap_unresolved",
  matrix_relative = file.path(
    "20_reusable", paste0(c("AN004960", "AN004962", "AN004963"), "_metabolite_matrix.tsv")
  ),
  annotation_relative = file.path(
    "20_reusable", paste0(c("AN004960", "AN004962", "AN004963"), "_metabolite_annotations.tsv")
  ),
  metadata_relative = file.path(
    "10_metadata",
    c("ST003025_AN004960_mwtab.txt", "ST003027_AN004962_mwtab.txt", "ST003027_AN004963_mwtab.txt")
  ),
  raw_zip_relative = c(
    "00_source/ST003025_Rawfiles.zip",
    "00_source/ST003027_Rawfiles.zip",
    "00_source/ST003027_Rawfiles.zip"
  ),
  source_unit_conflict = c(
    "source_units_field_mz_but_method_comments_describe_ug_per_g",
    "source_units_field_mz_but_method_comments_describe_nmol_per_g",
    "source_units_field_mz_but_method_comments_describe_nmol_per_g"
  ),
  source_upstream_imputation = c(
    "author-processed matrix includes minimum_x_random_0.1_to_0.5 imputation; seed and imputed cells unavailable",
    "not described in source metadata",
    "not described in source metadata"
  )
)

sample_map_relative <- "20_reusable/pr001876_analysis_sample_map.tsv"
study_summary_relative <- "20_reusable/pr001876_study_summary.tsv"
input_relative <- unique(c(
  analysis_spec$matrix_relative, analysis_spec$annotation_relative,
  analysis_spec$metadata_relative, analysis_spec$raw_zip_relative,
  sample_map_relative, study_summary_relative
))

manifest <- fread(manifest_path, colClasses = "character", na.strings = NULL)
manifest_inputs <- manifest[match(input_relative, relative_path)]
if (nrow(manifest_inputs) != length(input_relative) || any(is.na(manifest_inputs$relative_path)) ||
    !identical(manifest_inputs$relative_path, input_relative)) {
  stop("目标输入未完整登记到 PR001876 manifest。")
}
input_paths <- file.path(canonical_root, input_relative)
if (!all(file.exists(input_paths))) stop("已登记的目标输入文件缺失。")
input_size <- as.character(as.numeric(file.info(input_paths)$size))
input_sha256 <- vapply(
  input_paths, function(path) digest(path, algo = "sha256", file = TRUE, serialize = FALSE),
  character(1)
)
if (!all(input_size == manifest_inputs$size_bytes) || !all(input_sha256 == manifest_inputs$sha256)) {
  stop("目标输入文件大小或 SHA256 与 manifest 不一致。")
}

gc_metadata_text <- readLines(
  file.path(canonical_root, "10_metadata/ST003025_AN004960_mwtab.txt"),
  warn = FALSE, encoding = "UTF-8"
)
gc_imputation_disclosure_verified <- any(grepl(
  "minimum value by a random number between 0.1 and 0.5", gc_metadata_text,
  fixed = TRUE
))
if (!gc_imputation_disclosure_verified) stop("无法在 GC mwTab 中复核作者上游随机插补说明。")

sample_map <- fread(
  file.path(canonical_root, sample_map_relative),
  colClasses = "character", na.strings = NULL
)
target_sample_map <- sample_map[analysis_id %chin% analysis_spec$analysis_id]
subject_id_present <- target_sample_map[
  !is.na(source_subject_id) & nzchar(trimws(source_subject_id)), .N
]
if (subject_id_present != 0L) stop("目标矩阵出现与已冻结边界不一致的 source subject ID。")

atomic_fwrite <- function(object, path) {
  temp_path <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  fwrite(object, temp_path, sep = "\t", quote = FALSE, na = "")
  if (!file.rename(temp_path, path)) stop("无法原子发布文件：", path)
}

atomic_write_lines <- function(lines, path) {
  temp_path <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  on.exit(if (file.exists(temp_path)) unlink(temp_path), add = TRUE)
  writeLines(lines, temp_path, useBytes = TRUE)
  if (!file.rename(temp_path, path)) stop("无法原子发布文件：", path)
}

annotation_column <- function(annotation, column_name) {
  if (column_name %chin% names(annotation)) as.character(annotation[[column_name]]) else rep("", nrow(annotation))
}

normalize_metabolite_name <- function(value) tolower(gsub("[^[:alnum:]]", "", value))

safe_variance <- function(value) {
  value <- value[!is.na(value)]
  if (length(value) < 2L) return(NA_real_)
  var(value)
}

safe_wilcoxon <- function(escc, normal, minimum_n = 2L) {
  escc <- escc[!is.na(escc)]
  normal <- normal[!is.na(normal)]
  if (length(escc) < minimum_n || length(normal) < minimum_n) {
    return(c(statistic = NA_real_, p_value = NA_real_, rank_biserial = NA_real_))
  }
  test <- suppressWarnings(wilcox.test(
    escc, normal, alternative = "two.sided", exact = FALSE, correct = FALSE
  ))
  statistic <- unname(test$statistic)
  c(
    statistic = statistic,
    p_value = test$p.value,
    rank_biserial = 2 * statistic / (length(escc) * length(normal)) - 1
  )
}

safe_fisher <- function(positive_escc, observed_escc, positive_normal, observed_normal) {
  if (observed_escc < 1L || observed_normal < 1L) return(NA_real_)
  p_value <- unname(fisher.test(matrix(
    c(positive_escc, observed_escc - positive_escc,
      positive_normal, observed_normal - positive_normal),
    nrow = 2L, byrow = TRUE
  ))$p.value)
  pmin(1, pmax(0, p_value))
}

fit_main_limma <- function(matrix_value, analysis_group) {
  group_factor <- factor(analysis_group, levels = c("Normal_tissue", "Early_stage_ESCC"))
  design <- model.matrix(~0 + group_factor)
  colnames(design) <- c("Normal_tissue", "Early_stage_ESCC")
  fit <- lmFit(matrix_value, design)
  contrast <- makeContrasts(
    Early_stage_ESCC_minus_Normal_tissue = Early_stage_ESCC - Normal_tissue,
    levels = design
  )
  fit <- eBayes(contrasts.fit(fit, contrast), robust = TRUE)
  top <- as.data.table(topTable(
    fit, coef = "Early_stage_ESCC_minus_Normal_tissue", number = Inf,
    sort.by = "none", confint = TRUE
  ), keep.rownames = "feature_id")
  top
}

fit_covariate_limma <- function(matrix_value, analysis_group, covariate) {
  group_factor <- factor(analysis_group, levels = c("Normal_tissue", "Early_stage_ESCC"))
  design <- model.matrix(~ scale(covariate) + group_factor)
  coefficient <- "group_factorEarly_stage_ESCC"
  fit <- eBayes(lmFit(matrix_value, design), robust = TRUE)
  top <- as.data.table(topTable(
    fit, coef = coefficient, number = Inf, sort.by = "none"
  ), keep.rownames = "feature_id")
  top[, .(feature_id, effect = logFC, p = P.Value)]
}

fit_paired_sensitivity <- function(matrix_value, analysis_group, pairing_candidate_id) {
  pair_levels <- unique(pairing_candidate_id)
  escc_columns <- match(
    pair_levels, pairing_candidate_id[analysis_group == "Early_stage_ESCC"]
  )
  normal_columns <- match(
    pair_levels, pairing_candidate_id[analysis_group == "Normal_tissue"]
  )
  escc_indices <- which(analysis_group == "Early_stage_ESCC")[escc_columns]
  normal_indices <- which(analysis_group == "Normal_tissue")[normal_columns]
  if (anyNA(escc_indices) || anyNA(normal_indices) || length(pair_levels) != 16L) {
    stop("低置信编号配对敏感性无法组成 16 对候选编号。")
  }
  paired_difference <- matrix_value[, escc_indices, drop = FALSE] -
    matrix_value[, normal_indices, drop = FALSE]
  design <- matrix(1, nrow = ncol(paired_difference), ncol = 1L,
                   dimnames = list(NULL, "Early_stage_ESCC_minus_Normal_tissue"))
  fit <- eBayes(lmFit(paired_difference, design), robust = TRUE)
  top <- as.data.table(topTable(
    fit, coef = "Early_stage_ESCC_minus_Normal_tissue", number = Inf, sort.by = "none"
  ), keep.rownames = "feature_id")
  top[, .(
    feature_id, effect = logFC, p = P.Value,
    complete_candidate_pairs = rowSums(!is.na(paired_difference))
  )]
}

direction_from_effect <- function(value) fifelse(
  value > 0, "higher_in_early_ESCC",
  fifelse(value < 0, "lower_in_early_ESCC", "no_direction")
)

zip_entries <- function(zip_path) utils::unzip(zip_path, list = TRUE)$Name

parse_gc_acquisition <- function(zip_path, raw_file_names) {
  entries <- zip_entries(zip_path)
  extraction_dir <- tempfile("pr001876_gc_cdf_")
  dir.create(extraction_dir, recursive = TRUE)
  on.exit(unlink(extraction_dir, recursive = TRUE, force = TRUE), add = TRUE)
  output <- lapply(raw_file_names, function(raw_name) {
    entry <- entries[basename(entries) == raw_name]
    if (length(entry) != 1L) stop("GC raw ZIP 中无法唯一定位：", raw_name)
    utils::unzip(zip_path, files = entry, exdir = extraction_dir, junkpaths = TRUE)
    extracted_path <- file.path(extraction_dir, raw_name)
    strings_output <- system2(
      "/usr/bin/strings", c("-a", shQuote(extracted_path)), stdout = TRUE, stderr = TRUE
    )
    strings_trimmed <- trimws(strings_output)
    stamp_label <- which(strings_trimmed == "experiment_date_time_stamp")
    if (length(stamp_label) != 1L) stop("GC CDF 中采集时间字段不唯一：", raw_name)
    stamp_window <- strings_trimmed[seq.int(
      stamp_label + 1L, min(length(strings_trimmed), stamp_label + 4L)
    )]
    stamp <- grep("^[0-9]{14}[+-][0-9]{4}$", stamp_window, value = TRUE)
    if (length(stamp) < 1L) stop("GC CDF 中未解析到采集时间：", raw_name)
    data.table(
      raw_file_name = raw_name,
      acquisition_time = as.POSIXct(substr(stamp[1], 1L, 14L), format = "%Y%m%d%H%M%S", tz = "UTC")
    )
  })
  rbindlist(output)
}

parse_lc_acquisition <- function(zip_path, raw_file_names) {
  entries <- zip_entries(zip_path)
  output <- lapply(raw_file_names, function(raw_name) {
    entry <- entries[basename(entries) == raw_name]
    if (length(entry) != 1L) stop("LC raw ZIP 中无法唯一定位：", raw_name)
    connection <- unz(zip_path, entry, open = "r")
    lines <- readLines(connection, warn = FALSE)
    close(connection)
    run_line <- grep("startTimeStamp=", lines, value = TRUE)
    if (length(run_line) < 1L) stop("LC mzML 中未解析到采集时间：", raw_name)
    stamp <- sub('.*startTimeStamp="([^"]+)".*', "\\1", run_line[1])
    data.table(
      raw_file_name = raw_name,
      acquisition_time = as.POSIXct(stamp, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")
    )
  })
  rbindlist(output)
}

gc_raw_zip <- file.path(canonical_root, "00_source/ST003025_Rawfiles.zip")
lc_raw_zip <- file.path(canonical_root, "00_source/ST003027_Rawfiles.zip")
gc_file_names <- unique(target_sample_map[analysis_id == "AN004960", raw_file_name])
lc_file_names <- unique(target_sample_map[analysis_id == "AN004962", raw_file_name])
gc_acquisition <- parse_gc_acquisition(gc_raw_zip, gc_file_names)
lc_acquisition <- parse_lc_acquisition(lc_raw_zip, lc_file_names)
if (nrow(gc_acquisition) != 32L || nrow(lc_acquisition) != 32L ||
    anyNA(gc_acquisition$acquisition_time) || anyNA(lc_acquisition$acquisition_time)) {
  stop("原始采集时间解析未覆盖各研究的 32 个样本。")
}

feature_qc_list <- vector("list", nrow(analysis_spec))
differential_list <- vector("list", nrow(analysis_spec))
inventory_list <- vector("list", nrow(analysis_spec))
run_order_list <- vector("list", nrow(analysis_spec))
paired_list <- vector("list", nrow(analysis_spec))
scale_list <- list()
sample_qc_list <- vector("list", nrow(analysis_spec))
sample_alignment_checks <- logical(nrow(analysis_spec))
analysis_maps <- vector("list", nrow(analysis_spec))
names(analysis_maps) <- analysis_spec$analysis_id
asinh_scales <- c(1e-3, 1e-2, 1e-1, 1, 10, 100, 1000)

for (index in seq_len(nrow(analysis_spec))) {
  spec <- analysis_spec[index]
  abundance <- fread(
    file.path(canonical_root, spec$matrix_relative), header = TRUE,
    check.names = FALSE, na.strings = c("", "NA")
  )
  annotation <- fread(
    file.path(canonical_root, spec$annotation_relative), header = TRUE,
    check.names = FALSE, na.strings = NULL
  )
  if (ncol(abundance) != 33L || names(abundance)[1] != "metabolite_name") {
    stop(spec$analysis_id, " 矩阵不是 1 列代谢物加 32 列样本。")
  }
  if (nrow(annotation) != nrow(abundance) ||
      !identical(as.character(annotation[[1]]), as.character(abundance[[1]]))) {
    stop(spec$analysis_id, " 注释表与矩阵的代谢物顺序不一致。")
  }
  if (anyDuplicated(abundance[[1]])) stop(spec$analysis_id, " 存在未解析的重复代谢物名。")

  sample_columns <- names(abundance)[-1]
  current_map <- target_sample_map[analysis_id == spec$analysis_id]
  if (nrow(current_map) != 32L || uniqueN(current_map$sample_id) != 32L ||
      !setequal(sample_columns, current_map$sample_id)) {
    stop(spec$analysis_id, " 矩阵列名与 sample map 无法一对一对齐。")
  }
  current_map <- current_map[match(sample_columns, sample_id)]
  if (any(is.na(current_map$sample_id))) stop(spec$analysis_id, " 样本列名 join 后出现缺失。")
  current_map[, analysis_group := fifelse(
    source_factor == "Fator:Early stage ESCC", "Early_stage_ESCC",
    fifelse(source_factor == "Fator:Normal tissue", "Normal_tissue", NA_character_)
  )]
  if (anyNA(current_map$analysis_group) ||
      current_map[analysis_group == "Early_stage_ESCC", .N] != 16L ||
      current_map[analysis_group == "Normal_tissue", .N] != 16L) {
    stop(spec$analysis_id, " 分组不是早期 ESCC 16 例和正常组织 16 例。")
  }
  if (!all(current_map$paired_model_role == "unpaired_primary_only") ||
      !all(current_map$pairing_confidence == "low") ||
      uniqueN(current_map$pairing_candidate_id) != 16L ||
      any(current_map[, .N, by = pairing_candidate_id]$N != 2L)) {
    stop(spec$analysis_id, " 非配对主模型或低置信编号配对边界与 sample map 不一致。")
  }

  acquisition <- if (spec$study_id == "ST003025") gc_acquisition else lc_acquisition
  current_map[, acquisition_time := acquisition$acquisition_time[
    match(raw_file_name, acquisition$raw_file_name)
  ]]
  if (anyNA(current_map$acquisition_time)) stop(spec$analysis_id, " 样本采集时间 join 失败。")
  current_map[, acquisition_rank := rank(acquisition_time, ties.method = "first")]
  ordered <- order(current_map$acquisition_rank)
  block <- integer(nrow(current_map))
  block[ordered] <- rleid(current_map$analysis_group[ordered])
  current_map[, acquisition_group_block := block]
  escc_map <- current_map$analysis_group == "Early_stage_ESCC"
  normal_map <- current_map$analysis_group == "Normal_tissue"
  perfect_group_block <-
    max(current_map$acquisition_rank[escc_map]) < min(current_map$acquisition_rank[normal_map]) ||
    max(current_map$acquisition_rank[normal_map]) < min(current_map$acquisition_rank[escc_map])
  rank_group_correlation <- cor(
    current_map$acquisition_rank, as.integer(current_map$analysis_group == "Early_stage_ESCC")
  )
  current_map[, `:=`(
    perfect_group_block = perfect_group_block,
    rank_group_correlation = rank_group_correlation
  )]
  sample_alignment_checks[index] <- TRUE
  analysis_maps[[spec$analysis_id]] <- copy(current_map)
  sample_qc_list[[index]] <- current_map[, .(
    study_id, analysis_id, sample_id, raw_file_name, raw_index, physical_specimen_key,
    analysis_group, pairing_candidate_id, pairing_basis, pairing_confidence,
    primary_model_role = paired_model_role,
    acquisition_time_utc = format(acquisition_time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
    acquisition_rank, acquisition_group_block, perfect_group_block,
    rank_group_correlation
  )]

  raw_matrix <- as.matrix(abundance[, -1, with = FALSE])
  storage.mode(raw_matrix) <- "numeric"
  rownames(raw_matrix) <- paste0(spec$analysis_id, "__", seq_len(nrow(raw_matrix)))
  feature_id <- rownames(raw_matrix)
  metabolite_name <- as.character(abundance[[1]])
  escc_index <- current_map$analysis_group == "Early_stage_ESCC"
  normal_index <- current_map$analysis_group == "Normal_tissue"

  nonmissing_total <- rowSums(!is.na(raw_matrix))
  nonmissing_escc <- rowSums(!is.na(raw_matrix[, escc_index, drop = FALSE]))
  nonmissing_normal <- rowSums(!is.na(raw_matrix[, normal_index, drop = FALSE]))
  zero_total <- rowSums(raw_matrix == 0, na.rm = TRUE)
  negative_total <- rowSums(raw_matrix < 0, na.rm = TRUE)
  positive_escc <- rowSums(raw_matrix[, escc_index, drop = FALSE] > 0, na.rm = TRUE)
  positive_normal <- rowSums(raw_matrix[, normal_index, drop = FALSE] > 0, na.rm = TRUE)
  positive_total <- positive_escc + positive_normal
  positive_rate_escc <- positive_escc / pmax(nonmissing_escc, 1L)
  positive_rate_normal <- positive_normal / pmax(nonmissing_normal, 1L)
  all_na <- nonmissing_total == 0L
  all_zero <- nonmissing_total > 0L & zero_total == nonmissing_total
  raw_variance <- apply(raw_matrix, 1L, safe_variance)
  zero_variance <- !is.na(raw_variance) & raw_variance == 0
  label_independent_detection_pass <- positive_total >= 8L
  estimable <- nonmissing_escc >= 4L & nonmissing_normal >= 4L
  primary_test <- !all_na & !all_zero & !zero_variance & label_independent_detection_pass & estimable

  exclusion_reason <- fifelse(
    all_na, "all_na",
    fifelse(all_zero, "all_zero",
      fifelse(zero_variance, "zero_variance",
        fifelse(!estimable, "fewer_than_4_nonmissing_in_one_group",
          fifelse(!label_independent_detection_pass, "fewer_than_8_positive_values_across_32", "pass")
        )
      )
    )
  )

  feature_qc_list[[index]] <- data.table(
    study_id = spec$study_id, analysis_id = spec$analysis_id,
    platform = spec$platform, ion_mode = spec$ion_mode,
    evidence_family = spec$evidence_family,
    within_analysis_evidence_unit = spec$within_analysis_evidence_unit,
    gate_independence_group = spec$gate_independence_group,
    feature_id, metabolite_name,
    normalized_name_key = normalize_metabolite_name(metabolite_name),
    n_escc_nonmissing = nonmissing_escc, n_normal_nonmissing = nonmissing_normal,
    positive_count_total = positive_total,
    positive_detection_rate_escc = positive_rate_escc,
    positive_detection_rate_normal = positive_rate_normal,
    missing_fraction_all = 1 - nonmissing_total / ncol(raw_matrix),
    zero_count_all = zero_total, negative_count_all = negative_total,
    censoring_sensitive = zero_total > 0L | negative_total > 0L,
    all_na, all_zero, zero_variance,
    primary_detection_gate = "label-independent: >=8 positive values across 32; >=4 nonmissing per group",
    primary_test, exclusion_reason,
    source_unit_conflict = spec$source_unit_conflict,
    source_upstream_imputation = spec$source_upstream_imputation,
    additional_imputation_by_this_script = "none",
    absolute_concentration_claim_allowed = "no"
  )
  if (!any(primary_test)) stop(spec$analysis_id, " 没有特征通过标签无关主门禁。")

  transformed_matrix <- if (spec$transform == "log2") {
    if (any(raw_matrix[primary_test, , drop = FALSE] <= 0, na.rm = TRUE)) {
      stop(spec$analysis_id, " log2 矩阵出现非正值。")
    }
    log2(raw_matrix)
  } else {
    asinh(raw_matrix)
  }
  tested_indices <- which(primary_test)
  tested_feature_ids <- feature_id[tested_indices]
  tested_transformed <- transformed_matrix[primary_test, , drop = FALSE]
  main_top <- fit_main_limma(tested_transformed, current_map$analysis_group)
  if (!identical(main_top$feature_id, tested_feature_ids)) {
    stop(spec$analysis_id, " limma 输出顺序与输入特征不一致。")
  }

  wilcoxon_matrix <- t(vapply(tested_indices, function(row_index) {
    safe_wilcoxon(raw_matrix[row_index, escc_index], raw_matrix[row_index, normal_index])
  }, numeric(3)))
  positive_wilcoxon_matrix <- t(vapply(tested_indices, function(row_index) {
    safe_wilcoxon(
      raw_matrix[row_index, escc_index][raw_matrix[row_index, escc_index] > 0],
      raw_matrix[row_index, normal_index][raw_matrix[row_index, normal_index] > 0],
      minimum_n = 4L
    )
  }, numeric(3)))
  detection_fisher_p <- vapply(tested_indices, function(row_index) {
    safe_fisher(
      positive_escc[row_index], nonmissing_escc[row_index],
      positive_normal[row_index], nonmissing_normal[row_index]
    )
  }, numeric(1))

  row_summary <- function(matrix_value, column_index, fun) apply(
    matrix_value[, column_index, drop = FALSE], 1L,
    function(value) if (all(is.na(value))) NA_real_ else fun(value, na.rm = TRUE)
  )
  annotation_tested <- annotation[tested_indices]
  hmdb_id <- annotation_column(annotation_tested, "HMDB ID")
  kegg_id <- annotation_column(annotation_tested, "KEGG ID")
  pubchem_id <- annotation_column(annotation_tested, "PubChem ID")
  identity_status <- fifelse(
    nzchar(hmdb_id) | nzchar(kegg_id) | nzchar(pubchem_id),
    "source_identifier_present_unverified",
    fifelse(grepl("?", metabolite_name[tested_indices], fixed = TRUE),
      "source_label_uncertain_unmapped", "source_label_unmapped")
  )

  differential_list[[index]] <- data.table(
    study_id = spec$study_id, analysis_id = spec$analysis_id,
    platform = spec$platform, ion_mode = spec$ion_mode,
    evidence_family = spec$evidence_family,
    within_analysis_evidence_unit = spec$within_analysis_evidence_unit,
    gate_independence_group = spec$gate_independence_group,
    feature_id = tested_feature_ids,
    metabolite_name = metabolite_name[tested_indices],
    normalized_name_key = normalize_metabolite_name(metabolite_name[tested_indices]),
    cas = annotation_column(annotation_tested, "CAS"),
    molecular_formula = annotation_column(annotation_tested, "molecular formula"),
    quantitated_mz = annotation_column(annotation_tested, "quantitated m/z"),
    hmdb_id, kegg_id, pubchem_id,
    annotation_identity_status = identity_status,
    internal_standard = annotation_column(annotation_tested, "Internal standard name"),
    source_annotation_ion_mode = annotation_column(annotation_tested, "Ionization mode"),
    precursor_mz = annotation_column(annotation_tested, "Precursor ion m/z"),
    product_mz = annotation_column(annotation_tested, "Product ion m/z"),
    n_escc_nonmissing = nonmissing_escc[tested_indices],
    n_normal_nonmissing = nonmissing_normal[tested_indices],
    positive_count_total = positive_total[tested_indices],
    positive_detection_rate_escc = positive_rate_escc[tested_indices],
    positive_detection_rate_normal = positive_rate_normal[tested_indices],
    detection_rate_difference = positive_rate_escc[tested_indices] - positive_rate_normal[tested_indices],
    zero_count_all = zero_total[tested_indices], negative_count_all = negative_total[tested_indices],
    censoring_sensitive = zero_total[tested_indices] > 0L | negative_total[tested_indices] > 0L,
    detection_fisher_p,
    positive_abundance_n_escc = positive_escc[tested_indices],
    positive_abundance_n_normal = positive_normal[tested_indices],
    positive_abundance_wilcoxon_w = positive_wilcoxon_matrix[, "statistic"],
    positive_abundance_rank_biserial = positive_wilcoxon_matrix[, "rank_biserial"],
    positive_abundance_wilcoxon_p = positive_wilcoxon_matrix[, "p_value"],
    raw_mean_escc = row_summary(raw_matrix[tested_indices, , drop = FALSE], escc_index, mean),
    raw_mean_normal = row_summary(raw_matrix[tested_indices, , drop = FALSE], normal_index, mean),
    raw_median_escc = row_summary(raw_matrix[tested_indices, , drop = FALSE], escc_index, median),
    raw_median_normal = row_summary(raw_matrix[tested_indices, , drop = FALSE], normal_index, median),
    primary_transform = spec$transform,
    primary_asinh_cofactor = ifelse(spec$transform == "asinh", 1, NA_real_),
    limma_effect_transformed = main_top$logFC,
    limma_effect_scale = ifelse(
      spec$transform == "log2",
      "log2 difference: Early_stage_ESCC minus Normal_tissue",
      "asinh(x*1) difference: Early_stage_ESCC minus Normal_tissue"
    ),
    limma_ci95_low = main_top$CI.L, limma_ci95_high = main_top$CI.R,
    limma_t = main_top$t, limma_p = main_top$P.Value,
    limma_q_analysis = p.adjust(main_top$P.Value, method = "BH"),
    wilcoxon_w = wilcoxon_matrix[, "statistic"],
    wilcoxon_rank_biserial = wilcoxon_matrix[, "rank_biserial"],
    wilcoxon_p = wilcoxon_matrix[, "p_value"],
    source_unit_conflict = spec$source_unit_conflict,
    source_upstream_imputation = spec$source_upstream_imputation,
    additional_imputation_by_this_script = "none",
    nonpositive_value_semantics = ifelse(
      spec$transform == "asinh",
      "source zero/negative values retained in primary asinh; detection/censoring flags and positive-only abundance sensitivity reported",
      "source author-processed positive values used"
    ),
    absolute_concentration_claim_allowed = "no",
    primary_design = "unpaired",
    paired_model_used_as_primary = "no",
    host_multiomics_sample_level_link_allowed = "no"
  )

  run_order <- fit_covariate_limma(
    tested_transformed, current_map$analysis_group, current_map$acquisition_rank
  )
  run_order[, `:=`(
    study_id = spec$study_id, analysis_id = spec$analysis_id,
    platform = spec$platform, ion_mode = spec$ion_mode,
    evidence_family = spec$evidence_family,
    acquisition_source = "timestamp parsed from source raw ZIP",
    perfect_group_block = perfect_group_block,
    rank_group_correlation = rank_group_correlation,
    interpretation = ifelse(
      perfect_group_block,
      "run_order_confounded_sensitivity_only_not_corrected_truth",
      "run_order_sensitivity_only_not_corrected_truth"
    )
  )]
  run_order_list[[index]] <- run_order

  paired <- fit_paired_sensitivity(
    tested_transformed, current_map$analysis_group, current_map$pairing_candidate_id
  )
  paired[, `:=`(
    study_id = spec$study_id, analysis_id = spec$analysis_id,
    platform = spec$platform, ion_mode = spec$ion_mode,
    evidence_family = spec$evidence_family,
    pairing_confidence = "low",
    pairing_basis = unique(current_map$pairing_basis),
    interpretation = "number-offset pairing sensitivity only; unpaired model remains primary"
  )]
  paired[, estimable := !is.na(p)]
  paired_list[[index]] <- paired

  if (spec$transform == "asinh") {
    scale_rows <- lapply(asinh_scales, function(scale_value) {
      scale_top <- fit_main_limma(
        asinh(raw_matrix[primary_test, , drop = FALSE] * scale_value),
        current_map$analysis_group
      )
      data.table(
        study_id = spec$study_id, analysis_id = spec$analysis_id,
        platform = spec$platform, ion_mode = spec$ion_mode,
        evidence_family = spec$evidence_family,
        feature_id = scale_top$feature_id,
        asinh_scale_multiplier = scale_value,
        primary_scale = scale_value == 1,
        effect = scale_top$logFC, p = scale_top$P.Value,
        interpretation = "scale sensitivity; rank-based Wilcoxon is invariant and used for upgrade"
      )
    })
    scale_list[[length(scale_list) + 1L]] <- rbindlist(scale_rows)
  }

  inventory_list[[index]] <- data.table(
    study_id = spec$study_id, analysis_id = spec$analysis_id,
    platform = spec$platform, ion_mode = spec$ion_mode,
    transform = spec$transform, evidence_family = spec$evidence_family,
    within_analysis_evidence_unit = spec$within_analysis_evidence_unit,
    gate_independence_group = spec$gate_independence_group,
    n_escc = sum(escc_index), n_normal = sum(normal_index),
    source_subject_ids_present = 0L, primary_design = "unpaired",
    features_total = nrow(raw_matrix), features_all_na = sum(all_na),
    features_all_zero_non_all_na = sum(all_zero),
    features_with_any_zero = sum(zero_total > 0L),
    features_with_any_negative = sum(negative_total > 0L),
    features_primary_test = sum(primary_test),
    positive_detection_gate = "label-independent: >=8 positive values across 32; >=4 nonmissing per group",
    source_unit_conflict = spec$source_unit_conflict,
    source_upstream_imputation = spec$source_upstream_imputation,
    additional_imputation_by_this_script = "none",
    raw_acquisition_time_verified = TRUE,
    perfect_group_run_order_block = perfect_group_block,
    rank_group_correlation = rank_group_correlation,
    absolute_concentration_claim_allowed = "no",
    current_role = "independent metabolic-module calibration only"
  )
}

feature_qc <- rbindlist(feature_qc_list, use.names = TRUE, fill = TRUE)
differential <- rbindlist(differential_list, use.names = TRUE, fill = TRUE)
inventory <- rbindlist(inventory_list, use.names = TRUE, fill = TRUE)
sample_qc <- rbindlist(sample_qc_list, use.names = TRUE, fill = TRUE)
run_order_sensitivity <- rbindlist(run_order_list, use.names = TRUE, fill = TRUE)
paired_sensitivity <- rbindlist(paired_list, use.names = TRUE, fill = TRUE)
scale_sensitivity <- rbindlist(scale_list, use.names = TRUE, fill = TRUE)

positive_map <- analysis_maps[["AN004962"]][order(as.integer(raw_index))]
negative_map <- analysis_maps[["AN004963"]][order(as.integer(raw_index))]
positive_negative_same_specimens <-
  identical(positive_map$raw_index, negative_map$raw_index) &&
  identical(positive_map$sample_id, negative_map$sample_id) &&
  identical(positive_map$physical_specimen_key, negative_map$physical_specimen_key) &&
  identical(positive_map$analysis_group, negative_map$analysis_group)
if (!positive_negative_same_specimens) stop("ST003027 正负模式未能对齐为同一组 32 个物理标本及分组。")

cross_mode_name_overlap <- differential[
  study_id == "ST003027", uniqueN(analysis_id), by = normalized_name_key
][V1 > 1L, .N]

differential[, `:=`(
  limma_q_evidence_family = p.adjust(limma_p, method = "BH"),
  wilcoxon_q_evidence_family = p.adjust(wilcoxon_p, method = "BH"),
  detection_fisher_q_evidence_family = p.adjust(detection_fisher_p, method = "BH"),
  positive_abundance_wilcoxon_q_evidence_family = p.adjust(
    positive_abundance_wilcoxon_p, method = "BH"
  )
), by = evidence_family]
differential[, wilcoxon_q_analysis := p.adjust(wilcoxon_p, method = "BH"), by = analysis_id]
differential[, direction := direction_from_effect(limma_effect_transformed)]
differential[, wilcoxon_direction := direction_from_effect(wilcoxon_rank_biserial)]
differential[, positive_abundance_direction := direction_from_effect(positive_abundance_rank_biserial)]
differential[, direction_concordant := direction == wilcoxon_direction & direction != "no_direction"]
differential[, positive_abundance_direction_concordant :=
  direction == positive_abundance_direction & direction != "no_direction"]

run_order_sensitivity[, q_evidence_family := p.adjust(p, method = "BH"), by = evidence_family]
run_order_sensitivity[, direction := direction_from_effect(effect)]
paired_sensitivity[, q_evidence_family := p.adjust(p, method = "BH"), by = evidence_family]
paired_sensitivity[, direction := direction_from_effect(effect)]
scale_sensitivity[, q_evidence_family := p.adjust(p, method = "BH"),
                  by = .(evidence_family, asinh_scale_multiplier)]
scale_sensitivity[, direction := direction_from_effect(effect)]

run_join <- run_order_sensitivity[, .(
  analysis_id, feature_id,
  run_order_sensitivity_effect = effect,
  run_order_sensitivity_p = p,
  run_order_sensitivity_q_evidence_family = q_evidence_family,
  run_order_sensitivity_direction = direction,
  run_order_perfect_group_block = perfect_group_block,
  run_order_rank_group_correlation = rank_group_correlation
)]
paired_join <- paired_sensitivity[, .(
  analysis_id, feature_id,
  paired_sensitivity_effect = effect,
  paired_sensitivity_p = p,
  paired_sensitivity_q_evidence_family = q_evidence_family,
  paired_sensitivity_direction = direction,
  paired_sensitivity_estimable = estimable,
  paired_sensitivity_complete_candidate_pairs = complete_candidate_pairs
)]
scale_summary <- scale_sensitivity[, .(
  scale_sensitivity_min_q = min(q_evidence_family),
  scale_sensitivity_max_q = max(q_evidence_family),
  scale_sensitivity_all_q05 = all(q_evidence_family <= 0.05),
  scale_sensitivity_all_q10 = all(q_evidence_family <= 0.10),
  scale_sensitivity_direction_stable = uniqueN(direction) == 1L && direction[1] != "no_direction"
), by = .(analysis_id, feature_id)]

differential <- merge(differential, run_join, by = c("analysis_id", "feature_id"), all.x = TRUE, sort = FALSE)
differential <- merge(differential, paired_join, by = c("analysis_id", "feature_id"), all.x = TRUE, sort = FALSE)
differential <- merge(differential, scale_summary, by = c("analysis_id", "feature_id"), all.x = TRUE, sort = FALSE)
differential[, `:=`(
  run_order_sensitivity_direction_concordant =
    direction == run_order_sensitivity_direction & direction != "no_direction",
  paired_sensitivity_direction_concordant =
    direction == paired_sensitivity_direction & direction != "no_direction"
)]
differential[, detection_driven := censoring_sensitive &
  !is.na(detection_fisher_q_evidence_family) & detection_fisher_q_evidence_family <= 0.10 &
  (is.na(positive_abundance_wilcoxon_q_evidence_family) |
     positive_abundance_wilcoxon_q_evidence_family > 0.10 |
     !positive_abundance_direction_concordant)]

differential[, primary_candidate_tier := fifelse(
  limma_q_evidence_family <= 0.05, "limma_fdr_strong",
  fifelse(limma_q_evidence_family <= 0.10, "limma_fdr_exploratory",
    fifelse(limma_p <= 0.05 & wilcoxon_p <= 0.05 & direction_concordant,
      "nominal_concordant", "not_selected")
  )
)]
differential[, robust_lc_upgrade :=
  study_id == "ST003027" & primary_candidate_tier == "limma_fdr_strong" &
  wilcoxon_q_evidence_family <= 0.05 & direction_concordant &
  scale_sensitivity_all_q05 & scale_sensitivity_direction_stable &
  run_order_sensitivity_q_evidence_family <= 0.10 &
  run_order_sensitivity_direction_concordant & !detection_driven &
  (!censoring_sensitive |
     (positive_abundance_wilcoxon_q_evidence_family <= 0.10 &
        positive_abundance_direction_concordant))]
differential[, final_candidate_tier := fifelse(
  primary_candidate_tier == "not_selected", "not_selected",
  fifelse(study_id == "ST003025", "run_order_confounded_conditional",
    fifelse(robust_lc_upgrade, "robust_rank_and_scale_fdr",
      fifelse(primary_candidate_tier %chin% c("limma_fdr_strong", "limma_fdr_exploratory"),
        "conditional_fdr", "nominal_concordant_conditional")
    )
  )
)]
differential[, conditional_reason := fcase(
  final_candidate_tier == "not_selected", "not_selected",
  final_candidate_tier == "run_order_confounded_conditional",
    "GC disease group is confounded with acquisition block; source matrix also contains unreconstructable random imputation",
  final_candidate_tier == "robust_rank_and_scale_fdr",
    "all prespecified rank, scale, run-order and censoring upgrade criteria met",
  detection_driven,
    "detection/censoring component is not supported by positive-only abundance sensitivity",
  !is.na(scale_sensitivity_all_q05) & !scale_sensitivity_all_q05,
    "asinh scale sensitivity does not remain family q<=0.05 at all seven scales",
  wilcoxon_q_evidence_family > 0.05,
    "scale-invariant Wilcoxon does not reach family q<=0.05",
  !run_order_sensitivity_direction_concordant |
    run_order_sensitivity_q_evidence_family > 0.10,
    "acquisition-order sensitivity is not supportive at family q<=0.10",
  primary_candidate_tier == "limma_fdr_exploratory",
    "main unpaired limma reaches family q<=0.10 but not q<=0.05",
  primary_candidate_tier == "nominal_concordant",
    "nominal limma and Wilcoxon concordance only",
  default = "conditional because at least one prespecified upgrade criterion is unmet"
)]
differential[, evidence_claim_ceiling := paste(
  "PR001876 tissue metabolic feature or module calibration;",
  "not an independent patient-level host multiomics association"
)]

if (any(differential[study_id == "ST003025" & primary_candidate_tier != "not_selected",
                     final_candidate_tier != "run_order_confounded_conditional"])) {
  stop("GC 候选未全部执行顺序混杂降级。")
}

candidate_metabolites <- differential[final_candidate_tier != "not_selected"]
setorder(candidate_metabolites, evidence_family, final_candidate_tier,
         limma_q_evidence_family, limma_p, analysis_id)
candidate_metabolites[, candidate_rank := seq_len(.N), by = evidence_family]

candidate_counts <- candidate_metabolites[, .(
  candidates_total = .N,
  robust_rank_and_scale_fdr = sum(final_candidate_tier == "robust_rank_and_scale_fdr"),
  conditional_fdr = sum(final_candidate_tier == "conditional_fdr"),
  run_order_confounded_conditional = sum(final_candidate_tier == "run_order_confounded_conditional"),
  nominal_concordant_conditional = sum(final_candidate_tier == "nominal_concordant_conditional")
), by = analysis_id]
inventory <- merge(inventory, candidate_counts, by = "analysis_id", all.x = TRUE, sort = FALSE)
for (column_name in c(
  "candidates_total", "robust_rank_and_scale_fdr", "conditional_fdr",
  "run_order_confounded_conditional", "nominal_concordant_conditional"
)) set(inventory, which(is.na(inventory[[column_name]])), column_name, 0L)
inventory[, positive_negative_modes_same_32_specimens := fifelse(
  study_id == "ST003027", positive_negative_same_specimens, NA
)]
inventory[, source_name_overlap_between_lc_modes := fifelse(
  study_id == "ST003027", cross_mode_name_overlap, NA_integer_
)]
setorder(inventory, study_id, analysis_id)
setorder(feature_qc, study_id, analysis_id, feature_id)
setorder(differential, evidence_family, limma_q_evidence_family, limma_p, analysis_id)
setorder(run_order_sensitivity, evidence_family, q_evidence_family, analysis_id)
setorder(paired_sensitivity, evidence_family, q_evidence_family, analysis_id)
setorder(scale_sensitivity, evidence_family, asinh_scale_multiplier, q_evidence_family, analysis_id)
setorder(sample_qc, study_id, analysis_id, acquisition_rank)

# 运行时硬断言：边界、家族、敏感性和数值域。
if (!all(sample_alignment_checks) || subject_id_present != 0L ||
    nrow(sample_qc) != 96L || uniqueN(sample_qc[, .(analysis_id, sample_id)]) != 96L ||
    any(sample_qc[, .N, by = .(analysis_id, analysis_group)]$N != 16L)) {
  stop("样本身份或 16/16 分组硬断言失败。")
}
if (!positive_negative_same_specimens ||
    !all(sample_qc[analysis_id == "AN004960", perfect_group_block]) ||
    any(sample_qc[analysis_id %chin% c("AN004962", "AN004963"), perfect_group_block])) {
  stop("正负模式同标本或真实采集区块硬断言失败。")
}
expected_test_counts <- c(AN004960 = 11L, AN004962 = 171L, AN004963 = 150L)
observed_test_counts <- differential[, .N, by = analysis_id]
if (!identical(
  observed_test_counts[match(names(expected_test_counts), analysis_id), N],
  unname(expected_test_counts)
)) stop("标签无关过滤后的预期检验数硬断言失败。")
if (nrow(scale_sensitivity) != sum(expected_test_counts[c("AN004962", "AN004963")]) * length(asinh_scales) ||
    nrow(run_order_sensitivity) != nrow(differential) ||
    nrow(paired_sensitivity) != nrow(differential)) {
  stop("尺度、顺序或低置信配对敏感性行数硬断言失败。")
}
probability_columns <- c(
  "limma_p", "limma_q_analysis", "limma_q_evidence_family", "wilcoxon_p",
  "wilcoxon_q_analysis", "wilcoxon_q_evidence_family", "detection_fisher_p",
  "run_order_sensitivity_p", "run_order_sensitivity_q_evidence_family",
  "paired_sensitivity_p", "paired_sensitivity_q_evidence_family"
)
for (column_name in probability_columns) {
  value <- differential[[column_name]]
  na_allowed <- column_name == "detection_fisher_p" || grepl("^paired_sensitivity", column_name)
  if (!na_allowed && anyNA(value)) stop(column_name, " 存在缺失。")
  if (grepl("^paired_sensitivity", column_name) &&
      any(is.na(value) != !differential$paired_sensitivity_estimable)) {
    stop(column_name, " 的缺失与 paired_sensitivity_estimable 标志不一致。")
  }
  if (any(value < -1e-12 | value > 1 + 1e-12, na.rm = TRUE)) stop(column_name, " 超出 0–1。")
}
if (sum(inventory$candidates_total) != nrow(candidate_metabolites) ||
    sum(candidate_counts[, robust_rank_and_scale_fdr + conditional_fdr +
          run_order_confounded_conditional + nominal_concordant_conditional]) !=
      nrow(candidate_metabolites)) {
  stop("候选分层计数硬断言失败。")
}

inventory_path <- file.path(results_dir, "pr001876_targeted_ms_analysis_inventory.tsv")
sample_qc_path <- file.path(results_dir, "pr001876_targeted_ms_sample_qc.tsv")
feature_qc_path <- file.path(results_dir, "pr001876_targeted_ms_feature_qc.tsv")
differential_path <- file.path(results_dir, "pr001876_targeted_ms_differential.tsv")
candidate_path <- file.path(results_dir, "pr001876_targeted_ms_candidate_metabolites.tsv")
run_order_path <- file.path(results_dir, "pr001876_targeted_ms_run_order_sensitivity.tsv")
scale_path <- file.path(results_dir, "pr001876_targeted_ms_scale_sensitivity.tsv")
paired_path <- file.path(results_dir, "pr001876_targeted_ms_paired_sensitivity.tsv")
summary_path <- file.path(results_dir, "pr001876_targeted_ms_summary.md")
artifact_manifest_path <- file.path(results_dir, "pr001876_targeted_ms_artifact_manifest.tsv")
qa_path <- file.path(qa_dir, "pr001876_targeted_ms_analysis_qa_20260711.md")

format_candidate <- function(row) paste0(
  "- `", row$analysis_id, "` / ", row$metabolite_name, "：",
  ifelse(row$direction == "higher_in_early_ESCC", "早期 ESCC 较高", "早期 ESCC 较低"),
  "；主模型家族 q=", formatC(row$limma_q_evidence_family, digits = 3, format = "fg"),
  "，Wilcoxon 家族 q=", formatC(row$wilcoxon_q_evidence_family, digits = 3, format = "fg"),
  "，最终层级 `", row$final_candidate_tier, "`。"
)
top_candidates <- head(
  candidate_metabolites[final_candidate_tier == "robust_rank_and_scale_fdr"][
    order(limma_q_evidence_family, wilcoxon_q_evidence_family)
  ],
  12L
)
top_candidate_lines <- if (nrow(top_candidates)) vapply(
  seq_len(nrow(top_candidates)), function(i) format_candidate(top_candidates[i]), character(1)
) else "- 无特征进入宽松候选清单。"
inventory_lines <- vapply(seq_len(nrow(inventory)), function(i) {
  row <- inventory[i]
  paste0(
    "- `", row$analysis_id, "`：", row$n_escc, " 例早期 ESCC + ", row$n_normal,
    " 例正常组织；", row$features_total, " 个面板特征，", row$features_primary_test,
    " 个进入主检验，", row$candidates_total, " 个进入宽松候选清单。"
  )
}, character(1))
tier_counts <- candidate_metabolites[, .N, by = final_candidate_tier][order(final_candidate_tier)]
tier_lines <- paste0("- `", tier_counts$final_candidate_tier, "`：", tier_counts$N, " 个。")

summary_lines <- c(
  "# PR001876 靶向组织代谢组分析摘要",
  "", "## 分析范围", "",
  "本次分析包括 ST003025 靶向 GC-MS 及 ST003027 靶向 LC-MS 正、负离子模式。公开元数据没有可核验的 subject ID；三个主模型均为 16 例早期 ESCC 对 16 例正常组织的非配对比较。",
  "", inventory_lines,
  "", "## 主分析与宽松门禁", "",
  "- 主过滤与标签无关：32 个样本中至少 8 个正值，且每组至少 4 个非缺失值；不要求任一疾病组单独达到 50%。",
  "- GC-MS 使用 `log2`；LC-MS 主尺度使用 `asinh(x×1)`。本脚本不新增插补。",
  "- 主检验为非配对 `limma`；原始值 Wilcoxon 为尺度不变的优先升级证据。ST003027 正、负模式共用一个 BH 家族。",
  "- 宽松候选保留家族 q≤0.10，或 limma 与 Wilcoxon 均名义显著且方向一致的条目；门禁不要求所有敏感性完全一致，但最终层级明确标出稳健或条件性。",
  "", "## 最终候选层级", "", tier_lines,
  "", "## 关键稳健方向（按最终表排序）", "", top_candidate_lines,
  "", "## 敏感性与不可识别边界", "",
  "- 原始 CDF 时间戳显示 GC 的 1–16 与 17–32 是连续疾病分组区块，组别与采集顺序不可辨识。加入顺序的模型只用于展示不稳定性，不是可信的‘顺序校正真值’；所有 GC 候选统一降级为 `run_order_confounded_conditional`。",
  "- 原始 mzML 时间戳显示 LC 样本并非按疾病组完整成块；正式表仍报告实际顺序敏感性及方向一致性。",
  "- LC 对 `asinh` 乘数 1e-3、1e-2、1e-1、1、10、100、1000 全部重跑；只有同时获得尺度不变 Wilcoxon 支持者才能升级为 `robust_rank_and_scale_fdr`。",
  "- 编号偏移 16 的配对仅为低置信候选关系；配对模型只作为敏感性，非配对模型始终是主模型。",
  "- LC 的 0/负值保留在主 `asinh` 矩阵中，同时报告检出率 Fisher 检验、仅正值丰度 Wilcoxon、`censoring_sensitive` 与 `detection_driven` 标志。",
  "", "## 来源处理与证据边界", "",
  "- GC mwTab 明确说明作者处理矩阵已把原始缺失值替换为‘最小值×0.1–0.5 随机数’，但未提供随机种子和被插补单元。因此本脚本能复算固定发布矩阵，不能完整重建 raw→matrix 随机步骤；‘本脚本无额外插补’不等于‘来源矩阵未插补’。",
  "- 来源 mwTab 的单位字段与方法说明不一致；结果不声称经核验的绝对浓度，也不跨平台比较原始数值大小。",
  "- ST003025 与 ST003027 的受试者重叠无法确认，不计为两个独立患者验证队列；PR001876 与宿主多组学队列也不建立患者级关联。",
  "- LC 名称和源注释尚未完成人工 HMDB/KEGG 身份核验；带问号或无稳定 ID 的条目不能直接升级为精确化合物通路证据。"
)

atomic_fwrite(inventory, inventory_path)
atomic_fwrite(sample_qc, sample_qc_path)
atomic_fwrite(feature_qc, feature_qc_path)
atomic_fwrite(differential, differential_path)
atomic_fwrite(candidate_metabolites, candidate_path)
atomic_fwrite(run_order_sensitivity, run_order_path)
atomic_fwrite(scale_sensitivity, scale_path)
atomic_fwrite(paired_sensitivity, paired_path)
atomic_write_lines(summary_lines, summary_path)

formal_paths <- c(
  inventory_path, sample_qc_path, feature_qc_path, differential_path, candidate_path,
  run_order_path, scale_path, paired_path, summary_path
)
artifact_manifest <- data.table(
  artifact = tools::file_path_sans_ext(basename(formal_paths)),
  relative_path = file.path("results", basename(formal_paths)),
  file_size_bytes = as.character(as.numeric(file.info(formal_paths)$size)),
  sha256 = vapply(formal_paths, function(path) digest(
    path, algo = "sha256", file = TRUE, serialize = FALSE
  ), character(1)),
  generated_date = "2026-07-11",
  generation_script = "scripts/16_analyze_pr001876_metabolomics.R",
  source_dataset_key = target_dataset_key,
  status = "verified_after_atomic_publish"
)
atomic_fwrite(artifact_manifest, artifact_manifest_path)

readback_rows <- c(
  inventory = nrow(fread(inventory_path, na.strings = NULL)),
  sample_qc = nrow(fread(sample_qc_path, na.strings = NULL)),
  feature_qc = nrow(fread(feature_qc_path, na.strings = NULL)),
  differential = nrow(fread(differential_path, na.strings = NULL)),
  candidates = nrow(fread(candidate_path, na.strings = NULL)),
  run_order = nrow(fread(run_order_path, na.strings = NULL)),
  scale = nrow(fread(scale_path, na.strings = NULL)),
  paired = nrow(fread(paired_path, na.strings = NULL))
)
expected_rows <- c(
  inventory = nrow(inventory), sample_qc = nrow(sample_qc), feature_qc = nrow(feature_qc),
  differential = nrow(differential), candidates = nrow(candidate_metabolites),
  run_order = nrow(run_order_sensitivity), scale = nrow(scale_sensitivity),
  paired = nrow(paired_sensitivity)
)
if (!identical(readback_rows, expected_rows)) stop("正式输出回读行数不一致。")
manifest_readback <- fread(artifact_manifest_path, colClasses = "character", na.strings = NULL)
manifest_paths <- file.path(project_root, manifest_readback$relative_path)
manifest_sha_readback <- vapply(manifest_paths, function(path) digest(
  path, algo = "sha256", file = TRUE, serialize = FALSE
), character(1))
if (!all(file.exists(manifest_paths)) ||
    !all(as.character(as.numeric(file.info(manifest_paths)$size)) == manifest_readback$file_size_bytes) ||
    !all(manifest_sha_readback == manifest_readback$sha256)) {
  stop("artifact manifest 回读哈希验证失败。")
}

qa_lines <- c(
  "# PR001876 靶向代谢组分析 QA（2026-07-11）",
  "", "本文件是运行时检查记录，不是项目当前状态源。",
  "", "## 输入与身份", "",
  paste0("- CATALOG dataset_key：`", target_dataset_key, "`，状态 `verified`。"),
  paste0("- 输入文件：", length(input_relative), " 个；矩阵、注释、mwTab、raw ZIP、sample map 和汇总表的 manifest 大小及 SHA256 全部通过。"),
  paste0("- GC 作者上游随机插补说明复核：", gc_imputation_disclosure_verified, "；本脚本额外插补：无。"),
  paste0("- 三个矩阵列名与 sample map 一对一 join：", all(sample_alignment_checks), "；source subject ID 非空数：", subject_id_present, "。"),
  paste0("- ST003027 正负模式样本、物理标本和分组完全一致：", positive_negative_same_specimens, "。"),
  "", "## 硬断言", "",
  paste0("- 每个分析均为 16/16，96 个分析—样本键唯一：TRUE。"),
  paste0("- 标签无关过滤后检验数：", paste(names(expected_test_counts), expected_test_counts, sep = "=", collapse = "；"), "。"),
  paste0("- GC 为完整疾病组采集区块：", all(sample_qc[analysis_id == "AN004960", perfect_group_block]),
         "；LC 为完整疾病组采集区块：", any(sample_qc[analysis_id == "AN004962", perfect_group_block]), "。"),
  paste0("- 顺序敏感性 ", nrow(run_order_sensitivity), " 行；低置信编号配对敏感性 ", nrow(paired_sensitivity),
         " 行；LC 七尺度敏感性 ", nrow(scale_sensitivity), " 行。"),
  "- 顺序敏感性只用于不稳定性审查；不可识别的 GC 顺序模型未被解释为校正真值。",
  paste0("- 主 limma/Wilcoxon 与顺序敏感性 P/q 均非缺失且在 0–1 内；低置信配对中不可估计数为 ",
         sum(!differential$paired_sensitivity_estimable), "，已由 `paired_sensitivity_estimable` 显式标记；检出率 P 允许在不可估计时为空。"),
  paste0("- 候选总数：", nrow(candidate_metabolites), "；分层合计与候选行数一致：TRUE。"),
  paste0("- GC 入选条目全部为 `run_order_confounded_conditional`：",
         all(candidate_metabolites[study_id == "ST003025", final_candidate_tier == "run_order_confounded_conditional"]), "。"),
  "", "## 正式输出与 manifest", "",
  paste0("- 预期行数：", paste(names(expected_rows), expected_rows, sep = "=", collapse = "；"), "。"),
  paste0("- 回读行数：", paste(names(readback_rows), readback_rows, sep = "=", collapse = "；"), "。"),
  paste0("- artifact manifest：", nrow(artifact_manifest), " 个正式文件；大小和 SHA256 回读一致：TRUE。"),
  "", "## 软件", "",
  paste0("- R ", getRversion(), "；data.table ", packageVersion("data.table"),
         "；limma ", packageVersion("limma"), "；digest ", packageVersion("digest"), "。")
)
atomic_write_lines(qa_lines, qa_path)

message(
  "PR001876 靶向 MS 审查后分析完成：", nrow(differential), " 个可检验特征，",
  nrow(candidate_metabolites), " 个宽松候选；主模型均为非配对。"
)
