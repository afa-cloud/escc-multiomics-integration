#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(fs)
  library(mclust)
  library(MultiAssayExperiment)
  library(randomForest)
  library(SummarizedExperiment)
})

options(stringsAsFactors = FALSE)

# 使用固定作者 ECMS 表达型 randomForest 模型，将 TCGA-ESCC GDC TPM
# 投影到 ECMS1–4。锁定模型对作者仓库内置 78 例矩阵的预测始终作为 anchor；额外
# 16 例只在预锁定的重叠校准通过时进入 primary 关联。本脚本不把
# MOFA 聚类类似性冒充 ECMS 标签，也不宣称临床分类器性能。

formal_filenames <- c(
  "tcga_escc_ecms_patient_probabilities.tsv",
  "tcga_escc_ecms_projection_calibration.tsv",
  "tcga_escc_ecms_projection_qa.tsv",
  "tcga_escc_ecms_factor_associations.tsv",
  "tcga_escc_ecms_adjusted_factor_progeny_associations.tsv",
  "tcga_escc_ecms_projection_summary.md"
)
manifest_filename <- "tcga_escc_ecms_projection_artifact_manifest.tsv"

args <- commandArgs(trailingOnly = TRUE)
fields_only <- "--fields-only" %in% args
unknown_args <- setdiff(args, "--fields-only")
if (length(unknown_args)) {
  stop("未知参数：", paste(unknown_args, collapse = ", "), call. = FALSE)
}
if (fields_only) {
  cat(paste(c(formal_filenames, manifest_filename), collapse = "\n"), "\n")
  cat(
    "patient_key_fields\tpatient_id,official_anchor_label,gdc_projection_label,resolved_ecms_label,label_source,extension_status,eligible_for_primary_association,prob_ECMS1-4,margin_custom,low_margin_custom_flag\n",
    "calibration_gates\tagreement>=0.75;kappa>=0.60;ARI>=0.60;median_class_probability_spearman>=0.50\n",
    "margin_boundary\tmargin<0.10 is project-defined soft uncertainty only; never rejects a repository-78 anchor prediction\n",
    "evidence_ceiling\tnon-level factors max T2 because ECMS/MOFA/PROGENy share TCGA RNA; frozen level factors are T0 technical/background only\n",
    sep = ""
  )
  quit(save = "no", status = 0L)
}

fail_if <- function(condition, message) {
  if (length(condition) != 1L || is.na(condition) || isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

as_logical_strict <- function(values, label) {
  fail_if(any(!values %in% c("TRUE", "FALSE")),
          paste(label, "必须严格为 TRUE/FALSE。"))
  values == "TRUE"
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
fail_if(length(script_argument) != 1L,
        "无法从 --file 唯一定位 scripts/23_project_tcga_escc_ecms.R。")
script_path <- normalizePath(
  sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE
)
execution_script_sha256 <- digest(
  script_path, algo = "sha256", file = TRUE, serialize = FALSE
)
project_root <- normalizePath(
  file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE
)
results_dir <- file.path(project_root, "results")
work_intermediate_dir <- file.path(project_root, "_work", "intermediate")
data_root <- Sys.getenv("RESEARCH_DATA_ROOT", path.expand("~/ResearchDataHub"))
catalog_path <- file.path(data_root, "CATALOG.tsv")

dataset_key <- paste0(
  "GITHUB_CITYUHK_COMPUTATIONAL_BIOLOGY_ESCC_CMS_",
  "8a49db65f86f2c1e2a50d9bc4cb742ccc868e4c6"
)
locked_commit <- "8a49db65f86f2c1e2a50d9bc4cb742ccc868e4c6"
locked_model_size <- 1608173
locked_model_sha256 <-
  "4c05645d712ec630f8ed717cba347dced4228a999c1b479ebd760bf1436e1ee4"
locked_feature_sha256 <-
  "9e81299a2e93f75f4e8375f852ca2bdd156126271f628d3c118b015d48e4120d"
ecms_levels <- paste0("ECMS", 1:4)

required_project_inputs <- file.path(project_root, c(
  "PROJECT_INDEX.md",
  "data/datasets.tsv",
  "results/tcga_escc_multiassay_dr45.rds",
  "results/tcga_escc_multiassay_analysis_sets.tsv",
  "results/tcga_escc_multiassay_artifact_manifest.tsv",
  "results/tcga_escc_mofa_factor_scores.tsv",
  "results/tcga_escc_progeny_pathway_scores.tsv",
  "results/tcga_escc_mofa_level_factor_qc.tsv",
  "results/tcga_escc_heterogeneity_artifact_manifest.tsv"
))
fail_if(any(!file_exists(required_project_inputs)), paste(
  "缺少 ECMS 投影输入：",
  paste(required_project_inputs[!file_exists(required_project_inputs)], collapse = ";")
))
fail_if(!dir_exists(data_root) || !file_exists(catalog_path),
        paste("ResearchDataHub 或 CATALOG 不可读：", data_root))

dir_create(work_intermediate_dir, recurse = TRUE)
stage_dir <- tempfile(pattern = ".tcga_ecms_projection_", tmpdir = work_intermediate_dir)
dir_create(stage_dir)
on.exit(if (dir_exists(stage_dir)) dir_delete(stage_dir), add = TRUE)

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

atomic_fwrite <- function(object, path) {
  dir_create(dirname(path), recurse = TRUE)
  temp_path <- tempfile(pattern = paste0(".", basename(path), "."),
                        tmpdir = dirname(path))
  on.exit(if (file_exists(temp_path)) file_delete(temp_path), add = TRUE)
  fwrite(object, temp_path, sep = "\t", quote = FALSE, na = "")
  reread <- fread(temp_path, colClasses = "character", na.strings = NULL)
  fail_if(nrow(reread) != nrow(object), paste("原子 TSV 回读失败：", path))
  fail_if(!file.rename(temp_path, path), paste("无法原子更新：", path))
}

atomic_write_lines <- function(lines, path) {
  dir_create(dirname(path), recurse = TRUE)
  temp_path <- tempfile(pattern = paste0(".", basename(path), "."),
                        tmpdir = dirname(path))
  on.exit(if (file_exists(temp_path)) file_delete(temp_path), add = TRUE)
  writeLines(lines, temp_path, useBytes = TRUE)
  fail_if(!file.rename(temp_path, path), paste("无法原子更新：", path))
}

atomic_publish_file <- function(source, destination) {
  dir_create(dirname(destination), recurse = TRUE)
  temp_path <- tempfile(
    pattern = paste0(".", basename(destination), ".publishing."),
    tmpdir = dirname(destination)
  )
  on.exit(if (file_exists(temp_path)) file_delete(temp_path), add = TRUE)
  fail_if(!file.copy(source, temp_path, overwrite = FALSE, copy.mode = TRUE),
          paste("无法复制到发布临时文件：", destination))
  fail_if(as.numeric(file_info(source)$size) != as.numeric(file_info(temp_path)$size) ||
            sha256_file(source) != sha256_file(temp_path),
          paste("发布临时文件大小或 SHA256 不一致：", destination))
  fail_if(!file.rename(temp_path, destination),
          paste("无法原子发布：", destination))
  fail_if(sha256_file(source) != sha256_file(destination),
          paste("发布后 SHA256 不一致：", destination))
}

verify_manifest <- function(root, manifest_path) {
  manifest <- fread(manifest_path, colClasses = "character", na.strings = NULL)
  required <- c("relative_path", "size_bytes", "sha256", "file_status")
  fail_if(!all(required %in% names(manifest)) || !nrow(manifest) ||
            anyDuplicated(manifest$relative_path),
          "ECMS 模型 MANIFEST.tsv 字段、行数或唯一性失败。")
  paths <- file.path(root, manifest$relative_path)
  fail_if(any(!file_exists(paths)), "ECMS 模型 manifest 登记文件缺失。")
  observed_size <- as.character(as.numeric(file_info(paths)$size))
  observed_sha <- vapply(paths, sha256_file, character(1))
  fail_if(any(observed_size != manifest$size_bytes) || any(observed_sha != manifest$sha256),
          "ECMS 模型 manifest 全文件大小/SHA256 校验失败。")
  observed_relatives <- sort(as.character(path_rel(
    dir_ls(root, recurse = TRUE, type = "file", all = TRUE), start = root
  )))
  expected_relatives <- sort(c(
    manifest$relative_path, as.character(path_rel(manifest_path, start = root))
  ))
  fail_if(!identical(observed_relatives, expected_relatives),
          "ECMS 模型规范目录存在 manifest 未登记文件。")
  manifest
}

verify_project_artifact_manifest <- function(manifest_path, required_relatives) {
  manifest <- fread(manifest_path, colClasses = "character", na.strings = NULL)
  required_fields <- c(
    "relative_path", "file_size_bytes", "sha256", "status"
  )
  fail_if(!all(required_fields %in% names(manifest)) || !nrow(manifest) ||
            anyDuplicated(manifest$relative_path),
          paste("上游 artifact manifest 字段/唯一性失败：", manifest_path))
  fail_if(!all(required_relatives %in% manifest$relative_path),
          paste("上游 manifest 未覆盖锁定输入：", manifest_path))
  paths <- file.path(project_root, manifest$relative_path)
  fail_if(any(!file_exists(paths)),
          paste("上游 manifest 登记 artifact 缺失：", manifest_path))
  observed_size <- as.character(as.numeric(file_info(paths)$size))
  observed_sha <- vapply(paths, sha256_file, character(1))
  fail_if(any(observed_size != manifest$file_size_bytes) ||
            any(observed_sha != manifest$sha256) ||
            any(manifest$status != "verified"),
          paste("上游 manifest 全 artifact 大小/SHA256/status 失败：",
                manifest_path))
  manifest
}

feature_list_sha256 <- function(features) {
  payload <- charToRaw(paste0(paste(features, collapse = "\n"), "\n"))
  digest(payload, algo = "sha256", serialize = FALSE)
}

probability_margin <- function(probability_matrix) {
  sorted <- t(apply(probability_matrix, 1L, sort, decreasing = TRUE))
  margin <- as.numeric(sorted[, 1L] - sorted[, 2L])
  names(margin) <- rownames(probability_matrix)
  margin
}

cohen_kappa <- function(reference, prediction, levels) {
  reference <- factor(reference, levels = levels)
  prediction <- factor(prediction, levels = levels)
  tab <- table(reference, prediction)
  observed <- sum(diag(tab)) / sum(tab)
  expected <- sum(rowSums(tab) * colSums(tab)) / sum(tab)^2
  if (!is.finite(expected) || expected >= 1) return(NA_real_)
  (observed - expected) / (1 - expected)
}

safe_scale <- function(values, label) {
  fail_if(any(!is.finite(values)) || length(values) < 3L,
          paste("无法标准化：", label))
  standard_deviation <- sd(values)
  fail_if(!is.finite(standard_deviation) || standard_deviation <= 0,
          paste("零方差无法标准化：", label))
  as.numeric(scale(values))
}

multiassay_manifest_path <- file.path(
  results_dir, "tcga_escc_multiassay_artifact_manifest.tsv"
)
heterogeneity_manifest_path <- file.path(
  results_dir, "tcga_escc_heterogeneity_artifact_manifest.tsv"
)
level_factor_qc_path <- file.path(
  results_dir, "tcga_escc_mofa_level_factor_qc.tsv"
)
multiassay_upstream_manifest <- verify_project_artifact_manifest(
  multiassay_manifest_path,
  file.path("results", c(
    "tcga_escc_multiassay_dr45.rds",
    "tcga_escc_multiassay_analysis_sets.tsv"
  ))
)
heterogeneity_upstream_manifest <- verify_project_artifact_manifest(
  heterogeneity_manifest_path,
  file.path("results", c(
    "tcga_escc_mofa_factor_scores.tsv",
    "tcga_escc_progeny_pathway_scores.tsv",
    "tcga_escc_mofa_level_factor_qc.tsv"
  ))
)

message("[1/7] 验证 ResearchDataHub 锁定 ECMS 模型")
catalog <- fread(catalog_path, colClasses = "character", na.strings = NULL)
dataset_key_value <- dataset_key
catalog_row <- catalog[dataset_key == dataset_key_value]
fail_if(nrow(catalog_row) != 1L || catalog_row$status != "verified" ||
          catalog_row$version != paste0("commit_", locked_commit),
        "CATALOG 中缺少唯一 verified ECMS 锁定 commit 记录。")
model_root <- catalog_row$local_path[[1L]]
model_manifest_path <- catalog_row$manifest_path[[1L]]
fail_if(!dir_exists(model_root) || !file_exists(model_manifest_path),
        "ECMS 规范目录或 manifest 缺失。")
model_manifest <- verify_manifest(model_root, model_manifest_path)
model_path <- file.path(model_root, "00_source", "ECMS.model.rdata")
model_manifest_row <- model_manifest[relative_path == "00_source/ECMS.model.rdata"]
fail_if(nrow(model_manifest_row) != 1L ||
          as.numeric(file_info(model_path)$size) != locked_model_size ||
          sha256_file(model_path) != locked_model_sha256 ||
          model_manifest_row$sha256 != locked_model_sha256,
        "ECMS.model.rdata 未通过锁定大小/SHA256 校验。")

message("[2/7] 验证模型对象与作者仓库内置 TCGA 78 例矩阵预测")
model_environment <- new.env(parent = emptyenv())
loaded_objects <- load(model_path, envir = model_environment)
expected_objects <- c(
  "rf.cl", "gse53625.val.df", "gse45670.val.df", "tcga.val.df", "gene.features"
)
fail_if(!setequal(loaded_objects, expected_objects),
        "ECMS.model.rdata 对象集与锁定版本不一致。")
rf_classifier <- model_environment$rf.cl
tcga_official_matrix <- model_environment$tcga.val.df
gene_features <- as.character(model_environment$gene.features)
fail_if(!inherits(rf_classifier, "randomForest") || rf_classifier$ntree != 500L ||
          rf_classifier$mtry != 17L || !identical(rf_classifier$classes, ecms_levels) ||
          length(rf_classifier$forest$cutoff) != 4L ||
          any(abs(rf_classifier$forest$cutoff - 0.25) > 1e-12),
        "ECMS randomForest 类型、ntree、mtry、类别或 0.25 cutoff 不符。")
fail_if(length(gene_features) != 314L || anyNA(gene_features) ||
          any(!nzchar(gene_features)) || anyDuplicated(gene_features) ||
          feature_list_sha256(gene_features) != locked_feature_sha256,
        "ECMS 314 特征列表的数量、唯一性或 SHA256 不符。")
fail_if(!is.matrix(tcga_official_matrix) ||
          !identical(dim(tcga_official_matrix), c(78L, 314L)) ||
          !identical(colnames(tcga_official_matrix), gene_features) ||
          anyDuplicated(rownames(tcga_official_matrix)) ||
          any(!is.finite(tcga_official_matrix)),
        "内置 tcga.val.df 维度、特征顺序、患者唯一性或数值失败。")
tcga_official_column_means <- colMeans(tcga_official_matrix)
tcga_official_column_sds <- apply(tcga_official_matrix, 2L, sd)
fail_if(max(abs(tcga_official_column_means)) > 1e-6 ||
          max(abs(tcga_official_column_sds - 1)) > 1e-6,
        "作者仓库内置 tcga.val.df 未保持逐基因均值 0/SD 1。")

official_label <- predict(rf_classifier, tcga_official_matrix, type = "response")
official_probability <- predict(rf_classifier, tcga_official_matrix, type = "prob")
fail_if(!identical(colnames(official_probability), ecms_levels) ||
          !identical(rownames(official_probability), rownames(tcga_official_matrix)) ||
          any(abs(rowSums(official_probability) - 1) > 1e-12),
        "锁定模型对作者仓库内置 TCGA 矩阵的四类概率失败。")
official_counts <- table(factor(official_label, levels = ecms_levels))
fail_if(!identical(as.integer(official_counts), c(23L, 34L, 7L, 14L)),
        "锁定模型对作者仓库内置 TCGA 矩阵的预测计数未复现 23/34/7/14。")

message("[3/7] 用 GDC TPM 重建 94 例队列内 ECMS 输入")
mae_path <- file.path(results_dir, "tcga_escc_multiassay_dr45.rds")
analysis_sets_path <- file.path(results_dir, "tcga_escc_multiassay_analysis_sets.tsv")
factor_path <- file.path(results_dir, "tcga_escc_mofa_factor_scores.tsv")
progeny_path <- file.path(results_dir, "tcga_escc_progeny_pathway_scores.tsv")

mae <- readRDS(mae_path)
validObject(mae)
analysis_sets <- fread(analysis_sets_path, showProgress = FALSE)
patients <- analysis_sets[
  analysis_set == "five_layer_core" & included == TRUE, patient_id
]
fail_if(length(patients) != 94L || uniqueN(patients) != 94L,
        "five-layer core 必须为 94 位唯一患者。")
fail_if(!"RNA" %in% names(experiments(mae)), "MultiAssayExperiment 缺少 RNA 视图。")
rna_se <- experiments(mae)[["RNA"]]
fail_if(!"tpm" %in% assayNames(rna_se), "RNA 视图缺少 tpm assay。")
fail_if(any(!patients %in% colnames(rna_se)), "RNA TPM 未覆盖 94 例 five-layer core。")

rna_gene_ids <- rownames(rna_se)
stripped_gene_ids <- sub("[.][0-9]+$", "", rna_gene_ids)
fail_if(anyDuplicated(stripped_gene_ids),
        "RNA Ensembl ID 去版本后不唯一；不允许静默合并。")
feature_index <- match(gene_features, stripped_gene_ids)
fail_if(anyNA(feature_index) || length(unique(feature_index)) != 314L,
        "GDC TPM 未唯一覆盖全部 314 个 ECMS 特征。")

tpm_matrix <- assay(rna_se, "tpm")[feature_index, patients, drop = FALSE]
storage.mode(tpm_matrix) <- "double"
rownames(tpm_matrix) <- gene_features
fail_if(any(!is.finite(tpm_matrix)) || any(tpm_matrix < 0),
        "ECMS 特征 TPM 含非有限值或负值。")
log_tpm <- log2(tpm_matrix + 1)
feature_sd <- apply(log_tpm, 1L, sd)
fail_if(length(feature_sd) != 314L || any(!is.finite(feature_sd)) || any(feature_sd <= 0),
        "log2(TPM+1) 后存在零方差 ECMS 特征。")

gdc_input <- scale(t(log_tpm), center = TRUE, scale = TRUE)
fail_if(!identical(dim(gdc_input), c(94L, 314L)) ||
          !identical(rownames(gdc_input), patients) ||
          !identical(colnames(gdc_input), gene_features) ||
          any(!is.finite(gdc_input)),
        "GDC ECMS 队列内逐基因 Z 标准化失败。")

gdc_label <- predict(rf_classifier, gdc_input, type = "response")
gdc_probability <- predict(rf_classifier, gdc_input, type = "prob")
fail_if(!identical(colnames(gdc_probability), ecms_levels) ||
          !identical(rownames(gdc_probability), patients) ||
          any(abs(rowSums(gdc_probability) - 1) > 1e-12),
        "GDC TPM 重投影四类概率失败。")

message("[4/7] 在 78 例重叠患者上执行预锁定校准门禁")
official_ids <- rownames(tcga_official_matrix)
fail_if(length(official_ids) != 78L || !all(official_ids %in% patients),
        "作者仓库内置 TCGA 78 例矩阵 rownames 无法完整映射 five-layer core。")
gdc_overlap_label <- gdc_label[match(official_ids, rownames(gdc_input))]
gdc_overlap_probability <- gdc_probability[official_ids, , drop = FALSE]
official_label_character <- as.character(official_label)
gdc_overlap_character <- as.character(gdc_overlap_label)

agreement <- mean(official_label_character == gdc_overlap_character)
kappa <- cohen_kappa(official_label_character, gdc_overlap_character, ecms_levels)
ari <- adjustedRandIndex(official_label_character, gdc_overlap_character)
class_probability_spearman <- vapply(ecms_levels, function(class_name) {
  suppressWarnings(cor(
    official_probability[, class_name], gdc_overlap_probability[, class_name],
    method = "spearman", use = "complete.obs"
  ))
}, numeric(1))
median_probability_spearman <- median(class_probability_spearman, na.rm = TRUE)
official_margin <- probability_margin(official_probability)
gdc_overlap_margin <- probability_margin(gdc_overlap_probability)
margin_spearman <- suppressWarnings(cor(
  official_margin, gdc_overlap_margin, method = "spearman", use = "complete.obs"
))

gate_values <- c(
  exact_label_agreement = agreement,
  cohen_kappa = kappa,
  adjusted_rand_index = ari,
  median_class_probability_spearman = median_probability_spearman
)
gate_thresholds <- c(
  exact_label_agreement = 0.75,
  cohen_kappa = 0.60,
  adjusted_rand_index = 0.60,
  median_class_probability_spearman = 0.50
)
gate_pass <- is.finite(gate_values) & gate_values >= gate_thresholds
calibration_pass <- all(gate_pass)

calibration <- data.table(
  metric = names(gate_values),
  observed = as.numeric(gate_values),
  threshold = as.numeric(gate_thresholds),
  direction = ">=",
  gate_role = "prelocked_extension_gate",
  pass = as.logical(gate_pass),
  notes = c(
    "locked-model predictions on repository-bundled 78 matrix versus GDC TPM reprojection",
    "manual Cohen kappa on fixed ECMS1-4 levels",
    "mclust adjustedRandIndex",
    "median of four class-wise Spearman correlations across 78 patients"
  )
)

recall_rows <- rbindlist(lapply(ecms_levels, function(class_name) {
  in_class <- official_label_character == class_name
  data.table(
    metric = paste0("diagnostic_recall_", class_name),
    observed = mean(gdc_overlap_character[in_class] == class_name),
    threshold = NA_real_, direction = "diagnostic_only",
    gate_role = "diagnostic_not_veto", pass = NA,
    notes = paste0("repository-anchor class n=", sum(in_class),
                   "; small classes do not veto")
  )
}))
probability_rows <- data.table(
  metric = paste0("class_probability_spearman_", ecms_levels),
  observed = as.numeric(class_probability_spearman),
  threshold = NA_real_, direction = "diagnostic_only",
  gate_role = "component_of_median_gate", pass = NA,
  notes = "repository-78 anchor prediction versus GDC reprojection probability vector"
)
margin_rows <- data.table(
  metric = c(
    "official_margin_median", "gdc_reprojection_margin_median",
    "margin_spearman", "margin_median_absolute_difference",
    "overall_extension_calibration_pass"
  ),
  observed = c(
    median(official_margin), median(gdc_overlap_margin), margin_spearman,
    median(abs(official_margin - gdc_overlap_margin)), as.numeric(calibration_pass)
  ),
  threshold = c(NA, NA, NA, NA, 1),
  direction = c(rep("diagnostic_only", 4), "all_four_prelocked_gates"),
  gate_role = c(rep("diagnostic_not_veto", 4), "extension_decision"),
  pass = c(rep(NA, 4), calibration_pass),
  notes = c(
    "top1-top2 probability margin; no official threshold",
    "top1-top2 probability margin; no official threshold",
    "Spearman correlation of two project-reported margins",
    "absolute margin difference; descriptive only",
    "repository-78 locked-model anchor predictions are retained; only additional 16 depend on this decision"
  )
)
calibration <- rbindlist(
  list(calibration, recall_rows, probability_rows, margin_rows), use.names = TRUE
)

patient_probabilities <- data.table(patient_id = patients)
patient_probabilities[, in_official_78 := patient_id %in% official_ids]
patient_probabilities[, official_anchor_label := official_label_character[
  match(patient_id, official_ids)
]]
patient_probabilities[, gdc_projection_label := as.character(gdc_label)[
  match(patient_id, rownames(gdc_input))
]]

for (class_name in ecms_levels) {
  patient_probabilities[, (paste0("official_prob_", class_name)) :=
    official_probability[match(patient_id, official_ids), class_name]]
  patient_probabilities[, (paste0("gdc_prob_", class_name)) :=
    gdc_probability[patient_id, class_name]]
}
patient_probabilities[, official_margin_custom := official_margin[
  match(patient_id, official_ids)
]]
patient_probabilities[, gdc_margin_custom := probability_margin(gdc_probability)[patient_id]]
patient_probabilities[, gdc_agrees_with_official := fifelse(
  in_official_78,
  official_anchor_label == gdc_projection_label,
  NA
)]
patient_probabilities[, resolved_ecms_label := fifelse(
  in_official_78, official_anchor_label, gdc_projection_label
)]
patient_probabilities[, label_source := fifelse(
  in_official_78,
  "locked_model_prediction_on_repository_bundled_tcga78_anchor",
  fifelse(calibration_pass,
          "gdc_tpm_extension_after_overlap_calibration_pass",
          "gdc_tpm_extension_conditional_not_primary")
)]
patient_probabilities[, extension_status := fifelse(
  in_official_78, "not_applicable_official_anchor",
  fifelse(calibration_pass, "calibrated_extension", "conditional_extension_retained")
)]
patient_probabilities[, eligible_for_primary_association :=
  in_official_78 | calibration_pass]

for (class_name in ecms_levels) {
  official_column <- paste0("official_prob_", class_name)
  gdc_column <- paste0("gdc_prob_", class_name)
  resolved_column <- paste0("resolved_prob_", class_name)
  patient_probabilities[, (resolved_column) := fifelse(
    in_official_78, get(official_column), get(gdc_column)
  )]
}
patient_probabilities[, resolved_margin_custom := fifelse(
  in_official_78, official_margin_custom, gdc_margin_custom
)]
patient_probabilities[, low_margin_custom_flag := resolved_margin_custom < 0.10]
patient_probabilities[, soft_uncertainty_custom := fifelse(
  low_margin_custom_flag,
  "higher_uncertainty_margin_below_0.10",
  "not_flagged_by_project_margin_0.10"
)]
patient_probabilities[, `:=`(
  official_anchor_definition = paste0(
    "prediction from locked rf.cl on repository-bundled tcga.val.df; ",
    "not an author hand-curated label"
  ),
  margin_definition = "project-defined top1_probability minus top2_probability",
  margin_is_official_threshold = FALSE,
  low_margin_rejects_label = FALSE,
  normalization_scope = "94-patient GDC five-layer core for reprojection",
  single_sample_classifier_claim = FALSE,
  pseudo_label_generated = FALSE
)]
fail_if(patient_probabilities[in_official_78 == TRUE,
                              any(resolved_ecms_label != official_anchor_label)],
        "作者仓库内置 78 例矩阵的锁定模型 anchor 预测未被优先保留。")
fail_if(!calibration_pass &&
          patient_probabilities[in_official_78 == FALSE,
                                any(eligible_for_primary_association)],
        "校准失败时额外 16 例被误纳入 primary。")

message("[5/7] 计算 ECMS–Factor 和 ECMS 调整后 Factor–PROGENy 关联")
factor_scores <- fread(factor_path, showProgress = FALSE)
progeny_scores <- fread(progeny_path, showProgress = FALSE)
level_factor_qc <- fread(
  level_factor_qc_path,
  colClasses = "character",
  na.strings = NULL,
  showProgress = FALSE
)
factor_columns <- grep("^Factor[0-9]+$", names(factor_scores), value = TRUE)
pathway_columns <- setdiff(names(progeny_scores), "patient_id")
fail_if(length(factor_columns) != 8L || length(pathway_columns) != 14L ||
          uniqueN(factor_scores$patient_id) != 94L ||
          uniqueN(progeny_scores$patient_id) != 94L ||
          !setequal(factor_scores$patient_id, patients) ||
          !setequal(progeny_scores$patient_id, patients),
        "MOFA factor 或 PROGENy 输入维度/患者集不符。")

required_level_factor_fields <- c(
  "factor", "view", "level_factor_soft_flag"
)
expected_level_factor_views <- c("RNA", "miRNA", "HM450", "Mutation", "CNV")
fail_if(
  !all(required_level_factor_fields %in% names(level_factor_qc)) ||
    nrow(level_factor_qc) != 40L ||
    uniqueN(level_factor_qc[, .(factor, view)]) != 40L ||
    !setequal(level_factor_qc$factor, factor_columns) ||
    !setequal(level_factor_qc$view, expected_level_factor_views),
  "冻结 level-factor QC 字段、8×5 唯一结构或 factor/view 身份不符。"
)
level_factor_qc[, level_factor_soft_flag := as_logical_strict(
  level_factor_soft_flag, "level-factor QC flag"
)]
factor_level_flags <- level_factor_qc[, .(
  factor_level_soft_flag = any(level_factor_soft_flag)
), by = factor]
setorder(factor_level_flags, factor)
fail_if(
  nrow(factor_level_flags) != length(factor_columns) ||
    !setequal(factor_level_flags$factor, factor_columns) ||
    !factor_level_flags[factor == "Factor4", factor_level_soft_flag],
  "冻结 level-factor QC 未完整映射 8 个因子，或 Factor4 未被标记。"
)
factor_level_flag_map <- setNames(
  factor_level_flags$factor_level_soft_flag,
  factor_level_flags$factor
)

state_data <- merge(
  patient_probabilities[, .(patient_id, resolved_ecms_label)],
  factor_scores, by = "patient_id", all = FALSE
)
state_data <- merge(state_data, progeny_scores, by = "patient_id", all = FALSE)
fail_if(nrow(state_data) != 94L || any(!is.finite(as.matrix(
  state_data[, c(factor_columns, pathway_columns), with = FALSE]
))), "MOFA/PROGENy 数值不完整。")

association_scopes <- list(
  official_anchor_78 = list(
    ids = official_ids,
    role = if (calibration_pass) "anchor_calibration_secondary" else "primary",
    ceiling = "T2_explanatory_same_TCGA_RNA_shared_representation"
  ),
  hybrid_94_extension = list(
    ids = patients,
    role = if (calibration_pass) "calibrated_primary_extension" else
      "conditional_sensitivity_only",
    ceiling = if (calibration_pass)
      "T2_explanatory_same_TCGA_RNA_shared_representation" else
      "T1_conditional_extension_failed_overlap_calibration"
  )
)

factor_associations <- rbindlist(lapply(names(association_scopes), function(scope_name) {
  scope <- association_scopes[[scope_name]]
  scope_data <- state_data[match(scope$ids, patient_id)]
  fail_if(anyNA(scope_data$patient_id), paste("关联 scope 患者缺失：", scope_name))
  ecms <- factor(scope_data$resolved_ecms_label, levels = ecms_levels)
  fail_if(anyNA(ecms) || any(table(ecms) == 0L),
          paste("关联 scope 缺少 ECMS 类别：", scope_name))
  rbindlist(lapply(factor_columns, function(factor_name) {
    level_flag <- unname(factor_level_flag_map[[factor_name]])
    fail_if(length(level_flag) != 1L || is.na(level_flag),
            paste("level-factor flag 无法唯一映射：", factor_name))
    y <- scope_data[[factor_name]]
    model <- lm(y ~ ecms)
    anova_table <- anova(model)
    total_ss <- sum(anova_table$`Sum Sq`)
    fail_if(!is.finite(total_ss) || total_ss <= 0,
            paste("因子总离差平方和无效：", scope_name, factor_name))
    eta_squared <- anova_table$`Sum Sq`[[1L]] / total_ss
    kruskal <- kruskal.test(y ~ ecms)
    fail_if(any(!is.finite(c(
      anova_table$`F value`[[1L]], anova_table$`Pr(>F)`[[1L]], eta_squared,
      unname(kruskal$statistic), kruskal$p.value
    ))), paste("因子–ECMS 检验产生 NA/Inf：", scope_name, factor_name))
    class_means <- vapply(ecms_levels, function(class_name) {
      mean(y[ecms == class_name])
    }, numeric(1))
    class_n <- table(ecms)
    data.table(
      analysis_scope = scope_name,
      association_role = scope$role,
      is_primary_scope = scope$role %in% c("primary", "calibrated_primary_extension"),
      n_patients = length(y),
      factor = factor_name,
      ecms1_n = as.integer(class_n[["ECMS1"]]),
      ecms2_n = as.integer(class_n[["ECMS2"]]),
      ecms3_n = as.integer(class_n[["ECMS3"]]),
      ecms4_n = as.integer(class_n[["ECMS4"]]),
      ecms1_mean = class_means[["ECMS1"]],
      ecms2_mean = class_means[["ECMS2"]],
      ecms3_mean = class_means[["ECMS3"]],
      ecms4_mean = class_means[["ECMS4"]],
      anova_f = unname(anova_table$`F value`[[1L]]),
      anova_p_value = unname(anova_table$`Pr(>F)`[[1L]]),
      eta_squared = eta_squared,
      kruskal_h = unname(kruskal$statistic),
      kruskal_p_value = kruskal$p.value,
      extension_calibration_pass = calibration_pass,
      shared_rna_representation = TRUE,
      independent_validation = FALSE,
      factor_level_soft_flag = level_flag,
      evidence_ceiling = if (level_flag) {
        "T0_technical_background_level_factor"
      } else {
        scope$ceiling
      },
      interpretation = if (level_flag) {
        "level-factor soft QC; technical/background only, not a biological ECMS axis"
      } else {
        "same-TCGA explanatory association; not a new subtype or causal mechanism"
      }
    )
  }))
}), use.names = TRUE)
factor_associations[, anova_q_value := p.adjust(anova_p_value, method = "BH"),
                    by = analysis_scope]
factor_associations[, kruskal_q_value := p.adjust(kruskal_p_value, method = "BH"),
                    by = analysis_scope]
factor_numeric_fields <- c(
  "anova_f", "anova_p_value", "eta_squared", "kruskal_h",
  "kruskal_p_value", "anova_q_value", "kruskal_q_value"
)
fail_if(any(!is.finite(as.matrix(
  factor_associations[, ..factor_numeric_fields]
))), "ECMS–Factor 正式数值字段含 NA/Inf。")
setorder(factor_associations, analysis_scope, anova_q_value, -eta_squared)

adjusted_associations <- rbindlist(lapply(names(association_scopes), function(scope_name) {
  scope <- association_scopes[[scope_name]]
  scope_data <- state_data[match(scope$ids, patient_id)]
  ecms <- factor(scope_data$resolved_ecms_label, levels = ecms_levels)
  rbindlist(lapply(pathway_columns, function(pathway_name) {
    y <- safe_scale(scope_data[[pathway_name]], paste(scope_name, pathway_name))
    rbindlist(lapply(factor_columns, function(factor_name) {
      level_flag <- unname(factor_level_flag_map[[factor_name]])
      fail_if(length(level_flag) != 1L || is.na(level_flag),
              paste("level-factor flag 无法唯一映射：", factor_name))
      x <- safe_scale(scope_data[[factor_name]], paste(scope_name, factor_name))
      model_data <- data.frame(y = y, x = x, ecms = ecms)
      reduced <- lm(y ~ ecms, data = model_data)
      full <- lm(y ~ ecms + x, data = model_data)
      fail_if(
        reduced$rank < ncol(model.matrix(reduced)) ||
          full$rank < ncol(model.matrix(full)),
        paste("增量关联模型奇异/秩不足：", scope_name,
              pathway_name, factor_name)
      )
      comparison <- anova(reduced, full)
      coefficient <- summary(full)$coefficients["x", ]
      reduced_r2 <- summary(reduced)$r.squared
      full_r2 <- summary(full)$r.squared
      partial_r2 <- if (reduced_r2 < 1) {
        (full_r2 - reduced_r2) / (1 - reduced_r2)
      } else {
        NA_real_
      }
      critical_values <- c(
        coefficient, reduced_r2, full_r2, full_r2 - reduced_r2, partial_r2,
        comparison$F[[2L]], comparison$`Pr(>F)`[[2L]]
      )
      fail_if(any(!is.finite(critical_values)),
              paste("增量 Factor–PROGENy 数值含 NA/Inf：", scope_name,
                    pathway_name, factor_name))
      data.table(
        analysis_scope = scope_name,
        association_role = scope$role,
        is_primary_scope = scope$role %in% c("primary", "calibrated_primary_extension"),
        n_patients = nrow(model_data),
        pathway = pathway_name,
        factor = factor_name,
        reduced_model = "standardized_PROGENy ~ ECMS",
        full_model = "standardized_PROGENy ~ ECMS + standardized_MOFA_factor",
        factor_beta_standardized = unname(coefficient[["Estimate"]]),
        factor_standard_error = unname(coefficient[["Std. Error"]]),
        factor_t_value = unname(coefficient[["t value"]]),
        factor_coefficient_p_value = unname(coefficient[["Pr(>|t|)"]]),
        reduced_r_squared = reduced_r2,
        full_r_squared = full_r2,
        delta_r_squared = full_r2 - reduced_r2,
        partial_r_squared = partial_r2,
        incremental_f = comparison$F[[2L]],
        incremental_p_value = comparison$`Pr(>F)`[[2L]],
        extension_calibration_pass = calibration_pass,
        shared_rna_representation = TRUE,
        independent_validation = FALSE,
        factor_level_soft_flag = level_flag,
        evidence_ceiling = if (level_flag) {
          "T0_technical_background_level_factor"
        } else {
          scope$ceiling
        },
        interpretation = if (level_flag) {
          "level-factor soft QC; technical/background only, not a biological ECMS axis"
        } else {
          "incremental same-TCGA association after categorical ECMS adjustment; not causal"
        }
      )
    }))
  }))
}), use.names = TRUE)
adjusted_associations[, incremental_q_value := p.adjust(
  incremental_p_value, method = "BH"
), by = analysis_scope]
adjusted_associations[, factor_coefficient_q_value := p.adjust(
  factor_coefficient_p_value, method = "BH"
), by = analysis_scope]
adjusted_numeric_fields <- c(
  "factor_beta_standardized", "factor_standard_error", "factor_t_value",
  "factor_coefficient_p_value", "reduced_r_squared", "full_r_squared",
  "delta_r_squared", "partial_r_squared", "incremental_f",
  "incremental_p_value", "incremental_q_value", "factor_coefficient_q_value"
)
fail_if(any(!is.finite(as.matrix(
  adjusted_associations[, ..adjusted_numeric_fields]
))), "ECMS 调整后 Factor–PROGENy 正式数值字段含 NA/Inf。")
setorder(adjusted_associations, analysis_scope, incremental_q_value, -partial_r_squared)

verify_level_factor_propagation <- function(object, object_name) {
  observed <- unique(object[, .(factor, factor_level_soft_flag)])
  setorder(observed, factor)
  expected <- copy(factor_level_flags)
  setorder(expected, factor)
  fail_if(
    !identical(observed$factor, expected$factor) ||
      !identical(
        observed$factor_level_soft_flag,
        expected$factor_level_soft_flag
      ),
    paste(object_name, "的 level-factor flag 与冻结 QC 不一致。")
  )
  fail_if(
    any(object$factor_level_soft_flag &
          object$evidence_ceiling !=
            "T0_technical_background_level_factor"),
    paste(object_name, "存在触发 level-factor flag 但未降为 T0 的行。")
  )
}
verify_level_factor_propagation(factor_associations, "ECMS–Factor 表")
verify_level_factor_propagation(
  adjusted_associations, "ECMS 调整后 Factor–PROGENy 表"
)
factor_level_t0_rows <- factor_associations[
  factor_level_soft_flag == TRUE &
    evidence_ceiling == "T0_technical_background_level_factor", .N
]
adjusted_level_t0_rows <- adjusted_associations[
  factor_level_soft_flag == TRUE &
    evidence_ceiling == "T0_technical_background_level_factor", .N
]
fail_if(
  factor_associations[factor == "Factor4", .N] != 2L ||
    factor_associations[factor == "Factor4", sum(factor_level_soft_flag)] != 2L ||
    adjusted_associations[factor == "Factor4", .N] != 28L ||
    adjusted_associations[factor == "Factor4", sum(factor_level_soft_flag)] != 28L,
  "Factor4 必须在两个 scope 的 2 条 factor 行和 28 条 adjusted pathway 行中全部降为 T0。"
)

message("[6/7] 生成输入/模型 QA、边界摘要与 stage manifest")
qa <- data.table(
  check_id = c(
    "locked_commit", "model_bytes", "model_manifest_full_sha",
    "upstream_artifact_manifests", "model_object_set", "randomforest_structure",
    "feature_count",
    "feature_unique", "feature_list_sha256", "builtin_tcga_shape",
    "builtin_tcga_gene_scaling", "builtin_prediction_counts",
    "five_layer_core_count", "official_id_mapping",
    "gdc_tpm_assay", "ensembl_strip_unique", "feature_coverage",
    "feature_nonzero_variance", "gdc_projection_shape", "probability_row_sums",
    "official_anchor_precedence", "extension_gate", "margin_boundary",
    "pseudo_label_prohibition", "clinical_classifier_boundary"
  ),
  category = c(
    rep("model_provenance", 6), rep("feature_identity", 3),
    rep("repository_anchor_validation", 3), rep("patient_identity", 2),
    rep("gdc_input", 6), "label_resolution", "calibration", "uncertainty",
    "evidence_boundary", "evidence_boundary"
  ),
  hard_gate = c(rep(TRUE, 21), FALSE, TRUE, TRUE, TRUE),
  expected = c(
    locked_commit, paste0(locked_model_size, ";", locked_model_sha256),
    "all manifest-listed files match size and sha256",
    "multiassay and heterogeneity manifests cover locked inputs; all listed artifacts verified",
    paste(sort(expected_objects), collapse = ";"),
    "randomForest;ntree=500;mtry=17;classes=ECMS1-4;cutoff=0.25x4",
    "314", "314 unique",
    locked_feature_sha256, "78x314", "max_abs_mean<=1e-6;max_abs_sd_minus_1<=1e-6",
    "23/34/7/14", "94 unique", "78/78",
    "tpm", "all stripped IDs unique", "314/314", "314/314 nonzero variance",
    "94x314", "all row sums=1", "official labels retained for all 78",
    "four prelocked metrics", "margin<0.10 soft flag only",
    "no cluster/pathway pseudo-label", "no clinical classifier claim"
  ),
  observed = c(
    locked_commit, paste0(
      sprintf("%.0f", as.numeric(file_info(model_path)$size)),
      ";", sha256_file(model_path)
    ),
    paste0(nrow(model_manifest), " files verified"),
    paste0(nrow(multiassay_upstream_manifest), "+",
           nrow(heterogeneity_upstream_manifest), " artifacts verified"),
    paste(sort(loaded_objects), collapse = ";"),
    paste0(class(rf_classifier)[1], ";ntree=", rf_classifier$ntree,
           ";mtry=", rf_classifier$mtry, ";classes=",
           paste(rf_classifier$classes, collapse = ","), ";cutoff=",
           paste(rf_classifier$forest$cutoff, collapse = ",")),
    length(gene_features), uniqueN(gene_features), feature_list_sha256(gene_features),
    paste(dim(tcga_official_matrix), collapse = "x"),
    paste0("max_abs_mean=", max(abs(tcga_official_column_means)),
           ";max_abs_sd_minus_1=", max(abs(tcga_official_column_sds - 1))),
    paste(as.integer(official_counts), collapse = "/"), length(patients),
    paste0(sum(official_ids %in% patients), "/78"), "tpm",
    paste0(uniqueN(stripped_gene_ids), "/", length(stripped_gene_ids)),
    paste0(sum(!is.na(feature_index)), "/314"),
    paste0(sum(feature_sd > 0), "/314"), paste(dim(gdc_input), collapse = "x"),
    paste(range(rowSums(gdc_probability)), collapse = ";"),
    paste0(sum(patient_probabilities[
      in_official_78 == TRUE,
      resolved_ecms_label == official_anchor_label
    ]), "/78"),
    paste0(sum(gate_pass), "/4;overall=", calibration_pass),
    paste0(sum(patient_probabilities$low_margin_custom_flag),
           " flagged;0 rejected"),
    paste0(sum(patient_probabilities$pseudo_label_generated), " pseudo labels"),
    paste0(sum(patient_probabilities$single_sample_classifier_claim),
           " single-sample/clinical claims")
  ),
  pass = c(rep(TRUE, 21), calibration_pass, TRUE, TRUE, TRUE),
  notes = c(
    "official GitHub commit locked by script22", "model byte identity",
    paste0("manifest_sha256=", sha256_file(model_manifest_path)),
    paste0("multiassay_manifest_sha256=", sha256_file(multiassay_manifest_path),
           ";heterogeneity_manifest_sha256=", sha256_file(heterogeneity_manifest_path)),
    "loaded in isolated environment", "serialized fitted model application only",
    "model feature count", "no duplicated model feature", "UTF-8 lines with final newline",
    "repository-bundled matrix", "repository matrix is already gene-wise scaled",
    "locked-model prediction reproduction; not author hand-curated labels",
    "analysis set from formal project result", "all model rownames map to GDC patients",
    "uses GDC TPM, not Xena log2(count+1)", "no silent duplicate collapse",
    "full feature coverage", "checked after log2(TPM+1)",
    "within-94 cohort gene-wise scale", "four-class probabilities",
    "locked-model predictions on repository-bundled 78 matrix override GDC reprojection on overlap",
    "failure keeps extra 16 labels conditional and eligible=FALSE",
    "project-defined; not an author threshold; never rejects official labels",
    "all labels come from locked rf.cl", "batch research benchmark only"
  )
)
flagged_factors <- factor_level_flags[
  factor_level_soft_flag == TRUE, paste(sort(factor), collapse = ",")
]
qa <- rbind(
  qa,
  data.table(
    check_id = c(
      "level_factor_qc_identity",
      "level_factor_flag_propagation",
      "level_factor_evidence_ceiling"
    ),
    category = rep("level_factor_boundary", 3L),
    hard_gate = rep(TRUE, 3L),
    expected = c(
      "40 unique factor-view rows;Factor4 included among flagged factors",
      "factor and adjusted tables exactly match frozen factor flags",
      "all flagged factor rows are T0 technical/background"
    ),
    observed = c(
      paste0("40 unique factor-view rows;flagged factors=", flagged_factors),
      "both formal tables exactly match frozen factor flags",
      paste0(
        "factor table T0 rows=", factor_level_t0_rows,
        ";adjusted table T0 rows=", adjusted_level_t0_rows
      )
    ),
    pass = rep(TRUE, 3L),
    notes = c(
      paste0("level_factor_qc_sha256=", sha256_file(level_factor_qc_path)),
      "flags propagated before formal sorting and stage publication",
      "Factor4 is not interpreted as a biological heterogeneity axis"
    )
  ),
  use.names = TRUE
)
fail_if(anyDuplicated(qa$check_id), "ECMS 投影 QA check_id 不唯一。")
fail_if(any(qa$hard_gate & !qa$pass), "ECMS 投影硬 QA 未全部通过。")

primary_n <- sum(patient_probabilities$eligible_for_primary_association)
additional_status <- if (calibration_pass) {
  "额外 16 例已进入 primary 扩展"
} else {
  "额外 16 例只条件保留，未进入 primary"
}
summary_lines <- c(
  "# TCGA-ESCC ECMS 外部基准投影摘要",
  "",
  "## 模型与输入",
  "",
  paste0("- 作者 GitHub commit：`", locked_commit, "`。"),
  paste0("- `ECMS.model.rdata`：", locked_model_size, " B，SHA256 `",
         locked_model_sha256, "`。"),
  "- 输入为 94 例 five-layer core 的 GDC TPM；Ensembl ID 去版本后 314/314 特征唯一覆盖且均为非零方差。",
  "- 预处理复现作者 README：`log2(TPM+1)` 后在待投影 94 例队列内逐基因 Z 标准化。",
  "- 锁定模型对作者仓库内置 `tcga.val.df` 的预测计数复现为 ECMS1/2/3/4 = 23/34/7/14；这是模型预测 anchor，不是作者手工标签。",
  "",
  "## 标签解析与校准",
  "",
  paste0("- 78 例重叠患者始终采用锁定模型对作者仓库内置矩阵的预测 anchor；",
         additional_status, "。"),
  paste0("- 校准：agreement=", sprintf("%.3f", agreement),
         "，kappa=", sprintf("%.3f", kappa), "，ARI=", sprintf("%.3f", ari),
         "，四类概率向量 Spearman 中位数=",
         sprintf("%.3f", median_probability_spearman),
         "；整体门禁=`", calibration_pass, "`。"),
  paste0("- primary 关联可纳入患者数：", primary_n, "。"),
  "- top1–top2 概率 margin <0.10 只是项目自定义软不确定性标记；不是作者阈值，不拒绝仓库 78 例 anchor 预测。",
  "",
  "## 关联边界",
  "",
  paste0(
    "- 冻结 level-factor QC 标记：", flagged_factors,
    "；这些因子及其 ECMS 调整后通路行统一只作 T0 技术/背景项。"
  ),
  "- ECMS–MOFA 关联和控制 ECMS 后的 Factor–PROGENy 增量关联均在同一 TCGA 患者集中完成。",
  "- 未触发 level-factor flag 的 ECMS/MOFA/PROGENy 关联最高只是 T2 解释性关联；Factor4 不作为生物学异质性轴。",
  "- ECMS 和 PROGENy 直接依赖 RNA，MOFA 也含 RNA 视图；共享表示不能计作独立验证或因果调控轴。",
  "- 本投影是公共数据研究基准，不宣称单样本、临床诊断、预后或疗效分类器性能。",
  "- 没有从 MOFA 聚类数、通路相似性或手工规则生成 ECMS 伪标签。"
)

atomic_fwrite(patient_probabilities,
              file.path(stage_dir, "tcga_escc_ecms_patient_probabilities.tsv"))
atomic_fwrite(calibration,
              file.path(stage_dir, "tcga_escc_ecms_projection_calibration.tsv"))
atomic_fwrite(qa, file.path(stage_dir, "tcga_escc_ecms_projection_qa.tsv"))
atomic_fwrite(factor_associations,
              file.path(stage_dir, "tcga_escc_ecms_factor_associations.tsv"))
atomic_fwrite(
  adjusted_associations,
  file.path(stage_dir, "tcga_escc_ecms_adjusted_factor_progeny_associations.tsv")
)
atomic_write_lines(
  summary_lines, file.path(stage_dir, "tcga_escc_ecms_projection_summary.md")
)

input_hashes <- paste(
  c(
    paste0("model=", locked_model_sha256),
    paste0("model_manifest=", sha256_file(model_manifest_path)),
    paste0("multiassay_manifest=", sha256_file(multiassay_manifest_path)),
    paste0("heterogeneity_manifest=", sha256_file(heterogeneity_manifest_path)),
    paste0("multiassay=", sha256_file(mae_path)),
    paste0("analysis_sets=", sha256_file(analysis_sets_path)),
    paste0("factor_scores=", sha256_file(factor_path)),
    paste0("progeny_scores=", sha256_file(progeny_path)),
    paste0("level_factor_qc=", sha256_file(level_factor_qc_path))
  ),
  collapse = ";"
)
artifact_paths <- file.path(stage_dir, formal_filenames)
fail_if(any(!file_exists(artifact_paths)), "ECMS stage 正式 artifact 缺失。")
artifact_manifest <- data.table(
  artifact = formal_filenames,
  relative_path = file.path("results", formal_filenames),
  file_size_bytes = as.numeric(file_info(artifact_paths)$size),
  sha256 = vapply(artifact_paths, sha256_file, character(1)),
  generated_date = as.character(Sys.Date()),
  generation_script = "scripts/23_project_tcga_escc_ecms.R",
  execution_script_sha256 = execution_script_sha256,
  input_sha256 = input_hashes,
  model_commit = locked_commit,
  extension_calibration_pass = calibration_pass,
  primary_patient_count = primary_n,
  status = "verified"
)
atomic_fwrite(artifact_manifest, file.path(stage_dir, manifest_filename))

for (filename in formal_filenames) {
  path <- file.path(stage_dir, filename)
  if (grepl("[.]tsv$", filename)) {
    reread <- fread(path, colClasses = "character", na.strings = NULL)
    fail_if(!nrow(reread), paste("正式 TSV 为空：", filename))
  } else {
    fail_if(as.numeric(file_info(path)$size) <= 0,
            paste("正式摘要为空：", filename))
  }
}

message("[7/7] 原子发布正式结果，manifest 最后发布")
for (filename in formal_filenames) {
  atomic_publish_file(
    file.path(stage_dir, filename),
    file.path(results_dir, filename)
  )
}
for (i in seq_len(nrow(artifact_manifest))) {
  published <- file.path(project_root, artifact_manifest$relative_path[[i]])
  fail_if(!file_exists(published) ||
            as.numeric(file_info(published)$size) != artifact_manifest$file_size_bytes[[i]] ||
            sha256_file(published) != artifact_manifest$sha256[[i]],
          paste("发布 artifact 回读失败：", published))
}
atomic_publish_file(
  file.path(stage_dir, manifest_filename),
  file.path(results_dir, manifest_filename)
)
published_manifest <- fread(
  file.path(results_dir, manifest_filename),
  colClasses = "character", na.strings = NULL
)
fail_if(nrow(published_manifest) != nrow(artifact_manifest) ||
          !identical(published_manifest$relative_path, artifact_manifest$relative_path) ||
          any(published_manifest$sha256 != artifact_manifest$sha256),
        "ECMS 发布 manifest 与 stage 冻结版不一致。")

if (dir_exists(stage_dir)) dir_delete(stage_dir)
message(
  "完成：锁定模型对作者仓库内置 78 例矩阵的预测 anchor 已复现；GDC TPM 重投影校准=",
  calibration_pass, "；primary n=", primary_n,
  "；未触发 level-factor flag 的结论上限为同 TCGA RNA 共享表示的 T2 解释性关联；",
  flagged_factors, " 仅作 T0 技术/背景项。"
)
