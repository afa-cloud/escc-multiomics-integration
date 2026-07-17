#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(fs)
  library(glmnet)
  library(logistf)
  library(MOFA2)
  library(randomForest)
  library(readxl)
  library(R.utils)
  library(survival)
})

options(stringsAsFactors = FALSE)

# ESCC 连续状态第一阶段外部验证。
#
# 证据边界：
# 1) 以锁定 ECMS 资源中的 TCGA 78×314 矩阵训练 Factor1/Factor3 RNA
#    ridge proxy；它不是完整多组学因子的外部重建。
# 2) 模型与系数在读取 GSE53622/GSE53624 生存结局及 GSE45670 疗效标签前
#    冻结并写入 stage RDS；外部结局不得反向选择特征、方向、lambda 或截点。
# 3) GSE53622/GSE53624 分别重做队列内逐基因 Z 标准化；GSE53625 是两者
#    的精确超级系列，绝不作为第三个独立队列。
# 4) 生存是次要端点；阴性、PH 偏离或队列异质性保留，不作为 proxy
#    可计算性的机械失败门禁。
# 5) GSE45670 只作小样本探索性 pCR 关联，不宣称疗效预测器。

formal_filenames <- c(
  "escc_external_state_proxy_models.rds",
  "escc_external_state_proxy_definition.tsv",
  "escc_external_state_internal_cv.tsv",
  "escc_external_state_oof_predictions.tsv",
  "escc_external_state_patient_scores.tsv",
  "escc_external_state_survival_associations.tsv",
  "escc_external_state_ecms_increment.tsv",
  "escc_external_state_response_associations.tsv",
  "escc_external_validation_decision.tsv",
  "escc_external_validation_summary.md"
)
manifest_filename <- "escc_external_validation_artifact_manifest.tsv"

args <- commandArgs(trailingOnly = TRUE)
fields_only <- "--fields-only" %in% args
unknown_args <- setdiff(args, "--fields-only")
if (length(unknown_args)) {
  stop("未知参数：", paste(unknown_args, collapse = ", "), call. = FALSE)
}
if (fields_only) {
  cat(paste(c(formal_filenames, manifest_filename), collapse = "\n"), "\n")
  quit(save = "no", status = 0L)
}

fail_if <- function(condition, message) {
  if (length(condition) != 1L || is.na(condition) || isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
fail_if(length(script_argument) != 1L,
        "无法从 --file 唯一定位 scripts/25_validate_external_continuous_states.R。")
script_path <- normalizePath(
  sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE
)
project_root <- normalizePath(
  file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE
)
results_dir <- file.path(project_root, "results")
work_intermediate_dir <- file.path(project_root, "_work", "intermediate")
data_root <- Sys.getenv("RESEARCH_DATA_ROOT", path.expand("~/ResearchDataHub"))
catalog_path <- file.path(data_root, "CATALOG.tsv")
execution_script_sha256 <- digest(
  script_path, algo = "sha256", file = TRUE, serialize = FALSE
)

required_project_inputs <- file.path(project_root, c(
  "PROJECT_INDEX.md",
  "data/datasets.tsv",
  "results/tcga_escc_mofa_model.rds",
  "results/tcga_escc_mofa_factor_scores.tsv",
  "results/tcga_escc_heterogeneity_artifact_manifest.tsv"
))
fail_if(any(!file_exists(required_project_inputs)), paste(
  "缺少外部连续状态验证输入：",
  paste(required_project_inputs[!file_exists(required_project_inputs)], collapse = ";")
))
fail_if(!dir_exists(data_root) || !file_exists(catalog_path),
        paste("ResearchDataHub 或 CATALOG 不可读：", data_root))
fail_if(!dir_exists(results_dir), "results/ 不可读。")
dir_create(work_intermediate_dir, recurse = TRUE)
stage_dir <- tempfile(
  pattern = ".escc_external_continuous_states_",
  tmpdir = work_intermediate_dir
)
dir_create(stage_dir)
cleanup_active_stage <- function() {
  if (exists("stage_dir", inherits = TRUE) && dir_exists(stage_dir)) {
    dir_delete(stage_dir)
  }
  invisible(NULL)
}
previous_error_handler <- getOption("error")
options(error = function() {
  try(cleanup_active_stage(), silent = TRUE)
  if (is.function(previous_error_handler)) previous_error_handler()
})

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

feature_list_sha256 <- function(features) {
  payload <- charToRaw(paste0(paste(features, collapse = "\n"), "\n"))
  digest(payload, algo = "sha256", serialize = FALSE)
}

stage_fwrite <- function(object, filename) {
  path <- file.path(stage_dir, filename)
  fwrite(object, path, sep = "\t", quote = FALSE, na = "", logical01 = FALSE)
  reread <- fread(path, colClasses = "character", na.strings = NULL,
                  showProgress = FALSE)
  fail_if(nrow(reread) != nrow(object) || !identical(names(reread), names(object)),
          paste("stage TSV 回读失败：", filename))
  path
}

stage_write_lines <- function(lines, filename) {
  path <- file.path(stage_dir, filename)
  writeLines(lines, path, useBytes = TRUE)
  fail_if(!file_exists(path) || as.numeric(file_info(path)$size) <= 0,
          paste("stage 文本写入失败：", filename))
  path
}

atomic_publish_file <- function(source, destination) {
  fail_if(!file_exists(source), paste("发布源文件缺失：", source))
  dir_create(dirname(destination), recurse = TRUE)
  temp_path <- tempfile(
    pattern = paste0(".", basename(destination), ".publishing."),
    tmpdir = dirname(destination)
  )
  on.exit(if (file_exists(temp_path)) file_delete(temp_path), add = TRUE)
  copied <- file.copy(source, temp_path, overwrite = FALSE, copy.mode = TRUE)
  fail_if(!copied || !file_exists(temp_path) ||
            as.numeric(file_info(source)$size) != as.numeric(file_info(temp_path)$size) ||
            sha256_file(source) != sha256_file(temp_path),
          paste("发布临时复制或 SHA256 失败：", destination))
  renamed <- file.rename(temp_path, destination)
  fail_if(!renamed || !file_exists(destination) ||
            sha256_file(source) != sha256_file(destination),
          paste("原子发布失败：", destination))
  invisible(destination)
}

verify_hub_manifest <- function(dataset_root, manifest_path, dataset_key) {
  fail_if(!dir_exists(dataset_root) || !file_exists(manifest_path),
          paste("DataHub 数据集目录或 manifest 缺失：", dataset_key))
  manifest <- fread(manifest_path, colClasses = "character", na.strings = NULL,
                    showProgress = FALSE)
  required <- c("relative_path", "size_bytes", "sha256", "file_status")
  fail_if(!all(required %in% names(manifest)) || !nrow(manifest) ||
            anyDuplicated(manifest$relative_path),
          paste("DataHub manifest 字段、行数或唯一性失败：", dataset_key))
  paths <- file.path(dataset_root, manifest$relative_path)
  fail_if(any(!file_exists(paths)), paste("DataHub manifest 文件缺失：", dataset_key))
  observed_size <- as.character(as.numeric(file_info(paths)$size))
  observed_sha <- vapply(paths, sha256_file, character(1))
  fail_if(any(observed_size != manifest$size_bytes) ||
            any(observed_sha != manifest$sha256) ||
            any(!manifest$file_status %chin% c(
              "verified", "generated_verified", "verified_locked_git_blob",
              "verified_locked_commit_archive", "verified_api_metadata",
              "generated_verified_git_blob"
            )),
          paste("DataHub manifest 大小、SHA256 或状态失败：", dataset_key))
  manifest
}

verify_hub_dataset <- function(catalog, dataset_key) {
  key_value <- dataset_key
  row <- catalog[dataset_key == key_value]
  fail_if(nrow(row) != 1L || row$status != "verified",
          paste("CATALOG 缺少唯一 verified 数据集：", dataset_key))
  root <- row$local_path[[1L]]
  manifest_path <- row$manifest_path[[1L]]
  dataset_md <- file.path(root, "DATASET.md")
  fail_if(!file_exists(dataset_md) || as.numeric(file_info(dataset_md)$size) <= 0,
          paste("DATASET.md 缺失或为空：", dataset_key))
  manifest <- verify_hub_manifest(root, manifest_path, dataset_key)
  list(
    key = dataset_key,
    catalog_row = row,
    root = root,
    manifest_path = manifest_path,
    manifest = manifest,
    manifest_sha256 = sha256_file(manifest_path),
    dataset_md_sha256 = sha256_file(dataset_md)
  )
}

verify_project_manifest <- function(manifest_path, required_relatives) {
  manifest <- fread(manifest_path, colClasses = "character", na.strings = NULL,
                    showProgress = FALSE)
  required <- c("relative_path", "file_size_bytes", "sha256", "status",
                "generation_script")
  fail_if(!all(required %in% names(manifest)) || !nrow(manifest) ||
            anyDuplicated(manifest$relative_path) ||
            !all(required_relatives %in% manifest$relative_path),
          "异质性 artifact manifest 字段、唯一性或覆盖失败。")
  paths <- file.path(project_root, manifest$relative_path)
  fail_if(any(!file_exists(paths)), "异质性 artifact manifest 登记文件缺失。")
  observed_size <- as.character(as.numeric(file_info(paths)$size))
  observed_sha <- vapply(paths, sha256_file, character(1))
  fail_if(any(observed_size != manifest$file_size_bytes) ||
            any(observed_sha != manifest$sha256) || any(manifest$status != "verified") ||
            any(manifest$generation_script !=
                  "scripts/14_analyze_tcga_escc_heterogeneity.R"),
          "异质性 artifact manifest 大小、SHA256、状态或生成脚本失败。")
  manifest
}

hub_file <- function(dataset, relative_path) {
  fail_if(!relative_path %in% dataset$manifest$relative_path,
          paste("DataHub manifest 未登记：", dataset$key, relative_path))
  path <- file.path(dataset$root, relative_path)
  fail_if(!file_exists(path), paste("DataHub 文件缺失：", path))
  path
}

decision_rows <- list()
decision_counter <- 0L
add_decision <- function(
    decision_id, domain, factor = NA_character_, hard_gate = FALSE,
    metric, observed, rule, status, counts_toward_external_support,
    interpretation, boundary) {
  decision_counter <<- decision_counter + 1L
  decision_rows[[decision_counter]] <<- data.table(
    decision_id = decision_id,
    decision_domain = domain,
    factor = factor,
    hard_gate = hard_gate,
    metric = metric,
    observed = as.character(observed),
    threshold_or_rule = rule,
    status = status,
    counts_toward_external_support = counts_toward_external_support,
    interpretation = interpretation,
    boundary = boundary
  )
}

message("[1/10] 核验 ResearchDataHub CATALOG、DATASET、MANIFEST 与上游 artifact")
catalog <- fread(catalog_path, colClasses = "character", na.strings = NULL,
                 showProgress = FALSE)
dataset_keys <- c(
  ecms = paste0(
    "GITHUB_CITYUHK_COMPUTATIONAL_BIOLOGY_ESCC_CMS_",
    "8a49db65f86f2c1e2a50d9bc4cb742ccc868e4c6"
  ),
  gse53622 = "GEO_GSE53622_retrieved_20260628",
  gse53624 = "GEO_GSE53624_retrieved_20260628",
  gse53625 = "GEO_GSE53625_retrieved_20260625",
  gse45670 = "GEO_GSE45670_retrieved_20260625"
)
hub_datasets <- lapply(dataset_keys, function(key) verify_hub_dataset(catalog, key))
for (name in names(hub_datasets)) {
  add_decision(
    paste0("INPUT_HUB_", toupper(name)), "input_integrity", hard_gate = TRUE,
    metric = "CATALOG/DATASET/MANIFEST full SHA256",
    observed = paste0("verified;manifest_rows=", nrow(hub_datasets[[name]]$manifest)),
    rule = "unique CATALOG status=verified; DATASET.md present; all manifest bytes match",
    status = "PASS", counts_toward_external_support = FALSE,
    interpretation = "输入来源与本地字节完整性通过。",
    boundary = "完整性通过不等于生物学验证。"
  )
}

heterogeneity_manifest_path <- file.path(
  results_dir, "tcga_escc_heterogeneity_artifact_manifest.tsv"
)
heterogeneity_manifest <- verify_project_manifest(
  heterogeneity_manifest_path,
  file.path("results", c(
    "tcga_escc_mofa_model.rds", "tcga_escc_mofa_factor_scores.tsv"
  ))
)
add_decision(
  "INPUT_TCGA_MOFA_ARTIFACTS", "input_integrity", hard_gate = TRUE,
  metric = "heterogeneity artifact manifest",
  observed = paste0(nrow(heterogeneity_manifest), " artifacts verified"),
  rule = "factor scores and MOFA model listed; all artifact bytes/status verified",
  status = "PASS", counts_toward_external_support = FALSE,
  interpretation = "TCGA proxy target与敏感性权重来源已冻结。",
  boundary = "该模型仍来自同一TCGA发现队列。"
)

locked_model_size <- 1608173
locked_model_sha256 <-
  "4c05645d712ec630f8ed717cba347dced4228a999c1b479ebd760bf1436e1ee4"
locked_feature_sha256 <-
  "9e81299a2e93f75f4e8375f852ca2bdd156126271f628d3c118b015d48e4120d"
locked_commit <- "8a49db65f86f2c1e2a50d9bc4cb742ccc868e4c6"
model_path <- hub_file(hub_datasets$ecms, "00_source/ECMS.model.rdata")
fail_if(as.numeric(file_info(model_path)$size) != locked_model_size ||
          sha256_file(model_path) != locked_model_sha256,
        "ECMS.model.rdata 锁定大小或 SHA256 失败。")
add_decision(
  "INPUT_ECMS_MODEL_SHA", "input_integrity", hard_gate = TRUE,
  metric = "ECMS.model.rdata SHA256",
  observed = sha256_file(model_path), rule = locked_model_sha256,
  status = "PASS", counts_toward_external_support = FALSE,
  interpretation = "使用固定commit的官方模型字节。",
  boundary = "仓库无LICENSE；不把原始314基因矩阵复制进项目。"
)

message("[2/10] 验证锁定 ECMS 对象并构建 TCGA 78×314 proxy 训练面")
model_environment <- new.env(parent = emptyenv())
loaded_objects <- load(model_path, envir = model_environment)
expected_objects <- c(
  "rf.cl", "gse53625.val.df", "gse45670.val.df", "tcga.val.df",
  "gene.features"
)
fail_if(!setequal(loaded_objects, expected_objects),
        "ECMS.model.rdata 对象集与锁定版本不一致。")
rf_classifier <- model_environment$rf.cl
gene_features <- as.character(model_environment$gene.features)
tcga_matrix <- model_environment$tcga.val.df
gse53625_matrix_original <- model_environment$gse53625.val.df
gse45670_matrix_original <- model_environment$gse45670.val.df
ecms_levels <- paste0("ECMS", 1:4)
fail_if(!inherits(rf_classifier, "randomForest") || rf_classifier$ntree != 500L ||
          rf_classifier$mtry != 17L ||
          !identical(rf_classifier$classes, ecms_levels),
        "ECMS randomForest 结构失败。")

deterministic_rf_projection <- function(classifier, matrix, class_levels) {
  probability <- stats::predict(classifier, matrix, type = "prob")
  fail_if(!identical(colnames(probability), class_levels) ||
            !identical(rownames(probability), rownames(matrix)) ||
            any(!is.finite(probability)) ||
            any(abs(rowSums(probability) - 1) > 1e-12),
          "ECMS randomForest vote/probability矩阵失败。")
  maximum <- apply(probability, 1L, max)
  at_maximum <- abs(probability - maximum) <= 1e-12
  tie_count <- rowSums(at_maximum)
  selected_index <- apply(probability, 1L, which.max)
  selected_label <- class_levels[selected_index]
  tied_classes <- vapply(seq_len(nrow(probability)), function(index) {
    paste(class_levels[at_maximum[index, ]], collapse = ",")
  }, character(1))
  output <- data.table(
    gsm = rownames(probability),
    ecms_label = selected_label,
    ecms_max_vote_fraction = maximum,
    ecms_tie_count_at_max = as.integer(tie_count),
    ecms_tie_at_max = tie_count > 1L,
    ecms_tied_classes = tied_classes,
    ecms_label_rule = paste0(
      "argmax_randomForest_type_prob; exact ties resolved by first fixed class ",
      "level ECMS1<ECMS2<ECMS3<ECMS4; outcome-blinded"
    )
  )
  for (class_name in class_levels) {
    output[, (paste0("ecms_vote_", class_name)) := probability[, class_name]]
  }
  output
}
fail_if(length(gene_features) != 314L || anyDuplicated(gene_features) ||
          feature_list_sha256(gene_features) != locked_feature_sha256,
        "ECMS 314特征身份或 SHA256 失败。")
matrix_contract <- list(
  tcga.val.df = c(78L, 314L),
  gse53625.val.df = c(179L, 314L),
  gse45670.val.df = c(28L, 314L)
)
for (name in names(matrix_contract)) {
  x <- model_environment[[name]]
  fail_if(!is.matrix(x) || !identical(dim(x), matrix_contract[[name]]) ||
            !identical(colnames(x), gene_features) || anyDuplicated(rownames(x)) ||
            any(!is.finite(x)),
          paste("锁定模型矩阵维度、特征、行名或数值失败：", name))
}
fail_if(max(abs(colMeans(tcga_matrix))) > 1e-6 ||
          max(abs(apply(tcga_matrix, 2L, sd) - 1)) > 1e-6 ||
          max(abs(colMeans(gse53625_matrix_original))) > 1e-6 ||
          max(abs(apply(gse53625_matrix_original, 2L, sd) - 1)) > 1e-6 ||
          max(abs(colMeans(gse45670_matrix_original))) > 1e-6 ||
          max(abs(apply(gse45670_matrix_original, 2L, sd) - 1)) > 1e-6,
        "锁定TCGA/GSE矩阵未保持各自整体逐基因均值0/SD1。")

factor_scores_path <- file.path(results_dir, "tcga_escc_mofa_factor_scores.tsv")
factor_scores <- fread(factor_scores_path, showProgress = FALSE)
target_factors <- c("Factor1", "Factor3")
fail_if(!all(c("patient_id", target_factors) %in% names(factor_scores)) ||
          uniqueN(factor_scores$patient_id) != 94L ||
          !all(rownames(tcga_matrix) %in% factor_scores$patient_id),
        "TCGA MOFA factor scores 字段或78例映射失败。")
training_targets <- as.matrix(
  factor_scores[match(rownames(tcga_matrix), patient_id), ..target_factors]
)
rownames(training_targets) <- rownames(tcga_matrix)
fail_if(any(!is.finite(training_targets)), "TCGA Factor1/Factor3 target 含NA/Inf。")

mofa_model_path <- file.path(results_dir, "tcga_escc_mofa_model.rds")
mofa_model <- readRDS(mofa_model_path)
mofa_weights <- as.data.table(get_weights(mofa_model, as.data.frame = TRUE))
rna_weights <- mofa_weights[view == "RNA" & factor %chin% target_factors]
rna_weights[, canonical_feature := sub("_RNA$", "", feature)]
rna_weights[, gene_id := sub("[.][0-9]+$", "", canonical_feature)]
fail_if(rna_weights[, uniqueN(gene_id), by = factor][, any(V1 != 1500L)] ||
          rna_weights[, anyDuplicated(gene_id), by = factor][, any(V1 != 0L)],
        "MOFA RNA 权重去 _RNA/版本后不唯一或维度异常。")

sensitivity_weights <- setNames(vector("list", length(target_factors)), target_factors)
sensitivity_diagnostics <- list()
sensitivity_definition_rows <- list()
for (factor_name in target_factors) {
  weight_table <- rna_weights[factor == factor_name & gene_id %chin% gene_features]
  setkey(weight_table, gene_id)
  overlap_genes <- intersect(gene_features, weight_table$gene_id)
  weight_table <- weight_table[overlap_genes]
  fail_if(nrow(weight_table) != 51L || anyNA(weight_table$value),
          paste("MOFA固定权重与ECMS特征交集不等于51：", factor_name))
  normalized_weight <- weight_table$value / sum(abs(weight_table$value))
  names(normalized_weight) <- weight_table$gene_id
  sensitivity_weights[[factor_name]] <- normalized_weight
  observed <- training_targets[, factor_name]
  predicted <- as.vector(
    tcga_matrix[, names(normalized_weight), drop = FALSE] %*% normalized_weight
  )
  total_rna_weight_mass <- rna_weights[factor == factor_name, sum(value^2)]
  overlap_weight_mass <- sum(weight_table$value^2) / total_rna_weight_mass
  sensitivity_diagnostics[[factor_name]] <- data.table(
    record_type = "sensitivity_in_sample_diagnostic",
    model_variant = "mofa51_fixed_weight_sensitivity",
    factor = factor_name,
    repeat_id = NA_integer_, n = length(observed),
    spearman_rho = cor(observed, predicted, method = "spearman"),
    pearson_r = cor(observed, predicted),
    rmse = sqrt(mean((observed - predicted)^2)),
    median_absolute_error = median(abs(observed - predicted)),
    spearman_min = NA_real_, spearman_max = NA_real_,
    positive_orientation = cor(observed, predicted, method = "spearman") > 0,
    outer_folds = NA_integer_, inner_folds = NA_integer_,
    selection_rule = "fixed_original_MOFA_RNA_weights_no_outcome_fit",
    summary_statistic = "in_sample_transparency_diagnostic",
    gate_status = "sensitivity_only_not_a_gate",
    feature_count = 51L,
    rna_squared_weight_mass_fraction = overlap_weight_mass
  )
  sensitivity_definition_rows[[factor_name]] <- data.table(
    factor = factor_name,
    model_variant = "mofa51_fixed_weight_sensitivity",
    term = weight_table$gene_id,
    gene_id = weight_table$gene_id,
    coefficient = as.numeric(normalized_weight),
    source_weight = weight_table$value,
    selected_nonzero = normalized_weight != 0,
    feature_count = 51L,
    training_n = 78L,
    lambda_1se = NA_real_,
    alpha = NA_real_,
    preprocessing = paste(
      "ECMS matrix already cohort-wise gene Z; external GSE53622/24 separately",
      "re-Z; remove MOFA _RNA suffix then Ensembl version"
    ),
    model_role = "transparent_sensitivity_not_primary",
    rna_squared_weight_mass_fraction = overlap_weight_mass
  )
}
rm(mofa_model, mofa_weights)
gc(verbose = FALSE)

balanced_folds <- function(n, k, seed) {
  set.seed(seed)
  sample(rep(seq_len(k), length.out = n), n, replace = FALSE)
}

fit_repeated_nested_ridge <- function(x, y, factor_name, factor_index) {
  repeat_rows <- list()
  oof_rows <- list()
  for (repeat_id in seq_len(20L)) {
    outer_fold <- balanced_folds(
      nrow(x), 5L, 202607160L + factor_index * 10000L + repeat_id
    )
    prediction <- rep(NA_real_, nrow(x))
    lambda_by_patient <- rep(NA_real_, nrow(x))
    for (fold in seq_len(5L)) {
      test_index <- which(outer_fold == fold)
      train_index <- which(outer_fold != fold)
      inner_fold <- balanced_folds(
        length(train_index), 5L,
        202607160L + factor_index * 100000L + repeat_id * 100L + fold
      )
      inner_cv <- cv.glmnet(
        x = x[train_index, , drop = FALSE], y = y[train_index],
        family = "gaussian", alpha = 0, foldid = inner_fold,
        nfolds = 5L, type.measure = "mse", standardize = TRUE,
        intercept = TRUE, nlambda = 100L,
        control = list(maxit = 100000L),
        parallel = FALSE, keep = FALSE
      )
      lambda <- inner_cv$lambda.1se
      fail_if(length(lambda) != 1L || !is.finite(lambda) || lambda <= 0,
              paste("inner CV lambda.1se 无效：", factor_name, repeat_id, fold))
      prediction[test_index] <- as.numeric(stats::predict(
        inner_cv, newx = x[test_index, , drop = FALSE], s = "lambda.1se"
      ))
      lambda_by_patient[test_index] <- lambda
    }
    fail_if(any(!is.finite(prediction)) || any(!is.finite(lambda_by_patient)),
            paste("outer OOF prediction 不完整：", factor_name, repeat_id))
    rho <- cor(y, prediction, method = "spearman")
    repeat_rows[[repeat_id]] <- data.table(
      record_type = "outer_repeat",
      model_variant = "ridge314_primary",
      factor = factor_name,
      repeat_id = repeat_id,
      n = length(y),
      spearman_rho = rho,
      pearson_r = cor(y, prediction),
      rmse = sqrt(mean((y - prediction)^2)),
      median_absolute_error = median(abs(y - prediction)),
      spearman_min = NA_real_, spearman_max = NA_real_,
      positive_orientation = rho > 0,
      outer_folds = 5L, inner_folds = 5L,
      selection_rule = "nested_inner_5fold_lambda.1se",
      summary_statistic = "repeat",
      gate_status = if (rho >= 0.70) "GO" else if (rho >= 0.50) {
        "CONDITIONAL_GO"
      } else {
        "NO_GO"
      },
      feature_count = ncol(x),
      rna_squared_weight_mass_fraction = NA_real_
    )
    oof_rows[[repeat_id]] <- data.table(
      model_variant = "ridge314_primary",
      factor = factor_name,
      repeat_id = repeat_id,
      patient_id = rownames(x),
      outer_fold = outer_fold,
      observed_factor_score = y,
      oof_predicted_factor_score = prediction,
      residual = y - prediction,
      inner_selected_lambda_1se = lambda_by_patient,
      orientation_flipped_after_outcome = FALSE
    )
  }
  repeats <- rbindlist(repeat_rows)
  median_rho <- median(repeats$spearman_rho)
  summary_row <- data.table(
    record_type = "outer_repeat_summary",
    model_variant = "ridge314_primary",
    factor = factor_name,
    repeat_id = NA_integer_, n = length(y),
    spearman_rho = median_rho,
    pearson_r = median(repeats$pearson_r),
    rmse = median(repeats$rmse),
    median_absolute_error = median(repeats$median_absolute_error),
    spearman_min = min(repeats$spearman_rho),
    spearman_max = max(repeats$spearman_rho),
    positive_orientation = all(repeats$positive_orientation),
    outer_folds = 5L, inner_folds = 5L,
    selection_rule = "median_of_20_repeated_outer_5fold_CV",
    summary_statistic = "median_with_min_max",
    gate_status = if (median_rho >= 0.70) "GO" else if (median_rho >= 0.50) {
      "CONDITIONAL_GO"
    } else {
      "NO_GO"
    },
    feature_count = ncol(x),
    rna_squared_weight_mass_fraction = NA_real_
  )
  list(
    metrics = rbindlist(list(repeats, summary_row), use.names = TRUE),
    oof = rbindlist(oof_rows)
  )
}

message("[3/10] 执行 20×5 外层重复CV与内层5折 lambda.1se")
nested_results <- setNames(vector("list", length(target_factors)), target_factors)
for (factor_index in seq_along(target_factors)) {
  factor_name <- target_factors[[factor_index]]
  nested_results[[factor_name]] <- fit_repeated_nested_ridge(
    tcga_matrix, training_targets[, factor_name], factor_name, factor_index
  )
}
internal_cv <- rbindlist(lapply(nested_results, `[[`, "metrics"), fill = TRUE)
internal_cv <- rbindlist(
  c(list(internal_cv), sensitivity_diagnostics), fill = TRUE, use.names = TRUE
)
oof_predictions <- rbindlist(lapply(nested_results, `[[`, "oof"))
setorder(internal_cv, factor, model_variant, record_type, repeat_id)
setorder(oof_predictions, factor, repeat_id, patient_id)

final_models <- setNames(vector("list", length(target_factors)), target_factors)
ridge_definition_rows <- list()
for (factor_index in seq_along(target_factors)) {
  factor_name <- target_factors[[factor_index]]
  final_fold <- balanced_folds(
    nrow(tcga_matrix), 5L, 202607160L + factor_index * 1000000L
  )
  fit <- cv.glmnet(
    x = tcga_matrix, y = training_targets[, factor_name],
    family = "gaussian", alpha = 0, foldid = final_fold, nfolds = 5L,
    type.measure = "mse", standardize = TRUE, intercept = TRUE,
    nlambda = 100L, control = list(maxit = 100000L),
    parallel = FALSE, keep = FALSE
  )
  coefficient_matrix <- as.matrix(coef(fit, s = "lambda.1se"))
  fail_if(!identical(rownames(coefficient_matrix), c("(Intercept)", gene_features)) ||
            any(!is.finite(coefficient_matrix[, 1L])) ||
            !is.finite(fit$lambda.1se) || fit$lambda.1se <= 0,
          paste("最终 ridge 模型系数或 lambda.1se 失败：", factor_name))
  final_models[[factor_name]] <- list(
    cv_glmnet_fit = fit,
    lambda_1se = fit$lambda.1se,
    final_foldid = final_fold,
    coefficient = setNames(coefficient_matrix[, 1L], rownames(coefficient_matrix))
  )
  ridge_definition_rows[[factor_name]] <- data.table(
    factor = factor_name,
    model_variant = "ridge314_primary",
    term = rownames(coefficient_matrix),
    gene_id = fifelse(rownames(coefficient_matrix) == "(Intercept)",
                      NA_character_, rownames(coefficient_matrix)),
    coefficient = as.numeric(coefficient_matrix[, 1L]),
    source_weight = NA_real_,
    selected_nonzero = as.numeric(coefficient_matrix[, 1L]) != 0,
    feature_count = 314L,
    training_n = 78L,
    lambda_1se = fit$lambda.1se,
    alpha = 0,
    preprocessing = paste(
      "locked repository tcga.val.df; repeated nested CV; glmnet standardize=TRUE;",
      "external GSE53622/24 separately gene-wise re-Z"
    ),
    model_role = "primary_RNA_proxy_not_full_MOFA_reproduction",
    rna_squared_weight_mass_fraction = NA_real_
  )
}

proxy_definition <- rbindlist(
  c(ridge_definition_rows, sensitivity_definition_rows),
  use.names = TRUE, fill = TRUE
)
proxy_definition[, orientation_rule :=
  "retain original MOFA factor orientation; never flip after external outcome"]
proxy_definition[, external_outcomes_seen_during_definition := FALSE]
proxy_definition[, term_is_intercept := is.na(gene_id)]
setorder(proxy_definition, factor, model_variant, term_is_intercept, term)
proxy_definition[, term_is_intercept := NULL]

model_definition_payload <- list(
  feature_ids = gene_features,
  target_factors = target_factors,
  training_patient_ids = rownames(tcga_matrix),
  final_lambda = vapply(final_models, `[[`, numeric(1), "lambda_1se"),
  final_coefficients = lapply(final_models, `[[`, "coefficient"),
  sensitivity_weights = sensitivity_weights,
  preprocessing = paste(
    "ECMS locked cohort-wise gene Z; GSE53622/GSE53624 separately re-Z;",
    "GSE45670 locked tumor-only gene Z"
  ),
  outer_cv = "20 repeats x 5 folds; inner 5-fold lambda.1se",
  orientation = "fixed_original_MOFA_factor_orientation"
)
model_definition_sha256 <- digest(
  model_definition_payload, algo = "sha256", serialize = TRUE
)
proxy_definition[, model_definition_sha256 := model_definition_sha256]

proxy_models <- list(
  schema_version = "1.0",
  generated_date = as.character(Sys.Date()),
  model_definition_sha256 = model_definition_sha256,
  frozen_before_external_outcomes = TRUE,
  model_commit = locked_commit,
  model_sha256 = locked_model_sha256,
  feature_list_sha256 = locked_feature_sha256,
  feature_ids = gene_features,
  target_factors = target_factors,
  training_patient_ids = rownames(tcga_matrix),
  outer_cv_specification = "20 repeated outer 5-fold; inner 5-fold; lambda.1se",
  orientation_rule = "fixed original MOFA orientation; no outcome-directed flip",
  external_preprocessing = list(
    GSE53622 = "within-GSE53622 gene-wise Z",
    GSE53624 = "within-GSE53624 gene-wise Z",
    GSE45670 = "locked author tumor-only gene-wise Z"
  ),
  ridge314_primary = final_models,
  mofa51_fixed_weight_sensitivity = sensitivity_weights
)

message("[4/10] 在读取外部结局前冻结 proxy RDS、定义表和内部CV")
proxy_model_stage_path <- file.path(stage_dir, "escc_external_state_proxy_models.rds")
saveRDS(proxy_models, proxy_model_stage_path, compress = "gzip")
fail_if(!file_exists(proxy_model_stage_path) ||
          as.numeric(file_info(proxy_model_stage_path)$size) <= 0,
        "冻结 proxy RDS 写入失败。")
proxy_model_sha256 <- sha256_file(proxy_model_stage_path)
proxy_check <- readRDS(proxy_model_stage_path)
fail_if(!isTRUE(proxy_check$frozen_before_external_outcomes) ||
          proxy_check$model_definition_sha256 != model_definition_sha256 ||
          !identical(proxy_check$feature_ids, gene_features),
        "冻结 proxy RDS 回读身份失败。")
rm(proxy_models, final_models, nested_results)
gc(verbose = FALSE)
frozen_proxy <- readRDS(proxy_model_stage_path)

stage_fwrite(proxy_definition, "escc_external_state_proxy_definition.tsv")
stage_fwrite(internal_cv, "escc_external_state_internal_cv.tsv")
stage_fwrite(oof_predictions, "escc_external_state_oof_predictions.tsv")

for (factor_name in target_factors) {
  summary_row <- internal_cv[
    factor == factor_name & model_variant == "ridge314_primary" &
      record_type == "outer_repeat_summary"
  ]
  fail_if(nrow(summary_row) != 1L || summary_row$gate_status == "NO_GO" ||
            !summary_row$positive_orientation,
          paste("ridge proxy 内部CV未达到最低条件门禁：", factor_name))
  add_decision(
    paste0("PROXY_CV_", toupper(factor_name)), "proxy_fidelity",
    factor = factor_name, hard_gate = TRUE,
    metric = "median repeated outer-CV Spearman rho",
    observed = sprintf(
      "%.3f [min %.3f, max %.3f]",
      summary_row$spearman_rho, summary_row$spearman_min,
      summary_row$spearman_max
    ),
    rule = "GO>=0.70; CONDITIONAL_GO=0.50-0.69; NO_GO<0.50",
    status = summary_row$gate_status,
    counts_toward_external_support = TRUE,
    interpretation = "314基因RNA proxy在TCGA外层留出预测中保持原方向。",
    boundary = "这是同TCGA内proxy fidelity，不是外部生物学验证。"
  )
}
add_decision(
  "PROXY_MODEL_FROZEN", "leakage_control", hard_gate = TRUE,
  metric = "model frozen before external outcome parse",
  observed = paste0("TRUE;proxy_rds_sha256=", proxy_model_sha256),
  rule = "features, coefficients, lambda, direction and preprocessing frozen first",
  status = "PASS", counts_toward_external_support = FALSE,
  interpretation = "外部结局不能反向影响模型定义。",
  boundary = "外部分析仍属于公共数据关联研究。"
)

parse_tab_line <- function(line) {
  fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
  sub('"$', "", sub('^"', "", fields))
}

read_geo_header <- function(path) {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  lines <- character()
  for (index in seq_len(10000L)) {
    line <- readLines(connection, n = 1L, warn = FALSE)
    fail_if(!length(line), paste("GEO series matrix 缺少 table begin：", path))
    if (startsWith(line, "!series_matrix_table_begin")) break
    lines <- c(lines, line)
  }
  lines
}

normalize_field_name <- function(value) {
  value <- tolower(trimws(value))
  value <- gsub("[^a-z0-9]+", "_", value)
  value <- gsub("^_+|_+$", "", value)
  value[value == "tumor_loation"] <- "tumor_location"
  value
}

read_geo_sample_metadata <- function(path) {
  lines <- read_geo_header(path)
  get_single <- function(key) {
    candidates <- lines[startsWith(lines, paste0(key, "\t"))]
    fail_if(length(candidates) != 1L, paste("GEO header key不唯一：", key, path))
    parse_tab_line(candidates[[1L]])[-1L]
  }
  gsm <- get_single("!Sample_geo_accession")
  title <- get_single("!Sample_title")
  source_name <- get_single("!Sample_source_name_ch1")
  fail_if(length(gsm) != length(title) || length(gsm) != length(source_name) ||
            anyDuplicated(gsm),
          paste("GEO sample header维度或GSM唯一性失败：", path))
  output <- data.table(gsm = gsm, title = title, source_name = source_name)
  characteristic_lines <- lines[startsWith(lines, "!Sample_characteristics_ch1\t")]
  for (line in characteristic_lines) {
    values <- parse_tab_line(line)[-1L]
    fail_if(length(values) != length(gsm),
            paste("GEO characteristic长度不符：", path))
    key <- normalize_field_name(sub(":.*$", "", values[[1L]]))
    parsed <- trimws(sub("^[^:]+:[[:space:]]*", "", values))
    fail_if(!nzchar(key) || key %in% names(output),
            paste("GEO characteristic key空白或重复：", key, path))
    output[, (key) := parsed]
  }
  output
}

read_gzipped_excel <- function(path) {
  destination <- file.path(stage_dir, sub("[.]gz$", "", basename(path)))
  R.utils::gunzip(path, destname = destination, overwrite = TRUE, remove = FALSE)
  on.exit(if (file_exists(destination)) file_delete(destination), add = TRUE)
  sheets <- readxl::excel_sheets(destination)
  candidates <- lapply(sheets, function(sheet) {
    as.data.table(readxl::read_excel(destination, sheet = sheet))
  })
  row_counts <- vapply(candidates, nrow, integer(1))
  fail_if(!length(row_counts) || max(row_counts) <= 0,
          paste("Excel没有非空sheet：", path))
  candidates[[which.max(row_counts)]]
}

normalize_gse536_clinical <- function(object, cohort) {
  x <- copy(object)
  setnames(x, names(x), normalize_field_name(names(x)))
  required <- c(
    "patient_id", "age", "sex", "tobacco_use", "alcohol_use",
    "tumor_location", "tumor_grade", "t_stage", "n_stage", "tnm_stage",
    "arrhythmia", "pneumonia", "anastomotic_leak", "adjuvant_therapy",
    "death_at_fu", "survival_time_months"
  )
  fail_if(!setequal(names(x), required) || nrow(x) == 0L,
          paste("GSE536临床字段失败：", cohort))
  character_columns <- setdiff(required, c("age", "survival_time_months"))
  x[, (character_columns) := lapply(.SD, function(value) {
    tolower(trimws(as.character(value)))
  }), .SDcols = character_columns]
  x[, patient_id := tolower(trimws(patient_id))]
  x[, age := as.numeric(age)]
  x[, survival_time_months := as.numeric(survival_time_months)]
  fail_if(anyDuplicated(x$patient_id) || any(!is.finite(x$age)) ||
            any(!is.finite(x$survival_time_months)) ||
            any(x$survival_time_months <= 0) ||
            any(!x$death_at_fu %chin% c("yes", "no")) ||
            any(!x$sex %chin% c("female", "male")) ||
            any(!x$tnm_stage %chin% c("i", "ii", "iii")),
          paste("GSE536临床数值、编码或患者唯一性失败：", cohort))
  x[, os_event := as.integer(death_at_fu == "yes")]
  x[, tnm_stage_numeric := match(tnm_stage, c("i", "ii", "iii"))]
  x[, cohort := cohort]
  x
}

compare_geo_excel_clinical <- function(tumor_metadata, clinical, cohort) {
  geo <- copy(tumor_metadata)
  if ("tumor_loation" %in% names(geo) && !"tumor_location" %in% names(geo)) {
    setnames(geo, "tumor_loation", "tumor_location")
  }
  compare_fields <- c(
    "age", "sex", "tobacco_use", "alcohol_use", "tumor_location",
    "tumor_grade", "t_stage", "n_stage", "tnm_stage", "arrhythmia",
    "pneumonia", "anastomotic_leak", "adjuvant_therapy", "death_at_fu",
    "survival_time_months"
  )
  fail_if(!all(c("patient_id", compare_fields) %in% names(geo)),
          paste("GEO临床字段不完整：", cohort))
  geo[, patient_id := tolower(trimws(patient_id))]
  geo <- unique(geo[, c("patient_id", compare_fields), with = FALSE])
  fail_if(anyDuplicated(geo$patient_id) || !setequal(geo$patient_id, clinical$patient_id),
          paste("GEO tumor患者与Excel患者集合不一致：", cohort))
  geo <- geo[match(clinical$patient_id, patient_id)]
  mismatch <- 0L
  for (field in compare_fields) {
    if (field %chin% c("age", "survival_time_months")) {
      mismatch <- mismatch + sum(
        abs(as.numeric(geo[[field]]) - as.numeric(clinical[[field]])) > 1e-8
      )
    } else {
      mismatch <- mismatch + sum(
        tolower(trimws(as.character(geo[[field]]))) !=
          tolower(trimws(as.character(clinical[[field]])))
      )
    }
  }
  fail_if(mismatch != 0L, paste("GEO与Excel临床逐字段不一致：", cohort, mismatch))
  mismatch
}

scale_external_matrix <- function(x, cohort) {
  standard_deviation <- apply(x, 2L, sd)
  fail_if(any(!is.finite(standard_deviation)) || any(standard_deviation <= 0),
          paste("外部矩阵存在零/异常方差：", cohort))
  z <- scale(x, center = TRUE, scale = TRUE)
  fail_if(any(!is.finite(z)) || max(abs(colMeans(z))) > 1e-6 ||
            max(abs(apply(z, 2L, sd) - 1)) > 1e-6,
          paste("外部矩阵队列内重标化失败：", cohort))
  z
}

score_proxy_matrix <- function(x, frozen_model) {
  fail_if(!identical(colnames(x), frozen_model$feature_ids),
          "proxy scoring特征顺序与冻结模型不一致。")
  output <- data.table(gsm = rownames(x))
  for (factor_name in frozen_model$target_factors) {
    ridge_fit <- frozen_model$ridge314_primary[[factor_name]]$cv_glmnet_fit
    ridge_score <- as.numeric(stats::predict(
      ridge_fit, newx = x, s = "lambda.1se"
    ))
    weight <- frozen_model$mofa51_fixed_weight_sensitivity[[factor_name]]
    sensitivity_score <- as.numeric(
      x[, names(weight), drop = FALSE] %*% weight
    )
    prefix <- tolower(factor_name)
    output[, (paste0(prefix, "_ridge314_primary")) := ridge_score]
    output[, (paste0(prefix, "_mofa51_sensitivity")) := sensitivity_score]
  }
  output
}

add_within_cohort_score_z <- function(object) {
  x <- copy(object)
  raw_score_columns <- grep(
    "^factor[13]_(ridge314_primary|mofa51_sensitivity)$",
    names(x), value = TRUE
  )
  fail_if(length(raw_score_columns) != 4L,
          "患者分数表缺少4个固定raw score列。")
  for (column in raw_score_columns) {
    standard_deviation <- sd(x[[column]])
    fail_if(!is.finite(standard_deviation) || standard_deviation <= 0,
            paste("外部proxy score零方差：", x$subseries[[1L]], column))
    x[, (paste0(column, "_z")) := as.numeric(scale(get(column)))]
  }
  x
}

message("[5/10] 模型冻结后解析GEO样本、临床、生存和疗效标签")
gse53622_series <- hub_file(
  hub_datasets$gse53622, "00_source/GSE53622_series_matrix.txt.gz"
)
gse53624_series <- hub_file(
  hub_datasets$gse53624, "00_source/GSE53624_series_matrix.txt.gz"
)
gse53625_series <- hub_file(
  hub_datasets$gse53625, "00_source/GSE53625_series_matrix.txt.gz"
)
gse45670_series <- hub_file(
  hub_datasets$gse45670, "00_source/GSE45670_series_matrix.txt.gz"
)
gse53622_excel <- hub_file(
  hub_datasets$gse53622,
  "10_metadata/GSE53622_clinical_data_of_patients_independent_set.xlsx.gz"
)
gse53624_excel <- hub_file(
  hub_datasets$gse53624,
  "10_metadata/GSE53624_clinical_data_of_patients_orignial_set.xlsx.gz"
)

metadata_53622 <- read_geo_sample_metadata(gse53622_series)
metadata_53624 <- read_geo_sample_metadata(gse53624_series)
metadata_53625 <- read_geo_sample_metadata(gse53625_series)
metadata_45670 <- read_geo_sample_metadata(gse45670_series)
fail_if(nrow(metadata_53622) != 120L || nrow(metadata_53624) != 238L ||
          nrow(metadata_53625) != 358L || nrow(metadata_45670) != 38L,
        "GEO样本数未复现120/238/358/38。")
fail_if(!setequal(
  metadata_53625$gsm, c(metadata_53622$gsm, metadata_53624$gsm)
) || length(intersect(metadata_53622$gsm, metadata_53624$gsm)) != 0L,
"GSE53625不等于GSE53622+GSE53624的精确GSM并集。")

tumor_53622 <- metadata_53622[
  grepl("^cancer tissue", tolower(tissue)), .(
    gsm, patient_id = tolower(trimws(patient_id)), tissue
  )
]
tumor_53624 <- metadata_53624[
  grepl("^cancer tissue", tolower(tissue)), .(
    gsm, patient_id = tolower(trimws(patient_id)), tissue
  )
]
fail_if(nrow(tumor_53622) != 60L || nrow(tumor_53624) != 119L ||
          uniqueN(tumor_53622$patient_id) != 60L ||
          uniqueN(tumor_53624$patient_id) != 119L ||
          length(intersect(tumor_53622$patient_id, tumor_53624$patient_id)) != 0L,
        "GSE53622/24肿瘤患者数、唯一性或患者独立性失败。")
fail_if(!setequal(
  rownames(gse53625_matrix_original), c(tumor_53622$gsm, tumor_53624$gsm)
), "官方gse53625.val.df行名不等于两个子系列肿瘤GSM并集。")

clinical_53622 <- normalize_gse536_clinical(
  read_gzipped_excel(gse53622_excel), "GSE53622"
)
clinical_53624 <- normalize_gse536_clinical(
  read_gzipped_excel(gse53624_excel), "GSE53624"
)
compare_geo_excel_clinical(
  metadata_53622[grepl("^cancer tissue", tolower(tissue))],
  clinical_53622, "GSE53622"
)
compare_geo_excel_clinical(
  metadata_53624[grepl("^cancer tissue", tolower(tissue))],
  clinical_53624, "GSE53624"
)
fail_if(nrow(clinical_53622) != 60L || sum(clinical_53622$os_event) != 33L ||
          nrow(clinical_53624) != 119L || sum(clinical_53624$os_event) != 73L,
        "GSE53622/24临床样本或事件数未复现60/33与119/73。")

add_decision(
  "MAP_GSE53622", "patient_mapping", hard_gate = TRUE,
  metric = "tumor GSM -> patient_id -> clinical row",
  observed = "60/60;33 OS events;0 clinical mismatch",
  rule = "100% unique mapping; GEO and Excel 15 clinical fields agree",
  status = "PASS", counts_toward_external_support = TRUE,
  interpretation = "GSE53622可作为同研究家族的独立子系列估计。",
  boundary = "与GSE53624同属一个研究家族。"
)
add_decision(
  "MAP_GSE53624", "patient_mapping", hard_gate = TRUE,
  metric = "tumor GSM -> patient_id -> clinical row",
  observed = "119/119;73 OS events;0 clinical mismatch",
  rule = "100% unique mapping; GEO and Excel 15 clinical fields agree",
  status = "PASS", counts_toward_external_support = TRUE,
  interpretation = "GSE53624可作为同研究家族的另一个子系列估计。",
  boundary = "不是与GSE53622完全独立的研究来源。"
)
add_decision(
  "INDEPENDENCE_GSE536_FAMILY", "independence", hard_gate = TRUE,
  metric = "GSE53625 relationship",
  observed = "GSE53625=GSE53622(60 tumors)+GSE53624(119 tumors);0 patient overlap",
  rule = "GSE53625 never counted as a third cohort",
  status = "PASS", counts_toward_external_support = FALSE,
  interpretation = "分别报告子系列并用strata(subseries)合并。",
  boundary = "最多称一个外部研究家族、两个不重叠子系列。"
)

map_clinical_scores <- function(tumor_map, clinical, matrix_original, subseries) {
  mapping <- merge(tumor_map[, .(gsm, patient_id)], clinical,
                   by = "patient_id", all = FALSE, sort = FALSE)
  fail_if(nrow(mapping) != nrow(tumor_map) || anyDuplicated(mapping$gsm) ||
            !setequal(mapping$gsm, tumor_map$gsm),
          paste("肿瘤GSM与临床合并失败：", subseries))
  mapping <- mapping[match(tumor_map$gsm, gsm)]
  x <- matrix_original[mapping$gsm, , drop = FALSE]
  x <- scale_external_matrix(x, subseries)
  score <- score_proxy_matrix(x, frozen_proxy)
  mapping <- merge(mapping, score, by = "gsm", all = FALSE, sort = FALSE)
  mapping <- mapping[match(tumor_map$gsm, gsm)]
  mapping[, `:=`(
    dataset = subseries,
    independence_family = "GSE53625_family",
    subseries = subseries,
    sample_role = "ESCC_tumor",
    response_binary = NA_integer_,
    response_label = NA_character_,
    external_endpoint_role = "secondary_overall_survival"
  )]
  add_within_cohort_score_z(mapping)
}

patient_53622 <- map_clinical_scores(
  tumor_53622, clinical_53622, gse53625_matrix_original, "GSE53622"
)
patient_53624 <- map_clinical_scores(
  tumor_53624, clinical_53624, gse53625_matrix_original, "GSE53624"
)

# ECMS标签保留锁定模型在作者179矩阵上的投影；不对子系列重标化后重定义标签。
# randomForest含偶数棵树，response模式会随机打破并列。这里使用固定类水平顺序的
# 概率argmax，完全不读取结局，且把并列患者与规则显式写入患者表和decision。
gse53625_ecms <- deterministic_rf_projection(
  rf_classifier, gse53625_matrix_original, ecms_levels
)
attach_ecms_projection <- function(patient, projection, cohort) {
  patient_order <- patient$gsm
  output <- merge(patient, projection, by = "gsm", all = FALSE, sort = FALSE)
  output <- output[match(patient_order, gsm)]
  fail_if(nrow(output) != nrow(patient) || anyNA(output$ecms_label) ||
            anyDuplicated(output$gsm),
          paste("ECMS投影与患者表映射失败：", cohort))
  output
}
patient_53622 <- attach_ecms_projection(
  patient_53622, gse53625_ecms, "GSE53622"
)
patient_53624 <- attach_ecms_projection(
  patient_53624, gse53625_ecms, "GSE53624"
)
ecms_label_counts <- as.integer(table(factor(
  c(patient_53622$ecms_label, patient_53624$ecms_label),
  levels = ecms_levels
)))
gse53625_tie_patients <- gse53625_ecms[ecms_tie_at_max == TRUE, gsm]
fail_if(sum(ecms_label_counts) != 179L ||
          !identical(sort(c(patient_53622$gsm, patient_53624$gsm)),
                     sort(gse53625_ecms$gsm)),
        "GSE53625确定性ECMS投影未完整覆盖179例。")
add_decision(
  "ECMS_GSE53625_DETERMINISTIC_TIE", "model_projection",
  hard_gate = TRUE,
  metric = "fixed probability argmax and explicit RF vote ties",
  observed = paste0(
    "counts=", paste(ecms_label_counts, collapse = "/"),
    ";tie_count=", length(gse53625_tie_patients),
    ";tie_patients=", ifelse(
      length(gse53625_tie_patients),
      paste(gse53625_tie_patients, collapse = ","), "none"
    )
  ),
  rule = unique(gse53625_ecms$ecms_label_rule),
  status = "PASS", counts_toward_external_support = FALSE,
  interpretation = "锁定作者矩阵上的ECMS标签不再受RNG状态影响。",
  boundary = paste(
    "旧48/63/35/33为response随机打破并列的快照；",
    "确定性计数不以生存或疗效结局选择。"
  )
)

tumor_45670 <- metadata_45670[
  grepl("esophageal squamous cell carcinoma", tolower(tissue)),
  .(
    gsm,
    patient_id = gsm,
    age = as.numeric(age),
    sex = tolower(trimws(gender)),
    tumor_stage = trimws(tumor_stage),
    response_source_field = trimws(pathological_response_to_crp)
  )
]
tumor_45670[, response_binary := fcase(
  response_source_field == "pathological complete response", 1L,
  response_source_field == "not pathological complete response", 0L,
  default = NA_integer_
)]
fail_if(nrow(tumor_45670) != 28L || uniqueN(tumor_45670$gsm) != 28L ||
          !setequal(tumor_45670$gsm, rownames(gse45670_matrix_original)) ||
          sum(tumor_45670$response_binary == 1L) != 11L ||
          sum(tumor_45670$response_binary == 0L) != 17L ||
          any(!is.finite(tumor_45670$age)) || anyNA(tumor_45670$sex) ||
          anyNA(tumor_45670$tumor_stage),
        "GSE45670 28肿瘤/17非pCR/11pCR或临床映射失败。")
tumor_45670 <- tumor_45670[match(rownames(gse45670_matrix_original), gsm)]
score_45670 <- score_proxy_matrix(gse45670_matrix_original, frozen_proxy)
patient_45670 <- merge(
  tumor_45670, score_45670, by = "gsm", all = FALSE, sort = FALSE
)
patient_45670 <- patient_45670[match(rownames(gse45670_matrix_original), gsm)]
gse45670_ecms <- deterministic_rf_projection(
  rf_classifier, gse45670_matrix_original, ecms_levels
)
patient_45670 <- attach_ecms_projection(
  patient_45670, gse45670_ecms, "GSE45670"
)
patient_45670[, `:=`(
  dataset = "GSE45670",
  independence_family = "GSE45670",
  subseries = "GSE45670",
  sample_role = "ESCC_tumor",
  tumor_location = NA_character_,
  tumor_grade = NA_character_,
  t_stage = sub("N.*$", "", tumor_stage),
  n_stage = sub("^T[0-9]+", "", sub("M.*$", "", tumor_stage)),
  tnm_stage = NA_character_,
  tnm_stage_numeric = NA_integer_,
  tobacco_use = NA_character_,
  alcohol_use = NA_character_,
  arrhythmia = NA_character_,
  pneumonia = NA_character_,
  anastomotic_leak = NA_character_,
  adjuvant_therapy = NA_character_,
  death_at_fu = NA_character_,
  survival_time_months = NA_real_,
  os_event = NA_integer_,
  response_label = fifelse(
    response_binary == 1L, "pathological_complete_response",
    "not_pathological_complete_response"
  ),
  external_endpoint_role = "exploratory_pathological_response"
)]
patient_45670 <- add_within_cohort_score_z(patient_45670)
gse45670_tie_patients <- gse45670_ecms[ecms_tie_at_max == TRUE, gsm]
add_decision(
  "MAP_GSE45670", "patient_mapping", hard_gate = TRUE,
  metric = "official 28-row matrix -> tumor GSM -> response metadata",
  observed = paste0(
    "28/28;17 non-pCR;11 pCR;0 missing age/sex/stage/response;tie_count=",
    length(gse45670_tie_patients), ";tie_patients=",
    ifelse(length(gse45670_tie_patients),
           paste(gse45670_tie_patients, collapse = ","), "none")
  ),
  rule = paste(
    "100% tumor GSM mapping and exact response counts;",
    unique(gse45670_ecms$ecms_label_rule)
  ),
  status = "PASS", counts_toward_external_support = TRUE,
  interpretation = "可进行预设的探索性疗效关联。",
  boundary = "n=28；不能称治疗反应预测器。"
)

patient_scores <- rbindlist(
  list(patient_53622, patient_53624, patient_45670),
  use.names = TRUE, fill = TRUE
)
preferred_columns <- c(
  "dataset", "independence_family", "subseries", "gsm", "patient_id",
  "sample_role", "ecms_label", "ecms_max_vote_fraction",
  "ecms_tie_count_at_max", "ecms_tie_at_max", "ecms_tied_classes",
  "ecms_label_rule", paste0("ecms_vote_", ecms_levels),
  "factor1_ridge314_primary", "factor1_ridge314_primary_z",
  "factor3_ridge314_primary", "factor3_ridge314_primary_z",
  "factor1_mofa51_sensitivity", "factor1_mofa51_sensitivity_z",
  "factor3_mofa51_sensitivity", "factor3_mofa51_sensitivity_z",
  "age", "sex", "tnm_stage", "tnm_stage_numeric", "t_stage", "n_stage",
  "tumor_grade", "tumor_location", "tobacco_use", "alcohol_use",
  "adjuvant_therapy", "survival_time_months", "os_event",
  "response_binary", "response_label", "response_source_field",
  "external_endpoint_role"
)
fail_if(!all(preferred_columns %in% names(patient_scores)),
        "patient_scores缺少预锁定字段。")
setcolorder(
  patient_scores,
  c(preferred_columns, setdiff(names(patient_scores), preferred_columns))
)
setorder(patient_scores, dataset, patient_id)

score_specs <- data.table(
  factor = rep(target_factors, each = 2L),
  model_variant = rep(c(
    "ridge314_primary", "mofa51_fixed_weight_sensitivity"
  ), times = 2L),
  score_column = c(
    "factor1_ridge314_primary_z", "factor1_mofa51_sensitivity_z",
    "factor3_ridge314_primary_z", "factor3_mofa51_sensitivity_z"
  )
)

fit_cox_association <- function(
    data, cohort_scope, factor_name, model_variant, score_column,
    adjusted, stratified) {
  columns <- c("survival_time_months", "os_event", score_column)
  if (adjusted) columns <- c(columns, "age", "sex", "tnm_stage_numeric")
  if (stratified) columns <- c(columns, "subseries")
  d <- copy(data[complete.cases(data[, ..columns])])
  model_data <- data.frame(
    time = d$survival_time_months,
    event = d$os_event,
    score = d[[score_column]],
    age = d$age,
    sex = factor(d$sex, levels = c("female", "male")),
    stage = d$tnm_stage_numeric,
    subseries = factor(d$subseries)
  )
  formula_text <- if (adjusted && stratified) {
    "Surv(time,event) ~ score + age + sex + stage + strata(subseries)"
  } else if (!adjusted && stratified) {
    "Surv(time,event) ~ score + strata(subseries)"
  } else if (adjusted) {
    "Surv(time,event) ~ score + age + sex + stage"
  } else {
    "Surv(time,event) ~ score"
  }
  fit <- coxph(
    as.formula(formula_text), data = model_data,
    ties = "efron", x = TRUE, y = TRUE, model = TRUE
  )
  fit_summary <- summary(fit)
  coefficient <- fit_summary$coefficients["score", ]
  standard_error <- unname(coefficient[["se(coef)"]])
  zph <- cox.zph(fit, transform = "km")$table
  score_row <- grep("^score$", rownames(zph))
  global_row <- grep("^GLOBAL$", rownames(zph))
  fail_if(length(score_row) != 1L || length(global_row) != 1L,
          paste("cox.zph缺少score/GLOBAL：", cohort_scope, factor_name, model_variant))
  data.table(
    cohort_scope = cohort_scope,
    independence_family = "GSE53625_family",
    model_variant = model_variant,
    factor = factor_name,
    model_type = if (adjusted) "clinical_adjusted" else "unadjusted",
    stratified_by_subseries = stratified,
    model_formula = formula_text,
    n_patients = nrow(model_data),
    n_events = sum(model_data$event),
    patient_id_set_sha256 = feature_list_sha256(sort(as.character(d$patient_id))),
    subseries_levels = paste(sort(unique(as.character(d$subseries))), collapse = ","),
    hazard_ratio_per_1sd = exp(unname(coefficient[["coef"]])),
    ci_lower_95 = exp(unname(coefficient[["coef"]]) - 1.96 * standard_error),
    ci_upper_95 = exp(unname(coefficient[["coef"]]) + 1.96 * standard_error),
    standard_error_log_hr = standard_error,
    z_value = unname(coefficient[["z"]]),
    p_value = unname(coefficient[["Pr(>|z|)"]]),
    concordance = unname(fit_summary$concordance[[1L]]),
    ph_score_p_value = zph[score_row, "p"],
    ph_global_p_value = zph[global_row, "p"],
    ph_score_violation_0_05 = zph[score_row, "p"] < 0.05,
    ph_global_violation_0_05 = zph[global_row, "p"] < 0.05,
    endpoint_role = "secondary_not_proxy_hard_gate",
    optimized_cutpoint_used = FALSE
  )
}

message("[6/10] 完成单队列与GSE536家族分层合并生存分析及PH检查")
survival_data <- patient_scores[dataset %chin% c("GSE53622", "GSE53624")]
expected_survival_contract <- data.table(
  dataset = c("GSE53622", "GSE53624"),
  n_patients = c(60L, 119L),
  n_events = c(33L, 73L),
  expected_patient_id_set_sha256 = c(
    feature_list_sha256(sort(as.character(clinical_53622$patient_id))),
    feature_list_sha256(sort(as.character(clinical_53624$patient_id)))
  )
)
for (contract_index in seq_len(nrow(expected_survival_contract))) {
  contract <- expected_survival_contract[contract_index]
  cohort_value <- contract$dataset
  cohort_data <- survival_data[survival_data$dataset == cohort_value]
  fail_if(
    nrow(cohort_data) != contract$n_patients ||
      uniqueN(cohort_data$patient_id) != contract$n_patients ||
      sum(cohort_data$os_event) != contract$n_events ||
      anyDuplicated(cohort_data$gsm) ||
      feature_list_sha256(sort(as.character(cohort_data$patient_id))) !=
        contract$expected_patient_id_set_sha256,
    paste("单队列生存输入与clinical map合同不一致：", cohort_value)
  )
}
fail_if(
  nrow(survival_data) != 179L || uniqueN(survival_data$patient_id) != 179L ||
    sum(survival_data$os_event) != 106L ||
    !setequal(survival_data$patient_id, c(
      clinical_53622$patient_id, clinical_53624$patient_id
    )),
  "合并生存输入未严格复现179例/106事件/两个clinical map全集。"
)
survival_rows <- list()
survival_counter <- 0L
for (spec_index in seq_len(nrow(score_specs))) {
  spec <- score_specs[spec_index]
  for (cohort_value in c("GSE53622", "GSE53624")) {
    cohort_data <- survival_data[survival_data$dataset == cohort_value]
    for (adjusted in c(FALSE, TRUE)) {
      survival_counter <- survival_counter + 1L
      survival_rows[[survival_counter]] <- fit_cox_association(
        cohort_data, cohort_value,
        spec$factor, spec$model_variant, spec$score_column,
        adjusted = adjusted, stratified = FALSE
      )
    }
  }
  for (adjusted in c(FALSE, TRUE)) {
    survival_counter <- survival_counter + 1L
    survival_rows[[survival_counter]] <- fit_cox_association(
      survival_data, "GSE53622_GSE53624_stratified",
      spec$factor, spec$model_variant, spec$score_column,
      adjusted = adjusted, stratified = TRUE
    )
  }
}
survival_associations <- rbindlist(survival_rows)
survival_associations[, q_value := p.adjust(p_value, method = "BH")]
survival_associations[, direction := fifelse(
  hazard_ratio_per_1sd > 1, "higher_score_higher_hazard",
  "higher_score_lower_hazard"
)]
survival_associations[, interpretation_status := fcase(
  p_value <= 0.10 & !ph_score_violation_0_05, "nominal_secondary_support",
  p_value <= 0.10 & ph_score_violation_0_05, "conditional_PH_violation",
  default = "secondary_null_or_inconclusive"
)]
for (contract_index in seq_len(nrow(expected_survival_contract))) {
  contract <- expected_survival_contract[contract_index]
  cohort_rows <- survival_associations[cohort_scope == contract$dataset]
  fail_if(
    nrow(cohort_rows) != 8L ||
      any(cohort_rows$n_patients != contract$n_patients) ||
      any(cohort_rows$n_events != contract$n_events) ||
      any(cohort_rows$patient_id_set_sha256 !=
            contract$expected_patient_id_set_sha256) ||
      any(cohort_rows$stratified_by_subseries) ||
      any(grepl("strata\\(subseries\\)", cohort_rows$model_formula)),
    paste("单队列生存结果标签-数据合同失败：", contract$dataset)
  )
}
pooled_survival_rows <- survival_associations[
  cohort_scope == "GSE53622_GSE53624_stratified"
]
fail_if(
  nrow(pooled_survival_rows) != 8L ||
    any(pooled_survival_rows$n_patients != 179L) ||
    any(pooled_survival_rows$n_events != 106L) ||
    any(!pooled_survival_rows$stratified_by_subseries) ||
    any(!grepl("strata\\(subseries\\)", pooled_survival_rows$model_formula)) ||
    any(pooled_survival_rows$subseries_levels != "GSE53622,GSE53624") ||
    any(pooled_survival_rows$patient_id_set_sha256 != feature_list_sha256(
      sort(as.character(survival_data$patient_id))
    )),
  "合并生存结果未严格复现179/106且strata(subseries)。"
)
setorder(
  survival_associations, factor, model_variant, cohort_scope, model_type
)

cox_cindex <- function(fit) {
  as.numeric(summary(fit)$concordance[[1L]])
}

cox_cindex_on_data <- function(fit, data) {
  linear_predictor <- as.numeric(stats::predict(
    fit, newdata = data, type = "lp", reference = "zero"
  ))
  fail_if(length(linear_predictor) != nrow(data) ||
            any(!is.finite(linear_predictor)),
          "bootstrap Cox模型在评估数据上的线性预测值无效。")
  evaluation_data <- data.frame(
    time = data$time,
    event = data$event,
    subseries = data$subseries,
    linear_predictor = linear_predictor
  )
  as.numeric(concordance(
    Surv(time, event) ~ linear_predictor + strata(subseries),
    data = evaluation_data, reverse = TRUE
  )$concordance)
}

fit_ecms_increment <- function(
    data, factor_name, model_variant, score_column, bootstrap_seed) {
  required <- c(
    "survival_time_months", "os_event", "age", "sex", "tnm_stage_numeric",
    "subseries", "ecms_label", score_column
  )
  d <- copy(data[complete.cases(data[, ..required])])
  model_data <- data.frame(
    time = d$survival_time_months,
    event = d$os_event,
    score = d[[score_column]],
    age = d$age,
    sex = factor(d$sex, levels = c("female", "male")),
    stage = d$tnm_stage_numeric,
    subseries = factor(d$subseries),
    ecms = factor(d$ecms_label, levels = ecms_levels)
  )
  reduced_formula <-
    Surv(time, event) ~ age + sex + stage + ecms + strata(subseries)
  full_formula <-
    Surv(time, event) ~ age + sex + stage + ecms + score + strata(subseries)
  reduced <- coxph(
    reduced_formula, data = model_data, ties = "efron",
    x = TRUE, y = TRUE, model = TRUE
  )
  full <- coxph(
    full_formula, data = model_data, ties = "efron",
    x = TRUE, y = TRUE, model = TRUE
  )
  lrt_chisq <- 2 * (full$loglik[[2L]] - reduced$loglik[[2L]])
  lrt_df <- length(coef(full)) - length(coef(reduced))
  lrt_p <- pchisq(lrt_chisq, df = lrt_df, lower.tail = FALSE)
  coefficient <- summary(full)$coefficients["score", ]
  standard_error <- unname(coefficient[["se(coef)"]])
  reduced_c <- cox_cindex(reduced)
  full_c <- cox_cindex(full)
  apparent_delta_c <- full_c - reduced_c

  set.seed(bootstrap_seed)
  bootstrap_metrics <- matrix(
    NA_real_, nrow = 300L, ncol = 3L,
    dimnames = list(
      NULL,
      c("training_delta_cindex", "original_test_delta_cindex", "optimism")
    )
  )
  cohort_indices <- split(seq_len(nrow(model_data)), model_data$subseries)
  for (bootstrap_id in seq_len(300L)) {
    sampled <- unlist(lapply(cohort_indices, function(index) {
      sample(index, length(index), replace = TRUE)
    }), use.names = FALSE)
    boot_data <- model_data[sampled, , drop = FALSE]
    bootstrap_metrics[bootstrap_id, ] <- suppressWarnings(tryCatch({
      reduced_boot <- coxph(
        reduced_formula, data = boot_data, ties = "efron",
        x = TRUE, y = TRUE
      )
      full_boot <- coxph(
        full_formula, data = boot_data, ties = "efron",
        x = TRUE, y = TRUE
      )
      training_delta <-
        cox_cindex_on_data(full_boot, boot_data) -
        cox_cindex_on_data(reduced_boot, boot_data)
      original_test_delta <-
        cox_cindex_on_data(full_boot, model_data) -
        cox_cindex_on_data(reduced_boot, model_data)
      c(
        training_delta_cindex = training_delta,
        original_test_delta_cindex = original_test_delta,
        optimism = training_delta - original_test_delta
      )
    }, error = function(condition) rep(NA_real_, 3L)))
  }
  training_delta <- bootstrap_metrics[, "training_delta_cindex"]
  original_test_delta <- bootstrap_metrics[, "original_test_delta_cindex"]
  optimism <- bootstrap_metrics[, "optimism"]
  valid_joint <- apply(bootstrap_metrics, 1L, function(value) {
    all(is.finite(value))
  })
  valid_training <- training_delta[is.finite(training_delta)]
  valid_original_test <- original_test_delta[is.finite(original_test_delta)]
  valid_optimism <- optimism[valid_joint]
  fail_if(sum(valid_joint) < 240L,
          paste("分层bootstrap成对有效次数不足240/300：", factor_name, model_variant))
  corrected_delta_distribution <- apparent_delta_c - valid_optimism
  optimism_corrected_delta_c <- apparent_delta_c - mean(valid_optimism)
  zph <- cox.zph(full, transform = "km")$table
  score_row <- grep("^score$", rownames(zph))
  global_row <- grep("^GLOBAL$", rownames(zph))
  fail_if(length(score_row) != 1L || length(global_row) != 1L,
          paste("ECMS增量cox.zph缺少score/GLOBAL：", factor_name, model_variant))
  data.table(
    cohort_scope = "GSE53622_GSE53624_stratified",
    independence_family = "GSE53625_family",
    model_variant = model_variant,
    factor = factor_name,
    n_patients = nrow(model_data),
    n_events = sum(model_data$event),
    ecms1_n = sum(model_data$ecms == "ECMS1"),
    ecms2_n = sum(model_data$ecms == "ECMS2"),
    ecms3_n = sum(model_data$ecms == "ECMS3"),
    ecms4_n = sum(model_data$ecms == "ECMS4"),
    reduced_model = paste(deparse(reduced_formula), collapse = ""),
    full_model = paste(deparse(full_formula), collapse = ""),
    lrt_chisq = lrt_chisq,
    lrt_df = lrt_df,
    lrt_p_value = lrt_p,
    score_hazard_ratio_per_1sd = exp(unname(coefficient[["coef"]])),
    score_ci_lower_95 = exp(unname(coefficient[["coef"]]) - 1.96 * standard_error),
    score_ci_upper_95 = exp(unname(coefficient[["coef"]]) + 1.96 * standard_error),
    score_p_value = unname(coefficient[["Pr(>|z|)"]]),
    reduced_aic = AIC(reduced),
    full_aic = AIC(full),
    delta_aic_full_minus_reduced = AIC(full) - AIC(reduced),
    reduced_cindex = reduced_c,
    full_cindex = full_c,
    delta_cindex = apparent_delta_c,
    apparent_delta_cindex = apparent_delta_c,
    bootstrap_replicates_requested = 300L,
    bootstrap_replicates_valid = sum(valid_joint),
    bootstrap_training_replicates_valid = length(valid_training),
    bootstrap_original_test_replicates_valid = length(valid_original_test),
    bootstrap_optimism_replicates_valid = length(valid_optimism),
    bootstrap_training_delta_cindex_median = median(valid_training),
    bootstrap_training_delta_cindex_ci_lower_95 =
      quantile(valid_training, 0.025, names = FALSE),
    bootstrap_training_delta_cindex_ci_upper_95 =
      quantile(valid_training, 0.975, names = FALSE),
    bootstrap_original_test_delta_cindex_median = median(valid_original_test),
    bootstrap_original_test_delta_cindex_ci_lower_95 =
      quantile(valid_original_test, 0.025, names = FALSE),
    bootstrap_original_test_delta_cindex_ci_upper_95 =
      quantile(valid_original_test, 0.975, names = FALSE),
    bootstrap_optimism_mean = mean(valid_optimism),
    bootstrap_optimism_median = median(valid_optimism),
    optimism_corrected_delta_cindex = optimism_corrected_delta_c,
    optimism_corrected_delta_cindex_ci_lower_95 =
      quantile(corrected_delta_distribution, 0.025, names = FALSE),
    optimism_corrected_delta_cindex_ci_upper_95 =
      quantile(corrected_delta_distribution, 0.975, names = FALSE),
    ph_score_p_value = zph[score_row, "p"],
    ph_global_p_value = zph[global_row, "p"],
    endpoint_role = "conditional_increment_over_shared_314_gene_ECMS_context",
    independent_omics_validation = FALSE
  )
}

message("[7/10] 完成ECMS增量模型与300次分层bootstrap乐观偏差校正ΔC-index")
ecms_increment_rows <- lapply(seq_len(nrow(score_specs)), function(index) {
  spec <- score_specs[index]
  fit_ecms_increment(
    survival_data, spec$factor, spec$model_variant, spec$score_column,
    bootstrap_seed = 202607160L + index * 1000L
  )
})
ecms_increment <- rbindlist(ecms_increment_rows)
ecms_increment[, lrt_q_value := p.adjust(lrt_p_value, method = "BH")]
ecms_increment[, score_q_value := p.adjust(score_p_value, method = "BH")]
ecms_increment[, increment_status := fcase(
  lrt_p_value <= 0.10 & optimism_corrected_delta_cindex_ci_lower_95 > 0,
  "conditional_increment_supported",
  lrt_p_value <= 0.10,
  "nominal_increment_CI_includes_zero",
  optimism_corrected_delta_cindex_ci_lower_95 > 0,
  "Cindex_increment_without_LRT_support",
  default = "no_incremental_support"
)]
setorder(ecms_increment, factor, model_variant)

rank_biserial <- function(score, response_binary) {
  positive <- score[response_binary == 1L]
  negative <- score[response_binary == 0L]
  comparisons <- outer(positive, negative, "-")
  probability_superiority <-
    (sum(comparisons > 0) + 0.5 * sum(comparisons == 0)) / length(comparisons)
  2 * probability_superiority - 1
}

fit_response <- function(data, factor_name, model_variant, score_column) {
  d <- copy(data[complete.cases(response_binary, get(score_column))])
  model_data <- data.frame(
    response = d$response_binary,
    score = d[[score_column]]
  )
  wilcox <- wilcox.test(
    score ~ response, data = model_data,
    alternative = "two.sided", exact = FALSE, conf.int = FALSE
  )
  firth <- logistf::logistf(
    response ~ score, data = model_data, firth = TRUE, pl = TRUE
  )
  fail_if(!all(c("score") %in% names(firth$coefficients)) ||
            any(!is.finite(c(
              firth$coefficients[["score"]], firth$ci.lower[["score"]],
              firth$ci.upper[["score"]], firth$prob[["score"]]
            ))),
          paste("GSE45670 Firth logistic输出异常：", factor_name, model_variant))
  data.table(
    cohort = "GSE45670",
    model_variant = model_variant,
    factor = factor_name,
    n_total = nrow(model_data),
    n_pathological_complete_response = sum(model_data$response == 1L),
    n_not_pathological_complete_response = sum(model_data$response == 0L),
    median_score_pcr = median(model_data$score[model_data$response == 1L]),
    median_score_non_pcr = median(model_data$score[model_data$response == 0L]),
    wilcoxon_w = unname(wilcox$statistic),
    wilcoxon_p_value = wilcox$p.value,
    rank_biserial_pcr_higher_positive = rank_biserial(
      model_data$score, model_data$response
    ),
    firth_log_odds_per_1sd = firth$coefficients[["score"]],
    firth_odds_ratio_per_1sd = exp(firth$coefficients[["score"]]),
    firth_or_ci_lower_95 = exp(firth$ci.lower[["score"]]),
    firth_or_ci_upper_95 = exp(firth$ci.upper[["score"]]),
    firth_profile_likelihood_p_value = firth$prob[["score"]],
    firth_convergence_max_abs_score = unname(firth$conv[["max abs score"]]),
    endpoint_role = "exploratory_only",
    treatment_predictor_claim = FALSE,
    source_response_field = "pathological response to crp"
  )
}

message("[8/10] 完成GSE45670探索性Wilcoxon、rank-biserial与Firth分析")
response_associations <- rbindlist(lapply(
  seq_len(nrow(score_specs)), function(index) {
    spec <- score_specs[index]
    fit_response(
      patient_scores[dataset == "GSE45670"],
      spec$factor, spec$model_variant, spec$score_column
    )
  }
))
response_associations[, wilcoxon_q_value := p.adjust(
  wilcoxon_p_value, method = "BH"
)]
response_associations[, firth_q_value := p.adjust(
  firth_profile_likelihood_p_value, method = "BH"
)]
response_associations[, response_status := fcase(
  firth_profile_likelihood_p_value <= 0.10 & wilcoxon_p_value <= 0.10,
  "exploratory_concordant_signal",
  firth_profile_likelihood_p_value <= 0.10 | wilcoxon_p_value <= 0.10,
  "exploratory_single_method_signal",
  default = "exploratory_null_or_inconclusive"
)]
setorder(response_associations, factor, model_variant)

for (factor_name in target_factors) {
  main_rows <- survival_associations[
    factor == factor_name & model_variant == "ridge314_primary" &
      model_type == "clinical_adjusted"
  ]
  row_22 <- main_rows[cohort_scope == "GSE53622"]
  row_24 <- main_rows[cohort_scope == "GSE53624"]
  pooled <- main_rows[cohort_scope == "GSE53622_GSE53624_stratified"]
  fail_if(nrow(row_22) != 1L || nrow(row_24) != 1L || nrow(pooled) != 1L,
          paste("生存决策缺少主ridge调整模型：", factor_name))
  same_direction <- sign(log(row_22$hazard_ratio_per_1sd)) ==
    sign(log(row_24$hazard_ratio_per_1sd))
  survival_status <- if (same_direction && pooled$p_value <= 0.10 &&
                           !pooled$ph_score_violation_0_05) {
    "SECONDARY_SUPPORT"
  } else if (same_direction) {
    "CONDITIONAL_SAME_DIRECTION"
  } else {
    "HETEROGENEOUS_OR_NULL"
  }
  add_decision(
    paste0("SURVIVAL_", toupper(factor_name)), "external_secondary_endpoint",
    factor = factor_name, hard_gate = FALSE,
    metric = "adjusted HR direction in GSE53622/GSE53624 and pooled stratified p",
    observed = sprintf(
      "HR22=%.3f;HR24=%.3f;pooledHR=%.3f;p=%.3g;PHscoreP=%.3g",
      row_22$hazard_ratio_per_1sd, row_24$hazard_ratio_per_1sd,
      pooled$hazard_ratio_per_1sd, pooled$p_value, pooled$ph_score_p_value
    ),
    rule = paste(
      "support if both directions agree and pooled p<=0.10 without score PH failure;",
      "same direction otherwise conditional; survival never a proxy hard veto"
    ),
    status = survival_status,
    counts_toward_external_support = TRUE,
    interpretation = if (survival_status == "SECONDARY_SUPPORT") {
      "外部研究家族对次要生存端点提供方向一致支持。"
    } else if (survival_status == "CONDITIONAL_SAME_DIRECTION") {
      "方向一致但统计区间/功效不足，作为条件支持保留。"
    } else {
      "生存关联阴性或子系列异质；保留为端点边界。"
    },
    boundary = "TCGA内部OS原本不显著；外部OS阴性不否定状态proxy。"
  )

  increment_main <- ecms_increment[
    factor == factor_name & model_variant == "ridge314_primary"
  ]
  add_decision(
    paste0("ECMS_INCREMENT_", toupper(factor_name)), "conditional_increment",
    factor = factor_name, hard_gate = FALSE,
    metric = paste(
      "LRT and 300x stratified bootstrap apparent/training/original-test",
      "and optimism-corrected delta C-index"
    ),
    observed = sprintf(
      paste0(
        "LRTp=%.3g;apparent_deltaC=%.4f;optimism_corrected_deltaC=%.4f;",
        "corrected95%%CI=[%.4f,%.4f];paired_valid=%d/300"
      ),
      increment_main$lrt_p_value, increment_main$delta_cindex,
      increment_main$optimism_corrected_delta_cindex,
      increment_main$optimism_corrected_delta_cindex_ci_lower_95,
      increment_main$optimism_corrected_delta_cindex_ci_upper_95,
      increment_main$bootstrap_optimism_replicates_valid
    ),
    rule = "conditional evidence only; ECMS and proxy share the same 314-gene RNA representation",
    status = increment_main$increment_status,
    counts_toward_external_support = TRUE,
    interpretation = "评估proxy在临床变量和ECMS表达分型之上的条件增量。",
    boundary = "不是独立组学验证，也不等于临床模型可用。"
  )

  response_main <- response_associations[
    factor == factor_name & model_variant == "ridge314_primary"
  ]
  add_decision(
    paste0("RESPONSE_", toupper(factor_name)), "exploratory_response",
    factor = factor_name, hard_gate = FALSE,
    metric = "GSE45670 rank-biserial and Firth OR",
    observed = sprintf(
      "rank_biserial=%.3f;FirthOR=%.3f;profileP=%.3g",
      response_main$rank_biserial_pcr_higher_positive,
      response_main$firth_odds_ratio_per_1sd,
      response_main$firth_profile_likelihood_p_value
    ),
    rule = "exploratory only; no treatment-predictor claim",
    status = response_main$response_status,
    counts_toward_external_support = FALSE,
    interpretation = "小样本疗效语境提供探索性校准。",
    boundary = "28例、11个pCR；不进入总体GO硬门禁。"
  )
}

add_decision(
  "OVERALL_FIRST_STAGE_EXTERNAL_VALIDATION", "overall", hard_gate = FALSE,
  metric = "input/mapping hard gates + both ridge proxy fidelity gates",
  observed = "all hard identity gates PASS; Factor1/Factor3 ridge proxy CV at least CONDITIONAL_GO",
  rule = "complete first-stage external RNA-proxy calibration; survival/response are endpoint-specific secondary evidence",
  status = "GO_FIRST_STAGE_EXTERNAL_RNA_PROXY_COMPLETED",
  counts_toward_external_support = TRUE,
  interpretation = paste(
    "可以升级为独立表达队列中的连续状态RNA代理校准；",
    "阴性生存或疗效结果不被删除。"
  ),
  boundary = paste(
    "不能写成完整MOFA因子外部复现；9条突变/CNV事件-状态边仍是",
    "来源队列内部发现。"
  )
)

decisions <- rbindlist(decision_rows, use.names = TRUE, fill = TRUE)
fail_if(any(decisions$hard_gate & decisions$status != "PASS" &
              !decisions$status %chin% c("GO", "CONDITIONAL_GO")),
        "外部验证decision存在未通过硬门禁。")

message("[9/10] 写入正式结果、决策表与中文边界摘要")
stage_fwrite(patient_scores, "escc_external_state_patient_scores.tsv")
stage_fwrite(
  survival_associations, "escc_external_state_survival_associations.tsv"
)
stage_fwrite(ecms_increment, "escc_external_state_ecms_increment.tsv")
stage_fwrite(
  response_associations, "escc_external_state_response_associations.tsv"
)
stage_fwrite(decisions, "escc_external_validation_decision.tsv")

cv_factor1 <- internal_cv[
  factor == "Factor1" & model_variant == "ridge314_primary" &
    record_type == "outer_repeat_summary"
]
cv_factor3 <- internal_cv[
  factor == "Factor3" & model_variant == "ridge314_primary" &
    record_type == "outer_repeat_summary"
]
survival_main <- survival_associations[
  model_variant == "ridge314_primary" & model_type == "clinical_adjusted" &
    cohort_scope == "GSE53622_GSE53624_stratified"
]
response_main <- response_associations[model_variant == "ridge314_primary"]
summary_lines <- c(
  "# ESCC 连续状态第一阶段外部验证摘要",
  "",
  "## 完成范围",
  "",
  "- 输入字节、CATALOG、DATASET、MANIFEST、上游MOFA artifact和固定ECMS模型SHA256均已通过硬核验。",
  "- 在读取外部生存与疗效标签之前，已冻结Factor1/Factor3的314基因ridge proxy、lambda、系数、方向、预处理和51基因固定权重敏感性定义。",
  "- GSE53622与GSE53624分别重做队列内逐基因Z标准化；GSE53625只表示两者的超级系列，不计为第三队列。",
  "- 未使用最佳截点；生存以连续每1 SD分数进入Cox模型。",
  "",
  "## Proxy fidelity",
  "",
  sprintf(
    "- Factor1：20次外层5折CV的Spearman中位数 %.3f，范围 %.3f–%.3f，门禁 `%s`。",
    cv_factor1$spearman_rho, cv_factor1$spearman_min,
    cv_factor1$spearman_max, cv_factor1$gate_status
  ),
  sprintf(
    "- Factor3：20次外层5折CV的Spearman中位数 %.3f，范围 %.3f–%.3f，门禁 `%s`。",
    cv_factor3$spearman_rho, cv_factor3$spearman_min,
    cv_factor3$spearman_max, cv_factor3$gate_status
  ),
  "- 51基因原始MOFA RNA权重分数只作透明敏感性；它不替代314基因主proxy。",
  "",
  "## 外部生存（次要端点）",
  "",
  sprintf(
    "- Factor1：179例分层临床调整HR %.3f（95%%CI %.3f–%.3f），p=%.3g，q=%.3g；score PH p=%.3g。",
    survival_main[factor == "Factor1", hazard_ratio_per_1sd],
    survival_main[factor == "Factor1", ci_lower_95],
    survival_main[factor == "Factor1", ci_upper_95],
    survival_main[factor == "Factor1", p_value],
    survival_main[factor == "Factor1", q_value],
    survival_main[factor == "Factor1", ph_score_p_value]
  ),
  sprintf(
    "- Factor3：179例分层临床调整HR %.3f（95%%CI %.3f–%.3f），p=%.3g，q=%.3g；score PH p=%.3g。",
    survival_main[factor == "Factor3", hazard_ratio_per_1sd],
    survival_main[factor == "Factor3", ci_lower_95],
    survival_main[factor == "Factor3", ci_upper_95],
    survival_main[factor == "Factor3", p_value],
    survival_main[factor == "Factor3", q_value],
    survival_main[factor == "Factor3", ph_score_p_value]
  ),
  "- GSE53622与GSE53624来自同一研究家族；分队列效应和分层合并效应均保留，不能夸大成两个完全独立研究。",
  "- 生存阴性、置信区间跨1或PH偏离均属于端点边界，不构成proxy本身的机械失败。",
  "",
  "## ECMS条件增量",
  "",
  sprintf(
    paste0(
      "- 锁定作者179矩阵的确定性ECMS计数为 %s；vote并列 %d 例（%s）。",
      "并列按固定类水平顺序取首位，规则完全不读取结局。"
    ),
    paste(ecms_label_counts, collapse = "/"),
    length(gse53625_tie_patients),
    ifelse(length(gse53625_tie_patients),
           paste(gse53625_tie_patients, collapse = ","), "无")
  ),
  paste(
    "- 在179例分层模型中比较临床+ECMS与临床+ECMS+proxy；每次bootstrap同时记录",
    "training ΔC与原始样本test ΔC，并用两者差值进行乐观偏差校正。"
  ),
  sprintf(
    paste0(
      "- Factor1主proxy：apparent ΔC %.4f，optimism-corrected ΔC %.4f",
      "（95%% bootstrap区间 %.4f–%.4f；成对有效 %d/300）。"
    ),
    ecms_increment[factor == "Factor1" & model_variant == "ridge314_primary",
                   apparent_delta_cindex],
    ecms_increment[factor == "Factor1" & model_variant == "ridge314_primary",
                   optimism_corrected_delta_cindex],
    ecms_increment[factor == "Factor1" & model_variant == "ridge314_primary",
                   optimism_corrected_delta_cindex_ci_lower_95],
    ecms_increment[factor == "Factor1" & model_variant == "ridge314_primary",
                   optimism_corrected_delta_cindex_ci_upper_95],
    ecms_increment[factor == "Factor1" & model_variant == "ridge314_primary",
                   bootstrap_optimism_replicates_valid]
  ),
  sprintf(
    paste0(
      "- Factor3主proxy：apparent ΔC %.4f，optimism-corrected ΔC %.4f",
      "（95%% bootstrap区间 %.4f–%.4f；成对有效 %d/300）。"
    ),
    ecms_increment[factor == "Factor3" & model_variant == "ridge314_primary",
                   apparent_delta_cindex],
    ecms_increment[factor == "Factor3" & model_variant == "ridge314_primary",
                   optimism_corrected_delta_cindex],
    ecms_increment[factor == "Factor3" & model_variant == "ridge314_primary",
                   optimism_corrected_delta_cindex_ci_lower_95],
    ecms_increment[factor == "Factor3" & model_variant == "ridge314_primary",
                   optimism_corrected_delta_cindex_ci_upper_95],
    ecms_increment[factor == "Factor3" & model_variant == "ridge314_primary",
                   bootstrap_optimism_replicates_valid]
  ),
  "- ECMS与proxy共享同一314基因RNA表示；该结果最多是条件增量证据，不是独立组学验证。",
  "",
  "## GSE45670探索性疗效关联",
  "",
  sprintf(
    "- Factor1：rank-biserial %.3f，Firth OR/1SD %.3f（95%%CI %.3f–%.3f），profile p=%.3g。",
    response_main[factor == "Factor1", rank_biserial_pcr_higher_positive],
    response_main[factor == "Factor1", firth_odds_ratio_per_1sd],
    response_main[factor == "Factor1", firth_or_ci_lower_95],
    response_main[factor == "Factor1", firth_or_ci_upper_95],
    response_main[factor == "Factor1", firth_profile_likelihood_p_value]
  ),
  sprintf(
    "- Factor3：rank-biserial %.3f，Firth OR/1SD %.3f（95%%CI %.3f–%.3f），profile p=%.3g。",
    response_main[factor == "Factor3", rank_biserial_pcr_higher_positive],
    response_main[factor == "Factor3", firth_odds_ratio_per_1sd],
    response_main[factor == "Factor3", firth_or_ci_lower_95],
    response_main[factor == "Factor3", firth_or_ci_upper_95],
    response_main[factor == "Factor3", firth_profile_likelihood_p_value]
  ),
  "- GSE45670仅28例（17 non-pCR、11 pCR），因此只作探索性语境校准，不宣称疗效预测器。",
  "",
  "## 总体判定与证据上限",
  "",
  "- 总体状态：`GO_FIRST_STAGE_EXTERNAL_RNA_PROXY_COMPLETED`。",
  "- 可支持：Factor1/Factor3的冻结RNA代理已经完成独立表达队列的患者级临床/疗效校准。",
  "- 不可支持：完整多组学MOFA因子已经在外部队列重建，或9条突变/CNV事件—状态边已经独立复现。",
  "- 本研究只使用公共数据，未新采集人体样本、未开展动物实验或湿实验。"
)
stage_write_lines(summary_lines, "escc_external_validation_summary.md")

input_hashes <- c(
  paste0("catalog=", sha256_file(catalog_path)),
  paste0("heterogeneity_manifest=", sha256_file(heterogeneity_manifest_path)),
  paste0("mofa_model=", sha256_file(mofa_model_path)),
  paste0("factor_scores=", sha256_file(factor_scores_path)),
  paste0("ecms_model=", locked_model_sha256),
  unlist(lapply(hub_datasets, function(dataset) {
    c(
      paste0(dataset$key, "_manifest=", dataset$manifest_sha256),
      paste0(dataset$key, "_DATASET=", dataset$dataset_md_sha256)
    )
  }), use.names = FALSE)
)
input_hash_bundle <- paste(input_hashes, collapse = ";")

artifact_paths <- file.path(stage_dir, formal_filenames)
fail_if(any(!file_exists(artifact_paths)),
        paste("外部验证stage正式artifact缺失：",
              paste(formal_filenames[!file_exists(artifact_paths)], collapse = ";")))
artifact_manifest <- data.table(
  artifact = formal_filenames,
  relative_path = file.path("results", formal_filenames),
  file_size_bytes = as.numeric(file_info(artifact_paths)$size),
  sha256 = vapply(artifact_paths, sha256_file, character(1)),
  generated_date = as.character(Sys.Date()),
  generation_script = "scripts/25_validate_external_continuous_states.R",
  execution_script_sha256 = execution_script_sha256,
  input_sha256 = input_hash_bundle,
  proxy_model_sha256 = proxy_model_sha256,
  model_definition_sha256 = model_definition_sha256,
  outcome_parsed_after_model_freeze = TRUE,
  figure6_generated = FALSE,
  status = "verified"
)
stage_fwrite(artifact_manifest, manifest_filename)

message("[10/10] 原子发布10个artifact，manifest最后发布并清理stage")
for (filename in formal_filenames) {
  atomic_publish_file(
    file.path(stage_dir, filename), file.path(results_dir, filename)
  )
}
for (index in seq_len(nrow(artifact_manifest))) {
  published <- file.path(project_root, artifact_manifest$relative_path[[index]])
  fail_if(!file_exists(published) ||
            as.numeric(file_info(published)$size) !=
              artifact_manifest$file_size_bytes[[index]] ||
            sha256_file(published) != artifact_manifest$sha256[[index]],
          paste("发布artifact回读失败：", published))
}
atomic_publish_file(
  file.path(stage_dir, manifest_filename),
  file.path(results_dir, manifest_filename)
)
published_manifest <- fread(
  file.path(results_dir, manifest_filename),
  colClasses = "character", na.strings = NULL, showProgress = FALSE
)
fail_if(nrow(published_manifest) != length(formal_filenames) ||
          !identical(published_manifest$relative_path,
                     artifact_manifest$relative_path) ||
          any(published_manifest$sha256 != artifact_manifest$sha256),
        "外部验证发布manifest与stage冻结版不一致。")

if (dir_exists(stage_dir)) dir_delete(stage_dir)
message(
  "完成：Factor1/Factor3 314基因ridge proxy、51基因敏感性、20×5重复外层CV、",
  "GSE53622/24生存与ECMS增量、GSE45670探索性响应已发布；",
  "Figure6未生成，PROJECT_INDEX/data索引/ResearchDataHub未修改。"
)
