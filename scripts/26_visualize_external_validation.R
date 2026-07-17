#!/usr/bin/env Rscript

# ESCC 连续状态外部校准 Figure 6（R-only）。
#
# 证据边界：
# 1) 只读取 script 25 已发布的正式结果表，不重拟合 proxy、Cox、
#    ECMS 或 pCR 模型，不从旧图或摘要反推统计数值。
# 2) panel d 的主效应必须是表中真正的 bootstrap 乐观偏差校正
#    delta C-index；表观 delta C 只用空心菱形作对照。
# 3) GSE53622/GSE53624 是同一 GSE53625 研究家族的两个无患者
#    重叠子系列；GSE53625 不作第三个生存队列。
# 4) Figure 6 评估冻结 RNA proxy 的可计算性与临床语境校准，不复现
#    精确基因组事件—因子边，也不声称新的临床预测器。
# 5) 图件先在 _work/intermediate/ 临时渲染；4 个格式均通过结构、
#    尺寸、字体和 SHA256 检查后才发布，manifest 最后发布。

args <- commandArgs(trailingOnly = TRUE)
fields_only <- "--fields-only" %in% args
validate_only <- "--validate-only" %in% args
finalize_args <- grep("^--finalize-visual-qa=", args, value = TRUE)
if (length(finalize_args) > 1L) {
  stop("--finalize-visual-qa 只能指定一次。", call. = FALSE)
}
finalize_visual_qa <- length(finalize_args) == 1L
finalize_visual_qa_input <- if (finalize_visual_qa) {
  sub("^--finalize-visual-qa=", "", finalize_args[[1L]])
} else {
  NULL
}
unknown_args <- setdiff(
  args, c("--fields-only", "--validate-only", finalize_args)
)
if (length(unknown_args)) {
  stop("未知参数：", paste(unknown_args, collapse = ", "), call. = FALSE)
}
if (sum(c(fields_only, validate_only, finalize_visual_qa)) > 1L) {
  stop(
    "--fields-only、--validate-only 与 --finalize-visual-qa 不能同时使用。",
    call. = FALSE
  )
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) {
  stop("无法唯一定位当前脚本。", call. = FALSE)
}
script_path <- normalizePath(
  sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE
)
project_root <- dirname(dirname(script_path))
results_dir <- file.path(project_root, "results")
figures_dir <- file.path(project_root, "figures")
work_intermediate_dir <- file.path(project_root, "_work", "intermediate")
work_checks_dir <- file.path(project_root, "_work", "checks")

figure_id <- "Figure6"
figure_stem <- "escc_multiomics_figure6_external_validation"
figure_formats <- c("svg", "pdf", "tiff", "png")
figure_relative_paths <- file.path(
  "figures", paste0(figure_stem, ".", figure_formats)
)
figure_manifest_relative_path <- file.path(
  "results", "escc_external_validation_figure_artifact_manifest.tsv"
)
figure_contract_relative_path <- file.path(
  "_work", "checks", "escc_multiomics_figure6_contract_20260716.md"
)

if (fields_only) {
  cat(paste(figure_relative_paths, collapse = "\n"), "\n", sep = "")
  cat(figure_manifest_relative_path, "\n", sep = "")
  cat(
    "backend\tR-only\n",
    "size_mm\t183x215\n",
    "raster_contract\tTIFF and PNG at 600 dpi\n",
    "manifest_publish_order\t4 figure artifacts first; manifest last\n",
    "visual_qa_finalize\t--finalize-visual-qa=<qa.md>; requires Figure6: PASS\n",
    "panel_d_primary_statistic\toptimism_corrected_delta_cindex with corrected 95% CI\n",
    sep = ""
  )
  quit(save = "no", status = 0L)
}

required_packages <- c(
  "data.table", "digest", "fs", "ggplot2", "patchwork", "svglite",
  "ragg", "rsvg", "pdftools", "scales", "systemfonts"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "R-only Figure 6 依赖缺失：",
    paste(missing_packages, collapse = ", "),
    "。禁止切换 Python 代替渲染。",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(fs)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE, scipen = 999)

fail_if <- function(condition, message) {
  if (!is.logical(condition) || length(condition) != 1L || is.na(condition)) {
    stop(paste0(message, "（门禁为 NA 或非单一逻辑值）"), call. = FALSE)
  }
  if (isTRUE(condition)) stop(message, call. = FALSE)
}

require_columns <- function(object, columns, object_name) {
  missing <- setdiff(columns, names(object))
  fail_if(length(missing) > 0L, paste0(
    object_name, " 缺少字段：", paste(missing, collapse = ", ")
  ))
}

read_tsv <- function(path) {
  fail_if(!file_exists(path), paste("正式源表缺失：", path))
  fread(path, na.strings = c("", "NA"), showProgress = FALSE)
}

sha256_file <- function(path) {
  fail_if(!file_exists(path), paste("待校验文件缺失：", path))
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

as_logical_strict <- function(value, label, allow_na = FALSE) {
  input_na <- is.na(value)
  output <- rep(NA, length(value))
  if (is.logical(value)) {
    output <- value
  } else if (is.numeric(value)) {
    output[!input_na & value == 1] <- TRUE
    output[!input_na & value == 0] <- FALSE
  } else {
    normalized <- tolower(trimws(as.character(value)))
    output[!input_na & normalized %chin% c("true", "t", "1", "yes")] <- TRUE
    output[!input_na & normalized %chin% c("false", "f", "0", "no")] <- FALSE
  }
  fail_if(
    any(!input_na & is.na(output)) || (!allow_na && anyNA(output)),
    paste0(label, " 含不可识别值或不允许的 NA。")
  )
  output
}

format_p <- function(value) {
  output <- rep("NA", length(value))
  finite <- is.finite(value)
  output[finite & value >= 0.001] <- sprintf("%.3f", value[finite & value >= 0.001])
  output[finite & value < 0.001] <- format(
    value[finite & value < 0.001], scientific = TRUE, digits = 2
  )
  output
}

executed_script_sha256 <- sha256_file(script_path)
figure_contract_path <- file.path(project_root, figure_contract_relative_path)
fail_if(!file_exists(figure_contract_path), paste(
  "Figure 6 契约缺失：", figure_contract_path
))
figure_contract_sha256 <- sha256_file(figure_contract_path)

source_manifest_path <- file.path(
  results_dir, "escc_external_validation_artifact_manifest.tsv"
)
source_table_names <- c(
  internal_cv = "escc_external_state_internal_cv.tsv",
  oof_predictions = "escc_external_state_oof_predictions.tsv",
  patient_scores = "escc_external_state_patient_scores.tsv",
  survival_associations = "escc_external_state_survival_associations.tsv",
  ecms_increment = "escc_external_state_ecms_increment.tsv",
  response_associations = "escc_external_state_response_associations.tsv",
  validation_decision = "escc_external_validation_decision.tsv"
)
source_table_paths <- setNames(
  file.path(results_dir, unname(source_table_names)), names(source_table_names)
)

verify_source_family <- function() {
  fail_if(!file_exists(source_manifest_path), paste(
    "script 25 正式 manifest 缺失：", source_manifest_path
  ))
  manifest <- read_tsv(source_manifest_path)
  require_columns(
    manifest,
    c(
      "artifact", "relative_path", "file_size_bytes", "sha256",
      "generation_script", "execution_script_sha256",
      "outcome_parsed_after_model_freeze", "figure6_generated", "status"
    ),
    "external validation artifact manifest"
  )
  expected_artifacts <- c(
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
  fail_if(
    nrow(manifest) != 10L || uniqueN(manifest$artifact) != 10L ||
      !setequal(manifest$artifact, expected_artifacts) ||
      any(manifest$relative_path != file.path("results", manifest$artifact)) ||
      any(manifest$status != "verified") ||
      any(manifest$generation_script !=
            "scripts/25_validate_external_continuous_states.R"),
    "script 25 manifest 不是完整 10-artifact verified family。"
  )
  manifest[, outcome_parsed_after_model_freeze := as_logical_strict(
    outcome_parsed_after_model_freeze, "outcome parsed after model freeze"
  )]
  manifest[, figure6_generated := as_logical_strict(
    figure6_generated, "source manifest figure6 generated"
  )]
  fail_if(
    any(!manifest$outcome_parsed_after_model_freeze) ||
      any(manifest$figure6_generated),
    paste(
      "script 25 必须在结局解析前冻结模型，且不能冒充",
      "Figure 6 生成脚本。"
    )
  )
  artifact_paths <- file.path(project_root, manifest$relative_path)
  fail_if(any(!file_exists(artifact_paths)), "script 25 manifest 登记文件不完整。")
  observed_size <- as.numeric(file_info(artifact_paths)$size)
  observed_sha <- vapply(artifact_paths, sha256_file, character(1))
  fail_if(
    any(observed_size != as.numeric(manifest$file_size_bytes)) ||
      any(observed_sha != manifest$sha256),
    "script 25 artifact 当前大小/SHA256 与 manifest 不一致。"
  )
  script25_path <- file.path(
    project_root, "scripts", "25_validate_external_continuous_states.R"
  )
  fail_if(!file_exists(script25_path), "script 25 源脚本缺失。")
  manifest_script_sha <- unique(as.character(manifest$execution_script_sha256))
  fail_if(
    length(manifest_script_sha) != 1L ||
      manifest_script_sha != sha256_file(script25_path),
    "script 25 当前源码 SHA256 与已发布 manifest 不一致。"
  )
  manifest
}

source_manifest <- verify_source_family()
source_manifest_sha256 <- sha256_file(source_manifest_path)

inputs <- lapply(source_table_paths, read_tsv)

schema_contract <- list(
  internal_cv = c(
    "record_type", "model_variant", "factor", "repeat_id", "n",
    "spearman_rho", "pearson_r", "rmse", "median_absolute_error",
    "spearman_min", "spearman_max", "positive_orientation",
    "outer_folds", "inner_folds", "selection_rule", "summary_statistic",
    "gate_status", "feature_count", "rna_squared_weight_mass_fraction"
  ),
  oof_predictions = c(
    "model_variant", "factor", "repeat_id", "patient_id", "outer_fold",
    "observed_factor_score", "oof_predicted_factor_score", "residual",
    "inner_selected_lambda_1se", "orientation_flipped_after_outcome"
  ),
  patient_scores = c(
    "dataset", "independence_family", "subseries", "gsm", "patient_id",
    "sample_role", "ecms_label", "factor1_ridge314_primary_z",
    "factor3_ridge314_primary_z", "factor1_mofa51_sensitivity_z",
    "factor3_mofa51_sensitivity_z", "age", "sex", "tnm_stage_numeric",
    "survival_time_months", "os_event", "response_binary",
    "response_label", "external_endpoint_role"
  ),
  survival_associations = c(
    "cohort_scope", "independence_family", "model_variant", "factor",
    "model_type", "stratified_by_subseries", "n_patients", "n_events",
    "patient_id_set_sha256", "subseries_levels", "hazard_ratio_per_1sd",
    "ci_lower_95", "ci_upper_95", "p_value", "q_value",
    "ph_score_p_value", "ph_global_p_value", "ph_score_violation_0_05",
    "optimized_cutpoint_used", "interpretation_status"
  ),
  ecms_increment = c(
    "cohort_scope", "model_variant", "factor", "n_patients", "n_events",
    "lrt_p_value", "lrt_q_value", "delta_aic_full_minus_reduced",
    "apparent_delta_cindex", "bootstrap_replicates_requested",
    "bootstrap_optimism_replicates_valid",
    "bootstrap_training_delta_cindex_median",
    "bootstrap_original_test_delta_cindex_median",
    "bootstrap_optimism_mean", "optimism_corrected_delta_cindex",
    "optimism_corrected_delta_cindex_ci_lower_95",
    "optimism_corrected_delta_cindex_ci_upper_95",
    "independent_omics_validation", "increment_status"
  ),
  response_associations = c(
    "cohort", "model_variant", "factor", "n_total",
    "n_pathological_complete_response",
    "n_not_pathological_complete_response", "median_score_pcr",
    "median_score_non_pcr", "wilcoxon_p_value", "wilcoxon_q_value",
    "rank_biserial_pcr_higher_positive", "firth_odds_ratio_per_1sd",
    "firth_or_ci_lower_95", "firth_or_ci_upper_95",
    "firth_profile_likelihood_p_value", "firth_q_value", "endpoint_role",
    "treatment_predictor_claim", "response_status"
  ),
  validation_decision = c(
    "decision_id", "decision_domain", "factor", "hard_gate", "metric",
    "observed", "threshold_or_rule", "status",
    "counts_toward_external_support", "interpretation", "boundary"
  )
)

for (input_name in names(schema_contract)) {
  require_columns(inputs[[input_name]], schema_contract[[input_name]], input_name)
}

validate_semantics <- function(x) {
  target_factors <- c("Factor1", "Factor3")

  cv <- copy(x$internal_cv)
  cv[, positive_orientation := as_logical_strict(
    positive_orientation, "internal CV positive orientation"
  )]
  primary_repeats <- cv[
    model_variant == "ridge314_primary" & record_type == "outer_repeat"
  ]
  primary_summary <- cv[
    model_variant == "ridge314_primary" &
      record_type == "outer_repeat_summary"
  ]
  sensitivity <- cv[
    model_variant == "mofa51_fixed_weight_sensitivity" &
      record_type == "sensitivity_in_sample_diagnostic"
  ]
  fail_if(
    nrow(cv) != 44L || nrow(primary_repeats) != 40L ||
      nrow(primary_summary) != 2L || nrow(sensitivity) != 2L ||
      !setequal(cv$factor, target_factors),
    "internal CV 必须为 40 条外层 repeat + 2 条汇总 + 2 条敏感性。"
  )
  repeat_grid <- primary_repeats[, .(
    repeats = uniqueN(repeat_id), rows = .N, n_unique = uniqueN(n),
    n_value = unique(n), min_fold = min(outer_folds), max_fold = max(outer_folds)
  ), by = factor]
  fail_if(
    nrow(repeat_grid) != 2L || any(repeat_grid$repeats != 20L) ||
      any(repeat_grid$rows != 20L) || any(repeat_grid$n_unique != 1L) ||
      any(repeat_grid$n_value != 78L) || any(repeat_grid$min_fold != 5L) ||
      any(repeat_grid$max_fold != 5L) ||
      any(!is.finite(primary_repeats$spearman_rho)) ||
      any(!primary_repeats$positive_orientation),
    "20×5 外层 CV 记录的次数、样本数、折数或方向失败。"
  )
  fail_if(
    any(primary_summary$n != 78L) ||
      any(!primary_summary$gate_status %chin% c("GO", "CONDITIONAL_GO")) ||
      any(!is.finite(primary_summary$spearman_rho)) ||
      any(!is.finite(primary_summary$spearman_min)) ||
      any(!is.finite(primary_summary$spearman_max)),
    "Factor1/Factor3 ridge proxy CV 汇总未通过最低条件门禁。"
  )

  oof <- copy(x$oof_predictions)
  oof[, orientation_flipped_after_outcome := as_logical_strict(
    orientation_flipped_after_outcome, "OOF outcome orientation flip"
  )]
  fail_if(
    nrow(oof) != 3120L ||
      any(oof$model_variant != "ridge314_primary") ||
      !setequal(oof$factor, target_factors) ||
      any(oof$orientation_flipped_after_outcome) ||
      any(!is.finite(oof$observed_factor_score)) ||
      any(!is.finite(oof$oof_predicted_factor_score)) ||
      any(!is.finite(oof$inner_selected_lambda_1se)),
    "OOF 预测不是完整的 2×20×78 outcome-blinded family。"
  )
  oof_grid <- oof[, .(
    patients = uniqueN(patient_id), rows = .N,
    folds = uniqueN(outer_fold)
  ), by = .(factor, repeat_id)]
  oof_patient_consistency <- oof[, .(
    observed_values = uniqueN(observed_factor_score), repeats = .N
  ), by = .(factor, patient_id)]
  fail_if(
    nrow(oof_grid) != 40L || any(oof_grid$patients != 78L) ||
      any(oof_grid$rows != 78L) || any(oof_grid$folds != 5L) ||
      nrow(oof_patient_consistency) != 156L ||
      any(oof_patient_consistency$observed_values != 1L) ||
      any(oof_patient_consistency$repeats != 20L),
    "OOF factor×repeat 网格不是完整 78 例外层留出预测。"
  )

  patient <- copy(x$patient_scores)
  expected_counts <- c(GSE45670 = 28L, GSE53622 = 60L, GSE53624 = 119L)
  observed_counts <- patient[, .N, by = dataset]
  fail_if(
    nrow(patient) != 207L || uniqueN(patient$patient_id) != 207L ||
      !setequal(observed_counts$dataset, names(expected_counts)) ||
      any(observed_counts$N != expected_counts[observed_counts$dataset]),
    "patient scores 未保持 GSE45670=28、GSE53622=60、GSE53624=119。"
  )
  score_columns <- c(
    "factor1_ridge314_primary_z", "factor3_ridge314_primary_z",
    "factor1_mofa51_sensitivity_z", "factor3_mofa51_sensitivity_z"
  )
  fail_if(
    any(!is.finite(as.matrix(patient[, ..score_columns]))),
    "外部患者分数存在 NA/Inf。"
  )
  score_z_qc <- melt(
    patient[, c("dataset", score_columns), with = FALSE],
    id.vars = "dataset", variable.name = "score", value.name = "value"
  )[, .(mean = mean(value), sd = sd(value)), by = .(dataset, score)]
  fail_if(
    any(abs(score_z_qc$mean) > 1e-8) || any(abs(score_z_qc$sd - 1) > 1e-8),
    "外部分数 *_z 未按子队列保持均值 0/SD 1。"
  )
  response_patient <- patient[dataset == "GSE45670"]
  fail_if(
    nrow(response_patient) != 28L ||
      sum(response_patient$response_binary == 1L) != 11L ||
      sum(response_patient$response_binary == 0L) != 17L ||
      anyNA(response_patient$response_binary),
    "GSE45670 必须保持 11 pCR/17 non-pCR。"
  )

  survival <- copy(x$survival_associations)
  survival[, `:=`(
    stratified_by_subseries = as_logical_strict(
      stratified_by_subseries, "survival stratified flag"
    ),
    ph_score_violation_0_05 = as_logical_strict(
      ph_score_violation_0_05, "survival PH score flag"
    ),
    optimized_cutpoint_used = as_logical_strict(
      optimized_cutpoint_used, "survival optimized cutpoint"
    )
  )]
  fail_if(
    nrow(survival) != 24L ||
      !setequal(survival$factor, target_factors) ||
      any(survival$optimized_cutpoint_used) ||
      any(!is.finite(survival$hazard_ratio_per_1sd)) ||
      any(!is.finite(survival$ci_lower_95)) ||
      any(!is.finite(survival$ci_upper_95)) ||
      any(survival$ci_lower_95 <= 0 | survival$ci_upper_95 <= 0),
    "生存表结构、HR/CI 或预锁定连续分数边界失败。"
  )
  survival_primary <- survival[
    model_variant == "ridge314_primary" & model_type == "clinical_adjusted"
  ]
  expected_survival <- data.table(
    cohort_scope = c(
      "GSE53622", "GSE53624", "GSE53622_GSE53624_stratified"
    ),
    n_patients = c(60L, 119L, 179L), n_events = c(33L, 73L, 106L)
  )
  survival_check <- merge(
    unique(survival_primary[, .(cohort_scope, n_patients, n_events)]),
    expected_survival, by = "cohort_scope", suffixes = c("", "_expected")
  )
  fail_if(
    nrow(survival_primary) != 6L || nrow(survival_check) != 3L ||
      any(survival_check$n_patients != survival_check$n_patients_expected) ||
      any(survival_check$n_events != survival_check$n_events_expected),
    "主 ridge 调整生存行未保持 60/33、119/73、179/106。"
  )

  increment <- copy(x$ecms_increment)
  increment[, independent_omics_validation := as_logical_strict(
    independent_omics_validation, "ECMS independent omics validation"
  )]
  corrected_columns <- c(
    "optimism_corrected_delta_cindex",
    "optimism_corrected_delta_cindex_ci_lower_95",
    "optimism_corrected_delta_cindex_ci_upper_95"
  )
  fail_if(
    nrow(increment) != 4L || !setequal(increment$factor, target_factors) ||
      !setequal(increment$model_variant, c(
        "ridge314_primary", "mofa51_fixed_weight_sensitivity"
      )) ||
      any(increment$n_patients != 179L) || any(increment$n_events != 106L) ||
      any(increment$bootstrap_replicates_requested != 300L) ||
      any(increment$bootstrap_optimism_replicates_valid < 240L) ||
      any(!is.finite(as.matrix(increment[, ..corrected_columns]))) ||
      any(increment$optimism_corrected_delta_cindex_ci_lower_95 >
            increment$optimism_corrected_delta_cindex) ||
      any(increment$optimism_corrected_delta_cindex_ci_upper_95 <
            increment$optimism_corrected_delta_cindex) ||
      any(increment$independent_omics_validation),
    paste(
      "ECMS increment 必须有 4 条完整行、300 次分层 bootstrap",
      "与真正 optimism-corrected delta C/CI。"
    )
  )

  response <- copy(x$response_associations)
  response[, treatment_predictor_claim := as_logical_strict(
    treatment_predictor_claim, "response treatment predictor claim"
  )]
  fail_if(
    nrow(response) != 4L || !setequal(response$factor, target_factors) ||
      any(response$n_total != 28L) ||
      any(response$n_pathological_complete_response != 11L) ||
      any(response$n_not_pathological_complete_response != 17L) ||
      any(response$treatment_predictor_claim) ||
      any(!is.finite(response$firth_odds_ratio_per_1sd)) ||
      any(!is.finite(response$firth_or_ci_lower_95)) ||
      any(!is.finite(response$firth_or_ci_upper_95)),
    "GSE45670 探索性反应表的 n、OR/CI 或非预测器边界失败。"
  )

  decision <- copy(x$validation_decision)
  decision[, `:=`(
    hard_gate = as_logical_strict(hard_gate, "decision hard gate"),
    counts_toward_external_support = as_logical_strict(
      counts_toward_external_support, "decision external support"
    )
  )]
  required_decisions <- c(
    "PROXY_CV_FACTOR1", "PROXY_CV_FACTOR3", "PROXY_MODEL_FROZEN",
    "MAP_GSE53622", "MAP_GSE53624", "MAP_GSE45670",
    "SURVIVAL_FACTOR1", "SURVIVAL_FACTOR3",
    "ECMS_INCREMENT_FACTOR1", "ECMS_INCREMENT_FACTOR3",
    "RESPONSE_FACTOR1", "RESPONSE_FACTOR3",
    "OVERALL_FIRST_STAGE_EXTERNAL_VALIDATION"
  )
  fail_if(
    nrow(decision) != 22L || anyDuplicated(decision$decision_id) ||
      !all(required_decisions %in% decision$decision_id) ||
      any(decision$hard_gate & !decision$status %chin% c(
        "PASS", "GO", "CONDITIONAL_GO"
      )),
    "外部验证决策表缺少 Figure 6 所需门禁或存在硬门禁失败。"
  )
  overall <- decision[
    decision_id == "OVERALL_FIRST_STAGE_EXTERNAL_VALIDATION"
  ]
  fail_if(
    nrow(overall) != 1L ||
      overall$status != "GO_FIRST_STAGE_EXTERNAL_RNA_PROXY_COMPLETED" ||
      !grepl("9.*(mut|CNV|\u4e8b\u4ef6|event)", overall$boundary, ignore.case = TRUE),
    "总体决策未明确保留 9 条精确事件—状态边尚未外部复现的边界。"
  )

  list(
    internal_cv = cv,
    oof_predictions = oof,
    patient_scores = patient,
    survival_associations = survival,
    ecms_increment = increment,
    response_associations = response,
    validation_decision = decision
  )
}

inputs <- validate_semantics(inputs)

source_table_rows <- data.table(
  input_name = names(source_table_paths),
  relative_path = file.path("results", basename(source_table_paths)),
  sha256 = vapply(source_table_paths, sha256_file, character(1))
)
setorder(source_table_rows, input_name)
source_table_sha256_bundle <- paste0(
  source_table_rows$relative_path, "=", source_table_rows$sha256,
  collapse = ";"
)
source_bundle_signature <- digest(
  list(
    source_manifest_sha256 = source_manifest_sha256,
    source_tables = source_table_rows,
    figure_contract_sha256 = figure_contract_sha256,
    script_sha256 = executed_script_sha256
  ),
  algo = "sha256", serialize = TRUE
)

font_family <- "Arial"
base_size <- 6.2
minimum_text_pt <- 5.5
palette_contract <- c(
  neutral_dark = "#2B2B2B",
  neutral_mid = "#7C8288",
  neutral_light = "#D9DDE1",
  neutral_pale = "#F2F3F4",
  factor1 = "#4F81BD",
  factor1_light = "#AFC8E6",
  factor3 = "#6CAFA7",
  factor3_light = "#B8D7D2",
  conditional = "#D79A55",
  null = "#C76B67",
  boundary = "#B39A6A"
)
factor_palette <- c(
  Factor1 = palette_contract[["factor1"]],
  Factor3 = palette_contract[["factor3"]]
)

theme_nature_contract <- function(base_size = 6.2) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      axis.line = element_line(linewidth = 0.30, colour = "#303030"),
      axis.ticks = element_line(linewidth = 0.30, colour = "#303030"),
      axis.title = element_text(size = base_size, colour = "#202020"),
      axis.text = element_text(size = base_size - 0.6, colour = "#303030"),
      legend.title = element_text(size = base_size - 0.2, face = "bold"),
      legend.text = element_text(size = base_size - 0.7),
      strip.background = element_rect(fill = "#F2F3F4", colour = NA),
      strip.text = element_text(size = base_size - 0.2, face = "bold"),
      plot.title = element_text(size = base_size + 0.5, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.4, colour = "#555555"),
      plot.caption = element_text(size = 5.5, colour = "#666666", hjust = 0),
      plot.tag = element_text(size = 8, face = "bold"),
      panel.grid = element_blank(),
      legend.key.height = grid::unit(3.0, "mm"),
      plot.margin = margin(2.0, 2.0, 2.0, 2.0, unit = "mm")
    )
}
theme_set(theme_nature_contract())

build_panel_a <- function(x) {
  repeats <- copy(x$internal_cv)[
    record_type == "outer_repeat" & model_variant == "ridge314_primary"
  ]
  summary_rows <- copy(x$internal_cv)[
    record_type == "outer_repeat_summary" &
      model_variant == "ridge314_primary"
  ]
  repeats[, factor := factor(factor, levels = c("Factor1", "Factor3"))]
  summary_rows[, factor := factor(factor, levels = c("Factor1", "Factor3"))]
  y_min <- min(0.70, min(repeats$spearman_rho) - 0.02)
  y_max <- max(repeats$spearman_rho) + 0.035
  p_cv <- ggplot(repeats, aes(x = factor, y = spearman_rho, fill = factor)) +
    geom_hline(
      yintercept = 0.70, linetype = "dashed", linewidth = 0.30,
      colour = palette_contract[["neutral_mid"]]
    ) +
    geom_boxplot(
      width = 0.48, outlier.shape = NA, linewidth = 0.35, alpha = 0.40
    ) +
    geom_point(
      position = position_jitter(width = 0.10, height = 0, seed = 20260716),
      size = 0.90, shape = 21, stroke = 0.22, colour = "white", alpha = 0.88
    ) +
    geom_point(
      data = summary_rows, aes(y = spearman_rho),
      shape = 23, size = 2.15, fill = "white", colour = "#202020",
      stroke = 0.45
    ) +
    geom_text(
      data = summary_rows,
      aes(
        y = spearman_rho,
        label = sprintf("median %.3f\n[%.3f–%.3f]", spearman_rho,
                        spearman_min, spearman_max)
      ),
      family = font_family, size = 1.78, nudge_x = 0.31, hjust = 0
    ) +
    scale_fill_manual(values = factor_palette, guide = "none") +
    coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
    labs(
      tag = "a", title = "Nested-CV recoverability in TCGA",
      subtitle = "20 repeated outer 5-fold splits; diamond = frozen summary",
      x = NULL, y = "Outer-CV Spearman rho"
    ) +
    theme_nature_contract() +
    theme(plot.margin = margin(2.5, 6.5, 0.8, 2.0, unit = "mm"))

  oof <- copy(x$oof_predictions)
  consistency <- oof[, .(
    observed_n = uniqueN(observed_factor_score)
  ), by = .(factor, patient_id)]
  fail_if(any(consistency$observed_n != 1L),
          "OOF 同患者观察因子分数在 repeat 间不唯一。")
  oof_patient <- oof[, .(
    observed_factor_score = unique(observed_factor_score),
    median_oof_prediction = median(oof_predicted_factor_score),
    repeats = .N
  ), by = .(factor, patient_id)]
  fail_if(any(oof_patient$repeats != 20L) || nrow(oof_patient) != 156L,
          "OOF 患者中位预测不是 2×78 行。")
  oof_annotation <- oof_patient[, .(
    label = sprintf(
      "patient-level rho = %.3f\nn = %d; 20 OOF values/patient",
      cor(observed_factor_score, median_oof_prediction, method = "spearman"),
      .N
    )
  ), by = factor]
  oof_annotation[, `:=`(x = -Inf, y = Inf)]
  p_oof <- ggplot(
    oof_patient,
    aes(x = observed_factor_score, y = median_oof_prediction, colour = factor)
  ) +
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.28,
      colour = palette_contract[["neutral_mid"]]
    ) +
    geom_point(size = 0.92, alpha = 0.72) +
    geom_text(
      data = oof_annotation,
      aes(x = x, y = y, label = label), inherit.aes = FALSE,
      family = font_family, size = 1.70, hjust = -0.05, vjust = 1.10,
      colour = palette_contract[["neutral_dark"]]
    ) +
    facet_wrap(~factor, scales = "free", nrow = 1) +
    scale_colour_manual(values = factor_palette, guide = "none") +
    labs(
      title = "Median out-of-fold prediction per patient",
      subtitle = "Only held-out predictions; no training-set fitted values",
      x = "Observed MOFA factor score", y = "Median OOF proxy score"
    ) +
    theme_nature_contract() +
    theme(plot.margin = margin(0.8, 2.0, 2.0, 2.0, unit = "mm"))

  p_cv / p_oof + plot_layout(heights = c(0.43, 0.57))
}

build_panel_b <- function(x) {
  patient <- copy(x$patient_scores)
  relation <- rbindlist(list(
    patient[, .(
      dataset, patient_id, factor = "Factor1",
      ridge314_z = factor1_ridge314_primary_z,
      mofa51_z = factor1_mofa51_sensitivity_z
    )],
    patient[, .(
      dataset, patient_id, factor = "Factor3",
      ridge314_z = factor3_ridge314_primary_z,
      mofa51_z = factor3_mofa51_sensitivity_z
    )]
  ))
  relation[, factor := factor(factor, levels = c("Factor1", "Factor3"))]
  relation[, dataset := factor(
    dataset, levels = c("GSE53622", "GSE53624", "GSE45670")
  )]
  relation_annotation <- relation[, .(
    label = sprintf(
      "rho=%.2f; n=%d",
      cor(ridge314_z, mofa51_z, method = "spearman"), .N
    ),
    x = -Inf, y = Inf
  ), by = .(factor, dataset)]
  p_relation <- ggplot(
    relation,
    aes(x = ridge314_z, y = mofa51_z, colour = factor)
  ) +
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.24,
      colour = palette_contract[["neutral_mid"]]
    ) +
    geom_point(size = 0.60, alpha = 0.60) +
    geom_text(
      data = relation_annotation,
      aes(x = x, y = y, label = label), inherit.aes = FALSE,
      family = font_family, size = 1.45, hjust = -0.05, vjust = 1.08,
      colour = palette_contract[["neutral_dark"]]
    ) +
    facet_grid(factor ~ dataset, scales = "free") +
    scale_colour_manual(values = factor_palette, guide = "none") +
    labs(
      tag = "b", title = "Transported score sensitivity",
      subtitle = "314-gene ridge proxy versus fixed 51-gene MOFA-weight score",
      x = "Ridge proxy (within-cohort z)",
      y = "51-gene sensitivity score (z)"
    ) +
    theme_nature_contract() +
    theme(
      axis.text = element_text(size = 4.9),
      strip.text = element_text(size = 5.2, face = "bold"),
      plot.margin = margin(2.5, 2.0, 0.5, 2.0, unit = "mm")
    )

  score_long <- rbindlist(list(
    patient[, .(
      dataset, patient_id, factor = "Factor1", variant = "Ridge 314",
      score_z = factor1_ridge314_primary_z
    )],
    patient[, .(
      dataset, patient_id, factor = "Factor1", variant = "MOFA 51",
      score_z = factor1_mofa51_sensitivity_z
    )],
    patient[, .(
      dataset, patient_id, factor = "Factor3", variant = "Ridge 314",
      score_z = factor3_ridge314_primary_z
    )],
    patient[, .(
      dataset, patient_id, factor = "Factor3", variant = "MOFA 51",
      score_z = factor3_mofa51_sensitivity_z
    )]
  ))
  score_long[, `:=`(
    factor = factor(factor, levels = c("Factor1", "Factor3")),
    dataset = factor(dataset, levels = c("GSE53622", "GSE53624", "GSE45670")),
    variant = factor(variant, levels = c("Ridge 314", "MOFA 51"))
  )]
  p_distribution <- ggplot(
    score_long, aes(x = dataset, y = score_z, fill = variant)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.22, colour = "#A0A4A8") +
    geom_boxplot(
      width = 0.62, outlier.shape = NA, linewidth = 0.27,
      position = position_dodge(width = 0.68), alpha = 0.78
    ) +
    facet_wrap(~factor, nrow = 1) +
    scale_fill_manual(values = c(
      "Ridge 314" = palette_contract[["factor1"]],
      "MOFA 51" = palette_contract[["neutral_light"]]
    )) +
    scale_x_discrete(labels = c(
      GSE53622 = "GSE\n53622", GSE53624 = "GSE\n53624",
      GSE45670 = "GSE\n45670"
    )) +
    labs(
      title = "Within-cohort score distributions",
      x = NULL, y = "Score (z)", fill = "Proxy definition"
    ) +
    theme_nature_contract() +
    theme(
      legend.position = "bottom", legend.direction = "horizontal",
      legend.margin = margin(0, 0, 0, 0),
      axis.text.x = element_text(size = 5.0, lineheight = 0.88),
      plot.margin = margin(0.5, 2.0, 2.0, 2.0, unit = "mm")
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))

  p_relation / p_distribution + plot_layout(heights = c(0.67, 0.33))
}

build_panel_c <- function(x) {
  forest <- copy(x$survival_associations)[
    model_variant == "ridge314_primary" & model_type == "clinical_adjusted"
  ]
  cohort_labels <- c(
    GSE53622 = "GSE53622",
    GSE53624 = "GSE53624",
    GSE53622_GSE53624_stratified = "Stratified family"
  )
  forest[, row_label := sprintf(
    "%s  %s  (n=%d; events=%d)",
    factor, cohort_labels[cohort_scope], n_patients, n_events
  )]
  desired <- unlist(lapply(c("Factor1", "Factor3"), function(factor_name) {
    sprintf(
      "%s  %s  (n=%d; events=%d)", factor_name,
      cohort_labels[c(
        "GSE53622", "GSE53624", "GSE53622_GSE53624_stratified"
      )], c(60L, 119L, 179L), c(33L, 73L, 106L)
    )
  }))
  forest[, row_label := factor(row_label, levels = rev(desired))]
  forest[, annotation := sprintf(
    "q=%s; PH P=%s%s",
    format_p(q_value), format_p(ph_score_p_value),
    fifelse(ph_score_violation_0_05, "; PH flag", "")
  )]
  label_x <- max(forest$ci_upper_95) * 1.20
  x_limits <- c(min(forest$ci_lower_95) * 0.88, label_x * 1.85)
  ggplot(
    forest,
    aes(y = row_label, x = hazard_ratio_per_1sd, colour = factor)
  ) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.30,
               colour = palette_contract[["neutral_mid"]]) +
    geom_segment(
      aes(x = ci_lower_95, xend = ci_upper_95, yend = row_label),
      linewidth = 0.55
    ) +
    geom_point(shape = 21, size = 2.0, fill = "white", stroke = 0.55) +
    geom_text(
      aes(x = label_x, label = annotation), hjust = 0,
      family = font_family, size = 1.55, colour = palette_contract[["neutral_dark"]]
    ) +
    scale_colour_manual(values = factor_palette, guide = "none") +
    scale_x_log10(
      limits = x_limits,
      breaks = c(0.5, 0.75, 1, 1.5, 2),
      labels = scales::label_number(accuracy = 0.01)
    ) +
    labs(
      tag = "c", title = "External overall-survival associations",
      subtitle = "Adjusted HR per 1 SD; two subseries from one study family",
      x = "Hazard ratio (95% CI, log scale)", y = NULL
    ) +
    theme_nature_contract() +
    theme(
      axis.text.y = element_text(size = 5.0),
      plot.margin = margin(2.5, 2.0, 2.0, 2.0, unit = "mm")
    )
}

build_panel_d <- function(x) {
  increment <- copy(x$ecms_increment)
  increment[, variant_label := fifelse(
    model_variant == "ridge314_primary", "Ridge 314", "MOFA 51"
  )]
  increment[, factor_short := fifelse(factor == "Factor1", "F1", "F3")]
  increment[, row_label := paste(factor_short, variant_label, sep = "  ·  ")]
  order_rows <- c(
    "F1  ·  Ridge 314", "F1  ·  MOFA 51",
    "F3  ·  Ridge 314", "F3  ·  MOFA 51"
  )
  increment[, row_label := factor(row_label, levels = rev(order_rows))]
  increment[, y_value := as.numeric(row_label)]
  increment[, annotation := sprintf(
    "LRT P=%s\ndelta AIC=%+.2f\nvalid=%d/300",
    format_p(lrt_p_value), delta_aic_full_minus_reduced,
    bootstrap_optimism_replicates_valid
  )]
  range_values <- range(c(
    increment$optimism_corrected_delta_cindex_ci_lower_95,
    increment$optimism_corrected_delta_cindex_ci_upper_95,
    increment$apparent_delta_cindex
  ))
  span <- diff(range_values)
  if (!is.finite(span) || span <= 0) span <- 0.02
  x_limits <- c(range_values[[1L]] - 0.08 * span,
                range_values[[2L]] + 0.08 * span)
  p_effect <- ggplot(increment, aes(y = y_value, colour = factor)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.30,
               colour = palette_contract[["neutral_mid"]]) +
    geom_segment(
      aes(
        x = optimism_corrected_delta_cindex_ci_lower_95,
        xend = optimism_corrected_delta_cindex_ci_upper_95,
        yend = y_value
      ), linewidth = 0.62
    ) +
    geom_segment(
      aes(
        x = optimism_corrected_delta_cindex_ci_lower_95,
        xend = optimism_corrected_delta_cindex_ci_lower_95,
        y = y_value - 0.11, yend = y_value + 0.11
      ), linewidth = 0.45
    ) +
    geom_segment(
      aes(
        x = optimism_corrected_delta_cindex_ci_upper_95,
        xend = optimism_corrected_delta_cindex_ci_upper_95,
        y = y_value - 0.11, yend = y_value + 0.11
      ), linewidth = 0.45
    ) +
    geom_segment(
      aes(
        x = apparent_delta_cindex,
        xend = optimism_corrected_delta_cindex,
        yend = y_value
      ), linewidth = 0.26, colour = palette_contract[["neutral_mid"]]
    ) +
    geom_point(
      aes(x = apparent_delta_cindex), shape = 23, size = 1.85,
      fill = "white", colour = "#303030", stroke = 0.45
    ) +
    geom_point(
      aes(x = optimism_corrected_delta_cindex, fill = factor),
      shape = 21, size = 2.05, colour = "white", stroke = 0.38
    ) +
    scale_colour_manual(values = factor_palette, guide = "none") +
    scale_fill_manual(values = factor_palette, guide = "none") +
    scale_y_continuous(
      breaks = increment$y_value,
      labels = as.character(increment$row_label),
      expand = expansion(mult = c(0.10, 0.10))
    ) +
    scale_x_continuous(
      limits = x_limits,
      labels = scales::label_number(accuracy = 0.01)
    ) +
    labs(
      tag = "d", title = "ECMS-adjusted incremental value",
      subtitle = "Filled point/CI: optimism-corrected delta C; diamond: apparent delta C",
      x = "Change in Harrell C-index", y = NULL
    ) +
    theme_nature_contract() +
    theme(
      axis.text.y = element_text(size = 5.0),
      plot.subtitle = element_text(size = 5.2),
      plot.margin = margin(2.5, 0.8, 2.0, 2.0, unit = "mm")
    )

  p_text <- ggplot(increment, aes(x = 0, y = y_value, label = annotation)) +
    geom_text(
      family = font_family, size = 1.46, hjust = 0, lineheight = 0.88,
      colour = palette_contract[["neutral_dark"]]
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(
      limits = c(0.5, 4.5), breaks = NULL, expand = c(0, 0)
    ) +
    labs(title = "Model-fit statistics") +
    theme_void(base_family = font_family, base_size = base_size) +
    theme(
      plot.title = element_text(size = 5.6, face = "bold", hjust = 0),
      plot.margin = margin(13.0, 1.0, 2.0, 0.8, unit = "mm")
    )

  (p_effect | p_text) + plot_layout(widths = c(1.40, 0.60))
}

build_panel_e <- function(x) {
  patient <- copy(x$patient_scores)[dataset == "GSE45670"]
  score <- rbindlist(list(
    patient[, .(
      patient_id, factor = "Factor1", response_binary,
      score_z = factor1_ridge314_primary_z
    )],
    patient[, .(
      patient_id, factor = "Factor3", response_binary,
      score_z = factor3_ridge314_primary_z
    )]
  ))
  score[, `:=`(
    factor = factor(factor, levels = c("Factor1", "Factor3")),
    response = factor(
      fifelse(response_binary == 1L, "pCR", "non-pCR"),
      levels = c("non-pCR", "pCR")
    )
  )]
  stats <- copy(x$response_associations)[
    model_variant == "ridge314_primary"
  ]
  stats[, label := sprintf(
    paste0(
      "r_rb = %+.2f; P_W = %s\n",
      "Firth OR = %.2f [%.2f, %.2f]\nP_Firth = %s"
    ),
    rank_biserial_pcr_higher_positive, format_p(wilcoxon_p_value),
    firth_odds_ratio_per_1sd, firth_or_ci_lower_95,
    firth_or_ci_upper_95, format_p(firth_profile_likelihood_p_value)
  )]
  stats[, factor := factor(factor, levels = c("Factor1", "Factor3"))]
  p_stats <- ggplot(stats, aes(x = factor, y = 0.5, label = label)) +
    geom_tile(
      width = 0.94, height = 0.76,
      fill = palette_contract[["neutral_pale"]], colour = "white",
      linewidth = 0.35
    ) +
    geom_text(
      family = font_family, size = 1.43, lineheight = 0.90,
      colour = palette_contract[["neutral_dark"]]
    ) +
    coord_cartesian(ylim = c(0, 1), clip = "off") +
    labs(
      tag = "e", title = "Exploratory pathological response in GSE45670",
      subtitle = "Frozen ridge scores; 11 pCR and 17 non-pCR; no classifier or cutpoint"
    ) +
    theme_void(base_family = font_family, base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 0.5, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.4, colour = "#555555"),
      plot.tag = element_text(size = 8, face = "bold"),
      plot.margin = margin(2.5, 2.0, 0.3, 2.0, unit = "mm")
    )

  p_scores <- ggplot(score, aes(x = response, y = score_z, fill = response)) +
    geom_hline(yintercept = 0, linewidth = 0.24, colour = "#A0A4A8") +
    geom_boxplot(
      width = 0.48, outlier.shape = NA, linewidth = 0.35, alpha = 0.70
    ) +
    geom_point(
      position = position_jitter(width = 0.11, height = 0, seed = 20260716),
      shape = 21, size = 0.78, stroke = 0.20, colour = "white", alpha = 0.84
    ) +
    facet_wrap(~factor, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = c(
      "non-pCR" = palette_contract[["neutral_light"]],
      pCR = palette_contract[["factor1"]]
    ), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.06))) +
    labs(
      x = NULL, y = "Proxy score (within-cohort z)"
    ) +
    theme_nature_contract() +
    theme(plot.margin = margin(0.3, 2.0, 2.0, 2.0, unit = "mm"))

  p_stats / p_scores + plot_layout(heights = c(0.34, 0.66))
}

build_panel_f <- function(x) {
  decision <- copy(x$validation_decision)
  get_decision <- function(id) {
    row <- decision[decision_id == id]
    fail_if(nrow(row) != 1L, paste("决策行不唯一：", id))
    row
  }
  proxy <- rbindlist(lapply(c("Factor1", "Factor3"), function(factor_name) {
    row <- get_decision(paste0("PROXY_CV_", toupper(factor_name)))
    data.table(
      evidence_domain = "Proxy recoverability", scope = factor_name,
      status = row$status, cell_label = "criterion met"
    )
  }))
  map_ids <- c("MAP_GSE53622", "MAP_GSE53624", "MAP_GSE45670")
  map_rows <- rbindlist(lapply(map_ids, get_decision))
  fail_if(any(map_rows$status != "PASS") || any(!map_rows$hard_gate),
          "患者映射决策不是 3/3 硬门禁 PASS。")
  mapping_counts <- sub(";.*$", "", map_rows$observed)
  mapping <- data.table(
    evidence_domain = "Patient mapping", scope = "Combined", status = "PASS",
    cell_label = paste0(
      "mapping complete\n", mapping_counts[[1L]], " · ",
      mapping_counts[[2L]], "\n", mapping_counts[[3L]]
    )
  )
  survival <- rbindlist(lapply(c("Factor1", "Factor3"), function(factor_name) {
    row <- get_decision(paste0("SURVIVAL_", toupper(factor_name)))
    data.table(
      evidence_domain = "Overall survival", scope = factor_name,
      status = row$status,
      cell_label = fcase(
        row$status == "CONDITIONAL_SAME_DIRECTION",
          "same direction;\nnot significant",
        row$status == "HETEROGENEOUS_OR_NULL",
          "null or\nheterogeneous",
        default = row$status
      )
    )
  }))
  increment <- rbindlist(lapply(c("Factor1", "Factor3"), function(factor_name) {
    row <- get_decision(paste0("ECMS_INCREMENT_", toupper(factor_name)))
    data.table(
      evidence_domain = "Increment beyond ECMS", scope = factor_name,
      status = row$status,
      cell_label = fifelse(
        row$status == "no_incremental_support",
        "no incremental\nsupport", row$status
      )
    )
  }))
  response <- rbindlist(lapply(c("Factor1", "Factor3"), function(factor_name) {
    row <- get_decision(paste0("RESPONSE_", toupper(factor_name)))
    data.table(
      evidence_domain = "Pathological response", scope = factor_name,
      status = row$status,
      cell_label = fcase(
        row$status == "exploratory_concordant_signal",
          "exploratory\nassociation",
        row$status == "exploratory_single_method_signal",
          "exploratory;\none method",
        row$status == "exploratory_null_or_inconclusive",
          "exploratory;\ninconclusive",
        default = row$status
      )
    )
  }))
  overall_row <- get_decision("OVERALL_FIRST_STAGE_EXTERNAL_VALIDATION")
  exact_boundary <- data.table(
    evidence_domain = "Exact event-edge replication", scope = "Combined",
    status = "NOT_TESTED",
    cell_label = "not assessed\n9 TCGA-only edges"
  )
  overall <- data.table(
    evidence_domain = "External RNA transportability", scope = "Combined",
    status = overall_row$status,
    cell_label = "cohort-level\napplication feasible"
  )
  panel <- rbindlist(list(
    proxy, mapping, survival, increment, response, exact_boundary, overall
  ))
  domain_levels <- c(
    "Proxy recoverability", "Patient mapping", "Overall survival",
    "Increment beyond ECMS", "Pathological response",
    "Exact event-edge replication", "External RNA transportability"
  )
  panel[, `:=`(
    evidence_domain = factor(evidence_domain, levels = rev(domain_levels)),
    scope = factor(scope, levels = c("Factor1", "Factor3", "Combined"))
  )]
  status_palette <- c(
    PASS = palette_contract[["factor1_light"]],
    GO = palette_contract[["factor1"]],
    CONDITIONAL_GO = palette_contract[["conditional"]],
    CONDITIONAL_SAME_DIRECTION = palette_contract[["conditional"]],
    HETEROGENEOUS_OR_NULL = "#E3C0BE",
    no_incremental_support = palette_contract[["neutral_light"]],
    exploratory_concordant_signal = palette_contract[["factor3_light"]],
    exploratory_single_method_signal = "#D5E3E1",
    exploratory_null_or_inconclusive = palette_contract[["neutral_pale"]],
    NOT_TESTED = "#E8DFC9",
    GO_FIRST_STAGE_EXTERNAL_RNA_PROXY_COMPLETED = palette_contract[["factor1"]]
  )
  missing_status <- setdiff(unique(as.character(panel$status)), names(status_palette))
  fail_if(length(missing_status) > 0L, paste(
    "panel f 出现未预锁定状态：", paste(missing_status, collapse = ", ")
  ))
  ggplot(panel, aes(x = scope, y = evidence_domain, fill = status)) +
    geom_tile(
      width = 0.92, height = 0.82, linewidth = 0.42, colour = "white"
    ) +
    geom_text(
      aes(label = cell_label), family = font_family, size = 1.52,
      colour = palette_contract[["neutral_dark"]], lineheight = 0.88
    ) +
    scale_fill_manual(values = status_palette, guide = "none") +
    labs(
      tag = "f", title = "Summary of transportability evidence",
      subtitle = "Reproducibility, clinical associations and event-level limits.",
      x = NULL, y = NULL
    ) +
    theme_nature_contract() +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.y = element_text(size = 5.1),
      axis.text.x = element_text(size = 5.4, face = "bold"),
      plot.margin = margin(2.5, 2.0, 2.0, 2.0, unit = "mm")
    )
}

build_figure6 <- function(x) {
  p_a <- build_panel_a(x)
  p_b <- build_panel_b(x)
  p_c <- build_panel_c(x)
  p_d <- build_panel_d(x)
  p_e <- build_panel_e(x)
  p_f <- build_panel_f(x)
  top <- p_a | p_b
  middle <- p_c | p_d
  bottom <- p_e | p_f
  top / middle / bottom +
    plot_layout(heights = c(1.35, 0.86, 0.96)) &
    theme(plot.background = element_rect(fill = "white", colour = NA))
}

validate_rendered_file <- function(path, format) {
  fail_if(!file_exists(path), paste("渲染文件缺失：", path))
  size <- as.numeric(file_info(path)$size)
  fail_if(!is.finite(size) || size < 10240,
          paste("渲染文件过小或为空：", path))
  raw_prefix <- readBin(path, what = "raw", n = min(size, 65536L))
  if (format == "pdf") {
    fail_if(
      rawToChar(raw_prefix[seq_len(min(5L, length(raw_prefix)))]) != "%PDF-",
      paste("PDF magic header 异常：", path)
    )
    pdf_info <- pdftools::pdf_info(path)
    pdf_fonts <- pdftools::pdf_fonts(path)
    fail_if(
      pdf_info$pages != 1L || !nrow(pdf_fonts) || any(!pdf_fonts$embedded),
      paste("PDF 页数或字体嵌入检查失败：", path)
    )
    expected_width_pt <- 183 / 25.4 * 72
    expected_height_pt <- 215 / 25.4 * 72
    pdf_size <- pdftools::pdf_pagesize(path)
    fail_if(
      nrow(pdf_size) != 1L ||
        abs(pdf_size$width[[1L]] - expected_width_pt) > 1.0 ||
        abs(pdf_size$height[[1L]] - expected_height_pt) > 1.0,
      paste("PDF 尺寸不是 183×215 mm：", path)
    )
  } else if (format == "png") {
    png_magic <- as.raw(c(0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A))
    fail_if(
      length(raw_prefix) < 8L || !identical(raw_prefix[1:8], png_magic),
      paste("PNG magic header 异常：", path)
    )
  } else if (format == "tiff") {
    tiff_le <- as.raw(c(0x49, 0x49, 0x2A, 0x00))
    tiff_be <- as.raw(c(0x4D, 0x4D, 0x00, 0x2A))
    fail_if(
      length(raw_prefix) < 4L ||
        (!identical(raw_prefix[1:4], tiff_le) &&
           !identical(raw_prefix[1:4], tiff_be)),
      paste("TIFF magic header 异常：", path)
    )
  } else if (format == "svg") {
    svg_text <- rawToChar(raw_prefix)
    fail_if(!grepl("<svg", svg_text, fixed = TRUE),
            paste("SVG 根元素缺失：", path))
    fail_if(!grepl("<text", svg_text, fixed = TRUE),
            paste("SVG 未保留可编辑 text 节点：", path))
    fail_if(grepl("<image[^>]+data:image", svg_text),
            paste("SVG 疑似整体栅格化：", path))
  }
  invisible(TRUE)
}

verify_current_figure_manifest <- function() {
  manifest_path <- file.path(project_root, figure_manifest_relative_path)
  fail_if(!file_exists(manifest_path), paste(
    "Figure 6 manifest 缺失：", manifest_path
  ))
  manifest <- read_tsv(manifest_path)
  require_columns(
    manifest,
    c(
      "figure_id", "relative_path", "format", "width_mm", "height_mm",
      "dpi", "file_size_bytes", "sha256", "frozen_figure_contract_sha256",
      "panel_d_primary_statistic", "backend", "font_family",
      "generation_script", "generation_script_sha256",
      "source_manifest_sha256", "source_table_sha256",
      "source_bundle_signature", "structural_status", "visual_qa_status",
      "qa_path", "qa_sha256"
    ),
    "Figure 6 artifact manifest"
  )
  fail_if(
    nrow(manifest) != 4L || uniqueN(manifest$relative_path) != 4L ||
      any(manifest$figure_id != figure_id) ||
      !setequal(manifest$relative_path, figure_relative_paths) ||
      !setequal(manifest$format, figure_formats) ||
      any(as.numeric(manifest$width_mm) != 183) ||
      any(as.numeric(manifest$height_mm) != 215) ||
      any(manifest$structural_status !=
            "verified_after_staged_export_and_hash_check") ||
      any(manifest$backend != "R-only") ||
      any(manifest$panel_d_primary_statistic != paste(
        "optimism_corrected_delta_cindex with",
        "optimism_corrected_delta_cindex_ci_lower_95/upper_95"
      )) ||
      any(manifest$frozen_figure_contract_sha256 != figure_contract_sha256) ||
      any(manifest$generation_script !=
            "scripts/26_visualize_external_validation.R") ||
      any(manifest$generation_script_sha256 != executed_script_sha256) ||
      any(manifest$source_manifest_sha256 != source_manifest_sha256) ||
      any(manifest$source_table_sha256 != source_table_sha256_bundle) ||
      any(manifest$source_bundle_signature != source_bundle_signature),
    "Figure 6 manifest 不是当前源表/脚本对应的完整 1×4 family。"
  )
  artifact_paths <- file.path(project_root, manifest$relative_path)
  fail_if(any(!file_exists(artifact_paths)), "Figure 6 四格式文件不完整。")
  actual_size <- as.numeric(file_info(artifact_paths)$size)
  actual_sha <- vapply(artifact_paths, sha256_file, character(1))
  fail_if(
    any(actual_size != as.numeric(manifest$file_size_bytes)) ||
      any(actual_sha != manifest$sha256),
    "Figure 6 当前文件大小/SHA256 与 manifest 不一致。"
  )
  for (index in seq_len(nrow(manifest))) {
    validate_rendered_file(artifact_paths[[index]], manifest$format[[index]])
  }
  manifest
}

atomic_replace <- function(source, destination) {
  fail_if(!file_exists(source), paste("原子替换源文件缺失：", source))
  temp_path <- tempfile(
    pattern = paste0(".", basename(destination), ".publish-"),
    tmpdir = dirname(destination)
  )
  on.exit(if (file_exists(temp_path)) file_delete(temp_path), add = TRUE)
  copied <- file.copy(source, temp_path, overwrite = FALSE, copy.mode = TRUE)
  fail_if(
    !copied || !file_exists(temp_path) ||
      as.numeric(file_info(temp_path)$size) != as.numeric(file_info(source)$size) ||
      sha256_file(temp_path) != sha256_file(source),
    paste("发布临时复制/SHA256 失败：", destination)
  )
  moved <- file.rename(temp_path, destination)
  fail_if(
    !moved || !file_exists(destination) ||
      sha256_file(destination) != sha256_file(source),
    paste("原子替换失败：", destination)
  )
  invisible(destination)
}

finalize_visual_qa_manifest <- function(qa_input) {
  fail_if(is.null(qa_input) || !nzchar(trimws(qa_input)),
          "--finalize-visual-qa 未提供 QA Markdown 路径。")
  qa_candidate <- path_abs(qa_input, start = project_root)
  fail_if(!file_exists(qa_candidate), paste("视觉 QA 文件缺失：", qa_candidate))
  qa_abs <- path_real(qa_candidate)
  checks_abs <- path_real(work_checks_dir)
  fail_if(
    !startsWith(as.character(qa_abs), paste0(as.character(checks_abs), "/")) ||
      tolower(path_ext(qa_abs)) != "md",
    "Figure 6 视觉 QA 必须是 _work/checks/ 中的 Markdown。"
  )
  qa_lines <- readLines(qa_abs, warn = FALSE, encoding = "UTF-8")
  normalized <- tolower(trimws(gsub("：", ":", qa_lines, fixed = TRUE)))
  fail_if(
    any(grepl(
      "figure\\s*6.*(fail|failed|\u5931\u8d25|\u4e0d\u901a\u8fc7)", normalized,
      perl = TRUE
    )),
    "Figure 6 视觉 QA 仍含 FAIL/不通过，禁止 finalize。"
  )
  pass_pattern <- paste0(
    "^[-*]?[[:space:]]*figure[[:space:]]*6[[:space:]]*:",
    "[[:space:]]*(pass|passed|\u901a\u8fc7)[[:space:]]*$"
  )
  pass_lines <- normalized[grepl(pass_pattern, normalized, perl = TRUE)]
  fail_if(length(pass_lines) != 1L,
          "视觉 QA 必须各含且只含一条 `Figure6: PASS`。")

  lock_dir <- file.path(work_intermediate_dir, ".escc_figure6_publish.lock")
  lock_acquired <- dir.create(lock_dir, recursive = FALSE, showWarnings = FALSE)
  fail_if(!lock_acquired, paste("Figure 6 发布锁已存在：", lock_dir))
  on.exit(if (dir_exists(lock_dir)) try(dir_delete(lock_dir), silent = TRUE), add = TRUE)

  manifest <- verify_current_figure_manifest()
  qa_relative <- as.character(path_rel(qa_abs, start = project_root))
  qa_sha <- sha256_file(qa_abs)
  manifest[, `:=`(
    qa_path = as.character(qa_path), qa_sha256 = as.character(qa_sha256)
  )]
  manifest[, `:=`(
    visual_qa_status = "passed_reopened_review",
    qa_path = qa_relative, qa_sha256 = qa_sha
  )]
  manifest_path <- file.path(project_root, figure_manifest_relative_path)
  stage_path <- tempfile(
    pattern = paste0(".", basename(manifest_path), ".visualqa-"),
    tmpdir = results_dir
  )
  on.exit(if (file_exists(stage_path)) try(file_delete(stage_path), silent = TRUE),
          add = TRUE)
  fwrite(
    manifest, stage_path, sep = "\t", quote = FALSE, na = "",
    logical01 = FALSE
  )
  staged <- read_tsv(stage_path)
  fail_if(
    nrow(staged) != 4L ||
      any(staged$visual_qa_status != "passed_reopened_review") ||
      any(staged$qa_path != qa_relative) || any(staged$qa_sha256 != qa_sha),
    "Figure 6 视觉 QA manifest 临时回读失败。"
  )
  atomic_replace(stage_path, manifest_path)
  if (file_exists(stage_path)) file_delete(stage_path)
  finalized <- verify_current_figure_manifest()
  fail_if(
    any(finalized$visual_qa_status != "passed_reopened_review") ||
      any(finalized$qa_path != qa_relative) ||
      any(finalized$qa_sha256 != qa_sha),
    "Figure 6 视觉 QA 原子替换后验证失败。"
  )
  dir_delete(lock_dir)
  cat("FINALIZE_VISUAL_QA_OK\n")
  cat("qa_path\t", qa_relative, "\n", sep = "")
  cat("qa_sha256\t", qa_sha, "\n", sep = "")
  invisible(finalized)
}

if (finalize_visual_qa) {
  finalize_visual_qa_manifest(finalize_visual_qa_input)
  quit(save = "no", status = 0L)
}

if (validate_only) {
  validated <- verify_current_figure_manifest()
  cat("VALIDATION_OK\n")
  cat("figure_artifacts\t", nrow(validated), "\n", sep = "")
  cat("pdf_pages\t1\n")
  cat("fonts_embedded\tTRUE\n")
  cat(
    "visual_qa_status\t",
    paste(sort(unique(validated$visual_qa_status)), collapse = ";"),
    "\n", sep = ""
  )
  quit(save = "no", status = 0L)
}

run_generation <- function() {
message("[1/4] 构建 Figure 6 六面板（R-only，不重拟合上游模型）")
figure <- build_figure6(inputs)
fail_if(!inherits(figure, c("ggplot", "patchwork")),
        "Figure 6 对象不是 ggplot/patchwork。")

stage_dir <- tempfile(
  pattern = ".escc_external_validation_figure6_",
  tmpdir = work_intermediate_dir
)
dir_create(stage_dir, recurse = TRUE)
stage_active <- TRUE
on.exit({
  if (stage_active && dir_exists(stage_dir)) try(dir_delete(stage_dir), silent = TRUE)
}, add = TRUE)

render_plot <- function(plot, path, format, dpi = 600L) {
  width_in <- 183 / 25.4
  height_in <- 215 / 25.4
  device_open <- FALSE
  on.exit(if (device_open) try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  if (format == "svg") {
    svglite::svglite(
      path, width = width_in, height = height_in, bg = "white",
      system_fonts = list(sans = font_family)
    )
  } else if (format == "pdf") {
    svg_source <- sub("[.]pdf$", ".svg", path)
    fail_if(!file_exists(svg_source), paste(
      "PDF 转换缺少已验证 SVG：", svg_source
    ))
    rsvg::rsvg_pdf(svg_source, path)
    return(invisible(path))
  } else if (format == "tiff") {
    ragg::agg_tiff(
      path, width = width_in, height = height_in, units = "in", res = dpi,
      background = "white", compression = "lzw", scaling = 1
    )
  } else if (format == "png") {
    ragg::agg_png(
      path, width = width_in, height = height_in, units = "in", res = dpi,
      background = "white", scaling = 1
    )
  } else {
    stop("未支持图件格式：", format, call. = FALSE)
  }
  device_open <- TRUE
  print(plot)
  grDevices::dev.off()
  device_open <- FALSE
  invisible(path)
}

message("[2/4] 导出 SVG/PDF/600 dpi TIFF/PNG 并执行结构 QA")
stage_records <- vector("list", length(figure_formats))
for (index in seq_along(figure_formats)) {
  format <- figure_formats[[index]]
  filename <- paste0(figure_stem, ".", format)
  stage_path <- file.path(stage_dir, filename)
  render_plot(figure, stage_path, format, dpi = 600L)
  validate_rendered_file(stage_path, format)
  stage_records[[index]] <- data.table(
    figure_id = figure_id,
    artifact = filename,
    relative_path = file.path("figures", filename),
    format = format,
    width_mm = 183,
    height_mm = 215,
    dpi = if (format %chin% c("tiff", "png")) 600L else NA_integer_,
    file_size_bytes = as.numeric(file_info(stage_path)$size),
    sha256 = sha256_file(stage_path),
    claim = paste(
      "Frozen Factor1/Factor3 RNA proxies are computable in external ESCC",
      "expression cohorts; survival, ECMS increment and pCR are secondary or",
      "exploratory calibrations and do not replicate genomic event-state edges."
    ),
    archetype = "asymmetric mixed quantitative figure",
    panel_source_contract = paste(
      "a=internal_cv+OOF; b=patient_scores; c=survival_associations;",
      "d=ecms_increment true optimism correction;",
      "e=response_associations+patient_scores; f=validation_decision"
    ),
    statistical_unit = paste(
      "outer split/TCGA patient; external GSM/patient; survival patient;",
      "bootstrap replicate; decision category"
    ),
    sample_structure = paste(
      "TCGA 78 for proxy CV; GSE53622 60/33 events;",
      "GSE53624 119/73 events; stratified family 179/106 events;",
      "GSE45670 28 with 11 pCR and 17 non-pCR"
    ),
    statistical_encoding = paste(
      "formal source tables only; no manual p/q/stars; effect sizes and CIs",
      "shown for survival, ECMS increment and response; categorical evidence summary for f"
    ),
    panel_d_primary_statistic = paste(
      "optimism_corrected_delta_cindex with",
      "optimism_corrected_delta_cindex_ci_lower_95/upper_95"
    ),
    reviewer_risk = paste(
      "shared RNA representation with ECMS; one GSE53625 study family;",
      "small exploratory pCR cohort; no exact genomic event replication"
    ),
    frozen_figure_contract_sha256 = figure_contract_sha256,
    backend = "R-only",
    font_family = font_family,
    minimum_text_pt = minimum_text_pt,
    generation_script = "scripts/26_visualize_external_validation.R",
    generation_script_sha256 = executed_script_sha256,
    source_manifest_sha256 = source_manifest_sha256,
    source_table_sha256 = source_table_sha256_bundle,
    source_bundle_signature = source_bundle_signature,
    generated_date = as.character(Sys.Date()),
    structural_status = "verified_after_staged_export_and_hash_check",
    visual_qa_status = "pending_reopened_review",
    qa_path = NA_character_,
    qa_sha256 = NA_character_
  )
}
artifact_manifest <- rbindlist(stage_records)
fail_if(
  nrow(artifact_manifest) != 4L ||
    !setequal(artifact_manifest$relative_path, figure_relative_paths),
  "Figure 6 manifest 必须恰好登记 1×4 文件族。"
)
manifest_stage_path <- file.path(
  stage_dir, basename(figure_manifest_relative_path)
)
fwrite(
  artifact_manifest, manifest_stage_path, sep = "\t", quote = FALSE,
  na = "", logical01 = FALSE
)
manifest_reread <- read_tsv(manifest_stage_path)
fail_if(
  nrow(manifest_reread) != 4L ||
    !identical(names(manifest_reread), names(artifact_manifest)) ||
    any(manifest_reread$sha256 != artifact_manifest$sha256),
  "Figure 6 manifest 临时写入后回读失败。"
)

message("[3/4] 发布四格式 Figure 6，manifest 最后发布")
dir_create(figures_dir, recurse = TRUE)
dir_create(results_dir, recurse = TRUE)
lock_dir <- file.path(work_intermediate_dir, ".escc_figure6_publish.lock")
lock_acquired <- dir.create(lock_dir, recursive = FALSE, showWarnings = FALSE)
fail_if(!lock_acquired, paste("Figure 6 发布锁已存在：", lock_dir))
on.exit(if (dir_exists(lock_dir)) try(dir_delete(lock_dir), silent = TRUE), add = TRUE)

destinations <- c(
  file.path(project_root, artifact_manifest$relative_path),
  file.path(project_root, figure_manifest_relative_path)
)
sources <- c(
  file.path(stage_dir, artifact_manifest$artifact), manifest_stage_path
)
backup_dir <- file.path(stage_dir, "backup")
dir_create(backup_dir, recurse = TRUE)
original_exists <- file_exists(destinations)
backup_paths <- file.path(backup_dir, basename(destinations))
for (index in which(original_exists)) {
  copied <- file.copy(
    destinations[[index]], backup_paths[[index]], overwrite = FALSE,
    copy.mode = TRUE
  )
  fail_if(
    !copied || sha256_file(destinations[[index]]) != sha256_file(backup_paths[[index]]),
    paste("Figure 6 发布前备份失败：", destinations[[index]])
  )
}

published <- rep(FALSE, length(destinations))
publish_error <- NULL
publish_ok <- tryCatch({
  for (index in seq_along(destinations)) {
    atomic_replace(sources[[index]], destinations[[index]])
    published[[index]] <- TRUE
  }
  TRUE
}, error = function(condition) {
  publish_error <<- condition
  FALSE
}, interrupt = function(condition) {
  publish_error <<- condition
  FALSE
})

if (!publish_ok) {
  rollback_errors <- character()
  for (index in rev(seq_along(destinations))) {
    tryCatch({
      if (original_exists[[index]]) {
        atomic_replace(backup_paths[[index]], destinations[[index]])
      } else if (published[[index]] && file_exists(destinations[[index]])) {
        file_delete(destinations[[index]])
      }
    }, error = function(condition) {
      rollback_errors <<- c(
        rollback_errors,
        paste0(destinations[[index]], " => ", conditionMessage(condition))
      )
    })
  }
  if (dir_exists(lock_dir)) dir_delete(lock_dir)
  fail_if(length(rollback_errors) > 0L, paste(
    "Figure 6 发布失败且回滚不完整：",
    conditionMessage(publish_error), "; ", paste(rollback_errors, collapse = " | ")
  ))
  stop(
    paste0(
      "Figure 6 发布失败，已回滚到原始大小/SHA256：",
      conditionMessage(publish_error)
    ),
    call. = FALSE
  )
}

fail_if(!all(published), "Figure 6 bundle 未全部发布。")
if (dir_exists(lock_dir)) dir_delete(lock_dir)
validated_manifest <- verify_current_figure_manifest()
fail_if(any(validated_manifest$visual_qa_status != "pending_reopened_review"),
        "新发布 Figure 6 不得在人工重开前冒充视觉 QA 通过。")

message("[4/4] Figure 6 结构包发布完成；等待最终尺寸人工重开 QA")
if (dir_exists(stage_dir)) dir_delete(stage_dir)
stage_active <- FALSE
message(
  "完成：Figure 6 的 SVG/PDF/TIFF/PNG 和独立 manifest 已发布；",
  "visual_qa_status=pending_reopened_review。"
)
invisible(TRUE)
}

run_generation()
