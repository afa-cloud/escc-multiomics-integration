#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clue)
  library(data.table)
  library(digest)
  library(fs)
  library(matrixStats)
  library(MOFA2)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
})

options(stringsAsFactors = FALSE, scipen = 999)

# TCGA-ESCC driver-event -> MOFA state 表示泄漏审计。
#
# 本脚本重建 script 14 的冻结五视图输入，并通过多 seed、删除事件视图和
# 删除主叙事候选自身 Mutation/CNV 特征来审查 script 19 的内部关联。
# 所有关联仍来自同一批 TCGA 患者，泄漏控制不能升级为独立验证或因果证明。
# ECMS 只允许使用公开且锁定的投影权重；当前没有时只写 pending，不生成伪标签。

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("无法唯一定位当前脚本。")
script_path <- normalizePath(
  sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE
)
project_root <- dirname(dirname(script_path))
results_dir <- file.path(project_root, "results")
work_intermediate_dir <- file.path(project_root, "_work", "intermediate")
work_checks_dir <- file.path(project_root, "_work", "checks")

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(prefix, default = NULL) {
  hit <- args[startsWith(args, prefix)]
  if (!length(hit)) return(default)
  if (length(hit) != 1L) stop("参数重复：", prefix)
  sub(prefix, "", hit, fixed = TRUE)
}
fields_only <- "--fields-only" %in% args
validate_only <- "--validate-only" %in% args
force_retrain <- "--force-retrain" %in% args

seed_text <- arg_value("--seeds=", "20260711,20260712,20260713")
seeds <- suppressWarnings(as.integer(strsplit(seed_text, ",", fixed = TRUE)[[1L]]))
frozen_gate_seeds <- c(20260711L, 20260712L, 20260713L)
if (!identical(seeds, frozen_gate_seeds)) {
  stop("正式泄漏门禁已预锁定为 seeds=20260711,20260712,20260713，不允许运行时改写。")
}
maxiter <- suppressWarnings(as.integer(arg_value("--maxiter=", "1500")))
alignment_min_abs_correlation <- suppressWarnings(as.numeric(
  arg_value("--alignment-min-abs-correlation=", "0.35")
))
alignment_min_margin <- suppressWarnings(as.numeric(
  arg_value("--alignment-min-margin=", "0.05")
))
alignment_high_confidence_override <- suppressWarnings(as.numeric(
  arg_value("--alignment-high-confidence-override=", "0.70")
))
alignment_method <- arg_value("--alignment-method=", "spearman")

if (!is.finite(maxiter) || maxiter < 100L) stop("--maxiter 必须为不小于 100 的整数。")
if (!alignment_method %in% c("spearman", "pearson")) {
  stop("--alignment-method 仅支持 spearman 或 pearson。")
}
if (any(!is.finite(c(
  alignment_min_abs_correlation,
  alignment_min_margin,
  alignment_high_confidence_override
)))) stop("因子对齐阈值必须为有限数值。")
if (alignment_min_abs_correlation <= 0 || alignment_min_abs_correlation > 1 ||
    alignment_min_margin < 0 || alignment_min_margin > 1 ||
    alignment_high_confidence_override < alignment_min_abs_correlation ||
    alignment_high_confidence_override > 1) {
  stop("因子对齐阈值范围不合法。")
}
if (!isTRUE(all.equal(alignment_min_margin, 0.05, tolerance = 0))) {
  stop("正式泄漏门禁的 minimum margin 已预锁定为 0.05，不允许运行时改写。")
}

required_project_paths <- c(
  file.path(project_root, "PROJECT_INDEX.md"),
  results_dir,
  work_intermediate_dir,
  work_checks_dir
)
if (!all(file_exists(required_project_paths) | dir_exists(required_project_paths))) {
  stop("项目权威索引或规范目录缺失：", project_root)
}

fail_if <- function(condition, message) {
  if (!is.logical(condition) || length(condition) != 1L || is.na(condition)) {
    stop(
      "fail_if 收到 NA、非逻辑或非标量条件；按 fail-closed 中止：",
      message,
      call. = FALSE
    )
  }
  if (condition) stop(message, call. = FALSE)
}

require_columns <- function(object, columns, object_name) {
  missing <- setdiff(columns, names(object))
  fail_if(
    length(missing) > 0L,
    paste0(object_name, " 缺少字段：", paste(missing, collapse = ", "))
  )
}

read_tsv <- function(path) {
  fread(path, na.strings = c("", "NA"), showProgress = FALSE)
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
    paste0(label, " 含 NA 或不可识别逻辑值。")
  )
  output
}

collapse_unique <- function(value, empty = "") {
  value <- unique(as.character(value[!is.na(value) & nzchar(as.character(value))]))
  if (!length(value)) empty else paste(sort(value), collapse = ";")
}

median_finite <- function(value) {
  value <- value[is.finite(value)]
  if (!length(value)) NA_real_ else median(value)
}

min_finite <- function(value) {
  value <- value[is.finite(value)]
  if (!length(value)) NA_real_ else min(value)
}

row_zscore <- function(x) {
  means <- rowMeans(x, na.rm = TRUE)
  sds <- rowSds(x, na.rm = TRUE, useNames = FALSE)
  keep <- is.finite(sds) & sds > 0
  fail_if(!all(keep), "冻结特征进入 row_zscore 前出现非有限或零方差。")
  z <- sweep(x, 1L, means, "-")
  sweep(z, 1L, sds, "/")
}

impute_row_median <- function(x) {
  medians <- rowMedians(x, na.rm = TRUE, useNames = FALSE)
  fail_if(any(!is.finite(medians)), "冻结 HM450/CNV 特征存在无法行中位数填补的记录。")
  missing <- which(!is.finite(x), arr.ind = TRUE)
  if (nrow(missing)) x[missing] <- medians[missing[, 1L]]
  x
}

weighted_ploidy_proxy <- function(copy_number, segment_length) {
  keep <- is.finite(copy_number) & is.finite(segment_length) & segment_length > 0
  if (!any(keep)) return(NA_real_)
  weightedMedian(copy_number[keep], w = segment_length[keep], na.rm = TRUE)
}

binary_rank_biserial <- function(event, outcome) {
  keep <- !is.na(event) & is.finite(outcome)
  event <- as.logical(event[keep])
  outcome <- outcome[keep]
  n_event <- sum(event)
  n_reference <- sum(!event)
  if (n_event == 0L || n_reference == 0L) return(NA_real_)
  ranks <- rank(outcome, ties.method = "average")
  u <- sum(ranks[event]) - n_event * (n_event + 1) / 2
  2 * u / (n_event * n_reference) - 1
}

binary_test <- function(event, outcome) {
  keep <- !is.na(event) & is.finite(outcome)
  event <- as.logical(event[keep])
  outcome <- outcome[keep]
  n_event <- sum(event)
  n_reference <- sum(!event)
  p_value <- if (n_event >= 3L && n_reference >= 3L) {
    suppressWarnings(wilcox.test(
      outcome[event], outcome[!event], exact = FALSE
    )$p.value)
  } else {
    NA_real_
  }
  list(
    n_complete = length(outcome),
    n_event = n_event,
    n_reference = n_reference,
    effect = binary_rank_biserial(event, outcome),
    p_value = p_value
  )
}

continuous_test <- function(predictor, outcome) {
  keep <- is.finite(predictor) & is.finite(outcome)
  predictor <- predictor[keep]
  outcome <- outcome[keep]
  test <- if (length(predictor) >= 10L &&
              length(unique(predictor)) >= 3L &&
              length(unique(outcome)) >= 3L) {
    suppressWarnings(cor.test(
      predictor, outcome, method = "spearman", exact = FALSE
    ))
  } else {
    NULL
  }
  list(
    n_complete = length(predictor),
    n_event = NA_integer_,
    n_reference = NA_integer_,
    effect = if (is.null(test)) NA_real_ else unname(test$estimate),
    p_value = if (is.null(test)) NA_real_ else test$p.value
  )
}

query_mofapy2_version <- function() {
  fail_if(!requireNamespace("basilisk", quietly = TRUE) ||
            !requireNamespace("reticulate", quietly = TRUE),
          "无法通过 basilisk/reticulate 查询实际 mofapy2 版本。")
  mofa_environment <- get("mofa_env", envir = asNamespace("MOFA2"))
  process <- basilisk::basiliskStart(mofa_environment)
  on.exit(basilisk::basiliskStop(process), add = TRUE)
  version <- basilisk::basiliskRun(process, function() {
    module <- reticulate::import("mofapy2")
    as.character(module$version$`__version__`)
  })
  fail_if(length(version) != 1L || is.na(version) || !nzchar(version),
          "实际 mofapy2 版本查询失败。")
  version
}

extract_training_diagnostics <- function(fitted, configured_maxiter) {
  status <- as.character(fitted@status)
  elbo <- as.numeric(fitted@training_stats$elbo)
  elapsed_records <- as.numeric(fitted@training_stats$time)
  factor_records <- as.numeric(fitted@training_stats$number_factors)
  finite_elbo_indices <- which(is.finite(elbo))
  finite_elbo <- elbo[finite_elbo_indices]
  iteration_count <- max(0L, length(elapsed_records) - 1L)
  final_elbo <- if (length(finite_elbo)) tail(finite_elbo, 1L) else NA_real_
  last_relative_change <- if (length(finite_elbo) >= 2L) {
    previous <- finite_elbo[[length(finite_elbo) - 1L]]
    abs(final_elbo - previous) / max(1, abs(previous))
  } else {
    NA_real_
  }
  tail_elbo <- tail(finite_elbo, 6L)
  tail_non_decreasing_fraction <- if (length(tail_elbo) >= 2L) {
    mean(diff(tail_elbo) >= -1e-8)
  } else {
    NA_real_
  }
  finite_factor_records <- factor_records[is.finite(factor_records)]
  final_factor_count <- if (length(finite_factor_records)) {
    as.integer(tail(finite_factor_records, 1L))
  } else {
    NA_integer_
  }
  model_completed <- identical(status, "trained") &&
    iteration_count > 0L && length(finite_elbo) >= 2L && is.finite(final_elbo)
  model_converged <- model_completed &&
    is.finite(tail_non_decreasing_fraction) &&
    tail_non_decreasing_fraction >= 0.80 &&
    (
      iteration_count < configured_maxiter ||
        (is.finite(last_relative_change) && last_relative_change <= 1e-4)
    )
  list(
    mofa_status = status,
    model_completed = model_completed,
    model_converged = model_converged,
    training_iteration_count = as.integer(iteration_count),
    configured_maxiter = as.integer(configured_maxiter),
    elbo_finite_record_count = as.integer(length(finite_elbo)),
    final_elbo_record_index = if (length(finite_elbo_indices)) {
      as.integer(tail(finite_elbo_indices, 1L))
    } else {
      NA_integer_
    },
    final_elbo = final_elbo,
    last_relative_elbo_change = last_relative_change,
    tail_elbo_non_decreasing_fraction = tail_non_decreasing_fraction,
    final_factor_count_from_training_stats = final_factor_count,
    convergence_rule = paste(
      "status=trained; finite ELBO; >=80% non-decreasing among last six;",
      "early stop before maxiter or final relative ELBO change <=1e-4"
    )
  )
}

extract_event_blind_loadings <- function(
    fitted,
    frozen_feature_ids,
    common_views = c("RNA", "miRNA", "HM450")) {
  weights <- get_weights(
    fitted,
    views = common_views,
    factors = "all",
    as.data.frame = FALSE
  )
  fail_if(!identical(names(weights), common_views),
          "MOFA common non-event loading views 顺序异常。")
  standardized <- vector("list", length(common_views))
  names(standardized) <- common_views
  factor_names <- NULL
  for (view_name in common_views) {
    matrix <- as.matrix(weights[[view_name]])
    expected_ids <- frozen_feature_ids[[view_name]]
    fail_if(nrow(matrix) != length(expected_ids) || any(!is.finite(matrix)),
            paste("common loading 维度或有限性异常：", view_name))
    observed_ids <- rownames(matrix)
    fail_if(is.null(observed_ids) || anyNA(observed_ids) ||
              anyDuplicated(observed_ids) > 0L,
            paste("common loading 特征 ID 缺失或重复：", view_name))
    normalized_ids <- observed_ids
    view_suffix <- paste0("_", view_name)
    suffix_present <- endsWith(observed_ids, view_suffix)
    stripped_ids <- observed_ids
    stripped_ids[suffix_present] <- substr(
      observed_ids[suffix_present],
      1L,
      nchar(observed_ids[suffix_present]) - nchar(view_suffix)
    )
    # MOFA2 只会对跨视图重名 ID 追加 _<view>；仅当剥离后与同位
    # 冻结 ID 严格相等时才允许该标准化，不掩盖真实错序/错身份。
    suffix_is_mofa_namespace <- suffix_present & stripped_ids == expected_ids
    normalized_ids[suffix_is_mofa_namespace] <-
      stripped_ids[suffix_is_mofa_namespace]
    fail_if(!identical(normalized_ids, expected_ids),
            paste("common loading 与冻结 manifest 特征身份/顺序不一致：", view_name))
    if (is.null(factor_names)) {
      factor_names <- colnames(matrix)
    } else {
      fail_if(!identical(colnames(matrix), factor_names),
              "common loading 各视图因子顺序不一致。")
    }
    fail_if(is.null(factor_names) || anyDuplicated(factor_names) > 0L,
            "common loading 因子名缺失或重复。")
    # 删除/no-op 判定在 create_mofa 前使用 view-qualified 原始 ID；此处按冻结
    # 输入位置回填，避免 MOFA2 对跨视图重复 ID 自动加后缀导致错配。
    rownames(matrix) <- paste(view_name, expected_ids, sep = "::")
    centered <- sweep(matrix, 2L, colMeans(matrix), "-")
    norms <- sqrt(colSums(centered^2))
    fail_if(any(!is.finite(norms)) || any(norms <= 0),
            paste("common loading 无法按视图标准化：", view_name))
    standardized[[view_name]] <- sweep(centered, 2L, norms, "/")
  }
  combined <- do.call(rbind, standardized)
  fail_if(any(!is.finite(combined)), "拼接后的 event-blind common loadings 含非有限值。")
  combined
}

verify_manifest_inputs <- function(manifest_path, input_paths, full_sha = TRUE) {
  fail_if(!file_exists(manifest_path), paste("artifact manifest 缺失：", manifest_path))
  manifest <- read_tsv(manifest_path)
  require_columns(
    manifest,
    c("relative_path", "file_size_bytes", "sha256", "status"),
    basename(manifest_path)
  )
  fail_if(anyDuplicated(manifest$relative_path) > 0L,
          paste("manifest relative_path 重复：", basename(manifest_path)))
  expected_relative <- file.path("results", basename(input_paths))
  rows <- manifest[match(expected_relative, relative_path)]
  fail_if(
    nrow(rows) != length(input_paths) || anyNA(rows$relative_path),
    paste("manifest 未覆盖全部输入：", basename(manifest_path))
  )
  fail_if(any(rows$status != "verified"),
          paste("manifest 输入未标记 verified：", basename(manifest_path)))
  fail_if(any(!file_exists(input_paths)), "manifest 已登记输入文件缺失。")
  actual_size <- as.numeric(file_info(input_paths)$size)
  fail_if(any(actual_size != as.numeric(rows$file_size_bytes)),
          paste("输入大小与 manifest 不一致：", basename(manifest_path)))
  if (full_sha) {
    actual_sha <- vapply(
      input_paths,
      digest,
      FUN.VALUE = character(1),
      algo = "sha256",
      file = TRUE,
      serialize = FALSE
    )
    fail_if(any(actual_sha != rows$sha256),
            paste("输入 SHA256 与 manifest 不一致：", basename(manifest_path)))
  }
  rows[, .(relative_path, file_size_bytes, sha256, status)]
}

required_inputs <- list(
  mae = file.path(results_dir, "tcga_escc_multiassay_dr45.rds"),
  segments = file.path(results_dir, "tcga_escc_ascat2_segments_long.rds"),
  analysis_sets = file.path(results_dir, "tcga_escc_multiassay_analysis_sets.tsv"),
  driver_screen = file.path(results_dir, "tcga_escc_driver_candidate_screen.tsv"),
  feature_manifest = file.path(results_dir, "tcga_escc_heterogeneity_feature_manifest.tsv"),
  reference_factors = file.path(results_dir, "tcga_escc_mofa_factor_scores.tsv"),
  reference_model = file.path(results_dir, "tcga_escc_mofa_model.rds"),
  level_factor_qc = file.path(results_dir, "tcga_escc_mofa_level_factor_qc.tsv"),
  patient_events = file.path(results_dir, "tcga_escc_strong_driver_patient_events.tsv"),
  original_associations = file.path(results_dir, "tcga_escc_driver_state_associations.tsv"),
  original_network_edges = file.path(results_dir, "tcga_escc_driver_state_network_edges.tsv")
)
fail_if(
  any(!file_exists(unlist(required_inputs))),
  paste(
    "缺少泄漏审计输入：",
    paste(unlist(required_inputs)[!file_exists(unlist(required_inputs))], collapse = ", ")
  )
)

full_input_sha <- !fields_only
multiassay_manifest_rows <- verify_manifest_inputs(
  file.path(results_dir, "tcga_escc_multiassay_artifact_manifest.tsv"),
  unlist(required_inputs[c("mae", "segments", "analysis_sets")]),
  full_sha = full_input_sha
)
driver_manifest_rows <- verify_manifest_inputs(
  file.path(results_dir, "tcga_escc_driver_core_artifact_manifest.tsv"),
  required_inputs$driver_screen,
  full_sha = full_input_sha
)
heterogeneity_manifest_rows <- verify_manifest_inputs(
  file.path(results_dir, "tcga_escc_heterogeneity_artifact_manifest.tsv"),
  unlist(required_inputs[c(
    "feature_manifest", "reference_factors", "reference_model", "level_factor_qc"
  )]),
  full_sha = full_input_sha
)
state_manifest_rows <- verify_manifest_inputs(
  file.path(results_dir, "tcga_escc_driver_state_artifact_manifest.tsv"),
  unlist(required_inputs[c(
    "patient_events", "original_associations", "original_network_edges"
  )]),
  # network edge/association 是本审计的正式门禁源；即使 fields-only 也做 SHA256。
  full_sha = TRUE
)

# ECMS 是独立的下游外部坐标，不进入任何 MOFA 训练、因子匹配或
# driver-state 门禁。若 script23 已完整发布，则这里只校验其正式
# artifact 并把“可用但未参与门禁”的状态写入本审计；若尚未发布则保留 pending。
ecms_projection_inputs <- list(
  patient = file.path(results_dir, "tcga_escc_ecms_patient_probabilities.tsv"),
  calibration = file.path(results_dir, "tcga_escc_ecms_projection_calibration.tsv"),
  qa = file.path(results_dir, "tcga_escc_ecms_projection_qa.tsv"),
  factor_associations = file.path(
    results_dir, "tcga_escc_ecms_factor_associations.tsv"
  ),
  adjusted_factor_progeny = file.path(
    results_dir, "tcga_escc_ecms_adjusted_factor_progeny_associations.tsv"
  ),
  summary = file.path(results_dir, "tcga_escc_ecms_projection_summary.md")
)
ecms_projection_manifest_path <- file.path(
  results_dir, "tcga_escc_ecms_projection_artifact_manifest.tsv"
)
ecms_expected_relative_paths <- file.path(
  "results", basename(unlist(ecms_projection_inputs))
)
ecms_expected_family_paths <- c(
  unlist(ecms_projection_inputs), ecms_projection_manifest_path
)
ecms_projection_family_files <- dir_ls(
  results_dir,
  regexp = "tcga_escc_ecms_",
  type = "file",
  fail = FALSE
)
ecms_projection_family_files <- as.character(ecms_projection_family_files)
ecms_expected_family_paths <- as.character(ecms_expected_family_paths)
ecms_projection_any <- length(ecms_projection_family_files) > 0L
ecms_projection_available <- all(file_exists(ecms_expected_family_paths))
fail_if(
  ecms_projection_any && (
    !ecms_projection_available ||
      !setequal(ecms_projection_family_files, ecms_expected_family_paths)
  ),
  "ECMS 投影正式文件必须恰好为 script23 的 6 个 artifact 加 1 个 manifest。"
)
ecms_manifest_rows <- data.table()
ecms_projection_signature <- "pending_no_verified_ecms_projection"
ecms_extension_calibration_pass <- NA
ecms_primary_patient_count <- NA_integer_
ecms_anchor_patient_count <- NA_integer_
ecms_model_commit <- ""
if (ecms_projection_available) {
  ecms_manifest_rows <- verify_manifest_inputs(
    ecms_projection_manifest_path,
    unlist(ecms_projection_inputs),
    full_sha = TRUE
  )
  ecms_manifest <- read_tsv(ecms_projection_manifest_path)
  require_columns(
    ecms_manifest,
    c(
      "relative_path", "file_size_bytes", "sha256", "status",
      "generation_script", "execution_script_sha256", "model_commit",
      "extension_calibration_pass", "primary_patient_count"
    ),
    "tcga_escc_ecms_projection_artifact_manifest.tsv"
  )
  ecms_generator_path <- file.path(
    project_root, "scripts", "23_project_tcga_escc_ecms.R"
  )
  fail_if(!file_exists(ecms_generator_path), "ECMS 投影生成脚本缺失。")
  current_ecms_generator_sha <- digest(
    ecms_generator_path,
    algo = "sha256", file = TRUE, serialize = FALSE
  )
  fail_if(
    nrow(ecms_manifest) != length(ecms_expected_relative_paths) ||
      anyDuplicated(ecms_manifest$relative_path) > 0L ||
      !setequal(ecms_manifest$relative_path, ecms_expected_relative_paths) ||
      any(ecms_manifest$status != "verified") ||
      any(ecms_manifest$generation_script !=
            "scripts/23_project_tcga_escc_ecms.R") ||
      any(ecms_manifest$execution_script_sha256 != current_ecms_generator_sha) ||
      any(!grepl("^[0-9a-f]{64}$", ecms_manifest$sha256)) ||
      uniqueN(ecms_manifest$model_commit) != 1L ||
      ecms_manifest$model_commit[[1L]] !=
        "8a49db65f86f2c1e2a50d9bc4cb742ccc868e4c6" ||
      uniqueN(ecms_manifest$extension_calibration_pass) != 1L ||
      uniqueN(ecms_manifest$primary_patient_count) != 1L,
    "ECMS 投影 manifest 的完整路径集、生成脚本、SHA、commit 或状态失败。"
  )
  ecms_manifest_primary_numeric <- suppressWarnings(as.numeric(
    ecms_manifest$primary_patient_count
  ))
  fail_if(
    any(!is.finite(ecms_manifest_primary_numeric)) ||
      any(abs(ecms_manifest_primary_numeric -
                round(ecms_manifest_primary_numeric)) > 1e-12) ||
      uniqueN(ecms_manifest_primary_numeric) != 1L,
    "ECMS manifest primary_patient_count 必须是唯一整数。"
  )
  ecms_patient <- read_tsv(ecms_projection_inputs$patient)
  ecms_calibration <- read_tsv(ecms_projection_inputs$calibration)
  ecms_qa <- read_tsv(ecms_projection_inputs$qa)
  require_columns(
    ecms_patient,
    c(
      "patient_id", "in_official_78", "resolved_ecms_label",
      "official_anchor_label", "gdc_projection_label",
      "eligible_for_primary_association", "pseudo_label_generated",
      "single_sample_classifier_claim", "label_source",
      "resolved_margin_custom", "margin_is_official_threshold",
      "low_margin_rejects_label",
      paste0("official_prob_ECMS", 1:4),
      paste0("gdc_prob_ECMS", 1:4),
      paste0("resolved_prob_ECMS", 1:4)
    ),
    "tcga_escc_ecms_patient_probabilities.tsv"
  )
  require_columns(
    ecms_calibration,
    c("metric", "observed", "threshold", "direction", "gate_role", "pass"),
    "tcga_escc_ecms_projection_calibration.tsv"
  )
  require_columns(
    ecms_qa,
    c("check_id", "hard_gate", "pass"),
    "tcga_escc_ecms_projection_qa.tsv"
  )
  ecms_expected_patients_table <- read_tsv(required_inputs$analysis_sets)
  require_columns(
    ecms_expected_patients_table,
    c("analysis_set", "patient_id", "included"),
    "tcga_escc_multiassay_analysis_sets.tsv"
  )
  ecms_included_flag <- as_logical_strict(
    ecms_expected_patients_table$included,
    "analysis_sets included"
  )
  ecms_expected_patients <- ecms_expected_patients_table[
    analysis_set == "five_layer_core" & ecms_included_flag,
    patient_id
  ]
  fail_if(
    length(ecms_expected_patients) != 94L ||
      uniqueN(ecms_expected_patients) != 94L,
    "ECMS 校验无法冻结 94 位 five-layer core 患者。"
  )
  ecms_anchor_flag <- as_logical_strict(
    ecms_patient$in_official_78, "ECMS in_official_78"
  )
  ecms_eligible_flag <- as_logical_strict(
    ecms_patient$eligible_for_primary_association,
    "ECMS eligible_for_primary_association"
  )
  ecms_pseudo_flag <- as_logical_strict(
    ecms_patient$pseudo_label_generated,
    "ECMS pseudo_label_generated"
  )
  ecms_clinical_flag <- as_logical_strict(
    ecms_patient$single_sample_classifier_claim,
    "ECMS single_sample_classifier_claim"
  )
  ecms_margin_official_flag <- as_logical_strict(
    ecms_patient$margin_is_official_threshold,
    "ECMS margin_is_official_threshold"
  )
  ecms_margin_reject_flag <- as_logical_strict(
    ecms_patient$low_margin_rejects_label,
    "ECMS low_margin_rejects_label"
  )
  ecms_levels <- paste0("ECMS", 1:4)
  ecms_resolved_probability_columns <- paste0("resolved_prob_ECMS", 1:4)
  ecms_gdc_probability_columns <- paste0("gdc_prob_ECMS", 1:4)
  ecms_official_probability_columns <- paste0("official_prob_ECMS", 1:4)
  ecms_resolved_probabilities <- as.matrix(
    ecms_patient[, ..ecms_resolved_probability_columns]
  )
  ecms_gdc_probabilities <- as.matrix(
    ecms_patient[, ..ecms_gdc_probability_columns]
  )
  ecms_official_probabilities <- as.matrix(
    ecms_patient[, ..ecms_official_probability_columns]
  )
  storage.mode(ecms_resolved_probabilities) <- "double"
  storage.mode(ecms_gdc_probabilities) <- "double"
  storage.mode(ecms_official_probabilities) <- "double"
  ecms_resolved_label_index <- match(
    ecms_patient$resolved_ecms_label, ecms_levels
  )
  ecms_resolved_selected_probability <- ecms_resolved_probabilities[cbind(
    seq_len(nrow(ecms_resolved_probabilities)), ecms_resolved_label_index
  )]
  ecms_resolved_max_probability <- apply(
    ecms_resolved_probabilities, 1L, max
  )
  ecms_gdc_label_index <- match(
    ecms_patient$gdc_projection_label, ecms_levels
  )
  ecms_gdc_selected_probability <- ecms_gdc_probabilities[cbind(
    seq_len(nrow(ecms_gdc_probabilities)), ecms_gdc_label_index
  )]
  ecms_gdc_max_probability <- apply(ecms_gdc_probabilities, 1L, max)
  ecms_official_label_index <- match(
    ecms_patient$official_anchor_label, ecms_levels
  )
  ecms_official_selected_probability <- ecms_official_probabilities[cbind(
    seq_len(nrow(ecms_official_probabilities)), ecms_official_label_index
  )]
  ecms_official_max_probability <- rep(
    NA_real_, nrow(ecms_official_probabilities)
  )
  ecms_official_max_probability[ecms_anchor_flag] <- apply(
    ecms_official_probabilities[ecms_anchor_flag, , drop = FALSE],
    1L, max
  )
  ecms_resolved_sorted_probability <- t(apply(
    ecms_resolved_probabilities, 1L, sort, decreasing = TRUE
  ))
  ecms_resolved_margin_from_probability <-
    ecms_resolved_sorted_probability[, 1L] -
      ecms_resolved_sorted_probability[, 2L]
  ecms_allowed_label_sources <- c(
    "locked_model_prediction_on_repository_bundled_tcga78_anchor",
    "gdc_tpm_extension_after_overlap_calibration_pass",
    "gdc_tpm_extension_conditional_not_primary"
  )
  fail_if(
    nrow(ecms_patient) != 94L || uniqueN(ecms_patient$patient_id) != 94L ||
      !setequal(ecms_patient$patient_id, ecms_expected_patients) ||
      sum(ecms_anchor_flag) != 78L ||
      any(ecms_pseudo_flag) || any(ecms_clinical_flag) ||
      any(ecms_margin_official_flag) || any(ecms_margin_reject_flag) ||
      any(!ecms_patient$resolved_ecms_label %chin% ecms_levels) ||
      any(!ecms_patient$gdc_projection_label %chin% ecms_levels) ||
      any(ecms_anchor_flag &
            !ecms_patient$official_anchor_label %chin% ecms_levels) ||
      any(!ecms_patient$label_source %chin% ecms_allowed_label_sources) ||
      any(!is.finite(ecms_resolved_probabilities)) ||
      any(!is.finite(ecms_gdc_probabilities)) ||
      any(ecms_resolved_probabilities < 0 | ecms_resolved_probabilities > 1) ||
      any(ecms_gdc_probabilities < 0 | ecms_gdc_probabilities > 1) ||
      any(abs(rowSums(ecms_resolved_probabilities) - 1) > 1e-12) ||
      any(abs(rowSums(ecms_gdc_probabilities) - 1) > 1e-12) ||
      any(abs(
        ecms_resolved_selected_probability - ecms_resolved_max_probability
      ) > 1e-12) ||
      any(abs(ecms_gdc_selected_probability - ecms_gdc_max_probability) >
            1e-12) ||
      any(ecms_anchor_flag & abs(
        ecms_official_selected_probability - ecms_official_max_probability
      ) > 1e-12) ||
      any(!is.finite(ecms_patient$resolved_margin_custom)) ||
      any(ecms_patient$resolved_margin_custom < 0 |
            ecms_patient$resolved_margin_custom > 1) ||
      any(abs(
        ecms_patient$resolved_margin_custom -
          ecms_resolved_margin_from_probability
      ) > 1e-12) ||
      any(ecms_anchor_flag &
            ecms_patient$resolved_ecms_label != ecms_patient$official_anchor_label) ||
      any(!ecms_anchor_flag &
            ecms_patient$resolved_ecms_label != ecms_patient$gdc_projection_label) ||
      any(ecms_anchor_flag & !is.finite(ecms_official_probabilities)) ||
      any(ecms_anchor_flag & (
        ecms_official_probabilities < 0 | ecms_official_probabilities > 1
      )) ||
      any(abs(rowSums(
        ecms_official_probabilities[ecms_anchor_flag, , drop = FALSE]
      ) - 1) > 1e-12) ||
      any(abs(
        ecms_resolved_probabilities[ecms_anchor_flag, , drop = FALSE] -
          ecms_official_probabilities[ecms_anchor_flag, , drop = FALSE]
      ) > 1e-12) ||
      any(abs(
        ecms_resolved_probabilities[!ecms_anchor_flag, , drop = FALSE] -
          ecms_gdc_probabilities[!ecms_anchor_flag, , drop = FALSE]
      ) > 1e-12) ||
      any(!ecms_anchor_flag & !is.na(ecms_patient$official_anchor_label)) ||
      any(!ecms_anchor_flag & !is.na(ecms_official_probabilities)),
    "ECMS 患者集、标签来源、四类概率、anchor 或临床边界失败。"
  )
  ecms_manifest_extension <- as_logical_strict(
    ecms_manifest$extension_calibration_pass,
    "ECMS manifest extension_calibration_pass"
  )
  ecms_extension_calibration_pass <- unique(ecms_manifest_extension)
  ecms_gate_names <- c(
    "exact_label_agreement", "cohen_kappa", "adjusted_rand_index",
    "median_class_probability_spearman"
  )
  ecms_gate_thresholds <- c(0.75, 0.60, 0.60, 0.50)
  ecms_gate_rows <- ecms_calibration[
    metric %chin% ecms_gate_names & gate_role == "prelocked_extension_gate"
  ]
  ecms_overall <- ecms_calibration[
    metric == "overall_extension_calibration_pass" &
      gate_role == "extension_decision"
  ]
  fail_if(
    nrow(ecms_gate_rows) != 4L ||
      uniqueN(ecms_gate_rows$metric) != 4L ||
      !setequal(ecms_gate_rows$metric, ecms_gate_names) ||
      nrow(ecms_overall) != 1L,
    "ECMS calibration 缺少唯一的预锁定四门禁或 overall 行。"
  )
  ecms_gate_pass <- as_logical_strict(
    ecms_gate_rows$pass, "ECMS four prelocked gate pass"
  )
  ecms_overall_pass <- as_logical_strict(
    ecms_overall$pass, "ECMS overall extension pass"
  )
  ecms_gate_order <- match(ecms_gate_names, ecms_gate_rows$metric)
  ecms_gate_observed_ordered <- as.numeric(
    ecms_gate_rows$observed[ecms_gate_order]
  )
  ecms_gate_threshold_ordered <- as.numeric(
    ecms_gate_rows$threshold[ecms_gate_order]
  )
  ecms_gate_pass_ordered <- ecms_gate_pass[ecms_gate_order]
  ecms_gate_pass_derived <- is.finite(ecms_gate_observed_ordered) &
    ecms_gate_observed_ordered >= ecms_gate_thresholds
  ecms_overall_observed <- as.numeric(ecms_overall$observed[[1L]])
  ecms_overall_threshold <- as.numeric(ecms_overall$threshold[[1L]])
  fail_if(
    any(!is.finite(ecms_gate_observed_ordered)) ||
      any(ecms_gate_rows$direction != ">=") ||
      any(ecms_gate_threshold_ordered != ecms_gate_thresholds) ||
      !identical(ecms_gate_pass_ordered, ecms_gate_pass_derived) ||
      !is.finite(ecms_overall_observed) ||
      ecms_overall_observed != as.numeric(all(ecms_gate_pass_derived)) ||
      !is.finite(ecms_overall_threshold) || ecms_overall_threshold != 1 ||
      ecms_overall$direction[[1L]] != "all_four_prelocked_gates" ||
      ecms_extension_calibration_pass != all(ecms_gate_pass_derived) ||
      ecms_extension_calibration_pass != ecms_overall_pass,
    "ECMS calibration 缺少预锁定四门禁、overall 或与 manifest 不一致。"
  )
  ecms_expected_qa_ids <- c(
    "locked_commit", "model_bytes", "model_manifest_full_sha",
    "upstream_artifact_manifests", "model_object_set",
    "randomforest_structure", "feature_count",
    "feature_unique", "feature_list_sha256", "builtin_tcga_shape",
    "builtin_tcga_gene_scaling", "builtin_prediction_counts",
    "five_layer_core_count",
    "official_id_mapping", "gdc_tpm_assay", "ensembl_strip_unique",
    "feature_coverage", "feature_nonzero_variance", "gdc_projection_shape",
    "probability_row_sums", "official_anchor_precedence", "extension_gate",
    "margin_boundary", "pseudo_label_prohibition",
    "clinical_classifier_boundary", "level_factor_qc_identity",
    "level_factor_flag_propagation", "level_factor_evidence_ceiling"
  )
  ecms_qa_hard <- as_logical_strict(
    ecms_qa$hard_gate, "ECMS QA hard_gate"
  )
  ecms_qa_pass <- as_logical_strict(
    ecms_qa$pass, "ECMS QA pass"
  )
  fail_if(
    nrow(ecms_qa) != length(ecms_expected_qa_ids) ||
      anyDuplicated(ecms_qa$check_id) > 0L ||
      !setequal(ecms_qa$check_id, ecms_expected_qa_ids) ||
      sum(!ecms_qa_hard) != 1L ||
      ecms_qa$check_id[!ecms_qa_hard] != "extension_gate" ||
      any(ecms_qa_hard & !ecms_qa_pass) ||
      ecms_qa_pass[match("extension_gate", ecms_qa$check_id)] !=
        ecms_extension_calibration_pass,
    "ECMS QA 必需 check_id、hard gate 或 extension gate 状态失败。"
  )
  ecms_primary_patient_count <- as.integer(
    ecms_manifest_primary_numeric[[1L]]
  )
  ecms_anchor_patient_count <- sum(ecms_anchor_flag)
  expected_nonanchor_source <- if (ecms_extension_calibration_pass) {
    "gdc_tpm_extension_after_overlap_calibration_pass"
  } else {
    "gdc_tpm_extension_conditional_not_primary"
  }
  fail_if(
    !ecms_primary_patient_count %in% c(78L, 94L) ||
      ecms_primary_patient_count != sum(ecms_eligible_flag) ||
      (ecms_extension_calibration_pass && ecms_primary_patient_count != 94L) ||
      (!ecms_extension_calibration_pass && ecms_primary_patient_count != 78L),
    "ECMS 扩展门禁与 primary 患者数不一致。"
  )
  fail_if(
    any(ecms_patient$label_source[ecms_anchor_flag] !=
          "locked_model_prediction_on_repository_bundled_tcga78_anchor") ||
      any(ecms_patient$label_source[!ecms_anchor_flag] !=
            expected_nonanchor_source) ||
      any(ecms_eligible_flag[ecms_anchor_flag] != TRUE) ||
      any(ecms_eligible_flag[!ecms_anchor_flag] !=
            ecms_extension_calibration_pass),
    "ECMS anchor/扩展 label_source 与 primary eligibility 不一致。"
  )
  ecms_model_commit <- ecms_manifest$model_commit[[1L]]
  ecms_projection_signature <- digest(
    list(
      artifact_sha256 = ecms_manifest_rows$sha256,
      model_commit = ecms_model_commit,
      extension_calibration_pass = ecms_extension_calibration_pass,
      primary_patient_count = ecms_primary_patient_count
    ),
    algo = "sha256"
  )
}
source_manifest_rows <- rbindlist(list(
  multiassay_manifest_rows,
  driver_manifest_rows,
  heterogeneity_manifest_rows,
  state_manifest_rows
), use.names = TRUE)

message("[1/8] 读取冻结特征、参考因子和 driver-event 关联")
analysis_sets <- read_tsv(required_inputs$analysis_sets)
driver_screen <- read_tsv(required_inputs$driver_screen)
feature_manifest <- read_tsv(required_inputs$feature_manifest)
reference_factor_table <- read_tsv(required_inputs$reference_factors)
level_factor_qc <- read_tsv(required_inputs$level_factor_qc)
patient_events <- read_tsv(required_inputs$patient_events)
original_associations <- read_tsv(required_inputs$original_associations)
original_network_edges <- read_tsv(required_inputs$original_network_edges)

require_columns(analysis_sets, c(
  "patient_id", "analysis_set", "included"
), "tcga_escc_multiassay_analysis_sets.tsv")
require_columns(driver_screen, c(
  "gene_id", "gene_name", "decision", "primary_candidate_route"
), "tcga_escc_driver_candidate_screen.tsv")
require_columns(feature_manifest, c(
  "view", "feature_id", "feature_name", "selection_metric", "selection_value",
  "rank", "preprocessing"
), "tcga_escc_heterogeneity_feature_manifest.tsv")
require_columns(reference_factor_table, c(
  "patient_id", paste0("Factor", 1:8)
), "tcga_escc_mofa_factor_scores.tsv")
require_columns(level_factor_qc, c(
  "factor", "view", "level_factor_soft_flag"
), "tcga_escc_mofa_level_factor_qc.tsv")
require_columns(patient_events, c(
  "patient_id", "gene_id", "gene_name", "mutation", "relative_cnv",
  "amplification", "homozygous_deletion"
), "tcga_escc_strong_driver_patient_events.tsv")
require_columns(original_associations, c(
  "gene_id", "gene_name", "event_type", "target_type", "target",
  "n_complete", "n_event", "n_reference", "effect_measure", "effect",
  "p_value", "q_value", "analysis_eligible", "association_status",
  "level_factor_soft_flag", "predictor_feature_selected_in_target_model",
  "target_model_overlap_class"
), "tcga_escc_driver_state_associations.tsv")
require_columns(original_network_edges, c(
  "edge_id", "source_gene_id", "source_gene_name", "event_type",
  "target_node", "target_type", "effect_measure", "effect", "p_value",
  "q_value", "association_status", "level_factor_soft_flag"
), "tcga_escc_driver_state_network_edges.tsv")

patients <- reference_factor_table$patient_id
reference_factors <- as.matrix(reference_factor_table[, -"patient_id"])
rownames(reference_factors) <- patients
fail_if(length(patients) != 94L || uniqueN(patients) != 94L,
        "参考 MOFA 因子必须覆盖 94 位唯一患者。")
fail_if(!identical(colnames(reference_factors), paste0("Factor", 1:8)),
        "参考因子名称必须冻结为 Factor1–Factor8。")
fail_if(any(!is.finite(reference_factors)), "参考 Factor1–8 患者得分含非有限值。")
five_layer_patients <- analysis_sets[
  analysis_set == "five_layer_core" & included == TRUE,
  patient_id
]
fail_if(!identical(patients, five_layer_patients),
        "参考因子患者顺序与 script 14 five-layer core 冻结顺序不一致。")
fail_if(nrow(patient_events) != 94L * 12L ||
          uniqueN(patient_events$patient_id) != 94L ||
          uniqueN(patient_events$gene_id) != 12L,
        "冻结 patient-event 表必须为 94×12。")
fail_if(!setequal(patient_events$patient_id, patients),
        "冻结 patient-event 患者与参考因子患者不一致。")

expected_feature_counts <- c(
  RNA = 1500L,
  miRNA = 300L,
  HM450 = 2000L,
  Mutation = 50L,
  CNV = 550L
)
observed_feature_counts <- feature_manifest[, .N, by = view]
fail_if(
  !identical(
    observed_feature_counts$N[match(names(expected_feature_counts), observed_feature_counts$view)],
    unname(expected_feature_counts)
  ),
  "冻结五视图特征数与 script 14 不一致。"
)
fail_if(anyDuplicated(feature_manifest[, .(view, feature_id)]) > 0L,
        "冻结 feature manifest 存在 view-feature 重复。")
cross_view_duplicate_ids <- feature_manifest[, .(
  view_count = uniqueN(view),
  views = collapse_unique(view)
), by = feature_id][view_count > 1L]

frozen_associations <- original_associations[
  target_type == "MOFA_factor" & analysis_eligible == TRUE
]
fail_if(nrow(frozen_associations) != 200L,
        "eligible driver-event→MOFA 关系必须硬冻结为 200 条。")
fail_if(any(!frozen_associations$target %chin% colnames(reference_factors)),
        "冻结关联出现参考 Factor1–8 之外的目标。")
fail_if(anyDuplicated(frozen_associations[, .(gene_id, event_type, target)]) > 0L,
        "冻结关联 gene-event-factor 键重复。")

factor_level_flags <- level_factor_qc[, .(
  level_factor_soft_flag_frozen = any(as.logical(level_factor_soft_flag))
), by = .(factor)]
fail_if(nrow(factor_level_flags) != 8L ||
          !setequal(factor_level_flags$factor, colnames(reference_factors)),
        "level-factor QC 必须唯一覆盖 Factor1–Factor8。")
association_level_expected <- factor_level_flags$level_factor_soft_flag_frozen[
  match(
    original_associations[target_type == "MOFA_factor", target],
    factor_level_flags$factor
  )
]
fail_if(anyNA(association_level_expected) ||
          !identical(
            as.logical(original_associations[
              target_type == "MOFA_factor", level_factor_soft_flag
            ]),
            association_level_expected
          ),
        "original association 的 level-factor flag 与正式 QC 不一致。")

network_statuses <- c(
  "within_tcga_supported",
  "within_tcga_conditional",
  "within_tcga_conditional_level_factor_qc",
  "directional_exploratory",
  "directional_exploratory_level_factor_qc"
)
original_mofa_network <- original_network_edges[target_type == "MOFA_factor"]
expected_mofa_network <- original_associations[
  target_type == "MOFA_factor" & analysis_eligible == TRUE &
    association_status %chin% network_statuses
]
fail_if(nrow(original_mofa_network) != 41L || nrow(expected_mofa_network) != 41L,
        "正式 MOFA network edge 与 association 来源均必须为 41 条。")
fail_if(anyDuplicated(original_mofa_network$edge_id) > 0L ||
          anyDuplicated(original_mofa_network[, .(
            source_gene_id, event_type, target_node
          )]) > 0L,
        "正式 MOFA network edge ID 或关系键重复。")

numeric_vectors_equal <- function(x, y, tolerance = 1e-12) {
  if (length(x) != length(y)) return(FALSE)
  comparable <- (is.na(x) & is.na(y)) |
    (is.finite(x) & is.finite(y) & abs(x - y) <= tolerance)
  isTRUE(all(comparable))
}

network_check <- merge(
  expected_mofa_network[, .(
    gene_id,
    gene_name,
    event_type,
    target,
    association_effect_measure = effect_measure,
    association_effect = effect,
    association_p_value = p_value,
    association_q_value = q_value,
    association_status,
    association_level_factor_soft_flag = level_factor_soft_flag
  )],
  original_mofa_network[, .(
    original_edge_id = edge_id,
    gene_id = source_gene_id,
    gene_name = source_gene_name,
    event_type,
    target = target_node,
    network_effect_measure = effect_measure,
    network_effect = effect,
    network_p_value = p_value,
    network_q_value = q_value,
    network_status = association_status,
    network_level_factor_soft_flag = level_factor_soft_flag
  )],
  by = c("gene_id", "gene_name", "event_type", "target"),
  all = TRUE,
  sort = FALSE
)
fail_if(nrow(network_check) != 41L || anyNA(network_check$original_edge_id),
        "正式 MOFA network edge 与 association 关系键无法一一对应。")
fail_if(!identical(
  network_check$association_effect_measure,
  network_check$network_effect_measure
), "正式 MOFA network edge 的 effect_measure 与 association 不一致。")
fail_if(!numeric_vectors_equal(
  network_check$association_effect, network_check$network_effect
), "正式 MOFA network edge 的 effect 与 association 不一致。")
fail_if(!numeric_vectors_equal(
  network_check$association_p_value, network_check$network_p_value
), "正式 MOFA network edge 的 P 值与 association 不一致。")
fail_if(!numeric_vectors_equal(
  network_check$association_q_value, network_check$network_q_value
), "正式 MOFA network edge 的 q 值与 association 不一致。")
fail_if(!identical(network_check$association_status, network_check$network_status),
        "正式 MOFA network edge 的 status 与 association 不一致。")
fail_if(!identical(
  as.logical(network_check$association_level_factor_soft_flag),
  as.logical(network_check$network_level_factor_soft_flag)
), "正式 MOFA network edge 的 level-factor flag 与 association 不一致。")

frozen_associations <- merge(
  frozen_associations,
  network_check[, .(gene_id, event_type, target, original_edge_id)],
  by = c("gene_id", "event_type", "target"),
  all.x = TRUE,
  sort = FALSE
)

main_roles <- data.table(
  gene_name = c("PIK3CA", "NFE2L2", "GNAS", "ZNF750"),
  narrative_role = c(
    "known_positive_anchor",
    "known_positive_anchor",
    "exploratory_candidate",
    "boundary_context_reverse"
  ),
  role_interpretation = c(
    "已知阳性锚点仅用于检验泄漏控制能否保留，不预设必须阳性",
    "已知阳性锚点仅用于检验泄漏控制能否保留，不预设必须阳性",
    "探索性候选；未入模特征的 no-op 不得伪装成敏感性支持",
    "边界/情境反向候选；即使保留也不升级为普适机制"
  )
)
strong_driver <- driver_screen[decision == "strong_patient_level_candidate"]
main_roles <- merge(
  main_roles,
  strong_driver[, .(gene_id, gene_name, primary_candidate_route)],
  by = "gene_name",
  all.x = TRUE,
  sort = FALSE
)
fail_if(anyNA(main_roles$gene_id), "主叙事候选未全部属于当前 strong driver。")

mutation_feature_ids <- feature_manifest[view == "Mutation", feature_id]
cnv_feature_ids <- feature_manifest[view == "CNV", feature_id]
main_roles[, `:=`(
  mutation_feature_in_reference_model = gene_id %chin% mutation_feature_ids,
  cnv_feature_in_reference_model = gene_id %chin% cnv_feature_ids
)]

scenario_table <- rbindlist(list(
  data.table(
    scenario_order = 1L,
    scenario_id = "baseline",
    scenario_class = "baseline_multi_seed",
    target_gene_id = NA_character_,
    target_gene_name = NA_character_,
    drop_mutation_view = FALSE,
    drop_cnv_view = FALSE
  ),
  data.table(
    scenario_order = 2:4,
    scenario_id = c(
      "drop_mutation_view", "drop_cnv_view", "drop_both_event_views"
    ),
    scenario_class = c(
      "global_view_drop", "global_view_drop", "global_view_drop"
    ),
    target_gene_id = NA_character_,
    target_gene_name = NA_character_,
    drop_mutation_view = c(TRUE, FALSE, TRUE),
    drop_cnv_view = c(FALSE, TRUE, TRUE)
  ),
  main_roles[, .(
    scenario_order = 4L + seq_len(.N),
    scenario_id = paste0("drop_gene_event_features_", gene_name),
    scenario_class = "gene_event_feature_drop",
    target_gene_id = gene_id,
    target_gene_name = gene_name,
    drop_mutation_view = FALSE,
    drop_cnv_view = FALSE
  )]
), use.names = TRUE)

scenario_table <- merge(
  scenario_table,
  main_roles[, .(
    target_gene_id = gene_id,
    target_gene_name = gene_name,
    mutation_feature_in_reference_model,
    cnv_feature_in_reference_model
  )],
  by = c("target_gene_id", "target_gene_name"),
  all.x = TRUE,
  sort = FALSE
)
scenario_table[is.na(mutation_feature_in_reference_model),
               mutation_feature_in_reference_model := FALSE]
scenario_table[is.na(cnv_feature_in_reference_model),
               cnv_feature_in_reference_model := FALSE]
scenario_table[, mutation_action := fcase(
  drop_mutation_view, "view_removed",
  scenario_class == "gene_event_feature_drop" & mutation_feature_in_reference_model,
    "gene_feature_removed",
  scenario_class == "gene_event_feature_drop",
    "no_op_gene_feature_not_in_reference_model",
  default = "retained"
)]
scenario_table[, cnv_action := fcase(
  drop_cnv_view, "view_removed",
  scenario_class == "gene_event_feature_drop" & cnv_feature_in_reference_model,
    "gene_feature_removed",
  scenario_class == "gene_event_feature_drop",
    "no_op_gene_feature_not_in_reference_model",
  default = "retained"
)]
scenario_table[, scenario_no_op :=
  scenario_class == "gene_event_feature_drop" &
    !mutation_feature_in_reference_model & !cnv_feature_in_reference_model]
scenario_table[, expected_model_action := fifelse(
  scenario_no_op, "reuse_seed_matched_baseline_no_op", "fit_new_mofa_model"
)]
scenario_table[, requested_event_types := vapply(
  target_gene_name,
  function(gene) {
    if (is.na(gene)) return("all_analysis_eligible_events")
    collapse_unique(frozen_associations[gene_name == gene, event_type])
  },
  FUN.VALUE = character(1)
)]
scenario_table[, planned_views := vapply(seq_len(.N), function(i) {
  views <- names(expected_feature_counts)
  if (drop_mutation_view[[i]]) views <- setdiff(views, "Mutation")
  if (drop_cnv_view[[i]]) views <- setdiff(views, "CNV")
  paste(views, collapse = ";")
}, FUN.VALUE = character(1))]
scenario_table[, planned_feature_counts := vapply(seq_len(.N), function(i) {
  counts <- expected_feature_counts
  if (drop_mutation_view[[i]]) counts <- counts[names(counts) != "Mutation"]
  if (drop_cnv_view[[i]]) counts <- counts[names(counts) != "CNV"]
  if (scenario_class[[i]] == "gene_event_feature_drop") {
    if (mutation_feature_in_reference_model[[i]]) counts[["Mutation"]] <- counts[["Mutation"]] - 1L
    if (cnv_feature_in_reference_model[[i]]) counts[["CNV"]] <- counts[["CNV"]] - 1L
  }
  paste(names(counts), counts, sep = "=", collapse = ";")
}, FUN.VALUE = character(1))]
scenario_table[, feature_namespace_handling := paste(
  "view-qualified original feature IDs are removed before create_mofa;",
  "MOFA2 cross-view suffixes are never used for no-op detection"
)]
setorder(scenario_table, scenario_order)

model_plan <- rbindlist(lapply(seeds, function(current_seed) {
  output <- copy(scenario_table)
  output[, seed := current_seed]
  output
}), use.names = TRUE)
model_plan[, model_id := paste0(scenario_id, "__seed_", seed)]
setorder(model_plan, seed, scenario_order)
model_plan[, `:=`(
  model_fitted = FALSE,
  model_reused = FALSE,
  run_status = "planned_not_run",
  checkpoint_path = "",
  checkpoint_signature = "",
  factor_count = NA_integer_,
  mofa_status = "",
  model_completed = FALSE,
  model_converged = FALSE,
  training_iteration_count = NA_integer_,
  configured_maxiter = maxiter,
  elbo_finite_record_count = NA_integer_,
  final_elbo_record_index = NA_integer_,
  final_elbo = NA_real_,
  last_relative_elbo_change = NA_real_,
  tail_elbo_non_decreasing_fraction = NA_real_,
  final_factor_count_from_training_stats = NA_integer_,
  convergence_rule = "",
  runtime_seconds = NA_real_
)]

expected_planned_models <- nrow(model_plan)
expected_fitted_models <- model_plan[scenario_no_op == FALSE, .N]
expected_no_op_reuses <- model_plan[scenario_no_op == TRUE, .N]
fail_if(expected_planned_models != length(seeds) * 8L,
        "模型计划必须为每个 seed 8 个场景。")

if (fields_only) {
  message("FIELDS_ONLY_OK")
  message("冻结特征：", paste(names(expected_feature_counts), expected_feature_counts,
                               sep = "=", collapse = "；"))
  message("冻结可检验关联：", nrow(frozen_associations), " 条。")
  message("模型计划：", expected_planned_models, " 行；预计训练 ",
          expected_fitted_models, " 个；no-op 复用 ", expected_no_op_reuses, " 个。")
  message("主叙事特征状态：", paste0(
    main_roles$gene_name,
    "[Mutation=", main_roles$mutation_feature_in_reference_model,
    ",CNV=", main_roles$cnv_feature_in_reference_model, "]",
    collapse = "；"
  ))
  message("跨视图重复 feature_id：", nrow(cross_view_duplicate_ids),
          " 个；按 view-qualified manifest 处理。")
  quit(save = "no", status = 0L)
}

message("[2/8] 由正式 feature manifest 重建 script 14 冻结五视图")
mae <- readRDS(required_inputs$mae)
validObject(mae)
segments <- as.data.table(readRDS(required_inputs$segments))

extract_rows <- function(se, id_column, feature_ids, assay_name, view_name) {
  rd <- as.data.table(as.data.frame(rowData(se)))
  require_columns(rd, id_column, paste0(view_name, " rowData"))
  index <- match(feature_ids, rd[[id_column]])
  fail_if(anyNA(index), paste(view_name, " 冻结特征无法从 MultiAssayExperiment 回填。"))
  matrix <- assay(se, assay_name)[index, patients, drop = FALSE]
  rownames(matrix) <- feature_ids
  fail_if(!identical(colnames(matrix), patients), paste(view_name, " 患者顺序不一致。"))
  matrix
}

feature_ids_by_view <- lapply(names(expected_feature_counts), function(view_name) {
  feature_manifest[view == view_name][order(rank), feature_id]
})
names(feature_ids_by_view) <- names(expected_feature_counts)

rna_matrix <- extract_rows(
  experiments(mae)[["RNA"]], "gene_id", feature_ids_by_view$RNA,
  "tmm_logcpm", "RNA"
)
mirna_matrix <- extract_rows(
  experiments(mae)[["miRNA"]], "miRNA_ID", feature_ids_by_view$miRNA,
  "tmm_logcpm", "miRNA"
)
hm450_beta <- extract_rows(
  experiments(mae)[["HM450"]], "probe_id", feature_ids_by_view$HM450,
  "beta", "HM450"
)
hm450_beta <- impute_row_median(hm450_beta)
hm450_beta[hm450_beta < 0.001] <- 0.001
hm450_beta[hm450_beta > 0.999] <- 0.999
hm450_matrix <- log2(hm450_beta / (1 - hm450_beta))
mutation_matrix <- as.matrix(extract_rows(
  experiments(mae)[["Mutation"]], "gene_id", feature_ids_by_view$Mutation,
  "binary", "Mutation"
))

cnv_se <- experiments(mae)[["CNV_gene"]]
cnv_absolute <- extract_rows(
  cnv_se, "gene_id", feature_ids_by_view$CNV,
  "copy_number", "CNV"
)
require_columns(segments, c(
  "Chromosome", "Start", "End", "Copy_Number", "patient_id", "segment_length"
), "tcga_escc_ascat2_segments_long.rds")
segment_length_recomputed <- as.numeric(segments$End) -
  as.numeric(segments$Start) + 1
fail_if(anyNA(segments$segment_length) || anyNA(segment_length_recomputed) ||
          any(as.numeric(segments$segment_length) != segment_length_recomputed),
        "冻结 segment_length 与 End-Start+1 未逐行全等。")
segments <- segments[Chromosome %chin% paste0("chr", 1:22)]
ploidy <- segments[, .(
  ploidy_proxy = weighted_ploidy_proxy(Copy_Number, segment_length)
), by = patient_id]
ploidy_vector <- ploidy$ploidy_proxy[match(patients, ploidy$patient_id)]
fail_if(any(!is.finite(ploidy_vector)) || any(ploidy_vector <= 0),
        "five-layer core 存在无效 ploidy proxy。")
cnv_relative <- log2(
  (cnv_absolute + 0.5) /
    matrix(
      ploidy_vector + 0.5,
      nrow = nrow(cnv_absolute),
      ncol = ncol(cnv_absolute),
      byrow = TRUE
    )
)
cnv_relative <- impute_row_median(cnv_relative)

base_views <- list(
  RNA = row_zscore(rna_matrix),
  miRNA = row_zscore(mirna_matrix),
  HM450 = row_zscore(hm450_matrix),
  Mutation = mutation_matrix,
  CNV = row_zscore(cnv_relative)
)
fail_if(
  any(vapply(base_views, function(x) !identical(colnames(x), patients), logical(1))),
  "重建五视图患者顺序不一致。"
)
fail_if(
  !identical(vapply(base_views, nrow, integer(1)), expected_feature_counts),
  "重建五视图维度与冻结 manifest 不一致。"
)
fail_if(any(vapply(base_views, function(x) any(!is.finite(x)), logical(1))),
        "重建五视图仍含非有限值。")
fail_if(!all(base_views$Mutation %in% c(0, 1)), "Mutation 视图不是二元矩阵。")

frozen_feature_signature <- digest(
  list(
    patients = patients,
    features = feature_ids_by_view,
    dimensions = vapply(base_views, dim, integer(2))
  ),
  algo = "sha256"
)
base_view_numeric_hashes <- vapply(
  base_views,
  digest,
  FUN.VALUE = character(1),
  algo = "sha256",
  serialize = TRUE
)
base_views_numeric_hash <- digest(
  list(
    view_names = names(base_views),
    view_hashes = base_view_numeric_hashes
  ),
  algo = "sha256"
)
execution_script_sha256 <- digest(
  script_path, algo = "sha256", file = TRUE, serialize = FALSE
)
# 2026-07-12 首轮完整训练后只修正 model_plan 的 data.table 状态回写。
# 该修正不改变任何 MOFA 输入、参数或训练代码，因此模型检查点继续绑定到
# 首轮训练脚本的冻结 SHA；正式 artifact manifest 仍记录当前执行脚本 SHA。
model_training_contract_script_sha256 <-
  "6f4040a27d581ca9023ed2b85f622d1c0f15da09a58a0e085605a67b7d90c5d3"
r_version <- R.version.string
mofa2_version <- as.character(packageVersion("MOFA2"))
mofapy2_version <- query_mofapy2_version()

reference_model <- readRDS(required_inputs$reference_model)
fail_if(!inherits(reference_model, "MOFA"), "正式 reference MOFA RDS 类型异常。")
reference_model_scores <- get_factors(
  reference_model, factors = "all"
)[[1L]]
reference_model_scores <- as.matrix(reference_model_scores[patients, , drop = FALSE])
fail_if(!identical(dim(reference_model_scores), dim(reference_factors)) ||
          !identical(colnames(reference_model_scores), colnames(reference_factors)) ||
          !numeric_vectors_equal(
            as.numeric(reference_model_scores), as.numeric(reference_factors),
            tolerance = 1e-12
          ),
        "正式 reference MOFA RDS 得分与 factor_scores.tsv 不一致。")
reference_common_loadings <- extract_event_blind_loadings(
  reference_model,
  feature_ids_by_view
)
fail_if(!identical(colnames(reference_common_loadings), colnames(reference_factors)),
        "reference common loadings 未覆盖 Factor1–Factor8。")
reference_configured_maxiter <- suppressWarnings(as.integer(
  reference_model@training_options$maxiter
))
fail_if(length(reference_configured_maxiter) != 1L ||
          is.na(reference_configured_maxiter) ||
          reference_configured_maxiter < 1L,
        "正式 reference MOFA 缺少可追溯的自身 maxiter 训练配置。")
reference_training_diagnostics <- extract_training_diagnostics(
  reference_model, reference_configured_maxiter
)
fail_if(!isTRUE(reference_training_diagnostics$model_completed) ||
          !isTRUE(reference_training_diagnostics$model_converged),
        "正式 reference MOFA 未通过预锁定完成/收敛判定。")

model_input_signature <- digest(
  list(
    source_sha256 = source_manifest_rows$sha256,
    frozen_feature_signature = frozen_feature_signature,
    base_view_numeric_hashes = base_view_numeric_hashes,
    base_views_numeric_hash = base_views_numeric_hash,
    r_version = r_version,
    mofa2_version = mofa2_version,
    mofapy2_version = mofapy2_version,
    execution_script_sha256 = model_training_contract_script_sha256,
    mofa_options = list(
      num_factors = 10L,
      maxiter = maxiter,
      convergence_mode = "medium",
      drop_factor_threshold = 0.01,
      spikeslab_weights = TRUE,
      ard_weights = TRUE,
      scale_views = FALSE,
      center_groups = TRUE,
      use_float32 = TRUE
    )
  ),
  algo = "sha256"
)
input_signature <- digest(
  list(
    model_input_signature = model_input_signature,
    external_ecms_projection_signature = ecms_projection_signature,
    audit_options = list(
      seeds = seeds,
      alignment_method = alignment_method,
      alignment_min_abs_correlation = alignment_min_abs_correlation,
      alignment_min_margin = alignment_min_margin,
      alignment_high_confidence_override = alignment_high_confidence_override
    )
  ),
  algo = "sha256"
)
model_plan[, `:=`(
  model_input_signature = model_input_signature,
  audit_input_signature = input_signature,
  base_view_numeric_hashes = paste(
    names(base_view_numeric_hashes), base_view_numeric_hashes,
    sep = "=", collapse = ";"
  ),
  base_views_numeric_hash = base_views_numeric_hash,
  r_version = r_version,
  mofa2_version = mofa2_version,
  mofapy2_version = mofapy2_version,
  execution_script_sha256 = execution_script_sha256,
  model_training_contract_script_sha256 =
    model_training_contract_script_sha256
)]

if (validate_only) {
  message("VALIDATE_ONLY_OK")
  message("冻结视图已重建：", paste(
    names(base_views), vapply(base_views, nrow, integer(1)), sep = "=",
    collapse = "；"
  ))
  message("患者：", length(patients), "；输入 signature：", input_signature)
  message("base_views numeric hash：", base_views_numeric_hash)
  message("R/MOFA2/mofapy2：", r_version, " / ", mofa2_version,
          " / ", mofapy2_version)
  message("模型计划：", expected_planned_models, " 行；预计训练 ",
          expected_fitted_models, " 个；no-op 复用 ", expected_no_op_reuses, " 个。")
  quit(save = "no", status = 0L)
}

rm(
  mae, rna_matrix, mirna_matrix, hm450_beta, hm450_matrix,
  mutation_matrix, cnv_absolute, cnv_relative, cnv_se, segments, ploidy,
  reference_model, reference_model_scores
)
gc(verbose = FALSE)

dir_create(c(work_intermediate_dir, work_checks_dir), recurse = TRUE)
stage_dir <- tempfile(
  pattern = ".tcga_driver_state_leakage_",
  tmpdir = work_intermediate_dir
)
dir_create(stage_dir)
# 成功后显式删除；失败时保留在 _work/intermediate/ 供审计。
cache_dir <- file.path(
  work_intermediate_dir,
  "tcga_driver_state_leakage_cache"
)
dir_create(cache_dir, recurse = TRUE)

stage_tsv <- function(object, filename) {
  fwrite(object, file.path(stage_dir, filename), sep = "\t", quote = FALSE, na = "")
}

atomic_save_rds <- function(object, path) {
  temporary <- tempfile(
    pattern = paste0(".", path_file(path), "."),
    tmpdir = path_dir(path),
    fileext = ".tmp.rds"
  )
  keep_temporary <- TRUE
  on.exit({
    if (keep_temporary && file_exists(temporary)) file_delete(temporary)
  }, add = TRUE)
  saveRDS(object, temporary, compress = "xz")
  check <- readRDS(temporary)
  fail_if(is.null(check$checkpoint_signature), "checkpoint 临时文件回读失败。")
  fail_if(
    !identical(check$checkpoint_signature, object$checkpoint_signature),
    "checkpoint 临时文件身份回读不一致。"
  )
  renamed <- file.rename(temporary, path)
  fail_if(!renamed, paste0("checkpoint 原子替换失败：", path))
  keep_temporary <- FALSE
  invisible(path)
}

atomic_publish_file <- function(source, destination) {
  temporary <- tempfile(
    pattern = paste0(".", path_file(destination), "."),
    tmpdir = path_dir(destination),
    fileext = ".tmp"
  )
  keep_temporary <- TRUE
  on.exit({
    if (keep_temporary && file_exists(temporary)) file_delete(temporary)
  }, add = TRUE)
  file_copy(source, temporary, overwrite = TRUE)
  fail_if(
    as.numeric(file_info(source)$size) != as.numeric(file_info(temporary)$size) ||
      digest(source, algo = "sha256", file = TRUE, serialize = FALSE) !=
        digest(temporary, algo = "sha256", file = TRUE, serialize = FALSE),
    paste0("发布临时副本校验失败：", path_file(destination))
  )
  renamed <- file.rename(temporary, destination)
  fail_if(!renamed, paste0("正式 artifact 原子替换失败：", destination))
  keep_temporary <- FALSE
  invisible(destination)
}

make_scenario_views <- function(base, plan_row) {
  views <- lapply(base, function(x) x)
  if (isTRUE(plan_row$drop_mutation_view)) views$Mutation <- NULL
  if (isTRUE(plan_row$drop_cnv_view)) views$CNV <- NULL
  if (plan_row$scenario_class == "gene_event_feature_drop") {
    gene_id <- plan_row$target_gene_id
    if (!is.null(views$Mutation) && gene_id %in% rownames(views$Mutation)) {
      views$Mutation <- views$Mutation[rownames(views$Mutation) != gene_id, , drop = FALSE]
    }
    if (!is.null(views$CNV) && gene_id %in% rownames(views$CNV)) {
      views$CNV <- views$CNV[rownames(views$CNV) != gene_id, , drop = FALSE]
    }
  }
  fail_if(length(views) < 3L, "敏感性场景剩余视图少于 3 个。")
  fail_if(any(vapply(views, nrow, integer(1)) < 2L), "敏感性场景出现空或单特征视图。")
  views
}

fit_mofa_components <- function(views, seed, model_id) {
  mofa_object <- create_mofa(views)
  data_options <- get_default_data_options(mofa_object)
  data_options$scale_views <- FALSE
  data_options$center_groups <- TRUE
  data_options$use_float32 <- TRUE
  model_options <- get_default_model_options(mofa_object)
  model_options$num_factors <- 10L
  likelihood_map <- c(
    RNA = "gaussian",
    miRNA = "gaussian",
    HM450 = "gaussian",
    Mutation = "bernoulli",
    CNV = "gaussian"
  )
  model_options$likelihoods <- likelihood_map[names(views)]
  model_options$spikeslab_weights <- TRUE
  model_options$ard_weights <- TRUE
  training_options <- get_default_training_options(mofa_object)
  training_options$maxiter <- maxiter
  training_options$convergence_mode <- "medium"
  training_options$drop_factor_threshold <- 0.01
  training_options$seed <- as.integer(seed)
  training_options$verbose <- FALSE
  mofa_object <- prepare_mofa(
    mofa_object,
    data_options = data_options,
    model_options = model_options,
    training_options = training_options
  )
  hdf5_path <- file.path(stage_dir, paste0(model_id, ".hdf5"))
  if (file_exists(hdf5_path)) file_delete(hdf5_path)
  fitted <- run_mofa(
    mofa_object,
    outfile = hdf5_path,
    save_data = TRUE,
    use_basilisk = TRUE
  )
  scores <- get_factors(fitted, factors = "all")[[1L]]
  fail_if(nrow(scores) != length(patients) || !all(patients %in% rownames(scores)),
          paste("MOFA 患者维度异常：", model_id))
  scores <- as.matrix(scores[patients, , drop = FALSE])
  fail_if(ncol(scores) < 2L || any(!is.finite(scores)) ||
            is.null(colnames(scores)) || anyDuplicated(colnames(scores)) > 0L,
          paste("MOFA 因子数不足或含非有限值：", model_id))
  common_loadings <- extract_event_blind_loadings(
    fitted,
    feature_ids_by_view
  )
  fail_if(!identical(colnames(common_loadings), colnames(scores)),
          paste("MOFA common loadings 与 factor scores 因子顺序不一致：", model_id))
  training_diagnostics <- extract_training_diagnostics(fitted, maxiter)
  fail_if(!isTRUE(training_diagnostics$model_completed),
          paste("MOFA 未完成有效训练：", model_id))
  if (file_exists(hdf5_path)) file_delete(hdf5_path)
  rm(fitted, mofa_object)
  gc(verbose = FALSE)
  list(
    factor_scores = scores,
    common_loadings = common_loadings,
    training_diagnostics = training_diagnostics
  )
}

alignment_margin <- function(values, chosen_index) {
  chosen <- values[[chosen_index]]
  other <- values[-chosen_index]
  other <- other[is.finite(other)]
  if (!is.finite(chosen)) return(NA_real_)
  if (!length(other)) return(chosen)
  chosen - max(other)
}

align_factors <- function(
    reference_scores,
    candidate_scores,
    reference_loadings,
    candidate_loadings,
    model_id) {
  fail_if(!identical(rownames(reference_loadings), rownames(candidate_loadings)),
          paste("event-blind shared loading 特征顺序不一致：", model_id))
  fail_if(!identical(colnames(reference_loadings), colnames(reference_scores)) ||
            !identical(colnames(candidate_loadings), colnames(candidate_scores)),
          paste("loading 与 score 因子顺序不一致：", model_id))
  # 门禁匹配只使用共同的 RNA/miRNA/HM450 标准化载荷，完全不读取患者事件。
  loading_correlation <- suppressWarnings(cor(
    reference_loadings,
    candidate_loadings,
    use = "pairwise.complete.obs",
    method = "pearson"
  ))
  patient_score_correlation <- suppressWarnings(cor(
    reference_scores,
    candidate_scores,
    use = "pairwise.complete.obs",
    method = alignment_method
  ))
  score <- abs(loading_correlation)
  assignment_score <- score
  assignment_score[!is.finite(assignment_score)] <- 0
  n_reference <- nrow(assignment_score)
  n_candidate <- ncol(assignment_score)
  if (n_reference <= n_candidate) {
    assignment <- solve_LSAP(assignment_score, maximum = TRUE)
    pairs <- data.table(
      reference_index = seq_len(n_reference),
      candidate_index = as.integer(assignment)
    )
  } else {
    assignment <- solve_LSAP(t(assignment_score), maximum = TRUE)
    pairs <- data.table(
      reference_index = as.integer(assignment),
      candidate_index = seq_len(n_candidate)
    )
  }

  matched_rows <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    reference_index <- pairs$reference_index[[i]]
    candidate_index <- pairs$candidate_index[[i]]
    raw <- loading_correlation[reference_index, candidate_index]
    absolute <- abs(raw)
    patient_score_raw <- patient_score_correlation[
      reference_index, candidate_index
    ]
    reference_margin <- alignment_margin(
      score[reference_index, ], candidate_index
    )
    candidate_margin <- alignment_margin(
      score[, candidate_index], reference_index
    )
    minimum_margin <- min(reference_margin, candidate_margin, na.rm = TRUE)
    if (!is.finite(minimum_margin)) minimum_margin <- NA_real_
    low_margin_flag <- !is.finite(minimum_margin) ||
      minimum_margin < alignment_min_margin
    # high-correlation override 也必须通过同一预锁定 margin，不允许绕过。
    reliable <- is.finite(absolute) &&
      (absolute >= alignment_min_abs_correlation ||
         absolute >= alignment_high_confidence_override) &&
      !low_margin_flag
    status <- if (!is.finite(raw)) {
      "assigned_nonfinite_unmatched"
    } else if (absolute < alignment_min_abs_correlation) {
      "assigned_low_correlation_unmatched"
    } else if (low_margin_flag) {
      "assigned_low_margin_unmatched"
    } else if (!reliable) {
      "assigned_unreliable_unmatched"
    } else {
      "reliable_hungarian_match"
    }
    data.table(
      model_id = model_id,
      reference_factor = colnames(reference_scores)[[reference_index]],
      candidate_factor = colnames(candidate_scores)[[candidate_index]],
      raw_shared_loading_correlation = raw,
      absolute_shared_loading_correlation = absolute,
      patient_score_correlation_diagnostic_raw = patient_score_raw,
      patient_score_correlation_diagnostic_aligned =
        patient_score_raw * if (is.finite(raw) && raw < 0) -1 else 1,
      absolute_patient_score_correlation_diagnostic = abs(patient_score_raw),
      patient_score_diagnostic_is_gate = FALSE,
      sign_multiplier = if (is.finite(raw) && raw < 0) -1 else 1,
      sign_flipped = is.finite(raw) && raw < 0,
      reference_margin = reference_margin,
      candidate_margin = candidate_margin,
      minimum_match_margin = minimum_margin,
      low_margin_flag = low_margin_flag,
      match_reliable = reliable,
      matching_status = status,
      loading_alignment_is_event_blind = TRUE,
      loading_alignment_views = "RNA;miRNA;HM450",
      loading_standardization =
        "per-view per-factor center then unit-L2; concatenate shared features",
      reference_factor_count = ncol(reference_scores),
      candidate_factor_count = ncol(candidate_scores),
      loading_alignment_method = "Pearson correlation of standardized loadings",
      patient_score_diagnostic_method = alignment_method,
      minimum_absolute_correlation = alignment_min_abs_correlation,
      minimum_margin_required = alignment_min_margin,
      high_confidence_override = alignment_high_confidence_override,
      override_margin_requirement = alignment_min_margin
    )
  }))

  unmatched_reference <- setdiff(
    seq_len(ncol(reference_scores)), pairs$reference_index
  )
  if (length(unmatched_reference)) {
    matched_rows <- rbindlist(list(
      matched_rows,
      data.table(
        model_id = model_id,
        reference_factor = colnames(reference_scores)[unmatched_reference],
        candidate_factor = NA_character_,
        raw_shared_loading_correlation = NA_real_,
        absolute_shared_loading_correlation = NA_real_,
        patient_score_correlation_diagnostic_raw = NA_real_,
        patient_score_correlation_diagnostic_aligned = NA_real_,
        absolute_patient_score_correlation_diagnostic = NA_real_,
        patient_score_diagnostic_is_gate = FALSE,
        sign_multiplier = NA_real_,
        sign_flipped = NA,
        reference_margin = NA_real_,
        candidate_margin = NA_real_,
        minimum_match_margin = NA_real_,
        low_margin_flag = TRUE,
        match_reliable = FALSE,
        matching_status = "unmatched_reference_factor",
        loading_alignment_is_event_blind = TRUE,
        loading_alignment_views = "RNA;miRNA;HM450",
        loading_standardization =
          "per-view per-factor center then unit-L2; concatenate shared features",
        reference_factor_count = ncol(reference_scores),
        candidate_factor_count = ncol(candidate_scores),
        loading_alignment_method = "Pearson correlation of standardized loadings",
        patient_score_diagnostic_method = alignment_method,
        minimum_absolute_correlation = alignment_min_abs_correlation,
        minimum_margin_required = alignment_min_margin,
        high_confidence_override = alignment_high_confidence_override,
        override_margin_requirement = alignment_min_margin
      )
    ), use.names = TRUE)
  }

  unmatched_candidate <- setdiff(
    seq_len(ncol(candidate_scores)), pairs$candidate_index
  )
  extra_rows <- if (length(unmatched_candidate)) {
    data.table(
      model_id = model_id,
      reference_factor = NA_character_,
      candidate_factor = colnames(candidate_scores)[unmatched_candidate],
      raw_shared_loading_correlation = NA_real_,
      absolute_shared_loading_correlation = NA_real_,
      patient_score_correlation_diagnostic_raw = NA_real_,
      patient_score_correlation_diagnostic_aligned = NA_real_,
      absolute_patient_score_correlation_diagnostic = NA_real_,
      patient_score_diagnostic_is_gate = FALSE,
      sign_multiplier = NA_real_,
      sign_flipped = NA,
      reference_margin = NA_real_,
      candidate_margin = NA_real_,
      minimum_match_margin = NA_real_,
      low_margin_flag = TRUE,
      match_reliable = FALSE,
      matching_status = "unmatched_extra_candidate_factor",
      loading_alignment_is_event_blind = TRUE,
      loading_alignment_views = "RNA;miRNA;HM450",
      loading_standardization =
        "per-view per-factor center then unit-L2; concatenate shared features",
      reference_factor_count = ncol(reference_scores),
      candidate_factor_count = ncol(candidate_scores),
      loading_alignment_method = "Pearson correlation of standardized loadings",
      patient_score_diagnostic_method = alignment_method,
      minimum_absolute_correlation = alignment_min_abs_correlation,
      minimum_margin_required = alignment_min_margin,
      high_confidence_override = alignment_high_confidence_override,
      override_margin_requirement = alignment_min_margin
    )
  } else {
    matched_rows[0]
  }

  aligned <- matrix(
    NA_real_,
    nrow = nrow(reference_scores),
    ncol = ncol(reference_scores),
    dimnames = list(rownames(reference_scores), colnames(reference_scores))
  )
  reliable_rows <- matched_rows[
    !is.na(reference_factor) & match_reliable == TRUE
  ]
  for (i in seq_len(nrow(reliable_rows))) {
    aligned[, reliable_rows$reference_factor[[i]]] <-
      candidate_scores[, reliable_rows$candidate_factor[[i]]] *
      reliable_rows$sign_multiplier[[i]]
  }
  list(
    alignment = rbindlist(list(matched_rows, extra_rows), use.names = TRUE),
    aligned_scores = aligned
  )
}

message("[3/8] 训练或恢复多 seed / 视图删除 / 候选特征删除 MOFA")
model_results <- new.env(parent = emptyenv())
alignment_rows <- vector("list", nrow(model_plan))
association_rows <- vector("list", nrow(model_plan))

for (plan_index in seq_len(nrow(model_plan))) {
  plan <- model_plan[plan_index]
  model_id <- plan$model_id[[1L]]
  seed <- plan$seed[[1L]]
  start_time <- proc.time()[["elapsed"]]
  message(
    "  [", plan_index, "/", nrow(model_plan), "] ", model_id,
    if (plan$scenario_no_op[[1L]]) " (no-op baseline reuse)" else ""
  )

  if (plan$scenario_no_op[[1L]]) {
    baseline_id <- paste0("baseline__seed_", seed)
    fail_if(!exists(baseline_id, envir = model_results, inherits = FALSE),
            paste("no-op 场景找不到 seed-matched baseline：", baseline_id))
    result <- get(baseline_id, envir = model_results, inherits = FALSE)
    model_plan[plan_index, `:=`(
      model_fitted = FALSE,
      model_reused = TRUE,
      run_status = "reused_seed_matched_baseline_no_op",
      checkpoint_path = result$checkpoint_path,
      checkpoint_signature = result$checkpoint_signature,
      factor_count = ncol(result$factor_scores),
      runtime_seconds = proc.time()[["elapsed"]] - start_time
    )]
  } else {
    scenario_signature <- digest(
      list(
        model_input_signature = model_input_signature,
        seed = seed,
        scenario_id = plan$scenario_id[[1L]],
        scenario_class = plan$scenario_class[[1L]],
        target_gene_id = plan$target_gene_id[[1L]],
        mutation_action = plan$mutation_action[[1L]],
        cnv_action = plan$cnv_action[[1L]]
      ),
      algo = "sha256"
    )
    checkpoint_path <- file.path(
      cache_dir,
      paste0(gsub("[^A-Za-z0-9_.-]", "_", model_id), "__",
             substr(scenario_signature, 1L, 16L), ".rds")
    )
    expected_checkpoint_views <- strsplit(
      plan$planned_views[[1L]], ";", fixed = TRUE
    )[[1L]]
    expected_count_fields <- strsplit(
      plan$planned_feature_counts[[1L]], ";", fixed = TRUE
    )[[1L]]
    expected_checkpoint_feature_counts <- setNames(
      as.integer(sub("^[^=]+=", "", expected_count_fields)),
      sub("=.*$", "", expected_count_fields)
    )
    fail_if(anyNA(expected_checkpoint_feature_counts) ||
              !identical(
                names(expected_checkpoint_feature_counts),
                expected_checkpoint_views
              ),
            paste("无法从 model plan 解析预期视图/特征数：", model_id))
    result <- NULL
    if (file_exists(checkpoint_path) && !force_retrain) {
      checkpoint <- readRDS(checkpoint_path)
      valid_checkpoint <- is.list(checkpoint) &&
        identical(checkpoint$checkpoint_signature, scenario_signature) &&
        identical(checkpoint$model_input_signature, model_input_signature) &&
        identical(checkpoint$model_id, model_id) &&
        identical(as.integer(checkpoint$seed), as.integer(seed)) &&
        identical(checkpoint$scenario_id, plan$scenario_id[[1L]]) &&
        identical(checkpoint$view_names, expected_checkpoint_views) &&
        identical(
          checkpoint$feature_counts,
          expected_checkpoint_feature_counts
        ) &&
        identical(checkpoint$base_view_numeric_hashes, base_view_numeric_hashes) &&
        identical(checkpoint$base_views_numeric_hash, base_views_numeric_hash) &&
        identical(checkpoint$r_version, r_version) &&
        identical(checkpoint$mofa2_version, mofa2_version) &&
        identical(checkpoint$mofapy2_version, mofapy2_version) &&
        identical(
          checkpoint$execution_script_sha256,
          model_training_contract_script_sha256
        ) &&
        is.matrix(checkpoint$factor_scores) &&
        identical(rownames(checkpoint$factor_scores), patients) &&
        ncol(checkpoint$factor_scores) >= 2L &&
        !is.null(colnames(checkpoint$factor_scores)) &&
        anyDuplicated(colnames(checkpoint$factor_scores)) == 0L &&
        all(is.finite(checkpoint$factor_scores)) &&
        is.matrix(checkpoint$common_loadings) &&
        identical(colnames(checkpoint$common_loadings),
                  colnames(checkpoint$factor_scores)) &&
        identical(rownames(checkpoint$common_loadings),
                  rownames(reference_common_loadings)) &&
        all(is.finite(checkpoint$common_loadings)) &&
        is.list(checkpoint$training_diagnostics) &&
        isTRUE(checkpoint$training_diagnostics$model_completed)
      fail_if(!valid_checkpoint, paste("checkpoint 结构或签名异常：", checkpoint_path))
      result <- checkpoint
      current_run_status <- "reused_verified_checkpoint"
      current_model_fitted <- FALSE
      current_model_reused <- TRUE
    } else {
      views <- make_scenario_views(base_views, plan)
      components <- fit_mofa_components(views, seed, model_id)
      result <- list(
        checkpoint_signature = scenario_signature,
        model_input_signature = model_input_signature,
        base_view_numeric_hashes = base_view_numeric_hashes,
        base_views_numeric_hash = base_views_numeric_hash,
        r_version = r_version,
        mofa2_version = mofa2_version,
        mofapy2_version = mofapy2_version,
        execution_script_sha256 = model_training_contract_script_sha256,
        checkpoint_generation_script_sha256 = execution_script_sha256,
        model_id = model_id,
        seed = seed,
        scenario_id = plan$scenario_id[[1L]],
        factor_scores = components$factor_scores,
        common_loadings = components$common_loadings,
        training_diagnostics = components$training_diagnostics,
        view_names = names(views),
        feature_counts = vapply(views, nrow, integer(1)),
        generated_date = as.character(Sys.Date()),
        checkpoint_role = "reusable_intermediate_not_current_state"
      )
      atomic_save_rds(result, checkpoint_path)
      current_run_status <- "fitted_and_checkpointed"
      current_model_fitted <- TRUE
      current_model_reused <- FALSE
    }
    result$checkpoint_path <- checkpoint_path
    current_checkpoint_path <- checkpoint_path
    model_plan[plan_index, `:=`(
      model_fitted = current_model_fitted,
      model_reused = current_model_reused,
      run_status = current_run_status,
      checkpoint_path = current_checkpoint_path,
      checkpoint_signature = scenario_signature,
      factor_count = ncol(result$factor_scores),
      runtime_seconds = proc.time()[["elapsed"]] - start_time
    )]
  }

  diagnostics <- result$training_diagnostics
  fail_if(!is.list(diagnostics) || !isTRUE(diagnostics$model_completed),
          paste("模型完成诊断缺失或失败：", model_id))
  model_plan[plan_index, `:=`(
    mofa_status = diagnostics$mofa_status,
    model_completed = diagnostics$model_completed,
    model_converged = diagnostics$model_converged,
    training_iteration_count = diagnostics$training_iteration_count,
    configured_maxiter = diagnostics$configured_maxiter,
    elbo_finite_record_count = diagnostics$elbo_finite_record_count,
    final_elbo_record_index = diagnostics$final_elbo_record_index,
    final_elbo = diagnostics$final_elbo,
    last_relative_elbo_change = diagnostics$last_relative_elbo_change,
    tail_elbo_non_decreasing_fraction =
      diagnostics$tail_elbo_non_decreasing_fraction,
    final_factor_count_from_training_stats =
      diagnostics$final_factor_count_from_training_stats,
    convergence_rule = diagnostics$convergence_rule
  )]

  assign(model_id, result, envir = model_results)
  aligned <- align_factors(
    reference_factors,
    result$factor_scores,
    reference_common_loadings,
    result$common_loadings,
    model_id
  )
  alignment <- aligned$alignment
  alignment <- merge(
    alignment,
    model_plan[plan_index, .(
      model_id,
      seed,
      scenario_order,
      scenario_id,
      scenario_class,
      target_gene_id,
      target_gene_name,
      scenario_no_op,
      expected_model_action,
      mutation_action,
      cnv_action,
      run_status,
      model_completed,
      model_converged,
      training_iteration_count,
      final_elbo
    )],
    by = "model_id",
    all.x = TRUE,
    sort = FALSE
  )
  alignment[, counts_toward_sensitivity :=
    !scenario_no_op & model_completed & model_converged]
  alignment_rows[[plan_index]] <- alignment

  association_list <- vector("list", nrow(frozen_associations))
  for (association_index in seq_len(nrow(frozen_associations))) {
    frozen <- frozen_associations[association_index]
    match_row <- alignment[
      reference_factor == frozen$target[[1L]] & !is.na(reference_factor)
    ]
    fail_if(nrow(match_row) != 1L,
            paste("参考因子对齐行不唯一：", model_id, frozen$target[[1L]]))
    score <- aligned$aligned_scores[, frozen$target[[1L]]]
    gene_event <- patient_events[gene_id == frozen$gene_id[[1L]]]
    gene_event <- gene_event[match(patients, patient_id)]
    fail_if(anyNA(gene_event$patient_id) || !identical(gene_event$patient_id, patients),
            paste("患者事件顺序错误：", frozen$gene_name[[1L]]))
    predictor <- gene_event[[frozen$event_type[[1L]]]]
    test <- if (!match_row$match_reliable[[1L]]) {
      list(
        n_complete = 0L,
        n_event = if (frozen$event_type[[1L]] == "relative_cnv") NA_integer_ else sum(predictor > 0),
        n_reference = if (frozen$event_type[[1L]] == "relative_cnv") NA_integer_ else sum(predictor == 0),
        effect = NA_real_,
        p_value = NA_real_
      )
    } else if (frozen$event_type[[1L]] == "relative_cnv") {
      continuous_test(predictor, score)
    } else {
      binary_test(predictor, score)
    }
    association_list[[association_index]] <- data.table(
      model_id = model_id,
      seed = seed,
      scenario_order = plan$scenario_order[[1L]],
      scenario_id = plan$scenario_id[[1L]],
      scenario_class = plan$scenario_class[[1L]],
      scenario_target_gene_id = plan$target_gene_id[[1L]],
      scenario_target_gene_name = plan$target_gene_name[[1L]],
      scenario_no_op = plan$scenario_no_op[[1L]],
      mutation_action = plan$mutation_action[[1L]],
      cnv_action = plan$cnv_action[[1L]],
      model_completed = model_plan$model_completed[[plan_index]],
      model_converged = model_plan$model_converged[[plan_index]],
      training_iteration_count =
        model_plan$training_iteration_count[[plan_index]],
      final_elbo = model_plan$final_elbo[[plan_index]],
      original_edge_id = frozen$original_edge_id[[1L]],
      gene_id = frozen$gene_id[[1L]],
      gene_name = frozen$gene_name[[1L]],
      event_type = frozen$event_type[[1L]],
      reference_factor = frozen$target[[1L]],
      candidate_factor = match_row$candidate_factor[[1L]],
      factor_match_status = match_row$matching_status[[1L]],
      factor_match_reliable = match_row$match_reliable[[1L]],
      factor_match_abs_loading_correlation =
        match_row$absolute_shared_loading_correlation[[1L]],
      factor_match_minimum_margin = match_row$minimum_match_margin[[1L]],
      factor_match_low_margin_flag = match_row$low_margin_flag[[1L]],
      patient_score_correlation_diagnostic =
        match_row$patient_score_correlation_diagnostic_aligned[[1L]],
      patient_score_diagnostic_is_gate = FALSE,
      factor_sign_flipped = match_row$sign_flipped[[1L]],
      n_complete = test$n_complete,
      n_event = test$n_event,
      n_reference = test$n_reference,
      effect_measure = frozen$effect_measure[[1L]],
      aligned_effect = test$effect,
      p_value = test$p_value,
      p_q_scope = "TCGA_internal_descriptive_only",
      p_q_countable_as_independent_validation = FALSE,
      original_effect = frozen$effect[[1L]],
      original_p_value = frozen$p_value[[1L]],
      original_q_value = frozen$q_value[[1L]],
      original_association_status = frozen$association_status[[1L]],
      original_level_factor_soft_flag = frozen$level_factor_soft_flag[[1L]],
      predictor_feature_selected_in_reference_model =
        frozen$predictor_feature_selected_in_target_model[[1L]],
      original_overlap_class = frozen$target_model_overlap_class[[1L]]
    )
  }
  association_rows[[plan_index]] <- rbindlist(association_list)
}

factor_alignment <- rbindlist(alignment_rows, use.names = TRUE, fill = TRUE)
associations <- rbindlist(association_rows, use.names = TRUE, fill = TRUE)
setorder(
  factor_alignment,
  seed,
  scenario_order,
  reference_factor,
  candidate_factor,
  na.last = TRUE
)

message("[4/8] 重新计算固定 BH 家族并裁决方向/效应保留")
associations[, q_value := {
  output <- rep(NA_real_, .N)
  finite <- which(is.finite(p_value))
  if (length(finite)) output[finite] <- p.adjust(p_value[finite], method = "BH", n = .N)
  output
}, by = .(model_id, event_type)]
associations[, new_association_status_without_bootstrap := fcase(
  !factor_match_reliable, "factor_unmatched_or_unreliable",
  is.finite(q_value) & q_value <= 0.10 & abs(aligned_effect) >= 0.30,
    "within_tcga_supported_without_bootstrap",
  ((is.finite(q_value) & q_value <= 0.20) |
     (is.finite(p_value) & p_value <= 0.05)) & abs(aligned_effect) >= 0.20,
    "within_tcga_conditional_without_bootstrap",
  is.finite(p_value) & p_value <= 0.10 & abs(aligned_effect) >= 0.20,
    "directional_exploratory_without_bootstrap",
  default = "no_clear_internal_association"
)]
associations[, direction_retained := fifelse(
  is.finite(original_effect) & original_effect != 0 &
    is.finite(aligned_effect) & aligned_effect != 0,
  sign(aligned_effect) == sign(original_effect),
  NA
)]
associations[, absolute_effect_ratio := fifelse(
  is.finite(original_effect) & original_effect != 0 & is.finite(aligned_effect),
  abs(aligned_effect) / abs(original_effect),
  NA_real_
)]
associations[, magnitude_retained_50pct := fifelse(
  is.finite(absolute_effect_ratio), absolute_effect_ratio >= 0.50, NA
)]
associations[, original_network_edge := !is.na(original_edge_id)]
associations[, original_level_retained := fcase(
  !original_network_edge, NA,
  original_level_factor_soft_flag, NA,
  !factor_match_reliable | !(direction_retained %in% TRUE) |
    !(magnitude_retained_50pct %in% TRUE), FALSE,
  original_association_status == "within_tcga_supported",
    new_association_status_without_bootstrap == "within_tcga_supported_without_bootstrap",
  original_association_status %chin% c(
    "within_tcga_conditional", "within_tcga_conditional_level_factor_qc"
  ),
    new_association_status_without_bootstrap %chin% c(
      "within_tcga_supported_without_bootstrap",
      "within_tcga_conditional_without_bootstrap"
    ),
  original_association_status %chin% c(
    "directional_exploratory", "directional_exploratory_level_factor_qc"
  ),
    new_association_status_without_bootstrap %chin% c(
      "within_tcga_supported_without_bootstrap",
      "within_tcga_conditional_without_bootstrap",
      "directional_exploratory_without_bootstrap"
    ),
  default = FALSE
)]
associations[, gate_retained := fifelse(
  original_association_status == "within_tcga_supported" &
    original_level_factor_soft_flag == FALSE &
    original_level_retained %in% TRUE &
    new_association_status_without_bootstrap ==
      "within_tcga_supported_without_bootstrap",
  TRUE,
  FALSE
)]
associations[, counts_toward_retention :=
  !scenario_no_op & model_completed & model_converged &
    (
      scenario_class != "gene_event_feature_drop" |
        (!is.na(scenario_target_gene_id) & gene_id == scenario_target_gene_id)
    )]
associations[, no_op_not_counted_as_sensitivity := scenario_no_op]
associations <- merge(
  associations,
  main_roles[, .(gene_id, narrative_role, role_interpretation)],
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)
associations[is.na(narrative_role), `:=`(
  narrative_role = "other_strong_driver",
  role_interpretation = "其他 strong driver 的全局视图删除敏感性"
)]
associations[, conclusion_ceiling := paste(
  "TCGA 同患者表示泄漏敏感性；多 seed/删视图/删特征不能计作独立复现或因果证明"
)]
setorder(
  associations,
  seed,
  scenario_order,
  gene_name,
  event_type,
  reference_factor
)
fail_if(any(factor_alignment$patient_score_diagnostic_is_gate),
        "患者得分相关被错误用作因子匹配门禁。")
fail_if(any(!factor_alignment$loading_alignment_is_event_blind),
        "因子对齐未全部标记为 event-blind loading alignment。")
fail_if(any(factor_alignment$match_reliable &
              factor_alignment$low_margin_flag),
        "low-margin 因子匹配被错误标记为可靠。")
fail_if(any(factor_alignment$match_reliable &
              factor_alignment$minimum_match_margin < alignment_min_margin),
        "可靠因子匹配绕过了预锁定 margin 门槛。")
reference_alignment_counts <- factor_alignment[
  !is.na(reference_factor), .N, by = model_id
]
fail_if(nrow(reference_alignment_counts) != nrow(model_plan) ||
          any(reference_alignment_counts$N != ncol(reference_factors)),
        "每个模型必须恰好记录 Factor1–Factor8 八条参考对齐。")
association_counts <- associations[, .(
  eligible_relation_count = .N,
  original_network_edge_count = sum(!is.na(original_edge_id))
), by = model_id]
fail_if(nrow(association_counts) != nrow(model_plan) ||
          any(association_counts$eligible_relation_count != 200L) ||
          any(association_counts$original_network_edge_count != 41L),
        "每个模型必须重检 200 条 eligible 关系并标识 41 条正式原边。")
fail_if(any(associations$p_q_scope != "TCGA_internal_descriptive_only") ||
          any(associations$p_q_countable_as_independent_validation),
        "审计 p/q 被错误越界为独立验证证据。")

message("[5/8] 汇总跨 seed 因子匹配、方向与门禁保留率")
stability_source <- associations[
  original_network_edge == TRUE & counts_toward_retention == TRUE
]
edge_stability <- stability_source[, .(
  planned_seed_count = length(seeds),
  informative_model_count = .N,
  completed_converged_seed_count = uniqueN(seed),
  reliable_factor_match_count = sum(factor_match_reliable, na.rm = TRUE),
  reliable_factor_match_rate = mean(factor_match_reliable, na.rm = TRUE),
  direction_retained_count = sum(direction_retained %in% TRUE, na.rm = TRUE),
  direction_retention_rate = mean(direction_retained %in% TRUE, na.rm = TRUE),
  magnitude_retained_count = sum(magnitude_retained_50pct %in% TRUE, na.rm = TRUE),
  magnitude_retention_rate = mean(magnitude_retained_50pct %in% TRUE, na.rm = TRUE),
  gate_retained_count = sum(gate_retained %in% TRUE, na.rm = TRUE),
  gate_retention_rate = mean(gate_retained %in% TRUE, na.rm = TRUE),
  low_margin_count = sum(factor_match_low_margin_flag %in% TRUE),
  median_aligned_effect = median_finite(aligned_effect),
  median_q_value = median_finite(q_value),
  minimum_factor_match_abs_loading_correlation =
    min_finite(factor_match_abs_loading_correlation),
  median_factor_match_abs_loading_correlation =
    median_finite(factor_match_abs_loading_correlation),
  scenario_no_op = FALSE,
  no_op_reason = ""
), by = .(
  gene_id,
  original_edge_id,
  gene_name,
  narrative_role,
  event_type,
  reference_factor,
  original_effect,
  original_q_value,
  original_association_status,
  original_level_factor_soft_flag,
  scenario_id,
  scenario_class,
  scenario_target_gene_id,
  scenario_target_gene_name
)]

no_op_source <- associations[
  original_network_edge == TRUE &
    scenario_class == "gene_event_feature_drop" &
    scenario_no_op == TRUE &
    gene_id == scenario_target_gene_id
]
no_op_stability <- no_op_source[, .(
  planned_seed_count = length(seeds),
  informative_model_count = 0L,
  completed_converged_seed_count = 0L,
  reliable_factor_match_count = NA_integer_,
  reliable_factor_match_rate = NA_real_,
  direction_retained_count = NA_integer_,
  direction_retention_rate = NA_real_,
  magnitude_retained_count = NA_integer_,
  magnitude_retention_rate = NA_real_,
  gate_retained_count = NA_integer_,
  gate_retention_rate = NA_real_,
  low_margin_count = NA_integer_,
  median_aligned_effect = NA_real_,
  median_q_value = NA_real_,
  minimum_factor_match_abs_loading_correlation = NA_real_,
  median_factor_match_abs_loading_correlation = NA_real_,
  scenario_no_op = TRUE,
  no_op_reason = "candidate Mutation and CNV features were both absent from the reference MOFA input; baseline reuse is not sensitivity evidence"
), by = .(
  gene_id,
  original_edge_id,
  gene_name,
  narrative_role,
  event_type,
  reference_factor,
  original_effect,
  original_q_value,
  original_association_status,
  original_level_factor_soft_flag,
  scenario_id,
  scenario_class,
  scenario_target_gene_id,
  scenario_target_gene_name
)]
edge_stability <- rbindlist(
  list(edge_stability, no_op_stability),
  use.names = TRUE,
  fill = TRUE
)
edge_stability[, retention_decision := fcase(
  scenario_no_op, "not_applicable_no_op",
  informative_model_count < planned_seed_count,
    "incomplete_or_nonconverged_seed_models",
  original_level_factor_soft_flag, "level_factor_qc_not_primary",
  reliable_factor_match_rate < 2 / 3, "factor_alignment_unstable",
  direction_retention_rate >= 2 / 3 & gate_retention_rate >= 2 / 3 &
    low_margin_count > 0L,
    "retained_2of3_with_low_margin_seed_recorded",
  direction_retention_rate >= 2 / 3 & gate_retention_rate >= 2 / 3,
    "retained_across_seeds_no_low_margin",
  direction_retention_rate >= 2 / 3 & magnitude_retention_rate >= 2 / 3,
    "direction_and_magnitude_retained_statistical_attenuation",
  direction_retention_rate >= 2 / 3,
    "direction_retained_effect_attenuation",
  default = "not_retained_across_seeds"
)]
edge_stability[, conclusion_ceiling := fifelse(
  scenario_no_op,
  "no-op 只记录边界，不计作泄漏敏感性支持",
  "来源内泄漏控制稳定性；不是独立验证"
)]
setorder(
  edge_stability,
  narrative_role,
  gene_name,
  event_type,
  reference_factor,
  scenario_id
)

edge_keys <- frozen_associations[!is.na(original_edge_id), .(
  original_edge_id,
  gene_id,
  gene_name,
  event_type,
  reference_factor = target,
  original_effect = effect,
  original_q_value = q_value,
  original_association_status = association_status,
  original_level_factor_soft_flag = level_factor_soft_flag
)]
edge_keys <- merge(
  edge_keys,
  main_roles[, .(gene_id, narrative_role)],
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)
edge_keys[is.na(narrative_role), narrative_role := "other_strong_driver"]
fail_if(nrow(edge_keys) != 41L ||
          uniqueN(edge_keys$original_edge_id) != 41L ||
          anyDuplicated(edge_keys[, .(gene_id, event_type, reference_factor)]) > 0L,
        "candidate summary 的冻结原边键必须为 41 条唯一正式 MOFA 边。")

scenario_metrics <- function(
    key_gene_id,
    key_event_type,
    key_factor,
    key_scenario_id) {
  row <- edge_stability[
    gene_id == key_gene_id & event_type == key_event_type &
      reference_factor == key_factor & scenario_id == key_scenario_id
  ]
  if (!nrow(row)) {
    return(list(
      match_rate = NA_real_,
      direction_rate = NA_real_,
      magnitude_rate = NA_real_,
      gate_rate = NA_real_,
      completed_converged_seed_count = NA_integer_,
      low_margin_count = NA_integer_,
      decision = "not_scheduled",
      no_op = NA
    ))
  }
  list(
    match_rate = row$reliable_factor_match_rate[[1L]],
    direction_rate = row$direction_retention_rate[[1L]],
    magnitude_rate = row$magnitude_retention_rate[[1L]],
    gate_rate = row$gate_retention_rate[[1L]],
    completed_converged_seed_count =
      row$completed_converged_seed_count[[1L]],
    low_margin_count = row$low_margin_count[[1L]],
    decision = row$retention_decision[[1L]],
    no_op = row$scenario_no_op[[1L]]
  )
}

critical_scenario_ids <- c("baseline", "drop_both_event_views")
completed_converged_seed_count <- sum(vapply(seeds, function(current_seed) {
  rows <- model_plan[
    seed == current_seed & scenario_id %chin% critical_scenario_ids
  ]
  nrow(rows) == length(critical_scenario_ids) &&
    all(rows$model_completed) && all(rows$model_converged)
}, FUN.VALUE = logical(1)))

candidate_rows <- vector("list", nrow(edge_keys))
for (i in seq_len(nrow(edge_keys))) {
  key <- edge_keys[i]
  gene_id <- key$gene_id[[1L]]
  gene_name <- key$gene_name[[1L]]
  event_type <- key$event_type[[1L]]
  factor <- key$reference_factor[[1L]]
  baseline <- scenario_metrics(gene_id, event_type, factor, "baseline")
  drop_mutation <- scenario_metrics(gene_id, event_type, factor, "drop_mutation_view")
  drop_cnv <- scenario_metrics(gene_id, event_type, factor, "drop_cnv_view")
  drop_both <- scenario_metrics(gene_id, event_type, factor, "drop_both_event_views")
  gene_drop_id <- paste0("drop_gene_event_features_", gene_name)
  gene_drop <- scenario_metrics(gene_id, event_type, factor, gene_drop_id)

  original_supported <- identical(
    key$original_association_status[[1L]], "within_tcga_supported"
  )
  non_level_factor <- !isTRUE(key$original_level_factor_soft_flag[[1L]]) &&
    factor != "Factor4"
  gene_drop_scheduled <- gene_name %chin% main_roles$gene_name
  scheduled_gene_drop_no_op <- gene_drop_scheduled && isTRUE(gene_drop$no_op)
  # baseline 与 drop-both 是正式总门禁，若其一 no-op 则不可评估。
  # 候选自身事件特征本来就未进入参考模型时，gene-drop 属于“无对象可删”，
  # 只记录为 not applicable，既不冒充敏感性支持，也不追加否决。
  gate_required_scenario_no_op <- isTRUE(baseline$no_op) ||
    isTRUE(drop_both$no_op)
  all_critical_seeds_completed <-
    completed_converged_seed_count == length(seeds)
  leakage_gate_evaluable <- original_supported && non_level_factor &&
    isTRUE(reference_training_diagnostics$model_converged) &&
    all_critical_seeds_completed && !gate_required_scenario_no_op
  leakage_gate_pass <- leakage_gate_evaluable &&
    is.finite(baseline$match_rate) && baseline$match_rate >= 2 / 3 &&
    is.finite(drop_both$match_rate) && drop_both$match_rate >= 2 / 3 &&
    is.finite(drop_both$direction_rate) && drop_both$direction_rate >= 2 / 3 &&
    is.finite(drop_both$magnitude_rate) && drop_both$magnitude_rate >= 2 / 3 &&
    is.finite(drop_both$gate_rate) && drop_both$gate_rate >= 2 / 3

  failure_reasons <- character()
  if (!original_supported) failure_reasons <- c(
    failure_reasons, "original_edge_not_within_tcga_supported_no_upgrade_allowed"
  )
  if (isTRUE(key$original_level_factor_soft_flag[[1L]]) || factor == "Factor4") {
    failure_reasons <- c(failure_reasons, "level_factor_or_Factor4_not_primary")
  }
  if (!isTRUE(reference_training_diagnostics$model_converged)) {
    failure_reasons <- c(failure_reasons, "reference_model_not_converged")
  }
  if (!all_critical_seeds_completed) failure_reasons <- c(
    failure_reasons, "not_all_baseline_and_drop_both_seeds_completed_converged"
  )
  if (gate_required_scenario_no_op) failure_reasons <- c(
    failure_reasons, "required_baseline_or_drop_both_scenario_no_op"
  )
  if (!is.finite(baseline$match_rate) || baseline$match_rate < 2 / 3) {
    failure_reasons <- c(failure_reasons, "baseline_reliable_match_rate_below_2of3")
  }
  if (!is.finite(drop_both$match_rate) || drop_both$match_rate < 2 / 3) {
    failure_reasons <- c(failure_reasons, "drop_both_reliable_match_rate_below_2of3")
  }
  if (!is.finite(drop_both$direction_rate) || drop_both$direction_rate < 2 / 3) {
    failure_reasons <- c(failure_reasons, "drop_both_direction_retention_below_2of3")
  }
  if (!is.finite(drop_both$magnitude_rate) || drop_both$magnitude_rate < 2 / 3) {
    failure_reasons <- c(failure_reasons, "drop_both_magnitude_retention_below_2of3")
  }
  if (!is.finite(drop_both$gate_rate) || drop_both$gate_rate < 2 / 3) {
    failure_reasons <- c(failure_reasons, "drop_both_supported_gate_retention_below_2of3")
  }
  gate_failure_reason <- if (length(failure_reasons)) {
    paste(unique(failure_reasons), collapse = ";")
  } else {
    ""
  }
  overall <- if (leakage_gate_pass) {
    "original_supported_edge_retained_after_prelocked_leakage_gate"
  } else if (!original_supported) {
    "original_non_supported_edge_cannot_upgrade"
  } else if (!leakage_gate_evaluable) {
    "supported_edge_leakage_gate_not_evaluable"
  } else {
    "supported_edge_failed_prelocked_leakage_gate"
  }

  candidate_rows[[i]] <- data.table(
    original_edge_id = key$original_edge_id[[1L]],
    gene_id = gene_id,
    gene_name = gene_name,
    narrative_role = key$narrative_role[[1L]],
    event_type = event_type,
    reference_factor = factor,
    original_effect = key$original_effect[[1L]],
    original_q_value = key$original_q_value[[1L]],
    original_association_status = key$original_association_status[[1L]],
    original_level_factor_soft_flag = key$original_level_factor_soft_flag[[1L]],
    planned_seed_count = length(seeds),
    completed_converged_seed_count = completed_converged_seed_count,
    baseline_reliable_match_rate = baseline$match_rate,
    baseline_direction_retention_rate = baseline$direction_rate,
    baseline_gate_retention_rate = baseline$gate_rate,
    drop_mutation_match_rate = drop_mutation$match_rate,
    drop_mutation_direction_retention_rate = drop_mutation$direction_rate,
    drop_mutation_gate_retention_rate = drop_mutation$gate_rate,
    drop_cnv_match_rate = drop_cnv$match_rate,
    drop_cnv_direction_retention_rate = drop_cnv$direction_rate,
    drop_cnv_gate_retention_rate = drop_cnv$gate_rate,
    drop_both_reliable_match_rate = drop_both$match_rate,
    drop_both_direction_retention_rate = drop_both$direction_rate,
    drop_both_magnitude_retention_rate = drop_both$magnitude_rate,
    drop_both_gate_retention_rate = drop_both$gate_rate,
    drop_both_low_margin_count = drop_both$low_margin_count,
    drop_both_retention_decision = drop_both$decision,
    gene_feature_drop_scenario = gene_drop_id,
    gene_feature_drop_prespecified = gene_drop_scheduled,
    gene_feature_drop_no_op = gene_drop$no_op,
    scheduled_gene_feature_drop_no_op = scheduled_gene_drop_no_op,
    gate_required_scenario_no_op = gate_required_scenario_no_op,
    gene_feature_drop_match_rate = gene_drop$match_rate,
    gene_feature_drop_direction_retention_rate = gene_drop$direction_rate,
    gene_feature_drop_gate_retention_rate = gene_drop$gate_rate,
    leakage_gate_evaluable = leakage_gate_evaluable,
    leakage_gate_pass = leakage_gate_pass,
    countable_for_T3_T4 = leakage_gate_pass,
    gate_failure_reason = gate_failure_reason,
    post_audit_maximum_status = fifelse(
      leakage_gate_pass,
      key$original_association_status[[1L]],
      fifelse(
        original_supported,
        "supported_not_retained_for_T3_T4",
        key$original_association_status[[1L]]
      )
    ),
    leakage_control_decision = overall,
    exact_independent_validation = FALSE,
    p_q_scope = "TCGA_internal_descriptive_only",
    conclusion_ceiling = if (key$narrative_role[[1L]] == "boundary_context_reverse") {
      "情境/反向边界不因来源内泄漏控制而升级；仍需独立配对队列和功能验证"
    } else {
      "TCGA 来源内泄漏控制候选；不等于独立复现、因果驱动或治疗靶点"
    },
    required_next_validation = paste(
      "独立 ESCC 队列预锁定事件→状态投影；细胞组成校准；必要的扰动实验"
    )
  )
}
candidate_summary <- rbindlist(candidate_rows)
required_candidate_summary_columns <- c(
  "original_edge_id",
  "gene_id",
  "event_type",
  "reference_factor",
  "original_association_status",
  "original_level_factor_soft_flag",
  "planned_seed_count",
  "completed_converged_seed_count",
  "baseline_reliable_match_rate",
  "drop_both_reliable_match_rate",
  "drop_both_direction_retention_rate",
  "drop_both_magnitude_retention_rate",
  "drop_both_gate_retention_rate",
  "drop_both_low_margin_count",
  "drop_both_retention_decision",
  "leakage_gate_evaluable",
  "leakage_gate_pass",
  "countable_for_T3_T4",
  "gate_failure_reason"
)
require_columns(
  candidate_summary,
  required_candidate_summary_columns,
  "tcga_escc_driver_state_leakage_candidate_summary.tsv"
)
setcolorder(
  candidate_summary,
  c(
    required_candidate_summary_columns,
    setdiff(names(candidate_summary), required_candidate_summary_columns)
  )
)
setorder(
  candidate_summary,
  narrative_role,
  gene_name,
  event_type,
  reference_factor
)
fail_if(nrow(candidate_summary) != 41L ||
          uniqueN(candidate_summary$original_edge_id) != 41L ||
          !setequal(
            candidate_summary$original_edge_id,
            original_mofa_network$edge_id
          ),
        "candidate summary 必须每条正式 MOFA 原边恰好一行。")
fail_if(any(candidate_summary$countable_for_T3_T4 &
              (
                candidate_summary$original_association_status !=
                  "within_tcga_supported" |
                candidate_summary$original_level_factor_soft_flag |
                candidate_summary$reference_factor == "Factor4" |
                !candidate_summary$leakage_gate_evaluable |
                !candidate_summary$leakage_gate_pass |
                candidate_summary$completed_converged_seed_count !=
                  candidate_summary$planned_seed_count |
                candidate_summary$gate_required_scenario_no_op
              )),
        "countable_for_T3_T4 越过了原边层级、收敛或 no-op 门禁。")
fail_if(any(
  candidate_summary$original_association_status != "within_tcga_supported" &
    candidate_summary$countable_for_T3_T4
), "conditional/exploratory 原边被错误升级为 T3/T4 可计数边。")
fail_if(any(candidate_summary$leakage_gate_pass !=
              candidate_summary$countable_for_T3_T4),
        "leakage_gate_pass 与 countable_for_T3_T4 不一致。")
fail_if(any(candidate_summary$leakage_gate_pass &
              nzchar(candidate_summary$gate_failure_reason)),
        "门禁通过行仍含 gate_failure_reason。")
fail_if(any(!candidate_summary$leakage_gate_pass &
              !nzchar(candidate_summary$gate_failure_reason)),
        "门禁未通过行缺少 gate_failure_reason。")

message("[6/8] 冻结 ECMS 外部投影边界")
external_projection_status <- data.table(
  benchmark = "ECMS1-MET/ECMS2-CLS/ECMS3-IM/ECMS4-MES",
  benchmark_role = "external_consensus_subtype_reference_only",
  locked_public_projection_weights_available = ecms_projection_available,
  projection_artifact_manifest_verified = ecms_projection_available,
  projection_status = if (ecms_projection_available) {
    "available_downstream_verified_not_used_in_leakage_gate"
  } else {
    "pending_no_verified_projection_artifact_in_project_authority_surface"
  },
  model_commit = ecms_model_commit,
  repository_anchor_patient_count = ecms_anchor_patient_count,
  extension_calibration_pass = ecms_extension_calibration_pass,
  primary_patient_count = ecms_primary_patient_count,
  normalization_scope = if (ecms_projection_available) {
    "external_cohort_batch_scaling_not_single_sample"
  } else {
    "pending"
  },
  single_sample_valid = FALSE,
  used_in_mofa_training_or_leakage_gate = FALSE,
  verified_project_surfaces = if (ecms_projection_available) {
    paste(
      "results/tcga_escc_ecms_patient_probabilities.tsv;",
      "results/tcga_escc_ecms_projection_calibration.tsv;",
      "results/tcga_escc_ecms_projection_qa.tsv;",
      "results/tcga_escc_ecms_projection_artifact_manifest.tsv"
    )
  } else {
    "PROJECT_INDEX.md;results;data/datasets.tsv"
  },
  pseudo_labels_generated = FALSE,
  pseudo_label_prohibition = "cluster number or pathway resemblance cannot be converted into ECMS labels",
  required_resource = if (ecms_projection_available) {
    "satisfied by locked author model plus verified cohort-level preprocessing and calibration"
  } else {
    paste(
      "locked author classifier, feature preprocessing, gene identifiers,",
      "cohort-level calibration and artifact manifest"
    )
  },
  next_action = if (ecms_projection_available) {
    "use verified ECMS labels only as an external coordinate; never feed them back into MOFA training"
  } else {
    "run verified ECMS ingest/projection before external-coordinate interpretation"
  },
  conclusion_ceiling = if (ecms_projection_available) {
    "verified cohort-level external coordinate; not independent biological validation or clinical single-sample classifier"
  } else {
    "pending external benchmark; no ECMS assignment in this leakage analysis"
  }
)

message("[7/8] 生成摘要、QA、manifest 并执行结构回读")
stage_tsv(model_plan, "tcga_escc_driver_state_leakage_model_plan.tsv")
stage_tsv(factor_alignment, "tcga_escc_driver_state_leakage_factor_alignment.tsv")
stage_tsv(associations, "tcga_escc_driver_state_leakage_associations.tsv")
stage_tsv(edge_stability, "tcga_escc_driver_state_leakage_edge_stability.tsv")
stage_tsv(candidate_summary, "tcga_escc_driver_state_leakage_candidate_summary.tsv")
stage_tsv(external_projection_status, "tcga_escc_driver_state_external_projection_status.tsv")

format_metric <- function(value) {
  ifelse(is.finite(value), sprintf("%.2f", value), "NA")
}
main_candidate_lines <- candidate_summary[
  narrative_role != "other_strong_driver",
  paste0(
    "- `", gene_name, ":", event_type, "` → `", reference_factor,
    "`：", leakage_control_decision,
    "；baseline match=", format_metric(baseline_reliable_match_rate),
    "；drop-both direction=", format_metric(drop_both_direction_retention_rate),
    "；drop-both gate=", format_metric(drop_both_gate_retention_rate),
    "；gene-drop no-op=", gene_feature_drop_no_op, "。"
  )
]
if (!length(main_candidate_lines)) {
  main_candidate_lines <- "- 主叙事候选没有进入原始网络边；仅保留全关联表。"
}

summary_lines <- c(
  "# TCGA-ESCC 驱动事件—状态表示泄漏审计摘要",
  "",
  "## 分析范围",
  "",
  paste0(
    "- 使用 script 14 冻结的五视图：",
    paste(names(expected_feature_counts), expected_feature_counts, sep = "=", collapse = "；"),
    "；患者为 94 位 five-layer core。"
  ),
  paste0(
    "- 默认 seed：", paste(seeds, collapse = "、"),
    "；模型计划 ", expected_planned_models, " 行；实际拟合/恢复新场景 ",
    expected_fitted_models, " 个；no-op baseline 复用 ", expected_no_op_reuses, " 个。"
  ),
  paste0(
    "- 因子匹配只使用共享 RNA/miRNA/HM450 全载荷：每个视图内",
    "按因子中心化并单位 L2 标准化后拼接，以载荷 Pearson |r| 做",
    " Hungarian 一对一匹配；可靠门槛为 |r|≥",
    alignment_min_abs_correlation, " 且最小 margin≥", alignment_min_margin,
    "。高相关 override（|r|≥", alignment_high_confidence_override,
    "）仍必须满足 margin≥", alignment_min_margin, "。"
  ),
  paste0(
    "- 可靠匹配后按载荷相关符号翻转新因子；患者得分 ",
    alignment_method,
    " 相关仅是诊断字段，不参与分配、符号或门禁。"
  ),
  "- 逐原边门禁仅允许原本即为 within_tcga_supported 且非 Factor4/非 level-factor 的边；baseline 和 drop-both 全部 3 个预锁定 seed 必须完成并收敛，两者可靠对齐率均需≥2/3，drop-both 方向、≥50% 效应幅度和原 supported 级保留率均需≥2/3。low-margin 单次匹配不计作可靠匹配并完整留痕，但不在满足 2/3 复现门槛后追加一票否决。",
  "- p/q 按每个模型、每个 event type 的原始完整检验家族计算，只是 TCGA 同来源内部描述；不计作独立验证或新增支持层。",
  "",
  "## 主叙事候选",
  "",
  main_candidate_lines,
  "",
  "## no-op 与解释边界",
  "",
  paste0(
    "- 主叙事 Mutation/CNV 特征状态：",
    paste0(
      main_roles$gene_name,
      "[Mutation=", main_roles$mutation_feature_in_reference_model,
      ",CNV=", main_roles$cnv_feature_in_reference_model, "]",
      collapse = "；"
    ), "。"
  ),
  "- 若两种事件特征均未进入参考模型，gene-specific 场景复用 seed-matched baseline 并标记 no-op；这些行不进入保留率、不冒充支持，也不因本来无直接事件特征可删而追加否决。",
  paste0(
    "- 冻结 manifest 有 ", nrow(cross_view_duplicate_ids),
    " 个跨视图重复 feature_id；删除和 no-op 判定在 `create_mofa` 前按 view-qualified 原始 ID 完成，",
    "不依赖 MOFA2 自动添加的跨视图后缀。"
  ),
  "- 删除整个 Mutation/CNV view 控制的是事件层表示泄漏；即使边保留，也仍与同一 TCGA 患者共享数据，不能称独立验证。",
  "",
  "## ECMS 外部投影",
  "",
  paste0("- 状态：`", external_projection_status$projection_status, "`。"),
  if (ecms_projection_available) {
    paste0(
      "- 已验证锁定作者模型 commit `", ecms_model_commit,
      "`；仓库矩阵 anchor=", ecms_anchor_patient_count,
      "，扩展校准=", ecms_extension_calibration_pass,
      "，primary n=", ecms_primary_patient_count,
      "。该坐标未参与 MOFA 训练、因子匹配或泄漏门禁。"
    )
  } else {
    "- 未获得完整 verified 投影 artifact 前，不生成 ECMS1–4 伪标签，也不以聚类数相同冒充复现。"
  },
  "",
  "## 结论上限",
  "",
  "本审计用于判定当前 TCGA 内部 driver-state 边对随机种子和事件表示删除的敏感性。它能识别明显表示泄漏、符号不稳定和因子不可对齐，但不能替代独立 ESCC 队列、外部 ECMS 投影、细胞来源校准或功能实验。"
)
summary_filename <- "tcga_escc_driver_state_leakage_summary.md"
writeLines(summary_lines, file.path(stage_dir, summary_filename), useBytes = TRUE)

formal_filenames <- c(
  "tcga_escc_driver_state_leakage_model_plan.tsv",
  "tcga_escc_driver_state_leakage_factor_alignment.tsv",
  "tcga_escc_driver_state_leakage_associations.tsv",
  "tcga_escc_driver_state_leakage_edge_stability.tsv",
  "tcga_escc_driver_state_leakage_candidate_summary.tsv",
  "tcga_escc_driver_state_external_projection_status.tsv",
  summary_filename
)
formal_paths <- file.path(stage_dir, formal_filenames)
fail_if(any(!file_exists(formal_paths)), "泄漏审计正式 artifact 未完整生成。")

artifact_manifest <- data.table(
  artifact = formal_filenames,
  relative_path = file.path("results", formal_filenames),
  file_size_bytes = as.numeric(file_info(formal_paths)$size),
  sha256 = vapply(
    formal_paths,
    digest,
    FUN.VALUE = character(1),
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  ),
  generated_date = as.character(Sys.Date()),
  generation_script = "scripts/21_audit_tcga_escc_driver_state_leakage.R",
  generation_script_sha256 = execution_script_sha256,
  input_signature = input_signature,
  independence_group = "TCGA_ESCC_DR45_driver_state_shared_patients",
  status = "verified"
)
manifest_filename <- "tcga_escc_driver_state_leakage_artifact_manifest.tsv"
stage_tsv(artifact_manifest, manifest_filename)

for (filename in formal_filenames) {
  if (grepl("[.]tsv$", filename)) {
    reread <- read_tsv(file.path(stage_dir, filename))
    fail_if(!nrow(reread), paste("正式 TSV 为空：", filename))
  }
}
fail_if(nrow(model_plan) != expected_planned_models,
        "model plan 行数与冻结计划不一致。")
fail_if(model_plan[scenario_no_op == FALSE, .N] != expected_fitted_models,
        "实际拟合计划数异常。")
fail_if(any(model_plan$run_status == "planned_not_run"), "仍有未执行模型。")
fail_if(any(factor_alignment$matching_status == "reliable_hungarian_match" &
              !factor_alignment$match_reliable), "可靠匹配状态与布尔标志冲突。")
fail_if(any(associations$counts_toward_retention & associations$scenario_no_op),
        "no-op 被错误计入保留率。")
fail_if(any(candidate_summary$exact_independent_validation),
        "来源内泄漏审计被错误标记为独立验证。")
fail_if(any(external_projection_status$pseudo_labels_generated),
        "ECMS 伪标签禁止标志失败。")
fail_if(any(external_projection_status$used_in_mofa_training_or_leakage_gate),
        "ECMS 外部坐标被错误反馈进 MOFA 或泄漏门禁。")
fail_if(
  ecms_projection_available &&
    (!external_projection_status$projection_artifact_manifest_verified ||
       external_projection_status$projection_status !=
         "available_downstream_verified_not_used_in_leakage_gate"),
  "已验证 ECMS artifact 未正确写入外部坐标状态。"
)

qa_lines <- c(
  "# TCGA-ESCC driver-state 表示泄漏审计 QA",
  "",
  "本文件是运行时历史证据，不是项目当前状态源。",
  "",
  paste0("- 运行日期：", Sys.Date(), "。"),
  paste0("- 输入 manifest 大小与 SHA256：全部通过；input signature：", input_signature, "。"),
  paste0("- 五视图数值 SHA256：",
         paste(names(base_view_numeric_hashes), base_view_numeric_hashes,
               sep = "=", collapse = "；"), "；聚合 hash：",
         base_views_numeric_hash, "。"),
  paste0("- 执行脚本 SHA256：", execution_script_sha256, "。"),
  paste0(
    "- MOFA 模型训练契约脚本 SHA256：",
    model_training_contract_script_sha256,
    "；当前脚本仅修正状态回写并复用逐项验证通过的模型检查点。"
  ),
  paste0("- 冻结患者：", length(patients), "；冻结特征：",
         paste(names(expected_feature_counts), expected_feature_counts,
               sep = "=", collapse = "；"), "。"),
  paste0("- 跨视图重复 feature_id：", nrow(cross_view_duplicate_ids),
         "；按 view-qualified manifest 删除，未用 MOFA2 后缀判断 no-op。"),
  paste0("- 模型计划：", expected_planned_models, "；应拟合/恢复：",
         expected_fitted_models, "；no-op baseline 复用：", expected_no_op_reuses, "。"),
  paste0("- 实际 run status：", collapse_unique(paste0(
    model_plan[, .N, by = run_status]$run_status, "=",
    model_plan[, .N, by = run_status]$N
  )), "。"),
  paste0("- 完成且收敛的非 no-op 模型：",
         model_plan[scenario_no_op == FALSE & model_completed & model_converged, .N],
         "/", expected_fitted_models, "；参考模型 iterations=",
         reference_training_diagnostics$training_iteration_count,
         "，final ELBO=", reference_training_diagnostics$final_elbo, "。"),
  paste0("- 因子对齐行：", nrow(factor_alignment), "；可靠匹配：",
         sum(factor_alignment$match_reliable, na.rm = TRUE), "。"),
  paste0("- 重检关联：", nrow(associations), "；edge stability：",
         nrow(edge_stability), "；candidate summary：", nrow(candidate_summary), "。"),
  "- Hungarian 匹配与符号只使用 event-blind RNA/miRNA/HM450 标准化全载荷；患者得分相关仅作诊断。low-margin 不计作可靠匹配，但允许由其余至少 2/3 预锁定 seed 的一致结果通过弹性门禁。",
  "- no-op 行不计入方向、匹配或门禁保留率。",
  "- p/q 使用固定原始检验家族，仅为 TCGA 内部描述；跨 seed 保留率不冒充患者 bootstrap 或外部复现。",
  paste0(
    "- ECMS 外部坐标状态：",
    external_projection_status$projection_status,
    "；primary n=",
    ifelse(is.na(ecms_primary_patient_count), "NA", ecms_primary_patient_count),
    "；未进入 MOFA/泄漏门禁；伪标签生成数为 0。"
  ),
  paste0("- checkpoint：", cache_dir, "；仅可再生成中间层，不是当前状态源。"),
  paste0("- ", r_version, "；MOFA2 ", mofa2_version,
         "；mofapy2 ", mofapy2_version,
         "；clue ", packageVersion("clue"), "；data.table ",
         packageVersion("data.table"), "。"),
  "- ResearchDataHub 原始资料、项目上游结果、投稿包和 ECMS 标签均未修改或伪造。"
)
qa_filename <- "tcga_escc_driver_state_leakage_qa_20260711.md"
writeLines(qa_lines, file.path(stage_dir, qa_filename), useBytes = TRUE)

message("[8/8] 原子发布正式结果并逐项回读 SHA256")
for (filename in formal_filenames) {
  atomic_publish_file(
    file.path(stage_dir, filename),
    file.path(results_dir, filename)
  )
}
# 先逐项验证正式目标，再把 manifest 作为最后一个正式 artifact 原子发布。
# 这样即使目标验证失败，results/ 也不会出现宣称新版本已 verified 的 manifest。
for (i in seq_len(nrow(artifact_manifest))) {
  published <- file.path(project_root, artifact_manifest$relative_path[[i]])
  fail_if(!file_exists(published), paste("发布 artifact 缺失：", published))
  fail_if(as.numeric(file_info(published)$size) != artifact_manifest$file_size_bytes[[i]],
          paste("发布 artifact 大小不一致：", published))
  fail_if(
    digest(published, algo = "sha256", file = TRUE, serialize = FALSE) !=
      artifact_manifest$sha256[[i]],
    paste("发布 artifact SHA256 不一致：", published)
  )
}
atomic_publish_file(
  file.path(stage_dir, manifest_filename),
  file.path(results_dir, manifest_filename)
)
published_manifest <- read_tsv(file.path(results_dir, manifest_filename))
fail_if(
  nrow(published_manifest) != nrow(artifact_manifest) ||
    !setequal(published_manifest$relative_path, artifact_manifest$relative_path) ||
    any(published_manifest$sha256[
      match(artifact_manifest$relative_path, published_manifest$relative_path)
    ] != artifact_manifest$sha256),
  "results/ 发布 manifest 回读与 stage 冻结版不一致。"
)
atomic_publish_file(
  file.path(stage_dir, qa_filename),
  file.path(work_checks_dir, qa_filename)
)

if (dir_exists(stage_dir)) dir_delete(stage_dir)
message(
  "完成：计划 ", expected_planned_models, " 个模型场景行；实际训练/恢复 ",
  expected_fitted_models, " 个，no-op ", expected_no_op_reuses,
  " 个；关联 ", nrow(associations), " 行。"
)
