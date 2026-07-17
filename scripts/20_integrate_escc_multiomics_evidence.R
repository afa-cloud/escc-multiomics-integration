#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(fs)
})

options(stringsAsFactors = FALSE, scipen = 999)

# ESCC 多组学证据整合层。
#
# 本脚本只整合已经发布到 results/ 的正式结果，不重跑上游统计模型。
# 核心边界：
# 1) 不把不同队列拼成患者级关联；
# 2) 同一 TCGA 患者上的 driver、MOFA 和 PROGENy 关系不计为独立验证；
# 3) 缺失组学层记为不可比较，不记为阴性；
# 4) Cao 小样本方向性结果、情境特异反向和 PR001876/PRJNA766558
#    正交模块可条件保留，但结论上限随来源边界收窄；
# 5) T0-T4 是证据裁决主轴，数值分数只作排序辅助，不能越过硬门禁或结论上限。

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("无法唯一定位当前脚本。")
script_path <- normalizePath(
  sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE
)
executed_script_sha256 <- digest(
  script_path,
  algo = "sha256",
  file = TRUE,
  serialize = FALSE
)
project_root <- dirname(dirname(script_path))
results_dir <- file.path(project_root, "results")
work_intermediate_dir <- file.path(project_root, "_work", "intermediate")
work_checks_dir <- file.path(project_root, "_work", "checks")

args <- commandArgs(trailingOnly = TRUE)
fields_only <- "--fields-only" %in% args
self_test <- "--self-test" %in% args
unknown_args <- setdiff(args, c("--fields-only", "--self-test"))
if (length(unknown_args)) {
  stop("未知参数：", paste(unknown_args, collapse = ", "), call. = FALSE)
}
if (fields_only && self_test) {
  stop("--fields-only 与 --self-test 不能同时使用。", call. = FALSE)
}

integration_formal_filenames <- c(
  "escc_multiomics_integration_tier_definitions.tsv",
  "escc_multiomics_evidence_ledger.tsv",
  "escc_multiomics_integrated_driver_candidates.tsv",
  "escc_multiomics_integrated_axis_edges.tsv",
  "escc_multiomics_integrated_axis_summary.tsv",
  "escc_multiomics_heterogeneity_axes.tsv",
  "escc_multiomics_module_summaries.tsv",
  "escc_multiomics_integration_summary.md"
)
integration_manifest_filename <-
  "escc_multiomics_integration_artifact_manifest.tsv"

if (fields_only) {
  cat(
    paste(c(integration_formal_filenames, integration_manifest_filename),
          collapse = "\n"),
    "\n",
    sep = ""
  )
  cat(
    "required_upstream\tscript21 7+manifest; ECMS 6+manifest; ",
    "PRJNA766558 12+manifest when any formal family member exists\n",
    "ecms_primary\tofficial_anchor_78 when extension_calibration_pass=FALSE\n",
    "tier_caps\tECMS/PROGENy/PR001876/PRJNA766558 <= T2\n",
    "hard_boundaries\texact artifact family; strict logical fields; ",
    "no cross-cohort patient linkage; no representation leakage\n",
    sep = ""
  )
  quit(save = "no", status = 0L)
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

dir_create(c(work_intermediate_dir, work_checks_dir), recurse = TRUE)

fail_if <- function(condition, message) {
  if (!is.logical(condition) || length(condition) != 1L || is.na(condition)) {
    stop(
      paste0(message, "（门禁条件为 NA 或非单一逻辑值，按失败处理）"),
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

collapse_unique <- function(value, empty = "") {
  value <- unique(as.character(value[!is.na(value) & nzchar(as.character(value))]))
  if (!length(value)) empty else paste(sort(value), collapse = ";")
}

min_finite <- function(value) {
  value <- value[is.finite(value)]
  if (!length(value)) NA_real_ else min(value)
}

first_finite <- function(value) {
  value <- value[is.finite(value)]
  if (!length(value)) NA_real_ else value[[1L]]
}

as_logical_strict <- function(value, label = "logical field", allow_na = FALSE) {
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

# 只用于脚本内部生成、允许结构性 NA 的字段；外部门禁字段必须
# 显式调用 as_logical_strict(..., allow_na = FALSE)。
as_logical_safe <- function(value) {
  as_logical_strict(value, "internal logical field", allow_na = TRUE)
}

tier_rank <- c(T0 = 0L, T1 = 1L, T2 = 2L, T3 = 3L, T4 = 4L)
tier_score <- c(T0 = 0, T1 = 25, T2 = 55, T3 = 78, T4 = 90)

assign_state_edge_tier <- function(target_type, level_flag, leakage_countable,
                                   association_status) {
  fcase(
    target_type == "MOFA_factor" & level_flag, "T0",
    target_type == "MOFA_factor" & leakage_countable, "T3",
    target_type == "MOFA_factor" & association_status %chin% c(
      "within_tcga_supported", "within_tcga_conditional"
    ), "T2",
    target_type == "PROGENy_pathway" & association_status %chin% c(
      "within_tcga_supported", "within_tcga_conditional"
    ), "T2",
    default = "T1"
  )
}

leakage_gate_contract <- function(
    original_status,
    level_flag,
    evaluable,
    completed_seed_count,
    planned_seed_count,
    baseline_match_rate,
    drop_both_match_rate,
    drop_both_direction_rate,
    drop_both_magnitude_rate,
    drop_both_gate_rate) {
  original_status == "within_tcga_supported" &
    !level_flag & evaluable &
    completed_seed_count == planned_seed_count &
    baseline_match_rate >= 2 / 3 &
    drop_both_match_rate >= 2 / 3 &
    drop_both_direction_rate >= 2 / 3 &
    drop_both_magnitude_rate >= 2 / 3 &
    drop_both_gate_rate >= 2 / 3
}

is_forbidden_unpaired_host_source <- function(value) {
  grepl("PR001876|PRJNA766558", as.character(value))
}

artifact_family_exact <- function(expected_names, observed_names) {
  length(observed_names) == length(expected_names) &&
    !anyDuplicated(observed_names) &&
    setequal(observed_names, expected_names)
}

assign_ecms_factor_tier <- function(level_flag, anova_p, anova_q,
                                    kruskal_p, kruskal_q, eta_squared) {
  min_p <- pmin(anova_p, kruskal_p, na.rm = TRUE)
  min_q <- pmin(anova_q, kruskal_q, na.rm = TRUE)
  fcase(
    level_flag, "T0",
    min_q <= 0.10 | (min_p <= 0.05 & eta_squared >= 0.10), "T2",
    default = "T1"
  )
}

assign_ecms_adjusted_tier <- function(level_flag, incremental_p,
                                      incremental_q, partial_r_squared) {
  fcase(
    level_flag, "T0",
    (incremental_q <= 0.10 & partial_r_squared >= 0.05) |
      (incremental_q <= 0.20 & partial_r_squared >= 0.10), "T2",
    default = "T1"
  )
}

assign_heterogeneity_tier <- function(level_flag, leakage_t3_count,
                                      ecms_context_t2,
                                      adjusted_pathway_t2_count,
                                      conditional_driver_count = 0L) {
  fcase(
    level_flag, "T0",
    leakage_t3_count >= 1L & ecms_context_t2 &
      adjusted_pathway_t2_count >= 1L, "T3",
    leakage_t3_count >= 1L | ecms_context_t2 |
      adjusted_pathway_t2_count >= 1L | conditional_driver_count >= 1L, "T2",
    default = "T1"
  )
}

assign_axis_summary_tier <- function(distinct_event_units, leakage_t3_count,
                                     cao_rna_directional,
                                     ecms_contextualized_leakage_count,
                                     conditional_state_count = 0L,
                                     cao_conditional = FALSE,
                                     boundary_reverse = FALSE) {
  fcase(
    boundary_reverse, "T2",
    distinct_event_units >= 2L & leakage_t3_count >= 1L &
      cao_rna_directional & ecms_contextualized_leakage_count >= 1L, "T4",
    cao_conditional, "T2",
    leakage_t3_count >= 1L, "T3",
    conditional_state_count >= 1L, "T2",
    default = "T1"
  )
}

assign_metabolite_tier <- function(candidate_tier) {
  fcase(
    candidate_tier == "run_order_confounded_conditional", "T0",
    candidate_tier %chin% c("robust_rank_and_scale_fdr", "conditional_fdr"),
      "T2",
    default = "T1"
  )
}

assign_microbe_alpha_tier <- function(interpretation_scope, p_value, q_value,
                                      direction_consistency) {
  fcase(
    interpretation_scope == "paired_final_depth_diagnostic_not_diversity", "T0",
    interpretation_scope == "100_repeat_common_depth_rarefaction_sensitivity" &
      direction_consistency >= 0.80 & (q_value <= 0.10 | p_value <= 0.05), "T2",
    default = "T1"
  )
}

assign_microbe_beta_tier <- function(p_value, q_value, betadisper_p) {
  fifelse(
    betadisper_p > 0.05 & (q_value <= 0.10 | p_value <= 0.05),
    "T2", "T1"
  )
}

run_self_tests <- function() {
  # 1) 零 leakage gate：ECMS/PROGENy 再显著也不得生成综合 T3/T4。
  zero_state <- assign_state_edge_tier(
    rep("MOFA_factor", 41L), rep(FALSE, 41L), rep(FALSE, 41L),
    rep("within_tcga_supported", 41L)
  )
  fail_if(any(zero_state == "T3"), "SELFTEST1: zero leakage 仍产生 T3 edge。")
  zero_heterogeneity <- assign_heterogeneity_tier(
    FALSE, 0L, TRUE, 12L, 0L
  )
  zero_axis <- assign_axis_summary_tier(
    2L, 0L, TRUE, 0L, 0L, FALSE, FALSE
  )
  fail_if(zero_heterogeneity %chin% c("T3", "T4") ||
            zero_axis %chin% c("T3", "T4"),
          "SELFTEST1: zero leakage 仍产生 axis/heterogeneity T3/T4。")

  # 2) 3 seed 中 1 个 low-margin 不单票否决；2/3 已由 script21 冻结为 countable。
  low_margin_fixture <- data.table(
    original_status = "within_tcga_supported",
    level_flag = FALSE,
    evaluable = TRUE,
    completed = 3L,
    planned = 3L,
    baseline_match = 2 / 3,
    drop_match = 2 / 3,
    drop_direction = 2 / 3,
    drop_magnitude = 2 / 3,
    drop_gate = 2 / 3,
    low_margin_count = 1L
  )
  low_margin_countable <- leakage_gate_contract(
    low_margin_fixture$original_status,
    low_margin_fixture$level_flag,
    low_margin_fixture$evaluable,
    low_margin_fixture$completed,
    low_margin_fixture$planned,
    low_margin_fixture$baseline_match,
    low_margin_fixture$drop_match,
    low_margin_fixture$drop_direction,
    low_margin_fixture$drop_magnitude,
    low_margin_fixture$drop_gate
  )
  low_margin_state <- assign_state_edge_tier(
    "MOFA_factor", FALSE, low_margin_countable, "within_tcga_supported"
  )
  fail_if(low_margin_fixture$low_margin_count != 1L ||
            low_margin_state != "T3",
          "SELFTEST2: 2/3 稳定且 1 seed low-margin 被额外否决。")

  # 3) ECMS 校准 FALSE 时只能选官方 78。
  scope_fixture <- data.table(
    analysis_scope = c("official_anchor_78", "hybrid_94_extension"),
    is_primary_scope = c(TRUE, FALSE),
    n_patients = c(78L, 94L),
    extension_calibration_pass = c(FALSE, FALSE)
  )
  selected_scope <- scope_fixture[
    is_primary_scope & analysis_scope == "official_anchor_78"
  ]
  fail_if(nrow(selected_scope) != 1L || selected_scope$n_patients != 78L,
          "SELFTEST3: ECMS calibration FALSE 未唯一选中 official78。")

  # 4) 共享 RNA 表示即使极显著也最高 T2。
  ecms_tier <- assign_ecms_factor_tier(
    FALSE, 1e-30, 1e-25, 1e-20, 1e-18, 0.90
  )
  adjusted_tier <- assign_ecms_adjusted_tier(FALSE, 1e-30, 1e-25, 0.90)
  factor4_ecms_tier <- assign_ecms_factor_tier(
    TRUE, 1e-30, 1e-25, 1e-20, 1e-18, 0.90
  )
  factor4_adjusted_tier <- assign_ecms_adjusted_tier(
    TRUE, 1e-30, 1e-25, 0.90
  )
  fail_if(any(c(ecms_tier, adjusted_tier) != "T2") ||
            any(c(factor4_ecms_tier, factor4_adjusted_tier) != "T0"),
          "SELFTEST4: ECMS/PROGENy 共享 RNA 越过 T2 上限。")

  # 5) 非同患者代谢/微生物模块不得越过 T2 或进入宿主轴。
  metabolite_tier <- assign_metabolite_tier("robust_rank_and_scale_fdr")
  microbe_tier <- assign_microbe_beta_tier(1e-30, 1e-25, 0.90)
  axis_source_fixture <- c("PR001876", "PRJNA766558")
  fail_if(any(c(metabolite_tier, microbe_tier) %chin% c("T3", "T4")) ||
            any(!is_forbidden_unpaired_host_source(axis_source_fixture)),
          "SELFTEST5: 非同患者模块越级或进入宿主轴。")

  # 6) PRJNA 任一正式 artifact 缺失都不能视为完整发布。
  expected_prjna <- paste0("artifact_", sprintf("%02d", 1:12))
  observed_prjna <- expected_prjna[-12L]
  fail_if(!artifact_family_exact(expected_prjna, expected_prjna) ||
            artifact_family_exact(expected_prjna, observed_prjna),
          "SELFTEST6: PRJNA 部分发布被误认为完整。")

  # 7) 非法 logical 值必须 fail-closed。
  strict_failed <- FALSE
  tryCatch(
    as_logical_strict("FALSEE", "fixture_invalid_logical"),
    error = function(error) strict_failed <<- TRUE
  )
  fail_if(!strict_failed, "SELFTEST7: 非法 logical 值未 fail-closed。")

  data.table(
    fixture_id = sprintf("SELFTEST_%02d", 1:7),
    pass = TRUE,
    boundary = c(
      "zero_leakage_no_T3_T4",
      "one_low_margin_seed_no_extra_veto",
      "ecms_calibration_false_official78_only",
      "shared_rna_cap_T2",
      "unpaired_modules_cap_T2_no_host_axis",
      "prjna_partial_family_rejected",
      "invalid_logical_fail_closed"
    )
  )
}

if (self_test) {
  self_test_result <- run_self_tests()
  print(self_test_result)
  message("SELF_TEST_OK: 7/7")
  quit(save = "no", status = 0L)
}

stage_dir <- tempfile(
  pattern = ".escc_multiomics_integration_",
  tmpdir = work_intermediate_dir
)
dir_create(stage_dir)
# 成功发布后在脚本末尾显式删除；若运行失败则保留在 _work/intermediate/
# 供定位问题，不会污染正式 results/。

verify_manifest_inputs <- function(manifest_path, input_paths) {
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
  rows <- manifest[match(expected_relative, manifest$relative_path)]
  fail_if(
    nrow(rows) != length(input_paths) || anyNA(rows$relative_path),
    paste("manifest 未覆盖全部整合输入：", basename(manifest_path))
  )
  fail_if(
    any(!rows$status %chin% c("verified", "verified_after_atomic_publish")),
    paste("manifest 输入状态未验证：", basename(manifest_path))
  )
  fail_if(any(!file_exists(input_paths)), "manifest 已登记的整合输入文件缺失。")
  actual_size <- as.numeric(file_info(input_paths)$size)
  actual_sha <- vapply(
    input_paths,
    digest,
    FUN.VALUE = character(1),
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  )
  fail_if(
    any(actual_size != as.numeric(rows$file_size_bytes)) ||
      any(actual_sha != rows$sha256),
    paste("整合输入大小或 SHA256 与 manifest 不一致：", basename(manifest_path))
  )
  invisible(TRUE)
}

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

verify_exact_artifact_family <- function(
    manifest_path,
    formal_names,
    generation_script,
    generator_path,
    generator_sha_field,
    family_regex) {
  fail_if(!file_exists(manifest_path),
          paste("严格 artifact manifest 缺失：", manifest_path))
  manifest <- read_tsv(manifest_path)
  required_fields <- c(
    "relative_path", "file_size_bytes", "sha256", "status",
    "generation_script", generator_sha_field
  )
  require_columns(manifest, required_fields, basename(manifest_path))
  fail_if(anyDuplicated(manifest$relative_path) > 0L,
          paste("严格 manifest relative_path 重复：", basename(manifest_path)))

  expected_relative <- file.path("results", formal_names)
  fail_if(
    !artifact_family_exact(expected_relative, manifest$relative_path),
    paste("严格 manifest artifact family 不等于预锁定集合：",
          basename(manifest_path))
  )
  fail_if(any(manifest$status != "verified"),
          paste("严格 manifest 存在非 verified 行：", basename(manifest_path)))
  fail_if(any(manifest$generation_script != generation_script),
          paste("严格 manifest generation_script 不一致：",
                basename(manifest_path)))

  expected_generator_sha <- sha256_file(generator_path)
  observed_generator_sha <- as.character(manifest[[generator_sha_field]])
  fail_if(anyNA(observed_generator_sha) ||
            any(observed_generator_sha != expected_generator_sha),
          paste("严格 manifest 与当前上游脚本 SHA256 不一致：",
                basename(manifest_path)))

  artifact_paths <- file.path(project_root, manifest$relative_path)
  fail_if(any(!file_exists(artifact_paths)),
          paste("严格 manifest 登记 artifact 缺失：", basename(manifest_path)))
  observed_size <- as.numeric(file_info(artifact_paths)$size)
  observed_sha <- vapply(artifact_paths, sha256_file, character(1))
  fail_if(any(observed_size != as.numeric(manifest$file_size_bytes)) ||
            any(observed_sha != manifest$sha256),
          paste("严格 artifact family 大小/SHA256 不一致：",
                basename(manifest_path)))

  family_files <- basename(dir_ls(results_dir, type = "file", fail = FALSE))
  family_files <- family_files[grepl(family_regex, family_files, perl = TRUE)]
  expected_family_files <- c(formal_names, basename(manifest_path))
  fail_if(!artifact_family_exact(expected_family_files, family_files),
          paste("正式 results 中 artifact family 有缺失或 manifest 外额外文件：",
                basename(manifest_path)))

  manifest[match(expected_relative, manifest$relative_path)]
}

stage_tsv <- function(object, filename) {
  fwrite(object, file.path(stage_dir, filename), sep = "\t", quote = FALSE, na = "")
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

leakage_formal_names <- c(
  "tcga_escc_driver_state_leakage_model_plan.tsv",
  "tcga_escc_driver_state_leakage_factor_alignment.tsv",
  "tcga_escc_driver_state_leakage_associations.tsv",
  "tcga_escc_driver_state_leakage_edge_stability.tsv",
  "tcga_escc_driver_state_leakage_candidate_summary.tsv",
  "tcga_escc_driver_state_external_projection_status.tsv",
  "tcga_escc_driver_state_leakage_summary.md"
)
leakage_manifest_path <- file.path(
  results_dir, "tcga_escc_driver_state_leakage_artifact_manifest.tsv"
)
ecms_formal_names <- c(
  "tcga_escc_ecms_patient_probabilities.tsv",
  "tcga_escc_ecms_projection_calibration.tsv",
  "tcga_escc_ecms_projection_qa.tsv",
  "tcga_escc_ecms_factor_associations.tsv",
  "tcga_escc_ecms_adjusted_factor_progeny_associations.tsv",
  "tcga_escc_ecms_projection_summary.md"
)
ecms_manifest_path <- file.path(
  results_dir, "tcga_escc_ecms_projection_artifact_manifest.tsv"
)
prjna_formal_names <- c(
  "prjna766558_dada2_phyloseq.rds",
  "prjna766558_dada2_sample_qc.tsv",
  "prjna766558_dada2_asv_counts.tsv",
  "prjna766558_dada2_taxonomy.tsv",
  "prjna766558_dada2_asv_length_qc.tsv",
  "prjna766558_dada2_alpha_diversity.tsv",
  "prjna766558_dada2_alpha_paired_tests.tsv",
  "prjna766558_dada2_beta_tests.tsv",
  "prjna766558_dada2_genus_paired_differential.tsv",
  "prjna766558_dada2_ancombc2_sensitivity.tsv",
  "prjna766558_dada2_pipeline_sensitivity.tsv",
  "prjna766558_dada2_summary.md"
)
prjna_manifest_path <- file.path(
  results_dir, "prjna766558_dada2_artifact_manifest.tsv"
)

required_inputs <- list(
  driver = file.path(results_dir, "tcga_escc_driver_candidate_screen.tsv"),
  state_edges = file.path(results_dir, "tcga_escc_driver_state_network_edges.tsv"),
  leakage_edge_stability = file.path(
    results_dir, "tcga_escc_driver_state_leakage_edge_stability.tsv"
  ),
  leakage_candidate_summary = file.path(
    results_dir, "tcga_escc_driver_state_leakage_candidate_summary.tsv"
  ),
  leakage_external_projection = file.path(
    results_dir, "tcga_escc_driver_state_external_projection_status.tsv"
  ),
  ecms_patient_probabilities = file.path(
    results_dir, "tcga_escc_ecms_patient_probabilities.tsv"
  ),
  ecms_calibration = file.path(
    results_dir, "tcga_escc_ecms_projection_calibration.tsv"
  ),
  ecms_qa = file.path(results_dir, "tcga_escc_ecms_projection_qa.tsv"),
  ecms_factor_associations = file.path(
    results_dir, "tcga_escc_ecms_factor_associations.tsv"
  ),
  ecms_adjusted_associations = file.path(
    results_dir, "tcga_escc_ecms_adjusted_factor_progeny_associations.tsv"
  ),
  cao = file.path(results_dir, "cao2020_cross_layer_candidate_summary.tsv"),
  heterogeneity_cluster = file.path(
    results_dir, "tcga_escc_heterogeneity_cluster_evaluation.tsv"
  ),
  factor_pathways = file.path(results_dir, "tcga_escc_mofa_pathway_associations.tsv"),
  factor_weights = file.path(results_dir, "tcga_escc_mofa_top_weights.tsv"),
  factor_variance = file.path(results_dir, "tcga_escc_mofa_variance_explained.tsv"),
  factor_level_qc = file.path(results_dir, "tcga_escc_mofa_level_factor_qc.tsv"),
  metabolite_candidates = file.path(
    results_dir, "pr001876_targeted_ms_candidate_metabolites.tsv"
  ),
  metabolite_inventory = file.path(
    results_dir, "pr001876_targeted_ms_analysis_inventory.tsv"
  )
)
fail_if(
  any(!file_exists(unlist(required_inputs))),
  paste(
    "缺少正式整合输入：",
    paste(unlist(required_inputs)[!file_exists(unlist(required_inputs))], collapse = ", ")
  )
)

verify_manifest_inputs(
  file.path(results_dir, "tcga_escc_driver_core_artifact_manifest.tsv"),
  required_inputs$driver
)
verify_manifest_inputs(
  file.path(results_dir, "tcga_escc_driver_state_artifact_manifest.tsv"),
  required_inputs$state_edges
)
leakage_artifact_manifest <- verify_exact_artifact_family(
  leakage_manifest_path,
  leakage_formal_names,
  "scripts/21_audit_tcga_escc_driver_state_leakage.R",
  file.path(project_root, "scripts", "21_audit_tcga_escc_driver_state_leakage.R"),
  "generation_script_sha256",
  "^(tcga_escc_driver_state_leakage_.*|tcga_escc_driver_state_external_projection_status[.]tsv)$"
)
ecms_artifact_manifest <- verify_exact_artifact_family(
  ecms_manifest_path,
  ecms_formal_names,
  "scripts/23_project_tcga_escc_ecms.R",
  file.path(project_root, "scripts", "23_project_tcga_escc_ecms.R"),
  "execution_script_sha256",
  "^tcga_escc_ecms_.*"
)
verify_manifest_inputs(
  file.path(results_dir, "cao2020_cross_layer_artifact_manifest.tsv"),
  required_inputs$cao
)
verify_manifest_inputs(
  file.path(results_dir, "tcga_escc_heterogeneity_artifact_manifest.tsv"),
  unlist(required_inputs[c(
    "heterogeneity_cluster", "factor_pathways", "factor_weights",
    "factor_variance", "factor_level_qc"
  )])
)
verify_manifest_inputs(
  file.path(results_dir, "pr001876_targeted_ms_artifact_manifest.tsv"),
  unlist(required_inputs[c("metabolite_candidates", "metabolite_inventory")])
)

message("[1/7] 读取并冻结正式来源表")
driver <- read_tsv(required_inputs$driver)
state_edges <- read_tsv(required_inputs$state_edges)
leakage_edge_stability <- read_tsv(required_inputs$leakage_edge_stability)
leakage_candidate_summary <- read_tsv(required_inputs$leakage_candidate_summary)
leakage_external_projection <- read_tsv(
  required_inputs$leakage_external_projection
)
ecms_patient_probabilities <- read_tsv(
  required_inputs$ecms_patient_probabilities
)
ecms_calibration <- read_tsv(required_inputs$ecms_calibration)
ecms_qa <- read_tsv(required_inputs$ecms_qa)
ecms_factor_associations <- read_tsv(
  required_inputs$ecms_factor_associations
)
ecms_adjusted_associations <- read_tsv(
  required_inputs$ecms_adjusted_associations
)
cao <- read_tsv(required_inputs$cao)
cluster_evaluation <- read_tsv(required_inputs$heterogeneity_cluster)
factor_pathways <- read_tsv(required_inputs$factor_pathways)
factor_weights <- read_tsv(required_inputs$factor_weights)
factor_variance <- read_tsv(required_inputs$factor_variance)
factor_level_qc <- read_tsv(required_inputs$factor_level_qc)
metabolite_candidates <- read_tsv(required_inputs$metabolite_candidates)
metabolite_inventory <- read_tsv(required_inputs$metabolite_inventory)

require_columns(driver, c(
  "gene_id", "gene_name", "gene_type", "decision", "primary_candidate_route",
  "decision_basis", "conclusion_ceiling", "required_next_validation",
  "mutated_patients", "mutation_frequency_96", "high_confidence_mutation_pattern",
  "recurrent_mutation", "recurrent_high_level_cnv", "strong_dosage",
  "conditional_dosage", "precise_reverse_dosage", "multi_event_convergence",
  "evidence_unit_count", "spearman_rho_relative_cnv", "p_relative_cnv",
  "q_relative_cnv", "discordance_code", "test_eligible", "dosage_evidence",
  "amplification_frequency", "homozygous_deletion_frequency",
  "recurrent_cnv_class", "relative_log2_iqr"
), "tcga_escc_driver_candidate_screen.tsv")
require_columns(state_edges, c(
  "edge_id", "source_node", "source_gene_id", "source_gene_name",
  "event_type", "target_node",
  "target_type", "effect_measure", "effect", "effect_ci95_low",
  "effect_ci95_high", "p_value", "q_value", "association_status",
  "level_factor_soft_flag", "target_model_overlap_class", "independence_group",
  "countable_as_independent_validation", "conclusion_ceiling",
  "required_next_validation"
), "tcga_escc_driver_state_network_edges.tsv")
require_columns(leakage_edge_stability, c(
  "original_edge_id", "gene_id", "gene_name", "event_type",
  "reference_factor", "scenario_id", "planned_seed_count",
  "completed_converged_seed_count", "reliable_factor_match_rate",
  "direction_retention_rate", "magnitude_retention_rate",
  "gate_retention_rate", "low_margin_count", "retention_decision"
), "tcga_escc_driver_state_leakage_edge_stability.tsv")
require_columns(leakage_candidate_summary, c(
  "original_edge_id", "gene_id", "gene_name", "event_type",
  "reference_factor", "original_association_status",
  "original_level_factor_soft_flag", "planned_seed_count",
  "completed_converged_seed_count", "baseline_reliable_match_rate",
  "drop_both_reliable_match_rate", "drop_both_direction_retention_rate",
  "drop_both_magnitude_retention_rate", "drop_both_gate_retention_rate",
  "drop_both_low_margin_count", "drop_both_retention_decision",
  "leakage_gate_evaluable", "leakage_gate_pass",
  "countable_for_T3_T4", "gate_failure_reason",
  "post_audit_maximum_status", "leakage_control_decision",
  "p_q_scope", "conclusion_ceiling", "required_next_validation"
), "tcga_escc_driver_state_leakage_candidate_summary.tsv")
require_columns(leakage_external_projection, c(
  "benchmark", "benchmark_role", "locked_public_projection_weights_available",
  "projection_artifact_manifest_verified", "projection_status", "model_commit",
  "repository_anchor_patient_count", "extension_calibration_pass",
  "primary_patient_count", "normalization_scope", "single_sample_valid",
  "used_in_mofa_training_or_leakage_gate", "pseudo_labels_generated"
), "tcga_escc_driver_state_external_projection_status.tsv")
require_columns(ecms_patient_probabilities, c(
  "patient_id", "in_official_78", "official_anchor_label",
  "gdc_projection_label", "resolved_ecms_label", "label_source",
  "extension_status", "eligible_for_primary_association",
  "resolved_prob_ECMS1", "resolved_prob_ECMS2", "resolved_prob_ECMS3",
  "resolved_prob_ECMS4", "resolved_margin_custom", "low_margin_custom_flag",
  "margin_is_official_threshold", "low_margin_rejects_label",
  "single_sample_classifier_claim", "pseudo_label_generated"
), "tcga_escc_ecms_patient_probabilities.tsv")
require_columns(ecms_calibration, c(
  "metric", "observed", "threshold", "direction", "gate_role", "pass", "notes"
), "tcga_escc_ecms_projection_calibration.tsv")
require_columns(ecms_qa, c(
  "check_id", "category", "hard_gate", "expected", "observed", "pass", "notes"
), "tcga_escc_ecms_projection_qa.tsv")
require_columns(ecms_factor_associations, c(
  "analysis_scope", "association_role", "is_primary_scope", "n_patients",
  "factor", "ecms1_n", "ecms2_n", "ecms3_n", "ecms4_n",
  "ecms1_mean", "ecms2_mean", "ecms3_mean", "ecms4_mean",
  "anova_f", "anova_p_value", "eta_squared", "kruskal_h",
  "kruskal_p_value", "extension_calibration_pass",
  "shared_rna_representation", "independent_validation", "evidence_ceiling",
  "interpretation", "anova_q_value", "kruskal_q_value"
), "tcga_escc_ecms_factor_associations.tsv")
require_columns(ecms_adjusted_associations, c(
  "analysis_scope", "association_role", "is_primary_scope", "n_patients",
  "pathway", "factor", "reduced_model", "full_model",
  "factor_beta_standardized", "factor_standard_error", "factor_t_value",
  "factor_coefficient_p_value", "reduced_r_squared", "full_r_squared",
  "delta_r_squared", "partial_r_squared", "incremental_f",
  "incremental_p_value", "extension_calibration_pass",
  "shared_rna_representation", "independent_validation", "evidence_ceiling",
  "interpretation", "incremental_q_value", "factor_coefficient_q_value"
), "tcga_escc_ecms_adjusted_factor_progeny_associations.tsv")
require_columns(cao, c(
  "tcga_gene_id", "tcga_gene_name", "promoter_median_paired_delta_beta",
  "promoter_region_evaluation_status", "promoter_methylation_effect_class",
  "limma_log_fc_pc1", "limma_p_pc1", "limma_q_candidate_pc1",
  "rna_effect_class", "protein_measurement_status",
  "median_patient_log2_tumor_vs_normal", "promoter_rna_n",
  "promoter_rna_rho", "promoter_rna_p", "rna_protein_n",
  "rna_protein_rho", "rna_protein_p", "promoter_rna_relationship",
  "rna_protein_relationship", "cross_layer_class", "statistical_support_class",
  "available_layer_count", "fdr_supported_layer_count", "conclusion_ceiling",
  "required_follow_up"
), "cao2020_cross_layer_candidate_summary.tsv")
require_columns(cluster_evaluation, c(
  "k", "snf_mean_silhouette", "mofa_consensus_pac",
  "snf_mofa_adjusted_rand", "selected_for_description", "decision"
), "tcga_escc_heterogeneity_cluster_evaluation.tsv")
require_columns(factor_pathways, c(
  "factor", "pathway", "n", "spearman_rho", "p_value", "q_value",
  "abs_spearman_rho"
), "tcga_escc_mofa_pathway_associations.tsv")
require_columns(factor_weights, c(
  "view", "factor", "feature", "value", "abs_weight"
), "tcga_escc_mofa_top_weights.tsv")
require_columns(factor_variance, c(
  "r2_per_factor.view", "r2_per_factor.factor", "r2_per_factor.value"
), "tcga_escc_mofa_variance_explained.tsv")
require_columns(factor_level_qc, c(
  "factor", "view", "spearman_rho_with_view_mean", "q_value",
  "level_factor_soft_flag"
), "tcga_escc_mofa_level_factor_qc.tsv")
require_columns(metabolite_candidates, c(
  "analysis_id", "feature_id", "study_id", "platform", "ion_mode",
  "evidence_family", "gate_independence_group", "metabolite_name",
  "annotation_identity_status", "limma_effect_transformed", "limma_p",
  "limma_q_evidence_family", "direction", "final_candidate_tier",
  "conditional_reason", "evidence_claim_ceiling"
), "pr001876_targeted_ms_candidate_metabolites.tsv")
require_columns(metabolite_inventory, c(
  "analysis_id", "n_escc", "n_normal", "gate_independence_group",
  "perfect_group_run_order_block", "current_role", "candidates_total"
), "pr001876_targeted_ms_analysis_inventory.tsv")

# 所有会影响证据层级的外部 logical 字段必须在任何裁决前严格解析。
driver_logical_fields <- c(
  "high_confidence_mutation_pattern", "recurrent_mutation",
  "recurrent_high_level_cnv", "strong_dosage", "conditional_dosage",
  "precise_reverse_dosage", "multi_event_convergence", "test_eligible"
)
for (column in driver_logical_fields) {
  set(driver, j = column, value = as_logical_strict(
    driver[[column]], paste0("driver.", column)
  ))
}
state_edges[, level_factor_soft_flag := as_logical_strict(
  level_factor_soft_flag, "state_edges.level_factor_soft_flag"
)]
state_edges[, countable_as_independent_validation := as_logical_strict(
  countable_as_independent_validation,
  "state_edges.countable_as_independent_validation"
)]
for (column in c(
  "original_level_factor_soft_flag", "leakage_gate_evaluable",
  "leakage_gate_pass", "countable_for_T3_T4"
)) {
  set(leakage_candidate_summary, j = column, value = as_logical_strict(
    leakage_candidate_summary[[column]],
    paste0("leakage_candidate_summary.", column)
  ))
}
cluster_evaluation[, selected_for_description := as_logical_strict(
  selected_for_description, "cluster_evaluation.selected_for_description"
)]
factor_level_qc[, level_factor_soft_flag := as_logical_strict(
  level_factor_soft_flag, "factor_level_qc.level_factor_soft_flag"
)]
metabolite_inventory[, perfect_group_run_order_block := as_logical_strict(
  perfect_group_run_order_block,
  "metabolite_inventory.perfect_group_run_order_block"
)]

locked_ecms_commit <- "8a49db65f86f2c1e2a50d9bc4cb742ccc868e4c6"
fail_if(nrow(leakage_external_projection) != 1L,
        "script21 external projection status 必须恰好一行。")
for (column in c(
  "locked_public_projection_weights_available",
  "projection_artifact_manifest_verified", "extension_calibration_pass",
  "single_sample_valid", "used_in_mofa_training_or_leakage_gate",
  "pseudo_labels_generated"
)) {
  set(leakage_external_projection, j = column, value = as_logical_strict(
    leakage_external_projection[[column]],
    paste0("leakage_external_projection.", column)
  ))
}
fail_if(
  leakage_external_projection$model_commit[[1L]] != locked_ecms_commit ||
    leakage_external_projection$repository_anchor_patient_count[[1L]] != 78L ||
    leakage_external_projection$primary_patient_count[[1L]] != 78L ||
    leakage_external_projection$extension_calibration_pass[[1L]] ||
    !leakage_external_projection$locked_public_projection_weights_available[[1L]] ||
    !leakage_external_projection$projection_artifact_manifest_verified[[1L]] ||
    leakage_external_projection$used_in_mofa_training_or_leakage_gate[[1L]] ||
    leakage_external_projection$single_sample_valid[[1L]] ||
    leakage_external_projection$pseudo_labels_generated[[1L]] ||
    leakage_external_projection$projection_status[[1L]] !=
      "available_downstream_verified_not_used_in_leakage_gate",
  "script21 ECMS 下游可用但不进入 MOFA/泄漏门禁的边界不满足。"
)

require_columns(ecms_artifact_manifest, c(
  "model_commit", "extension_calibration_pass", "primary_patient_count"
), "tcga_escc_ecms_projection_artifact_manifest.tsv")
ecms_artifact_manifest[, extension_calibration_pass := as_logical_strict(
  extension_calibration_pass, "ecms_manifest.extension_calibration_pass"
)]
fail_if(any(ecms_artifact_manifest$model_commit != locked_ecms_commit) ||
          any(ecms_artifact_manifest$extension_calibration_pass) ||
          any(as.integer(ecms_artifact_manifest$primary_patient_count) != 78L),
        "ECMS manifest 必须锁定官方78主分析且 extension calibration=FALSE。")

for (column in c(
  "in_official_78", "eligible_for_primary_association",
  "low_margin_custom_flag", "margin_is_official_threshold",
  "low_margin_rejects_label", "single_sample_classifier_claim",
  "pseudo_label_generated"
)) {
  set(ecms_patient_probabilities, j = column, value = as_logical_strict(
    ecms_patient_probabilities[[column]],
    paste0("ecms_patient_probabilities.", column)
  ))
}
ecms_probability_columns <- paste0("resolved_prob_ECMS", 1:4)
ecms_probability_matrix <- as.matrix(
  ecms_patient_probabilities[, ..ecms_probability_columns]
)
fail_if(
  nrow(ecms_patient_probabilities) != 94L ||
    uniqueN(ecms_patient_probabilities$patient_id) != 94L ||
    sum(ecms_patient_probabilities$in_official_78) != 78L ||
    sum(ecms_patient_probabilities$eligible_for_primary_association) != 78L ||
    any(ecms_patient_probabilities[
      in_official_78 == FALSE, eligible_for_primary_association
    ]) ||
    any(ecms_patient_probabilities[
      in_official_78 == TRUE,
      label_source !=
        "locked_model_prediction_on_repository_bundled_tcga78_anchor" |
        resolved_ecms_label != official_anchor_label
    ]) ||
    any(ecms_patient_probabilities$margin_is_official_threshold) ||
    any(ecms_patient_probabilities$low_margin_rejects_label) ||
    any(ecms_patient_probabilities$single_sample_classifier_claim) ||
    any(ecms_patient_probabilities$pseudo_label_generated) ||
    any(!is.finite(ecms_probability_matrix)) ||
    any(abs(rowSums(ecms_probability_matrix) - 1) > 1e-8),
  "ECMS patient table 未严格满足 94/78 anchor、概率或伪标签边界。"
)

ecms_calibration[, pass_strict := as_logical_strict(
  pass, "ecms_calibration.pass", allow_na = TRUE
)]
calibration_gate_rows <- ecms_calibration[
  gate_role == "prelocked_extension_gate"
]
overall_calibration_row <- ecms_calibration[
  metric == "overall_extension_calibration_pass"
]
expected_calibration_metrics <- c(
  "exact_label_agreement", "cohen_kappa", "adjusted_rand_index",
  "median_class_probability_spearman"
)
fail_if(
  nrow(calibration_gate_rows) != 4L || nrow(overall_calibration_row) != 1L ||
    !setequal(calibration_gate_rows$metric, expected_calibration_metrics) ||
    any(calibration_gate_rows$direction != ">=") ||
    any(calibration_gate_rows$pass_strict !=
          (calibration_gate_rows$observed >= calibration_gate_rows$threshold)) ||
    overall_calibration_row$pass_strict[[1L]] ||
    overall_calibration_row$observed[[1L]] != 0 ||
    overall_calibration_row$threshold[[1L]] != 1 ||
    overall_calibration_row$direction[[1L]] != "all_four_prelocked_gates",
  "ECMS 预锁定校准门禁或 overall=FALSE 边界不一致。"
)
ecms_extension_calibration_pass <- FALSE
ecms_primary_patient_count <- 78L

ecms_qa[, hard_gate_strict := as_logical_strict(
  hard_gate, "ecms_qa.hard_gate"
)]
ecms_qa[, pass_strict := as_logical_strict(
  pass, "ecms_qa.pass"
)]
expected_ecms_qa_ids <- c(
  "locked_commit", "model_bytes", "model_manifest_full_sha",
  "upstream_artifact_manifests", "model_object_set", "randomforest_structure",
  "feature_count", "feature_unique", "feature_list_sha256",
  "builtin_tcga_shape", "builtin_tcga_gene_scaling",
  "builtin_prediction_counts", "five_layer_core_count", "official_id_mapping",
  "gdc_tpm_assay", "ensembl_strip_unique", "feature_coverage",
  "feature_nonzero_variance", "gdc_projection_shape", "probability_row_sums",
  "official_anchor_precedence", "extension_gate", "margin_boundary",
  "pseudo_label_prohibition", "clinical_classifier_boundary"
)
fail_if(
  uniqueN(ecms_qa$check_id) != nrow(ecms_qa) ||
    !all(expected_ecms_qa_ids %chin% ecms_qa$check_id) ||
    any(ecms_qa[hard_gate_strict == TRUE, !pass_strict]) ||
    ecms_qa[check_id == "extension_gate", .N] != 1L ||
    ecms_qa[check_id == "extension_gate", hard_gate_strict] ||
    ecms_qa[check_id == "extension_gate", pass_strict] ||
    any(ecms_qa[
      hard_gate_strict == FALSE & check_id != "extension_gate",
      !pass_strict
    ]),
  "ECMS QA 必须覆盖预锁定检查，所有硬门禁通过，且只允许软 extension gate 失败。"
)

for (object_name in c("ecms_factor_associations", "ecms_adjusted_associations")) {
  object <- get(object_name)
  for (column in c(
    "is_primary_scope", "extension_calibration_pass",
    "shared_rna_representation", "independent_validation"
  )) {
    set(object, j = column, value = as_logical_strict(
      object[[column]], paste0(object_name, ".", column)
    ))
  }
  if ("factor_level_soft_flag" %in% names(object)) {
    set(object, j = "factor_level_soft_flag", value = as_logical_strict(
      object$factor_level_soft_flag,
      paste0(object_name, ".factor_level_soft_flag")
    ))
  }
  assign(object_name, object)
}
expected_factors <- paste0("Factor", 1:8)
ecms_factor_primary <- ecms_factor_associations[
  analysis_scope == "official_anchor_78" & association_role == "primary" &
    is_primary_scope
]
ecms_factor_hybrid_boundary <- ecms_factor_associations[
  analysis_scope == "hybrid_94_extension"
]
ecms_adjusted_primary <- ecms_adjusted_associations[
  analysis_scope == "official_anchor_78" & association_role == "primary" &
    is_primary_scope
]
ecms_adjusted_hybrid_boundary <- ecms_adjusted_associations[
  analysis_scope == "hybrid_94_extension"
]
expected_pathways <- sort(unique(factor_pathways$pathway))
if ("factor_level_soft_flag" %in% names(ecms_factor_primary)) {
  ecms_factor_flag_check <- unique(ecms_factor_primary[, .(
    factor, ecms_output_level_flag = factor_level_soft_flag
  )])
  factor_level_flag_check <- merge(
    ecms_factor_flag_check, factor_level_qc[, .(
      project_level_flag = any(level_factor_soft_flag)
    ), by = factor],
    by = "factor", all = TRUE
  )
  fail_if(anyNA(factor_level_flag_check$ecms_output_level_flag) ||
            anyNA(factor_level_flag_check$project_level_flag) ||
            any(factor_level_flag_check$ecms_output_level_flag !=
                  factor_level_flag_check$project_level_flag),
          "script23 factor_level_soft_flag 与项目冻结 QC 不一致。")
}
if ("factor_level_soft_flag" %in% names(ecms_adjusted_primary)) {
  ecms_adjusted_flag_check <- unique(ecms_adjusted_primary[, .(
    factor, ecms_output_level_flag = factor_level_soft_flag
  )])
  factor_level_flag_check <- merge(
    ecms_adjusted_flag_check, factor_level_qc[, .(
      project_level_flag = any(level_factor_soft_flag)
    ), by = factor],
    by = "factor", all = TRUE
  )
  fail_if(anyNA(factor_level_flag_check$ecms_output_level_flag) ||
            anyNA(factor_level_flag_check$project_level_flag) ||
            any(factor_level_flag_check$ecms_output_level_flag !=
                  factor_level_flag_check$project_level_flag),
          "script23 adjusted factor level flag 与项目冻结 QC 不一致。")
}
fail_if(
  nrow(ecms_factor_primary) != 8L ||
    !setequal(ecms_factor_primary$factor, expected_factors) ||
    any(ecms_factor_primary$n_patients != 78L) ||
    nrow(ecms_factor_hybrid_boundary) != 8L ||
    any(ecms_factor_hybrid_boundary$is_primary_scope) ||
    any(ecms_factor_hybrid_boundary$n_patients != 94L) ||
    nrow(ecms_adjusted_primary) != 112L ||
    uniqueN(ecms_adjusted_primary, by = c("factor", "pathway")) != 112L ||
    !setequal(ecms_adjusted_primary$factor, expected_factors) ||
    !setequal(ecms_adjusted_primary$pathway, expected_pathways) ||
    any(ecms_adjusted_primary$n_patients != 78L) ||
    nrow(ecms_adjusted_hybrid_boundary) != 112L ||
    any(ecms_adjusted_hybrid_boundary$is_primary_scope) ||
    any(ecms_adjusted_hybrid_boundary$n_patients != 94L) ||
    any(ecms_factor_associations$extension_calibration_pass) ||
    any(ecms_adjusted_associations$extension_calibration_pass) ||
    any(!ecms_factor_associations$shared_rna_representation) ||
    any(!ecms_adjusted_associations$shared_rna_representation) ||
    any(ecms_factor_associations$independent_validation) ||
    any(ecms_adjusted_associations$independent_validation),
  "ECMS factor/adjusted association 未严格满足 official78 primary 和 hybrid94 boundary 契约。"
)

driver <- driver[decision %chin% c(
  "strong_patient_level_candidate", "conditional_candidate"
)]
fail_if(!nrow(driver), "没有可进入整合层的 strong/conditional driver 候选。")
fail_if(anyDuplicated(driver$gene_id) > 0L, "driver 候选 gene_id 重复。")
strong_pairs <- unique(driver[
  decision == "strong_patient_level_candidate",
  .(gene_id, gene_name)
])
cao_pairs <- unique(cao[, .(
  gene_id = tcga_gene_id,
  gene_name = tcga_gene_name
)])
fail_if(
  nrow(fsetdiff(strong_pairs, cao_pairs)) > 0L ||
    nrow(fsetdiff(cao_pairs, strong_pairs)) > 0L,
  "Cao 跨层校准 (gene_id, gene_name) 与当前 TCGA strong 候选不一致。"
)
state_pairs <- unique(state_edges[, .(
  gene_id = source_gene_id,
  gene_name = source_gene_name
)])
fail_if(
  nrow(fsetdiff(state_pairs, strong_pairs)) > 0L,
  "driver-state 网络出现当前 strong 候选之外的 (gene_id, gene_name)。"
)
strong_names <- strong_pairs$gene_name

# 泄漏审计只对 41 条 MOFA 原边做 event-blind 对齐与 drop-both
# 门禁；PROGENy 边仍可作同 TCGA 内部解释，但不得替代该门禁。
fail_if(
  state_edges[target_type == "MOFA_factor", .N] != 41L ||
    nrow(leakage_candidate_summary) != 41L ||
    uniqueN(leakage_candidate_summary$original_edge_id) != 41L,
  "MOFA 原边或泄漏门禁必须为 41 条唯一记录。"
)
fail_if(
  !setequal(
    state_edges[target_type == "MOFA_factor", edge_id],
    leakage_candidate_summary$original_edge_id
  ),
  "script19 的 41 条 MOFA 边与 script21 逐原边门禁不一致。"
)
leakage_gate_logical_fields <- c(
  "original_level_factor_soft_flag", "leakage_gate_evaluable",
  "leakage_gate_pass", "countable_for_T3_T4"
)
fail_if(
  any(vapply(
    leakage_gate_logical_fields,
    function(column) anyNA(leakage_candidate_summary[[column]]),
    FUN.VALUE = logical(1)
  )),
  "script21 逐原边门禁关键逻辑字段含 NA。"
)
countable_gate <- as_logical_safe(
  leakage_candidate_summary$countable_for_T3_T4
)
passed_gate <- as_logical_safe(leakage_candidate_summary$leakage_gate_pass)
evaluable_gate <- as_logical_safe(
  leakage_candidate_summary$leakage_gate_evaluable
)
original_level_flag <- as_logical_safe(
  leakage_candidate_summary$original_level_factor_soft_flag
)
fail_if(
  any(countable_gate != passed_gate),
  "countable_for_T3_T4 与 leakage_gate_pass 不一致。"
)
expected_countable_gate <- leakage_gate_contract(
  leakage_candidate_summary$original_association_status,
  original_level_flag,
  evaluable_gate,
  leakage_candidate_summary$completed_converged_seed_count,
  leakage_candidate_summary$planned_seed_count,
  leakage_candidate_summary$baseline_reliable_match_rate,
  leakage_candidate_summary$drop_both_reliable_match_rate,
  leakage_candidate_summary$drop_both_direction_retention_rate,
  leakage_candidate_summary$drop_both_magnitude_retention_rate,
  leakage_candidate_summary$drop_both_gate_retention_rate
)
fail_if(anyNA(expected_countable_gate) ||
          any(countable_gate != expected_countable_gate),
        "script21 countable gate 与预锁定 2/3 合同不一致。")
fail_if(
  any(countable_gate & (
    leakage_candidate_summary$original_association_status !=
      "within_tcga_supported" |
      original_level_flag |
      !evaluable_gate |
      leakage_candidate_summary$completed_converged_seed_count !=
        leakage_candidate_summary$planned_seed_count |
      !is.finite(leakage_candidate_summary$baseline_reliable_match_rate) |
      leakage_candidate_summary$baseline_reliable_match_rate < 2 / 3 |
      !is.finite(leakage_candidate_summary$drop_both_reliable_match_rate) |
      leakage_candidate_summary$drop_both_reliable_match_rate < 2 / 3 |
      !is.finite(leakage_candidate_summary$drop_both_direction_retention_rate) |
      leakage_candidate_summary$drop_both_direction_retention_rate < 2 / 3 |
      !is.finite(leakage_candidate_summary$drop_both_magnitude_retention_rate) |
      leakage_candidate_summary$drop_both_magnitude_retention_rate < 2 / 3 |
      !is.finite(leakage_candidate_summary$drop_both_gate_retention_rate) |
      leakage_candidate_summary$drop_both_gate_retention_rate < 2 / 3
  )),
  "countable_for_T3_T4 越过了原始层级、可评估性、收敛或 2/3 稳定性门禁。"
)
fail_if(
  any(
    !is.finite(leakage_candidate_summary$drop_both_low_margin_count) |
      leakage_candidate_summary$drop_both_low_margin_count < 0L |
      leakage_candidate_summary$drop_both_low_margin_count >
        leakage_candidate_summary$planned_seed_count
  ),
  "drop-both low-margin 计数超出预锁定 seed 范围。"
)
fail_if(
  any(leakage_candidate_summary$p_q_scope !=
        "TCGA_internal_descriptive_only"),
  "script21 p/q 证据边界不是 TCGA 内部描述。"
)
leakage_gate <- leakage_candidate_summary[, .(
  edge_id = original_edge_id,
  leakage_gene_id = gene_id,
  leakage_gene_name = gene_name,
  leakage_event_type = event_type,
  leakage_reference_factor = reference_factor,
  leakage_original_association_status = original_association_status,
  leakage_original_level_factor_soft_flag = original_level_factor_soft_flag,
  leakage_planned_seed_count = planned_seed_count,
  leakage_completed_converged_seed_count = completed_converged_seed_count,
  baseline_reliable_match_rate,
  drop_both_reliable_match_rate,
  drop_both_direction_retention_rate,
  drop_both_magnitude_retention_rate,
  drop_both_gate_retention_rate,
  drop_both_low_margin_count,
  drop_both_retention_decision,
  leakage_gate_evaluable,
  leakage_gate_pass,
  countable_for_T3_T4,
  leakage_gate_failure_reason = gate_failure_reason,
  post_audit_maximum_status,
  leakage_control_decision,
  leakage_p_q_scope = p_q_scope,
  leakage_conclusion_ceiling = conclusion_ceiling,
  leakage_required_next_validation = required_next_validation
)]
state_edges <- merge(
  state_edges, leakage_gate, by = "edge_id", all.x = TRUE, sort = FALSE
)
fail_if(
  any(state_edges[target_type == "MOFA_factor", is.na(leakage_gene_id)]) ||
    any(state_edges[target_type != "MOFA_factor", !is.na(leakage_gene_id)]),
  "script21 泄漏门禁未严格一对一映射 MOFA 原边。"
)
state_edges[, leakage_identity_mismatch :=
  target_type == "MOFA_factor" & (
    source_gene_id != leakage_gene_id |
      source_gene_name != leakage_gene_name |
      event_type != leakage_event_type |
      target_node != leakage_reference_factor |
      association_status != leakage_original_association_status |
      as_logical_safe(level_factor_soft_flag) !=
        as_logical_safe(leakage_original_level_factor_soft_flag)
  )]
fail_if(
  any(state_edges$leakage_identity_mismatch),
  "script19 与 script21 的 MOFA 原边身份/状态不一致。"
)
state_edges[, leakage_identity_mismatch := NULL]
state_edges[, leakage_controlled_support :=
  target_type == "MOFA_factor" & as_logical_safe(countable_for_T3_T4)]
state_edges[, integration_state_status := fcase(
  target_type == "MOFA_factor" & as_logical_safe(level_factor_soft_flag),
    "level_factor_background_not_gate",
  leakage_controlled_support,
    "leakage_controlled_within_tcga_supported",
  target_type == "MOFA_factor" &
    association_status == "within_tcga_supported",
    "raw_supported_leakage_not_retained_conditional",
  target_type == "MOFA_factor" &
    association_status == "within_tcga_conditional",
    "raw_conditional_not_upgradeable",
  target_type == "MOFA_factor",
    "raw_exploratory_not_upgradeable",
  target_type == "PROGENy_pathway" &
    association_status %chin% c("within_tcga_supported", "within_tcga_conditional"),
    "progeny_internal_interpretive_not_gate",
  default = "directional_exploratory"
)]
state_edges[, integration_state_tier := assign_state_edge_tier(
  target_type,
  level_factor_soft_flag,
  leakage_controlled_support,
  association_status
)]

tier_definitions <- data.table(
  evidence_tier = paste0("T", 0:4),
  chinese_label = c(
    "背景或不可评估", "保留探索", "条件候选", "当前证据范围内稳健候选",
    "高优先级验证对象"
  ),
  operational_definition = c(
    "技术不可评估、不可比较、层缺失或仅作背景；不把缺失写成阴性",
    "来源内方向线索或描述性关系，尚不足以形成正式候选轴",
    "单队列、同队列桥接、小样本方向性或可解释反向；带条件保留",
    "患者级严格统计或多种敏感性支持的来源内稳健候选；仍需独立复现",
    "跨层或跨来源收敛且已有可证伪验证路径；表示下一步优先级，不表示已验证因果"
  ),
  default_claim_ceiling = c(
    "仅背景或不可评估", "探索性线索", "条件性候选", "稳健候选而非因果机制",
    "高优先级验证候选而非已验证靶点"
  ),
  upgrade_requirement = c(
    "补齐身份、质量、可比层或可解释输入", "增加患者级稳定性或正交支持",
    "独立 ESCC 队列、患者级重抽样或明确同患者跨层支持",
    "独立事件复现、跨平台复核及必要的细胞来源定位",
    "独立队列复现并以扰动/功能实验检验预先定义方向"
  ),
  score_role = "数值分数仅作同层排序辅助；T 层级、硬门禁和结论上限优先"
)

message("[2/7] 构建 strong/conditional 驱动候选整合表")
state_gene_summary <- state_edges[, .(
  state_edge_count = .N,
  raw_state_supported_edges = sum(
    association_status == "within_tcga_supported"
  ),
  state_supported_edges = sum(leakage_controlled_support),
  state_conditional_edges = sum(association_status == "within_tcga_conditional"),
  state_exploratory_edges = sum(grepl("directional_exploratory", association_status)),
  supported_factor_targets = collapse_unique(
    target_node[leakage_controlled_support]
  ),
  leakage_not_retained_factor_targets = collapse_unique(
    target_node[target_type == "MOFA_factor" &
      association_status == "within_tcga_supported" &
      !leakage_controlled_support]
  ),
  conditional_factor_targets = collapse_unique(
    target_node[target_type == "MOFA_factor" &
                  association_status == "within_tcga_conditional"]
  ),
  conditional_pathway_targets = collapse_unique(
    target_node[target_type == "PROGENy_pathway" &
                  association_status == "within_tcga_conditional"]
  ),
  best_state_q = min_finite(q_value)
), by = .(gene_name = source_gene_name)]

state_best <- state_edges[
  order(-as.integer(leakage_controlled_support), q_value, p_value,
        na.last = TRUE), .SD[1L],
                          by = source_gene_name]
state_best <- state_best[, .(
  gene_name = source_gene_name,
  strongest_state_edge = paste(source_node, target_node, sep = " -> "),
  strongest_state_effect = effect,
  strongest_state_q = q_value,
  strongest_state_status = integration_state_status,
  strongest_state_raw_status = association_status,
  strongest_state_leakage_countable = leakage_controlled_support,
  strongest_state_leakage_failure_reason = leakage_gate_failure_reason
)]
state_gene_summary <- merge(
  state_gene_summary, state_best, by = "gene_name", all = TRUE, sort = FALSE
)

cao_driver <- cao[, .(
  gene_id = tcga_gene_id,
  gene_name = tcga_gene_name,
  cao_profile_status = "profiled_strong_candidate",
  cao_cross_layer_class = cross_layer_class,
  cao_statistical_support = statistical_support_class,
  cao_available_layer_count = available_layer_count,
  cao_fdr_supported_layer_count = fdr_supported_layer_count,
  cao_promoter_delta_beta = promoter_median_paired_delta_beta,
  cao_rna_log_fc = limma_log_fc_pc1,
  cao_rna_effect_class = rna_effect_class,
  cao_protein_status = protein_measurement_status,
  cao_protein_log2_tumor_vs_normal = median_patient_log2_tumor_vs_normal,
  cao_promoter_rna_relationship = promoter_rna_relationship,
  cao_rna_protein_relationship = rna_protein_relationship
)]

driver_integrated <- merge(
  driver, state_gene_summary, by = "gene_name", all.x = TRUE, sort = FALSE
)
driver_integrated <- merge(
  driver_integrated, cao_driver, by = c("gene_id", "gene_name"),
  all.x = TRUE, sort = FALSE
)

# 整合层重算“不同生物事件单元”，不直接沿用旧 screen 字段。
# 旧字段只将 >=5 例复发突变计为突变单元，会漏计少见但存在
# 多患者 hotspot/高置信功能模式的候选（如 GNAS、FBXW7）。CNV 复发与
# CNV–RNA dosage 是同一生物事件族，始终只计一个单元。
driver_integrated[, source_evidence_unit_count := as.integer(evidence_unit_count)]
driver_integrated[, mutation_event_unit := as.integer(
  as_logical_safe(high_confidence_mutation_pattern) |
    as_logical_safe(recurrent_mutation)
)]
driver_integrated[, cnv_expression_event_unit := as.integer(
  as_logical_safe(recurrent_high_level_cnv) &
    (as_logical_safe(strong_dosage) | as_logical_safe(conditional_dosage))
)]
driver_integrated[, distinct_event_unit_count :=
  mutation_event_unit + cnv_expression_event_unit]
# 保留原有下游语义入口，但其值明确更新为整合层重算结果。
driver_integrated[, evidence_unit_count := distinct_event_unit_count]
driver_integrated[is.na(cao_profile_status), `:=`(
  cao_profile_status = "not_profiled_current_strong_only",
  cao_cross_layer_class = "not_profiled_not_negative",
  cao_statistical_support = "not_profiled_not_negative",
  cao_protein_status = "not_profiled_not_negative"
)]
for (column in c(
  "state_edge_count", "raw_state_supported_edges", "state_supported_edges",
  "state_conditional_edges", "state_exploratory_edges"
)) {
  set(driver_integrated, which(is.na(driver_integrated[[column]])), column, 0L)
}

driver_integrated[, exact_driver_event_independent_validation_count := 0L]
driver_integrated[, independent_orthogonal_calibration_count := as.integer(
  cao_profile_status == "profiled_strong_candidate"
)]
driver_integrated[, quality_score_20 := fifelse(
  decision == "strong_patient_level_candidate", 18, 14
)]
driver_integrated[, within_layer_score_20 := fifelse(
  decision == "strong_patient_level_candidate",
  pmin(20, 16 + pmax(0, as.integer(evidence_unit_count) - 1L) * 4),
  11
)]
driver_integrated[, cross_omics_score_25 := fcase(
  cao_cross_layer_class == "three_layer_directional_hypothesis", 18,
  cao_cross_layer_class == "credible_reverse_retained", 13,
  cao_profile_status == "profiled_strong_candidate" &
    cao_rna_effect_class %chin% c("higher_in_tumor", "lower_in_tumor"), 8,
  cao_profile_status == "profiled_strong_candidate", 4,
  default = 0
)]
driver_integrated[, validation_score_20 := 0]
driver_integrated[, robustness_score_15 := pmin(
  15,
  3 * state_supported_edges +
    2 * state_conditional_edges +
    2 * as.integer(as_logical_safe(multi_event_convergence))
)]
driver_integrated[, auxiliary_score_100 := quality_score_20 +
                    within_layer_score_20 + cross_omics_score_25 +
                    validation_score_20 + robustness_score_15]
driver_integrated[, score_is_decision_gate := FALSE]

# T4 由可审计的结构触发，而不是由分数触发。当前要求：strong driver、
# 至少两个不同生物事件单元、来源内稳健状态边，以及 Cao 中可评估 RNA 方向。
driver_integrated[, high_priority_validation_rule :=
  decision == "strong_patient_level_candidate" &
    as.integer(evidence_unit_count) >= 2L &
    state_supported_edges >= 1L &
    cao_profile_status == "profiled_strong_candidate" &
    cao_rna_effect_class %chin% c("higher_in_tumor", "lower_in_tumor")]
driver_integrated[, integrated_tier := fcase(
  high_priority_validation_rule, "T4",
  decision == "strong_patient_level_candidate", "T3",
  decision == "conditional_candidate", "T2",
  default = "T1"
)]
driver_integrated[, integrated_decision := fcase(
  integrated_tier == "T4", "priority_for_independent_and_wetlab_validation",
  integrated_tier == "T3", "robust_source_candidate_retained",
  integrated_tier == "T2", "conditional_candidate_retained",
  default = "exploratory_retained"
)]
driver_integrated[, tier_assignment_basis := fcase(
  integrated_tier == "T4",
  paste(
    "TCGA strong driver + >=2 distinct event units + within-TCGA supported state edge",
    "+ Cao RNA direction; exact driver event still lacks independent validation"
  ),
  integrated_tier == "T3",
  "TCGA patient-level strong candidate; orthogonal layers calibrate but do not constitute exact event replication",
  integrated_tier == "T2",
  "TCGA conditional screen retained under the flexible gate; current Cao analysis covered strong candidates only",
  default = "exploratory evidence"
)]
driver_integrated[, integrated_conclusion_ceiling := fcase(
  integrated_tier == "T4",
  "高优先级独立队列/实验验证候选；不得写成已验证因果驱动轴或治疗靶点",
  integrated_tier == "T3",
  "当前公共数据范围内稳健驱动候选；无独立事件复现时不得称稳健复现",
  integrated_tier == "T2",
  "条件性驱动候选；不得进入主机制结论",
  default = "探索线索"
)]
driver_integrated[, integrated_required_validation := fcase(
  integrated_tier == "T4",
  "独立 ESCC 队列复现同一突变/CNV事件及方向；蛋白定量、细胞来源和扰动实验",
  integrated_tier == "T3",
  "独立 ESCC 患者队列复现精确事件；患者 bootstrap、跨平台与细胞来源定位",
  default = required_next_validation
)]
driver_integrated[, candidate_id := paste0("GENE:", gene_name)]
driver_integrated[, cross_cohort_patient_link_created := FALSE]
driver_integrated[, exact_independent_validation_present := FALSE]

driver_output_columns <- c(
  "candidate_id", "gene_id", "gene_name", "gene_type", "decision",
  "primary_candidate_route", "mutated_patients", "mutation_frequency_96",
  "high_confidence_mutation_pattern", "recurrent_mutation",
  "amplification_frequency", "homozygous_deletion_frequency",
  "recurrent_cnv_class", "spearman_rho_relative_cnv", "p_relative_cnv",
  "q_relative_cnv", "dosage_evidence", "discordance_code",
  "strong_dosage", "conditional_dosage", "precise_reverse_dosage",
  "multi_event_convergence", "source_evidence_unit_count",
  "mutation_event_unit", "cnv_expression_event_unit",
  "distinct_event_unit_count", "evidence_unit_count", "state_edge_count",
  "raw_state_supported_edges", "state_supported_edges",
  "state_conditional_edges", "state_exploratory_edges",
  "supported_factor_targets", "leakage_not_retained_factor_targets",
  "conditional_factor_targets",
  "conditional_pathway_targets", "strongest_state_edge",
  "strongest_state_effect", "strongest_state_q", "strongest_state_status",
  "strongest_state_raw_status", "strongest_state_leakage_countable",
  "strongest_state_leakage_failure_reason",
  "cao_profile_status", "cao_cross_layer_class", "cao_statistical_support",
  "cao_available_layer_count", "cao_fdr_supported_layer_count",
  "cao_promoter_delta_beta", "cao_rna_log_fc", "cao_rna_effect_class",
  "cao_protein_status", "cao_protein_log2_tumor_vs_normal",
  "cao_promoter_rna_relationship", "cao_rna_protein_relationship",
  "exact_driver_event_independent_validation_count",
  "independent_orthogonal_calibration_count", "quality_score_20",
  "within_layer_score_20", "cross_omics_score_25", "validation_score_20",
  "robustness_score_15", "auxiliary_score_100", "score_is_decision_gate",
  "high_priority_validation_rule", "integrated_tier", "integrated_decision",
  "tier_assignment_basis", "integrated_conclusion_ceiling",
  "integrated_required_validation", "cross_cohort_patient_link_created",
  "exact_independent_validation_present"
)
driver_candidates <- driver_integrated[, ..driver_output_columns]
driver_candidates[, tier_order_internal := unname(tier_rank[integrated_tier])]
setorder(
  driver_candidates,
  -tier_order_internal,
  -auxiliary_score_100,
  gene_name
)
driver_candidates[, tier_order_internal := NULL]

message("[3/7] 构建跨层轴边和候选轴摘要")
edge_template <- function(n) data.table(
  axis_candidate_id = rep(NA_character_, n),
  source_node = rep(NA_character_, n),
  source_layer = rep(NA_character_, n),
  target_node = rep(NA_character_, n),
  target_layer = rep(NA_character_, n),
  edge_class = rep(NA_character_, n),
  source_effect_measure = rep(NA_character_, n),
  source_effect = rep(NA_real_, n),
  target_effect_measure = rep(NA_character_, n),
  target_effect = rep(NA_real_, n),
  association_effect_measure = rep(NA_character_, n),
  association_effect = rep(NA_real_, n),
  p_value = rep(NA_real_, n),
  q_value = rep(NA_real_, n),
  n_complete = rep(NA_integer_, n),
  sample_design = rep(NA_character_, n),
  evidence_source = rep(NA_character_, n),
  independence_group = rep(NA_character_, n),
  same_patient_evidence = rep(FALSE, n),
  independent_from_tcga_discovery = rep(FALSE, n),
  countable_as_exact_driver_validation = rep(FALSE, n),
  representation_overlap = rep(NA_character_, n),
  original_state_edge_id = rep(NA_character_, n),
  leakage_gate_evaluable = rep(NA, n),
  leakage_gate_pass = rep(NA, n),
  countable_for_T3_T4 = rep(NA, n),
  leakage_gate_failure_reason = rep(NA_character_, n),
  leakage_control_decision = rep(NA_character_, n),
  baseline_reliable_match_rate = rep(NA_real_, n),
  drop_both_reliable_match_rate = rep(NA_real_, n),
  drop_both_direction_retention_rate = rep(NA_real_, n),
  drop_both_magnitude_retention_rate = rep(NA_real_, n),
  drop_both_gate_retention_rate = rep(NA_real_, n),
  drop_both_low_margin_count = rep(NA_integer_, n),
  support_status = rep(NA_character_, n),
  discordance_code = rep(NA_character_, n),
  evidence_tier = rep(NA_character_, n),
  layer_availability = rep(NA_character_, n),
  conclusion_ceiling = rep(NA_character_, n),
  required_next_validation = rep(NA_character_, n),
  source_file = rep(NA_character_, n),
  source_row_key = rep(NA_character_, n)
)

strong_driver <- driver_integrated[decision == "strong_patient_level_candidate"]
cnv_rna_edges <- edge_template(nrow(strong_driver))
cnv_rna_edges[, `:=`(
  axis_candidate_id = strong_driver$candidate_id,
  source_node = paste0(strong_driver$gene_name, ":CNV"),
  source_layer = "genome_CNV",
  target_node = paste0(strong_driver$gene_name, ":RNA"),
  target_layer = "transcriptome",
  edge_class = "TCGA_CNV_RNA_dosage",
  source_effect_measure = "relative_log2_CNV_IQR",
  source_effect = as.numeric(strong_driver$relative_log2_iqr),
  target_effect_measure = "not_separate_from_association",
  association_effect_measure = "Spearman rho",
  association_effect = as.numeric(strong_driver$spearman_rho_relative_cnv),
  p_value = as.numeric(strong_driver$p_relative_cnv),
  q_value = as.numeric(strong_driver$q_relative_cnv),
  n_complete = NA_integer_,
  sample_design = "TCGA-ESCC driver core; same-patient CNV and RNA",
  evidence_source = "TCGA_ESCC_DR45",
  independence_group = "TCGA_ESCC_DR45_driver_core_95",
  same_patient_evidence = TRUE,
  independent_from_tcga_discovery = FALSE,
  countable_as_exact_driver_validation = FALSE,
  representation_overlap = "direct_same_patient_CNV_RNA_no_external_replication",
  support_status = fcase(
    as_logical_safe(strong_driver$precise_reverse_dosage), "precise_reverse_axis_rejected",
    as_logical_safe(strong_driver$strong_dosage), "within_tcga_supported",
    as_logical_safe(strong_driver$conditional_dosage), "within_tcga_conditional",
    as_logical_safe(strong_driver$test_eligible), "directional_exploratory",
    default = "not_evaluable"
  ),
  discordance_code = fcase(
    as_logical_safe(strong_driver$precise_reverse_dosage), "comparable_precise_reverse",
    as_logical_safe(strong_driver$strong_dosage), "support",
    as_logical_safe(strong_driver$conditional_dosage),
      fifelse(is.na(strong_driver$discordance_code),
              "same_direction_nonsignificant", strong_driver$discordance_code),
    !as_logical_safe(strong_driver$test_eligible), "not_comparable",
    default = fifelse(is.na(strong_driver$discordance_code),
                      "same_direction_nonsignificant", strong_driver$discordance_code)
  ),
  evidence_tier = fcase(
    as_logical_safe(strong_driver$precise_reverse_dosage), "T0",
    as_logical_safe(strong_driver$strong_dosage), "T3",
    as_logical_safe(strong_driver$conditional_dosage), "T2",
    as_logical_safe(strong_driver$test_eligible), "T1",
    default = "T0"
  ),
  layer_availability = fifelse(
    as_logical_safe(strong_driver$test_eligible), "evaluable", "not_evaluable_not_negative"
  ),
  conclusion_ceiling = "TCGA 内部 CNV-RNA 剂量候选；不得写成独立复现或因果调控",
  required_next_validation = "独立 ESCC 队列复现 CNV-RNA 方向，并校准 focality、纯度和染色体臂背景",
  source_file = "results/tcga_escc_driver_candidate_screen.tsv",
  source_row_key = strong_driver$gene_id
)]

cao_promoter_rna <- edge_template(nrow(cao))
cao_promoter_rna[, `:=`(
  axis_candidate_id = paste0("GENE:", cao$tcga_gene_name),
  source_node = paste0(cao$tcga_gene_name, ":promoter_methylation"),
  source_layer = "epigenome_WGBS",
  target_node = paste0(cao$tcga_gene_name, ":RNA"),
  target_layer = "transcriptome",
  edge_class = "Cao2020_promoter_RNA",
  source_effect_measure = "median paired delta_beta tumor_minus_normal",
  source_effect = as.numeric(cao$promoter_median_paired_delta_beta),
  target_effect_measure = "limma logFC tumor_minus_normal pc1",
  target_effect = as.numeric(cao$limma_log_fc_pc1),
  association_effect_measure = "same-patient Spearman rho",
  association_effect = as.numeric(cao$promoter_rna_rho),
  p_value = as.numeric(cao$promoter_rna_p),
  q_value = NA_real_,
  n_complete = as.integer(cao$promoter_rna_n),
  sample_design = "Cao 2020; WGBS 9 pairs and RNA 10 pairs; overlapping patients only for correlation",
  evidence_source = "Cao2020_GSE149608_GSE149609",
  independence_group = "Cao2020_same_patient_cross_layer",
  same_patient_evidence = TRUE,
  independent_from_tcga_discovery = TRUE,
  countable_as_exact_driver_validation = FALSE,
  representation_overlap = "independent_cohort_orthogonal_calibration_not_exact_driver_event_replication",
  support_status = fcase(
    cao$promoter_region_evaluation_status == "technically_unavailable", "not_evaluable",
    cao$promoter_rna_relationship == "directionally_coherent_inverse", "directional_conditional",
    grepl("credible_reverse", cao$promoter_rna_relationship), "context_specific_reverse_retained",
    cao$promoter_rna_relationship %chin% c("weak_or_mixed", "rna_pseudocount_sensitive"),
      "directional_exploratory",
    default = "not_evaluable"
  ),
  discordance_code = fcase(
    cao$promoter_region_evaluation_status == "technically_unavailable", "not_comparable",
    cao$promoter_rna_relationship == "directionally_coherent_inverse", "support",
    grepl("credible_reverse", cao$promoter_rna_relationship), "context_specific_reverse",
    default = "same_direction_nonsignificant"
  ),
  evidence_tier = fcase(
    cao$promoter_region_evaluation_status == "technically_unavailable", "T0",
    cao$promoter_rna_relationship == "directionally_coherent_inverse", "T2",
    grepl("credible_reverse", cao$promoter_rna_relationship), "T2",
    default = "T1"
  ),
  layer_availability = fifelse(
    cao$promoter_region_evaluation_status == "technically_unavailable",
    "methylation_technically_unavailable_not_negative", "evaluable"
  ),
  conclusion_ceiling = cao$conclusion_ceiling,
  required_next_validation = cao$required_follow_up,
  source_file = "results/cao2020_cross_layer_candidate_summary.tsv",
  source_row_key = cao$tcga_gene_id
)]

cao_rna_protein <- edge_template(nrow(cao))
cao_rna_protein[, `:=`(
  axis_candidate_id = paste0("GENE:", cao$tcga_gene_name),
  source_node = paste0(cao$tcga_gene_name, ":RNA"),
  source_layer = "transcriptome",
  target_node = paste0(cao$tcga_gene_name, ":protein"),
  target_layer = "proteome",
  edge_class = "Cao2020_RNA_protein",
  source_effect_measure = "limma logFC tumor_minus_normal pc1",
  source_effect = as.numeric(cao$limma_log_fc_pc1),
  target_effect_measure = "median patient log2 tumor_vs_normal",
  target_effect = as.numeric(cao$median_patient_log2_tumor_vs_normal),
  association_effect_measure = "same-patient Spearman rho",
  association_effect = as.numeric(cao$rna_protein_rho),
  p_value = as.numeric(cao$rna_protein_p),
  q_value = NA_real_,
  n_complete = as.integer(cao$rna_protein_n),
  sample_design = "Cao 2020 RNA and proteome patient mapping; missing protein is not negative",
  evidence_source = "Cao2020_GSE149609_proteome_supplement",
  independence_group = "Cao2020_same_patient_cross_layer",
  same_patient_evidence = TRUE,
  independent_from_tcga_discovery = TRUE,
  countable_as_exact_driver_validation = FALSE,
  representation_overlap = "independent_cohort_orthogonal_calibration_not_exact_driver_event_replication",
  support_status = fcase(
    cao$protein_measurement_status == "quantified" &
      cao$rna_protein_relationship == "directionally_coherent", "directional_conditional",
    cao$protein_measurement_status == "quantified" &
      grepl("credible_reverse", cao$rna_protein_relationship),
      "context_specific_reverse_retained",
    cao$protein_measurement_status == "quantified", "directional_exploratory",
    default = "layer_missing_not_negative"
  ),
  discordance_code = fcase(
    cao$protein_measurement_status != "quantified", "not_comparable",
    cao$rna_protein_relationship == "directionally_coherent", "support",
    grepl("credible_reverse", cao$rna_protein_relationship), "context_specific_reverse",
    default = "same_direction_nonsignificant"
  ),
  evidence_tier = fcase(
    cao$protein_measurement_status != "quantified", "T0",
    cao$rna_protein_relationship == "directionally_coherent", "T2",
    grepl("credible_reverse", cao$rna_protein_relationship), "T2",
    default = "T1"
  ),
  layer_availability = fcase(
    cao$protein_measurement_status == "quantified", "evaluable",
    cao$protein_measurement_status == "identified_not_quantified",
      "identified_not_quantified_not_negative",
    default = "not_identified_not_negative"
  ),
  conclusion_ceiling = cao$conclusion_ceiling,
  required_next_validation = cao$required_follow_up,
  source_file = "results/cao2020_cross_layer_candidate_summary.tsv",
  source_row_key = cao$tcga_gene_id
)]

tcga_state_axis <- edge_template(nrow(state_edges))
tcga_state_axis[, `:=`(
  axis_candidate_id = paste0("GENE:", state_edges$source_gene_name),
  source_node = state_edges$source_node,
  source_layer = fifelse(
    state_edges$event_type == "mutation", "genome_mutation", "genome_CNV"
  ),
  target_node = state_edges$target_node,
  target_layer = fifelse(
    state_edges$target_type == "MOFA_factor", "multiomics_state",
    "RNA_derived_pathway_activity"
  ),
  edge_class = "TCGA_driver_state_internal",
  association_effect_measure = state_edges$effect_measure,
  association_effect = as.numeric(state_edges$effect),
  p_value = as.numeric(state_edges$p_value),
  q_value = as.numeric(state_edges$q_value),
  n_complete = NA_integer_,
  sample_design = "TCGA-ESCC five-layer core; same patients and partially shared model inputs",
  evidence_source = "TCGA_ESCC_DR45",
  independence_group = state_edges$independence_group,
  same_patient_evidence = TRUE,
  independent_from_tcga_discovery = FALSE,
  countable_as_exact_driver_validation = FALSE,
  representation_overlap = state_edges$target_model_overlap_class,
  original_state_edge_id = state_edges$edge_id,
  leakage_gate_evaluable = as_logical_safe(state_edges$leakage_gate_evaluable),
  leakage_gate_pass = as_logical_safe(state_edges$leakage_gate_pass),
  countable_for_T3_T4 = as_logical_safe(state_edges$countable_for_T3_T4),
  leakage_gate_failure_reason = state_edges$leakage_gate_failure_reason,
  leakage_control_decision = state_edges$leakage_control_decision,
  baseline_reliable_match_rate =
    as.numeric(state_edges$baseline_reliable_match_rate),
  drop_both_reliable_match_rate =
    as.numeric(state_edges$drop_both_reliable_match_rate),
  drop_both_direction_retention_rate =
    as.numeric(state_edges$drop_both_direction_retention_rate),
  drop_both_magnitude_retention_rate =
    as.numeric(state_edges$drop_both_magnitude_retention_rate),
  drop_both_gate_retention_rate =
    as.numeric(state_edges$drop_both_gate_retention_rate),
  drop_both_low_margin_count =
    as.integer(state_edges$drop_both_low_margin_count),
  support_status = state_edges$integration_state_status,
  discordance_code = fifelse(
    state_edges$leakage_controlled_support,
    "support", "same_direction_nonsignificant"
  ),
  evidence_tier = state_edges$integration_state_tier,
  layer_availability = fifelse(
    as_logical_safe(state_edges$level_factor_soft_flag),
    "evaluable_but_level_factor_soft_qc", "evaluable"
  ),
  conclusion_ceiling = fifelse(
    state_edges$target_type == "MOFA_factor" &
      !is.na(state_edges$leakage_conclusion_ceiling),
    state_edges$leakage_conclusion_ceiling,
    state_edges$conclusion_ceiling
  ),
  required_next_validation = fifelse(
    state_edges$target_type == "MOFA_factor" &
      !is.na(state_edges$leakage_required_next_validation),
    state_edges$leakage_required_next_validation,
    state_edges$required_next_validation
  ),
  source_file = fifelse(
    state_edges$target_type == "MOFA_factor",
    "results/tcga_escc_driver_state_leakage_candidate_summary.tsv",
    "results/tcga_escc_driver_state_network_edges.tsv"
  ),
  source_row_key = state_edges$edge_id
)]

factor_level_flag <- factor_level_qc[, .(
  level_factor_soft_flag = any(level_factor_soft_flag)
), by = factor]
fail_if(nrow(factor_level_flag) != 8L ||
          !setequal(factor_level_flag$factor, expected_factors) ||
          factor_level_flag[level_factor_soft_flag == TRUE, .N] != 1L ||
          factor_level_flag[level_factor_soft_flag == TRUE, factor] != "Factor4",
        "factor level 软门禁必须唯一标记 Factor4。")

ecms_factor_context <- merge(
  ecms_factor_primary, factor_level_flag,
  by = "factor", all.x = TRUE, sort = FALSE
)
fail_if(anyNA(ecms_factor_context$level_factor_soft_flag),
        "ECMS factor 无法一对一映射 level-factor QC。")
ecms_factor_context[, evidence_tier := assign_ecms_factor_tier(
  level_factor_soft_flag,
  anova_p_value,
  anova_q_value,
  kruskal_p_value,
  kruskal_q_value,
  eta_squared
)]
ecms_factor_context[, support_status := fcase(
  level_factor_soft_flag, "level_factor_background",
  evidence_tier == "T2", "ecms_context_explanatory_supported",
  default = "ecms_context_exploratory"
)]

ecms_factor_axis <- edge_template(nrow(ecms_factor_context))
ecms_factor_axis[, `:=`(
  axis_candidate_id = paste0("STATE:", ecms_factor_context$factor),
  source_node = "ECMS1-4 locked expression benchmark",
  source_layer = "external_expression_classifier_state",
  target_node = ecms_factor_context$factor,
  target_layer = "multiomics_state",
  edge_class = "TCGA_ECMS_factor_explanatory",
  association_effect_measure = "eta_squared categorical separation",
  association_effect = as.numeric(ecms_factor_context$eta_squared),
  p_value = pmin(
    ecms_factor_context$anova_p_value,
    ecms_factor_context$kruskal_p_value,
    na.rm = TRUE
  ),
  q_value = pmin(
    ecms_factor_context$anova_q_value,
    ecms_factor_context$kruskal_q_value,
    na.rm = TRUE
  ),
  n_complete = as.integer(ecms_factor_context$n_patients),
  sample_design = paste(
    "TCGA-ESCC official repository-anchor 78 patients;",
    "ECMS and MOFA share TCGA RNA representation"
  ),
  evidence_source = "TCGA_ESCC_DR45",
  independence_group = "TCGA_ESCC_ECMS_official78_shared_RNA",
  same_patient_evidence = TRUE,
  independent_from_tcga_discovery = FALSE,
  countable_as_exact_driver_validation = FALSE,
  representation_overlap = "ECMS_RNA_derived_and_MOFA_contains_RNA",
  support_status = ecms_factor_context$support_status,
  discordance_code = fcase(
    ecms_factor_context$level_factor_soft_flag,
    "technical_background_level_factor",
    ecms_factor_context$evidence_tier == "T2", "support",
    default = "same_direction_nonsignificant"
  ),
  evidence_tier = ecms_factor_context$evidence_tier,
  layer_availability = fifelse(
    ecms_factor_context$level_factor_soft_flag,
    "evaluable_but_level_factor_soft_qc", "evaluable"
  ),
  conclusion_ceiling = fcase(
    ecms_factor_context$evidence_tier == "T0",
    "官方78例同 TCGA RNA 共享表示的 T0 技术/背景 ECMS 语境；不作生物学轴",
    ecms_factor_context$evidence_tier == "T1",
    "官方78例同 TCGA RNA 共享表示的 T1 探索性 ECMS 语境；该证据类别最高 T2",
    ecms_factor_context$evidence_tier == "T2",
    "官方78例同 TCGA RNA 共享表示的 T2 解释性 ECMS 语境；不是新亚型或独立验证",
    default = "ECMS 行级证据层级异常"
  ),
  required_next_validation =
    "独立 ESCC 队列投影；去 RNA 表示敏感性；细胞来源和临床协变量校准",
  source_file = "results/tcga_escc_ecms_factor_associations.tsv",
  source_row_key = ecms_factor_context$factor
)]

ecms_adjusted_context <- merge(
  ecms_adjusted_primary, factor_level_flag,
  by = "factor", all.x = TRUE, sort = FALSE
)
fail_if(anyNA(ecms_adjusted_context$level_factor_soft_flag),
        "ECMS-adjusted Factor-PROGENy 无法映射 level-factor QC。")
ecms_adjusted_context[, evidence_tier := assign_ecms_adjusted_tier(
  level_factor_soft_flag,
  incremental_p_value,
  incremental_q_value,
  partial_r_squared
)]
ecms_adjusted_context[, support_status := fcase(
  level_factor_soft_flag, "level_factor_background",
  evidence_tier == "T2", "ecms_adjusted_explanatory_supported",
  default = "ecms_adjusted_directional_exploratory"
)]

ecms_adjusted_axis <- edge_template(nrow(ecms_adjusted_context))
ecms_adjusted_axis[, `:=`(
  axis_candidate_id = paste0("STATE:", ecms_adjusted_context$factor),
  source_node = ecms_adjusted_context$factor,
  source_layer = "multiomics_state",
  target_node = ecms_adjusted_context$pathway,
  target_layer = "RNA_derived_pathway_activity",
  edge_class = "TCGA_ECMS_adjusted_factor_PROGENy",
  source_effect_measure = "delta R2 after ECMS categorical baseline",
  source_effect = as.numeric(ecms_adjusted_context$delta_r_squared),
  target_effect_measure = "partial R2",
  target_effect = as.numeric(ecms_adjusted_context$partial_r_squared),
  association_effect_measure = "standardized factor beta after ECMS adjustment",
  association_effect = as.numeric(ecms_adjusted_context$factor_beta_standardized),
  p_value = as.numeric(ecms_adjusted_context$incremental_p_value),
  q_value = as.numeric(ecms_adjusted_context$incremental_q_value),
  n_complete = as.integer(ecms_adjusted_context$n_patients),
  sample_design = paste(
    "TCGA-ESCC official repository-anchor 78 patients;",
    "standardized PROGENy ~ ECMS + standardized MOFA factor"
  ),
  evidence_source = "TCGA_ESCC_DR45",
  independence_group = "TCGA_ESCC_ECMS_official78_shared_RNA",
  same_patient_evidence = TRUE,
  independent_from_tcga_discovery = FALSE,
  countable_as_exact_driver_validation = FALSE,
  representation_overlap =
    "ECMS_and_PROGENy_RNA_derived_MOFA_contains_RNA",
  support_status = ecms_adjusted_context$support_status,
  discordance_code = fcase(
    ecms_adjusted_context$level_factor_soft_flag,
    "technical_background_level_factor",
    ecms_adjusted_context$evidence_tier == "T2", "support",
    default = "same_direction_nonsignificant"
  ),
  evidence_tier = ecms_adjusted_context$evidence_tier,
  layer_availability = fifelse(
    ecms_adjusted_context$level_factor_soft_flag,
    "evaluable_but_level_factor_soft_qc", "evaluable"
  ),
  conclusion_ceiling = fcase(
    ecms_adjusted_context$evidence_tier == "T0",
    "ECMS 调整后的同 TCGA、共享 RNA 表示 T0 技术/背景边；不作生物学通路轴",
    ecms_adjusted_context$evidence_tier == "T1",
    "ECMS 调整后的同 TCGA、共享 RNA 表示 T1 探索边；该证据类别最高 T2",
    ecms_adjusted_context$evidence_tier == "T2",
    "ECMS 调整后的同 TCGA、共享 RNA 表示 T2 解释边；不是独立通路验证或因果链",
    default = "ECMS-adjusted 行级证据层级异常"
  ),
  required_next_validation =
    "去共享 RNA 表示敏感性、独立 ESCC 队列投影和通路活性功能验证",
  source_file =
    "results/tcga_escc_ecms_adjusted_factor_progeny_associations.tsv",
  source_row_key = paste(
    ecms_adjusted_context$factor,
    ecms_adjusted_context$pathway,
    sep = "::"
  )
)]

axis_edges <- rbindlist(list(
  cnv_rna_edges, cao_promoter_rna, cao_rna_protein,
  tcga_state_axis, ecms_factor_axis, ecms_adjusted_axis
), use.names = TRUE, fill = TRUE)
axis_edges[, edge_id := sprintf("ESCC_AXIS_%05d", seq_len(.N))]
setcolorder(axis_edges, c("edge_id", setdiff(names(axis_edges), "edge_id")))
axis_edges[, tier_order_internal := unname(tier_rank[evidence_tier])]
setorder(
  axis_edges,
  -tier_order_internal,
  axis_candidate_id,
  edge_class,
  q_value,
  p_value,
  na.last = TRUE
)
axis_edges[, tier_order_internal := NULL]

# 为每个 strong driver 形成一行候选轴裁决；不会把两跳路径写成直接边。
# 通路语境改用官方78例中 ECMS 调整后的 Factor-PROGENy T2 边；
# hybrid94 及未调整的 94 例相关不进入主组合路径。
ecms_adjusted_pathway_lookup <- ecms_adjusted_context[
  evidence_tier == "T2"
][order(incremental_q_value, -partial_r_squared), .(
  ecms_adjusted_pathways = collapse_unique(paste0(
    pathway, "(beta=", sprintf("%.2f", factor_beta_standardized),
    ",partialR2=", sprintf("%.2f", partial_r_squared), ")"
  )),
  ecms_adjusted_pathway_count = .N
), by = factor]
ecms_factor_context_lookup <- ecms_factor_context[, .(
  factor,
  ecms_context_tier = evidence_tier,
  ecms_context_status = support_status,
  ecms_eta_squared = eta_squared,
  ecms_context_q = pmin(anova_q_value, kruskal_q_value, na.rm = TRUE)
)]
state_factor_for_composition <- state_edges[
  target_type == "MOFA_factor" &
    integration_state_tier %chin% c("T3", "T2")
]
composed_rows <- merge(
  state_factor_for_composition[, .(
    original_state_edge_id = edge_id,
    gene_name = source_gene_name,
    event_type,
    factor = target_node,
    driver_factor_status = integration_state_status,
    driver_factor_tier = integration_state_tier
  )],
  ecms_factor_context_lookup,
  by = "factor",
  all.x = TRUE,
  sort = FALSE
)
composed_rows <- merge(
  composed_rows,
  ecms_adjusted_pathway_lookup,
  by = "factor",
  all.x = TRUE,
  sort = FALSE
)
composed_rows <- composed_rows[
  ecms_context_tier == "T2" &
    !is.na(ecms_adjusted_pathways) & nzchar(ecms_adjusted_pathways)
]
composed_summary <- composed_rows[, .(
  composed_state_pathway_hypotheses = collapse_unique(paste0(
    event_type, "->", factor, "->", ecms_adjusted_pathways,
    "[", driver_factor_status, ";", ecms_context_status, "]"
  )),
  ecms_contextualized_leakage_edge_count = uniqueN(
    original_state_edge_id[driver_factor_tier == "T3"]
  ),
  ecms_contextualized_factor_targets = collapse_unique(
    factor[driver_factor_tier == "T3"]
  ),
  ecms_adjusted_pathway_targets = collapse_unique(ecms_adjusted_pathways)
), by = gene_name]

axis_summary <- driver_integrated[
  decision == "strong_patient_level_candidate",
  .(
    candidate_id,
    gene_id,
    gene_name,
    tcga_driver_tier = integrated_tier,
    tcga_primary_route = primary_candidate_route,
    distinct_event_units = evidence_unit_count,
    cnv_rna_dosage_status = dosage_evidence,
    cnv_rna_rho = spearman_rho_relative_cnv,
    cnv_rna_q = q_relative_cnv,
    state_supported_edges,
    state_conditional_edges,
    supported_factor_targets,
    conditional_factor_targets,
    conditional_pathway_targets,
    cao_cross_layer_class,
    cao_statistical_support,
    cao_promoter_delta_beta,
    cao_rna_log_fc,
    cao_rna_effect_class,
    cao_protein_status,
    cao_protein_log2_tumor_vs_normal,
    cao_promoter_rna_relationship,
    cao_rna_protein_relationship,
    exact_driver_event_independent_validation_count
  )
]
axis_summary <- merge(axis_summary, composed_summary, by = "gene_name", all.x = TRUE)
axis_summary[is.na(composed_state_pathway_hypotheses),
             composed_state_pathway_hypotheses := ""]
axis_summary[is.na(ecms_contextualized_leakage_edge_count),
             ecms_contextualized_leakage_edge_count := 0L]
axis_summary[is.na(ecms_contextualized_factor_targets),
             ecms_contextualized_factor_targets := ""]
axis_summary[is.na(ecms_adjusted_pathway_targets),
             ecms_adjusted_pathway_targets := ""]
axis_summary[, `:=`(
  ecms_primary_patient_count = ecms_primary_patient_count,
  ecms_extension_calibration_pass = ecms_extension_calibration_pass,
  ecms_shared_rna_representation = TRUE,
  ecms_or_progeny_independent_validation = FALSE
)]
axis_summary[, composed_path_is_direct_evidence := FALSE]
axis_summary[, axis_class := fcase(
  cao_cross_layer_class == "credible_reverse_retained",
    "context_specific_reverse_axis_retained",
  as.integer(distinct_event_units) >= 2L & state_supported_edges >= 1L &
    cao_rna_effect_class %chin% c("higher_in_tumor", "lower_in_tumor") &
    ecms_contextualized_leakage_edge_count >= 1L,
    "driver_expression_ecms_state_bridge_priority",
  cao_cross_layer_class == "three_layer_directional_hypothesis",
    "three_layer_directional_axis_conditional",
  state_supported_edges >= 1L & ecms_contextualized_leakage_edge_count >= 1L,
    "driver_ecms_state_internal_bridge",
  state_supported_edges >= 1L, "driver_state_internal_bridge",
  state_conditional_edges >= 1L, "driver_state_conditional_bridge",
  default = "driver_candidate_axis_incomplete"
)]
axis_summary[, axis_tier := assign_axis_summary_tier(
  as.integer(distinct_event_units),
  as.integer(state_supported_edges),
  cao_rna_effect_class %chin% c("higher_in_tumor", "lower_in_tumor"),
  as.integer(ecms_contextualized_leakage_edge_count),
  as.integer(state_conditional_edges),
  cao_cross_layer_class %chin% c(
    "three_layer_directional_hypothesis", "credible_reverse_retained"
  ),
  cao_cross_layer_class == "credible_reverse_retained"
)]
axis_summary[, exception_or_boundary := fcase(
  axis_class == "context_specific_reverse_axis_retained",
    "情境特异反向保留；需复核组织/阶段/细胞组成，不以其他分数抵消",
  axis_class == "three_layer_directional_axis_conditional",
    "Cao 小样本三层方向性、候选层 q>0.05；条件保留",
  cao_protein_status != "quantified",
    "蛋白未鉴定或未定量，不计为阴性",
  default = "同 TCGA 状态边非独立验证"
)]
axis_summary[, conclusion_ceiling := fcase(
  axis_tier == "T4", "高优先级可证伪验证轴；不是已证实因果网络",
  axis_tier == "T3", "TCGA 来源内稳健 driver-state 候选轴；不是独立复现",
  axis_tier == "T2", "条件性、方向性或情境特异候选轴",
  default = "轴不完整，仅保留 driver 背景"
)]
axis_summary[, required_next_validation := fcase(
  axis_class == "driver_expression_ecms_state_bridge_priority",
    "独立 ESCC 复现精确 driver 事件→RNA/状态；补蛋白定量和扰动实验",
  axis_class == "three_layer_directional_axis_conditional",
    "扩大配对 WGBS-RNA-蛋白队列并做患者 bootstrap；独立复现方向",
  axis_class == "context_specific_reverse_axis_retained",
    "独立配对队列复核反向；以蛋白、细胞来源和功能实验区分脱耦",
  state_supported_edges + state_conditional_edges > 0L,
    "去除 predictor 后重训因子、跨 seed 对齐并投影独立 ESCC 队列",
  default = "补足可评估跨层读出"
)]

# driver candidate/ledger 的最终层级必须与最终集成轴严格同步。
# “patient-level strong screen”只是入口决策；若没有可计数的
# leakage-controlled state edge，不得因为初筛标签而在候选表中维持 T3。
# 这也使 EP300/FBXW7/NOTCH3 的条件轴、GNAS 的方向性轴以及
# ZNF750 的反向边界在结构化表与正文中保持一致。
strong_driver_candidate_index <- which(
  driver_candidates$decision == "strong_patient_level_candidate"
)
driver_axis_match <- match(
  driver_candidates$gene_name[strong_driver_candidate_index],
  axis_summary$gene_name
)
fail_if(anyNA(driver_axis_match),
        "strong driver candidate 缺少最终集成轴裁决，无法同步层级。")
driver_candidates[strong_driver_candidate_index,
  high_priority_validation_rule :=
    axis_summary$axis_tier[driver_axis_match] == "T4"]
driver_candidates[strong_driver_candidate_index,
  integrated_tier := axis_summary$axis_tier[driver_axis_match]]
driver_candidates[strong_driver_candidate_index,
  integrated_decision := fcase(
  integrated_tier == "T4", "priority_for_independent_and_wetlab_validation",
  integrated_tier == "T3", "robust_source_candidate_retained",
  axis_summary$axis_class[driver_axis_match] ==
    "three_layer_directional_axis_conditional",
    "directional_hypothesis_conditional",
  axis_summary$axis_class[driver_axis_match] ==
    "context_specific_reverse_axis_retained",
    "context_specific_reverse_retained",
  integrated_tier == "T2", "conditional_candidate_retained",
  default = "exploratory_retained"
)]
driver_candidates[strong_driver_candidate_index,
  tier_assignment_basis := fcase(
  integrated_tier == "T4",
  paste(
    "TCGA strong driver + >=2 event units + leakage-controlled state edge +",
    "Cao RNA direction + official78 ECMS context/adjusted pathway"
  ),
  integrated_tier == "T3",
  paste(
    "TCGA patient-level strong candidate + leakage-controlled state edge;",
    "orthogonal layers calibrate but do not constitute exact replication"
  ),
  axis_summary$axis_class[driver_axis_match] ==
    "three_layer_directional_axis_conditional",
    "Cao same-patient directional hypothesis without candidate-layer FDR support",
  axis_summary$axis_class[driver_axis_match] ==
    "context_specific_reverse_axis_retained",
    "context-specific reverse cross-layer boundary retained without score compensation",
  integrated_tier == "T2",
    "conditional state evidence without a leakage-controlled supported edge",
  default = "exploratory evidence"
)]
driver_candidates[strong_driver_candidate_index,
  integrated_conclusion_ceiling :=
  axis_summary$conclusion_ceiling[driver_axis_match]]
driver_candidates[strong_driver_candidate_index,
  integrated_required_validation :=
  axis_summary$required_next_validation[driver_axis_match]]
driver_candidates[, tier_order_internal := unname(tier_rank[integrated_tier])]
setorder(
  driver_candidates,
  -tier_order_internal,
  -auxiliary_score_100,
  gene_name
)
driver_candidates[, tier_order_internal := NULL]
axis_summary[, tcga_driver_tier := driver_candidates$integrated_tier[
  match(gene_name, driver_candidates$gene_name)
]]
fail_if(anyNA(axis_summary$tcga_driver_tier),
        "axis summary 无法回填最终 driver tier。")
axis_summary[, cross_cohort_patient_link_created := FALSE]
axis_summary[, exact_independent_validation_present := FALSE]
axis_summary[, tier_order_internal := unname(tier_rank[axis_tier])]
setorder(axis_summary, -tier_order_internal, gene_name)
axis_summary[, tier_order_internal := NULL]

message("[4/7] 汇总连续异质性轴，不命名未通过门禁的离散亚型")
variance_long <- factor_variance[, .(
  factor = `r2_per_factor.factor`,
  view = `r2_per_factor.view`,
  r2_percent = as.numeric(`r2_per_factor.value`)
)]
variance_profile <- variance_long[, .(
  view_r2_sum = sum(r2_percent, na.rm = TRUE),
  dominant_view = view[which.max(r2_percent)],
  dominant_view_r2_percent = max(r2_percent, na.rm = TRUE),
  view_r2_profile = collapse_unique(paste0(view, "=", sprintf("%.2f", r2_percent)))
), by = factor]
raw_pathway_factor_summary <- factor_pathways[
  order(q_value, -abs_spearman_rho), .(
  raw_top_pathway = pathway[[1L]],
  raw_top_pathway_rho = spearman_rho[[1L]],
  raw_top_pathway_q = q_value[[1L]],
  raw_q10_pathway_count = sum(q_value <= 0.10, na.rm = TRUE),
  raw_q10_pathways = collapse_unique(
    paste0(pathway[q_value <= 0.10], "(",
           sprintf("%.2f", spearman_rho[q_value <= 0.10]), ")")
  )
), by = factor]
ecms_adjusted_factor_summary <- ecms_adjusted_context[
  order(-as.integer(evidence_tier == "T2"),
        incremental_q_value, -partial_r_squared), .(
  top_pathway = pathway[[1L]],
  top_pathway_beta = factor_beta_standardized[[1L]],
  top_pathway_partial_r_squared = partial_r_squared[[1L]],
  top_pathway_p = incremental_p_value[[1L]],
  top_pathway_q = incremental_q_value[[1L]],
  ecms_adjusted_pathway_t2_count = sum(evidence_tier == "T2"),
  ecms_adjusted_pathways_t2 = collapse_unique(paste0(
    pathway[evidence_tier == "T2"], "(beta=",
    sprintf("%.2f", factor_beta_standardized[evidence_tier == "T2"]),
    ",partialR2=",
    sprintf("%.2f", partial_r_squared[evidence_tier == "T2"]), ")"
  ))
), by = factor]
ecms_factor_heterogeneity <- ecms_factor_context[, .(
  factor,
  ecms_context_tier = evidence_tier,
  ecms_context_status = support_status,
  ecms_eta_squared = eta_squared,
  ecms_context_q = pmin(anova_q_value, kruskal_q_value, na.rm = TRUE)
)]
weight_summary <- factor_weights[order(factor, -abs_weight), .(
  top_weight_features = collapse_unique(head(
    paste0(view, ":", feature, "(", sprintf("%.3f", value), ")"), 10L
  ))
), by = factor]
level_summary <- factor_level_qc[, .(
  level_factor_soft_flag = any(level_factor_soft_flag),
  strongest_view_mean_abs_rho = max(abs(spearman_rho_with_view_mean), na.rm = TRUE),
  strongest_view_mean_q = min_finite(q_value)
), by = factor]
driver_factor_summary <- state_edges[target_type == "MOFA_factor", .(
  raw_supported_driver_edge_count = sum(
    association_status == "within_tcga_supported"
  ),
  supported_driver_edge_count = sum(leakage_controlled_support),
  conditional_driver_edge_count = sum(integration_state_tier == "T2"),
  linked_driver_events = collapse_unique(paste0(
    source_node, "(", integration_state_status, ")"
  ))
), by = .(factor = target_node)]

factors <- data.table(factor = sort(unique(factor_pathways$factor)))
heterogeneity_axes <- Reduce(
  function(x, y) merge(x, y, by = "factor", all.x = TRUE, sort = FALSE),
  list(
    factors, variance_profile, raw_pathway_factor_summary,
    ecms_adjusted_factor_summary, ecms_factor_heterogeneity,
    weight_summary, level_summary, driver_factor_summary
  )
)
heterogeneity_axes[is.na(supported_driver_edge_count), supported_driver_edge_count := 0L]
heterogeneity_axes[is.na(raw_supported_driver_edge_count),
                   raw_supported_driver_edge_count := 0L]
heterogeneity_axes[is.na(conditional_driver_edge_count), conditional_driver_edge_count := 0L]
heterogeneity_axes[is.na(raw_q10_pathway_count), raw_q10_pathway_count := 0L]
heterogeneity_axes[is.na(ecms_adjusted_pathway_t2_count),
                   ecms_adjusted_pathway_t2_count := 0L]
heterogeneity_axes[, ecms_context_t2 := ecms_context_tier == "T2"]

selected_cluster <- cluster_evaluation[as_logical_safe(selected_for_description)]
fail_if(nrow(selected_cluster) != 1L, "描述性聚类选择必须唯一。")
heterogeneity_axes[, `:=`(
  patient_count = 94L,
  ecms_primary_patient_count = ecms_primary_patient_count,
  ecms_extension_calibration_pass = ecms_extension_calibration_pass,
  ecms_shared_rna_representation = TRUE,
  ecms_progeny_independent_validation = FALSE,
  discrete_subtype_support = FALSE,
  descriptive_cluster_k = selected_cluster$k[[1L]],
  descriptive_snf_silhouette = selected_cluster$snf_mean_silhouette[[1L]],
  descriptive_mofa_pac = selected_cluster$mofa_consensus_pac[[1L]],
  descriptive_snf_mofa_ari = selected_cluster$snf_mofa_adjusted_rand[[1L]],
  cluster_decision = "no_stable_discrete_clusters_continuous_factors_retained"
)]
heterogeneity_axes[, evidence_tier := assign_heterogeneity_tier(
  level_factor_soft_flag,
  as.integer(supported_driver_edge_count),
  ecms_context_t2,
  as.integer(ecms_adjusted_pathway_t2_count),
  as.integer(conditional_driver_edge_count)
)]
heterogeneity_axes[, t3_support_source := fifelse(
  evidence_tier == "T3",
  "leakage_controlled_driver_edge_plus_ECMS_context_and_adjusted_pathway",
  "not_applicable_no_composite_T3"
)]
heterogeneity_axes[, interpretation := fcase(
  level_factor_soft_flag,
    "组学均值 level-factor 软标记；保留解释但排除出聚类因子空间",
  ecms_adjusted_pathway_t2_count > 0L,
    paste0("连续状态轴；官方78例 ECMS 调整后优先通路：",
           ecms_adjusted_pathways_t2),
  raw_q10_pathway_count > 0L,
    paste0("连续状态轴；仅保留未调整通路敏感性：", raw_q10_pathways),
  default = "连续状态轴，当前缺少明确 ECMS-adjusted 通路解释"
)]
heterogeneity_axes[, conclusion_ceiling := fifelse(
  level_factor_soft_flag,
  "仅作水平因子/QC 背景，不用于亚型命名",
  paste(
    "TCGA 内部连续异质性假设；ECMS/PROGENy 共享 RNA 证据单独最高 T2；",
    "不得命名为稳定离散亚型或因果状态"
  )
)]
heterogeneity_axes[, required_next_validation := fifelse(
  level_factor_soft_flag,
  "检查 HM450 水平、纯度和技术协变量；独立模型不得以该因子驱动聚类",
  "跨 seed 因子对齐、独立队列投影、细胞来源和临床协变量校准"
)]
heterogeneity_axes[, tier_order_internal := unname(tier_rank[evidence_tier])]
setorder(heterogeneity_axes, -tier_order_internal, factor)
heterogeneity_axes[, tier_order_internal := NULL]

message("[5/7] 汇总代谢与可选微生物正交模块")
result_basenames <- basename(dir_ls(results_dir, type = "file", fail = FALSE))
microbe_any <- any(grepl("^prjna766558_dada2_", result_basenames))
microbe_available <- microbe_any
prjna_artifact_manifest <- data.table()

microbe_genus <- NULL
microbe_alpha <- NULL
microbe_beta <- NULL
microbe_pipeline <- NULL
microbe_ancombc2 <- NULL
if (microbe_available) {
  prjna_artifact_manifest <- verify_exact_artifact_family(
    prjna_manifest_path,
    prjna_formal_names,
    "scripts/18_analyze_prjna766558_dada2.R",
    file.path(project_root, "scripts", "18_analyze_prjna766558_dada2.R"),
    "executed_script_sha256",
    "^prjna766558_dada2_.*"
  )
  microbe_sample_qc <- read_tsv(file.path(
    results_dir, "prjna766558_dada2_sample_qc.tsv"
  ))
  microbe_alpha <- read_tsv(file.path(
    results_dir, "prjna766558_dada2_alpha_paired_tests.tsv"
  ))
  microbe_beta <- read_tsv(file.path(
    results_dir, "prjna766558_dada2_beta_tests.tsv"
  ))
  microbe_genus <- read_tsv(file.path(
    results_dir, "prjna766558_dada2_genus_paired_differential.tsv"
  ))
  microbe_pipeline <- read_tsv(file.path(
    results_dir, "prjna766558_dada2_pipeline_sensitivity.tsv"
  ))
  microbe_ancombc2 <- read_tsv(file.path(
    results_dir, "prjna766558_dada2_ancombc2_sensitivity.tsv"
  ))
  require_columns(microbe_sample_qc, c(
    "sample_name", "run_accession", "patient_pair_id", "tissue_role",
    "input_reads", "prokaryotic_reads_min_overlap12", "blank_control_status",
    "public_negative_control_run"
  ), "prjna766558_dada2_sample_qc.tsv")
  require_columns(microbe_alpha, c(
    "metric", "n_pairs", "median_paired_difference_tumor_minus_normal",
    "paired_wilcoxon_p", "paired_wilcoxon_q",
    "bootstrap_direction_consistency", "interpretation_scope"
  ), "prjna766558_dada2_alpha_paired_tests.tsv")
  require_columns(microbe_beta, c(
    "distance", "n_pairs", "permanova_r2", "restricted_permutation_p",
    "restricted_permutation_q", "betadisper_p"
  ), "prjna766558_dada2_beta_tests.tsv")
  require_columns(microbe_genus, c(
    "genus_label", "n_pairs", "median_clr_difference_tumor_minus_normal",
    "paired_wilcoxon_p", "paired_wilcoxon_q", "effect_direction",
    "candidate_tier", "pipeline_direction_stable",
    "prevalence_threshold_direction_stable", "bootstrap_direction_consistency",
    "blank_control_status", "contamination_risk", "independence_group",
    "host_sample_level_link_allowed", "conclusion_ceiling"
  ), "prjna766558_dada2_genus_paired_differential.tsv")
  require_columns(microbe_pipeline, c(
    "pipeline", "input_reads", "final_prokaryotic_reads", "final_taxa",
    "final_genera_or_deepest_labels", "bray_distance_spearman_vs_main",
    "differential_direction_concordance_vs_main", "final_read_retention"
  ), "prjna766558_dada2_pipeline_sensitivity.tsv")
  require_columns(microbe_ancombc2, c(
    "taxon", "analysis_status", "error_message"
  ), "prjna766558_dada2_ancombc2_sensitivity.tsv")

  for (column in c(
    "pipeline_direction_stable", "prevalence_threshold_direction_stable",
    "host_sample_level_link_allowed"
  )) {
    set(microbe_genus, j = column, value = as_logical_strict(
      microbe_genus[[column]], paste0("microbe_genus.", column)
    ))
  }
  pair_role_counts <- microbe_sample_qc[, .N, by = .(
    patient_pair_id, tissue_role
  )]
  expected_microbe_alpha_metrics <- c(
    "prokaryotic_reads_min_overlap12", "observed_asv", "shannon", "simpson",
    "rarefied_observed_asv_median", "rarefied_shannon_median",
    "rarefied_simpson_median"
  )
  expected_microbe_pipelines <- c(
    "main_min_overlap12", "min_overlap8", "forward_only",
    "main_plausible_v3v4_length_380_480nt"
  )
  allowed_microbe_tiers <- c(
    "paired_clr_fdr_supported_no_blank",
    "paired_clr_conditional_no_blank",
    "prevalence_shift_conditional_no_blank",
    "directional_exploratory_no_blank",
    "background_no_clear_paired_difference"
  )
  fail_if(
    nrow(microbe_sample_qc) != 42L ||
      uniqueN(microbe_sample_qc$sample_name) != 42L ||
      uniqueN(microbe_sample_qc$run_accession) != 42L ||
      uniqueN(microbe_sample_qc$patient_pair_id) != 21L ||
      nrow(pair_role_counts) != 42L || any(pair_role_counts$N != 1L) ||
      !setequal(pair_role_counts$tissue_role, c("tumor", "paired_non_tumor")) ||
      any(microbe_sample_qc$public_negative_control_run != "not_identified") ||
      any(microbe_sample_qc$blank_control_status !=
            "no_public_negative_control_or_extraction_blank") ||
      nrow(microbe_alpha) != 7L ||
      !setequal(microbe_alpha$metric, expected_microbe_alpha_metrics) ||
      any(microbe_alpha$n_pairs != 21L) ||
      nrow(microbe_beta) != 2L || any(microbe_beta$n_pairs != 21L) ||
      any(!microbe_genus$candidate_tier %chin% allowed_microbe_tiers) ||
      any(microbe_genus$n_pairs != 21L) ||
      any(microbe_genus$host_sample_level_link_allowed) ||
      any(microbe_genus$blank_control_status !=
            "no_public_negative_control_or_extraction_blank") ||
      nrow(microbe_pipeline) != 4L ||
      !setequal(microbe_pipeline$pipeline, expected_microbe_pipelines) ||
      any(!is.finite(microbe_pipeline$final_read_retention)) ||
      any(microbe_pipeline$final_read_retention <= 0),
    "PRJNA766558 完整包的 42/21 配对、无 blank、字段或四流程契约失败。"
  )
  microbe_alpha[, integration_evidence_tier := assign_microbe_alpha_tier(
    interpretation_scope,
    paired_wilcoxon_p,
    paired_wilcoxon_q,
    bootstrap_direction_consistency
  )]
  microbe_beta[, integration_evidence_tier := assign_microbe_beta_tier(
    restricted_permutation_p,
    restricted_permutation_q,
    betadisper_p
  )]
  microbe_genus[, integration_evidence_tier := fcase(
    candidate_tier %chin% c(
      "paired_clr_fdr_supported_no_blank",
      "paired_clr_conditional_no_blank",
      "prevalence_shift_conditional_no_blank"
    ), "T2",
    default = "T1"
  )]
}

metabolite_tier_counts <- metabolite_candidates[, .N, by = final_candidate_tier]
metabolite_tier_text <- collapse_unique(paste0(
  metabolite_tier_counts$final_candidate_tier, "=", metabolite_tier_counts$N
))
metabolite_unverified_identity <- sum(
  metabolite_candidates$annotation_identity_status != "verified_stable_identifier",
  na.rm = TRUE
)
gc_order_confounded <- sum(
  metabolite_inventory$perfect_group_run_order_block,
  na.rm = TRUE
)

microbe_module_tier <- "T0"
microbe_module_result <- "formal DADA2 outputs pending; module intentionally not interpreted"
microbe_module_ceiling <- "仅数据就绪背景；尚无微生物结果结论"
if (microbe_available) {
  microbe_counts <- microbe_genus[, .N, by = candidate_tier]
  microbe_retained_count <- sum(
    microbe_genus$candidate_tier != "background_no_clear_paired_difference"
  )
  microbe_module_signal <-
    any(microbe_genus$integration_evidence_tier == "T2") ||
    any(microbe_alpha$integration_evidence_tier == "T2") ||
    any(microbe_beta$integration_evidence_tier == "T2")
  microbe_module_tier <- if (microbe_module_signal) "T2" else "T1"
  microbe_module_result <- paste0(
    microbe_retained_count,
    " retained genus/deepest-label signals; ",
    collapse_unique(paste0(microbe_counts$candidate_tier, "=", microbe_counts$N)),
    "; no public blank"
  )
  microbe_module_ceiling <- paste(
    "配对分类生态信号；无 blank 限制下不建立因果微生物-宿主轴"
  )
}

tcga_driver_module_tier <- if (
  any(driver_candidates$integrated_tier %chin% c("T4", "T3"))
) {
  "T3"
} else if (any(driver_candidates$integrated_tier == "T2")) {
  "T2"
} else {
  "T1"
}
tcga_heterogeneity_module_tier <- if (
  any(heterogeneity_axes$evidence_tier == "T3")
) {
  "T3"
} else if (any(heterogeneity_axes$evidence_tier == "T2")) {
  "T2"
} else {
  "T1"
}
cao_module_tier <- if (
  any(cao$cross_layer_class %chin% c(
    "three_layer_directional_hypothesis", "credible_reverse_retained"
  ))
) "T2" else "T1"
metabolite_module_tier <- if (
  any(metabolite_candidates$final_candidate_tier %chin% c(
    "robust_rank_and_scale_fdr", "conditional_fdr"
  ))
) "T2" else "T1"
ecms_module_tier <- if (
  any(ecms_factor_context$evidence_tier == "T2") ||
    any(ecms_adjusted_context$evidence_tier == "T2")
) "T2" else "T1"

module_summaries <- data.table(
  module_id = c(
    "TCGA_DRIVER", "TCGA_HETEROGENEITY", "TCGA_ECMS_BENCHMARK",
    "CAO_CROSS_LAYER",
    "PR001876_METABOLOME", "PRJNA766558_MICROBIOME"
  ),
  omics_scope = c(
    "mutation+CNV+RNA", "RNA+miRNA+HM450+mutation+CNV+PROGENy",
    "locked 314-gene ECMS RNA classifier+MOFA+PROGENy",
    "WGBS+RNA+proteome", "targeted GC-MS/LC-MS", "16S tissue microbiome"
  ),
  dataset = c(
    "TCGA-ESCC GDC DR45", "TCGA-ESCC GDC DR45",
    "TCGA-ESCC official ECMS anchor", "Cao 2020",
    "PR001876", "PRJNA766558"
  ),
  design = c(
    "95-patient driver core; mutation mother set 96",
    "94-patient five-layer core",
    paste(
      "official repository-anchor 78 primary; hybrid94 excluded from primary",
      "because extension calibration failed"
    ),
    "9-pair WGBS; 10-pair RNA; mapped proteome",
    "16 early ESCC vs 16 normal per analysis; unpaired primary model",
    if (microbe_available) "21 paired tumor/non-tumor FFPE tissues" else
      "formal DADA2 result not available at this integration run"
  ),
  independence_group = c(
    "TCGA_ESCC_DR45", "TCGA_ESCC_DR45_shared_with_driver",
    "TCGA_ESCC_ECMS_official78_shared_RNA",
    "Cao2020_same_patient_cross_layer",
    "PR001876_early_stage_tissue_subject_overlap_unresolved",
    "PRJNA766558_FFPE_tissue_21_pairs"
  ),
  formal_analysis_available = c(
    TRUE, TRUE, TRUE, TRUE, TRUE, microbe_available
  ),
  host_sample_level_link_to_other_modules_allowed = FALSE,
  exact_driver_event_independent_validation = FALSE,
  missing_layer_interpretation = c(
    "RPPA missingness does not remove driver-core patients",
    "Factor4 level-factor is QC/background, not a negative factor",
    "hybrid94 calibration failure excludes 16 extension patients from primary; not a negative ECMS result",
    "protein not identified/not quantified is not negative",
    "metabolite identity unmapped is not host-axis evidence",
    if (microbe_available)
      "no public blank is a contamination boundary, not automatic rejection" else
      "pending formal DADA2 output; not a negative microbiome result"
  ),
  module_evidence_tier = c(
    tcga_driver_module_tier,
    tcga_heterogeneity_module_tier,
    ecms_module_tier,
    cao_module_tier,
    metabolite_module_tier,
    microbe_module_tier
  ),
  module_result = c(
    paste0(
      nrow(driver_candidates[decision == "strong_patient_level_candidate"]),
      " strong and ",
      nrow(driver_candidates[decision == "conditional_candidate"]),
      " conditional driver candidates; exact independent event validation=0"
    ),
    paste0(
      nrow(heterogeneity_axes[evidence_tier == "T3"]),
      " T3 continuous axes; descriptive k=", selected_cluster$k[[1L]],
      " failed discrete-cluster gate"
    ),
    paste0(
      "official78 primary; ECMS-context T2 factors=",
      sum(ecms_factor_context$evidence_tier == "T2"),
      "; ECMS-adjusted pathway T2 edges=",
      sum(ecms_adjusted_context$evidence_tier == "T2"),
      "; hybrid94 primary=0"
    ),
    paste0(
      sum(cao$cross_layer_class == "three_layer_directional_hypothesis"),
      " three-layer directional and ",
      sum(cao$cross_layer_class == "credible_reverse_retained"),
      " credible reverse; candidate-layer FDR support=",
      sum(cao$fdr_supported_layer_count > 0L), "/", nrow(cao)
    ),
    paste0(
      nrow(metabolite_candidates), " retained features; ",
      metabolite_tier_text, "; identity not fully verified=",
      metabolite_unverified_identity, "; GC perfect run-order blocks=", gc_order_confounded
    ),
    microbe_module_result
  ),
  conclusion_ceiling = c(
    "患者级高优先级候选，不是因果驱动证明",
    "连续异质性假设，不命名稳定离散亚型",
    "同 TCGA 且共享 RNA 表示的 ECMS 解释语境，最高 T2，不是独立验证",
    "小样本同患者方向性校准，不是因果轴",
    "非同患者正交代谢模块校准，不建立宿主患者级关联或精确通路轴",
    microbe_module_ceiling
  ),
  required_next_validation = c(
    "独立 ESCC 精确事件复现、跨平台和功能实验",
    "跨 seed 对齐、独立投影和细胞来源定位",
    "独立 ESCC 队列 ECMS 投影、去共享 RNA 敏感性和细胞来源校准",
    "扩大配对 WGBS-RNA-蛋白队列并复核反向",
    "HMDB/KEGG 身份核验、通路映射及独立队列",
    if (microbe_available)
      "独立配对队列、提取 blank/阴性对照、污染敏感性和功能层验证" else
      "等待正式 DADA2 结果完整发布并通过 artifact manifest"
  )
)

message("[6/7] 生成可追溯 evidence ledger")
driver_ledger <- driver_candidates[, .(
  candidate_id,
  candidate_label = gene_name,
  candidate_type = "driver_candidate",
  evidence_type = "TCGA_driver_composite",
  omics_layer = "genome+transcriptome",
  dataset = "TCGA-ESCC GDC DR45",
  analysis_unit = "unique patient",
  independence_group = "TCGA_ESCC_DR45_driver_core",
  comparison_or_relation = primary_candidate_route,
  effect_measure = fifelse(
    is.finite(spearman_rho_relative_cnv),
    "CNV-RNA Spearman rho", "mutated patient count"
  ),
  effect_value = fifelse(
    is.finite(spearman_rho_relative_cnv),
    spearman_rho_relative_cnv, as.numeric(mutated_patients)
  ),
  p_value = p_relative_cnv,
  q_value = q_relative_cnv,
  direction = fifelse(
    is.finite(spearman_rho_relative_cnv),
    fifelse(spearman_rho_relative_cnv > 0, "positive", "negative"),
    "event_recurrence"
  ),
  support_status = integrated_decision,
  discordance_code = fifelse(
    !is.na(discordance_code) & nzchar(discordance_code),
    discordance_code,
    fifelse(
      as_logical_safe(precise_reverse_dosage),
      "comparable_precise_reverse", "support"
    )
  ),
  layer_availability = fifelse(
    cao_profile_status == "profiled_strong_candidate",
    "TCGA_evaluable_Cao_profiled",
    "TCGA_evaluable_Cao_not_profiled_current_strong_only_not_negative"
  ),
  independent_from_discovery = FALSE,
  countable_as_exact_validation = FALSE,
  evidence_tier = integrated_tier,
  auxiliary_score_100,
  score_is_decision_gate = FALSE,
  evidence_details = paste0(
    "mutation_patients=", mutated_patients,
    ";amplification_frequency=", sprintf("%.4f", amplification_frequency),
    ";evidence_units=", evidence_unit_count,
    ";raw_state_supported_edges=", raw_state_supported_edges,
    ";leakage_gate_supported_edges=", state_supported_edges,
    ";cao_class=", cao_cross_layer_class
  ),
  conclusion_ceiling = integrated_conclusion_ceiling,
  required_next_validation = integrated_required_validation,
  source_file = "results/tcga_escc_driver_candidate_screen.tsv",
  source_row_key = gene_id
)]

axis_ledger <- axis_edges[, .(
  candidate_id = axis_candidate_id,
  candidate_label = sub("^(GENE:|STATE:)", "", axis_candidate_id),
  candidate_type = fifelse(grepl("^STATE:", axis_candidate_id),
                           "heterogeneity_state", "cross_layer_axis"),
  evidence_type = edge_class,
  omics_layer = paste(source_layer, target_layer, sep = "->"),
  dataset = evidence_source,
  analysis_unit = sample_design,
  independence_group,
  comparison_or_relation = paste(source_node, target_node, sep = " -> "),
  effect_measure = association_effect_measure,
  effect_value = association_effect,
  p_value,
  q_value,
  direction = fcase(
    edge_class == "TCGA_ECMS_factor_explanatory",
      "categorical_separation_no_pairwise_direction",
    is.finite(association_effect),
      fifelse(association_effect > 0, "positive", "negative"),
    default = fifelse(
      is.finite(source_effect) & is.finite(target_effect),
      paste0(
        "source_", fifelse(source_effect > 0, "positive", "negative"),
        "_target_", fifelse(target_effect > 0, "positive", "negative")
      ),
      "not_evaluable"
    )
  ),
  support_status,
  discordance_code,
  layer_availability,
  independent_from_discovery = independent_from_tcga_discovery,
  countable_as_exact_validation = countable_as_exact_driver_validation,
  evidence_tier,
  auxiliary_score_100 = unname(tier_score[evidence_tier]),
  score_is_decision_gate = FALSE,
  evidence_details = paste0(
    "source_effect=", ifelse(is.finite(source_effect), source_effect, "NA"),
    ";target_effect=", ifelse(is.finite(target_effect), target_effect, "NA"),
    ";representation_overlap=", representation_overlap,
    ";original_state_edge_id=", ifelse(
      is.na(original_state_edge_id), "NA", original_state_edge_id
    ),
    ";countable_for_T3_T4=", ifelse(
      is.na(countable_for_T3_T4), "NA", countable_for_T3_T4
    ),
    ";leakage_gate_failure_reason=", ifelse(
      is.na(leakage_gate_failure_reason), "NA", leakage_gate_failure_reason
    )
  ),
  conclusion_ceiling,
  required_next_validation,
  source_file,
  source_row_key
)]

metabolite_candidates[, integration_evidence_tier := assign_metabolite_tier(
  final_candidate_tier
)]
metabolite_ledger <- metabolite_candidates[, .(
  candidate_id = paste0("METABOLITE:", analysis_id, ":", feature_id),
  candidate_label = metabolite_name,
  candidate_type = "metabolic_feature",
  evidence_type = "PR001876_targeted_MS_feature",
  omics_layer = "metabolome",
  dataset = paste0("PR001876/", study_id, "/", analysis_id),
  analysis_unit = "unpaired tissue sample; 16 early ESCC vs 16 normal",
  independence_group = gate_independence_group,
  comparison_or_relation = "early_ESCC_vs_normal_tissue",
  effect_measure = "limma transformed-scale effect",
  effect_value = as.numeric(limma_effect_transformed),
  p_value = as.numeric(limma_p),
  q_value = as.numeric(limma_q_evidence_family),
  direction,
  support_status = final_candidate_tier,
  discordance_code = fifelse(
    final_candidate_tier == "run_order_confounded_conditional",
    "technical_failure", "support"
  ),
  layer_availability = "metabolic_feature_evaluable_identity_mapping_incomplete",
  independent_from_discovery = TRUE,
  countable_as_exact_validation = FALSE,
  evidence_tier = integration_evidence_tier,
  auxiliary_score_100 = unname(tier_score[integration_evidence_tier]),
  score_is_decision_gate = FALSE,
  evidence_details = paste0(
    "platform=", platform, ";ion_mode=", ion_mode,
    ";identity=", annotation_identity_status,
    ";conditional_reason=", conditional_reason
  ),
  conclusion_ceiling = evidence_claim_ceiling,
  required_next_validation = "核验化合物身份并在独立代谢组复现；无宿主患者映射时不得建立基因-代谢物患者级轴",
  source_file = "results/pr001876_targeted_ms_candidate_metabolites.tsv",
  source_row_key = feature_id
)]

heterogeneity_ledger <- heterogeneity_axes[, .(
  candidate_id = paste0("STATE:", factor),
  candidate_label = factor,
  candidate_type = "heterogeneity_state",
  evidence_type = "TCGA_continuous_MOFA_axis",
  omics_layer = "multiomics_continuous_state",
  dataset = "TCGA-ESCC GDC DR45",
  analysis_unit =
    "94-patient factor; pathway interpretation from official78 ECMS primary",
  independence_group = "TCGA_ESCC_factor94_ECMS78_shared_RNA",
  comparison_or_relation = paste0(factor, " continuous patient state"),
  effect_measure = "ECMS-adjusted standardized factor beta",
  effect_value = top_pathway_beta,
  p_value = top_pathway_p,
  q_value = top_pathway_q,
  direction = fifelse(top_pathway_beta > 0, "positive", "negative"),
  support_status = fifelse(
    level_factor_soft_flag,
    "level_factor_background", "continuous_axis_retained"
  ),
  discordance_code = fifelse(
    level_factor_soft_flag, "technical_failure", "support"
  ),
  layer_availability = fifelse(
    level_factor_soft_flag,
    "evaluable_but_level_factor_soft_qc", "evaluable"
  ),
  independent_from_discovery = FALSE,
  countable_as_exact_validation = FALSE,
  evidence_tier,
  auxiliary_score_100 = unname(tier_score[evidence_tier]),
  score_is_decision_gate = FALSE,
  evidence_details = paste0(
    "dominant_view=", dominant_view,
    ";top_pathway=", top_pathway,
    ";top_pathway_partial_R2=", top_pathway_partial_r_squared,
    ";driver_edges_supported=", supported_driver_edge_count,
    ";ecms_context_tier=", ecms_context_tier,
    ";ecms_adjusted_pathway_T2_count=", ecms_adjusted_pathway_t2_count,
    ";t3_support_source=", t3_support_source,
    ";discrete_subtype_support=FALSE"
  ),
  conclusion_ceiling,
  required_next_validation,
  source_file =
    "results/tcga_escc_ecms_adjusted_factor_progeny_associations.tsv",
  source_row_key = factor
)]

ecms_extension_boundary_ledger <- data.table(
  candidate_id = "MODULE:ECMS_HYBRID94_EXTENSION",
  candidate_label = "ECMS hybrid94 extension boundary",
  candidate_type = "module_status",
  evidence_type = "ECMS_extension_calibration_boundary",
  omics_layer = "RNA_classifier_benchmark",
  dataset = "TCGA-ESCC GDC DR45",
  analysis_unit = "94 projected patients; only repository-anchor 78 primary",
  independence_group = "TCGA_ESCC_ECMS_shared_RNA",
  comparison_or_relation = "official78_anchor_vs_GDC_reprojection_overlap_gate",
  effect_measure = "overall extension calibration pass",
  effect_value = 0,
  p_value = NA_real_,
  q_value = NA_real_,
  direction = "boundary_not_negative",
  support_status = "hybrid94_excluded_from_primary_calibration_failed",
  discordance_code = "not_comparable_as_primary_extension",
  layer_availability = "official78_primary_available_hybrid16_conditional_only",
  independent_from_discovery = FALSE,
  countable_as_exact_validation = FALSE,
  evidence_tier = "T0",
  auxiliary_score_100 = 0,
  score_is_decision_gate = FALSE,
  evidence_details = paste0(
    "extension_calibration_pass=FALSE;primary_n=78;hybrid_primary_n=0;",
    "shared_RNA_representation=TRUE"
  ),
  conclusion_ceiling =
    "额外16例只作条件敏感性；失败不否定官方78例 anchor",
  required_next_validation =
    "在独立 TPM 处理队列复核投影一致性；未通过前不扩展 primary",
  source_file = "results/tcga_escc_ecms_projection_calibration.tsv",
  source_row_key = "overall_extension_calibration_pass"
)

microbe_ledger <- if (microbe_available) {
  genus_candidates <- microbe_genus[
    candidate_tier != "background_no_clear_paired_difference"
  ]
  rbindlist(list(
    genus_candidates[, .(
      candidate_id = paste0("MICROBE:", genus_label),
      candidate_label = genus_label,
      candidate_type = "microbial_taxon",
      evidence_type = "PRJNA766558_paired_CLR",
      omics_layer = "microbiome_16S_taxonomy",
      dataset = "PRJNA766558",
      analysis_unit = "21 paired tumor/non-tumor tissues",
      independence_group,
      comparison_or_relation = "tumor_vs_paired_non_tumor",
      effect_measure = "median paired CLR difference",
      effect_value = median_clr_difference_tumor_minus_normal,
      p_value = paired_wilcoxon_p,
      q_value = paired_wilcoxon_q,
      direction = effect_direction,
      support_status = candidate_tier,
      discordance_code = "support",
      layer_availability = "taxonomic_ecology_evaluable_no_public_blank",
      independent_from_discovery = TRUE,
      countable_as_exact_validation = FALSE,
      evidence_tier = integration_evidence_tier,
      auxiliary_score_100 = unname(tier_score[integration_evidence_tier]),
      score_is_decision_gate = FALSE,
      evidence_details = paste0(
        "pipeline_direction_stable=", pipeline_direction_stable,
        ";prevalence_direction_stable=", prevalence_threshold_direction_stable,
        ";bootstrap_direction_consistency=", bootstrap_direction_consistency,
        ";blank=", blank_control_status,
        ";contamination=", contamination_risk
      ),
      conclusion_ceiling,
      required_next_validation = "独立配对队列、提取 blank/阴性对照、污染敏感性及功能层验证",
      source_file = "results/prjna766558_dada2_genus_paired_differential.tsv",
      source_row_key = genus_label
    )],
    microbe_alpha[, .(
      candidate_id = paste0("MICROBIOME_ALPHA:", metric),
      candidate_label = metric,
      candidate_type = "microbiome_module_metric",
      evidence_type = "PRJNA766558_paired_alpha_diversity",
      omics_layer = "microbiome_16S_alpha_diversity",
      dataset = "PRJNA766558",
      analysis_unit = paste0(n_pairs, " paired tissues"),
      independence_group = "PRJNA766558_FFPE_tissue_21_pairs",
      comparison_or_relation = "tumor_vs_paired_non_tumor",
      effect_measure = "median paired difference",
      effect_value = median_paired_difference_tumor_minus_normal,
      p_value = paired_wilcoxon_p,
      q_value = paired_wilcoxon_q,
      direction = fifelse(
        median_paired_difference_tumor_minus_normal > 0,
        "higher_in_tumor", "lower_in_tumor"
      ),
      support_status = fifelse(
        integration_evidence_tier == "T0",
        "depth_diagnostic_not_diversity",
        fifelse(
          integration_evidence_tier == "T2",
          "rarefied_alpha_conditional_supported_no_blank",
          "depth_sensitive_or_nonsignificant_description"
        )
      ),
      discordance_code = fifelse(
        integration_evidence_tier == "T2",
        "support",
        fifelse(integration_evidence_tier == "T0", "not_comparable",
                "same_direction_nonsignificant")
      ),
      layer_availability = interpretation_scope,
      independent_from_discovery = TRUE,
      countable_as_exact_validation = FALSE,
      evidence_tier = integration_evidence_tier,
      auxiliary_score_100 = unname(tier_score[integration_evidence_tier]),
      score_is_decision_gate = FALSE,
      evidence_details = paste0(
        "bootstrap_direction_consistency=", bootstrap_direction_consistency,
        ";interpretation_scope=", interpretation_scope,
        ";no_public_blank"
      ),
      conclusion_ceiling = fifelse(
        integration_evidence_tier == "T0",
        "测序深度 QC，不是 alpha 多样性证据",
        "配对 16S 群落描述；无 blank 时不建立因果微生物-宿主轴"
      ),
      required_next_validation = "独立配对队列和 blank/污染敏感性",
      source_file = "results/prjna766558_dada2_alpha_paired_tests.tsv",
      source_row_key = metric
    )],
    microbe_beta[, .(
      candidate_id = paste0("MICROBIOME_BETA:", gsub("[^A-Za-z0-9]+", "_", distance)),
      candidate_label = distance,
      candidate_type = "microbiome_module_metric",
      evidence_type = "PRJNA766558_restricted_PERMANOVA",
      omics_layer = "microbiome_16S_beta_diversity",
      dataset = "PRJNA766558",
      analysis_unit = paste0(n_pairs, " paired tissues"),
      independence_group = "PRJNA766558_FFPE_tissue_21_pairs",
      comparison_or_relation = "tumor_vs_paired_non_tumor",
      effect_measure = "restricted PERMANOVA R2",
      effect_value = permanova_r2,
      p_value = restricted_permutation_p,
      q_value = restricted_permutation_q,
      direction = "community_composition_difference",
      support_status = fifelse(
        integration_evidence_tier == "T2",
        "paired_permanova_supported_no_blank", "descriptive_or_dispersion_sensitive"
      ),
      discordance_code = fifelse(
        integration_evidence_tier == "T2",
        "support", "same_direction_nonsignificant"
      ),
      layer_availability = "evaluable_no_public_blank",
      independent_from_discovery = TRUE,
      countable_as_exact_validation = FALSE,
      evidence_tier = integration_evidence_tier,
      auxiliary_score_100 = unname(tier_score[integration_evidence_tier]),
      score_is_decision_gate = FALSE,
      evidence_details = paste0("betadisper_p=", betadisper_p, ";no_public_blank"),
      conclusion_ceiling = "配对 16S 群落差异；无 blank 时不建立因果微生物-宿主轴",
      required_next_validation = "独立配对队列和 blank/污染敏感性",
      source_file = "results/prjna766558_dada2_beta_tests.tsv",
      source_row_key = distance
    )]
  ), use.names = TRUE, fill = TRUE)
} else {
  data.table(
    candidate_id = "MODULE:PRJNA766558",
    candidate_label = "PRJNA766558 microbiome module",
    candidate_type = "module_status",
    evidence_type = "formal_result_pending",
    omics_layer = "microbiome_16S",
    dataset = "PRJNA766558",
    analysis_unit = "21 paired tissues registered; formal DADA2 output pending",
    independence_group = "PRJNA766558_FFPE_tissue_21_pairs",
    comparison_or_relation = "not_yet_integrated",
    effect_measure = NA_character_,
    effect_value = NA_real_,
    p_value = NA_real_,
    q_value = NA_real_,
    direction = "not_evaluable",
    support_status = "pending_not_negative",
    discordance_code = "not_comparable",
    layer_availability = "formal_DADA2_result_pending_not_negative",
    independent_from_discovery = TRUE,
    countable_as_exact_validation = FALSE,
    evidence_tier = "T0",
    auxiliary_score_100 = 0,
    score_is_decision_gate = FALSE,
    evidence_details = "No partial output was interpreted",
    conclusion_ceiling = "仅数据就绪背景；尚无微生物结果结论",
    required_next_validation = "等待正式 DADA2 结果完整发布并通过 artifact manifest",
    source_file = "results/prjna766558_dada2_artifact_manifest.tsv",
    source_row_key = "module_pending"
  )
}

evidence_ledger <- rbindlist(list(
  driver_ledger, axis_ledger, metabolite_ledger,
  heterogeneity_ledger, ecms_extension_boundary_ledger, microbe_ledger
), use.names = TRUE, fill = TRUE)
evidence_ledger[, cross_dataset_patient_link_allowed := FALSE]
fail_if(any(!evidence_ledger$evidence_tier %chin% names(tier_rank)),
        "evidence ledger 出现未知 T 层级。")
evidence_ledger[, evidence_id := sprintf("EVIDENCE_%06d", seq_len(.N))]
setcolorder(evidence_ledger, c("evidence_id", setdiff(names(evidence_ledger), "evidence_id")))
evidence_ledger[, tier_order_internal := unname(tier_rank[evidence_tier])]
setorder(
  evidence_ledger,
  -tier_order_internal,
  candidate_type,
  candidate_id,
  q_value,
  p_value,
  na.last = TRUE
)
evidence_ledger[, tier_order_internal := NULL]

message("[7/7] 写出正式整合包、摘要、QA 和 artifact manifest")
formal_tables <- list(
  escc_multiomics_integration_tier_definitions.tsv = tier_definitions,
  escc_multiomics_evidence_ledger.tsv = evidence_ledger,
  escc_multiomics_integrated_driver_candidates.tsv = driver_candidates,
  escc_multiomics_integrated_axis_edges.tsv = axis_edges,
  escc_multiomics_integrated_axis_summary.tsv = axis_summary,
  escc_multiomics_heterogeneity_axes.tsv = heterogeneity_axes,
  escc_multiomics_module_summaries.tsv = module_summaries
)
invisible(Map(stage_tsv, formal_tables, names(formal_tables)))

priority_axes <- axis_summary[axis_tier %chin% c("T4", "T3")]
priority_axis_lines <- if (nrow(priority_axes)) {
  paste0(
    "- `", priority_axes$gene_name, "`：", priority_axes$axis_class,
    "；", priority_axes$axis_tier, "；", priority_axes$conclusion_ceiling, "。"
  )
} else {
  "- 当前没有 T3/T4 候选轴。"
}

summary_lines <- c(
  "# ESCC 多组学证据整合摘要",
  "",
  "## 证据裁决",
  "",
  paste0(
    "- evidence ledger 共 ", nrow(evidence_ledger), " 条；T0-T4 计数：",
    collapse_unique(paste0(
      evidence_ledger[, .N, by = evidence_tier]$evidence_tier, "=",
      evidence_ledger[, .N, by = evidence_tier]$N
    )), "。"
  ),
  paste0(
    "- driver 候选：", nrow(driver_candidates), " 个（strong=",
    sum(driver_candidates$decision == "strong_patient_level_candidate"),
    "；conditional=", sum(driver_candidates$decision == "conditional_candidate"),
    "）；exploratory driver 未进入正式整合候选表。"
  ),
  paste0(
    "- T4 高优先级验证 driver：",
    if (any(driver_candidates$integrated_tier == "T4"))
      paste(driver_candidates[integrated_tier == "T4", gene_name], collapse = "、")
    else "无",
    "。T4 表示下一步优先级，不表示已验证因果。"
  ),
  paste0(
    "- driver–MOFA 表示泄漏门禁：41 条原边中 ",
    sum(as_logical_safe(leakage_candidate_summary$countable_for_T3_T4)),
    " 条可计入 T3/T4；其余仅条件保留或背景，",
    "PROGENy 边不替代 drop-both 门禁。"
  ),
  paste0(
    "- ECMS 基准：extension calibration=FALSE，仅官方 78 例为 primary；",
    "hybrid94 只作边界敏感性。ECMS-context T2 因子=",
    sum(ecms_factor_context$evidence_tier == "T2"),
    "；ECMS-adjusted Factor→PROGENy T2 边=",
    sum(ecms_adjusted_context$evidence_tier == "T2"), "。"
  ),
  "- ECMS、PROGENy 和 MOFA 共享 RNA 表示，相关边单独最高 T2；综合 T3 必须由 leakage-controlled driver 边提供。",
  "- 数值分数仅用于同层排序；T 层级由来源结构、独立性、冲突类型和结论上限裁决。",
  "",
  "## 关键候选轴",
  "",
  priority_axis_lines,
  "",
  "## 条件保留与反向",
  "",
  paste0(
    "- Cao 三层方向性候选：",
    paste(axis_summary[axis_class == "three_layer_directional_axis_conditional", gene_name],
          collapse = "、"),
    "；候选层 FDR 未通过，按 T2 条件保留。"
  ),
  paste0(
    "- 情境特异反向保留：",
    paste(axis_summary[axis_class == "context_specific_reverse_axis_retained", gene_name],
          collapse = "、"),
    "；反向不被其他分数抵消，需独立配对队列和蛋白/细胞来源复核。"
  ),
  if (microbe_available) {
    "- 蛋白未鉴定或未定量、代谢物身份未映射均不记为阴性；微生物无公开 blank 作为结论边界保留。"
  } else {
    "- 蛋白未鉴定或未定量、代谢物身份未映射、微生物正式结果待发布均不记为阴性。"
  },
  "",
  "## 疾病异质性",
  "",
  paste0(
    "- 当前保留 ", sum(!heterogeneity_axes$level_factor_soft_flag),
    " 个非 level-factor 连续轴；T3 轴 ",
    paste(heterogeneity_axes[evidence_tier == "T3", factor], collapse = "、"), "。"
  ),
  "- 连续异质性 T3 同时要求 leakage T3、官方78 ECMS-context T2 和 ECMS-adjusted pathway T2；Factor4 始终只作 level-factor QC。",
  paste0(
    "- 描述性 k=", selected_cluster$k[[1L]],
    "：SNF silhouette=", sprintf("%.3f", selected_cluster$snf_mean_silhouette[[1L]]),
    "，MOFA PAC=", sprintf("%.3f", selected_cluster$mofa_consensus_pac[[1L]]),
    "，ARI=", sprintf("%.3f", selected_cluster$snf_mofa_adjusted_rand[[1L]]),
    "；离散亚型门禁未通过，不命名新亚型。"
  ),
  "",
  "## ECMS、代谢与微生物模块",
  "",
  paste0("- ECMS：", module_summaries[
    module_id == "TCGA_ECMS_BENCHMARK", module_result
  ], "。"),
  paste0("- PR001876：", module_summaries[module_id == "PR001876_METABOLOME", module_result], "。"),
  paste0(
    "- PRJNA766558：",
    module_summaries[module_id == "PRJNA766558_MICROBIOME", module_result], "。"
  ),
  "- PR001876 和 PRJNA766558 均不与 TCGA/Cao 伪造患者级映射；所有代谢/微生物证据最高 T2，只作模块级正交校准。",
  "",
  "## 总结边界",
  "",
  if (microbe_available) {
    "本整合包揭示的是可证伪的公共数据候选网络、跨层轴和连续异质性框架。TCGA 同队列桥接、Cao 小样本方向性、PR001876 非配对代谢特征和 PRJNA766558 无 blank 的 16S 信号均不能单独证明因果机制或治疗靶点。"
  } else {
    "本整合包揭示的是可证伪的公共数据候选网络、跨层轴和连续异质性框架。TCGA 同队列桥接、Cao 小样本方向性和 PR001876 非配对代谢特征均不能单独证明因果机制或治疗靶点；PRJNA766558 正式 DADA2 结果尚未纳入，不能提前声称微生物信号。"
  }
)
summary_filename <- "escc_multiomics_integration_summary.md"
writeLines(summary_lines, file.path(stage_dir, summary_filename), useBytes = TRUE)

formal_filenames <- c(names(formal_tables), summary_filename)
fail_if(!identical(formal_filenames, integration_formal_filenames),
        "整合正式 artifact 列表与 --fields-only 契约不一致。")
artifact_paths <- file.path(stage_dir, formal_filenames)
fail_if(any(!file_exists(artifact_paths)), "正式整合 artifact 未完整生成。")
upstream_manifest_paths <- c(
  file.path(results_dir, "tcga_escc_driver_core_artifact_manifest.tsv"),
  file.path(results_dir, "tcga_escc_driver_state_artifact_manifest.tsv"),
  leakage_manifest_path,
  file.path(results_dir, "cao2020_cross_layer_artifact_manifest.tsv"),
  file.path(results_dir, "tcga_escc_heterogeneity_artifact_manifest.tsv"),
  file.path(results_dir, "pr001876_targeted_ms_artifact_manifest.tsv"),
  ecms_manifest_path
)
if (microbe_available) {
  upstream_manifest_paths <- c(upstream_manifest_paths, prjna_manifest_path)
}
fail_if(any(!file_exists(upstream_manifest_paths)),
        "整合 manifest 准备时上游 manifest 缺失。")
upstream_manifest_sha256 <- paste0(
  basename(upstream_manifest_paths), "=",
  vapply(upstream_manifest_paths, sha256_file, character(1)),
  collapse = ";"
)
artifact_manifest <- data.table(
  artifact = formal_filenames,
  relative_path = file.path("results", formal_filenames),
  file_size_bytes = as.numeric(file_info(artifact_paths)$size),
  sha256 = vapply(
    artifact_paths,
    digest,
    FUN.VALUE = character(1),
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  ),
  generated_date = as.character(Sys.Date()),
  generation_script = "scripts/20_integrate_escc_multiomics_evidence.R",
  executed_script_sha256 = executed_script_sha256,
  upstream_manifest_sha256 = upstream_manifest_sha256,
  ecms_model_commit = locked_ecms_commit,
  ecms_extension_calibration_pass = ecms_extension_calibration_pass,
  ecms_primary_patient_count = ecms_primary_patient_count,
  prjna766558_complete_family = microbe_available,
  source_boundary = if (microbe_available)
    paste(
      "TCGA+leakage_audit+ECMS official78 primary+Cao2020+PR001876+",
      "PRJNA766558 exact formal families; shared-RNA and unpaired-module caps applied"
    ) else paste(
      "TCGA+leakage_audit+ECMS official78 primary+Cao2020+PR001876 exact",
      "formal families; PRJNA766558 pending not negative"
    ),
  status = "verified"
)
manifest_filename <- integration_manifest_filename
stage_tsv(artifact_manifest, manifest_filename)

# 轻量回读与关键语义断言。
for (filename in names(formal_tables)) {
  reread <- read_tsv(file.path(stage_dir, filename))
  fail_if(!nrow(reread), paste("正式整合 TSV 为空：", filename))
}
fail_if(any(driver_candidates$decision == "exploratory_signal"),
        "探索性 driver 被错误写入正式候选表。")
fail_if(any(driver_candidates$exact_independent_validation_present),
        "当前结果不应出现精确 driver 事件独立验证。")
fail_if(any(axis_edges$countable_as_exact_driver_validation),
        "当前轴边不应计作精确 driver 独立验证。")
fail_if(
  any(axis_edges$edge_class == "TCGA_driver_state_internal" &
        axis_edges$evidence_tier == "T3" &
        !as_logical_safe(axis_edges$countable_for_T3_T4)),
  "TCGA driver-state T3 边越过了预锁定泄漏门禁。"
)
fail_if(
  any(axis_edges$edge_class == "TCGA_driver_state_internal" &
        axis_edges$target_layer == "RNA_derived_pathway_activity" &
        axis_edges$evidence_tier %chin% c("T3", "T4")),
  "PROGENy 同队列解释边被错误用于 T3/T4 泄漏门禁。"
)
fail_if(
  any(axis_edges$edge_class %chin% c(
        "TCGA_ECMS_factor_explanatory",
        "TCGA_ECMS_adjusted_factor_PROGENy"
      ) &
        axis_edges$evidence_tier %chin% c("T3", "T4")),
  "ECMS/PROGENy 共享 RNA 解释边被错误升级到 T3/T4。"
)
fail_if(
  ecms_factor_context[factor == "Factor4", any(evidence_tier != "T0")] ||
    ecms_adjusted_context[factor == "Factor4", any(evidence_tier != "T0")] ||
    heterogeneity_axes[factor == "Factor4", any(evidence_tier != "T0")],
  "Factor4 HM450 level-factor 未在 ECMS-factor、adjusted pathway 和 heterogeneity 三层同时锁定 T0。"
)
fail_if(
  any(axis_summary$axis_tier %chin% c("T3", "T4") &
        axis_summary$state_supported_edges < 1L),
  "driver 轴摘要 T3/T4 绕过了 leakage-controlled 状态边门禁。"
)
fail_if(
  any(heterogeneity_axes$evidence_tier == "T3" &
        (
          heterogeneity_axes$supported_driver_edge_count < 1L |
            !heterogeneity_axes$ecms_context_t2 |
            heterogeneity_axes$ecms_adjusted_pathway_t2_count < 1L
        )),
  "连续异质性 T3 轴缺少 leakage T3、ECMS-context T2 或 adjusted pathway T2。"
)
fail_if(
  !any(as_logical_safe(leakage_candidate_summary$countable_for_T3_T4)) &&
    (any(axis_summary$axis_tier %chin% c("T3", "T4")) ||
       any(heterogeneity_axes$evidence_tier == "T3")),
  "零条泄漏门禁通过时仍出现 driver-state/异质性 T3/T4。"
)
fail_if(any(axis_summary$cross_cohort_patient_link_created),
        "检测到禁止的跨队列患者级拼接标志。")
fail_if(any(module_summaries$host_sample_level_link_to_other_modules_allowed),
        "模块摘要错误允许了跨队列患者级连接。")
fail_if(any(evidence_ledger$cross_dataset_patient_link_allowed),
        "evidence ledger 错误允许跨队列患者级连接。")
fail_if(any(is_forbidden_unpaired_host_source(axis_edges$evidence_source)),
        "代谢或微生物来源被错误写成宿主轴边。")
fail_if(any(metabolite_ledger$evidence_tier %chin% c("T3", "T4")),
        "PR001876 非同患者证据越过 T2 上限。")
fail_if(any(microbe_ledger$evidence_tier %chin% c("T3", "T4")),
        "PRJNA766558 非同患者宿主证据越过 T2 上限。")
fail_if(any(module_summaries[
  module_id %chin% c(
    "TCGA_ECMS_BENCHMARK", "PR001876_METABOLOME",
    "PRJNA766558_MICROBIOME"
  ), module_evidence_tier
] %chin% c("T3", "T4")),
"共享 RNA 或非同患者模块越过 T2 上限。")
fail_if(any(axis_summary$axis_tier == "T4" &
              axis_summary$ecms_contextualized_leakage_edge_count < 1L),
        "T4 候选轴未在官方78 ECMS 语境中完成 leakage-controlled 桥接。")
fail_if(!setequal(
  driver_candidates[integrated_tier == "T4", gene_name],
  axis_summary[axis_tier == "T4", gene_name]
), "driver candidate T4 与最终集成轴 T4 未严格同步。")
fail_if(any(
  driver_candidates[
    decision == "strong_patient_level_candidate", integrated_tier
  ] !=
    axis_summary$axis_tier[match(
      driver_candidates[
        decision == "strong_patient_level_candidate", gene_name
      ],
      axis_summary$gene_name
    )]
), "driver candidate 与最终集成轴的 T0-T4 层级未逐基因同步。")
fail_if(any(axis_summary$axis_class == "context_specific_reverse_axis_retained" &
              axis_summary$axis_tier != "T2"),
        "情境特异反向未按 T2 条件保留。")
fail_if(microbe_available && !any(evidence_ledger$dataset == "PRJNA766558"),
        "DADA2 正式结果可用但未进入 evidence ledger。")
fail_if(microbe_available && microbe_alpha[
  interpretation_scope == "paired_final_depth_diagnostic_not_diversity",
  any(integration_evidence_tier != "T0")
], "测序深度诊断被错误写成 alpha 多样性候选。")
fail_if(!microbe_available &&
          !any(evidence_ledger$support_status == "pending_not_negative"),
        "DADA2 结果缺失未按 pending-not-negative 留痕。")
self_test_full_run <- run_self_tests()
fail_if(nrow(self_test_full_run) != 7L || any(!self_test_full_run$pass),
        "正式发布前 7 个纯内存 fixture 未全部通过。")

# 先原子发布正式 artifact，最后发布 manifest，避免中断时出现
# “新 manifest 指向旧/不完整文件”的窗口。
for (filename in formal_filenames) {
  atomic_publish_file(
    file.path(stage_dir, filename),
    file.path(results_dir, filename)
  )
}
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
published_manifest_path <- file.path(results_dir, manifest_filename)
fail_if(
  digest(
    published_manifest_path,
    algo = "sha256", file = TRUE, serialize = FALSE
  ) != digest(
    file.path(stage_dir, manifest_filename),
    algo = "sha256", file = TRUE, serialize = FALSE
  ),
  "发布后的整合 artifact manifest 与暂存版本不一致。"
)
published_manifest <- read_tsv(published_manifest_path)
fail_if(
  nrow(published_manifest) != nrow(artifact_manifest) ||
    any(published_manifest$sha256 != artifact_manifest$sha256) ||
    any(published_manifest$file_size_bytes != artifact_manifest$file_size_bytes),
  "发布后的整合 artifact manifest 回读不一致。"
)

qa_lines <- c(
  "# ESCC 多组学整合层 QA（2026-07-12）",
  "",
  "本文件是运行时检查记录，不是项目当前状态源。",
  "",
  "- script21 7+manifest 与 ECMS 6+manifest exact family/generator SHA：通过。",
  paste0(
    "- PRJNA766558 12+manifest exact family：",
    if (microbe_available) "通过" else "尚未发布，pending-not-negative",
    "。"
  ),
  paste0("- 上游 manifest SHA 冻结：", upstream_manifest_sha256, "。"),
  paste0("- strong+conditional driver：", nrow(driver_candidates), "；exploratory 纳入数：0。"),
  paste0("- evidence ledger：", nrow(evidence_ledger), " 条；axis edges：", nrow(axis_edges), " 条。"),
  paste0("- strong driver axis summary：", nrow(axis_summary), " 条。"),
  paste0("- 连续异质性轴：", nrow(heterogeneity_axes), "；离散亚型支持：FALSE。"),
  paste0("- PR001876 候选特征：", nrow(metabolite_candidates), "；仅模块级校准。"),
  paste0("- PRJNA766558 正式 DADA2 纳入：", microbe_available, "；若 FALSE 已记 pending-not-negative。"),
  paste0(
    "- ECMS：official78 primary；hybrid94 primary=0；context T2 factors=",
    sum(ecms_factor_context$evidence_tier == "T2"),
    "；adjusted pathway T2 edges=",
    sum(ecms_adjusted_context$evidence_tier == "T2"), "。"
  ),
  "- Factor4 HM450 level-factor：ECMS-factor、adjusted pathway、heterogeneity 全部 T0。",
  paste0(
    "- driver–MOFA 泄漏门禁：41 条原边中 ",
    sum(as_logical_safe(leakage_candidate_summary$countable_for_T3_T4)),
    " 条可计入 T3/T4；PROGENy 不替代 drop-both 门禁。"
  ),
  "- 跨队列患者级拼接：0；精确 driver 事件独立验证：0。",
  "- ECMS/PROGENy、PR001876、PRJNA766558 证据上限：T2；代谢/微生物宿主轴：0。",
  "- 7 个纯内存 fixture：7/7 通过，含 zero-gate 与 invalid-logical fail-closed。",
  "- 缺失蛋白、未映射代谢物身份及待发布微生物结果均未写成阴性。",
  "- 情境特异反向按 T2 保留，未被数值分数抵消。",
  paste0("- 正式 artifact：", nrow(artifact_manifest), " 个；发布后大小与 SHA256 回读一致。"),
  paste0("- R ", getRversion(), "；data.table ", packageVersion("data.table"),
         "；digest ", packageVersion("digest"), "；fs ", packageVersion("fs"), "。"),
  "- ResearchDataHub 原始资料、项目原始数据和投稿包均未修改或删除。"
)
qa_filename <- "escc_multiomics_integration_qa_20260712.md"
qa_stage_path <- file.path(stage_dir, qa_filename)
writeLines(qa_lines, qa_stage_path, useBytes = TRUE)
atomic_publish_file(qa_stage_path, file.path(work_checks_dir, qa_filename))

dir_delete(stage_dir)
message(
  "ESCC 多组学整合完成：", nrow(evidence_ledger), " 条 evidence；",
  nrow(driver_candidates), " 个 driver 候选；", nrow(axis_edges), " 条轴边；",
  "PRJNA766558 正式纳入=", microbe_available, "。"
)
