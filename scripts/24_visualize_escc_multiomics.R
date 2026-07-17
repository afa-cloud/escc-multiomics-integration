#!/usr/bin/env Rscript

# ESCC 多组学投稿级主图（R-only）。
#
# 设计边界：
# 1) 所有统计编码只读取 results/ 的正式源表；不从旧图片、摘要文字或
#    手工标注反推数值。
# 2) script 21、script 18、script 20 与 ECMS 投影必须是完整 exact artifact
#    family，且 manifest 中的生成脚本 SHA256 必须与当前脚本一致。
# 3) ECMS primary 只允许官方 78 例；Factor4 始终为 T0/HM450 level-factor。
# 4) Cao promoter、RNA、protein 使用独立数值尺度；缺失不画成 0。
# 5) PR001876 与 PRJNA766558 不建立患者级连接；整合网络只画无箭头关联线。
# 6) 图件先在 _work/intermediate/ 原子暂存，20 个图文件验证并发布后，
#    figure artifact manifest 才最后原子发布到 results/。

args <- commandArgs(trailingOnly = TRUE)
fields_only <- "--fields-only" %in% args
validate_only <- "--validate-only" %in% args
self_test <- "--self-test" %in% args
finalize_visual_qa_args <- grep(
  "^--finalize-visual-qa=", args, value = TRUE
)
if (length(finalize_visual_qa_args) > 1L) {
  stop("--finalize-visual-qa 只能指定一次。", call. = FALSE)
}
finalize_visual_qa <- length(finalize_visual_qa_args) == 1L
finalize_visual_qa_input <- if (finalize_visual_qa) {
  sub("^--finalize-visual-qa=", "", finalize_visual_qa_args[[1L]])
} else {
  NULL
}
unknown_args <- setdiff(
  args,
  c("--fields-only", "--validate-only", "--self-test", finalize_visual_qa_args)
)
if (length(unknown_args)) {
  stop("未知参数：", paste(unknown_args, collapse = ", "), call. = FALSE)
}
if (sum(c(fields_only, validate_only, self_test, finalize_visual_qa)) > 1L) {
  stop(paste(
    "--fields-only、--validate-only、--self-test 与",
    "--finalize-visual-qa 不能同时使用。"
  ),
       call. = FALSE)
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("无法唯一定位当前脚本。", call. = FALSE)
script_path <- normalizePath(
  sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE
)
project_root <- dirname(dirname(script_path))
results_dir <- file.path(project_root, "results")
figures_dir <- file.path(project_root, "figures")
work_intermediate_dir <- file.path(project_root, "_work", "intermediate")
work_checks_dir <- file.path(project_root, "_work", "checks")

figure_specs <- data.frame(
  figure_id = paste0("Figure", 1:5),
  file_stem = c(
    "escc_multiomics_figure1_driver_state_representation_overlap",
    "escc_multiomics_figure2_ecms_continuous_states",
    "escc_multiomics_figure3_cao_cross_layer_calibration",
    "escc_multiomics_figure4_orthogonal_metabolome_microbiome",
    "escc_multiomics_figure5_integrated_candidate_network"
  ),
  width_mm = rep(183, 5),
  height_mm = c(150, 155, 150, 175, 150),
  claim = c(
    paste(
      "Nine event-state associations remain after representation-overlap",
      "auditing, with the PIK3CA-Factor1 CNV family as the leading",
      "within-TCGA validation priority."
    ),
    paste(
      "Factor1 and Factor3 provide the most concentrated continuous-state",
      "signal across the official ECMS anchor and ECMS-adjusted pathways."
    ),
    paste(
      "Same-patient Cao data support a three-layer directional GNAS",
      "hypothesis and a falsifiable ZNF750 reverse pattern."
    ),
    paste(
      "PR001876 contains 83 rank-and-scale robust LC-MS features, while",
      "paired tissue 16S supports diversity and genus-level ecological shifts."
    ),
    paste(
      "A PIK3CA-Factor1-centered network prioritizes a compact set of",
      "within-source event-state hypotheses for independent validation."
    )
  ),
  archetype = c(
    "asymmetric mixed-modality figure",
    "quantitative grid",
    "asymmetric mixed-modality figure",
    "asymmetric mixed-modality figure",
    "asymmetric mixed-modality figure"
  ),
  panel_source_contract = c(
    paste(
      "a=the nine retained script21 candidate-summary edges;",
      "b=drop-both direction/magnitude/support retention;",
      "c=PIK3CA CNV-RNA and Factor1 evidence"
    ),
    paste(
      "a=official78 patient-level Factor1/Factor3 distributions;",
      "b=ECMS eta-squared and q; c=ECMS-adjusted pathway effects;",
      "d=compact continuous-axis evidence summary"
    ),
    paste(
      "a=GNAS patient-level WGBS/RNA/protein effects;",
      "b=ZNF750 patient-level promoter/RNA effects;",
      "c=two-gene layer-specific summary effects and q values"
    ),
    paste(
      "a=metabolite tier counts; b=top robust LC-MS effects;",
      "c=paired alpha-diversity effects; d=top stable paired CLR taxa"
    ),
    paste(
      "a=nine retained event-factor edges plus PIK3CA dosage and selected",
      "Factor1/Factor3 pathway context; b=T4/T3 validation-priority ladder"
    )
  ),
  filter_contract = c(
    paste(
      "only countable_for_T3_T4=TRUE edges in the main visual; all 41 edges",
      "and execution/alignment diagnostics remain in formal source tables"
    ),
    paste(
      "Factor1/Factor3 only; official78 primary patient distributions;",
      "formal ECMS-adjusted T2 pathway edges; extension16 and k diagnostics omitted"
    ),
    paste(
      "GNAS and ZNF750 only in the main visual; missing protein remains absent",
      "rather than zero; all 12 candidates remain in formal source tables"
    ),
    paste(
      "all tier counts; top ten robust LC-MS features by formal q/rank;",
      "three paired alpha metrics and top eight stable paired CLR taxa"
    ),
    paste(
      "nine retained event-factor edges; PIK3CA dosage; Factor1/Factor3",
      "ECMS context and selected adjusted pathways; T4/T3 genes only"
    )
  ),
  statistical_unit = c(
    "94 TCGA patients; retained original edge and three-seed drop-both retention",
    "78 official-anchor TCGA patients; factor/pathway association",
    "Cao patient pair and prespecified strong-driver candidate",
    "PR001876 analysis-feature and PRJNA766558 patient pair",
    "formal retained edge, pathway and candidate-priority rows; no new inferential test"
  ),
  sample_structure = c(
    "3 prespecified seeds; same TCGA patients; no independent validation",
    "official78 primary; extension16 excluded from the main positive visual",
    "9 paired WGBS; 10 paired RNA; mapped patient-level protein ratios",
    paste(
      "three PR001876 main analyses (AN004960/AN004962/AN004963), each 16 vs 16;",
      "96 analysis-sample keys; 21/41 published patient pairs available as an",
      "FFPE 16S reproducibility-only subset; no cross-dataset pairing"
    ),
    "TCGA event/factor/pathway rows remain source-internal"
  ),
  effect_definition = c(
    "original event-factor effect; direction, >=50% magnitude and support retention rates",
    "ECMS eta-squared; standardized adjusted beta; partial R-squared",
    "promoter median paired delta-beta; RNA limma logFC; protein median log2 T/N",
    "analysis-scale limma effect; paired alpha median difference; PERMANOVA R2; genus CLR difference",
    "categorical T3/T4 validation priority and source-table association effects"
  ),
  p_q_family = c(
    "per-model per-event-type BH family from script21; no manual stars",
    "formal ECMS factor and adjusted pathway BH families; cluster metrics descriptive",
    "candidate-layer BH fields from script15; no candidate passes are invented",
    "formal analysis/evidence-family BH; restricted PERMANOVA q; no manual stars",
    "no new p/q; inherits formal source-table families and categorical tier rules"
  ),
  reviewer_risk = c(
    "shared patients; three seeds; positive-only main display; full audit is supplemental",
    "shared RNA representation; continuous-state and subtype overclaim",
    "small n; candidate prescreen; missing protein; exact-event overclaim",
    "unverified metabolite labels; no public blank; unequal depth; 16S function",
    "causal-looking network; duplicated evidence; T4 overclaim; compact selection"
  ),
  stringsAsFactors = FALSE
)

figure_formats <- c("svg", "pdf", "tiff", "png")
formal_figure_relative_paths <- unlist(lapply(
  figure_specs$file_stem,
  function(stem) file.path("figures", paste0(stem, ".", figure_formats))
), use.names = FALSE)
figure_manifest_relative_path <-
  file.path("results", "escc_multiomics_figure_artifact_manifest.tsv")

if (fields_only) {
  cat(paste(formal_figure_relative_paths, collapse = "\n"), "\n", sep = "")
  cat(figure_manifest_relative_path, "\n", sep = "")
  cat(
    "backend\tR-only\n",
    "raster_contract\tTIFF and PNG at 600 dpi\n",
    "manifest_publish_order\t20 figure artifacts first; manifest last\n",
    "qa_status_contract\tstructural_status verified; visual_qa_status pending_reopened_review\n",
    "visual_qa_finalize\t--finalize-visual-qa=<qa.md>; requires explicit Figure1–Figure5 PASS\n",
    "hard_boundaries\tofficial78 primary; Factor4 T0; Cao free scales; ",
    "PR001876/PRJNA766558 disconnected; no causal arrows\n",
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
    "R-only 图件依赖缺失：", paste(missing_packages, collapse = ", "),
    "。禁止切换 Python 代替渲染。", call. = FALSE
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
executed_script_sha256 <- digest(
  script_path, algo = "sha256", file = TRUE, serialize = FALSE
)

fail_if <- function(condition, message) {
  if (!is.logical(condition) || length(condition) != 1L || is.na(condition)) {
    stop(paste0(message, "（门禁为 NA 或非单一逻辑值，按失败处理）"),
         call. = FALSE)
  }
  if (condition) stop(message, call. = FALSE)
}

require_columns <- function(object, columns, object_name) {
  missing <- setdiff(columns, names(object))
  fail_if(length(missing) > 0L, paste0(
    object_name, " 缺少字段：", paste(missing, collapse = ", ")
  ))
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
  fail_if(any(!input_na & is.na(output)) || (!allow_na && anyNA(output)),
          paste0(label, " 含不可识别值或不允许的 NA。"))
  output
}

read_tsv <- function(path) {
  fail_if(!file_exists(path), paste("正式源表缺失：", path))
  fread(path, na.strings = c("", "NA"), showProgress = FALSE)
}

sha256_file <- function(path) {
  fail_if(!file_exists(path), paste("待校验文件缺失：", path))
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
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
    "视觉 QA 必须是 _work/checks/ 内的 Markdown 文件。"
  )

  qa_lines <- readLines(qa_abs, warn = FALSE, encoding = "UTF-8")
  normalized <- tolower(trimws(gsub("：", ":", qa_lines, fixed = TRUE)))
  fail_if(
    any(grepl(
      "figure\\s*[1-5].*(fail|failed|失败|不通过)", normalized,
      perl = TRUE
    )),
    "视觉 QA 仍含 Figure1–Figure5 的 FAIL/不通过记录，禁止 finalize。"
  )
  pass_pattern <-
    "^[-*]?[[:space:]]*figure[[:space:]]*([1-5])[[:space:]]*:[[:space:]]*(pass|passed|通过)[[:space:]]*$"
  pass_matches <- regmatches(
    normalized, regexec(pass_pattern, normalized, perl = TRUE)
  )
  passed_ids <- as.integer(vapply(
    pass_matches[lengths(pass_matches) > 0L], `[[`, character(1), 2L
  ))
  fail_if(
    length(passed_ids) != 5L ||
      !identical(sort(passed_ids), 1:5),
    paste(
      "视觉 QA 必须各含且只含一条 Figure1: PASS 至",
      "Figure5: PASS（也接受中文“通过”）。"
    )
  )

  manifest_path <- file.path(project_root, figure_manifest_relative_path)
  lock_dir <- file.path(work_intermediate_dir, ".escc_figure_publish.lock")
  lock_acquired <- FALSE
  lock_release_allowed <- FALSE
  suspendInterrupts({
    lock_acquired <- dir.create(
      lock_dir, recursive = FALSE, showWarnings = FALSE
    )
    fail_if(!lock_acquired,
            paste("图件发布/QA finalize 锁已存在：", lock_dir))
    lock_release_allowed <- TRUE
    on.exit(suspendInterrupts({
      if (lock_acquired && lock_release_allowed && dir_exists(lock_dir)) {
        try(dir_delete(lock_dir), silent = TRUE)
      }
    }), add = TRUE, after = FALSE)
  })

  verify_current_manifest <- function() {
    fail_if(!file_exists(manifest_path),
            paste("figure artifact manifest 缺失：", manifest_path))
    manifest <- read_tsv(manifest_path)
    require_columns(
      manifest,
      c(
        "figure_id", "relative_path", "file_size_bytes", "sha256",
        "structural_status", "visual_qa_status"
      ),
      "figure artifact manifest"
    )
    fail_if(nrow(manifest) != 20L ||
              uniqueN(manifest$relative_path) != 20L ||
              !setequal(manifest$relative_path, formal_figure_relative_paths) ||
              !identical(sort(unique(manifest$figure_id)), paste0("Figure", 1:5)) ||
              any(manifest[, .N, by = figure_id]$N != 4L),
            "figure artifact manifest 不是完整 5×4 family。")
    fail_if(any(manifest$structural_status !=
                  "verified_after_staged_export_and_hash_check"),
            "figure artifact manifest 存在未通过的结构状态。")
    artifact_paths <- file.path(project_root, manifest$relative_path)
    fail_if(any(!file_exists(artifact_paths)),
            "figure artifact manifest 登记的 20 个图件不完整。")
    actual_size <- as.numeric(file_info(artifact_paths)$size)
    actual_sha <- vapply(artifact_paths, sha256_file, character(1))
    fail_if(any(actual_size != as.numeric(manifest$file_size_bytes)) ||
              any(actual_sha != manifest$sha256),
            "当前 20 个图件的大小/SHA256 与 manifest 不一致。")
    manifest
  }

  manifest <- verify_current_manifest()
  qa_relative <- as.character(path_rel(qa_abs, start = project_root))
  qa_sha <- sha256_file(qa_abs)
  # `fread()` 将全空的 QA 列推断为 logical；先显式转为 character，
  # 避免 finalize 时把路径和 SHA256 强制转换成 NA。
  manifest[, `:=`(
    qa_path = as.character(qa_path),
    qa_sha256 = as.character(qa_sha256)
  )]
  manifest[, `:=`(
    visual_qa_status = "passed_reopened_review",
    qa_path = qa_relative,
    qa_sha256 = qa_sha
  )]

  stage_name <- paste0(
    ".", basename(manifest_path), ".visualqa-", Sys.getpid(), "-",
    format(Sys.time(), "%Y%m%d%H%M%OS6")
  )
  stage_name <- gsub("[^A-Za-z0-9_.-]", "", stage_name)
  stage_path <- file.path(results_dir, stage_name)
  fail_if(file_exists(stage_path), paste("QA manifest 暂存路径已存在：", stage_path))
  stage_active <- TRUE
  on.exit({
    if (stage_active && file_exists(stage_path)) try(file_delete(stage_path), silent = TRUE)
  }, add = TRUE, after = FALSE)
  fwrite(
    manifest, stage_path, sep = "\t", quote = FALSE, na = "",
    logical01 = FALSE
  )
  staged <- read_tsv(stage_path)
  require_columns(staged, c("visual_qa_status", "qa_path", "qa_sha256"),
                  "staged visual QA manifest")
  fail_if(nrow(staged) != 20L ||
            any(staged$visual_qa_status != "passed_reopened_review") ||
            any(staged$qa_path != qa_relative) ||
            any(staged$qa_sha256 != qa_sha),
          "视觉 QA manifest 暂存内容未通过回读。")
  moved <- file.rename(stage_path, manifest_path)
  fail_if(!moved, "视觉 QA manifest 原子替换失败。")
  stage_active <- FALSE

  finalized <- verify_current_manifest()
  require_columns(finalized, c("qa_path", "qa_sha256"),
                  "finalized visual QA manifest")
  fail_if(any(finalized$visual_qa_status != "passed_reopened_review") ||
            any(finalized$qa_path != qa_relative) ||
            any(finalized$qa_sha256 != qa_sha),
          "视觉 QA manifest 原子替换后验证失败。")
  lock_release_allowed <- TRUE
  if (dir_exists(lock_dir)) dir_delete(lock_dir)
  fail_if(dir_exists(lock_dir), "视觉 QA finalize 后发布锁未清理。")
  lock_acquired <- FALSE
  cat("FINALIZE_VISUAL_QA_OK\n")
  cat("qa_path\t", qa_relative, "\n", sep = "")
  cat("qa_sha256\t", qa_sha, "\n", sep = "")
  invisible(finalized)
}

if (finalize_visual_qa) {
  finalize_visual_qa_manifest(finalize_visual_qa_input)
  quit(save = "no", status = 0L)
}

figure_contract_path <- file.path(
  project_root, "_work", "checks",
  "escc_multiomics_final_figure_contract_20260712.md"
)
expected_figure_contract_sha256 <-
  "10c968e11a997cef96a5ae71ecdf95cdf7887f749eb80835523eb7c2bce05a0f"
fail_if(!file_exists(figure_contract_path),
        "冻结的 ESCC 最终图件契约缺失。")
figure_contract_sha256 <- sha256_file(figure_contract_path)
fail_if(figure_contract_sha256 != expected_figure_contract_sha256,
        "冻结图件契约 SHA256 已漂移；需先显式审查并更新脚本绑定。")

artifact_family_exact <- function(expected, observed) {
  identical(sort(unique(expected)), sort(unique(observed))) &&
    length(expected) == length(observed) && !anyDuplicated(observed)
}

verify_artifact_family <- function(
    manifest_name,
    formal_names,
    generation_script,
    allowed_status,
    generator_sha_field = NULL,
    exact_current_generator = FALSE,
    family_regex = NULL) {
  manifest_path <- file.path(results_dir, manifest_name)
  fail_if(!file_exists(manifest_path),
          paste("artifact manifest 缺失：", manifest_name))
  manifest <- read_tsv(manifest_path)
  required <- c(
    "relative_path", "file_size_bytes", "sha256", "status",
    "generation_script"
  )
  if (exact_current_generator) required <- c(required, generator_sha_field)
  require_columns(manifest, required, manifest_name)
  fail_if(anyDuplicated(manifest$relative_path) > 0L,
          paste("manifest relative_path 重复：", manifest_name))
  observed_names <- basename(manifest$relative_path)
  fail_if(!artifact_family_exact(formal_names, observed_names), paste(
    "正式 artifact family 与预锁定集合不一致：", manifest_name,
    "；expected=", paste(formal_names, collapse = ";"),
    "；observed=", paste(observed_names, collapse = ";")
  ))
  fail_if(any(manifest$generation_script != generation_script),
          paste("manifest generation_script 不一致：", manifest_name))
  fail_if(any(!manifest$status %chin% allowed_status),
          paste("manifest 存在未允许状态：", manifest_name))

  artifact_paths <- file.path(project_root, manifest$relative_path)
  fail_if(any(!file_exists(artifact_paths)),
          paste("manifest 登记 artifact 缺失：", manifest_name))
  actual_size <- as.numeric(file_info(artifact_paths)$size)
  actual_sha <- vapply(artifact_paths, sha256_file, character(1))
  fail_if(
    any(actual_size != as.numeric(manifest$file_size_bytes)) ||
      any(actual_sha != manifest$sha256),
    paste("artifact 大小或 SHA256 与 manifest 不一致：", manifest_name)
  )

  if (exact_current_generator) {
    fail_if(is.null(generator_sha_field) || !nzchar(generator_sha_field),
            paste("未指定生成脚本 SHA 字段：", manifest_name))
    observed_generator_sha <- unique(as.character(manifest[[generator_sha_field]]))
    fail_if(length(observed_generator_sha) != 1L ||
              is.na(observed_generator_sha) ||
              !nzchar(observed_generator_sha),
            paste("manifest 生成脚本 SHA 不唯一：", manifest_name))
    generator_path <- file.path(project_root, generation_script)
    fail_if(!file_exists(generator_path),
            paste("当前生成脚本缺失：", generation_script))
    fail_if(sha256_file(generator_path) != observed_generator_sha, paste(
      "当前上游脚本 SHA256 与正式 manifest 不一致：", manifest_name,
      "。必须先由当前脚本重新发布完整 family。"
    ))
    fail_if(is.null(family_regex) || !nzchar(family_regex),
            paste("exact family 缺少 results 文件族正则：", manifest_name))
    results_members <- basename(dir_ls(results_dir, type = "file", fail = FALSE))
    results_members <- results_members[grepl(family_regex, results_members)]
    fail_if(
      !artifact_family_exact(c(formal_names, manifest_name), results_members),
      paste(
        "results/ 中 exact artifact family 有缺失或 manifest 外额外文件：",
        manifest_name, "；observed=", paste(results_members, collapse = ";")
      )
    )
  }

  list(
    manifest = manifest,
    manifest_path = manifest_path,
    manifest_sha256 = sha256_file(manifest_path),
    artifact_paths = artifact_paths,
    artifact_sha256 = actual_sha
  )
}

family_contracts <- list(
  leakage = list(
    manifest = "tcga_escc_driver_state_leakage_artifact_manifest.tsv",
    files = c(
      "tcga_escc_driver_state_leakage_model_plan.tsv",
      "tcga_escc_driver_state_leakage_factor_alignment.tsv",
      "tcga_escc_driver_state_leakage_associations.tsv",
      "tcga_escc_driver_state_leakage_edge_stability.tsv",
      "tcga_escc_driver_state_leakage_candidate_summary.tsv",
      "tcga_escc_driver_state_external_projection_status.tsv",
      "tcga_escc_driver_state_leakage_summary.md"
    ),
    script = "scripts/21_audit_tcga_escc_driver_state_leakage.R",
    status = "verified",
    sha_field = "generation_script_sha256",
    exact = TRUE,
    regex = paste0(
      "^tcga_escc_driver_state_(leakage_.*|external_projection_status[.]tsv)$"
    )
  ),
  microbiome = list(
    manifest = "prjna766558_dada2_artifact_manifest.tsv",
    files = c(
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
    ),
    script = "scripts/18_analyze_prjna766558_dada2.R",
    status = "verified",
    sha_field = "executed_script_sha256",
    exact = TRUE,
    regex = "^prjna766558_dada2_.*"
  ),
  integration = list(
    manifest = "escc_multiomics_integration_artifact_manifest.tsv",
    files = c(
      "escc_multiomics_integration_tier_definitions.tsv",
      "escc_multiomics_evidence_ledger.tsv",
      "escc_multiomics_integrated_driver_candidates.tsv",
      "escc_multiomics_integrated_axis_edges.tsv",
      "escc_multiomics_integrated_axis_summary.tsv",
      "escc_multiomics_heterogeneity_axes.tsv",
      "escc_multiomics_module_summaries.tsv",
      "escc_multiomics_integration_summary.md"
    ),
    script = "scripts/20_integrate_escc_multiomics_evidence.R",
    status = "verified",
    sha_field = "executed_script_sha256",
    exact = TRUE,
    regex = paste0(
      "^escc_multiomics_(integration_tier_definitions|evidence_ledger|",
      "integrated_driver_candidates|integrated_axis_edges|integrated_axis_summary|",
      "heterogeneity_axes|module_summaries|integration_summary|",
      "integration_artifact_manifest)[.]"
    )
  ),
  ecms = list(
    manifest = "tcga_escc_ecms_projection_artifact_manifest.tsv",
    files = c(
      "tcga_escc_ecms_patient_probabilities.tsv",
      "tcga_escc_ecms_projection_calibration.tsv",
      "tcga_escc_ecms_projection_qa.tsv",
      "tcga_escc_ecms_factor_associations.tsv",
      "tcga_escc_ecms_adjusted_factor_progeny_associations.tsv",
      "tcga_escc_ecms_projection_summary.md"
    ),
    script = "scripts/23_project_tcga_escc_ecms.R",
    status = "verified",
    sha_field = "execution_script_sha256",
    exact = TRUE,
    regex = "^tcga_escc_ecms_.*"
  ),
  heterogeneity = list(
    manifest = "tcga_escc_heterogeneity_artifact_manifest.tsv",
    files = c(
      "tcga_escc_heterogeneity_cluster_evaluation.tsv",
      "tcga_escc_heterogeneity_cluster_pathways.tsv",
      "tcga_escc_heterogeneity_feature_manifest.tsv",
      "tcga_escc_heterogeneity_patient_assignments.tsv",
      "tcga_escc_heterogeneity_summary.md",
      "tcga_escc_mofa_clinical_associations.tsv",
      "tcga_escc_mofa_factor_scores.tsv",
      "tcga_escc_mofa_level_factor_qc.tsv",
      "tcga_escc_mofa_model.rds",
      "tcga_escc_mofa_pathway_associations.tsv",
      "tcga_escc_mofa_top_weights.tsv",
      "tcga_escc_mofa_variance_explained.tsv",
      "tcga_escc_progeny_pathway_scores.tsv"
    ),
    script = "scripts/14_analyze_tcga_escc_heterogeneity.R",
    status = "verified",
    sha_field = NULL,
    exact = FALSE,
    regex = NULL
  ),
  cao = list(
    manifest = "cao2020_cross_layer_artifact_manifest.tsv",
    files = c(
      "cao2020_wgbs_candidate_pair_effects.tsv",
      "cao2020_wgbs_candidate_region_summary.tsv",
      "cao2020_rna_candidate_pair_effects.tsv",
      "cao2020_rna_candidate_summary.tsv",
      "cao2020_proteomics_candidate_pair_effects.tsv",
      "cao2020_proteomics_candidate_summary.tsv",
      "cao2020_cross_layer_patient_effects.tsv",
      "cao2020_cross_layer_candidate_summary.tsv",
      "cao2020_cross_layer_summary.md"
    ),
    script = "scripts/15_analyze_cao2020_cross_layer_axes.R",
    status = "verified",
    sha_field = NULL,
    exact = FALSE,
    regex = NULL
  ),
  metabolome = list(
    manifest = "pr001876_targeted_ms_artifact_manifest.tsv",
    files = c(
      "pr001876_targeted_ms_analysis_inventory.tsv",
      "pr001876_targeted_ms_sample_qc.tsv",
      "pr001876_targeted_ms_feature_qc.tsv",
      "pr001876_targeted_ms_differential.tsv",
      "pr001876_targeted_ms_candidate_metabolites.tsv",
      "pr001876_targeted_ms_run_order_sensitivity.tsv",
      "pr001876_targeted_ms_scale_sensitivity.tsv",
      "pr001876_targeted_ms_paired_sensitivity.tsv",
      "pr001876_targeted_ms_summary.md"
    ),
    script = "scripts/16_analyze_pr001876_metabolomics.R",
    status = "verified_after_atomic_publish",
    sha_field = NULL,
    exact = FALSE,
    regex = NULL
  )
)

verify_all_families <- function() {
  fail_if(!file_exists(file.path(project_root, "PROJECT_INDEX.md")),
          "PROJECT_INDEX.md 缺失。")
  fail_if(!dir_exists(results_dir) || !dir_exists(figures_dir) ||
            !dir_exists(work_intermediate_dir) || !dir_exists(work_checks_dir),
          "项目规范目录缺失；禁止静默改用其他输出目录。")
  output <- lapply(family_contracts, function(contract) {
    verify_artifact_family(
      manifest_name = contract$manifest,
      formal_names = contract$files,
      generation_script = contract$script,
      allowed_status = contract$status,
      generator_sha_field = contract$sha_field,
      exact_current_generator = contract$exact,
      family_regex = contract$regex
    )
  })
  names(output) <- names(family_contracts)
  output
}

if (self_test) {
  fail_if(!identical(length(formal_figure_relative_paths), 20L),
          "应有 5×4=20 个正式图件输出。")
  fail_if(anyDuplicated(formal_figure_relative_paths) > 0L,
          "正式图件路径重复。")
  fail_if(any(figure_specs$width_mm != 183), "主图宽度必须为 183 mm。")
  fail_if(!identical(figure_formats, c("svg", "pdf", "tiff", "png")),
          "导出格式契约漂移。")
  selftest_font <- systemfonts::match_fonts("Arial")
  fail_if(is.null(selftest_font$path) || !nzchar(selftest_font$path) ||
            !file_exists(selftest_font$path),
          "SELF TEST 未找到 Arial 字体。")
  script_text <- paste(readLines(script_path, warn = FALSE), collapse = "\n")
  required_contract_tokens <- c(
    paste0("publish_figure_bundle_", "transaction <- function"),
    paste0("rollback_", "verified"),
    paste0("manifest literal", "-last"),
    paste0("visual_qa_status = ", "\"pending_reopened_review\""),
    paste0("metabolite_analysis_", "inventory"),
    paste0("microbiome_sample_", "qc"),
    paste0("microbiome_", "ancombc2"),
    paste0("Patient-level tumour-minus-paired-non-tumour ", "differences"),
    paste0("effect_", "direction"),
    paste0("FFPE low ", "biomass; contamination unresolved"),
    paste0("reproducibility-only ", "subset"),
    paste0("--finalize-visual-", "qa="),
    paste0("passed_reopened_", "review"),
    paste0("qa_", "sha256"),
    paste0("build_figure1_", "positive"),
    paste0("build_figure2_", "positive"),
    paste0("build_figure3_", "positive"),
    paste0("build_figure4_", "positive"),
    paste0("build_figure5_", "positive"),
    paste0("Representation-overlap-audited ", "event–state bridges"),
    paste0("83 robust LC-MS ", "analysis-features"),
    paste0("PIK3CA leads the seven-gene T3/T4 ", "validation set")
  )
  fail_if(any(!vapply(
    required_contract_tokens,
    grepl,
    logical(1),
    x = script_text,
    fixed = TRUE
  )), "SELF TEST 缺少 bundle/边界/视觉QA 静态契约。")
  fail_if(grepl(paste0("atomic_publish_", "file <- function"),
                script_text, fixed = TRUE) ||
            grepl(paste0("six 16-vs-16 ", "PR001876 analyses"),
                  script_text, fixed = TRUE) ||
            grepl(paste0("tcga_escc_driver_state_", "artifact_manifest.tsv"), script_text,
                  fixed = TRUE),
          paste(
            "SELF TEST 发现旧逐文件发布、错误六分析或未使用的直接",
            "driver-state manifest 契约残留。"
          ))
  cat("SELF_TEST_OK\n")
  quit(save = "no", status = 0L)
}

verified_families <- verify_all_families()

input_paths <- list(
  leakage_model_plan = file.path(
    results_dir, "tcga_escc_driver_state_leakage_model_plan.tsv"
  ),
  leakage_factor_alignment = file.path(
    results_dir, "tcga_escc_driver_state_leakage_factor_alignment.tsv"
  ),
  leakage_associations = file.path(
    results_dir, "tcga_escc_driver_state_leakage_associations.tsv"
  ),
  leakage_edge_stability = file.path(
    results_dir, "tcga_escc_driver_state_leakage_edge_stability.tsv"
  ),
  leakage_candidate_summary = file.path(
    results_dir, "tcga_escc_driver_state_leakage_candidate_summary.tsv"
  ),
  ecms_probabilities = file.path(
    results_dir, "tcga_escc_ecms_patient_probabilities.tsv"
  ),
  ecms_calibration = file.path(
    results_dir, "tcga_escc_ecms_projection_calibration.tsv"
  ),
  ecms_factor = file.path(
    results_dir, "tcga_escc_ecms_factor_associations.tsv"
  ),
  ecms_adjusted = file.path(
    results_dir, "tcga_escc_ecms_adjusted_factor_progeny_associations.tsv"
  ),
  factor_scores = file.path(results_dir, "tcga_escc_mofa_factor_scores.tsv"),
  level_factor_qc = file.path(results_dir, "tcga_escc_mofa_level_factor_qc.tsv"),
  cluster_evaluation = file.path(
    results_dir, "tcga_escc_heterogeneity_cluster_evaluation.tsv"
  ),
  cao_summary = file.path(results_dir, "cao2020_cross_layer_candidate_summary.tsv"),
  cao_patient = file.path(results_dir, "cao2020_cross_layer_patient_effects.tsv"),
  metabolite_candidates = file.path(
    results_dir, "pr001876_targeted_ms_candidate_metabolites.tsv"
  ),
  metabolite_analysis_inventory = file.path(
    results_dir, "pr001876_targeted_ms_analysis_inventory.tsv"
  ),
  microbiome_sample_qc = file.path(
    results_dir, "prjna766558_dada2_sample_qc.tsv"
  ),
  microbiome_pipeline = file.path(
    results_dir, "prjna766558_dada2_pipeline_sensitivity.tsv"
  ),
  microbiome_alpha = file.path(
    results_dir, "prjna766558_dada2_alpha_paired_tests.tsv"
  ),
  microbiome_beta = file.path(
    results_dir, "prjna766558_dada2_beta_tests.tsv"
  ),
  microbiome_genus = file.path(
    results_dir, "prjna766558_dada2_genus_paired_differential.tsv"
  ),
  microbiome_ancombc2 = file.path(
    results_dir, "prjna766558_dada2_ancombc2_sensitivity.tsv"
  ),
  integrated_drivers = file.path(
    results_dir, "escc_multiomics_integrated_driver_candidates.tsv"
  ),
  integrated_axis_edges = file.path(
    results_dir, "escc_multiomics_integrated_axis_edges.tsv"
  ),
  integrated_axis_summary = file.path(
    results_dir, "escc_multiomics_integrated_axis_summary.tsv"
  ),
  heterogeneity_axes = file.path(
    results_dir, "escc_multiomics_heterogeneity_axes.tsv"
  ),
  module_summaries = file.path(
    results_dir, "escc_multiomics_module_summaries.tsv"
  )
)

inputs <- lapply(input_paths, read_tsv)

schema_contract <- list(
  leakage_model_plan = c(
    "model_id", "seed", "scenario_order", "scenario_id", "scenario_class",
    "scenario_no_op", "expected_model_action", "run_status", "model_reused",
    "model_completed", "model_converged", "checkpoint_path"
  ),
  leakage_factor_alignment = c(
    "model_id", "seed", "scenario_id", "scenario_order", "reference_factor",
    "candidate_factor", "absolute_shared_loading_correlation",
    "minimum_match_margin", "low_margin_flag", "match_reliable",
    "matching_status", "patient_score_correlation_diagnostic_aligned",
    "patient_score_diagnostic_is_gate", "loading_alignment_is_event_blind"
  ),
  leakage_associations = c(
    "model_id", "seed", "scenario_id", "scenario_class", "scenario_no_op",
    "original_edge_id", "gene_name", "event_type", "reference_factor",
    "aligned_effect", "q_value", "direction_retained",
    "magnitude_retained_50pct", "gate_retained", "counts_toward_retention"
  ),
  leakage_edge_stability = c(
    "original_edge_id", "gene_name", "event_type", "reference_factor",
    "scenario_id", "scenario_class", "planned_seed_count",
    "reliable_factor_match_rate", "direction_retention_rate",
    "magnitude_retention_rate", "gate_retention_rate", "low_margin_count",
    "scenario_no_op", "retention_decision"
  ),
  leakage_candidate_summary = c(
    "original_edge_id", "gene_name", "narrative_role", "event_type",
    "reference_factor", "original_effect", "original_association_status",
    "original_level_factor_soft_flag", "planned_seed_count",
    "completed_converged_seed_count", "baseline_reliable_match_rate",
    "drop_both_reliable_match_rate", "drop_both_direction_retention_rate",
    "drop_both_magnitude_retention_rate", "drop_both_gate_retention_rate",
    "drop_both_low_margin_count", "leakage_gate_evaluable",
    "leakage_gate_pass", "countable_for_T3_T4", "gate_failure_reason"
  ),
  ecms_probabilities = c(
    "patient_id", "in_official_78", "official_anchor_label",
    "extension_status", "eligible_for_primary_association",
    "resolved_ecms_label", "label_source", "low_margin_custom_flag",
    "single_sample_classifier_claim", "pseudo_label_generated"
  ),
  ecms_calibration = c(
    "metric", "observed", "threshold", "direction", "gate_role", "pass"
  ),
  ecms_factor = c(
    "analysis_scope", "association_role", "is_primary_scope", "n_patients",
    "factor", "ecms1_n", "ecms2_n", "ecms3_n", "ecms4_n",
    "ecms1_mean", "ecms2_mean", "ecms3_mean", "ecms4_mean",
    "eta_squared", "anova_q_value", "kruskal_q_value",
    "extension_calibration_pass", "shared_rna_representation",
    "independent_validation", "evidence_ceiling"
  ),
  ecms_adjusted = c(
    "analysis_scope", "association_role", "is_primary_scope", "n_patients",
    "pathway", "factor", "factor_beta_standardized", "partial_r_squared",
    "incremental_p_value", "incremental_q_value",
    "extension_calibration_pass", "shared_rna_representation",
    "independent_validation", "evidence_ceiling"
  ),
  factor_scores = c("patient_id", paste0("Factor", 1:8)),
  level_factor_qc = c(
    "factor", "view", "spearman_rho_with_view_mean", "q_value",
    "level_factor_soft_flag"
  ),
  cluster_evaluation = c(
    "k", "snf_mean_silhouette", "mofa_consensus_pac",
    "snf_mofa_adjusted_rand", "conditional_stability",
    "selected_for_description", "decision"
  ),
  cao_summary = c(
    "tcga_gene_name", "promoter_median_paired_delta_beta",
    "promoter_region_evaluation_status", "promoter_q_candidate_by_region",
    "limma_log_fc_pc1", "limma_q_candidate_pc1", "rna_effect_class",
    "protein_measurement_status", "q_candidate",
    "protein_effect_class", "median_patient_log2_tumor_vs_normal",
    "cross_layer_class", "statistical_support_class", "available_layer_count",
    "fdr_supported_layer_count"
  ),
  cao_patient = c(
    "tcga_gene_name", "paper_patient_id", "rna_delta_log2_tpm_pc1",
    "median_delta_beta_promoter", "coverage_gate_promoter",
    "protein_log2_tumor_vs_normal", "protein_patient_measurement_status"
  ),
  metabolite_candidates = c(
    "analysis_id", "platform", "primary_transform", "feature_id",
    "metabolite_name", "annotation_identity_status", "limma_effect_transformed",
    "limma_ci95_low", "limma_ci95_high", "limma_q_evidence_family",
    "wilcoxon_direction", "direction_concordant",
    "scale_sensitivity_direction_stable",
    "run_order_perfect_group_block",
    "run_order_sensitivity_direction_concordant",
    "paired_sensitivity_estimable", "paired_sensitivity_direction_concordant",
    "final_candidate_tier", "candidate_rank",
    "host_multiomics_sample_level_link_allowed"
  ),
  metabolite_analysis_inventory = c(
    "analysis_id", "study_id", "platform", "ion_mode", "transform",
    "n_escc", "n_normal", "candidates_total", "primary_design",
    "within_analysis_evidence_unit"
  ),
  microbiome_sample_qc = c(
    "patient_pair_id", "paper_pair_number", "sample_name", "tissue_role", "ffpe_status",
    "public_negative_control_run", "blank_control_status"
  ),
  microbiome_pipeline = c(
    "pipeline", "final_read_retention", "bray_distance_spearman_vs_main",
    "differential_direction_concordance_vs_main"
  ),
  microbiome_alpha = c(
    "metric", "n_pairs", "median_paired_difference_tumor_minus_normal",
    "effect_ci95_low", "effect_ci95_high", "paired_wilcoxon_p",
    "paired_wilcoxon_q", "interpretation_scope"
  ),
  microbiome_beta = c(
    "distance", "n_pairs", "permanova_r2", "restricted_permutation_p",
    "restricted_permutation_q", "betadisper_p",
    "paired_distance_to_centroid_p"
  ),
  microbiome_genus = c(
    "genus_label", "n_pairs", "median_clr_difference_tumor_minus_normal",
    "clr_effect_ci95_low", "clr_effect_ci95_high", "paired_wilcoxon_p",
    "paired_wilcoxon_q", "bootstrap_direction_consistency",
    "candidate_tier", "blank_control_status", "contamination_risk",
    "function_claim_scope", "host_sample_level_link_allowed"
  ),
  microbiome_ancombc2 = c(
    "taxon", "analysis_status", "error_message"
  ),
  integrated_drivers = c(
    "candidate_id", "gene_name", "decision", "raw_state_supported_edges",
    "state_supported_edges", "cao_cross_layer_class", "cao_statistical_support",
    "integrated_tier", "integrated_decision", "high_priority_validation_rule",
    "cross_cohort_patient_link_created", "exact_independent_validation_present"
  ),
  integrated_axis_edges = c(
    "edge_id", "axis_candidate_id", "source_node", "source_layer",
    "target_node", "target_layer", "edge_class", "association_effect",
    "association_effect_measure", "q_value", "evidence_tier", "support_status",
    "representation_overlap", "evidence_source", "independence_group",
    "independent_from_tcga_discovery", "conclusion_ceiling", "source_file",
    "source_row_key", "countable_as_exact_driver_validation"
  ),
  integrated_axis_summary = c(
    "gene_name", "candidate_id", "tcga_driver_tier", "state_supported_edges",
    "state_conditional_edges", "cao_cross_layer_class", "axis_class",
    "axis_tier", "exception_or_boundary", "cross_cohort_patient_link_created",
    "exact_independent_validation_present"
  ),
  heterogeneity_axes = c(
    "factor", "level_factor_soft_flag", "raw_supported_driver_edge_count",
    "supported_driver_edge_count", "conditional_driver_edge_count",
    "ecms_context_tier", "ecms_context_status", "ecms_eta_squared",
    "ecms_context_q", "ecms_adjusted_pathway_t2_count",
    "top_pathway", "top_pathway_beta", "top_pathway_partial_r_squared",
    "discrete_subtype_support", "evidence_tier", "cluster_decision"
  ),
  module_summaries = c(
    "module_id", "dataset", "design", "formal_analysis_available",
    "host_sample_level_link_to_other_modules_allowed",
    "exact_driver_event_independent_validation", "module_evidence_tier",
    "module_result", "conclusion_ceiling"
  )
)

for (input_name in names(schema_contract)) {
  require_columns(inputs[[input_name]], schema_contract[[input_name]], input_name)
}

validate_semantics <- function(x, families) {
  plan <- x$leakage_model_plan
  plan[, `:=`(
    scenario_no_op = as_logical_strict(scenario_no_op, "plan scenario_no_op"),
    model_reused = as_logical_strict(model_reused, "plan model_reused"),
    model_completed = as_logical_strict(model_completed, "plan model_completed"),
    model_converged = as_logical_strict(model_converged, "plan model_converged")
  )]
  fail_if(uniqueN(plan$seed) != 3L || nrow(plan) != 24L,
          "script21 模型计划必须为 3 seeds × 8 场景。")
  fail_if(any(!plan$model_completed) || any(!plan$model_converged),
          "script21 存在未完成或未收敛模型，禁止绘图。")
  fail_if(any(plan$run_status == "planned_not_run"),
          "script21 model plan 仍含 planned_not_run。")

  alignment <- x$leakage_factor_alignment
  alignment[, `:=`(
    low_margin_flag = as_logical_strict(
      low_margin_flag, "alignment low_margin_flag"
    ),
    match_reliable = as_logical_strict(match_reliable, "alignment match_reliable"),
    patient_score_diagnostic_is_gate = as_logical_strict(
      patient_score_diagnostic_is_gate, "alignment diagnostic gate"
    ),
    loading_alignment_is_event_blind = as_logical_strict(
      loading_alignment_is_event_blind, "alignment event blind"
    )
  )]
  fail_if(any(alignment$patient_score_diagnostic_is_gate) ||
            any(!alignment$loading_alignment_is_event_blind),
          "因子对齐使用了患者得分门禁或非 event-blind 载荷。")
  fail_if(any(alignment$match_reliable & alignment$low_margin_flag),
          "low-margin 匹配被标为可靠。")

  associations <- x$leakage_associations
  for (column in c(
    "scenario_no_op", "direction_retained", "magnitude_retained_50pct",
    "gate_retained", "counts_toward_retention"
  )) {
    set(associations, j = column, value = as_logical_strict(
      associations[[column]], paste("leakage associations", column),
      allow_na = column %chin% c(
        "direction_retained", "magnitude_retained_50pct", "gate_retained"
      )
    ))
  }
  fail_if(any(associations$counts_toward_retention & associations$scenario_no_op),
          "no-op 被错误计入关联保留率。")

  stability <- x$leakage_edge_stability
  stability[, scenario_no_op := as_logical_strict(
    scenario_no_op, "edge stability scenario_no_op"
  )]

  leakage <- x$leakage_candidate_summary
  logical_cols <- c(
    "original_level_factor_soft_flag", "leakage_gate_evaluable",
    "leakage_gate_pass", "countable_for_T3_T4"
  )
  for (column in logical_cols) {
    set(leakage, j = column, value = as_logical_strict(
      leakage[[column]], paste("leakage", column)
    ))
  }
  fail_if(nrow(leakage) != 41L || uniqueN(leakage$original_edge_id) != 41L,
          "泄漏审计必须保留全部 41 条原始 MOFA 网络边。")
  fail_if(any(leakage$leakage_gate_pass != leakage$countable_for_T3_T4),
          "leakage_gate_pass 与 countable_for_T3_T4 不一致。")
  fail_if(any(leakage$countable_for_T3_T4 &
                (leakage$reference_factor == "Factor4" |
                   leakage$original_level_factor_soft_flag)),
          "Factor4/level-factor 被错误允许进入 T3/T4。")

  probability <- x$ecms_probabilities
  probability[, `:=`(
    in_official_78 = as_logical_strict(in_official_78, "ECMS official78"),
    eligible_for_primary_association = as_logical_strict(
      eligible_for_primary_association, "ECMS primary eligibility"
    ),
    low_margin_custom_flag = as_logical_strict(
      low_margin_custom_flag, "ECMS low margin"
    ),
    single_sample_classifier_claim = as_logical_strict(
      single_sample_classifier_claim, "ECMS single sample claim"
    ),
    pseudo_label_generated = as_logical_strict(
      pseudo_label_generated, "ECMS pseudo label"
    )
  )]
  official <- probability[in_official_78 & eligible_for_primary_association]
  fail_if(nrow(official) != 78L || uniqueN(official$patient_id) != 78L,
          "ECMS primary 必须且只能是官方 78 例。")
  fail_if(nrow(probability) != 94L ||
            nrow(probability[in_official_78 == FALSE]) != 16L ||
            any(probability[
              in_official_78 == FALSE, eligible_for_primary_association
            ]),
          "ECMS 投影必须为 official78 primary + 16 extension boundary。")
  fail_if(anyNA(official$official_anchor_label) ||
            any(!official$official_anchor_label %chin% paste0("ECMS", 1:4)),
          "官方 78 例的 ECMS anchor label 不完整。")
  fail_if(any(probability$single_sample_classifier_claim) ||
            any(probability$pseudo_label_generated),
          "ECMS 被越界标为单样本分类器或生成伪标签。")
  ecms_manifest <- families$ecms$manifest
  require_columns(
    ecms_manifest,
    c("extension_calibration_pass", "primary_patient_count"),
    "ECMS manifest"
  )
  extension_pass <- unique(as_logical_strict(
    ecms_manifest$extension_calibration_pass,
    "ECMS manifest extension calibration"
  ))
  primary_count <- unique(as.integer(ecms_manifest$primary_patient_count))
  fail_if(length(extension_pass) != 1L || extension_pass,
          "当前图件契约要求 extension calibration=FALSE。")
  fail_if(length(primary_count) != 1L || primary_count != 78L,
          "ECMS manifest primary_patient_count 必须为 78。")
  calibration <- x$ecms_calibration
  calibration[, pass := as_logical_strict(
    pass, "ECMS calibration pass", allow_na = TRUE
  )]
  fail_if(sum(calibration$gate_role == "prelocked_extension_gate") != 4L ||
            sum(calibration$pass[calibration$gate_role ==
                                   "prelocked_extension_gate"]) != 2L,
          "ECMS 预锁定扩展校准必须完整记录 2/4 通过。")

  factor_primary <- x$ecms_factor
  factor_primary[, `:=`(
    is_primary_scope = as_logical_strict(is_primary_scope, "ECMS factor primary"),
    extension_calibration_pass = as_logical_strict(
      extension_calibration_pass, "ECMS factor extension"
    ),
    shared_rna_representation = as_logical_strict(
      shared_rna_representation, "ECMS factor shared RNA"
    ),
    independent_validation = as_logical_strict(
      independent_validation, "ECMS factor independence"
    )
  )]
  fail_if(nrow(factor_primary[is_primary_scope == TRUE]) != 8L ||
            any(factor_primary[is_primary_scope == TRUE, n_patients] != 78L),
          "ECMS factor primary 表必须为 official78 × 8 factors。")
  fail_if(any(!factor_primary[
    is_primary_scope == TRUE, shared_rna_representation
  ]) || any(factor_primary[is_primary_scope == TRUE, independent_validation]),
          "ECMS factor 共享 RNA/非独立边界漂移。")

  adjusted_primary <- x$ecms_adjusted
  adjusted_primary[, `:=`(
    is_primary_scope = as_logical_strict(
      is_primary_scope, "ECMS adjusted primary"
    ),
    extension_calibration_pass = as_logical_strict(
      extension_calibration_pass, "ECMS adjusted extension"
    ),
    shared_rna_representation = as_logical_strict(
      shared_rna_representation, "ECMS adjusted shared RNA"
    ),
    independent_validation = as_logical_strict(
      independent_validation, "ECMS adjusted independence"
    )
  )]
  fail_if(!nrow(adjusted_primary[is_primary_scope == TRUE]) ||
            any(adjusted_primary[is_primary_scope == TRUE, n_patients] != 78L),
          "ECMS-adjusted primary 表缺失或不是 official78。")
  fail_if(any(!adjusted_primary[
    is_primary_scope == TRUE, shared_rna_representation
  ]) || any(adjusted_primary[
    is_primary_scope == TRUE, independent_validation
  ]),
          "ECMS-adjusted 通路共享 RNA/非独立边界漂移。")

  level <- x$level_factor_qc
  level[, level_factor_soft_flag := as_logical_strict(
    level_factor_soft_flag, "level-factor QC"
  )]
  fail_if(!any(level$factor == "Factor4" & level$level_factor_soft_flag),
          "Factor4 未被正式 level-factor QC 标记。")

  cluster <- x$cluster_evaluation
  cluster[, `:=`(
    conditional_stability = as_logical_strict(
      conditional_stability, "cluster conditional stability"
    ),
    selected_for_description = as_logical_strict(
      selected_for_description, "cluster selected"
    )
  )]
  fail_if(any(cluster$conditional_stability) ||
            sum(cluster$selected_for_description) != 1L,
          "离散聚类门禁或唯一描述性 k 与冻结结论不一致。")

  cao <- x$cao_summary
  fail_if(nrow(cao) != 12L || uniqueN(cao$tcga_gene_name) != 12L,
          "Cao 校准必须覆盖全部 12 个 strong drivers。")

  metabolites <- x$metabolite_candidates
  metabolites[, host_multiomics_sample_level_link_allowed :=
    as_logical_strict(
      host_multiomics_sample_level_link_allowed,
      "PR001876 host sample link"
    )]
  fail_if(any(metabolites$host_multiomics_sample_level_link_allowed),
          "PR001876 被错误允许建立宿主患者级连接。")

  metabolite_inventory <- x$metabolite_analysis_inventory
  expected_analyses <- c("AN004960", "AN004962", "AN004963")
  fail_if(nrow(metabolite_inventory) != 3L ||
            !setequal(metabolite_inventory$analysis_id, expected_analyses) ||
            any(metabolite_inventory$n_escc != 16L) ||
            any(metabolite_inventory$n_normal != 16L) ||
            sum(metabolite_inventory$n_escc + metabolite_inventory$n_normal) != 96L,
          paste(
            "PR001876 权威分析必须为 AN004960/AN004962/AN004963 三个主模型，",
            "各16 vs16，共96个 analysis-sample keys。"
          ))
  candidate_counts <- metabolites[, .N, by = analysis_id]
  candidate_counts <- merge(
    metabolite_inventory[, .(analysis_id, candidates_total)],
    candidate_counts, by = "analysis_id", all.x = TRUE, sort = FALSE
  )
  candidate_counts[is.na(N), N := 0L]
  fail_if(any(candidate_counts$N != candidate_counts$candidates_total) ||
            !setequal(unique(metabolites$analysis_id), expected_analyses),
          "PR001876 candidate 表未完整覆盖三个主模型或与 analysis_inventory 计数不一致。")

  microbe_sample_qc <- x$microbiome_sample_qc
  microbe_sample_qc[, paper_pair_number := as.integer(paper_pair_number)]
  published_pair_map <- unique(
    microbe_sample_qc[, .(patient_pair_id, paper_pair_number)]
  )
  fail_if(nrow(microbe_sample_qc) != 42L ||
            uniqueN(microbe_sample_qc$patient_pair_id) != 21L ||
            nrow(published_pair_map) != 21L ||
            anyNA(published_pair_map$paper_pair_number) ||
            anyDuplicated(published_pair_map$paper_pair_number) > 0L ||
            !identical(
              sort(published_pair_map$paper_pair_number),
              seq.int(
                min(published_pair_map$paper_pair_number),
                max(published_pair_map$paper_pair_number), by = 2L
              )
            ) || min(published_pair_map$paper_pair_number) != 1L ||
            any(!microbe_sample_qc$public_negative_control_run %chin%
                  c("not_identified", "none")) ||
            any(!grepl("ffpe", microbe_sample_qc$ffpe_status, ignore.case = TRUE)) ||
            any(microbe_sample_qc$blank_control_status !=
                  "no_public_negative_control_or_extraction_blank"),
          paste(
            "PRJNA766558 必须为21对/42个 FFPE 样本，且没有公开 blank/",
            "extraction negative control。"
          ))

  genus <- x$microbiome_genus
  genus[, host_sample_level_link_allowed := as_logical_strict(
    host_sample_level_link_allowed, "PRJNA766558 host sample link"
  )]
  fail_if(any(genus$host_sample_level_link_allowed),
          "PRJNA766558 被错误允许建立宿主患者级连接。")
  fail_if(any(genus$blank_control_status !=
                "no_public_negative_control_or_extraction_blank"),
          "PRJNA766558 无公开 blank 边界未完整保留。")
  fail_if(nrow(genus) > 0L &&
            (any(!grepl("low_biomass_FFPE", genus$contamination_risk,
                        fixed = TRUE)) ||
               any(genus$function_claim_scope !=
                     "taxonomic_ecology_only_16S_function_potential_not_computed")),
          "PRJNA766558 FFPE低生物量污染边界或16S未测功能边界漂移。")

  ancombc2 <- x$microbiome_ancombc2
  fail_if(!nrow(ancombc2) ||
            uniqueN(ancombc2$analysis_status) != 1L ||
            !(unique(ancombc2$analysis_status) %chin% c(
              "completed_random_intercept_sensitivity",
              "failed_recorded_not_used_for_gate"
            )),
          "ANCOM-BC2 正交敏感性状态缺失或不可识别。")
  if (unique(ancombc2$analysis_status) ==
      "completed_random_intercept_sensitivity") {
    require_columns(
      ancombc2,
      c(
        "lfc_tissue_roletumor", "q_tissue_roletumor",
        "diff_robust_tissue_roletumor"
      ),
      "completed ANCOM-BC2 sensitivity"
    )
    fail_if(anyNA(ancombc2$taxon) || anyDuplicated(ancombc2$taxon) > 0L,
            "完成的 ANCOM-BC2 taxon 身份缺失或重复。")
  }

  drivers <- x$integrated_drivers
  drivers[, `:=`(
    cross_cohort_patient_link_created = as_logical_strict(
      cross_cohort_patient_link_created, "driver cross-cohort link"
    ),
    exact_independent_validation_present = as_logical_strict(
      exact_independent_validation_present, "driver exact validation"
    ),
    high_priority_validation_rule = as_logical_strict(
      high_priority_validation_rule, "driver high-priority rule"
    )
  )]
  strong <- drivers[decision == "strong_patient_level_candidate"]
  fail_if(nrow(strong) != 12L || uniqueN(strong$gene_name) != 12L,
          "正式 integrated-driver 源表必须完整保留 12 个 strong drivers。")
  fail_if(any(strong$cross_cohort_patient_link_created) ||
            any(strong$exact_independent_validation_present),
          "整合 driver 出现伪造跨队列患者连接或精确独立验证。")

  axes <- x$integrated_axis_edges
  axes[, countable_as_exact_driver_validation := as_logical_strict(
    countable_as_exact_driver_validation, "axis exact validation"
  )]
  axes[, independent_from_tcga_discovery := as_logical_strict(
    independent_from_tcga_discovery, "axis independence from TCGA"
  )]
  fail_if(any(axes$countable_as_exact_driver_validation),
          "当前轴边不应计为精确 driver 事件独立验证。")
  raw_progeny <- axes[
    edge_class == "TCGA_driver_state_internal" &
      target_layer == "RNA_derived_pathway_activity" & evidence_tier == "T2"
  ]
  adjusted_progeny <- axes[
    edge_class == "TCGA_ECMS_adjusted_factor_PROGENy" & evidence_tier == "T2"
  ]
  fail_if(nrow(raw_progeny) != 23L || uniqueN(raw_progeny$target_node) != 9L ||
            uniqueN(raw_progeny, by = c("axis_candidate_id", "target_node")) != 22L ||
            nrow(adjusted_progeny) != 31L ||
            uniqueN(adjusted_progeny$target_node) != 14L ||
            any(raw_progeny$independent_from_tcga_discovery) ||
            any(adjusted_progeny$independent_from_tcga_discovery),
          paste(
            "Figure 5 raw/ECMS-adjusted PROGENy 边的数量、节点或",
            "非独立证据边界漂移。"
          ))

  heterogeneity <- x$heterogeneity_axes
  heterogeneity[, `:=`(
    level_factor_soft_flag = as_logical_strict(
      level_factor_soft_flag, "integrated level-factor"
    ),
    discrete_subtype_support = as_logical_strict(
      discrete_subtype_support, "discrete subtype support"
    )
  )]
  fail_if(nrow(heterogeneity) != 8L || uniqueN(heterogeneity$factor) != 8L,
          "异质性表必须完整保留 Factor1–Factor8。")
  fail_if(nrow(heterogeneity[factor == "Factor4" &
                               level_factor_soft_flag & evidence_tier == "T0"]) != 1L,
          "Factor4 必须在最终整合表中冻结为 T0。")
  fail_if(any(heterogeneity$discrete_subtype_support),
          "连续异质性图不得出现已支持离散亚型。")

  modules <- x$module_summaries
  modules[, `:=`(
    formal_analysis_available = as_logical_strict(
      formal_analysis_available, "module availability"
    ),
    host_sample_level_link_to_other_modules_allowed = as_logical_strict(
      host_sample_level_link_to_other_modules_allowed, "module host link"
    ),
    exact_driver_event_independent_validation = as_logical_strict(
      exact_driver_event_independent_validation, "module exact validation"
    )
  )]
  orthogonal_modules <- modules[module_id %chin% c(
    "PR001876_METABOLOME", "PRJNA766558_MICROBIOME"
  )]
  fail_if(nrow(orthogonal_modules) != 2L ||
            any(!orthogonal_modules$formal_analysis_available) ||
            any(orthogonal_modules$host_sample_level_link_to_other_modules_allowed) ||
            any(orthogonal_modules$exact_driver_event_independent_validation),
          "代谢/微生物模块必须完整、断开且不冒充精确验证。")

  invisible(TRUE)
}

validate_semantics(inputs, verified_families)

font_match <- systemfonts::match_fonts("Arial")
fail_if(is.null(font_match$path) || !nzchar(font_match$path) ||
          !file_exists(font_match$path),
        "未找到 Arial 字体；R-only Nature 图件禁止静默替换字体。")

if (validate_only) {
  cat("VALIDATE_ONLY_OK\n")
  cat("source_families\t", paste(names(verified_families), collapse = ";"),
      "\n", sep = "")
  cat("official_ecms_primary\t78\n")
  cat("factor4\tT0_level_factor\n")
  cat("formal_outputs\t20 figures plus manifest-last\n")
  quit(save = "no", status = 0L)
}

font_family <- "Arial"
base_size <- 6.2
minimum_text_pt <- 5.5
palette_contract <- c(
  neutral_dark = "#2B2B2B",
  neutral_mid = "#7C8288",
  neutral_light = "#D9DDE1",
  neutral_pale = "#F2F3F4",
  signal_blue = "#4F81BD",
  signal_blue_light = "#AFC8E6",
  signal_teal = "#6CAFA7",
  accent_orange = "#D79A55",
  accent_red = "#C76B67",
  boundary = "#B39A6A"
)
tier_palette <- c(
  T0 = "#B8BDC2",
  T1 = "#E4E7E9",
  T2 = palette_contract[["signal_blue_light"]],
  T3 = palette_contract[["signal_blue"]],
  # T4 与 T3 同色；T4 只靠双外框表示湿实验优先级。
  T4 = palette_contract[["signal_blue"]]
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
      legend.key.height = grid::unit(3.2, "mm"),
      plot.margin = margin(2.2, 2.2, 2.2, 2.2, unit = "mm")
    )
}
theme_set(theme_nature_contract())

short_label <- function(value, width = 34L) {
  value <- as.character(value)
  needs_trim <- nchar(value) > width
  value[needs_trim] <- paste0(substr(value[needs_trim], 1L, width - 1L), "…")
  value
}

wrap_label <- function(value, width = 38L) {
  vapply(as.character(value), function(item) {
    paste(strwrap(item, width = width), collapse = "\n")
  }, character(1))
}

format_q <- function(value) {
  output <- rep("NA", length(value))
  finite <- is.finite(value)
  output[finite & value >= 0.001] <- sprintf("%.3f", value[finite & value >= 0.001])
  output[finite & value < 0.001] <- format(
    value[finite & value < 0.001], scientific = TRUE, digits = 2
  )
  output
}

empty_panel <- function(title, message, subtitle = NULL) {
  ggplot(data.table(x = 0.5, y = 0.5), aes(x, y)) +
    annotate(
      "rect", xmin = 0.05, xmax = 0.95, ymin = 0.15, ymax = 0.85,
      fill = palette_contract[["neutral_pale"]], colour = "#C8CDD1",
      linewidth = 0.35
    ) +
    annotate(
      "text", x = 0.5, y = 0.5, label = message,
      family = font_family, size = 2.0, colour = palette_contract[["neutral_dark"]],
      lineheight = 0.95
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_family = font_family, base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 0.4, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.4, colour = "#555555"),
      plot.margin = margin(2.2, 2.2, 2.2, 2.2, unit = "mm")
    )
}

panel_title_theme <- theme(
  plot.title.position = "plot",
  legend.position = "right"
)

build_figure1 <- function(x) {
  plan <- copy(x$leakage_model_plan)
  plan[, model_status := fcase(
    scenario_no_op & model_reused, "no-op / baseline reuse",
    model_completed & model_converged & grepl("checkpoint", run_status),
      "checkpoint / converged",
    model_completed & model_converged, "completed / converged",
    model_completed & !model_converged, "completed / not converged",
    default = "not completed"
  )]
  plan[, scenario_label := fcase(
    scenario_id == "baseline", "baseline",
    scenario_id == "drop_mutation_view", "− mutation view",
    scenario_id == "drop_cnv_view", "− CNV view",
    scenario_id == "drop_both_event_views", "− both event views",
    grepl("^drop_gene_event_features_", scenario_id),
      paste0("− ", sub("^drop_gene_event_features_", "", scenario_id), " features"),
    default = short_label(gsub("_", " ", scenario_id), 22L)
  )]
  plan[, scenario_label := factor(
    scenario_label,
    levels = rev(unique(scenario_label[order(scenario_order)]))
  )]
  plan[, seed_label := factor(
    paste0("s", sub("^2026", "", seed)),
    levels = paste0("s", sub("^2026", "", sort(unique(seed))))
  )]
  plan_palette <- c(
    "completed / converged" = palette_contract[["signal_blue"]],
    "checkpoint / converged" = palette_contract[["signal_teal"]],
    "no-op / baseline reuse" = palette_contract[["neutral_mid"]],
    "completed / not converged" = palette_contract[["accent_orange"]],
    "not completed" = palette_contract[["accent_red"]]
  )
  p_a <- ggplot(plan, aes(x = seed_label, y = scenario_label, fill = model_status)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(
      aes(label = ifelse(scenario_no_op, "N/A", "✓")),
      colour = "white", family = font_family, size = 1.95, fontface = "bold"
    ) +
    scale_fill_manual(values = plan_palette, drop = FALSE) +
    labs(
      title = "Model plan and execution",
      subtitle = "3 prespecified seeds × 8 scenarios; N/A means baseline reuse",
      x = NULL, y = NULL, fill = "Execution state"
    ) +
    theme_nature_contract() +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.y = element_text(size = 5.5), legend.position = "bottom",
      plot.margin = margin(3, 3, 3, 5, unit = "mm")
    ) +
    panel_title_theme

  alignment <- copy(x$leakage_factor_alignment)[!is.na(reference_factor)]
  alignment[, scenario_label := fcase(
    scenario_id == "baseline", "base",
    scenario_id == "drop_mutation_view", "−mut",
    scenario_id == "drop_cnv_view", "−CNV",
    scenario_id == "drop_both_event_views", "−both",
    grepl("^drop_gene_event_features_", scenario_id),
      paste0("−", sub("^drop_gene_event_features_", "", scenario_id)),
    default = short_label(gsub("_", " ", scenario_id), 12L)
  )]
  alignment[, model_label := factor(
    paste0("s", seed, " · ", scenario_label),
    levels = rev(unique(paste0("s", seed, " · ", scenario_label)[
      order(seed, scenario_order)
    ]))
  )]
  alignment[, reference_factor := factor(
    reference_factor, levels = paste0("Factor", 1:8)
  )]
  p_b <- ggplot(
    alignment,
    aes(
      x = reference_factor, y = model_label,
      fill = absolute_shared_loading_correlation,
      colour = match_reliable
    )
  ) +
    geom_tile(linewidth = 0.34) +
    geom_point(
      data = alignment[low_margin_flag == TRUE],
      shape = 16, size = 0.75, colour = palette_contract[["accent_red"]]
    ) +
    scale_fill_gradient(
      low = "#F2F3F4", high = palette_contract[["signal_blue"]],
      limits = c(0, 1), oob = scales::squish, na.value = "#D0D3D5"
    ) +
    scale_colour_manual(
      values = c(`TRUE` = "#202020", `FALSE` = "white"),
      labels = c(`TRUE` = "reliable outline", `FALSE` = "unreliable"),
      drop = FALSE
    ) +
    labs(
      title = "Event-blind factor alignment",
      subtitle = "Fill: |loading r|; outline: reliable; red dot: low margin",
      x = NULL, y = NULL, fill = "|loading r|", colour = "Match"
    ) +
    theme_nature_contract() +
    theme(
      axis.text.y = element_text(size = 5.5, lineheight = 0.9),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.line = element_blank(), axis.ticks = element_blank(),
      legend.position = "right", plot.margin = margin(3, 4, 3, 3, unit = "mm")
    ) +
    panel_title_theme

  stability <- copy(x$leakage_edge_stability)
  candidate <- copy(x$leakage_candidate_summary)
  candidate[, leakage_gate_pass := as_logical_strict(
    leakage_gate_pass, "Figure 1 leakage gate"
  )]
  event_short <- function(value) fcase(
    value == "relative_cnv", "relCNV",
    value == "amplification", "amp",
    value == "mutation", "mut",
    value == "deletion", "del",
    default = short_label(value, 8L)
  )
  stability[, edge_label := paste0(
    gene_name, " · ", event_short(event_type), "→",
    sub("^Factor", "F", reference_factor)
  )]
  candidate[, edge_label := paste0(
    gene_name, " · ", event_short(event_type), "→",
    sub("^Factor", "F", reference_factor)
  )]
  edge_order <- candidate[order(
    -as.integer(leakage_gate_pass), narrative_role, gene_name,
    event_type, reference_factor
  ), edge_label]
  fail_if(length(edge_order) != 41L || anyDuplicated(edge_order) > 0L,
          "Figure 1 全宽审计矩阵必须唯一覆盖 41 条 original edges。")

  stability[, audit_column := fcase(
    scenario_id == "baseline", "Base",
    scenario_id == "drop_mutation_view", "−Mut",
    scenario_id == "drop_cnv_view", "−CNV",
    scenario_id == "drop_both_event_views", "−Both",
    scenario_class == "gene_event_feature_drop", "Gene−",
    default = NA_character_
  )]
  stability[, cell_state := fcase(
    scenario_no_op, "N/A",
    grepl("^retained", retention_decision), "retained",
    grepl("direction_and_magnitude|direction_retained", retention_decision),
      "attenuated",
    grepl("unstable|incomplete", retention_decision), "unstable",
    default = "not retained"
  )]
  scenario_audit <- stability[!is.na(audit_column), .(
    cell_state = cell_state[[1L]], cell_text = ""
  ), by = .(edge_label, audit_column)]
  scenario_grid <- CJ(
    edge_label = edge_order,
    audit_column = c("Base", "−Mut", "−CNV", "−Both", "Gene−"),
    unique = TRUE
  )
  scenario_audit <- merge(
    scenario_grid, scenario_audit,
    by = c("edge_label", "audit_column"), all.x = TRUE, sort = FALSE
  )
  scenario_audit[is.na(cell_state), `:=`(cell_state = "N/A", cell_text = "")]

  rate_columns <- c(
    "drop_both_direction_retention_rate",
    "drop_both_magnitude_retention_rate",
    "drop_both_gate_retention_rate"
  )
  candidate[, (rate_columns) := lapply(.SD, as.numeric), .SDcols = rate_columns]
  rates <- melt(
    candidate,
    id.vars = "edge_label",
    measure.vars = rate_columns,
    variable.name = "metric", value.name = "rate"
  )
  rates[, audit_column := factor(
    metric, levels = rate_columns, labels = c("Dir", "Mag", "Support retention")
  )]
  rates[, `:=`(
    cell_state = fcase(
      !is.finite(rate), "N/A",
      rate >= 0.999, "retained",
      rate > 0, "attenuated",
      default = "not retained"
    ),
    cell_text = fifelse(is.finite(rate), sub("^0", "", sprintf("%.2f", rate)), "N/A")
  )]
  gate_audit <- candidate[, .(
    edge_label,
    audit_column = "Robust",
    cell_state = fifelse(leakage_gate_pass, "retained", "not retained"),
    cell_text = fifelse(leakage_gate_pass, "✓", "×")
  )]
  reason_audit <- candidate[, .(
    edge_label,
    audit_column = "Fail",
    cell_state = "annotation",
    cell_text = fcase(
      leakage_gate_pass, "pass",
      grepl("original_edge_not", gate_failure_reason), "orig",
      grepl("level_factor|Factor4", gate_failure_reason), "lvl",
      grepl("not_all|not_converged|no_op", gate_failure_reason), "run",
      grepl("match_rate", gate_failure_reason), "aln",
      grepl("direction_retention", gate_failure_reason), "dir",
      grepl("magnitude_retention", gate_failure_reason), "mag",
      grepl("supported_gate_retention", gate_failure_reason), "sup",
      default = "other"
    )
  )]
  audit_matrix <- rbindlist(list(
    scenario_audit,
    rates[, .(edge_label, audit_column = as.character(audit_column),
              cell_state, cell_text)],
    gate_audit,
    reason_audit
  ), use.names = TRUE)
  audit_columns <- c(
    "Base", "−Mut", "−CNV", "−Both", "Gene−",
    "Dir", "Mag", "Support retention", "Robust", "Fail"
  )
  audit_matrix[, `:=`(
    edge_label = factor(edge_label, levels = rev(edge_order)),
    audit_column = factor(audit_column, levels = audit_columns)
  )]
  fail_if(nrow(audit_matrix) != 41L * length(audit_columns),
          "Figure 1 全宽审计矩阵必须是 41×10 完整网格。")
  audit_palette <- c(
    retained = palette_contract[["signal_blue"]],
    attenuated = palette_contract[["accent_orange"]],
    unstable = palette_contract[["boundary"]],
    `not retained` = palette_contract[["neutral_light"]],
    `N/A` = "#F2F3F4", annotation = "#E8EAEC"
  )
  p_c <- ggplot(
    audit_matrix,
    aes(x = audit_column, y = edge_label, fill = cell_state)
  ) +
    geom_tile(colour = "white", linewidth = 0.16) +
    geom_text(
      aes(label = cell_text), family = font_family, size = 1.75,
      colour = "#303438", fontface = "plain"
    ) +
    scale_fill_manual(values = audit_palette, drop = FALSE) +
    labs(
      title = "All 41 original edges: retention rates and prespecified support criteria",
      subtitle = "Full-width audit; exact values and failure strings remain in the formal source tables",
      x = NULL, y = NULL, fill = "Audit state",
      caption = "Fail codes: orig original edge; lvl level-factor; run incomplete/no-op; aln alignment; dir direction; mag magnitude; sup support-retention criterion."
    ) +
    theme_nature_contract() +
    theme(
      axis.text.y = element_text(size = 5.5, lineheight = 0.88),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 5.5),
      axis.line = element_blank(), axis.ticks = element_blank(),
      legend.position = "bottom", legend.box = "horizontal",
      plot.caption = element_text(size = 5.5, hjust = 0),
      plot.margin = margin(3, 4, 3, 5, unit = "mm")
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))

  main_genes <- c("PIK3CA", "NFE2L2", "GNAS", "ZNF750")
  association <- copy(x$leakage_associations)[
    gene_name %chin% main_genes & !is.na(original_edge_id) &
      scenario_id %chin% c("baseline", "drop_both_event_views")
  ]
  if (!nrow(association)) {
    p_d <- empty_panel(
      "Narrative candidates across seeds",
      "N/A — none of the four prespecified genes occurs in the original edge set"
    )
  } else {
    association[, edge_label := paste0(
      gene_name, " · ", event_short(event_type), "→",
      sub("^Factor", "F", reference_factor)
    )]
    association[, scenario_label := factor(
      scenario_id,
      levels = c("baseline", "drop_both_event_views"),
      labels = c("baseline", "drop both event views")
    )]
    plotting <- association[is.finite(aligned_effect)]
    effect_order <- unique(plotting[order(
      match(gene_name, main_genes), event_type, reference_factor
    ), edge_label])
    plotting[, edge_label := factor(edge_label, levels = rev(effect_order))]
    p_d <- ggplot(
      plotting,
      aes(
        x = aligned_effect, y = edge_label,
        colour = scenario_label, shape = factor(seed)
      )
    ) +
      geom_vline(xintercept = 0, linewidth = 0.3, colour = "#888888") +
      geom_point(
        size = 1.25, alpha = 0.9,
        position = position_jitter(width = 0, height = 0.10, seed = 20260712)
      ) +
      scale_colour_manual(values = c(
        baseline = palette_contract[["neutral_dark"]],
        `drop both event views` = palette_contract[["signal_blue"]]
      )) +
      labs(
        title = "Prespecified genes: baseline versus event-view deletion",
        subtitle = "Per-seed effects remain TCGA-internal; all exact rows remain in the formal table",
        x = "Aligned association effect", y = NULL,
        colour = "Scenario", shape = "Seed"
      ) +
      theme_nature_contract() +
      theme(
        axis.text.y = element_text(size = 5.5), legend.position = "bottom",
        plot.margin = margin(3, 4, 3, 5, unit = "mm")
      ) +
      guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
  }

  top_row <- (p_a | p_b) + plot_layout(widths = c(0.92, 1.08))
  top_row /
    p_c /
    p_d +
    plot_layout(heights = c(0.92, 1.90, 0.72)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 8, face = "bold", family = font_family))
}

build_figure2 <- function(x) {
  probabilities <- copy(x$ecms_probabilities)
  primary <- probabilities[in_official_78 & eligible_for_primary_association]
  cohort_counts <- data.table(
    scope = factor(
      c("Official78", "Extension16"),
      levels = c("Official78", "Extension16")
    ),
    n = c(nrow(primary), nrow(probabilities[in_official_78 == FALSE])),
    role = c("primary", "conditional boundary")
  )
  p_a_counts <- ggplot(cohort_counts, aes(x = scope, y = n, fill = role)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = n), vjust = -0.25, family = font_family, size = 2.1) +
    scale_fill_manual(values = c(
      primary = palette_contract[["signal_blue"]],
      `conditional boundary` = palette_contract[["neutral_light"]]
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = "Anchor scope",
      subtitle = "Primary: Official78\nBoundary: Extension16",
      x = NULL, y = "Patients", fill = NULL
    ) +
    theme_nature_contract() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 5.5, angle = 35, hjust = 1)
    )

  calibration <- copy(x$ecms_calibration)[
    gate_role %chin% c("prelocked_extension_gate", "extension_decision")
  ]
  calibration[, gate_state := fifelse(pass, "pass", "not pass")]
  calibration[, metric_label := fcase(
    grepl("overall_extension", metric), "overall",
    grepl("median_class_probability", metric), "median probability",
    grepl("exact_label_agreement", metric), "label agreement",
    grepl("cohen_kappa", metric), "Cohen κ",
    grepl("adjusted_rand", metric), "ARI",
    default = short_label(gsub("_", " ", metric), 18L)
  )]
  p_a_gates <- ggplot(
    calibration, aes(x = "extension calibration", y = metric_label, fill = gate_state)
  ) +
    geom_tile(colour = "white", linewidth = 0.25) +
    geom_text(
      aes(label = paste0(
        ifelse(gate_state == "pass", "✓ ", "× "),
        ifelse(is.finite(observed), sprintf("%.3f", observed), "")
      )),
      family = font_family, size = 1.95
    ) +
    scale_fill_manual(values = c(
      pass = palette_contract[["signal_blue_light"]],
      `not pass` = palette_contract[["neutral_light"]]
    )) +
    labs(
      title = "Extension calibration",
      subtitle = "2/4 criteria met; not used in primary analysis",
      x = NULL, y = NULL
    ) +
    theme_nature_contract() +
    theme(
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      axis.line = element_blank(), legend.position = "none"
    )
  p_a <- (p_a_counts | p_a_gates) + plot_layout(widths = c(0.78, 1.22))

  scores <- copy(x$factor_scores)
  scores <- merge(
    scores,
    primary[, .(patient_id, official_anchor_label)],
    by = "patient_id", all = FALSE, sort = FALSE
  )
  fail_if(nrow(scores) != 78L,
          "Figure 2 的患者级因子分布未严格限制为官方 78 例。")
  score_long <- melt(
    scores,
    id.vars = c("patient_id", "official_anchor_label"),
    measure.vars = paste0("Factor", 1:8),
    variable.name = "factor", value.name = "score"
  )
  factor_means <- score_long[, .(
    mean_score = mean(score),
    n = .N
  ), by = .(factor, official_anchor_label)]
  factor_means[, factor := factor(factor, levels = rev(paste0("Factor", 1:8)))]
  factor_means[, official_anchor_label := factor(
    official_anchor_label, levels = paste0("ECMS", 1:4)
  )]
  factor_means[, display_score := fifelse(
    as.character(factor) == "Factor4", NA_real_, mean_score
  )]
  p_b_heat <- ggplot(
    factor_means,
    aes(x = official_anchor_label, y = factor, fill = display_score)
  ) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(
      data = factor_means[
        as.character(factor) == "Factor4" & official_anchor_label == "ECMS2"
      ],
      aes(label = "technical"), family = font_family, size = 1.75,
      colour = "#555555"
    ) +
    scale_fill_gradient2(
      low = "#C57A78", mid = "white", high = palette_contract[["signal_blue"]],
      midpoint = 0, na.value = tier_palette[["T0"]]
    ) +
    labs(
      title = "Official-anchor factor means",
      subtitle = "Factor4 is greyed as a technical/background level-factor",
      x = NULL, y = NULL, fill = "Mean score"
    ) +
    theme_nature_contract() +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 5.5)
    )

  factor_primary <- copy(x$ecms_factor)[is_primary_scope == TRUE]
  factor_primary[, q_min := pmin(anova_q_value, kruskal_q_value, na.rm = TRUE)]
  factor_primary[!is.finite(q_min), q_min := NA_real_]
  factor_primary[, q_short := fifelse(
    !is.finite(q_min), "NA",
    fifelse(
      q_min < 0.001,
      formatC(q_min, format = "e", digits = 1),
      sprintf("%.2f", q_min)
    )
  )]
  factor_primary[, factor := factor(factor, levels = rev(paste0("Factor", 1:8)))]
  p_b_stats <- ggplot(factor_primary, aes(x = 1, y = factor)) +
    geom_text(
      aes(label = paste0("η² ", sprintf("%.2f", eta_squared),
                         "\nq ", q_short)),
      family = font_family, size = 1.95, hjust = 0.5
    ) +
    scale_x_continuous(limits = c(0.5, 1.5)) +
    labs(title = "η² / q") +
    theme_void(base_family = font_family, base_size = base_size) +
    theme(plot.title = element_text(size = base_size + 0.4, face = "bold"))
  p_b <- (p_b_heat | p_b_stats) + plot_layout(widths = c(1, 0.48))

  heterogeneity <- copy(x$heterogeneity_axes)
  priority <- heterogeneity[
    !level_factor_soft_flag & ecms_context_tier == "T2"
  ][order(-ecms_eta_squared, ecms_context_q), head(.SD, 3L)]
  if (!nrow(priority)) {
    p_c <- empty_panel(
      "Priority factor distributions",
      "True empty state: no non-level factor met the interpretive ECMS rule"
    )
  } else {
    selected <- score_long[factor %chin% priority$factor]
    selected[, official_anchor_label := factor(
      official_anchor_label, levels = paste0("ECMS", 1:4)
    )]
    p_c <- ggplot(
      selected,
      aes(x = official_anchor_label, y = score, fill = official_anchor_label)
    ) +
      geom_violin(scale = "width", trim = FALSE, linewidth = 0.25, alpha = 0.55) +
      geom_boxplot(
        width = 0.20, outlier.shape = NA, linewidth = 0.3,
        fill = "white", alpha = 0.75
      ) +
      geom_point(
        position = position_jitter(width = 0.08, height = 0, seed = 20260712),
        size = 0.45, alpha = 0.45, colour = "#303030"
      ) +
      facet_wrap(~ factor, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = c(
        ECMS1 = "#AFC8E6", ECMS2 = "#BBD9D4",
        ECMS3 = "#D5C9E5", ECMS4 = "#E3C6A8"
      )) +
      labs(
        title = "Official78 distributions for ECMS-associated factors",
        subtitle = "Shared RNA representation; explanatory, not validation",
        x = NULL, y = "MOFA factor score"
      ) +
      theme_nature_contract() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))
  }

  adjusted <- copy(x$ecms_adjusted)[is_primary_scope == TRUE]
  formal_adjusted <- copy(x$integrated_axis_edges)[
    edge_class == "TCGA_ECMS_adjusted_factor_PROGENy"
  ]
  adjusted[, t2_key := paste(factor, pathway, sep = "::")]
  t2_keys <- formal_adjusted[evidence_tier == "T2", paste(
    source_node, target_node, sep = "::"
  )]
  adjusted[, formal_t2 := t2_key %chin% t2_keys]
  adjusted[, factor4 := factor == "Factor4"]
  adjusted[, display_beta := fifelse(factor4, NA_real_, factor_beta_standardized)]
  adjusted[, factor := factor(factor, levels = rev(paste0("Factor", 1:8)))]
  pathway_order <- adjusted[, .(
    best_q = min(incremental_q_value, na.rm = TRUE),
    max_r2 = max(partial_r_squared, na.rm = TRUE)
  ), by = pathway][order(best_q, -max_r2), pathway]
  adjusted[, pathway := factor(pathway, levels = pathway_order)]
  p_d <- ggplot(adjusted, aes(x = pathway, y = factor)) +
    geom_tile(aes(fill = display_beta), colour = "white", linewidth = 0.20) +
    geom_point(
      aes(size = partial_r_squared, shape = formal_t2),
      fill = "white", colour = "#202020", stroke = 0.35
    ) +
    geom_text(
      data = adjusted[factor4 == TRUE][1L], aes(label = "technical"),
      family = font_family, size = 1.75, colour = "#555555"
    ) +
    scale_fill_gradient2(
      low = "#C57A78", mid = "white", high = palette_contract[["signal_blue"]],
      midpoint = 0, na.value = tier_palette[["T0"]]
    ) +
    scale_shape_manual(
      values = c(`TRUE` = 21, `FALSE` = 1),
      labels = c(`TRUE` = "conditional", `FALSE` = "below rule")
    ) +
    scale_size_continuous(range = c(0.4, 2.2)) +
    labs(
      title = "ECMS-adjusted factor–pathway map",
      subtitle = "Fill β; size partial R²; shared RNA; interpretive only",
      x = NULL, y = NULL, fill = "β", size = "Partial R²",
      shape = "Status"
    ) +
    theme_nature_contract() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
      axis.line = element_blank(), axis.ticks = element_blank(),
      legend.position = "bottom"
    ) +
    guides(
      fill = guide_colourbar(
        title.position = "top", barwidth = grid::unit(22, "mm")
      ),
      size = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)
    )

  cluster <- copy(x$cluster_evaluation)
  cluster_long <- melt(
    cluster,
    id.vars = c("k", "conditional_stability", "selected_for_description"),
    measure.vars = c(
      "snf_mean_silhouette", "mofa_consensus_pac", "snf_mofa_adjusted_rand"
    ),
    variable.name = "metric", value.name = "value"
  )
  cluster_long[, metric := factor(
    metric,
    levels = c(
      "snf_mean_silhouette", "mofa_consensus_pac", "snf_mofa_adjusted_rand"
    ),
    labels = c("SNF silhouette", "MOFA PAC", "SNF–MOFA ARI")
  )]
  p_e <- ggplot(cluster_long, aes(x = k, y = value)) +
    geom_hline(yintercept = 0, colour = "#D0D3D5", linewidth = 0.25) +
    geom_line(colour = palette_contract[["neutral_mid"]], linewidth = 0.35) +
    geom_point(
      aes(fill = selected_for_description), shape = 21, size = 1.4,
      colour = "#303030", stroke = 0.35
    ) +
    facet_wrap(~ metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c(
      `TRUE` = palette_contract[["accent_orange"]], `FALSE` = "white"
    )) +
    scale_x_continuous(breaks = sort(unique(cluster$k))) +
    labs(
      title = "Discrete-cluster stability assessment",
      subtitle = "No k=2–6 solution met the prespecified stability criteria",
      x = "k", y = "Metric value"
    ) +
    theme_nature_contract() +
    theme(legend.position = "none")

  top_row <- (wrap_elements(full = p_a) | p_c) +
    plot_layout(widths = c(1.00, 1.48))
  middle_row <- (wrap_elements(full = p_b) | p_d) +
    plot_layout(widths = c(0.92, 1.56))
  top_row / middle_row / p_e +
    plot_layout(heights = c(0.82, 1.22, 0.66)) +
    plot_annotation(tag_levels = "a") &
    theme(
      plot.tag = element_text(size = 8, face = "bold", family = font_family),
      plot.margin = margin(3, 4, 3, 4, unit = "mm")
    )
}

build_figure3 <- function(x) {
  summary <- copy(x$cao_summary)
  driver_order <- x$integrated_drivers[
    decision == "strong_patient_level_candidate",
    gene_name
  ]
  summary[, gene := factor(tcga_gene_name, levels = rev(driver_order))]

  effect_panel <- function(value_column, title, x_label, colour,
                           availability_column = NULL) {
    plotting <- summary[, .(
      gene,
      value = as.numeric(get(value_column)),
      available = if (is.null(availability_column)) {
        is.finite(as.numeric(get(value_column)))
      } else {
        !grepl("not_identified|not_quantified|unavailable|missing",
               as.character(get(availability_column)), ignore.case = TRUE) &
          is.finite(as.numeric(get(value_column)))
      }
    )]
    ggplot(plotting[available == TRUE], aes(x = value, y = gene)) +
      geom_vline(xintercept = 0, linewidth = 0.30, colour = "#92979B") +
      geom_point(shape = 21, size = 1.7, fill = colour, colour = "#303030") +
      geom_text(
        data = plotting[available == FALSE],
        aes(x = 0, y = gene, label = "N/A"),
        inherit.aes = FALSE, family = font_family, size = 1.95,
        colour = palette_contract[["neutral_mid"]]
      ) +
      labs(title = title, x = x_label, y = NULL) +
      theme_nature_contract() +
      theme(axis.text.y = element_text(size = 5.5))
  }

  p_a_promoter <- effect_panel(
    "promoter_median_paired_delta_beta", "Promoter methylation",
    "Median paired Δβ", "#B8A5CD", "promoter_region_evaluation_status"
  )
  p_a_rna <- effect_panel(
    "limma_log_fc_pc1", "RNA",
    "limma logFC (T−N; pc=1)", palette_contract[["signal_blue"]]
  ) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_a_protein <- effect_panel(
    "median_patient_log2_tumor_vs_normal", "Protein",
    "Median log2(T/N)", palette_contract[["signal_teal"]],
    "protein_measurement_status"
  ) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p_a <- p_a_promoter + p_a_rna + p_a_protein +
    plot_layout(widths = c(1.15, 1, 1)) +
    plot_annotation(
      title = "Three independent effect scales",
      subtitle = "Promoter Δβ, RNA logFC and protein log2(T/N) never share an axis"
    ) &
    theme(
      plot.title = element_text(size = base_size + 0.5, face = "bold",
                                family = font_family),
      plot.subtitle = element_text(size = base_size - 0.4, colour = "#555555",
                                   family = font_family)
    )

  availability <- rbindlist(list(
    summary[, .(
      gene,
      layer = "promoter",
      state = fcase(
        promoter_region_evaluation_status == "evaluable", "evaluable",
        grepl("technical", promoter_region_evaluation_status),
          "technical N/A",
        default = "unavailable"
      )
    )],
    summary[, .(
      gene,
      layer = "RNA",
      state = fifelse(is.finite(limma_log_fc_pc1), "evaluable", "unavailable")
    )],
    summary[, .(
      gene,
      layer = "protein",
      state = fcase(
        protein_measurement_status == "quantified", "quantified",
        protein_measurement_status == "identified_not_quantified",
          "identified—not quantified",
        protein_measurement_status == "not_identified", "not identified",
        default = "unavailable"
      )
    )],
    summary[, .(
      gene,
      layer = "cross-layer",
      state = fcase(
        cross_layer_class == "three_layer_directional_hypothesis", "directional",
        cross_layer_class == "credible_reverse_retained", "credible reverse",
        default = "context/weak"
      )
    )]
  ))
  availability[, layer := factor(
    layer, levels = c("promoter", "RNA", "protein", "cross-layer")
  )]
  availability[, state_code := fcase(
    state == "evaluable", "E", state == "quantified", "Q",
    state == "technical N/A", "TECH", state == "unavailable", "NA",
    state == "identified—not quantified", "ID",
    state == "not identified", "NI", state == "directional", "DIR",
    state == "credible reverse", "REV", default = "CTX"
  )]
  availability_palette <- c(
    evaluable = "#DCE8F4", quantified = "#BBD9D4",
    `technical N/A` = "#ECEDEE", unavailable = "#E4E6E8",
    `identified—not quantified` = "#E4D8C9", `not identified` = "#D7D9DB",
    directional = "#AFC8E6", `credible reverse` = "#D7A6A3",
    `context/weak` = "#E2E4E6"
  )
  p_b <- ggplot(availability, aes(x = layer, y = gene, fill = state)) +
    geom_tile(colour = "white", linewidth = 0.28) +
    geom_text(aes(label = state_code), family = font_family, size = 1.75) +
    scale_fill_manual(values = availability_palette, drop = FALSE) +
    labs(
      title = "Layer availability",
      subtitle = "E/Q measured; ID/NI missing",
      x = NULL, y = NULL, fill = "State"
    ) +
    theme_nature_contract() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.line = element_blank(), axis.ticks = element_blank(),
      legend.position = "none"
    )

  patient_effect_panel <- function(gene_name, title) {
    patient <- copy(x$cao_patient)[tcga_gene_name == gene_name]
    if (!nrow(patient)) {
      return(empty_panel(title, paste("True empty state:", gene_name,
                                      "has no Cao patient-level rows")))
    }
    patient_long <- melt(
      patient,
      id.vars = c(
        "tcga_gene_name", "paper_patient_id", "coverage_gate_promoter",
        "protein_patient_measurement_status"
      ),
      measure.vars = c(
        "median_delta_beta_promoter", "rna_delta_log2_tpm_pc1",
        "protein_log2_tumor_vs_normal"
      ),
      variable.name = "layer", value.name = "effect"
    )
    patient_long[, layer := factor(
      layer,
      levels = c(
        "median_delta_beta_promoter", "rna_delta_log2_tpm_pc1",
        "protein_log2_tumor_vs_normal"
      ),
      labels = c("Promoter Δβ", "RNA Δlog2 TPM", "Protein log2 T/N")
    )]
    patient_long[, patient_order := as.integer(gsub("[^0-9]", "", paper_patient_id))]
    setorder(patient_long, patient_order)
    patient_long[, paper_patient_id := factor(
      paper_patient_id, levels = unique(paper_patient_id)
    )]
    p <- ggplot(
      patient_long[is.finite(effect)],
      aes(x = paper_patient_id, y = effect)
    ) +
      geom_hline(yintercept = 0, linewidth = 0.28, colour = "#979CA0") +
      geom_point(
        aes(fill = layer), shape = 21, size = 1.25,
        colour = "#303030", stroke = 0.3
      ) +
      facet_wrap(~ layer, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = c(
        `Promoter Δβ` = "#B8A5CD",
        `RNA Δlog2 TPM` = palette_contract[["signal_blue"]],
        `Protein log2 T/N` = palette_contract[["signal_teal"]]
      )) +
      labs(
        title = title,
        subtitle = paste(
          "Independent patient points; separate scales; no connecting lines",
          "or causal arrows"
        ),
        x = "Cao patient", y = "Patient-level effect"
      ) +
      theme_nature_contract() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5),
        legend.position = "none"
      )
    if (gene_name == "ZNF750" &&
        all(!is.finite(patient$protein_log2_tumor_vs_normal))) {
      p <- p + labs(caption = "Protein: not identified — not negative")
    }
    p
  }
  p_c <- patient_effect_panel("GNAS", "GNAS conditional three-layer pattern")
  p_d <- patient_effect_panel("ZNF750", "ZNF750 reverse-direction boundary")

  fdr_columns <- c(
    "promoter_q_candidate_by_region", "limma_q_candidate_pc1", "q_candidate"
  )
  fdr_matrix <- as.matrix(summary[, ..fdr_columns])
  fdr_supported <- rowSums(is.finite(fdr_matrix) & fdr_matrix <= 0.05)
  fail_if(any(fdr_supported != summary$fdr_supported_layer_count),
          "Cao 图中代码化 FDR 层计数与正式 summary 不一致。")
  class_rows <- rbindlist(list(
    summary[, .(
      gene,
      class_axis = "cross-layer class",
      class_value = cross_layer_class
    )],
    summary[, .(
      gene,
      class_axis = "statistical support",
      class_value = statistical_support_class
    )]
  ))
  class_rows[, class_axis := factor(
    class_axis, levels = c("cross-layer class", "statistical support")
  )]
  class_rows[, simplified := fcase(
    grepl("three_layer", class_value), "three-layer directional",
    grepl("credible_reverse", class_value), "credible reverse",
    grepl("no_candidate_layer_fdr", class_value), "0 layers at q≤0.05",
    grepl("context|weak", class_value), "context/weak",
    default = short_label(class_value, 24L)
  )]
  class_rows[, class_code := fcase(
    simplified == "three-layer directional", "DIR",
    simplified == "credible reverse", "REV",
    simplified == "0 layers at q≤0.05", "0q",
    default = "CTX"
  )]
  p_e <- ggplot(class_rows, aes(x = class_axis, y = gene, fill = simplified)) +
    geom_tile(colour = "white", linewidth = 0.28) +
    geom_text(aes(label = class_code), family = font_family, size = 1.75) +
    scale_fill_manual(values = c(
      `three-layer directional` = palette_contract[["signal_blue_light"]],
      `credible reverse` = "#D7A6A3",
      `0 layers at q≤0.05` = palette_contract[["neutral_light"]],
      `context/weak` = "#E7E8E9"
    )) +
    labs(
      title = "Cross-layer classes",
      subtitle = paste0(
        sum(fdr_supported > 0L), "/", nrow(summary),
        " candidate-layer q≤0.05"
      ),
      x = NULL, y = NULL, fill = "Class"
    ) +
    theme_nature_contract() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.line = element_blank(), axis.ticks = element_blank(),
      legend.position = "none"
    )

  left_column <- p_a / (p_b | p_e) + plot_layout(heights = c(1.1, 1))
  right_column <- p_c / p_d + plot_layout(heights = c(1, 1))
  figure <- left_column | right_column
  figure +
    plot_layout(widths = c(1.25, 1.20)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 8, face = "bold", family = font_family))
}

build_figure4 <- function(x) {
  metabolites <- copy(x$metabolite_candidates)
  analysis_inventory <- copy(x$metabolite_analysis_inventory)
  analysis_inventory[, analysis_order := match(
    analysis_id, c("AN004960", "AN004962", "AN004963")
  )]
  setorder(analysis_inventory, analysis_order)
  analysis_inventory[, analysis_label := fcase(
    analysis_id == "AN004960", "AN004960\nGC–MS",
    analysis_id == "AN004962", "AN004962\nLC–MS (+)",
    analysis_id == "AN004963", "AN004963\nLC–MS (−)",
    default = analysis_id
  )]
  analysis_label_lookup <- setNames(
    analysis_inventory$analysis_label, analysis_inventory$analysis_id
  )
  analysis_inventory[, analysis_order := NULL]
  logical_metabolite_columns <- c(
    "direction_concordant", "scale_sensitivity_direction_stable",
    "run_order_perfect_group_block",
    "run_order_sensitivity_direction_concordant",
    "paired_sensitivity_estimable", "paired_sensitivity_direction_concordant"
  )
  for (column in logical_metabolite_columns) {
    set(metabolites, j = column, value = as_logical_strict(
      metabolites[[column]], paste("metabolite", column), allow_na = TRUE
    ))
  }
  tier_order <- c(
    "robust_rank_and_scale_fdr", "conditional_fdr",
    "nominal_concordant_conditional", "run_order_confounded_conditional"
  )
  metabolites[, final_candidate_tier := factor(
    final_candidate_tier, levels = tier_order
  )]
  metabolite_counts <- metabolites[, .N, by = .(
    analysis_id, final_candidate_tier
  )]
  metabolite_counts <- merge(
    CJ(
      analysis_id = analysis_inventory$analysis_id,
      final_candidate_tier = tier_order,
      unique = TRUE
    ),
    metabolite_counts,
    by = c("analysis_id", "final_candidate_tier"),
    all.x = TRUE,
    sort = FALSE
  )
  metabolite_counts[is.na(N), N := 0L]
  metabolite_counts[, final_candidate_tier := factor(
    final_candidate_tier, levels = tier_order,
    labels = c("robust", "FDR cond.", "nominal", "run-order")
  )]
  metabolite_counts[, analysis_label := factor(
    unname(analysis_label_lookup[analysis_id]),
    levels = analysis_inventory$analysis_label
  )]
  p_a <- ggplot(
    metabolite_counts,
    aes(x = final_candidate_tier, y = analysis_label, fill = N)
  ) +
    geom_tile(colour = "white", linewidth = 0.28) +
    geom_text(aes(label = N), family = font_family, size = 1.95) +
    scale_fill_gradient(
      low = "#F1F2F3", high = palette_contract[["signal_blue"]]
    ) +
    labs(
      title = "PR001876 candidate tiers",
      subtitle = paste(
        "All three 16-vs-16 analyses; zero cells are true empty states"
      ),
      x = NULL, y = NULL, fill = "n"
    ) +
    theme_nature_contract() +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, size = 5.5),
      axis.line = element_blank(), axis.ticks = element_blank(),
      legend.position = "right"
    )

  display_metabolites <- metabolites[order(candidate_rank), head(.SD, 2L),
                                    by = analysis_id]
  display_metabolites <- unique(rbindlist(list(
    display_metabolites,
    metabolites[final_candidate_tier == "run_order_confounded_conditional"]
  ), fill = TRUE), by = c("analysis_id", "feature_id"))
  display_metabolites[, metabolite_label := paste0(
    short_label(metabolite_name, 22L), " · r", candidate_rank
  )]
  display_metabolites[, analysis_label := paste0(
    unname(analysis_label_lookup[analysis_id])
  )]
  display_metabolites[, analysis_label := factor(
    analysis_label, levels = analysis_inventory$analysis_label
  )]
  display_metabolites[, metabolite_label := factor(
    metabolite_label, levels = rev(unique(metabolite_label[order(candidate_rank)]))
  )]
  p_b <- ggplot(
    display_metabolites,
    aes(
      x = limma_effect_transformed, y = metabolite_label,
      colour = final_candidate_tier
    )
  ) +
    geom_vline(xintercept = 0, linewidth = 0.28, colour = "#979CA0") +
    geom_errorbar(
      aes(xmin = limma_ci95_low, xmax = limma_ci95_high),
      orientation = "y", width = 0, linewidth = 0.35
    ) +
    geom_point(size = 1.15) +
    facet_wrap(~ analysis_label, scales = "free", ncol = 1) +
    scale_colour_manual(values = c(
      robust_rank_and_scale_fdr = palette_contract[["signal_blue"]],
      conditional_fdr = palette_contract[["signal_blue_light"]],
      nominal_concordant_conditional = palette_contract[["neutral_mid"]],
      run_order_confounded_conditional = palette_contract[["accent_orange"]]
    ), drop = FALSE) +
    labs(
      title = "Analysis-specific transformed effects",
      subtitle = paste0(
        nrow(display_metabolites), "/", nrow(metabolites),
        " shown; all four GC run-order-confounded features retained"
      ),
      x = "limma effect (analysis-specific transform)", y = NULL,
      colour = "Tier"
    ) +
    theme_nature_contract() +
    theme(
      axis.text.y = element_text(size = 5.5),
      strip.text = element_text(size = 5.5),
      legend.position = "none"
    )

  robustness_source <- display_metabolites[, .(
    row_id = metabolite_label,
    `scale stable` = fifelse(
      scale_sensitivity_direction_stable %in% TRUE, "pass",
      fifelse(scale_sensitivity_direction_stable %in% FALSE, "not pass", "N/A")
    ),
    `Wilcoxon direction` = fifelse(
      direction_concordant %in% TRUE, "pass",
      fifelse(direction_concordant %in% FALSE, "not pass", "N/A")
    ),
    `paired sensitivity` = fifelse(
      !(paired_sensitivity_estimable %in% TRUE), "N/A",
      fifelse(paired_sensitivity_direction_concordant %in% TRUE, "pass", "not pass")
    ),
    `run-order` = fifelse(
      run_order_perfect_group_block %in% TRUE, "confounded",
      fifelse(
        run_order_sensitivity_direction_concordant %in% TRUE, "pass",
        fifelse(run_order_sensitivity_direction_concordant %in% FALSE,
                "not pass", "N/A")
      )
    ),
    `identity` = fifelse(
      grepl("verified", annotation_identity_status), "verified",
      "unverified/boundary"
    )
  )]
  robustness <- melt(
    robustness_source,
    id.vars = "row_id",
    variable.name = "check", value.name = "state"
  )
  robustness[, row_id := factor(row_id, levels = levels(display_metabolites$metabolite_label))]
  p_c <- ggplot(robustness, aes(x = check, y = row_id, fill = state)) +
    geom_tile(colour = "white", linewidth = 0.22) +
    scale_fill_manual(values = c(
      pass = palette_contract[["signal_blue_light"]],
      `not pass` = palette_contract[["neutral_mid"]],
      confounded = palette_contract[["accent_orange"]],
      verified = palette_contract[["signal_teal"]],
      `unverified/boundary` = palette_contract[["boundary"]],
      `N/A` = palette_contract[["neutral_pale"]]
    ), drop = FALSE) +
    labs(
      title = "Metabolite robustness matrix",
      subtitle = "Blue pass; orange run-order confounded; ochre ID unverified",
      x = NULL, y = NULL, fill = "State"
    ) +
    theme_nature_contract() +
    theme(
      axis.text.y = element_text(size = 5.5),
      axis.text.x = element_text(angle = 35, hjust = 1, size = 5.5),
      axis.line = element_blank(), axis.ticks = element_blank(),
      legend.position = "none"
    )

  pipeline <- copy(x$microbiome_pipeline)
  pipeline_long <- melt(
    pipeline,
    id.vars = "pipeline",
    measure.vars = c(
      "final_read_retention", "bray_distance_spearman_vs_main",
      "differential_direction_concordance_vs_main"
    ),
    variable.name = "metric", value.name = "value"
  )
  pipeline_long[, metric := factor(
    metric,
    levels = c(
      "final_read_retention", "bray_distance_spearman_vs_main",
      "differential_direction_concordance_vs_main"
    ),
    labels = c("Read retention", "Bray ρ vs main", "Direction concordance")
  )]
  pipeline_long[, pipeline := factor(
    pipeline, levels = unique(pipeline),
    labels = c("main 12", "min overlap 8", "forward only", "length filter")[
      seq_along(unique(pipeline))
    ]
  )]
  p_d <- ggplot(pipeline_long, aes(x = metric, y = pipeline, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.28) +
    geom_text(
      aes(label = sprintf("%.2f", value)), family = font_family, size = 1.85
    ) +
    scale_fill_gradient(
      low = "#F1F2F3", high = palette_contract[["signal_blue"]],
      limits = c(0, 1), oob = scales::squish
    ) +
    labs(
      title = "PRJNA766558 pipeline sensitivity",
      subtitle = "Cell text is the observed value",
      x = NULL, y = NULL, fill = "Value"
    ) +
    theme_nature_contract() +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1, size = 5.5),
      axis.line = element_blank(), axis.ticks = element_blank(),
      legend.position = "right"
    )

  alpha <- copy(x$microbiome_alpha)
  alpha[, metric_label := fcase(
    metric == "observed_asv", "Observed ASV",
    metric == "shannon", "Shannon",
    metric == "simpson", "Simpson",
    metric == "rarefied_observed_asv_median", "Rarefied observed",
    metric == "rarefied_shannon_median", "Rarefied Shannon",
    metric == "rarefied_simpson_median", "Rarefied Simpson",
    metric == "prokaryotic_reads_min_overlap12", "Sequencing depth",
    default = short_label(gsub("_", " ", metric), 22L)
  )]
  alpha_summary <- rbindlist(list(
    alpha[, .(metric_label, field = "Δ median",
              value = sprintf("%.2g", median_paired_difference_tumor_minus_normal))],
    alpha[, .(metric_label, field = "95% CI",
              value = paste0("[", sprintf("%.2g", effect_ci95_low), ", ",
                             sprintf("%.2g", effect_ci95_high), "]"))],
    alpha[, .(metric_label, field = "q",
              value = format_q(paired_wilcoxon_q))]
  ))
  p_e_alpha <- ggplot(alpha_summary, aes(x = field, y = metric_label)) +
    geom_tile(fill = "#EDF2F7", colour = "white", linewidth = 0.25) +
    geom_text(aes(label = value), family = font_family, size = 1.75) +
    labs(
      title = "Paired alpha effects",
      subtitle = "Original units; depth remains diagnostic only",
      x = NULL, y = NULL
    ) +
    theme_nature_contract() +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(size = 5.5)
    )

  beta <- copy(x$microbiome_beta)
  beta[, distance_label := short_label(gsub("_", " ", distance), 18L)]
  beta_summary <- rbindlist(list(
    beta[, .(distance_label, field = "R²", value = sprintf("%.3f", permanova_r2))],
    beta[, .(distance_label, field = "q", value = format_q(restricted_permutation_q))],
    beta[, .(distance_label, field = "disp. p", value = format_q(betadisper_p))],
    beta[, .(distance_label, field = "centroid p",
             value = format_q(paired_distance_to_centroid_p))]
  ))
  p_e_beta <- ggplot(beta_summary, aes(x = field, y = distance_label)) +
    geom_tile(fill = "#EDF2F7", colour = "white", linewidth = 0.25) +
    geom_text(aes(label = value), family = font_family, size = 1.75) +
    labs(
      title = "Paired beta-diversity checks",
      x = NULL, y = NULL
    ) +
    theme_nature_contract() +
    theme(axis.line = element_blank(), axis.ticks = element_blank())
  p_e <- p_e_alpha / p_e_beta + plot_layout(heights = c(1, 0.85))

  genus <- copy(x$microbiome_genus)
  genus_candidates <- genus[
    candidate_tier != "background_no_clear_paired_difference"
  ][order(
    factor(candidate_tier, levels = c(
      "paired_clr_fdr_supported_no_blank",
      "paired_clr_conditional_no_blank",
      "prevalence_shift_conditional_no_blank",
      "directional_exploratory_no_blank"
    )), paired_wilcoxon_q, paired_wilcoxon_p, -abs(median_clr_difference_tumor_minus_normal)
  )]
  total_genus_candidates <- nrow(genus_candidates)
  genus_display <- head(genus_candidates, 8L)
  ancombc2 <- copy(x$microbiome_ancombc2)
  ancom_status <- unique(ancombc2$analysis_status)
  fail_if(length(ancom_status) != 1L,
          "Figure 4 无法唯一读取 ANCOM-BC2 正交敏感性状态。")
  if (!nrow(genus_display)) {
    p_f <- empty_panel(
      "Genus-level paired CLR",
      "True empty state: no genus passed the formal non-background tiers",
      "All background rows remain in the formal source table"
    )
  } else {
    ancom_mapping_label <- "ANCOM-BC2 not completed"
    if (ancom_status == "completed_random_intercept_sensitivity") {
      ancom_map <- ancombc2[, .(
        genus_label = taxon,
        ancom_taxon_mapped = TRUE,
        ancom_effect = as.numeric(lfc_tissue_roletumor),
        ancom_q = as.numeric(q_tissue_roletumor),
        ancom_robust = as_logical_strict(
          diff_robust_tissue_roletumor,
          "ANCOM-BC2 diff_robust_tissue_roletumor"
        )
      )]
      genus_display <- merge(
        genus_display, ancom_map, by = "genus_label", all.x = TRUE,
        sort = FALSE
      )
      mapped_display_taxa <- sum(
        !is.na(genus_display$ancom_taxon_mapped) &
          genus_display$ancom_taxon_mapped
      )
      fail_if(mapped_display_taxa == 0L,
              paste(
                "ANCOM-BC2 已完成，但 Figure 4 展示 taxa 映射为 0/",
                nrow(genus_display), "；可能存在 taxon 名称漂移。"
              ))
      ancom_mapping_label <- paste0(
        "taxa mapped ", mapped_display_taxa, "/", nrow(genus_display)
      )
      genus_display[, ancom_state := fcase(
        is.na(ancom_taxon_mapped) | !ancom_taxon_mapped,
          "taxon unmapped",
        is.finite(ancom_q) & ancom_q <= 0.10 & ancom_robust &
          sign(ancom_effect) ==
            sign(median_clr_difference_tumor_minus_normal),
          "robust q≤0.10 concordant",
        is.finite(ancom_q) & ancom_q <= 0.10 & !ancom_robust,
          "q≤0.10 not pseudo-robust",
        is.finite(ancom_q) & ancom_q <= 0.10,
          "q≤0.10 discordant",
        default = "mapped; not q≤0.10 supported"
      )]
    } else {
      genus_display[, ancom_state := "ANCOM-BC2 unavailable"]
    }
    genus_display[, genus_plot_label := short_label(
      gsub("^.*g__", "", gsub("_", " ", genus_label)), 22L
    )]
    fail_if(anyDuplicated(genus_display$genus_plot_label) > 0L,
            "Figure 4 简化后的 genus 显示标签重复。")
    genus_display[, genus_plot_label := factor(
      genus_plot_label, levels = rev(genus_plot_label)
    )]
    p_f <- ggplot(
      genus_display,
      aes(
        x = median_clr_difference_tumor_minus_normal,
        y = genus_plot_label, colour = candidate_tier, shape = ancom_state
      )
    ) +
      geom_vline(xintercept = 0, linewidth = 0.28, colour = "#979CA0") +
      geom_errorbar(
        aes(xmin = clr_effect_ci95_low, xmax = clr_effect_ci95_high),
        orientation = "y", width = 0, linewidth = 0.38
      ) +
      geom_point(size = 1.25) +
      scale_colour_manual(values = c(
        paired_clr_fdr_supported_no_blank = palette_contract[["signal_blue"]],
        paired_clr_conditional_no_blank = palette_contract[["signal_blue_light"]],
        prevalence_shift_conditional_no_blank = palette_contract[["signal_teal"]],
        directional_exploratory_no_blank = palette_contract[["neutral_mid"]]
      )) +
      scale_shape_manual(values = c(
        `robust q≤0.10 concordant` = 16,
        `q≤0.10 not pseudo-robust` = 17,
        `q≤0.10 discordant` = 4,
        `mapped; not q≤0.10 supported` = 1,
        `taxon unmapped` = 8,
        `ANCOM-BC2 unavailable` = 3
      ), drop = FALSE) +
      labs(
        title = "Genus-level paired CLR candidates",
        subtitle = paste0(
          "Top ", nrow(genus_display), "/", total_genus_candidates,
          "; filled = ANCOM robust sensitivity"
        ),
        x = "Median CLR difference (tumor − paired non-tumor)", y = NULL,
        colour = "Formal tier", shape = "ANCOM-BC2 sensitivity"
      ) +
      theme_nature_contract() +
      theme(axis.text.y = element_text(size = 5.5), legend.position = "none")
  }

  microbe_sample_qc <- copy(x$microbiome_sample_qc)
  public_subset_pairs <- uniqueN(microbe_sample_qc$patient_pair_id)
  # 原发表队列上限来自正式 sample QC 的 paper_pair_number，不在图中硬编码。
  original_published_pairs <- max(
    as.integer(microbe_sample_qc$paper_pair_number), na.rm = TRUE
  )
  boundary_rows <- data.table(
    item = factor(
      c("Blank controls", "Material", "Function", "Coverage", "Sensitivity"),
      levels = rev(c(
        "Blank controls", "Material", "Function", "Coverage", "Sensitivity"
      ))
    ),
    boundary = c(
      "No public blank / negative control",
      "FFPE low biomass; contamination unresolved",
      "16S taxonomy only; function not measured",
      paste0(
        public_subset_pairs, "/", original_published_pairs,
        " published pairs; reproduction subset"
      ),
      "ANCOM sensitivity; supporting analysis only"
    )
  )
  p_microbe_boundary <- ggplot(
    boundary_rows, aes(x = 1, y = item, label = boundary)
  ) +
    geom_tile(
      fill = "#F0F1F2", colour = "white", linewidth = 0.30,
      width = 1, height = 0.92
    ) +
    geom_text(
      family = font_family, size = 1.95, hjust = 0.5,
      colour = "#3E4347"
    ) +
    scale_x_continuous(limits = c(0.5, 1.5)) +
    labs(
      title = "Microbiome evidence boundaries",
      subtitle = "Explicit limitations; none is an automatic veto",
      x = NULL, y = NULL
    ) +
    theme_nature_contract() +
    theme(
      axis.text.x = element_blank(), axis.ticks = element_blank(),
      axis.line = element_blank(),
      axis.text.y = element_text(size = 5.5, face = "bold")
    )

  metabolism_column <- p_a / p_b / p_c +
    plot_layout(heights = c(0.48, 1.42, 0.90))
  microbiome_column <- p_microbe_boundary / p_d / p_e / p_f +
    plot_layout(heights = c(0.56, 0.72, 1.12, 1.05))
  figure <- metabolism_column | plot_spacer() | microbiome_column
  figure +
    plot_layout(widths = c(1, 0.06, 1.08)) +
    plot_annotation(
      title = "Orthogonal disease-context modules",
      subtitle = paste(
        "PR001876 metabolome · no shared patient identifiers ·",
        "PRJNA766558 microbiome"
      ),
      tag_levels = "a"
    ) &
    theme(
      plot.title = element_text(size = 7.2, face = "bold", family = font_family),
      plot.subtitle = element_text(size = 5.8, colour = "#555555",
                                   family = font_family),
      plot.tag = element_text(size = 8, face = "bold", family = font_family)
    )
}

build_figure5 <- function(x) {
  drivers <- copy(x$integrated_drivers)[
    decision == "strong_patient_level_candidate"
  ]
  axes_all <- copy(x$integrated_axis_edges)
  axis_summary <- copy(x$integrated_axis_summary)
  main_genes <- axis_summary[axis_tier %chin% c("T3", "T4"), gene_name]
  boundary_genes <- c("GNAS", "ZNF750")
  focus_genes <- unique(c(main_genes, boundary_genes))

  t3_axes <- axes_all[
    evidence_tier == "T3" &
      sub("^GENE:", "", axis_candidate_id) %chin% focus_genes
  ]
  connected_factors <- unique(t3_axes[
    edge_class == "TCGA_driver_state_internal", target_node
  ])
  top_pathway_map <- copy(x$heterogeneity_axes)[
    factor %chin% connected_factors, .(source_node = factor, target_node = top_pathway)
  ]
  fail_if(anyNA(top_pathway_map) || anyDuplicated(top_pathway_map$source_node) > 0L,
          "Figure 5 连接 factor 的预冻结 top pathway 不完整。")
  adjusted_axes <- merge(
    axes_all[edge_class == "TCGA_ECMS_adjusted_factor_PROGENy" &
               evidence_tier == "T2"],
    top_pathway_map,
    by = c("source_node", "target_node"), all = FALSE, sort = FALSE
  )
  selected_pathways <- unique(top_pathway_map$target_node)
  raw_pathway_axes <- axes_all[
    edge_class == "TCGA_driver_state_internal" & evidence_tier == "T2" &
      target_layer == "RNA_derived_pathway_activity" &
      sub("^GENE:", "", axis_candidate_id) %chin% focus_genes &
      target_node %chin% selected_pathways
  ]
  boundary_axes <- axes_all[
    evidence_tier == "T2" &
      sub("^GENE:", "", axis_candidate_id) %chin% boundary_genes &
      !(edge_class == "TCGA_driver_state_internal" &
          target_layer == "RNA_derived_pathway_activity")
  ]
  axes <- unique(rbindlist(list(
    t3_axes, boundary_axes, raw_pathway_axes, adjusted_axes
  ), use.names = TRUE, fill = TRUE), by = "edge_id")
  fail_if(!nrow(axes) || any(!axes$evidence_tier %chin% c("T2", "T3")),
          "Figure 5 精简网络没有形成正式 T3 主轴/T2 边界集。")

  # 网络只显示正式 T3 主轴、GNAS/ZNF750 T2 边界和预冻结 top pathway。
  # 所有 geom_segment/geom_curve 均不设置 arrow 参数。
  driver_state <- axes[edge_class == "TCGA_driver_state_internal"]
  driver_state[, gene_name := sub("^GENE:", "", axis_candidate_id)]
  network_edges <- rbindlist(list(
    driver_state[, .(
      source = paste0("GENE:", gene_name),
      target = fifelse(
        target_layer == "RNA_derived_pathway_activity",
        paste0("PATH:", target_node), paste0("STATE:", target_node)
      ),
      edge_family = fifelse(
        target_layer == "RNA_derived_pathway_activity",
        "driver–pathway", "driver–factor"
      ),
      edge_provenance = fifelse(
        target_layer == "RNA_derived_pathway_activity",
        "raw driver–pathway (same TCGA)",
        "leakage-controlled driver–factor"
      ),
      tier = evidence_tier,
      effect = association_effect,
      q_value, association_effect_measure
    )],
    axes[edge_class == "TCGA_CNV_RNA_dosage", .(
      source = paste0("GENE:", sub("^GENE:", "", axis_candidate_id)),
      target = paste0("RNA:", sub("^GENE:", "", axis_candidate_id)),
      edge_family = "CNV–RNA",
      edge_provenance = "CNV–RNA dosage (same TCGA)",
      tier = evidence_tier,
      effect = association_effect,
      q_value, association_effect_measure
    )],
    axes[edge_class == "TCGA_ECMS_adjusted_factor_PROGENy", .(
      source = paste0("STATE:", source_node),
      target = paste0("PATH:", target_node),
      edge_family = "adjusted pathway",
      edge_provenance = "ECMS-adjusted factor–pathway",
      tier = evidence_tier,
      effect = association_effect,
      q_value, association_effect_measure
    )],
    axes[edge_class == "Cao2020_promoter_RNA", .(
      source = paste0("CAO_PROM:", sub(":promoter.*$", "", source_node)),
      target = paste0("CAO_RNA:", sub(":RNA.*$", "", target_node)),
      edge_family = "Cao promoter–RNA",
      edge_provenance = "Cao paired boundary",
      tier = evidence_tier,
      effect = association_effect,
      q_value, association_effect_measure
    )],
    axes[edge_class == "Cao2020_RNA_protein", .(
      source = paste0("CAO_RNA:", sub(":RNA.*$", "", source_node)),
      target = paste0("CAO_PROT:", sub(":protein.*$", "", target_node)),
      edge_family = "Cao RNA–protein",
      edge_provenance = "Cao paired boundary",
      tier = evidence_tier,
      effect = association_effect,
      q_value, association_effect_measure
    )]
  ), use.names = TRUE, fill = TRUE)
  # 不用 unique() 任意丢弃投影到同一 source–target 的多个正式事件。
  # 显式聚合规则：证据层级取上限、q 取最小有限值、effect 取中位数；
  # 正负号不一致时保留为 mixed，不被输入行顺序隐藏。
  network_edges[, raw_effect_direction := fcase(
    edge_family == "ECMS context", "unsigned / not directional",
    !is.finite(effect), "not estimated",
    effect > 0, "positive association",
    effect < 0, "negative association",
    default = "zero association"
  )]
  raw_network_edge_count <- nrow(network_edges)
  edge_tier_levels <- c("T2", "T3")
  network_edges <- network_edges[, {
    observed_directions <- sort(unique(raw_effect_direction))
    signed_directions <- intersect(
      observed_directions,
      c("positive association", "negative association")
    )
    aggregated_direction <- if (length(signed_directions) > 1L) {
      "mixed association signs"
    } else if (length(signed_directions) == 1L) {
      signed_directions[[1L]]
    } else if ("zero association" %chin% observed_directions) {
      "zero association"
    } else if ("unsigned / not directional" %chin% observed_directions) {
      "unsigned / not directional"
    } else {
      "not estimated"
    }
    finite_effect <- effect[is.finite(effect)]
    finite_q <- q_value[is.finite(q_value)]
    .(
      tier = edge_tier_levels[max(match(tier, edge_tier_levels))],
      effect = if (length(finite_effect)) median(finite_effect) else NA_real_,
      q_value = if (length(finite_q)) min(finite_q) else NA_real_,
      effect_direction = aggregated_direction,
      projected_event_count = .N,
      projected_tiers = paste(
        edge_tier_levels[edge_tier_levels %chin% unique(tier)],
        collapse = "+"
      )
    )
  }, by = .(
    source, target, edge_family, edge_provenance, association_effect_measure
  )]
  fail_if(sum(network_edges$projected_event_count) != raw_network_edge_count ||
            anyNA(network_edges$tier) || any(network_edges$projected_event_count < 1L),
          "Figure 5 并行投影边的确定性聚合不完整。")

  focus_driver_rows <- merge(
    drivers[gene_name %chin% focus_genes, .(gene_name)],
    axis_summary[, .(gene_name, axis_tier)],
    by = "gene_name", all.x = TRUE, sort = FALSE
  )
  base_nodes <- focus_driver_rows[, .(
    node_id = paste0("GENE:", gene_name),
    label = fcase(
      gene_name == "GNAS", "GNAS\nconditional",
      gene_name == "ZNF750", "ZNF750\nreverse",
      default = gene_name
    ),
    node_type = fcase(
      gene_name == "GNAS", "conditional boundary",
      gene_name == "ZNF750", "reverse boundary",
      default = "strong driver"
    ),
    tier = axis_tier
  )]
  edge_node_ids <- unique(c(network_edges$source, network_edges$target))
  derived_nodes <- data.table(node_id = setdiff(edge_node_ids, base_nodes$node_id))
  derived_nodes[, `:=`(
    label = fcase(
      node_id == "ECMS:official78", "ECMS official 78",
      grepl("^PATH:", node_id), sub("^PATH:", "", node_id),
      grepl("^RNA:", node_id), paste0(sub("^RNA:", "", node_id), " RNA"),
      grepl("^CAO_PROM:", node_id), paste0(sub("^CAO_PROM:", "", node_id),
                                            " promoter"),
      grepl("^CAO_RNA:", node_id), paste0(sub("^CAO_RNA:", "", node_id),
                                           " RNA"),
      grepl("^CAO_PROT:", node_id), paste0(sub("^CAO_PROT:", "", node_id),
                                            " protein"),
      grepl("^STATE:", node_id), sub("^STATE:", "", node_id),
      default = short_label(node_id, 25L)
    ),
    node_type = fcase(
      node_id == "ECMS:official78", "ECMS context",
      grepl("^PATH:", node_id), "pathway",
      grepl("^RNA:", node_id), "RNA",
      grepl("^CAO_PROM:", node_id), "Cao promoter",
      grepl("^CAO_RNA:", node_id), "Cao RNA",
      grepl("^CAO_PROT:", node_id), "Cao protein",
      grepl("^STATE:", node_id), "MOFA factor",
      default = "other"
    )
  )]
  tier_levels <- c("T0", "T1", "T2", "T3", "T4")
  edge_tier <- network_edges[, .(
    tier = tier_levels[max(match(tier, tier_levels), na.rm = TRUE)]
  ), by = .(node_id = target)]
  derived_nodes <- merge(
    derived_nodes, edge_tier, by = "node_id", all.x = TRUE, sort = FALSE
  )
  derived_nodes[is.na(tier), tier := "T2"]
  nodes <- rbindlist(list(base_nodes, derived_nodes), use.names = TRUE, fill = TRUE)

  x_lookup <- c(
    `strong driver` = 0.05,
    `conditional boundary` = 0.05,
    `reverse boundary` = 0.05,
    `Cao promoter` = 0.16,
    RNA = 0.24,
    `Cao RNA` = 0.29,
    `Cao protein` = 0.42,
    `MOFA factor` = 0.43,
    pathway = 0.84,
    other = 0.55
  )
  nodes[, x := unname(x_lookup[node_type])]
  nodes[is.na(x), x := 0.55]
  nodes[, y := {
    n <- .N
    if (n == 1L) 0.50 else seq(0.08, 0.92, length.out = n)
  }, by = node_type]
  nodes[node_type == "strong driver", y := {
    n <- .N
    if (n == 1L) 0.58 else seq(0.25, 0.92, length.out = n)
  }]
  nodes[node_type == "conditional boundary", y := 0.08]
  nodes[node_type == "reverse boundary", y := 0.16]
  # Cao 侧轨压到图底部，避免与主 driver→factor→pathway 布局伪装成同一路径。
  nodes[grepl("^Cao ", node_type), y := {
    n <- .N
    if (n == 1L) 0.06 else seq(0.03, 0.16, length.out = n)
  }, by = node_type]
  network_edges <- merge(
    network_edges,
    nodes[, .(source = node_id, x, y)],
    by = "source", all.x = TRUE, sort = FALSE
  )
  setnames(network_edges, c("x", "y"), c("x_source", "y_source"))
  network_edges <- merge(
    network_edges,
    nodes[, .(target = node_id, x, y)],
    by = "target", all.x = TRUE, sort = FALSE
  )
  setnames(network_edges, c("x", "y"), c("x_target", "y_target"))
  fail_if(anyNA(network_edges[, .(x_source, y_source, x_target, y_target)]),
          "整合网络存在无法定位的正式 T2/T3 节点。")

  edge_direction_palette <- c(
    `positive association` = palette_contract[["signal_blue"]],
    `negative association` = palette_contract[["accent_red"]],
    `mixed association signs` = "#7A5195",
    `zero association` = palette_contract[["neutral_mid"]],
    `unsigned / not directional` = palette_contract[["signal_teal"]],
    `not estimated` = palette_contract[["neutral_light"]]
  )
  edge_direction_linetype <- c(
    `positive association` = "solid",
    `negative association` = "dashed",
    `mixed association signs` = "twodash",
    `zero association` = "dotdash",
    `unsigned / not directional` = "dotted",
    `not estimated` = "longdash"
  )
  p_a <- ggplot() +
    geom_curve(
      data = network_edges,
      aes(
        x = x_source, y = y_source, xend = x_target, yend = y_target,
        colour = effect_direction, linetype = effect_direction,
        linewidth = tier, alpha = edge_provenance
      ),
      curvature = 0.08
    ) +
    geom_point(
      data = nodes[tier != "T4"],
      aes(x = x, y = y, fill = tier, shape = node_type),
      size = 2.35, colour = "#303030", stroke = 0.35
    ) +
    # T4 与 T3 使用同一填色，只增加双外框；不是“验证完成”。
    geom_point(
      data = nodes[tier == "T4"],
      aes(x = x, y = y), shape = 21, size = 3.25,
      fill = tier_palette[["T4"]], colour = "#202020", stroke = 0.45
    ) +
    geom_point(
      data = nodes[tier == "T4"],
      aes(x = x, y = y), shape = 21, size = 2.45,
      fill = tier_palette[["T4"]], colour = "white", stroke = 0.35
    ) +
    geom_text(
      data = nodes,
      aes(x = x, y = y, label = short_label(label, 20L)),
      family = font_family, size = 1.95, nudge_y = 0.022,
      check_overlap = TRUE
    ) +
    scale_fill_manual(values = tier_palette, drop = TRUE) +
    scale_colour_manual(values = edge_direction_palette, drop = TRUE) +
    scale_linetype_manual(values = edge_direction_linetype, drop = TRUE) +
    scale_linewidth_manual(values = c(T2 = 0.28, T3 = 0.50), drop = FALSE) +
    scale_alpha_manual(values = c(
      `raw driver–pathway (same TCGA)` = 0.48,
      `ECMS-adjusted factor–pathway` = 0.92,
      `leakage-controlled driver–factor` = 0.88,
      `CNV–RNA dosage (same TCGA)` = 0.78,
      `Cao paired boundary` = 0.72
    ), breaks = c(
      "Cao paired boundary", "CNV–RNA dosage (same TCGA)",
      "ECMS-adjusted factor–pathway",
      "leakage-controlled driver–factor",
      "raw driver–pathway (same TCGA)"
    ), labels = c(
      "Cao paired", "CNV–RNA", "ECMS-adjusted", "event–factor",
      "raw driver–pathway"
    ), drop = TRUE) +
    scale_shape_manual(values = c(
      `strong driver` = 21, `MOFA factor` = 22, pathway = 23,
      `conditional boundary` = 24, `reverse boundary` = 25,
      `ECMS context` = 24, RNA = 21, `Cao promoter` = 22,
      `Cao RNA` = 22, `Cao protein` = 22, other = 21
    ), breaks = c(
      "strong driver", "conditional boundary", "reverse boundary",
      "MOFA factor", "pathway", "RNA", "Cao promoter", "Cao RNA",
      "Cao protein"
    ), labels = c(
      "driver", "conditional", "reverse", "factor", "pathway", "RNA",
      "Cao promoter", "Cao RNA", "Cao protein"
    )) +
    coord_cartesian(xlim = c(0, 0.94), ylim = c(0, 1), clip = "off") +
    labs(
      title = paste(
        "Focused source-internal robust axes with explicit conditional",
        "boundaries"
      ),
      subtitle = paste0(
        "No arrows; raw and ECMS-adjusted pathway edges share\n",
        "TCGA/RNA and are not independent."
      ),
      colour = "Direction", linetype = "Direction",
      linewidth = "Evidence tier", alpha = "Edge source",
      fill = "Node role", shape = "Node"
    ) +
    theme_void(base_family = font_family, base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 0.5, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.4, colour = "#555555"),
      legend.position = "bottom", legend.box = "vertical",
      legend.title = element_text(size = base_size - 0.3, face = "bold"),
      legend.text = element_text(size = base_size - 0.7),
      plot.margin = margin(4, 5, 4, 5, unit = "mm")
    ) +
    guides(
      linewidth = "none",
      colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1),
      alpha = guide_legend(nrow = 2, byrow = TRUE),
      fill = "none", shape = guide_legend(nrow = 3, byrow = TRUE)
    )

  decision <- merge(
    drivers[, .(
      gene_name, integrated_tier, state_supported_edges,
      cao_cross_layer_class
    )],
    axis_summary[, .(gene_name, axis_tier, axis_class)],
    by = "gene_name", all.x = TRUE, sort = FALSE
  )
  decision[is.na(axis_tier), `:=`(
    axis_tier = "T1", axis_class = "no integrated axis"
  )]
  decision_rows <- rbindlist(list(
    decision[, .(
      gene_name, field = "driver role", state = integrated_tier
    )],
    decision[, .(
      gene_name, field = "robust bridge",
      state = fifelse(state_supported_edges > 0L, "pass", "no pass")
    )],
    decision[, .(
      gene_name, field = "Cao calibration",
      state = fcase(
        cao_cross_layer_class == "three_layer_directional_hypothesis",
          "directional",
        cao_cross_layer_class == "credible_reverse_retained", "reverse",
        default = "context/weak"
      )
    )],
    decision[, .(
      gene_name, field = "axis role", state = axis_tier
    )]
  ))
  decision_rows[, gene_name := factor(gene_name, levels = rev(drivers$gene_name))]
  decision_rows[, field := factor(
    field,
    levels = c(
      "driver role", "robust bridge", "Cao calibration", "axis role"
    )
  )]
  decision_palette <- c(
    tier_palette,
    pass = palette_contract[["signal_blue"]],
    `no pass` = palette_contract[["neutral_light"]],
    directional = palette_contract[["signal_blue_light"]],
    reverse = "#D7A6A3", `context/weak` = "#E7E8E9"
  )
  decision_rows[, state_code := fcase(
    state == "T4", "priority", state == "T3", "robust",
    state == "T2", "conditional", state == "T1", "context",
    state == "T0", "technical",
    state == "pass", "✓", state == "no pass", "×",
    state == "directional", "DIR", state == "reverse", "REV",
    default = "CTX"
  )]
  p_b <- ggplot(decision_rows, aes(x = field, y = gene_name, fill = state)) +
    geom_tile(colour = "white", linewidth = 0.28) +
    geom_text(
      aes(label = state_code), family = font_family, size = 1.70,
      colour = "#303438"
    ) +
    scale_fill_manual(values = decision_palette, drop = FALSE) +
    labs(
      title = "All 12 strong-driver decisions",
      subtitle = "Rows without a passed bridge remain visible",
      x = NULL, y = NULL, fill = "Categorical state",
      caption = "DIR directional; REV reverse; CTX context/weak."
    ) +
    theme_nature_contract() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.line = element_blank(), axis.ticks = element_blank(),
      legend.position = "none", plot.caption = element_text(size = 5.5, hjust = 0)
    )

  heterogeneity <- copy(x$heterogeneity_axes)
  factor_rows <- rbindlist(list(
    heterogeneity[, .(
      factor, field = "source bridge",
      state = fifelse(supported_driver_edge_count > 0L, "present", "absent")
    )],
    heterogeneity[, .(
      factor, field = "ECMS context",
      state = fifelse(ecms_context_tier == "T2", "present", "absent")
    )],
    heterogeneity[, .(
      factor, field = "adj. pathway",
      state = fifelse(ecms_adjusted_pathway_t2_count > 0L, "present", "absent")
    )],
    heterogeneity[, .(
      factor, field = "technical flag",
      state = fifelse(level_factor_soft_flag, "T0 level-factor", "not level-factor")
    )],
    heterogeneity[, .(
      factor, field = "combined role", state = evidence_tier
    )]
  ))
  factor_rows[, factor := factor(factor, levels = rev(paste0("Factor", 1:8)))]
  factor_rows[, field := factor(
    field,
    levels = c(
      "source bridge", "ECMS context", "adj. pathway", "technical flag",
      "combined role"
    )
  )]
  factor_rows[, state_code := fcase(
    state == "T4", "priority", state == "T3", "robust",
    state == "T2", "conditional", state == "T1", "context",
    state == "T0", "technical",
    state == "present", "✓", state == "absent", "×",
    state == "T0 level-factor", "technical", default = "—"
  )]
  p_c <- ggplot(factor_rows, aes(x = field, y = factor, fill = state)) +
    geom_tile(colour = "white", linewidth = 0.28) +
    geom_text(
      aes(label = state_code), family = font_family, size = 1.75,
      colour = "#404448"
    ) +
    scale_fill_manual(values = c(
      present = palette_contract[["signal_blue_light"]],
      absent = palette_contract[["neutral_light"]],
      `T0 level-factor` = tier_palette[["T0"]],
      `not level-factor` = "#F2F3F4",
      tier_palette
    ), drop = FALSE) +
    labs(
      title = "Factor-level composite evidence",
      subtitle = "Factor4 remains technical/background",
      x = NULL, y = NULL, fill = "Categorical state"
    ) +
    theme_nature_contract() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.line = element_blank(), axis.ticks = element_blank(),
      legend.position = "none"
    )

  modules <- copy(x$module_summaries)[module_id %chin% c(
    "PR001876_METABOLOME", "PRJNA766558_MICROBIOME"
  )]
  modules[, module_order := match(module_id, c(
    "PR001876_METABOLOME", "PRJNA766558_MICROBIOME"
  ))]
  setorder(modules, module_order)
  modules[, retained_n := as.integer(sub("^([0-9]+).*$", "\\1", module_result))]
  satellite_lines <- data.table(
    y = c(0.75, 0.50, 0.22),
    label = c(
      paste0(modules[1L, dataset], " · ", modules[1L, retained_n],
             " retained features · no host link"),
      paste0(modules[2L, dataset], " · ", modules[2L, retained_n],
             " retained genera · no host link"),
      "No shared patient identifiers · no connecting edge"
    ),
    fill = c("T2", "T2", "boundary")
  )
  p_d <- ggplot(satellite_lines, aes(x = 0.5, y = y, label = label)) +
    geom_tile(
      aes(fill = fill), width = 0.96, height = 0.20,
      colour = "white", linewidth = 0.30
    ) +
    geom_text(family = font_family, size = 1.82) +
    scale_fill_manual(values = c(
      T2 = tier_palette[["T2"]], boundary = palette_contract[["neutral_pale"]]
    )) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
    labs(
      title = "Disconnected orthogonal modules",
      subtitle = "No patient links; full boundaries in source table."
    ) +
    theme_void(base_family = font_family, base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 0.5, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.4, colour = "#555555"),
      legend.position = "none",
      plot.margin = margin(2.2, 2.2, 2.2, 2.2, unit = "mm")
    )

  right_column <- p_b / p_c / p_d +
    plot_layout(heights = c(1, 1, 0.72))
  figure <- p_a | right_column
  figure +
    plot_layout(widths = c(1.55, 1)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 8, face = "bold", family = font_family))
}

# Positive-story main figures. The full failure matrices and execution diagnostics
# remain auditable in the formal source tables; these builders only change the
# manuscript-facing visual hierarchy and never recompute upstream inference.
build_figure1_positive <- function(x) {
  retained <- copy(x$leakage_candidate_summary)[countable_for_T3_T4 == TRUE]
  fail_if(nrow(retained) != 9L || uniqueN(retained$original_edge_id) != 9L,
          "Figure 1 阳性主视觉必须恰好包含 9 条保留 event-state 边。")
  fail_if(any(retained$reference_factor == "Factor4") ||
            any(!retained$leakage_gate_pass),
          "Figure 1 阳性主视觉混入 Factor4 或未通过泄漏门禁的边。")

  retained[, event_label := fcase(
    event_type == "relative_cnv", "relative CNV",
    event_type == "amplification", "amplification",
    event_type == "mutation", "mutation",
    default = gsub("_", " ", event_type)
  )]
  retained[, effect_direction := fifelse(
    original_effect >= 0, "positive", "negative"
  )]
  retained[, edge_label := paste0(
    gene_name, " · ", event_label, " → ", reference_factor,
    " · q=", format_q(original_q_value)
  )]
  edge_order <- retained[order(-abs(original_effect)), edge_label]
  retained[, edge_label := factor(edge_label, levels = rev(edge_order))]

  sign_palette <- c(
    positive = palette_contract[["signal_blue"]],
    negative = palette_contract[["accent_red"]]
  )
  sign_shapes <- c(positive = 16, negative = 17)
  p_a <- ggplot(
    retained,
    aes(x = original_effect, y = edge_label, colour = effect_direction)
  ) +
    geom_vline(xintercept = 0, colour = "#AEB3B8", linewidth = 0.30) +
    geom_segment(
      aes(x = 0, xend = original_effect, yend = edge_label),
      linewidth = 0.55
    ) +
    geom_point(aes(shape = effect_direction), size = 2.15) +
    geom_point(
      data = retained[gene_name == "PIK3CA"],
      shape = 21, size = 3.15, fill = NA, colour = "#1F1F1F",
      stroke = 0.65
    ) +
    scale_colour_manual(values = sign_palette) +
    scale_shape_manual(values = sign_shapes) +
    coord_cartesian(xlim = c(-1.05, 0.86), clip = "off") +
    labs(
      title = "Nine retained event–state associations",
      subtitle = "PIK3CA–Factor1 CNV-family edges are double outlined",
      x = "Original source-table association effect", y = NULL,
      colour = "Direction", shape = "Direction"
    ) +
    theme_nature_contract() +
    theme(
      legend.position = "bottom",
      axis.text.y = element_text(size = 5.45),
      plot.margin = margin(2.2, 2.2, 2.2, 3.5, unit = "mm")
    ) +
    guides(colour = "none", shape = guide_legend(nrow = 1))

  retention_columns <- c(
    "drop_both_direction_retention_rate",
    "drop_both_magnitude_retention_rate",
    "drop_both_gate_retention_rate"
  )
  retained[, (retention_columns) := lapply(.SD, as.numeric),
           .SDcols = retention_columns]
  retention <- melt(
    retained,
    id.vars = c("original_edge_id", "edge_label", "planned_seed_count"),
    measure.vars = retention_columns,
    variable.name = "retention_metric", value.name = "retention_rate"
  )
  retention[, retention_metric := factor(
    retention_metric,
    levels = c(
      "drop_both_direction_retention_rate",
      "drop_both_magnitude_retention_rate",
      "drop_both_gate_retention_rate"
    ),
    labels = c("Direction", "≥50% magnitude", "Support")
  )]
  retention[, retained_seed_label := paste0(
    round(retention_rate * planned_seed_count), "/", planned_seed_count
  )]
  p_b <- ggplot(
    retention,
    aes(x = retention_metric, y = edge_label, fill = retention_rate)
  ) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_text(
      aes(label = retained_seed_label), family = font_family, size = 1.85,
      colour = "#202428"
    ) +
    scale_fill_gradient(
      low = "#E3E7EA", high = palette_contract[["signal_blue"]],
      limits = c(0, 1), breaks = c(2 / 3, 1),
      labels = c("2/3", "3/3")
    ) +
    labs(
      title = "Three-seed drop-both retention",
      subtitle = "Direction, magnitude and inferential support",
      x = NULL, y = NULL, fill = "Rate"
    ) +
    theme_nature_contract() +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.line = element_blank(), legend.position = "bottom"
    )

  pik_driver <- copy(x$integrated_drivers)[
    decision == "strong_patient_level_candidate" & gene_name == "PIK3CA"
  ]
  fail_if(nrow(pik_driver) != 1L ||
            anyNA(pik_driver[, .(
              spearman_rho_relative_cnv, q_relative_cnv
            )]),
          "Figure 1 缺少唯一 PIK3CA CNV–RNA dosage 行。")
  pik_edges <- retained[gene_name == "PIK3CA"]
  fail_if(nrow(pik_edges) != 2L ||
            !setequal(pik_edges$event_type, c("relative_cnv", "amplification")),
          "Figure 1 PIK3CA–Factor1 必须包含 relative CNV 与 amplification 两行。")
  pik_evidence <- rbindlist(list(
    pik_edges[, .(
      evidence = paste0(event_label, " → Factor1"),
      association_effect = original_effect,
      q_value = original_q_value,
      evidence_class = "event–state"
    )],
    pik_driver[, .(
      evidence = "relative CNV → PIK3CA RNA",
      association_effect = spearman_rho_relative_cnv,
      q_value = q_relative_cnv,
      evidence_class = "CNV–RNA dosage"
    )]
  ))
  pik_evidence[, evidence_label := paste0(
    evidence, "\nq=", format_q(q_value)
  )]
  pik_evidence[, evidence_label := factor(
    evidence_label, levels = rev(evidence_label[order(association_effect)])
  )]
  p_c <- ggplot(
    pik_evidence,
    aes(x = association_effect, y = evidence_label, colour = evidence_class)
  ) +
    geom_segment(
      aes(x = 0, xend = association_effect, yend = evidence_label),
      linewidth = 0.65
    ) +
    geom_point(size = 2.45) +
    geom_point(
      shape = 21, size = 3.35, fill = NA, colour = "#202020",
      stroke = 0.55
    ) +
    scale_colour_manual(values = c(
      `event–state` = palette_contract[["signal_blue"]],
      `CNV–RNA dosage` = palette_contract[["signal_teal"]]
    )) +
    coord_cartesian(xlim = c(0, 0.78), clip = "off") +
    labs(
      title = "PIK3CA–Factor1 is the leading within-TCGA priority",
      subtitle = paste(
        "Amplification and relative CNV provide complementary estimates",
        "from one CNV event family"
      ),
      x = "Association effect from the formal source row", y = NULL,
      colour = "Evidence"
    ) +
    theme_nature_contract() +
    theme(legend.position = "bottom", axis.text.y = element_text(size = 5.5)) +
    guides(colour = guide_legend(nrow = 1))

  (p_a | p_b) / p_c +
    plot_layout(heights = c(1.55, 0.72), widths = c(1.55, 0.82)) +
    plot_annotation(
      title = "Representation-overlap-audited event–state bridges",
      subtitle = paste(
        "Nine retained associations nominate PIK3CA–Factor1",
        "as the leading CNV-linked axis"
      ),
      tag_levels = "a"
    ) &
    theme(
      plot.title = element_text(size = 7.3, face = "bold", family = font_family),
      plot.subtitle = element_text(size = 5.7, colour = "#555555",
                                   family = font_family),
      plot.tag = element_text(size = 8, face = "bold", family = font_family)
    )
}

build_figure2_positive <- function(x) {
  focus_factors <- c("Factor1", "Factor3")
  official <- copy(x$ecms_probabilities)[
    in_official_78 & eligible_for_primary_association,
    .(patient_id, official_anchor_label)
  ]
  fail_if(nrow(official) != 78L || uniqueN(official$patient_id) != 78L,
          "Figure 2 主视觉必须且只能使用 official 78。")
  official <- merge(
    official,
    x$factor_scores[, c("patient_id", focus_factors), with = FALSE],
    by = "patient_id", all.x = TRUE, sort = FALSE
  )
  fail_if(anyNA(official[, ..focus_factors]),
          "Figure 2 official 78 的 Factor1/Factor3 分数不完整。")

  eta <- copy(x$ecms_factor)[
    is_primary_scope == TRUE & factor %chin% focus_factors
  ]
  fail_if(nrow(eta) != 2L || any(eta$n_patients != 78L),
          "Figure 2 official ECMS Factor1/Factor3 统计行不完整。")
  eta[, factor_panel := paste0(
    factor, "\nη²=", sprintf("%.3f", eta_squared),
    "; q=", format_q(anova_q_value)
  )]
  eta[, factor_panel := factor(factor_panel, levels = eta$factor_panel)]

  score_long <- melt(
    official,
    id.vars = c("patient_id", "official_anchor_label"),
    measure.vars = focus_factors,
    variable.name = "factor", value.name = "factor_score"
  )
  score_long <- merge(
    score_long, eta[, .(factor, factor_panel)],
    by = "factor", all.x = TRUE, sort = FALSE
  )
  score_long[, official_anchor_label := factor(
    official_anchor_label, levels = paste0("ECMS", 1:4)
  )]
  ecms_palette <- c(
    ECMS1 = "#9AB9D9", ECMS2 = "#6CAFA7",
    ECMS3 = "#D4A16A", ECMS4 = "#B49AC8"
  )
  p_a <- ggplot(
    score_long,
    aes(x = official_anchor_label, y = factor_score, fill = official_anchor_label)
  ) +
    geom_violin(
      scale = "width", trim = FALSE, alpha = 0.58,
      colour = "#6E7378", linewidth = 0.28
    ) +
    geom_boxplot(
      width = 0.18, outlier.shape = NA, fill = "white",
      colour = "#303438", linewidth = 0.32
    ) +
    geom_jitter(
      width = 0.10, height = 0, size = 0.52, alpha = 0.52,
      colour = "#303438"
    ) +
    facet_wrap(~factor_panel, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = ecms_palette) +
    labs(
      title = "Factor1 and Factor3 vary continuously across official ECMS anchors",
      subtitle = "Independent patient points; official anchor labels are explanatory context",
      x = NULL, y = "MOFA factor score", fill = "Official ECMS"
    ) +
    theme_nature_contract() +
    theme(legend.position = "none")

  eta[, factor := factor(factor, levels = rev(focus_factors))]
  p_b <- ggplot(eta, aes(x = eta_squared, y = factor)) +
    geom_segment(
      aes(x = 0, xend = eta_squared, yend = factor),
      colour = palette_contract[["signal_blue"]], linewidth = 0.70
    ) +
    geom_point(
      aes(size = -log10(anova_q_value)),
      colour = palette_contract[["signal_blue"]]
    ) +
    geom_text(
      aes(label = paste0("η²=", sprintf("%.3f", eta_squared),
                         "\nq=", format_q(anova_q_value))),
      hjust = -0.12, family = font_family, size = 1.72
    ) +
    scale_size_continuous(range = c(2.0, 3.2), guide = "none") +
    coord_cartesian(xlim = c(0, 0.43), clip = "off") +
    labs(
      title = "ECMS explanatory increment",
      x = "Eta-squared", y = NULL
    ) +
    theme_nature_contract() +
    theme(plot.margin = margin(2.2, 8, 2.2, 2.2, unit = "mm"))

  formal_t2 <- unique(x$integrated_axis_edges[
    edge_class == "TCGA_ECMS_adjusted_factor_PROGENy" &
      evidence_tier == "T2" &
      source_node %chin% focus_factors,
    .(factor = source_node, pathway = target_node)
  ])
  pathway <- merge(
    x$ecms_adjusted[
      is_primary_scope == TRUE & factor %chin% focus_factors
    ],
    formal_t2,
    by = c("factor", "pathway"), all = FALSE, sort = FALSE
  )
  fail_if(!nrow(pathway) || any(pathway$n_patients != 78L),
          "Figure 2 缺少 official78 的正式 T2 ECMS-adjusted pathway 行。")
  setorder(pathway, factor, incremental_q_value, -partial_r_squared, pathway)
  pathway <- pathway[, head(.SD, 5L), by = factor]
  pathway[, ci_low := factor_beta_standardized -
    1.96 * factor_standard_error]
  pathway[, ci_high := factor_beta_standardized +
    1.96 * factor_standard_error]
  pathway[, pathway_label := paste0(
    pathway, " · q=", format_q(incremental_q_value)
  )]
  pathway_label_order <- pathway[
    order(factor, factor_beta_standardized), unique(pathway_label)
  ]
  pathway[, pathway_label := factor(
    pathway_label,
    levels = rev(pathway_label_order)
  )]
  p_c <- ggplot(
    pathway,
    aes(
      x = factor_beta_standardized, y = pathway_label,
      colour = factor, shape = factor
    )
  ) +
    geom_vline(xintercept = 0, colour = "#AEB3B8", linewidth = 0.30) +
    geom_errorbar(
      aes(xmin = ci_low, xmax = ci_high),
      orientation = "y", width = 0.15, linewidth = 0.38
    ) +
    geom_point(aes(size = partial_r_squared), stroke = 0.35) +
    facet_wrap(~factor, scales = "free_y", nrow = 1) +
    scale_colour_manual(values = c(
      Factor1 = palette_contract[["signal_teal"]],
      Factor3 = palette_contract[["signal_blue"]]
    )) +
    scale_shape_manual(values = c(Factor1 = 16, Factor3 = 17)) +
    scale_size_continuous(range = c(1.7, 4.0), guide = "none") +
    labs(
      title = "ECMS-adjusted pathway associations",
      subtitle = "Top formal T2 rows per factor; point size = partial R²",
      x = "Standardized adjusted beta", y = NULL,
      colour = "Factor", shape = "Factor"
    ) +
    theme_nature_contract() +
    theme(
      legend.position = "bottom", axis.text.y = element_text(size = 5.1)
    ) +
    guides(colour = "none", shape = guide_legend(nrow = 1))

  summary_rows <- copy(x$heterogeneity_axes)[factor %chin% focus_factors]
  fail_if(nrow(summary_rows) != 2L,
          "Figure 2 缺少 Factor1/Factor3 最终异质性摘要。")
  summary_cards <- rbindlist(list(
    summary_rows[, .(
      factor, field = "Retained\nevent edges",
      value = as.character(supported_driver_edge_count)
    )],
    summary_rows[, .(
      factor, field = "ECMS\nη²",
      value = sprintf("%.3f", ecms_eta_squared)
    )],
    summary_rows[, .(
      factor, field = "Adjusted\nT2 pathways",
      value = as.character(ecms_adjusted_pathway_t2_count)
    )]
  ))
  summary_cards[, factor := factor(factor, levels = rev(focus_factors))]
  summary_cards[, field := factor(
    field,
    levels = c("Retained\nevent edges", "ECMS\nη²",
               "Adjusted\nT2 pathways")
  )]
  p_d <- ggplot(summary_cards, aes(x = field, y = factor)) +
    geom_tile(fill = "#EEF3F8", colour = "white", linewidth = 0.45) +
    geom_text(
      aes(label = value), family = font_family, size = 2.20,
      colour = "#1F405F", fontface = "bold"
    ) +
    labs(
      title = "Continuous-axis evidence summary",
      subtitle = "Counts and effects retain their original units",
      x = NULL, y = NULL
    ) +
    theme_nature_contract() +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(size = 5.2)
    )

  p_a / (p_b | p_c | p_d) +
    plot_layout(heights = c(1.18, 1), widths = c(0.72, 1.55, 0.78)) +
    plot_annotation(
      title = "Factor1 and Factor3 organize continuous ESCC states",
      subtitle = paste(
        "ECMS-anchored distributions, effect sizes and adjusted pathways",
        "in 78 tumors"
      ),
      tag_levels = "a"
    ) &
    theme(
      plot.title = element_text(size = 7.3, face = "bold", family = font_family),
      plot.subtitle = element_text(size = 5.7, colour = "#555555",
                                   family = font_family),
      plot.tag = element_text(size = 8, face = "bold", family = font_family)
    )
}

build_figure3_positive <- function(x) {
  focus_genes <- c("GNAS", "ZNF750")
  summary <- copy(x$cao_summary)[tcga_gene_name %chin% focus_genes]
  patient <- copy(x$cao_patient)[tcga_gene_name %chin% focus_genes]
  fail_if(nrow(summary) != 2L || !setequal(summary$tcga_gene_name, focus_genes),
          "Figure 3 必须包含 GNAS 与 ZNF750 两个正式 Cao 摘要行。")

  patient_strip <- function(gene, include_protein) {
    gene_patient <- patient[tcga_gene_name == gene]
    layers <- rbindlist(list(
      gene_patient[is.finite(median_delta_beta_promoter), .(
        patient_id = paper_patient_id,
        layer = "Promoter Δβ",
        effect = median_delta_beta_promoter
      )],
      gene_patient[is.finite(rna_delta_log2_tpm_pc1), .(
        patient_id = paper_patient_id,
        layer = "RNA Δlog2 TPM",
        effect = rna_delta_log2_tpm_pc1
      )],
      if (include_protein) {
        gene_patient[is.finite(protein_log2_tumor_vs_normal), .(
          patient_id = paper_patient_id,
          layer = "Protein log2 T/N",
          effect = protein_log2_tumor_vs_normal
        )]
      } else {
        NULL
      }
    ), use.names = TRUE, fill = TRUE)
    layer_summary <- layers[, .(
      n = .N, median_effect = median(effect)
    ), by = layer]
    layers <- merge(layers, layer_summary, by = "layer", sort = FALSE)
    layers[, panel := paste0(
      layer, "\nn=", n, "; median=", sprintf("%+.3f", median_effect)
    )]
    panel_order <- c("Promoter Δβ", "RNA Δlog2 TPM", "Protein log2 T/N")
    observed_order <- vapply(
      panel_order[panel_order %chin% layers$layer],
      function(layer_name) unique(layers[layer == layer_name, panel]),
      character(1)
    )
    layers[, panel := factor(panel, levels = observed_order)]
    ggplot(layers, aes(x = effect, y = 1)) +
      geom_vline(xintercept = 0, colour = "#AEB3B8", linewidth = 0.30) +
      geom_jitter(
        width = 0, height = 0.07, shape = 21, size = 1.5,
        fill = palette_contract[["signal_blue_light"]],
        colour = "#3F5368", stroke = 0.35, alpha = 0.82
      ) +
      geom_point(
        data = unique(layers[, .(panel, median_effect)]),
        aes(x = median_effect, y = 1), inherit.aes = FALSE,
        shape = 23, size = 2.7, fill = palette_contract[["signal_blue"]],
        colour = "#202020", stroke = 0.45
      ) +
      facet_wrap(~panel, nrow = 1, scales = "free_x") +
      scale_y_continuous(NULL, breaks = NULL) +
      labs(
        title = paste0(gene, " patient-paired layer effects"),
        subtitle = paste(
          "Patient-level tumour-minus-paired-non-tumour differences;",
          "each layer retains its native scale"
        ),
        x = "Tumour − paired non-tumour effect"
      ) +
      theme_nature_contract() +
      theme(panel.spacing = grid::unit(2.2, "mm"))
  }

  p_a <- patient_strip("GNAS", include_protein = TRUE)
  p_b <- patient_strip("ZNF750", include_protein = FALSE)

  summary_long <- rbindlist(list(
    summary[, .(
      gene = tcga_gene_name, layer = "Promoter",
      effect = promoter_median_paired_delta_beta,
      q_value = promoter_q_candidate_by_region,
      availability = promoter_region_evaluation_status
    )],
    summary[, .(
      gene = tcga_gene_name, layer = "RNA",
      effect = limma_log_fc_pc1,
      q_value = limma_q_candidate_pc1,
      availability = "evaluable"
    )],
    summary[, .(
      gene = tcga_gene_name, layer = "Protein",
      effect = median_patient_log2_tumor_vs_normal,
      q_value = q_candidate,
      availability = protein_measurement_status
    )]
  ))
  summary_long[, direction_class := fcase(
    !is.finite(effect), "unavailable",
    effect > 0, "positive",
    effect < 0, "negative",
    default = "zero"
  )]
  summary_long[, cell_label := fcase(
    !is.finite(effect), "not identified",
    default = paste0(
      sprintf("%+.3f", effect), "\nq=", format_q(q_value)
    )
  )]
  summary_long[, gene_label := fcase(
    gene == "GNAS", "GNAS · three-layer directional hypothesis",
    gene == "ZNF750", "ZNF750 · falsifiable reverse pattern"
  )]
  summary_long[, gene_label := factor(
    gene_label,
    levels = rev(c(
      "GNAS · three-layer directional hypothesis",
      "ZNF750 · falsifiable reverse pattern"
    ))
  )]
  summary_long[, layer := factor(layer, levels = c("Promoter", "RNA", "Protein"))]
  p_c <- ggplot(
    summary_long,
    aes(x = layer, y = gene_label, fill = direction_class)
  ) +
    geom_tile(colour = "white", linewidth = 0.55) +
    geom_text(
      aes(label = cell_label), family = font_family, size = 1.75,
      lineheight = 0.92, colour = "#25292D"
    ) +
    scale_fill_manual(values = c(
      positive = palette_contract[["signal_blue_light"]],
      negative = "#E2B2AF",
      zero = palette_contract[["neutral_light"]],
      unavailable = palette_contract[["neutral_pale"]]
    )) +
    labs(
      title = "Direction summary",
      subtitle = "Categorical direction only; values use layer-specific units",
      x = NULL, y = NULL, fill = "Direction"
    ) +
    theme_nature_contract() +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text.y = element_text(size = 5.1),
      legend.position = "bottom"
    ) +
    guides(fill = guide_legend(nrow = 1))

  p_a / (p_b | p_c) +
    plot_layout(heights = c(1, 1), widths = c(1.05, 1)) +
    plot_annotation(
      title = "Same-patient cross-layer calibration",
      subtitle = paste(
        "GNAS supplies a directional three-layer hypothesis;",
        "ZNF750 preserves an informative reverse pattern"
      ),
      tag_levels = "a"
    ) &
    theme(
      plot.title = element_text(size = 7.3, face = "bold", family = font_family),
      plot.subtitle = element_text(size = 5.7, colour = "#555555",
                                   family = font_family),
      plot.tag = element_text(size = 8, face = "bold", family = font_family)
    )
}

build_figure4_positive <- function(x) {
  metabolites <- copy(x$metabolite_candidates)
  all_expected_tiers <- c(
    "robust_rank_and_scale_fdr",
    "conditional_fdr",
    "nominal_concordant_conditional",
    "run_order_confounded_conditional"
  )
  tier_counts <- metabolites[, .N, by = final_candidate_tier]
  tier_counts <- merge(
    data.table(final_candidate_tier = all_expected_tiers),
    tier_counts, by = "final_candidate_tier", all.x = TRUE, sort = FALSE
  )
  tier_counts[is.na(N), N := 0L]
  expected_counts <- c(83L, 43L, 5L, 4L)
  fail_if(!identical(as.integer(tier_counts$N), expected_counts),
          "Figure 4 PR001876 候选层级计数不是 83/43/5/4。")
  display_tiers <- all_expected_tiers[1:3]
  tier_counts <- tier_counts[final_candidate_tier %chin% display_tiers]
  tier_counts[, tier_label := factor(
    final_candidate_tier,
    levels = rev(display_tiers),
    labels = rev(c(
      "Robust rank + scale FDR",
      "Conditional FDR",
      "Nominal concordant"
    ))
  )]
  p_a <- ggplot(
    tier_counts,
    aes(x = N, y = tier_label, fill = final_candidate_tier)
  ) +
    geom_col(width = 0.64) +
    geom_text(
      aes(label = N), hjust = -0.25, family = font_family,
      size = 2.05, fontface = "bold"
    ) +
    scale_fill_manual(values = c(
      robust_rank_and_scale_fdr = palette_contract[["signal_blue"]],
      conditional_fdr = palette_contract[["signal_blue_light"]],
      nominal_concordant_conditional = palette_contract[["neutral_mid"]]
    )) +
    coord_cartesian(xlim = c(0, 94), clip = "off") +
    labs(
      title = "PR001876 candidate evidence tiers",
      subtitle = "83 robust LC-MS analysis-features",
      x = "Retained analysis-features", y = NULL
    ) +
    theme_nature_contract() +
    theme(legend.position = "none")

  robust <- metabolites[
    final_candidate_tier == "robust_rank_and_scale_fdr"
  ][order(limma_q_evidence_family, candidate_rank)][1:10]
  fail_if(nrow(robust) != 10L ||
            any(robust$platform != "targeted LC-MS") ||
            any(!robust$scale_sensitivity_direction_stable),
          "Figure 4 top-10 robust LC-MS 行不完整或混入非稳健/非 LC-MS 行。")
  robust[, direction_class := fifelse(
    limma_effect_transformed >= 0, "higher in ESCC", "lower in ESCC"
  )]
  robust[, metabolite_label := paste0(
    short_label(metabolite_name, 27L),
    " · q=", format_q(limma_q_evidence_family)
  )]
  robust[, metabolite_label := factor(
    metabolite_label,
    levels = rev(robust[order(limma_effect_transformed), metabolite_label])
  )]
  direction_palette <- c(
    `higher in ESCC` = palette_contract[["signal_blue"]],
    `lower in ESCC` = palette_contract[["accent_red"]]
  )
  p_b <- ggplot(
    robust,
    aes(
      x = limma_effect_transformed, y = metabolite_label,
      colour = direction_class, shape = direction_class
    )
  ) +
    geom_vline(xintercept = 0, colour = "#AEB3B8", linewidth = 0.30) +
    geom_errorbar(
      aes(xmin = limma_ci95_low, xmax = limma_ci95_high),
      orientation = "y", width = 0.15, linewidth = 0.38
    ) +
    geom_point(size = 2.05) +
    scale_colour_manual(values = direction_palette) +
    scale_shape_manual(values = c(`higher in ESCC` = 16, `lower in ESCC` = 17)) +
    labs(
      title = "Top robust LC-MS effects",
      subtitle = "Primary-scale effects for source-labelled LC-MS features",
      x = "Limma effect on the primary transformed scale", y = NULL,
      colour = "Direction", shape = "Direction"
    ) +
    theme_nature_contract() +
    theme(
      legend.position = "bottom", axis.text.y = element_text(size = 5.0)
    ) +
    guides(colour = "none", shape = guide_legend(nrow = 1))

  alpha <- copy(x$microbiome_alpha)[
    metric %chin% c("observed_asv", "shannon", "simpson")
  ]
  fail_if(nrow(alpha) != 3L || any(alpha$n_pairs != 21L),
          "Figure 4 paired alpha-diversity 必须为 observed ASV/Shannon/Simpson 三行。")
  alpha[, metric_label := fcase(
    metric == "observed_asv", "Observed ASVs",
    metric == "shannon", "Shannon",
    metric == "simpson", "Simpson"
  )]
  alpha[, panel := paste0(
    metric_label, "\nq=", format_q(paired_wilcoxon_q)
  )]
  alpha[, panel := factor(
    panel,
    levels = alpha[match(c("observed_asv", "shannon", "simpson"), metric), panel]
  )]
  p_c <- ggplot(
    alpha,
    aes(x = median_paired_difference_tumor_minus_normal, y = 1)
  ) +
    geom_vline(xintercept = 0, colour = "#AEB3B8", linewidth = 0.30) +
    geom_errorbar(
      aes(xmin = effect_ci95_low, xmax = effect_ci95_high),
      orientation = "y", width = 0.16, linewidth = 0.45,
      colour = palette_contract[["accent_red"]]
    ) +
    geom_point(
      shape = 17, size = 2.25, colour = palette_contract[["accent_red"]]
    ) +
    facet_wrap(~panel, ncol = 1, scales = "free_x") +
    scale_y_continuous(NULL, breaks = NULL) +
    labs(
      title = "Paired diversity shifts",
      subtitle = "Tumour − paired non-tumour; n=21 pairs",
      x = "Median paired difference (metric-specific unit)"
    ) +
    theme_nature_contract() +
    theme(panel.spacing = grid::unit(1.5, "mm"))

  genus <- copy(x$microbiome_genus)[
    candidate_tier == "paired_clr_fdr_supported_no_blank"
  ][order(paired_wilcoxon_q, -bootstrap_direction_consistency, genus_label)][1:8]
  fail_if(nrow(genus) != 8L || any(genus$n_pairs != 21L),
          "Figure 4 top-8 paired CLR genus 行不完整。")
  genus[, direction_class := fifelse(
    median_clr_difference_tumor_minus_normal >= 0,
    "higher in tumour", "lower in tumour"
  )]
  genus[, display_label := paste0(
    short_label(gsub("^g__", "", genus_label), 24L),
    " · q=", format_q(paired_wilcoxon_q)
  )]
  genus[, display_label := factor(
    display_label,
    levels = rev(genus[order(
      median_clr_difference_tumor_minus_normal
    ), display_label])
  )]
  p_d <- ggplot(
    genus,
    aes(
      x = median_clr_difference_tumor_minus_normal,
      y = display_label, colour = direction_class, shape = direction_class
    )
  ) +
    geom_vline(xintercept = 0, colour = "#AEB3B8", linewidth = 0.30) +
    geom_errorbar(
      aes(xmin = clr_effect_ci95_low, xmax = clr_effect_ci95_high),
      orientation = "y", width = 0.15, linewidth = 0.38
    ) +
    geom_point(size = 2.05) +
    scale_colour_manual(values = c(
      `higher in tumour` = palette_contract[["signal_blue"]],
      `lower in tumour` = palette_contract[["accent_red"]]
    )) +
    scale_shape_manual(values = c(
      `higher in tumour` = 16, `lower in tumour` = 17
    )) +
    labs(
      title = "Stable paired genus-level CLR shifts",
      subtitle = "Top eight formal FDR rows; source labels retained",
      x = "Median CLR difference (tumour − paired non-tumour)", y = NULL,
      colour = "Direction", shape = "Direction"
    ) +
    theme_nature_contract() +
    theme(
      legend.position = "bottom", axis.text.y = element_text(size = 5.0)
    ) +
    guides(colour = "none", shape = guide_legend(nrow = 1))

  metabolism_column <- p_a / p_b + plot_layout(heights = c(0.62, 1.38))
  microbiome_column <- p_c / p_d + plot_layout(heights = c(0.88, 1.12))
  (metabolism_column | plot_spacer() | microbiome_column) +
    plot_layout(widths = c(1, 0.055, 1.06)) +
    plot_annotation(
      title = "Metabolomic and microbial disease-context signals",
      subtitle = paste(
        "PR001876: 83 robust LC-MS features;",
        "PRJNA766558: paired diversity and genus-level shifts"
      ),
      tag_levels = "a"
    ) &
    theme(
      plot.title = element_text(size = 7.3, face = "bold", family = font_family),
      plot.subtitle = element_text(size = 5.7, colour = "#555555",
                                   family = font_family),
      plot.tag = element_text(size = 8, face = "bold", family = font_family)
    )
}

build_figure5_positive <- function(x) {
  retained <- copy(x$leakage_candidate_summary)[countable_for_T3_T4 == TRUE]
  priority <- copy(x$integrated_drivers)[
    decision == "strong_patient_level_candidate" &
      integrated_tier %chin% c("T3", "T4")
  ]
  fail_if(nrow(retained) != 9L || uniqueN(retained$gene_name) != 7L,
          "Figure 5 中心网络必须由 9 条保留边和 7 个基因构成。")
  fail_if(nrow(priority) != 7L ||
            nrow(priority[integrated_tier == "T4"]) != 1L ||
            priority[integrated_tier == "T4", gene_name] != "PIK3CA",
          "Figure 5 T3/T4 优先级必须为 7 个基因且 PIK3CA 唯一 T4。")

  selected_pathways <- data.table(
    factor = c(
      rep("Factor1", 3L), rep("Factor3", 4L)
    ),
    pathway = c(
      "Estrogen", "NFkB", "TNFa",
      "VEGF", "EGFR", "TGFb", "MAPK"
    )
  )
  pathway <- merge(
    x$ecms_adjusted[is_primary_scope == TRUE],
    selected_pathways,
    by = c("factor", "pathway"), all = FALSE, sort = FALSE
  )
  formal_t2 <- unique(x$integrated_axis_edges[
    edge_class == "TCGA_ECMS_adjusted_factor_PROGENy" &
      evidence_tier == "T2",
    .(factor = source_node, pathway = target_node)
  ])
  pathway <- merge(
    pathway, formal_t2,
    by = c("factor", "pathway"), all = FALSE, sort = FALSE
  )
  fail_if(nrow(pathway) != 7L,
          "Figure 5 预选 Factor1/Factor3 adjusted pathway T2 行不完整。")

  nodes <- data.table(
    node_id = c(
      paste0("GENE:", c(
        "PIK3CA", "NFE2L2", "SMARCA4", "TP53",
        "PTCH1", "NOTCH1", "KMT2D"
      )),
      "RNA:PIK3CA",
      "STATE:Factor1", "STATE:Factor7", "STATE:Factor3",
      "ECMS:official78",
      paste0("PATH:", selected_pathways$pathway)
    ),
    label = c(
      "PIK3CA", "NFE2L2", "SMARCA4", "TP53",
      "PTCH1", "NOTCH1", "KMT2D",
      "PIK3CA RNA",
      "Factor1", "Factor7", "Factor3",
      "ECMS anchor\n(official 78)",
      selected_pathways$pathway
    ),
    node_type = c(
      rep("driver", 7L), "RNA",
      rep("MOFA factor", 3L), "ECMS context",
      rep("adjusted pathway", 7L)
    ),
    x = c(
      rep(0.06, 7L), 0.28,
      rep(0.46, 3L), 0.67,
      rep(0.91, 7L)
    ),
    y = c(
      0.86, 0.70, 0.58, 0.47, 0.31, 0.22, 0.10,
      0.95,
      0.70, 0.30, 0.10,
      0.52,
      0.92, 0.84, 0.76, 0.36, 0.29, 0.22, 0.15
    )
  )
  nodes[, tier := fcase(
    node_id == "GENE:PIK3CA", "T4",
    node_type == "driver", "T3",
    node_type == "MOFA factor", "T3",
    default = "T2"
  )]

  add_edge_coordinates <- function(edges) {
    output <- merge(
      edges,
      nodes[, .(source = node_id, x_source = x, y_source = y)],
      by = "source", all.x = TRUE, sort = FALSE
    )
    output <- merge(
      output,
      nodes[, .(target = node_id, x_target = x, y_target = y)],
      by = "target", all.x = TRUE, sort = FALSE
    )
    fail_if(anyNA(output[, .(x_source, y_source, x_target, y_target)]),
            "Figure 5 中心网络存在无法定位的节点。")
    output
  }

  event_edges <- retained[, .(
    source = paste0("GENE:", gene_name),
    target = paste0("STATE:", reference_factor),
    edge_class = "retained event–state",
    event_type,
    effect = original_effect,
    q_value = original_q_value,
    direction = fifelse(original_effect >= 0, "positive", "negative")
  )]
  event_edges <- add_edge_coordinates(event_edges)

  pik_driver <- priority[gene_name == "PIK3CA"]
  dosage_edge <- add_edge_coordinates(data.table(
    source = "GENE:PIK3CA",
    target = "RNA:PIK3CA",
    edge_class = "CNV–RNA dosage",
    event_type = "relative_cnv",
    effect = pik_driver$spearman_rho_relative_cnv,
    q_value = pik_driver$q_relative_cnv,
    direction = "positive"
  ))

  context_source <- copy(x$ecms_factor)[
    is_primary_scope == TRUE & factor %chin% c("Factor1", "Factor3")
  ]
  context_edges <- add_edge_coordinates(context_source[, .(
    source = paste0("STATE:", factor),
    target = "ECMS:official78",
    edge_class = "ECMS context",
    event_type = NA_character_,
    effect = eta_squared,
    q_value = anova_q_value,
    direction = "context"
  )])
  pathway_edges <- add_edge_coordinates(pathway[, .(
    source = paste0("STATE:", factor),
    target = paste0("PATH:", pathway),
    edge_class = "ECMS-adjusted pathway",
    event_type = NA_character_,
    effect = factor_beta_standardized,
    q_value = incremental_q_value,
    direction = fifelse(factor_beta_standardized >= 0, "positive", "negative"),
    partial_r_squared
  )])

  edge_palette <- c(
    positive = palette_contract[["signal_blue"]],
    negative = palette_contract[["accent_red"]],
    context = palette_contract[["signal_teal"]]
  )
  edge_linetype <- c(positive = "solid", negative = "dashed", context = "dotted")
  node_palette <- c(
    driver = tier_palette[["T3"]],
    RNA = "#BFD8D3",
    `MOFA factor` = "#AFC8E6",
    `ECMS context` = "#D3E5E2",
    `adjusted pathway` = "#E7EEF6"
  )
  node_shapes <- c(
    driver = 21, RNA = 21, `MOFA factor` = 22,
    `ECMS context` = 23, `adjusted pathway` = 24
  )

  p_a <- ggplot() +
    geom_segment(
      data = event_edges[!grepl("^GENE:PIK3CA$", source)],
      aes(
        x = x_source, y = y_source, xend = x_target, yend = y_target,
        colour = direction, linetype = direction
      ),
      linewidth = 0.78, alpha = 0.90
    ) +
    geom_curve(
      data = event_edges[
        source == "GENE:PIK3CA" & event_type == "amplification"
      ],
      aes(
        x = x_source, y = y_source, xend = x_target, yend = y_target,
        colour = direction, linetype = direction
      ),
      curvature = 0.17, linewidth = 1.12
    ) +
    geom_curve(
      data = event_edges[
        source == "GENE:PIK3CA" & event_type == "relative_cnv"
      ],
      aes(
        x = x_source, y = y_source, xend = x_target, yend = y_target,
        colour = direction, linetype = direction
      ),
      curvature = -0.17, linewidth = 1.12
    ) +
    geom_segment(
      data = dosage_edge,
      aes(
        x = x_source, y = y_source, xend = x_target, yend = y_target,
        colour = direction, linetype = direction
      ),
      linewidth = 0.92
    ) +
    geom_segment(
      data = context_edges,
      aes(
        x = x_source, y = y_source, xend = x_target, yend = y_target,
        colour = direction, linetype = direction
      ),
      linewidth = 0.58
    ) +
    geom_segment(
      data = pathway_edges,
      aes(
        x = x_source, y = y_source, xend = x_target, yend = y_target,
        colour = direction, linetype = direction,
        linewidth = partial_r_squared
      ),
      alpha = 0.78
    ) +
    geom_point(
      data = nodes,
      aes(x = x, y = y, fill = node_type, shape = node_type),
      size = 2.75, colour = "#303030", stroke = 0.40
    ) +
    geom_point(
      data = nodes[node_id == "GENE:PIK3CA"],
      aes(x = x, y = y), shape = 21, size = 4.15,
      fill = NA, colour = "#151515", stroke = 0.65
    ) +
    geom_text(
      data = nodes,
      aes(x = x, y = y, label = label),
      family = font_family, size = 1.83, nudge_y = 0.032,
      lineheight = 0.90, check_overlap = FALSE
    ) +
    annotate(
      "text", x = 0.245, y = 0.855,
      label = paste0(
        "relative CNV–F1: ", sprintf("%.3f", retained[
          gene_name == "PIK3CA" & event_type == "relative_cnv",
          original_effect
        ]), "; q=", format_q(retained[
          gene_name == "PIK3CA" & event_type == "relative_cnv",
          original_q_value
        ])
      ),
      family = font_family, size = 1.58, colour = "#274F78"
    ) +
    annotate(
      "text", x = 0.255, y = 0.785,
      label = paste0(
        "amplification–F1: ", sprintf("%.3f", retained[
          gene_name == "PIK3CA" & event_type == "amplification",
          original_effect
        ]), "; q=", format_q(retained[
          gene_name == "PIK3CA" & event_type == "amplification",
          original_q_value
        ])
      ),
      family = font_family, size = 1.58, colour = "#274F78"
    ) +
    annotate(
      "text", x = 0.135, y = 0.985,
      label = paste0(
        "CNV–RNA ρ=", sprintf("%.3f", pik_driver$spearman_rho_relative_cnv),
        "; q=", format_q(pik_driver$q_relative_cnv)
      ),
      family = font_family, size = 1.58, colour = "#376C68"
    ) +
    annotate(
      "text", x = 0.57, y = 0.65,
      label = paste0(
        "Factor1 ECMS η²=", sprintf("%.3f", context_source[
          factor == "Factor1", eta_squared
        ])
      ),
      family = font_family, size = 1.52, colour = "#376C68"
    ) +
    annotate(
      "text", x = 0.57, y = 0.36,
      label = paste0(
        "Factor3 ECMS η²=", sprintf("%.3f", context_source[
          factor == "Factor3", eta_squared
        ])
      ),
      family = font_family, size = 1.52, colour = "#376C68"
    ) +
    scale_colour_manual(values = edge_palette) +
    scale_linetype_manual(values = edge_linetype) +
    scale_linewidth_continuous(range = c(0.28, 0.88), guide = "none") +
    scale_fill_manual(values = node_palette) +
    scale_shape_manual(values = node_shapes) +
    coord_cartesian(xlim = c(0, 0.98), ylim = c(0.04, 1.02), clip = "off") +
    labs(
      title = "PIK3CA–Factor1-centred source-internal network",
      subtitle = paste(
        "Nine retained event–state edges connect dosage, ECMS context",
        "and selected adjusted pathways"
      ),
      colour = "Association", linetype = "Association",
      fill = "Node", shape = "Node"
    ) +
    theme_void(base_family = font_family, base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 0.5, face = "bold"),
      plot.subtitle = element_text(size = base_size - 0.4, colour = "#555555"),
      legend.position = "bottom", legend.box = "vertical",
      legend.title = element_text(size = base_size - 0.2, face = "bold"),
      legend.text = element_text(size = base_size - 0.7),
      plot.margin = margin(3.5, 4.5, 3.5, 3.5, unit = "mm")
    ) +
    guides(
      colour = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1),
      fill = "none", shape = "none"
    )

  priority[, tier_order := match(integrated_tier, c("T3", "T4"))]
  setorder(priority, -tier_order, -state_supported_edges, gene_name)
  priority[, gene_factor := factor(gene_name, levels = rev(gene_name))]
  priority[, tier_label := fifelse(
    integrated_tier == "T4", "T4 priority", "T3 retained"
  )]
  p_b <- ggplot(
    priority,
    aes(x = state_supported_edges, y = gene_factor, colour = integrated_tier)
  ) +
    geom_segment(
      aes(x = 0, xend = state_supported_edges, yend = gene_factor),
      linewidth = 0.68
    ) +
    geom_point(size = 2.45) +
    geom_point(
      data = priority[integrated_tier == "T4"],
      shape = 21, size = 3.65, fill = NA, colour = "#151515", stroke = 0.65
    ) +
    geom_text(
      aes(label = tier_label),
      hjust = -0.18, family = font_family, size = 1.78
    ) +
    scale_colour_manual(values = c(
      T3 = tier_palette[["T3"]], T4 = tier_palette[["T4"]]
    )) +
    scale_x_continuous(breaks = 0:2) +
    coord_cartesian(xlim = c(0, 2.95), clip = "off") +
    labs(
      title = "Validation-priority ladder",
      subtitle = "PIK3CA leads the seven-gene T3/T4 validation set",
      x = "Retained event–state edges", y = NULL, colour = "Tier"
    ) +
    theme_nature_contract() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 5.6, face = "bold")
    )

  (p_a | p_b) +
    plot_layout(widths = c(1.72, 0.63)) +
    plot_annotation(
      title = "PIK3CA–Factor1 anchors a compact validation-priority network",
      subtitle = paste(
        "Seven T3/T4 genes connect retained event–state associations",
        "with ECMS-aligned pathways"
      ),
      tag_levels = "a"
    ) &
    theme(
      plot.title = element_text(size = 7.3, face = "bold", family = font_family),
      plot.subtitle = element_text(size = 5.7, colour = "#555555",
                                   family = font_family),
      plot.tag = element_text(size = 8, face = "bold", family = font_family)
    )
}

run_figure_export_bundle <- function(inputs) {
message("[1/4] 构建 5 张 R-only 主图对象（不重算上游统计）")
figure_objects <- list(
  Figure1 = build_figure1_positive(inputs),
  Figure2 = build_figure2_positive(inputs),
  Figure3 = build_figure3_positive(inputs),
  Figure4 = build_figure4_positive(inputs),
  Figure5 = build_figure5_positive(inputs)
)
fail_if(!identical(names(figure_objects), figure_specs$figure_id),
        "图对象与冻结 figure_specs 顺序不一致。")
fail_if(any(!vapply(
  figure_objects,
  function(object) inherits(object, c("ggplot", "patchwork")),
  logical(1)
)), "至少一个主图对象不是 R ggplot/patchwork 对象。")

source_manifest_rows <- rbindlist(lapply(names(verified_families), function(name) {
  data.table(
    family = name,
    manifest = basename(verified_families[[name]]$manifest_path),
    manifest_sha256 = verified_families[[name]]$manifest_sha256
  )
}))
setorder(source_manifest_rows, family, manifest)
source_manifest_bundle <- paste0(
  source_manifest_rows$manifest, "=", source_manifest_rows$manifest_sha256,
  collapse = ";"
)
source_table_rows <- data.table(
  input_name = names(input_paths),
  relative_path = file.path("results", basename(unlist(input_paths))),
  sha256 = vapply(unlist(input_paths), sha256_file, character(1))
)
setorder(source_table_rows, input_name, relative_path)
source_table_bundle <- paste0(
  source_table_rows$relative_path, "=", source_table_rows$sha256,
  collapse = ";"
)
source_bundle_signature <- digest(
  list(
    source_manifests = source_manifest_rows,
    source_tables = source_table_rows,
    figure_contract = figure_specs,
    frozen_figure_contract_sha256 = figure_contract_sha256
  ),
  algo = "sha256", serialize = TRUE
)

stage_dir <- file.path(
  work_intermediate_dir,
  paste0(".escc_multiomics_figures_", substr(source_bundle_signature, 1L, 12L))
)
fail_if(dir_exists(stage_dir), paste(
  "图件 stage_dir 已存在；可能是未审查的中断运行：", stage_dir,
  "。请先人工核查，禁止自动覆盖。"
))
dir_create(stage_dir, recurse = TRUE)
stage_active <- TRUE
on.exit({
  if (stage_active && dir_exists(stage_dir)) dir_delete(stage_dir)
}, add = TRUE)

render_plot <- function(plot, path, format, width_mm, height_mm, dpi = 600L) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  device_open <- FALSE
  on.exit({
    if (device_open) try(grDevices::dev.off(), silent = TRUE)
  }, add = TRUE)
  if (format == "svg") {
    svglite::svglite(
      path, width = width_in, height = height_in,
      bg = "white", system_fonts = list(sans = font_family)
    )
  } else if (format == "pdf") {
    svg_source <- sub("[.]pdf$", ".svg", path)
    fail_if(
      !file_exists(svg_source),
      paste("PDF 转换缺少已验证的同图 SVG：", svg_source)
    )
    rsvg::rsvg_pdf(svg_source, path)
    return(invisible(path))
  } else if (format == "tiff") {
    ragg::agg_tiff(
      path, width = width_in, height = height_in, units = "in",
      res = dpi, background = "white", compression = "lzw",
      scaling = 1
    )
  } else if (format == "png") {
    ragg::agg_png(
      path, width = width_in, height = height_in, units = "in",
      res = dpi, background = "white", scaling = 1
    )
  } else {
    stop("未支持的图件格式：", format, call. = FALSE)
  }
  device_open <- TRUE
  print(plot)
  grDevices::dev.off()
  device_open <- FALSE
  invisible(path)
}

validate_rendered_file <- function(path, format) {
  fail_if(!file_exists(path), paste("渲染文件缺失：", path))
  size <- as.numeric(file_info(path)$size)
  fail_if(!is.finite(size) || size < 1024,
          paste("渲染文件过小或为空：", path))
  raw_prefix <- readBin(path, what = "raw", n = min(size, 4096L))
  if (format == "pdf") {
    fail_if(rawToChar(raw_prefix[seq_len(min(5L, length(raw_prefix)))]) != "%PDF-",
            paste("PDF magic header 异常：", path))
    pdf_info <- pdftools::pdf_info(path)
    pdf_fonts <- pdftools::pdf_fonts(path)
    fail_if(pdf_info$pages != 1L || !nrow(pdf_fonts) ||
              any(!pdf_fonts$embedded),
            paste("PDF 页数或字体嵌入检查失败：", path))
  } else if (format == "png") {
    png_magic <- as.raw(c(0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A))
    fail_if(length(raw_prefix) < 8L || !identical(raw_prefix[1:8], png_magic),
            paste("PNG magic header 异常：", path))
  } else if (format == "tiff") {
    tiff_le <- as.raw(c(0x49, 0x49, 0x2A, 0x00))
    tiff_be <- as.raw(c(0x4D, 0x4D, 0x00, 0x2A))
    fail_if(length(raw_prefix) < 4L ||
              (!identical(raw_prefix[1:4], tiff_le) &&
                 !identical(raw_prefix[1:4], tiff_be)),
            paste("TIFF magic header 异常：", path))
  } else if (format == "svg") {
    svg_text <- rawToChar(raw_prefix)
    fail_if(!grepl("<svg", svg_text, fixed = TRUE),
            paste("SVG 根元素缺失：", path))
    fail_if(grepl("<image[^>]+data:image", svg_text),
            paste("SVG 疑似整体栅格化，未保留可编辑矢量文本/图元：", path))
  }
  invisible(TRUE)
}

message("[2/4] 以 R 导出 SVG/PDF/600 dpi TIFF/PNG 并执行结构 QA")
stage_records <- vector("list", nrow(figure_specs) * length(figure_formats))
record_index <- 0L
for (i in seq_len(nrow(figure_specs))) {
  spec <- figure_specs[i, ]
  plot <- figure_objects[[spec$figure_id]]
  for (format in figure_formats) {
    record_index <- record_index + 1L
    filename <- paste0(spec$file_stem, ".", format)
    stage_path <- file.path(stage_dir, filename)
    render_plot(
      plot, stage_path, format,
      width_mm = spec$width_mm,
      height_mm = spec$height_mm,
      dpi = 600L
    )
    validate_rendered_file(stage_path, format)
    stage_records[[record_index]] <- data.table(
      figure_id = spec$figure_id,
      artifact = filename,
      relative_path = file.path("figures", filename),
      format = format,
      width_mm = spec$width_mm,
      height_mm = spec$height_mm,
      dpi = if (format %chin% c("tiff", "png")) 600L else NA_integer_,
      file_size_bytes = as.numeric(file_info(stage_path)$size),
      sha256 = sha256_file(stage_path),
      claim = spec$claim,
      archetype = spec$archetype,
      panel_source_contract = spec$panel_source_contract,
      filter_contract = spec$filter_contract,
      statistical_unit = spec$statistical_unit,
      sample_structure = spec$sample_structure,
      effect_definition = spec$effect_definition,
      p_q_family = spec$p_q_family,
      frozen_figure_contract_sha256 = figure_contract_sha256,
      statistical_encoding = paste(
        "formal source tables only; no manual p/q/stars;",
        "continuous effects separated from categorical evidence tiers"
      ),
      reviewer_risk = spec$reviewer_risk,
      backend = "R-only",
      font_family = font_family,
      minimum_text_pt = minimum_text_pt,
      generation_script = "scripts/24_visualize_escc_multiomics.R",
      generation_script_sha256 = executed_script_sha256,
      source_manifest_sha256 = source_manifest_bundle,
      source_table_sha256 = source_table_bundle,
      source_bundle_signature = source_bundle_signature,
      generated_date = as.character(Sys.Date()),
      structural_status = "verified_after_staged_export_and_hash_check",
      visual_qa_status = "pending_reopened_review",
      qa_path = NA_character_,
      qa_sha256 = NA_character_
    )
  }
}
artifact_manifest <- rbindlist(stage_records)
fail_if(nrow(artifact_manifest) != 20L ||
          uniqueN(artifact_manifest$relative_path) != 20L,
        "figure artifact manifest 必须恰好登记 20 个唯一图文件。")
fail_if(!setequal(artifact_manifest$relative_path, formal_figure_relative_paths),
        "figure artifact manifest 与 --fields-only 输出集合不一致。")
manifest_stage_path <- file.path(
  stage_dir, basename(figure_manifest_relative_path)
)
fwrite(
  artifact_manifest, manifest_stage_path,
  sep = "\t", quote = FALSE, na = "", logical01 = FALSE
)
fail_if(!file_exists(manifest_stage_path) ||
          as.numeric(file_info(manifest_stage_path)$size) < 1024,
        "figure artifact manifest 暂存失败。")

qa_lines <- c(
  "# ESCC 多组学主图生成结构 QA",
  "",
  "> 本文件属于 `_work/checks/` 过程性历史证据，不是项目当前状态源。",
  "",
  paste0("- 运行日期：", Sys.Date(), "。"),
  paste0("- 执行脚本 SHA256：", executed_script_sha256, "。"),
  paste0("- 上游 bundle signature：", source_bundle_signature, "。"),
  paste0("- 冻结图件契约 SHA256：", figure_contract_sha256, "。"),
  "- 后端：仅 R；未调用 Python 或其他后端生成、预览或导出图件。",
  "- 输出：5 张主图 × SVG/PDF/TIFF/PNG；TIFF 与 PNG 均为 600 dpi。",
  paste0("- 最小目标字号：", minimum_text_pt, " pt。"),
  "- SVG 根元素、PDF/PNG/TIFF magic header、文件非空和 stage SHA256 回读：通过。",
  "- 发布事务：20 个图件、过程 QA 和 manifest 使用同盘备份与整包回滚；manifest literal-last。",
  "- Figure 1：正文主视觉恰含 9 条 countable_for_T3_T4 边；完整 41 边及执行/对齐诊断留在正式源表。",
  "- Figure 2：仅 official 78 的 Factor1/Factor3；extension 16 与 k 稳定性诊断未进入正文主视觉。",
  "- Factor4：正式整合表继续冻结为 T0/HM450 level-factor；正文主视觉不展示。",
  "- Figure 3：GNAS/ZNF750 的 promoter、RNA、protein 独立尺度；缺失未填 0；无跨层连接线。",
  "- Figure 4：PR001876 层级计数为 83/43/5/4；PR001876 与 PRJNA766558 无共享患者连接。",
  "- Figure 5：中心网络含 9 条 retained event–factor 边与 7 个 T3/T4 基因；完整 12 行裁决留在 results/；所有关联线无箭头。",
  "- structural_status 已通过；visual_qa_status=pending_reopened_review。",
  "- 视觉 QA 尚需在最终尺寸重新打开 PNG/TIFF/SVG/PDF，由人工检查标签重叠、裁切和可读性；当前不宣称视觉完成。",
  "- 原始数据、ResearchDataHub、正式上游结果和投稿包未修改。"
)
qa_path <- file.path(
  work_checks_dir,
  paste0("escc_multiomics_figure_generation_qa_", format(Sys.Date(), "%Y%m%d"), ".md")
)
qa_stage_path <- file.path(stage_dir, basename(qa_path))
writeLines(qa_lines, qa_stage_path, useBytes = TRUE)
fail_if(!file_exists(qa_stage_path) ||
          as.numeric(file_info(qa_stage_path)$size) < 1024,
        "图件结构 QA 暂存失败。")

copy_verified <- function(source, destination) {
  fail_if(!file_exists(source), paste("复制源文件缺失：", source))
  fail_if(file_exists(destination), paste("复制目标已存在：", destination))
  dir_create(dirname(destination), recurse = TRUE)
  copied <- file.copy(
    source, destination, overwrite = FALSE, copy.mode = TRUE,
    copy.date = TRUE
  )
  fail_if(!copied || !file_exists(destination),
          paste("验证复制失败：", destination))
  fail_if(as.numeric(file_info(source)$size) !=
            as.numeric(file_info(destination)$size) ||
            sha256_file(source) != sha256_file(destination),
          paste("复制后大小/SHA256 不一致：", destination))
  invisible(destination)
}

publish_figure_bundle_transaction <- function(
    artifact_manifest,
    stage_dir,
    qa_stage_path,
    qa_destination,
    manifest_stage_path,
    manifest_destination) {
  figure_sources <- file.path(stage_dir, artifact_manifest$artifact)
  figure_destinations <- file.path(project_root, artifact_manifest$relative_path)
  sources <- c(figure_sources, qa_stage_path, manifest_stage_path)
  destinations <- c(
    figure_destinations, qa_destination, manifest_destination
  )
  relative_labels <- c(
    artifact_manifest$relative_path,
    file.path("_work", "checks", basename(qa_destination)),
    figure_manifest_relative_path
  )
  artifact_roles <- c(
    rep("figure", nrow(artifact_manifest)), "process_qa", "manifest"
  )
  fail_if(length(sources) != length(destinations) ||
            anyDuplicated(destinations) > 0L ||
            tail(artifact_roles, 1L) != "manifest" ||
            tail(destinations, 1L) != manifest_destination,
          "bundle transaction 路径、角色或 manifest-last 顺序异常。")
  fail_if(any(!file_exists(sources)),
          "bundle transaction 的 stage 输入不完整。")

  destination_dirs <- unique(dirname(destinations))
  fail_if(any(!dir_exists(destination_dirs)),
          "bundle transaction 目标目录缺失。")
  device_paths <- c(work_intermediate_dir, destination_dirs)
  device_ids <- as.character(file_info(device_paths)$device_id)
  fail_if(anyNA(device_ids) || uniqueN(device_ids) != 1L,
          paste(
            "stage/backup 与正式目标不在同一文件系统；禁止执行非原子 bundle",
            "transaction。"
          ))

  lock_dir <- file.path(work_intermediate_dir, ".escc_figure_publish.lock")
  lock_acquired <- FALSE
  lock_release_allowed <- FALSE
  # 锁获取与 handler 注册同处不可中断区，不留“有锁无清理”窗口。
  suspendInterrupts({
    lock_acquired <- dir.create(
      lock_dir, recursive = FALSE, showWarnings = FALSE
    )
    fail_if(!lock_acquired,
            paste("figure bundle 发布锁已存在，禁止并发发布：", lock_dir))
    lock_release_allowed <- TRUE
    on.exit(suspendInterrupts({
      if (lock_acquired && lock_release_allowed && dir_exists(lock_dir)) {
        try(dir_delete(lock_dir), silent = TRUE)
      }
    }), add = TRUE)
  })

  transaction_id <- paste0(
    Sys.getpid(), "-", format(Sys.time(), "%Y%m%d%H%M%OS6"), "-",
    substr(source_bundle_signature, 1L, 10L)
  )
  transaction_id <- gsub("[^A-Za-z0-9.-]", "", transaction_id)
  transaction_root <- file.path(
    work_intermediate_dir,
    paste0(".escc_figure_publish_bundle_", transaction_id)
  )
  hidden_publish <- character()
  hidden_rollback <- character()
  transaction_root_owned <- FALSE
  formal_targets_touched <- FALSE
  transaction_committed <- FALSE
  transaction_cleanup_verified <- FALSE
  rollback_verified <- FALSE
  rollback_cleanup_verified <- FALSE
  # transaction_root 创建与 handler 注册也不可被中断分开。
  suspendInterrupts({
    fail_if(dir_exists(transaction_root),
            paste("bundle transaction 目录已存在：", transaction_root))
    dir_create(file.path(transaction_root, "backup"), recurse = TRUE)
    transaction_root_owned <- TRUE
    on.exit(suspendInterrupts({
      for (path in c(hidden_publish, hidden_rollback)) {
        if (file_exists(path)) try(file_delete(path), silent = TRUE)
      }
      if (transaction_root_owned &&
          (!formal_targets_touched || transaction_cleanup_verified ||
             rollback_cleanup_verified) && dir_exists(transaction_root)) {
        try(dir_delete(transaction_root), silent = TRUE)
      }
    }), add = TRUE)
  })

  hidden_publish <- file.path(
    dirname(destinations),
    paste0(".", basename(destinations), ".publish-", transaction_id)
  )
  hidden_rollback <- file.path(
    dirname(destinations),
    paste0(".", basename(destinations), ".rollback-", transaction_id)
  )
  fail_if(any(file_exists(c(hidden_publish, hidden_rollback))),
          "当前 bundle transaction 的隐藏 publish/rollback 路径已存在。")

  original_exists <- file_exists(destinations)
  original_size <- rep(NA_real_, length(destinations))
  original_sha <- rep(NA_character_, length(destinations))
  original_size[original_exists] <- as.numeric(
    file_info(destinations[original_exists])$size
  )
  original_sha[original_exists] <- vapply(
    destinations[original_exists], sha256_file, character(1)
  )
  backup_paths <- file.path(transaction_root, "backup", relative_labels)
  published <- rep(FALSE, length(destinations))

  backup_error <- NULL
  backup_ok <- tryCatch({
    for (i in which(original_exists)) {
      copy_verified(destinations[[i]], backup_paths[[i]])
    }
    TRUE
  }, error = function(condition) {
    backup_error <<- condition
    FALSE
  }, interrupt = function(condition) {
    backup_error <<- condition
    FALSE
  })
  if (!backup_ok) {
    backup_cleanup_errors <- character()
    suspendInterrupts({
      rollback_verified <- TRUE  # 尚未替换任何正式目标。
      for (path in c(hidden_publish, hidden_rollback)) {
        tryCatch({
          if (file_exists(path)) file_delete(path)
          fail_if(file_exists(path), paste("中断后隐藏路径未清理：", path))
        }, error = function(condition) {
          backup_cleanup_errors <<- c(
            backup_cleanup_errors, conditionMessage(condition)
          )
        })
      }
      tryCatch({
        if (dir_exists(transaction_root)) dir_delete(transaction_root)
        fail_if(dir_exists(transaction_root),
                paste("备份中断后 transaction_root 未清理：",
                      transaction_root))
      }, error = function(condition) {
        backup_cleanup_errors <<- c(
          backup_cleanup_errors, conditionMessage(condition)
        )
      })
      rollback_cleanup_verified <- !length(backup_cleanup_errors)
      lock_release_allowed <- rollback_cleanup_verified
      if (rollback_cleanup_verified) {
        tryCatch({
          if (dir_exists(lock_dir)) dir_delete(lock_dir)
          fail_if(dir_exists(lock_dir), paste("备份中断后发布锁未清理：",
                                              lock_dir))
        }, error = function(condition) {
          backup_cleanup_errors <<- c(
            backup_cleanup_errors, conditionMessage(condition)
          )
        })
        rollback_cleanup_verified <- !length(backup_cleanup_errors)
        lock_release_allowed <- rollback_cleanup_verified
        if (rollback_cleanup_verified) lock_acquired <- FALSE
      }
    })
    if (!rollback_cleanup_verified) {
      stop(
        paste0(
          "bundle transaction 在备份阶段停止；正式目标未修改，",
          "但清理未完整：", conditionMessage(backup_error),
          "；cleanup errors=", paste(backup_cleanup_errors, collapse = " | "),
          "；发布锁保留待人工检查：", lock_dir
        ),
        call. = FALSE
      )
    }
    stop(
      paste0(
        "bundle transaction 在备份阶段停止，正式目标未修改且残留已清理：",
        conditionMessage(backup_error)
      ),
      call. = FALSE
    )
  }

  perform_bundle_rollback <- function() {
    rollback_errors_local <- character()
    rollback_cleanup_errors_local <- character()
    rollback_verified_local <- FALSE
    rollback_cleanup_verified_local <- FALSE
    suspendInterrupts({
      for (i in rev(seq_along(destinations))) {
        tryCatch({
          if (file_exists(hidden_publish[[i]])) file_delete(hidden_publish[[i]])
          if (original_exists[[i]]) {
            fail_if(!file_exists(backup_paths[[i]]),
                    paste("回滚备份缺失：", backup_paths[[i]]))
            if (file_exists(hidden_rollback[[i]])) file_delete(hidden_rollback[[i]])
            copy_verified(backup_paths[[i]], hidden_rollback[[i]])
            restored <- file.rename(hidden_rollback[[i]], destinations[[i]])
            if (!restored) {
              # 同盘 rename 极少数情况下不能覆盖；删除本事务新目标后重试。
              if (file_exists(destinations[[i]])) file_delete(destinations[[i]])
              restored <- file.rename(hidden_rollback[[i]], destinations[[i]])
            }
            if (!restored) {
              copied_back <- file.copy(
                backup_paths[[i]], destinations[[i]], overwrite = TRUE,
                copy.mode = TRUE, copy.date = TRUE
              )
              fail_if(!copied_back, paste("回滚复制失败：", destinations[[i]]))
            }
            fail_if(!file_exists(destinations[[i]]) ||
                      as.numeric(file_info(destinations[[i]])$size) !=
                        original_size[[i]] ||
                      sha256_file(destinations[[i]]) != original_sha[[i]],
                    paste("回滚 SHA256 验证失败：", destinations[[i]]))
          } else if (file_exists(destinations[[i]])) {
            file_delete(destinations[[i]])
            fail_if(file_exists(destinations[[i]]),
                    paste("回滚未能删除本事务新增目标：", destinations[[i]]))
          }
        }, error = function(condition) {
          rollback_errors_local <<- c(
            rollback_errors_local,
            paste0(destinations[[i]], " => ", conditionMessage(condition))
          )
        })
      }
      rollback_verified_local <- !length(rollback_errors_local)

      if (rollback_verified_local) {
        for (path in c(hidden_publish, hidden_rollback)) {
          tryCatch({
            if (file_exists(path)) file_delete(path)
            fail_if(file_exists(path), paste("回滚后隐藏路径未清理：", path))
          }, error = function(condition) {
            rollback_cleanup_errors_local <<- c(
              rollback_cleanup_errors_local, conditionMessage(condition)
            )
          })
        }
        tryCatch({
          if (dir_exists(transaction_root)) dir_delete(transaction_root)
          fail_if(dir_exists(transaction_root),
                  paste("回滚后 transaction_root 未清理：", transaction_root))
        }, error = function(condition) {
          rollback_cleanup_errors_local <<- c(
            rollback_cleanup_errors_local, conditionMessage(condition)
          )
        })
        rollback_cleanup_verified_local <-
          !length(rollback_cleanup_errors_local)
      }

      lock_release_allowed <<-
        rollback_verified_local && rollback_cleanup_verified_local
      if (lock_release_allowed) {
        tryCatch({
          if (dir_exists(lock_dir)) dir_delete(lock_dir)
          fail_if(dir_exists(lock_dir), paste("回滚后发布锁未清理：", lock_dir))
        }, error = function(condition) {
          rollback_cleanup_errors_local <<- c(
            rollback_cleanup_errors_local, conditionMessage(condition)
          )
        })
        rollback_cleanup_verified_local <-
          !length(rollback_cleanup_errors_local)
        lock_release_allowed <<-
          rollback_verified_local && rollback_cleanup_verified_local
        if (lock_release_allowed) lock_acquired <<- FALSE
      }
      rollback_verified <<- rollback_verified_local
      rollback_cleanup_verified <<- rollback_cleanup_verified_local
    })
    list(
      rollback_verified = rollback_verified_local,
      rollback_cleanup_verified = rollback_cleanup_verified_local,
      rollback_errors = rollback_errors_local,
      rollback_cleanup_errors = rollback_cleanup_errors_local
    )
  }

  rollback_result <- NULL
  publish_error <- NULL
  post_commit_error <- NULL
  published_manifest_final <- NULL
  handle_publish_condition <- function(condition) {
    if (transaction_committed) {
      post_commit_error <<- condition
      if (transaction_cleanup_verified && !dir_exists(lock_dir)) {
        "success"
      } else {
        "committed_cleanup_incomplete"
      }
    } else {
      publish_error <<- condition
      rollback_result <<- perform_bundle_rollback()
      if (!rollback_result$rollback_verified) {
        "rollback_incomplete"
      } else if (!rollback_result$rollback_cleanup_verified) {
        "rollback_cleanup_incomplete"
      } else {
        "rolled_back"
      }
    }
  }

  # 异常 handler 会在返回 publish_state 前直接回滚；本 on.exit 只防御
  # handler 内的非预期错误，且因为后注册而会先于普通清理 handler 执行。
  on.exit(suspendInterrupts({
    if (formal_targets_touched && !transaction_committed &&
        !rollback_verified) {
      try(perform_bundle_rollback(), silent = TRUE)
    }
  }), add = TRUE, after = FALSE)

  publish_state <- tryCatch({
    suspendInterrupts({
      formal_targets_touched <- TRUE
      lock_release_allowed <- FALSE
    })
    for (i in seq_along(sources)) {
      copy_verified(sources[[i]], hidden_publish[[i]])
      moved <- file.rename(hidden_publish[[i]], destinations[[i]])
      fail_if(!moved, paste(
        "bundle 原子替换失败：", destinations[[i]],
        "；manifest 尚未或最后替换，准备整包回滚。"
      ))
      published[[i]] <- TRUE
      fail_if(as.numeric(file_info(destinations[[i]])$size) !=
                as.numeric(file_info(sources[[i]])$size) ||
                sha256_file(destinations[[i]]) != sha256_file(sources[[i]]),
              paste("bundle 发布后回读失败：", destinations[[i]]))
    }
    fail_if(!all(published), "bundle transaction 未完成全部目标替换。")

    # manifest 是循环中的最后一个目标；到此才允许按 bundle 宣称结构发布。
    published_manifest_final <- read_tsv(manifest_destination)
    require_columns(
      published_manifest_final,
      c(
        "relative_path", "sha256", "structural_status", "visual_qa_status",
        "qa_path", "qa_sha256"
      ),
      "published figure manifest"
    )
    fail_if(nrow(published_manifest_final) != 20L ||
              !setequal(published_manifest_final$relative_path,
                        artifact_manifest$relative_path) ||
              any(published_manifest_final$sha256[
                match(artifact_manifest$relative_path,
                      published_manifest_final$relative_path)
              ] != artifact_manifest$sha256) ||
              any(published_manifest_final$structural_status !=
                    "verified_after_staged_export_and_hash_check") ||
              any(published_manifest_final$visual_qa_status !=
                    "pending_reopened_review") ||
              any(!is.na(published_manifest_final$qa_path)) ||
              any(!is.na(published_manifest_final$qa_sha256)),
            "最后发布的 figure manifest 与 stage/QA 状态契约不一致。")
    # commit 与成功清理同处一个不可中断区：待处理中断只能在
    # manifest 最后发布、验证、备份清理和锁释放后被投递。
    suspendInterrupts({
      transaction_committed <- TRUE
      for (path in c(hidden_publish, hidden_rollback)) {
        if (file_exists(path)) file_delete(path)
        fail_if(file_exists(path), paste("成功提交后隐藏路径未清理：", path))
      }
      if (dir_exists(transaction_root)) dir_delete(transaction_root)
      fail_if(dir_exists(transaction_root),
              paste("成功提交后备份包未清理：", transaction_root))
      transaction_cleanup_verified <- TRUE
      lock_release_allowed <- TRUE
      if (dir_exists(lock_dir)) dir_delete(lock_dir)
      fail_if(dir_exists(lock_dir), paste("成功提交后发布锁未清理：", lock_dir))
      lock_acquired <- FALSE
    })
    "success"
  }, error = function(condition) {
    handle_publish_condition(condition)
  }, interrupt = function(condition) {
    handle_publish_condition(condition)
  })

  if (identical(publish_state, "rollback_incomplete")) {
    stop(
      paste0(
        "bundle transaction 失败且回滚不完整：",
        conditionMessage(publish_error), "；rollback errors=",
        paste(rollback_result$rollback_errors, collapse = " | "),
        "；保留恢复包：", transaction_root,
        "；发布锁保留待人工恢复：", lock_dir
      ),
      call. = FALSE
    )
  }
  if (identical(publish_state, "rollback_cleanup_incomplete")) {
    stop(
      paste0(
        "bundle transaction 失败；全部正式目标已按原大小/SHA256回滚，",
        "但隐藏路径、备份包或发布锁清理未通过：",
        conditionMessage(publish_error), "；cleanup errors=",
        paste(rollback_result$rollback_cleanup_errors, collapse = " | "),
        "；发布锁保留待人工检查：", lock_dir
      ),
      call. = FALSE
    )
  }
  if (identical(publish_state, "rolled_back")) {
    stop(
      paste0(
        "bundle transaction 失败，全部正式目标已按原大小/SHA256回滚：",
        conditionMessage(publish_error)
      ),
      call. = FALSE
    )
  }

  if (identical(publish_state, "committed_cleanup_incomplete")) {
    stop(
      paste0(
        "bundle transaction 已在 manifest-last 后提交，不再自动回滚；",
        "但提交后备份/锁清理不完整：",
        conditionMessage(post_commit_error), "；恢复路径=", transaction_root,
        "；发布锁=", lock_dir
      ),
      call. = FALSE
    )
  }
  fail_if(!identical(publish_state, "success") ||
            !transaction_committed || !transaction_cleanup_verified ||
            dir_exists(transaction_root) || dir_exists(lock_dir),
          "bundle transaction 未达到已提交且备份/锁均清理的成功终态。")
  invisible(published_manifest_final)
}

message("[3/4] 执行可回滚 bundle transaction；manifest literal-last")
published_manifest <- publish_figure_bundle_transaction(
  artifact_manifest = artifact_manifest,
  stage_dir = stage_dir,
  qa_stage_path = qa_stage_path,
  qa_destination = qa_path,
  manifest_stage_path = manifest_stage_path,
  manifest_destination = file.path(project_root, figure_manifest_relative_path)
)
fail_if(nrow(published_manifest) != 20L,
        "bundle transaction 返回的 manifest 行数异常。")

message("[4/4] 图件结构包发布完成；等待最终尺寸视觉 QA")
if (dir_exists(stage_dir)) dir_delete(stage_dir)
stage_active <- FALSE
message(
  "完成：5 张主图，20 个图件 artifact；manifest 最后发布。",
  "视觉重开 QA 尚未由本脚本冒充完成。"
)
invisible(TRUE)
}

run_figure_export_bundle(inputs)
