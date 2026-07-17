#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("缺少 digest 包；请安装后重试。", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
allowed_args <- c("--validate-only")
unknown_args <- setdiff(args, allowed_args)
if (length(unknown_args) > 0L) {
  stop("未知参数：", paste(unknown_args, collapse = ", "), call. = FALSE)
}
validate_only <- "--validate-only" %in% args

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
path <- function(...) file.path(root, ...)
read_text <- function(file) {
  paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}
sha256_file <- function(file) digest::digest(file = file, algo = "sha256")

required_dirs <- c("scripts", "results", "figures", "manuscript", "_work/intermediate")
if (!all(dir.exists(path(required_dirs)))) {
  stop("项目目录结构不完整。", call. = FALSE)
}

master_file <- path("manuscript", "escc_multiomics_manuscript_master.md")
declarations_file <- path("manuscript", "escc_multiomics_submission_declarations.md")
references_file <- path("manuscript", "references.bib")
main_tables_file <- path("manuscript", "escc_multiomics_main_tables.md")
supp_tables_file <- path("manuscript", "escc_multiomics_supplementary_tables.md")
clean_supp_index_file <- path("results", "escc_submission_supplement_index.tsv")
clean_supp_manifest_file <- path("results", "escc_submission_supplement_artifact_manifest.tsv")
clean_supp_dir <- path("results", "submission_supplement")
human_completion_file <- path("results", "escc_human_completion_items_20260717.tsv")
submission_table_files <- c(
  path("results", "escc_manuscript_table1_cohorts.tsv"),
  path("results", "escc_manuscript_table2_event_state_associations.tsv"),
  path("results", "escc_manuscript_table3_external_validation.tsv"),
  path("results", "escc_manuscript_supplementary_table_index.tsv"),
  clean_supp_index_file
)

checks <- data.table(
  check_id = character(),
  domain = character(),
  requirement = character(),
  observed = character(),
  gate_class = character(),
  status = character(),
  blocks_target_formatting = logical(),
  blocks_final_submission = logical(),
  action = character()
)

add_check <- function(
    check_id,
    domain,
    requirement,
    observed,
    gate_class = "scientific",
    status = "PASS",
    blocks_target_formatting = TRUE,
    blocks_final_submission = TRUE,
    action = "") {
  checks <<- rbind(
    checks,
    data.table(
      check_id = check_id,
      domain = domain,
      requirement = requirement,
      observed = as.character(observed),
      gate_class = gate_class,
      status = status,
      blocks_target_formatting = blocks_target_formatting,
      blocks_final_submission = blocks_final_submission,
      action = action
    ),
    use.names = TRUE
  )
}

status_if <- function(ok, warn = FALSE) {
  if (isTRUE(ok)) "PASS" else if (isTRUE(warn)) "WARN" else "FAIL"
}

required_authority_files <- c(
  master_file,
  declarations_file,
  references_file,
  main_tables_file,
  supp_tables_file
)
authority_exists <- file.exists(required_authority_files)
add_check(
  "AUTHORITY_FILES",
  "投稿权威面",
  "英文母稿、声明、参考文献、主表和补充表均存在",
  paste(basename(required_authority_files), authority_exists, sep = "=", collapse = "; "),
  status = status_if(all(authority_exists)),
  action = "缺失文件必须先由正式组装脚本生成。"
)

master <- if (file.exists(master_file)) read_text(master_file) else ""
declarations <- if (file.exists(declarations_file)) read_text(declarations_file) else ""
main_tables <- if (file.exists(main_tables_file)) read_text(main_tables_file) else ""
supp_tables <- if (file.exists(supp_tables_file)) read_text(supp_tables_file) else ""
submission_table_text <- paste(vapply(
  submission_table_files[file.exists(submission_table_files)],
  read_text,
  character(1)
), collapse = "\n")
formal_text <- paste(
  master, main_tables, supp_tables, submission_table_text, declarations,
  sep = "\n"
)

science_placeholders <- gregexpr(
  "\\{\\{(?:EXTERNAL|SCIENTIFIC|RESULT|FIGURE|TABLE)_[^}]+\\}\\}|\\[待补文献\\]",
  formal_text,
  perl = TRUE
)
science_hits <- regmatches(formal_text, science_placeholders)[[1L]]
if (identical(science_hits, character(0)) || identical(science_hits, "")) science_hits <- character()
add_check(
  "NO_SCIENCE_PLACEHOLDERS",
  "投稿权威面",
  "无科学结果、外部验证、图表或待补文献占位",
  if (length(science_hits)) paste(unique(science_hits), collapse = "; ") else "0",
  status = status_if(length(science_hits) == 0L),
  action = "清除所有科学占位后才能进入目标期刊格式修订。"
)

admin_hits <- unique(unlist(regmatches(
  paste(master, declarations, sep = "\n"),
  gregexpr("\\{\\{[^}]+\\}\\}", paste(master, declarations, sep = "\n"), perl = TRUE)
)))
admin_hits <- admin_hits[nzchar(admin_hits)]
add_check(
  "ADMIN_FIELDS",
  "投稿行政信息",
  "作者、单位、基金、利益冲突、贡献、致谢和代码归档信息由作者确认",
  if (length(admin_hits)) paste(admin_hits, collapse = "; ") else "已全部填写",
  gate_class = "administrative",
  status = if (length(admin_hits)) "WARN" else "PASS",
  blocks_target_formatting = FALSE,
  blocks_final_submission = TRUE,
  action = "可先套目标期刊格式，但最终上传前必须由作者逐项确认。"
)

if (file.exists(human_completion_file)) {
  human_completion <- fread(human_completion_file, sep = "\t", encoding = "UTF-8")
  required_human_columns <- c("item_id", "blocking_level")
  if (all(required_human_columns %in% names(human_completion))) {
    open_human_items <- human_completion[
      blocking_level %chin% c("FINAL_UPLOAD_BLOCKER", "SUBMISSION_DAY_CHECK")
    ]
    human_observed <- if (nrow(open_human_items)) {
      paste0(
        open_human_items$item_id, "(", open_human_items$blocking_level, ")",
        collapse = "; "
      )
    } else {
      "所有阻断项与投稿当日检查均已完成"
    }
    human_status <- if (nrow(open_human_items)) "WARN" else "PASS"
  } else {
    human_observed <- paste0(
      "清单字段不完整：缺少",
      paste(setdiff(required_human_columns, names(human_completion)), collapse = ", ")
    )
    human_status <- "WARN"
  }
} else {
  human_observed <- "人工补齐清单缺失"
  human_status <- "WARN"
}
add_check(
  "HUMAN_COMPLETION_ITEMS",
  "投稿行政信息",
  "Elsevier声明、代码归档、投稿系统字段及投稿当日动态检查均已完成",
  human_observed,
  gate_class = "administrative",
  status = human_status,
  blocks_target_formatting = FALSE,
  blocks_final_submission = TRUE,
  action = "按人工补齐清单逐项完成并将阻断项更新为READY；OPTIONAL项不阻断投稿。"
)

target_journal <- Sys.getenv("TARGET_JOURNAL", unset = "")
add_check(
  "TARGET_JOURNAL",
  "期刊行政信息",
  "已冻结目标期刊及其最新作者指南",
  if (nzchar(target_journal)) target_journal else "未指定",
  gate_class = "target_journal",
  status = if (nzchar(target_journal)) "PASS" else "WARN",
  blocks_target_formatting = FALSE,
  blocks_final_submission = TRUE,
  action = "选刊后按最新作者指南调整章节、字数、图表和声明位置。"
)

local_leaks <- c("/Users/", "file://", "ResearchDataHub", "_work/")
leak_found <- local_leaks[vapply(local_leaks, grepl, logical(1), x = formal_text, fixed = TRUE)]
add_check(
  "NO_LOCAL_PATH_LEAKAGE",
  "投稿权威面",
  "投稿文本不暴露本机路径或过程目录",
  if (length(leak_found)) paste(leak_found, collapse = "; ") else "0",
  status = status_if(length(leak_found) == 0L),
  action = "用公共 accession、DOI 或代码仓库地址替代本机路径。"
)

internal_terms <- c(
  "formal upstream result", "quality gate", "claim ceiling", "non-negotiable",
  "pending external", "not yet generated", "final_final", "AI-generated",
  "artifact role", "manifest-verified", "authoritative", "upload instruction",
  "sha256", "results/", "artifact_role", "display_catalog_file",
  "source_files", "source_file", "source_manifest_files", "assembly_status",
  "upload_instruction", "data_rewritten", "source_row_key",
  "artifact manifest", "visual re-opening", "stage 1 go", "proxy fidelity go",
  "evidence gate", "extension gate", "discrete-cluster gate",
  "interpretation ceiling"
)
submission_science_text_lower <- tolower(paste(
  master, main_tables, supp_tables, submission_table_text
))
internal_found <- internal_terms[vapply(
  tolower(internal_terms),
  function(term) grepl(term, submission_science_text_lower, fixed = TRUE),
  logical(1)
)]
add_check(
  "NO_INTERNAL_PROCESS_LANGUAGE",
  "投稿权威面",
  "正文与投稿表格无内部流程词或项目演化痕迹",
  if (length(internal_found)) paste(internal_found, collapse = "; ") else "0",
  status = status_if(length(internal_found) == 0L),
  action = "改为常规科研语言；过程信息仅保留在 _work。"
)

abstract_text <- sub("(?s).*?## Abstract\\s*", "", master, perl = TRUE)
abstract_text <- sub("(?s)\\s*## Keywords.*", "", abstract_text, perl = TRUE)
abstract_words <- strsplit(gsub("[^[:alnum:]'–-]+", " ", abstract_text), "\\s+")[[1L]]
abstract_words <- abstract_words[nzchar(abstract_words)]
add_check(
  "ABSTRACT_LENGTH",
  "正文结构",
  "结构化套版前摘要不超过250词",
  length(abstract_words),
  status = status_if(length(abstract_words) <= 250L),
  action = "目标期刊若限制更严，再按其要求压缩。"
)

if (file.exists(references_file)) {
  citation_keys <- unique(unlist(regmatches(
    master,
    gregexpr("(?<![A-Za-z0-9._%+-])@[A-Za-z0-9_:.-]+", master, perl = TRUE)
  )))
  citation_keys <- sub("^@", "", citation_keys)
  bib <- read_text(references_file)
  bib_keys <- unique(sub(
    "^@[A-Za-z]+\\{([^,]+),.*$", "\\1",
    grep("^@[A-Za-z]+\\{[^,]+,", strsplit(bib, "\n", fixed = TRUE)[[1L]], value = TRUE)
  ))
  missing_bib <- setdiff(citation_keys, bib_keys)
  unused_bib <- setdiff(bib_keys, citation_keys)
  cite_ok <- length(citation_keys) > 0L && length(missing_bib) == 0L && length(unused_bib) == 0L
  cite_observed <- sprintf(
    "正文键=%d; BibTeX=%d; 缺失=%d; 未使用=%d",
    length(citation_keys), length(bib_keys), length(missing_bib), length(unused_bib)
  )
} else {
  cite_ok <- FALSE
  cite_observed <- "references.bib缺失"
}
add_check(
  "CITATION_COVERAGE",
  "参考文献",
  "正文引用键与BibTeX一一对应",
  cite_observed,
  status = status_if(cite_ok),
  action = "先修复缺失、重复或未使用键，再进行期刊引用样式转换。"
)

old_chinese_files <- c(
  path("manuscript", "escc_multiomics_results_and_discussion.md"),
  path("manuscript", "escc_multiomics_figure_legends.md")
)
old_authority_ok <- all(file.exists(old_chinese_files)) && all(vapply(
  old_chinese_files,
  function(file) grepl("不属于投稿包", read_text(file), fixed = TRUE),
  logical(1)
)) && !file.exists(path("manuscript", "escc_multiomics_supplementary_inventory.md"))
add_check(
  "AUTHORITY_CONVERGENCE",
  "投稿权威面",
  "旧中文稿明确排除，内部补充计划不留在投稿目录",
  if (old_authority_ok) "英文母稿为唯一正文权威面" else "权威面仍未收敛",
  status = status_if(old_authority_ok),
  action = "保留中文稿仅作核对；内部清单移入 _work/checks。"
)

external_manifest_file <- path("results", "escc_external_validation_artifact_manifest.tsv")
external_manifest_ok <- FALSE
external_manifest_observed <- "缺失"
if (file.exists(external_manifest_file)) {
  external_manifest <- fread(external_manifest_file)
  required_cols <- c("relative_path", "file_size_bytes", "sha256", "status")
  if (all(required_cols %in% names(external_manifest))) {
    actual_exists <- file.exists(path(external_manifest$relative_path))
    actual_sizes <- rep(NA_real_, nrow(external_manifest))
    actual_hashes <- rep(NA_character_, nrow(external_manifest))
    actual_sizes[actual_exists] <- file.info(path(external_manifest$relative_path[actual_exists]))$size
    actual_hashes[actual_exists] <- vapply(
      path(external_manifest$relative_path[actual_exists]), sha256_file, character(1)
    )
    external_manifest_ok <- all(actual_exists) &&
      all(actual_sizes == external_manifest$file_size_bytes) &&
      all(actual_hashes == external_manifest$sha256) &&
      all(external_manifest$status == "verified")
    external_manifest_observed <- sprintf(
      "%d artifacts; missing=%d; hash/size/status mismatch=%d",
      nrow(external_manifest), sum(!actual_exists),
      sum(!actual_exists | actual_sizes != external_manifest$file_size_bytes |
            actual_hashes != external_manifest$sha256 |
            external_manifest$status != "verified", na.rm = TRUE)
    )
  }
}
add_check(
  "EXTERNAL_ARTIFACTS",
  "外部验证",
  "外部验证正式产物、大小、SHA256和状态全部一致",
  external_manifest_observed,
  status = status_if(external_manifest_ok),
  action = "重跑 script25 并修复所有 manifest 不一致。"
)

external_tables <- c(
  cv = path("results", "escc_external_state_internal_cv.tsv"),
  survival = path("results", "escc_external_state_survival_associations.tsv"),
  increment = path("results", "escc_external_state_ecms_increment.tsv"),
  response = path("results", "escc_external_state_response_associations.tsv"),
  decision = path("results", "escc_external_validation_decision.tsv")
)
external_contract_ok <- all(file.exists(external_tables))
external_contract_observed <- "外部表缺失"
if (external_contract_ok) {
  cv <- fread(external_tables[["cv"]])
  survival <- fread(external_tables[["survival"]])
  increment <- fread(external_tables[["increment"]])
  response <- fread(external_tables[["response"]])
  decision <- fread(external_tables[["decision"]])
  cv_main <- cv[record_type == "outer_repeat_summary" & model_variant == "ridge314_primary"]
  survival_main <- survival[
    cohort_scope == "GSE53622_GSE53624_stratified" &
      model_variant == "ridge314_primary" & model_type == "clinical_adjusted"
  ]
  increment_main <- increment[model_variant == "ridge314_primary"]
  response_main <- response[model_variant == "ridge314_primary"]
  overall <- decision[decision_id == "OVERALL_FIRST_STAGE_EXTERNAL_VALIDATION"]
  external_contract_ok <- nrow(cv_main) == 2L && all(cv_main$gate_status == "GO") &&
    all(cv_main$n == 78L) &&
    nrow(survival_main) == 2L && all(survival_main$n_patients == 179L) &&
    all(survival_main$n_events == 106L) &&
    nrow(increment_main) == 2L &&
    all(increment_main$bootstrap_optimism_replicates_valid == 300L) &&
    all(increment_main$increment_status == "no_incremental_support") &&
    all(increment_main$optimism_corrected_delta_cindex_ci_lower_95 < 0) &&
    all(increment_main$optimism_corrected_delta_cindex_ci_upper_95 > 0) &&
    nrow(response_main) == 2L && all(response_main$n_total == 28L) &&
    all(response_main$n_pathological_complete_response == 11L) &&
    all(response_main$endpoint_role == "exploratory_only") &&
    nrow(overall) == 1L &&
    overall$status == "GO_FIRST_STAGE_EXTERNAL_RNA_PROXY_COMPLETED"
  external_contract_observed <- sprintf(
    "CV rho F1/F3=%.3f/%.3f; pooled=%d/%d events; bootstrap=%d/%d; pCR=%d/%d",
    cv_main[factor == "Factor1", spearman_rho],
    cv_main[factor == "Factor3", spearman_rho],
    survival_main[1L, n_patients], survival_main[1L, n_events],
    increment_main[factor == "Factor1", bootstrap_optimism_replicates_valid],
    increment_main[factor == "Factor3", bootstrap_optimism_replicates_valid],
    response_main[1L, n_pathological_complete_response], response_main[1L, n_total]
  )
}
add_check(
  "EXTERNAL_CONTRACT",
  "外部验证",
  "proxy fidelity、生存、真正乐观偏差校正、pCR和总体边界符合冻结合同",
  external_contract_observed,
  status = status_if(external_contract_ok),
  action = "不得把内部CV、共享RNA增量或28例pCR写成独立机制/临床验证。"
)

critical_manuscript_tokens <- c(
  "0.839", "0.896", "1.064", "1.073", "−0.0022", "0.0023",
  "0.455", "2.151", "GSE53622", "GSE53624", "GSE45670",
  "do not provide exact independent replication"
)
missing_tokens <- critical_manuscript_tokens[!vapply(
  critical_manuscript_tokens, grepl, logical(1), x = master, fixed = TRUE
)]
add_check(
  "EXTERNAL_MANUSCRIPT_MATCH",
  "外部验证",
  "母稿包含冻结的关键效应和证据边界",
  if (length(missing_tokens)) paste(missing_tokens, collapse = "; ") else "关键值与边界齐全",
  status = status_if(length(missing_tokens) == 0L),
  action = "从正式外部表重新填充摘要、方法、结果、讨论、Table 3和Figure 6图注。"
)

ledger_file <- path("results", "escc_multiomics_evidence_ledger.tsv")
ledger_ok <- FALSE
ledger_observed <- "缺失"
if (file.exists(ledger_file)) {
  ledger <- fread(ledger_file)
  tier_counts <- ledger[, .N, by = evidence_tier][order(evidence_tier)]
  expected <- data.table(
    evidence_tier = c("T0", "T1", "T2", "T3", "T4"),
    N = c(40L, 118L, 1508L, 22L, 1L)
  )
  ledger_ok <- identical(tier_counts$evidence_tier, expected$evidence_tier) &&
    identical(tier_counts$N, expected$N)
  ledger_observed <- paste(
    paste(tier_counts$evidence_tier, tier_counts$N, sep = "="),
    collapse = "; "
  )
}
add_check(
  "LEDGER_COUNTS",
  "整合结果",
  "证据台账层级与母稿冻结数字一致",
  ledger_observed,
  status = status_if(ledger_ok),
  action = "以正式 ledger 为准同步母稿、主表和索引。"
)

figure_manifest_specs <- list(
  list(
    file = path("results", "escc_multiomics_figure_artifact_manifest.tsv"),
    expected_rows = 20L,
    expected_visual = "passed_reopened_review",
    label = "Figure1-5"
  ),
  list(
    file = path("results", "escc_external_validation_figure_artifact_manifest.tsv"),
    expected_rows = 4L,
    expected_visual = "passed_reopened_review",
    label = "Figure6"
  )
)
figure_ok <- TRUE
figure_notes <- character()
for (spec in figure_manifest_specs) {
  ok <- FALSE
  note <- paste0(spec$label, " manifest缺失")
  if (file.exists(spec$file)) {
    manifest <- fread(spec$file)
    cols_ok <- all(c(
      "relative_path", "file_size_bytes", "sha256", "structural_status", "visual_qa_status"
    ) %in% names(manifest))
    if (cols_ok) {
      exists <- file.exists(path(manifest$relative_path))
      actual_size <- rep(NA_real_, nrow(manifest))
      actual_sha <- rep(NA_character_, nrow(manifest))
      actual_size[exists] <- file.info(path(manifest$relative_path[exists]))$size
      actual_sha[exists] <- vapply(path(manifest$relative_path[exists]), sha256_file, character(1))
      ok <- nrow(manifest) == spec$expected_rows && all(exists) &&
        all(actual_size == manifest$file_size_bytes) && all(actual_sha == manifest$sha256) &&
        all(grepl("verified", manifest$structural_status, fixed = TRUE)) &&
        all(manifest$visual_qa_status == spec$expected_visual) &&
        all(nzchar(manifest$qa_path)) && all(nzchar(manifest$qa_sha256))
      note <- sprintf(
        "%s rows=%d; visual=%s; missing=%d",
        spec$label, nrow(manifest), paste(unique(manifest$visual_qa_status), collapse = ","), sum(!exists)
      )
    }
  }
  figure_ok <- figure_ok && ok
  figure_notes <- c(figure_notes, note)
}
add_check(
  "FIGURE_PACKAGE",
  "图件",
  "Figure1-6四格式、SHA256、结构检查与逐图重开视觉QA全部通过",
  paste(figure_notes, collapse = "; "),
  status = status_if(figure_ok),
  action = "图件只能在人工重开通过后用 finalize 参数冻结 manifest。"
)

table_files <- c(
  submission_table_files,
  path("results", "escc_manuscript_table_artifact_manifest.tsv"),
  main_tables_file,
  supp_tables_file
)
tables_exist <- file.exists(table_files)
table_manifest_ok <- all(tables_exist)
table_observed <- paste(basename(table_files), tables_exist, sep = "=", collapse = "; ")
if (table_manifest_ok) {
  table_manifest <- fread(path("results", "escc_manuscript_table_artifact_manifest.tsv"))
  required_cols <- c("relative_path", "file_size_bytes", "sha256", "status")
  if (all(required_cols %in% names(table_manifest))) {
    exists <- file.exists(path(table_manifest$relative_path))
    actual_size <- rep(NA_real_, nrow(table_manifest))
    actual_sha <- rep(NA_character_, nrow(table_manifest))
    actual_size[exists] <- file.info(path(table_manifest$relative_path[exists]))$size
    actual_sha[exists] <- vapply(path(table_manifest$relative_path[exists]), sha256_file, character(1))
    table_manifest_ok <- all(exists) &&
      all(actual_size == table_manifest$file_size_bytes) &&
      all(actual_sha == table_manifest$sha256) &&
      all(table_manifest$status == "verified")
    table_observed <- sprintf(
      "%d manifest rows; missing=%d; bad_status=%d",
      nrow(table_manifest), sum(!exists), sum(table_manifest$status != "verified")
    )
  } else {
    table_manifest_ok <- FALSE
    table_observed <- "table manifest字段不完整"
  }
}
add_check(
  "TABLE_PACKAGE",
  "表格与补充材料",
  "主表、补充表、机器可读表和manifest均已组装并通过SHA256",
  table_observed,
  status = status_if(table_manifest_ok),
  action = "用 script27 从正式结果表重建，不手工复制数值。"
)

clean_supp_ok <- file.exists(clean_supp_index_file) &&
  file.exists(clean_supp_manifest_file) && dir.exists(clean_supp_dir)
clean_supp_observed <- "净化附件索引、内部 manifest 或附件目录缺失"
clean_supp_files <- character()
clean_supp_basenames <- character()
if (clean_supp_ok) {
  clean_index <- fread(clean_supp_index_file, encoding = "UTF-8")
  clean_manifest <- fread(clean_supp_manifest_file, encoding = "UTF-8")
  index_cols <- c(
    "supplementary_table_id", "title", "clean_attachment_basename", "rows",
    "columns", "source_basename", "interpretation_boundary"
  )
  manifest_cols <- c(
    "source_basename", "clean_attachment_basename", "relative_path",
    "source_sha256", "output_sha256", "file_size_bytes", "output_rows",
    "output_columns", "execution_script_sha256", "status"
  )
  clean_supp_ok <- all(index_cols %in% names(clean_index)) &&
    all(manifest_cols %in% names(clean_manifest)) &&
    nrow(clean_index) == 71L && nrow(clean_manifest) == 71L &&
    !anyDuplicated(clean_index$source_basename) &&
    !anyDuplicated(clean_index$clean_attachment_basename) &&
    !anyDuplicated(clean_manifest$relative_path)
  scan_hits <- character()
  missing_count <- 71L
  mismatch_count <- 71L
  if (clean_supp_ok) {
    setkey(clean_index, source_basename, clean_attachment_basename)
    setkey(clean_manifest, source_basename, clean_attachment_basename)
    clean_supp_ok <- identical(
      clean_index[, .(source_basename, clean_attachment_basename)],
      clean_manifest[, .(source_basename, clean_attachment_basename)]
    ) && all(clean_manifest$status == "verified_submission_clean") &&
      !any(clean_manifest$source_basename %in% c(
        "escc_external_validation_decision.tsv", "tcga_escc_ecms_projection_qa.tsv"
      ))
  }
  if (clean_supp_ok) {
    clean_supp_files <- path(clean_manifest$relative_path)
    clean_supp_basenames <- clean_manifest$clean_attachment_basename
    exists <- file.exists(clean_supp_files)
    missing_count <- sum(!exists)
    mismatch <- rep(TRUE, nrow(clean_manifest))
    banned_pattern <- paste0(
      "/Users/|ResearchDataHub|results/|_work/|sha256|manifest|artifact|",
      "source_file|source_row_key|(^|_)relative_path$|target_relative_path|",
      "(^|_)run_status$|checkpoint|(^|_)gate(_|$)|_ceiling($|_)|",
      "required[ _]next[ _]validation|decision[ _]basis|stage 1 go|",
      "proxy fidelity go|evidence gate|extension gate|discrete-cluster gate|",
      "countable|whitelist|prelocked|pre-locked|pre_locked"
    )
    if (any(exists)) {
      for (i in which(exists)) {
        f <- clean_supp_files[i]
        z <- fread(f, encoding = "UTF-8", showProgress = FALSE)
        chars <- unlist(z[, which(vapply(z, is.character, logical(1))), with = FALSE])
        tokens <- c(names(z), chars)
        han <- any(grepl("[一-龥]", tokens, perl = TRUE), na.rm = TRUE)
        internal <- any(grepl(banned_pattern, tokens, ignore.case = TRUE, perl = TRUE),
                        na.rm = TRUE)
        if (han) scan_hits <- c(scan_hits, paste0(basename(f), ":han"))
        if (internal) scan_hits <- c(scan_hits, paste0(basename(f), ":internal"))
        source_file <- path("results", clean_manifest$source_basename[i])
        mismatch[i] <- !file.exists(source_file) ||
          file.info(f)$size != clean_manifest$file_size_bytes[i] ||
          sha256_file(f) != clean_manifest$output_sha256[i] ||
          sha256_file(source_file) != clean_manifest$source_sha256[i] ||
          nrow(z) != clean_manifest$output_rows[i] ||
          ncol(z) != clean_manifest$output_columns[i] ||
          nrow(z) < 1L || ncol(z) < 1L
      }
    }
    mismatch_count <- sum(mismatch)
    actual_tsv <- list.files(clean_supp_dir, pattern = "\\.tsv$", full.names = TRUE)
    expected_paths <- normalizePath(clean_supp_files, winslash = "/", mustWork = FALSE)
    actual_paths <- normalizePath(actual_tsv, winslash = "/", mustWork = FALSE)
    directory_ok <- length(actual_tsv) == 71L && setequal(actual_paths, expected_paths)
    script29 <- path("scripts", "29_prepare_submission_supplement.R")
    script_ok <- file.exists(script29) &&
      all(clean_manifest$execution_script_sha256 == sha256_file(script29))
    index_dimension_ok <- all(clean_index$rows == clean_manifest$output_rows) &&
      all(clean_index$columns == clean_manifest$output_columns)
    clean_supp_ok <- all(exists) && mismatch_count == 0L &&
      length(scan_hits) == 0L && directory_ok && script_ok && index_dimension_ok
    clean_supp_observed <- sprintf(
      "attachments=%d; total_bytes=%d; missing=%d; hash/dimension mismatch=%d; QA hits=%d",
      nrow(clean_manifest), sum(clean_manifest$file_size_bytes), missing_count,
      mismatch_count, length(scan_hits)
    )
  }
}
add_check(
  "SUBMISSION_SUPPLEMENT_PACKAGE",
  "表格与补充材料",
  "71个英文净化附件均存在、哈希与script29 manifest一致，且无汉字、本机路径或内部裁决词",
  clean_supp_observed,
  gate_class = "submission_package",
  status = status_if(clean_supp_ok),
  blocks_target_formatting = FALSE,
  blocks_final_submission = TRUE,
  action = "可先继续目标期刊套版；最终上传前必须重跑script29并清零全部附件QA命中。"
)

novelty_files <- c(
  path("results", "escc_multiomics_novelty_gate.tsv"),
  path("results", "escc_multiomics_novelty_summary.md")
)
novelty_ok <- all(file.exists(novelty_files)) &&
  grepl("CONDITIONAL GO", read_text(novelty_files[[2L]]), fixed = TRUE)
add_check(
  "NOVELTY_GATE",
  "新颖性",
  "新颖性已冻结为可防守的条件性GO且无first-ever过度宣称",
  if (novelty_ok) "CONDITIONAL GO" else "新颖性文件缺失或状态未冻结",
  status = status_if(novelty_ok),
  action = "创新点限定为表示重叠审计、证据依赖架构与连续状态外部代理校准。"
)

gse151838_absent <- !grepl("GSE151838", paste(master, declarations), fixed = TRUE)
add_check(
  "DATA_AVAILABILITY_SCOPE",
  "数据可用性",
  "数据可用性不把未进入正式分析的GSE151838写成已分析数据",
  if (gse151838_absent) "GSE151838未列入已分析数据" else "发现GSE151838",
  status = status_if(gse151838_absent),
  action = "仅列GSE149608、GSE149609及正式使用的蛋白补充表。"
)

formatting_failures <- checks[
  blocks_target_formatting == TRUE & status == "FAIL", .N
]
submission_blockers <- checks[
  blocks_final_submission == TRUE & status != "PASS", .N
]
formatting_status <- if (formatting_failures == 0L) {
  "READY_FOR_TARGET_JOURNAL_FORMATTING"
} else {
  "NOT_READY_FOR_TARGET_JOURNAL_FORMATTING"
}
submission_status <- if (submission_blockers == 0L) {
  "READY_FOR_FINAL_SUBMISSION"
} else {
  "CONDITIONAL_NOT_READY_FOR_FINAL_UPLOAD"
}

add_check(
  "OVERALL_TARGET_FORMATTING",
  "总体判定",
  "科学结果、图表、引用和证据边界足以开始目标期刊格式修订",
  formatting_status,
  gate_class = "overall",
  status = if (formatting_failures == 0L) "PASS" else "FAIL",
  blocks_target_formatting = TRUE,
  blocks_final_submission = TRUE,
  action = "只有科学硬失败清零后才开始套版。"
)
add_check(
  "OVERALL_FINAL_UPLOAD",
  "总体判定",
  "作者行政字段、目标期刊和代码归档均完成后方可最终上传",
  submission_status,
  gate_class = "overall",
  status = if (submission_blockers == 0L) "PASS" else "WARN",
  blocks_target_formatting = FALSE,
  blocks_final_submission = TRUE,
  action = "格式修订可以先开始；最终上传前完成所有行政WARN。"
)

checks[, domain_order__ := match(domain, c(
  "总体判定", "投稿权威面", "外部验证", "整合结果", "图件",
  "表格与补充材料", "参考文献", "新颖性", "数据可用性",
  "正文结构", "投稿行政信息", "期刊行政信息"
))]
setorder(checks, domain_order__, check_id)
checks[, domain_order__ := NULL]

submission_files <- unique(c(
  required_authority_files,
  list.files(path("figures"), pattern = "^escc_multiomics_figure[1-6]_.*\\.(svg|pdf|tiff|png)$", full.names = TRUE),
  submission_table_files,
  clean_supp_files
))
submission_files <- submission_files[file.exists(submission_files)]
package_manifest <- data.table(
  artifact = basename(submission_files),
  relative_path = sub(paste0("^", root, "/"), "", submission_files),
  role = fifelse(
    basename(submission_files) %in% clean_supp_basenames, "supplementary_data_table",
    fifelse(
      grepl("^escc_multiomics_figure", basename(submission_files)), "main_figure",
      fifelse(
        grepl("supplementary", basename(submission_files)), "supplementary_material",
        fifelse(
          grepl("table[123]", basename(submission_files)), "main_table",
          fifelse(
            basename(submission_files) == "references.bib", "bibliography",
            fifelse(
              grepl("declarations", basename(submission_files)), "submission_declarations",
              "manuscript_or_manifest"
            )
          )
        )
      )
    )
  ),
  file_size_bytes = file.info(submission_files)$size,
  sha256 = vapply(submission_files, sha256_file, character(1)),
  generated_date = as.character(Sys.Date()),
  status = "verified_current_bytes"
)
setorder(package_manifest, role, artifact)

summary_lines <- c(
  "# ESCC 多组学投稿准备度终审",
  "",
  "## 总体结论",
  "",
  paste0("- 目标期刊格式修订：`", formatting_status, "`。"),
  paste0("- 最终投稿上传：`", submission_status, "`。"),
  paste0("- 科学/图表硬失败数：", formatting_failures, "。"),
  paste0("- 最终上传待处理门禁项：", submission_blockers, "。"),
  "",
  "格式修订门禁与最终上传门禁分开判定。已确认的作者、单位、基金、利益冲突、贡献和致谢，与仍待完成的代码归档、Elsevier声明、投稿系统字段及投稿当日检查均由行政门禁独立追踪；任何行政事项都不得掩盖科学结果、图表、引用或证据边界的缺口。",
  "",
  "## 已冻结的证据上限",
  "",
  "- 9 条事件—连续状态关联仍属于 TCGA 来源内、表示重叠敏感性审计后的候选关联。",
  "- Factor1/Factor3 的 314 基因 RNA 代理已完成第一阶段外部应用；生存和 ECMS 条件增量为阴性，28 例 pCR 仅探索性。",
  "- 外部表达队列不能重建完整 MOFA 因子，也不能复现精确突变/CNV 事件—状态边。",
  "- Cao、代谢组和微生物组模块保持各自分析单位与不匹配边界，不构造跨队列因果链。",
  "",
  "## 未通过或待作者确认的项目",
  ""
)
open_items <- checks[status != "PASS"]
if (nrow(open_items) == 0L) {
  summary_lines <- c(summary_lines, "- 无。")
} else {
  summary_lines <- c(
    summary_lines,
    paste0(
      "- `", open_items$check_id, "`（", open_items$status, "）：",
      open_items$observed, "；", open_items$action
    )
  )
}
summary_lines <- c(
  summary_lines,
  "",
  "## 权威入口",
  "",
  "- 英文母稿：`manuscript/escc_multiomics_manuscript_master.md`。",
  "- 主表：`manuscript/escc_multiomics_main_tables.md`。",
  "- 补充表：`manuscript/escc_multiomics_supplementary_tables.md`。",
  "- 投稿声明：`manuscript/escc_multiomics_submission_declarations.md`。",
  "- 逐项门禁：`results/escc_submission_readiness_gate.tsv`。",
  "- 投稿包文件清单：`results/escc_submission_package_manifest.tsv`。"
)

gate_output <- path("results", "escc_submission_readiness_gate.tsv")
package_output <- path("results", "escc_submission_package_manifest.tsv")
summary_output <- path("results", "escc_submission_readiness_summary.md")

if (!validate_only) {
  temp_gate <- tempfile(pattern = ".submission_readiness_gate_", tmpdir = path("_work/intermediate"))
  temp_package <- tempfile(pattern = ".submission_package_manifest_", tmpdir = path("_work/intermediate"))
  temp_summary <- tempfile(pattern = ".submission_readiness_summary_", tmpdir = path("_work/intermediate"))
  on.exit(unlink(c(temp_gate, temp_package, temp_summary), force = TRUE), add = TRUE)
  fwrite(checks, temp_gate, sep = "\t", quote = FALSE, na = "")
  fwrite(package_manifest, temp_package, sep = "\t", quote = FALSE, na = "")
  writeLines(summary_lines, temp_summary, useBytes = TRUE)
  if (!file.rename(temp_gate, gate_output) ||
      !file.rename(temp_package, package_output) ||
      !file.rename(temp_summary, summary_output)) {
    stop("投稿准备度结果原子发布失败。", call. = FALSE)
  }
}

cat("目标期刊格式修订：", formatting_status, "\n", sep = "")
cat("最终投稿上传：", submission_status, "\n", sep = "")
cat("科学/图表硬失败数：", formatting_failures, "\n", sep = "")
cat("最终上传待处理门禁项：", submission_blockers, "\n", sep = "")

if (validate_only && formatting_failures > 0L) quit(status = 1L)
