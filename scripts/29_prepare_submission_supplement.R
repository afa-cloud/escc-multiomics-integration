#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(fs)
})

# 从正式上游结果生成英文、无本机追溯信息的投稿补充附件。
# 上游 results 表只读；纯内部决策/QA 表明确排除。

OUT_DIR <- "results/submission_supplement"
OUT_INDEX <- "results/escc_submission_supplement_index.tsv"
OUT_MANIFEST <- "results/escc_submission_supplement_artifact_manifest.tsv"
SCRIPT_REL <- "scripts/29_prepare_submission_supplement.R"

args <- commandArgs(trailingOnly = TRUE)
allowed <- c("--validate-only", "--fields-only")
bad <- setdiff(args, allowed)
if (length(bad)) stop("未知参数: ", paste(bad, collapse = ", "), call. = FALSE)
if ("--fields-only" %in% args) {
  cat(OUT_DIR, OUT_INDEX, OUT_MANIFEST, sep = "\n")
  quit(save = "no", status = 0)
}
validate_only <- "--validate-only" %in% args

full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", full_args, value = TRUE)
script_abs <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), SCRIPT_REL), mustWork = TRUE)
}
root <- normalizePath(file.path(dirname(script_abs), ".."), mustWork = TRUE)
setwd(root)

fail <- function(...) stop(paste0(...), call. = FALSE)
ok <- function(x, msg) {
  if (length(x) != 1L || is.na(x) || !isTRUE(x)) fail(msg)
  invisible(TRUE)
}
need <- function(x, cols, label) {
  miss <- setdiff(cols, names(x))
  if (length(miss)) fail(label, " 缺少字段: ", paste(miss, collapse = ", "))
}
sha <- function(path) {
  ok(file_exists(path), paste0("文件不存在: ", path))
  digest(file = path, algo = "sha256", serialize = FALSE)
}
bytes <- function(path) as.numeric(file_info(path)$size)
has_han <- function(x) any(grepl("[一-龥]", as.character(x), perl = TRUE), na.rm = TRUE)

GROUP_SOURCES <- list(
  S1 = c("escc_multiomics_module_summaries.tsv", "tcga_escc_multiassay_analysis_sets.tsv", "tcga_escc_ecms_patient_probabilities.tsv", "escc_external_validation_decision.tsv"),
  S2 = c("tcga_escc_multiassay_sample_map.tsv", "tcga_escc_multiassay_analysis_sets.tsv", "tcga_escc_multiassay_qc.tsv"),
  S3 = c("tcga_escc_driver_candidate_screen.tsv", "tcga_escc_strong_driver_patient_events.tsv", "tcga_escc_driver_core_qc.tsv"),
  S4 = c("tcga_escc_cnv_expression_dosage.tsv", "tcga_escc_cnv_gene_summary.tsv", "tcga_escc_mutation_gene_summary.tsv", "tcga_escc_mutation_cluster_summary.tsv", "tcga_escc_mutation_pairwise_interactions.tsv"),
  S5 = c("tcga_escc_mofa_variance_explained.tsv", "tcga_escc_mofa_factor_scores.tsv", "tcga_escc_mofa_top_weights.tsv", "tcga_escc_mofa_pathway_associations.tsv", "tcga_escc_mofa_level_factor_qc.tsv", "tcga_escc_mofa_clinical_associations.tsv", "tcga_escc_progeny_pathway_scores.tsv"),
  S6 = c("tcga_escc_heterogeneity_cluster_evaluation.tsv", "tcga_escc_heterogeneity_cluster_pathways.tsv", "tcga_escc_heterogeneity_patient_assignments.tsv"),
  S7 = c("tcga_escc_driver_state_associations.tsv", "tcga_escc_driver_state_network_edges.tsv", "tcga_escc_driver_state_leakage_factor_alignment.tsv", "tcga_escc_driver_state_leakage_associations.tsv", "tcga_escc_driver_state_leakage_edge_stability.tsv", "tcga_escc_driver_state_leakage_candidate_summary.tsv"),
  S8 = c("tcga_escc_ecms_patient_probabilities.tsv", "tcga_escc_ecms_projection_calibration.tsv", "tcga_escc_ecms_projection_qa.tsv", "tcga_escc_ecms_factor_associations.tsv", "tcga_escc_ecms_adjusted_factor_progeny_associations.tsv"),
  S9 = c("cao2020_wgbs_candidate_pair_effects.tsv", "cao2020_wgbs_candidate_region_summary.tsv", "cao2020_rna_candidate_pair_effects.tsv", "cao2020_rna_candidate_summary.tsv", "cao2020_proteomics_candidate_pair_effects.tsv", "cao2020_proteomics_candidate_summary.tsv", "cao2020_cross_layer_patient_effects.tsv", "cao2020_cross_layer_candidate_summary.tsv"),
  S10 = c("pr001876_targeted_ms_analysis_inventory.tsv", "pr001876_targeted_ms_sample_qc.tsv", "pr001876_targeted_ms_feature_qc.tsv", "pr001876_targeted_ms_differential.tsv", "pr001876_targeted_ms_candidate_metabolites.tsv", "pr001876_targeted_ms_run_order_sensitivity.tsv", "pr001876_targeted_ms_scale_sensitivity.tsv", "pr001876_targeted_ms_paired_sensitivity.tsv"),
  S11 = c("prjna766558_dada2_sample_qc.tsv", "prjna766558_dada2_asv_counts.tsv", "prjna766558_dada2_taxonomy.tsv", "prjna766558_dada2_asv_length_qc.tsv", "prjna766558_dada2_alpha_diversity.tsv", "prjna766558_dada2_alpha_paired_tests.tsv", "prjna766558_dada2_beta_tests.tsv", "prjna766558_dada2_genus_paired_differential.tsv", "prjna766558_dada2_ancombc2_sensitivity.tsv", "prjna766558_dada2_pipeline_sensitivity.tsv"),
  S12 = c("escc_multiomics_integration_tier_definitions.tsv", "escc_multiomics_evidence_ledger.tsv", "escc_multiomics_integrated_driver_candidates.tsv", "escc_multiomics_integrated_axis_edges.tsv", "escc_multiomics_integrated_axis_summary.tsv", "escc_multiomics_heterogeneity_axes.tsv", "escc_multiomics_module_summaries.tsv"),
  S13 = c("escc_external_state_proxy_definition.tsv", "escc_external_state_internal_cv.tsv", "escc_external_state_oof_predictions.tsv", "escc_external_state_patient_scores.tsv", "escc_external_state_survival_associations.tsv", "escc_external_state_ecms_increment.tsv", "escc_external_state_response_associations.tsv", "escc_external_validation_decision.tsv")
)

GROUP_TITLES <- c(
  S1 = "Cohort inventory, analysis units, and independence families",
  S2 = "TCGA-ESCC multi-assay sample map and analysis sets",
  S3 = "TCGA-ESCC driver candidate screen and patient-level event matrix",
  S4 = "Somatic mutation and copy-number association details",
  S5 = "MOFA factor, pathway, and clinical annotation results",
  S6 = "Patient-state heterogeneity clustering results",
  S7 = "Representation-overlap-audited driver-state association analysis",
  S8 = "Locked ECMS projection, calibration, and factor associations",
  S9 = "Cao et al. paired cross-layer calibration",
  S10 = "PR001876 targeted metabolomics calibration",
  S11 = "PRJNA766558 paired 16S rRNA microbiome analysis",
  S12 = "Integrated evidence summary, axes, and candidate network",
  S13 = "External RNA-proxy and endpoint transportability calibration"
)

GROUP_BOUNDARIES <- c(
  S1 = "Cohort inventories define analysis units but are not biological validation.",
  S2 = "All patient subsets belong to one overlapping TCGA discovery family.",
  S3 = "Candidate priority is hypothesis-generating and does not establish causality.",
  S4 = "Correlated event encodings from the same TCGA patients are not independent validation.",
  S5 = "Latent-state annotations are association-based and do not establish mechanism.",
  S6 = "Cluster solutions are descriptive and are not independently validated subtypes.",
  S7 = "Omission sensitivity reduces representation overlap but is not independent replication.",
  S8 = "ECMS analyses share the TCGA RNA representation and are contextual calibration only.",
  S9 = "Cross-layer coverage is candidate- and assay-dependent; missingness is not negative evidence.",
  S10 = "Metabolomics results are unpaired and do not support absolute concentration claims.",
  S11 = "Microbiome results are compositional, FFPE-specific, and lack a public negative control.",
  S12 = "Evidence integration ranks hypotheses and does not convert association into mechanism.",
  S13 = "RNA-proxy transportability does not reproduce the full multi-omics model or event-level edges."
)

EXCLUDED <- c(
  "escc_external_validation_decision.tsv",
  "tcga_escc_ecms_projection_qa.tsv"
)

source_long <- rbindlist(lapply(names(GROUP_SOURCES), function(s) {
  data.table(supplementary_table_id = s, source_basename = GROUP_SOURCES[[s]])
}))
source_long[, source_order__ := seq_len(.N)]
source_map <- source_long[!source_basename %chin% EXCLUDED, .(
  supplementary_table_ids = paste(unique(supplementary_table_id), collapse = ";"),
  source_order = min(source_order__)
), by = source_basename]
setorder(source_map, source_order)
source_map[, output_basename := mapply(function(id_string, source) {
  ids <- strsplit(id_string, ";", fixed = TRUE)[[1]]
  prefix <- paste(sprintf("S%02d", as.integer(sub("^S", "", ids))), collapse = "-")
  display_source <- if ("S7" %chin% ids) {
    sub("_leakage_", "_representation_overlap_", source, fixed = TRUE)
  } else {
    source
  }
  paste0(prefix, "__", display_source)
}, supplementary_table_ids, source_basename, USE.NAMES = FALSE)]

ok(nrow(source_long) == 77L && uniqueN(source_long$source_basename) == 73L,
   "冻结的补充来源结构不再是 77 次引用/73 个唯一表。")
ok(nrow(source_map) == 71L && !anyDuplicated(source_map$source_basename) &&
     !anyDuplicated(source_map$output_basename),
   "投稿净化附件必须恰好对应 71 个唯一来源和 71 个唯一输出。")
ok(!any(source_map$source_basename %chin% EXCLUDED), "纯内部决策/QA 表未被排除。")

internal_name_pattern <- paste0(
  "(^|_)(sha256|sha|hash)(_|$)|manifest|artifact|",
  "source_file|source_row_key|source_manifest|",
  "(^|_)run_status$|",
  "countable|whitelist|prelocked|pre_locked|",
  "(^|_)relative_path$|target_relative_path|(^|_)path($|_)|_path$|",
  "(^|_)(local|central)_path($|_)|(^|_)(file|dir|directory)_path($|_)|",
  "(^|_)gate(_|$)|_ceiling($|_)|",
  "required_next_validation|required_validation|required_follow_up|",
  "decision_basis|upgrade_requirement"
)
internal_value_pattern <- paste0(
  "/Users/|file://|ResearchDataHub|results/|_work/|",
  "(^|[^[:alnum:]_])(sha256|sha|hash|manifest|artifact|gate)([^[:alnum:]_]|$)|",
  "(^|_)gate(_|$)|source_file|source_row_key|relative_path|target_relative_path|",
  "_ceiling($|_)|required next validation|decision basis|stage 1 go|proxy fidelity go|",
  "evidence gate|extension gate|discrete-cluster gate|checkpoint|",
  "countable|whitelist|prelocked|pre-locked|pre_locked"
)

protected_pattern <- paste0(
  "(^|_)(p_value|q_value|pval|qval|effect|estimate|hazard_ratio|odds_ratio|",
  "rho|correlation|beta|fold_change|log2fc|auc|cindex)(_|$)|",
  "^(sample|patient|subject|gene|feature|metabolite|asv|genus|taxon|entity|",
  "cohort|factor|event)_(id|name|key)$|^(accession|barcode)$|",
  "(^|_)n_(patients?|samples?|events?|genes?|features?|subjects?|responders?)(_|$)"
)
result_pattern <- paste0(
  "id|name|gene|feature|sample|patient|cohort|factor|pathway|event|effect|",
  "estimate|p_value|q_value|count|score|status|tier|class|label|metric|",
  "observed|abundance|diversity|taxon|taxonomy|asv|genus|metabolite|cnv|",
  "mutation|protein|rna|methyl|survival|hazard|odds|correlation|rho|beta|",
  "fold|auc|cindex|cluster|probability|definition"
)

clinical_boundary_map <- c(
  "患者级探索关联" = "Patient-level exploratory association.",
  "仅 31 个 OS 事件的单变量探索，不作预后定论" =
    "Univariable exploratory overall-survival association based on 31 events; not a definitive prognostic result."
)

column_rename_map <- c(
  countable_as_exact_driver_validation = "eligible_as_exact_driver_replication",
  countable_for_T3_T4 = "meets_T3_T4_retention_criteria",
  countable_as_independent_validation = "eligible_as_independent_validation",
  p_q_countable_as_independent_validation = "p_q_eligible_as_independent_validation",
  countable_as_exact_validation = "eligible_as_exact_replication",
  strongest_state_leakage_countable = "strongest_state_meets_retention_criteria",
  planned_seed_count = "prespecified_seed_count",
  narrative_role = "interpretation_role"
)

exact_value_map <- c(
  clinical_whitelist = "pathology_confirmed_cohort",
  escc_patient_whitelist = "pathology_confirmed_escc_cohort",
  all_four_prelocked_gates = "all_four_prespecified_criteria",
  source_file_count = "input_record_count"
)

tier_label <- c(
  T0 = "Background or not evaluable",
  T1 = "Exploratory evidence",
  T2 = "Conditional candidate",
  T3 = "Robust candidate within the current evidence scope",
  T4 = "High-priority validation candidate"
)
tier_definition <- c(
  T0 = "Technically unevaluable, non-comparable, missing at the relevant layer, or retained only as background.",
  T1 = "Within-source directional or descriptive evidence insufficient for a formal candidate axis.",
  T2 = "Single-cohort, within-cohort bridge, small-sample directional, or context-dependent evidence retained conditionally.",
  T3 = "Patient-level or sensitivity-supported evidence robust within its source but still requiring independent replication.",
  T4 = "Cross-layer or cross-source convergence with a falsifiable validation route; this tier denotes priority rather than established causality."
)

clean_one <- function(source_basename) {
  source_path <- file.path("results", source_basename)
  ok(file_exists(source_path), paste0("上游正式表不存在: ", source_path))
  x <- fread(source_path, encoding = "UTF-8", showProgress = FALSE)
  ok(nrow(x) > 0L && ncol(x) > 0L, paste0("上游表为空: ", source_basename))
  original_names <- names(x)
  removed <- data.table(column = character(), reason = character())
  translated <- character()
  renamed <- character()
  drop_cols <- function(cols, reason) {
    cols <- intersect(cols, names(x))
    if (!length(cols)) return(invisible(NULL))
    removed <<- rbind(removed, data.table(column = cols, reason = reason))
    x[, (cols) := NULL]
    invisible(NULL)
  }

  rename_hits <- intersect(names(column_rename_map), names(x))
  if (length(rename_hits)) {
    new_names <- unname(column_rename_map[rename_hits])
    conflicts <- intersect(new_names, setdiff(names(x), rename_hits))
    ok(!length(conflicts), paste0(
      source_basename, " 列重命名与既有字段冲突: ", paste(conflicts, collapse = ", ")
    ))
    setnames(x, rename_hits, new_names)
    renamed <- paste0(rename_hits, "->", new_names)
  }

  protected <- names(x)[
    grepl(protected_pattern, names(x), ignore.case = TRUE, perl = TRUE) &
      !grepl(internal_name_pattern, names(x), ignore.case = TRUE, perl = TRUE)
  ]

  if (identical(source_basename, "escc_multiomics_integration_tier_definitions.tsv")) {
    need(x, c("evidence_tier", "chinese_label", "operational_definition"), source_basename)
    ok(all(x$evidence_tier %chin% names(tier_label)), "证据层级表出现未冻结层级。")
    x[, evidence_label := unname(tier_label[evidence_tier])]
    x[, operational_definition_english := unname(tier_definition[evidence_tier])]
    x[, interpretation_boundary :=
        "Evidence tiers rank public-data support and do not establish causality or clinical utility."]
    translated <- c(translated, "chinese_label->evidence_label",
                    "operational_definition->operational_definition_english")
    drop_cols(c("chinese_label", "operational_definition"),
              "replaced_by_curated_english_fields")
  }

  if ("boundary" %chin% names(x) && has_han(x$boundary)) {
    vals <- unique(as.character(x$boundary[grepl("[一-龥]", x$boundary, perl = TRUE)]))
    if (length(vals) && all(vals %chin% names(clinical_boundary_map))) {
      hit <- as.character(x$boundary) %chin% names(clinical_boundary_map)
      x[hit, boundary := unname(clinical_boundary_map[as.character(boundary)])]
      translated <- c(translated, "boundary")
    }
  }

  for (column in names(x)[vapply(x, is.character, logical(1))]) {
    for (old_value in names(exact_value_map)) {
      hit <- !is.na(x[[column]]) & x[[column]] == old_value
      if (any(hit)) {
        set(x, which(hit), column, unname(exact_value_map[old_value]))
        translated <- c(
          translated,
          paste0(column, ":", old_value, "->", unname(exact_value_map[old_value]))
        )
      }
    }
  }

  name_internal <- names(x)[grepl(internal_name_pattern, names(x), ignore.case = TRUE, perl = TRUE)]
  drop_cols(name_internal, "internal_provenance_or_project_adjudication_column")

  han_cols <- names(x)[vapply(x, has_han, logical(1))]
  drop_cols(han_cols, "untranslated_chinese_explanatory_column")

  value_internal <- names(x)[vapply(x, function(v) {
    any(grepl(internal_value_pattern, as.character(v), ignore.case = TRUE, perl = TRUE), na.rm = TRUE)
  }, logical(1))]
  drop_cols(value_internal, "internal_path_or_process_content")

  lost_protected <- setdiff(protected, names(x))
  ok(!length(lost_protected), paste0(
    source_basename, " 净化会删除受保护的效应/P/q/标识列: ",
    paste(lost_protected, collapse = ", ")
  ))
  ok(nrow(x) > 0L && ncol(x) > 0L, paste0(source_basename, " 净化后为空。"))
  ok(any(grepl(result_pattern, names(x), ignore.case = TRUE, perl = TRUE)),
     paste0(source_basename, " 净化后无可识别的科学结果或标识列。"))
  ok(!has_han(names(x)) && !any(vapply(x, has_han, logical(1))),
     paste0(source_basename, " 净化后仍含汉字。"))
  ok(!any(grepl(internal_name_pattern, names(x), ignore.case = TRUE, perl = TRUE)),
     paste0(source_basename, " 净化后仍含内部列名。"))
  ok(!any(vapply(x, function(v) {
    any(grepl(internal_value_pattern, as.character(v), ignore.case = TRUE, perl = TRUE), na.rm = TRUE)
  }, logical(1))), paste0(source_basename, " 净化后仍含内部路径或过程词。"))

  if (nrow(removed)) removed <- unique(removed, by = "column")
  list(
    data = x,
    source_path = source_path,
    source_rows = nrow(fread(source_path, select = 1L, showProgress = FALSE)),
    source_columns = length(original_names),
    removed = removed,
    translated = unique(translated),
    renamed = unique(renamed)
  )
}

build_cleaned <- function() {
  lapply(source_map$source_basename, clean_one)
}

make_index <- function(cleaned) {
  rbindlist(lapply(seq_len(nrow(source_map)), function(i) {
    ids <- strsplit(source_map$supplementary_table_ids[i], ";", fixed = TRUE)[[1]]
    data.table(
      supplementary_table_id = paste(ids, collapse = ";"),
      title = paste(unname(GROUP_TITLES[ids]), collapse = " / "),
      clean_attachment_basename = source_map$output_basename[i],
      rows = nrow(cleaned[[i]]$data),
      columns = ncol(cleaned[[i]]$data),
      source_basename = source_map$source_basename[i],
      interpretation_boundary = paste(unique(unname(GROUP_BOUNDARIES[ids])), collapse = " ")
    )
  }))
}

make_manifest <- function(cleaned, output_base, output_is_stage = FALSE) {
  script_sha <- sha(SCRIPT_REL)
  rbindlist(lapply(seq_len(nrow(source_map)), function(i) {
    rel <- file.path(OUT_DIR, source_map$output_basename[i])
    out_path <- if (output_is_stage) file.path(output_base, rel) else rel
    z <- cleaned[[i]]
    rem <- z$removed
    data.table(
      supplementary_table_id = source_map$supplementary_table_ids[i],
      source_basename = source_map$source_basename[i],
      clean_attachment_basename = source_map$output_basename[i],
      relative_path = rel,
      source_sha256 = sha(z$source_path),
      output_sha256 = sha(out_path),
      file_size_bytes = bytes(out_path),
      source_rows = z$source_rows,
      source_columns = z$source_columns,
      output_rows = nrow(z$data),
      output_columns = ncol(z$data),
      removed_column_count = nrow(rem),
      removed_columns = if (nrow(rem)) paste(rem$column, collapse = ";") else "",
      removal_reasons = if (nrow(rem)) paste0(rem$column, "=", rem$reason, collapse = ";") else "none",
      translated_columns = if (length(z$translated)) paste(z$translated, collapse = ";") else "none",
      renamed_columns = if (length(z$renamed)) paste(z$renamed, collapse = ";") else "none",
      generation_script = SCRIPT_REL,
      execution_script_sha256 = script_sha,
      status = "verified_submission_clean"
    )
  }))
}

compare_dt <- function(expected, observed, label) {
  ok(identical(names(expected), names(observed)), paste0(label, " 字段顺序不一致。"))
  ok(nrow(expected) == nrow(observed), paste0(label, " 行数不一致。"))
  z <- all.equal(as.data.frame(expected), as.data.frame(observed),
                 check.attributes = FALSE, tolerance = 1e-12)
  ok(isTRUE(z), paste0(label, " 内容不一致: ", paste(z, collapse = "; ")))
}

validate_index <- function(index) {
  need(index, c(
    "supplementary_table_id", "title", "clean_attachment_basename", "rows",
    "columns", "source_basename", "interpretation_boundary"
  ), "投稿补充附件索引")
  ok(nrow(index) == 71L && !anyDuplicated(index$clean_attachment_basename) &&
       !anyDuplicated(index$source_basename), "投稿补充附件索引必须唯一登记 71 行。")
  vals <- unlist(index[, which(vapply(index, is.character, logical(1))), with = FALSE])
  ok(!has_han(vals), "投稿补充附件索引含汉字。")
  ok(!any(grepl("/Users/|ResearchDataHub|results/|_work/|sha256|manifest|artifact",
                vals, ignore.case = TRUE, perl = TRUE)),
     "投稿补充附件索引含路径、哈希或内部追溯词。")
  ok(!any(index$source_basename %chin% EXCLUDED), "投稿索引误纳入纯内部决策/QA 表。")
}

validate_manifest <- function(manifest, cleaned, base = root) {
  need(manifest, c(
    "supplementary_table_id", "source_basename", "clean_attachment_basename",
    "relative_path", "source_sha256", "output_sha256", "file_size_bytes",
    "source_rows", "source_columns", "output_rows", "output_columns",
    "removed_column_count", "removed_columns", "removal_reasons",
    "translated_columns", "renamed_columns", "generation_script",
    "execution_script_sha256", "status"
  ), "投稿补充附件内部 manifest")
  ok(nrow(manifest) == 71L && !anyDuplicated(manifest$relative_path),
     "投稿补充附件内部 manifest 必须唯一登记 71 行。")
  ok(all(manifest$status == "verified_submission_clean"), "投稿补充附件状态不合格。")
  for (i in seq_len(nrow(manifest))) {
    p <- file.path(base, manifest$relative_path[i])
    ok(file_exists(p), paste0("净化附件不存在: ", manifest$relative_path[i]))
    ok(identical(tolower(sha(p)), tolower(manifest$output_sha256[i])) &&
         identical(bytes(p), as.numeric(manifest$file_size_bytes[i])),
       paste0("净化附件哈希或大小不一致: ", manifest$relative_path[i]))
    ok(identical(tolower(sha(file.path("results", manifest$source_basename[i]))),
                 tolower(manifest$source_sha256[i])),
       paste0("上游来源哈希已变化: ", manifest$source_basename[i]))
    observed <- fread(p, encoding = "UTF-8", showProgress = FALSE)
    compare_dt(cleaned[[i]]$data, observed, paste0("净化附件 ", basename(p)))
  }
}

write_stage <- function(cleaned, index, stage) {
  for (i in seq_len(nrow(source_map))) {
    p <- file.path(stage, OUT_DIR, source_map$output_basename[i])
    dir_create(path_dir(p), recurse = TRUE)
    fwrite(cleaned[[i]]$data, p, sep = "\t", quote = FALSE, na = "")
    compare_dt(cleaned[[i]]$data, fread(p, encoding = "UTF-8", showProgress = FALSE),
               paste0("暂存附件 ", basename(p)))
  }
  ip <- file.path(stage, OUT_INDEX)
  dir_create(path_dir(ip), recurse = TRUE)
  fwrite(index, ip, sep = "\t", quote = FALSE, na = "")
  compare_dt(index, fread(ip, encoding = "UTF-8"), "暂存投稿补充附件索引")
}

atomic_publish <- function(stage, paths) {
  tmp <- paste0(paths, ".script29_tmp_", Sys.getpid())
  on.exit({
    left <- tmp[file_exists(tmp)]
    if (length(left)) file_delete(left)
  }, add = TRUE)
  for (i in seq_along(paths)) {
    dir_create(path_dir(paths[i]), recurse = TRUE)
    file_copy(file.path(stage, paths[i]), tmp[i], overwrite = TRUE)
    ok(identical(sha(tmp[i]), sha(file.path(stage, paths[i]))),
       paste0("发布临时副本哈希不一致: ", paths[i]))
  }
  for (i in seq_along(paths)) {
    ok(file.rename(tmp[i], paths[i]), paste0("原子发布失败: ", paths[i]))
  }
}

validate_published <- function(cleaned, expected_index = NULL) {
  ok(file_exists(OUT_INDEX) && file_exists(OUT_MANIFEST),
     "投稿补充附件索引或内部 manifest 尚未发布。")
  index <- fread(OUT_INDEX, encoding = "UTF-8")
  manifest <- fread(OUT_MANIFEST, encoding = "UTF-8")
  validate_index(index)
  if (!is.null(expected_index)) compare_dt(expected_index, index, "已发布投稿补充附件索引")
  validate_manifest(manifest, cleaned, root)
  expected_files <- file.path(OUT_DIR, source_map$output_basename)
  actual_files <- dir_ls(OUT_DIR, regexp = "\\.tsv$", type = "file", fail = FALSE)
  actual_rel <- sub(paste0("^", normalizePath(root, winslash = "/"), "/"), "",
                    normalizePath(actual_files, winslash = "/", mustWork = TRUE))
  ok(setequal(expected_files, actual_rel) && length(actual_rel) == 71L,
     "投稿附件目录含缺失或未登记的 TSV。")
  invisible(list(index = index, manifest = manifest))
}

main <- function() {
  message("正在净化 71 个唯一补充附件来源……")
  source_hash_before <- setNames(vapply(file.path("results", source_map$source_basename), sha, character(1)),
                                 source_map$source_basename)
  cleaned <- build_cleaned()
  index <- make_index(cleaned)
  validate_index(index)

  if (validate_only) {
    published <- validate_published(cleaned, index)
    source_hash_after <- vapply(file.path("results", source_map$source_basename), sha, character(1))
    ok(identical(unname(source_hash_before), unname(source_hash_after)), "上游正式结果在校验期间发生变化。")
    message(sprintf(
      "校验完成：71 个附件，合计 %.2f MiB，删除 %d 列，翻译 %d 项，重命名 %d 列；未写入正式文件。",
      sum(published$manifest$file_size_bytes) / 1024^2,
      sum(published$manifest$removed_column_count),
      sum(vapply(strsplit(published$manifest$translated_columns, ";", fixed = TRUE),
                 function(x) sum(x != "none"), integer(1))),
      sum(vapply(strsplit(published$manifest$renamed_columns, ";", fixed = TRUE),
                 function(x) sum(x != "none"), integer(1)))
    ))
    return(invisible(TRUE))
  }

  stage <- file.path("_work", "intermediate", paste0("script29_stage_", Sys.getpid()))
  if (dir_exists(stage)) dir_delete(stage)
  dir_create(stage, recurse = TRUE)
  on.exit(if (dir_exists(stage)) dir_delete(stage), add = TRUE)
  write_stage(cleaned, index, stage)
  manifest <- make_manifest(cleaned, stage, output_is_stage = TRUE)
  mp <- file.path(stage, OUT_MANIFEST)
  dir_create(path_dir(mp), recurse = TRUE)
  fwrite(manifest, mp, sep = "\t", quote = FALSE, na = "")
  validate_manifest(fread(mp, encoding = "UTF-8"), cleaned, stage)

  expected_paths <- c(
    file.path(OUT_DIR, source_map$output_basename),
    OUT_INDEX,
    OUT_MANIFEST
  )
  if (dir_exists(OUT_DIR)) {
    extras <- setdiff(
      path_file(dir_ls(OUT_DIR, regexp = "\\.tsv$", type = "file", fail = FALSE)),
      source_map$output_basename
    )
    ok(!length(extras), paste0(
      "投稿附件目录存在未登记 TSV，拒绝静默删除: ",
      paste(extras, collapse = ", ")
    ))
  }
  atomic_publish(stage, expected_paths)
  published <- validate_published(cleaned, index)
  source_hash_after <- vapply(file.path("results", source_map$source_basename), sha, character(1))
  ok(identical(unname(source_hash_before), unname(source_hash_after)), "上游正式结果在发布期间发生变化。")
  message(sprintf(
    "完成：71 个英文净化附件，合计 %.2f MiB，删除 %d 列，翻译 %d 项，重命名 %d 列；上游表未改写。",
    sum(published$manifest$file_size_bytes) / 1024^2,
    sum(published$manifest$removed_column_count),
    sum(vapply(strsplit(published$manifest$translated_columns, ";", fixed = TRUE),
               function(x) sum(x != "none"), integer(1))),
    sum(vapply(strsplit(published$manifest$renamed_columns, ";", fixed = TRUE),
               function(x) sum(x != "none"), integer(1)))
  ))
}

main()
