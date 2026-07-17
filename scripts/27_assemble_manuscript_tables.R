#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(fs)
})

# 只从正式 results 组装投稿表格；不复制或改写原始数据。
OUT <- c(
  table1 = "results/escc_manuscript_table1_cohorts.tsv",
  table2 = "results/escc_manuscript_table2_event_state_associations.tsv",
  table3 = "results/escc_manuscript_table3_external_validation.tsv",
  supp_index = "results/escc_manuscript_supplementary_table_index.tsv",
  manifest = "results/escc_manuscript_table_artifact_manifest.tsv",
  main_md = "manuscript/escc_multiomics_main_tables.md",
  supp_md = "manuscript/escc_multiomics_supplementary_tables.md"
)
CLEAN_SUPP_INDEX <- "results/escc_submission_supplement_index.tsv"
CLEAN_SUPP_MANIFEST <- "results/escc_submission_supplement_artifact_manifest.tsv"
CLEAN_SUPP_DIR <- "results/submission_supplement"
EXCLUDED_INTERNAL_SUPP <- c(
  "escc_external_validation_decision.tsv",
  "tcga_escc_ecms_projection_qa.tsv"
)

args <- commandArgs(trailingOnly = TRUE)
bad <- setdiff(args, c("--fields-only", "--validate-only"))
if (length(bad)) stop("未知参数: ", paste(bad, collapse = ", "), call. = FALSE)
if ("--fields-only" %in% args) {
  cat(paste(unname(OUT), collapse = "\n"), "\n", sep = "")
  quit(save = "no", status = 0)
}
validate_only <- "--validate-only" %in% args

full_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", full_args, value = TRUE)
script_abs <- if (length(script_arg)) {
  normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath(file.path(getwd(), "scripts", "27_assemble_manuscript_tables.R"), mustWork = TRUE)
}
root <- normalizePath(file.path(dirname(script_abs), ".."), mustWork = TRUE)
setwd(root)
script_rel <- "scripts/27_assemble_manuscript_tables.R"
script_sha <- digest(file = script_rel, algo = "sha256", serialize = FALSE)
bt <- intToUtf8(96)

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
safe_rel <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) &&
    grepl("^(results|manuscript)/", x) && !grepl("(^|/)\\.\\.(/|$)", x)
}

manifest_paths <- c(
  multiassay = "results/tcga_escc_multiassay_artifact_manifest.tsv",
  driver_core = "results/tcga_escc_driver_core_artifact_manifest.tsv",
  driver_state = "results/tcga_escc_driver_state_artifact_manifest.tsv",
  leakage = "results/tcga_escc_driver_state_leakage_artifact_manifest.tsv",
  heterogeneity = "results/tcga_escc_heterogeneity_artifact_manifest.tsv",
  ecms = "results/tcga_escc_ecms_projection_artifact_manifest.tsv",
  cao = "results/cao2020_cross_layer_artifact_manifest.tsv",
  metabolomics = "results/pr001876_targeted_ms_artifact_manifest.tsv",
  microbiome = "results/prjna766558_dada2_artifact_manifest.tsv",
  integration = "results/escc_multiomics_integration_artifact_manifest.tsv",
  external = "results/escc_external_validation_artifact_manifest.tsv"
)
load_manifest <- function(path) {
  ok(file_exists(path), paste0("正式 manifest 不存在: ", path))
  x <- fread(path)
  need(x, c("relative_path", "file_size_bytes", "sha256", "generation_script", "status"), path)
  ok(nrow(x) > 0L, paste0("正式 manifest 为空: ", path))
  ok(!anyDuplicated(x$relative_path), paste0("manifest 路径重复: ", path))
  ok(all(vapply(x$relative_path, safe_rel, logical(1))), paste0("manifest 含不安全路径: ", path))
  ok(all(grepl("^verified($|_)", x$status)),
     paste0("manifest 含非 verified 状态族条目: ", path))
  x[, manifest_file := path]
  x
}
manifests <- lapply(manifest_paths, load_manifest)
registry <- rbindlist(lapply(manifests, function(x) x[, .(
  relative_path,
  manifest_file,
  declared_bytes = as.numeric(file_size_bytes),
  declared_sha = tolower(as.character(sha256)),
  generation_script
)]))
ok(!anyDuplicated(registry$relative_path), "多个 manifest 重复声明同一正式源文件。")
manifest_sha <- setNames(vapply(unname(manifest_paths), sha, character(1)), unname(manifest_paths))
verified <- new.env(parent = emptyenv())

verify_source <- function(path) {
  if (exists(path, envir = verified, inherits = FALSE)) {
    return(get(path, envir = verified, inherits = FALSE))
  }
  hit <- registry[relative_path == path]
  ok(nrow(hit) == 1L, paste0("正式源文件未被唯一 manifest 覆盖: ", path))
  ok(file_exists(path), paste0("正式源文件不存在: ", path))
  ok(identical(bytes(path), hit$declared_bytes[[1]]), paste0("源文件大小不符: ", path))
  ok(identical(tolower(sha(path)), hit$declared_sha[[1]]), paste0("源文件 SHA256 不符: ", path))
  assign(path, hit, envir = verified)
  hit
}
verify_many <- function(paths) invisible(lapply(unique(paths), verify_source))
source_sha <- function(path) verify_source(path)$declared_sha[[1]]
hash_bundle <- function(paths) {
  paths <- unique(paths)
  verify_many(paths)
  paste0(paths, "=", vapply(paths, source_sha, character(1)), collapse = ";")
}
source_manifests <- function(paths) {
  sort(unique(vapply(unique(paths), function(p) verify_source(p)$manifest_file[[1]], character(1))))
}
manifest_bundle <- function(paths) {
  m <- source_manifests(paths)
  paste0(m, "=", unname(manifest_sha[m]), collapse = ";")
}
formal <- function(path, cols = character()) {
  verify_source(path)
  x <- fread(path)
  ok(nrow(x) > 0L, paste0("正式源表为空: ", path))
  if (length(cols)) need(x, cols, path)
  x
}

# Table 3 的不可放宽门禁。
external_gate <- function() {
  x <- manifests$external
  need(x, c("proxy_model_sha256", "model_definition_sha256",
            "outcome_parsed_after_model_freeze"), manifest_paths[["external"]])
  ok(all(x$generation_script == "scripts/25_validate_external_continuous_states.R"),
     "外部验证 manifest 不是由正式 script25 生成。")
  ok(all(x$outcome_parsed_after_model_freeze %in% TRUE),
     "外部验证未通过模型冻结后才解析结局的门禁。")
  p <- unique(x$proxy_model_sha256[!is.na(x$proxy_model_sha256) & nzchar(x$proxy_model_sha256)])
  d <- unique(x$model_definition_sha256[!is.na(x$model_definition_sha256) & nzchar(x$model_definition_sha256)])
  ok(length(p) == 1L && length(d) == 1L, "外部代理模型或定义 SHA256 不唯一。")
  verify_many(x$relative_path)
  dec <- formal("results/escc_external_validation_decision.tsv",
                c("decision_id", "decision_domain", "status"))
  row <- dec[decision_id == "OVERALL_FIRST_STAGE_EXTERNAL_VALIDATION" &
             decision_domain == "overall"]
  ok(nrow(row) == 1L, "外部验证总体决策行缺失或重复。")
  ok(row$status[[1]] == "GO_FIRST_STAGE_EXTERNAL_RNA_PROXY_COMPLETED",
     paste0("外部验证未达到正式 Table 3 门禁: ", row$status[[1]]))
  row$status[[1]]
}

t1_sources <- c(
  "results/tcga_escc_multiassay_analysis_sets.tsv",
  "results/tcga_escc_ecms_patient_probabilities.tsv",
  "results/cao2020_cross_layer_candidate_summary.tsv",
  "results/pr001876_targeted_ms_analysis_inventory.tsv",
  "results/prjna766558_dada2_sample_qc.tsv",
  "results/escc_external_state_patient_scores.tsv",
  "results/escc_multiomics_module_summaries.tsv"
)
t2_sources <- "results/tcga_escc_driver_state_leakage_candidate_summary.tsv"
t3_sources <- c(
  "results/escc_external_state_internal_cv.tsv",
  "results/escc_external_state_survival_associations.tsv",
  "results/escc_external_state_ecms_increment.tsv",
  "results/escc_external_state_response_associations.tsv",
  "results/escc_external_validation_decision.tsv"
)

supp <- list(
  list("S1", "Cohort inventory, analysis units, and independence families",
       "Cohort/resource and patient-level analysis sets",
       "Defines cohort composition, analysis roles, sample subsets, and non-independence boundaries.",
       "The cohort inventory alone is not biological validation; overlapping TCGA subsets and GSE53625 subseries are not independent.",
       c("results/escc_multiomics_module_summaries.tsv",
         "results/tcga_escc_multiassay_analysis_sets.tsv",
         "results/tcga_escc_ecms_patient_probabilities.tsv",
         "results/escc_external_validation_decision.tsv")),
  list("S2", "TCGA-ESCC multi-assay sample map and analysis sets",
       "TCGA-ESCC patient", "Documents assay availability and locked patient subsets.",
       "All subsets belong to one overlapping TCGA discovery family.",
       c("results/tcga_escc_multiassay_sample_map.tsv",
         "results/tcga_escc_multiassay_analysis_sets.tsv",
         "results/tcga_escc_multiassay_qc.tsv")),
  list("S3", "TCGA-ESCC driver candidate screen and patient-level event matrix",
       "Gene/event and TCGA-ESCC patient", "Reports the driver screen, patient event calls, and quality controls.",
       "Driver priority is hypothesis-generating, not causal or therapeutically actionable.",
       c("results/tcga_escc_driver_candidate_screen.tsv",
         "results/tcga_escc_strong_driver_patient_events.tsv",
         "results/tcga_escc_driver_core_qc.tsv")),
  list("S4", "Somatic mutation and copy-number association details",
       "Gene, event, patient, and event pair", "Reports mutation recurrence, copy-number dosage, clusters, and event-pair summaries.",
       "Correlated encodings from the same TCGA patients are not separate validation families.",
       c("results/tcga_escc_cnv_expression_dosage.tsv",
         "results/tcga_escc_cnv_gene_summary.tsv",
         "results/tcga_escc_mutation_gene_summary.tsv",
         "results/tcga_escc_mutation_cluster_summary.tsv",
         "results/tcga_escc_mutation_pairwise_interactions.tsv")),
  list("S5", "MOFA factor, pathway, and clinical annotation results",
       "Patient, factor, feature, pathway, and clinical variable", "Reports factor variance, scores, weights, pathway annotation, quality controls, and clinical associations.",
       "Latent factors and annotations do not by themselves prove mechanism.",
       c("results/tcga_escc_mofa_variance_explained.tsv",
         "results/tcga_escc_mofa_factor_scores.tsv",
         "results/tcga_escc_mofa_top_weights.tsv",
         "results/tcga_escc_mofa_pathway_associations.tsv",
         "results/tcga_escc_mofa_level_factor_qc.tsv",
         "results/tcga_escc_mofa_clinical_associations.tsv",
         "results/tcga_escc_progeny_pathway_scores.tsv")),
  list("S6", "Patient-state heterogeneity clustering results",
       "Patient and candidate cluster solution", "Reports cluster evaluation, pathway contrasts, and patient assignments.",
       "Cluster labels are descriptive summaries from the TCGA discovery cohort, not independently validated subtypes.",
       c("results/tcga_escc_heterogeneity_cluster_evaluation.tsv",
         "results/tcga_escc_heterogeneity_cluster_pathways.tsv",
         "results/tcga_escc_heterogeneity_patient_assignments.tsv")),
  list("S7", "Representation-overlap-audited driver-state association analysis",
       "TCGA-ESCC gene-event by latent-state association", "Reports original associations, network edges, factor alignment, stability, and retained associations from the pre-specified event-view-deletion analysis.",
       "Event-view deletion reduces representation circularity within TCGA but is not independent validation or causal proof.",
       c("results/tcga_escc_driver_state_associations.tsv",
         "results/tcga_escc_driver_state_network_edges.tsv",
         "results/tcga_escc_driver_state_leakage_factor_alignment.tsv",
         "results/tcga_escc_driver_state_leakage_associations.tsv",
         "results/tcga_escc_driver_state_leakage_edge_stability.tsv",
         "results/tcga_escc_driver_state_leakage_candidate_summary.tsv")),
  list("S8", "Locked ECMS projection, calibration, and factor associations",
       "Patient, ECMS probability, factor, and pathway", "Documents projection, calibration, quality checks, and continuous state associations.",
       "The 78-patient anchor is a TCGA subset; the 16 projected extensions are excluded from primary association inference.",
       c("results/tcga_escc_ecms_patient_probabilities.tsv",
         "results/tcga_escc_ecms_projection_calibration.tsv",
         "results/tcga_escc_ecms_projection_qa.tsv",
         "results/tcga_escc_ecms_factor_associations.tsv",
         "results/tcga_escc_ecms_adjusted_factor_progeny_associations.tsv")),
  list("S9", "Cao et al. paired cross-layer calibration",
       "Pre-specified candidate within paired tumor and non-tumor specimens", "Reports candidate-restricted WGBS, RNA, protein, and cross-layer calibration.",
       "Coverage is candidate- and layer-dependent; missing protein data are not negative evidence.",
       c("results/cao2020_wgbs_candidate_pair_effects.tsv",
         "results/cao2020_wgbs_candidate_region_summary.tsv",
         "results/cao2020_rna_candidate_pair_effects.tsv",
         "results/cao2020_rna_candidate_summary.tsv",
         "results/cao2020_proteomics_candidate_pair_effects.tsv",
         "results/cao2020_proteomics_candidate_summary.tsv",
         "results/cao2020_cross_layer_patient_effects.tsv",
         "results/cao2020_cross_layer_candidate_summary.tsv")),
  list("S10", "PR001876 targeted metabolomics calibration",
       "Metabolite feature within analysis/platform", "Reports inventory, quality control, differential results, candidates, and sensitivities.",
       "The design is unpaired; overlap is unresolved and absolute concentration claims are not allowed.",
       c("results/pr001876_targeted_ms_analysis_inventory.tsv",
         "results/pr001876_targeted_ms_sample_qc.tsv",
         "results/pr001876_targeted_ms_feature_qc.tsv",
         "results/pr001876_targeted_ms_differential.tsv",
         "results/pr001876_targeted_ms_candidate_metabolites.tsv",
         "results/pr001876_targeted_ms_run_order_sensitivity.tsv",
         "results/pr001876_targeted_ms_scale_sensitivity.tsv",
         "results/pr001876_targeted_ms_paired_sensitivity.tsv")),
  list("S11", "PRJNA766558 paired 16S rRNA microbiome analysis",
       "ASV, genus, and paired FFPE specimen", "Reports denoising, taxonomy, diversity, paired tests, and pipeline sensitivity.",
       "No public negative control was identified; results are compositional and FFPE-specific.",
       c("results/prjna766558_dada2_sample_qc.tsv",
         "results/prjna766558_dada2_asv_counts.tsv",
         "results/prjna766558_dada2_taxonomy.tsv",
         "results/prjna766558_dada2_asv_length_qc.tsv",
         "results/prjna766558_dada2_alpha_diversity.tsv",
         "results/prjna766558_dada2_alpha_paired_tests.tsv",
         "results/prjna766558_dada2_beta_tests.tsv",
         "results/prjna766558_dada2_genus_paired_differential.tsv",
         "results/prjna766558_dada2_ancombc2_sensitivity.tsv",
         "results/prjna766558_dada2_pipeline_sensitivity.tsv")),
  list("S12", "Integrated evidence summary, axes, and candidate network",
       "Evidence item, candidate, edge, and integrated axis", "Reports evidence tiers, the integrated evidence table, candidates, edges, axes, and module synthesis.",
       "Integration summarizes public-data convergence; it does not convert association into mechanism.",
       c("results/escc_multiomics_integration_tier_definitions.tsv",
         "results/escc_multiomics_evidence_ledger.tsv",
         "results/escc_multiomics_integrated_driver_candidates.tsv",
         "results/escc_multiomics_integrated_axis_edges.tsv",
         "results/escc_multiomics_integrated_axis_summary.tsv",
         "results/escc_multiomics_heterogeneity_axes.tsv",
         "results/escc_multiomics_module_summaries.tsv")),
  list("S13", "External RNA-proxy and endpoint transportability calibration",
       "Patient, repeated CV run, factor, cohort, and endpoint model", "Reports frozen proxy definitions, internal fidelity, patient-level scores, survival models, ECMS incremental analyses, and exploratory pathological-response associations.",
       "Internal proxy fidelity does not constitute external biological replication; GSE53622/24 are one family, ECMS shares the 314-gene representation, and response analyses are exploratory.",
       c("results/escc_external_state_proxy_definition.tsv",
         "results/escc_external_state_internal_cv.tsv",
         "results/escc_external_state_oof_predictions.tsv",
         "results/escc_external_state_patient_scores.tsv",
         "results/escc_external_state_survival_associations.tsv",
         "results/escc_external_state_ecms_increment.tsv",
         "results/escc_external_state_response_associations.tsv",
         "results/escc_external_validation_decision.tsv"))
)

t1row <- function(order, cohort, design, samples, layers, role, boundary, sources) {
  data.table(
    row_order = as.integer(order), cohort_or_resource = cohort,
    study_design_and_analysis_unit = design, sample_accounting = samples,
    omics_or_endpoints = layers, primary_manuscript_role = role,
    evidence_boundary = boundary
  )
}

build_table1 <- function() {
  sets <- formal("results/tcga_escc_multiassay_analysis_sets.tsv",
                 c("analysis_set", "patient_id", "included"))
  got <- sets[included %in% TRUE, .N, by = analysis_set]
  got <- setNames(got$N, got$analysis_set)
  exp <- c(clinical_whitelist = 96L, driver_core = 95L,
           five_layer_core = 94L, protein_deep_subset = 76L)
  ok(identical(as.integer(got[names(exp)]), as.integer(exp)),
     "TCGA 分析集不符合锁定的 96/95/94/76 结构。")

  ecms <- formal("results/tcga_escc_ecms_patient_probabilities.tsv",
                 c("in_official_78", "resolved_ecms_label", "eligible_for_primary_association"))
  ok(nrow(ecms) == 94L && sum(ecms$in_official_78 %in% TRUE) == 78L,
     "ECMS 94/78 样本结构不符。")
  ec <- ecms[in_official_78 %in% TRUE, .N, by = resolved_ecms_label]
  ec <- setNames(ec$N, ec$resolved_ecms_label)
  ok(identical(as.integer(ec[c("ECMS1", "ECMS2", "ECMS3", "ECMS4")]),
               c(23L, 34L, 7L, 14L)), "ECMS1-4 锁定计数不符。")

  cao <- formal("results/cao2020_cross_layer_candidate_summary.tsv",
                c("promoter_total_pairs", "paired_patients", "quantified_patients"))
  ok(nrow(cao) == 12L && max(cao$promoter_total_pairs) == 9L &&
       max(cao$paired_patients) == 10L &&
       identical(as.integer(range(cao$quantified_patients)), c(0L, 7L)),
     "Cao 候选或配对覆盖不符合正式结构。")
  met <- formal("results/pr001876_targeted_ms_analysis_inventory.tsv",
                c("analysis_id", "n_escc", "n_normal"))
  ok(nrow(met) == 3L && all(met$n_escc == 16L & met$n_normal == 16L),
     "PR001876 样本结构不符。")
  mic <- formal("results/prjna766558_dada2_sample_qc.tsv",
                c("patient_pair_id", "run_accession", "tissue_role"))
  ok(nrow(mic) == 42L && uniqueN(mic$patient_pair_id) == 21L,
     "PRJNA766558 42/21 样本配对结构不符。")
  ext <- formal("results/escc_external_state_patient_scores.tsv",
                c("dataset", "patient_id", "independence_family"))
  ex <- ext[, .N, by = dataset]
  ex <- setNames(ex$N, ex$dataset)
  ok(identical(as.integer(ex[c("GSE53622", "GSE53624", "GSE45670")]),
               c(60L, 119L, 28L)), "外部队列 60/119/28 结构不符。")
  formal("results/escc_multiomics_module_summaries.tsv", "module_id")

  rbindlist(list(
    t1row(1, "TCGA-ESCC multi-assay discovery cohort",
          "Primary ESCC tumors; patient-level matched multi-assay analysis",
          "96 clinical cohort; 95 driver-core; 94 five-layer core; 76 protein-deep",
          "Mutation, copy number, RNA, promoter methylation, and protein/RPPA",
          "Driver discovery, latent-state modeling, and heterogeneity mapping",
          "All TCGA-derived subsets belong to one overlapping discovery family.",
          c("results/tcga_escc_multiassay_analysis_sets.tsv",
            "results/escc_multiomics_module_summaries.tsv")),
    t1row(2, "Locked TCGA-ESCC ECMS anchor",
          "Locked classifier projection; patient-level probabilities and labels",
          "78 primary anchors: ECMS1 23, ECMS2 34, ECMS3 7, ECMS4 14; 16 extensions excluded from primary association",
          "RNA-derived ECMS probabilities with factor and pathway associations",
          "Continuous-state calibration and ECMS context",
          "This is a TCGA subset, not external validation; the 16 projected extensions are excluded from primary association inference.",
          "results/tcga_escc_ecms_patient_probabilities.tsv"),
    t1row(3, "Cao et al. paired ESCC cross-layer cohort",
          "Candidate-restricted paired tumor/non-tumor analysis",
          "12 candidates; 9 WGBS pairs, 10 RNA pairs, and candidate-dependent protein coverage of 0-7 patients",
          "WGBS, RNA, and proteomics", "Same-patient cross-layer calibration",
          "Coverage is layer-dependent; missing protein measurements are not negative evidence or genome-wide replication.",
          "results/cao2020_cross_layer_candidate_summary.tsv"),
    t1row(4, "PR001876 targeted metabolomics resource",
          "Unpaired early-stage ESCC versus normal tissue; feature-level analysis",
          "Three analyses; each includes 16 ESCC and 16 normal specimens",
          "Targeted GC-MS and positive/negative-mode LC-MS",
          "Orthogonal metabolic-module calibration",
          "Subject overlap is unresolved; LC modes share a family; absolute concentration claims are not allowed.",
          "results/pr001876_targeted_ms_analysis_inventory.tsv"),
    t1row(5, "PRJNA766558 paired FFPE microbiome resource",
          "Paired tumor/non-tumor 16S rRNA analysis; patient pair is the inferential unit",
          "42 specimens forming 21 patient pairs",
          "ASVs, taxonomy, diversity, and paired genus-level tests",
          "Orthogonal microbiome-module calibration",
          "No public negative control was identified; results are compositional and not patient-linked to TCGA.",
          "results/prjna766558_dada2_sample_qc.tsv"),
    t1row(6, "External expression cohorts",
          "Frozen continuous RNA-proxy projection with patient-level endpoints",
          "GSE53622: 60; GSE53624: 119; GSE45670: 28",
          "Bulk RNA proxy scores; overall survival; pathological response",
          "External transportability calibration of Factor1 and Factor3 RNA proxies",
          "GSE53622/24 are one GSE53625 family; GSE45670 is exploratory; neither reproduces event edges or full MOFA.",
          "results/escc_external_state_patient_scores.tsv")
  ))
}

build_table2 <- function() {
  p <- t2_sources[[1]]
  x <- formal(p, c("original_edge_id", "gene_id", "gene_name", "event_type",
                   "reference_factor", "original_association_status",
                   "original_effect", "original_q_value", "planned_seed_count",
                   "drop_both_gate_retention_rate", "countable_for_T3_T4",
                   "post_audit_maximum_status", "leakage_control_decision",
                   "exact_independent_validation", "p_q_scope"))
  x <- x[countable_for_T3_T4 %in% TRUE]
  ok(nrow(x) == 9L, "Table 2 必须完整保留 9 条 countable leakage-audit 边。")
  ok(all(x$exact_independent_validation %in% FALSE), "Table 2 不得标记独立验证。")
  ok(!anyDuplicated(x$original_edge_id), "Table 2 正式源边标识重复。")
  setorder(x, original_q_value, gene_name, reference_factor)
  x[, .(
    row_order = seq_len(.N), gene = gene_name, gene_id,
    event = fcase(event_type == "mutation", "Mutation",
                  event_type == "amplification", "Amplification",
                  event_type == "relative_cnv", "Relative copy-number state",
                  default = event_type),
    factor = reference_factor,
    effect_measure = fifelse(event_type == "relative_cnv", "Spearman rho",
                             "Rank-biserial correlation"),
    effect_estimate = as.numeric(original_effect),
    q_value = as.numeric(original_q_value),
    sensitivity_seed_count = as.integer(planned_seed_count),
    mutation_cnv_omission_retention_rate = as.numeric(drop_both_gate_retention_rate),
    audit_status = "Retained in the pre-specified within-TCGA representation-overlap analysis",
    event_family_note = fifelse(gene_name == "PIK3CA",
      "Relative-CNV and amplification encodings are non-independent views of the same TCGA copy-number family.",
      "One TCGA-derived event-state association."),
    evidence_boundary = "Within-TCGA descriptive support only; not independent replication, causal proof, or therapeutic-target validation."
  )]
}

t3row <- function(factor, component, scope, family, unit, n, nevent,
                  definition, measure, estimate, lo, hi, interval,
                  p, q, pcontext, secondary, sest, slo, shi,
                  status, ceiling, source, key) {
  data.table(
    factor, calibration_component = component, cohort_or_scope = scope,
    independence_family = family, analysis_unit = unit,
    n_total = as.integer(n), n_events_or_responders = as.integer(nevent),
    event_or_response_definition = definition, effect_measure = measure,
    effect_estimate = as.numeric(estimate), interval_lower = as.numeric(lo),
    interval_upper = as.numeric(hi), interval_type = interval,
    p_value = as.numeric(p), q_value = as.numeric(q), p_value_context = pcontext,
    secondary_measure = secondary, secondary_estimate = as.numeric(sest),
    secondary_interval_lower = as.numeric(slo),
    secondary_interval_upper = as.numeric(shi), analysis_status = status,
    evidence_boundary = ceiling, source_file = source, source_row_key = key
  )
}

build_table3 <- function() {
  cvp <- t3_sources[[1]]
  sp <- t3_sources[[2]]
  ip <- t3_sources[[3]]
  rp <- t3_sources[[4]]
  cv <- formal(cvp, c("record_type", "model_variant", "factor", "n",
                      "spearman_rho", "pearson_r", "spearman_min",
                      "spearman_max", "gate_status"))
  cv <- cv[record_type == "outer_repeat_summary" &
           model_variant == "ridge314_primary" &
           factor %chin% c("Factor1", "Factor3")]
  cv[, factor_order := match(factor, c("Factor1", "Factor3"))]
  setorder(cv, factor_order)
  ok(nrow(cv) == 2L && all(cv$gate_status == "GO"), "ridge314 CV summary 未完整通过。")
  cvout <- rbindlist(lapply(seq_len(nrow(cv)), function(i) t3row(
    cv$factor[i], "Internal proxy fidelity", "TCGA-ESCC locked ECMS anchor",
    "TCGA_ESCC_internal_nested_CV", "Patient; 20 repeated outer 5-fold cross-validations",
    cv$n[i], NA_integer_, "Not applicable",
    "Median Spearman rho across repeated outer CV", cv$spearman_rho[i],
    cv$spearman_min[i], cv$spearman_max[i], "Observed minimum-maximum across 20 repeats",
    NA_real_, NA_real_, "Not applicable", "Pearson correlation",
    cv$pearson_r[i], NA_real_, NA_real_, cv$gate_status[i],
    "Internal out-of-fold fidelity of an outcome-blind RNA proxy; not external biological replication.",
    cvp, paste(cv$record_type[i], cv$model_variant[i], cv$factor[i], sep = "|")
  )))

  sv <- formal(sp, c("cohort_scope", "independence_family", "model_variant",
                     "factor", "model_type", "n_patients", "n_events",
                     "hazard_ratio_per_1sd", "ci_lower_95", "ci_upper_95",
                     "p_value", "q_value", "ph_score_p_value",
                     "optimized_cutpoint_used", "interpretation_status"))
  sv <- sv[model_variant == "ridge314_primary" & model_type == "clinical_adjusted" &
           factor %chin% c("Factor1", "Factor3")]
  sv[, fo := match(factor, c("Factor1", "Factor3"))]
  sv[, co := match(cohort_scope, c("GSE53622", "GSE53624", "GSE53622_GSE53624_stratified"))]
  setorder(sv, fo, co)
  ok(nrow(sv) == 6L && all(!is.na(sv$co)) &&
       all(sv$optimized_cutpoint_used %in% FALSE),
     "ridge314 临床校正生存模型不是锁定的 6 行或使用了优化切点。")
  svout <- rbindlist(lapply(seq_len(nrow(sv)), function(i) t3row(
    sv$factor[i], "External overall-survival association", sv$cohort_scope[i],
    sv$independence_family[i], "Patient; clinical-adjusted Cox model",
    sv$n_patients[i], sv$n_events[i], "Overall-survival events",
    "Adjusted hazard ratio per 1-SD proxy score", sv$hazard_ratio_per_1sd[i],
    sv$ci_lower_95[i], sv$ci_upper_95[i], "95% confidence interval",
    sv$p_value[i], sv$q_value[i], "Wald test for the adjusted score coefficient",
    "Proportional-hazards score-test p value", sv$ph_score_p_value[i],
    NA_real_, NA_real_, sv$interpretation_status[i],
    "External association within one GSE53625 family; null or heterogeneous results are retained and no optimized cutpoint was used.",
    sp, paste(sv$cohort_scope[i], sv$model_variant[i], sv$factor[i], sv$model_type[i], sep = "|")
  )))

  inc <- formal(ip, c("cohort_scope", "independence_family", "model_variant",
                      "factor", "n_patients", "n_events", "lrt_p_value",
                      "lrt_q_value", "score_hazard_ratio_per_1sd",
                      "score_ci_lower_95", "score_ci_upper_95",
                      "optimism_corrected_delta_cindex",
                      "optimism_corrected_delta_cindex_ci_lower_95",
                      "optimism_corrected_delta_cindex_ci_upper_95",
                      "increment_status"))
  inc <- inc[model_variant == "ridge314_primary" & factor %chin% c("Factor1", "Factor3")]
  inc[, factor_order := match(factor, c("Factor1", "Factor3"))]
  setorder(inc, factor_order)
  ok(nrow(inc) == 2L, "ridge314 ECMS 增量结果不是 2 行。")
  iout <- rbindlist(lapply(seq_len(nrow(inc)), function(i) t3row(
    inc$factor[i], "Conditional increment over ECMS", inc$cohort_scope[i],
    inc$independence_family[i], "Patient; ECMS-adjusted stratified Cox comparison",
    inc$n_patients[i], inc$n_events[i], "Overall-survival events",
    "Optimism-corrected delta concordance index",
    inc$optimism_corrected_delta_cindex[i],
    inc$optimism_corrected_delta_cindex_ci_lower_95[i],
    inc$optimism_corrected_delta_cindex_ci_upper_95[i],
    "95% paired-bootstrap optimism-corrected interval",
    inc$lrt_p_value[i], inc$lrt_q_value[i],
    "Likelihood-ratio test for adding the score to the ECMS model",
    "Adjusted score hazard ratio per 1 SD", inc$score_hazard_ratio_per_1sd[i],
    inc$score_ci_lower_95[i], inc$score_ci_upper_95[i], inc$increment_status[i],
    "Conditional increment over ECMS built from the shared 314-gene representation; not independent omics validation.",
    ip, paste(inc$cohort_scope[i], inc$model_variant[i], inc$factor[i], "ecms_increment", sep = "|")
  )))

  rr <- formal(rp, c("cohort", "model_variant", "factor", "n_total",
                     "n_pathological_complete_response",
                     "n_not_pathological_complete_response",
                     "rank_biserial_pcr_higher_positive",
                     "firth_odds_ratio_per_1sd", "firth_or_ci_lower_95",
                     "firth_or_ci_upper_95", "firth_profile_likelihood_p_value",
                     "firth_q_value", "treatment_predictor_claim", "response_status"))
  rr <- rr[model_variant == "ridge314_primary" & factor %chin% c("Factor1", "Factor3")]
  rr[, factor_order := match(factor, c("Factor1", "Factor3"))]
  setorder(rr, factor_order)
  ok(nrow(rr) == 2L && all(rr$treatment_predictor_claim %in% FALSE),
     "ridge314 疗效探索结果不是 2 行或边界被升级。")
  rout <- rbindlist(lapply(seq_len(nrow(rr)), function(i) t3row(
    rr$factor[i], "Exploratory pathological-response association", rr$cohort[i],
    "GSE45670", "Patient; Firth logistic regression", rr$n_total[i],
    rr$n_pathological_complete_response[i],
    paste0(rr$n_pathological_complete_response[i], " pathological complete responses; ",
           rr$n_not_pathological_complete_response[i], " non-complete responses"),
    "Firth odds ratio for pathological complete response per 1 SD",
    rr$firth_odds_ratio_per_1sd[i], rr$firth_or_ci_lower_95[i],
    rr$firth_or_ci_upper_95[i], "95% profile-likelihood confidence interval",
    rr$firth_profile_likelihood_p_value[i], rr$firth_q_value[i],
    "Firth profile-likelihood test for the score coefficient",
    "Rank-biserial correlation; positive means higher score in pCR",
    rr$rank_biserial_pcr_higher_positive[i], NA_real_, NA_real_,
    rr$response_status[i],
    "Exploratory association only; no treatment-predictor, causal, or clinical-utility claim.",
    rp, paste(rr$cohort[i], rr$model_variant[i], rr$factor[i], "response", sep = "|")
  )))
  out <- rbindlist(list(cvout, svout, iout, rout))
  out[, analysis_status := fcase(
    analysis_status == "GO", "Proxy fidelity criterion met",
    analysis_status == "secondary_null_or_inconclusive", "Null or inconclusive",
    analysis_status == "no_incremental_support", "No incremental support",
    analysis_status == "exploratory_concordant_signal", "Exploratory concordant signal",
    analysis_status == "exploratory_single_method_signal", "Exploratory single-method signal",
    default = gsub("_", " ", analysis_status)
  )]
  out[, row_order := seq_len(.N)]
  setcolorder(out, c("row_order", setdiff(names(out), "row_order")))
  ok(nrow(out) == 12L && !anyDuplicated(out$source_row_key),
     "Table 3 必须无显著性筛选地保留 12 行且 source_row_key 唯一。")
  out[, c("source_file", "source_row_key") := NULL]
  out
}

load_clean_supp <- function() {
  ok(file_exists(CLEAN_SUPP_INDEX),
     "投稿净化附件索引不存在；请先运行 scripts/29_prepare_submission_supplement.R。")
  ok(file_exists(CLEAN_SUPP_MANIFEST),
     "投稿净化附件内部 manifest 不存在；请先运行 script29。")
  detail <- fread(CLEAN_SUPP_INDEX, encoding = "UTF-8")
  manifest <- fread(CLEAN_SUPP_MANIFEST, encoding = "UTF-8")
  need(detail, c(
    "supplementary_table_id", "title", "clean_attachment_basename", "rows",
    "columns", "source_basename", "interpretation_boundary"
  ), CLEAN_SUPP_INDEX)
  need(manifest, c(
    "source_basename", "clean_attachment_basename", "relative_path",
    "source_sha256", "output_sha256", "file_size_bytes", "output_rows",
    "output_columns", "execution_script_sha256", "status"
  ), CLEAN_SUPP_MANIFEST)
  ok(nrow(detail) == 71L && nrow(manifest) == 71L &&
       !anyDuplicated(detail$source_basename) &&
       !anyDuplicated(detail$clean_attachment_basename) &&
       !anyDuplicated(manifest$relative_path),
     "投稿净化附件索引/manifest 必须唯一登记 71 个附件。")
  expected_sources <- all_supp_sources()
  ok(setequal(detail$source_basename, path_file(expected_sources)) &&
       !any(detail$source_basename %chin% EXCLUDED_INTERNAL_SUPP),
     "投稿净化附件来源不等于冻结的 71 个可投稿上游表。")
  ok(all(manifest$status == "verified_submission_clean"),
     "投稿净化附件存在非 verified_submission_clean 状态。")
  ok(file_exists("scripts/29_prepare_submission_supplement.R") &&
       all(manifest$execution_script_sha256 ==
             sha("scripts/29_prepare_submission_supplement.R")),
     "script29 已变化或其 manifest 未更新。")
  setkey(detail, source_basename, clean_attachment_basename)
  setkey(manifest, source_basename, clean_attachment_basename)
  ok(identical(detail[, .(source_basename, clean_attachment_basename)],
               manifest[, .(source_basename, clean_attachment_basename)]),
     "投稿净化索引与内部 manifest 的来源—附件映射不一致。")
  for (i in seq_len(nrow(manifest))) {
    p <- manifest$relative_path[i]
    source <- file.path("results", manifest$source_basename[i])
    ok(file_exists(p) && file_exists(source), paste0("净化附件或来源缺失: ", p))
    ok(identical(bytes(p), as.numeric(manifest$file_size_bytes[i])) &&
         identical(tolower(sha(p)), tolower(manifest$output_sha256[i])) &&
         identical(tolower(sha(source)), tolower(manifest$source_sha256[i])),
       paste0("净化附件或来源哈希不一致: ", p))
    z <- fread(p, encoding = "UTF-8", showProgress = FALSE)
    ok(nrow(z) == manifest$output_rows[i] && ncol(z) == manifest$output_columns[i] &&
         nrow(z) > 0L && ncol(z) > 0L,
       paste0("净化附件维度不一致或为空: ", p))
    vals <- unlist(z[, which(vapply(z, is.character, logical(1))), with = FALSE])
    ok(!has_han(c(names(z), vals)), paste0("净化附件仍含汉字: ", p))
    ok(!any(grepl(
      "/Users/|ResearchDataHub|results/|_work/|sha256|manifest|artifact|source_file|source_row_key|(^|_)relative_path$|target_relative_path|(^|_)run_status$|checkpoint|(^|_)gate(_|$)|_ceiling($|_)|required[ _]next[ _]validation|decision[ _]basis",
      c(names(z), vals), ignore.case = TRUE, perl = TRUE
    )), paste0("净化附件仍含内部路径、哈希或裁决词: ", p))
  }
  detail[, source_order__ := match(source_basename, path_file(expected_sources))]
  setorder(detail, source_order__)
  detail[, source_order__ := NULL]
  detail
}

build_supp_index <- function(clean_detail) {
  supported_claims <- c(
    S1 = "Defines cohort composition and analysis-unit boundaries for interpreting overlap and evidence independence.",
    S2 = "Defines the locked TCGA patient sets used by each analysis module.",
    S3 = "Supports prioritization of recurrent driver candidates for downstream state analyses.",
    S4 = "Supports descriptive mutation, copy-number, and dosage relationships within TCGA-ESCC.",
    S5 = "Supports latent-state annotation and pathway and clinical characterization.",
    S6 = "Supports descriptive patient-state heterogeneity within TCGA-ESCC.",
    S7 = "Supports within-TCGA event-state associations that persist in pre-specified mutation/CNV omission sensitivity analyses.",
    S8 = "Supports ECMS-contextualized continuous-state associations in the locked TCGA anchor.",
    S9 = "Supports candidate-restricted same-patient cross-layer calibration.",
    S10 = "Supports orthogonal metabolic-module calibration under the stated design limitations.",
    S11 = "Supports paired FFPE microbiome contrasts under the stated compositional and quality-control limitations.",
    S12 = "Supports integrated ranking of hypotheses and cross-layer axes.",
    S13 = "Supports RNA-proxy transportability assessment and external endpoint calibration, not full MOFA or event-level replication."
  )
  x <- rbindlist(lapply(supp, function(d) {
    sources <- setdiff(path_file(d[[6]]), EXCLUDED_INTERNAL_SUPP)
    belongs <- vapply(
      strsplit(clean_detail$supplementary_table_id, ";", fixed = TRUE),
      function(ids) d[[1]] %chin% ids,
      logical(1)
    )
    attachments <- clean_detail[belongs]
    ok(setequal(attachments$source_basename, sources),
       paste0(d[[1]], " 的净化附件来源与冻结清单不一致。"))
    attachments[, source_order__ := match(source_basename, sources)]
    setorder(attachments, source_order__)
    attachments[, source_order__ := NULL]
    data.table(
      supplementary_table_id = d[[1]],
      title = d[[2]],
      analysis_unit = d[[3]],
      content_description = d[[4]],
      supported_claim = unname(supported_claims[d[[1]]]),
      interpretation_boundary = d[[5]],
      attachment_basenames = paste(attachments$clean_attachment_basename, collapse = ";")
    )
  }))
  ok(identical(x$supplementary_table_id, paste0("S", 1:13)),
     "补充表目录必须完整覆盖 S1-S13。")
  listed <- unlist(strsplit(x$attachment_basenames, ";", fixed = TRUE))
  ok(uniqueN(listed) == 71L && length(listed) == 74L,
     "submission-facing 补充附件应为 74 次分组引用/71 个唯一净化表。")
  x
}

md_escape <- function(x) {
  x <- ifelse(is.na(x), "—", as.character(x))
  x <- gsub("[\r\n]+", " ", x)
  gsub("\\|", "\\\\|", x)
}
md_table <- function(x) {
  ok(nrow(x) > 0L && ncol(x) > 0L, "Markdown 表不能为空。")
  y <- as.data.table(lapply(x, md_escape))
  c(
    paste0("| ", paste(names(y), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(y)), collapse = " | "), " |"),
    apply(y, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  )
}
fmt_p <- function(x) vapply(x, function(z) {
  if (is.na(z)) return("—")
  if (z < 0.001) return(formatC(z, format = "e", digits = 2))
  formatC(z, format = "f", digits = 3)
}, character(1))
status_label <- function(x) {
  map <- c(
    GO = "Proxy fidelity criterion met",
    secondary_null_or_inconclusive = "Null or inconclusive",
    no_incremental_support = "No incremental support",
    exploratory_concordant_signal = "Exploratory concordant signal",
    exploratory_single_method_signal = "Exploratory single-method signal"
  )
  z <- unname(map[x])
  z[is.na(z)] <- gsub("_", " ", x[is.na(z)])
  z
}
effect_text <- function(row) {
  e <- row$effect_estimate
  lo <- row$interval_lower
  hi <- row$interval_upper
  if (row$calibration_component == "Internal proxy fidelity")
    return(sprintf("rho = %.3f (repeat range %.3f to %.3f)", e, lo, hi))
  if (row$calibration_component == "External overall-survival association")
    return(sprintf("HR = %.3f (95%% CI %.3f to %.3f)", e, lo, hi))
  if (row$calibration_component == "Conditional increment over ECMS")
    return(sprintf("corrected delta C = %.4f (95%% bootstrap interval %.4f to %.4f)", e, lo, hi))
  sprintf("OR = %.3f (95%% profile-likelihood CI %.3f to %.3f)", e, lo, hi)
}

build_main_md <- function(t1, t2, t3) {
  d1 <- t1[, .(
    cohort_or_resource,
    study_design_and_analysis_unit,
    sample_accounting,
    omics_or_endpoints,
    primary_manuscript_role,
    evidence_boundary
  )]
  setnames(d1, c(
    "Cohort or resource", "Design and analysis unit", "Sample accounting",
    "Omics or endpoints", "Primary role", "Independence and interpretation boundary"
  ))
  d2 <- t2[, .(
    gene, event, factor, effect_measure,
    effect = sprintf("%.3f", effect_estimate),
    q = fmt_p(q_value),
    retention = sprintf("%.1f%%", 100 * mutation_cnv_omission_retention_rate),
    interpretation = audit_status
  )]
  setnames(d2, c(
    "Gene", "Event", "Factor", "Effect measure", "Effect", "q",
    "Drop-both retention", "Interpretation"
  ))
  ntext <- vapply(seq_len(nrow(t3)), function(i) {
    r <- t3[i]
    if (r$calibration_component == "Internal proxy fidelity") return(as.character(r$n_total))
    if (r$calibration_component == "Exploratory pathological-response association")
      return(sprintf("%d (%d pCR)", r$n_total, r$n_events_or_responders))
    sprintf("%d (%d events)", r$n_total, r$n_events_or_responders)
  }, character(1))
  etext <- vapply(seq_len(nrow(t3)), function(i) effect_text(as.list(t3[i])), character(1))
  d3 <- data.table(
    factor = t3$factor, component = t3$calibration_component,
    scope = t3$cohort_or_scope, n = ntext, effect = etext,
    p = fmt_p(t3$p_value), q = fmt_p(t3$q_value),
    status = t3$analysis_status, boundary = t3$evidence_boundary
  )
  setnames(d3, c(
    "Factor", "Component", "Cohort or scope", "n (events/responders)",
    "Effect (interval)", "P", "q", "Status", "Interpretation boundary"
  ))
  c(
    "# Main Tables",
    "",
    "## Table 1. Cohort and resource inventory with analysis-unit and independence boundaries",
    "", md_table(d1), "",
    paste0(
      "*Footnotes.* Counts are analysis-specific and may overlap. TCGA-derived rows represent one ",
      "discovery family. GSE53622 and GSE53624 are subseries of GSE53625 and form one external survival ",
      "family. TCGA, The Cancer Genome Atlas; ESCC, esophageal squamous cell carcinoma; ECMS, ESCC ",
      "molecular subtype; WGBS, whole-genome bisulfite sequencing; RPPA, reverse-phase protein array; ",
      "FFPE, formalin-fixed paraffin-embedded."
    ),
    "",
    "## Table 2. Representation-overlap-audited TCGA event-state associations retained in the pre-specified analysis",
    "", md_table(d2), "",
    paste0(
      "*Footnotes.* Binary mutation or amplification effects are rank-biserial correlations; relative ",
      "copy-number effects are Spearman correlations. q values use the original within-TCGA multiplicity ",
      "families. Drop-both retention is the proportion of three pre-specified seeds retaining the association ",
      "after both the mutation and copy-number event views were omitted from latent-state refitting. ",
      "The two PIK3CA encodings are non-independent. ",
      "Every row remains within-TCGA descriptive support and is not independent replication, causal proof, ",
      "or therapeutic-target validation."
    ),
    "",
    "## Table 3. Frozen RNA-proxy fidelity and external endpoint transportability calibration",
    "", md_table(d3), "",
    paste0(
      "*Footnotes.* The ridge314 proxy was frozen without external outcome fitting. Internal repeated nested ",
      "cross-validation quantifies fidelity and is not external biological replication. Survival hazard ratios ",
      "are adjusted for age, sex, and stage; the pooled model is stratified by GSE53622 versus GSE53624. No ",
      "optimized cutpoint was used. For ECMS increment rows, the displayed effect is the paired-bootstrap ",
      "optimism-corrected concordance change, whereas P and q refer to the likelihood-ratio test. GSE45670 ",
      "response analyses are exploratory and do not support a treatment-predictor claim. Null, heterogeneous, ",
      "and borderline results are retained without significance filtering. Event-level edges and the full ",
      "multi-omics factor model were not externally reproduced. HR, hazard ratio; CI, confidence interval; ",
      "pCR, pathological complete response; OR, odds ratio."
    )
  )
}

build_supp_md <- function(index) {
  ok(identical(index$supplementary_table_id, paste0("S", 1:13)),
     "补充表展示面必须完整覆盖 S1-S13。")
  lines <- c(
    "# Supplementary Tables", "",
    paste0(
      "*Abbreviations.* TCGA, The Cancer Genome Atlas; ESCC, esophageal squamous cell carcinoma; ",
      "ECMS, ESCC molecular subtype; MOFA, multi-omics factor analysis."
    )
  )
  for (i in seq_along(supp)) {
    d <- supp[[i]]
    row <- index[supplementary_table_id == d[[1]]]
    ok(nrow(row) == 1L, paste0("补充表展示面缺少 ", d[[1]], "。"))
    attachments <- strsplit(row$attachment_basenames, ";", fixed = TRUE)[[1]]
    lines <- c(
      lines, "", paste0("## Supplementary Table ", d[[1]], ". ", d[[2]]), "",
      paste0("**Contents.** ", d[[4]]), "",
      paste0("**Analysis unit.** ", d[[3]]), "",
      paste0("**Interpretation note.** ", d[[5]]), "",
      paste0("**Attachments.** ", paste0("`", attachments, "`", collapse = "; "))
    )
  }
  unname(lines)
}

has_han <- function(x) any(grepl("[\u4e00-\u9fff]", x, perl = TRUE), na.rm = TRUE)
validate_md <- function(lines, label) {
  ok(length(lines) > 10L && nzchar(trimws(paste(lines, collapse = "\n"))),
     paste0(label, " 为空或过短。"))
  text <- paste(lines, collapse = "\n")
  ok(!grepl("TODO|TBD|PLACEHOLDER|pending|planned|待补|待定", text, ignore.case = TRUE),
     paste0(label, " 含占位或未完成标记。"))
  ok(!grepl("/Users/", text, fixed = TRUE), paste0(label, " 泄露绝对路径。"))
  internal_terms <- paste(
    c(
      "Artifact role", "manifest-verified", "authoritative", "results/",
      "SHA256", "Upload instruction", "machine-readable", "display rendering",
      "display/provenance", "Assembly status", "ready_as_",
      "external endpoint validation"
    ),
    collapse = "|"
  )
  ok(!grepl(internal_terms, text, ignore.case = TRUE),
     paste0(label, " 含内部追溯或旧版 external endpoint validation 用语。"))
  ok(!has_han(lines), paste0(label, " 正式表题或脚注含中文。"))
}
validate_bundle <- function(b) {
  ok(nrow(b$table1) == 6L && nrow(b$table2) == 9L &&
       nrow(b$table3) == 12L && nrow(b$supp_index) == 13L,
     "表格行数门禁失败。")
  expected_columns <- list(
    table1 = c(
      "row_order", "cohort_or_resource", "study_design_and_analysis_unit",
      "sample_accounting", "omics_or_endpoints", "primary_manuscript_role",
      "evidence_boundary"
    ),
    table2 = c(
      "row_order", "gene", "gene_id", "event", "factor", "effect_measure",
      "effect_estimate", "q_value", "sensitivity_seed_count",
      "mutation_cnv_omission_retention_rate", "audit_status",
      "event_family_note", "evidence_boundary"
    ),
    table3 = c(
      "row_order", "factor", "calibration_component", "cohort_or_scope",
      "independence_family", "analysis_unit", "n_total",
      "n_events_or_responders", "event_or_response_definition",
      "effect_measure", "effect_estimate", "interval_lower", "interval_upper",
      "interval_type", "p_value", "q_value", "p_value_context",
      "secondary_measure", "secondary_estimate", "secondary_interval_lower",
      "secondary_interval_upper", "analysis_status", "evidence_boundary"
    ),
    supp_index = c(
      "supplementary_table_id", "title", "analysis_unit",
      "content_description", "supported_claim", "interpretation_boundary",
      "attachment_basenames"
    )
  )
  for (nm in names(expected_columns)) {
    ok(identical(names(b[[nm]]), expected_columns[[nm]]),
       paste0(nm, " 含非 submission-facing 字段或字段顺序不符。"))
  }
  ok(all(grepl("not independent replication", b$table2$evidence_boundary, fixed = TRUE)),
     "Table 2 独立证据边界被意外升级。")
  ok(any(b$table3$analysis_status == "Null or inconclusive") &&
       any(b$table3$analysis_status == "No incremental support"),
     "Table 3 未保留正式阴性结果。")
  ok(!any(grepl("/", b$supp_index$attachment_basenames, fixed = TRUE)),
     "补充附件列必须只含 basename。")
  attachment_names <- unlist(strsplit(
    b$supp_index$attachment_basenames, ";", fixed = TRUE
  ))
  ok(length(attachment_names) == 74L && uniqueN(attachment_names) == 71L &&
       all(grepl("^S[0-9]{2}(?:-S[0-9]{2})*__.+\\.tsv$", attachment_names, perl = TRUE)),
     "补充附件必须为 74 次分组引用/71 个唯一且可追踪的净化 TSV。")
  ok(!any(attachment_names %chin% path_file(all_supp_sources())),
     "补充目录仍直接列出上游正式结果表，而非投稿净化附件。")
  ok(!any(grepl("leakage_model_plan", b$supp_index$attachment_basenames, fixed = TRUE)),
     "S7 仍包含内部 leakage model plan。")
  for (nm in c("table1", "table2", "table3", "supp_index")) {
    x <- b[[nm]]
    vals <- unlist(x[, which(vapply(x, is.character, logical(1))), with = FALSE])
    ok(!has_han(vals), paste0(nm, " 正式表内容含中文。"))
    ok(!any(grepl("results/|/Users/|SHA256|manifest-verified|Artifact role|Upload instruction",
                  vals, ignore.case = TRUE)),
       paste0(nm, " 含内部路径、哈希或追溯语言。"))
  }
  validate_md(b$main_md, "主表 Markdown")
  validate_md(b$supp_md, "补充表 Markdown")
}

compare_dt <- function(expected, observed, label) {
  ok(identical(names(expected), names(observed)), paste0(label, " 字段顺序不一致。"))
  ok(nrow(expected) == nrow(observed), paste0(label, " 行数不一致。"))
  z <- all.equal(as.data.frame(expected), as.data.frame(observed),
                 check.attributes = FALSE, tolerance = 1e-12)
  ok(isTRUE(z), paste0(label, " 内容不一致: ", paste(z, collapse = "; ")))
}
write_stage <- function(b, stage) {
  dir_create(path(stage, "results"), recurse = TRUE)
  dir_create(path(stage, "manuscript"), recurse = TRUE)
  tables <- c("table1", "table2", "table3", "supp_index")
  for (nm in tables) {
    p <- path(stage, OUT[[nm]])
    fwrite(b[[nm]], p, sep = "\t", quote = FALSE, na = "")
    compare_dt(b[[nm]], fread(p, na.strings = ""), paste0("暂存 ", OUT[[nm]]))
  }
  writeLines(b$main_md, path(stage, OUT[["main_md"]]), useBytes = TRUE)
  writeLines(b$supp_md, path(stage, OUT[["supp_md"]]), useBytes = TRUE)
  ok(identical(readLines(path(stage, OUT[["main_md"]]), warn = FALSE), b$main_md),
     "主表 Markdown 暂存回读不一致。")
  ok(identical(readLines(path(stage, OUT[["supp_md"]]), warn = FALSE), b$supp_md),
     "补充表 Markdown 暂存回读不一致。")
}

all_supp_sources <- function() setdiff(
  unique(unlist(lapply(supp, function(d) d[[6]]))),
  file.path("results", EXCLUDED_INTERNAL_SUPP)
)
specs <- function() {
  main_up <- unique(c(t1_sources, t2_sources, t3_sources))
  list(
    list(OUT[["table1"]], "main_table_tsv", "authoritative_machine_readable_main_table",
         t1_sources, t1_sources),
    list(OUT[["table2"]], "main_table_tsv", "authoritative_machine_readable_main_table",
         t2_sources, t2_sources),
    list(OUT[["table3"]], "main_table_tsv", "authoritative_machine_readable_main_table",
         t3_sources, t3_sources),
    list(OUT[["supp_index"]], "supplementary_catalog_tsv", "authoritative_source_attachment_catalog",
         all_supp_sources(), all_supp_sources()),
    list(OUT[["main_md"]], "manuscript_display_markdown", "display_rendering_not_source_data",
         unname(OUT[c("table1", "table2", "table3")]), main_up),
    list(OUT[["supp_md"]], "manuscript_display_markdown", "display_rendering_not_source_data",
         OUT[["supp_index"]], all_supp_sources())
  )
}
count_content <- function(b, rel) {
  if (rel == OUT[["table1"]]) return(nrow(b$table1))
  if (rel == OUT[["table2"]]) return(nrow(b$table2))
  if (rel == OUT[["table3"]]) return(nrow(b$table3))
  if (rel == OUT[["supp_index"]]) return(nrow(b$supp_index))
  if (rel == OUT[["main_md"]]) return(length(b$main_md))
  if (rel == OUT[["supp_md"]]) return(length(b$supp_md))
  fail("未知制品: ", rel)
}
hash_any <- function(paths, base) paste0(
  paths, "=", vapply(paths, function(p) {
    if (p %chin% unname(OUT)) sha(path(base, p)) else source_sha(p)
  }, character(1)), collapse = ";"
)
make_manifest <- function(b, base, ext_status, date) {
  rbindlist(lapply(specs(), function(s) {
    rel <- s[[1]]
    p <- path(base, rel)
    ok(file_exists(p), paste0("待登记制品不存在: ", rel))
    data.table(
      artifact = path_file(rel), artifact_class = s[[2]],
      publication_role = s[[3]], relative_path = rel,
      row_or_line_count = as.integer(count_content(b, rel)),
      file_size_bytes = bytes(p), sha256 = sha(p),
      source_files = paste(s[[4]], collapse = ";"),
      source_sha256_bundle = hash_any(s[[4]], base),
      upstream_manifest_sha256_bundle = manifest_bundle(s[[5]]),
      generation_script = script_rel, execution_script_sha256 = script_sha,
      generated_date = as.character(date), external_validation_status = ext_status,
      raw_data_copied = FALSE, status = "verified"
    )
  }))
}
validate_manifest <- function(actual, expected, label) {
  required <- c(
    "artifact", "artifact_class", "publication_role", "relative_path",
    "row_or_line_count", "file_size_bytes", "sha256", "source_files",
    "source_sha256_bundle", "upstream_manifest_sha256_bundle",
    "generation_script", "execution_script_sha256", "generated_date",
    "external_validation_status", "raw_data_copied", "status"
  )
  need(actual, required, label)
  ok(nrow(actual) == 6L && !anyDuplicated(actual$relative_path),
     paste0(label, " 应唯一登记 6 个制品。"))
  ok(all(actual$status == "verified") && all(actual$raw_data_copied %in% FALSE),
     paste0(label, " status 或 raw_data_copied 不合格。"))
  a <- copy(actual)
  e <- copy(expected)
  a[, generated_date := as.character(generated_date)]
  e[, generated_date := as.character(generated_date)]
  a[, file_size_bytes := as.numeric(file_size_bytes)]
  e[, file_size_bytes := as.numeric(file_size_bytes)]
  a[, row_or_line_count := as.integer(row_or_line_count)]
  e[, row_or_line_count := as.integer(row_or_line_count)]
  setcolorder(a, names(e))
  compare_dt(e, a, label)
}
validate_content <- function(b, base, label) {
  for (nm in c("table1", "table2", "table3", "supp_index"))
    compare_dt(b[[nm]], fread(path(base, OUT[[nm]]), na.strings = ""),
               paste0(label, " ", OUT[[nm]]))
  ok(identical(readLines(path(base, OUT[["main_md"]]), warn = FALSE), b$main_md),
     paste0(label, " 主表 Markdown 不一致。"))
  ok(identical(readLines(path(base, OUT[["supp_md"]]), warn = FALSE), b$supp_md),
     paste0(label, " 补充表 Markdown 不一致。"))
}
validate_published <- function(b, ext_status, allow_absent = TRUE) {
  contents <- unname(OUT[names(OUT) != "manifest"])
  if (!file_exists(OUT[["manifest"]])) {
    ok(!any(file_exists(contents)), "发现无 manifest 的不完整投稿表格 bundle。")
    if (allow_absent) {
      message("尚无已发布投稿表格 bundle；内存组装校验已完成。")
      return(invisible(FALSE))
    }
    fail("投稿表格 manifest 未发布。")
  }
  ok(all(file_exists(contents)), "已发布 bundle 缺少制品。")
  m <- fread(OUT[["manifest"]])
  dates <- unique(as.character(m$generated_date))
  ok(length(dates) == 1L && nzchar(dates), "manifest generated_date 不唯一或为空。")
  expected <- make_manifest(b, root, ext_status, dates[[1]])
  validate_manifest(m, expected, "已发布制品 manifest")
  validate_content(b, root, "已发布制品")
  invisible(TRUE)
}
publish <- function(stage) {
  order <- c(unname(OUT[names(OUT) != "manifest"]), OUT[["manifest"]])
  tmp <- paste0(order, ".script27_tmp_", Sys.getpid())
  on.exit({
    left <- tmp[file_exists(tmp)]
    if (length(left)) file_delete(left)
  }, add = TRUE)
  for (i in seq_along(order)) {
    dir_create(path_dir(tmp[[i]]), recurse = TRUE)
    file_copy(path(stage, order[[i]]), tmp[[i]], overwrite = TRUE)
    ok(identical(sha(tmp[[i]]), sha(path(stage, order[[i]]))),
       paste0("发布临时副本 SHA 不一致: ", order[[i]]))
  }
  for (i in seq_along(order))
    ok(file.rename(tmp[[i]], order[[i]]), paste0("原子发布失败: ", order[[i]]))
}

main <- function() {
  message("正在核验外部验证硬门禁和正式 manifest……")
  ext_status <- external_gate()
  verify_many(unique(c(t1_sources, t2_sources, t3_sources, all_supp_sources())))
  clean_supp <- load_clean_supp()
  message("正在组装主表 1-3 与补充表 S1-S13 目录……")
  b <- list()
  b$table1 <- build_table1()
  b$table2 <- build_table2()
  b$table3 <- build_table3()
  b$supp_index <- build_supp_index(clean_supp)
  b$main_md <- build_main_md(b$table1, b$table2, b$table3)
  b$supp_md <- build_supp_md(b$supp_index)
  validate_bundle(b)

  if (validate_only) {
    validate_published(b, ext_status, allow_absent = TRUE)
    message("校验完成：未写入任何文件。")
    return(invisible(TRUE))
  }

  stage <- path("_work", "intermediate", paste0("script27_stage_", Sys.getpid()))
  if (dir_exists(stage)) dir_delete(stage)
  dir_create(stage, recurse = TRUE)
  on.exit(if (dir_exists(stage)) dir_delete(stage), add = TRUE)
  message("正在暂存并回读核验 6 个内容制品……")
  write_stage(b, stage)
  sm <- make_manifest(b, stage, ext_status, Sys.Date())
  mp <- path(stage, OUT[["manifest"]])
  dir_create(path_dir(mp), recurse = TRUE)
  fwrite(sm, mp, sep = "\t", quote = FALSE, na = "")
  observed <- fread(mp)
  expected <- make_manifest(b, stage, ext_status,
                            unique(as.character(observed$generated_date))[[1]])
  validate_manifest(observed, expected, "暂存制品 manifest")
  message("全部暂存制品通过，正在原子发布；manifest 最后发布……")
  publish(stage)
  validate_published(b, ext_status, allow_absent = FALSE)
  message("完成：主表、补充目录、Markdown 与制品 manifest 均已通过 SHA256 校验。")
}

main()
