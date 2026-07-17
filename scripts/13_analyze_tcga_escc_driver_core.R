#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(fs)
  library(maftools)
  library(matrixStats)
  library(Matrix)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
})

options(stringsAsFactors = FALSE)

# TCGA-ESCC 驱动核心首轮分析：
# 1) 非同义突变复发；2) ASCAT2 CNV 复发；3) CNV-RNA 剂量耦联；
# 4) 患者级候选收敛与突变共现/互斥。
#
# 本脚本不把公共数据关联写成因果驱动，不用 76 例完整个案替代
# 95 例 driver core，也不因单一 q 值或方向差异自动淘汰候选。

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required_inputs <- file.path(project_root, c(
  "PROJECT_INDEX.md",
  "results/tcga_escc_multiassay_dr45.rds",
  "results/tcga_escc_masked_maf_long.rds",
  "results/tcga_escc_ascat2_segments_long.rds",
  "results/tcga_escc_multiassay_analysis_sets.tsv"
))
if (any(!file_exists(required_inputs))) {
  stop("缺少 TCGA driver core 输入：",
       paste(required_inputs[!file_exists(required_inputs)], collapse = ", "))
}

results_dir <- file.path(project_root, "results")
work_intermediate_dir <- file.path(project_root, "_work", "intermediate")
dir_create(work_intermediate_dir, recurse = TRUE)
stage_dir <- tempfile(pattern = ".tcga_driver_core_", tmpdir = work_intermediate_dir)
dir_create(stage_dir)
on.exit({
  if (dir_exists(stage_dir)) dir_delete(stage_dir)
}, add = TRUE)

fail_if <- function(condition, message) {
  if (isTRUE(condition)) stop(message, call. = FALSE)
}

strip_ensembl_version <- function(x) sub("[.][0-9]+$", "", x)

representative_symbol <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[[1]]
}

extract_aa_position <- function(x) {
  x <- as.character(x)
  hit <- regexpr("[0-9]+", x)
  out <- rep(NA_real_, length(x))
  keep <- !is.na(hit) & hit > 0L
  out[keep] <- as.numeric(regmatches(x, hit))
  out
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) max(x) else NA_real_
}

pairwise_row_spearman <- function(x, y) {
  fail_if(!identical(dim(x), dim(y)), "相关矩阵维度不一致")
  rx <- rowRanks(
    x,
    ties.method = "average",
    na.last = "keep",
    preserveShape = TRUE,
    useNames = FALSE
  )
  ry <- rowRanks(
    y,
    ties.method = "average",
    na.last = "keep",
    preserveShape = TRUE,
    useNames = FALSE
  )
  valid <- is.finite(rx) & is.finite(ry)
  n <- rowSums(valid)
  rx[!valid] <- 0
  ry[!valid] <- 0
  mean_x <- rowSums(rx) / n
  mean_y <- rowSums(ry) / n
  dx <- sweep(rx, 1L, mean_x, "-")
  dy <- sweep(ry, 1L, mean_y, "-")
  dx[!valid] <- 0
  dy[!valid] <- 0
  numerator <- rowSums(dx * dy)
  denominator <- sqrt(rowSums(dx^2) * rowSums(dy^2))
  rho <- numerator / denominator
  rho[n < 4L | denominator == 0] <- NA_real_
  rho <- pmax(-1, pmin(1, rho))
  statistic <- rho * sqrt((n - 2) / pmax(1e-12, 1 - rho^2))
  p_value <- 2 * pt(-abs(statistic), df = pmax(1, n - 2))
  p_value[is.na(rho)] <- NA_real_
  data.table(n = n, rho = rho, p_value = p_value)
}

weighted_ploidy_proxy <- function(copy_number, segment_length) {
  keep <- is.finite(copy_number) & is.finite(segment_length) & segment_length > 0
  if (!any(keep)) return(NA_real_)
  weightedMedian(copy_number[keep], w = segment_length[keep], na.rm = TRUE)
}

mae <- readRDS(file.path(results_dir, "tcga_escc_multiassay_dr45.rds"))
validObject(mae)
maf_long <- as.data.table(readRDS(file.path(results_dir, "tcga_escc_masked_maf_long.rds")))
segments <- as.data.table(readRDS(file.path(results_dir, "tcga_escc_ascat2_segments_long.rds")))
analysis_sets <- fread(
  file.path(results_dir, "tcga_escc_multiassay_analysis_sets.tsv"),
  showProgress = FALSE
)

expected_experiments <- c("RNA", "miRNA", "HM450", "Mutation", "CNV_gene", "RPPA")
fail_if(!identical(sort(names(experiments(mae))), sort(expected_experiments)),
        "MultiAssayExperiment 实验视图不完整")

driver_patients <- analysis_sets[
  analysis_set == "driver_core" & included == TRUE,
  patient_id
]
fail_if(length(driver_patients) != 95L || uniqueN(driver_patients) != 95L,
        "driver core 必须为 95 位唯一患者")

rna_se <- experiments(mae)[["RNA"]]
cnv_se <- experiments(mae)[["CNV_gene"]]
mutation_se <- experiments(mae)[["Mutation"]]

rna_counts <- assay(rna_se, "counts")
rna_logcpm <- assay(rna_se, "tmm_logcpm")
cnv_total <- assay(cnv_se, "copy_number")
mutation_binary <- assay(mutation_se, "binary")

fail_if(!all(driver_patients %chin% colnames(rna_logcpm)),
        "driver core 患者未全部进入 RNA")
fail_if(!all(driver_patients %chin% colnames(cnv_total)),
        "driver core 患者未全部进入 CNV")
fail_if(!all(driver_patients %chin% colnames(mutation_binary)),
        "driver core 患者未全部进入 Mutation")

rna_rowdata <- as.data.table(as.data.frame(rowData(rna_se)))
cnv_rowdata <- as.data.table(as.data.frame(rowData(cnv_se)))
fail_if(anyDuplicated(rna_rowdata$gene_id) > 0L ||
          anyDuplicated(cnv_rowdata$gene_id) > 0L,
        "RNA 或 CNV gene_id 重复")

common_gene_ids <- intersect(rownames(rna_logcpm), rownames(cnv_total))
fail_if(length(common_gene_ids) != 60623L,
        "RNA-CNV 共同基因宇宙必须为 60,623 个 versioned Ensembl ID")
rna_idx <- match(common_gene_ids, rownames(rna_logcpm))
cnv_idx <- match(common_gene_ids, rownames(cnv_total))
cnv_total_common <- cnv_total[cnv_idx, , drop = FALSE]
fail_if(
  !identical(rownames(cnv_total_common), common_gene_ids),
  "ASCAT2 gene CNV 未按 RNA-CNV 共同基因顺序重排"
)
gene_annotation <- merge(
  cnv_rowdata,
  rna_rowdata[, .(gene_id, gene_type)],
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)
gene_annotation <- gene_annotation[match(common_gene_ids, gene_id)]

message("[1/5] 计算 segment 长度加权 ploidy proxy 与 CNV 复发")
segments_autosome <- segments[Chromosome %chin% paste0("chr", 1:22)]
ploidy <- segments_autosome[, .(
  ploidy_proxy = weighted_ploidy_proxy(Copy_Number, segment_length),
  autosomal_segment_count = .N,
  autosomal_covered_bp = sum(segment_length, na.rm = TRUE)
), by = patient_id]
fail_if(anyNA(ploidy$ploidy_proxy) || any(ploidy$ploidy_proxy <= 0),
        "ploidy proxy 计算失败")

all_cnv_patients <- colnames(cnv_total_common)
ploidy_vector <- ploidy$ploidy_proxy[match(all_cnv_patients, ploidy$patient_id)]
fail_if(anyNA(ploidy_vector), "部分 CNV 患者缺少 ploidy proxy")
relative_cnv <- log2(
  (cnv_total_common + 0.5) /
    matrix(
      ploidy_vector + 0.5,
      nrow = nrow(cnv_total_common),
      ncol = ncol(cnv_total_common),
      byrow = TRUE
    )
)

cnv_nonmissing <- rowSums(is.finite(cnv_total_common))
gain_frequency <- rowMeans(relative_cnv >= 0.3, na.rm = TRUE)
loss_frequency <- rowMeans(relative_cnv <= -0.3, na.rm = TRUE)
amplification_frequency <- rowMeans(
  relative_cnv >= 0.8 & cnv_total_common >= 4,
  na.rm = TRUE
)
homozygous_deletion_frequency <- rowMeans(cnv_total_common == 0, na.rm = TRUE)
cnv_summary <- data.table(
  gene_id = common_gene_ids,
  gene_name = gene_annotation$gene_name,
  chromosome = gene_annotation$chromosome,
  start = gene_annotation$start,
  end = gene_annotation$end,
  gene_type = gene_annotation$gene_type,
  nonmissing_patients = cnv_nonmissing,
  median_copy_number = rowMedians(
    cnv_total_common,
    na.rm = TRUE,
    useNames = FALSE
  ),
  median_relative_log2 = rowMedians(relative_cnv, na.rm = TRUE, useNames = FALSE),
  relative_log2_iqr = rowIQRs(relative_cnv, na.rm = TRUE, useNames = FALSE),
  gain_frequency = gain_frequency,
  loss_frequency = loss_frequency,
  amplification_frequency = amplification_frequency,
  homozygous_deletion_frequency = homozygous_deletion_frequency
)
cnv_summary[, dominant_cnv_event := fifelse(
  amplification_frequency >= homozygous_deletion_frequency &
    amplification_frequency > 0,
  "amplification",
  fifelse(homozygous_deletion_frequency > 0, "homozygous_deletion", "none")
)]
cnv_summary[, recurrent_cnv_class := fifelse(
  amplification_frequency >= 0.10, "recurrent_amplification",
  fifelse(homozygous_deletion_frequency >= 0.10, "recurrent_homozygous_deletion",
          fifelse(gain_frequency >= 0.20, "recurrent_gain",
                  fifelse(loss_frequency >= 0.20, "recurrent_loss", "not_recurrent")))
)]

message("[2/5] 计算 CNV-RNA 剂量耦联及参数敏感性")
driver_rna <- rna_logcpm[rna_idx, driver_patients, drop = FALSE]
driver_counts <- rna_counts[rna_idx, driver_patients, drop = FALSE]
driver_cnv_absolute <- cnv_total_common[, driver_patients, drop = FALSE]
driver_relative_cnv <- relative_cnv[, driver_patients, drop = FALSE]

expressed_filter <- rowSums(driver_counts >= 10, na.rm = TRUE) >= 10L
protein_coding_filter <- gene_annotation$gene_type == "protein_coding"
cnv_coverage_filter <- rowMeans(is.finite(driver_cnv_absolute)) >= 0.90
cnv_variation_filter <- rowIQRs(
  driver_relative_cnv,
  na.rm = TRUE,
  useNames = FALSE
) > 0.05
test_eligible <- expressed_filter & protein_coding_filter &
  cnv_coverage_filter & cnv_variation_filter

relative_cor <- data.table(
  n = rep(NA_integer_, length(common_gene_ids)),
  rho = rep(NA_real_, length(common_gene_ids)),
  p_value = rep(NA_real_, length(common_gene_ids))
)
absolute_cor <- copy(relative_cor)
relative_cor[test_eligible] <- pairwise_row_spearman(
  driver_relative_cnv[test_eligible, , drop = FALSE],
  driver_rna[test_eligible, , drop = FALSE]
)
absolute_cor[test_eligible] <- pairwise_row_spearman(
  driver_cnv_absolute[test_eligible, , drop = FALSE],
  driver_rna[test_eligible, , drop = FALSE]
)
relative_q <- rep(NA_real_, length(common_gene_ids))
absolute_q <- rep(NA_real_, length(common_gene_ids))
relative_q[test_eligible] <- p.adjust(relative_cor$p_value[test_eligible], method = "BH")
absolute_q[test_eligible] <- p.adjust(absolute_cor$p_value[test_eligible], method = "BH")

dosage <- data.table(
  gene_id = common_gene_ids,
  gene_name = gene_annotation$gene_name,
  gene_type = gene_annotation$gene_type,
  expressed_filter = expressed_filter,
  cnv_coverage_filter = cnv_coverage_filter,
  cnv_variation_filter = cnv_variation_filter,
  test_eligible = test_eligible,
  n_patients = relative_cor$n,
  spearman_rho_relative_cnv = relative_cor$rho,
  p_relative_cnv = relative_cor$p_value,
  q_relative_cnv = relative_q,
  spearman_rho_absolute_cnv = absolute_cor$rho,
  p_absolute_cnv = absolute_cor$p_value,
  q_absolute_cnv = absolute_q
)
dosage[, direction_stable_between_cnv_scales :=
         sign(spearman_rho_relative_cnv) == sign(spearman_rho_absolute_cnv)]
dosage[, discordance_code := fifelse(
  !test_eligible, "not_comparable",
  fifelse(
    spearman_rho_relative_cnv >= 0.30 & q_relative_cnv <= 0.10 &
      spearman_rho_absolute_cnv > 0,
    "support",
    fifelse(
      spearman_rho_relative_cnv >= 0.30 & p_relative_cnv <= 0.05 &
        direction_stable_between_cnv_scales,
      "same_direction_nonsignificant",
      fifelse(
        spearman_rho_relative_cnv <= -0.30 & q_relative_cnv <= 0.10 &
          spearman_rho_absolute_cnv < 0,
        "comparable_precise_reverse",
        fifelse(
          !direction_stable_between_cnv_scales &
            abs(spearman_rho_relative_cnv) >= 0.20,
          "opposite_imprecise",
          "same_direction_nonsignificant"
        )
      )
    )
  )
)]
dosage[, dosage_evidence := fifelse(
  discordance_code == "support", "strong_positive_dosage",
  fifelse(
    spearman_rho_relative_cnv >= 0.30 & p_relative_cnv <= 0.05 &
      direction_stable_between_cnv_scales,
    "conditional_positive_dosage",
    fifelse(discordance_code == "comparable_precise_reverse",
            "precise_reverse_for_dosage_axis", "no_strong_dosage_evidence")
  )
)]

message("[3/5] 汇总非同义突变并用 maftools 复核")
fail_if(!all(c("Hugo_Symbol", "Variant_Classification", "Tumor_Sample_Barcode") %chin%
               names(maf_long)),
        "MAF 长表缺少 maftools 核心字段")
maf_object <- read.maf(
  maf = as.data.frame(maf_long),
  removeDuplicatedVariants = TRUE,
  useAll = TRUE,
  verbose = FALSE
)
maftools_gene_summary <- as.data.table(getGeneSummary(maf_object))
maftools_sample_summary <- as.data.table(getSampleSummary(maf_object))
oncodrive_result <- suppressWarnings(oncodrive(
  maf_object,
  AACol = "HGVSp_Short",
  minMut = 5,
  pvalMethod = "zscore"
))
oncodrive_result <- as.data.table(oncodrive_result)

lof_classes <- c(
  "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site",
  "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation"
)
maf_long[, observed_aa_position := extract_aa_position(HGVSp_Short)]

mutation_summary <- maf_long[
  nonsynonymous == TRUE & !is.na(gene_id),
  .(
    Hugo_Symbol = representative_symbol(Hugo_Symbol),
    nonsynonymous_variants = .N,
    mutated_patients = uniqueN(patient_id),
    likely_lof_variants = sum(Variant_Classification %chin% lof_classes),
    likely_lof_patients = uniqueN(patient_id[Variant_Classification %chin% lof_classes]),
    hotspot_variants = sum(hotspot == "Y", na.rm = TRUE),
    max_observed_aa_position = safe_max(observed_aa_position),
    variant_classes = paste(sort(unique(Variant_Classification)), collapse = ";"),
    callers = paste(sort(unique(callers[!is.na(callers)])), collapse = ";")
  ),
  by = gene_id
]
mutation_summary <- merge(
  mutation_summary,
  oncodrive_result[, .(
    Hugo_Symbol,
    oncodrive_mutated_samples = MutatedSamples,
    oncodrive_clusters = clusters,
    oncodrive_cluster_score = clusterScores,
    oncodrive_mutations_in_clusters = muts_in_clusters,
    oncodrive_fraction_in_clusters = fract_muts_in_clusters,
    oncodrive_p = pval,
    oncodrive_fdr = fdr
  )],
  by = "Hugo_Symbol",
  all.x = TRUE,
  sort = FALSE
)
mutation_summary[, mutation_frequency_96 := mutated_patients / 96]
mutation_summary[, recurrent_mutation := mutated_patients >= 5L]
setorder(mutation_summary, -mutated_patients, -nonsynonymous_variants, Hugo_Symbol)

message("[4/5] 计算高频突变基因的患者级共现/互斥")
interaction_gene_ids <- head(mutation_summary$gene_id, 25L)
interaction_rows <- vector("list", choose(length(interaction_gene_ids), 2L))
counter <- 0L
if (length(interaction_gene_ids) >= 2L) {
  mut_mat <- as.matrix(mutation_binary[interaction_gene_ids, , drop = FALSE] > 0)
  for (i in seq_len(length(interaction_gene_ids) - 1L)) {
    for (j in seq.int(i + 1L, length(interaction_gene_ids))) {
      counter <- counter + 1L
      a <- as.logical(mut_mat[i, ])
      b <- as.logical(mut_mat[j, ])
      table_2x2 <- matrix(c(
        sum(a & b), sum(a & !b), sum(!a & b), sum(!a & !b)
      ), nrow = 2L, byrow = TRUE)
      ft <- fisher.test(table_2x2)
      interaction_rows[[counter]] <- data.table(
        gene1_id = interaction_gene_ids[[i]],
        gene1 = mutation_summary[gene_id == interaction_gene_ids[[i]], Hugo_Symbol][[1]],
        gene2_id = interaction_gene_ids[[j]],
        gene2 = mutation_summary[gene_id == interaction_gene_ids[[j]], Hugo_Symbol][[1]],
        both_mutated = table_2x2[1, 1],
        gene1_only = table_2x2[1, 2],
        gene2_only = table_2x2[2, 1],
        neither = table_2x2[2, 2],
        odds_ratio = unname(ft$estimate),
        p_value = ft$p.value
      )
    }
  }
}
interactions <- rbindlist(interaction_rows, fill = TRUE)
if (nrow(interactions)) {
  interactions[, q_value := p.adjust(p_value, method = "BH")]
  interactions[, interaction_class := fifelse(
    q_value <= 0.10 & odds_ratio > 1, "co_occurrence",
    fifelse(q_value <= 0.10 & odds_ratio < 1, "mutual_exclusivity", "descriptive_only")
  )]
  setorder(interactions, q_value, p_value)
}

message("[5/5] 构建驱动候选收敛表和柔性决策层级")
candidate <- merge(cnv_summary, dosage, by = c("gene_id", "gene_name", "gene_type"), all = TRUE)
candidate <- merge(
  candidate,
  mutation_summary[, .(
    gene_id,
    mutation_symbol = Hugo_Symbol,
    nonsynonymous_variants,
    mutated_patients,
    likely_lof_variants,
    likely_lof_patients,
    hotspot_variants,
    max_observed_aa_position,
    oncodrive_cluster_score,
    oncodrive_fraction_in_clusters,
    oncodrive_p,
    oncodrive_fdr,
    mutation_frequency_96,
    recurrent_mutation
  )],
  by = "gene_id",
  all.x = TRUE
)
candidate[is.na(mutated_patients), `:=`(
  nonsynonymous_variants = 0L,
  mutated_patients = 0L,
  likely_lof_variants = 0L,
  likely_lof_patients = 0L,
  hotspot_variants = 0L,
  max_observed_aa_position = NA_real_,
  mutation_frequency_96 = 0,
  recurrent_mutation = FALSE
)]
candidate[, gene_span_bp := end - start + 1]
candidate[, long_gene_passenger_flag :=
            gene_span_bp > 500000 & is.na(oncodrive_fdr) &
            hotspot_variants == 0 & likely_lof_patients < 5]
candidate[, large_protein_mutability_flag :=
            !is.na(max_observed_aa_position) & max_observed_aa_position > 2000 &
            (is.na(oncodrive_fdr) | oncodrive_fdr > 0.10) &
            hotspot_variants < 2 & likely_lof_patients < 5]
candidate[, oncodrive_support_only :=
            !is.na(oncodrive_fdr) & oncodrive_fdr <= 0.10]
candidate[, high_confidence_mutation_pattern :=
            (hotspot_variants >= 2 & mutated_patients >= 3) |
            (mutated_patients >= 10 & gene_span_bp <= 200000 &
               !large_protein_mutability_flag) |
            (likely_lof_patients >= 5 & mutated_patients >= 5)]
candidate[, recurrent_high_level_cnv :=
            pmax(amplification_frequency, homozygous_deletion_frequency, na.rm = TRUE) >= 0.10]
candidate[, recurrent_broad_cnv :=
            pmax(gain_frequency, loss_frequency, na.rm = TRUE) >= 0.20]
candidate[, strong_dosage :=
            !is.na(dosage_evidence) & dosage_evidence == "strong_positive_dosage"]
candidate[, conditional_dosage :=
            !is.na(dosage_evidence) & dosage_evidence == "conditional_positive_dosage"]
candidate[, precise_reverse_dosage :=
            !is.na(dosage_evidence) &
              dosage_evidence == "precise_reverse_for_dosage_axis"]
candidate[, multi_event_convergence := recurrent_mutation &
            recurrent_high_level_cnv & strong_dosage]
candidate[, multi_event_convergence_status := fifelse(
  multi_event_convergence,
  "conditional_pending_cnv_focality",
  "not_met"
)]
# CNV 复发与 CNV-RNA 剂量耦联来自同一 CNV 事件，不能当作两个独立证据单元。
candidate[, evidence_unit_count :=
            as.integer(recurrent_mutation) +
              as.integer(recurrent_high_level_cnv &
                           (strong_dosage | conditional_dosage))]
candidate[, candidate_model_eligible :=
            gene_type == "protein_coding" | recurrent_mutation]
candidate[, decision := fifelse(
  candidate_model_eligible & high_confidence_mutation_pattern,
  "strong_patient_level_candidate",
  fifelse(
    candidate_model_eligible & (recurrent_mutation |
      (recurrent_high_level_cnv & (strong_dosage | conditional_dosage))),
    "conditional_candidate",
    fifelse(
      candidate_model_eligible &
        (recurrent_high_level_cnv | strong_dosage | conditional_dosage),
      "exploratory_signal",
      "background_only"
    )
  )
)]
candidate[, primary_candidate_route := fifelse(
  high_confidence_mutation_pattern,
  "mutation_recurrence_hotspot_or_multihit_lof",
  fifelse(
    multi_event_convergence,
    "mutation_plus_cnv_dosage_conditional_pending_focality",
    fifelse(
      recurrent_mutation,
      "recurrent_mutation_conditional",
      fifelse(
        recurrent_high_level_cnv & (strong_dosage | conditional_dosage),
        "cnv_dosage_conditional_pending_focality",
        "broad_cnv_or_single_dosage_exploratory"
      )
    )
  )
)]
candidate[, decision_basis := fifelse(
  decision == "strong_patient_level_candidate",
  "高复发且非长基因频率假象，或存在多患者热点/截断模式；Oncodrive 不单独升级；剂量轴另行裁决",
  fifelse(
    decision == "conditional_candidate",
    "存在复发突变，或复发高水平 CNV 与方向稳定剂量证据；CNV 复发和剂量不重复计数，需独立校准和 focality",
    fifelse(
      decision == "exploratory_signal",
      "存在单层或宽阈值信号；保留探索，不进入主结论",
      "当前 TCGA 驱动核心未提供足够收敛证据"
    )
  )
)]
candidate[, conclusion_ceiling := fifelse(
  decision == "strong_patient_level_candidate",
  fifelse(
    high_confidence_mutation_pattern,
    "患者级高优先级突变候选；若剂量为精确反向，不得同步声称 CNV 剂量轴",
    "患者级高优先级多事件候选轴，待独立队列和同患者跨层验证"
  ),
  fifelse(
    decision == "conditional_candidate",
    "条件性候选，不得写成普适机制",
    fifelse(decision == "exploratory_signal", "单层探索线索", "背景项")
  )
)]
candidate[, required_next_validation := fifelse(
  decision == "strong_patient_level_candidate",
  "Cao 同患者多组学或独立 ESCC 队列复核；患者 bootstrap；细胞来源定位",
  fifelse(
    decision == "conditional_candidate",
    "独立患者队列、替代 CNV 尺度或情境分层至少一项",
    "仅在其他组学提供正交证据后重新评估"
  )
)]
fail_if(any(candidate$decision == "strong_patient_level_candidate" &
              !candidate$high_confidence_mutation_pattern),
        "强候选不得由 Oncodrive 或未校准 CNV 收敛单独升级")
fail_if(anyNA(candidate$evidence_unit_count),
        "证据单元计数不得含 NA")
candidate_screen <- candidate[decision != "background_only"]
decision_order <- c(
  strong_patient_level_candidate = 1L,
  conditional_candidate = 2L,
  exploratory_signal = 3L
)
candidate_screen[, decision_rank := decision_order[decision]]
setorder(
  candidate_screen,
  decision_rank,
  -high_confidence_mutation_pattern,
  oncodrive_fdr,
  -evidence_unit_count,
  -mutated_patients,
  -amplification_frequency,
  -homozygous_deletion_frequency,
  -spearman_rho_relative_cnv
)
candidate_screen[, decision_rank := NULL]

qc <- data.table(
  scope = c(
    "cohort", "matrix", "matrix", "matrix", "clinical", "method", "method",
    "method", "software", "software"
  ),
  metric = c(
    "driver_core_patients", "common_gene_universe", "dosage_test_eligible_genes",
    "maf_nonsynonymous_gene_count", "ploidy_proxy_range",
    "recurrent_mutation_definition", "recurrent_high_level_cnv_definition",
    "dosage_gate", "maftools_gene_summary_rows", "maftools_mutated_sample_rows"
  ),
  observed = c(
    length(driver_patients), length(common_gene_ids), sum(test_eligible),
    nrow(mutation_summary),
    paste0(sprintf("%.2f", min(ploidy$ploidy_proxy)), "-",
           sprintf("%.2f", max(ploidy$ploidy_proxy))),
    ">=5/96 patients", ">=10% amplification or homozygous deletion",
    "strong: rho>=0.30 and BH q<=0.10; conditional: rho>=0.30, p<=0.05 and scale-stable",
    nrow(maftools_gene_summary), nrow(maftools_sample_summary)
  ),
  expected = c(
    "95", "60623", "data-derived", "data-derived", "positive finite",
    "soft candidate threshold", "soft candidate threshold",
    "q is not the only gate", "package cross-check", "zero-event patients may be absent"
  ),
  gate_type = c("hard", "hard", rep("soft", 8)),
  status = c("pass", "pass", rep("pass", 8)),
  note = c(
    "不压缩到 76 例",
    "RNA 与 ASCAT2 gene CNV 精确交集",
    "仅 protein-coding、表达/覆盖/变异合格基因进入多重检验",
    "Ensembl 映射后非同义突变",
    "segment 长度加权中位数，仅作相对 CNV 尺度",
    "5% 左右为候选筛选，不是因果门禁",
    "阈值需在替代 CNV 尺度中复核",
    "相容性阴性与低功效反向可条件通过；精确反向不得用总分抵消",
    paste0("maftools 用于 MAF 结构与基因摘要复核；Oncodrive 聚簇基因 ",
           nrow(oncodrive_result), " 个，仅作候选证据"),
    "MAF 零事件是合法结果"
  )
)

decision_counts <- candidate_screen[, .N, by = decision]
get_decision_count <- function(label) {
  value <- decision_counts[decision == label, N]
  if (length(value)) value[[1]] else 0L
}
top_candidates <- head(candidate_screen, 20L)
top_lines <- if (nrow(top_candidates)) {
  paste0(
    seq_len(nrow(top_candidates)), ". `",
    fifelse(!is.na(top_candidates$gene_name) & top_candidates$gene_name != "",
            top_candidates$gene_name, top_candidates$gene_id),
    "`：", top_candidates$decision,
    "；突变患者 ", top_candidates$mutated_patients,
    "；扩增频率 ", sprintf("%.1f%%", 100 * top_candidates$amplification_frequency),
    "；纯合缺失频率 ", sprintf("%.1f%%", 100 * top_candidates$homozygous_deletion_frequency),
    "；CNV-RNA rho ", sprintf("%.3f", top_candidates$spearman_rho_relative_cnv)
  )
} else {
  "未形成候选"
}
summary_lines <- c(
  "# TCGA-ESCC 驱动核心首轮分析摘要",
  "",
  "## 分析边界",
  "",
  "- 主分析使用 95 例 RNA–突变–ASCAT2 CNV driver core；突变复发频率保留 96 例可用 MAF 母集。",
  "- segment 长度加权中位 copy number 仅作患者级 ploidy proxy；同时保留绝对 CNV 与相对 CNV 两种尺度。",
  "- CNV-RNA 剂量耦联只在 protein-coding、表达、覆盖和 CNV 变异合格的基因中做多重检验。",
  "- q 值不是唯一门禁：方向明确、名义效应且两种 CNV 尺度稳定者可条件保留；同语境精确反向不得由其他分数抵消。",
  "- maftools Oncodrive 只作蛋白位点聚簇线索；因本队列同义背景不足而使用预置背景，且该算法已被 OncodriveCLUSTL 取代，不能单独升级候选。",
  "- CNV 复发与 CNV-RNA 剂量耦联来自同一 CNV 事件，只计作一个 CNV-表达证据单元；在完成 focality/染色体臂背景校准前，多事件收敛最高为条件候选。",
  "",
  "## 候选层级",
  "",
  paste0(
    "- strong patient-level candidate：",
    get_decision_count("strong_patient_level_candidate"),
    " 个。"
  ),
  paste0(
    "- conditional candidate：",
    get_decision_count("conditional_candidate"),
    " 个。"
  ),
  paste0(
    "- exploratory signal：",
    get_decision_count("exploratory_signal"),
    " 个。"
  ),
  "",
  "## 排名前 20 的数据驱动候选",
  "",
  top_lines,
  "",
  "## 结论上限",
  "",
  "本表是 TCGA 患者级候选筛选，不是驱动因果证明。候选必须在 Cao 同患者多组学、独立 ESCC 队列、细胞来源或后续实验中继续升级。"
)

outputs <- list(
  tcga_escc_mutation_gene_summary.tsv = mutation_summary,
  tcga_escc_mutation_cluster_summary.tsv = oncodrive_result,
  tcga_escc_mutation_pairwise_interactions.tsv = interactions,
  tcga_escc_cnv_gene_summary.tsv = cnv_summary,
  tcga_escc_cnv_expression_dosage.tsv = dosage,
  tcga_escc_driver_candidate_screen.tsv = candidate_screen,
  tcga_escc_driver_core_qc.tsv = qc
)
for (name in names(outputs)) {
  fwrite(outputs[[name]], file.path(stage_dir, name), sep = "\t", na = "")
}
writeLines(
  summary_lines,
  file.path(stage_dir, "tcga_escc_driver_core_summary.md"),
  useBytes = TRUE
)

artifact_files <- dir_ls(stage_dir, type = "file")
artifact_manifest <- data.table(
  artifact = path_file(artifact_files),
  relative_path = file.path("results", path_file(artifact_files)),
  file_size_bytes = as.numeric(file_info(artifact_files)$size),
  sha256 = vapply(
    as.character(artifact_files), digest, character(1),
    algo = "sha256", file = TRUE
  ),
  generated_date = as.character(Sys.Date()),
  generation_script = "scripts/13_analyze_tcga_escc_driver_core.R",
  source_object = "results/tcga_escc_multiassay_dr45.rds",
  status = "verified"
)
fwrite(
  artifact_manifest,
  file.path(stage_dir, "tcga_escc_driver_core_artifact_manifest.tsv"),
  sep = "\t",
  na = ""
)

stage_files <- dir_ls(stage_dir, type = "file")
for (stage_file in stage_files) {
  destination <- file.path(results_dir, path_file(stage_file))
  file_copy(stage_file, destination, overwrite = TRUE)
  fail_if(
    digest(stage_file, algo = "sha256", file = TRUE) !=
      digest(destination, algo = "sha256", file = TRUE),
    paste0("发布后 SHA256 不一致：", path_file(stage_file))
  )
}

if (dir_exists(stage_dir)) dir_delete(stage_dir)

message("完成：TCGA-ESCC driver core 首轮结果已写入 results/")
message("候选数：", nrow(candidate_screen),
        "；剂量检验基因：", sum(test_eligible),
        "；突变基因：", nrow(mutation_summary))
