#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(fs)
  library(matrixStats)
  library(MultiAssayExperiment)
  library(SummarizedExperiment)
})

options(stringsAsFactors = FALSE)

# 将 TCGA-ESCC 的患者级强驱动候选与同一批患者的连续多组学状态相连。
#
# 这是同一 TCGA 队列内的“驱动事件 -> 状态”桥接，不是独立验证：
# - MOFA 因子来自同一 94 例患者，且包含 Mutation/CNV 等输入视图；
# - PROGENy 活性来自同一批患者的 RNA；
# - 任何显著关联都只能用于候选网络排序，不能证明因果方向。

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results")
work_intermediate_dir <- file.path(project_root, "_work", "intermediate")
required_inputs <- file.path(project_root, c(
  "PROJECT_INDEX.md",
  "results/tcga_escc_multiassay_dr45.rds",
  "results/tcga_escc_ascat2_segments_long.rds",
  "results/tcga_escc_driver_candidate_screen.tsv",
  "results/tcga_escc_mofa_factor_scores.tsv",
  "results/tcga_escc_progeny_pathway_scores.tsv",
  "results/tcga_escc_mofa_level_factor_qc.tsv",
  "results/tcga_escc_heterogeneity_feature_manifest.tsv"
))
if (any(!file_exists(required_inputs))) {
  stop(
    "缺少 TCGA 驱动—状态桥接输入：",
    paste(required_inputs[!file_exists(required_inputs)], collapse = ", ")
  )
}

dir_create(work_intermediate_dir, recurse = TRUE)
stage_dir <- tempfile(pattern = ".tcga_driver_state_", tmpdir = work_intermediate_dir)
dir_create(stage_dir)
on.exit({
  if (dir_exists(stage_dir)) dir_delete(stage_dir)
}, add = TRUE)

fail_if <- function(condition, message) {
  if (isTRUE(condition)) stop(message, call. = FALSE)
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
  mann_whitney_u <- sum(ranks[event]) - n_event * (n_event + 1) / 2
  2 * mann_whitney_u / (n_event * n_reference) - 1
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

bootstrap_effect <- function(predictor, outcome, event_type,
                             n_bootstrap = 500L, seed = 20260711L) {
  set.seed(seed)
  keep <- is.finite(predictor) & is.finite(outcome)
  predictor <- predictor[keep]
  outcome <- outcome[keep]
  estimates <- rep(NA_real_, n_bootstrap)

  if (event_type == "relative_cnv") {
    n <- length(predictor)
    for (b in seq_len(n_bootstrap)) {
      index <- sample.int(n, n, replace = TRUE)
      if (length(unique(predictor[index])) >= 3L &&
          length(unique(outcome[index])) >= 3L) {
        estimates[[b]] <- suppressWarnings(cor(
          predictor[index], outcome[index], method = "spearman"
        ))
      }
    }
  } else {
    event <- as.logical(predictor)
    event_index <- which(event)
    reference_index <- which(!event)
    for (b in seq_len(n_bootstrap)) {
      sampled_event <- sample(event_index, length(event_index), replace = TRUE)
      sampled_reference <- sample(
        reference_index, length(reference_index), replace = TRUE
      )
      sampled_outcome <- c(outcome[sampled_event], outcome[sampled_reference])
      sampled_group <- c(
        rep(TRUE, length(sampled_event)),
        rep(FALSE, length(sampled_reference))
      )
      estimates[[b]] <- binary_rank_biserial(sampled_group, sampled_outcome)
    }
  }

  estimates <- estimates[is.finite(estimates)]
  if (!length(estimates)) {
    return(list(
      valid_bootstraps = 0L,
      ci95_low = NA_real_,
      ci95_high = NA_real_,
      direction_consistency = NA_real_
    ))
  }
  list(
    valid_bootstraps = length(estimates),
    ci95_low = unname(quantile(estimates, 0.025, names = FALSE)),
    ci95_high = unname(quantile(estimates, 0.975, names = FALSE)),
    estimates = estimates
  )
}

stage_tsv <- function(object, filename) {
  fwrite(object, file.path(stage_dir, filename), sep = "\t", na = "")
}

message("[1/5] 读取 94 例 five-layer core 与 12 个强驱动候选")
mae <- readRDS(file.path(results_dir, "tcga_escc_multiassay_dr45.rds"))
validObject(mae)
segments <- as.data.table(readRDS(
  file.path(results_dir, "tcga_escc_ascat2_segments_long.rds")
))
driver_screen <- fread(
  file.path(results_dir, "tcga_escc_driver_candidate_screen.tsv"),
  na.strings = c("", "NA"),
  showProgress = FALSE
)
factor_scores <- fread(
  file.path(results_dir, "tcga_escc_mofa_factor_scores.tsv"),
  showProgress = FALSE
)
pathway_scores <- fread(
  file.path(results_dir, "tcga_escc_progeny_pathway_scores.tsv"),
  showProgress = FALSE
)
level_factor_qc <- fread(
  file.path(results_dir, "tcga_escc_mofa_level_factor_qc.tsv"),
  showProgress = FALSE
)
feature_manifest <- fread(
  file.path(results_dir, "tcga_escc_heterogeneity_feature_manifest.tsv"),
  showProgress = FALSE
)

strong <- driver_screen[decision == "strong_patient_level_candidate"]
patients <- factor_scores$patient_id
fail_if(nrow(strong) != 12L, "当前强驱动候选必须为 12 个")
fail_if(length(patients) != 94L || uniqueN(patients) != 94L,
        "MOFA 因子分数必须覆盖 94 位唯一患者")
fail_if(!identical(pathway_scores$patient_id, patients),
        "PROGENy 与 MOFA 患者顺序不一致")
fail_if(anyDuplicated(strong$gene_id) > 0L || anyDuplicated(strong$gene_name) > 0L,
        "强驱动候选基因 ID 或名称重复")

mutation_se <- experiments(mae)[["Mutation"]]
cnv_se <- experiments(mae)[["CNV_gene"]]
mutation_binary <- assay(mutation_se, "binary")
absolute_cnv <- assay(cnv_se, "copy_number")
fail_if(!all(strong$gene_id %chin% rownames(mutation_binary)),
        "强驱动候选未全部进入 Mutation 矩阵")
fail_if(!all(strong$gene_id %chin% rownames(absolute_cnv)),
        "强驱动候选未全部进入 CNV 矩阵")
fail_if(!all(patients %chin% colnames(mutation_binary)) ||
          !all(patients %chin% colnames(absolute_cnv)),
        "94 例患者未全部进入 Mutation/CNV 矩阵")

message("[2/5] 重建 ploidy-adjusted CNV 并冻结患者级事件表")
fail_if(!all(c("Chromosome", "Start", "End", "Copy_Number", "patient_id") %chin%
               names(segments)),
        "ASCAT2 segment 长表字段不完整")
segments <- segments[Chromosome %chin% paste0("chr", 1:22)]
segments[, segment_length := End - Start + 1]
ploidy <- segments[, .(
  ploidy_proxy = weighted_ploidy_proxy(Copy_Number, segment_length)
), by = patient_id]
ploidy_vector <- ploidy$ploidy_proxy[
  match(colnames(absolute_cnv), ploidy$patient_id)
]
fail_if(any(!is.finite(ploidy_vector)) || any(ploidy_vector <= 0),
        "部分 CNV 患者缺少有效 ploidy proxy")
relative_cnv <- log2(
  (absolute_cnv + 0.5) /
    matrix(
      ploidy_vector + 0.5,
      nrow = nrow(absolute_cnv),
      ncol = ncol(absolute_cnv),
      byrow = TRUE
    )
)

# 这里不沿用源表 evidence_unit_count 的窄语义：GNAS、FBXW7 的多患者热点
# 虽未达到 recurrent_mutation>=5，也应计为一个 mutation evidence unit。
strong[, mutation_evidence_unit :=
         high_confidence_mutation_pattern | recurrent_mutation]
strong[, cnv_expression_evidence_unit :=
         recurrent_high_level_cnv & (strong_dosage | conditional_dosage)]
strong[, distinct_event_unit_count_recomputed :=
         as.integer(mutation_evidence_unit) +
           as.integer(cnv_expression_evidence_unit)]
fail_if(any(!strong$mutation_evidence_unit),
        "12 个强驱动候选均应具有 mutation evidence unit")

event_rows <- vector("list", nrow(strong))
for (i in seq_len(nrow(strong))) {
  gene_id <- strong$gene_id[[i]]
  absolute <- as.numeric(absolute_cnv[gene_id, patients])
  relative <- as.numeric(relative_cnv[gene_id, patients])
  mutation <- as.integer(mutation_binary[gene_id, patients] > 0)
  event_rows[[i]] <- data.table(
    patient_id = patients,
    gene_id = gene_id,
    gene_name = strong$gene_name[[i]],
    tcga_driver_tier = strong$decision[[i]],
    primary_candidate_route = strong$primary_candidate_route[[i]],
    mutation = mutation,
    absolute_copy_number = absolute,
    relative_cnv = relative,
    amplification = as.integer(relative >= 0.8 & absolute >= 4),
    homozygous_deletion = as.integer(absolute == 0),
    mutation_evidence_unit = strong$mutation_evidence_unit[[i]],
    cnv_expression_evidence_unit = strong$cnv_expression_evidence_unit[[i]],
    distinct_event_unit_count_recomputed =
      strong$distinct_event_unit_count_recomputed[[i]],
    source_evidence_unit_count = strong$evidence_unit_count[[i]],
    independence_group = "TCGA_ESCC_DR45_five_layer_core_94",
    conclusion_ceiling = paste(
      "TCGA 内部患者级事件表；不等于外部验证或因果驱动证明"
    )
  )
}
patient_events <- rbindlist(event_rows)
fail_if(nrow(patient_events) != 94L * 12L,
        "患者级强驱动事件表应为 94×12 行")

message("[3/5] 检验驱动事件与 MOFA 因子/PROGENy 通路的患者级关联")
factor_columns <- setdiff(names(factor_scores), "patient_id")
pathway_columns <- setdiff(names(pathway_scores), "patient_id")
target_tables <- list(
  MOFA_factor = list(table = factor_scores, targets = factor_columns),
  PROGENy_pathway = list(table = pathway_scores, targets = pathway_columns)
)
event_types <- c("mutation", "relative_cnv", "amplification", "homozygous_deletion")
association_rows <- vector(
  "list",
  nrow(strong) * (length(factor_columns) + length(pathway_columns)) *
    length(event_types)
)
counter <- 0L
for (i in seq_len(nrow(strong))) {
  current_gene_id <- strong$gene_id[[i]]
  gene_name <- strong$gene_name[[i]]
  gene_events <- patient_events[gene_id == current_gene_id]
  gene_events <- gene_events[match(patients, patient_id)]
  fail_if(!identical(gene_events$patient_id, patients),
          paste("患者事件顺序错误：", gene_name))

  for (target_type in names(target_tables)) {
    target_table <- target_tables[[target_type]]$table
    for (target in target_tables[[target_type]]$targets) {
      outcome <- target_table[[target]]
      for (event_type in event_types) {
        predictor <- gene_events[[event_type]]
        test <- if (event_type == "relative_cnv") {
          continuous_test(predictor, outcome)
        } else {
          binary_test(predictor, outcome)
        }
        counter <- counter + 1L
        association_rows[[counter]] <- data.table(
          gene_id = current_gene_id,
          gene_name = gene_name,
          tcga_driver_tier = strong$decision[[i]],
          event_type = event_type,
          target_type = target_type,
          target = target,
          n_complete = test$n_complete,
          n_event = test$n_event,
          n_reference = test$n_reference,
          effect_measure = if (event_type == "relative_cnv") {
            "Spearman rho"
          } else {
            "rank-biserial correlation: event minus reference"
          },
          effect = test$effect,
          p_value = test$p_value,
          mutation_evidence_unit = strong$mutation_evidence_unit[[i]],
          cnv_expression_evidence_unit =
            strong$cnv_expression_evidence_unit[[i]],
          distinct_event_unit_count_recomputed =
            strong$distinct_event_unit_count_recomputed[[i]],
          independence_group = "TCGA_ESCC_DR45_driver_state_shared_patients",
          independent_validation = FALSE
        )
      }
    }
  }
}
associations <- rbindlist(association_rows[seq_len(counter)])
expected_associations <- nrow(strong) *
  (length(factor_columns) + length(pathway_columns)) * length(event_types)
fail_if(nrow(associations) != expected_associations,
        "驱动—状态关联行数错误")

associations[, analysis_eligible := fifelse(
  event_type == "relative_cnv",
  n_complete >= 80L,
  n_event >= 5L & n_reference >= 10L
)]
associations[, q_value := NA_real_]
associations[analysis_eligible == TRUE & is.finite(p_value),
  q_value := p.adjust(p_value, method = "BH"),
  by = .(target_type, event_type)
]

factor_level_flags <- level_factor_qc[, .(
  level_factor_soft_flag = any(level_factor_soft_flag)
), by = factor]
associations[, level_factor_soft_flag := FALSE]
associations[target_type == "MOFA_factor",
  level_factor_soft_flag := factor_level_flags$level_factor_soft_flag[
    match(target, factor_level_flags$factor)
  ]
]
fail_if(anyNA(associations$level_factor_soft_flag),
        "部分 MOFA 因子缺少 level-factor QC 状态")

mutation_features <- feature_manifest[
  view == "Mutation", unique(feature_id)
]
cnv_features <- feature_manifest[view == "CNV", unique(feature_id)]
associations[, predictor_feature_selected_in_target_model := FALSE]
associations[
  target_type == "MOFA_factor" & event_type == "mutation",
  predictor_feature_selected_in_target_model := gene_id %chin% mutation_features
]
associations[
  target_type == "MOFA_factor" &
    event_type %chin% c("relative_cnv", "amplification", "homozygous_deletion"),
  predictor_feature_selected_in_target_model := gene_id %chin% cnv_features
]
associations[, target_model_overlap_class := fifelse(
  target_type == "PROGENy_pathway",
  "same_patient_rna_derived_target",
  fifelse(
    predictor_feature_selected_in_target_model,
    "predictor_feature_in_mofa_input",
    "predictor_layer_in_mofa_input"
  )
)]

message("[4/5] 对合格关联执行 500 次患者 bootstrap 并柔性分层")
associations[, `:=`(
  bootstrap_iterations = 500L,
  valid_bootstraps = NA_integer_,
  effect_ci95_low = NA_real_,
  effect_ci95_high = NA_real_,
  bootstrap_direction_consistency = NA_real_
)]

eligible_indices <- which(associations$analysis_eligible)
for (row_index in eligible_indices) {
  row <- associations[row_index]
  gene_events <- patient_events[gene_id == row$gene_id]
  gene_events <- gene_events[match(patients, patient_id)]
  predictor <- gene_events[[row$event_type]]
  target_table <- target_tables[[row$target_type]]$table
  outcome <- target_table[[row$target]]
  boot <- bootstrap_effect(
    predictor,
    outcome,
    event_type = row$event_type,
    n_bootstrap = 500L,
    seed = 20260711L + row_index
  )
  direction_consistency <- if (
    is.finite(row$effect) && row$effect != 0 &&
      !is.null(boot$estimates) && length(boot$estimates)
  ) {
    mean(sign(boot$estimates) == sign(row$effect))
  } else {
    NA_real_
  }
  associations[row_index, `:=`(
    valid_bootstraps = boot$valid_bootstraps,
    effect_ci95_low = boot$ci95_low,
    effect_ci95_high = boot$ci95_high,
    bootstrap_direction_consistency = direction_consistency
  )]
}

associations[, association_status := fifelse(
  !analysis_eligible,
  "insufficient_event_count_descriptive",
  fifelse(
    is.finite(q_value) & q_value <= 0.10 & abs(effect) >= 0.30 &
      bootstrap_direction_consistency >= 0.80,
    "within_tcga_supported",
    fifelse(
      ((is.finite(q_value) & q_value <= 0.20) |
         (is.finite(p_value) & p_value <= 0.05)) &
        abs(effect) >= 0.20 & bootstrap_direction_consistency >= 0.70,
      "within_tcga_conditional",
      fifelse(
        is.finite(p_value) & p_value <= 0.10 & abs(effect) >= 0.20,
        "directional_exploratory",
        "no_clear_internal_association"
      )
    )
  )
)]
associations[
  level_factor_soft_flag == TRUE & association_status == "within_tcga_supported",
  association_status := "within_tcga_conditional_level_factor_qc"
]
associations[
  level_factor_soft_flag == TRUE & association_status == "within_tcga_conditional",
  association_status := "directional_exploratory_level_factor_qc"
]
associations[, discordance_code := fifelse(
  !analysis_eligible,
  "not_comparable",
  fifelse(
    association_status %chin% c(
      "within_tcga_supported",
      "within_tcga_conditional",
      "within_tcga_conditional_level_factor_qc"
    ),
    "support",
    fifelse(
      association_status %chin% c(
        "directional_exploratory",
        "directional_exploratory_level_factor_qc"
      ),
      "same_direction_nonsignificant",
      "same_direction_nonsignificant"
    )
  )
)]
associations[, countable_as_independent_validation := FALSE]
associations[, conclusion_ceiling := paste(
  "TCGA 内部患者级驱动事件—状态关联；MOFA/PROGENy 与预测变量共享患者",
  "或输入层，不能计作独立验证或因果边"
)]
associations[, required_next_validation := fifelse(
  target_type == "MOFA_factor",
  "去除预测事件后的因子敏感性、跨 seed 对齐及独立 ESCC 队列投影",
  "独立 ESCC 队列的预锁定通路活性验证及细胞组成敏感性"
)]

network_statuses <- c(
  "within_tcga_supported",
  "within_tcga_conditional",
  "within_tcga_conditional_level_factor_qc",
  "directional_exploratory",
  "directional_exploratory_level_factor_qc"
)
network_edges <- associations[
  analysis_eligible == TRUE & association_status %chin% network_statuses,
  .(
    edge_id = sprintf("TCGA_DRIVER_STATE_%04d", .I),
    source_node = paste(gene_name, event_type, sep = ":"),
    source_gene_id = gene_id,
    source_gene_name = gene_name,
    event_type,
    target_node = target,
    target_type,
    effect_measure,
    effect,
    effect_ci95_low,
    effect_ci95_high,
    p_value,
    q_value,
    bootstrap_direction_consistency,
    association_status,
    level_factor_soft_flag,
    target_model_overlap_class,
    predictor_feature_selected_in_target_model,
    independence_group,
    countable_as_independent_validation,
    conclusion_ceiling,
    required_next_validation
  )
]
setorder(network_edges, target_type, association_status, q_value, p_value)
setorder(
  associations,
  target_type,
  association_status,
  q_value,
  p_value,
  gene_name,
  event_type,
  target,
  na.last = TRUE
)

message("[5/5] 写出正式结果、摘要和 SHA256 manifest")
stage_tsv(patient_events, "tcga_escc_strong_driver_patient_events.tsv")
stage_tsv(associations, "tcga_escc_driver_state_associations.tsv")
stage_tsv(network_edges, "tcga_escc_driver_state_network_edges.tsv")

factor_primary <- network_edges[
  target_type == "MOFA_factor" &
    association_status %chin% c("within_tcga_supported", "within_tcga_conditional")
]
pathway_primary <- network_edges[
  target_type == "PROGENy_pathway" &
    association_status %chin% c("within_tcga_supported", "within_tcga_conditional")
]
top_edge_lines <- function(table, n = 10L) {
  if (!nrow(table)) return("- 无达到当前柔性内部关联门禁的边。")
  table <- copy(table)
  setorder(table, q_value, p_value, na.last = TRUE)
  table <- head(table, n)
  paste0(
    seq_len(nrow(table)), ". `", table$source_node, "` → `",
    table$target_node, "`：", table$effect_measure, "=",
    sprintf("%.3f", table$effect), "；P=", format(table$p_value, digits = 3),
    "；q=", format(table$q_value, digits = 3), "；",
    table$association_status
  )
}

summary_lines <- c(
  "# TCGA-ESCC 驱动事件—连续状态桥接摘要",
  "",
  "## 分析边界",
  "",
  "- 分析单位为 five-layer core 的 94 位唯一患者；输入为 12 个 TCGA 患者级强驱动候选。",
  "- 二分类事件要求至少 5 位事件患者和 10 位参照患者；不足者保留描述，但不进入关联升级。",
  "- 连续 CNV 使用与 driver core 相同的 segment 长度加权 ploidy proxy 和相对 log2 CNV。",
  "- 每个合格关联执行 500 次患者 bootstrap；q≤0.10 且效应/方向稳定者为来源内支持，q≤0.20 或 P≤0.05 且效应稳定者允许条件通过。",
  "- MOFA 因子来自单 seed 且包含 Mutation/CNV 输入，PROGENy 活性来自同一批 RNA；所有边均为同 TCGA 内部桥接，不能计作独立队列复现。",
  "- Factor4 是 HM450 level-factor 软标记，关联自动降级但不删除。",
  "",
  "## 事件单元语义",
  "",
  paste0(
    "- 12 个强候选均重算得到 1 个 mutation evidence unit；",
    "不再因 `recurrent_mutation>=5` 的窄定义漏掉 GNAS 或 FBXW7。"
  ),
  paste0(
    "- CNV recurrence 与 CNV–RNA dosage 始终合并为一个事件单元；",
    "只有 recurrent high-level CNV 且有 strong/conditional dosage 才新增该单元。"
  ),
  paste0(
    "- 当前 distinct event unit=2 的候选：",
    paste(strong[distinct_event_unit_count_recomputed == 2L, gene_name], collapse = "、"),
    "。"
  ),
  "",
  "## 主要因子边",
  "",
  top_edge_lines(factor_primary),
  "",
  "## 主要通路边",
  "",
  top_edge_lines(pathway_primary),
  "",
  "## 结论上限",
  "",
  "这些边用于把患者级驱动候选连接到连续状态和通路读出；它们没有时间顺序、干预证据或独立患者验证，不能写成已证实的驱动网络。预测事件直接进入 MOFA 输入的边还存在表示重叠，必须在去除该事件或独立队列投影后才能升级。"
)
writeLines(
  summary_lines,
  file.path(stage_dir, "tcga_escc_driver_state_summary.md"),
  useBytes = TRUE
)

artifact_names <- c(
  "tcga_escc_strong_driver_patient_events.tsv",
  "tcga_escc_driver_state_associations.tsv",
  "tcga_escc_driver_state_network_edges.tsv",
  "tcga_escc_driver_state_summary.md"
)
artifact_paths <- file.path(stage_dir, artifact_names)
fail_if(any(!file_exists(artifact_paths)), "驱动—状态正式产物未完整生成")
manifest <- data.table(
  artifact = artifact_names,
  relative_path = file.path("results", artifact_names),
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
  generation_script = "scripts/19_analyze_tcga_escc_driver_state_links.R",
  independence_group = "TCGA_ESCC_DR45_driver_state_shared_patients",
  status = "verified"
)
stage_tsv(manifest, "tcga_escc_driver_state_artifact_manifest.tsv")

for (artifact in c(artifact_names, "tcga_escc_driver_state_artifact_manifest.tsv")) {
  file_copy(
    file.path(stage_dir, artifact),
    file.path(results_dir, artifact),
    overwrite = TRUE
  )
}

dir_delete(stage_dir)
message(
  "完成：", nrow(patient_events), " 条患者事件、",
  nrow(associations), " 条关联、", nrow(network_edges), " 条候选网络边。"
)
