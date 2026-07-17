#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cluster)
  library(ConsensusClusterPlus)
  library(data.table)
  library(digest)
  library(fs)
  library(matrixStats)
  library(mclust)
  library(MOFA2)
  library(MultiAssayExperiment)
  library(progeny)
  library(SNFtool)
  library(survival)
  library(SummarizedExperiment)
})

options(stringsAsFactors = FALSE)

# TCGA-ESCC 五层异质性分析。
#
# 采用连续潜变量（MOFA2）与患者相似网络（SNF）两条路线；
# 不预设必须得到 4 个亚型，也不把低稳定性聚类强行命名为新亚型。

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results")
work_intermediate_dir <- file.path(project_root, "_work", "intermediate")
required_inputs <- file.path(project_root, c(
  "PROJECT_INDEX.md",
  "results/tcga_escc_multiassay_dr45.rds",
  "results/tcga_escc_ascat2_segments_long.rds",
  "results/tcga_escc_multiassay_analysis_sets.tsv",
  "results/tcga_escc_driver_candidate_screen.tsv"
))
if (any(!file_exists(required_inputs))) {
  stop("缺少五层异质性分析输入：",
       paste(required_inputs[!file_exists(required_inputs)], collapse = ", "))
}
dir_create(work_intermediate_dir, recurse = TRUE)
stage_dir <- tempfile(pattern = ".tcga_heterogeneity_", tmpdir = work_intermediate_dir)
dir_create(stage_dir)
on.exit({
  if (dir_exists(stage_dir)) dir_delete(stage_dir)
}, add = TRUE)

fail_if <- function(condition, message) {
  if (isTRUE(condition)) stop(message, call. = FALSE)
}

row_zscore <- function(x) {
  means <- rowMeans(x, na.rm = TRUE)
  sds <- rowSds(x, na.rm = TRUE, useNames = FALSE)
  keep <- is.finite(sds) & sds > 0
  fail_if(!all(keep), "进入 row_zscore 的特征必须有有限非零方差")
  z <- sweep(x, 1L, means, "-")
  sweep(z, 1L, sds, "/")
}

impute_row_median <- function(x) {
  medians <- rowMedians(x, na.rm = TRUE, useNames = FALSE)
  fail_if(any(!is.finite(medians)), "存在无法进行行中位数填补的特征")
  missing <- which(!is.finite(x), arr.ind = TRUE)
  if (nrow(missing)) x[missing] <- medians[missing[, 1L]]
  x
}

pac_score <- function(consensus_matrix, lower = 0.1, upper = 0.9) {
  values <- consensus_matrix[upper.tri(consensus_matrix)]
  mean(values > lower & values < upper, na.rm = TRUE)
}

stage_tsv <- function(object, filename) {
  fwrite(object, file.path(stage_dir, filename), sep = "\t", na = "")
}

mae <- readRDS(file.path(results_dir, "tcga_escc_multiassay_dr45.rds"))
validObject(mae)
segments <- as.data.table(readRDS(
  file.path(results_dir, "tcga_escc_ascat2_segments_long.rds")
))
analysis_sets <- fread(
  file.path(results_dir, "tcga_escc_multiassay_analysis_sets.tsv"),
  showProgress = FALSE
)
driver_screen <- fread(
  file.path(results_dir, "tcga_escc_driver_candidate_screen.tsv"),
  na.strings = c("", "NA"),
  showProgress = FALSE
)

patients <- analysis_sets[
  analysis_set == "five_layer_core" & included == TRUE,
  patient_id
]
fail_if(length(patients) != 94L || uniqueN(patients) != 94L,
        "five-layer core 必须为 94 位唯一患者")

rna_se <- experiments(mae)[["RNA"]]
mirna_se <- experiments(mae)[["miRNA"]]
meth_se <- experiments(mae)[["HM450"]]
mutation_se <- experiments(mae)[["Mutation"]]
cnv_se <- experiments(mae)[["CNV_gene"]]
patient_clinical <- as.data.table(as.data.frame(colData(mae)), keep.rownames = FALSE)
patient_clinical <- patient_clinical[match(patients, patient_id)]
fail_if(anyNA(patient_clinical$patient_id), "五层患者缺少临床 colData")

message("[1/6] 选择五个组学视图的固定特征")
rna_counts_full <- assay(rna_se, "counts")
rna_logcpm_full <- assay(rna_se, "tmm_logcpm")
rna_rd <- as.data.table(as.data.frame(rowData(rna_se)))
rna_var <- rowVars(rna_logcpm_full[, patients, drop = FALSE], useNames = FALSE)
rna_eligible <- rna_rd$gene_type == "protein_coding" &
  rowSums(rna_counts_full[, patients, drop = FALSE] >= 10) >= 10L &
  is.finite(rna_var) & rna_var > 0
rna_order <- order(rna_var, decreasing = TRUE, na.last = NA)
rna_selected <- head(rna_order[rna_eligible[rna_order]], 1500L)
rna_matrix <- rna_logcpm_full[rna_selected, patients, drop = FALSE]
rownames(rna_matrix) <- rna_rd$gene_id[rna_selected]

# PROGENy 使用唯一 gene symbol；重复 symbol 选择方差最高的 Ensembl 记录。
symbol_candidates <- data.table(
  row_index = seq_len(nrow(rna_rd)),
  gene_name = rna_rd$gene_name,
  gene_type = rna_rd$gene_type,
  variance = rna_var
)[gene_type == "protein_coding" & !is.na(gene_name) & nzchar(gene_name)]
setorder(symbol_candidates, gene_name, -variance, row_index)
symbol_candidates <- symbol_candidates[, .SD[1L], by = gene_name]
progeny_input <- rna_logcpm_full[
  symbol_candidates$row_index,
  patients,
  drop = FALSE
]
rownames(progeny_input) <- symbol_candidates$gene_name
progeny_scores <- as.data.table(
  progeny(
    progeny_input,
    scale = TRUE,
    organism = "Human",
    top = 500,
    perm = 1,
    verbose = FALSE,
    return_assay = FALSE
  ),
  keep.rownames = "patient_id"
)

mirna_counts_full <- assay(mirna_se, "counts")
mirna_logcpm_full <- assay(mirna_se, "tmm_logcpm")
mirna_rd <- as.data.table(as.data.frame(rowData(mirna_se)))
mirna_var <- rowVars(mirna_logcpm_full[, patients, drop = FALSE], useNames = FALSE)
mirna_eligible <- mirna_rd$cross_mapped_file_count <= 10L &
  rowSums(mirna_counts_full[, patients, drop = FALSE] >= 5) >= 10L &
  is.finite(mirna_var) & mirna_var > 0
mirna_order <- order(mirna_var, decreasing = TRUE, na.last = NA)
mirna_selected <- head(mirna_order[mirna_eligible[mirna_order]], 300L)
mirna_matrix <- mirna_logcpm_full[mirna_selected, patients, drop = FALSE]

meth_beta_full <- assay(meth_se, "beta")[, patients, drop = FALSE]
meth_rd <- as.data.table(as.data.frame(rowData(meth_se)))
meth_missing <- rowMeans(!is.finite(meth_beta_full))
meth_var <- rowVars(meth_beta_full, na.rm = TRUE, useNames = FALSE)
meth_eligible <- meth_rd$probe_class == "CpG" & meth_missing <= 0.20 &
  is.finite(meth_var) & meth_var > 0
meth_order <- order(meth_var, decreasing = TRUE, na.last = NA)
meth_selected <- head(meth_order[meth_eligible[meth_order]], 2000L)
meth_beta <- impute_row_median(meth_beta_full[meth_selected, , drop = FALSE])
meth_beta[meth_beta < 0.001] <- 0.001
meth_beta[meth_beta > 0.999] <- 0.999
meth_matrix <- log2(meth_beta / (1 - meth_beta))
rownames(meth_matrix) <- meth_rd$probe_id[meth_selected]

mutation_full <- assay(mutation_se, "binary")[, patients, drop = FALSE]
mutation_prevalence <- Matrix::rowMeans(mutation_full)
mutation_candidates <- which(mutation_prevalence >= 3 / 94)
mutation_selected <- head(
  mutation_candidates[order(mutation_prevalence[mutation_candidates], decreasing = TRUE)],
  50L
)
mutation_matrix <- as.matrix(mutation_full[mutation_selected, , drop = FALSE])
mutation_rd <- as.data.table(as.data.frame(rowData(mutation_se)))

cnv_total_full <- assay(cnv_se, "copy_number")
cnv_rd <- as.data.table(as.data.frame(rowData(cnv_se)))
segments_autosome <- segments[Chromosome %chin% paste0("chr", 1:22)]
ploidy <- segments_autosome[, .(
  ploidy_proxy = weightedMedian(Copy_Number, w = segment_length, na.rm = TRUE)
), by = patient_id]
ploidy_vector <- ploidy$ploidy_proxy[match(patients, ploidy$patient_id)]
fail_if(anyNA(ploidy_vector), "五层患者缺少 ploidy proxy")
cnv_absolute <- cnv_total_full[, patients, drop = FALSE]
cnv_relative <- log2(
  (cnv_absolute + 0.5) /
    matrix(ploidy_vector + 0.5, nrow = nrow(cnv_absolute),
           ncol = ncol(cnv_absolute), byrow = TRUE)
)
cnv_gene_type <- rna_rd$gene_type[match(rownames(cnv_relative), rna_rd$gene_id)]
cnv_var <- rowVars(cnv_relative, na.rm = TRUE, useNames = FALSE)
cnv_feature_table <- data.table(
  row_index = seq_len(nrow(cnv_rd)),
  gene_id = cnv_rd$gene_id,
  gene_name = cnv_rd$gene_name,
  chromosome = cnv_rd$chromosome,
  gene_type = cnv_gene_type,
  variance = cnv_var,
  coverage = rowMeans(is.finite(cnv_relative))
)[chromosome %chin% paste0("chr", 1:22) & gene_type == "protein_coding" &
    coverage >= 0.90 & is.finite(variance) & variance > 0]
setorder(cnv_feature_table, chromosome, -variance, row_index)
cnv_selected_table <- cnv_feature_table[, head(.SD, 25L), by = chromosome]
setorder(cnv_selected_table, chromosome, -variance)
cnv_selected <- cnv_selected_table$row_index
cnv_matrix <- impute_row_median(cnv_relative[cnv_selected, , drop = FALSE])

fail_if(length(rna_selected) != 1500L || length(mirna_selected) != 300L ||
          length(meth_selected) != 2000L || length(mutation_selected) < 10L ||
          length(cnv_selected) < 400L,
        "固定特征选择未达到预期维度")

mofa_views <- list(
  RNA = row_zscore(rna_matrix),
  miRNA = row_zscore(mirna_matrix),
  HM450 = row_zscore(meth_matrix),
  Mutation = mutation_matrix,
  CNV = row_zscore(cnv_matrix)
)
fail_if(any(vapply(mofa_views, function(x) !identical(colnames(x), patients), logical(1))),
        "五个视图的患者顺序不一致")

feature_manifest <- rbindlist(list(
  data.table(
    view = "RNA",
    feature_id = rna_rd$gene_id[rna_selected],
    feature_name = rna_rd$gene_name[rna_selected],
    selection_metric = "variance_tmm_logcpm",
    selection_value = rna_var[rna_selected],
    rank = seq_along(rna_selected),
    preprocessing = "protein_coding;count>=10_in_10;row_zscore"
  ),
  data.table(
    view = "miRNA",
    feature_id = mirna_rd$miRNA_ID[mirna_selected],
    feature_name = mirna_rd$miRNA_ID[mirna_selected],
    selection_metric = "variance_tmm_logcpm",
    selection_value = mirna_var[mirna_selected],
    rank = seq_along(mirna_selected),
    preprocessing = "cross_mapped_files<=10;count>=5_in_10;row_zscore"
  ),
  data.table(
    view = "HM450",
    feature_id = meth_rd$probe_id[meth_selected],
    feature_name = meth_rd$probe_id[meth_selected],
    selection_metric = "variance_beta",
    selection_value = meth_var[meth_selected],
    rank = seq_along(meth_selected),
    preprocessing = "CpG;missing<=20%;row_median_impute;M_value;row_zscore"
  ),
  data.table(
    view = "Mutation",
    feature_id = mutation_rd$gene_id[mutation_selected],
    feature_name = mutation_rd$gene_name[mutation_selected],
    selection_metric = "patient_prevalence",
    selection_value = mutation_prevalence[mutation_selected],
    rank = seq_along(mutation_selected),
    preprocessing = "binary;mutated_patients>=3"
  ),
  data.table(
    view = "CNV",
    feature_id = cnv_rd$gene_id[cnv_selected],
    feature_name = cnv_rd$gene_name[cnv_selected],
    selection_metric = "variance_ploidy_adjusted_log2",
    selection_value = cnv_var[cnv_selected],
    rank = seq_along(cnv_selected),
    preprocessing = "autosome;protein_coding;coverage>=90%;top25_per_chromosome;row_zscore"
  )
), use.names = TRUE)

rm(
  mae, rna_counts_full, rna_logcpm_full, mirna_counts_full, mirna_logcpm_full,
  meth_beta_full, meth_beta, mutation_full, cnv_total_full, cnv_absolute,
  cnv_relative, rna_se, mirna_se, meth_se, mutation_se, cnv_se
)
gc(verbose = FALSE)

message("[2/6] 训练 MOFA2 连续潜变量模型")
mofa_object <- create_mofa(mofa_views)
data_options <- get_default_data_options(mofa_object)
data_options$scale_views <- FALSE
data_options$center_groups <- TRUE
data_options$use_float32 <- TRUE
model_options <- get_default_model_options(mofa_object)
model_options$num_factors <- 10L
model_options$likelihoods <- c(
  RNA = "gaussian",
  miRNA = "gaussian",
  HM450 = "gaussian",
  Mutation = "bernoulli",
  CNV = "gaussian"
)
model_options$spikeslab_weights <- TRUE
model_options$ard_weights <- TRUE
training_options <- get_default_training_options(mofa_object)
training_options$maxiter <- 1500L
training_options$convergence_mode <- "medium"
training_options$drop_factor_threshold <- 0.01
training_options$seed <- 20260711L
training_options$verbose <- FALSE
mofa_object <- prepare_mofa(
  mofa_object,
  data_options = data_options,
  model_options = model_options,
  training_options = training_options
)
temporary_hdf5 <- file.path(stage_dir, "temporary_mofa_model.hdf5")
mofa_fit <- run_mofa(
  mofa_object,
  outfile = temporary_hdf5,
  save_data = TRUE,
  use_basilisk = TRUE
)
factor_matrix <- get_factors(mofa_fit, factors = "all")[[1]]
fail_if(nrow(factor_matrix) != 94L || !all(patients %chin% rownames(factor_matrix)),
        "MOFA2 factor scores 患者维度异常")
factor_matrix <- factor_matrix[patients, , drop = FALSE]
factor_scores <- as.data.table(factor_matrix, keep.rownames = "patient_id")
view_mean_levels <- vapply(
  mofa_views,
  function(x) colMeans(x, na.rm = TRUE),
  numeric(length(patients))
)
rownames(view_mean_levels) <- patients
level_factor_qc <- rbindlist(lapply(seq_len(ncol(factor_matrix)), function(i) {
  rbindlist(lapply(seq_len(ncol(view_mean_levels)), function(j) {
    test <- suppressWarnings(cor.test(
      factor_matrix[, i],
      view_mean_levels[, j],
      method = "spearman",
      exact = FALSE
    ))
    data.table(
      factor = colnames(factor_matrix)[[i]],
      view = colnames(view_mean_levels)[[j]],
      spearman_rho_with_view_mean = unname(test$estimate),
      p_value = test$p.value,
      level_factor_soft_flag = abs(unname(test$estimate)) >= 0.80
    )
  }))
}))
level_factor_qc[, q_value := p.adjust(p_value, method = "BH")]
level_flagged_factors <- unique(level_factor_qc[
  level_factor_soft_flag == TRUE,
  factor
])
factor_matrix_for_clustering <- factor_matrix[
  , !colnames(factor_matrix) %chin% level_flagged_factors,
  drop = FALSE
]
if (ncol(factor_matrix_for_clustering) < 2L) {
  factor_matrix_for_clustering <- factor_matrix
  level_flagged_factors <- character()
}
variance_explained <- as.data.table(
  get_variance_explained(mofa_fit, as.data.frame = TRUE)
)
weights_long <- as.data.table(get_weights(mofa_fit, as.data.frame = TRUE))
weight_value_col <- intersect(c("value", "weight"), names(weights_long))[[1]]
weights_long[, abs_weight := abs(get(weight_value_col))]
setorder(weights_long, view, factor, -abs_weight)
top_weights <- weights_long[, head(.SD, 20L), by = .(view, factor)]

message("[3/6] 计算 MOFA 因子与 PROGENy/临床的关联")
progeny_matrix <- as.matrix(progeny_scores[, -"patient_id"])
rownames(progeny_matrix) <- progeny_scores$patient_id
progeny_matrix <- progeny_matrix[patients, , drop = FALSE]
factor_pathway <- rbindlist(lapply(seq_len(ncol(factor_matrix)), function(i) {
  rbindlist(lapply(seq_len(ncol(progeny_matrix)), function(j) {
    test <- suppressWarnings(cor.test(
      factor_matrix[, i], progeny_matrix[, j],
      method = "spearman", exact = FALSE
    ))
    data.table(
      factor = colnames(factor_matrix)[[i]],
      pathway = colnames(progeny_matrix)[[j]],
      n = sum(complete.cases(factor_matrix[, i], progeny_matrix[, j])),
      spearman_rho = unname(test$estimate),
      p_value = test$p.value
    )
  }))
}))
factor_pathway[, q_value := p.adjust(p_value, method = "BH")]
factor_pathway[, abs_spearman_rho := abs(spearman_rho)]
setorder(factor_pathway, q_value, -abs_spearman_rho)

factor_clinical <- list()
fc_counter <- 0L
stage_text <- as.character(patient_clinical$ajcc_pathologic_stage)
stage_numeric <- fifelse(grepl("Stage IV", stage_text), 4,
                         fifelse(grepl("Stage III", stage_text), 3,
                                 fifelse(grepl("Stage II", stage_text), 2,
                                         fifelse(grepl("Stage I", stage_text), 1, NA_real_))))
for (i in seq_len(ncol(factor_matrix))) {
  factor_name <- colnames(factor_matrix)[[i]]
  score <- factor_matrix[, i]
  for (endpoint in c("age_at_index", "pathologic_stage")) {
    endpoint_value <- if (endpoint == "age_at_index") {
      patient_clinical$age_at_index
    } else {
      stage_numeric
    }
    keep <- complete.cases(score, endpoint_value)
    test <- suppressWarnings(cor.test(
      score[keep], endpoint_value[keep], method = "spearman", exact = FALSE
    ))
    fc_counter <- fc_counter + 1L
    factor_clinical[[fc_counter]] <- data.table(
      factor = factor_name,
      endpoint = endpoint,
      method = "spearman",
      n = sum(keep),
      effect = unname(test$estimate),
      p_value = test$p.value,
      effect_label = "rho"
    )
  }
  survival_keep <- complete.cases(
    score,
    patient_clinical$xena_os_time_days,
    patient_clinical$xena_os_event
  ) & patient_clinical$xena_os_time_days > 0
  if (sum(survival_keep) >= 30L && sum(patient_clinical$xena_os_event[survival_keep]) >= 10L) {
    fit <- coxph(
      Surv(patient_clinical$xena_os_time_days[survival_keep],
           patient_clinical$xena_os_event[survival_keep]) ~ scale(score[survival_keep])
    )
    coefficient <- summary(fit)$coefficients[1, ]
    fc_counter <- fc_counter + 1L
    factor_clinical[[fc_counter]] <- data.table(
      factor = factor_name,
      endpoint = "overall_survival",
      method = "univariable_cox",
      n = sum(survival_keep),
      effect = exp(coefficient[["coef"]]),
      p_value = coefficient[["Pr(>|z|)"]],
      effect_label = "hazard_ratio_per_1SD"
    )
  }
}
factor_clinical <- rbindlist(factor_clinical, fill = TRUE)
factor_clinical[, q_value := p.adjust(p_value, method = "BH")]
factor_clinical[, boundary := fifelse(
  endpoint == "overall_survival",
  "仅 31 个 OS 事件的单变量探索，不作预后定论",
  "患者级探索关联"
)]
setorder(factor_clinical, q_value, p_value)

message("[4/6] 构建五视图 SNF 患者相似网络")
snf_inputs <- lapply(mofa_views, function(x) {
  sample_feature <- t(x)
  sample_feature[!is.finite(sample_feature)] <- 0
  standardNormalization(sample_feature)
})
snf_affinities <- lapply(snf_inputs, function(x) {
  affinityMatrix(dist2(x, x), K = 20L, sigma = 0.5)
})
snf_network <- SNF(snf_affinities, K = 20L, t = 20L)
rownames(snf_network) <- colnames(snf_network) <- patients

# silhouette 必须在与 type=3 谱聚类一致的归一化 Laplacian 特征空间计算。
# 直接使用 1-affinity 会把 SNF 的稀疏小权重压缩到几乎相同的距离，
# 从而系统性低估 silhouette。
spectral_embedding_distance <- function(
    affinity,
    dimensions,
    exclude_trivial = FALSE,
    row_normalize = FALSE,
    column_scale = FALSE) {
  degree <- rowSums(affinity)
  degree[degree == 0] <- .Machine$double.eps
  laplacian <- diag(degree) - affinity
  degree_inv_sqrt <- diag(1 / sqrt(degree))
  normalized_laplacian <- degree_inv_sqrt %*% laplacian %*% degree_inv_sqrt
  eig <- eigen(normalized_laplacian, symmetric = TRUE)
  ordered <- order(abs(eig$values))
  offset <- if (exclude_trivial) 1L else 0L
  selected <- ordered[seq_len(dimensions) + offset]
  embedding <- eig$vectors[, selected, drop = FALSE]
  if (row_normalize) {
    norms <- sqrt(rowSums(embedding^2))
    norms[norms == 0] <- 1
    embedding <- embedding / norms
  }
  if (column_scale) embedding <- scale(embedding)
  rownames(embedding) <- patients
  dist(embedding)
}
fixed_snf_distance <- spectral_embedding_distance(
  snf_network,
  dimensions = 6L,
  exclude_trivial = TRUE,
  row_normalize = FALSE,
  column_scale = TRUE
)
snf_clusters <- lapply(2:6, function(k) {
  labels <- spectralClustering(snf_network, K = k, type = 3)
  names(labels) <- patients
  labels
})
names(snf_clusters) <- as.character(2:6)

message("[5/6] 进行 500 次 MOFA 因子共识聚类并与 SNF 三角验证")
temporary_consensus_pdf <- file.path(stage_dir, "temporary_consensus_diagnostics.pdf")
grDevices::pdf(temporary_consensus_pdf)
consensus <- ConsensusClusterPlus(
  d = t(factor_matrix_for_clustering),
  maxK = 6L,
  reps = 500L,
  pItem = 0.80,
  pFeature = 0.80,
  clusterAlg = "hc",
  innerLinkage = "average",
  finalLinkage = "average",
  distance = "pearson",
  seed = 20260711L,
  plot = NULL,
  writeTable = FALSE,
  verbose = FALSE
)
grDevices::dev.off()
if (file_exists(temporary_consensus_pdf)) file_delete(temporary_consensus_pdf)

cluster_eval <- rbindlist(lapply(2:6, function(k) {
  snf_labels <- snf_clusters[[as.character(k)]][patients]
  mofa_labels <- consensus[[k]]$consensusClass[patients]
  snf_silhouette <- mean(
    silhouette(snf_labels, fixed_snf_distance)[, "sil_width"]
  )
  snf_same_space_distance <- spectral_embedding_distance(
    snf_network,
    dimensions = k,
    exclude_trivial = FALSE,
    row_normalize = TRUE,
    column_scale = FALSE
  )
  snf_same_space_silhouette <- mean(
    silhouette(snf_labels, snf_same_space_distance)[, "sil_width"]
  )
  factor_distance <- dist(scale(factor_matrix_for_clustering))
  mofa_silhouette <- mean(silhouette(mofa_labels, factor_distance)[, "sil_width"])
  data.table(
    k = k,
    snf_mean_silhouette = snf_silhouette,
    snf_silhouette_space = "fixed_6d_normalized_laplacian_nontrivial",
    snf_same_space_silhouette_sensitivity = snf_same_space_silhouette,
    snf_min_cluster_size = min(table(snf_labels)),
    mofa_consensus_pac = pac_score(consensus[[k]]$consensusMatrix),
    mofa_mean_silhouette = mofa_silhouette,
    mofa_min_cluster_size = min(table(mofa_labels)),
    snf_mofa_adjusted_rand = adjustedRandIndex(snf_labels, mofa_labels)
  )
}))
cluster_eval[, strong_stability :=
               snf_mean_silhouette >= 0.15 & mofa_consensus_pac <= 0.20 &
               snf_mofa_adjusted_rand >= 0.40 & snf_min_cluster_size >= 8L &
               mofa_min_cluster_size >= 8L]
cluster_eval[, conditional_stability :=
               snf_mean_silhouette >= 0.05 & mofa_consensus_pac <= 0.35 &
               snf_mofa_adjusted_rand >= 0.20 & snf_min_cluster_size >= 6L &
               mofa_min_cluster_size >= 6L]
cluster_eval[, composite_rank :=
               frank(-snf_mean_silhouette, ties.method = "average") +
               frank(mofa_consensus_pac, ties.method = "average") +
               frank(-snf_mofa_adjusted_rand, ties.method = "average")]

if (cluster_eval[strong_stability == TRUE, .N]) {
  selected <- cluster_eval[strong_stability == TRUE][which.min(composite_rank)]
  heterogeneity_status <- "stable_supported_subtype_hypothesis"
} else if (cluster_eval[conditional_stability == TRUE, .N]) {
  selected <- cluster_eval[conditional_stability == TRUE][which.min(composite_rank)]
  heterogeneity_status <- "conditional_subtype_hypothesis"
} else {
  selected <- cluster_eval[which.min(composite_rank)]
  heterogeneity_status <- "no_stable_discrete_clusters_continuous_factors_retained"
}
selected_k <- selected$k[[1]]
cluster_eval[, selected_for_description := k == selected_k]
cluster_eval[, decision := fifelse(
  strong_stability, "stable",
  fifelse(conditional_stability, "conditional", "descriptive_only")
)]

snf_selected <- snf_clusters[[as.character(selected_k)]][patients]
mofa_selected <- consensus[[selected_k]]$consensusClass[patients]
patient_assignments <- merge(
  data.table(
    patient_id = patients,
    selected_k = selected_k,
    heterogeneity_status = heterogeneity_status,
    snf_cluster = as.integer(snf_selected),
    mofa_consensus_cluster = as.integer(mofa_selected)
  ),
  factor_scores,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)
patient_assignments <- merge(
  patient_assignments,
  patient_clinical[, .(
    patient_id,
    age_at_index,
    sex_at_birth,
    ajcc_pathologic_stage,
    xena_os_time_days,
    xena_os_event
  )],
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)

cluster_pathway <- rbindlist(lapply(colnames(progeny_matrix), function(pathway) {
  values <- progeny_matrix[, pathway]
  test <- kruskal.test(values ~ factor(snf_selected))
  k <- length(unique(snf_selected))
  epsilon_squared <- max(0, (unname(test$statistic) - k + 1) / (length(values) - k))
  data.table(
    pathway = pathway,
    selected_k = selected_k,
    heterogeneity_status = heterogeneity_status,
    kruskal_h = unname(test$statistic),
    epsilon_squared = epsilon_squared,
    p_value = test$p.value
  )
}))
cluster_pathway[, q_value := p.adjust(p_value, method = "BH")]
setorder(cluster_pathway, q_value, -epsilon_squared)

message("[6/6] 校验并发布异质性结果")
model_rds <- file.path(stage_dir, "tcga_escc_mofa_model.rds")
saveRDS(mofa_fit, model_rds, compress = "gzip")
if (file_exists(temporary_hdf5)) file_delete(temporary_hdf5)
mofa_check <- readRDS(model_rds)
fail_if(nrow(get_factors(mofa_check)[[1]]) != 94L, "MOFA RDS 重新读取失败")
rm(mofa_check)

top_pathway_lines <- head(factor_pathway, 10L)
top_pathway_text <- if (nrow(top_pathway_lines)) {
  paste0(
    seq_len(nrow(top_pathway_lines)), ". `",
    top_pathway_lines$factor, "`–`", top_pathway_lines$pathway,
    "`：rho=", sprintf("%.3f", top_pathway_lines$spearman_rho),
    "，q=", format(top_pathway_lines$q_value, digits = 3, scientific = TRUE)
  )
} else {
  "无可报告关联"
}
summary_lines <- c(
  "# TCGA-ESCC 五层异质性首轮分析摘要",
  "",
  "## 输入与方法",
  "",
  "- 分析母集为 94 例 five-layer core，不因缺 RPPA 压缩到 76 例。",
  paste0("- MOFA2 输入：RNA ", nrow(mofa_views$RNA), "、miRNA ",
         nrow(mofa_views$miRNA), "、HM450 ", nrow(mofa_views$HM450),
         "、Mutation ", nrow(mofa_views$Mutation), "、CNV ",
         nrow(mofa_views$CNV), " 个特征。"),
  paste0("- 当前单一 seed 的 MOFA2 模型保留 ", ncol(factor_matrix),
         " 个因子；这表示本次模型收敛，不等同于因子数已跨 seed 稳定；SNF 使用 K=20、迭代 20 次。"),
  paste0("- 与任一组学视图样本均值 |rho|≥0.80 的 level-factor 软标记为 ",
         length(level_flagged_factors), " 个；保留用于解释，但从聚类因子空间中排除。"),
  "- MOFA 因子空间进行 500 次患者/特征重采样共识聚类，并与 SNF 谱聚类比较。",
  "",
  "## 异质性门禁结果",
  "",
  paste0("- 描述性最优 k=", selected_k, "；状态：`", heterogeneity_status, "`。"),
  paste0("- SNF silhouette=", sprintf("%.3f", selected$snf_mean_silhouette),
         "，MOFA PAC=", sprintf("%.3f", selected$mofa_consensus_pac),
         "，两路线 ARI=", sprintf("%.3f", selected$snf_mofa_adjusted_rand), "。"),
  paste0("- SNF silhouette 使用所有 k 共用的固定 6 维非平凡归一化 Laplacian 距离；",
         "同 k 聚类空间的内部敏感性值为 ",
         sprintf("%.3f", selected$snf_same_space_silhouette_sensitivity),
         "，只作内部拟合参考，不参与跨 k 排名。"),
  "- 当前参数下未检出稳定离散聚类；保留连续因子作为描述框架，但这不是‘疾病必然连续’的证明，不命名新亚型。",
  "- ECMS1–4 仅作为后续外部基准，当前不会把聚类数相同当作创新或复现证明。",
  "",
  "## 最强的因子–通路关联",
  "",
  top_pathway_text,
  "",
  "## 证据上限",
  "",
  "MOFA 因子、SNF 聚类和 PROGENy 活性均为公共数据中的患者异质性假设；未经过独立队列和同患者多组学验证前，不写成固定分子亚型或因果网络。"
)

stage_tsv(feature_manifest, "tcga_escc_heterogeneity_feature_manifest.tsv")
stage_tsv(progeny_scores, "tcga_escc_progeny_pathway_scores.tsv")
stage_tsv(factor_scores, "tcga_escc_mofa_factor_scores.tsv")
stage_tsv(variance_explained, "tcga_escc_mofa_variance_explained.tsv")
stage_tsv(top_weights, "tcga_escc_mofa_top_weights.tsv")
stage_tsv(level_factor_qc, "tcga_escc_mofa_level_factor_qc.tsv")
stage_tsv(factor_pathway, "tcga_escc_mofa_pathway_associations.tsv")
stage_tsv(factor_clinical, "tcga_escc_mofa_clinical_associations.tsv")
stage_tsv(cluster_eval, "tcga_escc_heterogeneity_cluster_evaluation.tsv")
stage_tsv(patient_assignments, "tcga_escc_heterogeneity_patient_assignments.tsv")
stage_tsv(cluster_pathway, "tcga_escc_heterogeneity_cluster_pathways.tsv")
writeLines(
  summary_lines,
  file.path(stage_dir, "tcga_escc_heterogeneity_summary.md"),
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
  generation_script = "scripts/14_analyze_tcga_escc_heterogeneity.R",
  source_object = "results/tcga_escc_multiassay_dr45.rds",
  status = "verified"
)
stage_tsv(artifact_manifest, "tcga_escc_heterogeneity_artifact_manifest.tsv")

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

message("完成：TCGA-ESCC 五层异质性结果已写入 results/")
message("MOFA factors：", ncol(factor_matrix),
        "；selected k：", selected_k,
        "；status：", heterogeneity_status)
