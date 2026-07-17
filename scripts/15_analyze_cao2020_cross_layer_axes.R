#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(fs)
  library(limma)
})

options(stringsAsFactors = FALSE)

# Cao 2020 同患者多组学校准：
# 1) 9 对 WGBS 共同 CpG 的 promoter / gene-body 配对效应；
# 2) 10 对 RNA TPM 的患者固定效应 limma 与伪计数敏感性；
# 3) Supplementary Data 6 蛋白组比值的 protein-group 身份传播与 P7 批次折叠；
# 4) 结果仅用于公共数据假设生成和跨层校准，不证明因果驱动或治疗靶点。

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results")
work_intermediate_dir <- file.path(project_root, "_work", "intermediate")
work_checks_dir <- file.path(project_root, "_work", "checks")
dir_create(c(results_dir, work_intermediate_dir, work_checks_dir), recurse = TRUE)

stage_dir <- tempfile(pattern = ".cao2020_cross_layer_", tmpdir = work_intermediate_dir)
dir_create(stage_dir)
# 失败时保留 stage_dir 供审计；只在全部发布成功后显式删除。

fail_if <- function(condition, message) {
  if (isTRUE(condition)) stop(message, call. = FALSE)
}

safe_wilcox_p <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3L || all(x == 0)) return(NA_real_)
  tryCatch(
    wilcox.test(x, mu = 0, exact = FALSE, conf.int = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  n <- sum(keep)
  if (n < 4L || uniqueN(x[keep]) < 2L || uniqueN(y[keep]) < 2L) {
    return(data.table(n = n, rho = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  data.table(n = n, rho = unname(test$estimate), p_value = test$p.value)
}

strip_ensembl_version <- function(x) sub("[.][0-9]+$", "", x)

extract_gtf_attribute <- function(x, key) {
  pattern <- paste0(".*", key, " \\\"([^\\\"]+)\\\".*")
  hit <- grepl(pattern, x)
  out <- rep(NA_character_, length(x))
  out[hit] <- sub(pattern, "\\1", x[hit])
  out
}

extract_gn <- function(x) {
  matched <- regexec("(?:^| )GN=([^ ]+)", as.character(x), perl = TRUE)
  values <- regmatches(as.character(x), matched)
  vapply(values, function(z) if (length(z) >= 2L) z[[2L]] else NA_character_, character(1))
}

coverage_gate <- function(region_type, n_common_cpg) {
  if (region_type == "promoter") {
    if (n_common_cpg >= 5L) return("primary_coverage")
    if (n_common_cpg >= 3L) return("conditional_coverage")
    return("technically_unavailable")
  }
  if (n_common_cpg >= 20L) return("primary_coverage")
  if (n_common_cpg >= 10L) return("conditional_coverage")
  "technically_unavailable"
}

region_effect_direction <- function(x, threshold = 0.05) {
  fifelse(
    !is.finite(x), "not_evaluable",
    fifelse(x >= threshold, "hypermethylated_in_tumor",
            fifelse(x <= -threshold, "hypomethylated_in_tumor", "small_or_mixed_shift"))
  )
}

rna_effect_direction <- function(log_fc_pc1, direction_stable, threshold = 0.30) {
  fifelse(
    !is.finite(log_fc_pc1), "not_evaluable",
    fifelse(!direction_stable, "pseudocount_sensitive",
            fifelse(log_fc_pc1 >= threshold, "higher_in_tumor",
                    fifelse(log_fc_pc1 <= -threshold, "lower_in_tumor", "small_or_mixed_shift")))
  )
}

protein_effect_direction <- function(status, n_patients, median_effect, threshold = 0.20) {
  if (status == "not_identified") return("layer_missing_not_negative")
  if (status == "identified_not_quantified") return("identified_not_quantified_not_negative")
  if (!is.finite(n_patients) || n_patients < 3L) return("quantified_but_sparse")
  if (!is.finite(median_effect)) return("quantified_but_not_evaluable")
  if (median_effect >= threshold) return("higher_in_tumor")
  if (median_effect <= -threshold) return("lower_in_tumor")
  "small_or_mixed_shift"
}

input_paths <- c(
  driver_screen = file.path(results_dir, "tcga_escc_driver_candidate_screen.tsv")
)
fail_if(any(!file_exists(input_paths)), paste("TCGA 候选输入缺失：", paste(input_paths[!file_exists(input_paths)], collapse = ", ")))

data_root <- Sys.getenv("RESEARCH_DATA_ROOT", unset = path.expand("~/ResearchDataHub"))
fail_if(!dir_exists(data_root), paste("ResearchDataHub 不可访问：", data_root))

wgbs_root <- file.path(data_root, "datasets", "public", "GEO", "GSE149608", "retrieved_20260711_author_processed")
rna_root <- file.path(data_root, "datasets", "public", "GEO", "GSE149609", "retrieved_20260711_author_processed")
protein_root <- file.path(
  data_root, "datasets", "public", "Nature_Communications",
  "doi_10.1038_s41467-020-17227-z", "supplement_data6_snapshot_20231204"
)
gtf_path <- file.path(
  data_root, "references", "Homo_sapiens", "GRCh37", "GENCODE", "release_19",
  "00_source", "gencode.v19.annotation.gtf.gz"
)
bedtools <- Sys.getenv(
  "BEDTOOLS",
  path.expand("~/.local/share/mamba/envs/escc-amplicon/bin/bedtools")
)

required_canonical <- c(
  file.path(wgbs_root, "DATASET.md"),
  file.path(wgbs_root, "10_metadata", "sample_file_map.tsv"),
  file.path(rna_root, "DATASET.md"),
  file.path(rna_root, "10_metadata", "sample_file_map.tsv"),
  file.path(protein_root, "DATASET.md"),
  file.path(protein_root, "20_reusable", "proteomics_protein_identification.tsv"),
  file.path(protein_root, "20_reusable", "proteomics_ratio_long.tsv"),
  gtf_path,
  bedtools
)
fail_if(any(!file_exists(required_canonical)), paste(
  "Cao 2020 或 GENCODE v19 输入缺失：",
  paste(required_canonical[!file_exists(required_canonical)], collapse = ", ")
))
fail_if(file.access(bedtools, mode = 1L) != 0L, paste("bedtools 不可执行：", bedtools))

driver_screen <- fread(input_paths[["driver_screen"]], na.strings = c("", "NA"), showProgress = FALSE)
candidates <- driver_screen[
  decision == "strong_patient_level_candidate",
  .(
    tcga_gene_id = gene_id,
    tcga_gene_name = gene_name,
    tcga_primary_candidate_route = primary_candidate_route,
    tcga_decision_basis = decision_basis,
    tcga_conclusion_ceiling = conclusion_ceiling,
    tcga_required_next_validation = required_next_validation
  )
]
fail_if(!nrow(candidates), "当前 TCGA 结果中没有 strong_patient_level_candidate")
fail_if(anyDuplicated(candidates$tcga_gene_id) > 0L || anyDuplicated(candidates$tcga_gene_name) > 0L,
        "TCGA 强候选存在重复 gene_id 或 gene_name")
candidates[, candidate_rank := seq_len(.N)]
candidates[, stable_gene_id := strip_ensembl_version(tcga_gene_id)]

message("[1/7] 使用 GENCODE v19/GRCh37 冻结候选区域")
gtf_cmd <- sprintf(
  "gzip -cd %s | awk -F '\\t' '$3 == \"gene\"'",
  shQuote(gtf_path)
)
gtf_genes <- fread(
  cmd = gtf_cmd,
  sep = "\t",
  header = FALSE,
  quote = "",
  fill = TRUE,
  showProgress = FALSE
)
fail_if(ncol(gtf_genes) < 9L, "GENCODE v19 GTF gene 记录解析失败")
setnames(gtf_genes, paste0("V", seq_len(ncol(gtf_genes))))
gtf_genes <- gtf_genes[, .(
  chromosome = V1,
  feature = V3,
  gtf_start = as.integer(V4),
  gtf_end = as.integer(V5),
  strand = V7,
  attributes = V9
)]
gtf_genes[, gencode19_gene_id := extract_gtf_attribute(attributes, "gene_id")]
gtf_genes[, stable_gene_id := strip_ensembl_version(gencode19_gene_id)]
gtf_genes[, gencode19_gene_name := extract_gtf_attribute(attributes, "gene_name")]
gtf_candidates <- merge(
  candidates,
  gtf_genes[, .(
    stable_gene_id, gencode19_gene_id, gencode19_gene_name,
    chromosome, gtf_start, gtf_end, strand
  )],
  by = "stable_gene_id",
  all.x = TRUE,
  sort = FALSE
)
setorder(gtf_candidates, candidate_rank)
fail_if(anyNA(gtf_candidates$gencode19_gene_id), paste(
  "强候选无法映射到 GENCODE v19：",
  paste(gtf_candidates[is.na(gencode19_gene_id), tcga_gene_name], collapse = ", ")
))
fail_if(anyDuplicated(gtf_candidates$stable_gene_id) > 0L, "GENCODE v19 候选 gene 记录不唯一")

gene_body <- gtf_candidates[, .(
  chromosome,
  bed_start = pmax(0L, gtf_start - 1L),
  bed_end = gtf_end,
  tcga_gene_id,
  tcga_gene_name,
  gencode19_gene_name,
  region_type = "gene_body",
  strand,
  candidate_rank
)]
promoter <- gtf_candidates[, .(
  chromosome,
  bed_start = ifelse(
    strand == "+",
    pmax(0L, gtf_start - 2001L),
    pmax(0L, gtf_end - 500L)
  ),
  bed_end = ifelse(
    strand == "+",
    gtf_start + 499L,
    gtf_end + 2000L
  ),
  tcga_gene_id,
  tcga_gene_name,
  gencode19_gene_name,
  region_type = "promoter",
  strand,
  candidate_rank
)]
regions <- rbindlist(list(promoter, gene_body), use.names = TRUE)
regions[, chromosome_order := fifelse(
  grepl("^chr[0-9]+$", chromosome),
  as.integer(sub("^chr", "", chromosome)),
  1000L
)]
setorder(regions, chromosome_order, bed_start, bed_end, candidate_rank, region_type)
regions_bed <- file.path(stage_dir, "candidate_regions_hg19.bed")
fwrite(
  regions[, .(
    chromosome, bed_start, bed_end, tcga_gene_id, tcga_gene_name,
    gencode19_gene_name, region_type, strand
  )],
  regions_bed,
  sep = "\t",
  col.names = FALSE
)

message("[2/7] 以 bedtools 抽取 9 对 WGBS 候选区域并取配对共同 CpG")
wgbs_map <- fread(file.path(wgbs_root, "10_metadata", "sample_file_map.tsv"), showProgress = FALSE)
wgbs_map <- wgbs_map[analysis_eligibility == "paired_patient_analysis"]
fail_if(nrow(wgbs_map) != 18L || uniqueN(wgbs_map$paper_patient_id) != 9L,
        "WGBS 严格配对必须为 9 位患者/18 个文件")
fail_if(any(wgbs_map[, .N, by = .(paper_patient_id, condition)]$N != 1L),
        "WGBS 患者-条件映射不唯一")

extract_wgbs_sample <- function(i) {
  sample_row <- wgbs_map[i]
  source_file <- file.path(wgbs_root, sample_row$filename)
  fail_if(!file_exists(source_file), paste("WGBS 文件不存在：", source_file))
  out_file <- file.path(stage_dir, paste0("wgbs_extract_", sample_row$gsm, ".tsv"))
  err_file <- file.path(stage_dir, paste0("wgbs_extract_", sample_row$gsm, ".stderr.log"))
  status <- system2(
    bedtools,
    args = c("intersect", "-a", source_file, "-b", regions_bed, "-wa", "-wb"),
    stdout = out_file,
    stderr = err_file
  )
  fail_if(status != 0L, paste("bedtools WGBS 抽取失败：", sample_row$gsm, "；见", err_file))
  x <- fread(out_file, sep = "\t", header = FALSE, fill = TRUE, showProgress = FALSE)
  fail_if(ncol(x) != 12L, paste("WGBS bedtools 输出列数异常：", sample_row$gsm, ncol(x)))
  setnames(x, c(
    "cpg_chr", "cpg_start", "cpg_end", "methylation_percent",
    "region_chr", "region_start", "region_end", "tcga_gene_id",
    "tcga_gene_name", "gencode19_gene_name", "region_type", "strand"
  ))
  x[, methylation_percent := as.numeric(methylation_percent)]
  fail_if(any(!is.finite(x$methylation_percent) | x$methylation_percent < 0 | x$methylation_percent > 100),
          paste("WGBS 甲基化百分比越界：", sample_row$gsm))
  x[, `:=`(
    gsm = sample_row$gsm,
    paper_patient_id = sample_row$paper_patient_id,
    condition = sample_row$condition
  )]
  x[]
}

wgbs_extracted <- rbindlist(lapply(seq_len(nrow(wgbs_map)), extract_wgbs_sample), use.names = TRUE)
fail_if(!nrow(wgbs_extracted), "WGBS 候选区域未抽取到 CpG")
fail_if(anyDuplicated(wgbs_extracted[, .(
  gsm, tcga_gene_id, region_type, cpg_chr, cpg_start, cpg_end
)]) > 0L, "WGBS 单样本-候选区域存在重复 CpG 键")

wgbs_common <- dcast(
  wgbs_extracted,
  tcga_gene_id + tcga_gene_name + gencode19_gene_name + paper_patient_id +
    region_type + cpg_chr + cpg_start + cpg_end ~ condition,
  value.var = "methylation_percent"
)
fail_if(!all(c("N", "T") %chin% names(wgbs_common)), "WGBS 配对展宽缺少 N/T 列")
wgbs_common <- wgbs_common[is.finite(N) & is.finite(T)]
wgbs_common[, delta_beta := (T - N) / 100]
wgbs_pair_effects <- wgbs_common[, .(
  n_common_cpg = .N,
  median_normal_beta = median(N / 100),
  median_tumor_beta = median(T / 100),
  median_delta_beta = median(delta_beta),
  mean_delta_beta = mean(delta_beta),
  delta_beta_iqr = IQR(delta_beta),
  fraction_cpg_higher_in_tumor = mean(delta_beta > 0),
  fraction_cpg_lower_in_tumor = mean(delta_beta < 0)
), by = .(
  tcga_gene_id, tcga_gene_name, gencode19_gene_name,
  paper_patient_id, region_type
)]
wgbs_pair_effects[, coverage_gate := mapply(coverage_gate, region_type, n_common_cpg)]
wgbs_pair_effects <- merge(
  wgbs_pair_effects,
  candidates[, .(tcga_gene_id, candidate_rank)],
  by = "tcga_gene_id",
  all.x = TRUE,
  sort = FALSE
)
setorder(wgbs_pair_effects, candidate_rank, region_type, paper_patient_id)

wgbs_region_summary <- wgbs_pair_effects[, {
  evaluable <- coverage_gate != "technically_unavailable" & is.finite(median_delta_beta)
  effects <- median_delta_beta[evaluable]
  .(
    total_pairs = .N,
    primary_coverage_pairs = sum(coverage_gate == "primary_coverage"),
    conditional_coverage_pairs = sum(coverage_gate == "conditional_coverage"),
    technically_unavailable_pairs = sum(coverage_gate == "technically_unavailable"),
    evaluable_pairs = sum(evaluable),
    median_common_cpg = median(n_common_cpg),
    min_common_cpg = min(n_common_cpg),
    median_paired_delta_beta = if (length(effects)) median(effects) else NA_real_,
    paired_delta_beta_iqr = if (length(effects)) IQR(effects) else NA_real_,
    positive_pairs = sum(effects > 0),
    negative_pairs = sum(effects < 0),
    p_value = safe_wilcox_p(effects),
    region_evaluation_status = if (sum(evaluable) >= 5L) "evaluable" else if (sum(evaluable) >= 3L) "conditional_evaluable" else "technically_unavailable"
  )
}, by = .(tcga_gene_id, tcga_gene_name, gencode19_gene_name, candidate_rank, region_type)]
wgbs_region_summary[, q_candidate_by_region := p.adjust(p_value, method = "BH"), by = region_type]
wgbs_region_summary[, methylation_effect_class := region_effect_direction(median_paired_delta_beta)]
setorder(wgbs_region_summary, candidate_rank, region_type)

message("[3/7] 构建 10 对 RNA TPM 患者固定效应 limma")
rna_map <- fread(file.path(rna_root, "10_metadata", "sample_file_map.tsv"), showProgress = FALSE)
rna_map <- rna_map[analysis_eligibility == "paired_patient_analysis"]
fail_if(nrow(rna_map) != 20L || uniqueN(rna_map$paper_patient_id) != 10L,
        "RNA 配对必须为 10 位患者/20 个文件")
fail_if(any(rna_map[, .N, by = .(paper_patient_id, condition)]$N != 1L),
        "RNA 患者-条件映射不唯一")

read_rna_sample <- function(i) {
  sample_row <- rna_map[i]
  source_file <- file.path(rna_root, sample_row$filename)
  fail_if(!file_exists(source_file), paste("RNA 文件不存在：", source_file))
  x <- fread(
    cmd = paste("gzip -cd", shQuote(source_file)),
    select = c("gene_id", "TPM"),
    showProgress = FALSE
  )
  fail_if(anyDuplicated(x$gene_id) > 0L, paste("RNA gene_id 不唯一：", sample_row$gsm))
  fail_if(any(!is.finite(x$TPM) | x$TPM < 0), paste("RNA TPM 非有限或小于 0：", sample_row$gsm))
  x[, `:=`(
    gsm = sample_row$gsm,
    paper_patient_id = sample_row$paper_patient_id,
    condition = sample_row$condition
  )]
  x[]
}

rna_long <- rbindlist(lapply(seq_len(nrow(rna_map)), read_rna_sample), use.names = TRUE)
rna_matrix_dt <- dcast(rna_long, gene_id ~ gsm, value.var = "TPM")
sample_order <- rna_map$gsm
fail_if(!all(sample_order %chin% names(rna_matrix_dt)), "RNA 矩阵缺少样本")
rna_gene_ids <- rna_matrix_dt$gene_id
rna_tpm <- as.matrix(rna_matrix_dt[, ..sample_order])
rownames(rna_tpm) <- rna_gene_ids
storage.mode(rna_tpm) <- "double"

rna_design_data <- rna_map[match(sample_order, gsm), .(
  gsm,
  patient = factor(paper_patient_id),
  condition = factor(condition, levels = c("N", "T"))
)]
design <- model.matrix(~ patient + condition, data = rna_design_data)
fail_if(!"conditionT" %chin% colnames(design), "RNA limma 设计未生成 conditionT 系数")

fit_rna_pc <- function(pc) {
  fit <- lmFit(log2(rna_tpm + pc), design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  coefficient <- match("conditionT", colnames(fit$coefficients))
  data.table(
    gene_id = rownames(rna_tpm),
    log_fc = fit$coefficients[, coefficient],
    t_statistic = fit$t[, coefficient],
    p_value = fit$p.value[, coefficient],
    q_transcriptome = p.adjust(fit$p.value[, coefficient], method = "BH")
  )
}

rna_fit_pc01 <- fit_rna_pc(0.1)
rna_fit_pc1 <- fit_rna_pc(1)
gtf_candidates[, rna_source_gene_id := paste0(gencode19_gene_name, "_", gencode19_gene_name)]
missing_rna_ids <- setdiff(gtf_candidates$rna_source_gene_id, rna_gene_ids)
fail_if(length(missing_rna_ids) > 0L, paste(
  "GENCODE v19 候选无法在 Cao RNA gene_id 中定位：",
  paste(missing_rna_ids, collapse = ", ")
))

rna_candidate_long <- merge(
  rna_long[gene_id %chin% gtf_candidates$rna_source_gene_id],
  gtf_candidates[, .(
    tcga_gene_id, tcga_gene_name, gencode19_gene_name,
    rna_source_gene_id, candidate_rank
  )],
  by.x = "gene_id",
  by.y = "rna_source_gene_id",
  all.x = TRUE,
  sort = FALSE
)
rna_pair_effects <- dcast(
  rna_candidate_long,
  tcga_gene_id + tcga_gene_name + gencode19_gene_name + candidate_rank +
    paper_patient_id ~ condition,
  value.var = "TPM"
)
fail_if(anyNA(rna_pair_effects[, .(N, T)]), "RNA 候选配对存在缺失 TPM")
rna_pair_effects[, `:=`(
  delta_log2_tpm_pc01 = log2(T + 0.1) - log2(N + 0.1),
  delta_log2_tpm_pc1 = log2(T + 1) - log2(N + 1)
)]
setorder(rna_pair_effects, candidate_rank, paper_patient_id)

rna_summary <- rna_pair_effects[, .(
  paired_patients = .N,
  median_normal_tpm = median(N),
  median_tumor_tpm = median(T),
  median_delta_log2_tpm_pc01 = median(delta_log2_tpm_pc01),
  median_delta_log2_tpm_pc1 = median(delta_log2_tpm_pc1),
  positive_pairs_pc1 = sum(delta_log2_tpm_pc1 > 0),
  negative_pairs_pc1 = sum(delta_log2_tpm_pc1 < 0),
  paired_wilcox_p_pc1 = safe_wilcox_p(delta_log2_tpm_pc1)
), by = .(tcga_gene_id, tcga_gene_name, gencode19_gene_name, candidate_rank)]
rna_summary <- merge(
  rna_summary,
  gtf_candidates[, .(tcga_gene_id, rna_source_gene_id)],
  by = "tcga_gene_id",
  all.x = TRUE,
  sort = FALSE
)
rna_summary <- merge(
  rna_summary,
  rna_fit_pc01[, .(
    rna_source_gene_id = gene_id,
    limma_log_fc_pc01 = log_fc,
    limma_t_pc01 = t_statistic,
    limma_p_pc01 = p_value,
    limma_q_transcriptome_pc01 = q_transcriptome
  )],
  by = "rna_source_gene_id",
  all.x = TRUE,
  sort = FALSE
)
rna_summary <- merge(
  rna_summary,
  rna_fit_pc1[, .(
    rna_source_gene_id = gene_id,
    limma_log_fc_pc1 = log_fc,
    limma_t_pc1 = t_statistic,
    limma_p_pc1 = p_value,
    limma_q_transcriptome_pc1 = q_transcriptome
  )],
  by = "rna_source_gene_id",
  all.x = TRUE,
  sort = FALSE
)
rna_summary[, `:=`(
  limma_q_candidate_pc01 = p.adjust(limma_p_pc01, method = "BH"),
  limma_q_candidate_pc1 = p.adjust(limma_p_pc1, method = "BH"),
  paired_wilcox_q_candidate_pc1 = p.adjust(paired_wilcox_p_pc1, method = "BH")
)]
rna_summary[, pseudocount_direction_stable :=
              sign(limma_log_fc_pc01) == sign(limma_log_fc_pc1)]
rna_summary[, rna_effect_class := rna_effect_direction(
  limma_log_fc_pc1,
  pseudocount_direction_stable
)]
setorder(rna_summary, candidate_rank)

message("[4/7] 按 protein group 传播 GN 并将 P7 批次折叠为一位患者")
protein_id <- fread(
  file.path(protein_root, "20_reusable", "proteomics_protein_identification.tsv"),
  showProgress = FALSE
)
protein_ratio <- fread(
  file.path(protein_root, "20_reusable", "proteomics_ratio_long.tsv"),
  na.strings = c("", "NA"),
  showProgress = FALSE
)
protein_id[, source_gene_symbol := extract_gn(description)]
candidate_aliases <- unique(rbindlist(list(
  gtf_candidates[, .(tcga_gene_id, tcga_gene_name, candidate_rank, source_gene_symbol = tcga_gene_name)],
  gtf_candidates[, .(tcga_gene_id, tcga_gene_name, candidate_rank, source_gene_symbol = gencode19_gene_name)]
)))
candidate_aliases <- candidate_aliases[!is.na(source_gene_symbol) & nzchar(source_gene_symbol)]

protein_group_candidates <- unique(merge(
  protein_id[!is.na(source_gene_symbol), .(
    batch, sheet_name, group_id, source_gene_symbol
  )],
  candidate_aliases,
  by = "source_gene_symbol",
  all = FALSE,
  allow.cartesian = TRUE
))

protein_group_rows <- merge(
  protein_id[, .(
    batch, sheet_name, source_row, group_id, hit_number,
    description, uniprot_swissprot_accession
  )],
  protein_group_candidates[, .(
    batch, sheet_name, group_id, tcga_gene_id, tcga_gene_name,
    candidate_rank, source_gene_symbol
  )],
  by = c("batch", "sheet_name", "group_id"),
  all = FALSE,
  allow.cartesian = TRUE
)
protein_ratio_mapped <- merge(
  protein_ratio,
  protein_group_rows,
  by = c("batch", "sheet_name", "source_row"),
  all = FALSE,
  allow.cartesian = TRUE
)

protein_group_patient_batch <- protein_ratio_mapped[, {
  numeric_values <- log2_tumor_vs_normal[is.finite(log2_tumor_vs_normal)]
  .(
    n_source_ratio_rows = .N,
    n_numeric_source_rows = length(numeric_values),
    log2_tumor_vs_normal = if (length(numeric_values)) median(numeric_values) else NA_real_,
    raw_ratio_values = paste(sort(unique(as.character(normal_vs_tumor_ratio))), collapse = ";"),
    source_rows = paste(sort(unique(source_row)), collapse = ";"),
    uniprot_swissprot_accessions = paste(sort(unique(
      uniprot_swissprot_accession[!is.na(uniprot_swissprot_accession) & nzchar(uniprot_swissprot_accession)]
    )), collapse = ";")
  )
}, by = .(
  tcga_gene_id, tcga_gene_name, candidate_rank, source_gene_symbol,
  batch, sheet_name, group_id, paper_patient_id
)]

protein_batch_patient <- protein_group_patient_batch[, {
  numeric_values <- log2_tumor_vs_normal[is.finite(log2_tumor_vs_normal)]
  .(
    protein_groups = uniqueN(group_id),
    numeric_protein_groups = sum(is.finite(log2_tumor_vs_normal)),
    log2_tumor_vs_normal = if (length(numeric_values)) median(numeric_values) else NA_real_,
    raw_ratio_values = paste(sort(unique(raw_ratio_values)), collapse = " | "),
    group_ids = paste(sort(unique(group_id)), collapse = ";"),
    uniprot_swissprot_accessions = paste(sort(unique(uniprot_swissprot_accessions[nzchar(uniprot_swissprot_accessions)])), collapse = ";")
  )
}, by = .(
  tcga_gene_id, tcga_gene_name, candidate_rank,
  batch, sheet_name, paper_patient_id
)]

protein_pair_effects <- protein_batch_patient[, {
  numeric_values <- log2_tumor_vs_normal[is.finite(log2_tumor_vs_normal)]
  .(
    identified_batches = .N,
    quantified_batches = length(numeric_values),
    log2_tumor_vs_normal = if (length(numeric_values)) median(numeric_values) else NA_real_,
    raw_batch_log2_tumor_vs_normal = paste(
      paste0("batch", batch, ":", fifelse(is.finite(log2_tumor_vs_normal), format(round(log2_tumor_vs_normal, 6), trim = TRUE), "NA")),
      collapse = ";"
    ),
    group_ids = paste(sort(unique(group_ids)), collapse = " | "),
    uniprot_swissprot_accessions = paste(sort(unique(uniprot_swissprot_accessions[nzchar(uniprot_swissprot_accessions)])), collapse = ";"),
    measurement_status = if (length(numeric_values)) "quantified" else "identified_not_quantified"
  )
}, by = .(tcga_gene_id, tcga_gene_name, candidate_rank, paper_patient_id)]
setorder(protein_pair_effects, candidate_rank, paper_patient_id)

protein_group_summary <- protein_group_candidates[, .(
  identified_protein_groups = uniqueN(paste(batch, sheet_name, group_id)),
  identified_batches = uniqueN(batch)
), by = .(tcga_gene_id, tcga_gene_name, candidate_rank)]
protein_quant_summary <- protein_pair_effects[, {
  values <- log2_tumor_vs_normal[is.finite(log2_tumor_vs_normal)]
  .(
    patient_records = .N,
    quantified_patients = length(values),
    numeric_batch_values = sum(quantified_batches),
    p7_quantified_batches = sum(quantified_batches[paper_patient_id == "P7"]),
    median_patient_log2_tumor_vs_normal = if (length(values)) median(values) else NA_real_,
    protein_effect_iqr = if (length(values)) IQR(values) else NA_real_,
    positive_patients = sum(values > 0),
    negative_patients = sum(values < 0),
    p_value = safe_wilcox_p(values)
  )
}, by = .(tcga_gene_id, tcga_gene_name, candidate_rank)]
protein_summary <- merge(
  candidates[, .(tcga_gene_id, tcga_gene_name, candidate_rank)],
  protein_group_summary,
  by = c("tcga_gene_id", "tcga_gene_name", "candidate_rank"),
  all.x = TRUE,
  sort = FALSE
)
protein_summary <- merge(
  protein_summary,
  protein_quant_summary,
  by = c("tcga_gene_id", "tcga_gene_name", "candidate_rank"),
  all.x = TRUE,
  sort = FALSE
)
for (column in c(
  "identified_protein_groups", "identified_batches", "patient_records",
  "quantified_patients", "numeric_batch_values", "p7_quantified_batches",
  "positive_patients", "negative_patients"
)) {
  set(protein_summary, which(is.na(protein_summary[[column]])), column, 0L)
}
protein_summary[, protein_measurement_status := fifelse(
  identified_protein_groups == 0L,
  "not_identified",
  fifelse(quantified_patients == 0L, "identified_not_quantified", "quantified")
)]
protein_summary[, q_candidate := p.adjust(p_value, method = "BH")]
protein_summary[, protein_effect_class := mapply(
  protein_effect_direction,
  protein_measurement_status,
  quantified_patients,
  median_patient_log2_tumor_vs_normal
)]
setorder(protein_summary, candidate_rank)

message("[5/7] 整合患者级 WGBS-RNA-蛋白效应并保留可信反向")
wgbs_patient_wide <- dcast(
  wgbs_pair_effects,
  tcga_gene_id + tcga_gene_name + candidate_rank + paper_patient_id ~ region_type,
  value.var = c("median_delta_beta", "n_common_cpg", "coverage_gate")
)
cross_layer_patient <- CJ(
  tcga_gene_id = candidates$tcga_gene_id,
  paper_patient_id = sort(unique(rna_map$paper_patient_id)),
  unique = TRUE
)
cross_layer_patient <- merge(
  cross_layer_patient,
  candidates[, .(tcga_gene_id, tcga_gene_name, candidate_rank)],
  by = "tcga_gene_id",
  all.x = TRUE,
  sort = FALSE
)
cross_layer_patient <- merge(
  cross_layer_patient,
  rna_pair_effects[, .(
    tcga_gene_id, paper_patient_id,
    normal_tpm = N,
    tumor_tpm = T,
    rna_delta_log2_tpm_pc01 = delta_log2_tpm_pc01,
    rna_delta_log2_tpm_pc1 = delta_log2_tpm_pc1
  )],
  by = c("tcga_gene_id", "paper_patient_id"),
  all.x = TRUE,
  sort = FALSE
)
cross_layer_patient <- merge(
  cross_layer_patient,
  wgbs_patient_wide,
  by = c("tcga_gene_id", "tcga_gene_name", "candidate_rank", "paper_patient_id"),
  all.x = TRUE,
  sort = FALSE
)
cross_layer_patient <- merge(
  cross_layer_patient,
  protein_pair_effects[, .(
    tcga_gene_id, paper_patient_id,
    protein_log2_tumor_vs_normal = log2_tumor_vs_normal,
    protein_identified_batches = identified_batches,
    protein_quantified_batches = quantified_batches,
    protein_patient_measurement_status = measurement_status
  )],
  by = c("tcga_gene_id", "paper_patient_id"),
  all.x = TRUE,
  sort = FALSE
)
setorder(cross_layer_patient, candidate_rank, paper_patient_id)

cross_correlations <- cross_layer_patient[, {
  pr <- safe_spearman(median_delta_beta_promoter, rna_delta_log2_tpm_pc1)
  rp <- safe_spearman(rna_delta_log2_tpm_pc1, protein_log2_tumor_vs_normal)
  mp <- safe_spearman(median_delta_beta_promoter, protein_log2_tumor_vs_normal)
  .(
    promoter_rna_n = pr$n,
    promoter_rna_rho = pr$rho,
    promoter_rna_p = pr$p_value,
    rna_protein_n = rp$n,
    rna_protein_rho = rp$rho,
    rna_protein_p = rp$p_value,
    promoter_protein_n = mp$n,
    promoter_protein_rho = mp$rho,
    promoter_protein_p = mp$p_value
  )
}, by = .(tcga_gene_id, tcga_gene_name, candidate_rank)]

promoter_summary <- copy(wgbs_region_summary[region_type == "promoter"])
setnames(
  promoter_summary,
  setdiff(names(promoter_summary), c("tcga_gene_id", "tcga_gene_name", "gencode19_gene_name", "candidate_rank", "region_type")),
  paste0("promoter_", setdiff(names(promoter_summary), c("tcga_gene_id", "tcga_gene_name", "gencode19_gene_name", "candidate_rank", "region_type")))
)
body_summary <- copy(wgbs_region_summary[region_type == "gene_body"])
setnames(
  body_summary,
  setdiff(names(body_summary), c("tcga_gene_id", "tcga_gene_name", "gencode19_gene_name", "candidate_rank", "region_type")),
  paste0("gene_body_", setdiff(names(body_summary), c("tcga_gene_id", "tcga_gene_name", "gencode19_gene_name", "candidate_rank", "region_type")))
)
promoter_summary[, region_type := NULL]
body_summary[, region_type := NULL]

cross_layer_summary <- merge(
  candidates,
  promoter_summary,
  by = c("tcga_gene_id", "tcga_gene_name", "candidate_rank"),
  all.x = TRUE,
  sort = FALSE
)
cross_layer_summary <- merge(
  cross_layer_summary,
  body_summary[, setdiff(names(body_summary), "gencode19_gene_name"), with = FALSE],
  by = c("tcga_gene_id", "tcga_gene_name", "candidate_rank"),
  all.x = TRUE,
  sort = FALSE
)
cross_layer_summary <- merge(
  cross_layer_summary,
  rna_summary,
  by = c("tcga_gene_id", "tcga_gene_name", "candidate_rank", "gencode19_gene_name"),
  all.x = TRUE,
  sort = FALSE
)
cross_layer_summary <- merge(
  cross_layer_summary,
  protein_summary,
  by = c("tcga_gene_id", "tcga_gene_name", "candidate_rank"),
  all.x = TRUE,
  sort = FALSE
)
cross_layer_summary <- merge(
  cross_layer_summary,
  cross_correlations,
  by = c("tcga_gene_id", "tcga_gene_name", "candidate_rank"),
  all.x = TRUE,
  sort = FALSE
)

cross_layer_summary[, promoter_rna_relationship := {
  methylation <- promoter_methylation_effect_class
  rna <- rna_effect_class
  fifelse(
    promoter_region_evaluation_status == "technically_unavailable",
    "methylation_technically_unavailable",
    fifelse(
      methylation == "small_or_mixed_shift" | rna == "small_or_mixed_shift",
      "weak_or_mixed",
      fifelse(
        rna == "pseudocount_sensitive",
        "rna_pseudocount_sensitive",
        fifelse(
          (methylation == "hypermethylated_in_tumor" & rna == "lower_in_tumor") |
            (methylation == "hypomethylated_in_tumor" & rna == "higher_in_tumor"),
          "directionally_coherent_inverse",
          fifelse(
            (methylation == "hypermethylated_in_tumor" & rna == "higher_in_tumor") |
              (methylation == "hypomethylated_in_tumor" & rna == "lower_in_tumor"),
            "credible_reverse_same_direction_retained",
            "not_evaluable"
          )
        )
      )
    )
  )
}]

cross_layer_summary[, rna_protein_relationship := {
  rna <- rna_effect_class
  protein <- protein_effect_class
  fifelse(
    protein_measurement_status == "not_identified",
    "protein_missing_not_negative",
    fifelse(
      protein_measurement_status == "identified_not_quantified",
      "protein_identified_not_quantified",
      fifelse(
        protein %chin% c("quantified_but_sparse", "small_or_mixed_shift") |
          rna %chin% c("small_or_mixed_shift", "pseudocount_sensitive"),
        "weak_sparse_or_mixed",
        fifelse(
          (rna == "higher_in_tumor" & protein == "higher_in_tumor") |
            (rna == "lower_in_tumor" & protein == "lower_in_tumor"),
          "directionally_coherent",
          fifelse(
            (rna == "higher_in_tumor" & protein == "lower_in_tumor") |
              (rna == "lower_in_tumor" & protein == "higher_in_tumor"),
            "credible_reverse_retained",
            "not_evaluable"
          )
        )
      )
    )
  )
}]

cross_layer_summary[, cross_layer_class := fifelse(
  grepl("credible_reverse", promoter_rna_relationship) |
    grepl("credible_reverse", rna_protein_relationship),
  "credible_reverse_retained",
  fifelse(
    promoter_rna_relationship == "directionally_coherent_inverse" &
      rna_protein_relationship == "directionally_coherent",
    "three_layer_directional_hypothesis",
    fifelse(
      promoter_rna_relationship == "directionally_coherent_inverse" &
        rna_protein_relationship %chin% c("protein_missing_not_negative", "protein_identified_not_quantified"),
      "two_layer_directional_hypothesis_protein_missing_not_negative",
      fifelse(
        promoter_rna_relationship == "directionally_coherent_inverse",
        "two_layer_directional_hypothesis",
        "context_dependent_or_weak"
      )
    )
  )
)]
cross_layer_summary[, fdr_supported_layer_count :=
                      as.integer(is.finite(promoter_q_candidate_by_region) & promoter_q_candidate_by_region <= 0.05) +
                      as.integer(is.finite(limma_q_candidate_pc1) & limma_q_candidate_pc1 <= 0.05) +
                      as.integer(is.finite(q_candidate) & q_candidate <= 0.05)]
cross_layer_summary[, statistical_support_class := fifelse(
  fdr_supported_layer_count >= 2L,
  "multi_layer_candidate_fdr_support",
  fifelse(
    fdr_supported_layer_count == 1L,
    "single_layer_candidate_fdr_support",
    "directional_only_no_candidate_layer_fdr"
  )
)]
cross_layer_summary[, available_layer_count :=
                      1L +
                      as.integer(promoter_region_evaluation_status != "technically_unavailable") +
                      as.integer(protein_measurement_status == "quantified")]
cross_layer_summary[, conclusion_ceiling := "公共数据假设生成与跨层校准；不证明因果机制或治疗靶点"]
cross_layer_summary[, required_follow_up := fifelse(
  cross_layer_class == "credible_reverse_retained",
  "在独立 ESCC 配对队列中复核反向，并以蛋白/功能实验区分技术与生物学脱耦",
  fifelse(
    protein_measurement_status != "quantified",
    "补充候选蛋白定量与独立 ESCC 配对队列；缺失蛋白层不视为阴性",
    "在独立 ESCC 队列与患者级重抽样中复核，再进入机制实验"
  )
)]
setorder(cross_layer_summary, candidate_rank)

message("[6/7] 生成正式结果、摘要、artifact manifest 和单一 QA 记录")
stage_tsv <- function(object, filename) {
  fwrite(object, file.path(stage_dir, filename), sep = "\t", na = "")
}

formal_tables <- list(
  cao2020_wgbs_candidate_pair_effects.tsv = wgbs_pair_effects[, candidate_rank := NULL][],
  cao2020_wgbs_candidate_region_summary.tsv = wgbs_region_summary[, candidate_rank := NULL][],
  cao2020_rna_candidate_pair_effects.tsv = rna_pair_effects[, candidate_rank := NULL][],
  cao2020_rna_candidate_summary.tsv = rna_summary[, candidate_rank := NULL][],
  cao2020_proteomics_candidate_pair_effects.tsv = protein_pair_effects[, candidate_rank := NULL][],
  cao2020_proteomics_candidate_summary.tsv = protein_summary[, candidate_rank := NULL][],
  cao2020_cross_layer_patient_effects.tsv = cross_layer_patient[, candidate_rank := NULL][],
  cao2020_cross_layer_candidate_summary.tsv = cross_layer_summary[, candidate_rank := NULL][]
)
invisible(Map(stage_tsv, formal_tables, names(formal_tables)))

format_number <- function(x, digits = 3L) {
  ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "NA")
}

axis_lines <- cross_layer_summary[, paste0(
  "- `", tcga_gene_name, "`：", cross_layer_class,
  "；promoter median Δβ=", format_number(promoter_median_paired_delta_beta),
  "；RNA limma logFC(T−N, pc=1)=", format_number(limma_log_fc_pc1),
  "；蛋白状态=", protein_measurement_status,
  fifelse(is.finite(median_patient_log2_tumor_vs_normal),
          paste0("，median log2(T/N)=", format_number(median_patient_log2_tumor_vs_normal)), ""),
  "；统计层级=", statistical_support_class,
  "。"
)]
reverse_genes <- cross_layer_summary[cross_layer_class == "credible_reverse_retained", tcga_gene_name]
protein_missing_genes <- cross_layer_summary[protein_measurement_status != "quantified", tcga_gene_name]
summary_lines <- c(
  "# Cao 2020 同患者跨层候选轴校准摘要",
  "",
  "## 结论边界",
  "",
  "本分析是公共数据假设生成与跨层校准，不证明因果机制、治疗靶点或已完成湿实验验证。",
  "",
  "## 分析母集",
  "",
  paste0("- TCGA `strong_patient_level_candidate`动态读取：", nrow(candidates), " 个，未固定旧候选名单。"),
  "- WGBS：9 对患者，GENCODE v19/GRCh37，promoter 定义为 TSS 上游 2 kb/下游 500 bp；每位患者先取 N/T 共同 CpG 再计算 Δβ。",
  "- RNA：10 对 RSEM TPM，以患者固定效应 limma 估计 T−N，并比较 pseudocount 0.1 与 1。",
  "- 蛋白：按 batch+sheet+group_id 传播 GN；P7 先折叠批次为一位患者；缺失层不视为阴性。",
  paste0("- 方向性分类不等于统计显著；本轮至少一个候选层 q≤0.05 的候选为 ", sum(cross_layer_summary$fdr_supported_layer_count >= 1L), "/", nrow(cross_layer_summary), "。"),
  "",
  "## 候选轴状态",
  "",
  axis_lines,
  "",
  "## 冲突与缺失",
  "",
  paste0("- 保留的可信反向：", if (length(reverse_genes)) paste(reverse_genes, collapse = "、") else "无", "。"),
  paste0("- 蛋白未定量或未鉴定（不当作阴性）：", if (length(protein_missing_genes)) paste(protein_missing_genes, collapse = "、") else "无", "。"),
  "- promoter–RNA 反向符合只表示方向性校准；gene-body 甲基化、蛋白脱耦与小样本反向必须单独保留。",
  "",
  "## 下一步",
  "",
  "优先复核具有方向性收敛或稳定反向的候选，但仍需独立 ESCC 配对队列、患者级重抽样、蛋白定量和机制实验后才能升级结论。"
)
summary_filename <- "cao2020_cross_layer_summary.md"
writeLines(summary_lines, file.path(stage_dir, summary_filename), useBytes = TRUE)

formal_filenames <- c(names(formal_tables), summary_filename)
artifact_manifest <- rbindlist(lapply(formal_filenames, function(filename) {
  path <- file.path(stage_dir, filename)
  data.table(
    artifact = filename,
    relative_path = file.path("results", filename),
    file_size_bytes = file_size(path),
    sha256 = digest(path, algo = "sha256", file = TRUE),
    generated_date = as.character(Sys.Date()),
    generation_script = "scripts/15_analyze_cao2020_cross_layer_axes.R",
    source_object = "TCGA strong candidates; Cao 2020 GSE149608/GSE149609/Supplementary Data 6; GENCODE v19",
    status = "verified"
  )
}))
manifest_filename <- "cao2020_cross_layer_artifact_manifest.tsv"
stage_tsv(artifact_manifest, manifest_filename)

bedtools_version <- system2(bedtools, "--version", stdout = TRUE, stderr = TRUE)
qa_filename <- "cao2020_cross_layer_analysis_20260711.md"
qa_lines <- c(
  "# Cao 2020 跨层候选轴分析 QA 记录",
  "",
  "此文件是运行时历史证据，不是项目当前状态源。",
  "",
  paste0("- 运行日期：", Sys.Date(), "。"),
  paste0("- R：", R.version.string, "。"),
  paste0("- data.table：", as.character(packageVersion("data.table")), "。"),
  paste0("- limma：", as.character(packageVersion("limma")), "。"),
  paste0("- bedtools：", paste(bedtools_version, collapse = " "), "。"),
  paste0("- 动态 TCGA 强候选：", nrow(candidates), " 个。"),
  "- GENCODE v19/GRCh37 候选映射：全部唯一通过。",
  paste0("- WGBS：", uniqueN(wgbs_pair_effects$paper_patient_id), " 对患者，", nrow(wgbs_pair_effects), " 个患者-候选-区域效应。"),
  paste0("- RNA：", uniqueN(rna_pair_effects$paper_patient_id), " 对患者，pc=0.1/1 方向稳定候选 ", sum(rna_summary$pseudocount_direction_stable), "/", nrow(rna_summary), "。"),
  paste0("- 蛋白：定量候选 ", sum(protein_summary$protein_measurement_status == "quantified"),
         "，仅鉴定未定量 ", sum(protein_summary$protein_measurement_status == "identified_not_quantified"),
         "，未鉴定 ", sum(protein_summary$protein_measurement_status == "not_identified"), "。"),
  paste0("- P7 最大定量批次数：", max(protein_pair_effects$quantified_batches, na.rm = TRUE), "；患者级输出中只保留一行。"),
  paste0("- 可信反向保留：", sum(cross_layer_summary$cross_layer_class == "credible_reverse_retained"), " 个。"),
  "- 缺失 WGBS/蛋白层均保留为技术不可评估或缺失，没有写成阴性。",
  "- 公共数据结论上限已固定为假设生成/校准，没有写成因果机制。",
  paste0("- 正式 artifact：", nrow(artifact_manifest), " 个，已生成大小与 SHA256。"),
  "- 成功发布后由脚本显式删除 stage_dir；失败时保留供追溯。",
  "- ResearchDataHub canonical 原始资料、项目原始数据和投稿包均未修改或删除。"
)
writeLines(qa_lines, file.path(stage_dir, qa_filename), useBytes = TRUE)

message("[7/7] 重读、验证并原子化发布")
for (filename in names(formal_tables)) {
  reread <- fread(file.path(stage_dir, filename), showProgress = FALSE)
  fail_if(!nrow(reread), paste("正式 TSV 为空：", filename))
}
fail_if(nrow(fread(file.path(stage_dir, "cao2020_cross_layer_candidate_summary.tsv"))) != nrow(candidates),
        "跨层候选摘要行数与动态 TCGA 强候选数不一致")
fail_if(nrow(fread(file.path(stage_dir, "cao2020_rna_candidate_summary.tsv"))) != nrow(candidates),
        "RNA 候选摘要行数异常")
fail_if(nrow(fread(file.path(stage_dir, "cao2020_proteomics_candidate_summary.tsv"))) != nrow(candidates),
        "蛋白候选摘要行数异常")
fail_if(any(artifact_manifest$file_size_bytes <= 0), "artifact 存在零字节文件")

publish_files <- c(formal_filenames, manifest_filename)
for (filename in publish_files) {
  file_copy(file.path(stage_dir, filename), file.path(results_dir, filename), overwrite = TRUE)
}
file_copy(
  file.path(stage_dir, qa_filename),
  file.path(work_checks_dir, qa_filename),
  overwrite = TRUE
)

for (i in seq_len(nrow(artifact_manifest))) {
  published <- file.path(project_root, artifact_manifest$relative_path[i])
  fail_if(!file_exists(published), paste("发布 artifact 不存在：", published))
  fail_if(as.numeric(file_size(published)) != as.numeric(artifact_manifest$file_size_bytes[i]),
          paste("发布 artifact 大小不一致：", published))
  fail_if(digest(published, algo = "sha256", file = TRUE) != artifact_manifest$sha256[i],
          paste("发布 artifact SHA256 不一致：", published))
}

dir_delete(stage_dir)
message("完成：Cao 2020 跨层候选轴分析已发布，stage_dir 已显式删除。")
