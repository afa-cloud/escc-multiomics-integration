#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(edgeR)
  library(fs)
  library(jsonlite)
  library(Matrix)
  library(MultiAssayExperiment)
  library(S4Vectors)
  library(SummarizedExperiment)
})

options(stringsAsFactors = FALSE)

# TCGA-ESCC GDC DR45 开放多组学患者级矩阵构建。
#
# 只读取 ResearchDataHub 中已校验的 canonical 文件，不联网、不下载、
# 不重算 653 个源文件的全量 SHA256。脚本只完成患者映射、矩阵化、
# 结构/质量门禁和 MultiAssayExperiment 构建，不做聚类或生物学建模。

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required_project_files <- c(
  "PROJECT_INDEX.md",
  "data/datasets.tsv",
  "results/tcga_escc_patient_whitelist.tsv"
)
missing_project_files <- required_project_files[
  !file_exists(file.path(project_root, required_project_files))
]
if (length(missing_project_files)) {
  stop(
    "当前工作目录不是完整项目根目录，缺少：",
    paste(missing_project_files, collapse = ", ")
  )
}

data_root <- Sys.getenv("RESEARCH_DATA_ROOT", path.expand("~/ResearchDataHub"))
if (!dir_exists(data_root)) {
  stop("ResearchDataHub 不存在或未挂载：", data_root)
}

results_dir <- file.path(project_root, "results")
work_intermediate_dir <- file.path(project_root, "_work", "intermediate")
work_checks_dir <- file.path(project_root, "_work", "checks")
dir_create(c(results_dir, work_intermediate_dir, work_checks_dir), recurse = TRUE)

target_dataset_key <- "GDC_TCGA_ESCA_open_multiomics_gdc_dr45_retrieved_20260711"
target_clinical_key <- "GDC_TCGA_ESCA_clinical_gdc_dr45_retrieved_20260711"
project_datasets <- fread(
  file.path(project_root, "data", "datasets.tsv"),
  na.strings = c("", "NA"),
  showProgress = FALSE
)
if (project_datasets[project_datasets[["dataset_key"]] == target_dataset_key, .N] != 1L) {
  stop("data/datasets.tsv 中目标开放多组学 dataset_key 必须唯一")
}
if (project_datasets[project_datasets[["dataset_key"]] == target_clinical_key, .N] != 1L) {
  stop("data/datasets.tsv 中目标临床 dataset_key 必须唯一")
}

canonical_root <- project_datasets[
  project_datasets[["dataset_key"]] == target_dataset_key,
  central_path
]
clinical_root <- project_datasets[
  project_datasets[["dataset_key"]] == target_clinical_key,
  central_path
]
if (!dir_exists(canonical_root) || !dir_exists(clinical_root)) {
  stop("TCGA canonical 或临床 canonical 路径不可访问")
}

manifest_path <- file.path(
  canonical_root,
  "10_metadata",
  "tcga_escc_open_multiomics_file_manifest.tsv"
)
entity_map_path <- file.path(
  canonical_root,
  "10_metadata",
  "tcga_escc_open_multiomics_entity_map.tsv"
)
availability_path <- file.path(
  canonical_root,
  "20_reusable",
  "tcga_escc_open_multiomics_patient_availability.tsv"
)
intersections_path <- file.path(
  canonical_root,
  "20_reusable",
  "tcga_escc_open_multiomics_intersections.tsv"
)
clinical_json_path <- file.path(
  clinical_root,
  "10_metadata",
  "tcga_esca_gdc_cases_dr45.json"
)

required_inputs <- c(
  manifest_path,
  entity_map_path,
  availability_path,
  intersections_path,
  clinical_json_path,
  file.path(canonical_root, "DATASET.md"),
  file.path(canonical_root, "90_manifests", "MANIFEST.tsv")
)
if (any(!file_exists(required_inputs))) {
  stop(
    "TCGA canonical 输入不完整：",
    paste(required_inputs[!file_exists(required_inputs)], collapse = ", ")
  )
}

fail_if <- function(condition, message) {
  if (isTRUE(condition)) stop(message, call. = FALSE)
}

as_flag <- function(x) {
  if (is.logical(x)) return(x)
  tolower(as.character(x)) %chin% c("true", "t", "1", "yes")
}

strip_ensembl_version <- function(x) sub("[.][0-9]+$", "", x)

align_values <- function(reference_ids, observed_ids, values, label) {
  observed_ids <- as.character(observed_ids)
  fail_if(anyDuplicated(observed_ids) > 0L, paste0(label, " feature ID 重复"))
  if (identical(reference_ids, observed_ids)) return(values)
  idx <- match(reference_ids, observed_ids)
  fail_if(
    length(observed_ids) != length(reference_ids) || anyNA(idx),
    paste0(label, " feature 集合与首文件不一致")
  )
  values[idx]
}

make_measurement_coldata <- function(file_selection, layer, patients) {
  dt <- file_selection[
    layer_id == layer & selected_for_patient_matrix,
    .(
      patient_id,
      file_id,
      file_name,
      sample_submitter_id,
      portion_submitter_id,
      analyte_submitter_id,
      aliquot_submitter_id,
      matched_normal_submitter_id,
      target_relative_path
    )
  ]
  setkey(dt, patient_id)
  dt <- dt[patients]
  fail_if(anyNA(dt$file_id), paste0(layer, " measurement colData 存在缺失文件"))
  df <- as.data.frame(dt)
  rownames(df) <- df$patient_id
  S4Vectors::DataFrame(df, check.names = FALSE)
}

robust_outlier <- function(x, cutoff = 3) {
  x <- as.numeric(x)
  med <- median(x, na.rm = TRUE)
  scale <- mad(x, center = med, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(scale) || scale == 0) return(rep(FALSE, length(x)))
  abs(x - med) / scale > cutoff
}

qc_rows <- data.table(
  scope = character(),
  layer = character(),
  entity_id = character(),
  metric = character(),
  observed = character(),
  expected = character(),
  gate_type = character(),
  status = character(),
  action = character(),
  note = character()
)

add_qc <- function(
    scope,
    layer,
    entity_id,
    metric,
    observed,
    expected,
    gate_type,
    status,
    action,
    note = "") {
  qc_rows <<- rbind(
    qc_rows,
    data.table(
      scope = as.character(scope),
      layer = as.character(layer),
      entity_id = as.character(entity_id),
      metric = as.character(metric),
      observed = as.character(observed),
      expected = as.character(expected),
      gate_type = as.character(gate_type),
      status = as.character(status),
      action = as.character(action),
      note = as.character(note)
    )
  )
}

add_soft_metric <- function(layer, patient_ids, metric, values, note = "") {
  flags <- robust_outlier(values)
  for (i in seq_along(patient_ids)) {
    add_qc(
      scope = "sample",
      layer = layer,
      entity_id = patient_ids[[i]],
      metric = metric,
      observed = format(values[[i]], scientific = FALSE, trim = TRUE),
      expected = "按队列分布复核；3 MAD 仅作软标记",
      gate_type = "soft",
      status = if (isTRUE(flags[[i]])) "soft_flag" else "pass",
      action = if (isTRUE(flags[[i]])) {
        "保留样本；进入层内 QC、PCA 和敏感性分析"
      } else {
        "保留"
      },
      note = note
    )
  }
}

manifest <- fread(manifest_path, na.strings = c("", "NA"), showProgress = FALSE)
entity_map <- fread(entity_map_path, na.strings = c("", "NA"), showProgress = FALSE)
availability <- fread(availability_path, na.strings = c("", "NA"), showProgress = FALSE)
intersections <- fread(intersections_path, na.strings = c("", "NA"), showProgress = FALSE)
whitelist <- fread(
  file.path(project_root, "results", "tcga_escc_patient_whitelist.tsv"),
  na.strings = c("", "NA"),
  showProgress = FALSE
)

expected_file_counts <- c(
  rna_star = 95L,
  mirna_bcgsc = 96L,
  methylation_hm450 = 96L,
  masked_maf = 96L,
  cnv_gene_ascat2 = 96L,
  cnv_segment_ascat2 = 96L,
  rppa = 78L
)
actual_file_counts <- manifest[, .N, by = layer_id]
actual_file_counts <- setNames(actual_file_counts$N, actual_file_counts$layer_id)
fail_if(
  !identical(as.integer(actual_file_counts[names(expected_file_counts)]),
             as.integer(expected_file_counts)),
  "七层文件数不符合 95/96/96/96/96/96/78 固定门禁"
)
fail_if(nrow(manifest) != 653L, "开放多组学 manifest 必须为 653 个文件")
fail_if(any(manifest$download_status != "verified"), "存在未 verified 的开放多组学文件")
fail_if(any(!as_flag(entity_map$hierarchy_match)), "entity map 存在层级映射失败")
fail_if(nrow(whitelist) != 96L || uniqueN(whitelist$patient) != 96L,
        "ESCC 患者白名单必须为 96 位唯一患者")

source_paths <- file.path(canonical_root, manifest$target_relative_path)
fail_if(any(!file_exists(source_paths)), "manifest 中存在缺失源文件")
observed_sizes <- as.numeric(file_info(source_paths)$size)
fail_if(any(observed_sizes != as.numeric(manifest$local_size_bytes)),
        "源文件大小与 canonical manifest 不一致")

single_layers <- c("rna_star", "mirna_bcgsc", "methylation_hm450", "rppa")
pair_layers <- c("masked_maf", "cnv_gene_ascat2", "cnv_segment_ascat2")
single_checks <- entity_map[layer_id %chin% single_layers, .(
  primary_n = sum(analysis_role == "primary_tumor_measurement")
), by = .(layer_id, file_id)]
pair_checks <- entity_map[layer_id %chin% pair_layers, .(
  tumor_n = sum(analysis_role == "tumor_associated_input"),
  normal_n = sum(analysis_role == "matched_normal_associated_input")
), by = .(layer_id, file_id)]
fail_if(any(single_checks$primary_n != 1L),
        "RNA/miRNA/HM450/RPPA 必须每文件映射一个原发肿瘤实体")
fail_if(any(pair_checks$tumor_n != 1L | pair_checks$normal_n != 1L),
        "MAF/CNV 必须每文件映射一个肿瘤和一个匹配正常实体")

tumor_entities <- entity_map[
  analysis_role %chin% c("primary_tumor_measurement", "tumor_associated_input"),
  .(
    layer_id,
    file_id,
    mapped_patient_id = case_submitter_id,
    sample_submitter_id,
    sample_type,
    portion_submitter_id,
    analyte_submitter_id,
    aliquot_id,
    aliquot_submitter_id,
    analysis_role
  )
]
normal_entities <- entity_map[
  analysis_role == "matched_normal_associated_input",
  .(
    layer_id,
    file_id,
    matched_normal_submitter_id = aliquot_submitter_id,
    matched_normal_sample_type = sample_type
  )
]

file_selection <- manifest[, .(
  layer_id,
  file_id,
  patient_id,
  file_name,
  target_relative_path,
  file_size = as.numeric(file_size),
  local_sha256,
  workflow_type,
  workflow_version
)]
file_selection <- merge(
  file_selection,
  tumor_entities,
  by = c("layer_id", "file_id"),
  all.x = TRUE,
  sort = FALSE
)
file_selection <- merge(
  file_selection,
  normal_entities,
  by = c("layer_id", "file_id"),
  all.x = TRUE,
  sort = FALSE
)
fail_if(anyNA(file_selection$mapped_patient_id), "文件缺少肿瘤实体映射")
fail_if(any(file_selection$patient_id != file_selection$mapped_patient_id),
        "file manifest 与 entity map 的患者映射不一致")
file_selection[, `:=`(
  selected_for_patient_matrix = TRUE,
  selection_rule = "该患者该层唯一文件",
  feature_record_count = NA_real_,
  total_read_count = NA_real_,
  detected_feature_count = NA_real_,
  missing_value_count = NA_real_,
  missing_fraction = NA_real_
)]

for (layer in names(expected_file_counts)) {
  add_qc(
    "layer",
    layer,
    "all",
    "source_file_count",
    actual_file_counts[[layer]],
    expected_file_counts[[layer]],
    "hard",
    "pass",
    "继续",
    "固定 DR45 文件数"
  )
}
add_qc(
  "cohort", "all", "all", "escc_patient_whitelist",
  nrow(whitelist), 96L, "hard", "pass", "继续",
  "仅纳入 DR45 病理确认 ESCC"
)

# 补充患者级人口学信息；原始 JSON 只作同版本临床协变量来源。
clinical_hits <- fromJSON(clinical_json_path, simplifyVector = FALSE)$data$hits
demographics <- rbindlist(lapply(clinical_hits, function(hit) {
  d <- hit$demographic
  data.table(
    patient_id = hit$submitter_id %||% NA_character_,
    age_at_index = d$age_at_index %||% NA_real_,
    sex_at_birth = d$sex_at_birth %||% NA_character_,
    race = d$race %||% NA_character_,
    ethnicity = d$ethnicity %||% NA_character_,
    country_of_residence_at_enrollment =
      d$country_of_residence_at_enrollment %||% NA_character_
  )
}), fill = TRUE)
fail_if(anyDuplicated(demographics$patient_id) > 0L, "临床人口学患者 ID 重复")

setnames(whitelist, "patient", "patient_id")
patient_coldata <- merge(
  whitelist,
  demographics,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)
patient_coldata <- merge(
  patient_coldata,
  availability,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("", "_availability")
)
setorder(patient_coldata, patient_id)
fail_if(nrow(patient_coldata) != 96L, "患者级 colData 构建后不是 96 行")

stage_dir <- tempfile(pattern = ".tcga_multiassay_", tmpdir = work_intermediate_dir)
dir_create(stage_dir)
on.exit({
  if (dir_exists(stage_dir)) dir_delete(stage_dir)
}, add = TRUE)

message("[1/7] 构建 RNA STAR counts/TPM/TMM-logCPM")
rna_files <- file_selection[layer_id == "rna_star"]
setorder(rna_files, patient_id)
read_rna <- function(path) {
  fread(
    path,
    skip = "gene_id",
    select = c(
      "gene_id", "gene_name", "gene_type", "unstranded", "tpm_unstranded"
    ),
    na.strings = c("", "NA"),
    showProgress = FALSE
  )
}
rna_first <- read_rna(file.path(canonical_root, rna_files$target_relative_path[[1]]))
rna_first <- rna_first[!startsWith(gene_id, "N_")]
rna_ids <- as.character(rna_first$gene_id)
fail_if(length(rna_ids) != 60660L, "RNA 去除四个 STAR 汇总行后应为 60,660 个基因")
fail_if(anyDuplicated(rna_ids) > 0L, "RNA gene_id 重复")
rna_counts <- matrix(
  0L,
  nrow = length(rna_ids),
  ncol = nrow(rna_files),
  dimnames = list(rna_ids, rna_files$patient_id)
)
rna_tpm <- matrix(
  NA_real_,
  nrow = length(rna_ids),
  ncol = nrow(rna_files),
  dimnames = list(rna_ids, rna_files$patient_id)
)
rna_library <- numeric(nrow(rna_files))
rna_detected <- integer(nrow(rna_files))
for (j in seq_len(nrow(rna_files))) {
  dt <- read_rna(file.path(canonical_root, rna_files$target_relative_path[[j]]))
  dt <- dt[!startsWith(gene_id, "N_")]
  counts_j <- align_values(rna_ids, dt$gene_id, dt$unstranded, "RNA")
  tpm_j <- align_values(rna_ids, dt$gene_id, dt$tpm_unstranded, "RNA")
  fail_if(anyNA(counts_j) || any(counts_j < 0) || any(counts_j != floor(counts_j)),
          "RNA unstranded counts 必须为非负整数")
  rna_counts[, j] <- as.integer(counts_j)
  rna_tpm[, j] <- as.numeric(tpm_j)
  rna_library[[j]] <- sum(counts_j)
  rna_detected[[j]] <- sum(counts_j > 0)
}
rna_dge <- DGEList(counts = rna_counts)
rna_dge <- calcNormFactors(rna_dge, method = "TMM")
rna_logcpm <- cpm(rna_dge, log = TRUE, prior.count = 2)
rna_annotation <- rna_first[, .(gene_id, gene_name, gene_type)]
rna_coldata <- make_measurement_coldata(file_selection, "rna_star", rna_files$patient_id)
rna_se <- SummarizedExperiment(
  assays = list(counts = rna_counts, tpm = rna_tpm, tmm_logcpm = rna_logcpm),
  rowData = S4Vectors::DataFrame(as.data.frame(rna_annotation), check.names = FALSE),
  colData = rna_coldata
)
rna_metrics <- data.table(
  file_id = rna_files$file_id,
  feature_record_count = nrow(rna_counts),
  total_read_count = rna_library,
  detected_feature_count = rna_detected,
  missing_value_count = 0,
  missing_fraction = 0
)
file_selection[rna_metrics, on = "file_id", `:=`(
  feature_record_count = i.feature_record_count,
  total_read_count = i.total_read_count,
  detected_feature_count = i.detected_feature_count,
  missing_value_count = i.missing_value_count,
  missing_fraction = i.missing_fraction
)]
add_soft_metric("RNA", rna_files$patient_id, "library_size", rna_library,
                "TMM 归一化前总 unstranded counts")
add_soft_metric("RNA", rna_files$patient_id, "detected_gene_count", rna_detected,
                "count > 0 的基因数")
add_qc("layer", "RNA", "all", "matrix_dimension",
       paste(dim(rna_counts), collapse = "x"), "60660x95", "hard", "pass", "继续")

message("[2/7] 构建 miRNA counts/RPM，并冻结重复 aliquot")
mirna_files <- file_selection[layer_id == "mirna_bcgsc"]
setorder(mirna_files, patient_id, file_id)
read_mirna <- function(path) {
  dt <- fread(path, na.strings = c("", "NA"), showProgress = FALSE)
  setnames(dt, "cross-mapped", "cross_mapped")
  dt
}
mirna_first <- read_mirna(file.path(canonical_root, mirna_files$target_relative_path[[1]]))
mirna_ids <- as.character(mirna_first$miRNA_ID)
fail_if(length(mirna_ids) != 1881L || anyDuplicated(mirna_ids) > 0L,
        "miRNA feature 集必须为 1,881 个唯一 miRBase 21 ID")
mirna_counts_all <- matrix(
  0L,
  nrow = length(mirna_ids),
  ncol = nrow(mirna_files),
  dimnames = list(mirna_ids, mirna_files$file_id)
)
mirna_rpm_all <- matrix(
  NA_real_,
  nrow = length(mirna_ids),
  ncol = nrow(mirna_files),
  dimnames = list(mirna_ids, mirna_files$file_id)
)
mirna_crossmap_n <- integer(length(mirna_ids))
for (j in seq_len(nrow(mirna_files))) {
  dt <- read_mirna(file.path(canonical_root, mirna_files$target_relative_path[[j]]))
  counts_j <- align_values(mirna_ids, dt$miRNA_ID, dt$read_count, "miRNA")
  rpm_j <- align_values(
    mirna_ids,
    dt$miRNA_ID,
    dt$reads_per_million_miRNA_mapped,
    "miRNA"
  )
  cross_j <- align_values(mirna_ids, dt$miRNA_ID, dt$cross_mapped, "miRNA")
  fail_if(anyNA(counts_j) || any(counts_j < 0) || any(counts_j != floor(counts_j)),
          "miRNA read_count 必须为非负整数")
  mirna_counts_all[, j] <- as.integer(counts_j)
  mirna_rpm_all[, j] <- as.numeric(rpm_j)
  mirna_crossmap_n <- mirna_crossmap_n + as.integer(cross_j == "Y")
}
mirna_metrics <- data.table(
  file_id = mirna_files$file_id,
  patient_id = mirna_files$patient_id,
  detected_feature_count = colSums(mirna_counts_all > 0),
  total_read_count = colSums(mirna_counts_all)
)
setorder(
  mirna_metrics,
  patient_id,
  -detected_feature_count,
  -total_read_count,
  file_id
)
mirna_metrics[, selected := seq_len(.N) == 1L, by = patient_id]
selected_mirna_ids <- mirna_metrics[selected == TRUE, file_id]
file_selection[layer_id == "mirna_bcgsc", `:=`(
  selected_for_patient_matrix = file_id %chin% selected_mirna_ids,
  selection_rule = fifelse(
    file_id %chin% selected_mirna_ids,
    "同患者优先检测 miRNA 数较多者，其次总 reads 较高者；唯一文件直接入选",
    "同患者重复 aliquot 未入主矩阵，保留作敏感性分析"
  )
)]
mirna_metrics[, `:=`(
  feature_record_count = length(mirna_ids),
  missing_value_count = 0,
  missing_fraction = 0
)]
file_selection[mirna_metrics, on = "file_id", `:=`(
  feature_record_count = i.feature_record_count,
  total_read_count = i.total_read_count,
  detected_feature_count = i.detected_feature_count,
  missing_value_count = i.missing_value_count,
  missing_fraction = i.missing_fraction
)]
mirna_selected_files <- file_selection[
  layer_id == "mirna_bcgsc" & selected_for_patient_matrix
]
setorder(mirna_selected_files, patient_id)
mirna_selected_idx <- match(mirna_selected_files$file_id, colnames(mirna_counts_all))
mirna_counts <- mirna_counts_all[, mirna_selected_idx, drop = FALSE]
mirna_rpm <- mirna_rpm_all[, mirna_selected_idx, drop = FALSE]
colnames(mirna_counts) <- mirna_selected_files$patient_id
colnames(mirna_rpm) <- mirna_selected_files$patient_id
mirna_dge <- DGEList(counts = mirna_counts)
mirna_dge <- calcNormFactors(mirna_dge, method = "TMM")
mirna_logcpm <- cpm(mirna_dge, log = TRUE, prior.count = 2)
mirna_annotation <- data.table(
  miRNA_ID = mirna_ids,
  cross_mapped_file_count = mirna_crossmap_n
)
mirna_coldata <- make_measurement_coldata(
  file_selection,
  "mirna_bcgsc",
  mirna_selected_files$patient_id
)
mirna_se <- SummarizedExperiment(
  assays = list(counts = mirna_counts, rpm = mirna_rpm, tmm_logcpm = mirna_logcpm),
  rowData = S4Vectors::DataFrame(as.data.frame(mirna_annotation), check.names = FALSE),
  colData = mirna_coldata
)
duplicate_patient <- mirna_metrics[, .N, by = patient_id][N > 1L, patient_id]
fail_if(!identical(duplicate_patient, "TCGA-IG-A3YB"),
        "miRNA 重复患者不符合冻结预期")
duplicate_file_ids <- mirna_metrics[patient_id == duplicate_patient, file_id]
duplicate_cor <- cor(
  log2(mirna_rpm_all[, duplicate_file_ids[[1]]] + 1),
  log2(mirna_rpm_all[, duplicate_file_ids[[2]]] + 1),
  method = "spearman",
  use = "pairwise.complete.obs"
)
add_qc(
  "sample", "miRNA", duplicate_patient, "duplicate_aliquot_log2rpm_spearman",
  sprintf("%.6f", duplicate_cor), "仅作敏感性；不平均 raw count", "soft",
  "soft_flag", "主 aliquot 已冻结；另一 aliquot 保留作替换敏感性",
  paste("主文件", mirna_selected_files[patient_id == duplicate_patient, file_id])
)
add_soft_metric("miRNA", mirna_selected_files$patient_id, "library_size",
                colSums(mirna_counts), "主 aliquot 的总 read_count")
add_soft_metric("miRNA", mirna_selected_files$patient_id, "detected_mirna_count",
                colSums(mirna_counts > 0), "主 aliquot 中 read_count > 0")
add_qc("layer", "miRNA", "all", "matrix_dimension",
       paste(dim(mirna_counts), collapse = "x"), "1881x95", "hard", "pass", "继续")

rm(mirna_counts_all, mirna_rpm_all)
gc(verbose = FALSE)

message("[3/7] 构建 HM450 全 probe beta 矩阵")
meth_files <- file_selection[layer_id == "methylation_hm450"]
setorder(meth_files, patient_id)
read_methylation <- function(path) {
  fread(
    path,
    header = FALSE,
    col.names = c("probe_id", "beta"),
    na.strings = c("", "NA", "NaN"),
    showProgress = FALSE
  )
}
meth_first <- read_methylation(file.path(canonical_root, meth_files$target_relative_path[[1]]))
probe_ids <- as.character(meth_first$probe_id)
fail_if(length(probe_ids) != 486427L || anyDuplicated(probe_ids) > 0L,
        "HM450 必须为 486,427 个唯一 probe")
meth_beta <- matrix(
  NA_real_,
  nrow = length(probe_ids),
  ncol = nrow(meth_files),
  dimnames = list(probe_ids, meth_files$patient_id)
)
meth_nonmissing <- integer(nrow(meth_files))
meth_missing <- integer(nrow(meth_files))
for (j in seq_len(nrow(meth_files))) {
  dt <- read_methylation(file.path(canonical_root, meth_files$target_relative_path[[j]]))
  beta_j <- align_values(probe_ids, dt$probe_id, dt$beta, "HM450")
  fail_if(any(beta_j[!is.na(beta_j)] < 0 | beta_j[!is.na(beta_j)] > 1),
          "HM450 beta 必须在 [0,1] 或 NA")
  meth_beta[, j] <- as.numeric(beta_j)
  meth_nonmissing[[j]] <- sum(!is.na(beta_j))
  meth_missing[[j]] <- sum(is.na(beta_j))
}
probe_class <- fifelse(
  grepl("^cg[0-9]+$", probe_ids), "CpG",
  fifelse(grepl("^ch[0-9]+$", probe_ids), "non_CpG",
          fifelse(grepl("^rs[0-9]+$", probe_ids), "SNP", "other"))
)
meth_annotation <- data.table(probe_id = probe_ids, probe_class = probe_class)
meth_coldata <- make_measurement_coldata(
  file_selection,
  "methylation_hm450",
  meth_files$patient_id
)
meth_se <- SummarizedExperiment(
  assays = list(beta = meth_beta),
  rowData = S4Vectors::DataFrame(as.data.frame(meth_annotation), check.names = FALSE),
  colData = meth_coldata
)
meth_metrics <- data.table(
  file_id = meth_files$file_id,
  feature_record_count = nrow(meth_beta),
  detected_feature_count = meth_nonmissing,
  missing_value_count = meth_missing,
  missing_fraction = meth_missing / nrow(meth_beta)
)
file_selection[meth_metrics, on = "file_id", `:=`(
  feature_record_count = i.feature_record_count,
  detected_feature_count = i.detected_feature_count,
  missing_value_count = i.missing_value_count,
  missing_fraction = i.missing_fraction
)]
add_soft_metric("HM450", meth_files$patient_id, "missing_fraction",
                meth_missing / nrow(meth_beta),
                "NA 不补 0；后续按 probe 覆盖率过滤")
add_qc("layer", "HM450", "all", "matrix_dimension",
       paste(dim(meth_beta), collapse = "x"), "486427x96", "hard", "pass", "继续")

message("[4/7] 构建 ASCAT2 gene CNV 三个 assay 与 segment 长表")
cnv_files <- file_selection[layer_id == "cnv_gene_ascat2"]
setorder(cnv_files, patient_id)
read_cnv_gene <- function(path) {
  fread(path, na.strings = c("", "NA", "NaN"), showProgress = FALSE)
}
cnv_first <- read_cnv_gene(file.path(canonical_root, cnv_files$target_relative_path[[1]]))
cnv_ids <- as.character(cnv_first$gene_id)
fail_if(length(cnv_ids) != 60623L || anyDuplicated(cnv_ids) > 0L,
        "ASCAT2 gene CNV 必须为 60,623 个唯一基因")
cnv_total <- matrix(
  NA_real_, nrow = length(cnv_ids), ncol = nrow(cnv_files),
  dimnames = list(cnv_ids, cnv_files$patient_id)
)
cnv_min <- cnv_total
cnv_max <- cnv_total
cnv_missing <- integer(nrow(cnv_files))
for (j in seq_len(nrow(cnv_files))) {
  dt <- read_cnv_gene(file.path(canonical_root, cnv_files$target_relative_path[[j]]))
  total_j <- align_values(cnv_ids, dt$gene_id, dt$copy_number, "ASCAT2 gene CNV")
  min_j <- align_values(cnv_ids, dt$gene_id, dt$min_copy_number, "ASCAT2 gene CNV")
  max_j <- align_values(cnv_ids, dt$gene_id, dt$max_copy_number, "ASCAT2 gene CNV")
  fail_if(any(total_j[!is.na(total_j)] < 0) ||
            any(min_j[!is.na(min_j)] < 0) ||
            any(max_j[!is.na(max_j)] < 0),
          "ASCAT2 copy number 不能为负数")
  valid_triplet <- !is.na(total_j) & !is.na(min_j) & !is.na(max_j)
  fail_if(any(min_j[valid_triplet] > total_j[valid_triplet]) ||
            any(total_j[valid_triplet] > max_j[valid_triplet]),
          "ASCAT2 min/total/max 关系异常")
  cnv_total[, j] <- as.numeric(total_j)
  cnv_min[, j] <- as.numeric(min_j)
  cnv_max[, j] <- as.numeric(max_j)
  cnv_missing[[j]] <- sum(is.na(total_j))
}
cnv_annotation <- cnv_first[, .(gene_id, gene_name, chromosome, start, end)]
cnv_coldata <- make_measurement_coldata(
  file_selection,
  "cnv_gene_ascat2",
  cnv_files$patient_id
)
cnv_se <- SummarizedExperiment(
  assays = list(
    copy_number = cnv_total,
    min_copy_number = cnv_min,
    max_copy_number = cnv_max
  ),
  rowData = S4Vectors::DataFrame(as.data.frame(cnv_annotation), check.names = FALSE),
  colData = cnv_coldata
)
cnv_metrics <- data.table(
  file_id = cnv_files$file_id,
  feature_record_count = nrow(cnv_total),
  detected_feature_count = nrow(cnv_total) - cnv_missing,
  missing_value_count = cnv_missing,
  missing_fraction = cnv_missing / nrow(cnv_total)
)
file_selection[cnv_metrics, on = "file_id", `:=`(
  feature_record_count = i.feature_record_count,
  detected_feature_count = i.detected_feature_count,
  missing_value_count = i.missing_value_count,
  missing_fraction = i.missing_fraction
)]
add_soft_metric("CNV_gene", cnv_files$patient_id, "missing_fraction",
                cnv_missing / nrow(cnv_total),
                "空 copy number 保留 NA，不补二倍体或 0")
add_qc("layer", "CNV_gene", "all", "matrix_dimension",
       paste(dim(cnv_total), collapse = "x"), "60623x96", "hard", "pass", "继续")

segment_files <- file_selection[layer_id == "cnv_segment_ascat2"]
setorder(segment_files, patient_id)
segment_list <- vector("list", nrow(segment_files))
segment_counts <- integer(nrow(segment_files))
for (j in seq_len(nrow(segment_files))) {
  dt <- fread(
    file.path(canonical_root, segment_files$target_relative_path[[j]]),
    na.strings = c("", "NA", "NaN"),
    showProgress = FALSE
  )
  tumor_aliquot_id <- segment_files$aliquot_id[[j]]
  fail_if(uniqueN(dt$GDC_Aliquot) != 1L || unique(dt$GDC_Aliquot) != tumor_aliquot_id,
          "ASCAT2 segment 的 GDC_Aliquot 与肿瘤实体不一致")
  fail_if(any(dt$Start > dt$End), "ASCAT2 segment 存在 Start > End")
  fail_if(any(dt$Copy_Number < 0 | dt$Major_Copy_Number < 0 |
                dt$Minor_Copy_Number < 0, na.rm = TRUE),
          "ASCAT2 segment copy number 不能为负数")
  dt[, `:=`(
    patient_id = segment_files$patient_id[[j]],
    file_id = segment_files$file_id[[j]],
    tumor_aliquot_id = tumor_aliquot_id,
    segment_length = End - Start + 1
  )]
  segment_list[[j]] <- dt
  segment_counts[[j]] <- nrow(dt)
}
segments_long <- rbindlist(segment_list, use.names = TRUE, fill = TRUE)
segment_metrics <- data.table(
  file_id = segment_files$file_id,
  feature_record_count = segment_counts,
  detected_feature_count = segment_counts,
  missing_value_count = 0,
  missing_fraction = 0
)
file_selection[segment_metrics, on = "file_id", `:=`(
  feature_record_count = i.feature_record_count,
  detected_feature_count = i.detected_feature_count,
  missing_value_count = i.missing_value_count,
  missing_fraction = i.missing_fraction
)]
add_soft_metric("CNV_segment", segment_files$patient_id, "segment_count",
                segment_counts, "与 gene CNV 为同一 ASCAT2 证据的区间表示")

message("[5/7] 合并 masked MAF 并构建基因×患者稀疏突变矩阵")
maf_files <- file_selection[layer_id == "masked_maf"]
setorder(maf_files, patient_id)
maf_columns <- c(
  "Hugo_Symbol", "Entrez_Gene_Id", "NCBI_Build", "Chromosome",
  "Start_Position", "End_Position", "Variant_Classification", "Variant_Type",
  "Reference_Allele", "Tumor_Seq_Allele1", "Tumor_Seq_Allele2",
  "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
  "Tumor_Sample_UUID", "Matched_Norm_Sample_UUID", "HGVSc", "HGVSp",
  "HGVSp_Short", "Transcript_ID", "t_depth", "t_ref_count", "t_alt_count",
  "n_depth", "n_ref_count", "n_alt_count", "Gene", "Feature", "Consequence",
  "IMPACT", "PICK", "GDC_FILTER", "hotspot", "callers", "Mutation_Status"
)
maf_list <- vector("list", nrow(maf_files))
maf_counts_per_file <- integer(nrow(maf_files))
for (j in seq_len(nrow(maf_files))) {
  dt <- fread(
    file.path(canonical_root, maf_files$target_relative_path[[j]]),
    skip = "Hugo_Symbol",
    select = maf_columns,
    na.strings = c("", "NA", "NaN"),
    showProgress = FALSE
  )
  if (nrow(dt)) {
    fail_if(any(dt$NCBI_Build != "GRCh38", na.rm = TRUE), "MAF 必须为 GRCh38")
    fail_if(any(dt$Mutation_Status != "Somatic", na.rm = TRUE), "MAF 必须为 Somatic")
    patient_from_barcode <- substr(dt$Tumor_Sample_Barcode, 1L, 12L)
    fail_if(any(patient_from_barcode != maf_files$patient_id[[j]]),
            "MAF 肿瘤条形码与 manifest 患者不一致")
  }
  dt[, `:=`(
    patient_id = maf_files$patient_id[[j]],
    file_id = maf_files$file_id[[j]],
    tumor_aliquot_submitter_id = maf_files$aliquot_submitter_id[[j]],
    matched_normal_submitter_id = maf_files$matched_normal_submitter_id[[j]]
  )]
  maf_list[[j]] <- dt
  maf_counts_per_file[[j]] <- nrow(dt)
}
maf_long <- rbindlist(maf_list, use.names = TRUE, fill = TRUE)
maf_variant_key <- c(
  "patient_id", "Chromosome", "Start_Position", "End_Position",
  "Reference_Allele", "Tumor_Seq_Allele2", "Gene"
)
maf_duplicates <- nrow(maf_long) - uniqueN(maf_long, by = maf_variant_key)
if (maf_duplicates > 0L) {
  maf_long <- unique(maf_long, by = maf_variant_key)
}
non_synonymous_classes <- c(
  "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Splice_Region",
  "Translation_Start_Site", "Nonsense_Mutation", "Nonstop_Mutation",
  "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation"
)
maf_long[, nonsynonymous := Variant_Classification %chin% non_synonymous_classes]
cnv_gene_key <- strip_ensembl_version(cnv_ids)
fail_if(anyDuplicated(cnv_gene_key) > 0L,
        "ASCAT2 去版本 Ensembl ID 存在重复，不能构建突变共同宇宙")
gene_lookup <- data.table(gene_key = cnv_gene_key, gene_id = cnv_ids)
maf_long[, gene_key := strip_ensembl_version(Gene)]
maf_long <- merge(maf_long, gene_lookup, by = "gene_key", all.x = TRUE, sort = FALSE)
maf_nonsyn <- maf_long[nonsynonymous & !is.na(gene_id)]
mutation_counts_long <- maf_nonsyn[, .N, by = .(gene_id, patient_id)]
mutation_count <- sparseMatrix(
  i = match(mutation_counts_long$gene_id, cnv_ids),
  j = match(mutation_counts_long$patient_id, patient_coldata$patient_id),
  x = mutation_counts_long$N,
  dims = c(length(cnv_ids), nrow(patient_coldata)),
  dimnames = list(cnv_ids, patient_coldata$patient_id)
)
mutation_binary <- mutation_count
mutation_binary@x[] <- 1
mutation_coldata <- make_measurement_coldata(
  file_selection,
  "masked_maf",
  patient_coldata$patient_id
)
mutation_se <- SummarizedExperiment(
  assays = list(count = mutation_count, binary = mutation_binary),
  rowData = S4Vectors::DataFrame(as.data.frame(cnv_annotation), check.names = FALSE),
  colData = mutation_coldata
)
maf_patient_nonsyn <- data.table(patient_id = patient_coldata$patient_id)
maf_patient_nonsyn <- maf_nonsyn[, .N, by = patient_id][maf_patient_nonsyn, on = "patient_id"]
maf_patient_nonsyn[is.na(N), N := 0L]
maf_metrics <- data.table(
  file_id = maf_files$file_id,
  feature_record_count = maf_counts_per_file,
  detected_feature_count = maf_counts_per_file,
  missing_value_count = 0,
  missing_fraction = 0
)
file_selection[maf_metrics, on = "file_id", `:=`(
  feature_record_count = i.feature_record_count,
  detected_feature_count = i.detected_feature_count,
  missing_value_count = i.missing_value_count,
  missing_fraction = i.missing_fraction
)]
for (i in seq_len(nrow(maf_patient_nonsyn))) {
  add_qc(
    "sample", "Mutation", maf_patient_nonsyn$patient_id[[i]],
    "nonsynonymous_variant_count", maf_patient_nonsyn$N[[i]],
    "零事件为合法结果；仅作队列分布复核", "soft", "pass", "保留",
    "不按零事件自动失败"
  )
}
add_qc(
  "layer", "Mutation", "all", "unmapped_nonsynonymous_variants",
  maf_long[nonsynonymous & is.na(gene_id), .N],
  "保留在长表，不强行映射到共同基因宇宙", "soft", "pass", "记录边界"
)
add_qc(
  "layer", "Mutation", "all", "exact_duplicate_variant_rows_removed",
  maf_duplicates, "0 为理想；若存在则按患者-坐标-等位基因-基因去重",
  "soft", if (maf_duplicates > 0L) "soft_flag" else "pass",
  "长表保留唯一变异事件"
)

message("[6/7] 构建 GDC RPPA 来源矩阵")
rppa_files <- file_selection[layer_id == "rppa"]
setorder(rppa_files, patient_id)
read_rppa <- function(path) {
  fread(path, na.strings = c("", "NA", "NaN"), showProgress = FALSE)
}
rppa_first <- read_rppa(file.path(canonical_root, rppa_files$target_relative_path[[1]]))
rppa_feature_key <- paste(
  rppa_first$AGID,
  rppa_first$peptide_target,
  rppa_first$set_id,
  sep = "|"
)
fail_if(length(rppa_feature_key) != 487L || anyDuplicated(rppa_feature_key) > 0L,
        "RPPA 必须为 487 个唯一 AGID|peptide_target|set_id 特征")
rppa_matrix <- matrix(
  NA_real_,
  nrow = length(rppa_feature_key),
  ncol = nrow(rppa_files),
  dimnames = list(rppa_feature_key, rppa_files$patient_id)
)
rppa_missing <- integer(nrow(rppa_files))
for (j in seq_len(nrow(rppa_files))) {
  dt <- read_rppa(file.path(canonical_root, rppa_files$target_relative_path[[j]]))
  key_j <- paste(dt$AGID, dt$peptide_target, dt$set_id, sep = "|")
  expression_j <- align_values(
    rppa_feature_key,
    key_j,
    dt$protein_expression,
    "RPPA"
  )
  rppa_matrix[, j] <- as.numeric(expression_j)
  rppa_missing[[j]] <- sum(is.na(expression_j))
}
rppa_annotation <- rppa_first[, .(
  feature_id = rppa_feature_key,
  AGID,
  lab_id,
  catalog_number,
  set_id,
  peptide_target
)]
rppa_coldata <- make_measurement_coldata(
  file_selection,
  "rppa",
  rppa_files$patient_id
)
rppa_se <- SummarizedExperiment(
  assays = list(protein_expression = rppa_matrix),
  rowData = S4Vectors::DataFrame(as.data.frame(rppa_annotation), check.names = FALSE),
  colData = rppa_coldata
)
rppa_metrics <- data.table(
  file_id = rppa_files$file_id,
  feature_record_count = nrow(rppa_matrix),
  detected_feature_count = nrow(rppa_matrix) - rppa_missing,
  missing_value_count = rppa_missing,
  missing_fraction = rppa_missing / nrow(rppa_matrix)
)
file_selection[rppa_metrics, on = "file_id", `:=`(
  feature_record_count = i.feature_record_count,
  detected_feature_count = i.detected_feature_count,
  missing_value_count = i.missing_value_count,
  missing_fraction = i.missing_fraction
)]
add_soft_metric("RPPA", rppa_files$patient_id, "missing_fraction",
                rppa_missing / nrow(rppa_matrix),
                "GDC RPPA 为来源层；正式蛋白解释优先 TCPA Level 4 校准")
add_qc("layer", "RPPA", "all", "matrix_dimension",
       paste(dim(rppa_matrix), collapse = "x"), "487x78", "hard", "pass", "继续")

message("[7/7] 组装 MultiAssayExperiment、分析集与正式输出")
patient_df <- as.data.frame(patient_coldata)
rownames(patient_df) <- patient_df$patient_id
mae_coldata <- S4Vectors::DataFrame(patient_df, check.names = FALSE)
experiments <- list(
  RNA = rna_se,
  miRNA = mirna_se,
  HM450 = meth_se,
  Mutation = mutation_se,
  CNV_gene = cnv_se,
  RPPA = rppa_se
)
mae_sample_map <- rbindlist(lapply(names(experiments), function(assay_name) {
  ids <- colnames(experiments[[assay_name]])
  data.table(assay = assay_name, primary = ids, colname = ids)
}))
mae <- MultiAssayExperiment(
  experiments = ExperimentList(experiments),
  colData = mae_coldata,
  sampleMap = S4Vectors::DataFrame(as.data.frame(mae_sample_map))
)
validObject(mae)

actual_sets <- data.table(
  analysis_set = c(
    "clinical_whitelist", "driver_core", "five_layer_core", "protein_deep_subset"
  ),
  patient_count = c(
    nrow(patient_coldata),
    sum(as_flag(patient_coldata$driver_core)),
    sum(as_flag(patient_coldata$five_layer_core)),
    sum(as_flag(patient_coldata$protein_deep_subset))
  ),
  expected_count = c(96L, 95L, 94L, 76L)
)
fail_if(any(actual_sets$patient_count != actual_sets$expected_count),
        "MultiAssayExperiment 分析集不符合 96/95/94/76")
for (i in seq_len(nrow(actual_sets))) {
  add_qc(
    "analysis_set", "all", actual_sets$analysis_set[[i]], "patient_count",
    actual_sets$patient_count[[i]], actual_sets$expected_count[[i]],
    "hard", "pass", "继续",
    "嵌套分析集，不以 76 例替代全部单层样本"
  )
}
add_qc(
  "object", "all", "tcga_escc_multiassay_dr45", "validObject",
  "TRUE", "TRUE", "hard", "pass", "继续",
  "MultiAssayExperiment 结构有效"
)

required_by_set <- list(
  clinical_whitelist = character(),
  driver_core = c("rna_star_available", "masked_maf_available", "cnv_gene_ascat2_available"),
  five_layer_core = c(
    "rna_star_available", "mirna_bcgsc_available", "methylation_hm450_available",
    "masked_maf_available", "cnv_gene_ascat2_available"
  ),
  protein_deep_subset = c(
    "rna_star_available", "mirna_bcgsc_available", "methylation_hm450_available",
    "masked_maf_available", "cnv_gene_ascat2_available", "rppa_available"
  )
)
analysis_sets_long <- rbindlist(lapply(names(required_by_set), function(set_name) {
  required_cols <- required_by_set[[set_name]]
  if (!length(required_cols)) {
    included <- rep(TRUE, nrow(patient_coldata))
    missing_layers <- rep("", nrow(patient_coldata))
  } else {
    avail_matrix <- as.data.frame(patient_coldata[, ..required_cols])
    avail_matrix[] <- lapply(avail_matrix, as_flag)
    included <- apply(avail_matrix, 1L, all)
    missing_layers <- apply(avail_matrix, 1L, function(x) {
      paste(sub("_available$", "", required_cols[!x]), collapse = ";")
    })
  }
  data.table(
    analysis_set = set_name,
    patient_id = patient_coldata$patient_id,
    included = included,
    exclusion_reason = fifelse(included, "", paste0("missing:", missing_layers)),
    independent_inference_unit = "patient"
  )
}))

setorder(file_selection, layer_id, patient_id, file_id)
file_selection[, mapped_patient_id := NULL]
sample_map_out <- file_selection[, .(
  layer_id,
  file_id,
  patient_id,
  file_name,
  selected_for_patient_matrix,
  selection_rule,
  sample_submitter_id,
  sample_type,
  portion_submitter_id,
  analyte_submitter_id,
  aliquot_id,
  aliquot_submitter_id,
  matched_normal_submitter_id,
  matched_normal_sample_type,
  analysis_role,
  feature_record_count,
  total_read_count,
  detected_feature_count,
  missing_value_count,
  missing_fraction,
  target_relative_path,
  local_sha256
)]

soft_flag_count <- qc_rows[gate_type == "soft" & status == "soft_flag", .N]
hard_fail_count <- qc_rows[gate_type == "hard" & status != "pass", .N]
fail_if(hard_fail_count > 0L, "存在未通过的硬门禁，不发布正式对象")

summary_lines <- c(
  "# TCGA-ESCC GDC DR45 多组学对象构建摘要",
  "",
  "## 构建结果",
  "",
  "- 已从统一仓库中 653 个已验证开放文件构建一个患者并集 MultiAssayExperiment。",
  "- 患者母集为 96 例病理确认 ESCC；driver core、five-layer core 和 protein-deep subset 分别为 95、94 和 76 例。",
  "- RNA：60,660 基因 × 95 例，保留整数 counts、TPM 和 TMM-logCPM。",
  "- miRNA：1,881 个 miRBase 21 特征 × 95 例；TCGA-IG-A3YB 的两个 aliquot 未平均，主矩阵选择检测特征更多且总 reads 更高的 A360-13。",
  sprintf("- 两个 miRNA aliquot 的 log2(RPM+1) Spearman 相关为 %.6f；未入选文件保留在样本映射中用于敏感性分析。", duplicate_cor),
  "- HM450：486,427 个 probe × 96 例；beta 的 NA 保留，未补 0。",
  "- masked MAF：长表保留全部已选字段；非同义突变按 60,623 个共同 Ensembl 基因构建稀疏矩阵。",
  "- ASCAT2 gene CNV：60,623 基因 × 96 例，分别保存 total/min/max copy number；空值保留 NA。",
  sprintf("- ASCAT2 segment：%s 个区间，单独保存长表，不作为独立组学层重复计分。", format(nrow(segments_long), big.mark = ",")),
  "- GDC RPPA：487 个抗体特征 × 78 例；正式蛋白解释仍需 TCPA Level 4 校准。",
  "",
  "## 门禁边界",
  "",
  "- 样本身份、病理、患者去重、文件完整性、实体层级和矩阵结构是硬门禁。",
  sprintf("- 本次硬门禁失败数为 %d。", hard_fail_count),
  sprintf("- 软异常标记数为 %d；软标记不自动删除患者，进入层内 QC、PCA、替代 aliquot 或敏感性分析。", soft_flag_count),
  "- MAF 零事件、跨平台方向差异或单层缺失不自动判为失败；必须保留并在后续候选门禁表中分类。",
  "- 不把 76 例完整个案设为全项目硬门禁，各组学先用最大有效样本，再使用嵌套交集整合。",
  "",
  "## 证据上限",
  "",
  "该对象用于公共数据假设生成、患者级校准和跨层候选轴筛选，不构成因果机制或治疗靶点证明。"
)

stage_outputs <- file.path(stage_dir, c(
  "tcga_escc_multiassay_dr45.rds",
  "tcga_escc_masked_maf_long.rds",
  "tcga_escc_ascat2_segments_long.rds",
  "tcga_escc_multiassay_sample_map.tsv",
  "tcga_escc_multiassay_analysis_sets.tsv",
  "tcga_escc_multiassay_qc.tsv",
  "tcga_escc_multiassay_build_summary.md"
))
names(stage_outputs) <- path_file(stage_outputs)

saveRDS(mae, stage_outputs[["tcga_escc_multiassay_dr45.rds"]], compress = "gzip")
saveRDS(maf_long, stage_outputs[["tcga_escc_masked_maf_long.rds"]], compress = "gzip")
saveRDS(segments_long, stage_outputs[["tcga_escc_ascat2_segments_long.rds"]], compress = "gzip")
fwrite(sample_map_out, stage_outputs[["tcga_escc_multiassay_sample_map.tsv"]], sep = "\t", na = "")
fwrite(analysis_sets_long, stage_outputs[["tcga_escc_multiassay_analysis_sets.tsv"]], sep = "\t", na = "")
fwrite(qc_rows, stage_outputs[["tcga_escc_multiassay_qc.tsv"]], sep = "\t", na = "")
writeLines(summary_lines, stage_outputs[["tcga_escc_multiassay_build_summary.md"]], useBytes = TRUE)

# 重新读取正式对象候选，避免发布不可反序列化的 RDS。
rm(
  experiments, rna_se, mirna_se, meth_se, mutation_se, cnv_se, rppa_se,
  rna_counts, rna_tpm, rna_logcpm, mirna_counts, mirna_rpm, mirna_logcpm,
  meth_beta, mutation_count, mutation_binary, cnv_total, cnv_min, cnv_max,
  rppa_matrix, mae
)
gc(verbose = FALSE)
mae_check <- readRDS(stage_outputs[["tcga_escc_multiassay_dr45.rds"]])
validObject(mae_check)
fail_if(length(experiments(mae_check)) != 6L, "重新读取后 MAE assay 数异常")
rm(mae_check)
gc(verbose = FALSE)

artifact_manifest <- data.table(
  artifact = names(stage_outputs),
  relative_path = file.path("results", names(stage_outputs)),
  file_size_bytes = as.numeric(file_info(unname(stage_outputs))$size),
  sha256 = vapply(
    unname(stage_outputs),
    digest,
    character(1),
    algo = "sha256",
    file = TRUE
  ),
  generated_date = as.character(Sys.Date()),
  generation_script = "scripts/12_build_tcga_escc_multiassay.R",
  source_dataset_key = target_dataset_key,
  status = "verified"
)
artifact_manifest_path <- file.path(stage_dir, "tcga_escc_multiassay_artifact_manifest.tsv")
fwrite(artifact_manifest, artifact_manifest_path, sep = "\t", na = "")
stage_outputs <- c(stage_outputs, tcga_escc_multiassay_artifact_manifest.tsv = artifact_manifest_path)

for (name in names(stage_outputs)) {
  destination <- file.path(results_dir, name)
  file_copy(stage_outputs[[name]], destination, overwrite = TRUE)
  fail_if(
    digest(destination, algo = "sha256", file = TRUE) !=
      digest(stage_outputs[[name]], algo = "sha256", file = TRUE),
    paste0("发布后 SHA256 不一致：", name)
  )
}

session_lines <- c(
  "# TCGA-ESCC 多组学对象构建环境记录",
  "",
  "此文件是运行时历史证据，不是项目当前状态源。",
  "",
  paste0("- 构建日期：", Sys.Date()),
  paste0("- R：", R.version.string),
  paste0("- MultiAssayExperiment：", packageVersion("MultiAssayExperiment")),
  paste0("- SummarizedExperiment：", packageVersion("SummarizedExperiment")),
  paste0("- edgeR：", packageVersion("edgeR")),
  paste0("- data.table：", packageVersion("data.table")),
  paste0("- Matrix：", packageVersion("Matrix")),
  "- 源数据：GDC TCGA-ESCA DR45，ESCC 病理白名单 96 例。",
  "- 发布前重新读取 RDS 并通过 validObject()。"
)
writeLines(
  session_lines,
  file.path(work_checks_dir, "tcga_escc_multiassay_build_session_20260711.md"),
  useBytes = TRUE
)

if (dir_exists(stage_dir)) dir_delete(stage_dir)

message("完成：TCGA-ESCC MultiAssayExperiment 已发布到 results/")
message("分析集：96 / 95 / 94 / 76；硬门禁失败：", hard_fail_count,
        "；软标记：", soft_flag_count)
