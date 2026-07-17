#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

project_root <- normalizePath(getwd(), mustWork = TRUE)
registry_path <- file.path(project_root, "data", "datasets.tsv")
registry <- fread(registry_path, sep = "\t", colClasses = "character", na.strings = NULL)

require_one_path <- function(target_name) {
  value <- registry[registry$logical_name == target_name, central_path]
  if (length(value) != 1L || !dir.exists(value)) {
    stop("无法从 data/datasets.tsv 唯一定位：", target_name)
  }
  value
}

expression_root <- require_one_path("tcga_escc_bulk_expression")
clinical_root <- require_one_path("tcga_escc_clinical")

expression_path <- file.path(
  expression_root, "00_source", "TCGA-ESCA.star_counts.tsv.gz"
)
survival_path <- file.path(
  expression_root, "00_source", "TCGA-ESCA.survival.tsv.gz"
)
clinical_path <- file.path(
  clinical_root, "20_reusable", "tcga_esca_gdc_primary_diagnosis_flat_dr45.tsv"
)

required_files <- c(expression_path, survival_path, clinical_path)
if (!all(file.exists(required_files))) {
  stop("TCGA-ESCA 白名单输入缺失：", paste(required_files[!file.exists(required_files)], collapse = ";"))
}

clinical <- fread(clinical_path, colClasses = "character", na.strings = c("", "NA"))
rename_if_present <- function(old_name, new_name) {
  if (old_name %in% names(clinical) && !new_name %in% names(clinical)) {
    setnames(clinical, old_name, new_name)
  }
}
rename_if_present("case_submitter_id", "patient")
rename_if_present("vital_status", "gdc_vital_status")
rename_if_present("days_to_death", "gdc_days_to_death")
rename_if_present("days_to_last_follow_up", "gdc_days_to_last_follow_up")
if (!"tumor_stage" %in% names(clinical)) {
  clinical[, tumor_stage := NA_character_]
}
required_clinical <- c(
  "patient", "primary_diagnosis", "morphology", "is_escc",
  "gdc_vital_status", "gdc_days_to_death", "gdc_days_to_last_follow_up",
  "tumor_stage", "ajcc_pathologic_stage"
)
if (!all(required_clinical %in% names(clinical))) {
  stop("GDC 临床扁平表字段不完整。")
}
if (anyDuplicated(clinical$patient)) {
  stop("GDC 临床扁平表存在重复 patient。")
}

clinical[, is_escc_flag := tolower(as.character(is_escc)) == "true"]
clinical_escc <- clinical[is_escc_flag == TRUE]
if (nrow(clinical_escc) == 0L) {
  stop("没有识别到 ESCC 病例。")
}

expression_header <- names(fread(
  expression_path,
  sep = "\t",
  nrows = 0L,
  check.names = FALSE
))
if (length(expression_header) < 2L || expression_header[1] != "Ensembl_ID") {
  stop("STAR 表达矩阵表头不符合预期。")
}

sample_map <- data.table(sample = expression_header[-1L])
sample_map[, patient := substr(sample, 1L, 12L)]
sample_map[, sample_type_code := substr(sample, 14L, 15L)]
sample_map[, sample_type := fcase(
  sample_type_code == "01", "Primary Tumor",
  sample_type_code == "06", "Metastatic",
  sample_type_code == "11", "Solid Tissue Normal",
  default = "Other"
)]
sample_map <- sample_map[patient %in% clinical_escc$patient]
setorder(sample_map, patient, sample_type_code, sample)
sample_map[, patient_sample_type_count := .N, by = .(patient, sample_type_code)]
sample_map[, selected_primary_tumor :=
  sample_type_code == "01" & seq_len(.N) == 1L,
  by = patient
]

duplicate_primary <- sample_map[
  sample_type_code == "01",
  .N,
  by = patient
][N > 1L]
if (nrow(duplicate_primary) > 0L) {
  stop("同一 ESCC 患者存在多个处理后 primary tumor 样本，需人工冻结规则。")
}

survival <- fread(survival_path, sep = "\t", na.strings = c("", "NA"))
required_survival <- c("sample", "OS.time", "OS", "_PATIENT")
if (!all(required_survival %in% names(survival))) {
  stop("Xena 生存表字段不完整。")
}
survival_primary <- survival[
  `_PATIENT` %in% clinical_escc$patient & substr(sample, 14L, 15L) == "01"
]
if (anyDuplicated(survival_primary$`_PATIENT`)) {
  stop("Xena 生存表同一患者出现多个 primary tumor 生存记录。")
}
survival_patient <- survival_primary[, .(
  patient = `_PATIENT`,
  survival_source_sample = sample,
  xena_os_time_days = as.numeric(OS.time),
  xena_os_event = as.integer(OS)
)]

collapse_samples <- function(x) {
  x <- sort(unique(x))
  if (length(x) == 0L) "" else paste(x, collapse = ";")
}

patient_sample_summary <- sample_map[, .(
  primary_tumor_count = sum(sample_type_code == "01"),
  normal_count = sum(sample_type_code == "11"),
  metastatic_count = sum(sample_type_code == "06"),
  primary_tumor_samples = collapse_samples(sample[sample_type_code == "01"]),
  normal_samples = collapse_samples(sample[sample_type_code == "11"]),
  metastatic_samples = collapse_samples(sample[sample_type_code == "06"])
), by = patient]

patient_whitelist <- merge(
  clinical_escc,
  patient_sample_summary,
  by = "patient",
  all.x = TRUE,
  sort = FALSE
)
patient_whitelist <- merge(
  patient_whitelist,
  survival_patient,
  by = "patient",
  all.x = TRUE,
  sort = FALSE
)

count_columns <- c("primary_tumor_count", "normal_count", "metastatic_count")
for (column in count_columns) {
  set(patient_whitelist, which(is.na(patient_whitelist[[column]])), column, 0L)
}
sample_list_columns <- c("primary_tumor_samples", "normal_samples", "metastatic_samples")
for (column in sample_list_columns) {
  set(patient_whitelist, which(is.na(patient_whitelist[[column]])), column, "")
}

patient_whitelist[, main_primary_expression_eligible := primary_tumor_count == 1L]
patient_whitelist[, survival_model_eligible :=
  main_primary_expression_eligible & is.finite(xena_os_time_days) &
    xena_os_event %in% c(0L, 1L)
]
patient_whitelist[, whitelist_status := fcase(
  primary_tumor_count == 1L, "include_primary_expression",
  primary_tumor_count == 0L, "exclude_no_primary_expression",
  default = "review_multiple_primary_expression"
)]
patient_whitelist[, expression_unit := "log2(count+1)"]
patient_whitelist[, provider_same_sample_records_averaged := TRUE]

clinical_for_samples <- clinical_escc[, .(
  patient,
  primary_diagnosis,
  morphology,
  gdc_vital_status,
  tumor_stage,
  ajcc_pathologic_stage
)]
sample_map <- merge(sample_map, clinical_for_samples, by = "patient", all.x = TRUE, sort = FALSE)
sample_map <- merge(sample_map, survival_patient, by = "patient", all.x = TRUE, sort = FALSE)
sample_map[, analysis_role := fcase(
  selected_primary_tumor, "main_primary_expression",
  sample_type_code == "11", "context_normal_expression",
  sample_type_code == "06", "context_metastatic_expression",
  default = "exclude_other_sample_type"
)]
sample_map[, independent_inference_unit := "patient"]
sample_map[, expression_unit := "log2(count+1)"]
sample_map[, provider_same_sample_records_averaged := TRUE]

observed <- c(
  clinical_escc_patients = nrow(clinical_escc),
  expression_samples_escc = nrow(sample_map),
  primary_tumor_samples = sum(sample_map$sample_type_code == "01"),
  normal_samples = sum(sample_map$sample_type_code == "11"),
  metastatic_samples = sum(sample_map$sample_type_code == "06"),
  primary_expression_patients = sum(patient_whitelist$main_primary_expression_eligible),
  survival_eligible_patients = sum(patient_whitelist$survival_model_eligible),
  survival_events = sum(patient_whitelist$xena_os_event[patient_whitelist$survival_model_eligible] == 1L)
)
expected <- c(
  clinical_escc_patients = 96L,
  expression_samples_escc = 99L,
  primary_tumor_samples = 95L,
  normal_samples = 3L,
  metastatic_samples = 1L,
  primary_expression_patients = 95L,
  survival_eligible_patients = 94L,
  survival_events = 31L
)
summary_table <- data.table(
  metric = names(expected),
  expected = as.integer(expected),
  observed = as.integer(observed[names(expected)])
)
summary_table[, status := fifelse(expected == observed, "passed", "failed")]

patient_output <- file.path(project_root, "results", "tcga_escc_patient_whitelist.tsv")
sample_output <- file.path(project_root, "results", "tcga_escc_expression_sample_map.tsv")
summary_output <- file.path(project_root, "results", "tcga_escc_whitelist_summary.tsv")

patient_columns <- c(
  "patient", "primary_diagnosis", "morphology", "gdc_vital_status",
  "gdc_days_to_death", "gdc_days_to_last_follow_up", "tumor_stage",
  "ajcc_pathologic_stage", "primary_tumor_count", "normal_count",
  "metastatic_count", "primary_tumor_samples", "normal_samples",
  "metastatic_samples", "survival_source_sample", "xena_os_time_days",
  "xena_os_event", "main_primary_expression_eligible",
  "survival_model_eligible", "whitelist_status", "expression_unit",
  "provider_same_sample_records_averaged"
)
sample_columns <- c(
  "patient", "sample", "sample_type_code", "sample_type",
  "patient_sample_type_count", "selected_primary_tumor", "analysis_role",
  "primary_diagnosis", "morphology", "gdc_vital_status", "tumor_stage",
  "ajcc_pathologic_stage", "xena_os_time_days", "xena_os_event",
  "independent_inference_unit", "expression_unit",
  "provider_same_sample_records_averaged"
)

fwrite(patient_whitelist[, ..patient_columns], patient_output, sep = "\t", na = "NA", quote = FALSE)
fwrite(sample_map[, ..sample_columns], sample_output, sep = "\t", na = "NA", quote = FALSE)
fwrite(summary_table, summary_output, sep = "\t", na = "NA", quote = FALSE)

message("ESCC 临床患者：", observed[["clinical_escc_patients"]])
message("主分析原发肿瘤表达患者：", observed[["primary_expression_patients"]])
message("ESCC 正常表达样本：", observed[["normal_samples"]])
message("ESCC 转移表达样本：", observed[["metastatic_samples"]])
message("生存可用患者/事件：", observed[["survival_eligible_patients"]], "/", observed[["survival_events"]])

if (any(summary_table$status == "failed")) {
  stop("当前 TCGA 固定快照与预期样本计数不一致；请查看 tcga_escc_whitelist_summary.tsv。")
}
