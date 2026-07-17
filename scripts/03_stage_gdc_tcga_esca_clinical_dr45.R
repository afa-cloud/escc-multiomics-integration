#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(httr2)
  library(jsonlite)
})

project_root <- normalizePath(getwd(), mustWork = TRUE)
data_root <- Sys.getenv("RESEARCH_DATA_ROOT", unset = path.expand("~/ResearchDataHub"))
data_root <- normalizePath(data_root, mustWork = TRUE)
snapshot_version <- "gdc_dr45_retrieved_20260711"
dataset_key <- "GDC_TCGA_ESCA_clinical_gdc_dr45_retrieved_20260711"
retrieval_date <- "2026-07-11"
stage_root <- file.path(
  data_root, "_incoming", "GDC", "TCGA-ESCA-clinical", snapshot_version
)

if (dir.exists(stage_root) && length(list.files(stage_root, recursive = TRUE, all.files = TRUE)) > 0L) {
  stop("暂存目录非空，为避免覆盖请先人工核查：", stage_root)
}

metadata_dir <- file.path(stage_root, "10_metadata")
reusable_dir <- file.path(stage_root, "20_reusable")
manifest_dir <- file.path(stage_root, "90_manifests")
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(reusable_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

status_url <- "https://api.gdc.cancer.gov/status"
cases_url <- "https://api.gdc.cancer.gov/cases"

status_response <- request(status_url) |>
  req_user_agent("ESCC-multiomics-project/1.0") |>
  req_retry(max_tries = 5L) |>
  req_perform()
status_raw <- resp_body_raw(status_response)
status_parsed <- fromJSON(rawToChar(status_raw), simplifyVector = TRUE)

if (as.integer(status_parsed$data_release_version$major) != 45L ||
    as.integer(status_parsed$data_release_version$minor) != 0L ||
    !identical(status_parsed$data_release_version$release_date, "2025-12-04")) {
  stop("GDC live data release 已变化；不得沿用 DR45 版本名。")
}

case_filter <- list(
  op = "in",
  content = list(field = "project.project_id", value = list("TCGA-ESCA"))
)
case_fields <- c(
  "case_id",
  "submitter_id",
  "project.project_id",
  "demographic.vital_status",
  "demographic.days_to_death",
  "diagnoses.diagnosis_id",
  "diagnoses.diagnosis_is_primary_disease",
  "diagnoses.classification_of_tumor",
  "diagnoses.primary_diagnosis",
  "diagnoses.morphology",
  "diagnoses.ajcc_pathologic_stage",
  "diagnoses.ajcc_clinical_stage",
  "diagnoses.ajcc_pathologic_t",
  "diagnoses.ajcc_pathologic_n",
  "diagnoses.ajcc_pathologic_m",
  "diagnoses.tissue_or_organ_of_origin",
  "diagnoses.site_of_resection_or_biopsy",
  "diagnoses.days_to_diagnosis",
  "diagnoses.days_to_last_follow_up",
  "diagnoses.year_of_diagnosis",
  "diagnoses.tumor_grade"
)
query_parameters <- list(
  filters = case_filter,
  fields = case_fields,
  expand = c("diagnoses", "demographic", "project"),
  size = 1000L,
  format = "JSON"
)

cases_response <- request(cases_url) |>
  req_user_agent("ESCC-multiomics-project/1.0") |>
  req_url_query(
    filters = toJSON(case_filter, auto_unbox = TRUE),
    fields = paste(case_fields, collapse = ","),
    expand = paste(query_parameters$expand, collapse = ","),
    size = query_parameters$size,
    format = query_parameters$format
  ) |>
  req_retry(max_tries = 5L) |>
  req_perform()
cases_raw <- resp_body_raw(cases_response)
cases_parsed <- fromJSON(rawToChar(cases_raw), simplifyVector = FALSE)

status_path <- file.path(metadata_dir, "gdc_status_dr45.json")
response_path <- file.path(metadata_dir, "tcga_esca_gdc_cases_dr45.json")
query_path <- file.path(metadata_dir, "tcga_esca_gdc_cases_query_dr45.json")
writeBin(status_raw, status_path)
writeBin(cases_raw, response_path)
writeLines(
  toJSON(
    list(
      dataset_key = dataset_key,
      retrieval_date = retrieval_date,
      endpoint = cases_url,
      data_release = status_parsed$data_release,
      query = query_parameters
    ),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  ),
  query_path,
  useBytes = TRUE
)

hits <- cases_parsed$data$hits
if (length(hits) != 185L) {
  stop("TCGA-ESCA case 数不等于 DR45 预期的 185。")
}

scalar_value <- function(x, field) {
  value <- x[[field]]
  if (is.null(value) || length(value) == 0L) return(NA_character_)
  as.character(value[[1L]])
}

diagnosis_rows <- rbindlist(lapply(hits, function(hit) {
  diagnoses <- hit$diagnoses
  if (is.null(diagnoses) || length(diagnoses) == 0L) {
    diagnoses <- list(list())
  }
  rbindlist(lapply(seq_along(diagnoses), function(index) {
    diagnosis <- diagnoses[[index]]
    data.table(
      case_id = scalar_value(hit, "case_id"),
      case_submitter_id = scalar_value(hit, "submitter_id"),
      diagnosis_index = index,
      diagnosis_id = scalar_value(diagnosis, "diagnosis_id"),
      diagnosis_is_primary_disease = tolower(scalar_value(diagnosis, "diagnosis_is_primary_disease")) == "true",
      classification_of_tumor = scalar_value(diagnosis, "classification_of_tumor"),
      primary_diagnosis = scalar_value(diagnosis, "primary_diagnosis"),
      morphology = scalar_value(diagnosis, "morphology"),
      ajcc_pathologic_stage = scalar_value(diagnosis, "ajcc_pathologic_stage"),
      ajcc_clinical_stage = scalar_value(diagnosis, "ajcc_clinical_stage"),
      ajcc_pathologic_t = scalar_value(diagnosis, "ajcc_pathologic_t"),
      ajcc_pathologic_n = scalar_value(diagnosis, "ajcc_pathologic_n"),
      ajcc_pathologic_m = scalar_value(diagnosis, "ajcc_pathologic_m"),
      tissue_or_organ_of_origin = scalar_value(diagnosis, "tissue_or_organ_of_origin"),
      site_of_resection_or_biopsy = scalar_value(diagnosis, "site_of_resection_or_biopsy"),
      days_to_diagnosis = suppressWarnings(as.numeric(scalar_value(diagnosis, "days_to_diagnosis"))),
      days_to_last_follow_up = suppressWarnings(as.numeric(scalar_value(diagnosis, "days_to_last_follow_up"))),
      year_of_diagnosis = suppressWarnings(as.integer(scalar_value(diagnosis, "year_of_diagnosis"))),
      tumor_grade = scalar_value(diagnosis, "tumor_grade"),
      vital_status = scalar_value(hit$demographic, "vital_status"),
      days_to_death = suppressWarnings(as.numeric(scalar_value(hit$demographic, "days_to_death")))
    )
  }), fill = TRUE)
}), fill = TRUE)

diagnosis_rows[, selected_primary_diagnosis :=
  diagnosis_is_primary_disease & tolower(classification_of_tumor) == "primary"
]
primary_counts <- diagnosis_rows[, .(
  selected_primary_diagnosis_count = sum(selected_primary_diagnosis, na.rm = TRUE)
), by = .(case_id, case_submitter_id)]
if (any(primary_counts$selected_primary_diagnosis_count != 1L)) {
  stop("DR45 中存在无法唯一选择 primary disease diagnosis 的 case。")
}

primary_flat <- diagnosis_rows[selected_primary_diagnosis == TRUE]
primary_flat[, is_escc := morphology %in% c("8070/3", "8071/3", "8083/3")]
primary_flat[, escc_exclusion_reason := fifelse(
  is_escc,
  "",
  paste0("non_escc_morphology:", morphology)
)]
setorder(primary_flat, case_submitter_id)
setorder(diagnosis_rows, case_submitter_id, diagnosis_index)

if (nrow(primary_flat) != 185L || sum(primary_flat$is_escc) != 96L) {
  stop("DR45 primary diagnosis 或 ESCC 计数不符合 185/96 预期。")
}
morphology_check <- primary_flat[is_escc == TRUE, .N, by = morphology][order(morphology)]
expected_morphology <- data.table(
  morphology = c("8070/3", "8071/3", "8083/3"),
  N = c(90L, 5L, 1L)
)
if (!identical(morphology_check, expected_morphology)) {
  stop("DR45 ESCC morphology 计数不符合 90/5/1 预期。")
}

diagnoses_path <- file.path(reusable_dir, "tcga_esca_gdc_diagnoses_long_dr45.tsv")
flat_path <- file.path(reusable_dir, "tcga_esca_gdc_primary_diagnosis_flat_dr45.tsv")
summary_path <- file.path(reusable_dir, "tcga_esca_gdc_histology_summary_dr45.tsv")
provenance_path <- file.path(reusable_dir, "PROVENANCE.tsv")
source_metadata_path <- file.path(metadata_dir, "source_metadata.tsv")

fwrite(diagnosis_rows, diagnoses_path, sep = "\t", na = "NA", quote = FALSE)
fwrite(primary_flat, flat_path, sep = "\t", na = "NA", quote = FALSE)
fwrite(
  primary_flat[, .(cases = .N), by = .(is_escc, morphology, primary_diagnosis)][order(-is_escc, morphology)],
  summary_path,
  sep = "\t",
  na = "NA",
  quote = FALSE
)

source_sha <- digest(response_path, algo = "sha256", file = TRUE, serialize = FALSE)
source_metadata <- data.table(
  dataset_key = dataset_key,
  endpoint = cases_url,
  retrieval_date = retrieval_date,
  data_release = status_parsed$data_release,
  data_release_major = status_parsed$data_release_version$major,
  data_release_minor = status_parsed$data_release_version$minor,
  release_date = status_parsed$data_release_version$release_date,
  total_cases = nrow(primary_flat),
  escc_cases = sum(primary_flat$is_escc),
  source_response_sha256 = source_sha
)
fwrite(source_metadata, source_metadata_path, sep = "\t", na = "NA", quote = FALSE)

provenance <- data.table(
  output_file = c(
    "20_reusable/tcga_esca_gdc_diagnoses_long_dr45.tsv",
    "20_reusable/tcga_esca_gdc_primary_diagnosis_flat_dr45.tsv",
    "20_reusable/tcga_esca_gdc_histology_summary_dr45.tsv"
  ),
  input_file = "10_metadata/tcga_esca_gdc_cases_dr45.json",
  input_sha256 = source_sha,
  generation_script = file.path(project_root, "scripts", "03_stage_gdc_tcga_esca_clinical_dr45.R"),
  parameters = c(
    "expand all diagnoses; no row selection",
    "diagnosis_is_primary_disease=true; classification_of_tumor=Primary; ESCC morphology=8070/3|8071/3|8083/3",
    "aggregate selected primary diagnosis by morphology and text"
  ),
  software = paste0("R ", getRversion(), "; data.table; jsonlite; httr2"),
  generation_date = retrieval_date,
  reproducible = "yes"
)
fwrite(provenance, provenance_path, sep = "\t", na = "NA", quote = FALSE)

manifest_files <- c(
  status_path,
  response_path,
  query_path,
  source_metadata_path,
  diagnoses_path,
  flat_path,
  summary_path,
  provenance_path
)
relative_paths <- substring(manifest_files, nchar(stage_root) + 2L)
levels <- dirname(relative_paths)
source_urls <- c(
  status_url,
  cases_url,
  "generated query definition",
  "generated from GDC status and cases response",
  cases_url,
  cases_url,
  cases_url,
  "generated provenance record"
)
file_status <- c(
  "verified", "verified", "generated_verified", "generated_verified",
  "generated_verified", "generated_verified", "generated_verified", "generated_verified"
)
corresponding_source <- c(
  "", "", "", "10_metadata/gdc_status_dr45.json;10_metadata/tcga_esca_gdc_cases_dr45.json",
  "10_metadata/tcga_esca_gdc_cases_dr45.json",
  "10_metadata/tcga_esca_gdc_cases_dr45.json",
  "10_metadata/tcga_esca_gdc_cases_dr45.json",
  "10_metadata/tcga_esca_gdc_cases_dr45.json"
)
generation_method <- c(
  "source bytes retrieved without modification",
  "source bytes retrieved without modification",
  "script-encoded API query",
  "03_stage_gdc_tcga_esca_clinical_dr45.R",
  "03_stage_gdc_tcga_esca_clinical_dr45.R",
  "03_stage_gdc_tcga_esca_clinical_dr45.R",
  "03_stage_gdc_tcga_esca_clinical_dr45.R",
  "03_stage_gdc_tcga_esca_clinical_dr45.R"
)
manifest <- data.table(
  relative_path = relative_paths,
  file_level = levels,
  size_bytes = file.info(manifest_files)$size,
  sha256 = vapply(
    manifest_files,
    digest,
    character(1),
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  ),
  source_url = source_urls,
  download_date = retrieval_date,
  file_status = file_status,
  corresponding_source_file = corresponding_source,
  generation_method = generation_method,
  notes = c(
    "GDC live status snapshot",
    "185 TCGA-ESCA cases; all diagnosis objects preserved",
    "exact filter and fields",
    "release and source response summary",
    "all diagnoses; primary selection flag retained",
    "one primary disease diagnosis per case; 96 ESCC",
    "ESCC morphology distribution 90/5/1",
    "L2 generation provenance"
  )
)
fwrite(
  manifest,
  file.path(manifest_dir, "MANIFEST.tsv"),
  sep = "\t",
  na = "NA",
  quote = FALSE
)

message("暂存目录：", stage_root)
message("TCGA-ESCA cases：", nrow(primary_flat))
message("ESCC cases：", sum(primary_flat$is_escc))
message("ESCC morphology：", paste(morphology_check$morphology, morphology_check$N, collapse = "; "))
