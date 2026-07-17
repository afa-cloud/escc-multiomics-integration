#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(fs)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_root <- Sys.getenv("RESEARCH_DATA_ROOT", path.expand("~/ResearchDataHub"))
catalog_path <- file.path(data_root, "CATALOG.tsv")
stage_root <- file.path(
  data_root, "_incoming", "Metabolomics_Workbench", "PR001876",
  "retrieved_20260711"
)
canonical_root <- file.path(
  data_root, "datasets", "public", "Metabolomics_Workbench", "PR001876",
  "retrieved_20260711"
)
download_date <- "2026-07-11"
dataset_key <- "METABOLOMICS_WORKBENCH_PR001876_retrieved_20260711"

stopifnot(
  dir_exists(data_root),
  file_access(catalog_path, "read"),
  dir_exists(stage_root)
)
if (dir_exists(canonical_root)) stop("规范目录已存在，拒绝覆盖：", canonical_root)

partials <- dir_ls(
  stage_root,
  recurse = TRUE,
  type = "file",
  regexp = "\\.(aria2|part)$|\\.invalid\\."
)
if (length(partials)) stop("暂存目录仍有未完成或无效下载：", paste(partials, collapse = ";"))

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

write_tsv_atomic <- function(x, path) {
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = dirname(path))
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "")
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("无法原子更新文件：", path)
  }
}

raw_definitions <- data.table(
  study_id = c("ST003013", "ST003014", "ST003015", "ST003025", "ST003027"),
  filename = paste0(c("ST003013", "ST003014", "ST003015", "ST003025", "ST003027"), "_Rawfiles.zip"),
  expected_size = c(111464557, 88718718, 58144183, 637059, 46783597),
  expected_md5 = c(
    "044e1a26d894560a825c616773e8f580",
    "cf0a3c8ba0a5fa6dfcfaca877f6da5cf",
    "73c2626b9c28385c008d7dc9f3195d8f",
    "d444ca5d7b5dcf72a37b2f232596aa4f",
    "6b11b4d20538ecff1b808c41785b8b24"
  ),
  specimen = c("serum", "urine", "tissue", "tissue", "tissue"),
  platform = c("NMR", "NMR", "NMR", "targeted GC-MS", "targeted LC-MS")
)
raw_definitions[, path := file.path(stage_root, "00_source", filename)]
if (!all(file_exists(raw_definitions$path))) stop("原始 ZIP 不完整")

observed_size <- as.numeric(file_info(raw_definitions$path)$size)
observed_md5 <- unname(tools::md5sum(raw_definitions$path))
if (!all(observed_size == raw_definitions$expected_size)) stop("原始 ZIP 大小校验失败")
if (!all(tolower(observed_md5) == raw_definitions$expected_md5)) stop("原始 ZIP 官方 MD5 校验失败")

archive_inventory <- rbindlist(lapply(seq_len(nrow(raw_definitions)), function(i) {
  zip_path <- raw_definitions$path[i]
  if (!identical(system2("/usr/bin/unzip", c("-tq", shQuote(zip_path))), 0L)) {
    stop("ZIP 完整性校验失败：", zip_path)
  }
  members <- as.data.table(utils::unzip(zip_path, list = TRUE))
  setnames(members, c("Name", "Length", "Date"), c("member_name", "uncompressed_size_bytes", "member_date"))
  unsafe <- startsWith(members$member_name, "/") |
    vapply(strsplit(members$member_name, "/", fixed = TRUE), function(parts) {
      any(parts == "..")
    }, logical(1))
  if (any(unsafe)) stop("ZIP 含不安全成员路径：", zip_path)
  if (anyDuplicated(members$member_name)) stop("ZIP 含重复成员名：", zip_path)
  members[, `:=`(
    study_id = raw_definitions$study_id[i],
    archive_relative_path = file.path("00_source", raw_definitions$filename[i]),
    archive_member_status = "listed_from_verified_zip"
  )]
  setcolorder(members, c(
    "study_id", "archive_relative_path", "member_name",
    "uncompressed_size_bytes", "member_date", "archive_member_status"
  ))
  members
}), use.names = TRUE, fill = TRUE)
fwrite(
  archive_inventory,
  file.path(stage_root, "10_metadata", "archive_members.tsv"),
  sep = "\t", quote = FALSE, na = ""
)

analysis_definitions <- data.table(
  study_id = c("ST003013", "ST003014", "ST003015", "ST003025", "ST003027", "ST003027"),
  analysis_id = c("AN004946", "AN004947", "AN004948", "AN004960", "AN004962", "AN004963"),
  filename = c(
    "ST003013_AN004946_mwtab.txt",
    "ST003014_AN004947_mwtab.txt",
    "ST003015_AN004948_mwtab.txt",
    "ST003025_AN004960_mwtab.txt",
    "ST003027_AN004962_mwtab.txt",
    "ST003027_AN004963_mwtab.txt"
  ),
  specimen = c("serum", "urine", "tissue", "tissue", "tissue", "tissue"),
  platform = c("NMR", "NMR", "NMR", "targeted GC-MS", "targeted LC-MS", "targeted LC-MS"),
  ion_mode = c(NA, NA, NA, "source_recorded", "positive", "negative")
)
analysis_definitions[, path := file.path(stage_root, "10_metadata", filename)]
if (!all(file_exists(analysis_definitions$path))) stop("mwTab 元数据不完整")

extract_section <- function(lines, start_marker, end_marker) {
  start <- which(lines == start_marker)
  end <- which(lines == end_marker)
  if (length(start) != 1L || length(end) != 1L || end <= start) return(character())
  lines[(start + 1L):(end - 1L)]
}

parse_sample_factors <- function(path, study_id, analysis_id, specimen, platform, ion_mode) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  header <- lines[1L]
  expected_identity <- c(
    paste0("STUDY_ID:", study_id),
    paste0("ANALYSIS_ID:", analysis_id),
    "PROJECT_ID:PR001876"
  )
  if (!all(vapply(expected_identity, grepl, logical(1), x = header, fixed = TRUE))) {
    stop("mwTab 身份字段不符：", path)
  }

  factor_lines <- lines[startsWith(trimws(lines), "SUBJECT_SAMPLE_FACTORS")]
  split <- tstrsplit(factor_lines, "\t", fixed = TRUE, fill = "")
  factors <- data.table(
    tag = trimws(split[[1L]]),
    source_subject_id = trimws(split[[2L]]),
    sample_id = trimws(split[[3L]]),
    source_factor = trimws(split[[4L]]),
    raw_field = trimws(split[[5L]])
  )
  factors <- factors[tag == "SUBJECT_SAMPLE_FACTORS"]
  factors[, raw_file_name := sub("^RAW_FILE_NAME=", "", raw_field)]
  factors[, raw_index := suppressWarnings(as.integer(sub(
    ".*(?:Sample-|\\.)([0-9]+)(?:\\.fid|\\.mzML|\\.CDF)?$",
    "\\1", raw_file_name, perl = TRUE
  )))]
  factors[is.na(raw_index), raw_index := suppressWarnings(as.integer(sample_id))]
  factors[, `:=`(
    study_id = study_id,
    analysis_id = analysis_id,
    specimen = specimen,
    platform = platform,
    ion_mode = ion_mode,
    physical_specimen_key = paste(study_id, raw_file_name, sep = "::"),
    pairing_candidate_id = NA_character_,
    pairing_basis = "none",
    pairing_confidence = "none",
    paired_model_role = "not_paired"
  )]

  if (study_id == "ST003013") {
    factors[grepl("Pre-operation", source_factor), `:=`(
      pairing_candidate_id = sprintf("ST003013_ESCC_%02d", raw_index),
      pairing_basis = "pre/post numbering offset 54; source subject ID absent",
      pairing_confidence = "medium",
      paired_model_role = "sensitivity_only"
    )]
    factors[grepl("Post-operation", source_factor), `:=`(
      pairing_candidate_id = sprintf("ST003013_ESCC_%02d", raw_index - 54L),
      pairing_basis = "pre/post numbering offset 54; source subject ID absent",
      pairing_confidence = "medium",
      paired_model_role = "sensitivity_only"
    )]
  } else if (study_id == "ST003014") {
    factors[grepl("Pre-operation", source_factor), `:=`(
      pairing_candidate_id = sprintf("ST003014_ESCC_%02d", raw_index),
      pairing_basis = "pre/post numbering offset 54; source subject ID absent",
      pairing_confidence = "medium",
      paired_model_role = "sensitivity_only"
    )]
    factors[grepl("Post-operation", source_factor), `:=`(
      pairing_candidate_id = sprintf("ST003014_ESCC_%02d", raw_index - 54L),
      pairing_basis = "pre/post numbering offset 54; source subject ID absent",
      pairing_confidence = "medium",
      paired_model_role = "sensitivity_only"
    )]
  } else if (study_id == "ST003015") {
    factors[grepl("Tumor tissue", source_factor), `:=`(
      pairing_candidate_id = sprintf("ST003015_ESCC_%02d", raw_index),
      pairing_basis = "tumor/normal numbering offset 54; source subject ID absent",
      pairing_confidence = "medium",
      paired_model_role = "sensitivity_only"
    )]
    factors[grepl("Normal tissue", source_factor), `:=`(
      pairing_candidate_id = sprintf("ST003015_ESCC_%02d", raw_index - 54L),
      pairing_basis = "tumor/normal numbering offset 54; source subject ID absent",
      pairing_confidence = "medium",
      paired_model_role = "sensitivity_only"
    )]
  } else if (study_id %in% c("ST003025", "ST003027")) {
    factors[grepl("Early stage ESCC", source_factor), `:=`(
      pairing_candidate_id = sprintf(paste0(study_id, "_candidate_%02d"), raw_index),
      pairing_basis = "case/control numbering offset 16; participant mapping not supplied",
      pairing_confidence = "low",
      paired_model_role = "unpaired_primary_only"
    )]
    factors[grepl("Normal tissue", source_factor), `:=`(
      pairing_candidate_id = sprintf(paste0(study_id, "_candidate_%02d"), raw_index - 16L),
      pairing_basis = "case/control numbering offset 16; participant mapping not supplied",
      pairing_confidence = "low",
      paired_model_role = "unpaired_primary_only"
    )]
  }

  factors[, source_subject_id := fifelse(source_subject_id == "-", NA_character_, source_subject_id)]
  factors[, raw_field := NULL]
  setcolorder(factors, c(
    "study_id", "analysis_id", "source_subject_id", "sample_id", "source_factor",
    "raw_file_name", "raw_index", "specimen", "platform", "ion_mode",
    "physical_specimen_key", "pairing_candidate_id", "pairing_basis",
    "pairing_confidence", "paired_model_role"
  ))
  list(lines = lines, factors = factors)
}

dir_create(file.path(stage_root, "20_reusable"), recurse = TRUE)
sample_maps <- list()
generated_matrix_files <- character()

for (i in seq_len(nrow(analysis_definitions))) {
  definition <- analysis_definitions[i]
  parsed <- parse_sample_factors(
    definition$path,
    definition$study_id,
    definition$analysis_id,
    definition$specimen,
    definition$platform,
    definition$ion_mode
  )
  sample_maps[[definition$analysis_id]] <- parsed$factors

  data_lines <- extract_section(
    parsed$lines,
    "MS_METABOLITE_DATA_START",
    "MS_METABOLITE_DATA_END"
  )
  metabolite_lines <- extract_section(
    parsed$lines,
    "METABOLITES_START",
    "METABOLITES_END"
  )

  if (length(data_lines)) {
    header <- strsplit(data_lines[1L], "\t", fixed = TRUE)[[1L]]
    if (length(data_lines) < 3L || !startsWith(data_lines[2L], "Factors\t")) {
      stop("MS_METABOLITE_DATA 缺少 Factors 行：", definition$path)
    }
    matrix <- fread(
      text = paste(data_lines[-c(1L, 2L)], collapse = "\n"),
      sep = "\t", header = FALSE, fill = TRUE, quote = "",
      na.strings = c("NA", ""), encoding = "UTF-8"
    )
    if (ncol(matrix) != length(header)) stop("代谢物矩阵列数与表头不一致")
    setnames(matrix, header)
    setnames(matrix, 1L, "metabolite_name")
    matrix_path <- file.path(
      stage_root, "20_reusable",
      paste0(definition$analysis_id, "_metabolite_matrix.tsv")
    )
    fwrite(matrix, matrix_path, sep = "\t", quote = FALSE, na = "NA")
    generated_matrix_files <- c(generated_matrix_files, matrix_path)
  }

  if (length(metabolite_lines)) {
    metabolites <- fread(
      text = paste(metabolite_lines, collapse = "\n"),
      sep = "\t", header = TRUE, fill = TRUE, quote = "",
      na.strings = c("NA", ""), encoding = "UTF-8"
    )
    metabolite_path <- file.path(
      stage_root, "20_reusable",
      paste0(definition$analysis_id, "_metabolite_annotations.tsv")
    )
    fwrite(metabolites, metabolite_path, sep = "\t", quote = FALSE, na = "")
    generated_matrix_files <- c(generated_matrix_files, metabolite_path)
  }
}

sample_map <- rbindlist(sample_maps, use.names = TRUE, fill = TRUE)
sample_map_path <- file.path(stage_root, "20_reusable", "pr001876_analysis_sample_map.tsv")
fwrite(sample_map, sample_map_path, sep = "\t", quote = FALSE, na = "")

study_summary <- sample_map[, .(
  analysis_rows = .N,
  physical_specimens = uniqueN(physical_specimen_key),
  source_subject_ids_present = sum(!is.na(source_subject_id)),
  group_summary = paste(names(table(source_factor)), as.integer(table(source_factor)), collapse = ";"),
  pairing_confidence = paste(sort(unique(pairing_confidence)), collapse = ";"),
  paired_model_role = paste(sort(unique(paired_model_role)), collapse = ";")
), by = .(study_id, analysis_id, specimen, platform, ion_mode)]
study_summary[, processed_matrix_available := analysis_id %in% c("AN004960", "AN004962", "AN004963")]
study_summary_path <- file.path(stage_root, "20_reusable", "pr001876_study_summary.tsv")
fwrite(study_summary, study_summary_path, sep = "\t", quote = FALSE, na = "")

provenance <- data.table(
  output_scope = c(
    "10_metadata/archive_members.tsv",
    "20_reusable/pr001876_analysis_sample_map.tsv",
    "20_reusable/pr001876_study_summary.tsv",
    "20_reusable/AN004960_metabolite_matrix.tsv and annotations",
    "20_reusable/AN004962/AN004963 metabolite matrices and annotations"
  ),
  corresponding_source_file = c(
    "00_source/*_Rawfiles.zip",
    "10_metadata/*_mwtab.txt",
    "10_metadata/*_mwtab.txt",
    "10_metadata/ST003025_AN004960_mwtab.txt",
    "10_metadata/ST003027_AN004962_mwtab.txt;10_metadata/ST003027_AN004963_mwtab.txt"
  ),
  generating_script = "scripts/07_ingest_pr001876_metabolomics.R",
  key_parameters = c(
    "list verified ZIP members without extraction",
    "preserve source IDs; inferred numbering pairs marked sensitivity-only or unpaired",
    "summarize actual analysis rows and physical specimen keys",
    "extract source-provided wide matrix and metabolite annotations",
    "keep positive and negative analyses as two modes of the same ST003027 sample set"
  ),
  software = paste0("R ", getRversion(), "; data.table ", packageVersion("data.table")),
  generated_date = download_date,
  regenerable = "yes"
)
provenance_path <- file.path(stage_root, "20_reusable", "PROVENANCE.tsv")
fwrite(provenance, provenance_path, sep = "\t", quote = FALSE, na = "")

dataset_md <- c(
  "# PR001876 ESCC 代谢组数据说明",
  "",
  "## 基本信息",
  "",
  paste0("- `dataset_key`：`", dataset_key, "`。"),
  "- 来源：美国国立卫生研究院代谢组学工作台（Metabolomics Workbench）。",
  "- 项目：`PR001876`；子研究：`ST003013`、`ST003014`、`ST003015`、`ST003025`、`ST003027`。",
  "- 本地冻结版本：`retrieved_20260711`；获取日期：2026-07-11。",
  "- 物种：人（Homo sapiens）；疾病：食管鳞状细胞癌（ESCC）。",
  "- 组学层：血清、尿液和组织 NMR；早期 ESCC 组织靶向 GC-MS/LC-MS。",
  "- 当前完整性状态：`verified`。五个原始 ZIP 已通过官方大小、MD5、ZIP 容器完整性和本地 SHA256；六个 mwTab 已核对项目、研究和分析 ID。",
  "",
  "## 目录与主要文件",
  "",
  "- `00_source/*_Rawfiles.zip`：五个子研究的来源原始包，不原位解压或改写。",
  "- `10_metadata/*_mwtab.txt`：六个官方 mwTab 快照。",
  "- `10_metadata/archive_members.tsv`：五个 ZIP 的成员路径、解压后大小与安全路径审计。",
  "- `20_reusable/pr001876_analysis_sample_map.tsv`：分析级样本、体液/组织、分组和配对边界。",
  "- `20_reusable/pr001876_study_summary.tsv`：各 study/analysis 的实际样本行、物理样本数和可用矩阵状态。",
  "- `20_reusable/AN004960/AN004962/AN004963_*`：mwTab 中公开的靶向 MS 代谢物矩阵与注释。",
  "- `20_reusable/PROVENANCE.tsv`：生成关系和参数。",
  "- `90_manifests/MANIFEST.tsv`：逐文件大小、SHA256、来源与完整性状态。",
  "",
  "## 推荐使用场景",
  "",
  "- 作为 ESCC 代谢通路、体液方向和早期组织代谢表型的独立正交证据。",
  "- ST003025 与 ST003027 的来源处理矩阵可用于组间差异、正负离子模式一致性和通路层校准。",
  "- NMR 原始包可用于后续可复现预处理；当前 mwTab 未提供 NMR 定量矩阵。",
  "",
  "## 配对、纳入排除与证据边界",
  "",
  "- mwTab 的 `SUBJECT` 字段均为空；不能把编号直接升级为已确认参与者 ID。",
  "- ST003013/ST003014 的 54 例术前/54 例术后以及 ST003015 的 54 例肿瘤/54 例正常可按编号偏移 54 建立候选配对，但在补充材料确认前只用于配对敏感性分析。",
  "- ST003025/ST003027 的 16 例早期 ESCC 与 16 例正常编号偏移 16 仅是低置信候选关系；主分析使用非配对组间模型。",
  "- ST003027 的 AN004962 与 AN004963 共享同一组 32 个物理样本，是正/负模式，不计为两个独立队列。",
  "- 血清和尿液健康对照数量不同，且无参与者 ID；不得假设跨体液同人。",
  "- 本队列与基因组、转录组和蛋白组主队列不是同一患者，只在代谢物、反应、通路或模块层连接。",
  "",
  "## 来源、引用与许可",
  "",
  "- 项目页面：https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR001876",
  "- 使用时应引用 Metabolomics Workbench 项目页与对应原始研究，并遵守来源条款。",
  paste0("- 使用项目：`", project_root, "`。"),
  "- `verified` 仅表示当前快照字节与来源身份完整，不证明样本配对或生物学因果。"
)
writeLines(dataset_md, file.path(stage_root, "DATASET.md"), useBytes = TRUE)

source_files <- c(raw_definitions$path, analysis_definitions$path)
generated_files <- c(
  file.path(stage_root, "10_metadata", "archive_members.tsv"),
  sample_map_path,
  study_summary_path,
  generated_matrix_files,
  provenance_path
)
manifest_files <- c(source_files, generated_files)
relative_paths <- path_rel(manifest_files, start = stage_root)
manifest <- data.table(
  relative_path = relative_paths,
  file_level = sub("/.*$", "", relative_paths),
  size_bytes = as.numeric(file_info(manifest_files)$size),
  sha256 = vapply(manifest_files, sha256_file, character(1)),
  source_url = "generated from verified Metabolomics Workbench files",
  download_date = download_date,
  file_status = "generated_verified",
  corresponding_source_file = "10_metadata/*_mwtab.txt",
  generation_method = "scripts/07_ingest_pr001876_metabolomics.R",
  notes = ""
)

for (i in seq_len(nrow(raw_definitions))) {
  relative <- file.path("00_source", raw_definitions$filename[i])
  manifest[relative_path == relative, `:=`(
    source_url = paste0(
      "https://www.metabolomicsworkbench.org/studydownload/",
      raw_definitions$filename[i]
    ),
    file_status = "verified",
    corresponding_source_file = "",
    generation_method = "source bytes retrieved without modification",
    notes = paste0("official_md5=", raw_definitions$expected_md5[i], ";zip_integrity_ok")
  )]
}
for (i in seq_len(nrow(analysis_definitions))) {
  relative <- file.path("10_metadata", analysis_definitions$filename[i])
  manifest[relative_path == relative, `:=`(
    source_url = paste0(
      "https://www.metabolomicsworkbench.org/data/study_textformat_view.php?ANALYSIS_ID=",
      analysis_definitions$analysis_id[i],
      "&MODE=d&STUDY_ID=", analysis_definitions$study_id[i]
    ),
    file_status = "verified",
    corresponding_source_file = "",
    generation_method = "source bytes retrieved without modification",
    notes = "project study and analysis identifiers verified"
  )]
}
manifest[relative_path == "10_metadata/archive_members.tsv", `:=`(
  corresponding_source_file = "00_source/*_Rawfiles.zip",
  notes = "member inventory only; raw archives not extracted"
)]
setorder(manifest, relative_path)
fwrite(
  manifest,
  file.path(stage_root, "90_manifests", "MANIFEST.tsv"),
  sep = "\t", quote = FALSE, na = ""
)

dir_create(dirname(canonical_root), recurse = TRUE)
dir_copy(stage_root, canonical_root, overwrite = FALSE)
canonical_manifest <- fread(
  file.path(canonical_root, "90_manifests", "MANIFEST.tsv"),
  colClasses = "character"
)
canonical_files <- file.path(canonical_root, canonical_manifest$relative_path)
if (!all(file_exists(canonical_files))) stop("规范路径复制后存在缺失文件")
observed_sizes <- as.character(as.numeric(file_info(canonical_files)$size))
observed_sha <- vapply(canonical_files, sha256_file, character(1))
if (!all(observed_sizes == canonical_manifest$size_bytes) ||
    !all(observed_sha == canonical_manifest$sha256)) {
  stop("规范路径复制后大小或 SHA256 校验失败")
}

catalog <- fread(catalog_path, colClasses = "character", na.strings = NULL)
if (dataset_key %in% catalog$dataset_key) stop("CATALOG.tsv 已存在 dataset_key")
catalog_row <- data.table(
  dataset_key = dataset_key,
  record_type = "multi_biofluid_metabolomics_dataset",
  access_level = "public",
  source = "Metabolomics Workbench",
  accession = "PR001876;ST003013;ST003014;ST003015;ST003025;ST003027",
  version = "retrieved_20260711",
  species = "Homo sapiens",
  disease = "esophageal squamous cell carcinoma",
  tissue = "serum;urine;primary tumor tissue;normal tissue",
  assay = "NMR;targeted GC-MS;targeted LC-MS",
  sample_summary = "mwTab 实际记录 195 血清、147 尿液、108 NMR 组织样本；靶向 MS 为 32 个组织样本，ST003027 两种离子模式共享样本",
  available_levels = "00_source;10_metadata;20_reusable",
  status = "verified",
  local_path = canonical_root,
  manifest_path = file.path(canonical_root, "90_manifests", "MANIFEST.tsv"),
  source_url = "https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR001876",
  download_date = download_date,
  last_verified = download_date,
  license_or_access = "public; cite Metabolomics Workbench and original study; follow source terms",
  projects_using = project_root,
  recommended_use = "ESCC 代谢通路和多体液方向校准；靶向 MS 组织矩阵用于独立组间验证",
  limitations = "source subject IDs absent; numbering-derived pairs are sensitivity-only or unpaired; NMR mwTab lacks quantitative matrix; not matched to genomic cohorts",
  notes = "ST003027 positive and negative analyses are two modes of one 32-sample cohort and cannot be counted independently."
)
catalog <- rbindlist(list(catalog, catalog_row), use.names = TRUE, fill = FALSE)
write_tsv_atomic(catalog, catalog_path)

project_datasets_path <- file.path(project_root, "data", "datasets.tsv")
project_datasets <- fread(project_datasets_path, colClasses = "character", na.strings = NULL)
if (dataset_key %in% project_datasets$dataset_key) stop("data/datasets.tsv 已存在 dataset_key")
project_row <- data.table(
  logical_name = "escc_metabolomics_pr001876",
  dataset_key = dataset_key,
  accession = "PR001876;ST003013;ST003014;ST003015;ST003025;ST003027",
  version = "retrieved_20260711",
  data_level = "00_source;10_metadata;20_reusable",
  central_path = canonical_root,
  project_purpose = "代谢通路、多体液方向和早期组织靶向代谢验证",
  inclusion_exclusion = "靶向 MS 主分析按非配对组间设计；编号推定配对仅作敏感性分析；正负模式不重复计为独立队列",
  limitations = "无可确认 subject ID；NMR 尚需从原始数据重处理；与主分子队列无同患者映射",
  local_link = "",
  last_verified = download_date
)
project_datasets <- rbindlist(
  list(project_datasets, project_row),
  use.names = TRUE,
  fill = FALSE
)
write_tsv_atomic(project_datasets, project_datasets_path)

message("已完成 PR001876 原始 ZIP/mwTab 审计、靶向 MS 矩阵提取、规范路径复制、全 SHA 校验和索引登记。")
