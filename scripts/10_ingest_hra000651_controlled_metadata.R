#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(filelock)
  library(fs)
  library(jsonlite)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_root <- Sys.getenv("RESEARCH_DATA_ROOT", path.expand("~/ResearchDataHub"))
catalog_path <- file.path(data_root, "CATALOG.tsv")
project_datasets_path <- file.path(project_root, "data", "datasets.tsv")
download_date <- "2026-07-11"
dataset_key <- "GSA_HUMAN_HRA000651_retrieved_20260711"
stage_root <- file.path(
  data_root, "_incoming", "GSA-Human", "HRA000651", "retrieved_20260711"
)
canonical_root <- file.path(
  data_root, "datasets", "controlled", "GSA-Human", "HRA000651",
  "retrieved_20260711"
)

if (!file_exists(file.path(project_root, "PROJECT_INDEX.md")) ||
    !file_exists(project_datasets_path)) {
  stop("必须从项目根目录运行，且 PROJECT_INDEX.md 与 data/datasets.tsv 必须存在")
}
if (format(Sys.Date(), "%Y-%m-%d") != download_date) {
  stop("快照日期与实际运行日期不一致；请先显式升级版本目录和 download_date")
}
stopifnot(dir_exists(data_root), file_access(catalog_path, "read"))
if (dir_exists(canonical_root)) stop("规范目录已存在，拒绝覆盖：", canonical_root)

catalog <- fread(catalog_path, colClasses = "character", na.strings = NULL)
if (dataset_key %in% catalog$dataset_key) stop("CATALOG.tsv 已存在 dataset_key")
project_datasets_preflight <- fread(
  project_datasets_path, colClasses = "character", na.strings = NULL
)
required_project_columns <- c(
  "logical_name", "dataset_key", "accession", "version", "data_level",
  "central_path", "project_purpose", "inclusion_exclusion", "limitations",
  "local_link", "last_verified"
)
if (!identical(names(project_datasets_preflight), required_project_columns)) {
  stop("data/datasets.tsv schema 不符合预期")
}
if (dataset_key %in% project_datasets_preflight$dataset_key) {
  stop("data/datasets.tsv 已存在 dataset_key")
}

for (relative in c("00_source", "10_metadata", "20_reusable", "90_manifests")) {
  dir_create(file.path(stage_root, relative), recurse = TRUE)
}

individual_url <- paste0(
  "https://ngdc.cncb.ac.cn/gsa-human/ajaxb/indinstudy?",
  "accession=HRA000651&pageNo=1&pageSize=200"
)
run_url <- paste0(
  "https://ngdc.cncb.ac.cn/gsa-human/ajaxb/runinstudy?",
  "accession=HRA000651&pageNo=1&pageSize=200"
)
study_url <- "https://ngdc.cncb.ac.cn/gsa-human/browse/HRA000651"
paper_url <- "https://pmc.ncbi.nlm.nih.gov/articles/PMC7930383/"

download_atomic <- function(url, output_path) {
  part_path <- paste0(output_path, ".part")
  if (file_exists(part_path)) file_delete(part_path)
  status <- tryCatch(
    download.file(url, part_path, mode = "wb", method = "libcurl", quiet = TRUE),
    error = function(e) 1L
  )
  if (!identical(status, 0L) || !file_exists(part_path) || file_info(part_path)$size < 100) {
    stop("元数据下载失败：", url)
  }
  if (file_exists(output_path)) file_delete(output_path)
  file_move(part_path, output_path)
  if (!file_exists(output_path)) stop("无法发布元数据文件：", output_path)
}

individual_json_path <- file.path(
  stage_root, "10_metadata", "gsa_human_individuals_HRA000651_20260711.json"
)
run_json_path <- file.path(
  stage_root, "10_metadata", "gsa_human_runs_HRA000651_20260711.json"
)
download_atomic(individual_url, individual_json_path)
download_atomic(run_url, run_json_path)

individual_response <- fromJSON(individual_json_path, simplifyDataFrame = TRUE)
run_response <- fromJSON(run_json_path, simplifyDataFrame = TRUE)
individuals <- as.data.table(individual_response$individualViews)
runs <- as.data.table(run_response$runViews)

if (nrow(individuals) != 53L || uniqueN(individuals$accession) != 53L) {
  stop("HRA000651 不是 53 个唯一个体")
}
if (nrow(runs) != 106L || uniqueN(runs$runAcc) != 53L) {
  stop("HRA000651 不是 53 个 run/106 个文件行")
}
run_file_counts <- runs[, .N, by = runAcc]
if (!all(run_file_counts$N == 2L)) stop("每个 run 不是两个方向文件")
if (!all(runs$platform == "Illumina HiSeq 2500")) stop("平台不一致")

individuals[, tissue_role := fifelse(
  grepl("PN$", sampleName),
  "physiological_normal",
  fifelse(grepl("T$", sampleName), "escc_tumor", "unresolved")
)]
if (sum(individuals$tissue_role == "escc_tumor") != 38L ||
    sum(individuals$tissue_role == "physiological_normal") != 15L ||
    any(individuals$tissue_role == "unresolved")) {
  stop("肿瘤/生理正常构成不是 38/15")
}

individual_map <- individuals[, .(
  individual_accession = accession,
  individual_identifier = name,
  gender,
  sample_accession = sampleAcc,
  biosample_accession = biosampleAcc,
  sample_name = sampleName,
  sample_title = sampleTitle,
  tissue_role,
  pairing_status = "unpaired_independent_individual",
  cohort_assignment = "not_resolved_from_public_repository_metadata",
  sample_description = sampleDesc
)]
setorder(individual_map, tissue_role, sample_name)

run_map <- merge(
  runs,
  individual_map[, .(
    individual_accession, individual_identifier, sample_accession,
    biosample_accession, sample_name, tissue_role, pairing_status
  )],
  by.x = c("sampleAcc", "biosampleAcc"),
  by.y = c("sample_accession", "biosample_accession"),
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(run_map$individual_accession)) stop("run 到个体映射不完整")
run_map <- run_map[, .(
  individual_accession,
  individual_identifier,
  sample_accession = sampleAcc,
  biosample_accession = biosampleAcc,
  sample_name,
  tissue_role,
  pairing_status,
  experiment_accession = expName,
  run_accession = runAcc,
  run_title = runTitle,
  platform,
  file_name = runFileName,
  repository_reported_file_size = runFileSize,
  repository_reported_file_size_unit = "GB",
  raw_access_status = "controlled_pending"
)]
setorder(run_map, run_accession, file_name)

individual_map_path <- file.path(
  stage_root, "10_metadata", "hra000651_individual_sample_map.tsv"
)
run_map_path <- file.path(stage_root, "10_metadata", "hra000651_run_file_map.tsv")
fwrite(individual_map, individual_map_path, sep = "\t", quote = FALSE, na = "")
fwrite(run_map, run_map_path, sep = "\t", quote = FALSE, na = "")

study_context <- data.table(
  accession = "HRA000651",
  repository_access_status = "controlled_pending",
  repository_observation_date = download_date,
  repository_url = study_url,
  metadata_endpoints_public = "yes",
  raw_download_attempted = "no",
  access_boundary = paste(
    "当前 GSA-Human 仓库页面标记受控访问；",
    "论文早期的 publicly accessible 表述不替代当前仓库状态"
  ),
  dac_status = "not_applied",
  project_use_status = "metadata_only_not_analysis_eligible"
)
study_context_path <- file.path(
  stage_root, "10_metadata", "hra000651_access_context.tsv"
)
fwrite(study_context, study_context_path, sep = "\t", quote = FALSE, na = "")

publication_design <- data.table(
  cohort = c("discovery", "validation"),
  escc_tumor_individuals = c(18L, 20L),
  physiological_normal_individuals = c(11L, 4L),
  sample_level_cohort_mapping = "not_resolved_from_public_repository_metadata",
  marker_region = "16S rRNA V4",
  primer_forward = "515F GTGCCAGCMGCCGCGGTAA",
  primer_reverse = "806R GGACTACHVGGGTWTCTAAT",
  sequencing = "Illumina HiSeq 2500 2x250 bp",
  negative_control_protocol = "H2O negative control reported for sequencing in each sample",
  negative_control_raw_available = "no_independent_blank_run_identified",
  source_url = paper_url
)
publication_design_path <- file.path(
  stage_root, "10_metadata", "hra000651_publication_design.tsv"
)
fwrite(publication_design, publication_design_path, sep = "\t", quote = FALSE, na = "")

provenance <- data.table(
  output_scope = c(
    "10_metadata/hra000651_individual_sample_map.tsv",
    "10_metadata/hra000651_run_file_map.tsv",
    "10_metadata/hra000651_access_context.tsv",
    "10_metadata/hra000651_publication_design.tsv"
  ),
  corresponding_source_file = c(
    "10_metadata/gsa_human_individuals_HRA000651_20260711.json",
    "10_metadata/gsa_human_runs_HRA000651_20260711.json",
    study_url,
    paper_url
  ),
  generating_script = "scripts/10_ingest_hra000651_controlled_metadata.R",
  key_parameters = c(
    "53 unique individuals; 38 T and 15 PN; no pairing inferred from shared numeric labels",
    "53 unique HRR runs; 106 file rows; two file rows per run",
    "metadata-only registration; raw access controlled_pending; no raw download attempted",
    "discovery 18T+11PN; validation 20T+4PN; cohort-to-sample map unresolved"
  ),
  software = paste0(
    "R ", getRversion(), "; data.table ", packageVersion("data.table"),
    "; jsonlite ", packageVersion("jsonlite")
  ),
  generated_date = download_date,
  regenerable = "yes"
)
provenance_path <- file.path(stage_root, "20_reusable", "PROVENANCE.tsv")
fwrite(provenance, provenance_path, sep = "\t", quote = FALSE, na = "")

dataset_md <- c(
  "# HRA000651 ESCC 食管组织 16S 受控数据说明",
  "",
  "## 基本信息",
  "",
  paste0("- `dataset_key`：`", dataset_key, "`。"),
  "- 来源：GSA-Human；accession `HRA000651`。",
  "- 本地元数据冻结版本：`retrieved_20260711`。",
  "- 当前状态：`incomplete` / `controlled_pending`；未获批、未下载原始 FASTQ，不可进入正式分析。",
  "- 公开元数据：53 位独立个体、53 个 HRR run、106 个文件行；38 例 ESCC 肿瘤、15 例生理正常。",
  "- 设计为非配对队列；数字编号相同的 T 与 PN 也是不同 HRI 个体，不得强行配对。",
  "",
  "## 目录与文件",
  "",
  "- `00_source/`：受控批准前保持为空。",
  "- `10_metadata/gsa_human_*json`：GSA-Human 公开个体与 run 接口原始快照。",
  "- `10_metadata/hra000651_individual_sample_map.tsv`：个体、样本、角色与非配对边界。",
  "- `10_metadata/hra000651_run_file_map.tsv`：HRX/HRR 和 106 个受控文件行。",
  "- `10_metadata/hra000651_access_context.tsv`：当前访问边界和 DAC 状态。",
  "- `10_metadata/hra000651_publication_design.tsv`：论文队列、引物、平台和阴性对照说明。",
  "",
  "## 推荐使用与限制",
  "",
  "- 获批后在本队列内独立完成 V4 ASV 与污染敏感性分析。",
  "- 与 PRJNA766558 的 V3-V4/MiSeq/FFPE 数据不直接拼接 ASV 表；只在属、科或功能模块层校准。",
  "- 论文报告 H2O 阴性对照，但当前 53 个公开元数据 run 中未识别独立 blank run。",
  "- 与主分子队列不是同患者，不做伪样本级相关。",
  "",
  "## 访问和引用",
  "",
  paste0("- 当前仓库页：", study_url),
  paste0("- 原始论文：", paper_url),
  "- 原始 FASTQ 须经当前 GSA-Human DAC 流程批准；未经授权不尝试获取。"
)
dataset_md_path <- file.path(stage_root, "DATASET.md")
writeLines(dataset_md, dataset_md_path, useBytes = TRUE)

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}
manifest_files <- c(
  individual_json_path, run_json_path, individual_map_path, run_map_path,
  study_context_path, publication_design_path, provenance_path, dataset_md_path
)
manifest_relative <- path_rel(manifest_files, start = stage_root)
manifest <- data.table(
  relative_path = manifest_relative,
  file_level = sub("/.*$", "", manifest_relative),
  size_bytes = as.numeric(file_info(manifest_files)$size),
  sha256 = vapply(manifest_files, sha256_file, character(1)),
  source_url = c(
    individual_url, run_url, individual_url, run_url,
    study_url, paper_url, "generated from registered source metadata",
    "generated from registered source metadata"
  ),
  download_date = download_date,
  file_status = c("verified", "verified", rep("generated_verified", 6L)),
  corresponding_source_file = c(
    "", "",
    "10_metadata/gsa_human_individuals_HRA000651_20260711.json",
    "10_metadata/gsa_human_runs_HRA000651_20260711.json",
    study_url, paper_url,
    "10_metadata/gsa_human_individuals_HRA000651_20260711.json;10_metadata/gsa_human_runs_HRA000651_20260711.json",
    "10_metadata/hra000651_access_context.tsv;10_metadata/hra000651_publication_design.tsv"
  ),
  generation_method = c(
    "GSA-Human public metadata API raw response",
    "GSA-Human public metadata API raw response",
    rep("scripts/10_ingest_hra000651_controlled_metadata.R", 6L)
  ),
  notes = c(
    "public metadata only; not raw sequence",
    "public metadata only; file rows remain controlled",
    "53 unique independent individuals",
    "53 runs and 106 controlled file rows",
    "current repository access boundary; no raw download attempted",
    "publication-level design; sample-to-cohort assignment unresolved",
    "",
    "controlled metadata dataset documentation"
  )
)
setorder(manifest, relative_path)
manifest_path <- file.path(stage_root, "90_manifests", "MANIFEST.tsv")
fwrite(manifest, manifest_path, sep = "\t", quote = FALSE, na = "")

expected_stage_files <- sort(as.character(c(
  path_rel(manifest_files, start = stage_root),
  "90_manifests/MANIFEST.tsv"
)))
observed_stage_files <- sort(as.character(path_rel(
  dir_ls(stage_root, recurse = TRUE, type = "file", all = TRUE),
  start = stage_root
)))
if (!identical(observed_stage_files, expected_stage_files)) {
  stop(
    "_incoming 存在未登记或缺失文件：",
    paste(setdiff(union(observed_stage_files, expected_stage_files),
                  intersect(observed_stage_files, expected_stage_files)), collapse = ";")
  )
}

dir_create(dirname(canonical_root), recurse = TRUE)
publish_root <- paste0(canonical_root, ".publishing-", Sys.getpid())
if (dir_exists(publish_root)) stop("临时发布目录已存在：", publish_root)
on.exit({
  if (dir_exists(publish_root)) dir_delete(publish_root)
}, add = TRUE)
dir_copy(stage_root, publish_root, overwrite = FALSE)
canonical_manifest <- fread(
  file.path(publish_root, "90_manifests", "MANIFEST.tsv"),
  colClasses = "character"
)
canonical_files <- file.path(publish_root, canonical_manifest$relative_path)
if (!all(file_exists(canonical_files))) stop("规范路径复制后存在缺失文件")
if (!all(as.character(as.numeric(file_info(canonical_files)$size)) == canonical_manifest$size_bytes) ||
    !all(vapply(canonical_files, sha256_file, character(1)) == canonical_manifest$sha256)) {
  stop("规范路径复制后大小或 SHA256 校验失败")
}
observed_publish_files <- sort(as.character(path_rel(
  dir_ls(publish_root, recurse = TRUE, type = "file", all = TRUE),
  start = publish_root
)))
if (!identical(observed_publish_files, expected_stage_files)) {
  stop("临时发布目录文件集合与登记白名单不一致")
}
file_move(publish_root, canonical_root)

catalog_row <- data.table(
  dataset_key = dataset_key,
  record_type = "controlled_amplicon_metadata",
  access_level = "controlled",
  source = "GSA-Human",
  accession = "HRA000651",
  version = "retrieved_20260711",
  species = "Homo sapiens",
  disease = "esophageal squamous cell carcinoma",
  tissue = "esophageal tumor; physiological normal esophagus",
  assay = "16S rRNA V4; Illumina HiSeq 2500 paired-end 250 bp",
  sample_summary = "53 independent individuals; 38 ESCC tumor and 15 physiological normal; 53 HRR runs; metadata only",
  available_levels = "10_metadata",
  status = "incomplete",
  local_path = canonical_root,
  manifest_path = file.path(canonical_root, "90_manifests", "MANIFEST.tsv"),
  source_url = study_url,
  download_date = download_date,
  last_verified = download_date,
  license_or_access = "controlled access; DAC approval required before raw download",
  projects_using = project_root,
  recommended_use = "metadata planning only until controlled access approval",
  limitations = "raw FASTQ unavailable locally; unpaired design; cohort-to-sample map unresolved; no independent blank run identified",
  notes = "Paper described public accessibility, but current repository access status governs; no unauthorized raw retrieval attempted."
)
project_row <- data.table(
  logical_name = "escc_tissue_microbiome_hra000651_controlled",
  dataset_key = dataset_key,
  accession = "HRA000651",
  version = "retrieved_20260711",
  data_level = "10_metadata",
  central_path = canonical_root,
  project_purpose = "受控微生物队列元数据规划与未来独立方向校准",
  inclusion_exclusion = "当前仅元数据；DAC 批准和原始完整性校验前不进分析",
  limitations = "原始 FASTQ 受控未取得；非配对设计；发现/验证到样本映射未解析；无独立 blank run",
  local_link = "",
  last_verified = download_date
)

# 全局目录和项目引用必须串行更新；持锁期间重读两份当前文件，避免丢失其他入库记录。
registry_lock <- filelock::lock(
  file.path(data_root, ".registry.lock"),
  timeout = 60000
)
if (is.null(registry_lock)) stop("无法获取 ResearchDataHub 登记锁")
on.exit(filelock::unlock(registry_lock), add = TRUE)

catalog_before <- fread(catalog_path, colClasses = "character", na.strings = NULL)
project_before <- fread(
  project_datasets_path, colClasses = "character", na.strings = NULL
)
if (dataset_key %in% catalog_before$dataset_key) stop("CATALOG.tsv 已存在 dataset_key")
if (dataset_key %in% project_before$dataset_key) {
  stop("data/datasets.tsv 已存在 dataset_key")
}

catalog_after <- rbindlist(
  list(catalog_before, catalog_row), use.names = TRUE, fill = FALSE
)
project_after <- rbindlist(
  list(project_before, project_row), use.names = TRUE, fill = FALSE
)
temp_catalog <- tempfile(pattern = "CATALOG.tsv.", tmpdir = dirname(catalog_path))
temp_project <- tempfile(pattern = "datasets.tsv.", tmpdir = dirname(project_datasets_path))
fwrite(catalog_after, temp_catalog, sep = "\t", quote = FALSE, na = "")
fwrite(project_after, temp_project, sep = "\t", quote = FALSE, na = "")

if (!file.rename(temp_catalog, catalog_path)) stop("无法原子更新 CATALOG.tsv")
if (!file.rename(temp_project, project_datasets_path)) {
  rollback_catalog <- tempfile(
    pattern = "CATALOG.tsv.rollback.", tmpdir = dirname(catalog_path)
  )
  fwrite(catalog_before, rollback_catalog, sep = "\t", quote = FALSE, na = "")
  rollback_ok <- file.rename(rollback_catalog, catalog_path)
  if (!rollback_ok) {
    stop("项目数据引用更新失败，且 CATALOG.tsv 自动回滚失败；需人工对账")
  }
  stop("无法原子更新 data/datasets.tsv；CATALOG.tsv 已回滚")
}

message("已冻结 HRA000651 的 53 个体/106 文件行元数据，并以 incomplete/controlled_pending 登记；未下载原始 FASTQ。")
