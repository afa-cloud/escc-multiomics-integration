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
  data_root, "_incoming", "NCBI_SRA", "PRJNA766558", "retrieved_20260711"
)
canonical_root <- file.path(
  data_root, "datasets", "public", "NCBI_SRA", "PRJNA766558",
  "retrieved_20260711"
)
env_prefix <- Sys.getenv(
  "AMP_ENV_PREFIX",
  path.expand("~/.local/share/mamba/envs/escc-amplicon")
)
prefetch_bin <- file.path(env_prefix, "bin", "prefetch")
validate_bin <- file.path(env_prefix, "bin", "vdb-validate")
download_date <- "2026-07-11"
dataset_key <- "NCBI_SRA_PRJNA766558_retrieved_20260711"

stopifnot(
  dir_exists(data_root),
  file_access(catalog_path, "read"),
  dir_exists(stage_root),
  file_exists(prefetch_bin),
  file_exists(validate_bin)
)
if (dir_exists(canonical_root)) stop("规范目录已存在，拒绝覆盖：", canonical_root)

catalog <- fread(catalog_path, colClasses = "character", na.strings = NULL)
if (dataset_key %in% catalog$dataset_key) stop("CATALOG.tsv 已存在 dataset_key")

runinfo_path <- file.path(
  stage_root, "10_metadata", "ncbi_sra_runinfo_PRJNA766558_20260711.csv"
)
required_metadata <- c(
  file.path(stage_root, "10_metadata", "ncbi_sra_esearch_PRJNA766558_20260711.json"),
  runinfo_path,
  file.path(stage_root, "10_metadata", "ncbi_sra_experiment_package_PRJNA766558_20260711.xml"),
  file.path(stage_root, "10_metadata", "ncbi_bioproject_esearch_PRJNA766558_20260711.json"),
  file.path(stage_root, "10_metadata", "ncbi_bioproject_PRJNA766558_20260711.xml")
)
if (!all(file_exists(required_metadata))) {
  stop("NCBI 元数据尚未完整冻结：", paste(required_metadata[!file_exists(required_metadata)], collapse = ";"))
}

runinfo <- fread(runinfo_path, colClasses = "character", na.strings = NULL)
required_columns <- c(
  "Run", "ReleaseDate", "spots", "bases", "spots_with_mates", "avgLength",
  "size_MB", "download_path", "Experiment", "LibraryStrategy",
  "LibrarySelection", "LibrarySource", "LibraryLayout", "Platform", "Model",
  "SRAStudy", "BioProject", "Sample", "BioSample", "SampleName", "RunHash",
  "ReadHash"
)
if (!all(required_columns %in% names(runinfo))) stop("RunInfo 缺少关键字段")
if (nrow(runinfo) != 42L || uniqueN(runinfo$Run) != 42L) stop("PRJNA766558 不是 42 个唯一 run")

run_numbers <- sort(as.integer(sub("^SRR", "", runinfo$Run)))
if (!identical(run_numbers, 16095367:16095408)) {
  stop("SRR 范围不是预期的连续 SRR16095367-SRR16095408")
}
if (!all(runinfo$LibraryStrategy == "AMPLICON") ||
    !all(runinfo$LibrarySelection == "PCR") ||
    !all(runinfo$LibraryLayout == "PAIRED") ||
    !all(runinfo$Platform == "ILLUMINA") ||
    !all(runinfo$Model == "Illumina MiSeq")) {
  stop("RunInfo 的扩增子、双端或平台字段不一致")
}
if (!all(grepl("^(NS|TS)[0-9]+$", runinfo$SampleName))) {
  stop("SampleName 不能解析为 NS/TS 配对编号")
}

pairing <- runinfo[, .(
  run_accession = Run,
  experiment_accession = Experiment,
  biosample_accession = BioSample,
  sra_sample_accession = Sample,
  sample_name = SampleName,
  paper_pair_number = as.integer(sub("^(NS|TS)", "", SampleName)),
  tissue_role = fifelse(startsWith(SampleName, "TS"), "tumor", "paired_non_tumor"),
  library_strategy = LibraryStrategy,
  library_selection = LibrarySelection,
  library_source = LibrarySource,
  layout = LibraryLayout,
  platform = Platform,
  instrument_model = Model,
  spots = as.numeric(spots),
  bases = as.numeric(bases),
  spots_with_mates = as.numeric(spots_with_mates),
  avg_length = as.numeric(avgLength),
  source_size_mb = as.numeric(size_MB),
  release_date = ReleaseDate,
  run_hash = RunHash,
  read_hash = ReadHash,
  source_download_path = download_path
)]
pairing[, patient_pair_id := sprintf("FFPE_pair_%02d", paper_pair_number)]
pair_check <- pairing[, .(
  sample_count = .N,
  tumor_count = sum(tissue_role == "tumor"),
  non_tumor_count = sum(tissue_role == "paired_non_tumor")
), by = patient_pair_id]
if (nrow(pair_check) != 21L ||
    !all(pair_check$sample_count == 2L) ||
    !all(pair_check$tumor_count == 1L) ||
    !all(pair_check$non_tumor_count == 1L)) {
  stop("21 对 TS/NS 映射校验失败")
}

pairing[, `:=`(
  marker = "16S rRNA V3-V4",
  primer_forward_name = "336F",
  primer_forward_sequence = "GTACTCCTACGGGAGGCAGCA",
  primer_reverse_name = "806R",
  primer_reverse_sequence = "GTGGACTACHVGGGTWTCTAAT",
  read_configuration = "Illumina MiSeq 2x250 bp paired-end",
  ffpe_status = "FFPE tissue block",
  public_negative_control_run = "not_identified",
  pairing_status = "confirmed_by_NS_TS_sample_name",
  analysis_role = "paired_patient_amplicon_analysis"
)]
setcolorder(pairing, c(
  "patient_pair_id", "paper_pair_number", "tissue_role", "sample_name",
  "run_accession", "experiment_accession", "biosample_accession",
  "sra_sample_accession", "library_strategy", "library_selection",
  "library_source", "layout", "platform", "instrument_model", "spots",
  "bases", "spots_with_mates", "avg_length", "source_size_mb", "release_date",
  "run_hash", "read_hash", "marker", "primer_forward_name",
  "primer_forward_sequence", "primer_reverse_name", "primer_reverse_sequence",
  "read_configuration", "ffpe_status", "public_negative_control_run",
  "pairing_status", "analysis_role", "source_download_path"
))
setorder(pairing, paper_pair_number, tissue_role)

pairing_path <- file.path(stage_root, "10_metadata", "prjna766558_pairing.tsv")
fwrite(pairing, pairing_path, sep = "\t", quote = FALSE, na = "")
run_list_path <- file.path(stage_root, "10_metadata", "prjna766558_run_accessions.txt")
writeLines(sort(pairing$run_accession), run_list_path, useBytes = TRUE)

sra_root <- file.path(stage_root, "00_source", "sra")
dir_create(sra_root, recurse = TRUE)
prefetch_status <- system2(
  prefetch_bin,
  c(
    "--option-file", shQuote(run_list_path),
    "--output-directory", shQuote(sra_root),
    "--max-size", "u"
  )
)
if (!identical(prefetch_status, 0L)) stop("prefetch 下载失败；可原路径断点续传")

sra_files <- file.path(sra_root, pairing$run_accession, paste0(pairing$run_accession, ".sra"))
if (!all(file_exists(sra_files))) {
  stop("prefetch 完成后仍缺少 SRA 文件：", paste(pairing$run_accession[!file_exists(sra_files)], collapse = ";"))
}

validation_status <- vapply(sra_files, function(path) {
  output <- system2(validate_bin, shQuote(path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  is.null(status) || identical(status, 0L)
}, logical(1))
if (!all(validation_status)) {
  failed <- pairing$run_accession[!validation_status]
  stop("vdb-validate 失败：", paste(failed, collapse = ";"))
}

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}
source_inventory <- data.table(
  run_accession = pairing$run_accession,
  relative_path = path_rel(sra_files, start = stage_root),
  size_bytes = as.numeric(file_info(sra_files)$size),
  sha256 = vapply(sra_files, sha256_file, character(1)),
  vdb_validate_status = "passed",
  source_url = paste0("https://www.ncbi.nlm.nih.gov/sra/", pairing$run_accession)
)
setorder(source_inventory, run_accession)
source_inventory_path <- file.path(stage_root, "10_metadata", "sra_source_files.tsv")
fwrite(source_inventory, source_inventory_path, sep = "\t", quote = FALSE, na = "")

provenance <- data.table(
  output_scope = c(
    "10_metadata/prjna766558_pairing.tsv",
    "00_source/sra/*/*.sra",
    "10_metadata/sra_source_files.tsv"
  ),
  corresponding_source_file = c(
    "10_metadata/ncbi_sra_runinfo_PRJNA766558_20260711.csv",
    "NCBI SRA via prefetch",
    "00_source/sra/*/*.sra"
  ),
  generating_script = "scripts/08_download_ingest_prjna766558_sra.R",
  key_parameters = c(
    "42 unique runs; NS/TS suffix pairing; 21 pairs; exact primer sequences from source publication",
    "prefetch --max-size u; preserve full SRA quality values; no lite files",
    "vdb-validate each run; local SHA256"
  ),
  software = paste0(
    "R ", getRversion(), "; data.table ", packageVersion("data.table"),
    "; sra-tools from ", env_prefix
  ),
  generated_date = download_date,
  regenerable = "yes"
)
provenance_path <- file.path(stage_root, "20_reusable", "PROVENANCE.tsv")
fwrite(provenance, provenance_path, sep = "\t", quote = FALSE, na = "")

dataset_md <- c(
  "# PRJNA766558 ESCC 组织 16S 数据说明",
  "",
  "## 基本信息",
  "",
  paste0("- `dataset_key`：`", dataset_key, "`。"),
  "- 来源：NCBI BioProject/SRA；BioProject `PRJNA766558`，SRA Study `SRP339053`。",
  "- 本地冻结版本：`retrieved_20260711`；获取日期：2026-07-11。",
  "- 物种：人（Homo sapiens）；疾病：食管鳞状细胞癌（ESCC）。",
  "- 样本：21 对 FFPE 肿瘤/配对癌旁组织，共 42 个 paired-end 16S run。",
  "- 标记区：16S rRNA V3-V4；336F `GTACTCCTACGGGAGGCAGCA`；806R `GTGGACTACHVGGGTWTCTAAT`。",
  "- 测序：Illumina MiSeq 2×250 bp。",
  "- 当前完整性状态：`verified`；42 个完整 `.sra` 均通过 `vdb-validate` 与本地 SHA256。",
  "",
  "## 目录与主要文件",
  "",
  "- `00_source/sra/<SRR>/<SRR>.sra`：通过 SRA Toolkit `prefetch` 取得的完整 SRA，不使用 `.lite`。",
  "- `10_metadata/ncbi_sra_*`：NCBI Entrez ESearch、RunInfo 和 Experiment Package 原始快照。",
  "- `10_metadata/ncbi_bioproject_*`：NCBI BioProject 查询与 XML 快照。",
  "- `10_metadata/prjna766558_pairing.tsv`：21 对 TS/NS、BioSample、run、引物和 FFPE 角色映射。",
  "- `10_metadata/sra_source_files.tsv`：逐 run 大小、SHA256 和 `vdb-validate` 状态。",
  "- `20_reusable/PROVENANCE.tsv`：元数据映射和下载生成关系。",
  "- `90_manifests/MANIFEST.tsv`：逐文件大小、SHA256、来源和完整性状态。",
  "",
  "## 推荐使用场景",
  "",
  "- 以患者为推断单位的 21 对肿瘤/癌旁 ASV、α/β 多样性和差异丰度分析。",
  "- 与宿主多组学在属、科、功能潜力或通路模块层做独立生态校准。",
  "",
  "## 不适用场景与限制",
  "",
  "- FFPE、低生物量组织和口腔来源污染风险需要单独敏感性分析。",
  "- 当前公开 RunInfo 中未识别到独立阴性对照或提取空白 run；该缺失作为软性降分，不直接硬拒绝队列。",
  "- 公开 16S 仅覆盖论文 41 对组织中的 21 对；WES 需向作者申请，不得写成已公开取得。",
  "- 在查看 read quality 前不得预设 DADA2 截断长度。",
  "- 与主基因组/转录组队列不是同一患者，不做虚假样本级相关。",
  "- 微生物功能预测只表示潜力，不等同于实测宿主代谢或通路活动。",
  "",
  "## 来源、引用与许可",
  "",
  "- BioProject：https://www.ncbi.nlm.nih.gov/bioproject/PRJNA766558",
  "- 原始论文：https://link.springer.com/article/10.1186/s12866-021-02352-6",
  "- 使用时应引用 BioProject、SRA 和原始论文，并遵守来源条款。",
  paste0("- 使用项目：`", project_root, "`。")
)
writeLines(dataset_md, file.path(stage_root, "DATASET.md"), useBytes = TRUE)

manifest_files <- c(
  required_metadata,
  pairing_path,
  run_list_path,
  source_inventory_path,
  sra_files,
  provenance_path
)
relative_paths <- path_rel(manifest_files, start = stage_root)
manifest <- data.table(
  relative_path = relative_paths,
  file_level = sub("/.*$", "", relative_paths),
  size_bytes = as.numeric(file_info(manifest_files)$size),
  sha256 = vapply(manifest_files, sha256_file, character(1)),
  source_url = "generated from NCBI SRA metadata",
  download_date = download_date,
  file_status = "generated_verified",
  corresponding_source_file = "10_metadata/ncbi_sra_runinfo_PRJNA766558_20260711.csv",
  generation_method = "scripts/08_download_ingest_prjna766558_sra.R",
  notes = ""
)
manifest[relative_path %in% path_rel(required_metadata, start = stage_root), `:=`(
  source_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
  file_status = "verified",
  corresponding_source_file = "",
  generation_method = "NCBI Entrez skill raw response",
  notes = "official NCBI metadata snapshot"
)]
for (i in seq_len(nrow(source_inventory))) {
  relative <- source_inventory$relative_path[i]
  manifest[relative_path == relative, `:=`(
    source_url = source_inventory$source_url[i],
    file_status = "verified",
    corresponding_source_file = "",
    generation_method = "SRA Toolkit prefetch; source bytes retained",
    notes = "full_sra;vdb_validate_passed"
  )]
}
setorder(manifest, relative_path)
manifest_path <- file.path(stage_root, "90_manifests", "MANIFEST.tsv")
fwrite(manifest, manifest_path, sep = "\t", quote = FALSE, na = "")

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

catalog_row <- data.table(
  dataset_key = dataset_key,
  record_type = "paired_tissue_amplicon_dataset",
  access_level = "public",
  source = "NCBI SRA",
  accession = "PRJNA766558;SRP339053;SRR16095367-SRR16095408",
  version = "retrieved_20260711",
  species = "Homo sapiens",
  disease = "esophageal squamous cell carcinoma",
  tissue = "FFPE primary tumor; paired adjacent non-tumor esophagus",
  assay = "16S rRNA V3-V4 amplicon; Illumina MiSeq paired-end 250 bp",
  sample_summary = "42 complete SRA runs; 21 confirmed TS/NS patient pairs",
  available_levels = "00_source;10_metadata",
  status = "verified",
  local_path = canonical_root,
  manifest_path = file.path(canonical_root, "90_manifests", "MANIFEST.tsv"),
  source_url = "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA766558",
  download_date = download_date,
  last_verified = download_date,
  license_or_access = "public; cite NCBI SRA and original study; follow source terms",
  projects_using = project_root,
  recommended_use = "paired ESCC tumor/non-tumor 16S ASV and independent microbial ecology calibration",
  limitations = "FFPE low-biomass contamination risk; no public negative-control run identified; 21-pair public subset only; WES unavailable",
  notes = "Do not use RunInfo .lite download paths for DADA2; full prefetch SRA files are retained."
)
# 下载与全量 SHA 可能持续较久；登记前重读目录，避免用启动时快照覆盖其他串行入库结果。
catalog <- fread(catalog_path, colClasses = "character", na.strings = NULL)
if (dataset_key %in% catalog$dataset_key) stop("CATALOG.tsv 已存在 dataset_key")
catalog <- rbindlist(list(catalog, catalog_row), use.names = TRUE, fill = FALSE)
tmp_catalog <- tempfile(pattern = "CATALOG.tsv.", tmpdir = dirname(catalog_path))
fwrite(catalog, tmp_catalog, sep = "\t", quote = FALSE, na = "")
if (!file.rename(tmp_catalog, catalog_path)) stop("无法原子更新 CATALOG.tsv")

project_datasets_path <- file.path(project_root, "data", "datasets.tsv")
project_datasets <- fread(project_datasets_path, colClasses = "character", na.strings = NULL)
if (dataset_key %in% project_datasets$dataset_key) stop("data/datasets.tsv 已存在 dataset_key")
project_row <- data.table(
  logical_name = "escc_tissue_microbiome_prjna766558",
  dataset_key = dataset_key,
  accession = "PRJNA766558;SRP339053;SRR16095367-SRR16095408",
  version = "retrieved_20260711",
  data_level = "00_source;10_metadata",
  central_path = canonical_root,
  project_purpose = "21 对 ESCC 肿瘤/癌旁 16S 生态差异及宿主通路模块校准",
  inclusion_exclusion = "保留患者配对；查看质量分布后再定截断；无阴性对照时执行污染敏感性分析",
  limitations = "FFPE/低生物量；无公开阴性对照 run；仅 21 对公开 16S；与宿主主队列无同患者映射",
  local_link = "",
  last_verified = download_date
)
project_datasets <- rbindlist(
  list(project_datasets, project_row),
  use.names = TRUE,
  fill = FALSE
)
tmp_project <- tempfile(pattern = "datasets.tsv.", tmpdir = dirname(project_datasets_path))
fwrite(project_datasets, tmp_project, sep = "\t", quote = FALSE, na = "")
if (!file.rename(tmp_project, project_datasets_path)) stop("无法原子更新 data/datasets.tsv")

message("已完成 PRJNA766558 42 个完整 SRA 下载、21 对映射、vdb-validate、全 SHA 校验和索引登记。")
