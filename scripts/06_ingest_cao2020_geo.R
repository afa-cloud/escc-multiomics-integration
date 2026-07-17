#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(fs)
  library(GEOquery)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_root <- Sys.getenv("RESEARCH_DATA_ROOT", path.expand("~/ResearchDataHub"))
catalog_path <- file.path(data_root, "CATALOG.tsv")
download_date <- "2026-07-11"
version <- "retrieved_20260711_author_processed"

stopifnot(
  dir_exists(data_root),
  file_access(catalog_path, "read")
)

cross_map_path <- file.path(project_root, "results", "cao2020_cross_omics_sample_map.tsv")
stopifnot(file_exists(cross_map_path))
cross_map <- fread(cross_map_path, colClasses = "character", na.strings = "")

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

safe_tar_members <- function(tar_path, expected_count, expected_pattern) {
  members <- system2("/usr/bin/tar", c("-tf", shQuote(tar_path)), stdout = TRUE)
  if (!length(members)) stop("tar 无成员：", tar_path)

  unsafe <- vapply(strsplit(members, "/", fixed = TRUE), function(parts) {
    any(parts == "..")
  }, logical(1)) | startsWith(members, "/")
  if (any(unsafe)) stop("tar 含不安全路径：", paste(members[unsafe], collapse = ";"))
  if (anyDuplicated(members)) stop("tar 含重复成员名：", tar_path)

  verbose <- system2("/usr/bin/tar", c("-tvf", shQuote(tar_path)), stdout = TRUE)
  if (any(substr(verbose, 1L, 1L) %in% c("l", "h"))) {
    stop("tar 含符号链接或硬链接：", tar_path)
  }
  if (length(members) != expected_count) {
    stop("tar 成员数不符：", tar_path, " expected=", expected_count,
         " observed=", length(members))
  }
  if (!all(grepl(expected_pattern, basename(members)))) {
    stop("tar 成员扩展名或命名不符：", tar_path)
  }
  members
}

meta_value <- function(gsm_object, key) {
  value <- GEOquery::Meta(gsm_object)[[key]]
  if (is.null(value) || !length(value)) return(NA_character_)
  paste(as.character(value), collapse = "; ")
}

extract_accession <- function(text, pattern) {
  if (is.na(text) || !nzchar(text)) return(NA_character_)
  hits <- unique(unlist(regmatches(text, gregexpr(pattern, text, perl = TRUE))))
  if (!length(hits)) NA_character_ else paste(hits, collapse = ";")
}

geo_metadata_table <- function(soft_path) {
  gse <- GEOquery::getGEO(
    filename = soft_path,
    GSEMatrix = FALSE,
    AnnotGPL = FALSE
  )
  gsm_list <- GEOquery::GSMList(gse)
  rbindlist(lapply(names(gsm_list), function(gsm) {
    obj <- gsm_list[[gsm]]
    title <- meta_value(obj, "title")
    source <- meta_value(obj, "source_name_ch1")
    characteristics <- meta_value(obj, "characteristics_ch1")
    relation <- meta_value(obj, "relation")
    combined <- paste(title, source, characteristics, relation, sep = "; ")
    data.table(
      gsm = gsm,
      sample_title = title,
      source_name = source,
      characteristics = characteristics,
      biosample = extract_accession(combined, "SAMN[0-9]+"),
      sra_accession = extract_accession(combined, "SRR[0-9]+")
    )
  }), fill = TRUE)
}

make_sample_map <- function(accession, members, soft_path) {
  samples <- data.table(
    gsm = sub("^(GSM[0-9]+).*$", "\\1", basename(members)),
    filename = file.path("20_reusable", "author_processed", basename(members))
  )
  samples <- merge(samples, geo_metadata_table(soft_path), by = "gsm", all.x = TRUE)

  if (accession == "GSE149608") {
    mapping <- cross_map[, .(
      gsm = wgbseq_gsm,
      paper_patient_id,
      condition,
      mapping_note = notes
    )]
    samples <- merge(samples, mapping, by = "gsm", all.x = TRUE)
    if (anyNA(samples$paper_patient_id)) stop("GSE149608 存在未映射 GSM")
    samples[, `:=`(
      sample_role = fifelse(
        paper_patient_id == "P15" & condition == "T",
        "unpaired_tumor_sensitivity_only",
        fifelse(condition == "T", "paired_primary_tumor", "paired_adjacent_normal")
      ),
      analysis_eligibility = fifelse(
        paper_patient_id == "P15" & condition == "T",
        "unpaired_sensitivity_only",
        "paired_patient_analysis"
      ),
      mapping_evidence = "GEO SOFT + Cao 2020 published sample mapping + project cross-omics audit",
      mapping_confidence = "high"
    )]
  } else if (accession == "GSE149609") {
    mapping <- cross_map[, .(
      gsm = rna_gsm,
      paper_patient_id,
      condition,
      mapping_note = notes
    )]
    samples <- merge(samples, mapping, by = "gsm", all.x = TRUE)
    samples[, `:=`(
      sample_role = fifelse(
        !is.na(paper_patient_id) & condition == "T",
        "paired_primary_tumor",
        fifelse(!is.na(paper_patient_id), "paired_adjacent_normal", "cell_line_perturbation")
      ),
      analysis_eligibility = fifelse(
        !is.na(paper_patient_id),
        "paired_patient_analysis",
        "cell_line_mechanism_calibration"
      ),
      mapping_evidence = fifelse(
        !is.na(paper_patient_id),
        "GEO SOFT + Cao 2020 published sample mapping + project cross-omics audit",
        "GEO SOFT metadata and author filename"
      ),
      mapping_confidence = fifelse(!is.na(paper_patient_id), "high", "medium")
    )]
    if (sum(!is.na(samples$paper_patient_id)) != 20L) {
      stop("GSE149609 患者组织文件数不是 20")
    }
  } else if (accession == "GSE151838") {
    samples[, `:=`(
      paper_patient_id = NA_character_,
      condition = NA_character_,
      mapping_note = "Het-1A/EC109 EZH2 ChIP-seq；仅细胞系机制校准",
      sample_role = "cell_line_chipseq",
      analysis_eligibility = "cell_line_mechanism_calibration",
      mapping_evidence = "GEO SOFT metadata and author filename",
      mapping_confidence = "high"
    )]
  }

  samples[, genome_build := "hg19"]
  setcolorder(samples, c(
    "gsm", "paper_patient_id", "condition", "sample_role", "sample_title",
    "source_name", "characteristics", "biosample", "sra_accession", "filename",
    "genome_build", "analysis_eligibility", "mapping_evidence",
    "mapping_confidence", "mapping_note"
  ))
  setorder(samples, gsm)
  samples
}

dataset_definitions <- list(
  GSE149608 = list(
    accession = "GSE149608",
    expected_count = 19L,
    expected_pattern = "\\.bed\\.gz$",
    assay = "whole-genome bisulfite sequencing; author-processed CpG methylation",
    tissue = "primary tumor; adjacent normal esophagus",
    summary = "19 个作者处理 WGBS 样本；9 个完整肿瘤-癌旁配对，另有未配对 T15",
    use = "同患者甲基化—转录跨层锚点；严格配对分析使用 9 对",
    limitations = "N15 原始测序约 1.58x 且无作者处理 BED；T15 只可用于非配对敏感性分析；hg19/MOABS 处理值",
    notes = "verified 表示来源快照和生成文件完整，不表示 19 个样本均可进入严格配对分析。"
  ),
  GSE149609 = list(
    accession = "GSE149609",
    expected_count = 35L,
    expected_pattern = "\\.txt\\.gz$",
    assay = "RNA-seq; author-processed RSEM TPM",
    tissue = "primary tumor; adjacent normal esophagus; ESCC-related cell lines",
    summary = "35 个作者处理 TPM 文件；20 个患者组织和 15 个细胞系扰动样本",
    use = "10 对患者 RNA 配对层及细胞系机制校准",
    limitations = "TPM 不是原始整数计数，不能直接输入 DESeq2；15 个细胞系文件不得计入患者多组学交集",
    notes = "患者组织与细胞系已在 sample_file_map.tsv 中明确拆分。"
  ),
  GSE151838 = list(
    accession = "GSE151838",
    expected_count = 2L,
    expected_pattern = "\\.bed\\.gz$",
    assay = "EZH2 ChIP-seq; author-processed MACS2 peaks",
    tissue = "Het-1A and EC109 cell lines",
    summary = "2 个细胞系 EZH2 ChIP-seq peak 文件",
    use = "EZH2 结合位点的细胞系机制校准",
    limitations = "无患者组织，不能充当患者级表观遗传证据或独立队列",
    notes = "该记录仅进入 cell_line_mechanism_calibration 角色。"
  )
)

catalog_rows <- list()
project_rows <- list()

for (definition in dataset_definitions) {
  accession <- definition$accession
  stage_root <- file.path(
    data_root, "_incoming", "GEO", accession,
    "retrieved_20260711_author_processed"
  )
  canonical_root <- file.path(
    data_root, "datasets", "public", "GEO", accession,
    "retrieved_20260711_author_processed"
  )
  tar_name <- paste0(accession, "_RAW.tar")
  soft_name <- paste0(accession, "_family.soft.gz")
  tar_path <- file.path(stage_root, "00_source", tar_name)
  soft_path <- file.path(stage_root, "10_metadata", soft_name)
  filelist_path <- file.path(stage_root, "10_metadata", "filelist.txt")

  required <- c(tar_path, soft_path, filelist_path)
  if (!all(file_exists(required))) {
    stop("下载未完成：", paste(required[!file_exists(required)], collapse = ";"))
  }
  partials <- dir_ls(stage_root, recurse = TRUE, type = "file", regexp = "\\.(aria2|part)$")
  if (length(partials)) stop("暂存目录仍有未完成下载：", paste(partials, collapse = ";"))
  if (dir_exists(canonical_root)) stop("规范目录已存在，拒绝覆盖：", canonical_root)

  members <- safe_tar_members(
    tar_path,
    definition$expected_count,
    definition$expected_pattern
  )
  filelist_text <- readLines(filelist_path, warn = FALSE, encoding = "UTF-8")
  present_in_filelist <- vapply(
    basename(members),
    function(name) any(grepl(name, filelist_text, fixed = TRUE)),
    logical(1)
  )
  if (!all(present_in_filelist)) {
    stop(accession, " tar 成员未全部出现在 GEO filelist.txt")
  }

  extract_dir <- file.path(stage_root, "20_reusable", "author_processed")
  if (dir_exists(extract_dir)) dir_delete(extract_dir)
  dir_create(extract_dir, recurse = TRUE)
  utils::untar(tar_path, exdir = extract_dir)

  extracted <- file.path(extract_dir, basename(members))
  if (!all(file_exists(extracted))) stop(accession, " 解包成员缺失")
  gzip_status <- vapply(extracted, function(path) {
    identical(system2("/usr/bin/gzip", c("-t", shQuote(path))), 0L)
  }, logical(1))
  if (!all(gzip_status)) stop(accession, " 内部 gzip 成员校验失败")

  tar_members <- data.table(
    member_name = members,
    member_size_bytes = as.numeric(file_info(extracted)$size),
    extracted_relative_path = file.path(
      "20_reusable", "author_processed", basename(extracted)
    ),
    sha256 = vapply(extracted, sha256_file, character(1)),
    file_type = sub("^.*\\.", "", basename(extracted)),
    integrity_status = "gzip_verified"
  )
  fwrite(
    tar_members,
    file.path(stage_root, "10_metadata", "tar_members.tsv"),
    sep = "\t", quote = FALSE, na = ""
  )

  sample_map <- make_sample_map(accession, members, soft_path)
  fwrite(
    sample_map,
    file.path(stage_root, "10_metadata", "sample_file_map.tsv"),
    sep = "\t", quote = FALSE, na = ""
  )

  provenance <- data.table(
    output_scope = "20_reusable/author_processed/*",
    corresponding_source_file = file.path("00_source", tar_name),
    input_sha256 = sha256_file(tar_path),
    generating_script = "scripts/06_ingest_cao2020_geo.R",
    key_parameters = paste0(
      "safe member audit; expected_count=", definition$expected_count,
      "; preserve author filenames; gzip -t each member"
    ),
    software = paste0(
      "R ", getRversion(), "; GEOquery ", packageVersion("GEOquery"),
      "; data.table ", packageVersion("data.table")
    ),
    generated_date = download_date,
    regenerable = "yes"
  )
  fwrite(
    provenance,
    file.path(stage_root, "20_reusable", "PROVENANCE.tsv"),
    sep = "\t", quote = FALSE, na = ""
  )

  dataset_key <- paste0("GEO_", accession, "_", version)
  dataset_md <- c(
    paste0("# ", accession, " 作者处理数据说明"),
    "",
    "## 基本信息",
    "",
    paste0("- `dataset_key`：`", dataset_key, "`。"),
    paste0("- 来源：NCBI Gene Expression Omnibus（GEO），accession `", accession, "`。"),
    paste0("- 本地冻结版本：`", version, "`；获取日期：", download_date, "。"),
    "- 物种：人（Homo sapiens）。",
    "- 疾病：食管鳞状细胞癌（ESCC）。",
    paste0("- 检测类型：", definition$assay, "。"),
    paste0("- 样本摘要：", definition$summary, "。"),
    "- 当前完整性状态：`verified`；该状态表示来源字节与生成文件完整，不等于全部样本均可进入患者级整合。",
    "",
    "## 目录与主要文件",
    "",
    paste0("- `00_source/", tar_name, "`：GEO 作者处理数据原始 tar，保留原始字节。"),
    paste0("- `10_metadata/", soft_name, "`：GEO family SOFT 元数据。"),
    "- `10_metadata/filelist.txt`：GEO 官方补充文件清单。",
    "- `10_metadata/tar_members.tsv`：tar 成员、大小、SHA256 和逐成员 gzip 校验。",
    "- `10_metadata/sample_file_map.tsv`：GSM、患者/细胞系角色、组织条件和分析资格。",
    "- `20_reusable/author_processed/`：从原 tar 安全解包、文件名和内容未转换的作者处理文件。",
    "- `20_reusable/PROVENANCE.tsv`：解包与映射生成关系。",
    "- `90_manifests/MANIFEST.tsv`：逐文件大小、SHA256、来源和完整性状态。",
    "",
    "## 推荐使用场景",
    "",
    paste0("- ", definition$use, "。"),
    "",
    "## 不适用场景与限制",
    "",
    paste0("- ", definition$limitations, "。"),
    "- 公共数据关联只用于假设生成，不证明因果机制或治疗靶点。",
    "",
    "## 来源、生成关系与复用",
    "",
    paste0("- GEO 页面：https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", accession),
    "- `20_reusable/author_processed/` 可由 `00_source/` 原 tar 重新生成；未修改成员内容。",
    "- 项目脚本通过 `RESEARCH_DATA_ROOT` 与 `data/datasets.tsv` 定位本记录。",
    paste0("- 使用项目：`", project_root, "`。"),
    "",
    "## 引用、许可与完整性",
    "",
    "- 使用时应引用 GEO 记录和 Cao 等 2020 原始论文，并遵守来源条款。",
    "- SHA256 清单及逐文件状态见 `90_manifests/MANIFEST.tsv`。",
    paste0("- 特别说明：", definition$notes)
  )
  writeLines(dataset_md, file.path(stage_root, "DATASET.md"), useBytes = TRUE)

  tar_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    substr(accession, 1L, 6L), "nnn/", accession,
    "/suppl/", tar_name
  )
  soft_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    substr(accession, 1L, 6L), "nnn/", accession,
    "/soft/", soft_name
  )
  filelist_url <- paste0(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/",
    substr(accession, 1L, 6L), "nnn/", accession,
    "/suppl/filelist.txt"
  )

  manifest_files <- c(
    tar_path,
    soft_path,
    filelist_path,
    file.path(stage_root, "10_metadata", "tar_members.tsv"),
    file.path(stage_root, "10_metadata", "sample_file_map.tsv"),
    extracted,
    file.path(stage_root, "20_reusable", "PROVENANCE.tsv")
  )
  relative_paths <- path_rel(manifest_files, start = stage_root)
  manifest <- data.table(
    relative_path = relative_paths,
    file_level = sub("/.*$", "", relative_paths),
    size_bytes = as.numeric(file_info(manifest_files)$size),
    sha256 = vapply(manifest_files, sha256_file, character(1)),
    source_url = "generated from GEO source files",
    download_date = download_date,
    file_status = "generated_verified",
    corresponding_source_file = file.path("00_source", tar_name),
    generation_method = "scripts/06_ingest_cao2020_geo.R",
    notes = ""
  )
  manifest[relative_path == file.path("00_source", tar_name), `:=`(
    source_url = tar_url,
    file_status = "verified",
    corresponding_source_file = "",
    generation_method = "source bytes retrieved without modification",
    notes = "tar container and all internal gzip members verified"
  )]
  manifest[relative_path == file.path("10_metadata", soft_name), `:=`(
    source_url = soft_url,
    file_status = "verified",
    corresponding_source_file = "",
    generation_method = "source bytes retrieved without modification",
    notes = "gzip stream verified"
  )]
  manifest[relative_path == "10_metadata/filelist.txt", `:=`(
    source_url = filelist_url,
    file_status = "verified",
    corresponding_source_file = "",
    generation_method = "source bytes retrieved without modification",
    notes = "GEO official supplementary file list"
  )]
  manifest[startsWith(relative_path, "20_reusable/author_processed/"), `:=`(
    source_url = tar_url,
    corresponding_source_file = file.path("00_source", tar_name),
    generation_method = "safe extraction from GEO RAW tar; content unchanged",
    notes = "gzip_verified"
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
  if (!all(file_exists(canonical_files))) stop(accession, " 规范路径文件缺失")
  observed_sizes <- as.character(as.numeric(file_info(canonical_files)$size))
  observed_sha <- vapply(canonical_files, sha256_file, character(1))
  if (!all(observed_sizes == canonical_manifest$size_bytes) ||
      !all(observed_sha == canonical_manifest$sha256)) {
    stop(accession, " 规范路径复制后校验失败")
  }

  record_type <- switch(
    accession,
    GSE149608 = "paired_epigenome_dataset",
    GSE149609 = "paired_transcriptome_dataset",
    GSE151838 = "cell_line_regulatory_dataset"
  )
  catalog_rows[[accession]] <- data.table(
    dataset_key = dataset_key,
    record_type = record_type,
    access_level = "public",
    source = "GEO",
    accession = accession,
    version = version,
    species = "Homo sapiens",
    disease = "esophageal squamous cell carcinoma",
    tissue = definition$tissue,
    assay = definition$assay,
    sample_summary = definition$summary,
    available_levels = "00_source;10_metadata;20_reusable",
    status = "verified",
    local_path = canonical_root,
    manifest_path = file.path(canonical_root, "90_manifests", "MANIFEST.tsv"),
    source_url = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", accession),
    download_date = download_date,
    last_verified = download_date,
    license_or_access = "public; cite GEO record and Cao et al. 2020; follow source terms",
    projects_using = project_root,
    recommended_use = definition$use,
    limitations = definition$limitations,
    notes = definition$notes
  )

  logical_name <- switch(
    accession,
    GSE149608 = "cao2020_wgbs_author_processed",
    GSE149609 = "cao2020_rna_tpm_author_processed",
    GSE151838 = "cao2020_ezh2_chipseq_cell_line"
  )
  inclusion <- switch(
    accession,
    GSE149608 = "严格患者级配对分析仅纳入 9 对；T15 只进入非配对敏感性分析",
    GSE149609 = "患者组织 10 对进入配对层；15 个细胞系样本单列为机制校准",
    GSE151838 = "仅作 Het-1A/EC109 细胞系机制校准，不进入患者多组学交集"
  )
  project_rows[[accession]] <- data.table(
    logical_name = logical_name,
    dataset_key = dataset_key,
    accession = accession,
    version = version,
    data_level = "00_source;10_metadata;20_reusable",
    central_path = canonical_root,
    project_purpose = definition$use,
    inclusion_exclusion = inclusion,
    limitations = definition$limitations,
    local_link = "",
    last_verified = download_date
  )
}

catalog <- fread(catalog_path, colClasses = "character", na.strings = NULL)
new_catalog <- rbindlist(catalog_rows, use.names = TRUE, fill = FALSE)
if (any(new_catalog$dataset_key %in% catalog$dataset_key)) {
  stop("CATALOG.tsv 已存在待写入 dataset_key，拒绝重复登记")
}
catalog <- rbindlist(list(catalog, new_catalog), use.names = TRUE, fill = FALSE)
write_tsv_atomic(catalog, catalog_path)

project_datasets_path <- file.path(project_root, "data", "datasets.tsv")
project_datasets <- fread(
  project_datasets_path,
  colClasses = "character",
  na.strings = NULL
)
new_project_rows <- rbindlist(project_rows, use.names = TRUE, fill = FALSE)
if (any(new_project_rows$dataset_key %in% project_datasets$dataset_key)) {
  stop("data/datasets.tsv 已存在待写入 dataset_key，拒绝重复登记")
}
project_datasets <- rbindlist(
  list(project_datasets, new_project_rows),
  use.names = TRUE,
  fill = FALSE
)
write_tsv_atomic(project_datasets, project_datasets_path)

message("已完成 Cao 2020 三个 GEO 作者处理数据集的安全解包、规范路径复制、全 SHA 校验和索引登记。")
