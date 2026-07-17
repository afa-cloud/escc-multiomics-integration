#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(fs)
  library(jsonlite)
  library(parallel)
})

options(stringsAsFactors = FALSE, scipen = 999)

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stop("无法唯一定位当前脚本。")
script_path <- normalizePath(
  sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE
)
project_root <- dirname(dirname(script_path))
results_dir <- file.path(project_root, "results")
checks_dir <- file.path(project_root, "_work", "checks")
if (!file.exists(file.path(project_root, "PROJECT_INDEX.md")) ||
    !dir.exists(results_dir) || !dir.exists(checks_dir)) {
  stop("项目根目录或规范输出目录缺失。")
}

data_root <- Sys.getenv("RESEARCH_DATA_ROOT", path.expand("~/ResearchDataHub"))
target_dataset_key <- "NCBI_SRA_PRJNA766558_retrieved_20260711"
canonical_root <- file.path(
  data_root, "datasets", "public", "NCBI_SRA", "PRJNA766558",
  "retrieved_20260711"
)
catalog_path <- file.path(data_root, "CATALOG.tsv")
manifest_path <- file.path(canonical_root, "90_manifests", "MANIFEST.tsv")
pairing_path <- file.path(canonical_root, "10_metadata", "prjna766558_pairing.tsv")
provenance_path <- file.path(canonical_root, "20_reusable", "PROVENANCE.tsv")
fastq_qc_path <- file.path(canonical_root, "10_metadata", "fastq_conversion_qc.tsv")
final_fastq_dir <- file.path(canonical_root, "20_reusable", "fastq")
incoming_dir <- file.path(
  data_root, "_incoming", "PRJNA766558_fastq_retrieved_20260711"
)

required_paths <- c(
  canonical_root, catalog_path, manifest_path, pairing_path, provenance_path
)
if (!all(file.exists(required_paths) | dir.exists(required_paths))) {
  stop("PRJNA766558 规范输入、CATALOG 或 manifest 缺失。")
}

amplicon_env <- Sys.getenv(
  "ESCC_AMPLICON_ENV",
  path.expand("~/.local/share/mamba/envs/escc-amplicon")
)
tool_path <- function(name) file.path(amplicon_env, "bin", name)
tools <- c(
  fasterq_dump = tool_path("fasterq-dump"),
  pigz = tool_path("pigz"),
  seqkit = tool_path("seqkit"),
  cutadapt = tool_path("cutadapt")
)
if (any(!file.exists(tools))) {
  stop("扩增子环境工具缺失：", paste(names(tools)[!file.exists(tools)], collapse = ", "))
}

fail_if <- function(condition, message) {
  if (isTRUE(condition)) stop(message, call. = FALSE)
}

atomic_fwrite <- function(object, path) {
  dir_create(dirname(path), recurse = TRUE)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  fwrite(object, temporary, sep = "\t", quote = FALSE, na = "")
  fail_if(!file.rename(temporary, path), paste0("无法原子发布：", path))
}

catalog <- fread(catalog_path, colClasses = "character", na.strings = NULL)
catalog_row <- catalog[catalog[["dataset_key"]] == target_dataset_key]
fail_if(
  nrow(catalog_row) != 1L || catalog_row$status != "verified" ||
    normalizePath(catalog_row$local_path, winslash = "/", mustWork = TRUE) !=
      normalizePath(canonical_root, winslash = "/", mustWork = TRUE),
  "PRJNA766558 CATALOG 记录不唯一、未验证或路径不一致"
)

pairing <- fread(pairing_path, colClasses = "character", na.strings = NULL)
pairing[, spots := as.integer(spots)]
fail_if(
  nrow(pairing) != 42L || uniqueN(pairing$run_accession) != 42L ||
    uniqueN(pairing$patient_pair_id) != 21L ||
    any(pairing[, .N, by = patient_pair_id]$N != 2L),
  "必须为 42 个唯一 run 和 21 个完整患者对"
)
fail_if(
  !all(pairing$layout == "PAIRED") ||
    !all(pairing$analysis_role == "paired_patient_amplicon_analysis"),
  "存在非 paired-end 或非正式配对分析 run"
)

manifest <- fread(manifest_path, colClasses = "character", na.strings = NULL)
sra_relative <- file.path(
  "00_source", "sra", pairing$run_accession,
  paste0(pairing$run_accession, ".sra")
)
sra_manifest <- manifest[match(sra_relative, relative_path)]
fail_if(
  nrow(sra_manifest) != 42L || anyNA(sra_manifest$relative_path) ||
    !all(sra_manifest$file_status == "verified"),
  "42 个 SRA 未完整登记为 verified"
)
sra_paths <- file.path(canonical_root, sra_relative)
fail_if(any(!file.exists(sra_paths)), "已登记 SRA 文件缺失")
sra_size_ok <- as.character(as.numeric(file.info(sra_paths)$size)) ==
  sra_manifest$size_bytes
sra_sha <- vapply(
  sra_paths, digest, character(1), algo = "sha256", file = TRUE,
  serialize = FALSE
)
fail_if(
  !all(sra_size_ok) || !identical(unname(sra_sha), sra_manifest$sha256),
  "SRA 大小或 SHA256 与 manifest 不一致"
)

dir_create(incoming_dir, recurse = TRUE)
dir_create(final_fastq_dir, recurse = TRUE)

# SRA Toolkit 3.4.1 在包含中文字符的当前工作目录下可触发
# rcConverting/rcBuffer/rcInsufficient。所有 fasterq 子进程固定从 /tmp 启动，
# 输入输出仍使用绝对路径。
old_workdir <- getwd()
setwd("/tmp")
on.exit(setwd(old_workdir), add = TRUE)

convert_one <- function(index) {
  run <- pairing$run_accession[[index]]
  sra <- sra_paths[[index]]
  stage_r1 <- file.path(incoming_dir, paste0(run, "_1.fastq.gz"))
  stage_r2 <- file.path(incoming_dir, paste0(run, "_2.fastq.gz"))
  if (file.exists(stage_r1) && file.exists(stage_r2)) {
    gzip_ok <- system2(tools[["pigz"]], c("-t", stage_r1, stage_r2),
                       stdout = FALSE, stderr = FALSE) == 0L
    if (gzip_ok) return(data.table(run_accession = run, conversion = "reused_incoming"))
  }

  temporary <- file.path("/tmp", paste0("prjna766558_fasterq_", run))
  if (dir.exists(temporary)) dir_delete(temporary)
  dir_create(temporary, recurse = TRUE)
  on.exit(if (dir.exists(temporary)) dir_delete(temporary), add = TRUE)

  log <- system2(
    tools[["fasterq_dump"]],
    c(
      "--split-files", "--threads", "2", "--temp", temporary,
      "--outdir", temporary, sra
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  exit_status <- attr(log, "status")
  if (is.null(exit_status)) exit_status <- 0L
  fail_if(exit_status != 0L, paste(run, "fasterq-dump 失败：", paste(log, collapse = " | ")))

  raw_r1 <- file.path(temporary, paste0(run, "_1.fastq"))
  raw_r2 <- file.path(temporary, paste0(run, "_2.fastq"))
  fail_if(!file.exists(raw_r1) || !file.exists(raw_r2), paste(run, "缺少成对 FASTQ"))
  pigz_status <- system2(
    tools[["pigz"]], c("-p", "2", "-f", raw_r1, raw_r2),
    stdout = FALSE, stderr = FALSE
  )
  fail_if(pigz_status != 0L, paste(run, "pigz 压缩失败"))
  file_copy(paste0(raw_r1, ".gz"), stage_r1, overwrite = TRUE)
  file_copy(paste0(raw_r2, ".gz"), stage_r2, overwrite = TRUE)
  fail_if(
    system2(tools[["pigz"]], c("-t", stage_r1, stage_r2),
            stdout = FALSE, stderr = FALSE) != 0L,
    paste(run, "gzip 完整性失败")
  )
  data.table(run_accession = run, conversion = "generated_from_verified_sra")
}

message("并行转换 42 个 SRA 为成对 gzip FASTQ")
conversion <- rbindlist(
  mclapply(seq_len(nrow(pairing)), convert_one, mc.cores = 4L),
  use.names = TRUE
)
fail_if(nrow(conversion) != 42L, "FASTQ 转换结果行数错误")

stage_fastq <- unlist(lapply(pairing$run_accession, function(run) {
  file.path(incoming_dir, paste0(run, c("_1.fastq.gz", "_2.fastq.gz")))
}), use.names = FALSE)
fail_if(length(stage_fastq) != 84L || any(!file.exists(stage_fastq)),
        "未生成 84 个 FASTQ")

message("计算 FASTQ 统计与引物残留")
stats_text <- system2(
  tools[["seqkit"]], c("stats", "-T", "-a", stage_fastq),
  stdout = TRUE, stderr = FALSE
)
stats_status <- attr(stats_text, "status")
if (is.null(stats_status)) stats_status <- 0L
fail_if(stats_status != 0L, "seqkit stats 失败")
fastq_stats <- fread(text = paste(stats_text, collapse = "\n"), check.names = TRUE)
setnames(
  fastq_stats,
  old = c("num_seqs", "sum_len", "min_len", "avg_len", "max_len",
          "Q20...", "Q30...", "AvgQual", "GC...", "sum_n"),
  new = c("read_count", "base_count", "min_length", "mean_length", "max_length",
          "q20_percent", "q30_percent", "mean_quality", "gc_percent", "n_bases"),
  skip_absent = TRUE
)
required_stat_columns <- c(
  "file", "read_count", "base_count", "min_length", "mean_length", "max_length",
  "q20_percent", "q30_percent", "mean_quality", "gc_percent", "n_bases"
)
fail_if(!all(required_stat_columns %chin% names(fastq_stats)),
        "seqkit stats 表头与预期不一致")
fastq_stats[, run_accession := sub("_[12][.]fastq[.]gz$", "", basename(file))]
fastq_stats[, read_direction := fifelse(
  grepl("_1[.]fastq[.]gz$", file), "R1", "R2"
)]

primer_forward <- unique(pairing$primer_forward_sequence)
primer_reverse <- unique(pairing$primer_reverse_sequence)
fail_if(length(primer_forward) != 1L || length(primer_reverse) != 1L,
        "引物序列不唯一")

primer_one <- function(run) {
  r1 <- file.path(incoming_dir, paste0(run, "_1.fastq.gz"))
  r2 <- file.path(incoming_dir, paste0(run, "_2.fastq.gz"))
  json_path <- file.path("/tmp", paste0("prjna766558_primer_", run, ".json"))
  on.exit(if (file.exists(json_path)) unlink(json_path), add = TRUE)
  log <- system2(
    tools[["cutadapt"]],
    c(
      "-g", paste0("^", primer_forward),
      "-G", paste0("^", primer_reverse),
      "--action=none", "--json", json_path,
      "-o", "/dev/null", "-p", "/dev/null", r1, r2
    ),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(log, "status")
  if (is.null(status)) status <- 0L
  fail_if(status != 0L || !file.exists(json_path), paste(run, "cutadapt 引物审计失败"))
  report <- fromJSON(json_path, simplifyVector = FALSE)
  data.table(
    run_accession = run,
    primer_r1_anchored_matches = as.integer(report$read_counts$read1_with_adapter),
    primer_r2_anchored_matches = as.integer(report$read_counts$read2_with_adapter)
  )
}
primer_audit <- rbindlist(
  mclapply(pairing$run_accession, primer_one, mc.cores = 4L),
  use.names = TRUE
)

fastq_qc <- merge(
  fastq_stats[, ..required_stat_columns][, `:=`(
    run_accession = sub("_[12][.]fastq[.]gz$", "", basename(file)),
    read_direction = fifelse(grepl("_1[.]fastq[.]gz$", file), "R1", "R2")
  )],
  pairing[, .(
    run_accession, patient_pair_id, paper_pair_number, tissue_role, sample_name,
    expected_spots = spots, marker, primer_forward_sequence,
    primer_reverse_sequence, ffpe_status, public_negative_control_run
  )],
  by = "run_accession", all.x = TRUE, sort = FALSE
)
fastq_qc <- merge(fastq_qc, primer_audit, by = "run_accession", all.x = TRUE, sort = FALSE)
fastq_qc[, `:=`(
  relative_path = file.path("20_reusable", "fastq", basename(file)),
  read_count_matches_sra_spots = read_count == expected_spots,
  primer_status = fifelse(
    primer_r1_anchored_matches == 0L & primer_r2_anchored_matches == 0L,
    "no_anchored_source_primer_detected_already_trimmed_likely",
    "anchored_source_primer_detected"
  )
)]
setorder(fastq_qc, patient_pair_id, tissue_role, read_direction)

fail_if(nrow(fastq_qc) != 84L || anyNA(fastq_qc$patient_pair_id),
        "FASTQ QC 与样本表无法完整合并")
fail_if(!all(fastq_qc$read_count_matches_sra_spots),
        "FASTQ read 数与 SRA spots 不一致")
fail_if(any(fastq_qc$min_length != fastq_qc$max_length),
        "本冻结版本预期每个 FASTQ 内 read 长度固定")
fail_if(
  !all(fastq_qc[read_direction == "R1", mean_length] == 221) ||
    !all(fastq_qc[read_direction == "R2", mean_length] == 220),
  "FASTQ 固定长度不是预期的 R1=221/R2=220"
)

message("复制校验后的 FASTQ 到中央仓库 L2")
final_fastq <- file.path(final_fastq_dir, basename(stage_fastq))
stage_sha <- vapply(
  stage_fastq, digest, character(1), algo = "sha256", file = TRUE,
  serialize = FALSE
)
for (i in seq_along(stage_fastq)) {
  if (file.exists(final_fastq[[i]])) {
    existing_sha <- digest(
      final_fastq[[i]], algo = "sha256", file = TRUE, serialize = FALSE
    )
    fail_if(existing_sha != stage_sha[[i]],
            paste0("目标 FASTQ 已存在但 SHA256 不同：", final_fastq[[i]]))
  } else {
    file_copy(stage_fastq[[i]], final_fastq[[i]], overwrite = FALSE)
  }
  fail_if(
    digest(final_fastq[[i]], algo = "sha256", file = TRUE, serialize = FALSE) !=
      stage_sha[[i]],
    paste0("FASTQ 发布后 SHA256 不一致：", final_fastq[[i]])
  )
}

fastq_qc[, `:=`(
  file = basename(file),
  file_size_bytes = as.numeric(file.info(final_fastq_dir |>
    file.path(basename(relative_path)))$size),
  sha256 = stage_sha[match(basename(relative_path), basename(stage_fastq))]
)]

provenance <- fread(provenance_path, colClasses = "character", na.strings = NULL)
new_provenance <- data.table(
  output_scope = c(
    "20_reusable/fastq/*_1.fastq.gz;20_reusable/fastq/*_2.fastq.gz",
    "10_metadata/fastq_conversion_qc.tsv"
  ),
  corresponding_source_file = c(
    "00_source/sra/*/*.sra",
    "20_reusable/fastq/*.fastq.gz;10_metadata/prjna766558_pairing.tsv"
  ),
  generating_script = "scripts/17_convert_prjna766558_sra_to_fastq.R",
  key_parameters = c(
    "fasterq-dump --split-files --threads 2 from /tmp; pigz -p 2; preserve full 221/220 bp reads",
    "seqkit full statistics; exact read-count check; cutadapt anchored 336F/806R detection with action=none"
  ),
  software = paste0(
    "R ", getRversion(), "; SRA Toolkit 3.4.1; pigz 2.8; seqkit 2.13.0; cutadapt 5.2"
  ),
  generated_date = as.character(Sys.Date()),
  regenerable = "yes"
)
provenance <- rbind(
  provenance[generating_script != "scripts/17_convert_prjna766558_sra_to_fastq.R"],
  new_provenance,
  fill = TRUE
)
atomic_fwrite(provenance, provenance_path)
atomic_fwrite(fastq_qc, fastq_qc_path)

fastq_manifest <- data.table(
  relative_path = file.path("20_reusable", "fastq", basename(final_fastq)),
  file_level = "20_reusable",
  size_bytes = as.character(as.numeric(file.info(final_fastq)$size)),
  sha256 = stage_sha,
  source_url = "generated from local verified full SRA",
  download_date = as.character(Sys.Date()),
  file_status = "generated_verified",
  corresponding_source_file = file.path(
    "00_source", "sra",
    sub("_[12][.]fastq[.]gz$", "", basename(final_fastq)),
    paste0(sub("_[12][.]fastq[.]gz$", "", basename(final_fastq)), ".sra")
  ),
  generation_method = "fasterq-dump 3.4.1 --split-files; pigz 2.8",
  notes = "read_count_matches_sra_spots;fixed_length_R1_221_R2_220;gzip_test_passed"
)
metadata_paths <- c(fastq_qc_path, provenance_path)
metadata_manifest <- data.table(
  relative_path = c(
    "10_metadata/fastq_conversion_qc.tsv",
    "20_reusable/PROVENANCE.tsv"
  ),
  file_level = c("10_metadata", "20_reusable"),
  size_bytes = as.character(as.numeric(file.info(metadata_paths)$size)),
  sha256 = vapply(
    metadata_paths, digest, character(1), algo = "sha256", file = TRUE,
    serialize = FALSE
  ),
  source_url = "generated from verified local dataset",
  download_date = as.character(Sys.Date()),
  file_status = "generated_verified",
  corresponding_source_file = c(
    "00_source/sra/*/*.sra;10_metadata/prjna766558_pairing.tsv",
    "00_source/sra/*/*.sra"
  ),
  generation_method = "scripts/17_convert_prjna766558_sra_to_fastq.R",
  notes = c("84 FASTQ rows; full read and primer audit", "updated L2 provenance")
)
manifest <- manifest[
  !relative_path %chin% c(fastq_manifest$relative_path, metadata_manifest$relative_path)
]
manifest <- rbindlist(list(manifest, fastq_manifest, metadata_manifest), fill = TRUE)
setorder(manifest, file_level, relative_path)
atomic_fwrite(manifest, manifest_path)

published_manifest <- fread(manifest_path, colClasses = "character", na.strings = NULL)
published_rows <- published_manifest[
  relative_path %chin% c(fastq_manifest$relative_path, metadata_manifest$relative_path)
]
published_paths <- file.path(canonical_root, published_rows$relative_path)
manifest_ok <-
  file.exists(published_paths) &
  as.character(as.numeric(file.info(published_paths)$size)) == published_rows$size_bytes &
  vapply(
    published_paths, digest, character(1), algo = "sha256", file = TRUE,
    serialize = FALSE
  ) == published_rows$sha256
fail_if(!all(manifest_ok) || nrow(published_rows) != 86L,
        "发布后的 84 FASTQ 和 2 个元数据文件未通过 manifest 复核")

qa_lines <- c(
  "# PRJNA766558 FASTQ 转换 QA（2026-07-11）",
  "",
  "本文件是运行时历史证据，不是项目当前状态源。",
  "",
  paste0("- verified SRA 输入：", nrow(pairing), " 个；完整 SHA256 复核通过。"),
  paste0("- FASTQ 输出：", nrow(fastq_qc), " 个文件；", uniqueN(fastq_qc$run_accession), " 个 run。"),
  "- fasterq-dump 固定从 ASCII 路径 `/tmp` 启动，规避中文工作目录下的 SRA Toolkit buffer 错误。",
  "- 每个 run 的 R1/R2 read 数均与 SRA spots 一致；R1 固定 221 bp，R2 固定 220 bp。",
  paste0("- R1 锚定 336F 匹配总数：", sum(fastq_qc$primer_r1_anchored_matches),
         "；R2 锚定 806R 匹配总数：", sum(fastq_qc$primer_r2_anchored_matches), "。"),
  "- 未再次切除引物；正式 DADA2 保留全长以保护 V3–V4 合并重叠。",
  paste0("- Q30 范围：", sprintf("%.2f", min(fastq_qc$q30_percent)), "%–",
         sprintf("%.2f", max(fastq_qc$q30_percent)), "%；无 N 碱基文件数：",
         sum(fastq_qc$n_bases == 0), "/84。"),
  "- 84 个 FASTQ、转换 QC 和 PROVENANCE 已写入中央 manifest 并逐文件复核。",
  "- 原始 SRA 未修改或删除；本项目目录未保存 FASTQ 副本。"
)
writeLines(
  qa_lines,
  file.path(checks_dir, "prjna766558_fastq_conversion_qa_20260711.md"),
  useBytes = TRUE
)

if (dir.exists(incoming_dir)) dir_delete(incoming_dir)
message("PRJNA766558 FASTQ 转换完成：84 个文件，全部通过 read 数、gzip 和 manifest 校验。")
