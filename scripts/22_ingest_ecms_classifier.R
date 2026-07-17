#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(fs)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)

# 将 ECMS 表达型随机森林分类器锁定到作者 GitHub 的固定 commit，
# 先进入 ResearchDataHub/_incoming，完成字节、Git blob、归档成员和许可边界
# 校验后才以新版本目录发布。本脚本不更新项目 data/datasets.tsv。

args <- commandArgs(trailingOnly = TRUE)
plan_only <- "--plan-only" %in% args
unknown_args <- setdiff(args, "--plan-only")
if (length(unknown_args)) {
  stop("未知参数：", paste(unknown_args, collapse = ", "), call. = FALSE)
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
fail_script_location <- function(message) stop(message, call. = FALSE)
if (length(script_argument) != 1L) {
  fail_script_location("无法从 --file 唯一定位 scripts/22_ingest_ecms_classifier.R。")
}
script_path <- normalizePath(
  sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE
)
project_root <- normalizePath(
  file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE
)
project_required <- file.path(project_root, c("PROJECT_INDEX.md", "data/datasets.tsv"))
if (any(!file.exists(project_required))) {
  fail_script_location(paste(
    "脚本定位到的项目根目录缺少权威索引：",
    paste(project_required[!file.exists(project_required)], collapse = ";")
  ))
}
data_root <- Sys.getenv("RESEARCH_DATA_ROOT", path.expand("~/ResearchDataHub"))
catalog_path <- file.path(data_root, "CATALOG.tsv")

repository <- "CityUHK-CompBio/ESCC_CMS"
commit_sha <- "8a49db65f86f2c1e2a50d9bc4cb742ccc868e4c6"
commit_tree_sha <- "5726b305a38632c85e07c7fe41ebb9215ad84116"
version <- paste0("commit_", commit_sha)
dataset_key <- paste0("GITHUB_CITYUHK_COMPUTATIONAL_BIOLOGY_ESCC_CMS_", commit_sha)
retrieval_date <- as.character(Sys.Date())

stage_root <- file.path(
  data_root, "_incoming", "GitHub", "CityUHK-CompBio", "ESCC_CMS", version
)
canonical_root <- file.path(
  data_root, "datasets", "public", "GitHub", "CityUHK-CompBio_ESCC_CMS", version
)

archive_name <- paste0("ESCC_CMS-", commit_sha, ".tar.gz")
archive_url <- paste0(
  "https://github.com/", repository, "/archive/", commit_sha, ".tar.gz"
)
commit_api_url <- paste0(
  "https://api.github.com/repos/", repository, "/commits/", commit_sha
)
tree_api_url <- paste0(
  "https://api.github.com/repos/", repository, "/git/trees/", commit_sha,
  "?recursive=1"
)
raw_base_url <- paste0(
  "https://raw.githubusercontent.com/", repository, "/", commit_sha, "/"
)

tree_expected <- data.table(
  repository_path = c(
    "ECMS.Rmd", "ECMS.model.rdata", "ECMS_2024.png", "README.md",
    "imECMS.Rmd", "rst_final_model.rdata"
  ),
  git_blob_sha1 = c(
    "421e045c7ead11a94db1e1fa303dd5ccf2323e03",
    "b794fd61457ede021b37ab0dc5d4b108ca10101e",
    "d6819f92d5eca636ccf1a57c1114d0ce8ef8b363",
    "e6197e794f334faa24d182fcbc767719a4263c48",
    "9d2b8f92bdd81614bea165d05516a0ba0c85ef12",
    "3bd8e55895f9a0528af8cdefe1f196cb36ae2efc"
  ),
  size_bytes = c(16354, 1608173, 1487584, 4268, 23143, 439717)
)

key_files <- data.table(
  repository_path = c("ECMS.model.rdata", "README.md", "ECMS.Rmd"),
  relative_path = file.path(
    "00_source", c("ECMS.model.rdata", "README.md", "ECMS.Rmd")
  ),
  expected_size_bytes = c(1608173, 4268, 16354),
  expected_sha256 = c(
    "4c05645d712ec630f8ed717cba347dced4228a999c1b479ebd760bf1436e1ee4",
    NA_character_, NA_character_
  ),
  expected_git_blob_sha1 = c(
    "b794fd61457ede021b37ab0dc5d4b108ca10101e",
    "e6197e794f334faa24d182fcbc767719a4263c48",
    "421e045c7ead11a94db1e1fa303dd5ccf2323e03"
  )
)

if (plan_only) {
  cat(paste0(
    "dataset_key\t", dataset_key, "\n",
    "commit\t", commit_sha, "\n",
    "stage_root\t", stage_root, "\n",
    "canonical_root\t", canonical_root, "\n",
    "model_size_bytes\t1608173\n",
    "model_sha256\t4c05645d712ec630f8ed717cba347dced4228a999c1b479ebd760bf1436e1ee4\n",
    "catalog_update\tlatest_snapshot_under_global_lock\n",
    "project_data_index_update\tno\n",
    "license_boundary\tno_LICENSE_file_at_locked_commit\n"
  ))
  quit(save = "no", status = 0L)
}

fail_if <- function(condition, message) {
  if (length(condition) != 1L || is.na(condition) || isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

git_blob_sha1 <- function(path) {
  size <- as.numeric(file_info(path)$size)
  bytes <- readBin(path, what = "raw", n = size)
  fail_if(length(bytes) != size, paste("无法完整读取 Git blob：", path))
  payload <- c(charToRaw(paste0("blob ", size)), as.raw(0L), bytes)
  digest(payload, algo = "sha1", serialize = FALSE)
}

atomic_fwrite <- function(object, path) {
  dir_create(dirname(path), recurse = TRUE)
  temp_path <- tempfile(
    pattern = paste0(".", basename(path), "."), tmpdir = dirname(path)
  )
  on.exit(if (file_exists(temp_path)) file_delete(temp_path), add = TRUE)
  fwrite(object, temp_path, sep = "\t", quote = FALSE, na = "")
  reread <- fread(temp_path, colClasses = "character", na.strings = NULL)
  fail_if(nrow(reread) != nrow(object), paste("原子 TSV 回读失败：", path))
  fail_if(!file.rename(temp_path, path), paste("无法原子更新：", path))
}

atomic_write_lines <- function(lines, path) {
  dir_create(dirname(path), recurse = TRUE)
  temp_path <- tempfile(
    pattern = paste0(".", basename(path), "."), tmpdir = dirname(path)
  )
  on.exit(if (file_exists(temp_path)) file_delete(temp_path), add = TRUE)
  writeLines(lines, temp_path, useBytes = TRUE)
  fail_if(!file.rename(temp_path, path), paste("无法原子更新：", path))
}

atomic_download <- function(url, destination) {
  dir_create(dirname(destination), recurse = TRUE)
  if (file_exists(destination)) return(invisible(destination))
  part_path <- paste0(destination, ".part")
  curl_bin <- Sys.which("curl")
  fail_if(!nzchar(curl_bin), "缺少 curl 命令。")
  status <- system2(
    curl_bin,
    c(
      "-fL", "--retry", "5", "--retry-delay", "2",
      "--connect-timeout", "30", "--continue-at", "-",
      "--output", shQuote(part_path), shQuote(url)
    )
  )
  fail_if(!identical(status, 0L), paste("下载失败：", url))
  fail_if(!file_exists(part_path) || as.numeric(file_info(part_path)$size) <= 0,
          paste("下载产物为空：", url))
  fail_if(!file.rename(part_path, destination),
          paste("无法原子提升下载文件：", destination))
  invisible(destination)
}

copy_verified_atomic <- function(source, destination) {
  dir_create(dirname(destination), recurse = TRUE)
  fail_if(file_exists(destination), paste("目标已存在，拒绝覆盖：", destination))
  temp_path <- tempfile(
    pattern = paste0(".", basename(destination), ".copy."),
    tmpdir = dirname(destination)
  )
  on.exit(if (file_exists(temp_path)) file_delete(temp_path), add = TRUE)
  fail_if(!file.copy(source, temp_path, overwrite = FALSE, copy.mode = TRUE),
          paste("复制失败：", source))
  fail_if(
    as.numeric(file_info(source)$size) != as.numeric(file_info(temp_path)$size) ||
      sha256_file(source) != sha256_file(temp_path),
    paste("复制后大小或 SHA256 不一致：", source)
  )
  fail_if(!file.rename(temp_path, destination),
          paste("无法原子提升复制文件：", destination))
}

verify_manifest_root <- function(root) {
  manifest_path <- file.path(root, "90_manifests", "MANIFEST.tsv")
  fail_if(!file_exists(manifest_path), paste("缺少 manifest：", manifest_path))
  manifest <- fread(manifest_path, colClasses = "character", na.strings = NULL)
  required_fields <- c("relative_path", "size_bytes", "sha256", "file_status")
  fail_if(!all(required_fields %in% names(manifest)), "MANIFEST.tsv 字段不完整。")
  fail_if(!nrow(manifest) || anyDuplicated(manifest$relative_path),
          "MANIFEST.tsv 为空或相对路径重复。")
  paths <- file.path(root, manifest$relative_path)
  fail_if(any(!file_exists(paths)), "manifest 登记文件缺失。")
  observed_sizes <- as.character(as.numeric(file_info(paths)$size))
  observed_sha <- vapply(paths, sha256_file, character(1))
  fail_if(any(observed_sizes != manifest$size_bytes) || any(observed_sha != manifest$sha256),
          "manifest 登记文件大小或 SHA256 失败。")
  observed_relatives <- sort(as.character(path_rel(
    dir_ls(root, recurse = TRUE, type = "file", all = TRUE), start = root
  )))
  expected_relatives <- sort(c(
    manifest$relative_path, "90_manifests/MANIFEST.tsv"
  ))
  fail_if(!identical(observed_relatives, expected_relatives),
          "规范目录存在 manifest 未登记文件或登记文件缺失。")
  model_row <- manifest[relative_path == "00_source/ECMS.model.rdata"]
  fail_if(nrow(model_row) != 1L || model_row$size_bytes != "1608173" ||
            model_row$sha256 !=
              "4c05645d712ec630f8ed717cba347dced4228a999c1b479ebd760bf1436e1ee4",
          "规范目录 ECMS.model.rdata 未通过锁定校验。")
  invisible(manifest)
}

fail_if(!dir_exists(data_root), paste("ResearchDataHub 不可访问：", data_root))
fail_if(!file_exists(catalog_path) || !file_access(catalog_path, "read"),
        paste("ResearchDataHub/CATALOG.tsv 不可读：", catalog_path))

catalog_start_sha256 <- sha256_file(catalog_path)
catalog_start <- fread(catalog_path, colClasses = "character", na.strings = NULL)
fail_if(!"dataset_key" %in% names(catalog_start) || anyDuplicated(catalog_start$dataset_key),
        "CATALOG.tsv 缺少 dataset_key 或已存在重复键。")

dataset_key_value <- dataset_key
existing_at_start <- catalog_start[dataset_key == dataset_key_value]
if (nrow(existing_at_start)) {
  fail_if(
    nrow(existing_at_start) != 1L || existing_at_start$status != "verified" ||
      existing_at_start$local_path != canonical_root,
    "CATALOG.tsv 已存在同 dataset_key，但不是同一 verified 规范记录。"
  )
  verify_manifest_root(canonical_root)
  if (dir_exists(stage_root)) dir_delete(stage_root)
  message("ECMS 固定 commit 已验证入库；未覆盖、未更新项目索引。")
  quit(save = "no", status = 0L)
}

if (!dir_exists(canonical_root)) {
  dir_create(file.path(stage_root, "00_source"), recurse = TRUE)
  dir_create(file.path(stage_root, "10_metadata"), recurse = TRUE)
  dir_create(file.path(stage_root, "20_reusable"), recurse = TRUE)
  dir_create(file.path(stage_root, "90_manifests"), recurse = TRUE)

  message("[1/6] 下载固定 commit 归档、关键原文件与 GitHub 元数据")
  archive_path <- file.path(stage_root, "00_source", archive_name)
  atomic_download(archive_url, archive_path)

  key_files[, path := file.path(stage_root, relative_path)]
  for (i in seq_len(nrow(key_files))) {
    atomic_download(
      paste0(raw_base_url, key_files$repository_path[[i]]),
      key_files$path[[i]]
    )
  }

  commit_json_path <- file.path(stage_root, "10_metadata", "github_commit.json")
  tree_json_path <- file.path(stage_root, "10_metadata", "github_tree.json")
  atomic_download(commit_api_url, commit_json_path)
  atomic_download(tree_api_url, tree_json_path)

  message("[2/6] 校验 commit、tree、Git blob 与锁定模型字节")
  commit_json <- fromJSON(commit_json_path, simplifyVector = FALSE)
  fail_if(!identical(commit_json$sha, commit_sha), "GitHub commit JSON 的 SHA 不符。")
  fail_if(!identical(commit_json$commit$tree$sha, commit_tree_sha),
          "GitHub commit tree SHA 不符。")

  tree_json <- fromJSON(tree_json_path, simplifyVector = TRUE)
  tree_observed <- as.data.table(tree_json$tree)
  fail_if(isTRUE(tree_json$truncated), "GitHub tree API 结果被截断。")
  fail_if(nrow(tree_observed) != nrow(tree_expected),
          "固定 commit 的仓库文件数与预期不符。")
  tree_check <- merge(
    tree_expected, tree_observed[, .(repository_path = path, type, sha, size)],
    by = "repository_path", all = TRUE
  )
  fail_if(
    nrow(tree_check) != nrow(tree_expected) || anyNA(tree_check$git_blob_sha1) ||
      any(tree_check$type != "blob") ||
      any(tree_check$git_blob_sha1 != tree_check$sha) ||
      any(tree_check$size_bytes != tree_check$size),
    "GitHub tree 的路径、blob SHA 或大小不符。"
  )
  fail_if(any(grepl("^(LICENSE|COPYING)([.]|$)", basename(tree_observed$path),
                       ignore.case = TRUE)),
          "锁定 commit 的许可文件状态与预先审计不一致。")

  for (i in seq_len(nrow(key_files))) {
    path <- key_files$path[[i]]
    fail_if(as.numeric(file_info(path)$size) != key_files$expected_size_bytes[[i]],
            paste("关键原文件大小不符：", key_files$repository_path[[i]]))
    fail_if(git_blob_sha1(path) != key_files$expected_git_blob_sha1[[i]],
            paste("关键原文件 Git blob SHA 不符：",
                  key_files$repository_path[[i]]))
    if (!is.na(key_files$expected_sha256[[i]])) {
      fail_if(sha256_file(path) != key_files$expected_sha256[[i]],
              paste("关键原文件 SHA256 不符：",
                    key_files$repository_path[[i]]))
    }
  }

  message("[3/6] 安全校验并解开 commit archive")
  members <- system2("/usr/bin/tar", c("-tzf", shQuote(archive_path)), stdout = TRUE)
  verbose_members <- system2(
    "/usr/bin/tar", c("-tvzf", shQuote(archive_path)), stdout = TRUE
  )
  fail_if(!length(members) || anyDuplicated(members), "commit archive 成员为空或重复。")
  fail_if(any(startsWith(members, "/")) || any(vapply(
    strsplit(members, "/", fixed = TRUE),
    function(parts) any(parts == ".."), logical(1)
  )), "commit archive 包含不安全路径。")
  fail_if(length(verbose_members) != length(members) ||
            any(!substr(verbose_members, 1L, 1L) %in% c("-", "d")),
          "commit archive 包含非常规文件、符号链接或硬链接。")
  archive_roots <- unique(sub("/.*$", "", members))
  fail_if(length(archive_roots) != 1L, "commit archive 顶层目录不唯一。")
  expected_members <- c(
    paste0(archive_roots, "/"),
    file.path(archive_roots, tree_expected$repository_path)
  )
  fail_if(!setequal(members, expected_members), "commit archive 成员集与 Git tree 不一致。")

  snapshot_dir <- file.path(stage_root, "20_reusable", "repository_snapshot")
  if (dir_exists(snapshot_dir)) dir_delete(snapshot_dir)
  extract_temp <- file.path(
    stage_root, "20_reusable", paste0(".extract-", Sys.getpid())
  )
  if (dir_exists(extract_temp)) dir_delete(extract_temp)
  dir_create(extract_temp)
  utils::untar(archive_path, exdir = extract_temp)
  extracted_root <- file.path(extract_temp, archive_roots)
  fail_if(!dir_exists(extracted_root), "commit archive 解包根目录缺失。")
  fail_if(!file.rename(extracted_root, snapshot_dir),
          "无法原子提升 repository snapshot。")
  if (dir_exists(extract_temp)) dir_delete(extract_temp)

  snapshot_files <- file.path(snapshot_dir, tree_expected$repository_path)
  fail_if(any(!file_exists(snapshot_files)), "repository snapshot 文件缺失。")
  snapshot_inventory <- copy(tree_expected)
  snapshot_inventory[, `:=`(
    relative_path = file.path("20_reusable", "repository_snapshot", repository_path),
    observed_size_bytes = as.numeric(file_info(snapshot_files)$size),
    observed_git_blob_sha1 = vapply(snapshot_files, git_blob_sha1, character(1)),
    sha256 = vapply(snapshot_files, sha256_file, character(1)),
    integrity_status = "verified_against_locked_git_tree"
  )]
  fail_if(
    any(snapshot_inventory$observed_size_bytes != snapshot_inventory$size_bytes) ||
      any(snapshot_inventory$observed_git_blob_sha1 != snapshot_inventory$git_blob_sha1),
    "repository snapshot 与锁定 Git tree 不一致。"
  )

  model_snapshot_path <- file.path(snapshot_dir, "ECMS.model.rdata")
  fail_if(
    sha256_file(model_snapshot_path) !=
      "4c05645d712ec630f8ed717cba347dced4228a999c1b479ebd760bf1436e1ee4",
    "archive 内 ECMS.model.rdata SHA256 不符。"
  )

  message("[4/6] 生成中文数据说明、生成关系与 manifest")
  inventory_path <- file.path(
    stage_root, "10_metadata", "repository_file_inventory.tsv"
  )
  atomic_fwrite(snapshot_inventory, inventory_path)

  provenance <- data.table(
    output_scope = c(
      "20_reusable/repository_snapshot/*",
      "10_metadata/repository_file_inventory.tsv"
    ),
    corresponding_source_file = c(
      file.path("00_source", archive_name),
      "10_metadata/github_tree.json"
    ),
    input_sha256 = c(sha256_file(archive_path), sha256_file(tree_json_path)),
    generating_script = "scripts/22_ingest_ecms_classifier.R",
    key_parameters = c(
      paste0("locked_commit=", commit_sha,
             ";exact Git tree;safe tar paths;no links;Git blob verification"),
      paste0("locked_commit=", commit_sha, ";recursive Git tree API")
    ),
    software = paste0(
      "R ", getRversion(), "; data.table ", packageVersion("data.table"),
      "; digest ", packageVersion("digest"), "; jsonlite ",
      packageVersion("jsonlite")
    ),
    generated_date = retrieval_date,
    regenerable = "yes",
    catalog_sha256_at_script_start = catalog_start_sha256
  )
  provenance_path <- file.path(stage_root, "20_reusable", "PROVENANCE.tsv")
  atomic_fwrite(provenance, provenance_path)

  dataset_md <- c(
    "# ESCC ECMS 表达型分类器固定版本说明",
    "",
    "## 基本信息",
    "",
    paste0("- `dataset_key`：`", dataset_key, "`。"),
    paste0("- 来源：GitHub `", repository, "`。"),
    paste0("- 固定 commit：`", commit_sha, "`；Git tree：`",
           commit_tree_sha, "`。"),
    "- 资源类型：ESCC 共识分子分型（ECMS）表达型 randomForest 模型与作者代码快照。",
    "- 原始论文 DOI：`10.1038/s41392-026-02577-9`。",
    paste0("- 固定仓库链接：https://github.com/", repository, "/tree/", commit_sha),
    "- 模型训练队列 152 例；模型包内置 TCGA 78、GSE53625 179 和 GSE45670 28 例验证矩阵。",
    paste0("- 当前使用项目：`", project_root, "`。"),
    paste0("- 获取日期：", retrieval_date, "；完整性状态：`verified`。"),
    "",
    "## 目录与生成关系",
    "",
    paste0("- `00_source/", archive_name, "`：固定 commit 的 GitHub archive。"),
    "- `00_source/ECMS.model.rdata`：表达型模型原始字节，大小 1,608,173 B，SHA256 `4c05645d712ec630f8ed717cba347dced4228a999c1b479ebd760bf1436e1ee4`。",
    "- `00_source/README.md` 与 `00_source/ECMS.Rmd`：作者使用说明和分析代码原文件。",
    "- `10_metadata/github_commit.json` 与 `github_tree.json`：固定 commit/tree 的 GitHub API 元数据。",
    "- `20_reusable/repository_snapshot/`：从已验证 archive 安全解开的完整仓库快照；每个文件与 Git blob SHA-1 对照。",
    "- `90_manifests/MANIFEST.tsv`：逐文件大小、SHA256、来源和生成方式。",
    "",
    "## 推荐用途与限制",
    "",
    "- 用于多样本 ESCC 队列的 ECMS1–4 表达型外部基准投影；输入必须按作者 README 先 `log2(TPM+1)`，再在待投影队列内逐基因 Z 标准化。",
    "- 模型不包含训练集固定均值/标准差，不是独立单样本分类器；队列组成可影响标签。",
    "- 模型没有官方最低概率、margin 或拒绝分类阈值；项目中的 margin 只能作自定义软不确定性。",
    "- 该资源用于公共数据假设生成与分型校准，不构成临床诊断、疗效预测或因果证明。",
    "",
    "## 许可、引用与分发边界",
    "",
    "- 锁定 commit 的 Git tree 中没有 `LICENSE` 或 `COPYING` 文件，GitHub 仓库的公开可访问不等于授予软件再分发许可。",
    "- 使用时应引用原始 ECMS 论文与作者仓库；对外再分发、派生部署或商业使用前应向权利人确认授权。",
    "- 论文正文/补充材料的许可不自动推定完整覆盖仓库代码与序列化模型。",
    "",
    "## 完整性",
    "",
    "- 完整固定信息见 `10_metadata/repository_file_inventory.tsv`。",
    "- 生成关系见 `20_reusable/PROVENANCE.tsv`。",
    "- 本数据集不包含论文全文，不触发文献全文入 Zotero。"
  )
  dataset_md_path <- file.path(stage_root, "DATASET.md")
  atomic_write_lines(dataset_md, dataset_md_path)

  manifest_files <- c(
    archive_path,
    key_files$path,
    commit_json_path,
    tree_json_path,
    inventory_path,
    snapshot_files,
    provenance_path,
    dataset_md_path
  )
  relative_paths <- as.character(path_rel(manifest_files, start = stage_root))
  manifest <- data.table(
    relative_path = relative_paths,
    file_level = fifelse(grepl("/", relative_paths), sub("/.*$", "", relative_paths), "root"),
    size_bytes = as.numeric(file_info(manifest_files)$size),
    sha256 = vapply(manifest_files, sha256_file, character(1)),
    source_url = "generated within locked repository ingest",
    download_date = retrieval_date,
    file_status = "generated_verified",
    corresponding_source_file = file.path("00_source", archive_name),
    generation_method = "scripts/22_ingest_ecms_classifier.R",
    notes = ""
  )
  manifest[relative_path == file.path("00_source", archive_name), `:=`(
    source_url = archive_url,
    file_status = "verified_locked_commit_archive",
    corresponding_source_file = "",
    generation_method = "GitHub commit archive source bytes",
    notes = paste0("commit=", commit_sha, ";members_match_locked_git_tree")
  )]
  for (i in seq_len(nrow(key_files))) {
    manifest[relative_path == key_files$relative_path[[i]], `:=`(
      source_url = paste0(raw_base_url, key_files$repository_path[[i]]),
      file_status = "verified_locked_git_blob",
      corresponding_source_file = "",
      generation_method = "GitHub raw source bytes",
      notes = paste0("git_blob_sha1=", key_files$expected_git_blob_sha1[[i]])
    )]
  }
  manifest[relative_path == "10_metadata/github_commit.json", `:=`(
    source_url = commit_api_url,
    file_status = "verified_api_metadata",
    corresponding_source_file = "",
    generation_method = "GitHub REST API response",
    notes = paste0("commit=", commit_sha)
  )]
  manifest[relative_path == "10_metadata/github_tree.json", `:=`(
    source_url = tree_api_url,
    file_status = "verified_api_metadata",
    corresponding_source_file = "",
    generation_method = "GitHub REST API response",
    notes = paste0("tree=", commit_tree_sha, ";truncated=false")
  )]
  manifest[startsWith(relative_path, "20_reusable/repository_snapshot/"), `:=`(
    source_url = archive_url,
    file_status = "generated_verified_git_blob",
    corresponding_source_file = file.path("00_source", archive_name),
    generation_method = "safe extraction; content unchanged; Git blob verified",
    notes = "no symlink or hardlink"
  )]
  setorder(manifest, relative_path)
  manifest_path <- file.path(stage_root, "90_manifests", "MANIFEST.tsv")
  atomic_fwrite(manifest, manifest_path)

  observed_stage_files <- sort(as.character(path_rel(
    dir_ls(stage_root, recurse = TRUE, type = "file", all = TRUE),
    start = stage_root
  )))
  expected_stage_files <- sort(c(manifest$relative_path, "90_manifests/MANIFEST.tsv"))
  fail_if(!identical(observed_stage_files, expected_stage_files),
          "_incoming 中存在未登记、缺失或残留 .part 文件。")

  message("[5/6] 以同目录临时路径发布新规范版本")
  canonical_parent <- dirname(canonical_root)
  dir_create(canonical_parent, recurse = TRUE)
  publish_temp <- file.path(
    canonical_parent, paste0(".", basename(canonical_root), ".publishing-", Sys.getpid())
  )
  fail_if(dir_exists(publish_temp), paste("发布临时目录已存在：", publish_temp))
  dir_create(publish_temp)
  publish_failed <- TRUE
  on.exit({
    if (publish_failed && dir_exists(publish_temp)) dir_delete(publish_temp)
  }, add = TRUE)

  publish_relatives <- expected_stage_files
  for (relative in publish_relatives) {
    copy_verified_atomic(
      file.path(stage_root, relative),
      file.path(publish_temp, relative)
    )
  }
  verify_manifest_root(publish_temp)
  fail_if(dir_exists(canonical_root), paste("规范版本目录已存在：", canonical_root))
  fail_if(!file.rename(publish_temp, canonical_root),
          "无法原子提升新规范版本目录。")
  publish_failed <- FALSE
} else {
  message("发现未登记的同版本规范目录；先执行全 manifest 验证。")
  verify_manifest_root(canonical_root)
}

message("[6/6] 在全局锁内重读最新 CATALOG 并唯一登记")
canonical_manifest_path <- file.path(canonical_root, "90_manifests", "MANIFEST.tsv")
verify_manifest_root(canonical_root)

catalog_row <- data.table(
  dataset_key = dataset_key,
  record_type = "locked_public_classifier_model_and_code_snapshot",
  access_level = "public",
  source = "GitHub",
  accession = repository,
  version = version,
  species = "Homo sapiens",
  disease = "esophageal squamous cell carcinoma",
  tissue = "bulk tumor transcriptome classifier",
  assay = "314-gene expression randomForest ECMS classifier; author code snapshot",
  sample_summary = "serialized randomForest model trained on 152 ESCC samples; bundled TCGA 78, GSE53625 179 and GSE45670 28 validation matrices",
  available_levels = "00_source;10_metadata;20_reusable",
  status = "verified",
  local_path = canonical_root,
  manifest_path = canonical_manifest_path,
  source_url = paste0("https://github.com/", repository, "/tree/", commit_sha),
  download_date = retrieval_date,
  last_verified = retrieval_date,
  license_or_access = "publicly accessible GitHub repository; no LICENSE/COPYING file at locked commit; cite paper/repository and obtain permission before redistribution or derivative deployment",
  projects_using = project_root,
  recommended_use = "batch-cohort ECMS1-4 transcriptomic benchmark after log2(TPM+1) and within-cohort gene-wise scaling",
  limitations = "no training-set centering/scaling parameters; not valid as an independent single-sample classifier; no official probability/margin rejection threshold; no complete retraining environment; no repository LICENSE",
  notes = paste0(
    "commit=", commit_sha,
    ";model_size=1608173;model_sha256=4c05645d712ec630f8ed717cba347dced4228a999c1b479ebd760bf1436e1ee4;project data/datasets.tsv intentionally not updated by ingest script"
  )
)

catalog_lock <- paste0(catalog_path, ".update.lock")
# 顶层 Rscript 中的 on.exit() 不提供可靠的脚本退出清理。把锁生命周期
# 放进立即执行函数，使 on.exit 真正绑定函数栈；任何门禁失败或中断都会
# 清理空锁目录，避免后续合法入库被遗留锁永久阻塞。
catalog_lock_result <- (function() {
  fail_if(!dir.create(catalog_lock, showWarnings = FALSE),
          paste("CATALOG.tsv 正在被其他进程更新：", catalog_lock))
  on.exit(unlink(catalog_lock, recursive = TRUE, force = TRUE), add = TRUE)

  catalog_preupdate_sha256_local <- sha256_file(catalog_path)
  catalog_latest <- fread(
    catalog_path, colClasses = "character", na.strings = NULL
  )
  fail_if(!identical(names(catalog_latest), names(catalog_row)),
          "CATALOG.tsv 字段顺序与预期 schema 不一致。")
  fail_if(anyDuplicated(catalog_latest$dataset_key),
          "最新 CATALOG.tsv 已存在重复 dataset_key。")

  existing_latest <- catalog_latest[dataset_key == dataset_key_value]
  if (!nrow(existing_latest)) {
    merged_catalog <- rbindlist(
      list(catalog_latest, catalog_row), use.names = TRUE, fill = FALSE
    )
    fail_if(anyDuplicated(merged_catalog$dataset_key),
            "合并最新 CATALOG 后 dataset_key 不唯一。")
    atomic_fwrite(merged_catalog, catalog_path)
  } else {
    fail_if(
      nrow(existing_latest) != 1L || existing_latest$status != "verified" ||
        existing_latest$local_path != canonical_root ||
        existing_latest$manifest_path != canonical_manifest_path,
      "并发更新后 CATALOG 已出现同 dataset_key 的冲突记录。"
    )
  }

  catalog_final <- fread(
    catalog_path, colClasses = "character", na.strings = NULL
  )
  catalog_final_row <- catalog_final[dataset_key == dataset_key_value]
  fail_if(
    nrow(catalog_final_row) != 1L ||
      catalog_final_row$status != "verified" ||
      catalog_final_row$local_path != canonical_root,
    "CATALOG.tsv 原子更新后回读失败。"
  )
  fail_if(anyDuplicated(catalog_final$dataset_key),
          "CATALOG.tsv 发布后出现重复 dataset_key。")
  list(
    preupdate_sha256 = catalog_preupdate_sha256_local,
    final_row = catalog_final_row
  )
})()
catalog_preupdate_sha256 <- catalog_lock_result$preupdate_sha256

# _incoming 只是入库暂存区。仅在规范目录全 manifest 与 CATALOG
# 回读均通过后，清理本固定 commit 的唯一 stage；不触碰任何其他原始资料。
verify_manifest_root(canonical_root)
if (dir_exists(stage_root)) dir_delete(stage_root)

if (!identical(catalog_start_sha256, catalog_preupdate_sha256)) {
  message(
    "CATALOG.tsv 在脚本运行期间发生变化；已在全局锁内基于最新表合并，",
    "未以起始快照覆盖并发新增行。"
  )
}

message(
  "完成：ECMS 表达型分类器固定到 commit ", commit_sha,
  "，规范路径全文件 SHA256 已验证，CATALOG 已并发安全登记。",
  "未更新 data/datasets.tsv、PROJECT_INDEX.md 或 results/。"
)
