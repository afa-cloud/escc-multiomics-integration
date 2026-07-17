#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ANCOMBC)
  library(Biostrings)
  library(dada2)
  library(data.table)
  library(digest)
  library(fs)
  library(permute)
  library(phyloseq)
  library(vegan)
})

options(stringsAsFactors = FALSE)

# PRJNA766558：21 对 ESCC 肿瘤/配对癌旁 FFPE 组织 16S V3-V4。
#
# 主流程：不重复切引物、不预截断、maxEE=2/2、minOverlap=12；
# 敏感性：minOverlap=8 与 forward-only。患者对是唯一推断单位。
# 无公开阴性对照，因此不能声称已排除低生物量/FFPE 污染。

script_argument <- grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)
if (length(script_argument) != 1L) {
  stop("无法唯一定位当前 DADA2 脚本。", call. = FALSE)
}
script_path <- normalizePath(
  sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE
)
executed_script_sha256 <- digest(
  script_path,
  algo = "sha256",
  file = TRUE,
  serialize = FALSE
)
pipeline_schema_version <- "PRJNA766558_DADA2_v2_bound_signature_atomic_publish"
species_execution_schema_version <-
  "PRJNA766558_DADA2_species_length_grouped_isolated_chunks_v2"
project_root <- dirname(dirname(script_path))
results_dir <- file.path(project_root, "results")
work_intermediate_dir <- file.path(project_root, "_work", "intermediate")
datasets_file <- file.path(project_root, "data", "datasets.tsv")

fail_if <- function(condition, message) {
  if (!is.logical(condition) || length(condition) != 1L || is.na(condition)) {
    stop(
      paste0(message, "（门禁条件为 NA 或非单一逻辑值，按失败处理）"),
      call. = FALSE
    )
  }
  if (condition) stop(message, call. = FALSE)
}

require_columns <- function(object, columns, object_name) {
  missing_columns <- setdiff(columns, names(object))
  fail_if(
    length(missing_columns) > 0L,
    paste0(object_name, " 缺少字段：", paste(missing_columns, collapse = ", "))
  )
}

require_list_fields <- function(object, fields, object_name) {
  fail_if(!is.list(object), paste0(object_name, " 不是 list"))
  missing_fields <- setdiff(fields, names(object))
  fail_if(
    length(missing_fields) > 0L,
    paste0(object_name, " 缺少字段：", paste(missing_fields, collapse = ", "))
  )
}

stage_tsv <- function(object, filename) {
  fwrite(object, file.path(stage_dir, filename), sep = "\t", na = "")
}

atomic_save_rds <- function(object, path, validator, compress = "xz") {
  dir_create(dirname(path), recurse = TRUE)
  temporary_path <- tempfile(
    pattern = paste0(".", basename(path), "_"),
    tmpdir = dirname(path),
    fileext = ".tmp.rds"
  )
  keep_temporary <- TRUE
  on.exit({
    if (keep_temporary && file_exists(temporary_path)) {
      file_delete(temporary_path)
    }
  }, add = TRUE)
  validator(object)
  saveRDS(object, temporary_path, compress = compress)
  validator(readRDS(temporary_path))
  renamed <- file.rename(temporary_path, path)
  fail_if(!renamed, paste0("checkpoint 原子替换失败：", path))
  keep_temporary <- FALSE
  invisible(path)
}

atomic_publish_file <- function(source, destination) {
  dir_create(dirname(destination), recurse = TRUE)
  temporary_path <- tempfile(
    pattern = paste0(".", basename(destination), "_"),
    tmpdir = dirname(destination),
    fileext = ".tmp"
  )
  keep_temporary <- TRUE
  on.exit({
    if (keep_temporary && file_exists(temporary_path)) {
      file_delete(temporary_path)
    }
  }, add = TRUE)
  file_copy(source, temporary_path, overwrite = TRUE)
  fail_if(
    as.numeric(file_info(source)$size) != as.numeric(file_info(temporary_path)$size) ||
      digest(source, algo = "sha256", file = TRUE, serialize = FALSE) !=
        digest(temporary_path, algo = "sha256", file = TRUE, serialize = FALSE),
    paste0("发布临时副本校验失败：", basename(destination))
  )
  renamed <- file.rename(temporary_path, destination)
  fail_if(!renamed, paste0("结果原子替换失败：", destination))
  keep_temporary <- FALSE
  invisible(destination)
}

get_n <- function(x) sum(getUniques(x))

paired_rank_biserial <- function(differences) {
  differences <- differences[is.finite(differences) & differences != 0]
  if (!length(differences)) return(NA_real_)
  ranks <- rank(abs(differences), ties.method = "average")
  sum(sign(differences) * ranks) / sum(ranks)
}

bootstrap_paired_median <- function(differences, n_bootstrap = 1000L,
                                    seed = 20260711L) {
  differences <- differences[is.finite(differences)]
  if (!length(differences)) {
    return(list(
      valid = 0L,
      ci_low = NA_real_,
      ci_high = NA_real_,
      direction_consistency = NA_real_
    ))
  }
  set.seed(seed)
  estimates <- replicate(
    n_bootstrap,
    median(sample(differences, length(differences), replace = TRUE))
  )
  observed <- median(differences)
  list(
    valid = length(estimates),
    ci_low = unname(quantile(estimates, 0.025, names = FALSE)),
    ci_high = unname(quantile(estimates, 0.975, names = FALSE)),
    direction_consistency = if (observed == 0) {
      NA_real_
    } else {
      mean(sign(estimates) == sign(observed))
    }
  )
}

taxonomy_label <- function(taxonomy) {
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  taxonomy <- as.matrix(taxonomy)
  missing_ranks <- setdiff(ranks, colnames(taxonomy))
  if (length(missing_ranks)) {
    taxonomy <- cbind(
      taxonomy,
      matrix(
        NA_character_,
        nrow = nrow(taxonomy),
        ncol = length(missing_ranks),
        dimnames = list(rownames(taxonomy), missing_ranks)
      )
    )
  }
  taxonomy <- taxonomy[, ranks, drop = FALSE]
  labels <- character(nrow(taxonomy))
  for (i in seq_len(nrow(taxonomy))) {
    if (!is.na(taxonomy[i, "Genus"]) && nzchar(taxonomy[i, "Genus"])) {
      labels[[i]] <- paste0("g__", taxonomy[i, "Genus"])
    } else {
      available <- ranks[!is.na(taxonomy[i, ranks]) &
                           nzchar(taxonomy[i, ranks])]
      if (length(available)) {
        deepest <- tail(available, 1L)
        labels[[i]] <- paste0(
          substr(deepest, 1L, 1L), "__", taxonomy[i, deepest],
          "_unclassified"
        )
      } else {
        labels[[i]] <- "unclassified_prokaryote"
      }
    }
  }
  labels
}

filter_prokaryotic_taxa <- function(sequence_table, taxonomy) {
  taxonomy <- as.matrix(taxonomy)
  taxonomy_text <- apply(taxonomy, 1L, function(x) {
    paste(x[!is.na(x)], collapse = ";")
  })
  kingdom <- if ("Kingdom" %in% colnames(taxonomy)) {
    taxonomy[, "Kingdom"]
  } else {
    rep(NA_character_, nrow(taxonomy))
  }
  prokaryote <- kingdom %in% c("Bacteria", "Archaea")
  organelle <- grepl("Chloroplast|Mitochondria", taxonomy_text, ignore.case = TRUE)
  keep <- prokaryote & !organelle
  fail_if(!any(keep), "分类后没有可保留的细菌/古菌 ASV")
  list(
    sequence_table = sequence_table[, rownames(taxonomy)[keep], drop = FALSE],
    taxonomy = taxonomy[keep, , drop = FALSE],
    assigned_prokaryote_taxa = sum(prokaryote),
    removed_non_prokaryote_taxa = sum(!prokaryote),
    removed_organelle_taxa = sum(prokaryote & organelle)
  )
}

aggregate_genus <- function(sequence_table, taxonomy) {
  fail_if(!identical(colnames(sequence_table), rownames(taxonomy)),
          "ASV 表与 taxonomy 顺序不一致")
  labels <- taxonomy_label(taxonomy)
  genus_names <- sort(unique(labels))
  genus_table <- vapply(genus_names, function(label) {
    rowSums(sequence_table[, labels == label, drop = FALSE])
  }, FUN.VALUE = numeric(nrow(sequence_table)))
  if (is.null(dim(genus_table))) {
    genus_table <- matrix(
      genus_table,
      ncol = 1L,
      dimnames = list(rownames(sequence_table), genus_names)
    )
  } else {
    rownames(genus_table) <- rownames(sequence_table)
    colnames(genus_table) <- genus_names
  }
  genus_table
}

assign_taxonomy_chunked <- function(sequences, reference,
                                    chunk_size = 1000L,
                                    threads = 8L) {
  sequences <- unique(as.character(sequences))
  fail_if(!length(sequences), "没有可分类序列")
  chunks <- split(
    seq_along(sequences),
    ceiling(seq_along(sequences) / chunk_size)
  )
  output <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    message(
      "  taxonomy chunk ", i, "/", length(chunks),
      "；sequences=", length(chunks[[i]])
    )
    output[[i]] <- assignTaxonomy(
      sequences[chunks[[i]]],
      reference,
      minBoot = 50,
      tryRC = TRUE,
      multithread = threads,
      verbose = FALSE
    )
    gc(verbose = FALSE)
  }
  taxonomy <- do.call(rbind, output)
  fail_if(nrow(taxonomy) != length(sequences) ||
            !identical(rownames(taxonomy), sequences),
          "分块 taxonomy 回填顺序错误")
  taxonomy
}

run_species_worker_isolated <- function(taxonomy_chunk, reference, inner_n) {
  worker <- NULL
  worker_reaped <- FALSE
  on.exit({
    if (!is.null(worker) && !worker_reaped) {
      suspendInterrupts({
        try(tools::pskill(worker$pid, signal = tools::SIGTERM), silent = TRUE)
        reaped <- NULL
        deadline <- Sys.time() + 5
        while (is.null(reaped) && Sys.time() < deadline) {
          candidate <- try(
            parallel::mccollect(worker, wait = FALSE),
            silent = TRUE
          )
          if (!inherits(candidate, "try-error") && !is.null(candidate)) {
            reaped <- candidate
          } else {
            Sys.sleep(0.1)
          }
        }
        if (is.null(reaped)) {
          try(tools::pskill(worker$pid, signal = tools::SIGKILL), silent = TRUE)
          try(parallel::mccollect(worker, wait = TRUE), silent = TRUE)
        }
      })
    }
  }, add = TRUE)
  # 清理保护先注册，再创建 fork，消除 spawn 与 cleanup 注册之间的窗口。
  worker <- parallel::mcparallel(
    addSpecies(
      taxonomy_chunk,
      reference,
      allowMultiple = TRUE,
      tryRC = TRUE,
      n = inner_n,
      verbose = FALSE
    ),
    mc.set.seed = FALSE,
    silent = TRUE
  )
  collected <- parallel::mccollect(worker, wait = TRUE)
  worker_reaped <- !is.null(collected)
  collected
}

add_species_chunked_isolated <- function(
    taxonomy, reference, cache_directory, pipeline_signature,
    main_sequence_hash, main_taxonomy_hash, species_execution_signature,
    outer_chunk_size = 1000L, inner_n = 250L) {
  fail_if(!is.matrix(taxonomy) || nrow(taxonomy) < 1L,
          "species 输入 taxonomy 不是非空 matrix")
  fail_if(is.null(rownames(taxonomy)) || anyDuplicated(rownames(taxonomy)) > 0L,
          "species 输入 taxonomy 序列名缺失或重复")
  fail_if(!file_exists(reference), "species 参考文件不存在")
  fail_if(
    !is.numeric(outer_chunk_size) || length(outer_chunk_size) != 1L ||
      !is.finite(outer_chunk_size) || outer_chunk_size < 1L ||
      outer_chunk_size != as.integer(outer_chunk_size),
    "species 外层分块大小无效"
  )
  fail_if(
    !is.numeric(inner_n) || length(inner_n) != 1L ||
      !is.finite(inner_n) || inner_n < 1L || inner_n != as.integer(inner_n),
    "species assignSpecies 内层 n 无效"
  )
  outer_chunk_size <- as.integer(outer_chunk_size)
  inner_n <- as.integer(inner_n)
  species_chunk_dir <- file.path(
    cache_directory, "species_chunks", species_execution_signature
  )
  dir_create(species_chunk_dir, recurse = TRUE)
  # DADA2 assignSpecies 会按 query 长度逐组扫描完整参考。若按原行顺序
  # 简单切块，同一长度会散落到多个块而重复扫描。这里把每个长度组保持
  # 在同一块内，再按 outer_chunk_size 打包；只有单个长度组自身超过上限
  # 时才不可避免地拆分，使总参考扫描次数接近唯一长度数。
  sequence_length <- nchar(rownames(taxonomy))
  length_groups <- split(
    seq_len(nrow(taxonomy)),
    factor(sequence_length, levels = sort(unique(sequence_length)))
  )
  chunks <- list()
  current_chunk <- integer()
  for (length_group in length_groups) {
    if (length(length_group) > outer_chunk_size) {
      if (length(current_chunk)) {
        chunks[[length(chunks) + 1L]] <- current_chunk
        current_chunk <- integer()
      }
      split_large_group <- split(
        length_group,
        ceiling(seq_along(length_group) / outer_chunk_size)
      )
      chunks <- c(chunks, unname(split_large_group))
    } else if (
      length(current_chunk) > 0L &&
        length(current_chunk) + length(length_group) > outer_chunk_size
    ) {
      chunks[[length(chunks) + 1L]] <- current_chunk
      current_chunk <- length_group
    } else {
      current_chunk <- c(current_chunk, length_group)
    }
  }
  if (length(current_chunk)) chunks[[length(chunks) + 1L]] <- current_chunk
  fail_if(
    !identical(sort(unlist(chunks, use.names = FALSE)), seq_len(nrow(taxonomy))) ||
      any(vapply(chunks, length, integer(1)) > outer_chunk_size),
    "species 长度分组打包未覆盖全部序列或超过块上限"
  )
  output <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    taxonomy_chunk <- taxonomy[chunks[[i]], , drop = FALSE]
    chunk_sequence_hash <- digest(
      paste(rownames(taxonomy_chunk), collapse = "|"),
      algo = "sha256",
      serialize = FALSE
    )
    chunk_taxonomy_hash <- digest(
      taxonomy_chunk,
      algo = "sha256",
      serialize = TRUE
    )
    checkpoint_file <- file.path(
      species_chunk_dir,
      sprintf(
        "species_chunk_%03d_of_%03d.rds",
        i, length(chunks)
      )
    )
    validate_chunk_checkpoint <- local({
      expected_index <- as.integer(i)
      expected_total <- as.integer(length(chunks))
      expected_sequences <- rownames(taxonomy_chunk)
      expected_taxonomy <- taxonomy_chunk
      expected_sequence_hash <- chunk_sequence_hash
      expected_taxonomy_hash <- chunk_taxonomy_hash
      expected_columns <- c(colnames(taxonomy_chunk), "Species")
      expected_outer_chunk_size <- outer_chunk_size
      expected_inner_n <- inner_n
      function(checkpoint) {
        require_list_fields(
          checkpoint,
          c(
            "pipeline_signature", "species_execution_signature",
            "main_sequence_hash", "main_taxonomy_hash", "chunk_index",
            "chunk_total", "chunk_sequence_hash", "chunk_taxonomy_hash",
            "outer_chunk_size", "inner_n", "species_output_hash",
            "taxonomy_chunk"
          ),
          "species chunk checkpoint"
        )
        fail_if(
          !identical(checkpoint$pipeline_signature, pipeline_signature) ||
            !identical(
              checkpoint$species_execution_signature,
              species_execution_signature
            ) ||
            !identical(checkpoint$main_sequence_hash, main_sequence_hash) ||
            !identical(checkpoint$main_taxonomy_hash, main_taxonomy_hash) ||
            !identical(as.integer(checkpoint$chunk_index), expected_index) ||
            !identical(as.integer(checkpoint$chunk_total), expected_total) ||
            !identical(
              checkpoint$chunk_sequence_hash,
              expected_sequence_hash
            ) ||
            !identical(
              checkpoint$chunk_taxonomy_hash,
              expected_taxonomy_hash
            ) ||
            !identical(
              as.integer(checkpoint$outer_chunk_size),
              expected_outer_chunk_size
            ) ||
            !identical(as.integer(checkpoint$inner_n), expected_inner_n),
          paste0("species chunk ", expected_index, " 签名/索引不一致")
        )
        value <- checkpoint$taxonomy_chunk
        fail_if(!is.matrix(value),
                paste0("species chunk ", expected_index, " 不是 matrix"))
        fail_if(
          !identical(rownames(value), expected_sequences) ||
            anyDuplicated(rownames(value)) > 0L ||
            !identical(colnames(value), expected_columns),
          paste0("species chunk ", expected_index, " 行列结构不一致")
        )
        fail_if(
          !identical(
            value[, colnames(expected_taxonomy), drop = FALSE],
            expected_taxonomy
          ),
          paste0("species chunk ", expected_index, " 原 taxonomy 被改变")
        )
        fail_if(
          !identical(
            checkpoint$species_output_hash,
            digest(value, algo = "sha256", serialize = TRUE)
          ) || any(!is.na(value[, "Species"]) & !nzchar(value[, "Species"])),
          paste0("species chunk ", expected_index, " 输出 hash 或 Species 值异常")
        )
        invisible(TRUE)
      }
    })
    message(
      "  species chunk ", i, "/", length(chunks),
      "；sequences=", length(chunks[[i]]),
      if (file_exists(checkpoint_file)) "；读取 checkpoint" else "；隔离子进程"
    )
    if (file_exists(checkpoint_file)) {
      chunk_checkpoint <- readRDS(checkpoint_file)
      validate_chunk_checkpoint(chunk_checkpoint)
    } else {
      collected <- tryCatch(
        run_species_worker_isolated(taxonomy_chunk, reference, inner_n),
        error = function(condition) condition
      )
      fail_if(
        inherits(collected, "error") || is.null(collected) ||
          length(collected) != 1L,
        paste0(
          "species chunk ", i,
          " 子进程未返回唯一结果",
          if (inherits(collected, "error")) {
            paste0("：", conditionMessage(collected))
          } else {
            ""
          }
        )
      )
      chunk_result <- collected[[1L]]
      fail_if(
        inherits(chunk_result, "try-error"),
        paste0(
          "species chunk ", i, " 子进程失败：",
          paste(as.character(chunk_result), collapse = " ")
        )
      )
      species_output_hash <- digest(
        chunk_result,
        algo = "sha256",
        serialize = TRUE
      )
      chunk_checkpoint <- list(
        pipeline_signature = pipeline_signature,
        species_execution_signature = species_execution_signature,
        main_sequence_hash = main_sequence_hash,
        main_taxonomy_hash = main_taxonomy_hash,
        chunk_index = as.integer(i),
        chunk_total = as.integer(length(chunks)),
        chunk_sequence_hash = chunk_sequence_hash,
        chunk_taxonomy_hash = chunk_taxonomy_hash,
        outer_chunk_size = outer_chunk_size,
        inner_n = inner_n,
        species_output_hash = species_output_hash,
        taxonomy_chunk = chunk_result,
        generated_date = as.character(Sys.Date()),
        cache_role = "reusable intermediate; not current project state"
      )
      atomic_save_rds(
        chunk_checkpoint,
        checkpoint_file,
        validate_chunk_checkpoint,
        compress = "xz"
      )
    }
    output[[i]] <- chunk_checkpoint$taxonomy_chunk
    rm(taxonomy_chunk, chunk_checkpoint)
    gc(verbose = FALSE)
  }
  output <- do.call(rbind, output)
  fail_if(
    nrow(output) != nrow(taxonomy) || anyDuplicated(rownames(output)) > 0L ||
      !setequal(rownames(output), rownames(taxonomy)),
    "分块 species 回填序列集合错误"
  )
  output <- output[rownames(taxonomy), , drop = FALSE]
  fail_if(!identical(rownames(output), rownames(taxonomy)),
          "分块 species 回填顺序错误")
  fail_if(
    !identical(colnames(output), c(colnames(taxonomy), "Species")) ||
      !identical(output[, colnames(taxonomy), drop = FALSE], taxonomy),
    "分块 species 合并后原 taxonomy 或列结构被改变"
  )
  output
}

learn_errors_checked <- function(files, nbases, seed, label) {
  set.seed(seed)
  nonconvergence_message <- FALSE
  errors <- withCallingHandlers(
    learnErrors(
      files,
      nbases = nbases,
      randomize = TRUE,
      MAX_CONSIST = 20,
      multithread = 8,
      verbose = TRUE
    ),
    message = function(message_condition) {
      if (grepl(
        "terminated before convergence",
        conditionMessage(message_condition),
        fixed = TRUE
      )) {
        nonconvergence_message <<- TRUE
      }
    }
  )
  fail_if(
    nonconvergence_message,
    paste0(label, " 误差模型在 20 次自洽迭代内仍未收敛")
  )
  errors
}

paired_clr_test <- function(genus_table, metadata, prevalence_cut = 0.20,
                            n_bootstrap = 0L, seed = 20260711L) {
  fail_if(!identical(rownames(genus_table), metadata$sample_name),
          "属表与样本元数据顺序不一致")
  prevalence <- colMeans(genus_table > 0)
  total_count <- colSums(genus_table)
  keep <- prevalence >= prevalence_cut & total_count >= 10L
  fail_if(sum(keep) < 2L, "属水平 prevalence 过滤后特征少于 2 个")
  filtered <- genus_table[, keep, drop = FALSE]
  log_counts <- log(filtered + 0.5)
  clr <- sweep(log_counts, 1L, rowMeans(log_counts), "-")
  relative <- genus_table / rowSums(genus_table)

  pair_ids <- sort(unique(metadata$patient_pair_id))
  pair_index <- lapply(pair_ids, function(pair_id) {
    tumor <- which(
      metadata$patient_pair_id == pair_id & metadata$tissue_role == "tumor"
    )
    normal <- which(
      metadata$patient_pair_id == pair_id &
        metadata$tissue_role == "paired_non_tumor"
    )
    fail_if(length(tumor) != 1L || length(normal) != 1L,
            paste("患者对不完整：", pair_id))
    c(tumor = tumor, normal = normal)
  })

  rows <- vector("list", ncol(filtered))
  for (j in seq_len(ncol(filtered))) {
    clr_differences <- vapply(pair_index, function(index) {
      clr[index[["tumor"]], j] - clr[index[["normal"]], j]
    }, FUN.VALUE = numeric(1))
    relative_differences <- vapply(pair_index, function(index) {
      relative[index[["tumor"]], colnames(filtered)[[j]]] -
        relative[index[["normal"]], colnames(filtered)[[j]]]
    }, FUN.VALUE = numeric(1))
    tumor_present <- vapply(pair_index, function(index) {
      filtered[index[["tumor"]], j] > 0
    }, FUN.VALUE = logical(1))
    normal_present <- vapply(pair_index, function(index) {
      filtered[index[["normal"]], j] > 0
    }, FUN.VALUE = logical(1))
    tumor_only <- sum(tumor_present & !normal_present)
    normal_only <- sum(!tumor_present & normal_present)
    discordant_presence <- tumor_only + normal_only
    presence_p <- if (discordant_presence > 0L) {
      binom.test(tumor_only, discordant_presence, p = 0.5)$p.value
    } else {
      1
    }
    wilcox_p <- suppressWarnings(wilcox.test(
      clr_differences, mu = 0, exact = FALSE
    )$p.value)
    bootstrap <- if (n_bootstrap > 0L) {
      bootstrap_paired_median(
        clr_differences,
        n_bootstrap = n_bootstrap,
        seed = seed + j
      )
    } else {
      list(
        valid = NA_integer_, ci_low = NA_real_, ci_high = NA_real_,
        direction_consistency = NA_real_
      )
    }
    rows[[j]] <- data.table(
      genus_label = colnames(filtered)[[j]],
      prevalence_cut = prevalence_cut,
      n_pairs = length(pair_ids),
      prevalence_all = prevalence[[colnames(filtered)[[j]]]],
      prevalence_tumor = mean(tumor_present),
      prevalence_paired_non_tumor = mean(normal_present),
      tumor_only_presence_pairs = tumor_only,
      normal_only_presence_pairs = normal_only,
      presence_exact_p = presence_p,
      median_clr_difference_tumor_minus_normal = median(clr_differences),
      mean_clr_difference_tumor_minus_normal = mean(clr_differences),
      paired_rank_biserial = paired_rank_biserial(clr_differences),
      paired_wilcoxon_p = wilcox_p,
      median_relative_abundance_difference = median(relative_differences),
      bootstrap_iterations = n_bootstrap,
      bootstrap_valid = bootstrap$valid,
      clr_effect_ci95_low = bootstrap$ci_low,
      clr_effect_ci95_high = bootstrap$ci_high,
      bootstrap_direction_consistency = bootstrap$direction_consistency
    )
  }
  result <- rbindlist(rows)
  result[, paired_wilcoxon_q := p.adjust(paired_wilcoxon_p, method = "BH")]
  result[, presence_exact_q := p.adjust(presence_exact_p, method = "BH")]
  setorder(result, paired_wilcoxon_q, paired_wilcoxon_p, genus_label)
  list(result = result, filtered_counts = filtered, clr = clr)
}

paired_alpha_test <- function(alpha_table, metric, seed) {
  tumor <- alpha_table[tissue_role == "tumor", .(
    patient_pair_id,
    tumor_value = get(metric)
  )]
  normal <- alpha_table[tissue_role == "paired_non_tumor", .(
    patient_pair_id,
    normal_value = get(metric)
  )]
  paired <- merge(tumor, normal, by = "patient_pair_id", all = FALSE)
  differences <- paired$tumor_value - paired$normal_value
  test <- suppressWarnings(wilcox.test(
    paired$tumor_value,
    paired$normal_value,
    paired = TRUE,
    exact = FALSE
  ))
  bootstrap <- bootstrap_paired_median(
    differences,
    n_bootstrap = 1000L,
    seed = seed
  )
  data.table(
    metric = metric,
    n_pairs = nrow(paired),
    median_tumor = median(paired$tumor_value),
    median_paired_non_tumor = median(paired$normal_value),
    median_paired_difference_tumor_minus_normal = median(differences),
    mean_paired_difference_tumor_minus_normal = mean(differences),
    paired_rank_biserial = paired_rank_biserial(differences),
    paired_wilcoxon_p = test$p.value,
    effect_ci95_low = bootstrap$ci_low,
    effect_ci95_high = bootstrap$ci_high,
    bootstrap_direction_consistency = bootstrap$direction_consistency
  )
}

paired_permanova <- function(distance_object, metadata, distance_name, seed) {
  fail_if(
    !identical(attr(distance_object, "Labels"), metadata$sample_name),
    paste0(distance_name, " 距离对象标签/顺序与 metadata 不一致")
  )
  permutation_design <- how(nperm = 999)
  setBlocks(permutation_design) <- metadata$patient_pair_id
  setWithin(permutation_design) <- Within(type = "free")
  set.seed(seed)
  adonis <- adonis2(
    distance_object ~ tissue_role,
    data = as.data.frame(metadata),
    permutations = permutation_design
  )
  dispersion <- betadisper(distance_object, group = metadata$tissue_role)
  set.seed(seed + 1L)
  dispersion_test <- permutest(
    dispersion,
    permutations = permutation_design
  )
  dispersion_values <- data.table(
    sample_name = metadata$sample_name,
    patient_pair_id = metadata$patient_pair_id,
    tissue_role = metadata$tissue_role,
    distance_to_centroid = dispersion$distances
  )
  tumor <- dispersion_values[tissue_role == "tumor", .(
    patient_pair_id,
    tumor_distance = distance_to_centroid
  )]
  normal <- dispersion_values[tissue_role == "paired_non_tumor", .(
    patient_pair_id,
    normal_distance = distance_to_centroid
  )]
  dispersion_pairs <- merge(tumor, normal, by = "patient_pair_id")
  paired_dispersion_p <- suppressWarnings(wilcox.test(
    dispersion_pairs$tumor_distance,
    dispersion_pairs$normal_distance,
    paired = TRUE,
    exact = FALSE
  )$p.value)
  data.table(
    distance = distance_name,
    n_samples = nrow(metadata),
    n_pairs = uniqueN(metadata$patient_pair_id),
    permanova_f = adonis$F[[1]],
    permanova_r2 = adonis$R2[[1]],
    restricted_permutation_p = adonis$`Pr(>F)`[[1]],
    permutation_scheme = "999 permutations restricted within patient_pair_id",
    permutation_seed = seed,
    betadisper_f = dispersion_test$tab$F[[1]],
    betadisper_p = dispersion_test$tab$`Pr(>F)`[[1]],
    betadisper_permutation_scheme =
      "999 tissue-label permutations restricted within patient_pair_id",
    betadisper_seed = seed + 1L,
    paired_distance_to_centroid_p = paired_dispersion_p
  )
}

bray_distance_correlation <- function(table_a, table_b) {
  common_samples <- intersect(rownames(table_a), rownames(table_b))
  all_taxa <- union(colnames(table_a), colnames(table_b))
  align <- function(x) {
    output <- matrix(
      0,
      nrow = length(common_samples),
      ncol = length(all_taxa),
      dimnames = list(common_samples, all_taxa)
    )
    output[, colnames(x)] <- x[common_samples, , drop = FALSE]
    output / rowSums(output)
  }
  distance_a <- as.numeric(vegdist(align(table_a), method = "bray"))
  distance_b <- as.numeric(vegdist(align(table_b), method = "bray"))
  suppressWarnings(cor(distance_a, distance_b, method = "spearman"))
}

fail_if(!file_exists(datasets_file), "缺少 data/datasets.tsv")
datasets <- fread(datasets_file, showProgress = FALSE)
microbiome_row <- datasets[grepl("PRJNA766558", dataset_key)]
silva_row <- datasets[grepl("SILVA.*138.2", dataset_key)]
fail_if(nrow(microbiome_row) != 1L, "PRJNA766558 数据引用必须唯一")
fail_if(nrow(silva_row) != 1L, "SILVA v138.2 数据引用必须唯一")

microbiome_root <- microbiome_row$central_path[[1]]
silva_root <- silva_row$central_path[[1]]
metadata_file <- file.path(
  microbiome_root, "10_metadata", "prjna766558_pairing.tsv"
)
fastq_qc_file <- file.path(
  microbiome_root, "10_metadata", "fastq_conversion_qc.tsv"
)
fastq_dir <- file.path(microbiome_root, "20_reusable", "fastq")
taxonomy_reference <- file.path(
  silva_root, "00_source", "silva_nr99_v138.2_toGenus_trainset.fa.gz"
)
species_reference <- file.path(
  silva_root, "00_source", "silva_v138.2_assignSpecies.fa.gz"
)
required_inputs <- c(
  metadata_file, fastq_qc_file, taxonomy_reference, species_reference
)
fail_if(any(!file_exists(required_inputs)), paste(
  "缺少 PRJNA766558 DADA2 输入：",
  paste(required_inputs[!file_exists(required_inputs)], collapse = ", ")
))

dir_create(work_intermediate_dir, recurse = TRUE)
stage_dir <- tempfile(pattern = ".prjna766558_dada2_", tmpdir = work_intermediate_dir)
cache_dir <- file.path(work_intermediate_dir, "prjna766558_dada2_cache")
sequence_checkpoint_file <- file.path(cache_dir, "sequence_checkpoint.rds")
taxonomy_checkpoint_file <- file.path(cache_dir, "taxonomy_checkpoint.rds")
species_checkpoint_dir <- file.path(cache_dir, "species_checkpoints")
filtered_dir <- file.path(stage_dir, "filtered_fastq")
dir_create(filtered_dir, recurse = TRUE)
dir_create(cache_dir, recurse = TRUE)
on.exit({
  if (dir_exists(stage_dir)) dir_delete(stage_dir)
}, add = TRUE)

message("[1/8] 校验 21 对样本映射与 84 个 L2 FASTQ")
metadata <- fread(metadata_file, showProgress = FALSE)
fastq_qc <- fread(fastq_qc_file, showProgress = FALSE)
metadata_identity_columns <- c(
  "patient_pair_id", "paper_pair_number", "tissue_role", "sample_name",
  "run_accession", "ffpe_status", "public_negative_control_run"
)
fastq_identity_columns <- c(
  "run_accession", "read_direction", "patient_pair_id", "paper_pair_number",
  "tissue_role", "sample_name", "ffpe_status", "public_negative_control_run",
  "relative_path", "file_size_bytes", "sha256",
  "read_count_matches_sra_spots", "read_count", "n_bases", "base_count",
  "primer_r1_anchored_matches", "primer_r2_anchored_matches"
)
require_columns(metadata, metadata_identity_columns, "prjna766558_pairing.tsv")
require_columns(fastq_qc, fastq_identity_columns, "fastq_conversion_qc.tsv")
fail_if(
  anyNA(metadata[, ..metadata_identity_columns]),
  "pairing metadata 的患者/样本身份字段存在 NA"
)
fail_if(
  anyNA(fastq_qc[, ..fastq_identity_columns]),
  "FASTQ QC 的映射、大小、SHA 或完整性字段存在 NA"
)
metadata[, tissue_order := match(
  tissue_role,
  c("paired_non_tumor", "tumor")
)]
setorder(metadata, paper_pair_number, tissue_order)
metadata[, tissue_order := NULL]
metadata[, tissue_role := factor(
  tissue_role,
  levels = c("paired_non_tumor", "tumor")
)]
fail_if(nrow(metadata) != 42L || uniqueN(metadata$patient_pair_id) != 21L,
        "PRJNA766558 必须为 21 对/42 个样本")
fail_if(anyDuplicated(metadata$sample_name) > 0L ||
          anyDuplicated(metadata$run_accession) > 0L,
        "sample_name 或 run_accession 重复")
pair_counts <- metadata[, .N, by = .(patient_pair_id, tissue_role)]
fail_if(nrow(pair_counts) != 42L || any(pair_counts$N != 1L),
        "患者对的肿瘤/癌旁结构不完整")
fail_if(nrow(fastq_qc) != 84L || any(!fastq_qc$read_count_matches_sra_spots),
        "84 个 FASTQ 的 read 数核查未通过")
ambiguous_base_fraction <- sum(fastq_qc$n_bases) / sum(fastq_qc$base_count)
fail_if(!is.finite(ambiguous_base_fraction) || ambiguous_base_fraction > 1e-4,
        "FASTQ 模糊碱基比例超过预设软容忍范围；需先检查来源质量")
anchored_primer_matches <- sum(
  fastq_qc$primer_r1_anchored_matches +
    fastq_qc$primer_r2_anchored_matches
)
total_fastq_reads <- sum(fastq_qc$read_count)
anchored_primer_fraction <- anchored_primer_matches / total_fastq_reads
primer_second_trim_review_threshold <- 1e-3
fail_if(
  !is.finite(anchored_primer_fraction) ||
    anchored_primer_fraction > primer_second_trim_review_threshold,
  "残余锚定引物比例超过 0.1% 软审查阈值；不能继续执行“不二次切引物”路线"
)
silva_manifest <- fread(
  file.path(silva_root, "90_manifests", "MANIFEST.tsv"),
  showProgress = FALSE
)
require_columns(
  silva_manifest,
  c("relative_path", "size_bytes", "sha256", "file_status"),
  "SILVA MANIFEST.tsv"
)
reference_relative_paths <- c(
  "00_source/silva_nr99_v138.2_toGenus_trainset.fa.gz",
  "00_source/silva_v138.2_assignSpecies.fa.gz"
)
reference_records <- silva_manifest[
  match(reference_relative_paths, relative_path),
  .(relative_path, size_bytes, sha256, file_status)
]
fail_if(
  nrow(reference_records) != 2L || anyNA(reference_records$relative_path) ||
    any(!nzchar(reference_records$sha256)) ||
    any(reference_records$file_status != "verified"),
        "SILVA manifest 缺少两个参考文件 SHA256")

# 在签名前锁定 FASTQ 字节与患者/方向映射，避免“同一批字节但
# 样本或 R1/R2 标签互换”时错误复用 checkpoint。
fail_if(
  nrow(fastq_qc) != 84L ||
    !setequal(unique(fastq_qc$read_direction), c("R1", "R2")),
  "FASTQ QC 必须包含 42 run 各一个 R1/R2"
)
fastq_direction_counts <- fastq_qc[, .N, by = .(run_accession, read_direction)]
fail_if(
  nrow(fastq_direction_counts) != 84L || any(fastq_direction_counts$N != 1L) ||
    uniqueN(fastq_qc$run_accession) != 42L,
  "FASTQ QC 中 run–read_direction 映射不是 42×2 唯一结构"
)
fastq_sample_map <- unique(fastq_qc[, .(
  run_accession, sample_name, patient_pair_id, paper_pair_number,
  tissue_role = as.character(tissue_role), ffpe_status,
  public_negative_control_run
)])
metadata_sample_map <- metadata[, .(
  run_accession, sample_name, patient_pair_id, paper_pair_number,
  tissue_role = as.character(tissue_role), ffpe_status,
  public_negative_control_run
)]
setorder(fastq_sample_map, run_accession)
setorder(metadata_sample_map, run_accession)
fail_if(
  nrow(fastq_sample_map) != 42L ||
    !identical(fastq_sample_map, metadata_sample_map),
  "FASTQ QC 与 pairing metadata 的 run/样本/患者/组织映射不一致"
)
fastq_signature_rows <- fastq_qc[, .(
  sample_name,
  run_accession,
  read_direction,
  relative_path,
  file_size_bytes,
  sha256
)]
setorder(
  fastq_signature_rows,
  sample_name, run_accession, read_direction, relative_path
)
fastq_signature_rows[, absolute_path := file.path(
  microbiome_root, relative_path
)]
fail_if(any(!file_exists(fastq_signature_rows$absolute_path)),
        "FASTQ QC 登记的 L2 文件缺失")
expected_fastq_paths <- file.path(
  fastq_dir,
  paste0(fastq_signature_rows$run_accession, "_",
         fifelse(fastq_signature_rows$read_direction == "R1", "1", "2"),
         ".fastq.gz")
)
fail_if(
  !identical(
    normalizePath(fastq_signature_rows$absolute_path, winslash = "/"),
    normalizePath(expected_fastq_paths, winslash = "/")
  ),
  "FASTQ relative_path 与 run/read_direction 期望路径不一致"
)
actual_fastq_size <- as.numeric(file_info(fastq_signature_rows$absolute_path)$size)
actual_fastq_sha <- vapply(
  fastq_signature_rows$absolute_path,
  digest,
  FUN.VALUE = character(1),
  algo = "sha256",
  file = TRUE,
  serialize = FALSE
)
fail_if(
  any(actual_fastq_size != as.numeric(fastq_signature_rows$file_size_bytes)) ||
    any(actual_fastq_sha != fastq_signature_rows$sha256),
  "84 个 FASTQ 的当前字节大小/SHA256 与冻结 QC 不一致"
)
reference_records[, absolute_path := file.path(silva_root, relative_path)]
fail_if(any(!file_exists(reference_records$absolute_path)),
        "SILVA 参考实体文件缺失")
actual_reference_size <- as.numeric(file_info(reference_records$absolute_path)$size)
actual_reference_sha <- vapply(
  reference_records$absolute_path,
  digest,
  FUN.VALUE = character(1),
  algo = "sha256",
  file = TRUE,
  serialize = FALSE
)
fail_if(
  any(actual_reference_size != as.numeric(reference_records$size_bytes)) ||
    any(actual_reference_sha != reference_records$sha256),
  "SILVA 参考当前字节大小/SHA256 与 manifest 不一致"
)
signature_components <- c(
  fastq_signature_rows[, paste(
    "fastq", paste0("sample=", sample_name),
    paste0("run=", run_accession), paste0("direction=", read_direction),
    paste0("relative_path=", relative_path),
    paste0("size=", file_size_bytes), paste0("sha256=", sha256),
    sep = "|"
  )],
  reference_records[, paste(
    "reference", paste0("relative_path=", relative_path),
    paste0("size=", size_bytes), paste0("sha256=", sha256), sep = "|"
  )],
  metadata_sample_map[, paste(
    "metadata", paste0("sample=", sample_name),
    paste0("run=", run_accession), paste0("pair=", patient_pair_id),
    paste0("paper_pair_number=", paper_pair_number),
    paste0("tissue_role=", tissue_role), paste0("ffpe_status=", ffpe_status),
    paste0("public_negative_control_run=", public_negative_control_run),
    sep = "|"
  )],
  paste0("pipeline_schema_version=", pipeline_schema_version),
  paste0("dada2=", packageVersion("dada2")),
  "truncLen=0/0;trimLeft=0/0;maxN=0;maxEE=2/2;truncQ=2;minLen=100;rm.phix=TRUE",
  paste0(
    "primer_second_trim=none;anchored_fraction=",
    format(anchored_primer_fraction, scientific = TRUE, digits = 12),
    ";review_threshold=", primer_second_trim_review_threshold
  ),
  "error_nbases=5e8;randomize=TRUE;MAX_CONSIST=20;pool=FALSE",
  "merge_main=12;merge_sensitivity=8;maxMismatch=0;chimera=consensus",
  "taxonomy=SILVA138.2;minBoot=50;tryRC=TRUE;chunk_size=1000",
  # 下面的 n=2000 是已冻结上游 sequence/taxonomy checkpoint 的历史
  # signature component。DADA2 的 n 只控制 species 精确匹配的计算分块，
  # 不改变匹配定义；当前实际分块/worker 参数另由 species_execution_signature
  # 和正式 manifest 字段精确记录，避免为内存实现修复重跑去噪与属级分类。
  "species=SILVA138.2;allowMultiple=TRUE;tryRC=TRUE;n=2000"
)
# digest(..., serialize = FALSE) 对 character vector 只会消费首个元素，
# 因此必须先折叠为唯一标量 payload，确保 84 个 FASTQ SHA、
# 两个参考 SHA 和全部锁定参数都进入检查点签名。
signature_payload <- paste(signature_components, collapse = "|")
pipeline_signature <- digest(
  signature_payload,
  algo = "sha256",
  serialize = FALSE
)
species_outer_chunk_size <- 1000L
species_inner_n <- 250L
species_execution_components <- c(
  paste0("pipeline_signature=", pipeline_signature),
  paste0("species_execution_schema_version=", species_execution_schema_version),
  paste0("outer_chunk_size=", species_outer_chunk_size),
  paste0("inner_n=", species_inner_n),
  "scheduler=length_group_preserving_binpack",
  "worker=parallel_mcparallel_sequential",
  "checkpoint=signature_partitioned_per_chunk_atomic_xz_with_output_hash",
  "result_semantics=addSpecies_allowMultiple_TRUE_tryRC_TRUE"
)
species_execution_signature <- digest(
  paste(species_execution_components, collapse = "|"),
  algo = "sha256",
  serialize = FALSE
)
dir_create(species_checkpoint_dir, recurse = TRUE)
species_checkpoint_file <- file.path(
  species_checkpoint_dir,
  paste0("species_checkpoint_", species_execution_signature, ".rds")
)

metadata[, forward_fastq := file.path(
  fastq_dir, paste0(run_accession, "_1.fastq.gz")
)]
metadata[, reverse_fastq := file.path(
  fastq_dir, paste0(run_accession, "_2.fastq.gz")
)]
fail_if(any(!file_exists(metadata$forward_fastq)) ||
          any(!file_exists(metadata$reverse_fastq)),
        "部分 L2 FASTQ 不可读取")
metadata[, filtered_forward := file.path(
  filtered_dir, paste0(sample_name, "_F_filt.fastq.gz")
)]
metadata[, filtered_reverse := file.path(
  filtered_dir, paste0(sample_name, "_R_filt.fastq.gz")
)]

validate_count_vector <- function(value, label) {
  fail_if(!is.numeric(value), paste0(label, " 不是数值向量"))
  fail_if(
    !identical(names(value), metadata$sample_name),
    paste0(label, " 的样本 names/顺序与 metadata 不一致")
  )
  fail_if(
    any(!is.finite(value)) || any(value < 0),
    paste0(label, " 存在 NA/非有限或负值")
  )
  invisible(TRUE)
}

validate_sequence_table <- function(value, label) {
  fail_if(!is.matrix(value), paste0(label, " 不是 matrix"))
  fail_if(
    !identical(rownames(value), metadata$sample_name),
    paste0(label, " 样本行/顺序与 metadata 不一致")
  )
  fail_if(
    ncol(value) < 1L || is.null(colnames(value)) ||
      anyDuplicated(colnames(value)) > 0L || any(!nzchar(colnames(value))),
    paste0(label, " 序列列名缺失或重复")
  )
  fail_if(
    !is.numeric(value) || any(!is.finite(value)) || any(value < 0),
    paste0(label, " 计数存在 NA/非有限或负值")
  )
  invisible(TRUE)
}

validate_sequence_checkpoint <- function(checkpoint) {
  required_fields <- c(
    "pipeline_signature", "filter_tracking", "denoised_forward_counts",
    "denoised_reverse_counts", "merged_main_counts",
    "merged_overlap8_counts", "sequence_main_nochim",
    "sequence_overlap8_nochim", "sequence_forward_nochim"
  )
  require_list_fields(checkpoint, required_fields, "sequence checkpoint")
  fail_if(
    !identical(checkpoint$pipeline_signature, pipeline_signature),
    "sequence checkpoint 与当前输入/参数签名不一致"
  )
  tracking <- checkpoint$filter_tracking
  fail_if(!is.matrix(tracking), "sequence checkpoint filter_tracking 不是 matrix")
  fail_if(
    !identical(rownames(tracking), metadata$sample_name) ||
      !all(c("reads.in", "reads.out") %in% colnames(tracking)),
    "sequence checkpoint filter_tracking 样本/字段不完整"
  )
  fail_if(
    any(!is.finite(tracking[, c("reads.in", "reads.out"), drop = FALSE])) ||
      any(tracking[, c("reads.in", "reads.out"), drop = FALSE] < 0),
    "sequence checkpoint filter_tracking 存在 NA/非有限或负值"
  )
  validate_count_vector(checkpoint$denoised_forward_counts,
                        "denoised_forward_counts")
  validate_count_vector(checkpoint$denoised_reverse_counts,
                        "denoised_reverse_counts")
  validate_count_vector(checkpoint$merged_main_counts,
                        "merged_main_counts")
  validate_count_vector(checkpoint$merged_overlap8_counts,
                        "merged_overlap8_counts")
  validate_sequence_table(checkpoint$sequence_main_nochim,
                          "sequence_main_nochim")
  validate_sequence_table(checkpoint$sequence_overlap8_nochim,
                          "sequence_overlap8_nochim")
  validate_sequence_table(checkpoint$sequence_forward_nochim,
                          "sequence_forward_nochim")
  reads_in <- tracking[, "reads.in"]
  reads_out <- tracking[, "reads.out"]
  fail_if(any(reads_out > reads_in),
          "sequence checkpoint 出现 filtered reads > input reads")
  fail_if(
    any(checkpoint$denoised_forward_counts > reads_out) ||
      any(checkpoint$denoised_reverse_counts > reads_out),
    "sequence checkpoint 出现 denoised reads > filtered reads"
  )
  fail_if(
    any(checkpoint$merged_main_counts > checkpoint$denoised_forward_counts) ||
      any(checkpoint$merged_main_counts > checkpoint$denoised_reverse_counts) ||
      any(checkpoint$merged_overlap8_counts > checkpoint$denoised_forward_counts) ||
      any(checkpoint$merged_overlap8_counts > checkpoint$denoised_reverse_counts) ||
      any(checkpoint$merged_main_counts > checkpoint$merged_overlap8_counts),
    "sequence checkpoint 合并阶段 reads 账目关系异常"
  )
  fail_if(
    any(rowSums(checkpoint$sequence_main_nochim) >
          checkpoint$merged_main_counts) ||
      any(rowSums(checkpoint$sequence_overlap8_nochim) >
          checkpoint$merged_overlap8_counts) ||
      any(rowSums(checkpoint$sequence_forward_nochim) >
          checkpoint$denoised_forward_counts),
    "sequence checkpoint 出现 nonchim reads 超过上一阶段"
  )
  invisible(TRUE)
}

if (file_exists(sequence_checkpoint_file)) {
  message("[2/8] 读取已验证的去嵌合序列表 checkpoint")
  sequence_checkpoint <- readRDS(sequence_checkpoint_file)
  validate_sequence_checkpoint(sequence_checkpoint)
  filter_tracking <- sequence_checkpoint$filter_tracking
  denoised_forward_counts <- sequence_checkpoint$denoised_forward_counts
  denoised_reverse_counts <- sequence_checkpoint$denoised_reverse_counts
  merged_main_counts <- sequence_checkpoint$merged_main_counts
  merged_overlap8_counts <- sequence_checkpoint$merged_overlap8_counts
  sequence_main_nochim <- sequence_checkpoint$sequence_main_nochim
  sequence_overlap8_nochim <- sequence_checkpoint$sequence_overlap8_nochim
  sequence_forward_nochim <- sequence_checkpoint$sequence_forward_nochim
} else {
  message("[2/8] expected-error 过滤：不切引物、不预截断、maxEE=2/2")
  filter_tracking <- filterAndTrim(
  metadata$forward_fastq,
  metadata$filtered_forward,
  metadata$reverse_fastq,
  metadata$filtered_reverse,
  truncLen = c(0, 0),
  trimLeft = c(0, 0),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  minLen = 100,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = 8,
  verbose = TRUE
  )
  rownames(filter_tracking) <- metadata$sample_name
  fail_if(any(filter_tracking[, "reads.out"] == 0L),
          "至少一个样本过滤后无 reads")

  message("[3/8] 学习误差、去噪并建立主/敏感性序列表")
  err_forward <- learn_errors_checked(
  metadata$filtered_forward,
  nbases = 5e8,
  seed = 20260731L,
  label = "forward"
  )
  err_reverse <- learn_errors_checked(
  metadata$filtered_reverse,
  nbases = 5e8,
  seed = 20260732L,
  label = "reverse"
  )
  derep_forward <- derepFastq(metadata$filtered_forward, verbose = FALSE)
  derep_reverse <- derepFastq(metadata$filtered_reverse, verbose = FALSE)
  names(derep_forward) <- metadata$sample_name
  names(derep_reverse) <- metadata$sample_name
  dada_forward <- dada(
  derep_forward, err = err_forward, pool = FALSE,
  multithread = 8, verbose = TRUE
  )
  dada_reverse <- dada(
  derep_reverse, err = err_reverse, pool = FALSE,
  multithread = 8, verbose = TRUE
  )

  merge_main <- mergePairs(
  dada_forward, derep_forward, dada_reverse, derep_reverse,
  minOverlap = 12,
  maxMismatch = 0,
  returnRejects = FALSE,
  verbose = TRUE
  )
  merge_overlap8 <- mergePairs(
  dada_forward, derep_forward, dada_reverse, derep_reverse,
  minOverlap = 8,
  maxMismatch = 0,
  returnRejects = FALSE,
  verbose = TRUE
  )
  sequence_main <- makeSequenceTable(merge_main)
  sequence_overlap8 <- makeSequenceTable(merge_overlap8)
  sequence_forward <- makeSequenceTable(dada_forward)
  sequence_main_nochim <- removeBimeraDenovo(
  sequence_main, method = "consensus", multithread = 8, verbose = TRUE
  )
  sequence_overlap8_nochim <- removeBimeraDenovo(
  sequence_overlap8, method = "consensus", multithread = 8, verbose = TRUE
  )
  sequence_forward_nochim <- removeBimeraDenovo(
  sequence_forward, method = "consensus", multithread = 8, verbose = TRUE
  )
  fail_if(ncol(sequence_main_nochim) == 0L, "主流程去嵌合后无 ASV")

  denoised_forward_counts <- vapply(dada_forward, get_n, numeric(1))
  denoised_reverse_counts <- vapply(dada_reverse, get_n, numeric(1))
  merged_main_counts <- vapply(
    merge_main, function(x) sum(x$abundance), numeric(1)
  )
  merged_overlap8_counts <- vapply(
    merge_overlap8, function(x) sum(x$abundance), numeric(1)
  )
  sequence_checkpoint <- list(
    pipeline_signature = pipeline_signature,
    filter_tracking = filter_tracking,
    denoised_forward_counts = denoised_forward_counts,
    denoised_reverse_counts = denoised_reverse_counts,
    merged_main_counts = merged_main_counts,
    merged_overlap8_counts = merged_overlap8_counts,
    sequence_main_nochim = sequence_main_nochim,
    sequence_overlap8_nochim = sequence_overlap8_nochim,
    sequence_forward_nochim = sequence_forward_nochim,
    generated_date = as.character(Sys.Date()),
    cache_role = "reusable intermediate; not current project state"
  )
  atomic_save_rds(
    sequence_checkpoint,
    sequence_checkpoint_file,
    validate_sequence_checkpoint,
    compress = "xz"
  )
  rm(
    err_forward, err_reverse, derep_forward, derep_reverse,
    dada_forward, dada_reverse, merge_main, merge_overlap8,
    sequence_main, sequence_overlap8, sequence_forward
  )
  gc(verbose = FALSE)
}

message("[4/8] SILVA v138.2 分块分类并去除非细菌/古菌及细胞器序列")
taxonomy_sequences <- unique(c(
  colnames(sequence_main_nochim),
  colnames(sequence_overlap8_nochim),
  colnames(sequence_forward_nochim)
))
taxonomy_sequence_hash <- digest(
  paste(taxonomy_sequences, collapse = "|"),
  algo = "sha256",
  serialize = FALSE
)
validate_taxonomy_checkpoint <- function(checkpoint) {
  require_list_fields(
    checkpoint,
    c("pipeline_signature", "taxonomy_sequence_hash", "taxonomy_union"),
    "taxonomy checkpoint"
  )
  fail_if(
    !identical(checkpoint$pipeline_signature, pipeline_signature) ||
      !identical(checkpoint$taxonomy_sequence_hash, taxonomy_sequence_hash),
    "taxonomy checkpoint 与当前序列表签名不一致"
  )
  taxonomy_value <- checkpoint$taxonomy_union
  fail_if(!is.matrix(taxonomy_value), "taxonomy checkpoint 分类表不是 matrix")
  fail_if(
    !identical(rownames(taxonomy_value), taxonomy_sequences) ||
      anyDuplicated(rownames(taxonomy_value)) > 0L,
    "taxonomy checkpoint 序列行/顺序与当前序列集不一致"
  )
  fail_if(
    !all(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus") %in%
           colnames(taxonomy_value)),
    "taxonomy checkpoint 缺少必要分类层级"
  )
  invisible(TRUE)
}
if (file_exists(taxonomy_checkpoint_file)) {
  taxonomy_checkpoint <- readRDS(taxonomy_checkpoint_file)
  validate_taxonomy_checkpoint(taxonomy_checkpoint)
  taxonomy_union <- taxonomy_checkpoint$taxonomy_union
  message("  已读取 taxonomy checkpoint：", nrow(taxonomy_union), " 条序列")
} else {
  taxonomy_union <- assign_taxonomy_chunked(
    taxonomy_sequences,
    taxonomy_reference,
    chunk_size = 1000L,
    threads = 8L
  )
  taxonomy_checkpoint <- list(
    pipeline_signature = pipeline_signature,
    taxonomy_sequence_hash = taxonomy_sequence_hash,
    taxonomy_union = taxonomy_union,
    generated_date = as.character(Sys.Date()),
    cache_role = "reusable intermediate; not current project state"
  )
  atomic_save_rds(
    taxonomy_checkpoint,
    taxonomy_checkpoint_file,
    validate_taxonomy_checkpoint,
    compress = "xz"
  )
}
taxonomy_main <- taxonomy_union[
  colnames(sequence_main_nochim), , drop = FALSE
]
main_sequence_hash <- digest(
  paste(rownames(taxonomy_main), collapse = "|"),
  algo = "sha256",
  serialize = FALSE
)
main_taxonomy_hash <- digest(
  taxonomy_main,
  algo = "sha256",
  serialize = TRUE
)
taxonomy_main_without_species <- taxonomy_main
validate_species_checkpoint <- function(checkpoint) {
  require_list_fields(
    checkpoint,
    c(
      "pipeline_signature", "species_execution_signature",
      "main_sequence_hash", "main_taxonomy_hash", "species_output_hash",
      "taxonomy_main"
    ),
    "species checkpoint"
  )
  fail_if(
    !identical(checkpoint$pipeline_signature, pipeline_signature) ||
      !identical(
        checkpoint$species_execution_signature,
        species_execution_signature
      ) ||
      !identical(checkpoint$main_sequence_hash, main_sequence_hash) ||
      !identical(checkpoint$main_taxonomy_hash, main_taxonomy_hash),
    "species checkpoint 与当前主序列表签名不一致"
  )
  taxonomy_value <- checkpoint$taxonomy_main
  fail_if(!is.matrix(taxonomy_value), "species checkpoint 分类表不是 matrix")
  fail_if(
    !identical(rownames(taxonomy_value), colnames(sequence_main_nochim)) ||
      anyDuplicated(rownames(taxonomy_value)) > 0L,
    "species checkpoint 序列行/顺序与主序列表不一致"
  )
  fail_if(
    !identical(
      colnames(taxonomy_value),
      c(colnames(taxonomy_main_without_species), "Species")
    ) || !identical(
      taxonomy_value[, colnames(taxonomy_main_without_species), drop = FALSE],
      taxonomy_main_without_species
    ),
    "species checkpoint 分类层级或原 taxonomy 被改变"
  )
  fail_if(
    !identical(
      checkpoint$species_output_hash,
      digest(taxonomy_value, algo = "sha256", serialize = TRUE)
    ) || any(
      !is.na(taxonomy_value[, "Species"]) &
        !nzchar(taxonomy_value[, "Species"])
    ),
    "species checkpoint 输出 hash 或 Species 值异常"
  )
  invisible(TRUE)
}
if (file_exists(species_checkpoint_file)) {
  species_checkpoint <- readRDS(species_checkpoint_file)
  validate_species_checkpoint(species_checkpoint)
  taxonomy_main <- species_checkpoint$taxonomy_main
  message("  已读取 species checkpoint：", nrow(taxonomy_main), " 条主序列")
} else {
  # assignSpecies 会为每条 query 暂存相对于完整 species 参考的命中向量；
  # 18,337 条主 ASV 一次执行会超过本机 24 GiB R vector heap 上限。
  # 外层 1,000 条分块并逐块在独立子进程中运行，使每块结束后释放参考与
  # 命中矩阵；每块均绑定输入 taxonomy/hash 并原子保存，可安全断点续跑。
  taxonomy_main <- add_species_chunked_isolated(
    taxonomy_main,
    reference = species_reference,
    cache_directory = cache_dir,
    pipeline_signature = pipeline_signature,
    main_sequence_hash = main_sequence_hash,
    main_taxonomy_hash = main_taxonomy_hash,
    species_execution_signature = species_execution_signature,
    outer_chunk_size = species_outer_chunk_size,
    inner_n = species_inner_n
  )
  species_output_hash <- digest(
    taxonomy_main,
    algo = "sha256",
    serialize = TRUE
  )
  species_checkpoint <- list(
    pipeline_signature = pipeline_signature,
    species_execution_signature = species_execution_signature,
    main_sequence_hash = main_sequence_hash,
    main_taxonomy_hash = main_taxonomy_hash,
    species_output_hash = species_output_hash,
    taxonomy_main = taxonomy_main,
    generated_date = as.character(Sys.Date()),
    cache_role = "reusable intermediate; not current project state"
  )
  atomic_save_rds(
    species_checkpoint,
    species_checkpoint_file,
    validate_species_checkpoint,
    compress = "xz"
  )
}
taxonomy_overlap8 <- taxonomy_union[
  colnames(sequence_overlap8_nochim), , drop = FALSE
]
taxonomy_forward <- taxonomy_union[
  colnames(sequence_forward_nochim), , drop = FALSE
]

filtered_main <- filter_prokaryotic_taxa(
  sequence_main_nochim, taxonomy_main
)
filtered_overlap8 <- filter_prokaryotic_taxa(
  sequence_overlap8_nochim, taxonomy_overlap8
)
filtered_forward <- filter_prokaryotic_taxa(
  sequence_forward_nochim, taxonomy_forward
)
sequence_main_final <- filtered_main$sequence_table
taxonomy_main_final <- filtered_main$taxonomy
sequence_overlap8_final <- filtered_overlap8$sequence_table
taxonomy_overlap8_final <- filtered_overlap8$taxonomy
sequence_forward_final <- filtered_forward$sequence_table
taxonomy_forward_final <- filtered_forward$taxonomy
fail_if(
  any(rowSums(sequence_main_final) > rowSums(sequence_main_nochim)) ||
    any(rowSums(sequence_overlap8_final) >
          rowSums(sequence_overlap8_nochim)) ||
    any(rowSums(sequence_forward_final) > rowSums(sequence_forward_nochim)),
  "分类过滤后的 prokaryotic reads 超过 nonchim reads"
)
fail_if(
  any(rowSums(sequence_main_final) == 0L) ||
    any(rowSums(sequence_overlap8_final) == 0L) ||
    any(rowSums(sequence_forward_final) == 0L),
  "至少一个样本在主流程或敏感性分类过滤后无 reads"
)

# 336F–806R V3–V4 的理论长度会随物种及是否保留引物而波动。
# 使用 380–480 nt 作为宽松的“合理长度”敏感性子集，异常长度
# ASV 仍保留在主结果，不机械删除。
plausible_v3v4_length_min <- 380L
plausible_v3v4_length_max <- 480L
main_asv_length <- nchar(colnames(sequence_main_final))
plausible_length_keep <-
  main_asv_length >= plausible_v3v4_length_min &
  main_asv_length <= plausible_v3v4_length_max
fail_if(!any(plausible_length_keep),
        "主流程没有 380–480 nt 的合理 V3–V4 长度 ASV")
sequence_main_plausible_length <- sequence_main_final[
  , plausible_length_keep, drop = FALSE
]
taxonomy_main_plausible_length <- taxonomy_main_final[
  colnames(sequence_main_plausible_length), , drop = FALSE
]
fail_if(
  any(rowSums(sequence_main_plausible_length) == 0L),
  "至少一个样本在 380–480 nt 合理长度敏感性子集中无 reads"
)

genus_main <- aggregate_genus(sequence_main_final, taxonomy_main_final)
genus_plausible_length <- aggregate_genus(
  sequence_main_plausible_length,
  taxonomy_main_plausible_length
)
genus_overlap8 <- aggregate_genus(
  sequence_overlap8_final, taxonomy_overlap8_final
)
genus_forward <- aggregate_genus(
  sequence_forward_final, taxonomy_forward_final
)
metadata <- metadata[match(rownames(sequence_main_final), sample_name)]
fail_if(anyNA(metadata$sample_name), "序列表样本无法回填元数据")

message("[5/8] α/β 多样性与患者对内限制置换")
sample_qc <- metadata[, .(
  patient_pair_id,
  paper_pair_number,
  sample_name,
  run_accession,
  tissue_role = as.character(tissue_role),
  ffpe_status,
  public_negative_control_run
)]
sample_qc[, `:=`(
  input_reads = filter_tracking[sample_name, "reads.in"],
  filtered_reads = filter_tracking[sample_name, "reads.out"],
  denoised_forward_reads = denoised_forward_counts[sample_name],
  denoised_reverse_reads = denoised_reverse_counts[sample_name],
  merged_reads_min_overlap12 = merged_main_counts[sample_name],
  merged_reads_min_overlap8 = merged_overlap8_counts[sample_name],
  nonchim_reads_min_overlap12 = rowSums(sequence_main_nochim)[sample_name],
  prokaryotic_reads_min_overlap12 = rowSums(sequence_main_final)[sample_name],
  prokaryotic_reads_plausible_v3v4_length =
    rowSums(sequence_main_plausible_length)[sample_name],
  prokaryotic_reads_min_overlap8 = rowSums(sequence_overlap8_final)[sample_name],
  prokaryotic_reads_forward_only = rowSums(sequence_forward_final)[sample_name]
)]
sample_qc[, `:=`(
  filter_retention = filtered_reads / input_reads,
  main_final_retention = prokaryotic_reads_min_overlap12 / input_reads,
  plausible_v3v4_length_read_fraction =
    prokaryotic_reads_plausible_v3v4_length /
      prokaryotic_reads_min_overlap12,
  overlap8_final_retention = prokaryotic_reads_min_overlap8 / input_reads,
  forward_only_final_retention = prokaryotic_reads_forward_only / input_reads,
  primer_handling = "no_second_primer_cut_source_reads_already_trimmed_likely",
  anchored_primer_matches_all_fastq = anchored_primer_matches,
  anchored_primer_fraction_all_fastq = anchored_primer_fraction,
  primer_second_trim_review_threshold = primer_second_trim_review_threshold,
  plausible_v3v4_length_range = paste0(
    plausible_v3v4_length_min, "-", plausible_v3v4_length_max, " nt"
  ),
  main_filter = "truncLen=0/0;maxN=0;maxEE=2/2;truncQ=2;minLen=100",
  blank_control_status = "no_public_negative_control_or_extraction_blank"
)]

common_rarefaction_depth <- min(rowSums(sequence_main_final))
fail_if(
  !is.finite(common_rarefaction_depth) || common_rarefaction_depth < 1L,
  "主流程全样本稀释深度不可用"
)
rarefaction_repeats <- 100L
set.seed(20260726L)
rarefied_metrics <- lapply(seq_len(rarefaction_repeats), function(iteration) {
  rarefied <- suppressWarnings(rrarefy(
    sequence_main_final,
    sample = common_rarefaction_depth
  ))
  data.table(
    sample_name = rownames(rarefied),
    observed_asv = rowSums(rarefied > 0),
    shannon = diversity(rarefied, index = "shannon"),
    simpson = diversity(rarefied, index = "simpson")
  )
})
rarefied_long <- rbindlist(
  lapply(seq_along(rarefied_metrics), function(iteration) {
    output <- copy(rarefied_metrics[[iteration]])
    output[, rarefaction_iteration := iteration]
    output
  })
)
rarefied_summary <- rarefied_long[, .(
  rarefied_observed_asv_median = median(observed_asv),
  rarefied_shannon_median = median(shannon),
  rarefied_simpson_median = median(simpson)
), by = sample_name]

alpha_diversity <- data.table(
  sample_name = rownames(sequence_main_final),
  observed_asv = rowSums(sequence_main_final > 0),
  shannon = diversity(sequence_main_final, index = "shannon"),
  simpson = diversity(sequence_main_final, index = "simpson")
)
alpha_diversity <- merge(
  sample_qc[, .(
    sample_name, patient_pair_id, paper_pair_number, tissue_role,
    prokaryotic_reads_min_overlap12
  )],
  alpha_diversity,
  by = "sample_name",
  sort = FALSE
)
alpha_diversity <- merge(
  alpha_diversity,
  rarefied_summary,
  by = "sample_name",
  all.x = TRUE,
  sort = FALSE
)
alpha_diversity[, `:=`(
  rarefaction_depth = common_rarefaction_depth,
  rarefaction_repeats = rarefaction_repeats,
  rarefaction_seed = 20260726L
)]
alpha_diversity[, tissue_order := match(
  tissue_role,
  c("paired_non_tumor", "tumor")
)]
setorder(alpha_diversity, paper_pair_number, tissue_order)
alpha_diversity[, tissue_order := NULL]
alpha_tests <- rbindlist(list(
  paired_alpha_test(
    alpha_diversity, "prokaryotic_reads_min_overlap12", 20260720L
  ),
  paired_alpha_test(alpha_diversity, "observed_asv", 20260721L),
  paired_alpha_test(alpha_diversity, "shannon", 20260722L),
  paired_alpha_test(alpha_diversity, "simpson", 20260723L),
  paired_alpha_test(
    alpha_diversity, "rarefied_observed_asv_median", 20260727L
  ),
  paired_alpha_test(
    alpha_diversity, "rarefied_shannon_median", 20260728L
  ),
  paired_alpha_test(
    alpha_diversity, "rarefied_simpson_median", 20260729L
  )
))
alpha_tests[, paired_wilcoxon_q := p.adjust(
  paired_wilcoxon_p, method = "BH"
)]
alpha_tests[, interpretation_scope := fcase(
  metric == "prokaryotic_reads_min_overlap12",
  "paired_final_depth_diagnostic_not_diversity",
  grepl("^rarefied_", metric),
  "100_repeat_common_depth_rarefaction_sensitivity",
  default = "raw_count_table_depth_sensitive_description"
)]

relative_main <- sequence_main_final / rowSums(sequence_main_final)
bray_distance <- vegdist(relative_main, method = "bray")
aitchison_keep <- colMeans(sequence_main_final > 0) >= 0.10 &
  colSums(sequence_main_final) >= 10L
fail_if(sum(aitchison_keep) < 2L,
        "Aitchison 距离过滤后 ASV 少于 2 个")
aitchison_log <- log(sequence_main_final[, aitchison_keep, drop = FALSE] + 0.5)
aitchison_clr <- sweep(aitchison_log, 1L, rowMeans(aitchison_log), "-")
aitchison_distance <- dist(aitchison_clr, method = "euclidean")
metadata_beta <- metadata[, .(
  sample_name,
  patient_pair_id,
  tissue_role = factor(tissue_role,
                       levels = c("paired_non_tumor", "tumor"))
)]
beta_tests <- rbindlist(list(
  paired_permanova(
    bray_distance, metadata_beta,
    "Bray-Curtis ASV relative abundance", 20260724L
  ),
  paired_permanova(
    aitchison_distance, metadata_beta,
    "Aitchison ASV CLR", 20260725L
  )
))
beta_tests[, restricted_permutation_q := p.adjust(
  restricted_permutation_p, method = "BH"
)]

message("[6/8] 属水平配对 CLR、presence exact 与流程敏感性")
genus_primary <- paired_clr_test(
  genus_main, metadata,
  prevalence_cut = 0.20,
  n_bootstrap = 1000L,
  seed = 20260800L
)
genus_main_p10 <- paired_clr_test(
  genus_main, metadata, prevalence_cut = 0.10, n_bootstrap = 0L
)$result
genus_main_p30 <- paired_clr_test(
  genus_main, metadata, prevalence_cut = 0.30, n_bootstrap = 0L
)$result
genus_overlap8_p20 <- paired_clr_test(
  genus_overlap8, metadata, prevalence_cut = 0.20, n_bootstrap = 0L
)$result
genus_forward_p20 <- paired_clr_test(
  genus_forward, metadata, prevalence_cut = 0.20, n_bootstrap = 0L
)$result
genus_plausible_length_p20 <- paired_clr_test(
  genus_plausible_length, metadata,
  prevalence_cut = 0.20, n_bootstrap = 0L
)$result

genus_differential <- copy(genus_primary$result)
genus_differential <- merge(
  genus_differential,
  genus_main_p10[, .(
    genus_label,
    clr_effect_prevalence10 = median_clr_difference_tumor_minus_normal
  )],
  by = "genus_label", all.x = TRUE, sort = FALSE
)
genus_differential <- merge(
  genus_differential,
  genus_main_p30[, .(
    genus_label,
    clr_effect_prevalence30 = median_clr_difference_tumor_minus_normal
  )],
  by = "genus_label", all.x = TRUE, sort = FALSE
)
genus_differential <- merge(
  genus_differential,
  genus_overlap8_p20[, .(
    genus_label,
    clr_effect_min_overlap8 = median_clr_difference_tumor_minus_normal
  )],
  by = "genus_label", all.x = TRUE, sort = FALSE
)
genus_differential <- merge(
  genus_differential,
  genus_forward_p20[, .(
    genus_label,
    clr_effect_forward_only = median_clr_difference_tumor_minus_normal
  )],
  by = "genus_label", all.x = TRUE, sort = FALSE
)
genus_differential <- merge(
  genus_differential,
  genus_plausible_length_p20[, .(
    genus_label,
    clr_effect_plausible_v3v4_length =
      median_clr_difference_tumor_minus_normal
  )],
  by = "genus_label", all.x = TRUE, sort = FALSE
)
genus_differential[, prevalence_threshold_direction_stable :=
  is.finite(clr_effect_prevalence10) & is.finite(clr_effect_prevalence30) &
    sign(median_clr_difference_tumor_minus_normal) ==
      sign(clr_effect_prevalence10) &
    sign(median_clr_difference_tumor_minus_normal) ==
      sign(clr_effect_prevalence30)]
genus_differential[, pipeline_direction_stable :=
  is.finite(clr_effect_min_overlap8) & is.finite(clr_effect_forward_only) &
    sign(median_clr_difference_tumor_minus_normal) ==
      sign(clr_effect_min_overlap8) &
    sign(median_clr_difference_tumor_minus_normal) ==
      sign(clr_effect_forward_only)]
genus_differential[, min_overlap8_direction_concordant :=
  is.finite(clr_effect_min_overlap8) &
    sign(median_clr_difference_tumor_minus_normal) ==
      sign(clr_effect_min_overlap8)]
genus_differential[, forward_only_direction_concordant :=
  is.finite(clr_effect_forward_only) &
    sign(median_clr_difference_tumor_minus_normal) ==
      sign(clr_effect_forward_only)]
genus_differential[, plausible_length_direction_concordant :=
  is.finite(clr_effect_plausible_v3v4_length) &
    sign(median_clr_difference_tumor_minus_normal) ==
      sign(clr_effect_plausible_v3v4_length)]
genus_differential[, effect_direction := fifelse(
  median_clr_difference_tumor_minus_normal > 0,
  "higher_in_tumor",
  fifelse(
    median_clr_difference_tumor_minus_normal < 0,
    "lower_in_tumor",
    "no_direction"
  )
)]
genus_differential[, candidate_tier := fifelse(
  paired_wilcoxon_q <= 0.10 &
    abs(median_clr_difference_tumor_minus_normal) >= 0.30 &
    pipeline_direction_stable & prevalence_threshold_direction_stable &
    plausible_length_direction_concordant &
    bootstrap_direction_consistency >= 0.80,
  "paired_clr_fdr_supported_no_blank",
  fifelse(
    (paired_wilcoxon_q <= 0.20 | paired_wilcoxon_p <= 0.05) &
      abs(median_clr_difference_tumor_minus_normal) >= 0.20 &
      pipeline_direction_stable & bootstrap_direction_consistency >= 0.70,
    "paired_clr_conditional_no_blank",
    fifelse(
      (presence_exact_q <= 0.20 | presence_exact_p <= 0.05) &
        (tumor_only_presence_pairs + normal_only_presence_pairs) >= 5L,
      "prevalence_shift_conditional_no_blank",
      fifelse(
        paired_wilcoxon_p <= 0.10 &
          abs(median_clr_difference_tumor_minus_normal) >= 0.20,
        "directional_exploratory_no_blank",
        "background_no_clear_paired_difference"
      )
    )
  )
)]
genus_differential[, `:=`(
  blank_control_status = "no_public_negative_control_or_extraction_blank",
  contamination_risk = "unresolved_low_biomass_FFPE_reagent_or_oral_carryover",
  decontam_status = "not_run_no_blank_and_no_DNA_concentration",
  function_claim_scope = "taxonomic_ecology_only_16S_function_potential_not_computed",
  host_sample_level_link_allowed = FALSE,
  independence_group = "PRJNA766558_FFPE_tissue_21_pairs",
  conclusion_ceiling = paste(
    "paired taxonomic ecology signal conditional on no-blank limitation;",
    "not a causal microbe-host axis"
  )
)]
setorder(
  genus_differential,
  candidate_tier,
  paired_wilcoxon_q,
  paired_wilcoxon_p,
  genus_label
)

message("[7/8] ANCOM-BC2 配对随机截距正交敏感性")
ancom_counts <- t(genus_primary$filtered_counts)
ancom_metadata <- as.data.frame(metadata[, .(
  sample_name,
  tissue_role = factor(tissue_role,
                       levels = c("paired_non_tumor", "tumor")),
  patient_pair_id = factor(patient_pair_id)
)])
rownames(ancom_metadata) <- ancom_metadata$sample_name
ancom_metadata$sample_name <- NULL
set.seed(20260901L)
ancombc2_condition_log <- character()
record_ancombc2_condition <- function(condition) {
  condition_text <- conditionMessage(condition)
  if (grepl("boundary (singular) fit", condition_text, fixed = TRUE)) {
    ancombc2_condition_log <<- c(ancombc2_condition_log, condition_text)
  }
}
ancom_output <- tryCatch(
  withCallingHandlers(
    ancombc2(
      data = ancom_counts,
      taxa_are_rows = TRUE,
      meta_data = ancom_metadata,
      fix_formula = "tissue_role",
      rand_formula = "(1 | patient_pair_id)",
      p_adj_method = "BH",
      pseudo = 0,
      pseudo_sens = TRUE,
      prv_cut = 0,
      lib_cut = 0,
      s0_perc = 0.05,
      group = NULL,
      struc_zero = FALSE,
      neg_lb = FALSE,
      alpha = 0.10,
      n_cl = 2,
      verbose = FALSE,
      global = FALSE,
      pairwise = FALSE,
      dunnet = FALSE,
      trend = FALSE
    ),
    message = record_ancombc2_condition,
    warning = record_ancombc2_condition
  ),
  error = function(error) error
)
ancombc2_singular_condition_count <- length(ancombc2_condition_log)
if (inherits(ancom_output, "error")) {
  ancombc2_status <- "failed_recorded_not_used_for_gate"
  ancombc2_fit_status <- "not_available_analysis_failed"
  ancombc2_differential_count <- NA_integer_
  ancombc2_diff_robust_count <- NA_integer_
  ancombc2_error <- conditionMessage(ancom_output)
  ancombc2_results <- data.table(
    taxon = NA_character_,
    analysis_status = ancombc2_status,
    mixed_model_fit_status = ancombc2_fit_status,
    boundary_singular_condition_count = ancombc2_singular_condition_count,
    differential_taxa_count = ancombc2_differential_count,
    robust_differential_taxa_count = ancombc2_diff_robust_count,
    error_message = ancombc2_error
  )
} else {
  ancombc2_status <- "completed_random_intercept_sensitivity"
  ancombc2_fit_status <- if (ancombc2_singular_condition_count > 0L) {
    "completed_with_boundary_singular_fit_condition"
  } else {
    "completed_without_boundary_singular_fit_condition"
  }
  ancombc2_error <- ""
  if ("taxon" %in% names(ancom_output$res)) {
    ancombc2_results <- as.data.table(ancom_output$res)
  } else {
    ancombc2_results <- as.data.table(
      ancom_output$res,
      keep.rownames = "taxon"
    )
  }
  fail_if(
    anyNA(ancombc2_results$taxon) ||
      any(!nzchar(ancombc2_results$taxon)) ||
      anyDuplicated(ancombc2_results$taxon) > 0L,
    "ANCOM-BC2 结果 taxon 身份缺失或重复"
  )
  required_ancombc2_result_fields <- c(
    "diff_tissue_roletumor", "passed_ss_tissue_roletumor",
    "diff_robust_tissue_roletumor"
  )
  fail_if(
    !all(required_ancombc2_result_fields %in% names(ancombc2_results)),
    "ANCOM-BC2 结果缺少 tissue_role 差异/敏感性字段"
  )
  ancombc2_differential_count <- sum(
    ancombc2_results$diff_tissue_roletumor %in% TRUE
  )
  ancombc2_diff_robust_count <- sum(
    ancombc2_results$diff_robust_tissue_roletumor %in% TRUE
  )
  ancombc2_results[, `:=`(
    analysis_status = ancombc2_status,
    mixed_model_fit_status = ancombc2_fit_status,
    boundary_singular_condition_count = ancombc2_singular_condition_count,
    differential_taxa_count = ancombc2_differential_count,
    robust_differential_taxa_count = ancombc2_diff_robust_count,
    error_message = ""
  )]
  setorder(ancombc2_results, taxon)
}

pipeline_sensitivity <- data.table(
  pipeline = c(
    "main_min_overlap12", "min_overlap8", "forward_only",
    "main_plausible_v3v4_length_380_480nt"
  ),
  input_reads = sum(sample_qc$input_reads),
  final_prokaryotic_reads = c(
    sum(sample_qc$prokaryotic_reads_min_overlap12),
    sum(sample_qc$prokaryotic_reads_min_overlap8),
    sum(sample_qc$prokaryotic_reads_forward_only),
    sum(sample_qc$prokaryotic_reads_plausible_v3v4_length)
  ),
  final_taxa = c(
    ncol(sequence_main_final),
    ncol(sequence_overlap8_final),
    ncol(sequence_forward_final),
    ncol(sequence_main_plausible_length)
  ),
  final_genera_or_deepest_labels = c(
    ncol(genus_main), ncol(genus_overlap8), ncol(genus_forward),
    ncol(genus_plausible_length)
  ),
  bray_distance_spearman_vs_main = c(
    1,
    bray_distance_correlation(genus_main, genus_overlap8),
    bray_distance_correlation(genus_main, genus_forward),
    bray_distance_correlation(genus_main, genus_plausible_length)
  ),
  differential_direction_concordance_vs_main = c(
    1,
    mean(genus_differential$min_overlap8_direction_concordant,
         na.rm = TRUE),
    mean(genus_differential$forward_only_direction_concordant,
         na.rm = TRUE),
    mean(genus_differential$plausible_length_direction_concordant,
         na.rm = TRUE)
  )
)
pipeline_sensitivity[, final_read_retention :=
  final_prokaryotic_reads / input_reads]

message("[8/8] 写出正式对象、结构化结果、摘要与 manifest")
sequence_order <- order(
  -colSums(sequence_main_final),
  colnames(sequence_main_final)
)
sequence_main_final <- sequence_main_final[, sequence_order, drop = FALSE]
taxonomy_main_final <- taxonomy_main_final[
  colnames(sequence_main_final), , drop = FALSE
]
asv_ids <- sprintf("ASV%06d", seq_len(ncol(sequence_main_final)))
asv_counts <- as.data.table(t(sequence_main_final), keep.rownames = "sequence")
asv_counts[, asv_id := asv_ids]
setcolorder(asv_counts, c("asv_id", "sequence", rownames(sequence_main_final)))
taxonomy_table <- as.data.table(taxonomy_main_final, keep.rownames = "sequence")
taxonomy_table[, `:=`(
  asv_id = asv_ids,
  total_count = colSums(sequence_main_final),
  prevalence = colMeans(sequence_main_final > 0),
  genus_label = taxonomy_label(taxonomy_main_final),
  sequence_length_nt = nchar(sequence),
  plausible_v3v4_length_380_480nt =
    nchar(sequence) >= plausible_v3v4_length_min &
      nchar(sequence) <= plausible_v3v4_length_max
)]
setcolorder(
  taxonomy_table,
  c(
    "asv_id", "sequence", "sequence_length_nt",
    "plausible_v3v4_length_380_480nt",
    "total_count", "prevalence", "genus_label",
    setdiff(names(taxonomy_table),
            c(
              "asv_id", "sequence", "sequence_length_nt",
              "plausible_v3v4_length_380_480nt",
              "total_count", "prevalence", "genus_label"
            ))
  )
)
asv_length_qc <- taxonomy_table[, .(
  asv_count = .N,
  total_reads = sum(total_count),
  median_prevalence = median(prevalence)
), by = .(
  sequence_length_nt,
  plausible_v3v4_length_380_480nt
)]
asv_length_qc[, `:=`(
  asv_fraction = asv_count / sum(asv_count),
  read_fraction = total_reads / sum(total_reads),
  sensitivity_role = fifelse(
    plausible_v3v4_length_380_480nt,
    "broad_expected_v3v4_length_subset",
    "retained_atypical_length_not_mechanically_removed"
  )
)]
setorder(asv_length_qc, sequence_length_nt)

phyloseq_object <- phyloseq(
  otu_table(t(sequence_main_final), taxa_are_rows = TRUE),
  tax_table(taxonomy_main_final),
  sample_data(data.frame(
    metadata[, .(
      patient_pair_id,
      paper_pair_number,
      tissue_role = as.character(tissue_role),
      run_accession,
      ffpe_status,
      public_negative_control_run
    )],
    row.names = metadata$sample_name,
    check.names = FALSE
  ))
)
taxa_names(phyloseq_object) <- asv_ids
saveRDS(
  phyloseq_object,
  file.path(stage_dir, "prjna766558_dada2_phyloseq.rds"),
  compress = "xz"
)

stage_tsv(sample_qc, "prjna766558_dada2_sample_qc.tsv")
stage_tsv(asv_counts, "prjna766558_dada2_asv_counts.tsv")
stage_tsv(taxonomy_table, "prjna766558_dada2_taxonomy.tsv")
stage_tsv(asv_length_qc, "prjna766558_dada2_asv_length_qc.tsv")
stage_tsv(alpha_diversity, "prjna766558_dada2_alpha_diversity.tsv")
stage_tsv(alpha_tests, "prjna766558_dada2_alpha_paired_tests.tsv")
stage_tsv(beta_tests, "prjna766558_dada2_beta_tests.tsv")
stage_tsv(genus_differential, "prjna766558_dada2_genus_paired_differential.tsv")
stage_tsv(ancombc2_results, "prjna766558_dada2_ancombc2_sensitivity.tsv")
stage_tsv(pipeline_sensitivity, "prjna766558_dada2_pipeline_sensitivity.tsv")

candidate_counts <- genus_differential[, .N, by = candidate_tier]
setorder(candidate_counts, candidate_tier)
top_genera <- genus_differential[
  candidate_tier != "background_no_clear_paired_difference"
][order(paired_wilcoxon_q, paired_wilcoxon_p)]
top_lines <- if (nrow(top_genera)) {
  paste0(
    seq_len(min(15L, nrow(top_genera))), ". `",
    head(top_genera$genus_label, 15L), "`：CLR 配对中位差=",
    sprintf("%.3f", head(
      top_genera$median_clr_difference_tumor_minus_normal, 15L
    )), "；P=", format(head(top_genera$paired_wilcoxon_p, 15L), digits = 3),
    "；q=", format(head(top_genera$paired_wilcoxon_q, 15L), digits = 3),
    "；", head(top_genera$candidate_tier, 15L)
  )
} else {
  "- 当前柔性门禁下无属水平候选；完整结果仍保留。"
}

summary_lines <- c(
  "# PRJNA766558 配对 FFPE 组织 16S DADA2 分析摘要",
  "",
  "## 数据与流程",
  "",
  "- 21 对 ESCC 肿瘤/配对癌旁 FFPE 组织，共 42 个患者级样本。",
  "- 主流程不重复切引物、不预截断，使用 maxEE=2/2、truncQ=2、minOverlap=12；另保留 minOverlap=8 和 forward-only 敏感性。",
  paste0(
    "- 原始 FASTQ 锚定残余引物命中 ", anchored_primer_matches,
    "/", total_fastq_reads, "（",
    sprintf("%.5f", 100 * anchored_primer_fraction),
    "%），低于 0.1% 复审阈值，支持不二次切引物。"
  ),
  paste0(
    "- 宽松 380–480 nt V3–V4 长度子集保留 ",
    sum(taxonomy_table$plausible_v3v4_length_380_480nt), "/",
    nrow(taxonomy_table),
    " 个 ASV，覆盖 ", sprintf(
      "%.2f",
      100 * sum(taxonomy_table[
        plausible_v3v4_length_380_480nt == TRUE, total_count
      ]) / sum(taxonomy_table$total_count)
    ), "% 主流程 reads；其他长度仍保留并作软敏感性边界。"
  ),
  "- SILVA v138.2 DADA2 参考用于 Kingdom–Genus 分类和 V3–V4 片段的精确 species 参考匹配；该匹配不等于全长 16S 或物种级确定。非细菌/古菌及叶绿体/线粒体序列已从正式表排除。",
  paste0(
    "- species 精确匹配采用 ", species_outer_chunk_size,
    " 条外层分块、DADA2 `n=", species_inner_n,
    "` 的顺序隔离子进程；每块绑定 taxonomy/hash 并原子保存，",
    "仅改变内存执行方式，不改变匹配定义。"
  ),
  paste0(
    "- 主流程获得 ", ncol(sequence_main_final), " 个 ASV、",
    ncol(genus_main), " 个属或最深可分类标签；总最终 reads ",
    format(sum(sample_qc$prokaryotic_reads_min_overlap12), big.mark = ","),
    "。"
  ),
  "",
  "## α 多样性",
  "",
  paste0(
    "- 原始 alpha 指标作测序深度敏感描述；同时报告最终 reads 配对诊断和 ",
    rarefaction_repeats, " 次统一深度稀释敏感性。稀释深度=",
    common_rarefaction_depth, "；最小/中位最终 reads=",
    min(sample_qc$prokaryotic_reads_min_overlap12), "/",
    sprintf("%.0f", median(sample_qc$prokaryotic_reads_min_overlap12)), "。"
  ),
  paste0(
    "- ", alpha_tests$metric, "：肿瘤−癌旁配对中位差=",
    sprintf("%.3f", alpha_tests$median_paired_difference_tumor_minus_normal),
    "；P=", format(alpha_tests$paired_wilcoxon_p, digits = 3),
    "；q=", format(alpha_tests$paired_wilcoxon_q, digits = 3), "。"
  ),
  "",
  "## β 多样性",
  "",
  paste0(
    "- ", beta_tests$distance, "：患者对内限制置换 PERMANOVA R²=",
    sprintf("%.3f", beta_tests$permanova_r2),
    "；P=", format(beta_tests$restricted_permutation_p, digits = 3),
    "；q=", format(beta_tests$restricted_permutation_q, digits = 3),
    "；betadisper P=", format(beta_tests$betadisper_p, digits = 3), "。"
  ),
  "",
  "## 属水平候选",
  "",
  paste0(
    "- 分层计数：",
    paste0(candidate_counts$candidate_tier, "=", candidate_counts$N,
           collapse = "；"), "。"
  ),
  top_lines,
  "",
  "## 敏感性与证据边界",
  "",
  paste0(
    "- ANCOM-BC2 随机截距敏感性状态：`", ancombc2_status, "`",
    if (nzchar(ancombc2_error)) paste0("；错误：", ancombc2_error) else "。"
  ),
  paste0(
    "- ANCOM-BC2 混合模型拟合状态：`", ancombc2_fit_status,
    "`；运行时记录 boundary singular condition ",
    ancombc2_singular_condition_count, " 次；初始差异条目 ",
    ifelse(is.na(ancombc2_differential_count), "NA", ancombc2_differential_count),
    " 个，pseudo-count 敏感性后 `diff_robust` ",
    ifelse(is.na(ancombc2_diff_robust_count), "NA", ancombc2_diff_robust_count),
    " 个。该模型只作正交敏感性，不参与属级候选门禁。"
  ),
  "- 当前没有公开阴性对照、提取空白或 DNA 浓度，未运行 decontam；流程稳定不能替代污染排除。",
  "- 所有微生物候选最高为带无空白限制的配对生态信号；16S 功能潜力本轮未计算，不能写成实测功能。",
  "- PRJNA766558 与 TCGA、Cao 或 PR001876 不是同一批患者，不建立虚假菌群—宿主患者级相关或因果边。"
)
writeLines(
  summary_lines,
  file.path(stage_dir, "prjna766558_dada2_summary.md"),
  useBytes = TRUE
)

artifact_names <- c(
  "prjna766558_dada2_phyloseq.rds",
  "prjna766558_dada2_sample_qc.tsv",
  "prjna766558_dada2_asv_counts.tsv",
  "prjna766558_dada2_taxonomy.tsv",
  "prjna766558_dada2_asv_length_qc.tsv",
  "prjna766558_dada2_alpha_diversity.tsv",
  "prjna766558_dada2_alpha_paired_tests.tsv",
  "prjna766558_dada2_beta_tests.tsv",
  "prjna766558_dada2_genus_paired_differential.tsv",
  "prjna766558_dada2_ancombc2_sensitivity.tsv",
  "prjna766558_dada2_pipeline_sensitivity.tsv",
  "prjna766558_dada2_summary.md"
)
artifact_paths <- file.path(stage_dir, artifact_names)
fail_if(any(!file_exists(artifact_paths)), "PRJNA766558 正式产物未完整生成")
manifest <- data.table(
  artifact = artifact_names,
  relative_path = file.path("results", artifact_names),
  file_size_bytes = as.numeric(file_info(artifact_paths)$size),
  sha256 = vapply(
    artifact_paths,
    digest,
    FUN.VALUE = character(1),
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  ),
  generated_date = as.character(Sys.Date()),
  generation_script = "scripts/18_analyze_prjna766558_dada2.R",
  executed_script_sha256 = executed_script_sha256,
  pipeline_schema_version = pipeline_schema_version,
  pipeline_signature_sha256 = pipeline_signature,
  species_execution_schema_version = species_execution_schema_version,
  species_execution_signature_sha256 = species_execution_signature,
  species_outer_chunk_size = species_outer_chunk_size,
  species_inner_n = species_inner_n,
  signature_component_count = length(signature_components),
  independence_group = "PRJNA766558_FFPE_tissue_21_pairs",
  status = "verified"
)
stage_tsv(manifest, "prjna766558_dada2_artifact_manifest.tsv")

for (artifact in artifact_names) {
  atomic_publish_file(
    file.path(stage_dir, artifact),
    file.path(results_dir, artifact)
  )
}
# manifest 最后原子替换：若中途中断，旧 manifest 会与已部分更新的
# artifact 不一致而使下游校验失败，不会把混合版本误认为 verified。
atomic_publish_file(
  file.path(stage_dir, "prjna766558_dada2_artifact_manifest.tsv"),
  file.path(results_dir, "prjna766558_dada2_artifact_manifest.tsv")
)
published_paths <- file.path(results_dir, artifact_names)
published_size <- as.numeric(file_info(published_paths)$size)
published_sha <- vapply(
  published_paths,
  digest,
  FUN.VALUE = character(1),
  algo = "sha256",
  file = TRUE,
  serialize = FALSE
)
fail_if(
  any(published_size != manifest$file_size_bytes) ||
    any(published_sha != manifest$sha256),
  "results/ 正式 artifact 发布后大小/SHA256 回读不一致"
)
published_manifest <- fread(
  file.path(results_dir, "prjna766558_dada2_artifact_manifest.tsv"),
  colClasses = "character",
  na.strings = NULL,
  showProgress = FALSE
)
stage_manifest_reread <- fread(
  file.path(stage_dir, "prjna766558_dada2_artifact_manifest.tsv"),
  colClasses = "character",
  na.strings = NULL,
  showProgress = FALSE
)
fail_if(
  !identical(names(published_manifest), names(stage_manifest_reread)) ||
    nrow(published_manifest) != nrow(stage_manifest_reread) ||
    any(!vapply(
      names(stage_manifest_reread),
      function(column_name) identical(
        published_manifest[[column_name]],
        stage_manifest_reread[[column_name]]
      ),
      logical(1)
    )),
  "results/ 正式 manifest 与 stage 冻结版不一致"
)

dir_delete(stage_dir)
message(
  "完成：", ncol(sequence_main_final), " 个主流程 ASV，",
  ncol(genus_main), " 个属/最深标签，",
  nrow(genus_differential), " 个 prevalence>=20% 配对检验特征。"
)
