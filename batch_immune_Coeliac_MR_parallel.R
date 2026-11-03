#!/usr/bin/env Rscript

# ==============================================================================
# 批量 two-sample MR 分析脚本
# Exposure：GWAS Catalog harmonised 免疫表型（*.h.tsv.gz）
# Outcome ：FinnGen R13 K11_COELIAC（Coeliac disease）
#
# 该脚本根据已经完成的单个 Exposure 测试流程编写，并修复以下常见问题：
# 1. 不完整读入 FinnGen 文件；汇总全部 clumped SNP 后只流式扫描一次。
# 2. harmonise_data(action = 3)，删除全部回文 SNP。
# 3. 对每个 Exposure 使用 tryCatch，单个文件报错不会终止整个批量流程。
# 4. 按最终 SNP 数选择 Wald ratio、IVW、weighted median 和 MR-Egger。
# 5. 自动汇总异质性、Egger 截距、leave-one-out 和 MR-PRESSO。
# 6. 自动进行全局 BH-FDR 校正。


#
# 推荐运行方式：
# conda activate mr_env
# Rscript batch_immune_Coeliac_MR_parallel.R
# ==============================================================================


# ==============================================================================
# 0. 加载依赖包
# ==============================================================================

required_packages <- c("data.table", "TwoSampleMR", "ieugwasr", "MRPRESSO")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "缺少 R 包：",
    paste(missing_packages, collapse = ", "),
    "\n请先安装后再运行。"
  )
}

library(data.table)


# ==============================================================================
# 1. 用户参数：请首先修改这里
# ==============================================================================

# 所有 Exposure *.h.tsv.gz 文件所在目录
exposure_dir <- path.expand("./MR_project/data/exposure")

# Exposure 元数据文件
# 必须包含：
# exposure_id, exposure_name, file_name, group, include
metadata_file <- path.expand(
  "./MR_project/config/exposure_metadata_all_731_extended.tsv"
)


# FinnGen outcome 文件
outcome_file <- path.expand(
  "./MR_project/data/outcome/summary_stats_finngen_R13_K11_COELIAC.gz"
)

# 1000 Genomes EUR PLINK 文件共同前缀
# 例如该路径应对应 EUR.bed、EUR.bim、EUR.fam
ld_reference <- path.expand("./MR_project/reference/1000G/EUR")

# PLINK 路径；若已加入 PATH，一般可自动识别
plink_bin <- Sys.which("plink")

# 输出目录
result_root <- path.expand(
  "./MR_project/results/batch_all_731_Coeliac_MR_parallel_NonMHC"
)

# 工具变量筛选参数
p_threshold <- 1e-5
f_threshold <- 10

# KSD 参考框架的 LD clumping 参数
clump_r2 <- 0.1
clump_kb <- 500

# 是否排除扩展 MHC 区域
exclude_mhc <- FALSE

# 扩展 MHC 区域定义
mhc_chr <- "6"
mhc_start <- 25000000
mhc_end <- 34000000

# 删除全部 palindromic SNP
harmonise_action <- 3

# MR-PRESSO 初筛重复数
# 全部批量初筛设为 1000
presso_nb_distribution <- 1000

# Outcome 信息
outcome_id <- "FinnGen_R13_K11_COELIAC"

outcome_name <- "Coeliac disease"

# -----------------------------并行计算参数--------------------------------------
# 第一阶段需要同时解压大型Exposure文件并调用PLINK，内存和磁盘压力最大。
# 建议拥有32 GB内存时使用2；64 GB及以上可尝试3–4。
clump_workers <- 3L

# 第三阶段仅处理clumping后的较小数据，并运行MR和MR-PRESSO。
# 建议先使用6–8；不要直接设置为全部逻辑核心数。
mr_workers <- 8L

# 是否复用已经完成的单Exposure检查点。
resume_completed <- TRUE

# 上次状态为error的任务是否重新运行。
rerun_error_tasks <- TRUE

# MR-PRESSO包含随机抽样；使用固定种子保证单个Exposure可重复。
parallel_seed <- 20260722L

# 本并行版本使用parallel::mclapply，仅支持Linux/macOS。
if (.Platform$OS.type != "unix") {
  stop("该并行脚本使用mclapply，仅支持Linux/macOS。")
}

# 防止每个并行子进程再次调用多个data.table/BLAS线程，造成CPU过度订阅。
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

data.table::setDTthreads(1L)

physical_cores <- parallel::detectCores(logical = FALSE)
if (is.na(physical_cores) || physical_cores < 1L) {
  physical_cores <- 1L
}

# 至少保留一个物理核心给系统；同时不能超过Exposure任务数。
max_safe_workers <- max(1L, physical_cores - 1L)
clump_workers <- max(1L, min(clump_workers, max_safe_workers))
mr_workers <- max(1L, min(mr_workers, max_safe_workers))


# ==============================================================================
# 2. 输出目录
# ==============================================================================

dir_clump <- file.path(result_root, "01_clumped_exposure")
dir_harmonised <- file.path(result_root, "02_harmonised")
dir_detail <- file.path(result_root, "03_MR_details")
dir_log <- file.path(result_root, "04_logs")
dir_summary <- file.path(result_root, "05_summary")
dir_cache_clump <- file.path(result_root, "06_cache_clump")
dir_cache_mr <- file.path(result_root, "07_cache_MR")

for (x in c(
  result_root, dir_clump, dir_harmonised,
  dir_detail, dir_log, dir_summary,
  dir_cache_clump, dir_cache_mr
)) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}


# ==============================================================================
# 3. 输入文件和软件检查
# ==============================================================================

if (!dir.exists(exposure_dir)) {
  stop("Exposure 目录不存在：", exposure_dir)
}

if (!file.exists(metadata_file)) {
  stop("Exposure 元数据文件不存在：", metadata_file)
}

if (!file.exists(outcome_file)) {
  stop("FinnGen outcome 文件不存在：", outcome_file)
}

for (suffix in c(".bed", ".bim", ".fam")) {
  if (!file.exists(paste0(ld_reference, suffix))) {
    stop("缺少 LD reference 文件：", paste0(ld_reference, suffix))
  }
}

if (identical(plink_bin, "") || is.na(plink_bin)) {
  stop(
    "未找到 PLINK。请确认终端可运行 plink，",
    "或在脚本中手动设置 plink_bin。"
  )
}


# ==============================================================================
# 4. 读取 Exposure 元数据
# ==============================================================================

metadata <- fread(metadata_file)

required_meta <- c(
  "exposure_id",
  "exposure_name",
  "file_name",
  "group",
  "include"
)

missing_meta <- setdiff(required_meta, names(metadata))

if (length(missing_meta) > 0) {
  stop(
    "元数据表缺少以下列：",
    paste(missing_meta, collapse = ", ")
  )
}

# include 支持 TRUE/FALSE、1/0、yes/no
metadata[, include := tolower(trimws(as.character(include))) %in%
           c("true", "1", "yes", "y")]

metadata <- metadata[include == TRUE]

if (nrow(metadata) == 0) {
  stop("元数据中没有 include=TRUE 的 Exposure。")
}

metadata[, exposure_file := file.path(exposure_dir, file_name)]

missing_files <- metadata[!file.exists(exposure_file)]

if (nrow(missing_files) > 0) {
  fwrite(
    missing_files,
    file.path(dir_log, "missing_exposure_files.tsv"),
    sep = "\t"
  )
  stop("存在缺失 Exposure 文件，请查看 missing_exposure_files.tsv。")
}


# ==============================================================================
# 5. 通用辅助函数
# ==============================================================================

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

empty_presso_summary <- function() {
  data.table(
    presso_global_p = NA_character_,
    presso_n_outliers = NA_integer_,
    presso_outlier_snps = NA_character_,
    presso_raw_beta = NA_real_,
    presso_raw_se = NA_real_,
    presso_raw_p = NA_real_,
    presso_corrected_beta = NA_real_,
    presso_corrected_se = NA_real_,
    presso_corrected_p = NA_real_,
    presso_distortion_p = NA_character_
  )
}

empty_mr_summary <- function(
    exposure_id,
    exposure_name,
    group,
    status,
    error_message = NA_character_) {

  data.table(
    exposure_id = exposure_id,
    exposure_name = exposure_name,
    group = group,
    n_matched_outcome = NA_integer_,
    n_harmonised = NA_integer_,
    status = status,
    error_message = error_message,

    wald_beta = NA_real_,
    wald_se = NA_real_,
    wald_p = NA_real_,
    wald_or = NA_real_,
    wald_lci = NA_real_,
    wald_uci = NA_real_,

    ivw_beta = NA_real_,
    ivw_se = NA_real_,
    ivw_p = NA_real_,
    ivw_or = NA_real_,
    ivw_lci = NA_real_,
    ivw_uci = NA_real_,

    wm_beta = NA_real_,
    wm_se = NA_real_,
    wm_p = NA_real_,
    wm_or = NA_real_,
    wm_lci = NA_real_,
    wm_uci = NA_real_,

    egger_beta = NA_real_,
    egger_se = NA_real_,
    egger_p = NA_real_,
    egger_or = NA_real_,
    egger_lci = NA_real_,
    egger_uci = NA_real_,

    ivw_Q = NA_real_,
    ivw_Q_df = NA_real_,
    ivw_Q_p = NA_real_,
    egger_Q = NA_real_,
    egger_Q_df = NA_real_,
    egger_Q_p = NA_real_,

    egger_intercept = NA_real_,
    egger_intercept_se = NA_real_,
    egger_intercept_p = NA_real_,

    presso_error = NA_character_
  )
}


# ==============================================================================
# 6. 处理单个 Exposure：读取、QC、F值、MHC排除和LD clumping
# ==============================================================================

prepare_one_exposure <- function(
    exposure_file,
    exposure_id,
    exposure_name,
    group) {

  raw <- fread(exposure_file)

  # 根据已经测试通过的 GWAS Catalog harmonised 文件列名进行映射
  required_cols <- c(
    "hm_rsid",
    "hm_chrom",
    "hm_pos",
    "hm_other_allele",
    "hm_effect_allele",
    "hm_beta",
    "hm_effect_allele_frequency",
    "standard_error",
    "p_value",
    "n"
  )

  missing_cols <- setdiff(required_cols, names(raw))

  if (length(missing_cols) > 0) {
    stop(
      "Exposure 文件缺少列：",
      paste(missing_cols, collapse = ", ")
    )
  }

  # 必须将 harmonised 等位基因与 harmonised beta 配套使用
  exp <- raw[, .(
    SNP = trimws(as.character(hm_rsid)),
    chr = as.character(hm_chrom),
    pos = safe_numeric(hm_pos),
    effect_allele = toupper(trimws(as.character(hm_effect_allele))),
    other_allele = toupper(trimws(as.character(hm_other_allele))),
    beta = safe_numeric(hm_beta),
    se = safe_numeric(standard_error),
    eaf = safe_numeric(hm_effect_allele_frequency),
    pval = safe_numeric(p_value),
    samplesize = safe_numeric(n)
  )]

  # 仅保留标准 rsID 和关键字段完整的数据
  exp <- exp[
    !is.na(SNP) &
      grepl("^rs[0-9]+$", SNP) &
      !is.na(beta) &
      !is.na(se) &
      se > 0 &
      !is.na(pval) &
      !is.na(effect_allele) &
      effect_allele != "" &
      !is.na(other_allele) &
      other_allele != ""
  ]

  # 同一 rsID 若重复，保留 P 值最小的一行
  setorder(exp, SNP, pval)
  exp <- exp[!duplicated(SNP)]

  n_valid <- nrow(exp)

  # P < 1e-5
  exp <- exp[pval < p_threshold]
  n_p_threshold <- nrow(exp)

  if (n_p_threshold == 0) {
    return(list(
      status = "no_p_threshold_instruments",
      qc = data.table(
        exposure_id, exposure_name, group,
        n_valid,
        n_p_threshold = 0L,
        n_F_ge_10 = 0L,
        n_MHC_removed = 0L,
        n_after_MHC = 0L,
        n_after_clump = 0L
      ),
      exposure_dat = NULL,
      clump_result = NULL
    ))
  }

  # F = (beta/se)^2
  exp[, F_stat := (beta / se)^2]

  exp <- exp[F_stat >= f_threshold]
  n_F_ge_10 <- nrow(exp)

  if (n_F_ge_10 == 0) {
    return(list(
      status = "no_strong_instruments",
      qc = data.table(
        exposure_id, exposure_name, group,
        n_valid,
        n_p_threshold,
        n_F_ge_10 = 0L,
        n_MHC_removed = 0L,
        n_after_MHC = 0L,
        n_after_clump = 0L
      ),
      exposure_dat = NULL,
      clump_result = NULL
    ))
  }

  # MHC 排除应在 clumping 前完成
  n_MHC_removed <- 0L

  if (exclude_mhc) {
    mhc_flag <- (
      exp$chr == mhc_chr &
        exp$pos >= mhc_start &
        exp$pos <= mhc_end
    )

    mhc_flag[is.na(mhc_flag)] <- FALSE
    n_MHC_removed <- sum(mhc_flag)
    exp <- exp[!mhc_flag]
  }

  n_after_MHC <- nrow(exp)

  if (n_after_MHC == 0) {
    return(list(
      status = "no_nonMHC_instruments",
      qc = data.table(
        exposure_id, exposure_name, group,
        n_valid,
        n_p_threshold,
        n_F_ge_10,
        n_MHC_removed,
        n_after_MHC = 0L,
        n_after_clump = 0L
      ),
      exposure_dat = NULL,
      clump_result = NULL
    ))
  }

  # 修复 phenotype_col=NULL 报错：创建真实 phenotype 列
  exp[, phenotype := exposure_name]

  # 修复 data.table 与 format_data 的兼容问题：强制 as.data.frame()
  exposure_dat <- TwoSampleMR::format_data(
    dat = as.data.frame(exp),
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "pval",
    chr_col = "chr",
    pos_col = "pos",
    samplesize_col = "samplesize",
    phenotype_col = "phenotype"
  )

  exposure_dat$id.exposure <- exposure_id
  exposure_dat$exposure <- exposure_name

  # 不使用 transmute()，直接构建普通 data.frame
  clump_input <- data.frame(
    rsid = as.character(exposure_dat$SNP),
    pval = as.numeric(exposure_dat$pval.exposure),
    id = as.character(exposure_dat$id.exposure),
    stringsAsFactors = FALSE
  )

  clump_input <- clump_input[
    !is.na(clump_input$rsid) &
      !is.na(clump_input$pval) &
      !duplicated(clump_input$rsid),
    ,
    drop = FALSE
  ]

  clump_result <- ieugwasr::ld_clump_local(
    dat = clump_input,
    clump_kb = clump_kb,
    clump_r2 = clump_r2,
    clump_p = p_threshold,
    bfile = ld_reference,
    plink_bin = plink_bin
  )

  n_after_clump <- nrow(clump_result)

  if (n_after_clump == 0) {
    return(list(
      status = "no_instruments_after_clump",
      qc = data.table(
        exposure_id, exposure_name, group,
        n_valid,
        n_p_threshold,
        n_F_ge_10,
        n_MHC_removed,
        n_after_MHC,
        n_after_clump = 0L
      ),
      exposure_dat = NULL,
      clump_result = NULL
    ))
  }

  # clump_result 只包含 rsid/pval/id，需要回连完整 Exposure 信息
  exposure_clumped <- exposure_dat[
    exposure_dat$SNP %in% clump_result$rsid,
    ,
    drop = FALSE
  ]

  # 按 clump 结果顺序排列，方便追踪
  exposure_clumped <- exposure_clumped[
    match(clump_result$rsid, exposure_clumped$SNP),
    ,
    drop = FALSE
  ]

  exposure_clumped <- exposure_clumped[
    !is.na(exposure_clumped$SNP),
    ,
    drop = FALSE
  ]

  list(
    status = "success",
    qc = data.table(
      exposure_id, exposure_name, group,
      n_valid,
      n_p_threshold,
      n_F_ge_10,
      n_MHC_removed,
      n_after_MHC,
      n_after_clump = nrow(exposure_clumped)
    ),
    exposure_dat = exposure_clumped,
    clump_result = clump_result
  )
}


# ==============================================================================
# 7. 一次性流式提取 FinnGen 中的所有目标 SNP
# ==============================================================================

extract_finngen_target_snps <- function(
    outcome_file,
    target_snps,
    output_file) {

  target_snps <- unique(trimws(as.character(target_snps)))
  target_snps <- target_snps[grepl("^rs[0-9]+$", target_snps)]

  if (length(target_snps) == 0) {
    stop("没有标准 rsID 可用于 FinnGen 提取。")
  }

  # 只读取表头，不完整载入 outcome
  outcome_header <- fread(outcome_file, nrows = 0)
  rsid_col <- match("rsids", names(outcome_header))

  if (is.na(rsid_col)) {
    stop("FinnGen 文件中不存在 rsids 列。")
  }

  snp_file <- tempfile(fileext = ".txt")
  writeLines(target_snps, snp_file)

  # awk 指定 rsids 列精确匹配，避免在其他字段误命中
  awk_code <- sprintf(
    paste0(
      "BEGIN{FS=OFS=\"\\t\"} ",
      "NR==FNR{wanted[$1]=1; next} ",
      "FNR==1{print; next} ",
      "$%d in wanted{print}"
    ),
    rsid_col
  )

  shell_inner <- sprintf(
    "awk '%s' %s <(gzip -dc %s) > %s",
    awk_code,
    shQuote(snp_file),
    shQuote(outcome_file),
    shQuote(output_file)
  )

  status <- system2(
    "bash",
    c("-c", shQuote(shell_inner))
  )

  if (status != 0 || !file.exists(output_file)) {
    stop("FinnGen 流式提取失败。")
  }

  extracted <- fread(output_file)

  if (nrow(extracted) == 0) {
    stop("FinnGen 中未匹配到任何目标 SNP。")
  }

  extracted
}


# ==============================================================================
# 8. 将 FinnGen 子集整理为 TwoSampleMR outcome
# ==============================================================================

prepare_outcome <- function(outcome_subset) {

  required_cols <- c(
    "#chrom",
    "pos",
    "ref",
    "alt",
    "rsids",
    "pval",
    "beta",
    "sebeta",
    "af_alt"
  )

  missing_cols <- setdiff(required_cols, names(outcome_subset))

  if (length(missing_cols) > 0) {
    stop(
      "FinnGen 子集缺少列：",
      paste(missing_cols, collapse = ", ")
    )
  }

  outcome_clean <- data.table(
    SNP = trimws(as.character(outcome_subset$rsids)),
    beta = safe_numeric(outcome_subset$beta),
    se = safe_numeric(outcome_subset$sebeta),
    effect_allele = toupper(trimws(as.character(outcome_subset$alt))),
    other_allele = toupper(trimws(as.character(outcome_subset$ref))),
    eaf = safe_numeric(outcome_subset$af_alt),
    pval = safe_numeric(outcome_subset$pval),
    chr = as.character(outcome_subset[["#chrom"]]),
    pos = safe_numeric(outcome_subset$pos),
    phenotype = outcome_name
  )

  outcome_clean <- outcome_clean[
    !is.na(SNP) &
      grepl("^rs[0-9]+$", SNP) &
      !is.na(beta) &
      !is.na(se) &
      se > 0 &
      !is.na(effect_allele) &
      effect_allele != "" &
      !is.na(other_allele) &
      other_allele != ""
  ]

  setorder(outcome_clean, SNP, pval)
  outcome_clean <- outcome_clean[!duplicated(SNP)]

  # 不使用 phenotype_col=NULL
  outcome_dat <- TwoSampleMR::format_data(
    dat = as.data.frame(outcome_clean),
    type = "outcome",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "pval",
    chr_col = "chr",
    pos_col = "pos",
    phenotype_col = "phenotype"
  )

  outcome_dat$id.outcome <- outcome_id
  outcome_dat$outcome <- outcome_name

  outcome_dat
}


# ==============================================================================
# 9. 提取某个 MR 方法的结果
# ==============================================================================

extract_method_result <- function(mr_or_dt, method_name) {

  x <- mr_or_dt[method == method_name]

  if (nrow(x) == 0) {
    return(list(
      beta = NA_real_,
      se = NA_real_,
      p = NA_real_,
      or = NA_real_,
      lci = NA_real_,
      uci = NA_real_
    ))
  }

  list(
    beta = x$b[1],
    se = x$se[1],
    p = x$pval[1],
    or = x$or[1],
    lci = x$or_lci95[1],
    uci = x$or_uci95[1]
  )
}


# ==============================================================================
# 10. 提取 MR-PRESSO 汇总结果
# ==============================================================================

extract_presso_result <- function(presso_result, harmonised_keep) {

  ans <- empty_presso_summary()

  if (is.null(presso_result)) {
    return(ans)
  }

  main_results <- presso_result[["Main MR results"]]

  if (!is.null(main_results) && nrow(main_results) > 0) {

    raw_row <- main_results[
      main_results$`MR Analysis` == "Raw",
      ,
      drop = FALSE
    ]

    if (nrow(raw_row) > 0) {
      ans$presso_raw_beta <- raw_row$`Causal Estimate`[1]
      ans$presso_raw_se <- raw_row$Sd[1]
      ans$presso_raw_p <- raw_row$`P-value`[1]
    }

    corrected_row <- main_results[
      main_results$`MR Analysis` == "Outlier-corrected",
      ,
      drop = FALSE
    ]

    if (nrow(corrected_row) > 0) {
      ans$presso_corrected_beta <- corrected_row$`Causal Estimate`[1]
      ans$presso_corrected_se <- corrected_row$Sd[1]
      ans$presso_corrected_p <- corrected_row$`P-value`[1]
    }
  }

  presso_block <- presso_result[["MR-PRESSO results"]]

  if (is.null(presso_block)) {
    return(ans)
  }

  global_test <- presso_block[["Global Test"]]

  if (!is.null(global_test)) {
    ans$presso_global_p <- as.character(global_test$Pvalue)
  }

  distortion_test <- presso_block[["Distortion Test"]]

  if (!is.null(distortion_test)) {

    ans$presso_distortion_p <- as.character(
      distortion_test$Pvalue
    )

    indices <- distortion_test[["Outliers Indices"]]

    if (!is.null(indices) &&
        length(indices) > 0 &&
        !all(is.na(indices))) {

      indices <- suppressWarnings(as.integer(indices))
      indices <- indices[
        !is.na(indices) &
          indices >= 1 &
          indices <= nrow(harmonised_keep)
      ]

      ans$presso_n_outliers <- length(indices)

      if (length(indices) > 0) {
        ans$presso_outlier_snps <- paste(
          harmonised_keep$SNP[indices],
          collapse = ";"
        )
      }
    } else {
      ans$presso_n_outliers <- 0L
    }
  }

  ans
}


# ==============================================================================
# 11. 对单个 Exposure 运行 harmonisation、MR 和敏感性分析
# ==============================================================================

run_one_mr <- function(
    exposure_dat,
    outcome_dat,
    exposure_id,
    exposure_name,
    group) {

  out_dir <- file.path(dir_detail, exposure_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # 只取当前 Exposure 的 SNP
  outcome_sub <- outcome_dat[
    outcome_dat$SNP %in% exposure_dat$SNP,
    ,
    drop = FALSE
  ]

  n_matched_outcome <- length(unique(outcome_sub$SNP))

  if (n_matched_outcome == 0) {
    ans <- empty_mr_summary(
      exposure_id,
      exposure_name,
      group,
      status = "no_outcome_match"
    )
    ans$n_matched_outcome <- 0L
    ans$n_harmonised <- 0L
    return(ans)
  }

  harmonised_all <- TwoSampleMR::harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat = outcome_sub,
    action = harmonise_action
  )

  harmonised_keep <- harmonised_all[
    harmonised_all$mr_keep,
    ,
    drop = FALSE
  ]

  # 保险性去重
  harmonised_keep <- harmonised_keep[
    !duplicated(harmonised_keep$SNP),
    ,
    drop = FALSE
  ]

  n_harmonised <- nrow(harmonised_keep)

  fwrite(
    as.data.table(harmonised_all),
    file.path(
      dir_harmonised,
      paste0(exposure_id, "_harmonised_all.tsv")
    ),
    sep = "\t"
  )

  fwrite(
    as.data.table(harmonised_keep),
    file.path(
      dir_harmonised,
      paste0(exposure_id, "_harmonised_keep.tsv")
    ),
    sep = "\t"
  )

  if (n_harmonised == 0) {
    ans <- empty_mr_summary(
      exposure_id,
      exposure_name,
      group,
      status = "no_harmonised_instruments"
    )
    ans$n_matched_outcome <- n_matched_outcome
    ans$n_harmonised <- 0L
    return(ans)
  }

  # SNP 数量决定 MR 方法
  if (n_harmonised == 1) {
    method_list <- "mr_wald_ratio"
  } else if (n_harmonised == 2) {
    method_list <- "mr_ivw"
  } else {
    method_list <- c(
      "mr_ivw",
      "mr_weighted_median",
      "mr_egger_regression"
    )
  }

  mr_result <- TwoSampleMR::mr(
    harmonised_keep,
    method_list = method_list
  )

  mr_or_dt <- as.data.table(
    TwoSampleMR::generate_odds_ratios(mr_result)
  )

  fwrite(
    mr_or_dt,
    file.path(out_dir, "MR_OR.tsv"),
    sep = "\t"
  )

  wald <- extract_method_result(mr_or_dt, "Wald ratio")
  ivw <- extract_method_result(mr_or_dt, "Inverse variance weighted")
  wm <- extract_method_result(mr_or_dt, "Weighted median")
  egger <- extract_method_result(mr_or_dt, "MR Egger")

  # -------------------------异质性检验-------------------------
  ivw_Q <- ivw_Q_df <- ivw_Q_p <- NA_real_
  egger_Q <- egger_Q_df <- egger_Q_p <- NA_real_

  if (n_harmonised >= 2) {

    heterogeneity_result <- tryCatch(
      TwoSampleMR::mr_heterogeneity(harmonised_keep),
      error = function(e) NULL
    )

    if (!is.null(heterogeneity_result) &&
        nrow(heterogeneity_result) > 0) {

      heterogeneity_dt <- as.data.table(heterogeneity_result)

      fwrite(
        heterogeneity_dt,
        file.path(out_dir, "heterogeneity.tsv"),
        sep = "\t"
      )

      x <- heterogeneity_dt[
        method == "Inverse variance weighted"
      ]

      if (nrow(x) > 0) {
        ivw_Q <- x$Q[1]
        ivw_Q_df <- x$Q_df[1]
        ivw_Q_p <- x$Q_pval[1]
      }

      x <- heterogeneity_dt[
        method == "MR Egger"
      ]

      if (nrow(x) > 0) {
        egger_Q <- x$Q[1]
        egger_Q_df <- x$Q_df[1]
        egger_Q_p <- x$Q_pval[1]
      }
    }
  }

  # -------------------------Egger 截距-------------------------
  egger_intercept <- NA_real_
  egger_intercept_se <- NA_real_
  egger_intercept_p <- NA_real_

  if (n_harmonised >= 3) {

    pleiotropy_result <- tryCatch(
      TwoSampleMR::mr_pleiotropy_test(harmonised_keep),
      error = function(e) NULL
    )

    if (!is.null(pleiotropy_result) &&
        nrow(pleiotropy_result) > 0) {

      pleiotropy_dt <- as.data.table(pleiotropy_result)

      fwrite(
        pleiotropy_dt,
        file.path(out_dir, "egger_intercept.tsv"),
        sep = "\t"
      )

      egger_intercept <- pleiotropy_dt$egger_intercept[1]
      egger_intercept_se <- pleiotropy_dt$se[1]
      egger_intercept_p <- pleiotropy_dt$pval[1]
    }
  }

  # -------------------------Leave-one-out-----------------------
  if (n_harmonised >= 3) {

    loo_result <- tryCatch(
      TwoSampleMR::mr_leaveoneout(harmonised_keep),
      error = function(e) NULL
    )

    if (!is.null(loo_result)) {
      fwrite(
        as.data.table(loo_result),
        file.path(out_dir, "leave_one_out.tsv"),
        sep = "\t"
      )
    }
  }

  # -------------------------单 SNP 分析-------------------------
  single_result <- tryCatch(
    TwoSampleMR::mr_singlesnp(harmonised_keep),
    error = function(e) NULL
  )

  if (!is.null(single_result)) {
    fwrite(
      as.data.table(single_result),
      file.path(out_dir, "single_SNP.tsv"),
      sep = "\t"
    )
  }

  # -------------------------MR-PRESSO---------------------------
  presso_summary <- empty_presso_summary()
  presso_error <- NA_character_

  # 至少 4 SNP 才尝试 MR-PRESSO
  if (n_harmonised >= 4) {

    presso_result <- tryCatch(
      MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = harmonised_keep,
        NbDistribution = presso_nb_distribution,
        SignifThreshold = 0.05
      ),
      error = function(e) {
        presso_error <<- conditionMessage(e)
        NULL
      }
    )

    if (!is.null(presso_result)) {

      saveRDS(
        presso_result,
        file.path(out_dir, "MR_PRESSO_result.rds")
      )

      capture.output(
        print(presso_result),
        file = file.path(out_dir, "MR_PRESSO_result.txt")
      )

      presso_summary <- extract_presso_result(
        presso_result,
        harmonised_keep
      )
    }
  }

  ans <- data.table(
    exposure_id = exposure_id,
    exposure_name = exposure_name,
    group = group,
    n_matched_outcome = n_matched_outcome,
    n_harmonised = n_harmonised,
    status = "success",
    error_message = NA_character_,

    wald_beta = wald$beta,
    wald_se = wald$se,
    wald_p = wald$p,
    wald_or = wald$or,
    wald_lci = wald$lci,
    wald_uci = wald$uci,

    ivw_beta = ivw$beta,
    ivw_se = ivw$se,
    ivw_p = ivw$p,
    ivw_or = ivw$or,
    ivw_lci = ivw$lci,
    ivw_uci = ivw$uci,

    wm_beta = wm$beta,
    wm_se = wm$se,
    wm_p = wm$p,
    wm_or = wm$or,
    wm_lci = wm$lci,
    wm_uci = wm$uci,

    egger_beta = egger$beta,
    egger_se = egger$se,
    egger_p = egger$p,
    egger_or = egger$or,
    egger_lci = egger$lci,
    egger_uci = egger$uci,

    ivw_Q = ivw_Q,
    ivw_Q_df = ivw_Q_df,
    ivw_Q_p = ivw_Q_p,
    egger_Q = egger_Q,
    egger_Q_df = egger_Q_df,
    egger_Q_p = egger_Q_p,

    egger_intercept = egger_intercept,
    egger_intercept_se = egger_intercept_se,
    egger_intercept_p = egger_intercept_p,

    presso_error = presso_error
  )

  cbind(ans, presso_summary)
}


# ==============================================================================
# 12. 第一阶段：并行读取Exposure并完成本地LD clumping
# ==============================================================================

message("阶段 1/3：并行处理Exposure并运行本地LD clumping。")
message("clumping并发数：", clump_workers)

# 单个clumping任务。每个任务只写自己的文件，不会发生并发覆盖。
run_one_clump_task <- function(i) {

  data.table::setDTthreads(1L)

  exposure_id_i <- metadata$exposure_id[i]
  exposure_name_i <- metadata$exposure_name[i]
  group_i <- metadata$group[i]
  exposure_file_i <- metadata$exposure_file[i]

  exposure_cache_file <- file.path(
    dir_cache_clump,
    paste0(exposure_id_i, "_exposure_dat.rds")
  )

  task_cache_file <- file.path(
    dir_cache_clump,
    paste0(exposure_id_i, "_task_result.rds")
  )

  # 只复用完整成功的clumping任务。错误任务默认重新运行。
  if (resume_completed &&
      file.exists(task_cache_file) &&
      file.exists(exposure_cache_file)) {

    cached <- tryCatch(readRDS(task_cache_file), error = function(e) NULL)

    if (!is.null(cached) && isTRUE(cached$success)) {
      message(sprintf("[clump %d/%d] %s：读取检查点", i, nrow(metadata), exposure_id_i))
      return(cached)
    }
  }

  message(sprintf("[clump %d/%d] %s | %s", i, nrow(metadata), exposure_id_i, exposure_name_i))

  result_i <- tryCatch(
    prepare_one_exposure(
      exposure_file = exposure_file_i,
      exposure_id = exposure_id_i,
      exposure_name = exposure_name_i,
      group = group_i
    ),
    error = function(e) {
      list(
        status = "error",
        qc = data.table(
          exposure_id = exposure_id_i,
          exposure_name = exposure_name_i,
          group = group_i,
          n_valid = NA_integer_,
          n_p_threshold = NA_integer_,
          n_F_ge_10 = NA_integer_,
          n_MHC_removed = NA_integer_,
          n_after_MHC = NA_integer_,
          n_after_clump = NA_integer_
        ),
        exposure_dat = NULL,
        clump_result = NULL,
        error_message = conditionMessage(e)
      )
    }
  )

  status_i <- data.table(
    exposure_id = exposure_id_i,
    exposure_name = exposure_name_i,
    stage = "clumping",
    status = result_i$status,
    error_message = if (!is.null(result_i$error_message)) {
      result_i$error_message
    } else {
      NA_character_
    }
  )

  success_i <- identical(result_i$status, "success")

  if (success_i) {

    # Exposure对象保存为RDS，第三阶段直接读取，不把大型对象长期保存在主进程内存中。
    saveRDS(result_i$exposure_dat, exposure_cache_file)

    fwrite(
      as.data.table(result_i$exposure_dat),
      file.path(dir_clump, paste0(exposure_id_i, "_clumped_exposure.tsv")),
      sep = "\t"
    )

    fwrite(
      as.data.table(result_i$clump_result),
      file.path(dir_clump, paste0(exposure_id_i, "_clump_result.tsv")),
      sep = "\t"
    )
  }

  task_result <- list(
    exposure_id = exposure_id_i,
    success = success_i,
    qc = result_i$qc,
    status = status_i
  )

  # 使用临时文件加原子重命名，避免任务被中断时留下半写入RDS。
  task_tmp <- paste0(task_cache_file, ".tmp.", Sys.getpid())
  saveRDS(task_result, task_tmp)
  file.rename(task_tmp, task_cache_file)

  task_result
}

# mc.preschedule=FALSE使731个任务动态分配，避免少数大文件拖住固定worker。
clump_task_results <- parallel::mclapply(
  X = seq_len(nrow(metadata)),
  FUN = run_one_clump_task,
  mc.cores = min(clump_workers, nrow(metadata)),
  mc.preschedule = FALSE,
  mc.set.seed = FALSE,
  mc.cleanup = TRUE
)

clump_qc <- rbindlist(
  lapply(clump_task_results, function(x) x$qc),
  fill = TRUE
)

status_list <- lapply(
  clump_task_results,
  function(x) x$status
)

successful_ids <- vapply(
  clump_task_results,
  function(x) if (isTRUE(x$success)) x$exposure_id else NA_character_,
  character(1)
)

successful_ids <- successful_ids[!is.na(successful_ids)]

fwrite(
  clump_qc,
  file.path(dir_summary, "clumping_QC_summary.tsv"),
  sep = "\t"
)

if (length(successful_ids) == 0) {
  stop("没有Exposure成功完成clumping。")
}


# ==============================================================================
# 13. 第二阶段：汇总全部SNP，一次性扫描FinnGen
# ==============================================================================

message("阶段 2/3：汇总全部clumped SNP，并流式提取FinnGen outcome。")

exposure_cache_files <- setNames(
  file.path(dir_cache_clump, paste0(successful_ids, "_exposure_dat.rds")),
  successful_ids
)

all_clumped_snps <- unique(
  unlist(
    lapply(
      exposure_cache_files,
      function(f) as.character(readRDS(f)$SNP)
    ),
    use.names = FALSE
  )
)

all_clumped_file <- file.path(dir_summary, "all_clumped_rsids.txt")
outcome_subset_file <- file.path(
  result_root,
  "FinnGen_K11_COELIAC_target_SNPs.tsv"
)

# 只有目标SNP集合完全相同，才复用上次提取的FinnGen子集。
reuse_outcome_subset <- FALSE

if (resume_completed &&
    file.exists(all_clumped_file) &&
    file.exists(outcome_subset_file)) {

  old_targets <- unique(trimws(readLines(all_clumped_file, warn = FALSE)))
  old_targets <- old_targets[old_targets != ""]

  reuse_outcome_subset <- setequal(old_targets, all_clumped_snps)
}

writeLines(all_clumped_snps, all_clumped_file)

if (reuse_outcome_subset) {
  message("目标SNP集合未改变，复用已经提取的FinnGen子集。")
  outcome_subset <- fread(outcome_subset_file)
} else {
  outcome_subset <- extract_finngen_target_snps(
    outcome_file = outcome_file,
    target_snps = all_clumped_snps,
    output_file = outcome_subset_file
  )
}

matched_outcome_snps <- unique(trimws(as.character(outcome_subset$rsids)))
missing_outcome_snps <- setdiff(all_clumped_snps, matched_outcome_snps)

writeLines(
  missing_outcome_snps,
  file.path(dir_log, "missing_in_FinnGen_rsids.txt")
)

message("全部唯一clumped SNP数：", length(all_clumped_snps))
message("FinnGen匹配SNP数：", length(matched_outcome_snps))
message("FinnGen未匹配SNP数：", length(missing_outcome_snps))

outcome_dat <- prepare_outcome(outcome_subset)


# ==============================================================================
# 14. 第三阶段：并行harmonisation、MR和敏感性分析
# ==============================================================================

message("阶段 3/3：并行运行harmonisation、MR和敏感性分析。")
message("MR并发数：", mr_workers)

run_one_mr_task <- function(i) {

  data.table::setDTthreads(1L)

  exposure_id_i <- successful_ids[i]
  meta_i <- metadata[exposure_id == exposure_id_i][1]

  exposure_cache_file <- file.path(
    dir_cache_clump,
    paste0(exposure_id_i, "_exposure_dat.rds")
  )

  mr_cache_file <- file.path(
    dir_cache_mr,
    paste0(exposure_id_i, "_MR_task_result.rds")
  )

  # 读取已完成检查点。若上次为error且rerun_error_tasks=TRUE，则重新运行。
  if (resume_completed && file.exists(mr_cache_file)) {

    cached <- tryCatch(readRDS(mr_cache_file), error = function(e) NULL)

    if (!is.null(cached)) {
      cached_status <- cached$summary$status[1]

      if (!(rerun_error_tasks && identical(cached_status, "error"))) {
        message(sprintf("[MR %d/%d] %s：读取检查点", i, length(successful_ids), exposure_id_i))
        return(cached)
      }
    }
  }

  message(sprintf("[MR %d/%d] %s | %s", i, length(successful_ids), exposure_id_i, meta_i$exposure_name))

  # 每个Exposure使用固定但不同的随机种子，确保MR-PRESSO可重复且不受调度顺序影响。
  set.seed(parallel_seed + i)

  exposure_dat_i <- readRDS(exposure_cache_file)

  mr_i <- tryCatch(
    run_one_mr(
      exposure_dat = exposure_dat_i,
      outcome_dat = outcome_dat,
      exposure_id = exposure_id_i,
      exposure_name = meta_i$exposure_name,
      group = meta_i$group
    ),
    error = function(e) {
      empty_mr_summary(
        exposure_id = exposure_id_i,
        exposure_name = meta_i$exposure_name,
        group = meta_i$group,
        status = "error",
        error_message = conditionMessage(e)
      )
    }
  )

  status_i <- data.table(
    exposure_id = exposure_id_i,
    exposure_name = meta_i$exposure_name,
    stage = "MR",
    status = mr_i$status[1],
    error_message = mr_i$error_message[1]
  )

  task_result <- list(
    summary = mr_i,
    status = status_i
  )

  mr_tmp <- paste0(mr_cache_file, ".tmp.", Sys.getpid())
  saveRDS(task_result, mr_tmp)
  file.rename(mr_tmp, mr_cache_file)

  # 每个Exposure保存独立TSV，主脚本中断后也不会丢失已完成结果。
  fwrite(
    mr_i,
    file.path(dir_cache_mr, paste0(exposure_id_i, "_MR_summary.tsv")),
    sep = "\t"
  )

  task_result
}

RNGkind("L'Ecuyer-CMRG")
set.seed(parallel_seed)

mr_task_results <- parallel::mclapply(
  X = seq_along(successful_ids),
  FUN = run_one_mr_task,
  mc.cores = min(mr_workers, length(successful_ids)),
  mc.preschedule = FALSE,
  mc.set.seed = TRUE,
  mc.cleanup = TRUE
)

mr_summary_list <- lapply(
  mr_task_results,
  function(x) x$summary
)

status_list <- c(
  status_list,
  lapply(mr_task_results, function(x) x$status)
)

# 汇总阶段性结果。真正的断点续跑依赖每个Exposure独立的RDS/TSV检查点。
interim <- rbindlist(mr_summary_list, fill = TRUE)

fwrite(
  interim,
  file.path(dir_summary, "MR_summary_interim.tsv"),
  sep = "\t"
)

# ==============================================================================
# 15. 汇总所有结果并进行 FDR 校正
# ==============================================================================

mr_summary <- rbindlist(mr_summary_list, fill = TRUE)

final_summary <- merge(
  clump_qc,
  mr_summary,
  by = c("exposure_id", "exposure_name", "group"),
  all = TRUE
)

# 将extended metadata中的panel、trait_type等注释合并进最终结果。
optional_meta_cols <- intersect(
  c("panel", "trait_type", "analysis_tier", "statistical_name"),
  names(metadata)
)

if (length(optional_meta_cols) > 0) {
  meta_annotation <- unique(
    metadata[, c(
      "exposure_id", "exposure_name", "group", optional_meta_cols
    ), with = FALSE]
  )

  final_summary <- merge(
    meta_annotation,
    final_summary,
    by = c("exposure_id", "exposure_name", "group"),
    all.y = TRUE
  )
}

# 单 SNP 使用 Wald ratio；2个及以上 SNP 使用 IVW 为主分析
final_summary[, primary_method := fifelse(
  !is.na(ivw_p),
  "Inverse variance weighted",
  fifelse(!is.na(wald_p), "Wald ratio", NA_character_)
)]

final_summary[, primary_beta := fifelse(
  primary_method == "Inverse variance weighted",
  ivw_beta,
  fifelse(primary_method == "Wald ratio", wald_beta, NA_real_)
)]

final_summary[, primary_p := fifelse(
  primary_method == "Inverse variance weighted",
  ivw_p,
  fifelse(primary_method == "Wald ratio", wald_p, NA_real_)
)]

final_summary[, primary_or := fifelse(
  primary_method == "Inverse variance weighted",
  ivw_or,
  fifelse(primary_method == "Wald ratio", wald_or, NA_real_)
)]

final_summary[, primary_lci := fifelse(
  primary_method == "Inverse variance weighted",
  ivw_lci,
  fifelse(primary_method == "Wald ratio", wald_lci, NA_real_)
)]

final_summary[, primary_uci := fifelse(
  primary_method == "Inverse variance weighted",
  ivw_uci,
  fifelse(primary_method == "Wald ratio", wald_uci, NA_real_)
)]

# 对所有纳入 Exposure 的主分析 P 值进行全局 BH-FDR
final_summary[, primary_fdr := p.adjust(primary_p, method = "BH")]

# 同时计算分组内FDR，作为补充结果
final_summary[, group_fdr := p.adjust(primary_p, method = "BH"), by = group]

# extended metadata存在panel列时，同时计算panel内FDR。
if ("panel" %in% names(final_summary)) {
  final_summary[, panel_fdr := p.adjust(primary_p, method = "BH"), by = panel]
}

# 常用结果标记
final_summary[, nominal_significant :=
                !is.na(primary_p) & primary_p < 0.05]

final_summary[, global_fdr_significant :=
                !is.na(primary_fdr) & primary_fdr < 0.10]

final_summary[, fdr_significant :=
                !is.na(panel_fdr) & panel_fdr < 0.10]

final_summary[, heterogeneity_flag :=
                !is.na(ivw_Q_p) & ivw_Q_p < 0.05]

final_summary[, egger_pleiotropy_flag :=
                !is.na(egger_intercept_p) &
                egger_intercept_p < 0.05]

final_summary[, wm_direction_consistent :=
                !is.na(ivw_beta) &
                !is.na(wm_beta) &
                sign(ivw_beta) == sign(wm_beta)]

final_summary[, egger_direction_consistent :=
                !is.na(ivw_beta) &
                !is.na(egger_beta) &
                sign(ivw_beta) == sign(egger_beta)]

# 标记 MHC 和 LD 参数，方便后续与其他参数组合合并比较
final_summary[, mhc_setting := ifelse(
  exclude_mhc,
  "MHC_excluded",
  "MHC_retained"
)]

final_summary[, clump_r2_setting := clump_r2]
final_summary[, clump_kb_setting := clump_kb]
final_summary[, p_threshold_setting := p_threshold]
final_summary[, f_threshold_setting := f_threshold]

setorder(final_summary, primary_fdr, primary_p)

fwrite(
  final_summary,
  file.path(dir_summary, "batch_MR_full_summary.tsv"),
  sep = "\t"
)

fwrite(
  final_summary[nominal_significant == TRUE],
  file.path(dir_summary, "batch_MR_nominal_significant.tsv"),
  sep = "\t"
)

fwrite(
  final_summary[fdr_significant == TRUE],
  file.path(dir_summary, "batch_MR_FDR_significant.tsv"),
  sep = "\t"
)

status_table <- rbindlist(status_list, fill = TRUE)

fwrite(
  status_table,
  file.path(dir_log, "batch_status_log.tsv"),
  sep = "\t"
)


# ==============================================================================
# 16. 运行结束摘要
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("批量 MR 分析完成\n")
cat("============================================================\n")
cat("元数据中纳入Exposure数：", nrow(metadata), "\n")
cat("clumping并发数：", clump_workers, "\n")
cat("MR并发数：", mr_workers, "\n")
cat("成功完成clumping数：", length(successful_ids), "\n")
cat(
  "获得主 MR P 值的 Exposure 数：",
  sum(!is.na(final_summary$primary_p)),
  "\n"
)
cat(
  "名义显著 Exposure 数：",
  sum(final_summary$nominal_significant, na.rm = TRUE),
  "\n"
)
cat(
  "FDR 显著 Exposure 数：",
  sum(final_summary$fdr_significant, na.rm = TRUE),
  "\n"
)
cat(
  "完整汇总表：",
  file.path(dir_summary, "batch_MR_full_summary.tsv"),
  "\n"
)
cat(
  "运行状态日志：",
  file.path(dir_log, "batch_status_log.tsv"),
  "\n"
)
cat("============================================================\n")
