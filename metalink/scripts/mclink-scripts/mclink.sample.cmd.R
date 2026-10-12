###
#' @Date: 2022-11-02 23:40:20
#' @LastEditors: LiuyangLi
#' @LastEditTime: 2024-07-12 14:03:49
#' @FilePath: /liliuyang/software/mcLink/mclink.sample.cmd.R
#' @Description: transfer the count to relative abundance

suppressMessages(library(optparse))
rm(list = ls(all.names = TRUE))
package = c("mclink","dplyr")
repo = "http://mirrors.tuna.tsinghua.edu.cn/CRAN/"
ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    # Try GitHub first for mclink package
    if ("mclink" %in% new.pkg) {
      tryCatch({
        if (!requireNamespace("remotes", quietly = TRUE)) {
          install.packages("remotes", repos = repo)
        }
        if (!require("devtools")) install.packages("devtools")
        devtools::install_github(paste0("LiuyangLee/", "mclink"))
        new.pkg <- new.pkg[new.pkg != "mclink"]
      }, error = function(e) {
        cat("GitHub install failed for mclink, falling back to CRAN\n")
      })
    }
    # Install remaining packages from CRAN
    if (length(new.pkg)) {
      install.packages(new.pkg, dependencies = TRUE, repos = repo)
    }
  }
  invisible(sapply(pkg, require, character.only = TRUE, quietly = TRUE))
  #cat(paste0("Loading packages: ", pkg, "\n"))
}
suppressMessages(ipak(package))

option_list <- list(
  make_option(c("--in_KO_pathway_info"), type="character", default=NULL, help="in_KO_pathway_info"),
  make_option(c("--in_KO_Sample_wide_RelativeAbundance"), type="character", default=NULL, help="in_KO_Sample_wide_RelativeAbundance: primary input matrix (RA scale); the Keep filter mask is derived from it"),
  make_option(c("--in_KO_Sample_wide_SumCoverage"), type="character", default=NULL, help="optional second abundance matrix (SumCoverage scale); filtered with the Keep mask derived from the primary input (--in_KO_Sample_wide_RelativeAbundance, RA scale); when table_feature=abundance, mclink is additionally run on this matrix with out_dir=<out_DIR>/SumCoverage"),
  make_option(c("--table_feature"), type="character", default=NULL, help="table_feature, 'abundance' or 'completeness'"),
  make_option(c("--plus_scale_method"), type="character", default=NULL, help="plus_scale_method, e.g., 'mean', 'max', or 'min'"),
  make_option(c("--comma_scale_method"), type="character", default=NULL, help="comma_scale_method, e.g., 'max' for presence/absence data, or 'sum' for abundance data"),
  make_option(c("--minsample"), type="character", default=NULL, help="minsample: minimum number of samples in which a KO must be detected; 0 or empty means no filtering"),
  make_option(c("--prevalence"), type="character", default=NULL, help="prevalence: minimum detection rate (0-1) of a KO across all samples; 0 or empty means no filtering"),
  make_option(c("--abundance"), type="character", default=NULL, help="abundance: minimum mean relative abundance of a KO across all samples; 0 or empty means no filtering"),
  make_option(c("--out_DIR"), type="character", default=NULL, help="output directory"),
  make_option(c("--split_by_pathway"), type="character", default=NULL, help="TURE or FALSE"),
  make_option(c("--n_cores"), type="character", default="1", help="n_cores: number of cores for mclink module-parallel processing (default: 1, sequential)"),
  make_option(c("--verbose"), action="store_true", default=FALSE, help="print messages")
)
opt <- parse_args(OptionParser(option_list=option_list))

in_KO_pathway_info <- opt$in_KO_pathway_info
in_KO_Sample_wide_RelativeAbundance <- opt$in_KO_Sample_wide_RelativeAbundance
in_KO_Sample_wide_SumCoverage <- opt$in_KO_Sample_wide_SumCoverage
table_feature <- opt$table_feature
plus_scale_method <- opt$plus_scale_method
comma_scale_method <- opt$comma_scale_method
minsample <- opt$minsample
prevalence <- opt$prevalence
abundance <- opt$abundance
out_DIR <- opt$out_DIR
split_by_pathway <- opt$split_by_pathway
n_cores <- opt$n_cores
verbose <- opt$verbose

# empty/invalid n_cores falls back to sequential processing
if (is.null(n_cores) || is.na(suppressWarnings(as.integer(n_cores))) || as.integer(n_cores) < 1) n_cores <- 1

# empty thresholds mean no KO filtering (equivalent to keeping all KOs)
if (is.null(minsample) || minsample == "") minsample <- 0
if (is.null(prevalence) || prevalence == "") prevalence <- 0
if (is.null(abundance) || abundance == "") abundance <- 0

suppressMessages(library(mclink))


ra_cal = function(data = MAGs_rpkm){
  # relative abundance Calculate
  # MAGs_rpkm: row,KO; col,samples.
  tax_df = data %>%
    {.[rowSums(., na.rm = TRUE)>0,, drop = F]} %>%
    {.[,colSums(., na.rm = TRUE)>0, drop = F]} %>%
    {as.data.frame(t(.),stringsAsFactors = F)} %>%
    {.*1e6/rowSums(., na.rm = TRUE)} %>%
    {as.data.frame(t(.),stringsAsFactors = F)}
  return(tax_df)
}

pmr_cal = function(data = MAGs_rpkm
                  ){
  # prevalence and mean relative abundance Calculation
  # MAGs_rpkm: row,KO; col,samples.
  MAGs_stat = data %>%
    {ra_cal(.)} %>%
    {data.frame(
      Min_sample = apply(., 1, function (x) {sum(ifelse(x>0, 1, 0))}),
      Prevalence = apply(., 1, function (x) {sum(ifelse(x>0, 1, 0)/dim(.)[2])}),
      Mean_RA = apply(., 1, mean)
    )}
  return(MAGs_stat)
}

pmr_sub = function(data = MAGs_rpkm,
                   pre= 0.25,
                   ra = 0.000001,
                   min_sample = 0
                  ){
  # subset the result from pmr_cal
  # MAGs_rpkm: row,KO; col,samples.
  MAGs_stat = pmr_cal(data) %>%
    {mutate(., Keep = (.$Prevalence >= as.numeric(pre)) &
      (.$Mean_RA >= as.numeric(ra)) &
      (.$Min_sample >= as.numeric(min_sample))
      )}
  print(paste0('INFO: ', sum(MAGs_stat$Keep), ' KOs are kept'))
  data.scale = data %>%
    {.[rownames(.) %in% rownames(MAGs_stat)[MAGs_stat$Keep %in% 'TRUE'], , drop = FALSE]}
  return(data.scale)
}

#data(KO_pathway_ref)
#data(KO_Sample_wide)
#table_feature = "completeness"
#plus_scale_method = "min"
#comma_scale_method = "max"
#out_dir = "output_completeness_mean_max"
#split_by_pathway = TRUE

para_KO_pathway_ref = read.csv(in_KO_pathway_info, sep='\t',header = TRUE, check.names = F)
para_KO_Sample_wide = data.table::fread(in_KO_Sample_wide_RelativeAbundance, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
para_plus_scale_method = as.character(plus_scale_method) # mean, max, min
para_comma_scale_method = as.character(comma_scale_method) # max, sum
para_table_feature = as.character(table_feature) # abundance, completeness"
para_out_dir = as.character(out_DIR)
para_split_by_pathway = split_by_pathway

df_start <- para_KO_Sample_wide
rownames(df_start) <- df_start[[1]]
ko_col <- df_start[[1]]
n_na <- sum(is.na(as.matrix(df_start[, -1, drop = FALSE])))
if (n_na > 0) {
  cat(sprintf("[%s] WARNING: input table contains %d NA value(s); replaced with 0, please check the upstream data\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), n_na))
  df_start[, -1][is.na(df_start[, -1])] <- 0
}
cat(sprintf("[%s] Original table: %d KO (rows) with %d samples (columns)\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            nrow(df_start), ncol(df_start) - 1))
dt_filt <- df_start[, -1, drop = FALSE] %>%
  as.data.frame() %>%
  pmr_sub(., pre = prevalence, ra = abundance, min_sample = minsample)
rownames(dt_filt) <- ko_col[match(rownames(dt_filt), ko_col)]
df_final <- cbind(KO = rownames(dt_filt), dt_filt)
rownames(df_final) <- NULL
cat(sprintf("[%s] After filtering: %d KO (rows) with %d samples (columns)\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            nrow(df_final), ncol(df_final) - 1))

# abundance results go to the RelativeAbundance subdirectory; completeness results go to the Completeness subdirectory (shared, computed only once)
if (para_table_feature == "abundance") {
    para_out_dir_track <- file.path(para_out_dir, "RelativeAbundance")
} else {
    para_out_dir_track <- file.path(para_out_dir, "Completeness")
}
if (!dir.exists(para_out_dir_track)) dir.create(para_out_dir_track, recursive = TRUE)

mclink_args <- list(ref = para_KO_pathway_ref,
                    data = df_final,
                    table_feature = para_table_feature,
                    plus_scale_method = para_plus_scale_method,
                    comma_scale_method = para_comma_scale_method,
                    out_dir = para_out_dir_track,
                    split_by_pathway = para_split_by_pathway)
# n_cores only exists in mclink >= 1.3.0; older versions run sequentially
if ("n_cores" %in% names(formals(mclink::mclink))) {
    mclink_args$n_cores <- as.integer(n_cores)
} else {
    cat(sprintf("[%s] mclink version without n_cores support; running sequentially.\n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}
mc_list = do.call(mclink::mclink, mclink_args)

# optional second abundance matrix (SumCoverage scale): filtered with the Keep mask derived from the RA-scale primary input
if (para_table_feature == "abundance" && !is.null(in_KO_Sample_wide_SumCoverage) && in_KO_Sample_wide_SumCoverage != "") {
  para_KO_Sample_wide_SumCoverage = data.table::fread(in_KO_Sample_wide_SumCoverage, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
  df_cov <- para_KO_Sample_wide_SumCoverage
  rownames(df_cov) <- df_cov[[1]]
  n_na <- sum(is.na(as.matrix(df_cov[, -1, drop = FALSE])))
  if (n_na > 0) {
    cat(sprintf("[%s] WARNING: SumCoverage table contains %d NA value(s); replaced with 0, please check the upstream data\n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S"), n_na))
    df_cov[, -1][is.na(df_cov[, -1])] <- 0
  }
  cat(sprintf("[%s] SumCoverage table: %d KO (rows) with %d samples (columns)\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              nrow(df_cov), ncol(df_cov) - 1))
  # apply the RA-derived Keep mask to the SumCoverage matrix (subset by KO row names)
  kept_kos <- df_final$KO
  n_missing_kos <- sum(!kept_kos %in% rownames(df_cov))
  if (n_missing_kos > 0) {
    cat(sprintf("[%s] WARNING: %d KO(s) kept by the RA mask are missing from the SumCoverage matrix; skipped\n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S"), n_missing_kos))
  }
  kept_kos <- kept_kos[kept_kos %in% rownames(df_cov)]
  kept_samples <- colnames(df_final)[-1]
  n_missing_samples <- sum(!kept_samples %in% colnames(df_cov))
  if (n_missing_samples > 0) {
    cat(sprintf("[%s] WARNING: %d sample column(s) of the primary input are missing from the SumCoverage matrix; skipped\n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%S"), n_missing_samples))
  }
  kept_samples <- kept_samples[kept_samples %in% colnames(df_cov)]
  dt_cov <- df_cov[kept_kos, kept_samples, drop = FALSE]
  df_cov_final <- cbind(KO = rownames(dt_cov), dt_cov)
  rownames(df_cov_final) <- NULL
  cat(sprintf("[%s] SumCoverage after filtering: %d KO (rows) with %d samples (columns)\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              nrow(df_cov_final), ncol(df_cov_final) - 1))
  para_out_dir_cov <- file.path(para_out_dir, "SumCoverage")
  if (!dir.exists(para_out_dir_cov)) dir.create(para_out_dir_cov, recursive = TRUE)
  mclink_args_cov <- mclink_args
  mclink_args_cov$data <- df_cov_final
  mclink_args_cov$out_dir <- para_out_dir_cov
  mc_list_cov = do.call(mclink::mclink, mclink_args_cov)
}
