###
#' @Date: 2023-06-14 16:48:43
#' @LastEditors: LiuyangLi
#' @LastEditTime: 2024-10-16 23:29:28
#' @FilePath: /liliuyang/R/Rscript.genome/combine_taxonomy_function_depth.vmcr.cmd.R
#' @Description:
###

rm(list = ls(all.names = TRUE))
packages=c("reshape2","dplyr","optparse","tidyverse","data.table")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  invisible(sapply(pkg, require, character.only = TRUE))
  #cat(paste0("Loading packages: ", pkg, "\n"))
}
suppressMessages(ipak(packages))

#MAGs_min_completeness1 = 99
#MAGs_max_contamination1 = 0.1
#MAGs_min_completeness2 = 99
#MAGs_max_contamination2 = 0.1
#MAGs_taxonomy_classify_method = "BAT"
#source("C:\\Users\\Administrator\\Desktop\\Script\\Rscript\\R-20241204\\RLib\\R\\RLib.r")
##source("C:\\Users\\liliuyang\\Desktop\\Script\\Rscript\\R-20241204\\RLib\\R\\RLib.r")
#in_MAGs_tax_detail="../results/MetaG-MAGs-sample/09_summary/marmgCAMI2_2018.08.15_09.49.32_sample_8/mMAG_summary_info.tsv"
#out_MAGs_tax_detail = "out_MAGs_tax_detail.tsv"

script_path <- sub("--file=", "", commandArgs()[grepl("--file=", commandArgs())])
if (length(script_path) == 0 || script_path == "") {
  script_path <- sys.frame(1)$ofile
}
source(file.path(dirname(script_path), "RLib.r"))
option_list <- list(
  make_option(c("--MAGs_taxonomy_classify_method"), type="character", default=NULL, help="gtdbtk or BAT"),
  make_option(c("--in_MAGs_tax_detail"), type="character", default=NULL, help="in_MAGs_tax_detail.tsv"),
  make_option(c("--MAGs_min_completeness1"), type="double", default=0, help="Minimum completeness threshold (%) for MAGs quality filtering (default: 0)"),
  make_option(c("--MAGs_max_contamination1"), type="double", default=100, help="Maximum contamination threshold (%) for MAGs quality filtering (default: 100)"),
  make_option(c("--MAGs_min_completeness2"), type="double", default=0, help="Minimum completeness threshold (%) for MAGs quality filtering (default: 0)"),
  make_option(c("--MAGs_max_contamination2"), type="double", default=100, help="Maximum contamination threshold (%) for MAGs quality filtering (default: 100)"),
  make_option(c("--out_MAGs_tax_detail"), type="character", default=NULL, help="out_MAGs_tax_detail.tsv"),
  make_option(c("--verbose"), action="store_true", default=FALSE, help="print messages")
)
opt <- parse_args(OptionParser(option_list=option_list))
MAGs_taxonomy_classify_method <- opt$MAGs_taxonomy_classify_method
in_MAGs_tax_detail <- opt$in_MAGs_tax_detail
MAGs_min_completeness1 <- opt$MAGs_min_completeness1
MAGs_max_contamination1 <- opt$MAGs_max_contamination1
MAGs_min_completeness2 <- opt$MAGs_min_completeness2
MAGs_max_contamination2 <- opt$MAGs_max_contamination2
out_MAGs_tax_detail <- opt$out_MAGs_tax_detail
verbose <- opt$verbose
#################################################################
########################    check     ###########################
#################################################################
if (is.null(MAGs_taxonomy_classify_method) || !MAGs_taxonomy_classify_method %in% c("gtdbtk","BAT")) {
  stop("--MAGs_taxonomy_classify_method must be 'gtdbtk' or 'BAT'")
}
validate_MAGs_param <- function(param_value, param_name, min_val = 0, max_val = 100) {
  if (is.null(param_value) || is.na(param_value) || is.nan(param_value)) {
    stop(paste0("--", param_name, " must be a valid number"))
  }
  if (param_value < min_val || param_value > max_val) {
    stop(paste0("--", param_name, " must be between ", min_val, " and ", max_val))
  }
  return(param_value)
}
MAGs_min_completeness1 <- validate_MAGs_param(MAGs_min_completeness1, "MAGs_min_completeness1")
MAGs_max_contamination1 <- validate_MAGs_param(MAGs_max_contamination1, "MAGs_max_contamination1")
MAGs_min_completeness2 <- validate_MAGs_param(MAGs_min_completeness2, "MAGs_min_completeness2")
MAGs_max_contamination2 <- validate_MAGs_param(MAGs_max_contamination2, "MAGs_max_contamination2")

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --MAGs_taxonomy_classify_method:',MAGs_taxonomy_classify_method,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --MAGs_min_completeness1:',MAGs_min_completeness1,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --MAGs_max_contamination1:',MAGs_max_contamination1,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --MAGs_min_completeness2:',MAGs_min_completeness2,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --MAGs_max_contamination2:',MAGs_max_contamination2,'\n')
##########################################################################
########################     MAG taxononmy     ###########################
##########################################################################

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract MAGs taxononmy from:',in_MAGs_tax_detail,'\n')
MAG_tax_detail = fread(in_MAGs_tax_detail,  sep = "\t", header = T, data.table = F, fill = T)
colnames(MAG_tax_detail)[1] = "genome"
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total MAGs extracted:',nrow(MAG_tax_detail), '\n')

quality_cols <- c("Completeness1","Contamination1","Completeness2","Contamination2")
missing_quality_cols <- setdiff(quality_cols, colnames(MAG_tax_detail))
if (length(missing_quality_cols) > 0) {
  stop(paste0("ERROR: column(s) not found in ", in_MAGs_tax_detail, ": ", paste(missing_quality_cols, collapse = ", ")))
}
n_na_quality <- sum(is.na(as.matrix(MAG_tax_detail[, quality_cols, drop = FALSE])))
if (n_na_quality > 0) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING:',in_MAGs_tax_detail,'contains',n_na_quality,'NA value(s) in Completeness/Contamination columns; the corresponding MAGs will be removed by quality filtering, please check the upstream data','\n')
}

MAG_tax_detail <- MAG_tax_detail %>%
  filter(Completeness1 >= MAGs_min_completeness1,
         Contamination1 <= MAGs_max_contamination1,
         Completeness2 >= MAGs_min_completeness2,
         Contamination2 <= MAGs_max_contamination2)

if (nrow(MAG_tax_detail) > 0) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'MAGs after quality filtering:',nrow(MAG_tax_detail), '\n')
  MAG_tax_detail %>%
    {fwrite(., out_MAGs_tax_detail, row.names = FALSE, col.names = T,
            sep = '\t', quote = FALSE)}
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'No MAGs passed quality filtering!\n')
  file.create(out_MAGs_tax_detail)
}
