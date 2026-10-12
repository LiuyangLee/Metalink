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

#checkv_quality = "Low-quality"
#source("C:\\Users\\Administrator\\Desktop\\Script\\Rscript\\R-20241204\\RLib\\R\\RLib.r")
#source("C:\\Users\\liliuyang\\Desktop\\Script\\Rscript\\R-20241204\\RLib\\R\\RLib.r")
#in_vMAGs_tax_detail = '../results/MetaG-virus-sample/10_summary/Y36_17/vMAG_summary_info.tsv'
#out_vMAGs_tax_detail = '../results/MetaG-virus-sample/10_summary/marmgCAMI2_2018.08.15_09.49.32_sample_8/vMAG_summary_info.filter.tsv'
script_path <- sub("--file=", "", commandArgs()[grepl("--file=", commandArgs())])
if (length(script_path) == 0 || script_path == "") {
  script_path <- sys.frame(1)$ofile
}
source(file.path(dirname(script_path), "RLib.r"))
in_vMAGs_tax_detail = NULL

option_list <- list(
  make_option(c("--in_vMAGs_tax_detail"), type="character", default=NULL, help="vMAG_summary_info.tsv"),
  make_option(c("--checkv_quality"), type="character", default=NULL, help="checkv quality filter: Complete, High-quality, Medium-quality, Low-quality, Not-determined (default: Low-quality)"),
  make_option(c("--out_vMAGs_tax_detail"), type="character", default=NULL, help="vMAG_summary_info.filter.tsv"),
  make_option(c("--verbose"), action="store_true", default=FALSE, help="print messages")
)
parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

in_vMAGs_tax_detail <- opt$in_vMAGs_tax_detail
checkv_quality <- opt$checkv_quality
out_vMAGs_tax_detail <- opt$out_vMAGs_tax_detail
verbose <- opt$verbose
#################################################################
########################    check     ###########################
#################################################################
if (is.null(checkv_quality) || checkv_quality == "") {
  checkv_quality <- "Low-quality"
}
valid_qualities <- c("Complete", "High-quality", "Medium-quality", "Low-quality", "Not-determined")
if (!checkv_quality %in% valid_qualities) {
  stop("--checkv_quality must be one of: Complete, High-quality, Medium-quality, Low-quality, Not-determined")
}

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --checkv_quality:',checkv_quality,'\n')

##########################################################################
########################    virus MAGs taxononmy   ###########################
##########################################################################

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract virus MAGs taxononmy from:',in_vMAGs_tax_detail,'\n')
virus_tax_detail = fread(in_vMAGs_tax_detail,  sep = "\t", #select = c(1,2,4),
                         header = T, data.table = T, fill = T)
colnames(virus_tax_detail)[1] = 'genome'
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total virus MAGs extracted:',nrow(virus_tax_detail), '\n')

quality_levels <- c("Complete"=5, "High-quality"=4, "Medium-quality"=3, "Low-quality"=2, "Not-determined"=1)
if (!"checkv_quality" %in% colnames(virus_tax_detail)) {
  stop(paste0("ERROR: column 'checkv_quality' not found in ", in_vMAGs_tax_detail))
}
n_bad_quality <- sum(is.na(virus_tax_detail$checkv_quality) | !(virus_tax_detail$checkv_quality %in% names(quality_levels)))
if (n_bad_quality > 0) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING:',in_vMAGs_tax_detail,'contains',n_bad_quality,'NA or unrecognized value(s) in checkv_quality column; the corresponding rows were removed by quality filtering, please check the upstream data','\n')
}
target_level <- quality_levels[checkv_quality]
current_levels <- quality_levels[virus_tax_detail$checkv_quality]
virus_tax_detail <- virus_tax_detail[current_levels >= target_level]

# vMAGs are confirmed viruses by the identification chain; an empty superkingdom
# means "no close reference in the database", not "not a virus" -- promote the
# bare Domain (same convention as the combine_taxonomy_function_depth scripts)
if ("superkingdom" %in% colnames(virus_tax_detail)) {
  virus_tax_detail[is.na(superkingdom) | superkingdom == "", superkingdom := "d__Viruses"]
}

if (nrow(virus_tax_detail) > 0) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Viruses after quality filtering:',nrow(virus_tax_detail), '\n')
  virus_tax_detail %>%
    {fwrite(., out_vMAGs_tax_detail, row.names = FALSE, col.names = T,
            sep = '\t', quote = FALSE)}
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'No Viruses passed quality filtering!\n')
  file.create(out_vMAGs_tax_detail)
}
