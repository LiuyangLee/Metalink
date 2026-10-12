###
#' @Date: 2023-06-14 16:48:43
#' @LastEditors: LiuyangLi
#' @LastEditTime: 2026-09-07
#' @FilePath: /liliuyang/R/Rscript.genome/combine_taxonomy_function_depth.vc.cmd.R
#' @Description:
###

rm(list = ls(all.names = TRUE))
packages=c("dplyr","optparse","stringr","data.table")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  invisible(sapply(pkg, require, character.only = TRUE))
  #cat(paste0("Loading packages: ", pkg, "\n"))
}
suppressMessages(ipak(packages))

#CAMI_format_output = T
#sample_name = 'Y36_17'
#KO_mode='ignore_multiple_KOs'
#contigs_taxonomy_classify_method = "CAT_pack"
#min_KO_num_in_one_taxon = 100
#minORFLen=50
#source("C:\\Users\\Administrator\\Desktop\\Script\\Rscript\\R-20241204\\RLib\\R\\RLib.r")
#in_read_length = "../results/MetaG-base-sample/01_trim_stats/Y36_17_trimmed_1.fastq.stats"
#in_taxid_tax="../taxid/taxid_tax_lineage.tsv.gz"
#in_taxid_fun="../taxid/taxid_fun_lineage.tsv.gz"
#in_contig_fun_detail = "../results/MetaG-contigs-sample/02_function/emapper/annotations/Y36_17.emapper.annotations.gz"
#in_contig_ffn_length = "../results/MetaG-base-sample/07_contigs/contigs_ffn_length/Y36_17"
#in_contig_depth = "../results/MetaG-base-sample/07_contigs/contigs_filter_depth/Y36_17"
#in_contig_tax_detail="../results/MetaG-contigs-sample/01_taxonomy/CAT_pack/Y36_17.CAT_pack.contig2classification.tsv"
#in_contigs_ORF_loci = '../results/MetaG-base-sample/07_contigs/contigs_ORF_loci/Y36_17'
#in_vMAGs_tax_detail = '../results/MetaG-virus-sample/10_summary/Y36_17/vMAG_summary_info.tsv'
#in_vMAGs_ffn_length = '../results/MetaG-virus-sample/10_summary/Y36_17/Virus_annotation_results.tsv'
#in_vMAGs_contig_length = "../results/MetaG-virus-sample/10_summary/Y36_17/final_virus.contig_length.txt"
#detail = "all.test.detail.tsv"
#out_Taxa_abundance = "test.detail.2.tsv"
#out_vc_consistency = "test.vc_consistency.tsv"
#out_CAMI_format = "test.out_CAMI_format.tsv"

script_path <- sub("--file=", "", commandArgs()[grepl("--file=", commandArgs())])
if (length(script_path) == 0 || script_path == "") {
  script_path <- sys.frame(1)$ofile
}
source(file.path(dirname(script_path), "RLib.r"))
# GForce-optimized gsum rounds differently from base R sum() (which dplyr uses
# per group); fall back to base evaluation so grouped sums are bit-identical to
# the dplyr implementation
options(datatable.optimize = 0L)
in_contig_fun_detail = NULL
in_contig_ffn_length = NULL
in_contig_depth = NULL
in_contig_tax_detail=NULL
detail = NULL
out_Taxa_abundance = NULL
in_contigs_ORF_loci = NULL
in_vMAGs_tax_detail = NULL
in_vMAGs_ffn_length = NULL
in_vMAGs_contig_length = NULL
# ===================== optparse parameter definition =====================
option_list <- list(
  make_option("--in_taxid_tax", type="character", default=NULL,
              help="taxid_tax table, e.g., taxid_tax.tsv"),
  make_option("--in_taxid_fun", type="character", default=NULL,
              help="taxid_tax table, e.g., taxid_fun_scale.tsv"),
  make_option("--contigs_taxonomy_classify_method", type="character", default=NULL,
              help="kraken2 or CAT_pack"),
  make_option("--contigs_function_classify_method", type="character", default=NULL,
              help="kraken2 or emapper"),
  make_option("--in_contig_tax_detail", type="character", default=NULL,
              help="in_contig_tax_detail.tsv"),
  make_option("--in_contig_fun_detail", type="character", default=NULL,
              help="in_contigs.emapper.annotations.gz"),
  make_option("--in_read_length", type="character", default=NULL,
              help="sample_trimmed_1.fastq.stats"),
  make_option("--in_contig_ffn_length", type="character", default=NULL,
              help="in_contig_ffn_length.tsv"),
  make_option("--in_contig_depth", type="character", default=NULL,
              help="in_contig_depth.tsv"),
  make_option("--in_contigs_ORF_loci", type="character", default=NULL,
              help="contigs_ORF_loci.tsv"),
  make_option("--in_vMAGs_tax_detail", type="character", default=NULL,
              help="vMAG_summary_info.tsv"),
  make_option("--in_vMAGs_ffn_length", type="character", default=NULL,
              help="Virus_annotation_results.tsv"),
  make_option("--in_vMAGs_contig_length", type="character", default=NULL,
              help="final_virus.contig_length.txt"),
  make_option("--sample_name", type="character", default=NULL,
              help="sample_name, e.g., 'SRR000001' (default: test)"),
  make_option("--KO_mode", type="character", default="ignore_multiple_KOs",
              help="KO handling mode: 'keep_multiple_KOs', 'ignore_multiple_KOs', or 'split_multiple_KOs' (default: ignore_multiple_KOs)"),
  make_option("--min_KO_num_in_one_taxon", type="numeric", default=0,
              help="min_KO_num_in_one_taxon: min number of KO genes (rows) a taxon (column) must contain"),
  make_option("--minORFLen", type="numeric", default=0,
              help="Min threshhold for ORF length"),
  make_option("--summary_Tax", type="character", default=NULL,
              help="output summary_Tax"),
  make_option("--summary_ORF_KO", type="character", default=NULL,
              help="output summary_ORF_KO"),
  make_option("--summary_ORF_Tax_KO", type="character", default=NULL,
              help="output summary"),
  make_option("--summary_ORF_Source", type="character", default=NULL,
              help="output summary_ORF_Source"),
  make_option("--summary_ORF_Tax_Source", type="character", default=NULL,
              help="output summary"),
  make_option("--out_Taxa_abundance", type="character", default=NULL,
              help="out_Taxa_abundance"),
  make_option("--out_KO_eggnog_depth", type="character", default=NULL,
              help="output out_KO_eggnog_depth"),
  make_option("--out_KO_Taxa_SumCoverage_wide", type="character", default=NULL,
              help="out_KO_Taxa_SumCoverage_wide"),
  make_option("--out_KO_Taxa_RelativeAbundance_wide", type="character", default=NULL,
              help="out_KO_Taxa_RelativeAbundance_wide"),
  make_option("--out_vc_consistency", type="character", default=NULL,
              help="out_vc_consistency"),
  make_option("--detail", type="character", default=NULL,
              help="output detail"),
  make_option("--CAMI_format_output", type="character", default="FALSE",
              help="Enable CAMI output format. Disabled when MAG mode uses GTDB-TK (default: FALSE)"),
  make_option("--out_CAMI_format", type="character", default=NULL,
              help="Output file path for CAMI format (only used when --CAMI_format_output is set)"),
  make_option("--rank", type="character", default="Domain,Phylum,Class,Order,Family,Genus,Species",
              help="Comma-separated taxonomic ranks to output KO_Taxa tables for (default: Domain,Phylum,Class,Order,Family,Genus,Species)"),
  make_option("--output_detail", type="character", default="FALSE",
              help="Whether to output the detail file (default: TRUE)"),
  make_option("--verbose", action="store_true", default=FALSE,
              help="print messages"),
  make_option("--out_mode_flag", type="character", default=NULL,
              help="Output file path to write the execution mode flag (default: NULL)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Assign to variables with same names as before for minimal code change
in_taxid_tax <- opt$in_taxid_tax
in_taxid_fun <- opt$in_taxid_fun
contigs_taxonomy_classify_method <- opt$contigs_taxonomy_classify_method
contigs_function_classify_method <- opt$contigs_function_classify_method
in_contig_tax_detail <- opt$in_contig_tax_detail
in_contig_fun_detail <- opt$in_contig_fun_detail
in_read_length <- opt$in_read_length
in_contig_ffn_length <- opt$in_contig_ffn_length
in_contig_depth <- opt$in_contig_depth
in_contigs_ORF_loci <- opt$in_contigs_ORF_loci
in_vMAGs_tax_detail <- opt$in_vMAGs_tax_detail
in_vMAGs_ffn_length <- opt$in_vMAGs_ffn_length
in_vMAGs_contig_length <- opt$in_vMAGs_contig_length
sample_name <- opt$sample_name
KO_mode <- opt$KO_mode
min_KO_num_in_one_taxon <- opt$min_KO_num_in_one_taxon
minORFLen <- opt$minORFLen
summary_Tax <- opt$summary_Tax
summary_ORF_KO <- opt$summary_ORF_KO
summary_ORF_Tax_KO <- opt$summary_ORF_Tax_KO
summary_ORF_Source <- opt$summary_ORF_Source
summary_ORF_Tax_Source <- opt$summary_ORF_Tax_Source
out_Taxa_abundance <- opt$out_Taxa_abundance
out_KO_eggnog_depth <- opt$out_KO_eggnog_depth
out_KO_Taxa_SumCoverage_wide <- opt$out_KO_Taxa_SumCoverage_wide
out_KO_Taxa_RelativeAbundance_wide <- opt$out_KO_Taxa_RelativeAbundance_wide
out_vc_consistency <- opt$out_vc_consistency
detail <- opt$detail
# make_KO_Taxa_tables() writes the per-rank long tables into dirname(detail);
# the directory is no longer auto-created via a declared output
if (!is.null(detail)) dir.create(dirname(detail), showWarnings = FALSE, recursive = TRUE)
CAMI_format_output <- opt$CAMI_format_output
out_CAMI_format <- opt$out_CAMI_format
rank <- opt$rank
output_detail <- opt$output_detail
verbose <- opt$verbose
out_mode_flag <- opt$out_mode_flag
#################################################################
########################    check     ###########################
#################################################################
if (is.null(sample_name)) {
  stop("--sample_name must be specified")
}

if (is.null(KO_mode) || !KO_mode %in% c("keep_multiple_KOs","ignore_multiple_KOs","split_multiple_KOs")) {
  stop("--KO_mode must be 'keep_multiple_KOs', 'ignore_multiple_KOs', or 'split_multiple_KOs'")
}

if (is.null(contigs_taxonomy_classify_method) || !contigs_taxonomy_classify_method %in% c("kraken2","CAT_pack")) {
  stop("--contigs_taxonomy_classify_method must be 'kraken2' or 'CAT_pack'")
}

if (is.null(contigs_function_classify_method) || !contigs_function_classify_method %in% c("kraken2","emapper")) {
  stop("--contigs_function_classify_method must be 'kraken2' or 'emapper'")
}

if (is.null(minORFLen) ||
    is.na(suppressWarnings(as.numeric(minORFLen))) ||
    is.nan(suppressWarnings(as.numeric(minORFLen)))) {
  stop("--minORFLen must be a valid number")
}
if (is.null(min_KO_num_in_one_taxon) ||
    is.na(suppressWarnings(as.numeric(min_KO_num_in_one_taxon))) ||
    is.nan(suppressWarnings(as.numeric(min_KO_num_in_one_taxon))) ||
    as.numeric(min_KO_num_in_one_taxon) < 0) {
  stop("--min_KO_num_in_one_taxon must be a valid non-negative number")
}

str2bool <- function(x, default = FALSE) {
  if (is.null(x) || x == "") return(default)
  x <- toupper(trimws(x))
  if (x %in% c("TRUE", "YES", "1", "T", "Y")) return(TRUE)
  if (x %in% c("FALSE", "NO", "0", "F", "N")) return(FALSE)
  stop(paste("Invalid boolean value:", x, ". Use TRUE/FALSE, yes/no, or 1/0"))
}
CAMI_format_output <- str2bool(CAMI_format_output, default = FALSE)
output_detail <- str2bool(output_detail, default = FALSE)

# Parse rank string into vector
rank_vec <- unlist(strsplit(rank, ","))
rank_vec <- trimws(rank_vec)
valid_ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
invalid_ranks <- setdiff(rank_vec, valid_ranks)
if (length(invalid_ranks) > 0) {
  stop(paste("Invalid rank(s):", paste(invalid_ranks, collapse=", "),
             ". Valid ranks:", paste(valid_ranks, collapse=", ")))
}

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --sample_name:',sample_name,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --min_KO_num_in_one_taxon:',min_KO_num_in_one_taxon,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --contigs_taxonomy_classify_method:',contigs_taxonomy_classify_method,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --contigs_function_classify_method:',contigs_function_classify_method,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --KO_mode:',KO_mode,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --minORFLen:',minORFLen,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --CAMI_format_output:',CAMI_format_output,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --output_detail:',output_detail,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --rank:',rank,'\n')

################################################################
####################     taxid_taxonomy      ######################
################################################################
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract mapping relationship of taxid-taxonomy from:',in_taxid_tax,'\n')
taxid_taxonomy = fread(in_taxid_tax, sep = "\t", colClasses = "character", col.names = c('taxid','Taxonomy','Lineage'), header = FALSE, data.table = T)

rank_cols    <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")
taxid_cols   <- c("taxid_Root","taxid_Domain", "taxid_Phylum", "taxid_Class", "taxid_Order", "taxid_Family", "taxid_Genus", "taxid_Species", "taxid_Strain")
lineage_cols <- c("Lineage_Root","Lineage_Domain", "Lineage_Phylum", "Lineage_Class", "Lineage_Order", "Lineage_Family", "Lineage_Genus", "Lineage_Species", "Lineage_Strain")

# Expand taxonomy strings into the rank/NewTax (and CAMI lineage) columns.
# The expansion is a pure function of Taxonomy (and Lineage), so it is applied
# once to the unique strings and mapped back by index instead of row by row.
expand_taxonomy_columns <- function(dt, id_col, CAMI_format_output) {
  # rename_GTDB_taxonomy_dataframe() converts NA to "" via data[is.na(data)] <- "";
  # do the same on the source columns before the unique mapping
  set(dt, which(is.na(dt[["Taxonomy"]])), "Taxonomy", "")
  if (CAMI_format_output) {
    set(dt, which(is.na(dt[["Lineage"]])), "Lineage", "")
    uniq <- unique(dt[, .(Taxonomy, Lineage)])
    uniq[, (taxid_cols) := split_lineage_fill(Lineage, n = 9)]
    uniq <- uniq %>%
      as.data.frame(.) %>%
      mutate(., Lineage_Strain = (paste0(.$taxid_Root,";",.$taxid_Domain,";",.$taxid_Phylum,";",.$taxid_Class,";",.$taxid_Order,";",.$taxid_Family,";",.$taxid_Genus,";",.$taxid_Species,";",.$taxid_Strain))) %>%
      mutate(., Lineage_Species= (paste0(.$taxid_Root,";",.$taxid_Domain,";",.$taxid_Phylum,";",.$taxid_Class,";",.$taxid_Order,";",.$taxid_Family,";",.$taxid_Genus,";",.$taxid_Species))) %>%
      mutate(., Lineage_Genus  = (paste0(.$taxid_Root,";",.$taxid_Domain,";",.$taxid_Phylum,";",.$taxid_Class,";",.$taxid_Order,";",.$taxid_Family,";",.$taxid_Genus))) %>%
      mutate(., Lineage_Family = (paste0(.$taxid_Root,";",.$taxid_Domain,";",.$taxid_Phylum,";",.$taxid_Class,";",.$taxid_Order,";",.$taxid_Family))) %>%
      mutate(., Lineage_Order  = (paste0(.$taxid_Root,";",.$taxid_Domain,";",.$taxid_Phylum,";",.$taxid_Class,";",.$taxid_Order))) %>%
      mutate(., Lineage_Class  = (paste0(.$taxid_Root,";",.$taxid_Domain,";",.$taxid_Phylum,";",.$taxid_Class))) %>%
      mutate(., Lineage_Phylum = (paste0(.$taxid_Root,";",.$taxid_Domain,";",.$taxid_Phylum))) %>%
      mutate(., Lineage_Domain = (paste0(.$taxid_Root,";",.$taxid_Domain))) %>%
      mutate(., Lineage_Root   = (paste0(.$taxid_Root))) %>%
      mutate(., Lineage_Strain = sub(";+$", "", .$Lineage_Strain)) %>%
      mutate(., Lineage_Species= sub(";+$", "", .$Lineage_Species)) %>%
      mutate(., Lineage_Genus  = sub(";+$", "", .$Lineage_Genus)) %>%
      mutate(., Lineage_Family = sub(";+$", "", .$Lineage_Family)) %>%
      mutate(., Lineage_Order  = sub(";+$", "", .$Lineage_Order)) %>%
      mutate(., Lineage_Class  = sub(";+$", "", .$Lineage_Class)) %>%
      mutate(., Lineage_Phylum = sub(";+$", "", .$Lineage_Phylum)) %>%
      mutate(., Lineage_Domain = sub(";+$", "", .$Lineage_Domain)) %>%
      mutate(., Lineage_Root   = sub(";+$", "", .$Lineage_Root))
    uniq <- as.data.table(uniq)
    uniq[, (rank_cols) := split_lineage_fill(Taxonomy, n = 8)]
    uniq <- uniq %>%
      as.data.frame(.) %>%
      rename_GTDB_taxonomy_dataframe(.) %>%
      rename_strain_dataframe(.) %>%
      mutate(., Strain = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus,";",.$Species,";",.$Strain))) %>%
      mutate(., Species= (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus,";",.$Species))) %>%
      mutate(., Genus  = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus))) %>%
      mutate(., Family = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family))) %>%
      mutate(., Order  = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order))) %>%
      mutate(., Class  = (paste0(.$Domain,";",.$Phylum,";",.$Class))) %>%
      mutate(., Phylum = (paste0(.$Domain,";",.$Phylum)))
    uniq <- as.data.table(uniq)
    idx <- match(paste(dt$Taxonomy, dt$Lineage, sep = "\x01"),
                 paste(uniq$Taxonomy, uniq$Lineage, sep = "\x01"))
    # only columns used downstream are mapped back (taxid_Domain..taxid_Strain
    # are just intermediates for the Lineage_* strings)
    for (cl in c("taxid_Root", lineage_cols, rank_cols, "NewTax"))
      set(dt, j = cl, value = uniq[[cl]][idx])
    setcolorder(dt, c(id_col, "taxid", "Lineage", "Taxonomy"))
  } else {
    uniq <- unique(dt[, .(Taxonomy)])
    uniq[, (rank_cols) := split_lineage_fill(Taxonomy, n = 8)]
    uniq <- uniq %>%
      as.data.frame(.) %>%
      rename_GTDB_taxonomy_dataframe(.) %>%
      rename_strain_dataframe(.) %>%
      mutate(., Strain = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus,";",.$Species,";",.$Strain))) %>%
      mutate(., Species= (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus,";",.$Species))) %>%
      mutate(., Genus  = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus))) %>%
      mutate(., Family = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family))) %>%
      mutate(., Order  = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order))) %>%
      mutate(., Class  = (paste0(.$Domain,";",.$Phylum,";",.$Class))) %>%
      mutate(., Phylum = (paste0(.$Domain,";",.$Phylum)))
    uniq <- as.data.table(uniq)
    idx <- match(dt$Taxonomy, uniq$Taxonomy)
    for (cl in c(rank_cols, "NewTax"))
      set(dt, j = cl, value = uniq[[cl]][idx])
    dt[, c("taxid", "Lineage") := NULL]
  }
  invisible(dt)
}

make_vc_depth <- function(data2, group_col, level, ReadLength) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Calculating depth at the level:',group_col,'\n')
  d <- data2[, .(Taxa_ContigNum = .N,
                 Taxa_SeqBases = sum(SeqBases),
                 Taxa_RelativeAbundance = sum(ContigRA_Bases)),
             keyby = c("sample", group_col)]
  d[, Taxa_ReadCount := round(Taxa_SeqBases / ReadLength)]
  d[, Level := level]
  setnames(d, 2L, "Taxonomy")
  d[]
}

# Same, grouped by (<rank>, <Lineage_rank>) for the CAMI profile tables.
make_vc_lineage_depth <- function(data2, group_cols, level, ReadLength, min_semicolons = NULL) {
  d <- data2[, .(Taxa_ContigNum = .N,
                 Taxa_SeqBases = sum(SeqBases),
                 Taxa_RelativeAbundance = sum(ContigRA_Bases)),
             keyby = group_cols]
  d[, Taxa_ReadCount := round(Taxa_SeqBases / ReadLength)]
  d[, Level := level]
  if (!is.null(min_semicolons)) d <- d[str_count(d[[group_cols[2]]], ";") >= min_semicolons]
  setnames(d, c(1L, 2L), c("Taxonomy", "Lineage"))
  d[]
}

# group counts plus an extra "All" level on the first group column, reproducing
# the original bind_rows(group_by(cols) summarise, group_by(cols[1]="All", ...)
# summarise) followed by per-group totals and proportions
summary_with_all <- function(dt, count_cols, count_name, total_name, total_by, prop_name = "Proportion") {
  rest_cols <- count_cols[-1]
  s1 <- dt[, .(N = .N), keyby = count_cols]
  s2 <- dt[, .(N = .N), keyby = rest_cols]
  s2[, (count_cols[1]) := "All"]
  setcolorder(s2, c(count_cols[1], count_cols[-1]))
  out <- rbindlist(list(s1, s2))
  setnames(out, "N", count_name)
  out[, (total_name) := sum(get(count_name)), by = total_by]
  out[, (prop_name) := get(count_name) / get(total_name)]
  out[]
}

##########################################################################
########################    contig taxononmy   ###########################
##########################################################################
if (contigs_taxonomy_classify_method == "CAT_pack"){
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract contigs taxononmy from:',in_contig_tax_detail,'\n')
  seqid_taxid_detail = fread(in_contig_tax_detail, sep = "\t", header = T, data.table = T, fill = T)
  colnames(seqid_taxid_detail) = c("sample_Contig", "classification","reason", "Lineage", "lineage_scores","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  seqid_taxid_detail[, Lineage := gsub("\\*", "", Lineage)]
  seqid_taxid_detail[, taxid := gsub(".*;", "", Lineage)]
  seqid_taxid_detail[is.na(taxid)| taxid == "", taxid := "0"]
  seqid_taxid_detail = seqid_taxid_detail %>%
    select(., c("sample_Contig", "taxid"))
} else if (contigs_taxonomy_classify_method == "kraken2"){
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract contigs taxononmy from:',in_contig_tax_detail,'\n')
  seqid_taxid_detail = fread(in_contig_tax_detail, sep = "\t", select = c(2,3), header = FALSE, data.table = T)
  colnames(seqid_taxid_detail) = c("sample_Contig", "taxid")
  seqid_taxid_detail$taxid = seqid_taxid_detail$taxid %>%
    gsub(".* \\(taxid |\\)", "", .)
}

contig_seqid_tax = seqid_taxid_detail %>%
  merge.data.table(., taxid_taxonomy, by = 'taxid', all.x = T)

if (CAMI_format_output) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Keep taxid and lineage column for contig-level annotation.\n')
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Do not keep taxid and lineage column for contig-level annotation.\n')
}
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total contigs extracted:',nrow(seqid_taxid_detail), '\n')

# split_lineage_fill + rename_* + lineage paste chain applied once per unique
# Taxonomy (and Lineage) string and mapped back by index
expand_taxonomy_columns(contig_seqid_tax, "sample_Contig", CAMI_format_output)

#table(contig_seqid_tax$Domain)
#table(contig_seqid_tax$NewTax)

##########################################################################
########################    virus MAGs taxononmy   ###########################
##########################################################################
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract virus MAGs taxononmy from:',in_vMAGs_tax_detail,'\n')
virus_tax_detail = fread(in_vMAGs_tax_detail,  sep = "\t", #select = c(1,2,4),
                         header = T, data.table = T, fill = T)
colnames(virus_tax_detail)[1] = 'genome'
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'The number of virus MAGs used for analysis',nrow(virus_tax_detail), '\n')

virus_tax_detail[, lineage := gsub("\\*", "", lineage)]
virus_tax_detail[, taxid := gsub(".*;", "", lineage)]
virus_tax_detail[is.na(taxid)| taxid == "", taxid := "0"]

virus_tax_detail = virus_tax_detail %>%
  merge.data.table(., taxid_taxonomy, by = 'taxid', all.x = T)

# vMAGs are confirmed viruses by the identification chain; a Taxonomy left
# unclassified means "no close reference in the database", not "not a virus" --
# promote the bare Domain so Domain-level tables count them as d__Viruses
virus_tax_detail[is.na(Taxonomy) | Taxonomy == "" | grepl("^d__unclassified", Taxonomy),
                 Taxonomy := "d__Viruses"]

if (CAMI_format_output) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Keep taxid and lineage column for vMAG-level annotation.\n')
  virus_seqid_tax = as.data.table(virus_tax_detail)
  virus_seqid_tax[, c("taxid_Root","taxid_Domain", "taxid_Phylum", "taxid_Class", "taxid_Order", "taxid_Family", "taxid_Genus", "taxid_Species", "taxid_Strain") := split_lineage_fill(Lineage, n = 9)]
  virus_seqid_tax = virus_seqid_tax %>%
    as.data.frame(.) %>%
    mutate(., Lineage_Strain = (paste0(.$taxid_Root,";",.$taxid_Domain,";",.$taxid_Phylum,";",.$taxid_Class,";",.$taxid_Order,";",.$taxid_Family,";",.$taxid_Genus,";",.$taxid_Species,";",.$taxid_Strain))) %>%
    mutate(., Lineage_Species= (paste0(.$taxid_Root,";",.$taxid_Domain,";",.$taxid_Phylum,";",.$taxid_Class,";",.$taxid_Order,";",.$taxid_Family,";",.$taxid_Genus,";",.$taxid_Species))) %>%
    mutate(., Lineage_Genus  = (paste0(.$taxid_Root,";",.$taxid_Domain,";",.$taxid_Phylum,";",.$taxid_Class,";",.$taxid_Order,";",.$taxid_Family,";",.$taxid_Genus))) %>%
    mutate(., Lineage_Family = (paste0(.$taxid_Root,";",.$taxid_Domain,";",.$taxid_Phylum,";",.$taxid_Class,";",.$taxid_Order,";",.$taxid_Family))) %>%
    mutate(., Lineage_Order  = (paste0(.$taxid_Root,";",.$taxid_Domain,";",.$taxid_Phylum,";",.$taxid_Class,";",.$taxid_Order))) %>%
    mutate(., Lineage_Class  = (paste0(.$taxid_Root,";",.$taxid_Domain,";",.$taxid_Phylum,";",.$taxid_Class))) %>%
    mutate(., Lineage_Phylum = (paste0(.$taxid_Root,";",.$taxid_Domain,";",.$taxid_Phylum))) %>%
    mutate(., Lineage_Domain = (paste0(.$taxid_Root,";",.$taxid_Domain))) %>%
    mutate(., Lineage_Root   = (paste0(.$taxid_Root))) %>%
    mutate(., Lineage_Strain = sub(";+$", "", .$Lineage_Strain)) %>%
    mutate(., Lineage_Species= sub(";+$", "", .$Lineage_Species)) %>%
    mutate(., Lineage_Genus  = sub(";+$", "", .$Lineage_Genus)) %>%
    mutate(., Lineage_Family = sub(";+$", "", .$Lineage_Family)) %>%
    mutate(., Lineage_Order  = sub(";+$", "", .$Lineage_Order)) %>%
    mutate(., Lineage_Class  = sub(";+$", "", .$Lineage_Class)) %>%
    mutate(., Lineage_Phylum = sub(";+$", "", .$Lineage_Phylum)) %>%
    mutate(., Lineage_Domain = sub(";+$", "", .$Lineage_Domain)) %>%
    mutate(., Lineage_Root   = sub(";+$", "", .$Lineage_Root))
  virus_seqid_tax = virus_seqid_tax %>%
    dplyr::select(., c("genome","taxid","Lineage","Taxonomy",
                       "taxid_Root","taxid_Domain", "taxid_Phylum", "taxid_Class", "taxid_Order", "taxid_Family", "taxid_Genus", "taxid_Species", "taxid_Strain",
                       "Lineage_Root","Lineage_Domain", "Lineage_Phylum", "Lineage_Class", "Lineage_Order", "Lineage_Family", "Lineage_Genus", "Lineage_Species", "Lineage_Strain")) %>%
    as.data.table(.)
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Do not keep taxid and lineage column for vMAG-level annotation.\n')
  virus_seqid_tax = virus_tax_detail %>%
    select(., c("genome", "taxid","Taxonomy"))
}


vMAGs_depth <- virus_tax_detail[
  , .(scaffold = unlist(tstrsplit(scaffold, ",", fixed = TRUE)))
  , by = .(genome, lytic_state, coverm_abundance)
]
vMAGs_depth[, sample_Contig := tstrsplit(scaffold, "\\|\\|", fixed = FALSE, type.convert = TRUE)[[1]]
][, sample_Contig := tstrsplit(sample_Contig, "\\|provirus_", fixed = FALSE, type.convert = TRUE)[[1]]
][, sample_Contig := tstrsplit(sample_Contig, "\\|vs2region_", fixed = FALSE, type.convert = TRUE)[[1]]]
vMAGs_depth[, c("sample", "virus_contig") := tstrsplit(scaffold, "---", fill = "", keep = 1:2)]
vMAGs_depth[, c("sample", "contig") := tstrsplit(sample_Contig, "---", fill = "", keep = 1:2)]
vMAGs_depth[, genome_Contig := paste0(genome, "---", contig)]

########################  vMAGs_contig_length  ###########################
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract virus contig length from:',in_vMAGs_contig_length,'\n')
vMAGs_contig_length = in_vMAGs_contig_length %>% {fread(., sep = '\t', header = T, data.table = T)}
colnames(vMAGs_contig_length)=c("scaffold","ContigLen")
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Virus contigs number:',nrow(vMAGs_contig_length), '\n')
vMAGs_contig_length <- vMAGs_contig_length[scaffold %in% vMAGs_depth$scaffold]

vMAGs_depth <- vMAGs_depth %>%
  merge(., vMAGs_contig_length, by = 'scaffold', all.x = T)

colnames(vMAGs_depth) = c("scaffold","genome", "lytic_state","coverm_abundance",
                          "sample_Contig","sample","virus_contig",
                          "contig","genome_Contig", "ContigLen")
vMAGs_depth[, SeqBases := ContigLen * coverm_abundance]
vMAGs_depth.agg <- vMAGs_depth[
  , .(
    scaffold = paste(scaffold, collapse = "; "),
    Integrated_Virus_ContigLen = sum(ContigLen),
    Integrated_Virus_SeqBases = sum(SeqBases)
  ),
  by = c("sample_Contig")
]

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total virus contigs extracted:',nrow(vMAGs_depth),'\n')
virus_seqid_tax[, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain") := split_lineage_fill(Taxonomy, n = 8)]

virus_seqid_tax = virus_seqid_tax %>%
  as.data.frame(.) %>%
  rename_GTDB_taxonomy_dataframe(.) %>%
  rename_strain_dataframe(.) %>%
  #dplyr::select(., c("genome","Taxonomy","NewTax","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")) %>%
  mutate(., Strain = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus,";",.$Species,";",.$Strain))) %>%
  mutate(., Species= (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus,";",.$Species))) %>%
  mutate(., Genus  = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus))) %>%
  mutate(., Family = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family))) %>%
  mutate(., Order  = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order))) %>%
  mutate(., Class  = (paste0(.$Domain,";",.$Phylum,";",.$Class))) %>%
  mutate(., Phylum = (paste0(.$Domain,";",.$Phylum)))

#table(virus_seqid_tax$Domain)
#table(virus_seqid_tax$NewTax)
#####

##########################################################################
########## contig seq id combine from vMAG and contig taxonomy ############
##########################################################################
####

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract ORF length of vMAGs from:',in_vMAGs_ffn_length,'\n')
vMAGs_ffn_infor = in_vMAGs_ffn_length %>% {fread(., sep = '\t', header = T, data.table = T)}
colnames(vMAGs_ffn_infor)[1:2] = c("genome","ORF")
vMAGs_ffn_infor[, ORFLen := `ORF end` - `ORF start` + 1]
n_na_orflen <- sum(is.na(vMAGs_ffn_infor$ORFLen))
if (n_na_orflen > 0) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING:',in_vMAGs_ffn_length,'contains',n_na_orflen,'NA value(s) in ORF start/ORF end; ORFLen replaced with 0, please check the upstream data','\n')
  vMAGs_ffn_infor[is.na(ORFLen), ORFLen := 0]
}

vMAGs_ffn_infor[
  , sample_Contig := tstrsplit(scaffold, "\\|\\|", fixed = FALSE, type.convert = TRUE)[[1]]
][
  , sample_Contig := tstrsplit(sample_Contig, "\\|provirus_", fixed = FALSE, type.convert = TRUE)[[1]]
][
  , sample_Contig := tstrsplit(sample_Contig, "\\|vs2region_", fixed = FALSE, type.convert = TRUE)[[1]]
]
vMAGs_ffn_infor[, c("sample", "contig") := tstrsplit(sample_Contig, "---", fill = "", keep = 1:2)]
vMAGs_ffn_infor[, c("sample", "ORF") := tstrsplit(ORF, "---", fill = "", keep = 1:2)]
vMAGs_ffn_infor[, qaccver := paste0(genome, "---", ORF)]
vMAGs_ffn_infor[, Source := "vMAG"]
vMAGs_ffn_infor <- vMAGs_ffn_infor[scaffold %in% vMAGs_depth$scaffold]

cols <- c("qaccver", "ORFLen", "sample_Contig", "sample","contig", "ORF", "Source",'scaffold')
vMAGs_ffn_length <- vMAGs_ffn_infor[, ..cols]

vMAGs_ffn_length = vMAGs_ffn_length %>%
  merge(., vMAGs_depth[,c('genome','scaffold')], by = 'scaffold', all.x = T) %>%
  mutate(genome_Contig = paste0(genome,'---', contig))
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total ORFs at the vMAGs level:',nrow(vMAGs_ffn_length), '\n')
#####
####

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract ORF length of contigs from:',in_contig_ffn_length,'\n')
Contigs_ffn_length = in_contig_ffn_length %>% {fread(., sep = '\t', header = T, data.table = T)}
colnames(Contigs_ffn_length) = c('qaccver','ORFLen')
n_na_orflen <- sum(is.na(Contigs_ffn_length$ORFLen))
if (n_na_orflen > 0) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING:',in_contig_ffn_length,'contains',n_na_orflen,'NA value(s) in ORFLen; replaced with 0, please check the upstream data','\n')
  Contigs_ffn_length[is.na(ORFLen), ORFLen := 0]
}
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total ORFs at the contig level:',nrow(Contigs_ffn_length), '\n')
Contigs_ffn_length[, c("sample_Contig") := tstrsplit(qaccver, "_\\d*$")]
Contigs_ffn_length[, c("sample", "contig") := tstrsplit(sample_Contig, "---", fill = "", keep = 1:2)]
Contigs_ffn_length[, c("genome", "ORF") := tstrsplit(qaccver, "---", fill = "", keep = 1:2)]
Contigs_ffn_length[, Source := "Contig"]

Contigs_ffn_length = Contigs_ffn_length %>%
  mutate(genome_Contig = paste0(genome,'---', contig))

Contigs_contig_genome = Contigs_ffn_length %>%
  mutate(genome_Contig = paste0(genome,'---', contig)) %>%
  select(c("sample_Contig", "genome", "genome_Contig")) %>%
  unique()

# 导入 contigs ORF 位置
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract ORF loci of contigs from:',in_contigs_ORF_loci,'\n')
contigs_ORF_loci = in_contigs_ORF_loci %>% {fread(., sep = '\t', header = T, data.table = T)}
colnames(contigs_ORF_loci)[1] = "qaccver"
contigs_ORF_loci[, c("sample_Contig") := tstrsplit(qaccver, "_\\d*$")]
contigs_ORF_loci[, c("sample", "contig") := tstrsplit(sample_Contig, "---", fill = "", keep = 1:2)]
contigs_ORF_loci[, c("genome", "ORF") := tstrsplit(qaccver, "---", fill = "", keep = 1:2)]

contigs_ORF_loci = contigs_ORF_loci %>%
  mutate(genome_Contig = paste0(genome,'---', contig)) %>%
  select(.,c("qaccver", "start","end", "strand"))


##################################
###### integrated and non_integrated #######
############################################
# c分成了3个部分，整合噬菌体，独立contigs被认为是病毒，独立contigs被认为是微生物
# c_ffn_length_virus 导出一个文件
# c中的整合噬菌体若要纳入ffn_length计算系统，
# 由于genomad调用的pyprodigal和virsorter2调用的prodigal是两个软件，产生的ORF编号和直接用的prodigal产生的宿主ORF编号不完全一致，所以用范围来排除
# 在这里按照病毒边界扣掉contig对应病毒区域的序列，让它们分别对应病毒分类和宿主分类，整合型噬菌体的ORF（进行了重新编号）和宿主的ORF这两套编号系统整合到一个ORF的文件中，
# 需要计算vMAGs中每个ORF的起点和终点
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Combine ORF length of virus...','\n')

Contigs_ffn_length_loci = Contigs_ffn_length %>%
  mutate(., ORF_name = paste0(sample, '---', ORF)) %>%
  merge(., contigs_ORF_loci, all.x = T, by = 'qaccver')
Contigs_ffn_length_loci[, `:=`(mORF_start = start, mORF_end = end)]

# 导入病毒 ORF 位置
vMAGs_ffn_infor[, `:=`(
  `ORF start` = as.integer(`ORF start`),
  `ORF end`   = as.integer(`ORF end`)
)]

# |vs2region_S_E scaffolds carry parent-contig coords in the ID, but their ORF
# coords are fragment-relative (prodigal ran on the extracted fragment); shift
# ORFs back to the parent-contig coordinate system (the host-ORF exclusion
# range below is taken from the embedded boundaries where present, so this
# shift only matters for scaffolds without embedded coordinates)
is_vs2region = grepl("\\|vs2region_\\d+_\\d+$", vMAGs_ffn_infor$scaffold)
if (any(is_vs2region)) {
  vs2_offset = as.integer(sub(".*\\|vs2region_(\\d+)_(\\d+)$", "\\1",
                              vMAGs_ffn_infor$scaffold[is_vs2region])) - 1L
  vMAGs_ffn_infor[is_vs2region, `:=`(`ORF start` = `ORF start` + vs2_offset,
                                        `ORF end`   = `ORF end`   + vs2_offset)]
}
vMAGs_scaffold_range <- vMAGs_ffn_infor[
  , .(vORF_start = min(`ORF start`, na.rm = TRUE),
      vORF_end   = max(`ORF end`,   na.rm = TRUE)),
  by = c("scaffold", "sample_Contig")
]
# where the scaffold ID embeds parent-contig boundaries (|provirus_S_E /
# |vs2region_S_E), use them as the viral region instead of the outermost viral
# ORF span: the host-side and viral-side prodigal runs often disagree at
# fragment ends, leaving host ORFs in the flank outside the outermost viral
# ORFs (same semantics as pmr_vm_mask_viral.py)
has_embedded_boundary = grepl("\\|(?:provirus|vs2region)_\\d+_\\d+$", vMAGs_scaffold_range$scaffold)
vMAGs_scaffold_range[has_embedded_boundary, `:=`(
  vORF_start = as.integer(sub(".*\\|(?:provirus|vs2region)_(\\d+)_(\\d+)$", "\\1", scaffold)),
  vORF_end   = as.integer(sub(".*\\|(?:provirus|vs2region)_(\\d+)_(\\d+)$", "\\2", scaffold)))]
vMAGs_scaffold_range[, `:=`(vORF_s = vORF_start, vORF_e = vORF_end)]

# 1. 宿主ORF 完全落于病毒区域，宿主ORF可在整个病毒区间内滑动
inside_virus <- vMAGs_scaffold_range[
  Contigs_ffn_length_loci,
  on = .(sample_Contig,
         vORF_s <= start,
         vORF_e >= end)
]
inside_virus[, `:=`(vORF_s = NULL, vORF_e = NULL)]
setcolorder(inside_virus,
            c("scaffold", "sample_Contig",
              "vORF_start", "vORF_end", "mORF_start", "mORF_end"))
inside_virus <- inside_virus[!is.na(scaffold)]

# 2. 宿主ORF 部分落于病毒区域，宿主ORF可在整个病毒区间内滑动以及允许部分区域在外面
host_virus_overlap <- vMAGs_scaffold_range[
  Contigs_ffn_length_loci,
  on = .(sample_Contig,
         vORF_s <= end,
         vORF_e >= start
  )
]
host_virus_overlap[, `:=`(vORF_s = NULL, vORF_e = NULL)]
setcolorder(host_virus_overlap,
            c("scaffold", "sample_Contig",
              "vORF_start", "vORF_end", "mORF_start", "mORF_end"))
host_virus_overlap <- host_virus_overlap[!is.na(scaffold)]

# 3. 宿主单条ORF 大于病毒区域
inside_host <- vMAGs_scaffold_range[
  Contigs_ffn_length_loci,
  on = .(sample_Contig,
         vORF_s >= start,
         vORF_e <= end)
]
inside_host[, `:=`(vORF_s = NULL, vORF_e = NULL)]
setcolorder(inside_host,
            c("scaffold", "sample_Contig",
              "vORF_start", "vORF_end", "mORF_start", "mORF_end"))
inside_host <- inside_host[!is.na(scaffold)]

# 4. 宿主和病毒边界存在分歧的列表  （按照病毒计算，从宿主中去除ORF）
# 先去掉宿主ORF完全在病毒内部的    （按照病毒计算，从宿主中去除ORF）
# 再去掉宿主ORF包裹了整个病毒区域的（按照病毒计算，从宿主中去除ORF）
host_viral_boundary_conflict.list =
  host_virus_overlap[!ORF_name %in% c(inside_virus$ORF_name,
                                      inside_host$ORF_name)]

host_viral_boundary_conflict = host_virus_overlap[qaccver %in% host_viral_boundary_conflict.list$qaccver]

colnames(host_viral_boundary_conflict) =
  c("virus_scaffold", "sample_Contig","virus_start","virus_end",
    "mORF_start","mORF_end",
    "mORF_qaccver","mORFLen","sample", "contig", "genome", "ORF", "Source",
    "genome_Contig", "ORF_name", "mORF_strand")
host_viral_boundary_conflict = host_viral_boundary_conflict %>%
  select(., c("virus_scaffold", "virus_start","virus_end",
              "mORF_start","mORF_end", "mORF_strand",
              "mORF_qaccver","Source"))

################################
cat(sprintf('[%s] Detected %d raw-contig ORFs overlapping viral regions in the raw contig annotation.\n',
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            nrow(host_virus_overlap)))
cat(sprintf('[%s] Detected %d viral ORFs from viral contigs; will use these virus-based predictions instead of overlapping raw-contig annotations.\n',
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            nrow(vMAGs_ffn_length)))


vc_ffn_length.hosts = Contigs_ffn_length %>%
  {.[!.$qaccver %in% host_virus_overlap$qaccver, ]} %>%
  select("qaccver","ORFLen", "sample_Contig", "sample", "contig", "genome", "ORF", "Source", "genome_Contig")
colnames(vc_ffn_length.hosts) = c("qaccver","ORFLen", "sample_Contig", "sample", "contig", "genome", "ORF", "Source", "genome_Contig")


vc_ffn_length.virus = vMAGs_ffn_length %>%
  select("qaccver","ORFLen", "sample_Contig", "sample", "contig", "genome", "ORF", "Source", "genome_Contig")
colnames(vc_ffn_length.virus) = c("qaccver","ORFLen", "sample_Contig", "sample", "contig", "genome", "ORF", "Source", "genome_Contig")
vc_ffn_length.virus[, Source := "vMAG"]

vc_ffn_length = rbindlist(list(vc_ffn_length.virus, vc_ffn_length.hosts), use.names = TRUE, fill = TRUE)


##########################################################################
###########################  contig depth  ##############################
##########################################################################
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract depth of contigs from:',in_contig_depth,'\n')
depth = in_contig_depth %>% {fread(., sep = '\t', header = T, data.table = T)}
colnames(depth)=c("sample_Contig","ContigLen","Avg_fold")
n_na_depth <- sum(is.na(depth$ContigLen)) + sum(is.na(depth$Avg_fold))
if (n_na_depth > 0) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING:',in_contig_depth,'contains',n_na_depth,'NA value(s) in ContigLen/Avg_fold; replaced with 0, please check the upstream data','\n')
  depth[is.na(ContigLen), ContigLen := 0]
  depth[is.na(Avg_fold), Avg_fold := 0]
}
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total contigs at the contig level:',nrow(depth), '\n')
# 使用 tstrsplit() 将 sample_Contig 列拆分成两列
depth[, c("sample", "contig") := tstrsplit(sample_Contig, "---", fill = "", keep = 1:2)]

depth = depth %>% select(-c(contig))

##############################################################
########### coverm trimmed mean depth based vMAGs   ###############
##############################################################
Contigs_depth = depth %>%
  merge(., Contigs_contig_genome, by = "sample_Contig", all.x = T) %>%
  {mutate(., ContigBases = ContigLen*Avg_fold)} %>%
  {mutate(., genome_Contig = ifelse(is.na(.$genome_Contig), .$sample_Contig, .$genome_Contig))} %>%
  {mutate(., genome = ifelse(is.na(.$genome), .$sample, .$genome))} %>%
  select(-c("sample", "genome"))
colnames(Contigs_depth) = c("sample_Contig","SeqLen","Avg_fold","genome_Contig","SeqBases")

vc_depth.virus = vMAGs_depth %>%
  select(c("sample_Contig","ContigLen","coverm_abundance","genome_Contig","SeqBases"))
colnames(vc_depth.virus) = c("sample_Contig","SeqLen","Avg_fold","genome_Contig","SeqBases")


vc_depth.hosts = Contigs_depth[!sample_Contig %in% vMAGs_depth.agg$sample_Contig]

vc_depth.overlap = Contigs_depth[sample_Contig %in% vMAGs_depth$sample_Contig] %>%
  merge(., vMAGs_depth.agg, by = "sample_Contig", all.x = T)

vc_depth.overlap[, Remove_Virus_SeqLen := SeqLen - Integrated_Virus_ContigLen]
vc_depth.overlap[, Remove_Virus_SeqBases := SeqBases - Integrated_Virus_SeqBases]
vc_depth.overlap[, Remove_Virus_Avg_fold := fifelse(Remove_Virus_SeqLen == 0,0,Remove_Virus_SeqBases / Remove_Virus_SeqLen)]

fully_viral_contigs <- vc_depth.overlap[is.na(Remove_Virus_SeqLen) | Remove_Virus_SeqLen <= 0, sample_Contig]
n_fully_viral_host_ORFs <- sum(vc_ffn_length.hosts$sample_Contig %in% fully_viral_contigs)
vc_ffn_length.hosts <- vc_ffn_length.hosts[!sample_Contig %in% fully_viral_contigs]
vc_ffn_length <- vc_ffn_length[!(Source == "Contig" & sample_Contig %in% fully_viral_contigs)]
cat(sprintf('[%s] Detected %d contigs whose entire depth was assigned to the virus branch (total viral region length >= contig length); excluded their %d raw-contig ORFs from the Contig branch (virus-based ORF predictions are used for these contigs).\n',
            format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            length(fully_viral_contigs), n_fully_viral_host_ORFs))
vc_depth.overlap <- vc_depth.overlap[Remove_Virus_SeqLen > 0]
vc_depth.overlap = vc_depth.overlap %>%
  select(c("sample_Contig","Remove_Virus_SeqLen","Remove_Virus_Avg_fold","genome_Contig","Remove_Virus_SeqBases"))
colnames(vc_depth.overlap) = c("sample_Contig","SeqLen","Avg_fold","genome_Contig","SeqBases")

vc_depth = rbindlist(list(vc_depth.overlap, vc_depth.virus, vc_depth.hosts), use.names = TRUE, fill = TRUE)
#length(unique(vc_depth$genome_Contig))

##########################################################################
###########   contig taxonomy  for virus  ##############整个模块移动了All_Taxa_abundance位置之后，其他模式在All_Taxa_abundance位置之前
##########################################################################

contig_seqid_tax[, genome_Contig := sample_Contig]

contig_taxonomy = contig_seqid_tax %>%
  select(-c("sample_Contig"))

v_bin_contig_genome = vc_ffn_length.virus %>%
  select(c("sample_Contig", "genome", "genome_Contig")) %>%
  unique()

if (CAMI_format_output) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Keep taxid and lineage column for contig-level annotation.\n')
  v_bin_seqid_tax = vMAGs_depth %>%
    merge(., virus_seqid_tax, all.x = T, by = "genome") %>%
    dplyr::select(., c("sample_Contig","taxid","Lineage","Taxonomy",
                       "taxid_Root","taxid_Domain", "taxid_Phylum", "taxid_Class", "taxid_Order", "taxid_Family", "taxid_Genus", "taxid_Species", "taxid_Strain",
                       "Lineage_Root","Lineage_Domain", "Lineage_Phylum", "Lineage_Class", "Lineage_Order", "Lineage_Family", "Lineage_Genus", "Lineage_Species", "Lineage_Strain",
                       "Domain","Phylum", "Class", "Order","Family", "Genus","Species", "Strain","NewTax","genome_Contig"))
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Do not keep taxid and lineage column for contig-level annotation.\n')
  v_bin_seqid_tax = vMAGs_depth %>%
    merge(., virus_seqid_tax, all.x = T, by = "genome") %>%
    select(c("sample_Contig","Taxonomy","NewTax", "Domain","Phylum", "Class", "Order","Family", "Genus","Species", "Strain","genome_Contig"))
}
# vc_depth.overlap，这些sample_Contig同时有两个分类，其他的宿主contig和病毒contig都只有一个分类
vc_seqid_tax.virus = v_bin_seqid_tax

vc_seqid_tax.hosts = contig_seqid_tax %>%
  {.[!.$sample_Contig %in% vMAGs_depth.agg$sample_Contig,]}

vc_seqid_tax.overlap = contig_seqid_tax %>%
  {.[.$sample_Contig %in% vc_depth.overlap$sample_Contig,]}

vc_seqid_tax = rbindlist(list(vc_seqid_tax.overlap, vc_seqid_tax.virus, vc_seqid_tax.hosts), use.names = TRUE, fill = TRUE)

vc_taxonomy = vc_seqid_tax[, setdiff(names(vc_seqid_tax), "sample_Contig"), with = FALSE]

contig_seqid_tax.inter = contig_seqid_tax %>%
  select(c("sample_Contig","Domain","Phylum", "Class", "Order", "Family","Genus","Species"))
v_bin_seqid_tax.inter = v_bin_seqid_tax %>%
  select(c("sample_Contig","genome_Contig","Domain","Phylum", "Class", "Order", "Family","Genus","Species"))


colnames(contig_seqid_tax.inter) = c("sample_Contig", "Domain_contig","Phylum_contig", "Class_contig", "Order_contig", "Family_contig","Genus_contig","Species_contig")
colnames(v_bin_seqid_tax.inter) = c("sample_Contig","genome_Contig","Domain_vMAG","Phylum_vMAG", "Class_vMAG", "Order_vMAG", "Family_vMAG","Genus_vMAG","Species_vMAG")

vc_consistency <- contig_seqid_tax.inter %>%
  {.[.$sample_Contig %in% v_bin_seqid_tax.inter$sample_Contig,]} %>%
  as.data.frame(.) %>%
  merge(., as.data.frame(v_bin_seqid_tax.inter), by = "sample_Contig", all.x = T) %>%
  mutate(
    Domain_consistency  = Domain_contig  == Domain_vMAG,
    Phylum_consistency  = Phylum_contig  == Phylum_vMAG,
    Class_consistency   = Class_contig   == Class_vMAG,
    Order_consistency   = Order_contig   == Order_vMAG,
    Family_consistency  = Family_contig  == Family_vMAG,
    Genus_consistency   = Genus_contig   == Genus_vMAG,
    Species_consistency = Species_contig == Species_vMAG
  ) %>%
  select(c("sample_Contig", "genome_Contig",
           "Domain_contig", "Domain_vMAG", "Domain_consistency",
           "Phylum_contig", "Phylum_vMAG", "Phylum_consistency",
           "Class_contig",  "Class_vMAG",  "Class_consistency",
           "Order_contig",  "Order_vMAG",  "Order_consistency",
           "Family_contig", "Family_vMAG", "Family_consistency",
           "Genus_contig",  "Genus_vMAG",  "Genus_consistency",
           "Species_contig","Species_vMAG","Species_consistency"))

vc_seqid_tax.inter.D =  vc_consistency %>%
  {.[.$sample_Contig %in% vMAGs_depth$sample_Contig,]} %>%
  {.[.$Domain_consistency %in% "FALSE",]}
cat(sprintf(
  '[%s] %d contigs contain partial viral signatures: they were classified as microbial/root/unclassified in metalink-c mode, but sections of the same sequence were flagged as viral in metalink-v mode. We annotate those viral segments with the corresponding viral taxonomy while keeping the original classification for the remaining regions.\n',
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  nrow(vc_seqid_tax.inter.D)
))

vc_seqid_tax.inter.S =  vc_consistency %>%
  {.[.$sample_Contig %in% vMAGs_depth$sample_Contig,]} %>%
  {.[.$Domain_consistency %in% "TRUE",]} %>%
  {.[.$Species_consistency %in% "FALSE",]}
cat(sprintf(
  '[%s] %d viral sequences received discordant species-level assignments between metalink-c mode and metalink-v mode; the metalink-v mode classification is adopted.\n',
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  nrow(vc_seqid_tax.inter.S)
))

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output virus and MAGs/contigs consistency to:',out_vc_consistency,'\n')
vc_consistency %>%
  {fwrite(., out_vc_consistency, row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}



################################################################
####################     read length      ######################
################################################################
read_length = in_read_length %>% {fread(., sep = '\t', header = T, data.table = T)}
ReadLength = max(read_length$avg_len, na.rm = TRUE)


vc_depth_sum = sum(vc_depth$SeqBases)

# left join of vc_taxonomy (x) with per-contig depth (y) on genome_Contig,
# equivalent to merge(vc_taxonomy, vc_depth, all.x = TRUE, by = "genome_Contig")
vc_seqid_tax_depth <- vc_depth[, .(genome_Contig, SeqLen, Avg_fold, SeqBases, sample = sample_name)][
  vc_taxonomy, on = "genome_Contig"]
vc_seqid_tax_depth[, Sum_Bases_vc := vc_depth_sum]
vc_seqid_tax_depth[, ContigRA_Bases := 1e2 * SeqBases/Sum_Bases_vc]
colnames(vc_seqid_tax_depth)[1] = "SeqID"

Domain_depth  = make_vc_depth(vc_seqid_tax_depth, "Domain",  "D", ReadLength)
Phylum_depth  = make_vc_depth(vc_seqid_tax_depth, "Phylum",  "P", ReadLength)
Class_depth   = make_vc_depth(vc_seqid_tax_depth, "Class",   "C", ReadLength)
Order_depth   = make_vc_depth(vc_seqid_tax_depth, "Order",   "O", ReadLength)
Family_depth  = make_vc_depth(vc_seqid_tax_depth, "Family",  "F", ReadLength)
Genus_depth   = make_vc_depth(vc_seqid_tax_depth, "Genus",   "G", ReadLength)
Species_depth = make_vc_depth(vc_seqid_tax_depth, "Species", "S", ReadLength)
Strain_depth  = make_vc_depth(vc_seqid_tax_depth, "Strain",  "T", ReadLength)

All_Taxa_abundance = rbindlist(list(
  Domain_depth,  Phylum_depth,
  Class_depth,   Order_depth,
  Family_depth,  Genus_depth,
  Species_depth, Strain_depth))

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the integrated abundance table of vMAGs, mMAGs and contigs to:',out_Taxa_abundance,'\n')
All_Taxa_abundance[, .(sample, Level, Taxonomy,
                   Taxa_ContigNum, Taxa_ReadCount, Taxa_RelativeAbundance)] %>%
  {fwrite(., paste0(out_Taxa_abundance), row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

#############################################
############### out_CAMI_format ##############
#############################################

if (CAMI_format_output) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the integrated abundance table of contigs with CAMI format to:',out_CAMI_format,'\n')

  vc_Lineage_Strain_depth = make_vc_lineage_depth(vc_seqid_tax_depth, c("Strain","Lineage_Strain"), "strain", ReadLength, 8)
  vc_Lineage_Species_depth = make_vc_lineage_depth(vc_seqid_tax_depth, c("Species","Lineage_Species"), "species", ReadLength, 7)
  vc_Lineage_Genus_depth = make_vc_lineage_depth(vc_seqid_tax_depth, c("Genus","Lineage_Genus"), "genus", ReadLength, 6)
  vc_Lineage_Family_depth = make_vc_lineage_depth(vc_seqid_tax_depth, c("Family","Lineage_Family"), "family", ReadLength, 5)
  vc_Lineage_Order_depth = make_vc_lineage_depth(vc_seqid_tax_depth, c("Order","Lineage_Order"), "order", ReadLength, 4)
  vc_Lineage_Class_depth = make_vc_lineage_depth(vc_seqid_tax_depth, c("Class","Lineage_Class"), "class", ReadLength, 3)
  vc_Lineage_Phylum_depth = make_vc_lineage_depth(vc_seqid_tax_depth, c("Phylum","Lineage_Phylum"), "phylum", ReadLength, 2)
  vc_Lineage_Domain_depth = make_vc_lineage_depth(vc_seqid_tax_depth, c("Domain","Lineage_Domain"), "superkingdom", ReadLength, 1)
  vc_Lineage_Root_depth = make_vc_lineage_depth(vc_seqid_tax_depth, c("taxid_Root","Lineage_Root"), "root", ReadLength)
  vc_Lineage_Root_depth[, Taxonomy := gsub("1", "root", gsub("0", "unclassified", Taxonomy))]

  All_Taxa_Lineage_depth = rbindlist(list(
    vc_Lineage_Root_depth,
    vc_Lineage_Domain_depth,  vc_Lineage_Phylum_depth,
    vc_Lineage_Class_depth,   vc_Lineage_Order_depth,
    vc_Lineage_Family_depth,  vc_Lineage_Genus_depth,
    vc_Lineage_Species_depth, vc_Lineage_Strain_depth))
  All_Taxa_Lineage_depth[, taxid := gsub(".*;", "", Lineage)]
  All_Taxa_Lineage_depth[, Lineage := gsub(";", "|", Lineage)]
  All_Taxa_Lineage_depth[, Lineage := sub("^1\\|", "", Lineage)]
  All_Taxa_Lineage_depth[, Taxonomy := gsub(";", "|", Taxonomy)]
  All_Taxa_Lineage_depth[, Taxonomy := sub("^d__", "", Taxonomy)]
  All_Taxa_Lineage_depth[, Taxonomy := sub("\\|p__", "\\|", Taxonomy)]
  All_Taxa_Lineage_depth[, Taxonomy := sub("\\|c__", "\\|", Taxonomy)]
  All_Taxa_Lineage_depth[, Taxonomy := sub("\\|o__", "\\|", Taxonomy)]
  All_Taxa_Lineage_depth[, Taxonomy := sub("\\|f__", "\\|", Taxonomy)]
  All_Taxa_Lineage_depth[, Taxonomy := sub("\\|g__", "\\|", Taxonomy)]
  All_Taxa_Lineage_depth[, Taxonomy := sub("\\|s__", "\\|", Taxonomy)]
  All_Taxa_Lineage_depth[, Taxonomy := sub("\\|t__", "\\|", Taxonomy)]
  #All_Taxa_Lineage_depth[, Taxonomy := mapply(sync_empty_levels, Lineage, Taxonomy)]
  #All_Taxa_Lineage_depth[, Taxonomy := gsub("^\\|+", "", Taxonomy)]
  #All_Taxa_Lineage_depth[, Taxonomy := gsub("\\|+$", "", Taxonomy)]
  All_Taxa_Lineage_depth = All_Taxa_Lineage_depth %>%
    as.data.frame(.)%>%
    select(c("taxid","Level","Lineage","Taxonomy", "Taxa_RelativeAbundance"))
  colnames(All_Taxa_Lineage_depth) = c("@@TAXID","RANK","TAXPATH","TAXPATHSN","PERCENTAGE")

  header_lines <- c(
    "# Taxonomic Profiling Output",
    paste0("@SampleID:",sample_name),
    paste0("@Version:", Sys.getenv("METALINK_VERSION", "unknown")),
    "@Ranks:superkingdom|phylum|class|order|family|genus|species|strain",
    paste0("@TaxonomyID:",in_taxid_tax),
    "@__program__:metalink vc",
    "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE"
  )
  #out_CAMI_format <- "output_cami_profile.txt"
  writeLines(header_lines, out_CAMI_format)
  write.table(All_Taxa_Lineage_depth, out_CAMI_format,sep = "\t",row.names = F,col.names = F,append = T,quote = F)
}
#############################################
######## ORF length for vc ##########
#############################################
vc_ffn_length_depth <- vc_depth[, .(genome_Contig, sample_Contig, Avg_fold)][
  vc_ffn_length[, .(genome_Contig, qaccver, ORFLen, sample, Source)], on = "genome_Contig"]
vc_ffn_length_depth <- vc_ffn_length_depth[, .(sample_Contig, qaccver, ORFLen, Avg_fold,
                                               genome_Contig, sample, Source)]
# vc_taxonomy already carries genome_Contig (from the rbindlist of the three
# vc_seqid_tax parts); no column-rename safeguard needed here


##########################################################################
###########   ORF function of contigs   ##############
##########################################################################
if (contigs_function_classify_method == "emapper"){
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract ORF function of contigs from:',in_contig_fun_detail,'\n')
  contigs_fun_detail <- fread(in_contig_fun_detail, sep = "\t", header = T, skip = 1, data.table = T)
  contigs_fun_detail <- contigs_fun_detail[, .SD, .SDcols = c(1, 12)]
  colnames(contigs_fun_detail) = c('qaccver','KO')
  contigs_fun_detail[, KO := gsub("ko:", "", KO)]
  contigs_fun_detail[, KO := gsub(",", "_", KO)]
} else if (contigs_function_classify_method == "kraken2"){
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract ORF function of contigs from:',in_contig_fun_detail,'\n')
  contigs_fun_detail = fread(in_contig_fun_detail, sep = "\t", select = c(2,3), header = FALSE, data.table = T)
  colnames(contigs_fun_detail) = c("qaccver", "taxid")
  contigs_fun_detail$taxid = contigs_fun_detail$taxid %>%
    gsub(".* \\(taxid |\\)", "", .)
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract mapping relationship of taxid-function from:',in_taxid_fun,'\n')
  taxid_function = fread(in_taxid_fun, sep = "\t", colClasses = "character", col.names = c('taxid','Fun_Taxonomy','Lineage'), header = FALSE, data.table = T)
  taxid_function[, Strain := tstrsplit(Fun_Taxonomy, ";", fixed = TRUE)[[8]]]
  taxid_function[, Strain := gsub("d__unclassified strain", "Unclassified", Strain)]
  taxid_function[, Strain := gsub("root strain", "Root", Strain)]
  # Vectorized equivalent of the original by-row loop:
  #   vec = strsplit(Strain, " "); if (length(vec) >= 3) sub("^[a-z]__", "", vec[2])
  #   else sub("^[a-z]__", "", Strain)
  # (strsplit drops trailing empty tokens while tstrsplit keeps them, hence the
  #  trailing-space guard in has3.)
  tok <- tstrsplit(taxid_function$Strain, " ", fixed = TRUE)
  tok3 <- if (length(tok) >= 3) tok[[3]] else rep(NA_character_, nrow(taxid_function))
  has3 <- !is.na(tok3) & !(tok3 == "" & grepl(" $", taxid_function$Strain))
  taxid_function[, KO := ifelse(has3, sub("^[a-z]__", "", tok[[2]]), sub("^[a-z]__", "", Strain))]
  # only taxid and KO are used downstream (Fun_Taxonomy/Lineage/Strain never
  # reach any output)
  taxid_function = taxid_function[, .(taxid, KO)]
  rm(tok, tok3, has3)
  contigs_fun_detail = contigs_fun_detail %>%
    merge.data.table(., taxid_function, by = 'taxid', all.x = T)
  rm(taxid_function)
}

contigs_fun_detail[, c("sample_Contig") := tstrsplit(qaccver, "_\\d*$")]
contigs_fun_detail = contigs_fun_detail[!KO %in% '-']

# Fun_Taxonomy (kraken2 branch) never reaches any output; keep qaccver and KO only
contigs_fun_detail = contigs_fun_detail[, .(qaccver, KO)]

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total annotated ORFs at the contig level:',nrow(contigs_fun_detail), '\n')
##########################################################################
###########   ORF function of virus   ##############
##########################################################################
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract ORF function of viral contigs from:',in_vMAGs_ffn_length,'\n')
virus_fun_detail = vMAGs_ffn_infor %>%
  as.data.frame(.) %>%
  select(c("qaccver","KO")) %>%
  filter(!KO %in% '')

##########################################################
###########  merge ORF function of v and c   ############
##########################################################
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Combine ORF function of virus and contigs...','\n')
vc_fun_detail.hosts = contigs_fun_detail %>%
  {.[.$qaccver %in% unique(vc_ffn_length.hosts$qaccver), ]}

vc_fun_detail.virus = virus_fun_detail

vc_fun_detail = rbindlist(list(vc_fun_detail.virus, vc_fun_detail.hosts), use.names = TRUE, fill = TRUE)
#############################################
######## function abundance for c ##########
#############################################

seperate_emapper_ko_depth <- function(data) {
  setDT(data)
  data[, .(KO = unlist(strsplit(KO, "_"))), by = .(ORF, Avg_fold)
  ][, Duplicate_Annotations_num := .N, by = ORF]
}

if (KO_mode == "split_multiple_KOs") {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'For MAGs and contigs: Split the depth of an ORF with multiple KO annotations by evenly distributing the total value across all assigned KOs.','\n')
  vc_ORF_KO_depth = vc_fun_detail[vc_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  vc_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  vc_ORF_KO_depth = seperate_emapper_ko_depth(vc_ORF_KO_depth[, .(ORF = qaccver, Avg_fold, KO)])
  setnames(vc_ORF_KO_depth, "ORF", "qaccver")
} else if (KO_mode == "keep_multiple_KOs") {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'For MAGs and contigs: Keep original multiple KO annotations.','\n')
  vc_ORF_KO_depth = vc_fun_detail[vc_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  vc_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  vc_ORF_KO_depth[, Duplicate_Annotations_num := 1]
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'For MAGs and contigs: Ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.','\n')
  vc_ORF_KO_depth = vc_fun_detail[vc_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  vc_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  vc_ORF_KO_depth[, Duplicate_Annotations_num := 1]
  vc_ORF_KO_depth = vc_ORF_KO_depth[!str_detect(KO, "_")]
}
vc_ORF_KO_depth[, Scaled_Avg_fold := Avg_fold / Duplicate_Annotations_num]
vc_ORF_KO_depth[vc_ffn_length_depth, on = "qaccver",
                `:=`(sample_Contig = i.sample_Contig,
                     ORFLen = i.ORFLen,
                     sample = i.sample,
                     genome_Contig = i.genome_Contig,
                     Source = i.Source)]
setorderv(vc_ORF_KO_depth, c("sample_Contig", "qaccver"), na.last = TRUE)

setDT(vc_ORF_KO_depth)
need_tax_cols <- unique(c("Taxonomy", "NewTax", "Domain", rank_vec))
if (CAMI_format_output) need_tax_cols <- c(need_tax_cols, "taxid", "Lineage")
vc_ORF_KO_depth[vc_taxonomy, on = "genome_Contig",
                (need_tax_cols) := mget(paste0("i.", need_tax_cols))]

#############################################
######## function abundance for vc #########
#############################################
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Combine ORFs KO annotations of vMAGs and contigs...','\n')

vc_ORF_KO_depth[, scaleOrfBases := ORFLen*Scaled_Avg_fold]
vc_ORF_KO_depth <- vc_ORF_KO_depth[ORFLen >= as.numeric(minORFLen)]

ORF_depth_sum = vc_ORF_KO_depth[, .(sumOrfBases = sum(scaleOrfBases, na.rm = TRUE),
                                    sumOrfDepth = sum(Scaled_Avg_fold, na.rm = TRUE)), keyby = sample]

ORF_eggnog_depth = vc_ORF_KO_depth
ORF_eggnog_depth[ORF_depth_sum, on = "sample",
                 `:=`(sumOrfBases = i.sumOrfBases, sumOrfDepth = i.sumOrfDepth)]
ORF_eggnog_depth[, ORFRA := 1e6 * Scaled_Avg_fold/sumOrfDepth]
ORF_eggnog_depth[, ORFRA_Bases := 1e6 * scaleOrfBases/sumOrfBases]
ORF_eggnog_depth[, NewTax := gsub("unclassified", "Unclassified", NewTax)]
ORF_eggnog_depth[, NewTax := gsub("New", "Unclassified", NewTax)]

#############################################
############### summary for KO ##############
#############################################

KO_eggnog_depth = ORF_eggnog_depth[, .(KO_ORFNum = .N,
                                       KO_SumCoverage = sum(Scaled_Avg_fold, na.rm = TRUE),
                                       KO_RelativeAbundance = sum(ORFRA, na.rm = TRUE),
                                       ORFLen_median = floor(median(ORFLen, na.rm = TRUE))), keyby = .(sample, KO)]

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the table of KOs mean coverage and relative abundance to:',out_KO_eggnog_depth,'\n')
KO_eggnog_depth %>%
  {fwrite(., paste0(out_KO_eggnog_depth), row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

#############################################
############ summary for KO_Taxa ############
#############################################

KO_Taxa_SumCoverage = ORF_eggnog_depth[, .(Taxa_KO_ORF_Num = .N,
                                     Taxa_KO_SumCoverage = sum(Scaled_Avg_fold, na.rm = TRUE),
                                     Taxa_KO_RelativeAbundance = sum(ORFRA, na.rm = TRUE)),
                                 keyby = .(sample, KO, Taxonomy)]
KO_Taxa_SumCoverage[, sTaxonomy := paste0(sample,"---",Taxonomy)]

KO_Taxa_SumCoverage <- KO_Taxa_SumCoverage[, if (.N >= as.numeric(min_KO_num_in_one_taxon)) .SD, by = sTaxonomy]
if (nrow(KO_Taxa_SumCoverage) > 0) {
  KO_Taxa_SumCoverage_wide <- data.table::dcast(
    KO_Taxa_SumCoverage,
    KO ~ sTaxonomy,
    value.var = "Taxa_KO_SumCoverage",
    fun.aggregate = sum,
    fill = 0
  )
  KO_Taxa_RelativeAbundance_wide <- data.table::dcast(
    KO_Taxa_SumCoverage,
    KO ~ sTaxonomy,
    value.var = "Taxa_KO_RelativeAbundance",
    fun.aggregate = sum,
    fill = 0
  )
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output mean coverage table of KO-Taxa to:',out_KO_Taxa_SumCoverage_wide,'\n')
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output relative abundance table of KO-Taxa to:',out_KO_Taxa_RelativeAbundance_wide,'\n')
  KO_Taxa_SumCoverage_wide %>%
    {fwrite(., out_KO_Taxa_SumCoverage_wide, row.names = FALSE, col.names = T,
            sep = '\t', quote = FALSE)}
  KO_Taxa_RelativeAbundance_wide %>%
    {fwrite(., out_KO_Taxa_RelativeAbundance_wide, row.names = FALSE, col.names = T,
            sep = '\t', quote = FALSE)}
} else {
  warning("No taxon were retained, maybe min_KO_num_in_one_taxon is too large.")
  fwrite(KO_Taxa_SumCoverage, out_KO_Taxa_SumCoverage_wide, row.names = FALSE, col.names = T,sep = '\t', quote = FALSE)
  fwrite(KO_Taxa_SumCoverage, out_KO_Taxa_RelativeAbundance_wide, row.names = FALSE, col.names = T,sep = '\t', quote = FALSE)
}

#' Summarise ORF_eggnog_depth by a chosen taxonomic rank and write two wide tables to disk
#'
#' @param rank Character string specifying the taxonomic rank to group by,
#'   e.g. "Phylum", "Class", "Order", etc.
#' @param ORF_eggnog_depth Input data frame that must contain the columns:
#'   sample, KO, Scaled_Avg_fold, and ORFRA.
#' @param out_KO_Taxa_SumCoverage_wide File path (without extension) for the
#'   depth-wide output TSV.
#' @param out_KO_Taxa_RelativeAbundance_wide File path (without extension) for the
#'   relative-abundance-wide output TSV.
#'
#' @return Invisible NULL. Two TSV files are written to disk:
#'   one for depth and one for relative abundance.
#' @export
make_KO_Taxa_tables <- function(rank,
                                ORF_eggnog_depth,
                                out_KO_Taxa_SumCoverage_wide,
                                out_KO_Taxa_RelativeAbundance_wide,
                                min_KO_num_in_one_taxon = 100) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Input rank for make_KO_Taxa_tables function:',rank,'\n')
  if (!rank %in% colnames(ORF_eggnog_depth))
    stop(sprintf("input dataframe lack col: %s", rank))
  KO_Taxa_SumCoverage_out <- ORF_eggnog_depth[, .(ORF_Num = .N,
                                            Taxa_KO_SumCoverage = sum(Scaled_Avg_fold, na.rm = TRUE),
                                            Taxa_KO_RelativeAbundance = sum(ORFRA, na.rm = TRUE)),
                                        keyby = c("sample", "KO", rank)]
  out_KO_Taxa_long = sprintf("%s_%s.tsv", sub("_seq$", "", sub("\\.tsv$", "", detail)), rank)
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the table of aggregated taxonomy, function, and abundance to:',out_KO_Taxa_long,'\n')
  fwrite(KO_Taxa_SumCoverage_out, out_KO_Taxa_long, row.names = FALSE, col.names = T,sep = '\t', quote = FALSE)

  KO_Taxa_SumCoverage = KO_Taxa_SumCoverage_out[, srank := paste(sample, get(rank), sep = "---")]

  KO_Taxa_SumCoverage <- KO_Taxa_SumCoverage[, if (.N >= as.numeric(min_KO_num_in_one_taxon)) .SD, by = srank]
  if (nrow(KO_Taxa_SumCoverage) > 0) {
    KO_Taxa_SumCoverage_wide <- data.table::dcast(KO_Taxa_SumCoverage, KO ~ srank, value.var = "Taxa_KO_SumCoverage", fun.aggregate = sum, fill = 0)
    KO_Taxa_RelativeAbundance_wide <- data.table::dcast(KO_Taxa_SumCoverage, KO ~ srank, value.var = "Taxa_KO_RelativeAbundance", fun.aggregate = sum, fill = 0)
    out_KO_Taxa_SumCoverage_wide = sprintf("%s.%s.tsv", sub("\\.tsv$", "", out_KO_Taxa_SumCoverage_wide), rank)
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output mean coverage table of KO-Taxa to:',out_KO_Taxa_SumCoverage_wide,'\n')
    fwrite(KO_Taxa_SumCoverage_wide, out_KO_Taxa_SumCoverage_wide, row.names = FALSE, col.names = T,sep = '\t', quote = FALSE)
    out_KO_Taxa_RelativeAbundance_wide = sprintf("%s.%s.tsv", sub("\\.tsv$", "", out_KO_Taxa_RelativeAbundance_wide), rank)
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output relative abundance table of KO-Taxa to:',out_KO_Taxa_RelativeAbundance_wide,'\n')
    fwrite(KO_Taxa_RelativeAbundance_wide, out_KO_Taxa_RelativeAbundance_wide, row.names = FALSE, col.names = T,sep = '\t', quote = FALSE)
  } else {
    warning("No taxon were retained, maybe min_KO_num_in_one_taxon is too large.")
    fwrite(KO_Taxa_SumCoverage, sprintf("%s.%s.tsv", sub("\\.tsv$", "", out_KO_Taxa_SumCoverage_wide), rank), row.names = FALSE, col.names = T,sep = '\t', quote = FALSE)
    fwrite(KO_Taxa_SumCoverage, sprintf("%s.%s.tsv", sub("\\.tsv$", "", out_KO_Taxa_RelativeAbundance_wide), rank), row.names = FALSE, col.names = T,sep = '\t', quote = FALSE)
  }
  return(KO_Taxa_SumCoverage_out)
}
for (rk in rank_vec) {
  make_KO_Taxa_tables(rk, ORF_eggnog_depth, out_KO_Taxa_SumCoverage_wide, out_KO_Taxa_RelativeAbundance_wide, as.numeric(min_KO_num_in_one_taxon))
}

#############################################
########### summary for ORF_Tax_KO ##########
#############################################

# summary
ORF_eggnog_depth[, KO_annotate := ifelse(KO %in% c("Unclassified"),"Unclassified","Classified")]
ORF_Tax_KO_summary <- summary_with_all(ORF_eggnog_depth,
                                       c("Domain", "NewTax", "KO_annotate"),
                                       "ORF_Num", "Total_ORF", c("Domain", "NewTax"))

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the summary table of ORFs taxonomy and function to:',summary_ORF_Tax_KO,'\n')
#summary = "{sample}.contig_Tax_summary.tsv
ORF_Tax_KO_summary %>%
  {fwrite(., summary_ORF_Tax_KO, row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

#############################################
############# summary for ORF KO ############
#############################################

ORF_KO_summary <- summary_with_all(ORF_eggnog_depth,
                                   c("Domain", "KO_annotate"),
                                   "ORF_Num", "Total_ORF", "Domain")

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the summary table of ORFs KO annotations to:',summary_ORF_KO,'\n')
#summary = "{sample}.contig_Tax_summary.tsv"
ORF_KO_summary %>%
  {fwrite(., summary_ORF_KO, row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

#############################################
########### summary for ORF_Tax_Source ##########
#############################################
ORF_Tax_Source_summary <- summary_with_all(ORF_eggnog_depth,
                                           c("Domain", "NewTax", "Source"),
                                           "ORF_Num", "Total_ORF", c("Domain", "NewTax"))

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the summary table of ORFs taxonomy and source to:',summary_ORF_Tax_Source,'\n')
#summary = "{sample}.contig_Tax_summary.tsv
ORF_Tax_Source_summary %>%
  {fwrite(., summary_ORF_Tax_Source, row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

#############################################
############# summary for ORF Source ############
#############################################
ORF_Source_summary <- summary_with_all(ORF_eggnog_depth,
                                       c("Domain", "Source"),
                                       "ORF_Num", "Total_ORF", "Domain")

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the summary table of ORFs source to:',summary_ORF_Source,'\n')
#summary = "{sample}.contig_Tax_summary.tsv"
ORF_Source_summary %>%
  {fwrite(., summary_ORF_Source, row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

#############################################
########### detail output (write last to ########
########### avoid the subset copy)     ########
#############################################
if (CAMI_format_output) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Keep taxid and lineage column for',detail,'\n')
  detail_cols = c("qaccver", "KO",  "Taxonomy","NewTax","taxid","Lineage",
                  "Scaled_Avg_fold","ORFLen", "ORFRA","ORFRA_Bases","Source")
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Do not keep taxid and lineage column for',detail,'\n')
  detail_cols = c("qaccver", "KO",  "Taxonomy","NewTax",
                  "Scaled_Avg_fold","ORFLen", "ORFRA","ORFRA_Bases","Source")
}
#sum(ORF_eggnog_depth$scaleOrfBases)
if (output_detail) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the table of ORFs taxonomy, function, and abundance to:',detail,'\n')
  # write the whole table after dropping non-output columns: avoids the
  # N x length(detail_cols) pointer copy made by dt[, ..detail_cols]
  ORF_eggnog_depth[, setdiff(colnames(ORF_eggnog_depth), detail_cols) := NULL]
  setcolorder(ORF_eggnog_depth, detail_cols)
  fwrite(ORF_eggnog_depth, paste0(detail), row.names = FALSE, col.names = T,
         sep = '\t', quote = FALSE)
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Skip output of detail file (--output_detail is FALSE)\n')
}

rm(ORF_eggnog_depth)
#gc()

#############################################
############## summary for Taxa #############
#############################################
c_taxonomy = copy(contig_seqid_tax)
c_taxonomy[, NewTax := gsub("unclassified", "Unclassified", NewTax)]
c_taxonomy[, NewTax := gsub("New", "Unclassified", NewTax)]

c_Tax_summary <- summary_with_all(c_taxonomy,
                                  c("Domain", "NewTax"),
                                  "Contig_Num_c", "Total_Contig_c", "Domain",
                                  "Proportion_Contig_c")

v_taxonomy = as.data.table(copy(virus_seqid_tax))
v_taxonomy[, NewTax := gsub("unclassified", "Unclassified", NewTax)]
v_taxonomy[, NewTax := gsub("New", "Unclassified", NewTax)]

v_Tax_summary <- summary_with_all(v_taxonomy,
                                  c("Domain", "NewTax"),
                                  "vMAG_Num", "Total_vMAG", "Domain",
                                  "Proportion_v")

vc_taxonomy[, NewTax := gsub("unclassified", "Unclassified", NewTax)]
vc_taxonomy[, NewTax := gsub("New", "Unclassified", NewTax)]

vc_Tax_summary <- summary_with_all(vc_taxonomy,
                                   c("Domain", "NewTax"),
                                   "Contig_Num_vc", "Total_Contig_vc", "Domain",
                                   "Proportion_Contig_vc")

# base merge (sort=TRUE) reproduces the original locale-sorted row order;
# inputs are converted to data.frame so merge.data.frame (not merge.data.table)
# is dispatched
all_Tax_summary = as.data.frame(vc_Tax_summary) %>%
  merge(., as.data.frame(c_Tax_summary), all = T, by = c("Domain","NewTax")) %>%
  merge(., as.data.frame(v_Tax_summary), all = T, by = c("Domain","NewTax"))
# mutate(across(where(is.numeric), ~replace_na(., 0))), preserving column types
setDT(all_Tax_summary)
for (cl in names(all_Tax_summary)) {
  if (is.numeric(all_Tax_summary[[cl]]))
    set(all_Tax_summary, which(is.na(all_Tax_summary[[cl]])), cl,
        if (is.integer(all_Tax_summary[[cl]])) 0L else 0)
}


cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the summary table of vMAGs and contigs taxonomy to:',summary_Tax,'\n')
#summary = "{sample}.cr_Tax_summary.tsv"
all_Tax_summary %>%
  {fwrite(., summary_Tax, row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

mode <- "vc"
if (!is.null(out_mode_flag)) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Exporting execution mode:',mode,'to',out_mode_flag,'\n')
  dir.create(dirname(out_mode_flag), showWarnings = FALSE, recursive = TRUE)
  writeLines(mode, out_mode_flag)
}
