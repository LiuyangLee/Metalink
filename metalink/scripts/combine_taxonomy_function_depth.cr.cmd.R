###
#' @Date: 2023-06-14 16:48:43
#' @LastEditors: LiuyangLi
#' @LastEditTime: 2026-09-09
#' @FilePath: /liliuyang/R/Rscript.genome/combine_taxonomy_function_depth.cr.cmd.R
#' @Description:
#'   Since Metalink-0.3.0: read-side ORF depth is normalized by KO reference
#'   length (--in_KO_length_ref): each read-derived ORF contributes
#'   ORFLen / (3 * len_median) instead of 1, so that read-side KO depths share
#'   the coverage scale of the contig side. ORFs whose KO is absent from the
#'   reference table fall back to the global median of len_median.
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

#sample_name = "test"
#KO_mode='ignore_multiple_KOs'
#contigs_taxonomy_classify_method = "CAT_pack"
#reads_function_classify_method = "emapper"
#read_type='paired'
#read_type='single'
#source("C:\\Users\\Administrator\\Desktop\\Script\\Rscript\\R-20241204\\RLib\\R\\RLib.r")
#min_KO_num_in_one_taxon = 100
#minORFLen=50
#in_read_paired_tax_detail = '../data/02_kraken2_unmapped_reads_paired/Y36_17.kraken2.read2classification.tsv.gz'
#in_read_singleton_tax_detail = '../data/02_kraken2_unmapped_reads_singleton/Y36_17.kraken2.read2classification.tsv.gz'
#in_read_singleton_tax_detail = '../data/02_kraken2_unmapped_reads_singleton/SRR10325488_top10.kraken2.read2classification.tsv.gz'
#in_read_length = "../data/Y36_17_trimmed_1.fastq.stats"
#in_read_ffn_length = "../data/reads_ffn_length/Y36_17_unmapped_both_merge_unmerge.ffn_length.gz"
#in_read_fun_detail = '../data/03_eggnog_function/Y36_17.reads_cr.emapper.annotations.gz'
#in_taxid_tax="../taxid/taxid_tax_lineage.tsv.gz"
#in_taxid_fun="../taxid/taxid_fun_lineage.tsv.gz"
#in_contig_fun_detail = "../data/Y36_17.emapper.annotations.gz"
#in_contig_ffn_length = "../data/contigs_ffn_length/Y36_17"
#in_contig_depth = "../data/contigs_metabat_depth/Y36_17"
#in_contig_tax_detail="../data/Y36_17.CAT_pack.contig2classification.tsv"
#detail = "all.test.detail.tsv"
#out_Taxa_abundance = "test.detail.2.tsv"
script_path <- sub("--file=", "", commandArgs()[grepl("--file=", commandArgs())])
if (length(script_path) == 0 || script_path == "") {
  script_path <- sys.frame(1)$ofile
}
source(file.path(dirname(script_path), "RLib.r"))
# GForce-optimized gsum rounds differently from base R sum() (which dplyr uses
# per group); fall back to base evaluation so grouped sums are bit-identical to
# the dplyr implementation
options(datatable.optimize = 0L)
in_read_paired_tax_detail = NULL
in_read_singleton_tax_detail = NULL
in_read_ffn_length = NULL
in_read_fna_length = NULL
in_taxid_tax=NULL
in_taxid_fun=NULL
in_contig_fun_detail = NULL
in_read_fun_detail = NULL
in_contig_ffn_length = NULL
in_contig_depth = NULL
in_contig_tax_detail=NULL
detail = NULL
out_Taxa_abundance = NULL

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
  make_option("--reads_function_classify_method", type="character", default=NULL,
              help="kraken2 or emapper"),
  make_option("--in_read_paired_tax_detail", type="character", default=NULL,
              help="in_read_paired_tax_detail.tsv"),
  make_option("--in_read_singleton_tax_detail", type="character", default=NULL,
              help="in_read_singleton_tax_detail.tsv"),
  make_option("--in_read_fun_detail", type="character", default=NULL,
              help="Y36_17.reads_protein.emapper.annotations.gz"),
  make_option("--in_contig_ffn_length", type="character", default=NULL,
              help="in_contig_ffn_length.tsv"),
  make_option("--in_read_ffn_length", type="character", default=NULL,
              help="in_read_ffn_length.tsv"),
  make_option("--in_read_fna_length", type="character", default=NULL,
              help="in_read_fna_length.tsv (per-read length for long reads)"),
  make_option("--in_contig_depth", type="character", default=NULL,
              help="in_contig_depth.tsv"),
  make_option("--in_read_length", type="character", default=NULL,
              help="sample_trimmed_1.fastq.stats"),
  make_option("--sample_name", type="character", default=NULL,
              help="sample_name, e.g., 'SRR000001' (default: test)"),
  make_option("--read_type", type="character", default=NULL,
              help="read type: 'single' or 'paired' (default: paired)"),
  make_option("--KO_mode", type="character", default="ignore_multiple_KOs",
              help="KO handling mode: 'keep_multiple_KOs', 'ignore_multiple_KOs', or 'split_multiple_KOs' (default: ignore_multiple_KOs)"),
  make_option("--in_KO_length_ref", type="character", default=NULL,
              help="KO reference length table with accession and len_median columns (amino acids), e.g., KO_faa_seq_length_summary.tsv.gz; used to normalize read-side ORF depth"),
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
reads_function_classify_method <- opt$reads_function_classify_method
in_read_paired_tax_detail <- opt$in_read_paired_tax_detail
in_read_singleton_tax_detail <- opt$in_read_singleton_tax_detail
in_read_fun_detail <- opt$in_read_fun_detail
in_contig_ffn_length <- opt$in_contig_ffn_length
in_read_ffn_length <- opt$in_read_ffn_length
in_read_fna_length <- opt$in_read_fna_length
in_contig_depth <- opt$in_contig_depth
in_read_length <- opt$in_read_length
sample_name <- opt$sample_name
read_type <- opt$read_type
KO_mode <- opt$KO_mode
in_KO_length_ref <- opt$in_KO_length_ref
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

if (is.null(read_type) || !read_type %in% c('paired','single',
                                            'paired-rna','single-rna',
                                            'pacbio-raw', 'pacbio-corr', 'pacbio-hifi',
                                            'nano-raw', 'nano-corr', 'nano-hq')) {
  stop("--read_type must be one of 'paired','single','paired-rna','single-rna', 'pacbio-raw', 'pacbio-corr', 'pacbio-hifi','nano-raw', 'nano-corr', 'nano-hq'")
}

if (is.null(KO_mode) || !KO_mode %in% c("keep_multiple_KOs","ignore_multiple_KOs","split_multiple_KOs")) {
  stop("--KO_mode must be 'keep_multiple_KOs', 'ignore_multiple_KOs', or 'split_multiple_KOs'")
}

if (is.null(in_KO_length_ref) || !file.exists(in_KO_length_ref)) {
  stop("--in_KO_length_ref must be specified and point to an existing KO reference length table (accession + len_median columns)")
}

if (is.null(contigs_taxonomy_classify_method) || !contigs_taxonomy_classify_method %in% c("kraken2","CAT_pack")) {
  stop("--contigs_taxonomy_classify_method must be 'kraken2' or 'CAT_pack'")
}

if (is.null(contigs_function_classify_method) || !contigs_function_classify_method %in% c("kraken2","emapper")) {
  stop("--contigs_function_classify_method must be 'kraken2' or 'emapper'")
}

if (is.null(reads_function_classify_method) || !reads_function_classify_method %in% c("kraken2","emapper")) {
  stop("--reads_function_classify_method must be 'kraken2' or 'emapper'")
}

if (read_type %in% c('paired','paired-rna')) {
  if (is.null(in_read_paired_tax_detail) || is.null(in_read_singleton_tax_detail)) {
    stop("In paired mode, both --in_read_paired_tax_detail and --in_read_singleton_tax_detail must be specified")
  }
} else { # single
  if (!is.null(in_read_paired_tax_detail)) {
    stop("In single mode, --in_read_paired_tax_detail must not be specified")
  }
  if (is.null(in_read_singleton_tax_detail)) {
    stop("In single mode, --in_read_singleton_tax_detail must be specified")
  }
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
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --read_type:',read_type,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --contigs_taxonomy_classify_method:',contigs_taxonomy_classify_method,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --contigs_function_classify_method:',contigs_function_classify_method,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --reads_function_classify_method:',reads_function_classify_method,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --KO_mode:',KO_mode,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --in_KO_length_ref:',in_KO_length_ref,'\n')
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
expand_taxonomy_columns <- function(dt, CAMI_format_output, id_col) {
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

# Vectorized equivalent of the original by-row loop over taxid_function:
#   vec = strsplit(Strain, " "); if (length(vec) >= 3) sub("^[a-z]__", "", vec[2])
#   else sub("^[a-z]__", "", Strain)
# (strsplit drops trailing empty tokens while tstrsplit keeps them, hence the
#  trailing-space guard in has3.)
# Only taxid and KO are used downstream (Fun_Taxonomy/Lineage/Strain never
# reach any output). Loaded once and shared by the contig/read branches.
taxid_function_cache <- NULL
get_taxid_function <- function(in_taxid_fun) {
  if (is.null(taxid_function_cache)) {
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract mapping relationship of taxid-function from:',in_taxid_fun,'\n')
    taxid_function = fread(in_taxid_fun, sep = "\t", colClasses = "character", col.names = c('taxid','Fun_Taxonomy','Lineage'), header = FALSE, data.table = T)
    taxid_function[, Strain := tstrsplit(Fun_Taxonomy, ";", fixed = TRUE)[[8]]]
    taxid_function[, Strain := gsub("d__unclassified strain", "Unclassified", Strain)]
    taxid_function[, Strain := gsub("root strain", "Root", Strain)]
    tok <- tstrsplit(taxid_function$Strain, " ", fixed = TRUE)
    tok3 <- if (length(tok) >= 3) tok[[3]] else rep(NA_character_, nrow(taxid_function))
    has3 <- !is.na(tok3) & !(tok3 == "" & grepl(" $", taxid_function$Strain))
    taxid_function[, KO := ifelse(has3, sub("^[a-z]__", "", tok[[2]]), sub("^[a-z]__", "", Strain))]
    taxid_function_cache <<- taxid_function[, .(taxid, KO)]
  }
  taxid_function_cache
}

##########################################################################
########################    contig taxononmy   ###########################
##########################################################################
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Contigs taxonomy classify method:',contigs_taxonomy_classify_method,'\n')

if (contigs_taxonomy_classify_method == "CAT_pack"){
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract contigs taxononmy from:',in_contig_tax_detail,'\n')
  seqid_taxid_detail = fread(in_contig_tax_detail, sep = "\t", header = T, data.table = T, fill = T)
  colnames(seqid_taxid_detail) = c("sample_Contig", "classification","reason", "Lineage", "lineage_scores","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  seqid_taxid_detail[, Lineage := gsub("\\*", "", Lineage)]
  seqid_taxid_detail[, taxid := gsub(".*;", "", Lineage)]
  seqid_taxid_detail[is.na(taxid)| taxid == "", taxid := "0"]
  seqid_taxid_detail = seqid_taxid_detail[, .(sample_Contig, taxid)]
} else if (contigs_taxonomy_classify_method == "kraken2"){
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract contigs taxononmy from:',in_contig_tax_detail,'\n')
  seqid_taxid_detail = fread(in_contig_tax_detail, sep = "\t", select = c(2,3), header = FALSE, data.table = T)
  colnames(seqid_taxid_detail) = c("sample_Contig", "taxid")
  seqid_taxid_detail[, taxid := gsub(".* \\(taxid |\\)", "", taxid)]
}

contig_seqid_tax = merge.data.table(seqid_taxid_detail, taxid_taxonomy, by = 'taxid', all.x = T)

if (CAMI_format_output) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Keep taxid and lineage column for contig-level annotation.\n')
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Do not keep taxid and lineage column for contig-level annotation.\n')
}
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total contigs extracted:',nrow(seqid_taxid_detail), '\n')
rm(seqid_taxid_detail)

expand_taxonomy_columns(contig_seqid_tax, CAMI_format_output, "sample_Contig")

#table(contig_seqid_tax$Domain)
#table(contig_seqid_tax$NewTax)

################################################################
####################     read length      ######################
################################################################
read_length = in_read_length %>% {fread(., sep = '\t', header = T, data.table = T)}
ReadLength = max(read_length$avg_len, na.rm = TRUE)
TotalReads = max(read_length$num_seqs, na.rm = TRUE)
TotalBases = max(read_length$sum_len, na.rm = TRUE)

# For long reads, derive separate average lengths for mapped (contig side) and unmapped (read side) reads.
ReadLength_read = ReadLength
ReadLength_contig = ReadLength
if (!is.null(in_read_fna_length) && in_read_fna_length != "" && file.exists(in_read_fna_length)) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Derive mapped/unmapped read lengths from fna_length for long reads.\n')
  read_fna_length_all = fread(in_read_fna_length, sep = '\t', header = T, data.table = T)[, c(1,2), with = FALSE]
  colnames(read_fna_length_all) = c('sample_Read','SeqLen')
  if (nrow(read_fna_length_all) > 0) {
    ReadLength_read = sum(read_fna_length_all$SeqLen) / nrow(read_fna_length_all)
    mapped_reads = TotalReads - nrow(read_fna_length_all)
    if (mapped_reads > 0) {
      ReadLength_contig = (TotalBases - sum(read_fna_length_all$SeqLen)) / mapped_reads
    }
  }
}
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Mean read length for unmapped read side:',ReadLength_read,'; for contig side:',ReadLength_contig,'\n')

##########################################################################
####################     read taxononmy paired     #######################
##########################################################################
empty_result <- data.frame(
  sample_Read = character(),
  Taxonomy    = character(),
  NewTax      = character(),
  Domain      = character(),
  Phylum      = character(),
  Class       = character(),
  Order       = character(),
  Family      = character(),
  Genus       = character(),
  Species     = character(),
  Strain      = character(),
  SeqLen      = integer()
)

process_read_tax_detail <- function(in_file, seq_len, CAMI_format_output, taxid_taxonomy, label) {
  tryCatch({
    dt <- fread(in_file, sep = "\t", select = c(2,3), header = FALSE, data.table = TRUE)
    setnames(dt, c("sample_Read", "taxid"))
    dt[, taxid := gsub(".* \\(taxid |\\)", "", taxid)]
    dt <- merge.data.table(dt, taxid_taxonomy, by = 'taxid', all.x = TRUE)
    if (CAMI_format_output) {
      cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'), paste0('Keep taxid and lineage column for ', label, ' annotation.\n'))
    } else {
      cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'), paste0('Do not keep taxid and lineage column for ', label, ' annotation.\n'))
    }
    expand_taxonomy_columns(dt, CAMI_format_output, "sample_Read")
    dt[, SeqLen := seq_len]
    dt
  }, error = function(e)
    empty_result
  )
}

if (read_type %in% c('paired','paired-rna')) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract paired reads taxonomy from:',in_read_paired_tax_detail,'\n')
  read_paired_seqid_tax = process_read_tax_detail(in_read_paired_tax_detail, ReadLength * 2, CAMI_format_output, taxid_taxonomy, "paired read-level")
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total paired reads extracted:',nrow(read_paired_seqid_tax), '\n')
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract single reads taxonomy from:',in_read_singleton_tax_detail,'\n')
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract reads taxonomy from:',in_read_singleton_tax_detail,'\n')
}
##########################################################################
####################     read taxononmy single    ########################
##########################################################################
read_singleton_seqid_tax =
  tryCatch({
    read_singleton_tax_detail = fread(in_read_singleton_tax_detail, sep = "\t", select = c(2,3), header = FALSE, data.table = T)
    colnames(read_singleton_tax_detail) = c("sample_Read", "taxid")
    read_singleton_tax_detail[, taxid := gsub(".* \\(taxid |\\)", "", taxid)]

    read_singleton_seqid_tax = merge.data.table(read_singleton_tax_detail, taxid_taxonomy, by = 'taxid', all.x = T)

    if (CAMI_format_output) {
      cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Keep taxid and lineage column for singleton read-level annotation.\n')
    } else {
      cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Do not keep taxid and lineage column for singleton read-level annotation.\n')
    }
    expand_taxonomy_columns(read_singleton_seqid_tax, CAMI_format_output, "sample_Read")

    if (read_type %in% c('pacbio-raw','pacbio-corr','pacbio-hifi','nano-raw','nano-corr','nano-hq')) {
      cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Merge per-read fna_length for long-read singleton reads.\n')
      read_fna_length = fread(in_read_fna_length, sep = '\t', header = T, data.table = T)[, c(1,2), with = FALSE]
      colnames(read_fna_length) = c('sample_Read','SeqLen')
      read_singleton_seqid_tax[read_fna_length, on = "sample_Read", SeqLen := i.SeqLen]
    } else {
      read_singleton_seqid_tax[, SeqLen := ReadLength]
    }
    read_singleton_seqid_tax
  }, error = function(e)
    empty_result
  )

##########################################################################
#############    read taxononmy merge single and paired    ###############
##########################################################################

if (read_type %in% c('paired','paired-rna')) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Combine both mapped and unmapped reads taxonomy...','\n')
  read_seqid_tax = rbindlist(list(read_paired_seqid_tax, read_singleton_seqid_tax), use.names = TRUE, fill = TRUE)
  read_seqid_tax[, Avg_fold := 1]
  read_seqid_tax[, SeqBases := Avg_fold * SeqLen]

  rm(read_paired_seqid_tax, read_singleton_seqid_tax)
  #gc()
} else { # single
  read_seqid_tax = read_singleton_seqid_tax
  setDT(read_seqid_tax)
  read_seqid_tax[, Avg_fold := 1]
  read_seqid_tax[, SeqBases := Avg_fold * SeqLen]

  rm(read_singleton_seqid_tax)
  #gc()
}
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total reads extracted:',nrow(read_seqid_tax), '\n')

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
  contigs_fun_detail[, taxid := gsub(".* \\(taxid |\\)", "", taxid)]
  taxid_function = get_taxid_function(in_taxid_fun)
  contigs_fun_detail = merge.data.table(contigs_fun_detail, taxid_function, by = 'taxid', all.x = T)
}

contigs_fun_detail = contigs_fun_detail[!KO %in% '-']

c_fun_detail = contigs_fun_detail[, .(qaccver, KO)]
rm(contigs_fun_detail)

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total annotated ORFs at the contig level:',nrow(c_fun_detail), '\n')
###################################################
###########   ORF function of read   ##############
###################################################

if (reads_function_classify_method == "emapper"){
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract ORF function of reads from:',in_read_fun_detail,'\n')
  read_ORF_KO <- fread(in_read_fun_detail, sep = "\t", header = T, skip = 1, data.table = T)
  read_ORF_KO <- read_ORF_KO[, .SD, .SDcols = c(1, 12)]
  colnames(read_ORF_KO) = c('qaccver','KO')
  read_ORF_KO[, KO := gsub("ko:", "", KO)]
  read_ORF_KO[, KO := gsub(",", "_", KO)]
  read_ORF_KO = read_ORF_KO[!KO %in% '-']
} else if (reads_function_classify_method == "kraken2"){
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract ORF function of reads from:',in_read_fun_detail,'\n')
  read_ORF_KO = fread(in_read_fun_detail, sep = "\t", select = c(2,3), header = FALSE, data.table = T)
  colnames(read_ORF_KO) = c("qaccver", "taxid")
  read_ORF_KO[, taxid := gsub(".* \\(taxid |\\)", "", taxid)]
  taxid_function = get_taxid_function(in_taxid_fun)
  read_ORF_KO = merge.data.table(read_ORF_KO, taxid_function, by = 'taxid', all.x = T)
}

read_ORF_KO = read_ORF_KO[, .(qaccver, KO)]

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total annotated ORFs at the read level:',nrow(read_ORF_KO), '\n')
##########################################################################
########## contig seq id combine from bin and contig taxonomy ############
##########################################################################
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
Contigs_ffn_length[, c("sample", "ORF") := tstrsplit(qaccver, "---", fill = "", keep = 1:2)]
Contigs_ffn_length[, Source := "Contig"]

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract ORF length of reads from:',in_read_ffn_length,'\n')
read_ffn_length = in_read_ffn_length %>% {fread(., sep = '\t', header = T, data.table = T)}
colnames(read_ffn_length) = c('qaccver','ORFLen')
n_na_orflen <- sum(is.na(read_ffn_length$ORFLen))
if (n_na_orflen > 0) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING:',in_read_ffn_length,'contains',n_na_orflen,'NA value(s) in ORFLen; replaced with 0, please check the upstream data','\n')
  read_ffn_length[is.na(ORFLen), ORFLen := 0]
}
# equivalent to strcapture("^(.*)_(\\d+)_(\\d+)_(\\+|-)$", qaccver, ...) keeping
# only sample_Contig (start/end/sign are dropped before any output)
read_ffn_length[, sample_Contig := sub("^(F---|R---)", "", sub("_\\d+_\\d+_[+-]$", "", qaccver))]
read_ffn_length[, sample_Contig := sub("/[12]$", "", sample_Contig)]
read_ffn_length[, Source := "Read"]

##########################################################################
############ taxa abundance including unclassified species ###############
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

c_depth = depth[, .(sample_Contig, SeqLen = ContigLen, Avg_fold, SeqBases = ContigLen * Avg_fold)]
rm(depth)

##############################################################
########### coverm trimmed mean depth based    ###############
##############################################################
c_seqid_tax_depth = c_depth[contig_seqid_tax, on = "sample_Contig"]

setnames(read_seqid_tax, 1, "SeqID")

#colnames(read_seqid_tax)
#colnames(c_seqid_tax_depth)

c_depth_sum = sum(c_depth$SeqBases)
Read_depth_sum = sum(read_seqid_tax$SeqBases)
cr_depth_sum = c_depth_sum + Read_depth_sum

#' Merge contigs and read depth tables and compute relative abundances
#' @param data1      c_seqid_tax_depth, data.table
#' @param data2      read_seqid_tax, data.table
#' @param group_col  character, grouping column name(s), e.g. "Species", "Genus"
#' @param level      character, level abbreviation, e.g. "S", "G"
#' @param c_sum     numeric, total contigs bases (c_depth_sum)
#' @param read_sum   numeric, total read bases (Read_depth_sum)
#' @param read_length    numeric, ReadLength
#'
#' @return data.table with Taxa_Bases_*, Taxa_RelativeAbundance_* columns.
make_cr_depth <- function(data1 = c_seqid_tax_depth,
                          data2 = read_seqid_tax,
                          group_col,
                          level,
                          c_sum   = c_depth_sum,
                          read_sum = Read_depth_sum,
                          read_length  = ReadLength_read,
                          contig_read_length = ReadLength_contig) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Calculating depth at the level:',group_col,'\n')
  # ----- contigs side -----
  c_depth_g <- data1[, .(Taxa_Bases_contig = sum(SeqBases)), keyby = group_col]
  # ----- read side -----
  read_depth_g <- data2[, .(Taxa_Bases_read = sum(SeqBases)), keyby = group_col]

  # ----- combine, widen, add relative abundances -----
  # pivot_wider keeps group keys in order of first appearance in
  # bind_rows(c_depth, read_depth); unique() on the stacked keys reproduces it
  cr_Taxa_abundance <- unique(rbind(c_depth_g[, ..group_col], read_depth_g[, ..group_col]))
  cr_Taxa_abundance[, `:=`(Taxa_Bases_contig = NA_real_, Taxa_Bases_read = NA_real_)]
  cr_Taxa_abundance[c_depth_g, Taxa_Bases_contig := i.Taxa_Bases_contig, on = group_col]
  cr_Taxa_abundance[read_depth_g, Taxa_Bases_read := i.Taxa_Bases_read, on = group_col]
  cr_Taxa_abundance[is.na(Taxa_Bases_contig), Taxa_Bases_contig := 0]
  cr_Taxa_abundance[is.na(Taxa_Bases_read), Taxa_Bases_read := 0]
  cr_Taxa_abundance[, `:=`(
    Taxa_Bases_cr = Taxa_Bases_contig + Taxa_Bases_read,
    Sum_Bases_contig = c_sum,
    Sum_Bases_read = read_sum,
    Sum_Bases_cr = c_sum + read_sum,
    read_length = read_length
  )]
  cr_Taxa_abundance[, Taxa_ReadCount_contig := round(Taxa_Bases_contig / contig_read_length)]
  cr_Taxa_abundance[, Taxa_ReadCount_read := round(Taxa_Bases_read / read_length)]
  cr_Taxa_abundance[, Taxa_ReadCount_cr := Taxa_ReadCount_contig + Taxa_ReadCount_read]
  cr_Taxa_abundance[, `:=`(
    Taxa_RelativeAbundance_contig = 1e2 * Taxa_Bases_contig / Sum_Bases_contig,
    Taxa_RelativeAbundance_read = 1e2 * Taxa_Bases_read / Sum_Bases_read,
    Taxa_RelativeAbundance_cr = 1e2 * Taxa_Bases_cr / Sum_Bases_cr,
    Level = level
  )]
  setnames(cr_Taxa_abundance, 1, "Taxonomy")
  return(cr_Taxa_abundance)
}
cr_Strain_depth  = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, "Strain", "T", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig)
cr_Species_depth = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, "Species","S", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig)
cr_Genus_depth   = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, "Genus",  "G", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig)
cr_Family_depth  = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, "Family", "F", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig)
cr_Order_depth   = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, "Order",  "O", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig)
cr_Class_depth   = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, "Class",  "C", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig)
cr_Phylum_depth  = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, "Phylum", "P", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig)
cr_Domain_depth  = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, "Domain", "D", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig)

All_Taxa_abundance = rbindlist(list(
  cr_Domain_depth,  cr_Phylum_depth,
  cr_Class_depth,   cr_Order_depth,
  cr_Family_depth,  cr_Genus_depth,
  cr_Species_depth, cr_Strain_depth))

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the integrated abundance table of contigs and reads to:',out_Taxa_abundance,'\n')
All_Taxa_abundance[, .(sample = sample_name, Level, Taxonomy,
                   Taxa_ReadCount_contig, Taxa_ReadCount_read, Taxa_ReadCount_cr,
                   Taxa_RelativeAbundance_contig, Taxa_RelativeAbundance_read, Taxa_RelativeAbundance_cr)
][, `:=`(Taxa_ReadCount = Taxa_ReadCount_cr,
         Taxa_RelativeAbundance = Taxa_RelativeAbundance_cr)] %>%
  {fwrite(., paste0(out_Taxa_abundance), row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

rm(cr_Strain_depth, cr_Species_depth, cr_Genus_depth, cr_Family_depth, cr_Order_depth, cr_Class_depth, cr_Phylum_depth, cr_Domain_depth, All_Taxa_abundance)
#gc()

#############################################
############### out_CAMI_format ##############
#############################################

if (CAMI_format_output) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the integrated abundance table of contigs and reads with CAMI format to:',out_CAMI_format,'\n')
  cr_Lineage_Strain_depth  = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, c("Strain","Lineage_Strain"), "strain", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig) %>%
    filter(str_count(Lineage_Strain, ";") >= 8)
  colnames(cr_Lineage_Strain_depth)[2] = "Lineage"
  cr_Lineage_Species_depth = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, c("Species","Lineage_Species"),"species", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig) %>%
    filter(str_count(Lineage_Species, ";") >= 7)
  colnames(cr_Lineage_Species_depth)[2] = "Lineage"
  cr_Lineage_Genus_depth   = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, c("Genus","Lineage_Genus"),  "genus", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig) %>%
    filter(str_count(Lineage_Genus, ";") >= 6)
  colnames(cr_Lineage_Genus_depth)[2] = "Lineage"
  cr_Lineage_Family_depth  = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, c("Family","Lineage_Family"), "family", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig) %>%
    filter(str_count(Lineage_Family, ";") >= 5)
  colnames(cr_Lineage_Family_depth)[2] = "Lineage"
  cr_Lineage_Order_depth   = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, c("Order","Lineage_Order"),  "order", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig) %>%
    filter(str_count(Lineage_Order, ";") >= 4)
  colnames(cr_Lineage_Order_depth)[2] = "Lineage"
  cr_Lineage_Class_depth   = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, c("Class","Lineage_Class"),  "class", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig) %>%
    filter(str_count(Lineage_Class, ";") >= 3)
  colnames(cr_Lineage_Class_depth)[2] = "Lineage"
  cr_Lineage_Phylum_depth  = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, c("Phylum","Lineage_Phylum"), "phylum", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig) %>%
    filter(str_count(Lineage_Phylum, ";") >= 2)
  colnames(cr_Lineage_Phylum_depth)[2] = "Lineage"
  cr_Lineage_Domain_depth  = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, c("Domain","Lineage_Domain"), "superkingdom", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig) %>%
    filter(str_count(Lineage_Domain, ";") >= 1)
  colnames(cr_Lineage_Domain_depth)[2] = "Lineage"
  cr_Lineage_Root_depth    = make_cr_depth(c_seqid_tax_depth, read_seqid_tax, c("taxid_Root","Lineage_Root"), "root", c_depth_sum, Read_depth_sum, ReadLength_read, ReadLength_contig)
  cr_Lineage_Root_depth[, Taxonomy := gsub("1", "root", gsub("0", "unclassified", Taxonomy))]
  colnames(cr_Lineage_Root_depth)[2] = "Lineage"
  All_Taxa_Lineage_depth = rbindlist(list(
    cr_Lineage_Root_depth,
    cr_Lineage_Domain_depth,  cr_Lineage_Phylum_depth,
    cr_Lineage_Class_depth,   cr_Lineage_Order_depth,
    cr_Lineage_Family_depth,  cr_Lineage_Genus_depth,
    cr_Lineage_Species_depth, cr_Lineage_Strain_depth), use.names = TRUE, fill = TRUE)
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
  All_Taxa_Lineage_depth = All_Taxa_Lineage_depth[, .(taxid, Level, Lineage, Taxonomy, Taxa_RelativeAbundance_cr)]
  colnames(All_Taxa_Lineage_depth) = c("@@TAXID","RANK","TAXPATH","TAXPATHSN","PERCENTAGE")

  header_lines <- c(
    "# Taxonomic Profiling Output",
    paste0("@SampleID:",sample_name),
    paste0("@Version:", Sys.getenv("METALINK_VERSION", "unknown")),
    "@Ranks:superkingdom|phylum|class|order|family|genus|species|strain",
    paste0("@TaxonomyID:",in_taxid_tax),
    "@__program__:metalink cr",
    "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE"
  )
  #out_CAMI_format <- "output_cami_profile.txt"
  writeLines(header_lines, out_CAMI_format)
  write.table(All_Taxa_Lineage_depth, out_CAMI_format,sep = "\t",row.names = F,col.names = F,append = T,quote = F)

  rm(cr_Lineage_Strain_depth, cr_Lineage_Species_depth, cr_Lineage_Genus_depth, cr_Lineage_Family_depth, cr_Lineage_Order_depth, cr_Lineage_Class_depth, cr_Lineage_Phylum_depth, cr_Lineage_Domain_depth, cr_Lineage_Root_depth, All_Taxa_Lineage_depth)
  # the CAMI lineage helper columns are not used downstream of this block
  read_seqid_tax[, c("taxid_Root", lineage_cols) := NULL]
  c_seqid_tax_depth[, c("taxid_Root", lineage_cols) := NULL]
  #gc()
}
#############################################
######## ORF length for contigs and read ##########
#############################################

c_ffn_length_depth = c_depth[Contigs_ffn_length, on = "sample_Contig"]
c_ffn_length_depth = c_ffn_length_depth[, .(sample_Contig, qaccver, ORFLen, Avg_fold, sample, Source)]

read_ffn_length_depth = read_ffn_length[, .(sample_Contig, qaccver, ORFLen, Source)]
read_ffn_length_depth[, `:=`(Avg_fold = 1, sample = sample_name)]
setcolorder(read_ffn_length_depth, c("sample_Contig","qaccver","ORFLen","Avg_fold","sample","Source"))

rm(c_depth, read_ffn_length, Contigs_ffn_length)
#gc()

#colnames(c_taxonomy)
#colnames(read_taxonomy)
#############################################
######## function abundance for contigs ##########
#############################################

seperate_emapper_ko_depth <- function(data) {
  setDT(data)
  data[, .(KO = unlist(strsplit(KO, "_"))), by = .(ORF, Avg_fold)
  ][, Duplicate_Annotations_num := .N, by = ORF]
}


if (KO_mode == "split_multiple_KOs") {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'For contigs: Split the depth of an ORF with multiple KO annotations by evenly distributing the total value across all assigned KOs.','\n')
  c_ORF_KO_depth = c_fun_detail[, .(qaccver, KO)][c_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  c_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  c_ORF_KO_depth = seperate_emapper_ko_depth(c_ORF_KO_depth[, .(ORF = qaccver, Avg_fold, KO)])
  setnames(c_ORF_KO_depth, "ORF", "qaccver")
} else if (KO_mode == "keep_multiple_KOs") {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'For contigs: Ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.','\n')
  c_ORF_KO_depth = c_fun_detail[, .(qaccver, KO)][c_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  c_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  c_ORF_KO_depth[, Duplicate_Annotations_num := 1]
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'For contigs: Ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.','\n')
  c_ORF_KO_depth = c_fun_detail[, .(qaccver, KO)][c_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  c_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  c_ORF_KO_depth[, Duplicate_Annotations_num := 1]
  c_ORF_KO_depth = c_ORF_KO_depth[!str_detect(KO, "_")]
}
c_ORF_KO_depth[, Scaled_Avg_fold := Avg_fold / Duplicate_Annotations_num]
c_ORF_KO_depth[c_ffn_length_depth, on = "qaccver",
               `:=`(sample_Contig = i.sample_Contig,
                    ORFLen = i.ORFLen,
                    sample = i.sample,
                    Source = i.Source)]

#############################################
######## function abundance for read ########
#############################################

if (KO_mode == "split_multiple_KOs") {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'For reads: Split the depth of an ORF with multiple KO annotations by evenly distributing the total value across all assigned KOs.','\n')
  read_ORF_KO_depth = read_ORF_KO[, .(qaccver, KO)][read_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  read_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  read_ORF_KO_depth = seperate_emapper_ko_depth(read_ORF_KO_depth[, .(ORF = qaccver, Avg_fold, KO)])
  setnames(read_ORF_KO_depth, "ORF", "qaccver")
} else if (KO_mode == "keep_multiple_KOs") {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'For reads: Keep original multiple KO annotations.','\n')
  read_ORF_KO_depth = read_ORF_KO[, .(qaccver, KO)][read_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  read_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  read_ORF_KO_depth[, Duplicate_Annotations_num := 1]
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'For reads: Ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.','\n')
  read_ORF_KO_depth = read_ORF_KO[, .(qaccver, KO)][read_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  read_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  read_ORF_KO_depth[, Duplicate_Annotations_num := 1]
  read_ORF_KO_depth = read_ORF_KO_depth[!str_detect(KO, "_")]
}
read_ORF_KO_depth[read_ffn_length_depth, on = "qaccver",
                  `:=`(sample_Contig = i.sample_Contig,
                       ORFLen = i.ORFLen,
                       sample = i.sample,
                       Source = i.Source)]

# Normalize read-side ORF depth by KO reference length: a read-derived ORF is a
# fragment of the full-length KO gene, so counting it as 1 mixes gene length
# into the abundance. Weighting by ORFLen / (3 * len_median) (len_median is in
# amino acids) puts the sum on the same coverage scale as the contig side
# (covered bases / gene length). KOs absent from the reference table (e.g.
# emapper-only annotations or on-the-fly multi-KO strings) fall back to the
# global median of len_median over all accessions.
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract mapping relationship of accession-KO reference length from:',in_KO_length_ref,'\n')
ko_len_ref = fread(in_KO_length_ref, sep = "\t", header = TRUE,
                   select = c("accession", "len_median"), data.table = TRUE)
ko_len_ref = ko_len_ref[!is.na(len_median) & len_median > 0]
global_len_median = median(ko_len_ref$len_median)
read_ORF_KO_depth[ko_len_ref, on = c(KO = "accession"), KO_ref_len := i.len_median]
n_fallback = read_ORF_KO_depth[, sum(is.na(KO_ref_len))]
read_ORF_KO_depth[is.na(KO_ref_len), KO_ref_len := global_len_median]
read_ORF_KO_depth[, Scaled_Avg_fold := ORFLen / (3 * KO_ref_len) / Duplicate_Annotations_num]
read_ORF_KO_depth[, KO_ref_len := NULL]
if (nrow(read_ORF_KO_depth) > 0) {
  cat(sprintf('[%s] Read-side ORF depth normalized by KO reference length; %d of %d read-ORF row(s) (%.2f%%) have no accession match in %s and fall back to the global median len_median = %.1f aa.\n',
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              n_fallback, nrow(read_ORF_KO_depth),
              1e2 * n_fallback / nrow(read_ORF_KO_depth),
              in_KO_length_ref, global_len_median))
  if (n_fallback / nrow(read_ORF_KO_depth) > 0.9) {
    warning(sprintf("%.2f%% of read-ORF rows fell back to the global median KO length; check that --in_KO_length_ref matches the version of the function annotation database.",
                    1e2 * n_fallback / nrow(read_ORF_KO_depth)))
  }
}

# attach the taxonomy columns (only those used downstream) by sample_Contig
need_tax_cols <- unique(c("Taxonomy", "NewTax", "Domain", rank_vec))
if (CAMI_format_output) need_tax_cols <- c(need_tax_cols, "taxid", "Lineage")
c_ORF_KO_depth[contig_seqid_tax, on = "sample_Contig",
               (need_tax_cols) := mget(paste0("i.", need_tax_cols))]
read_ORF_KO_depth[read_seqid_tax, on = .(sample_Contig = SeqID),
                  (need_tax_cols) := mget(paste0("i.", need_tax_cols))]

rm(c_ffn_length_depth, read_ffn_length_depth, c_fun_detail, read_ORF_KO)
#gc()

#colnames(c_ORF_KO_depth)
#colnames(read_ORF_KO_depth)
#dplyr::all_equal(c_ORF_KO_depth_2, c_ORF_KO_depth)
#############################################
######## function abundance for cr #########
#############################################
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Combine ORFs KO annotations of contigs and reads...','\n')

# Row order is canonicalized to C-locale radix order. Grouped sums of
# Scaled_Avg_fold / ORFRA can differ from the pre-0.2.1 implementation in the
# last printed digit only, because floating-point addition is not associative;
# all other digits are identical. Per-row values (scaleOrfBases, ORFRA,
# ORFRA_Bases) are computed with the same expressions and are row-order
# independent.
cr_ORF_KO_depth = rbindlist(list(c_ORF_KO_depth, read_ORF_KO_depth), use.names = TRUE, fill = TRUE)
cr_ORF_KO_depth[, scaleOrfBases := ORFLen * Avg_fold / Duplicate_Annotations_num]
cr_ORF_KO_depth <- cr_ORF_KO_depth[ORFLen >= as.numeric(minORFLen)]

ORF_depth_sum = cr_ORF_KO_depth[, .(sumOrfBases = sum(scaleOrfBases),
                                    sumOrfDepth = sum(Scaled_Avg_fold)), keyby = sample]

ORF_eggnog_depth = cr_ORF_KO_depth
ORF_eggnog_depth[ORF_depth_sum, on = "sample",
                 `:=`(sumOrfBases = i.sumOrfBases, sumOrfDepth = i.sumOrfDepth)]
ORF_eggnog_depth[, ORFRA := 1e6 * Scaled_Avg_fold/sumOrfDepth]
ORF_eggnog_depth[, ORFRA_Bases := 1e6 * scaleOrfBases/sumOrfBases]
ORF_eggnog_depth[, NewTax := gsub("unclassified", "Unclassified", NewTax)]
ORF_eggnog_depth[, NewTax := gsub("New", "Unclassified", NewTax)]
setorderv(ORF_eggnog_depth, c("sample_Contig", "qaccver"), na.last = TRUE)

rm(cr_ORF_KO_depth, ORF_depth_sum)
#gc()

#############################################
############### summary for KO ##############
#############################################

KO_eggnog_depth = ORF_eggnog_depth[, .(KO_ORFNum = .N,
                                       KO_SumCoverage = sum(Scaled_Avg_fold),
                                       KO_RelativeAbundance = sum(ORFRA),
                                       ORFLen_median = floor(median(ORFLen, na.rm = TRUE))), keyby = .(sample, KO)]

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the table of KOs mean coverage and relative abundance to:',out_KO_eggnog_depth,'\n')
KO_eggnog_depth %>%
  {fwrite(., paste0(out_KO_eggnog_depth), row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

rm(KO_eggnog_depth)

#############################################
############ summary for KO_Taxa ############
#############################################

KO_Taxa_SumCoverage = ORF_eggnog_depth[, .(Taxa_KO_ORF_Num = .N,
                                     Taxa_KO_SumCoverage = sum(Scaled_Avg_fold),
                                     Taxa_KO_RelativeAbundance = sum(ORFRA)),
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
  rm(KO_Taxa_SumCoverage_wide, KO_Taxa_RelativeAbundance_wide)
} else {
  warning("No taxon were retained, maybe min_KO_num_in_one_taxon is too large.")
  fwrite(KO_Taxa_SumCoverage, out_KO_Taxa_SumCoverage_wide, row.names = FALSE, col.names = T,sep = '\t', quote = FALSE)
  fwrite(KO_Taxa_SumCoverage, out_KO_Taxa_RelativeAbundance_wide, row.names = FALSE, col.names = T,sep = '\t', quote = FALSE)
}

rm(KO_Taxa_SumCoverage)

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
                                            Taxa_KO_SumCoverage = sum(Scaled_Avg_fold),
                                            Taxa_KO_RelativeAbundance = sum(ORFRA)),
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
    rm(KO_Taxa_SumCoverage_wide, KO_Taxa_RelativeAbundance_wide)
  } else {
    warning("No taxon were retained, maybe min_KO_num_in_one_taxon is too large.")
    fwrite(KO_Taxa_SumCoverage, sprintf("%s.%s.tsv", sub("\\.tsv$", "", out_KO_Taxa_SumCoverage_wide), rank), row.names = FALSE, col.names = T,sep = '\t', quote = FALSE)
    fwrite(KO_Taxa_SumCoverage, sprintf("%s.%s.tsv", sub("\\.tsv$", "", out_KO_Taxa_RelativeAbundance_wide), rank), row.names = FALSE, col.names = T,sep = '\t', quote = FALSE)
  }
  rm(KO_Taxa_SumCoverage)
  return(KO_Taxa_SumCoverage_out)
}
for (rk in rank_vec) {
  make_KO_Taxa_tables(rk, ORF_eggnog_depth, out_KO_Taxa_SumCoverage_wide, out_KO_Taxa_RelativeAbundance_wide, as.numeric(min_KO_num_in_one_taxon))
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

#############################################
########### summary for ORF_Tax_KO ##########
#############################################

# summary
ORF_eggnog_depth[, KO_annotate := ifelse(KO %in% c("Unclassified"),"Unclassified","Classified")]
ORF_Tax_KO_summary <- summary_with_all(ORF_eggnog_depth,
                                       c("Domain", "NewTax", "KO_annotate"),
                                       "ORF_Num", "Total_KO", c("Domain", "NewTax"))

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the summary table of ORFs taxonomy and function to:',summary_ORF_Tax_KO,'\n')
#summary = "{sample}.contig_Tax_summary.tsv
ORF_Tax_KO_summary %>%
  {fwrite(., summary_ORF_Tax_KO, row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

rm(ORF_Tax_KO_summary)

#############################################
############# summary for ORF KO ############
#############################################

ORF_KO_summary <- summary_with_all(ORF_eggnog_depth,
                                   c("Domain", "KO_annotate"),
                                   "ORF_Num", "Total_KO", "Domain")

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the summary table of ORFs KO annotations to:',summary_ORF_KO,'\n')
#summary = "{sample}.contig_Tax_summary.tsv"
ORF_KO_summary %>%
  {fwrite(., summary_ORF_KO, row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

rm(ORF_KO_summary)
#gc()

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

rm(ORF_Tax_Source_summary)
#gc()

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

rm(ORF_Source_summary)

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

c_tax_for_sum = contig_seqid_tax[, .(Domain, NewTax)]
c_tax_for_sum[, NewTax := gsub("unclassified", "Unclassified", NewTax)]
c_tax_for_sum[, NewTax := gsub("New", "Unclassified", NewTax)]

c_Tax_summary <- summary_with_all(c_tax_for_sum,
                                  c("Domain", "NewTax"),
                                  "Contig_Num_c", "Total_Contig_c", "Domain", "Proportion_Contig_c")

rm(c_tax_for_sum)

r_tax_for_sum = read_seqid_tax[, .(Domain, NewTax)]
r_tax_for_sum[, NewTax := gsub("unclassified", "Unclassified", NewTax)]
r_tax_for_sum[, NewTax := gsub("New", "Unclassified", NewTax)]

r_Tax_summary <- summary_with_all(r_tax_for_sum,
                                  c("Domain", "NewTax"),
                                  "Read_Num_r", "Total_Read_r", "Domain", "Proportion_r")

rm(r_tax_for_sum)

##########################################################################
########### cr Tax summary with unified SeqBases scale ###############
##########################################################################
c_Seq_Bases_by_tax <- c_seqid_tax_depth[, .(Seq_Bases_c = sum(SeqBases)), keyby = .(Domain, NewTax)]
read_Seq_Bases_by_tax <- read_seqid_tax[, .(Seq_Bases_read = sum(SeqBases)), keyby = .(Domain, NewTax)]
cr_tax_by_SeqBases <- full_join(
  c_Seq_Bases_by_tax,
  read_Seq_Bases_by_tax,
  by = c("Domain", "NewTax")
) %>%
  as.data.frame()
# mutate(across(starts_with("Seq_Bases"), ~replace_na(., 0)))
for (cl in c("Seq_Bases_c", "Seq_Bases_read")) {
  v <- cr_tax_by_SeqBases[[cl]]
  v[is.na(v)] <- 0
  cr_tax_by_SeqBases[[cl]] <- v
}
cr_tax_by_SeqBases$Seq_Bases_cr <- cr_tax_by_SeqBases$Seq_Bases_c + cr_tax_by_SeqBases$Seq_Bases_read

cr_Tax_summary <- cr_tax_by_SeqBases %>%
  bind_rows(
    cr_tax_by_SeqBases %>%
      mutate(Domain = "All") %>%
      group_by(Domain, NewTax) %>%
      summarise(
        Seq_Bases_cr = sum(Seq_Bases_cr),
        Seq_Bases_c = sum(Seq_Bases_c),
        Seq_Bases_read = sum(Seq_Bases_read),
        .groups = "drop"
      )
  ) %>%
  group_by(Domain) %>%
  mutate(
    Total_Seq_cr = sum(Seq_Bases_cr),
    Proportion_cr = ifelse(Total_Seq_cr == 0, 0, Seq_Bases_cr / Total_Seq_cr)
  ) %>%
  ungroup() %>%
  as.data.frame()
cr_Tax_summary = cr_Tax_summary %>%
  {mutate(., NewTax = gsub("unclassified", "Unclassified", .$NewTax))} %>%
  {mutate(., NewTax = gsub("New", "Unclassified", .$NewTax))}

all_Tax_summary = cr_Tax_summary %>%
  merge(., c_Tax_summary, all = T, by = c("Domain","NewTax")) %>%
  merge(., r_Tax_summary, all = T, by = c("Domain","NewTax"))
# mutate(across(where(is.numeric), ~replace_na(., 0)))
for (cl in names(all_Tax_summary)[vapply(all_Tax_summary, is.numeric, logical(1))]) {
  v <- all_Tax_summary[[cl]]
  v[is.na(v)] <- 0
  all_Tax_summary[[cl]] <- v
}


cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the summary table of contigs and reads taxonomy to:',summary_Tax,'\n')
#summary = "{sample}.cr_Tax_summary.tsv"
all_Tax_summary %>%
  {fwrite(., summary_Tax, row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

mode <- "cr"
if (!is.null(out_mode_flag)) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Exporting execution mode:',mode,'to',out_mode_flag,'\n')
  dir.create(dirname(out_mode_flag), showWarnings = FALSE, recursive = TRUE)
  writeLines(mode, out_mode_flag)
}
