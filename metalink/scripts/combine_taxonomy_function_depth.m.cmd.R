###
#' @Date: 2023-06-14 16:48:43
#' @LastEditors: LiuyangLi
#' @LastEditTime: 2026-09-07
#' @FilePath: /liliuyang/R/Rscript.genome/combine_taxonomy_function_depth.m.cmd.R
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

MAGs_function_classify_method = "kraken2"
script_path <- sub("--file=", "", commandArgs()[grepl("--file=", commandArgs())])
if (length(script_path) == 0 || script_path == "") {
  script_path <- sys.frame(1)$ofile
}
source(file.path(dirname(script_path), "RLib.r"))
# GForce-optimized gsum rounds differently from base R sum() (which dplyr uses
# per group); fall back to base evaluation so grouped sums are bit-identical to
# the dplyr implementation
options(datatable.optimize = 0L)
in_taxid_fun=NULL

# ===================== optparse parameter definition =====================
option_list <- list(
  make_option("--in_taxid_fun", type="character", default=NULL,
              help="taxid_tax table, e.g., taxid_fun_scale.tsv"),
  make_option("--MAGs_taxonomy_classify_method", type="character", default=NULL,
              help="gtdbtk or BAT"),
  make_option("--MAGs_function_classify_method", type="character", default=NULL,
              help="kraken2 or emapper"),
  make_option("--in_MAGs_tax_detail", type="character", default=NULL,
              help="in_MAGs_tax_detail.tsv"),
  make_option("--in_MAGs_fun_detail", type="character", default=NULL,
              help="in_MAGs_fun_detail.tsv"),
  make_option("--in_MAGs_ffn_length", type="character", default=NULL,
              help="in_MAGs_ffn_length.tsv"),
  make_option("--in_contig_depth", type="character", default=NULL,
              help="in_contig_depth.tsv"),
  make_option("--minORFLen", type="numeric", default=0,
              help="Min threshhold for ORF length"),
  make_option("--sample_name", type="character", default=NULL,
              help="sample_name, e.g., 'SRR000001' (default: test)"),
  make_option("--KO_mode", type="character", default="ignore_multiple_KOs",
              help="KO handling mode: 'keep_multiple_KOs', 'ignore_multiple_KOs', or 'split_multiple_KOs' (default: ignore_multiple_KOs)"),
  make_option("--min_KO_num_in_one_taxon", type="character", default=NULL,
              help="min_KO_num_in_one_taxon: min number of KO genes (rows) a taxon (column) must contain"),
  make_option("--detail", type="character", default=NULL,
              help="output detail"),
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
  make_option("--out_KO_eggnog_depth", type="character", default=NULL,
              help="output out_KO_eggnog_depth"),
  make_option("--out_Taxa_abundance", type="character", default=NULL,
              help="out_Taxa_abundance"),
  make_option("--out_KO_Taxa_SumCoverage_wide", type="character", default=NULL,
              help="out_KO_Taxa_SumCoverage_wide"),
  make_option("--out_KO_Taxa_RelativeAbundance_wide", type="character", default=NULL,
              help="out_KO_Taxa_RelativeAbundance_wide"),
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
in_taxid_fun <- opt$in_taxid_fun
MAGs_taxonomy_classify_method <- opt$MAGs_taxonomy_classify_method
MAGs_function_classify_method <- opt$MAGs_function_classify_method
in_MAGs_tax_detail <- opt$in_MAGs_tax_detail
in_MAGs_fun_detail <- opt$in_MAGs_fun_detail
in_MAGs_ffn_length <- opt$in_MAGs_ffn_length
in_contig_depth <- opt$in_contig_depth
minORFLen <- opt$minORFLen
sample_name <- opt$sample_name
KO_mode <- opt$KO_mode
min_KO_num_in_one_taxon <- opt$min_KO_num_in_one_taxon
detail <- opt$detail
# make_KO_Taxa_tables() writes the per-rank long tables into dirname(detail);
# the directory is no longer auto-created via a declared output
if (!is.null(detail)) dir.create(dirname(detail), showWarnings = FALSE, recursive = TRUE)
summary_Tax <- opt$summary_Tax
summary_ORF_KO <- opt$summary_ORF_KO
summary_ORF_Tax_KO <- opt$summary_ORF_Tax_KO
summary_ORF_Source <- opt$summary_ORF_Source
summary_ORF_Tax_Source <- opt$summary_ORF_Tax_Source
out_KO_eggnog_depth <- opt$out_KO_eggnog_depth
out_Taxa_abundance <- opt$out_Taxa_abundance
out_KO_Taxa_SumCoverage_wide <- opt$out_KO_Taxa_SumCoverage_wide
out_KO_Taxa_RelativeAbundance_wide <- opt$out_KO_Taxa_RelativeAbundance_wide
CAMI_format_output <- opt$CAMI_format_output
out_CAMI_format <- opt$out_CAMI_format
rank <- opt$rank
output_detail <- opt$output_detail
verbose <- opt$verbose
out_mode_flag <- opt$out_mode_flag

if (is.null(sample_name)) {
  stop("--sample_name must be specified")
}

if (is.null(KO_mode) || !KO_mode %in% c("keep_multiple_KOs","ignore_multiple_KOs","split_multiple_KOs")) {
  stop("--KO_mode must be 'keep_multiple_KOs', 'ignore_multiple_KOs', or 'split_multiple_KOs'")
}

if (is.null(MAGs_taxonomy_classify_method) || !MAGs_taxonomy_classify_method %in% c("gtdbtk","BAT")) {
  stop("--MAGs_taxonomy_classify_method must be 'gtdbtk' or 'BAT'")
}

if (is.null(MAGs_function_classify_method) || !MAGs_function_classify_method %in% c("kraken2","emapper")) {
  stop("--MAGs_function_classify_method must be 'kraken2' or 'emapper'")
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
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --MAGs_taxonomy_classify_method:',MAGs_taxonomy_classify_method,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --MAGs_function_classify_method:',MAGs_function_classify_method,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --KO_mode:',KO_mode,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --minORFLen:',minORFLen,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --CAMI_format_output:',CAMI_format_output,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --output_detail:',output_detail,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --rank:',rank,'\n')

##########################################################################
########################     MAG taxononmy     ###########################
##########################################################################
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract MAGs taxononmy from:',in_MAGs_tax_detail,'\n')
MAG_tax_detail = fread(in_MAGs_tax_detail,  sep = "\t", header = T, data.table = F, fill = T)
colnames(MAG_tax_detail)[1] = "genome"
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'The number of MAGs used for analysis',nrow(MAG_tax_detail), '\n')

# The MAG-level taxonomy table is small (one row per genome), so the original
# dplyr chains are kept as-is here.
MAG_tax_detail = MAG_tax_detail %>%
  rename_GTDB_taxonomy_dataframe(.) %>%
  rename_strain_dataframe(.)

gtdbtk_detail = MAG_tax_detail %>%
  as.data.frame(.) %>%
  mutate(., Strain = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus,";",.$Species,";",.$Strain))) %>%
  mutate(., Species= (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus,";",.$Species))) %>%
  mutate(., Genus  = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus))) %>%
  mutate(., Family = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family))) %>%
  mutate(., Order  = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order))) %>%
  mutate(., Class  = (paste0(.$Domain,";",.$Phylum,";",.$Class))) %>%
  mutate(., Phylum = (paste0(.$Domain,";",.$Phylum))) %>%
  mutate(., relative_abundance = 1e2 * relative_abundance) %>%
  mutate(., Taxonomy = Strain)

if (CAMI_format_output) {
  if (MAGs_taxonomy_classify_method == "gtdbtk"){
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'It is unavailable to keep taxid and lineage column when MAGs_taxonomy_classify_method is set as gtdbtk.\n')
    MAG_seqid_tax = gtdbtk_detail %>%
      dplyr::select(., c("genome","Taxonomy","NewTax","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"))
    genome_tax_depth = gtdbtk_detail %>%
      mutate(., sample = sample_name) %>%
      select(c("sample", "Taxonomy", "Total scaffolds","coverm_abundance","relative_abundance",
              "Domain", "Phylum", "Class", "Order", "Family", "Genus","Species", "Strain", "genome"))
    colnames(genome_tax_depth)[1:5] = c("sample", "Taxonomy","Taxa_ContigNum","coverm_abundance","relative_abundance")
  } else if (MAGs_taxonomy_classify_method == "BAT"){
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Keep taxid and lineage column for MAG-level annotation.\n')
    gtdbtk_detail = as.data.table(gtdbtk_detail)
    gtdbtk_detail$Lineage = gtdbtk_detail$Lineage
    gtdbtk_detail$taxid = as.character(gtdbtk_detail$taxid)
    gtdbtk_detail[, c("taxid_Root","taxid_Domain", "taxid_Phylum", "taxid_Class", "taxid_Order", "taxid_Family", "taxid_Genus", "taxid_Species", "taxid_Strain") := split_lineage_fill(Lineage, n = 9)]
    gtdbtk_detail = gtdbtk_detail %>%
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
    MAG_seqid_tax = gtdbtk_detail %>%
      dplyr::select(., c("genome","taxid","Lineage","Taxonomy","NewTax",
                         "taxid_Root","taxid_Domain", "taxid_Phylum", "taxid_Class", "taxid_Order", "taxid_Family", "taxid_Genus", "taxid_Species", "taxid_Strain",
                         "Lineage_Root","Lineage_Domain", "Lineage_Phylum", "Lineage_Class", "Lineage_Order", "Lineage_Family", "Lineage_Genus", "Lineage_Species", "Lineage_Strain",
                         "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"))
    genome_tax_depth = gtdbtk_detail %>%
      mutate(., sample = sample_name) %>%
      select(c("sample", "Taxonomy", "Total scaffolds","coverm_abundance","relative_abundance",
               "taxid_Root","taxid_Domain", "taxid_Phylum", "taxid_Class", "taxid_Order", "taxid_Family", "taxid_Genus", "taxid_Species", "taxid_Strain",
               "Lineage_Root","Lineage_Domain", "Lineage_Phylum", "Lineage_Class", "Lineage_Order", "Lineage_Family", "Lineage_Genus", "Lineage_Species", "Lineage_Strain",
               "Domain", "Phylum", "Class", "Order", "Family", "Genus","Species", "Strain", "genome"))
    colnames(genome_tax_depth)[1:5] = c("sample", "Taxonomy","Taxa_ContigNum","coverm_abundance","relative_abundance")
  }
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Do not keep taxid and lineage column for MAG-level annotation.\n')
  MAG_seqid_tax = gtdbtk_detail %>%
    dplyr::select(., c("genome","Taxonomy","NewTax","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"))
  genome_tax_depth = gtdbtk_detail %>%
    mutate(., sample = sample_name) %>%
    select(c("sample", "Taxonomy", "Total scaffolds","coverm_abundance","relative_abundance",
            "Domain", "Phylum", "Class", "Order", "Family", "Genus","Species", "Strain", "genome"))
  colnames(genome_tax_depth)[1:5] = c("sample", "Taxonomy","Taxa_ContigNum","coverm_abundance","relative_abundance")
}
setDT(MAG_seqid_tax)
setDT(genome_tax_depth)

##########################################################################
###########   ORF function of MAG    ##############
##########################################################################

if (MAGs_function_classify_method == "emapper"){
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract ORF function of MAGs from:',in_MAGs_fun_detail,'\n')
  MAGs_fun_detail <- fread(in_MAGs_fun_detail, sep = "\t", header = T, skip = 1, data.table = T)
  MAGs_fun_detail <- MAGs_fun_detail[, .SD, .SDcols = c(1, 12)]
  colnames(MAGs_fun_detail) = c('qaccver','KO')
  MAGs_fun_detail[, KO := gsub("ko:", "", KO)]
  MAGs_fun_detail[, KO := gsub(",", "_", KO)]
} else if (MAGs_function_classify_method == "kraken2"){
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract ORF function of MAGs from:',in_MAGs_fun_detail,'\n')
  MAGs_fun_detail = fread(in_MAGs_fun_detail, sep = "\t", select = c(2,3), header = FALSE, data.table = T)
  colnames(MAGs_fun_detail) = c("qaccver", "taxid")
  MAGs_fun_detail[, taxid := gsub(".* \\(taxid |\\)", "", taxid)]
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
  MAGs_fun_detail = MAGs_fun_detail %>%
    merge.data.table(., taxid_function, by = 'taxid', all.x = T)
  rm(taxid_function)
}

MAGs_fun_detail[, c("sample_Contig") := tstrsplit(qaccver, "_\\d*$")]
MAGs_fun_detail[, contig := tstrsplit(sample_Contig, "---")[[2]]]
MAGs_fun_detail[, sample := sample_name]
MAGs_fun_detail[, sample_Contig := paste0(sample_name, "---", contig)]
MAGs_fun_detail[, c("genome", "ORF") := tstrsplit(qaccver, "---", fill = "", keep = 1:2)]
MAGs_fun_detail[, Source := "mMAG"]
MAGs_fun_detail <- MAGs_fun_detail[genome %in% MAG_seqid_tax$genome]

MAGs_fun_detail <- MAGs_fun_detail[!KO %in% '-']

m_fun_detail = MAGs_fun_detail[, .(qaccver, KO)]

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total annotated ORFs at the MAG level:',nrow(m_fun_detail), '\n')
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
# 使用 tstrsplit() 将 sample_Contig 列拆分成两列
depth[, c("sample", "contig") := tstrsplit(sample_Contig, "---", fill = "", keep = 1:2)]

depth[, c("sample", "contig") := NULL]

##############################################################
########### coverm trimmed mean depth based    ###############
##############################################################

#' Summarise genome_tax_depth by a chosen taxonomic rank
#' @param data2      genome_tax_depth, data.table
#' @param group_col  character, grouping column name, e.g. "Species", "Genus"
#' @param level      character, level abbreviation, e.g. "S", "G"
#' @return data.table with Taxa_ContigNum/Taxa_SumCoverage/Taxa_RelativeAbundance columns.
make_mag_depth <- function(data2 = genome_tax_depth, group_col, level) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Calculating depth at the level:',group_col,'\n')
  mag_Taxa_abundance <- data2[, .(Taxa_ContigNum = sum(Taxa_ContigNum),
                              Taxa_SumCoverage = sum(coverm_abundance),
                              Taxa_RelativeAbundance = sum(relative_abundance)),
                          keyby = c("sample", group_col)]
  mag_Taxa_abundance[, Level := level]
  setnames(mag_Taxa_abundance, 2, "Taxonomy")
  return(mag_Taxa_abundance)
}
Domain_depth_genome  = make_mag_depth(genome_tax_depth, "Domain",  "D")
Phylum_depth_genome  = make_mag_depth(genome_tax_depth, "Phylum",  "P")
Class_depth_genome   = make_mag_depth(genome_tax_depth, "Class",   "C")
Order_depth_genome   = make_mag_depth(genome_tax_depth, "Order",   "O")
Family_depth_genome  = make_mag_depth(genome_tax_depth, "Family",  "F")
Genus_depth_genome   = make_mag_depth(genome_tax_depth, "Genus",   "G")
Species_depth_genome = make_mag_depth(genome_tax_depth, "Species", "S")
Strain_depth_genome  = make_mag_depth(genome_tax_depth, "Strain",  "T")

All_Taxa_abundance_genome = rbindlist(list(
  Domain_depth_genome,Phylum_depth_genome,
  Class_depth_genome,Order_depth_genome,
  Family_depth_genome,Genus_depth_genome,
  Species_depth_genome,Strain_depth_genome))

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the integrated abundance table of MAGs to:',out_Taxa_abundance,'\n')
All_Taxa_abundance_genome[, .(sample, Level, Taxonomy,
                          Taxa_ContigNum,Taxa_SumCoverage,Taxa_RelativeAbundance)] %>%
  {fwrite(., paste0(out_Taxa_abundance), row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

rm(Domain_depth_genome, Phylum_depth_genome, Class_depth_genome, Order_depth_genome,
   Family_depth_genome, Genus_depth_genome, Species_depth_genome, Strain_depth_genome,
   All_Taxa_abundance_genome)

#############################################
############### out_CAMI_format ##############
#############################################

if (CAMI_format_output) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the integrated abundance table of contigs with CAMI format to:',out_CAMI_format,'\n')

  m_Lineage_Strain_depth = genome_tax_depth[, .(Taxa_ContigNum = sum(Taxa_ContigNum),
                                                Taxa_SumCoverage = sum(coverm_abundance),
                                                Taxa_RelativeAbundance = sum(relative_abundance)),
                                            keyby = .(Strain, Lineage_Strain)]
  m_Lineage_Strain_depth[, Level := "strain"]
  m_Lineage_Strain_depth = m_Lineage_Strain_depth[str_count(Lineage_Strain, ";") >= 8]
  colnames(m_Lineage_Strain_depth)[2] = "Lineage"

  m_Lineage_Species_depth = genome_tax_depth[, .(Taxa_ContigNum = sum(Taxa_ContigNum),
                                                 Taxa_SumCoverage = sum(coverm_abundance),
                                                 Taxa_RelativeAbundance = sum(relative_abundance)),
                                             keyby = .(Species, Lineage_Species)]
  m_Lineage_Species_depth[, Level := "species"]
  m_Lineage_Species_depth = m_Lineage_Species_depth[str_count(Lineage_Species, ";") >= 7]
  colnames(m_Lineage_Species_depth)[2] = "Lineage"

  m_Lineage_Genus_depth = genome_tax_depth[, .(Taxa_ContigNum = sum(Taxa_ContigNum),
                                               Taxa_SumCoverage = sum(coverm_abundance),
                                               Taxa_RelativeAbundance = sum(relative_abundance)),
                                           keyby = .(Genus, Lineage_Genus)]
  m_Lineage_Genus_depth[, Level := "genus"]
  m_Lineage_Genus_depth = m_Lineage_Genus_depth[str_count(Lineage_Genus, ";") >= 6]
  colnames(m_Lineage_Genus_depth)[2] = "Lineage"

  m_Lineage_Family_depth = genome_tax_depth[, .(Taxa_ContigNum = sum(Taxa_ContigNum),
                                                Taxa_SumCoverage = sum(coverm_abundance),
                                                Taxa_RelativeAbundance = sum(relative_abundance)),
                                            keyby = .(Family, Lineage_Family)]
  m_Lineage_Family_depth[, Level := "family"]
  m_Lineage_Family_depth = m_Lineage_Family_depth[str_count(Lineage_Family, ";") >= 5]
  colnames(m_Lineage_Family_depth)[2] = "Lineage"

  m_Lineage_Order_depth = genome_tax_depth[, .(Taxa_ContigNum = sum(Taxa_ContigNum),
                                               Taxa_SumCoverage = sum(coverm_abundance),
                                               Taxa_RelativeAbundance = sum(relative_abundance)),
                                           keyby = .(Order, Lineage_Order)]
  m_Lineage_Order_depth[, Level := "order"]
  m_Lineage_Order_depth = m_Lineage_Order_depth[str_count(Lineage_Order, ";") >= 4]
  colnames(m_Lineage_Order_depth)[2] = "Lineage"

  m_Lineage_Class_depth = genome_tax_depth[, .(Taxa_ContigNum = sum(Taxa_ContigNum),
                                               Taxa_SumCoverage = sum(coverm_abundance),
                                               Taxa_RelativeAbundance = sum(relative_abundance)),
                                           keyby = .(Class, Lineage_Class)]
  m_Lineage_Class_depth[, Level := "class"]
  m_Lineage_Class_depth = m_Lineage_Class_depth[str_count(Lineage_Class, ";") >= 3]
  colnames(m_Lineage_Class_depth)[2] = "Lineage"

  m_Lineage_Phylum_depth = genome_tax_depth[, .(Taxa_ContigNum = sum(Taxa_ContigNum),
                                                Taxa_SumCoverage = sum(coverm_abundance),
                                                Taxa_RelativeAbundance = sum(relative_abundance)),
                                            keyby = .(Phylum, Lineage_Phylum)]
  m_Lineage_Phylum_depth[, Level := "phylum"]
  m_Lineage_Phylum_depth = m_Lineage_Phylum_depth[str_count(Lineage_Phylum, ";") >= 2]
  colnames(m_Lineage_Phylum_depth)[2] = "Lineage"

  m_Lineage_Domain_depth = genome_tax_depth[, .(Taxa_ContigNum = sum(Taxa_ContigNum),
                                                Taxa_SumCoverage = sum(coverm_abundance),
                                                Taxa_RelativeAbundance = sum(relative_abundance)),
                                            keyby = .(Domain, Lineage_Domain)]
  m_Lineage_Domain_depth[, Level := "superkingdom"]
  m_Lineage_Domain_depth = m_Lineage_Domain_depth[str_count(Lineage_Domain, ";") >= 1]
  colnames(m_Lineage_Domain_depth)[2] = "Lineage"

  m_Lineage_Root_depth = genome_tax_depth[, .(Taxa_ContigNum = sum(Taxa_ContigNum),
                                              Taxa_SumCoverage = sum(coverm_abundance),
                                              Taxa_RelativeAbundance = sum(relative_abundance)),
                                          keyby = .(taxid_Root, Lineage_Root)]
  m_Lineage_Root_depth[, Level := "root"]
  colnames(m_Lineage_Root_depth)[2] = "Lineage"

  colnames(m_Lineage_Root_depth)[1] = "Taxonomy"
  colnames(m_Lineage_Domain_depth)[1] = "Taxonomy"
  colnames(m_Lineage_Phylum_depth)[1] = "Taxonomy"
  colnames(m_Lineage_Class_depth)[1] = "Taxonomy"
  colnames(m_Lineage_Order_depth)[1] = "Taxonomy"
  colnames(m_Lineage_Family_depth)[1] = "Taxonomy"
  colnames(m_Lineage_Genus_depth)[1] = "Taxonomy"
  colnames(m_Lineage_Species_depth)[1] = "Taxonomy"
  colnames(m_Lineage_Strain_depth)[1] = "Taxonomy"

  m_Lineage_Root_depth[, Taxonomy := gsub("1", "root", gsub("0", "unclassified", Taxonomy))]

  All_Taxa_Lineage_depth = rbindlist(list(
    m_Lineage_Root_depth,
    m_Lineage_Domain_depth,  m_Lineage_Phylum_depth,
    m_Lineage_Class_depth,   m_Lineage_Order_depth,
    m_Lineage_Family_depth,  m_Lineage_Genus_depth,
    m_Lineage_Species_depth, m_Lineage_Strain_depth))
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
  All_Taxa_Lineage_depth = All_Taxa_Lineage_depth[, .(taxid, Level, Lineage, Taxonomy, Taxa_RelativeAbundance)]
  colnames(All_Taxa_Lineage_depth) = c("@@TAXID","RANK","TAXPATH","TAXPATHSN","PERCENTAGE")

  header_lines <- c(
    "# Taxonomic Profiling Output",
    paste0("@SampleID:",sample_name),
    paste0("@Version:", Sys.getenv("METALINK_VERSION", "unknown")),
    "@Ranks:superkingdom|phylum|class|order|family|genus|species|strain",
    paste0("@TaxonomyID:",in_MAGs_tax_detail),
    "@__program__:metalink m",
    "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE"
  )
  #out_CAMI_format <- "output_cami_profile.txt"
  writeLines(header_lines, out_CAMI_format)
  write.table(All_Taxa_Lineage_depth, out_CAMI_format,sep = "\t",row.names = F,col.names = F,append = T,quote = F)

  rm(m_Lineage_Strain_depth, m_Lineage_Species_depth, m_Lineage_Genus_depth,
     m_Lineage_Family_depth, m_Lineage_Order_depth, m_Lineage_Class_depth,
     m_Lineage_Phylum_depth, m_Lineage_Domain_depth, m_Lineage_Root_depth,
     All_Taxa_Lineage_depth)
}
#############################################
############ function abundance ###############
#############################################

seperate_emapper_ko_depth <- function(data) {
  setDT(data)
  data[, .(KO = unlist(strsplit(KO, "_"))), by = .(ORF, Avg_fold)
  ][, Duplicate_Annotations_num := .N, by = ORF]
}

m_taxonomy = MAG_seqid_tax

mags_ffn_length = in_MAGs_ffn_length %>% {fread(., sep = '\t', header = T, data.table = T)}
colnames(mags_ffn_length) = c('qaccver','ORFLen')
n_na_orflen <- sum(is.na(mags_ffn_length$ORFLen))
if (n_na_orflen > 0) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING:',in_MAGs_ffn_length,'contains',n_na_orflen,'NA value(s) in ORFLen; replaced with 0, please check the upstream data','\n')
  mags_ffn_length[is.na(ORFLen), ORFLen := 0]
}
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total ORFs at the MAG level:',nrow(mags_ffn_length), '\n')
mags_ffn_length[, c("sample_Contig") := tstrsplit(qaccver, "_\\d*$")]
mags_ffn_length[, contig := tstrsplit(sample_Contig, "---")[[2]]]
mags_ffn_length[, sample := sample_name]
mags_ffn_length[, sample_Contig := paste0(sample_name, "---", contig)]
mags_ffn_length[, c("genome", "ORF") := tstrsplit(qaccver, "---", fill = "", keep = 1:2)]
mags_ffn_length <- mags_ffn_length[genome %in% MAG_seqid_tax$genome]
m_ffn_length_depth = mags_ffn_length[, .(sample_Contig, qaccver, ORFLen, sample, genome)]
m_ffn_length_depth[depth, on = "sample_Contig", Avg_fold := i.Avg_fold]

# Row order is canonicalized to C-locale radix order. In ignore/keep modes
# every row has Scaled_Avg_fold == Avg_fold, so all downstream grouped sums are
# exactly order-independent. In split mode, grouped sums of Scaled_Avg_fold
# (= Avg_fold / N) and ORFRA can differ from the pre-0.2.1 implementation in
# the last printed digit only, because floating-point addition is not
# associative; all other digits are identical.

if (KO_mode == "split_multiple_KOs") {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Split the depth of an ORF with multiple KO annotations by evenly distributing the total value across all assigned KOs.','\n')
  m_ORF_KO_depth = m_fun_detail[m_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  m_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  m_ORF_KO_depth = seperate_emapper_ko_depth(m_ORF_KO_depth[, .(ORF = qaccver, Avg_fold, KO)])
  setnames(m_ORF_KO_depth, "ORF", "qaccver")
} else if (KO_mode == "keep_multiple_KOs") {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Keep original multiple KO annotations.','\n')
  m_ORF_KO_depth = m_fun_detail[m_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  m_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  m_ORF_KO_depth[, Duplicate_Annotations_num := 1]
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.','\n')
  m_ORF_KO_depth = m_fun_detail[m_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  m_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  m_ORF_KO_depth[, Duplicate_Annotations_num := 1]
  m_ORF_KO_depth = m_ORF_KO_depth[!str_detect(KO, "_")]
}
m_ORF_KO_depth[, Scaled_Avg_fold := Avg_fold / Duplicate_Annotations_num]
m_ORF_KO_depth[m_ffn_length_depth, on = "qaccver",
               `:=`(sample_Contig = i.sample_Contig,
                    ORFLen = i.ORFLen,
                    sample = i.sample,
                    genome = i.genome)]
setDT(m_ORF_KO_depth)
need_tax_cols <- unique(c("Taxonomy", "NewTax", "Domain", rank_vec))
if (CAMI_format_output) need_tax_cols <- c(need_tax_cols, "taxid", "Lineage")
m_ORF_KO_depth[m_taxonomy, on = "genome",
               (need_tax_cols) := mget(paste0("i.", need_tax_cols))]
setorderv(m_ORF_KO_depth, c("sample_Contig", "qaccver"), na.last = TRUE)

rm(mags_ffn_length, m_ffn_length_depth, m_fun_detail, MAGs_fun_detail)

#############################################
######## function abundance for mcr #########
#############################################

m_ORF_KO_depth[, scaleOrfBases := ORFLen*Scaled_Avg_fold]
m_ORF_KO_depth <- m_ORF_KO_depth[ORFLen >= as.numeric(minORFLen)]

ORF_depth_sum = m_ORF_KO_depth[, .(sumOrfBases = sum(scaleOrfBases),
                                   sumOrfDepth = sum(Scaled_Avg_fold)), keyby = sample]

ORF_eggnog_depth = m_ORF_KO_depth
ORF_eggnog_depth[ORF_depth_sum, on = "sample",
                 `:=`(sumOrfBases = i.sumOrfBases, sumOrfDepth = i.sumOrfDepth)]
ORF_eggnog_depth[, ORFRA := 1e6 * Scaled_Avg_fold/sumOrfDepth]
ORF_eggnog_depth[, ORFRA_Bases := 1e6 * scaleOrfBases/sumOrfBases]
ORF_eggnog_depth[, Source := "mMAG"]
ORF_eggnog_depth[, NewTax := gsub("unclassified", "Unclassified", NewTax)]
ORF_eggnog_depth[, NewTax := gsub("New", "Unclassified", NewTax)]

rm(m_ORF_KO_depth, ORF_depth_sum)

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
summary_with_all <- function(dt, count_cols, count_name, total_name, total_by) {
  rest_cols <- count_cols[-1]
  s1 <- dt[, .(N = .N), keyby = count_cols]
  s2 <- dt[, .(N = .N), keyby = rest_cols]
  s2[, (count_cols[1]) := "All"]
  setcolorder(s2, c(count_cols[1], count_cols[-1]))
  out <- rbindlist(list(s1, s2))
  setnames(out, "N", count_name)
  out[, (total_name) := sum(get(count_name)), by = total_by]
  out[, Proportion := get(count_name) / get(total_name)]
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
#summary = "{sample}.MAG_Tax_summary.tsv"
ORF_KO_summary %>%
  {fwrite(., summary_ORF_KO, row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

rm(ORF_KO_summary)

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

#############################################
############## summary for Taxa #############
#############################################

MAG_seqid_tax[, NewTax := gsub("New", "Unclassified", NewTax)]

MAG_Tax_summary <- summary_with_all(MAG_seqid_tax,
                                    c("Domain", "NewTax"),
                                    "MAG_Num", "Total_MAG", "Domain")

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the summary table of MAGs taxonomy to:',summary_Tax,'\n')
#summary = "{sample}.MAG_Tax_summary.tsv"
MAG_Tax_summary %>%
  {fwrite(., summary_Tax, row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

mode <- "m"
if (!is.null(out_mode_flag)) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Exporting execution mode:',mode,'to',out_mode_flag,'\n')
  dir.create(dirname(out_mode_flag), showWarnings = FALSE, recursive = TRUE)
  writeLines(mode, out_mode_flag)
}
