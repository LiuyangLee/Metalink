###
#' @Date: 2026-09-20
#' @LastEditors: LiuyangLi
#' @LastEditTime: 2026-09-20
#' @FilePath: /liliuyang/R/Rscript.genome/combine_taxonomy_function_depth.pmr_v.cmd.R
#' @Description: integrate vMAG taxonomy (CAT/BAT lineage via taxid), function
#'   (Virus_annotation_results KO), and project-level coverm depth for pmr --mags v.
#'   Derived from combine_taxonomy_function_depth.pmr_m.cmd.R; the vMAG-side
#'   column handling mirrors combine_taxonomy_function_depth.vm.cmd.R.
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

script_path <- sub("--file=", "", commandArgs()[grepl("--file=", commandArgs())])
if (length(script_path) == 0 || script_path == "") {
  script_path <- sys.frame(1)$ofile
}
source(file.path(dirname(script_path), "RLib.r"))
# GForce-optimized gsum rounds differently from base R sum() (which dplyr uses
# per group); fall back to base evaluation so grouped sums are bit-identical to
# the dplyr implementation
options(datatable.optimize = 0L)

# ===================== optparse parameter definition =====================
option_list <- list(
  make_option("--in_taxid_tax", type="character", default=NULL,
              help="taxid_tax table, e.g., taxid_tax_scale.tsv"),
  make_option("--in_vMAGs_tax_detail", type="character", default=NULL,
              help="CheckV-filtered vMAG_summary_info.filter.tsv"),
  make_option("--in_vMAGs_ffn_length", type="character", default=NULL,
              help="project-level Virus_annotation_results.tsv (KO + ORF start/end)"),
  make_option("--in_MAGs_depth", type="character", default=NULL,
              help="in_MAGs_depth.tsv (coverm genome output)"),
  make_option("--in_contig_depth", type="character", default=NULL,
              help="in_contig_depth.tsv (coverm contig dense output)"),
  make_option("--minORFLen", type="numeric", default=0,
              help="Min threshhold for ORF length"),
  make_option("--sample_name", type="character", default=NULL,
              help="sample_name, e.g., 'SRR000001' (default: test)"),
  make_option("--KO_mode", type="character", default="ignore_multiple_KOs",
              help="KO handling mode: 'keep_multiple_KOs', 'ignore_multiple_KOs', or 'split_multiple_KOs' (default: ignore_multiple_KOs)"),
  make_option("--min_KO_num_in_one_taxon", type="numeric", default=0,
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
  make_option("--out_MAG_tax_depth", type="character", default=NULL,
              help="out_MAG_tax_depth"),
  make_option("--CAMI_format_output", type="character", default="FALSE",
              help="Enable CAMI output format (default: FALSE)"),
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
in_vMAGs_tax_detail <- opt$in_vMAGs_tax_detail
in_vMAGs_ffn_length <- opt$in_vMAGs_ffn_length
in_MAGs_depth <- opt$in_MAGs_depth
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
out_MAG_tax_depth <- opt$out_MAG_tax_depth
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
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --KO_mode:',KO_mode,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --minORFLen:',minORFLen,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --CAMI_format_output:',CAMI_format_output,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --output_detail:',output_detail,'\n')
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Parameter --rank:',rank,'\n')

##########################################################################
########################     vMAG abundance     ##########################
##########################################################################

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract vMAGs abundance from:',in_MAGs_depth,'\n')
Ref_coverm_trimmed_mean = fread(in_MAGs_depth,  sep = "\t", header = T, data.table = F, fill = T)
colnames(Ref_coverm_trimmed_mean) = c("genome","coverm_abundance")
n_na_coverm <- sum(is.na(Ref_coverm_trimmed_mean$coverm_abundance))
if (n_na_coverm > 0) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING:',in_MAGs_depth,'contains',n_na_coverm,'NA value(s) in coverm_abundance; replaced with 0, please check the upstream data','\n')
  Ref_coverm_trimmed_mean$coverm_abundance[is.na(Ref_coverm_trimmed_mean$coverm_abundance)] <- 0
}
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total vMAGs extracted:',nrow(Ref_coverm_trimmed_mean), '\n')

Ref_coverm_trimmed_mean = Ref_coverm_trimmed_mean %>%
  {.[!(.$genome %in% 'unmapped'), ]} %>%
  mutate(coverm_abundance = as.numeric(coverm_abundance)) %>%
  mutate(total_abundance = sum(coverm_abundance, na.rm = TRUE)) %>%
  mutate(
    relative_abundance = if_else(
      total_abundance == 0,
      0,         # 注意：必须与另一侧同类型
      coverm_abundance / total_abundance
    )
  ) %>%
  select(-total_abundance)

##########################################################################
########################     taxid taxonomy     ##########################
##########################################################################

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract mapping relationship of taxid-taxonomy from:',in_taxid_tax,'\n')
taxid_taxonomy = fread(in_taxid_tax, sep = "\t", colClasses = "character",
                       col.names = c('taxid','Taxonomy','Lineage'), header = FALSE, data.table = T)

##########################################################################
########################    virus MAGs taxononmy   #######################
##########################################################################

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract virus MAGs taxononmy from:',in_vMAGs_tax_detail,'\n')
virus_tax_detail = fread(in_vMAGs_tax_detail,  sep = "\t", #select = c(1,2,4),
                         header = T, data.table = T, fill = T)
colnames(virus_tax_detail)[1] = 'genome'
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'The number of virus MAGs used for analysis',nrow(virus_tax_detail), '\n')

# the sample-level coverm_abundance/relative_abundance columns carried over by
# the project-level vMAG_summary_info.tsv are stale; fresh depth comes from
# the pmr_v coverm output joined below
virus_tax_detail = virus_tax_detail %>%
  select(-any_of(c("coverm_abundance","relative_abundance")))

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

# raw vMAG taxonomy + coverm depth (kept for out_MAG_tax_depth): intersection
# of the taxonomy table and the coverm genome output, the same genome-set
# semantics as the mMAG side (MAG_tax_depth_raw) and pmr_vm
virus_tax_depth_raw <- Ref_coverm_trimmed_mean %>%
  inner_join(virus_tax_detail, by = 'genome')

if (!is.null(out_MAG_tax_depth) && out_MAG_tax_depth != "") {
  virus_tax_depth_raw %>%
    {fwrite(., out_MAG_tax_depth, row.names = FALSE, col.names = T,
            sep = '\t', quote = FALSE)}
}

############

virus_tax_depth_raw = as.data.table(virus_tax_depth_raw)
virus_tax_depth_raw[, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain") := split_lineage_fill(Taxonomy, n = 8)]

virus_tax_depth_raw = virus_tax_depth_raw %>%
  as.data.frame(.) %>%
  rename_GTDB_taxonomy_dataframe(.) %>%
  rename_strain_dataframe(.)

virus_seqid_tax = virus_tax_depth_raw %>%
  as.data.frame(.) %>%
  mutate(., Strain = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus,";",.$Species,";",.$Strain))) %>%
  mutate(., Species= (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus,";",.$Species))) %>%
  mutate(., Genus  = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family,";",.$Genus))) %>%
  mutate(., Family = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order,";",.$Family))) %>%
  mutate(., Order  = (paste0(.$Domain,";",.$Phylum,";",.$Class,";",.$Order))) %>%
  mutate(., Class  = (paste0(.$Domain,";",.$Phylum,";",.$Class))) %>%
  mutate(., Phylum = (paste0(.$Domain,";",.$Phylum))) %>%
  mutate(., relative_abundance = 1e2 * relative_abundance) %>%
  mutate(., Taxonomy = Strain) %>%
  mutate(., sample = sample_name)   # MODIFIED: use sample_name parameter directly

# number of member scaffolds of each vMAG cluster (Taxa_ContigNum), counted by
# splitting the comma-separated scaffold column of vMAG_summary_info
if ("scaffold" %in% colnames(virus_tax_detail)) {
  v_scaf_num = virus_tax_detail[, .(Taxa_ContigNum = lengths(strsplit(scaffold, ",", fixed = TRUE))), by = .(genome)]
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING: column scaffold not found in',in_vMAGs_tax_detail,'; Taxa_ContigNum of vMAGs set to 0','\n')
  v_scaf_num = data.table(genome = virus_tax_detail$genome, Taxa_ContigNum = 0L)
}

if (CAMI_format_output) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Keep taxid and lineage column for vMAG-level annotation.\n')
  virus_seqid_tax = as.data.table(virus_seqid_tax)
  virus_seqid_tax$Lineage = virus_seqid_tax$Lineage
  virus_seqid_tax$taxid = as.character(virus_seqid_tax$taxid)
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
  v_genome_tax_depth = virus_seqid_tax %>%
    mutate(., sample = sample_name) %>%   # MODIFIED: use sample_name parameter directly
    left_join(v_scaf_num %>% as.data.frame(.), by = "genome") %>%
    select(c("sample", "Taxonomy", "Taxa_ContigNum","coverm_abundance","relative_abundance",
             "taxid_Root","taxid_Domain", "taxid_Phylum", "taxid_Class", "taxid_Order", "taxid_Family", "taxid_Genus", "taxid_Species", "taxid_Strain",
             "Lineage_Root","Lineage_Domain", "Lineage_Phylum", "Lineage_Class", "Lineage_Order", "Lineage_Family", "Lineage_Genus", "Lineage_Species", "Lineage_Strain",
             "Domain", "Phylum", "Class", "Order", "Family", "Genus","Species", "Strain", "genome"))
  virus_seqid_tax = virus_seqid_tax %>%
    dplyr::select(., c("genome","taxid","Lineage","Taxonomy","NewTax",
                       "taxid_Root","taxid_Domain", "taxid_Phylum", "taxid_Class", "taxid_Order", "taxid_Family", "taxid_Genus", "taxid_Species", "taxid_Strain",
                       "Lineage_Root","Lineage_Domain", "Lineage_Phylum", "Lineage_Class", "Lineage_Order", "Lineage_Family", "Lineage_Genus", "Lineage_Species", "Lineage_Strain",
                       "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"))
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Do not keep taxid and lineage column for vMAG-level annotation.\n')
  v_genome_tax_depth = virus_seqid_tax %>%
    mutate(., sample = sample_name) %>%   # MODIFIED: use sample_name parameter directly
    left_join(v_scaf_num %>% as.data.frame(.), by = "genome") %>%
    select(c("sample", "Taxonomy", "Taxa_ContigNum","coverm_abundance","relative_abundance",
            "Domain", "Phylum", "Class", "Order", "Family", "Genus","Species", "Strain", "genome"))
  virus_seqid_tax = virus_seqid_tax %>%
    dplyr::select(., c("genome","Taxonomy","NewTax","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"))
}

##########################################################################
###########   ORF function and length of vMAG    ##############
##########################################################################

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract ORF length and function of vMAGs from:',in_vMAGs_ffn_length,'\n')
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
vMAGs_ffn_infor[, sample := sample_name]
vMAGs_ffn_infor[, Source := "vMAG"]
vMAGs_ffn_infor <- vMAGs_ffn_infor[genome %in% virus_seqid_tax$genome]

v_ffn_length = vMAGs_ffn_infor[, .(qaccver, ORFLen, sample_Contig, sample, contig, genome, Source)]
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total ORFs at the vMAG level:',nrow(v_ffn_length), '\n')

virus_fun_detail = vMAGs_ffn_infor[!KO %in% '', .(qaccver, KO)]
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total annotated ORFs at the vMAG level:',nrow(virus_fun_detail), '\n')

##########################################################################
############ taxa abundance including unclassified species ###############
##########################################################################
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract depth of contigs from:',in_contig_depth,'\n')
depth = in_contig_depth %>% {fread(., sep = '\t', header = T, data.table = T)}
colnames(depth)=c("genome_Contig","Avg_fold")
n_na_depth <- sum(is.na(depth$Avg_fold))
if (n_na_depth > 0) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING:',in_contig_depth,'contains',n_na_depth,'NA value(s) in Avg_fold; replaced with 0, please check the upstream data','\n')
  depth[is.na(Avg_fold), Avg_fold := 0]
}
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total contigs extracted:',nrow(depth), '\n')
# vMAG contig headers can carry viral boundary marks (||..., |provirus_...);
# strip them so that ORF and depth tables share the sample---contig key
depth[, genome_Contig := tstrsplit(genome_Contig, "\\|\\|", fixed = FALSE, type.convert = TRUE)[[1]]]
depth[, genome_Contig := tstrsplit(genome_Contig, "\\|provirus_", fixed = FALSE, type.convert = TRUE)[[1]]]
depth[, genome_Contig := tstrsplit(genome_Contig, "\\|vs2region_", fixed = FALSE, type.convert = TRUE)[[1]]]
depth[, contig := tstrsplit(genome_Contig, "---")[[2]]]
depth[, sample_Contig := paste0(sample_name, "---", contig)]

depth[, c("contig", "genome_Contig") := NULL]

######## filter ORFs to contigs present in the depth table #######
v_ffn_length = v_ffn_length[sample_Contig %in% depth$sample_Contig]
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Total ORFs for vMAGs with depth:',nrow(v_ffn_length), '\n')

# left join of v_ffn_length (x) onto depth (y) by sample_Contig,
# replacing merge(., depth, all.x = T, by = "sample_Contig")
v_ffn_length_depth = depth[
  v_ffn_length[, .(sample_Contig, qaccver, ORFLen, sample, genome)],
  on = "sample_Contig"
][, .(sample_Contig, qaccver, ORFLen, sample, Avg_fold, genome)]

rm(vMAGs_ffn_infor, v_ffn_length)

seperate_emapper_ko_depth <- function(data) {
  setDT(data)
  data[, .(KO = unlist(strsplit(KO, "_"))), by = .(ORF, Avg_fold)
  ][, Duplicate_Annotations_num := .N, by = ORF]
}

if (KO_mode == "split_multiple_KOs") {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Split the depth of an ORF with multiple KO annotations by evenly distributing the total value across all assigned KOs.','\n')
  v_ORF_KO_depth = virus_fun_detail[v_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  v_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  v_ORF_KO_depth = seperate_emapper_ko_depth(v_ORF_KO_depth[, .(ORF = qaccver, Avg_fold, KO)])
  setnames(v_ORF_KO_depth, "ORF", "qaccver")
} else if (KO_mode == "keep_multiple_KOs") {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Keep original multiple KO annotations.','\n')
  v_ORF_KO_depth = virus_fun_detail[v_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  v_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  v_ORF_KO_depth[, Duplicate_Annotations_num := 1]
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Ignore the ORFs with multiple KO annotations and only use the ORFs with single KO annotation.','\n')
  v_ORF_KO_depth = virus_fun_detail[v_ffn_length_depth[, .(qaccver, Avg_fold)], on = "qaccver"]
  v_ORF_KO_depth[is.na(KO), KO := "Unclassified"]
  v_ORF_KO_depth[, Duplicate_Annotations_num := 1]
  v_ORF_KO_depth = v_ORF_KO_depth[!str_detect(KO, "_")]
}
v_ORF_KO_depth[, Scaled_Avg_fold := Avg_fold / Duplicate_Annotations_num]
# attach the remaining per-ORF columns, then the vMAG taxonomy columns
v_ORF_KO_depth[v_ffn_length_depth, on = "qaccver",
               `:=`(sample_Contig = i.sample_Contig,
                    ORFLen = i.ORFLen,
                    sample = i.sample,
                    genome = i.genome)]
setorderv(v_ORF_KO_depth, c("sample_Contig", "qaccver"), na.last = TRUE)

v_taxonomy = as.data.table(virus_seqid_tax)
need_tax_cols <- unique(c("Taxonomy", "NewTax", "Domain", rank_vec))
if (CAMI_format_output) need_tax_cols <- c(need_tax_cols, "taxid", "Lineage")
v_ORF_KO_depth[v_taxonomy, on = "genome",
               (need_tax_cols) := mget(paste0("i.", need_tax_cols))]

rm(v_ffn_length_depth, virus_fun_detail)
#gc()

##########################################################################
############## coverm trimmed mean depth based    ########################
##########################################################################

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Calculating depth at the level: Domain','\n')
Domain_depth_genome = v_genome_tax_depth %>%
  group_by(sample, Domain) %>%
  summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
            Taxa_SumCoverage = sum(coverm_abundance),
            Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
  ) %>%
  as.data.frame(.) %>%
  mutate(., Level = "D")

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Calculating depth at the level: Phylum','\n')
Phylum_depth_genome = v_genome_tax_depth %>%
  group_by(sample, Phylum) %>%
  summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
            Taxa_SumCoverage = sum(coverm_abundance),
            Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
  ) %>%
  as.data.frame(.) %>%
  mutate(., Level = "P")

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Calculating depth at the level: Class','\n')
Class_depth_genome = v_genome_tax_depth %>%
  group_by(sample, Class) %>%
  summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
            Taxa_SumCoverage = sum(coverm_abundance),
            Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
  ) %>%
  as.data.frame(.) %>%
  mutate(., Level = "C")

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Calculating depth at the level: Order','\n')
Order_depth_genome = v_genome_tax_depth %>%
  group_by(sample, Order) %>%
  summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
            Taxa_SumCoverage = sum(coverm_abundance),
            Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
  ) %>%
  as.data.frame(.) %>%
  mutate(., Level = "O")

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Calculating depth at the level: Family','\n')
Family_depth_genome = v_genome_tax_depth %>%
  group_by(sample, Family) %>%
  summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
            Taxa_SumCoverage = sum(coverm_abundance),
            Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
  ) %>%
  as.data.frame(.) %>%
  mutate(., Level = "F")

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Calculating depth at the level: Genus','\n')
Genus_depth_genome = v_genome_tax_depth %>%
  group_by(sample, Genus) %>%
  summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
            Taxa_SumCoverage = sum(coverm_abundance),
            Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
  ) %>%
  as.data.frame(.) %>%
  mutate(., Level = "G")

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Calculating depth at the level: Species','\n')
Species_depth_genome = v_genome_tax_depth %>%
  group_by(sample, Species) %>%
  summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
            Taxa_SumCoverage = sum(coverm_abundance),
            Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
  ) %>%
  as.data.frame(.) %>%
  mutate(., Level = "S")

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Calculating depth at the level: Strain','\n')
Strain_depth_genome = v_genome_tax_depth %>%
  group_by(sample, Strain) %>%
  summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
            Taxa_SumCoverage = sum(coverm_abundance),
            Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
  ) %>%
  as.data.frame(.) %>%
  mutate(., Level = "T")

colnames(Domain_depth_genome)[2] = "Taxonomy"
colnames(Phylum_depth_genome)[2] = "Taxonomy"
colnames(Class_depth_genome)[2] = "Taxonomy"
colnames(Order_depth_genome)[2] = "Taxonomy"
colnames(Family_depth_genome)[2] = "Taxonomy"
colnames(Genus_depth_genome)[2] = "Taxonomy"
colnames(Species_depth_genome)[2] = "Taxonomy"
colnames(Strain_depth_genome)[2] = "Taxonomy"

All_Taxa_abundance_genome = rbind(
  Domain_depth_genome,Phylum_depth_genome,
  Class_depth_genome,Order_depth_genome,
  Family_depth_genome,Genus_depth_genome,
  Species_depth_genome,Strain_depth_genome)

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the integrated abundance table of vMAGs to:',out_Taxa_abundance,'\n')
All_Taxa_abundance_genome %>%
  select(.,c("sample", "Level", "Taxonomy",
             "Taxa_ContigNum","Taxa_SumCoverage","Taxa_RelativeAbundance")) %>%
  {fwrite(., paste0(out_Taxa_abundance), row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

#############################################
############### out_CAMI_format ##############
#############################################

if (CAMI_format_output) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the integrated abundance table of contigs with CAMI format to:',out_CAMI_format,'\n')

  v_Lineage_Strain_depth = v_genome_tax_depth %>%
    group_by(across(all_of(c("Strain","Lineage_Strain")))) %>%
    summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
              Taxa_SumCoverage = sum(coverm_abundance),
              Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
    ) %>%
    as.data.frame(.) %>%
    mutate(., Level = "strain") %>%
    filter(str_count(Lineage_Strain, ";") >= 8)
  colnames(v_Lineage_Strain_depth)[2] = "Lineage"

  v_Lineage_Species_depth = v_genome_tax_depth %>%
    group_by(across(all_of(c("Species","Lineage_Species")))) %>%
    summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
              Taxa_SumCoverage = sum(coverm_abundance),
              Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
    ) %>%
    as.data.frame(.) %>%
    mutate(., Level = "species") %>%
    filter(str_count(Lineage_Species, ";") >= 7)
  colnames(v_Lineage_Species_depth)[2] = "Lineage"

  v_Lineage_Genus_depth = v_genome_tax_depth %>%
    group_by(across(all_of(c("Genus","Lineage_Genus")))) %>%
    summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
              Taxa_SumCoverage = sum(coverm_abundance),
              Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
    ) %>%
    as.data.frame(.) %>%
    mutate(., Level = "genus") %>%
    filter(str_count(Lineage_Genus, ";") >= 6)
  colnames(v_Lineage_Genus_depth)[2] = "Lineage"

  v_Lineage_Family_depth = v_genome_tax_depth %>%
    group_by(across(all_of(c("Family","Lineage_Family")))) %>%
    summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
              Taxa_SumCoverage = sum(coverm_abundance),
              Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
    ) %>%
    as.data.frame(.) %>%
    mutate(., Level = "family") %>%
    filter(str_count(Lineage_Family, ";") >= 5)
  colnames(v_Lineage_Family_depth)[2] = "Lineage"

  v_Lineage_Order_depth = v_genome_tax_depth %>%
    group_by(across(all_of(c("Order","Lineage_Order")))) %>%
    summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
              Taxa_SumCoverage = sum(coverm_abundance),
              Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
    ) %>%
    as.data.frame(.) %>%
    mutate(., Level = "order") %>%
    filter(str_count(Lineage_Order, ";") >= 4)
  colnames(v_Lineage_Order_depth)[2] = "Lineage"

  v_Lineage_Class_depth = v_genome_tax_depth %>%
    group_by(across(all_of(c("Class","Lineage_Class")))) %>%
    summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
              Taxa_SumCoverage = sum(coverm_abundance),
              Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
    ) %>%
    as.data.frame(.) %>%
    mutate(., Level = "class") %>%
    filter(str_count(Lineage_Class, ";") >= 3)
  colnames(v_Lineage_Class_depth)[2] = "Lineage"

  v_Lineage_Phylum_depth = v_genome_tax_depth %>%
    group_by(across(all_of(c("Phylum","Lineage_Phylum")))) %>%
    summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
              Taxa_SumCoverage = sum(coverm_abundance),
              Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
    ) %>%
    as.data.frame(.) %>%
    mutate(., Level = "phylum") %>%
    filter(str_count(Lineage_Phylum, ";") >= 2)
  colnames(v_Lineage_Phylum_depth)[2] = "Lineage"

  v_Lineage_Domain_depth = v_genome_tax_depth %>%
    group_by(across(all_of(c("Domain","Lineage_Domain")))) %>%
    summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
              Taxa_SumCoverage = sum(coverm_abundance),
              Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
    ) %>%
    as.data.frame(.) %>%
    mutate(., Level = "superkingdom") %>%
    filter(str_count(Lineage_Domain, ";") >= 1)
  colnames(v_Lineage_Domain_depth)[2] = "Lineage"

  v_Lineage_Root_depth = v_genome_tax_depth %>%
    group_by(across(all_of(c("taxid_Root","Lineage_Root")))) %>%
    summarise(Taxa_ContigNum = sum(Taxa_ContigNum),
              Taxa_SumCoverage = sum(coverm_abundance),
              Taxa_RelativeAbundance = sum(relative_abundance), .groups = "drop"
    ) %>%
    as.data.frame(.) %>%
    mutate(., Level = "root")
  colnames(v_Lineage_Root_depth)[2] = "Lineage"

  colnames(v_Lineage_Root_depth)[1] = "Taxonomy"
  colnames(v_Lineage_Domain_depth)[1] = "Taxonomy"
  colnames(v_Lineage_Phylum_depth)[1] = "Taxonomy"
  colnames(v_Lineage_Class_depth)[1] = "Taxonomy"
  colnames(v_Lineage_Order_depth)[1] = "Taxonomy"
  colnames(v_Lineage_Family_depth)[1] = "Taxonomy"
  colnames(v_Lineage_Genus_depth)[1] = "Taxonomy"
  colnames(v_Lineage_Species_depth)[1] = "Taxonomy"
  colnames(v_Lineage_Strain_depth)[1] = "Taxonomy"

  v_Lineage_Root_depth = v_Lineage_Root_depth %>%
    {mutate(., Taxonomy = gsub("0","unclassified",.$Taxonomy))} %>%
    {mutate(., Taxonomy = gsub("1","root",.$Taxonomy))}

  All_Taxa_Lineage_depth = bind_rows(
    v_Lineage_Root_depth,
    v_Lineage_Domain_depth,  v_Lineage_Phylum_depth,
    v_Lineage_Class_depth,   v_Lineage_Order_depth,
    v_Lineage_Family_depth,  v_Lineage_Genus_depth,
    v_Lineage_Species_depth, v_Lineage_Strain_depth) %>%
    as.data.table()
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
    "@__program__:metalink pmr",
    "@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE"
  )
  #out_CAMI_format <- "output_cami_profile.txt"
  writeLines(header_lines, out_CAMI_format)
  write.table(All_Taxa_Lineage_depth, out_CAMI_format,sep = "\t",row.names = F,col.names = F,append = T,quote = F)
}
#############################################
############ function abundance ###############
#############################################

v_ORF_KO_depth[, scaleOrfBases := ORFLen*Scaled_Avg_fold]
v_ORF_KO_depth <- v_ORF_KO_depth[ORFLen >= as.numeric(minORFLen)]

ORF_depth_sum = v_ORF_KO_depth[, .(sumOrfBases = sum(scaleOrfBases),
                                   sumOrfDepth = sum(Scaled_Avg_fold)), keyby = sample]

ORF_eggnog_depth = v_ORF_KO_depth
ORF_eggnog_depth[ORF_depth_sum, on = "sample",
                 `:=`(sumOrfBases = i.sumOrfBases, sumOrfDepth = i.sumOrfDepth)]
ORF_eggnog_depth[, ORFRA := 1e6 * Scaled_Avg_fold/sumOrfDepth]
ORF_eggnog_depth[, ORFRA_Bases := 1e6 * scaleOrfBases/sumOrfBases]
ORF_eggnog_depth[, Source := "vMAG"]
ORF_eggnog_depth[, NewTax := gsub("unclassified", "Unclassified", NewTax)]
ORF_eggnog_depth[, NewTax := gsub("New", "Unclassified", NewTax)]

rm(v_ORF_KO_depth, ORF_depth_sum)
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

virus_seqid_tax = as.data.table(virus_seqid_tax)
virus_seqid_tax[, NewTax := gsub("unclassified", "Unclassified", NewTax)]
virus_seqid_tax[, NewTax := gsub("New", "Unclassified", NewTax)]

MAG_Tax_summary <- summary_with_all(virus_seqid_tax,
                                    c("Domain", "NewTax"),
                                    "MAG_Num", "Total_MAG", "Domain")

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Output the summary table of vMAGs taxonomy to:',summary_Tax,'\n')
#summary = "{sample}.MAG_Tax_summary.tsv"
MAG_Tax_summary %>%
  {fwrite(., summary_Tax, row.names = FALSE, col.names = T,
          sep = '\t', quote = FALSE)}

mode <- "pmr_v"
if (!is.null(out_mode_flag)) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Exporting execution mode:',mode,'to',out_mode_flag,'\n')
  dir.create(dirname(out_mode_flag), showWarnings = FALSE, recursive = TRUE)
  writeLines(mode, out_mode_flag)
}
