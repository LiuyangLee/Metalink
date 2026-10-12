###
#' @Date: 2022-11-02 23:40:20
#' @LastEditors: LiuyangLi
#' @LastEditTime: 2024-07-20 13:04:22
#' @FilePath: /liliuyang/R/Rscript.genome/combine_metadata_MAGs.cmd.R
#' @Description: 新增了 protine_stats 的内容

#rm(list = ls(all.names = TRUE))
packages=c("dplyr","optparse","tidyverse","data.table")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg))
    install.packages(new.pkg, dependencies = TRUE, repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  invisible(sapply(pkg, require, character.only = TRUE, quietly = TRUE))
  #cat(paste0("Loading packages: ", pkg, "\n"))
}
suppressMessages(ipak(packages))

#source("C:\\Users\\Administrator\\Desktop\\Script\\Rscript\\R-20241204\\RLib\\R\\RLib.r")
#in_taxid_tax="../taxid/taxid_tax_scale.tsv.gz"
#taxonomy_classify_method = "BAT" #gtdbtk

script_path <- sub("--file=", "", commandArgs()[grepl("--file=", commandArgs())])
if (length(script_path) == 0 || script_path == "") {
  script_path <- sys.frame(1)$ofile
}
source(file.path(dirname(script_path), "RLib.r"))
in_taxid_tax=NULL

option_list <- list(
  make_option(c("--checkM2"), type="character", default=NULL,
              help="Taxonomy table, colnames would be transformed to Genome_id, Domain,Phylum, Class,
                Order, Family,Genus, Species, e.g., meer_drep99_all.summary.tsv"),
  make_option(c("--in_taxid_tax"), type="character", default=NULL,
              help="taxid_tax table, e.g., taxid_tax.tsv"),
  make_option(c("--taxonomy_classify_method"), type="character", default=NULL,
              help="gtdbtk or BAT"),
  make_option(c("--taxonomy"), type="character", default=NULL,
              help="Taxonomy table, colnames would be transformed to Genome_id, Domain,Phylum, Class,
                Order, Family,Genus, Species, e.g., meer_drep99_all.summary.tsv"),
  make_option(c("--checkM"), type="character", default=NULL,
              help="checkM table, colnames would be transformed to completeness, contamination, e.g., genomeInformation.csv"),
  make_option(c("--genome_stat"), type="character", default=NULL,
              help="Genome statis table, colnames would be transformed to Genome, GC,..., e.g., genome_stat.tsv"),
  make_option(c("--protein_stat"), type="character", default=NULL,
              help="Protein statis table, colnames would be transformed to Genome, GC,..., e.g., protein_stat.tsv"),
  make_option(c("--in_coverm_trimmed_mean"), type="character", default=NULL,
              help="sum_coverm_trimmed_mean.tsv table, without colnames"),
  make_option(c("--drep_99_Cdb"), type="character", default=NULL,
              help="drep_99_Cdb table, colnames would be transformed to primary_cluster, secondary_cluster, e.g., Cdb.csv"),
  make_option(c("--drep_95_Cdb"), type="character", default=NULL,
              help="drep_95_Cdb table, colnames would be transformed to primary_cluster, secondary_cluster, e.g., Cdb.csv"),
  make_option(c("--drep_95_Wdb"), type="character", dest="drep_95_Wdb", default=NULL,
              help="drep_95 table, colnames would be transformed to primary_cluster, secondary_cluster, e.g., drep_95_Wdb.csv"),
  make_option(c("--drep_99_Wdb"), type="character", dest="drep_99_Wdb", default=NULL,
              help="Wdb.csv from the dRep 99 run (columns: genome,cluster,score; winners only). Used to pick the dRep_99 cluster representative. If missing/invalid, falls back to lexicographically largest ID"),
  make_option(c("--out_genome"), type="character", default=NULL,
              help="output file for mMAG_summary_info, e.g., mMAG_summary_info.tsv"),
  make_option(c("--out_species"), type="character", default=NULL,
              help="output file for mMAG_species_summary_info, e.g., mMAG_species_summary_info.tsv"),
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="print messages")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

checkM2 <- opt$checkM2
in_taxid_tax <- opt$in_taxid_tax
taxonomy_classify_method <- opt$taxonomy_classify_method
taxonomy <- opt$taxonomy
checkM <- opt$checkM
genome_stat <- opt$genome_stat
protein_stat <- opt$protein_stat
in_coverm_trimmed_mean <- opt$in_coverm_trimmed_mean
drep_99_Cdb <- opt$drep_99_Cdb
drep_95_Cdb <- opt$drep_95_Cdb
drep_95_Wdb <- opt$drep_95_Wdb
drep_99_Wdb <- opt$drep_99_Wdb
out_genome <- opt$out_genome
out_species <- opt$out_species
verbose <- opt$verbose
##########################      Data input      ##########################

########### checkM2 ########
Ref_checkM2 = read.csv(checkM2, sep = '\t', quote='', header = T, check.names=FALSE) %>%
  {.[,1:5]}
colnames(Ref_checkM2)=c(
  'ID',	'Completeness2','Contamination2',
  'Completeness_Model_Used','Translation_Table_Used'
)

################################################################
####################     taxid_taxonomy      ######################
################################################################
cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract mapping relationship of taxid-taxonomy from:',in_taxid_tax,'\n')
taxid_taxonomy = fread(in_taxid_tax, sep = "\t", colClasses = "character", col.names = c('taxid','Taxonomy','Lineage'), header = FALSE, data.table = T)

########### taxonomy ########

cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'taxonomy_classify_method:',taxonomy_classify_method,'\n')
if (taxonomy_classify_method == "gtdbtk"){
  #taxonomy="../results/MetaG-MAGs-sample/02_taxonomy/gtdbtk.bin2classification.tsv"
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract MAGs taxononmy from:',taxonomy,'\n')
  Ref_gtdbtk = read.csv(taxonomy, sep = '\t', quote='', header = T, check.names=FALSE) %>%
    {.[,1:3]}
  colnames(Ref_gtdbtk)=c('ID', 'Taxonomy', 'closest_genome_reference')
  Ref_gtdbtk$Raw_Taxonomy = Ref_gtdbtk$Taxonomy
  seqid_taxid_detail = Ref_gtdbtk %>%
    gtdb_tax() %>%
    #mutate(.,  Species =  gsub(" ", "_", .$Species)) %>%
    mutate(.,  Strain =  paste0("t__", .$ID)) %>%
    {.[,colnames(.)%in%  c("ID", "Domain","Phylum", "Class",
                           "Order", "Family","Genus", "Species","Strain","closest_genome_reference", "Raw_Taxonomy"
    )]}
} else if (taxonomy_classify_method == "BAT"){
  #taxonomy="../results/MetaG-MAGs-sample/02_taxonomy/BAT.bin2classification.tsv"
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Extract MAGs taxononmy from:',taxonomy,'\n')
  seqid_taxid_detail = fread(taxonomy, sep = "\t", header = T, data.table = T, fill = T)
  colnames(seqid_taxid_detail) = c("ID", "classification","reason", "lineage", "lineage_scores","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  seqid_taxid_detail[, lineage := gsub("\\*", "", lineage)]
  seqid_taxid_detail[, taxid := gsub(".*;", "", lineage)]
  seqid_taxid_detail[is.na(taxid)| taxid == "", taxid := "0"]
  seqid_taxid_detail[, Raw_Taxonomy := paste(Domain, Phylum, Class, Order, Family, Genus, Species, sep = ";")]
  seqid_taxid_detail[, Raw_Taxonomy := gsub(":\\s*[0-9.]+", "", Raw_Taxonomy)]
  seqid_taxid_detail = seqid_taxid_detail %>%
    merge.data.table(., taxid_taxonomy, by = 'taxid', all.x = T)
  seqid_taxid_detail[, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain") := tstrsplit(Taxonomy, ";")]
  seqid_taxid_detail = seqid_taxid_detail %>%
    select(., c("ID", "taxid", "Lineage","Raw_Taxonomy","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species","Strain"))
  seqid_taxid_detail$ID = seqid_taxid_detail$ID %>% {gsub("\\.fna$", "", .)}
  seqid_taxid_detail$ID = seqid_taxid_detail$ID %>% {gsub("\\.faa$", "", .)}
  seqid_taxid_detail$ID = seqid_taxid_detail$ID %>% {gsub("\\.fa$", "", .)}
}

Ref_taxonomy = seqid_taxid_detail %>%
  as.data.frame() %>%
  rename_GTDB_taxonomy_dataframe(.) %>%
  {mutate(., NewTax = gsub("New", "Unclassified", .$NewTax))}
########### checkM from dRep_99 ########
Ref_checkM = read.csv(checkM, sep = ',', quote='', header = T, check.names=FALSE) %>%
  {.[,c(1:4)]}
colnames(Ref_checkM)=c(
  'ID',	'Completeness1','Contamination1','strain_heterogeneity'
)
Ref_checkM$ID = Ref_checkM$ID %>% {gsub("\\.fna$", "", .)}
Ref_checkM$ID = Ref_checkM$ID %>% {gsub("\\.faa$", "", .)}
Ref_checkM$ID = Ref_checkM$ID %>% {gsub("\\.fa$", "", .)}

########### genome_stat & protein_stat ########
Ref_genome_stat = read.csv(genome_stat, sep = '\t', quote='', header = T, check.names=FALSE)
Ref_genome_stat$ID = Ref_genome_stat[,1]

Ref_protein_stat = read.csv(protein_stat, sep = '\t', quote='', header = T, check.names=FALSE)
Ref_protein_stat$ID = Ref_protein_stat[,1]

########### dRep_99 ########
Ref_dRep_99 = read.csv(drep_99_Cdb, sep = ',', quote='', header = T, check.names=FALSE) %>%
  {.[,c(1,2,6)]}
colnames(Ref_dRep_99)=c(
  'ID',	'dRep_99_cluster','dRep_90_cluster'
)
Ref_dRep_99$ID = Ref_dRep_99$ID %>% {gsub("\\.fna$", "", .)}
Ref_dRep_99$ID = Ref_dRep_99$ID %>% {gsub("\\.faa$", "", .)}
Ref_dRep_99$ID = Ref_dRep_99$ID %>% {gsub("\\.fa$", "", .)}

########### dRep_95 ########
Ref_dRep_95 = read.csv(drep_95_Cdb, sep = ',', quote='', header = T, check.names=FALSE) %>%
  {.[,c(1,2)]}
colnames(Ref_dRep_95)=c(
  'ID',	'dRep_95_cluster'
)
Ref_dRep_95$ID = Ref_dRep_95$ID %>% {gsub("\\.fna$", "", .)}
Ref_dRep_95$ID = Ref_dRep_95$ID %>% {gsub("\\.faa$", "", .)}
Ref_dRep_95$ID = Ref_dRep_95$ID %>% {gsub("\\.fa$", "", .)}

########### drep_95_Wdb.csv ########
drep_95_Wdb_df = read.csv(drep_95_Wdb, sep = ',', quote='', header = T, check.names=FALSE) %>%
  {.[,c(1,3)]}
colnames(drep_95_Wdb_df)=c(
  'ID',	'score'
)
drep_95_Wdb_df$ID = drep_95_Wdb_df$ID %>% {gsub("\\.fna$", "", .)}
drep_95_Wdb_df$ID = drep_95_Wdb_df$ID %>% {gsub("\\.faa$", "", .)}
drep_95_Wdb_df$ID = drep_95_Wdb_df$ID %>% {gsub("\\.fa$", "", .)}

########### drep_99_Wdb.csv ########
drep_99_Wdb_df = NULL
if (!is.null(drep_99_Wdb) && file.exists(drep_99_Wdb) && file.size(drep_99_Wdb) > 0) {
  drep_99_Wdb_df = tryCatch(
    read.csv(drep_99_Wdb, sep = ',', quote='', header = T, check.names=FALSE),
    error = function(e) NULL
  )
  if (is.null(drep_99_Wdb_df) || ncol(drep_99_Wdb_df) < 3 || nrow(drep_99_Wdb_df) == 0) {
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING: drep_99_Wdb file',drep_99_Wdb,'is invalid; dRep_99 cluster representatives will fall back to lexicographically largest ID','\n')
    drep_99_Wdb_df = NULL
  } else {
    drep_99_Wdb_df = drep_99_Wdb_df %>%
      {.[,c(1,3)]}
    colnames(drep_99_Wdb_df)=c(
      'ID',	'score_99'
    )
    drep_99_Wdb_df$ID = drep_99_Wdb_df$ID %>% {gsub("\\.fna$", "", .)}
    drep_99_Wdb_df$ID = drep_99_Wdb_df$ID %>% {gsub("\\.faa$", "", .)}
    drep_99_Wdb_df$ID = drep_99_Wdb_df$ID %>% {gsub("\\.fa$", "", .)}
  }
} else {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING: drep_99_Wdb is not provided, missing or empty (',drep_99_Wdb,'); dRep_99 cluster representatives will fall back to lexicographically largest ID','\n')
}

########### get dRep_99_representative ########
dRep_representative = Ref_dRep_99 %>%
  {merge(., Ref_dRep_95, by = 'ID', all.x = T)} %>%
  {merge(., drep_95_Wdb_df, by = 'ID', all.x = T)} %>%
  {mutate(., dRep_99_representative = .$ID)}
if (!is.null(drep_99_Wdb_df)) {
  dRep_representative = dRep_representative %>%
    {merge(., drep_99_Wdb_df, by = 'ID', all.x = T)}
}

duplicate_values_99 <- unique(dRep_representative[duplicated(dRep_representative$dRep_99_cluster), 'dRep_99_cluster'])
dRep_representative_sub_non_NA_99 = dRep_representative %>%
  {.[.$dRep_99_cluster %in% duplicate_values_99, ]} %>%
  {.[!is.na(.$dRep_95_cluster),]}
for (variable in duplicate_values_99) {
  sub_rows_99 = dRep_representative_sub_non_NA_99[dRep_representative_sub_non_NA_99$dRep_99_cluster %in% variable, ]
  if (nrow(sub_rows_99) == 0) {
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING: dRep_99_cluster ',variable,' has no non-NA dRep_95_cluster; its dRep_95_cluster/dRep_99_representative are left unchanged, please check the dRep_95 inputs','\n')
    next
  }
  dRep_representative[dRep_representative$dRep_99_cluster %in% variable, "dRep_95_cluster"] =
    max(sub_rows_99$dRep_95_cluster)
  rep99 = max(sub_rows_99$ID)  # fallback: lexicographic max
  if (!is.null(drep_99_Wdb_df)) {
    winners_99 = sub_rows_99[!is.na(sub_rows_99$score_99), ]
    if (nrow(winners_99) > 0) {
      rep99 = winners_99$ID[which.max(winners_99$score_99)]
    } else {
      cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING: dRep_99_cluster ',variable,' has no winner in drep_99_Wdb; falling back to lexicographically largest ID','\n')
    }
  }
  dRep_representative[dRep_representative$dRep_99_cluster %in% variable, "dRep_99_representative"] =
    rep99
}
#dRep_representative$dRep_95_cluster = paste0(dRep_representative$dRep_90_cluster,"_",
#                                       sapply(as.character(dRep_representative$dRep_95_cluster),
#                                              function(x) unlist(strsplit(x, "_"))[2]))
#
########### get dRep_95_representative ########
dRep_representative = dRep_representative %>%
  {mutate(., dRep_95_representative = .$dRep_99_representative)}
duplicate_values_95 <- unique(dRep_representative[duplicated(dRep_representative$dRep_95_cluster), 'dRep_95_cluster'])
#duplicate_values_95 = duplicate_values_95[!is.na(duplicate_values_95)]
dRep_representative_sub_non_NA_95 = dRep_representative %>%
  {.[.$dRep_95_cluster %in% duplicate_values_95, ]} %>%
  {.[!is.na(.$score),]}
for (variable in duplicate_values_95) {
  sub_rows_95 = dRep_representative_sub_non_NA_95[dRep_representative_sub_non_NA_95$dRep_95_cluster %in% variable, ]
  if (nrow(sub_rows_95) == 0) {
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING: dRep_95_cluster ',variable,' has no non-NA score; its dRep_95_representative is left unchanged, please check the drep_95_Wdb input','\n')
    next
  }
  dRep_representative[dRep_representative$dRep_95_cluster %in% variable, "dRep_95_representative"] =
    max(sub_rows_95$dRep_99_representative)
}

dRep_representative = dRep_representative %>%
  {dplyr::select(., -any_of(c("dRep_99_cluster","dRep_90_cluster","dRep_95_cluster","score","score_99")))}

### coverm ###
Ref_coverm_trimmed_mean = read.csv(in_coverm_trimmed_mean, sep = '\t', quote='', header = T, check.names=FALSE)
colnames(Ref_coverm_trimmed_mean)=c('ID',	'coverm_abundance')
n_na_coverm <- sum(is.na(Ref_coverm_trimmed_mean$coverm_abundance))
if (n_na_coverm > 0) {
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'WARNING:',in_coverm_trimmed_mean,'contains',n_na_coverm,'NA value(s) in coverm_abundance; replaced with 0, please check the upstream data','\n')
  Ref_coverm_trimmed_mean$coverm_abundance[is.na(Ref_coverm_trimmed_mean$coverm_abundance)] <- 0
}
Ref_coverm_trimmed_mean = Ref_coverm_trimmed_mean %>%
  {.[!(.$ID %in% 'unmapped'), ]} %>%
  mutate(coverm_abundance = as.numeric(coverm_abundance)) %>%
  mutate(
    total_abundance = sum(coverm_abundance, na.rm = TRUE),
    relative_abundance = if_else(
      total_abundance == 0,
      0,
      coverm_abundance / total_abundance
    )
  ) %>%
  select(-total_abundance)

dRep_coverm_genome = dRep_representative %>%
  {merge(., Ref_coverm_trimmed_mean, by = 'ID', all.x = T)}
dRep_coverm_genome[is.na(dRep_coverm_genome$dRep_99_representative),"dRep_99_representative"] = "CheckM_filter"
dRep_coverm_genome[is.na(dRep_coverm_genome$dRep_95_representative),"dRep_95_representative"] = "CheckM_filter"

Ref_metadata_genome = dRep_coverm_genome %>%
  {merge(., Ref_taxonomy, by = 'ID', all.x = T)} %>%
  {merge(., Ref_checkM2, by = 'ID', all.x = T)} %>%
  {merge(., Ref_checkM, by = 'ID', all.x = T)} %>%
  {merge(., Ref_genome_stat, by = 'ID', all.x = T)} %>%
  {merge(., Ref_protein_stat, by = 'ID', all.x = T)} %>%
  {mutate(., `Coding density (%)` = 100 * (.$`Coding bases`)/(.$`Total base (bp)`))} %>%
  {mutate(., `Average gene length` = (.$`Coding bases`)/(.$`Protein count`))}

colnames(Ref_metadata_genome)[1] = 'Prokaryotes'
Ref_metadata_genome %>%
  {write.table(., out_genome,
               sep = "\t", row.names = F,  quote = FALSE,  fileEncoding="UTF-8")}


dRep_coverm_species = dRep_coverm_genome %>%
  group_by(dRep_95_representative) %>%
  summarise(
    sum_coverm_abundance = sum(coverm_abundance, na.rm = TRUE),
    sum_relative_abundance = sum(relative_abundance, na.rm = TRUE)
  ) %>%
  mutate(ID = dRep_95_representative)

Ref_metadata_species = dRep_coverm_species %>%
  {merge(., Ref_taxonomy, by = 'ID', all.x = T)} %>%
  {merge(., Ref_checkM2, by = 'ID', all.x = T)} %>%
  {merge(., Ref_checkM, by = 'ID', all.x = T)} %>%
  {merge(., Ref_genome_stat, by = 'ID', all.x = T)} %>%
  {merge(., Ref_protein_stat, by = 'ID', all.x = T)} %>%
  {mutate(., `Coding density (%)` = 100 * (.$`Coding bases`)/(.$`Total base (bp)`))} %>%
  {mutate(., `Average gene length` = (.$`Coding bases`)/(.$`Protein count`))}

colnames(Ref_metadata_species)[1] = 'Prokaryotes'
Ref_metadata_species %>%
  {write.table(., out_species,
               sep = "\t", row.names = F,  quote = FALSE,  fileEncoding="UTF-8")}
