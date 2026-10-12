###
#' @Date: 2022-11-02 23:40:20
#' @LastEditors: LiuyangLi
#' @LastEditTime: 2024-06-18 21:18:51
#' @FilePath: /liliuyang/R/Rscript.genome/combine_metadata_MAGs.cmd.R
#' @Description: merge protine_stats

rm(list = ls(all.names = TRUE))
packages=c("dplyr","optparse","tidyverse")
ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE, repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    invisible(sapply(pkg, require, character.only = TRUE, quietly = TRUE))
    #cat(paste0("Loading packages: ", pkg, "\n"))
}
suppressMessages(ipak(packages))

option_list = list(
  make_option(c("--metadata"), type="character", dest="metadata", default=NULL,
              help="metadata table, colnames may include abundance, Domain,Phylum, Class,
                Order, Family,Genus, Species, e.g., mMAG_summary_info.tsv"),
  make_option(c("--drep_99_Cdb"), type="character", dest="drep_99_Cdb", default=NULL,
              help="drep_99_Cdb table, colnames would be transformed to primary_cluster, secondary_cluster, e.g., Cdb.csv"),
  make_option(c("--drep_95_Cdb"), type="character", dest="drep_95_Cdb", default=NULL,
              help="drep_95_Cdb table, colnames would be transformed to primary_cluster, secondary_cluster, e.g., Cdb.csv"),
  make_option(c("--drep_95_Wdb"), type="character", dest="drep_95_Wdb", default=NULL,
              help="drep_95 table, colnames would be transformed to primary_cluster, secondary_cluster, e.g., drep_95_Wdb.csv"),
  make_option(c("--drep_99_Wdb"), type="character", dest="drep_99_Wdb", default=NULL,
              help="Wdb.csv from the dRep 99 run (columns: genome,cluster,score; winners only). Used to pick the dRep_99 cluster representative. If missing/invalid, falls back to lexicographically largest ID"),
  make_option(c("-o", "--Output"), type="character", dest="Output", default=NULL,
              help="output Dir and prefix, e.g., meer.gtdbtk.tsv"),
  make_option(c("--verbose"), action="store_true", dest="verbose", default=FALSE,
              help="print messages")
)
opt = parse_args(OptionParser(option_list=option_list))
metadata = opt$metadata
drep_99_Cdb = opt$drep_99_Cdb
drep_95_Cdb = opt$drep_95_Cdb
drep_95_Wdb = opt$drep_95_Wdb
drep_99_Wdb = opt$drep_99_Wdb
Output = opt$Output
verbose = opt$verbose
##########################      Data input      ##########################

########### metadata ########
in_metadata = read.csv(metadata, sep = '\t', quote='', header = T, check.names=FALSE)
colnames(in_metadata)[1:3]=c('ID','dRep_99_representative_by_sample','dRep_95_representative_by_sample')
in_metadata = in_metadata %>%
  select(-c('dRep_99_representative_by_sample','dRep_95_representative_by_sample')) %>%
  unique()

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

########### merge ########
Ref_metadata = dRep_representative %>%
  {merge(., in_metadata, by = 'ID', all.y = T)}

Ref_metadata[is.na(Ref_metadata$dRep_99_representative),"dRep_99_representative"] = "CheckM_filter"
Ref_metadata[is.na(Ref_metadata$dRep_95_representative),"dRep_95_representative"] = "CheckM_filter"

Ref_metadata %>%
  {write.table(., Output,
               sep = "\t", row.names = F,  quote = FALSE,  fileEncoding="UTF-8")}
