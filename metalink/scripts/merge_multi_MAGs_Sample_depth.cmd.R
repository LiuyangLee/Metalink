###
#' @Date: 2023-08-18 20:24:02
#' @LastEditors: LiuyangLi
#' @LastEditTime: 2024-06-20 14:15:07
#' @FilePath: /liliuyang/R/Rscript.MASH/202308-v7/M0S2.Merge_multi_Sample_KO_depth.cmd.R
#' @Description:
###

suppressMessages(library(optparse))
rm(list = ls(all.names = TRUE))
packages=c("data.table")
ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE, repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    invisible(sapply(pkg, require, character.only = TRUE, quietly = TRUE))
    #cat(paste0("Loading packages: ", pkg, "\n"))
}
suppressMessages(ipak(packages))
option_list = list(
  make_option(c("--user_Sample_in_DIR", "-u"), type="character", default=NULL,
              help="user_Sample directory", metavar="character"),
  make_option(c("--out_file", "-o"), type="character", default=NULL,
              help="output file", metavar="character"),
  make_option(c("--tax_table"), type="character", default=NULL,
              help="optional comma-separated taxonomy table(s) (e.g. 02_quality_gate filter.tsv); only the seven rank columns are left-joined onto Genome, with vMAG-side superkingdom..species renamed to Domain..Species", metavar="character"),
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="print messages")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
user_Sample_in_DIR <- opt$user_Sample_in_DIR
out_file <- opt$out_file
tax_table <- opt$tax_table
verbose <- opt$verbose


count = 1
user_file <- dir(user_Sample_in_DIR)
if (length(user_file) == 0) {
  stop(paste0(
    "ERROR: no per-sample MAG abundance tables found in ", user_Sample_in_DIR, "\n",
    "This means every sample's MAG coverm table was missing or empty, i.e. no MAG had any mapped reads in any sample.\n",
    "Check upstream: (1) did any MAG pass the project quality gate (see 02_quality_gate), and ",
    "(2) did reads actually map to the MAG reference (see 03_mapping / 03_coverm_MAGs logs)?"
  ))
}
for (file in user_file) {
  print(file)
  user.tmp <- data.table(fread(paste(user_Sample_in_DIR, file, sep = '/'),
                               quote='',sep = '\t',header = T, check.names = FALSE))
  setnames(user.tmp, 1, 'Genome', skip_absent =T)
  n_na <- sum(is.na(user.tmp))
  if (n_na > 0) {
    cat(paste0("WARNING: ", file, " contains ", n_na, " NA value(s); replaced with 0 in the merged output, please check the upstream data\n"))
  }
  if (count == 1){
    result = user.tmp
  } else {
    result = merge.data.table(result, user.tmp, all = T, by = 'Genome')
  }
  count = count + 1
}

result[is.na(result)] <- 0

if (!is.null(tax_table) && tax_table != "") {
  tax_files <- strsplit(tax_table, ",", fixed = TRUE)[[1]]
  rank_from <- c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
  rank_to   <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  sample_cols <- setdiff(colnames(result), "Genome")
  # normalize rank columns with the same rename pipeline used by the combine
  # scripts (idempotent on already-normalized mMAG tables); sources RLib.r for
  # rename_GTDB_taxonomy_dataframe / rename_strain_dataframe
  script_path <- sub("--file=", "", commandArgs()[grepl("--file=", commandArgs())])
  source(file.path(dirname(script_path), "RLib.r"))
  normalize_tax <- function(tax, keep_cols) {
    if (!("Strain" %in% colnames(tax))) tax[, Strain := ""]
    tax <- tax %>%
      as.data.frame(.) %>%
      rename_GTDB_taxonomy_dataframe(.) %>%
      rename_strain_dataframe(.) %>%
      as.data.table(.)
    tax[, NewTax := gsub("New", "Unclassified", NewTax)]
    tax[, ..keep_cols]
  }
  if (length(tax_files) == 1) {
    # single-type matrix (pmr --mags m or v): self-contained species-level
    # summary -- carry ALL columns of the filter table; when ranks use CAT/BAT
    # names (superkingdom..species), additionally prepend canonical
    # Domain..Species with trailing confidence scores stripped
    dt <- data.table(fread(tax_files, quote = '', sep = '\t', header = TRUE, check.names = FALSE))
    setnames(dt, 1, 'Genome', skip_absent = TRUE)
    for (i in seq_along(rank_from))
      if (rank_from[i] %in% colnames(dt) && !(rank_to[i] %in% colnames(dt)))
        dt[, (rank_to[i]) := sub(":.*$", "", dt[[rank_from[i]]])]
    tax <- normalize_tax(dt, colnames(dt))
    front <- intersect(rank_to, colnames(tax))
    tax_cols <- c(front, setdiff(colnames(tax), c("Genome", front)))
  } else {
    # combined mMAG+vMAG matrix (pmr --mags vm): keep only the seven rank
    # columns per table so both sides align on the same columns (no sparse
    # union); CAT/BAT rank names are renamed and scores stripped
    tax_list <- lapply(tax_files, function(f) {
      dt <- data.table(fread(f, quote = '', sep = '\t', header = TRUE, check.names = FALSE))
      setnames(dt, 1, 'Genome', skip_absent = TRUE)
      for (i in seq_along(rank_from))
        if (rank_from[i] %in% colnames(dt) && !(rank_to[i] %in% colnames(dt)))
          setnames(dt, rank_from[i], rank_to[i])
      keep <- c("Genome", intersect(rank_to, colnames(dt)))
      dt <- dt[, ..keep]
      for (cl in keep[-1]) if (is.character(dt[[cl]])) set(dt, NULL, cl, sub(":.*$", "", dt[[cl]]))
      dt
    })
    tax_raw <- rbindlist(tax_list, use.names = TRUE, fill = TRUE)
    tax <- normalize_tax(tax_raw, c("Genome", rank_to))
    tax_cols <- intersect(rank_to, colnames(tax))
  }
  result <- merge.data.table(result, tax, by = 'Genome', all.x = TRUE, sort = FALSE)
  setcolorder(result, c("Genome", tax_cols, sample_cols))
  # group rows by taxonomy (Domain -> Species); within a group, natural-sort
  # genomes by the trailing number (bin.2 before bin.10), names without a
  # trailing number fall back to plain alphabetical order at the end
  result[, .natkey := suppressWarnings(as.integer(sub(".*?(\\d+)$", "\\1", Genome)))]
  setorderv(result, c(intersect(rank_to, colnames(result)), ".natkey", "Genome"), na.last = TRUE)
  result[, .natkey := NULL]
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Joined',length(tax_cols),'taxonomy/annotation columns from:',tax_table,'\n')
} else {
  result[, .natkey := suppressWarnings(as.integer(sub(".*?(\\d+)$", "\\1", Genome)))]
  setorderv(result, c(".natkey", "Genome"), na.last = TRUE)
  result[, .natkey := NULL]
}

fwrite(result, out_file, row.names = FALSE, col.names = T,
            sep = '\t', quote = FALSE)

print("merge_multi_MAGs_Sample_depth.cmd.R finished!")
