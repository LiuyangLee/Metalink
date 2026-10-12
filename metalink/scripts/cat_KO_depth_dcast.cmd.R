#!/usr/bin/env Rscript
suppressMessages(library(data.table))

suppressMessages(library(optparse))


# ------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------
option_list <- list(
  make_option("--sample_tsv_list", type = "character", default = NULL,
              help = "text file that lists full paths to per-sample tsv files, one per line"),
  make_option("--out_long", type = "character", default = "Merged.long.tsv",
              help = "output long-format table (before casting)"),
  make_option("--out_wide", type = "character", default = "Merged.wide.tsv",
              help = "output wide-format table; when multiple value_var are supplied this path is used as prefix"),
  make_option("--threads", type = "integer", default = 1,
              help = "number of CPU threads for data.table"),
  make_option("--id_var", type = "character", default = "KO",
              help = "column name(s) for the left-hand side of the casting formula"),
  make_option("--value_var", type = "character", default = "KO_RA",
              help = "column name(s) to be cast into wide format; comma-separated for several columns")
)

args <- parse_args(OptionParser(option_list = option_list))
setDTthreads(args$threads)

# ------------------------------------------------------------------
# Read list of input files
# ------------------------------------------------------------------
files <- readLines(args$sample_tsv_list, encoding = "UTF-8")
if (length(files) == 0) stop("ERROR: sample list file is empty.")

# ------------------------------------------------------------------
# Stack all per-sample tables and add a 'sample' column
# ------------------------------------------------------------------
dt_long <- rbindlist(
  lapply(files, function(f) {
    sample_id <- sub("(.*)\\..*\\..*$", "\\1", basename(f))
    fread(f)[, sample := sample_id]
  }),
  use.names = TRUE, fill = TRUE
)

# ------------------------------------------------------------------
# Save the complete long table (useful for debugging)
# ------------------------------------------------------------------
fwrite(dt_long, args$out_long, sep = "\t", quote = FALSE)

# ------------------------------------------------------------------
# Parse the value columns that must be cast
# ------------------------------------------------------------------
value_vars <- unlist(strsplit(args$value_var, ","))
missing    <- setdiff(value_vars, names(dt_long))
if (length(missing))
  stop("ERROR: value_var column(s) not found: ", paste(missing, collapse = ", "))

n_na <- sum(is.na(as.matrix(dt_long[, ..value_vars])))
if (n_na > 0) {
  cat(sprintf("WARNING: input tables contain %d NA value(s) in value column(s) (%s); replaced with 0, please check the upstream data\n",
              n_na, paste(value_vars, collapse = ", ")))
  for (v in value_vars) set(dt_long, which(is.na(dt_long[[v]])), v, 0)
}

# ------------------------------------------------------------------
# Build left-hand side of the formula
# ------------------------------------------------------------------
lhs_cols <- args$id_var
if (args$id_var == "Taxonomy") {
  if (!"Level" %in% names(dt_long))
    stop("Level column is required when id_var = 'Taxonomy'")
  lhs_cols <- c(args$id_var, "Level")
}

# ------------------------------------------------------------------
# Cast each value column separately and write the corresponding file
# ------------------------------------------------------------------
for (val in value_vars) {
  dt_wide <- dcast(dt_long,
                   as.formula(sprintf("%s ~ sample", paste(lhs_cols, collapse = " + "))),
                   value.var   = val,
                   fun.aggregate = sum,
                   fill = 0)
  # Output file name: replace the last ".wide.tsv" with ".{val}.wide.tsv"
  out_file <- sub("\\.wide\\.tsv$", sprintf(".%s.wide.tsv", val), args$out_wide)
  fwrite(dt_wide, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
}
