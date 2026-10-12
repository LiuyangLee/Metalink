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
              help = "output long-format table (before aggregation)"),
  make_option("--out_wide", type = "character", default = "Merged.wide.tsv",
              help = "output wide-format table prefix"),
  make_option("--threads", type = "integer", default = 1,
              help = "number of CPU threads for data.table"),
  make_option("--id_var", type = "character", default = "KO,Phylum",
              help = "column name(s) for grouping (comma-separated); First = row, Second = column in wide format"),
  make_option("--value_var", type = "character", default = "ORF_Num,Taxa_KO_SumCoverage,Taxa_KO_RelativeAbundance",
              help = "column name(s) to aggregate and cast; comma-separated"),
  make_option("--min_abundance", type = "double", default = 0,
              help = "minimum mean relative abundance threshold on the per-million (GPM) scale of Taxa_KO_RelativeAbundance (default 0 = no filter; e.g., 1 = 1 GPM = 1/1e6)"),
  make_option("--min_prevalence", type = "double", default = 0,
              help = "minimum prevalence threshold (default 0 = no filter; e.g., 0.1 = 10%%)")
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
# Parse columns
# ------------------------------------------------------------------
id_vars <- unlist(strsplit(args$id_var, ","))
if (length(id_vars) != 2)
  stop("ERROR: id_var must contain exactly 2 columns (row_var ~ col_var), e.g., 'KO,Phylum'")

row_var <- id_vars[1]
col_var <- id_vars[2]
value_vars <- unlist(strsplit(args$value_var, ","))

missing_id <- setdiff(id_vars, names(dt_long))
if (length(missing_id))
  stop("ERROR: id_var column(s) not found: ", paste(missing_id, collapse = ", "))

missing_val <- setdiff(value_vars, names(dt_long))
if (length(missing_val))
  stop("ERROR: value_var column(s) not found: ", paste(missing_val, collapse = ", "))

n_na <- sum(is.na(as.matrix(dt_long[, ..value_vars])))
if (n_na > 0) {
  cat(sprintf("WARNING: input tables contain %d NA value(s) in value column(s) (%s); replaced with 0, please check the upstream data\n",
              n_na, paste(value_vars, collapse = ", ")))
  for (v in value_vars) set(dt_long, which(is.na(dt_long[[v]])), v, 0)
}

# Total number of samples
n_samples <- length(files)

# ------------------------------------------------------------------
# 1. Deduplicate within the same sample (sum merge duplicate rows)
# ------------------------------------------------------------------
dt_long[, (value_vars) := lapply(.SD, sum),
        by = c("sample", id_vars),
        .SDcols = value_vars]

dt_long <- unique(dt_long, by = c("sample", id_vars))

# Keep only needed columns
keep_cols <- c("sample", id_vars, value_vars)
dt_long <- dt_long[, ..keep_cols]

# ------------------------------------------------------------------
# 2. Calculate cross-sample mean: sum(non-NA) / total samples (missing treated as 0)
# ------------------------------------------------------------------
dt_mean <- dt_long[, lapply(.SD, function(x) sum(x, na.rm = TRUE) / n_samples),
                   by = id_vars,
                   .SDcols = value_vars]

fwrite(dt_mean, sub("\\.tsv$", ".tsv", args$out_long), sep = "\t", quote = FALSE)

# ------------------------------------------------------------------
# 3. Calculate prevalence (proportion of samples where the combination appears)
# ------------------------------------------------------------------
prevalence_dt <- dt_long[, .(prevalence = .N / n_samples), by = id_vars]
dt_mean <- merge(dt_mean, prevalence_dt, by = id_vars)

# ------------------------------------------------------------------
# 4. Optional: apply abundance and prevalence filtering
# ------------------------------------------------------------------
cat(sprintf("\n=== Pre-filter statistics ===\n"))
cat(sprintf("Total samples: %d\n", n_samples))
cat(sprintf("Total long-format rows: %d\n", nrow(dt_long)))
cat(sprintf("Unique combinations (KO;%s): %d\n", col_var, nrow(dt_mean)))

rel_abund_col <- "Taxa_KO_RelativeAbundance"
if (rel_abund_col %in% names(dt_mean)) {
  cat(sprintf("Mean relative abundance range: %.6f ~ %.6f\n",
              min(dt_mean[[rel_abund_col]], na.rm = TRUE),
              max(dt_mean[[rel_abund_col]], na.rm = TRUE)))
}
cat(sprintf("Prevalence range: %.2f%% ~ %.2f%%\n",
            min(dt_mean$prevalence, na.rm = TRUE) * 100,
            max(dt_mean$prevalence, na.rm = TRUE) * 100))

# Determine whether to apply filtering: enabled when min_abundance > 0 or min_prevalence > 0
do_filter <- args$min_abundance > 0 || args$min_prevalence > 0

if (do_filter) {
  cat(sprintf("\n=== Applying filter ===\n"))
  if (args$min_abundance > 0) {
    cat(sprintf("Minimum mean relative abundance: %.4f%% (%.4f)\n", args$min_abundance * 100, args$min_abundance))
  }
  if (args$min_prevalence > 0) {
    cat(sprintf("Minimum prevalence: %.1f%% (%.2f)\n", args$min_prevalence * 100, args$min_prevalence))
  }

  # Build filter condition
  keep_idx <- rep(TRUE, nrow(dt_mean))

  if (args$min_abundance > 0 && rel_abund_col %in% names(dt_mean)) {
    keep_idx <- keep_idx & (dt_mean[[rel_abund_col]] >= args$min_abundance)
  } else if (args$min_abundance > 0) {
    cat("WARNING: Taxa_KO_RelativeAbundance column not found, skipping abundance filter\n")
  }

  if (args$min_prevalence > 0) {
    keep_idx <- keep_idx & (dt_mean$prevalence >= args$min_prevalence)
  }

  cat(sprintf("Combinations passing filter: %d / %d (%.1f%%)\n",
              sum(keep_idx), nrow(dt_mean), sum(keep_idx) / nrow(dt_mean) * 100))

  # Get filtered id_vars combinations
  filtered_ids <- dt_mean[keep_idx, ..id_vars]

  # Filter dt_mean and dt_long
  dt_mean <- dt_mean[keep_idx]
  dt_long <- dt_long[filtered_ids, on = id_vars]

  cat(sprintf("Long-format rows after filtering: %d\n", nrow(dt_long)))
  cat(sprintf("Unique combinations after filtering: %d\n", nrow(dt_mean)))
} else {
  cat("\n=== No filter applied (min_abundance=0, min_prevalence=0) ===\n")
  cat("Tip: If the output matrix is too large, set filter thresholds, e.g.:\n")
  cat("  --min_abundance 0.0001 --min_prevalence 0.1\n")
}

# ------------------------------------------------------------------
# 5. Cast: KO ~ Phylum (mean value matrix)
# ------------------------------------------------------------------
for (val in value_vars) {
  dt_wide <- dcast(dt_mean,
                   as.formula(sprintf("%s ~ %s", row_var, col_var)),
                   value.var     = val,
                   fun.aggregate = sum,
                   fill          = 0)

  out_file <- sub("\\.wide\\.tsv$", sprintf(".%s.mean.wide.tsv", val), args$out_wide)
  fwrite(dt_wide, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Written (mean): %s (%d rows x %d cols)\n", out_file, nrow(dt_wide), ncol(dt_wide)))
}

# ------------------------------------------------------------------
# 6. Per-sample wide matrix (KO;Phylum ~ sample)
# ------------------------------------------------------------------
#fwrite(dt_long, sub("\\.mean.long.tsv$", ".per_sample.filter.long.tsv", args$out_long), sep = "\t", quote = FALSE)
dt_long[, KO_Phylum := paste(get(row_var), get(col_var), sep = ";")]

for (val in value_vars) {
  n_long_rows <- nrow(dt_long)
  n_unique_ko_phylum <- length(unique(dt_long$KO_Phylum))
  n_unique_samples <- length(unique(dt_long$sample))
  n_wide_rows <- n_unique_ko_phylum
  n_wide_cols <- n_unique_samples + 1

  cat(sprintf("\n=== Generating per-sample wide matrix: %s ===\n", val))
  cat(sprintf("Long-format rows: %d\n", n_long_rows))
  cat(sprintf("Will cast to wide matrix: %d rows x %d cols (KO;Phylum ~ sample)\n", n_wide_rows, n_wide_cols))

  est_memory_mb <- (n_wide_rows * n_wide_cols * 8) / (1024 * 1024)
  cat(sprintf("Estimated memory usage: ~%.1f MB\n", est_memory_mb))

  if (n_wide_rows > 10000 || n_wide_cols > 500) {
    cat("WARNING: Matrix size is large!\n")
    cat("  Suggestion: set stricter filter thresholds\n")
    cat(sprintf("  e.g.: --min_abundance 1 (sum: 1,000,000) and/or --min_prevalence 0.3 (sum: 1)\n"))
  }

  dt_wide_per_sample <- dcast(dt_long,
                              KO_Phylum ~ sample,
                              value.var = val,
                              fun.aggregate = sum,
                              fill = 0)

  out_file <- sub("\\.wide\\.tsv$", sprintf(".%s.per_sample.filter.wide.tsv", val), args$out_wide)

  fwrite(dt_wide_per_sample, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Written (per-sample): %s (%d rows x %d cols)\n", out_file, nrow(dt_wide_per_sample), ncol(dt_wide_per_sample)))
}

cat("\nDone.\n")
