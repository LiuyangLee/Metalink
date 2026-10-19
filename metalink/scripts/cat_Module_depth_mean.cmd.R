#!/usr/bin/env Rscript
# ==============================================================================
# cat_Module_depth_mean.cmd.R
# Project-level Module abundance merge (analogous to cat_KO_depth_mean.cmd.R)
# Input: per-sample mc_abundance.tsv files (wide format with Module metadata)
# Output: merged long table (cross-sample mean), mean wide matrix, per-sample wide matrix
# ==============================================================================
suppressMessages(library(data.table))

suppressMessages(library(optparse))


# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------
option_list <- list(
  make_option("--sample_tsv_list", type = "character", default = NULL,
              help = "text file that lists full paths to per-sample mc_abundance.tsv files, one per line"),
  make_option("--out_long", type = "character", default = "Merged.long.tsv",
              help = "output long-format mean table (cross-sample mean)"),
  make_option("--out_wide", type = "character", default = "Merged.wide.tsv",
              help = "output wide-format table prefix"),
  make_option("--threads", type = "integer", default = 1,
              help = "number of CPU threads for data.table"),
  make_option("--id_var", type = "character", default = "Module_Entry,Phylum",
              help = "column name(s) for grouping (comma-separated); First = row, Second = column in wide format"),
  make_option("--value_var", type = "character", default = "Taxa_Module_RelativeAbundance",
              help = "column name(s) to aggregate and cast; comma-separated"),
  make_option("--min_abundance", type = "double", default = 0,
              help = "minimum mean relative abundance threshold (default 0 = no filter)"),
  make_option("--min_prevalence", type = "double", default = 0,
              help = "minimum prevalence threshold (default 0 = no filter)")
)

args <- parse_args(OptionParser(option_list = option_list))
setDTthreads(args$threads)

# ------------------------------------------------------------------------------
# Module metadata columns (fixed, always preserved)
# ------------------------------------------------------------------------------
module_meta_cols <- c("Level_2", "Level_3", "Module_Name")

# ------------------------------------------------------------------------------
# Extract sample name from path: .../{sample}/{rank}/mc_abundance.tsv
# ------------------------------------------------------------------------------
extract_sample_id <- function(file_path) {
  fp <- normalizePath(file_path, mustWork = FALSE)
  parts <- strsplit(fp, "/")[[1]]
  mc_idx <- which(parts == "mc_abundance.tsv")
  if (length(mc_idx) == 1 && mc_idx > 2) {
    return(parts[mc_idx - 2])
  }
  return(sub("(.*)\\..*$", "\\1", basename(file_path)))
}

# ------------------------------------------------------------------------------
# Read list of input files
# ------------------------------------------------------------------------------
files <- readLines(args$sample_tsv_list, encoding = "UTF-8")
if (length(files) == 0) {
  # Not an error: every per-sample input was missing/empty upstream (e.g.
  # nothing classified, or no taxon had enough KOs). Write header-only
  # outputs so the workflow continues; downstream steps skip empty tables.
  id_vars_e <- unlist(strsplit(args$id_var, ","))
  value_vars_e <- unlist(strsplit(args$value_var, ","))
  cat("NOTE: sample list file is empty - no non-empty per-sample tables were\n")
  cat("available (all samples were skipped upstream, e.g. nothing classified or\n")
  cat("no taxon had enough KOs). Writing empty outputs; downstream steps skip them.\n")
  dir.create(dirname(args$out_long), showWarnings = FALSE, recursive = TRUE)
  writeLines(paste(c(id_vars_e, value_vars_e), collapse = "\t"), args$out_long)
  for (val in value_vars_e) {
    writeLines(id_vars_e[1], sub("\\.wide\\.tsv$", sprintf(".%s.mean.wide.tsv", val), args$out_wide))
    writeLines(id_vars_e[1], sub("\\.wide\\.tsv$", sprintf(".%s.per_sample.filter.wide.tsv", val), args$out_wide))
  }
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# Parse columns
# ------------------------------------------------------------------------------
id_vars <- unlist(strsplit(args$id_var, ","))
if (length(id_vars) != 2)
  stop("ERROR: id_var must contain exactly 2 columns (row_var ~ col_var), e.g., 'Module_Entry,Phylum'")

row_var <- id_vars[1]
col_var <- id_vars[2]
value_vars <- unlist(strsplit(args$value_var, ","))

# ------------------------------------------------------------------------------
# Stack all per-sample tables: convert wide to long
# ------------------------------------------------------------------------------
all_meta_list <- list()

dt_long <- rbindlist(
  lapply(files, function(f) {
    sample_id <- extract_sample_id(f)
    cat(sprintf("Reading: %s (sample: %s)\n", f, sample_id))

    # Read wide format
    dt_wide <- fread(f, showProgress = FALSE)

    # Validate row_var exists
    if (!(row_var %in% names(dt_wide)))
      stop(sprintf("ERROR: row_var '%s' not found in %s", row_var, f))

    # Identify metadata and species columns
    meta_cols_in_file <- intersect(names(dt_wide), c(row_var, module_meta_cols))
    species_cols <- setdiff(names(dt_wide), meta_cols_in_file)

    if (length(species_cols) == 0)
      stop(sprintf("ERROR: No species columns found in %s", f))

    # Extract and store Module metadata
    if (row_var %in% names(dt_wide)) {
      meta_unique <- unique(dt_wide[, ..meta_cols_in_file])
      all_meta_list[[length(all_meta_list) + 1]] <<- meta_unique
    }

    # Melt to long format
    dt_long_sample <- melt(dt_wide,
                          id.vars = meta_cols_in_file,
                          measure.vars = species_cols,
                          variable.name = "Species",
                          value.name = value_vars,
                          na.rm = TRUE)

    # Ensure value columns are numeric
    for (v in value_vars) {
      dt_long_sample[, (v) := as.numeric(get(v))]
    }

    # Remove NA and zeros
    for (v in value_vars) {
      dt_long_sample <- dt_long_sample[!is.na(get(v)) & get(v) > 0]
    }

    if (nrow(dt_long_sample) == 0) {
      cat(sprintf("  WARNING: %s has no non-zero values after melt\n", sample_id))
      return(NULL)
    }

    # Extract rank from Species column (taxonomy path)
    dt_long_sample[, (col_var) := sapply(as.character(Species), function(s) {
      parts <- strsplit(s, ";")[[1]]
      rank_prefix <- switch(col_var,
        "Domain"  = "d__", "Phylum"  = "p__", "Class"   = "c__",
        "Order"   = "o__", "Family"  = "f__", "Genus"   = "g__", "Species" = "s__",
        stop("Unknown rank: ", col_var))
      for (part in parts) {
        if (startsWith(part, rank_prefix)) {
          return(substr(part, nchar(rank_prefix) + 1, nchar(part)))
        }
      }
      return(NA_character_)
    })]

    # Add sample ID
    dt_long_sample[, sample := sample_id]

    # Remove rows with NA rank
    dt_long_sample <- dt_long_sample[!is.na(get(col_var))]

    # Keep needed columns
    keep_cols <- c("sample", row_var, col_var, value_vars)
    keep_cols <- c(meta_cols_in_file, keep_cols)
    keep_cols <- unique(keep_cols)
    dt_long_sample <- dt_long_sample[, ..keep_cols]

    cat(sprintf("  -> %d rows after processing\n", nrow(dt_long_sample)))
    return(dt_long_sample)
  }),
  use.names = TRUE, fill = TRUE
)

if (nrow(dt_long) == 0) {
  # All per-sample tables were read but contained no usable (non-zero) values:
  # a valid empty result, not an error. Write header-only outputs and exit.
  cat("NOTE: no usable (non-zero) values in any per-sample table; writing\n")
  cat("empty outputs so the workflow continues; downstream steps skip them.\n")
  dir.create(dirname(args$out_long), showWarnings = FALSE, recursive = TRUE)
  writeLines(paste(c(id_vars, value_vars), collapse = "\t"), args$out_long)
  for (val in value_vars) {
    writeLines(row_var, sub("\\.wide\\.tsv$", sprintf(".%s.mean.wide.tsv", val), args$out_wide))
    writeLines(row_var, sub("\\.wide\\.tsv$", sprintf(".%s.per_sample.filter.wide.tsv", val), args$out_wide))
  }
  quit(save = "no", status = 0)
}

# Validate columns
missing_id <- setdiff(id_vars, names(dt_long))
if (length(missing_id))
  stop("ERROR: id_var column(s) not found: ", paste(missing_id, collapse = ", "))

missing_val <- setdiff(value_vars, names(dt_long))
if (length(missing_val))
  stop("ERROR: value_var column(s) not found: ", paste(missing_val, collapse = ", "))

n_samples <- length(unique(dt_long$sample))
cat(sprintf("\n=== Stacked %d files, %d unique samples, total rows: %d ===\n",
            length(files), n_samples, nrow(dt_long)))

# Combine all metadata
all_meta_dt <- unique(rbindlist(all_meta_list, fill = TRUE))
cat(sprintf("Unique Module metadata entries: %d\n", nrow(all_meta_dt)))

# ------------------------------------------------------------------------------
# 1. Deduplicate within the same sample (sum merge duplicate Module-rank pairs)
# ------------------------------------------------------------------------------
dedup_by_cols <- c("sample", id_vars)
if (row_var %in% names(dt_long)) {
  dedup_by_cols <- c(row_var, dedup_by_cols)
  dedup_by_cols <- unique(dedup_by_cols)
}

dt_long[, (value_vars) := lapply(.SD, sum),
        by = dedup_by_cols,
        .SDcols = value_vars]

dt_long <- unique(dt_long, by = dedup_by_cols)

# Keep only needed columns
keep_cols <- c("sample", id_vars, value_vars)
if (row_var %in% names(dt_long)) {
  keep_cols <- c(row_var, keep_cols)
  keep_cols <- unique(keep_cols)
}
# Also keep metadata columns if present
meta_present <- intersect(module_meta_cols, names(dt_long))
keep_cols <- unique(c(meta_present, keep_cols))
dt_long <- dt_long[, ..keep_cols]

cat(sprintf("After dedup: %d rows\n", nrow(dt_long)))

# ------------------------------------------------------------------------------
# 2. Calculate cross-sample mean
# ------------------------------------------------------------------------------
dt_mean <- dt_long[, lapply(.SD, function(x) sum(x, na.rm = TRUE) / n_samples),
                   by = id_vars,
                   .SDcols = value_vars]

# Merge Module metadata to mean table
if (nrow(all_meta_dt) > 0 && row_var %in% names(all_meta_dt)) {
  dt_mean <- merge(dt_mean, all_meta_dt, by = row_var, all.x = TRUE)
  # Reorder: metadata first
  meta_present <- intersect(module_meta_cols, names(dt_mean))
  other_cols <- setdiff(names(dt_mean), c(row_var, meta_present))
  setcolorder(dt_mean, c(row_var, meta_present, other_cols))
}

fwrite(dt_mean, args$out_long, sep = "\t", quote = FALSE)
cat(sprintf("Mean table saved: %d rows\n", nrow(dt_mean)))

# ------------------------------------------------------------------------------
# 3. Calculate prevalence
# ------------------------------------------------------------------------------
prevalence_dt <- dt_long[, .(prevalence = .N / n_samples), by = id_vars]
dt_mean <- merge(dt_mean, prevalence_dt, by = id_vars)

# ------------------------------------------------------------------------------
# 4. Optional: apply abundance and prevalence filtering
# ------------------------------------------------------------------------------
cat(sprintf("\n=== Pre-filter statistics ===\n"))
cat(sprintf("Total samples: %d\n", n_samples))
cat(sprintf("Total long-format rows: %d\n", nrow(dt_long)))
cat(sprintf("Unique combinations (Module;%s): %d\n", col_var, nrow(dt_mean)))

# Use first value_var for abundance stats
rel_abund_col <- value_vars[1]
if (rel_abund_col %in% names(dt_mean)) {
  cat(sprintf("Mean %s range: %.6f ~ %.6f\n",
              rel_abund_col,
              min(dt_mean[[rel_abund_col]], na.rm = TRUE),
              max(dt_mean[[rel_abund_col]], na.rm = TRUE)))
}
cat(sprintf("Prevalence range: %.2f%% ~ %.2f%%\n",
            min(dt_mean$prevalence, na.rm = TRUE) * 100,
            max(dt_mean$prevalence, na.rm = TRUE) * 100))

do_filter <- args$min_abundance > 0 || args$min_prevalence > 0

if (do_filter) {
  cat(sprintf("\n=== Applying filter ===\n"))
  if (args$min_abundance > 0) {
    cat(sprintf("Minimum mean relative abundance: %.4f%% (%.4f)\n",
                args$min_abundance * 100, args$min_abundance))
  }
  if (args$min_prevalence > 0) {
    cat(sprintf("Minimum prevalence: %.1f%% (%.2f)\n",
                args$min_prevalence * 100, args$min_prevalence))
  }

  keep_idx <- rep(TRUE, nrow(dt_mean))

  if (args$min_abundance > 0 && rel_abund_col %in% names(dt_mean)) {
    keep_idx <- keep_idx & (dt_mean[[rel_abund_col]] >= args$min_abundance)
  } else if (args$min_abundance > 0) {
    cat("WARNING: Abundance column not found, skipping abundance filter\n")
  }

  if (args$min_prevalence > 0) {
    keep_idx <- keep_idx & (dt_mean$prevalence >= args$min_prevalence)
  }

  cat(sprintf("Combinations passing filter: %d / %d (%.1f%%)\n",
              sum(keep_idx), nrow(dt_mean), sum(keep_idx) / nrow(dt_mean) * 100))

  filtered_ids <- dt_mean[keep_idx, ..id_vars]
  dt_mean <- dt_mean[keep_idx]
  dt_long <- dt_long[filtered_ids, on = id_vars]

  cat(sprintf("Long-format rows after filtering: %d\n", nrow(dt_long)))
  cat(sprintf("Unique combinations after filtering: %d\n", nrow(dt_mean)))
} else {
  cat("\n=== No filter applied (min_abundance=0, min_prevalence=0) ===\n")
  cat("Tip: If the output matrix is too large, set filter thresholds, e.g.:\n")
  cat("  --min_abundance 0.0001 --min_prevalence 0.1\n")
}

# ------------------------------------------------------------------------------
# 5. Cast: Module ~ Phylum (mean value matrix)
# ------------------------------------------------------------------------------
for (val in value_vars) {
  dt_wide <- dcast(dt_mean,
                   as.formula(sprintf("%s ~ %s", row_var, col_var)),
                   value.var     = val,
                   fun.aggregate = sum,
                   fill          = 0)

  # Merge metadata to wide
  if (nrow(all_meta_dt) > 0 && row_var %in% names(all_meta_dt)) {
    meta_present <- intersect(module_meta_cols, names(all_meta_dt))
    if (length(meta_present) > 0) {
      dt_wide <- merge(dt_wide, all_meta_dt[, c(row_var, meta_present), with = FALSE],
                       by = row_var, all.x = TRUE)
      other_cols <- setdiff(names(dt_wide), c(row_var, meta_present))
      setcolorder(dt_wide, c(row_var, meta_present, other_cols))
    }
  }

  out_file <- sub("\\.wide\\.tsv$", sprintf(".%s.mean.wide.tsv", val), args$out_wide)
  fwrite(dt_wide, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Written (mean): %s (%d rows x %d cols)\n", out_file, nrow(dt_wide), ncol(dt_wide)))
}

# ------------------------------------------------------------------------------
# 6. Per-sample wide matrix (Module;Taxa ~ sample)
# ------------------------------------------------------------------------------
dt_long[, Module_Taxa := paste(get(row_var), get(col_var), sep = ";")]

for (val in value_vars) {
  n_long_rows <- nrow(dt_long)
  n_unique_pairs <- length(unique(dt_long$Module_Taxa))
  n_unique_samples <- length(unique(dt_long$sample))
  n_wide_rows <- n_unique_pairs
  n_wide_cols <- n_unique_samples + 1

  cat(sprintf("\n=== Generating per-sample wide matrix: %s ===\n", val))
  cat(sprintf("Long-format rows: %d\n", n_long_rows))
  cat(sprintf("Will cast to wide matrix: %d rows x %d cols (Module;Taxa ~ sample)\n", n_wide_rows, n_wide_cols))

  est_memory_mb <- (n_wide_rows * n_wide_cols * 8) / (1024 * 1024)
  cat(sprintf("Estimated memory usage: ~%.1f MB\n", est_memory_mb))

  if (n_wide_rows > 10000 || n_wide_cols > 500) {
    cat("WARNING: Matrix size is large!\n")
    cat("  Suggestion: set stricter filter thresholds\n")
    cat(sprintf("  e.g.: --min_abundance 1 and/or --min_prevalence 0.3\n"))
  }

  dt_wide_per_sample <- dcast(dt_long,
                              Module_Taxa ~ sample,
                              value.var = val,
                              fun.aggregate = sum,
                              fill = 0)

  out_file <- sub("\\.wide\\.tsv$", sprintf(".%s.per_sample.filter.wide.tsv", val), args$out_wide)

  fwrite(dt_wide_per_sample, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Written (per-sample): %s (%d rows x %d cols)\n", out_file, nrow(dt_wide_per_sample), ncol(dt_wide_per_sample)))
}

cat("\nDone.\n")
