#!/usr/bin/env Rscript
# Filter module abundance by completeness threshold
# Usage: Rscript filter_by_completeness.R [options]
# Example: Rscript filter_by_completeness.R -a mc_abundance.tsv -c mc_completeness.tsv -o filtered.tsv -t 0.5

# ============================================================================
# 0. LOAD REQUIRED PACKAGES
# ============================================================================

# Check and install optparse if not available
if (!requireNamespace("optparse", quietly = TRUE)) {
  cat("Package 'optparse' not found. Installing from CRAN...\n")
  install.packages("optparse", repos = "https://cloud.r-project.org/")
}

suppressMessages(library(optparse))
suppressMessages(library(data.table))

# ============================================================================
# 1. COMMAND-LINE ARGUMENT PARSING
# ============================================================================

option_list <- list(
  make_option(c("-a", "--abundance_RelativeAbundance"),
              type = "character",
              dest = "abundance_RelativeAbundance",
              default = NA_character_,
              metavar = "FILE",
              help = "Input abundance table of the primary track (RelativeAbundance scale, tsv format). Required."),

  make_option("--abundance_SumCoverage",
              type = "character",
              dest = "abundance_SumCoverage",
              default = NA_character_,
              metavar = "FILE",
              help = "Optional second abundance table (SumCoverage scale, tsv format); filtered with the same completeness table. Must be used together with --output_SumCoverage."),

  make_option(c("-c", "--completeness"),
              type = "character",
              dest = "completeness",
              default = NA_character_,
              metavar = "FILE",
              help = "Input completeness table (tsv format). Required. Shared by all abundance tracks."),

  make_option(c("-o", "--output_RelativeAbundance"),
              type = "character",
              dest = "output_RelativeAbundance",
              default = "mc_abundance_filtered.tsv",
              metavar = "FILE",
              help = "Output filtered abundance table for the RelativeAbundance track [default: %default]"),

  make_option("--output_SumCoverage",
              type = "character",
              dest = "output_SumCoverage",
              default = NA_character_,
              metavar = "FILE",
              help = "Output filtered abundance table for the SumCoverage track. Must be used together with --abundance_SumCoverage."),

  make_option(c("-t", "--completeness-threshold"),
              type = "double",
              dest = "completeness_threshold",
              default = 0.5,
              metavar = "FLOAT",
              help = "Completeness threshold (0-1). Abundance values with completeness below this are set to zero [default: %default]"),

  make_option("--meta-cols",
              type = "character",
              dest = "meta_cols",
              default = "Module_Entry,Level_2,Level_3,Module_Name,Definition",
              metavar = "STRING",
              help = "Comma-separated metadata column names [default: %default]"),

  make_option("--strategy",
              type = "character",
              dest = "strategy",
              default = "zero",
              metavar = "STRING",
              help = "Filtering strategy: 'zero' (set to 0) or 'na' (set to NA) [default: %default]"),

  make_option("--header",
              action = "store_true",
              dest = "header",
              default = TRUE,
              help = "Input tsv files have a header row [default: TRUE]"),

  make_option(c("-q", "--quiet"),
              action = "store_true",
              dest = "quiet",
              default = FALSE,
              help = "Suppress non-essential output [default: FALSE]")
)

parser <- OptionParser(
  description = "Filter module abundance by completeness threshold. Values in the abundance table whose corresponding completeness is below the threshold are set to zero.",
  epilogue = "Example: Rscript filter_by_completeness.R -a abundance.tsv -c completeness.tsv -o filtered.tsv -t 0.5",
  option_list = option_list
)

args <- parse_args(parser, positional_arguments = FALSE)

# Manual required argument checks (optparse has no built-in required)
if (is.na(args$abundance_RelativeAbundance)) {
  stop("Missing required argument: -a/--abundance_RelativeAbundance")
}
if (is.na(args$completeness)) {
  stop("Missing required argument: -c/--completeness")
}
if (xor(!is.na(args$abundance_SumCoverage), !is.na(args$output_SumCoverage))) {
  stop("--abundance_SumCoverage and --output_SumCoverage must be provided together")
}

# Manual choices check (optparse has no built-in choices)
if (!args$strategy %in% c("zero", "na")) {
  stop(sprintf("Invalid strategy '%s'. Must be one of: zero, na", args$strategy))
}

# ============================================================================
# 2. UNPACK ARGUMENTS
# ============================================================================

abundance_RA_file   <- args$abundance_RelativeAbundance
abundance_SC_file   <- args$abundance_SumCoverage
completeness_file   <- args$completeness
output_RA_file      <- args$output_RelativeAbundance
output_SC_file      <- args$output_SumCoverage
completeness_threshold         <- args$completeness_threshold
strategy            <- args$strategy
has_header          <- args$header
quiet_mode          <- args$quiet

# Parse metadata columns
meta_cols <- strsplit(args$meta_cols, ",")[[1]]
meta_cols <- trimws(meta_cols)  # Remove whitespace

# ============================================================================
# 3. HELPER: CONDITIONAL PRINTING
# ============================================================================

log_msg <- function(...) {
  if (!quiet_mode) {
    msg <- sprintf(...)
    if (nchar(msg) > 0) {
      ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      cat(sprintf("[%s] %s", ts, msg), "\n")
    } else {
      cat("\n")
    }
  }
}

# ============================================================================
# 4. READ COMPLETENESS TABLE (shared by all abundance tracks)
# ============================================================================

log_msg("Loading completeness table: %s", completeness_file)
dt_comp  <- fread(completeness_file, sep = "\t", header = has_header, stringsAsFactors = FALSE)

# ============================================================================
# 5. FILTERING FUNCTION (shared by all abundance scales)
# ============================================================================

filter_abundance_by_completeness <- function(abundance_file, output_file) {

  log_msg("Loading abundance table: %s", abundance_file)
  dt_abund <- fread(abundance_file, sep = "\t", header = has_header, stringsAsFactors = FALSE)

  # Verify structure consistency
  if (!identical(names(dt_abund), names(dt_comp))) {
    stop(sprintf("Column names mismatch!\nAbundance: %s\nCompleteness: %s",
                 paste(names(dt_abund), collapse = ", "),
                 paste(names(dt_comp), collapse = ", ")))
  }

  if (nrow(dt_abund) != nrow(dt_comp)) {
    stop(sprintf("Row count mismatch! Abundance: %d, Completeness: %d",
                 nrow(dt_abund), nrow(dt_comp)))
  }

  missing_meta <- setdiff(meta_cols, names(dt_abund))
  if (length(missing_meta) > 0) {
    stop(sprintf("Metadata columns not found in input: %s",
                 paste(missing_meta, collapse = ", ")))
  }

  log_msg("Tables loaded: %d rows x %d columns", nrow(dt_abund), ncol(dt_abund))
  log_msg("Metadata columns: %s", paste(meta_cols, collapse = ", "))
  log_msg("Sample columns: %d", ncol(dt_abund) - length(meta_cols))

  # Extract sample columns
  sample_cols <- setdiff(names(dt_abund), meta_cols)

  if (length(sample_cols) == 0) {
    stop("No sample columns found after removing metadata columns. Check --meta-cols.")
  }

  # Extract numeric matrices
  mat_abund <- as.matrix(dt_abund[, ..sample_cols])
  mat_comp  <- as.matrix(dt_comp[, ..sample_cols])

  mode(mat_abund) <- "numeric"
  mode(mat_comp)  <- "numeric"

  n_na_abund <- sum(is.na(mat_abund))
  n_na_comp  <- sum(is.na(mat_comp))
  if (n_na_abund > 0) {
    log_msg("WARNING: abundance table %s contains %d NA value(s), please check the upstream data", abundance_file, n_na_abund)
  }
  if (n_na_comp > 0) {
    log_msg("WARNING: completeness table %s contains %d NA value(s); the corresponding abundance values are kept as-is (not filtered), please check the upstream data", completeness_file, n_na_comp)
  }

  # Filtering
  mat_abund_filtered <- mat_abund
  mask_low <- mat_comp < completeness_threshold
  mask_low[is.na(mask_low)] <- FALSE

  n_total   <- length(mat_abund)
  n_filtered <- sum(mask_low, na.rm = TRUE)

  log_msg("Completeness threshold: %.3f", completeness_threshold)
  log_msg("Strategy: %s", strategy)
  log_msg("Total values: %d", n_total)
  log_msg("Values below completeness threshold: %d (%.2f%%)", n_filtered, 100 * n_filtered / n_total)

  # Apply filtering
  if (strategy == "zero") {
    mat_abund_filtered[mask_low] <- 0
  } else if (strategy == "na") {
    mat_abund_filtered[mask_low] <- NA_real_
  }

  # Construct output
  dt_out <- cbind(
    dt_abund[, ..meta_cols],
    as.data.table(mat_abund_filtered)
  )
  setnames(dt_out, c(meta_cols, sample_cols))

  # Save output
  fwrite(dt_out, output_file, sep = "\t", row.names = FALSE)
  log_msg("Output saved: %s", output_file)
  log_msg("Dimensions: %d rows x %d columns", nrow(dt_out), ncol(dt_out))

  # Summary statistics
  if (!quiet_mode) {
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] --- Filtering Summary ---\n", ts))
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] Original total abundance:  %.2f\n", ts, sum(mat_abund, na.rm = TRUE)))
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] Filtered total abundance:  %.2f\n", ts, sum(mat_abund_filtered, na.rm = TRUE)))
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] Original non-zero values:  %d\n", ts, sum(mat_abund > 0, na.rm = TRUE)))
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] Filtered non-zero values:  %d\n", ts, sum(mat_abund_filtered > 0, na.rm = TRUE)))

    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] --- Per-sample retention ---\n", ts))
    for (sc in sample_cols) {
      orig <- sum(dt_abund[[sc]] > 0, na.rm = TRUE)
      filt <- sum(dt_out[[sc]] > 0, na.rm = TRUE)
      pct  <- ifelse(orig > 0, 100 * filt / orig, 0)
      ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      cat(sprintf("[%s]   %s: %d -> %d (%.1f%%)\n", ts, sc, orig, filt, pct))
    }
  }

  invisible(dt_out)
}

# ============================================================================
# 6. RUN FILTERING
# ============================================================================

filter_abundance_by_completeness(abundance_RA_file, output_RA_file)

if (!is.na(abundance_SC_file)) {
  log_msg("Filtering the SumCoverage-scale abundance table with the same completeness table")
  filter_abundance_by_completeness(abundance_SC_file, output_SC_file)
}
