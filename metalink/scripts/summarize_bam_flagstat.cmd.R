#!/usr/bin/env Rscript
# ==============================================================================
# summarize_bam_flagstat.cmd.R
# Aggregate per-sample samtools flagstat outputs into one project summary table
# (one row per sample) and a mapping-rate summary figure (PDF).
#
# Usage:
#   Rscript summarize_bam_flagstat.cmd.R \
#     --input_dir <dir with per-sample flagstat tsv files> \
#     --out_tsv <project>.bam_mapping_summary.tsv \
#     --out_pdf <project>.bam_mapping_summary.pdf
# ==============================================================================
suppressMessages(library(optparse))

option_list <- list(
  make_option("--input_dir", type = "character", default = NULL,
              help = "directory containing per-sample samtools flagstat .tsv files"),
  make_option("--out_tsv", type = "character", default = NULL,
              help = "output summary table (one row per sample)"),
  make_option("--out_pdf", type = "character", default = NULL,
              help = "output mapping-rate summary figure (PDF)")
)
opt <- parse_args(OptionParser(option_list = option_list))

files <- list.files(opt$input_dir, pattern = "\\.tsv$", full.names = TRUE)
if (length(files) == 0) {
  stop(paste0("ERROR: no flagstat .tsv files found in ", opt$input_dir))
}

# flagstat line format: "N + 0 metric (P% : N/A)" or "N + 0 metric"
metric_keys <- c("in total", "primary", "secondary", "supplementary",
                 "duplicates", "primary duplicates",
                 "mapped", "primary mapped",
                 "paired in sequencing", "read1", "read2",
                 "properly paired", "with itself and mate mapped", "singletons",
                 "with mate mapped to a different chr",
                 "with mate mapped to a different chr (mapQ>=5)")
col_names <- c("total_reads", "primary", "secondary", "supplementary",
               "duplicates", "primary_duplicates",
               "mapped_reads", "primary_mapped",
               "paired_in_sequencing", "read1", "read2",
               "properly_paired", "itself_and_mate_mapped", "singletons",
               "mate_mapped_diff_chr", "mate_mapped_diff_chr_mapq5")

parse_flagstat <- function(f) {
  lines <- readLines(f, warn = FALSE)
  vals <- setNames(rep(NA_real_, length(metric_keys)), col_names)
  rates <- c(mapping_rate = NA_real_, properly_paired_rate = NA_real_,
             singleton_rate = NA_real_)
  for (ln in lines) {
    # "(mapQ>=5)" is part of the metric name, not a trailing parenthetical
    is_mapq5 <- grepl("(mapQ>=5)", ln, fixed = TRUE)
    m <- regmatches(ln, regexec("^([0-9]+) \\+ [0-9]+ (.*?)( \\((.*)\\))?$", ln))[[1]]
    if (length(m) == 0) next
    cnt <- as.numeric(m[2])
    metric <- if (is_mapq5) "with mate mapped to a different chr (mapQ>=5)" else m[3]
    paren <- m[5]
    pct <- if (!is.na(paren) && grepl("^[0-9.]+%", paren)) {
      suppressWarnings(as.numeric(sub("%.*", "", paren)))
    } else {
      NA_real_
    }
    idx <- match(metric, metric_keys)
    if (!is.na(idx)) vals[idx] <- cnt
    if (metric == "mapped" && !is.na(pct)) rates["mapping_rate"] <- pct
    if (metric == "properly paired" && !is.na(pct)) rates["properly_paired_rate"] <- pct
    if (metric == "singletons" && !is.na(pct)) rates["singleton_rate"] <- pct
  }
  nm <- sub("\\.tsv$", "", basename(f))
  nm <- sub("\\.contigs_bam_stats$|\\.bam_stats$", "", nm)
  c(sample = nm, vals, rates)
}

rows <- lapply(files, parse_flagstat)
df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
num_cols <- setdiff(colnames(df), "sample")
df[num_cols] <- lapply(df[num_cols], as.numeric)
df <- df[order(df$sample), ]

dir.create(dirname(opt$out_tsv), showWarnings = FALSE, recursive = TRUE)
write.table(df, opt$out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(sprintf("Summary table written: %s (%d samples)\n", opt$out_tsv, nrow(df)))

# ---- figure ----
dir.create(dirname(opt$out_pdf), showWarnings = FALSE, recursive = TRUE)
pdf(opt$out_pdf, width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 2), cex = 1.3)

rate <- df$mapping_rate
rate_plot <- rate[!is.na(rate)]
thr <- c(50, 80, 95)
thr_col <- c("#EF5350", "#FF9800", "#26A69A")

# left: sorted mapping-rate curve (a line needs >= 2 points; draw points
# for very small sample sets, a text note when there is no data at all)
if (length(rate_plot) == 0) {
  plot.new()
  text(0.5, 0.5, "no mapping-rate data available")
} else {
  plot(sort(rate_plot),
       type = if (length(rate_plot) > 1) "l" else "p",
       lwd = 1.5, pch = 19, xlab = "Samples (sorted by mapping rate)",
       ylab = "Mapping rate (%)", ylim = c(0, 100),
       main = sprintf("Mapping rate across %d sample(s)", length(rate_plot)))
  if (length(rate_plot) == 1) points(1, rate_plot, pch = 19)
  abline(h = thr, col = thr_col, lty = 2, lwd = 2)
  legend("bottomright", legend = paste0(thr, "%"), col = thr_col, lty = 2, bty = "n")
}

# right: sample counts per bin
bins <- cut(rate_plot, breaks = c(-Inf, 50, 80, 95, Inf),
            labels = c("<50%", "50-80%", "80-95%", ">=95%"))
counts <- table(bins)
bp <- barplot(counts, col = adjustcolor(c("#EF5350", "#FF9800", "#26A69A", "#42A5F5"), alpha.f = 0.7),
              ylab = "Number of samples", xlab = "Mapping rate",
              main = "Samples per mapping-rate bin", ylim = c(0, max(counts) * 1.2 + 1))
text(bp, counts, labels = counts, pos = 3, cex = 1.1)

dev.off()
cat(sprintf("Summary figure written: %s\n", opt$out_pdf))
