#!/usr/bin/env Rscript
# ==============================================================================
# summarize_MAG_quality.cmd.R
# Aggregate per-sample checkm stats into one project-level MAG quality table
# (one row per MAG) and a completeness-vs-contamination figure with marginal
# histograms, colored by MIMAG quality tier.
#
# Usage:
#   Rscript summarize_MAG_quality.cmd.R \
#     --input_dir <dir with per-sample checkm_stats csv files> \
#     --out_tsv <project>.MAG_quality_summary.tsv \
#     --out_pdf <project>.MAG_quality_summary.pdf
# ==============================================================================
suppressMessages(library(optparse))

option_list <- list(
  make_option("--input_dir", type = "character", default = NULL,
              help = "directory containing per-sample checkm_stats .csv files"),
  make_option("--out_tsv", type = "character", default = NULL,
              help = "output summary table (one row per MAG)"),
  make_option("--out_pdf", type = "character", default = NULL,
              help = "output quality-distribution figure (PDF)")
)
opt <- parse_args(OptionParser(option_list = option_list))

files <- list.files(opt$input_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  stop(paste0("ERROR: no checkm_stats .csv files found in ", opt$input_dir))
}

dfs <- lapply(files, function(f) {
  d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(d) == 0) {
    cat(sprintf("NOTE: %s has no MAG rows (sample passed no bins); skipped\n", basename(f)))
    return(NULL)
  }
  d$sample <- sub("\\.csv$", "", basename(f))
  d
})
dfs <- dfs[!vapply(dfs, is.null, logical(1))]
if (length(dfs) == 0) {
  stop(paste0("ERROR: all checkm_stats files in ", opt$input_dir, " are empty (no MAGs in any sample)"))
}
df <- do.call(rbind, dfs)
df <- df[, c("genome", "sample", "completeness", "contamination",
             "strain_heterogeneity", "length", "N50", "QS")]

# MIMAG tiers (Medium contamination boundary extended to <= 10)
df$quality_tier <- ifelse(df$completeness >= 90 & df$contamination <= 5, "High",
                   ifelse(df$completeness >= 50 & df$contamination <= 10, "Medium",
                          "Partial"))

dir.create(dirname(opt$out_tsv), showWarnings = FALSE, recursive = TRUE)
write.table(df, opt$out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Summary table written: %s (%d MAGs from %d samples)\n",
            opt$out_tsv, nrow(df), length(files)))

# ---- figure: scatter + marginal histograms (base R, no extra packages) ----
tier_levels <- c("High", "Medium", "Partial")
tier_col <- c(High = "#F88800", Medium = "#6898E8", Partial = "#307020")
# unname(): a NAMED color vector silently breaks barplot/layout rendering
tier_col_a <- unname(adjustcolor(tier_col, alpha.f = 0.5))
tier <- factor(df$quality_tier, levels = tier_levels)
counts <- table(tier)
legend_txt <- sprintf("%s quality (%d genomes; %.2f%%)", tier_levels,
                      as.integer(counts), 100 * as.numeric(counts) / nrow(df))

dir.create(dirname(opt$out_pdf), showWarnings = FALSE, recursive = TRUE)
pdf(opt$out_pdf, width = 9, height = 8)
layout(matrix(c(1, 0, 2, 3), 2, 2, byrow = TRUE), widths = c(4, 1.1),
       heights = c(1.1, 4))

xlim <- c(50, 100)
ylim <- c(0, max(10, ceiling(max(df$contamination))))

# top: completeness histogram
par(mar = c(0, 5, 2, 0))
hh <- hist(df$completeness, breaks = seq(xlim[1], xlim[2], by = 2.5),
           plot = FALSE)
counts_x <- lapply(tier_levels, function(t)
  hist(df$completeness[tier == t], breaks = hh$breaks, plot = FALSE)$counts)
counts_x <- do.call(rbind, counts_x)
# barplot() draws borders even for zero-count segments (as stray lines in the
# wrong color), so draw fills borderless and add opaque borders manually
mp <- barplot(counts_x, beside = FALSE, col = tier_col_a,
              border = NA, axes = FALSE, space = 0)
csum <- rbind(0, apply(counts_x, 2, cumsum))
for (i in seq_along(mp)) {
  for (j in seq_along(tier_levels)) {
    if (counts_x[j, i] > 0)
      rect(mp[i] - 0.5, csum[j, i], mp[i] + 0.5, csum[j + 1, i],
           border = unname(tier_col)[j])
  }
}

# bottom-left: scatter (horizontal tick labels via las = 1; axis titles unchanged)
par(mar = c(6.5, 6.5, 0, 0), cex.lab = 2.0, cex.axis = 1.6, las = 1)
plot(df$completeness, df$contamination, col = tier_col_a[as.integer(tier)], pch = 16,
     cex = 0.7, xlim = xlim, ylim = ylim, xlab = "Completeness (%)",
     ylab = "Contamination (%)")
abline(v = c(90, 50), h = c(5, 10), lty = 2, col = "grey50")
legend("topleft", legend = legend_txt, bty = "n", cex = 1.6)

# right: contamination histogram (horizontal)
par(mar = c(5, 0, 0, 2))
hv <- hist(df$contamination, breaks = seq(ylim[1], ylim[2], length.out = 25),
           plot = FALSE)
counts_y <- lapply(tier_levels, function(t)
  hist(df$contamination[tier == t], breaks = hv$breaks, plot = FALSE)$counts)
counts_y <- do.call(rbind, counts_y)
mp <- barplot(counts_y, beside = FALSE, col = tier_col_a,
              border = NA, axes = FALSE, space = 0, horiz = TRUE)
csum <- rbind(0, apply(counts_y, 2, cumsum))
for (i in seq_along(mp)) {
  for (j in seq_along(tier_levels)) {
    if (counts_y[j, i] > 0)
      rect(csum[j, i], mp[i] - 0.5, csum[j + 1, i], mp[i] + 0.5,
           border = unname(tier_col)[j])
  }
}

dev.off()
cat(sprintf("Summary figure written: %s\n", opt$out_pdf))
