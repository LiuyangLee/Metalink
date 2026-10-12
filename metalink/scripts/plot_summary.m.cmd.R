#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(tidyverse)
  library(dplyr)
  library(ggrepel)
  library(data.table)
})

# -------------------------- command-line parser --------------------------
option_list <- list(
  make_option(c("-t", "--Tax"),  type = "character", default = "MAG_Tax_summary.tsv",
              help = "path to MAG_Tax_summary.tsv [default %default]"),
  make_option(c("-k", "--ORFKO"),      type = "character", default = "ORF_KO_summary.tsv",
              help = "path to ORF_KO_summary.tsv [default %default]"),
  make_option(c("-a", "--ORFTaxKO"),   type = "character", default = "ORF_Tax_KO_summary.tsv",
              help = "path to ORF_Tax_KO_summary.tsv [default %default]"),
  make_option(c("-o", "--outPrefix"),  type = "character", default = "plot",
              help = "output file prefix [default %default]"),
  make_option(c("-f", "--format"),     type = "character", default = "png",
              help = "output image format: png|pdf|svg [default %default]"),
  make_option(c("-w", "--width"),      type = "numeric",   default = 8,
              help = "plot width in inches [default %default]"),
  make_option(c("-h", "--height"),     type = "numeric",   default = 5,
              help = "plot height in inches [default %default]"),
  make_option(c("-s", "--fontSize"),   type = "numeric",   default = 13,
              help = "base font size [default %default]"),
  make_option(c("-c", "--colorSet"),   type = "character",
              default = "#0288D1,#FF9800,#26A69A,#E6EE9C,#00BCD4,#EF5350,#9CCC65,#C5CCFF,#9E9E9E",
              help = "Color specification: either 'brewer:PaletteName' (e.g., brewer:Set3) or comma-separated hex colors. Default matches Set3 palette [default %default]"),
  make_option(c("--alpha"),             type = "numeric",   default = 0.7,
              help = "bar transparency (0-1) [default %default]"),
  make_option(c("-n", "--noPlot"),     action = "store_true", default = FALSE,
              help = "skip plotting, only prepare data")
)
#42A5F5,
#FF9800,
#26A69A,
#E6EE9C,
#EF5350,
#00BCD4,
#9CCC65,
#C5CCFF,
#9E9E9E
opt <- parse_args(OptionParser(option_list = option_list,add_help_option=FALSE))

# -------------------------- helper: save plot ----------------------------
save_plot <- function(p, suffix, width_factor, height_factor = 1) {
  if (opt$noPlot) return(invisible(NULL))
  outfile <- paste0(opt$outPrefix, "_", suffix, ".", opt$format)
  ggplot2::ggsave(outfile, plot = p, width = opt$width * width_factor, height = opt$height * height_factor,create.dir = TRUE)
  message("Saved: ", outfile)
}

# -------------------------- helper: color scale ----------------------------
get_fill_scale <- function(color_set, n = 9) {
  if (grepl("^brewer:", color_set)) {
    palette_name <- sub("^brewer:", "", color_set)
    return(scale_fill_brewer(palette = palette_name))
  } else {
    colors <- strsplit(color_set, ",")[[1]]
    colors <- trimws(colors)
    if (length(colors) < n) {
      warning(sprintf("Only %d colors provided, need %d. Recycling...", length(colors), n))
      colors <- rep(colors, length.out = n)
    }
    return(scale_fill_manual(values = colors[1:n]))
  }
}

# -------------------------- load data ------------------------------------
MAG_Tax_summary <- fread(opt$Tax) %>% as_tibble()
ORF_KO_summary     <- fread(opt$ORFKO)     %>% as_tibble()
ORF_Tax_KO_summary <- fread(opt$ORFTaxKO)  %>% as_tibble()

# -------------------------- fixed levels ---------------------------------
levels_Domain      <- c("All", "d__Bacteria", "d__Archaea")
levels_NewTax      <- c("Known species", "Unclassified species", "Unclassified genus",
                        "Unclassified family", "Unclassified order", "Unclassified class",
                        "Unclassified phylum", "Unclassified")
levels_KO_annotate <- c("Classified", "Unclassified")

# -------------------------- MAG_Tax_summary plots ---------------------
MAG_Tax_summary <- MAG_Tax_summary %>%
  mutate(Domain = factor(Domain, levels = levels_Domain),
         NewTax = factor(NewTax, levels = levels_NewTax)) %>%
  complete(Domain, NewTax, fill = list(MAG_Num = 0, Proportion = 0))


p1 <- ggplot(MAG_Tax_summary,
             aes(x = Domain, y = MAG_Num, fill = NewTax)) +
  geom_col(position = position_dodge(preserve = "single"), alpha = opt$alpha) +
  geom_text(aes(label = paste0(sprintf("%.0f",MAG_Num),"")),
            position = position_dodge(width = 0.9),color = "black", vjust = -0.5, size=2, check_overlap = T) +
  get_fill_scale(opt$colorSet) +
  labs(title = "MAG Count by Domain and NewTax",
       x = "Domain", y = "MAG Count") +
  theme_minimal(base_size = opt$fontSize) +
  theme(axis.ticks.length = unit(0.2, "cm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

save_plot(p1, "MAG_Num", width_factor = 1.6)


p2 <- ggplot(MAG_Tax_summary,
             aes(x = Domain, y = Proportion, fill = NewTax)) +
  geom_col(position = position_dodge(preserve = "single"), alpha = opt$alpha) +
  geom_text(aes(label = sprintf("%.0f%%", Proportion * 100)),
            position = position_dodge(width = 0.9),color = "black", vjust = -0.5, size=2, check_overlap = T) +
  get_fill_scale(opt$colorSet) +
  labs(title = "MAG Proportion by Domain and NewTax",
       x = "Domain", y = "MAG Proportion") +
  theme_minimal(base_size = opt$fontSize) +
  theme(axis.ticks.length = unit(0.2, "cm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

save_plot(p2, "MAG_prop", width_factor = 1.6)

# -------------------------- ORF_KO_summary plots -------------------------
ORF_KO_summary <- ORF_KO_summary %>%
  mutate(Domain = factor(Domain, levels = levels_Domain),
         KO_annotate = factor(KO_annotate, levels = levels_KO_annotate)) %>%
  complete(Domain, KO_annotate, fill = list(ORF_Num = 0, Proportion = 0))

p3 <- ggplot(ORF_KO_summary,
             aes(x = Domain, y = ORF_Num, fill = KO_annotate)) +
  geom_col(position = position_dodge(preserve = "single"), alpha = opt$alpha) +
  geom_text(aes(label = paste0(sprintf("%.0f",ORF_Num),"")),
            position = position_dodge(width = 0.9),color = "black", vjust = -0.5, size=2, check_overlap = T) +
  get_fill_scale(opt$colorSet, n = 2) +
  labs(title = "ORF Count by Domain and KO_annotate",
       x = "Domain", y = "ORF Count") +
  theme_minimal(base_size = opt$fontSize) +
  theme(axis.ticks.length = unit(0.2, "cm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

save_plot(p3, "ORFKO_count", width_factor = 1)

p4 <- ggplot(ORF_KO_summary,
             aes(x = Domain, y = Proportion, fill = KO_annotate)) +
  geom_col(position = position_dodge(preserve = "single"), alpha = opt$alpha) +
  geom_text(aes(label = sprintf("%.0f%%", Proportion * 100)),
            position = position_dodge(width = 0.9),color = "black", vjust = -0.5, size=2, check_overlap = T) +
  get_fill_scale(opt$colorSet, n = 2) +
  labs(title = "ORF Proportion by Domain and KO_annotate",
       x = "Domain", y = "ORF Proportion") +
  theme_minimal(base_size = opt$fontSize) +
  theme(axis.ticks.length = unit(0.2, "cm"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

save_plot(p4, "ORFKO_prop", width_factor = 1)

# -------------------------- ORF_Tax_KO_summary plots ---------------------
levels_Domain_NewTax = vector()
for (x in levels_Domain[1:3]){
  for (y in levels_NewTax){
    levels_Domain_NewTax = c(levels_Domain_NewTax, paste0(x,"_", y))
  }
}

ORF_Tax_KO_summary <- ORF_Tax_KO_summary %>%
  {.[.$Domain %in% levels_Domain[1:3],]} %>%
  mutate(Domain = factor(Domain, levels = levels_Domain[1:3]),
         NewTax = factor(NewTax, levels = levels_NewTax),
         KO_annotate = factor(KO_annotate, levels = levels_KO_annotate)) %>%
  complete(Domain, NewTax, KO_annotate,fill = list(ORF_Num = 0, Proportion = 0)) %>%
  mutate(Domain_NewTax = factor(paste(Domain, NewTax, sep = "_"),
                                levels = levels_Domain_NewTax))

p5 <- ggplot(ORF_Tax_KO_summary,
             aes(x = Domain_NewTax, y = ORF_Num, fill = KO_annotate)) +
  geom_col(position = position_dodge(preserve = "single"), alpha = opt$alpha) +
  geom_text(aes(label = paste0(sprintf("%.0f",ORF_Num),"")),
            position = position_dodge(width = 0.9),color = "black", vjust = -0.5, size=2, check_overlap = T) +
  get_fill_scale(opt$colorSet, n = 2) +
  labs(title = "ORF Count by Domain_NewTax and KO_annotate",
       x = "Domain_NewTax", y = "ORF Count") +
  theme_minimal(base_size = opt$fontSize) +
  theme(axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

save_plot(p5, "ORFTaxKO_count",  width_factor = 2.4, height_factor = 1.2)

p6 <- ggplot(ORF_Tax_KO_summary,
             aes(x = Domain_NewTax, y = Proportion, fill = KO_annotate)) +
  geom_col(position = position_dodge(preserve = "single"), alpha = opt$alpha) +
  geom_text(aes(label = sprintf("%.0f%%", Proportion * 100)),
            position = position_dodge(width = 0.9),color = "black", vjust = -0.5, size=2, check_overlap = T) +
  get_fill_scale(opt$colorSet, n = 2) +
  labs(title = "ORF Proportion by Domain_NewTax and KO_annotate",
       x = "Domain_NewTax", y = "ORF Proportion") +
  theme_minimal(base_size = opt$fontSize) +
  theme(axis.ticks.length = unit(0.2, "cm"),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))

save_plot(p6, "ORFTaxKO_prop",  width_factor = 2.4, height_factor = 1.2)

message("All plots saved.")
