#!/usr/bin/env Rscript plot_summary.cmd.R
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
  make_option(c("-t", "--Tax"),  type = "character", default = "contig_Tax_summary.tsv",
              help = "path to contig_Tax_summary.tsv [default %default]"),
  make_option(c("--TaxFormat"),  type = "character", default = "c",
              help = "Tax file format: c|read|pmr|mc|cr|vc|vm|mcr|vmc|vcr|vmcr [default %default]"),
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

opt <- parse_args(OptionParser(option_list = option_list,add_help_option=FALSE))

# -------------------------- helper: save plot ----------------------------
save_plot <- function(p, suffix, width_factor, height_factor = 1) {
  if (opt$noPlot) return(invisible(NULL))
  outfile <- paste0(opt$outPrefix, "_", suffix, ".", opt$format)
  ggplot2::ggsave(outfile, plot = p, width = opt$width * width_factor, height = opt$height * height_factor,create.dir = TRUE)
  message("Saved: ", outfile)
}

# -------------------------- helper: format big numbers -------------------
format_big_num <- function(x) {
  ifelse(x >= 1e9, sprintf("%.0fG", x / 1e9),
         ifelse(x >= 1e6, sprintf("%.0fM", x / 1e6),
                ifelse(x >= 1e3, sprintf("%.0fk", x / 1e3),
                       sprintf("%.0f", x))))
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

# -------------------------- format definitions ---------------------------
format_defs <- list(
  c = list(
    c = c("Contig_Num", "Total_Contig", "Proportion")
  ),
  read = list(
    read = c("Read_Num", "Total_Read", "Proportion")
  ),
  pmr = list(
    MAG = c("MAG_Num", "Total_MAG", "Proportion")
  ),
  mc = list(
    mc  = c("Contig_Num_mc", "Total_Contig_mc", "Proportion_Contig_mc"),
    c   = c("Contig_Num_c", "Total_Contig_c", "Proportion_Contig_c"),
    MAG = c("MAG_Num", "Total_MAG", "Proportion_m")
  ),
  cr = list(
    cr = c("Seq_Bases_cr", "Total_Seq_cr", "Proportion_cr"),
    r  = c("Read_Num_r", "Total_Read_r", "Proportion_r"),
    c  = c("Contig_Num_c", "Total_Contig_c", "Proportion_Contig_c")
  ),
  vc = list(
    vc   = c("Contig_Num_vc", "Total_Contig_vc", "Proportion_Contig_vc"),
    c    = c("Contig_Num_c", "Total_Contig_c", "Proportion_Contig_c"),
    vMAG = c("vMAG_Num", "Total_vMAG", "Proportion_v")
  ),
  vm = list(
    vm   = c("Contig_Num_vm", "Total_Contig_vm", "Proportion_Contig_vm"),
    mMAG = c("mMAG_Num", "Total_mMAG", "Proportion_m"),
    vMAG = c("vMAG_Num", "Total_vMAG", "Proportion_v")
  ),
  mcr = list(
    mcr  = c("Seq_Bases_mcr", "Total_Seq_mcr", "Proportion_mcr"),
    r    = c("Read_Num_r", "Total_Read_r", "Proportion_r"),
    mc   = c("Contig_Num_mc", "Total_Contig_mc", "Proportion_Contig_mc"),
    c    = c("Contig_Num_c", "Total_Contig_c", "Proportion_Contig_c"),
    MAG  = c("MAG_Num", "Total_MAG", "Proportion_m")
  ),
  vmc = list(
    vmc  = c("Contig_Num_vmc", "Total_Contig_vmc", "Proportion_Contig_vmc"),
    mc   = c("Contig_Num_mc", "Total_Contig_mc", "Proportion_Contig_mc"),
    c    = c("Contig_Num_c", "Total_Contig_c", "Proportion_Contig_c"),
    mMAG = c("mMAG_Num", "Total_mMAG", "Proportion_m"),
    vMAG = c("vMAG_Num", "Total_vMAG", "Proportion_v")
  ),
  vcr = list(
    vcr  = c("Seq_Bases_vcr", "Total_Seq_vcr", "Proportion_vcr"),
    r    = c("Read_Num_r", "Total_Read_r", "Proportion_r"),
    vc   = c("Contig_Num_vc", "Total_Contig_vc", "Proportion_Contig_vc"),
    c    = c("Contig_Num_c", "Total_Contig_c", "Proportion_Contig_c"),
    vMAG = c("vMAG_Num", "Total_vMAG", "Proportion_v")
  ),
  vmcr = list(
    vmcr = c("Seq_Bases_vmcr", "Total_Seq_vmcr", "Proportion_vmcr"),
    r    = c("Read_Num_r", "Total_Read_r", "Proportion_r"),
    vmc  = c("Contig_Num_vmc", "Total_Contig_vmc", "Proportion_Contig_vmc"),
    mc   = c("Contig_Num_mc", "Total_Contig_mc", "Proportion_Contig_mc"),
    c    = c("Contig_Num_c", "Total_Contig_c", "Proportion_Contig_c"),
    mMAG = c("mMAG_Num", "Total_mMAG", "Proportion_m"),
    vMAG = c("vMAG_Num", "Total_vMAG", "Proportion_v")
  )
)

if (!(opt$TaxFormat %in% names(format_defs))) {
  stop("Unknown TaxFormat: ", opt$TaxFormat, ". Supported: ", paste(names(format_defs), collapse = ", "))
}

source_map <- format_defs[[opt$TaxFormat]]
is_multi_source <- length(source_map) > 1

# -------------------------- load data ------------------------------------
contig_Tax_summary <- fread(opt$Tax) %>% as_tibble() %>%
  filter(!is.na(NewTax), !is.na(Domain))  # filter out NA rows before factorization

ORF_KO_summary     <- fread(opt$ORFKO)     %>% as_tibble()
ORF_Tax_KO_summary <- fread(opt$ORFTaxKO)  %>% as_tibble()

# -------------------------- fixed levels ---------------------------------
levels_Domain      <- c("All", "d__Bacteria", "d__Archaea", "d__Viruses",
                        "d__Eukaryota", "d__unclassified", "d__root domain")
levels_NewTax      <- c("Known species", "Unclassified species", "Unclassified genus",
                        "Unclassified family", "Unclassified order", "Unclassified class",
                        "Unclassified phylum", "Root", "Unclassified")
levels_KO_annotate <- c("Classified", "Unclassified")

# -------------------------- helper: infer count type from column name ----
infer_count_type <- function(num_col) {
  if (grepl("Read_Num", num_col)) {
    return("Read_Num")
  } else if (grepl("Contig_Num", num_col)) {
    return("Contig_Num")
  } else if (num_col == "MAG_Num") {
    return("MAG_Num")
  } else if (num_col == "mMAG_Num") {
    return("mMAG_Num")
  } else if (num_col == "vMAG_Num") {
    return("vMAG_Num")
  } else if (grepl("Seq_Bases", num_col)) {
    return("Seq_Bases")
  } else {
    return("Count")
  }
}

# -------------------------- reshape data ---------------------------------
if (is_multi_source) {
  # Multi-source: convert to long format
  long_list <- list()
  source_types <- list()
  for (src in names(source_map)) {
    cols <- source_map[[src]]
    if (all(cols %in% colnames(contig_Tax_summary))) {
      count_type <- infer_count_type(cols[1])
      source_types[[src]] <- count_type
      df_src <- contig_Tax_summary %>%
        select(Domain, NewTax, all_of(cols)) %>%
        rename(Count = !!sym(cols[1]),
               Total = !!sym(cols[2]),
               Proportion = !!sym(cols[3])) %>%
        mutate(Source = src, Type = count_type)
      long_list[[src]] <- df_src
    } else {
      missing <- setdiff(cols, colnames(contig_Tax_summary))
      warning("Source '", src, "' skipped: missing columns ", paste(missing, collapse = ", "))
    }
  }

  if (length(long_list) == 0) {
    stop("No valid sources found in the input file for format '", opt$TaxFormat, "'")
  }

  contig_Tax_summary <- bind_rows(long_list)

  # Only use sources that actually exist in the data
  actual_sources <- names(source_types)
  contig_Tax_summary$Source <- factor(contig_Tax_summary$Source, levels = actual_sources)

  # Apply factor levels and filtering
  contig_Tax_summary <- contig_Tax_summary %>%
    mutate(Domain = factor(Domain, levels = levels_Domain),
           NewTax = factor(NewTax, levels = levels_NewTax)) %>%
    complete(Domain, NewTax, Source, fill = list(Count = 0, Proportion = 0)) %>%
    filter(!is.na(NewTax)) %>%  # double-check: remove any NA NewTax after complete
    {.[.$Domain %in% levels_Domain[1:5], ]}

  # Create facet labels with type info: "Source (Type)"
  facet_labels <- sapply(actual_sources, function(s) {
    paste0(s, " (", source_types[[s]], ")")
  })
  contig_Tax_summary$FacetLabel <- factor(
    contig_Tax_summary$Source,
    levels = actual_sources,
    labels = facet_labels
  )

  # Remove any rows with NA FacetLabel (shouldn't happen, but just in case)
  contig_Tax_summary <- contig_Tax_summary %>% filter(!is.na(FacetLabel))

  # Determine labels for multi-source plot
  count_ylabel <- "Count"
  count_title  <- "Count by Domain and NewTax"
  prop_ylabel  <- "Proportion"
  prop_title   <- "Proportion by Domain and NewTax"

} else {
  # Standard format
  cols <- source_map[[1]]
  count_type <- infer_count_type(cols[1])
  contig_Tax_summary <- contig_Tax_summary %>%
    rename(Count = !!sym(cols[1]),
           Total = !!sym(cols[2]),
           Proportion = !!sym(cols[3])) %>%
    mutate(Domain = factor(Domain, levels = levels_Domain),
           NewTax = factor(NewTax, levels = levels_NewTax)) %>%
    complete(Domain, NewTax, fill = list(Count = 0, Proportion = 0)) %>%
    filter(!is.na(NewTax)) %>%  # double-check: remove any NA NewTax after complete
    {.[.$Domain %in% levels_Domain[1:5], ]}

  count_ylabel <- paste0(count_type, " Count")
  count_title  <- paste0(count_type, " Count by Domain and NewTax")
  prop_ylabel  <- paste0(count_type, " Proportion")
  prop_title   <- paste0(count_type, " Proportion by Domain and NewTax")
}

# -------------------------- contig_Tax_summary plots ---------------------
if (is_multi_source) {
  n_sources <- length(actual_sources)
  ncol_facet <- ifelse(n_sources <= 3, n_sources, ifelse(n_sources <= 9, 3, 4))

  p1 <- ggplot(contig_Tax_summary,
               aes(x = Domain, y = Count, fill = NewTax)) +
    geom_col(position = position_dodge(preserve = "single"), alpha = opt$alpha) +
    geom_text(aes(label = format_big_num(Count)),
              position = position_dodge(width = 0.9),color = "black", vjust = -0.5, size=2, check_overlap = T) +
    get_fill_scale(opt$colorSet) +
    labs(title = count_title,
         x = "Domain", y = count_ylabel) +
    facet_wrap(~FacetLabel, scales = "free_y", ncol = ncol_facet) +
    theme_minimal(base_size = opt$fontSize) +
    theme(axis.ticks.length = unit(0.2, "cm"),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_rect(fill = "grey90", colour = "black"),
          strip.text = element_text(face = "bold"))

  wf <- max(1.4, ncol_facet * 0.8)
  hf <- max(1, ceiling(n_sources / ncol_facet) * 0.7)
  save_plot(p1, "contig_count", width_factor = wf, height_factor = hf)

  p2 <- ggplot(contig_Tax_summary,
               aes(x = Domain, y = Proportion, fill = NewTax)) +
    geom_col(position = position_dodge(preserve = "single"), alpha = opt$alpha) +
    geom_text(aes(label = sprintf("%.0f%%", Proportion * 100)),
              position = position_dodge(width = 0.9),color = "black", vjust = -0.5, size=2, check_overlap = T) +
    get_fill_scale(opt$colorSet) +
    labs(title = prop_title,
         x = "Domain", y = prop_ylabel) +
    facet_wrap(~FacetLabel, scales = "free_y", ncol = ncol_facet) +
    theme_minimal(base_size = opt$fontSize) +
    theme(axis.ticks.length = unit(0.2, "cm"),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_rect(fill = "grey90", colour = "black"),
          strip.text = element_text(face = "bold"))

  save_plot(p2, "contig_prop", width_factor = wf, height_factor = hf)

} else {
  p1 <- ggplot(contig_Tax_summary,
               aes(x = Domain, y = Count, fill = NewTax)) +
    geom_col(position = position_dodge(preserve = "single"), alpha = opt$alpha) +
    geom_text(aes(label = format_big_num(Count)),
              position = position_dodge(width = 0.9),color = "black", vjust = -0.5, size=2, check_overlap = T) +
    get_fill_scale(opt$colorSet) +
    labs(title = count_title,
         x = "Domain", y = count_ylabel) +
    theme_minimal(base_size = opt$fontSize) +
    theme(axis.ticks.length = unit(0.2, "cm"),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))

  save_plot(p1, "contig_count", width_factor = 1.4)

  p2 <- ggplot(contig_Tax_summary,
               aes(x = Domain, y = Proportion, fill = NewTax)) +
    geom_col(position = position_dodge(preserve = "single"), alpha = opt$alpha) +
    geom_text(aes(label = sprintf("%.0f%%", Proportion * 100)),
              position = position_dodge(width = 0.9),color = "black", vjust = -0.5, size=2, check_overlap = T) +
    get_fill_scale(opt$colorSet) +
    labs(title = prop_title,
         x = "Domain", y = prop_ylabel) +
    theme_minimal(base_size = opt$fontSize) +
    theme(axis.ticks.length = unit(0.2, "cm"),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))

  save_plot(p2, "contig_prop", width_factor = 1.4)
}

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

save_plot(p3, "ORFKO_count", width_factor = 1.2)

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

save_plot(p4, "ORFKO_prop", width_factor = 1.2)

# -------------------------- ORF_Tax_KO_summary plots -----------------------
levels_Domain_NewTax = vector()
for (x in levels_Domain[1:5]){
  for (y in levels_NewTax[1:7]){
    levels_Domain_NewTax = c(levels_Domain_NewTax, paste0(x,"_", y))
  }
}

ORF_Tax_KO_summary <- ORF_Tax_KO_summary %>%
  {.[.$Domain %in% levels_Domain[2:5],]} %>%
  {.[.$NewTax %in% levels_NewTax[1:7],]} %>%
  mutate(Domain = factor(Domain, levels = levels_Domain[2:5]),
         NewTax = factor(NewTax, levels = levels_NewTax[1:7]),
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

save_plot(p5, "ORFTaxKO_count",  width_factor = 1.35, height_factor = 1.2)

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

save_plot(p6, "ORFTaxKO_prop",  width_factor = 1.35, height_factor = 1.2)

message("All plots saved.")
