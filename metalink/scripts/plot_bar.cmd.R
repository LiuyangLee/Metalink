###
#' @Date: 2022-11-02 23:40:20
#' @LastEditors: LiuyangLi
#' @LastEditTime: 2026-07-18
#' @FilePath: plot_bar.cmd.R
#' @Description: Stacked barplot for Metalink wide tables (Taxa / KO / Module)
###
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(plotly))
suppressMessages(library(htmlwidgets))
suppressMessages(library(reshape2))
rm(list = ls(all.names = TRUE))

option_list = list(
  make_option("--Input", type = "character", default = NULL,
              help = "input wide TSV file, e.g., Taxa_Sample.Taxa_ReadCount.wide.tsv (Taxa_Sample.Taxa_SumCoverage.wide.tsv for m/pmr modes), KO_Sample.KO_SumCoverage.wide.tsv or mc_abundance_filtered_by_completeness.tsv"),
  make_option("--Type", type = "character", default = NULL,
              help = "input table type: Taxa (Taxa_Sample.Taxa_ReadCount.wide.tsv, or Taxa_Sample.Taxa_SumCoverage.wide.tsv for m/pmr modes, plot D/P/C/O/F/G/S levels), KO (KO_Sample.KO_SumCoverage.wide.tsv, plot KO directly), Module (mc_abundance_filtered_by_completeness.tsv, plot Level_2 and Module_Name)"),
  make_option("--TopN", type = "integer", default = 29,
              help = "number of top categories to keep, the rest are merged into Others, maximum 29 [default %default]"),
  make_option("--Output", type = "character", default = NULL,
              help = "output Dir for Barplot figure, e.g., Barplot_Figure"),
  make_option("--Alpha", type = "double", default = 0.7,
              help = "bar fill transparency, 0=fully transparent, 1=fully opaque [default %default]"),
  make_option("--verbose", action = "store_true", default = FALSE,
              help = "print messages")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

Input  = opt$Input
Type   = opt$Type
TopN   = opt$TopN
Output = opt$Output
Alpha  = opt$Alpha
if (is.null(Input) | is.null(Type) | is.null(Output)){
  print_help(opt_parser)
  stop("Input, Type and Output are required.", call. = FALSE)
}

TopN = as.integer(TopN)
if (TopN > 29){
  message("TopN is capped at 29 (color theme holds 30 colors: 29 categories + Others).")
  TopN = 29
}
if (TopN < 1) stop("TopN must be a positive integer.")
if (! Type %in% c('Taxa', 'KO', 'Module')) stop("Type must be one of: Taxa, KO, Module")
Alpha = as.numeric(Alpha)
if (Alpha < 0 || Alpha > 1) stop("Alpha must be between 0 and 1.")

################################## for ##################################
color_theme = c(
  "#0288D1","#FF9800","#F44336","#26A69A","#FFCA28","#26C6DA",
  "#FFEE58","#009688","#8BC34A","#AB47BC","#CDDC39","#FFC107",
  "#E91E63","#9CCC65","#795548","#9C27B0","#3F51B5","#42A5F5",
  "#EF5350","#00BCD4","#66BB6A","#FF5722","#E6EE9C","#C5CCFF",
  "#FFEB3B","#D4E157","#673AB7","#4CAF50","#EC407A","#9E9E9E"
)
if(! dir.exists(Output)){
  dir.create(Output)
}
# ----- read the wide input table -----
raw_table = read.delim(Input, sep='\t', header = TRUE, check.names = FALSE,
                       stringsAsFactors = FALSE, comment.char = "")
names(raw_table) = sub("\r$", "", names(raw_table)) # tolerate CRLF line endings

# ----- type-specific settings -----
if (Type == 'Taxa'){
  levels_to_plot = c('D', 'P', 'C', 'O', 'F', 'G', 'S')
  meta_cols = c('Taxonomy', 'Level')
} else if (Type == 'KO'){
  levels_to_plot = c('KO')
  meta_cols = c('KO')
} else { # Module
  levels_to_plot = c('Level_2', 'Level_3', 'Module_Name')
  meta_cols = c('Module_Entry', 'Level_2', 'Level_3', 'Module_Name', 'Definition')
}
sample_cols = setdiff(names(raw_table), meta_cols)

n_na = sum(is.na(as.matrix(raw_table[, sample_cols, drop = FALSE])))
if (n_na > 0){
  message('WARNING: ', Input, ' contains ', n_na, ' NA/empty value(s) in sample columns; they are plotted as missing, please check the upstream data')
}

# ----- legend columns per level -----
legend_col_map = list(
  'Taxa'   = c('D' = 6, 'P' = 5, 'C' = 2, 'O' = 1, 'F' = 1, 'G' = 1, 'S' = 1),
  'KO'     = c('KO' = 6),
  'Module' = c('Level_2' = 3, 'Level_3' = 3, 'Module_Name' = 1)
)

for (level in levels_to_plot) {
  print(paste0('Type: ', Type))
  print(paste0('level: ', level))

  legend_col = unname(legend_col_map[[Type]][[level]])

  # ----- build the label vector for this level -----
  if (Type == 'Taxa'){
    sub_table = raw_table[raw_table[['Level']] == level, , drop = FALSE]
    #label_vec = sub('.*;', '', sub_table[['Taxonomy']]) # keep the last rank, e.g. g__WIAN01
    label_vec = sub_table[['Taxonomy']]
  } else if (Type == 'KO'){
    sub_table = raw_table
    label_vec = sub_table[['KO']]
  } else { # Module, level is 'Level_2' or 'Module_Name'
    sub_table = raw_table
    label_vec = sub_table[[level]]
  }

  # ----- wide to long -----
  long_table = data.frame(Taxonomy = label_vec, sub_table[, sample_cols, drop = FALSE],
                          check.names = FALSE, stringsAsFactors = FALSE) %>%
    pivot_longer(cols = all_of(sample_cols), names_to = 'sample_name', values_to = 'RawValue') %>%
    mutate(RawValue = as.numeric(RawValue))

  # ----- relative abundance (%) within each sample -----
  long_table = long_table %>%
    group_by(sample_name) %>%
    mutate(Abudance = if (sum(RawValue, na.rm = TRUE) == 0) 0
                      else RawValue / sum(RawValue, na.rm = TRUE) * 100) %>%
    ungroup()

  # ----- top N categories + Others -----
  cat_total = sort(tapply(long_table$Abudance, long_table$Taxonomy, sum), decreasing = TRUE)
  zero_cats = names(cat_total)[cat_total <= 0] # categories absent from all samples
  if (length(zero_cats) > 0){
    message('Dropped all-zero categories: ', paste(zero_cats, collapse = ', '))
    long_table = long_table[!(long_table$Taxonomy %in% zero_cats), , drop = FALSE]
    cat_total = cat_total[cat_total > 0]
  }
  if (nrow(long_table) == 0){
    message('Skipping level ', level, ': all categories are zero.')
    next
  }
  top_cats = names(cat_total)[seq_len(min(TopN, length(cat_total)))]
  has_others = length(cat_total) > TopN
  long_table$Taxonomy[!(long_table$Taxonomy %in% top_cats)] = 'Others'
  data = aggregate(Abudance ~ sample_name + Taxonomy, data = long_table, FUN = sum)
  levels_taxonomy = if (has_others) c(top_cats, 'Others') else top_cats

  # ----- named fill colors, Others forced to grey -----
  fill_colors = setNames(color_theme[seq_along(levels_taxonomy)], levels_taxonomy)
  if ('Others' %in% levels_taxonomy) fill_colors[['Others']] = '#9E9E9E'

  # ----- sample order (from input table column order) -----
  levels_sample_names = sample_cols[sample_cols %in% unique(data$sample_name)]

  data = data %>%
    {mutate(., Taxonomy = factor(.$Taxonomy, levels = (levels_taxonomy)))} %>%
    {mutate(., sample_name = factor((.$sample_name), levels = (levels_sample_names)))}

  # ----- export the plotting table (wide format, relative abundance) -----
  write.table(data %>%
              {pivot_wider(select(., sample_name, Taxonomy, Abudance),
                           names_from = sample_name, values_from = Abudance)} %>%
              {arrange(., Taxonomy)},
            file = paste0(Output, "/Barplot.", Type, ".", level, ".tsv"),
            row.names = FALSE, quote = FALSE, sep = "\t")

  num_sample = length(unique(levels_sample_names))
  print(paste0("num of sample: ", num_sample))

  figure_width = 24 + 0.1*num_sample

  figure_vjust = 1.5 + 0.005*num_sample
  if (figure_vjust > 2){
    figure_vjust = 2
  }
  figure_text_size = (15 - length(levels_sample_names)/20)
  if (figure_text_size < 5){
    figure_text_size = 5
  }

  figure_text_size = figure_text_size + 5
  figure_vjust = 1.2


  n_taxa  <- nlevels(data$Taxonomy)
  n_row   <- ceiling(n_taxa / legend_col)

  legend_h_in <- (n_row * 15 + 2) / 48

  panel_h_in  <- 9
  top_margin  <- 0.5
  bottom_margin <- 0.5
  total_h_in  <- legend_h_in + panel_h_in + top_margin + bottom_margin

  p = data %>%
    {ggplot(.,aes(x=sample_name,y=Abudance,fill=Taxonomy))+
        geom_bar(stat = "identity", width = 0.7, position = position_stack(reverse = TRUE), alpha = Alpha) +
        scale_fill_manual(values = fill_colors)+
        theme(plot.title = element_text(hjust=0.5),
              title=element_text(size=20),
              axis.text.x = element_text(angle = 90, hjust = 1, vjust = figure_vjust, size= figure_text_size),
              axis.text.y = element_text(hjust = 1, size= figure_text_size),
              axis.ticks.x = element_blank(),
              legend.title=element_blank(),
              legend.position = 'bottom',
              legend.text = element_text(size = 15),
              legend.margin = margin(0,0,0,0),
              legend.spacing = unit(0, "pt"),
              strip.text = element_text(size = 24),
              panel.background = element_rect(fill=NA),
              panel.grid = element_blank()
              ) +
        labs(x=NULL,y='Relative Abundance (%)') +
        guides(fill=guide_legend(ncol = legend_col, byrow = T, reverse = F, label.position = 'right'))
    }
  #saveWidget(ggplotly(p), file = paste0(Output,"/Barplot.", Type, ".", level, ".html"))
  ggsave(paste0("Barplot.", Type, ".", level, ".pdf"), p,path = Output,width = figure_width, height = total_h_in, limitsize = F)
  print(paste("Barplot completed:", Type, level))

}


# ----- Module mode: Energy metabolism Level_3 breakdown -----
if (Type == 'Module') {
  energy_modules = raw_table[raw_table[['Level_2']] == 'Energy metabolism', , drop = FALSE]
  if (nrow(energy_modules) > 0) {
    energy_level3s = unique(energy_modules[['Level_3']])
    for (lvl3 in energy_level3s) {
      print(paste0('Processing Energy metabolism Level_3: ', lvl3))

      sub_energy = energy_modules[energy_modules[['Level_3']] == lvl3, , drop = FALSE]
      label_vec = sub_energy[['Module_Name']]

      # wide to long
      long_table = data.frame(Taxonomy = label_vec, sub_energy[, sample_cols, drop = FALSE],
                              check.names = FALSE, stringsAsFactors = FALSE) %>%
        pivot_longer(cols = all_of(sample_cols), names_to = 'sample_name', values_to = 'RawValue') %>%
        mutate(RawValue = as.numeric(RawValue))

      # relative abundance (%)
      long_table = long_table %>%
        group_by(sample_name) %>%
        mutate(Abudance = if (sum(RawValue, na.rm = TRUE) == 0) 0
                          else RawValue / sum(RawValue, na.rm = TRUE) * 100) %>%
        ungroup()

      # top N + Others
      cat_total = sort(tapply(long_table$Abudance, long_table$Taxonomy, sum), decreasing = TRUE)
      zero_cats = names(cat_total)[cat_total <= 0]
      if (length(zero_cats) > 0){
        message('Dropped all-zero categories: ', paste(zero_cats, collapse = ', '))
        long_table = long_table[!(long_table$Taxonomy %in% zero_cats), , drop = FALSE]
        cat_total = cat_total[cat_total > 0]
      }
      if (nrow(long_table) == 0){
        message('Skipping Energy metabolism Level_3 ', lvl3, ': all categories are zero.')
        next
      }
      top_cats = names(cat_total)[seq_len(min(TopN, length(cat_total)))]
      has_others = length(cat_total) > TopN
      long_table$Taxonomy[!(long_table$Taxonomy %in% top_cats)] = 'Others'
      data = aggregate(Abudance ~ sample_name + Taxonomy, data = long_table, FUN = sum)
      levels_taxonomy = if (has_others) c(top_cats, 'Others') else top_cats

      # colors
      fill_colors = setNames(color_theme[seq_along(levels_taxonomy)], levels_taxonomy)
      if ('Others' %in% levels_taxonomy) fill_colors[['Others']] = '#9E9E9E'

      # sample order
      levels_sample_names = sample_cols[sample_cols %in% unique(data$sample_name)]

      data = data %>%
        {mutate(., Taxonomy = factor(.$Taxonomy, levels = (levels_taxonomy)))} %>%
        {mutate(., sample_name = factor((.$sample_name), levels = (levels_sample_names)))}

      # export CSV
      write.table(data %>%
                  {pivot_wider(select(., sample_name, Taxonomy, Abudance),
                               names_from = sample_name, values_from = Abudance)} %>%
                  {arrange(., Taxonomy)},
                file = paste0(Output, "/Barplot.Module.Energy.", lvl3, ".tsv"),
                row.names = FALSE, quote = FALSE, sep = "\t")

      num_sample = length(unique(levels_sample_names))
      figure_width = 24 + 0.1*num_sample

      figure_vjust = 1.5 + 0.005*num_sample
      if (figure_vjust > 2) figure_vjust = 2
      figure_text_size = (15 - length(levels_sample_names)/20)
      if (figure_text_size < 5) figure_text_size = 5
      figure_text_size = figure_text_size + 5
      figure_vjust = 1.2

      n_taxa  <- nlevels(data$Taxonomy)
      n_row   <- ceiling(n_taxa / 1)
      legend_h_in <- (n_row * 15 + 2) / 48
      panel_h_in  <- 9
      top_margin  <- 0.5
      bottom_margin <- 0.5
      total_h_in  <- legend_h_in + panel_h_in + top_margin + bottom_margin

      p = data %>%
        {ggplot(.,aes(x=sample_name,y=Abudance,fill=Taxonomy))+
            geom_bar(stat = "identity", width = 0.7, position = position_stack(reverse = TRUE), alpha = Alpha) +
            scale_fill_manual(values = fill_colors)+
            theme(plot.title = element_text(hjust=0.5),
                  title=element_text(size=20),
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = figure_vjust, size= figure_text_size),
                  axis.text.y = element_text(hjust = 1, size= figure_text_size),
                  axis.ticks.x = element_blank(),
                  legend.title=element_blank(),
                  legend.position = 'bottom',
                  legend.text = element_text(size = 15),
                  legend.margin = margin(0,0,0,0),
                  legend.spacing = unit(0, "pt"),
                  strip.text = element_text(size = 24),
                  panel.background = element_rect(fill=NA),
                  panel.grid = element_blank()
                  ) +
            labs(x=NULL, y='Relative Abundance (%)',
                 title = paste('Energy metabolism:', lvl3)) +
            guides(fill=guide_legend(ncol = 1, byrow = T, reverse = F, label.position = 'right'))
        }

      ggsave(paste0("Barplot.Module.Energy.", lvl3, ".pdf"), p,
             path = Output, width = figure_width, height = total_h_in, limitsize = F)
      print(paste("Energy metabolism barplot completed:", lvl3))
    }
  }
}

print("plot_bar.cmd.R finished!")
