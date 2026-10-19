#' @Date: 2022-11-02 23:40:20
#' @LastEditors: LiuyangLi
#' @LastEditTime: 2024-03-20 20:50:08
#' @FilePath: /liliuyang/R/Rscript.MASH/202308-v7/M1S5.Heatmap.cmd.R
#' @Description: 删除分组信息

suppressMessages(library(optparse))
suppressMessages(library(pheatmap))
suppressMessages(library(tidyverse))

# setwd("E:/deepsea/R/metadata.20221101/R.5group")
rm(list = ls(all.names = TRUE))

option_list <- list(
  make_option(c("--user_Sample_in_DIR"), type="character", default=NULL, help="user_Sample directory"),
  make_option(c("--out_DIR"), type="character", default=NULL, help="output directory for Heatmap figure"),
  make_option(c("--verbose"), action="store_true", default=FALSE, help="print messages")
)
opt <- parse_args(OptionParser(option_list=option_list))

user_Sample_in_DIR <- opt$user_Sample_in_DIR
out_DIR <- opt$out_DIR
verbose <- opt$verbose

if(! dir.exists(out_DIR)){
  dir.create(out_DIR)
}

draw_heatmap <- function(sub_Sample_KO) {
  min_val <- min(sub_Sample_KO, na.rm = TRUE)
  max_val <- max(sub_Sample_KO, na.rm = TRUE)

  # FIX: 单样品时无法聚类，强制关闭列聚类
  do_cluster_cols <- ncol(sub_Sample_KO) >= 2

  # 动态生成 colors，使0精确对应白色
  n_colors <- 100  # 减少颜色数量，避免color key绘制问题

  if (min_val < 0 && max_val > 0) {
    # 数据跨越0：计算0在色板中的位置，精确分配颜色
    zero_pos <- (0 - min_val) / (max_val - min_val) * n_colors
    n_neg <- floor(zero_pos)   # 负数段颜色数
    n_pos <- n_colors - n_neg  # 正数段颜色数

    colors <- c(colorRampPalette(c("#4575B4", "white"))(n_neg),
                colorRampPalette(c("white", "#D73027"))(n_pos))

  } else if (min_val == max_val) {
    colors <- rep("white", n_colors)
  } else if (min_val >= 0) {
    colors <- colorRampPalette(c("white", "#D73027"))(n_colors)
  } else {
    colors <- colorRampPalette(c("#4575B4", "white"))(n_colors)
  }

  # 不指定breaks，让pheatmap自动生成
  # 自动生成breaks: seq(min_val, max_val, length.out = n_colors + 1)
  # 0会落在白色边界上

  legend_breaks <- pretty(c(min_val, max_val), n = 5)

  tryCatch({
    # FIX: 当 min_val == max_val 时，pheatmap 的自动 breaks 会全为 0，导致 cut.default 报错
    # 解决方案：手动构造不重复的 breaks，legend 只显示一个值
    if (min_val == max_val) {
      # 构造一个对称的微小区间，使 breaks 不重复
      eps <- 1e-10
      breaks <- seq(min_val - eps, max_val + eps, length.out = n_colors + 1)
      # FIX: legend_breaks 不能重复，只显示一个值 0
      legend_breaks_single <- min_val
      legend_labels_single <- as.character(min_val)

      p = sub_Sample_KO %>%
        {pheatmap(.,
                  method = 'average',
                  na_col='grey',
                  color = colors,
                  breaks = breaks,  # 手动指定不重复的 breaks
                  legend_breaks = legend_breaks_single,
                  legend_labels = legend_labels_single,
                  cluster_rows = F,
                  cluster_cols = do_cluster_cols,
                  show_rownames = T,show_colnames = T,
                  cellwidth = 6,cellheight = 6,
                  display_numbers = F,
                  fontsize_row = 5,
                  fontsize_col	= 5,
                  angle_col = 45,
                  silent = TRUE
        )}
    } else {
      p = sub_Sample_KO %>%
        {pheatmap(.,
                  method = 'average',
                  na_col='grey',
                  color = colors,
                  legend_breaks = legend_breaks,
                  cluster_rows = F,
                  cluster_cols = do_cluster_cols,
                  show_rownames = T,show_colnames = T,
                  cellwidth = 6,cellheight = 6,
                  display_numbers = F,
                  fontsize_row = 5,
                  fontsize_col	= 5,
                  angle_col = 45,
                  silent = TRUE
        )}
    }
    return(p)
  }, error = function(e) {
    message("Error in pheatmap: ", e)
    colors <- rep("white", times = 100)
    p <- pheatmap(sub_Sample_KO,
                  na_col = 'grey',
                  scale = 'none',
                  color = colors,
                  cluster_rows = F,
                  cluster_cols = F,
                  show_rownames = T,show_colnames = T,
                  cellwidth = 6,cellheight = 6,
                  display_numbers = F,
                  fontsize_row = 5,
                  fontsize_col = 5,
                  angle_col = 45,
                  silent = TRUE
    )
    return(p)
  })
}

user_file <- dir(user_Sample_in_DIR)
user_file <- user_file[grepl(".tsv", user_file)]
if (length(user_file) == 0) {
  message("No module tables (.tsv) found in ", user_Sample_in_DIR,
          "; skipping heatmap generation (upstream mclink step was skipped or produced no output).")
  quit(save = "no", status = 0)
}
for (file in user_file) {
  if (grepl(".tsv", file)){
    file_name = gsub(".tsv", "", file)
    print(file)
    sub_Sample_KO <- read.csv(paste(user_Sample_in_DIR, file, sep = '/'),
                                quote='',sep = '\t',header = T, row.names = 1, check.names = FALSE)

    n_na <- sum(is.na(as.matrix(sub_Sample_KO)))
    if (n_na > 0) {
      message(paste0("WARNING: ", file, " contains ", n_na, " NA value(s); plotted as grey cells, please check the upstream data"))
    }

    num_sample = length(unique(colnames(sub_Sample_KO)))
    print(paste0("num of sample: ", num_sample))
    figure_width=7+0.12*num_sample

    num_KO = length(unique(rownames(sub_Sample_KO)))
    print(paste0("num of KO: ", num_KO))
    figure_height=5+0.12*num_KO

    p =  draw_heatmap(sub_Sample_KO)
    ggsave(paste0(file_name, ".pdf"),
                  path = out_DIR, p, limitsize = FALSE,
                  width = figure_width, height = figure_height)
  }
}
print("M1S5.Heatmap.cmd.R finished!")
