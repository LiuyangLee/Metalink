#' @Date: 2022-11-02 23:40:20
#' @LastEditors: LiuyangLi
#' @LastEditTime: 2024-03-20 20:52:01
#' @FilePath: /liliuyang/R/Rscript.MASH/202308-v7/M1S4.Heatmap.cmd.R
#' @Description: 删除分组信息

suppressMessages(library(optparse))
suppressMessages(library(pheatmap))
suppressMessages(library(tidyverse))

# setwd("E:/deepsea/R/metadata.20221101/R.5group")
rm(list = ls(all.names = TRUE))

option_list <- list(
  make_option(c("--user_Sample_in_DIR"), type="character", default=NULL,
              help="user_Sample directory"),
  make_option(c("--out_DIR"), type="character", default=NULL,
              help="output directory for Heatmap figure"),
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help="print messages")
)
opt <- parse_args(OptionParser(option_list=option_list))
user_Sample_in_DIR <- opt$user_Sample_in_DIR
out_DIR <- opt$out_DIR
verbose <- opt$verbose

draw_heatmap <- function(sub_Sample_KO) {
    #  'sub_Sample_KO' 是热图数据
  # FIX: 单样品时无法聚类，强制关闭列聚类
  do_cluster_cols <- ncol(sub_Sample_KO) >= 2
  tryCatch({
    p <- pheatmap(sub_Sample_KO,
                  method = 'average',
                  na_col = 'grey',
                  cluster_rows = F,
                  cluster_cols = do_cluster_cols,###是否对行列聚类
                  show_rownames = T,show_colnames = T,###是否显示行名列名
                  cellwidth = 6,cellheight = 6,###固定每个cell的宽和高
                  display_numbers = F,####是否显示每个单元格的值
                  fontsize_row = 5,
                  fontsize_col = 5,
                  angle_col = 45,
                  silent = TRUE
    )
    return(p)
  }, error = function(e) {
    message("Error in pheatmap: ", e)
    min_val <- min(sub_Sample_KO, na.rm = TRUE)
    max_val <- max(sub_Sample_KO, na.rm = TRUE)

    # FIX: 当所有值相同时（如全为0），构造不重复的 breaks，
    # 颜色全蓝（与正常图0值对应的颜色一致），legend 只显示一个值 0
    if (min_val == max_val) {
      # pheatmap 要求 breaks 长度 = colors 长度 + 1，且严格递增
      eps <- 1e-10
      breaks <- seq(min_val - eps, max_val + eps, length.out = 101)
      # 全蓝色，与正常热图 color key 最底端颜色一致
      colors <- rep("#4575B4", times = 100)
      # legend 只显示一个刻度：0（或实际值），无梯度
      legend_breaks <- min_val
      legend_labels <- as.character(min_val)
    } else {
      breaks <- NULL
      colors <- NULL
      legend_breaks <- NULL
      legend_labels <- NULL
    }

    p <- pheatmap(sub_Sample_KO,
                  na_col = 'grey',
                  scale = 'none',
                  color = colors,
                  breaks = breaks,
                  legend_breaks = legend_breaks,
                  legend_labels = legend_labels,
                  cluster_rows = F,
                  cluster_cols = F,###当出现错误时，不进行列聚类
                  show_rownames = T,show_colnames = T,###是否显示行名列名
                  cellwidth = 6,cellheight = 6,###固定每个cell的宽和高
                  display_numbers = F,####是否显示每个单元格的值
                  fontsize_row = 5,
                  fontsize_col = 5,
                  angle_col = 45,
                  silent = TRUE
    )
    return(p)
  })
}


if(! dir.exists(out_DIR)){
  dir.create(out_DIR)
}

user_file <- dir(user_Sample_in_DIR)
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
    # 按照线性关系width=3+0.5*num_group
    figure_width=7+0.12*num_sample

    num_KO = length(unique(rownames(sub_Sample_KO)))
    print(paste0("num of KO: ", num_KO))
    # 按照线性关系width=3+0.5*num_group
    figure_height=5+0.12*num_KO
    ## 热图
    p =  draw_heatmap(sub_Sample_KO)
    ggsave(paste0(file_name, ".pdf"),
                  path = out_DIR, p, limitsize = FALSE,
                  width = figure_width, height = figure_height)
  }
}

print("M1S4.Heatmap.cmd.R finished!")
