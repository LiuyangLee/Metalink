#!/usr/bin/env Rscript
#' @Date: 2026-07-16
#' @Description: mclink.M1S5.Heatmap.html.cmd.R (Optimized for large matrices)

`%+%` <- function(a, b) paste0(a, b)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(pheatmap)
})

# ==================== Arguments ====================
option_list <- list(
  make_option(c("--mc_completeness_DIR"), type = "character", default = NULL,
              help = "Directory containing module completeness/abundance TSV files", metavar = "DIR"),
  make_option(c("--mc_detected_KOs"), type = "character", default = NULL,
              help = "Path to detected KOs TSV file", metavar = "FILE"),
  make_option(c("--out_DIR"), type = "character", default = NULL,
              help = "Output directory for HTML files", metavar = "DIR")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$mc_completeness_DIR) || is.null(opt$mc_detected_KOs) || is.null(opt$out_DIR)) {
  print_help(opt_parser)
  stop("--mc_completeness_DIR, --mc_detected_KOs, and --out_DIR are required!", call. = FALSE)
}

# ==================== Settings ====================
CELL_WIDTH <- 24
CELL_HEIGHT <- 24
FONT_SIZE <- 9
DENDROGRAM_HEIGHT <- 80
COLORBAR_MARGIN <- 60
X_AXIS_GAP <- 10

# ==================== Color Function (M1S5: #4575B4 -> white -> #D73027) ====================
STOPS_M1S5 <- matrix(c(
  0.000,  69, 117, 180,
  0.125, 114, 150, 198,
  0.250, 163, 187, 218,
  0.375, 208, 220, 236,
  0.500, 255, 253, 253,
  0.625, 245, 203, 200,
  0.750, 235, 153, 148,
  0.875, 225,  98,  91,
  1.000, 215,  48,  39
), ncol = 4, byrow = TRUE)

val_to_color <- function(v, zmin, zmax, zmid) {
  if (is.na(v)) return("#CCCCCC")
  if (zmax == zmin || v == zmid) {
    t <- 0.5
  } else if (v < zmid) {
    t <- (v - zmin) / (zmid - zmin) * 0.5
  } else {
    t <- 0.5 + (v - zmid) / (zmax - zmid) * 0.5
  }
  t <- max(0, min(1, t))
  r <- round(approx(STOPS_M1S5[,1], STOPS_M1S5[,2], xout = t, rule = 2)$y)
  g <- round(approx(STOPS_M1S5[,1], STOPS_M1S5[,3], xout = t, rule = 2)$y)
  b <- round(approx(STOPS_M1S5[,1], STOPS_M1S5[,4], xout = t, rule = 2)$y)
  return(sprintf("#%02X%02X%02X", r, g, b))
}

val_to_color_vec <- function(v_vec, zmin, zmax, zmid) {
  t <- rep(0.5, length(v_vec))
  not_eq <- !(zmax == zmin | v_vec == zmid)
  lt_zmid <- not_eq & (v_vec < zmid)
  gt_zmid <- not_eq & (v_vec > zmid)
  t[lt_zmid] <- (v_vec[lt_zmid] - zmin) / (zmid - zmin) * 0.5
  t[gt_zmid] <- 0.5 + (v_vec[gt_zmid] - zmid) / (zmax - zmid) * 0.5
  t <- pmax(0, pmin(1, t))
  r <- round(approx(STOPS_M1S5[,1], STOPS_M1S5[,2], xout = t, rule = 2)$y)
  g <- round(approx(STOPS_M1S5[,1], STOPS_M1S5[,3], xout = t, rule = 2)$y)
  b <- round(approx(STOPS_M1S5[,1], STOPS_M1S5[,4], xout = t, rule = 2)$y)
  res <- sprintf("#%02X%02X%02X", r, g, b)
  res[is.na(v_vec)] <- "#CCCCCC"
  res
}

# ==================== Load detected KOs (shared across all files) ====================
message("Loading: ", opt$mc_detected_KOs)
mc_kos <- fread(opt$mc_detected_KOs, sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Build module info lookup (shared)
module_info <- list()

kos_cols <- names(mc_kos)
meta_cols <- intersect(c("Module_Entry", "Level_2", "Level_3", "Module_Name", "Definition"), kos_cols)
kos_sample_cols <- kos_cols[!kos_cols %in% meta_cols]

for (i in seq_len(nrow(mc_kos))) {
  mod_name <- as.character(mc_kos$Module_Name[i])
  all_kos <- character(0)
  for (sc in kos_sample_cols) {
    if (sc %in% names(mc_kos)) {
      v <- as.character(mc_kos[[sc]][i])
      if (!is.na(v) && nchar(v) > 0 && v != "0" && v != "-") {
        k <- unlist(strsplit(v, "\\s+", perl = TRUE))
        k <- k[nchar(k) > 0]
        all_kos <- c(all_kos, k)
      }
    }
  }
  dkos <- paste(unique(all_kos), collapse = " ")
  if (dkos == "") dkos <- "None"
  def <- as.character(mc_kos$Definition[i])
  if (is.na(def)) def <- "N/A"
  entry <- as.character(mc_kos$Module_Entry[i])
  if (is.na(entry)) entry <- "N/A"
  module_info[[mod_name]] <- list(
    entry = entry,
    definition = def,
    detected_kos = dkos
  )
}

message(sprintf("Module info lookup: %d entries", length(module_info)))

# ==================== Generate Column Dendrogram SVG ====================
generate_dendrogram_svg <- function(hclust_obj, width, height) {
  n <- length(hclust_obj$labels)
  if (n <= 1) return("")
  merge <- hclust_obj$merge
  heights <- hclust_obj$height
  order_vec <- hclust_obj$order
  display_pos <- rep(0, n)
  for (disp_idx in 1:n) {
    orig_idx <- order_vec[disp_idx]
    display_pos[orig_idx] <- disp_idx
  }
  x_pos <- rep(0, n + nrow(merge))
  y_pos <- rep(0, n + nrow(merge))
  for (i in 1:n) {
    x_pos[i] <- display_pos[i] - 0.5
    y_pos[i] <- 0
  }
  for (i in 1:nrow(merge)) {
    left <- merge[i, 1]
    right <- merge[i, 2]
    left_idx <- ifelse(left < 0, -left, n + left)
    right_idx <- ifelse(right < 0, -right, n + right)
    x_pos[n + i] <- (x_pos[left_idx] + x_pos[right_idx]) / 2
    y_pos[n + i] <- heights[i]
  }
  max_h <- max(heights)
  if (max_h == 0) max_h <- 1
  x_scale <- width / n
  y_scale <- height / max_h
  paths <- character(0)
  for (i in 1:nrow(merge)) {
    left <- merge[i, 1]
    right <- merge[i, 2]
    left_idx <- ifelse(left < 0, -left, n + left)
    right_idx <- ifelse(right < 0, -right, n + right)
    x1 <- x_pos[left_idx] * x_scale
    y1 <- height - y_pos[left_idx] * y_scale
    x2 <- x_pos[n + i] * x_scale
    y2 <- height - y_pos[n + i] * y_scale
    x3 <- x_pos[right_idx] * x_scale
    y3 <- height - y_pos[right_idx] * y_scale
    path <- sprintf("M%.2f,%.2f L%.2f,%.2f L%.2f,%.2f L%.2f,%.2f L%.2f,%.2f",
                    x1, y1, x1, y2, x2, y2, x3, y2, x3, y3)
    paths <- c(paths, path)
  }
  return(paste(paths, collapse = " "))
}

# ==================== Process each file in mc_completeness_DIR ====================
if (!dir.exists(opt$out_DIR)) {
  dir.create(opt$out_DIR, recursive = TRUE)
}

user_files <- dir(opt$mc_completeness_DIR)
user_files <- user_files[grepl(".tsv$", user_files)]
if (length(user_files) == 0) {
  message("No module tables (.tsv) found in ", opt$mc_completeness_DIR,
          "; skipping HTML heatmap generation (upstream mclink step was skipped or produced no output).")
  quit(save = "no", status = 0)
}

for (file in user_files) {
  if (!grepl(".tsv$", file)) {
    next
  }

  tryCatch({

  file_name <- gsub(".tsv$", "", file)
  file_path <- file.path(opt$mc_completeness_DIR, file)
  out_html_path <- file.path(opt$out_DIR, file_name %+% ".html")

  message("\n========================================")
  message("Processing: ", file)
  message("Output: ", out_html_path)

  # Load data
  sub_Sample_KO <- read.csv(file_path, quote = '', sep = '\t', header = TRUE,
                            row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

  nr <- nrow(sub_Sample_KO)
  nc <- ncol(sub_Sample_KO)
  message(sprintf("Matrix: %d modules x %d samples", nr, nc))

  if (nr == 0 || nc == 0) {
    message("Skipping empty matrix.")
    next
  }

  # FIX: Single sample skip
  if (nc == 1) {
    message("Single sample detected, skipping heatmap generation (use SingleSample script for single-sample heatmaps).")
    next
  }

  # Column Clustering via pheatmap
  message("Performing column clustering via pheatmap (method = average)...")

  min_val <- min(sub_Sample_KO, na.rm = TRUE)
  max_val <- max(sub_Sample_KO, na.rm = TRUE)

  # 数据完全相同时 pheatmap 的自动 breaks 会重复，导致 cut.default 报错
  is_constant <- !is.finite(min_val) || !is.finite(max_val) || (max_val - min_val < 1e-12)
  do_cluster_cols <- (ncol(sub_Sample_KO) >= 2) && !is_constant

  n_colors <- 100
  if (min_val < 0 && max_val > 0) {
    zero_pos <- (0 - min_val) / (max_val - min_val) * n_colors
    n_neg <- floor(zero_pos)
    n_pos <- n_colors - n_neg
    colors <- c(colorRampPalette(c("#4575B4", "white"))(n_neg),
                colorRampPalette(c("white", "#D73027"))(n_pos))
  } else if (min_val == max_val) {
    colors <- rep("white", n_colors)
  } else if (min_val >= 0) {
    colors <- colorRampPalette(c("white", "#D73027"))(n_colors)
  } else {
    colors <- colorRampPalette(c("#4575B4", "white"))(n_colors)
  }

  if (do_cluster_cols) {
    ph <- pheatmap(sub_Sample_KO,
                   method = 'average',
                   na_col = 'grey',
                   color = colors,
                   cluster_rows = FALSE,
                   cluster_cols = TRUE,
                   show_rownames = TRUE,
                   show_colnames = TRUE,
                   cellwidth = 6,
                   cellheight = 6,
                   display_numbers = FALSE,
                   fontsize_row = 5,
                   fontsize_col = 5,
                   angle_col = 45,
                   silent = TRUE
    )
    col_order <- ph$tree_col$order
    col_hclust <- ph$tree_col
  } else {
    col_order <- seq_len(ncol(sub_Sample_KO))
    col_hclust <- NULL
  }

  message(sprintf("Column order from pheatmap: %s", paste(col_order, collapse = ", ")))

  mat <- as.matrix(sub_Sample_KO)
  mat <- mat[, col_order, drop = FALSE]
  sample_cols <- colnames(sub_Sample_KO)[col_order]
  nc <- length(sample_cols)

  message(sprintf("Column names after pheatmap clustering: %s", paste(sample_cols, collapse = ", ")))

  # Generate dendrogram SVG
  dendrogram_svg <- generate_dendrogram_svg(col_hclust, nc * CELL_WIDTH, DENDROGRAM_HEIGHT)

  # [MODIFIED] Color Range: mirror the PDF version's palette logic
  zmin <- min(mat, na.rm = TRUE)
  zmax <- max(mat, na.rm = TRUE)

  if (zmin < 0 && zmax > 0) {
    zmid <- 0
  } else if (zmax <= 0) {
    zmid <- zmax
  } else {
    zmid <- zmin
  }

  message(sprintf("Color range: [%.3f, %.3f, %.3f]", zmin, zmid, zmax))

  # [NEW] Build colorbar gradient by sampling val_to_color so that
  # the colorbar always matches the cell colors exactly
  if (zmax > zmin) {
    cb_vals <- sort(unique(c(seq(zmin, zmax, length.out = 41), zmid)))
    cb_pcts <- (cb_vals - zmin) / (zmax - zmin) * 100
  } else {
    cb_vals <- zmin
    cb_pcts <- 0
  }
  cb_colors <- val_to_color_vec(cb_vals, zmin, zmax, zmid)
  cb_stops <- sprintf("%s %.2f%%", cb_colors, cb_pcts)
  colorbar_gradient <- paste(cb_stops, collapse = ", ")

  # [OPTIMIZED] Build Color Matrix using vectorized function
  color_mat <- matrix(val_to_color_vec(as.vector(mat), zmin, zmax, zmid), nrow = nr, ncol = nc)

  # Calculate Dimensions
  max_row_label_len <- max(nchar(rownames(mat)))
  ROW_LABEL_WIDTH <- max(400, round(max_row_label_len * 5.8) + 20)
  message(sprintf("Row label width: %d px (max length: %d)", ROW_LABEL_WIDTH, max_row_label_len))

  max_col_label_len <- max(nchar(colnames(mat)))
  X_AXIS_LABEL_FONT_SIZE <- 9
  X_AXIS_CHAR_WIDTH <- 6
  label_vertical_ext <- round(max_col_label_len * X_AXIS_CHAR_WIDTH * 0.707)
  HEADER_HEIGHT <- max(20, label_vertical_ext + 10)

  message(sprintf("Max column label length: %d chars", max_col_label_len))
  message(sprintf("X-axis gap: %d px, Header height: %d px", X_AXIS_GAP, HEADER_HEIGHT))

  COLORBAR_WIDTH <- 80
  COLORBAR_HEIGHT <- 200
  total_width <- COLORBAR_WIDTH + COLORBAR_MARGIN + nc * CELL_WIDTH + ROW_LABEL_WIDTH + 50
  message(sprintf("Canvas width: %d pixels", total_width))

  # Colorbar Ticks
  cb_breaks <- pretty(c(zmin, zmax), n = 5)
  cb_breaks <- cb_breaks[cb_breaks >= zmin - 1e-9 & cb_breaks <= zmax + 1e-9]
  if (length(cb_breaks) == 0) {
    cb_breaks <- c(zmin, zmax)
  }
  message(sprintf("Colorbar breaks (pretty): %s", paste(sprintf("%.2f", cb_breaks), collapse = ", ")))

  # [OPTIMIZED] Generate Tooltip Data with pre-allocation
  n_cells <- nr * nc
  tooltip_data <- character(n_cells)
  idx <- 0
  for (i in seq_len(nr)) {
    mod <- rownames(mat)[i]
    info <- module_info[[mod]]
    if (is.null(info)) info <- list(entry = "N/A", definition = "N/A", detected_kos = "N/A")

    for (j in seq_len(nc)) {
      samp <- colnames(mat)[j]
      v <- mat[i, j]
      vstr <- ifelse(is.na(v), "N/A", sprintf("%.4f", v))

      tip <- paste0(
        '<div class="tt-row"><div class="tt-label">Module Entry</div><div class="tt-sep">:</div><div class="tt-value">', info$entry, '</div></div>',
        '<div class="tt-row"><div class="tt-label">Module Name</div><div class="tt-sep">:</div><div class="tt-value">', mod, '</div></div>',
        '<div class="tt-row"><div class="tt-label">Module Definition</div><div class="tt-sep">:</div><div class="tt-value def-val">', info$definition, '</div></div>',
        '<div class="tt-row"><div class="tt-label">Detect KOs</div><div class="tt-sep">:</div><div class="tt-value">', info$detected_kos, '</div></div>',
        '<div class="tt-row"><div class="tt-label">Module Abundance</div><div class="tt-sep">:</div><div class="tt-value">', vstr, '</div></div>',
        '<div class="tt-row"><div class="tt-label">Sample/Group</div><div class="tt-sep">:</div><div class="tt-value">', samp, '</div></div>'
      )

      tip_js <- gsub('\\', '\\\\', tip, fixed = TRUE)
      tip_js <- gsub('"', '\\"', tip_js, fixed = TRUE)
      tip_js <- gsub('\n', '\\n', tip_js, fixed = TRUE)

      tooltip_data[idx + 1] <- '    "' %+% idx %+% '": "' %+% tip_js %+% '"'
      idx <- idx + 1
    }
  }

  tooltip_json <- paste(tooltip_data, collapse = ",\n")

  # [OPTIMIZED] Stream HTML to file instead of building a huge vector
  con <- file(out_html_path, "w", encoding = "UTF-8")

  w <- function(...) {
    cat(..., "\n", file = con, sep = "")
  }

  w('<!DOCTYPE html>')
  w('<html><head>')
  w('<meta charset="UTF-8">')
  w('<title>Module Heatmap</title>')
  w('<style>')

  cat('
* { box-sizing: border-box; margin: 0; padding: 0; }

body {
  font-family: Arial, sans-serif;
  background: #fff;
  padding: 20px;
  min-width: ', total_width, 'px;
}

.title {
  font-size: 16px;
  font-weight: bold;
  color: #333;
  margin-bottom: 10px;
  padding-left: ', (COLORBAR_WIDTH + COLORBAR_MARGIN), 'px;
}

/* Main container: colorbar | data | row labels */
.heatmap-wrapper {
  display: flex;
  flex-direction: row;
  align-items: flex-start;
}

/* Left: colorbar */
.colorbar-col {
  width: ', COLORBAR_WIDTH, 'px;
  flex-shrink: 0;
  padding-top: ', DENDROGRAM_HEIGHT, 'px;
  margin-right: ', COLORBAR_MARGIN, 'px;
}

.colorbar-area {
  position: relative;
  width: ', COLORBAR_WIDTH, 'px;
  height: ', (nr * CELL_HEIGHT), 'px;
}

.colorbar {
  position: absolute;
  right: 0;
  top: 0;
  width: 16px;
  height: ', COLORBAR_HEIGHT, 'px;
  background: linear-gradient(to top, ', colorbar_gradient, ');
  border: 1px solid #999;
}

.colorbar-labels {
  position: absolute;
  right: 22px;
  top: 0;
  width: 50px;
  height: ', COLORBAR_HEIGHT, 'px;
}

.colorbar-label {
  position: absolute;
  right: 0;
  font-size: 10px;
  color: #333;
  line-height: 1;
  transform: translateY(-50%);
  white-space: nowrap;
}

/* Middle: data area */
.data-area {
  display: flex;
  flex-direction: column;
  position: relative;
}

/* Dendrogram area */
.dendrogram-row {
  height: ', DENDROGRAM_HEIGHT, 'px;
}

.dendrogram-svg {
  width: ', (nc * CELL_WIDTH), 'px;
  height: ', DENDROGRAM_HEIGHT, 'px;
}

/* Data rows: black outer border */
.data-rows {
  display: flex;
  flex-direction: column;
  border: 2px solid #000;
}

/* Data row: no bottom border (cell borders handle it) */
.data-row {
  display: flex;
  flex-direction: row;
  height: ', CELL_HEIGHT, 'px;
}

/* Cell: gray inner border */
.cell {
  width: ', CELL_WIDTH, 'px;
  height: ', CELL_HEIGHT, 'px;
  flex-shrink: 0;
  cursor: crosshair;
  position: relative;
  border: 1px solid #999;
}

.cell-inner {
  width: 100%;
  height: 100%;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 10px;
  color: transparent;
  pointer-events: none;
}

.cell:hover .cell-inner {
  color: #333;
  text-shadow: 0 0 2px #fff, 0 0 3px #fff, 0 0 4px #fff;
}

/* Column headers: positioned below heatmap with gap */
.col-headers-row {
  position: absolute;
  bottom: -', X_AXIS_GAP, 'px;
  left: 0;
  height: ', HEADER_HEIGHT, 'px;
  display: flex;
  flex-direction: row;
  white-space: nowrap;
  overflow: visible;
  pointer-events: none;
}

.col-header {
  width: ', CELL_WIDTH, 'px;
  height: ', HEADER_HEIGHT, 'px;
  text-align: left;
  font-size: ', FONT_SIZE, 'px;
  color: #333;
  position: relative;
  flex-shrink: 0;
  overflow: visible;
}

/* X-axis labels: at bottom of header, extend left and down (away from heatmap) */
.col-header span {
  position: absolute;
  left: 50%;
  bottom: 0;
  transform: translateX(-100%) rotate(-45deg);
  transform-origin: right center;
  white-space: nowrap;
  font-size: ', X_AXIS_LABEL_FONT_SIZE, 'px;
}

/* Right: row labels */
.row-labels-col {
  width: ', ROW_LABEL_WIDTH, 'px;
  flex-shrink: 0;
  padding-top: ', DENDROGRAM_HEIGHT, 'px;
}

.row-label {
  height: ', CELL_HEIGHT, 'px;
  line-height: ', CELL_HEIGHT, 'px;
  font-size: ', FONT_SIZE, 'px;
  text-align: left;
  padding-left: 8px;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  color: #333;
  border-bottom: 1px solid #eee;
}

/* Tooltip */
.tooltip {
  display: none;
  position: fixed;
  background: white;
  border: 1.5px solid #CC3333;
  border-radius: 4px;
  padding: 10px 14px;
  font-size: 12px;
  line-height: 1.5;
  color: #333;
  z-index: 10000;
  max-width: 500px;
  box-shadow: 2px 2px 8px rgba(0,0,0,0.15);
  pointer-events: none;
  font-family: Arial, sans-serif;
}

.tt-row {
  display: grid;
  grid-template-columns: 130px 10px 1fr;
  gap: 2px;
  margin-bottom: 2px;
}

.tt-label {
  text-align: right;
  font-weight: bold;
  color: #555;
  white-space: nowrap;
}

.tt-sep {
  text-align: center;
  font-weight: bold;
  color: #555;
}

.tt-value {
  text-align: left;
  word-break: break-all;
  overflow-wrap: break-word;
}

.tt-value.def-val {
  max-width: 350px;
}
', file = con, sep = "")

  w('</style>')
  w('</head><body>')

  w('<div class="title">Module Completeness/Abundance Heatmap</div>')

  w('<div class="heatmap-wrapper">')

  # Colorbar (LEFT)
  w('<div class="colorbar-col">')
  w('<div class="colorbar-area">')
  w('<div class="colorbar"></div>')
  w('<div class="colorbar-labels">')

  for (b in cb_breaks) {
    if (zmax == zmin) {
      tpos <- 0.5
    } else {
      tpos <- (b - zmin) / (zmax - zmin)
    }
    tpos <- max(0, min(1, tpos))
    label_top <- COLORBAR_HEIGHT * (1 - tpos)

    if (abs(b - round(b)) < 1e-10) {
      bstr <- sprintf("%.0f", b)
    } else {
      bstr <- sprintf("%.2f", b)
    }

    w('<div class="colorbar-label" style="top:', label_top, 'px">', bstr, '</div>')
  }

  w('</div>')  # colorbar-labels
  w('</div>')  # colorbar-area
  w('</div>')  # colorbar-col

  # Data area (MIDDLE)
  w('<div class="data-area">')

  # Dendrogram
  w('<div class="dendrogram-row">')
  w('<svg class="dendrogram-svg" viewBox="0 0 ', (nc * CELL_WIDTH), ' ', DENDROGRAM_HEIGHT, '" preserveAspectRatio="none">')
  w('<path d="', dendrogram_svg, '" stroke="#333" stroke-width="1.5" fill="none"/>')
  w('</svg>')
  w('</div>')

  # Data rows
  w('<div class="data-rows">')
  idx <- 0
  for (i in seq_len(nr)) {
    w('<div class="data-row">')
    for (j in seq_len(nc)) {
      v <- mat[i, j]
      color <- color_mat[i, j]
      cell_text <- ifelse(is.na(v), "NA", sprintf("%.2f", v))
      w('<div class="cell" data-idx="', idx, '" style="background-color:', color, '">',
        '<div class="cell-inner">', cell_text, '</div>',
        '</div>')
      idx <- idx + 1
    }
    w('</div>')
  }
  w('</div>')  # data-rows

  # Column headers (BOTTOM, with gap)
  w('<div class="col-headers-row">')
  for (j in seq_len(nc)) {
    sname <- colnames(mat)[j]
    w('<div class="col-header" title="', sname, '"><span>', sname, '</span></div>')
  }
  w('</div>')

  w('</div>')  # data-area

  # Row labels (RIGHT)
  w('<div class="row-labels-col">')
  for (i in seq_len(nr)) {
    mod_name <- rownames(mat)[i]
    w('<div class="row-label" title="', mod_name, '">', mod_name, '</div>')
  }
  w('</div>')

  w('</div>')  # heatmap-wrapper

  # Tooltip
  w('<div id="tooltip" class="tooltip"></div>')

  # JavaScript
  w('<script>')
  w('(function() {')
  w('  var tooltipData = {')
  cat(tooltip_json, "\n", file = con, sep = "")
  w('  };')
  w('  var tooltip = document.getElementById("tooltip");')
  w('  var cells = document.querySelectorAll(".cell");')
  w('  cells.forEach(function(cell) {')
  w('    cell.addEventListener("mouseenter", function(e) {')
  w('      var idx = this.getAttribute("data-idx");')
  w('      var tip = tooltipData[idx];')
  w('      if (!tip) return;')
  w('      tooltip.innerHTML = tip;')
  w('      tooltip.style.display = "block";')
  w('      var rect = this.getBoundingClientRect();')
  w('      var tipRect = tooltip.getBoundingClientRect();')
  w('      var left = rect.right + 10;')
  w('      var top = rect.top;')
  w('      if (left + tipRect.width > window.innerWidth) {')
  w('        left = rect.left - tipRect.width - 10;')
  w('      }')
  w('      if (top + tipRect.height > window.innerHeight) {')
  w('        top = window.innerHeight - tipRect.height - 10;')
  w('      }')
  w('      if (top < 0) top = 10;')
  w('      tooltip.style.left = left + "px";')
  w('      tooltip.style.top = top + "px";')
  w('    });')
  w('    cell.addEventListener("mouseleave", function() {')
  w('      tooltip.style.display = "none";')
  w('    });')
  w('  });')
  w('})();')
  w('</script>')

  w('</body></html>')

  close(con)

  message("Done: ", out_html_path)
  message(sprintf("Rows: %d, Cols: %d, Width: %d px", nr, nc, total_width))
  }, error = function(e) {
    if (exists("con", inherits = FALSE) && inherits(con, "connection") && isOpen(con)) {
      close(con)
    }
    message("ERROR processing ", file, ": ", conditionMessage(e))
    message("Skipping to next file.")
  })
}

message("\n========================================")
message("M1S5 HTML heatmap generation finished!")
