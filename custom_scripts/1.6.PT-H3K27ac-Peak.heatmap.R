#' @description identification of cell type-group peaks
#' @Author: Zhilin Long
#' @date 2026-10-20

# ~/AgingMouse/PairedTag/3.histoneLevel/H3K27ac/peak/mergePeak_zhuojie/cellTypeSpecificPeak/lineage
library(tidyverse)
library(GenomicRanges)
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(cowplot)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(matrixStats)
library(grid)
set.seed(101)
set.resolutions <- seq(0.1, 1.2, by = 0.1)
source("~/Functions/colorAssign.R")
PT.info <- readRDS("~/AgingMouse/PairedTag/1.dataPreprocessing/PT.info.rds")
setwd("~/AgingMouse/PairedTag/3.histoneLevel/H3K27ac/peak/mergePeak_zhuojie/cellTypeSpecificPeak/lineage")


write_tab <- function(x, name) {
  data.table::fwrite(x, file.path(out, name), sep = "\t")
}
msg <- function(...) {
  cat(format(Sys.time(), "%F %T"), ..., "\n")
  flush.console()
}

#### ---- cluster ---- ####
out <- "~/AgingMouse/PairedTag/3.histoneLevel/H3K27ac/peak/mergePeak_zhuojie/cellTypeSpecificPeak/lineage/final_peak_CPM_20261020"
# 2. Load the selection result and check matrix/metadata alignment
# m and p retain the original selection indices; mm and pp will hold display order.
a <- readRDS(file.path(out, "hierarchical_analysis.rds"))
m <- a$map
p <- a$primary
cats <- a$parameters$category_order
groups <- a$groups
gcats <- a$group_categories
stopifnot(all(a$memberships$level == "Main enriched"), !any(p$category_hit), !any(p$local_hit), !anyDuplicated(p$peak_id),
  identical(rownames(a$log_sample), p$peak_id), identical(colnames(a$log_sample), m$tissue_cellType))

# Pearson distance is 1 - r. Clamp negative round-off to zero and force zero diagonals.
# sample_r compares profiles across all selected log peaks; group_r is saved for reference.
distance <- function(r) {
  d <- pmax(1 - r, 0)
  diag(d) <- 0
  as.dist(d)
}
group_r <- cor(a$log_group)
sample_r <- cor(a$log_sample)
stopifnot(all(is.finite(group_r)), all(is.finite(sample_r)))

# 3. Fix category/main block order before within-group clustering
# Rank mains by decreasing original profile count within each category. Resolve ties
# lexicographically using all represented tissues, then their first metadata row.
# Block order is prescribed and does not imply a global unsupervised relationship.
tissue_order <- names(color.tissue)
stopifnot(setequal(unique(m$tissue), tissue_order))
count_order <- rbindlist(lapply(groups, function(g) {
  jj <- which(m$main_cellType == g)
  ti <- sort(match(unique(m$tissue[jj]), tissue_order))
  data.frame(main_cellType = g, category = unname(gcats[g]), category_rank = match(gcats[g], cats),
    n_celltypes = length(jj), n_distinct_cellType_labels = length(unique(m$cellType[jj])), first_tissue = tissue_order[ti[1]],
    first_tissue_rank = ti[1], all_tissues = paste(tissue_order[ti], collapse = ";"), tissue_order_key = paste(sprintf("%02d",
      c(ti, rep(99L, length(tissue_order) - length(ti)))), collapse = "-"), map_first_row = min(jj))
}))
# Compare complete tissue lists after the first tissue; exact ties keep mapping order.
setorderv(count_order, c("category_rank", "n_celltypes", "tissue_order_key", "map_first_row"), c(1,
  -1, 1, 1))
group_order <- count_order$main_cellType
stopifnot(sum(count_order$n_celltypes) == nrow(m), setequal(group_order, groups))
for (ct in cats) {
  z <- count_order[category == ct]
  stopifnot(all(diff(z$n_celltypes) <= 0))
  for (n in unique(z$n_celltypes)) {
    zz <- z[n_celltypes == n]
    stopifnot(identical(order(zz$tissue_order_key, zz$map_first_row), seq_len(nrow(zz))))
  }
}
write_tab(data.frame(order = seq_len(nrow(count_order)), count_order), "main_count_tissue_order.tsv")

# 4. Cluster profiles independently within each main group
# Use Pearson distance across all selected peaks and average linkage.
sample_order <- unlist(lapply(group_order, function(g) {
  jj <- which(m$main_cellType == g)
  if (length(jj) > 1L)
    jj <- jj[hclust(distance(sample_r[jj, jj, drop = FALSE]), method = "average")$order]
  jj
}), use.names = FALSE)
mm <- m[sample_order, , drop = FALSE]
mm <- left_join(mm, unique(PT.info[, c("tissue_cellType", "tissue_new", "colors")]), by = "tissue_cellType")

# 5. Build main-specific peak blocks and cluster peaks within each block
# Peak correlations use all original profiles. No peak is dropped for display.
group_cat <- unname(gcats[group_order])
p$display_block <- ifelse(p$level == "Category shared", paste(p$category, "Shared", sep = "::"),
  paste(p$category, p$main_cellType, p$level, sep = "::"))
block_order <- unlist(lapply(cats, function(ct) {
  g <- group_order[gcats[group_order] == ct]
  candidates <- unlist(lapply(g, function(z) paste(ct, z, "Main enriched", sep = "::")), use.names = FALSE)
  candidates[candidates %in% p$display_block]
}), use.names = FALSE)
msg("Clustering", nrow(p), "distinct peaks in", length(block_order), "hierarchical blocks")
peak_order <- unlist(lapply(block_order, function(b) {
  ii <- which(p$display_block == b)
  if (length(ii) > 1L) {
    r <- cor(t(a$log_sample[ii, , drop = FALSE]))
    ii <- ii[hclust(distance(r), method = "average")$order]
  }
  ii
}), use.names = FALSE)
stopifnot(setequal(peak_order, seq_len(nrow(p))), setequal(sample_order, seq_len(nrow(m))))

# Prepare ordered annotations and fixed color mappings.
pp <- p[peak_order, , drop = FALSE]
category_cols <- color.category[cats]
# Save the prepared display indices before rendering.
write_tab(data.frame(order = seq_len(nrow(mm)), mm), "all195_cluster_order.tsv")
write_tab(data.frame(order = seq_len(nrow(pp)), pp), "peak_display_order.tsv")
saveRDS(list(row_indices = sample_order, peak_indices = peak_order), file.path(out, "prepared_display_layout.rds"))
# 6. Category and tissue annotation tracks
tissue_cols <- color.tissue
lgp <- list(title_gp = gpar(fontsize = 8, fontface = "bold"), labels_gp = gpar(fontsize = 7))
top_tracks <- function() HeatmapAnnotation(
  Category = factor(pp$category, levels = cats),
  col = list(Category = category_cols),
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 7),
  simple_anno_size = unit(2.4, "mm"),
  annotation_legend_param = list(Category = c(list(title = "Category", at = cats), lgp))
)
cellType_cols <- unique(mm$colors)
names(cellType_cols) <- unique(mm$cellType)
left_tracks <- function() rowAnnotation(
  Category = factor(mm$category, levels = cats),
  Tissue = factor(mm$tissue, levels = tissue_order),
  CellType = factor(mm$cellType, levels = unique(mm$cellType)),
  col = list(Category = category_cols, Tissue = tissue_cols, CellType = cellType_cols),
  simple_anno_size = unit(2.5, "mm"),
  show_legend = c(TRUE, TRUE, FALSE),
  show_annotation_name = FALSE,
  annotation_name_gp = gpar(fontsize = 7),
  annotation_legend_param = list(Tissue = c(list(title = "Tissue", at = tissue_order), lgp))
)

# Reuse profile correlations computed from log2CPM for both heatmap versions.
# This keeps their row dendrograms and final orders identical.
cluster_members <- function(mat) {
  ids <- rownames(mat)
  if (length(ids) == 1L)
    return(structure(1L, label = ids, members = 1L, height = 0, leaf = TRUE, class = "dendrogram"))
  hclust(distance(sample_r[ids, ids, drop = FALSE]), method = "average")
}

# 7. Display matrices and heatmap construction
# Matrices here are peaks x profiles. Z-score each log peak across all profiles
# using the sample standard deviation; transpose only when passing data to Heatmap.
# Color-map endpoints saturate displayed colors without clipping matrix values.
log_mat <- a$log_sample
z_mat <- (log_mat - rowMeans(log_mat)) / rowSds(log_mat)
stopifnot(max(abs(rowMeans(z_mat))) < 1e-10, max(abs(rowSds(z_mat) - 1)) < 1e-10)
log_col <- colorRamp2(c(0, 1.5, 3, 4.5, 6), c("#FFFFFF", "#FEE8C8", "#edb78b", "#f76a54", "#a50417"))
z_col <- colorRamp2(c(-2, 0, 2, 4), c("#2166AC", "#F7F7F7", "#D6604D", "#B2182B"))

# Respect precomputed main/peak blocks. Only within-main row dendrograms are shown;
# column clustering has already been calculated explicitly in peak_order.
make_heatmap <- function(mat, name, cols, legend) Heatmap(
  t(mat[peak_order, sample_order, drop = FALSE]),
  name = name,
  col = cols,
  # Set the matrix body size independently of annotations and the output page.
  width = unit(8, "cm"),
  height = unit(14, "cm"),
  row_split = factor(mm$main_cellType, levels = group_order),
  cluster_rows = cluster_members,
  row_dend_reorder = FALSE,
  cluster_row_slices = FALSE,
  row_gap = unit(0, "mm"),
  cluster_columns = FALSE,
  column_split = factor(pp$display_block, levels = block_order),
  cluster_column_slices = FALSE,
  column_gap = unit(0, "mm"),
  column_title = NULL,
  # Full profile and group labels remain in the ordered annotation tables.
  row_title = NULL,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 7),
  row_dend_width = unit(18, "mm"),
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 6),
  row_names_max_width = unit(90, "mm"),
  left_annotation = left_tracks(),
  top_annotation = top_tracks(),
  use_raster = TRUE,
  raster_device = "png",
  raster_resize_mat = TRUE,
  raster_quality = 2,
  heatmap_legend_param = c(legend, lgp)
)
# 8. Render the same layout as PDF and PNG
# Signal cells are rasterized; PDF text and annotation lines remain vector objects.
# Measure rendered slice viewports so the requested body size is verified.
rendered_geometry <- list()
measure_body <- function(obj, suffix, device) {
  heatmap_name <- obj@ht_list[[1L]]@name
  slice_size <- function(i, j, dimension) {
    seekViewport(paste0(heatmap_name, "_heatmap_body_", i, "_", j))
    size <- if (dimension == "width") {
      convertWidth(unit(1, "npc"), "cm", valueOnly = TRUE)
    } else {
      convertHeight(unit(1, "npc"), "cm", valueOnly = TRUE)
    }
    upViewport(0)
    size
  }
  widths <- vapply(seq_along(block_order), function(j) slice_size(1L, j, "width"), numeric(1))
  heights <- vapply(seq_along(group_order), function(i) slice_size(i, 1L, "height"), numeric(1))
  stopifnot(abs(sum(widths) - 8) < 1e-5, abs(sum(heights) - 14) < 1e-5)
  rendered_geometry[[paste(suffix, device, sep = "_")]] <<- list(
    body_width_cm = sum(widths), body_height_cm = sum(heights),
    row_gap_mm = 0, column_gap_mm = 0, page_width_cm = 20, page_height_cm = 17
  )
}

render_output <- function(ht, suffix) {
  stem <- paste0("Final_", suffix)
  render <- function(device) {
    obj <- draw(ht, background = "transparent", merge_legends = TRUE, heatmap_legend_side = "right",
      annotation_legend_side = "right",
      padding = unit(c(5, 5, 15, 5), "mm"))
    grid.text(sprintf("H3K27ac | %d profiles | %s peaks",
      nrow(m), format(nrow(p), big.mark = ",")),
      x = 0.5, y = unit(1, "npc") - unit(5, "mm"),
      gp = gpar(fontsize = 9, fontface = "bold"))
    grid.text(if (suffix == "log2CPM") "log2(CPM + 1)" else "Peak Z-score of log2(CPM + 1)",
      x = 0.5, y = unit(1, "npc") - unit(9, "mm"), gp = gpar(fontsize = 7))
    measure_body(obj, suffix, device)
    invisible(obj)
  }
  pdf(file.path(out, paste0(stem, ".pdf")), width = 34, height = 30,
    family = "Helvetica", useDingbats = FALSE)
  obj <- render()
  dev.off()
  ragg::agg_png(file.path(out, paste0(stem, ".png")), width = 20, height = 17,
    units = "cm", res = 300, background = "white")
  render()
  dev.off()
  invisible(obj)
}

# 9. Draw both heatmaps and validate the actual displayed orders
msg("Rendering log2CPM heatmap")
hlog <- render_output(make_heatmap(log_mat, "log2(CPM+1)", log_col, list(at = 0:6, labels = c(as.character(0:5),
  "6+"))), "log2CPM")
msg("Rendering Z-score heatmap")
hz <- render_output(make_heatmap(z_mat, "Peak Z-score", z_col, list(at = c(-2, 0, 2, 4), labels = c("-2 or less",
  "0", "2", "4 or more"))), "Zscore")

# Query ComplexHeatmap after drawing rather than assuming the prepared order survived.
rl <- unlist(row_order(hlog), use.names = FALSE)
rz <- unlist(row_order(hz), use.names = FALSE)
cl <- unlist(column_order(hlog), use.names = FALSE)
cz <- unlist(column_order(hz), use.names = FALSE)
stopifnot(identical(rl, rz), identical(cl, cz), identical(rl, seq_len(nrow(mm))), identical(cl, seq_len(nrow(pp))),
  identical(unique(mm$category[rl]), cats), !anyDuplicated(pp$peak_id[cl]))
write_tab(data.frame(order = seq_len(nrow(mm)), mm[rl, , drop = FALSE]), "all195_cluster_order.tsv")
write_tab(data.frame(order = seq_along(group_order), category = group_cat, main_cellType = group_order),
  "main_display_order.tsv")
write_tab(data.frame(order = seq_len(nrow(pp)), pp[cl, , drop = FALSE]), "peak_display_order.tsv")
write_tab(data.frame(tissue_cellType = rownames(sample_r), sample_r, check.names = FALSE), "all195_Pearson_correlation.tsv")
saveRDS(list(row_indices = sample_order[rl], peak_indices = peak_order[cl], group_order = group_order,
  block_order = block_order, sample_r = sample_r, group_r = group_r), file.path(out, "display_layout.rds"))


#### ---- produce png
block_cats <- pp$category[match(block_order, pp$display_block)]
cgaps <- ifelse(head(block_cats, -1) == tail(block_cats, -1), 0.5, 2)
rgaps <- ifelse(head(group_cat, -1) == tail(group_cat, -1), 1, 3)
make_heatmap <- function(mat, name, cols, legend) Heatmap(
  t(mat[peak_order, sample_order, drop = FALSE]),
  name = name,
  col = cols,
  row_split = factor(mm$main_cellType, levels = group_order),
  cluster_rows = cluster_members,
  show_row_dend = FALSE,
  row_dend_reorder = FALSE,
  cluster_row_slices = FALSE,
  row_gap = unit(rgaps, "mm"),
  cluster_columns = FALSE,
  column_split = factor(pp$display_block, levels = block_order),
  cluster_column_slices = FALSE,
  column_gap = unit(cgaps, "mm"),
  column_title = NULL,
  row_title = NULL,
  show_column_names = FALSE,
  show_row_names = FALSE,
  row_names_max_width = unit(90, "mm"),
  left_annotation = left_tracks(),
  top_annotation = top_tracks(),
  use_raster = TRUE,
  raster_device = "png",
  raster_resize_mat = TRUE,
  raster_quality = 2,
  heatmap_legend_param = c(legend, lgp)
)
render_output <- function(ht, suffix, subtitle) {
  stem <- paste0("H3K27ac_big", suffix)
  render <- function() {
    obj <- draw(ht, merge_legends = TRUE, heatmap_legend_side = "right", annotation_legend_side = "right", padding = unit(c(6, 6, 25, 6), "mm"))
    # Draw one border around the entire heatmap body.
    heatmap_name <- obj@ht_list[[1L]]@name
    seekViewport(paste0(heatmap_name, "_heatmap_body_wrap"))
    grid.rect(gp = gpar(fill = NA, col = "black", lwd = 0.8))
    upViewport(0)
    invisible(obj)
  }
  pdf(file.path(out, paste0(stem, ".pdf")), width = 20, height = 30, family = "Helvetica", useDingbats = FALSE)
  obj <- render()
  dev.off()
  invisible(obj)
}

# 10. Draw both heatmaps and validate the actual displayed orders
msg("Rendering log2CPM heatmap")
hlog <- render_output(make_heatmap(log_mat, "log2(CPM+1)", log_col, list(at = 0:6, labels = c(as.character(0:5),
  "6+"))), "log2CPM", "Unstandardized log2(CPM + 1) | main enriched only; category shared and local selection disabled")
# msg("Rendering Z-score heatmap")
# hz <- render_output(make_heatmap(z_mat, "Peak Z-score", z_col, list(at = c(-2, 0, 2, 4), labels = c("-2 or less",
#   "0", "2", "4 or more"))), "Zscore", sprintf("Peak Z-score of log2(CPM + 1) across all %d original profiles | same layout as log2CPM",
#   nrow(m)))