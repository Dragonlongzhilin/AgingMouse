#' @description The H3K27ac & TE associated with H3K9me3-C1 peaks
#' @Author: Zhilin Long
#' @date 2026-8-21

library(tidyverse)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(GenomicRanges)
library(ggrepel)
library(rtracklayer)
library(cowplot)
library(gridExtra)
library(ggrepel)
set.seed(101)
source("~/longzhilin/ProjectCode/AgingMouse/Functions/colorAssign.R")
source("~/longzhilin/ProjectCode/AgingProject/plotDiffOverlap.R")
source("~/longzhilin/ProjectCode/AgingMouse/Functions/utilities.R")
setwd("~/longzhilin/Projects/AgingMouse/PairedTag/3.histoneLevel/H3K9me3/peak")
PT.info <- readRDS("~/longzhilin/Projects/AgingMouse/PairedTag/1.dataPreprocessing/PT.info.rds")
PT.info.H3K9me3.cellType <- PT.info %>% dplyr::select(tissue_cellType, cellType, tissue, tissue_new, tissue_cellType_new, category, H3K9me3.depth, colors) %>% group_by(tissue_cellType) %>% summarise(cellType = unique(cellType), read = sum(H3K9me3.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), tissue_cellType_new = unique(tissue_cellType_new), category = unique(category), colors = unique(colors), .groups = "drop")
PT.info.H3K27ac.cellType <- PT.info %>% dplyr::select(tissue_cellType, cellType, tissue, tissue_new, tissue_cellType_new, category, H3K27ac.depth, colors) %>% group_by(tissue_cellType) %>% summarise(cellType = unique(cellType), read = sum(H3K27ac.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), tissue_cellType_new = unique(tissue_cellType_new), category = unique(category), colors = unique(colors), .groups = "drop")
PT.info.RNA.cellType <- PT.info %>% dplyr::select(tissue_cellType, cellType, tissue, tissue_new, tissue_cellType_new, category, RNA.depth, colors) %>% group_by(tissue_cellType) %>% summarise(cellType = unique(cellType), read = sum(RNA.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), tissue_cellType_new = unique(tissue_cellType_new), category = unique(category), colors = unique(colors), .groups = "drop")

# calculate the average H3K9me3 signal change in C1 peaks for each tissue_cellType
cpms.fc <- readRDS("differentialPeakAnalysis/peak_cpms.fc.tissue.rds")
heatmap.order <- readRDS("heatmap.order.tissue.final.rds")
rank_cluster <- cpms.fc[which(rownames(cpms.fc) %in% heatmap.order$row_order$peak[which(heatmap.order$row_order$cluster == "C1")]),]
rank_cluster[is.na(rank_cluster)] <- 0
rank_order <- colMeans(rank_cluster, na.rm = T)

## TE information
mm10_rmsk_TE <- import("~/longzhilin/Softwares/genomes/mouse/mm10/TE/mm10_rmsk_TE.gtf") # TEtranscripts, 1-base, 3725827 regions
mm10_rmsk_TE.class <- split(mm10_rmsk_TE, mm10_rmsk_TE$class_id)
for (class_name in names(mm10_rmsk_TE.class)) {
  gr <- mm10_rmsk_TE.class[[class_name]]
  safe_name <- gsub("[^A-Za-z0-9]", "_", class_name)
  output_file <- file.path("~/longzhilin/Softwares/genomes/mouse/mm10/TE", paste0("mm10_TE_", safe_name, ".gtf"))
  export(gr, con = output_file, format = "gtf")
  cat(sprintf("Export: %s (total %d regions)\n", class_name, length(gr)))
}

#### ---- K27ac-TE in H3K9me3-C1 ---- ####
## 27ac_ov_H3K9me3-C1
H3K9me3_C1_27ac <- readRDS("H3K27ac_TE/H3K9me3.C1_H3K27ac.rds")
plot.data <- H3K9me3_C1_27ac
## keep the cell type in H3K9me3-heatmap
# plot.data <- subset(plot.data, tissue_cellType %in% names(heatmap.order$column_order))
plot.data$log2_odds_ratio <- log2(plot.data$fisher_odds_ratio)
plot.data$category <- factor(plot.data$category, levels = names(color.category))
plot.data$category_pos <- length(color.category) + 1 - as.integer(plot.data$category)
# plot.data$point_size <- factor(
#   ifelse(plot.data$fisher_pvalue < 0.05, "Significant", "Not significant"),
#   levels = c("Not significant", "Significant")
# )

# Store point offsets so labels connect to the exact points shown in the plot.
set.seed(123)
plot.data$point_y <- plot.data$category_pos + runif(nrow(plot.data), -0.19, 0.19)
label.data <- tail(plot.data[order(plot.data$fisher_odds_ratio), , drop = FALSE], 5)
axis_max <- max(plot.data$log2_odds_ratio)
label.data$label_x <- axis_max + 0.25
label.data$label_y <- seq(9.9, 11.9, length.out = nrow(label.data))

p1 <- ggplot(plot.data, aes(x = log2_odds_ratio, y = category_pos)) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
  geom_boxplot(aes(group = category, color = category), orientation = "y",
               width = 0.5, alpha = 0.8, outlier.shape = NA, linewidth = 0.5) +
  geom_point(
    aes(y = point_y, fill = category, color = after_scale(fill)),
    shape = 21,
    stroke = 0.2,
    size = 1,
    alpha = 0.8
  ) +
  geom_segment(data = label.data,
               aes(x = log2_odds_ratio, y = point_y, xend = label_x, yend = label_y),
               inherit.aes = FALSE, color = "grey35", linewidth = 0.3) +
  geom_text(data = label.data, aes(x = label_x, y = label_y, label = tissue_cellType_new),
            inherit.aes = FALSE, hjust = 0, size = 3, color = "black") +
  scale_y_continuous(breaks = seq_along(color.category), labels = rev(names(color.category))) +
  scale_fill_manual(values = color.category) +
  scale_color_manual(values = color.category) +
  expand_limits(x = c(min(plot.data$log2_odds_ratio), axis_max + 4)) +
  labs(
    x = "log2(Odds ratio)",
    y = "Category"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

#### significant TE number: edgeR analysis
TE.diff <- readRDS("~/suzhuojie/projects/pairedtag_27ac_4me1/result/H3K27ac/te_locus_featurecounts/diffTE.subfamily.C1.list.rds")
TE.diff.up.number <- lapply(TE.diff, function(x) sum(x$fdr < 0.05 & x$logFC > 0))
TE.diff.up.number <- data.frame(tissue_cellType = names(TE.diff.up.number), upregulated_TE_count = unlist(TE.diff.up.number))
plot.data <- TE.diff.up.number %>% left_join(PT.info.H3K9me3.cellType, by = "tissue_cellType")
# plot.data <- subset(plot.data, tissue_cellType %in% names(heatmap.order$column_order)) # keep the cell type in H3K9me3-heatmap
plot.data$category <- factor(plot.data$category, levels = names(color.category))
plot.data$category_pos <- length(color.category) + 1 - as.integer(plot.data$category)

# Store jittered positions so labels connect to the displayed points.
set.seed(123)
plot.data$point_y <- plot.data$category_pos + runif(nrow(plot.data), -0.19, 0.19)
# Select the five cell types with the highest TE counts.
label.data <- tail(plot.data[order(plot.data$upregulated_TE_count), , drop = FALSE], 5)
p2 <- ggplot(plot.data, aes(x = upregulated_TE_count, y = category_pos)) +
  geom_boxplot(
    aes(group = category, color = category),
    orientation = "y",
    width = 0.5,
    alpha = 0.8,
    outlier.shape = NA,
    linewidth = 0.5
  ) +
  geom_point(
    aes(y = point_y, color = category, fill = category),
    alpha = 0.8, size = 1, shape = 21, stroke = 0.2) +
  geom_text_repel(
    data = label.data,
    aes(
      x = upregulated_TE_count,
      y = point_y,
      label = tissue_cellType_new
    ),
    inherit.aes = FALSE,
    nudge_x = 0.8,
    direction = "y",
    hjust = 0,
    size = 1,
    color = "black",
    segment.color = "grey35",
    segment.size = 0.3,
    min.segment.length = 0,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf,
    seed = 123
  ) +
  scale_y_continuous(
    breaks = seq_along(color.category),
    labels = rev(names(color.category))
  ) +
  scale_color_manual(
    values = color.category,
    guide = "none"
  ) +
  scale_fill_manual(
    values = color.category,
    guide = "none"
  ) +
  labs(
    x = "The number of upregulated TE subfamilies",
    y = "Category"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right")

plot.data$category <- factor(plot.data$category, levels = rev(names(color.category)))
p22 <- ggboxplot(plot.data, x = "category", y = "upregulated_TE_count", color = "black", add = "jitter", add.params = list(size = 1, alpha = 0.8)) +
  labs(
    x = "",
    y = "The number of upregulated TE subfamilies"
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "none") +
  coord_flip()

pdf("H3K27ac_TE/H3K9me3-C1_K27ac_TE_boxplot_all.pdf", height = 3.5, width = 2.5)
print(p1)
print(p1 + theme(legend.position = "none"))
# print(p11)
print(p2)
print(p22)
dev.off()

#### plot the number of TE subfamily
# 1. Load the heatmap order and TE subfamily counts.
heatmap_order <- readRDS("heatmap.order.tissue.final.rds")
selected_cell_type <- names(heatmap_order$column_order)

te_subfamily_by_celltype_counts <- read.table(
  "~/suzhuojie/projects/pairedtag_27ac_4me1/result/H3K27ac/C1_up_peak_TE_subfamily_counts/peak_center_TE_number/tables/counts/TE_subfamily_by_celltype_counts.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

# 2. Map cell types to categories and select the heatmap cell types.
groups <- PT.info.H3K27ac.cellType$category
names(groups) <- PT.info.H3K27ac.cellType$tissue_cellType
missing_cell_type <- setdiff(selected_cell_type, colnames(te_subfamily_by_celltype_counts))
if (length(missing_cell_type) > 0) {
  stop(
    "Selected cell types are missing from the count table: ",
    paste(missing_cell_type, collapse = ", ")
  )
}

selected_counts <- te_subfamily_by_celltype_counts[, selected_cell_type, drop = FALSE]
# selected_counts <- te_subfamily_by_celltype_counts
selected_groups <- as.character(groups[colnames(selected_counts)])

if (anyNA(selected_groups)) {
  stop("Some selected cell types have no category assignment.")
}

# 3. Sum counts by category for each TE subfamily.
category_levels <- names(color.category)
group_sums <- do.call(
  cbind,
  lapply(category_levels, function(category_name) {
    rowSums(
      selected_counts[
        ,
        selected_groups == category_name,
        drop = FALSE
      ]
    )
  })
)

rownames(group_sums) <- rownames(selected_counts)
colnames(group_sums) <- category_levels

# 4. Rank TE subfamilies by total count and keep the top 30.
total_counts <- rowSums(group_sums)
row_order <- order(total_counts, decreasing = TRUE)

group_sums <- group_sums[row_order, , drop = FALSE]
total_counts <- total_counts[row_order]

top_n <- min(30L, nrow(group_sums))
top_te <- rownames(group_sums)[seq_len(top_n)]

plot_data <- data.frame(
  te_subfamily = rep(rownames(group_sums), times = ncol(group_sums)),
  category = rep(colnames(group_sums), each = nrow(group_sums)),
  count = as.numeric(group_sums),
  total_count = rep(total_counts, times = ncol(group_sums)),
  stringsAsFactors = FALSE
)

plot_data_top30 <- plot_data[plot_data$te_subfamily %in% top_te,]
plot_data_top30$te_subfamily <- factor(plot_data_top30$te_subfamily, levels = top_te)
total_labels <- unique(plot_data_top30[c("te_subfamily", "total_count")])
plot_data_top30$category <- factor(plot_data_top30$category, levels = names(color.category))
plot_data_top30$category <- droplevels(plot_data_top30$category)

# 5. Draw the stacked bar plot.
p <- ggplot2::ggplot(
  plot_data_top30,
  ggplot2::aes(
    x = te_subfamily,
    y = count,
    fill = category
  )
) +
  ggplot2::geom_col(width = 0.8) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.08))
  ) +
  ggplot2::scale_fill_manual(
    values = color.category,
    breaks = category_levels,
    drop = FALSE
  ) +
  ggplot2::labs(
    title = "Top 30 TE subfamilies by category",
    x = "TE subfamily",
    y = "Total count",
    fill = "Category"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, size = 10),
    panel.grid = ggplot2::element_blank(),
    text = ggplot2::element_text(colour = "black"),
    axis.text = ggplot2::element_text(colour = "black"),
    axis.title = ggplot2::element_text(colour = "black"),
    legend.text = ggplot2::element_text(colour = "black"),
    legend.title = ggplot2::element_text(colour = "black"),
    plot.title = ggplot2::element_text(colour = "black")
  )
pdf("H3K27ac_TE/TE_subfamily_count_top30.pdf", height = 5, width = 10)
print(p)
dev.off()
