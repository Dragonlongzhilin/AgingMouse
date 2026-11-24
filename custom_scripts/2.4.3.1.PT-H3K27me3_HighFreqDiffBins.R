#' @description read change in bins among different age
#' @Author: Zhilin Long
#' @date 2026-8-24

library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(gridExtra)
library(GenomicRanges)
library(GenomeInfoDb)
set.seed(101)
source("~/Functions/colorAssign.R")
source("~/Functions/plotDiffOverlap.R")
source("~/Functions/utilities.R")
setwd("~/longzhilin/Projects/AgingMouse/PairedTag/3.histoneLevel/H3K27me3")

histone <- "H3K27me3"
log2FC.threshold <- log2(1.2)
sig.label <- "PValue"

# summarise by tissue, age and cell type
PT.info <- readRDS("~/longzhilin/Projects/AgingMouse/PairedTag/1.dataPreprocessing/PT.info.rds")
PT.info.H3K27me3.cellType <- PT.info %>% select(tissue_cellType, tissue_cellType_new, cellType, tissue, tissue_new, category, H3K27me3.depth) %>% group_by(tissue_cellType) %>% summarise(tissue_cellType_new = unique(tissue_cellType_new), cellType = unique(cellType), read = sum(H3K27me3.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), category = unique(category), .groups = "drop")
PT.info.H3K27me3.tissue <- PT.info %>% select(tissue, H3K27me3.depth, tissue_new) %>% group_by(tissue) %>% summarise(read = sum(H3K27me3.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), .groups = "drop")
PT.info.H3K27me3.sample <- PT.info %>% select(tissue_cellType, tissue_cellType_new, sample, age_group, category, H3K27me3.depth) %>% group_by(tissue_cellType, sample) %>% summarise(read = sum(H3K27me3.depth, na.rm = T), age_group = unique(age_group), category = unique(category), tissue_cellType_new = unique(tissue_cellType_new), .groups = "drop") %>% filter(read > 0)
PT.info.H3K27me3.cellType_age <- PT.info %>% dplyr::select(tissue_cellType, tissue_cellType_new, cellType, age_group, tissue, tissue_new, category, H3K27me3.depth, colors) %>% group_by(tissue_cellType, age_group) %>% summarise(tissue_cellType_new = unique(tissue_cellType_new), age_group = unique(age_group), cellType = unique(cellType), read = sum(H3K27me3.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), category = unique(category), colors = unique(colors), .groups = "drop")

# set sample palette
sample.info <- readRDS("~/longzhilin/Projects/AgingMouse/PairedTag/1.dataPreprocessing/sample.info.rds")
sample.info <- sample.info[sample.info$histone == histone,]
info.Old <- unique(sample.info$replication[which(sample.info$age_group == "Old")])
info.Young <- unique(sample.info$replication[which(sample.info$age_group == "Young")])
Old.sample.colors <- colorRampPalette(color.olds)(length(info.Old))
Young.sample.colors <- colorRampPalette(color.youngs)(length(info.Young))
color.samples <- c(Young.sample.colors, Old.sample.colors)
names(color.samples) <- c(info.Young, info.Old)

#### plot read depths
plot.data <- unique(PT.info.H3K27me3.sample[,c("tissue_cellType", "sample", "age_group", "read")])
plot.data$tissue <- gsub(":.*", "", plot.data$tissue_cellType)
plot.data <- plot.data %>% left_join(sample.info, by = c("sample", "age_group", "tissue"))
idx <- which(plot.data$age_group == "Young")
plot.data$read[idx] <- plot.data$read[idx] * (-1)
plot.data <- tidyr::pivot_wider(plot.data, id_cols = c(tissue_cellType, replication), names_from = age_group, values_from = read)
plot.data_sum <- plot.data %>% group_by(tissue_cellType) %>% summarise(Old_sum = sum(Old, na.rm = TRUE), Young_sum = sum(Young, na.rm = TRUE))
total_sum <- tapply(PT.info.H3K27me3.sample$read, PT.info.H3K27me3.sample$tissue_cellType, sum)
plot.data$tissue_cellType <- factor(plot.data$tissue_cellType, levels = names(sort(total_sum, decreasing = T)))
symmetric_limit <- max(abs(na.omit(plot.data_sum$Old_sum)), abs(na.omit(plot.data_sum$Young_sum)))
pdf("H3K27me3.depths.pdf", height = 5, width = 30)
ggplot(plot.data, aes(x = tissue_cellType)) +
geom_bar(aes(y = Old, fill = replication), stat = "identity", width = 0.5) + 
geom_bar(aes(y = Young, fill = replication), stat = "identity", width = 0.5) +
scale_fill_manual(name = "replication", values = c(color.samples)) + 
geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
scale_y_continuous(limits = c(-symmetric_limit, symmetric_limit)) + 
geom_text(data = plot.data_sum, aes(x = tissue_cellType, y = Old_sum, label = round(Old_sum, 2)), color = "black", size = 3, angle = 90, hjust = -0.5) +
geom_text(data = plot.data_sum, aes(x = tissue_cellType, y = Young_sum, label = abs(round(Young_sum, 2))), color = "black", angle = 90, size = 3, hjust = 1.5) +
labs(x = "", y = "Read depths(M)") + theme_classic(base_size = 12) + theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(color = "black"), axis.ticks = element_line(color = "black"), axis.line = element_line(color = "black"), legend.position = "top", axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

## select cell types with at least 0.5M reads in each age group
selected.tissue_cellType <- PT.info.H3K27me3.cellType_age %>% group_by(tissue_cellType) %>% filter(all(read >= 0.5)) %>% ungroup() %>% dplyr::select(tissue_cellType, tissue) %>% distinct() # 115 cell types

#### ---- differential bins overlapping young peaks ---- ####
cellType_Young_peakList <- readRDS("peak/cellType_Young_peakList.rds")
names(cellType_Young_peakList) <- gsub("_Young", "", names(cellType_Young_peakList))
diffBin.list <- readRDS("bin10kb/differentialBinAnalysis/diffBin.list.rds")
binOverlapingPeak <- lapply(names(diffBin.list), function(tissue_cellType){
  cat("Analysis ", tissue_cellType, "\n")
  diffBin <- diffBin.list[[tissue_cellType]]
  peak <- cellType_Young_peakList[[tissue_cellType]]
  diffBin.gr <- GRanges(seqnames = gsub(":.*", "", diffBin$Region), ranges = IRanges(start = as.numeric(gsub("-.*", "", gsub(".*:", "", diffBin$Region))), end = as.numeric(gsub(".*-", "", diffBin$Region))), strand = "*")
  idx1 <- which(diffBin[, sig.label] < 0.05 & diffBin$logFC > log2FC.threshold)
  if(length(idx1) > 0){
    diffBin_up.gr <- diffBin.gr[idx1]
    up.idx <- findOverlaps(diffBin_up.gr, peak)
    up.number <- length(queryHits(up.idx))
    up.ratio <- up.number / length(diffBin_up.gr)
  }else{
    up.number <- 0
    up.ratio <- 0
  }

  idx2 <- which(diffBin[, sig.label] < 0.05 & diffBin$logFC < (log2FC.threshold * -1))
  if(length(idx2) > 0){
    diffBin_down.gr <- diffBin.gr[idx2]
    down.idx <- findOverlaps(diffBin_down.gr, peak)
    down.number <- length(queryHits(down.idx))
    down.ratio <- down.number / length(diffBin_down.gr)
  }else{
    down.number <- 0
    down.ratio <- 0
  }
  info <- data.frame(tissue_cellType = tissue_cellType, Peak = length(peak), diffBin_UP = length(idx1), diffBin_UP.overlappingPeak = up.number, diffBin_UP.overlappingPeak.ratio = up.ratio, diffBin_DOWN = length(idx2), diffBin_DOWN.overlappingPeak = down.number, diffBin_DOWN.overlappingPeak.ratio = down.ratio)
  return(info)
})
binOverlapingPeak <- do.call(rbind, binOverlapingPeak)
binOverlapingPeak <- binOverlapingPeak %>% left_join(PT.info.H3K27me3.cellType, by = "tissue_cellType")
saveRDS(binOverlapingPeak, "bin10kb/differentialBinAnalysis/binOverlapingPeak.rds")

## plot ratio
binOverlapingPeak$tissue_new <- factor(binOverlapingPeak$tissue_new, levels = names(color.tissue2))
p_scatter <- ggscatter(binOverlapingPeak, x = "diffBin_UP.overlappingPeak.ratio", y = "diffBin_DOWN.overlappingPeak.ratio", xlab = "Proportion of upregulated bins overlapping with peaks", ylab = "Proportion of downregulated bins overlapping with peaks", color = "tissue_new", palette = color.tissue2, size = "read", alpha = 0.8) + scale_size_continuous(range = c(1, 4), name = "Read(M)") + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") + scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0, 1)) + theme(aspect.ratio = 1, legend.position = "right")

p_scatter_up <- ggscatter(binOverlapingPeak, x = "diffBin_UP", y = "diffBin_UP.overlappingPeak.ratio", ylab = "Proportion of upregulated bins overlapping with peaks", xlab = "Number of upregulated bins", color = "tissue_new", palette = color.tissue2, size = "Peak", alpha = 0.8) + scale_size_continuous(range = c(1, 4), name = "Peak number") + theme(aspect.ratio = 1, legend.position = "right") + geom_smooth(method = "lm", se = TRUE, color = "#CE3D32", linetype = "dashed", aes(group = 1)) + stat_cor(method = "spearman", label.x = max(binOverlapingPeak$diffBin_UP, na.rm = TRUE) * 0.05, label.y = 0.95, size = 4,digits = 3, aes(group = 1))
p_scatter_down <- ggscatter(binOverlapingPeak, x = "diffBin_DOWN", y = "diffBin_DOWN.overlappingPeak.ratio", ylab = "Proportion of downregulated bins overlapping with peaks", xlab = "Number of downregulated bins", color = "tissue_new", palette = color.tissue2, size = "Peak", alpha = 0.8) + scale_size_continuous(range = c(1, 4), name = "Peak number") + theme(aspect.ratio = 1, legend.position = "right") + geom_smooth(method = "lm", se = TRUE, color = "#3A6CAD", linetype = "dashed", aes(group = 1)) + stat_cor(method = "spearman", label.x = max(binOverlapingPeak$diffBin_UP, na.rm = TRUE) * 0.05, label.y = 0.95, size = 4,digits = 3, aes(group = 1))

p_test <- ggpaired(binOverlapingPeak, cond1 = "diffBin_UP.overlappingPeak.ratio", cond2 = "diffBin_DOWN.overlappingPeak.ratio", color = "condition", line.color = "#D1D1D1", line.size = 0.1, palette = c("#CE3D32", "#3A6CAD"), add.params = list(size = 0.5), xlab = "", ylab = "Proportion of dyregulated bins overlapping with peaks") + stat_compare_means(paired = TRUE, method = "wilcox.test", label = "p.signif", label.x = 1.5) + theme(legend.position = "none", aspect.ratio = 1.5)

## violin plot 
# Keep complete pairs for both plotting and statistical testing.
paired_data <- binOverlapingPeak %>%
  transmute(
    Increased = as.numeric(diffBin_UP.overlappingPeak.ratio),
    Decreased = as.numeric(diffBin_DOWN.overlappingPeak.ratio)
  ) %>%
  filter(complete.cases(.))

# Run a paired Wilcoxon signed-rank test.
paired_test <- wilcox.test(
  paired_data$Increased,
  paired_data$Decreased,
  paired = TRUE,
  exact = FALSE,
  conf.int = TRUE
)

# Summarize paired differences as Increased minus Decreased.
paired_difference <- paired_data$Increased - paired_data$Decreased

difference_summary <- data.frame(
  n_pairs = nrow(paired_data),
  median_difference = median(paired_difference),
  mean_difference = mean(paired_difference),
  Hodges_Lehmann_estimate = unname(paired_test$estimate),
  CI_lower = paired_test$conf.int[1],
  CI_upper = paired_test$conf.int[2]
)

# Reshape the paired measurements into long format for plotting.
plot_data <- paired_data %>%
  pivot_longer(
    cols = c(Increased, Decreased),
    names_to = "Regulation",
    values_to = "Proportion"
  ) %>%
  mutate(
    Regulation = factor(Regulation, levels = c("Increased", "Decreased"))
  )

p_annotation <- data.frame(
  group1 = "Increased",
  group2 = "Decreased",
  y.position = 1.08,
  label = paired_test$p.value
)

# Plot the proportions and display the paired-test p-value.
p_violin <- ggviolin(
  plot_data,
  x = "Regulation",
  y = "Proportion",
  color = "Regulation",
  palette = c("#CE3D32", "#3A6CAD"),
  add = "boxplot",
  add.params = list(fill = "white"),
  trim = TRUE
) +
  stat_pvalue_manual(
    p_annotation,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 4
  ) +
  labs(
    x = NULL,
    y = "Proportion of dysregulated bins overlapping with peaks"
  ) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1.15)) +
  theme_classic() +
  theme(legend.position = "none")

## boxplot: diffBin_DOWN.overlappingPeak.ratio by tissue (sorted by median)
rank.order <- tapply(binOverlapingPeak$diffBin_DOWN.overlappingPeak.ratio, binOverlapingPeak$tissue_new, median)
binOverlapingPeak$tissue_new <- factor(binOverlapingPeak$tissue_new, levels = names(sort(rank.order, decreasing = T)))
p_tissue <- ggboxplot(binOverlapingPeak, x = "tissue_new", y = "diffBin_DOWN.overlappingPeak.ratio", add = "jitter", color = "tissue_new", palette = color.tissue2, xlab = "", ylab = "Proportion of downregulated bins overlapping with peaks") + ggtitle("diffBin_DOWN ratio by tissue") + theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), plot.title = element_text(hjust = 0.5), aspect.ratio = 1)

## boxplot: diffBin_DOWN.overlappingPeak.ratio by category (sorted by median)
rank.order <- tapply(binOverlapingPeak$diffBin_DOWN.overlappingPeak.ratio, binOverlapingPeak$category, median)
binOverlapingPeak$category <- factor(binOverlapingPeak$category, levels = names(sort(rank.order, decreasing = T)))
p_category <- ggboxplot(binOverlapingPeak, x = "category", y = "diffBin_DOWN.overlappingPeak.ratio", add = "jitter", color = "category", palette = color.category, xlab = "", ylab = "Proportion of downregulated bins overlapping with peaks") + ggtitle("diffBin_DOWN ratio by category") + theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1), plot.title = element_text(hjust = 0.5), aspect.ratio = 1)

pdf("bin10kb/differentialBinAnalysis/binOverlapingPeak.pdf", height = 5, width = 5)
print(p_scatter)
print(ggarrange(p_test, p_violin, nrow = 2, ncol = 2))
print(p_scatter_up)
print(p_scatter_down)
print(p_tissue)
print(p_category)
dev.off()

#### ---- 1.analysis the high frequency bins ---- ####
diffBin.Frequency <- readRDS("bin10kb/differentialBinAnalysis/diffBin.Frequency.rds")
cpm.log2fc <- readRDS("bin10kb/differentialBinAnalysis/cpm.log2fc.rds")
n_cellTypes <- 40
#### category distribution of differential bins
up_highFreq <- diffBin.Frequency[which(diffBin.Frequency$n_sig_up_cellType >= n_cellTypes),] # 17108 bins
down_highFreq <- diffBin.Frequency[which(diffBin.Frequency$n_sig_down_cellType >= n_cellTypes),] # 3313 bins
diffBin.highFreq <- c(as.character(up_highFreq$Region), as.character(down_highFreq$Region))
diffBin.highFreq.log2fc <- cpm.log2fc[diffBin.highFreq,]

## plot the high-frequency of different peak
plot.data.count <- diffBin.Frequency
plot.data.count$group <- "lowFreq"
plot.data.count$group[which(plot.data.count$n_sig_up_cellType >= n_cellTypes)] <- "up-highFreq"
plot.data.count$group[which(plot.data.count$n_sig_down_cellType >= n_cellTypes)] <- "down-highFreq"
plot.data.count$group[which(plot.data.count$n_sig_down_cellType >= n_cellTypes & plot.data.count$n_sig_up_cellType >= n_cellTypes)] <- "up-down-highFreq"
xy_limit <- max(plot.data[,c("n_sig_up_cellType", "n_sig_down_cellType")], na.rm = T)
p <- ggplot(plot.data.count, aes(x = n_sig_up_cellType, y = n_sig_down_cellType, color = group)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(name = "group", values = c("lowFreq" = "#D1D1D1", "up-highFreq" = "#CE3D32", "down-highFreq" = "#3A6CAD", "up-down-highFreq" = "black")) +
  xlim(0, xy_limit) + ylim(0, xy_limit) +
  geom_vline(xintercept = n_cellTypes, linetype = "dashed", color = "#BEBEBE", linewidth = 0.5) +
  geom_hline(yintercept = n_cellTypes, linetype = "dashed", color = "#BEBEBE", linewidth = 0.5) +
  annotate("text", x = xy_limit * 0.85, y = xy_limit * 0.05,
           label = paste0("up-highFreq: ", nrow(up_highFreq), " bins\n(n_sig_up >= 40)"),
           color = "#CE3D32", size = 3, hjust = 1, vjust = 0) +
  annotate("text", x = xy_limit * 0.05, y = xy_limit * 0.85,
           label = paste0("down-highFreq: ", nrow(down_highFreq), " bins\n(n_sig_down >= 40)"),
           color = "#3A6CAD", size = 3, hjust = 0.75, vjust = 0.5, angle = 90) +
  ggtitle(paste0("Total differential bins: ", nrow(plot.data))) +
  xlab("# of cell types with increased H3K27me3") +
  ylab("# of cell types with decreased H3K27me3") +
  theme_classic() +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(hjust = 0.5)
  )
pdf("bin10kb/diffBin.HighFrequency.bubblePlot.pdf", width = 5, height = 5)
print(p)
dev.off()

## overlap with H3K9me3 peaks
H3K9me3_tissue_peak_C1 <- read.table("~/longzhilin/public_html/AgingMouse/PairedTag_new/MergeIGV/H3K9me3/H3K9me3_tissue_peak_C1.bed")
H3K9me3_tissue_peak_C1 <- makeGRangesFromDataFrame(H3K9me3_tissue_peak_C1, keep.extra.columns = T, seqnames.field = "V1", start.field = "V2", end.field = "V3")
H3K9me3_tissue_peak_C1$type <- "C1"
H3K9me3_tissue_peak_C2 <- read.table("~/longzhilin/public_html/AgingMouse/PairedTag_new/MergeIGV/H3K9me3/H3K9me3_tissue_peak_C2.bed")
H3K9me3_tissue_peak_C2 <- makeGRangesFromDataFrame(H3K9me3_tissue_peak_C2, keep.extra.columns = T, seqnames.field = "V1", start.field = "V2", end.field = "V3")
H3K9me3_tissue_peak_C2$type <- "C2"
H3K9me3_tissue_domain <- c(H3K9me3_tissue_peak_C1, H3K9me3_tissue_peak_C2)

up_highFreq.gr <- GRanges(up_highFreq$Region)
down_highFreq.gr <- GRanges(down_highFreq$Region)
other.gr <- GRanges(setdiff(diffBin.Frequency$Region, c(up_highFreq$Region, down_highFreq$Region)))
bin_by_group <- list(up_highFreq.bin = up_highFreq.gr, down_highFreq.bin = down_highFreq.gr, other.bin = other.gr)
status_levels <- c("C1 only", "C2 only", "C1 + C2", "No overlap")
plot_data <- do.call(rbind, lapply(names(bin_by_group), function(group_name) {
  bins <- bin_by_group[[group_name]]

  overlap_C1 <- overlapsAny(bins, H3K9me3_tissue_peak_C1)
  overlap_C2 <- overlapsAny(bins, H3K9me3_tissue_peak_C2)

  status <- ifelse(
    overlap_C1 & overlap_C2, "C1 + C2",
    ifelse(overlap_C1, "C1 only",
           ifelse(overlap_C2, "C2 only", "No overlap"))
  )

  counts <- table(factor(status, levels = status_levels))
  total_bins <- length(bins)

  data.frame(
    group = group_name,
    status = names(counts),
    bin_count = as.integer(counts),
    total_bins = total_bins,
    proportion = if (total_bins > 0) {
      as.integer(counts) / total_bins
    } else {
      NA_real_
    }
  )
}))
print(plot_data, row.names = FALSE)
saveRDS(plot_data, "bin10kb/diffBin.highFreq.H3K9me3_peakOverlap.rds")
plot_data <- subset(plot_data, bin_count > 0)
plot_data$group <- factor(plot_data$group, levels = c("up_highFreq.bin", "down_highFreq.bin", "other.bin"))
plot_data <- subset(plot_data, status != "No overlap")
plot_data$status <- droplevels(factor(plot_data$status, levels = status_levels))
p <- ggplot(plot_data, aes(x = group, y = proportion, fill = status)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", proportion * 100)), position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(
    values = c(
      "C1 only" = "#E18727",
      "C2 only" = "#00A087",
      "C1 + C2" = "#9B59B6",
      "No overlap" = "#BDBDBD"
    ),
    breaks = c("C1 only", "C2 only", "C1 + C2", "No overlap"),
    drop = TRUE
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = function(x) paste0(x * 100, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "Proportion", fill = NULL) +
  theme_classic(base_size = 13)
pdf("bin10kb/diffBin.highFreq.H3K9me3_peakOverlap.pdf", height = 2, width = 3)
print(p)
dev.off()

#### ---- 2.quantify the changed levels of high frequency bins ---- ####
counts <- readRDS("bin10kb/differentialBinAnalysis/counts.rds")
cpms <- readRDS("bin10kb/differentialBinAnalysis/cpms.rds")
condition <- readRDS("bin10kb/differentialBinAnalysis/condition.rds")
cpms.up <- cpms[as.character(up_highFreq$Region),]
cpms.down <- cpms[as.character(down_highFreq$Region),]
condition$group <- paste0(condition$tissue_cellType, "_", condition$age_group)
up.group_means <- t(apply(cpms.up, 1, function(row) {tapply(row, condition$group, mean)}))
down.group_means <- t(apply(cpms.down, 1, function(row) {tapply(row, condition$group, mean)}))

up.values <- colSums(up.group_means)
down.values <- colSums(down.group_means)
up.cpm.FC <- up.values[grep("_Old", names(up.values))] / up.values[grep("_Young", names(up.values))] # change foldchange
down.cpm.FC <- down.values[grep("_Old", names(down.values))] / down.values[grep("_Young", names(down.values))] # change foldchange
xylimit <- 1

plot.data <- data.frame(tissue_cellType = names(up.cpm.FC), up = up.cpm.FC, down = down.cpm.FC[names(up.cpm.FC)])
plot.data$tissue_cellType <- gsub("_Old|_Young", "", plot.data$tissue_cellType)
plot.data <- plot.data %>% left_join(PT.info.H3K27me3.cellType, by = "tissue_cellType")
plot.data <- na.omit(plot.data)
plot.data <- plot.data[order(plot.data$up, decreasing = T),]
saveRDS(plot.data, "bin10kb/diffBin.highFreq.FC.rds")
write.csv(plot.data, "bin10kb/diffBin.highFreq.FC.csv", row.names = F)

#### boxplot: log2FC[CPM(Old/Young)] by tissue (sorted by median)
plot.data$up <- log2(plot.data$up)
plot.data$down <- log2(plot.data$down)
ymax <- max(abs(c(plot.data$up, plot.data$down)), na.rm = TRUE)
ybreaks <- pretty(c(0, ymax), n = 3)
ymax <- max(ybreaks)
## sort tissue by median of up
rank.order <- tapply(plot.data$up, plot.data$tissue, median)
plot.data$tissue <- factor(plot.data$tissue, levels = names(sort(rank.order, decreasing = T)))
p1 <- ggboxplot(plot.data, x = "tissue", y = "up", xlab = "", ylab = "Mean-log2FC[CPM(Old/Young)]", color = "tissue", palette = color.tissue, add = "jitter", add.params = list(size = 0.5)) + ggtitle("Differential up bins (high-frequency)") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + scale_y_continuous(breaks = ybreaks, expand = c(0, 0)) + coord_cartesian(ylim = c(0, ymax), expand = FALSE)
p2 <- ggboxplot(plot.data, x = "tissue", y = "down", xlab = "", ylab = "Mean-log2FC[CPM(Old/Young)]", color = "tissue", palette = color.tissue, add = "jitter", add.params = list(size = 0.5)) + ggtitle("Differential down bins (high-frequency)") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + scale_y_continuous(breaks = -ybreaks, expand = c(0, 0)) + coord_cartesian(ylim = c(-ymax, max(plot.data$down)), expand = FALSE)
## sort category by median of up
rank.order <- tapply(plot.data$up, plot.data$category, median)
plot.data$category <- factor(plot.data$category, levels = names(sort(rank.order, decreasing = T)))
p3 <- ggboxplot(plot.data, x = "category", y = "up", xlab = "", ylab = "Mean-log2FC[CPM(Old/Young)]", color = "category", palette = color.category, add = "jitter", add.params = list(size = 0.5)) + ggtitle("Differential up bins (high-frequency)") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + scale_y_continuous(breaks = ybreaks, expand = c(0, 0)) + coord_cartesian(ylim = c(0, ymax), expand = FALSE)
p4 <- ggboxplot(plot.data, x = "category", y = "down", xlab = "", ylab = "Mean-log2FC[CPM(Old/Young)]", color = "category", palette = color.category, add = "jitter", add.params = list(size = 0.5)) + ggtitle("Differential down bins (high-frequency)") + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + scale_y_continuous(breaks = -ybreaks, expand = c(0, 0)) + coord_cartesian(ylim = c(-ymax, max(plot.data$down)), expand = FALSE)


pdf("bin10kb/diffBin.highFreq.log2FC.boxPlot.pdf", height = 6, width = 6)
ggarrange(p1, p2, p3, p4, nrow = 2, ncol = 2, align = "hv", byrow = T)
dev.off()

