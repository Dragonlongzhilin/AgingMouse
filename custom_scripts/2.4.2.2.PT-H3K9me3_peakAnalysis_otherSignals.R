#' @description H3K9me3 peak for other signals
#' @Author: Zhilin Long
#' @date 2026-6-23

library(Seurat)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(GenomicRanges)
library(rtracklayer)
library(ggrepel)
library(Rsubread)
library(cowplot)
library(gridExtra)
set.seed(101)
source("~/longzhilin/ProjectCode/AgingMouse/Functions/colorAssign.R")
source("~/longzhilin/ProjectCode/AgingProject/plotDiffOverlap.R")
source("~/longzhilin/ProjectCode/AgingMouse/Functions/utilities.R")
setwd("~/longzhilin/Projects/AgingMouse/PairedTag/3.histoneLevel/H3K9me3/peak")
PT.info <- readRDS("~/longzhilin/Projects/AgingMouse/PairedTag/1.dataPreprocessing/PT.info.rds")

# summarise by tissue, age and cell type
PT.info.H3K9me3.cellType <- PT.info %>% select(tissue_cellType, cellType, tissue, tissue_new, category, H3K9me3.depth) %>% group_by(tissue_cellType) %>% summarise(cellType = unique(cellType), read = sum(H3K9me3.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), category = unique(category), .groups = "drop")

# load peaks
H3K9me3_peak <- read.table("~/longzhilin/public_html/AgingMouse/PairedTag_new/MergeIGV/H3K9me3/H3K9me3_tissue_peak.bed", sep = "\t", header = F, stringsAsFactors = F)
H3K9me3_peak.gr <- makeGRangesFromDataFrame(H3K9me3_peak, keep.extra.columns = T, seqnames.field = "V1", start.field = "V2", end.field = "V3")
H3K9me3_peak.gr$gene_id <- paste0(seqnames(H3K9me3_peak.gr), ":", start(H3K9me3_peak.gr), "-", end(H3K9me3_peak.gr))
peak.SAF <- data.frame(GeneID = H3K9me3_peak.gr$gene_id, Chr = seqnames(H3K9me3_peak.gr), Start = start(H3K9me3_peak.gr), End = end(H3K9me3_peak.gr), Strand = "*")

other.histones <- names(color.histone)[grep("H3K9me3|ATAC", names(color.histone), invert = T)]
#### ---- quantify the signal of other histones in peaks ---- ####
dir.create("otherSignals", showWarnings = F, recursive = T)
for(histone in other.histones){
  #### construct the bam files
  cellType.bam.files <- lapply(names(color.tissue), function(tissue){
    bam.files <- list.files("~/longzhilin/public_html/AgingMouse/PairedTag_new", pattern = ".bam$", full.names = T, recursive = T)
    bam.files <- bam.files[grep(paste0("/", tissue, "/"), bam.files)]
    bam.files <- bam.files[grep(histone, bam.files)]
    bam.files <- bam.files[grep("cellType_age_replication/bam", bam.files)]
    file.names <- gsub(".bam", "", basename(bam.files))
    file.names <- gsub(paste0(histone, "."), "", basename(file.names))
    file.names <- paste0(tissue, ":", file.names)
    info <- data.frame(fileName = file.names, tissue = tissue, bamFile = bam.files)
    return(info)
  })
  cellType.bam.files <- do.call(rbind, cellType.bam.files)

  counts <- featureCounts(files = cellType.bam.files$bamFile, annot.inbuilt = "mm10", isGTFAnnotationFile = F, annot.ext = peak.SAF, nthreads = 20, fraction = T, countMultiMappingReads = T)
  counts <- counts$counts
  colnames(counts) <- cellType.bam.files$fileName
  saveRDS(counts, paste0("otherSignals/", histone, ".tissue.counts.rds"))
  cpms <- edgeR::cpm(counts)
  saveRDS(cpms, paste0("otherSignals/", histone, ".tissue.cpms.rds"))
}

#### ---- plot the heatmap of histone-log2FC(cpm[Old/Young]) ---- ####
color.chr <- ggsci::pal_igv("default")(26)
names(color.chr) <- c(paste0("chr", 1:22), "chrM", "chrX", "chrY", "chrUn")
heatmap.order <- readRDS("heatmap.order.tissue.final.rds")
clusters <- unique(heatmap.order$row_order$cluster)
cluster.color <- ggsci::pal_npg("nrc")(length(clusters))
names(cluster.color) <- clusters

other.histones <- c("H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3")
cluster1_mean.histone.list <- lapply(other.histones, function(histone){
  cat("Analyzing ", histone, "\n")
  # summarise by tissue, age and cell type
  idx <- grep(paste0(histone, ".depth"), colnames(PT.info))
  PT.info.subset <- PT.info[,c(23,5,26,21,idx)]
  colnames(PT.info.subset)[5] <- "depth"
  PT.info.cellType <- PT.info.subset %>% group_by(tissue_cellType) %>% summarise(cellType = unique(cellType), read = sum(depth, na.rm = T), tissue = unique(tissue_new), category = unique(category), .groups = "drop")

  cpms <- readRDS(paste0("otherSignals/", histone, ".tissue.cpms.rds"))
  cpms.peak <- cpms[grep("chr", rownames(cpms)),]
  groups <- gsub("_mouse.*", "", colnames(cpms.peak))
  cpms.mean.peak <- t(apply(cpms.peak, 1, function(x) {
    tapply(x, groups, mean)
  }))
  cpms_Old <- cpms.mean.peak[,grep("Old", colnames(cpms.mean.peak))] + 0.01
  cpms_Young <- cpms.mean.peak[,grep("Young", colnames(cpms.mean.peak))] + 0.01
  cpms.fc <- log2(cpms_Old / cpms_Young)
  colnames(cpms.fc) <- gsub("_Old|_Young", "", colnames(cpms.fc))
  saveRDS(cpms.fc, paste0("otherSignals/", histone, ".peak_cpms.fc.rds"))

  plot.data <- cpms.fc
  peak.info <- heatmap.order$row_order[na.omit(match(rownames(plot.data), heatmap.order$row_order$peak)),]
  peak.info$chr <- gsub(":.*", "", peak.info$peak)

  column_info <- data.frame(tissue_cellType = colnames(plot.data))
  column_info <- column_info %>% left_join(PT.info.cellType, by = "tissue_cellType")
  # column_info <- subset(column_info, read >= 2) # remove the cell type with less than 2M reads
  plot.data <- plot.data[,column_info$tissue_cellType]

  pdf(paste0("otherSignals/", histone,".CPM_log2FC.pdf"), height = 15, width = 40)
  ## rank by mean log2FC of total peak of cluster 1
  column_order <- sort(heatmap.order$column_order)
  idx <- match(names(column_order), column_info$tissue_cellType)
  column_info.subset <- column_info[na.omit(idx),]
  plot.data.subset <- plot.data[,na.omit(idx)]
  row_order <- heatmap.order$row_order
  # row_order$cluster <- factor(row_order$cluster, levels = c("C1", "C3", "C2"))
  # row_order <- arrange(row_order, cluster)
  peak.info.subset <- peak.info[match(row_order$peak, peak.info$peak),]
  plot.data.subset <- plot.data.subset[row_order$peak,]
  scores <- as.numeric(unique(plot.data.subset))
  col_fun <- colorRamp2(c(quantile(scores, 0.01), 0, quantile(scores, 0.99)), c("blue", "white", "red"))
  H3K9me3_col_fun <- colorRamp2(c(quantile(column_order, 0.01), 0, quantile(column_order, 0.99)), c("blue", "white", "red"))
  column_ha <- HeatmapAnnotation(category = column_info.subset$category, tissue = column_info.subset$tissue, H3K9me3_change = column_order[which(!is.na(idx))], col = list(category = color.category, tissue = color.tissue2, H3K9me3_change = H3K9me3_col_fun))
  left_annotation <- rowAnnotation(chr = peak.info.subset$chr, cluster = peak.info.subset$cluster, col = list(chr = color.chr[unique(peak.info.subset$chr)], cluster = cluster.color))
  p <- Heatmap(plot.data.subset, col = col_fun, name = "CPM:log2(Old/Young)", column_title = paste0("The changes of ", histone," signal in H3K9me3 peaks(", nrow(plot.data.subset), " peaks)"), top_annotation = column_ha, show_column_dend = F, show_row_dend = F, show_row_names = F, cluster_rows = F, left_annotation = left_annotation, cluster_columns = F, column_title_gp = gpar(fontsize = 14), width = unit(60, "cm"), row_split = peak.info.subset$cluster, height = unit(15, "cm"), column_names_rot = 60, cluster_row_slices = FALSE, row_gap = unit(2, "mm"))
  draw(p)
  p <- Heatmap(plot.data.subset, col = col_fun, name = "CPM:log2(Old/Young)", column_title = paste0("The changes of ", histone," signal in H3K9me3 peaks(", nrow(plot.data.subset), " peaks)"), top_annotation = column_ha, show_column_dend = F, show_row_dend = F, show_row_names = F, cluster_rows = F, left_annotation = left_annotation, cluster_columns = F, column_title_gp = gpar(fontsize = 14), width = unit(60, "cm"), row_split = peak.info.subset$cluster, column_split = column_info.subset$category, height = unit(15, "cm"), column_names_rot = 60, cluster_row_slices = FALSE, row_gap = unit(2, "mm"))
  draw(p)

  dev.off()
  return(p)
})

#### optimal heatmap
heatmap.order <- readRDS("heatmap.order.tissue.final.rds")
clusters <- unique(heatmap.order$row_order$cluster)
cluster.color <- ggsci::pal_npg("nrc")(length(clusters))
names(cluster.color) <- clusters
pdf("otherSignals/CPM_log2FC.tissue.heatmap.pdf", height = 12, width = 12)
cluster1_mean.histone.list <- lapply(other.histones, function(histone){
  cat("Analyzing ", histone, "\n")
  # summarise by tissue, age and cell type
  idx <- grep(paste0(histone, ".depth"), colnames(PT.info))
  PT.info.subset <- PT.info[,c(23,5,26,21,idx)]
  colnames(PT.info.subset)[5] <- "depth"
  PT.info.cellType <- PT.info.subset %>% group_by(tissue_cellType) %>% summarise(cellType = unique(cellType), read = sum(depth, na.rm = T), tissue = unique(tissue_new), category = unique(category), .groups = "drop")

  cpms <- readRDS(paste0("otherSignals/", histone, ".tissue.cpms.rds"))
  cpms.peak <- cpms[grep("chr", rownames(cpms)),]
  groups <- gsub("_mouse.*", "", colnames(cpms.peak))
  cpms.mean.peak <- t(apply(cpms.peak, 1, function(x) {
    tapply(x, groups, mean)
  }))
  cpms_Old <- cpms.mean.peak[,grep("Old", colnames(cpms.mean.peak))] + 0.01
  cpms_Young <- cpms.mean.peak[,grep("Young", colnames(cpms.mean.peak))] + 0.01
  cpms.fc <- log2(cpms_Old / cpms_Young)
  colnames(cpms.fc) <- gsub("_Old|_Young", "", colnames(cpms.fc))
  saveRDS(cpms.fc, paste0("otherSignals/", histone, ".peak_cpms.fc.rds"))

  plot.data <- cpms.fc
  row_order <- heatmap.order$row_order
  peak.info <- row_order[na.omit(match(rownames(plot.data), row_order$peak)),]
  peak.info$chr <- gsub(":.*", "", peak.info$peak)

  column_info <- data.frame(tissue_cellType = colnames(plot.data))
  column_info <- column_info %>% left_join(PT.info.cellType, by = "tissue_cellType")
  # column_info <- subset(column_info, read >= 2) # remove the cell type with less than 2M reads
  plot.data <- plot.data[,column_info$tissue_cellType]

  ## rank by mean H3K9me3 log2FC of total peak of cluster 1
  column_order <- sort(heatmap.order$column_order)
  idx <- match(names(column_order), column_info$tissue_cellType)
  column_info.subset <- column_info[na.omit(idx),]
  plot.data.subset <- plot.data[,na.omit(idx)]
  peak.info.subset <- peak.info[match(row_order$peak, peak.info$peak),]
  plot.data.subset <- plot.data.subset[row_order$peak,]
  scores <- as.numeric(unique(plot.data.subset))
  if(histone != "H3K27me3"){
    # col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
    col_fun <- colorRamp2(c(quantile(scores, 0.01), 0, quantile(scores, 0.99)), c("#002fff", "white", "#f94433"))
  }else{
    col_fun <- colorRamp2(c(-5, 0, 5), c("blue", "white", "red"))
  }
  column_ha <- HeatmapAnnotation(category = column_info.subset$category, tissue = column_info.subset$tissue, col = list(category = color.category, tissue = color.tissue2))
  left_annotation <- rowAnnotation(cluster = peak.info.subset$cluster, col = list(cluster = cluster.color))
  p <- Heatmap(plot.data.subset, col = col_fun, name = "CPM:log2(Old/Young)", column_title = paste0(histone," signal in H3K9me3 peaks(", nrow(plot.data.subset), " peaks)"), top_annotation = column_ha, show_column_dend = F, show_row_dend = F, show_row_names = F, cluster_rows = F, left_annotation = left_annotation, cluster_columns = F, show_column_names = F, column_title_gp = gpar(fontsize = 14), width = unit(8, "cm"), row_split = peak.info.subset$cluster, height = unit(6, "cm"), column_names_rot = 60, cluster_row_slices = FALSE, row_gap = unit(0.5, "mm"), column_gap = unit(0.5, "mm"))
  draw(p)
  p <- Heatmap(plot.data.subset, col = col_fun, name = "CPM:log2(Old/Young)", column_title = paste0(histone," signal in H3K9me3 peaks(", nrow(plot.data.subset), " peaks)"), top_annotation = column_ha, show_column_dend = F, show_row_dend = F, show_row_names = F, cluster_rows = F, left_annotation = left_annotation, cluster_columns = F, show_column_names = F, column_title_gp = gpar(fontsize = 14), width = unit(8, "cm"), row_split = peak.info.subset$cluster, column_split = column_info.subset$category, height = unit(6, "cm"), column_names_rot = 60, cluster_row_slices = FALSE, row_gap = unit(0.5, "mm"), column_gap = unit(0.5, "mm"))
  draw(p)

  # calculate the mean log2FC of C1 and C2
  peak_cluster <- peak.info.subset
  cluster.mean.log2FC <- lapply(unique(peak_cluster$cluster), function(i){
    mean_km <- cpms.fc[peak_cluster$peak[which(peak_cluster$cluster == i)],column_info.subset$tissue_cellType]
    mean_km <- colMeans(mean_km, na.rm = T)
    mean_km <- data.frame(mean = mean_km, tissue_cellType = names(mean_km))
    colnames(mean_km)[1] <- paste0(i, "_mean")
    return(mean_km)
  })
  cluster.mean.log2FC <- Reduce(function(x,y) merge(x,y,all=T,by = "tissue_cellType"), cluster.mean.log2FC)
  cluster.mean.log2FC <- cluster.mean.log2FC %>% left_join(PT.info.H3K9me3.cellType, by = "tissue_cellType")
  return(cluster.mean.log2FC)
})
dev.off()
names(cluster1_mean.histone.list) <- other.histones
saveRDS(cluster1_mean.histone.list, "otherSignals/cluster1_mean.histone.tissue.list.rds")
