#' @description differential peak analysis for H3K9me3
#' @Author: Zhilin Long
#' @date 2026-6-23

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
source("~/Functions/colorAssign.R")
source("~/Functions/plotDiffOverlap.R")
source("~/Functions/utilities.R")
setwd("~/longzhilin/Projects/AgingMouse/PairedTag/3.histoneLevel/H3K9me3/peak")
PT.info <- readRDS("~/longzhilin/Projects/AgingMouse/PairedTag/1.dataPreprocessing/PT.info.rds")

# summarise by tissue, age and cell type
PT.info.H3K9me3.cellType <- PT.info %>% dplyr::select(tissue_cellType, cellType, tissue, tissue_new, tissue_cellType_new, category, H3K9me3.depth, colors) %>% group_by(tissue_cellType) %>% summarise(cellType = unique(cellType), read = sum(H3K9me3.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), tissue_cellType_new = unique(tissue_cellType_new), category = unique(category), colors = unique(colors), .groups = "drop")
PT.info.RNA.cellType <- PT.info %>% dplyr::select(tissue_cellType, cellType, tissue, tissue_new, tissue_cellType_new, category, RNA.depth, colors) %>% group_by(tissue_cellType) %>% summarise(cellType = unique(cellType), read = sum(RNA.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), tissue_cellType_new = unique(tissue_cellType_new), category = unique(category), colors = unique(colors), .groups = "drop")
PT.info.H3K9me3.cellType_age <- PT.info %>% dplyr::select(tissue_cellType, cellType, age_group, tissue, tissue_new, tissue_cellType_new, category, H3K9me3.depth, colors) %>% group_by(tissue_cellType, age_group) %>% summarise(age_group = unique(age_group), cellType = unique(cellType), read = sum(H3K9me3.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), tissue_cellType_new = unique(tissue_cellType_new), category = unique(category), colors = unique(colors), .groups = "drop")

# load PT-peak
H3K9me3_peak <- read.table("~/longzhilin/public_html/AgingMouse/PairedTag_new/MergeIGV/H3K9me3/H3K9me3_tissue_peak.bed", sep = "\t", header = F, stringsAsFactors = F) # 21666 peak, 0-based
log2FC.threshold <- 0
sig.label <- "fdr"

#### ---- differential analysis ---- ####
dir.create("differentialPeakAnalysis", recursive = TRUE)
bam.files <- list.files("~/longzhilin/public_html/AgingMouse/PairedTag_new", pattern = ".bam$", full.names = T, recursive = T)
bam.files <- bam.files[grep("H3K9me3/cellType_age_replication/bam", bam.files)]

info <- data.frame(label = bam.files)
info$tissue_cellType <- gsub("~/longzhilin/public_html/AgingMouse/PairedTag_new/", "", info$label)
info$tissue_cellType <- gsub("/H3K9me3/cellType_age_replication/bam/H3K9me3.", ":", info$tissue_cellType)
info$sample <- gsub(".bam", "", info$tissue_cellType)
info$batch <- gsub(".*_", "", info$sample)
info$age <- gsub(".*_Old_.*", "Old", info$tissue_cellType)
info$age <- gsub(".*_Young_.*", "Young", info$age)
info$tissue_cellType <- gsub("_Old_.*|_Young_.*", "", info$tissue_cellType)
info <- info %>% left_join(PT.info.H3K9me3.cellType, by = "tissue_cellType")
saveRDS(info, "differentialPeakAnalysis/condition.rds")

#### feature counts
library(Rsubread)
H3K9me3_peak.gr <- makeGRangesFromDataFrame(H3K9me3_peak, keep.extra.columns = T, seqnames.field = "V1", start.field = "V2", end.field = "V3")
H3K9me3_peak.gr$gene_id <- paste0(seqnames(H3K9me3_peak.gr), ":", start(H3K9me3_peak.gr), "-", end(H3K9me3_peak.gr))
peak.SAF <- data.frame(GeneID = H3K9me3_peak.gr$gene_id, Chr = seqnames(H3K9me3_peak.gr), Start = start(H3K9me3_peak.gr), End = end(H3K9me3_peak.gr), Strand = "*")
peak.count <- featureCounts(files = bam.files, annot.inbuilt = "mm10", isGTFAnnotationFile = F, annot.ext = peak.SAF, nthreads = 20, fraction = T, countMultiMappingReads = T)
counts <- peak.count$counts
colnames(counts) <- info$sample
saveRDS(counts, "differentialPeakAnalysis/counts.tissue_peak.rds")
cpms <- edgeR::cpm(counts)
saveRDS(cpms, "differentialPeakAnalysis/cpms.tissue_peak.rds")

#### edgeR analysis
source("~/longzhilin/ProjectCode/AgingProject/edgeR_diffRegions.R")
diffPeak.list <- lapply(unique(info$tissue_cellType), function(tissue_cellType){
  cat("   ", tissue_cellType, "\n")
  info.subset <- info[which(info$tissue_cellType == tissue_cellType),]
  count.subset <- counts[,info.subset$sample]
  cellType <- gsub(".*:", "", tissue_cellType)
  tissue <- gsub(":.*", "", tissue_cellType)
  if(length(unique(info.subset$batch)) > 1){
    batch <- info.subset$batch
  }else{
    batch <- NULL
  }
  outDir <- file.path("differentialPeakAnalysis", tissue)
  if(!dir.exists(outDir)){
    dir.create(outDir, recursive = TRUE)
  }
  group <- factor(info.subset$age, levels = c("Young", "Old"))
  title <- paste0(tissue_cellType, "(depth[", unique(info.subset$read), "M])")
  diffPeaks <- edgeR.diffRegions(countTables = count.subset, group = group, batch = batch, outPath = outDir, type = cellType, title = title, logFC = 1, PValue = 0.05, ContinuousComparison = FALSE, cpm.filter = TRUE, min.regions = 10)
  res <- diffPeaks$GroupComparison.Old.vs.Young
  return(res)
})
names(diffPeak.list) <- unique(info$tissue_cellType)
saveRDS(diffPeak.list, "differentialPeakAnalysis/diffTissuePeak.list.rds")

#### ---- counting & plot differential Peaks ---- ####
H3K9me3_peak.gr$peakLength <- width(H3K9me3_peak.gr)
selected.peaks <- H3K9me3_peak.gr$gene_id[which(H3K9me3_peak.gr$peakLength >= 200000)] # >= 200 kb; 913 peaks
color.chr <- ggsci::pal_igv("default")(26)
names(color.chr) <- c(paste0("chr", 1:22), "chrM", "chrX", "chrY", "chrUn")
diffPeak.list <- readRDS("differentialPeakAnalysis/diffTissuePeak.list.rds")

#### The celltype frequency of differential Peaks
source("~/longzhilin/AnalysisCode/IDConvert.R")
library(ChIPseeker)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
diffPeaks <- diffPeak.list
groups <- c("All", names(color.category))
groups <- groups[grep("Neural|Others|Germline", groups, invert = T)]
pdf("differentialPeakAnalysis/diffTissuePeakFrequency.pdf", width = 15, height = 8)
diffPeaks.Freq.list <- lapply(groups, function(group){
  cat("Analysis group:", group, "\n")
  if(group == "All"){
    cellTypes <- PT.info.H3K9me3.cellType$tissue_cellType
    tissues <- unique(as.character(PT.info.H3K9me3.cellType$tissue))
  }else{
    cellTypes <- PT.info.H3K9me3.cellType$tissue_cellType[PT.info.H3K9me3.cellType$category == group]
    tissues <- unique(as.character(PT.info.H3K9me3.cellType$tissue[PT.info.H3K9me3.cellType$category == group]))
  }
  diffPeaks.subset <- diffPeaks[which(names(diffPeaks) %in% cellTypes)]

  ## The frequency of up-regulated peaks
  diffPeaks.cellType.up <- lapply(diffPeaks.subset, function(x){
    x <- x[which(x$Region %in% selected.peaks),]
    x <- x$Region[x$logFC > log2FC.threshold & x[,sig.label] < 0.05]
    return(x)
  })
  diffPeaks.cellType.up <- unlist(diffPeaks.cellType.up)
  diffPeaks.cellType.up.Freq <- as.data.frame(table(diffPeaks.cellType.up))
  diffPeaks.cellType.up.Freq <- diffPeaks.cellType.up.Freq[order(diffPeaks.cellType.up.Freq$Freq, decreasing = T),]
  plot.data <- diffPeaks.cellType.up.Freq
  plot.data$category_var <- as.character(plot.data$Freq)
  plot.data <- as.data.frame(table(plot.data$category_var))
  plot.data$Var1 <- factor(plot.data$Var1, levels = as.character(sort(as.numeric(levels(plot.data$Var1)))))
  p1 <- ggbarplot(plot.data, x = "Var1", y = "Freq", ylab = "Peak number", xlab = "Cell type", label = T, fill = "grey") + ggtitle("Differential-up Peaks(Old/Young)") + theme(plot.title = element_text(hjust = 0.5))

  diffPeaks.tissue.up <- data.frame(label = names(diffPeaks.cellType.up), Region = diffPeaks.cellType.up)
  rownames(diffPeaks.tissue.up) <- NULL
  diffPeaks.tissue.up$tissue <- gsub(":.*", "", diffPeaks.tissue.up$label)
  diffPeaks.tissue.up <- unique(diffPeaks.tissue.up[,c("tissue", "Region")])
  diffPeaks.tissue.up.Freq <- as.data.frame(table(diffPeaks.tissue.up$Region))
  diffPeaks.tissue.up.Freq <- diffPeaks.tissue.up.Freq[order(diffPeaks.tissue.up.Freq$Freq, decreasing = T),]

  ## The frequency of down-regulated peaks
  diffPeaks.cellType.down <- lapply(diffPeaks.subset, function(x){
    x <- x[which(x$Region %in% selected.peaks),]
    x <- x$Region[x$logFC < (log2FC.threshold * -1) & x[,sig.label] < 0.05]
    return(x)
  })
  diffPeaks.cellType.down <- unlist(diffPeaks.cellType.down)
  diffPeaks.cellType.down.Freq <- as.data.frame(table(diffPeaks.cellType.down))
  diffPeaks.cellType.down.Freq <- diffPeaks.cellType.down.Freq[order(diffPeaks.cellType.down.Freq$Freq, decreasing = T),]
  plot.data <- diffPeaks.cellType.down.Freq
  plot.data$category_var <- as.character(plot.data$Freq)
  plot.data <- as.data.frame(table(plot.data$category_var))
  plot.data$Var1 <- factor(plot.data$Var1, levels = as.character(sort(as.numeric(levels(plot.data$Var1)))))
  p2 <- ggbarplot(plot.data, x = "Var1", y = "Freq", ylab = "Peak number", xlab = "Cell type", label = T, fill = "grey") + ggtitle("Differential-down Peaks(Old/Young)") + theme(plot.title = element_text(hjust = 0.5))
  if(group == "All"){
    p <- ggarrange(p1, p2, ncol = 1, nrow = 2)
  }else{
    p <- ggarrange(p1, NULL, p2, NULL, ncol = 2, nrow = 2)
  }
  annotated_combined <- annotate_figure(p, top = text_grob(group, color = "black", face = "bold", size = 16))
  print(annotated_combined)

  diffPeaks.tissue.down <- data.frame(label = names(diffPeaks.cellType.down), Region = diffPeaks.cellType.down)
  rownames(diffPeaks.tissue.down) <- NULL
  diffPeaks.tissue.down$tissue <- gsub(":.*", "", diffPeaks.tissue.down$label)
  diffPeaks.tissue.down <- unique(diffPeaks.tissue.down[,c("tissue", "Region")])
  diffPeaks.tissue.down.Freq <- as.data.frame(table(diffPeaks.tissue.down$Region))
  diffPeaks.tissue.down.Freq <- diffPeaks.tissue.down.Freq[order(diffPeaks.tissue.down.Freq$Freq, decreasing = T),]

  ## merge to output
  colnames(diffPeaks.cellType.up.Freq) <- c("Peak", "n_sig_cellType_up")
  colnames(diffPeaks.cellType.down.Freq) <- c("Peak", "n_sig_cellType_down")
  colnames(diffPeaks.tissue.up.Freq) <- c("Peak", "sig_cellType_up_cover_n_tissue")
  colnames(diffPeaks.tissue.down.Freq) <- c("Peak", "sig_cellType_down_cover_n_tissue")
  df_list <- list(diffPeaks.cellType.up.Freq, diffPeaks.cellType.down.Freq, diffPeaks.tissue.up.Freq,diffPeaks.tissue.down.Freq)
  res <- Reduce(function(x,y) merge(x,y,all=T,by = "Peak"), df_list)
  res[is.na(res)] <- 0
  res$n_sig_cellType_total <- rowSums(res[,2:3], na.rm = T)
  res <- res[order(res$n_sig_cellType_total, decreasing = T),]
  res$sig_cellType_cover_total_tissue <- rowSums(res[,4:5], na.rm = T)
  res$n_cellType <- length(cellTypes)
  res$n_tissue <- length(tissues)

  ## annotation peak
  peak.gr <- GRanges(seqnames = gsub(":.*", "", res$Peak), ranges = IRanges(start = as.numeric(gsub("-.*", "", gsub(".*:", "", res$Peak))), end = as.numeric(gsub(".*-", "", res$Peak))), strand = "*")
  peak.anno <- annotatePeak(peak.gr, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Mm.eg.db")
  peak.anno <- as.data.frame(peak.anno)
  peak.anno$annotation <- gsub(" \\(.*", "", peak.anno$annotation)
  peak.anno$Peak <- paste0(peak.anno$seqnames, ":", peak.anno$start, "-", peak.anno$end)
  res <- res %>% left_join(peak.anno[,c("Peak", "geneId", "annotation")], by = "Peak")
  gene_symbol <- IDConvert(peak.anno$geneId, from = "ENTREZID", to = "SYMBOL", species = "Mus musculus", method = "clusterProfiler")
  colnames(gene_symbol)[1] <- "geneId"
  res <- res %>% left_join(gene_symbol, by = "geneId")
  res <- res[,c("Peak", "geneId", "annotation", "SYMBOL", "n_sig_cellType_up", "n_sig_cellType_down", "n_sig_cellType_total", "n_cellType", "sig_cellType_up_cover_n_tissue", "sig_cellType_down_cover_n_tissue", "sig_cellType_cover_total_tissue", "n_tissue")]
  return(res)
})
dev.off()
names(diffPeaks.Freq.list) <- groups
saveRDS(diffPeaks.Freq.list, "differentialPeakAnalysis/diffTissuePeaks.Freq.rds")
writexl::write_xlsx(diffPeaks.Freq.list, path = "differentialPeakAnalysis/diffTissuePeaks.Freq.xlsx")

## plot the frequency of different peak
# Count overlapping points at each (x, y) coordinate for bubble plot
plot.data <- diffPeaks.Freq.list$All
plot.data.count <- plot.data %>% group_by(n_sig_cellType_up, n_sig_cellType_down) %>% summarise(n = n(), .groups = "drop") %>% mutate(log10_n = log10(n))
top_peaks <- plot.data %>% slice_max(n_sig_cellType_up, n = 5) %>% bind_rows(plot.data %>% slice_max(n_sig_cellType_down, n = 5)) %>% distinct(Peak)
xy_limit <- max(plot.data[,c("n_sig_cellType_up", "n_sig_cellType_down")], na.rm = T)
points_to_label <- subset(plot.data, Peak %in% top_peaks$Peak)
p <- ggplot(plot.data.count, aes(x = n_sig_cellType_up, y = n_sig_cellType_down)) +
  geom_count(aes(size = log10_n, color = log10_n), alpha = 0.7) +
  scale_size_continuous(name = "log10 (# of peaks)", range = c(1, 3)) +
  geom_text_repel(data = points_to_label, 
                  aes(x = n_sig_cellType_up, y = n_sig_cellType_down, label = Peak), 
                  size = 3, box.padding = 0.5, point.padding = 0.3, 
                  force = 2, max.overlaps = 20) +
  scale_color_viridis_c(option = "D", name = "log10 (# of peaks)") + 
  xlim(0, xy_limit) + ylim(0, xy_limit) +
  ggtitle(paste0("Total differential peaks: ", nrow(plot.data))) + 
  xlab("# of cell types with increased H3K9me3") + 
  ylab("# of cell types with decreased H3K9me3") +
  theme_classic() + 
  theme(
    aspect.ratio = 1, 
    plot.title = element_text(hjust = 0.5)
  )
pdf("differentialPeakAnalysis/differentialPeaks_bubblePlot.pdf", width = 5, height = 5)
print(p)
dev.off()

#### ---- plot log2FC(CPM) heatmap for differential Peak ---- ####
selected.tissue_cellType <- PT.info.H3K9me3.cellType_age %>% group_by(tissue_cellType) %>% filter(all(read >= 0.5)) %>% ungroup() %>% dplyr::select(tissue_cellType, tissue) %>% distinct() # 100 cell types
plot.data <- as.data.frame(table(selected.tissue_cellType$tissue))
ggbarplot(plot.data, x = "Var1", y = "Freq", fill = "Var1", color = "Var1", palette = color.tissue, ylab = "Number of cell types", xlab = "") + geom_text(aes(label = Freq), angle = 90, hjust = -0.25, vjust = 0.5, size = 3.5, color = "black") + theme(legend.position = "none", axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))

#### extract differential peaks
diffPeak.list <- readRDS("differentialPeakAnalysis/diffTissuePeak.list.rds")
info <- readRDS("differentialPeakAnalysis/condition.rds")
info$group <- paste0(info$tissue_cellType, "_", info$age)
diffPeak.sig <- lapply(diffPeak.list, function(y){
  y <- y[which(y$Region %in% selected.peaks),]
  up.peaks <- y$Region[y$logFC > log2FC.threshold & y[,sig.label] < 0.05]
  down.peaks <- y$Region[y$logFC < (log2FC.threshold * -1) & y[,sig.label] < 0.05]
  diffPeak <- unique(c(up.peaks, down.peaks))
  return(unique(diffPeak))
})
diffPeak.sig <- unique(unlist(diffPeak.sig)) # 890 peaks
filtered_type <- ".tissue"

#### calculate log2FC(CPM-Old-mean/CPM-Young-mean)
cpms.peak <- readRDS("differentialPeakAnalysis/cpms.tissue_peak.rds")
cpms.mean.peak <- t(apply(cpms.peak, 1, function(x) {
  tapply(x, info$group, mean)
}))
cpms_Old <- cpms.mean.peak[,grep("Old", colnames(cpms.mean.peak))] + 0.01
cpms_Young <- cpms.mean.peak[,grep("Young", colnames(cpms.mean.peak))] + 0.01
cpms.fc <- log2(cpms_Old / cpms_Young)
colnames(cpms.fc) <- gsub("_Old|_Young", "", colnames(cpms.fc))
saveRDS(cpms.fc, paste0("differentialPeakAnalysis/peak_cpms.fc", filtered_type, ".rds"))

#### heatmap configure
type <- "cpm" # cpm; edgeR
heatmap.name <- "CPM:log2(Old/Young)" # CPM:log2(Old/Young); edgeR:logFC
plot.data <- cpms.fc[diffPeak.sig, selected.tissue_cellType$tissue_cellType]
# plot.data <- edgeR.logFC[, selected.tissue_cellType$tissue_cellType]
peak.info <- data.frame(peak = rownames(plot.data), chr = gsub(":.*", "", rownames(plot.data)))
peak.info.gr <- GRanges(seqnames = gsub(":.*", "", peak.info$peak), ranges = IRanges(start = as.numeric(gsub("-.*", "", gsub(".*:", "", peak.info$peak))), end = as.numeric(gsub(".*-", "", peak.info$peak))), strand = "*")
column_info <- PT.info.H3K9me3.cellType[match(colnames(plot.data), PT.info.H3K9me3.cellType$tissue_cellType),]

#### peak clustering with hierarchical clustering (euclidean + ward.D2) ####
cluster.method <- "hc_cluster"
X <- plot.data
X[is.na(X)] <- 0
row_dist <- dist(X, method = "euclidean")
row_hc <- hclust(row_dist, method = "ward.D2")
row_clusters <- cutree(row_hc, k = 2)

# The cluster with fewer peaks is defined as C1, and the other cluster as C2.
raw_cl <- row_clusters[row_hc$order]
c1_cluster <- names(which.min(table(raw_cl))) # raw cluster id with fewer peaks

# Explicitly put C1 (fewer peaks) rows first, keeping the dendrogram order within each cluster.
# (Do NOT rely on rev(as.dendrogram()) + as.hclust() to flip the tree: it does not reliably
# reverse $order, which is why C1 could still end up at the bottom.)
ord <- row_hc$order
ord <- c(ord[as.character(raw_cl) == c1_cluster], ord[as.character(raw_cl) != c1_cluster])

peak_order <- data.frame(peak = row_hc$labels[ord], cluster = ifelse(as.character(row_clusters[ord]) == c1_cluster, "C1", "C2"))
row_split <- factor(peak_order$cluster, levels = c("C1", "C2"))
peak.info <- peak.info[match(peak_order$peak, peak.info$peak),]
plot.data <- plot.data[peak_order$peak,]

# rank columns by mean log2FC of C1
rank_cluster <- plot.data[peak_order$cluster == "C1", , drop = FALSE]
rank_cluster[is.na(rank_cluster)] <- 0
rank_order <- colMeans(rank_cluster, na.rm = T)
plot.data <- plot.data[,order(rank_order)]

## row/column annotation (built after final row/column order to stay aligned)
column_info <- column_info[match(colnames(plot.data), column_info$tissue_cellType),]
column_ha <- HeatmapAnnotation(category = column_info$category, tissue = column_info$tissue_new, col = list(category = color.category, tissue = color.tissue2))

dot.colors <- PT.info.H3K9me3.cellType$colors
names(dot.colors) <- PT.info.H3K9me3.cellType$tissue_cellType
dot.colors <- dot.colors[colnames(plot.data)]
point_legend <- Legend(labels = names(dot.colors), title = "cellType", type = "points", legend_gp = gpar(col = dot.colors, pch = 21, cex = 2), ncol = 20)
bottom_annotation <- HeatmapAnnotation(cellType = anno_points(rep(1, ncol(plot.data)), annotation_height = unit(1, "cm"), gp = gpar(col = dot.colors, cex = 0.8)), show_annotation_name = TRUE)

scores <- as.numeric(plot.data)
col_fun <- colorRamp2(c(quantile(scores, 0.01), 0, quantile(scores, 0.99)), c("blue", "white", "red"))
pdf(paste0("peak_", type, "_", cluster.method, filtered_type, ".pdf"), height = 30, width = 60)
## show cluster
res <- X[,names(sort(rank_order))]
res_info <- column_info[match(colnames(res), column_info$tissue_cellType),]
res_ha <- HeatmapAnnotation(category = res_info$category, tissue = res_info$tissue, col = list(category = color.category, tissue = color.tissue))
left_info <- peak.info[match(rownames(res), peak.info$peak),]
p <- Heatmap(res, col = col_fun, name = heatmap.name, column_title = paste0("The changes of H3K9me3 signal in peaks(", nrow(plot.data), " peaks): age-read >= 0.5M"), top_annotation = res_ha, bottom_annotation = bottom_annotation, show_column_dend = F, show_row_dend = T, show_row_names = F, show_column_names = F, cluster_columns = F, column_title_gp = gpar(fontsize = 14), width = unit(15, "cm"), height = unit(6, "cm"), column_names_rot = 60, gap = unit(0.5, "mm"), row_split = row_clusters, cluster_row_slices = F)
p <- draw(p, annotation_legend_list = list(point_legend), annotation_legend_side = "bottom")

## cluster row for each cluster
p <- Heatmap(plot.data, col = col_fun, name = heatmap.name, column_title = paste0("The changes of H3K9me3 signal in peaks(", nrow(plot.data), " peaks): age-read >= 0.5M"), top_annotation = column_ha, bottom_annotation = bottom_annotation, show_column_dend = F, show_row_dend = F, show_row_names = F, show_column_names = F, cluster_columns = F, column_title_gp = gpar(fontsize = 14), width = unit(15, "cm"), height = unit(6, "cm"), column_names_rot = 60, gap = unit(0.5, "mm"), cluster_rows = T, cluster_row_slices = F, row_split = row_split, clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2")
p <- draw(p, annotation_legend_list = list(point_legend), annotation_legend_side = "bottom")

## column split by category
p1 <- Heatmap(plot.data, col = col_fun, name = heatmap.name, column_title = paste0("The changes of H3K9me3 signal in peaks(", nrow(plot.data), " peaks): age-read >= 0.5M"), top_annotation = column_ha, bottom_annotation = bottom_annotation, show_column_dend = F, show_row_dend = F, show_row_names = F, cluster_columns = F, column_title_gp = gpar(fontsize = 14), width = unit(15, "cm"), column_split = column_info$category, height = unit(6, "cm"), column_names_rot = 90, gap = unit(0.5, "mm"), cluster_rows = T, cluster_row_slices = F, row_split = row_split, clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2")
p1 <- draw(p1, annotation_legend_list = list(point_legend), annotation_legend_side = "bottom")
dev.off()

heatmap_row_order <- row_order(p1)
peak_order <- lapply(names(heatmap_row_order), function(x){
  res <- heatmap_row_order[[x]]
  res <- data.frame(peak = rownames(plot.data)[res], cluster = x)
  return(res)
})
peak_order <- do.call(rbind, peak_order)
heatmap.order <- list(row_order = peak_order, column_order = sort(rank_order))
saveRDS(heatmap.order, paste0("heatmap.order", filtered_type, ".final.rds"))

#### optimal heatmap
scores <- as.numeric(unique(unlist(plot.data)))
# col_fun <- colorRamp2(c(quantile(scores, 0.01), 0, quantile(scores, 0.99)), c("blue", "white", "red"))
col_fun <- colorRamp2(c(quantile(scores, 0.01), 0, quantile(scores, 0.99)), c("#002fff", "white", "#f94433"))
pdf(paste0("peak_", type, "_", cluster.method, filtered_type, ".final2.pdf"), height = 20, width = 25)
## plot the selected cell types
plot.data <- plot.data[heatmap.order$row_order$peak, names(heatmap.order$column_order)]
column_labels <- seq_len(ncol(plot.data))
row_split <- heatmap.order$row_order$cluster
clusters <- unique(heatmap.order$row_order$cluster)
cluster.color <- ggsci::pal_npg("nrc")(length(clusters)) 
names(cluster.color) <- clusters
column_info <- column_info[match(colnames(plot.data), column_info$tissue_cellType),]
column_ha <- HeatmapAnnotation(category = column_info$category, tissue = column_info$tissue_new, col = list(category = color.category, tissue = color.tissue2))
p <- Heatmap(plot.data, col = col_fun, name = heatmap.name, column_title = paste0("H3K9me3 signal in ", nrow(plot.data), " peaks"), top_annotation = column_ha, show_column_dend = F, show_row_dend = F, show_row_names = F, cluster_columns = F, column_title_gp = gpar(fontsize = 14), width = unit(8, "cm"), height = unit(6, "cm"), column_names_rot = 90, cluster_rows = F, cluster_row_slices = F, row_split = row_split, row_gap = unit(0.5, "mm"), column_gap = unit(0.5, "mm"), show_column_names = TRUE, column_names_side = "bottom", column_labels = column_labels, column_names_gp = gpar(fontsize = 6))
p <- draw(p)

## column split by category
column_order_new <- PT.info.H3K9me3.cellType[match(names(heatmap.order$column_order), PT.info.H3K9me3.cellType$tissue_cellType),]
column_order_new <- arrange(column_order_new, category)
idx <- match(column_order_new$tissue_cellType, names(heatmap.order$column_order))
column_labels <- column_labels[idx]
plot.data <- plot.data[heatmap.order$row_order$peak, column_order_new$tissue_cellType]
dot.colors <- column_order_new$colors
names(dot.colors) <- column_order_new$cellType
cellType.colors <- column_order_new$colors
names(cellType.colors) <- column_order_new$cellType
point_legend <- Legend(labels = names(dot.colors), title = "cellType", type = "points", legend_gp = gpar(col = dot.colors, pch = 21, cex = 1), ncol = 4)
column_ha <- HeatmapAnnotation(category = column_order_new$category, tissue = column_order_new$tissue_new, col = list(category = color.category, tissue = color.tissue2))
bottom_annotation <- HeatmapAnnotation(cellType_new = column_order_new$cellType, col = list(cellType_new = cellType.colors), cellType = anno_points(rep(1, ncol(plot.data)), annotation_height = unit(0.25, "cm"), gp = gpar(col = dot.colors, cex = 0.5)), show_annotation_name = TRUE)
p1 <- Heatmap(plot.data, col = col_fun, name = heatmap.name, column_title = paste0("H3K9me3 signal in ", nrow(plot.data), " peaks"), top_annotation = column_ha, show_column_dend = F, show_row_dend = F, show_row_names = F, bottom_annotation = bottom_annotation, cluster_columns = F, column_title_gp = gpar(fontsize = 14), width = unit(10, "cm"), column_split = column_order_new$category, height = unit(6, "cm"), column_names_rot = 90, cluster_rows = F, cluster_row_slices = F, row_split = row_split, row_gap = unit(0.5, "mm"), column_gap = unit(0.5, "mm"), show_column_names = TRUE, column_names_side = "bottom", column_labels = column_labels, column_names_gp = gpar(fontsize = 6))
p1 <- draw(p1, annotation_legend_list = list(point_legend), annotation_legend_side = "right")

## plot the all cell types
plot.data.raw <- cpms.fc[heatmap.order$row_order$peak, ]
rank_cluster <- plot.data.raw[heatmap.order$row_order$cluster == "C1", , drop = FALSE]
rank_cluster[is.na(rank_cluster)] <- 0
rank_order <- colMeans(rank_cluster, na.rm = T)
plot.data.raw <- plot.data.raw[,order(rank_order)]
column_order_new <- PT.info.H3K9me3.cellType[match(colnames(plot.data.raw), PT.info.H3K9me3.cellType$tissue_cellType),]
column_order_new <- arrange(column_order_new, category)
plot.data.raw <- plot.data.raw[heatmap.order$row_order$peak, column_order_new$tissue_cellType]
dot.colors <- column_order_new$colors
names(dot.colors) <- column_order_new$cellType
cellType.colors <- column_order_new$colors
names(cellType.colors) <- column_order_new$cellType
point_legend <- Legend(labels = names(dot.colors), title = "cellType", type = "points", legend_gp = gpar(col = dot.colors, pch = 21, cex = 1), ncol = 4)
column_ha <- HeatmapAnnotation(category = column_order_new$category, tissue = column_order_new$tissue_new, col = list(category = color.category, tissue = color.tissue2))
bottom_annotation <- HeatmapAnnotation(cellType_new = column_order_new$cellType, col = list(cellType_new = cellType.colors), cellType = anno_points(rep(1, ncol(plot.data.raw)), annotation_height = unit(0.25, "cm"), gp = gpar(col = dot.colors, cex = 0.5)), show_annotation_name = TRUE)
p1 <- Heatmap(plot.data.raw, col = col_fun, name = heatmap.name, column_title = paste0("H3K9me3 signal in ", nrow(plot.data.raw), " peaks"), top_annotation = column_ha, show_column_dend = F, show_row_dend = F, show_row_names = F, show_column_names = F, bottom_annotation = bottom_annotation, cluster_columns = F, column_title_gp = gpar(fontsize = 14), width = unit(12, "cm"), column_split = column_order_new$category, height = unit(8, "cm"), column_names_rot = 90, cluster_rows = F, cluster_row_slices = F, row_split = row_split, row_gap = unit(0.5, "mm"), column_gap = unit(0.5, "mm"))
p1 <- draw(p1, annotation_legend_list = list(point_legend), annotation_legend_side = "right")
dev.off()

peak.gr <- GRanges(seqnames = gsub(":.*", "", heatmap.order$row_order$peak), ranges = IRanges(start = as.numeric(gsub("-.*", "", gsub(".*:", "", heatmap.order$row_order$peak))), end = as.numeric(gsub(".*-", "", heatmap.order$row_order$peak))), strand = "*")
peak.gr.list <- split(peak.gr, heatmap.order$row_order$cluster)
write.table(as.data.frame(peak.gr.list$C1)[,1:3], "~/longzhilin/public_html/AgingMouse/PairedTag_new/MergeIGV/H3K9me3/H3K9me3_tissue_peak_C1.bed", sep = "\t", quote = F, row.names = F, col.names = F)
write.table(as.data.frame(peak.gr.list$C2)[,1:3], "~/longzhilin/public_html/AgingMouse/PairedTag_new/MergeIGV/H3K9me3/H3K9me3_tissue_peak_C2.bed", sep = "\t", quote = F, row.names = F, col.names = F)

#### plot boxplot
peak_cluster <- heatmap.order$row_order
plot.data <- lapply(unique(peak_cluster$cluster), function(i){
  mean_km <- cpms.fc[peak_cluster$peak[which(peak_cluster$cluster == i)],selected.tissue_cellType$tissue_cellType]
  mean_km <- colMeans(mean_km, na.rm = T)
  mean_km <- data.frame(mean = mean_km, tissue_cellType = names(mean_km))
  colnames(mean_km)[1] <- paste0(i, "_mean")
  return(mean_km)
})
plot.data <- Reduce(function(x,y) merge(x,y,all=T,by = "tissue_cellType"), plot.data)
plot.data <- plot.data %>% left_join(PT.info.H3K9me3.cellType, by = "tissue_cellType")
order_data <- plot.data %>% group_by(tissue) %>% summarise(median_val = median(C1_mean, na.rm = TRUE)) %>% arrange(median_val)
plot.data$tissue <- factor(plot.data$tissue, levels = order_data$tissue)
p1 <- ggboxplot(plot.data, x = "tissue", y = "C1_mean", xlab = "", ylab = "Mean-log2FC[CPM(Old/Young)]", color = "tissue", palette = color.tissue, outlier.shape = NA) + ggtitle("H3K9me3 signal in C1 peaks") + theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1), plot.title = element_text(hjust = 0.5))
order_data <- plot.data %>% group_by(tissue) %>% summarise(median_val = median(C2_mean, na.rm = TRUE)) %>% arrange(median_val)
plot.data$tissue <- factor(plot.data$tissue, levels = order_data$tissue)
p2 <- ggboxplot(plot.data, x = "tissue", y = "C2_mean", xlab = "", ylab = "Mean-log2FC[CPM(Old/Young)]", color = "tissue", palette = color.tissue, outlier.shape = NA) + ggtitle("H3K9me3 signal in C2 peaks") + theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1), plot.title = element_text(hjust = 0.5))
order_data <- plot.data %>% group_by(category) %>% summarise(median_val = median(C1_mean, na.rm = TRUE)) %>% arrange(median_val)
plot.data$category <- factor(plot.data$category, levels = order_data$category)
p11 <- ggboxplot(plot.data, x = "category", y = "C1_mean", xlab = "", ylab = "Mean-log2FC[CPM(Old/Young)]", color = "category", palette = color.category, outlier.shape = NA) + ylim(c(-0.7,0.25)) + ggtitle("H3K9me3 signal in C1 peaks") + theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1), plot.title = element_text(hjust = 0.5)) # note the ylim
order_data <- plot.data %>% group_by(category) %>% summarise(median_val = median(C2_mean, na.rm = TRUE)) %>% arrange(median_val)
plot.data$category <- factor(plot.data$category, levels = order_data$category)
p22 <- ggboxplot(plot.data, x = "category", y = "C2_mean", xlab = "", ylab = "Mean-log2FC[CPM(Old/Young)]", color = "category", palette = color.category, outlier.shape = NA) + ggtitle("H3K9me3 signal in C2 peaks") + theme(legend.position = "none", axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1), plot.title = element_text(hjust = 0.5))
pdf(paste0("peak_cpm_hc_cluster", filtered_type, ".mean.pdf"), height = 5, width = 8)
p <- ggarrange(p1, p2, p11, p22, ncol = 2, nrow = 2, align = "hv")
print(p)
print(ggarrange(p22, NULL, NULL, NULL, ncol = 2, nrow = 2, align = "hv", heights = c(2,1)))
dev.off()

#### output the mean-logFC for each peak
peak.cpms.fc <- as.data.frame(cpms.fc[heatmap.order$row_order$peak, names(sort(heatmap.order$column_order))])
peak.cpms.fc$meanLog2FC <- rowMeans(peak.cpms.fc)
write.csv(peak.cpms.fc, paste0("peak_cpm_hc_cluster", filtered_type, ".log2FC.csv"), row.names = T)

cluster.mean.log2FC <- lapply(unique(heatmap.order$row_order$cluster), function(i){
  peaks <- heatmap.order$row_order$peak[heatmap.order$row_order$cluster == i]
  cluster.log2FC <- as.data.frame(cpms.fc[peaks, ])
  cluster.mean.log2FC <- colMeans(cluster.log2FC, na.rm = T)
  cluster.mean.log2FC <- data.frame(tissue_cellType = names(cluster.mean.log2FC), meanLog2FC = cluster.mean.log2FC)
  colnames(cluster.mean.log2FC)[2] <- paste0(i, "_", colnames(cluster.mean.log2FC)[2])
  return(cluster.mean.log2FC)
})
cluster.mean.log2FC <- merge(cluster.mean.log2FC[[1]], cluster.mean.log2FC[[2]], by = "tissue_cellType")
cluster.mean.log2FC <- cluster.mean.log2FC[order(cluster.mean.log2FC$C1_meanLog2FC),]
saveRDS(cluster.mean.log2FC, paste0("peak_cpm_hc_cluster", filtered_type, ".meanLog2FC.rds"))

#### ---- Pathway enrichment ---- ####
heatmap.order <- readRDS("heatmap.order.tissue.final.rds")
peak.gr <- GRanges(seqnames = gsub(":.*", "", heatmap.order$row_order$peak), ranges = IRanges(start = as.numeric(gsub("-.*", "", gsub(".*:", "", heatmap.order$row_order$peak))), end = as.numeric(gsub(".*-", "", heatmap.order$row_order$peak))), strand = "*")
peak.gr.list <- split(peak.gr, heatmap.order$row_order$cluster)
#### rGREAT method
source("~/longzhilin/AnalysisCode/Epigenomics/Great.enrichment.R")
source("~/longzhilin/AnalysisCode/PathwayEnrichment/pathway.dotplot.R")
enrichedPathways <- lapply(peak.gr.list, function(x){
  p <- Great.enrichment(gr = x, bg = NULL, model = "local", species = "Mus musculus", local.genesets = "GO:BP", use.MsigDB = F, plot = F, title = NULL, simplify.GO = F)
  return(p)
})
writexl::write_xlsx(enrichedPathways, "peak.cluster.enrichedPathways.xlsx")
p <- pathway.dotplot(enrichedResult.lists = enrichedPathways, title = "Enriched pathways for H3K9me3 peaks", items = c("description", "fold_enrichment", "p_adjust"), threshold = 0.05, top.n = 30, GeneRatio = F, max.dot.size = 4, sub.words = 80, y.text.size = 11)
print(p)

#### ---- explore the correlation between gene expression/change and H3K9me3 change of C1 peaks ---- ####
heatmap.order <- readRDS("heatmap.order.tissue.final.rds")
C1.Mean_log2FC <- data.frame(tissue_cellType = names(heatmap.order$column_order), Mean_log2FC = as.numeric(heatmap.order$column_order))
C1.Mean_log2FC <- C1.Mean_log2FC %>% left_join(PT.info.H3K9me3.cellType, by = "tissue_cellType")

#### correlation between edgeR-FC(Old/Young) and H3K9me3 change of C1 peaks
RNA.edgeR.list <- readRDS("~/longzhilin/Projects/AgingMouse/PairedTag/2.rnaLevel/differentialGeneAnalysis/RNA.edgeR.list.rds")
RNA.edgeR.logFC <- lapply(names(RNA.edgeR.list), function(tissue){
  edgeR.FC <- RNA.edgeR.list[[tissue]]
  logFC <- lapply(names(edgeR.FC), function(cellType){
    res <- edgeR.FC[[cellType]]
    res$tissue_cellType <- paste0(tissue, ":", gsub(" |-", "_", cellType))
    return(res[,c("tissue_cellType", "gene", "logFC")])
  })
  logFC <- do.call(rbind, logFC)
  return(logFC)
})
RNA.edgeR.logFC <- do.call(rbind, RNA.edgeR.logFC)
RNA.edgeR.logFC <- RNA.edgeR.logFC %>% pivot_wider(names_from = gene, values_from = logFC, values_fill = NA)
RNA.edgeR.logFC <- as.data.frame(RNA.edgeR.logFC)
idx <- match(C1.Mean_log2FC$tissue_cellType, RNA.edgeR.logFC$tissue_cellType)
RNA.edgeR.logFC <- RNA.edgeR.logFC[na.omit(idx),]
rownames(RNA.edgeR.logFC) <- RNA.edgeR.logFC$tissue_cellType
RNA.edgeR.logFC <- RNA.edgeR.logFC[,-1]
RNA.edgeR.logFC[is.na(RNA.edgeR.logFC)] <- 0

edgeRlogFC.cor <- apply(RNA.edgeR.logFC, 2, function(x){
  res.raw <- data.frame(RNA_logFC = x, H3K9me3_log2FC = C1.Mean_log2FC$Mean_log2FC[which(!is.na(idx))])
  res <- na.omit(res.raw)
  cat("Matrix size:", nrow(res), "\n")
  if(nrow(res) > nrow(res.raw)/3){
    up.cor <- cor.test(res$RNA_logFC, res$H3K9me3_log2FC, method = "spearman")
    res <- data.frame(RNA_edgeRlogFC.H3K9me3_C1.cor = up.cor$estimate, RNA_edgeRlogFC.H3K9me3_C1.p = up.cor$p.value)
  }else{
    res <- data.frame(RNA_edgeRlogFC.H3K9me3_C1.cor = NA, RNA_edgeRlogFC.H3K9me3_C1.p = NA)
  }
  return(res)
})
edgeRlogFC.cor <- do.call(rbind, edgeRlogFC.cor)
edgeRlogFC.cor <- na.omit(edgeRlogFC.cor)
edgeRlogFC.cor$gene <- rownames(edgeRlogFC.cor)
edgeRlogFC.cor <- edgeRlogFC.cor[order(edgeRlogFC.cor$RNA_edgeRlogFC.H3K9me3_C1.cor, decreasing = T),]
edgeRlogFC.cor$RNA_edgeRlogFC.H3K9me3_C1.fdr <- p.adjust(edgeRlogFC.cor$RNA_edgeRlogFC.H3K9me3_C1.p, method = "fdr")
edgeRlogFC.cor <- edgeRlogFC.cor[, c("gene", "RNA_edgeRlogFC.H3K9me3_C1.cor", "RNA_edgeRlogFC.H3K9me3_C1.p", "RNA_edgeRlogFC.H3K9me3_C1.fdr")]

#### scatter plot
genes <- c("Pcgf5", "H2-Q7", "H2-Q6", "H2-Q5", "H2-Q4", "Psmb9", "Ranbp17", "Col16a1", "Eln", "Col6a3", "Lama2")
cor_col <- "RNA_edgeRlogFC.H3K9me3_C1.cor"
fdr_col <- "RNA_edgeRlogFC.H3K9me3_C1.fdr"
plot.data <- edgeRlogFC.cor[, c("gene", cor_col, fdr_col)]
plot.data <- plot.data[order(plot.data[[cor_col]], na.last = NA), ]
plot.data$Rank <- seq_len(nrow(plot.data))
point_labels <- subset(plot.data, gene %in% genes)
point_labels$direction <- ifelse(
  point_labels[[cor_col]] > 0, "Positive",
  ifelse(point_labels[[cor_col]] < 0, "Negative", "Zero")
)
n <- nrow(plot.data)
pdf("CorrelationWithGene/H3K9me3_C1.GeneExpression.correlation.scatter.pdf",
    width = 6, height = 6)
ggplot(plot.data, aes(x = Rank, y = .data[[cor_col]])) +
  geom_point(color = "grey75", size = 1, alpha = 0.8) +
  geom_point(
    data = point_labels,
    aes(color = direction),
    size = 2
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = point_labels,
    aes(label = gene, color = direction),
    size = 2,
    segment.size = 0.5,
    point.padding = 0.3,
    max.overlaps = 100,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name = "Correlation",
    values = c(
      Negative = "#0859FD",
      Positive = "#d30808",
      Zero = "grey50"
    )
  ) +
  scale_x_continuous(
    breaks = c(1, round((n + 1) / 2), n),
    labels = c("Most negative", "Intermediate", "Most positive")
  ) +
  labs(
    x = "Genes ranked by Spearman correlation",
    y = "Spearman correlation",
    title = "Spearman correlation between gene expression and H3K9me3-C1"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

#### pathway enrichment
source("~/Functions/clusterProfiler.GSEA.R")
source("~/Functions/clusterProfiler.enricher.R")
source("~/Functions/pathway.dotplot.R")
## GSEA method: rank by sign(cor) × -log10(p)
GSEA.rank <- edgeRlogFC.cor
p_min_nonzero <- min(GSEA.rank$RNA_edgeRlogFC.H3K9me3_C1.p[is.finite(GSEA.rank$RNA_edgeRlogFC.H3K9me3_C1.p) & GSEA.rank$RNA_edgeRlogFC.H3K9me3_C1.p > 0], na.rm = TRUE)
GSEA.rank$RNA_edgeRlogFC.H3K9me3_C1.p[which(GSEA.rank$RNA_edgeRlogFC.H3K9me3_C1.p == 0)] <- p_min_nonzero * 0.01
H3K9me3.gene.change.list <- list(RNA_CPM.H3K9me3_C1 = sign(GSEA.rank$RNA_CPM.H3K9me3_C1.cor) * -log10(GSEA.rank$RNA_CPM.H3K9me3_C1.p), RNA_log2FC.H3K9me3_C1 = sign(GSEA.rank$RNA_log2FC.H3K9me3_C1.cor) * -log10(GSEA.rank$RNA_log2FC.H3K9me3_C1.p), RNA_edgeRlogFC.H3K9me3_C1 = sign(GSEA.rank$RNA_edgeRlogFC.H3K9me3_C1.cor) * -log10(GSEA.rank$RNA_edgeRlogFC.H3K9me3_C1.p))
GSEA.enrichedPathways <- lapply(H3K9me3.gene.change.list, function(x){
  names(x) <- GSEA.rank$gene
  geneRank <- x[!is.na(x)]
  GSEA.res <- clusterProfiler.GSEA(geneRank, geneType = "SYMBOL", db.type = "GO", GO.ont = "BP", pvalueCutoff = 1, species = "Mus musculus")
  GSEA.res.up <- GSEA.res %>% filter(NES > 0) %>% arrange(p.adjust)
  GSEA.res.down <- GSEA.res %>% filter(NES < 0) %>% arrange(p.adjust)
  GSEA.res <- list(positive = GSEA.res.up, negative = GSEA.res.down)
  return(GSEA.res)
})
names(GSEA.enrichedPathways) <- names(H3K9me3.gene.change.list)
GSEA.enrichedPathways <- unlist(GSEA.enrichedPathways, recursive = F)
p <- pathway.dotplot(enrichedResult.lists = GSEA.enrichedPathways[1:2], title = "", items = c("Description", "setSize", "p.adjust"), threshold = 0.05, top.n = 15, GeneRatio = F, max.dot.size = 4, sub.words = 80, y.text.size = 11)
print(p)

Correlation.enrichedPathways <- list(Hypergeometric = Hypergeometric.enrichedPathways, GSEA = GSEA.enrichedPathways)
saveRDS(Correlation.enrichedPathways, "CorrelationWithGene/H3K9me3_C1.GeneExpression.enrichedPathways.rds")
