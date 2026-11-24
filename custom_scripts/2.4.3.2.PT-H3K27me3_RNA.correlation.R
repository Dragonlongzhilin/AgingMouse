#' @description correlation between H3K27me3 and gene expression
#' @Author: Zhilin Long
#' @date 2026-8-25

library(Seurat)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(gridExtra)
library(ggrepel)
library(GenomicRanges)
library(GenomeInfoDb)
set.seed(101)
source("~/Functions/colorAssign.R")
source("~/Functions/utilities.R")
setwd("~/longzhilin/Projects/AgingMouse/PairedTag")

histone <- "H3K27me3"
RNA.cpm.mats <- readRDS("1.dataPreprocessing/RNA.cpm.mats.rds")
PT.info <- readRDS("1.dataPreprocessing/PT.info.rds")
PT.info.H3K27me3.cellType <- PT.info %>% select(tissue_cellType, tissue_cellType_new, cellType, tissue, tissue_new, category, H3K27me3.depth) %>% group_by(tissue_cellType) %>% summarise(tissue_cellType_new = unique(tissue_cellType_new), cellType = unique(cellType), read = sum(H3K27me3.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), category = unique(category), .groups = "drop")
PT.info.H3K27me3.tissue <- PT.info %>% select(tissue, H3K27me3.depth, tissue_new) %>% group_by(tissue) %>% summarise(read = sum(H3K27me3.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), .groups = "drop")
PT.info.H3K27me3.cellType_age <- PT.info %>% dplyr::select(tissue_cellType, tissue_cellType_new, cellType, age_group, tissue, tissue_new, category, H3K27me3.depth, colors) %>% group_by(tissue_cellType, age_group) %>% summarise(tissue_cellType_new = unique(tissue_cellType_new), age_group = unique(age_group), cellType = unique(cellType), read = sum(H3K27me3.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), category = unique(category), colors = unique(colors), .groups = "drop")

PT.info.RNA.cellType <- PT.info %>% dplyr::select(tissue_cellType, tissue_cellType_new, cellType, tissue, tissue_new, category, RNA.depth) %>% group_by(tissue_cellType) %>% summarise(tissue_cellType_new = unique(tissue_cellType_new), cellType = unique(cellType), read = sum(RNA.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), category = unique(category), .groups = "drop")
PT.info.RNA.cellType_age <- PT.info %>% dplyr::select(tissue_cellType, tissue_cellType_new, cellType, age_group, tissue, tissue_new, category, RNA.depth, colors) %>% group_by(tissue_cellType, age_group) %>% summarise(tissue_cellType_new = unique(tissue_cellType_new), age_group = unique(age_group), cellType = unique(cellType), read = sum(RNA.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), category = unique(category), colors = unique(colors), .groups = "drop")

## select cell types with at least 0.5M reads in each age group
selected.tissue_cellType <- PT.info.H3K27me3.cellType_age %>% group_by(tissue_cellType) %>% filter(all(read >= 0.5)) %>% ungroup() %>% dplyr::select(tissue_cellType, tissue) %>% distinct() # 115 cell types

#### ---- quantify the histone signal in TSS+/- 5kb ---- ####
library(rtracklayer)
options(scipen = 999)
gencodeGTF <- import("~/longzhilin/Softwares/genomes/mouse/mm10/gencode.vM25.basic.annotation.gtf") # 1-based
genomics_region <- read.table("~/longzhilin/Softwares/genomes/mouse/mm10/genomics_region.bed", header = F, sep = "\t", stringsAsFactors = F) # 0-base
transcripts <- gencodeGTF[gencodeGTF$type == "transcript", ]
tss_pos <- ifelse(strand(transcripts) == "+", start(transcripts), end(transcripts))
tss_granges <- GRanges(
  seqnames = seqnames(transcripts),
  ranges = IRanges(start = tss_pos, width = 1),
  strand = strand(transcripts),
  source = transcripts$source,
  type = "exon",
  score = transcripts$score,
  phase = transcripts$phase,
  gene_id = transcripts$gene_id,
  transcript_id = transcripts$transcript_id,
  gene_name = transcripts$gene_name,
  exon_number = transcripts$exon_number,
  exon_id = transcripts$exon_id
)
tss_extended_5kb <- promoters(tss_granges, upstream = 5000, downstream = 5000)
optimized_fix_bounds <- function(regions, genomicRange) {
  chr_max <- setNames(genomicRange$V3, genomicRange$V1)
  current_chrs <- as.character(seqnames(regions))
  max_lens <- chr_max[current_chrs]
  start_fixed <- pmax(1, start(regions))
  end_fixed <- pmin(max_lens, end(regions))
  idx <- which(is.na(end_fixed))
  end_fixed[idx] <- end(regions)[idx]
  start_fixed <- pmin(start_fixed, end_fixed)
  ranges(regions) <- IRanges(start = start_fixed, end = end_fixed)
  return(regions)
}
tss_extended_5kb <- optimized_fix_bounds(tss_extended_5kb, genomics_region)
export(tss_extended_5kb, "~/longzhilin/Softwares/genomes/mouse/mm10/mm10.gencode_tss_extended_5kb.gtf", format = "GTF")

## run featureCounts
cellType.bam.files <- lapply(names(color.tissue), function(tissue){
  bam.files <- list.files("~/longzhilin/public_html/AgingMouse/PairedTag_new", pattern = ".bam$", full.names = T, recursive = T)
  bam.files <- bam.files[grep(paste0("/", tissue, "/"), bam.files)]
  bam.files <- bam.files[grep(histone, bam.files)]
  bam.files <- bam.files[grep("cellType_age_replication/bam", bam.files)]
  file.names <- gsub(".bam", "", basename(bam.files))
  file.names <- gsub(paste0(histone, "."), "", file.names)
  file.names <- paste0(tissue, ":", file.names)
  info <- data.frame(fileName = file.names, tissue = tissue, bamFile = bam.files)
  return(info)
})
cellType.bam.files <- do.call(rbind, cellType.bam.files)

TSS_5kb_gtf <- "~/longzhilin/Softwares/genomes/mouse/mm10/mm10.gencode_tss_extended_5kb.gtf"
outDir <- "~/longzhilin/Projects/AgingMouse/PairedTag/3.histoneLevel/H3K27me3/TSS_5kb"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
library(Rsubread)
counts <- featureCounts(files = cellType.bam.files$bamFile, annot.inbuilt = "mm10", isGTFAnnotationFile = T, annot.ext = TSS_5kb_gtf, nthreads = 20, fraction = T, countMultiMappingReads = T, GTF.attrType = "gene_name")
counts <- counts$counts
colnames(counts) <- cellType.bam.files$fileName
saveRDS(counts, "3.histoneLevel/H3K27me3/TSS_5kb/counts.rds")

#### ---- edgeR analysis for H3K27me3 signal (TSS+/-5kb)---- ####
source("~/longzhilin/ProjectCode/AgingProject/edgeR_gene.R")
library(edgeR)
library(EnhancedVolcano)

counts_bin10kb <- readRDS("3.histoneLevel/H3K27me3/bin10kb/differentialBinAnalysis/counts.rds")
condition <- readRDS("3.histoneLevel/H3K27me3/bin10kb/differentialBinAnalysis/condition.rds")
DNA_depth <- colSums(counts_bin10kb)
counts_TSS_5kb <- readRDS("3.histoneLevel/H3K27me3/TSS_5kb/counts.rds")
counts_TSS_5kb <- counts_TSS_5kb[,colnames(counts_bin10kb)]
edgeR_figures_outdir <- "3.histoneLevel/H3K27me3/TSS_5kb/edgeR_figures"
dir.create(edgeR_figures_outdir, recursive = TRUE, showWarnings = FALSE)
TSS_5kb.edgeR.list <- lapply(unique(condition$tissue_cellType), function(x){
  cat("Analysis ", x, "...\n")
  idx <- which(condition$tissue_cellType == x)
  counts_TSS_5kb.subset <- counts_TSS_5kb[,idx]
  condition.subset <- condition[idx,]
  DNA_depth.subset <- DNA_depth[condition.subset$label]
  cpms_TSS_5kb.subset <- edgeR::cpm(counts_TSS_5kb.subset, lib.size = DNA_depth.subset) 
  group <- factor(condition.subset$age, levels = c("Young", "Old"))

  if(length(unique(condition.subset$batch)) > 1){
    batch <- factor(condition.subset$batch)
    design <- model.matrix(~0+group+batch)
    lvls <- levels(factor(group))
    colnames(design)[1:length(lvls)] <- lvls
  }else{
    batch <- NULL
    design <- model.matrix(~0+group)
    colnames(design) <- levels(factor(group))
  }
  cat("Batch: ", batch, "\n")

  myContrasts <- "Old-Young"
  contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"), as.list(myContrasts),levels=list(design))))
  obj <- DGEList(counts = counts_TSS_5kb.subset, group = group, lib.size = DNA_depth.subset)

  cat("Filtering method: cpm >= 1 in at least 2 samples\n")
  count_gt <- apply(cpms_TSS_5kb.subset, 1, function(x) sum(x >= 1, na.rm = TRUE))
  keep <- which(count_gt >= 2)
  genes.left <- length(keep)
  cat(nrow(counts_TSS_5kb.subset), " genes in raw matrix\n")
  cat(genes.left, " genes in filted matrix\n")
  if(genes.left > 10){
    obj <- obj[keep,, keep.lib.sizes=T]
    obj <- calcNormFactors(obj)
    obj <- estimateDisp(obj, design, robust=TRUE)
    fit_tag <- glmFit(obj, design)
    lrt <- glmLRT(fit_tag, contrast = contrast.matrix[,1])
    fdr <- p.adjust(lrt$table$PValue, method = "fdr")
    out <- cbind(cpms_TSS_5kb.subset[rownames(lrt$table),], lrt$table, fdr)
    out <- as.data.frame(out)
    out$gene <- rownames(lrt$table)
    out <- out[order(out$fdr),]

    pdf(paste0(edgeR_figures_outdir, "/", x, ".pdf"), width = 6, height = 6)
    p <- edgeR_plot(edgeR_result = out, PValue = T, threshold = 0.05, top.n = 5, size = 3.5, title = x)
    print(p)
    dev.off()
  }else{
    cat("    less than 10 genes in ", x, "\n")
    out <- NULL
  }
  return(out)
})
names(TSS_5kb.edgeR.list) <- unique(condition$tissue_cellType)
saveRDS(TSS_5kb.edgeR.list, "3.histoneLevel/H3K27me3/TSS_5kb/TSS_5kb.edgeR.list.rds")

#### ---- The correlation between H3K27me3 change and Gene expression ---- ####
type.colors <- c("#E64B35", "#00468B", "#D1D1D1")
names(type.colors) <- c("up", "down", "ns")
log2FC.threshold <- log2(1.2)
sig.label <- "PValue"

TSS_5kb.edgeR.list <- readRDS("3.histoneLevel/H3K27me3/TSS_5kb/TSS_5kb.edgeR.list.rds")
RNA.edgeR.list <- readRDS("2.rnaLevel/differentialGeneAnalysis/RNA.edgeR.list.rds")
RNA.edgeR.list <- unlist(RNA.edgeR.list, recursive = F)
names(RNA.edgeR.list) <- gsub("\\.", ":", names(RNA.edgeR.list))
tissue_cellTypes <- intersect(names(TSS_5kb.edgeR.list), names(RNA.edgeR.list))
tissues <- unique(gsub(":.*", "", tissue_cellTypes))

## TSS extended 5kb overlapping with H3K27me3 peak
tss_extended_5kb <- rtracklayer::import("~/longzhilin/Softwares/genomes/mouse/mm10/mm10.gencode_tss_extended_5kb.gtf")
H3K27me3_peak <- readRDS("3.histoneLevel/H3K27me3/peak/H3K27me3_tissue_peak.rds") # 13954 peaks
idx <- findOverlaps(H3K27me3_peak, tss_extended_5kb)
tss_extended_5kb$type <- "ns"
tss_extended_5kb$type[subjectHits(idx)] <- "OverlapPeak"
tss_extended_5kb <- as.data.frame(tss_extended_5kb)
gene_peak <- unique(tss_extended_5kb[,c("gene_name", "type")])
colnames(gene_peak) <- c("gene", "class")
gene_peak_unique <- gene_peak %>% group_by(gene) %>% summarise(class = ifelse("OverlapPeak" %in% class, "OverlapPeak", "ns")) %>% ungroup()

####  Correlation between H3K27me3 change and Gene expression
pdf("3.histoneLevel/H3K27me3/TSS_5kb/H3K27me3_expression.correlation.pdf", width = 30, height = 20)
plot.correlation <- lapply(tissues, function(tissue){
  cat("Analysis ", tissue, "...\n")
  tissue_cellType_subset <- tissue_cellTypes[grep(paste0(tissue, ":"), tissue_cellTypes)]
  results <- lapply(tissue_cellType_subset, function(tissue_cellType){
    cellType <- gsub(".*:", "", tissue_cellType)
    cat("   Analysis ", cellType, "...\n")
    TSS.data <- TSS_5kb.edgeR.list[[tissue_cellType]]
    RNA.data <- RNA.edgeR.list[[tissue_cellType]]

    RNA.data$type <- "ns"
    RNA.data$type[which(RNA.data$fdr < 0.05 & RNA.data$logFC > 0)] <- "up"
    RNA.data$type[which(RNA.data$fdr < 0.05 & RNA.data$logFC < 0)] <- "down"
    RNA.data$type <- factor(RNA.data$type, levels = c("up", "down", "ns"))

    TSS.data$type <- "ns"
    TSS.data$type[which(TSS.data[,sig.label] < 0.05 & TSS.data$logFC > log2FC.threshold)] <- "up"
    TSS.data$type[which(TSS.data[,sig.label] < 0.05 & TSS.data$logFC < -log2FC.threshold)] <- "down"
    TSS.data$type <- factor(TSS.data$type, levels = c("up", "down", "ns"))
    TSS.data <- TSS.data %>% left_join(gene_peak_unique, by = "gene")

    plot.data <- merge(RNA.data[,c("gene", "logFC", "fdr", "type")], TSS.data[,c("gene", "logFC", sig.label, "type", "class")], by = "gene", suffixes = c(".RNA", ".TSS"))
    plot.data <- plot.data[order(plot.data[,3]),]
    plot.data <- plot.data[which(plot.data$type.TSS != "ns"),]
    # plot.data <- plot.data[which(plot.data$type.TSS != "ns" & plot.data$class == "OverlapPeak"),]

    label.up <- if(any(plot.data$type.RNA == "up")) head(plot.data$gene[plot.data$type.RNA == "up"], 5)
    label.down <- if(any(plot.data$type.RNA == "down")) head(plot.data$gene[plot.data$type.RNA == "down"], 5)
    if(is.null(label.up) & is.null(label.down)){
      show.label <- NULL
    }else{
      show.label <- plot.data$gene
    }

    # split into 4 quadrants
    plot.data.subset <- plot.data[which(plot.data$type.RNA != "ns" & plot.data$type.TSS != "ns"),]
    Q1 <- sum(plot.data.subset$logFC.RNA > 0 & plot.data.subset$logFC.TSS > 0) # gene up, H3K27me3 up
    Q2 <- sum(plot.data.subset$logFC.RNA > 0 & plot.data.subset$logFC.TSS < 0) # gene up, H3K27me3 down
    Q3 <- sum(plot.data.subset$logFC.RNA < 0 & plot.data.subset$logFC.TSS < 0) # gene down, H3K27me3 down
    Q4 <- sum(plot.data.subset$logFC.RNA < 0 & plot.data.subset$logFC.TSS > 0) # gene down, H3K27me3 up
    total <- nrow(plot.data.subset)
    labels <- c(
      sprintf("n=%d (%.1f%%)", Q1, Q1/total*100),
      sprintf("n=%d (%.1f%%)", Q2, Q2/total*100),
      sprintf("n=%d (%.1f%%)", Q3, Q3/total*100),
      sprintf("n=%d (%.1f%%)", Q4, Q4/total*100)
    )

    # Fisher test
    data <- matrix(c(Q2, Q1, Q3, Q4), nrow = 2, ncol = 2, byrow = FALSE)
    res <- fisher.test(data)
    subtitle <- paste0("Fisher's Exact Test: ", p_to_stars(res$p.value))
    info <- data.frame(tissue_cellType = tissue_cellType, Gene.up_H3K27me3.up = Q1, Gene.up_H3K27me3.down = Q2, Gene.down_H3K27me3.down = Q3, Gene.down_H3K27me3.up = Q4, total = total, p.value = res$p.value)

    p <- ggscatter(plot.data, x = "logFC.RNA", y = "logFC.TSS", xlab = "RNA: logfoldchanges", ylab = "H3K27me3: logfoldchanges", color = "type.RNA", palette = type.colors, size = 0.5, label = show.label, repel = T, label.select = unique(c(label.up, label.down)), font.label = c(10, "plain")) + geom_hline(yintercept = 0, lty = "dashed", color = "black") + geom_vline(xintercept = 0, lty = "dashed", color = "black") + annotate("text", x = c(max(plot.data$logFC.RNA)*0.8, max(plot.data$logFC.RNA)*0.8, min(plot.data$logFC.RNA)*0.5, min(plot.data$logFC.RNA)*0.5), y = c(max(plot.data$logFC.TSS)*0.8, min(plot.data$logFC.TSS)*0.8, min(plot.data$logFC.TSS)*0.8, max(plot.data$logFC.TSS)*0.8), label = labels, size = c(3.5, 4, 3.5, 4), color = c("black", "red", "black", "red")) + ggtitle(cellType) + labs(subtitle = subtitle) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5, size = 11), axis.text = element_text(size = 10), axis.title = element_text(size = 11), plot.subtitle = element_text(size = 10))
    return(list(info = info, p = p))
  })
  plots <- lapply(results, function(x){x$p})
  p <- patchwork::wrap_plots(plots, ncol = 6, guides = "collect") + patchwork::plot_annotation(title = tissue, theme = theme(plot.title = element_text(size = 12, hjust = 0.5), legend.position = "top"))
  print(p)

  info <- lapply(results, function(x){x$info})
  info <- do.call(rbind, info)
  return(info)
})
dev.off()
plot.correlation <- do.call(rbind, plot.correlation)
plot.correlation <- plot.correlation[order(plot.correlation$p.value),]
saveRDS(plot.correlation, "3.histoneLevel/H3K27me3/TSS_5kb/H3K27me3_expression_correlation.rds")

#### Gene ratio with H3K27me3 change
source("~/longzhilin/ProjectCode/AgingProject/plotDiffOverlap.R")
Expression.TSS_5kb.change <- lapply(tissue_cellTypes, function(tissue_cellType){
  cat("Analysis ", tissue_cellType, "...\n")
  TSS.data <- TSS_5kb.edgeR.list[[tissue_cellType]]
  RNA.data <- RNA.edgeR.list[[tissue_cellType]]

  RNA.data$type <- "ns"
  RNA.data$type[which(RNA.data$fdr < 0.05 & RNA.data$logFC > 0)] <- "up"
  RNA.data$type[which(RNA.data$fdr < 0.05 & RNA.data$logFC < 0)] <- "down"
  RNA.data$type <- factor(RNA.data$type, levels = c("up", "down", "ns"))
  total.RNA.up <- length(which(RNA.data$type == "up"))
  total.RNA.down <- length(which(RNA.data$type == "down"))

  TSS.data$type <- "ns"
  TSS.data$type[which(TSS.data[,sig.label] < 0.05 & TSS.data$logFC > log2FC.threshold)] <- "up"
  TSS.data$type[which(TSS.data[,sig.label] < 0.05 & TSS.data$logFC < -log2FC.threshold)] <- "down"
  TSS.data$type <- factor(TSS.data$type, levels = c("up", "down", "ns"))
  TSS.data <- TSS.data %>% left_join(gene_peak_unique, by = "gene")

  plot.data <- merge(RNA.data[,c("gene", "logFC", "fdr", "type")], TSS.data[,c("gene", "logFC", sig.label, "type", "class")], by = "gene", suffixes = c(".RNA", ".TSS"))
  # plot.data <- plot.data %>% filter(class == "OverlapPeak")
  Q1 <- sum(plot.data$type.RNA == "up" & plot.data$type.TSS == "down") # gene up, H3K27me3 down
  Q2 <- sum(plot.data$type.RNA == "down" & plot.data$type.TSS == "up") # gene down, H3K27me3 up
  Q3 <- total.RNA.up - Q1 # gene up with no/other H3K27me3 change
  Q4 <- total.RNA.down - Q2 # gene down with no/other H3K27me3 change
  up.ratio <- Q1 / total.RNA.up * 100
  down.ratio <- Q2 / total.RNA.down * 100

  # Fisher test
  data <- matrix(c(Q1, Q2, Q3, Q4), nrow = 2, ncol = 2, byrow = FALSE)
  res <- fisher.test(data)
  info <- data.frame(tissue_cellType = tissue_cellType, up.ratio = up.ratio, down.ratio = down.ratio, p.value = res$p.value)

  ## the percentage of RNA change with H3K27me3 change
  # upregulated genes with decreased H3K27me3
  RNA.up <- data.frame(tissue_cellType = tissue_cellType, type = c("total upregulated genes", "upregulated genes with decreased H3K27me3"), number = c(total.RNA.up, Q1))
  # downregulated genes with increased H3K27me3
  RNA.down <- data.frame(tissue_cellType = tissue_cellType, type = c("total downregulated genes", "downregulated genes with increased H3K27me3"), number = c(total.RNA.down, Q2))
  RNA.up$direction <- "Up"
  RNA.down$direction <- "Down"
  RNA.change <- rbind(RNA.up, RNA.down)

  ## the percentage of H3K27me3 change overlapping peak
  total.TSS.up <- length(which(TSS.data$type == "up"))
  total.TSS.down <- length(which(TSS.data$type == "down"))
  TSS_5kb.change <- data.frame(tissue_cellType = tissue_cellType, type = c("total upregulated H3K27me3 signls", "total downregulated H3K27me3 signls"), number = c(total.TSS.up, total.TSS.down), direction = c("Up", "Down"))
  return(list(info = info, RNA.change = RNA.change, TSS_5kb.change = TSS_5kb.change))
})
saveRDS(Expression.TSS_5kb.change, "3.histoneLevel/H3K27me3/TSS_5kb/Expression.TSS_5kb.change.rds")

## H3K27me3 & RNA change
Expression.TSS_5kb.correlation <- lapply(Expression.TSS_5kb.change, function(x){x$info})
Expression.TSS_5kb.correlation <- do.call(rbind, Expression.TSS_5kb.correlation)
Expression.TSS_5kb.correlation <- Expression.TSS_5kb.correlation %>% left_join(PT.info.H3K27me3.cellType, by = "tissue_cellType")
Expression.TSS_5kb.correlation <- Expression.TSS_5kb.correlation[order(Expression.TSS_5kb.correlation$up.ratio, decreasing = TRUE),]
Expression.TSS_5kb.correlation$label <- sapply(Expression.TSS_5kb.correlation$p.value, p_to_stars)
saveRDS(Expression.TSS_5kb.correlation, "3.histoneLevel/H3K27me3/TSS_5kb/Expression.TSS_5kb.correlation.rds")
writexl::write_xlsx(Expression.TSS_5kb.correlation, "3.histoneLevel/H3K27me3/TSS_5kb/Expression.TSS_5kb.correlation.xlsx")

## RNA change
RNA.change <- lapply(Expression.TSS_5kb.change, function(x){x$RNA.change})
RNA.change <- do.call(rbind, RNA.change)
RNA.change <- RNA.change %>% left_join(PT.info.RNA.cellType, by = "tissue_cellType")
RNA.change$cols <- as.character(plyr::mapvalues(RNA.change$category, from = names(color.category), to = color.category))
cols.info <- unique(RNA.change[,c("tissue_cellType", "cols")])
cols <- cols.info$cols
names(cols) <- cols.info$tissue_cellType
RNA.change.up <- RNA.change[which(RNA.change$direction == "Up"),]
plot.RNA.up <- plotDiffOverlap2(RNA.change.up[,c("tissue_cellType", "type", "number")], group = "tissue_cellType", cols = cols)
plot.RNA.up <- plot.RNA.up %>% left_join(PT.info.RNA.cellType[,c("tissue_cellType", "tissue_cellType_new", "tissue_new", "category")], by = "tissue_cellType")
RNA.change.down <- RNA.change[which(RNA.change$direction == "Down"),]
plot.RNA.down <- plotDiffOverlap2(RNA.change.down[,c("tissue_cellType", "type", "number")], group = "tissue_cellType", cols = cols)
plot.RNA.down <- plot.RNA.down %>% left_join(PT.info.RNA.cellType[,c("tissue_cellType", "tissue_cellType_new", "tissue_new", "category")], by = "tissue_cellType")
plot.RNA.down$number <- plot.RNA.down$number * -1
# select cell type with >= 250
cellType.up.DEG_filter <- plot.RNA.up[which(plot.RNA.up$type == "total upregulated genes"),]
cellType.up.DEG_filter <- cellType.up.DEG_filter$tissue_cellType[which(cellType.up.DEG_filter$number >= 250)]
cellType.down.DEG_filter <- plot.RNA.down[which(plot.RNA.down$type == "total downregulated genes"),]
cellType.down.DEG_filter <- cellType.down.DEG_filter$tissue_cellType[which(abs(cellType.down.DEG_filter$number) >= 250)]
cellType.DEG_filter <- unique(c(cellType.up.DEG_filter, cellType.down.DEG_filter)) # 84 cell types

## TSS_5kb change
TSS_5kb.change <- lapply(Expression.TSS_5kb.change, function(x){x$TSS_5kb.change})
TSS_5kb.change <- do.call(rbind, TSS_5kb.change)
TSS_5kb.change <- TSS_5kb.change %>% left_join(PT.info.H3K27me3.cellType, by = "tissue_cellType")
TSS_5kb.change$cols <- as.character(plyr::mapvalues(TSS_5kb.change$category, from = names(color.category), to = color.category))
cols.info <- unique(TSS_5kb.change[,c("tissue_cellType", "cols")])
cols <- cols.info$cols
names(cols) <- cols.info$tissue_cellType
TSS_5kb.change.up <- TSS_5kb.change[which(TSS_5kb.change$direction == "Up"),]
plot.TSS_5kb.up <- plotDiffOverlap2(TSS_5kb.change.up[,c("tissue_cellType", "type", "number")], group = "tissue_cellType", cols = cols)
plot.TSS_5kb.up <- plot.TSS_5kb.up %>% left_join(PT.info.H3K27me3.cellType[,c("tissue_cellType", "tissue_cellType_new", "tissue_new", "category")], by = "tissue_cellType")
TSS_5kb.change.down <- TSS_5kb.change[which(TSS_5kb.change$direction == "Down"),]
plot.TSS_5kb.down <- plotDiffOverlap2(TSS_5kb.change.down[,c("tissue_cellType", "type", "number")], group = "tissue_cellType", cols = cols)
plot.TSS_5kb.down <- plot.TSS_5kb.down %>% left_join(PT.info.H3K27me3.cellType[,c("tissue_cellType", "tissue_cellType_new", "tissue_new", "category")], by = "tissue_cellType")
plot.TSS_5kb.down$number <- plot.TSS_5kb.down$number * -1

plot.data <- Expression.TSS_5kb.correlation
plot.data <- plot.data[order(plot.data$up.ratio, decreasing = T),]
plot.data <- plot.data %>% dplyr::rename("H3K27me3.reads" = read) %>% left_join(PT.info.RNA.cellType[, c("tissue_cellType", "read")] %>% dplyr::rename("RNA.reads" = read), by = "tissue_cellType")

# filter low read samples
selected.cellTypes <- intersect(cellType.DEG_filter, selected.tissue_cellType$tissue_cellType)
plot.data <- subset(plot.data, tissue_cellType %in% selected.cellTypes)

plot.data$tissue_cellType_new <- factor(plot.data$tissue_cellType_new, levels = plot.data$tissue_cellType_new)
plot.RNA.up <- plot.RNA.up[which(plot.RNA.up$tissue_cellType_new %in% plot.data$tissue_cellType_new),]
plot.RNA.up$tissue_cellType_new <- factor(plot.RNA.up$tissue_cellType_new, levels = plot.data$tissue_cellType_new)
plot.RNA.up <- arrange(plot.RNA.up, tissue_cellType_new)
plot.RNA.down <- plot.RNA.down[which(plot.RNA.down$tissue_cellType_new %in% plot.data$tissue_cellType_new),]
plot.RNA.down$tissue_cellType_new <- factor(plot.RNA.down$tissue_cellType_new, levels = plot.data$tissue_cellType_new)
plot.RNA.down <- arrange(plot.RNA.down, tissue_cellType_new)
plot.TSS_5kb.up <- plot.TSS_5kb.up[which(plot.TSS_5kb.up$tissue_cellType_new %in% plot.data$tissue_cellType_new),]
plot.TSS_5kb.up$tissue_cellType_new <- factor(plot.TSS_5kb.up$tissue_cellType_new, levels = plot.data$tissue_cellType_new)
plot.TSS_5kb.up <- arrange(plot.TSS_5kb.up, tissue_cellType_new)
plot.TSS_5kb.down <- plot.TSS_5kb.down[which(plot.TSS_5kb.down$tissue_cellType_new %in% plot.data$tissue_cellType_new),]
plot.TSS_5kb.down$tissue_cellType_new <- factor(plot.TSS_5kb.down$tissue_cellType_new, levels = plot.data$tissue_cellType_new)
plot.TSS_5kb.down <- arrange(plot.TSS_5kb.down, tissue_cellType_new)

symmetric_limit <- max(abs(na.omit(plot.data$up.ratio)), abs(na.omit(plot.data$down.ratio)))
## annother style
plot.data <- plot.data %>% mutate(down.abs = abs(down.ratio), neglog10_p = -log10(pmax(p.value, .Machine$double.xmin)))
label.data <- plot.data %>% arrange(desc(up.ratio)) %>% slice_head(n = 3)
max_value <- max(c(plot.data$down.abs, plot.data$up.ratio), na.rm = TRUE)
axis_breaks <- pretty(c(0, max_value), n = 5)
axis_max <- max(axis_breaks)
p1 <- ggplot(plot.data, aes(x = down.abs, y = up.ratio, size = neglog10_p)) +
  geom_point(alpha = 0.85, color = "#B93F2B") +
  geom_text_repel(data = label.data, aes(x = down.abs, y = up.ratio, label = tissue_cellType_new), inherit.aes = FALSE, color = "black", size = 3.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_x_continuous(limits = c(0, axis_max), breaks = axis_breaks, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, axis_max), breaks = axis_breaks, expand = c(0, 0)) +
  scale_size_continuous(range = c(1, 5), name = expression(-log[10](p.value))) +
  coord_fixed(ratio = 1) +
  labs(
    x = "Downregulated genes with increased H3K27me3 (%)",
    y = "Upregulated genes with decreased H3K27me3 (%)"
  ) +
  theme_classic() + theme(aspect.ratio = 1)

plot.data.rank.wide <- plot.data %>% arrange(desc(up.ratio)) %>% mutate(rank = row_number(), down.abs = abs(down.ratio))
plot.data.rank <- plot.data.rank.wide %>% pivot_longer(cols = c(up.ratio, down.abs), names_to = "direction", values_to = "ratio")
label.data <- plot.data.rank.wide %>% slice_head(n = 3) %>% mutate(label_y = pmax(up.ratio, down.abs))
p2 <- ggplot(plot.data.rank, aes(x = rank, y = ratio, color = direction)) +
  geom_segment(data = plot.data.rank.wide, aes(x = rank, xend = rank, y = up.ratio, yend = down.abs), inherit.aes = FALSE, color = "grey60", linetype = "dashed",linewidth = 0.5) +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_text_repel(data = label.data, aes(x = rank, y = label_y, label = tissue_cellType_new), inherit.aes = FALSE, color = "black", size = 3.5, nudge_y = 0.5) +
  scale_color_manual(values = c(up.ratio = "#E64B35", down.abs = "#00468B"),
    labels = c(
      up.ratio = "Upregulated genes with decreased H3K27me3",
      down.abs = "Downregulated genes with increased H3K27me3"
    )
  ) +
  scale_x_continuous(breaks = pretty(plot.data.rank$rank)) +
  labs(x = "Ranked tissue-cell type", y = "Ratio (%)", color = NULL) +
  theme_classic()
pdf("3.histoneLevel/H3K27me3/TSS_5kb/Expression.TSS_5kb.relationship.filter2.pdf", height = 5, width = 12)
print(p1)
print(p2)
dev.off()

#### output the correlation scatter plot data
RNA.change.Up <- subset(RNA.change, type == "total upregulated genes")
RNA.change.Down <- subset(RNA.change, type == "total downregulated genes")
TSS_5kb.change.Up <- subset(TSS_5kb.change, type == "total upregulated H3K27me3 signls")
TSS_5kb.change.Down <- subset(TSS_5kb.change, type == "total downregulated H3K27me3 signls")
plot.correlation <- left_join(plot.correlation, RNA.change.Up[,c("tissue_cellType", "read", "number")] %>% dplyr::rename("RNA.read" = read, "RNA.up" = number), by = c("tissue_cellType"))
plot.correlation <- left_join(plot.correlation, RNA.change.Down[,c("tissue_cellType", "number")] %>% dplyr::rename("RNA.down" = number), by = c("tissue_cellType"))
plot.correlation <- left_join(plot.correlation, TSS_5kb.change.Up[,c("tissue_cellType", "read", "number")] %>% dplyr::rename("H3K27me3.read" = read, "H3K27me3.up" = number), by = c("tissue_cellType"))
plot.correlation <- left_join(plot.correlation, TSS_5kb.change.Down[,c("tissue_cellType", "number")] %>% dplyr::rename("H3K27me3.down" = number), by = c("tissue_cellType"))
writexl::write_xlsx(plot.correlation, "3.histoneLevel/H3K27me3/TSS_5kb/H3K27me3_expression_correlation.xlsx")
saveRDS(plot.correlation, "3.histoneLevel/H3K27me3/TSS_5kb/H3K27me3_expression_correlation.rds")

#### ---- Frequent differential signals in different tissues and cell types ---- ####
source("~/Functions/clusterProfiler.enricher.R")
source("~/Functions/pathway.dotplot.R")
log2FC.threshold <- log2(1.2)
sig.label <- "PValue"

#### H3K27me3
TSS_5kb.edgeR.list <- readRDS("3.histoneLevel/H3K27me3/TSS_5kb/TSS_5kb.edgeR.list.rds")
## Up-regulated bins: extract significant bins with cell type information
TSS_5kb.up.df <- lapply(names(TSS_5kb.edgeR.list), function(ct){
  x <- TSS_5kb.edgeR.list[[ct]]
  idx <- which(x[, sig.label] < 0.05 & x$logFC > log2FC.threshold)
  if(length(idx) > 0){
    return(data.frame(gene = x$gene[idx], cellType = ct, stringsAsFactors = FALSE))
  }else{
    return(NULL)
  }
})
TSS_5kb.up.df <- do.call(rbind, TSS_5kb.up.df)
TSS_5kb.up.Freq <- TSS_5kb.up.df %>% group_by(gene) %>% summarise(n_sig_up_cellType = n(), sig_up_cellType = paste(cellType, collapse = ";"), .groups = "drop") %>% arrange(desc(n_sig_up_cellType))

## Down-regulated bins: extract significant bins with cell type information
TSS_5kb.down.df <- lapply(names(TSS_5kb.edgeR.list), function(ct){
  x <- TSS_5kb.edgeR.list[[ct]]
  idx <- which(x[, sig.label] < 0.05 & x$logFC < (log2FC.threshold * -1))
  if(length(idx) > 0){
    return(data.frame(gene = x$gene[idx], cellType = ct, stringsAsFactors = FALSE))
  }else{
    return(NULL)
  }
})
TSS_5kb.down.df <- do.call(rbind, TSS_5kb.down.df)
TSS_5kb.down.Freq <- TSS_5kb.down.df %>% group_by(gene) %>% summarise(n_sig_down_cellType = n(), sig_down_cellType = paste(cellType, collapse = ";"), .groups = "drop") %>% arrange(desc(n_sig_down_cellType))

TSS_5kb.Frequency <- merge(TSS_5kb.up.Freq, TSS_5kb.down.Freq, all = T, by = "gene")
TSS_5kb.Frequency$sig_up_cellType[is.na(TSS_5kb.Frequency$sig_up_cellType)] <- ""
TSS_5kb.Frequency$sig_down_cellType[is.na(TSS_5kb.Frequency$sig_down_cellType)] <- ""
TSS_5kb.Frequency[is.na(TSS_5kb.Frequency)] <- 0
TSS_5kb.Frequency$n_sig_cellType <- rowSums(TSS_5kb.Frequency[,c("n_sig_up_cellType", "n_sig_down_cellType")])
TSS_5kb.Frequency <- TSS_5kb.Frequency[,c("gene", "n_sig_cellType", "n_sig_up_cellType", "sig_up_cellType", "n_sig_down_cellType", "sig_down_cellType")]
TSS_5kb.Frequency <- TSS_5kb.Frequency[order(TSS_5kb.Frequency$n_sig_up_cellType, decreasing = T),]

symmetric_limit <- max(max(TSS_5kb.up.Freq$n_sig_up_cellType), max(TSS_5kb.down.Freq$n_sig_down_cellType))
p1 <- gghistogram(TSS_5kb.up.Freq, x = "n_sig_up_cellType", fill = "lightgray", bins = 50, xlab = "The number of cell type", ylab = "The number of differential bin") + ggtitle("Differential-up TSS_5kb(Old/Young)") + xlim(0, symmetric_limit) + theme(plot.title = element_text(hjust = 0.5))
p2 <- gghistogram(TSS_5kb.down.Freq, x = "n_sig_down_cellType", fill = "lightgray", bins = 50, xlab = "The number of cell type", ylab = "The number of differential bin") + ggtitle("Differential-down TSS_5kb(Old/Young)") + xlim(0, symmetric_limit) + theme(plot.title = element_text(hjust = 0.5))
pdf("3.histoneLevel/H3K27me3/TSS_5kb/TSS_5kb.Frequency.pdf")
ggarrange(p1, NULL, p2, NULL, ncol = 2, nrow = 2, align = "hv")
dev.off()

## pathway enrichment for top 100
plot.merge.up <- unique(as.data.frame(TSS_5kb.up.Freq[,c("gene", "n_sig_up_cellType")]))
select.up.genes <- as.character(plot.merge.up[order(plot.merge.up$n_sig_up_cellType, decreasing = T),][1:200, "gene"])
plot.merge.down <- unique(as.data.frame(TSS_5kb.down.Freq[,c("gene", "n_sig_down_cellType")]))
select.down.genes <- as.character(plot.merge.down[order(plot.merge.down$n_sig_down_cellType, decreasing = T),][1:200, "gene"])
H3K27me3.list <- list(Increased = select.up.genes, Decreased = select.down.genes)
enrichedPathways <- lapply(H3K27me3.list, function(x){
  res <- clusterProfiler.enricher(gene = x, geneType = "SYMBOL", db.type = "GO", GO.ont = "BP", qvalueCutoff = 0.05, pvalueCutoff = 0.05, species = "Mus musculus")
  return(res)
})
p <- pathway.dotplot(enrichedResult.lists = enrichedPathways, title = "", items = c("Description", "GeneRatio", "p.adjust"), threshold = 0.05, top.n = 15, GeneRatio = F, max.dot.size = 4, sub.words = 80, y.text.size = 11)
print(p)

## add gene expression (CPM: Old - Young)
RNA.cpm.mats <- readRDS("1.dataPreprocessing/RNA.cpm.mats.rds")
RNA.cpm.mats[is.na(RNA.cpm.mats)] <- 0
condition <- data.frame(label = colnames(RNA.cpm.mats))
condition$age <- gsub("Old--.*", "Old", condition$label)
condition$age <- gsub("Young--.*", "Young", condition$age)
condition$tissue_cellType <- gsub("Old--|Young--", "", condition$label)
condition$tissue_cellType <- gsub("--", ":", condition$tissue_cellType)
condition$tissue_cellType <- gsub(" |-", "_", condition$tissue_cellType)
cpm.Old <- RNA.cpm.mats[,which(condition$age == "Old")]
cpm.Young <- RNA.cpm.mats[,which(condition$age == "Young")]
cpm.diff <- cpm.Old - cpm.Young
cpm.mean <- rowMeans(cpm.diff)
pos_count <- rowSums(cpm.diff > 0)
neg_count <- rowSums(cpm.diff < 0)
cpm.diff.info <- data.frame(gene = rownames(cpm.diff), RNA.cpm.mean = cpm.mean, RNA.pos_count = pos_count, RNA.neg_count = neg_count)
TSS_5kb.Frequency <- left_join(TSS_5kb.Frequency, cpm.diff.info, by = "gene")
saveRDS(TSS_5kb.Frequency, "3.histoneLevel/H3K27me3/TSS_5kb/TSS_5kb.Frequency.rds")
writexl::write_xlsx(TSS_5kb.Frequency, "3.histoneLevel/H3K27me3/TSS_5kb/TSS_5kb.Frequency.xlsx")

#### H3K27me3 & gene expression (edgeR)
## Up/Down-regulated genes
RNA.edgeR.list <- readRDS("2.rnaLevel/differentialGeneAnalysis/RNA.edgeR.list.rds")
RNA.edgeR.list <- unlist(RNA.edgeR.list, recursive = F)
names(RNA.edgeR.list) <- gsub("\\.", ":", names(RNA.edgeR.list))
DEG.up.df <- lapply(names(RNA.edgeR.list), function(ct){
  x <- RNA.edgeR.list[[ct]]
  idx <- which(x$fdr < 0.05 & x$logFC > 0)
  if(length(idx) > 0){
    return(data.frame(gene = x$gene[idx], cellType = ct, stringsAsFactors = FALSE))
  }else{
    return(NULL)
  }
})
DEG.up.df <- do.call(rbind, DEG.up.df)

DEG.down.df <- lapply(names(RNA.edgeR.list), function(ct){
  x <- RNA.edgeR.list[[ct]]
  idx <- which(x$fdr < 0.05 & x$logFC < 0)
  if(length(idx) > 0){
    return(data.frame(gene = x$gene[idx], cellType = ct, stringsAsFactors = FALSE))
  }else{
    return(NULL)
  }
})
DEG.down.df <- do.call(rbind, DEG.down.df)

# gene up with decreased H3K27me3
DEG_up.TSS_5kb_down <- merge(DEG.up.df, TSS_5kb.down.df, by = c("gene", "cellType"))
DEG_up.TSS_5kb_down.Freq <- DEG_up.TSS_5kb_down %>% group_by(gene) %>% summarise(n_sig_DEG.up.H3K27me3.down_cellType = n(), sig_DEG.up.H3K27me3.down_cellType = paste(cellType, collapse = ";"), .groups = "drop") %>% arrange(desc(n_sig_DEG.up.H3K27me3.down_cellType))

# gene down with increased H3K27me3
DEG_down.TSS_5kb_up <- merge(DEG.down.df, TSS_5kb.up.df, by = c("gene", "cellType"))
DEG_down.TSS_5kb_up.Freq <- DEG_down.TSS_5kb_up %>% group_by(gene) %>% summarise(n_sig_DEG.down.H3K27me3.up_cellType = n(), sig_DEG.down.H3K27me3.up_cellType = paste(cellType, collapse = ";"), .groups = "drop") %>% arrange(desc(n_sig_DEG.down.H3K27me3.up_cellType))

# plot distribution of frequency
symmetric_limit <- max(max(DEG_up.TSS_5kb_down.Freq$n_sig_DEG.up.H3K27me3.down_cellType), max(DEG_down.TSS_5kb_up.Freq$n_sig_DEG.down.H3K27me3.up_cellType))
p1 <- gghistogram(DEG_up.TSS_5kb_down.Freq[,c("gene", "n_sig_DEG.up.H3K27me3.down_cellType")], x = "n_sig_DEG.up.H3K27me3.down_cellType", fill = "lightgray", bins = 50, xlab = "The number of cell type", ylab = "The number of gene") + ggtitle("Gene up with decreased H3K27me3") + xlim(0, symmetric_limit) + theme(plot.title = element_text(hjust = 0.5))
p2 <- gghistogram(DEG_down.TSS_5kb_up.Freq[,c("gene", "n_sig_DEG.down.H3K27me3.up_cellType")], x = "n_sig_DEG.down.H3K27me3.up_cellType", fill = "lightgray", bins = 50, xlab = "The number of cell type", ylab = "The number of gene") + ggtitle("Gene down with increased H3K27me3") + xlim(0, symmetric_limit) + theme(plot.title = element_text(hjust = 0.5))
pdf("3.histoneLevel/H3K27me3/TSS_5kb/DEG.H3K27me3.Frequency.pdf")
ggarrange(p1, NULL, p2, NULL, ncol = 2, nrow = 2, align = "hv")
dev.off()

DEG.H3K27me3.Frequency <- merge(DEG_up.TSS_5kb_down.Freq, DEG_down.TSS_5kb_up.Freq, all = T, by = "gene")
DEG.H3K27me3.Frequency$sig_DEG.up.H3K27me3.down_cellType[is.na(DEG.H3K27me3.Frequency$sig_DEG.up.H3K27me3.down_cellType)] <- ""
DEG.H3K27me3.Frequency$sig_DEG.down.H3K27me3.up_cellType[is.na(DEG.H3K27me3.Frequency$sig_DEG.down.H3K27me3.up_cellType)] <- ""
DEG.H3K27me3.Frequency[is.na(DEG.H3K27me3.Frequency)] <- 0
DEG.H3K27me3.Frequency$n_sig_cellType <- rowSums(DEG.H3K27me3.Frequency[,c("n_sig_DEG.up.H3K27me3.down_cellType", "n_sig_DEG.down.H3K27me3.up_cellType")])
DEG.H3K27me3.Frequency <- DEG.H3K27me3.Frequency[,c("gene", "n_sig_cellType", "n_sig_DEG.up.H3K27me3.down_cellType", "sig_DEG.up.H3K27me3.down_cellType", "n_sig_DEG.down.H3K27me3.up_cellType", "sig_DEG.down.H3K27me3.up_cellType")]
DEG.H3K27me3.Frequency <- DEG.H3K27me3.Frequency[order(DEG.H3K27me3.Frequency$n_sig_DEG.up.H3K27me3.down_cellType, decreasing = T),]
saveRDS(DEG.H3K27me3.Frequency, "3.histoneLevel/H3K27me3/TSS_5kb/DEG.H3K27me3.Frequency.rds")
writexl::write_xlsx(DEG.H3K27me3.Frequency, "3.histoneLevel/H3K27me3/TSS_5kb/DEG.H3K27me3.Frequency.xlsx")

## plot high-up-frequency
plot.data <- DEG.H3K27me3.Frequency[which(DEG.H3K27me3.Frequency$n_sig_DEG.up.H3K27me3.down_cellType > 6),]
mean.logFC <- apply(plot.data, 1, function(x){
  cellTypes <- unlist(strsplit(x[4], ";"))
  cellTypes <- cellTypes[which(cellTypes != "")]
  # mean RNA-logFC
  mean.RNA.logFC <- lapply(cellTypes, function(ct){
    RNA.edgeR.list[[ct]]$logFC[which(RNA.edgeR.list[[ct]]$gene == x[1])]
  }) %>% unlist()
  mean.RNA.logFC <- mean(mean.RNA.logFC)

  # mean H3K27me3-logFC
  mean.H3K27me3.logFC <- lapply(cellTypes, function(ct){
    TSS_5kb.edgeR.list[[ct]]$logFC[which(TSS_5kb.edgeR.list[[ct]]$gene == x[1])]
  }) %>% unlist()
  mean.H3K27me3.logFC <- mean(mean.H3K27me3.logFC)
  mean.logFC <- data.frame(gene = x[1], mean.RNA.logFC = mean.RNA.logFC, mean.H3K27me3.logFC = mean.H3K27me3.logFC, Freq = as.numeric(x[3]))
  return(mean.logFC)
})
mean.logFC <- do.call(rbind, mean.logFC)

## scatter plot: mean.RNA.logFC vs mean.H3K27me3.logFC, point size by Freq, label top 10 genes
library(ggrepel)
genes <- c("Zc3hav1", "Tgfbr3", "H2-K1", "H2-D1", "H2-Ab1", "Nlrc5", "Runx1", "Ezh2", "Bach2", "Nrxn3", "Ptprz1", "Evl")
mean.logFC <- mean.logFC[order(mean.logFC$Freq, decreasing = TRUE),]
range(mean.logFC$mean.H3K27me3.logFC)
top10_genes <- unique(c(head(mean.logFC$gene, 5), genes))
points_to_label <- subset(mean.logFC, gene %in% top10_genes)

p1 <- ggplot(mean.logFC, aes(x = Freq, y = mean.RNA.logFC)) + geom_point(aes(color = mean.H3K27me3.logFC), color = "#E64B35", alpha = 0.8, size = 3) +
  # scale_color_viridis_c(option = "magma", direction = 1, name = "H3K27me3:\nMean logFC") +
  geom_text_repel(data = points_to_label, aes(x = Freq, y = mean.RNA.logFC, label = gene), size = 4, box.padding = 0.5, point.padding = 0.3, force = 2, max.overlaps = Inf, min.segment.length = 0) +
  labs(x = "# of cell types with increased RNA & decreased H3K27me3", y = "RNA: mean logFC", title = "") +
  theme_classic() + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
p2 <- ggplot(mean.logFC, aes(x = Freq, y = mean.H3K27me3.logFC)) + geom_point(aes(color = mean.RNA.logFC), color = "#0874C0", alpha = 0.8, size = 3) +
  # scale_color_viridis_c(option = "plasma", direction = 1, name = "RNA:\nMean logFC") +
  geom_text_repel(data = points_to_label, aes(x = Freq, y = mean.H3K27me3.logFC, label = gene), size = 4, box.padding = 0.5, point.padding = 0.3, force = 2, max.overlaps = Inf, min.segment.length = 0) +
  labs(x = "# of cell types with increased RNA & decreased H3K27me3", y = "H3K27me3: mean logFC", title = "") +
  theme_classic() + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))

## plot high-down-frequency
plot.data <- DEG.H3K27me3.Frequency[which(DEG.H3K27me3.Frequency$n_sig_DEG.down.H3K27me3.up_cellType > 0),]
mean.logFC <- apply(plot.data, 1, function(x){
  cellTypes <- unlist(strsplit(x[6], ";"))
  cellTypes <- cellTypes[which(cellTypes != "")]
  # mean RNA-logFC
  mean.RNA.logFC <- lapply(cellTypes, function(ct){
    RNA.edgeR.list[[ct]]$logFC[which(RNA.edgeR.list[[ct]]$gene == x[1])]
  }) %>% unlist()
  mean.RNA.logFC <- mean(mean.RNA.logFC)

  # mean H3K27me3-logFC
  mean.H3K27me3.logFC <- lapply(cellTypes, function(ct){
    TSS_5kb.edgeR.list[[ct]]$logFC[which(TSS_5kb.edgeR.list[[ct]]$gene == x[1])]
  }) %>% unlist()
  mean.H3K27me3.logFC <- mean(mean.H3K27me3.logFC)
  mean.logFC <- data.frame(gene = x[1], mean.RNA.logFC = mean.RNA.logFC, mean.H3K27me3.logFC = mean.H3K27me3.logFC, Freq = as.numeric(x[5]))
  return(mean.logFC)
})
mean.logFC <- do.call(rbind, mean.logFC)

## scatter plot: mean.RNA.logFC vs mean.H3K27me3.logFC, point size by Freq, label top 10 genes
library(ggrepel)
mean.logFC <- mean.logFC[order(mean.logFC$Freq, decreasing = TRUE),]
top10_genes <- unique(c(head(mean.logFC$gene, 5), genes))
points_to_label <- subset(mean.logFC, gene %in% top10_genes)

p11 <- ggplot(mean.logFC, aes(x = Freq, y = mean.RNA.logFC)) + geom_point(aes(color = mean.H3K27me3.logFC), alpha = 0.8, color = "#E64B35", size = 3) +
  # scale_color_viridis_c(option = "magma", direction = 1, name = "H3K27me3:\nMean logFC") +
  geom_text_repel(data = points_to_label, aes(x = Freq, y = mean.RNA.logFC, label = gene), size = 4, box.padding = 0.5, point.padding = 0.3, force = 2, max.overlaps = Inf, min.segment.length = 0) +
  labs(x = "# of cell types with decreased RNA & increased H3K27me3", y = "RNA: mean logFC", title = "") +
  theme_classic() + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
p22 <- ggplot(mean.logFC, aes(x = Freq, y = mean.H3K27me3.logFC)) + geom_point(aes(color = mean.RNA.logFC), alpha = 0.8, color = "#0874C0", size = 3) +
  # scale_color_viridis_c(option = "plasma", direction = 1, name = "RNA:\nMean logFC") +
  geom_text_repel(data = points_to_label, aes(x = Freq, y = mean.H3K27me3.logFC, label = gene), size = 4, box.padding = 0.5, point.padding = 0.3, force = 2, max.overlaps = Inf, min.segment.length = 0) +
  labs(x = "# of cell types with decreased RNA & increased H3K27me3", y = "H3K27me3: mean logFC", title = "") +
  theme_classic() + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))

pdf("3.histoneLevel/H3K27me3/TSS_5kb/Frequency.mean.logFC.scatter.pdf", width = 10, height = 10)
print(ggarrange(p1, p2, p11, p22, ncol = 2, nrow = 2, align = "hv"))
dev.off()

## pathway enrichment for freq > 6
select.up.genes <- as.character(DEG.H3K27me3.Frequency$gene[which(DEG.H3K27me3.Frequency$n_sig_DEG.up.H3K27me3.down_cellType > 6)])
select.down.genes <- as.character(DEG.H3K27me3.Frequency$gene[which(DEG.H3K27me3.Frequency$n_sig_DEG.down.H3K27me3.up_cellType > 0)])
RNA.list <- list(DEG.up.H3K27me3.down = select.up.genes, DEG.down.H3K27me3.up = select.down.genes)
enrichedPathways <- lapply(RNA.list, function(x){
  res <- clusterProfiler.enricher(gene = x, geneType = "SYMBOL", db.type = "GO", GO.ont = "BP", qvalueCutoff = 0.05, pvalueCutoff = 0.05, species = "Mus musculus", return.clusterProfiler = T)
  return(res)
})
saveRDS(enrichedPathways, "3.histoneLevel/H3K27me3/TSS_5kb/DEG.H3K27me3.FrequencyOver6.enrichedPathways.rds")

## calculate the mean-logFC across category
genes <- unlist(strsplit(enrichedPathways$DEG.up.H3K27me3.down$geneID[1], "/"))
plot.data <- DEG.H3K27me3.Frequency[which(DEG.H3K27me3.Frequency$gene %in% genes),]
pdf("3.histoneLevel/H3K27me3/TSS_5kb/Pathway.meanLogFC.pdf", width = 5, height = 5)
mean.logFC <- apply(plot.data, 1, function(x){
  cellTypes <- unlist(strsplit(x[4], ";"))
  cellTypes <- cellTypes[which(cellTypes != "")]
  # mean RNA-logFC
  RNA.logFC <- lapply(cellTypes, function(ct){
    RNA.edgeR.list[[ct]]$logFC[which(RNA.edgeR.list[[ct]]$gene == x[1])]
  }) %>% unlist()
  RNA.logFC <- data.frame(tissue_cellType = cellTypes, RNA.logFC = RNA.logFC)
  RNA.logFC <- PT.info.H3K27me3.cellType[,c("tissue_cellType", "category")] %>% left_join(RNA.logFC, by = "tissue_cellType")
  RNA.logFC$RNA.logFC[is.na(RNA.logFC$RNA.logFC)] <- 0
  mean.RNA.logFC <- tapply(RNA.logFC$RNA.logFC, RNA.logFC$category, mean)

  # mean H3K27me3-logFC
  H3K27me3.logFC <- lapply(cellTypes, function(ct){
    TSS_5kb.edgeR.list[[ct]]$logFC[which(TSS_5kb.edgeR.list[[ct]]$gene == x[1])]
  }) %>% unlist()
  H3K27me3.logFC <- data.frame(tissue_cellType = cellTypes, H3K27me3.logFC = H3K27me3.logFC)
  H3K27me3.logFC <- PT.info.H3K27me3.cellType[,c("tissue_cellType", "category")] %>% left_join(H3K27me3.logFC, by = "tissue_cellType")
  H3K27me3.logFC$H3K27me3.logFC[is.na(H3K27me3.logFC$H3K27me3.logFC)] <- 0
  mean.H3K27me3.logFC <- tapply(H3K27me3.logFC$H3K27me3.logFC, H3K27me3.logFC$category, mean)

  RNA.H3K27me3.data <- RNA.logFC %>% left_join(H3K27me3.logFC[,c("tissue_cellType", "H3K27me3.logFC")], by = "tissue_cellType")
  Freq <- table(RNA.H3K27me3.data$category[RNA.H3K27me3.data$RNA.logFC > 0])
  mean.logFC <- data.frame(gene = x[1], Freq = as.numeric(Freq[names(mean.RNA.logFC)]), category = names(mean.RNA.logFC), mean.RNA.logFC = mean.RNA.logFC, mean.H3K27me3.logFC = mean.H3K27me3.logFC[names(mean.RNA.logFC)], row.names = NULL)
  mean.logFC$category <- factor(mean.logFC$category, levels = names(color.category))
  p <- ggscatter(mean.logFC, y = "mean.RNA.logFC", x = "mean.H3K27me3.logFC", color = "category", palette = color.category, size = "Freq", ylab = "RNA:mean logFC", xlab = "H3K27me3:mean logFC") + ggtitle(x[1]) + scale_size(range = c(1, 6), name = "n.sig") + scale_x_reverse() + theme(plot.title = element_text(hjust = 0.5), legend.position = "right", aspect.ratio = 1)
  print(p)
  return(mean.logFC)
})
dev.off()
mean.logFC <- do.call(rbind, mean.logFC)
mean.logFC$category <- factor(mean.logFC$category, levels = names(color.category))
