#' @description differential gene expression
#' @Author: Zhilin Long
#' @date 2026-6-18

library(Seurat)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggrepel)
library(gridExtra)
set.seed(101)
source("~/Functions/colorAssign.R")
source("~/Functions/utilities.R")
setwd("~/AgingMouse/PairedTag")

RNA.cpm.mats <- readRDS("1.dataPreprocessing/RNA.cpm.mats.rds")
PT.info <- readRDS("1.dataPreprocessing/PT.info.rds")
PT.info.RNA.cellType <- PT.info %>% select(tissue_cellType, tissue_cellType_new, cellType, tissue, tissue_new, category, RNA.depth) %>% group_by(tissue_cellType) %>% summarise(tissue_cellType_new = unique(tissue_cellType_new), cellType = unique(cellType), read = sum(RNA.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), category = unique(category), .groups = "drop")
RNA_age_sample_replication <- readRDS("1.dataPreprocessing/RNA_age_sample_replication.rds")

tissue.seurat.final <- readRDS("1.dataPreprocessing/tissue.seurat.final.rds")
#### ---- edgeR analysis ---- ####
source("/storage/zhangyanxiaoLab/longzhilin/ProjectCode/AgingProject/edgeR_gene.R")
library(edgeR)
dir.create("2.rnaLevel/differentialGeneAnalysis", showWarnings = FALSE, recursive = TRUE)
RNA.edgeR.list <- lapply(tissue.seurat.final, function(seuratObj){
  cat(seuratObj$tissue[1], "\n")
  seuratObj$batch <- gsub("-.*", "", seuratObj$assay)
  seuratObj$batch <- factor(seuratObj$batch, levels = unique(seuratObj$batch))
  seuratObj$group <- paste0(seuratObj$cell_type, "_", seuratObj$age_group, "_", seuratObj$sample, "_", seuratObj$batch) # add histone as replication
  counts <- AggregateExpression(seuratObj, assays = "RNA", group.by = "group")$RNA

  condition <- data.frame(label = colnames(counts))
  condition$tissue <- seuratObj$tissue[1]
  condition$cellType <- gsub("-mouse.*", "", condition$label)
  condition$age_group <- gsub(".*-", "", condition$cellType)
  condition$cellType <- gsub("-Old|-Young", "", condition$cellType)
  condition$cellType <- gsub(" |-", "_", condition$cellType)
  condition$sample <- gsub(".*-mouse", "mouse", condition$label)
  condition$batch <- gsub(".*-", "", condition$sample)
  condition$sample <- gsub("-.*", "", condition$sample)
  counts <- counts[, condition$label]
  pdf(paste0("2.rnaLevel/differentialGeneAnalysis/", seuratObj$tissue[1], ".edgeR.pdf"), width = 6, height = 6)
  cellType.DEGs <- lapply(unique(condition$cellType), function(cellType){
    cat("    ", cellType, "\n")
    counts.subset <- counts[,which(condition$cellType == cellType)]
    cpms.subset <- edgeR::cpm(counts.subset)
    condition.subset <- condition[which(condition$cellType == cellType),]
    group <- factor(condition.subset$age_group, levels = names(color.age_group))
    batch <- factor(condition.subset$batch)
    myContrasts <- "Old-Young"
    out <- edgeR_gene_simple(counts = counts.subset, group = group, myContrasts = myContrasts, batch = batch, cpm.threshold = 1, cpm.filter = TRUE)
    p <- edgeR_plot(edgeR_result = out, PValue = FALSE, threshold = 0.05, top.n = 5, size = 3.5, title = paste0(seuratObj$tissue[1], ": ", cellType))
    print(p)
    return(out)
  })
  dev.off()
  names(cellType.DEGs) <- unique(condition$cellType)
  return(cellType.DEGs)
})
saveRDS(RNA.edgeR.list, "2.rnaLevel/differentialGeneAnalysis/RNA.edgeR.list.rds") # output RNA.condition.csv

#### ---- Frequent differential genes in different tissues and cell types ---- ####
source("~/Functions/utilities.R")
log2FC.threshold <- 0
sig.label <- "fdr"
RNA.edgeR.list <- readRDS("2.rnaLevel/differentialGeneAnalysis/RNA.edgeR.list.rds")
RNA.edgeR.list <- unlist(RNA.edgeR.list, recursive = F)
names(RNA.edgeR.list) <- gsub("\\.", ":", names(RNA.edgeR.list))

#### The celltype frequency of differential genes
## Up-regulated genes: extract significant genes with cell type information
DEG.up.df <- lapply(names(RNA.edgeR.list), function(ct){
  x <- RNA.edgeR.list[[ct]]
  idx <- which(x[, sig.label] < 0.05 & x$logFC > log2FC.threshold)
  if(length(idx) > 0){
    return(data.frame(gene = x$gene[idx], cellType = ct, stringsAsFactors = FALSE))
  }else{
    return(NULL)
  }
})
DEG.up.df <- do.call(rbind, DEG.up.df)
DEG.up.Freq <- DEG.up.df %>% group_by(gene) %>% summarise(n_sig_up_cellType = n(), sig_up_cellType = paste(cellType, collapse = ";"), .groups = "drop") %>% arrange(desc(n_sig_up_cellType))

## Down-regulated genes: extract significant genes with cell type information
DEG.down.df <- lapply(names(RNA.edgeR.list), function(ct){
  x <- RNA.edgeR.list[[ct]]
  idx <- which(x[, sig.label] < 0.05 & x$logFC < (log2FC.threshold * -1))
  if(length(idx) > 0){
    return(data.frame(gene = x$gene[idx], cellType = ct, stringsAsFactors = FALSE))
  }else{
    return(NULL)
  }
})
DEG.down.df <- do.call(rbind, DEG.down.df)
DEG.down.Freq <- DEG.down.df %>% group_by(gene) %>% summarise(n_sig_down_cellType = n(), sig_down_cellType = paste(cellType, collapse = ";"), .groups = "drop") %>% arrange(desc(n_sig_down_cellType))

symmetric_limit <- max(max(DEG.up.Freq$n_sig_up_cellType), max(DEG.down.Freq$n_sig_down_cellType))
p1 <- gghistogram(DEG.up.Freq[,c("gene", "n_sig_up_cellType")], x = "n_sig_up_cellType", fill = "lightgray", bins = 50, xlab = "The number of cell type", ylab = "The number of differential gene") + ggtitle("Differential-up genes(old/young)") + xlim(0, symmetric_limit) + theme(plot.title = element_text(hjust = 0.5))
p2 <- gghistogram(DEG.down.Freq[,c("gene", "n_sig_down_cellType")], x = "n_sig_down_cellType", fill = "lightgray", bins = 50, xlab = "The number of cell type", ylab = "The number of differential gene") + ggtitle("Differential-down genes(old/young)") + xlim(0, symmetric_limit) + theme(plot.title = element_text(hjust = 0.5))

#### The tissue frequency of differential genes
DEG.tissue.up <- DEG.up.df
DEG.tissue.up$tissue <- gsub(":.*", "", DEG.tissue.up$cellType)
DEG.tissue.up <- unique(DEG.tissue.up[,c("gene", "tissue")])
DEG.tissue.up.Freq <- DEG.tissue.up %>% distinct(gene, tissue) %>% count(gene, name = "sig_up_cellType_cover_n_tissue") %>% arrange(desc(sig_up_cellType_cover_n_tissue))
plot.up.Freq <- as.data.frame(table(DEG.tissue.up.Freq$sig_up_cellType_cover_n_tissue))
colnames(plot.up.Freq)[1] <- "category_var"
plot.up.Freq$category_var <- factor(plot.up.Freq$category_var, levels = as.character(sort(as.numeric(levels(plot.up.Freq$category_var)))))

DEG.tissue.down <- DEG.down.df
DEG.tissue.down$tissue <- gsub(":.*", "", DEG.tissue.down$cellType)
DEG.tissue.down <- unique(DEG.tissue.down[,c("gene", "tissue")])
DEG.tissue.down.Freq <- DEG.tissue.down %>% distinct(gene, tissue) %>% count(gene, name = "sig_down_cellType_cover_n_tissue") %>% arrange(desc(sig_down_cellType_cover_n_tissue))
plot.down.Freq <- as.data.frame(table(DEG.tissue.down.Freq$sig_down_cellType_cover_n_tissue))
colnames(plot.down.Freq)[1] <- "category_var"
plot.down.Freq$category_var <- factor(plot.down.Freq$category_var, levels = as.character(sort(as.numeric(levels(plot.down.Freq$category_var)))))

p11 <- ggbarplot(plot.up.Freq, x = "category_var", y = "Freq", xlab = "The number of tissue", ylab = "The number of differential gene", label = T, fill = "grey") + ggtitle("Differential-up genes(old/young)") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
p22 <- ggbarplot(plot.down.Freq, x = "category_var", y = "Freq", xlab = "The number of tissue", ylab = "The number of differential gene", label = T, fill = "grey") + ggtitle("Differential-down genes(old/young)") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)

pdf("2.rnaLevel/differentialGeneAnalysis/DEG.Frequency.pdf")
ggarrange(p1, p11, p2, p22, ncol = 2, nrow = 2, align = "hv")
dev.off()

#### plot top 30
pdf("2.rnaLevel/differentialGeneAnalysis/RNA.edgeR.DEGs.top30.pdf", width = 50, height = 15)
# plot LFC heatmap
Freq.mat <- unique(DEG.up.Freq[,c("gene", "n_sig_up_cellType")])[1:30,]
select.genes <- as.character(Freq.mat$gene)
p11 <- edgeR.LFC.heatmap(edgeR.result.list = RNA.edgeR.list, select.genes = select.genes, assay.info = PT.info.RNA.cellType, color.tissue = color.tissue, color.category = color.category, Freq.mat = Freq.mat, hgt = 15, wth = 80, title = "Top 30 frequent increased genes(old/young)", setNA = 0)
Freq.mat <- unique(DEG.down.Freq[,c("gene", "n_sig_down_cellType")])[1:30,]
select.genes <- as.character(Freq.mat$gene)
p22 <- edgeR.LFC.heatmap(edgeR.result.list = RNA.edgeR.list, select.genes = select.genes, assay.info = PT.info.RNA.cellType, color.tissue = color.tissue, color.category = color.category, Freq.mat = Freq.mat, hgt = 15, wth = 80, title = "Top 30 frequent decreased genes(old/young)", setNA = 0)

# split category
categorys <- as.character(unique(PT.info.RNA.cellType$category))
categorys <- categorys[grep("Others|Neural", categorys, invert = T)]
categorys.Freq <- lapply(categorys, function(x){
  cat("Analsis category:", x, "\n")
  up.res <- DEG.up.df
  names(up.res)[2] <- "tissue_cellType"
  select.category.info <- unique(PT.info.RNA.cellType[,c("tissue_cellType", "category")])
  select.category.info <- select.category.info[which(select.category.info$category == x),]
  up.res <- up.res[which(up.res$tissue_cellType %in% select.category.info$tissue_cellType),]
  up.genes <- unique(unlist(up.res))
  up.Freq <- up.res %>% distinct(gene, tissue_cellType) %>% group_by(gene) %>% summarise(n_sig_up_cellType = n(), sig_up_cellType = paste(tissue_cellType, collapse = ";"), .groups = "drop") %>% arrange(desc(n_sig_up_cellType))
  up.Freq.mat <- unique(up.Freq[,c("gene", "n_sig_up_cellType")])[1:30,]
  select.genes <- as.character(up.Freq.mat$gene)
  p111 <- edgeR.LFC.heatmap(edgeR.result.list = RNA.edgeR.list, select.genes = select.genes, assay.info = PT.info.RNA.cellType, color.tissue = color.tissue, color.category = color.category, Freq.mat = up.Freq.mat, hgt = 15, title = paste0(x, " :Top 30 frequent increased genes(old/young)"), setNA = 0)

  down.res <- DEG.down.df
  names(down.res)[2] <- "tissue_cellType"
  select.category.info <- unique(PT.info.RNA.cellType[,c("tissue_cellType", "category")])
  select.category.info <- select.category.info[which(select.category.info$category == x),]
  down.res <- down.res[which(down.res$tissue_cellType %in% select.category.info$tissue_cellType),]
  down.genes <- unique(unlist(down.res))
  down.Freq <- down.res %>% distinct(gene, tissue_cellType) %>% group_by(gene) %>% summarise(n_sig_down_cellType = n(), sig_down_cellType = paste(tissue_cellType, collapse = ";"), .groups = "drop") %>% arrange(desc(n_sig_down_cellType))
  down.Freq.mat <- unique(down.Freq[,c("gene", "n_sig_down_cellType")])[1:30,]
  select.genes <- as.character(down.Freq.mat$gene)
  p222 <- edgeR.LFC.heatmap(edgeR.result.list = RNA.edgeR.list, select.genes = select.genes, assay.info = PT.info.RNA.cellType, color.tissue = color.tissue, color.category = color.category, Freq.mat = down.Freq.mat, hgt = 15, title = paste0(x, " :Top 30 frequent decreased genes(old/young)"), setNA = 0)
  DEG.Freq <- merge(up.Freq, down.Freq, by = "gene", all = T)
  DEG.Freq$sig_up_cellType[is.na(DEG.Freq$sig_up_cellType)] <- ""
  DEG.Freq$sig_down_cellType[is.na(DEG.Freq$sig_down_cellType)] <- ""
  DEG.Freq[is.na(DEG.Freq)] <- 0
  DEG.Freq$n_sig_cellType <- rowSums(DEG.Freq[,c("n_sig_up_cellType", "n_sig_down_cellType")])
  DEG.Freq <- DEG.Freq[,c("gene", "n_sig_cellType", "n_sig_up_cellType", "sig_up_cellType", "n_sig_down_cellType", "sig_down_cellType")]
  DEG.Freq <- DEG.Freq[order(DEG.Freq$n_sig_up_cellType, decreasing = T),]
  return(DEG.Freq)
})
names(categorys.Freq) <- categorys
dev.off()

## output the result
DEG.Freq <- merge(DEG.up.Freq, DEG.down.Freq, by = "gene", all = T)
DEG.Freq$sig_up_cellType[is.na(DEG.Freq$sig_up_cellType)] <- ""
DEG.Freq$sig_down_cellType[is.na(DEG.Freq$sig_down_cellType)] <- ""
DEG.Freq[is.na(DEG.Freq)] <- 0
DEG.Freq$n_sig_cellType <- rowSums(DEG.Freq[,c("n_sig_up_cellType", "n_sig_down_cellType")])
DEG.Freq <- DEG.Freq[,c("gene", "n_sig_cellType", "n_sig_up_cellType", "sig_up_cellType", "n_sig_down_cellType", "sig_down_cellType")]
DEG.Freq <- DEG.Freq[order(DEG.Freq$n_sig_up_cellType, decreasing = T),]
DEG.Freq.list <- c(list(All = DEG.Freq), categorys.Freq)
saveRDS(DEG.Freq.list, "2.rnaLevel/differentialGeneAnalysis/RNA.edgeR.DEGs.Freq.rds")
writexl::write_xlsx(DEG.Freq.list, "2.rnaLevel/differentialGeneAnalysis/RNA.edgeR.DEGs.Freq.xlsx")
