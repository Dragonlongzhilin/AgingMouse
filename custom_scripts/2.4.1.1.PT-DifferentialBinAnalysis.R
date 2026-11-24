#' @description read change in bins among different age
#' @Author: Zhilin Long
#' @date 2026-6-5

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

histone <- "H3K4me1"
binsize <- "10kb"
bin.SAF <- paste0("~/longzhilin/Softwares/genomes/mouse/mm10/mm10_", binsize, "_bins.saf")
mm10.bin <- read.table(bin.SAF, sep = "\t", header = T, stringsAsFactors = F)
outDir <- paste0("~/longzhilin/Projects/AgingMouse/PairedTag/3.histoneLevel/", histone, "/bin", binsize, "/differentialBinAnalysis")
dir.create(outDir, recursive = T, showWarnings = F)
setwd(outDir)
log2FC.threshold <- log2(1.2)
sig.label <- "PValue"

# summarise by tissue, age and cell type
PT.info <- readRDS("~/longzhilin/Projects/AgingMouse/PairedTag/1.dataPreprocessing/PT.info.rds")
PT.info.cellType <- PT.info %>% dplyr::select(tissue_cellType, tissue_cellType_new, cellType, tissue, tissue_new, category, H3K4me1.depth) %>% group_by(tissue_cellType) %>% summarise(tissue_cellType_new = unique(tissue_cellType_new), cellType = unique(cellType), read = sum(H3K4me1.depth, na.rm = T), tissue = unique(tissue), tissue_new = unique(tissue_new), category = unique(category), .groups = "drop")

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

#### ---- 1.differential bin analysis ---- ####
#### feature counts
bin.count <- featureCounts(files = cellType.bam.files$bamFile, annot.inbuilt = "mm10", isGTFAnnotationFile = F, annot.ext = bin.SAF, nthreads = 20, fraction = T, countMultiMappingReads = T)
counts <- bin.count$counts
colnames(counts) <- cellType.bam.files$fileName
idx <- match(mm10.bin$GeneID, rownames(counts))
rownames(counts)[idx] <- paste0(mm10.bin$Chr, ":", mm10.bin$Start, "-", mm10.bin$End)
saveRDS(counts, "counts.rds")
cpms <- edgeR::cpm(counts)
saveRDS(cpms, "cpms.rds")

#### condition information
info <- data.frame(label = colnames(counts))
info$tissue <- gsub(":.*", "", info$label)
info$sample <- gsub(".*_mouse", "mouse", info$label)
info$sample <- gsub("_.*", "", info$sample)
info$batch <- gsub(".*_", "", info$label)
info$age_group <- gsub(".*_Young_.*", "Young", info$label)
info$age_group <- gsub(".*_Old_.*", "Old", info$age_group)
info$cellType <- gsub(".*:", "", info$label)
info$cellType <- gsub("_age.*|_Young.*|_Old.*", "", info$cellType)
info$tissue_cellType <- paste0(info$tissue, ":", info$cellType)
info <- info %>% left_join(PT.info.cellType[,c("tissue_cellType", "tissue_cellType_new", "tissue_new", "category", "read")], by = "tissue_cellType")
saveRDS(info, "condition.rds")

#### log2fc[cpm:Old/Young]
info$group <- paste0(info$tissue_cellType, "_", info$age_group)
cpms.mean.bin <- t(apply(cpms, 1, function(x) {
  tapply(x, info$group, mean)
}))
min_value <- min(cpms.mean.bin[cpms.mean.bin > 0], na.rm = TRUE)
cat("    min_value: ", min_value, "\n")
cpms_Old <- cpms.mean.bin[,grep("Old", colnames(cpms.mean.bin))] + 0.001
cpms_Young <- cpms.mean.bin[,grep("Young", colnames(cpms.mean.bin))] + 0.001
cpms.fc <- log2(cpms_Old / cpms_Young)
colnames(cpms.fc) <- gsub("_Old|_Young", "", colnames(cpms.fc))
saveRDS(cpm.log2fc, "cpm.log2fc.rds")

#### edgeR analysis
source("~/Functions/edgeR_diffRegions.R")
diffBin.list <- lapply(unique(info$tissue_cellType), function(tissue_cellType){
  cat("   ", tissue_cellType, "\n")
  info.subset <- info[which(info$tissue_cellType == tissue_cellType),]
  count.subset <- counts[,info.subset$label]
  cellType <- unique(info.subset$cellType)
  tissue <- unique(info.subset$tissue)
  if(length(unique(info.subset$batch)) > 1){
    batch <- info.subset$batch
  }else{
    batch <- NULL
  }
  outDir <- file.path("./", tissue)
  if(!dir.exists(outDir)){
    dir.create(outDir, recursive = TRUE)
  }
  group <- factor(info.subset$age, levels = c("Young", "Old"))
  title <- paste0(tissue, ":", cellType, "(depth[", unique(info.subset$read), "M])")
  diffBins <- edgeR.diffRegions(countTables = count.subset, group = group, batch = batch, outPath = outDir, type = cellType, title = title, logFC = log2FC.threshold, PValue = 0.05, ContinuousComparison = FALSE, cpm.filter = TRUE, min.regions = 10)
  diffBins <- diffBins$GroupComparison.Old.vs.Young
  if(!is.null(diffBins)){
    diffBins <- diffBins[grep("chrM:", diffBins$Region, invert = TRUE),] # remove chrM

    # output significant bins
    idx <- which(diffBins[,sig.label] < 0.05 & abs(diffBins$logFC) > log2FC.threshold)
    if(length(idx) > 0){
      diffBins.sig <- diffBins[idx,]
    }else{
      diffBins.sig <- diffBins
    }
    write.csv(diffBins.sig, paste0(outDir, "/", tissue_cellType, ".sigDiffBin.csv"), row.names = FALSE)
  }
  return(diffBins)
})
names(diffBin.list) <- unique(info$tissue_cellType)
saveRDS(diffBin.list, "diffBin.list.rds")

#### ---- Frequency differential bins ---- ####
diffBin.list <- readRDS("diffBin.list.rds")
#### The celltype frequency of differential bins
## Up-regulated bins: extract significant bins with cell type information
diffBin.up.df <- lapply(names(diffBin.list), function(ct){
  x <- diffBin.list[[ct]]
  idx <- which(x[, sig.label] < 0.05 & x$logFC > log2FC.threshold)
  if(length(idx) > 0){
    return(data.frame(Region = x$Region[idx], cellType = ct, stringsAsFactors = FALSE))
  }else{
    return(NULL)
  }
})
diffBin.up.df <- do.call(rbind, diffBin.up.df)
diffBin.up.Freq <- diffBin.up.df %>% group_by(Region) %>% summarise(n_sig_up_cellType = n(), sig_up_cellType = paste(cellType, collapse = ";"), .groups = "drop") %>% arrange(desc(n_sig_up_cellType))

## Down-regulated bins: extract significant bins with cell type information
diffBin.down.df <- lapply(names(diffBin.list), function(ct){
  x <- diffBin.list[[ct]]
  idx <- which(x[, sig.label] < 0.05 & x$logFC < (log2FC.threshold * -1))
  if(length(idx) > 0){
    return(data.frame(Region = x$Region[idx], cellType = ct, stringsAsFactors = FALSE))
  }else{
    return(NULL)
  }
})
diffBin.down.df <- do.call(rbind, diffBin.down.df)
diffBin.down.Freq <- diffBin.down.df %>% group_by(Region) %>% summarise(n_sig_down_cellType = n(), sig_down_cellType = paste(cellType, collapse = ";"), .groups = "drop") %>% arrange(desc(n_sig_down_cellType))

symmetric_limit <- max(max(diffBin.up.Freq$n_sig_up_cellType), max(diffBin.down.Freq$n_sig_down_cellType))
p1 <- gghistogram(diffBin.up.Freq, x = "n_sig_up_cellType", fill = "lightgray", bins = 50, xlab = "The number of cell type", ylab = "The number of differential bin") + ggtitle("Differential-up bins(Old/Young)") + xlim(0, symmetric_limit) + theme(plot.title = element_text(hjust = 0.5))
p2 <- gghistogram(diffBin.down.Freq, x = "n_sig_down_cellType", fill = "lightgray", bins = 50, xlab = "The number of cell type", ylab = "The number of differential bin") + ggtitle("Differential-down bins(Old/Young)") + xlim(0, symmetric_limit) + theme(plot.title = element_text(hjust = 0.5))

#### The tissue frequency of differential bins
diffBin.tissue.up <- lapply(names(diffBin.list), function(y){
  x <- diffBin.list[[y]]
  idx <- which(x[, sig.label] < 0.05 & x$logFC > log2FC.threshold)
  if(length(idx) > 0){
    data.frame(Region = x$Region[idx], tissue_cellType = y, stringsAsFactors = FALSE)
  }else{
    NULL
  }
})
diffBin.tissue.up <- do.call(rbind, diffBin.tissue.up)
diffBin.tissue.up$tissue <- gsub(":.*", "", diffBin.tissue.up$tissue_cellType)
diffBin.tissue.up.Freq <- diffBin.tissue.up %>% distinct(Region, tissue) %>% count(Region, name = "Freq") %>% arrange(desc(Freq))
plot.up.Freq <- diffBin.tissue.up.Freq %>% count(Freq) %>% select(category_var = Freq, Freq = n) %>% mutate(category_var = factor(category_var, levels = as.character(sort(unique(category_var)))))

diffBin.tissue.down <- lapply(names(diffBin.list), function(y){
  x <- diffBin.list[[y]]
  idx <- which(x[, sig.label] < 0.05 & x$logFC < (log2FC.threshold * -1))
  if(length(idx) > 0){
    data.frame(Region = x$Region[idx], tissue_cellType = y, stringsAsFactors = FALSE)
  }else{
    NULL
  }
})
diffBin.tissue.down <- do.call(rbind, diffBin.tissue.down)
diffBin.tissue.down$tissue <- gsub(":.*", "", diffBin.tissue.down$tissue_cellType)
diffBin.tissue.down.Freq <- diffBin.tissue.down %>% distinct(Region, tissue) %>% count(Region, name = "Freq") %>% arrange(desc(Freq))
plot.down.Freq <- diffBin.tissue.down.Freq %>% count(Freq) %>% select(category_var = Freq, Freq = n) %>% mutate(category_var = factor(category_var, levels = as.character(sort(unique(category_var)))))

p11 <- ggbarplot(plot.up.Freq, x = "category_var", y = "Freq", xlab = "The number of tissue", ylab = "The number of differential bin", label = T, fill = "grey") + ggtitle("Differential-up bins(Old/Young)") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
p22 <- ggbarplot(plot.down.Freq, x = "category_var", y = "Freq", xlab = "The number of tissue", ylab = "The number of differential bin", label = T, fill = "grey") + ggtitle("Differential-down bins(Old/Young)") + theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)

pdf("diffBin.Frequency.pdf")
ggarrange(p1, p11, p2, p22, ncol = 2, nrow = 2, align = "hv")
dev.off()

colnames(diffBin.tissue.up.Freq) <- c("Region", "sig_up_cellType_cover_n_tissue")
colnames(diffBin.tissue.down.Freq) <- c("Region", "sig_down_cellType_cover_n_tissue")
diffBin.Frequency <- list(cellType.up = diffBin.up.Freq, cellType.down = diffBin.down.Freq, tissue.up = diffBin.tissue.up.Freq, tissue.down = diffBin.tissue.down.Freq)
diffBin.Frequency <- Reduce(function(x,y) merge(x, y, all = T, by = "Region"), diffBin.Frequency)
diffBin.Frequency$sig_up_cellType[is.na(diffBin.Frequency$sig_up_cellType)] <- ""
diffBin.Frequency$sig_down_cellType[is.na(diffBin.Frequency$sig_down_cellType)] <- ""
diffBin.Frequency[is.na(diffBin.Frequency)] <- 0
diffBin.Frequency$n_sig_cellType <- rowSums(diffBin.Frequency[,c("n_sig_up_cellType", "n_sig_down_cellType")])
diffBin.Frequency$sig_cellType_cover_n_tissue <- rowSums(diffBin.Frequency[,c("sig_up_cellType_cover_n_tissue", "sig_down_cellType_cover_n_tissue")])
diffBin.Frequency <- diffBin.Frequency[,c("Region", "n_sig_cellType", "n_sig_up_cellType", "sig_up_cellType", "n_sig_down_cellType", "sig_down_cellType", "sig_cellType_cover_n_tissue", "sig_up_cellType_cover_n_tissue", "sig_down_cellType_cover_n_tissue")]
diffBin.Frequency <- diffBin.Frequency[order(diffBin.Frequency$n_sig_up_cellType, decreasing = T),]
saveRDS(diffBin.Frequency, "diffBin.Frequency.rds")
writexl::write_xlsx(diffBin.Frequency, "diffBin.Frequency.xlsx")

## plot the frequency of different peak
# Count overlapping points at each (x, y) coordinate for bubble plot
plot.data <- diffBin.Frequency
plot.data.count <- plot.data %>% group_by(n_sig_up_cellType, n_sig_down_cellType) %>% summarise(n = n(), .groups = "drop") %>% mutate(log10_n = log10(n))
top_bins <- plot.data %>% slice_max(n_sig_up_cellType, n = 5, with_ties = FALSE) %>% bind_rows(plot.data %>% slice_max(n_sig_down_cellType, n = 5, with_ties = FALSE)) %>% distinct(Region)
xy_limit <- max(plot.data[,c("n_sig_up_cellType", "n_sig_down_cellType")], na.rm = T)
points_to_label <- subset(plot.data, Region %in% top_bins$Region)
p <- ggplot(plot.data.count, aes(x = n_sig_up_cellType, y = n_sig_down_cellType)) +
  geom_count(aes(size = log10_n, color = log10_n), alpha = 0.7) +
  scale_size_continuous(name = "log10 (# of bins)", range = c(1, 3)) +
  geom_text_repel(data = points_to_label, 
                  aes(x = n_sig_up_cellType, y = n_sig_down_cellType, label = Region), 
                  size = 2, box.padding = 0.5, point.padding = 0.3, 
                  force = 2, max.overlaps = 50) +
  scale_color_viridis_c(option = "D", name = "log10 (# of bins)") + 
  xlim(0, xy_limit) + ylim(0, xy_limit) +
  ggtitle(paste0("Total differential bins: ", nrow(plot.data))) + 
  xlab(paste0("# of cell types with increased ", histone)) + 
  ylab(paste0("# of cell types with decreased ", histone)) +
  theme_classic() + 
  theme(
    aspect.ratio = 1, 
    plot.title = element_text(hjust = 0.5)
  )
pdf("diffBin.Frequency.bubblePlot.pdf", width = 5, height = 5)
print(p)
dev.off()

