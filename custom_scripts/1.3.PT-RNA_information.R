#' @description cell type annotation for brain tissue
#' @Author: Zhilin Long
#' @date 2026-6-10

library(Seurat)
library(BPCells)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(anndataR)
set.seed(101)
require(future)
plan("multicore", workers = 5)
options(future.globals.maxSize = 100000 * 1024^2)
set.resolutions <- seq(0.1, 1.2, by = 0.1)
source("~/Functions/colorAssign.R")
setwd("~/AgingMouse/PairedTag/1.dataPreprocessing")
tissue.seurat.final <- readRDS("tissue.seurat.final.rds")

#### ---- read depth and split category ---- ####
#### count the depth of RNA and DNA
depth.list <- lapply(names(tissue.seurat.final), function(tissue){
  x <- tissue.seurat.final[[tissue]]
  x$batch <- gsub("-.*", "", x$assay)
  x$group <- paste0(x$cellType, "_", x$age_group, "_", x$sample, "_", x$batch)
  groups <- unique(x$group)
  RNA.depth <- round(tapply(x$nCount_RNA, x$group, sum) / 1000000, 2)
  RNA.cells <- table(x$group)
  depth.info <- data.frame(tissue = tissue, group = groups, RNA.depth = as.numeric(RNA.depth[groups]), RNA.cells = as.numeric(RNA.cells[groups]))
  depth.info$cellType <- gsub("_Old.*|_Young.*", "", depth.info$group)
  depth.info$sample <- gsub(".*_mouse", "mouse", depth.info$group)
  depth.info$batch <- gsub(".*_", "", depth.info$sample)
  depth.info$sample <- gsub("_.*", "", depth.info$sample)
  depth.info$age_group <- gsub(".*_Old.*", "Old", depth.info$group)
  depth.info$age_group <- gsub(".*_Young.*", "Young", depth.info$age_group)

  histones <- as.character(unique(x$histone))
  DNA.depths <- lapply(histones, function(y){
    histone.metadata <- x@meta.data[which(x$histone == y),]
    histone.depth <- round(tapply(histone.metadata$nCount_DNA, histone.metadata$group, sum) / 1000000, 2)
    histone.cells <- table(histone.metadata$group)
    histone.info <- data.frame(tissue = tissue, group = names(histone.cells), depth = as.numeric(histone.depth[names(histone.cells)]), cells = as.numeric(histone.cells))
    colnames(histone.info) <- c("tissue", "group", paste0(y, ".depth"), paste0(y, ".cells"))
    return(histone.info)
  })
  DNA.depths <- Reduce(function(x,y) merge(x,y, by = c("tissue", "group"), all = TRUE), DNA.depths)
  depth.info <- merge(depth.info, DNA.depths, by = c("tissue", "group"), all = TRUE)
  depth.info <- depth.info[order(depth.info$RNA.depth, decreasing = T),]
  return(depth.info)
})
PT.info <- do.call(rbind, depth.list)
PT.info$group <- paste0(PT.info$tissue, ":", PT.info$group)

#### add category
# ref: /mnt/transposon1/zhangyanxiaoLab/longzhilin/Data/AgingMouse/singlecell/Altas/Zhang.2024.Science.EasySci/Table S2 - Metadata for main cell types annotated in each organ or tissue.csv
PT.info$category <- "Others" # Choroid_plexus/Enteric_neurons
PT.info$category[which(PT.info$cellType %in% c("Oligodendrocytes", "Astrocytes", "Bergmann glia", "OPC", "Ependymal cells", "Enteric glia", "Schwann cells"))] <- "Glial"
PT.info$category[grep(" Gaba", PT.info$cellType)] <- "GABAergic"
PT.info$category[grep(" Glut", PT.info$cellType)] <- "Glutamatergic"
PT.info$category[which(PT.info$cellType %in% c("B cells", "T cells", "Plasma cells"))] <- "Lymphoid"
PT.info$category[which(PT.info$cellType %in% c("Macrophages", "Microglia", "Kupffer cells", "Alveolar macrophages", "Interstitial macrophages", "Monocytes", "Myeloid cells", "Dendritic cells", "Neutrophils", "Microglia/BAM", "Erythroid cells", "Langerhans cells"))] <- "Myeloid"
PT.info$category[which(PT.info$cellType %in% c("Type I alveolar epithelial cells", "Type II alveolar epithelial cells", "Secretory cells", "Ciliated cells", "Urothelial cells", "PT-S1", "PT-S2", "PT-S3", "CD-PC", "CD-IC", "LOH-AL", "LOH-DL", "DCT", "Podocytes", "Colonic epithelial cells", "Goblet cells", "Enteroendocrine cells", "Tuft cells", "Cholangiocytes", "Hepatocytes", "Mesothelial cells", "Juxtaglomerular cells", "Sertoli cells", "Choroid plexus", "IFE basal", "IFE suprabasal", "Hair follicle stem cells", "Upper hair follicle cells", "IFE cycling basal", "Sebocytes"))] <- "Epithelial"
PT.info$category[which(PT.info$cellType %in% c("Vascular endothelial cells", "Lymphatic endothelial cells", "Aerocytes", "General capillary cells", "Endocardial endothelial cells"))] <- "Endothelial"
PT.info$category[which(PT.info$cellType %in% c("Fibroblasts", "Mural cells", "Mfge8 positive cells", "Ngf positive cells", "Tenocytes", "Fibro-adipogenic progenitors", "Asic2 positive cells", "Hepatic stellate cells", "Pericytes", "Ncam1 positive cells", "VLMC", "ABC", "Leydig cells", "Myoid cells", "Telocytes", "Dermal papilla cells", "Dermal sheath cells"))] <- "Stromal"
PT.info$category[which(PT.info$cellType %in% c("Neuromuscular junction myonuclei", "Myotendinous junction myonuclei", "Skeletal muscle cells", "Ventricular cardiomyocytes", "Muscle satellite cells", "Myonuclei"))] <- "Muscle"
PT.info$category[which(PT.info$cellType %in% c("Adipocytes", "Brown adipocytes", "ASPC"))] <- "Adipose"
PT.info$category[which(PT.info$cellType %in% c("Spermatogonia", "Spermatocytes", "Early round spermatids", "Late round spermatids", "Elongating spermatids"))] <- "Germline"
PT.info$cell_type <- gsub(" |-", "_", PT.info$cellType)
PT.info$cell_type <- gsub("\\/", "", PT.info$cell_type)
PT.info$category <- factor(PT.info$category, levels = names(color.category))
PT.info$tissue <- factor(PT.info$tissue, levels = names(color.tissue))
PT.info$tissue_cellType <- paste0(PT.info$tissue, ":", PT.info$cell_type)

## add color for cell type
load("~/PanSci_Color_Palette.RData")
color_palette <- unique(unlist(main_cluster_palette_list))
set.seed(123)
color_palette <- sample(color_palette)
total_cellTypes <- unique(PT.info$cell_type)
cellType_colors <- color_palette[1:length(total_cellTypes)]
PT.info$colors <- plyr::mapvalues(PT.info$cell_type, from = total_cellTypes, to = cellType_colors)
PT.info$category.colors <- color.category[as.character(PT.info$category)]
PT.info$tissue_new <- recode(as.character(PT.info$tissue), "brainCB" = "Cerebellum", "brainFC" = "Frontal cortex", "brainHip" = "Hippocampus", "heart" = "Heart","liver" = "Liver", "lung" = "Lung", "kidney" = "Kidney", "BAT" = "BAT", "muscle" = "Muscle", "spleen" = "Spleen", "colon" = "Colon", "skin" = "Skin", "testis" = "Testis")
PT.info$tissue_new <- factor(PT.info$tissue_new, levels = names(color.tissue2))

# order cell type with tissue and RNA depth
PT.info.RNA.cellType <- PT.info %>% dplyr::select(tissue_cellType, cellType, tissue, category, RNA.depth) %>% group_by(tissue_cellType) %>% summarise(cellType = unique(cellType), read = sum(RNA.depth, na.rm = T), tissue = unique(tissue), category = unique(category), .groups = "drop")
PT.info.RNA.cellType$cellType <- as.character(PT.info.RNA.cellType$cellType)
PT.info.RNA.cellType <- arrange(PT.info.RNA.cellType, category, tissue, cellType, desc(read))
cellType.order <- unique(PT.info.RNA.cellType$cellType)
PT.info$cellType <- factor(PT.info$cellType, levels = cellType.order)
cellType_new <- paste0(seq_along(cellType.order), ".", cellType.order)
PT.info$cellType_new <- plyr::mapvalues(PT.info$cellType, from = cellType.order, to = cellType_new)
PT.info$cellType_new <- factor(PT.info$cellType_new, levels = cellType_new)
PT.info$tissue_cellType_new <- paste0(PT.info$tissue_new, ":", PT.info$cell_type)
saveRDS(PT.info, "PT.info.rds")

#### update the read depth: 2026.8.15
# after remove duplicated reads
bam_read_counts <- read.table("~/AgingMouse/PairedTag_new/bam_read_counts.txt", header = T, sep = "\t", stringsAsFactors = F)
bam_read_counts$dir_type <- "cellType"
bam_read_counts$dir_type[grep("cellType_age/", bam_read_counts$bam_file)] <- "cellType_age"
bam_read_counts$dir_type[grep("cellType_age/RNA", bam_read_counts$bam_file)] <- "cellType_age_RNA"
bam_read_counts$dir_type[grep("cellType_age_replication", bam_read_counts$bam_file)] <- "cellType_age_replication"
bam_read_counts$label <- basename(bam_read_counts$bam_file)
bam_read_counts$label <- gsub(".bam", "", bam_read_counts$label)
bam_read_counts$histone <- gsub("\\..*", "", bam_read_counts$label)
bam_read_counts$cell_type <- gsub(".*\\.", "", bam_read_counts$label)
bam_read_counts$bam_file <- gsub("~/AgingMouse/PairedTag_new/", "", bam_read_counts$bam_file)
bam_read_counts$tissue <- gsub("/.*", "", bam_read_counts$bam_file)
saveRDS(bam_read_counts, "bam_read_counts.rds")

# add to PT.info
bam_read_counts <- readRDS("bam_read_counts.rds")
DNA_read_counts <- bam_read_counts[which(bam_read_counts$dir_type == "cellType_age_replication"),]
DNA_read_counts$cell_type <- gsub("_Old.*|_Young.*", "", gsub(".*\\.", "", DNA_read_counts$label))
DNA_read_counts$age <- gsub(".*_Old_.*", "Old", DNA_read_counts$label)
DNA_read_counts$age <- gsub(".*_Young_.*", "Young", DNA_read_counts$age)
DNA_read_counts$sample <- gsub(".*mouse", "mouse", DNA_read_counts$label)
DNA_read_counts$batch <- gsub(".*_", "", DNA_read_counts$sample)
DNA_read_counts$sample <- gsub("_.*", "", DNA_read_counts$sample)
DNA_read_counts$depth <- round(DNA_read_counts$read_count / 1e6, 2)
DNA_read_counts$key <- paste(DNA_read_counts$tissue, DNA_read_counts$cell_type, DNA_read_counts$age, DNA_read_counts$sample, DNA_read_counts$batch, sep = "|")
PT.info$key <- paste(as.character(PT.info$tissue), PT.info$cell_type, PT.info$age_group, PT.info$sample, PT.info$batch, sep = "|")
for (h in unique(DNA_read_counts$histone)) {
  info <- DNA_read_counts[which(DNA_read_counts$histone == h),]
  PT.info[[paste0(h, ".depth")]] <- info$depth[match(PT.info$key, info$key)]
}
PT.info$key <- NULL
saveRDS(PT.info, "PT.info.rds")

#### output the celltype information
PT.info.cellType <- PT.info %>% group_by(tissue_cellType_new) %>% summarise(tissue = unique(tissue_new), cellType = unique(cellType), category = unique(category), RNA.cell = sum(RNA.cells, na.rm = T), RNA.read = sum(RNA.depth, na.rm = T), H3K9me3.cell = sum(H3K9me3.cells, na.rm = T), H3K9me3.depth = sum(H3K9me3.depth, na.rm = T), H3K27me3.cell = sum(H3K27me3.cells, na.rm = T), H3K27me3.depth = sum(H3K27me3.depth, na.rm = T), H3K36me3.cell = sum(H3K36me3.cells, na.rm = T), H3K36me3.depth = sum(H3K36me3.depth, na.rm = T), H3K27ac.cell = sum(H3K27ac.cells, na.rm = T), H3K27ac.depth = sum(H3K27ac.depth, na.rm = T), H3K4me1.cell = sum(H3K4me1.cells, na.rm = T), H3K4me1.depth = sum(H3K4me1.depth, na.rm = T), H3K4me3.cell = sum(H3K4me3.cells, na.rm = T), H3K4me3.depth = sum(H3K4me3.depth, na.rm = T), colors = unique(colors), category.colors = unique(category.colors), .groups = "drop")
PT.info.cellType <- arrange(PT.info.cellType, tissue, desc(RNA.read))
saveRDS(PT.info.cellType, "PT.info.cellType.rds")
writexl::write_xlsx(PT.info.cellType, "PT.info.cellType.xlsx")

#### ---- plot Cell number for histone mark of each tissue ---- ####
PT.info.sample <- PT.info %>% group_by(tissue, sample, batch) %>% summarise(tissue_new = unique(tissue_new), age_group = unique(age_group), RNA.cell = sum(RNA.cells, na.rm = T), RNA.read = sum(RNA.depth, na.rm = T), H3K9me3.cell = sum(H3K9me3.cells, na.rm = T), H3K9me3.depth = sum(H3K9me3.depth, na.rm = T), H3K27me3.cell = sum(H3K27me3.cells, na.rm = T), H3K27me3.depth = sum(H3K27me3.depth, na.rm = T), H3K36me3.cell = sum(H3K36me3.cells, na.rm = T), H3K36me3.depth = sum(H3K36me3.depth, na.rm = T), H3K27ac.cell = sum(H3K27ac.cells, na.rm = T), H3K27ac.depth = sum(H3K27ac.depth, na.rm = T), H3K4me1.cell = sum(H3K4me1.cells, na.rm = T), H3K4me1.depth = sum(H3K4me1.depth, na.rm = T), H3K4me3.cell = sum(H3K4me3.cells, na.rm = T), H3K4me3.depth = sum(H3K4me3.depth, na.rm = T), .groups = "drop")

histone_age_sample_replication <- lapply(unique(PT.info.sample$tissue), function(tissue){
  info.subset <- PT.info.sample[which(PT.info.sample$tissue == tissue),]
  tissue_new <- unique(info.subset$tissue_new)
  histone.info <- lapply(names(color.histone), function(histone){
    info <- unique(na.omit(info.subset[,c("sample", "batch", "age_group", paste0(histone, ".cell"))]))
    info <- info[which(info[,4] > 0),]
    info <- info %>% group_by(age_group) %>% mutate(replication = paste0("rep", row_number())) %>% ungroup()
    colnames(info)[4] <- "cellNumber"
    info$histone <- histone
    return(info)
  })
  histone.info <- do.call(rbind, histone.info)
  histone.info$tissue <- tissue
  histone.info$tissue_new <- tissue_new
  return(histone.info)
})
histone_age_sample_replication <- do.call(rbind, histone_age_sample_replication)
histone_age_sample_replication$replication <- paste0(histone_age_sample_replication$age_group, "_", histone_age_sample_replication$replication)
histone_age_sample_replication$age_group <- factor(histone_age_sample_replication$age_group, levels = names(color.age_group))
histone_age_sample_replication$histone <- factor(histone_age_sample_replication$histone, levels = names(color.histone))
histone_age_sample_replication$tissue <- factor(histone_age_sample_replication$tissue, levels = names(color.tissue))
histone_age_sample_replication$tissue_new <- factor(histone_age_sample_replication$tissue_new, levels = names(color.tissue2))
histone_age_sample_replication <- arrange(histone_age_sample_replication, age_group, histone, replication)
histone_age_sample_replication$replication <- factor(histone_age_sample_replication$replication, levels = unique(histone_age_sample_replication$replication))
saveRDS(histone_age_sample_replication, "histone_age_sample_replication.rds")

## output the sample information
sample.information <- histone_age_sample_replication[,c("sample", "tissue_new", "age_group", "histone", "cellNumber", "batch")]
sample.information$tissue_new <- factor(sample.information$tissue_new, levels = names(color.tissue2))
sample.information <- arrange(sample.information, tissue_new, histone, age_group, batch)
sample.information$tissue_histone <- paste0(sample.information$tissue_new, ":", sample.information$histone)
sample.information$age_sample <- paste0(sample.information$age_group, ":", sample.information$sample, "_", sample.information$batch)
sample.information <- sample.information %>% dplyr::select(tissue_histone, age_sample, cellNumber) %>% pivot_wider(names_from = age_sample, values_from = cellNumber, values_fill = NA)
write.csv(sample.information, "PT.info.sample.csv", na = "NA", row.names = FALSE, quote = FALSE)

#### plots
pdf("tissue_sample_cellNumber3.pdf", height = 4, width = 5)
plot.data$age_group <- factor(plot.data$age_group, levels = names(color.age_group))
ggballoonplot(
  plot.data,
  x = "age_group", y = "tissue_new",
  color = "age_group", fill = "age_group", palette = color.age_group,
  size = "cellNumber", size.range = c(3, 5)
) +
  facet_wrap(~ histone, nrow = 1) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "gray90", color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.spacing = unit(0, "mm")
  ) +
  scale_size_continuous(
    breaks = pretty(plot.data$cellNumber, n = 4),
    range  = c(3, 5),
    limits = range(plot.data$cellNumber),
  )
dev.off()

#### ---- batch information ---- ####
histone_age_sample_replication <- readRDS("histone_age_sample_replication.rds")
sample.info <- lapply(names(tissue.seurat.final), function(tissue){
  seuratObj <- tissue.seurat.final[[tissue]]
  seuratObj$batch <- gsub("-.*", "", seuratObj$assay)
  condition <- unique(seuratObj@meta.data[,c("tissue", "assay", "batch", "sample", "age", "age_group", "histone")])
  rownames(condition) <- NULL
  return(condition)
})
sample.info <- do.call(rbind, sample.info)
sample.info <- left_join(sample.info, histone_age_sample_replication, by = c("tissue", "sample", "age_group", "histone", "batch"))
sample.info$assay <- gsub("scPT-Natural-Aging-", "", sample.info$assay)
sample.info$tissue <- factor(sample.info$tissue, levels = names(color.tissue))
sample.info$tissue_new <- factor(sample.info$tissue_new, levels = names(color.tissue2))
sample.info$histone <- factor(sample.info$histone, levels = names(color.histone))
sample.info <- arrange(sample.info, tissue, histone, age_group, replication)
sample.info$sample <- gsub("\\..*", "", sample.info$sample)
saveRDS(sample.info, file = "sample.info.rds")
write.csv(sample.info, file = "sample.info.csv", row.names = FALSE, quote = FALSE)

#### ---- UMAP plot for each tissue ---- ####
PT.info.cellTypeMarker <- readxl::read_xlsx("PT.info.cellTypeMarker.xlsx")
plot.list <- lapply(names(tissue.seurat.final), function(tissue){
  PairedTag_data <- tissue.seurat.final[[tissue]]
  info.subset <- unique(PT.info[which(PT.info$tissue == tissue),c("cellType", "colors", "tissue_new")])
  tissue_new <- unique(info.subset$tissue_new)
  color.cellType <- info.subset$colors
  names(color.cellType) <- info.subset$cellType
  p1 <- DimPlot(PairedTag_data, reduction = "umap", group.by = "cellType", cols = color.cellType, label = T, repel = T, label.size = 5, raster=FALSE) + xlab("UMAP_1") + ylab("UMAP_2") + ggtitle(paste0(tissue_new, ": ", format(ncol(PairedTag_data), big.mark = ","), " cells")) + guides(color = guide_legend(title = "Cell type", ncol = 1, override.aes = list(size = 3))) + theme(aspect.ratio = 1, legend.position = "right", plot.title = element_text(size = 16, hjust = 0.5, face = "plain"), axis.text = element_text(size = 12), axis.title = element_text(size = 14))

  cellType.order <- sort(levels(PairedTag_data$cellType))
  geneMarkers <- PT.info.cellTypeMarker[which(PT.info.cellTypeMarker$tissue == tissue_new),]
  geneMarkers <- geneMarkers[match(cellType.order, geneMarkers$cellType),]
  genes <- unique(unlist(strsplit(geneMarkers$GeneMarkers, ", ")))
  PairedTag_data$cellType <- factor(PairedTag_data$cellType, levels = cellType.order)
  p2 <- DotPlot(PairedTag_data, features = genes, group.by = "cellType", cols = c("#1e90ff", "#ff5a36"), dot.scale = 4) + xlab("") + ylab("") + ggtitle(tissue_new) + theme(axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5), plot.title = element_text(hjust = 0.5, size = 14))
  return(list(p1 = p1, p2 = p2))
})
names(plot.list) <- names(tissue.seurat.final)
png("Final.tissue_umap.png", width = 30, height = 35, units = "in", res = 300)
ggarrange(plot.list[[1]][[1]], plot.list[[2]][[1]], plot.list[[3]][[1]], plot.list[[4]][[1]], plot.list[[5]][[1]], plot.list[[6]][[1]], plot.list[[7]][[1]], plot.list[[8]][[1]], plot.list[[9]][[1]], plot.list[[10]][[1]], plot.list[[11]][[1]], plot.list[[13]][[1]], plot.list[[12]][[1]], nrow = 5, ncol = 3, align = "hv")
dev.off()
pdf("Final.cellMarker.DotPlot.pdf", height = 8, width = 20)
res <- lapply(plot.list, function(x) print(x[[2]]))
dev.off()

#### ---- plot merge depth ---- ####
PT.info <- readRDS("PT.info.rds")
histone_age_sample_replication <- readRDS("histone_age_sample_replication.rds")

#### tissue
## RNA
PT.info.tissue <- PT.info %>% group_by(tissue, sample, age_group, batch) %>%
  summarise(tissue_new = unique(tissue_new), RNA.depth = sum(RNA.depth, na.rm = T), RNA.cells = sum(RNA.cells, na.rm = T), H3K27ac.depth = sum(H3K27ac.depth, na.rm = T), H3K27me3.depth = sum(H3K27me3.depth, na.rm = T), H3K4me1.depth = sum(H3K4me1.depth, na.rm = T), H3K4me3.depth = sum(H3K4me3.depth, na.rm = T), H3K36me3.depth = sum(H3K36me3.depth, na.rm = T), H3K9me3.depth = sum(H3K9me3.depth, na.rm = T), .groups = "drop")
PT.info.tissue[PT.info.tissue == 0] <- NA
# set y axis scale
DNA.volume <- PT.info.tissue[, grep("depth", colnames(PT.info.tissue))]
DNA.volume <- DNA.volume[,grep("RNA.depth", colnames(DNA.volume), invert = T)]
y.max <- rowsum(DNA.volume, group = PT.info.tissue$tissue, na.rm = TRUE)
y.max <- max(unlist(y.max))

PT.info.tissue.RNA <- PT.info.tissue[,c("tissue", "sample", "age_group", "batch", "tissue_new", "RNA.depth", "RNA.cells")]
PT.info.tissue.RNA <- PT.info.tissue.RNA %>% left_join(histone_age_sample_replication[,c("tissue", "sample", "age_group", "batch", "histone")], by = c("tissue", "sample", "age_group", "batch")) %>% group_by(tissue, sample, age_group, batch, tissue_new, RNA.depth, RNA.cells) %>% summarise(histone = paste(histone, collapse = "-"), .groups = "drop")
PT.info.tissue.RNA$age_group <- factor(PT.info.tissue.RNA$age_group, levels = names(color.age_group))
PT.info.tissue.RNA <- arrange(PT.info.tissue.RNA, tissue, batch, age_group)
PT.info.tissue.RNA <- PT.info.tissue.RNA %>% group_by(tissue, age_group) %>% mutate(replication = paste0("rep", row_number())) %>% ungroup()
PT.info.tissue.RNA$replication <- paste0(PT.info.tissue.RNA$age_group, "_", PT.info.tissue.RNA$replication)

# set replication palette
info.old.rep <- gtools::mixedsort(as.character(unique(PT.info.tissue.RNA$replication[which(PT.info.tissue.RNA$age_group == "Old")])))
info.young.rep <- gtools::mixedsort(as.character(unique(PT.info.tissue.RNA$replication[which(PT.info.tissue.RNA$age_group == "Young")])))
old.rep.colors <- colorRampPalette(color.olds)(length(info.old.rep))
young.rep.colors <- colorRampPalette(color.youngs)(length(info.young.rep))
color.replications <- c(young.rep.colors, old.rep.colors)
names(color.replications) <- c(info.young.rep, info.old.rep)
PT.info.tissue.RNA$replication <- factor(PT.info.tissue.RNA$replication, levels = names(color.replications))
saveRDS(PT.info.tissue.RNA, "RNA_age_sample_replication.rds")

library(cowplot)
library(gridExtra)
PT.info.tissue.RNA$tissue_new <- factor(PT.info.tissue.RNA$tissue_new, levels = rev(names(color.tissue2)))
total_cells <- PT.info.tissue.RNA %>% group_by(tissue_new) %>% summarise(total = sum(RNA.cells, na.rm = TRUE))
p1 <- ggbarplot(PT.info.tissue.RNA, x = "tissue_new", y = "RNA.cells", xlab = "", ylab = "Cell number", color = "replication", fill = "replication", palette = color.replications, label = F, orientation = "horiz", position = position_stack(reverse = TRUE)) + ggtitle("Cell number") + theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5)) + geom_text(data = total_cells, aes(x = tissue_new, y = total, label = total), inherit.aes = FALSE, hjust = -0.2, size = 3.5)
total_depth <- PT.info.tissue.RNA %>% group_by(tissue_new) %>% summarise(total = sum(RNA.depth, na.rm = TRUE))
p2 <- ggbarplot(PT.info.tissue.RNA, x = "tissue_new", y = "RNA.depth", xlab = "", ylab = "Reads (M)", color = "replication", fill = "replication", palette = color.replications, label = F, orientation = "horiz", position = position_stack(reverse = TRUE)) + ggtitle("RNA") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position = "right", plot.title = element_text(size = 12, hjust = 0.5)) + geom_text(data = total_depth, aes(x = tissue_new, y = total, label = paste0(round(total, 2), "M")), inherit.aes = FALSE, hjust = -0.2, size = 3.5)
aligned_plots1 <- align_plots(p1, p2, align = "h",axis = "tb") 

## DNA
histone_age_sample_replication <- readRDS("histone_age_sample_replication.rds")

# set replication palette
info.old.rep <- unique(histone_age_sample_replication$replication[which(histone_age_sample_replication$age_group == "Old")])
info.young.rep <- unique(histone_age_sample_replication$replication[which(histone_age_sample_replication$age_group == "Young")])
old.rep.colors <- colorRampPalette(color.olds)(length(info.old.rep))
young.rep.colors <- colorRampPalette(color.youngs)(length(info.young.rep))
color.replications <- c(young.rep.colors, old.rep.colors)
names(color.replications) <- c(sort(info.young.rep), sort(info.old.rep))
plot.data <- PT.info.tissue %>% dplyr::select(-RNA.depth, -RNA.cells) %>% pivot_longer(cols = starts_with("H3K"), names_to = "histone", values_to = "depth") %>% filter(!is.na(depth))
plot.data$histone <- gsub(".depth", "", plot.data$histone)
plot.data$histone <- factor(plot.data$histone, levels = rev(names(color.histone)))
plot.data$histone <- droplevels(plot.data$histone)
plot.data <- plot.data %>% left_join(histone_age_sample_replication, by = c("histone", "age_group", "sample", "tissue", "batch", "tissue_new"))
plot.data.subset <- plot.data[which(plot.data$tissue == names(color.tissue)[1]),]
plot.data.subset$replication <- factor(plot.data.subset$replication, levels = names(color.replications))
total_depth <- plot.data.subset %>% group_by(histone) %>% summarise(total = sum(depth, na.rm = TRUE))
p3 <- ggbarplot(plot.data.subset, x = "histone", y = "depth", xlab = "", ylab = "Reads (M)", color = "replication", fill = "replication", palette = color.replications, label = F, orientation = "horiz", position = position_stack(reverse = TRUE)) + scale_y_continuous(limits = c(0, y.max)) + ggtitle(unique(plot.data.subset$tissue_new)) + theme(legend.position = "right", plot.title = element_text(size = 12, hjust = 0.5)) + geom_hline(yintercept = 50, lty = "dashed", lwd = 0.5, color = "black") + geom_text(data = total_depth, aes(x = histone, y = total, label = paste0(round(total, 2), "M")), inherit.aes = FALSE, hjust = -0.2, size = 3.5)
plot.list <- lapply(names(color.tissue)[-1], function(x){
  plot.data.subset <- plot.data[which(plot.data$tissue == x),]
  total_depth <- plot.data.subset %>% group_by(histone) %>% summarise(total = sum(depth, na.rm = TRUE))
  p <- ggbarplot(plot.data.subset, x = "histone", y = "depth", xlab = "", ylab = "Reads (M)", color = "replication", fill = "replication", palette = color.replications, label = F, orientation = "horiz", position = position_stack(reverse = TRUE)) + scale_y_continuous(limits = c(0, y.max)) + ggtitle(unique(plot.data.subset$tissue_new)) + theme(legend.position = "right", plot.title = element_text(size = 12, hjust = 0.5), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) + geom_hline(yintercept = 50, lty = "dashed", lwd = 0.5, color = "black") + geom_text(data = total_depth, aes(x = histone, y = total, label = paste0(round(total, 2), "M")), inherit.aes = FALSE, hjust = -0.2, size = 3.5)
  return(p)
})
aligned_plots2 <- align_plots(p3, plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]], plot.list[[5]], plot.list[[6]], plot.list[[7]], plot.list[[8]], plot.list[[9]], plot.list[[10]], plot.list[[11]], plot.list[[12]], align = "h",axis = "tb")
pdf("CellNumberAndVolume.tissue2.pdf", width = 45, height = 5)
grid.arrange(grobs = aligned_plots1, nrow = 1, ncol = 13, widths = c(4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3))
grid.arrange(grobs = aligned_plots2, nrow = 2, widths = c(4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3))
dev.off()

#### cell type
order.label <- tapply(PT.info$RNA.depth, paste0(PT.info$tissue, ":", PT.info$cell_type), sum)
order.label <- data.frame(label = names(order.label), depth = order.label)
order.label$tissue <- gsub(":.*", "", order.label$label)
order.label$cellType <- gsub(".*:", "", order.label$label)
order.label$tissue <- factor(order.label$tissue, levels = rev(names(color.tissue)))
order.label <- arrange(order.label, tissue, depth)
PT.info$cellType <- paste0(PT.info$tissue, ":", PT.info$cell_type)
PT.info$cellType <- factor(PT.info$cellType, levels = order.label$label)

# set sample palette
info <- unique(PT.info[, c("age_group", "sample")])
info.old <- info[which(info$age_group == "Old"),]
info.young <- info[which(info$age_group == "Young"),]
old.sample.colors <- colorRampPalette(color.olds)(nrow(info.old))
young.sample.colors <- colorRampPalette(color.youngs)(nrow(info.young))
color.samples <- c(young.sample.colors, old.sample.colors)
names(color.samples) <- c(info.young$sample, info.old$sample)
PT.info$sample <- factor(PT.info$sample, levels = names(color.samples))

# set y axis scale
DNA.volume <- PT.info[, grep("depth", colnames(PT.info))]
DNA.volume <- DNA.volume[,grep("RNA.depth", colnames(DNA.volume), invert = T)]
y.max <- rowsum(DNA.volume, group = PT.info$cellType, na.rm = TRUE)
y.max <- max(unlist(y.max))

library(cowplot)
library(gridExtra)
p1 <- ggbarplot(PT.info, x = "cellType", y = "RNA.cells", xlab = "", ylab = "Cell number", color = "sample", fill = "sample", palette = color.samples, label = F, orientation = "horiz", position = position_stack(reverse = TRUE)) + ggtitle("Cell number") + theme(legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5))
p2 <- ggbarplot(PT.info, x = "cellType", y = "RNA.depth", xlab = "", ylab = "Reads (M)", color = "sample", fill = "sample", palette = color.samples, label = F, orientation = "horiz", position = position_stack(reverse = TRUE)) + ggtitle("RNA") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5))
p3 <- ggbarplot(PT.info, x = "cellType", y = "H3K9me3.depth", xlab = "", ylab = "Reads (M)", color = "sample", fill = "sample", palette = color.samples, label = F, orientation = "horiz", position = position_stack(reverse = TRUE)) + ggtitle("H3K9me3") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5)) + geom_hline(yintercept = c(10), lty = "dashed", lwd = 0.5, color = "black")
p4 <- ggbarplot(PT.info, x = "cellType", y = "H3K27me3.depth", xlab = "", ylab = "Reads (M)", color = "sample", fill = "sample", palette = color.samples, label = F, orientation = "horiz", position = position_stack(reverse = TRUE)) + ggtitle("H3K27me3") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5)) + geom_hline(yintercept = c(10), lty = "dashed", lwd = 0.5, color = "black")
p5 <- ggbarplot(PT.info, x = "cellType", y = "H3K4me1.depth", xlab = "", ylab = "Reads (M)", color = "sample", fill = "sample", palette = color.samples, label = F, orientation = "horiz", position = position_stack(reverse = TRUE)) + ggtitle("H3K4me1") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5)) + geom_hline(yintercept = c(10), lty = "dashed", lwd = 0.5, color = "black")
p6 <- ggbarplot(PT.info, x = "cellType", y = "H3K27ac.depth", xlab = "", ylab = "Reads (M)", color = "sample", fill = "sample", palette = color.samples, label = F, orientation = "horiz", position = position_stack(reverse = TRUE)) + ggtitle("H3K27ac") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5)) + geom_hline(yintercept = c(10), lty = "dashed", lwd = 0.5, color = "black")
p7 <- ggbarplot(PT.info, x = "cellType", y = "H3K4me3.depth", xlab = "", ylab = "Reads (M)", color = "sample", fill = "sample", palette = color.samples, label = F, orientation = "horiz", position = position_stack(reverse = TRUE)) + ggtitle("H3K4me3") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5)) + geom_hline(yintercept = c(10), lty = "dashed", lwd = 0.5, color = "black")
p8 <- ggbarplot(PT.info, x = "cellType", y = "H3K36me3.depth", xlab = "", ylab = "Reads (M)", color = "sample", fill = "sample", palette = color.samples, label = F, orientation = "horiz", position = position_stack(reverse = TRUE))  + ggtitle("H3K36me3") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), legend.position = "none", plot.title = element_text(size = 12, hjust = 0.5)) + geom_hline(yintercept = c(10), lty = "dashed", lwd = 0.5, color = "black")
aligned_plots <- align_plots(p1, p2, p3, p4, p5, p6, p7, p8, align = "h",axis = "tb")
pdf("CellNumberAndVolume.cellType.pdf", width = 30, height = 40)
grid.arrange(grobs = aligned_plots, nrow = 1, widths = c(4, 2, 2, 2, 2, 2, 2, 2))
dev.off()

#### ---- General RNA/CPM matrix ---- ####
#### age level
count.mat.list <- lapply(names(tissue.seurat.final), function(tissue){
  seurat.obj <- tissue.seurat.final[[tissue]]
  seurat.obj$group <- paste0(seurat.obj$age_group, "__", tissue, "__", seurat.obj$cell_type)
  counts <- AggregateExpression(seurat.obj, group.by = "group", assays = "RNA", slot = "counts")$RNA
  counts <- as.data.frame(counts)
  counts$gene <- rownames(counts)
  return(counts)
})
count.mats <- Reduce(function(x,y) merge(x, y, all=T, by = "gene"), count.mat.list)
rownames(count.mats) <- count.mats$gene
count.mats <- count.mats[,-1]
saveRDS(count.mats, "RNA.count.mats.rds")

count.mats[is.na(count.mats)] <- 0
cpm.mats <- edgeR::cpm(count.mats)
saveRDS(cpm.mats, "RNA.cpm.mats.rds")

#### sample level
count.sample.mat.list <- lapply(names(tissue.seurat.final), function(tissue){
  seurat.obj <- tissue.seurat.final[[tissue]]
  seurat.obj$batch <- gsub("-.*", "", seurat.obj$assay)
  seurat.obj$group <- paste0(seurat.obj$age_group, "__", tissue, "__", seurat.obj$cell_type, "__", seurat.obj$sample, "__", seurat.obj$batch)
  counts <- AggregateExpression(seurat.obj, group.by = "group", assays = "RNA", slot = "counts")$RNA
  counts <- as.data.frame(counts)
  counts$gene <- rownames(counts)
  return(counts)
})
count.sample.mat.list <- Reduce(function(x,y) merge(x, y, all=T, by = "gene"), count.sample.mat.list)
rownames(count.sample.mat.list) <- count.sample.mat.list$gene
count.sample.mat.list <- count.sample.mat.list[,-1]
saveRDS(count.sample.mat.list, "RNA.count.sample.mats.rds")

count.sample.mat.list[is.na(count.sample.mat.list)] <- 0
cpm.sample.mat.list <- edgeR::cpm(count.sample.mat.list)
saveRDS(cpm.sample.mat.list, "RNA.cpm.sample.mats.rds")

# median cpm for each age group
info <- data.frame(label = colnames(cpm.sample.mat.list))
info <- info %>% separate_wider_delim(label, delim = "--", names = c("age_group", "tissue", "celltype", "sample", "batch"))
info$label <- colnames(cpm.sample.mat.list)
saveRDS(info, "RNA.sample.mats.rds")

info$group <- paste0(info$tissue, ":", info$celltype, "_", info$age_group)
medianCPM.mats <- apply(cpm.sample.mat.list, 1, function(x){
  medianCPM <- tapply(x, info$group, median, na.rm = TRUE)
  return(medianCPM)
})
medianCPM.mats <- t(medianCPM.mats)
colnames(medianCPM.mats) <- gsub("-", "_", colnames(medianCPM.mats))
saveRDS(medianCPM.mats, "RNA.medianCPM.mats.rds")
