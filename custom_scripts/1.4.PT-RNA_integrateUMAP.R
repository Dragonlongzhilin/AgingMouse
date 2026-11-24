#' @description integrate all tissue to plot UMAP
#' @Author: Zhilin Long
#' @date 2026-6-10

# ulimit -n 200000 # set the maximum number of open files

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
PT.info <- readRDS("PT.info.rds")
colors.info <- unique(PT.info[,c("cellType", "colors")])
color.cellType <- colors.info$colors
names(color.cellType) <- colors.info$cellType
color.cellType <- color.cellType[levels(PT.info$cellType)]

tissue.seurat.final <- readRDS("tissue.seurat.final.rds")
dir.create("./MergeTissue/BP_matrix", recursive = TRUE)
#### BPCells merge all cells
seurat.BP.list <- lapply(names(tissue.seurat.final), function(tissue){
  cat("Analysis ", tissue, "\n")
  seurat_obj <- tissue.seurat.final[[tissue]]
  metadata <- seurat_obj@meta.data
  data_dir <- paste0("./MergeTissue/BP_matrix/", tissue)
  write_matrix_dir(mat = seurat_obj[["RNA"]]$counts, dir = data_dir, overwrite = T)
  counts_mat <- open_matrix_dir(dir = data_dir)
  seurat_obj <- CreateSeuratObject(counts = counts_mat, meta.data = metadata)
  return(seurat_obj)
})
names(seurat.BP.list) <- names(tissue.seurat.final)
seurat.BP.merge <- merge(x = seurat.BP.list[[1]], y = seurat.BP.list[-1], add.cell.ids = names(seurat.BP.list))
seurat.BP.merge[["RNA"]] <- JoinLayers(seurat.BP.merge[["RNA"]])
write_matrix_dir(mat = seurat.BP.merge[["RNA"]]$counts, dir = './MergeTissue/BP_matrix/AllTissue', overwrite = T)
metadata <- seurat.BP.merge@meta.data
saveRDS(metadata, "./MergeTissue/BP_matrix/metadata.rds")

#### ---- Seurat integration ---- ####
#### subset 5M cells for each tissue
tissue.seurat.subset <- lapply(tissue.seurat.final, function(x){
  x <- DietSeurat(x, layers = "RNA")
  x <- NormalizeData(x)
  x[["RNA"]] <- split(x[["RNA"]], f = x$histone_sample_batch)
  x <- FindVariableFeatures(x, verbose = FALSE)
  seuratObj <- SketchData(object = x, ncells = 2000, method = "LeverageScore", sketched.assay = "sketch")
  seuratObj.subset <- subset(x, cells = colnames(seuratObj[["sketch"]]))
  return(seuratObj.subset)
})
saveRDS(tissue.seurat.subset, "MergeTissue/tissue.seurat.final.subset.rds")

#### BPCells merge cellsubset
seurat.BP.subset <- lapply(names(tissue.seurat.subset), function(tissue){
  cat("Analysis ", tissue, "\n")
  seurat_obj <- tissue.seurat.subset[[tissue]]
  metadata <- seurat_obj@meta.data
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
  data_dir <- paste0("./MergeTissue/BP_subset_matrix/", tissue)
  write_matrix_dir(mat = seurat_obj[["RNA"]]$counts, dir = data_dir, overwrite = T)
  counts_mat <- open_matrix_dir(dir = data_dir)
  seurat_obj <- CreateSeuratObject(counts = counts_mat, meta.data = metadata)
  return(seurat_obj)
})
names(seurat.BP.subset) <- names(tissue.seurat.subset)

#### batch correction & cluster
tissue.seurat.sketch <- merge(x = seurat.BP.subset[[1]], y = seurat.BP.subset[-1], add.cell.ids = names(seurat.BP.subset))
tissue.seurat.sketch[["RNA"]] <- JoinLayers(tissue.seurat.sketch[["RNA"]])
DefaultAssay(tissue.seurat.sketch) <- "RNA"
tissue.seurat.sketch[["RNA"]] <- split(tissue.seurat.sketch[["RNA"]], f = tissue.seurat.sketch$histone_sample_batch)
tissue.seurat.sketch <- NormalizeData(tissue.seurat.sketch)
tissue.seurat.sketch <- FindVariableFeatures(tissue.seurat.sketch, verbose = FALSE)
tissue.seurat.sketch <- ScaleData(tissue.seurat.sketch, vars.to.regress = c("mt_ratio_RNA", "nCount_RNA"), verbose = F)
tissue.seurat.sketch <- RunPCA(tissue.seurat.sketch, verbose = F)
tissue.seurat.sketch <- IntegrateLayers(object = tissue.seurat.sketch, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
tissue.seurat.sketch <- FindNeighbors(tissue.seurat.sketch, reduction = "harmony", dims = 1:50, verbose = FALSE)
tissue.seurat.sketch <- FindClusters(tissue.seurat.sketch, resolution = set.resolutions, verbose = FALSE)
tissue.seurat.sketch <- RunUMAP(tissue.seurat.sketch, dims = 1:50, n.neighbors = 50, min.dist = 0.15, reduction = "harmony", verbose = FALSE, return.model = T)
saveRDS(tissue.seurat.sketch, file = "MergeTissue/tissue.seurat.sketch.rds")

tissue.seurat.sketch$tissue_cellType <- paste0(tissue.seurat.sketch$tissue, ":", tissue.seurat.sketch$cell_type)
info.subset <- unique(PT.info[,c("tissue_cellType", "category")])
category.info <- plyr::mapvalues(tissue.seurat.sketch$tissue_cellType, from = info.subset$tissue_cellType, to = as.character(info.subset$category))
tissue.seurat.sketch$category <- as.character(category.info)
tissue.seurat.sketch$seurat_clusters <- tissue.seurat.sketch$RNA_snn_res.0.4

source("~/Functions/refAnnotationRatio.R")
assign.ratio <- refAnnotationRatio(tissue.seurat.sketch@meta.data, predictedLabels = "category", width = unit(6, "cm"), height = unit(10, "cm"), plot = T)
assign.ratio$cluster <- as.character(assign.ratio$cluster)
assign.ratio$category.annotation <- as.character(assign.ratio$category.annotation)
category.annotation <- plyr::mapvalues(tissue.seurat.sketch$seurat_clusters, from = as.character(assign.ratio$cluster), to = as.character(assign.ratio$category.annotation))
tissue.seurat.sketch <- AddMetaData(tissue.seurat.sketch, category.annotation, "category.annotation")
tissue.seurat.sketch$category.annotation <- factor(tissue.seurat.sketch$category.annotation, levels = names(color.category))
tissue.seurat.sketch$cellType <- factor(tissue.seurat.sketch$cellType, levels = levels(PT.info$cellType))

pdf("MergeTissue/tissue.seurat.sketch.pdf", width = 30, height = 8)
DimPlot(tissue.seurat.sketch, reduction = "umap", group.by = "seurat_clusters", cols = color.sample, label = T, repel = T, label.size = 4, raster=FALSE) + theme(aspect.ratio = 1)
DimPlot(tissue.seurat.sketch, reduction = "umap", group.by = "cellType", cols = color.cellType, label = T, repel = T, label.size = 3, raster=FALSE) + ggtitle(paste0("Merge tissue: ", ncol(tissue.seurat.sketch), " cells")) + theme(aspect.ratio = 1, plot.title = element_text(size = 16, hjust = 0.5))
DimPlot(tissue.seurat.sketch, reduction = "umap", group.by = "tissue", cols = color.tissue, label.size = 3, raster=FALSE) + ggtitle(paste0("Merge tissue: ", ncol(tissue.seurat.sketch), " cells")) + theme(aspect.ratio = 1, plot.title = element_text(size = 16, hjust = 0.5))
DimPlot(tissue.seurat.sketch, reduction = "umap", group.by = "category", cols = color.category, label.size = 3, raster=FALSE) + ggtitle(paste0("Merge tissue: ", ncol(tissue.seurat.sketch), " cells")) + theme(aspect.ratio = 1, plot.title = element_text(size = 16, hjust = 0.5))
DimPlot(tissue.seurat.sketch, reduction = "umap", group.by = "category.annotation", cols = color.category, label.size = 3, raster=FALSE) + ggtitle(paste0("Merge tissue: ", ncol(tissue.seurat.sketch), " cells")) + theme(aspect.ratio = 1, plot.title = element_text(size = 16, hjust = 0.5))
dev.off()

#### Extend sketch results to the full datasets
metadata_all <- readRDS("./MergeTissue/BP_matrix/metadata.rds")
counts_all <- open_matrix_dir(dir = "./MergeTissue/BP_matrix/AllTissue")
tissue.seurat <- CreateSeuratObject(counts = counts_all, meta.data = metadata_all)
DefaultAssay(tissue.seurat) <- "RNA"
tissue.seurat <- NormalizeData(tissue.seurat)
tissue.seurat[["RNA"]] <- split(tissue.seurat[["RNA"]], f = tissue.seurat$histone_sample_batch)

tissue.seurat[["sketch"]] <- tissue.seurat.sketch[["RNA"]]
DefaultAssay(tissue.seurat) <- "sketch"
tissue.seurat[["pca"]] <- tissue.seurat.sketch[["pca"]]
tissue.seurat@reductions$pca@assay.used <- "sketch"
tissue.seurat[["harmony"]] <- tissue.seurat.sketch[["harmony"]]
tissue.seurat@reductions$harmony@assay.used <- "sketch"
tissue.seurat[["umap"]] <- tissue.seurat.sketch[["umap"]]
tissue.seurat@reductions$umap@assay.used <- "sketch"
# tissue.seurat[["sketch"]] <- split(tissue.seurat[["sketch"]], f = tissue.seurat$histone_sample_batch)
tissue.seurat <- ProjectIntegration(object = tissue.seurat, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
tissue.seurat <- ProjectData(object = tissue.seurat, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "harmony", full.reduction = "harmony.full", dims = 1:50, refdata = list(cellType.full = "cellType"))
tissue.seurat <- RunUMAP(tissue.seurat, reduction = "harmony.full", dims = 1:50, reduction.name = "umap.full", reduction.key = "UMAP_full_", n.neighbors = 40, min.dist = 0.15)
tissue.seurat[["RNA"]] <- JoinLayers(tissue.seurat[["RNA"]])
saveRDS(tissue.seurat, file = "MergeTissue/tissue.seurat.rds")

DefaultAssay(tissue.seurat) <- "RNA"
tissue.seurat$tissue_cellType <- paste0(tissue.seurat$tissue, ":", tissue.seurat$cell_type)
info.subset <- unique(PT.info[,c("tissue_cellType", "category")])
category.info <- plyr::mapvalues(tissue.seurat$tissue_cellType, from = info.subset$tissue_cellType, to = as.character(info.subset$category))
tissue.seurat$category <- as.character(category.info)
tissue.seurat$tissue_new <- recode(as.character(tissue.seurat$tissue), "brainCB" = "Cerebellum", "brainFC" = "Frontal cortex", "brainHip" = "Hippocampus", "heart" = "Heart","liver" = "Liver", "lung" = "Lung", "kidney" = "Kidney", "BAT" = "BAT", "muscle" = "Muscle", "spleen" = "Spleen", "colon" = "Colon", "skin" = "Skin", "testis" = "Testis")
tissue.seurat$tissue_new <- factor(tissue.seurat$tissue_new, levels = names(color.tissue2))
tissue.seurat$category <- factor(tissue.seurat$category, levels = names(color.category))
tissue.seurat$histone <- factor(tissue.seurat$histone, levels = names(color.histone))
tissue.seurat$age_group <- factor(tissue.seurat$age_group, levels = names(color.age_group))
tissue.seurat$cellType <- factor(tissue.seurat$cellType, levels = levels(PT.info$cellType))
pdf("MergeTissue/tissue.seurat.pdf", width = 30, height = 8)
DimPlot(tissue.seurat, reduction = "umap.full", group.by = "cellType", cols = color.cellType, label = T, repel = T, label.size = 3) + ggtitle(paste0("Merge tissue: ", ncol(tissue.seurat), " cells")) + theme(aspect.ratio = 1, plot.title = element_text(size = 16, hjust = 0.5))
DimPlot(tissue.seurat, reduction = "umap.full", group.by = "tissue", cols = color.tissue, label.size = 3) + ggtitle(paste0("Merge tissue: ", ncol(tissue.seurat), " cells")) + theme(aspect.ratio = 1, plot.title = element_text(size = 16, hjust = 0.5))
DimPlot(tissue.seurat, reduction = "umap.full", group.by = "category", cols = color.category, label.size = 3) + ggtitle(paste0("Merge tissue: ", ncol(tissue.seurat), " cells")) + theme(aspect.ratio = 1, plot.title = element_text(size = 16, hjust = 0.5))
dev.off()
png("MergeTissue/allcellType.umap.png", width = 30, height = 8, units = "in", res = 300)
DimPlot(tissue.seurat, reduction = "umap.full", group.by = "cellType", cols = color.cellType) + xlab("UMAP_1") + ylab("UMAP_2") + ggtitle("") + theme(aspect.ratio = 1, plot.title = element_text(size = 16, hjust = 0.5, face = "plain")) + guides(color = guide_legend(title = "Cell type", ncol = 3, override.aes = list(size = 2.5))) + scale_color_manual(values = color.cellType[levels(PT.info$cellType)], labels = levels(PT.info$cellType_new))
dev.off()
pdf("MergeTissue/allcellType.umap.pdf", width = 30, height = 8)
DimPlot(tissue.seurat, reduction = "umap.full", group.by = "cellType", cols = color.cellType) + xlab("UMAP_1") + ylab("UMAP_2") + ggtitle("") + theme(aspect.ratio = 1, plot.title = element_text(size = 16, hjust = 0.5, face = "plain")) + guides(color = guide_legend(title = "Cell type", ncol = 3, override.aes = list(size = 2.5))) + scale_color_manual(values = color.cellType[levels(PT.info$cellType)], labels = levels(PT.info$cellType_new))
dev.off()
png("MergeTissue/allTissue.umap.png", width = 10, height = 8, units = "in", res = 300)
DimPlot(tissue.seurat, reduction = "umap.full", group.by = "tissue_new", cols = color.tissue2) + xlab("UMAP_1") + ylab("UMAP_2") + ggtitle("") + theme(aspect.ratio = 1, plot.title = element_text(size = 16, hjust = 0.5)) + guides(color = guide_legend(title = "Tissue", override.aes = list(size = 4)))
dev.off()
png("MergeTissue/allCategory.umap.png", width = 10, height = 8, units = "in", res = 300)
DimPlot(tissue.seurat, reduction = "umap.full", group.by = "category", cols = color.category) + xlab("UMAP_1") + ylab("UMAP_2") + ggtitle("") + theme(aspect.ratio = 1, plot.title = element_text(size = 16, hjust = 0.5, face = "plain")) + guides(color = guide_legend(title = "Category", override.aes = list(size = 4))) # text = element_text(size = 20), axis.text = element_text(size = 15), axis.title = element_text(size = 16)
dev.off()
png("MergeTissue/allHistone.umap.png", width = 10, height = 8, units = "in", res = 300)
DimPlot(tissue.seurat, reduction = "umap.full", group.by = "histone", cols = color.histone) + ggtitle("") + xlab("UMAP_1") + ylab("UMAP_2") + theme(aspect.ratio = 1, plot.title = element_text(size = 16, hjust = 0.5)) + guides(color = guide_legend(title = "Histone", override.aes = list(size = 4)))
dev.off()
png("MergeTissue/allAge.umap.png", width = 10, height = 8, units = "in", res = 300)
DimPlot(tissue.seurat, reduction = "umap.full", group.by = "age_group", cols = color.age_group) + ggtitle("") + xlab("UMAP_1") + ylab("UMAP_2") + theme(aspect.ratio = 1, plot.title = element_text(size = 16, hjust = 0.5)) + guides(color = guide_legend(title = "Age", override.aes = list(size = 4)))
dev.off()

RNA_age_sample_replication <- readRDS("RNA_age_sample_replication.rds")
info.Old.rep <- gtools::mixedsort(as.character(unique(RNA_age_sample_replication$replication[which(RNA_age_sample_replication$age_group == "Old")])))
info.Young.rep <- gtools::mixedsort(as.character(unique(RNA_age_sample_replication$replication[which(RNA_age_sample_replication$age_group == "Young")])))
Old.rep.colors <- colorRampPalette(color.olds)(length(info.Old.rep))
Young.rep.colors <- colorRampPalette(color.youngs)(length(info.Young.rep))
color.replications <- c(Young.rep.colors, Old.rep.colors)
names(color.replications) <- c(info.Young.rep, info.Old.rep)
meta.data <- tissue.seurat@meta.data %>% left_join(RNA_age_sample_replication[,c("tissue", "age_group", "sample", "batch", "replication")], by = c("tissue", "age_group", "sample", "batch"))
tissue.seurat$replication <- factor(meta.data$replication.y, levels = names(color.replications))
png("MergeTissue/allReplication.umap.png", width = 10, height = 8, units = "in", res = 300)
DimPlot(tissue.seurat, reduction = "umap.full", group.by = "replication", cols = color.replications, label.size = 3) + ggtitle("") + xlab("UMAP_1") + ylab("UMAP_2") + theme(aspect.ratio = 1, plot.title = element_text(size = 16, hjust = 0.5)) + guides(color = guide_legend(title = "Sample", override.aes = list(size = 4)))
dev.off()

#### cell marker expression
tissue.seurat <- readRDS("MergeTissue/tissue.seurat.rds")

PT.info.cellTypeMarker <- readxl::read_xlsx("PT.info.cellTypeMarker.xlsx")
PT.info.cellTypeMarker$category <- factor(PT.info.cellTypeMarker$category, levels = names(color.category))
cellTypeMarker.list <- split(PT.info.cellTypeMarker$GeneMarkers, PT.info.cellTypeMarker$category)
cellTypeMarker.list <- lapply(cellTypeMarker.list, function(x){return(unique(unlist(strsplit(x, ", "))))})
cellTypeMarker.df <- stack(cellTypeMarker.list)
colnames(cellTypeMarker.df) <- c("gene", "category")
cellTypeMarker.df <- cellTypeMarker.df[, c("category", "gene")]
cellTypeMarker.df <- cellTypeMarker.df[!duplicated(cellTypeMarker.df$gene),]

features.use <- intersect(unique(as.character(cellTypeMarker.df$gene)), rownames(tissue.seurat[["RNA"]]))
# tissue_cellType.exp <- AverageExpression(tissue.seurat, features = features.use, group.by = "tissue_cellType", layer = "data")
# plot.data <- tissue_cellType.exp$RNA

# Manual pseudobulk: bypass buggy Seurat AverageExpression internal comparison
data.mat <- GetAssayData(tissue.seurat, assay = "RNA", layer = "data")
data.mat <- data.mat[features.use, , drop = FALSE]
groups <- as.character(tissue.seurat$tissue_cellType)
group_levels <- unique(groups)
plot.data <- do.call(cbind, lapply(group_levels, function(g) {
  cells <- colnames(tissue.seurat)[groups == g]
  if (length(cells) > 1) {
    Matrix::rowMeans(data.mat[, cells, drop = FALSE])
  } else {
    data.mat[, cells]
  }
}))
colnames(plot.data) <- group_levels

# Match gene-to-category annotation
gene_anno <- cellTypeMarker.df[match(features.use, cellTypeMarker.df$gene), ]
gene_anno$category <- factor(gene_anno$category, levels = names(color.category))
gene_anno <- gene_anno[order(gene_anno$category), ]
plot.data <- plot.data[gene_anno$gene, ]

# Z-score scaling across tissue_cellType
plot.data.scaled <- t(scale(t(plot.data)))

# Row annotation: category colors
row_ha <- rowAnnotation(
  Category = gene_anno$category,
  col = list(Category = color.category),
  show_annotation_name = TRUE,
  annotation_legend_param = list(
    Category = list(title = "Category", title_position = "topcenter")
  )
)
row_split <- gene_anno$category

# Column annotation: tissue colors based on tissue_cellType
col_anno <- data.frame(tissue_cellType = colnames(plot.data.scaled))
col_anno <- col_anno %>% left_join(unique(PT.info[,c("tissue_cellType", "tissue_cellType_new", "tissue_new", "cellType", "category")]), by = "tissue_cellType")
col_anno$category <- factor(col_anno$category, levels = names(color.category))
col_ha <- HeatmapAnnotation(
  Category = col_anno$category,
  Tissue = col_anno$tissue_new,
  col = list(Tissue = color.tissue2, Category = color.category))
col_split <- col_anno$category
colnames(plot.data.scaled) <- col_anno$tissue_cellType_new
colnames(plot.data) <- col_anno$tissue_cellType_new

# Heatmap color palette
heatmap_colors <- colorRamp2(c(quantile(plot.data.scaled, 0.05), 0, quantile(plot.data.scaled, 0.95)), c("#2166AC", "white", "#B2182B"))

# Draw heatmap with row_split by category
ht <- Heatmap(
  plot.data.scaled,
  name = "Z-score",
  col = heatmap_colors,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 8),
  row_split = row_split,
  column_split = col_split,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 11),
  top_annotation = col_ha,
  right_annotation = row_ha,
  row_title = "%s",
  use_raster = FALSE,
  border = TRUE
)

pdf("cellTypeMarker.heatmap.pdf", width = 25, height = 50)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

#### plot heatmap for each tissue
tissue.seurat <- readRDS("MergeTissue/tissue.seurat.rds")
PT.info.cellTypeMarker <- readxl::read_xlsx("PT.info.cellTypeMarker.xlsx")
PT.info.cellTypeMarker$category <- factor(PT.info.cellTypeMarker$category, levels = names(color.category))

pdf("cellTypeMarker.tissue.heatmap.pdf", width = 20, height = 20)
res <- lapply(names(color.tissue2), function(x){
  info.subset <- PT.info.cellTypeMarker[PT.info.cellTypeMarker$tissue == x, ]
  info.subset <- arrange(info.subset, category)
  cellType_order <- unique(info.subset$tissue_cellType)
  genes <- unique(unlist(strsplit(info.subset$GeneMarkers, ", ")))
  idx <- match(genes, rownames(plot.data))
  plot.data.subset <- plot.data[na.omit(idx), cellType_order]
  colnames(plot.data.subset) <- gsub(".*:", "", colnames(plot.data.subset))

  heatmap_colors <- colorRamp2(c(0, quantile(plot.data.subset, 0.95)), c("white", "#B2182B"))
  ht <- Heatmap(as.matrix(plot.data.subset), name = "Normalized expression", col = heatmap_colors, cluster_rows = F, cluster_columns = F, show_row_names = TRUE, show_column_names = TRUE, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), row_title_rot = 0, row_title_gp = gpar(fontsize = 11), use_raster = FALSE, border = TRUE, column_title = x, width = unit(0.4 * ncol(plot.data.subset), "cm"), height = unit(0.4 * nrow(plot.data.subset), "cm"))
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")

  plot.data.subset.scaled <- t(scale(t(plot.data.subset)))
  heatmap_colors <- colorRamp2(c(min(plot.data.subset.scaled), 0, max(plot.data.subset.scaled)), c("#2166AC", "white", "#B2182B"))
  ht <- Heatmap(as.matrix(plot.data.subset.scaled), name = "Z-score", col = heatmap_colors, cluster_rows = F, cluster_columns = F, show_row_names = TRUE, show_column_names = TRUE, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), row_title_rot = 0, row_title_gp = gpar(fontsize = 11), use_raster = FALSE, border = TRUE, column_title = x, width = unit(0.4 * ncol(plot.data.subset.scaled), "cm"), height = unit(0.4 * nrow(plot.data.subset.scaled), "cm"))
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
})
dev.off()