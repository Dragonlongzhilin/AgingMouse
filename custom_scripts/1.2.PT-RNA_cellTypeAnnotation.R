#' @description validation cell type annotation with reference dataset
#' @Author: Zhilin Long
#' @date 2026-7-28

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
source("~/Functions/refAnnotationRatio.R")
setwd("~/AgingMouse/PairedTag/1.dataPreprocessing")
tissue.seurat.final <- readRDS("tissue.seurat.final.rds")
PT.info <- readRDS("PT.info.rds")
colors.info <- unique(PT.info[,c("cellType", "colors")])
color.cellType <- colors.info$colors
names(color.cellType) <- colors.info$cellType

#### ---- brainCB ---- ####
## paper: A transcriptomic atlas of mouse cerebellar cortex comprehensively defines cell types | Nature 2021
Kozareva.2021.Nature <- readRDS("~/AgingMouse/singlecell/brain/Kozareva.2021.Nature/brainCB_BP.rds")
DefaultAssay(Kozareva.2021.Nature) <- "sketch"
Yao.2023.Nature.CB <- readRDS("~/AgingMouse/singlecell/brain/Yao.2023.Nature/CB.rds") # cerebellum
DefaultAssay(Yao.2023.Nature.CB) <- "RNA"

PairedTag_data <- tissue.seurat.final$brainCB
anchors <- FindTransferAnchors(
  reference = Yao.2023.Nature.CB,
  query = PairedTag_data,
  reference.reduction = "pca",
  dims = 1:30
)
predictions <- TransferData(anchorset = anchors, refdata = Yao.2023.Nature.CB$subclass, dims = 1:30)
PairedTag_data <- AddMetaData(PairedTag_data, metadata = predictions$predicted.id, col.name = "Yao.2023.Nature.CB")

anchors <- FindTransferAnchors(
  reference = Kozareva.2021.Nature,
  query = PairedTag_data,
  reference.reduction = "pca",
  dims = 1:30
)
ref_labels <- as.character(Kozareva.2021.Nature$cluster[colnames(Kozareva.2021.Nature@assays$sketch)])
predictions <- TransferData(anchorset = anchors, refdata = ref_labels, dims = 1:30)
PairedTag_data <- AddMetaData(PairedTag_data, metadata = predictions$predicted.id, col.name = "Kozareva.2021.Nature")
tissue.seurat.final$brainCB <- PairedTag_data

seuratObj <- tissue.seurat.final$brainCB
pdf("brainCB/referenceMap.pdf")
assign.ratio <- refAnnotationRatio(seuratObj@meta.data, predictedLabels = "Yao.2023.Nature.CB", width = unit(8, "cm"), height = unit(7, "cm"), plot = T, group.by = "cellType", row_title = "Our dataset:Cerebellum")
assign.ratio <- refAnnotationRatio(seuratObj@meta.data, predictedLabels = "Kozareva.2021.Nature", width = unit(6, "cm"), height = unit(6, "cm"), plot = T, group.by = "cellType", row_title = "Our dataset:Cerebellum")
dev.off()

#### ---- brainFC ---- ####
## paper: A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation
Yao.2021.Cell <- readRDS("~/AgingMouse/singlecell/brain/Yao.2021.Cell/10X/seuratObj_10X_BP.rds")
DefaultAssay(Yao.2021.Cell) <- "sketch"
Yao.2021.Cell.FC <- readRDS("~/AgingMouse/singlecell/brain/Yao.2021.Cell/10X/seuratObj_10X_BP.FC.rds")
DefaultAssay(Yao.2021.Cell.FC) <- "RNA"
## paper: Continuous cell-type diversification in mouse visual cortex development | Nature 2025
Gao.2025.Nature <- readRDS("~/AgingMouse/singlecell/brain/Gao.2025.Nature/multiome_snRNA.rds")
DefaultAssay(Gao.2025.Nature) <- "RNA"
## paper: A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain | Nature 2023
Yao.2023.Nature.Isocortex <- readRDS("~/AgingMouse/singlecell/brain/Yao.2023.Nature/Isocortex.rds") # Isocortex
DefaultAssay(Yao.2023.Nature.Isocortex) <- "RNA"

PairedTag_data <- tissue.seurat.final$brainFC
anchors <- FindTransferAnchors(
  reference = Yao.2021.Cell.FC,
  query = PairedTag_data,
  reference.reduction = "pca",
  dims = 1:30
)
predictions <- TransferData(anchorset = anchors, refdata = Yao.2021.Cell.FC$subclass_label, dims = 1:30)
PairedTag_data <- AddMetaData(PairedTag_data, metadata = predictions$predicted.id, col.name = "Yao.2021.Cell.FC")

anchors <- FindTransferAnchors(
  reference = Gao.2025.Nature,
  query = PairedTag_data,
  reference.reduction = "pca",
  dims = 1:30
)
predictions <- TransferData(anchorset = anchors, refdata = Gao.2025.Nature$subclass_label, dims = 1:30)
PairedTag_data <- AddMetaData(PairedTag_data, metadata = predictions$predicted.id, col.name = "Gao.2025.Nature")

anchors <- FindTransferAnchors(
  reference = Yao.2023.Nature.Isocortex,
  query = PairedTag_data,
  reference.reduction = "pca",
  dims = 1:30
)
predictions <- TransferData(anchorset = anchors, refdata = Yao.2023.Nature.Isocortex$subclass, dims = 1:30)
PairedTag_data <- AddMetaData(PairedTag_data, metadata = predictions$predicted.id, col.name = "Yao.2023.Nature.Isocortex")
tissue.seurat.final$brainFC <- PairedTag_data

seuratObj <- tissue.seurat.final$brainFC
pdf("brainFC/referenceMap.pdf", width = 20, height = 12)
assign.ratio <- refAnnotationRatio(seuratObj@meta.data, predictedLabels = "Yao.2021.Cell.FC", width = unit(9, "cm"), height = unit(8, "cm"), plot = T, group.by = "cellType", row_title = "Our dataset:Frontal cortex")
assign.ratio <- refAnnotationRatio(seuratObj@meta.data, predictedLabels = "Yao.2023.Nature.Isocortex", width = unit(25, "cm"), height = unit(8, "cm"), plot = T, group.by = "cellType", row_title = "Our dataset:Frontal cortex")
assign.ratio <- refAnnotationRatio(seuratObj@meta.data, predictedLabels = "Gao.2025.Nature", width = unit(12, "cm"), height = unit(8, "cm"), plot = T, group.by = "cellType", row_title = "Our dataset:Frontal cortex")
dev.off()

#### ---- brainHip ---- ####
Yao.2021.Cell.Hip <- readRDS("~/AgingMouse/singlecell/brain/Yao.2021.Cell/10X/seuratObj_10X_BP.Hip.rds")
DefaultAssay(Yao.2021.Cell.Hip) <- "RNA"
Yao.2023.Nature.HPF <- readRDS("~/AgingMouse/singlecell/brain/Yao.2023.Nature/HPF.rds") # hippocampal formation (HPF)
DefaultAssay(Yao.2023.Nature.HPF) <- "RNA"

PairedTag_data <- tissue.seurat.final$brainHip
anchors <- FindTransferAnchors(
  reference = Yao.2021.Cell.Hip,
  query = PairedTag_data,
  reference.reduction = "pca",
  dims = 1:30
)
predictions <- TransferData(anchorset = anchors, refdata = Yao.2021.Cell.Hip$subclass_label, dims = 1:30)
PairedTag_data <- AddMetaData(PairedTag_data, metadata = predictions$predicted.id, col.name = "Yao.2021.Cell.Hip")

anchors <- FindTransferAnchors(
  reference = Yao.2023.Nature.HPF,
  query = PairedTag_data,
  reference.reduction = "pca",
  dims = 1:30
)
predictions <- TransferData(anchorset = anchors, refdata = Yao.2023.Nature.HPF$subclass, dims = 1:30)
PairedTag_data <- AddMetaData(PairedTag_data, metadata = predictions$predicted.id, col.name = "Yao.2023.Nature.HPF")
tissue.seurat.final$brainHip <- PairedTag_data

seuratObj <- tissue.seurat.final$brainHip
pdf("brainHip/referenceMap.pdf", width = 12, height = 10)
assign.ratio <- refAnnotationRatio(seuratObj@meta.data, predictedLabels = "Yao.2021.Cell.Hip", width = unit(12, "cm"), height = unit(8, "cm"), plot = T, group.by = "cellType", row_title = "Our dataset:Hippocampus")
assign.ratio <- refAnnotationRatio(seuratObj@meta.data, predictedLabels = "Yao.2023.Nature.HPF", width = unit(20, "cm"), height = unit(8, "cm"), plot = T, group.by = "cellType", row_title = "Our dataset:Hippocampus")
dev.off()

#### ---- liver/heart/lung/colon/kidney/BAT/muscle ---- ####
tissues <- c("colon", "muscle", "BAT", "liver", "lung", "heart", "kidney")
for(tissue in tissues){
	cat("Analysis ", tissue, "\n")
	PairedTag_data <- tissue.seurat.final[[tissue]]
	## Zhang.2024.Science.EasySci
	Zhang.2024.Science.EasySci <- readRDS(paste0("~/AgingMouse/singlecell/Altas/Zhang.2024.Science.EasySci/", tissue, "/", tissue, ".rds"))
	anchors <- FindTransferAnchors(
		reference = Zhang.2024.Science.EasySci,
		query = PairedTag_data,
		reference.reduction = "pca",
		dims = 1:30
	)
	predictions <- TransferData(anchorset = anchors, refdata = Zhang.2024.Science.EasySci$Main_cell_type, dims = 1:30)
	PairedTag_data <- AddMetaData(PairedTag_data, metadata = predictions$predicted.id, col.name = "Zhang.2024.Science.EasySci")
	tissue.seurat.final[[tissue]] <- PairedTag_data
}
saveRDS(tissue.seurat.final, "tissue.seurat.final.rds")

seuratObj <- tissue.seurat.final$colon
pdf("colon/referenceMap.pdf")
assign.ratio <- refAnnotationRatio(seuratObj@meta.data, predictedLabels = "Zhang.2024.Science.EasySci", width = unit(7, "cm"), height = unit(6, "cm"), plot = T, group.by = "cellType", row_title = "Our dataset:Colon", Unpredictable_cellType = "Vascular endothelial cells")
dev.off()

#### ---- kidney ---- ####
PairedTag_data <- tissue.seurat.final$kidney
Chen.2025.NatGenet <- readRDS("~/AgingMouse/singlecell/kidney/Chen.2025.NatGenet/multiome_10x.rds")
anchors <- FindTransferAnchors(
  reference = Chen.2025.NatGenet,
  query = PairedTag_data,
  reference.reduction = "pca",
  dims = 1:30
)
predictions <- TransferData(anchorset = anchors, refdata = Chen.2025.NatGenet$celltype_final, dims = 1:30)
PairedTag_data <- AddMetaData(PairedTag_data, metadata = predictions$predicted.id, col.name = "Chen.2025.NatGenet")
tissue.seurat.final$kidney <- PairedTag_data

seuratObj <- tissue.seurat.final$kidney
pdf("kidney/referenceMap2.pdf")
assign.ratio <- refAnnotationRatio(seuratObj@meta.data, predictedLabels = "Chen.2025.NatGenet", width = unit(9, "cm"), height = unit(6, "cm"), plot = T, group.by = "cellType", row_title = "Our dataset:Kidney", Unpredictable_cellType = c("Adipocytes", "Juxtaglomerular cells", "Ncam1 positive cells"))
dev.off()

#### ---- spleen ---- ####
PairedTag_data <- tissue.seurat.final$spleen
Zhang.2024.Science.EasySci.Immune <- readRDS("~/AgingMouse/singlecell/Altas/Zhang.2024.Science.EasySci/Immune/Immune.subset.rds")
anchors <- FindTransferAnchors(
  reference = Zhang.2024.Science.EasySci.Immune,
  query = PairedTag_data,
  reference.reduction = "pca",
  dims = 1:30
)
predictions <- TransferData(anchorset = anchors, refdata = Zhang.2024.Science.EasySci.Immune$Main_cell_type, dims = 1:30)
PairedTag_data <- AddMetaData(PairedTag_data, metadata = predictions$predicted.id, col.name = "Zhang.2024.Science.EasySci")
tissue.seurat.final$spleen <- PairedTag_data

seuratObj <- tissue.seurat.final$spleen
pdf("spleen/referenceMap.pdf")
assign.ratio <- refAnnotationRatio(seuratObj@meta.data, predictedLabels = "Zhang.2024.Science.EasySci", width = unit(7, "cm"), height = unit(6, "cm"), plot = T, group.by = "cellType", row_title = "Our dataset:Spleen")
dev.off()

saveRDS(tissue.seurat.final, "tissue.seurat.final.rds")
