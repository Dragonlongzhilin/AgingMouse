#' @description data processing
#' @Author: Zhilin Long
#' @date 2025-6-30
# conda activate PairedTagAnalysis

source("~/Functions/QualityControlSingle.R")
source("~/Functions/QualityControlMerge.R")
#### ---- BAT ---- ####
#### H3K27ac-H3K27me3
H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20241112-scPT-Natural-Aging-BAT-H3K27ac-H3K27me3", tissueName = "BAT", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K4me1-H3K4me3
H3K4me1_H3K4me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20241203-scPT-Natural-Aging-BAT-H3K4me1-H3K4me3", tissueName = "BAT", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K9me3-H3K36me3
H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20241125-scPT-Natural-Aging-BAT-H3K9me3-H3K36me3", tissueName = "BAT", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### Merge all asssay
MergeAllAssay <- QualityControlMerge(output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing/BAT", dataName = "PairedTag_data", tissueName = "BAT", merge_dirname = "MergeAllAssay", adjustedVariable = "histone_sample_batch")

#### ---- kidney ---- ####
#### H3K27ac-H3K27me3
H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20240606-scPT-Natural-Aging-kidney-H3K27ac-H3K27me3", tissueName = "kidney", batch = "20240606", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20250619-scPT-Natural-Aging-kidney-H3K27ac-H3K27me3", tissueName = "kidney", batch = "20250619", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K9me3-H3K36me3
H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20240715-scPT-Natural-Aging-kidney-H3K9me3-H3K36me3", tissueName = "kidney", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K4me1-H3K4me3
H3K4me1_H3K4me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20240708-scPT-Natural-Aging-kidney-H3K4me1-H3K4me3", tissueName = "kidney", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### Merge all asssay
MergeAllAssay <- QualityControlMerge(output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing/kidney", dataName = "PairedTag_data", tissueName = "kidney", merge_dirname = "MergeAllAssay", adjustedVariable = "histone_sample_batch", rmHistoneSample = c("H3K27me3:mouse138", "H3K27ac:mouse143", "H3K27ac:mouse144", "H3K27ac:mouse264", "H3K27ac:mouse265")) # remove 20240606-H3K27me3:mouse138 (low cells); 20250619-H3K27ac assay (low DNA counts)

#### ---- colon ---- ####
#### H3K27ac-H3K27me3
H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20240627-scPT-Natural-Aging-colon-H3K27ac-H3K27me3", tissueName = "colon", batch = "20240627", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20251103-scPT-Natural-Aging-colon-H3K27me3-H3K27me3", tissueName = "colon", batch = "20251103", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K9me3-H3K36me3
H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20240713-scPT-Natural-Aging-colon-H3K9me3-H3K36me3", tissueName = "colon", batch = "20240713", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20241121-scPT-Natural-Aging-colon-H3K9me3-H3K36me3", tissueName = "colon", batch = "20241121", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K4me1-H3K4me3
H3K4me1_H3K4me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20240823-scPT-Natural-Aging-colon-H3K4me1-H3K4me3", tissueName = "colon", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### Merge all asssay
MergeAllAssay <- QualityControlMerge(output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing/colon", dataName = "PairedTag_data", tissueName = "colon", merge_dirname = "MergeAllAssay", adjustedVariable = "histone_sample_batch", rmHistoneSample = c("mouse137", "mouse128")) # remove mouse137 with low cells and genes; remove mouse128 with low genes

#### ---- muscle ---- ####
#### H3K27ac-H3K27me3
H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20241230-scPT-Natural-Aging-muscle-H3K27ac-H3K27me3", tissueName = "muscle", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K9me3-H3K36me3
H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20250107-scPT-Natural-Aging-muscle-H3K9me3-H3K36me3", tissueName = "muscle", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K4me1-H3K4me3
H3K4me1_H3K4me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20250113-scPT-Natural-Aging-muscle-H3K4me1-H3K4me3", tissueName = "muscle", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### Merge all asssay
MergeAllAssay <- QualityControlMerge(output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing/muscle", dataName = "PairedTag_data", tissueName = "muscle", merge_dirname = "MergeAllAssay", adjustedVariable = "histone_sample_batch")

#### ---- spleen ---- ####
#### H3K27ac-H3K27me3
H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20251212-scPT-Natural-Aging-spleen-H3K27ac-H3K27me3", tissueName = "spleen", batch = "20251212", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1, nCount_RNA_limit = 100, nFeature_RNA_limit = 100, nCount_DNA_limit = 50)

#### H3K9me3-H3K36me3
H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20251013-scPT-Natural-Aging-spleen-H3K9me3-H3K36me3", tissueName = "spleen", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1, nCount_RNA_limit = 100, nFeature_RNA_limit = 100, nCount_DNA_limit = 50)

#### H3K4me1-H3K4me3
H3K4me1_H3K4me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20251113-scPT-Natural-Aging-spleen-H3K4me1-H3K4me3", tissueName = "spleen", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1, nCount_RNA_limit = 100, nFeature_RNA_limit = 100, nCount_DNA_limit = 50)

#### Merge all asssay
MergeAllAssay <- QualityControlMerge(output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing/spleen", dataName = "PairedTag_data", tissueName = "spleen", merge_dirname = "MergeAllAssay", adjustedVariable = "histone_sample_batch", rmHistoneSample = c("mouse126", "H3K27me3:mouse143", "H3K27me3:mouse144", "H3K27me3:mouse126", "H3K27me3:mouse127")) # mouse126: Imbalance in cell type ratios (B cell); 20251121-H3K27me3: low read depths

#### ---- liver ---- ####
#### H3K27ac-H3K27me3
H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20240620-scPT-Natural-Aging-liver-H3K27ac-H3K27me3", tissueName = "liver", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K9me3-H3K36me3
H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20240704-scPT-Natural-Aging-liver-H3K9me3-H3K36me3", tissueName = "liver", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K4me1-H3K4me3
H3K4me1_H3K4me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20240808-scPT-Natural-Aging-liver-H3K4me1-H3K4me3", tissueName = "liver", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### Merge all asssay
MergeAllAssay <- QualityControlMerge(output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing/liver", dataName = "PairedTag_data", tissueName = "liver", merge_dirname = "MergeAllAssay", adjustedVariable = "histone_sample_batch")

#### ---- heart ---- ####
#### H3K27ac-H3K27me3
H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20240930-scPT-Natural-Aging-heart-H3K27ac-H3K27me3", tissueName = "heart", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K9me3-H3K36me3
H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20241012-scPT-Natural-Aging-heart-H3K9me3-H3K36me3", tissueName = "heart", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K4me1-H3K4me3
H3K4me1_H3K4me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20241013-scPT-Natural-Aging-heart-H3K4me1-H3K4me3", tissueName = "heart", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### Merge all asssay
MergeAllAssay <- QualityControlMerge(output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing/heart", dataName = "PairedTag_data", tissueName = "heart", merge_dirname = "MergeAllAssay", adjustedVariable = "histone_sample_batch")

#### ---- lung ---- ####
#### H3K27ac-H3K27me3
H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20240808-scPT-Natural-Aging-lung-H3K27ac-H3K27me3", tissueName = "lung", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K9me3-H3K36me3
H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20240819-scPT-Natural-Aging-lung-H3K9me3-H3K36me3", tissueName = "lung", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K4me1-H3K4me3
H3K4me1_H3K4me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20240827-scPT-Natural-Aging-lung-H3K4me1-H3K4me3", tissueName = "lung", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### Merge all asssay
MergeAllAssay <- QualityControlMerge(output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing/lung", dataName = "PairedTag_data", tissueName = "lung", merge_dirname = "MergeAllAssay", adjustedVariable = "histone_sample_batch")

#### ---- brainCB ---- ####
#### H3K27ac-H3K27me3
H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20241130-scPT-Natural-Aging-brainCB-H3K27ac-H3K27me3", tissueName = "brainCB", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K9me3-H3K36me3
H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20241117-scPT-Natural-Aging-brainCB-H3K9me3-H3K36me3", tissueName = "brainCB", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K4me1-H3K4me3
H3K4me1_H3K4me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20250113-scPT-Natural-Aging-brainCB-H3K4me1-H3K4me3", tissueName = "brainCB", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### Merge all asssay
MergeAllAssay <- QualityControlMerge(output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing/brainCB", dataName = "PairedTag_data", tissueName = "brainCB", merge_dirname = "MergeAllAssay", adjustedVariable = "histone_sample_batch")

#### ---- brainFC ---- ####
#### H3K27ac-H3K27me3
H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20241207-scPT-Natural-Aging-brainFC-H3K27ac-H3K27me3", tissueName = "brainFC", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K9me3-H3K36me3
H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20240507-scPT-Natural-Aging-brainFC-H3K9me3-H3K36me3", tissueName = "brainFC", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K4me1-H3K4me3
H3K4me1_H3K4me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20241215-scPT-Natural-Aging-brainFC-H3K4me1-H3K4me3", tissueName = "brainFC", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### Merge all asssay
MergeAllAssay <- QualityControlMerge(output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing/brainFC", dataName = "PairedTag_data", tissueName = "brainFC", merge_dirname = "MergeAllAssay", adjustedVariable = "histone_sample_batch")

#### ---- brainHip ---- ####
#### H3K27ac-H3K27me3
H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20241110-scPT-Natural-Aging-brainHip-H3K27ac-H3K27me3", tissueName = "brainHip", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K9me3-H3K36me3
H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20241027-scPT-Natural-Aging-brainHip-H3K9me3-H3K36me3", tissueName = "brainHip", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K4me1-H3K4me3
H3K4me1_H3K4me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20241222-scPT-Natural-Aging-brainHip-H3K4me1-H3K4me3", tissueName = "brainHip", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### Merge all asssay
MergeAllAssay <- QualityControlMerge(output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing/brainHip", dataName = "PairedTag_data", tissueName = "brainHip", merge_dirname = "MergeAllAssay", adjustedVariable = "histone_sample_batch")

#### ---- testis ---- ####
#### H3K27ac-H3K27me3
H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20250524-scPT-Natural-Aging-testis-H3K27ac-H3K27me3", tissueName = "testis", batch = "20250524", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20260207-scPT-Natural-Aging-testis-H3K27ac-H3K27ac", tissueName = "testis", batch = "20260207", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K9me3-H3K36me3
H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20250410-scPT-Natural-Aging-testis-H3K9me3-H3K36me3", tissueName = "testis", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### H3K4me1-H3K4me3
H3K4me1_H3K4me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20250418-scPT-Natural-Aging-testis-H3K4me1-H3K4me3", tissueName = "testis", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1)

#### Merge all asssay
MergeAllAssay <- QualityControlMerge(output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing/testis", dataName = "PairedTag_data", tissueName = "testis", merge_dirname = "MergeAllAssay", adjustedVariable = "histone_sample_batch")

#### ---- skin ---- ####
#### H3K27ac-H3K27me3
H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20260121-scPT-Natural-Aging-skin-H3K27ac-H3K27me3", tissueName = "skin", batch = "20260121", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1, nCount_RNA_limit = 100, nFeature_RNA_limit = 100, nCount_DNA_limit = 50)

H3K27ac_H3K27me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20260315-scPT-Natural-Aging-skin-H3K27ac-H3K27me3", tissueName = "skin", batch = "20260315", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1, nCount_RNA_limit = 100, nFeature_RNA_limit = 100, nCount_DNA_limit = 50)

#### H3K9me3-H3K36me3
H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20251224-scPT-Natural-Aging-skin-H3K9me3-H3K36me3", tissueName = "skin", batch = "20251224", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1, nCount_RNA_limit = 100, nFeature_RNA_limit = 100, nCount_DNA_limit = 50)

H3K9me3_H3K36me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20260329-scPT-Natural-Aging-skin-H3K9me3-H3K36me3", tissueName = "skin", batch = "20260329", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1, nCount_RNA_limit = 100, nFeature_RNA_limit = 100, nCount_DNA_limit = 50)

#### H3K4me1-H3K4me3
H3K4me1_H3K4me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20250430-scPT-Natural-Aging-skin-H3K4me1-H3K4me3", tissueName = "skin", batch = "20250430", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1, nCount_RNA_limit = 100, nFeature_RNA_limit = 100, nCount_DNA_limit = 50)

H3K4me1_H3K4me3 <- QualityControlSingle(rawData.path = "~/paired_seq_tag/data/merge/20260412-scPT-Natural-Aging-skin-H3K4me1-H3K4me3", tissueName = "skin", batch = "20260412", output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing", doublet.rate = 0.1, nCount_RNA_limit = 100, nFeature_RNA_limit = 100, nCount_DNA_limit = 50)

#### Merge all asssay
source("~/Functions/QualityControlMerge.R")
MergeAllAssay <- QualityControlMerge(output.path = "~/AgingMouse/PairedTag/1.dataPreprocessing/skin", dataName = "PairedTag_data", tissueName = "skin", merge_dirname = "MergeAllAssay", adjustedVariable = "histone_sample_batch")

#### ---- Merge QC control ---- ####
library(Seurat)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(dplyr)
source("~/Functions/colorAssign.R")
setwd("~/AgingMouse/PairedTag/1.dataPreprocessing")
kidney <- readRDS("~/AgingMouse/PairedTag/1.dataPreprocessing/kidney/MergeAllAssay/PairedTag_data.rds")
BAT <- readRDS("~/AgingMouse/PairedTag/1.dataPreprocessing/BAT/MergeAllAssay/PairedTag_data.rds")
colon <- readRDS("~/AgingMouse/PairedTag/1.dataPreprocessing/colon/MergeAllAssay/PairedTag_data.rds")
muscle <- readRDS("~/AgingMouse/PairedTag/1.dataPreprocessing/muscle/MergeAllAssay/PairedTag_data.rds")
tissue.seurat.list <- list(kidney = kidney, BAT = BAT, colon = colon, muscle = muscle)
color.tissue <- ggsci::pal_igv("default")(length(names(tissue.seurat.list)))
names(color.tissue) <- names(tissue.seurat.list)
QC.info <- lapply(tissue.seurat.list, function(x){
  return(x@meta.data[,c("tissue", "age", "histone", "histone_age", "histone_age_sample", "nCount_RNA", "nFeature_RNA", "mt_ratio_RNA", "rp_ratio_RNA", "nCount_DNA", "cellType")])
})
QC.info <- do.call(rbind, QC.info)
QC.info <- plyr::arrange(QC.info, histone, tissue, age, histone_age_sample)
QC.info$label <- paste0(QC.info$tissue, ":", QC.info$histone_age_sample)
QC.info$label <- factor(QC.info$label, levels = unique(QC.info$label))

#### qc for each sample
p1 <- ggviolin(QC.info, x = "label", y = "nCount_RNA", color = "age", palette = color.age, xlab = "", ylab = "nCount_RNA", add = "median_iqr") + yscale("log10", .format = TRUE) + xlab("") + theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
p2 <- ggviolin(QC.info, x = "label", y = "nCount_DNA", color = "age", palette = color.age, xlab = "", ylab = "nCount_DNA", add = "median_iqr") + yscale("log10", .format = TRUE) + xlab("") + theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
p3 <- ggviolin(QC.info, x = "label", y = "nFeature_RNA", color = "age", palette = color.age, xlab = "", ylab = "nFeature_RNA", add = "median_iqr") + xlab("") + theme(legend.position="right", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
p4 <- ggviolin(QC.info, x = "label", y = "mt_ratio_RNA", color = "age", palette = color.age, xlab = "", ylab = "mt_ratio_RNA") + xlab("") + theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
p5 <- ggviolin(QC.info, x = "label", y = "rp_ratio_RNA", color = "age", palette = color.age, xlab = "", ylab = "rp_ratio_RNA", add = "median_iqr", add.params = list(size = 0.25)) + xlab("") + theme(legend.position="none", axis.text.x = element_text(color = "black", size = 10, angle = 60, hjust = 1, vjust = 1))
pdf("QC_control.pdf", height = 12, width = 25)
aligned_plots <- align_plots(p1, p2, p3, p4, p5, align = "v", axis = "lr")
grid.arrange(grobs = aligned_plots, ncol = 1, heights = c(1, 1, 1, 1, 3))
dev.off()

#### reads for each sample
DNA.reads <- tapply(QC.info$nCount_DNA, QC.info$label, sum)
RNA.reads <- tapply(QC.info$nCount_RNA, QC.info$label, sum)
plot.reads <- data.frame(label = names(DNA.reads), DNA.reads = DNA.reads, RNA.reads = RNA.reads)
plot.reads$sample <- gsub(".*_mouse", "mouse", plot.reads$label)
plot.reads$group <- gsub("_mouse.*", "", plot.reads$label)
plot.reads$DNA.reads <- round(plot.reads$DNA.reads / 1000000, 2)
plot.reads$RNA.reads <- round(plot.reads$RNA.reads / 1000000, 2)
plot.reads$cellNumber <- as.numeric(table(QC.info$label))
p1 <- ggbarplot(plot.reads, x = "group", y = "cellNumber", fill = "sample", color = "sample", palette = color.sample, xlab = "", ylab = "Cell number") + theme(legend.position="none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
p2 <- ggbarplot(plot.reads, x = "group", y = "RNA.reads", fill = "sample", color = "sample", palette = color.sample, xlab = "", ylab = "RNA reads (millions)") + theme(legend.position="right", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
p3 <- ggbarplot(plot.reads, x = "group", y = "DNA.reads", fill = "sample", color = "sample", palette = color.sample, xlab = "", ylab = "DNA reads (millions)") + theme(legend.position="none", axis.text.x = element_text(color = "black", size = 10, angle = 60, hjust = 1, vjust = 1))
pdf("QC_reads.pdf", height = 10, width = 25)
aligned_plots <- align_plots(p1, p2, p3, align = "v", axis = "lr")
grid.arrange(grobs = aligned_plots, ncol = 1, heights = c(1, 1, 2))
dev.off()

#### reads for each cell types
QC.info$label <- paste0(QC.info$tissue, ":", QC.info$cellType, "_", QC.info$histone)
rank.label <- as.data.frame(table(QC.info$label))
rank.label$tissue <- gsub(":.*", "", rank.label$Var1)
rank.label$histone <- gsub(".*_", "", rank.label$Var1)
rank.label <- plyr::arrange(rank.label, tissue, histone, Freq)
QC.info$label <- factor(QC.info$label, levels = rank.label$Var1)

DNA.reads <- tapply(QC.info$nCount_DNA, QC.info$label, sum)
RNA.reads <- tapply(QC.info$nCount_RNA, QC.info$label, sum)
plot.reads <- data.frame(label = names(DNA.reads), DNA.reads = DNA.reads, RNA.reads = RNA.reads)
plot.reads$cellType <- gsub("_.*", "", plot.reads$label)
plot.reads$histone <- gsub(".*_", "", plot.reads$label)
plot.reads$tissue <- gsub(":.*", "", plot.reads$label)
plot.reads$DNA.reads <- round(plot.reads$DNA.reads / 1000000, 2)
plot.reads$RNA.reads <- round(plot.reads$RNA.reads / 1000000, 2)
plot.reads$cellNumber <- as.numeric(table(QC.info$label))

p3 <- ggbarplot(plot.reads, x = "cellType", y = "DNA.reads", fill = "tissue", color = "tissue", palette = color.tissue, xlab = "", ylab = "DNA reads (millions)", orientation = "horiz", facet.by = "histone") + geom_vline(aes(xintercept = 10), lty = "dashed", lwd = 2, color = "red") + theme(legend.position="none", axis.text.x = element_text(color = "black", size = 10, angle = 60, hjust = 1, vjust = 1)) + facet_wrap(~ histone, nrow = 1, scales = "fixed")
pdf("QC_cellType_reads.pdf", height = 10, width = 20)
p3
dev.off()
