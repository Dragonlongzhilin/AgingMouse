#' @description merge multiple samples of Paired-Tag-RNA data
#' @Author: Zhilin Long
#' @date 2025-6-12

#' @param rmHistoneSample: histone + sample or histone or sample, e.g. "H3K9me3:mouse140", "H3K9me3", "mouse140"

QualityControlMerge <- function(output.path, dataName, tissueName, merge_dirname = "MergeAllAssay", adjustedVariable = "histone_sample_batch", rmPattern = "keepDoublets|depths.info", rmHistoneSample = NULL, nPCs = 30){
  require(Seurat)
  require(tidyverse)
  require(ggpubr)
  require(ggsci)
  require(clustree)
  require(sctransform)
  require(scCustomize)
  require(ComplexHeatmap)
  require(circlize)
  set.seed(101)
  require(future)
  plan("multicore", workers = 5)
  options(future.globals.maxSize = 110000 * 1024^2)
  set.resolutions <- seq(0.1, 1.2, by = 0.1)
  source("~/Functions/colorAssign.R")
  merge_dir <- file.path(output.path, merge_dirname)
  if (!dir.exists(merge_dir)) {
    dir.create(merge_dir, recursive = T, showWarnings = F)
  }
  setwd(merge_dir)

  #### ---- Merge all assay ---- ####
  Histone.datasets <- list.files(output.path, ".rds", full.names = T, recursive = T)
  Histone.datasets <- Histone.datasets[grep(rmPattern, Histone.datasets, invert = T)]
  Histone.datasets <- Histone.datasets[grep("H3K27me3|H3K4me3|H3K4me1|H3K9me3|H3K36me3|H3K27ac", Histone.datasets)]
  cat("Merge assay: ", Histone.datasets, "\n")
  Histone.datasets.list <- lapply(Histone.datasets, function(x){
    Histone.dataset <- readRDS(x)
    return(Histone.dataset)
  })
  samples <- gsub(".rds", "", basename(Histone.datasets))
  PairedTag_data <- merge(x = Histone.datasets.list[[1]], y = Histone.datasets.list[-1], add.cell.ids = samples)

  # remove low quality assay
  if(!is.null(rmHistoneSample)){
    rmCells <- lapply(rmHistoneSample, function(x){
      rmHistoneSample.split <- unlist(strsplit(x, ":")) # histone + sample
      if(length(rmHistoneSample.split) == 1){
        if(rmHistoneSample.split %in% names(color.histone)){
          rmCell <- colnames(PairedTag_data)[which(PairedTag_data$histone == rmHistoneSample.split)]
        }else{
          rmCell <- colnames(PairedTag_data)[which(PairedTag_data$sample == rmHistoneSample.split)]
        }
      }else{
        rmCell <- colnames(PairedTag_data)[which(PairedTag_data$histone == rmHistoneSample.split[1] & PairedTag_data$sample == rmHistoneSample.split[2])]
      }
      return(rmCell)
    })
    rmCells <- unique(unlist(rmCells))
    PairedTag_data <- subset(PairedTag_data, cells = setdiff(colnames(PairedTag_data), rmCells))
  }
  PairedTag_data$age <- factor(PairedTag_data$age, levels = names(color.age))
  PairedTag_data$age <- droplevels(PairedTag_data$age)
  PairedTag_data$histone <- factor(PairedTag_data$histone, levels = names(color.histone))
  PairedTag_data$histone <- droplevels(PairedTag_data$histone)
  PairedTag_data$barcode <- factor(PairedTag_data$barcode, levels = as.character(sort(as.numeric(unique(PairedTag_data$barcode)))))
  PairedTag_data$barcode <- droplevels(PairedTag_data$barcode)

  # remove top 1% discrete cells for each barcode in each assay
  PairedTag_data <- subset(PairedTag_data, subset = discrete.group == "Pass")

  histone_age_sample.num <- as.data.frame(table(PairedTag_data$histone_age_sample))
  histone_age_sample.num$histone_age <- gsub("_mouse.*", "", histone_age_sample.num$Var1)
  histone_age_sample.num$histone_age <- factor(histone_age_sample.num$histone_age, levels = names(color.histone_age))
  histone_age_sample.num$histone_age <- droplevels(histone_age_sample.num$histone_age)
  histone_age_sample.num$Var1 <- factor(histone_age_sample.num$Var1, levels = histone_age_sample.num$Var1)
  p1 <- ggbarplot(histone_age_sample.num, x = "Var1", y = "Freq", fill = "histone_age", lab.pos = "out", label = T, palette = color.histone_age, xlab = "", ylab = "Cell number", lab.size = 3) + ggtitle(paste0(tissueName, ": ", ncol(PairedTag_data), " cells")) + rotate_x_text(90) + theme(legend.position="right", plot.title = element_text(hjust = 0.5))
  pdf("1.2.CellNumber.pdf", height = 6, width = 10)
  print(p1)
  dev.off()

  group <- names(table(PairedTag_data$histone_age_sample))
  group <- gsub("_mouse.*", "", group)
  colors <- color.histone_age[match(group, names(color.histone_age))]
  names(colors) <- names(table(PairedTag_data$histone_age_sample))
  pdf("1.2.AfterFilter.pdf", height = 5)
  print(VlnPlot(
    object = PairedTag_data,
    features = "nCount_RNA",
    ncol = 1,
    cols = colors,
    pt.size = 0,
    group.by = 'histone_age_sample'
  ) & xlab("") & theme(legend.position="none", title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  print(VlnPlot(
    object = PairedTag_data,
    features = "nCount_DNA",
    ncol = 1,
    cols = colors,
    pt.size = 0,
    group.by = 'histone_age_sample'
  ) & xlab("") & theme(legend.position="none", title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  print(VlnPlot(
    object = PairedTag_data,
    features = "nFeature_RNA",
    ncol = 1,
    cols = colors,
    pt.size = 0,
    group.by = 'histone_age_sample'
  ) & xlab("") & theme(legend.position="none", title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  print(VlnPlot(
    object = PairedTag_data,
    features = "mt_ratio_RNA",
    ncol = 1,
    cols = colors,
    pt.size = 0,
    group.by = 'histone_age_sample'
  ) & xlab("") & theme(legend.position="none", title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  print(VlnPlot(
    object = PairedTag_data,
    features = "rp_ratio_RNA",
    ncol = 1,
    cols = colors,
    pt.size = 0,
    group.by = 'histone_age_sample'
  ) & xlab("") & theme(legend.position="none", title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  print(VlnPlot(
    object = PairedTag_data,
    features = "nCount_RNA",
    ncol = 1,
    cols = colors,
    pt.size = 0,
    group.by = 'histone_age_sample'
  ) & yscale("log10", .format = TRUE) & xlab("") & theme(legend.position="none", title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  print(VlnPlot(
    object = PairedTag_data,
    features = "nCount_DNA",
    ncol = 1,
    cols = colors,
    pt.size = 0,
    group.by = 'histone_age_sample'
  ) & yscale("log10", .format = TRUE) & xlab("") & theme(legend.position="none", title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  print(VlnPlot(
    object = PairedTag_data,
    features = "nFeature_RNA",
    ncol = 1,
    cols = colors,
    pt.size = 0,
    group.by = 'histone_age_sample'
  ) & yscale("log10", .format = TRUE) & xlab("") & theme(legend.position="none", title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  dev.off()

  p1 <- FeatureScatter(PairedTag_data, group.by = "histone_age_sample", feature1 = "nCount_RNA", feature2 = "mt_ratio_RNA", cols = color.sample, raster=FALSE) + theme(aspect.ratio = 1)
  p2 <- FeatureScatter(PairedTag_data, group.by = "histone_age_sample", feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = color.sample, raster=FALSE) + theme(aspect.ratio = 1)
  p3 <- FeatureScatter(PairedTag_data, group.by = "histone_age_sample", feature1 = "nCount_RNA", feature2 = "nCount_DNA", cols = color.sample, raster=FALSE) + theme(aspect.ratio = 1)
  pdf("1.2.FeatureScatter.pdf", height = 10, width = 10)
  print(p1)
  print(p2)
  print(p3)
  dev.off()

  #### ---- Cluster and batch correction ---- ####
  PairedTag_data <- JoinLayers(PairedTag_data)
  PairedTag_data <- DietSeurat(PairedTag_data, assays = "RNA", layers = c("counts"))
  PairedTag_data[["RNA"]] <- split(PairedTag_data[["RNA"]], f = PairedTag_data@meta.data[,adjustedVariable])
  PairedTag_data <- NormalizeData(PairedTag_data, verbose = FALSE)
  PairedTag_data <- FindVariableFeatures(PairedTag_data, verbose = FALSE)
  PairedTag_data <- ScaleData(PairedTag_data, vars.to.regress = c("mt_ratio_RNA", "nCount_RNA"), verbose = FALSE)
  PairedTag_data <- RunPCA(PairedTag_data, verbose = FALSE)
  PairedTag_data <- RunUMAP(PairedTag_data, dims = 1:nPCs, reduction = "pca", reduction.name = "umap.uncorrected", verbose = FALSE)
  PairedTag_data <- IntegrateLayers(object = PairedTag_data, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
  # re-join layers before performing any differential expression analysis
  PairedTag_data[["RNA"]] <- JoinLayers(PairedTag_data[["RNA"]])
  PairedTag_data <- FindNeighbors(PairedTag_data, reduction = "harmony", dims = 1:nPCs, verbose = FALSE)
  PairedTag_data <- FindClusters(PairedTag_data, resolution = set.resolutions, verbose = FALSE)
  PairedTag_data <- RunUMAP(PairedTag_data, dims = 1:nPCs, reduction = "harmony", verbose = FALSE)

  pdf("1.3.RNA_ElbowPlotClustree.pdf", height = 8)
  print(ElbowPlot(PairedTag_data, ndims = 50))
  p <- clustree(PairedTag_data)
  print(p)
  dev.off()

  pdf("1.3.RNA_ClusterPhenotype.pdf", width = 20, height = 15)
  p1 <- DimPlot(PairedTag_data, reduction = "umap.uncorrected", group.by = "age", cols = color.age, raster=FALSE) + theme(aspect.ratio = 1)
  p2 <- DimPlot(PairedTag_data, reduction = "umap.uncorrected", group.by = "histone", cols = color.histone, raster=FALSE) + theme(aspect.ratio = 1)
  p3 <- DimPlot(PairedTag_data, reduction = "umap.uncorrected", group.by = "sample", cols = color.sample, raster=FALSE) + theme(aspect.ratio = 1)
  p4 <- DimPlot(PairedTag_data, reduction = "umap.uncorrected", group.by = "barcode", cols = color.cellType, raster=FALSE) + theme(aspect.ratio = 1)
  print(ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, align = "hv"))

  p1 <- DimPlot(PairedTag_data, reduction = "umap", group.by = "age", cols = color.age, raster=FALSE) + theme(aspect.ratio = 1)
  p2 <- DimPlot(PairedTag_data, reduction = "umap", group.by = "histone", cols = color.histone, raster=FALSE) + theme(aspect.ratio = 1)
  p3 <- DimPlot(PairedTag_data, reduction = "umap", group.by = "sample", cols = color.sample, raster=FALSE) + theme(aspect.ratio = 1)
  p4 <- DimPlot(PairedTag_data, reduction = "umap", group.by = "barcode", cols = color.cellType, raster=FALSE) + theme(aspect.ratio = 1)
  print(ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, align = "hv"))

  p1 <- DimPlot(PairedTag_data, reduction = "umap.uncorrected", group.by = "histone_age_sample", cols = color.sample, raster=FALSE) + theme(aspect.ratio = 1)
  p2 <- DimPlot(PairedTag_data, reduction = "umap", group.by = "histone_age_sample", cols = color.sample, raster=FALSE) + theme(aspect.ratio = 1)
  print(ggarrange(p1, p2, ncol = 1, nrow = 2, align = "hv"))
  dev.off()

  pdf("1.3.RNA_Cluster_QCmetrics.pdf", height = 8, width = 10)
  p <- FeaturePlot(PairedTag_data, reduction = "umap", features = c("nCount_RNA", "nCount_DNA", "nFeature_RNA", "mt_ratio_RNA", "rp_ratio_RNA"), ncol = 3, cols = c("lightgrey", "red")) & theme(plot.title = element_text(size = 11), axis.title = element_text(size = 10), axis.text = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 10), aspect.ratio = 1)
  print(p)
  p <- FeaturePlot_scCustom(PairedTag_data, reduction = "umap", features = c("nCount_RNA", "nCount_DNA", "nFeature_RNA", "mt_ratio_RNA", "rp_ratio_RNA"), num_columns = 3) & theme(plot.title = element_text(size = 11), axis.title = element_text(size = 10), axis.text = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 10), aspect.ratio = 1)
  print(p)
  dev.off()
  saveRDS(PairedTag_data, paste0(dataName, ".rds"))

  #### ---- compare the depth between pre- and after-filtering ---- ####
  ## load data
  depths.info.list <- list.files(output.path, "^depths.info.rds", full.names = T, recursive = T)
  depths.info.list <- depths.info.list[grep(merge_dirname, depths.info.list, invert = T)]
  depths.info <- lapply(depths.info.list, function(x) readRDS(x))
  depths.info <- do.call("rbind", depths.info)

  ## after-filtering
  RNA.sample.depths3 <- tapply(PairedTag_data$nCount_RNA, PairedTag_data$histone_age_sample, sum)
  DNA.sample.depths3 <- tapply(PairedTag_data$nCount_DNA, PairedTag_data$histone_age_sample, sum)
  Cell_number <- table(PairedTag_data$histone_age_sample)
  idx <- match(depths.info$histone_age_sample, names(Cell_number))
  depths.info$Cell_number.filter2 <- as.numeric(Cell_number[idx])
  depths.info$RNA.reads.filter2 <- RNA.sample.depths3[idx]
  depths.info$DNA.reads.filter2 <- DNA.sample.depths3[idx]
  saveRDS(depths.info, "merge.depths.info.rds")
  return(PairedTag_data)
}
