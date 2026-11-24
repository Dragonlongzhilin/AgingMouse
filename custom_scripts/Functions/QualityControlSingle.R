#' @description quality control of Paired-Tag-RNA data for single sample
#' @Author: Zhilin Long
#' @date 2025-3-26

# # settings
# rawData.path <- "~/qihongjian/projects/paired_seq_tag/data/merge/20241125-scPT-Natural-Aging-BAT-H3K9me3-H3K36me3"
# samples <- c("mouse137", "mouse136", "mouse122", "mouse123") # barcode represent the mouse sample (order)
# histoneMarks <- c("H3K27ac", "H3K27me3") # barcode represent the histone mark (order)
# tissueName <- "BAT"
# output.path <- "~/longzhilin/Projects/AgingMouse/PairedTag"

# meta file including the barcode, sample and histone information, for example
# bc1	histone	mouse	age
# 1	H3K9me3	mouse140	age3M
# 2	H3K9me3	mouse140	age3M
# 3	H3K9me3	mouse141	age3M
# 4	H3K9me3	mouse141	age3M
# 5	H3K36me3	mouse140	age3M
# 6	H3K36me3	mouse140	age3M
# 7	H3K36me3	mouse141	age3M
# 8	H3K36me3	mouse141	age3M
# 9	H3K9me3	mouse129	age27M
# 10	H3K9me3	mouse129	age27M
# 11	H3K9me3	mouse130	age27M
# 12	H3K9me3	mouse130	age27M
# 13	H3K36me3	mouse129	age27M
# 14	H3K36me3	mouse129	age27M
# 15	H3K36me3	mouse130	age27M
# 16	H3K36me3	mouse130	age27M


QualityControlSingle <- function(rawData.path, tissueName, output.path, batch = NULL, doublet.rate = 0.05, nPCs = 30, nCount_RNA_limit = 200, nFeature_RNA_limit = 200, nCount_DNA_limit = 50){
  require(Seurat)
  require(tidyverse)
  require(ggpubr)
  require(clustree)
  require(sctransform)
  require(ComplexHeatmap)
  require(circlize)
  require(scCustomize)
  require(reticulate)
  set.seed(101)
  require(future)
  plan("multicore", workers = 5)
  options(future.globals.maxSize = 100000 * 1024^2)
  set.resolutions <- seq(0.1, 1.2, by = 0.1)
  use_condaenv("~/anaconda3/envs/sctk-reticulate/bin/python", required = TRUE)
  scrublet <- import("scrublet")
  source("~/Functions/colorAssign.R")

  meta <- file.path("~/paired_seq_tag/code/machine-experiments/pipeline/config", paste0(basename(rawData.path), ".metadata.csv"))
  #### Reading meta file
  cat("1.1.Reading ", meta, "...\n")
  if (file.exists(meta)) {
    meta <- read.csv(meta, stringsAsFactors = F)
  } else {
    stop(meta, " file not exist!")
  }
  meta$bc1 <- as.character(meta$bc1)
  meta$age <- gsub("age3M", "m3", meta$age)
  meta$age <- gsub("age27M", "m27", meta$age)
  samples <- unique(meta$mouse)
  histoneMarks <- unique(meta$histone)

  if(length(histoneMarks) == 1){
    fileName <- histoneMarks
  }else{
    fileName <- paste0(histoneMarks[1], "_", histoneMarks[2])
  }

  if(is.null(batch)){
    workPath <- file.path(output.path, tissueName, fileName)
  }else{
    workPath <- file.path(output.path, tissueName, paste0(fileName, "_", batch))
  }

  if(!dir.exists(workPath)){
    dir.create(workPath, recursive = T, showWarnings = F)
  }
  setwd(workPath)
  
  #### load data
  cat("1.1.Loading ", tissueName, " data...\n")
  snRNA <- Read10X(data.dir = file.path(rawData.path, "rna_bam/RNA_matrix"))
  snRNA <- CreateSeuratObject(counts = snRNA, project = tissueName, min.cells = 3)
  snRNA[["mt_ratio_RNA"]] <- PercentageFeatureSet(snRNA, pattern = "^mt-")
  snRNA[["rp_ratio_RNA"]] <- PercentageFeatureSet(snRNA, pattern = "^Rpl|^Rps")
  snRNA$tissue <- tissueName
  snRNA$assay <- gsub("-scPT-Natural-Aging", "", basename(rawData.path))
  snRNA$batch <- gsub("-.*", "", snRNA$assay)

  barcode <- as.character(as.numeric(gsub(".*\\:", "", colnames(snRNA))))
  snRNA$barcode <- barcode
  # split histone
  snRNA$histone <- plyr::mapvalues(barcode, from = meta$bc1, to = meta$histone)
  snRNA$histone <- factor(snRNA$histone, levels = names(color.histone))
  snRNA$histone <- droplevels(snRNA$histone)

  # split age
  snRNA$age <- plyr::mapvalues(barcode, from = meta$bc1, to = meta$age)
  snRNA$age <- factor(snRNA$age, levels = names(color.age))

  # split sample
  snRNA$sample <- plyr::mapvalues(barcode, from = meta$bc1, to = meta$mouse)
  snRNA$sample <- factor(snRNA$sample, levels = samples)

  # add information
  snRNA$histone_sample <- paste0(snRNA$histone, "_", snRNA$sample)
  snRNA$histone_age <- paste0(snRNA$histone, "_", snRNA$age)
  snRNA$histone_age_sample <- paste0(snRNA$histone, "_", snRNA$age, "_", snRNA$sample)
  # add batch
  snRNA$histone_sample_batch <- paste0(snRNA$histone, "_", snRNA$sample, "_", snRNA$batch)
  snRNA$sample_batch <- paste0(snRNA$sample, "_", snRNA$batch)

  # add DNA counts
  DNA_dirs <- list.dirs(path = "~/paired_seq_tag/code/machine-experiments/DNA/counts/data", full.names = T, recursive = FALSE)
  DNA_dirs <- grep(paste0(tissueName, "-"), DNA_dirs, value = TRUE)
  DNA_dirs <- grep("_bins$", DNA_dirs, value = TRUE)
  DNA_dirs <- grep(paste0(as.character(unique(snRNA$histone)), collapse="|"), DNA_dirs, value = TRUE)
  if(length(DNA_dirs) > 0){
    info <- paste0("There are ", length(DNA_dirs), " directories in ~/paired_seq_tag/code/machine-experiments/DNA/counts/data.\n")
    cat(info)
    DNA_counts <- lapply(DNA_dirs, function(DNA_dir){
      DNA_matrix <- Read10X(data.dir = DNA_dir)
      DNA_count <- colSums(DNA_matrix)
      return(DNA_count)
    })
    DNA_counts <- unlist(DNA_counts)
    snRNA <- AddMetaData(snRNA, DNA_counts, col.name = "nCount_DNA")
  }else{
    info <- paste0("There are no directories in ~/paired_seq_tag/code/machine-experiments/DNA/counts/data. Please check!!!\n")
    cat(info)
    return(NULL)
  }

  #### ---- Quality control ---- ####
  pdf("1.1.RNA_PreFilter.pdf", height = 5, width = 12)
  print(VlnPlot(
    object = snRNA,
    features = c("nCount_RNA", "nCount_DNA", "nFeature_RNA", "mt_ratio_RNA", "rp_ratio_RNA"),
    ncol = 5,
    cols = color.sample,
    pt.size = 0,
    group.by = 'histone_age_sample'
  ) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  print(VlnPlot(
    object = snRNA,
    features = c("nCount_RNA", "nCount_DNA", "nFeature_RNA"),
    ncol = 5,
    cols = color.sample,
    pt.size = 0,
    group.by = 'histone_age_sample'
  ) & yscale("log10", .format = TRUE) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  dev.off()

  ## filtering cells
  cat("1.2.Filtering cells with nCount_RNA >", nCount_RNA_limit, " & nFeature_RNA >", nFeature_RNA_limit, " & nCount_DNA >", nCount_DNA_limit, " & mt_ratio_RNA < 5 ...\n")
  snRNA.pro <- subset(snRNA, subset = nCount_RNA > nCount_RNA_limit & nFeature_RNA > nFeature_RNA_limit & nCount_DNA > nCount_DNA_limit & mt_ratio_RNA < 5)

  ## remove doublets
  source("~/Functions/doubletDetect.R")
  cat("1.2.Remove doublets...\n")
  groups <- unique(snRNA.pro$barcode)
  pdf("1.2.RNA_Doublets.pdf", height = 10, width = 15)
  test.list <- lapply(groups, function(group){
    cat("  ", group, "...\n")
    test <- subset(snRNA.pro, subset = barcode == group)
    test <- NormalizeData(test, verbose = FALSE)
    test <- FindVariableFeatures(test, verbose = FALSE)
    test <- ScaleData(test, vars.to.regress = c("mt_ratio_RNA", "nCount_RNA"), verbose = FALSE)
    test <- RunPCA(test, verbose = FALSE)
    test <- FindNeighbors(test, reduction = "pca", dims = 1:nPCs, verbose = FALSE)
    test <- FindClusters(test, resolution = 0.3, verbose = FALSE)
    test <- RunUMAP(test, dims = 1:nPCs, reduction = "pca", verbose = FALSE)
    ElbowPlot(test, ndims = 50)

    # remove top 1% discrete cells
    nCount_RNA.top <- quantile(test$nCount_RNA, probs = 0.99)
    nFeature_RNA.top <- quantile(test$nFeature_RNA, probs = 0.99)
    nCount_DNA.top <- quantile(test$nCount_DNA, probs = 0.99)
    discrete.group <- rep("Pass", ncol(test))
    discrete.group[test$nCount_RNA > nCount_RNA.top | test$nFeature_RNA > nFeature_RNA.top | test$nCount_DNA > nCount_DNA.top] <- "Fail"
    test$discrete.group <- discrete.group

    # DoubletFinder method
    test <- FunctionDoubletFinder(Seurat.object = test, PCs = 1:nPCs, doublet.rate = doublet.rate, annotation = "seurat_clusters", num.cores = 10)
    test$DoubletFinder.prediction <- factor(test$DoubletFinder.prediction, levels = c("Singlet", "Doublet"))
    DoubletFinder.res <- length(which(test$DoubletFinder.prediction == "Doublet"))
    DoubletFinder.res <- paste0(DoubletFinder.res, " in ", ncol(test), " cells", " (", round(DoubletFinder.res / ncol(test) * 100, 2), "%)")

    # scrublet method
    counts_matrix <- GetAssayData(test, slot = "counts")
    scrub <- scrublet$Scrublet(t(counts_matrix), expected_doublet_rate = doublet.rate)
    scrublet.score <- scrub$scrub_doublets(n_prin_comps = as.integer(nPCs))
    scrublet.prediction <- rep("Singlet", ncol(counts_matrix))
    scrublet.prediction[scrublet.score[[2]]] <- "Doublet"
    test$scrublet.prediction <- factor(scrublet.prediction, levels = c("Singlet", "Doublet"))
    test$scrublet.score <- scrublet.score[[1]]
    scrublet.res <- length(which(test$scrublet.prediction == "Doublet"))
    scrublet.res <- paste0(scrublet.res, " in ", ncol(test), " cells", " (", round(scrublet.res / ncol(test) * 100, 2), "%)")

    p1 <- DimPlot_scCustom(test, group.by = "seurat_clusters", label = T, raster=FALSE) + theme(aspect.ratio = 1)
    p2 <- FeaturePlot_scCustom(test, features = "nCount_RNA") + theme(aspect.ratio = 1)
    p3 <- FeaturePlot_scCustom(test, features = "nCount_DNA") + theme(aspect.ratio = 1)
    p4 <- DimPlot_scCustom(test, group.by = "DoubletFinder.prediction", colors_use = c("#E7E7E7", "#B41008"), raster=FALSE) + ggtitle(paste0("Doublets (DoubletFinder)\n", DoubletFinder.res)) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))
    p5 <- DimPlot_scCustom(test, group.by = "scrublet.prediction", colors_use = c("#E7E7E7", "#B41008"), raster=FALSE) + ggtitle(paste0("Doublets (scrublet)\n", scrublet.res)) + theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5))

    p11 <- VlnPlot_scCustom(seurat_object = test, features = "nCount_RNA", group.by = "DoubletFinder.prediction", plot_boxplot = TRUE, colors_use = c("#E7E7E7", "#B41008")) & yscale("log10", .format = TRUE) & NoLegend() & xlab("")
    p22 <- VlnPlot_scCustom(seurat_object = test, features = "nCount_DNA", group.by = "DoubletFinder.prediction", plot_boxplot = TRUE, colors_use = c("#E7E7E7", "#B41008")) & yscale("log10", .format = TRUE) & NoLegend() & xlab("")
    p33 <- VlnPlot_scCustom(seurat_object = test, features = "nCount_RNA", group.by = "scrublet.prediction", plot_boxplot = TRUE, colors_use = c("#E7E7E7", "#B41008")) & yscale("log10", .format = TRUE) & NoLegend() & xlab("")
    p44 <- VlnPlot_scCustom(seurat_object = test, features = "nCount_DNA", group.by = "scrublet.prediction", plot_boxplot = TRUE, colors_use = c("#E7E7E7", "#B41008")) & yscale("log10", .format = TRUE) & NoLegend() & xlab("")
    p6 <- ggarrange(p11, p22, p33, p44, nrow = 2, ncol = 2, align = "hv")

    p <- ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2, align = "hv")
    annotated_p <- annotate_figure(p, top = text_grob(group, size = 14))
    print(annotated_p)
    return(test)
  })
  dev.off()
  test <- merge(x = test.list[[1]], y = test.list[-1])

  snRNA.pro <- AddMetaData(snRNA.pro, metadata = test$DoubletFinder.prediction, col.name = "DoubletFinder.prediction")
  snRNA.pro <- AddMetaData(snRNA.pro, metadata = test$discrete.group, col.name = "discrete.group")
  snRNA.pro <- AddMetaData(snRNA.pro, metadata = test$scrublet.prediction, col.name = "scrublet.prediction")
  snRNA.pro <- AddMetaData(snRNA.pro, metadata = test$scrublet.score, col.name = "scrublet.score")

  pdf("1.2.RNA_AfterFilter.pdf", height = 5, width = 12)
  print(VlnPlot(
    object = snRNA.pro,
    features = c("nCount_RNA", "nCount_DNA", "nFeature_RNA", "mt_ratio_RNA", "rp_ratio_RNA"),
    ncol = 5,
    cols = color.sample,
    pt.size = 0,
    group.by = 'histone_age_sample'
  ) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  print(VlnPlot(
    object = snRNA.pro,
    features = c("nCount_RNA", "nCount_DNA", "nFeature_RNA"),
    ncol = 5,
    cols = color.sample,
    pt.size = 0,
    group.by = 'histone_age_sample'
  ) & yscale("log10", .format = TRUE) & xlab("") & theme(title = element_text(size=10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)))
  dev.off()

  # plot cell number
  histone_age_sample.num <- as.data.frame(table(snRNA$histone_age_sample))
  histone_age_sample.num$histone_age <- gsub("_mouse.*", "", histone_age_sample.num$Var1)
  histone_age_sample.num$histone_age <- factor(histone_age_sample.num$histone_age, levels = names(color.histone_age))
  histone_age_sample.num$histone_age <- droplevels(histone_age_sample.num$histone_age)
  histone_age_sample.num$Var1 <- factor(histone_age_sample.num$Var1, levels = histone_age_sample.num$Var1)
  p1 <- ggbarplot(histone_age_sample.num, x = "Var1", y = "Freq", fill = "histone_age", lab.pos = "out", label = T, palette = color.histone_age, xlab = "", ylab = "Cell number") + ggtitle(paste0(tissueName, ": ", ncol(snRNA), " cells")) + rotate_x_text(90) + theme(legend.position="none", plot.title = element_text(hjust = 0.5))

  histone_age_sample.num <- as.data.frame(table(snRNA.pro$histone_age_sample))
  histone_age_sample.num$histone_age <- gsub("_mouse.*", "", histone_age_sample.num$Var1)
  histone_age_sample.num$histone_age <- factor(histone_age_sample.num$histone_age, levels = names(color.histone_age))
  histone_age_sample.num$histone_age <- droplevels(histone_age_sample.num$histone_age)
  histone_age_sample.num$Var1 <- factor(histone_age_sample.num$Var1, levels = histone_age_sample.num$Var1)
  p2 <- ggbarplot(histone_age_sample.num, x = "Var1", y = "Freq", fill = "histone_age", lab.pos = "out", label = T, palette = color.histone_age, xlab = "", ylab = "Cell number") + ggtitle(paste0(tissueName, ": ", ncol(snRNA.pro), " cells")) + rotate_x_text(90) + theme(legend.position="none", plot.title = element_text(hjust = 0.5))

  pdf("1.2.RNA_CellNumber.pdf", height = 6)
  print(ggarrange(p1, p2, nrow = 1, ncol = 2, align = "hv"))
  dev.off()

  p1 <- FeatureScatter(snRNA.pro, group.by = "histone_age_sample", feature1 = "nCount_RNA", feature2 = "mt_ratio_RNA", cols = color.sample, raster=FALSE) + theme(aspect.ratio = 1)
  p2 <- FeatureScatter(snRNA.pro, group.by = "histone_age_sample", feature1 = "nCount_RNA", feature2 = "nFeature_RNA", cols = color.sample, raster=FALSE) + theme(aspect.ratio = 1)
  p3 <- FeatureScatter(snRNA.pro, group.by = "histone_age_sample", feature1 = "nCount_RNA", feature2 = "nCount_DNA", cols = color.sample, raster=FALSE) + theme(aspect.ratio = 1)
  pdf("1.2.RNA_FeatureScatter.pdf")
  print(p1)
  print(p2)
  print(p3)
  dev.off()

  #### ---- Cluster and batch correction ---- ####
  # split datasets and process without integration
  cat("1.3.Clustering and batch correction...\n")
  snRNA.pro[["RNA"]] <- split(snRNA.pro[["RNA"]], f = snRNA.pro$histone_sample_batch)
  snRNA.pro <- NormalizeData(snRNA.pro, verbose = FALSE)
  snRNA.pro <- FindVariableFeatures(snRNA.pro, verbose = FALSE)
  snRNA.pro <- ScaleData(snRNA.pro, vars.to.regress = c("mt_ratio_RNA", "nCount_RNA"), verbose = FALSE)
  snRNA.pro <- RunPCA(snRNA.pro, verbose = FALSE)
  snRNA.pro <- RunUMAP(snRNA.pro, dims = 1:nPCs, reduction = "pca", reduction.name = "umap.uncorrected", verbose = FALSE)
  snRNA.pro <- IntegrateLayers(object = snRNA.pro, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", npcs = nPCs, verbose = FALSE)
  # re-join layers before performing any differential expression analysis
  snRNA.pro[["RNA"]] <- JoinLayers(snRNA.pro[["RNA"]])
  snRNA.pro <- FindNeighbors(snRNA.pro, reduction = "harmony", dims = 1:nPCs, verbose = FALSE)
  snRNA.pro <- FindClusters(snRNA.pro, resolution = set.resolutions, verbose = FALSE)
  snRNA.pro <- RunUMAP(snRNA.pro, dims = 1:nPCs, reduction = "harmony", verbose = FALSE)
  pdf("1.3.RNA_ElbowPlotClustree.pdf", height = 8)
  ElbowPlot(snRNA.pro, ndims = 50)
  p <- clustree(snRNA.pro)
  print(p)
  dev.off()

  pdf("1.3.RNA_ClusterPhenotype.pdf", width = 15, height = 10)
  p1 <- DimPlot(snRNA.pro, reduction = "umap.uncorrected", group.by = "age", cols = color.age, raster=FALSE) + theme(aspect.ratio = 1)
  p2 <- DimPlot(snRNA.pro, reduction = "umap.uncorrected", group.by = "histone", cols = color.histone, raster=FALSE) + theme(aspect.ratio = 1)
  p3 <- DimPlot(snRNA.pro, reduction = "umap.uncorrected", group.by = "sample", cols = color.sample, raster=FALSE) + theme(aspect.ratio = 1)
  p4 <- DimPlot(snRNA.pro, reduction = "umap.uncorrected", group.by = "barcode", cols = color.cellType, raster=FALSE) + theme(aspect.ratio = 1)
  print(ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, align = "hv"))

  p1 <- DimPlot(snRNA.pro, reduction = "umap", group.by = "age", cols = color.age, raster=FALSE) + theme(aspect.ratio = 1)
  p2 <- DimPlot(snRNA.pro, reduction = "umap", group.by = "histone", cols = color.histone, raster=FALSE) + theme(aspect.ratio = 1)
  p3 <- DimPlot(snRNA.pro, reduction = "umap", group.by = "sample", cols = color.sample, raster=FALSE) + theme(aspect.ratio = 1)
  p4 <- DimPlot(snRNA.pro, reduction = "umap", group.by = "barcode", cols = color.cellType, raster=FALSE) + theme(aspect.ratio = 1)
  print(ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, align = "hv"))

  p1 <- DimPlot(snRNA.pro, reduction = "umap.uncorrected", group.by = "histone_sample_batch", cols = color.sample, raster=FALSE) + theme(aspect.ratio = 1)
  p2 <- DimPlot(snRNA.pro, reduction = "umap", group.by = "histone_sample_batch", cols = color.sample, raster=FALSE) + theme(aspect.ratio = 1)
  print(ggarrange(p1, p2, ncol = 1, nrow = 2, align = "hv"))
  dev.off()

  pdf("1.3.RNA_Cluster_QCmetrics.pdf", height = 8, width = 8)
  p <- FeaturePlot(snRNA.pro, reduction = "umap", features = c("nCount_RNA", "nCount_DNA", "nFeature_RNA", "mt_ratio_RNA", "rp_ratio_RNA"), ncol = 2, cols = c("lightgrey", "red")) & theme(plot.title = element_text(size = 11), axis.title = element_text(size = 10), axis.text = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 10), aspect.ratio = 1)
  print(p)
  p <- FeaturePlot_scCustom(snRNA.pro, reduction = "umap", features = c("nCount_RNA", "nCount_DNA", "nFeature_RNA", "mt_ratio_RNA", "rp_ratio_RNA"), num_columns = 2) & theme(plot.title = element_text(size = 11), axis.title = element_text(size = 10), axis.text = element_text(size = 10), legend.text = element_text(size = 10), legend.title = element_text(size = 10), aspect.ratio = 1)
  print(p)
  dev.off()

  #### ---- compare the depth between pre- and after-filtering ---- ####
  ## pre-filtering
  depths.info <- as.data.frame(table(snRNA@meta.data[,c("histone_age_sample", "assay")]))
  colnames(depths.info) <- c("histone_age_sample", "assay", "Cell_number")
  depths.info <- depths.info[which(depths.info$Cell_number > 0),]
  RNA.sample.depths <- tapply(snRNA$nCount_RNA, snRNA$histone_age_sample, sum)
  DNA.sample.depths <- tapply(snRNA$nCount_DNA, snRNA$histone_age_sample, sum)
  depths.info$RNA.reads <- RNA.sample.depths[depths.info$histone_age_sample]
  depths.info$DNA.reads <- DNA.sample.depths[depths.info$histone_age_sample]

  ## after-filtering
  RNA.sample.depths2 <- tapply(snRNA.pro$nCount_RNA, snRNA.pro$histone_age_sample, sum)
  DNA.sample.depths2 <- tapply(snRNA.pro$nCount_DNA, snRNA.pro$histone_age_sample, sum)
  depths.info$Cell_number.filter1 <- as.numeric(table(snRNA.pro$histone_age_sample)[depths.info$histone_age_sample])
  depths.info$RNA.reads.filter1 <- RNA.sample.depths2[depths.info$histone_age_sample]
  depths.info$DNA.reads.filter1 <- DNA.sample.depths2[depths.info$histone_age_sample]
  saveRDS(depths.info, "depths.info.rds")

  saveRDS(snRNA.pro, paste0(fileName,".rds"))
  return(snRNA.pro)
}
