#' @description HUVEC experiment
#' @Author: Zhilin Long
#' @date 2026-8-25

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
setwd("~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment")

#### ---- quantify the expression ---- ####
## conda activate BigDataProcessing
#### replicative senescence
cd ~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/RNAseq/replicative_senescence
nohup bash run_rnaseq.sh > run_rnaseq.log 2>&1 &

#### GSK126 treatment
cd ~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/RNAseq/GSK126_treatment
nohup bash run_rnaseq.sh > run_rnaseq.log 2>&1 &

#### convert the ENSG to gene name
library(readr)
library(tidyverse)
library(rtracklayer)
gene_mapping <- readRDS("~/longzhilin/AnalysisCode/Data/ID_Mapping/Human/hg38.gencode.v38.ENSG_geneName.rds") # ~/longzhilin/AnalysisCode/Data/ID_Mapping/Human/ENSG_geneSymbol.R
setwd("~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/RNAseq/GSK126_treatment")
gene_counts <- read.table("combined.counts", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gene_counts$Geneid <- sub("\\..*$", "", gene_counts$Geneid)
gene_counts$gene_length_kb <- gene_counts$Length / 1000
sample_cols <- grep(".sorted.bam", colnames(gene_counts))

## remove the gene with zero counts
gene_counts <- gene_counts[rowSums(gene_counts[,sample_cols]) > 0,]

## ID mapping
gene_counts <- gene_counts %>% dplyr::rename(gene_id = Geneid) %>% filter(gene_id %in% gene_mapping$gene_id)
counts_annotated <- gene_counts %>% left_join(gene_mapping %>% dplyr::select(gene_id, gene_name), by = "gene_id")

## calculate TPM
rpk <- sweep(as.matrix(counts_annotated[, sample_cols]), 1, counts_annotated$gene_length_kb, "/")
tpm_matrix <- sweep(rpk, 2, colSums(rpk), "/") * 1e6
tpm <- data.frame(gene_id = counts_annotated$gene_id, gene_name = counts_annotated$gene_name)
tpm_result <- cbind(
  tpm %>% dplyr::select(gene_id, gene_name),
  as.data.frame(tpm_matrix, check.names = FALSE)
)
gene_counts <- counts_annotated[,grep("HM0", colnames(counts_annotated))]
rownames(gene_counts) <- counts_annotated$gene_name
gene_counts <- gene_counts[rowSums(gene_counts) > 0,]
gene_tpm <- tpm_result[,grep("HM0", colnames(tpm_result))]
rownames(gene_tpm) <- tpm_result$gene_name
gene_tpm <- gene_tpm[rownames(gene_counts),]
saveRDS(gene_counts, "gene_counts.rds")
saveRDS(gene_tpm, "gene_tpm.rds")

#### ---- differential analysis for replicative senescence ---- ####
# conda activate SingleCellAnalysis
library(edgeR)
setwd("~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/RNAseq/replicative_senescence")
source("~/longzhilin/ProjectCode/AgingProject/edgeR_gene.R")
counts <- readRDS("gene_counts.rds")
condition <- data.frame(sample = colnames(counts), label = c("P11_1", "P11_2", "P24_1", "P24_2", "P39_1", "P39_2"))
condition$group <- gsub("_.*", "", condition$label)
condition$group <- factor(condition$group, levels = c("P11", "P24", "P39"))
compare.groups <- list(P39.vs.P24 = c("P39_1", "P39_2", "P24_1", "P24_2"), P24.vs.P11 = c("P24_1", "P24_2", "P11_1", "P11_2"), P39.vs.P11 = c("P39_1", "P39_2", "P11_1", "P11_2"))
dir.create("4.differentialGeneAnalysis", showWarnings = FALSE, recursive = TRUE)
pdf("4.differentialGeneAnalysis/edgeR.pdf", width = 6, height = 6)
group.DEGs <- lapply(compare.groups, function(label){
  cat("    ", label, "\n")
  condition.subset <- condition[which(condition$label %in% label),]
  counts.subset <- counts[,condition.subset$sample]
  cpms.subset <- edgeR::cpm(counts.subset)
  group <- droplevels(condition.subset$group)
  cat("group: ", group, "\n")
  out <- edgeR_gene_simple(counts = counts.subset, group = group, myContrasts = paste0(rev(levels(group)), collapse = "-"), cpm.filter = FALSE)
  p <- edgeR_plot(edgeR_result = out, PValue = FALSE, threshold = 0.05, top.n = 5, size = 3.5, title = paste0(rev(levels(group)), collapse = " vs "))
  print(p)
  write.csv(out, paste0("4.differentialGeneAnalysis/", paste0(rev(levels(group)), collapse = ".vs."), ".csv"))
  return(out)
})
dev.off()
saveRDS(group.DEGs, "4.differentialGeneAnalysis/group.DEGs.rds")

#### pathway enrichment
## GSEA method
source("~/Funtions/clusterProfiler.GSEA.R")
GSEA.res <- lapply(group.DEGs, function(x){
  geneRank <- sign(x$logFC) * sqrt(x$LR)
  names(geneRank) <- x$gene
  GSEA.res <- clusterProfiler.GSEA(geneRank, geneType = "SYMBOL", db.type = "GO", GO.ont = "BP", pvalueCutoff = 1, species = "Homo sapiens")
  up <- GSEA.res[which(GSEA.res$NES > 0), ]
  down <- GSEA.res[which(GSEA.res$NES < 0), ]
  return(list(up = up, down = down))
})
GSEA.res <- unlist(GSEA.res, recursive = FALSE)
pathwayEnrichment <- list(GSEA = GSEA.res)
pathwayEnrichment <- unlist(pathwayEnrichment, recursive = FALSE)
saveRDS(pathwayEnrichment, "4.differentialGeneAnalysis/pathwayEnrichment.rds")
writexl::write_xlsx(pathwayEnrichment, "4.differentialGeneAnalysis/pathwayEnrichment.xlsx")

selected.pathway.up <- c("regulation of immune effector process", "extracellular matrix organization", "negative regulation of cell population proliferation", "positive regulation of T cell mediated immunity", "positive regulation of lymphocyte mediated immunity", "positive regulation of leukocyte adhesion to vascular endothelial cell", "regulation of endothelial cell apoptotic process", "regulation of interleukin-6 production", "regulation of inflammatory response", "antigen processing and presentation of peptide antigen", "antigen processing and presentation via MHC class Ib", "positive regulation of leukocyte mediated cytotoxicity")
selected.pathway.down <- c("cytoplasmic translation", "ATP synthesis coupled electron transport", "oxidative phosphorylation", "mitotic sister chromatid segregation", "positive regulation of DNA metabolic process", "positive regulation of DNA repair", "positive regulation of telomere maintenance", "positive regulation of double-strand break repair", "DNA-templated DNA replication", "positive regulation of vascular endothelial cell proliferation", "ribosome biogenesis", "mitotic cell cycle phase transition")
source("~/Funtions/pathway.barplot.R")
pdf("4.differentialGeneAnalysis/enrichedPathways.barplot.pdf", width = 3, height = 4)
pathway.barplot(enrichedResult = pathwayEnrichment$GSEA.P39.vs.P11.up, selected.pathways = selected.pathway.up, label.value = "Description", rank.value = "p.adjust", title = "P39 vs P11: NES > 0", color = "#ec513c", threshold = 0.05)
pathway.barplot(enrichedResult = pathwayEnrichment$GSEA.P39.vs.P11.down, selected.pathways = selected.pathway.down, label.value = "Description", rank.value = "p.adjust", title = "P39 vs P11: NES < 0", color = "#035AA6", threshold = 0.05)

pathway.barplot(enrichedResult = pathwayEnrichment$GSEA.P39.vs.P24.up, selected.pathways = selected.pathway.up, label.value = "Description", rank.value = "p.adjust", title = "P39 vs P24: NES > 0", color = "#ec513c", threshold = 0.05)
pathway.barplot(enrichedResult = pathwayEnrichment$GSEA.P39.vs.P24.down, selected.pathways = unique(c(selected.pathway.down, "DNA strand elongation", "positive regulation of cell cycle")), label.value = "Description", rank.value = "p.adjust", title = "P39 vs P24: NES < 0", color = "#035AA6", threshold = 0.05)
dev.off()

#### ---- Gene expression analysis: P11, P24, P39 ---- ####
library(ComplexHeatmap)
library(circlize)
library(ggpubr)
setwd("~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/RNAseq/replicative_senescence")
group.DEGs <- readRDS("4.differentialGeneAnalysis/group.DEGs.rds")
HUVEC.tpm.mat <- readRDS("gene_tpm.rds")
colnames(HUVEC.tpm.mat) <- c("P11_1", "P11_2", "P24_1", "P24_2", "P39_1", "P39_2")

#### plot SASP genes
SASP.genes <- c("CDKN1A", "IL1A", "IL1B", "IL6", "MMP10", "CXCL2", "CCL2", "ICAM1", "AREG", "IGFBP2")
idx <- match(SASP.genes, rownames(HUVEC.tpm.mat))
plot.data <- log2(HUVEC.tpm.mat[na.omit(idx),] + 1)
plot.data <- plot.data[!duplicated(rownames(plot.data)), , drop = FALSE]
p1 <- Heatmap(t(plot.data), name = "Log2(TPM+1)", height = unit(3, "cm"), width = unit(6, "cm"), , cluster_columns = F, show_row_dend = F, cluster_rows = F, show_column_dend = F)

plot.data <- as.matrix(plot.data)
row_norm <- sqrt(rowSums(plot.data^2, na.rm = TRUE))
row_norm[row_norm == 0] <- 1  # 全零行保持为零
plot.data_l2 <- sweep(plot.data, 1, row_norm, FUN = "/")

p11 <- Heatmap(t(plot.data_l2), name = "L2-Normalized expression", height = unit(3, "cm"), width = unit(6, "cm"), cluster_columns = F, show_row_dend = F, cluster_rows = F, show_column_dend = F, heatmap_legend_param = list(title = "L2-Normalized expression"))
pdf("HUVEC.SASP.genes.pdf", width = 6, height = 6)
draw(p1)
draw(p11)
dev.off()

#### ---- differential analysis for GSK126_treatment ---- ####
# conda activate SingleCellAnalysis
library(edgeR)
setwd("~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/RNAseq/GSK126_treatment")
source("~/Funtions/edgeR_gene.R")
counts <- readRDS("gene_counts.rds")
condition <- data.frame(sample = colnames(counts), label = c("DMSO_1", "DMSO_2", "GSK126_1", "GSK126_2", "GSK126_3"))
condition$group <- gsub("_.*", "", condition$label)
condition$group <- factor(condition$group, levels = c("DMSO", "GSK126"))
compare.groups <- list(GSK126.vs.DMSO = c("GSK126_1", "GSK126_2", "GSK126_3", "DMSO_1", "DMSO_2"))
dir.create("4.differentialGeneAnalysis", showWarnings = FALSE, recursive = TRUE)
pdf("4.differentialGeneAnalysis/edgeR.pdf", width = 6, height = 6)
group.DEGs <- lapply(compare.groups, function(label){
  cat("    ", label, "\n")
  condition.subset <- condition[which(condition$label %in% label),]
  counts.subset <- counts[,condition.subset$sample]
  cpms.subset <- edgeR::cpm(counts.subset)
  group <- droplevels(condition.subset$group)
  cat("group: ", group, "\n")
  out <- edgeR_gene_simple(counts = counts.subset, group = group, myContrasts = paste0(rev(levels(group)), collapse = "-"), cpm.filter = FALSE)
  p <- edgeR_plot(edgeR_result = out, PValue = FALSE, threshold = 0.05, top.n = 5, size = 3.5, title = paste0(rev(levels(group)), collapse = " vs "))
  print(p)
  write.csv(out, paste0("4.differentialGeneAnalysis/", paste0(rev(levels(group)), collapse = ".vs."), ".csv"))
  return(out)
})
dev.off()
saveRDS(group.DEGs, "4.differentialGeneAnalysis/group.DEGs.rds")

#### pathway enrichment
## GSEA method
source("~/Funtions/clusterProfiler.GSEA.R")
GSEA.res <- lapply(group.DEGs, function(x){
  geneRank <- sign(x$logFC) * sqrt(x$LR)
  names(geneRank) <- x$gene
  GSEA.res <- clusterProfiler.GSEA(geneRank, geneType = "SYMBOL", db.type = "GO", GO.ont = "BP", pvalueCutoff = 1, species = "Homo sapiens")
  up <- GSEA.res[which(GSEA.res$NES > 0), ]
  down <- GSEA.res[which(GSEA.res$NES < 0), ]
  return(list(up = up, down = down))
})
GSEA.res <- unlist(GSEA.res, recursive = FALSE)
pathwayEnrichment <- list(GSEA = GSEA.res)
pathwayEnrichment <- unlist(pathwayEnrichment, recursive = FALSE)
saveRDS(pathwayEnrichment, "4.differentialGeneAnalysis/pathwayEnrichment.rds")
writexl::write_xlsx(pathwayEnrichment, "4.differentialGeneAnalysis/pathwayEnrichment.xlsx")

## plot enriched pathways
source("~/Funtions/pathway.dotplot.R")
p <- pathway.dotplot(enrichedResult.lists = pathwayEnrichment[3:4], title = "", items = c("Description", "GeneRatio", "p.adjust"), threshold = 0.05, top.n = 15, GeneRatio = T, max.dot.size = 4, sub.words = 80, y.text.size = 11)
print(p)
p <- pathway.dotplot(enrichedResult.lists = pathwayEnrichment[1:2], title = "", items = c("Description", "NES", "p.adjust"), threshold = 0.05, top.n = 15, GeneRatio = F, max.dot.size = 4, sub.words = 80, y.text.size = 11)
print(p)

source("~/Funtions/pathway.barplot.R")
pdf("4.differentialGeneAnalysis/enrichedPathways.barplot.pdf", width = 3, height = 4)
pathway.barplot(enrichedResult = pathwayEnrichment$GSEA.GSK126.vs.DMSO.up, selected.pathways = unique(c(selected.pathway.up, "neuron projection guidance", "circulatory system process", "cellular response to BMP stimulus", "positive regulation of apoptotic process", "tissue morphogenesis", "negative regulation of endothelial cell proliferation", "antigen processing and presentation via MHC class Ib")), label.value = "Description", rank.value = "p.adjust", title = "GSK126 vs DMSO: NES > 0", color = "#ec513c", threshold = 0.05)
pathway.barplot(enrichedResult = pathwayEnrichment$GSEA.GSK126.vs.DMSO.down, selected.pathways = unique(c(selected.pathway.up, "mRNA processing", "negative regulation of telomere maintenance", "DNA-templated DNA replication", "regulation of mitotic cell cycle phase transition", "regulation of endothelial cell development", "RNA export from nucleus", "facultative heterochromatin formation", "double-strand break repair", "DNA replication", "positive regulation of mRNA catabolic process", "regulation of chromatin organization", "regulation of endothelial cell development")), label.value = "Description", rank.value = "p.adjust", title = "GSK126 vs DMSO: NES < 0", color = "#035AA6", threshold = 0.05)
dev.off()

#### ---- GSK126_treatment vs replicative_senescence ---- ####
setwd("~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/RNAseq")
GSK126_treatment.tpm <- readRDS("GSK126_treatment/gene_tpm.rds")
colnames(GSK126_treatment.tpm) <- c("DMSO_1", "DMSO_2", "GSK126_1", "GSK126_2", "GSK126_3")
replicative_senescence.tpm <- readRDS("replicative_senescence/gene_tpm.rds")
colnames(replicative_senescence.tpm) <- c("P11_1", "P11_2", "P24_1", "P24_2", "P39_1", "P39_2")
common.genes <- intersect(rownames(GSK126_treatment.tpm), rownames(replicative_senescence.tpm))
GSK126_treatment.tpm <- GSK126_treatment.tpm[common.genes,]
replicative_senescence.tpm <- replicative_senescence.tpm[common.genes,]
tpm.mat <- cbind(replicative_senescence.tpm, GSK126_treatment.tpm)

#### hypergeometric-test
SASP.signature <- read.table("~/longzhilin/AnalysisCode/Data/signatureGeneSet/senescence/cellularSenescenceSignature.Zhang.2022.CellRes.txt", header = T, stringsAsFactors = F, sep = "\t")
up_overlap <- list(P39_vs_P24 = P39_vs_P24.up, GSK126_vs_DMSO = plot.data$gene[which(plot.data$logFC > 1 & plot.data$fdr < 0.05)])
up_P39 <- up_overlap$P39_vs_P24
up_GSK126 <- up_overlap$GSK126_vs_DMSO

library(hypeR)
mhyp <- hypeR(up_GSK126, list(P39_vs_P24 = up_P39), test="hypergeometric")
pval_txt <- format(mhyp$data$pval, scientific = TRUE, digits = 2)

pdf("GSK126_DMSO.P39_vs_P24.up.pdf", width = 4, height = 4)
ggvenn::ggvenn(up_overlap, show_percentage = F, auto_scale = FALSE, set_name_size = 5, text_size = 5, fill_color = c("#039E73", "#d30808")) + ggtitle(bquote("Hypergeometric test: " ~ italic(p) == .(pval_txt)))
mhyp <- hypeR(intersect(up_P39, up_GSK126), list(SASP = SASP.signature$Gene), test="hypergeometric")
pval_txt <- format(mhyp$data$pval, scientific = TRUE, digits = 2)
ggvenn::ggvenn(list(Overlap = intersect(up_P39, up_GSK126), SASP = SASP.signature$Gene), show_percentage = F, auto_scale = FALSE, set_name_size = 5, text_size = 5, fill_color = c("#039E73", "#d30808")) + ggtitle(bquote("Hypergeometric test: " ~ italic(p) == .(pval_txt)))
dev.off()

#### Correlation of logFC
library(ggplot2)
library(ggrepel)
SASP.genes <- c("CDKN1A", "IL1A", "IL1B", "IL6", "MMP10", "CXCL2", "CCL2", "ICAM1", "AREG", "IGFBP2")
genes <- unique(na.omit(c(point_labels$gene, SASP.genes)))
HUVEC_GSK <- GSK.DEGs$GSK126.vs.DMSO[, c("gene", "fdr", "logFC")]
HUVEC_GSK <- HUVEC_GSK[!duplicated(HUVEC_GSK$gene), ]
colnames(HUVEC_GSK)[2:3] <- c("GSK126_vs_DMSO.fdr", "GSK126_vs_DMSO.logFC")
pdf("logFC.correlation.pdf", width = 6, height = 6)
logFC.plot <- lapply(names(HUVEC.DEGs), function(group) {
  HUVEC.res <- HUVEC.DEGs[[group]][, c("gene", "fdr", "logFC")]
  HUVEC.res <- HUVEC.res[!duplicated(HUVEC.res$gene), ]
  colnames(HUVEC.res)[2:3] <- paste0(group, c(".fdr", ".logFC"))
  plot.data <- merge(HUVEC.res, HUVEC_GSK, by = "gene")
  
  lfc_x <- plot.data[["GSK126_vs_DMSO.logFC"]]
  lfc_y <- plot.data[[paste0(group, ".logFC")]]
  fdr_x <- plot.data[["GSK126_vs_DMSO.fdr"]]
  fdr_y <- plot.data[[paste0(group, ".fdr")]]

  keep <- complete.cases(plot.data[, 2:5]) & is.finite(lfc_x) & is.finite(lfc_y) & is.finite(fdr_x) & is.finite(fdr_y)
  plot.data <- plot.data[keep, ]
  lfc_x <- plot.data[["GSK126_vs_DMSO.logFC"]]
  lfc_y <- plot.data[[paste0(group, ".logFC")]]
  fdr_x <- plot.data[["GSK126_vs_DMSO.fdr"]]
  fdr_y <- plot.data[[paste0(group, ".fdr")]]

  plot.data$lfc_x <- lfc_x
  plot.data$lfc_y <- lfc_y

  x_status <- ifelse(
    fdr_x < 0.05 & lfc_x > 1, "Up",
    ifelse(fdr_x < 0.05 & lfc_x < -1, "Down", "NS")
  )

  y_status <- ifelse(
    fdr_y < 0.05 & lfc_y > 1, "Up",
    ifelse(fdr_y < 0.05 & lfc_y < -1, "Down", "NS")
  )
  status_levels <- c("Up", "Down", "NS")
  x_status <- factor(x_status, levels = status_levels)
  y_status <- factor(y_status, levels = status_levels)

  plot.data$quadrant <- factor(
    ifelse(x_status == "Up" & y_status == "Up", "Q1",
      ifelse(x_status == "Down" & y_status == "Up", "Q2",
        ifelse(x_status == "Down" & y_status == "Down", "Q3",
          ifelse(x_status == "Up" & y_status == "Down", "Q4", "NS")
        )
      )
    ),
    levels = c("Q1", "Q2", "Q3", "Q4", "NS")
  )

  status_table <- table(
    factor(x_status, levels = status_levels),
    factor(y_status, levels = status_levels)
  )

  chisq_result <- suppressWarnings(chisq.test(status_table))
  chisq_p <- chisq_result$p.value

  subtitle_text <- paste0(
    "3-state association: Chi-squared p = ",
    format.pval(chisq_p, digits = 3, eps = 1e-3)
  )

  quadrant_count <- table(
    factor(
      plot.data$quadrant,
      levels = c("Q1", "Q2", "Q3", "Q4", "NS")
    )
  )

  q_count <- as.integer(quadrant_count[c("Q1", "Q2", "Q3", "Q4")])
  sig.genes <- plot.data[
    plot.data$gene %in% genes & plot.data$quadrant != "NS",
    ,
    drop = FALSE
  ]

  xlim <- range(lfc_x, na.rm = TRUE)
  ylim <- range(lfc_y, na.rm = TRUE)
  if (diff(xlim) == 0) xlim <- xlim + c(-1, 1)
  if (diff(ylim) == 0) ylim <- ylim + c(-1, 1)

  label_pos <- data.frame(
    quadrant = factor(
      c("Q1", "Q2", "Q3", "Q4"),
      levels = c("Q1", "Q2", "Q3", "Q4")
    ),
    x = c(xlim[2], xlim[1], xlim[1], xlim[2]) * 0.9,
    y = c(ylim[2], ylim[2], ylim[1], ylim[1]) * 0.9,
    n = q_count
  )

  label_pos$label <- paste0(
    label_pos$quadrant,
    "\nn = ",
    label_pos$n
  )

  cor_data <- plot.data[is.finite(plot.data$lfc_x) & is.finite(plot.data$lfc_y),]
  cor_res <- cor.test(cor_data$lfc_x, cor_data$lfc_y, method = "spearman")
  cor_text <- sprintf(
    "Spearman r = %.2f; p = %s; n = %d",
    unname(cor_res$estimate),
    format.pval(cor_res$p.value, digits = 2, eps = 1e-16),
    nrow(cor_data)
  )

  cols <- c(Q1 = "#E64B35", Q2 = "#4DBBD5", Q3 = "#00A087", Q4 = "#3C5488")
  p_plot <- ggplot(
    plot.data,
    aes(x = lfc_x, y = lfc_y)
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey60"
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey60"
    ) +
    geom_smooth(
      data = cor_data,
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "blue",
      linewidth = 0.8
    ) +
    geom_point(
      data = subset(plot.data, quadrant == "NS"),
      color = "grey80",
      size = 1.2,
      alpha = 0.6
    ) +
    geom_point(
      data = subset(plot.data, quadrant != "NS"),
      aes(color = quadrant),
      size = 1.6,
      alpha = 0.85
    ) +
    scale_color_manual(
      values = cols,
      name = NULL
    ) +
    geom_text(
      data = label_pos,
      aes(
        x = x,
        y = y,
        label = label,
        color = quadrant
      ),
      fontface = "bold",
      size = 3.2,
      show.legend = FALSE
    ) +
    labs(
      x = "GSK126 vs DMSO logFC",
      y = paste0(group, " logFC"),
      title = paste0("logFC correlation: GSK126 vs ", group),
      subtitle = paste(subtitle_text, cor_text, sep = "\n")
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      aspect.ratio = 1,
      legend.position = "none"
    )

  if (nrow(sig.genes) > 0) {
    p_plot <- p_plot +
      geom_text_repel(
        data = sig.genes,
        aes(x = lfc_x, y = lfc_y, label = gene),
        color = "black",
        size = 3,
        segment.size = 0.8,
        min.segment.length = 0,
        segment.color = "grey40",
        max.overlaps = Inf,
        show.legend = FALSE
      )
  }
  print(p_plot)
  return(p_plot)
})
names(logFC.plot) <- names(HUVEC.DEGs)
dev.off()
