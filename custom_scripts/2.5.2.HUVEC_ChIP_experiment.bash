#!/bin/bash
#' @description CUT&Tag for IMR90 H3K9me3
#' @Author: Zhilin Long
#' @date 2026-09-08

#### ---- Data processing ---- ####
# conda activate BigDataProcessing
#### H3K27me3
cd ~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/H3K27me3
nohup bash run_chipseq.sh > run_chipseq.log 2>&1 &

#### H3K9me3
cd ~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/H3K9me3
nohup bash run_chipseq.sh > run_chipseq.log 2>&1 &

#### ---- re-smooth bw ---- ####
cd ~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/H3K27me3
bamCoverage -p 16 --binSize 50 --smoothLength 1000 -b "bam/DYQ302-1_L03.nodup.bam" -o "bigWig/DYQ302-1_L03.nodup.smooth.bw" --outFileFormat bigwig --normalizeUsing RPKM
bamCoverage -p 16 --binSize 50 --smoothLength 1000 -b "bam/DYQ302-2_L03.nodup.bam" -o "bigWig/DYQ302-2_L03.nodup.smooth.bw" --outFileFormat bigwig --normalizeUsing RPKM
bamCoverage -p 16 --binSize 50 --smoothLength 1000 -b "bam/DYQ302-5_L03.nodup.bam" -o "bigWig/DYQ302-5_L03.nodup.smooth.bw" --outFileFormat bigwig --normalizeUsing RPKM
bamCoverage -p 16 --binSize 50 --smoothLength 1000 -b "bam/DYQ302-6_L03.nodup.bam" -o "bigWig/DYQ302-6_L03.nodup.smooth.bw" --outFileFormat bigwig --normalizeUsing RPKM

#### ---- call peak for P11 ---- ####
## H3K9me3
export PYTHONPATH="~/longzhilin/Softwares/Epigenomics/SICER2:${PYTHONPATH:-}"
sicer="~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer"
workPath="~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/H3K9me3"
peak_dir="${workPath}/peak_sicer2"
mkdir -p "${peak_dir}"
cd "${peak_dir}"
# P11
merged_bam="$peak_dir/P11_H3K9me3.bam"
samtools merge -@ 12 -f "$merged_bam" "${workPath}/bam/DYQ303-1_L04.nodup.bam" "${workPath}/bam/DYQ303-2_L04.nodup.bam"
samtools index -@ 12 "$merged_bam"
~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t "$merged_bam" -s hg38 -w 5000 -g 10000 -cpu 16

# P24
sicer="~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer"
merged_bam="$peak_dir/P24_H3K9me3.bam"
samtools merge -@ 12 -f "$merged_bam" "${workPath}/bam/DYQ303-3_L04.nodup.bam" "${workPath}/bam/DYQ303-4_L04.nodup.bam"
samtools index -@ 12 "$merged_bam"
~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t "$merged_bam" -s hg38 -w 5000 -g 10000 -cpu 16

## H3K27me3
export PYTHONPATH="~/longzhilin/Softwares/Epigenomics/SICER2:${PYTHONPATH:-}"
sicer="~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer"
workPath="~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/H3K27me3"
peak_dir="${workPath}/peak_sicer2"
mkdir -p "${peak_dir}"
cd "${peak_dir}"

# P11
merged_bam="$peak_dir/P11_H3K27me3.bam"
samtools merge -@ 12 -f "$merged_bam" "${workPath}/bam/DYQ302-1_L03.nodup.bam" "${workPath}/bam/DYQ302-2_L03.nodup.bam"
samtools index -@ 12 "$merged_bam"
~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t "$merged_bam" -s hg38 -w 5000 -g 10000 -cpu 16

# P24
merged_bam="$peak_dir/P24_H3K27me3.bam"
samtools merge -@ 12 -f "$merged_bam" "${workPath}/bam/DYQ302-3_L03.nodup.bam" "${workPath}/bam/DYQ302-4_L03.nodup.bam"
samtools index -@ 12 "$merged_bam"
~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t "$merged_bam" -s hg38 -w 5000 -g 10000 -cpu 16

#### ---- quantify the 10kb bins ---- ####
## make 10kb bins
cd ~/longzhilin/Softwares/genomes/human/hg38
genome_sizes="~/longzhilin/Softwares/genomes/human/hg38/hg38.chrom.sizes"
bedtools  makewindows -g $genome_sizes  -w 10000 > hg38_10kb_bins.bed;
awk -v OFS='\t' '{print "bin" NR, $0}' hg38_10kb_bins.bed > hg38_10kb_bins.tmp.saf;
awk 'BEGIN{FS=OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"}{print $1, $2, $3+1, $4, "."}' hg38_10kb_bins.tmp.saf > hg38_10kb_bins.saf;
rm hg38_10kb_bins.tmp.saf;

## quantify
workPath="~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/H3K27me3"
bam_files=$(find $workPath -name "*.nodup.bam" | sort -V)
mkdir -p "${workPath}/bin10kb"
cd "${workPath}/bin10kb"
featureCounts -T 16 -p -M --fraction -F SAF \
    -a ~/longzhilin/Softwares/genomes/human/hg38/hg38_10kb_bins.saf \
    -o counts_10kb_bins.txt ${bam_files}

workPath="~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/H3K9me3"
bam_files=$(find $workPath -name "*.nodup.bam" | sort -V)
mkdir -p "${workPath}/bin10kb"
cd "${workPath}/bin10kb"
featureCounts -T 16 -p -M --fraction -F SAF \
    -a ~/longzhilin/Softwares/genomes/human/hg38/hg38_10kb_bins.saf \
    -o counts_10kb_bins.txt ${bam_files}

## differential bins (R code)
setwd("~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment")
source("~/longzhilin/ProjectCode/AgingProject/edgeR_diffRegions.R")
histones <- c("H3K9me3", "H3K27me3")
sampleNames <- c("P11_1", "P11_2", "P24_1", "P24_2", "P39_1", "P39_2")

diffBin.list <- lapply(histones, function(histone){
  bin10kb <- read.table(paste0(histone, "/bin10kb/counts_10kb_bins.txt"), header = T, sep = "\t", stringsAsFactors = F)
  peak <- paste0(bin10kb$Chr, ":", bin10kb$Start, "-", bin10kb$End)
  bin10kb.count <- bin10kb[,7:12]
  colnames(bin10kb.count) <- sampleNames
  rownames(bin10kb.count) <- peak
  count.subset <- bin10kb.count[, c("P24_1", "P24_2", "P39_1", "P39_2")]
  group <- factor(c(rep("P24", 2), rep("P39", 2)), levels = c("P24", "P39"))
  diffBins <- edgeR.diffRegions(countTables = count.subset, group = group, outPath = ".", type = "P39_vs_P24", title = NULL, logFC = 1, PValue = 0.05, ContinuousComparison = FALSE, cpm.filter = TRUE, min.regions = 10)
  diffBins <- diffBins$GroupComparison.P39.vs.P24
  write.csv(diffBins, file = paste0(histone, "_diffBins.csv"), row.names = F)
  return(diffBins)
})
names(diffBin.list) <- histones
saveRDS(diffBin.list, "diffBin.list.P39.vs.P24.rds")

#### ---- MA plot
library(ggplot2)
library(dplyr)
plot_data <- diffBin.list$H3K27me3 %>%
  mutate(
    group = case_when(
      !is.na(fdr) & fdr < 0.05 & logFC > 0 ~ "Upregulated",
      !is.na(fdr) & fdr < 0.05 & logFC < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )
counts <- table(factor(
  plot_data$group,
  levels = c("Upregulated", "Downregulated", "Not significant")
))
cat("显著上调数量：", as.integer(counts["Upregulated"]), "\n")
cat("显著下调数量：", as.integer(counts["Downregulated"]), "\n")
ggplot(plot_data, aes(x = logCPM, y = logFC, color = group)) +
  geom_point(size = 1.5, alpha = 0.75) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  annotate(
    "label",
    x = Inf, y = Inf,
    label = paste0("Upregulated: ", as.integer(counts["Upregulated"])),
    hjust = 1.05, vjust = 1.2,
    color = "#E64B35",
    fill = "white",
    label.size = NA
  ) +
  annotate(
    "label",
    x = Inf, y = Inf,
    label = paste0("Downregulated: ", as.integer(counts["Downregulated"])),
    hjust = 1.05, vjust = 2.2,
    color = "#3C5488",
    fill = "white",
    label.size = NA
  ) +
  scale_color_manual(
    values = c(
      "Upregulated" = "#E64B35",
      "Downregulated" = "#3C5488",
      "Not significant" = "grey75"
    )
  ) +
  labs(
    x = "logCPM",
    y = "logFC",
    color = NULL
  ) +
  theme_classic()

#### ---- bin overlapping peak ---- ####
library(GenomicRanges)
library(ggpubr)
diffBin.list <- readRDS("diffBin.list.P39.vs.P24.rds")
#### H3K27me3
H3K27me3.peak <- read.table("~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/H3K27me3/peak_sicer2/P24_H3K27me3-W5000-G10000.scoreisland", header = F, stringsAsFactors = F, sep = "\t")
write.table(H3K27me3.peak[,1:3], "~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/H3K27me3/peak_sicer2/P24_H3K27me3-W5000-G10000.bed", col.names = F, row.names = F, sep = "\t", quote = F)
H3K27me3.peak <- makeGRangesFromDataFrame(H3K27me3.peak, seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = T)
x <- diffBin.list$H3K27me3
region <- x$Region
region <- data.frame(Chr = gsub(":.*", "", region), Start = gsub(".*:", "", gsub("-.*", "", region)), End = gsub(".*-", "", region))
region <- makeGRangesFromDataFrame(region, seqnames.field = "Chr", start.field = "Start", end.field = "End")
idx <- findOverlaps(region, H3K27me3.peak)
region$type <- "Not"
region$type[unique(queryHits(idx))] <- "overlappingPeak"
region.up <- region[which(x$fdr < 0.05 & x$logFC > 0),]
region.down <- region[which(x$fdr < 0.05 & x$logFC < 0),]
H3K27me3.diffBin.sig <- list(up = region.up, down = region.down)
saveRDS(H3K27me3.diffBin.sig, "H3K27me3/H3K27me3.P39.vs.P24.diffBin.sig.rds")

#### H3K9me3
H3K9me3.peak <- read.table("~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/H3K9me3/peak_sicer2/P24_H3K9me3-W5000-G10000.scoreisland", header = F, stringsAsFactors = F, sep = "\t")
write.table(H3K9me3.peak[,1:3], "~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/H3K9me3/peak_sicer2/P24_H3K9me3-W5000-G10000.bed", col.names = F, row.names = F, sep = "\t", quote = F)
H3K9me3.peak <- makeGRangesFromDataFrame(H3K9me3.peak, seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = T)
x <- diffBin.list$H3K9me3
region <- x$Region
region <- data.frame(Chr = gsub(":.*", "", region), Start = gsub(".*:", "", gsub("-.*", "", region)), End = gsub(".*-", "", region))
region <- makeGRangesFromDataFrame(region, seqnames.field = "Chr", start.field = "Start", end.field = "End")
idx <- findOverlaps(region, H3K9me3.peak)
region$type <- "Not"
region$type[unique(queryHits(idx))] <- "overlappingPeak"
region.up <- region[which(x$fdr < 0.05 & x$logFC > 0),]
region.down <- region[which(x$fdr < 0.05 & x$logFC < 0),]
H3K9me3.diffBin.sig <- list(up = region.up, down = region.down)
saveRDS(H3K9me3.diffBin.sig, "H3K9me3/H3K9me3.P39.vs.P24.diffBin.sig.rds")

#### plot number
source("~/longzhilin/ProjectCode/AgingProject/plotDiffOverlap.R")
diffBin.sig.number <- lapply(names(H3K9me3.diffBin.sig), function(type){
y <- H3K9me3.diffBin.sig[[type]]
total.bins <- length(y)
peak.bins <- length(y[which(y$type == "overlappingPeak")])
data.frame(group = "H3K9me3_P39_vs_P24", label = type, type = c("overlappingPeak", "totalPeak"), number = c(peak.bins, total.bins))
})
diffBin.sig.number <- do.call(rbind, diffBin.sig.number)

# plot split up/down
# Fisher's exact test: down overlap proportion > up
up_overlap_count <- diffBin.sig.number$number[which(diffBin.sig.number$label == "up" & diffBin.sig.number$type == "overlappingPeak")]
up_total <- diffBin.sig.number$number[which(diffBin.sig.number$label == "up" & diffBin.sig.number$type == "totalPeak")]
down_overlap_count <- diffBin.sig.number$number[which(diffBin.sig.number$label == "down" & diffBin.sig.number$type == "overlappingPeak")]
down_total <- diffBin.sig.number$number[which(diffBin.sig.number$label == "down" & diffBin.sig.number$type == "totalPeak")]
overlap_mat <- matrix(c(down_overlap_count, down_total - down_overlap_count,
                        up_overlap_count, up_total - up_overlap_count),
                      nrow = 2, byrow = TRUE,
                      dimnames = list(c("down", "up"), c("Overlap", "NonOverlap")))
fisher_res <- fisher.test(overlap_mat, alternative = "greater")
cat("Fisher's exact test (down > up) p-value:", fisher_res$p.value, "\n")
cat("Odds ratio:", unname(fisher_res$estimate), "\n")
cat("95% CI:", fisher_res$conf.int, "\n")

diffBin.sig.number$group_label <- paste0(diffBin.sig.number$group, "_", diffBin.sig.number$label)
diffBin.sig.number$group_label <- factor(diffBin.sig.number$group_label, levels = c("H3K9me3_P39_vs_P24_up", "H3K9me3_P39_vs_P24_down"))
diffBin.sig.number$type <- factor(diffBin.sig.number$type, levels = c("totalPeak", "overlappingPeak"))
diffBin.sig.number <- diffBin.sig.number[order(diffBin.sig.number$group_label, 
                                                diffBin.sig.number$type), ]
p11 <- ggplot(diffBin.sig.number, aes(x = group_label, y = number, fill = label, alpha = type)) +
  geom_col(position = "identity", width = 0.7, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  coord_flip() +
  scale_fill_manual(values = c(up = "#d30808", down = "#0859FD"), name = "Direction") +
  scale_alpha_manual(values = c(totalPeak = 0.4, overlappingPeak = 0.8), name = "Type") +
  labs(title = "H3K9me3", x = "", y = "Number of differential bins", subtitle = paste0("Fisher p(down>up) = ", signif(fisher_res$p.value, 3))) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        aspect.ratio = 0.5,
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"))

## plot number: diverging barplot (up positive/red on right, down negative/blue on left)
diffBin.sig.number$number[diffBin.sig.number$label == "down"] <- -diffBin.sig.number$number[diffBin.sig.number$label == "down"]
diffBin.sig.number$label <- factor(diffBin.sig.number$label, levels = c("up", "down"))
diffBin.sig.number$type <- factor(diffBin.sig.number$type, levels = c("totalPeak", "overlappingPeak"))
diffBin.sig.number <- diffBin.sig.number[order(diffBin.sig.number$group, diffBin.sig.number$label, diffBin.sig.number$type), ]
p1 <- ggplot(diffBin.sig.number, aes(x = group, y = number, fill = label, alpha = type)) +
  geom_col(position = "identity", width = 0.7, color = NA) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  coord_flip() +
  scale_fill_manual(values = c(up = "#d30808", down = "#0859FD"), name = "Direction") +
  scale_alpha_manual(values = c(overlappingPeak = 0.8, totalPeak = 0.4), name = "Type") +
  labs(title = "H3K9me3", x = "", y = "Differential bin numbers") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        aspect.ratio = 0.5,
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"))

pdf("H3K9me3/H3K9me3_binOverlapPeak2.pdf", width = 7, height = 5)
ggarrange(p1, p11, ncol = 1, nrow = 3, align = "hv")
dev.off()

#### ---- bw compare ---- ####
set -euo pipefail
cd ~/longzhilin/Projects/AgingMouse/PairedTag/5.otherData/HUVEC_experiment/H3K27me3/bigWig

# P39 vs P11
bigwigCompare -b1 DYQ302-1_L03.nodup.bw -b2 DYQ302-2_L03.nodup.bw \
  --operation add --binSize 50 --numberOfProcessors 8 \
  -o "P11.bw"

bigwigCompare -b1 DYQ302-5_L03.nodup.bw -b2 DYQ302-6_L03.nodup.bw \
  --operation add --binSize 50 --numberOfProcessors 8 \
  -o "P39.bw"

bigwigCompare -b1 "P39.bw" -b2 "P11.bw" \
  --operation log2 --pseudocount 1 1 --skipZeroOverZero \
  --binSize 50 --numberOfProcessors 8 \
  -o P39.vs.P11.log2ratio.nodup.bw

# P39 vs P24
bigwigCompare -b1 DYQ302-3_L03.nodup.bw -b2 DYQ302-4_L03.nodup.bw \
  --operation add --binSize 50 --numberOfProcessors 8 \
  -o "P24.bw"

bigwigCompare -b1 "P39.bw" -b2 "P24.bw" \
  --operation log2 --pseudocount 1 1 --skipZeroOverZero \
  --binSize 50 --numberOfProcessors 8 \
  -o P39.vs.P24.log2ratio.nodup.bw

# P24 vs P11
bigwigCompare -b1 "P24.bw" -b2 "P11.bw" \
  --operation log2 --pseudocount 1 1 --skipZeroOverZero \
  --binSize 50 --numberOfProcessors 8 \
  -o P24.vs.P11.log2ratio.nodup.bw
