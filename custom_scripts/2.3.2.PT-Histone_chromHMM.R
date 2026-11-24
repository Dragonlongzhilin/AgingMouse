#!/bin/bash
#' @description running ChromHMM method for each cell type in Paired-Tag data
#' @author Zhilin Long
#' @date 2026-3-6
# conda activate PairedTagAnalysis

library(ggpubr)
library(ComplexHeatmap)
library(circlize)
set.seed(101)
workPath <- "~/longzhilin/Projects/AgingMouse/PairedTag/4.celltypeLevel/ChromHMM"
setwd(workPath)

#### read depths
depths.files <- list.files(file.path(workPath, "State15"), "input.rds", full.names = T, recursive = T)
depths <- lapply(depths.files, function(depths.file){
  depths <- readRDS(depths.file)
  depths$tissue <- gsub(":.*", "", depths$sample)
  return(depths[,c("tissue", "cellType", "age", "sample", "histone", "reads", "file")])
})
depths <- do.call("rbind", depths)
saveRDS(depths, "DNA.depths.rds")

#### ---- emissions 15 ---- ####
emissions_15 <- read.table("State15/MergeAllTissue/learnModel/emissions_15.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)
emissions_15 <- emissions_15[,c("H3K27me3", "H3K9me3", "H3K36me3", "H3K27ac", "H3K4me1", "H3K4me3")]
rownames(emissions_15) <- paste0("E", rownames(emissions_15))
# col_fun <- circlize::colorRamp2(
#   c(0, 0.25, 0.5, 0.75, 1),
#   c("white", "#FCDE6B", "#F89E48", "#F36A2D", "#E31A1C")
# )
col_fun <- circlize::colorRamp2(
  c(0, 1),
  c("white", "#E31A1C")
)
p1 <- Heatmap(as.matrix(emissions_15), col = col_fun, name = "Emission\nprobability", show_row_names = T, show_column_names = T, cluster_columns = T, cluster_rows = T, height = unit(6, "cm"), width = unit(3, "cm"))

state_label <- c(
  E1 = "Active (Pr-A)",
  E2 = "Upstream (Pr-Up)",
  E3 = "Bivalent (Pr-B)",
  E4 = "Strong (Tr-S)",
  E5 = "Permissive (Tr-P)",
  E6 = "Initiation (Tr-I)",
  E7 = "Strong, TSS-dist(En-Sd)",
  E8 = "Strong, TSS-prox(En-Sp)",
  E9 = "Weak (En-W)",
  E10 = "Poised (En-P)",
  E11 = "Polycomb (Hc-P)",
  E12 = "Polycomb weak (Hc-Pw)",
  E13 = "H3K9me3 (Hc-H)",
  E14 = "H3K9me3 weak (Hc-Hw)",
  E15 = "No signal (Ns)"
)
mat <- as.matrix(emissions_15)
mat <- mat[names(state_label),]
row_labels <- paste(rownames(mat), state_label[rownames(mat)])
p2 <- Heatmap(mat, col = col_fun, row_labels = row_labels, name = "Emission\nprobability", rect_gp = grid::gpar(col = "black", lwd = 1), border = TRUE, cluster_columns = F, cluster_rows = F, border_gp = grid::gpar(col = "black", lwd = 1), height = unit(6, "cm"), width = unit(3, "cm"))
p2_scale <- Heatmap(t(scale(t(mat))), col = col_fun, row_labels = row_labels, name = "Emission\nprobability", rect_gp = grid::gpar(col = "black", lwd = 1), border = TRUE, cluster_columns = F, cluster_rows = F, border_gp = grid::gpar(col = "black", lwd = 1), height = unit(6, "cm"), width = unit(3, "cm"))

#### proportions of each state in whole genome
library(data.table)
library(ggplot2)
model_dir <- "~/longzhilin/Projects/AgingMouse/PairedTag/4.celltypeLevel/ChromHMM/State15/MergeAllTissue"
result_dir <- file.path(model_dir, "VisualizationResults_20261014")
manifest_file <- file.path(model_dir, "VisualizationRedesign_20261014/sample_mark_manifest.tsv")
segment_dir <- file.path(model_dir, "learnModel")
chromsizes_file <- "~/longzhilin/Softwares/Epigenomics/ChromHMM/CHROMSIZES/mm10.txt"

# Read all sequence lengths in mm10.txt and use their sum as the fixed denominator.
# Do not filter chromosomes or exclude gap or blacklist regions.
chromsizes <- fread(chromsizes_file, header = FALSE, col.names = c("chrom", "length_bp"))
if(nrow(chromsizes) == 0L || anyNA(chromsizes) || any(chromsizes$length_bp <= 0) || anyDuplicated(chromsizes$chrom)) {
  stop("Invalid mm10 chromosome-size file: ", chromsizes_file)
}
reference_genome_bp <- sum(as.numeric(chromsizes$length_bp))
reference_contigs <- chromsizes$chrom

# Read sample metadata.
manifest <- unique(fread(manifest_file, select = c("sample", "tissue", "celltype", "age", "composite_rep")))
if (nrow(manifest) == 0L || anyDuplicated(manifest$sample)) {
  stop("Sample manifest is empty or contains duplicated samples: ", manifest_file)
}
setorder(manifest, sample)

# Define the expected ChromHMM states and allocate an in-memory results list.
state_ids <- 1:15
profile_results <- vector("list", nrow(manifest))

# Read each profile BED file and calculate base-pair coverage by state.
for (i in seq_len(nrow(manifest))) {
  sample_id <- manifest$sample[i]
  bed_file <- file.path(segment_dir, paste0(sample_id, "_15_segments.bed.gz"))

  if (!file.exists(bed_file)) {
    stop("Missing BED file: ", bed_file)
  }
  seg <- fread(bed_file, header = FALSE, col.names = c("chrom", "start", "end", "state"))

  if (nrow(seg) == 0L || anyNA(seg) || any(seg$start < 0) || any(seg$end <= seg$start)) {
    stop("Invalid or empty BED file: ", bed_file)
  }

  # Validate sequence names against mm10.txt without filtering any sequences.
  unknown_contigs <- setdiff(unique(seg$chrom), reference_contigs)

  if (length(unknown_contigs) > 0L) {
    stop("BED contains contigs absent from mm10.txt in ", bed_file, ": ", paste(unknown_contigs, collapse = ", "))
  }

  # Convert state labels E1-E15 to integer state IDs.
  seg[, state_id := as.integer(sub("^E", "", as.character(state)))]

  if (anyNA(seg$state_id) || any(!seg$state_id %in% state_ids)) {
    stop("Expected states E1-E15 in BED file: ", bed_file)
  }

  # Check that intervals tile each represented sequence continuously from 0.
  setorder(seg, chrom, start, end)
  seg[, previous_end := shift(end), by = chrom]

  if (any(seg[!is.na(previous_end), start != previous_end])) {
    stop("BED contains a gap or overlap: ", bed_file)
  }

  first_segments <- seg[, .(first_start = start[1L]), by = chrom]

  if (any(first_segments$first_start != 0)) {
    stop("A represented contig does not start at coordinate 0: ", bed_file)
  }

  # Sum BED interval widths by state; no chromosome or interval masking is applied.
  state_bp <- seg[
    ,
    .(bp = sum(as.numeric(end) - as.numeric(start))),
    by = state_id
  ]

  # Add rows for states with no intervals in a profile.
  full_states <- merge(
    data.table(state_id = state_ids),
    state_bp,
    by = "state_id",
    all.x = TRUE
  )

  full_states[is.na(bp), bp := 0]

  # Calculate the total number of bases represented in this profile's BED file.
  segment_genome_bp <- sum(full_states$bp)

  if (segment_genome_bp > reference_genome_bp) {
    stop("Segmented BED length exceeds mm10 reference length: ", bed_file)
  }

  # Calculate each state's percentage of the full mm10 reference genome.
  full_states[, sample := sample_id]
  full_states[, tissue := manifest$tissue[i]]
  full_states[, celltype := manifest$celltype[i]]
  full_states[, age := manifest$age[i]]
  full_states[, composite_rep := manifest$composite_rep[i]]
  full_states[, segment_genome_bp := segment_genome_bp]
  full_states[, reference_genome_bp := reference_genome_bp]
  full_states[, genome_percent := 100 * bp / reference_genome_bp]

  profile_results[[i]] <- full_states

  if (i %% 20L == 0L || i == nrow(manifest)) {
    message("Processed ", i, "/", nrow(manifest), " profiles")
  }
}

# Combine profiles in memory; do not write an intermediate coverage table.
coverage <- rbindlist(profile_results, use.names = TRUE)
if (any(coverage$bp < 0) || any(coverage$genome_percent < 0)) {
  stop("Calculated state coverage contains negative values")
}

coverage[, state_label := factor(
  paste0("E", state_id),
  levels = paste0("E", 15:1)
)]
coverage$state_label <- factor(coverage$state_label, levels = rev(names(state_label)))
saveRDS(coverage, file.path(model_dir, "emissions.coverage.rds"))

# Plot the genome-coverage distribution for each ChromHMM state.
p11 <- ggplot(
  coverage,
  aes(x = genome_percent, y = state_label)
) +
  geom_boxplot(
    width = 0.62,
    fill = "#D9D9D9",
    colour = "#333333",
    linewidth = 0.38,
    # outlier.shape = 21,
    # outlier.size = 0.8,
    # outlier.stroke = 0.25,
    # outlier.fill = "#555555",
    outlier.shape = NA
  ) +
  scale_x_continuous(
    # Pseudo-log scaling preserves zero values while showing a wide data range.
    trans = scales::transform_pseudo_log(sigma = 0.001),
    breaks = c(0, 0.001, 0.01, 0.1, 1, 10, 50),
    labels = c("0", "0.001", "0.01", "0.1", "1", "10", "50"),
    expand = expansion(mult = c(0.025, 0.04))
  ) +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = "Genome coverage (%)",
    y = "ChromHMM state",
    title = "Genome coverage by ChromHMM state",
    subtitle = paste0(
      "Each box summarizes ",
      uniqueN(coverage$sample),
      " profiles; denominator = mm10 reference genome (",
      format(reference_genome_bp, big.mark = ","),
      " bp)"
    )
  ) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(
      size = 7.5,
      colour = "#555555",
      margin = margin(b = 7)
    ),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7.5, colour = "#222222"),
    axis.text.y = element_text(margin = margin(r = 4)),
    axis.line = element_line(linewidth = 0.35, colour = "#222222"),
    axis.ticks = element_line(linewidth = 0.3, colour = "#222222"),
    axis.ticks.length = grid::unit(1.5, "mm"),
    plot.margin = margin(7, 8, 5, 5, unit = "mm")
  )

# Plot total covered base pairs for each state across profiles.
p22 <- ggplot(
  coverage,
  aes(x = bp, y = state_label)
) +
  geom_boxplot(
    width = 0.62,
    fill = "#D9D9D9",
    colour = "#333333",
    linewidth = 0.38,
    # outlier.size = 0.8,
    # outlier.stroke = 0.25,
    # outlier.fill = "#555555",
    outlier.shape = NA
  ) +
  scale_x_continuous(
    # Pseudo-log scaling preserves zero values and accommodates a wide range.
    trans = scales::transform_pseudo_log(sigma = 1000),
    breaks = c(1e5, 1e6, 1e7, 1e8, 1e9),
    labels = c("100 kb", "1 Mb", "10 Mb", "100 Mb", "1 Gb"),
    expand = expansion(mult = c(0.025, 0.04))
  ) +
  coord_cartesian(
  xlim = c(1e5, max(coverage$bp, na.rm = TRUE)),
  expand = FALSE
) +
  scale_y_discrete(drop = FALSE) +
  labs(
    x = "Coverage length (bp)",
    y = "ChromHMM state",
    title = "Coverage length by ChromHMM state",
    subtitle = paste0(
      "Each box summarizes ",
      uniqueN(coverage$sample),
      " profiles; values are summed BED interval lengths"
    )
  ) +
  theme_classic(base_size = 8, base_family = "sans") +
  theme(
    plot.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(
      size = 7.5,
      colour = "#555555",
      margin = margin(b = 7)
    ),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7.5, colour = "#222222"),
    axis.text.y = element_text(margin = margin(r = 4)),
    axis.line = element_line(linewidth = 0.35, colour = "#222222"),
    axis.ticks = element_line(linewidth = 0.3, colour = "#222222"),
    axis.ticks.length = grid::unit(1.5, "mm"),
    plot.margin = margin(7, 8, 5, 5, unit = "mm")
  )

pdf(file.path(model_dir, "emissions.pdf"), width = 6, height = 3.5)
print(p1)
print(p2)
print(p2_scale)
print(ggarrange(p11, p22, ncol = 4, nrow = 1, align = "hv"))
dev.off()

#### ---- state of genomic annotation ---- ####
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPseeker)
library(data.table)
workPath <- "~/longzhilin/Projects/AgingMouse/PairedTag/4.celltypeLevel/ChromHMM/State15/MergeAllTissue"
outDir <- file.path(workPath, "AnalysisResults")
dir.create(outDir, recursive = T, showWarnings = F)
setwd(outDir)
segments_1k_files <- list.files(file.path(workPath, "learnModel/split_1k"), "segments_1k.bed", full.names = T, recursive = T)
info <- data.frame(file = segments_1k_files)
info$sample <- gsub("_15_segments_1k.bed", "", basename(info$file))
info$tissue <- gsub(":.*", "", info$sample)
info$age <- gsub(".*_Old_.*", "Old", info$sample)
info$age <- gsub(".*_Young_.*", "Young", info$age)
info$cellType <- gsub("_Old_.*|_Young_.*", "", info$sample)

segments_1k_anno <- lapply(segments_1k_files, function(segments_1k_file){
  segments_1k <- fread(segments_1k_file, header = F, col.names = c("chrom", "start", "end", "state"))
  segments_1k.gr <- makeGRangesFromDataFrame(segments_1k, keep.extra.columns = T, starts.in.df.are.0based = TRUE) # 1-based start coordinate
  segments_1k_anno <- annotatePeak(segments_1k.gr, TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, tssRegion = c(-3000, 3000), verbose = FALSE)
  segments_1k_anno <- as.data.frame(segments_1k_anno)
  segments_1k_anno$annotationType <- gsub(" \\(.*\\)", "", segments_1k_anno$annotation)
  idx <- which(segments_1k_anno$annotationType == "Promoter")
  segments_1k_anno$annotationType[idx] <- segments_1k_anno$annotation[idx]
  segments_1k_anno$sample <- gsub("_15_segments_1k.bed", "", basename(segments_1k_file))
  return(segments_1k_anno)
})
segments_1k_anno <- do.call("rbind", segments_1k_anno)
segments_1k_anno <- as.data.table(segments_1k_anno)
saveRDS(segments_1k_anno, "segments_1k_anno.rds")
state_annotation_ratio <- segments_1k_anno[
  ,
  .(n_regions = .N),
  by = .(sample, state, annotationType)
][
  ,
  `:=`(
    proportion = n_regions / sum(n_regions),
    percent = 100 * n_regions / sum(n_regions)
  ),
  by = .(sample, state)
][
  order(sample, state, annotationType)
]
saveRDS(state_annotation_ratio, "state_annotation_ratio.rds")

#### plot the proportion of each state in different genomic annotation
library(data.table)
library(ggplot2)
plot_data <- as.data.table(copy(state_annotation_ratio))
# Extract age from sample names such as BAT:ASPC_Old_rep_1.
plot_data[
  ,
  age := fifelse(
    grepl("_Young_", sample),
    "Young",
    fifelse(grepl("_Old_", sample), "Old", NA_character_)
  )
]

if (anyNA(plot_data$age)) {
  stop("Could not determine age for one or more samples.")
}

annotation_levels <- c(
  "Promoter (<=1kb)",
  "Promoter (1-2kb)",
  "Promoter (2-3kb)",
  "Exon",
  "Intron",
  "Downstream",
  "Distal Intergenic",
  "3' UTR",
  "5' UTR"
)

unexpected_annotations <- setdiff(
  unique(plot_data$annotationType),
  annotation_levels
)

if (length(unexpected_annotations) > 0L) {
  stop(
    "Unexpected annotation categories: ",
    paste(unexpected_annotations, collapse = ", ")
  )
}

# Sum region counts across samples within each age, state, and annotation.
age_state_counts <- plot_data[
  ,
  .(n_regions = sum(n_regions)),
  by = .(age, state, annotationType)
]

# Complete absent annotation categories with zero counts.
state_grid <- unique(age_state_counts[, .(age, state)])
plot_grid <- state_grid[rep(seq_len(nrow(state_grid)), each = length(annotation_levels))]

plot_grid[
  ,
  annotationType := rep(annotation_levels, times = nrow(state_grid))
]

age_state_counts <- merge(
  plot_grid,
  age_state_counts,
  by = c("age", "state", "annotationType"),
  all.x = TRUE,
  sort = FALSE
)
age_state_counts[is.na(n_regions), n_regions := 0L]

# Calculate annotation percentages within each age and state.
age_state_counts[
  ,
  total_regions := sum(n_regions),
  by = .(age, state)
]

age_state_counts[
  ,
  percent := fifelse(total_regions > 0, 100 * n_regions / total_regions, 0)
]

# Order states numerically.
state_levels <- names(state_label)
state_numbers <- as.integer(gsub("[^0-9]", "", state_levels))
if (anyNA(state_numbers)) {
  stop("State labels must contain numeric IDs, for example E1, E2, ..., E15.")
}
age_state_counts[, age := factor(age, levels = c("Young", "Old"))]
age_state_counts[, state := factor(state, levels = rev(state_levels))]
age_state_counts[, annotationType := factor(annotationType, levels = rev(annotation_levels))]
age_state_counts$state <- factor(age_state_counts$state, levels = rev(state_levels))

annotation_colors <- c(
  "Promoter (<=1kb)" = "#B6385B",
  "Promoter (1-2kb)" = "#F28A5B",
  "Promoter (2-3kb)" = "#F5E75D",
  "Exon" = "#159A9B",
  "Intron" = "#2ED5BF",
  "Downstream" = "#183B61",
  "Distal Intergenic" = "#477DB0",
  "3' UTR" = "#C69BE0",
  "5' UTR" = "#B7A5D8"
)

legend_breaks <- rev(annotation_levels)
legend_labels <- rev(annotation_levels)
p <- ggplot(
  age_state_counts,
  aes(x = percent, y = state, fill = annotationType)
) +
  geom_col(
    width = 0.82,
    color = "black",
    linewidth = 0.25
  ) +
  facet_grid(. ~ age) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 25),
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = annotation_colors,
    breaks = legend_breaks,
    labels = legend_labels,
    drop = FALSE
  ) +
  labs(
    x = "% of regions",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.y = element_text(face = "bold", color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(face = "bold"),
    legend.position = "right",
    panel.spacing.x = grid::unit(0.8, "lines")
  )

pdf("chromHMM_state_genomic_annotations.pdf", width = 5, height = 5)
print(p)
dev.off()
