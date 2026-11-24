#!/bin/bash
#' @description ChromHMM state stranfer 
#' @author Zhilin Long
#' @date 2026-3-12
# conda activate PairedTagAnalysis

library(dplyr)
library(tidyr)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
set.seed(101)
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
state_levels <- names(state_label)
bam_information <- readRDS("~/longzhilin/Projects/AgingMouse/PairedTag/4.celltypeLevel/ChromHMM/DNA.depths.rds")
bam_information$cellType <- paste0(bam_information$tissue, ":", bam_information$cellType)
bam_information <- bam_information[grep("_rep_1|_rep_2", bam_information$sample),]
workPath <- "~/longzhilin/Projects/AgingMouse/PairedTag/4.celltypeLevel/ChromHMM/State15/MergeAllTissue"
outDir <- file.path(workPath, "AnalysisResults")
dir.create(outDir, recursive = T, showWarnings = F)
setwd(outDir)

#### ---- calculate_enrichment_score_cross_comparison ---- ####
segments_1k_files <- list.files(file.path(workPath, "learnModel/split_1k"), "segments_1k.bed", full.names = T, recursive = T)
info <- data.frame(file = segments_1k_files)
info$sample <- gsub("_15_segments_1k.bed", "", basename(info$file))
info$tissue <- gsub(":.*", "", info$sample)
info$age <- gsub(".*_Old_.*", "Old", info$sample)
info$age <- gsub(".*_Young_.*", "Young", info$age)
info$cellType <- gsub("_Old_.*|_Young_.*", "", info$sample)

#### prepare the background transfer matrix
cellType.state.background <- lapply(unique(info$cellType), function(x){
  info.Young <- info[info$cellType == x & info$age == "Young",]
  info.Old <- info[info$cellType == x & info$age == "Old",]
  cat(x, ": Young(", nrow(info.Young), "), Old(", nrow(info.Old), "), only select top 2 replicates\n")
  if(nrow(info.Young) < 1 | nrow(info.Old) < 1){
    cat(x, " does have less than 2 replication in Young or Old group!\n")
    return(NULL)
  }else{
    info.Young <- info.Young[grep("rep_1|rep_2", info.Young$sample),]
    info.Old <- info.Old[grep("rep_1|rep_2", info.Old$sample),]
    # Young
    Young_1 <- read.delim(info.Young$file[1], header = F, stringsAsFactors = F)
    Young_1$label <- paste0(Young_1$V1, "-", Young_1$V2 ,"-", Young_1$V3)
    states <- sort(unique(Young_1$V4))
    Young_background_transfer_matrix <- lapply(states, function(state){
      state_region <- Young_1[which(Young_1$V4 == state),,drop=F]
      Young_state_list <- lapply(2:nrow(info.Young), function(i){
        Young_others <- read.delim(info.Young$file[i], header = F, stringsAsFactors = F)
        Young_others$label <- paste0(Young_others$V1, "-", Young_others$V2, "-", Young_others$V3)
        transfer_state <- Young_others[which(Young_others$label %in% state_region$label),,drop=F]
        t_transfer_matrix <- as.data.frame(table(transfer_state$V4))
        colnames(t_transfer_matrix) <- paste0("Young_", i, c("_state", "_Freq"))
        t_transfer_matrix <- data.frame(Young_1_state = state, t_transfer_matrix)
        t_transfer_matrix$transfer_state <- paste0(t_transfer_matrix[,1], "-", t_transfer_matrix[,2])
        return(t_transfer_matrix)
      })
      Young_background_list <- Reduce(function(x,y) merge(x, y, all = T, by = c("transfer_state", "Young_1_state")), Young_state_list)
      return(Young_background_list)
    })
    Young_background_transfer_matrix <- Reduce(function(x,y) rbind(x, y), Young_background_transfer_matrix)

    # Old
    Old_1 <- read.delim(info.Old$file[1], header = F, stringsAsFactors = F)
    Old_1$label <- paste0(Old_1$V1, "-", Old_1$V2 ,"-", Old_1$V3)
    states <- sort(unique(Old_1$V4))
    Old_background_transfer_matrix <- lapply(states, function(state){
      state_region <- Old_1[which(Old_1$V4 == state),,drop=F]
      Old_state_list <- lapply(2:nrow(info.Old), function(i){
        Old_others <- read.delim(info.Old$file[i], header = F, stringsAsFactors = F)
        Old_others$label <- paste0(Old_others$V1, "-", Old_others$V2, "-", Old_others$V3)
        transfer_state <- Old_others[which(Old_others$label %in% state_region$label),,drop=F]
        t_transfer_matrix <- as.data.frame(table(transfer_state$V4))
        colnames(t_transfer_matrix) <- paste0("Old_", i, c("_state", "_Freq"))
        t_transfer_matrix <- data.frame(Old_1_state = state, t_transfer_matrix)
        t_transfer_matrix$transfer_state <- paste0(t_transfer_matrix[,1], "-", t_transfer_matrix[,2])
        return(t_transfer_matrix)
      })
      Old_background_list <- Reduce(function(x,y) merge(x, y, all = T, by = c("transfer_state", "Old_1_state")), Old_state_list)
      return(Old_background_list)
    })
    Old_background_transfer_matrix <- Reduce(function(x,y) rbind(x, y), Old_background_transfer_matrix)

    t_background_transfer_matrix <- merge(Young_background_transfer_matrix, Old_background_transfer_matrix, by = "transfer_state", all = T)
    # Complete the state-pair grid so an unobserved transition in an available comparison contributes zero.
    background_pair_grid <- expand.grid(rep1_state = state_levels, rep2_state = state_levels, stringsAsFactors = FALSE)
    background_pair_grid$transfer_state <- paste(background_pair_grid$rep1_state, background_pair_grid$rep2_state, sep = "-")
    t_background_transfer_matrix <- merge(background_pair_grid["transfer_state"], t_background_transfer_matrix, by = "transfer_state", all.x = TRUE, sort = FALSE)
    background_freq_columns <- grep("_Freq$", colnames(t_background_transfer_matrix), value = TRUE)
    t_background_transfer_matrix[background_freq_columns] <- lapply(t_background_transfer_matrix[background_freq_columns], function(x) {
      x[is.na(x)] <- 0
      x
    })
    t_background_transfer_matrix$Freq <- rowMeans(t_background_transfer_matrix[, background_freq_columns, drop = FALSE])
    t_background_transfer_matrix <- t_background_transfer_matrix %>% separate(col = transfer_state, into = c("rep1_state", "rep2_state"), sep = "-", remove = TRUE, convert = FALSE)
    t_background_transfer_matrix$tissue_cellType <- x
    return(t_background_transfer_matrix)
  }
})
names(cellType.state.background) <- unique(info$cellType)
cellType.state.background <- Reduce(function(x,y) rbind(x,y), cellType.state.background)
saveRDS(cellType.state.background, "cellType.state.background.rds")

#### state: young to old
cellType.state.aging <- lapply(unique(info$cellType), function(x){
  info.Young <- info[info$cellType == x & info$age == "Young",]
  info.Old <- info[info$cellType == x & info$age == "Old",]
  
  cat(x, ": Young(", nrow(info.Young), "), Old(", nrow(info.Old), "), only select top 2 replicates\n")
  if(nrow(info.Young) < 1 | nrow(info.Old) < 1){
    cat(x, " does have less than 2 replication in Young or Old group!\n")
    return(NULL)
  }else{
    info.Young <- info.Young[grep("rep_1|rep_2", info.Young$sample),]
    info.Old <- info.Old[grep("rep_1|rep_2", info.Old$sample),]

    transfer_matrix_list <- lapply(1:nrow(info.Young), function(i){
      Young <- read.delim(info.Young$file[i], header = F, stringsAsFactors = F)
      Young$label <- paste0(Young$V1, "-", Young$V2 ,"-", Young$V3)
      Old <- read.delim(info.Old$file[i], header = F, stringsAsFactors = F)
      Old$label <- paste0(Old$V1, "-", Old$V2 ,"-", Old$V3)
      Young <- Young[which(Young$label %in% Old$label),,drop=F]
      Old <-  Old[which(Old$label %in% Young$label),,drop=F]
      states <- sort(unique(Young$V4))

      state.list <- lapply(states, function(state){
        state_region <- Young[which(Young$V4 == state),,drop=F]
        transfer_state <- Old[which(Old$label %in% state_region$label),]
        transfer_matrix <- as.data.frame(table(transfer_state$V4))
        colnames(transfer_matrix) <- c("Old_state", "Freq")
        t_transfer_matrix <- data.frame(tissue_cellType = x, Young_state = state, transfer_matrix)
        t_transfer_matrix$label <- paste0(t_transfer_matrix$Young_state, "-", t_transfer_matrix$Old_state)
        return(t_transfer_matrix)
      })
      state.list <- Reduce(function(x,y) rbind(x, y), state.list)
      return(state.list)
    })
    extract_columns <- lapply(transfer_matrix_list, function(df) {select(df, label, Freq)})
    transfer_matrix_summary <- Reduce(function(x,y) merge(x, y, all = T, by = "label"), extract_columns)
    # Complete all 15 x 15 state pairs; absent pairs in an observed comparison are zero counts.
    aging_pair_grid <- expand.grid(Young_state = state_levels, Old_state = state_levels, stringsAsFactors = FALSE)
    aging_pair_grid$label <- paste(aging_pair_grid$Young_state, aging_pair_grid$Old_state, sep = "-")
    transfer_matrix_summary <- merge(aging_pair_grid["label"], transfer_matrix_summary, by = "label", all.x = TRUE, sort = FALSE)
    aging_freq_columns <- setdiff(colnames(transfer_matrix_summary), "label")
    transfer_matrix_summary[aging_freq_columns] <- lapply(transfer_matrix_summary[aging_freq_columns], function(x) {
      x[is.na(x)] <- 0
      x
    })

    transfer_matrix_summary$average_Freq <- rowMeans(transfer_matrix_summary[, aging_freq_columns, drop = FALSE])
    transfer_matrix_summary <- transfer_matrix_summary[,c("label","average_Freq")]

    t_aging_transfer_matrix <- transfer_matrix_summary %>% separate(col = label, into = c("Young_state", "Old_state"), sep = "-", remove = TRUE, convert = FALSE)
    t_aging_transfer_matrix$tissue_cellType <- x
    colnames(t_aging_transfer_matrix)[3] <- "Freq"
    t_aging_transfer_matrix <- t_aging_transfer_matrix[,c("tissue_cellType", "Young_state", "Old_state", "Freq")]
    return(t_aging_transfer_matrix)
  }
})
names(cellType.state.aging) <- unique(info$cellType)
cellType.state.aging <- Reduce(function(x,y) rbind(x,y), cellType.state.aging)
saveRDS(cellType.state.aging, "cellType.state.aging.rds")

#### genome coverage
genomeCoverage.list <- lapply(1:nrow(info), function(i){
  file <- info$file[i]
  segments_1k <- read.delim(file, header = F, stringsAsFactors = F)
  segments_1k$width <- as.numeric(segments_1k$V3) - as.numeric(segments_1k$V2)
  genomeCoverage <- segments_1k %>% group_by(V4) %>% summarise(total_width = sum(width), .groups = 'drop')
  genomeCoverage$tissue_cellType <- info$cellType[i]
  genomeCoverage$sample <- info$sample[i]

  file <- gsub("/split_1k", "", file)
  file <- gsub("_segments_1k.bed", "_overlap.txt", file)
  genomeOverlap <- read.delim(file, header = T, stringsAsFactors = F)
  genomeOverlap <- genomeOverlap[1:15,]
  genomeOverlap$V4 <- paste0("E", genomeOverlap[,1])
  idx <- match(genomeCoverage$V4, genomeOverlap$V4)
  genomeCoverage$genomePercentage <- as.numeric(genomeOverlap[idx,2])
  return(genomeCoverage)
})
genomeCoverage <- Reduce(function(x,y) rbind(x, y), genomeCoverage.list)
genomeCoverage <- as.data.frame(genomeCoverage)
saveRDS(genomeCoverage, "genomeCoverage.rds")

#### analysis different levels
bam_depths <- bam_information[which(bam_information$cellType %in% info$cellType),]
bam_depths <- tapply(bam_depths$reads, bam_depths$histone, sum)
bam_depths <- round(bam_depths / 1000000, 2)

#### plot heatmap for merge all celltypes
cellType.state.aging$Young_state <- factor(cellType.state.aging$Young_state, levels = paste0("E", 1:15))
cellType.state.aging$Old_state <- factor(cellType.state.aging$Old_state, levels = paste0("E", 1:15))
cellType.state.background$rep1_state <- factor(cellType.state.background$rep1_state, levels = paste0("E", 1:15))
cellType.state.background$rep2_state <- factor(cellType.state.background$rep2_state, levels = paste0("E", 1:15))

aging_transfer_matrix <- cellType.state.aging %>% group_by(Young_state, Old_state) %>% summarise(mean_freq = mean(Freq))
aging_transfer_matrix2 <- aging_transfer_matrix 
aging_transfer_matrix2$percent_freq <- aging_transfer_matrix2$mean_freq / sum(aging_transfer_matrix2$mean_freq) * 100
aging_transfer_matrix_to_plot <- reshape2::dcast(aging_transfer_matrix2, formula = Young_state~Old_state, value.var = "percent_freq")
rownames(aging_transfer_matrix_to_plot) <- aging_transfer_matrix_to_plot[,1]
aging_transfer_matrix_to_plot <- aging_transfer_matrix_to_plot[,-1]

background_transfer_matrix <- cellType.state.background %>% group_by(rep1_state, rep2_state) %>% summarise(mean_freq = mean(Freq))
background_transfer_matrix2 <- background_transfer_matrix 
background_transfer_matrix2$percent_freq <- background_transfer_matrix2$mean_freq / sum(background_transfer_matrix2$mean_freq) * 100
background_transfer_matrix_to_plot <- reshape2::dcast(background_transfer_matrix2, formula = rep1_state~rep2_state, value.var = "percent_freq")
rownames(background_transfer_matrix_to_plot) <- background_transfer_matrix_to_plot[,1]
background_transfer_matrix_to_plot <- background_transfer_matrix_to_plot[,-1]

aging <- as.matrix(aging_transfer_matrix_to_plot)
background <- as.matrix(background_transfer_matrix_to_plot)
stopifnot(identical(dimnames(aging), dimnames(background)))
result_to_plot <- matrix(NA_real_, nrow = nrow(aging),ncol = ncol(aging), dimnames = dimnames(aging))
valid <- is.finite(aging) & is.finite(background) & background > 0
result_to_plot[valid] <- aging[valid] / background[valid]

result_to_plot2 <- matrix(NA_real_, nrow = nrow(result_to_plot), ncol = ncol(result_to_plot), dimnames = dimnames(result_to_plot))
valid_log <- is.finite(result_to_plot) & result_to_plot > 0
result_to_plot2[valid_log] <- log2(result_to_plot[valid_log])
diag(result_to_plot2) <- NA

# just consider the off-diagonal elements for bubble plot
aging_transfer_matrix2 <- aging_transfer_matrix 
aging_transfer_matrix2$mean_freq[which(aging_transfer_matrix2$Young_state == aging_transfer_matrix2$Old_state)] <- NA
aging_transfer_matrix2$percent_freq <- aging_transfer_matrix2$mean_freq / sum(aging_transfer_matrix2$mean_freq, na.rm = T) * 100
bubble_aging_transfer_matrix_to_plot <- reshape2::dcast(aging_transfer_matrix2, formula = Young_state~Old_state, value.var = "percent_freq") 
rownames(bubble_aging_transfer_matrix_to_plot) <- bubble_aging_transfer_matrix_to_plot$Young_state
bubble_aging_transfer_matrix_to_plot <- bubble_aging_transfer_matrix_to_plot[,-1]

# change order
bubble_aging_transfer_matrix_to_plot <- bubble_aging_transfer_matrix_to_plot[state_levels, state_levels]
result_to_plot2 <- result_to_plot2[state_levels, state_levels]

p1 <- Heatmap(result_to_plot2, name = "log2(aging/background)", col = col_fun, cluster_rows = F, cluster_columns = F, height = unit(10, "cm"), width = unit(10, "cm"), cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.2f", result_to_plot2[i, j]), x, y, gp = gpar(fontsize = 7))}, row_names_side = "left", column_names_side = "top")

library(bubbleHeatmap)
min.value <- floor(min(na.omit(as.numeric(unlist(result_to_plot2)))))
max.value <- floor(max(na.omit(as.numeric(unlist(result_to_plot2)))))
max.value2 <- floor(max(na.omit(as.numeric(unlist(bubble_aging_transfer_matrix_to_plot)))))
col_fun <- colorRamp2(c(min.value, 0, max.value), c("blue", "white", "red"))
customColors <- col_fun(seq(min.value, max.value, length.out = 100))
tree <- bubbleHeatmap(as.matrix(result_to_plot2), sizeMat = as.matrix(bubble_aging_transfer_matrix_to_plot), colorLim = c(min.value, max.value), sizeLim = c(0, max.value2), legendTitles = c("Proportion", "log2(Enrichment score)"), colorSeq = customColors, xTitle = "Old state(Old)",yTitle = "Young state(Young)")
grDevices::cairo_pdf("MergeTissue.pdf")
print(p1)
grid.newpage()
grid.draw(tree)
dev.off()

#### ---- Enrichment of chromatin states relative to annotated genes ---- ####
# https://www.nature.com/articles/s41586-020-2093-3/figures/11
library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(EnrichedHeatmap)
workPath <- "~/longzhilin/Projects/AgingMouse/PairedTag/4.celltypeLevel/ChromHMM/State15/MergeAllTissue"
inputDir <- file.path(workPath, "learnModel")
outDir <- file.path(workPath, "AnalysisResults")

if (!dir.exists(inputDir)) {
  stop("ChromHMM input directory does not exist: ", inputDir)
}
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Define 20 bins: five upstream, 10 gene-body, and five downstream bins.
extendBp <- 5000L
windowWidth <- 1000L
targetRatio <- 0.5
nExtendedWindows <- 2 * extendBp / windowWidth
nTargetWindows <- nExtendedWindows * targetRatio / (1 - targetRatio)
nPositions <- nExtendedWindows + nTargetWindows

if (any(abs(c(nExtendedWindows, nTargetWindows, nPositions) -
           round(c(nExtendedWindows, nTargetWindows, nPositions))) > 1e-8) ||
    nPositions != 20L) {
  stop("Window settings must produce 20 integer bins.")
}
nExtendedWindows <- as.integer(nExtendedWindows)
nTargetWindows <- as.integer(nTargetWindows)
nPositions <- as.integer(nPositions)

# Represent each gene by its longest annotated transcript.
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
tx <- transcripts(txdb, columns = "gene_id")
txGeneIdList <- mcols(tx)$gene_id
txGeneIdCount <- S4Vectors::elementNROWS(txGeneIdList)
txGeneId <- rep(NA_character_, length(tx))
hasSingleGeneId <- txGeneIdCount == 1L
txGeneId[hasSingleGeneId] <- as.character(
  unlist(txGeneIdList[hasSingleGeneId], use.names = FALSE)
)

txTable <- data.table(
  seqnames = as.character(seqnames(tx)),
  start = start(tx),
  end = end(tx),
  strand = as.character(strand(tx)),
  gene_id = txGeneId
)
txTable <- txTable[!is.na(gene_id) & nzchar(gene_id)]
txTable[, transcript_width := end - start + 1L]
setorder(txTable, gene_id, -transcript_width, seqnames, start, end)
geneTable <- txTable[!duplicated(gene_id)]

geneGr <- makeGRangesFromDataFrame(
  geneTable,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand",
  keep.extra.columns = TRUE
)
if (length(geneGr) == 0L) {
  stop("No gene intervals were obtained from the mm10 TxDb.")
}
rm(tx, txGeneIdList, txGeneIdCount, txGeneId, txTable, geneTable)

segmentFiles <- list.files(
  inputDir,
  pattern = "_15_segments\\.bed(\\.gz)?$",
  full.names = TRUE,
  recursive = TRUE
)
if (length(segmentFiles) == 0L) {
  stop("No merged ChromHMM segment files were found under: ", inputDir)
}

sampleInfo <- data.table(file = segmentFiles)
sampleInfo[, sample := sub(
  "_15_segments\\.bed(\\.gz)?$",
  "",
  basename(file)
)]
sampleInfo[, age := fifelse(
  grepl("_Young_", sample), "Young",
  fifelse(grepl("_Old_", sample), "Old", NA_character_)
)]
if (anyNA(sampleInfo$age)) {
  stop("Could not determine age for: ",
       paste(sampleInfo[is.na(age), sample], collapse = ", "))
}
if (anyDuplicated(sampleInfo$sample)) {
  stop("Sample names are not unique after parsing segment filenames.")
}

# Include the window settings in the cache directory to prevent stale reuse.
cacheDir <- file.path(
  outDir,
  sprintf(
    "profile_cache_parallel_20bins_ext%d_w%d_ratio%s",
    extendBp,
    windowWidth,
    format(targetRatio, trim = TRUE)
  )
)
dir.create(cacheDir, recursive = TRUE, showWarnings = FALSE)

cacheNames <- gsub("[^[:alnum:]_.-]", "_", sampleInfo$sample)
if (anyDuplicated(cacheNames)) {
  stop("Sample names collide after cache filename sanitization.")
}
cacheFiles <- file.path(cacheDir, paste0(cacheNames, ".rds"))

profileOneSample <- function(file, sample, age) {
  segments <- fread(
    file,
    header = FALSE,
    select = 1:4,
    col.names = c("chrom", "start", "end", "state"),
    showProgress = FALSE
  )
  if (nrow(segments) == 0L) {
    stop("Segment file is empty: ", file)
  }

  segments[, state := as.character(state)]
  unknownStates <- setdiff(unique(segments$state), state_levels)
  if (length(unknownStates) > 0L) {
    stop(
      "Unexpected state labels in ", file, ": ",
      paste(unknownStates, collapse = ", ")
    )
  }

  # BED starts are zero-based and half-open; convert them to GRanges coordinates.
  segmentGr <- makeGRangesFromDataFrame(
    segments,
    seqnames.field = "chrom",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE,
    starts.in.df.are.0based = TRUE
  )
  segmentsByState <- split(
    segmentGr,
    factor(mcols(segmentGr)$state, levels = state_levels)
  )

  profileList <- lapply(state_levels, function(currentState) {
    stateGr <- segmentsByState[[currentState]]

    if (length(stateGr) == 0L) {
      profile <- numeric(nPositions)
    } else {
      # Merge adjacent intervals with the same state before calculating coverage.
      stateGr <- reduce(stateGr, ignore.strand = TRUE)
      stateMatrix <- normalizeToMatrix(
        signal = stateGr,
        target = geneGr,
        extend = extendBp,
        w = windowWidth,
        target_ratio = targetRatio,
        include_target = TRUE,
        mean_mode = "coverage",
        background = 0,
        smooth = FALSE,
        verbose = FALSE
      )

      if (ncol(stateMatrix) != nPositions) {
        stop(
          "Unexpected number of bins for state ", currentState,
          " in sample ", sample, ": ", ncol(stateMatrix)
        )
      }

      profile <- colMeans(stateMatrix, na.rm = TRUE)
      if (any(!is.finite(profile))) {
        stop("Non-finite profile values for state ", currentState,
             " in sample ", sample)
      }
    }

    data.table(
      age = age,
      sample = sample,
      state = currentState,
      position = seq_len(nPositions),
      coverage = profile
    )
  })

  rbindlist(profileList, use.names = TRUE)
}

cacheMetadata <- function(i) {
  inputInfo <- file.info(sampleInfo$file[[i]])
  list(
    cache_version = 1L,
    input_file = normalizePath(sampleInfo$file[[i]], mustWork = TRUE),
    input_size = unname(inputInfo$size),
    input_mtime = as.numeric(inputInfo$mtime),
    age = sampleInfo$age[[i]],
    sample = sampleInfo$sample[[i]],
    extend_bp = extendBp,
    window_width = windowWidth,
    target_ratio = targetRatio,
    state_levels = state_levels,
    gene_model = paste0(
      "TxDb.Mmusculus.UCSC.mm10.knownGene ",
      as.character(packageVersion("TxDb.Mmusculus.UCSC.mm10.knownGene"))
    )
  )
}

processOneSample <- function(i) {
  cacheFile <- cacheFiles[[i]]
  expectedMetadata <- cacheMetadata(i)

  if (file.exists(cacheFile)) {
    cached <- readRDS(cacheFile)
    if (is.list(cached) &&
        identical(cached$metadata, expectedMetadata) &&
        is.data.table(cached$profile)) {
      return(cached$profile)
    }
  }

  result <- profileOneSample(
    file = sampleInfo$file[[i]],
    sample = sampleInfo$sample[[i]],
    age = sampleInfo$age[[i]]
  )
  tempFile <- tempfile(
    pattern = paste0(".", cacheNames[[i]], "_"),
    tmpdir = cacheDir,
    fileext = ".tmp"
  )
  on.exit(unlink(tempFile), add = TRUE)
  saveRDS(
    list(metadata = expectedMetadata, profile = result),
    tempFile
  )
  if (!file.rename(tempFile, cacheFile)) {
    stop("Failed to write cache for sample: ", sampleInfo$sample[[i]])
  }
  result
}

workersRequested <- suppressWarnings(
  as.integer(Sys.getenv("PROFILE_WORKERS", unset = "4"))
)
if (is.na(workersRequested) || workersRequested < 1L) {
  stop("PROFILE_WORKERS must be a positive integer.")
}
if (.Platform$OS.type == "windows") {
  workers <- 1L
} else {
  workers <- min(workersRequested, nrow(sampleInfo))
}

message("Samples: ", nrow(sampleInfo), "; workers: ", workers)
sampleIndices <- seq_len(nrow(sampleInfo))
if (workers > 1L) {
  profileResults <- parallel::mclapply(
    sampleIndices,
    processOneSample,
    mc.cores = workers,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
} else {
  profileResults <- lapply(sampleIndices, processOneSample)
}

failed <- vapply(profileResults, inherits, logical(1), "try-error")
if (any(failed)) {
  stop(
    "Profile calculation failed for sample(s): ",
    paste(sampleInfo$sample[failed], collapse = ", ")
  )
}

sampleProfiles <- rbindlist(profileResults, use.names = TRUE)
if (anyDuplicated(sampleProfiles[, .(sample, state, position)])) {
  stop("Duplicate sample-state-position rows were generated.")
}
if (any(sampleProfiles$position < 1L | sampleProfiles$position > 20L)) {
  stop("Profile positions fall outside the expected 1-20 range.")
}

saveRDS(
  sampleProfiles,
  file.path(outDir, "sample_profiles_20bins_parallel.rds")
)
fwrite(
  sampleProfiles,
  file.path(outDir, "sample_profiles_20bins_parallel.tsv"),
  sep = "\t"
)

#### plot heatmap
library(data.table)
library(ggplot2)
library(scales)
analysisDir <- file.path(
  "~/longzhilin/Projects/AgingMouse/PairedTag",
  "4.celltypeLevel/ChromHMM/State15/MergeAllTissue/AnalysisResults"
)
cacheDir <- file.path(
  analysisDir,
  "profile_cache_parallel_20bins_ext5000_w1000_ratio0.5"
)

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
stateOrder <- names(state_label)
expectedPositions <- 20L

cacheFiles <- list.files(cacheDir, pattern = "\\.rds$", full.names = TRUE)
if (length(cacheFiles) == 0L) {
  stop("No per-sample cache files found in: ", cacheDir)
}

# The parallel profile script stores metadata and the sample profile together.
cacheObjects <- lapply(cacheFiles, readRDS)
hasProfile <- vapply(
  cacheObjects,
  function(x) is.list(x) && !is.null(x$profile),
  logical(1)
)
if (!all(hasProfile)) {
  stop("One or more cache files do not contain a profile entry.")
}
profileData <- rbindlist(
  lapply(cacheObjects, function(x) as.data.table(x$profile)),
  use.names = TRUE
)
rm(cacheObjects)

requiredColumns <- c("age", "sample", "state", "position", "coverage")
if (!all(requiredColumns %in% names(profileData))) {
  stop("Profile cache is missing required columns.")
}
if (anyNA(profileData[, ..requiredColumns])) {
  stop("Profile cache contains missing values in required fields.")
}
if (any(!is.finite(profileData$coverage)) || any(profileData$coverage < 0)) {
  stop("Coverage must be finite and non-negative.")
}
if (!all(profileData$age %in% c("Young", "Old"))) {
  stop("Unexpected age labels are present in the profile cache.")
}
if (!setequal(unique(as.character(profileData$state)), stateOrder)) {
  stop("Observed state labels do not match the expected 15-state order.")
}
if (!setequal(unique(as.integer(profileData$position)), seq_len(expectedPositions))) {
  stop("Observed positions do not match the expected 1-20 bins.")
}

# Retain cell types represented at both ages for a paired age comparison.
profileData[, tissue := sub(":.*$", "", sample)]
profileData[, cellType := sub("^[^:]+:", "", sample)]
profileData[, cellType := sub("_(Young|Old)_rep_[0-9]+$", "", cellType)]
sampleMetadata <- unique(profileData[, .(sample, age, tissue, cellType)])
if (anyDuplicated(sampleMetadata$sample)) {
  stop("A sample is assigned to more than one age or cell type.")
}
pairedCellTypes <- sampleMetadata[
  , .(n_ages = uniqueN(age)),
  by = .(tissue, cellType)
][n_ages == 2L, .(tissue, cellType)]
if (nrow(pairedCellTypes) == 0L) {
  stop("No cell types have profiles at both ages.")
}
profileData <- profileData[
  pairedCellTypes,
  on = .(tissue, cellType),
  nomatch = 0L
]

# Confirm that every retained sample contributes all 15 states and 20 bins.
sampleDimensions <- profileData[
  , .(n_rows = .N, n_states = uniqueN(state), n_positions = uniqueN(position)),
  by = .(sample, age, tissue, cellType)
]
if (any(sampleDimensions$n_rows != length(stateOrder) * expectedPositions) ||
    any(sampleDimensions$n_states != length(stateOrder)) ||
    any(sampleDimensions$n_positions != expectedPositions)) {
  stop("At least one sample has incomplete state-by-position profiles.")
}

# Average replicates within each cell type, then weight cell types equally.
cellTypeProfiles <- profileData[
  , .(celltype_coverage = mean(coverage)),
  by = .(tissue, cellType, age, state, position)
]
ageProfiles <- cellTypeProfiles[
  , .(
    mean_coverage = mean(celltype_coverage),
    sd_celltype = sd(celltype_coverage),
    n_celltypes = .N
  ),
  by = .(age, state, position)
]

## Apply row-wise min-max normalization to each age-state profile across 20 bins.
colorUpperLimit <- 1
outDir <- file.path(
  analysisDir,
  "profile_heatmap_tss_tes_20bins_minmax"
)
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
ageProfiles[
  , minmax_normalized_overlap := {
    rangeValue <- range(mean_coverage)
    if (any(!is.finite(rangeValue)) || diff(rangeValue) == 0) {
      rep(NA_real_, .N)
    } else {
      (mean_coverage - rangeValue[[1L]]) / diff(rangeValue)
    }
  },
  by = .(age, state)
]
if (anyNA(ageProfiles$minmax_normalized_overlap)) {
  stop("At least one age-state profile has a zero or invalid range.")
}

minmaxCheck <- ageProfiles[
  , .(row_min = min(minmax_normalized_overlap), row_max = max(minmax_normalized_overlap)),
  by = .(age, state)
]
if (nrow(minmaxCheck) != 30L || any(abs(minmaxCheck$row_min) > 1e-10) || any(abs(minmaxCheck$row_max - 1) > 1e-10)) {
  stop("Each non-constant age-state profile must span 0 to 1 after min-max normalization.")
}

ageProfiles[, state_label := unname(state_label[state])]
ageProfiles[, age := factor(age, levels = c("Young", "Old"))]
plotData <- copy(ageProfiles)
plotData[, state := factor(state, levels = rev(stateOrder))]

observedNormalizedMax <- max(plotData$minmax_normalized_overlap)
fractionAboveLimit <- mean(plotData$minmax_normalized_overlap > colorUpperLimit)
overlapColors <- c(
  "#f5f3db",
  "#F6E99A",
  "#D4E89A",
  "#83CFA7",
  "#36AEB5",
  "#2D79B8",
  "#1D2B7A"
)

# Positions 1-5 are upstream, 6-15 are the gene body, and 16-20 are downstream.
p <- ggplot(
  plotData,
  aes(x = position, y = state, fill = minmax_normalized_overlap)
) +
  geom_tile(width = 1, height = 0.98) +
  geom_vline(
    xintercept = c(5.5, 15.5),
    color = "black",
    linewidth = 0.55
  ) +
  geom_hline(
    yintercept = c(1.5, 3.5, 5.5, 9.5, 12.5),
    color = "black",
    linewidth = 0.55
  ) +
  facet_grid(. ~ age, drop = FALSE) +
  scale_x_continuous(
    position = "top",
    breaks = c(1, 5.5, 15.5, 20),
    labels = c("-5 kb", "TSS", "TTS", "+5 kb"),
    limits = c(0.5, 20.5),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    labels = function(x) unname(state_label[x]),
    expand = c(0, 0),
    drop = FALSE
  ) +
  scale_fill_gradientn(
    colours = overlapColors,
    limits = c(0, colorUpperLimit),
    breaks = seq(0, colorUpperLimit, by = 0.2),
    labels = label_number(accuracy = 0.01),
    oob = squish,
    name = "Min-max normalized\nOverlap (0-1)",
    na.value = "grey80"
  ) +
  labs(
    title = "ChromHMM state distribution around gene bodies",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(color = "black", face = "bold", size = 9),
    axis.text.y = element_text(color = "black", size = 8.5),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    panel.spacing.x = grid::unit(2.2, "lines")
  )

orderedProfiles <- ageProfiles[
  order(age, match(state, stateOrder), position)
]
saveRDS(
  orderedProfiles,
  file.path(outDir, "age_mean_profile_minmax_normalized.rds")
)
ggsave(
  filename = file.path(outDir, "chromHMM_state_TSS_TTS_heatmap_minmax.pdf"),
  plot = p,
  width = 10,
  height = 4,
  units = "in"
)

## L2-normalize each age-state profile across the 20 displayed bins
colorUpperLimit <- 0.35
outDir <- file.path(
  analysisDir,
  "profile_heatmap_tss_tes_20bins_l2"
)
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

ageProfiles[
  , l2_normalized_overlap := {
    normValue <- sqrt(sum(mean_coverage^2))
    if (!is.finite(normValue) || normValue == 0) {
      rep(NA_real_, .N)
    } else {
      mean_coverage / normValue
    }
  },
  by = .(age, state)
]
if (anyNA(ageProfiles$l2_normalized_overlap)) {
  stop("At least one age-state profile has zero or invalid L2 norm.")
}

l2NormCheck <- ageProfiles[
  , .(l2_norm = sqrt(sum(l2_normalized_overlap^2))),
  by = .(age, state)
]
if (nrow(l2NormCheck) != 30L || any(abs(l2NormCheck$l2_norm - 1) > 1e-10)) {
  stop("Each age-state profile must have unit L2 norm across 20 bins.")
}

ageProfiles[, state_label := unname(state_label[state])]
ageProfiles[, age := factor(age, levels = c("Young", "Old"))]
plotData <- copy(ageProfiles)
plotData[, state := factor(state, levels = rev(stateOrder))]

observedL2Max <- max(plotData$l2_normalized_overlap)
fractionAboveLimit <- mean(plotData$l2_normalized_overlap > colorUpperLimit)
overlapColors <- c(
  "#FFFDF2",
  "#F6E99A",
  "#D4E89A",
  "#83CFA7",
  "#36AEB5",
  "#2D79B8",
  "#1D2B7A"
)

# Positions 1-5 are upstream, 6-15 are the gene body, and 16-20 are downstream.
p <- ggplot(
  plotData,
  aes(x = position, y = state, fill = l2_normalized_overlap)
) +
  geom_tile(width = 1, height = 0.98) +
  geom_vline(
    xintercept = c(5.5, 15.5),
    color = "black",
    linewidth = 0.55
  ) +
  geom_hline(
    yintercept = c(1.5, 3.5, 5.5, 9.5, 12.5),
    color = "black",
    linewidth = 0.55
  ) +
  facet_grid(. ~ age, drop = FALSE) +
  scale_x_continuous(
    position = "top",
    breaks = c(1, 5.5, 15.5, 20),
    labels = c("-5 kb", "TSS", "TTS", "+5 kb"),
    limits = c(0.5, 20.5),
    expand = c(0, 0)
  ) +
  scale_y_discrete(
    labels = function(x) unname(state_label[x]),
    expand = c(0, 0),
    drop = FALSE
  ) +
  scale_fill_gradientn(
    colours = overlapColors,
    limits = c(0, colorUpperLimit),
    breaks = seq(0, colorUpperLimit, by = 0.1),
    labels = label_number(accuracy = 0.01),
    oob = squish,
    name = "L2-normalized\nOverlap",
    na.value = "grey80"
  ) +
  labs(
    title = "ChromHMM state distribution around gene bodies",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(color = "black", face = "bold", size = 9),
    axis.text.y = element_text(color = "black", size = 8.5),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    panel.spacing.x = grid::unit(2.2, "lines")
  )

orderedProfiles <- ageProfiles[
  order(age, match(state, stateOrder), position)
]
saveRDS(
  orderedProfiles,
  file.path(outDir, "age_mean_profile_l2_normalized.rds")
)
ggsave(
  filename = file.path(outDir, "chromHMM_state_TSS_TTS_heatmap_l2.pdf"),
  plot = p,
  width = 10,
  height = 4,
  units = "in"
)