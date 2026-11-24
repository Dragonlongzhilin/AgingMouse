#!/bin/bash
#' @description running ChromHMM method for each cell type in Paired-Tag data
#' @author Zhilin Long
#' @date 2026-3-6
# conda activate PairedTagAnalysis

ChromHMM_path=~/longzhilin/Softwares/Epigenomics/ChromHMM
workPath=~/longzhilin/Projects/AgingMouse/PairedTag/4.celltypeLevel/ChromHMM
cd $workPath

#### create the input file
# R code
source("~/Functions/ChromHMM_input.R")
tissues <- c("brainCB", "brainFC", "brainHip", "heart", "liver", "lung", "kidney", "BAT", "muscle", "spleen", "colon", "skin", "testis")
bam_read_counts <- read.table("~/longzhilin/public_html/AgingMouse/PairedTag_new/bam_read_counts.txt", header = T, sep = "\t", stringsAsFactors = F)
for(i in tissues){
  cat("Analysis ", i, "\n")
  res <- ChromHMM_input(bamPath = paste0("~/longzhilin/public_html/AgingMouse/PairedTag_new/", i), outPath = "~/longzhilin/Projects/AgingMouse/PairedTag/4.celltypeLevel/ChromHMM/State15", pattern = "cellType_age_replication", suffix = ".bam$", bam_read_counts = bam_read_counts)
}


#### ---- Merge all tissues ---- ####
#### create the input file
# R code
library(tidyverse)
input.list <- list.files("~/longzhilin/Projects/AgingMouse/PairedTag/4.celltypeLevel/ChromHMM/State15", "input.rds", full.names = T, recursive = T)
merge.cellmarkfile <- lapply(input.list, function(x){
  input <- readRDS(x)
  tissue <- basename(dirname(x))
  input$file <- gsub("\\./", paste0("./", tissue, "/"), input$file)
  input_filtered <- input %>% group_by(cellType) %>% filter(sum(reads > 1000000) / n() >= 0.5) %>% ungroup() # at least 50% of replicates have more than 1M reads
  return(input_filtered)
})
merge.cellmarkfile <- do.call(rbind, merge.cellmarkfile)
write.table(merge.cellmarkfile[,c("sample", "histone", "file")], file.path("~/longzhilin/Projects/AgingMouse/PairedTag/4.celltypeLevel/ChromHMM/State15/MergeAllTissue", "cellmarkfile.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

#### run ChromHMM
bamPath="~/longzhilin/public_html/AgingMouse/PairedTag_new"
outDir="~/longzhilin/Projects/AgingMouse/PairedTag/4.celltypeLevel/ChromHMM/State15/MergeAllTissue"
# BinarizeBam
mkdir -p $outDir/binarizeBam;
java -mx8000M -jar $ChromHMM_path/ChromHMM.jar BinarizeBam -gzip -b 1000 $ChromHMM_path/CHROMSIZES/mm10.txt $bamPath $outDir/cellmarkfile.txt $outDir/binarizeBam;
# LearnModel
mkdir -p $outDir/learnModel;
java -mx8000M -jar $ChromHMM_path/ChromHMM.jar LearnModel -b 1000 -gzip -p 24 -holdcolumnorder $outDir/binarizeBam $outDir/learnModel 15 mm10;
