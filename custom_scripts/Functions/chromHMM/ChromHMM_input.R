#' @description construct the input file for ChromHMM
#' @author Zhilin Long
#' @date 2026-3-6

ChromHMM_input <- function(bamPath, outPath, pattern = "cellType_age_replication", suffix = ".bam$", bam_read_counts = NULL){
  require(dplyr)
  bamFiles <- list.files(bamPath, pattern = suffix, full.names = T, recursive = T)
  bamFiles <- bamFiles[grep(pattern, bamFiles)]

  # count reads
  if(is.null(bam_read_counts)){
    reads <- lapply(bamFiles, function(i){
      cmd <- paste0("samtools view -c -@ 8 ", i)
      read <- system(cmd, intern = T)
      return(read)
    })
    reads <- as.numeric(unlist(reads))
  }else{
    reads <- bam_read_counts$read_count[match(bamFiles, bam_read_counts$bam_file)]
  }

  bamFiles <- gsub(bamPath, ".", bamFiles)
  names(reads) <- bamFiles
  fileName <- gsub(suffix, "", basename(bamFiles))
  age <- gsub(".*_Old_.*", "Old", fileName)
  age <- gsub(".*_Young_.*", "Young", age)
  cellType <- gsub("_Old.*|_Young.*", "", fileName)
  cellType <- gsub(".*\\.", "", cellType)
  tissue <- basename(bamPath)
  label <- paste0(tissue, ":", cellType)
  input <- data.frame(cellType = cellType, age = age, sample = gsub("mouse.*", "rep", fileName), histone = gsub("\\..*", "", fileName), file = bamFiles)
  input$sample <- paste0(label, "_", input$age)
  input$reads <- reads[input$file]
  input <- arrange(input, cellType, histone, age, desc(reads)) # rank by read depths
  input <- input %>% group_by(sample, histone) %>% mutate(sample = paste0(sample, "_rep_", row_number())) %>% ungroup()
  cat("Select the top two replicates with the highest sequencing depths for each age group.\n")
  input <- input[grep("_rep_1|_rep_2", input$sample),] # only keep the first two replicates

  for(i in unique(cellType)){
    input.subset <- input[input$cellType == i,]
    paths <- file.path(outPath, tissue, i, "Old_Young")
    dir.create(paths, recursive = TRUE, showWarnings = F)
    write.table(input.subset[,c("sample", "histone", "file")], file.path(paths, "cellmarkfile.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

    input.subset <- input[input$cellType == i & input$age == "Old",]
    paths <- file.path(outPath, tissue, i, "Old")
    dir.create(paths, recursive = TRUE, showWarnings = F)
    write.table(input.subset[,c("sample", "histone", "file")], file.path(paths, "cellmarkfile.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

    input.subset <- input[input$cellType == i & input$age == "Young",]
    paths <- file.path(outPath, tissue, i, "Young")
    dir.create(paths, recursive = TRUE, showWarnings = F)
    write.table(input.subset[,c("sample", "histone", "file")], file.path(paths, "cellmarkfile.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
  }
  saveRDS(input, file.path(outPath, tissue, "input.rds"))
  return(input)
}