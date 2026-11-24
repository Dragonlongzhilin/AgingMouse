#' @description make outfile for extract_sample_barcode_from_bam.py
#' @Author: Zhilin Long
#' @date 2025-3-29

#' @param PairedTag
#' @param workPath: public html
#' @param group_name: age/cellType
#' @param split.by: histone

makeGroupForPseudoBulk <- function(meta.data, workPath, group_name, split.by = NULL){
  # meta.data <- PairedTag@meta.data
  if(!is.null(split.by)){
    split.group <- unique(meta.data[,split.by])
    for(i in split.group){
      pseudo_bulk.path <- file.path(workPath, i, group_name) 
      if (!dir.exists(pseudo_bulk.path)) {
        dir.create(pseudo_bulk.path, recursive = TRUE)
      }
      meta.data.subset <- meta.data[which(meta.data[,split.by] == i),]
      cat("split ", i, " in ", group_name, ": ", nrow(meta.data.subset), " cells\n")
      rownames(meta.data.subset) <- gsub(".*_", "", rownames(meta.data.subset))
      group <- data.frame(barcode = rownames(meta.data.subset), sample = meta.data.subset[,group_name])
      group$sample <- gsub(" |-", "_", group$sample)
      write.table(group, file.path(pseudo_bulk.path, paste0(group_name, ".txt")), row.names = F, col.names = T, quote = F, sep = "\t")
    }
  }else{
    pseudo_bulk.path <- file.path(workPath, group_name) 
    if (!dir.exists(pseudo_bulk.path)) {
      dir.create(pseudo_bulk.path, recursive = TRUE)
    }
    rownames(meta.data) <- gsub(".*_", "", rownames(meta.data))
    cat(group_name, ": ", nrow(meta.data), " cells\n")
    group <- data.frame(barcode = rownames(meta.data), sample = meta.data[,group_name])
    group$sample <- gsub(" |-", "_", group$sample)
    write.table(group, file.path(pseudo_bulk.path, paste0(group_name, ".txt")), row.names = F, col.names = T, quote = F, sep = "\t")
  }
  return(group)
}
