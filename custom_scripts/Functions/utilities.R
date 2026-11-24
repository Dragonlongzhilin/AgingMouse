#' @description convert p-value to stars
p_to_stars <- function(p_value) {
  if (is.na(p_value)) {
    return(NA)
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

#' @description count the frequency across tissues
count.freq <- function(genes, sig.list){
  Freq <- lapply(genes, function(gene){
    flags <- lapply(sig.list, function(x){
      if(gene %in% x){
        return(T)
      }else{
        return(F)
      }
    })
    flags <- names(flags)[which(unlist(flags))]
    # tissue <- str_extract(flags, "^[^\\.]+")
    tissue <- gsub(":.*", "", flags)
    info <- as.data.frame(table(tissue))
    info$tissue <- as.character(info$tissue)
    info$gene <- gene
    info[nrow(info)+1,] <- c("Merge", sum(info$Freq), gene)
    return(info)
  })
  Freq <- Reduce(function(x,y) rbind(x,y), Freq)
  Freq$Freq <- as.numeric(Freq$Freq)
  return(Freq)
}

#' @description count the frequency across tissues + cellType, and category it
count.freq2 <- function(genes, sig.list, assay.info){
  require(stringr)
  require(dplyr)
  Freq <- lapply(genes, function(gene){
    flags <- lapply(sig.list, function(x){
      if(gene %in% x){
        return(T)
      }else{
        return(F)
      }
    })
    flags <- names(flags)[which(unlist(flags))]
    flags <- gsub(".Up$|.Down$", "", flags)
    flags <- gsub(".up$|.down$", "", flags)
    flags <- gsub("\\.", ":", flags)
    # tissue + cellType
    info <- data.frame(gene = gene, tissue_cellType = flags, Freq = length(flags))
    res <- info %>% dplyr::left_join(unique(assay.info), by = "tissue_cellType")
    return(res)
  })
  Freq <- Reduce(function(x,y) rbind(x,y), Freq)

  # count the total frequency
  merge.Freq <- unique(Freq[,c("gene", "Freq")])
  merge.Freq <- merge.Freq[order(merge.Freq$Freq, decreasing = T),]
  Freq$gene <- factor(Freq$gene, levels = merge.Freq$gene)
  Freq <- dplyr::arrange(Freq, gene)
  return(Freq)
}

#' @description plot the 0-1 heatmap of interested/top30 high frequency genes/bins across tissues + cell type
Freq.plot <- function(Freq.data, color.tissue, color.category, assay.info, top = 30, select.genes = NULL, title = "", hgt = 15, wth = NULL, colorLegend = NULL, select.category = NULL, addMissing = FALSE, row_split = NULL, col.order = NULL, row.order = NULL){
  require(ComplexHeatmap)
  require(circlize)
  require(tidyr)
  Freq.data$Freq <- 1

  if(is.null(select.category)){
    order.genes <- sort(table(Freq.data$gene), decreasing = T)
  }else{
    order.genes <- sort(table(Freq.data$gene[which(Freq.data$category %in% select.category)]), decreasing = T)
  }

  if(is.null(select.genes)){
    genes <- names(order.genes)[1:top]
  }else{
    idx <- match(select.genes, names(order.genes))
    genes <- order.genes[na.omit(idx)]
    genes <- genes[order(genes, decreasing = T)]
    genes <- names(genes)
  }
  genes <- na.omit(genes)

  if(!is.null(row.order)){
    genes <- genes[na.omit(match(row.order, genes))]
  }

  plot.data <- Freq.data[which(Freq.data$gene %in% genes),]
  if(!is.null(colorLegend)){
    idx <- match(colorLegend, colnames(plot.data))
    colnames(plot.data)[idx] <- "P.value"
    mat_table <- plot.data %>% select(gene, P.value, tissue_cellType) %>% pivot_wider(names_from = tissue_cellType, values_from = P.value, values_fill = NA) %>% as.data.frame()
    rownames(mat_table) <- mat_table$gene
    mat_table <- mat_table[,-1]
    mat_table[mat_table == "NULL"] <- NA
    a <- apply(mat_table, 2, function(x){as.numeric(as.character(x))})
    rownames(a) <- rownames(mat_table)
    mat_table <- a
  }else{
    mat_table <- as.data.frame.matrix(table(plot.data$gene, plot.data$tissue_cellType))
    mat_table[mat_table > 0] <- 1
  }

  mat_table <- mat_table[genes,]
  idx <- match(colnames(mat_table), assay.info$tissue_cellType)
  if(addMissing){
    na.idx <- setdiff(1:nrow(assay.info), idx)
    na_table <- matrix(NA, nrow = nrow(mat_table), ncol = length(na.idx))
    rownames(na_table) <- rownames(mat_table)
    colnames(na_table) <- assay.info$tissue_cellType[na.idx]
    mat_table <- cbind(mat_table, na_table)
  }

  # rank the column
  if(!is.null(col.order)){
    mat_table <- mat_table[,col.order]
  }else{
    rank.mat <- mat_table
    rank.mat[is.na(rank.mat)] <- 0
    rank.order <- colMeans(rank.mat)
    rank.order <- sort(rank.order, decreasing = T)
    mat_table <- mat_table[,names(rank.order)]
  }

  # add information
  idx <- match(colnames(mat_table), assay.info$tissue_cellType)
  phenotype.info <- assay.info[idx,]
  depth <- round(phenotype.info$read, 2)
  col_ha <- HeatmapAnnotation(tissue = phenotype.info$tissue, category = phenotype.info$category, read = anno_barplot(depth, bar_width = 1, gp = gpar(col = "white", fill = "#35B199"), height = unit(2, "cm"), add_numbers = T, numbers_rot = 90), col = list(tissue = color.tissue, category = color.category))
  bar.info <- data.frame(gene = genes, Freq = as.numeric(order.genes[genes]))
  bar.info <- na.omit(bar.info)
  bar.info <- bar.info[match(genes, bar.info$gene),]

  if(length(unique(bar.info$Freq)) == 1){
    left_annotation <- NULL
  }else{
    left_annotation <- rowAnnotation(Freq = anno_numeric(bar.info$Freq, align_to = "right", bg_gp = gpar(fill = "grey", col = "grey")), annotation_name_rot = 0)
  }

  if(!is.null(colorLegend)){
    cols <- colorRamp2(c(0.05, min(na.omit(as.numeric(mat_table)))), c("#E7C76F", "#C30812"))
    name <- colorLegend
  }else{
    mat_table_chr <- apply(mat_table, 2, as.character)
    colnames(mat_table_chr) <- colnames(mat_table)
    rownames(mat_table_chr) <- rownames(mat_table)
    mat_table <- mat_table_chr
    cols <- c("0" = "#BEBEBE", "1" = "#C30812")
    name <- "Appearance"
  }

  if(is.null(wth)){
    wth <- ncol(mat_table) * 0.5
  }
  ## split category
  column_split <- phenotype.info$category

  if(!is.null(row_split)){
    idx <- match(rownames(mat_table), row_split$Gene)
    row_split <- row_split$GroupSet[idx]
    row_gap = unit(2, "mm")
  }else{
    row_gap = unit(1, "mm")
  }
  p <- Heatmap(mat_table, name = name, col = cols, column_title = title, cluster_rows = F, cluster_columns = F, column_names_rot = 60, show_row_dend = F, column_split = column_split, left_annotation = left_annotation, top_annotation = col_ha, height = unit(hgt, "cm"), width = unit(wth, "cm"), row_gap = row_gap, rect_gp = gpar(col = "white", lwd = 2), row_split = row_split)
  return(list(Freq.top = bar.info, p = p))
}

#' @description plot the logFC heatmap of select high frequency genes/bins across tissues + cell type
## edgeR: LFC & fdr/PValue
edgeR.LFC.heatmap <- function(edgeR.result.list, select.genes, assay.info, color.tissue, color.category, Freq.mat = NULL, title = "", hgt = 12, wth = 70, select.category = NULL, cluster_rows = TRUE, plot = TRUE, fontsize = 12, setNA = NULL, sig.label = "fdr"){
  require(ComplexHeatmap)
  require(circlize)
  require(dplyr)
  if(!is.list(edgeR.result.list[[1]][[1]])){
    LFC_FDR.list <- lapply(names(edgeR.result.list), function(tissue_cellType){
      y <- edgeR.result.list[[tissue_cellType]]
      idx <- match(select.genes, y$gene)
      if(all(is.na(idx))){
        LFC_FDR <- data.frame(tissue_cellType = tissue_cellType, tissue = gsub(":.*", "", tissue_cellType), cellType = gsub(".*:", "", tissue_cellType), gene = select.genes, LFC = NA, fdr = NA)
      }else{
        LFC_FDR <- data.frame(tissue_cellType = tissue_cellType, tissue = gsub(":.*", "", tissue_cellType), cellType = gsub(".*:", "", tissue_cellType), gene = select.genes, LFC = y$logFC[idx], fdr = y[idx, sig.label])
      }
      return(LFC_FDR)
    })
    LFC_FDR.list <- bind_rows(LFC_FDR.list)
  }else{
    LFC_FDR.list <- lapply(names(edgeR.result.list), function(tissue){
      edgeR.res <- edgeR.result.list[[tissue]]
      LFC_FDR.res <- lapply(names(edgeR.res), function(cellType){
        y <- edgeR.res[[cellType]]
        idx <- match(select.genes, y$gene)
        if(all(is.na(idx))){
          LFC_FDR <- data.frame(tissue = tissue, cellType = gsub("-| ", "_", cellType), gene = select.genes, LFC = NA, fdr = NA)
        }else{
          LFC_FDR <- data.frame(tissue = tissue, cellType = gsub("-| ", "_", cellType), gene = select.genes, LFC = y$logFC[idx], fdr = y[idx, sig.label])
        }
        return(LFC_FDR)
      })
      LFC_FDR.res <- bind_rows(LFC_FDR.res)
      return(LFC_FDR.res)
    })
    LFC_FDR.list <- bind_rows(LFC_FDR.list)
    LFC_FDR.list$tissue_cellType <- paste0(LFC_FDR.list$tissue, ":", LFC_FDR.list$cellType)
  }
  LFC.matrix <- tidyr::spread(LFC_FDR.list[,c("gene", "tissue_cellType", "LFC")], tissue_cellType, LFC)
  rownames(LFC.matrix) <- LFC.matrix[,1]
  LFC.matrix <- LFC.matrix[,-1]

  FDR.matrix <- tidyr::spread(LFC_FDR.list[,c("gene", "tissue_cellType", "fdr")], tissue_cellType, fdr)
  rownames(FDR.matrix) <- FDR.matrix[,1]
  FDR.matrix <- FDR.matrix[,-1]
  FDR.matrix[is.na(FDR.matrix)] <- 1

  info.subset <- data.frame(tissue_cellType = colnames(LFC.matrix))
  info.subset <- info.subset %>% left_join(unique(assay.info), by = "tissue_cellType")
  info.subset <- unique(arrange(info.subset, category, tissue_cellType))

  # rank by mean logFC
  if(!is.null(setNA)){
    cat("set edgeR LFC NA --> ", setNA, ", use to calculate mean logFC\n")
    a <- LFC.matrix
    a[is.na(a)] <- setNA
    rank_order <- colMeans(a, na.rm = T)
  }else{
    rank_order <- colMeans(LFC.matrix, na.rm = T)
  }
  rank_order <- sort(rank_order, decreasing = T)
  info.subset <- info.subset[match(names(rank_order), info.subset$tissue_cellType),]

  if(!is.null(select.category)){
    info.subset <- info.subset[info.subset$category %in% select.category,]
  }

  plot.data <- LFC.matrix[select.genes, info.subset$tissue_cellType]
  plot.data2 <- FDR.matrix[select.genes, info.subset$tissue_cellType]

  if(plot){
    if(!is.null(Freq.mat)){
      Freq.mat <- as.data.frame(Freq.mat)
      Freq.mat <- Freq.mat[match(select.genes, Freq.mat[,1]),]
      left_annotation <- rowAnnotation(Freq = anno_numeric(as.numeric(Freq.mat[,2]), align_to = "right", bg_gp = gpar(fill = "grey", col = "grey")), annotation_name_rot = 0)
    }else{
      left_annotation <- NULL
    }
    if(is.null(wth)){
      wth <- ncol(plot.data) * 0.5
    }

    column_ha <- HeatmapAnnotation(category = info.subset$category, tissue = info.subset$tissue, read = anno_barplot(as.numeric(info.subset$read), bar_width = 1, gp = gpar(col = "white", fill = "#35B199"), height = unit(2, "cm"), add_numbers = T, numbers_rot = 90), col = list(category = color.category, tissue = color.tissue))
    col_fun <- colorRamp2(c(quantile(na.omit(as.numeric(unlist(plot.data))), 0.001), 0, quantile(na.omit(as.numeric(unlist(plot.data))), 0.999)), c("blue", "white", "red"))
    p <- Heatmap(plot.data, col = col_fun, name = "edgeR:LFC", column_title = title, top_annotation = column_ha, show_column_dend = F, show_row_dend = F, cluster_rows = cluster_rows, left_annotation = left_annotation, clustering_method_rows = "average", cluster_columns = F, column_names_gp = gpar(fontsize = fontsize), row_names_gp = gpar(fontsize = fontsize), column_title_gp = gpar(fontsize = floor(fontsize * 1.2)), clustering_method_columns = "average", width = unit(wth, "cm"), height = unit(hgt, "cm"), column_names_rot = 60, cell_fun = function(j, i, x, y, width, height, fill) {if(!is.na(plot.data2[i, j]) && plot.data2[i, j] < 0.05) grid.text("*", x, y, gp = gpar(fontsize = 8, fontface = "bold"))})
    draw(p)

    column_split <- info.subset$category
    p <- Heatmap(plot.data, col = col_fun, name = "edgeR:LFC", column_title = title, top_annotation = column_ha, show_column_dend = F, show_row_dend = F, cluster_rows = cluster_rows, left_annotation = left_annotation, cluster_columns = F, clustering_method_rows = "average", column_names_gp = gpar(fontsize = fontsize), row_names_gp = gpar(fontsize = fontsize), column_title_gp = gpar(fontsize = floor(fontsize * 1.2)), column_split = column_split, width = unit(wth, "cm"), height = unit(hgt, "cm"), column_names_rot = 60, cell_fun = function(j, i, x, y, width, height, fill) {if(!is.na(plot.data2[i, j]) && plot.data2[i, j] < 0.05) grid.text("*", x, y, gp = gpar(fontsize = 8, fontface = "bold"))})
    draw(p)
  }

  return(list(LFC = plot.data, FDR = plot.data2, info = info.subset))
}

#' @description plot the score heatmap of select pathways across tissues + cell type
## e.g., GSEA: NES & FDR
Pathway.Heatmap <- function(PathwayMatrix, select.pathways, assay.info, select.cols = c("Description", "NES", "p.adjust"), color.tissue, color.category, name = "GSEA:NES", Freq.mat = NULL, title = "", hgt = 12, wth = 70, select.category = NULL, cluster_rows = TRUE, plot = TRUE, show_row_names = TRUE, fontsize = 12, row_names_side = "left", setNA = NULL, cell_fun = FALSE, addMissing = FALSE){
  require(ComplexHeatmap)
  require(circlize)
  require(dplyr)

  if(addMissing){
    select.pathways.info <- lapply(names(PathwayMatrix), function(y){
      x <- PathwayMatrix[[y]]
      target_pathways <- select.pathways
      if(is.null(x)){
        # Create a dataframe with all target pathways, score=0, p.adjust=1
        info.sig <- data.frame(Description = target_pathways, score = setNA, p.adjust = 1, tissue_cellType = y)
        return(info.sig)
      }else{
        # Find which target pathways exist in x
        existing_idx <- which(x[,select.cols[1]] %in% target_pathways)
        if(length(existing_idx) == 0){
          # No pathways found - create dataframe with all target pathways
          info.sig <- data.frame(Description = target_pathways, score = setNA, p.adjust = 1, tissue_cellType = y)
          return(info.sig)
        }else{
          # Get info for existing pathways
          info.sig <- x[existing_idx, select.cols]
          colnames(info.sig) <- c("Description", "score", "p.adjust")
          info.sig$tissue_cellType <- y
          
          # Check for missing pathways and add them with score=0, p.adjust=1
          found_pathways <- info.sig$Description
          missing_pathways <- setdiff(target_pathways, found_pathways)
          if(length(missing_pathways) > 0){
            missing_df <- data.frame(Description = missing_pathways,score = setNA, p.adjust = 1, tissue_cellType = y
            )
            info.sig <- rbind(info.sig, missing_df)
          }
          return(info.sig)
        }
      }
    })
  }else{
    select.pathways.info <- lapply(names(PathwayMatrix), function(y){
      x <- PathwayMatrix[[y]]
      if(is.null(x)){
        return(NULL)
      }else{
        info.sig <- x[which(x[,select.cols[1]] %in% select.pathways), select.cols]
        if(nrow(info.sig) == 0){
          return(NULL)
        }else{
          colnames(info.sig) <- c("Description", "score", "p.adjust")
          info.sig$tissue_cellType <- y
          return(info.sig)
        }
      }
    })
  }
  select.pathways.info <- select.pathways.info[!sapply(select.pathways.info,is.null)]
  select.pathways.info <- do.call("rbind", select.pathways.info)
  select.pathways.info$tissue_cellType <- gsub("\\.", ":", select.pathways.info$tissue_cellType)
  select.pathways.info$tissue_cellType <- gsub("-| ", "_", select.pathways.info$tissue_cellType)
  score.data <- select.pathways.info %>%  dplyr::select(Description, score, tissue_cellType) %>% pivot_wider(names_from = tissue_cellType, values_from = score, values_fill = NA)
  FDR.data <- select.pathways.info %>%  dplyr::select(Description, p.adjust, tissue_cellType) %>% pivot_wider(names_from = tissue_cellType, values_from = p.adjust, values_fill = NA)
  score.data <- as.data.frame(score.data)
  FDR.data <- as.data.frame(FDR.data)
  rownames(score.data) <- score.data$Description
  score.data <- score.data[,-1]
  score.data <- as.matrix(score.data)
  rownames(FDR.data) <- FDR.data$Description
  FDR.data <- FDR.data[,-1]
  FDR.data <- as.matrix(FDR.data)
  FDR.data[is.na(FDR.data)] <- 1

  info.subset <- data.frame(tissue_cellType = colnames(score.data))
  info.subset <- info.subset %>% left_join(unique(assay.info), by = "tissue_cellType")
  info.subset <- unique(arrange(info.subset, category, tissue_cellType))

  if(cluster_rows){
    cat("set NA --> ", setNA, ", use to cluster and calculate mean Rank\n")
    score.data[is.na(score.data)] <- setNA
    rank_order <- colMeans(score.data, na.rm = T)
  }else{
    # rank by mean Rank
    if(!is.null(setNA)){
      cat("set NA --> ", setNA, ", use to calculate mean Rank\n")
      a <- score.data
      a[is.na(a)] <- setNA
      rank_order <- colMeans(a, na.rm = T)
    }else{
      rank_order <- colMeans(score.data, na.rm = T)
    }
  }
  rank_order <- sort(rank_order, decreasing = T)
  info.subset <- info.subset[match(names(rank_order), info.subset$tissue_cellType),]

  if(!is.null(select.category)){
    info.subset <- info.subset[info.subset$category %in% select.category,]
  }

  plot.data <- score.data[select.pathways, info.subset$tissue_cellType]
  plot.data2 <- FDR.data[select.pathways, info.subset$tissue_cellType]

  if(plot){
    if(!is.null(Freq.mat)){
      Freq.mat <- as.data.frame(Freq.mat)
      Freq.mat <- Freq.mat[match(select.pathways, Freq.mat[,1]),]
      left_annotation <- rowAnnotation(Freq = anno_numeric(as.numeric(Freq.mat[,2]), align_to = "left", bg_gp = gpar(fill = "grey", col = "grey")), annotation_name_rot = 0)
    }else{
      left_annotation <- NULL
    }
    if(is.null(wth)){
      wth <- ncol(plot.data) * 0.5
    }

    idx <- grep("PeakNumber", colnames(info.subset))
    if(length(idx) > 0){
      column_ha <- HeatmapAnnotation(category = info.subset$category, tissue = info.subset$tissue, read = anno_barplot(as.numeric(info.subset$read), bar_width = 1, gp = gpar(col = "white", fill = "#35B199"), height = unit(2, "cm"), add_numbers = T, numbers_rot = 90), PeakNumber = anno_barplot(as.numeric(info.subset$PeakNumber), bar_width = 1, gp = gpar(col = "white", fill = "#D3765A"), height = unit(2, "cm"), add_numbers = T, numbers_rot = 90), col = list(category = color.category, tissue = color.tissue))
    }else{
      column_ha <- HeatmapAnnotation(category = info.subset$category, tissue = info.subset$tissue, read = anno_barplot(as.numeric(info.subset$read), bar_width = 1, gp = gpar(col = "white", fill = "#35B199"), height = unit(2, "cm"), add_numbers = T, numbers_rot = 90), col = list(category = color.category, tissue = color.tissue))
    }

    scores <- unique(na.omit(as.numeric(unlist(plot.data))))
    if(all(scores >= 0, na.rm = TRUE)){
      col_fun <- colorRamp2(c(0, quantile(scores, 0.9)), c("white", "red"))
    }else if(all(scores < 0, na.rm = TRUE)){
      col_fun <- colorRamp2(c(quantile(scores, 0.1), 0), c("blue", "white"))
    }else{
      col_fun <- colorRamp2(c(quantile(scores, 0.01), 0, quantile(scores, 0.99)), c("blue", "white", "red"))
    }

    if(cell_fun){
      p <- Heatmap(plot.data, col = col_fun, name = name, column_title = title, top_annotation = column_ha, show_column_dend = F, show_row_dend = F, show_row_names = show_row_names, cluster_rows = cluster_rows, right_annotation = left_annotation, clustering_method_rows = "average", row_names_side = row_names_side, cluster_columns = F, column_names_gp = gpar(fontsize = fontsize), row_names_gp = gpar(fontsize = fontsize), column_title_gp = gpar(fontsize = floor(fontsize * 1.2)), clustering_method_columns = "average", width = unit(wth, "cm"), height = unit(hgt, "cm"), column_names_rot = 60, cell_fun = function(j, i, x, y, width, height, fill) {if(!is.na(plot.data2[i, j]) && plot.data2[i, j] < 0.05) grid.text("*", x, y, gp = gpar(fontsize = 8, fontface = "bold"))})
      draw(p)

      column_split <- info.subset$category
      p <- Heatmap(plot.data, col = col_fun, name = name, column_title = title, top_annotation = column_ha, show_column_dend = F, show_row_dend = F, show_row_names = show_row_names, cluster_rows = cluster_rows, right_annotation = left_annotation, cluster_columns = F, clustering_method_rows = "average", column_names_gp = gpar(fontsize = fontsize), row_names_gp = gpar(fontsize = fontsize), row_names_side = row_names_side, column_title_gp = gpar(fontsize = floor(fontsize * 1.2)), column_split = column_split, width = unit(wth, "cm"), height = unit(hgt, "cm"), column_names_rot = 60, cell_fun = function(j, i, x, y, width, height, fill) {if(!is.na(plot.data2[i, j]) && plot.data2[i, j] < 0.05) grid.text("*", x, y, gp = gpar(fontsize = 8, fontface = "bold"))})
      draw(p)
    }else{
      p <- Heatmap(plot.data, col = col_fun, name = name, column_title = title, top_annotation = column_ha, show_column_dend = F, show_row_dend = F, show_row_names = show_row_names, cluster_rows = cluster_rows, right_annotation = left_annotation, clustering_method_rows = "average", row_names_side = row_names_side, cluster_columns = F, column_names_gp = gpar(fontsize = fontsize), row_names_gp = gpar(fontsize = fontsize), column_title_gp = gpar(fontsize = floor(fontsize * 1.2)), clustering_method_columns = "average", width = unit(wth, "cm"), height = unit(hgt, "cm"), column_names_rot = 60)
      draw(p)

      column_split <- info.subset$category
      p <- Heatmap(plot.data, col = col_fun, name = name, column_title = title, top_annotation = column_ha, show_column_dend = F, show_row_dend = F, show_row_names = show_row_names, cluster_rows = cluster_rows, right_annotation = left_annotation, cluster_columns = F, clustering_method_rows = "average", column_names_gp = gpar(fontsize = fontsize), row_names_gp = gpar(fontsize = fontsize), row_names_side = row_names_side, column_title_gp = gpar(fontsize = floor(fontsize * 1.2)), column_split = column_split, width = unit(wth, "cm"), height = unit(hgt, "cm"), column_names_rot = 60)
      draw(p)
    }
  }
  return(list(score = plot.data, FDR = plot.data2, info = info.subset))
}

#' @description group the TFs across different category
TF.Group <- function(TF.list){
  df <- TF.list %>% enframe(name = "Group", value = "Gene") %>% unnest(cols = c(Gene)) %>% distinct()

  # 查看每个基因出现在哪些组
  gene_groups <- df %>% group_by(Gene) %>% summarise(Groups = list(unique(Group)), .groups = "drop")

  # 按“出现在相同分组集合”归类
  gene_clusters <- gene_groups %>% mutate(GroupSet = sapply(Groups, function(x) paste(sort(x), collapse = "|"))) %>% group_by(GroupSet) %>% mutate(Cluster = cur_group_id()) %>% ungroup()

  # 查看每个簇（基因组合）对应的分组集合
  gene_clusters %>% select(GroupSet, Cluster, Gene) %>% group_by(Cluster, GroupSet) %>% summarise(Genes = paste(Gene, collapse = ", "), .groups = "drop")
  gene_clusters <- as.data.frame(gene_clusters)
  return(gene_clusters)
}
