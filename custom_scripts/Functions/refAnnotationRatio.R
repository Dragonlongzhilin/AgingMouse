#' @description plot the heatmap and cluster ratio for reference annotation
#' @param PairedTag_data
#' @return
#' @export

refAnnotationRatio <- function(meta.data, group.by = "seurat_clusters", predictedLabels, width = unit(6, "cm"), height = unit(8, "cm"), font_size = 10, plot = T, row_title = character(0), Unpredictable_cellType = NULL){
  require(tidyverse)
  require(ComplexHeatmap)
  require(circlize)
  assign.ratio <- lapply(predictedLabels, function(x){
    df <- data.frame(cluster = meta.data[,group.by], annotation = meta.data[,x])
    most_frequent_annotations <- df %>%
      group_by(cluster, annotation) %>%
      summarise(count = n()) %>%
      arrange(cluster, desc(count)) %>%
      group_by(cluster) %>%
      slice_max(count, n = 1, with_ties = FALSE) %>%
      ungroup()
    most_frequent_annotations <- as.data.frame(most_frequent_annotations)
    counts <- as.data.frame(table(df$cluster))
    most_frequent_annotations$ratio <- round(most_frequent_annotations$count / counts$Freq, 4)

    # ratio < 0.5 --- mix
    most_frequent_annotations$annotationHQ <- most_frequent_annotations$annotation
    idx <- which(most_frequent_annotations$ratio < 0.5)
    if(length(idx) > 0){
      most_frequent_annotations$annotationHQ[idx] <- paste0("Mix_", 1:length(idx))
    }
    colnames(most_frequent_annotations)[2:5] <- paste0(x, ".", colnames(most_frequent_annotations)[2:5])

    # plot heatmap
    col_fun <- colorRamp2(c(0, 0.5, 1), c("#441C5A", "white", "#F1E91D"))
    plot.data <- as.data.frame(table(df))
    plot.data$Total <- rep(counts$Freq, times = length(unique(df$annotation)))
    plot.data$Ratio <- plot.data$Freq / plot.data$Total
    plot.data <- spread(plot.data[,c(1,2,5)], annotation, Ratio)
    rownames(plot.data) <- plot.data$cluster
    plot.data <- plot.data[,-1]

    # reorder columns
    power <- 10  # >1 emphasizes peak values; higher = sharper band (tune 5~20)
    col_scores <- apply(plot.data, 2, function(x) {
      weighted.mean(1:nrow(plot.data), w = x^power, na.rm = TRUE)
    })
    col_levels <- colnames(plot.data)[order(col_scores)]
    plot.data <- plot.data[, col_levels]

    # remove Unpredictable_cellType
    if(!is.null(Unpredictable_cellType)){
      idx <- which(rownames(plot.data) %in% Unpredictable_cellType)
      plot.data[idx,] <- NA
    }

    p <- Heatmap(plot.data, name = "Ratio", width = width, height = height, column_title = x, row_title = row_title, col = col_fun, show_row_dend = F, show_column_dend = F, cluster_rows = F, cluster_columns = F, row_names_gp = gpar(fontsize = font_size), column_names_gp = gpar(fontsize = font_size), cell_fun = function(j, i, x, y, width, height, fill) {if(!is.na(plot.data[i, j]) && plot.data[i, j] > 0.5) grid.text("*", x, y, gp = gpar(fontsize = 10))})

    if(plot){
      print(p)
    }
    return(most_frequent_annotations)
  })
  assign.ratio <- Reduce(function(x,y) merge(x,y,by = "cluster"), assign.ratio)
  return(assign.ratio)
}