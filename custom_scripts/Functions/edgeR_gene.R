#' @description pseudobulk groupserential analysis with edgeR method
#' @Author: Zhilin Long
#' @date 2025-4-29

#' @param group factor, e.g., c("d5", "d30", "d50", "d70")

edgeR_gene <- function(counts, group, label, batch = NULL, use.PValue = FALSE, top.n = 5, TransposableElement = NULL, ContinuousComparison = FALSE, cpm.threshold = 1, cpm.filter = FALSE, sample.label = NULL){
  require(edgeR)
  require(EnhancedVolcano)
  require(ggpubr)
  require(dplyr)

  if(!is.factor(group)){
    group <- as.factor(group)
  }
  lvls <- levels(group)
  if(!is.null(batch)){
    if (!is.factor(batch)) {
      batch <- as.factor(batch)
    }
    design <- model.matrix(~0+group+batch)
    colnames(design)[1:length(lvls)] <- lvls
  }else{
    design <- model.matrix(~0+group)
    colnames(design) <- lvls
  }
  if(is.null(batch)){
    cat("Batch: none\n")
  }else{
    cat("Batch: ", paste(levels(batch), collapse = ", "), "\n")
  }
  cpms <- edgeR::cpm(counts)

  cat(nrow(counts), " genes in raw matrix\n")
  y <- DGEList(counts = counts, group = group)
  if(cpm.filter){
    # cat("Filtering method: cpm >= ", cpm.threshold, " in at least 2 samples in each group\n")
    # # cpm >= 1 in at least 2 samples in each group
    # group_counts <- lapply(levels(group), function(x){
    #   idx <- which(group == x)
    #   count_gt1 <- apply(cpms[, idx], 1, function(x) sum(x >= cpm.threshold, na.rm = TRUE))
    #   return(which(count_gt1 >= 2))
    # })
    # keep <- unique(unlist(group_counts))
    # genes.left <- length(keep)
    # filter.ratio <- genes.left / nrow(counts)

    # cpm >= 1 in at least 2 samples in each group
    cat("Filtering method: cpm >= ", cpm.threshold, " in at least 2 samples\n")
    count_gt <- apply(cpms, 1, function(x) sum(x >= cpm.threshold, na.rm = TRUE))
    keep <- which(count_gt >= 2)
    genes.left <- length(keep)
  }else{
    cat("Filtering method: edgeR default\n")
    keep <- filterByExpr(y)
    genes.left <- length(which(keep))
  }
  filter.ratio <- genes.left / nrow(counts)
  cat(genes.left, " genes in filted matrix: ", round(filter.ratio, 2), "\n")

  y <- y[keep,, keep.lib.sizes=F]
  y <-  calcNormFactors(y)
  my_colors <- RColorBrewer::brewer.pal(n = nlevels(group), name = "Set1")
  if(!is.null(sample.label)){
    plotMDS(y, col = my_colors[as.numeric(group)], pch = 16, cex = 1)
    mds_data <- plotMDS(y, col = my_colors[as.numeric(group)], cex = 0.8, plot = FALSE)
    text(x = mds_data$x, y = mds_data$y, labels = sample.label, pos = 3, cex = 0.8, col = my_colors[as.numeric(group)])
    legend("topright", legend = levels(group), col = my_colors, pch = 16)
  }else{
    plotMDS(y, col = my_colors[as.numeric(group)], pch = 16, cex = 1)
    mds_data <- plotMDS(y, col = my_colors[as.numeric(group)], cex = 0.8, plot = FALSE)
    text(x = mds_data$x, y = mds_data$y, labels = colnames(counts), pos = 3, cex = 0.8, col = my_colors[as.numeric(group)])
    legend("topright", legend = levels(group), col = my_colors, pch = 16)
  }
  y <- estimateDisp(y, design, robust=TRUE)
  bcv <- sqrt(y$tagwise.dispersion)
  plotBCV(y)
  fit_tag <- glmFit(y, design)

  # Continuous comparison
  if(ContinuousComparison){
    cat("edgeR: continuous comparison \n")
    len <- length(lvls)
    myContrasts <- paste(lvls[2:len], lvls[1:(len-1)], sep = '-')
    contrast.matrix <- makeContrasts(contrasts = myContrasts, levels = design)
    lrt <- glmLRT(fit_tag, contrast = contrast.matrix)
    fdr <- p.adjust(lrt$table$PValue, method = "fdr")
    ContinuousComparison <- cbind(cpms[rownames(lrt$table),], lrt$table, fdr)
    ContinuousComparison <- as.data.frame(ContinuousComparison)
    ContinuousComparison$gene <- rownames(lrt$table)
    ContinuousComparison$bcv <- bcv
    ContinuousComparison <- ContinuousComparison[order(ContinuousComparison$fdr),]
  }else{
    ContinuousComparison <- NULL
  }

  # group comparison
  combinations <- combn(levels(group), 2)
  myContrasts <- paste0(combinations[2,], "-", combinations[1,])
  contrast.matrix <- makeContrasts(contrasts = myContrasts, levels = design)
  GroupComparison <- lapply(1:length(myContrasts), function(i){
    cat("edgeR: group comparison (", myContrasts[i], ")\n")
    lrt <- glmLRT(fit_tag, contrast = contrast.matrix[,i])
    plotMD(lrt)
    fdr <- p.adjust(lrt$table$PValue, method = "fdr")
    out <- cbind(cpms[rownames(lrt$table),], lrt$table, fdr)
    out <- as.data.frame(out)
    out$gene <- rownames(lrt$table)
    out$bcv <- bcv
    out <- out[order(out$fdr),]
    titles <- gsub("-", ".vs.", myContrasts[i])

    plot.data <- out
    if(!is.null(TransposableElement)){
      plot.data <- plot.data[plot.data$gene %in% TransposableElement,]
    }
    plot.data$groups <- "NS"

    if(use.PValue){
      plot.data$groups[plot.data$logFC > 0 & plot.data$PValue < 0.05] <- "UP"
      plot.data$groups[plot.data$logFC < 0 & plot.data$PValue < 0.05] <- "DOWN"
      y.label <- "PValue"
    }else{
      plot.data$groups[plot.data$logFC > 0 & plot.data$fdr < 0.05] <- "UP"
      plot.data$groups[plot.data$logFC < 0 & plot.data$fdr < 0.05] <- "DOWN"
      y.label <- "fdr"
    }

    up.num <- sum(plot.data$groups == "UP")
    down.num <- sum(plot.data$groups == "DOWN")
    keyvals <- ifelse(
      plot.data$groups == 'UP', 'red3',
      ifelse(plot.data$groups == 'DOWN', 'royalblue', 'grey30')
    )
    names(keyvals) <- plot.data$groups

    if(top.n > 0){
      top.genes <- plot.data %>% filter(groups %in% c("UP", "DOWN")) %>% group_by(groups) %>% arrange(fdr, .by_group = TRUE) %>% slice_head(n = top.n) %>% pull(gene)
      lab <- rep(NA, nrow(plot.data))
      lab[match(top.genes, plot.data$gene)] <- top.genes
    }else{
      lab <- ""
    }

    p <- EnhancedVolcano(plot.data, lab = lab, x = "logFC", y = y.label, xlab = "LogFC", ylab = paste0("-Log10(", y.label, ")"), pCutoff = 0.05, labSize = 4, FCcutoff = 0, title = paste0(label, ":", titles), subtitle = paste0("Up: ", up.num, "; Down: ", down.num, ", Not sig: ", nrow(plot.data) - up.num - down.num), legendPosition = 'none', colCustom = keyvals) + theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
    print(p)
    return(out)
  })
  names(GroupComparison) <- gsub("-", ".vs.", myContrasts)

  DEG.list <- list(ContinuousComparison = ContinuousComparison)
  DEG.list <- c(DEG.list, GroupComparison)
  return(DEG.list)
}

edgeR_plot <- function(edgeR_result, genes = NULL, PValue = FALSE, threshold = 0.05, top.n = 10, size = 4, title = NULL){
  require(ggpubr)
  require(dplyr)
  require(ggrepel)
  require(ggtext)
  plot.data <- edgeR_result
  if(PValue){
    plot.data$fdr <- plot.data$PValue
    sig.label <- "PValue"
  }else{
    sig.label <- "fdr"
  }

  plot.data$groups <- "NS"
  plot.data$groups[plot.data$logFC > 0 & plot.data$fdr < threshold] <- "UP"
  plot.data$groups[plot.data$logFC < 0 & plot.data$fdr < threshold] <- "DOWN"
  plot.data$groups <- factor(plot.data$groups, levels = c("UP", "DOWN", "NS"))
  up.num <- sum(plot.data$groups == "UP")
  down.num <- sum(plot.data$groups == "DOWN")

  cols <- c("#E64B35", "#4DBBD5", "#D3D3D3")
  names(cols) <- c("UP", "DOWN", "NS")
  # plot.data$point_size <- ifelse(plot.data$groups %in% c("UP", "DOWN"), 1, 1)

  p <- ggscatter(plot.data, x = "logCPM", y = "logFC", color = "groups", xlab = "Average logCPM", ylab = "LogFC", palette = cols, alpha = 0.8, size = 1) + labs(title = title, subtitle = paste0(sig.label, " < ", threshold, " --> Up: ", up.num, "; Down: ", down.num, ", Not sig: ", nrow(plot.data) - up.num - down.num)) + theme(plot.title = element_text(hjust = 0.5, size = 12), aspect.ratio = 1)

  # select to label
  if(top.n > 0){
    top_genes <- plot.data %>% filter(groups %in% c("UP", "DOWN")) %>% group_by(groups) %>% arrange(fdr, .by_group = TRUE) %>% slice_head(n = top.n) %>% pull(gene)
    top_genes <- plot.data %>% filter(gene %in% top_genes)
    p <- p + geom_text_repel(data = top_genes, aes(label = gene, color = groups), size = size, box.padding = 0.5, arrow = arrow(length = unit(0.01, "npc"), type = "closed"), max.overlaps = Inf)
  }

  # label intersect genes
  if(!is.null(genes)){
    if(is.list(genes)){
      # UP
      interest_genes <- plot.data %>% filter(gene %in% genes$UP)
      p <- p + geom_text_repel(data = interest_genes, aes(label = gene), color = "#5F1B74", size = size, box.padding = 0.5, arrow = arrow(length = unit(0.01, "npc"), type = "closed"), max.overlaps = Inf)
      # DOWN
      interest_genes <- plot.data %>% filter(gene %in% genes$DOWN)
      p <- p + geom_text_repel(data = interest_genes, aes(label = gene), color = "#00A087", size = size, box.padding = 0.5, arrow = arrow(length = unit(0.01, "npc"), type = "closed"), max.overlaps = Inf)

      # add legend
      text_data <- data.frame(x = Inf, y = min(plot.data$logFC), label = c("DEGs from other data: ", "UP", "DOWN"), color = c("black", "#5F1B74", "#00A087"), hjust = c(1.5, 3.5, 1))
      p <- p + geom_text(data = text_data, aes(x = x, y = y, label = label, hjust = hjust), color = text_data$color, size = 3.5)
    }else{
      interest_genes <- plot.data %>% filter(gene %in% genes)
      p <- p + geom_text_repel(data = interest_genes, aes(label = gene), color = "#00A087", size = size, box.padding = 0.5, arrow = arrow(length = unit(0.01, "npc"), type = "closed"), max.overlaps = Inf)
      text_data <- data.frame(x = Inf, y = min(plot.data$logFC), label = c("DEGs from other data"), color = c("#00A087"), hjust = 1.5)
      p <- p + geom_text(data = text_data, aes(x = x, y = y, label = label, hjust = hjust), color = text_data$color, size = 3.5)
    }
  }
  return(p)
}

#' @param counts
#' @param group: require order consistent with counts
#' @param myContrasts: m27-m3, d50-d5
#' @param cpm.threshold: filter genes with low counts
edgeR_gene_simple <- function(counts, group, myContrasts, batch = NULL, cpm.threshold = 1, cpm.filter = FALSE, sample.label = NULL){
  require(edgeR)
  cpms <- edgeR::cpm(counts)
  if(!is.factor(group)){
    group <- as.factor(group)
  }
  lvls <- levels(group)
  if(!is.null(batch)){
    if (!is.factor(batch)) {
      batch <- as.factor(batch)
    }
    design <- model.matrix(~0+group+batch)
    colnames(design)[1:length(lvls)] <- lvls
  }else{
    design <- model.matrix(~0+group)
    colnames(design) <- lvls
  }
  if(is.null(batch)){
    cat("Batch: none\n")
  }else{
    cat("Batch: ", paste(levels(batch), collapse = ", "), "\n")
  }
  contrast.matrix <- makeContrasts(contrasts = myContrasts, levels = design)

  cat(nrow(counts), " genes in raw matrix\n")
  obj <- DGEList(counts = counts, group = group)
  if(cpm.filter){
    # cat("Filtering method: cpm >= ", cpm.threshold, " in at least 2 samples in each group\n")
    # # cpm >= 1 in at least 2 samples in each group
    # group_counts <- lapply(levels(group), function(x){
    #   idx <- which(group == x)
    #   count_gt1 <- apply(cpms[, idx], 1, function(x) sum(x >= cpm.threshold, na.rm = TRUE))
    #   return(which(count_gt1 >= 2))
    # })
    # keep <- unique(unlist(group_counts))

    cat("Filtering method: cpm >= ", cpm.threshold, " in at least 2 samples\n")
    count_gt <- apply(cpms, 1, function(x) sum(x >= cpm.threshold, na.rm = TRUE))
    keep <- which(count_gt >= 2)
    genes.left <- length(keep)
  }else{
    cat("Filtering method: edgeR default\n")
    keep <- filterByExpr(obj)
    genes.left <- length(which(keep))
  }
  filter.ratio <- genes.left / nrow(counts)
  cat(genes.left, " genes in filted matrix: ", round(filter.ratio, 2), "\n")

  obj <- obj[keep,, keep.lib.sizes=F]
  obj <-  calcNormFactors(obj)
  my_colors <- RColorBrewer::brewer.pal(n = nlevels(group), name = "Set1")
  if(!is.null(sample.label)){
    plotMDS(obj, col = my_colors[as.numeric(group)], pch = 16, cex = 1)
    mds_data <- plotMDS(obj, col = my_colors[as.numeric(group)], cex = 0.8, plot = FALSE)
    text(x = mds_data$x, y = mds_data$y, labels = sample.label, pos = 3, cex = 0.8, col = my_colors[as.numeric(group)])
    legend("topright", legend = levels(group), col = my_colors, pch = 16)
  }else{
    plotMDS(obj, col = my_colors[as.numeric(group)], pch = 16, cex = 1)
    mds_data <- plotMDS(obj, col = my_colors[as.numeric(group)], cex = 0.8, plot = FALSE)
    text(x = mds_data$x, y = mds_data$y, labels = colnames(counts), pos = 3, cex = 0.8, col = my_colors[as.numeric(group)])
    legend("topright", legend = levels(group), col = my_colors, pch = 16)
  }

  obj <- estimateDisp(obj, design, robust=TRUE)
  bcv <- sqrt(obj$tagwise.dispersion)
  plotBCV(obj)
  fit_tag <- glmFit(obj, design)
  lrt <- glmLRT(fit_tag, contrast = contrast.matrix[,1])
  fdr <- p.adjust(lrt$table$PValue, method = "fdr")
  out <- cbind(cpms[rownames(lrt$table),], lrt$table, fdr)
  out <- as.data.frame(out)
  out$bcv <- bcv
  out$gene <- rownames(lrt$table)
  out <- out[order(out$fdr),]
  return(out)
}
