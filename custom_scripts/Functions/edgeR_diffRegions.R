#' @description differential regions analysis among different age with edgeR method
#' @Author: Zhilin Long
#' @date 2025-7-18

#' @param group factor, e.g., c("d5", "d30", "d50", "d70")
source("/storage/zhangyanxiaoLab/longzhilin/ProjectCode/AgingProject/Bin2GRange.R")
edgeR.diffRegions <- function(countTables, group, outPath, type, batch = NULL, logFC = 1, PValue = 0.05, title = NULL, ContinuousComparison = TRUE, min.regions = 1000, cpm.filter = FALSE, cpm.threshold = 1, sample.label = NULL){
  require(edgeR)
  require(EnhancedVolcano)
  if(!is.factor(group)){
    group <- as.factor(group)
  }
  lvls <- levels(group)
  pdf(paste0(outPath, "/", type, ".pdf"))
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
  cat("Batch: ", batch, "\n")
  cpms <- edgeR::cpm(countTables)

  cat(nrow(countTables), " regions in raw matrix\n")
  y <- DGEList(counts = countTables, group = group)
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
    regions.left <- length(keep)
  }else{
    cat("Filtering method: edgeR default\n")
    keep <- filterByExpr(y)
    regions.left <- length(which(keep))
  }
  filter.ratio <- regions.left / nrow(countTables)
  cat(regions.left, " regions in filted matrix: ", round(filter.ratio, 2), "\n")

  if(regions.left > min.regions){
    y <- y[keep,, keep.lib.sizes=F]
    y <-  calcNormFactors(y)
    n.groups <- nlevels(group)
    my_colors <- if(n.groups <= 9) RColorBrewer::brewer.pal(n = n.groups, name = "Set1") else colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(n.groups)
    if(!is.null(sample.label)){
      plotMDS(y, col = my_colors[as.numeric(group)], pch = 16, cex = 1, labels = rep("", ncol(countTables)))
      mds_data <- plotMDS(y, col = my_colors[as.numeric(group)], cex = 0.8, plot = FALSE)
      text(x = mds_data$x, y = mds_data$y, labels = sample.label, pos = 3, cex = 0.8, col = my_colors[as.numeric(group)])
      legend("topright", legend = levels(group), col = my_colors, pch = 16)
    }else{
      plotMDS(y, col = my_colors[as.numeric(group)], pch = 16, cex = 1, labels = rep("", ncol(countTables)))
      mds_data <- plotMDS(y, col = my_colors[as.numeric(group)], cex = 0.8, plot = FALSE)
      text(x = mds_data$x, y = mds_data$y, labels = colnames(countTables), pos = 3, cex = 0.8, col = my_colors[as.numeric(group)])
      legend("topright", legend = levels(group), col = my_colors, pch = 16)
    }

    y <- tryCatch({
      res <- estimateDisp(y, design, robust = TRUE)
      }, error = function(e) {
      message("estimateDisp error: ", e$message)
      return(NULL)
    })
    # y <- estimateDisp(y, design, robust = TRUE)

    if(is.null(y)){
      diffRegions <- NULL
    }else{
      bcv <- sqrt(y$tagwise.dispersion)
      plotBCV(y)
      fit_tag <- glmFit(y, design)
      
      # Continuous comparison
      if(ContinuousComparison){
        cat("edgeR: continuous comparison \n")
        len <- length(lvls)
        myContrasts <- paste(lvls[2:len],lvls[1:(len-1)],sep='-')
        contrast.matrix <- makeContrasts(contrasts = myContrasts, levels = design)
        lrt <- glmLRT(fit_tag, contrast = contrast.matrix)
        fdr <- p.adjust(lrt$table$PValue, method="fdr")
        ContinuousComparison.res <- cbind(cpms[rownames(lrt$table),],lrt$table,fdr)
        ContinuousComparison.res <- as.data.frame(ContinuousComparison.res)
        ContinuousComparison.res$Region <- rownames(lrt$table)
        ContinuousComparison.res$bcv <- bcv
        ContinuousComparison.res <- ContinuousComparison.res[order(ContinuousComparison.res$PValue),]

        # regions <- as.data.frame(Bin2GRange(ContinuousComparison$Region, "_"))
        # logFC.idx <- grep("logFC.", colnames(ContinuousComparison))
        # idx <- which((abs(ContinuousComparison$logFC.d30.d05) > logFC | abs(ContinuousComparison$logFC.d50.d30) > logFC) & ContinuousComparison$PValue < PValue)
        # if(length(idx) > 0){
        #   ContinuousComparison.bed <- data.frame(Chr = Regions[idx,1], Start = Regions[idx,2], End = Regions[idx,3], Name = ContinuousComparison$Region[idx], PValue = ContinuousComparison$PValue[idx], Strand = ".", logFC.d30.d05 = ContinuousComparison$logFC.d30.d05[idx], logFC.d50.d30 = ContinuousComparison$logFC.d50.d30[idx])
        #   write.table(ContinuousComparison.bed, paste0(outPath, "/", type, "_ContinuousComparison.bed"), sep = "\t", row.names = F, col.names = F, quote = F)
        # }
      }else{
        ContinuousComparison.res <- NULL
      }

      # Group comparison
      combinations <- combn(levels(group), 2)
      myContrasts <- paste0(combinations[2,], "-", combinations[1,])
      contrast.matrix <- makeContrasts(contrasts = myContrasts, levels = design)
      GroupComparison <- lapply(1:length(myContrasts), function(i){
        lrt <- glmLRT(fit_tag, contrast = contrast.matrix[,i])
        plotMD(lrt)
        fdr <- p.adjust(lrt$table$PValue, method="fdr")
        out <- cbind(cpms[rownames(lrt$table),], lrt$table, fdr)
        out <- as.data.frame(out)
        out$Region <- rownames(lrt$table)
        out$bcv <- bcv
        out <- out[order(out$PValue),]

        if(is.null(title)){
          title <- paste0(type, ":", gsub("-", ".vs.", myContrasts[i]))
        }

        # Regions <- as.data.frame(Bin2GRange(out$Region, "_"))
        # up.idx <- which(out$logFC > logFC & out$PValue < PValue)
        # down.idx <- which(out$logFC < (logFC*-1) & out$PValue < PValue)
        # if(length(up.idx) > 0){
        #   GroupComparison.up <- data.frame(Chr = Regions[up.idx,1], Start = Regions[up.idx,2], End = Regions[up.idx,3], logFC = out$logFC[up.idx], PValue = out$PValue[up.idx], Strand = ".")
        #   write.table(GroupComparison.up, paste0(outPath, "/", type, "_", title, "_up.bed"), sep = "\t", row.names = F, col.names = F, quote = F)
        # }
        # if(length(down.idx) > 0){
        #   GroupComparison.down <- data.frame(Chr = Regions[down.idx,1], Start = Regions[down.idx,2], End = Regions[down.idx,3], logFC = out$logFC[down.idx], PValue = out$PValue[down.idx], Strand = ".")
        #   write.table(GroupComparison.down, paste0(outPath, "/", type, "_", title, "_down.bed"), sep = "\t", row.names = F, col.names = F, quote = F)
        # }

        plot.data <- out
        plot.data$diff <- "NS"
        plot.data$diff[plot.data$logFC > logFC & plot.data$PValue < PValue] <- "UP"
        plot.data$diff[plot.data$logFC < (logFC*-1) & plot.data$PValue < PValue] <- "DOWN"
        up.num <- sum(plot.data$diff == "UP")
        down.num <- sum(plot.data$diff == "DOWN")
        keyvals <- ifelse(
          plot.data$diff == 'UP', 'red3',
          ifelse(plot.data$diff == 'DOWN', 'royalblue', 'grey30')
        )
        names(keyvals) <- plot.data$diff
        p <- EnhancedVolcano(plot.data, lab = "", x = "logFC", y = "PValue", xlab = "LogFC", ylab = "-Log10(P-value)", pCutoff = PValue, FCcutoff = logFC, title = title, subtitle = paste0("Up: ", up.num, "; Down: ", down.num, ", Not sig: ", nrow(plot.data) - up.num - down.num), legendPosition = 'none', colCustom = keyvals) + theme(plot.title = element_text(hjust = 0.5), aspect.ratio = 1)
        print(p)
        return(out)
      })
      names(GroupComparison) <- paste0("GroupComparison.", gsub("-", ".vs.", myContrasts))

      diffRegions <- list(ContinuousComparison = ContinuousComparison.res)
      diffRegions <- c(diffRegions, GroupComparison)
    }
  }else{
    diffRegions <- NULL
  }
  dev.off()
  return(diffRegions)
}
