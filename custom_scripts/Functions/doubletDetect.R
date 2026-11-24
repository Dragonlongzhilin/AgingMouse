########## ---- scRNA-seq ---- ##########

####Using DoubletFinder function for Doublet detection
# The input Seurat object, the request has been filtered, PCA, FindClusters, etc.

#' @param annotation seurat result
#' @param GT GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
#' @param doublet.rate the doublet rate, doubletDetect.doubletRate.txt or automatically calculated based on 10X

#' @author Zhilin Long

FunctionDoubletFinder <- function(Seurat.object, PCs, doublet.rate = "10X", annotation = NULL, pN_value = 0.25, GT = FALSE, sct = FALSE, seed = 12345, num.cores = 10){
    set.seed(seed)
    require(DoubletFinder) #Require cleanup of low-quality cells in advance; v2.0.4
    patient <- paramSweep(Seurat.object, PCs = PCs, sct = sct, num.cores = num.cores)
    gt.calls <- summarizeSweep(patient, GT = GT)
    
    if(!is.numeric(doublet.rate)){
        if(doublet.rate == "10X"){
            dbr <- 0.008 * (nrow(Seurat.object@meta.data) / 1000) # 10X genomics
        }else{
            doubletRate <- read.table(doublet.rate, header = T, sep = "\t", stringsAsFactors = F)
            cell.ranges <- round(nrow(Seurat.object@meta.data)/1000, 0)*1000
            dbr <- doubletRate[which(doubletRate[, 2] == cell.ranges), 3]/100
        }
    }else{
        dbr <- doublet.rate
    }
    cat("    doublet.rate: ", dbr, "\n")

    #Computes and visualizes the mean-variance normalized bimodality coefficient (BCmvn) score for each pK value tested during doubletFinder_ParamSweep. 
    #Optimal pK for any scRNA-seq data can be manually discerned as maxima in BCmvn distributions. 
    #If ground-truth doublet classifications are available, BCmvn is plotted along with mean ROC AUC for each pK.
    sweep <- find.pK(gt.calls)
    pK_value <- as.numeric(as.character(sweep$pK[sweep$BCmetric == max(sweep$BCmetric)])) # optimal pK

    par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
    plot(x = as.numeric(as.character(sweep$pK)), y = sweep$BCmetric, pch = 16,type="b", col = "blue",lty=1)
    abline(v=pK_value,lwd=2,col='red',lty=2)
    title("The BCmvn distributions")
    text(pK_value,max(sweep$BCmetric),as.character(pK_value),pos = 4,col = "red")

    #potential doublet rate，
    nExp_poi <- round(dbr*nrow(Seurat.object@meta.data))
    if(is.null(annotation)){
        Seurat.object <- doubletFinder(Seurat.object, PCs = PCs, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = NULL, sct = sct)
        label <- paste0("DF.classifications_", pN_value, "_", pK_value,'_', nExp_poi)
    }else{
        annotations <- Seurat.object@meta.data[, annotation]
        homotypic.prop <- modelHomotypic(annotations)
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        pANN_value <- paste0("pANN_", pN_value, "_", pK_value, '_', nExp_poi)

        # https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/225#issuecomment-2786505997
        Seurat.object <- doubletFinder(Seurat.object, PCs = PCs, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = NULL, sct = sct)
        Seurat.object <- doubletFinder(Seurat.object, PCs = PCs, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value, sct = sct)
        label <- paste0("DF.classifications_", pN_value, "_", pK_value,'_', nExp_poi.adj)
    }
    Seurat.object@meta.data$DoubletFinder.prediction <- Seurat.object@meta.data[, label]
    
    # remove label information
    idx <- grep("pANN|DF.classifications", colnames(Seurat.object@meta.data))
    Seurat.object@meta.data <- Seurat.object@meta.data[,-idx]

    return(Seurat.object)
}

####Using scDblFinder function for Doublet detection
# https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html

#' @param Seurat.object
#' @param clusters, cluster-based approach
#' @param samples, preferable to look for doublets separately for each sample
#' @param n.core, threads
#' @param dbr The expected doublet rate. By default this is assumed to be 1% per thousand cells captured (so 4% among 4000 thousand cells), which is appropriate for 10x datasets. 

#' @note Of note, if you have very large differences in number of cells between samples the scores will not be directly comparable.
# We are working on improving this, but in the meantime it would be preferable to stratify similar samples and threshold the sets separately.
FunctionScDblFinder.RNA <- function(Seurat.object, clusters = NULL, samples = NULL, doublet.rate = NULL, nfeatures = 1352, dims = 20, n.core = 10, seed = 12345){
    set.seed(seed)
    require(scDblFinder)
    require(BiocParallel)

    if(!is.numeric(doublet.rate)){
        if(is.null(doublet.rate)){
            dbr <- NULL # default value, Unless manually given, the expected number of doublets (e) is specified by e = n2/10−5 (where n is the number of cells captured).
            cat("    doublet.rate estimated by scDblFinder method \n")
        }else{
            dbr <- 0.008 * (nrow(Seurat.object@meta.data) / 1000) # 10X genomics
            cat("    doublet.rate estimated by 10X method \n")
        }
    }else{
        dbr <- doublet.rate
        cat("    doublet.rate: ", dbr, "\n")
    }

    sce <- as.SingleCellExperiment(Seurat.object)
    sce <- scDblFinder(sce, clusters = clusters, samples = samples, nfeatures = nfeatures, dims = dims, dbr = dbr, BPPARAM = MulticoreParam(n.core))
    Seurat.object <- AddMetaData(Seurat.object, sce$scDblFinder.class, "scDblFinder.prediction")

    # keep the label consistent with DoubleFinder method
    Seurat.object$scDblFinder.prediction <- gsub("doublet", "Doublet", Seurat.object$scDblFinder.prediction) 
    Seurat.object$scDblFinder.prediction <- gsub("singlet", "Singlet", Seurat.object$scDblFinder.prediction)
    Seurat.object <- AddMetaData(Seurat.object, sce$scDblFinder.score, "scDblFinder.score")
    # Seurat.object <- AddMetaData(Seurat.object, sce$scDblFinder.ratio, "scDblFinder.ratio")
    # Seurat.object <- AddMetaData(Seurat.object, sce$scDblFinder.cluster, "scDblFinder.cluster")
    return(Seurat.object)
}
# sce$scDblFinder.ratio : the proportion of artificial doublets among the neighborhood (the higher, the more chances that the cell is a doublet)
# sce$scDblFinder.weighted : the proportion of artificial doublets among the neighborhood, weighted by distance
# sce$scDblFinder.score : the final doublet score
# sce$scDblFinder.class : the classification (doublet or singlet)

####Using Scrublet function for Doublet detection
# https://rdrr.io/github/compbiomed/singleCellTK/man/runScrublet.html



########## ---- scATAC-seq ---- ##########

####Using scDblFinder function for Doublet detection
# Since working with all features (i.e. tiles or peaks) is computationally very expensive, an alternative approach is to begin by reducing the size of the dataset, not through selection (as in scRNAseq), but by aggregating correlated features into a relatively small set. 
# This has the advantage of using all information, as well as making the count data more continuous.

FunctionScDblFinder.ATAC <- function(Seurat.object, clusters = NULL, samples = NULL, nfeatures = 25, dims = 20, dbr = NULL, n.core = 10, seed = 12345){
    set.seed(seed)
    require(scDblFinder)
    require(BiocParallel)
    sce <- as.SingleCellExperiment(Seurat.object)
    sce <- scDblFinder(sce, clusters = clusters, samples = samples, aggregateFeatures=TRUE, nfeatures=nfeatures, processing="normFeatures", dims = dims, BPPARAM = MulticoreParam(n.core))
    Seurat.object <- AddMetaData(Seurat.object, sce$scDblFinder.class, "scDblFinder.prediction")
    Seurat.object <- AddMetaData(Seurat.object, sce$scDblFinder.score, "scDblFinder.score")
    # Seurat.object <- AddMetaData(Seurat.object, sce$scDblFinder.ratio, "scDblFinder.ratio")
    # Seurat.object <- AddMetaData(Seurat.object, sce$scDblFinder.cluster, "scDblFinder.cluster")
    return(Seurat.object)
}

#### AMULET method: AMULET: a novel read count-based method for effective multiplet detection from single nucleus ATAC-seq data. 2021 genome biology
# cd /data/ExtraDisk/sdd/longzhilin/softwares/AMULET-1.1
# ./AMULET.sh \
# /data/ExtraDisk/sdd/longzhilin/Project/Thalamus.20211022/Data/cellRangerResults/scATAC/aggr/fragments.tsv.gz \
# /data/ExtraDisk/sdd/longzhilin/Project/Thalamus.20211022/Data/cellRangerResults/scATAC/aggr/singlecell.csv \
# /data/ExtraDisk/sdd/longzhilin/softwares/AMULET-1.1/human_autosomes.txt \
# /data/ExtraDisk/sdd/longzhilin/softwares/AMULET-1.1/knownRepetitiveElements/RepeatFilterFiles/blacklist_repeats_segdups_rmsk_hg38.bed \
# /data/ExtraDisk/sdd/longzhilin/Project/Thalamus.20211022/Result/1.QualityControl/scATAC/AMULET/ \
# /data/ExtraDisk/sdd/longzhilin/softwares/AMULET-1.1/

FunctionAmulet.ATAC <- function(Seurat.object, AMULET.resultFile, threshold = 0.01){
    AMULET.result <- read.table(file = AMULET.resultFile, header = T, stringsAsFactors = F, sep = "\t")
    cells <- Cells(Seurat.object)
    cells <- gsub(".*_", "", cells)
    AMULET.Pvalue <- AMULET.result$p.value
    names(AMULET.Pvalue) <- AMULET.result$cell_id
    AMULET.Qvalue <- AMULET.result$q.value
    names(AMULET.Qvalue) <- AMULET.result$cell_id
    Seurat.object <- AddMetaData(Seurat.object, AMULET.Pvalue, "AMULET.Pvalue")
    Seurat.object <- AddMetaData(Seurat.object, AMULET.Qvalue, "AMULET.Qvalue")    
    Seurat.object$AMULET.class <- rep("doublet", length(Seurat.object$AMULET.Qvalue))
    Seurat.object$AMULET.class[which(Seurat.object$AMULET.Qvalue > threshold)] <- "singlet" # more larger --- more doublet
    return(Seurat.object)
}

# The Amulet method tends to perform best with datasets that have homotypic doublets and where cells have a high library size (i.e. median library size per cell of 10-15k reads), 
# while the scDblFinder-based approach works better for heterotypic doublets. 
# Until an optimal solution is found, we recommend using multiple approaches to inform decisions, 
# in particular using the p-value combination method below.
# ref: https://bioconductor.org/packages/devel/bioc/vignettes/scDblFinder/inst/doc/scATAC.html

FunctionCombined.AMULET_scDblFinder <- function(Seurat.object, AMULET.resultFile, threshold = 0.05){
    AMULET.result <- read.table(file = AMULET.resultFile, header = T, stringsAsFactors = F, sep = "\t")
    cells <- Cells(Seurat.object)
    cells <- gsub(".*_", "", cells)
    idx <- match(cells, AMULET.result$cell_id)
    Seurat.object$AMULET.Pvalue <- AMULET.result$p.value[idx]
    Seurat.object$AMULET.Qvalue <- AMULET.result$q.value[idx]
    Seurat.object$AMULET.class <- rep("doublet", length(Seurat.object$AMULET.Qvalue))
    Seurat.object$AMULET.class[which(Seurat.object$AMULET.Qvalue > 0.01)] <- "singlet" # more larger --- more doublet
    Seurat.object$scDblFinder.Pvalue <- 1 - Seurat.object$scDblFinder.score

    scDblFinder.min.value <- sort(unique(Seurat.object$scDblFinder.Pvalue))
    if(scDblFinder.min.value[1] == 0){
        scDblFinder.min.value <- scDblFinder.min.value[2]
    }else{
        scDblFinder.min.value <- scDblFinder.min.value[1]
    }
    AMULET.min.value <- sort(unique(Seurat.object$AMULET.Pvalue))
    if(AMULET.min.value[1] == 0){
        AMULET.min.value <- AMULET.min.value[2]
    }else{
        AMULET.min.value <- AMULET.min.value[1]
    }

    combined <- apply(Seurat.object@meta.data[, c("scDblFinder.Pvalue", "AMULET.Pvalue")], 1, function(x){
        if(x[1] == 0){
            x[1] <- scDblFinder.min.value
        }
        if(x[2] == 0){
            x[2] <- AMULET.min.value
        }
        suppressWarnings(aggregation::fisher(x))
    })
    Seurat.object$combined.Pvalue <- combined
    Seurat.object$Multiplet <- rep("singlet", length(Seurat.object$combined.Pvalue))
    Seurat.object$Multiplet[which(Seurat.object$combined.Pvalue < threshold)] <- "doublet"
    return(Seurat.object)
}
