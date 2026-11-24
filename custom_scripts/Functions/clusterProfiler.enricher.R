#' @description pathway enrichment by clusterProfiler
#' @author Zhilin Long
#' @param gene: character，input genes. Such as: FUT7 AIRE IFI35
#' @param geneType: character, type of input gene. ENTREZID or SYMBOL or ENSEMBL
#' @param db.type: select the database type
#' @param saveDir: save plot files
#' @param pAdjustMethod: the method to adjust p values
#' @param MsigDB.category: set the category of MsigDB for enrichment
#' MsigDB.pathway, msigdbr(species = "Homo sapiens", category = i)， set i is H, C1, C2, C3, C4, C5, C6, C7, C8

#' @param GO.ont specifies the entry type used when GSEA analyzes the GO database, supporting "BP", "MF", and "CC", as well as "ALL"

# BgRatio = unique(pathway genes) / unique(union genes from all selected pathways)
# GeneRatio = intersection genes between DEGs and pathway genes / intersection genes between DEGs and union genes from all selected pathways

#' The principle of functional enrichment analysis is based on hypergeometric test
#' clusterProfiler mainly performs analysis based on EntrezID
#' The universe in the analysis process, that is, the number of background genes, is related to the subcategory in the db you choose, and is not fixed
clusterProfiler.enricher <- function(gene, geneType = c("ENTREZID", "SYMBOL", "ENSEMBL"),
                                     db.type = c("MsigDB", "GO"),
                                     # saveDir = NULL, title = NULL,
                                     pAdjustMethod = "BH",
                                     qvalueCutoff = 0.2,
                                     pvalueCutoff = 0.05,
                                     species = c("Homo sapiens", "Mus musculus", "Drosophila melanogaster"),
                                     MsigDB.category = list(H = c("All"), C2 = c("KEGG", "REACTOME"), C5 = c("BP")),
                                     GO.ont = "BP", simplify = F, return.clusterProfiler = F,
                                     minGSSize = 10, maxGSSize = 500, showCategory = 20
                                     ){
  require(enrichplot)
  require(clusterProfiler)
  require(tidyverse)

  geneType <- match.arg(geneType)
  db.type <- match.arg(db.type)
  species <- match.arg(species)

  cat("Input the gene set belong to", species, "\n")
  if(length(gene) < minGSSize){
    cat("Input the gene set is too small to do the pathway enrichment!\n")
    return(NULL)
  }else{
    #### MSigDb over-presentaton analysis
    if(db.type == "MsigDB"){
      library(msigdbr)
      gene_type <- switch(geneType, "ENTREZID" = "entrez_gene", "SYMBOL" = "gene_symbol", "ENSEMBL" = "ensembl_gene", stop("Invalid geneType! Must be 'ENTREZID', 'SYMBOL', or 'ENSEMBL'"))
      # extract the pathway gene sets
      MsigDB.pathway <- tibble(gs_name = NULL, gene_type = NULL)
      if(class(MsigDB.category) == "character"){
        # for example: H, C1, C2
        MsigDB.pathway <- msigdbr(species = species, category = MsigDB.category)
      }else{
        for(i in names(MsigDB.category)){
          if(MsigDB.category[[i]][1] == "All"){
            MsigDB.pathway <- rbind(MsigDB.pathway, msigdbr(species = species, category = i))
          }else{
            for(j in MsigDB.category[[i]]){
              MsigDB.pathway <- rbind(MsigDB.pathway, msigdbr(species = species, category = i, subcategory = j))
            }
          }
        }
      }
      MsigDB.pathway <- MsigDB.pathway[,c("gs_name", gene_type)]
      em.res <- enricher(gene, TERM2GENE = MsigDB.pathway, qvalueCutoff = qvalueCutoff, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize)
    }
    
    if(species == "Mus musculus"){
      require(org.Mm.eg.db)
      db <- org.Mm.eg.db
    }else if(species == "Drosophila melanogaster"){
      require(org.Dm.eg.db)
      db <- org.Dm.eg.db
    }else{
      require(org.Hs.eg.db)
      db <- org.Hs.eg.db
    }

    #### GO over-representation analysis
    if(db.type == "GO"){
      em.res <- enrichGO(gene = gene,
                         OrgDb = db, keyType = geneType,
                         ont = GO.ont, minGSSize = minGSSize, maxGSSize = maxGSSize,
                         pAdjustMethod = pAdjustMethod,
                         pvalueCutoff  = pvalueCutoff,
                         qvalueCutoff  = qvalueCutoff)
      # simplify method, we can remove redundant terms.The function internally called GOSemSim (Yu et al. 2010) to calculate semantic similarity among GO terms and remove those highly similar terms by keeping one representative term. 
      if(simplify){
        cat("simplify GO enrichment!\n")
        em.res <- simplify(em.res, cutoff = 0.7, by="p.adjust", select_fun = min)
      }
    }
    
    if(return.clusterProfiler){
      return(em.res)
    }else{
      # convert geneSymbol
      if(geneType == "ENTREZID"){
        em.res <- setReadable(em.res, db, keyType = "ENTREZID")
      }
      return(em.res@result)
    }

    # abandoned, due to convert ID error
    # #### KEGG over-representation analysis
    # # https://www.genome.jp/kegg/catalog/org_list.html
    # if(db.type == "KEGG"){
    #   # If the gene name type is not ENTREZID, ID conversion is required
    #   if(geneType == "ENTREZID"){
    #     cat("Do not convert gene id\n")
    #     geneList <- gene
    #   }
    #   if(geneType %in% c("SYMBOL", "ENSEMBL")){
    #     converted_genes <- bitr(geneID = gene, fromType = geneType, toType = "ENTREZID", OrgDb = OrgDb)
    #     cat("Before convert: ", length(gene), "\n")
    #     index <- match(gene, converted_genes[,1])
    #     geneList <- as.character(converted_genes[na.omit(index),2])
    #     cat("After convert: ", length(geneList), "\n")    
    #   }
    #   em.res <- enrichKEGG(gene = geneList, organism = KEGG.organism, keyType = "kegg",
    #                        pAdjustMethod = pAdjustMethod, minGSSize = minGSSize, maxGSSize = maxGSSize,
    #                        pvalueCutoff  = pvalueCutoff,
    #                        qvalueCutoff  = qvalueCutoff)
    # }

    # # convert geneSymbol
    # em.res.geneSymbol <- setReadable(em.res, db, keyType = "ENTREZID")
    
    # if(nrow(em.res@result) > 0){
    #   write.csv(em.res.geneSymbol, file = file.path(saveDir, paste0(title, "_", "enricherResult.csv")))
    #   tryCatch({
    #     p1 <- barplot(em.res.geneSymbol, showCategory = showCategory, main = title) + theme(axis.text.x = element_text(size=4), axis.text.y = element_text(size=4))
    #     print(p1)
    #     p2 <- dotplot(em.res.geneSymbol, showCategory = showCategory, title = title) + theme(axis.text.x = element_text(size=4), axis.text.y = element_text(size=4))
    #     print(p2)
    #     p3 <- heatplot(em.res.geneSymbol, showCategory = showCategory) + theme(axis.text.x = element_text(size=4), axis.text.y = element_text(size=4))
    #     print(p3)
    #   }, warning = function(w){
    #     cat("Warning in plot!\n")
    #   }, error = function(e){
    #     cat("Error in plot!\n")
    #   },finally = {
    #     cat("Plot complete!\n")
    #   })
    # }
  }
}