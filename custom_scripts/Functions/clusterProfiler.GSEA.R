#### 使用clusterProfiler进行通路富集分析
#' @author Zhilin Long
#' @date 2025-5-16

#' @param geneRank 要求输入一个已经拍好序的rank数值向量, name属性为对应的gene id。为了保证后续正确，函数内部我增加了排序的功能，所以可以支持non-rank的数据
#' 如
#' A  B C
#' 5.2  1.2 0.8
#' @param geneType 字符型向量，指明输入的gene中name属性类型，即基因id的类型
#' @param db.type 字符型向量，指明GSEA针对的数据库
#' @param saveDir gsea输入目录
#' @param pAdjustMethod, pvalueCutoff p值矫正方法和阈值
#' @param MsigDB.category 指定特定的MsigDB中的数据库进行GSEA分析
#' MsigDB.pathway, msigdbr(species = "Homo sapiens", category = i)， 指定i为H, C1, C2, C3, C4, C5, C6, C7, C8
# 通过gs_cat gs_subcat，来选择库
#' @param GO.ont 指定GSEA分析GO数据库时使用的条目类型，支持"BP", "MF", and "CC"，以及"ALL"
#' @param simplify 逻辑值，指明是否需要采用simplify方法过滤冗余的go term结果：https://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
# GO分析默认为BP
#' minGSSsize默认是10-500
clusterProfiler.GSEA <- function(geneRank, geneType = c("ENTREZID", "SYMBOL", "ENSEMBL"),
                                 db.type = c("MsigDB", "GO"),
                                 species = c("Homo sapiens", "Mus musculus", "Drosophila melanogaster"),
                                 # saveDir,
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05,
                                 MsigDB.category = list(H = c("All"), C2 = c("KEGG", "REACTOME"), C5 = c("BP")),
                                 GO.ont = "BP", simplify = F, seed = 101,
                                 minGSSize = 10, maxGSSize = 500, showCategory = 20, return.gsea = F){
  require(enrichplot)
  require(clusterProfiler)
  require(tidyverse)
  #pdf(file = file.path(saveDir, "GSEAplot.pdf"))
  geneType <- match.arg(geneType)
  db.type <- match.arg(db.type)  
  species <- match.arg(species)

  cat("Input the gene set belong to", species, "\n")

  # GSEA requires the input gene, which is after rank
  geneRank <- sort(geneRank, decreasing = TRUE)
 
  #### GSEA analysis of MsigDB database
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
    # run GSEA
    gsea.res <- GSEA(geneRank, TERM2GENE = MsigDB.pathway, pvalueCutoff = pvalueCutoff, pAdjustMethod = pAdjustMethod, verbose = FALSE, maxGSSize = maxGSSize, minGSSize = minGSSize, seed = seed)
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

  #### GSEA analysis of Gene Ontology
  if(db.type == "GO"){
    gsea.res <- gseGO(gene = geneRank, OrgDb = db, ont = GO.ont, keyType = geneType,
                      pvalueCutoff = pvalueCutoff, maxGSSize = maxGSSize, minGSSize = minGSSize,
                      verbose = FALSE, seed = seed)
    # simplify method, we can remove redundant terms.The function internally called GOSemSim (Yu et al. 2010) to calculate semantic similarity among GO terms and remove those highly similar terms by keeping one representative term. 
    if(simplify){
      cat("simplify GO enrichment!\n")
      gsea.res <- simplify(gsea.res, cutoff = 0.7, by="p.adjust", select_fun = min)
    }
  }

  # #### GSEA analysis of KEGG database
  # if(db.type == "KEGG"){
  #   gsea.res <- gseKEGG(gene = geneRank,
  #                       organism = KEGG.organism, keyType = "kegg",
  #                       pvalueCutoff = pvalueCutoff, minGSSize = kegg.minGSSize,
  #                       verbose = FALSE, seed = seed)
  # }

  # #转换为geneSymbol
  # gsea.res.geneSymbol <- setReadable(gsea.res, org.Hs.eg.db, keyType = "ENTREZID")
  # if(nrow(gsea.res@result) > 0){
  #   write.csv(gsea.res.geneSymbol, file = file.path(saveDir, "GSEAResult.csv"))

  #   p1 <- dotplot(gsea.res.geneSymbol, showCategory = showCategory)
  #   print(p1)

  #   p2.list <- sapply(1:nrow(gsea.res.geneSymbol@result), function(x){
  #     p2 <- gseaplot2(gsea.res.geneSymbol, geneSetID = x, title = gsea.res.geneSymbol$Description[x])
  #     print(p2)
  #   })
  # }
  if(return.gsea){
    return(gsea.res)
  }else{
    return(gsea.res@result)
  }
}
