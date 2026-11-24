#' @description assign colors
#' @author Zhilin Long
#' @date 2025-6-11

library(ggsci)
color.age <- c("#2B2BF9", "#eb1818")
names(color.age) <- c("m3", "m27")

color.age_group <- c("#2B2BF9", "#eb1818")
names(color.age_group) <- c("Young", "Old")

color.olds <- c("#E54B35", "#f9c2b8")
color.youngs <- c("#0A74B2", "#aadefa")

color.sample <- ggsci::pal_igv("default")(51)

replication.color <- ggsci::pal_nejm("default")(2)

color.histone <- c("#2C4B9B", "#512B85", "#EC9022", "#2FB5B3", "#E54B37", "#268550") # "#4DBBD5"
names(color.histone) <- c("H3K9me3", "H3K27me3", "H3K27ac", "H3K4me1", "H3K4me3", "H3K36me3") # "ATAC"

color.histone_age <-ggsci::pal_d3("category20")(14)
names(color.histone_age) <- as.vector(t(outer(names(color.histone), names(color.age), paste, sep = "_")))

color.subCellType <- pal_simpsons("springfield")(16)
names(color.subCellType) <- paste0("C", 1:12)

color.cellType <- c(ggsci::pal_d3("category20")(20), ggsci::pal_d3("category20b")(20), ggsci::pal_d3("category20c")(20))

color.tissue <- c("#1F78B4", "#FC8D62", "#36dbac", "#e04a42", "#94e873", "#724ba7", "#e7e972", "#EC639A", "#247b42", "#8baaf2", "#8C564B",  "#19b9d2", "#a7489e")
names(color.tissue) <- c("brainCB", "brainFC", "brainHip", "heart", "liver", "lung", "kidney", "BAT", "muscle", "spleen", "colon", "skin", "testis")
color.tissue2 <- color.tissue
names(color.tissue2) <- c("Cerebellum", "Frontal cortex", "Hippocampus", "Heart", "Liver", "Lung", "Kidney", "BAT", "Muscle", "Spleen", "Colon", "Skin", "Testis")

color.category <- c("#D46026", "#C684B8", "#46732E", "#E6A023", "#039E73", "#197EC0", "#CE3D32", "#4DBBD5", "#3b5ba1", "#ae2766", "#036473", "#80796B")
names(color.category) <- c("Glutamatergic", "GABAergic", "Glial", "Myeloid", "Lymphoid", "Endothelial", "Epithelial", "Stromal", "Adipose", "Muscle", "Germline", "Others")
