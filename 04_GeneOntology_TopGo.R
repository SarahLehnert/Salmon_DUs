#R script for running Gene Ontology analysis using TopGo
#Script similar to used in previous publication
#Adapted from scripts provided by Wellband et al. 2019 (Mol Ecol)

#set working directory
setwd("~/Desktop/Sarah/Salmon/GeneOntology/DU2_220K/")

#Libraries
library(topGO)

#all gene ontology for salmo
#this uses the original salmon assembly
go.anno <- read.csv("~/Desktop/Sarah/Salmon/GeneOntology/Ssal_ICSASG_v2_GOAccession.csv", header = T, stringsAsFactors = F, sep="\t")
go.anno[,1] <- gsub("\\.t[0-9]*", "", go.anno[,1])
go.anno=go.anno[,c(1,4)]


#Outliers (intersected) - DU2 outliers
outliers <- read.table("outliers_DU2_220Kdata.accnos", stringsAsFactors = F)[,1]
outliers <- gsub("\\.t[0-9]*", "", outliers)


#Get outliers in annotated gene ontology
outliers.w.go <- outliers[outliers %in% go.anno[,1]]

#all genes -all loci used in PCA 
all.WGS.genes <- read.table("all_loci_DU2_220Kdata.accnos", stringsAsFactors = F)[,1]
all.WGS.genes.w.go <- all.WGS.genes[all.WGS.genes %in% go.anno[,1]]


#Create function for topGo
CreateGene2GOMapping <- function(x) {
  nr <- nrow(x)
  map <- list()
  for(r in 1:nr) {
    map[[as.character(x[r, 1])]] <- append(map[[as.character(x[r, 1])]], x[r,2])
  }
  return(map)
}


# create GO mapping, takes about 20 minutes
gene2GO.map <- CreateGene2GOMapping(go.anno)

gene2GO.map <- gene2GO.map[names(gene2GO.map) %in% all.WGS.genes.w.go]
gene.list <- factor(as.integer(all.WGS.genes.w.go %in% outliers.w.go))
names(gene.list) <- all.WGS.genes.w.go

outlier_GOdata <- new("topGOdata",
                      description = "Outliers DU2 Labrador", ontology = "BP",
                      allGenes = gene.list,
                      nodeSize = 5,
                      annot = annFUN.gene2GO, gene2GO = gene2GO.map)

outlier.result.Fisher <- runTest(outlier_GOdata, algorithm = "weight01", statistic = "fisher")

geneData(outlier.result.Fisher)
hist(score(outlier.result.Fisher), 50, xlab = "p-values")

getwd()
write.table(GenTable(outlier_GOdata, outlier.result.Fisher, topNodes = length(outlier.result.Fisher@score)), 
            "GeneOntology_DU2_Labrador_220K.txt", quote = F, row.names = F , sep="\t", col.names = T)


