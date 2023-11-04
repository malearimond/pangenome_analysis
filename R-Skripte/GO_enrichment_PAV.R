install.packages("topGO")
if (!requireNamespace("BiocManager", quietly = TRUE))install.packages("BiocManager")
BiocManager::install("topGO")
BiocManager::install("Rgraphviz")

library(topGO)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("topGO")
#BiocManager::install("Rgraphviz")
library(Rgraphviz)
# BiocManager::install("AnnotationDbi", force = TRUE)
BPterms <- ls(GOBPTerm)

################################
setwd("~/Studium/Köln_Biological Science/Master Thesis/Results/")
####################################GO analysis with Version 5.1

###bring the GOassignment into the right format
#Bastiaan
#GOassignment<-read.table(file = "~/Studium/Köln_Biological Science/Master Thesis/Results/eggnog_geneID2GO.txt", sep="&", stringsAsFactors = F)
#new approach
#gene2GO_new<- readMappings(file = "~/Studium/Köln_Biological Science/Master Thesis/Results/eggnog_geneID2GO_wo_comma.txt", sep = "\t",IDsep = " ")

####################################################################################
#Read in eggnog annotation file
#(remove commas in original file)
gene2GO_new<- readMappings(file = "~/Studium/Köln_Biological Science/Master Thesis/Results/eggnog_geneID2GO_wo_comma.txt", sep = "\t",IDsep = " ")
geneNames <- names(gene2GO_new)

#### core_genes ####
#read in the gene file you want to analyse
genes<-read.table("~/Studium/Köln_Biological Science/Master Thesis/Results/core_genes_wo_contig_2.txt", sep="\t")

genes<-genes[order(genes$V1, decreasing = F),]

#look at overlap between your genes and overall genes
#OTHER WAY AROUND
geneList <- factor(as.integer(geneNames %in% genes))
names(geneList) <- geneNames
str(geneList)

#run GO in BP (Biological processes): BUILD TopTO OBJECT
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = gene2GO_new)

sg <- sigGenes(GOdata)
str(sg)
numSigGenes(GOdata)


####enrichment test

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher

resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
resultKS.elim


par(cex=0.2)


showSigOfNodes(GOdata, score(resultKS.elim),firstSigNodes = 4, useInfo ='def',.NO.CHAR = 15, plotFunction = GOplot)
title(main="GO enrichment core genes", cex.main=5)


allRes <- GenTable(GOdata, classicFisher = resultFisher,classicKS = resultKS, elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 20)

allRes$Term

write.table(allRes, file = "~/Studium/Köln_Biological Science/Master Thesis/Results/core_GO.csv", sep=";", col.names = T,row.names = F)





### SOFTCORE ####
genes_soft<-read.table("~/Studium/Köln_Biological Science/Master Thesis/Results/soft_genes_wo_contig_2.txt", sep="\t")

genes_soft<-genes_soft[order(genes_soft$V1, decreasing = F),]


geneList_soft <- factor(as.integer(geneNames %in% genes_soft))
names(geneList_soft) <- geneNames
str(geneList_soft)

GOdata_soft <- new("topGOdata", ontology = "BP", allGenes = geneList_soft, annot = annFUN.gene2GO, gene2GO = gene2GO_new)

####enrichment test

resultFisher_soft <- runTest(GOdata_soft, algorithm = "classic", statistic = "fisher")
resultFisher_soft

resultKS_soft <- runTest(GOdata_soft, algorithm = "classic", statistic = "ks")
resultKS_soft
resultKS.elim_soft <- runTest(GOdata_soft, algorithm = "elim", statistic = "ks")
resultKS.elim_soft


par(cex=0.2)
showSigOfNodes(GOdata_soft, score(resultKS.elim_soft), firstSigNodes = 5, useInfo ='def',.NO.CHAR = 20, plotFunction = GOplot)

title(main="GO enrichment soft-core genes", cex.main=5)

allRes_soft <- GenTable(GOdata_soft, classicFisher = resultFisher_soft,classicKS = resultKS_soft, elimKS = resultKS.elim_soft, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 20)

allRes_soft$Term

write.table(allRes_soft, file = "~/Studium/Köln_Biological Science/Master Thesis/Results/softcore_GO.csv", sep=";", col.names = T,row.names = F)


#### PRIVATE ####
genes_private <-read.table("~/Studium/Köln_Biological Science/Master Thesis/Results/private_wo_contig_2.txt", sep="\t")

genes_private<-genes_private[order(genes_private$V1, decreasing = F),]

geneList_private <- factor(as.integer(geneNames %in% genes_private))
names(geneList_private) <- geneNames
str(geneList_private)

GOdata_private <- new("topGOdata", ontology = "BP", allGenes = geneList_private, annot = annFUN.gene2GO, gene2GO = gene2GO_new)
summary(GOdata_private)
####enrichment test

resultFisher_private <- runTest(GOdata_private, algorithm = "classic", statistic = "fisher")
resultFisher_private

resultKS_private <- runTest(GOdata_private, algorithm = "classic", statistic = "ks")
resultKS_private 
resultKS.elim_private <- runTest(GOdata_private, algorithm = "elim", statistic = "ks")
resultKS.elim_private

par(cex=0.2)
showSigOfNodes(GOdata_private, score(resultKS_private), firstSigNodes = 5, useInfo ='def',.NO.CHAR = 20, plotFunction = GOplot)

title(main="GO enrichment private genes", cex.main=5)

allRes_private <- GenTable(GOdata_private, classicFisher = resultFisher_private,classicKS = resultKS_private, elimKS = resultKS.elim_private, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)

allRes_private$Term

write.table(allRes_private, file = "~/Studium/Köln_Biological Science/Master Thesis/Results/private_GO.csv", sep=";", col.names = T,row.names = F)

### DISPENSABLE ####
genes_dispensable <-read.table("~/Studium/Köln_Biological Science/Master Thesis/Results/dispensable_genes_wo_contig_2.txt", sep="\t")

genes_dispensable<-genes_dispensable[order(genes_dispensable$V1, decreasing = F),]

geneList_dispensable <- factor(as.integer(geneNames %in% genes_dispensable))
names(geneList_dispensable) <- geneNames
str(geneList_dispensable)

GOdata_dispensable <- new("topGOdata", ontology = "BP", allGenes = geneList_dispensable, annot = annFUN.gene2GO, gene2GO = gene2GO_new)
summary(GOdata_dispensable)
####enrichment test

resultFisher_dispensable <- runTest(GOdata_dispensable, algorithm = "classic", statistic = "fisher")
resultFisher_dispensable

resultKS_dispensable <- runTest(GOdata_dispensable, algorithm = "classic", statistic = "ks")
resultKS_dispensable 
resultKS.elim_dispensable <- runTest(GOdata_dispensable, algorithm = "elim", statistic = "ks")
resultKS.elim_dispensable

par(cex=0.2)
showSigOfNodes(GOdata_dispensable, score(resultKS_dispensable), firstSigNodes = 5, useInfo ='def',.NO.CHAR = 20, plotFunction = GOplot)

title(main="GO enrichment dispensable genes", cex.main=5)

allRes_dispensable <- GenTable(GOdata_dispensable, classicFisher = resultFisher_dispensable,classicKS = resultKS_dispensable, elimKS = resultKS.elim_dispensable, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)

allRes_dispensable$Term

write.table(allRes_dispensable, file = "~/Studium/Köln_Biological Science/Master Thesis/Results/dispensable_GO.csv", sep=";", col.names = T,row.names = F)

####################################################################################
#### MERGED PRIVATE + DISPENSABLE ####

genes_merged <-read.table("~/Studium/Köln_Biological Science/Master Thesis/Results/private_dispensable_genes_wo_contig_2.txt", sep="\t")

genes_merged<-genes_merged[order(genes_merged$V1, decreasing = F),]

#OTHER WAY AROUND
geneList_merged <- factor(as.integer(geneNames %in% genes_merged))
names(geneList_merged) <- geneNames
str(geneList)

GOdata_merged <- new("topGOdata", ontology = "BP", allGenes = geneList_merged, annot = annFUN.gene2GO, gene2GO = gene2GO_new)



####enrichment test

resultFisher_merged <- runTest(GOdata_merged, algorithm = "classic", statistic = "fisher")
resultFisher_merged

resultKS_merged <- runTest(GOdata_merged, algorithm = "classic", statistic = "ks")
resultKS_merged 
resultKS.elim_merged <- runTest(GOdata_merged, algorithm = "elim", statistic = "ks")
resultKS.elim_merged

par(cex=0.2)
showSigOfNodes(GOdata_merged, score(resultKS_merged), firstSigNodes = 5, useInfo ='def',.NO.CHAR = 20, plotFunction = GOplot)
title(main="GO enrichment merged genes", cex.main=5)

allRes_merged <- GenTable(GOdata_merged, classicKS = resultKS_merged,orderBy = "KS", ranksOf = "classicFisher", topNodes = 30)

allRes_merged$Term

write.table(allRes_merged, file = "~/Studium/Köln_Biological Science/Master Thesis/Results/merged_GO.csv", sep=";", col.names = T,row.names = F)


pValue.classic <- score(resultKS_merged)
pValue_disp <- as.data.frame(pValue.classic)
pValue_disp$Orthogroup <-  row.names(pValue_disp)
# Assuming your data frame is named pValue_spain
pValue_disp <- pValue_disp[order(pValue_disp$pValue.classic), ]
pValue_disp <- subset(pValue_disp, select = -Orthogroup)
write.table(pValue_spain,"GO_p_value_disp.txt")
head(pValue_spain)
install.packages("DOSE")
library(DOSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de)

edo2 <- gseDO(geneList)
dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")

library(enrichplot)
barplot(edo, showCategory=20) 

#### SIGNIFICANT OG####
allRes$Term
allRes_private$Term
allRes_dispensable$Term
allRes_merged$GO.ID

##### PLOTS ####
par(mar = c(0, 0, 0, 0))
par(mfrow=c(1,2),cex=0.35)
par(cex=0.2)
showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 4, useInfo ='def',.NO.CHAR = 38, plotFunction = GOplot)
title(main="GO enrichment core genes", cex.main=3.5)

showSigOfNodes(GOdata_soft, score(resultKS_soft), firstSigNodes = 4, useInfo ='def',.NO.CHAR = 20, plotFunction = GOplot)
title(main="GO enrichment softcore genes", cex.main=5)

showSigOfNodes(GOdata_dispensable, score(resultKS_dispensable), firstSigNodes = 4, useInfo ='def',.NO.CHAR = 20, plotFunction = GOplot)
title(main="GO enrichment dispensable genes", cex.main=2)

showSigOfNodes(GOdata_private, score(resultKS_private), firstSigNodes = 4, useInfo ='def',.NO.CHAR = 38, plotFunction = GOplot)
title(main="GO enrichment variable genes", cex.main=5)
allRes_private$Term

showSigOfNodes(GOdata_merged, score(resultKS_merged), firstSigNodes = 4, useInfo ='def',.NO.CHAR = 150, plotFunction = GOplot)
title(main="GO enrichment merged genes", cex.main=5)

tiff("GO_enrichment_core.png")
showSigOfNodes(GOdata, score(resultKS), firstSigNodes = 4, useInfo ='def',.NO.CHAR = 30, plotFunction = GOplot)
title(main="GO enrichment core genes", cex.main=5)
dev.off()


printGraph(GOdata,  score(resultKS), firstSigNodes = 4, fn.prefix = "tGO",pdfSW = TRUE)
printGraph(GOdata, score(resultKS), firstSigNodes = 4)








### 1.1.3) GO enrichment for unqiue gene sets per region####
#Scandinavia
scandinavia_genes <-read.table("~/Studium/Köln_Biological Science/Master Thesis/Results/scandinavia_genes_unique.txt", sep="\t")
scandinavia_genes<-scandinavia_genes[order(scandinavia_genes$V1, decreasing = F),]

#OTHER WAY AROUND
geneList_scandinavia <- factor(as.integer(geneNames %in% scandinavia_genes))
names(geneList_scandinavia) <- geneNames
str(geneList)

GOdata_scandinavia <- new("topGOdata", ontology = "BP", allGenes = geneList_scandinavia, annot = annFUN.gene2GO, gene2GO = gene2GO_new)



####enrichment test

resultFisher_scandinavia <- runTest(GOdata_scandinavia, algorithm = "classic", statistic = "fisher")
resultFisher_scandinavia

resultKS_scandinavia <- runTest(GOdata_scandinavia, algorithm = "classic", statistic = "ks")
resultKS_scandinavia 
resultKS.elim_scandinavia <- runTest(GOdata_scandinavia, algorithm = "elim", statistic = "ks")
resultKS.elim_scandinavia

par
par(cex=0.3)
showSigOfNodes(GOdata_scandinavia, score(resultKS_scandinavia), firstSigNodes = 5, useInfo ='def',.NO.CHAR = 20, plotFunction = GOplot)
title(main="GO enrichment scandinavian genes", cex=8)

allRes_scandinvia <- GenTable(GOdata_scandinavia, classicFisher = resultFisher_scandinavia,classicKS = resultKS_scandinavia, elimKS = resultKS.elim_scandinavia, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 20)
allRes_scandinvia$Term
write.table(allRes_scandinvia, file = "~/Studium/Köln_Biological Science/Master Thesis/Results/Scandinavian_GO.csv", sep=";", col.names = T,row.names = F)



##SPAIN GO ENRICHMENT 
spain_genes <-read.table("~/Studium/Köln_Biological Science/Master Thesis/Results/spain_genes.txt", sep="\t")

spain_genes<-spain_genes[order(spain_genes$V1, decreasing = F),]
head(scandinavia_genes)
#OTHER WAY AROUND
geneList_spain <- factor(as.integer(geneNames %in% spain_genes))
names(geneList_spain) <- geneNames


GOdata_spain <- new("topGOdata", ontology = "BP", allGenes = geneList_spain, annot = annFUN.gene2GO, gene2GO = gene2GO_new)



####enrichment test

resultFisher_spain <- runTest(GOdata_spain, algorithm = "classic", statistic = "fisher")
resultFisher_spain

resultKS_spain <- runTest(GOdata_spain, algorithm = "classic", statistic = "ks")
resultKS_spain 
resultKS.elim_spain <- runTest(GOdata_spain, algorithm = "elim", statistic = "ks")
resultKS.elim_spain

par
par(cex=0.4)
showSigOfNodes(GOdata_spain, score(resultKS_spain), firstSigNodes = 5, useInfo ='def',.NO.CHAR = 10, plotFunction = GOplot)
title(main="GO enrichment spanish genes", cex=5)

allRes_spain <- GenTable(GOdata_scandinavia, classicFisher = resultFisher_spain,classicKS = resultKS_spain, elimKS = resultKS.elim_spain, orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 20)
allRes_spain$Term
write.table(allRes_spain, file = "~/Studium/Köln_Biological Science/Master Thesis/Results/Spain_GO.csv", sep=";", col.names = T,row.names = F)

pValue.classic <- score(resultKS_spain)
pValue_spain <- as.data.frame(pValue.classic,row.names = F)
pValue_spain$Orthogroup <-  row.names(pValue_spain)
# Assuming your data frame is named pValue_spain
pValue_spain <- pValue_spain[order(pValue_spain$pValue.classic), ]
pValue_spain <- subset(pValue_spain, select = -Orthogroup)
write.table(pValue_spain,"GO_p_value_spain.txt")
head(pValue_spain)
