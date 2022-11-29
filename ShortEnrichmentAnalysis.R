# Get Genelist and expression marix
topTabAvsB <- read.table ("datasets/Top_AvsB.csv2", 
                          head=T, sep=";", dec=",", row.names=1)
expresAvsB <- read.table ("datasets/expres_AvsB.csv2", 
                          head=T, sep=";", dec=",", row.names=1)

# Define Gene list using arbitrary though reasonable cutoffs
probesUniverse <- rownames(topTabAvsB)
whichGenesInTop<- topTab["adj.P.Val"]<0.05 & topTab["logFC"] > 1

# Annotate Gene Universe and Gene list
entrezUniverse<- select(hgu133a.db, probesUniverse, "ENTREZID")
entrezUniverse <- entrezUniverse$ENTREZID
topGenes <-   entrezUniverse[whichGenesInTop]
head(topGenes)

# Remove possible duplicates

topGenes <- topGenes[!duplicated(topGenes)]
entrezUniverse <- entrezUniverse[!duplicated(entrezUniverse)]

# ORA
library(GOstats)
GOparams = new("GOHyperGParams",
               geneIds=topGenes, universeGeneIds=entrezUniverse,
               annotation="hgu133a.db", ontology="BP",
               pvalueCutoff=0.001)
GOhyper = hyperGTest(GOparams)
head(summary(GOhyper))
