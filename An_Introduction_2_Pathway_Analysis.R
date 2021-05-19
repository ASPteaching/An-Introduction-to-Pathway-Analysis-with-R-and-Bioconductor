## ----global_options, include=FALSE-------------------------------------------------------------
knitr::opts_chunk$set(fig.width=12, fig.height=8, cache=FALSE,
                      echo=TRUE, warning=FALSE, message=FALSE, results ='markup')
options(warn=-1)


## ----installPackages, eval=FALSE---------------------------------------------------------------
## installifnot <- function (packageName){
##  if (!(require(packageName, character.only=TRUE))) {
##     install.packages(packageName)
##   }else{
##     detach(paste ("package", packageName, sep=":"), character.only=TRUE)
##   }
## }
## bioCifnot <- function (packageName){
##  if (!(require(packageName, character.only=TRUE))) {
##     BiocManager::install(packageName)
##  }else{
##   detach(paste ("package", packageName, sep=":"), character.only=TRUE)
## }
## }
## installifnot("knitr")
## installifnot("XML") # May yield problems if some libraries (xml2-config) not available in linux
## bioCifnot ("org.Hs.eg.db")
## bioCifnot ("hgu133a.db")
## bioCifnot ("GO.db")
## bioCifnot ("GOstats")
## bioCifnot ("topGO")
## bioCifnot ("annotate")
## bioCifnot ("Rgraphviz")
## bioCifnot ("clusterProfiler")


## ----readData1---------------------------------------------------------------------------------
inputDir="datasets"
topTabAvsB <- read.table (file.path(inputDir, "Top_AvsB.csv2"), head=T, sep=";", dec=",", row.names=1)
expresAvsB <- read.table (file.path(inputDir, "expres_AvsB.csv2"), head=T, sep=";", dec=",", row.names=1)
comparisonName <- "AvsB"
dim(topTabAvsB); head(topTabAvsB)
dim(expresAvsB); head(expresAvsB)


## ----probes------------------------------------------------------------------------------------
myProbes <- rownames(expresAvsB)
head(myProbes)


## ----mappings----------------------------------------------------------------------------------
library(hgu133a.db)
keytypes(hgu133a.db)


## ----------------------------------------------------------------------------------------------
geneAnots <- select(hgu133a.db, myProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
head(geneAnots)


## ----------------------------------------------------------------------------------------------
selected<- topTabAvsB[,"adj.P.Val"]<0.05 & topTabAvsB[,"logFC"] > 1
sum(selected)
selectedTopTab <- topTabAvsB[selected,]
selectedProbes <- rownames(selectedTopTab)
selectedAnots <-  select(hgu133a.db, selectedProbes, c("SYMBOL", "ENTREZID", "GENENAME"))
selectedInfo <- cbind(selectedAnots, selectedTopTab)
write.csv2(selectedInfo, file="selectedTopTab_AvsB.csv2")


## ----------------------------------------------------------------------------------------------
require(GOstats)
require(GO.db)
require(hgu133a.db); 
require(annotate) # Loads the required libraries.


## ----top25-------------------------------------------------------------------------------------
probes <- rownames(expresAvsB)[1:5]


## ----------------------------------------------------------------------------------------------
charFromSelect <- function(df,posNames=1, posValues=2){
  res <- df[,posValues]
  names(res) <- df[,posNames]
  return(res)
}
  
require(annotate)
geneIDs <-  select(hgu133a.db, probes, c("ENTREZID", "SYMBOL"))
entrezs <-  charFromSelect(geneIDs, 1, 2)
simbols <-  charFromSelect(geneIDs, 1, 3) 


## ----------------------------------------------------------------------------------------------
GOAcc<-mget(probes,env=hgu133aGO)
GOAcc[[1]][1:5]


## ----GOtable-----------------------------------------------------------------------------------
# % WANING
#  The previous chunk can be substituted by the followink code chunk, shorter and more efficient.
library(hgu133a.db)
keytypes(hgu133a.db)
res <- select(hgu133a.db, keys=probes, keytype = "PROBEID", columns = c("ENTREZID", "SYMBOL","ONTOLOGY"))
res1 <- select(hgu133a.db, keys=probes, keytype = "PROBEID",  columns = c("ENTREZID", "SYMBOL","GO"))


## ----echo=TRUE---------------------------------------------------------------------------------
print(head(res, n=10))
print(head(res1, n=10))


## ----------------------------------------------------------------------------------------------
oneTerm <- "GO:0006338"
oneParent<-get(oneTerm, GOBPPARENTS)   # the vector of its parent terms in the BP ontology.
oneParent
oneChildren<-get(oneTerm, GOBPCHILDREN) # the vector of its children terms in the BP ontolog
oneChildren
oneOffspring<-get(oneTerm, GOBPOFFSPRING) # the vector of its offspring terms in the BP ontology.
oneOffspring
oneChildren %in% oneOffspring


## ----------------------------------------------------------------------------------------------
require(org.Hs.eg.db) # loads the library
myEIDs3 <-entrezs[1:3] # Create vecotor of input Entrez IDs
myGO <- unlist(org.Hs.egGO[[as.character(myEIDs3[1])]])
myGO_All <- mget(myEIDs3, org.Hs.egGO)
GOgenes <- org.Hs.egGO2ALLEGS[[myGO[1]]]
GOgenes_All <- mget(myGO[1], org.Hs.egGO2ALLEGS)


## ----------------------------------------------------------------------------------------------
require(hgu133a.db)
topTab <- topTabAvsB 
entrezUniverse = unlist(mget(rownames(topTab), hgu133aENTREZID, ifnotfound=NA)) 
whichGenes<- topTab["adj.P.Val"]<0.05 & topTab["logFC"] > 1
sum(whichGenes)
topGenes <-   entrezUniverse[whichGenes]
allMyGenes <-topTab$adj.P.Val
names(allMyGenes)<-rownames(topTab)


## ----------------------------------------------------------------------------------------------
require(topGO)
data(geneList) # adds function "topDiffGenes"
myGOData <- new("topGOdata", ontology="BP", 
                allGenes=allMyGenes,
                geneSel=topDiffGenes, nodeSize=10,  
                annot= annFUN.db, affyLib="hgu133a.db")

Myenrichment_Fisher <- runTest(myGOData, algorithm= "classic", statistic="fisher")
Myenrichment_Fisher

head(score(Myenrichment_Fisher), 25) # Displays p values for every GO term
geneData(Myenrichment_Fisher) # A table showing Medata data for enrichment


## ----------------------------------------------------------------------------------------------
Myenrichment_KS <- runTest(myGOData, algorithm= "classic", statistic="ks")


## ----------------------------------------------------------------------------------------------
enrich_table <-GenTable(myGOData, classicFisher=Myenrichment_Fisher,topNodes = 20)
adjustedEnrichP <- cbind(enrich_table, adjP=p.adjust(enrich_table$classicFisher, method = "BH"))
head(adjustedEnrichP, n=25) # get the enrichment results as table


## ----------------------------------------------------------------------------------------------
showSigOfNodes(myGOData, score(Myenrichment_Fisher), firstSigNodes=5, useInfo="all") # Plot the enrichment GO graph
gostat <- termStat(myGOData, names(score(Myenrichment_Fisher)))
plot(score(Myenrichment_Fisher), score(Myenrichment_KS)[names(score(Myenrichment_Fisher))], xlab="P values Fisher test", ylab="P values KS test", cex=(gostat$Annotated/max(gostat$Annotated))*4, col=heat.colors(gostat$Significant))
print(showGroupDensity(myGOData, enrich_table[1, "GO.ID"], ranks=TRUE))


## ----------------------------------------------------------------------------------------------
entrezIDs <- AnnotationDbi::select(hgu133a.db, rownames(topTabAvsB), c("ENTREZID"))
entrezIDs <- charFromSelect(entrezIDs)
geneList <- cbind(topTabAvsB, ENTREZ =entrezIDs)

#ordenem per logFC absolut per eliminar els duplicats amb menor logFC absolut
geneList <- geneList[order(abs(geneList$logFC), decreasing=T),]
geneList <- geneList[ !duplicated(geneList$ENTREZ), ]  ### Keep highest
    #tornem a ordenar per logFC per fer el GSEA
geneList <- geneList[order(geneList$logFC, decreasing=T),]
genesVector <- geneList$logFC
names(genesVector) <- geneList$ENTREZ


#fixem seed per reproduibilitat dels resultats
#  set.seed(123)
library(clusterProfiler)
gseResulti <- gseKEGG(geneList = genesVector,
                      organism = "hsa",
                      keyType = "kegg",
                      exponent = 1,
                      minGSSize = 10,maxGSSize = 500,
                      pvalueCutoff = 0.05,pAdjustMethod = "BH",
                      # nPerm = 10000, #augmentem permutacions a 10000
                      verbose = TRUE,
                      use_internal_data = FALSE,
                      seed = TRUE,
                      eps=0,
                      by = "fgsea"
                )

# keggResultsList[[i]] <- gseResulti


## ----results='asis'----------------------------------------------------------------------------
library(kableExtra)
gsea.result <- setReadable(gseResulti, OrgDb = org.Hs.eg.db, keyType ="ENTREZID" )

gsea.result.df <- as.data.frame(gsea.result)
print(kable(gsea.result.df[,c("Description","setSize","NES","p.adjust")])%>% scroll_box(height = "500px"))
  


## ----eval=TRUE---------------------------------------------------------------------------------
library(ggplot2)
# for (i in 1:length(files)){
#   cat("\nComparison: ", namesC[i],"\n")
   cat("DOTPLOT\n")
#   if(nrow(keggResultsList[[i]]) > 0){
 if(nrow(gseResulti) > 0){
   p<- dotplot(gseResulti, showCategory = 20, font.size = 15,
            title =paste("Enriched Pathways\n", comparisonName ,
            split=".sign") + facet_grid(.~.sign))
   plot(p)
   cat("\nENRICHMENT MAP\n")
   em<- emapplot(gseResulti)
   plot(em)
   #guardem en pdf
   pdf(file = paste0("KEGGplots.",comparisonName,".pdf"), 
                        width = 14, height = 14)
   print(p)
   print(em)
   dev.off()
   }else{
      cat("\nNo enriched terms found\n")
 }


