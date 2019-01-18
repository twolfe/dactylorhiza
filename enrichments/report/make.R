################ TOPGO ####################################################
library("topGO")
library("GOplot")
#library(ALL)
library("xtable")
library("drake")
pkgconfig::set_config("drake::strings_in_dots" = "literals")

## if you need to rerun everything
## drake::clean(destroy = TRUE)
source("/home/thomas/Documents/phd/scripts/enrichments/src/functions.R")
## create softlinks for data
## ---------------------------

#normSet <- readRDS("/home/thomas/Documents/phd/scripts/enrichments/data/normSet.rds")
### FULL SET OF NON FILTERED TRANSCRIPTS ##########################################
#geneID2GO <- readMappings(file = "/data/phdData/orchis/Orchis.transcriptome.go.annotation.topGO.map")
### REDUCED SET OF FILTERED TRANSCRIPTS ##########################################
#geneID2GO <- readMappings(file = "/data/phdData/orchis/counts/de.featureCounts/retained.after.filtering.go.map.txt")
geneID2GO <- readMappings(file = "/home/thomas/Documents/phd/scripts/enrichments/data/retained.genes.justGO.trinotate.txt")
deGenes <- read.table("/home/thomas/Documents/phd/research/HEB/HEB/all.HEB.res", header = T)
deGenes <- deGenes$gene[deGenes$color == "orange"]

deGenes <- read.table("/home/thomas/Documents/phd/scripts/enrichments/data/both.transcripts.de.txt")
both.go <- read.table("/home/thomas/Documents/phd/scripts/enrichments/data/both.blast2go.enrichments.txt"
                      , header = T, sep = "\t")
resEdgeRtopTags <- readRDS("/home/thomas/Documents/phd/scripts/enrichments/data/all.results.rds")
dat.filtered <- read.table("/home/thomas/Documents/phd/scripts/enrichments/data/datFilteredCounts.txt")
deGenes <- as.vector(deGenes$V1)
### ON THE RUN ANALYSIS
#deGenes <- FDRres[FDRres$FDR <= 0.05,]
#env_gen_analysis<-readRDS("/home/botanik/Documents/phd/research/environmental_study/deEnviroCorrelations/2018/env_gen_analysis.RDS")
#finalenvmatrixrep <- readRDS("/home/botanik/Documents/phd/research/environmental_study/deEnviroCorrelations/2018/finalenvmatrixrep.RDS")
#treeCover <- readRDS("/home/botanik/Documents/phd/research/environmental_study/deEnviroCorrelations/2018/treecoverageDactpolyploids.RDS")
GO2geneID <- inverseList(geneID2GO)
geneNames <- names(geneID2GO)
head(geneNames)
#resEdgeRtopTags is the table found after running cleanDE.R
geneList <- factor(as.integer(geneNames %in% rownames(myInterestingGenes)))

plan <- drake_plan(
  name_list = c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch"),
  table = geneNames %in% deGenes,
  int_table = as.integer(table),
  int_fac_table = factor(int_table),
  fac_table = rename(table = int_fac_table, geneNames = geneNames),
  
  GOdata.BP = new("topGOdata", ontology = "BP", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO),
  GOdata.MF = new("topGOdata", ontology = "MF", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO),
  GOdata.CC = new("topGOdata", ontology = "CC", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO),
  resultWeight01.BP = runTest(GOdata.BP, statistic = "fisher"),
  resultWeight01.MF = runTest(GOdata.MF, statistic = "fisher"),
  resultWeight01.CC = runTest(GOdata.CC, statistic = "fisher"),
  
  allRes.BP1 = GenTable(GOdata.BP, weight01_pval=resultWeight01.BP, orderBy = "weight01", ranksOf = "weight01",topNodes = 100),
  allRes.BP2 = cbind(allRes.BP1,"BP"),
  allRes.BP = change_names(data = allRes.BP2, name_list = name_list),

  allRes.MF1 = GenTable(GOdata.MF, weight01_pval=resultWeight01.MF, orderBy = "weight01", ranksOf = "weight01",topNodes = 100),
  allRes.MF2 = cbind(allRes.MF1,"MF"),
  allRes.MF = change_names(data = allRes.MF2, name_list = name_list),
  
  allRes.CC1 = GenTable(GOdata.CC, weight01_pval=resultWeight01.CC, orderBy = "weight01", ranksOf = "weight01",topNodes = 100),
  allRes.CC2 = cbind(allRes.CC1,"CC"),
  allRes.CC = change_names(data = allRes.CC2, name_list = name_list),
  #colnames(allRes.CC) = c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch"),
  
  allRes1 = rbind(allRes.BP,allRes.MF),
  allRes = rbind(allRes1, allRes.CC),
  
  allGO.BP = genesInTerm(GOdata.BP),
  allGO.MF = genesInTerm(GOdata.MF),
  allGO.CC = genesInTerm(GOdata.CC),
  allGO = c(allGO.BP, allGO.MF, allGO.CC),
  
  SAM_ANOTATION = lapply(allGO,function(x) x[x %in%  deGenes]),
  enriched_go_with_my_genes = lapply(SAM_ANOTATION[allRes[,1]], paste0, collapse = ", "),
  enriched_go_with_my_genes.list = attach_enriched_go_genes(enriched_go_with_my_genes),
  go.dataframe = data.frame("Category" = allRes$branch, "ID" = allRes$GO.ID, "Term" = allRes$Term, "Genes" = as.vector(enriched_go_with_my_genes.list), "adj_pval" = as.numeric(sub(",", ".", allRes$weight01_pval, fixed = TRUE))),
  EC.genelist = data.frame("ID" = rownames(resEdgeRtopTags$table), "logFC" = resEdgeRtopTags$table$logFC, "AveExpr" = rowMeans(dat.filtered[rownames(resEdgeRtopTags$table),]), "P.Value" = resEdgeRtopTags$table$PValue, "adj.P.Val" = resEdgeRtopTags$table$F),
  circ = circle_dat(go.dataframe, EC.genelist)
  #BP = topGoAnalysis("BP", fac_table, geneID2GO, 100),
  #MF = topGoAnalysis("MF", fac_table, geneID2GO, 100),
  #CC = topGoAnalysis("CC", fac_table, geneID2GO, 100),
  #allRes1 = rbind(BP, MF),
  #allRes2 = rbind(allRes1, CC),
  #allRes_scientific = format(allRes2, scientific = TRUE),
  #allRes_scientific_weight01_pval = allRes_scientific$weight01_pval,
  #allRes1_weight01_pval = as.double(allRes_scientific_weight01_pval),
  #allRes_sign = allRes1[allRes1$weight01_pval <= 0.05,]
)
make(plan)

enrichments <- list()
enrichmentsDE <- list()
myInterestingDeGenes <- list()
for (i in names(env_gen_analysis)){
  print(i[1])
  tryCatch({
  myInterestingGenes <- env_gen_analysis[[i]]
  #myInterestingGenes <- rownames(myInterestingGenes$table[myInterestingGenes$table$FDR <= 0.05,])
  ########## ENRICHMENTS FOR NON_DE GENES ################################################
  myInterestingGenes <- myInterestingGenes$table[myInterestingGenes$table$FDR <= 0.05,]
  geneList <- factor(as.integer(geneNames %in% rownames(myInterestingGenes)))
  names(geneList) <- geneNames
  str(geneList)
  GOdata.MF <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  GOdata.BP <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  GOdata.CC <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  resultWeight01.BP <- runTest(GOdata.BP, statistic = "fisher")
  resultWeight01.MF <- runTest(GOdata.MF, statistic = "fisher")
  resultWeight01.CC <- runTest(GOdata.CC, statistic = "fisher")
  allRes.BP <- GenTable(GOdata.BP, weight01_pval=resultWeight01.BP, orderBy = "weight01", ranksOf = "weight01",topNodes = 100)
  allRes.BP <- cbind(allRes.BP,"BP")
  colnames(allRes.BP) <- c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")
  
  allRes.MF <- GenTable(GOdata.MF, weight01_pval=resultWeight01.MF, orderBy = "weight01", ranksOf = "weight01",topNodes = 100)
  allRes.MF <- cbind(allRes.MF,"MF")
  colnames(allRes.MF) <- c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")
  
  allRes.CC <- GenTable(GOdata.CC, weight01_pval=resultWeight01.CC, orderBy = "weight01", ranksOf = "weight01",topNodes = 100)
  allRes.CC <- cbind(allRes.CC,"CC")
  colnames(allRes.CC) <- c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")
  
  allRes <- rbind(allRes.BP,allRes.MF)
  allRes <- rbind(allRes, allRes.CC)
  enrichments[[i]] <- allRes
  
  ############### ENRICHMENTS FOR DE GENES ##########################
  deCorrGenes <- myInterestingGenes[rownames(myInterestingGenes) %in% deGenes,]
  myInterestingDeGenes[[i]] <- deCorrGenes
  geneList <- factor(as.integer(geneNames %in% rownames(deCorrGenes)))
  names(geneList) <- geneNames
  str(geneList)
  
  geneList <- factor(as.integer(geneNames %in% deGenes))
  names(geneList) <- geneNames
  topGoAnalysis("BP", geneList, geneID2GO, 100)
  
  GOdata.MF <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  GOdata.BP <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  GOdata.CC <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
  resultWeight01.BP <- runTest(GOdata.BP, statistic = "fisher")
  resultWeight01.MF <- runTest(GOdata.MF, statistic = "fisher")
  resultWeight01.CC <- runTest(GOdata.CC, statistic = "fisher")
  allRes.BP <- GenTable(GOdata.BP, weight01_pval=resultWeight01.BP, orderBy = "weight01", ranksOf = "weight01",topNodes = 100)
  allRes.BP <- cbind(allRes.BP,"BP")
  colnames(allRes.BP) <- c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")
  
  allRes.MF <- GenTable(GOdata.MF, weight01_pval=resultWeight01.MF, orderBy = "weight01", ranksOf = "weight01",topNodes = 100)
  allRes.MF <- cbind(allRes.MF,"MF")
  colnames(allRes.MF) <- c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")
  
  allRes.CC <- GenTable(GOdata.CC, weight01_pval=resultWeight01.CC, orderBy = "weight01", ranksOf = "weight01",topNodes = 100)
  allRes.CC <- cbind(allRes.CC,"CC")
  colnames(allRes.CC) <- c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")
  
  allRes <- rbind(allRes.BP,allRes.MF)
  allRes <- rbind(allRes, allRes.CC)
  
  enrichmentsDE[[i]] <- allRes
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
myInterestingGenes <- env_gen_analysis$N.NH4[env_gen_analysis$N.NH4$table$FDR <= 0.05,]

GO2geneID.de <- inverseList(geneID2GO.de)

geneList <- factor(as.integer(geneNames %in% myInterestingGenes))

names(geneList) <- geneNames
str(geneList)

GOdata.MF <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.BP <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata.CC <- new("topGOdata", ontology = "CC", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)

#TESTS
#resultFis.BP <- runTest(GOdata.BP, algorithm = "classic", statistic = "fisher")
#resultKS.BP <- runTest(GOdata.BP, algorithm = "elim", statistic = "fisher")
#resultWeight.BP <- runTest(GOdata.BP, algorithm = "weight", statistic = "fisher")
resultWeight01.BP <- runTest(GOdata.BP, statistic = "fisher")
resultWeight01.MF <- runTest(GOdata.MF, statistic = "fisher")
resultWeight01.CC <- runTest(GOdata.CC, statistic = "fisher")

allRes.BP <- GenTable(GOdata.BP, weight01_pval=resultWeight01.BP, orderBy = "weight01", ranksOf = "weight01",topNodes = 100)
allRes.BP <- cbind(allRes.BP,"BP")
colnames(allRes.BP) <- c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")

allRes.MF <- GenTable(GOdata.MF, weight01_pval=resultWeight01.MF, orderBy = "weight01", ranksOf = "weight01",topNodes = 100)
allRes.MF <- cbind(allRes.MF,"MF")
colnames(allRes.MF) <- c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")

allRes.CC <- GenTable(GOdata.CC, weight01_pval=resultWeight01.CC, orderBy = "weight01", ranksOf = "weight01",topNodes = 100)
allRes.CC <- cbind(allRes.CC,"CC")
colnames(allRes.CC) <- c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")

allRes <- rbind(allRes.BP,allRes.MF)
allRes <- rbind(allRes, allRes.CC)

head(allRes, n=100)
x.big <- xtable(allRes, caption = "A \\code{longtable} spanning several pages")
print(x.big, hline.after=c(-1, 0), tabular.environment = "longtable", scalebox = 0.7)
term.genes <- genes(GOdata, allRes$GO.ID)

### Extract significant results 
allRes.scientific <- format(allRes, scientific = TRUE)
allRes$weight01_pval <- as.double(allRes$weight01_pval)
allRes.sign <- allRes[allRes$weight01_pval <= 0.05,]

##### link here
allGO.BP = genesInTerm(GOdata.BP)
allGO.MF = genesInTerm(GOdata.MF)
allGO.CC = genesInTerm(GOdata.CC)
allGO = c(allGO.BP, allGO.MF, allGO.CC)

SAM_ANOTATION = lapply(allRes,function(x) x[x %in% myInterestingGenes])
enriched_go_with_my_genes <- lapply(SAM_ANOTATION[allRes[,1]], paste0, collapse = ", ")

enriched_go_with_my_genes.list <- c()
for (i in 1:length(enriched_go_with_my_genes)){
  enriched_go_with_my_genes.list <- c(enriched_go_with_my_genes.list, enriched_go_with_my_genes[[i]])
}
############### GOPLOT #################################################
library(splitstackshape)

########################################################################
### very important: change this in the function before running it : ####
########################################################################

#go.dataframe <- data.frame("Category" = both.go$GO.Category, "ID" = both.go$GO.ID, "Term" = both.go$GO.Name, "Genes" = both.go$TestSet.Sequences, "adj_pval" = both.go$FDR)
go.dataframe <- data.frame("Category" = allRes$branch, "ID" = allRes$GO.ID, "Term" = allRes$Term, "Genes" = as.vector(enriched_go_with_my_genes.list), "adj_pval" = as.numeric(allRes$weight01_pval))
EC.genelist <- data.frame("ID" = rownames(resEdgeRtopTags$table), "logFC" = resEdgeRtopTags$table$logFC, "AveExpr" = rowMeans(dat.filtered[rownames(resEdgeRtopTags$table),]), "P.Value" = resEdgeRtopTags$table$PValue, "adj.P.Val" = resEdgeRtopTags$table$F)
#resEdgeRtopTags$table[myInterestingGenes,]

circ <- circle_dat(go.dataframe, EC.genelist)
GOplot::GOBar(circ, zsc.col = c("blue", "white", "red"), display = 'multiple')

GOplot::GOCircle(subset(circ, category == 'BP' | category == 'MF'), lfc.col = c("blue", "red"))


circ.BP.CC.MF <- rbind(subset(circ, category == "BP"), subset(circ, category == "CC"), subset(circ, category == "MF"))

reduced_circ <- reduce_overlap(circ.BP.CC.MF, overlap = 0.6)
GOBubble(reduced_circ, labels = 4.5, display = 'multiple', table.legend = T)
reduced_circ[c(1,2,3,4,5,6,7,24,36,38,44,45,54,55,56),]

write.table(circ.day, "circ.day.txt", sep = "\t")
write.table(circ.night, "circ.night.txt", sep = "\t")
