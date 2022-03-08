################ TOPGO ####################################################
library("topGO")
library("GOplot")
#library(ALL)
library("xtable")
library("drake")
pkgconfig::set_config("drake::strings_in_dots" = "literals")

## if you need to rerun everything
## drake::clean(destroy = TRUE)

## Should be modified to correspond to your system
source("src/functions.R")

#GENOME STUFF
geneID2GO <- readMappings(file = "data/all_annotations_justGO.txt")
GO2geneID <- inverseList(geneID2GO)
geneNames <- names(geneID2GO)
deGenes <- read.table("data/both.transcripts.de.txt")
resEdgeRtopTags <- read.table("data/counts52KPS_de_results.txt")
dat.filtered <- read.table("data/datFileredCounts.txt", header = T)
head(geneNames)

plan <- drake_plan(
  name_list = c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch"),
  table = as.factor(geneNames) %in% deGenes$V1,
  int_table = as.integer(table),
  int_fac_table = factor(int_table),
  fac_table = rename(table = int_fac_table, geneNames = geneNames),
  
  GOdata.BP = new("topGOdata", ontology = "BP", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO),
  GOdata.MF = new("topGOdata", ontology = "MF", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO),
  GOdata.CC = new("topGOdata", ontology = "CC", allGenes = fac_table, annot = annFUN.gene2GO, gene2GO = geneID2GO),
  resultWeight01.BP = runTest(GOdata.BP, statistic = "fisher"),
  resultWeight01.MF = runTest(GOdata.MF, statistic = "fisher"),
  resultWeight01.CC = runTest(GOdata.CC, statistic = "fisher"),
  
  allRes.BP1 = GenTable(GOdata.BP, weight01_pval=resultWeight01.BP, orderBy = "weight01", ranksOf = "weight01",topNodes = 500),
  allRes.BP2 = cbind(allRes.BP1,"BP"),
  allRes.BP = change_names(data = allRes.BP2, name_list = name_list),

  allRes.MF1 = GenTable(GOdata.MF, weight01_pval=resultWeight01.MF, orderBy = "weight01", ranksOf = "weight01",topNodes = 500),
  allRes.MF2 = cbind(allRes.MF1,"MF"),
  allRes.MF = change_names(data = allRes.MF2, name_list = name_list),
  
  allRes.CC1 = GenTable(GOdata.CC, weight01_pval=resultWeight01.CC, orderBy = "weight01", ranksOf = "weight01",topNodes = 500),
  allRes.CC2 = cbind(allRes.CC1,"CC"),
  allRes.CC = change_names(data = allRes.CC2, name_list = name_list),
  #colnames(allRes.CC) = c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch"),
  
  allRes1 = rbind(allRes.BP,allRes.MF),
  allRes2 = rbind(allRes1, allRes.CC),
  allRes = allRes2[allRes2$weight01_pval<=0.1,],
  
  
  allGO.BP = genesInTerm(GOdata.BP),
  allGO.MF = genesInTerm(GOdata.MF),
  allGO.CC = genesInTerm(GOdata.CC),
  allGO = c(allGO.BP, allGO.MF, allGO.CC),
  
  SAM_ANOTATION = lapply(allGO,function(x) x[x %in%  deGenes$V1]),
  enriched_go_with_my_genes = lapply(SAM_ANOTATION[allRes[,1]], paste0, collapse = ", "),
  enriched_go_with_my_genes.list = attach_enriched_go_genes(enriched_go_with_my_genes),
  go.dataframe = data.frame("Category" = allRes$branch, "ID" = allRes$GO.ID, "Term" = allRes$Term, "Genes" = as.vector(enriched_go_with_my_genes.list), "adj_pval" = as.numeric(sub(",", ".", allRes$weight01_pval, fixed = TRUE))),
  EC.genelist = data.frame("ID" = rownames(resEdgeRtopTags), "logFC" = resEdgeRtopTags$logFC, "AveExpr" = rowMeans(dat.filtered[rownames(resEdgeRtopTags),]), "P.Value" = resEdgeRtopTags$PValue, "adj.P.Val" = resEdgeRtopTags$F),
  circ = circle_dat(go.dataframe, EC.genelist)
  
)

make(plan)
# Now, your data is saved in the .drake folder as .rds files. 
# Like this, you will only change the variables needed when the script is modified. Usefull for large pipelines.

# Retrieve the circ variable from the .drake folder.
circ <- readd(circ)

# Work with the datframe.
reduced_circ <- reduce_overlap(circ, overlap = 1.00)
reduced_circ_counts <- subset(reduced_circ, count > 2)

GOplot::GOBar(subset(reduced_circ, category == 'BP'), zsc.col = c("blue", "white", "red"))
GOplot::GOBar(subset(reduced_circ_counts, category == 'BP'), zsc.col = c("blue", "white", "red"))

GOplot::GOCircle(subset(circ, category == 'BP' | category == 'MF'), lfc.col = c("blue", "red"), zsc.col = c("blue", "white", "red"))

reduced_circBP <- reduce_overlap(subset(circ, category == "BP"), overlap = 0.6)
reduced_circMF <- reduce_overlap(subset(circ, category == "MF"), overlap = 0.6)
reduced_circCC <- reduce_overlap(subset(circ, category == "CC"), overlap = 0.6)


reduced_circ <- reduce_overlap(circ.BP.CC.MF, overlap = 0.8)
GOBubble(reduced_circ, labels = 4.5, display = 'multiple', table.legend = T)

