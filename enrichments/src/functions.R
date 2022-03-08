topGoAnalysis <- function(ontology, allGenes, gene2GO, nbNodes){
  GOdata <- new("topGOdata", ontology = ontology, allGenes = allGenes, annot = annFUN.gene2GO, gene2GO = gene2GO)
  resultWeight01 <- runTest(GOdata, statistic = "fisher")
  allRes <- GenTable(GOdata, weight01_pval=resultWeight01, orderBy = "weight01", ranksOf = "weight01",topNodes = nbNodes)
  allRes <- cbind(allRes,"BP")
  colnames(allRes) <- c("GO.ID","Term","Annotated","Significant","Expected","weight01_pval", "branch")
  return(allRes)
}

change_names <- function(data, name_list){
  colnames(data) <- name_list
  return(data)
}

rename <- function(table, geneNames){
  names(table) <- geneNames
  return(table)
}

attach_enriched_go_genes <- function(enriched_go_with_my_genes){
  enriched_go_with_my_genes.list = c()
  for (i in 1:length(enriched_go_with_my_genes)){
    enriched_go_with_my_genes.list = c(enriched_go_with_my_genes.list, enriched_go_with_my_genes[[i]])
  }
  return(enriched_go_with_my_genes.list)
}
circle_dat <- function(terms, genes){
  
  colnames(terms) <- tolower(colnames(terms))
  terms$genes <- toupper(terms$genes)
  genes$ID <- toupper(genes$ID)
  tgenes <- strsplit(as.vector(terms$genes), ', ')
  if (length(tgenes[[1]]) == 1) tgenes <- strsplit(as.vector(terms$genes), ',')
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  logFC <- sapply(unlist(tgenes), function(x) genes$logFC[match(x, genes$ID)])
  if(class(logFC) == 'factor'){
    logFC <- gsub(",", ".", gsub("\\.", "", logFC))
    logFC <- as.numeric(logFC)
  }
  s <- 1; zsc <- c()
  for (c in 1:length(count)){
    value <- 0
    e <- s + count[c] - 1
    value <- logFC[s:e]
    #value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 1, -1))
    zsc <- c(zsc, sum(value, na.rm = F) / sqrt(count[c])) #### HERE : na.rm = TRUE, takes all the genes into account !
    s <- e + 1
  }
  if (is.null(terms$id)){
    df <- data.frame(category = rep(as.character(terms$category), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }else{
    df <- data.frame(category = rep(as.character(terms$category), count), ID = rep(as.character(terms$id), count), term = rep(as.character(terms$term), count),
                     count = rep(count, count), genes = as.character(unlist(tgenes)), logFC = logFC, adj_pval = rep(terms$adj_pval, count),
                     zscore = rep(zsc, count), stringsAsFactors = FALSE)
  }
  return(df)
}

#SAM_ANOTATION = lapply(allGO,function(x) x[x %in% myInterestingGenes] )
#enriched_go_with_my_genes <- lapply(SAM_ANOTATION[allRes.sign[,1]], paste0, collapse = ", ")

#enriched_go_with_my_genes.list <- c()
#for (i in 1:length(enriched_go_with_my_genes)){
#  enriched_go_with_my_genes.list <- c(enriched_go_with_my_genes.list, enriched_go_with_my_genes[[i]])
#}