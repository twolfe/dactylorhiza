species <- c()

remove_id <- function(data, remove_list){
  for(ind in remove_list){
    data[ind] <- NULL
  }
  return(data)
}

change_names <- function(data, name_list){
  colnames(data) <- name_list
  return(data)
}

combineRep = function(data){
  data = transform(data,  tB1830=tB1830_1+tB1830_2, tB1833=tB1833_1+tB1833_2, ## ebudensis
                   fP1707=fP1707_1+fP1707_2, ##fuchsii
                   mA1567=mA1567_1+mA1567_2,  mP1722= mP1722_1+mP1722_2, mP1748=mP1748_1+mP1748_2, mS1765=mS1765_1+mS1765_2, ##majalis
                   tB1798=tB1798_1+tB1798_2, tB1805=tB1805_1+tB1805_2, tB1812=tB1812_1+tB1812_2, tS1920=tS1920_1+tS1920_2) ##traunteineri
  data = data[, -grep("_", colnames(data))]
  return(data)
}

# ----+ functions: makeVectors
makeSpeciesVector = function(data){
  indm = with(data, grepl("m", colnames(data))); indm[indm==1] = 2
  indt = with(data, grepl("t", colnames(data))); indt[indt==1] = 3
  inde = with(data, grepl("e", colnames(data))); inde[inde==1] = 3
  indf = with(data, grepl("f", colnames(data))); indf[indf==1] = 5
  indi = with(data, grepl("i", colnames(data))); indi[indi==1] = 6 #modifie here if you want to make ebu => traun
  ind = indm+indt+inde+indf+indi; ind = sub(2, "majalis", ind); ind = sub(3, "traunsteineri", ind); ind = sub(5, "fuchsii", ind); ind = sub(6, "incarnata", ind)
  #species <<- ind
  return(ind)
}
makeGeographiesVector = function(data){
  indA = with(data, grepl("A", colnames(data))); indA[indA==1] = 2
  indP = with(data, grepl("P", colnames(data))); indP[indP==1] = 3
  indB = with(data, grepl("B", colnames(data))); indB[indB==1] = 4
  indS = with(data, grepl("S", colnames(data))); indS[indS==1] = 5
  ind = indA+indP+indB+indS; ind = sub(2, "alps", ind); ind = sub(3, "pyrenees", ind); ind = sub(4, "britain", ind); ind = sub(5, "scandinavia", ind)
  #geographies <<- ind
  return(ind)
}
# Homoeolog function helper
makeOriginsVector = function(data){
  indf = with(data, grepl("P1\\.", colnames(data))); indf[indf==1] = 1
  indi = with(data, grepl("P2\\.", colnames(data))); indi[indi==1] = 2
  ind = indf+indi; ind = sub(1, "fuchsiiSubGenome", ind); ind = sub(2, "incarnataSubGenome", ind)
  return(ind)
}
makeBatchMatrix = function(data){
  batchIndexVector = with(data, grepl("_", colnames(data)))
  batchVec = which(batchIndexVector == TRUE)
  return(matrix(batchVec, byrow = TRUE, ncol = 2))
}
#---- -

## MAPS
geoMap = c("alps"=17, "britain" = 19, "pyrenees" = 15, "scandinavia" = 18)
colorMap = c("ebudensis"="darkgreen","fuchsii"="purple","incarnata"="orange","majalis"="red","traunsteineri"="blue")
##Could add a Batch Map
##ploidies = c(rep("4",4), rep("2", 10), rep("4", 27))
##ploidyMap = c("2"="blue", "4"="red")



## HELPER FUNCTIONS

# ----+ functions: selectFac
selectSpeciesFac = function(ind){
  return(as.factor(species[ind]))
}
selectGeographiesFac = function(ind){
  return(as.factor(geographies[ind]))
}
selectOriginsFac = function(ind){
  return(as.factor(origins[ind]))
}
##selectPloidiesFac = function(ind){
##  return(as.factor(ploidies[ind]))
##}

## FUNCTIONS

# ----+ function:  selectSpecies
selectSpecies = function(data, spec1, spec2, spec3, spec4, spec5){
  index1 = with(data, grepl(spec1, colnames(data)))
  index2 = with(data, grepl(spec2, colnames(data)))
  index3 = with(data, grepl(spec3, colnames(data)))
  index4 = with(data, grepl(spec4, colnames(data)))
  index5 = with(data, grepl(spec5, colnames(data)))
  ind = (index1 | index2 | index3 | index4 | index5) # To remove duplicated columns when there is for instance a comparaison between "m" and "A"
  #species <<- makeSpeciesVector(data)
  #geographies <<- makeGeographiesVector(data)
  #origins <<- makeOriginsVector(data)
  #specFac <<- selectSpeciesFac(ind) # I made sure the order of the species fac & geographies fac are correct
  #geoFac <<- selectGeographiesFac(ind)
  #oriFac <<- selectOriginsFac(ind)
  #geoNum <<- geoMap[as.character(geoFac)]
  #specieColors <<- colorMap[as.character(specFac)]
  #batchMatrix <<- makeBatchMatrix(data = data[ind])
  #ploidyFac <<- selectPloidiesFac(ind)
  return(data[ind])
  #return(ind)
}
#---- -
## spec1 & spec2 are the first character of the chosen species
## eg. fuchsii = "f" or "fP" or "A", one spec should be " " if only two comparaisons
## (Should be put in the variable dat)


## RUV (REMOVAL OF UNWANTED VARIANCE)

# ----+ function: makeRUVset
makeRUVset = function (dat){
  #plotRLE(as.matrix(dat),outline=FALSE,ylim=c(-4,4),col=specieColors, main = "Unwanted Variance not removed")
  uq = betweenLaneNormalization(as.matrix(dat), which = "full")
  set = newSeqExpressionSet(uq)
  #plotRLE(set,outline=FALSE,ylim=c(-4,4),col=specieColors, main = "Unwanted Variance not removed")
  
  #plotPCA(set, col = specieColors, k=2, cex = 1.2, main = "First two PCs'")#, labels = T)#pch = geoNum
  
  return(set)
}

plot_pca_set <- function (set){
    plotRLE(set,outline=FALSE,ylim=c(-4,4),col=specieColors, main = "Unwanted Variance not removed")
    plotPCA(set, col = specieColors, k=2, cex = 1.2, main = "First two PCs'")
}
#---- -
## Quick analysis without any normalisation
## (Should be put in variable set)


## TODO STUFF, NORMALISATION WITH REPLICATE SAMPLE
# ----+ function: makeRUVnormalisedSet
## Normalise with genes, the "empirical" ones found with findNonDEgenes
makeRUVnormalisedSet = function(set, empircalNonDEgenes, numFactors){
  normalisedSet = RUVg(as.matrix(counts(set)), empircalNonDEgenes, k=numFactors)
  return(normalisedSet)
}
#---- -
## (Should be put in variable normSet)

#Using replicate samples to remove variance
makeRUVrepNormalisedSet = function(set, empircalNonDEgenes, numFactors, replicates){
  normalisedSet = RUVs(as.matrix(counts(set)), cIdx=empircalNonDEgenes, k=numFactors, scIdx=replicates)
  #plotRLE(normalisedSet$normalizedCounts,outline=F, col=specieColors, cex = 0.5, main = "Unwanted Variance removed")
  #plotPCA(normalisedSet$normalizedCounts, k=2, col = specieColors, main = "First two PCs'")#pch=geoNum
  return(normalisedSet)
}