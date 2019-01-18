library("gplots")
library("RUVSeq")
library("edgeR")
library("HTSFilter")
library("RColorBrewer")
library("adegenet")
library("statmod")
library("cgwtools")                                        #---- -

#red <- transp("red",.4)
red = rgb(255, 0, 0, max = 255, alpha = 125, names = "red")
black <- rgb(0, 0, 0, max = 255, alpha = 125, names = "black")
## GLOBAL VARIABLES & FUNCTION FOR DATA ORGANISATION
#data = read.table("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/modelsContrasts/trimmed.counts.txt/trdata.txt", header = T, sep = " ")
data <- read.table("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/19.01.2017.counts.txt", header = TRUE, row.names = 1)
data <- read.table("/home/botanik/Documents/phd/research/siRNA/differential.expression/counts_sRNA20-22.txt", header = TRUE, row.names = 1)
data.orchis <- read.table("/data/phdData/orchis/counts/counts.corset.txt", header = TRUE, row.names = 1)
data.hylite <- read.table("/data/phdData/orchis/hylite/HyLiTE_analysis/HyLiTE_analysis.expression.txt", header = TRUE, row.names = 1)

load("/data/phdData/orchis/counts/count.featureCounts.RData")
data.featureCounts <- as.data.frame(fc_PE$counts)
colnames(data.featureCounts) <- c("tB1830_1","tB1830_2","tB1833_1","tB1833_2",
                            "fB1804","fB1855","fP1001","fP1707_1","fP1707_2",
                            "iA1586","iB1176","iB1870","iS1904","iS1908",
                            "mA1567_1","mA1567_2","mA1568","mA1573","mA1661","mA1775","mP1722_1","mP1722_2","mP1744","mP1748_1","mP1748_2","mS1757","mS1765_1","mS1765_2",
                            "tA1553","tA1641","tA1670","tB1798_1","tB1798_2","tB1805_1","tB1805_2","tB1812_1","tB1812_2","tS1901","tS1902","tS1920_1","tS1920_2")

combineRep = function(data){
  data = transform(data,  tB1830=tB1830_1+tB1830_2, tB1833=tB1833_1+tB1833_2, ## ebudensis
                   fP1707=fP1707_1+fP1707_2, ##fuchsii
                   mA1567=mA1567_1+mA1567_2,  mP1722= mP1722_1+mP1722_2, mP1748=mP1748_1+mP1748_2, mS1765=mS1765_1+mS1765_2, ##majalis
                   tB1798=tB1798_1+tB1798_2, tB1805=tB1805_1+tB1805_2, tB1812=tB1812_1+tB1812_2, tS1920=tS1920_1+tS1920_2) ##traunteineri
  data = data[, -grep("_", colnames(data))]
  return(data)
}

species <- c()
#dataC = dataC[,-grep("_", colnames(dataC))]

# ----+ functions: makeVectors
makeSpeciesVector = function(data){
  indm = with(data, grepl("m", colnames(data))); indm[indm==1] = 2
  indt = with(data, grepl("t", colnames(data))); indt[indt==1] = 3
  inde = with(data, grepl("e", colnames(data))); inde[inde==1] = 3
  indf = with(data, grepl("f", colnames(data))); indf[indf==1] = 5
  indi = with(data, grepl("i", colnames(data))); indi[indi==1] = 6 #modifie here if you want to make ebu => traun
  ind = indm+indt+inde+indf+indi; ind = sub(2, "majalis", ind); ind = sub(3, "traunsteineri", ind); ind = sub(5, "fuchsii", ind); ind = sub(6, "incarnata", ind)
  species <<- ind
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
  species <<- makeSpeciesVector(data)
  geographies <<- makeGeographiesVector(data)
  origins <<- makeOriginsVector(data)
  specFac <<- selectSpeciesFac(ind) # I made sure the order of the species fac & geographies fac are correct
  geoFac <<- selectGeographiesFac(ind)
  oriFac <<- selectOriginsFac(ind)
  geoNum <<- geoMap[as.character(geoFac)]
  specieColors <<- colorMap[as.character(specFac)]
  batchMatrix <<- makeBatchMatrix(data = data[ind])
  #ploidyFac <<- selectPloidiesFac(ind)
  return(data[ind])
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
  plotRLE(set,outline=FALSE,ylim=c(-4,4),col=specieColors, main = "Unwanted Variance not removed")
  
  plotPCA(set, col = specieColors, k=2, cex = 1.2, main = "First two PCs'")#, labels = T)#pch = geoNum
  
  return(set)
}
#---- -
## Quick analysis without any normalisation
## (Should be put in variable set)


## TODO STUFF, NORMALISATION WITH REPLICATE SAMPLE
# ----+ function: makeRUVnormalisedSet
## Normalise with genes, the "empirical" ones found with findNonDEgenes
makeRUVnormalisedSet = function(set, empircalNonDEgenes, numFactors){
  normalisedSet = RUVg(as.matrix(counts(set)), empircalNonDEgenes, k=numFactors)
  plotRLE(normalisedSet$normalizedCounts,outline=F, col=specieColors, cex = 0.5, main = "Unwanted Variance removed")
  plotPCA(normalisedSet$normalizedCounts, k=2, pch=geoNum, col = specieColors, main = "First two PCs'")
  return(normalisedSet)
}
#---- -
## (Should be put in variable normSet)

#Using replicate samples to remove variance
makeRUVrepNormalisedSet = function(set, empircalNonDEgenes, numFactors, replicates){
  normalisedSet = RUVs(as.matrix(counts(set)), cIdx=empircalNonDEgenes, k=numFactors, scIdx=replicates)
  plotRLE(normalisedSet$normalizedCounts,outline=F, col=specieColors, cex = 0.5, main = "Unwanted Variance removed")
  plotPCA(normalisedSet$normalizedCounts, k=2, col = specieColors, main = "First two PCs'")#pch=geoNum
  return(normalisedSet)
}
data.featureCounts[colnames(data.featureCounts) == "mA1567_1" | colnames(data.featureCounts) == "mA1567_2"] <- NULL
data.featureCounts[colnames(data.featureCounts) == "mP1748_1" | colnames(data.featureCounts) == "mP1748_2"] <- NULL

chloro.mito <- read.table("/data/phdData/orchis/counts/de.featureCounts/HEB/chloroplastic.and.mitochondrial.genes.txt")
#REMOVE TRANSCRIPTS WITH NO ANNOTATION FULL=.
dat <- selectSpecies(data.featureCounts, "P1.t", "P2.t", " "," "," ")

myvars <- rownames(dat) %in% chloro.mito$V1
dat.no.chloro <- dat[!myvars,]
#check filtering step
keep <- rowSums(cpm(dat)>0.08) >= 10
summary(keep)
dat.filtered <- dat[keep,]
set <- makeRUVset(dat = dat.filtered)
normSet <- makeRUVrepNormalisedSet(set, rownames(counts(set)), 3, batchMatrix)
design <- model.matrix(~specFac+normSet$W-1)
design <- model.matrix(~geoFac+normSet$W-1)
colnames(design) <- c("majalis", "traunsteineri", "nW1", "nW2")#, "nW3")
y <- DGEList(counts = counts(set), group = specFac) # my stuff: counts(set), group = group) # Don't forget, model is taking into account the RUV normalised matrix of counts
#get an idea of how to filter
cpm(5, mean(y$samples$lib.size))
y <- calcNormFactors(y)
y <- estimateDisp(y, design = design, robust = T)
plotBCV(y)
fit <- glmQLFit(y, design, robust = T)
hist(fit$coefficients[,1], breaks = 100)
summary(fit$df.prior)
plotQLDisp(fit)
con <- makeContrasts(traunsteineri - majalis, levels = design)
Group <- factor(c(rep("maj.fuch", 10), rep("trau.fuch", 17), rep("maj.inc", 10), rep("trau.inc", 17)))
design <- model.matrix(~Group+normSet$W-1)
colnames(design) <- c("maj.fuch", "maj.inc", "trau.fuch", "trau.inc", "nW1", "nW2", "nW3")
my.contrasts <- makeContrasts(majHEB = maj.inc - maj.fuch,
                              trauHEB = trau.inc - trau.fuch,
                              fuchAE = trau.fuch - maj.fuch,
                              incAE = trau.inc - maj.inc, levels = design)
resEdge <- glmQLFTest(fit, contrast = con)
#res <- glmQLFTest(fit, coef = 2)
de = decideTestsDGE(res, adjust.method="fdr", p.value = 0.05)
#de.genes = rownames(res)[as.logical(de)]
#plotSmear(res, de.tags = de.genes, cex=0.5, col=c("blue","red"),main="MA-plot (p-value <= 0.05)")
is.de <- decideTestsDGE(resEdge)
plotMD(resEdge, status=is.de, values=c(1,0,-1), col=c("blue","purple", "red"), legend="topright", pch=19, cex=c(0.5,0.1))
resEdgeRtopTags <- topTags(resEdge, n=nrow(set))
FDRres = topTags(resEdgeRtopTags, n=nrow(dat))$table
write.table(FDRres[FDRres$FDR <= 0.05,], "/data/phdData/orchis/counts/de.featureCounts/both.de.results.txt", sep="\t", quote = F)



#TEST
par(mfrow = c(1,2))
retained = c()
stripround = c()
de = c()
cpmval = seq(0,2, by = 0.25)
dat = selectSpecies(data, "f", "i", "m","t"," ")
for(j in c(0,2,5))
{
  for(i in cpmval)
    #for(k in c(1,3,5,7))
  {
    keep <- rowSums(cpm(dat)>i) >= j
    retained = c(retained, sum(keep))
    datS <- dat[keep,]
    #dat = selectSpecies(data[rowSums(cpm(data[,c(5,6,7,8,9)])>i) >= j | rowSums(cpm(data[,c(10,11,12,13,14)])>i) >= j,], "f", "i", "m","t"," ")
    set <- makeRUVset(dat = datS)
    normSet <- makeRUVrepNormalisedSet(set, rownames(counts(set)), 3, batchMatrix)
    y <- DGEList(counts = counts(set), group = specFac) # my stuff: counts(set), group = group) # Don't forget, model is taking into account the RUV normalised matrix of counts
    design <- model.matrix(~specFac+normSet$W-1)
    colnames(design) <- c("fuchsii", "incarnata", "majalis", "traunsteineri", "nW1", "nW2", "nW3")
    y <- DGEList(counts = counts(set), group = specFac) # my stuff: counts(set), group = group) # Don't forget, model is taking into account the RUV normalised matrix of counts
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design = design, robust = T)
    fit <- glmQLFit(y, design, robust = T)
    #hist(fit$coefficients[,1], breaks = 100)
    con <- makeContrasts(traunsteineri - majalis, levels = design)
    res <- glmQLFTest(fit, contrast = con)
    resEdgeRtopTags <- topTags(res, n=nrow(set))
    de = c(de, nrow(resEdgeRtopTags$table[resEdgeRtopTags$table$FDR <= 0.05,]))
    #stripround = c(stripround, length(rownames(fit$coefficients[fit$coefficients[,1] < -15 | fit$coefficients[,2] < -15 | fit$coefficients[,3] < -15 | fit$coefficients[,4] < -15,])))
  }
  plot(cpmval, retained)
  #plot(cpmval, stripround)
  plot(cpmval, de)
  retained = c()
  #stripround = c()
  de = c()
}

#LOOP FOR COMPARISONS
dat <- selectSpecies(dat, "f", "i", "m","t"," ")
dat[colnames(dat) == "mA1567_1" | colnames(dat) == "mA1567_2" | colnames(dat) == "mP1748_1" | colnames(dat) == "mP1748_2"] <- NULL
y <- DGEList(counts = dat, group = specFac)
cpm(5, mean(y$samples$lib.size))
keep <- rowSums(cpm(dat)>0.1) >= 5
dat.filtered <- dat[keep,]
set <- makeRUVset(dat = dat.filtered)
normSet <- makeRUVrepNormalisedSet(set, rownames(counts(set)), 3, batchMatrix)
design <- model.matrix(~specFac+normSet$W-1)
colnames(design) <- c("fuchsii", "incarnata", "majalis", "traunsteineri", "nW1", "nW2", "nW3")
y <- DGEList(counts = counts(set), group = specFac) 
y <- calcNormFactors(y)
y <- estimateDisp(y, design = design, robust = T)
plotBCV(y)
fit <- glmQLFit(y, design, robust = T)
hist(fit$coefficients[,1], breaks = 100)
summary(fit$df.prior)
plotQLDisp(fit)

contrast = list(c(-1,1,0,0,0,0,0), c(-1,0,1,0,0,0,0), c(-1,0,0,1,0,0,0), c(0,-1,1,0,0,0,0),c(0,-1,0,1,0,0,0),c(0,0,-1,1,0,0,0))
for(i in seq(1,6,by=1))
{
  res <- glmQLFTest(fit, contrast = as.vector(contrast[[i]]))
  resEdgeRtopTags <- topTags(res, n=nrow(set))
  FDRres = topTags(resEdgeRtopTags, n=nrow(dat))$table
  write.table(FDRres[FDRres$FDR <= 0.05,],  paste("/data/phdData/orchis/counts/de.featureCounts/contrasts/contrasts", contrast[i], "txt", sep="."), sep="\t", quote = F)
  save(resEdgeRtopTags, file = paste("/data/phdData/orchis/counts/de.featureCounts/contrasts/contrasts", contrast[i], "RData", sep="."))
}

########################### PLOT NICE PCA ############################################
load("/home/botanik/Documents/phd/research/environmental_study/deEnviroCorrelations/2018/normSet.RData")
fac <- factor(c(rep("traunsteineri",4), rep("majalis",10), rep("traunsteineri",3), rep("traunsteineri",6), rep("traunsteineri",4)))
fac1 <- factor(c(rep("traunsteineri Britain",4), rep("majalis",10), rep("traunsteineri continent",3), rep("traunsteineri Britain",6), rep("traunsteineri continent",4)))
pca <- dudi.pca(t(normSet$normalizedCounts), center=T, scale=T)
s.class(pca$li, fac=fac1,
        col=c("red", "turquoise", "lightblue"),
        axesel=F, cstar=0, cpoint=1, xax = 1, yax = 2, clabel=0.6, add.plot = T,pch=1)

s.class(pca$li, fac=fac,
        col=c("red", "blue"),
        axesel=T, cstar=1, cpoint=1, xax = 1, yax = 2, clabel=0.9)
