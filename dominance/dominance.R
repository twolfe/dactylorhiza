library("RColorBrewer")
library("edgeR")
library("corrplot")
library("lmPerm")
library("ggplot2")
library("d3radarR")
library("ggradar")
library("doBy")
library("dplyr")
library("plotly")
library("plyr")
library("fmsb")

#install_github("jerryzhujian9/ezmisc")
suppressPackageStartupMessages(library(dplyr))
library("scales")

add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
                     function(x) 
                       rgb(x[1], x[2], x[3], alpha=alpha))  
}
map.category.color = c("TRY1"="black", "TRY2"="black", "TRY3"="black", "TRY4"="black", "TRY5"="black", "TRY6"="black", "TRY7"="black", "TRY8"="black", "TRY9"="black", "TRY10"="black", "TRY11"="black", "TRY12"="black", "TRY13"="black", "TRY14"="black", 
                       "I"="magenta", "II"="purple", "III"="darkgreen", "IV"="orange", "V"="darkgreen", "VI"="darkgreen", "VII"="steelblue", "VIII"="steelblue", "IX"="orange", "X"="darkgreen", "XI"="purple", "XII"="magenta", "XIII"="black")


symdiff <- function(x, y) { setdiff( union(x, y), intersect(x, y))}

red = add.alpha("red", 0.7)
blue = add.alpha("blue",0.7)
darkgreen = add.alpha("darkgreen",0.5)
green = add.alpha("green",0.5)
purple = add.alpha("purple",0.1)
magenta = add.alpha("magenta",0.5)
steelblue = add.alpha("steelblue",0.6)
yellow = add.alpha("yellow",0.5)
darkblue = add.alpha("darkblue",0.3)
black = add.alpha("black",0.01)
orange = add.alpha("orange",0.1)
grey = add.alpha("grey",0.01)

#DIPLOIDS
load("/data/phdData/orchis/counts/de.featureCounts/contrasts/contrasts.c(-1, 1, 0, 0, 0, 0, 0).RData")
#load("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/DEanalysisObjects/fuchsii.incarnata.RData")
setFI = set
DEfi= resEdgeRtopTags
#DEmtNoAlps005= DEfi[DEfi$table$FDR <= 0.05,]
GLMfi = resEdgeRglm
# (<0) => Up regulated in incarnata, (>0) => Up regulated in traunsteineri

#MAJALIS
load("/data/phdData/orchis/counts/de.featureCounts/contrasts/contrasts.c(-1, 0, 1, 0, 0, 0, 0).RData")
#load("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/DEanalysisObjects/fuchsii.majalis.RData")
setFM = set
DEfm= resEdgeRtopTags
#DEmtNoAlps005= DEfi[DEfi$table$FDR <= 0.05,]
GLMfm = resEdgeRglm
# (<0) => Up regulated in incarnata, (>0) => Up regulated in traunsteineri

load("/data/phdData/orchis/counts/de.featureCounts/contrasts/contrasts.c(0, -1, 1, 0, 0, 0, 0).RData")
#load("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/DEanalysisObjects/incarnata.majalis.RData")
setIM = set
DEim= resEdgeRtopTags
#DEmtNoAlps005= DEfi[DEfi$table$FDR <= 0.05,]
GLMim = resEdgeRglm
# (<0) => Up regulated in incarnata, (>0) => Up regulated in traunsteineri

#TRAUNSTEINERI
load("/data/phdData/orchis/counts/de.featureCounts/contrasts/contrasts.c(-1, 0, 0, 1, 0, 0, 0).RData")

#load("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/DEanalysisObjects/fuchsii.traunsteineri.RData")
setFT = set
DEft= resEdgeRtopTags
#DEmtNoAlps005= DEfi[DEfi$table$FDR <= 0.05,]
GLMft = resEdgeRglm
# (<0) => Up regulated in incarnata, (>0) => Up regulated in traunsteineri

load("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/DEanalysisObjects/fuchsii.traunsteineri.britain.RData")
setFTb = set
DEftB= resEdgeRtopTags
#DEmtNoAlps005= DEfi[DEfi$table$FDR <= 0.05,]
GLMftB = resEdgeRglm

load("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/DEanalysisObjects/fuchsii.traunsteineri.continent.RData")
setFTc = set
DEftC= resEdgeRtopTags
#DEmtNoAlps005= DEfi[DEfi$table$FDR <= 0.05,]
GLMftC = resEdgeRglm

load("/data/phdData/orchis/counts/de.featureCounts/contrasts/contrasts.c(0, -1, 0, 1, 0, 0, 0).RData")
#load("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/DEanalysisObjects/incarnata.traunsteineri.RData")
setIT = set
DEit= resEdgeRtopTags
#DEmtNoAlps005= DEfi[DEfi$table$FDR <= 0.05,]
GLMit = resEdgeRglm
# (<0) => Up regulated in incarnata, (>0) => Up regulated in traunsteineri

load("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/DEanalysisObjects/incarnata.traunsteineri.britain.RData")
setITb = set
DEitB= resEdgeRtopTags
#DEmtNoAlps005= DEfi[DEfi$table$FDR <= 0.05,]
GLMitB = resEdgeRglm

load("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/DEanalysisObjects/incarnata.traunsteineri.continent.RData")
setITc = set
DEitC= resEdgeRtopTags
#DEmtNoAlps005= DEfi[DEfi$table$FDR <= 0.05,]
GLMitC = resEdgeRglm

#POLYPLOIDS
load("/data/phdData/orchis/counts/de.featureCounts/contrasts/contrasts.c(0, 0, -1, 1, 0, 0, 0).RData")
#load("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/DEanalysisObjects/majalis.traunsteineri.RData")
setMT = set
DEmt = resEdgeRtopTags
#DEmtNoAlps005= DEfi[DEfi$table$FDR <= 0.05,]
GLMmt = resEdgeRglm
# (<0) => Up regulated in incarnata, (>0) => Up regulated in traunsteineri

load("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/DEanalysisObjects/majalis.traunsteineri.britain.RData")
setMTb = set
DEmtB = resEdgeRtopTags
#DEmtNoAlps005= DEfi[DEfi$table$FDR <= 0.05,]
GLMmtB = resEdgeRglm

load("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/DEanalysisObjects/majalis.traunsteineri.continent.RData")
setMTc = set
DEmtC = resEdgeRtopTags
#DEmtNoAlps005= DEfi[DEfi$table$FDR <= 0.05,]
GLMmtC = resEdgeRglm

load("/home/botanik/Documents/phd/research/differential.expression/2016.11/withExtra2RNAlibraries/without.mP1748.mA1567/DEanalysisObjects/tB.tC.continent.RData")
setTbTc = set
DEtBtC = resEdgeRtopTags
#DEmtNoAlps005= DEfi[DEfi$table$FDR <= 0.05,]
GLMtBtC = resEdgeRglm

dominanceFun = function(topTags1, topTags2, topTags3, FDRval, spec){
  tab1 <<- topTags1
  tab2 <<- topTags2
  tab3 <<- topTags3
  DE1 <<- rownames(topTags1[topTags1$table$FDR <= FDRval,])
  DE2 <<- rownames(topTags2[topTags2$table$FDR <= FDRval,])
  DE3 <<- rownames(topTags3[topTags3$table$FDR <= FDRval,])
  nonDE1 <<- rownames(topTags1[topTags1$table$FDR > FDRval,])
  nonDE2 <<- rownames(topTags2[topTags2$table$FDR > FDRval,])
  nonDE3 <<- rownames(topTags3[topTags3$table$FDR > FDRval,])

  UpFu1 <<- rownames(tab1$table[tab1$table$logFC <= 0.0 & tab1$table$FDR <= FDRval,])#++
  print(paste("UpFu1 size", as.character(length(UpFu1)), sep = " = "))
  
  UpIn1 <<- rownames(tab1$table[tab1$table$logFC > 0.0 & tab1$table$FDR <= FDRval,])#**
  print(paste("UpIn1 size", as.character(length(UpIn1)), sep = " = "))
  
  UpFuPoly <<- rownames(tab2$table[tab2$table$logFC <= 0.0 & tab2$table$FDR <= FDRval,])#++
  print(paste("UpFuPoly size", as.character(length(UpFuPoly)), sep = " = "))
  
  UpPolyFu <<- rownames(tab2$table[tab2$table$logFC > 0.0 & tab2$table$FDR <= FDRval,])#**
  print(paste("UpPolyFu size", as.character(length(UpPolyFu)), sep = " = "))
  
  UpInPoly <<- rownames(tab3$table[tab3$table$logFC <= 0.0 & tab3$table$FDR <= FDRval,])#**
  print(paste("UpInPoly size", as.character(length(UpInPoly)), sep = " = "))
  
  UpPolyIn <<- rownames(tab3$table[tab3$table$logFC > 0.0 & tab3$table$FDR <= FDRval,])#++
  print(paste("UpPolyIn size", as.character(length(UpPolyIn)), sep = " = "))
  
  
  
  #DE1, DE2...
  XIII = Reduce(intersect, list(nonDE1, nonDE2, nonDE3))#V
  #assign(paste(as.character(spec), "XIII", sep = "."), XIII, envir = .GlobalEnv)
  #plot(tab2$table[XIII,]$logFC, tab3$table[XIII,]$logFC, pch=20, cex=0.4, col=black, xlim=c(-15,15),ylim=c(-15,15), 
  #    xlab = paste(spec, "-fuchsii logFC", sep=""), ylab = paste(spec, "-incarnata logFC", sep=""))
  #p <- ggplot(d,aes(x+0.05,y+0.05))+geom_tile(aes(fill=enrichment)) + scale_fill_gradient(low="white", high="steelblue", limits=c(0,120))
  
  TRY1 = Reduce(intersect, list(nonDE1, nonDE2, UpInPoly))#XIII
  points(tab2$table[TRY1,]$logFC, tab3$table[TRY1,]$logFC, pch=20, cex=0.4, col=black)
  
  TRY2 = Reduce(intersect, list(nonDE1, nonDE2, UpPolyIn))#XIII
  points(tab2$table[TRY2,]$logFC, tab3$table[TRY2,]$logFC, pch=20, cex=0.4, col=black)
  
  #DE1, UpFuPoly...
  TRY3 = Reduce(intersect, list(nonDE1, UpFuPoly, nonDE3))#XIII
  points(tab2$table[TRY3,]$logFC, tab3$table[TRY3,]$logFC, pch=20, cex=0.4, col=black)
  
  VII = Reduce(intersect, list(nonDE1, UpFuPoly, UpInPoly))
  assign(paste(as.character(spec), "VII", sep = "."), VII, envir = .GlobalEnv)
  points(tab2$table[VII,]$logFC, tab3$table[VII,]$logFC, pch=20, cex=0.4, col=steelblue)
  
  TRY4 = Reduce(intersect, list(nonDE1, UpFuPoly, UpPolyIn))#XIII
  points(tab2$table[TRY4,]$logFC, tab3$table[TRY4,]$logFC, pch=20, cex=0.4, col=black)
  
  #DE1, UpPolyFu...
  TRY5 = Reduce(intersect, list(nonDE1, UpPolyFu, nonDE3))#XIII
  points(tab2$table[TRY5,]$logFC, tab3$table[TRY5,]$logFC, pch=20, cex=0.4, col=black)
  
  TRY6 = Reduce(intersect, list(nonDE1, UpPolyFu, UpInPoly))#XIII
  points(tab2$table[TRY6,]$logFC, tab3$table[TRY6,]$logFC, pch=20, cex=0.4, col=black)
  
  VIII = Reduce(intersect, list(nonDE1, UpPolyFu, UpPolyIn))#V
  assign(paste(as.character(spec), "VIII", sep = "."), VIII, envir = .GlobalEnv)
  points(tab2$table[VIII,]$logFC, tab3$table[VIII,]$logFC, pch=20, cex=0.4, col=steelblue)
  
  #XIII <- c(XIII, TRY1, TRY2, TRY3, TRY4, TRY5, TRY6)
  assign(paste(as.character(spec), "XIII", sep = "."), XIII, envir = .GlobalEnv)
  plot(tab2$table[XIII,]$logFC, tab3$table[XIII,]$logFC, pch=20, cex=0.4, col=black, xlim=c(-5,5),ylim=c(-5,5), 
       xlab = paste(spec, "-fuchsii logFC", sep=""), ylab = paste(spec, "-incarnata logFC", sep=""))
  print(XIII)
  #UpFu1,DE2...
  TRY7 = Reduce(intersect, list(UpFu1, nonDE2, nonDE3))#I
  points(tab2$table[TRY7,]$logFC, tab3$table[TRY7,]$logFC, pch=20, cex=0.4, col=black)
  
  TRY8 = Reduce(intersect, list(UpFu1, nonDE2, UpInPoly))#TBA :D 
  points(tab2$table[TRY8,]$logFC, tab3$table[TRY8,]$logFC, pch=20, cex=0.4, col=black)
  XIII <- c(XIII, TRY1, TRY2, TRY3, TRY4, TRY5, TRY6, TRY8)
  
  II = Reduce(intersect, list(UpFu1, nonDE2, UpPolyIn))
  assign(paste(as.character(spec), "II", sep = "."), II, envir = .GlobalEnv)
  points(tab2$table[II,]$logFC, tab3$table[II,]$logFC, pch=20, cex=0.4, col=purple)
  
  #UpFu1, UpFuPoly...
  IX = Reduce(intersect, list(UpFu1, UpFuPoly, nonDE3))
  assign(paste(as.character(spec), "IX", sep = "."), IX, envir = .GlobalEnv)
  points(tab2$table[IX,]$logFC, tab3$table[IX,]$logFC, pch=20, cex=0.4, col=orange)
  
  III = Reduce(intersect, list(UpFu1, UpFuPoly, UpInPoly))
  assign(paste(as.character(spec), "III", sep = "."), III, envir = .GlobalEnv)
  points(tab2$table[III,]$logFC, tab3$table[III,]$logFC, pch=20, cex=0.4, col=darkgreen)
  
  
  I = Reduce(intersect, list(UpFu1, UpFuPoly, UpPolyIn))#
  #assign(paste(as.character(spec), "I", sep = "."), I, envir = .GlobalEnv) 
  #points(tab2$table[I,]$logFC, tab3$table[I,]$logFC, pch=20, cex=0.4, col=magenta)
  I <- c(I, TRY7)
  assign(paste(as.character(spec), "I", sep = "."), I, envir = .GlobalEnv)
  points(tab2$table[I,]$logFC, tab3$table[I,]$logFC, pch=20, cex=0.4, col=magenta)
  
  #UpFu1, UpPolyFu...
  TRY9 = Reduce(intersect, list(UpFu1, UpPolyFu, nonDE3))
  points(tab2$table[TRY9,]$logFC, tab3$table[TRY9,]$logFC, pch=20, cex=0.4, col=black)
  
  TRY10 = Reduce(intersect, list(UpFu1, UpPolyFu, UpInPoly))
  points(tab2$table[TRY10,]$logFC, tab3$table[TRY10,]$logFC, pch=20, cex=0.4, col=black)
  
  V = Reduce(intersect, list(UpFu1, UpPolyFu, UpPolyIn))#
  assign(paste(as.character(spec), "V", sep = "."), V, envir = .GlobalEnv)
  points(tab2$table[V,]$logFC, tab3$table[V,]$logFC, pch=20, cex=0.4, col=darkgreen)
  
  #UpIn1, DE2...
  TRY11 = Reduce(intersect, list(UpIn1, nonDE2, nonDE3))#XII
  points(tab2$table[TRY11,]$logFC, tab3$table[TRY11,]$logFC, pch=20, cex=0.4, col=black)
  
  XI = Reduce(intersect, list(UpIn1, nonDE2, UpInPoly))#
  assign(paste(as.character(spec), "XI", sep = "."), XI, envir = .GlobalEnv)
  points(tab2$table[XI,]$logFC, tab3$table[XI,]$logFC, pch=20, cex=0.4, col=purple)
  
  TRY12 = Reduce(intersect, list(UpIn1, nonDE2, UpPolyIn))
  points(tab2$table[TRY12,]$logFC, tab3$table[TRY12,]$logFC, pch=20, cex=0.4, col=black)
  
  #UpIn1, UpFuPoly...
  TRY13 = Reduce(intersect, list(UpIn1, UpFuPoly, nonDE3))
  points(tab2$table[TRY13,]$logFC, tab3$table[TRY13,]$logFC, pch=20, cex=0.4, col=black)
  
  X = Reduce(intersect, list(UpIn1, UpFuPoly, UpInPoly))#
  assign(paste(as.character(spec), "X", sep = "."), X, envir = .GlobalEnv)
  points(tab2$table[X,]$logFC, tab3$table[X,]$logFC, pch=20, cex=0.4, col=darkgreen)
  
  TRY14 = Reduce(intersect, list(UpIn1, UpFuPoly, UpPolyIn))
  points(tab2$table[TRY14,]$logFC, tab3$table[TRY14,]$logFC, pch=20, cex=0.4, col=black)
  
  
  #UpIn1, UpPolyFu...
  IV = Reduce(intersect, list(UpIn1, UpPolyFu, nonDE3))#
  assign(paste(as.character(spec), "IV", sep = "."), IV, envir = .GlobalEnv)
  #points(tab2$table[IV,]$logFC, tab3$table[IV,]$logFC, pch=20, cex=0.4, col=orange)
  
  XII = Reduce(intersect, list(UpIn1, UpPolyFu, UpInPoly))#
  XII <- c(XII, TRY11)
  assign(paste(as.character(spec), "XII", sep = "."), XII, envir = .GlobalEnv)
  points(tab2$table[XII,]$logFC, tab3$table[XII,]$logFC, pch=20, cex=0.4, col=magenta)
  
  VI = Reduce(intersect, list(UpIn1, UpPolyFu, UpPolyIn))#
  assign(paste(as.character(spec), "VI", sep = "."), VI, envir = .GlobalEnv)
  points(tab2$table[VI,]$logFC, tab3$table[VI,]$logFC, pch=20, cex=0.4, col=darkgreen)
  
  XIII <- c(XIII, TRY13)
  # EXTRA MISSING POINTS TRY

  allnames <<- c(TRY1, TRY2, TRY3, TRY4, TRY5, TRY6, TRY7,TRY8, TRY9, TRY10, TRY11, TRY12, TRY13, TRY14, I, II, III, IV, V, VI, VII, VIII, IX, X, XI, XII, XIII)
  namnam <<- c(I, II, III, IV, V, VI, VII, VIII, IX, X, XI, XII, XIII)
  notIn.names <<- symdiff(rownames(DEfi), namnam)
  #print(length(notIn.names))
  print(length(interaction(namnam)))
  points(tab2$table[notIn.names,]$logFC, tab3$table[notIn.names,]$logFC, pch=20, cex=4, col=steelblue)
  
  #PLOTING NICE
  points(tab2$table[II,]$logFC, tab3$table[II,]$logFC, pch=20, cex=0.4, col=orange)
  points(tab2$table[IV,]$logFC, tab3$table[IV,]$logFC, pch=20, cex=0.4, col=purple)
  points(tab2$table[IX,]$logFC, tab3$table[IX,]$logFC, pch=20, cex=0.4, col=purple)
  points(tab2$table[XI,]$logFC, tab3$table[XI,]$logFC, pch=20, cex=0.4, col=orange)
  
  abline(h=0, v=0)
  
  dataframe <- data.frame("geneID" = namnam, "logFC.fuchsii" = tab2$table[namnam,]$logFC, "logFC.inacarnata" = tab3$table[namnam,]$logFC,
                          "category" = c( #rep("TRY1", length(TRY1)), rep("TRY2", length(TRY2)), rep("TRY3", length(TRY3)), rep("TRY4", length(TRY4)),
                                          #rep("TRY5", length(TRY5)), rep("TRY6", length(TRY6)), rep("TRY7", length(TRY7)), rep("TRY8", length(TRY8)),
                                          #rep("TRY9", length(TRY9)), rep("TRY10", length(TRY10)), rep("TRY11", length(TRY11)), rep("TRY12", length(TRY12)),
                                          #rep("TRY13", length(TRY13)), rep("TRY14", length(TRY14)),
                                          rep("I", length(I)), rep("II", length(II)), rep("III", length(III)), 
                                          rep("IV", length(IV)), rep("V", length(V)), rep("VI", length(VI)), 
                                          rep("VII", length(VII)), rep("VIII", length(VIII)), rep("IX", length(IX)),
                                          rep("X", length(X)), rep("XI", length(XI)), rep("XII", length(XII)),
                                          rep("XIII", length(XIII)) )
                          )
  
  dataframe <- data.frame(dataframe, "color" = map.category.color[as.character(dataframe$category)])
  
  assign(paste(as.character(spec), "dataframe", sep = "."), dataframe, envir = .GlobalEnv)
}

# write.table(file="traunsteineri.up.transgressive.txt", c(traunsteineri.V,traunsteineri.VI,traunsteineri.VIII))
# write.table(file="traunsteineri.up.transgressive.txt", c(traunsteineri.V,traunsteineri.VI,traunsteineri.VIII),quote=F,row.names=F,col.names=F)
# write.table(file="traunsteineri.down.transgressive.txt", c(traunsteineri.III,traunsteineri.VII,traunsteineri.X),quote=F,row.names=F,col.names=F)

#+RIP
################## POLAR COORDINATES ####################################
getSteps <- function(x,y) {
  d <- diff(complex(real = x, imaginary = y))
  data.frame(size = Mod(d), 
             angle = c(NA, diff(Arg(d)) %% (2*pi)) * 360/(2*pi))
}
####### Category comparaison / permutation / correlations / variance using technical replicates ############################
m = aovp(getSteps(tab2.maj$table[symdiff(allnames.maj,allnames.trau),]$logFC,tab3.maj$table[symdiff(allnames.maj,allnames.trau),]$logFC)$angle~getSteps(tab2.trau$table[symdiff(allnames.maj,allnames.trau),]$logFC,tab3.trau$table[symdiff(allnames.maj,allnames.trau),]$logFC)$angle, perm = "Exact")
summary.aovp(m)
plot(m)
plot(getSteps(tab2.maj$table[intersect(allnames.maj,allnames.trau),]$logFC,tab3.maj$table[intersect(allnames.maj,allnames.trau),]$logFC)$angle, getSteps(tab2.trau$table[intersect(allnames.maj,allnames.trau),]$logFC,tab3.trau$table[intersect(allnames.maj,allnames.trau),]$logFC)$angle, pch=20, cex=0.6, col=red)
points(getSteps(tab2.maj$table[symdiff(allnames.maj,allnames.trau),]$logFC,tab3.maj$table[symdiff(allnames.maj,allnames.trau),]$logFC)$angle, getSteps(tab2.trau$table[symdiff(allnames.maj,allnames.trau),]$logFC,tab3.trau$table[symdiff(allnames.maj,allnames.trau),]$logFC)$angle, pch=20, cex=0.6, col=green)

#########################################################
#            PLOTTING GGPLOT                            #
#########################################################
dominanceFun(DEfi,DEfm,DEim,0.05,"majalis")
dominanceFun(DEfi,DEft,DEit,0.05,"traunsteineri")
p <- ggplot(traunsteineri.dataframe[traunsteineri.dataframe$color != "black",], aes(x = logFC.fuchsii, y = logFC.inacarnata)) + geom_point(colour = traunsteineri.dataframe[traunsteineri.dataframe$color != "black",]$color, alpha=1/10)
p
p <- ggplot(majalis.dataframe[majalis.dataframe$color != "black",], aes(x = logFC.fuchsii, y = logFC.inacarnata)) + geom_point(colour = majalis.dataframe[majalis.dataframe$color != "black",]$color, alpha=1/10)
p

# RADAR PLOT
colors_border=c( rgb(1.,0.,0.,0.9), rgb(0.,0.,1.,0.9) )
colors_in=c( rgb(1.,0.,0.,0.4), rgb(0.,0.,1.,0.4) )

x = count(majalis.dataframe, 'category')
x["percentage"] <- x$freq/sum(x$freq)*100

y = count(traunsteineri.dataframe, 'category')
y["percentage"] <- y$freq/sum(y$freq)*100

x.good <- x[x$category %in% c("I","II","III","IV","IX","V","VI","VII","VIII","X","XI","XII"),]
y.good <- y[y$category %in% c("I","II","III","IV","IX","V","VI","VII","VIII","X","XI","XII"),]
# http://www.r-graph-gallery.com/143-spider-chart-with-saveral-individuals/

radar.dataframe <- data.frame("majalis" = x.good$percentage, "traunsteineri" = y.good$percentage)
rownames(radar.dataframe) <- c("I","II","III","IV","IX","V","VI","VII","VIII","X","XI","XII")
radar.dataframe <- as.data.frame(rbind(rep(8,12) , rep(0,12) , t(radar.dataframe)))

radar.dataframe.large <- radar.dataframe[c("IV","IX","VII","II","XI","I","XII")]

radar.dataframe.small <- radar.dataframe[c("III","V","VI","X","VIII")]
radar.dataframe.small <- radar.dataframe.small[-c(1,2),]
radar.dataframe.small <- as.data.frame(rbind(rep(0.13,5) , rep(0,5) , radar.dataframe.small))

radarchart( radar.dataframe.large[,sample(ncol(radar.dataframe.large))]  , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=1 , plty=1, pty = 20,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(0,8,2), cglwd=0.8,
            #custom labels
            vlcex=0.8 
)

radarchart( radar.dataframe.small[,sample(ncol(radar.dataframe.small))]  , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=1 , plty=1, pty = 20,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(0,0.13,0.03), cglwd=0.8,
            #custom labels
            vlcex=0.8 
)

#HEATMAP
majalis.dataframe.ordered <- majalis.dataframe[order(majalis.dataframe$geneID),]
traunsteineri.dataframe.ordered <- traunsteineri.dataframe[order(traunsteineri.dataframe$geneID),]
#ADD truth column to know if transcript is DE between the polyploids
majalis.dataframe.ordered <- data.frame(majalis.dataframe.ordered, "polyploidDE" = majalis.dataframe.ordered$geneID %in% rownames(DEmt$table[DEmt$table$FDR<=0.05,]))
traunsteineri.dataframe.ordered <- data.frame(traunsteineri.dataframe.ordered, "polyploidDE" = traunsteineri.dataframe.ordered$geneID %in% rownames(DEmt$table[DEmt$table$FDR<=0.05,]))

combined.ordered.categories <- data.frame("majalis" = majalis.dataframe.ordered$category, "traunsteineri" = traunsteineri.dataframe.ordered$category)
combined.ordered.polyploidDE <- data.frame("majalis" = majalis.dataframe.ordered$polyploidDE, "traunsteineri" = traunsteineri.dataframe.ordered$polyploidDE)

frequency.table <- table(combined.ordered.categories$majalis, combined.ordered.categories$traunsteineri)
frequency.table.polyploidDE <- table(combined.ordered.categories$majalis[combined.ordered.polyploidDE == "TRUE"], combined.ordered.categories$traunsteineri[combined.ordered.polyploidDE == "TRUE"])

frequency.dataframe <- as.data.frame(frequency.table[c("II","XI","IV","IX","V","VI","VIII","III","X","VII","I","XII","XIII"),c("II","XI","IV","IX","V","VI","VIII","III","X","VII","I","XII","XIII")])
frequency.dataframe.polyploidDE <- as.data.frame(frequency.table.polyploidDE[c("II","XI","IV","IX","V","VI","VIII","III","X","VII","I","XII","XIII"),c("II","XI","IV","IX","V","VI","VIII","III","X","VII","I","XII","XIII")])

#frequency.dataframe <- as.data.frame(frequency.table[c("II","XI","IV","V","VI","VIII","III","X","VII","I","XII","XIII"),c("II","XI","IV","V","VI","VIII","III","X","VII","I","XII","XIII")])#without cat "IX"
#frequency.dataframe.polyploidDE <- as.data.frame(frequency.table.polyploidDE[c("II","XI","IV","V","VI","VIII","III","X","VII","I","XII","XIII"),c("II","XI","IV","V","VI","VIII","III","X","VII","I","XII","XIII")])#without cat "IX"

frequency.dataframe.combined <- data.frame(frequency.dataframe, "Freq.polyploidDE" = frequency.dataframe.polyploidDE$Freq, 
                                           "percentage" = frequency.dataframe.polyploidDE$Freq/64150)

# You can also change the colorscale using the "colorscale parameter", fill = Freq
ggheatmap <- ggplot(frequency.dataframe.combined[,], aes(Var2, Var1, fill = frequency.dataframe.polyploidDE$Freq))+
  geom_tile(color = "white")+
  scale_fill_gradient(limits=c(0, max(frequency.dataframe.polyploidDE$Freq)),low = "white", high = "orange",space = "Lab", 
                       name="Number of\ntranscripts") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1))+
  coord_fixed()
ggheatmap + geom_text(aes(Var2, Var1, label = Freq), color = "black", size = 3)

col4 <- colorRampPalette(c("white","yellow", "purple"), space = "rgb")
D <- frequency.table[c("II","XI","IV","IX","V","VI","VIII","III","X","VII","I","XII"),c("II","XI","IV","IX","V","VI","VIII","III","X","VII","I","XII")]
P <- frequency.table.polyploidDE[c("II","XI","IV","IX","V","VI","VIII","III","X","VII","I","XII"),c("II","XI","IV","IX","V","VI","VIII","III","X","VII","I","XII")]

corrplot(D[1:12,1:12], is.corr = F, addrect=3, col=col4(50))#, method = "pie")
corrplot(P[1:12,1:12], is.corr = F, addrect=3, col=col4(50), method = "pie")

outputD <- data.frame(x1 = rowSums(D[,1:2]), x2 = rowSums(D[,3:4]), x3 = rowSums(D[,5:7]), x4 = rowSums(D[,8:10]), x5 = rowSums(D[,11:12]))
outputD2 <- data.frame(x1 = colSums(outputD[1:2,]), x2 = colSums(outputD[3:4,]), x3 = colSums(outputD[5:7,]), x4 = colSums(outputD[8:10,]), x5 = colSums(outputD[11:12,]))
corrplot(as.matrix(outputD2), is.corr = F, addrect=3, col=col4(50))#, method = "pie")

outputP <- data.frame(x1 = rowSums(P[,1:2]), x2 = rowSums(P[,3:4]), x3 = rowSums(P[,5:7]), x4 = rowSums(P[,8:10]), x5 = rowSums(P[,11:12]))
outputP2 <- data.frame(x1 = colSums(outputP[1:2,]), x2 = colSums(outputP[3:4,]), x3 = colSums(outputP[5:7,]), x4 = colSums(outputP[8:10,]), x5 = colSums(outputP[11:12,]))
corrplot(as.matrix(outputP2), is.corr = F, addrect=3, col=col4(100))#, method = "pie")

x <- outputP2/outputD2
x[is.na(x)] <- 0
corrplot(as.matrix(x), is.corr = F, addrect=3, col=col4(100))#, method = "pie")

# VIOL PLOT OF THE COUNTS IN DIFFERENT CATEGORIES
vioplot(rowSums(data[majalis.II,5:9])/5, rowSums(data[majalis.II, 15:25])/ rowSums(data[majalis.II, 10:14])/5)
#########################################################
#           PERMUTATION TESTS                           #
#########################################################
# SAVE CLUSTERS
for(i in c("II","XI","IV","IX","V","VI","VIII","III","X","VII","I","XII","XIII"))
{
  for(j in c("II","XI","IV","IX","V","VI","VIII","III","X","VII","I","XII","XIII"))
  {
      write.table(intersect(majalis.dataframe.ordered$geneID[majalis.dataframe.ordered$category==i], traunsteineri.dataframe.ordered$geneID[traunsteineri.dataframe.ordered$category==j]), paste("/home/botanik/Documents/phd/research/dominance/2017.03/enrichments/", i, j, "txt", sep="."), sep="\n", row.names = FALSE, col.names=FALSE, quote = F)
  }
}
