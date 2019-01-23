# ----+ LIBRARIES
library("RUVSeq")
library("edgeR")
library("HTSFilter")
library("RColorBrewer")
library("adegenet")
library("statmod")
library("cgwtools")   

design = model.matrix(~specFac+normSet$W-1)
y = DGEList(counts = counts(set), group = specFac)
y = calcNormFactors(y, method = "upperquartile")
y = estimateGLMCommonDisp(y, design)
y = estimateGLMTagwiseDisp(y, design)

fitQL <- glmQLFit(y, design, robust = T)

confimt<-makeContrasts(
  fvsi=c(-1,1,0,0,0,0),
  fvsm=c(-1,0,1,0,0,0),
  fvst=c(-1,0,0,1,0,0),
  ivsm=c(0,-1,1,0,0,0),
  ivst=c(0,-1,0,1,0,0),
  mvst=c(0,0,-1,1,0,0),
  levels=design)

anov <- glmQLFTest(fitQL, contrast = confimt)

ModeldominanceFun = function(topTags1, FDRval, spec){ #(threewaycomparison, parents, fuchsiidiploid)
  tab1 <<- topTags1
  DE1 <<- topTags1$table[topTags1$table$FDR <= FDRval,]
  nonDE1 <<- rownames(topTags1[topTags1$table$FDR > FDRval,])
  
  
  UpFu1 <<- rownames(tab2$table[tab2$table$logFC <= 0.0 & tab2$table$FDR <= FDRval,])#++
  print(paste("UpFu1 size", as.character(length(UpFu1)), sep = " = "))
  
  UpIn1 <<- rownames(tab2$table[tab2$table$logFC > 0.0 & tab2$table$FDR <= FDRval,])#**
  print(paste("UpIn1 size", as.character(length(UpIn1)), sep = " = "))
  
  UpFuPoly <<- rownames(tab3$table[tab3$table$logFC <= 0.0 & tab3$table$FDR <= FDRval,])#++
  print(paste("UpFuPoly size", as.character(length(UpFuPoly)), sep = " = "))
  
  UpPolyFu <<- rownames(tab3$table[tab3$table$logFC > 0.0 & tab3$table$FDR <= FDRval,])#**
  print(paste("UpPolyFu size", as.character(length(UpPolyFu)), sep = " = "))
  
  #MAJALIS CATEGORIES
  I <- rownames(DE1[DE1[1] < 0 & DE1[2] < 0 & DE1[2]-DE1[1] > 0,])# _-' (inc, poly, fuc)
  assign(paste(as.character(spec), "I", sep = "."), I, envir = .GlobalEnv)
  
  III <- rownames(DE1[DE1[1] < 0 & DE1[2] < 0 & DE1[2]-DE1[1] < 0,])# -_'
  assign(paste(as.character(spec), "III", sep = "."), III, envir = .GlobalEnv)
  
  V <- rownames(DE1[DE1[1] < 0 & DE1[2] > 0 & DE1[2]-DE1[1] > 0,])# _'-
  assign(paste(as.character(spec), "V", sep = "."), V, envir = .GlobalEnv)
  
  VI <- rownames(DE1[DE1[1] > 0 & DE1[2] > 0 & DE1[2]-DE1[1] > 0,])# -'_
  assign(paste(as.character(spec), "VI", sep = "."), VI, envir = .GlobalEnv)
  
  X <- rownames(DE1[DE1[1] > 0 & DE1[2] < 0 & DE1[2]-DE1[1] < 0,])# '_-
  assign(paste(as.character(spec), "X", sep = "."), X, envir = .GlobalEnv)
  
  XII <- rownames(DE1[DE1[1] > 0 & DE1[2] > 0 & DE1[2]-DE1[1] < 0,])# '-_
  assign(paste(as.character(spec), "XII", sep = "."), XII, envir = .GlobalEnv)
  
  II <- Reduce(intersect, list(nonDE1, UpFu1, nonDE3))# _--
  assign(paste(as.character(spec), "II", sep = "."), II, envir = .GlobalEnv)
  
  IV <- Reduce(intersect, list(nonDE1, UpIn1, UpPolyFu))# --_
  assign(paste(as.character(spec), "IV", sep = "."), IV, envir = .GlobalEnv)
  
  
  VII <- Reduce(intersect, list(nonDE1, nonDE2, UpFuPoly))# -_-
  assign(paste(as.character(spec), "VII", sep = "."), VII, envir = .GlobalEnv)
  
  VIII <- Reduce(intersect, list(nonDE1,nonDE2, UpPolyFu))# _-_
  assign(paste(as.character(spec), "VIII", sep = "."), VIII, envir = .GlobalEnv)
  
  IX <- Reduce(intersect, list(nonDE1, UpFu1, UpFuPoly))# __-
  assign(paste(as.character(spec), "IX", sep = "."), IX, envir = .GlobalEnv)
  
  XI <- Reduce(intersect, list(nonDE1,UpIn1, nonDE3))# -__
  assign(paste(as.character(spec), "XI", sep = "."), XI, envir = .GlobalEnv)
  
  # TRAUNSTEINERI CATEGORIES
  
  XIII <- nonDE1# ---
  assign(paste(as.character(spec), "XIII", sep = "."), XIII, envir = .GlobalEnv)
  
  cats <- c(I, II, III, IV, V, VI, VII, VIII, IX, X, XI, XII, XIII)
  
  print(length(c(I, II, III, IV, V, VI, VII, VIII, IX, X, XI, XII, XIII)))
  
  dataframe <- data.frame("geneID" = c(I, II, III, IV, V, VI, VII, VIII, IX, X, XI, XII, XIII),
                          "logFC.fuchsii" = tab1$table[cats,2], "logFC.inacarnata" = tab1$table[cats,2]-tab1$table[cats,1],
                          "category" = c(
                            rep("I", length(I)), rep("II", length(II)), rep("III", length(III)), 
                            rep("IV", length(IV)), rep("V", length(V)), rep("VI", length(VI)), 
                            rep("VII", length(VII)), rep("VIII", length(VIII)), rep("IX", length(IX)),
                            rep("X", length(X)), rep("XI", length(XI)), rep("XII", length(XII)),
                            rep("XIII", length(XIII)) )
  )
  
  dataframe <- data.frame(dataframe, "color" = map.category.color[as.character(dataframe$category)])
  
  assign(paste(as.character(spec), "dataframe", sep = "."), dataframe, envir = .GlobalEnv)
}
