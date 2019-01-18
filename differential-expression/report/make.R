library("gplots")
library("tidyverse")
library("broom")
library("here")
library("data.table")
library("RUVSeq")
library("edgeR")
library("statmod")
library("RColorBrewer")
library("drake")
pkgconfig::set_config("drake::strings_in_dots" = "literals")

## if you need to rerun everything
## drake::clean(destroy = TRUE)

## create softlinks for data
## ---------------------------

#fs::link_create("/home/botanik/Documents/phd/scripts/Rscripts/differential-expression/data/count.featureCounts.txt", "count.featureCounts.txt")
## GLOBAL VARIABLES & FUNCTION FOR DATA ORGANISATION
#data <- read.table("/home/botanik/Documents/phd/research/siRNA/differential.expression/counts_sRNA20-22.txt", header = TRUE, row.names = 1)
#data.orchis <- read.table("/data/phdData/orchis/counts/counts.corset.txt", header = TRUE, row.names = 1)
#data.hylite <- read.table("/data/phdData/orchis/hylite/HyLiTE_analysis/HyLiTE_analysis.expression.txt", header = TRUE, row.names = 1)
#chloro.mito <- read.table("/data/phdData/orchis/counts/de.featureCounts/HEB/chloroplastic.and.mitochondrial.genes.txt")

#fs::link_create(here("..", "data", "genotyping_data_subset_train.raw"), "genotyping_data_subset_train.raw")

## check supporting files
## -------------------------

file.exists(here(path = "..", "src", "functions.R"))
file.exists("report.Rmd")

# Your custom code is a bunch of functions.
## -------------------------
source("/home/thomas/Documents/phd/scripts/differential-expression/src/functions.R")
source(here(path = "..", "src", "functions.R"))

# The workflow plan data frame outlines what you are going to do.

plan <- drake_plan(
  name_list = c("tB1830_1","tB1830_2","tB1833_1","tB1833_2", 
                 "fB1804","fB1855","fP1001","fP1707_1","fP1707_2", 
                 "iA1586","iB1176","iB1870","iS1904","iS1908", 
                 "mA1567_1","mA1567_2","mA1568","mA1573","mA1661","mA1775","mP1722_1","mP1722_2","mP1744","mP1748_1","mP1748_2","mS1757","mS1765_1","mS1765_2", 
                 "tA1553","tA1641","tA1670","tB1798_1","tB1798_2","tB1805_1","tB1805_2","tB1812_1","tB1812_2","tS1901","tS1902","tS1920_1","tS1920_2"),
  
  remove_list = c("mA1567_1", "mA1567_2", "mP1748_1", "mP1748_2"),
  data.featureCounts = read.table("count.featureCounts.txt"),
  data.featureCounts.names = change_names(data = data.featureCounts, name_list = name_list),
  data.featureCounts.clean = remove_id(data.featureCounts.names, remove_list),
  dat = selectSpecies(data.featureCounts.clean, "t", "m", " "," "," "),
  species = makeSpeciesVector(dat),
  specFac = as.factor(species),
  batchMatrix = makeBatchMatrix(data = dat),
  #check filtering step
  keep = rowSums(cpm(dat)>0.08) >= 10,
  s = summary(keep),
  dat.filtered = dat[keep,],
  set = makeRUVset(dat = dat.filtered),
  counts.set = rownames(counts(set)),
  normSet = makeRUVrepNormalisedSet(set, counts.set, 3, batchMatrix),
  des = data.frame(specFac, normSet$W),
  #formula = "~ specFac + normSet$W",
  #formul = as.formula("~species+normSetW-1"),
  design = model.matrix(~specFac + W_1 + W_2 + W_3 - 1, data = des),
  #design = model.matrix(~geoFac+normSet$W-1)
  des.names = c("majalis", "traunsteineri", "nW1", "nW2", "nW3"),
  design.names = change_names(data = design, name_list = des.names),
  y0 = DGEList(counts = counts(set), group = specFac), # my stuff: counts(set), group = group) # Don't forget, model is taking into account the RUV normalised matrix of counts
  #get an idea of how to filter
  #cpm = cpm(5, mean(y0$samples$lib.size)),
  y1 = calcNormFactors(y0),
  y2 = estimateDisp(y1, design = design, robust = T),
  pBCV = plotBCV(y2),
  fit = glmQLFit(y2, design, robust = T),
  hist = hist(fit$coefficients[,1], breaks = 100),
  sFit = summary(fit$df.prior),
  pQLDisp = plotQLDisp(fit),
  con = c(-1,1,0,0,0), #Where the contrasts are chosen
  resEdge = glmQLFTest(fit, contrast = con),
  de = decideTestsDGE(resEdge, adjust.method="fdr", p.value = 0.05),
  is.de = decideTestsDGE(resEdge),
  pMD = plotMD(resEdge, status=is.de, values=c(1,0,-1), col=c("blue","black", "red"), legend="topright", pch=19, cex=c(0.5,0.1)),
  resEdgeRtopTags = topTags(resEdge, n=nrow(set)),
  FDRres = topTags(resEdgeRtopTags, n=nrow(dat))$table
  #con <- makeContrasts(traunsteineri - majalis, levels = design)
  )
make(plan)
