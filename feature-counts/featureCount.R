library("Rsubread", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
load("/home/botanik/Documents/phd/scripts/Rscripts/DEclean.R")

filenamesP1 <- list.files(path="/data/phdData/orchis/alignments/bamFiles/sorted/P1",
                          pattern="*.clean.RG.P1.bam", recursive = T, full.names = T)
filenamesP2 <- list.files(path="/data/phdData/orchis/alignments/bamFiles/sorted/P2",
                          pattern="*.clean.RG.P2.bam", recursive = T,  full.names = T)


filenames <- c(filenamesP1, filenamesP2)

fc_HEB <- featureCounts(filenames,annot.ext="/data/phdData/orchis/blast2go.annotation.Orchis.gff",isGTFAnnotationFile=T,GTF.featureType="CDS",isPairedEnd=TRUE,genome="/data/phdData/orchis/Oita_infl_unigenes.unleaved.fasta", GTF.attrType="ID")
data.featureCounts <- as.data.frame(fc_HEB$counts)

P1.names <- c("P1.tB1830_1","P1.tB1830_2","P1.tB1833_1","P1.tB1833_2",
              "P1.fB1804","P1.fB1855","P1.fP1001","P1.fP1707_1","P1.fP1707_2",
              "P1.iA1586","P1.iB1176","P1.iB1870","P1.iS1904","P1.iS1908",
              "P1.mA1567_1","P1.mA1567_2","P1.mA1568","P1.mA1573","P1.mA1661","P1.mA1775","P1.mP1722_1","P1.mP1722_2","P1.mP1744","P1.mP1748_1","P1.mP1748_2","P1.mS1757","P1.mS1765_1","P1.mS1765_2",
              "P1.tA1553","P1.tA1641","P1.tA1670","P1.tB1798_1","P1.tB1798_2","P1.tB1805_1","P1.tB1805_2","P1.tB1812_1","P1.tB1812_2","P1.tS1901","P1.tS1902","P1.tS1920_1","P1.tS1920_2")

P2.names <- c("P2.tB1830_1","P2.tB1830_2","P2.tB1833_1","P2.tB1833_2",
              "P2.fB1804","P2.fB1855","P2.fP1001","P2.fP1707_1","P2.fP1707_2",
              "P2.iA1586","P2.iB1176","P2.iB1870","P2.iS1904","P2.iS1908",
              "P2.mA1567_1","P2.mA1567_2","P2.mA1568","P2.mA1573","P2.mA1661","P2.mA1775","P2.mP1722_1","P2.mP1722_2","P2.mP1744","P2.mP1748_1","P2.mP1748_2","P2.mS1757","P2.mS1765_1","P2.mS1765_2",
              "P2.tA1553","P2.tA1641","P2.tA1670","P2.tB1798_1","P2.tB1798_2","P2.tB1805_1","P2.tB1805_2","P2.tB1812_1","P2.tB1812_2","P2.tS1901","P2.tS1902","P2.tS1920_1","P2.tS1920_2")

nam <- c(P1.names, P2.names)

colnames(data.featureCounts) <- nam
data.featureCounts[colnames(data.featureCounts) == "P1.mA1567_1" | colnames(data.featureCounts) == "P1.mA1567_2" | colnames(data.featureCounts) == "P2.mA1567_1" | colnames(data.featureCounts) == "P2.mA1567_2"] <- NULL
data.featureCounts[colnames(data.featureCounts) == "P1.mP1748_1" | colnames(data.featureCounts) == "P1.mP1748_2" | colnames(data.featureCounts) == "P2.mP1748_1" | colnames(data.featureCounts) == "P2.mP1748_2"] <- NULL

data.featureCounts <- data.featureCounts[order(colnames(data.featureCounts))]

# SELECT COMPARISONS OF INTEREST
dat <- selectSpecies(data.featureCounts, "P1.m", "P1.t", " "," "," ")


