name_list <- c(
  "tB1830_1", "tB1830_2", "tB1833_1", "tB1833_2", "fB1804", "fB1855",
  "fP1001", "fP1707_1", "fP1707_2", "iA1586", "iB1176", "iB1870",
  "iS1904", "iS1908", "mA1567_1", "mA1567_2", "mA1568", "mA1573",
  "mA1661", "mA1775", "mP1722_1", "mP1722_2", "mP1744", "mP1748_1",
  "mP1748_2", "mS1757", "mS1765_1", "mS1765_2", "tA1553", "tA1641",
  "tA1670", "tB1798_1", "tB1798_2", "tB1805_1", "tB1805_2",
  "tB1812_1", "tB1812_2", "tS1901", "tS1902", "tS1920_1", "tS1920_2"
)
remove_list <- c("mA1567_1", "mA1567_2", "mP1748_1", "mP1748_2")
data.featureCounts <- read.table("count.featureCounts.txt")
des.names <- c("majalis", "traunsteineri", "nW1", "nW2", "nW3")
data.featureCounts.names <- change_names(data = data.featureCounts, name_list = name_list)
data.featureCounts.clean <- remove_id(data.featureCounts.names, remove_list)
dat <- selectSpecies(data.featureCounts.clean, "t", "m", " ", " ", " ")
species <- makeSpeciesVector(dat)
keep <- rowSums(cpm(dat) > 0.08) >= 10
batchMatrix <- makeBatchMatrix(data = dat)
specFac <- as.factor(species)
dat.filtered <- dat[keep, ]
s <- summary(keep)
set <- makeRUVset(dat = dat.filtered)
counts.set <- rownames(counts(set))
y0 <- DGEList(counts = counts(set), group = specFac)
normSet <- makeRUVrepNormalisedSet(set, counts.set, 3, batchMatrix)
y1 <- calcNormFactors(y0)
des <- data.frame(specFac, normSet$W)
design <- model.matrix(~specFac + W_1 + W_2 + W_3 - 1, data = des)
y2 <- estimateDisp(y1, design = design, robust = T)
design.names <- change_names(data = design, name_list = des.names)
fit <- glmQLFit(y2, design, robust = T)
pBCV <- plotBCV(y2)
hist <- hist(fit$coefficients[, 1], breaks = 100)
sFit <- summary(fit$df.prior)
pQLDisp <- plotQLDisp(fit)
