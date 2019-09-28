load("Chromium.rda")
tpm <- tpm[apply(tpm, 1, sd) > 0, ]
#data preparation

mds_tpm <- cmdscale(as.dist(1 - cor(tpm)))
#expression-level dimension reduction

library(viper)
load("../../MSigDB-regulon/Immuno-Hs-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset1 <- gset[len >= 190 & len <= 210]
gset1 <- gset1[grep("_UP", names(gset1))]
load("../../MSigDB-regulon/GO-BP-Hs-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset2 <- gset[len >= 50 & len <= 100]
gset <- c(gset1, gset2)
nes <- aREA(tpm, gset)$nes
#BPA

mds_nes <- cmdscale(as.dist(viperSimilarity(nes)))
#BPA-level dimension reduction

table <- unique(unlist(lapply(gset, function(x) names(x$tfmode))))
gsetRandom <- lapply(gset, function(x, table){
  gene <- sample(setdiff(table, names(x$tfmode)), length(x$tfmode))
  list(tfmode=structure(rep(1, length(gene)), names=gene), likelihood=rep(1, length(gene)))
}, table=table)
nes_r <- aREA(tpm, gsetRandom)$nes
#BPA on random genesets

library(dbscan)
dbscan <- dbscan(mds_nes, eps = 0.02)
#clustering analysis

b <- which(dbscan$cluster == 1)
t <- which(dbscan$cluster == 2)
mono <- which(dbscan$cluster == 3)
stat <- cbind(b=names(sort(rowMeans(nes[, b]) - rowMeans(nes[, -1*b]), decreasing = T)),
              t=names(sort(rowMeans(nes[, t]) - rowMeans(nes[, -1*t]), decreasing = T)),
              mono=names(sort(rowMeans(nes[, mono]) - rowMeans(nes[, -1*mono]), decreasing = T)))
write.csv(stat, file = "stat.csv", quote = F)
#differentially activated BPs

save(mds_tpm, nes, mds_nes, nes_r, dbscan, file = "GO.rda")
