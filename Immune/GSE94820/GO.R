load("GSE94820.rda")
load("anno.rda")
ind <- colSums(expmat > 0) > 2000 & colSums(expmat > 0) < 8000
expmat <- expmat[, ind]
anno <- anno[ind]
#data preparation

library(viper)
load("../../MSigDB-regulon/Immuno-Hs-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset1 <- gset[len >= 190 & len <= 210]
gset1 <- gset1[grep("_UP", names(gset1))]
load("../../MSigDB-regulon/GO-BP-Hs-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset2 <- gset[len >= 50 & len <= 100]
gset <- c(gset1, gset2)
nes <- aREA(expmat, gset)$nes
#BPA

library(made4)
mds_nes <- cmdscale(as.dist(viperSimilarity(nes)))
#BPA-level dimension reduction

table <- unique(unlist(lapply(gset, function(x) names(x$tfmode))))
gsetRandom <- lapply(gset, function(x, table){
  gene <- sample(setdiff(table, names(x$tfmode)), length(x$tfmode))
  list(tfmode=structure(rep(1, length(gene)), names=gene), likelihood=rep(1, length(gene)))
}, table=table)
nes_r <- aREA(expmat, gsetRandom)$nes
#BPA on random genesets

mono <- grep("Mono", anno)
pdc <- grep("pDC", anno)
cdc <- c(grep("CD141", anno), grep("CD1C", anno), grep("DoubleNeg", anno))
stat <- cbind(mono=names(sort(rowMeans(nes[, mono]) - rowMeans(nes[, -1*mono]), decreasing = T)),
              pdc=names(sort(rowMeans(nes[, pdc]) - rowMeans(nes[, -1*pdc]), decreasing = T)),
              cdc=names(sort(rowMeans(nes[, cdc]) - rowMeans(nes[, -1*cdc]), decreasing = T)))
write.csv(stat, file = "stat.csv", quote = F)
#differentially activated BPs

save(nes, mds_nes, nes_r, file = "GO.rda")
