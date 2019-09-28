c1 <- get(load("Spleen_cpm.rda"))
colnames(c1) <- paste("Spleen", 1:ncol(c1), sep = "_")
c2 <- get(load("Thymus_cpm.rda"))
colnames(c2) <- paste("Thymus", 1:ncol(c2), sep = "_")
cpm <- cbind(c1, c2)
cpm <- cpm[apply(cpm, 1, sd) > 0, ]
#data preparation

mds_cpm <- cmdscale(as.dist(1 - cor(cpm)))
#expression-level dimension reduction

library(viper)
load("../../MSigDB-regulon/Immuno-Mm-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset1 <- gset[len >= 190 & len <= 210]
gset1 <- gset1[grep("_UP", names(gset1))]
load("../../MSigDB-regulon/GO-BP-Mm-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset2 <- gset[len >= 50 & len <= 100]
gset <- c(gset1, gset2)
nes <- aREA(cpm, gset)$nes
#BPA

mds_nes <- cmdscale(as.dist(viperSimilarity(nes)))
#BPA-level dimension reduction

table <- unique(unlist(lapply(gset, function(x) names(x$tfmode))))
gsetRandom <- lapply(gset, function(x, table){
  gene <- sample(setdiff(table, names(x$tfmode)), length(x$tfmode))
  list(tfmode=structure(rep(1, length(gene)), names=gene), likelihood=rep(1, length(gene)))
}, table=table)
nes_r <- aREA(cpm, gsetRandom)$nes
#BPA on random genesets

library(dbscan)
dbscan <- dbscan(mds_nes, eps = 0.0125)
#clustering analysis

b <- which(dbscan$cluster == 1)
t <- which(dbscan$cluster == 2)
stat <- cbind(b=names(sort(rowMeans(nes[, b]) - rowMeans(nes[, t]), decreasing = T)),
              t=names(sort(rowMeans(nes[, t]) - rowMeans(nes[, b]), decreasing = T)))
write.csv(stat, file = "stat.csv", quote = F)
#differentially activated BPs

save(mds_cpm, nes, mds_nes, nes_r, dbscan, file = "GO.rda")
