library(viper)
load("../Embryo/E-MTAB-3929/rpkm.rda")
rpkm <- rpkm[apply(rpkm, 1, sd) > 0, ]
load("../MSigDB-regulon/GO-BP-Hs-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset <- gset[len >= 50 & len <= 100]
nes1 <- viper(rpkm[, 1:20], gset, method = "none", pleiotropy = T)
nes2 <- viper(rpkm[, 1:20], gset, method = "none", pleiotropy = F)
save(nes1, nes2, file = "shadow.rda")
pdf("shadow.pdf")
plot(rowMeans(nes1), rowMeans(nes2), xlab = "Shadow=T", ylab = "Shadow=F", pch = 16, main = "Average Activity of 20 Ramdom Cells")
abline(0, 1, NULL, NULL, col = 2, lty = 4)
dev.off()
#shadow analysis

pdf("gsetSize.pdf", width = 9, height = 6)
par(mfrow = c(2, 3), mar = c(4, 4, 4, 2))
load("../MSigDB-regulon/GO-BP-Hs-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
plot(density(log10(len)), xlab = "log10(#Genes/Geneset)", main = "GO-BP-Hs-MSigDB")
abline(NULL, NULL, NULL, log10(c(50, 100)), col = 2, lty = 4)
legend("topright", legend = c(paste("Total:", length(len)), paste("Remaining:", sum(len >= 50 & len <= 100))), bty = "n")
load("../MSigDB-regulon/GO-BP-Mm-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
plot(density(log10(len)), xlab = "log10(#Genes/Geneset)", main = "GO-BP-Mm-MSigDB")
abline(NULL, NULL, NULL, log10(c(50, 100)), col = 2, lty = 4)
legend("topright", legend = c(paste("Total:", length(len)), paste("Remaining:", sum(len >= 50 & len <= 100))), bty = "n")
load("../MSigDB-regulon/GO-BP-Dr-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
plot(density(log10(len)), xlab = "log10(#Genes/Geneset)", main = "GO-BP-Dr-MSigDB")
abline(NULL, NULL, NULL, log10(c(50, 100)), col = 2, lty = 4)
legend("topright", legend = c(paste("Total:", length(len)), paste("Remaining:", sum(len >= 50 & len <= 100))), bty = "n")
#BP
load("../MSigDB-regulon/Immuno-Hs-MSigDB-regulon.rda")
len <- unlist(lapply(gset[grep("_UP", names(gset))], function(x) length(x$tfmode)))
plot(density(log10(len)), xlab = "log10(#Genes/Geneset)", main = "Immuno-Hs-MSigDB (UP)")
abline(NULL, NULL, NULL, log10(c(190, 210)), col = 2, lty = 4)
legend("topleft", legend = c(paste("Total:", length(len)), paste("Remaining:", sum(len >= 190 & len <= 210))), bty = "n")
load("../MSigDB-regulon/Immuno-Mm-MSigDB-regulon.rda")
len <- unlist(lapply(gset[grep("_UP", names(gset))], function(x) length(x$tfmode)))
plot(density(log10(len)), xlab = "log10(#Genes/Geneset)", main = "Immuno-Mm-MSigDB (UP)")
abline(NULL, NULL, NULL, log10(c(190, 210)), col = 2, lty = 4)
legend("topleft", legend = c(paste("Total:", length(len)), paste("Remaining:", sum(len >= 190 & len <= 210))), bty = "n")
#Immuno
dev.off()
#gene set size bias

library(viper)
load("../Embryo/GSE65525/tpm.rda")
load("../MSigDB-regulon/GO-BP-Hs-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
nes <- aREA(tpm, gset)$nes
pdf("gsetFilter.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4, 4, 2))
plot(log10(len[rownames(nes)]), abs(rowMeans(nes)), pch = 16, cex = 0.5,
     xlab = "log10(Size, Geneset)", ylab = "Average Abs(Activity)", main = "Mouse ES")
plot(log10(len[rownames(nes)]), log10(apply(nes, 1, sd)/abs(rowMeans(nes))), pch = 16, cex = 0.5,
     xlab = "log10(Size, Geneset)", ylab = "log10(CV, Activity)", main = "Mouse ES")
dev.off()
#gene set size selection

load("../Immune/Chromium/GO.rda")
pdf("gsetRandom.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))
plot(density(apply(nes, 1, sd)), col = 1, ylim = c(0, 8), main = "Density Distribution of SD/Term")
lines(density(apply(nes_r, 1, sd)), col = 4)
plot(density(rowMeans(nes)), col = 1, ylim = c(0, 1.5), xlab = "", main = "Density Distribution of Mean/Term")
lines(density(rowMeans(nes_r)), col = 4)
legend("topright", legend = c("Original", "Random"), fill = c(1, 4),  bty = "n", border = NA)
dev.off()
#BPA with random gene sets

library(viper)
library(Rtsne)
lung <- get(load("../Dropout/GTEx/Lung.rda"))
eso <- get(load("../Dropout/GTEx/Esophagus.rda"))
dset1 <- cbind(lung[, 1:10], eso[, 1:90])
dset2 <- cbind(lung[, 101:150], eso[, 101:150])
load("../MSigDB-regulon/GO-BP-Hs-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset <- gset[len >= 50 & len <= 100]

pdf("normalization.pdf", width = 8, height = 4)
par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
rank <- apply(dset1, 2, rank)
nes1 <- aREA((rank - apply(rank, 1, median))/apply(rank, 1, mad), gset)$nes
rank <- apply(dset2, 2, rank)
nes2 <- aREA((rank - apply(rank, 1, median))/apply(rank, 1, mad), gset)$nes
image(cor(cbind(nes1, nes2)), col = cm.colors(100), axes = F, main = "Internal Normalization")
axis(side = 1, line = -0.5, at = seq(0, 1, length.out = 200)[c(5, 55, 125, 175)], labels = c("Lung", "Eso", "Lung'", "Eso'"), tick = F, las = 2)
axis(side = 2, line = -0.5, at = seq(0, 1, length.out = 200)[c(5, 55, 125, 175)], labels = c("Lung", "Eso", "Lung'", "Eso'"), tick = F, las = 2)
abline(NULL, NULL, c(0, 0.5, 1), c(0, 0.5, 1), lty = 1)
abline(NULL, NULL, c(0.05, 0.75), c(0.05, 0.75), lty = 4)

nes1 <- aREA(dset1, gset)$nes
nes2 <- aREA(dset2, gset)$nes
image(cor(cbind(nes1, nes2)), col = cm.colors(100), axes = F, main = "Raw")
axis(side = 1, line = -0.5, at = seq(0, 1, length.out = 200)[c(5, 55, 125, 175)], labels = c("Lung", "Eso", "Lung'", "Eso'"), tick = F, las = 2)
axis(side = 2, line = -0.5, at = seq(0, 1, length.out = 200)[c(5, 55, 125, 175)], labels = c("Lung", "Eso", "Lung'", "Eso'"), tick = F, las = 2)
abline(NULL, NULL, c(0, 0.5, 1), c(0, 0.5, 1), lty = 1)
abline(NULL, NULL, c(0.05, 0.75), c(0.05, 0.75), lty = 4)
dev.off()
#"absolute" and "relative" activity

load("../Immune/Chromium/Chromium.rda")
col <- tpm[c("CD3E", "CD14", "MS4A1"), ]
col <- (col - apply(col, 1, min))/apply(col, 1, sd)
col <- apply(col, 2, prop.table)
col[!is.finite(col)] <- 0
col <- rgb(col[1, ], col[2, ], col[3, ])
pdf("batchCorrect.pdf", width = 9, height = 3)
par(mfrow = c(1, 3), mar = c(4, 4, 4, 2))
load("./mnnCorrect_tsne_coords.rda")
plot(joint_tsne$Y, col = col, pch = c(rep(3, 500), rep(1, 500)), cex = 0.5, xlab = "Dim1", ylab = "Dim2", main = "MNN", cex.main = 2)
load("./scanorama_tsne_coords.rda")
plot(joint_tsne$Y, col = col, pch = c(rep(3, 500), rep(1, 500)), cex = 0.5, xlab = "Dim1", ylab = "Dim2", main = "scanorama", cex.main = 2)
load("./scMerge_tsne_coords.rda")
plot(joint_tsne$Y, col = col, pch = c(rep(3, 500), rep(1, 500)), cex = 0.5, xlab = "Dim1", ylab = "Dim2", main = "scMerge", cex.main = 2)
legend("bottomright", legend = c("V1", "V2"), pch = c(3, 1), cex = 1, bty = "n", border = NA)
legend("bottomleft", legend = c("T", "Mono", "B"), fill = c("Red", "Green", "Blue"), cex = 1, bty = "n", border = NA)
dev.off()
#batch correction algorithms

library(made4)
load("./cca_tsne_coords.rda")
pdf("seurat.pdf")
plot(tsne_coords, col = getcol(21)[as.factor(meta$species)], pch = 16, cex = 0.5, xlab = "Dim1", ylab = "Dim2", main = "Seurat", cex.main = 2)
legend("bottomright", legend = levels(as.factor(meta$species)), fill = getcol(21), cex = 1, bty = "n", border = NA)
dev.off()
#seurat

library(viper)
library(Rtsne)
gset1 <- get(load("../MSigDB-regulon/Immuno-Hs-MSigDB-regulon.rda"))
len <- unlist(lapply(gset1, function(x) length(x$tfmode)))
gset1 <- gset1[len >= 190 & len <= 210]
gset1 <- gset1[grep("_UP", names(gset1))]
gset2 <- get(load("../MSigDB-regulon/GO-BP-Hs-MSigDB-regulon.rda"))
len <- unlist(lapply(gset2, function(x) length(x$tfmode)))
gset2 <- gset2[len >= 50 & len <= 100]
gset <- c(gset1, gset2)

load("../Immune/Chromium/Chromium.rda")
col <- tpm[c("CD3E", "CD14", "MS4A1"), ]
col <- (col - apply(col, 1, min))/apply(col, 1, sd)
col <- apply(col, 2, prop.table)
col[!is.finite(col)] <- 0
col <- rgb(col[1, ], col[2, ], col[3, ])
nes <- aREA(tpm, gset)$nes

tsne_all <- Rtsne(as.dist(viperSimilarity(nes)), is_distance = T, pca = F, max_iter = 500, perplexity = 50)$Y
tsne_500 <- Rtsne(as.dist(viperSimilarity(nes[order(apply(nes, 1, sd), decreasing = T)[1:500], ])), is_distance = T, pca = F, max_iter = 500, perplexity = 50)$Y
tsne_100 <- Rtsne(as.dist(viperSimilarity(nes[order(apply(nes, 1, sd), decreasing = T)[1:100], ])), is_distance = T, pca = F, max_iter = 500, perplexity = 50)$Y
save(tsne_all, tsne_500, tsne_100, file = "variableTerms.rda")

pdf("variableTerms.pdf", width = 8, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 4, 1))
plot(tsne_all, col = col, pch = 16, cex = 0.75, xlab = "Dim1", ylab = "Dim2", main = "All Terms")
plot(tsne_500, col = col, pch = 16, cex = 0.75, xlab = "Dim1", ylab = "Dim2", main = "Top 500 Variable Terms")
plot(tsne_100, col = col, pch = 16, cex = 0.75, xlab = "Dim1", ylab = "Dim2", main = "Top 100 Variable Terms")
legend("bottomleft", legend = c("T", "Mono", "B"), fill = c("Red", "Green", "Blue"), bty = "n", border = NA)
plot(density(apply(nes[1:length(gset1), ], 1, sd)), ylim = c(0, 2), xlab = "", main = "SD/Term")
lines(density(apply(nes[(length(gset1)+1):(length(gset1)+length(gset2)), ], 1, sd)), col = 2)
legend("topright", legend = c("Immuno", "GO-BP"), fill = 1:2, bty = "n", border = NA)
dev.off()
#variable terms
