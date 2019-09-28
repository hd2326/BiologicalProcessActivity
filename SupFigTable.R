system("cp ./Benchmark/batchCorrect.pdf ./SupFig1.pdf")
#SupFig1

load("./Dropout/GTEx/GO.rda")
load("./Immune/Chromium/GO.rda")
load("./Immune/Chromium/Chromium.rda")
load("./Embryo/E-MTAB-3929/GO.rda")
pdf("SupFig2.pdf", width = 9, height = 6)
par(mfcol = c(2, 3))
plot(rowMeans(rpkm_drop[, 1:20]), rowSums(rpkm_drop[, 1:20] > 0)/20, xlab = "Average Expression", ylab = "Expression Rate", pch = 16, main = "Esophagus")
legend("bottomright", legend = paste("Overall Drop-out(%)", signif(sum(rpkm_drop[, 1:20] == 0)/length(rpkm_drop[, 1:20])*100, 4)), bty = "n", border = NA)
plot(rowMeans(rpkm_drop[, 21:40]), rowSums(rpkm_drop[, 21:40] > 0)/20, xlab = "Average Expression", ylab = "Expression Rate", pch = 16, main = "Lung")
legend("bottomright", legend = paste("Overall Drop-out(%)", signif(sum(rpkm_drop[, 21:40] == 0)/length(rpkm_drop[, 21:40])*100, 4)), bty = "n", border = NA)
plot(rowMeans(tpm[, dbscan$cluster == 1]), rowSums(tpm[, dbscan$cluster == 1] > 0)/sum(dbscan$cluster == 1),
     xlab = "Average Expression", ylab = "Expression Rate", pch = 16, main = "Chromium B-Cell")
legend("bottomright", legend = paste("Overall Drop-out(%)", signif(sum(tpm[, dbscan$cluster == 1] == 0)/length(tpm[, dbscan$cluster == 1])*100, 4)), bty = "n", border = NA)
plot(rowMeans(tpm[, dbscan$cluster == 2]), rowSums(tpm[, dbscan$cluster == 2] > 0)/sum(dbscan$cluster == 2),
     xlab = "Average Expression", ylab = "Expression Rate", pch = 16, main = "Chromium T-Cell")
legend("bottomright", legend = paste("Overall Drop-out(%)", signif(sum(tpm[, dbscan$cluster == 2] == 0)/length(tpm[, dbscan$cluster == 2])*100, 4)), bty = "n", border = NA)
plot(rowMeans(rpkm[, colnames(rpkm) == "E6"]), rowSums(rpkm[, colnames(rpkm) == "E6"] > 0)/sum(colnames(rpkm) == "E6"),
     xlab = "Average Expression", ylab = "Expression Rate", pch = 16, main = "E-MTAB-3929 E6")
legend("bottomright", legend = paste("Overall Drop-out(%)", signif(sum(rpkm[, colnames(rpkm) == "E6"] == 0)/length(rpkm[, colnames(rpkm) == "E6"])*100, 4)), bty = "n", border = NA)
plot(rowMeans(rpkm[, colnames(rpkm) == "E7"]), rowSums(rpkm[, colnames(rpkm) == "E7"] > 0)/sum(colnames(rpkm) == "E7"),
     xlab = "Average Expression", ylab = "Expression Rate", pch = 16, main = "E-MTAB-3929 E7")
legend("bottomright", legend = paste("Overall Drop-out(%)", signif(sum(rpkm[, colnames(rpkm) == "E7"] == 0)/length(rpkm[, colnames(rpkm) == "E7"])*100, 4)), bty = "n", border = NA)
dev.off()
#SupFig2

source("./ColorGradient.R")
pdf("SupFig3.pdf", width = 12, height = 9)
par(mfrow = c(3, 4))
load("./Immune/Chromium/GO.rda")
load("./Immune/Chromium/Chromium.rda")
plot(mds_nes, col = c(1, 4, 5, 6)[dbscan$cluster+1], pch = 16, cex = 0.75, xlab = "Dim1", ylab = "Dim2", main = "Chromium")
legend("topright", legend = c("Outlier", "B", "T", "Mono"), fill = c(1, 4, 5, 6),  bty = "n", border = NA)
plot(mds_nes, col = expColor("CD3E", tpm), pch = 16, cex = 0.75, xlab = "Dim1", ylab = "Dim2", main = "T-cell, CD3E")
plot(mds_nes, col = expColor("CD14", tpm), pch = 16, cex = 0.75, xlab = "Dim1", ylab = "Dim2", main = "Monocyte, CD14")
plot(mds_nes, col = expColor("MS4A1", tpm), pch = 16, cex = 0.75, xlab = "Dim1", ylab = "Dim2", main = "B-cell, CD20")
#Chromium
load("./Immune/QuakeMCA/GO.rda")
c1 <- get(load("./Immune/QuakeMCA/Spleen_cpm.rda"))
c2 <- get(load("./Immune/QuakeMCA/Thymus_cpm.rda"))
cpm <- cbind(c1, c2)
plot(mds_nes, col = c(1, 4, 5, rep(1, 5))[dbscan$cluster+1], pch = 16, cex = 0.75, xlab = "Dim1", ylab = "Dim2", main = "Tabula Muris")
legend("bottomright", legend = c("Outlier", "B", "T"), fill = c(1, 4, 5),  bty = "n", border = NA)
plot(mds_nes, col = expColor("Cd3e", cpm), pch = 16, cex = 0.75, xlab = "Dim1", ylab = "Dim2", main = "T-cell, Cd3e")
plot(mds_nes, col = expColor("Cd14", cpm), pch = 16, cex = 0.75, xlab = "Dim1", ylab = "Dim2", main = "Monocyte, Cd14")
plot(mds_nes, col = expColor("Ms4a1", cpm), pch = 16, cex = 0.75, xlab = "Dim1", ylab = "Dim2", main = "B-cell, Cd20")
#QuakeMCA
load("./Immune/GSE94820/GO.rda")
anno <- sapply(strsplit(colnames(nes), split = "_"), function(x) x[1])
plot(mds_nes, col = getcol(21)[as.factor(anno)], pch = 16, cex = 0.75, xlab = "Dim1", ylab = "Dim2", xlim = c(-0.3, 0.3), ylim = c(-0.2, 0.2), main = "GSE94820")
legend("bottomleft", legend = levels(as.factor(anno)), fill = getcol(21), bty = "n", border = NA)
#GSE94820
dev.off()
#SupFig3

library(princurve)
load("./Embryo/E-MTAB-3929/GO.rda")
pdf("SupFig4.pdf", width = 12, height = 4)
par(mfcol = c(1, 3), mar = c(5, 5, 5, 1))
pc_rpkm <- princurve::principal.curve(tsne_rpkm, smoother = "lowess")
plot(tsne_rpkm, col = rainbow(5)[as.factor(colnames(nes))], cex = 0.75, pch = 16, xlab = "Dim1", ylab = "Dim2", cex.lab = 2, main = "Expression", cex.main = 2)
lines(pc_rpkm, lwd = 2)
pc_nes <- princurve::principal.curve(tsne_nes, smoother = "lowess")
plot(tsne_nes, col = rainbow(5)[as.factor(colnames(nes))], cex = 0.75, pch = 16, xlab = "Dim1", ylab = "Dim2", cex.lab = 2, main = "Activity", cex.main = 2)
lines(pc_nes, lwd = 2)
plot(pc_rpkm$lambda, pc_nes$lambda, xlab = "Expression", ylab = "Activity", cex.lab = 2, pch = 16, cex = 0.75, col = rainbow(5)[as.factor(colnames(nes))], main = "Lineage", cex.main = 2)
legend("bottomright", legend = levels(as.factor(colnames(nes))), fill = rainbow(5), cex = 1, bty = "n", border = NA)
dev.off()
#SupFig4

library(dtw)
load("./ES/GO.rda")
nes <- cbind(nes_Hs[intersect(rownames(nes_Hs), rownames(nes_Mm)), ],
             nes_Mm[intersect(rownames(nes_Hs), rownames(nes_Mm)), ])
dist <- cor(nes)
dist <- dist[1:ncol(nes_Hs), (ncol(nes_Hs)+1):ncol(nes)]
dtw <- dtw(1 - dist)
pdf("SupFig5.pdf")
par(mar = c(6, 6, 6, 2))
plot(dtw$index1, dtw$index2, axes = F, xlab = "", ylab = "", main = "Dynamic Time Warping", cex.main = 2)
axis(side = 1, tick = F, las = 2, at = 1:28, col.axis = 1,
     labels = c(rep("oocyte", 3), rep("pronuclear", 3), rep("zygote", 2),
                rep("2cell", 3), rep("4cell", 4), rep("8cell", 10), rep("morula", 3)))
axis(side = 2, tick = F, las = 2, at = 1:17, col.axis = 4,
     labels = c(rep("oocyte", 2), rep("pronuclear", 3), rep("2cell", 3), 
                rep("4cell", 3), rep("8cell", 3), rep("morula", 3)))
legend("bottomright", legend = c("Hs", "Mm"), text.col = c(1, 4), bty = "n", cex = 2)
dev.off()
#SupFig5

system("cp ./Immune/Chromium-GSE94820-QuakeMCA/hclust.pdf ./SupFig6.pdf")
#SupFig6

system("cp ./Embryo/E-MTAB-3929-GSE66688-GSE65525/E-MTAB-3929-GSE66688-GSE65525.pdf ./SupFig7.pdf")
#SupFig7

system("cp ./Benchmark/seurat.pdf ./SupFig8.pdf")
#SupFig8

system("cp ./Benchmark/GSE116272/bulk.pdf ./SupFig9.pdf")
#SupFig9

system("cp ./Benchmark/gsetFilter.pdf ./SupFig10.pdf")
#SupFig10

system("cp ./Benchmark/gsetSize.pdf ./SupFig11.pdf")
#SupFig11

system("cp ./Benchmark/normalization.pdf ./SupFig12.pdf")
#SupFig12

system("cp ./Benchmark/shadow.pdf ./SupFig13.pdf")
#SupFig13

system("cp ./Benchmark/gsetRandom.pdf ./SupFig14.pdf")
#SupFig14

system("cp ./Benchmark/variableTerms.pdf ./SupFig15.pdf")
#SupFig15

system("cp ./Immune/Chromium/stat.csv ./SupTable1.csv")
#SupTable1

system("cp ./Immune/QuakeMCA/stat.csv ./SupTable2.csv")
#SupTable2

system("cp ./Immune/GSE94820/stat.csv ./SupTable3.csv")
#SupTable3
