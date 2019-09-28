load("../E-MTAB-3929/GO.rda")
nes_1 <- nes
load("../GSE66688/GO.rda")
nes_2 <- nes
load("../GSE65525/GO.rda")
nes_3 <- nes
#data preparation

library(viper)
library(Rtsne)
library(made4)
go <- Reduce(intersect, list(rownames(nes_1), rownames(nes_2), rownames(nes_3)))
nes <- cbind(nes_1[go, ], nes_2[go, ], nes_3[go, ])
tsne <- Rtsne(as.dist(viperSimilarity(nes)), is_distance = T, pca = F, max_iter = 500, perplexity = 50)$Y
anno <- c(colnames(nes_1), rep("Dr", ncol(nes_2)), rep("Mm", ncol(nes_3)))
pdf("E-MTAB-3929-GSE66688-GSE65525.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))
plot(tsne, col = c(rep(1, ncol(nes_1)), rep(4, ncol(nes_2)), rep(5, ncol(nes_3))), pch = 16, cex = 0.5, xlab = "Dim1", ylab = "Dim2", main = "GO-BP-GSEA")
legend("topleft", legend = c("Hs", "Dr", "Mm"), fill = c(1, 4, 5),  bty = "n", border = NA)
plot(tsne[1:ncol(nes_1), ], col = getcol(21)[as.factor(colnames(nes_1))], xlim = range(tsne[, 1]), ylim = range(tsne[, 2]), pch = 16, cex = 0.5, xlab = "Dim1", ylab = "Dim2", main = "E-MTAB-3929")
legend("topleft", legend = levels(as.factor(colnames(nes_1))), fill = getcol(21),  bty = "n", border = NA)
plot(tsne[(ncol(nes_1)+1):(ncol(nes_1)+ncol(nes_2)), ], col = 1, xlim = range(tsne[, 1]), ylim = range(tsne[, 2]), pch = 16, cex = 0.5, xlab = "Dim1", ylab = "Dim2", main = "GSE66688")
plot(tsne[(ncol(nes_1)+ncol(nes_2)+1):(ncol(nes_1)+ncol(nes_2)+ncol(nes_3)), ], col = 1, xlim = range(tsne[, 1]), ylim = range(tsne[, 2]), pch = 16, cex = 0.5, xlab = "Dim1", ylab = "Dim2", main = "GSE65525")
dev.off()

save(nes, anno, tsne, file = "GO.rda")
