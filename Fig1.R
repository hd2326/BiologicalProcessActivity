library(princurve)
pdf("Fig1.pdf", width = 7, height = 6)
layout(matrix(c(1, 4, 5,
                2, 4, 5,
                3, 4, 5,
                6, 9, 10,
                7, 9, 10,
                8, 9, 10), 6, 3, byrow = T), widths = c(2, 6, 6), heights = rep(2, 6))

par(mar = c(3, 1, 3, 1))
image(matrix(1:20, 20, 1, byrow = T), col = colorRampPalette(c("Black", "Red"))(20), axes = F, main = "CD3E, T")
axis(side = 1, line = 0, at = c(0, 1), labels = c("-", "+"), tick = F, cex.axis = 2)
image(matrix(1:20, 20, 1, byrow = T), col = colorRampPalette(c("Black", "Green"))(20), axes = F, main = "CD14, Mono")
axis(side = 1, line = 0, at = c(0, 1), labels = c("-", "+"), tick = F, cex.axis = 2)
image(matrix(1:20, 20, 1, byrow = T), col = colorRampPalette(c("Black", "Blue"))(20), axes = F, main = "CD20, B")
axis(side = 1, line = 0, at = c(0, 1), labels = c("-", "+"), tick = F, cex.axis = 2)

par(mar = c(5, 5, 3, 3))
load("./Immune/Chromium/GO.rda")
load("./Immune/Chromium/Chromium.rda")
col <- tpm[c("CD3E", "CD14", "MS4A1"), ]
col <- (col - apply(col, 1, min))/apply(col, 1, sd)
col <- apply(col, 2, prop.table)
col[!is.finite(col)] <- 0
col <- rgb(col[1, ], col[2, ], col[3, ])
plot(mds_tpm, col = col, pch = c(rep(3, 500), rep(1, 500)), cex = 0.5, xlab = "Dim1", ylab = "Dim2", main = "Expression", cex.main = 2)
plot(mds_nes, col = col, pch = c(rep(3, 500), rep(1, 500)), cex = 0.5, xlab = "Dim1", ylab = "Dim2", main = "Cell Type", cex.main = 2)
legend("topright", legend = c("V1", "V2"), pch = c(3, 1), cex = 1, bty = "n", border = NA)

par(mar = c(3, 1, 3, 1))
plot.new()
image(matrix(1:20, 20, 1, byrow = T), col = colorRampPalette(c("White", "Yellow", "Red"))(20), axes = F, main = "r(Pearson)")
axis(side = 1, line = 0, at = c(0, 0.5, 1), labels = c("0.7", "0.85", "1"), tick = F, cex.axis = 1)
plot.new()

par(mar = c(6, 6, 2, 2))
load("./Dropout/GTEx/GO.rda")
r1 <- t(cor(cbind(rpkm, rpkm_drop)))
r1[upper.tri(r1)] <- 0
r2 <- t(cor(cbind(nes, nes_drop)))
r2[lower.tri(r2)] <- 0
image(t(r1 + r2), col = colorRampPalette(c("White", "Yellow", "Red"))(20), breaks = seq(0.7, 1, length.out = 21), axes = F, main = "", cex.main = 2)
axis(side = 1, line = -1, at = seq(0.125, 0.875, length.out = 4), labels = c("E-Ori", "L-Ori", "E-Dro", "L-Dro"), tick = F, las = 2, cex.axis = 2)
axis(side = 2, line = -1, at = seq(0.125, 0.875, length.out = 4), labels = c("E-Ori", "L-Ori", "E-Dro", "L-Dro"), tick = F, las = 2, cex.axis = 2)
axis(side = 3, line = 0, at = 0.5, labels = c("Expression"), tick = F, cex.axis = 2)
axis(side = 4, line = 0, at = 0.5, labels = c("Activity"), tick = F, cex.axis = 2)

par(mar = c(5, 5, 3, 3))
load("./Embryo/E-MTAB-3929/GO.rda")
pc <- princurve::principal.curve(tsne_nes, smoother = "lowess")
plot(tsne_nes, col = rainbow(5)[as.factor(colnames(nes))], cex = 0.75, pch = 16, xlab = "Dim1", ylab = "Dim2", main = "Cell Transition", cex.main = 2)
lines(pc, lwd = 2)
legend("topright", legend = levels(as.factor(colnames(nes))), fill = rainbow(5), cex = 1, bty = "n", border = NA)
legend("bottomright", legend = c(paste("r(Lin,GATA3)=", signif(cor(pc$lambda, rpkm["GATA3", ]), 2), sep = ""),
                                 paste("r(Lin,DPPA5)=", signif(cor(pc$lambda, rpkm["DPPA5", ]), 2), sep = "")), bty = "n", border = NA)
dev.off()
