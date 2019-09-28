library(made4)
pdf("Fig2.pdf", width = 12, height = 8)
layout(matrix(c(1, 3, 6,
                1, 3, 6,
                1, 3, 6,
                1, 3, 6,
                1, 3, 7,
                1, 3, 7,
                2, 4, 7,
                2, 5, 7), 8, 3, byrow = T), widths = c(6, 2, 4), heights = rep(1, 8))

load("./ES/GO.rda")
nes <- cbind(nes_Hs[intersect(rownames(nes_Hs), rownames(nes_Mm)), ],
             nes_Mm[intersect(rownames(nes_Hs), rownames(nes_Mm)), ])
dist <- cor(nes)
dist[1:ncol(nes_Hs), (ncol(nes_Hs)+1):ncol(nes)] <- 0
for (i in 1:ncol(nes_Hs)){
  for (j in 1:ncol(nes_Hs)){
    if (i > j) dist[i, j] <- 0}}
for (i in (ncol(nes_Hs)+1):ncol(nes)){
  for (j in (ncol(nes_Hs)+1):ncol(nes)){
    if (i > j) dist[i, j] <- 0}}
#r(Pearson)
bp <- nes[c("GO_NEGATIVE_REGULATION_OF_INTRACELLULAR_PROTEIN_TRANSPORT",
            "GO_ACTIVATION_OF_GTPASE_ACTIVITY",
            "GO_POSITIVE_REGULATION_OF_TRANSCRIPTION_FACTOR_IMPORT_INTO_NUCLEUS",
            "GO_RNA_3_END_PROCESSING",
            "GO_RIBOSOME_ASSEMBLY",
            "GO_REGULATION_OF_TRANSLATIONAL_INITIATION",
            "GO_REGULATION_OF_PROTEIN_TARGETING_TO_MITOCHONDRION"), ]
bp <- (bp - rowMeans(bp))/apply(bp, 1, sd)
#BP-Activity
par(mar = c(5, 5, 5, 5))
image(abs(dist)^0.5*sign(dist), col = colorRampPalette(c("Blue", "White", "Red"))(10000), breaks = seq(-1, 1, length.out = 10001), axes = F)
for (i in c(2, 4)){
  axis(side = i, line = -1, at = 0:2/44, tick = F, labels = rep("oocyte", 3), las = 2, cex = 0.5, col.axis = getcol(21)[1])
  axis(side = i, line = -1, at = 3:5/44, tick = F, labels = rep("pronuclear", 3), las = 2, cex = 0.5, col.axis = getcol(21)[2])
  axis(side = i, line = -1, at = 6:7/44, tick = F, labels = rep("zygote", 2), las = 2, cex = 0.5, col.axis = getcol(21)[3])
  axis(side = i, line = -1, at = 8:10/44, tick = F, labels = rep("2cell", 3), las = 2, cex = 0.5, col.axis = getcol(21)[4])
  axis(side = i, line = -1, at = 11:14/44, tick = F, labels = rep("4cell", 4), las = 2, cex = 0.5, col.axis = getcol(21)[5])
  axis(side = i, line = -1, at = 15:24/44, tick = F, labels = rep("8cell", 10), las = 2, cex = 0.5, col.axis = getcol(21)[6])
  axis(side = i, line = -1, at = 25:27/44, tick = F, labels = rep("morula", 3), las = 2, cex = 0.5, col.axis = getcol(21)[7])
}
for (i in c(1, 3)){
  axis(side = i, line = -1, at = 28:29/44, tick = F, labels = rep("oocyte", 2), las = 2, cex = 0.5, col.axis = getcol(21)[1])
  axis(side = i, line = -1, at = 30:32/44, tick = F, labels = rep("pronuclear", 3), las = 2, cex = 0.5, col.axis = getcol(21)[2])
  axis(side = i, line = -1, at = 33:35/44, tick = F, labels = rep("2cell", 3), las = 2, cex = 0.5, col.axis = getcol(21)[4])
  axis(side = i, line = -1, at = 36:38/44, tick = F, labels = rep("4cell", 3), las = 2, cex = 0.5, col.axis = getcol(21)[5])
  axis(side = i, line = -1, at = 39:41/44, tick = F, labels = rep("8cell", 3), las = 2, cex = 0.5, col.axis = getcol(21)[6])
  axis(side = i, line = -1, at = 42:44/44, tick = F, labels = rep("morula", 3), las = 2, cex = 0.5, col.axis = getcol(21)[7])
}
axis(side = 2, line = -26, at = 14/44, tick = F, labels = "Hs", cex.axis = 4)
axis(side = 3, line = -14, at = 36/44, tick = F, labels = "Mm", cex.axis = 4)
legend("topleft", legend = c("oocyte", "pronuclear", "zygote", "2cell", "4cell", "8cell", "morula"), fill = getcol(21), bty = "n", border = NA, cex = 1.5)
#r(Pearson)
par(mar = c(2, 5, 0, 5))
bp1 <- bp
bp1[bp1 > 2] <- 2
bp1[bp1 < -2] <- -2
bp1[, 1:ncol(nes_Hs)] <- 0
image(t(bp1), col = colorRampPalette(c("Purple4", "White", "Orange2"))(10000), breaks = seq(-2, 2, length.out = 10001), axes = F)
axis(side = 2, line = -22, at = 0:6/6, labels = c("ProteinTransport", "GTPaseSignaling",
                                                  "TranscriptionRegulation", "RNAProcessing", "RibosomeBiogenesis",
                                                  "Translation", "Mitochondrion"), las = 2, cex.axis = 2, col = 1, tick = F)
#Mm activity
par(mar = c(5, 0, 5, 2))
bp1 <- bp
bp1[, (ncol(nes_Hs)+1):ncol(nes)] <- 0
bp1[bp1 > 2] <- 2
bp1[bp1 < -2] <- -2
image(bp1, col = colorRampPalette(c("Purple4", "White", "Orange2"))(10000), breaks = seq(-2, 2, length.out = 10001), axes = F)
axis(side = 3, line = -14, at = 0:6/6, labels = c("ProteinTransport", "GTPaseSignaling",
                                                  "TranscriptionRegulation", "RNAProcessing", "RibosomeBiogenesis",
                                                  "Translation", "Mitochondrion"), las = 2, cex.axis = 2, col = 1, tick = F)
#Hs activity
par(mar = c(3, 1, 3, 1))
image(matrix(1:20, 20, 1, byrow = T), col = colorRampPalette(c("Blue", "White", "Red"))(20), axes = F, main = "r(Pearson)", cex.main = 2)
axis(side = 1, line = 0, at = c(0, 0.5, 1), labels = c(-1, 0, 1), tick = F, cex.axis = 2)
image(matrix(1:20, 20, 1, byrow = T), col = colorRampPalette(c("Purple4", "White", "Orange2"))(20), axes = F, main = "Z-Score", cex.main = 2)
axis(side = 1, line = 0, at = c(0, 0.5, 1), labels = c(-2, 0, 2), tick = F, cex.axis = 2)
#scale
#GSE44183

load("./Immune/Chromium/GO.rda")
anno_1 <- c("Out", "B", "T", "Mono")[dbscan$cluster+1]
load("./Immune/QuakeMCA/GO.rda")
anno_2 <- c("Out", "B", "T", rep("Out", 5))[dbscan$cluster+1]
load("./Immune/GSE94820/GO.rda")
anno_3 <- sapply(strsplit(colnames(nes), split = "_"), function(x) x[1])
load("./Immune/Chromium-GSE94820-QuakeMCA/GO.rda")
par(mar = c(3, 3, 4, 1))
plot(mds, col = c(rep("darkred", length(anno_1)), rep("darksalmon", length(anno_2)), rep("darkorchid", length(anno_3))),
     pch = 16, cex = 0.75, xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5), xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "Datasets", cex.main = 2)
axis(side = 1, line = 0, at = 0, tick = F, labels = "Dim1", cex.axis = 2)
axis(side = 2, line = 0, at = 0, tick = F, labels = "Dim2", cex.axis = 2)
legend("topright", legend = c("Chromium (Hs)", "Tabula Muris (Mm)", "GSE94820 (Hs)"), fill = c("darkred", "darksalmon", "darkorchid"), cex = 1.5, bty = "n", border = NA)
plot(mds, col = c(c(4, 6, 1, 5)[as.factor(anno_1)], c(4, 1, 5)[as.factor(anno_2)], c(7, 7, 7, 6, 8)[as.factor(anno_3)]),
     pch = 16, cex = 0.75, xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5), xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "Cell Type", cex.main = 2)
axis(side = 1, line = 0, at = 0, tick = F, labels = "Dim1", cex.axis = 2)
axis(side = 2, line = 0, at = 0, tick = F, labels = "Dim2", cex.axis = 2)
legend("topright", legend = c("Outlier", "B", "T", "Mono", "cDC", "pDC"), fill = c(1, 4, 5, 6, 7, 8), cex = 1.5, bty = "n", border = NA)
#Immuno
dev.off()
