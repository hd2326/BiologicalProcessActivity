load("../Chromium/GO.rda")
nes_1 <- nes
anno_1 <- c("Out", "B", "T", "Mono")[dbscan$cluster+1]
load("../QuakeMCA/GO.rda")
nes_2 <- nes
anno_2 <- c("Out", "B", "T", rep("Out", 5))[dbscan$cluster+1]
load("../GSE94820/GO.rda")
nes_3 <- nes
anno_3 <- sapply(strsplit(colnames(nes_3), split = "_"), function(x) x[1])
load("../Chromium/Chromium.rda")
c1 <- get(load("../QuakeMCA/Spleen_cpm.rda"))
colnames(c1) <- paste("Spleen", 1:ncol(c1), sep = "_")
c2 <- get(load("../QuakeMCA/Thymus_cpm.rda"))
colnames(c2) <- paste("Thymus", 1:ncol(c2), sep = "_")
cpm <- cbind(c1, c2)
cpm <- cpm[apply(cpm, 1, sd) > 0, ]
load("../GSE94820/GSE94820.rda")
expmat <- expmat[, colSums(expmat > 0) > 2000 & colSums(expmat > 0) < 8000]
#data preparation

library(viper)
library(made4)
go <- Reduce(intersect, list(rownames(nes_1), rownames(nes_2), rownames(nes_3)))
nes <- cbind(nes_1[go, ], nes_2[go, ], nes_3[go, ])
mds <- cmdscale(as.dist(viperSimilarity(nes[apply(nes, 1, sd) > 1.5, ])))
save(mds, file = "GO.rda")
#BPA-level dimension reduction

library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
table <- getBM(attributes = c("mmusculus_homolog_associated_gene_name", "external_gene_name"),
               filters = "external_gene_name",
               values = intersect(rownames(tpm), rownames(expmat)),
               mart = ensembl)
table <- table[match(intersect(rownames(cpm), table$mmusculus_homolog_associated_gene_name), table$mmusculus_homolog_associated_gene_name), ]
exp <- cbind(tpm[table$external_gene_name, ], cpm[table$mmusculus_homolog_associated_gene_name, ], expmat[table$external_gene_name, ])
#ortholog gene mapping

library(viper)
library(dendextend)
hc1 <- hclust(as.dist(viperSimilarity(nes[apply(nes, 1, sd) > 1.5, ])), method = "ward.D2")
dend1 <- as.dendrogram(hc1)
labels(dend1) <- ""
hc2 <- hclust(as.dist(1 - cor(exp[apply(exp, 1, sd) > 0, ])), method = "ward.D2")
dend2 <- as.dendrogram(hc2)
labels(dend2) <- ""
pdf("hclust.pdf", width = 10, height = 7)
layout(matrix(1:6, 3, 2, byrow = F), heights = c(5, 1, 1))
par(mar = c(1, 5, 5, 5))
plot(dend1, main = "Activity")
legend("topright", legend = c("Chromium", "QuakeMCA", "GSE94820"), fill = c(1, 4, 5),  bty = "n", border = NA)
legend("topleft", legend = c("Outlier", "B", "T", "Mono", "DC"), fill = c(1, 4, 6, 5, 8),  bty = "n", border = NA)
par(mar = c(2, 5, 2, 5))
barplot(rep(1, ncol(nes)), col = c(rep(1, ncol(nes_1)),
                                   rep(4, ncol(nes_2)),
                                   rep(5, ncol(nes_3)))[hc1$order], space = 0, border = NA, axes = F, main = "Dataset")
barplot(rep(1, ncol(nes)), col = c(c(4, 5, 1, 6)[as.factor(anno_1)],
                                   c(4, 1, 6)[as.factor(anno_2)],
                                   c(8, 5)[grepl("Mono", anno_3)+1])[hc1$order], space = 0, border = NA, axes = F, main = "Cell")

par(mar = c(1, 5, 5, 5))
plot(dend2, main = "Expression")
par(mar = c(2, 5, 2, 5))
barplot(rep(1, ncol(exp)), col = c(rep(1, ncol(tpm)),
                                   rep(4, ncol(cpm)),
                                   rep(5, ncol(expmat)))[hc2$order], space = 0, border = NA, axes = F, main = "Dataset")
barplot(rep(1, ncol(exp)), col = c(c(4, 5, 1, 6)[as.factor(anno_1)],
                                   c(4, 1, 6)[as.factor(anno_2)],
                                   c(8, 5)[grepl("Mono", anno_3)+1])[hc2$order], space = 0, border = NA, axes = F, main = "Cell")
dev.off()
#dendrogram

cell <- c(anno_1, anno_2, anno_3)
cell[cell == "CD141" | cell == "CD1C" | cell == "DoubleNeg" | cell == "pDC" | cell == "Mono"] <- "Mono-DC"
batch <- c(rep("B1", length(anno_1)), rep("B2", length(anno_2)), rep("B3", length(anno_3)))
d1 <- structure(lapply(unique(cell), function(c, cell, nes){
  dist <- viperSimilarity(cbind(nes[, cell == c], mean=rowMeans(nes[, cell == c])))
  dist[, "mean"]
}, cell=cell, nes=nes), names=unique(cell))
d2 <- structure(lapply(unique(batch), function(b, batch, cell, nes){
  dist <- viperSimilarity(cbind(nes[, batch == b & cell != "Out"], mean=rowMeans(nes[, batch == b & cell != "Out"])))
  dist[, "mean"]
}, batch=batch, cell=cell, nes=nes), names=unique(batch))
# > var.test(unlist(d1[1:3]), unlist(d2))
# 
# F test to compare two variances
# 
# data:  unlist(d1[1:3]) and unlist(d2)
# F = 0.81256, num df = 4437, denom df = 4437, p-value = 4.974e-12
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   0.7661200 0.8618234
# sample estimates:
#   ratio of variances 
# 0.812564 

d1 <- structure(lapply(unique(cell), function(c, cell, exp){
  dist <- cbind(exp[, cell == c], mean=rowMeans(exp[, cell == c]))
  dist <- 1 - cor(dist[apply(dist, 1, sd) > 0, ])
  dist[, "mean"]
}, cell=cell, exp=exp), names=unique(cell))
d2 <- structure(lapply(unique(batch), function(b, batch, cell, exp){
  dist <- cbind(exp[, batch == b & cell != "Out"], mean=rowMeans(exp[, batch == b & cell != "Out"]))
  dist <- 1 - cor(dist[apply(dist, 1, sd) > 0, ])
  dist[, "mean"]
}, batch=batch, cell=cell, exp=exp), names=unique(batch))
# > var.test(unlist(d1[1:3]), unlist(d2))
# 
# F test to compare two variances
# 
# data:  unlist(d1[1:3]) and unlist(d2)
# F = 1.0913, num df = 4437, denom df = 4437, p-value = 0.003639
# alternative hypothesis: true ratio of variances is not equal to 1
# 95 percent confidence interval:
#   1.028881 1.157409
# sample estimates:
#   ratio of variances 
# 1.091254 
#F-test
