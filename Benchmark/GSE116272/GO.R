file <- list.files(pattern = ".counts.txt")[1:3]
gene <- unique(unlist(lapply(file, function(f){
  counts <- read.delim(f)
  as.character(counts$gene)
})))
rpm_Mm <- matrix(0, length(gene), length(file), dimnames = list(gene, file))
for (f in file){
  message(f)
  counts <- read.delim(f)
  counts <- counts[!duplicated(counts$gene), ]
  rpm_Mm[as.character(counts$gene), f] <- log2(counts[, 2]/(sum(counts[, 2])/1e6) + 1)}

file <- list.files(pattern = ".counts.txt")[8:11]
gene <- unique(unlist(lapply(file, function(f){
  counts <- read.delim(f)
  as.character(counts$gene)
})))
rpm_Hs <- matrix(0, length(gene), length(file), dimnames = list(gene, file))
for (f in file){
  message(f)
  counts <- read.delim(f)
  counts <- counts[!duplicated(counts$gene), ]
  rpm_Hs[as.character(counts$gene), f] <- log2(counts[, 2]/(sum(counts[, 2])/1e6) + 1)}

library(viper)
load("../../Dropout/GTEx/GO.rda")

load("../../MSigDB-regulon/GO-BP-Hs-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset <- gset[len >= 50 & len <= 100]
nes_Hs <- aREA(rpm_Hs, gset)$nes

load("/../../MSigDB-regulon/GO-BP-Mm-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset <- gset[len >= 50 & len <= 100]
nes_Mm <- aREA(rpm_Mm, gset)$nes

nes <- cbind(nes_Hs[Reduce(intersect, list(rownames(nes_Hs), rownames(nes_Mm), rownames(nes))), ],
             nes_Mm[Reduce(intersect, list(rownames(nes_Hs), rownames(nes_Mm), rownames(nes))), ],
             nes[Reduce(intersect, list(rownames(nes_Hs), rownames(nes_Mm), rownames(nes))), 21:40])

pdf("bulk.pdf", width = 5, height = 5)
par(mar = c(6, 6, 1, 1))
image(cor(nes), col = cm.colors(100), axes = F)
axis(side = 1, at = seq(0, 1, length.out = ncol(nes)), labels = c(rep("Eso_Epi_Mm", 3), rep("Eso_Epi_Hs", 4), rep("Lung_Whole_Hs", 20)), tick = F, las = 2, cex.axis = 0.5)
axis(side = 2, at = seq(0, 1, length.out = ncol(nes)), labels = c(rep("Eso_Epi_Mm", 3), rep("Eso_Epi_Hs", 4), rep("Lung_Whole_Hs", 20)), tick = F, las = 2, cex.axis = 0.5)
dev.off()
