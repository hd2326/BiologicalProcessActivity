load("rpkm.rda")
rpkm <- rpkm[apply(rpkm, 1, sd) > 0, ]
#data preparation

library(viper)
load("/private/groups/stuartlab/hding16/GO/MSigDB-regulon/GO-BP-Hs-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset <- gset[len >= 50 & len <= 100]
nes <- aREA(rpkm, gset)$nes
#BPA

library(Rtsne)
tsne_rpkm <- Rtsne(as.dist(1 - cor(rpkm)), is_distance = T, pca = F, max_iter = 500, perplexity = 150)$Y
tsne_nes <- Rtsne(as.dist(viperSimilarity(nes)), is_distance = T, pca = F, max_iter = 500, perplexity = 150)$Y
#dimension reduction

save(tsne_rpkm, tsne_nes, rpkm, nes, file = "GO.rda")
