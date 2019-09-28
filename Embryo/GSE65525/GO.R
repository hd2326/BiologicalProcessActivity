load("tpm.rda")
tpm <- tpm[apply(tpm, 1, sd) > 0, ]
#data preparation

library(viper)
load("/private/groups/stuartlab/hding16/GO/MSigDB-regulon/GO-BP-Dr-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset <- gset[len >= 50 & len <= 100]
nes <- aREA(tpm, gset)$nes
#BPA

library(Rtsne)
tsne_tpm <- Rtsne(as.dist(1 - cor(tpm)), is_distance = T, pca = F, max_iter = 500, perplexity = 150)$Y
tsne_nes <- Rtsne(as.dist(viperSimilarity(nes)), is_distance = T, pca = F, max_iter = 500, perplexity = 150)$Y
#dimension reduction

save(tsne_tpm, tsne_nes, tpm, nes, file = "GO.rda")
