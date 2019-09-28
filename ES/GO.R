library(made4)
load("GSE44183.rda")
Hs <- log2(Hs+1)
Mm <- log2(Mm+1)
#data preparation

library(viper)
load("../MSigDB-regulon/GO-BP-Hs-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset <- gset[len >= 50 & len <= 100]
rank <- apply(Hs, 2, rank)
rank <- (rank - apply(rank, 1, median))/apply(rank, 1, mad)
nes_Hs <- aREA(rank, gset)$nes
load("../MSigDB-regulon/GO-BP-Mm-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset <- gset[len >= 50 & len <= 100]
rank <- apply(Mm, 2, rank)
rank <- (rank - apply(rank, 1, median))/apply(rank, 1, mad)
nes_Mm <- aREA(rank, gset)$nes
save(nes_Hs, nes_Mm, Hs, Mm, file = "GO.rda")
#BPA
