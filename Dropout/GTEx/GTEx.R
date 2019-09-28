library(made4)
library(viper)
load("./dset.rda")
load("../MSigDB-regulon/GO-BP-Hs-MSigDB-regulon.rda")
len <- unlist(lapply(gset, function(x) length(x$tfmode)))
gset <- gset[len >= 50 & len <= 100]
nes <- aREA(rpkm, gset)$nes
#regular BPA

dropout <- function(x){
  x[sample(1:length(x), floor(length(x)*(1 - min(mean(x)/5, 1))))] <- 0 
  return(x)
}
rpkm_drop <- cbind(t(apply(rpkm[, 1:20], 1, dropout)), t(apply(rpkm[, 21:40], 1, dropout)))
nes_drop <- aREA(rpkm_drop, gset)$nes
#simulated drop-out BPA

save(rpkm, rpkm_drop, nes, nes_drop, file = "GO.rda")
