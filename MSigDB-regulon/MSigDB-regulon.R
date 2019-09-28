library(msigdbr)

table <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")#GO, BP
len <- table(table$gs_name)
gset <- lapply(names(len), function(x, table){
  gene <- table$gene_symbol[grep(x, table$gs_name)]
  list(tfmode=structure(rep(1, length(gene)), names=gene), likelihood=rep(1, length(gene)))
}, table=table)
names(gset) <- names(len)
save(gset, file = "GO-BP-Hs-MSigDB-regulon.rda")
#generate GO-BP-Hs-MSigDB regulon

table <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")#GO, BP
len <- table(table$gs_name)
gset <- lapply(names(len), function(x, table){
  gene <- table$gene_symbol[grep(x, table$gs_name)]
  list(tfmode=structure(rep(1, length(gene)), names=gene), likelihood=rep(1, length(gene)))
}, table=table)
names(gset) <- names(len)
save(gset, file = "GO-BP-Mm-MSigDB-regulon.rda")
#generate GO-BP-Mm-MSigDB regulon

table <- msigdbr(species = "Danio rerio", category = "C5", subcategory = "BP")#GO, BP
len <- table(table$gs_name)
gset <- lapply(names(len), function(x, table){
  gene <- table$gene_symbol[grep(x, table$gs_name)]
  list(tfmode=structure(rep(1, length(gene)), names=gene), likelihood=rep(1, length(gene)))
}, table=table)
names(gset) <- names(len)
save(gset, file = "GO-BP-Dr-MSigDB-regulon.rda")
#generate GO-BP-Dr-MSigDB regulon

table <- msigdbr(species = "Homo sapiens", category = "C7")#Immuno
len <- table(table$gs_name)
gset <- lapply(names(len), function(x, table){
  gene <- table$gene_symbol[grep(x, table$gs_name)]
  list(tfmode=structure(rep(1, length(gene)), names=gene), likelihood=rep(1, length(gene)))
}, table=table)
names(gset) <- names(len)
save(gset, file = "Immuno-Hs-MSigDB-regulon.rda")
#generate Immuno-Hs-MSigDB regulon

table <- msigdbr(species = "Mus musculus", category = "C7")#Immuno
len <- table(table$gs_name)
gset <- lapply(names(len), function(x, table){
  gene <- table$gene_symbol[grep(x, table$gs_name)]
  list(tfmode=structure(rep(1, length(gene)), names=gene), likelihood=rep(1, length(gene)))
}, table=table)
names(gset) <- names(len)
save(gset, file = "Immuno-Mm-MSigDB-regulon.rda")
#generate Immuno-Mm-MSigDB regulon


