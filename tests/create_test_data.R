### Vectors
v1 = c(1,2,3,4,5)
saveRDS(v1,"atomic_vector_num.rds")

v2 = c("a","b","c","d","e")
saveRDS(v2,"atomic_vector_char.rds")

v3 = c(TRUE,TRUE,FALSE,FALSE,TRUE)
saveRDS(v3,"atomic_vector_bool.rds")

### Regular sequence
rseq = 1:5
saveRDS(rseq,"regular_sequence.rds")

### Data Frame
data_frame <- data.frame (
  CharVec = c("a",   "b",  "c" ),
  NumVec  = c( 1,     2,    3  ),
  BoolVec = c(TRUE, FALSE, TRUE)
)
saveRDS(data_frame,"data.frame.rds")

### Seurat (S4)
download.file("https://cf.10xgenomics.com/samples/cell-vdj/4.0.0/sc5p_v2_hs_PBMC_1k/sc5p_v2_hs_PBMC_1k_filtered_feature_bc_matrix.tar.gz","PBMC.tar.gz")
untar("PBMC.tar.gz", exdir="10xPBMC")
library(Seurat)
pbmc = Read10X("10xPBMC/filtered_feature_bc_matrix")
seurat_object = CreateSeuratObject(counts = pbmc$`Gene Expression`)
saveRDS(seurat_object, "seurat.rds")
