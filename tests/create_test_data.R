# Vectors ##################
v1 = c(1,2,3,4,5)
saveRDS(v1,"atomic_vector_num.rds")

v2 = c("a","b","c","d","e")
saveRDS(v2,"atomic_vector_char.rds")

v3 = c(TRUE,TRUE,FALSE,FALSE,TRUE)
saveRDS(v3,"atomic_vector_bool.rds")

# Regular sequence ##########
rseq = 1:5
saveRDS(rseq,"regular_sequence.rds")

# Data Frame ###############
data_frame <- data.frame (
  CharVec = c(NA, "a",   "b",  "c" ),
  NumVec  = c( 1,     2,    3  , NA),
  BoolVec = c(TRUE, FALSE, NA, TRUE),
  Factor  = factor(c('A',NA,'B','C')),
  row.names = c('r1','r2','r3','r4')
)
saveRDS(data_frame,"data.frame_with_rownames.rds")
rownames(data_frame) = NULL

saveRDS(data_frame,"data.frame_without_rownames.rds")

# environment
a = new.env()
a$var = 1
a$foo = 'a'
a$bar = 'b'
saveRDS(a,'environment.rds')


# Seurat  ########
# with sketch
library(Seurat)
curl::curl_download('https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz',destfile = 'pbmc3k_filtered_gene_bc_matrices.tar.gz')
system("tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz")

pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
# Initialize the Seurat object with the raw (non-normalized data).
obj <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
obj

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- SketchData(
  object = obj,
  ncells = 500,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
obj

DefaultAssay(obj) <- "sketch"
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- FindNeighbors(obj, dims = 1:50)
obj <- FindClusters(obj, resolution = 2)

obj@reductions$pca@assay.used
obj@assays$sketch
saveRDS(obj,'seu_sketch.rds')

rownames(obj@meta.data) = NULL

saveRDS(obj,'seu_sketch_no_cellnames.rds')
