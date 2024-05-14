show_rds <- function(object, header = TRUE,print=TRUE) {
  rds <- object |>
    serialize(connection = NULL, ascii = TRUE) |>
    rawToChar() |>
    strsplit(split = "\n") |>
    unlist()
  if(header == FALSE) rds <- rds[-(1:6)]
  if(print)
    cat(rds, sep = "\n")
  invisible(rds)
}

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


# create more exotic types ######
# environments, parsing is not implemented yet
a = new.env()
a$var = 1
a$foo = 'a'
a$bar = 'b'
# show_rds(a) # strange irregular stretches of 254....
saveRDS(a,'env.rds')


# pairlists are used to store attributes and expressions
pl = pairlist(name1='val1',name2='val1',name2='val1')
# the names of pairlists are added to seqpool and all next occurences are recorded as index=(i+1)/2^8 - 2 in the seqpool
# show_rds(pl)
saveRDS(pl,'pairlist.rds')

# symbols looks similar to chars
sym = rlang::sym('abc')
# show_rds(sym)
saveRDS(sym,'symbol.rds')

# but contrary to chars symbols also stored in seqpools
# show_rds(list('abc','abc',rlang::sym('abc'),rlang::sym('abc'),rlang::sym('abc1')))

# which joined with pairlists makes things difficult
# show_rds(pairlist(name1='val1',name2='val1',name2=rlang::sym('abc'),name2=rlang::sym('abc')))