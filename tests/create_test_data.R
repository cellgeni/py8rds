v1 = c(1,2,3,4,5)
saveRDS(v1,"atomic_vector_num.rds")

v2 = c("a","b","c","d","e")
saveRDS(v2,"atomic_vector_char.rds")

v3 = c(TRUE,TRUE,FALSE,FALSE,TRUE)
saveRDS(v2,"atomic_vector_bool.rds")

rseq = 1:5
saveRDS(rseq,"regular_sequence.rds")

data_frame <- data.frame (
  CharVec = c("a",   "b",  "c" ),
  NumVec  = c( 1,     2,    3  ),
  BoolVec = c(TRUE, FALSE, TRUE)
)
saveRDS(data_frame,"data.frame.rds")
