#!usr/bin/R
library(rbenchmark)
n_clust=NULL
n_init=100
kfunc="kmeans"
maxneighbor=100
seed=42
ncenters=11
nstarts=100

assay="HTO"
ab <- read.table("/Users/antonogorodnikov/Documents/Work/DataSci/SCITO-seq/tests/ab_raw_200K_matix.csv", sep=",", header=TRUE)
ab = subset(ab, select = -c(X) )
#colnames(ab)<-seq(1:dim(ab)[2])
rownames(ab)<-seq(1:dim(ab)[1])
ab_hashtag <- read.table("/Users/antonogorodnikov/Documents/Work/DataSci/SCITO-seq/tests/ab_merge_batch_200K_matrix.txt", sep=",", header=F)
ab_hashtag <-subset(ab_hashtag, select = -c(V1) )
colnames(ab_hashtag)<-colnames(ab)
rownames(ab_hashtag)<-paste("Batch",seq(1:10),sep="")

ab_hashtag_obj_in <- CreateSeuratObject(counts = ab)

# Add HTO data as a new assay independent from RNA
ab_hashtag_obj_in[["HTO"]] <- CreateAssayObject(counts = ab_hashtag)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation


# benchmark the original Seurat function
result <- benchmark(
  {
    ab_hashtag_obj <- NormalizeData(ab_hashtag_obj_in, assay = "HTO", normalization.method = "LogNormalize")
    rst <- HTODemux(ab_hashtag_obj, assay = "HTO", positive.quantile = 0.99,  kfunc = "kmeans")
  },
  replications = 2
)
write.table(result, file="")


  
  





