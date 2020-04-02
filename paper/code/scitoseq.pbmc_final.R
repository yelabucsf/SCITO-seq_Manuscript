
## Deconvoluting multiplets using R code edited by Jimmie Ye and Byungjin Hwang
## Last modified 2020.4.1 by Byungjin Hwang (bjhwang113@gmail.com)
## let's make this code generalizable to a lot of batches
library(Seurat);
library(mixtools);
library(ggplot2);
library(ggridges)
library(parallelDist)
##library(wordspace)
library(dplyr)
library(rlang)
library(cluster)
library(fitdistrplus)
library(scater);
##library(loom)
library(Matrix);

setwd("~/Box/Papers/In Preparation/2019.Hwang.scitoseq/Code")


MaxN <- function(x, N = 2){
  len <- length(x)
  if (N > len) {
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x, partial = len - N + 1)[len - N + 1]
}

HTODemux <- function (object, assay = "HTO", positive.quantile = 0.99, init = NULL, 
                      nstarts = 100, kfunc = "clara", nsamples = 100, seed = 42, 
                      verbose = TRUE) {
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay <- assay %||% DefaultAssay(object = object)
  data <- GetAssayData(object = object, assay = assay)
  counts <- GetAssayData(object = object, assay = assay, slot = "counts")[, 
                                                                          colnames(x = object)]
  counts <- as.matrix(x = counts)
  ncenters <- init %||% (nrow(x = data) + 1)
  switch(EXPR = kfunc, kmeans = {
    init.clusters <- kmeans(x = t(x = GetAssayData(object = object, 
                                                   assay = assay)), centers = ncenters, nstart = nstarts)
    Idents(object = object, cells = names(x = init.clusters$cluster)) <- init.clusters$cluster
  }, clara = {
    init.clusters <- clara(x = t(x = GetAssayData(object = object, 
                                                  assay = assay)), k = ncenters, samples = nsamples)
    Idents(object = object, cells = names(x = init.clusters$clustering), 
           drop = TRUE) <- init.clusters$clustering
  }, stop("Unknown k-means function ", kfunc, ", please choose from 'kmeans' or 'clara'"))
  average.expression <- AverageExpression(object = object, 
                                          assays = assay, verbose = FALSE)[[assay]]
  if (sum(average.expression == 0) > 0) {
    stop("Cells with zero counts exist as a cluster.")
  }
  discrete <- GetAssayData(object = object, assay = assay)
  discrete[discrete > 0] <- 0
  for (iter in rownames(x = data)) {
    values <- counts[iter, colnames(object)]
    values.use <- values[WhichCells(object = object, idents = levels(x = Idents(object = object))[[which.min(x = average.expression[iter, 
                                                                                                                                    ])]])]
    fit <- suppressWarnings(expr = fitdist(data = values.use, 
                                           distr = "nbinom"))
    cutoff <- as.numeric(x = quantile(x = fit, probs = positive.quantile)$quantiles[1])
    discrete[iter, names(x = which(x = values > cutoff))] <- 1
    if (verbose) {
      message(paste0("Cutoff for ", iter, " : ", cutoff, 
                     " reads"))
    }
  }
  npositive <- Matrix::colSums(x = discrete)
  classification.global <- npositive
  classification.global[npositive == 0] <- "Negative"
  classification.global[npositive == 1] <- "Singlet"
  classification.global[npositive == 2] <- "Doublet"
  classification.global[npositive == 3] <- "Triplet"
  classification.global[npositive == 4] <- "Quadruplet"
  classification.global[npositive == 5] <- "Quintuplet"
  classification.global[npositive == 6] <- "Sextuplet"
  classification.global[npositive == 7] <- "Septuplet"
  donor.id = rownames(x = data)
  hash.max <- apply(X = data, MARGIN = 2, FUN = max)
  hash.maxID <- apply(X = data, MARGIN = 2, FUN = which.max)
  hash.second <- apply(X = data, MARGIN = 2, FUN = MaxN, N = 2)
  hash.maxID <- as.character(x = donor.id[sapply(X = 1:ncol(x = data), 
                                                 FUN = function(x) {
                                                   return(which(x = data[, x] == hash.max[x])[1])
                                                 })])
  hash.secondID <- as.character(x = donor.id[sapply(X = 1:ncol(x = data), 
                                                    FUN = function(x) {
                                                      return(which(x = data[, x] == hash.second[x])[1])
                                                    })])
  ##browser();
  hash.margin <- hash.max - hash.second
  doublet_id <- sapply(X = 1:length(x = hash.maxID), FUN = function(x) {
    return(paste(sort(x = c(hash.maxID[x], hash.secondID[x])), 
                 collapse = "_"))
  })
  classification <- classification.global
  classification[classification.global == "Negative"] <- "Negative"
  classification[classification.global == "Singlet"] <- hash.maxID[which(x = classification.global == 
                                                                           "Singlet")]
  classification[classification.global == "Doublet"] <- doublet_id[which(x = classification.global == 
                                                                           "Doublet")]
  classification.metadata <- data.frame(hash.maxID, hash.secondID, 
                                        hash.margin, classification, classification.global)
  colnames(x = classification.metadata) <- paste(assay, c("maxID", 
                                                          "secondID", "margin", "classification", "classification.global"), 
                                                 sep = "_")
  object <- AddMetaData(object = object, metadata = classification.metadata)
  Idents(object) <- paste0(assay, "_classification")
  doublets <- rownames(x = object[[]])[which(object[[paste0(assay, 
                                                            "_classification.global")]] == "Doublet")]
  Idents(object = object, cells = doublets) <- "Doublet"
  object$hash.ID <- Idents(object = object)
  ##object$discrete <- discrete;
  return(list(object,discrete))
}

## specify your h5 input file here ##
experiment <- "200k_PBMC_29ab"
raw <- Read10X_h5(paste("../Text/", experiment, "/filtered_feature_bc_matrix.h5",sep=""))
tx <- raw$`Gene Expression`[,which(Matrix::colSums(raw$`Gene Expression`)>0)]
#ab <- raw$`Antibody Capture`
ab <- raw$`Antibody Capture`[,which(Matrix::colSums(raw$`Gene Expression`)>0)];
## ccd = cell_containing_drops
ccd_indices <- intersect(which(Matrix::colSums(tx)>0), which(Matrix::colSums(ab)>0))
ab <- ab[,ccd_indices];
tx <- tx[,ccd_indices];
ab_colSums = Matrix::colSums(ab);
tx_colSums = Matrix::colSums(tx);

## Load vireo output here (can also be demuxlet or freemuxlet ouput)
donor <- read.table("../Text/100K_PBMC_29ab/donor_ids.tsv", header=TRUE)
#IND<-as.data.frame(donor$donor_id)
IND_filter <- donor[ccd_indices,]
IND <- as.data.frame(IND_filter$donor_id)
DBC_real <- as.data.frame(IND_filter$cell)

colnames(IND)<-"Demuxlet"
rownames(IND)<-colnames(ab)

colnames(DBC_real)<-"DBC_real"
rownames(DBC_real)<-colnames(ab)

aux <- cbind(IND,DBC_real)
#donor_ids <- ReadH5AD("/Users/yechun/Downloads/100k.h5ad", assay="ADT")
###################################################
## first detect which batch each droplet contains
## input matrix here is still under optimization, we recommend currently using colSum or (normalization to mean/median of CD45 expression) approach.
###################################################

#ab_hashtag <- rbind(Matrix::colSums(ab[1:28,]),
#                    Matrix::colSums(ab[29:56,]),
#                    Matrix::colSums(ab[57:84,]),
#                    Matrix::colSums(ab[85:112,]),
#                    Matrix::colSums(ab[113:140,]),
#                    Matrix::colSums(ab[141:168,]),
#                    Matrix::colSums(ab[169:196,]),
#                    Matrix::colSums(ab[197:224,]),
#                    Matrix::colSums(ab[225:252,]),
#                    Matrix::colSums(ab[253:280,]))

# median
#ab_hashtag <- rbind(round(colSums(t(t(ab[1:28,])/median(ab[1,]))),0),
#                    round(colSums(t(t(ab[29:56,])/median(ab[29,]))),0),
#                    round(colSums(t(t(ab[57:84,])/median(ab[57,]))),0),
#                    round(colSums(t(t(ab[85:112,])/median(ab[85,]))),0),
#                    round(colSums(t(t(ab[113:140,])/median(ab[113,]))),0),
#                    round(colSums(t(t(ab[141:168,])/median(ab[141,]))),0),
#                    round(colSums(t(t(ab[169:196,])/median(ab[169,]))),0),
#                    round(colSums(t(t(ab[197:224,])/median(ab[197,]))),0),
#                    round(colSums(t(t(ab[225:252,])/median(ab[225,]))),0),
#                    round(colSums(t(t(ab[253:280,])/median(ab[253,]))),0))
# mean no round-up
#ab_hashtag <- rbind(colSums(t(t(ab[1:28,])/mean(ab[1,]))),
#                    colSums(t(t(ab[29:56,])/mean(ab[29,]))),
#                    colSums(t(t(ab[57:84,])/mean(ab[57,]))),
#                    colSums(t(t(ab[85:112,])/mean(ab[85,]))),
#                    colSums(t(t(ab[113:140,])/mean(ab[113,]))),
#                    colSums(t(t(ab[141:168,])/mean(ab[141,]))),
#                    colSums(t(t(ab[169:196,])/mean(ab[169,]))),
#                    colSums(t(t(ab[197:224,])/mean(ab[197,]))),
#                    colSums(t(t(ab[225:252,])/mean(ab[225,]))),
#                    colSums(t(t(ab[253:280,])/mean(ab[253,]))))

#mean round-up
ab_hashtag <- rbind(round(colSums(t(t(ab[1:28,])/mean(ab[1,]))),0),
                    round(colSums(t(t(ab[29:56,])/mean(ab[29,]))),0),
                    round(colSums(t(t(ab[57:84,])/mean(ab[57,]))),0),
                    round(colSums(t(t(ab[85:112,])/mean(ab[85,]))),0),
                    round(colSums(t(t(ab[113:140,])/mean(ab[113,]))),0),
                    round(colSums(t(t(ab[141:168,])/mean(ab[141,]))),0),
                    round(colSums(t(t(ab[169:196,])/mean(ab[169,]))),0),
                    round(colSums(t(t(ab[197:224,])/mean(ab[197,]))),0),
                    round(colSums(t(t(ab[225:252,])/mean(ab[225,]))),0),
                    round(colSums(t(t(ab[253:280,])/mean(ab[253,]))),0))


#mean *100
#ab_hashtag <- rbind(round(colSums(t(t(ab[1:28,])/mean(ab[1,]))*100),0),
#                    round(colSums(t(t(ab[29:56,])/mean(ab[29,]))*100),0),
#                    round(colSums(t(t(ab[57:84,])/mean(ab[57,]))*100),0),
#                    round(colSums(t(t(ab[85:112,])/mean(ab[85,]))*100),0),
#                    round(colSums(t(t(ab[113:140,])/mean(ab[113,]))*100),0),
#                    round(colSums(t(t(ab[141:168,])/mean(ab[141,]))*100),0),
#                    round(colSums(t(t(ab[169:196,])/mean(ab[169,]))*100),0),
#                    round(colSums(t(t(ab[197:224,])/mean(ab[197,]))*100),0),
#                    round(colSums(t(t(ab[225:252,])/mean(ab[225,]))*100),0),
#                    round(colSums(t(t(ab[253:280,])/mean(ab[253,]))*100),0))

#geormetric mean
#ab_hashtag <- rbind(round(colSums(t(t(ab[1:28,])/exp(mean(log(ab[1,]+1))))),0),
#                    round(colSums(t(t(ab[29:56,])/exp(mean(log(ab[29,]+1))))),0),
#                    round(colSums(t(t(ab[57:84,])/exp(mean(log(ab[57,]+1))))),0),
#                    round(colSums(t(t(ab[85:112,])/exp(mean(log(ab[85,]+1))))),0),
#                    round(colSums(t(t(ab[113:140,])/exp(mean(log(ab[113,]+1))))),0),
#                    round(colSums(t(t(ab[141:168,])/exp(mean(log(ab[141,]+1))))),0),
#                    round(colSums(t(t(ab[169:196,])/exp(mean(log(ab[169,]+1))))),0),
#                    round(colSums(t(t(ab[197:224,])/exp(mean(log(ab[197,]+1))))),0),
#                    round(colSums(t(t(ab[225:252,])/exp(mean(log(ab[225,]+1))))),0),
#                    round(colSums(t(t(ab[253:280,])/exp(mean(log(ab[253,]+1))))),0))


rownames(ab_hashtag) = paste("Batch",seq(1:10),sep="")

# Setup Seurat object
ab_hashtag_obj <- CreateSeuratObject(counts = ab, meta.data=aux)

# Add HTO data as a new assay independent from RNA
ab_hashtag_obj[["HTO"]] <- CreateAssayObject(counts = ab_hashtag)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
ab_hashtag_obj <- NormalizeData(ab_hashtag_obj, assay = "HTO", normalization.method = "CLR")

rst <- HTODemux(ab_hashtag_obj, assay = "HTO", positive.quantile = 0.99)
ab_hashtag_obj <- rst[[1]];
discrete <- rst[[2]];

table(ab_hashtag_obj$HTO_classification.global)


############################################################
# (Optional step)
# The code in this section could be performed as a second round of classification in case you didn' classify much singlet or multiplets as expected 
# get only negative barcodes and negative HTO classification
negative_class <- ab_hashtag_obj$HTO_classification[ab_hashtag_obj$HTO_classification=='Negative']
negative_dbc <- attributes(negative_class)$names

library(tidyverse)
#IND_neg <- donor[donor$cell==negative_dbc]
IND_neg <- subset(donor, donor[,1] %in% negative_dbc)
ab_neg <-ab[,negative_dbc]

IND <- as.data.frame(IND_neg$donor_id)
DBC_real <- as.data.frame(IND_neg$cell)

colnames(IND)<-"Demuxlet"
rownames(IND)<-colnames(ab_neg)

colnames(DBC_real)<-"DBC_real"
rownames(DBC_real)<-colnames(ab_neg)

aux <- cbind(IND,DBC_real)
ab_hashtag <- rbind(colSums(t(t(ab_neg[1:28,])/mean(ab_neg[1,]))),
                    colSums(t(t(ab_neg[29:56,])/mean(ab_neg))),
                    colSums(t(t(ab_neg[57:84,])/mean(ab_neg[57,]))),
                    colSums(t(t(ab_neg[85:112,])/mean(ab_neg[85,]))),
                    colSums(t(t(ab_neg[113:140,])/mean(ab_neg[113,]))),
                    colSums(t(t(ab_neg[141:168,])/mean(ab_neg[141,]))),
                    colSums(t(t(ab_neg[169:196,])/mean(ab_neg[169,]))),
                    colSums(t(t(ab_neg[197:224,])/mean(ab_neg[197,]))),
                    colSums(t(t(ab_neg[225:252,])/mean(ab_neg[225,]))),
                    colSums(t(t(ab_neg[253:280,])/mean(ab_neg[253,]))))

rownames(ab_hashtag) = paste("Batch",seq(1:10),sep="")

# Setup Seurat object
ab_hashtag_obj <- CreateSeuratObject(counts = ab_neg, meta.data=aux)

# Add HTO data as a new assay independent from RNA
ab_hashtag_obj[["HTO"]] <- CreateAssayObject(counts = ab_hashtag)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
ab_hashtag_obj <- NormalizeData(ab_hashtag_obj, assay = "HTO", normalization.method = "CLR")

rst <- HTODemux(ab_hashtag_obj, assay = "HTO", positive.quantile = 0.9999)
ab_hashtag_obj <- rst[[1]];
discrete <- rst[[2]];

table(ab_hashtag_obj$HTO_classification.global)
singlet_rescue <- ab_hashtag_obj$HTO_classification[(ab_hashtag_obj$HTO_classification.global) =='Singlet']
singlet_dbc <- attributes(singlet_rescue)$names
multiplet_rescue <- ab_hashtag_obj$HTO_classification[ab_hashtag_obj$HTO_classification.global !='Singlet'& ab_hashtag_obj$HTO_classification.global !='Negative']
multiplet_dbc <- attributes(multiplet_rescue)$names

write.table(singlet_dbc,"100K_rescue_negative_singlet_bc.txt", col.names=c("barcodes"), sep="," )

write.table(multiplet_dbc, "100K_rescue_negative_multiplet_bc.txt", col.names=c("barcodes"), sep="," )
############################################################


Idents(ab_hashtag_obj) <- "HTO_maxID"

# check if batch classification is in the expected range
RidgePlot(ab_hashtag_obj, assay = "HTO", features = rownames(ab_hashtag_obj[["HTO"]])[1:4], ncol = 2)

#################################################################
## now only pick out the batch/droplet barcode ids that are real
#################################################################

ab_resolved <- cbind(ab[1:28,which(discrete[1,]==1)],
                     ab[29:56,which(discrete[2,]==1)],
                     ab[57:84,which(discrete[3,]==1)],
                     ab[85:112,which(discrete[4,]==1)],
                     ab[113:140,which(discrete[5,]==1)],
                     ab[141:168,which(discrete[6,]==1)],
                     ab[169:196,which(discrete[7,]==1)],
                     ab[197:224,which(discrete[8,]==1)],
                     ab[225:252,which(discrete[9,]==1)],
                     ab[253:280,which(discrete[10,]==1)])


colnames(ab_resolved) <- make.names(colnames(ab_resolved),unique=TRUE)

classification <- as.character(ab_hashtag_obj$HTO_classification.global)
DBC_real <- as.character(ab_hashtag_obj@meta.data$DBC_real)
Demuxlet <- as.character(ab_hashtag_obj@meta.data$Demuxlet)

prop.table(table(data.frame(classification, Demuxlet)),margin=2)

meta <- data.frame(Batch=c(rep("Batch1",sum(discrete[1,])),
                           rep("Batch2",sum(discrete[2,])),
                           rep("Batch3",sum(discrete[3,])),
                           rep("Batch4",sum(discrete[4,])),
                           rep("Batch5",sum(discrete[5,])),
                           rep("Batch6",sum(discrete[6,])),
                           rep("Batch7",sum(discrete[7,])),
                           rep("Batch8",sum(discrete[8,])),
                           rep("Batch9",sum(discrete[9,])),
                           rep("Batch10",sum(discrete[10,]))),
                   Type=c(classification[which(discrete[1,]==1)],
                          classification[which(discrete[2,]==1)],
                          classification[which(discrete[3,]==1)],
                          classification[which(discrete[4,]==1)],
                          classification[which(discrete[5,]==1)],
                          classification[which(discrete[6,]==1)],
                          classification[which(discrete[7,]==1)],
                          classification[which(discrete[8,]==1)],
                          classification[which(discrete[9,]==1)],
                          classification[which(discrete[10,]==1)]),
                    Demuxlet=c(Demuxlet[which(discrete[1,]==1)],
                               Demuxlet[which(discrete[2,]==1)],
                               Demuxlet[which(discrete[3,]==1)],
                               Demuxlet[which(discrete[4,]==1)],
                               Demuxlet[which(discrete[5,]==1)],
                               Demuxlet[which(discrete[6,]==1)],
                               Demuxlet[which(discrete[7,]==1)],
                               Demuxlet[which(discrete[8,]==1)],
                               Demuxlet[which(discrete[9,]==1)],
                               Demuxlet[which(discrete[10,]==1)]),
                    DBC=c(DBC_real[which(discrete[1,]==1)],
                          DBC_real[which(discrete[2,]==1)],
                          DBC_real[which(discrete[3,]==1)],
                          DBC_real[which(discrete[4,]==1)],
                          DBC_real[which(discrete[5,]==1)],
                          DBC_real[which(discrete[6,]==1)],
                          DBC_real[which(discrete[7,]==1)],
                          DBC_real[which(discrete[8,]==1)],
                          DBC_real[which(discrete[9,]==1)],
                          DBC_real[which(discrete[10,]==1)]))
                         
rownames(meta) <- colnames(ab_resolved)


ab_resolved_obj <- CreateSeuratObject(counts = ab_resolved, meta.data=meta)
ab_resolved_obj[["ADT"]] <- CreateAssayObject(counts = ab_resolved)


DefaultAssay(ab_resolved_obj) <- "ADT";

sampled <- sample(1:ncol(ab_resolved), ncol(ab_resolved))

ab_resolved_obj <- ab_resolved_obj[,sampled];

selected_ab_cell_matrix <-ab_resolved_obj[['ADT']]@counts
write.table(selected_ab_cell_matrix,"100K_selected_ab_cell_matrix_for_david_0.99_cd45mean_negative.txt", sep=",")
meta_ab_cell_matrix <-ab_resolved_obj@meta.data
write.table(meta_ab_cell_matrix,"100K_meta_ab_cell_matrix.txt_for_david_0.99_cd45mean_negative", sep=",")

## from here you can process the data further with scanpy or either Seurat (below)

## remove barcodes that are demuxlet doublets but ab singlets
ab_resolved_obj <- ab_resolved_obj[,-intersect(which(ab_resolved_obj$Demuxlet=="doublet"), which(ab_resolved_obj$Type=="Singlet"))]

ab_resolved_obj <- NormalizeData(ab_resolved_obj, assay = "ADT", normalization.method = "CLR")
ab_resolved_obj <- ScaleData(ab_resolved_obj, split.by = "Batch", assay = "ADT")##, vars.to.regress = c("Batch"))
ab_resolved_obj <- RunPCA(ab_resolved_obj, features=rownames(ab_resolved_obj))



## use top pcs
##ab_resolved_obj <- RunTSNE(ab_resolved_obj, dims=1:15, reduction="pca", assay = "ADT", reduction.key= "adtTSNE_")
ab_resolved_obj <- RunUMAP(ab_resolved_obj, dims=1:15, reduction="pca", assay = "ADT", reduction.key= "adtUMAP_")
ab_resolved_obj <- FindNeighbors(ab_resolved_obj, reduction="pca")
ab_resolved_obj <- FindClusters(ab_resolved_obj, resolution = 0.2)



plot_df <- cbind(ab_resolved_obj@reductions$umap@cell.embeddings, ab_resolved_obj@meta.data, t(GetAssayData(ab_resolved_obj,slot="data")))
colnames(plot_df) <- make.names(colnames(plot_df))

ggplot(plot_df[plot_df$Type=="Singlet",], aes(CD4.barcode1,fill=Batch))+geom_histogram()+facet_grid(Batch~.)
ggplot(plot_df, aes(adtUMAP_1, adtUMAP_2, color=Batch))+geom_point()

##plot_df_sampled <- rbind(plot_df[sample(which(plot_df$Type=="Doublet"),19082, replace=F),], plot_df[sample(which(plot_df$Type=="Singlet"),19082, replace=F),])

ggplot(plot_df, aes(adtUMAP_1, adtUMAP_2, color=seurat_clusters))+geom_point()+facet_grid(Type~.)
ggplot(plot_df[plot_df$seurat_clusters==4,], aes(adtUMAP_1, adtUMAP_2, color=seurat_clusters))+geom_point()+facet_grid(Type~.)
ggplot(plot_df[plot_df$Demuxlet=="doublet",], aes(adtUMAP_1, adtUMAP_2, color=seurat_clusters))+geom_point()+facet_grid(Type~.)

ggplot(plot_df, aes(adtUMAP_1, adtUMAP_2, color=Demuxlet))+geom_point(size=0.5)+scale_color_manual(values=c("gray","gray","red","gray"))+facet_grid(Type~.)

all_singlets_tab <- table(ab_resolved_obj$seurat_clusters[ab_resolved_obj$Type=="Singlet"]);
all_multiplets_tab <- table(ab_resolved_obj$seurat_clusters[ab_resolved_obj$Type!="Singlet"]);

all_singlets_freq <- as.vector(prop.table(table(ab_resolved_obj$seurat_clusters[ab_resolved_obj$Type=="Singlet"])))
all_multiplets_freq <- as.vector(prop.table(table(ab_resolved_obj$seurat_clusters[ab_resolved_obj$Type!="Singlet"])))

collision_cluster <- 6;
all_singlets_nc_freq <- as.vector(prop.table(all_singlets_tab[-collision_cluster]))
all_multiplets_nc_freq <- as.vector(prop.table(all_multiplets_tab[-collision_cluster]))

cor(all_singlets_freq, all_multiplets_freq)

## let's look at different cell types for their hashtags
RidgePlot(ab_hashtag_obj[,match(ab_resolved_obj$DBC[which(ab_resolved_obj$Type=="Singlet" &
                                                             ab_resolved_obj$seurat_clusters==0)], ab_hashtag_obj$DBC_real)], assay = "HTO", features = rownames(ab_hashtag_obj[["HTO"]])[1:4], ncol = 2)

RidgePlot(ab_hashtag_obj[,unique(match(ab_resolved_obj$DBC[which(ab_resolved_obj$Type!="Singlet" &
                                                            ab_resolved_obj$seurat_clusters==0)], ab_hashtag_obj$DBC_real))], assay = "HTO", features = rownames(ab_hashtag_obj[["HTO"]])[1:4], ncol = 2)

## let's look at different cell types for their hashtags
RidgePlot(ab_hashtag_obj[,match(ab_resolved_obj$DBC[which(ab_resolved_obj$Type=="Singlet" &
                                                            ab_resolved_obj$seurat_clusters==1)], ab_hashtag_obj$DBC_real)], assay = "HTO", features = rownames(ab_hashtag_obj[["HTO"]])[1:4], ncol = 2)

RidgePlot(ab_hashtag_obj[,unique(match(ab_resolved_obj$DBC[which(ab_resolved_obj$Type!="Singlet" &
                                                                   ab_resolved_obj$seurat_clusters==1)], ab_hashtag_obj$DBC_real))], assay = "HTO", features = rownames(ab_hashtag_obj[["HTO"]])[1:4], ncol = 2)

RidgePlot(ab_hashtag_obj[,match(ab_resolved_obj$DBC[which(ab_resolved_obj$Type=="Singlet" &
                                                            ab_resolved_obj$seurat_clusters==5)], ab_hashtag_obj$DBC_real)], assay = "HTO", features = rownames(ab_hashtag_obj[["HTO"]])[1:4], ncol = 2)
RidgePlot(ab_hashtag_obj[,unique(match(ab_resolved_obj$DBC[which(ab_resolved_obj$Type!="Singlet" &
                                                                   ab_resolved_obj$seurat_clusters==5)], ab_hashtag_obj$DBC_real))], assay = "HTO", features = rownames(ab_hashtag_obj[["HTO"]])[1:4], ncol = 2)

##RidgePlot(ab_resolved_obj, assay = "HTO", features = rownames(ab_resolved_obj[["HTO"]])[1:4], ncol = 2)


## find singlets
d0_demux_dbc <- IND_filter$cell[IND_filter$donor_id=="donor0"]
d1_demux_dbc <- IND_filter$cell[IND_filter$donor_id=="donor1"]


ab_resolved_filtered_obj <- ab_resolved_obj[-intersect(ab_resolved_obj$Demuxlet=="doublet", ab_resolved_obj$Type=="Singlet"),]

d0_singlets <- ab_resolved_filtered_obj[,intersect(which(!is.na(match(ab_resolved_filtered_obj$DBC, d0_demux_dbc))), which(ab_resolved_filtered_obj$Type=="Singlet"))];
d0_singlets_freq <- as.vector(table(d0_singlets$seurat_clusters)/sum(table(d0_singlets$seurat_clusters)))
d0_singlets_batch_freq <- prop.table(table(data.frame(d0_singlets$seurat_clusters, d0_singlets$Batch)),margin=2)

d1_singlets <- ab_resolved_filtered_obj[,intersect(which(!is.na(match(ab_resolved_filtered_obj$DBC, d1_demux_dbc))), which(ab_resolved_filtered_obj$Type=="Singlet"))]
d1_singlets_freq <- as.vector(table(d1_singlets$seurat_clusters)/sum(table(d1_singlets$seurat_clusters)))
d1_singlets_batch_freq <- prop.table(table(data.frame(d1_singlets$seurat_clusters, d1_singlets$Batch)),margin=2)

## find multiplets
d0_multiplets <- ab_resolved_filtered_obj[,intersect(which(!is.na(match(ab_resolved_filtered_obj$DBC, d0_demux_dbc))), which(ab_resolved_filtered_obj$Type!="Singlet"))];
d0_multiplets_freq <- as.vector(table(d0_multiplets$seurat_clusters)/sum(table(d0_multiplets$seurat_clusters)))
d0_multiplets_batch_freq <- prop.table(table(data.frame(d0_multiplets$seurat_clusters, d0_multiplets$Batch)),margin=2)

d1_multiplets <- ab_resolved_filtered_obj[,intersect(which(!is.na(match(ab_resolved_filtered_obj$DBC, d1_demux_dbc))), which(ab_resolved_filtered_obj$Type!="Singlet"))];
d1_multiplets_freq <- as.vector(table(d1_multiplets$seurat_clusters)/sum(table(d1_multiplets$seurat_clusters)))
d1_multiplets_batch_freq <- prop.table(table(data.frame(d1_multiplets$seurat_clusters, d1_multiplets$Batch)),margin=2)

cor(d0_singlets_freq, d0_multiplets_freq)
cor(d1_singlets_freq, d1_multiplets_freq)
