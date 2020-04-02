library(Matrix);
library(dplyr)
library(data.table)
library(rlang)
library(MASS)
library(fitdistrplus)

drop_id_to_type <- function(drop_id, npools) {
  drop_type = rep("singlet",length(drop_id))
  drop_type[which(drop_id=="")] <- "unstained"
  drop_type[grep("_",drop_id)] <- "multiplet_cross_ct"
  for(i in grep("_",drop_id)) {
    id = strsplit(drop_id[i],"_")[[1]];
    ##print(id)
    if((length(na.omit(match(1:npools, id)))==length(id)) || (length(na.omit(match((npools+1):(2*npools), id)))==length(id))) {
         ##print("here i am")
         drop_type[i] <- "multiplet_same_ct";
    }
  }
  
  ##length(ab_droplet_type[ab_droplet_type == "doublet_cross_ct"])
  ##length(ab_droplet_type[ab_droplet_type == "doublet_same_ct"])
  
  return(drop_type);
}

MaxN <- function(x, N = 2){
  len <- length(x)
  if (N > len) {
    warning('N greater than length(x).  Setting N=length(x)')
    N <- length(x)
  }
  sort(x, partial = len - N + 1)[len - N + 1]
}

HTODemux <- function (object, assay = "HTO", positive.quantile = 0.90, collision.quantile = 0.99, init = NULL, 
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
  classification.global[npositive > 1] <- "Doublet"
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

htodemux_threshold <- function(ab, threshold=0.95) {
  npools <- nrow(ab)/2;
  
  ab_obj <- CreateSeuratObject(counts = ab)
  ab_obj[["HTO"]] <- CreateAssayObject(counts = ab)
  
  ab_obj <- NormalizeData(ab_obj, assay = "HTO", normalization.method = "LogNormalize")
  
  rst <- HTODemux(ab_obj, assay = "HTO", kfunc="clara", positive.quantile = threshold)
  ab_obj <- rst[[1]];
  discrete <- rst[[2]];
  
  ##  discrete2 <- fread("~/Downloads/TB_200K_discrete_values_data_frame_raw.csv",header=T)
  
  table(ab_obj$HTO_classification.global)
  
  ab_resolved_hto <- NULL;
  
  for(i in 1:npools) {
    ct1_index <- i;
    ct2_index <- i+npools;
    ab_resolved_hto <- cbind(ab_resolved_hto,discrete[c(ct1_index, ct2_index),]);
  }
  
  cell_type <- apply(ab_resolved_hto,2,paste,collapse="_")
  
  drop_id <- apply(discrete, 2, function(x) {
    paste(which(x==1),collapse="_")
  })
  
  drop_type <- drop_id_to_type(drop_id,npools);
  
  return(cbind(cell_type, drop_id=rep(drop_id,npools), drop_type=rep(drop_type,npools)))
}

auto1d_threshold <- function(ab_log_norm, sigmas0) {
  rst <- NULL
  
  npools <- nrow(ab_log_norm)/2;

  thresholds <- sapply(1:nrow(ab_log_norm), function(i) {
    rst[[i]] <<- normalmixEM(ab_log_norm[i,],sigma=sigma0s)
    cat("mus: ", rst[[i]]$mu, 
        " sigmas: ", rst[[i]]$sigma, 
        " mean mus: ", mean(c(rst[[i]]$mu[1], rst[[i]]$mu[2])),
        "\n")
    mean(c(rst[[i]]$mu[1], rst[[i]]$mu[2]))
  })
  cat(thresholds, "\n")

  
  cell_type = NULL;

  for(i in 1:npools) {
    print(i);
  
    cell_type <- c(cell_type, apply(ab_log_norm[c(i,i+npools),],2,function(x){
      paste(sort(unlist(which(x>thresholds[c(i,i+npools)]))),collapse="_")
    }));
  }  

  ## ok, we'll filter mouse and human separately
  drop_id <- apply(ab_log_norm,2,function(x){
    ##paste(sort(unlist(which(x>0))),collapse="_")
  
    paste(sort(unlist(which(x>thresholds))),collapse="_")
  })
  
  drop_type <- drop_id_to_type(drop_id,npools);
  
  return(cbind(cell_type, drop_id=rep(drop_id,npools), drop_type=rep(drop_type,npools)))
  
}

manual1d_threshold <- function(ab_log, manual_thresholds) {
  cell_type = NULL;
  
  npools <- nrow(ab_log)/2;
  
  for(i in 1:npools) {
    print(i);
    
    cell_type <- c(cell_type, apply(ab_log[c(i,i+npools),],2,function(x){
      paste(sort(unlist(which(x>manual_thresholds[c(i,i+npools)]))),collapse="_")
    }));
  }  
  
  drop_id <- apply(ab_log,2,function(x){
    paste(sort(unlist(which(x>manual_thresholds))),collapse="_")
  });
  
  drop_type <- drop_id_to_type(drop_id,npools);
  
  ##browser();
  return(cbind(cell_type, drop_id=rep(drop_id,npools), drop_type=rep(drop_type,npools)))
}

#centers <- rbind(c(-2.5,-2.5),
#                  c(2.5,0),
#                  c(0,2.5),
#                  c(2.5,2.5));
# ab_droplet_ids <- rep("", ncol(ab_log_norm));
# kmeans_rst <- NULL;
# for(i in 1:npools) {
#   ##clara_rst <- clara(t(ab_log_norm[c(i,i+npools),]),4,metric="manhattan");
#   kmeans_rst[[i]] <- kmeans(t(ab_log_norm[c(i,i+npools),]),4,nstart=30,centers=centers);
#   
#   print(kmeans_rst[[i]]$centers)
#   
#   ## identifying centers for first ab
#   # ab1_indices <- which(kmeans_rst$cluster == which((kmeans_rst$centers[,1] > 0) & (kmeans_rst$centers[,2] < 0)))
#   # ab2_indices <- which(kmeans_rst$cluster == which((kmeans_rst$centers[,2] > 0) & (kmeans_rst$centers[,1] < 0)))
#   # ab1_2_indices <- which(kmeans_rst$cluster == which((kmeans_rst$centers[,1] > 0) & (kmeans_rst$centers[,2] > 0)))
#   ab1_indices <- which(kmeans_rst[[i]]$cluster == 2);
#   ab2_indices <- which(kmeans_rst[[i]]$cluster == 3);
#   ab1_2_indices <- which(kmeans_rst[[i]]$cluster == 4);
#   
#   ab_droplet_ids[ab1_indices] = sapply(ab_droplet_ids[ab1_indices], function(x) {
#     if(x == "") {
#       i
#     } else{
#       paste(x, i, sep="_");
#     }
# 
#     if(x == "") {
#       i+npools
#     } else {
#       paste(x, i+npools, sep="_");
#     }
#   })
# 
#   ab_droplet_ids[ab1_2_indices] = sapply(ab_droplet_ids[ab1_2_indices], function(x) {
#     if(x == "") {
#       paste(i, i+npools, sep="_");
#     } else {
#       paste(x, i, i+npools, sep="_");
#     }
#   })
#}
# 
# 
# ## 2d thresholding
# total_ab_log_norm_2d = total_ab_log_norm[total_ab_log_norm$ab_cell_type!="1_2",]
# total_ab_log_norm_2d = total_ab_log_norm_2d[total_ab_log_norm_2d$demux_droplet_type!="doublet",]
# total_ab_log_norm_2d = total_ab_log_norm_2d[total_ab_log_norm_2d$demux_droplet_type!="unassigned",]
# total_ab_log_norm_2d = total_ab_log_norm_2d[total_ab_log_norm_2d$ab_droplet_type!="unstained",]
# 
# summary(as.factor(total_ab_log_norm_2d$ab_cell_type[total_ab_log_norm_2d$ab_droplet_type=="singlet" & 
#                                                       total_ab_log_norm_2d$demux_droplet_type=="donor0"]))
# 
# summary(as.factor(total_ab_log_norm_2d$ab_cell_type[total_ab_log_norm_2d$ab_droplet_type=="singlet" & 
#                                                       total_ab_log_norm_2d$demux_droplet_type=="donor1"]))
# 
# summary(as.factor(total_ab_log_norm_2d$ab_cell_type[total_ab_log_norm_2d$ab_droplet_type!="singlet" & 
#                                                       total_ab_log_norm_2d$demux_droplet_type=="donor0"]))
# 
# summary(as.factor(total_ab_log_norm_2d$ab_cell_type[total_ab_log_norm_2d$ab_droplet_type!="singlet" & 
#                                                       total_ab_log_norm_2d$demux_droplet_type=="donor1"]))
# 
# summary(total_ab_log_norm_2d$tx_type[total_ab_log_norm_2d$demux_droplet_type=="donor0"&total_ab_log_norm_2d$ab_droplet_type=="singlet"])
# 
# summary(total_ab_log_norm_2d$tx_type[total_ab_log_norm_2d$demux_droplet_type=="donor1"&total_ab_log_norm_2d$ab_droplet_type=="singlet"])
# ##ggplot(total_ab_log_norm_2d, aes(ct1, ct2, color=ab_cell_type))+geom_point()+facet_grid(ab_droplet_type~demux_droplet_type)
# 

