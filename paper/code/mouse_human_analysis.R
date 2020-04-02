## let's make this code generalizable to a lot of batches
library(Seurat);
library(mixtools);
library(Matrix);
library(ggplot2);
library(ggridges)
library(data.table);
library(cluster)
library(dplyr)

source("threshold.R")
source("correlate_rates.R");

setwd("~/Box/Papers/In Preparation/2019.Hwang.scitoseq/Code")

### open experiment files (h5 here)
##experiment <- "200k_human_human_TB_5batch"
experiment <- "200k_human_human_TB_5batch"
##experiment <- "100k_univ"
dir.create(experiment, showWarnings = TRUE, recursive = FALSE, mode = "0777")
raw <- Read10X_h5(paste("../Text/", experiment, "/filtered_feature_bc_matrix.h5",sep=""))
tx <- raw$`Gene Expression`[,which(colSums(raw$`Gene Expression`)>0)]
ab <- raw$`Antibody Capture`[,which(colSums(raw$`Gene Expression`)>0)];

## ccd = cell_containing_drops
##ccd_indices <- intersect(which(colSums(tx)>0), which(colSums(ab)>0))
##ab <- ab[,ccd_indices];
##tx <- tx[,ccd_indices];

ab_colSums = colSums(ab);
tx_colSums = colSums(tx);

if(experiment %in% c("20k_human_human", "50k_human_human", "100k_human_human_TB_5batch",  "200k_human_human_TB_5batch")) {
  donor_ids <- fread(paste("../Text/", experiment, "/donor_ids.tsv",sep=""))
}

##tx <- tx;##[,colSums(tx)>1000]

if(experiment %in% c("20k_univ", "20k_direct", "100k_univ")) {
  ## call transcriptome cell types
  tx_hg_perc = colSums(tx[grep("hg19",rownames(tx)),])/(colSums(tx[grep("hg19",rownames(tx)),])
                                                       +colSums(tx[grep("mm10",rownames(tx)),]))
  tx_mm_perc = colSums(tx[grep("mm10",rownames(tx)),])/(colSums(tx[grep("hg19",rownames(tx)),])
                                                       +colSums(tx[grep("mm10",rownames(tx)),]))
  tx_type = rep("doublet_cross_ct", length(tx_hg_perc))
  tx_type[tx_hg_perc > .92] = "hg19";
  tx_type[tx_mm_perc > .92] = "mm10";
  tx_hg19_nsinglets <- sum(tx_type=="hg19");
  tx_mm10_nsinglets <- sum(tx_type=="mm10");
  tx_nmultiplets <- sum(tx_type=="doublet_cross_ct");
  ##tx_nsinglets <- tx_hg19_nsinglets + tx_mm10_nsinglets;
  tx_nsinglets <- c(tx_hg19_nsinglets, tx_mm10_nsinglets);
  names(tx_nsinglets) <- c("human","mouse")
  ndrops <- length(tx_type);
  mr_tx <- (tx_nmultiplets/ndrops)/(1-(tx_hg19_nsinglets^2+tx_mm10_nsinglets^2)/sum(tx_nsinglets)^2)
}

if(experiment %in% c("20k_human_human", "50k_human_human", "100k_human_human_TB_5batch",  "200k_human_human_TB_5batch")) {
  ## call transcriptome cell types
  tx_b_perc = tx[grep("^MS4A1$",rownames(tx)),]/(tx[grep("^MS4A1$",rownames(tx)),]
                                                        +tx[grep("^CD3E$",rownames(tx)),])
  tx_t_perc = tx[grep("^CD3E$",rownames(tx)),]/(tx[grep("^MS4A1$",rownames(tx)),]
                                                        +tx[grep("^CD3E$",rownames(tx)),])
  tx_type = rep("doublet_cross_ct", length(tx_b_perc))
  tx_type[tx_b_perc > .9] = "b";
  tx_type[tx_t_perc > .9] = "t";
  tx_type[is.nan(tx_t_perc) & is.nan(tx_b_perc)] = "unknown"
  tx_b_nsinglets <- sum(tx_type=="b");
  tx_t_nsinglets <- sum(tx_type=="t");
  tx_nmultiplets <- sum(tx_type=="doublet_cross_ct");
  ##tx_nsinglets <- tx_b_nsinglets + tx_t_nsinglets;
  tx_nsinglets <- c(tx_b_nsinglets, tx_t_nsinglets);
  names(tx_nsinglets) <- c("b","t")
  ndrops <- length(tx_type)
  mr_tx <- (tx_nmultiplets/(tx_nmultiplets+sum(tx_nsinglets)))/(1-(tx_b_nsinglets^2+tx_t_nsinglets^2)/sum(tx_nsinglets)^2)
}

npools <- nrow(ab)/2;

## tx count normalization + log normal normalization 
ab_norm <- ab*(mean(tx_colSums)/tx_colSums);
ab_ab_norm <- ab*(mean(ab_colSums)/ab_colSums);
ab_log_norm = log1p(ab_norm)
ab_log <- log1p(ab);

# priors for the gaussian mixture model and manual thresholds for single and double positive populations
if(experiment %in% c("20k_human_human")) {
    mu0s <- rbind(c(2.5,5.5),
                c(2.5,5.5),
                c(2.5,5.5),
                c(2.5,5.5))
    ## fixing sigmals works better for this
    sigma0s <- c(0.1);
    manual_thresholds <- rbind(c(3.8,3.8),
                               c(4.3,4))
    
}
    
if(experiment %in% c("50k_human_human")) {
      mu0s <- rbind(c(3,5.5),
                       c(4,5.5),
                       c(3,5.5),
                       c(4,5.5))
         sigma0s <- c(0.1);
         manual_thresholds <- rbind(c(6,5),
                                    c(6,4.2))
         manual_thresholds <- rbind(c(4.3,4.3),
                                    c(4.7,4.7))
         
}

if(experiment %in% c("200k_human_human_TB_5batch")) {
  mu0s <- rbind(c(4,6),
                c(2.5,5.5),
                c(2.5,5),
                c(3.5,6),
                 c(3.5,6),
                 c(4,6),
                 c(4,6),
                 c(2,5),
                  c(4,6),
                  c(4, 6))
  sigma0s <- c(0.1);
  manual_thresholds <- rbind(c(4.8,5),
                             c(5.1,5.2),
                             c(3.8,3.8),
                             c(4.8,5),
                             c(4.6,5));
  
  
  
}

if(experiment %in% c("100k_human_human_TB_5batch")) {
  mu0s <- rbind(c(4,6),
                c(2.5,5.5),
                c(2.5,5),
                c(3.5,6),
                c(3.5,6),
                c(4,6),
                c(4,6),
                c(2,5),
                c(4,6),
                c(4, 6))
  manual_thresholds <- rbind(c(4.4,4.8),
                             c(4.8,4.9),
                             c(3.8,3),
                             c(4.5,4.4),
                             c(4.4,4.4));
  manual_ab_norm_thresholds <- rbind(c(4.5,5),
                                     c(5,5),
                                     c(3.8,3.5),
                                     c(4.6,4.8),
                                     c(4.4,4.8));

  sigma0s <- c(0.1);
}

if(experiment %in% c("20k_univ","20k_direct")) {
  mu0s <- rbind(c(3,5.5),
                c(4,5.5),
                c(3,5.5),
                c(4,5.5))
  sigma0s <- c(0.1,0.1);
  manual_thresholds <- rbind(c(6,5),
                             c(6,4.2))

}

if(experiment %in% c("100k_univ")) {
  mu0s <- rbind(c(2.5,6.5),
                c(3,6.5),
                c(2.5,6.5),
                c(4.5,7.5),
                c(4.5,7.5),
                c(3.5,6),
                c(3.5,5.5),
                c(3,5.5),
                c(3.5,6),
                c(4,7))
  manual_thresholds <- rbind(c(4,4.4),
                             c(4.5,4.2),
                             c(4.5,4),
                             c(6,4.4),
                             c(5.5,5.1));
  sigma0s <- c(0.1,0.1);
}
    
mus <- NULL;
sigmas <- NULL;
for(i in 1:nrow(ab_log_norm)) {
  x <- ab_log_norm[i,];
  rst = normalmixEM(x, mu=mu0s[i,], sigma=sigma0s,maxit=1000000,epsilon=1e-11);
  mus <<- rbind(mus, data.frame(neg=rst$mu[1],pos=rst$mu[2]))
  sigmas <<- rbind(sigmas, data.frame(neg=rst$sigma[1],pos=rst$sigma[2]))
  
}
mus_mean <- rbind(colMeans(mus[c(1,3),]),colMeans(mus[c(1,3),]),colMeans(mus[c(2,4),]),colMeans(mus[c(2,4),]))
sigmas_mean <- rbind(colMeans(sigmas[c(1,3),]),colMeans(sigmas[c(1,3),]),colMeans(sigmas[c(2,4),]),colMeans(sigmas[c(2,4),]))

# this is normalizing to the difference
ab_log_norm = t(sapply(1:nrow(ab_log_norm), function(i) {
  x <- ab_log_norm[i,];
  rst = normalmixEM(x, mu=mu0s[i,], sigma=sigma0s, maxit=100000,epsilon=1e-9)
  return_x <- x;
  threshold <- rst$x[which.min(abs(rst$posterior[,2]-rst$posterior[,1]))]
  print(threshold);

    return_x <- (x-threshold)/rst$sigma[1]
  
  print(rst$mu);
  print(rst$sigma)

  return(return_x)
}))


##########################################################################
## HTODemux threhsolding using raw data
##########################################################################

htodemux_rst <- htodemux_threshold(ab);
##htodemux_cell_type <- htodemux_rst[[1]]
##htodemux_drop_type <- htodemux_rst[[2]]

## standardize
ab_log_norm <- t(apply(ab_log_norm, 1, function(x) {(x-mean(x))/sd(x)}))

##################################################################
## 1-d automatic thresholding
##################################################################

##let's put in some automatic filtering code
auto1d_rst <- auto1d_threshold(ab_log_norm, sigma0s)
##auto1d_cell_type <- auto1d_rst[[1]];
##auto1d_drop_type <- auto1d_rst[[2]];

##################################################################
## 1-d manual thresholding
##################################################################
manual1d_rst <- manual1d_threshold(ab_log, manual_thresholds)

##manual1d_rst <- manual1d_threshold(log1p(ab_ab_norm), manual_ab_norm_thresholds)

##manual1d_cell_type <- manual1d_rst[[1]];
##manual1d_drop_type <- manual1d_rst[[2]];

###################################################################
## 2-d thresholding per pool
###################################################################
# 

htodemux_exp_obs_df <- correlate_rates(ab, htodemux_rst[,2])
auto1d_exp_obs_df <- correlate_rates(ab, auto1d_rst[,2])
manual1d_exp_obs_df <- correlate_rates(ab, manual1d_rst[,2])

ggsave(ggplot(htodemux_exp_obs_df, aes(exp, obs))+geom_point()+theme_bw()+geom_abline(intercept=0,slope=1), file=paste(experiment,"/htodemux_exp_obs.png",sep=""),width=3,height=3)

ggsave(ggplot(auto1d_exp_obs_df, aes(exp, obs))+geom_point()+theme_bw()+geom_abline(intercept=0,slope=1), file=paste(experiment,"/auto1d_exp_obs.png",sep=""),width=3,height=3)

ggsave(ggplot(manual1d_exp_obs_df, aes(exp, obs))+geom_point()+theme_bw()+geom_abline(intercept=0,slope=1), file=paste(experiment,"/manual1d_exp_obs.png",sep=""),width=3,height=3)

# if(experiment == "100k_human_human_TB_5batch"){
#   ggsave(ggplot(exp_obs_df, aes(exp, obs))+geom_point()+theme_bw()+scale_x_continuous(limits=c(0,0.01))+scale_y_continuous(limits=c(0,0.01))+geom_abline(intercept=0,slope=1), file=paste(experiment,"/exp_obs.png",sep=""),width=3,height=3)
# } else if(experiment == "200k_human_human_TB_5batch") {
#   ggsave(ggplot(exp_obs_df, aes(exp, obs))+geom_point()+theme_bw()+scale_x_continuous(limits=c(0,0.01))+scale_y_continuous(limits=c(0,0.01))+geom_abline(intercept=0,slope=1), file=paste(experiment,"/exp_obs.png",sep=""),width=3,height=3)
#   
# }


##########################################################################
## HTODemux threhsolding by separating out all the batches
##########################################################################
# ab_resolved <- NULL;
# for(i in 1:npools) {
#   cd4_index <- i;
#   cd20_index <- i+npools;
#   ab_resolved <- cbind(ab_resolved,ab[c(cd4_index, cd20_index),]);
# }
# # Setup Seurat object
# ab_resolved_obj <- CreateSeuratObject(counts = ab_resolved)
# 
# # Add HTO data as a new assay independent from RNA
# ab_resolved_obj[["HTO"]] <- CreateAssayObject(counts = ab_resolved)
# 
# # Normalize HTO data, here we use centered log-ratio (CLR) transformation
# ab_resolved_obj <- NormalizeData(ab_resolved_obj, assay = "HTO", normalization.method = "LogNormalize")
# 
# rst_resolved <- HTODemux(ab_resolved_obj, assay = "HTO", kfunc="kmeans", positive.quantile = 0.95)
# ab_resolved_obj <- rst_resolved[[1]];
# discrete_resolved <- rst_resolved[[2]];
# 
# table(ab_resolved_obj$HTO_classification.global)
# 
# ab_resolved_hto <- NULL;
# 
# for(i in 1:npools) {
#   cd4_index <- i;
#   cd20_index <- i+npools;
#   ab_resolved_hto <- cbind(ab_resolved_hto,discrete[c(cd4_index, cd20_index),]);
#   ##  ab_resolved_hto <- cbind(ab_resolved_hto, ab[start:end,which(discrete[i,]==i)]),
# }
# 
# labels <- apply(ab_resolved_hto,2,paste,collapse="_")
# 
# ab_droplet_ids <- apply(discrete, 2, function(x) {
#   paste(which(x==1),collapse="_")
# })

# 
# ########################################################################
# ## threshold first and call droplets first
# ########################################################################
# 
# # Setup Seurat object
# ab_obj <- CreateSeuratObject(counts = ab)
# 
# # Add HTO data as a new assay independent from RNA
# ab_hto <- NULL;
# 
# for(i in 1:npools) {
#   cd4_index <- i;
#   cd20_index <- i+npools;
#   ab_hto <- cbind(ab_hto, Matrix::colSums(ab[c(cd4_index,cd20_index,])))
# }
# 
# ab_obj[["HTO"]] <- CreateAssayObject(counts = ab_hto)
# 
# # Normalize HTO data, here we use centered log-ratio (CLR) transformation
# ab_obj <- NormalizeData(ab_obj, assay = "HTO", normalization.method = "LogNormalize")
# 
# rst <- HTODemux(ab_obj, assay = "HTO", positive.quantile = 0.99)
# ab_obj <- rst[[1]];
# discrete <- rst[[2]];
# 
# table(ab_obj$HTO_classification.global)
# 
# ab_resolved_hto <- NULL;
# 
# for(i in 1:npools) {
#   cd4_index <- i;
#   cd20_index <- i+npools;
#   ab_resolved_hto <- cbind(ab_resolved_hto, discrete[c(cd4_index, cd20_index),]);
#   ##  ab_resolved_hto <- cbind(ab_resolved_hto, ab[start:end,which(discrete[i,]==i)]),
# }
# 
# labels <- apply(ab_resolved_hto,2,paste,collapse="_")

# 
# centers <- rbind(c(-2.5,-2.5),
#                  c(2.5,0),
#                  c(0,2.5),
#                  c(2.5,2.5));
# ab_droplet_ids <- rep("", ncol(ab_log_norm));
# gmm_rst <- NULL;
# for(i in 1:npools) {
#   dat = center_scale(t(ab_log_norm[c(i,i+npools),]), mean_center = T, sd_scale = T)  # centering and scaling the data
#   
#   gmm <- GMM(dat, gaussian_comps=4, dist_mode = "eucl_dist", seed_mode = "random_subset", km_iter = 1000,
#              em_iter = 1000, verbose = F)##, dist_mode = "maha_dist");
# 
#   gmm_rst[[i]] = predict_GMM(dat, gmm$centroids, gmm$covariance_matrices, gmm$weights)    
#   
#   print(gmm$centroids)
#   
#   ## identifying centers for first ab
#   # ab1_indices <- which(kmeans_rst$cluster == which((kmeans_rst$centers[,1] > 0) & (kmeans_rst$centers[,2] < 0)))
#   # ab2_indices <- which(kmeans_rst$cluster == which((kmeans_rst$centers[,2] > 0) & (kmeans_rst$centers[,1] < 0)))
#   # ab1_2_indices <- which(kmeans_rst$cluster == which((kmeans_rst$centers[,1] > 0) & (kmeans_rst$centers[,2] > 0)))
#   ab1_indices <- which(gmm_rst[[i]]$cluster_labels == 2);
#   ab2_indices <- which(gmm_rst[[i]]$cluster_labels == 3);
#   ab1_2_indices <- which(gmm_rst[[i]]$cluster_labels == 4);
#   
#   ab_droplet_ids[ab1_indices] = sapply(ab_droplet_ids[ab1_indices], function(x) {
#     if(x == "") {
#       i
#     } else{
#       paste(x, i, sep="_");
#     }
#   })
#   
#   ab_droplet_ids[ab2_indices] = sapply(ab_droplet_ids[ab2_indices], function(x) {
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
# }

## discrepancies
# disc = rep(0,ndrops);
# disc[setdiff(tx_doublets_indices,ab_cross_species_doublets_indices)] <- 1;

cell_drop_type <- data.frame(htodemux=htodemux_rst,
                             manual1d=manual1d_rst,
                             auto1d=auto1d_rst)

## some plotting functions 1
total_ab_tx_norm_log = NULL;
for(i in 1:npools) {
  print(i);
  total_ab_tx_norm_log <- rbind(total_ab_tx_norm_log,
                        data.frame(ct1=log1p(ab_norm[i,]),
                                   ct2=log1p(ab_norm[i+npools,]),
                                   batch=paste("batch",i,sep="")))#,
  ##                                   ab_droplet_type=ab_droplet_type,
  ##                                   disc=disc))
}

total_ab_tx_norm_log <- cbind(total_ab_tx_norm_log, cell_drop_type)

p = ggplot(total_ab_tx_norm_log, aes(ct1,ct2,color=batch))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none")
#+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
##ggplot(total_ab_tx_norm_log, aes(ct1,ct2,color=batch))+geom_point()
ggsave(p, file=paste(experiment, "/total_ab_tx_norm_log_batch.png",sep=""),width=3,height=3)
p

for(i in 1:npools) {
  p = ggplot(total_ab_tx_norm_log[total_ab_tx_norm_log$batch==paste("batch",i,sep=""),], aes(ct1,ct2,color=htodemux.cell_type))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none");##+facet_grid(batch~.)
  #+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
  ggsave(p, file=paste(experiment,"/total_ab_tx_norm_log_htodemux_cell_type_batch",i,".png",sep=""),width=3,height=3)
  p
}

for(i in 1:npools) {
  p = ggplot(total_ab_tx_norm_log[total_ab_tx_norm_log$batch==paste("batch",i,sep=""),], aes(ct1,ct2,color=auto1d.cell_type))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none");##+facet_grid(batch~.)
  #+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
  ggsave(p, file=paste(experiment,"/total_ab_tx_norm_log_auto1d_cell_type_batch",i,".png",sep=""),width=3,height=3)
  p
}

for(i in 1:npools) {
  p = ggplot(total_ab_tx_norm_log[total_ab_tx_norm_log$batch==paste("batch",i,sep=""),], aes(ct1,ct2,color=manual1d.cell_type))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none");##+facet_grid(batch~.)
  #+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
  ggsave(p, file=paste(experiment,"/total_ab_tx_norm_log_manual1d_cell_type_batch",i,".png",sep=""),width=3,height=3)
  p
}



## some plotting functions 2
total_ab_ab_norm_log = NULL;
for(i in 1:npools) {
  print(i);
  total_ab_ab_norm_log <- rbind(total_ab_ab_norm_log,
                                data.frame(ct1=log1p(ab_ab_norm[i,]),
                                           ct2=log1p(ab_ab_norm[i+npools,]),
                                           batch=paste("batch",i,sep="")))#,
  ##                                   ab_droplet_type=ab_droplet_type,
  ##                                   disc=disc))
}

total_ab_ab_norm_log <- cbind(total_ab_ab_norm_log, cell_drop_type)

p = ggplot(total_ab_ab_norm_log, aes(ct1,ct2,color=batch))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none")
#+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
##ggplot(total_ab_ab_norm_log, aes(ct1,ct2,color=batch))+geom_point()
ggsave(p, file=paste(experiment, "/total_ab_ab_norm_log_batch.png",sep=""),width=3,height=3)
p

for(i in 1:npools) {
  p = ggplot(total_ab_ab_norm_log[total_ab_ab_norm_log$batch==paste("batch",i,sep=""),], aes(ct1,ct2,color=htodemux.cell_type))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none");##+facet_grid(batch~.)
  #+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
  ggsave(p, file=paste(experiment,"/total_ab_ab_norm_log_htodemux_cell_type_batch",i,".png",sep=""),width=3,height=3)
  p
}

for(i in 1:npools) {
  p = ggplot(total_ab_ab_norm_log[total_ab_ab_norm_log$batch==paste("batch",i,sep=""),], aes(ct1,ct2,color=auto1d.cell_type))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none");##+facet_grid(batch~.)
  #+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
  ggsave(p, file=paste(experiment,"/total_ab_ab_norm_log_auto1d_cell_type_batch",i,".png",sep=""),width=3,height=3)
  p
}

for(i in 1:npools) {
  p = ggplot(total_ab_ab_norm_log[total_ab_ab_norm_log$batch==paste("batch",i,sep=""),], aes(ct1,ct2,color=manual1d.cell_type))+geom_point(size=0.02)+theme_bw()+theme(legend.position = "none");##+facet_grid(batch~.)
  #+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
  ggsave(p, file=paste(experiment,"/total_ab_ab_norm_log_manual1d_cell_type_batch",i,".png",sep=""),width=3,height=3)
  p
}


## some plotting functions 3
total_ab_log = NULL;
for(i in 1:npools) {
  print(i);
  total_ab_log <- rbind(total_ab_log,
                        data.frame(ct1=ab_log[i,],
                                   ct2=ab_log[i+npools,],
                                   batch=paste("batch",i,sep="")))#,
##                                   ab_droplet_type=ab_droplet_type,
##                                   disc=disc))
}

total_ab_log <- cbind(total_ab_log, cell_drop_type)

p = ggplot(total_ab_log, aes(ct1,ct2,color=batch))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none")
#+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
##ggplot(total_ab_log, aes(ct1,ct2,color=batch))+geom_point()
ggsave(p, file=paste(experiment, "/total_ab_log_batch.png",sep=""),width=3,height=3)
p

for(i in 1:npools) {
  p = ggplot(total_ab_log[total_ab_log$batch==paste("batch",i,sep=""),], aes(ct1,ct2,color=htodemux.cell_type))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none");##+facet_grid(batch~.)
  #+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
  ggsave(p, file=paste(experiment,"/total_ab_log_htodemux_cell_type_batch",i,".png",sep=""),width=3,height=3)
  p
}

for(i in 1:npools) {
  p = ggplot(total_ab_log[total_ab_log$batch==paste("batch",i,sep=""),], aes(ct1,ct2,color=auto1d.cell_type))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none");##+facet_grid(batch~.)
  #+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
  ggsave(p, file=paste(experiment,"/total_ab_log_auto1d_cell_type_batch",i,".png",sep=""),width=3,height=3)
  p
}

for(i in 1:npools) {
  p = ggplot(total_ab_log[total_ab_log$batch==paste("batch",i,sep=""),], aes(ct1,ct2,color=manual1d.cell_type))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none");##+facet_grid(batch~.)
  #+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
  ggsave(p, file=paste(experiment,"/total_ab_log_manual1d_cell_type_batch",i,".png",sep=""),width=3,height=3)
  p
}

## this will identify mouse/human doublets with batch on top
ggplot(total_ab_log, aes(ct1,ct2,color=batch))+geom_point()
##ggplot(total_ab_log, aes(ct1,ct2,color=ab_droplet_type))+geom_point()

## some plotting functions 4
## first plot the resolved data
total_ab_log_norm = NULL;
for(i in 1:npools) {
  print(i);

  total_ab_log_norm <- rbind(total_ab_log_norm,
        data.frame(ct1=ab_log_norm[i,],
                   ct2=ab_log_norm[i+npools,],
                   batch=paste("batch",i,sep=""),
                   tx_type=tx_type))
  ##                 disc=disc))
}


total_ab_log_norm$demux_droplet_type=donor_ids$donor_id[match(sapply(rownames(total_ab_log_norm),function(x){substr(x,1,18)}), donor_ids$cell)]
total_ab_log_norm <- cbind(total_ab_log_norm, cell_drop_type)
####### write the output file for further plotting using python (generating the resolved scatter and density plots are done with python jupyter notebook)
write.table(total_ab_log_norm, file=paste(experiment,"/total_ab_log_norm.tab",sep=""), row.names=T, col.names=T, quote=F, sep="\t")

htodemux_summary <- table(total_ab_log_norm[,c("batch","htodemux.cell_type")])[,-1]
htodemux_summary <- rbind(htodemux_summary, colSums(htodemux_summary))
htodemux_summary_rat <- t(apply(htodemux_summary, 1, function(x) {x/sum(x)}))

auto1d_summary <- table(total_ab_log_norm[,c("batch","auto1d.cell_type")])[,-1]
auto1d_summary <- rbind(auto1d_summary, colSums(auto1d_summary))
auto1d_summary_rat <- t(apply(auto1d_summary, 1, function(x) {x/sum(x)}))

manual1d_summary <- table(total_ab_log_norm[,c("batch","manual1d.cell_type")])[,-1]
manual1d_summary <- rbind(manual1d_summary, colSums(manual1d_summary))
manual1d_summary_rat <- t(apply(manual1d_summary, 1, function(x) {x/sum(x)}))


p = ggplot(total_ab_log_norm, aes(ct1,ct2,color=batch))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none")
#+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
##ggplot(total_ab_log_norm, aes(ct1,ct2,color=batch))+geom_point()
ggsave(p, file=paste(experiment, "/total_ab_log_norm_batch.png",sep=""),width=3,height=3)
p

p = ggplot(total_ab_log_norm, aes(ct1,ct2,color=htodemux.cell_type))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none")
#+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
ggsave(p, file=paste(experiment,"/total_ab_log_norm_htodemux_cell_type.png",sep=""),width=3,height=3)
p

p = ggplot(total_ab_log_norm, aes(ct1,ct2,color=manual1d.cell_type))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none")
#+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
ggsave(p, file=paste(experiment,"/total_ab_log_norm_manual1d_cell_type.png",sep=""),width=3,height=3)
p

p = ggplot(total_ab_log_norm, aes(ct1,ct2,color=auto1d.cell_type))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none")
#+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
ggsave(p, file=paste(experiment,"/total_ab_log_norm_auto1d_cell_type.png",sep=""),width=3,height=3)
p


for(i in 1:npools) {
  p = ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",i,sep=""),], aes(ct1,ct2,color=htodemux.cell_type))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none");##+facet_grid(batch~.)
  #+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
  ggsave(p, file=paste(experiment,"/total_ab_log_norm_htodemux_cell_type_batch",i,".png",sep=""),width=3,height=3)
  p
}

for(i in 1:npools) {
  p = ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",i,sep=""),], aes(ct1,ct2,color=auto1d.cell_type))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none");##+facet_grid(batch~.)
  #+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
  ggsave(p, file=paste(experiment,"/total_ab_log_norm_auto1d_cell_type_batch",i,".png",sep=""),width=3,height=3)
  p
}

for(i in 1:npools) {
  p = ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",i,sep=""),], aes(ct1,ct2,color=manual1d.cell_type))+geom_point(size=0.2)+theme_bw()+theme(legend.position = "none");##+facet_grid(batch~.)
  #+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
  ggsave(p, file=paste(experiment,"/total_ab_log_norm_manual1d_cell_type_batch",i,".png",sep=""),width=3,height=3)
  p
}


# ggsave(ggplot(total_ab_log_norm, aes(ct1,ct2,color=ab_droplet_type))
#        +geom_point()
#        +theme_bw()
#        +scale_color_manual(values=c("blue","red","purple","gray")), file=paste(experiment,"/total_ab_log_norm_type.png",sep=""),width=3,height=3)
# ##ggplot(total_ab_log_norm, aes(ct1,ct2,color=batch))+geom_bin2d()+facet_grid(batch~.)


p = ggplot(total_ab_log_norm, aes(ct1,ct2,color=batch))+stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+scale_fill_distiller(palette=4, direction=-1)+theme_bw()+theme(legend.position = "none")
#+scale_x_continuous(limits=c(-6,8))+scale_y_continuous(limits=c(-6,8))
ggsave(p, file=paste(experiment,"/total_ab_log_norm_density.png",sep=""),width=3,height=3)
p

##ggplot(total_ab_log_norm, aes(ct1,ct2,color=ab_cell_type))+geom_point()+facet_grid(batch~.)
##ggplot(total_ab_log_norm[total_ab_log_norm$ab_droplet_type=="doublet_same_ct",], aes(ct1,ct2,color=ab_droplet_type))+geom_point()
##ggplot(total_ab_log_norm, aes(human,mouse,color=as.factor(disc)))+geom_point()

ab_merged <- rbind(colSums(ab[1:npools,]), colSums(ab[(npools+1):(2*npools),]))
ab_merged <- ab_merged*(mean(tx_colSums)/tx_colSums);
ab_merged_log = log1p(ab_merged)
ab_merged_log_norm = t(apply(ab_merged_log,1, function(x) {
   (x-mean(x))/sd(x);
}))
total_ab_merged_log_norm <- data.frame(ct1=ab_merged_log_norm[1,],
                              ct2=ab_merged_log_norm[2,],
                              tx_type=tx_type)
total_ab_merged_log_norm <- cbind(total_ab_merged_log_norm, cell_drop_type[1:ncol(ab),c(3,6,9)])
                              #disc=disc)

#source("ridge_plot.R")


##ggplot(total_ab_merged_log_norm, aes(ct1,ct2,color=tx_type))+geom_point()
# 
# ggsave(ggplot(total_ab_merged_log_norm, aes(ct1,ct2,color=ab_droplet_type))
#        +geom_point(size=0.2)
#        +theme_bw()
#        +scale_x_continuous(limits=c(-2.5,3.5))
#        +scale_y_continuous(limits=c(-2.5,3.5))
#        +scale_color_manual(values=c("blue","red","purple","gray"))
#        +theme(legend.position = "none"), file=paste(experiment,"/total_ab_merged_log_norm_type.png",sep=""),width=3,height=3)
# 
# ggsave(ggplot(total_ab_merged_log_norm, aes(ct1,ct2,color=ab_droplet_type))
#        +stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)
#        +scale_fill_distiller(palette=4, direction=-1)
#        +theme_bw()
#        +scale_x_continuous(limits=c(-2.5,3.5))
#        +scale_y_continuous(limits=c(-2.5,3.5))
#        +theme(legend.position = "none"), file=paste(experiment,"total_ab_merged_log_norm_density.png",sep=""),width=3,height=3)
# ##ggplot(total_ab_merged_log_norm, aes(ct1,ct2,color=ab_droplet_type))+geom_point()

## plot genetic singlets but ab doublets (cells from same donor encapsulated together but can be resolved by barcodes)

# 
# ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",1,sep=""),], aes(ct1,ct2,color=factor(kmeans_rst[[1]]$cluster)))+geom_point()
# ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",2,sep=""),], aes(ct1,ct2,color=factor(kmeans_rst[[2]]$cluster)))+geom_point()
# ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",3,sep=""),], aes(ct1,ct2,color=factor(kmeans_rst[[3]]$cluster)))+geom_point()
# ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",4,sep=""),], aes(ct1,ct2,color=factor(kmeans_rst[[4]]$cluster)))+geom_point()
# ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",5,sep=""),], aes(ct1,ct2,color=factor(kmeans_rst[[5]]$cluster)))+geom_point()

# batch1_indices <- which(total_ab_log_norm$batch==paste("batch",1,sep=""));
# batch1_background_cor <- cor(total_ab_log_norm[intersect(batch1_indices, kmeans_rst[[1]]$cluster==1),c("ct1","ct2")])
# total_ab_log_norm[batch1_indices, c("ct2")] = lm(total_ab_log_norm[batch1_indices, c("ct2")]~total_ab_log_norm[batch1_indices, c("ct1")])$residual
# 
# 
# batch2_indices <- which(total_ab_log_norm$batch==paste("batch",2,sep=""));
# batch2_background_cor <- cor(total_ab_log_norm[intersect(batch2_indices, kmeans_rst[[1]]$cluster==2),c("ct1","ct2")])
# total_ab_log_norm[batch2_indices, c("ct2")] = lm(total_ab_log_norm[batch2_indices, c("ct2")]~total_ab_log_norm[batch2_indices, c("ct1")])$residual
# batch2_kmeans <- kmeans(total_ab_log_norm[batch2_indices, c("ct1","ct2")],4,nstart=30,centers=centers);

# total_ab_log_norm$kmeans[total_ab_log_norm$batch==paste("batch",1,sep="")] = kmeans_rst[[1]]$cluster;
# total_ab_log_norm$kmeans[total_ab_log_norm$batch==paste("batch",2,sep="")] = kmeans_rst[[2]]$cluster;
# total_ab_log_norm$kmeans[total_ab_log_norm$batch==paste("batch",3,sep="")] = kmeans_rst[[3]]$cluster;
# total_ab_log_norm$kmeans[total_ab_log_norm$batch==paste("batch",4,sep="")] = kmeans_rst[[4]]$cluster;
# total_ab_log_norm$kmeans[total_ab_log_norm$batch==paste("batch",5,sep="")] = kmeans_rst[[5]]$cluster;

# centers <- rbind(c(-.5,-.5),
#                  c(2,0),
#                  c(0,2),
#                  c(2,2));
# 
# ##kmeans_rst_together <- kmeans(total_ab_log_norm[, c("ct1","ct2")],5,nstart=100,centers=centers);
# kmeans_rst_together <- kmeans(total_ab_log_norm[, c("ct1","ct2")],4,nstart=100,centers=centers);
# total_ab_log_norm$kmeans_together = kmeans_rst_together$cluster;
# ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",1,sep=""),], aes(ct1,ct2,color=factor(kmeans_together)))+geom_point(size=0.2)

# library(ggExtra)
# 
# ##total_ab_log$hto_demux = labels;
# p = ggplot(total_ab_log[total_ab_log$batch==paste("batch",2,sep=""),], aes(ct1,ct2,color=factor(htodemux.cell_type)))+geom_point(size=0.2)
# ggMarginal(p, type="density",groupColour=TRUE)
# 
# total_ab_log_norm$hto_demux = labels;
# p = ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",2,sep=""),], aes(ct1,ct2,color=factor(hto_demux)))+geom_point(size=0.2)
# ggMarginal(p, type="density",groupColour=TRUE)
# 
# summary(factor(total_ab_log_norm$hto_demux))


# clara <- clara(total_ab_log_norm[, c("ct1","ct2")], 4, metric="jaccard", samples=1000) ##, dist_mode = "maha_dist", seed_mode =
# total_ab_log_norm$clara_together = clara$clustering
# ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",2,sep=""),], aes(ct1,ct2,color=factor(clara_together)))+geom_point(size=0.2)


# library(mclust);
# mod <- Mclust(total_ab_log_norm[, c("ct1","ct2")], 4) ##, dist_mode = "maha_dist", seed_mode =
# total_ab_log_norm$mclust_together = mod$classification
# ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",2,sep=""),], aes(ct1,ct2,color=factor(mclust_together)))+geom_point(size=0.2)
# 

# if(experiment == "100k_univ") {
#   jimmie_thresholds <- rbind(c(4,4.4),
#                              c(4.5,4.2),
#                              c(4.5,4),
#                              c(6,4.4),
#                              c(5.5,5.1));
# 
#   total_ab_log_norm$jimmie = "0_0";
#   for (i in 1:npools) {
#     print(i);
#     pool_index <- which(total_ab_log_norm$batch==paste("batch",i,sep=""))
#     total_ab_log_norm$jimmie[pool_index][total_ab_log2[pool_index,"ct1"] > jimmie_thresholds[i,1]] <- "0_1";
#     total_ab_log_norm$jimmie[pool_index][total_ab_log2[pool_index,"ct2"] > jimmie_thresholds[i,2]] <- "1_0";
#     total_ab_log_norm$jimmie[pool_index][intersect(which(total_ab_log2[pool_index,"ct1"] > jimmie_thresholds[i,1]),
#                                                    which(total_ab_log2[pool_index,"ct2"] > jimmie_thresholds[i,2]))] <- "1_1";
#   }
#   total_ab_log2$jimmie = total_ab_log_norm$jimmie;
# }
# 
# ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",2,sep=""),], aes(ct1,ct2,color=factor(jimmie)))+geom_point(size=0.2)
# ggplot(total_ab_log2[total_ab_log2$batch==paste("batch",2,sep=""),], aes(ct1,ct2,color=factor(jimmie)))+geom_point(size=0.2)
# summary(factor(total_ab_log_norm$jimmie))
# 

# 
# ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",1,sep=""),], aes(ct1,ct2,color=factor(gmm_rst[[1]]$cluster_labels)))+geom_point()
# ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",2,sep=""),], aes(ct1,ct2,color=factor(gmm_rst[[2]]$cluster_labels)))+geom_point()
# ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",3,sep=""),], aes(ct1,ct2,color=factor(gmm_rst[[3]]$cluster_labels)))+geom_point()
# ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",4,sep=""),], aes(ct1,ct2,color=factor(gmm_rst[[4]]$cluster_labels)))+geom_point()
# ggplot(total_ab_log_norm[total_ab_log_norm$batch==paste("batch",5,sep=""),], aes(ct1,ct2,color=factor(gmm_rst[[5]]$cluster_labels)))+geom_point()
# 
# 
# total_ab_log_norm$gmm[total_ab_log_norm$batch==paste("batch",1,sep="")] = gmm_rst[[1]]$cluster_labels;
# total_ab_log_norm$gmm[total_ab_log_norm$batch==paste("batch",2,sep="")] = gmm_rst[[2]]$cluster_labels;
# total_ab_log_norm$gmm[total_ab_log_norm$batch==paste("batch",3,sep="")] = gmm_rst[[3]]$cluster_labels;
# total_ab_log_norm$gmm[total_ab_log_norm$batch==paste("batch",4,sep="")] = gmm_rst[[4]]$cluster_labels;
# total_ab_log_norm$gmm[total_ab_log_norm$batch==paste("batch",5,sep="")] = gmm_rst[[5]]$cluster_labels;

# 
# ## some plotting functions
# total_ab_log_norm_batch = NULL;
# for(i in 1:npools) {
#   print(i);
#   total_ab_log_norm_batch <- rbind(total_ab_log_norm_batch,
#                              data.frame(ct1=ab_log_norm_batch[i,],
#                                         ct2=ab_log_norm_batch[i+npools,],
#                                         batch=paste("batch",i,sep=""),
#                                         ab_droplet_type=ab_droplet_type,
#                                         disc=disc))
# }
# 
# ## this will identify mouse/human doublets with batch on top
# ggplot(total_ab_log_norm_batch, aes(ct1,ct2,color=batch))+geom_point()
# ggplot(total_ab_log_norm_batch, aes(ct1,ct2,color=ab_droplet_type))+geom_point()
# ggplot(total_ab_log_norm_batch[total_ab_log_norm_batch$ab_droplet_type=="doublet_same_species",], aes(ct1,ct2,color=ab_droplet_type))+geom_point()
# ##ggplot(total_ab_log_norm_batch, aes(human,mouse,color=as.factor(disc)))+geom_point()



# 
# total_ab_500_log_norm_mouse = data.frame(batch1=ab_log_norm[3,],batch2=ab_log_norm[4,]
#                                          ,ab_droplet_type=ab_droplet_type,disc=disc)
# total_ab_500_log_norm_human = data.frame(batch1=ab_log_norm[1,],batch2=ab_log_norm[2,]
#                                          ,ab_droplet_type=ab_droplet_type,disc=disc)
# ggplot(total_ab_500_log_norm_mouse, aes(batch1,batch2,color=ab_droplet_type))+geom_point()
# ggplot(total_ab_500_log_norm_human, aes(batch1,batch2,color=ab_droplet_type))+geom_point()
# 
# 
# total_ab_500_log_norm_merged = rbind(data.frame(mouse=ab_log_norm[3,]
#                                                 +ab_log_norm[4,], 
#                                                 human=ab_log_norm[1,]
#                                                 +ab_log_norm[2,], ab_droplet_type=ab_droplet_type))
# ggplot(total_ab_500_log_norm_merged, aes(human,mouse,color=ab_droplet_type))+geom_point()
# 
# ## let's look at the transcriptome
# total_tx_500 = data.frame(mouse=colSums(tx[grep("mm10",rownames(tx)),]),
#                                human=colSums(tx[grep("hg19",rownames(tx)),]),
#                           ab_droplet_type=ab_droplet_type,tx_type=tx_type,disc=disc)
#                               
# ggplot(total_tx_500, aes(log10(human),log10(mouse),color=ab_droplet_type))+geom_point()
# ggplot(total_tx_500, aes(log10(human),log10(mouse),color=as.factor(disc)))+geom_point()

## let's threshold for doublets here

##total_ab_500 = rbind(data.frame(mouse=ab_log_norm_batch_adj[2,], human=ab_log_norm_batch_adj[1,],batch="batch1"),data.frame(mouse=ab_log_norm_batch_adj[4,],human=ab_log_norm_batch_adj[3,],batch="batch2"))
# 
# ab_merged <- rbind(ab[1,]+ab[2,], ab[3,]+ab[4,])
# ab_merged_log = log1p(ab_merged)
# ab_merged_log_norm = t(apply(ab_merged_log,1, function(x) {
#   (x-mean(x))/sd(x);
# }))
# 
# total_ab_500_merged_log_norm = data.frame(mouse=ab_merged_log_norm[1,], 
#                                          human=ab_merged_log_norm[2,])
