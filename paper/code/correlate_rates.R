perm = function(n, x) {
  factorial(n) / factorial(n-x)
}


correlate_rates <- function(ab, drop_id) {
  ndrops <- ncol(ab);
  npools <- nrow(ab)/2;
 
  drop_id <- drop_id[1:ndrops];
  # 
  # ## compare tx with ab doublets
  # overlap_doublets_indices <- intersect(grep("doublet_cross_ct", tx_type), grep("doublet_cross_ct", ab_droplet_type));
  # tx_doublets_indices <- grep("doublet", tx_type);
  # ab_doublets_indices <- grep("doublet", ab_droplet_type);
  # ab_cross_species_doublets_indices <- grep("doublet_cross_ct",ab_droplet_type);
  # cat("tx doublets: ", length(tx_doublets_indices), " prop tx overlap: ", length(overlap_doublets_indices)/length(tx_doublets_indices),"\n")
  # cat("ab doublets: ", length(ab_cross_species_doublets_indices), " prop ab overlap: ", length(overlap_doublets_indices)/length(ab_cross_species_doublets_indices),"\n")

  ## now let's compare the expected vs observed doublets of various types
  ab_ct1_nsinglets <- sapply(1:npools, function(i) {
    length(which(drop_id==i));
  })
  ab_ct2_nsinglets <- sapply(1:npools, function(i) {
    length(which(drop_id==npools+i));
  })
  ab_nsinglets <- cbind(ab_ct1_nsinglets, ab_ct2_nsinglets);
  names(ab_nsinglets) <- rownames(ab)
  
  ab_nmultiplets <- length(grep("_", drop_id))
  ab_nunstained <- length(which(drop_id==""))

    ## these are the expected doublet rates after accounting for doublets we can't detect
  ##mr_ab <- ab_nmultiplets/(ab_nsinglets+ab_nmultiplets)/(1-(h1_ncells^2+h2_ncells^2+m1_ncells^2+m2_ncells^2)/nsinglets^2);
  
  ## new code, here we rowsum first which is the number of cells of a particular batch
  ##mr_ab <- (ab_nmultiplets/(sum(ab_nsinglets)+ab_nmultiplets))/(1-sum(sapply(rowSums(ab_nsinglets), function(x){x^2}))/sum(ab_nsinglets)^2);
  ## old code...here we assume that batch barcode and antibody can separate all cells
  mr_ab <- (ab_nmultiplets/(sum(ab_nsinglets)+ab_nmultiplets))/(1-sum(sapply(ab_nsinglets, function(x){x^2}))/sum(ab_nsinglets)^2);
  
  multiplet_index = grep("_",drop_id);
  multiplet_array = sapply(drop_id[multiplet_index], function(x) {strsplit(x,"_")})
  
  ##browser()
  ##multiplet_index = grep("_",ab_droplet_ids);
  ##multiplet_array = sapply(ab_droplet_ids[multiplet_index], function(x) {strsplit(x,"_")})
  
  ## matrix of doublet rates
  ########################################################################################
  ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ## this code really only works for doublets, if lots of multiplets, what's calculated is
  ## the co-occurrence frequencies
  ## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ########################################################################################
  
  exp_doublet_rates <- matrix(NA, nrow=npools*2, ncol=npools*2);
  for(i in 1:(2*npools)) {
    for(j in i:(2*npools)) {
      exp_doublet_rates[i,j] <- exp_doublet_rates[j,i] <- ab_nsinglets[i]*ab_nsinglets[j]/sum(ab_nsinglets)^2*mr_ab;
     }
  }
  
  ## we need to normalize by the frequency that two anbibodies from the same
  ## pool cooccur
  obs_doublet_rates <- matrix(NA,nrow=npools*2,ncol=npools*2);
  for(i in 1:(2*npools)) {
    for(j in i:(2*npools)) {
      ## code to try to deal with droplets with > 2 cells
      if(i != j) {
        matched_multiplets <- unlist(sapply(multiplet_array, function(x) {
          if(length(na.omit(match(c(i,j),x)))==2) {
            ##return(1/length(x));
            ##return(1);
            return(1/perm(length(x),2))
            ##return(1/choose(length(x),2)/2)
          }
        }))
        
        obs_doublet_rates[i,j] <- obs_doublet_rates[j,i] <- sum(matched_multiplets)/ndrops
      } else {
        obs_doublet_rates[i,j] <- obs_doublet_rates[j,i] <- 0
      }
     
      # if(i+npools == j) {
      #   obs_doublet_rates[i,j] <- obs_doublet_rates[j,i] <- obs_doublet_rates[i,j]*exp_doublet_rates[i,j]
      # }
      
      ## old code for dealing with doublets
      ##obs_doublet_rates[i,j] <- obs_doublet_rates[j,i] <- length(grep(paste("^",i,"_",j,"$",sep=""),ab_droplet_ids))/ndrops/2
      
      ## new code with dealing with doublets, only count singlets and doublets
      ##ab_doublets <- which(sapply(multiplet_array, function(x) {length(x)==2}))
      ##obs_doublet_rates[i,j] <- obs_doublet_rates[j,i] <- length(grep(paste("^",i,"_",j,"$",sep=""),ab_droplet_ids))/(sum(ab_nsinglets)+length(ab_doublets))/2
    }
  }
  
  
  # sim_doublet_rates <- matrix(NA, nrow=npools*2, ncol=npools*2);
  # ## let's simulate
  # sim_ndrops <- 1e6;
  # sim_singlets <- unlist(sapply(1:(npools*2), function(i) {
  #   rep(i, ab_nsinglets[i]*sim_ndrops/ndrops)
  # }))
  # sim_doublets_first <- sample(sim_singlets, sim_ndrops*mr_ab, replace=T)
  # sim_doublets_second <- sample(sim_singlets, sim_ndrops*mr_ab, replace=T)
  # sim_doublet_ids <- sapply(1:length(sim_doublets_first), function(i) {
  #   paste(sort(c(sim_doublets_first[i],sim_doublets_second[i])),collapse="_");
  # })
  # for(i in 1:(2*npools)) {
  #   for(j in i:(2*npools)) {
  #     sim_doublet_rates[i,j] <- sim_doublet_rates[j,i] <- length(grep(paste("^",i,"_",j,"$",sep=""),sim_doublet_ids))/ndrops
  #     ##print(length(grep(paste("^",i,"_",j,"$",sep=""),ab_droplet_ids))/ndrops)
  #     ##print(obs_doublet_rates)
  #   }
  # }
  
  # remove the entries for same batch
  for(i in 1:npools) {
    exp_doublet_rates[i,i+npools] <- exp_doublet_rates[i+npools, i] <- NA
    obs_doublet_rates[i,i+npools] <- obs_doublet_rates[i+npools, i] <- NA
  }

  
  print(cor(exp_doublet_rates[upper.tri(exp_doublet_rates)], obs_doublet_rates[upper.tri(obs_doublet_rates)],use='complete.obs'))
  exp_obs_df <- data.frame(exp=exp_doublet_rates[upper.tri(exp_doublet_rates)], obs=obs_doublet_rates[upper.tri(obs_doublet_rates)])
  
  return(exp_obs_df);
  
  # cat("tx singlets:\n")
  # print(tx_nsinglets)
  # cat("ab singlets:\n");
  # print(ab_nsinglets)
  # cat("npools:",npools, "mr_ab_resolve:", mr_ab-sum(obs_doublet_rates), "mr_ab:",mr_ab," mr_tx:",mr_tx,'\n')
  
}