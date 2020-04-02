cr_sim_singlet_multiplet <- function(cells_loaded, pools, crr) {
##print(cells_loaded)
##cells_loaded <- 200000
##pools = 10
##crr = 0.56


# sampling from poisson
cells_per_drop = rpois(100000, lambda=cells_loaded/100000);
##sim = sapply(sim, function(x) rbinom(1,x,prob=.6));
cells_per_drop <- cells_per_drop[sample(1:length(cells_per_drop), length(cells_per_drop)*crr)];
drops <- length(which(cells_per_drop>0)); 
cells <-  sum(cells_per_drop);
multiplets <- length(which(cells_per_drop>1))

## code to calculate collision rate as a function of pool barcodes

## droplet containing cells
cells_per_dcc = cells_per_drop[which(cells_per_drop>0)];

unique_barcodes_per_dcc <- NULL;
bc2_sim <- NULL;
bc0_sim <- NULL;
bc1_sim <- NULL;

genetic_singlets <- NULL;
genetic_doublets <- NULL;
ab_singlets <- NULL;
ab_multiplets <- NULL;
ab_cells <- NULL;

# cells_per_dcc - cells per droplet-containing-cells (dcc)
for(x in cells_per_dcc) {
  tmp <- sample(1:pools,x,replace=T)
  genetic <- sample(1:2,x,replace=T)
  unique_barcodes_per_dcc <- c(unique_barcodes_per_dcc,length(unique(tmp)));
  tab <- tabulate(tmp,nbins=pools);
  genetic_tab <- tabulate(genetic, nbins=2);

  bc2_sim <- c(bc2_sim, sum(tab>1))
  bc0_sim <- c(bc0_sim, sum(tab==0))
  bc1_sim <- c(bc1_sim, sum(tab==1))
  
  genetic_singlets <- c(genetic_singlets, sum(genetic_tab>0)==1)
  genetic_doublets <- c(genetic_doublets, sum(genetic_tab>0)==2)
  
  ab_singlets <- c(ab_singlets, sum(tab>0)==1)
  ab_multiplets <- c(ab_multiplets, sum(tab>0)>1)
  ab_cells <- c(ab_cells, sum(tab>0))
}

cr2_genetic_singlets <- sum(bc2_sim[which(genetic_singlets)])/sum(unique_barcodes_per_dcc[which(genetic_singlets)])
cr2_genetic_doublets <- sum(bc2_sim[which(genetic_doublets)])/sum(unique_barcodes_per_dcc[which(genetic_doublets)])

cr2_singlets <- sum(bc2_sim[which(ab_singlets)])/sum(unique_barcodes_per_dcc[which(ab_singlets)])
cr2_multiplets <- sum(bc2_sim[which(ab_multiplets)])/sum(unique_barcodes_per_dcc[which(ab_multiplets)])

cr2 <- sum(bc2_sim)/sum(unique_barcodes_per_dcc)

return(ab_cells);
}