library(ggplot2);
library(plyr)
setwd("~/Box/Papers/In Preparation/2019.Hwang.scitoseq/Code")

## cost per read assuming $7083 per 3 billion reads
## cost.seq.base <- 7083/3e9
## assuming $7/gigabase and we read 100 bases
##cost.seq.per.read <- 5/1e9*100
cost.seq.per.read <- 28332/12e9;
## ucsf discount
cost.seq.per.read <- 2.4/1e6

## cost per well
cost.prep.per.well <- 1482;
## cost per antibody
cost.ab.per.ab <- 7.5;
## number of markers
n.markers <- 50;

## simulate cid-seq
## assume 500 (tn5) + 934 rt + 1150 master mix + 300 = 2884
## cost per drop = 657/50000+(300)/50000 = 0.03

## crr cell recovery rate
crr <- 0.6;

pools_list <- c(1,2,3,4,5,8,10,11,12,14,16,20,22,23,24,25,32,36,48,64,96,128,144,160,176,192,256,320, 352,384);

## estimates cells_loaded for a specific cell_recovered
cells_recovered <- function(cells_loaded, cells_recovered_target, crr=0.6) {
  lambda <- cells_loaded/100000;
  
  cells_per_drop <- rpois(100000,lambda);
  ## this is because the capture rate is ~50%
  cells_per_drop <- cells_per_drop[sample(1:length(cells_per_drop), length(cells_per_drop)*crr)];

  cells_recovered <- sum(cells_per_drop)
    
  return(abs(cells_recovered-cells_recovered_target));
}

## estimates cells_loaded for a specific cell_recovered
drops_recovered <- function(cells_loaded, drops_recovered_target, crr=0.6) {
  lambda <- cells_loaded/100000;
  
  cells_per_drop <- rpois(100000,lambda);
  ## this is because the capture rate is ~50%
  cells_per_drop <- cells_per_drop[sample(1:length(cells_per_drop), length(cells_per_drop)*crr)];
  
  drops_recovered <- sum(cells_per_drop>0)
  
  return(abs(drops_recovered-drops_recovered_target));
}

cr.func.exact <- function(cells_loaded, pools) {
  numerator <- 1-exp(-cells_loaded/(pools*100000))*(1+cells_loaded/(pools*100000))
  denominator <- 1-exp(-cells_loaded/(pools*100000))
  return(numerator/denominator)
  
}
#cr.func2(optimize(cells_recovered, interval=c(1e3,1e6), cells_recovered_target=35371)$minimum,1)
  
##lambda.cid <- 1.6;
##cr.func <- function(cells_loaded, lambda, pools, cr_target=0.05) {
cr.func2 <- function(cells_loaded, pools, crr) {
  ## cells_loaded
  lambda = cells_loaded/100000;
  ##cells_per_drop<-rpois(100000,lambda);
  cells_per_drop<-rpois(100000,lambda);
  ## this is because the capture rate is ~50%
  cells_per_drop <- cells_per_drop[sample(1:length(cells_per_drop), length(cells_per_drop)*crr)];
  
  ## number of pool barcodes not tagging a cell
  bc0 <- pools * (1-1/pools)^cells_per_drop;
  
  ## number of pool barcodes tagging exactly one cell
  bc1 <- cells_per_drop*(1-1/pools)^(cells_per_drop-1);
  
  ## number of pool barcodes tagging > 1 cell
  bc2 <- pools - bc0 - bc1;
  
  tmp <- bc2/(bc1+bc2)*cells_per_drop;
  
  cr <- sum(na.omit(tmp))/sum(cells_per_drop[!is.na(tmp)]);
  
  cr3 <- mean(bc2/(bc1+bc2),na.rm=T);
  
  ##return(abs(cr-cr_target))
  return(list(sum(cells_per_drop),sum(bc1+bc2),cr3))
}

cr.func <- function(cells_loaded, pools, cr_target=0.05) {
  ## cells_loaded
  lambda = cells_loaded/100000;
  ##cells_per_drop<-rpois(100000,lambda);
  cells_per_drop<-rpois(100000,lambda);
  ## this is because the capture rate is ~50%
  cells_per_drop <- cells_per_drop[sample(1:length(cells_per_drop), length(cells_per_drop)*crr)];
  
  ## number of pool barcodes not tagging a cell
  bc0 <- pools * (1-1/pools)^cells_per_drop;
  
  ## number of pool barcodes tagging exactly one cell
  bc1 <- cells_per_drop*(1-1/pools)^(cells_per_drop-1);
  
  ## number of pool barcodes tagging > 1 cell
  bc2 <- pools - bc0 - bc1;
  
  tmp <- bc2/(bc1+bc2)*cells_per_drop;
  
  cr <- sum(na.omit(tmp))/sum(cells_per_drop[!is.na(tmp)]);

  cr3 <- mean(bc2/(bc1+bc2),na.rm=T);
  
  ##return(abs(cr-cr_target))
  return(abs(cr3-cr_target))
}

sim <- function(cr_target = 0.05) {
  cells <- NULL;
  crs <- NULL;
  crs2 <- NULL;
  crs3 <- NULL;
  cost.prep <- NULL;
  cost.seq <- NULL;
  cost.ab <- NULL;
  singlets <- NULL;
  collisions <- NULL;
  lambdas <- NULL;
  empties2 <- NULL;
  singlets2 <- NULL;
  collisions2 <- NULL;
  drops <- NULL;
  
  ## simulate for collision rate of 0.05
  for(pools in pools_list) {
    ## use two different sets of intervals to avoid some weird boundary conditions
    if(pools < 10) {
      interval = c(0,10*100000);
    } else {
      interval = c(0,100*100000);
    }
    cells_loaded <- optimize(cr.func,interval=interval,pools=pools,tol=1e-7, cr_target = cr_target)$minimum
    lambda.cid <- cells_loaded/100000
    lambdas <- c(lambdas, lambda.cid);
    ##cells_per_drop<-rpois(cells_loaded/2,lambda.cid);
    cells_per_drop<-rpois(100000,lambda.cid);
    cells_per_drop <- cells_per_drop[sample(1:length(cells_per_drop), length(cells_per_drop)*crr)];
    
    cells <- c(cells, sum(cells_per_drop));
    bc0 <- pools * (1-1/pools)^cells_per_drop;
    bc1 <- cells_per_drop*(1-1/pools)^(cells_per_drop-1);
    bc2 <- pools - bc0 - bc1;
    tmp <- bc2/(bc1+bc2)*cells_per_drop;
    cr <- sum(na.omit(tmp))/sum(cells_per_drop[!is.na(tmp)]);
    crs <- c(crs, cr);
    
    ## definition from vijay ramani
    crs3 <- c(crs3, mean(bc2/(bc1+bc2),na.rm=T))
  
    singlet <- sum(cells_per_drop)*(1-cr)
    collision <- sum(cells_per_drop)*cr
    singlets <- c(singlets, singlet);
    collisions <- c(collisions, collision);
    
    empty2 <-  length(which(cells_per_drop==0));
    singlet2 <- length(which(cells_per_drop==1));
    collision2 <- length(which(cells_per_drop>1));
    empties2 <-  c(empties2, empty2);
    singlets2 <- c(singlets2, singlet2);
    collisions2 <- c(collisions2, collision2);

    crs2 <- c(crs2, collision2/(collision2+singlet2))
    
    drops <- c(drops, length(which(cells_per_drop>0)))
  }

  cost.prep <- c(cost.prep, cost.prep.per.well/singlets);##*100000
  ## assuming 10 reads per marker
  cost.seq <- rep(c(cost.seq, 40*n.markers*cost.seq.per.read)/(1-cr_target))##/(1-cr))*100000,length(singlets));
  ## assuming 50 abs, $1k per antibody and each conjugation can be used for 100 reactions
  cost.ab <- c(cost.ab, cost.ab.per.ab*n.markers/singlets);
    
  ## plot everything
  plot.cost.df <- rbind(data.frame(cost=cost.prep+cost.seq+cost.ab, pools=pools_list, method="scito-seq",type="total",lambda=lambdas),
                   data.frame(cost=cost.prep, pools=pools_list, method="scito-seq", type="prep",lambda=lambdas),
                   data.frame(cost=cost.seq, pools=pools_list, method="scito-seq", type="seq",lambda=lambdas),
                   data.frame(cost=cost.ab, pools=pools_list, method="scito-seq", type="ab",lambda=lambdas))
  
  plot.cells.df <- rbind(data.frame(drops=drops, cells=cells, singlets=singlets, collisions=collisions, empties2=empties2, singlets2=singlets2, collisions2=collisions2, cr=crs, cr2=crs2, cr3 = crs3, pools=pools_list, method="scito-seq"))
  
  
  return(list(plot.cost.df, plot.cells.df));
}

sim.5 <- sim(0.05);
sim.1 <- sim(0.01);
sim.10 <- sim(0.10);
sim.20 <- sim(0.20);

## plot cost vs pools
plot.cost.final.df <- rbind(##data.frame(sim.1[[1]],cr=0.01),
                            data.frame(sim.5[[1]],cr_target=0.05))#,
                            ##data.frame(sim.10[[1]],cr=0.1))
##ggsave("cost.vs.pools.pdf", ggplot(aes(pools, cost, linetype=type, color=as.factor(cr)),data=plot.cost.final.df)+geom_line()+theme_bw(),width=3,height=3);##+theme(legend.position='none')
##ggsave("cost.vs.pools.log.pdf", ggplot(aes(pools, cost, linetype=type, color=as.factor(cr)),data=plot.cost.final.df)+geom_line()+theme_bw()+scale_y_log10(breaks=c(1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1)),width=3,height=3);##+theme(legend.position='none')
ggsave("cost.vs.pools.5.pdf", ggplot(aes(pools, cost, color=as.factor(type)),data=plot.cost.final.df)+geom_line()+theme_bw()+theme(legend.position='none'),width=3,height=3)
plot.cost.final.df$type <- factor(plot.cost.final.df$type, levels=c("seq","prep","ab","total"))
plot.cost.final.df <- plot.cost.final.df[order(plot.cost.final.df$type),]
ggsave("cost.vs.pools.log.5.pdf", ggplot(aes(pools, cost, color=type),data=plot.cost.final.df)+geom_line()+theme_bw()+theme(legend.position='none')+scale_color_manual(values=c("#7CAE00","#F8766D","#00BFC4","#C77CFF"))+scale_x_log10(breaks=c(1,2,6,11,24,48,96,160,384))+scale_y_log10(breaks=c(1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1)),width=3,height=3);##

ggsave("prep.cost.vs.pools.5.pdf", ggplot(aes(pools, cost, color=as.factor(type)),data=plot.cost.final.df[plot.cost.final.df$type!="total" & plot.cost.final.df$type!="seq",])+geom_line()+theme_bw()+theme(legend.position='none'),width=3,height=3)
ggsave("prep.cost.vs.pools.log.5.pdf", ggplot(aes(pools, cost, color=as.factor(type)),data=plot.cost.final.df[plot.cost.final.df$type!="total" & plot.cost.final.df$type!="seq",])+geom_line()+theme_bw()+theme(legend.position='none')+scale_x_log10(breaks=c(1,2,6,11,24,48,96,160,384))+scale_y_log10(breaks=c(1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1)),width=3,height=3);##+theme(legend.position='none')


## plot cells vs pools
plot.cells.final.df <- rbind(data.frame(sim.1[[2]],cr_target=0.01),
                             data.frame(sim.5[[2]],cr_target=0.05),
                             data.frame(sim.10[[2]],cr_target=0.1))
ggsave("cells.vs.pools.pdf", ggplot(aes(pools, cells, color=as.factor(cr_target)),data=plot.cells.final.df)+geom_line()+theme_bw()+theme(legend.position='none'),width=3,height=3);
ggsave("cells.vs.pools.log.pdf", ggplot(aes(pools, cells, color=as.factor(cr_target)),data=plot.cells.final.df)+geom_line()+theme_bw()+theme(legend.position='none')+scale_x_log10(breaks=c(1,2,6,11,24,48,96,160,384))+scale_y_log10(breaks=c(1000,2000,5000,10000,20000,50000,100000,200000,500000,1e6,2e6, 5e6)),width=3,height=3);
##ggsave("cells.vs.pools.5.pdf", ggplot(aes(pools, cells, color=as.factor(cr)),data=plot.cells.final.df)+geom_line()+theme_bw()+theme(legend.position='none'),width=3,height=3);
##ggsave("cells.vs.pools.log.5.pdf", ggplot(aes(pools, cells, color=as.factor(cr)),data=plot.cells.final.df)+geom_line()+theme_bw()+theme(legend.position='none')+scale_y_log10(breaks=c(1000,2000,5000,10000,20000,50000,100000,200000,500000,1e6)),width=3,height=3);

## plot collision rates vs cells loaded

# ## let's first simulate pretending we have only one batch and set batch = 0
# ## this should match batch = 1 results below
# drops_10x <- NULL; cells_10x <- NULL; crs_10x <- NULL;
# for (lambda in seq(0.1,3,0.1)) {
#   cells_per_drop = rpois(100000, lambda=lambda);
#   ##sim = sapply(sim, function(x) rbinom(1,x,prob=.6));
#   cells_per_drop <- cells_per_drop[sample(1:length(cells_per_drop), length(cells_per_drop)*.6)];
#   drops_10x <<- c(drops_10x,length(which(cells_per_drop>0)));
#   cells_10x <<- c(cells_10x, sum(cells_per_drop));
#   crs_10x <<- c(crs_10x,length(which(cells_per_drop>1))/length(which(cells_per_drop>0)))
# }
# 
# plot.cells_loaded.crs.df <- data.frame(drops=drops_10x, cells_loaded=cells_10x/.6, cells_recovered=cells_10x, crs=crs_10x, type="expected", pools=0);
plot.cells_loaded.crs.df <- NULL;

## now let's simulate multiple pools
for (pools in c(1,2,5,10,20,48,96,384)) {
  drops <- NULL; multiplets <- NULL; cells <- NULL; crs2 <- NULL;
  crs <- NULL; crs3 <- NULL; cells_loaded_list <- 10^seq(3, 6, 0.5);

  for (cells_loaded in cells_loaded_list) {
    print(cells_loaded)
    cells_per_drop = rpois(100000, lambda=cells_loaded/100000);
    ##sim = sapply(sim, function(x) rbinom(1,x,prob=.6));
    cells_per_drop <- cells_per_drop[sample(1:length(cells_per_drop), length(cells_per_drop)*crr)];
    drops <<- c(drops,length(which(cells_per_drop>0))); 
    cells <<- c(cells, sum(cells_per_drop));
    multiplets <<- c(multiplets, length(which(cells_per_drop>1)))
  
    ## code to calculate collision rate as a function of pool barcodes
    
    ## droplet containing cells
    cells_per_dcc = cells_per_drop[which(cells_per_drop>0)];
    
    ## number of pool barcodes not tagging a cell in a droplet
    bc0 <- pools * (1-1/pools)^cells_per_dcc;
  
    ## number of pool barcodes tagging exactly one cell in a dropet
    ## bc1 <- cells_per_dcc*(1-1/pools)^(cells_per_dcc-1);
    bc1 <- cells_per_dcc * (1-1/pools)^(cells_per_dcc-1);
    
    ## number of pool barcodes tagging > 1 cell in a droplet
    bc2 <- pools - bc0 - bc1;
    
    ##tmp <- bc2/(bc1+bc2)*cells_per_dcc;
    ##cr <- sum(na.omit(tmp))/sum(cells_per_dcc[!is.na(tmp)]);

    ##cr3 <- mean(bc2/(bc1+bc2),na.rm=T);
    cr <- sum(bc2)/sum(bc1+bc2);
    
    crs <<- c(crs,cr)
    
    unique_barcodes_per_dcc <- NULL;
    bc2_sim <- NULL;
    bc0_sim <- NULL;
    bc1_sim <- NULL;

    # for(x in cells_per_dcc) {
    #   bc2_sim_iter <- NULL;
    #   bc0_sim_iter <- NULL;
    #   bc1_sim_iter <- NULL;
    #   
    #   for(iter in 1:100) {
    #     ##print(x)
    #     tmp <- sample(1:pools,x,replace=T)
    #     tab <- tabulate(tmp,nbins=pools);
    #     
    #     unique_barcodes_per_dcc <- c(unique_barcodes_per_dcc,length(unique(tmp)));
    #     bc2_sim_iter <- c(bc2_sim_iter, sum(tab>1))
    #     bc0_sim_iter <- c(bc0_sim_iter, sum(tab==0))
    #     bc1_sim_iter <- c(bc1_sim_iter, sum(tab==1))
    #   }
    #   bc2_sim <- c(bc2_sim, mean(bc2_sim_iter));
    #   bc0_sim <- c(bc0_sim, mean(bc0_sim_iter));
    #   bc1_sim <- c(bc1_sim, mean(bc1_sim_iter));
    #   
    #   
    # }

    for(x in cells_per_dcc) {
        tmp <- sample(1:pools,x,replace=T)
        unique_barcodes_per_dcc <- c(unique_barcodes_per_dcc,length(unique(tmp)));
        tab <- tabulate(tmp,nbins=pools);
        bc2_sim <- c(bc2_sim, sum(tab>1))
        bc0_sim <- c(bc0_sim, sum(tab==0))
        bc1_sim <- c(bc1_sim, sum(tab==1))
    }

    cr2 <- sum(bc2_sim)/sum(unique_barcodes_per_dcc)
    cr3 <- sum(unique_barcodes_per_dcc!=cells_per_dcc)/length(cells_per_dcc)
    
    crs2 <<- c(crs2, cr2)
    crs3 <<- c(crs3, cr3)
  }

  plot.cells_loaded.crs.df <- rbind(plot.cells_loaded.crs.df, 
                                    data.frame(drops=drops,
                                               cells_loaded=cells_loaded_list, 
                                               cells_recovered = cells,
                                               multiplets = multiplets,
                                               singlets=drops-multiplets,
                                               empties=100000*crr-drops,
                                               crs=crs,
                                               crs_drops = crs3,
                                               pools=pools,
                                               type="expected"),
                                    data.frame(drops=drops,
                                               cells_loaded=cells_loaded_list, 
                                               cells_recovered = cells,
                                               multiplets=multiplets,
                                               singlets=drops-multiplets,
                                               empties=100000*crr-drops,
                                               crs=crs2,
                                               crs_drops = crs3,
                                               pools=pools,
                                               type="simulated"));
                                    ##data.frame(drops=drops,cells_loaded=cells_loaded_list, cells_recovered = cells,crs=crs,pools=pools,type="expected"));
}


##observed_drops <- c(8765,9067,7357,17730,34549,38631,41391,55404);
##observed_drs <- c(.1409,.1777,.1049,0.1679,0.3804,0.4552,0.483438429,0.807838784)

observed_drops <- c(8765,7357,17730,34549,38674)##,41391,59433);

## add the 100k and 200k t/b experiments
observed_drops <- c(observed_drops, c(42386, 59433))
##observed_crs <- c(.1409,.1049,0.1679,0.3804,0.4552,0.483438429,0.807838784)
##observed_crs <- c(.0845,.0837,0.1679,0.3804,0.4552,0.483438429,0.807838784)

## this is based on tx
observed_crs <- c(.0906,.1064,0.1679,0.3804,0.4552,0.483438429,0.8403718)
## this is based on ab
observed_crs <- c(.0819,.0883,0.1717,0.3930,0.3865)

## add the 100k and 200k t/b experiment
observed_crs <- c(observed_crs, c(0.483438429,0.8403718))

## this is the 2 batch results
observed_drops <- c(observed_drops, 8765,7357,17730, 34549, 38674);
observed_crs <- c(observed_crs, .0250,.0270, 0.06, 0.14, 0.045)

## add in the t/b experiments
observed_drops <- c(observed_drops, c(42386, 59433));
observed_crs <- c(observed_crs, c(0.101,0.15)) ## not getting 6.3% yet 

## drops to cells conversion
observed_cells_loaded_pred <- NULL;
observed_drops_pred <- NULL;
observed_cells_pred <- NULL;
observed_multiplets_pred <- NULL;

for(observed_cr in observed_crs) {
  observed_cells_loaded_pred_i <- optimize(cr.func,interval=c(0,10*100000),pools=1,tol=1e-7, cr_target = observed_cr)$minimum;
  observed_cells_loaded_pred <<- c(observed_cells_loaded_pred,observed_cells_loaded_pred_i);
  
  observed_cells_per_drop_pred_i = rpois(100000, lambda=observed_cells_loaded_pred_i/100000);
  observed_cells_per_drop_pred_i <- observed_cells_per_drop_pred_i[sample(1:length(observed_cells_per_drop_pred_i), length(observed_cells_per_drop_pred_i)*crr)];
  observed_drops_pred <<- c(observed_drops_pred,length(which(observed_cells_per_drop_pred_i>0))); 
  observed_cells_pred <<- c(observed_cells_pred, sum(observed_cells_per_drop_pred_i));
  observed_multiplets_pred <<- c(observed_multiplets_pred,length(which(observed_cells_per_drop_pred_i>1))); 
  
}

plot.cells_loaded.crs.df <- rbind(data.frame(drops=observed_drops,
                                             cells_loaded=observed_cells_loaded_pred,
                                             cells_recovered=observed_cells_pred,
                                             multiplets=observed_multiplets_pred,
                                             singlets=observed_drops-observed_multiplets_pred,
                                             empties=100000*crr-observed_drops,
                                             crs=observed_crs,
                                             crs_drops=observed_crs,type="observed",pools=1),
             plot.cells_loaded.crs.df)

ggsave("cells_loaded.vs.crs.pdf",ggplot(data=plot.cells_loaded.crs.df, aes(cells_loaded, crs, color=as.factor(pools)))+
  geom_point(aes(alpha=as.factor(type)))+geom_line(aes(linetype=as.factor(type)),position = position_dodge(width = 0.05))+
  scale_alpha_manual(values=c(1,0,0))+scale_linetype_manual(values=c("blank","dashed","solid"))+
  theme_bw()+scale_x_log10(breaks=c(1e3,2e3,5e3,1e4,2e4,5e4,1e5,2e5,5e5,1e6))+
  scale_y_log10(breaks=c(1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1)), width=7, height=7)

plot.cells_loaded.crs.df <- NULL;

## now let's simulate multiple pools
for (pools in c(1,2,5,10,20,48,96,384)) {
  drops <- NULL; multiplets <- NULL; cells <- NULL; crs2 <- NULL;
  crs <- NULL; crs3 <- NULL; cells_loaded_list <- 10^seq(3, 6, 0.5);
  
  for (cells_loaded in cells_loaded_list) {
    print(cells_loaded)
    cells_per_drop = rpois(100000, lambda=cells_loaded/100000);
    ##sim = sapply(sim, function(x) rbinom(1,x,prob=.6));
    cells_per_drop <- cells_per_drop[sample(1:length(cells_per_drop), length(cells_per_drop)*crr)];
    drops <<- c(drops,length(which(cells_per_drop>0))); 
    cells <<- c(cells, sum(cells_per_drop));
    multiplets <<- c(multiplets, length(which(cells_per_drop>1)))
    
    ## code to calculate collision rate as a function of pool barcodes
    
    ## droplet containing cells
    cells_per_dcc = cells_per_drop[which(cells_per_drop>0)];
    
    ## number of pool barcodes not tagging a cell in a droplet
    bc0 <- pools * (1-1/pools)^cells_per_dcc;
    
    ## number of pool barcodes tagging exactly one cell in a dropet
    ## bc1 <- cells_per_dcc*(1-1/pools)^(cells_per_dcc-1);
    bc1 <- cells_per_dcc * (1-1/pools)^(cells_per_dcc-1);
    
    ## number of pool barcodes tagging > 1 cell in a droplet
    bc2 <- pools - bc0 - bc1;
    
    ##tmp <- bc2/(bc1+bc2)*cells_per_dcc;
    ##cr <- sum(na.omit(tmp))/sum(cells_per_dcc[!is.na(tmp)]);
    
    ##cr3 <- mean(bc2/(bc1+bc2),na.rm=T);
    cr <- sum(bc2)/sum(bc1+bc2);
    
    crs <<- c(crs,cr)
    
    unique_barcodes_per_dcc <- NULL;
    bc2_sim <- NULL;
    bc0_sim <- NULL;
    bc1_sim <- NULL;
    
    # for(x in cells_per_dcc) {
    #   bc2_sim_iter <- NULL;
    #   bc0_sim_iter <- NULL;
    #   bc1_sim_iter <- NULL;
    #   
    #   for(iter in 1:100) {
    #     ##print(x)
    #     tmp <- sample(1:pools,x,replace=T)
    #     tab <- tabulate(tmp,nbins=pools);
    #     
    #     unique_barcodes_per_dcc <- c(unique_barcodes_per_dcc,length(unique(tmp)));
    #     bc2_sim_iter <- c(bc2_sim_iter, sum(tab>1))
    #     bc0_sim_iter <- c(bc0_sim_iter, sum(tab==0))
    #     bc1_sim_iter <- c(bc1_sim_iter, sum(tab==1))
    #   }
    #   bc2_sim <- c(bc2_sim, mean(bc2_sim_iter));
    #   bc0_sim <- c(bc0_sim, mean(bc0_sim_iter));
    #   bc1_sim <- c(bc1_sim, mean(bc1_sim_iter));
    #   
    #   
    # }
    
    for(x in cells_per_dcc) {
      tmp <- sample(1:pools,x,replace=T)
      unique_barcodes_per_dcc <- c(unique_barcodes_per_dcc,length(unique(tmp)));
      tab <- tabulate(tmp,nbins=pools);
      bc2_sim <- c(bc2_sim, sum(tab>1))
      bc0_sim <- c(bc0_sim, sum(tab==0))
      bc1_sim <- c(bc1_sim, sum(tab==1))
    }
    
    cr2 <- sum(bc2_sim)/sum(unique_barcodes_per_dcc)
    cr3 <- sum(unique_barcodes_per_dcc!=cells_per_dcc)/length(cells_per_dcc)
    
    crs2 <<- c(crs2, cr2)
    crs3 <<- c(crs3, cr3)
  }
  
  plot.cells_loaded.crs.df <- rbind(plot.cells_loaded.crs.df, 
                                    data.frame(drops=drops,
                                               cells_loaded=cells_loaded_list, 
                                               cells_recovered = cells,
                                               multiplets = multiplets,
                                               singlets=drops-multiplets,
                                               empties=100000*crr-drops,
                                               crs=crs,
                                               crs_drops = crs3,
                                               pools=pools,
                                               type="expected"),
                                    data.frame(drops=drops,
                                               cells_loaded=cells_loaded_list, 
                                               cells_recovered = cells,
                                               multiplets=multiplets,
                                               singlets=drops-multiplets,
                                               empties=100000*crr-drops,
                                               crs=crs2,
                                               crs_drops = crs3,
                                               pools=pools,
                                               type="simulated"));
  ##data.frame(drops=drops,cells_loaded=cells_loaded_list, cells_recovered = cells,crs=crs,pools=pools,type="expected"));
}

plot.cells_loaded.crs.expected.df = plot.cells_loaded.crs.df[(plot.cells_loaded.crs.df$type=="expected") & (plot.cells_loaded.crs.df$batch==1),]

plot.cells_loaded.drops.df <- rbind(data.frame(cells_loaded=plot.cells_loaded.crs.expected.df$cells_loaded, drops=plot.cells_loaded.crs.expected.df$multiplets/(100000*crr), type="multiplet"),
                                    data.frame(cells_loaded=plot.cells_loaded.crs.expected.df$cells_loaded, drops=plot.cells_loaded.crs.expected.df$singlets/(100000*crr), type="singlet"),
                                    data.frame(cells_loaded=plot.cells_loaded.crs.expected.df$cells_loaded, drops=plot.cells_loaded.crs.expected.df$empties/(100000*crr), type="empties"))


ggsave("cells_loaded.vs.drops.pdf", ggplot(data=plot.cells_loaded.drops.df, aes(cells_loaded, drops, color=as.factor(type)))+
         geom_point()+geom_line()+theme_bw()+scale_x_log10(breaks=c(1e3,2e3,5e3,1e4,2e4,5e4,1e5,2e5,5e5,1e6)), width=7, height=7)


## plot density as a function of loading concentrations and lambda
cells <- NULL;
crs <- NULL;
crs2 <- NULL;
crs3 <- NULL;
cost.prep <- NULL;
cost.seq <- NULL;
cost.ab <- NULL;
singlets <- NULL;
collisions <- NULL;
lambdas <- NULL;
empties2 <- NULL;
singlets2 <- NULL;
collisions2 <- NULL;
drops <- NULL;

plot.cells_loaded.dist.df <- NULL

for(cells_loaded in c(1e3,2e3,5e3,1e4,2e4,5e4,1e5,1.81e5,2e5,5e5,1e6)) {
  lambda.cid <- cells_loaded/100000
  lambdas <- c(lambdas, lambda.cid);
  ##cells_per_drop<-rpois(cells_loaded/2,lambda.cid);
  cells_per_drop<-rpois(100000,lambda.cid);
  cells_per_drop <- cells_per_drop[sample(1:length(cells_per_drop), length(cells_per_drop)*crr)];

  cells <- c(cells, sum(cells_per_drop));
  bc0 <- pools * (1-1/pools)^cells_per_drop;
  bc1 <- cells_per_drop*(1-1/pools)^(cells_per_drop-1);
  bc2 <- pools - bc0 - bc1;
  tmp <- bc2/(bc1+bc2)*cells_per_drop;
  cr <- sum(na.omit(tmp))/sum(cells_per_drop[!is.na(tmp)]);
  crs <- c(crs, cr);
  
  ## definition from vijay ramani
  crs3 <- c(crs3, mean(bc2/(bc1+bc2),na.rm=T))
  
  singlet <- sum(cells_per_drop)*(1-cr)
  collision <- sum(cells_per_drop)*cr
  singlets <- c(singlets, singlet);
  collisions <- c(collisions, collision);
  
  empty2 <-  length(which(cells_per_drop==0));
  singlet2 <- length(which(cells_per_drop==1));
  collision2 <- length(which(cells_per_drop>1));
  empties2 <-  c(empties2, empty2);
  singlets2 <- c(singlets2, singlet2);
  collisions2 <- c(collisions2, collision2);
  
  plot.cells_loaded.dist.df <- rbind(plot.cells_loaded.dist.df,
                                     data.frame(cells_per_drop=cells_per_drop,
                                                cells_loaded=cells_loaded))
                                     
                                     
}

ggsave("cells.per.drop.dist.pdf", ggplot(data=plot.cells_loaded.dist.df, aes(cells_per_drop,fill=factor(cells_loaded)))+
  scale_y_log10()+
  scale_x_continuous(breaks=1:25)+
  geom_histogram(bins=20)+
  facet_grid(factor(cells_loaded)~.)+
  theme_bw(),width=7,height=7)
