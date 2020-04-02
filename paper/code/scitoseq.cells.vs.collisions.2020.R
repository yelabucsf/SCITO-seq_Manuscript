library(ggplot2);

## simulate scito-seq
scitoseq <- function(w1) {
  
  cells <- NULL;
  ecrs <- NULL;
  acrs <- NULL;
  drops <- NULL;
  
  lambdas <- c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.5,1.9,3,10)
  for(lambda in lambdas) {
    n2<-rpois(100000,lambda*1.9);
    drops <- c(drops, sum(n2>0))
    cells <- c(cells, sum(n2));
    bc0 <- w1 * (1-1/w1)^n2;
    bc1 <- n2*(1-1/w1)^(n2-1);
    bc2 <- w1 - bc0 - bc1;
    tmp <- bc2/(bc1+bc2)*n2;
    ecr <- sum(na.omit(tmp))/sum(n2[!is.na(tmp)]);
    ecrs <- c(ecrs, ecr);
    acr <- sum(which(n2>1))/sum(which(n2>0))
    acrs <- c(acrs, acr);
  }
  return(list(drops=drops, cells=cells, ecrs=ecrs, acrs=acrs, lambdas = lambdas))
}

scitoseq.1.rst <- scitoseq(1);
scitoseq.2.rst <- scitoseq(2);
scitoseq.4.rst <- scitoseq(4);
scitoseq.8.rst <- scitoseq(8);
scitoseq.16.rst <- scitoseq(16);
scitoseq.32.rst <- scitoseq(32);
scitoseq.48.rst <- scitoseq(48);
scitoseq.96.rst <- scitoseq(96);
##scitoseq.384.rst <- scitoseq(384);

dropseq <- function() {
  ## simluate drop-seq
  cells <- NULL;
  crs <- NULL;
  
  cells2 <- NULL;
  crs2 <- NULL;
  
  lambdas <-  c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,1.1,1.2);
  for(lambda in lambdas) {
    n2<-rpois(100000,lambda*1.9);
    
    # for(lambda in c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,0.2,0.3,0.4)) {
    #     n2<-rpois(50000,lambda);
    
    cells <- c(cells, length(which(n2>0)));
    cr <- length(which(n2>1))/length(which(n2>0));
    crs <- c(crs, cr);
    
    cells2 <- c(cells2, sum(n2))
    cr2 <- sum(n2[which(n2>1)])/sum(n2);
    crs2 <- c(crs2, cr2);
    
  }
  return(list(cells=cells, cells2=cells2, crs=crs, crs2=crs2))
}

dropseq.rst <- dropseq();

## plot everything
plot.df <- rbind(data.frame(cells=scitoseq.1.rst$cells, collision.rate=scitoseq.1.rst$ecrs, method="dsc-seq"),
  data.frame(cells=scitoseq.2.rst$cells, collision.rate=scitoseq.2.rst$ecrs, method="scito-seq.2"),
  data.frame(cells=scitoseq.4.rst$cells, collision.rate=scitoseq.4.rst$ecrs, method="scito-seq.4"),
  data.frame(cells=scitoseq.8.rst$cells, collision.rate=scitoseq.8.rst$ecrs, method="scito-seq.8"),
  data.frame(cells=scitoseq.16.rst$cells, collision.rate=scitoseq.16.rst$ecrs, method="scito-seq.16"),
  data.frame(cells=scitoseq.32.rst$cells, collision.rate=scitoseq.32.rst$ecrs, method="scito-seq.32"),
  data.frame(cells=scitoseq.48.rst$cells, collision.rate=scitoseq.48.rst$ecrs, method="scito-seq.48"),
  data.frame(cells=scitoseq.96.rst$cells, collision.rate=scitoseq.96.rst$ecrs, method="scito-seq.96")
)

ggsave("collision.vs.cells.pdf",width=5, height=5, ggplot(aes(collision.rate, cells, color=method),data=plot.df)+geom_point()+geom_smooth(method=glm, formula=y~x, family = gaussian(link = 'log'), se=FALSE)+scale_y_continuous(limits=c(1000,1000000), trans="log2")+scale_x_continuous(limits=c(1e-5,0.1),breaks=seq(0,0.1,0.01),trans="log10")+theme_bw())
#+scale_y_log10(breaks=seq(0,0.2,0.02));
