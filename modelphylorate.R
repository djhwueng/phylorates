rm(list=ls())
library(ape)
library(geiger)
library(extraDistr)
library(LaplacesDemon)

#################################################
#########Draw proposal from Prior ###############
#################################################

priorsample<-function(prmodel="halfnorm",prtausq="invgamma"){
  #change to density so you have prior
  if(prmodel=="halfnorm"){kappa<-LaplacesDemon::rhalfnorm(1,scale=sqrt(1000*pi/2));alpha<-LaplacesDemon::rhalfnorm(1,scale=sqrt(1000*pi/2))}
  if(prmodel=="hcauchy"){kappa<-rhcauchy(1, sigma = .0010);alpha<-rhcauchy(1, sigma = .0010)}
  if(prmodel=="halft"){kappa<-rhalft(1, scale=0.01, nu=10);alpha<-rhalft(1, scale=0.01, nu=10)}
  if(prmodel=="lnorm"){kappa<-rlnorm(1, meanlog = -5, sdlog = 0.2);alpha<-rlnorm(1, meanlog = -5, sdlog = 0.2)}
  if(prmodel=="invgamma"){kappa<-LaplacesDemon::rinvgamma(1, shape=2, scale=0.01);alpha<-LaplacesDemon::rinvgamma(1, shape=2, scale=0.01)}
  ## --------------------
  if(prtausq=="invgamma"){tausq<-LaplacesDemon::rinvgamma(1,shape=2, scale=0.01)}
  if(prtausq=="lnorm"){tausq<-rlnorm(1,meanlog = -5, sdlog = 0.2)}
  params<-c(kappa,alpha,tausq)
  names(params)<-c("kappa","alpha","tausq")
  return(params)
}

priorsample(prmodel="halfnorm",prtausq="invgamma")
priorsample(prmodel="hcauchy",prtausq="invgamma")
priorsample(prmodel="halft",prtausq="invgamma")
priorsample(prmodel="lnorm",prtausq="invgamma")
priorsample(prmodel="invgamma",prtausq="invgamma")

priorsample(prmodel="halfnorm",prtausq="lnorm")
priorsample(prmodel="hcauchy",prtausq="lnorm")
priorsample(prmodel="halft",prtausq="lnorm")
priorsample(prmodel="lnorm",prtausq="lnorm")
priorsample(prmodel="invgamma",prtausq="lnorm")

################################################
#             The prior function               #
################################################

priordensity<-function(params=params,prmodel="halfnorm",prtausq="invgamma"){
  kappa<-params["kappa"]# for Garch(1,1) just use kappa = 0
  alpha<-params["alpha"]
  tausq<-params["tausq"]
  #change to density so you have prior
  if(prmodel=="halfnorm"){prm.kappa<-LaplacesDemon::dhalfnorm(kappa,scale=sqrt(1000*pi/2),log=T);prm.alpha<-LaplacesDemon::dhalfnorm(alpha,scale=sqrt(1000*pi/2),log=T)}
  if(prmodel=="hcauchy"){prm.kappa<-dhcauchy(kappa,sigma = .001,log=T);prm.alpha<-dhcauchy(alpha,sigma = .001,log=T)}
  if(prmodel=="halft"){prm.kappa<-dhalft(kappa,scale=0.01, nu=10,log=T);prm.alpha<-dhalft(alpha,scale=0.01, nu=10,log=T)}
  if(prmodel=="lnorm"){prm.kappa<-dlnorm(kappa,meanlog = -5, sdlog = 0.2,log=T);prm.alpha<-dlnorm(alpha,meanlog = -5, sdlog = 0.2,log=T)}
  if(prmodel=="invgamma"){prm.kappa<-LaplacesDemon::dinvgamma(kappa,shape=2, scale=0.01,log=T);prm.alpha<-LaplacesDemon::dinvgamma(alpha,shape=2, scale=0.01,log=T)}
  ## --------------------
  if(prtausq=="invgamma"){prm.tausq<-LaplacesDemon::dinvgamma(tausq,shape=2, scale=0.01,log=T)}
  if(prtausq=="lnorm"){prm.tausq<-dlnorm(tausq,meanlog = -5, sdlog = 0.2,log=T)}
  return(list(prm.kappa=prm.kappa,prm.alpha=prm.alpha,prm.tausq=prm.tausq))
}


params<-priorsample(prmodel="halfnorm",prtausq="invgamma")
priordensity(params=params,prmodel="halfnorm",prtau="invgamma")
priordensity(params=params,prmodel="hcauchy",prtausq="invgamma")
priordensity(params=params,prmodel="halft",prtausq="invgamma")
priordensity(params=params,prmodel="lnorm",prtausq="invgamma")
priordensity(params=params,prmodel="invgamma",prtausq="invgamma")

priordensity(params=params,prmodel="halfnorm",prtausq="lnorm")
priordensity(params=params,prmodel="hcauchy",prtausq="lnorm")
priordensity(params=params,prmodel="halft",prtausq="lnorm")
priordensity(params=params,prmodel="lnorm",prtausq="lnorm")
priordensity(params=params,prmodel="invgamma",prtausq="lnorm")



################################################
#         The likelihood function              #
################################################

like10<-function(params=params,trait=trait,tree=tree){
  kappa<-params["kappa"]# for Garch(1,1) just use kappa = 0
  alpha<-params["alpha"]
  tausq<-params["tausq"]
  
  y<-trait
  makeL(tree)->L
  makeL1(tree)->L1
  
  optL <- function(lambda){
    #lambda<-0.001
    y <- scale(y)#(y-mean(y))/sd(y)
    betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%t(L)) %*% (as.matrix(y) - rootV)
    y.hat <- (L %*% betas) + rootV
    Rvar <- array()
    for (i in 1:Ntip(tree)) {
      #    i<-2
      ace.tip <- betas[match(names(which(L[i, ] != 0)),rownames(betas)), ]
      mat = as.matrix(dist(ace.tip))
      Rvar[i] <- sum(mat[row(mat) == col(mat) + 1])
    }
    #  the rate variation within clades
    abs(1 - (var(Rvar)))  #+ (mean(as.matrix(y))/mean(y.hat))) #y is scaled so zero mean
  }
  
  u <- data.frame(trait, (1/diag(vcv(tree))^2))
  u <- u[order(u[, ncol(u)], decreasing = TRUE), ]
  u1 <- u[1:(nrow(u) * 0.1), , drop = FALSE]
  rootV <- c(apply(u1[, 1:(ncol(u1) - 1), drop = FALSE],2, function(x) weighted.mean(x, u1[, dim(u1)[2]])))# differ 
  
  ntaxa<-length(y)
  tree<-reorder(tree,"postorder")
  tree$root.edge<-0
  N<-dim(tree$edge)[1]
  anc<-tree$edge[,1]
  des<-tree$edge[,2]
  treelength<-tree$edge.length
  
  h <- mle(optL, start = list(lambda = 1), method = "L-BFGS-B",upper = 10, lower = 0.001)
  lambda <- h@coef
  lambda
  betas.rrphylo <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%t(L)) %*% (as.matrix(y) - rootV)
  betas.rrphylo
  betas.rrphylo<-betas.rrphylo[c(ntaxa:(2*ntaxa-1),1:(ntaxa-1)),]
  betas.rrphylo
  
  
  betas.archrate<-array(NA,c(2*ntaxa-1))#this is Garch(1,0)
  names(betas.archrate)<-colnames(L)
  betas.archrate<-betas.archrate[c(ntaxa:(2*ntaxa-1),1:(ntaxa-1))]
  betas.archrate[names(betas.archrate)==ntaxa+1] <- betas.rrphylo[ntaxa+1]#sqrt(bm.fit$opt$sigsq)
  betas.archrate
  
  for(index in N:1){
    betas.archrate[des[index]]<- betas.rrphylo[des[index]] + alpha/treelength[index]*betas.archrate[anc[index]]
  }  
  betas.rrphylo<-betas.rrphylo[c((ntaxa+1):(2*ntaxa-1),1:ntaxa)]
  betas.archrate<-betas.archrate[c((ntaxa+1):(2*ntaxa-1),1:ntaxa)]
  
  #this is prediction
  trait
  y.hat.rrphylo <- (L %*% betas.rrphylo) + rootV
  ace.rrphylo <- (L1 %*% betas.rrphylo[1:Nnode(tree) ]) + rootV
  y.hat.archrate <- (L %*% betas.archrate) + rootV
  ace.archrate <- (L1 %*% betas.archrate[1:Nnode(tree) ]) + rootV
  ###
  
  loglike<--ntaxa/2*tausq-t(trait-L%*%betas.archrate)%*%(trait-L%*%betas.archrate)/tausq/2  
  return(loglike)
}


ntaxa<-5
tree<-rtree(ntaxa)
phen<-fastBM(tree,internal=T)
trait<-phen[1:ntaxa]
trait
like10(params=params,trait=trait,tree=tree)

################################################
#            The posterior distribution        #
################################################

postg10 <- function(params,trait=trait,tree=tree,prmodel="halfnorm",prtau="invgamma"){
  prden<-priordensity(params=params,prmodel="halfnorm",prtau="invgamma")
  return(like10(params=params,trait=trait,tree=tree) + prden$prm.kappa+prden$prm.alpha+prden$prm.tausq)
}



###################################################
#            main program                         #
###################################################
ntaxa<-5
rtree(ntaxa)->tree
fastBM(tree,internal=T)->phen
phen[1:ntaxa]->trait
trait

prmodel="halfnorm"
prtausq="invgamma"

postg10(params,trait=trait,tree=tree,prmodel="halfnorm",prtau="invgamma")



