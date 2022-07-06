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


################################################
#         The likelihood function              #
################################################
#here your can actually specify the three model garch01,garch11, and the null model garch00 

LikelihoodAll<-function(params=params,trait=trait,tree=tree,modeltype="garch00"){
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
  
  # compute the rootV
  u <- data.frame(trait, (1/diag(vcv(tree))^2))
  u <- u[order(u[, ncol(u)], decreasing = TRUE), ]
  u1 <- u[1:(nrow(u) * 0.1), , drop = FALSE]
  rootV <- c(apply(u1[, 1:(ncol(u1) - 1), drop = FALSE],2, function(x) weighted.mean(x, u1[, dim(u1)[2]])))# differ 
  
  ### make poster order for ease of computing the likelihood
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
  
  ##here is the null model likelihood, no prior for garch but only the tausq 
  betas.rrphylo <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%t(L)) %*% (as.matrix(y) - rootV)
  betas.rrphylo
  
  
  ###########################
  # garch00,garch10, garch11#
  ###########################
  
  if(modeltype=="garch00"){loglike<--ntaxa/2*tausq-t(trait-L%*%betas.rrphylo)%*%(trait-L%*%betas.rrphylo)/tausq/2
  }
  
  if(modeltype=="garch10"){
    ##reorder the taxa and node 
    betas.rrphylo<-betas.rrphylo[c(ntaxa:(2*ntaxa-1),1:(ntaxa-1)),]
    betas.archrate<-array(NA,c(2*ntaxa-1))#this is Garch(1,0)
    names(betas.archrate)<-colnames(L)
    betas.archrate<-betas.archrate[c(ntaxa:(2*ntaxa-1),1:(ntaxa-1))]
    betas.archrate[names(betas.archrate)==ntaxa+1] <- betas.rrphylo[ntaxa+1]
    for(index in N:1){
      betas.archrate[des[index]]<- betas.rrphylo[des[index]] + alpha/treelength[index]*betas.archrate[anc[index]]
    }  
    betas.archrate<-betas.archrate[c((ntaxa+1):(2*ntaxa-1),1:ntaxa)]
    loglike<--ntaxa/2*tausq-t(trait-L%*%betas.archrate)%*%(trait-L%*%betas.archrate)/tausq/2  
  }
  
  if(modeltype=="garch11"){
    ##reorder the taxa and node 
    betas.rrphylo<-betas.rrphylo[c(ntaxa:(2*ntaxa-1),1:(ntaxa-1)),]
    betas.garch11rate<-array(NA,c(2*ntaxa-1))#this is Garch(1,1)
    names(betas.garch11rate)<-colnames(L)
    betas.garch11rate<-betas.garch11rate[c(ntaxa:(2*ntaxa-1),1:(ntaxa-1))]
    betas.garch11rate[names(betas.garch11rate)==ntaxa+1] <- betas.rrphylo[ntaxa+1]
    for(index in N:1){
      betas.garch11rate[des[index]]<- betas.rrphylo[des[index]]+ kappa/treelength[index]*betas.garch11rate[anc[index]]*rnorm(1,0,treelength[index])^2  + alpha/treelength[index]*betas.garch11rate[anc[index]]
    }  
    betas.garch11rate<-betas.garch11rate[c((ntaxa+1):(2*ntaxa-1),1:ntaxa)]
    loglike<--ntaxa/2*tausq-t(trait-L%*%betas.garch11rate)%*%(trait-L%*%betas.garch11rate)/tausq/2  
  }
  return(loglike)
}

#this is prediction
# trait
# y.hat.rrphylo <- (L %*% betas.rrphylo) + rootV
# ace.rrphylo <- (L1 %*% betas.rrphylo[1:Nnode(tree) ]) + rootV
# y.hat.archrate <- (L %*% betas.archrate) + rootV
# ace.archrate <- (L1 %*% betas.archrate[1:Nnode(tree) ]) + rootV
###

################################################
#            The posterior distribution        #
################################################

posteriorAll <- function(params,trait=trait,tree=tree,prmodel="halfnorm",prtau="invgamma",modeltype="garch00"){
  prden<-priordensity(params=params,prmodel="halfnorm",prtau="invgamma")
  return(LikelihoodAll(params=params,trait=trait,tree=tree,modeltype=modeltype) + prden$prm.kappa+prden$prm.alpha+prden$prm.tausq)
}

################################################
#                 proposalfunction             #       
################################################
proposalfunction <- function(param){
  props<-rnorm(3,mean = param, sd= c(0.01,0.01,0.01))
  names(props)<-c("kappa","alpha","tausq")
  return(props)
}
