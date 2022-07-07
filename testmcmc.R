rm(list=ls())
source("~/Dropbox/Conference/1Prof/202207ECAISBA/rcode/modelphylorate.R")
###################################################
#            main program                         #
###################################################
# ntaxa<-5
# tree<-rtree(ntaxa)
# phen<-fastBM(tree,internal=T)
# trait<-phen[1:ntaxa]
# trait
# dim(trait)
# typeof(trait)
# tree$tip.label
###################################################
#            sim trait  using true parameter      #
###################################################

ntaxa<-5
kappa<-0.002
alpha<-0.006
tausq<-1
betas<-rnorm(2*ntaxa-1)
tree<-rtree(ntaxa)
tree1<-tree
makeL(tree)->L
makeL1(tree)->L1
names(betas)<-colnames(L)
trait<-array(NA,c(2*ntaxa-1)) #this is Garch(1,1)
names(trait)<-colnames(L)
betas<-betas[c(ntaxa:(2*ntaxa-1),1:(ntaxa-1))]
trait<-trait[c(ntaxa:(2*ntaxa-1),1:(ntaxa-1))]
trait[names(trait)==ntaxa+1] <- betas[ntaxa+1]
betas
trait
tree<-reorder(tree,"postorder")
tree$root.edge<-0
N<-dim(tree$edge)[1]
anc<-tree$edge[,1]
des<-tree$edge[,2]
treelength<-tree$edge.length

for(index in N:1){
  trait[des[index]]<- betas[des[index]]+ kappa/treelength[index]*betas[anc[index]]*rnorm(1,0,treelength[index])^2  + alpha/treelength[index]*betas[anc[index]] + rnorm(1,0,tausq)
}  
trait<-c(trait[1:ntaxa])
dim(trait)
trait
tree<-tree1#no postorder
tree$tip.label
plot(tree)
####

params<-priorsample(prmodel="halfnorm",prtausq="invgamma")
LikelihoodAll(params=params,trait=trait,tree=tree,modeltype="garch00")
LikelihoodAll(params=params,trait=trait,tree=tree,modeltype="garch10")
LikelihoodAll(params=params,trait=trait,tree=tree,modeltype="garch11")

posteriorAll(params,trait=trait,tree=tree,prmodel="halfnorm",prtau="invgamma",modeltype="garch00")
posteriorAll(params,trait=trait,tree=tree,prmodel="halfnorm",prtau="invgamma",modeltype="garch10")
posteriorAll(params,trait=trait,tree=tree,prmodel="halfnorm",prtau="invgamma",modeltype="garch11")

#reselect startvalues
startvalue<-c(0.1,0.1,1)
names(startvalue)<-c("kappa","alpha","tausq")
posteriorAll(startvalue,trait=trait,tree=tree,prmodel="halfnorm",prtau="invgamma",modeltype="garch11")

# next can do mcmc shall be work better now
#https://khayatrayen.github.io/MCMC.html

iter<-1e5
chain<-array(dim=c(iter+1, length(startvalue)))
colnames(chain)<-names(params)
chain[1,]<-startvalue
head(chain)
for (i in 1:iter){
  if(i%%1e3==0){print(i)}
#  i<-1
  #print(i)
  #Tony you shall use a symmetric proposal, normal distribution shall be nice that avoud negative. This is equal to 1 if the proposal density is symmetric. 
  proposal = proposalfunction(chain[i,])
  #print(round(proposal,4))
  new<-posteriorAll(proposal,trait=trait,tree=tree,prmodel="halfnorm",prtau="invgamma",modeltype="garch00")
  old<-posteriorAll(chain[i,],trait=trait,tree=tree,prmodel="halfnorm",prtau="invgamma",modeltype="garch00")
  probab= exp(  new  -  old   )
  if(runif(1)<probab){
    chain[i+1,]<-proposal
  }else{
    chain[i+1,]<-chain[i,]
  }
}

burnIn<-50000
par(mfrow=c(1,3))
plot(chain[-(1:burnIn),"kappa"],main="kappa",type="l")
plot(chain[-(1:burnIn),"alpha"],main="alpha",type="l")
plot(chain[-(1:burnIn),"tausq"],main="tausq",type="l")

acceptance=1-mean(duplicated(chain[-(1:burnIn),]))
acceptance

par(mfrow = c(2,3))
hist(chain[-(1:burnIn),"kappa"],nclass=30, , main="Posterior of kappa", xlab="True value = red line" )
abline(v = median(chain[-(1:burnIn),"kappa"]),col="blue",lwd=3)
abline(v = 0, col="red",lty=2 )

hist(chain[-(1:burnIn),"alpha"],nclass=30, main="Posterior of alpha", xlab="True value = red line")
abline(v = median(chain[-(1:burnIn),"alpha"]),col="blue",lwd=3)
abline(v = 0, col="red",lty=2 )

hist(chain[-(1:burnIn),"tausq"],nclass=30, main="Posterior of tausq", xlab="")
abline(v = median(chain[-(1:burnIn),"tausq"]) ,col="blue",lwd=2)
#abline(v = 0.5, col="red" )

plot(chain[-(1:burnIn),"kappa"], type = "l", xlab="True value = red line" , main = "Chain values of kappa", )
abline(h = 0, col="red" )

plot(chain[-(1:burnIn),"alpha"], type = "l", xlab="True value = red line" , main = "Chain values of alpha", )
abline(h = 0, col="red" )

plot(chain[-(1:burnIn),"tausq"], type = "l", xlab="True value = red line" , main = "Chain values of tausq", )
#abline(h = 0.5, col="red" )





