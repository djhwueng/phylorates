# https://cran.r-project.org/web/packages/RRphylo/vignettes/RRphylo.html
rm(list=ls())
library(geiger)
#optL, which is written as to minimize the rate variation within clades, thereby acting conservatively in terms of the chance to find rate shifts and introducing phylogenetic autocorrelation in evolutionary rates (Sakamoto & Venditti, Eastmann et al. 2011) .
optL <- function(lambda){
  #lambda<-0.001
  y <- scale(y)#(y-mean(y))/sd(y)
  betas <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%t(L)) %*% (as.matrix(y) - rootV)
  y.hat <- (L %*% betas) + rootV
  Rvar <- array()
  for (i in 1:Ntip(t)) {
    #    i<-2
    ace.tip <- betas[match(names(which(L[i, ] != 0)),rownames(betas)), ]
    mat = as.matrix(dist(ace.tip))
    Rvar[i] <- sum(mat[row(mat) == col(mat) + 1])
  }
  #  the rate variation within clades
  abs(1 - (var(Rvar)))  #+ (mean(as.matrix(y))/mean(y.hat))) #y is scaled so zero mean
}

ntaxa<-30
rtree(ntaxa)->tree
fastBM(tree,internal=T)->phen
phen[1:ntaxa]->y
y

makeL(tree)->L
makeL1(tree)->L1

#find root
t <- tree
toriginal <- t
yoriginal <- y <- treedataMatch(tree, y)[[1]]
Loriginal <- L <- makeL(t)
L1original <- L1 <- makeL1(t)
u <- data.frame(yoriginal, (1/diag(vcv(toriginal))^2))
u <- u[order(u[, ncol(u)], decreasing = TRUE), ]
u1 <- u[1:(nrow(u) * 0.1), , drop = FALSE]
rootV <- c(apply(u1[, 1:(ncol(u1) - 1), drop = FALSE],2, function(x) weighted.mean(x, u1[, dim(u1)[2]])))# differ from mle, i think they use weighted average
#rep(1,5)%*%solve(vcv(tree))%*%yoriginal/rep(1,5)%*%solve(vcv(tree))%*%matrix(rep(1,5),ncol=1)

tree<-reorder(tree,"postorder")
tree$root.edge<-0
N<-dim(tree$edge)[1]
anc<-tree$edge[,1]
des<-tree$edge[,2]
treelength<-tree$edge.length

## only you figure out the mechanism of the arch and garch 
## So you put the prior for the two 
## Now this seems work better it just need to make the tree  you put small alpha
# archratephy<-function(alpha,tree=tree,yoriginal=yoriginal,rootV=rootV){

#install.packages("fdrtool")
#  you shall find those who can give more zero closed priors to test
# install.packages("extraDistr")
## You have 5 priors 

## you can share them how you find the zero means by sampling and check 

library(extraDistr)
library("LaplacesDemon")

# the mean for half normal sqrt(2m/pi)= 0.7978846
# the bigger the m1, the smaller the mean to zero

#################################################
###Possible prior for the parameters#############
#################################################

sqrt(2/pi)
m1 = 10000
s1<-LaplacesDemon::rhalfnorm(1000,scale=sqrt(m1*pi/2)) 
mean(s1)

# the smaller the sigma, the smaller the mean
s2<-rhcauchy(1000, sigma = .0010) # + half cauchy
mean(s2)

## the smaller the scale the smaller the mean. scale seems not matter much
s3<-rhalft(1000, scale=0.01, nu=10)# + half t $\mathcal{Ht}$
mean(s3)

## mean exp(mu + sigma2/2)
s4<-rlnorm(1000, meanlog = -5, sdlog = 0.2)# + log normal
mean(s4)

## mean scale/(shape -1)
s5<-LaplacesDemon::rinvgamma(1000, shape=2, scale=0.01)
mean(s5)

## error likelihood use   tau ~ invgamma
s6<-rhcauchy(1, sigma = 0.1) # sd(y) + half cauchy
mean(s2)

s7<-rlnorm(1, meanlog = mean(y), sdlog = sd(y)) # + log normal
mean(s7)


#?rhalfnorm
alpha<-0.128
kappa<-0.236
n<-Ntip(tree)
N<-dim(tree$edge)[1]
betas.archrate<-array(NA,c(2*n-1))#this is Garch(1,0)
betas.garch11rate<-array(NA,c(2*n-1))#this is Garch(1,1)

plot(tree)
nodelabels()
tiplabels()

#bm.fit<-fitContinuous(tree,yoriginal,model="BM")
#omega<-sqrt(bm.fit$opt$sigsq)
##what RRphylo did, use their initial 
##  it could be hard to make the same array though but shall be a method to approach it.

h <- mle(optL, start = list(lambda = 1), method = "L-BFGS-B",upper = 10, lower = 0.001)
lambda <- h@coef
lambda
betas.rrphylo <- (solve(t(L) %*% L + lambda * diag(ncol(L))) %*%t(L)) %*% (as.matrix(y) - rootV)
betas.rrphylo
y.hat.rrphylo <- (L %*% betas.rrphylo) + rootV
ace.rrphylo <- (L1 %*% betas.rrphylo[1:Nnode(t) ]) + rootV

names(betas.archrate)<-colnames(L)
betas.archrate<-betas.archrate[c(n:(2*n-1),1:(n-1))]
betas.archrate

names(betas.garch11rate)<-colnames(L)
betas.garch11rate<-betas.garch11rate[c(n:(2*n-1),1:(n-1))]
betas.garch11rate

betas.rrphylo<-betas.rrphylo[c(n:(2*n-1),1:(n-1)),]
betas.rrphylo

betas.archrate[names(betas.archrate)==n+1] <- betas.rrphylo[n+1]#sqrt(bm.fit$opt$sigsq)
betas.archrate

betas.garch11rate[names(betas.garch11rate)==n+1] <- betas.rrphylo[n+1]#sqrt(bm.fit$opt$sigsq)
betas.garch11rate

anc<-tree$edge[,1]
des<-tree$edge[,2]
treelength<-tree$edge.length
cbind(tree$edge,tree$edge.length)

for(index in N:1){
  #  index<-N
  #        betas.archrate[des[index]]<- betas.rrphylo[rownames(betas.rrphylo)==des[index]] + alpha/treelength[index]*betas.archrate[anc[index]]#this would not go to
  betas.archrate[des[index]]<- betas.rrphylo[des[index]] + alpha/treelength[index]*betas.archrate[anc[index]]
  
  ## Here this is garch(1,1)
  betas.garch11rate[des[index]]<- betas.rrphylo[des[index]]+ kappa/treelength[index]*betas.garch11rate[anc[index]]*rnorm(1,0,treelength[index])^2  + alpha/treelength[index]*betas.archrate[anc[index]]
  
  #c(betas.archrate[anc[index]], betas.rrphylo[anc[index],])
  # c(betas.archrate[des[index]], betas.rrphylo[des[index],])
  # #c(betas.archrate[names(betas.archrate)==anc[index]], betas.rrphylo[rownames(betas.rrphylo)==anc[index],])
  # c(betas.archrate[names(betas.archrate)==des[index]], betas.rrphylo[rownames(betas.rrphylo)==des[index],])
}
sim.df<-cbind(betas.rrphylo,betas.archrate,betas.garch11rate)
colnames(sim.df)<-c("ridge","garch10","garch11")
head(sim.df)

betas.rrphylo<-betas.rrphylo[c((n+1):(2*n-1),1:n)]
betas.archrate<-betas.archrate[c((n+1):(2*n-1),1:n)]
betas.garch11rate<-betas.garch11rate[c((n+1):(2*n-1),1:n)]
cbind(betas.archrate,betas.rrphylo)

y.hat.rrphylo <- (L %*% betas.rrphylo) + rootV
ace.rrphylo <- (L1 %*% betas.rrphylo[1:Nnode(t) ]) + rootV

y.hat.archrate <- (L %*% betas.archrate) + rootV
ace.archrate <- (L1 %*% betas.archrate[1:Nnode(t) ]) + rootV

y.hat.garch11rate <- (L %*% betas.garch11rate) + rootV
ace.garch11rate <- (L1 %*% betas.garch11rate[1:Nnode(t) ]) + rootV

### Here is prediction to the trait
par(mfrow=c(1,3))
plot(yoriginal,y.hat.rrphylo,xlim=range(yoriginal,y.hat.rrphylo),ylim=range(yoriginal,y.hat.rrphylo),main= paste("RRphylo, lambda=", round(lambda,3), sep=""))
abline(a=0,b=1)

plot(yoriginal,y.hat.archrate,xlim=range(yoriginal,y.hat.archrate),ylim=range(yoriginal,y.hat.archrate),main=paste("Garch(1,0) Rate, alpha =", alpha ,sep=""))
abline(a=0,b=1)


plot(yoriginal,y.hat.garch11rate,xlim=range(yoriginal,y.hat.garch11rate),ylim=range(yoriginal,y.hat.garch11rate),main=paste("Garch(1,1) Rate, alpha =", alpha, " kappa =",kappa,sep=""))
abline(a=0,b=1)

### Here is ancestral estimation

par(mfrow=c(1,3))
makeL1(tree)->L1
ace <- L1 %*% betas.rrphylo[1:Nnode(tree)]
# plot(c(ace.rrphylo,yoriginal),betas.rrphylo,bg=c(rep("red",ntaxa-1),rep("green",ntaxa)),pch=21,cex=1.5,mgp=c(1.8,0.5,0),xlab="phenotypes at nodes and tips",ylab="rates",main="rates vs simulated phenotypes")
# legend("topright",legend=c("nodes","tips"),fill=c("red","green"),bty="n")
# compute ancestral character estimates at nodes as within RRphylo

plot(ace.rrphylo,ace,pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),xlab="simulated",ylab="estimated",main="RRphylo: phenotype at nodes (aces)",xlim=range(c(ace.rrphylo,ace)),ylim=range(c(ace.rrphylo,ace)) )
abline(a=0,b=1)

plot(ace.archrate,ace,pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),xlab="simulated",ylab="estimated",main="Garch(1,0): phenotype at nodes (aces)",xlim=range(c(ace.archrate,ace)),ylim=range(c(ace.archrate,ace)) )
abline(a=0,b=1)

plot(ace.garch11rate,ace,pch=21,bg="gray80",cex=1.5,mgp=c(1.8,0.5,0),xlab="simulated",ylab="estimated",main="Garch(1,1): phenotype at nodes (aces)",xlim=range(c(ace.garch11rate,ace)),ylim=range(c(ace.garch11rate,ace)) )
abline(a=0,b=1)

# Linear Regression using bayesian statistics Metropolis-Hastings MCMC in R
# https://khayatrayen.github.io/MCMC.html

