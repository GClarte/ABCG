library(ggplot2)
library(ggpubr)
library(smfsb)

#generation of the dataset
a=2
var=1
sigm=1
Nk=10
thetreel=c(2.774149)
obs=c(2.5666002)
# thetreel=rnorm(20,a,var)
# obs=numeric(20)
# for (i in 1:20){
#   obs[i]=mean(rnorm(10,thetreel[i],sigm))
# }

#parameter and hyperparameter update

#data = dataset
#par, hyper values at the previous step
#sigm = variance of the observations
#var = variance of the parameter
#qq=nombre of observations in each class
gibbsparam <- function(data, hyper, var, sigm, nbeps, qq) {
  p = length(data)
  thetc = numeric(p)
  for (i in 1:p) { 
    thettest = rnorm(nbeps, hyper, var)
    test = rowMeans(matrix(rnorm(qq * nbeps, thettest, sigm), nrow = nbeps))
    dist = abs(test - data[i])
    thetc[i] = thettest[which.min(dist)]
  }
  return(thetc)
}


#uniform prior on [0,400] for the parameters

gibbshyper <- function(thet, nbeps2, var) {
  res = runif(nbeps2, -4, 4)
  test = rowMeans(matrix(rnorm(length(thet) * nbeps2, res, var), ncol=length(thet)))
  dist = abs(test - mean(thet))
  return(res[which.min(dist)])
}

gibbstot <- function(data, thetini, hyperini, sigm, var, nbeps1, nbeps2, nbpts, qq) {
  reshyper = rep(NA, nbpts + 1)
  resparam = matrix(NA, ncol = nbpts + 1, nrow = length(thetini))
  reshyper[1] = hyperini
  resparam[,1] = thetini
  for (i in 2:(nbpts + 1)) {
    resparam[,i] = gibbsparam(data, reshyper[i - 1], var, sigm, nbeps1, qq)
    reshyper[i] = gibbshyper(resparam[, i], nbeps2, var)
  }
  return(list(reshyper,resparam))
}


#for comparison, simple ABC
#here, ntot=nbpts*(qq*nbeps1 + nbeps2)
ABCsimple <- function(data, var, sigm, nbpts, ntot, qq) {  
  reshyper = rep(NA, ntot)
  resparam = matrix(NA, ncol = ntot, nrow = length(data))
  dist = rep(NA, ntot)
  for (i in 1:ntot) {
    reshyper[i] = runif(1, -4, 4)
    resparam[,i] = rnorm(length(data), reshyper[i], var)
    dist[i] = sum(abs(data - rowMeans(
      matrix(rnorm(length(data) * qq, rep(resparam[,i], qq), sigm), nrow = length(data))
    )))
  }
  # compute the distance in the end
  v = order(dist)[1:nbpts]
  return(list(reshyper[v], resparam[,v]))
}


#functions used in ABC SMC


rdist <- function(par){
  sum(abs(Statstar - rowMeans(matrix(rnorm( 10, rep(par[2:length(par)], 10), sigm), nrow = 1))))
}

rprior <- function(){
  P=numeric(2)
  P[1]=runif(1,-4,4)
  P[2]=rnorm(1,P[1],var)
  return(P)
}

dprior <- function(P,...){
  dunif(P[1],-4,4,...)+sum(dnorm(P[2],P[1],var,...))
}

rperturb <- function(p){
  p+rnorm(2,0,.1)
}

dperturb <- function(p,pancien,...){
  sum(dnorm(p-pancien,0,.1,...))
}

Statstar=obs

#Gibbs exact changer
gibbsexact <- function(data,k,sigm,var,qq)
{
  R1=matrix(NA,ncol=20,nrow=k)
  R2=rep(NA,k)
  R2[1]=rnorm(1,-4,4)
  R1[1,]=rnorm(20,(qq*data/var + R2[1]/sigm)/(1/sigm + qq/var),1/(1/sigm + qq/var))
  for (i in 2:k)
  {
    cand = runif(1,-4,4)
    u=runif(1,0,1)
    if(u<(prod(dnorm(R1[i-1,],cand,var)))/(prod(dnorm(R1[i-1,],R2[i-1],var)))){
      R2[i]=cand
    } else {
      R2[i]=R2[i-1]
    }
    R1[i,]=rnorm(20,(qq*data/var + R2[i]/sigm)/(1/sigm + qq/var),1/(1/sigm + qq/var))
  }
  return(list(R1,R2))
}

#plusieurs densités superposées
N <- 100

hyper = numeric()
par1 = numeric()
nbr=numeric()
type=numeric()
exact=gibbsexact(obs,1000,sigm,var,10)



N=100
for (i in 1:N){
  A=gibbstot(obs,rep(0,20),0,sigm,var,30,30,floor(1000/30),10)
  B=ABCsimple(obs,var,sigm,floor(1000/30),1000,10)
  C=abcSmc(1000, rprior, dprior, rdist, rperturb,dperturb, verb=TRUE, steps=5, factor=6)
  hyper=c(hyper,A[[1]][5:34],B[[1]],C[,1])
  par1=c(par1,A[[2]][1,5:34],B[[2]],C[,2])
  nbr=c(nbr,rep(i,length(A[[1]][5:34])+length(B[[1]])+length(C[,1])))
  type=c(type,c(rep("ABC-Gibbs",length(A[[1]][5:34])),rep("Simple ABC",length(B[[1]])),rep("SMC-ABC",length(C[,1]))))
}


gib_curves30 <- data.frame(Group = nbr,
                      mu1 = par1,
                      hyperparameter = hyper,
                      Method=type)



g=ggplot(data=gib_curves30,aes(x=mu1))+theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
for (i in 1:N){
  g=g+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="ABC-Gibbs"),],stat = "density",alpha=.1)
  g=g+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="Simple ABC"),],stat = "density",alpha=.1)
  g=g+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="SMC-ABC"),],stat = "density",alpha=.1)
}
g=g+facet_grid(Method~.,scales = "free" )+xlab(expression(mu[1]))
g=g+geom_line(data=data.frame(mu1=exact[[1]][,1],Type=rep("Exact",1000)),stat = "density",linetype="dashed")+ scale_colour_viridis_d()
g
ggsave("mu1comp1D.pdf",height=15,width=7,unit="cm")


f=ggplot(data=gib_curves30,aes(x=hyperparameter))+theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
for (i in 1:N){
  f=f+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="ABC-Gibbs"),],stat = "density",alpha=0.2)
  f=f+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="Simple ABC"),],stat = "density",alpha=0.2)
  f=f+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="SMC-ABC"),],stat = "density",alpha=0.2)
}
f=f+facet_grid(Method~.,scales = "free" )+xlab(expression(alpha))
f=f+geom_line(data=data.frame(hyperparameter=exact[[2]],Type=rep("Exact",1000)),stat = "density",linetype="dashed")+xlim(c(0,5))+ scale_colour_viridis_d()
f
ggsave("hypercomp1D.pdf",height=15,width=7,unit="cm")


