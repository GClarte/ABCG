library(gk)
library(smfsb)
library(ggplot2)

ABCAbas <- function(N,x,hyper,B,g,k){
  par=rnorm(N,hyper,1)
  dist=sapply(par,function(y){Z=rgk(length(x),y,B,g,k);return(sum(abs(quantile(Z,(0:8)/8)-quantile(x,(0:8)/8))))})
  return(list(par[which.min(dist)],min(dist)))
}


ABCAhaut <- function(N,As){
  par=runif(N,-10,10)
  l=length(As)
  dist=sapply(par,function(y){abs(mean(rnorm(l,y,1))-mean(As))})
  return(par[which.min(dist)])
}


gibbs <- function(M,x,K){
  n=ncol(X)
  l=nrow(X)
  As=matrix(NA,ncol=n,nrow=K)
  hyper=numeric(K)
  hyper[1]=runif(1,-10,10)
  As[1,]=rnorm(n,hyper[1],1)
  distacc=matrix(ncol=K-1,nrow=n)
  for (i in 2:K){
    for (j in 1:n){
      VV=ABCAbas(M[1],X[,j],hyper[i-1],B,g,k)
      As[i,j]=VV[[1]]
      distacc[j,i-1]=VV[[2]]
    }
    hyper[i]=ABCAhaut(M[2],As[i,])
  }
  return(list(cbind(hyper,As),distacc))
}

#nouveau ABCSMC

myabcSmc <- function(N, rprior, dprior, rdist, rperturb, dperturb, factor=10,
                   steps=15, verb=FALSE) {
  priorLW = log(rep(1/N, N))
  priorSample = mclapply(as.list(priorLW), function(x) {rprior()})
  quant=numeric(steps)
  for (i in steps:1) {
    if (verb) message(paste(i,""), appendLF=FALSE)
    out = myabcSmcStep(dprior, priorSample, priorLW, rdist, rperturb,
                     dperturb, factor)
    priorSample = out[[1]]
    priorLW = out[[2]]
    quant[i]=out[[3]]
  }
  if (verb) message("")
  return(list(t(sapply(sample(priorSample, N, replace=TRUE, prob=exp(priorLW)), identity)),quant))
}

myabcSmcStep <- function(dprior, priorSample, priorLW, rdist, rperturb,
                       dperturb, factor=10) {
  n = length(priorSample)
  mx = max(priorLW)
  rw = exp(priorLW - mx)
  prior = sample(priorSample, n*factor, replace=TRUE, prob=rw)
  prop = mcMap(rperturb, prior)
  dist = mcMap(rdist, prop) # forward simulation in parallel
  qCut = quantile(unlist(dist), 1/factor)
  new = prop[dist < qCut]
  lw = mcMap( function(th) {
    terms = priorLW + sapply(priorSample,
                             function(x){dperturb(th, x, log=TRUE)})
    mt = max(terms)
    denom = mt + log(sum(exp(terms - mt)))
    dprior(th, log=TRUE) - denom
  } , new)
  lw = unlist(lw)
  mx = max(lw)
  lw = lw - mx
  nlw = log(exp(lw)/sum(exp(lw)))
  list(sample = new, lw = nlw,qCut)
}


#versions avec smfsb
#paramètres de la forme par=c(hyper,As,B,g,k)

rdist <- function(par){
  stat=numeric()
  for (i in 1:50){
    stat=c(stat,quantile(rgk(20,par[1+i],B,g,k),(0:8)/8))
  }
  return(sum(abs(stat-Statstar)))
}

rprior <- function(){
  P=numeric(51)
  P[1]=runif(1,-10,10)
  P[2:51]=rnorm(50,P[1],1)
  return(P)
}

dprior <- function(P,...){
  dunif(P[1],-10,10,...)+sum(dnorm(P[2:51],P[1],1,...))
}

rperturb <- function(p){
  p+rnorm(51,0,.1)
}

dperturb <- function(p,pancien,...){
  sum(dnorm(p-pancien,0,.1,...))
}


#version vanille
abcvan <- function(k,N,x){
  Par=matrix(NA,ncol=51,nrow=N)
  dist=numeric(N)
  for (i in 1:N){
    Par[i,]=rprior()
    dist[i]=rdist(Par[i,])
  }
  V=order(dist)[1:k]
  return(list(Par[V,],sort(dist)[k]))
}

#simulations

hyper=runif(1,-10,10)
As=rnorm(50,hyper,1)
B=runif(1,0,1)
g=runif(1,0,1)
k=runif(1,0,1)

X=matrix(NA,ncol=50,nrow=20)
for (i in 1:50){
  X[,i]=rgk(20,As[i],B,g,k)
}



Statstar=numeric()
for (i in 1:50){
  Statstar=c(Statstar,quantile(X[,i],(0:8)/8))
}



#comparaison

outSMC = myabcSmc(1000, rprior, dprior, rdist, rperturb,dperturb, verb=TRUE, steps=10, factor=10)
outVan = abcvan(1000,100*1000,X)
outgib = gibbs(c(100,50),X,1000)

#comparaison des distances ABC

mean(apply(outVan[[1]],1,rdist))
sd(apply(outVan[[1]],1,rdist))

mean(apply(outSMC[[1]],1,rdist))
sd(apply(outSMC[[1]],1,rdist))

mean(apply(outgib[[1]],1,rdist))
sd(apply(outgib[[1]][10:1000,],1,rdist))

mean(apply(matrix(c(0,As),nrow=100,ncol=51,byrow=T),1,rdist))
sd(apply(matrix(c(0,As),nrow=100,ncol=51,byrow=T),1,rdist))

Dat <- data.frame(
  value = c(outSMC[,1], outVan[,1], outgib[,1], 
            outSMC[,2], outVan[,2], outgib[,2],
            outSMC[,3], outVan[,3], outgib[,3],
            outSMC[,4], outVan[,4], outgib[,4], 
            outSMC[,5], outVan[,5], outgib[,5]),
  Method = rep(c(rep("ABC-SMC",1000), rep("vanilla ABC",1000), rep("ABC Gibbs",1000)),5),
  Parameter = c(rep("Hyperparameter",3000),
                rep("mu1",3000), rep("mu2",3000), rep("mu3",3000),rep("mu4",3000)))

Dat$Parameter <- factor(Dat$Parameter, 
                        levels = c("Hyperparameter", "mu1", "mu2", "mu3", "mu4"), 
                        labels = c('Hyperparameter' = expression(alpha),
                                   'mu1' = expression(mu[1]),
                                   'mu2' = expression(mu[2]),
                                   'mu3' = expression(mu[3]),
                                   'mu4' = expression(mu[4]))
)

theta_star <- data.frame(value = c(hyper,As[1:4]),
                         Parameter = levels(Dat$Parameter))


f <- ggplot(data = Dat) + stat_density(aes(x = value, linetype = Method), 
                                       geom="line", position="identity") + 
  geom_vline(data = theta_star, aes(xintercept = value)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank()) + coord_cartesian(xlim = c(-10,-6), ylim = c(0,5), expand = TRUE,
                                                        default = FALSE, clip = "on")
# geom_density(data=Dat,aes(x=value,color=Method)) 

f <- f + facet_wrap(Parameter~., labeller = label_parsed, ncol = 3) + scale_colour_viridis_d()
f

ggsave("gksimple_nouv.pdf",height=10,width=15,units = 'cm')


Dat=data.frame(value=c(outSMC[,1],outVan[,1],outgib[,1],outSMC[,2],outVan[,2],outgib[,2],outSMC[,3],outVan[,3],outgib[,3],outSMC[,4],outVan[,4],outgib[,4],outSMC[,5],outVan[,5],outgib[,5]),Method=rep(c(rep("ABC-SMC",1000),rep("vanilla ABC",1000),rep("ABC Gibbs",1000)),5),Parameter=c(rep("Hyperparameter",3000),rep("1",3000),rep("2",3000),rep("3",3000),rep("4",3000)))
f=ggplot(data = Dat) + geom_density(data=Dat,aes(x=value,color=Method)) +geom_vline(data=data.frame(Truth=c(hyper,As[1:4]),Parameter=c("Hyperparameter","1","2","3","4")),mapping = aes(xintercept = Truth))
f=f+facet_grid(Parameter~.)+theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")+scale_colour_viridis_d() + xlim(c(-2.5,3))
f
ggsave("gksimple191012.pdf",height=7,width=25,units = 'cm')


