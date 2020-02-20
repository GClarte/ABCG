library(gk)
library(smfsb)
library(ggplot2)

ABCAbas <- function(N,x,hyper,B,g,k){
  par=rnorm(N,hyper,1)
  dist=sapply(par,function(y){Z=rgk(length(x),y,B,g,k);return(sum(abs(quantile(Z,(0:8)/8)-quantile(x,(0:8)/8))))})
  return(par[which.min(dist)])
}


ABCAhaut <- function(N,As){
  par=runif(N,-10,10)
  l=length(As)
  dist=sapply(par,function(y){abs(mean(rnorm(l,y,1))-mean(As))})
  return(par[which.min(dist)])
}


ABCB <- function(N,X,As,g,k){
  par=runif(N,0,1)
  n=ncol(X)
  l=nrow(X)
  dist=numeric(N)
  for (i in 1:N){
    dist[i]=sum(sapply(1:n,function(z){K=rgk(l,As[z],par[i],g,k);return(sum(abs(quantile(K,(0:8)/8)-quantile(X[,z],(0:8)/8))))}))
  }
  return(par[which.min(dist)])
}


ABCg <- function(N,X,As,B,k){
  par=runif(N,0,1)
  n=ncol(X)
  l=nrow(X)
  dist=numeric(N)
  for (i in 1:N){
    dist[i]=sum(sapply(1:n,function(z){K=rgk(l,As[z],B,par[i],k);return(sum(abs(quantile(K,(0:8)/8)-quantile(X[,z],(0:8)/8))))}))
  }
  return(par[which.min(dist)])
}


ABCk <- function(N,X,As,B,g){
  par=runif(N,0,1)
  n=ncol(X)
  l=nrow(X)
  dist=numeric(N)
  for (i in 1:N){
    dist[i]=sum(sapply(1:n,function(z){K=rgk(l,As[z],B,g,par[i]);return(sum(abs(quantile(K,(0:8)/8)-quantile(X[,z],(0:8)/8))))}))
  }
  return(par[which.min(dist)])
}

gibbs <- function(M,x,k){
  n=ncol(X)
  l=nrow(X)
  As=matrix(NA,ncol=n,nrow=k)
  B=numeric(k)
  g=numeric(k)
  K=numeric(k)
  hyper=numeric(k)
  hyper[1]=runif(1,-10,10)
  As[1,]=rnorm(n,hyper[1],1)
  B[1]=runif(1,0,1)
  g[1]=runif(1,0,1)
  K[1]=runif(1,0,1)
  for (i in 2:k){
    for (j in 1:n){
      As[i,j]=ABCAbas(M[1],X[,j],hyper[i-1],B[i-1],g[i-1],K[i-1])
    }
    hyper[i]=ABCAhaut(M[2],As[i,])
    B[i]=ABCB(M[3],X,As[i,],g[i-1],K[i-1])
    g[i]=ABCg(M[4],X,As[i,],B[i],K[i-1])
    K[i]=ABCk(M[5],X,As[i,],B[i],g[i])
  }
  return(cbind(hyper,As,B,g,K))
}


#versions avec smfsb
#paramètres de la forme par=c(hyper,As,B,g,k)

rdist <- function(par){
  stat=numeric()
  for (i in 1:50){
    stat=c(stat,quantile(rgk(20,par[1+i],par[52],par[53],par[54]),(0:8)/8))
  }
  return(sum(abs(stat-Statstar)))
}

rprior <- function(){
  P=numeric(54)
  P[1]=runif(1,-10,10)
  P[52:54]=runif(3,0,1)
  P[2:51]=rnorm(50,P[1],1)
  return(P)
}

dprior <- function(P,...){
  dunif(P[1],-10,10,...)+sum(dunif(P[52:54],0,1,...))+sum(dnorm(P[2:51],P[1],1,...))
}

rperturb <- function(p){
  p+rnorm(54,0,.1)
}

dperturb <- function(p,pancien,...){
  sum(dnorm(p-pancien,0,.1,...))
}


#version vanille
abcvan <- function(k,N,x){
  Par=matrix(NA,ncol=54,nrow=N)
  dist=numeric(N)
  for (i in 1:N){
    Par[i,]=rprior()
    dist[i]=rdist(Par[i,])
  }
  V=order(dist)[1:k]
  return(Par[V,])
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

outSMC = abcSmc(1000, rprior, dprior, rdist, rperturb,dperturb, verb=TRUE, steps=10, factor=10)
outVan = abcvan(1000,100*1000,X)
outgib = gibbs(c(100,50,50,50,50),X,1000)

#comparaison des distances ABC

mean(apply(outVan,1,rdist))
sd(apply(outVan,1,rdist))

mean(apply(outSMC,1,rdist))
sd(apply(outSMC,1,rdist))

mean(apply(outgib[10:1000,],1,rdist))
sd(apply(outgib[10:1000,],1,rdist))

#plots pour le paramètres
Dat=data.frame(parameter=c(outSMC[,2],outVan[,2],outgib[,2]),type=c(rep("ABC-SMC",1000),rep("vanilla ABC",1000),rep("ABC Gibbs",1000)))
f=ggplot(data = Dat) + geom_density(data=Dat,aes(x=parameter,color=type)) + geom_vline(xintercept = As[1])
f

colMeans(outSMC)-c(hyper,As,B,g,k)
colMeans(outVan)-c(hyper,As,B,g,k)
colMeans(outgib)-c(hyper,As,B,g,k)

apply(outSMC,2,sd)
apply(outVan,2,sd)
apply(outgib,2,sd)


Dat <- data.frame(
  value = c(outSMC[,2], outVan[,2], outgib[,2], 
            outSMC[,3], outVan[,3], outgib[,3], 
            outSMC[,4], outVan[,4], outgib[,4],
            outSMC[,5], outVan[,5], outgib[,5]),
  Method = rep(c(rep("ABC-SMC", 1000), rep("vanilla ABC",1000),
                 rep("ABC Gibbs",1000)), 4), 
  Parameter = c(rep("mu1", 3000), rep("mu2", 3000), rep("mu3", 3000), rep("mu4", 3000))
  )

Dat$Parameter <- factor(Dat$Parameter, 
                             levels = c("mu1", "mu2", "mu3", "mu4"), 
                             labels = c('mu1' = expression(mu[1]),
                                        'mu2' = expression(mu[2]),
                                        'mu3' = expression(mu[3]),
                                        'mu4' = expression(mu[4]))
)

theta_star <- data.frame(value = As[1:4],
                         Parameter = levels(Dat$Parameter))

f <- ggplot(data = Dat) + stat_density(aes(x = value, linetype = Method), 
                                       geom="line", position="identity") + 
  geom_vline(data = theta_star, aes(xintercept = value)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
f <- f + coord_cartesian(xlim = c(-4,-0), ylim = c(0,5), expand = TRUE,default = FALSE, clip = "on")
  # geom_density(data=Dat,aes(x=value,color=Method)) 
    
f <- f + facet_grid(.~Parameter, labeller = label_parsed) + scale_colour_viridis_d()
f
ggsave("4parbis_nouv.pdf",height=7,width=15,units = 'cm')

#plots pour les hyperparamètres

Dat <- data.frame(value = c(outSMC[,1], outVan[,1], outgib[,1],
                            outSMC[,52], outVan[,52], outgib[,52],
                            outSMC[,53], outVan[,53], outgib[,53],
                            outSMC[,54], outVan[,54], outgib[,54]),
               Method = rep(c(rep("ABC-SMC",1000), rep("vanilla ABC",1000), 
                              rep("ABC Gibbs",1000)),4),
               Parameter = c(rep("hyperparameter",3000), rep("B",3000), rep("g",3000), rep('k',3000)))

Dat$Parameter <- factor(Dat$Parameter, 
                        levels = c("hyperparameter", "B", "g", "k"), 
                        labels = c('hyperparameter' = expression(alpha),
                                   'B' = expression(B),
                                   'g' = expression(g),
                                   'k' = expression(k))
)

theta_star <- data.frame(value = c(hyper,B,g,k),
                         Parameter = levels(Dat$Parameter))

f <- ggplot(data = Dat) + stat_density(aes(x = value, linetype = Method),
                                       geom="line", position="identity") +
  geom_vline(data = theta_star, aes(xintercept = value)) +
  theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())
# geom_density(data=Dat,aes(x=value,color=Method)) 
f <- f + coord_cartesian(ylim = c(0,10), expand = TRUE,default = FALSE, clip = "on")
f <- f + facet_grid(.~Parameter, labeller = label_parsed,scale="free") + scale_colour_viridis_d()
f

ggsave("autresparbis_nouv.pdf",height=7,width=15,units = 'cm')

Dat=data.frame(parameter=c(outSMC[,1],outVan[,1],outgib[,1]),type=c(rep("ABC-SMC",1000),rep("vanilla ABC",1000),rep("ABC Gibbs",1000)))
f=ggplot(data = Dat) + geom_density(data=Dat,aes(x=parameter,color=type))
f

Dat=data.frame(parameter=c(outSMC[,52],outVan[,52],outgib[,52]),type=c(rep("ABC-SMC",1000),rep("vanilla ABC",1000),rep("ABC Gibbs",1000)))
f=ggplot(data = Dat) + geom_density(data=Dat,aes(x=parameter,color=type))
f

Dat=data.frame(parameter=c(outSMC[,53],outVan[,53],outgib[,53]),type=c(rep("ABC-SMC",1000),rep("vanilla ABC",1000),rep("ABC Gibbs",1000)))
f=ggplot(data = Dat) + geom_density(data=Dat,aes(x=parameter,color=type))
f
Dat=data.frame(parameter=c(outSMC[,54],outVan[,54],outgib[,54]),type=c(rep("ABC-SMC",1000),rep("vanilla ABC",1000),rep("ABC Gibbs",1000)))
f=ggplot(data = Dat) + geom_density(data=Dat,aes(x=parameter,color=type))
f
