library(gk)
library(smfsb)
library(ggplot2)
library(mvtnorm)


source("functions_gk_simple.R")




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

#définition du modèle pour SMCmaison
rdist <- function(par,Statstar){
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


model=rdist

monprior =list(density=dprior,simu=rprior)



outVan = abcvan(1000,100*1000,Statstar)
outgib = gibbs(c(100,50,50,50,50),X,1000)
outSMC = SMCmaison(1000,5,500,0,.9,Statstar,monprior,model)



Dat <- data.frame(
  value = c(outSMC[[1]][1,], outVan[,1], outgib[[1]][,1], 
            outSMC[[1]][2,], outVan[,2], outgib[[1]][,2],
            outSMC[[1]][3,], outVan[,3], outgib[[1]][,3],
            outSMC[[1]][4,], outVan[,4], outgib[[1]][,4], 
            outSMC[[1]][5,], outVan[,5], outgib[[1]][,5]),
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
        axis.title.y=element_blank()) + coord_cartesian(xlim = c(8,12), ylim = c(0,10), expand = TRUE,
                                                        default = FALSE, clip = "on")
# geom_density(data=Dat,aes(x=value,color=Method)) 

f <- f + facet_wrap(Parameter~., labeller = label_parsed, ncol = 3) + scale_colour_viridis_d()
f

ggsave("gksimpleSMC_nouv.pdf",height=10,width=15,units = 'cm')
  