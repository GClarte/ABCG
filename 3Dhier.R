


library(ggplot2)
library(ggpubr)
library(smfsb)
library(EasyABC)
library(mvtnorm)


gibbsparam <- function(data, hyper, var, sigm, nbeps, qq) {
  #gibbs step for the parmeter, prior induced by the hyperparameter
  p = length(data)
  thetc = numeric(p)
  dists=0
  for (i in 1:p) { 
    thettest = rnorm(nbeps, hyper, var)
    test = rowMeans(matrix(rnorm(qq * nbeps, thettest, sigm), nrow = nbeps))
    dist = abs(test - data[i])
    thetc[i] = thettest[which.min(dist)]
    dists=dists+min(dist)
  }
  return(list(thetc,dists))
}



gibbshyper <- function(thet, nbeps2, var) {
  #gibbs step for the hyperparameter, with uniform prior
  res = runif(nbeps2, -4, 4)
  test = rowMeans(matrix(rnorm(length(thet) * nbeps2, res, var), ncol=length(thet)))
  dist = abs(test - mean(thet))
  return(list(res[which.min(dist)],min(dist)))
}

gibbstot <- function(data, thetini, hyperini, sigm, var, nbeps1, nbeps2, nbpts, qq) {
  #full function
  reshyper = rep(NA, nbpts + 1)
  resparam = matrix(NA, ncol = nbpts + 1, nrow = length(thetini))
  reshyper[1] = hyperini
  resparam[,1] = thetini
  resdist = rep(NA,nbpts)
  for (i in 2:(nbpts + 1)) {
    resdist[i-1]= 0
    VV=gibbsparam(data, reshyper[i - 1], var, sigm, nbeps1, qq)
    resparam[,i] = VV[[1]]
    resdist[i-1]= resdist[i-1]+VV[[2]]
    WW=gibbshyper(resparam[, i], nbeps2, var)
    reshyper[i] = WW[[1]]
  }
  return(list(reshyper,resparam,resdist))
}


#for comparison, simple ABC
#here, ntot=nbpts*(qq*nbeps1 + nbeps2)
ABCsimple <- function(data, vari, sigma, nbpts, ntot, qq) {  
  reshyper = rep(NA, ntot)
  resparam = matrix(NA, ncol = ntot, nrow = length(data))
  dist = matrix(NA, nrow=length(data),ncol=ntot)
  for (i in 1:ntot) {
    reshyper[i] = runif(1, -4, 4)
    resparam[,i] = rnorm(length(data), reshyper[i], vari)
    dist[i] = sum(abs(data - rowMeans(
      matrix(rnorm(length(data) * qq, rep(resparam[,i], qq), sigma), nrow = length(data))
    )))
  }
  # compute the distance in the end
  v = order(dist)[1:nbpts]
  return(list(reshyper[v], resparam[,v],dist[v[nbpts]]))
}

SMCmaison <- function(npart,M,itermax,epstarget,alph,data,prior,model){
  par=sapply(1:npart,function(x){prior$simu()})
  stats=sapply(1:npart,function(x){sapply(1:M,function(y){model(par[,x],data)})})
  j=1
  pds=rep(1/npart,npart)
  ESS=npart
  eps=max(stats)
  histeps=eps
  histpart=matrix(ncol=npart,nrow=2*itermax)
  histpds=matrix(ncol=npart,nrow=itermax)
  histpds[1,]=pds
  while(j< itermax && eps>epstarget){
    VV=changementeps(stats,eps,pds,alph)
    eps=VV[[1]]
    pds=VV[[2]]
    ESS=VV[[3]]
    histeps=c(histeps,eps)
    histpart[2*(j-1)+1,]=par[1,]
    if (ESS<npart/2){
      kk=1
      sd=2*cov(t(par[,which(pds!=0)]))
      if(any(diag(sd)==0)){
        sd=sd+diag(.1,nrow(par))
      }
      part=par
      statst=stats
      for (i in 1:npart){
        NON=TRUE
        while(NON & kk<1000000){
          qui=sample(1:npart,1,prob=pds,rep=TRUE)
          para=par[,qui]
          statsa=stats[,qui]
          part=para+rmvnorm(1,sigma=sd)
          statst=sapply(1:length(statsa),function(x){model(part,data)})
          u=runif(1)
          kk=kk+1
          if( u < (sum(statst<=eps)/sum(statsa<=eps)*prior$density(part)/prior$density(para))){
            par[,i]=part
            stats[,i]=statst
            NON=FALSE
          }
        }
        if (kk==1000000){
          print("n'y arrive pas")
          return(list(par,pds,stats,histeps,histpart))
        }
      }
      pds=rep(1/npart,npart)
      ESS=npart
      histpart[2*j,]=par[1,]
    }
    VV=pasnoyau(par,eps,stats,pds,data,prior,model)
    par=VV[[1]]
    histpds[1,]=pds
    stats=VV[[2]]
    print(paste(j,eps))
    j=j+1
  }
  #on finit par un resampling pour avoir quelque chose de propre
  qui=sample(1:npart,prob=pds,rep=TRUE)
  par=par[,qui]
  pds=rep(1/npart,npart)
  ESS=npart
  stats=stats[,qui]
  return(list(par,pds,stats,histeps,histpart,histpds))
}


#par : matrice des paramètres des particules actuelles, en colonne les coordonnées des paramètres
#simu : matrice des distances entre obs et simus
pasnoyau <- function(par,eps,stats,pds,data,prior,model){
  sd=2*cov(t(par[,which(par[,i]!=0)]))
  if(any(diag(sd)==0)){
    sd=sd+diag(.1,nrow(par))
  }
  part=par
  statst=stats
  for (i in 1:ncol(par)){
    if(pds[i]>0){
      VV=chgtpar(par[,i],stats[,i],sd,eps,data,prior,model)
      part[,i]=VV[[1]]
      statst[,i]=VV[[2]]
    }
    
  }
  return(list(part,statst))
}

changementeps <- function(stats,eps,pds,alph){
  posseps = sort(stats[stats<=eps],decreasing = TRUE)
  pdst=pds
  epst=eps
  test=FALSE
  k=1
  newpds=pds
  ESS=1/(sum(pds^2))
  ESSt=ESS
  while (k<(length(posseps)-1) && !test){
    neweps=posseps[k]
    for (i in 1:length(pds)){
      if (newpds[i]!=0){
        newpds[i]=pds[i]*sum(stats[,i]<=neweps)/sum(stats[,i]<=eps)
      }
    }
    newpds=newpds/sum(newpds)
    newESS=1/sum(newpds^2)
    if (newESS<alph*ESS){
      test=TRUE
      pdst=newpds
      epst=neweps
      ESSt=newESS
    }
    k=k+1
  }
  return(list(epst,pdst,ESSt))
}


chgtpar <- function(par,stats,sd,eps,data,prior,model){
  part=par+rmvnorm(1,sigma=sd)
  statst=sapply(1:length(stats),function(x){model(part,data)})
  u=runif(1)
  if( u < (sum(statst<=eps)/sum(stats<=eps)*prior$density(part)/prior$density(par))){
    return(list(part,statst))
  } else {
    return(list(par,stats))
  }
}



gibbsexact <- function(data,k,sigm,var,qq)
{
  R1=matrix(NA,ncol=2,nrow=k)
  R2=rep(NA,k)
  R2[1]=rnorm(1,-4,4)
  R1[1,]=rnorm(2,(qq*data/var + R2[1]/sigm)/(1/sigm + qq/var),1/(1/sigm + qq/var))
  for (i in 2:k)
  {
    cand = runif(1,-4,4)
    u=runif(1,0,1)
    if(u<(prod(dnorm(R1[i-1,],cand,var)))/(prod(dnorm(R1[i-1,],R2[i-1],var)))){
      R2[i]=cand
    } else {
      R2[i]=R2[i-1]
    }
    R1[i,]=rnorm(2,(qq*data/var + R2[i]/sigm)/(1/sigm + qq/var),1/(1/sigm + qq/var))
  }
  return(list(R1,R2))
}

#generation of the dataset
a=2
#on code les variances en dur pour éviter les problèmes
vari=1
sigma=1
Nk=10
thetreel=c(2.774149,2.438900)
obs=c(2.5666002,2.5808230)

#pour SMCmaison

rpriortot <- function(){
  alpha=runif(1,-4,4)
  thet=rnorm(2,alpha)
  return(c(alpha,thet))
}

dpriortot <- function(x){
  if (abs(x[1])<4){
    return(1/4*prod(dnorm(x[2:3],x[1])))
  } else {
    return(0)
  }
}

monprior=list(density=dpriortot,simu=rpriortot)

model <- function(pars,data){
  return(sum(abs(rowMeans(matrix(rnorm(2 * 10, rep(pars[2:3], 10), 1), nrow = 2))-data)))
}

VV=SMCmaison(1000,10,100,0,.9,obs,monprior,model)
resS100=ABCsimple(obs,vari,sigma,1000,1000*1000,10)



hyper = numeric()
par1 = numeric()
par2 = numeric()
nbr=numeric()
type=numeric()
epsfin=matrix(ncol=3,nrow=10)
exact=gibbsexact(obs,1000,sigma,vari,10)


N=10
for (i in 1:N){
  A=gibbstot(obs,rep(0,2),0,sigma,vari,30,30,1000,10)
  B=ABCsimple(obs,vari,sigma,1000,30000,10)
  C=SMCmaison(1000,5,30,A[[3]][length(A[[3]])],.9,obs,monprior,model)
  hyper=c(hyper,A[[1]][5:50],B[[1]],C[[1]][1,])
  par1=c(par1,A[[2]][1,5:50],B[[2]][1,],C[[1]][2,])
  par2=c(par2,A[[2]][2,5:50],B[[2]][2,],C[[1]][3,])
  epsfin[1,]=c(A[[3]][length(A[[3]])],B[[3]],C[[3]][length(C[[3]])])
  nbr=c(nbr,rep(i,length(A[[1]][5:50])+length(B[[1]])+length(C[[1]][1,])))
  type=c(type,c(rep("ABC-Gibbs",length(A[[1]][5:50])),rep("Simple ABC",length(B[[1]])),rep("SMC-ABC",length(C[[1]][1,]))))
}


gib_curves30 <- data.frame(Group = nbr,
                           mu1 = par1,
                           mu2 = par2,
                           hyperparameter = hyper,
                           Method=type)



g=ggplot(data=gib_curves30,aes(x=mu1))+theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
for (i in 1:N){
  g=g+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="ABC-Gibbs"),],stat = "density",alpha=.2)
  g=g+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="Simple ABC"),],stat = "density",alpha=.2)
  g=g+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="SMC-ABC"),],stat = "density",alpha=.2)
}
g=g+facet_grid(Method~.,scales = "free_y" )+xlab(expression(mu[1]))
g=g+geom_line(data=data.frame(mu1=exact[[1]][,1],Type=rep("Exact",1000)),stat = "density",linetype="dashed")+ scale_colour_viridis_d()
g
ggsave("mu1SMC3d.pdf",height=15,width=7,unit="cm")

f=ggplot(data=gib_curves30,aes(x=hyperparameter))+theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
for (i in 1:N){
  f=f+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="ABC-Gibbs"),],stat = "density",alpha=0.2)
  f=f+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="Simple ABC"),],stat = "density",alpha=0.2)
  f=f+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="SMC-ABC"),],stat = "density",alpha=0.2)
}
f=f+facet_grid(Method~.,scales = "free" )+xlab(expression(alpha))
f=f+geom_line(data=data.frame(hyperparameter=exact[[2]],Type=rep("Exact",1000)),stat = "density",linetype="dashed")+xlim(c(0,5))+ scale_colour_viridis_d()
f
ggsave("hyperSMC3d.pdf",height=15,width=7,unit="cm")
