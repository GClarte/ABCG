

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
  dist=numeric(k)
  for (i in 2:k){
    for (j in 1:n){
      As[i,j]=ABCAbas(M[1],X[,j],hyper[i-1],B[i-1],g[i-1],K[i-1])
    }
    hyper[i]=ABCAhaut(M[2],As[i,])
    B[i]=ABCB(M[3],X,As[i,],g[i-1],K[i-1])
    g[i]=ABCg(M[4],X,As[i,],B[i],K[i-1])
    K[i]=ABCk(M[5],X,As[i,],B[i],g[i])
  }
  return(list(cbind(hyper,As,B,g,K),))
}


abcvan <- function(k,N,Statstar){
  Par=matrix(NA,ncol=54,nrow=N)
  dist=numeric(N)
  for (i in 1:N){
    Par[i,]=rprior()
    dist[i]=rdist(Par[i,],Statstar)
  }
  V=order(dist)[1:k]
  return(list(Par[V,],max(dist(V))))
}



SMCmaison <- function(npart,M,itermax,epstarget,alph,data,prior,model){
  par=sapply(1:npart,function(x){prior$simu()})
  stats=sapply(1:npart,function(x){sapply(1:M,function(y){model(par[,x],data)})})
  j=1
  pds=rep(1/npart,npart)
  ESS=npart
  eps=max(stats)
  histeps=eps
  while(j< itermax && eps>epstarget){
    VV=changementeps(stats,eps,pds,alph)
    eps=VV[[1]]
    pds=VV[[2]]
    ESS=VV[[3]]
    histeps=c(histeps,eps)
    if (ESS<npart/2){
      qui=sample(1:npart,prob=pds,rep=TRUE)
      par=par[,qui]
      pds=rep(1/npart,npart)
      ESS=npart
      stats=stats[,qui]
    }
    VV=pasnoyau(par,eps,stats,pds,data,prior,model)
    par=VV[[1]]
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
  return(list(par,pds,stats,histeps))
}



#par : matrice des paramètres des particules actuelles, en colonne les coordonnées des paramètres
#simu : matrice des distances entre obs et simus
pasnoyau <- function(par,eps,stats,pds,data,prior,model){
  sd=2*cov(t(par))
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
