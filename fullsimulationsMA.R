#samples n times a set of parameters for MA model for a given value of the hyperparameter
prior <- function(n,hyper){
  X_1=rbeta(n,hyper[1],hyper[2]+hyper[3])
  Temp=rbeta(n,hyper[2],hyper[3])
  X_2=Temp*(1-X_1)
  R_1=X_1-X_2
  R_2=2*(X_1+X_2)-1
  return(cbind(R_1,R_2))
}

#simulates a trajectory of the MA model with tht1 and tht2 the paramters and sigm the variance
simu <- function(tht1,tht2,sigm,n)
{
  A=rnorm(n+2,0,sigm)
  U=A[3:(n+2)]
  V=A[2:(n+1)]
  W=A[1:n]
  return(U+tht1*V+tht2*W)
}

#computes autocorrelation distance between two trajectories
distautcor <- function(Z,X)
{
  p=length(Z)
  a2=sum(Z[3:p]*Z[1:(p-2)])
  b2=sum(X[3:p]*X[1:(p-2)])
  a1=sum(Z[2:p]*Z[1:(p-1)])
  b1=sum(X[2:p]*X[1:(p-1)])
  return(sqrt((a2-b2)^2+(a1-b1)^2))
}

#ABC step to infer the parameter of a single trajectory X for a given hyperparameter hyper
#and a given variance sigm. N is the number of simulations.
ABCbas <- function(X,N,sigm,hyper)
{
  res=prior(N,hyper)
  dist=apply(res,1,function(x){distautcor(X-mean(X),simu(x[1],x[2],sigm,length(X)))})
  return(res[which.min(dist),])
}

#computes the sufficient statistic for Dirichlet, U is a matrix with two colums, for each
#parameter, and each row corresponds to one set of parameters
statdir <- function(U){
  q=(U[,2]+2*U[,1]+1)/4
  return(c(sum(log(q)),sum(log(abs(q-U[,1])))))
}

#ABC step for the inference of hyperparameters, Q1 and Q2 are the vectors of the first
#and second parameter respectively. L is the number of simulations.
ABChaut <- function(Q1,Q2,L){
  Q=cbind(Q1,Q2)
  hyper=matrix(rexp(3*L,1),ncol=3,nrow=L)
  dist=apply(hyper,1,function(x){U=prior(length(Q1),x);sum((statdir(U)-statdir(Q))^2)})
  return(hyper[which.min(dist),])
}

#ABC step for sigma, Dat is the whole dataset, j the index of the inferred sigma,
#parh1 and parh2 are the vectors of first and second parameter of MA,
#parsigm1 and parsigm2 are the hyperparameters on sigma
#N the number of simulations
ABCsigm <- function(Dat,j,parh1,parh2,parsigm1,parsigm2,N){
  dist=numeric(N)
  n=nrow(Dat)
  w=seq(1,n,3)
  si=1/(rgamma(N,parsigm1,scale=parsigm2))
  dist=lapply(si,function(x){abs(var(simu(parh1[j],parh2[j],x,n)[w])- var(Dat[w,j]))})
  return(si[which.min(dist)])
}

#ABC step for sigma hyperparameter.
#sigm vector of sigmas, Q number of simulations
ABCsigmhaut <- function(sigm,Q){
  dist=numeric(Q)
  n=length(sigm)
  pri=matrix(abs(rcauchy(2*Q,1)),ncol=2,nrow=Q)
  dist=apply(pri,1,function(h){y=1/(rgamma(n,h[1],scale=h[2])) ;return(abs(sum(log(y))-sum(log(sigm)))+abs(sum(y)-sum(sigm)))})
  return(pri[which.min(dist),])
}


#Simulates a chain according to ABCGibbs,
#X is the dataset, N the simulations to infer each parameter
#M simulations for the hyperparameters, P simulations for sigma
#P2 simulations for hypersigma, Npts number of desired points.
#returns a list of matrices : hyperparameters, first parameter, second, sigma, hypersigma
#for parameters and sigma, the columns corresponds to the time series.
Gibbs <- function(X,N,M,P,P2,Npts){
  hyper=matrix(ncol=3,nrow=Npts)
  hyper[1,]=rexp(3,1)
  par1=matrix(NA,ncol=ncol(X),nrow=Npts)
  par2=matrix(NA,ncol=ncol(X),nrow=Npts)
  sigm=matrix(NA,ncol=ncol(X),nrow=Npts)
  hypersigm=matrix(NA,ncol=2,nrow=Npts)
  hypersigm[1,]=abs(rcauchy(2,1))
  sigm[1,]=1/(rgamma(ncol(X),hypersigm[1,1],scale=hypersigm[1,2]))
  d=prior(ncol(X),hyper[1,])
  par1[1,]=d[,1]
  par2[1,]=d[,2]
  m=ncol(X)
  V=X-matrix(colMeans(X),ncol=m,nrow=nrow(X),byrow=T)
  for (i in 2:Npts){
    for (j in 1:ncol(X)){
      U=ABCbas(X[,j],N,sigm[i-1],hyper[i-1,])
      par1[i,j]=U[1]
      par2[i,j]=U[2]
      sigm[i,j]=ABCsigm(X,j,par1[i,],par2[i,],hypersigm[i-1,1],hypersigm[i-1,2],P)
    }
    hyper[i,]=ABChaut(par1[i,],par2[i,],M)
    hypersigm[i,]=ABCsigmhaut(sigm[i,],P2)
  }
  return(list(hyper,par1,par2,sigm,hypersigm))
}

#Simple ABC corresponding, data is the dataset, Npts the number of desired points,
#Ntot number of simulations of the whole hierarchy. JJ renormalisation factor for
#autocorrelation distance in the same order as the colums of data.
#II same for the variance.
#same output as Gibbs
ABCsimple <- function(data,Npts,Ntot,JJ,II){
  hyper=matrix(rexp(3*Ntot,1),ncol=3)
  par1=matrix(NA,ncol=ncol(data),nrow=Ntot)
  par2=matrix(NA,ncol=ncol(data),nrow=Ntot)
  hypersig=matrix(abs(rcauchy(2*Ntot)),ncol=2)
  sig=t(apply(hypersig,1,function(x){1/(rgamma(ncol(data),x[1],x[2]))}))
  dist=numeric(Ntot)
  w=seq(1,nrow(data),3)
  for (i in 1:Ntot){
    V=prior(ncol(data),hyper[i,])
    par1[i,]=V[,1]
    par2[i,]=V[,2]
    u=0
    for (j in 1:ncol(data)){
      z=simu(par1[i,j],par2[i,j],sig[i,j],nrow(data))
      u=u+distautcor(data[,j],z)/JJ[[j]]+abs(var(z[w])- var(data[w,j]))/II[[j]]
    }
    dist[i]=u
  }
  H=order(dist)[1:Npts]
  return(list(hyper[H,],par1[H,],par2[H,],sig[H,],hypersig[H,]))
}

#computes a reference table of size Ntot and returns a reference table, same output as Gibbs
#returns also two vectors dist1 and dist2 that are sample of the a priori distances
#autocorrelation and variance
table <- function(data,Ntot){
  hyper=matrix(rexp(3*Ntot,1),ncol=3)
  par1=matrix(NA,ncol=ncol(data),nrow=Ntot)
  par2=matrix(NA,ncol=ncol(data),nrow=Ntot)
  hypersig=matrix(abs(rcauchy(2*Ntot)),ncol=2)
  sig=t(apply(hypersig,1,function(x){1/(rgamma(ncol(data),x[1],x[2]))}))
  dist1=matrix(NA,nrow=Ntot,ncol=ncol(data))
  dist2=matrix(NA,nrow=Ntot,ncol=ncol(data))
  w=seq(1,nrow(data),3)
  for (i in 1:Ntot){
    V=prior(ncol(data),hyper[i,])
    par1[i,]=V[,1]
    par2[i,]=V[,2]
    u=0
    for (j in 1:ncol(data)){
      z=simu(par1[i,j],par2[i,j],sig[i,j],nrow(data))
      dist1[i,j]=distautcor(data[,j],z)
      dist2[i,j]=abs(var(z[w])- var(data[w,j]))
    }
  }
  return(list(hyper,par1,par2,sig,hypersig,dist1,dist2))
}

#compares two output. BB and CC are the renormalisation constants (same as JJ and II previous)
#U and V are the output, data the dataset.
#returns two matrices, for each the value [i,j] corresponds to the value of the total distance
#for the time serie j with the parameters i.
comparaison <- function(U,V,data,BB,CC){
  n=nrow(U[[1]])
  m=nrow(V[[1]])
  k=ncol(U[[2]])
  d1=matrix(0,nrow=n,ncol=k)
  d2=matrix(0,nrow=n,ncol=k)
  w=seq(1,nrow(data),3)
  for (i in 1:n){
    for (j in 1:k){
      z=simu(U[[2]][i,j],U[[3]][i,j],U[[4]][i,j],nrow(data))
      d1[i,j] = distautcor(data[,j],z)/BB[[j]]+abs(var(z[w])- var(data[w,j]))/CC[[j]]
    }
  }
  for (i in 1:m){
    for (j in 1:k){
      z=simu(V[[2]][i,j],V[[3]][i,j],V[[4]][i,j],nrow(data))
      d2[i,j] =distautcor(data[,j],z)/BB[[j]]+abs(var(z[w])- var(data[w,j]))/CC[[j]]
    }
  }
  return(list(d1,d2))
}

#example on real dataset
YY=read.table("fluxcorrige8GHz.csv",header=T,sep=",")
Tab=table(YY,100000)
Dist=Tab[[6]][which(rowSums(Tab[[6]])<1e30),]
JJ=apply(Dist,2,function(x){quantile(x,0.001)})
Dist2=Tab[[7]][which(rowSums(Tab[[7]])<1e30),]
II=apply(Dist2,2,function(x){quantile(x,0.001)})

flux8G=Gibbs(YY,500,100,100,100,1000)
flux8S=ABCsimple(YY,1000,500000,JJ,II)

Xreel=comparaison(flux8G,flux8S,YY,JJ,II)

#example on toy dataset
V=read.table("toydata.csv")
Tab=table(V,100000)
Dist=Tab[[6]][which(rowSums(Tab[[6]])<1e30),]
JJ=apply(Dist,2,function(x){quantile(x,0.001)})
Dist2=Tab[[7]][which(rowSums(Tab[[7]])<1e30),]
II=apply(Dist2,2,function(x){quantile(x,0.001)})

toyG=Gibbs(V,1000,100,100,100,1000)
toyS=ABCsimple(V,1000,1000000,JJ,II)

Xtoy=comparaison(toyG,toyS,V,JJ,II)


