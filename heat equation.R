
#a coeff, u état avant, dt pas de temps
pas <- function(a,u,dt)
{
  n=length(u)
  b=1/dt*c(1/(3*n),1/(6*n),rep(0,n-3),1/(6*n))
  Q=suppressWarnings(matrix(b[matrix(1:n,n+1,n+1,byrow=T)[c(1,n:2),1:n]],n,n))
  Q[1,n]=1/(6*n*dt)
  Q[n,1]=1/(6*n*dt)
  c=numeric(n^2)
  c[seq(2,n^2,n+1)]=-a[2:n]
  c[seq(n+1,n^2,n+1)]=-a[2:n]
  K=matrix(c,ncol=n,byrow=T)
  K[1,n]=-a[1]
  K[n,1]=-a[1]
  diag(K)=-rowSums(K)
  return(solve(Q+K,Q%*%u))
}

#simule une trajectoire
simu <- function(a,u0,dt,k){
  n=length(a)
  R=matrix(NA,ncol=k,nrow=n)
  R[,1]=u0
  for (i in 2:10){
    R[,i]=pas(a,R[,i-1],0.1)
  }
  return(R+rnorm(n*k,0,0.1))
}

gibbsparam <-function(data,t0,k,j,a,dt,neps,obs)
{
  dist=rep(NA,neps)
  par=runif(neps,1,2)
  for (i in 1:neps)
  {
    ac=a
    ac[j]=par[i]
    dist[i]=sum(abs(simu(ac,t0,dt,k)[obs,]-data[obs,]))
  }
  return(par[which.min(dist)])
}

#gibbsabc
gibbstot <- function(data,t0,k,dt,neps,npts)
{
  p=nrow(data)
  A=matrix(nrow=npts,ncol=p)
  A[1,]=runif(p,1,2)
  At=A[1,]
  for (i in 2:npts)
  {
    for (j in 1:p)
    {
      At[j]=gibbsparam(data,t0,k,j,At,dt,neps,(j-2):(j+1)%%(p)+1)
    }
    A[i,]=At
  }
  return(A)
}

#simpleABC
abcsimple <- function(data,t0,k,dt,ntot,npts)
{
  dist=rep(NA,ntot)
  p=nrow(data)
  A=matrix(runif(ntot*p,1,2),nrow=ntot,ncol=p)
  for (i in 1:ntot)
  {
    dist[i]=sum(abs(simu(A[i,],t0,dt,k)-data))
  }
  return(A[order(dist)[1:npts],])
}


A1=c(0.7543740135,0.4924927489,0.6960668999,0.6353256113,0.5189776772,0.0439311250,0.3174562373,0.4521930572,0.6808836150,0.5748584443,0.1353280602,0.1143711742,0.252031299,0.2166931028,0.0001339563,0.7680274749,0.0338563418,0.9298740495,0.0814381025,0.6435299462)
t0=rep(c(1,5,10,14,7),4)
data1=simu(A1,t0,.1,10)

estimation <- function(data,neps,npts,dt,qq)
{
  R1=matrix(NA,ncol=nrow(data),nrow=qq)
  R2=matrix(NA,ncol=nrow(data),nrow=qq)
  t0=data[,1]
  k=ncol(data)
  for (i in 1:qq)
  {
    R1[i,]=colMeans(gibbstot(data,t0,k,dt,neps,npts))
    R2[i,]=colMeans(abcsimple(data,t0,k,dt,20*neps*npts,npts))
  }
  return(list(apply(R1,2,mean),apply(R2,2,mean),apply(R1,2,var),apply(R2,2,var)))
}

#testing mean and variance of the a posteriori estimates
Neps=seq(1,50,5)
Ntot=floor(500/Neps)*Neps
Npts=floor(500/Neps)

resmean1G=matrix(NA,ncol=10,nrow=20)
resmean1S=matrix(NA,ncol=10,nrow=20)
resvar1G=matrix(NA,ncol=10,nrow=20)
resvar1S=matrix(NA,ncol=10,nrow=20)

for (i in 1:10)
{
  U=estimation(data1,Neps[i],Npts[i],.1,100)
  resmean1G[,i]=U[[1]]
  resmean1S[,i]=U[[2]]
  resvar1G[,i]=U[[3]]
  resvar1S[,i]=U[[4]]
  print(i)
}

#comparing density
U=gibbstot(data1,t0,10,.1,40,1000)
V=abcsimple(data1,t0,10,.1,800000,1000)

#save.image("190523chaleur.RData")
#===============================================
#plotting the results
#================================================


res=data.frame(Estimator=c(resmean1G[1,],resmean1S[1,]),Variance=c(resvar1G[1,],resvar1S[1,]),Neps=rep(Neps,2),Method=c(rep("ABC Gibbs",10),rep("Simple ABC",10)))
q1=ggplot(data=res)
q1 <- q1+geom_line(aes(y=Estimator,x=Neps,linetype=Method))+ theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
q1 <- q1+geom_hline(yintercept = A1[1])
q2=ggplot(data=res)
q2 <- q2+geom_line(aes(y=Variance,x=Neps,linetype=Method))+ theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
ggarrange(q1,q2,common.legend = T,legend="bottom")
ggsave("nonhierevol.pdf",width = 20,height=10,unit="cm")

histo=data.frame(Parameter=c(U[,1],V[,1]),Method=c(rep("ABC Gibbs",1000),rep("Simple ABC",1000)))
p <- ggplot(data=histo,aes(x=Parameter))
p<-p+geom_histogram(position="dodge",aes(y=..density..,alpha=.7))
p <- p+geom_vline(xintercept=A1[1])+theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
p <- p + theme(legend.position = "none",
               axis.title.x=element_blank(),
               axis.title.y=element_blank(),)
p = p +facet_grid(cols=vars(Method))
print(p)
ggsave("nonhiercompdensity.pdf",width=20,height=10,unit="cm")

