library(ggplot2)
library(ggpubr)

#generation of the dataset
a=2
var=1
sigm=1
Nk=10
thetreel=c(2.774149,2.438900,3.538051,3.134760,3.192536,2.134164,2.460961,2.785255,1.821516,1.999794,2.147105,2.627481,2.273351,1.829135,1.095750,2.192478,1.951274,2.781893,1.795393,2.563236)
obs=c(2.5666002,2.5808230,3.4319618,3.4865059,3.5939341,1.7859019,2.3653864,2.7560831,2.2540691,2.0810668,2.3309626,2.1923769,2.8244147,2.3624706,0.8032528,1.7128994,2.1605168,3.5554632,1.4499120,2.3010078)
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


#estimation of the variance of the estimator
resgibbs <- function(data,thetini,hyperini,sigm,var,nbeps1,nbeps2,nbpts,nbvar,qq)
{
  res1=rep(NA,nbvar)
  res2=matrix(NA,ncol=nbvar,nrow=length(thetini))
  for(i in 1:nbvar)
  {
    U=gibbstot(data,thetini,hyperini,sigm,var,nbeps1,nbeps2,nbpts,qq)
    res1[i]=mean(U[[1]][2:(nbpts+1)])
    res2[,i]=rowSums(U[[2]][,2:(nbpts+1)])/nbpts
  }
  return(list(var(res1),mean(res1),apply(res2,1,var),apply(res2,1,mean)))
}
ressimple <- function(data,sigm,var,nbpts,ntot,nbvar,qq)
{
  res1=rep(NA,nbvar)
  res2=matrix(NA,ncol=nbvar,nrow=length(data))
  for(i in 1:nbvar)
  {
    U=ABCsimple(data,var,sigm,nbpts,ntot,qq)
    res1[i]=mean(U[[1]])
    res2[,i]=rowSums(U[[2]])/nbpts
  }
  return(list(var(res1),mean(res1),apply(res2,1,var),apply(res2,1,mean)))
}

#number of iterations are all integers...
Ntot=floor(1000/(1:400))*1:400
Neps=1:1000
Npts=floor(1000/(1:400))

#variation of Neps
res1simplemeancomp=rep(NA,400)
res2simplemeancomp=matrix(NA,ncol=400,nrow=20)
res1simplevarcomp=rep(NA,400)
res2simplevarcomp=matrix(NA,ncol=400,nrow=20)
for (i in 1:400)
{
  V=ressimple(obs,sigm,var,Npts[i],Ntot[i],100,3)
  res1simplemeancomp[i]=V[[2]]
  res2simplemeancomp[,i]=V[[4]]
  res1simplevarcomp[i]=V[[1]]
  res2simplevarcomp[,i]=V[[3]]
}


res1meancomp=rep(NA,400)
res2meancomp=matrix(NA,ncol=400,nrow=20)
res1varcomp=rep(NA,400)
res2varcomp=matrix(NA,ncol=400,nrow=20)
for (i in 1:400)
{
  V=resgibbs(obs,rep(0,20),0,sigm,var,Neps[i],Neps[i],Npts[i],100,10)
  res1meancomp[i]=V[[2]]
  res2meancomp[,i]=V[[4]]
  res1varcomp[i]=V[[1]]
  res2varcomp[,i]=V[[3]]
}

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


#pour les histogrammes
resG100=gibbstot(obs,rep(0,20),0,sigm,var,100,100,1000,10)
resS100=ABCsimple(obs,var,sigm,1000,100*100,10)


resG1=gibbstot(obs,rep(0,20),0,sigm,var,1,1,1000,10)
resG5=gibbstot(obs,rep(0,20),0,sigm,var,5,5,1000,10)
resG10=gibbstot(obs,rep(0,20),0,sigm,var,10,10,1000,10)
resG1000=gibbstot(obs,rep(0,20),0,sigm,var,1000,1000,1000,10)

exact=gibbsexact(obs,10000,sigm,var,10)

#plusieurs densités superposées

N <- 100

hyper = numeric()
par1 = numeric()
par2 = numeric()
par3 = numeric()
nbr=numeric()
type=numeric()
exact=gibbsexact(obs,1000,sigm,var,10)



for (i in 1:N){
  A=gibbstot(obs,rep(0,20),0,sigm,var,30,30,floor(1000/30),10)
  B=ABCsimple(obs,var,sigm,floor(1000/30),1000,10)
  hyper=c(hyper,A[[1]][5:34],B[[1]])
  par1=c(par1,A[[2]][1,5:34],B[[2]][1,])
  par2=c(par2,A[[2]][2,5:34],B[[2]][2,])
  par3=c(par3,A[[2]][3,5:34],B[[2]][3,])
  nbr=c(nbr,rep(i,length(A[[1]][5:34])+length(B[[1]])))
  type=c(type,c(rep("ABC-Gibbs",length(A[[1]][5:34])),rep("Simple ABC",length(B[[1]]))))
}


gib_curves30 <- data.frame(Group = nbr,
                      mu1 = par1,
                      mu2 = par2,
                      mu3 = par3,
                      hyperparameter = hyper,
                      Method=type)

g=ggplot(data=gib_curves30,aes(x=mu1,color=Method))+theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
for (i in 1:N){
  g=g+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="ABC-Gibbs"),],stat = "density",alpha=0.1)
  g=g+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="Simple ABC"),],stat = "density",alpha=0.1)
}
g=g+geom_line(data=data.frame(mu1=exact[[1]][,1],Method=rep("Exact",1000)),stat = "density")
g
ggsave("mu1comp.pdf",height=10,width=10,unit="cm")

f=ggplot(data=gib_curves30,aes(x=hyperparameter,color=Method))+theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
for (i in 1:N){
  f=f+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="ABC-Gibbs"),],stat = "density",alpha=0.1)
  f=f+geom_line(data=gib_curves30[which(gib_curves30$Group==i & gib_curves30$Method=="Simple ABC"),],stat = "density",alpha=0.1)
}
f=f+geom_line(data=data.frame(hyperparameter=exact[[2]],Method=rep("Exact",1000)),stat = "density")+xlim(c(0,5))
f
ggsave("hypercomp.pdf",height=10,width=10,unit="cm")



#ploting

comparaisonmoyennehyper <- data.frame(value=c(res1meancomp,res1simplemeancomp),Neps=c(1:400,1:400),Method=c(rep("ABC Gibbs",400),rep("Vanilla ABC",400)))
comparaisonvariancehyper <- data.frame(value=c(res1varcomp,res1simplevarcomp),Neps=c(1:400,1:400),Method=c(rep("ABC Gibbs",400),rep("Vanilla ABC",400)))
comparaisonmoyennepar <- data.frame(value=c(res2meancomp[4,],res2simplemeancomp[4,]),Neps=c(1:400,1:400),Method=c(rep("ABC Gibbs",400),rep("Vanilla ABC",400)))
comparaisonvariancepar <- data.frame(value=c(res2varcomp[4,],res2simplevarcomp[4,]),Neps=c(1:400,1:400),Method=c(rep("ABC Gibbs",400),rep("Vanilla ABC",400)))

comptot <- data.frame(value=c(res1meancomp,res1simplemeancomp,res1varcomp,res1simplevarcomp,res2meancomp[1,],res2simplemeancomp[1,],res2varcomp[1,],res2simplevarcomp[1,]),
                      Neps=rep(1:400,8),
                      Method=rep(c(rep("ABC Gibbs",400),rep("Vanilla ABC",400)),4),
                      Obj1=c(rep("Hyperparameter",1600),rep("Parameter",1600)),
                      Obj2=rep(c(rep("Mean",800),rep("Variance",800)),2)
                      )

p =ggplot(comptot)
p=p+theme(axis.title.x=element_blank(),
          axis.title.y=element_blank())
p=p+geom_line(aes(y=value,x=Neps,linetype=Method))+ theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
p=p+facet_grid(cols=vars(Obj1),rows=vars(Obj2),scales = "free")

p1 <- ggplot(data=comparaisonmoyennehyper)+geom_line(aes(y=value,x=Neps,linetype=Method)) + theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
p1 <- p1 + theme(axis.title.x=element_blank(),
                    axis.title.y=element_blank())
p1 <- p1 + ylim(1.5,3.5) + geom_hline(yintercept=mean(thetreel),color='red',linetype=1)+ geom_hline(yintercept=2,color='red',linetype=2)
p2 <- ggplot(data=comparaisonvariancehyper)+geom_line(aes(y=value,x=Neps,linetype=Method))+ theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
p2 <- p2 + theme(axis.title.x=element_blank(),
                 axis.title.y=element_blank())
p3 <- ggplot(data=comparaisonmoyennepar)+geom_line(aes(y=value,x=Neps,linetype=Method))+ theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
p3 <- p3 + theme(axis.title.x=element_blank(),
                 axis.title.y=element_blank())
p3 <- p3 + ylim(2.4,4.5) + geom_hline(yintercept=obs[4],color='red',linetype=1)+ geom_hline(yintercept=thetreel[4],color='red',linetype=2)
p4 <- ggplot(data=comparaisonvariancepar)+geom_line(aes(y=value,x=Neps,linetype=Method))+ theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
p4 <- p4 + theme(axis.title.x=element_blank(),
                 axis.title.y=element_blank())

ggarrange(p1, p2,p3,p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom",align="h", labels="AUTO",hjust=-5,vjust=2)
ggsave("hiercomparaisontoy.pdf",height=20,width=20,unit="cm")

exactp <- data.frame(Parameter=exact[[1]][,1],Hyperparameter=exact[[2]])

compneps = data.frame(Hyperparameter=c(resG1[[1]],resG5[[1]],resG100[[1]],resG100[[1]],resG1000[[1]]),Parameter=c(resG1[[2]][1,],resG5[[2]][1,],resG10[[2]][1,],resG100[[2]][1,],resG1000[[2]][1,]), Neps=c(rep("1",1001),rep("5",1001),rep("10",1001),rep("100",1001),rep("1000",1001)))
compneps$Neps <- factor(compneps$Neps, levels = c("1", "5", "10","100","1000"))
q1 <- ggplot(data=compneps)
q1 <- q1+geom_density(data=compneps,adjust=0.7,aes(x=Parameter,colour=Neps))
q1 <- q1 +theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
q1 <- q1 +geom_density(data=exactp,aes(x=Parameter),linetype=2,show.legend = F)
q1 <- q1 + theme(axis.title.x=element_blank(),
                 axis.title.y=element_blank())
q2 <- ggplot(data=compneps)
q2 <- q2+geom_density(data=compneps,adjust=0.7,aes(x=Hyperparameter,colour=Neps))
q2 <- q2+theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
q2 <- q2 +geom_density(data=exactp,aes(x=Hyperparameter),linetype=2,show.legend = F)
q2 <- q2 + theme(axis.title.x=element_blank(),
                 axis.title.y=element_blank())

ggarrange(q2,q1,ncol=2,nrow=1,common.legend=T,legend="bottom")
ggsave("varnepstoy.pdf",width = 30,height=15,unit="cm")

compdens<-data.frame(Hyperparameter=c(resG100[[1]],resS100[[1]]),Parameter=c(resG100[[2]][1,],resS100[[2]][1,]),Method=c(rep("ABC Gibbs",1001),rep("Simple ABC",1000)))
r1 <- ggplot(data=compdens,aes(x=Parameter))
r1 <- r1+geom_density(aes(x=Parameter,linetype=Method))
#r1<-r1+geom_histogram(position="dodge",aes(y=..density..,fill=Method))
#r1 <- r1+geom_vline(xintercept=thetreel[1])
r1 <- r1+ theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
r1 <- r1 +geom_density(data=exactp,aes(x=Parameter),linetype="dotdash",show.legend = F) 
r1 <- r1 + theme(legend.position = "none",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank())

r2 <- ggplot(data=compdens,aes(x=Hyperparameter))
#r2 <- r2+geom_histogram(position="identity",aes(y=..density..,fill=Method,alpha=.07))
r2 <- r2+geom_density(aes(x=Hyperparameter,linetype=Method))
#r2 <- r2+geom_vline(xintercept=2)
r2 = r2 + theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")
r2 <- r2 +geom_density(data=exactp,aes(x=Hyperparameter),linetype="dotdash",show.legend = F)
r2 <- r2 + theme(legend.position = "none",
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank())


ggarrange(r1,r2,ncol=2,nrow=1,common.legend = T,legend="bottom")
ggsave("compdensity.pdf",width=30,height=15,unit="cm")

#plotte pour plusieurs Nmu des densités
N=c(5,10,30,50,100)

hyper = numeric()
par1 = numeric()
par2 = numeric()
par3 = numeric()
Nb=numeric()
type=numeric()
exact=gibbsexact(obs,1000,sigm,var,10)


for (i in 1:5){
  A=gibbstot(obs,rep(0,20),0,sigm,var,N[i],N[i],floor(1000/N[i]),10)
  B=ABCsimple(obs,var,sigm,floor(1000/N[i]),1000,10)
  hyper=c(hyper,A[[1]],B[[1]])
  par1=c(par1,A[[2]][1,],B[[2]][1,])
  par2=c(par2,A[[2]][2,],B[[2]][2,])
  par3=c(par3,A[[2]][3,],B[[2]][3,])
  Nb=c(Nb,rep(N[i],length(A[[1]])+length(B[[1]])))
  type=c(type,c(rep("ABC-Gibbs",length(A[[1]])),rep("Simple ABC",length(B[[1]]))))
}

pourplot <- data.frame(Value=c(hyper,par1,par2,par3),
                       parameter=c(rep("hyperparameter",length(hyper)),rep("mu1",length(par1)),rep("mu2",length(par2)),rep("mu3",length(par3))),
                       Npts=rep(Nb,4),
                       Method=rep(type,4))

h=ggplot(data=pourplot)
h+geom_density(aes(x=Value,linetype=Method))+geom_vline(data=data.frame(value=c(a,thetreel[1:3]),parameter=c("hyperparameter","mu1","mu2","mu3")),aes(xintercept=value))+
  facet_grid(Npts~parameter)+theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")+
  theme(axis.title.x=element_blank(),
          axis.title.y=element_blank())

ggsave("pleindetrucs.pdf",height=23,width=15,unit="cm")

#version avec nombre de points constant
N=c(5,10,30,50,100)

hyper = numeric()
par1 = numeric()
par2 = numeric()
par3 = numeric()
Nb=numeric()
type=numeric()
exact=gibbsexact(obs,1000,sigm,var,10)


for (i in 1:5){
  A=gibbstot(obs,rep(0,20),0,sigm,var,N[i],N[i],1000,10)
  B=ABCsimple(obs,var,sigm,1000,1000*N[i],10)
  hyper=c(hyper,A[[1]],B[[1]])
  par1=c(par1,A[[2]][1,],B[[2]][1,])
  par2=c(par2,A[[2]][2,],B[[2]][2,])
  par3=c(par3,A[[2]][3,],B[[2]][3,])
  Nb=c(Nb,rep(N[i],length(A[[1]])+length(B[[1]])))
  type=c(type,c(rep("ABC-Gibbs",length(A[[1]])),rep("Simple ABC",length(B[[1]]))))
}

pourplot <- data.frame(Value=c(hyper,par1,par2,par3),
                       parameter=c(rep("hyperparameter",length(hyper)),rep("mu1",length(par1)),rep("mu2",length(par2)),rep("mu3",length(par3))),
                       Npts=rep(Nb,4),
                       Method=rep(type,4))

h=ggplot(data=pourplot)
h+geom_density(aes(x=Value,linetype=Method))+geom_vline(data=data.frame(value=c(a,thetreel[1:3]),parameter=c("hyperparameter","mu1","mu2","mu3")),aes(xintercept=value))+
  facet_grid(Npts~parameter)+theme(panel.background = element_rect(fill = 'white', colour = 'grey'),legend.position="bottom")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

ggsave("pleindetrucsbis.pdf",height=23,width=15,unit="cm")
