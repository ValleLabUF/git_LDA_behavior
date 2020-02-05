rm(list=ls(all=TRUE))
library('MCMCpack')
set.seed(3)

#basic settings
ntsegm=20*40 #overall number of time segments
#be careful with large nbehavior because it is easy to:
#- generate phi1 and phi2 in such a way that things are not clearly separable
#- generate theta in such a way that we don't have approx. pure time segments
nbehavior=3
ncat.data=c(6,10,2,8,5)
ndata.types=length(ncat.data)
nobs.tsegm=100 #number of observations per time segment

#get parameters
phi=list()
for (i in 1:ndata.types){
  phi[[i]]=rdirichlet(nbehavior,alpha=rep(0.1,ncat.data[i]))
}
theta.true=theta=rdirichlet(ntsegm,alpha=rep(0.1,nbehavior))
phi.true=phi

#look at these parameters
boxplot(theta)
apply(theta>0.95,2,mean)

par(mfrow=c(ceiling(nbehavior/2),2),mar=rep(1,4))
for (i in 1:nbehavior) plot(phi.true[[1]][i,],type='h',main=i)
for (i in 1:nbehavior) plot(phi.true[[2]][i,],type='h',main=i)

#create storage space
y.disagg.true=z.true=list()
for (i in 1:ndata.types){
  z.true[[i]]=matrix(0,ntsegm,nbehavior)
  y.disagg.true[[i]]=array(NA,dim=c(ntsegm,ncat.data[i],nbehavior))
}

#generate disaggregated data (already partitioned into different behaviors)
for (j in 1:ndata.types){
  for (t in 1:ntsegm){
    z=rmultinom(1,size=nobs.tsegm,prob=theta[t,]) #number of obs from each behavior
    z.true[[j]][t,]=z
    for (k in 1:nbehavior){
      y.disagg.true[[j]][t,,k]=rmultinom(1,size=z[k],prob=phi[[j]][k,]) #number of obs in each categ from each behavior
    }
  }
}

#aggregate data
y=numeric()
for (j in 1:ndata.types){
  aux=apply(y.disagg.true[[j]],1:2,sum) #sum across behaviors
  nomes=paste0('y',j,'.',1:ncat.data[j])
  colnames(aux)=nomes
  y=cbind(y,aux) 
}

#export fake data
write.csv(y,'fake data.csv',row.names=F)
