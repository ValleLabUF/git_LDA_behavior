plot(res$loglikel[100:ngibbs],type='l')

compare=function(true1,estim1){
  rango=range(c(true1,estim1))
  plot(true1,estim1,ylim=rango,xlim=rango)
  lines(rango,rango,col='red')
}

par(mfrow=c(1,1))
theta.estim=matrix(res$theta[ngibbs,],nrow(dat),nmaxclust)
boxplot(theta.estim)

nclust=3
z1.tmp=apply(res$z.agg[[1]],c(1,3),sum)[,1:nclust] 

#find best order
ordem=numeric()
for (i in 1:ncol(z.true[[1]])){
  tmp=rep(NA,ncol(z.true[[1]]))
  for (j in 1:ncol(z1.tmp)){
    tmp[j]=cor(cbind(z1.tmp[,j],z.true[[1]][,i]))[1,2]
  }
  ind=which(tmp==max(tmp))
  ordem=c(ordem,ind)
}

#look at z's
par(mfrow=c(ceiling(ndata.types/2),2),mar=rep(1,4))
for (j in 1:ndata.types){
  z1.tmp=apply(res$z.agg[[j]],c(1,3),sum)[,1:nclust] 
  compare(z.true[[j]],z1.tmp[,ordem])
}

#look at theta's
compare(theta.true,theta.estim[,ordem])

#look at phi's
par(mfrow=c(ceiling(ndata.types/2),2),mar=rep(1,4))
for (j in 1:ndata.types){
  phi1.estim=matrix(res$phi[[j]][ngibbs,],nmaxclust,ncat.data[j])
  compare(phi.true[[j]],phi1.estim[ordem,])
}
