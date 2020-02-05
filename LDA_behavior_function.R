sample.z=function(ntsegm,nbins,y,nmaxclust,phi,ltheta,zeroes,ndata.types){
  z.agg=list()
  for (i in 1:ndata.types){
    tmp=SampleZAgg(ntsegm=ntsegm,b1=nbins[i],y1=y[[i]],nmaxclust=nmaxclust,
                   lphi1=log(phi[[i]]), ltheta=ltheta,zeroes=zeroes[[i]])  
    z.agg[[i]]=tmp$Z1Agg
  }
  z.agg
}
#-----------------------------------
sample.v=function(z.agg,gamma1,ntsegm,ndata.types,nmaxclust){
  soma.fim=matrix(0,ntsegm,nmaxclust)
  cumsum.fim=matrix(0,ntsegm,nmaxclust-1)
  for (i in 1:ndata.types){
    soma=apply(z.agg[[i]],c(1,3),sum)
    soma.fim=soma.fim+soma
    tmp=CumSumInv(ntsegm=ntsegm,nmaxclust=nmaxclust,z=soma)
    cumsum.fim=cumsum.fim+tmp[,-1]
  }
  
  v=matrix(NA,ntsegm,nmaxclust-1)
  for (i in 1:(nmaxclust-1)){
    v[,i]=rbeta(ntsegm,soma.fim[,i]+1,cumsum.fim[,i]+gamma1)
  }
  cbind(v,1)
}
#-----------------------------------
get.theta=function(v,nmaxclust,ntsegm){
  theta=matrix(NA,ntsegm,nmaxclust)
  theta[,1]=v[,1]
  prod=rep(1,ntsegm)
  for (i in 2:nmaxclust){
    prod=prod*(1-v[,i-1])
    theta[,i]=v[,i]*prod
  }
  theta
}
#-----------------------------------
sample.phi=function(z.agg,alpha,nmaxclust,nbins,ndata.types){
  phi=list()
  for (j in 1:ndata.types){
    soma=t(apply(z.agg[[j]],2:3,sum))  
    tmp=matrix(NA,nmaxclust,nbins[j])
    for (i in 1:nmaxclust){
      tmp[i,]=rdirichlet(1,soma[i,]+alpha)
    }
    phi[[j]]=tmp
  }
  phi
}

