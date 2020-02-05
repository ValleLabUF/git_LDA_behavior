LDA_behavior_gibbs=function(dat,gamma1,alpha,ngibbs,nmaxclust,nburn,ndata.types){
  ntsegm=nrow(dat)
  
  #separate variables
  y=list()
  ncat.dat=rep(NA,ndata.types)
  for (i in 1:ndata.types){
    nome=paste0('y',i)
    ind=grep(nome,colnames(dat))
    y[[i]]=data.matrix(dat[,ind])
    ncat.dat[i]=length(ind)
  }
  
  #initial values
  phi=z.agg=list()
  for (i in 1:ndata.types){
    phi[[i]]=matrix(1/ncat.dat[i],nmaxclust,ncat.dat[i])
    z.agg[[i]]=array(NA,dim=c(ntsegm,ncat.dat[i],nmaxclust))
  }
  theta=matrix(1/nmaxclust,ntsegm,nmaxclust)
  
  for (j in 1:ndata.types){
    for (i in 1:ntsegm){
      for (k in 1:ncat.dat[j]){
        z.agg[[j]][i,k,]=rmultinom(1,size=y[[j]][i,k],prob=rep(1/nmaxclust,nmaxclust))
      }
    }  
  }
  
  #prepare for gibbs
  store.phi=zeroes=list()
  for (i in 1:ndata.types){
    store.phi[[i]]=matrix(NA,ngibbs,nmaxclust*ncat.dat[i])
    zeroes[[i]]=array(0,c(ntsegm,ncat.dat[i],nmaxclust))
  }
  store.theta=matrix(NA,ngibbs,ntsegm*nmaxclust)
  store.loglikel=rep(NA,1)
  
  #run gibbs sampler
  for (i in 1:ngibbs){
    print(i)
    
    #re-order clusters
    if (i < nburn & i%%50==0){
      med=apply(theta,2,mean)
      ordem=order(med,decreasing=T)
      theta=theta[,ordem]
      
      for (j in 1:ndata.types){
        z.agg[[j]]=z.agg[[j]][,,ordem]
        phi[[j]]=phi[[j]][ordem,]
      }
    }
    
    #sample from FCD's 
    z.agg=sample.z(ntsegm=ntsegm,ncat.dat=ncat.dat,y=y, nmaxclust=nmaxclust,
                   phi=phi,ltheta=log(theta),zeroes=zeroes,ndata.types=ndata.types)
    
    v=sample.v(z.agg=z.agg,gamma1=gamma1,
               ntsegm=ntsegm,ndata.types=ndata.types,nmaxclust=nmaxclust)
    theta=get.theta(v=v,nmaxclust=nmaxclust,ntsegm=ntsegm)
    # theta=theta.true
    
    phi=sample.phi(z.agg=z.agg,alpha=alpha,nmaxclust=nmaxclust,
                   ncat.dat=ncat.dat,ndata.types=ndata.types)
    
    #calculate log-likelihood
    p1=0
    for (j in 1:ndata.types){
      prob1=theta%*%phi[[j]]  
      p1=p1+sum(y[[j]]*log(prob1))
    }
    
    #store results
    for (j in 1:ndata.types){
      store.phi[[j]][i,]=phi[[j]]
    }
    store.theta[i,]=theta
    store.loglikel[i]=p1
  }
  
  list(phi=store.phi,theta=store.theta,
       loglikel=store.loglikel,z.agg=z.agg)  
}