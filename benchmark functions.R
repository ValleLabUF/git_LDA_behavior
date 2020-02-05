library(microbenchmark)

microbenchmark(
  #sample from FCD's
  z1agg=SampleZ1Agg(nobs=nobs,b1=b1,y1=y1, nmaxclust=nmaxclust,
                  lphi1=log(phi1),ltheta=log(theta),zeroes=zeroes1),
  z2agg=SampleZ2Agg(nobs=nobs,b2=b2,y2=y2, nmaxclust=nmaxclust,
                  lphi2=log(phi2),ltheta=log(theta),zeroes=zeroes2),
  v=sample.v(z1.agg=z1.agg,z2.agg=z2.agg,gamma1=gamma1,
             nobs=nobs,nbehav=nmaxclust),
  theta=get.theta(v=v,nbehav=nmaxclust,nobs=nobs),
  phi1=sample.phi1(z1.agg=z1.agg,alpha=alpha,nbehav=nmaxclust,b1=b1),
  phi2=sample.phi2(z2.agg=z2.agg,alpha=alpha,nbehav=nmaxclust,b2=b2)
)
