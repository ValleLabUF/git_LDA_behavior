set.seed(1)
tmp=SampleZ2Agg(nobs=nobs,b2=b2,y2=y2, nmaxclust=nmaxclust,
                lphi2=log(phi2),ltheta=log(theta))
z2.agg=tmp$Z2Agg


set.seed(1)
z2.agg.denis=sample.z2.agg(lphi2=log(phi2),ltheta=log(theta),y2=y2,
                     nobs=nobs,b2=b2,nbehav=nmaxclust)

hist(unique(z2.agg-z2.agg.denis))
range(z2.agg-z2.agg.denis)