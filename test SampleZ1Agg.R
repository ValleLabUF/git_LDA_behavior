set.seed(1)
z1.agg=sample.z1.agg(lphi1=log(phi1),ltheta=log(theta),y1=y1,
                     nobs=nobs,b1=b1,nbehav=nmaxclust)

set.seed(1)
tmp=SampleZ1Agg(nobs=nobs,b1=b1,y1=y1, nmaxclust=nmaxclust,
                lphi1=log(phi1),ltheta=log(theta), zeroes=array(0,c(nobs,b1,nmaxclust)))
z1.agg.cpp=tmp$Z1Agg

hist(unique(z1.agg-z1.agg.cpp))
range(z1.agg-z1.agg.cpp)
diff=z1.agg-z1.agg.cpp

for (i in 1:nobs){
  diff=z1.agg[i,,]-z1.agg.cpp[i,,]
  cond=diff < -1
  if (sum(cond)>0) break
}
diff=z1.agg[i,,]-z1.agg.cpp[i,,]
diff=z1.agg[10,1,]-z1.agg.cpp[10,1,]
y1[10,1]