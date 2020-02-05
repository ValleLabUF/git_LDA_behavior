soma1=apply(z1.agg,c(1,3),sum)
soma2=apply(z2.agg,c(1,3),sum)
soma.fim=soma1+soma2

cumsum1=t(apply(soma1[,nbehav:1],1,cumsum))
cumsum1=cumsum1[,nbehav:1]

cumsum.cpp=CumSumInv(nobs=nobs,nmaxclust=nmaxclust,z=soma1)
unique(cumsum1-cumsum.cpp)