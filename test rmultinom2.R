nclass=5
n=1000
# prob=runif(nclass)
# prob=prob/sum(prob)
prob=c(0,0,0,0.01,0.99)

z=rmultinom2(prob=prob,n=n,randu=rep(1,n),nmaxclust=nclass)

prob.estim=z/sum(z)
rango=range(c(prob,prob.estim))
plot(prob,prob.estim,xlim=rango,ylim=rango)
lines(rango,rango)
#----------------------------------
#Do I get the same set of number numbers?
set.seed(1)
z=rmultinom2(prob=prob,n=n,randu=runif(n),nmaxclust=nclass)
set.seed(1)
z1=rmultinom(1,size=n,prob=prob)