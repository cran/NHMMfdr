forwardbackward1.indep.kernel <-
function(x, ptheta, f0, f1x)
{


## Initialize

NUM<-length(x)

## Densities

f0x<-dnorm(x, f0[1], f0[2])

gamma<-matrix(1:(NUM*2), NUM, 2, byrow=TRUE)

gamma[,1] <- ptheta[1]*f0x/(ptheta[1]*f0x + ptheta[2]*f1x)
gamma[,2] <- 1 - gamma[,1]

lfdr <- gamma[,1]

c0 <- f0x*ptheta[1] + f1x*ptheta[2]

forwardbackward.var<-list(lf=lfdr, pr=gamma,rescale=c0)
return(forwardbackward.var)
  
}

