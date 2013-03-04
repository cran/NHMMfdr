forwardbackward.indep <-
function(x, ptheta, pc, f0, f1)
{


## Initialize

NUM<-length(x)
L<-length(pc)

## Densities

f0x<-dnorm(x, f0[1], f0[2])
f1x<-rep(0, NUM)
for (c in 1:L)
{
  f1x<-f1x+pc[c]*dnorm(x, f1[c, 1], f1[c, 2])
}

gamma<-matrix(1:(NUM*2), NUM, 2, byrow=TRUE)

gamma[,1] <- ptheta[1]*f0x/(ptheta[1]*f0x + ptheta[2]*f1x)
gamma[,2] <- 1 - gamma[,1]

lfdr <- gamma[,1]
omega<-matrix(rep(0, NUM*L), NUM, L, byrow=TRUE)

 for (c in 1:L)
  { 
    f1c<-dnorm(x, f1[c, 1], f1[c, 2])
    omega[, c]<-gamma[, 2]*pc[c]*f1c/f1x
  }


c0 <- f0x*ptheta[1] + f1x*ptheta[2]

forwardbackward.var <- list(lf=lfdr, pr=gamma, wt=omega,rescale=c0)
return(forwardbackward.var)
  
}
