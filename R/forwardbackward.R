forwardbackward <-
function(x, pii, A, pc, f0, f1)
{

NUM<-length(x)
L<-length(pc)


f0x<-dnorm(x, f0[1], f0[2])
f1x<-rep(0, NUM)
for (c in 1:L)
{
  f1x<-f1x+pc[c]*dnorm(x, f1[c, 1], f1[c, 2])
}


alpha<-matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)
c0<-rep(0, NUM)

alpha[1, 1]<-pii[1]*f0x[1]
alpha[1, 2]<-pii[2]*f1x[1]
c0[1]<-1/sum(alpha[1, ])
alpha[1, ]<-c0[1]*alpha[1, ]

alpha.tmp <- .C('calAlpha',alpha=as.numeric(alpha),c0=as.numeric(c0),as.numeric(A),as.numeric(f0x),as.numeric(f1x),as.integer(NUM))

alpha <- alpha.tmp$alpha
dim(alpha) <- c(NUM,2)

c0 <- alpha.tmp$c0

beta<-matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)

beta[NUM, 1]<-c0[NUM]
beta[NUM, 2]<-c0[NUM]

beta.tmp <- .C('calBeta',beta=as.numeric(beta),as.numeric(c0),as.numeric(A),as.numeric(f0x),as.numeric(f1x),as.integer(NUM))

beta <- beta.tmp$beta
dim(beta) <- c(NUM,2)

lfdr<-rep(0, NUM)

lfdr.tmp <- .C('calLfdr',as.numeric(alpha),as.numeric(beta),lfdr=as.numeric(lfdr),as.integer(NUM))

lfdr <- lfdr.tmp$lfdr

gamma<-matrix(1:(NUM*2), NUM, 2, byrow=TRUE)
gamma[NUM, ]<-c(lfdr[NUM], 1-lfdr[NUM])
dgamma<-array(rep(0, (NUM-1)*4), c(2, 2, (NUM-1)))

gamma.tmp <- .C('calGamma',as.numeric(alpha),as.numeric(beta),as.numeric(A),as.numeric(f0x),as.numeric(f1x),gamma=as.numeric(gamma),dgamma=as.numeric(dgamma),as.integer(NUM))

gamma <- gamma.tmp$gamma
dgamma <- gamma.tmp$dgamma
dim(gamma) <- c(NUM,2)
dim(dgamma) <- c(2, 2, (NUM-1))

omega<-matrix(rep(0, NUM*L), NUM, L, byrow=TRUE)

 for (c in 1:L)
  { 
    f1c<-dnorm(x, f1[c, 1], f1[c, 2])
    omega[, c]<-gamma[, 2]*pc[c]*f1c/f1x
  }


forwardbackward.var<-list(bw=alpha, fw=beta, lf=lfdr, pr=gamma, ts=dgamma, wt=omega, rescale=c0)
return(forwardbackward.var)
  
}

