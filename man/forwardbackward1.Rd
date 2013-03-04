\name{forwardbackward1}
\alias{forwardbackward1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intermediate function
}
\description{
Intermediate function
}
\usage{
forwardbackward1(x, pii, A, f0, f1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{pii}{
%%     ~~Describe \code{pii} here~~
}
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{f0}{
%%     ~~Describe \code{f0} here~~
}
  \item{f1}{
%%     ~~Describe \code{f1} here~~
}
}
\references{
P.F. Kuan and D.Y. Chiang (2012). Integrating Prior Knowledge in Multiple Testing Under 
Dependence with Applications in Detecting Differential DNA Methylation. Biometrics, 
doi: 10.1111/j.1541-0420.2011.01730.x.

B. Efron (2004). Large-scale simultaneous hypothesis testing: the choice of a null hypothesis. 
Journal of the American Stat. Assoc. 99, 96-104.

W. Sun and T. Cai (2009). Large-scale multiple testing under dependence. 
J. R.Stat. Soc B. 71, 393-424.
}
\author{Pei Fen Kuan}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x, pii, A, f0, f1)
{

## Initialize

NUM<-length(x)

## Densities

f0x<-dnorm(x, f0[1], f0[2])
f1x<-dnorm(x, f1[1], f1[2])

## the backward-forward procedure

alpha<-matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)

c0<-rep(0, NUM)

alpha[1, 1]<-pii[1]*f0x[1]
alpha[1, 2]<-pii[2]*f1x[1]

c0[1]<-1/sum(alpha[1, ])
alpha[1, ]<-c0[1]*alpha[1, ]

alpha.tmp <- .C('calAlpha',alpha = as.numeric(alpha),c0=as.numeric(c0),as.numeric(A),
as.numeric(f0x),as.numeric(f1x),as.integer(NUM))

alpha <- alpha.tmp$alpha
dim(alpha) <- c(NUM,2)

c0 <- alpha.tmp$c0

beta<-matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)

beta[NUM, 1]<-c0[NUM]
beta[NUM, 2]<-c0[NUM]

beta.tmp <- .C('calBeta',beta=as.numeric(beta),as.numeric(c0),as.numeric(A),
as.numeric(f0x),as.numeric(f1x),as.integer(NUM))

beta <- beta.tmp$beta
dim(beta) <- c(NUM,2)

lfdr<-rep(0, NUM)

lfdr.tmp <- .C('calLfdr',as.numeric(alpha),as.numeric(beta),lfdr = as.numeric(lfdr),
as.integer(NUM))

lfdr <- lfdr.tmp$lfdr
 
gamma<-matrix(1:(NUM*2), NUM, 2, byrow=TRUE)
gamma[NUM, ]<-c(lfdr[NUM], 1-lfdr[NUM])
dgamma<-array(rep(0, (NUM-1)*4), c(2, 2, (NUM-1)))

gamma.tmp <- .C('calGamma',as.numeric(alpha),as.numeric(beta),as.numeric(A),as.numeric(f0x),
as.numeric(f1x),gamma=as.numeric(gamma),dgamma=as.numeric(dgamma),as.integer(NUM))

gamma <- gamma.tmp$gamma
dgamma <- gamma.tmp$dgamma
dim(gamma) <- c(NUM,2)
dim(dgamma) <- c(2, 2, (NUM-1))

forwardbackward.var<-list(bw=alpha, fw=beta, lf=lfdr, pr=gamma, ts=dgamma, rescale=c0)
return(forwardbackward.var)
  
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
