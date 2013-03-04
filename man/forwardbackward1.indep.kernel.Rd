\name{forwardbackward1.indep.kernel}
\alias{forwardbackward1.indep.kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intermediate function
}
\description{
Intermediate function
}
\usage{
forwardbackward1.indep.kernel(x, ptheta, f0, f1x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{ptheta}{
%%     ~~Describe \code{ptheta} here~~
}
  \item{f0}{
%%     ~~Describe \code{f0} here~~
}
  \item{f1x}{
%%     ~~Describe \code{f1x} here~~
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
