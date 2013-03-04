\name{simdata1.nhmm}
\alias{simdata1.nhmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intermediate function
}
\description{
Intermediate function
}
\usage{
simdata1.nhmm(NUM, pii, A, f0, f1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{NUM}{
%%     ~~Describe \code{NUM} here~~
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
Dependence with Applications in Detecting Differential DNA Methylation. Biometrics, doi: 10.1111/j.1541-0420.2011.01730.x.

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
function(NUM, pii, A, f0, f1)
{


theta<-rep(0, NUM)
x<-rep(0, NUM)

## generating the states
 # initial state
theta[1]<-rbinom(1, 1, pii[2])
 # other states
for (i in 2:NUM)
{
  if (theta[i-1]==0)
     theta[i]<-rbinom(1, 1, A[1, 2, i-1])
  else
     theta[i]<-rbinom(1, 1, A[2, 2, i-1])
}

## generating the observations
for (i in 1:NUM)
{
  if (theta[i]==0)
  {
    x[i]<-rnorm(1, mean=f0[1], sd=f0[2])
  }
  else
  { 
    x[i]<-rnorm(1, mean=f1[1], sd=f1[2])
  }
}
data<-list(s=theta, o=x)
return (data)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
