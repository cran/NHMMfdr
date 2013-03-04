\name{compute.gradient}
\alias{compute.gradient}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intermediate function
}
\description{
Intermediate function
}
\usage{
compute.gradient(Z, dgamma, gamma, trans.par1, trans.par2, dist.included = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{dgamma}{
%%     ~~Describe \code{dgamma} here~~
}
  \item{gamma}{
%%     ~~Describe \code{gamma} here~~
}
  \item{trans.par1}{
%%     ~~Describe \code{trans.par1} here~~
}
  \item{trans.par2}{
%%     ~~Describe \code{trans.par2} here~~
}
\item{dist.included}{
%%     ~~Describe \code{dist.included} here~~
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
function(Z, dgamma, gamma, trans.par1, trans.par2, dist.included=TRUE){

	gradient <- rep(0,3+dim(Z)[2])

	tmp.trans.prob <- compute.A.nhmm(Z, trans.par1, trans.par2, dist.included=dist.included)

	pii <- tmp.trans.prob$pii
	A <- tmp.trans.prob$A
	
	gradient[1] <- gamma[1,2] - pii[2]
	gradient[2] <- sum(dgamma[1,2,] - gamma[-dim(gamma)[1],1]*A[1,2,])
	gradient[3] <- sum(dgamma[2,2,] - gamma[-dim(gamma)[1],2]*A[2,2,])
		
	tmp <- rep(0,dim(dgamma)[3])
        
	for(i in 1:2){
		tmp <- tmp + gamma[-dim(Z)[1],i]*A[i,2,] 
	}
		gradient[-c(1:3)] <- (gamma[1,2] - pii[2])*Z[1,] 
		+ apply(matrix((gamma[-1,2] - tmp)*Z[-1,],ncol=dim(Z)[2]),2,sum)

        if(dist.included==TRUE){
          gradient[4] <- (gamma[1,2] - pii[2])*Z[1,1] + sum(dgamma[1,2,]*Z[-1,1]) 
          - sum(dgamma[2,2,]*Z[-1,1]) - sum(gamma[-dim(Z)[1],1]*A[1,2,]*Z[-1,1]) 
          + sum(gamma[-dim(Z)[1],2]*A[2,2,]*Z[-1,1])
        }
        
	return(-gradient)

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
