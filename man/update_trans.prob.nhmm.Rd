\name{update_trans.prob.nhmm}
\alias{update_trans.prob.nhmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intermediate function
}
\description{
Intermediate function
}
\usage{
update_trans.prob.nhmm(Z, dgamma, gamma, trans.par1, trans.par2, iter.conj.grad = 10, 
dist.included = TRUE)
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
  \item{iter.conj.grad}{
%%     ~~Describe \code{iter.conj.grad} here~~
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
function(Z, dgamma, gamma, trans.par1, trans.par2,iter.conj.grad=10, dist.included=TRUE){

	trans.par1.old <- trans.par1
	trans.par2.old <- trans.par2

	gradient.old <- compute.gradient(Z, dgamma, gamma, trans.par1.old, trans.par2.old,
	dist.included=dist.included)
	delta2.old <- -gradient.old

	n.conj.grad = 0
	while(n.conj.grad <= iter.conj.grad){

	tmp <- try(find.nu(Z, dgamma, gamma, c(0,0,0,0),delta2.old,trans.par1.old, trans.par2.old,
	dist.included=dist.included))

	if(is.na(tmp$nu)){
		break		
	}
	trans.par1.new <- tmp$trans.par1
	trans.par2.new <- tmp$trans.par2

	gradient.new <- compute.gradient(Z, dgamma, gamma,trans.par1.new, trans.par2.new,
	dist.included=dist.included)
	delta2.new <- update_delta2(Z, dgamma, gamma, delta2.old, gradient.new, gradient.old)

	trans.par1.old <- trans.par1.new
	trans.par2.old <- trans.par2.new
        if(dist.included==TRUE&trans.par2.new[4]<0)trans.par2.new[4] <- 0.5
	gradient.old <- gradient.new
	delta2.old <- delta2.new
	n.conj.grad <- n.conj.grad + 1
	
	}

	if(!is.na(tmp$nu)){
		
	tmp.trans.prob <- compute.A.nhmm(Z, trans.par1.new, trans.par2.new,
	dist.included=dist.included)
	pii.new <- tmp.trans.prob$pii
	A.new <- tmp.trans.prob$A
	return(list(A = A.new, pii = pii.new, trans.par1 = trans.par1.new, 
	trans.par2 = trans.par2.new))
	}else return(1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
