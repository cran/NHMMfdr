\name{update_delta2}
\alias{update_delta2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intermediate function
}
\description{
Intermediate function
}
\usage{
update_delta2(Z, dgamma, gamma, delta.old, gradient.new, gradient.old)
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
  \item{delta.old}{
%%     ~~Describe \code{delta.old} here~~
}
  \item{gradient.new}{
%%     ~~Describe \code{gradient.new} here~~
}
  \item{gradient.old}{
%%     ~~Describe \code{gradient.old} here~~
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
function(Z, dgamma, gamma, delta.old, gradient.new, gradient.old){

	num.fac <- sum((gradient.new - gradient.old)*gradient.new)
	denom.fac <- sum(gradient.old^2)
	fac <- num.fac/denom.fac
	if(fac<0) fac <- 0
	delta2 <- -gradient.new + fac*delta.old
	return(delta2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
