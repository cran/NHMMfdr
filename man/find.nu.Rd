\name{find.nu}
\alias{find.nu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intermediate function
}
\description{
Intermediate function
}
\usage{
find.nu(Z, dgamma, gamma, delta1, delta2, trans.par1, trans.par2,dist.included = TRUE)
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
  \item{delta1}{
%%     ~~Describe \code{delta1} here~~
}
  \item{delta2}{
%%     ~~Describe \code{delta2} here~~
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
function(Z, dgamma, gamma, delta1, delta2, trans.par1, trans.par2, dist.included=TRUE){

	p <- dim(Z)[2]
	# delta1 = (lambda1_H, sigma11_H, sigma21_H, rho1_H)
	delta_H <- array(0,c(2,2,c(dim(Z)[1]-1)))
	delta1_H0 <- delta1[1] + sum(delta1[-c(1:3)]*Z[1,])
	delta2_H0 <- delta2[1] + sum(delta2[-c(1:3)]*Z[1,])
    
	tmp11 <- t(delta1[-c(1:3)]*t(Z[-1,]))
	tmp21 <- t(delta1[-c(1:3)]*t(Z[-1,]))
	tmp12 <- t(delta2[-c(1:3)]*t(Z[-1,]))
	tmp22 <- t(delta2[-c(1:3)]*t(Z[-1,]))

        if(dist.included==TRUE){
          delta1_H0 <- delta1[1] + sum(delta1[-c(1:4)]*Z[1,-1])
          delta2_H0 <- delta2[1] + sum(delta2[-c(1:4)]*Z[1,-1])
    
          #tmp11 <- t(c(delta1[4],delta1[-c(1:4)])*t(Z[-1,]))
          tmp21 <- t(c(-delta1[4],delta1[-c(1:4)])*t(Z[-1,]))
          #tmp12 <- t(c(delta2[4],delta2[-c(1:4)])*t(Z[-1,]))
          tmp22 <- t(c(-delta2[4],delta2[-c(1:4)])*t(Z[-1,]))
        }
        
	denom11 = denom21 = denom12 = denom22 <- rep(0,dim(Z)[1]-1)

	for(i in 1:p){
		denom11 <- denom11 + tmp11[,i]
		denom21 <- denom21 + tmp21[,i]
		denom12 <- denom12 + tmp12[,i]
		denom22 <- denom22 + tmp22[,i]
	}
	
	delta_H[1,1,] <- delta1[2] + denom11
	delta_H[2,1,] <- delta1[3] + denom21
	delta_H[1,2,] <- delta2[2] + denom12
	delta_H[2,2,] <- delta2[3] + denom22	

	
	nu.new <- 0
	ptol <- 1
	iter <- 0

while(ptol > 1e-6){	
	iter <- iter + 1
	
	trans.par1.new <- trans.par1 + nu.new*delta1
	trans.par2.new <- trans.par2 + nu.new*delta2

	tmp.trans.prob <- compute.A.nhmm (Z, trans.par1.new, trans.par2.new,
	dist.included=dist.included)
	pii <- tmp.trans.prob$pii
	A <- tmp.trans.prob$A
				
	f <- sum(c(delta1_H0,delta2_H0)*(gamma[1,] - pii))
	for(i in 1:2){
		for(j in 1:2){	
			f <- f + sum(delta_H[i,j,]*(dgamma[i,j,] - gamma[-dim(gamma)[1],i]*A[i,j,]))
		}
	}

     fprime <- -sum(c(delta1_H0,delta2_H0)^2*pii*(1 - pii))
	for(i in 1:2){
		for(j in 1:2){
			fprime <- fprime - 
			sum(delta_H[i,j,]^2*gamma[-dim(gamma)[1],i]*A[i,j,]*(1 - A[i,j,]))
		}
	}
	nu.old <- nu.new
	nu.new <- nu.old - f/fprime
	ptol <- abs(f)
	if(is.na(ptol)||iter > 100){
		nu.new <- NaN
		break
	}
}
	trans.par1.new <- trans.par1 + nu.new*delta1
	trans.par2.new <- trans.par2 + nu.new*delta2
	
	return(list(nu=nu.new,trans.par1=trans.par1.new,trans.par2= trans.par2.new))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
