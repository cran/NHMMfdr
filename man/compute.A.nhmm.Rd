\name{compute.A.nhmm}
\alias{compute.A.nhmm}
\title{
Compute transition probability (internal function)
}
\description{
An internal function to be used by \code{\link{fdr.nhmm}}.
}
\usage{
compute.A.nhmm(Z, trans.par1, trans.par2, dist.included = TRUE)
}
\arguments{
  \item{Z}{
A matrix of covariates
}
  \item{trans.par1}{
Transition parameters for State 0
}
  \item{trans.par2}{
Transition parameters for State 1
}
\item{dist.included}{
Logical value. \code{dist.included = TRUE} if probe spacing/distance is included in the first column of \code{Z}.
}
}
\value{
Returns transition probability matrix and initial probabilities. 
}
\references{
P.F. Kuan and D.Y. Chiang (2012). Integrating Prior Knowledge in Multiple Testing Under Dependence with Applications in Detecting Differential DNA Methylation. Biometrics, doi: 10.1111/j.1541-0420.2011.01730.x.}
\author{Pei Fen Kuan}
\seealso{ \code{\link{fdr.nhmm}}, \code{\link{LIS.adjust}} }
\examples{

library(NHMMfdr)

################
# Simulate data
################

### simulate covariate and transition prob
NUM1 <- 1000
Z <- rnorm(NUM1)
Z <- matrix(Z,ncol=1)

Z <- apply(Z,2,scale)
trans.par1.true <- c(0,0,0,0)

trans.par2.true <- rnorm(3+dim(Z)[2])

print(trans.par2.true[-1])

A.true <- compute.A.nhmm(Z, trans.par1.true, trans.par2.true, 
	dist.included=FALSE)$A
pii.true <- compute.A.nhmm(Z, trans.par1.true, trans.par2.true,
 	dist.included=FALSE)$pii

### the null distribution
f0 <- c(0, 1)

### the alternative distribution
f1 <- c(3, 1)

### the NHMM data

simdat <- simdata.nhmm(NUM1, pii.true, A.true, f0, 1, f1)

### the observed values
x <- simdat$o

### the unobserved true states
theta1 <- simdat$s

################
# Model fitting
################

fit.nhmm <- fdr.nhmm(x, Z, dist = NULL, log.transform.dist = FALSE,
	alttype ='mixnormal', L=1, maxiter = 100, nulltype = 0, modeltype = 'NHMM', epsilon=1e-4)
	
### checking estimated parameters
print(fit.nhmm$trans.par2[-1])

################
# Adjust LIS
################

LIS.adjust <- LIS.adjust(fit.nhmm$LIS, fdr = 0.1, adjust = TRUE)

### tests which are statistically significant

sig.test <- which(LIS.adjust$States == 1)
length(sig.test)

}
