\name{LIS.adjust}
\alias{LIS.adjust}
\title{
Determine statistically significant tests at a user-specified FDR level.
}
\description{
Compute adjusted local index of significance (LIS) and determine the tests which are statistically significant at a user-specified FDR level.

}
\usage{
LIS.adjust(lis, fdr = 0.001, adjust = TRUE)
}
\arguments{
  \item{lis}{
LIS scores from \code{\link{fdr.nhmm}}.
}
  \item{fdr}{
Desirable FDR level. Usually 0.05. 
}
  \item{adjust}{
Logical value. \code{TRUE} if adjusted LIS is to be computed. See \code{details}.
}
}
\details{
If \code{adjust = TRUE}, the adjusted LIS will be computed. This will be useful if the user wants to control at different FDR level. The significant tests corresponds to those adjusted LIS which are less than the FDR level.
}
\value{
\item{States}{A binary valued vector. Takes value 1 if the test or observation is declared to be significant.}
\item{aLIS}{Adjusted LIS.}
}
\references{
P.F. Kuan and D.Y. Chiang (2012). Integrating Prior Knowledge in Multiple Testing Under Dependence with Applications in Detecting Differential DNA Methylation. Biometrics, doi: 10.1111/j.1541-0420.2011.01730.x.}
\author{Pei Fen Kuan}
\seealso{ \code{\link{fdr.nhmm}}}
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

