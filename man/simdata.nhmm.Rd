\name{simdata.nhmm}
\alias{simdata.nhmm}
\title{
Simulate from a non-homogeneous hidden Markov Model (NHMM).
}
\description{
Function to simulate observations from a NHMM. 
}
\usage{
simdata.nhmm(NUM, pii, A, f0, pc, f1)
}
\arguments{
  \item{NUM}{
Number of observations.
}
  \item{pii}{
Initial probabilities.
}
  \item{A}{
Transition probability matrix.
}
  \item{f0}{
Null distribution.
}
  \item{pc}{
Mixture proportion for alternative distribution.
}
  \item{f1}{
Alternative distribution.
}
}
\value{
\item{s}{Simulated hidden states.}
\item{o}{Simulated observations.}
}
\references{
P.F. Kuan and D.Y. Chiang (2011). Integrating Prior Knowledge in Multiple Testing Under Dependence with Applications in Detecting Differential DNA Methylation. Biometrics, doi: 10.1111/j.1541-0420.2011.01730.x.}
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

A.true <- compute.A.nhmm(Z, trans.par1.true, trans.par2.true)$A
pii.true <- compute.A.nhmm(Z, trans.par1.true, trans.par2.true)$pii

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

fit.nhmm <- fdr.nhmm(x, Z, alttype ='mixnormal', L=1, maxiter = 100, nulltype = 0, 
modeltype = 'NHMM', epsilon=1e-4)

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

