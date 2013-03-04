\name{fdr.nhmm}
\alias{fdr.nhmm}
\title{
FDR control based for multiple hypothesis tests under dependence
}
\description{
Main function which fits the NHMM and HMM based fdr control. It offers numerous options for FDR control as described in Kuan et al. (2011).

}
\usage{
fdr.nhmm(x, Z = NULL, dist = NULL, log.transform.dist = TRUE, alttype = "kernel", L = 2, 
maxiter = 1000, nulltype = 0, modeltype = "NHMM", symmetric = FALSE,epsilon =1e-04)
}
\arguments{
  \item{x}{
A vector of the observed data. It is assumed that the data has been ordered by time or position.
}
\item{Z}{
A matrix of covariates EXCLUDING spacing/distance between probes. Each column corresponds to one covariate. If there is only one covariate, define \code{Z} to be a matrix with one column. If fitting a homogeneous HMM, i.e., no covariate is required, set \code{Z = NULL}.
}
\item{dist}{
A vector of spacing/distance between probes. If fitting a homogeneous HMM, i.e., no covariate is required, set \code{dist = NULL}.
}
\item{log.transform.dist}{
Logical value. \code{TRUE} is log2(Distance+2) is to be applied to the \code{dist} covariate. Recommended for numerical stability.
}
\item{alttype}{
Type of estimating for alternative hypothesis f1(x). Available choices are "\code{mixnormal}" for gaussian mixtures or "\code{kernel}" for non-parametric kernel density estimates.
}
  \item{L}{
Number of mixture component for \code{alttype = "mixnormal"}. If  \code{alttype = "kernel"}, this is irrelevant.
}
\item{maxiter}{
Maximum iterations in the EM algorithm to speed up computation.
}
\item{nulltype}{
Type of null hypothesis assumed in estimating f0(z), for use in the fdr calculations. \code{0} is the theoretical null N(0,1), \code{1} is maximum likelihood estimation, \code{2} is central matching estimation. This method is imported from R package \code{locfdr}.\\ NOTE: Recommended to use nulltype = 0, i.e., theoretical null to avoid double correction of correlation structure.
}
\item{modeltype}{
Types of dependence structure. Available choices are "\code{Indep}", "\code{HMM}" or "\code{NHMM}". See \code{details}.
}
\item{symmetric}{
Logical value. \code{TRUE} if the alternative hypothesis f1(x) is assumed to be symmetrical. Note that this option is only available for \code{alttype = "kernel"}.
}
\item{epsilon}{
Convergence control in the EM algorithm.
}
}
\details{
\code{modeltype = "Indep"} assumes the tests are independent. \code{modeltype = "HMM"} assumes the dependence structure follows a homogeneous HMM based on Sun and Cai (2009). 
}
\value{

\item{LIS}{Local index of significance (analog of p-values). This will be used for FDR control in \code{\link{LIS.adjust}}.}
\item{BIC}{Bayesian information criterion (BIC).}
\item{pii}{Initial probabilities.}
\item{A}{Transition probability matrix.}
\item{f0}{Null hypothesis.}
\item{f1}{Alternative hypothesis.}
\item{logL}{Log likelihood.}
\item{trans.par2}{Transition parameters for State 1.}
}
\references{
P.F. Kuan and D.Y. Chiang (2012). Integrating Prior Knowledge in Multiple Testing Under Dependence with Applications in Detecting Differential DNA Methylation. Biometrics, doi: 10.1111/j.1541-0420.2011.01730.x.

B. Efron (2004). Large-scale simultaneous hypothesis testing: the choice of a null hypothesis. Journal of the American Stat. Assoc. 99, 96-104.

W. Sun and T. Cai (2009). Large-scale multiple testing under dependence. J. R.Stat. Soc B. 71, 393-424.
}
\author{Pei Fen Kuan}
\seealso{ \code{\link{LIS.adjust}} }
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

