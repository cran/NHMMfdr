\name{em.indep}
\alias{em.indep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intermediate function
}
\description{
Intermediate function
}
\usage{
em.indep(x, alttype = "mixnormal", L = 2, maxiter = 1000, nulltype = 2, 
symmetric = FALSE,epsilon = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{alttype}{
%%     ~~Describe \code{alttype} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{maxiter}{
%%     ~~Describe \code{maxiter} here~~
}
  \item{nulltype}{
%%     ~~Describe \code{nulltype} here~~
}
  \item{symmetric}{
%%     ~~Describe \code{symmetric} here~~
}
\item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
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
function(x, alttype='mixnormal', L=2, maxiter=1000, nulltype = 2, symmetric = FALSE, 
epsilon=1e-4)
{

NUM<-length(x)

ptol<-epsilon
niter<-0

# Assuming it will converge 
converged=TRUE

### initializing model parameters
if(alttype == 'kernel'){

f0.new<-c(2, 1)

locfdr_p0 <- locfdr(x,plot=0)
if(nulltype == 0){
	f0.new <- c(0,1)
}
if(nulltype == 1){
	f0.new <- c(locfdr_p0$fp0[3,1],locfdr_p0$fp0[3,2])
}
if(nulltype == 2){
	f0.new <- c(locfdr_p0$fp0[5,1],locfdr_p0$fp0[5,2])
}

f1.new <- 0.5*dnorm(x,4,1)+0.5*dnorm(x,-4,1)
ptheta.new <- c(0.5,0.5)

diff<-1

### The E-M Algorithm

while(diff>ptol && niter<maxiter)
{

niter<-niter+1

ptheta.old <- ptheta.new
f0.old <- f0.new
f1.old <- f1.new

## updating the weights and probabilities of hidden states

forwardbackward.res <- forwardbackward1.indep.kernel(x, ptheta.old, f0.old, f1.old)

gamma <- forwardbackward.res$pr
c0 <- forwardbackward.res$rescale

## updating the parameter estimates

ptheta.new <- apply(gamma,2,sum)/NUM

q5 <- sum(gamma[, 1]*x)
mu0 <- q5/sum(gamma[, 1])

q6 <- sum(gamma[, 1]*(x-mu0)*(x-mu0))
sd0 <- sqrt(q6/sum(gamma[, 1]))

f0.new<-c(mu0, sd0)

if(nulltype == 0){
	f0.new <- c(0,1)
}
if(nulltype == 1){
	f0.new <- c(locfdr_p0$fp0[3,1],locfdr_p0$fp0[3,2])
}
if(nulltype == 2){
	f0.new <- c(locfdr_p0$fp0[5,1],locfdr_p0$fp0[5,2])
}

if(symmetric == FALSE){
	kern.f1 <- density(x,weights=gamma[,2]/sum(gamma[,2]))
	f1.new <- approx(kern.f1$x, kern.f1$y, x, rule = 2, ties="ordered")$y
}

if(symmetric == TRUE){
	kern.f1 <- density(c(x,2*f0.new[1]-x),
	weights=c(gamma[,2],gamma[,2])/sum(c(gamma[,2],gamma[,2])))
	f1.new <- approx(kern.f1$x, kern.f1$y, x, rule = 2, ties="ordered")$y
}

df2 <- abs(f1.old-f1.new)
diff <- max(df2)

if (is.na(diff)) {
   converged=FALSE;
   break;
}

}

 lfdr <- gamma[, 1]
if (converged) {
 logL<- sum(log(c0))
if (nulltype > 0) {
	BIC <- logL-(3+2-2)*log(NUM)/2 
} else {
	BIC <- logL-(3-2)*log(NUM)/2 
}

 em.var<-list(ptheta=ptheta.new, f0=f0.new, f1=kern.f1, LIS=lfdr, logL=logL, BIC=BIC, 
 ni=niter, converged=converged) 
} else {
 BIC<- logL<- (-Inf)
 em.var<-list(ptheta=ptheta.old, f0=f0.old, f1=kern.f1, LIS=lfdr, logL=logL, BIC=BIC, 
 ni=niter, converged=converged)
}

}

if(alttype == 'mixnormal'){

########
# L=1
########

if (L==1)
{

f0.new<-c(2, 1)

locfdr_p0 <- locfdr(x,plot=0)
if(nulltype == 0){
	f0.new <- c(0,1)
}
if(nulltype == 1){
	f0.new <- c(locfdr_p0$fp0[3,1],locfdr_p0$fp0[3,2])
}
if(nulltype == 2){
	f0.new <- c(locfdr_p0$fp0[5,1],locfdr_p0$fp0[5,2])
}

f1.new<-c(4, 1)
ptheta.new <- c(0.5,0.5)

diff<-1

### The E-M Algorithm

while(diff>ptol && niter<maxiter)
{

niter<-niter+1

ptheta.old <- ptheta.new
f0.old <- f0.new
f1.old <- f1.new

## updating the weights and probabilities of hidden states

forwardbackward.res <- forwardbackward1.indep(x, ptheta.old, f0.old, f1.old)

gamma <- forwardbackward.res$pr
c0 <- forwardbackward.res$rescale

## updating the parameter estimates

ptheta.new <- apply(gamma,2,sum)/NUM

q5 <- sum(gamma[, 1]*x)
mu0 <- q5/sum(gamma[, 1])

q6 <- sum(gamma[, 1]*(x-mu0)*(x-mu0))
sd0 <- sqrt(q6/sum(gamma[, 1]))

f0.new<-c(mu0, sd0)

if(nulltype == 0){
	f0.new <- c(0,1)
}
if(nulltype == 1){
	f0.new <- c(locfdr_p0$fp0[3,1],locfdr_p0$fp0[3,2])
}
if(nulltype == 2){
	f0.new <- c(locfdr_p0$fp0[5,1],locfdr_p0$fp0[5,2])
}

q1 <- sum(gamma[, 2])
q2 <- sum(gamma[, 2]*x)
mu1 <- q2/q1
q3 <- sum(gamma[, 2]*(x-mu1)*(x-mu1))
sd1 <- sqrt(q3/q1)
f1.new <- c(mu1, sd1)

df2 <- abs(f1.old-f1.new)
diff <- max(df2)

if (is.na(diff)) {
   converged=FALSE;
   break;
}

}

 lfdr <- gamma[, 1]
if (converged) {
 logL<- sum(log(c0))
if (nulltype > 0) {
	BIC<-logL-(3*L+2)*log(NUM)/2 
} else {
	BIC<-logL-(3*L)*log(NUM)/2 
}

 em.var<-list(ptheta=ptheta.new, f0=f0.new, f1=f1.new, LIS=lfdr, logL=logL, BIC=BIC, 
 ni=niter, converged=converged) 
} else {
 BIC<- logL<- (-Inf)
 em.var<-list(ptheta=ptheta.old, f0=f0.old, f1=f1.old, LIS=lfdr, logL=logL, BIC=BIC, 
 ni=niter, converged=converged)
}

}

#######
# L>1
#######

else if (L>1)
{

ptheta.new <- c(0.5, 0.5)
pc.new <- rep(1, L)/L
mus <- seq(from=-1, by=1.5, length=L)
sds <- rep(1, L)
f0.new <- c(2, 1)

locfdr_p0 <- locfdr(x,plot=0)
if(nulltype == 0){
	f0.new <- c(0,1)
}
if(nulltype == 1){
	f0.new <- c(locfdr_p0$fp0[3,1],locfdr_p0$fp0[3,2])
}
if(nulltype == 2){
	f0.new <- c(locfdr_p0$fp0[5,1],locfdr_p0$fp0[5,2])
}


f1.new<-cbind(mus, sds)

diff<-1

### The E-M Algorithm

while(diff>ptol && niter<maxiter)
{

niter <- niter+1

ptheta.old <- ptheta.new
pc.old <- pc.new
f0.old <- f0.new
f1.old <- f1.new

## updating the weights and probabilities of hidden states

forwardbackward.res <- forwardbackward.indep(x, ptheta.old, pc.old, f0.old, f1.old)

gamma <- forwardbackward.res$pr
omega <- forwardbackward.res$wt
c0 <- forwardbackward.res$rescale

## updating the parameter estimates

ptheta.new <- apply(gamma,2,sum)/NUM

q5 <- sum(gamma[, 1]*x)
mu0 <- q5/sum(gamma[, 1])

q6 <- sum(gamma[, 1]*(x-mu0)*(x-mu0))
sd0 <- sqrt(q6/sum(gamma[, 1]))

f0.new <- c(mu0, sd0)

if(nulltype == 0){
	f0.new <- c(0,1)
}
if(nulltype == 1){
	f0.new <- c(locfdr_p0$fp0[3,1],locfdr_p0$fp0[3,2])
}
if(nulltype == 2){
	f0.new <- c(locfdr_p0$fp0[5,1],locfdr_p0$fp0[5,2])
}


mus <- 1:L
sds <- 1:L

for (c in 1:L)
{

  q1 <- sum(omega[, c])
  q2 <- sum(gamma[, 2])
  pc.new[c] <- q1/q2
  
  q3 <- sum(omega[, c]*x)
  mus[c] <- q3/q1

  q4 <- sum(omega[, c]*(x-mus[c])*(x-mus[c]))
  sds[c] <- sqrt(q4/q1)

}

f1.new <- cbind(mus, sds)

df2 <- abs(f1.old-f1.new)
diff <- max(df2)

if (is.na(diff)) {
   converged=FALSE;
   break;
}

}

lfdr <- gamma[, 1]
if (converged) {
 logL <- sum(log(c0))
if (nulltype > 0) {
	BIC <- logL-(3*L+2)*log(NUM)/2 
} else {
	BIC <- logL-(3*L)*log(NUM)/2 
}
 em.var<-list(ptheta=ptheta.new, pc=pc.new, f0=f0.new, f1=f1.new, LIS=lfdr, logL=logL, 
 BIC=BIC, ni=niter, converged=converged) 
} else {
 logL <- (-Inf)
 BIC <- logL<- (-Inf)
 em.var <-list(ptheta=ptheta.old, pc=pc.old, f0=f0.old, f1=f1.old, LIS=lfdr, logL=logL, 
 BIC=BIC, ni=niter, converged=converged)
}

}
}

return (em.var)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
