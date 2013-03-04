\name{em.nhmm}
\alias{em.nhmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intermediate function
}
\description{
Intermediate function
}
\usage{
em.nhmm(x, Z, dist, dist.included = TRUE, alttype = "mixnormal", L = 2, maxiter = 1000, 
nulltype = 2,symmetric = FALSE, epsilon = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
\item{dist}{
%%     ~~Describe \code{dist} here~~
}
\item{dist.included}{
%%     ~~Describe \code{dist.included} here~~
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
function(x, Z, dist, dist.included=TRUE, alttype='mixnormal', L=2, maxiter=1000, 
nulltype=2, symmetric=FALSE, epsilon=1e-4)
{

NUM<-length(x)

ptol<-epsilon

niter<-0

if(length(dist)>0&length(Z)>0){
  z1 <- length(dist)
  if(is.vector(Z)==TRUE){
    z2 <- length(Z)
  }else z2 <- dim(Z)[1]

  if(z1!=z2){
    cat('Error: x and Z not compatible','\n')
    stop
  }
}

if(length(Z)>0){
if(is.vector(Z)==TRUE) Z <- matrix(Z,ncol=1)
if(dim(Z)[1]!=NUM){
	cat('Error: x and Z not compatible','\n')
	stop
}
if(dim(Z)[2] == 1) Z <- scale(Z)
if(dim(Z)[2] > 1) Z <- apply(Z,2,scale)
cat('Scaling covariates','\n')
}

Z <- cbind(dist,Z)
Z_all <- matrix(0,ncol=dim(Z)[2],nrow=(NUM-1)*2)
for(i in 1:dim(Z)[2]){
	Z_all[,i] <- rep(Z[-1,i],each=2)
}

# Assuming it will converge 
converged=TRUE

if(alttype == 'kernel'){

tmp.trans.update <- 0
diff<-1
logL.iter <- 0

while(diff>ptol && niter<maxiter)
{

### control for error in update_trans.prob.nhmm ###

if(length(tmp.trans.update) == 1){

diff<-1
logL.iter <- 0
niter <- 0

trans.par1.new <- c(0,0,0,0)
trans.par2.new <- rnorm(3+dim(Z)[2])
trans.par2.new[4] <- abs(trans.par2.new[4])

tmp.trans.prob <- compute.A.nhmm (Z, trans.par1.new, trans.par2.new, 
dist.included=dist.included)
pii.new <- tmp.trans.prob$pii
A.new <- tmp.trans.prob$A

f0.new<-c(0, 1)
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
}

### The E-M Algorithm

niter<-niter+1

trans.par1.old <- trans.par1.new
trans.par2.old <- trans.par2.new
pii.old<-pii.new
A.old<-A.new
f0.old<-f0.new
f1.old<-f1.new

## updating the weights and probabilities of hidden states

forwardbackward.res <- forwardbackward1.kernel(x, pii.old, A.old, f0.old, f1.old)

gamma <- forwardbackward.res$pr
dgamma <- forwardbackward.res$ts
c0 <- forwardbackward.res$rescale

## updating the parameter estimates

tmp.trans.update <- try(update_trans.prob.nhmm(Z, dgamma, gamma, trans.par1.old, 
trans.par2.old,iter.conj.grad=10, dist.included=dist.included))

if(length(tmp.trans.update) > 2){

pii.new <- tmp.trans.update$pii
A.new <- tmp.trans.update$A
trans.par1.new <- tmp.trans.update$trans.par1
trans.par2.new <- tmp.trans.update$trans.par2

npar.A <- dim(Z)[2] + 2

q5<-sum(gamma[, 1]*x)
mu0<-q5/sum(gamma[, 1])

q6<-sum(gamma[, 1]*(x-mu0)*(x-mu0))
sd0<-sqrt(q6/sum(gamma[, 1]))

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

logL.iter <- c(logL.iter,-sum(log(c0)))

df1<-abs(trans.par2.old[-1]-trans.par2.new[-1])
df2<-abs(f1.old-f1.new)
df3 <- abs(logL.iter[niter+1] - logL.iter[niter])

diff<-max(df1, df2, df3)

}

if (is.na(diff)) {
   converged=FALSE;
   break;
}

}

lfdr<-gamma[, 1]
if (converged) {
 	logL<--sum(log(c0))
if (nulltype > 0) {
	BIC <- logL - (3 + npar.A + 2 - 2)*log(NUM)/2 
} else {
	BIC <- logL - (3 + npar.A - 2)*log(NUM)/2 
}

em.var<-list(pii=pii.new, A=A.new, f0=f0.new, f1=kern.f1, LIS=lfdr, logL=logL, BIC=BIC, 
ni=niter, trans.par2 = trans.par2.new, converged=converged,logL.iter=logL.iter[-1]) 
} else {
BIC<- logL<- (-Inf)
em.var<-list(pii=pii.old, A=A.old, f0=f0.old, f1=kern.f1, LIS=lfdr, logL=logL, BIC= BIC, 
ni=niter, trans.par2 = trans.par2.new, converged=converged,logL.iter=logL.iter[-1])
}

}


if(alttype == 'mixnormal'){

########
# L=1
########

if (L==1)
{

tmp.trans.update <- 0
diff<-1
logL.iter <- 0

while(diff>ptol && niter<maxiter)
{

### control for error in update_trans.prob.nhmm ###

if(length(tmp.trans.update) == 1){

diff<-1
logL.iter <- 0
niter <- 0

trans.par1.new <- c(0,0,0,0)
trans.par2.new <- rnorm(3+dim(Z)[2])
trans.par2.new[4] <- abs(trans.par2.new[4])

tmp.trans.prob <- compute.A.nhmm (Z, trans.par1.new, trans.par2.new,
dist.included=dist.included)
pii.new <- tmp.trans.prob$pii
A.new <- tmp.trans.prob$A

f0.new<-c(0, 1)
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
f1.new <- c(4, 1)
}

### The E-M Algorithm

niter<-niter+1

trans.par1.old <- trans.par1.new
trans.par2.old <- trans.par2.new
pii.old<-pii.new
A.old<-A.new
f0.old<-f0.new
f1.old<-f1.new

## updating the weights and probabilities of hidden states

forwardbackward.res <- forwardbackward1(x, pii.old, A.old, f0.old, f1.old)

gamma <- forwardbackward.res$pr
dgamma <- forwardbackward.res$ts
c0 <- forwardbackward.res$rescale

## updating the parameter estimates

tmp.trans.update <- try(update_trans.prob.nhmm(Z, dgamma, gamma, trans.par1.old, 
trans.par2.old,iter.conj.grad=10,dist.included=dist.included))

if(length(tmp.trans.update) > 2){

pii.new <- tmp.trans.update$pii
A.new <- tmp.trans.update$A
trans.par1.new <- tmp.trans.update$trans.par1
trans.par2.new <- tmp.trans.update$trans.par2

npar.A <- dim(Z)[2] + 2

q5<-sum(gamma[, 1]*x)
mu0<-q5/sum(gamma[, 1])

q6<-sum(gamma[, 1]*(x-mu0)*(x-mu0))
sd0<-sqrt(q6/sum(gamma[, 1]))

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

q1<-sum(gamma[, 2])
q2<-sum(gamma[, 2]*x)
mu1<-q2/q1
q3<-sum(gamma[, 2]*(x-mu1)*(x-mu1))
sd1<-sqrt(q3/q1)
f1.new<-c(mu1, sd1)

logL.iter <- c(logL.iter,-sum(log(c0)))

df1<-abs(trans.par2.old[-1]-trans.par2.new[-1])
df2<-abs(f1.old-f1.new)
df3 <- abs(logL.iter[niter+1] - logL.iter[niter])

diff<-max(df1, df2, df3)

}

if (is.na(diff)) {
   converged=FALSE;
   break;
}

}

lfdr<-gamma[, 1]
if (converged) {
 	logL<--sum(log(c0))
if (nulltype > 0) {
	BIC <- logL - (3*L + npar.A + 2)*log(NUM)/2 
} else {
	BIC <- logL - (3*L + npar.A)*log(NUM)/2 
}

em.var<-list(pii=pii.new, A=A.new, f0=f0.new, f1=f1.new, LIS=lfdr, logL=logL, BIC=BIC, 
ni=niter, trans.par2 = trans.par2.new, converged=converged,logL.iter=logL.iter[-1]) 
} else {
BIC<- logL<- (-Inf)
em.var<-list(pii=pii.old, A=A.old, f0=f0.old, f1=f1.old, LIS=lfdr, logL=logL, BIC= BIC, 
ni=niter, trans.par2 = trans.par2.new, converged=converged,logL.iter=logL.iter[-1])
}

}

#######
# L>1
#######

else if (L>1)
{

tmp.trans.update <- 0
diff<-1
logL.iter <- 0

while(diff>ptol && niter<maxiter)
{

### control for error in update_trans.prob.nhmm ###

if(length(tmp.trans.update) == 1){

diff<-1
logL.iter <- 0
niter <- 0

trans.par1.new <- c(0,0,0,0)
trans.par2.new <- rnorm(3+dim(Z)[2])
trans.par2.new[4] <- abs(trans.par2.new[4])

tmp.trans.prob <- compute.A.nhmm (Z, trans.par1.new, trans.par2.new,dist.included=dist.included)
pii.new <- tmp.trans.prob$pii
A.new <- tmp.trans.prob$A

pc.new<-rep(1, L)/L
mus<-seq(from=-1, by=1.5, length=L)
sds<-rep(1, L)
f0.new<-c(0, 1)
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

}

### The E-M Algorithm

niter<-niter+1

trans.par1.old <- trans.par1.new
trans.par2.old <- trans.par2.new
pii.old<-pii.new
A.old<-A.new
pc.old <- pc.new
f0.old<-f0.new
f1.old<-f1.new

## updating the weights and probabilities of hidden states

forwardbackward.res <- forwardbackward(x, pii.old, A.old, pc.old, f0.old, f1.old)

gamma <- forwardbackward.res$pr
dgamma <- forwardbackward.res$ts
omega <- forwardbackward.res$wt
c0 <- forwardbackward.res$rescale

## updating the parameter estimates

tmp.trans.update <- try(update_trans.prob.nhmm(Z, dgamma, gamma, trans.par1.old, 
trans.par2.old,iter.conj.grad=10,dist.included=dist.included))

if(length(tmp.trans.update) > 2){

pii.new <- tmp.trans.update$pii
A.new <- tmp.trans.update$A
trans.par1.new <- tmp.trans.update$trans.par1
trans.par2.new <- tmp.trans.update$trans.par2

npar.A <- dim(Z)[2] + 2

q5<-sum(gamma[, 1]*x)
mu0<-q5/sum(gamma[, 1])

q6<-sum(gamma[, 1]*(x-mu0)*(x-mu0))
sd0<-sqrt(q6/sum(gamma[, 1]))

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

mus<-1:L
sds<-1:L

for (c in 1:L)
{

 q1<-sum(omega[, c])
 q2<-sum(gamma[, 2])
 pc.new[c]<-q1/q2
  
 q3<-sum(omega[, c]*x)
 mus[c]<-q3/q1

 q4<-sum(omega[, c]*(x-mus[c])*(x-mus[c]))
 sds[c]<-sqrt(q4/q1)

}
f1.new<-cbind(mus, sds)

logL.iter <- c(logL.iter,-sum(log(c0)))


df1<-abs(trans.par2.old[-1]-trans.par2.new[-1])
df2<-abs(f1.old-f1.new)
df3 <- abs(logL.iter[niter+1] - logL.iter[niter])

diff<-max(df1, df2, df3)

}

if (is.na(diff)) {
   converged=FALSE;
   break;
}

}

 lfdr<-gamma[, 1]
if (converged) {
 logL <- -sum(log(c0))
if (nulltype>0) {
	BIC <- logL - (3*L + npar.A + 2)*log(NUM)/2 
} else {
	BIC <- logL - (3*L + npar.A)*log(NUM)/2 
}
em.var<-list(pii=pii.new, A=A.new, pc=pc.new, f0=f0.new, f1=f1.new, LIS=lfdr, logL=logL, 
BIC=BIC, ni=niter, trans.par2 = trans.par2.new, converged=converged,
logL.iter=logL.iter[-1])
} else {
 logL<- (-Inf)
 BIC<- logL<- (-Inf)
 em.var <- list(pii=pii.old, A=A.old, pc=pc.old, f0=f0.old, f1=f1.old, LIS=lfdr, logL=logL, 
 BIC=BIC, ni=niter, trans.par2 = trans.par2.new, converged=converged,
 logL.iter=logL.iter[-1])
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
