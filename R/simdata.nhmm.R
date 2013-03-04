simdata.nhmm <-
function(NUM, pii, A, f0, pc, f1)
{


theta<-rep(0, NUM)
x<-rep(0, NUM)
nc<-length(pc)

if (nc==1)
{ 
  data <- simdata1.nhmm (NUM, pii, A, f0, f1)
}    

else
{
## generating the states
 # initial state
theta[1]<-rbinom(1, 1, pii[2])
 # other states
for (i in 2:NUM)
{
  if (theta[i-1]==0)
     theta[i]<-rbinom(1, 1, A[1, 2, i-1])
  else
     theta[i]<-rbinom(1, 1, A[2, 2, i-1])
}

## generating the observations
for (i in 1:NUM)
{
  if (theta[i]==0)
  {
    mu0<-f0[1]
    sd0<-f0[2]
    x[i]<-rnorm(1, mean=mu0, sd=sd0)
  }
  else
  { 
    c<-sample(1:nc, 1, prob=pc)
    mu1<-f1[c, 1]
    sd1<-f1[c, 2]
    x[i]<-rnorm(1, mean=mu1, sd=sd1)
  }
}
data<-list(s=theta, o=x)
}

return (data)

}

