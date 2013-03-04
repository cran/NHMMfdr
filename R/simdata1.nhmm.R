simdata1.nhmm <-
function(NUM, pii, A, f0, f1)
{


theta<-rep(0, NUM)
x<-rep(0, NUM)

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
    x[i]<-rnorm(1, mean=f0[1], sd=f0[2])
  }
  else
  { 
    x[i]<-rnorm(1, mean=f1[1], sd=f1[2])
  }
}
data<-list(s=theta, o=x)
return (data)

}

