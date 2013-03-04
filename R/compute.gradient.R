compute.gradient <-
function(Z, dgamma, gamma, trans.par1, trans.par2, dist.included=TRUE){

	gradient <- rep(0,3+dim(Z)[2])

	tmp.trans.prob <- compute.A.nhmm(Z, trans.par1, trans.par2, dist.included=dist.included)

	pii <- tmp.trans.prob$pii
	A <- tmp.trans.prob$A
	
	gradient[1] <- gamma[1,2] - pii[2]
	gradient[2] <- sum(dgamma[1,2,] - gamma[-dim(gamma)[1],1]*A[1,2,])
	gradient[3] <- sum(dgamma[2,2,] - gamma[-dim(gamma)[1],2]*A[2,2,])
		
	tmp <- rep(0,dim(dgamma)[3])
        
	for(i in 1:2){
		tmp <- tmp + gamma[-dim(Z)[1],i]*A[i,2,] 
	}
		gradient[-c(1:3)] <- (gamma[1,2] - pii[2])*Z[1,] + apply(matrix((gamma[-1,2] - tmp)*Z[-1,],ncol=dim(Z)[2]),2,sum)

        if(dist.included==TRUE){
          gradient[4] <- (gamma[1,2] - pii[2])*Z[1,1] + sum(dgamma[1,2,]*Z[-1,1]) - sum(dgamma[2,2,]*Z[-1,1]) - sum(gamma[-dim(Z)[1],1]*A[1,2,]*Z[-1,1]) + sum(gamma[-dim(Z)[1],2]*A[2,2,]*Z[-1,1])
        }
        
	return(-gradient)

}

