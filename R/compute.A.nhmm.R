compute.A.nhmm <-
function(Z, trans.par1, trans.par2, dist.included=TRUE){

	# Z is a matrix of size m x p
	p <- dim(Z)[2]

	pii <- rep(0,2)
	A <- array(0,c(2,2,dim(Z)[1]-1))

	pii[1] <- sum(exp(trans.par1[1] + sum(trans.par1[-c(1:3)]*Z[1,])))/(sum(exp(trans.par1[1] + sum(trans.par1[-c(1:3)]*Z[1,]))) + sum(exp(trans.par2[1] + sum(trans.par2[-c(1:3)]*Z[1,]))))
	pii[2] = 1 - pii[1]
	
	tmp12 <- t(trans.par2[-c(1:3)]*t(Z[-1,]))
	tmp22 <- t(trans.par2[-c(1:3)]*t(Z[-1,]))

        if(dist.included==TRUE){
          #print('#####YES########')
        
          tmp12 <- t(c(trans.par2[4],trans.par2[-c(1:4)])*t(Z[-1,]))
          tmp22 <- t(c(-trans.par2[4],trans.par2[-c(1:4)])*t(Z[-1,]))
        }
        
	denom12 = denom22 <- rep(0,dim(Z)[1]-1)

	for(i in 1:p){
		denom12 <- denom12 + tmp12[,i]
		denom22 <- denom22 + tmp22[,i]
	}
	
	num11 <- exp(trans.par1[2])
	num21 <- exp(trans.par1[3])
	num12 <- exp(trans.par2[2] + denom12)
	num22 <- exp(trans.par2[3] + denom22)
	
	
	A[1,1,] <- num11/(num11 + num12)
	A[1,2,] <- num12/(num11 + num12)
	A[2,1,] <- num21/(num21 + num22)
	A[2,2,] <- num22/(num21 + num22)

	return(list(A = A, pii = pii))
}

