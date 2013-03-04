find.nu <-
function(Z, dgamma, gamma, delta1, delta2, trans.par1, trans.par2, dist.included=TRUE){

	p <- dim(Z)[2]
	# delta1 = (lambda1_H, sigma11_H, sigma21_H, rho1_H)
	delta_H <- array(0,c(2,2,c(dim(Z)[1]-1)))
	delta1_H0 <- delta1[1] + sum(delta1[-c(1:3)]*Z[1,])
	delta2_H0 <- delta2[1] + sum(delta2[-c(1:3)]*Z[1,])
    
	tmp11 <- t(delta1[-c(1:3)]*t(Z[-1,]))
	tmp21 <- t(delta1[-c(1:3)]*t(Z[-1,]))
	tmp12 <- t(delta2[-c(1:3)]*t(Z[-1,]))
	tmp22 <- t(delta2[-c(1:3)]*t(Z[-1,]))

        if(dist.included==TRUE){
          delta1_H0 <- delta1[1] + sum(delta1[-c(1:4)]*Z[1,-1])
          delta2_H0 <- delta2[1] + sum(delta2[-c(1:4)]*Z[1,-1])
    
          #tmp11 <- t(c(delta1[4],delta1[-c(1:4)])*t(Z[-1,]))
          tmp21 <- t(c(-delta1[4],delta1[-c(1:4)])*t(Z[-1,]))
          #tmp12 <- t(c(delta2[4],delta2[-c(1:4)])*t(Z[-1,]))
          tmp22 <- t(c(-delta2[4],delta2[-c(1:4)])*t(Z[-1,]))
        }
        
	denom11 = denom21 = denom12 = denom22 <- rep(0,dim(Z)[1]-1)

	for(i in 1:p){
		denom11 <- denom11 + tmp11[,i]
		denom21 <- denom21 + tmp21[,i]
		denom12 <- denom12 + tmp12[,i]
		denom22 <- denom22 + tmp22[,i]
	}
	
	delta_H[1,1,] <- delta1[2] + denom11
	delta_H[2,1,] <- delta1[3] + denom21
	delta_H[1,2,] <- delta2[2] + denom12
	delta_H[2,2,] <- delta2[3] + denom22	

	
	nu.new <- 0
	ptol <- 1
	iter <- 0

while(ptol > 1e-6){	
	iter <- iter + 1
	
	trans.par1.new <- trans.par1 + nu.new*delta1
	trans.par2.new <- trans.par2 + nu.new*delta2

	tmp.trans.prob <- compute.A.nhmm (Z, trans.par1.new, trans.par2.new,dist.included=dist.included)
	pii <- tmp.trans.prob$pii
	A <- tmp.trans.prob$A
				
	f <- sum(c(delta1_H0,delta2_H0)*(gamma[1,] - pii))
	for(i in 1:2){
		for(j in 1:2){	
			f <- f + sum(delta_H[i,j,]*(dgamma[i,j,] - gamma[-dim(gamma)[1],i]*A[i,j,]))
		}
	}

     fprime <- -sum(c(delta1_H0,delta2_H0)^2*pii*(1 - pii))
	for(i in 1:2){
		for(j in 1:2){
			fprime <- fprime - sum(delta_H[i,j,]^2*gamma[-dim(gamma)[1],i]*A[i,j,]*(1 - A[i,j,]))
		}
	}
	nu.old <- nu.new
	nu.new <- nu.old - f/fprime
	ptol <- abs(f)
	if(is.na(ptol)||iter > 100){
		nu.new <- NaN
		break
	}
}
	trans.par1.new <- trans.par1 + nu.new*delta1
	trans.par2.new <- trans.par2 + nu.new*delta2
	
	return(list(nu=nu.new,trans.par1=trans.par1.new,trans.par2= trans.par2.new))
}

