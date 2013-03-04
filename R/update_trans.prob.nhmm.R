update_trans.prob.nhmm <-
function(Z, dgamma, gamma, trans.par1, trans.par2,iter.conj.grad=10, dist.included=TRUE){

	trans.par1.old <- trans.par1
	trans.par2.old <- trans.par2

	gradient.old <- compute.gradient(Z, dgamma, gamma, trans.par1.old, trans.par2.old,dist.included=dist.included)
	delta2.old <- -gradient.old

	n.conj.grad = 0
	while(n.conj.grad <= iter.conj.grad){

	tmp <- try(find.nu(Z, dgamma, gamma, c(0,0,0,0),delta2.old,trans.par1.old, trans.par2.old,dist.included=dist.included))

	if(is.na(tmp$nu)){
		break		
	}
	trans.par1.new <- tmp$trans.par1
	trans.par2.new <- tmp$trans.par2

	gradient.new <- compute.gradient(Z, dgamma, gamma,trans.par1.new, trans.par2.new,dist.included=dist.included)
	delta2.new <- update_delta2(Z, dgamma, gamma, delta2.old, gradient.new, gradient.old)

	trans.par1.old <- trans.par1.new
	trans.par2.old <- trans.par2.new
        if(dist.included==TRUE&trans.par2.new[4]<0)trans.par2.new[4] <- 0.5
	gradient.old <- gradient.new
	delta2.old <- delta2.new
	n.conj.grad <- n.conj.grad + 1
	
	}

	if(!is.na(tmp$nu)){
		
	tmp.trans.prob <- compute.A.nhmm(Z, trans.par1.new, trans.par2.new,dist.included=dist.included)
	pii.new <- tmp.trans.prob$pii
	A.new <- tmp.trans.prob$A
	return(list(A = A.new, pii = pii.new, trans.par1 = trans.par1.new, trans.par2 = trans.par2.new))
	}else return(1)
}

