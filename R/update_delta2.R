update_delta2 <-
function(Z, dgamma, gamma, delta.old, gradient.new, gradient.old){

	num.fac <- sum((gradient.new - gradient.old)*gradient.new)
	denom.fac <- sum(gradient.old^2)
	fac <- num.fac/denom.fac
	if(fac<0) fac <- 0
	delta2 <- -gradient.new + fac*delta.old
	return(delta2)
}

