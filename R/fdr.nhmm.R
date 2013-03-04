fdr.nhmm <-
function(x, Z = NULL, dist = NULL, log.transform.dist = TRUE, alttype = 'kernel', L=2, maxiter=1000, nulltype=0, modeltype = 'NHMM', symmetric=FALSE, epsilon=1e-4)
{
	cat('NOTE: symmetric option is only available for alttype kernel','\n')
	if(modeltype!='NHMM'&modeltype!='HMM'&modeltype!='Indep'){
		cat('Error: modeltype must be NHMM, HMM or Indep','\n')
		return(0)
	}

	if(alttype!='kernel'&alttype!='mixnormal'){
		cat('Error: alttype must be kernel or mixnormal','\n')
		return(0)
	}

	if(alttype=='mixnormal'&is.null(L)==T){
		cat('Warning: alttype is mixnorm but L missing','\n')
		cat('Running with alttype kernel','\n')
		alttype = 'kernel'
		#return(0)
	}
	if(alttype=='kernel'){
		L = 1
	}
	if(nulltype >3){
		cat('Error: nulltype must be 0, 1, 2 or 3','\n')
		return(0)
	}

	fdr_res <- 0
	n.try <- 1
        dist.included=FALSE
        if(length(dist)>0) {
          dist.included=TRUE
          if(length(which(dist<0))>0){
            cat('Error: Distance must be positive!','\n')
            return(0)
          }
          if(log.transform.dist == TRUE) {
            cat('Transforming distance','\n')
            dist <- log2(dist+2)
          }
        }
	if(modeltype == 'NHMM' & length(Z)==0 &length(dist)==0){
		cat('Warning: Z missing...run as HMM','\n')

		while(length(fdr_res)==1&n.try<3){
		fdr_res <- try(em.hmm(x=x, alttype=alttype, L=L, maxiter=maxiter, nulltype=nulltype, symmetric=symmetric, epsilon=epsilon))
		if(length(fdr_res)==1)cat('Numerical ERROR: Rerunning...','\n')
		n.try <- n.try + 1
		}
	}
	if(modeltype == 'HMM'){
		cat('Running with alltype ',alttype,', nulltype', nulltype,', modeltype ',modeltype,'...','\n')

		while(length(fdr_res)==1&n.try<3){
		fdr_res <- try(em.hmm(x=x, alttype=alttype, L=L, maxiter=maxiter, nulltype=nulltype, symmetric=symmetric, epsilon=epsilon))
		if(length(fdr_res)==1)cat('Numerical ERROR: Rerunning...','\n')
		n.try <- n.try + 1
		}

	}
	if(modeltype == 'NHMM' & (length(Z) >0|length(dist)>0)){
		cat('Running with alltype ',alttype,', nulltype', nulltype,', modeltype ',modeltype,'...','\n')

		while(length(fdr_res)==1&n.try<3){
		fdr_res <- try(em.nhmm(x=x, Z=Z, dist, dist.included=dist.included, alttype=alttype, L=L, maxiter=maxiter, nulltype=nulltype, symmetric=symmetric, epsilon=epsilon))
		if(length(fdr_res)==1)cat('Numerical ERROR: Rerunning...','\n')
		n.try <- n.try + 1
		}

	}
	if(modeltype == 'Indep'){
		cat('Running with alltype ',alttype,', nulltype', nulltype,', modeltype ',modeltype,'...','\n')
		
		while(length(fdr_res)==1&n.try<3){
		fdr_res <- try(em.indep(x=x, alttype=alttype, L=L, maxiter=maxiter, nulltype=nulltype, symmetric=symmetric, epsilon=epsilon))
		if(length(fdr_res)==1)cat('Numerical ERROR: Rerunning...','\n')
		n.try <- n.try + 1
		}

	}
	if(length(fdr_res)==1)cat('ERROR: Try with different parameter settings','\n')
	if(length(fdr_res)>1)cat('DONE!','\n')
	return(fdr_res)
}

