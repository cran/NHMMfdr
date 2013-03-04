LIS.adjust <-
function (lis, fdr = 0.001, adjust = TRUE) 
 {
     n = length(lis)
     s.lis = sort(lis)
     rank.lis = rank(lis)
 
     mean_s.lis <- cumsum(s.lis)/c(1:n)
     i <- min(which(mean_s.lis>fdr))
 
     nNonNull = i - 1
     States = rep(0, n)
     if (nNonNull > 0) 
         States[lis <= s.lis[nNonNull]] = 1
     if (adjust) {
         #aLIS = sapply(lis, function(cut) mean(lis[which(lis <= 
            # cut)]))
    		s.aLIS = cumsum(s.lis)/1:n
    		aLIS = s.aLIS[rank.lis]
         return(list(States = States, aLIS = aLIS))
     }
     else {
         return(list(States = States))
     }
 }

