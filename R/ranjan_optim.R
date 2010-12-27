ranjan_optim <-
function(x, T, method.param, model, type="UK"){

	krig  <- predict.km(model, newdata=as.data.frame(t(x)), type)
	mk    <- krig$mean
	sk    <- krig$sd
	alpha <- method.param

	p <- 1 - pnorm(T,mk,sk)

	n.points <- length(mk)
	G <- mk

	for (i in 1:n.points){
		if (sk[i]>0){
			t <- (mk[i]-T)/sk[i]
			tplus <- t + alpha		
			tminus <- t - alpha
			G[i] <- (alpha^2 - 1 - t^2)*(pnorm(tplus)-pnorm(tminus))- 2*t*(dnorm(tplus) - dnorm(tminus))+ tplus*dnorm(tplus) - tminus*dnorm(tminus)
		}
		else G[i] <- 0
	}
	
	crit <- G*(sk^2)
	return(crit)
}

