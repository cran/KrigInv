bichon_optim <-
function(x, T, method.param, model, type="UK"){

	krig  <- predict.km(model, newdata=as.data.frame(t(x)), type)
	mk    <- krig$mean
	sk    <- krig$sd
	alpha <- method.param

	t <- (mk-T)/sk
	tplus <- t + alpha
	tminus <- t - alpha

	G <- alpha*(pnorm(tplus)-pnorm(tminus))
		- t*(2*pnorm(t) - pnorm(tminus) - pnorm(tminus))
		- (2*dnorm(t) - dnorm(tplus) - dnorm(tminus))

	crit <- G*sk
	return(crit)
}

