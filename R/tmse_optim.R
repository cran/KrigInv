tmse_optim <-
function(x, T, method.param=0, model, type="UK"){

	krig  <- predict.km(model, newdata=as.data.frame(t(x)), type)
	mk    <- krig$mean
	sk    <- krig$sd
	epsilon <- method.param

	W <- 1/sqrt(2*pi*(sk^2+epsilon^2))*exp(-0.5*((mk-T)/sqrt(sk^2+epsilon^2))^2)

	tmse <- W*sk^2
	return(tmse)
}