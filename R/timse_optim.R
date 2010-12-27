timse_optim <- function(x, weight=1, integration.points, model, new.noise.var=0, type="UK"){

	mindist <- dist(rbind(model@X[1,],t(x)))
	for (i in 1:model@n)  mindist <- min(mindist, dist(rbind(model@X[i,],t(x))))

	if (mindist > 1e-4 || new.noise.var > 0){
		# Update kriging with new point
		X.new <- t(x)
		y.new <- mean(model@y)		#the algo does not depend of this value but we need to provide one for updating the kriging model
		model <- update.km(model=model,NewX=X.new,NewY=y.new,CovReEstimate=FALSE,new.noise.var=new.noise.var)
	}

	krig  <- predict.km(model, newdata=as.data.frame(integration.points), type)
	sk    <- krig$sd
	timse <- mean(weight*sk^2)

	return(timse)
}
