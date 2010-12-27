sur_optim <- function(x, y.integration=NULL, integration.points, model, T, new.noise.var=0, type="UK"){

	mindist <- dist(rbind(model@X[1,],t(x)))
	
	if (is.null(y.integration)){
		n.int.y <- 10
		y.temp <- seq(0,1-1/n.int.y,1/n.int.y)
		y.temp <- .5/(n.int.y) + y.temp
		y.integration <- qnorm(y.temp)
    }
	
	integration.points <- as.matrix(integration.points)
	if (!identical(colnames(integration.points), colnames(model@X))) colnames(integration.points) <- colnames(model@X)

	
	for (i in 1:model@n)  mindist <- min(mindist, dist(rbind(model@X[i,],t(x))))

	if (mindist > 1e-4 || new.noise.var > 0){
		# Compute simulated values of y.new
		krig  <- predict.km(model, newdata=as.data.frame(t(x)), type)
		mk <- krig$mean
		sk <- krig$sd
		y.new <- mk + sqrt(sk^2+new.noise.var)*y.integration
		n.y <- length(y.new)
		
		# Update kriging with new point
		X.new <- t(x)
		ynew <- y.new[1]
		model <- update.km(model=model,NewX=X.new,NewY=ynew,CovReEstimate=FALSE,new.noise.var=new.noise.var)
		
		p.all <- seq(1,n.y)
		
		for (i in 1:n.y){
			#modify the value of the last point
			ynew <- y.new[i]
			
			#solution 1 (a bit slower)
			#model <- update.km(model=model,NewX=model@X[model@n],NewY=ynew,NewX_AllreadyExist=TRUE,CovReEstimate=FALSE,new.noise.var=new.noise.var)
			
			#solution 2 (a bit faster because these are no full copy of the km object)
			model@y[model@n] <- ynew
			model@z <- as.numeric(backsolve(t(model@T), model@y-model@F%*%as.matrix(model@trend.coef), upper.tri=FALSE))

			#calculation of mk and sk, the kriging mean an standard-dev
			if (i==1){
				#this call to predict.km is expensive. However we do it only once
				krig  <- predict.km(model, newdata=as.data.frame(integration.points), type)
				mk    <- krig$mean	
				sk    <- krig$sd	#calculated only once
				
				Tinv.c.newdata=krig$Tinv.c	#save this matrix to accelerate the calculation of mk for the next values of i
				F.newdata <- model.matrix(model@trend.formula, data=data.frame(integration.points))
				y.predict.trend <- F.newdata%*%model@trend.coef		#save this array to accelerate the calculation of mk for the next values of i
			}
			else{
				y.predict.complement <- t(Tinv.c.newdata)%*%model@z			#model@z has been updated in the update.km function
				mk <- as.numeric(y.predict.trend + y.predict.complement)	#quick calculation of the kriging mean. And the kriging sd is unchanged !!!
			}
			
			p.all[i] = mean(1 - pnorm(abs(T-mk), 0, sk))
		}
		crit <- mean(p.all)
	}
	else{
		krig  <- predict.km(model, newdata=as.data.frame(integration.points), type)
		mk    <- krig$mean
		sk    <- krig$sd
		p = 1 - pnorm(abs(T-mk), 0, sk)
		crit <- mean(p)
	}
	
	return(crit)
}













