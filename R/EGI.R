EGI <-function(T, model, method.param=NULL, sampling.method=NULL, discrete.X=NULL,
               new.noise.var=0, fun, iter, lower, upper, integration.points=NULL,
               parinit=NULL, control=NULL, kmcontrol=NULL) {

	n <- nrow(model@X)

	if (is.null(kmcontrol$penalty)) kmcontrol$penalty <- model@penalty
	if (length(model@penalty==0)) kmcontrol$penalty <- NULL 

	if (is.null(kmcontrol$optim.method)) kmcontrol$optim.method <- model@optim.method 
	if (is.null(kmcontrol$parinit)) kmcontrol$parinit <- model@parinit
	if (is.null(kmcontrol$control)) kmcontrol$control <- model@control
	
	if (is.null(sampling.method)) sampling.method <- "ranjan"

	for (i in 1:iter) {

        if (sampling.method=="timse"){
			oEGI <- max_timse(lower=lower, upper=upper,
							parinit=parinit, control=control,
							T=T,epsilon=method.param,discrete.X=discrete.X,integration.points=integration.points,
							model=model, new.noise.var=new.noise.var, type="UK")
		}
		else if (sampling.method=="sur"){
			oEGI <- max_sur(lower=lower, upper=upper,
							parinit=parinit, control=control,
							T=T,model=model,discrete.X=discrete.X,integration.points=integration.points,
							new.noise.var=new.noise.var, type="UK")
		}
        else{
			oEGI <- max_infill_criterion(lower=lower, upper=upper,
										parinit=parinit, control=control,
										T=T, sampling.method=sampling.method, method.param=method.param,discrete.X=discrete.X,
										model=model,type="UK")
		}
		
		print("New point")
        print(oEGI$par)
		
		#call to update.km
		X.new <- oEGI$par
		y.new <- fun(oEGI$par) + rnorm(1, 0, new.noise.var)
		model <- update.km(model=model,NewX=X.new,NewY=y.new,CovReEstimate=TRUE,new.noise.var=new.noise.var,kmcontrol=kmcontrol)
	}
	
	return(list(
				par=model@X[(n+1):(n+iter),, drop=FALSE], 
				value=model@y[(n+1):(n+iter),, drop=FALSE], 
				npoints=1, 
				nsteps=iter, 
				lastmodel=model
				)
			)

}
