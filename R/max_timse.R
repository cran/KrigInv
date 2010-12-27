max_timse <- function(lower, upper, parinit=NULL, control=NULL, discrete.X=NULL, integration.points=NULL,
                      T, epsilon=0, model, new.noise.var=0, type="UK"){
	
	if(is.null(epsilon)) epsilon <- 0
	
	if (is.null(discrete.X)){
		d <- model@d
		
		# INTEGRATION POINTS
		# case 0: integration points defined by user
		
		if (is.null(integration.points)){
			# case 1: integration points generated using default values
			n.int.points <- 100*d
			integration.points <- maximinLHS(n.int.points, d, 100)
		}
		
		else if (length(integration.points) == 1){
			# case 2: number of integration points defined by user
			n.int.points <- integration.points
			integration.points <- maximinLHS(n.int.points, d, 100)
		}
		
		# default OPTIMIZATION PARAMETERS when they are not provided
		if (is.null(control$pop.size))  control$pop.size <- 10
		if (is.null(control$max.generations))  control$max.generations <- 10*d#100*d
		if (is.null(control$wait.generations))  control$wait.generations <- 2
		if (is.null(control$BFGSburnin)) control$BFGSburnin <- 0#10#0
		if (is.null(parinit))  parinit <- lower + runif(d) * (upper - lower)
		
		domaine <- cbind(lower, upper)
		
		# COMPUTE WEIGHT BEFORE OPTIMIZATION
		krig  <- predict.km(model, newdata=as.data.frame(integration.points), type)
		mk <- krig$mean
		sk    <- krig$sd
		weight <- 1/sqrt(2*pi*(sk^2+epsilon^2))*exp(-0.5*((mk-T)/sqrt(sk^2+epsilon^2))^2)
		weight[is.nan(weight)] <- 0

		
		# OPTIMIZATION
		o <- genoud(fn=timse_optim, nvars=d, max=FALSE, pop.size=control$pop.size,
				max.generations=control$max.generations,
				wait.generations=control$wait.generations,
				hard.generation.limit=TRUE, starting.values=parinit, MemoryMatrix=TRUE,
				Domains=domaine, default.domains=10, solution.tolerance=0.01,
				#gr=gr,
				boundary.enforcement=2, lexical=FALSE, gradient.check=FALSE, BFGS=TRUE,
				data.type.int=FALSE, hessian=TRUE, unif.seed=812821, int.seed=53058,
				print.level=2, share.type=0, instance.number=0,
				output.path="stdout", output.append=FALSE, project.path=NULL,
				P1=50, P2=50, P3=50, P4=50, P5=50, P6=50, P7=50, P8=50, P9=0,
				P9mix=NULL, BFGSburnin=control$BFGSburnin,BFGSfn=NULL, BFGShelp=NULL,
				cluster=FALSE, balance=FALSE, debug=FALSE,
				model=model, weight=weight, integration.points=integration.points, 
				new.noise.var=new.noise.var, type="SK")
	}

	else{
		# COMPUTE WEIGHT BEFORE OPTIMIZATION
		krig  <- predict.km(model, newdata=as.data.frame(discrete.X), type)
		mk <- krig$mean
		sk    <- krig$sd
		
		weight <- 1/sqrt(2*pi*(sk^2+epsilon^2))*exp(-0.5*((mk-T)/sqrt(sk^2+epsilon^2))^2)
		weight[is.nan(weight)] <- 0

		all.crit <- seq(1,nrow(discrete.X))
		for (i in 1:nrow(discrete.X)){
			all.crit[i] <- timse_optim(x=t(discrete.X[i,]), integration.points=discrete.X,
			weight=weight, model=model, new.noise.var=new.noise.var, type=type)
		}
		
		ibest <- which.min(all.crit)
		o <- list(2)
		o$par <- t(discrete.X[ibest,])
		o$value <- min(all.crit)
		
		#scatterplot3d(discrete.X[,1],discrete.X[,2],all.crit)
		#print("ibest")
		#print(ibest)
		#print(o$value)
	}

	o$par
	o$par <- t(as.matrix(o$par))
	colnames(o$par) <- colnames(model@X)
	o$value <- as.matrix(o$value)
	colnames(o$value) <- colnames(model@y)   
	return(list(par=o$par, value=o$value)) 
}