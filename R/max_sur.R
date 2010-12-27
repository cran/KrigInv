max_sur <- function(lower, upper, parinit=NULL, control=NULL, discrete.X=NULL, n.int.y=10, integration.points=NULL,
                    T, model, new.noise.var=0, type="UK"){

	# COMPUTE Y VALUES BEFORE OPTIMIZATION
	y.temp <- seq(0,1-1/n.int.y,1/n.int.y)
	y.temp <- .5/(n.int.y) + y.temp
	y.integration <- qnorm(y.temp)

	d <- model@d

	if (is.null(discrete.X)){
		# INTEGRATION POINTS
		
		# case 0: integration points defined by user (nothing to do then)
		
		if (is.null(integration.points)){
			# case 1: integration points generated using default values
			#n.points.by.vertice <- round((100*d)^(1/d))
			n.int.points <- 100*d
			integration.points <- maximinLHS(n.int.points, d, 100)
		}
		
		else if (length(integration.points) == 1){
			# case 2: number of integration points defined by user
			n.int.points <- integration.points
			integration.points <- maximinLHS(n.int.points, d, 100)
		}
	
		# default OPTIMIZATION PARAMETERS when they are not provided
		if (is.null(control$pop.size))  control$pop.size <- 10#floor(4 + 3 * log(d))
		if (is.null(control$max.generations))  control$max.generations <- 10*d#100*d
		if (is.null(control$wait.generations))  control$wait.generations <- 10#2
		if (is.null(control$BFGSburnin)) control$BFGSburnin <- 0#10#0
		if (is.null(parinit))  parinit <- lower + runif(d) * (upper - lower)
	
		domaine <- cbind(lower, upper)
		
		# OPTIMIZATION
		o <- genoud(fn=sur_optim, nvars=d, max=FALSE, pop.size=control$pop.size,
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
			model=model, T=T, y.integration=y.integration, integration.points=integration.points, 
	        new.noise.var=new.noise.var, type="SK")
	}

	else{
		all.crit <- seq(1,nrow(discrete.X))
		
		for (i in 1:nrow(discrete.X)){
			all.crit[i] <- sur_optim(x=t(discrete.X[i,]), integration.points=discrete.X,
			T=T, y.integration=y.integration, model=model, new.noise.var=new.noise.var, type="SK")
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

	o$par <- t(as.matrix(o$par))
	colnames(o$par) <- colnames(model@X)
	o$value <- as.matrix(o$value)
	colnames(o$value) <- colnames(model@y)   
	return(list(par=o$par, value=o$value)) 
}