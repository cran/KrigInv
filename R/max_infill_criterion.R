max_infill_criterion <- function(lower, upper, parinit=NULL, sampling.method, method.param=0, control=NULL, 
                                 discrete.X=NULL, T, model, type="UK"){

	d <- model@d

	if (sampling.method=="tmse")        funk.optim <- tmse_optim
	else if (sampling.method=="ranjan") funk.optim <- ranjan_optim
	else if (sampling.method=="bichon") funk.optim <- bichon_optim
	else{
		funk.optim <- ranjan_optim
        print("Unknown sampling criterion - switched to Ranjan EI")
	}

	if (is.null(discrete.X)){
		#default OPTIMIZATION PARAMETERS when they are not provided
		if (is.null(control$pop.size))  control$pop.size <- 10#floor(4 + 3 * log(d))
		if (is.null(control$max.generations))  control$max.generations <- 10*d#100*d
		if (is.null(control$wait.generations))  control$wait.generations <- 2#10#2
		if (is.null(control$BFGSburnin)) control$BFGSburnin <- 0#10#0
		if (is.null(parinit))  parinit <- lower + runif(d) * (upper - lower)

		domaine <- cbind(lower, upper)

		o <- genoud(fn=funk.optim, nvars=d, max=TRUE, pop.size=control$pop.size,
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
				model=model, T=T, method.param=method.param, type=type)
	}
	else{
        x <- t(discrete.X)
		all.crit <- funk.optim(x=x, T=T, method.param=method.param, model, type="SK")
		ibest <- which.max(all.crit)
		o <- list(2)
		o$par <- t(discrete.X[ibest,])
		o$value <- max(all.crit)
		
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
