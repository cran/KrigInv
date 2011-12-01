sur_optim_parallel <- function(x, integration.points,integration.weights=NULL,
						intpoints.oldmean,intpoints.oldsd,precalc.data,
						model, T, new.noise.var=NULL,batchsize,current.sur){
	
	if(!is.null(new.noise.var)){
		if(new.noise.var == 0) {
			new.noise.var <- NULL
		}
	}
	#x is a vector of size d * batchsize
	d <- model@d
	n <- model@n
	X.new <- matrix(x,nrow=d)
	mindist <- Inf
		
	tp1 <- c(as.numeric(t(model@X)),x)
	for (i in 1:batchsize){
		#distance between the i^th point and all other points (in the DOE or in the batch)
		xx <- X.new[,i]
		tp2<-matrix(tp1-as.numeric(xx),ncol=d,byrow=TRUE)^2
		mysums <- sqrt(rowSums(tp2))
		mysums[n+i] <- Inf #because this one is usually equal to zero...
		mindist <- min(mindist,mysums)		
	}
		
	if (!identical(colnames(integration.points), colnames(model@X))) colnames(integration.points) <- colnames(model@X)
	
	if ((mindist > 1e-5) || (!is.null(new.noise.var))){
		X.new <- t(X.new)	
		krig  <- predict_nobias_km(object=model, newdata=as.data.frame(X.new), 
								type="UK",se.compute=TRUE, cov.compute=TRUE) 
	
		mk <- krig$mean ; sk <- krig$sd ; newXvar <- sk*sk
		F.newdata <- krig$F.newdata ; c.newdata <- krig$c;Sigma.r <- krig$cov
		Sigma.r.norm <- krig$cov / crossprod(t(sk))
	
		kn = computeQuickKrigcov(model,integration.points,X.new,precalc.data, F.newdata , c.newdata) 
	
		krig2  <- predict_update_km_parallel (newXmean=mk,newXvar=newXvar,newXvalue=mk, 
				Sigma.r=Sigma.r,newdata.oldmean=intpoints.oldmean,newdata.oldsd=intpoints.oldsd,kn=kn)
		if(!is.null(krig2$error)) return(current.sur)
		
		sk.new <- krig2$sd;lambda <- krig2$lambda	
		a <- (intpoints.oldmean-T) / sk.new
	
		tmp <- crossprod(t(1/sk.new),t(sk))
		b <- lambda * tmp
	
		a[a==Inf]<- 1000 ;a[a== -Inf] <- 1000;b[b==Inf]<-0;b[b== -Inf]<-0
		a[is.nan(a)] <- 1000; b[is.nan(b)]<-0
	
		bt.sigma.r.b <- as.numeric(rowSums(b*crossprod(t(b),Sigma.r.norm)))
		a.new <- as.numeric(a/sqrt(1+bt.sigma.r.b))
		b.new <- as.numeric(-bt.sigma.r.b/(1+bt.sigma.r.b))
		Phi.biv.a.b <- pbivnorm(a.new,-a.new,b.new)  #c.d.f of the bivariate gaussian distribution
	
		if (is.null(integration.weights)) {crit <- mean(Phi.biv.a.b)
		}else crit <- sum(Phi.biv.a.b*integration.weights)
	}else crit <- current.sur + 0.01	
	
	return(crit)
}

