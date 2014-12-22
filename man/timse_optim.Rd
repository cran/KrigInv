\name{timse_optim}
\Rdversion{1.3}
\alias{timse_optim}

\title{ Targeted IMSE criterion }

\description{
Evaluation of the \code{"timse"} criterion for a candidate point. To be used in optimization routines, like in \code{\link{max_timse}}. 
To avoid numerical instabilities, the new point is evaluated only if it is not too close to an existing observation, or if there is some observation noise. 
The criterion is the integral of the posterior timse uncertainty.
}
\usage{
timse_optim(x, integration.points, integration.weights = NULL, 
intpoints.oldmean = NULL, intpoints.oldsd = NULL, 
precalc.data, model, T, new.noise.var = 0, weight = NULL)
}

\arguments{
  	\item{x}{Input vector of size d at which one wants to evaluate the criterion.}
  	\item{integration.points}{p*d matrix of points for numerical integration in the X space.}
	\item{integration.weights}{Vector of size p corresponding to the weights of these integration points.}
	\item{intpoints.oldmean}{Vector of size p corresponding to the kriging mean at the integration points before adding \code{x} to the design of experiments.}
	\item{intpoints.oldsd}{Vector of size p corresponding to the kriging standard deviation at the integration points before adding \code{x} to the design of experiments.}
	\item{precalc.data}{List containing useful data to compute quickly the updated kriging variance. This list can be generated using the \code{\link{precomputeUpdateData}} function.}
 	\item{model}{Object of class \code{\link[DiceKriging]{km}} (Kriging model).}
	\item{T}{Target value (scalar)}
	\item{new.noise.var}{Optional scalar value of the noise variance for the new observations.}
	\item{weight}{
Vector of weight function (length must be equal to the number of lines of the matrix integration.points). 
If nothing is set, the imse criterion is used instead if timse. It corresponds to equal weights.
}
}

\value{targeted imse value}
\references{
Picheny, V., Ginsbourger, D., Roustant, O., Haftka, R.T., Adaptive designs of experiments for accurate approximation of a target region, J. Mech. Des. - July 2010 - Volume 132, Issue 7, \url{http://dx.doi.org/10.1115/1.4001873}

Picheny V., Improving accuracy and compensating for uncertainty in surrogate modeling, Ph.D. thesis, University of Florida and Ecole Nationale Superieure des Mines de Saint-Etienne

Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2011), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set} ,\url{http://hal.archives-ouvertes.fr/hal-00641108/}
}

\author{
Victor Picheny (CERFACS, Toulouse, France) 

David Ginsbourger (IMSV, University of Bern, Switzerland) 
 
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

\seealso{\code{\link{EGI}}, \code{\link{max_timse}}}

\examples{
#timse_optim

set.seed(8)
N <- 9 #number of observations
T <- 80 #threshold
testfun <- branin

#a 9 points initial design
design <- data.frame( matrix(runif(2*N),ncol=2) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

###we need to compute some additional arguments:
#integration points, and current kriging means and variances at these points
integcontrol <- list(n.points=50,distrib="timse",init.distrib="MC")
obj <- integration_design(integcontrol=integcontrol,lower=c(0,0),
upper=c(1,1),model=model,T=T)

integration.points <- obj$integration.points
integration.weights <- obj$integration.weights
pred <- predict_nobias_km(object=model,newdata=integration.points,
type="UK",se.compute=TRUE)
intpoints.oldmean <- pred$mean ; intpoints.oldsd<-pred$sd

#another precomputation
precalc.data <- precomputeUpdateData(model,integration.points)

#we also need to compute weights. Otherwise the (more simple) 
#imse criterion will be evaluated
weight <- 1/sqrt(2*pi*intpoints.oldsd^2) * 
exp(-0.5*((intpoints.oldmean-T)/sqrt(intpoints.oldsd^2))^2)
weight[is.nan(weight)] <- 0

x <- c(0.5,0.4)#one evaluation of the timse criterion
timse_optim(x=x,integration.points=integration.points,
         integration.weights=integration.weights,
         intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
         precalc.data=precalc.data,T=T,model=model,weight=weight)

n.grid <- 20 #you can run it with 100
x.grid <- y.grid <- seq(0,1,length=n.grid)
x <- expand.grid(x.grid, y.grid)
timse.grid <- apply(X=x,FUN=timse_optim,MARGIN=1,integration.points=integration.points,
         integration.weights=integration.weights,
         intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
         precalc.data=precalc.data,T=T,model=model,weight=weight)

z.grid <- matrix(timse.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
image(x=x.grid,y=y.grid,z=z.grid,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid,15,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)

i.best <- which.min(timse.grid)
points(x[i.best,], col="blue", pch=17, lwd=4,cex=3)

#plots the real (unknown in practice) curve f(x)=T
testfun.grid <- apply(x,1,testfun)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE,lwd=5)
title("Contour lines of timse criterion (black) and of f(x)=T (blue)")
}
