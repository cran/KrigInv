\name{jn_optim}
\Rdversion{1.3}
\alias{jn_optim}

\title{ jn criterion optimization }

\description{
Evaluation of the \code{"jn"} criterion for a candidate point. 
To be used in the optimization routines \code{max_sur} with the argument \code{real.volume.variance=TRUE}. 
To avoid numerical instabilities, a new point is added to the design of experiments only if it is not too close to an existing observation, or if there is some observation noise. 
The criterion is the integral of the posterior expected \code{jn} uncertainty, which is the posterior expected variance of the excursion set's volume.
}

\usage{
jn_optim(x, integration.points, integration.weights = NULL, 
  intpoints.oldmean, intpoints.oldsd, 
  precalc.data, model, T, new.noise.var = NULL,current.sur=0)
}

\arguments{
  	\item{x}{Input vector of size d at which the criterion is evaluated.}
  	\item{integration.points}{p*d matrix of points for numerical integration in the X space.}
	\item{integration.weights}{Vector of size p corresponding to the weights of these integration points. If not provided, uniform weight is used.}
	\item{intpoints.oldmean}{Vector of size p corresponding to the kriging mean at the integration points before adding \code{x} to the design of experiments.}
	\item{intpoints.oldsd}{Vector of size p corresponding to the kriging standard deviation at the integration points before adding \code{x} to the design of experiments.}
	\item{precalc.data}{List containing useful data to compute quickly the updated kriging variance. This list can be generated using the \code{\link{precomputeUpdateData}} function.}
 	\item{model}{Object of class \code{\link[DiceKriging]{km}} (Kriging model).}
	\item{T}{Target value (scalar).}
	\item{new.noise.var}{Optional scalar value of the noise variance for the new observations.}
    \item{current.sur}{Arbitrary value given to the \code{"jn"} criterion for candidate points that are too close to existing observations. This argument applies only if the noise variance is zero.}
}

\value{jn value}
\references{
Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2011), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set} ,\url{http://hal.archives-ouvertes.fr/hal-00641108/}

Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing, pp.1-21, 2011, \url{http://arxiv.org/abs/1009.5177}
}

\author{ 
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

\seealso{\code{\link{EGI}}, \code{\link{max_sur}}, \code{\link{sur_optim}}}

\examples{
#jn_optim

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
integcontrol <- list(n.points=50,distrib="jn",init.distrib="MC")
obj <- integration_design(integcontrol=integcontrol,lower=c(0,0),upper=c(1,1),
       model=model,T=T)

integration.points <- obj$integration.points
integration.weights <- obj$integration.weights
pred <- predict_nobias_km(object=model,newdata=integration.points,
                          type="UK",se.compute=TRUE)
intpoints.oldmean <- pred$mean ; intpoints.oldsd<-pred$sd

#another precomputation
precalc.data <- precomputeUpdateData(model,integration.points)

x <- c(0.5,0.4)#one evaluation of the jn criterion
jn_optim(x=x,integration.points=integration.points,
         integration.weights=integration.weights,
         intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
         precalc.data=precalc.data,T=T,model=model)
#returns a negative value
#this is normal as a positive (constant) part of the criterion is not computed
#this constant part does not depend on x

n.grid <- 20 #you can run it with 100
x.grid <- y.grid <- seq(0,1,length=n.grid)
x <- expand.grid(x.grid, y.grid)
jn.grid <- apply(X=x,FUN=jn_optim,MARGIN=1,integration.points=integration.points,
         integration.weights=integration.weights,
         intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
         precalc.data=precalc.data,T=T,model=model)
z.grid <- matrix(jn.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
image(x=x.grid,y=y.grid,z=z.grid,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid,15,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)

i.best <- which.min(jn.grid)
points(x[i.best,], col="blue", pch=17, lwd=4,cex=3)

#plots the real (unknown in practice) curve f(x)=T
testfun.grid <- apply(x,1,testfun)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE,lwd=5)
title("Contour lines of jn criterion (black) and of f(x)=T (blue)")

}
