\name{sur_optim_parallel}
\Rdversion{1.2}
\alias{sur_optim_parallel}

\title{ Parallel sur criterion }
\description{
Evaluation of the parallel sur criterion for some candidate points. To be used in optimization routines, like in \code{\link{max_sur_parallel}}. 
To avoid numerical instabilities, the new points are evaluated only if they are not too close to an existing observation, or if there is some observation noise. 
The criterion is the integral of the posterior sur uncertainty.
}
\usage{
sur_optim_parallel(x, integration.points, integration.weights = NULL, 
intpoints.oldmean, intpoints.oldsd, 
precalc.data, model, T, 
new.noise.var = NULL, batchsize, current.sur)
}

\arguments{
  	\item{x}{Input vector of size batchsize*d at which one wants to evaluate the criterion. This argument is NOT a matrix.}
  	\item{integration.points}{p*d matrix of points for numerical integration in the X space.}
	\item{integration.weights}{Vector of size p corresponding to the weights of these integration points.}
  	\item{intpoints.oldmean}{Vector of size p corresponding to the kriging mean at the integration points before adding the batchsize points \code{x} to the design of experiments.}
	\item{intpoints.oldsd}{Vector of size p corresponding to the kriging standard deviation at the integration points before adding  the batchsize points \code{x} to the design of experiments.}
  	\item{precalc.data}{List containing useful data to compute quickly the updated kriging variance. This list can be generated using the \code{\link{precomputeUpdateData}} function.}
 	\item{model}{Object of class \code{\link[DiceKriging]{km}} (Kriging model).}
	\item{T}{Target value (scalar).}
  	\item{new.noise.var}{Optional scalar value of the noise variance for the new observations.}
  	\item{batchsize}{Number of points to sample simultaneously. The sampling criterion will return batchsize points at a time for sampling.}
  	\item{current.sur}{Current value of the sur criterion (before adding new observations)}
}
\details{
The first argument \code{x} has been chosen to be a vector of size batchsize*d (and not a matrix with batchsize rows and d columns) so that an optimizer like genoud can optimize it easily.
For example if d=2, batchsize=3 and \code{x=c(0.1,0.2,0.3,0.4,0.5,0.6)}, we will evaluate the parallel criterion at the three points (0.1,0.2),(0.3,0.4) and (0.5,0.6).
The last argument \code{current.sur} is used as a default value for the sur criterion when the new points \code{x} are too close to existing observations.
}

\value{Parallel sur value}
\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing, pp.1-21, 2011, \url{http://arxiv.org/abs/1009.5177}

Vazquez, E., Bect, J.: A sequential Bayesian algorithm to estimate a probability of failure. In: Proceedings of the 15th IFAC Symposium on System Identification, (SYSID 2009), Saint-Malo, France (2009)

Chevalier C., Bect J., Ginsbourger D., Picheny V., Richet Y., Vazquez E. (2011), \emph{Fast kriging-based stepwise uncertainty reduction with application to the identification of an excursion set},   ,\url{http://hal.archives-ouvertes.fr/docs/00/64/11/08/PDF/FastKrigingInversion.pdf}
}

\author{
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

\seealso{\code{\link{EGIparallel}}, \code{\link{max_sur_parallel}} }

\examples{
#sur_optim_parallel

set.seed(8)
N <- 9 #number of observations
T <- 80 #threshold
testfun <- branin

#a 9 points initial design (LHS in 2 dimensions)
design <- data.frame(maximinLHS(N,k=2))
response <- testfun(design)

#km object with matern3_2 covariance (params estimated by ML from the observations)
model <- km(formula=~., design = design, response = response,covtype="matern3_2")

###we need to compute some additional arguments:
#integration points, and current kriging means and variances at these points
integcontrol <- list(n.points=50,distrib="sur")
obj <- integration_design(integcontrol=integcontrol,
lower=c(0,0),upper=c(1,1),model=model,T=T)

integration.points <- obj$integration.points
integration.weights <- obj$integration.weights
pred <- predict_nobias_km(object=model,newdata=integration.points,
type="UK",se.compute=TRUE)
intpoints.oldmean <- pred$mean ; intpoints.oldsd<-pred$sd

#another precomputation
precalc.data <- precomputeUpdateData(model,integration.points)

batchsize <- 4
x <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
#one evaluation of the sur_optim_parallel criterion
#we calculate the expectation of the future "sur" uncertainty 
#when 4 points are added to the doe
#the 4 points are (0.1,0.2) , (0.3,0.4), (0.5,0.6), (0.7,0.8)
sur_optim_parallel(x=x,integration.points=integration.points,
          integration.weights=integration.weights,
          intpoints.oldmean=intpoints.oldmean,intpoints.oldsd=intpoints.oldsd,
          precalc.data=precalc.data,T=T,model=model,
          batchsize=batchsize,current.sur=Inf)
#return 0.07882718

#the function max_sur_parallel will help to find the optimum: 
#ie: the batch of 4 minimizing the expectation of the future uncertainty
}
