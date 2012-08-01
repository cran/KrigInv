\name{max_sur_parallel}
\Rdversion{1.3}
\alias{max_sur_parallel}

\title{ Minimizer of the parallel sur criterion}

\description{
Minimization, based on the package rgenoud (or on exhaustive search on a discrete set), of the sur criterion for a batch of candidate sampling points.  
}

\usage{
max_sur_parallel(lower, upper, optimcontrol = NULL, 
batchsize, integration.param, T,
model, new.noise.var = 0)
}

\arguments{
  \item{lower}{
Vector containing the lower bounds of the design space.
}
  \item{upper}{
Vector containing the upper bounds of the design space.
}
  \item{optimcontrol}{
Optional list of control parameters for the optimization of the sampling criterion. The field \code{method} defines which optimization method is used: it can be either \code{"genoud"} (default) for an optimisation using the genoud algorithm, or \code{"discrete"} for an optimisation over a specified discrete set. 
If the field \code{method} is set to \code{"genoud"}, one can set some parameters of this algorithm: 
\code{pop.size}  (default : 50*d),  \code{max.generations} (10*d), \code{wait.generations} (2),  \code{BFGSburnin} (2) and the mutations \code{P1}, \code{P2}, up to \code{P9} (see \code{\link[rgenoud]{genoud}}). Numbers into brackets are the default values.
If the field \code{method} is set to \code{"discrete"}, one can set the field \code{optim.points}: p * d matrix corresponding to the p points where the criterion will be evaluated. If nothing is specified, 100*d points are chosen randomly. 
Finally, one can control the field {optim.option} in order to decide how to optimize the sampling criterion. 
If \code{optim.option} is set to 2 (default), batchsize sequential optimizations in dimension d are performed to find the optimum.
If \code{optim.option} is set to 1, only one optimization in dimension batchsize*d is performed. This option is only available with \code{"genoud"}. This option might provide more global and accurate solutions, but is a lot more expensive.
}
  \item{batchsize}{
Number of points to sample simultaneously. The sampling criterion will return batchsize points at a time for sampling.
}
  \item{integration.param}{
Optional list of control parameter for the computation of integrals, containing the fields \code{integration.points}: a p*d matrix corresponding to p integrations points and \code{integration.weights}: a vector of size p corresponding to the weights of these integration points. 
If nothing is specified, default values are used (see: function \code{\link{integration_design}} for more details). 
}
  \item{T}{
Target value (scalar).
}
  \item{model}{
A Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{new.noise.var}{
Optional scalar value of the noise variance of the new observations.
}
}

\value{A list with components:
\item{par}{the best set of points found.}
\item{value}{the value of the sur criterion at par.}
\item{allvalues}{If an optimization on a discrete set of points is chosen, the value of the criterion at all these points.}
}

\references{
Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2011), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set} ,\url{http://hal.archives-ouvertes.fr/hal-00641108/}

Chevalier C., Ginsbourger D. (2012), \emph{Corrected Kriging update formulae for batch-sequential data assimilation} ,\url{http://arxiv.org/pdf/1203.6452.pdf}
}

\author{
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

\seealso{\code{\link{EGIparallel}},\code{\link{sur_optim_parallel}}}

\examples{
#max_sur_parallel


set.seed(8)
N <- 9 #number of observations
T <- 80 #threshold
testfun <- branin
lower <- c(0,0)
upper <- c(1,1)

#a 9 points initial design (LHS in 2 dimensions)
design <- data.frame(maximinLHS(N,k=2))
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

optimcontrol <- list(method="genoud",pop.size=50,optim.option=1)
integcontrol <- list(distrib="sur",n.points=50)
integration.param <- integration_design(integcontrol=integcontrol,d=2,
                                            lower=lower,upper=upper,model=model,
                                            T=T)
batchsize <- 5 #number of new points

\dontrun{
obj <- max_sur_parallel(lower=lower,upper=upper,optimcontrol=optimcontrol,
                        batchsize=batchsize,T=T,model=model,
                        integration.param=integration.param)
                        #one optim in dimension 5*2 !

obj$par;obj$value #optimum in 5 new points
new.model <- update_km(model=model,NewX=obj$par,NewY=apply(obj$par,1,testfun),
                       CovReEstimate=TRUE)


par(mfrow=c(1,2))
print_uncertainty(model=model,T=T,type="pn",lower=lower,upper=upper,
cex.points=2.5,main="probability of excursion")

print_uncertainty(model=new.model,T=T,type="pn",lower=lower,upper=upper,
new.points=batchsize,col.points.end="red",cex.points=2.5,
main="updated probability of excursion")
}
}
