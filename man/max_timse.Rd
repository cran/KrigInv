\name{max_timse}
\Rdversion{1.3}
\alias{max_timse}

\title{ Minimizer of the IMSE or targeted IMSE criterion}

\description{
Minimization, based on the package rgenoud (or on exhaustive search on a discrete set), of the targeted imse (or imse) criterion.
}

\usage{
max_timse(lower, upper, optimcontrol = NULL, 
integration.param = NULL, T, model, 
new.noise.var = 0, epsilon = 0, imse = FALSE)
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
Optional scalar value of the noise variance of the new observation.
}
  \item{epsilon}{
Optional tolerance value (a real positive number). Default value is 0.
}
  \item{imse}{
Optional boolean to decide if the \code{"imse"} criterion should be used instead of \code{"timse"}. Default: \code{FALSE}.
}
}

\value{A list with components:
\item{par}{the best set of parameters found.}
\item{value}{the value of the criterion at par.}
\item{allvalues}{if an optimization on a discrete set of points is chosen, the value of the criterion at all these points}
}

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

\seealso{\code{\link{EGI}},\code{\link{timse_optim}}}

\examples{
#max_timse

set.seed(8)
N <- 9 #number of observations
T <- 80 #threshold
testfun <- branin
lower <- c(0,0)
upper <- c(1,1)

#a 9 points initial design
design <- data.frame( matrix(runif(2*N),ncol=2) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

optimcontrol <- list(method="genoud",pop.size=50)
integcontrol <- list(distrib="timse",n.points=50,init.distrib="MC")
integration.param <- integration_design(integcontrol=integcontrol,d=2,
                                            lower=lower,upper=upper,model=model,
                                            T=T)

\dontrun{
obj <- max_timse(lower=lower,upper=upper,optimcontrol=optimcontrol,T=T,
                model=model,integration.param=integration.param)

obj$par;obj$value
new.model <- update_km(model=model,NewX=obj$par,NewY=testfun(obj$par),
                       CovReEstimate=TRUE)

par(mfrow=c(1,2))
print_uncertainty(model=model,T=T,type="pn",lower=lower,upper=upper,
cex.points=2.5,main="probability of excursion")

print_uncertainty(model=new.model,T=T,type="pn",lower=lower,upper=upper,
new.points=1,col.points.end="red",cex.points=2.5,
main="updated probability of excursion")
}
}
