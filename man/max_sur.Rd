\name{max_sur}
\Rdversion{1.2}
\alias{max_sur}

\title{ Minimizer of the sur criterion}

\description{
Minimization, based on the package rgenoud (or on exhaustive search on a discrete set), of the sur criterion.  
}

\usage{
max_sur(lower, upper, optimcontrol = NULL, 
integration.param = NULL, T, model, 
new.noise.var = 0, real.volume.variance = FALSE, real.volume.constant = FALSE)
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
Optional scalar value of the noise variance at the new observation.
}
  \item{real.volume.variance}{
Boolean to decide if the \code{"jn"} sampling criterion should be used instead of the \code{"sur"} criterion. 
When it is equal to \code{FALSE} (default), the \code{"sur"} criterion is used (faster computation).
}
  \item{real.volume.constant}{
When the\code{"jn"} criterion is chosen, this argument decides whether a constant part of the \code{"jn"} criterion should be computed or not. 
Computing this constant does NOT change the optimum of the \code{"jn"} criterion. Default value: \code{FALSE}.
This argument is ignored if the argument \code{real.volume.variance} is set to \code{FALSE}.
}
}
\value{A list with components:
\item{par}{Best set of parameters found.}
\item{value}{Value of the criterion at par.}
\item{allvalues}{If an optimization on a discrete set of points is chosen, the value of the criterion at all these points}
\item{variance.volume}{If the arguments \code{real.volume.variance} and \code{real.volume.constant} are both set to \code{TRUE},the value of the computed constant}
}

\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing, pp.1-21, 2011, \url{http://arxiv.org/abs/1009.5177}

Vazquez, E., Bect, J.: A sequential Bayesian algorithm to estimate a probability of failure. In: Proceedings of the 15th IFAC Symposium on System Identification, (SYSID 2009), Saint-Malo, France (2009)

Chevalier C., Bect J., Ginsbourger D., Picheny V., Richet Y., Vazquez E. (2011), \emph{Fast kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}  ,\url{http://hal.archives-ouvertes.fr/docs/00/64/11/08/PDF/FastKrigingInversion.pdf}
}

\author{
Victor Picheny (CERFACS, Toulouse, France) 

David Ginsbourger (IMSV, University of Bern, Switzerland) 
 
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

\seealso{\code{\link{EGI}},\code{\link{sur_optim}}}

\examples{
#max_sur

set.seed(8)
N <- 9 #number of observations
T <- 80 #threshold
testfun <- branin
lower <- c(0,0)
upper <- c(1,1)

#a 9 points initial design (LHS in 2 dimensions)
design <- data.frame(maximinLHS(N,k=2))
response <- testfun(design)

#km object with matern3_2 covariance (params estimated by ML from the observations)
model <- km(formula=~., design = design, response = response,covtype="matern3_2")

optimcontrol <- list(method="genoud",pop.size=50)
integcontrol <- list(distrib="sur",n.points=50)
integration.param <- integration_design(integcontrol=integcontrol,d=2,
                                        lower=lower,upper=upper,model=model,
                                        T=T)

\dontrun{
obj <- max_sur(lower=lower,upper=upper,optimcontrol=optimcontrol,T=T,
                model=model,integration.param=integration.param)

obj$par;obj$value #optimum in (0.4550191,0.8614298)
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
