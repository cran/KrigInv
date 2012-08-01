
\name{computeRealVolumeConstant}
\Rdversion{1.3}
\alias{computeRealVolumeConstant}

\title{A constant used to calculate the expected excursion set's volume variance}

\description{
This function computes a constant used to calculate exactly the value of the \code{"jn"} criterion at one point. 
Computing this constant does NOT change the optimum of the \code{"jn"} criterion. 
Therefore, its calculation is indicative only and is only necessary to know exactly (in expectation) the excursion set's volume variance. 
}

\usage{
computeRealVolumeConstant(model,integration.points,
integration.weights=NULL,T)
}

\arguments{

  \item{model}{A Kriging model of \code{\link[DiceKriging]{km}} class.}
	\item{integration.points}{p*d matrix of points for numerical integration in the X space.}
	\item{integration.weights}{(Optional) Vector of size p corresponding to the weights of these integration points. If not provided, all weights are set to 1.}
 	\item{T}{ the targeted (scalar) output value }
 	
}

\value{a scalar
}

\references{
Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2011), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set} ,\url{http://hal.archives-ouvertes.fr/hal-00641108/}

Chevalier C., Ginsbourger D. (2012), \emph{Corrected Kriging update formulae for batch-sequential data assimilation} ,\url{http://arxiv.org/pdf/1203.6452.pdf}
}

\author{
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

\seealso{\code{\link{precomputeUpdateData}}, \code{\link{predict_nobias_km}}, \code{\link{predict_update_km}}}

\examples{
#computeRealVolumeConstant

set.seed(8)
N <- 9 #number of observations
testfun <- branin
T <- 80

#a 9 points initial design (LHS in 2 dimensions)
design <- data.frame(maximinLHS(N,k=2))
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

integcontrol <- list(n.points=50,distrib="sur")
obj <- integration_design(integcontrol=integcontrol,
	lower=c(0,0),upper=c(1,1),
model=model,T=T)
integration.points <- obj$integration.points
integration.weights <- obj$integration.weights

computeRealVolumeConstant(model=model,
	integration.points=integration.points,
        integration.weights=integration.weights,T=T)

}
