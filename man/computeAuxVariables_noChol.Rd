\name{computeAuxVariables_noChol}
\Rdversion{1.3}
\alias{computeAuxVariables_noChol}

\title{Auxiliary variables for kriging}


\description{
Computes or updates some auxiliary variables used for kriging (see below). 
This function is a copy of the computeAuxVariables function from the DiceKriging package, except that the calculation of the Cholesky decomposition is not performed, for cpu time savings.
} 

\usage{
computeAuxVariables_noChol(model)
}
\arguments{
  \item{model}{An object of class \code{\link[DiceKriging]{km}} with missing or non updated fields. }
}

\value{
An updated \code{\link[DiceKriging]{km}} objet, where the changes concern the following field:
\item{z}{Vector equal to t(T)^(-1)*(y - F*beta), where y, F and Beta are respectively the vector of response, the experimental matrix and the trend coefficients specified in \code{model}.
}
}
  
\references{ 
J.-S. Park and J. Baek (2001), Efficient computation of maximum likelihood estimators in a spatial linear model with power exponential covariogram, \emph{Computer Geosciences}, \bold{27} no. 1, 1-7. 
}
                      
\author{O. Roustant, D. Ginsbourger, Ecole des Mines de St-Etienne}

\seealso{ \code{\link[base]{chol}}, \code{\link[base]{backsolve}}. }

\keyword{models}

\examples{
#computeAuxVariables_noChol
set.seed(8)
N <- 9 #number of observations
testfun <- branin

#a 9 points initial design
design <- data.frame( matrix(runif(2*N),ncol=2) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

model2 <- computeAuxVariables_noChol(model)
#in this example model2 and model are the same!
}
