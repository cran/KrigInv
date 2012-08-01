\name{computeAuxVariables_update}
\Rdversion{1.3}
\alias{computeAuxVariables_update}

\title{
Auxiliary variables in an update_km procedure
}
\description{
Function similar to the \code{\link[DiceKriging]{computeAuxVariables}} of the DiceKriging package, with a quicker implementation.
}
\usage{
computeAuxVariables_update(model)
}

\arguments{
  \item{model}{
A Kriging model of \code{\link[DiceKriging]{km}} class.
}
}
\details{
This function was introduced to optimize the calculation time of some expensive to evaluate integral criteria.
}
\value{
An updated \code{\link[DiceKriging]{km}} model
}

\author{
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

\seealso{
\code{\link[DiceKriging]{computeAuxVariables}}
}

\examples{
#computeAuxVariables_update

set.seed(8)
N <- 9 #number of observations
testfun <- branin

#a 9 points initial design (LHS in 2 dimensions)
design <- data.frame(maximinLHS(N,k=2))
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

model2 <- computeAuxVariables_update(model)
#in this example model2 and model are the same!
}
