\name{vorob_threshold}
\Rdversion{1.3}
\alias{vorob_threshold}

\title{ Calculation of the Vorob'ev threshold }
\description{
Evaluation of the Vorob'ev threshold given an excursion probability vector. This threshold is such that the volume of the set (x : pn(x) > threshold) is equal to the integral of pn. 
}
\usage{
vorob_threshold(pn)
}

\arguments{
  	\item{pn}{Input vector of arbitrary size containing the excursion probabilities pn(x).}
}
\details{
In this function, all the points x are supposed tobe equaly weighted.
}

\value{a scalar: the Vorob'ev thresold}
\references{
Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2011), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set} ,\url{http://hal.archives-ouvertes.fr/hal-00641108/}

Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing, pp.1-21, 2011, \url{http://arxiv.org/abs/1009.5177}
}

\author{
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

\seealso{\code{\link{max_vorob_parallel}}, \code{\link{vorob_optim_parallel}} }

\examples{
#vorob_threshold

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

\dontrun{
###we need to compute some additional arguments:
#integration points, and current kriging means and variances at these points
integcontrol <- list(n.points=50,distrib="sobol")
obj <- integration_design(integcontrol=integcontrol,
lower=c(0,0),upper=c(1,1),model=model,T=T)

integration.points <- obj$integration.points

pred <- predict_nobias_km(object=model,newdata=integration.points,
type="UK",se.compute=TRUE)
pn <- pnorm((pred$mean-T)/pred$sd)

vorob_threshold(pn)
}
}
