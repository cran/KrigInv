\name{computeQuickKrigcov}
\Rdversion{1.3}
\alias{computeQuickKrigcov}

\title{Quick computation of kriging covariances}

\description{
Computes kriging covariances between one new point and many integration points, using precomputed data.
}

\usage{
computeQuickKrigcov(model,integration.points,X.new,
precalc.data, F.newdata , c.newdata)
}

\arguments{
 
  \item{model}{A Kriging model of \code{\link[DiceKriging]{km}} class.}
	\item{integration.points}{p*d matrix of points for numerical integration in the X space.}
 	\item{X.new}{The new point where we calculate kriging covariances. The calculated covariances are the covariances between this new point and all the integration points.}
 	\item{precalc.data}{List containing precalculated data. This list is generated using the function \code{\link{precomputeUpdateData}} }
 	\item{F.newdata}{The value of the kriging trend basis function at point X.new}
 	\item{c.newdata}{The (unconditional) covariance between X.new and the design points}

}

\details{
This function requires to use another function in order to generate the proper arguments. 
The argument \code{precalc.data} can be generated using \code{\link{precomputeUpdateData}}.
The arguments \code{F.newdata} and \code{c.newdata} can be obtained using \code{\link{predict_nobias_km}}, which returns a field \code{F.newdata} and a field \code{c}.
}

\value{A vector containing kriging covariances
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
#computeQuickKrigcov

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

#the points where we want to compute prediction 
#if a point new.x is added to the doe
n.grid <- 20 #you can run it with 100
x.grid <- y.grid <- seq(0,1,length=n.grid)
newdata <- expand.grid(x.grid,y.grid)

#precalculation
precalc.data <- precomputeUpdateData(model=model,
		integration.points=newdata)

#now we can compute very quickly kriging covariances 
#between these data and any other points
other.x <- matrix(c(0.6,0.6),ncol=2)
pred <- predict_nobias_km(object=model,
	newdata=other.x,type="UK",se.compute=TRUE)

kn <- computeQuickKrigcov(model=model,
	integration.points=newdata,X.new=other.x,
        precalc.data=precalc.data,
	F.newdata=pred$F.newdata,
        c.newdata=pred$c)

z.grid <- matrix(kn, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
image(x=x.grid,y=y.grid,z=z.grid,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid,15,add=TRUE)
contour(x=x.grid,y=y.grid,z=z.grid,levels=0,add=TRUE,col="blue",lwd=5)
points(design, col="black", pch=17, lwd=4,cex=2)
points(other.x, col="red", pch=17, lwd=4,cex=3)
title("Kriging covariances with the point (0.6,0.6), in red")
}
