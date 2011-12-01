
\name{precomputeUpdateData}
\Rdversion{1.2}
\alias{precomputeUpdateData}

\title{Useful data to quickly update kriging mean and variance}

\description{
This function is used in combination with \code{\link{computeQuickKrigcov}} and computes an output list that serves as input in that function.
}

\usage{
precomputeUpdateData(model, integration.points)
}

\arguments{
  \item{model}{
A Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{integration.points}{
p*d matrix of points for numerical integration in the X space.
}
}

\value{A list with components:
\item{Kinv.c.olddata}{Matrix equal to K^(-1)*c where K is the non conditional covariance matrix at the design points and c is the non conditional covariances between the design points and the integration points.}
\item{Kinv.F}{Matrix equal to K^(-1)*F where F is a matrix with the values of the trend functions at the design points.}
\item{first.member}{Matrix with a complicated expression.}
}
\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing, pp.1-21, 2011, \url{http://arxiv.org/abs/1009.5177}

Chevalier C., Bect J., Ginsbourger D., Picheny V., Richet Y., Vazquez E. (2011), \emph{Fast kriging-based stepwise uncertainty reduction with application to the identification of an excursion set} ,\url{http://hal.archives-ouvertes.fr/docs/00/64/11/08/PDF/FastKrigingInversion.pdf}
}

\author{
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

\seealso{\code{\link{computeQuickKrigcov}}, \code{\link{predict_nobias_km}}, \code{\link{predict_update_km}}}

\examples{
#precomputeUpdateData

set.seed(8)
N <- 9 #number of observations
testfun <- branin

#a 9 points initial design (LHS in 2 dimensions)
design <- data.frame(maximinLHS(N,k=2))
response <- testfun(design)

#km object with matern3_2 covariance (params estimated by ML from the observations)
model <- km(formula=~., design = design, response = response,covtype="matern3_2")

#the points where we want to compute prediction (if a point new.x is added to the doe)
n.grid <- 20 #you can run it with 100
x.grid <- y.grid <- seq(0,1,length=n.grid)
newdata <- expand.grid(x.grid,y.grid)
precalc.data <- precomputeUpdateData(model=model,integration.points=newdata)

#now we can compute very quickly kriging covariances 
#between these data and any other points
other.x <- matrix(c(0.6,0.6),ncol=2)
pred <- predict_nobias_km(object=model,newdata=other.x,type="UK",se.compute=TRUE)

kn <- computeQuickKrigcov(model=model,integration.points=newdata,X.new=other.x,
                    precalc.data=precalc.data,F.newdata=pred$F.newdata,
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