\name{predict_update_km}
\Rdversion{1.3}
\alias{predict_update_km}

\title{
Quick computation of updated kriging means and variances.
}

\description{
This function uses kriging update formula to quickly compute kriging mean and variances at points \code{newdata}, when a point \code{newX} is added. 
The required additional informations are the old kriging mean and variance at point \code{newX}, the output value f(newX), the old kriging mean and variances at points \code{newdata} and the kriging covariance between \code{newX} and {newdata}.
}
\usage{
predict_update_km(newXmean, newXvar, newXvalue, 
newdata.oldmean, newdata.oldsd, kn)
}

\arguments{
  \item{newXmean}{
Scalar: old kriging mean at \code{newX} (before adding \code{newX} to the design).
}
  \item{newXvar}{
Scalar: old kriging variance at \code{newX} (before adding \code{newX} to the design).
}
  \item{newXvalue}{
Scalar: value of the objective function at \code{newX}, f(\code{newX}).
}
  \item{newdata.oldmean}{
Vector: old kriging mean at the points \code{newdata} (before adding \code{newX} to the design)
}
  \item{newdata.oldsd}{
Vector: old kriging standard deviations at the points \code{newdata} (before adding \code{newX} to the design)
}
  \item{kn}{
Kriging covariances between the points \code{newdata} and \code{newX}. These covariances can be computed using the function \code{\link{computeQuickKrigcov}}
}
}

\value{A list with the following fields:
\item{mean}{Updated kriging mean at points \code{newdata}}
\item{sd}{Updated kriging standard deviation at points \code{newdata}}
\item{lambda}{New kriging weight of \code{newX} for the prediction at points \code{newdata}}
}

\references{
Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2011), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set} ,\url{http://hal.archives-ouvertes.fr/hal-00641108/}

Chevalier C., Ginsbourger D. (2012), \emph{Corrected Kriging update formulae for batch-sequential data assimilation} ,\url{http://arxiv.org/pdf/1203.6452.pdf}
}

\author{
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

\seealso{\code{\link{precomputeUpdateData}}, \code{\link{predict_nobias_km}}, \code{\link{computeQuickKrigcov}}}

\examples{
#predict_update_km

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

#points where we want to compute prediction (if a point new.x is added to the doe)
n.grid <- 20 #you can run it with 100
x.grid <- y.grid <- seq(0,1,length=n.grid)
newdata <- expand.grid(x.grid,y.grid)
precalc.data <- precomputeUpdateData(model=model,integration.points=newdata)
pred2 <- predict_nobias_km(object=model,newdata=newdata,type="UK",se.compute=TRUE)
newdata.oldmean <- pred2$mean; newdata.oldsd <- pred2$sd

new.x <- matrix(c(0.6,0.6),ncol=2)   #the point that we are going to add
pred1 <- predict_nobias_km(object=model,newdata=new.x,type="UK",se.compute=TRUE)
newXmean <- pred1$mean; newXvar <- pred1$sd^2; newXvalue <- pred1$mean + 2*pred1$sd

kn <- computeQuickKrigcov(model=model,integration.points=newdata,X.new=new.x,
                    precalc.data=precalc.data,F.newdata=pred1$F.newdata,
                    c.newdata=pred1$c)

updated.predictions <- predict_update_km(newXmean=newXmean,newXvar=newXvar,
newXvalue=newXvalue,newdata.oldmean=newdata.oldmean,
newdata.oldsd=newdata.oldsd,kn=kn)

#the new kriging variance is usually lower than the old one
updated.predictions$sd - newdata.oldsd 

z.grid1 <- matrix(newdata.oldsd, n.grid, n.grid)
z.grid2 <- matrix(updated.predictions$sd, n.grid, n.grid)

par(mfrow=c(1,2))

#plots: contour of the criterion, doe points and new point
image(x=x.grid,y=y.grid,z=z.grid1,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid1,15,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)
title("Kriging standard deviation")

image(x=x.grid,y=y.grid,z=z.grid2,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid2,15,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)
points(new.x, col="red", pch=17, lwd=4,cex=2)
title("updated Kriging standard deviation")

}
