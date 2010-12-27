\name{timse_optim}
\Rdversion{1.1}
\alias{timse_optim}

\title{ Targeted IMSE criterion }

\description{
Evaluation of the Targeted IMSE criterion for a candidate point. To be used in optimization routines, like in max.timse. The weight function is defined outside before the evaluation of this function. The new point is added to the design of experiments only if it is not too close to an existing observation, or if there is some observation noise.}

\usage{
timse_optim(x, weight=1, integration.points, model, new.noise.var=0, 
type="UK")
}

\arguments{
  \item{x}{ the input vector at which one wants to evaluate the criterion }
  \item{weight}{ the vector of weight function (length must be equal to the number of lines of integration.points } 
  \item{model}{ An object of class \code{\link[DiceKriging]{km}} (Kriging model).}
  \item{type}{ Kriging type (string): "SK" or "UK" (default) }
  \item{integration.points}{ Matrix of points for numerical integration (at which the weight vector is defined) }
  \item{new.noise.var}{ optional scalar value of the noise variance for the new observation }
}
%\details{}
\value{targeted IMSE value}
\references{
Picheny, V., Ginsbourger, D., Roustant, O., Haftka, R.T., Adaptive designs of experiments for accurate approximation of a target region, J. Mech. Des. - July 2010 - Volume 132, Issue 7, \url{http://dx.doi.org/10.1115/1.4001873}

Picheny V., Improving accuracy and compensating for uncertainty in surrogate modeling, Ph.D. thesis, University of Florida and Ecole Nationale Superieure des Mines de Saint-Etienne
}
\author{ 
Victor Picheny

David Ginsbourger }
%\note{}

\seealso{\code{\link{EGI}}, \code{\link{max_timse}}}

\examples{
####################################################################
#a 9-point full factorial initial design
design.fact <- expand.grid(seq(0,1,length=3), seq(0,1,length=3))

design.fact <- data.frame(design.fact)
names(design.fact) <- c ( "x1","x2")
testfun <- camelback2			#our test function

#the response
response <- testfun(design.fact)

#the initial km model
model <- km(formula=~., design = design.fact, response = response, 
covtype="matern5_2")

#the integration points
n.grid <- 20
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)

#compute weight on the grid
T <- 0
epsilon <- 0.01
krig  <- predict.km(model, newdata=as.data.frame(design.grid), "UK")
mk <- krig$mean
sk    <- krig$sd
weight <- 1/sqrt(2*pi*(sk^2+epsilon^2))*exp(-0.5*((mk-T)/sqrt(sk^2
+epsilon^2))^2)
weight[is.nan(weight)] <- 0

#evaluate criterion on the grid		
timse.EI.grid <- apply(design.grid, 1, timse_optim, weight=weight, 
model=model, integration.points=design.grid)
z.grid <- matrix(timse.EI.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
contour(x.grid,y.grid,z.grid,25)
points(design.fact, col="black", pch=20, lwd=4)

i.best <- which.min(timse.EI.grid)
points(design.grid[i.best,], col="blue", pch=20, lwd=4)

#plots: contour of the actual function at threshold
testfun.grid <- testfun(design.grid)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE)
title("Contour lines of timse criterion (black) and of f(x)=T 
(blue)")
####################################################################
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. \keyword{}
% \keyword{internal}
