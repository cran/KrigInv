\name{tmse_optim}
\Rdversion{1.1}
\alias{tmse_optim}

\title{ Targeted MSE criterion }

\description{
Evaluation of the Targeted MSE criterion. To be used in optimization routines, like in max.infill.criterion}

\usage{
tmse_optim(x, T, method.param=0, model, type="UK")
}

\arguments{
  \item{x}{ the input vector at which one wants to evaluate the criterion }
  \item{T}{ the targeted (scalar) output value } 
  \item{method.param}{ the (scalar) tolerance around the target T }
  \item{model}{ An object of class \code{\link[DiceKriging]{km}} (Kriging model)}
  \item{type}{ Kriging type (string): "SK" or "UK" (default) }
}
%\details{}
\value{targeted MSE value}
\references{
Picheny, V., Ginsbourger, D., Roustant, O., Haftka, R.T., Adaptive designs of experiments for accurate approximation of a target region, J. Mech. Des. - July 2010 - Volume 132, Issue 7, \url{http://dx.doi.org/10.1115/1.4001873}

Picheny V., Improving accuracy and compensating for uncertainty in surrogate modeling, Ph.D. thesis, University of Florida and Ecole Nationale Superieure des Mines de Saint-Etienne
}
\author{ Victor Picheny

David Ginsbourger }
%\note{}

\seealso{\code{\link{EGI}}, \code{\link{max_infill_criterion}}}

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
n.grid <- 30
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)

#evaluate criterion on the grid
T <- 0
tmse.EI.grid <- apply(design.grid, 1, tmse_optim, T=T, method.param=
0.1, model=model)
z.grid <- matrix(tmse.EI.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
contour(x.grid,y.grid,z.grid,25)
points(design.fact, col="black", pch=20, lwd=4)

i.best <- which.max(tmse.EI.grid)
points(design.grid[i.best,], col="blue", pch=20, lwd=4)

#plots: contour of the actual function at threshold
testfun.grid <- testfun(design.grid)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE)
title("Contour lines of tmse criterion (black) and of 
f(x)=T (blue)")
#####################################################################
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. \keyword{}
%\keyword{internal}
