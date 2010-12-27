\name{sur_optim}
\Rdversion{1.1}
\alias{sur_optim}

\title{ SUR criterion }

\description{
Evaluation of the SUR criterion for a candidate point. To be used in optimization routines, like in max.sur. The new point is added to the design of experiments only if it is not too close to an existing observation, or if there is some observation noise. The criterion is the integral of the posterior probability of misclassification.}

\usage{
sur_optim(x, y.integration=NULL, integration.points, model, T, 
new.noise.var=0, type="UK")
}

\arguments{
  \item{x}{ the input vector at which one wants to evaluate the criterion }
  \item{T}{ the targeted (scalar) output value } 
  \item{y.integration}{ the vector of integration points in the y space; if not provided, 10 points using importance sampling are used} 
  \item{model}{ An object of class \code{\link[DiceKriging]{km}} (Kriging model)}
  \item{type}{ Kriging type (string): "SK" or "UK" (default) }
  \item{integration.points}{ Matrix of points for numerical integration in the X space}
  \item{new.noise.var}{ optional scalar value of the noise variance for the new observation (default value is 0)}
}
%\details{}
\value{SUR value}
\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, accepted with minor revisions to the Journal of Statistics and Computing, \url{http://arxiv.org/abs/1009.5177}

Vazquez, E., Bect, J.: A sequential Bayesian algorithm to estimate a probability of failure. In: Proceedings of the 15th IFAC Symposium on System Identification, (SYSID 2009), Saint-Malo, France (2009)}
\author{ 
Victor Picheny

David Ginsbourger }
%\note{}

\seealso{\code{\link{EGI}}, \code{\link{max_sur}}}

\examples{
#####################################################################
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

#evaluate criterion on the grid
T <- 0
		
sur.grid <- apply(design.grid, 1, sur_optim, T=T, model=model, 
integration.points=design.grid)
z.grid <- matrix(sur.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
contour(x.grid,y.grid,z.grid,25)
points(design.fact, col="black", pch=20, lwd=4)

i.best <- which.min(sur.grid)
points(design.grid[i.best,], col="blue", pch=20, lwd=4)

#plots: contour of the actual function at threshold
testfun.grid <- testfun(design.grid)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE)
title("Contour lines of sur criterion (black) and of f(x)=T (blue)")
#####################################################################
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. \keyword{}
% \keyword{internal}
