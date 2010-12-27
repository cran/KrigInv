\name{bichon_optim}
\Rdversion{1.1}
\alias{bichon_optim}

\title{ Bichon et al.'s Expected Feasibility criterion }

\description{
Evaluation of Bichon's Expected Feasibility criterion. To be used in optimization routines, like in max_infill_criterion}

\usage{
bichon_optim(x, T, method.param, model, type="UK")
}

\arguments{
  \item{x}{ the input vector at which one wants to evaluate the criterion }
  \item{T}{ the targeted (scalar) output value } 
  \item{method.param}{ the (scalar) tolerance around the target T }
  \item{model}{ An object of class \code{\link[DiceKriging]{km}} (Kriging model)}
  \item{type}{ Kriging type (string): "SK" or "UK" (default) }
}
%\details{}
\value{Bichon EF criterion}
\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, accepted with minor revisions to the Journal of Statistics and Computing, \url{http://arxiv.org/abs/1009.5177}

Bichon, B.J., Eldred, M.S., Swiler, L.P., Mahadevan, S., McFarland, J.M.: Efficient global reliability analysis for nonlinear implicit performance functions. 
AIAA Journal 46 (10), 2459-2468 (2008)
}
\author{ 
V. Picheny (Ecole Centrale Paris)

D. Ginsbourger (IMSV, University of Bern, Switzerland)
}
%\note{}

\seealso{\code{\link{EGI}}, \code{\link{max_infill_criterion}}}

\examples{
##################################################################
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
bichon.EI.grid <- apply(design.grid, 1, bichon_optim, T=T, 
method.param=0.1, model=model)
z.grid <- matrix(bichon.EI.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
contour(x.grid,y.grid,z.grid,25)
points(design.fact, col="black", pch=20, lwd=4)

i.best <- which.max(bichon.EI.grid)
points(design.grid[i.best,], col="blue", pch=20, lwd=4)

#plots: contour of the actual function at threshold
testfun.grid <- testfun(design.grid)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE)
title("Contour lines of Bichon criterion (black) and of f(x)=T (blue)")
####################################################################
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. \keyword{}
%\keyword{}
