\name{ranjan_optim}
\Rdversion{1.1}
\alias{ranjan_optim}

\title{ Ranjan et al.'s Expected Improvement criterion }

\description{
Evaluation of Ranjan's Expected Improvement criterion. To be used in optimization routines, like in max.infill.criterion}

\usage{
ranjan_optim(x, T, method.param, model, type="UK")
}

\arguments{
  \item{x}{ the input vector at which one wants to evaluate the criterion }
  \item{T}{ the targeted (scalar) output value } 
  \item{method.param}{ the (scalar) tolerance around the target T }
  \item{model}{ An object of class \code{\link[DiceKriging]{km}} (Kriging model).}
  \item{type}{ Kriging type (string): "SK" or "UK" (default) }
}
%\details{}
\value{Ranjan EI criterion}
\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, accepted with minor revisions to the Journal of Statistics and Computing, \url{http://arxiv.org/abs/1009.5177}

Ranjan, P., Bingham, D., Michailidis, G.: Sequential experiment design for contour estimation from complex computer codes. Technometrics 50(4), 527-541 (2008)}
\author{
Victor Picheny 

David Ginsbourger }
%\note{}

\seealso{\code{\link{EGI}}, \code{\link{max_infill_criterion}}}

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
n.grid <- 30
x.grid <- y.grid <- seq(0,1,length=n.grid)
design.grid <- expand.grid(x.grid, y.grid)

#evaluate criterion on the grid
T <- 0
ranjan.EI.grid <- apply(design.grid, 1, ranjan_optim, T=T, 
method.param=0.1, model=model)
z.grid <- matrix(ranjan.EI.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
contour(x.grid,y.grid,z.grid,25)
points(design.fact, col="black", pch=20, lwd=4)

i.best <- which.max(ranjan.EI.grid)
points(design.grid[i.best,], col="blue", pch=20, lwd=4)

#plots: contour of the actual function at threshold
testfun.grid <- testfun(design.grid)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE)
title("Contour lines of Ranjan criterion (black) and of f(x)=T (blue)")
#####################################################################
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. \keyword{}
% \keyword{}
