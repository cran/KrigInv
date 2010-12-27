\name{max_sur}
\Rdversion{1.1}
\alias{max_sur}

\title{ Minimizer of the SUR criterion}

\description{
Minimization, based on the package rgenoud (or on exhaustive search on a discrete set), of the SUR criterion. If not provided, the integration points in X space are generated before optimization by latin hypercube sampling. 
The integration points in Y space are fixed. The default number of integration points in X is 100 times the dimension of the problem.}

\usage{
max_sur(lower, upper, parinit=NULL, control=NULL, discrete.X=NULL, 
n.int.y=10, integration.points=NULL, T, model, new.noise.var=0, type="UK")
}

\arguments{
  \item{lower}{ vector containing the lower bounds of the variables to be optimized over }

  \item{upper}{ vector containing the upper bounds of the variables to be optimized over }

  \item{parinit}{ optional vector containing the initial values for the variables to be optimized over }

   \item{control}{  optional list of control parameters for optimization. 
One can control  \code{"pop.size"}  (default : [4+3*log(nb of variables)]),  \code{"max.generations"} (5),  
\code{"wait.generations"} (2) and  \code{"BFGSburnin"} (0) of function  \code{"genoud"} (see \code{\link[rgenoud]{genoud}}). 
Numbers into brackets are the default values }  

  \item{T}{target value (a real number)} 

  \item{model}{ An object of class \code{\link[DiceKriging]{km}} (Kriging model)}

  \item{new.noise.var}{ optional scalar value of the noise variance for the new observations }
  
  \item{type}{ Kriging type (string): "SK" or "UK" (default) }

  \item{n.int.y}{ optional scalar value of the the number of integration points in the y space. }
  
  \item{discrete.X}{ optional matrix of candidate points. If provided, the search for new observations is made on this discrete set instead of running the continuous optimisation, and it is also used as integration points }

  \item{integration.points}{ optional scalar or matrix. If it is scalar, it defines the number of integration points, then generated as an LHS design; if it is a matrix, it defines the integration points directly. 
                                      If discrete.X is provided, this input is not taken into account. }

}


\value{A list with components:
\item{par}{the best set of parameters found.}
\item{value}{the value of the SUR criterion at par.}
}
\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, accepted with minor revisions to the Journal of Statistics and Computing, \url{http://arxiv.org/abs/1009.5177}

Vazquez, E., Bect, J.: A sequential Bayesian algorithm to estimate a probability of failure. In: Proceedings of the 15th IFAC Symposium on System Identification, (SYSID 2009), Saint-Malo, France (2009)}
\author{ Victor Picheny, David Ginsbourger}
%\note{}

\seealso{\code{\link{EGI}},\code{\link{sur_optim}}}

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

#plots: contour of the actual function at threshold
testfun.grid <- testfun(design.grid)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE)
title("Contour lines of sur criterion (black) and of f(x)=T (blue)")

#search best point with Genoud
opt <- max_sur(lower=c(0,0), upper=c(1,1), T=T, model=model, 
integration.points=design.grid)
points(opt$par, col="blue", pch=20, lwd=4)
              
###################################################################
}
