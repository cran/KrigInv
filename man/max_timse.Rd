\name{max_timse}
\Rdversion{1.1}
\alias{max_timse}

\title{ Minimizer of the targeted IMSE criterion}

\description{
Minimization, based on the package rgenoud (or on exhaustive search on a discrete set), of the targeted IMSE criterion. If not provided, the integration points (for the IMSE evaluation) are generated before optimization by latin hypercube sampling. The default number of integration points is 100 times the dimension of the problem.}

\usage{
max_timse(lower, upper, parinit=NULL, control=NULL, discrete.X=NULL, 
integration.points=NULL, T, epsilon=0, model, new.noise.var=0, type="UK")
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

  \item{epsilon}{ optional tolerance value (a real positive number); default value is 0.}

  \item{model}{ An object of class \code{\link[DiceKriging]{km}} (Kriging model)}
  
  \item{type}{ Kriging type (string): "SK" or "UK" (default) }
  
  \item{new.noise.var}{ optional scalar value of the noise variance for the new observations }

  \item{discrete.X}{ optional matrix of candidate points. If provided, the search for new observations is made on this discrete set instead of running the continuous optimisation, and it is also used as integration points }

  \item{integration.points}{ optional scalar or matrix. If it is scalar, it defines the number of integration points, then generated as an LHS design; if it is a matrix, it defines the integration points directly.
 					 If discrete.X is provided, this input is not taken into account.  }

}


\value{A list with components:
\item{par}{the best set of parameters found.}
\item{value}{the value tIMSE at par.}
}
\references{
Picheny, V., Ginsbourger, D., Roustant, O., Haftka, R.T., Adaptive designs of experiments for accurate approximation of a target region, J. Mech. Des. - July 2010 - Volume 132, Issue 7, http://dx.doi.org/10.1115/1.4001873

Picheny V., Improving accuracy and compensating for uncertainty in surrogate modeling, Ph.D. thesis, University of Florida and Ecole Nationale Superieure des Mines de Saint-Etienne
}
\author{ 
Victor Picheny 

David Ginsbourger}
%\note{}

\seealso{\code{\link{EGI}},\code{\link{timse_optim}}}

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

#compute weight on the grid
T <- 0
epsilon <- 0.01
krig  <- predict.km(model, newdata=as.data.frame(design.grid), "UK")
mk <- krig$mean
sk    <- krig$sd
weight <- 1/sqrt(2*pi*(sk^2+epsilon^2))*exp(-0.5*((mk-T)/sqrt(sk^2+epsilon^2))^2)
weight[is.nan(weight)] <- 0

#evaluate criterion on the grid		
timse.EI.grid <- apply(design.grid, 1, timse_optim, weight=weight, 
model=model, integration.points=design.grid)
z.grid <- matrix(timse.EI.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
contour(x.grid,y.grid,z.grid,25)
points(design.fact, col="black", pch=20, lwd=4)

#plots: contour of the actual function at threshold
testfun.grid <- testfun(design.grid)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE)
title("Contour lines of timse criterion (black) and of f(x)=T (blue)")

#search best point with Genoud
opt <- max_timse(lower=c(0,0), upper=c(1,1), T=T, epsilon=epsilon, 
model=model, integration.points=design.grid)
points(opt$par, col="blue", pch=20, lwd=4)

                      
#####################################################################
}
