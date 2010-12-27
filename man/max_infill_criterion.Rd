\name{max_infill_criterion}
\alias{max_infill_criterion}
\title{Optimizer for the infill criteria}
\Rdversion{1.1}

\description{Global optimization, based on the package rgenoud (or on exhaustive search on a discrete set), of the chosen infill criterion (maximization or minimization, depending on the case)}


\usage{
max_infill_criterion(lower, upper, parinit=NULL, sampling.method, 
method.param=0, control=NULL, discrete.X=NULL, T, model, type="UK")
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

  \item{method.param}{optional tolerance value (a real number)}

  \item{model}{ An object of class \code{\link[DiceKriging]{km}} (Kriging model)}
  
  \item{type}{ Kriging type (string): "SK" or "UK" (default) }

  \item{sampling.method}{ criterion used for choosing observations: "ranjan" (default) , "bichon", or "tmse" }

  \item{discrete.X}{ optional matrix of candidate points. If provided, the search for new observations is made on this discrete set instead of running the continuous optimisation }

}


\value{A list with components:
\item{par}{the best set of parameters found.}
\item{value}{the value of the chosen criterion at par.}
}
\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, accepted with minor revisions to the Journal of Statistics and Computing, \url{http://arxiv.org/abs/1009.5177}

Picheny, V., Ginsbourger, D., Roustant, O., Haftka, R.T., Adaptive designs of experiments for accurate approximation of a target region, J. Mech. Des. - July 2010 - Volume 132, Issue 7, http://dx.doi.org/10.1115/1.4001873

Bichon, B.J., Eldred, M.S., Swiler, L.P., Mahadevan, S., McFarland, J.M.: Efficient global reliability analysis for nonlinear implicit performance functions. AIAA Journal 46 (10), 2459-2468 (2008)

Ranjan, P., Bingham, D., Michailidis, G.: Sequential experiment design for contour estimation from complex computer codes. Technometrics 50(4), 527-541 (2008)
}
\author{ Victor Picheny, David Ginsbourger}
%\note{}

\seealso{\code{\link{EGI}},\code{\link{ranjan_optim}},\code{\link{tmse_optim}},\code{\link{bichon_optim}}}

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

#evaluate Bichon's criterion on the grid
T <- 0
bichon.EI.grid <- apply(design.grid, 1, bichon_optim, T=T, 
method.param=0.1, model=model)
z.grid <- matrix(bichon.EI.grid, n.grid, n.grid)

#plots: contour of the criterion, doe points and new point
contour(x.grid,y.grid,z.grid,25)
points(design.fact, col="black", pch=20, lwd=4)

#plots: contour of the actual function at threshold
testfun.grid <- testfun(design.grid)
z.grid.2 <- matrix(testfun.grid, n.grid, n.grid)
contour(x.grid,y.grid,z.grid.2,levels=T,col="blue",add=TRUE)
title("Contour lines of Bichon criterion (black) and of f(x)=T (blue)")

#search best point with Genoud
opt <- max_infill_criterion(lower=c(0,0), upper=c(1,1), 
sampling.method="bichon", method.param=0.1, T=T, model=model)
points(opt$par, col="blue", pch=20, lwd=4)

##################################################################
}
