
\name{EGI}
\Rdversion{1.1}
\alias{EGI}

\title{Efficient Global Inversion: sequential inversion algorithm based on Kriging}

\description{
Sequential sampling based on the optimization of a kriging-based criterion, with model update after each evaluation. Five criteria are available for selecting experiments, three inexpensive ("bichon", "ranjan", and "tmse") and two expensive ones ("timse" and "sur", which require numerical integration). }

\usage{
EGI(T, model, method.param, sampling.method, discrete.X,
new.noise.var, fun, iter, lower, upper, integration.points,
parinit, control, kmcontrol)
}


\arguments{
  \item{T}{target value (a real number). The sampling algorithm and the underlying kriging model aim to find the points below (resp. over) T.} 

  \item{model}{ a Kriging model of "km" class }

  \item{method.param}{ optional tolerance value (a real number) for methods "ranjan", "bichon", "tmse" and "timse"}

  \item{sampling.method}{ Criterion used for choosing observations: "ranjan" (default) , "bichon", "tmse", "timse" or "sur".}

  \item{discrete.X}{ optional matrix of candidate points. If provided, the search of the optimum for new observations is made on this discrete set instead of running the continuous optimisation. For "sur" and "timse", discrete.X also serves as integration points (and then the argument "integration.points" is ignored)}

  \item{new.noise.var}{ optional scalar value of the noise variance for the new observations }

  \item{fun}{ objective function }

  \item{iter}{ Number of iterations }

  \item{lower}{ vector containing the lower bounds of the variables to be optimized over }

  \item{upper}{ vector containing the upper bounds of the variables to be optimized over }

  \item{integration.points}{ optional scalar or matrix for "sur" and "timse" criteria. If it is scalar, it defines the number of integration points; if it is a matrix, it defines the integration points}

  \item{parinit}{ optional vector of initial values for the variables to be optimized over }

   \item{control}{  optional list of control parameters for the optimization of the variables. 
One can control  \code{"pop.size"}  (default : [4+3*log(nb of variables)]),  \code{"max.generations"} (5),  
\code{"wait.generations"} (2) and  \code{"BFGSburnin"} (0) of function  \code{"genoud"} (see \code{\link[rgenoud]{genoud}}). 
Numbers into brackets are the default values }  

    \item{kmcontrol}{ an optional list representing the control variables for the re-estimation of the kriging model once new points are sampled. 
		    The items are the same as in \code{\link[DiceKriging]{km}} : 

		\code{penalty}, \code{optim.method}, \code{parinit}, \code{control}. 

		The default values are those contained in \code{model}, typically corresponding to the variables 
		used in \code{\link[DiceKriging]{km}} to estimate a kriging model from the initial design points.}
}


\value{A list with components:
\item{par}{ The added observations}
\item{value}{the added observation values}
\item{nsteps}{the number of added observations}
\item{lastmodel}{the current (last) kriging model of "km" class}
}
%\references{}
\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, accepted with minor revisions to the Journal of Statistics and Computing, \url{http://arxiv.org/abs/1009.5177}
}

\author{
Victor Picheny (Ecole Centrale Paris, France) 

David Ginsbourger (IMSV, University of Bern, Switzerland) 
 
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

%\note{}

\seealso{\code{\link{max_sur}}, \code{\link{max_timse}}, \code{\link{max_infill_criterion}}}

\examples{
####################################################################
#a 8 points initial design  (a LHS in 2d)
design.fact <- maximinLHS(8,k=2)		#dimension = 2
design.fact <- data.frame(design.fact)
names(design.fact) <- c ( "x1","x2")
testfun <- camelback2			#our test function

#the response
response <- testfun(design.fact)

#the initial km model
model <- km(formula=~., design = design.fact, response = response, covtype="matern5_2")

#the integration points
intpoints <- expand.grid(seq(0,1,length=20), seq(0,1,length=20))

#the inversion threshold
T <- 0

#the number of iterations
iter <- 1


#inversion
res1  <- EGI(T=T,model=model,iter=iter,sampling.method="sur",
lower=c(0,0),upper=c(1,1),fun=testfun,discrete.X=intpoints)
res2  <- EGI(T=T,model=model,iter=iter,sampling.method="sur",
lower=c(0,0),upper=c(1,1),fun=testfun)
res3  <- EGI(T=T,model=model,iter=iter,sampling.method="timse",
lower=c(0,0),upper=c(1,1),fun=testfun,discrete.X=intpoints)	
res4  <- EGI(T=T,model=model,iter=iter,sampling.method="ranjan",
lower=c(0,0),upper=c(1,1),fun=testfun,method.param=0.5)	
res5  <- EGI(T=T,model=model,iter=iter,sampling.method="bichon",
lower=c(0,0),upper=c(1,1),fun=testfun,method.param=0.5)	

#plot of the function on a 50 x 50 grid
resolution <- 50
s <- seq(0,1,,resolution)
full_design <- expand.grid(s,s)
resp <- testfun(full_design)

par(mfrow=c(1,1))
contour(s,s,matrix(resp, nrow=resolution),20)

#initial points
points(design.fact, col="black", pch=20, lwd=4)

#new points
points(res1$par, col="red", pch=20, lwd=4)
text( res1$par[,1],res1$par[,2], seq(1,nrow(res1$par)), pos=3)

#new points
points(res2$par, col="blue", pch=20, lwd=4)
text( res2$par[,1],res2$par[,2], seq(1,nrow(res2$par)), pos=3)

#new points
points(res3$par, col="grey", pch=20, lwd=4)
text( res3$par[,1],res3$par[,2], seq(1,nrow(res3$par)), pos=3)

#new points
points(res4$par, col="green", pch=20, lwd=4)
text( res4$par[,1],res4$par[,2], seq(1,nrow(res4$par)), pos=3)

#new points
points(res5$par, col="purple", pch=20, lwd=4)
text( res5$par[,1],res5$par[,2], seq(1,nrow(res5$par)), pos=3)
##################################################################
#same example with noisy initial observations and noisy new 
#observations:

#library("lhs")
#a 8 points initial design  (a LHS in 2d)
design.fact <- maximinLHS(8,k=2)		#dimension = 2
design.fact <- data.frame(design.fact)
names(design.fact) <- c ( "x1","x2")
testfun <- camelback2			#our test function

#the response
response <- testfun(design.fact)+rnorm(n=8,mean=0,sd=0.25)

#the noise variance vector
noise.var <- rep(x=0.05,times=8)

#the initial km model
model <- km(formula=~1, design = design.fact, response = response,noise.var = noise.var )

#the integration points
intpoints <- expand.grid(seq(0,1,length=20), seq(0,1,length=20))

#the inversion threshold
T <- -1

#the number of iterations
iter <- 1

#the new observations variance
new.var <- 0.01

#EGI call
res  <- EGI(T=T,model=model,iter=iter,sampling.method="timse",
lower=c(0,0),upper=c(1,1),fun=testfun,discrete.X=intpoints,
new.noise.var = new.var)

#plot of the function on a 50 x 50 grid
resolution <- 50
s <- seq(0,1,,resolution)
full_design <- expand.grid(s,s)
resp <- testfun(full_design)

par(mfrow=c(1,1))
contour(s,s,matrix(resp, nrow=resolution),20)

#initial points
points(design.fact, col="black", pch=20, lwd=4)

#new points
points(res$par, col="red", pch=20, lwd=4)
text( res$par[,1],res$par[,2], seq(1,nrow(res$par)), pos=3)
##################################################################
}