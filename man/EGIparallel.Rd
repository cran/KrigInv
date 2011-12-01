\name{EGIparallel}
\Rdversion{1.2}
\alias{EGIparallel}

\title{
Efficient Global Inversion: Parallel version.
}
\description{
Sequential sampling based on the optimization of a kriging-based criterion, with model update after each evaluation. 
Similar to EGI, except that samples can be chosen in batches instead of one at a time. Two criteria are available for now.
}
\usage{
EGIparallel(T, model, method = NULL, method.param=NULL, 
fun, iter, batchsize = 1, 
lower, upper, new.noise.var = 0, 
optimcontrol = NULL, kmcontrol = NULL, integcontrol = NULL, ...)
}

\arguments{
  \item{T}{
Target value (scalar). The sampling algorithm and the underlying kriging model at finding the points for which the output is close to T.
}
  \item{model}{
A Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{method}{
Criterion used for choosing observations. The two sampling methods available in parallel version are \code{"sur"} and \code{"timse"}.
}
  \item{method.param}{
Optional tolerance value (a real number) for method \code{"timse"}.
}
  \item{batchsize}{
Number of points to sample simultaneously. The sampling criterion will return \code{batchsize} points at a time for sampling.
}
  \item{new.noise.var}{
Optional scalar value of the noise variance of the new observations.
}
  \item{fun}{
Objective function.
}
  \item{iter}{
Number of iterations.
}
  \item{lower}{
Vector containing the lower bounds of the variables to be optimized over.
}
  \item{upper}{
Vector containing the upper bounds of the variables to be optimized over.
}
  \item{optimcontrol}{
Optional list of control parameters for the optimization of the sampling criterion. The field \code{method} defines which optimization method is used: it can be either \code{"genoud"} (default) for an optimisation using the genoud algorithm, or \code{"discrete"} for an optimisation over a specified discrete set. 
If the field \code{method} is set to \code{"genoud"}, one can set some parameters of this algorithm: 
\code{pop.size}  (default : 50*d),  \code{max.generations} (10*d), \code{wait.generations} (2),  \code{BFGSburnin} (2) and the mutations \code{P1}, \code{P2}, up to \code{P9} (see \code{\link[rgenoud]{genoud}}). Numbers into brackets are the default values.
If the field \code{method} is set to \code{"discrete"}, one can set the field \code{optim.points}: p * d matrix corresponding to the p points where the criterion will be evaluated. If nothing is specified, 100*d points are chosen randomly. 
Finally, one can control the field {optim.option} in order to decide how to optimize the sampling criterion. 
If \code{optim.option} is set to 2 (default), batchsize sequential optimizations in dimension d are performed to find the optimum.
If \code{optim.option} is set to 1, only one optimization in dimension batchsize*d is performed. This option is only available with \code{"genoud"}. This option might provide more global and accurate solutions, but is a lot more expensive.
}
  \item{kmcontrol}{
Optional list representing the control variables for the re-estimation of the kriging model once new points are sampled. 
The items are the same as in \code{\link[DiceKriging]{km}}.
}
  \item{integcontrol}{
Optional list specifying the procedure to build the integration points and weights. Many options are possible.
A) If nothing is specified, 100*d points are chosen using the Sobol sequence. 
B) One can directly set the field \code{integration.points} (a p * d matrix) for prespecified integration points. In this case these integration points and the corresponding vector \code{integration.weights} will be used for all the iterations of the algorithm. 
C) If the field \code{integration.points} is not set then the integration points are renewed at each iteration. 
In that case one can control the number of integration points \code{n.points} (default: 100*d) and a specific distribution \code{distrib}. Possible values for \code{distrib} are: 
\code{"sobol"}, \code{"MC"}, \code{"timse"}, \code{"imse"}, \code{"sur"} and \code{"jn"} (default: \code{"sobol"}). 
C.1) The choice \code{"sobol"} corresponds to integration points chosen with the Sobol sequence in dimension d (uniform weight). 
C.2) The choice \code{"MC"} corresponds to points chosen randomly, uniformly on the domain. 
C.3) The choices \code{"timse"}, \code{"imse"}, \code{"sur"} and \code{"jn"} correspond to importance sampling distributions (unequal weights). It is strongly recommended to use the importance sampling distribution corresponding to the chosen sampling criterion. 
When important sampling procedures are chosen, \code{n.points} points are chosen using importance sampling among a discrete set of \code{n.candidates} points (default: \code{n.points*10}) which are distributed according to a distribution \code{init.distrib} (default: \code{"sobol"}). Possible values for \code{init.distrib} are the space filling distributions \code{"sobol"} and \code{"MC"} or an user defined distribution \code{"spec"}. 
The \code{"sobol"} and \code{"MC"} choices correspond to quasi random and random points in the domain. 
If the \code{"spec"} value is chosen the user must fill in manually the field \code{init.distrib.spec} to specify himself a n.candidates * d matrix of points in dimension d.
}
  \item{\dots}{
Other arguments of the target function \code{fun}.
}
}
\details{
The function used to build the integration points and weights (based on the options specified in \code{integcontrol}) is the function \code{\link{integration_design}}
}
\value{A list with components:
\item{par}{The added observations ((ite*batchsize) * d matrix)}
\item{value}{The value of the function \code{fun} at the added observations (vector of size ite*batchsize)}
\item{nsteps}{The number of added observations (=ite*batchsize).}
\item{lastmodel}{The current (last) kriging model of \code{\link[DiceKriging]{km}} class.}
\item{lastvalue}{The value of the criterion at the last added batch of points.}
\item{allvalues}{If an optimization on a discrete set of points is chosen, the value of the criterion at all these points, for the last iteration, for the last point of the batch.}
}

\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing, pp.1-21, 2011, \url{http://arxiv.org/abs/1009.5177}

Chevalier C., Bect J., Ginsbourger D., Picheny V., Richet Y., Vazquez E. (2011), \emph{Fast kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}  ,\url{http://hal.archives-ouvertes.fr/docs/00/64/11/08/PDF/FastKrigingInversion.pdf}
}

\author{
Clement Chevalier (IMSV, Switzerland, and IRSN, France)

Victor Picheny (CERFACS, Toulouse, France) 
}

\seealso{\code{\link{EGI}}, \code{\link{max_sur_parallel}}, \code{\link{sur_optim_parallel}}}

\examples{
#EGIparallel

set.seed(8)
N <- 9 #number of observations
T <- 80 #threshold
testfun <- branin
lower <- c(0,0)
upper <- c(1,1)

#a 9 points initial design (LHS in 2 dimensions)
design <- data.frame(maximinLHS(N,k=2))
response <- testfun(design)

#km object with matern3_2 covariance (params estimated by ML from the observations)
model <- km(formula=~., design = design, response = response,covtype="matern3_2")

optimcontrol <- list(method="genoud",pop.size=50)
integcontrol <- list(distrib="sur",n.points=50)
iter <- 1
batchsize <- 6

\dontrun{
obj <- EGIparallel(T=T,model=model,method="sur",batchsize=batchsize,
                    fun=testfun,iter=iter,lower=lower,upper=upper,
                    optimcontrol=optimcontrol,integcontrol=integcontrol)

par(mfrow=c(1,2))
print_uncertainty_2d(model=model,T=T,main="probability of excursion",
type="pn",new.points=0,cex.points=2)

print_uncertainty_2d(model=obj$lastmodel,T=T,
main="probability of excursion, parallel sur sampling",
type="pn",new.points=iter*batchsize,col.points.end="red",cex.points=2)
}

##############
#same example with noisy initial observations and noisy new observations
branin.noise <- function(x) return(branin(x)+rnorm(n=1,sd=30))

set.seed(8)
N <- 9;T <- 80
testfun <- branin.noise
lower <- c(0,0);upper <- c(1,1)

design <- data.frame(maximinLHS(N,k=2))
response.noise <- apply(design,1,testfun)
response.noise - response

model.noise <- km(formula=~., design = design, response = response.noise,
covtype="matern3_2",noise.var=rep(30*30,times=N))

optimcontrol <- list(method="genoud",pop.size=50)
integcontrol <- list(distrib="sur",n.points=50)
iter <- 1
batchsize <- 6

\dontrun{
obj <- EGIparallel(T=T,model=model.noise,method="sur",batchsize=batchsize,
                   fun=testfun,iter=iter,lower=lower,upper=upper,
                   optimcontrol=optimcontrol,integcontrol=integcontrol,
                   new.noise.var=30*30)


par(mfrow=c(1,2))
print_uncertainty_2d(model=model.noise,T=T,
main="probability of excursion, noisy obs.",
type="pn",new.points=0,cex.points=2)

print_uncertainty_2d(model=obj$lastmodel,T=T,
main="probability of excursion, parallel sur sampling, noisy obs.",
type="pn",new.points=iter*batchsize,col.points.end="red",cex.points=2)
}
}
