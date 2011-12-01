\name{print_uncertainty}
\Rdversion{1.2}
\alias{print_uncertainty}

\title{
Prints a measure of uncertainty for a function of any dimension.
}
\description{
This function prints in the whole input domain the value of a given measure of uncertainty. 
Possible measures are \code{"pn"} (the probability of excursion) and measures specific to a sampling criterion: \code{"sur"}, \code{"timse"} and \code{"imse"}.
This function can be used to print relevant outputs after having used the function \code{\link{EGI}}.
}
\usage{
print_uncertainty(model, T, type = "pn", ...)
}

\arguments{
  \item{model}{
Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{T}{
Target value (scalar).
}
  \item{type}{
Type of uncertainty that the user wants to print. 
	Possible values are \code{"pn"} (probability of excursion), or 
	\code{"sur"}, \code{"imse"}, \code{"timse"} if we print a measure of uncertainty corresponding to one criterion.
}
  \item{\dots}{
Other arguments of the function \code{print_uncertainty_1d}, \code{print_uncertainty_2d} or \code{print_uncertainty_nd}.
}
}

\value{
the integrated uncertainty
}

\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing, pp.1-21, 2011, \url{http://arxiv.org/abs/1009.5177}
}

\author{ 
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

\seealso{\code{\link{EGI}}}

\examples{
#print_uncertainty

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

#par(mfrow=c(2,2)) #you could do many plots, but only one is run here
print_uncertainty(model=model,T=T,main="probability of excursion",type="pn")
#print_uncertainty(model=model,T=T,main="imse uncertainty",type="imse")
#print_uncertainty(model=model,T=T,main="timse uncertainty",type="timse")
#print_uncertainty(model=model,T=T,main="sur uncertainty",type="sur")
}
