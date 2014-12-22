\name{integration_design}
\Rdversion{1.3}
\alias{integration_design}

\title{
Construction of a sample of integration points and weights
}
\description{
Generic function to build integration points for some sampling criterion. 
Available important sampling schemes are \code{"sur"}, \code{"jn"}, \code{"timse"}, \code{"vorob"} and \code{"imse"}. 
Each of them corresponds to a sampling criterion. 
}
\usage{
integration_design(integcontrol = NULL, d = NULL, 
lower, upper, model = NULL, T = NULL,min.prob=0.001)
}

\arguments{
  \item{integcontrol}{
Optional list specifying the procedure to build the integration points and weights, relevant only for the sampling criteria based on numerical integration: 
(\code{"imse"}, \code{"timse"}, \code{"sur"} or \code{"jn"}). Many options are possible.
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
  \item{d}{
The dimension of the input set. If not provided d is set equal to the length of \code{lower}.
}
  \item{lower}{
Vector containing the lower bounds of the design space.
}
  \item{upper}{
Vector containing the upper bounds of the design space.
}
  \item{model}{
A Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{T}{
Target value (scalar). The sampling algorithm and the underlying kriging model aim at finding the points for which the output is close to T.
}
  \item{min.prob}{
  This argument applies only when importance sampling distributions are chosen (to compute integral criteria like \code{"sur"} or \code{"timse"}). For numerical reasons we give a minimum probability for a point to belong to the importance sample. This avoids probabilities equal to zero and importance sampling weights equal to infinity. In an importance sample of M points, the maximum weight becomes \code{1/min.prob * 1/M}.
  }
}
\details{
The important sampling aims at improving the accuracy of the calculation of numerical integrals present in these criteria.
}
\value{A list with components:
\item{integration.points}{ p x d matrix of p points used for the numerical calculation of integrals}
\item{integration.weights}{ a vector of size p corresponding to the weight of each point. If all the points are equally weighted, \code{integration.weights} is set to \code{NULL}}
\item{alpha}{ if the \code{"vorob"} important sampling schemes is chosen, the function also returns a scalar, alpha, being the calculated Vorob'ev threshold}
}
\references{
Chevalier C., Picheny V., Ginsbourger D. (2012), \emph{The KrigInv package: An efficient and user-friendly R implementation of Kriging-based inversion algorithms} ,
\url{http://hal.archives-ouvertes.fr/hal-00713537/}

Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2011), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set} ,\url{http://hal.archives-ouvertes.fr/hal-00641108/}

Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing, pp.1-21, 2011, \url{http://arxiv.org/abs/1009.5177}
}

\author{
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

\seealso{\code{\link{EGI}}, \code{\link{max_timse}}, \code{\link{max_sur}}}

\examples{
#integration_design

\dontrun{
set.seed(8)
#when nothing is specified: integration points 
#are chosen with the sobol sequence
integ.param <- integration_design(lower=c(0,0),upper=c(1,1))
plot(integ.param$integration.points)
}

#an example with pure random integration points
integcontrol <- list(distrib="MC",n.points=50)
integ.param <- integration_design(integcontrol=integcontrol,
		lower=c(0,0),upper=c(1,1))
plot(integ.param$integration.points)

#an example with important sampling distributions
#these distributions are used to compute integral criterion like
#"sur","timse" or "imse"

#for these, we need a kriging model
N <- 14;testfun <- branin; T <- 80
lower <- c(0,0);upper <- c(1,1)
design <- data.frame( matrix(runif(2*N),ncol=2) )
response <- testfun(design)
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

integcontrol <- list(distrib="timse",n.points=200,n.candidates=5000,init.distrib="MC")
integ.param <- integration_design(integcontrol=integcontrol,
		lower=c(0,0),upper=c(1,1), model=model,T=T)

print_uncertainty_2d(model=model,T=T,type="timse",
col.points.init="red",cex.points=2,
main="timse uncertainty and one sample of integration points")

points(integ.param$integration.points,pch=17,cex=1)

}
