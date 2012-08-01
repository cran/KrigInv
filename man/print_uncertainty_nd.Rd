\name{print_uncertainty_nd}
\Rdversion{1.3}
\alias{print_uncertainty_nd}

\title{
Print a measure of uncertainty for functions with dimension d strictly higher than 2.
}
\description{
This function draws projections on various plans of a given measure of uncertainty. 
Possible measures are \code{"pn"} (probability of excursion) and measures specific to a sampling criterion: \code{"sur"}, \code{"timse"} and \code{"imse"}.
This function can be used to print relevant outputs after having used the function \code{\link{EGI}}.
}

\usage{
print_uncertainty_nd(model,T,type="pn",lower=NULL,upper=NULL,
  		resolution=20, nintegpoints=400,main="",
			cex.main=1,cex.lab=1,cex.contourlab=1,cex.axis=1,
			nlevels=10,levels=NULL,
			xdecal=3,ydecal=3, option="mean")
}

\arguments{
  \item{model}{
Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{T}{
Target value (a real number). The sampling algorithm and the underlying kriging model aim to find the points below (resp. over) T.
}
  \item{type}{
Type of uncertainty that the user wants to print. 
	Possible values are \code{"pn"} (probability of excursion), or 
	\code{"sur"}, \code{"imse"}, \code{"timse"} if we print a measure of uncertainty corresponding to one criterion.
}
  \item{lower}{
Vector containing the lower bounds of the input domain. If nothing is set we use a vector of 0.
}
  \item{upper}{
Vector containing the upper bounds of the input domain. If nothing is set we use a vector of 1.
}
  \item{resolution}{
Number of points to discretize a plan included in the domain. For the moment, we cannot use values higher than 40.
}
  \item{nintegpoints}{
to do
}
  \item{main}{
Title of the graph.
}
  \item{cex.main}{
Multiplicative factor for the size of the title.
}
  \item{cex.lab}{
Multiplicative factor for the size of titles of the axis.
}
  \item{cex.contourlab}{
Multiplicative factor for the size of labels of the contour plot.
}
  \item{cex.axis}{
Multiplicative factor for the size of the axis graduations.
}
  \item{nlevels}{
Integer corresponding to the number of levels of the contour plot.
}
  \item{levels}{
Array: one can directly set the levels of the contour plot.
}
  \item{xdecal}{
Optional position shifting of the titles of the x axis.
}
  \item{ydecal}{
Optional position shifting of the titles of the y axis.
}
  \item{option}{
Optional argument (a string). The 3 possible values are \code{"mean"} (default), \code{"max"} and \code{"min"}.
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
#print_uncertainty_nd

set.seed(8)
N <- 25 #number of observations
T <- -1 #threshold
testfun <- hartman3
#The hartman3 function is defined over the domain [0,1]^3. 

hartman3(runif(3))

lower <- rep(0,times=3)
upper <- rep(1,times=3)

#a 9 points initial design (LHS in 2 dimensions)
design <- data.frame(maximinLHS(N,k=3))
response <- apply(design,1,testfun)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

\dontrun{
print_uncertainty_nd(model=model,T=T,main="average probability of excursion",type="pn",
                    option="mean")

print_uncertainty_nd(model=model,T=T,main="maximum probability of excursion",type="pn",
                     option="max")

}
}
