\name{update_km}
\Rdversion{1.3}
\alias{update_km}

\title{Update of a kriging model}

\description{
Update a \code{\link[DiceKriging]{km}} object when one or many new observations are added. Many, but not all, fields of the \code{\link[DiceKriging]{km}} object need to be recalculated when new observations are added. 
It is also possible to modify the k last (existing) observations.
}

\usage{
update_km(model, NewX, NewY, NewX_AllreadyExist = FALSE, 
CovReEstimate = FALSE, new.noise.var = NULL, 
kmcontrol = NULL, F.newdata = NULL)
}

\arguments{
	\item{model}{Kriging model of \code{\link[DiceKriging]{km}} class.}

	\item{NewX}{Matrix with \code{model@d} columns and r rows corresponding to the r locations of the observations to be updated. These locations can be new locations or existing ones.}

	\item{NewY}{Matrix with one column and r rows corresponding to the r responses at the r locations NewX.}

	\item{NewX_AllreadyExist}{Boolean: indicate whether the locations NewX are all news or not.}

	\item{CovReEstimate}{Boolean to decide whether the covariance parameters of the \code{\link[DiceKriging]{km}} object should be re-estimated.}

	\item{new.noise.var}{Vector containing the noise variance at each new observations.}
	  
	\item{kmcontrol}{Optional list representing the control variables for the re-estimation of the kriging model once new points are sampled. The items are the same as in \code{\link[DiceKriging]{km}}}
	
	\item{F.newdata}{Optional matrix containing the value of the trend at the new locations. Setting this argument avoids a call to an expensive function.}
}

\value{
Updated km object
}

\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing, pp.1-21, 2011, \url{http://arxiv.org/abs/1009.5177}

Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2011), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set} ,\url{http://hal.archives-ouvertes.fr/hal-00641108/}
}
\author{ 
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

\seealso{\code{\link[DiceKriging]{km}}}

\examples{
#update_km

set.seed(8)
N <- 9 #number of observations
testfun <- branin

#a 9 points initial design
design <- data.frame( matrix(runif(2*N),ncol=2) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")
model@covariance

new.x <- matrix(c(0.4,0.5),ncol=2)#the point that we are going to add in the km object
new.y <- testfun(new.x)
new.model <- update_km(model=model,NewX=new.x,NewY=new.y,CovReEstimate=TRUE)
new.model@covariance
}
