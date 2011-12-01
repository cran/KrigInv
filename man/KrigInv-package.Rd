\name{KrigInv-package}
\alias{KrigInv}
\docType{package}

\title{ Kriging-based inversion of deterministic and stochastic computer codes }

\description{
Sequential algorithms based on Kriging for computer experiments, meant to explore the subset of input parameters corresponding to a prescribed level of the output.
}

\details{
\tabular{ll}{
Package: \tab KrigInv\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2011-12-01\cr
License: \tab GPL version 3\cr
LazyLoad: \tab yes\cr
}
}

\author{
C. Chevalier (IMSV, University of Bern, Switzerland and IRSN)

V. Picheny (Ecole Centrale Paris)

D. Ginsbourger (IMSV, University of Bern, Switzerland)

with contributions from Yann Richet (IRSN)

Maintainer: C. Chevalier (clement.chevalier@stat.unibe.ch)
}

\note{
A first prototype of this package was originally developed by D. Ginsbourger in the frame of a collaboration with IRSN (Institut de Radioprotection et de Surete Nucleaire), acting through Yann Richet.

The three main authors thank IRSN for sponsoring open source research, and allowing them to spread the present package and publish it on CRAN.

They also would like to warmly thank Yann Richet for numerous discussions concerning this package, and more!

Package \code{rgenoud} >=5.3.3. is recommended.

Important function (key to all proposed methods):
\tabular{ll}{
\code{EGI} \tab Sequential Kriging-based inversion \cr

}

}


\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing, pp.1-21, 2011, \url{http://arxiv.org/abs/1009.5177}

Chevalier C., Bect J., Ginsbourger D., Picheny V., Richet Y., Vazquez E. (2011), \emph{Fast kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}  ,\url{http://hal.archives-ouvertes.fr/docs/00/64/11/08/PDF/FastKrigingInversion.pdf}

Picheny, V., Ginsbourger, D., Roustant, O., Haftka, R.T., Adaptive designs of experiments for accurate approximation of a target region, J. Mech. Des. - July 2010 - Volume 132, Issue 7, http://dx.doi.org/10.1115/1.4001873

Picheny V., Improving accuracy and compensating for uncertainty in surrogate modeling, Ph.D. thesis, University of Florida and Ecole Nationale Superieure des Mines de Saint-Etienne

Vazquez, E., Bect, J.: A sequential Bayesian algorithm to estimate a probability of failure. In: Proceedings of the 15th IFAC Symposium on System Identification, (SYSID 2009), Saint-Malo, France (2009)

Bichon, B.J., Eldred, M.S., Swiler, L.P., Mahadevan, S., McFarland, J.M.: Efficient global reliability analysis for nonlinear implicit performance functions. AIAA Journal 46 (10), 2459-2468 (2008)

Ranjan, P., Bingham, D., Michailidis, G.: Sequential experiment design for contour estimation from complex computer codes. Technometrics 50(4), 527-541 (2008)

Rasmussen C.E., Williams C.K.I. (2006), \emph{Gaussian Processes for Machine Learning}, the MIT Press, \url{www.GaussianProcess.org/gpml}
}

\keyword{package}