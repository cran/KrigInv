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
Version: \tab 1.1\cr
Date: \tab 2010-12-27\cr
License: \tab GPL version 3\cr
LazyLoad: \tab yes\cr
}
}

\author{
V. Picheny (Ecole Centrale Paris)

D. Ginsbourger (IMSV, University of Bern, Switzerland)

with contributions from Clement Chevalier (IMSV and IRSN) and Yann Richet (IRSN)

Maintainer: D. Ginsbourger (david.ginsbourger@stat.unibe.ch)
}

\note{
A first prototype of this package was originally developed by D. Ginsbourger in the frame of a collaboration with IRSN (Institut de Radioprotection et de Surete Nucleaire), acting through Yann Richet.

The two main authors thank IRSN for sponsoring open source research, and allowing them to spread the present package and publish it on CRAN.

They also would like to warmly thank Clement Chevalier for his precious contributions and advices in programming, as well as Yann Richet for numerous discussions concerning this package, and more!

Package \code{rgenoud} >=5.3.3. is recommended.

Important function (key to all proposed methods):
\tabular{ll}{
\code{EGI} \tab Sequential Kriging-based inversion \cr

%\code{qEI} \tab q-points noise-free EI criterion\cr
}

}


\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, accepted with minor revisions to the Journal of Statistics and Computing, \url{http://arxiv.org/abs/1009.5177}

Picheny, V., Ginsbourger, D., Roustant, O., Haftka, R.T., Adaptive designs of experiments for accurate approximation of a target region, J. Mech. Des. - July 2010 - Volume 132, Issue 7, \url{http://dx.doi.org/10.1115/1.4001873}

Vazquez, E., Bect, J.: A sequential Bayesian algorithm to estimate a probability of failure. In: Proceedings of the 15th IFAC Symposium on System Identification, (SYSID 2009), Saint-Malo, France (2009)

Bichon, B.J., Eldred, M.S., Swiler, L.P., Mahadevan, S., McFarland, J.M.: Efficient global reliability analysis for nonlinear implicit performance functions. AIAA Journal 46 (10), 2459-2468 (2008)

Ranjan, P., Bingham, D., Michailidis, G.: Sequential experiment design for contour estimation from complex computer codes. Technometrics 50(4), 527-541 (2008)

Rasmussen C.E., Williams C.K.I. (2006), \emph{Gaussian Processes for Machine Learning}, the MIT Press, \url{www.GaussianProcess.org/gpml}
}

%\examples{}

\keyword{package}