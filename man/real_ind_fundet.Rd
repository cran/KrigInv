\name{real_ind_fundet}
\Rdversion{1.1}
\alias{real_ind_fundet}

\title{ Indicator function for the condition: fundet lies in the band ]T-epsilon, T+epsilon[ }

\description{
Binary function I[T-epsilon < fundet(x) < T+epsilon] indicating when the test function "fundet" lies in a prescribed band}

\usage{
real_ind_fundet(x,T,epsilon)
}

\arguments{
  \item{x}{1-dimensional input}
  \item{T}{ Target value (a real number)}
  \item{epsilon}{Tolerance value (a real number)}
}
%\details{}
\value{0 or 1 depending on the validity of the condition T-epsilon < fundet(x) < T+epsilon}
%\references{}
\author{ David Ginsbourger }
%\note{}

%\seealso{}

\examples{
level <- 0.5
epsilon <- 0.1
t <- seq(0,1,by=0.01)
plot(t, real_ind_fundet(t,T=level,epsilon=epsilon), type="l", col="red", lwd=2)
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory. 
\keyword{internal}
