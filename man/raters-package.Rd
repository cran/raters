\name{raters-package}
\alias{raters-package}
\alias{raters}
\title{
Inter-rater agreement among a set of raters

}
\description{
Computes a statistic as an index of inter-rater agreement among a set of raters.This procedure is based on a statistic not affected by Kappa paradoxes.
It is also possible to evaluate if the agreement is nil using the test argument.
The p value can be approximated using the Normal, Chi-squared distribution or 
using Monte Carlo algorithm. Fleiss' Kappa is also shown.

}
\details{
\tabular{ll}{
Package: \tab raters\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2014-03-20\cr
License: \tab GPL-2\cr
}

}
\author{
Daniele Giardiello, Piero Quatto and Stefano Vigliani

Maintainer: Daniele Giardiello <daniele.giardiello1@gmail.com>

}
\references{

Fleiss, J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bulletin \bold{76, 378-382}.

Falotico, R. Quatto, P. (2010). On avoiding paradoxes in assessing inter-rater agreement. Italian Journal of Applied Statistics \bold{22, 151-160.} 

}

\keyword{datasets}

\examples{
data(diagnostic)
concordance(diagnostic,test="Normal")
}