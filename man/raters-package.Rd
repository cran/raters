\name{raters-package}
\alias{raters-package}
\alias{raters}
\title{
A Modification of Fleiss' Kappa in case of Nominal and Ordinal Variables
}
\description{
Computes a statistic as an index of inter-rater agreement among a set of raters in case of nominal or ordinal data.This procedure is based on a statistic not affected by Kappa paradoxes.
In case of ordinal data, the weighted versions of the statistic has been developed using a matrix of linear or quadratic weights. 
The percentile Boostrap confidence interval is computed and the test argument allows to perform if the agreement is nil.
The p value can be approximated using the Normal, Chi-squared distribution or using Monte Carlo algorithm in case of nominal data. Otherwise, the  approximation and the Monte Carlo algorithm is computed.
Fleiss' Kappa index is also shown in case of nominal data.
In a nutshell, the function \code{concordance} can be used
in case of nominal scale while the functions \code{wlin.conc} and \code{wquad.conc} can be used in case of ordinal 
data using linear or quadratic weights, respectively.


}
\details{
\tabular{ll}{
Package: \tab raters\cr
Type: \tab Package\cr
Version: \tab 2.0.1\cr
Date: \tab 2014-12-29\cr
License: \tab GPL(>=2)\cr
}

}
\author{
Daniele Giardiello, Piero Quatto, Enrico Ripamonti and Stefano Vigliani

Maintainer: Daniele Giardiello <daniele.giardiello1@gmail.com>

}
\references{

Fleiss, J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bulletin \bold{76, 378-382}.

Falotico, R. Quatto, P. (2010). On avoiding paradoxes in assessing inter-rater agreement. Italian Journal of Applied Statistics \bold{22, 151-160.} 

Falotico, R., Quatto, P. (2014). Fleiss' kappa statistic without paradoxes. Quality & Quantity, \bold{1-8}.

Marasini, D. Quatto, P. Ripamonti, E. (2014).  Assessing the inter-rater agreement for ordinal data through weighted indexes. Statistical methods in medical research.
}

\keyword{datasets}

\examples{
# Nominal data
data(diagnostic)
concordance(diagnostic,test="Normal")

# Ordinal data with linear weights
data(winetable)
set.seed(12345)
wlin.conc(winetable,test="MC")

# Ordinal data with quadratic weights
data(winetable)
set.seed(12345)
wquad.conc(winetable,test="MC")
}