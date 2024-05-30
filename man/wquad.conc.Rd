\name{wquad.conc}
\alias{wquad.conc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Inter-rater agreement among a set of raters for ordinal data using quadratic weights
%%  ~~function to do ... ~~
}
\description{
Computes a statistic as an index of inter-rater agreement among a set of raters in case of ordinal data using quadratic weights.
The matrix of quadratic weights is defined inside the function.
This procedure is based on a statistic not affected by Kappa paradoxes.
It is also possible to get the confidence interval at level alpha using the percentile Bootstrap and to evaluate if the agreement is nil using the Monte Carlo algorithm.
Fleiss' Kappa cannot be used in case of ordinal data.
It is advisable to use \code{set.seed} to get the same replications for Bootstrap confidence limits and Montecarlo test.


%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wquad.conc(db, test = "Default", B = 1000, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{db}{n*c matrix or data frame, n subjects c categories. The numbers inside the matrix or data frame indicate how many raters chose a specific category for a given subject.  A sum of row indicates the total number of raters who evaluated a given subject. In case of ordinal data, the c categories can be sorted according to a specific scale.}
  \item{test}{Statistical test to evaluate if the raters make random assignment regardless of the characteristic of each subject. Under null hypothesis, it corresponds to a high percentage of assignment errors. Thus, the expected agreement is weak. If this argument is not specified the p value are not being computed.
}
  \item{B}{ Number of iterations for the percentile Bootstrap and for Monte Carlo test.}
  \item{alpha}{ Level of significance for Bootstrap confidence interval and for Monte Carlo algorithm if it is specified}

}


\value{
  A list containing the following components:

  \item{$Statistic}{A list with the index of inter-rater agreement not affected by Kappa paradoxes for ordinal data and the percentile Bootstrap confidence interval. If the test argument is specified the p value is also shown. }



}
\references{

Fleiss, J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bulletin \bold{76, 378-382}

Falotico, R. Quatto, P. (2010). On avoiding paradoxes in assessing inter-rater agreement. Italian Journal of Applied Statistics \bold{22, 151-160}

Marasini, D. Quatto, P. Ripamonti, E. (2014).  Assessing the inter-rater agreement for ordinal data through weighted indexes. Statistical methods in medical research.

}
\author{
Piero Quatto \email{piero.quatto@unimib.it},
Daniele Giardiello \email{daniele.giardiello1@gmail.com},
Stefano Vigliani
}

\examples{
data(uterine)
set.seed(12345)
wquad.conc(uterine, test = "MC", B = 25)
}

\keyword{datasets}
