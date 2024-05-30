\name{concordance}
\alias{concordance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Inter-rater agreement among a set of raters for nominal data
%%  ~~function to do ... ~~
}
\description{
Computes a statistic as an index of inter-rater agreement among a set of raters in case of nominal data.  This procedure is based on a statistic not affected by paradoxes of Kappa.
It is also possible to get the confidence interval at level alpha using the percentile Bootstrap and to evaluate if the agreement is nil using the test argument.
The p value can be approximated using the Normal, Chi squared distribution or using Monte Carlo algorithm.
Normal approximation and Monte Carlo procedure can be calculated
even though the number of observers is not the same for each evaluated subject.
Fleiss Kappa is also shown and its confidence interval, standard error and pvalue using Normal approximation are
available when the number of observes is the same for each classified subject and the test argument is specified.
The functions \code{wlin.conc} and \code{wquad.conc} can be used in case of ordinal data using linear or quadratic weight matrix, respectively.


%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
concordance(db, test = "Default", B = 1000, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{db}{n*c matrix or data frame, n subjects c categories. The numbers inside the matrix or data frame indicate how many raters chose a specific category for a given subject.  A sum of row indicates the total number of raters who evaluated a given subject. }
  \item{test}{Statistical test to evaluate if the raters make random assignment regardless of the characteristic of each subject. Under null hypothesis, it corresponds to a high percentage of assignment errors. Thus, the expected agreement is weak.
  Normal approximation is advisable when the number of subject is pretty large while a Chi square approximation when the number of raters is large.
  Monte Carlo test is useful for small samples even though the higher number of simulations, the more time the procedure takes.
  If test is not mentioned, no test will be computed.
}
  \item{B}{ Number of iterations for the percentile Bootstrap and for Monte Carlo test.}
  \item{alpha}{ Level of significance for Bootstrap confidence interval}

}


\value{
  A list containing the following components:
  \item{$Fleiss}{A list with Kappa of Fleiss index. When the number of raters is the same for every evaluated subject, the standard deviation, the Z Wald test and the p value are also shown.}

  \item{$Statistic}{A list with the index of inter-rater agreement not affected by Kappa paradoxes and the percentile Bootstrap confidence interval. If the test argument is specified the p value is also shown. }



}
\references{

Fleiss, J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bulletin \bold{76, 378-382}

Falotico, R. Quatto, P. (2010). On avoiding paradoxes in assessing inter-rater agreement. Italian Journal of Applied Statistics \bold{22, 151-160}

}
\author{
Piero Quatto \email{piero.quatto@unimib.it},
Daniele Giardiello \email{daniele.giardiello1@gmail.com},
Stefano Vigliani
}

\examples{
data(diagnostic)
concordance(diagnostic, test = "Chisq")
concordance(diagnostic, test = "Normal")
concordance(diagnostic, test = "MC", B = 100)
}

\keyword{datasets}
