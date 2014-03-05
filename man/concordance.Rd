\name{concordance}
\alias{concordance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Inter rater agreement among a set of raters
%%  ~~function to do ... ~~
}
\description{
Computes a statistic as an index of inter rater agreement among a set of raters.  This procedure is based on a statistic not affected by \dQuote{Kappa} paradoxes.  
It is also possible to evaluate if the agreement is nil using the test argument.
The p value can be approximated using the Normal, Chi squared distribution or 
using Monte Carlo algorithm.
Normal approximation and Monte Carlo procedure can be calculated
even though the number of observers is not the same for each evaluated subject.
\sQuote{Fleiss} Kappa is also shown and its confidence interval is
available when the number of observes is the same for each classified subject.


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
  \item{B}{ Number of iterations for Monte Carlo test.}   
  \item{alpha}{ Level of significance.}

}


\value{
  A list containing the following components:
  \item{$Fleiss}{ \sQuote{Fleiss} Kappa index. }
  \item{$S}{ Computed index of inter rater agreement not affected by Kappa paradoxes. }
  \item{$Lower}{ Lower limit of \sQuote{Fleiss} Kappa confidence interval. It is computed only when the number of raters is the same for all subjects. }
  \item{$Upper}{ Upper limit of \sQuote{Fleiss} Kappa confidence interval. It is computed only when the number of raters is the same for all subjects. }
  \item{$pvalue}{ P value for the statistical test. }
  
  
}
\references{
  
Fleiss, J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bulletin \bold{76, 378-382}

Falotico, R. Quatto, P. (2010). On avoiding paradoxes in assessing inter-rater agreement. Italian Journal of Applied Statistics \bold{22, 151-160}
}
\author{
Piero Quatto \email{piero.quatto@unimib.it},
Daniele Giardiello \email{daniele.giardiello1@gmail.com}
}

\examples{
data(diagnostic)
concordance(diagnostic,test="Chisq")
concordance(diagnostic,test="Normal")
concordance(diagnostic,test="MC",B=100)
}

\keyword{datasets}
