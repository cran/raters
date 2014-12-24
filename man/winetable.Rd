\name{winetable}
\alias{winetable}
\docType{data}
\title{
Sensory wine evaluation
}
\description{ The data wine in ordinal package represent a factorial experiment on factors determining the bitterness of wine
with 1 as \dQuote{Least bitter}and 5 as \dQuote{Most bitter}.
In this case, we supposed that eight different bottles of wine were evaluated by nine judges according to the bitterness scale described above.
It is possible to get this data using the dataframe wine included in ordinal package using \code{table(wine$bottle,wine$rating)}
}
\usage{data(winetable)}
\format{
  A data frame with 8 observations on 5 variables.
}

\references{
%%  ~~ possibly secondary sources and usages ~~
Randall J.H.(1989) The Analysis of Sensory Data by Generalized Linear Model. Biometrical Journal
\bold{vol 31,issue 7, pp 781-793}
}

\keyword{datasets}
