\name{USArmyRoofs}
\alias{USArmyRoofs}
\title{Roof Quality in US Army Bases}
\usage{data(USArmyRoofs)}
\description{
  The \code{USArmyRoofs} data frame has 153 observations of roof sections of US
  Army bases and 2 columns, \code{age} and \code{fci}.  This is Example 2
  of He & Ng (1999).
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{age}{numeric vector giving the roof's age in years.}
    \item{fci}{numeric, giving the FCI, the flash condition index, i.e.,
      the percentage of flashing which is in good condition.}
  }
}
\source{
  From shar file available from \url{https://www2.nau.edu/PinNg/cobs.html}
}
\references{
  He, X. and Ng, P. (1999)
  COBS: Qualitatively Constrained Smoothing via Linear Programming;
  \emph{Computational Statistics} \bold{14}, 315--337.
}
\examples{
data(USArmyRoofs)
plot(fci ~ age, data = USArmyRoofs, main = "US Army Roofs data")

%% currently He & Ng "Example 2" is --> ../tests/roof.R
}
\keyword{datasets}
