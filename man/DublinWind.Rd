\name{DublinWind}
\alias{DublinWind}
\non_function{}
\title{Daily Wind Speeds in Dublin}
\usage{data(DublinWind)}
\description{
  The \code{DublinWind} data frame is basically the time series of daily
  average wind speeds from 1961 to 1978, measured in Dublin, Ireland.
  These are 6574 observations (18 full years among which four leap years).
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{speed}{numeric vector of average daily wind speed in knots}
    \item{day}{an integer vector giving the day number of the year,
      i.e., 1:366.}
  }
}
\details{
  The periodic pattern along the 18 years measured and the
  autocorrelation are to be taken into account for analysis, see the
  references.  This is Example 3 of the COBS paper.
}
\source{
  From shar file available from \url{http://www.cba.nau.edu/pin-ng/cobs.html}

  Also available from \url{.......STATLIB..........}
}
\references{
  Haslett, J. and Raftery, A. (1989)
  Space-Time Modelling with Long-Memory Dependence: {A}ssessing
  {I}reland's Wind Power Resource.
  .......... % FIXME -- no journal !
  \bold{38}, 1--21.

  COBS: Qualitatively Constrained Smoothing via Linear Programming;
  \emph{Computational Statistics} \bold{14}, 315--337.

  He, X. and Ng, P. (1999)
  COBS: Qualitatively Constrained Smoothing via Linear Programming;
  \emph{Computational Statistics} \bold{14}, 315--337.
}
\examples{
data(DublinWind)
str(DublinWind)
plot(speed ~ day, data = DublinWind, type = "l")
attach(DublinWind)
plot(wSpeed <- ts(speed))
%%## FIXME: use correct time axis!--> need leap year correctness
%%##wSpeed <- ts(speed, start = 1961, end = 1978.999, periodicity= (365 + 1/4))

%% currently He & Ng "Example 3" is --> ../tests/wind.R
}
\keyword{datasets}
