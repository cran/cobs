\name{globtemp}
\alias{globtemp}
\non_function{}
\title{Annual Average Global Surface Temperature}
\usage{data(globtemp)}
\description{
  Time Series of length 113 of annual average global surface temperature
  deviations from 1880 to 1992.
%% deviations from what?  The mean is not 0, but -0.084; median = -0.06
}
\details{
  This is Example 1 of the COBS paper, where the hypothesis of a
  monotonely increasing trend is considered; Koenker and Schorfheide
  (1994) consider modeling the autocorrelations.
}
\source{
  \file{temp.data} in file \file{cobs.shar} available from
  \url{http://www.cba.nau.edu/pin-ng/cobs.html}.
}
\references{
  He, X. and Ng, P. (1999)
  COBS: Qualitatively Constrained Smoothing via Linear Programming;
  \emph{Computational Statistics} \bold{14}, 315--337.

  Koenker, R. and Schorfheide F. (1994)
  Quantile Spline Models for Global Temperature Change;
  \emph{Climate Change} \bold{28}, 395--404.
}
\examples{
data(globtemp)
plot(globtemp, main = "Annual Global Temperature Deviations")

%% currently He & Ng "Example 1" is --> ../tests/temp.R
}
\keyword{datasets}
