\name{globtemp}
\alias{globtemp}
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
str(globtemp)
## forget about time-series, just use numeric vectors:
year <- as.vector(time(globtemp))
temp <- as.vector(globtemp)

%% --- even much more is in ../tests/temp.R
%%                          ~~~~~~~~~~~~~~~
##---- Code for Figure 1a of He and Ng (1999) ----------

a50 <- cobs(year, temp, knots.add = TRUE, degree = 1, constraint = "increase")
summary(a50)
## As suggested in the warning message, we increase the number of knots to 9
a50 <- cobs(year, temp, nknots = 9, knots.add = TRUE, degree = 1,
            constraint = "increase")
summary(a50)
## Here, we use the same knots sequence chosen for the 50th percentile
a10 <- cobs(year, temp, nknots = length(a50$knots), knots = a50$knot,
            degree = 1, tau = 0.1, constraint = "increase")
summary(a10)
a90 <- cobs(year, temp, nknots = length(a50$knots), knots = a50$knot,
            degree = 1, tau = 0.9, constraint = "increase")
summary(a90)

which(hot.idx  <- temp >= a90$fit)
which(cold.idx <- temp <= a10$fit)
normal.idx <- !hot.idx & !cold.idx

plot(year, temp, type = "n", ylab = "Temperature (C)", ylim = c(-.7,.6))
lines(predict(a50, year, interval = "both"), col = 2)
lines(predict(a10, year, interval = "both"), col = 3)
lines(predict(a90, year, interval = "both"), col = 3)
points(year, temp, pch = c(1,3)[2 - normal.idx])

## label the "hot" and "cold" days
text(year[hot.idx], temp[hot.idx] + .03, labels = year[hot.idx])
text(year[cold.idx],temp[cold.idx]- .03, labels = year[cold.idx])
}
\keyword{datasets}
