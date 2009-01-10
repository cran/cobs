\name{DublinWind}
\alias{DublinWind}
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
      i.e., one of 1:366.}
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
  Space-Time Modelling with Long-Memory Dependence: Assessing
  Ireland's Wind Power Resource.
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
plot(speed ~ day, data = DublinWind)# not so nice; want time scale

## transform 'day' to correct "Date" object; and then plot
Dday <- seq(from = as.Date("1961-01-01"), by = 1,
            length = nrow(DublinWind))
plot(speed ~ Dday, data = DublinWind, type = "l",
     main = paste("DublinWind speed daily data, n=",
                  nrow(DublinWind)))

##--- ~ He & Ng "Example 3" %% much more is in ../tests/wind.R
co.o50 <-
 with(DublinWind, ## use nknots > (default) 6 :
  cobs(day, speed, knots.add = TRUE, constraint= "periodic", nknots = 10,
       tau = .5, method = "uniform"))
summary(co.o50)
lines(Dday, fitted(co.o50), col=2, lwd = 2)

## the periodic "smooth" function - in one period
plot(predict(co.o50), type = "o", cex = 0.2, col=2,
     xlab = "day", ylim = c(0,20))
points(speed ~ day, data = DublinWind, pch = ".")
}
\keyword{datasets}
