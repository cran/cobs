\name{interpSplineCon}
\alias{interpSplineCon}
\alias{isIsplineCon}
\title{(Cubic) Interpolation Spline from "conreg"}
\description{
  From a \code{"\link{conreg}"} object representing a \emph{linear}
  spline,
  \describe{
    \item{\code{interpSplineCon()}}{produces the corresponding (cubic)
      spline (via package \pkg{splines}' \code{\link{interpSpline}()})
      by interpolating at the knots, thus \dQuote{smoothing the kinks}.}
    \item{\code{isIsplineCon()}}{determines if the spline fulfills the
      same convexity / concavity constraints as the underlying
      \code{"\link{conreg}"} object.}
  }
}
\usage{
interpSplineCon(object, \dots)
isIsplineCon(object, isp, \dots)
}
\arguments{
  \item{object}{an \R object as resulting from \code{\link{conreg}()}.}
  \item{isp}{optionally, the result of \code{interpSplineCon(object,
      \dots)}; useful if that is already available in the caller.}
  \item{\dots}{optional further arguments passed to
    \code{\link{interpSpline}()}.}
}
%% \details{
%% }
\value{
  \describe{
    \item{\code{interpSplineCon()}}{returns the
      \code{\link{interpSpline}()} interpolation spline object.}
    \item{\code{isIsplineCon()}}{is \code{TRUE} (or \code{FALSE}),
      indicating if the convexity/concavity constraints are fulfilled (in
      knot intervals).}
  }
}
\author{Martin Maechler}
\seealso{
 \code{\link{conreg}}, \code{\link{interpSpline}}.
}
\examples{
cc <- conreg(cars[,"speed"], cars[,"dist"], convex=TRUE)
iS <- interpSplineCon(cc)
(isC <- isIsplineCon(cc)) # FALSE: not strictly convex
## Passing the interpolation spline --- if you have it anyway ---
## is more efficient (faster) :
stopifnot(identical(isC,
                    isIsplineCon(cc, isp = iS)))

## the interpolation spline is not quite convex:
plot(cc)
with(cars, points(dist ~ speed, col = adjustcolor(1, 1/2)))
lines(predict(iS, seq(1,28, by=1/4)),
      col = adjustcolor("forest green", 3/4), lwd=2)
}
\keyword{smooth}
\keyword{utilities}


