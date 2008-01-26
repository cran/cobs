\name{predict.cobs}
\alias{predict.cobs}
\title{Predict method for COBS Fits}
\description{
  Compute predicted values and simultaneous or pointwise confidence
  bounds for \code{\link{cobs}} objects.
}
\usage{
\method{predict}{cobs}(object,
        z, minz = knots[1], maxz = knots[nknots], nz = 100,
        interval = c("none", "confidence", "simultaneous", "both"),
        level = 0.95, \dots)
}
\arguments{
  \item{object}{object of class \code{cobs}.}
  \item{z}{vector of grid points at which the fitted values are
    evaluated; defaults to an equally spaced grid with \code{nz} grid
    points between \code{minz} and \code{maxz}.  Note that now \code{z}
    may lie outside of the knots interval which was not allowed originally.}
  \item{minz}{numeric needed if \code{z} is not specified; defaults to
    \code{min(x)} or the first knot if \code{knots} are given.}
  \item{maxz}{analogous to \code{minz}; defaults to \code{max(x)} or the
    last knot if \code{knots} are given.}
  \item{nz}{number of grid points in \code{z} if that is not given;
    defaults to 100.}
  \item{interval}{type of interval calculation, see below}
  \item{level}{confidence level}
  \item{\dots}{further arguments passed to and from methods.}
}
% \details{

% }
\value{
  a matrix of predictions and bounds if \code{interval} is set (not
  "none").  The columns are named \code{z}, \code{fit}, further
  \code{cb.lo} and \code{cb.up} for the \code{simultaneous} confidence
  band, and \code{ci.lo} and \code{ci.up} the pointwise
  \code{confidence} intervals according to specified \code{level}.

  If \code{z} has been specified, it is unchanged in the result.
}
\author{Martin Maechler, based on He and Ng's code in \code{cobs()}.}
\seealso{\code{\link{cobs}} the model fitting function.}
\examples{
example(cobs) # continuing :
(pRbs <- predict(Rbs))
#str(pSbs <- predict(Sbs, xx, interval = "both"))
str(pSbs <- predict(Sbs, xx, interval = "none"))

plot(x,y, xlim = range(xx), ylim = range(y, pSbs[,2], finite = TRUE),
     main = "COBS Median smoothing spline, automatical lambda")
lines(pSbs, col = "red")
lines(spline(x,f.true), col = "gray40")
#matlines(pSbs[,1], pSbs[,-(1:2)],
#         col= rep(c("green","blue"),c(2,2)), lty=2)
}
\keyword{regression}
