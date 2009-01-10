\name{conreg-methods}
\title{Summary Methods for 'conreg' Objects}
\alias{fitted.conreg}
\alias{residuals.conreg}
\alias{knots.conreg}
\alias{lines.conreg}
\alias{plot.conreg}
\alias{predict.conreg}
\alias{print.conreg}
%
\description{
  Methods for \code{\link{conreg}} objects
}
\usage{
\method{fitted}{conreg}(object, ...)
\method{residuals}{conreg}(object, ...)
\method{knots}{conreg}(Fn, ...)

\method{lines}{conreg}(x, type = "l", col = 2, lwd = 1.5, show.knots = TRUE,
      add.iSpline = TRUE, force.iSpl = FALSE, ...)

\method{plot}{conreg}(x, type = "l", col = 2, lwd = 1.5, show.knots = TRUE,
    add.iSpline = TRUE, force.iSpl = FALSE,
    xlab = "x", ylab = expression(s[c](x)),
    sub = "simple concave regression", col.sub = col, ...)

\method{predict}{conreg}(object, x, deriv = 0, ...)

\method{print}{conreg}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
  \item{object, Fn, x}{an \R object of class \code{conreg}, i.e., typically
    the result of \code{\link{conreg}(..)}.   For \code{predict()},
    \code{x} is a numeric vector of abscissa values at which to evaluate
    the concave/convex spline function.}
  \item{type, col, lwd, xlab, ylab, sub, col.sub}{plotting arguments as in
    \code{\link{plot.default}}.}
  \item{show.knots}{logical indicating the spline knots should be marked
    additionally.}
  \item{add.iSpline}{logical indicating if an \emph{interpolation}
    spline should be considered for plotting. This is only used when it
    is itself concave/convex, \bold{unless} \code{force.iSpl} is
    \code{TRUE}.}
  \item{force.iSpl}{logical indicating if an interpolating spline is drawn
    even when it is not convex/concave.}

  \item{deriv}{for \code{predict}, integer specifying the derivate to be
    computed; currently must be \code{0} or \code{1}.}
  \item{digits}{number of significant digits for printing.}
  \item{\dots}{further arguments, potentially passed to methods.}
}
% \details{
%% FIXME
% }
% \value{
%% FIXME really
% }
\author{Martin Maechler}
\seealso{\code{\link{conreg}}, ....
}
\examples{
example(conreg, echo = FALSE)
class(rc) # "conreg"
rc # calls the print method
knots(rc)
plot(rc)
xx <- seq(-0.1, 1.1, length=201) # slightly extrapolate
yx <- predict(rc, xx)
plot(xx, yx, type = "l",
     main="plot(xx,  predict( conreg(.), xx))")
}
\keyword{print}
\keyword{models}
