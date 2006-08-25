\name{cobs-methods}
\alias{coef.cobs}
\alias{fitted.cobs}
\alias{knots.cobs}
\alias{print.cobs}
\alias{residuals.cobs}
\alias{summary.cobs}
\title{Methods for COBS Objects}
\description{
  Print, summary and other methods for \code{\link{cobs}} objects.
}
\usage{
\method{print}{cobs}(x, digits = getOption("digits"), \dots)
\method{summary}{cobs}(object, digits = getOption("digits"), \dots)

\method{coef}{cobs}(object, \dots)
\method{fitted}{cobs}(object, \dots)
\method{knots}{cobs}(Fn, \dots)
\method{residuals}{cobs}(object, \dots)
}
\arguments{
  \item{x, object, Fn}{object of class \code{cobs}.}
  \item{digits}{number of digits to use for printing.}
  \item{\dots}{further arguments passed from and to methods.}
}
\value{
  \code{print.cobs()} returns its argument invisibly.
  The \code{coef()}, \code{fitted()}, \code{knots()}, and
  \code{residuals()} methods return a numeric vector.
}
\details{
  These are methods for fitted COBS objects, as computed by
  \code{\link{cobs}}.
}
\author{Martin Maechler}
\seealso{
  \code{\link{predict.cobs}} for the \code{\link{predict}} method,
  \code{\link{plot.cobs}} for the \code{\link{plot}} method,
  and \code{\link{cobs}} for examples.
}
\examples{
example(cobs)% ./cobs.Rd
Sbs # uses print.*

summary(Sbs)

coef(Sbs)
knots(Sbs)
}
\keyword{print}
