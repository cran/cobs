\name{cobs-methods}
\alias{print.cobs}
\alias{summary.cobs}
\alias{fitted.cobs}
\alias{residuals.cobs}
\title{Methods for COBS Objects}
\description{
  Print, summary and other methods for \code{\link{cobs}} objects.
}
\usage{
\method{print}{cobs}(x, digits = getOption("digits"), \dots)
\method{summary}{cobs}(object, digits = getOption("digits"), \dots)

\method{fitted}{cobs}(object, \dots)
\method{residuals}{cobs}(object, \dots)

}
\arguments{
  \item{x, object}{object of class \code{cobs}.}
  \item{digits}{number of digits to use for printing.}
  \item{\dots}{further arguments passed from and to methods.}
}
\value{
  \code{print.cobs()} returns its argument invisibly.
}
\details{
  These are methods for fitted COBS objects, as computed by
  \code{\link{cobs}}.
}
\author{Martin Maechler}
\seealso{
  \code{\link{predict.cobs}} for the \code{\link{predict}} method, and
  \code{\link{cobs}} for examples.
}
\examples{
example(cobs)% ./cobs.Rd
summary(Sbs)

plot(fitted(Sbs), resid(Sbs), main = deparse(Sbs$call))
}
\keyword{print}
