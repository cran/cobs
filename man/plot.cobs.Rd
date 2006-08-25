\name{plot.cobs}
\alias{plot.cobs}
\title{Plot Method for COBS Objects}
\description{
  The \code{\link{plot}} method for \code{\link{cobs}} objects.
  If there was \code{lambda} selection, it provides two plots by default.
}
\usage{
\method{plot}{cobs}(x, which = if(x$select.lambda) 1:2 else 2,
     show.k = TRUE,
     col = par("col"), l.col = c("red","pink"), k.col = gray(c(0.6, 0.8)),
     lwd = 2, cex = 0.4,
     xlab = deparse(x$call[[2]]),
     ylab = deparse(x$call[[3]]),
     main = paste(deparse(x$call, width.cut = 100), collapse="\n"),
     subtit= c("choosing lambda", "data & spline curve") , \dots)
}
\arguments{
  \item{x}{object of class \code{cobs}.}
  \item{which}{integer vector specifying \emph{which} plots should be
    drawn;} %% 1: SIC ~ lambda  2: data + smooth curve
  \item{show.k}{logical indicating if the \dQuote{effective
      dimensionality} \eqn{k} should also be shown.  Only applicable when
    \code{which} contains 1.}
  \item{col, l.col, k.col}{colors for plotting; \code{k.col} only
    applies when \code{show.k} is true in the first plot (\code{which
      == 1}) where \code{l.col[2]} and \code{k.col[2]} are only used as well,
    for denoting \dQuote{doubtful} points; \code{col} is only used for
    the 2nd plot (\code{which == 2}).}
  \item{lwd, cex}{line width and point size for the 2nd plot
    (i.e. \code{which == 2}).}
  \item{xlab, ylab, main}{axis annotation; see also \code{\link{axis}}.}
  \item{subtit}{a vector of length 2, specifying a sub title for each
    plot (according to \code{which}).}
  \item{\dots}{further arguments passed and to internal
    \code{\link{plot}} methods.}
}
\details{
  \code{plot(.)} produces two side-by-side plots in case there was a
  search for the optimal lambda(\code{which = 1:2}), and only the
  (second) data plus spline curve plot otherwise (\code{which = 2}).

  %% FIXME
}
\author{Martin Maechler}
\seealso{
  There are several other methods for COBS objects, see, e.g.
  \code{\link{summary.cobs}} or \code{\link{predict.cobs}}.

  \code{\link{cobs}} for examples.
}
\examples{
example(cobs)% ./cobs.Rd

plot(Sbs)
plot(fitted(Sbs), resid(Sbs),
     main = "Tukey-Anscombe plot for cobs()",
     sub = deparse(Sbs$call))
}
\keyword{print}
