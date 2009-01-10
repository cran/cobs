\name{conreg}
\alias{conreg}
\title{Convex / Concave Regression}
\description{
  Compute a univariate concave or convex regression, i.e.,
  for given vectors, \eqn{x,y,w} in \eqn{R^n}, where \eqn{x} has to be
  strictly sorted (\eqn{x_1 < x_2 < \ldots < x_n}), compute an
  \eqn{n}-vector \eqn{m} minimizing the weighted sum of squares
  \eqn{\sum_{i=1}^n {w_i (y_i - m_i)^2}}{sum(i=1..n; w_i * (y_i - m_i)^2)}
  under the constraints
  \deqn{(m_i - m_{i-1})/(x_i - x_{i-1}) \ge (m_{i+1} - m_i)/(x_{i+1} - x_i),}{%
    (m[i] - m[i-1])/(x[i] - x[i-1]) >= (m[i+1] - m[i])/(x[i+1] - x[i]),}
  for \eqn{1 \le i \le n}{1 <= i <= n} and
  \eqn{m_0 := m_{n+1} := -\infty}{m[0] := m[n+1] := -Inf},
  for concavity.
  For convexity (\code{convex=TRUE}), replace \eqn{\ge}{>=} by
  \eqn{\le}{<=} and \eqn{-\infty}{-Inf} by \eqn{+\infty}{+Inf}.
}
\usage{
conreg(x, y = NULL, w = NULL, convex = FALSE,
       tol = 1e-07, maxit = c(200, 20), adjTol = TRUE, verbose = FALSE)
}

\arguments{
  \item{x, y}{numeric vectors giving the values of the predictor and
    response variable.  Alternatively a single \dQuote{plotting}
    structure (two-column matrix / y-values only / list, etc) can be
    specified: see \code{\link{xy.coords}}.}
  \item{w}{optional vector of weights of the same length as \code{x};
    defaults to all 1.}
  \item{convex}{logical indicating if convex or concave regression is desired.}
  \item{tol}{convergence tolerance; do not make this too small!}
  \item{maxit}{maximal number of (outer and inner) iterations of knot
    selection.}
  \item{adjTol}{logical indicating if the convergence test tolerance is
    to be adjusted (increased) in some cases.}
  \item{verbose}{logical indicating if knot placement iterations should
    be \dQuote{reported}.}
}
\details{
  The algorithm is an active-set method, needing some numerical
  tolerance because of rounding errors in computation of finite
  difference rations.

%   ## TODO {this was 'matlab'}:
%   ## --- dessen Zwischenschritte auch graphisch illustriert werden.
}
\value{
  an object of class \code{conreg} which is basically a list with components
  \item{x}{sorted (and possibly aggregated) abscissa values \code{x}.}
  \item{y}{corresponding y values.}
  \item{yf}{corresponding fitted values.}
  \item{iKnots}{integer vector giving indices of the \emph{knots},
    i.e. locations where the fitted curve has kinks.
    Formally, these are exactly those indices where the constraint is
    fulfilled strictly, i.e., those \eqn{i} where
    \deqn{(m_i - m_{i-1})/(x_i-x_{i-1}) > (m_{i+1} - m_i)/(x_{i+1}-x_i).}{%
      (m[i] - m[i-1])/(x[i] - x[i-1]) > (m[i+1] - m[i])/(x[i+1] - x[i]).}
  }
  \item{call}{the \code{\link{call}} to \code{conreg()} used.}

  \item{..}{... FIXME ...}
  \item{..}{... FIXME ...}

  Note that there are several methods defined for \code{conreg} objects,
  see \code{\link{predict.conreg}}.
  Notably \code{print} and \code{plot}; also
  \code{\link{predict}}, \code{\link{residuals}}, \code{\link{fitted}},
  \code{\link{knots}}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{Lutz Duembgen programmed the original Matlab code in July 2006;
  Martin Maechler ported it to \R, tested, catch infinite loops, added
  more options, improved tolerance, etc;  from April 27, 2007.
%  Ported to R and enhanced:
%  - work for unordered, even duplicated 'x'
%  - made slightly faster; more options;
%  - detect infinite loop; auto-adjust tol;
%  - changed tol (prec) to 1e-7
%  - new arg.  convex == FALSE   <==>  concave == TRUE
%  - define "class" and many methods, plot, predict, ...
%
%  Martin Maechler, 27.-28. Apr 2007
}

\seealso{\code{\link[stats]{isoreg}} for isotone (monotone) regression;
  CRAN packages \pkg{ftnonpar}, \pkg{cobs}, \pkg{logcondens}.
}
\examples{

## Generated data :
N <- 100
f <- function(X) 4*X*(1 - X)

xx <- seq(0,1, length=501)# for plotting true f()
set.seed(1)
x <- sort(runif(N))
y <- f(x) + 0.2 * rnorm(N)
plot(x,y, cex = 0.6)
lines(xx, f(xx), col = "blue", lty=2)

rc <- conreg(x,y)
lines(rc, col = 2)
title("Concave Regression in R")

## Trivial cases work too:
r.1 <- conreg(1,7)
r.2 <- conreg(1:2,7:6)
r.3  <- conreg(1:3,c(4:5,1))
r.3. <- conreg(1:3,c(4:5,8))
stopifnot(resid(r.1) == 0,
          resid(r.2) == 0,
          resid(r.3) == 0,
          all.equal(fitted(r.3.),
                    c(11,17,23)/3, tol=1e-12))
}
\keyword{regression}
\keyword{smooth}
