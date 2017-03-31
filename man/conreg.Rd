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
       method = c("Duembgen06_R", "SR"),
       tol = c(1e-10, 1e-7), maxit = c(500, 20),
       adjTol = TRUE, verbose = FALSE)
}

\arguments{
  \item{x, y}{numeric vectors giving the values of the predictor and
    response variable.  Alternatively a single \dQuote{plotting}
    structure (two-column matrix / y-values only / list, etc) can be
    specified: see \code{\link{xy.coords}}.}
  \item{w}{for \code{method} \code{"Duembgen06_R"} only: optional vector
    of weights of the same length as \code{x}; defaults to all 1.}
  \item{convex}{logical indicating if convex or concave regression is
    desired.}
  \item{method}{a character string indicating the method used,
    \describe{
      \item{"Duembgen06_R"}{is an active set method written by Lutz
	Duembgen (University of Berne, CH) in Matlab in 2006 and
	translated to \R by Martin Maechler.}
      \item{"SR"}{is an \R interface to the C code of a \bold{S}upport
	\bold{R}eduction algorithm written by Piet Groeneboom (TU Delft,
	NL) and donated to the \pkg{cobs} package in July 2012.}
    }}
  \item{tol}{convergence tolerance(s); do not make this too small!}
  \item{maxit}{maximal number of (outer and inner) iterations of knot
    selection.}
  \item{adjTol}{(for \code{"Duembgen06_R"} only:) logical indicating if
    the convergence test tolerance is to be adjusted (increased) in some cases.}
  \item{verbose}{logical or integer indicating if (and how much) knot
    placement and fitting iterations should be \dQuote{reported}.}
}
\details{
  Both algorithms need some numerical tolerances because of rounding
  errors in computation of finite difference ratios.
  The active-set \code{"Duembgen06_R"} method notably has two different
  such tolerances which were both \code{1e-7}\eqn{= 10^{7}} up to March
  2016.

  The two default tolerances (and the exact convergence checks) may
  change in the future, possibly to become more adaptive.
%   ## TODO {this was 'matlab'}:
%   ## --- dessen Zwischenschritte auch graphisch illustriert werden.
}
\value{
  an object of class \code{conreg} which is basically a list with components
  \item{x}{sorted (and possibly aggregated) abscissa values \code{x}.}
  \item{y}{corresponding y values.}
  \item{w}{corresponding weights, \bold{only} for \code{"Duembgen06_R"}.}
  \item{yf}{corresponding fitted values.}
  \item{convex}{logical indicating if a convex or a concave fit has been
    computed.}
  \item{iKnots}{integer vector giving indices of the \emph{knots},
    i.e. locations where the fitted curve has kinks.
    Formally, these are exactly those indices where the constraint is
    fulfilled strictly, i.e., those \eqn{i} where
    \deqn{(m_i - m_{i-1})/(x_i-x_{i-1}) > (m_{i+1} - m_i)/(x_{i+1}-x_i).}{%
      (m[i] - m[i-1])/(x[i] - x[i-1]) > (m[i+1] - m[i])/(x[i+1] - x[i]).}
  }
  \item{call}{the \code{\link{call}} to \code{conreg()} used.}
  \item{iter}{integer (vector of length one or two) with the number of
    iterations used (in the outer and inner loop for \code{"Duembgen06_R"}).}

  %%--- these are not yet okay  for "SR" -- "SR" can also give more ('r','R',..)
  %% \item{deriv.loc}{... FIXME ...}
  %% \item{conv.loc}{... FIXME ...}

  Note that there are several methods defined for \code{conreg} objects,
  see \code{\link{predict.conreg}} or \code{\link{methods}(class = "conreg")}.

  Notably \code{print} and \code{plot}; also
  \code{\link{predict}}, \code{\link{residuals}}, \code{\link{fitted}},
  \code{\link{knots}}.

  Also, \code{\link{interpSplineCon}()} to construct a more smooth
  (\emph{cubic}) spline, and \code{\link{isIsplineCon}()} which checks
  if the int is strictly concave or convex the same as the
  \code{conreg()} result from which it was constructed.
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
set.seed(1)# -> conreg does not give convex cubic

x <- sort(runif(N))
y <- f(x) + 0.2 * rnorm(N)
plot(x,y, cex = 0.6)
lines(xx, f(xx), col = "blue", lty=2)
rc <- conreg(x,y)
lines(rc, col = 2, force.iSpl = TRUE)
 # 'force.iSpl': force the drawing of the cubic spline through the kinks
title("Concave Regression in R")

y2 <- y


## Trivial cases work too:
(r.1 <- conreg(1,7))
(r.2 <- conreg(1:2,7:6))
(r.3  <- conreg(1:3,c(4:5,1)))
}
\keyword{regression}
\keyword{smooth}
