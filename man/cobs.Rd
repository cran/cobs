\name{cobs}
\alias{cobs}
\alias{n1000cut}
\title{COnstrained B-Splines Nonparametric Regression Quantiles}
\description{
  Computes constrained quantile curves using linear or
  quadratic splines.  The median spline (\eqn{L_1} loss) is a robust
  (constrained) smoother.
}
\usage{
cobs(x, y, constraint = c("none", "increase", "decrease",
                          "convex", "concave", "periodic"),
     knots, nknots, method = "quantile",
     degree = 2, tau = 0.5, lambda = 0, ic = "aic",
     n.sub = n1000cut(n),
     knots.add = FALSE, pointwise,
     print.warn = TRUE, print.mesg = TRUE, trace = print.mesg,
     coef = rep(0,nvar), w = rep(1,n),
     maxiter = 20*n, lstart = 7872, toler.kn = 1e-6,
     eps = .Machine$double.eps, factor = 1)

n1000cut(n)
}
\arguments{
  \item{x}{vector of covariate; missing values are omitted.}
  \item{y}{vector of response variable.  It must have the same length as
    \code{x}.}
  \item{constraint}{character (string) specifying the kind of
    constraint; must be one of the values in the default list, above;
    may be abbreviated.}
  \item{knots}{vector of locations of the knot mesh; if missing,
    \code{nknots} number of \code{knots} will be created using the
    specified \code{method} and automatic knot selection will be carried
    out for regression B-spline (lambda=0); if not missing and
    \code{length(knots)==nknots}, the provided knot mesh will be used in
    the fit and no automatic knot selection will be performed;
    otherwise, automatic knots selection will be performed on the
    provided \code{knots}.}
  \item{nknots}{maximum number of knots; defaults to 6 for regression
    B-spline, 20 for smoothing B-spline.}
  \item{method}{character specifying the method for generating
    \code{nknots} number of \code{knots} when \code{knots} is not provided;
    "quantile" (equally spaced in percentile levels)
    or "uniform" (equally spaced knots); defaults to "quantile".}
  \item{degree}{degree of the splines; 1 for linear spline and 2 for quadratic
    spline; defaults to 2.}
  \item{tau}{desired quantile level; defaults to 0.5 (median).}
  \item{lambda}{penalty parameter;
    lambda = 0: no penalty (regression B-spline);\cr
    lambda > 0: smoothing B-spline with the given lambda;\cr
    lambda < 0: smoothing B-spline with lambda chosen by an information
    criterion, see \code{ic}.}
  \item{ic}{information criterion used in knot deletion and addition for
    regression B-spline method when lambda=0; "aic" (Akaike-type) or
    "sic" (Schwarz-type); default to "aic".}
  \item{n.sub}{integer, not larger than sample size \code{n}; the
    default has \code{n.sub == n} as long as \code{n} is less than 1000.}
  \item{knots.add}{logical indicating if an additional step of stepwise
    knot addition should be performed for regression B-splines.}
  \item{pointwise}{an optional three-column matrix with each row
    specifies one of the following constraints:
    \describe{
      \item{\code{( 1,xi,yi)}:}{fitted value at xi will be \eqn{\ge}{>=} yi;}
      \item{\code{(-1,xi,yi)}:}{fitted value at xi will be \eqn{\le}{<=} yi;}
      \item{\code{( 0,xi,yi)}:}{fitted value at xi will be \eqn{=} yi;}
      \item{\code{( 2,xi,yi)}:}{derivative of the fitted function at xi
	will be yi.}
      }
  }
  \item{print.warn}{
    logical flag for printing of interactive warning messages; true by default;
    probably needs to be set to \code{FALSE} if performing simulation.}
  \item{print.mesg}{
    logical flag or integer for printing of intermediate messages; true
    by default.  Probably needs to be set to \code{FALSE} in simulations.}
  \item{trace}{integer \eqn{\ge 0}{>= 0} indicating how much the Fortran
    routine \code{drqssbc} should print intermediate messages;
    defaults to \code{print.mesg}.}
  \item{coef}{initial guess of the B-spline coefficients; default to a
    vector of zeros.}
  \item{w}{
    vector of weights the same length as x (y) assigned to both x and y;
    default to uniform weights adding up to one; using normalized weights
    that add up to one will speed up computation.
  }
  \item{maxiter}{upper bound of the number of iteration; default to 20*n.}
  \item{lstart}{starting value for lambda when performing parametric
    programming in lambda if \code{lambda < 0}; defaults to
    \code{log(big)^2}.}
  \item{toler.kn}{numeric tolerance for shifting the boundary knots
    outside; default 1e-6 used to be built in.}
  \item{eps}{tolerance passed to \code{\link{qbsks}} and \code{\link{drqssbc}}.}
  \item{factor}{
    determines how big a step to the next smaller lambda should be while
    performing parametric programming in lambda; the default \code{1}
    will give all unique lambda's; use of a bigger factor \eqn{(> 1 &  < 4)}
    will save time for big problems.}
  %
  \item{n}{integer, the sample size.}
}
\value{
  an object of class \code{cobs}, a list with components
  \item{call}{the \code{cobs(..)} call used for creation.}
  \item{tau, degree}{as input}
  \item{constraint}{as input (but no more abbreviated).}
  \item{call}{the \code{cobs(..)} call used for creation.}
  \item{coef}{B-spline coefficients.}
  \item{knots}{the final set of knots used in the computation.}
  \item{ifl}{exit code:
    \describe{
      \item{1}{-- ok;}
      \item{2}{-- problem is infeasible, check specification of the
	\code{pointwise} argument;}
      \item{3}{-- maxiter is reached before finding a solution,
	either increase maxiter and restart the program with \code{coef}
	set to the value upon previous exit or use a smaller
	\code{lstart} value when lambda<0 or use a smaller \code{lambda}
	value when lambda>0;}
      \item{4}{-- program aborted, numerical difficulties due to
	ill-conditioning.}
    }
  }
  \item{icyc}{length 2: number of cycles taken to achieve convergence
    for final lambda, and total number of cycles for all lambdas.}
  \item{k}{the effective dimensionality of the final fit.}
  \item{k0}{(usually the same)}
  \item{x.ps}{the pseudo design matrix \eqn{X} (as returned by
    \code{\link{qbsks}}).}
  \item{resid}{vector of residuals from the fit.}
  \item{fitted}{vector of fitted values from the fit.}
  \item{SSy}{the sum of squares around centered \code{y} (e.g. for
    computation of \eqn{R^2}.)}
  \item{lambda}{the penalty parameter used in the final fit.}
  \item{pp.lambda}{vector of all unique lambda's obtained from
    parametric programming when \code{lambda} < 0 on input.}
  \item{sic}{vector of Schwarz information criteria evaluated at
    \code{pp.lambda}.}
}%end{value}
\details{
  \code{cobs()} computes the constraint quantile smoothing B-spline with
  penalty when lambda is not zero.\cr
  If lambda < 0, an optimal lambda will be chosen using Schwarz type
  information criterion. \cr
  If lambda > 0, the supplied lambda will be used.\cr
  If lambda = 0, cobs computes the constraint quantile regression B-spline
  with no penalty using the provided knots or those selected by Akaike or
  Schwarz information criterion.
}
\references{
  He, X. and Ng, P. (1999)
  COBS: Qualitatively Constrained Smoothing via Linear Programming;
  \emph{Computational Statistics} \bold{14}, 315--337.

  Koenker, R. and Ng, P. (1996)
  A Remark on Bartels and Conn's Linearly Constrained L1 Algorithm,
  \emph{ACM Transaction on Mathematical Software} \bold{22}, 493--495.

  Ng, P. (1996)
  An Algorithm for Quantile Smoothing Splines,
  \emph{Computational Statistics & Data Analysis} \bold{22}, 99--118.

  Bartels, R. and Conn A. (1980)
  Linearly Constrained Discrete \eqn{L_1} Problems,
  \emph{ACM Transaction on Mathematical Software} \bold{6}, 594--608.

  A postscript version of the paper that describes the details of COBS
  can be downloaded from \url{http://www.cba.nau.edu/pin-ng/cobs.html}
}
\section{Warning}{
  This is still a beta version, and we do appreciate comments and
  suggestions; \code{library(help = cobs)} shows the authors.
}
\seealso{\code{\link[modreg]{smooth.spline}} for unconstrained smoothing
    splines; \code{\link[splines]{bs}} for unconstrained (regression)
    B-splines.}
}
\examples{
x <- seq(-1,1,,50)
y <- (f.true <- pnorm(2*x)) + rnorm(50)/10
## specify pointwise constraints (boundary conditions)
con <- rbind(c( 1,min(x),0), # f(min(x)) >= 0
             c(-1,max(x),1), # f(max(x)) <= 1
             c(0,  0,   0.5))# f(0)      = 0.5

## obtain the median regression B-spline using automatically selected knots
Rbs <- cobs(x,y,constraint="increase",pointwise=con)
Rbs

plot(x,y)
lines(predict(Rbs), col = 2, lwd = 1.5)
lines(spline(x,f.true), col = "gray40")

Rbsub <- cobs(x,y,constraint="increase",pointwise=con, n.sub = 45)
%%--> Current newly introduced BUG -- ../R/cobs.R.~21~ *does* work
summary(Rbsub)
lines(predict(Rbsub), col = 4, lwd = 1)


## compute the median smoothing B-spline using automatically chosen lambda
Sbs <- cobs(x,y,constraint="increase",pointwise=con,lambda=-1)
Sbs
plot(Sbs$pp.lambda[-1], Sbs$sic[-1], log = "x",
     main = "SIC ~ lambda", xlab = expression(lambda), ylab = "SIC")
axis(1, at = Sbs$lambda, label = expression(hat(lambda)),
     col.axis = 2, mgp = c(3, 0.5, 0))

Sb1 <- cobs(x,y,constraint="increase",pointwise=con,lambda=-1, degree=1)
summary(Sb1)
pxx <- predict(Sb1, xx <- seq(-1.2, 1.2, len = 201), interval = "both")
plot(x,y, main = deparse(Sb1$call),
     xlim = range(xx), ylim = range(y, pxx[,"fit"]))
lines(pxx, col = 2)
rug(Sb1$knots, col = 4, lwd = 1.6)# (too many knots)
matlines(pxx[,1], pxx[,-(1:2)], col= rep(c("green","blue"),c(2,2)), lty=2)
}
% Converted by Sd2Rd version 1.21.
\keyword{smooth}
\keyword{regression}
