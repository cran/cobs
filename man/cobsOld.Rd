\name{cobsOld}
\alias{cobsOld}
\title{COnstrained B-Splines Nonparametric Regression Quantiles}
}
\description{
  Computes constrained quantile curves using linear or
  quadratic splines.  The median spline (\eqn{L_1} loss) is a robust
  (constrained) smoother.
}
\usage{
cobsOld(x, y, constraint = c("none", "increase", "decrease",
                          "convex", "concave", "periodic"),
     z, minz = knots[1], maxz = knots[nknots], nz = 100,
     knots, nknots, method = "quantile",
     degree = 2, tau = 0.5, lambda = 0, ic = "aic",
     knots.add = FALSE, alpha = 0.1, pointwise,
     print.warn = TRUE, print.mesg = TRUE, trace = print.mesg,
     coef = rep(0,nvar), w = rep(1,n),
     maxiter = 20*n, lstart = log(big)^2, toler = 1e-6,
     factor = 1)
}
\arguments{
  \item{x}{vector of covariate.}
  \item{y}{vector of response variable. It must have the same length as x.}
  \item{constraint}{character (string) specifying the kind of constraint;
    must be one of "increase", "decrease", "convex", "concave",
    "periodic" or "none".}
  \item{z}{vector of grid points at which the fitted values are
    evaluated; default to an equally spaced grid with \code{nz} grid
    points between \code{minz} and \code{maxz}.  If the fitted values at
    x are desired, use \code{z = unique(x)}.}
  \item{minz}{numeric needed if \code{z} is not specified; defaults to
    \code{min(x)} or the first knot if \code{knots} are given.}
  \item{maxz}{analogous to \code{minz}; defaults to \code{max(x)} or the
    last knot if \code{knots} are given.}
  \item{nz}{number of grid points in \code{z} if that is not given;
    defaults to 100.}
  \item{knots}{vector of locations of the knot mesh; if missing,
    \code{nknots} number of \code{knots} will be created using the
    specified \code{method} and automatic knot selection will be carried
    out for regression B-spline (lambda=0); if not missing and
    \code{length(knots)==nknots}, the provided knot mesh will be used in
    the fit and no automatic knot selection will be performed;
    otherwise, automatic knots selection will be performed on the
    provided \code{knots}.}
  \item{nknots}{maximum number of knots; defaults to 6 for regression
    B-spline, 10 for smoothing B-spline.}
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
    lambda < 0: smoothing B-spline with lambda chosen by a Schwarz-type
    information criterion.}
  \item{ic}{information criterion used in knot deletion and addition for
    regression B-spline method when lambda=0; "aic" (Akaike-type) or
    "sic" (Schwarz-type); default to "aic".}
  \item{knots.add}{logical indicating if an additional step of stepwise
    knot addition should be performed for regression B-splines.}
  \item{alpha}{level of significance for the confidence band.}
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
    logical flag for printing of interactive warning messages; default to T;
    probably needs to be set to F if performing monte carlo simulation.
  }
  \item{print.mesg}{
    logical flag for printing of intermediate messages; default to T;
    probably needs to be set to F if performing monte carlo simulation.
  }
  \item{trace}{integer \eqn{\ge 0}{>= 0} indicating how much the Fortran
    routine \code{drqssbc} should print intermediate messages;
    defaults to \code{print.mesg}, i.e. 1 (or 0).}
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
  \item{toler}{numeric tolerance for ???? ; default 1e-6 used to be builtin.}
  \item{factor}{
    determines how big a step to the next smaller lambda should be while
    performing parametric programming in lambda; the default \code{1}
    will give all unique lambda's; use of a bigger factor \eqn{(> 1 &  < 4)}
    will save time for big problems.}
}
\value{
  a list with components
  \item{coef}{B-spline coefficients.}
  \item{fit}{fitted value at z.}
  \item{resid}{vector of residuals from the fit.}
  \item{z}{as in input.}
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
  \item{icyc}{number of cycles taken to achieve convergence.}
  \item{k}{the effective dimensionality of the final fit.}
  \item{lambda}{the penalty parameter used in the final fit.}
  \item{pp.lambda}{vector of all unique lambda's obtained from
    parametric programming when lambda < 0 on input.}
  \item{sic}{vector of Schwarz information criteria evaluated at pp.lambda.}
  \item{cb.lo}{lower bound of the confidence band}
  \item{cb.up}{upper bound of the confidence band}
  \item{ci.lo}{lower bound of the pointwise confidence interval}
  \item{ci.up}{upper bound of the pointwise confidence interval}
}%end{value}
\details{
  \code{cobsOld()} computes the constraint quantile smoothing B-spline with
  penalty when lambda is not zero.\cr
  If lambda < 0, an optimal lambda will be chosen using Schwarz type
  information criterion. \cr
  If lambda > 0, the supplied lambda will be used.\cr
  If lambda = 0, cobsOld computes the constraint quantile regression B-spline
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
cobsOld(x,y,constraint="increase",pointwise=con)->cobs.o

plot(x,y)
lines(cobs.o$z, cobs.o$fit, col = 2, lwd = 1.5)
lines(spline(x,f.true), col = "gray40")
lines(cobs.o$z, cobs.o$cb.lo,lty=2, col = 3)
lines(cobs.o$z, cobs.o$cb.up,lty=2, col = 3)

## compute the median smoothing B-spline using automatically chosen lambda
cobsOld(x,y,constraint="increase",pointwise=con,lambda=-1)->cobs.oo
plot(x,y, main = "COBS Median smoothing spline, automatical lambda")
lines(spline(x,f.true), col = "gray40")
lines(cobs.oo$z, cobs.oo$fit)
lines(cobs.oo$z, cobs.oo$cb.lo,lty=2)
lines(cobs.oo$z, cobs.oo$cb.up,lty=2)

plot(cobs.oo$pp.lambda[-1], cobs.oo$sic[-1], log = "x",
     main = "SIC ~ lambda", xlab = expression(lambda), ylab = "SIC")
axis(1, at = cobs.oo$lambda, label = expression(hat(lambda)),
     col.axis = 2, mgp = c(3, 0.5, 0))
}
% Converted by Sd2Rd version 1.21.
\keyword{smooth}
\keyword{regression}
