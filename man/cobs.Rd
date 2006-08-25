\name{cobs}
\alias{cobs}
\alias{cut00}
\title{COnstrained B-Splines Nonparametric Regression Quantiles}
\description{
  Computes constrained quantile curves using linear or
  quadratic splines.  The median spline (\eqn{L_1} loss) is a robust
  (constrained) smoother.
}
\usage{
cobs(x, y, constraint = c("none", "increase", "decrease",
                          "convex", "concave", "periodic"),
     w = rep(1,n),
     knots, nknots = if(lambda == 0) 6 else 20,
     method = "quantile", degree = 2, tau = 0.5,
     lambda = 0, ic = c("sic", "aic", "bic", "SIC", "AIC", "BIC"),
     knots.add = FALSE, repeat.delete.add = FALSE, pointwise = NULL,
     keep.data = TRUE, keep.x.ps = TRUE,
     print.warn = TRUE, print.mesg = TRUE, trace = print.mesg,
     lambdaSet = exp(seq(log(lambda.lo), log(lambda.hi), length = lambda.length)),
     lambda.lo = f.lambda*1e-4, lambda.hi = f.lambda*1e3, lambda.length = 25,
     maxiter = 100,
     rq.tol = 1e-8, toler.kn = 1e-6, tol.0res = 1e-6,
     eps, n.sub, coef, lstart, factor)
}
\arguments{
  \item{x}{vector of covariate; missing values are omitted.}
  \item{y}{vector of response variable.  It must have the same length as
    \code{x}.}
  \item{constraint}{character (string) specifying the kind of
    constraint; must be one of the values in the default list above;
    may be abbreviated.  More flexible constraints can be specified via
    the \code{pointwise} argument (below).}
  \item{w}{
    vector of weights the same length as \code{x} (y) assigned to both x and y;
    default to uniform weights adding up to one; using normalized weights
    that add up to one will speed up computation.
  }
  \item{knots}{vector of locations of the knot mesh; if missing,
    \code{nknots} number of \code{knots} will be created using the
    specified \code{method} and automatic knot selection will be carried
    out for regression B-spline (\code{lambda=0}); if not missing and
    \code{length(knots)==nknots}, the provided knot mesh will be used in
    the fit and no automatic knot selection will be performed;
    otherwise, automatic knots selection will be performed on the
    provided \code{knots}.}
  \item{nknots}{maximum number of knots; defaults to 6 for regression
    B-splines, 20 for smoothing B-splines.}
  \item{method}{character string specifying the method for generating
    \code{nknots} number of \code{knots} when \code{knots} is not provided;
    either \code{"quantile"} (equally spaced in percentile levels)
    or \code{"uniform"} (equally spaced knots); defaults to "quantile".}
  \item{degree}{degree of the splines; 1 for linear spline and 2 for quadratic
    spline; defaults to 2.}
  \item{tau}{desired quantile level; defaults to 0.5 (median).}
  \item{lambda}{penalty parameter \eqn{\lambda}\cr
    \eqn{\lambda = 0}: no penalty (regression B-spline);\cr
    \eqn{\lambda > 0}: smoothing B-spline with the given \eqn{\lambda};\cr
    \eqn{\lambda < 0}: smoothing B-spline with \eqn{\lambda} chosen by a
    Schwarz-type information criterion.}
  \item{ic}{string indicating the information  criterion used in knot
    deletion and addition for regression B-spline method when
    \code{lambda == 0};\cr
    \code{"aic"} (Akaike-type, equivalently \code{"AIC"}) or\cr
    \code{"sic"} (Schwarz-type, equivalently \code{"BIC"}, \code{"SIC"}
    or \code{"bic"}).  Defaults to \code{"sic"}.}
  \item{knots.add}{logical indicating if an additional step of stepwise
    knot addition should be performed for regression B-splines.}
  \item{repeat.delete.add}{logical indicating if an additional step of stepwise
    knot deletion should be performed for regression B-splines.}
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
  \item{keep.data}{logical indicating if the \code{x} and \code{y} input
    vectors \bold{after} removing \code{\link{NA}}s should be kept in
    the result.}
  \item{keep.x.ps}{logical indicating if the pseudo design matrix
    \eqn{\tilde{X}}{X\~} should be returned (as \emph{sparse} matrix).
    That is needed for interval prediction, \code{\link{predict.cobs}(*,
      interval=..)}.}
  \item{print.warn}{flag for printing of interactive warning messages;
    true by default; set to \code{FALSE} if performing simulation.}
  \item{print.mesg}{
    logical flag or integer for printing of intermediate messages; true
    by default.  Probably needs to be set to \code{FALSE} in simulations.}
  \item{trace}{integer \eqn{\ge 0}{>= 0} indicating how much diagnostics
    the low-level code in \code{\link{drqssbc2}} should print;  defaults
    to \code{print.mesg}.}
  \item{lambdaSet}{numeric vector of lambda values to use for grid search;
    in that case, defaults to a geometric sequence (a \dQuote{grid in
      log scale}) from \code{lambda.lo} to \code{lambda.hi} of length
    \code{lambda.length}.}
  \item{lambda.lo, lambda.hi}{lower and upper bound of the grid search
    for lambda (when \code{lambda < 0}).  The defaults are meant to keep
    everything scale-equivariant and are hence using the factor
%%    \eqn{f = \sigma_y\cdot\sigma_x^d}{f = s[y] * s[x]^d}, i.e.,
%%    \code{f.lambda <- sd(y) * sd(x)^degree}.
    \eqn{f = \sigma_x^d}{f = s[x]^d}, i.e.,
    \code{f.lambda <- sd(x)^degree}.
    %%% FIXME ?!? in quantreg
    Note however that currently the underlying algorithms in package
    \pkg{quantreg} are \emph{not} scale equivariant yet.}
  \item{lambda.length}{number of grid points in the grid search for
    optimal lambda.}
  \item{maxiter}{upper bound of the number of iterations; defaults to 100.}
  \item{rq.tol}{numeric convergence tolerance for the interior point
    algorithm called from \code{\link[quantreg]{rq.fit.sfnc}()} or
    \code{\link[quantreg]{rq.fit.sfn}()}.}
  \item{toler.kn}{numeric tolerance for shifting the boundary knots
    outside; defaults to \eqn{10^{-6}}{10^(-6)}.}
  \item{tol.0res}{tolerance for testing \eqn{|r_i| = 0},  passed to
    \code{\link{qbsks2}} and \code{\link{drqssbc2}}.}
  \item{eps, n.sub, coef, lstart, factor}{unused argument of the 1999/2002 version
    of \code{cobs()}, currently still allowed (with a
    \code{\link{warning}}) for back compatibility reasons.}
}
\value{
  an object of class \code{cobs}, a list with components
%%% FIXME: needs to be updated!
  \item{call}{the \code{cobs(..)} call used for creation.}
  \item{tau, degree}{as input}
  \item{constraint}{as input (but no more abbreviated).}
  \item{pointwise}{as input.}
  \item{coef}{B-spline coefficients.}
  \item{knots}{the final set of knots used in the computation.}
  \item{ifl}{exit codes: see \code{\link[quantreg]{rq.fit.sfnc}}.}
  \item{icyc}{length 2: number of cycles taken to achieve convergence
    for final lambda, and total number of cycles for all lambdas.}
  \item{k}{the effective dimensionality of the final fit.}
  \item{k0}{(usually the same)}
  \item{x.ps}{the pseudo design matrix \eqn{X} (as returned by
    \code{\link{qbsks2}}).}
  \item{resid}{vector of residuals from the fit.}
  \item{fitted}{vector of fitted values from the fit.}
  \item{SSy}{the sum of squares around centered \code{y} (e.g. for
    computation of \eqn{R^2}.)}
  \item{lambda}{the penalty parameter used in the final fit.}
  \item{pp.lambda}{vector of all lambdas used for
    lambda search when \code{lambda} < 0 on input.}
  \item{pp.sic}{vector of Schwarz information criteria evaluated at
    \code{pp.lambda}; note that it is not quite sure how good these are
    for determining an optimal \code{lambda}.}
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
  Koenker, R. and Ng, P. (2005)
  Inequality Constrained Quantile Regression,
  \emph{Sankhya, The Indian Journal of Statistics} \bold{67}, 418--440.

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
% \section{Warning}{
%   This is still a beta version, and we do appreciate comments and
%   suggestions; \code{library(help = cobs)} shows the authors.
% }
\seealso{\code{\link{smooth.spline}} for unconstrained smoothing
    splines; \code{\link[splines]{bs}} for unconstrained (regression)
    B-splines.}
}
\examples{
x <- seq(-1,3,,150)
y <- (f.true <- pnorm(2*x)) + rnorm(150)/10
## specify pointwise constraints (boundary conditions)
con <- rbind(c( 1,min(x),0), # f(min(x)) >= 0
             c(-1,max(x),1), # f(max(x)) <= 1
             c(0,  0,   0.5))# f(0)      = 0.5

## obtain the median  REGRESSION  B-spline using automatically selected knots
Rbs <- cobs(x,y, constraint= "increase", pointwise = con)
Rbs
plot(Rbs, lwd = 2.5)
lines(spline(x, f.true), col = "gray40")
lines(predict(cobs(x,y)), col = "blue")
mtext("cobs(x,y)   # completely unconstrained", 3, col= "blue")

## compute the median  SMOOTHING  B-spline using automatically chosen lambda
Sbs <- cobs(x,y, constraint="increase", pointwise= con, lambda= -1)
summary(Sbs)
plot(Sbs) ## by default  includes  SIC ~ lambda

Sb1 <- cobs(x,y, constraint="increase", pointwise= con, lambda= -1,
            degree = 1)
summary(Sb1)
plot(Sb1)

plot(Sb1, which = 2) # only the  data + smooth
rug(Sb1$knots, col = 4, lwd = 1.6)# (too many knots)
mtext(" + pointwise and simultaneous 95\% - confidence intervals")
xx <- seq(min(x) - .2, max(x)+ .2, len = 201)
pxx <- predict(Sb1, xx, interval = "both")
lines(pxx, col = 2)
matlines(pxx[,1], pxx[,-(1:2)], col= rep(c("green3","blue"), c(2,2)), lty=2)
}
\keyword{smooth}
\keyword{regression}
