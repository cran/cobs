\name{drqssbc}
\alias{drqssbc}
\title{Regression Quantile Smoothing Spline with Constraints}
\description{
  Estimate the B-spline coefficients for a regression quantile
  \emph{smoothing} spline with optional constraints, using Ng(1996)'s
  algorithm.
}
\usage{
drqssbc(x, y, w, pw, knots, degree, Tlambda, constraint,
        n.sub = n1000cut(n),
        equal, smaller, greater, gradient, coef, maxiter,
        trace, n.equal, n.smaller, n.greater, n.gradient,
        nrq, nl1, neqc, niqc, nvar, nj0,
        tau, lam, tmin, kmax, lstart, factor,
        eps = .Machine$double.eps, print.warn)
}
\arguments{
  \item{x}{numeric vector, sorted increasingly, the abscissa values}
  \item{y}{numeric, same length as \code{x}, the observations.}
  \item{w}{numeric vector of weights, same length as \code{x}, as in
    \code{\link{cobs}}.}
  \item{pw}{penalty weights vector passed to \code{\link{l1.design}} or
    \code{\link{loo.design}}.}
  \item{knots}{ ~~Describe \code{knots} here~~ }
  \item{degree}{integer, must be 1 or 2.}
  \item{Tlambda}{ ~~Describe \code{Tlambda} here~~ }
  \item{constraint}{see \code{\link{cobs}} (but cannot be abbreviated here).}
  \item{n.sub}{integer, not larger than sample size \code{n}; the
    default has \code{n.sub == n} as long as \code{n} is less than 1000.}
  \item{equal,smaller, greater}{3-column matrices specifying the
    respective constraints.  The has 0 zeros if there no constraints of
    the corresponding kind.}
  \item{gradient}{3-column matrix for gradient constraints.}
  \item{coef}{numeric vector, the initial guess for the B-spline coefficients.}
  \item{maxiter}{upper bound of the number of iteration; default to 20*n.}
  \item{trace}{ ~~Describe \code{trace} here~~ }
  \item{n.equal,n.smaller,n.greater,n.gradient}{
    ~~Describe \code{n.gradient} here~~ }
  \item{nrq}{ ~~Describe \code{nrq} here~~ }
  \item{nl1}{ ~~Describe \code{nl1} here~~ }
  \item{neqc}{integer giving the number of equations.}
  \item{niqc}{integer giving the number of \bold{i}ne\bold{q}uality
    \bold{c}onstraints.}
  \item{nvar}{integer giving the number of equations \emph{and} constraints.}
  \item{nj0}{ ~~Describe \code{nj0} here~~ }
  \item{tau}{desired quantile level; defaults to 0.5 (median).}
  \item{lam}{ ~~Describe \code{lam} here~~ }
  \item{tmin}{ ~~Describe \code{tmin} here~~ }
  \item{kmax}{ ~~Describe \code{kmax} here~~ }
  \item{lstart}{number, see \code{\link{cobs}.}}
  \item{factor}{number in [1,4], see \code{\link{cobs}}.}
  \item{eps}{tolerance used in the fortran code in many different contexts.}
  \item{print.warn}{logical indicating if warnings should be printed,
    when the algorithm seems to have behaved somewhat unexpectedly.}
}
\details{
  This is an auxiliary function for \code{\link{cobs}}, possibly
  interesting on its own.  This documentation is currently sparse; read
  the source code!
}
\value{
  a list with components
  \item{comp1 }{Description of `comp1'}
  \item{comp2 }{Description of `comp2'}
  ...
}
\references{
  Ng, P. (1996)
  An Algorithm for Quantile Smoothing Splines,
  \emph{Computational Statistics \& Data Analysis} \bold{22}, 99--118.
}
\author{Pin Ng; this help page: Martin Maechler.}
\seealso{
  the main function \code{\link{cobs}} and its auxiliary
  \code{\link{qbsks}} which calls \code{drqssbc()} repeatedly.
}
\examples{
set.seed(1243)
x  <- 1:32
fx <- (x-5)*(x-15)^2*(x-21)
y  <- fx + round(rnorm(x,s = 0.25),2)
## FAILS  drqssbc(x,y,nrq=32,lam=1,degree=1,knots=c(1,5,15,32))
}
\keyword{smooth}
\keyword{utilities}
