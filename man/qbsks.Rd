\name{qbsks}
\alias{qbsks}
\title{Quantile B-Spline with Fixed Knots}
\description{
  Compute B-spline coefficients for regression quantile B-spline with
  stepwise knots selection and quantile B-spline with fixed knots
  \bold{regression spline}, using Ng (1996)'s algorithm.
}
\usage{
qbsks(x,y,w,pw, knots,nknots, degree,Tlambda, constraint,
      n.sub = n1000cut(n), equal,smaller, greater,gradient, coef,maxiter,
      trace, n.equal,n.smaller,n.greater,n.gradient,
      nrq,nl1, neqc, nj0, tau,lam,tmin,kmax,lstart,
      ks,mk.flag, knots.add, ic, print.mesg,
      factor, tol.kn = 1e-6, eps = .Machine$double.eps, print.warn)
}
\arguments{
  \item{x}{numeric vector, sorted increasingly, the abscissa values}
  \item{y}{numeric, same length as \code{x}, the observations.}
  \item{w}{numeric vector of weights, same length as \code{x}, as in
    \code{\link{cobs}}.}
  \item{pw}{penalty weights vector ... ...}
  \item{knots}{ ~~Describe \code{knots} here~~ }
  \item{nknots}{ ~~Describe \code{nknots} here~~ }
  \item{degree}{integer specifying polynomial degree; must be 1 or 2.}
  \item{Tlambda}{ ~~Describe \code{Tlambda} here~~ }
  \item{constraint}{ ~~Describe \code{constraint} here~~ }
  \item{n.sub}{integer, not larger than sample size \code{n}; the
    default has \code{n.sub == n} as long as \code{n} is less than 1000.}
  \item{equal}{ ~~Describe \code{equal} here~~ }
  \item{smaller}{ ~~Describe \code{smaller} here~~ }
  \item{greater}{ ~~Describe \code{greater} here~~ }
  \item{gradient}{ ~~Describe \code{gradient} here~~ }
  \item{coef}{ ~~Describe \code{coef} here~~ }
  \item{maxiter}{ ~~Describe \code{maxiter} here~~ }
  \item{trace}{ ~~Describe \code{trace} here~~ }
  \item{n.equal}{ ~~Describe \code{n.equal} here~~ }
  \item{n.smaller}{ ~~Describe \code{n.smaller} here~~ }
  \item{n.greater}{ ~~Describe \code{n.greater} here~~ }
  \item{n.gradient}{ ~~Describe \code{n.gradient} here~~ }
  \item{nrq}{ ~~Describe \code{nrq} here~~ }
  \item{nl1}{ ~~Describe \code{nl1} here~~ }
  \item{neqc}{ ~~Describe \code{neqc} here~~ }
  \item{nj0}{ ~~Describe \code{nj0} here~~ }
  \item{tau}{ ~~Describe \code{tau} here~~ }
  \item{lam}{ ~~Describe \code{lam} here~~ }
  \item{tmin}{ ~~Describe \code{tmin} here~~ }
  \item{kmax}{ ~~Describe \code{kmax} here~~ }
  \item{lstart}{ ~~Describe \code{lstart} here~~ }
  \item{ks}{ ~~Describe \code{ks} here~~ }
  \item{mk.flag}{ ~~Describe \code{mk.flag} here~~ }
  \item{knots.add}{ ~~Describe \code{knots.add} here~~ }
  \item{ic}{ ~~Describe \code{ic} here~~ }
  \item{print.mesg}{ ~~Describe \code{print.mesg} here~~ }
  \item{factor}{ ~~Describe \code{factor} here~~ }
  \item{tol.kn}{``tolerance'' for shifting the outer knots.}
  \item{eps}{tolerance passed to \code{\link{drqssbc}}.}
  \item{print.warn}{flag indicating if and how much warnings and
    information is to be printed;  currently just passed to \code{\link{drqssbc}}.}
}
\details{
  This is an auxiliary function for \code{\link{cobs}(*, lambda = 0)},
  possibly interesting on its own.  This documentation is currently sparse; read
  the source code!
}
\value{
  a list with components
  \item{coef}{..}
  \item{fidel}{..}
  \item{k}{dimensionality of model fit.}
  \item{ifl}{integer ``flag''; the return code.}
  \item{icyc}{integer of length 2, see \code{\link{cobs}}.}
  \item{knots}{the vector of inner knots.}
  \item{nknots}{the number of inner knots.}
  \item{nvar}{the number of ``variables'', i.e. unknowns including
    constraints.}
  \item{lambda}{the penalty factor, chosen or given.}
  \item{pseudo.x}{the pseudo design matrix \eqn{X}, as returned from
    \code{\link{drqssbc}}.}
}
\references{
  Ng, P. (1996)
  An Algorithm for Quantile Smoothing Splines,
  \emph{Computational Statistics \& Data Analysis} \bold{22}, 99--118.
}
\author{Pin Ng; this help page: Martin Maechler.}

\seealso{the main function \code{\link{cobs}}; further
  \code{\link{drqssbc}} which is called from \code{qbsks()}.
}
\examples{
}
\keyword{smooth}
\keyword{utilities}
