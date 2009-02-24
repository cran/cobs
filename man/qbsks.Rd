\name{qbsks2}
\alias{qbsks2}
\title{Quantile B-Spline with Fixed Knots}
\description{
  Compute B-spline coefficients for regression quantile B-spline with
  stepwise knots selection and quantile B-spline with fixed knots
  \bold{regression spline}, using Ng (1996)'s algorithm.
}
\usage{
qbsks2(x,y,w,pw, knots,nknots, degree,Tlambda, constraint, ptConstr,
      maxiter, trace, nrq,nl1, neqc, tau, select.lambda,
      ks, do.select, knots.add, repeat.delete.add, ic, print.mesg,
      give.pseudo.x = TRUE,
      rq.tol = 1e-8, tol.kn = 1e-6, tol.0res = 1e-6, print.warn, nk.start)
}
\arguments{
  \item{x}{numeric vector, sorted increasingly, the abscissa values}
  \item{y}{numeric, same length as \code{x}, the observations.}
  \item{w}{numeric vector of weights, same length as \code{x}, as in
    \code{\link{cobs}}.}
  \item{pw}{penalty weights vector ... ...}%% << FIXME, currently unused
  \item{knots}{numeric vector of knots of which \code{nknots} will be used.}
  \item{nknots}{number of \code{knots} to be used.}% FIXME: Say more in \details
  \item{degree}{integer specifying polynomial degree; must be 1 or 2.}
  \item{Tlambda}{ ~~Describe \code{Tlambda} here~~ }
  \item{constraint}{string (or empty) specifying the global constraints;
    see \code{\link{cobs}}.}
  \item{ptConstr}{\code{\link{list}} of \bold{p}oin\bold{t}wise constraints.}
  \item{maxiter}{ ~~Describe \code{maxiter} here~~ }
  \item{trace}{integer or logical indicating the tracing level of the
    \emph{underlying} algorithms; not implemented (due to lack of trace
    in quantreg ...)}
  \item{nrq,nl1,neqc}{integers specifying dimensionalities, directly
    passed to \code{\link{drqssbc2}}, see there.}
  \item{tau}{desired quantile level (in interval \eqn{(0,1)}).}
  \item{select.lambda}{ ~~Describe \code{lam} here~~ }
  \item{ks}{number used as offset in SIC/AIC/BIC.}
  \item{do.select}{logical indicating if knots shall be selected
    (instead of used as specified).}
  \item{knots.add, repeat.delete.add}{logicals, see \code{\link{cobs}}.}
  \item{ic}{information criterion to use, see \code{\link{cobs}}.}
  \item{print.mesg}{an integer indicating how \code{qbsks2()} should
    print message about its current stages.}
  \item{give.pseudo.x}{logical indicating if the pseudo design matrix
    \eqn{\tilde{X}}{X\~} should be returned (as \emph{sparse} matrix).}
  \item{rq.tol}{numeric convergence tolerance for the interior point
    algorithm called from \code{\link[quantreg]{rq.fit.sfnc}()} or
    \code{\link[quantreg]{rq.fit.sfn}()}.}
  \item{tol.kn}{\dQuote{tolerance} for shifting the outer knots.}
  \item{tol.0res}{tolerance passed to \code{\link{drqssbc2}}.}
  \item{print.warn}{flag indicating if and how much warnings and
    information is to be printed; currently just passed to
    \code{\link{drqssbc2}}.}
  \item{nk.start}{number of starting knots used in automatic knot selection.}
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
    \code{\link{drqssbc2}}.}
}
\references{
  Ng, P. (1996)
  An Algorithm for Quantile Smoothing Splines,
  \emph{Computational Statistics \& Data Analysis} \bold{22}, 99--118.

  See also the \emph{references} in \code{\link{cobs}}.
}
\author{Pin Ng; this help page: Martin Maechler.}

\seealso{the main function \code{\link{cobs}}; further
  \code{\link{drqssbc2}} which is called from \code{qbsks2()}.
}
% \examples{
% }
\keyword{smooth}
\keyword{utilities}
