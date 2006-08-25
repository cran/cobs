\name{drqssbc2}
\alias{drqssbc2}
\title{Regression Quantile Smoothing Spline with Constraints}
\description{
  Estimate the B-spline coefficients for a regression quantile
  \emph{smoothing} spline with optional constraints, using Ng(1996)'s
  algorithm.
}
\usage{
drqssbc2(x, y, w = rep.int(1,n), pw, knots, degree, Tlambda,
        constraint, ptConstr, maxiter = 100, trace = 0,
        nrq = length(x), nl1, neqc, niqc, nvar,
        tau = 0.5, select.lambda, give.pseudo.x = FALSE,
        rq.tol = 1e-8 * sc.y, tol.0res = 1e-6,
        print.warn = TRUE, rq.print.warn = FALSE)
}
\arguments{
  \item{x}{numeric vector, sorted increasingly, the abscissa values}
  \item{y}{numeric, same length as \code{x}, the observations.}
  \item{w}{numeric vector of weights, same length as \code{x}, as in
    \code{\link{cobs}}.}
  \item{pw}{penalty weights vector passed to \code{\link{l1.design2}} or
    \code{\link{loo.design2}}. \bold{FIXME}: This is currently unused.}
  \item{knots}{numeric vector of knots for the splines.}
  \item{degree}{integer, must be 1 or 2.}
  \item{Tlambda}{vector of smoothing parameter values \eqn{\lambda}; if
    it is longer than one, an \dQuote{optimal} value will be selected
    from these.}
  \item{constraint}{see \code{\link{cobs}} (but cannot be abbreviated here).}
  \item{ptConstr}{\code{\link{list}} of \bold{p}oin\bold{t}wise constraints;
    notably \code{equal}, \code{smaller}, \code{greater} and \code{gradient}
    are 3-column matrices specifying the respective constraints.  May have
    0 rows if there are no constraints of the corresponding kind.}
  \item{maxiter}{maximal number of iterations; defaults to 100.}
  \item{trace}{integer or logical indicating the tracing level of the
    \emph{underlying} algorithms; not much implemented (due to lack of trace
    in quantreg ...)}
  \item{nrq}{ ~~Describe \code{nrq} here~~ }
  \item{nl1}{ ~~Describe \code{nl1} here~~ }
  \item{neqc}{integer giving the number of equations.}
  \item{niqc}{integer giving the number of \bold{i}ne\bold{q}uality
    \bold{c}onstraints; of the same length as \code{constraint}.}
  \item{nvar}{integer giving the number of equations \emph{and} constraints.}
  \item{tau}{desired quantile level; defaults to 0.5 (median).}
  \item{select.lambda}{logical indicating if an optimal lambda should be
    selected from the vector of \code{Tlambda}.}
  \item{give.pseudo.x}{logical indicating if the pseudo design matrix
    \eqn{\tilde{X}}{X\~} should be returned (as \emph{sparse} matrix).}
  \item{rq.tol}{numeric convergence tolerance for the interior point
    algorithm called from \code{\link[quantreg]{rq.fit.sfnc}()} or
    \code{\link[quantreg]{rq.fit.sfn}()}.  Note that (for scale
    invariance) this has to be in units of \code{y}, which the default
    makes use of.}
  \item{tol.0res}{tolerance used to check for zero residuals, i.e.,
    \eqn{|r_i| < tol * mean(|r_i|)}.}
  \item{print.warn}{logical indicating if warnings should be printed,
    when the algorithm seems to have behaved somewhat unexpectedly.}
  \item{rq.print.warn}{logical indicating if warnings should be printed
    from inside the \code{rq.*} function calls, see below.}
}
\details{
  This is an auxiliary function for \code{\link{cobs}}, possibly
  interesting on its own.  Depending on \code{degree}, either
  \code{\link{l1.design2}} or \code{\link{loo.design2}} are
  called for construction of the sparse design matrix.

  Subsequently, either \code{\link[quantreg]{rq.fit.sfnc}} or
  \code{\link[quantreg]{rq.fit.sfn}} is called as the main \dQuote{work horse}.

  This documentation is currently sparse; read the source code!
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
  The main function \code{\link{cobs}} and its auxiliary
  \code{\link{qbsks2}} which calls \code{drqssbc2()} repeatedly.

  \code{\link{l1.design2}} and \code{\link{loo.design2}};
  further \code{\link[quantreg]{rq.fit.sfnc}} and
  \code{\link[quantreg]{rq.fit.sfn}} from package \pkg{quantreg}.
}
\examples{
set.seed(1243)
x  <- 1:32
fx <- (x-5)*(x-15)^2*(x-21)
y  <- fx + round(rnorm(x,s = 0.25),2)
%% FAILS --- FIXME
%% drqssbc2(x,y,nrq=32,lam=1,degree=1,knots=c(1,5,15,32))
}
\keyword{smooth}
\keyword{utilities}
