\name{cobs-internal}
\title{Internal COBS functions}
%      -----------------------
%%-- maybe more of these really should be documented and ``exported''
\alias{getdim}
\alias{l1.design}
\alias{loo.design}
\alias{dn}% two simple bandwidth rules
\alias{shat}% sparsity estimate via resid.
\alias{.splBasis}
\alias{.splValue}

\alias{cut00}
%
\description{
  Internal cobs functions.
}
\usage{
dn(p, n, hs = FALSE, alpha)

getdim(degree, nknots, constraint)

l1.design(x, w, constraint, equal, smaller, greater, gradient,
          knots, pw, n.equal, n.smaller, n.greater, n.gradient, nrq,
          nl1, neqc, niqc, nvar, lambda)
loo.design(x, w, constraint, equal, smaller, greater, gradient,
          knots, pw, n.equal, n.smaller, n.greater, n.gradient, nrq,
          nl1, neqc, niqc, nvar, lambda)

shat(residual, tau, alpha, hs)

.splValue(degree, knots, coef, xo)
.splBasis(ord, knots, ncoef, xo, derivs)
}
\details{
  These are not (yet?) to be called by the user and have not been
  documented by the original COBS authors.
}
\keyword{internal}
