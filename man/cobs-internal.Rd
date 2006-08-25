\name{scobs-internal}
%%-- maybe more of these really should be documented and ``exported''
\alias{getdim2}
\alias{l1.design2}
\alias{loo.design2}
\alias{dn}% two simple bandwidth rules
\alias{shat}% sparsity estimate via resid.
\alias{.splBasis}
\alias{.splValue}
% The following should be eliminated (MM thinks)
%\alias{mesg}
%\alias{warn}
\title{Internal COBS functions}
\description{
  Internal scobs functions.
}
\usage{
dn(p, n, hs = FALSE, alpha)

getdim2(degree, nknots, constraint)

l1.design2(x, w, constraint, ptConstr,
          knots, pw, nrq, nl1, neqc, niqc, nvar, lambda)
loo.design2(x, w, constraint, ptConstr,
          knots, pw, nrq, nl1, neqc, niqc, nvar, lambda)

shat(residual, tau, alpha, hs)

.splValue(degree, knots, coef, xo)
.splBasis(ord, knots, ncoef, xo, derivs)

%mesg(number, \dots)
% warn(number, \dots)
}
\details{
  These are not (yet?) to be called by the user and have not been
  documented by the original COBS authors.
}
\keyword{internal}
