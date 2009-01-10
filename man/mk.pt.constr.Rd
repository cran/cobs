\name{mk.pt.constr}
\alias{mk.pt.constr}
\title{COBS auxiliary for constructing pointwise constraint specifications}
\description{
  COBS (\code{\link{cobs}}) auxiliary function for constructing the
  pointwise constraint specification list from the \code{pointwise}
  3-column matrix (as used as argument in \code{\link{cobs}}).
}
\usage{
mk.pt.constr(pointwise)
}
\arguments{
  \item{pointwise}{numeric 3-column matrix, see \code{pointwise} in \code{\link{cobs}}.}
}
\value{
  A list with components
  \item{n.equal}{number of equality constraints}
  \item{n.greater}{number of ``greater'' constraints}
  \item{n.smaller}{number of ``smaller'' constraints}
  \item{n.gradient}{number of gradient constraints}

  Unless the input \code{pointwise} was \code{NULL}, the result also
  has corresponding components:
  \item{equal}{3-column matrix ofequality constraints}
  \item{greater}{3-column matrix of``greater'' constraints}
  \item{smaller}{3-column matrix of``smaller'' constraints}
  \item{gradient}{3-column matrix ofgradient constraints}
}
\author{Martin Maechler}
\examples{
 ## from  ?cobs:
 x <- seq(-1,3,,150)
 con <- rbind(c( 1,min(x),0), # f(min(x)) >= 0
              c(-1,max(x),1), # f(max(x)) <= 1
              c(0,  0,   0.5))# f(0)      = 0.5
 r <- mk.pt.constr(con)
 str(r)
}
\keyword{misc}

