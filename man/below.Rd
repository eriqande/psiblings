\name{below}
\alias{below}
\title{return logical vector picking out the lower triangle of an n x n square matrix}
\usage{
below(n, strictly = F)
}
\arguments{
  \item{n}{dimension of a square matrix}

  \item{strictly}{TRUE means do not include the diagonal
  and FALSE means do.}
}
\description{
A function that delivers the same functionality as
lower.tri(matrix(1,n,n))
}
\note{
There is similar functionality in base R, but Hans calls
this in places so we include it here.
}
\author{
Hans J. Skaug \email{Hans.Skaug@math.uib.no}
}

