\name{par_read}
\alias{par_read}
\title{Read from a .par file produced by admb}
\usage{
par_read(name, ncols = list())
}
\arguments{
  \item{name}{path of file to read from}

  \item{ncols}{named list with the number of columns for
  matrix parameters}
}
\description{
Read from a .par file produced by admb or a file of
identical structure. This reads the values of the objective
function and the named parameters. If some of the
parameters are matrices, then their number of columns
should be supplied in a named list. (See example)
}
\author{
Hans J. Skaug \email{Hans.Skaug@math.uib.no}
}

