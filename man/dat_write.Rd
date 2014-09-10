\name{dat_write}
\alias{dat_write}
\title{Write to .dat file for admb}
\usage{
dat_write(name, L)
}
\arguments{
  \item{name}{path of file to write to}

  \item{L}{list of values to write (ought to be named for
  clarity)}
}
\description{
Writes a named list to a .dat file that admb can read.
Basically just write the values of vectors or matrices to
the file and puts a comment above each telling which
parameter it is.
}
\author{
Hans J. Skaug \email{Hans.Skaug@math.uib.no}
}

