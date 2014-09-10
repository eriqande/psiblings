\name{pin_write}
\alias{pin_write}
\title{Write to .pin file for admb}
\usage{
pin_write(name, L)
}
\arguments{
  \item{name}{path of file to write to}

  \item{L}{list of values to write (ought to be named for
  clarity)}
}
\description{
Writes a named list to a .pin file that admb can read.
Probably just write the values of vectors or matrices to
the file and puts a comment above each telling which
parameter it is.
}
\author{
Hans J. Skaug \email{Hans.Skaug@math.uib.no}
}

