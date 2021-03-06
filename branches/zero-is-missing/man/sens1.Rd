\name{sens1}
\alias{sens1}
\title{ Calculate Mean Sensitivity  }
\description{
  This function calculates mean sensitivity of a detrended ring-width
  series.
}
\usage{
  sens1(x, na.rm = FALSE)
}
\arguments{
  \item{x}{ a \code{numeric} vector containing the series }
  \item{na.rm}{ a \code{logical} flag.  If \code{TRUE}, the sensitivity
    \code{NA} values are removed.  See \sQuote{Details}.  If
    \code{FALSE}, \code{NA} values will propagate to the result. }
}
\details{
  This calculates mean sensitivity according to Eq. 1 in Biondi and
  Qeadan (2008).  This is the standard measure of sensitivity in
  dendrochronology and is typically calculated on detrended series.
  Users unfamiliar with sensitivity in tree-ring data should see Cook
  and Kairiukstis (1990) and Fritts (2001) for further details.

  Sensitivity is computed using terms involving each pair of two
  consecutive values in \code{\var{x}}.  If \code{\var{na.rm}} is
  \code{TRUE}, pairs with at least one \code{NA} value are dropped. The
  number of pairs used is taken into account.
}
\value{ the mean sensitivity.
}
\references{

  Biondi, F. and Qeadan, F. (2008) Inequality in Paleorecords.
  \emph{Ecology}, 89(4):1056\enc{–}{--}1067.

  Cook, E. R. and Kairiukstis, L. A. (1990) \emph{Methods of
  Dendrochronology: Applications in the Environmental Sciences}.
  Springer.  \acronym{ISBN-13}: 978-0-7923-0586-6.

  Fritts, H. C. (2001) \emph{Tree Rings and Climate}.  Blackburn.
  \acronym{ISBN-13}: 978-1-930665-39-2.
  
}
\author{ Mikko Korpela, based on original by Andy Bunn }
\seealso{ \code{\link{sens2}}, \code{\link{rwl.stats}} }
\examples{data(ca533)
ca533.rwi <- detrend(rwl = ca533, method = "ModNegExp")
sens1(ca533.rwi[, 1])
}
\keyword{ univar }
