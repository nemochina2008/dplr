\name{rwl.stats}
\alias{rwl.stats}
\title{ Calculate Descriptive Statistics on Ring-Width Series }
\description{
  This function calculates descriptive statistics on a \code{data.frame}
  of raw or detrended ring-width series.
}
\usage{
rwl.stats(rwl, na.rm.or.pass = TRUE, round.decimals = 3)
}
\arguments{
  \item{rwl}{ a \code{data.frame} with (usually) raw ring-width series
    as columns and years as rows such as that produced by
    \code{\link{read.rwl}}.  It is sometimes desirable to run this on
    detrended (e.g., rwi) data particularly for the sensitivity measures
    (\code{\link{sens1}}, \code{\link{sens2}}, \code{\link{gini.coef}}). }

  \item{na.rm.or.pass}{ a \code{logical} flag.  This is used as the
    \code{var{na.rm}} argument for various functions computing the
    statistics, and if \code{TRUE}, the first order autocorrelation is
    computed with \code{\var{na.action}=na.pass} (see \code{\link{acf}})
    instead of \code{\var{na.action}=na.fail} used otherwise. }

  \item{round.decimals}{ a \code{numeric} value giving the number of
    decimals the results will be rounded to.  A negative value can be
    used to disable rounding. }
}
\details{
  This calculates a variety of descriptive statistics commonly used in
  dendrochronology (see below). Users unfamiliar with these should see
  Cook and Kairiukstis (1990) and Fritts (2001) for further details.
}
\value{
  A \code{data.frame} containing descriptive stats on each
  \code{"series"}. These are the first and last year of the series as
  well as the length of the series (\code{"first"}, \code{"last"},
  \code{"year"}).  The mean, median, standard deviation are given
  (\code{"mean"}, \code{"median"}, \code{"stdev"}) as are the skewness,
  two measures of sensitivity, the Gini coefficient, and first order
  autocorrelation (\code{"skew"}, \code{"\link{sens1}"},
  \code{"\link{sens2}"}, \code{"\link{gini.coef}"}, \code{"ar1"}).
}
\references{

  Cook, E. R. and Kairiukstis, L.A. (1990) \emph{Methods of
  Dendrochronology: Applications in the Environmental Sciences}.
  Springer.  \acronym{ISBN-13}: 978-0-7923-0586-6.

  Fritts, H. C. (2001) \emph{Tree Rings and Climate}.  Blackburn.
  \acronym{ISBN-13}: 978-1-930665-39-2.

}
\author{ Andy Bunn. Slightly improved by Mikko Korpela. }
\seealso{ \code{\link{rwi.stats}}, \code{\link{read.rwl}} }
\examples{data(ca533)
rwl.stats(ca533)
}
\keyword{ misc }
