\name{chron}
\alias{chron}
\title{Build Mean Value Chronology}
\description{
  This function builds a mean value chronology, typically from a
  \code{data.frame} of detrended ring widths as produced by
  \code{\link{detrend}}.
}
\usage{
chron(x, prefix = "xxx", biweight = TRUE, prewhiten = FALSE,
      ids = NULL, x.out = FALSE)
}
\arguments{
  \item{x}{a \code{data.frame} of ring widths with
    \code{rownames(\var{x})} containing years and \code{colnames(x)}
    containing each series \acronym{ID} such as produced by
    \code{\link{read.rwl}}}
  \item{prefix}{a \code{character} string with less than 4
    characters. Defaults to \code{"xxx"}}
  \item{biweight}{\code{logical} flag. If \code{TRUE} then a robust
    mean is calculated using \code{\link{tbrm}.}}
  \item{prewhiten}{\code{logical} flag. If \code{TRUE} each series is
    whitened using \code{\link{ar}} prior to averaging.}
  \item{ids}{ an optional \code{data.frame} with column one named
    \code{"tree"} giving a \code{numeric} \acronym{ID} for each tree and
    column two named \code{"core"} giving a \code{numeric} \acronym{ID}
    for each core. Defaults to one core per tree as
    \code{data.frame(tree = 1:ncol(\var{x}), core = rep(1,
      ncol(\var{x})))}. }
  \item{x.out}{ a \code{logical} flag.  If \code{TRUE}, one
    (\code{\var{prewhiten}==FALSE}) or two
    (\code{\var{prewhiten}==TRUE}) data.frames are included in the
    returned list in addition to its usual components.  See
    \sQuote{Value}. }
}
\details{
  This either averages the rows of the \code{data.frame} using a mean or
  a robust mean (the so-called standard chronology) or can do so from
  the residuals of an ar process (the residual chronology).
}
\value{
  A \code{data.frame} with the standard chronology \code{"std"},
  residual chronology \code{"res"} (if
  prewhitening was performed), and the sample depths for both standard
  and residual chronologies (\code{"samps"} and optionally
  \code{"samps.res"}, respectively).  If \code{\var{x.out}} is
  \code{TRUE}, the result prior to the final averaging step is also
  returned (\code{"x.std"} and optionally \code{"x.res"}).
}
\references{
  Cook, E. R. and Kairiukstis, L. A. (1990) \emph{Methods of
    Dendrochronology: Applications in the Environmental Sciences.}
  Springer.  \acronym{ISBN-13}: 978-0-7923-0586-6.

  Fritts, H. C. (2001) \emph{Tree Rings and Climate.}  Blackburn.
  \acronym{ISBN-13}: 978-1-930665-39-2.
}
\author{ Andy Bunn. Patched and improved by Mikko Korpela. }
\seealso{ \code{\link{read.rwl}}, \code{\link{detrend}},
  \code{\link{ar}}, \code{\link{chron}}
}
\examples{data(ca533)
ca533.rwi <- detrend(rwl = ca533, method = "ModNegExp")
ca533.crn <- chron(ca533.rwi, prefix = "CAM")
## With residual chron
ca533.crn <- chron(ca533.rwi, prefix = "CAM", prewhiten = TRUE)
}
\keyword{ manip }
