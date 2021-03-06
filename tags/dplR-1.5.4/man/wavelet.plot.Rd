\name{wavelet.plot}
\alias{wavelet.plot}
\title{ Plot a Continuous Wavelet Transform }
\description{
  This function creates a \code{filled.contour} plot of a continuous
  wavelet transform as output from \code{\link{morlet}}.
}
\usage{
wavelet.plot(wave.list,
             wavelet.levels = quantile(wave.list$Power,
                                       probs = seq(from=0, to=1, by=0.1)),
             add.coi = TRUE, add.sig = TRUE, x.lab = gettext("Time"),
             period.lab = gettext("Period"), crn.lab = gettext("RWI"),
             key.cols = rev(rainbow(length(wavelet.levels)-1)),
             key.lab = parse(text=paste0("\"", gettext("Power"), "\"^2")),
             add.spline = FALSE, f = 0.5, nyrs = NULL,
             crn.col = "black", crn.lwd = 1,
             crn.ylim = range(wave.list$y)*1.1, side.by.side = FALSE)
}
\arguments{
  \item{wave.list}{A \code{list}. Output from \code{\link{morlet}}.}
  \item{wavelet.levels}{A \code{numeric} vector. Values for levels of
    the filled contours for the wavelet plot.}
  \item{add.coi}{A \code{logical} flag. Add cone of influence?}
  \item{add.sig}{A \code{logical} flag. Add contour lines for significance?}
  \item{x.lab}{X-axis label.}
  \item{period.lab}{Y-axis label for the wavelet plot.}
  \item{crn.lab}{Y-axis label for the time-series plot.}
  \item{key.cols}{A vector of colors for the wavelets and the key.}
  \item{key.lab}{Label for key.}
  \item{add.spline}{A \code{logical} flag. Add a spline to the
    time-series plot using \code{\link{ffcsaps}}?. }
  \item{nyrs}{ A number giving the rigidity of the smoothing spline, defaults to
    0.33 of series length if nyrs is \code{NULL}. }
  \item{f}{ A number between 0 and 1 giving the frequency response or
    wavelength cutoff for the smoothing spline. Defaults to 0.5. }
  \item{crn.col}{Line color for the time-series plot.}
  \item{crn.lwd}{Line width for the time-series plot.}
  \item{crn.ylim}{Axis limits for the time-series plot.}
  \item{side.by.side}{A \code{logical} flag. Plots will be in one row if
    \code{TRUE}. }
}
\details{
  This produces a plot of a continuous wavelet transform and plots the
  original time series. Contours are added for significance and a cone of
  influence polygon can be added as well. Anything within the cone of
  influence should not be interpreted.

  The time series can be plotted with a smoothing spline as well.
}
\value{
  None. This function is invoked for its side effect, which is to produce a
  plot.
}
\references{
  Torrence, C. and Compo, G. P. (1998) A practical guide to wavelet
  analysis.  \emph{Bulletin of the American Meteorological Society},
  79:61--78.
}
\note{ The function \code{morlet} is a port of Torrence's
  \acronym{IDL} code available at
  \url{http://atoc.colorado.edu/research/wavelets/software.html}
}
\author{ Andy Bunn. Patched and improved by Mikko Korpela. }
\seealso{ \code{\link{morlet}}, \code{\link{ffcsaps}}}
\examples{data(ca533)
ca533.rwi <- detrend(rwl = ca533, method = "ModNegExp")
ca533.crn <- chron(ca533.rwi, prefix = "CAM", prewhiten = FALSE)
Years <- as.numeric(rownames(ca533.crn))
CAMstd <- ca533.crn[, 1]
out.wave <- morlet(y1 = CAMstd, x1 = Years, p2 = 9, dj = 0.1,
                   siglvl = 0.99)
wavelet.plot(out.wave)
levs <- quantile(out.wave$Power, probs = c(0, 0.5, 0.75, 0.9, 0.99))
wavelet.plot(out.wave, wavelet.levels = levs, add.sig = FALSE,
             key.cols = c("white", "green", "blue", "red"))
}
\keyword{ hplot }
