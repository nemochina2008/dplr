\name{read.compact}
\alias{read.compact}
\title{ Read \acronym{DPL} Compact Format Ring Width File }
\description{  
  This function reads in a \acronym{DPL} compact format file of ring
  widths.
}
\usage{
read.compact(fname, zero.as.na = TRUE)
}
\arguments{
  \item{fname}{ a \code{character} vector giving the file name of the
    rwl file. }
  \item{zero.as.na}{ a \code{logical} flag.  If \code{TRUE} (the
    default), zeros in the input are interpreted as missing values
    (\code{NA}). }
}
\details{
  This function should be able to read files written by the
  Dendrochronology Program Library (\acronym{DPL}) in its compact
  format.
}
\value{
  A \code{data.frame} with the series in columns and the years as
  rows. The series \acronym{ID}s are the column names and the years are
  the row names.
}
\author{ Mikko Korpela }
\seealso{ \code{\link{read.rwl}}, \code{\link{read.tucson}},
  \code{\link{read.tridas}}, \code{\link{read.fh}},
  \code{\link{write.compact}}}
\keyword{ IO }
