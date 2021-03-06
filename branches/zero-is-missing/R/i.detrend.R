`i.detrend` <- function(rwl, y.name=names(rwl), nyrs = NULL, f = 0.5,
                        pos.slope = FALSE)
{
    if (!is.data.frame(rwl)) {
        stop("'rwl' must be a data.frame")
    }
    out <- rwl
    n.col <- ncol(rwl)
    fmt <- gettext("Detrend series %d of %d\n", domain="R-dplR")
    for (i in seq_len(n.col)) {
        cat(sprintf(fmt, i, n.col))
        fits <- i.detrend.series(rwl[[i]], y.name=y.name[i], nyrs = nyrs,
                                 f = f, pos.slope = pos.slope)
        out[, i] <- fits
    }
    out
}
