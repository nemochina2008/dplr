`sens1` <- function(x, na.rm = FALSE)
{
    check.flags(na.rm)
    idx.good <- which(!is.na(x))
    n.good <- length(idx.good)
    if (n.good > 0) {
        first.good <- idx.good[1]
        last.good <- idx.good[n.good]
        y <- as.double(x[first.good:last.good])
    } else {
        y <- vector(mode="numeric", length=0)
    }
    .C(dplR.sens1,
       y, as.integer(length(y)), as.integer(na.rm), result=NaN,
       NAOK=TRUE, DUP=FALSE)$result
}
