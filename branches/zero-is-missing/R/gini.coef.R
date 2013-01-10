`gini.coef` <- function(x, na.rm = TRUE)
{
    check.flags(na.rm)
    if (na.rm) {
        y <- as.double(x[!is.na(x)])
    } else {
        idx.good <- which(!is.na(x))
        n.good <- length(idx.good)
        if (n.good > 0) {
            y <- as.double(x[idx.good[1]:idx.good[n.good]])
        } else {
            y <- vector(mode="numeric", length=0)
        }
    }
    .C(dplR.gini,
       y, as.integer(length(y)), result=NaN, NAOK=TRUE, DUP=FALSE)$result
}
