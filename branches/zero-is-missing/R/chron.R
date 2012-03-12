`chron` <-
    function(x, prefix="xxx", biweight=TRUE, prewhiten=FALSE)
{
    if(!is.data.frame(x)) {
        stop("'x' must be a data.frame")
    }
    check.flags(prewhiten, biweight)
    prefix.str <- as.character(prefix)
    if (length(prefix.str) != 1 || nchar(prefix.str) > 3) {
        stop("'prefix' must be a character string with less than 4 characters")
    }
    x2 <- x
    if (!all(diff(as.numeric(row.names(x))) == 1)) {
        x2 <- complete.rwl.df(x)
    }
    samps <- rowSums(!is.na(x2))
    if (!biweight) {
        std <- rowMeans(x2, na.rm=TRUE)
    } else {
        std <- apply(x2, 1, tbrm, C=9)
    }
    if (prewhiten) {
        x.ar <- apply(x2, 2, ar.func)
        if (!biweight) {
            res <- rowMeans(x.ar, na.rm=TRUE)
        } else {
            res <- apply(x.ar, 1, tbrm, C=9)
        }
        res[is.nan(res)] <- NA
        out <- data.frame(std, res, samps)
        names(out) <- c(paste0(prefix.str, "std"),
                        paste0(prefix.str, "res"),
                        "samp.depth")
    } else {
        out <- data.frame(std, samps)
        names(out) <- c(paste0(prefix.str, "std"), "samp.depth")
    }
    row.names(out) <- row.names(x2)
    out
}
