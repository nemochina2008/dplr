bai.out <- function(rwl, diam = NULL, warn.na = TRUE) {

    check.flags(warn.na)
    if(!is.data.frame(rwl)) {
        stop("'rwl' must be a data.frame")
    }
    if(!is.null(diam)) {
        if(ncol(rwl) != nrow(diam))
            stop("dimension problem: ", "'ncol(rw)' != 'nrow(diam)'")
        if(!all(diam[, 1] %in% names(rwl)))
            stop("series ids in 'diam' and 'rwl' do not match")
        diam.vec <- diam[, 2]
    }
    rwl2 <- rwl
    if (!all(diff(as.numeric(row.names(rwl))) == 1)) {
        rwl2 <- complete.rwl.df(rwl)
    }

    out <- rwl2
    ## vector of years
    n.vec <- seq_len(nrow(rwl2))
    for(i in seq_len(ncol(rwl2))){
        ## series to work with
        dat <- rwl2[[i]]
        ## strip out data from NA
        idx.good <- which(!is.na(dat))
        n.good <- length(idx.good)
        if (n.good > 0) {
            first.good <- idx.good[1]
            last.good <- idx.good[n.good]
            idx.seq <- first.good:last.good
            dat <- dat[idx.seq]
            ## get diameter if not given
            if (is.null(diam)) {
                d <- sum(dat)*2
            } else {
                d <- diam.vec[i]
            }
            ## get ring area
            r0 <- d/2 - c(0, cumsum(rev(dat)))
            bai <- -pi*rev(diff(r0*r0))
            if (warn.na && any(is.na(bai))) {
                warning(gettextf("NA values in series %s",
                                 names(rwl2)[i]), domain=NA)
            }
            ## write result
            out[idx.seq, i] <- bai
        }
    }
    ## return result
    out
}
