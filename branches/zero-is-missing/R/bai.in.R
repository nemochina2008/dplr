bai.in <- function(rwl, d2pith = NULL, warn.na = TRUE) {

    check.flags(warn.na)
    if(!is.data.frame(rwl)) {
        stop("'rwl' must be a data.frame")
    }
    if(!is.null(d2pith)) {
        if(ncol(rwl) != nrow(d2pith))
            stop("dimension problem: ", "'ncol(rw)' != 'nrow(d2pith)'")
        if(!all(d2pith[, 1] %in% names(rwl)))
            stop("series ids in 'd2pith' and 'rwl' do not match")
        d2pith.vec <- d2pith[, 2]
    } else {
        ## distance offset if not given
        d2pith.vec <- rep(0, ncol(rwl))
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
            ## get ring area
            bai <- pi*dat*(dat+2*(cumsum(dat) + d2pith.vec[i] - dat))
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
