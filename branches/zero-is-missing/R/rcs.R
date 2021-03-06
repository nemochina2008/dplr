rcs <- function(rwl, po, nyrs=NULL, f=0.5, biweight=TRUE, ratios=TRUE,
                rc.out=FALSE, make.plot=TRUE, ...) {
    if (!is.data.frame(rwl)) {
        stop("'rwl' must be a data.frame")
    }
    n.col <- ncol(rwl)
    if (n.col == 0) {
        return(rwl)
    }
    if (n.col != nrow(po)) {
        stop("dimension problem: ", "'ncol(rw)' != 'nrow(po)'")
    }
    col.names <- names(rwl)
    if (is.null(col.names) || anyDuplicated(col.names) ||
        any(is.na(col.names))) {
        stop("'rwl' must have unique, non-NA names")
    }
    if (!all(sort(po[, 1]) == sort(col.names))) {
        stop("series ids in 'po' and 'rwl' do not match")
    }
    if (any(po[, 2] < 1)) {
        stop("minimum 'po' is 1")
    }
    if (!all(is.int(po[, 2]))) {
        stop("each value in 'po' must be an integer")
    }
    check.flags(make.plot, rc.out, ratios, biweight)

    seq.cols <- seq_len(n.col)
    rwl2 <- rwl
    rownames(rwl2) <- rownames(rwl2) # guard against NULL names funniness
    yrs <- as.numeric(row.names(rwl2))
    stopifnot(diff(yrs) == 1)

    rwl.ord <- apply(rwl2, 2, sortByIndex)
    rwca <- matrix(NA,
                   ncol = n.col,
                   nrow = nrow(rwl.ord) + max(po[, 2]))
    nrow.m1 <- nrow(rwl.ord) - 1
    yrs2pith <- po[match(col.names, po[, 1]), 2]
    for (i in seq.cols) {
        rwca[yrs2pith[i]:(yrs2pith[i] + nrow.m1), i] <- rwl.ord[, i]
    }

    if (biweight) {
        ca.m <- apply(rwca, 1, tbrm, C = 9)
    } else {
        ca.m <- rowMeans(rwca, na.rm=TRUE)
    }

    ## spline follows B&Q 2008 as 10% of the RC length
    good.idx <- which(!is.na(ca.m))
    if (is.null(nyrs)) {
        n.good <- length(good.idx)
        if (n.good >= 1) {
            rc.extent <- good.idx[n.good] - good.idx[1] + 1
        } else {
            rc.extent <- 0
        }
        if (rc.extent < 10) {
            nyrs2 <- 0
        } else {
            nyrs2 <- floor(rc.extent * 0.1)
        }
    } else {
        nyrs2 <- nyrs
    }
    tmp <- ffcsaps(y=ca.m[good.idx], x=good.idx, nyrs=nyrs2, f=f)
    rc <- rep(NA, nrow(rwca))
    rc[good.idx] <- tmp
    ## calculate indices as ratios or differences
    if (ratios) {
        rwica <- rwca / rc
    } else {
        rwica <- rwca - rc
    }
    ## and restore to cal years
    rwi <- rwl2
    for (i in seq.cols) {
        series.yrs <- yr.range(rwl2[[i]], yr.vec=yrs)
        first <- series.yrs[1]
        last <- series.yrs[2]
        seq.yrs <- seq_len(last - first + 1)
        rwi[[i]][seq.yrs + (first - yrs[1])] <-
            rwica[seq.yrs + (yrs2pith[i] - 1), i]
    }
    if (make.plot) {
        par(mar = c(4, 4, 4, 4) + 0.1, mgp = c(1.25, 0.25, 0), tcl = 0.25)
        plot(rwca[, 1], ylim=range(rwca, na.rm=TRUE), type="n", ylab="mm",
             xlab=gettext("Cambial Age (Years)", domain="R-dplR"), ...)
        for (i in seq.cols) {
            lines(rwca[, i], col="grey")
        }
        lines(ca.m, lwd=1.5, col="black")
        lines(rc, lwd=2, col="red")
    }
    if (rc.out) {
        list(rwi=rwi, rc=rc)
    } else {
        rwi
    }
}
