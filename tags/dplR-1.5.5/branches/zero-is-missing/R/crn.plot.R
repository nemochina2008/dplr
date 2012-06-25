`crn.plot` <- function(crn, add.spline=FALSE, nyrs=NULL, f=0.5, ...){
    if (!is.data.frame(crn)) {
        stop("'crn' must be a data.frame")
    }
    check.flags(add.spline)

    op <- par(no.readonly=TRUE) # Save par
    on.exit(par(op))            # Reset par on exit
    par(mar=c(3, 3, 3, 3), mgp=c(1.25, 0.25, 0), tcl=0.25)

    yr.vec <- as.numeric(row.names(crn))
    crn2 <- crn
    if (!all(diff(yr.vec) == 1)) {
        crn2 <- complete.rwl.df(crn)
        yr.vec <- as.numeric(row.names(crn2))
    }
    crn.names <- names(crn2)
    nCrn <- ncol(crn2)
    ## Check to see if the crn has sample depth
    sd.exist <- crn.names[nCrn] == "samp.depth"
    if (sd.exist) {
        samp.depth <- crn2[[nCrn]]
        nCrn <- nCrn-1
    }
    if (nCrn > 1) {
        layout(matrix(seq_len(nCrn), nrow=nCrn, ncol=1))
    }
    text.years <- gettext("Years", domain="R-dplR")
    text.rwi <- gettext("RWI", domain="R-dplR")
    text.samp <- gettext("Sample Depth", domain="R-dplR")
    nyrs2 <- nyrs
    for (i in seq_len(nCrn)) {
        this.crn <- crn2[[i]]
        plot(yr.vec, this.crn, type="l",
             xlab=text.years, ylab=text.rwi, main=crn.names[i], ...)
        idx.good <- which(!is.na(this.crn))
        n.good <- length(idx.good)
        if (n.good > 0) {
            this.crn <- this.crn[idx.good[1]:idx.good[n.good]]
        } else {
            this.crn <- numeric(0)
        }
        ## Only possibly NULL in the first round of the for loop
        if (is.null(nyrs2) || nyrs2 == 0) {
            nyrs2 <- length(this.crn) * 0.33
        }
        idx.good2 <- which(!is.na(this.crn))
        n.good2 <- length(idx.good2)
        if (n.good2 >= 3) {
            spl <- rep(NA_real_, length(yr.vec))
            spl[idx.good2 + (idx.good[1] - 1)] <-
                ffcsaps(y=this.crn[idx.good2], x=idx.good2, nyrs=nyrs2, f=f)
            if (add.spline) {
                lines(yr.vec, spl, col="red", lwd=2)
            }
        }
        abline(h=1)
        if (sd.exist) {
            par(new=TRUE)
            plot(yr.vec, samp.depth, type="l", lty="dashed",
                 xlab="", ylab="", axes=FALSE)
            axis(4, at=pretty(samp.depth))
            mtext(text.samp, side=4, line=1.25)
        }
    }
}
