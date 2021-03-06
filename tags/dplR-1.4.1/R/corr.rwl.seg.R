corr.rwl.seg <- function(rwl, seg.length=50, bin.floor=100, n=NULL,
                         prewhiten = TRUE, pcrit=0.05, biweight=TRUE,
                         make.plot = TRUE, label.cex=1, ...){


    ## helper function
    yr.range <- function(x) {
        if(any(mask <- !is.na(x))) range(as.numeric(names(x))[mask])
        else c(NA, NA)
    }



    ## run error checks
    qa.xdate(rwl, seg.length, n, bin.floor)

    ## turn off warnings for this function
    ## The sig test for spearman's rho often produces warnings.
    w <- options("warn")
    on.exit(options(w))
    options(warn = -1)

    rnames <- rownames(rwl)
    cnames <- colnames(rwl)
    seg.lag <- seg.length/2
    nseries <- ncol(rwl)
    if(nseries < 2) stop("At least 2 series are needed in rwl")
    yrs <- as.numeric(rnames)
    nyrs <- length(yrs)
    min.yr <- min(yrs)
    max.yr <- max(yrs)
    if(is.null(bin.floor) || bin.floor == 0) min.bin <- min.yr
    else min.bin <- ceiling(min.yr/bin.floor)*bin.floor
    bins <- seq(from=min.bin, to=max.yr-seg.length+1, by=seg.lag)
    bins <- cbind(bins, bins+(seg.length-1))
    nbins <- nrow(bins)
    bin.names <- paste(bins[, 1], ".", bins[, 2], sep="")
    ## structures for results
    res.cor <- matrix(NA, nseries, nbins)
    rownames(res.cor) <- cnames
    colnames(res.cor) <- bin.names

    res.pval <- matrix(NA, nseries, nbins)
    rownames(res.pval) <- cnames
    colnames(res.pval) <- bin.names

    overall.cor <- matrix(NA, nseries, 2)
    rownames(overall.cor) <- cnames
    colnames(overall.cor) <- c("rho", "p-val")

    segavg.cor <- rep(NA, nbins)
    names(segavg.cor) <- bin.names

    ## normalize all series
    norm.one <- normalize1(rwl, n, prewhiten)
    ## rwi for segments altered by normalizing
    rwi <- norm.one$master
    idx.good <- norm.one$idx.good

    ## loop through series
    for(i in 1:nseries){
        idx.noti <- rep(TRUE, nseries)
        idx.noti[i] <- FALSE
        master.norm <- rwi[, idx.good & idx.noti, drop=FALSE]

        ## compute master series by normal mean or robust mean
        master <- vector(mode="numeric", length=nyrs)
        if (!biweight){
            for (j in 1:nyrs)
                master[j] <- exactmean(master.norm[j, ])
        } else {
            ## surprisingly, for loop is faster than apply
            for (j in 1:nyrs)
                master[j] <- tbrm(master.norm[j, ], C=9)
        }
        series <- rwi[, i]
        ## loop through bins
        for(j in 1:nbins){
            mask <- yrs%in%seq(from=bins[j, 1], to=bins[j, 2])
            ## cor is NA if there is not complete overlap
            if(!any(mask) ||
               any(is.na(series[mask])) ||
               any(is.na(master[mask]))){
                bin.cor <- NA
                bin.pval <- NA
            }
            else {
                tmp <- cor.test(series[mask], master[mask],
                                method = "spearman", alternative = "g")
                bin.cor <- tmp$estimate
                bin.pval <- tmp$p.val
            }
            res.cor[i, j] <- bin.cor
            res.pval[i, j] <- bin.pval
        }
        ## overall correlation
        tmp <- cor.test(series, master, method = "spearman", alternative = "g")
        overall.cor[i, 1] <- tmp$estimate
        overall.cor[i, 2] <- tmp$p.val
    }
    ## avg seg correlation
    segavg.cor <- colMeans(res.cor, na.rm=T)

    ## make a list of problem segments
    seg.flags <- rep(NA, nrow(res.pval))
    names(seg.flags) <- rownames(res.pval)
    flag.logical <- res.pval >= pcrit
    flag.logical[is.na(flag.logical)] <- FALSE
    for(i in 1:length(seg.flags))
        seg.flags[i] <- paste(names(flag.logical[i, flag.logical[i, ]]),
                              collapse = ", ")
    seg.flags <- seg.flags[seg.flags != ""]
    ## plot
    if(make.plot){
        p.val <- res.pval
        ## odd segs
        segs <- rwi
        this.seq <- seq(from=1, to=nrow(bins), by=2)
        odd.bins <- bins[this.seq, , drop=FALSE]
        odd.p.val <- p.val[, this.seq, drop=FALSE]
        nodd.bins <- nrow(odd.bins)
        com.segs <- as.data.frame(matrix(1, ncol=nseries, nrow=nyrs))
        rownames(com.segs) <- rnames
        colnames(com.segs) <- cnames
        flag.segs <- as.data.frame(matrix(NA, ncol=nseries, nrow=nyrs))
        rownames(flag.segs) <- rnames
        colnames(flag.segs) <- cnames
        ## loop through odd.bins
        tmp <- odd.p.val > pcrit
        for(i in 1:nseries){
            for(j in 1:nodd.bins){
                ## minus 1 deals with edge in segment graphing
                mask <- yrs%in%seq(from=odd.bins[j, 1], to=odd.bins[j, 2]-1)
                mask2 <- yrs%in%seq(from=odd.bins[j, 1], to=odd.bins[j, 2])
                ## note lack of complete overlap
                if(any(is.na(segs[mask, i]))) com.segs[mask, i] <- NA
                if(is.na(tmp[i, j])) com.segs[mask, i] <- NA
                if(!is.na(tmp[i, j]) && tmp[i, j]) flag.segs[mask2, i] <- 1
            }
            ## make sure any segment past year range is NA for
            ## com.segs this is a fix because there was an errant 1 in
            ## floating series that was not getting tagged because of
            ## the mask indexing: odd.bins[j, 2]-1 that deals with the
            ## bin overlap issue. This should fix it.
            series.yrs <- yr.range(segs[, i])
            mask3 <- yrs%in%seq(series.yrs[1], series.yrs[2])
            com.segs[!mask3, i] <- NA
        }
        idx.small <- yrs < min(odd.bins)
        idx.large <- yrs > max(odd.bins)
        com.segs[idx.small, ] <- NA
        com.segs[idx.large, ] <- NA
        flag.segs[idx.small, ] <- NA
        flag.segs[idx.large, ] <- NA

        extreme.year <- apply(segs, 2, yr.range)
        first.year <- extreme.year[1, ]
        rsult <- sort.int(first.year, decreasing=FALSE, index.return=TRUE)
        neworder <- rsult$ix
        segs <- segs[, neworder, drop=FALSE]
        com.segs <- com.segs[, neworder, drop=FALSE]
        flag.segs <- flag.segs[, neworder, drop=FALSE]




        segs.df <- data.frame(t(extreme.year[, neworder]))
        names(segs.df) <- c("first.yr", "last.yr")
        com.segs.df <- data.frame(t(apply(com.segs, 2, yr.range)))
        names(com.segs.df) <- c("first.yr", "last.yr")
        flag.segs.df <- data.frame(t(apply(flag.segs, 2, yr.range)))
        names(flag.segs.df) <- c("first.yr", "last.yr")

        op <- par(no.readonly=TRUE)
        col.pal <- c("#E41A1C", "#377EB8", "#4DAF4A")
        par(mar=c(4, 5, 4, 5) + 0.1, mgp=c(1.25, 0.25, 0), tcl=0.25)
        plot(yrs, segs[, 1], type="n", ylim=c(0, ncol(segs)),
             axes=FALSE, ylab="", xlab="Year",
             sub=paste("Segments: length=", seg.length, ",lag=",
             seg.lag, sep=""), ...)
        ## bounding poly for even series
        for(i in seq(from=1, to=nseries, by=2)){
            xx <- c(min.yr-100, max.yr+100)
            xx <- c(xx, rev(xx))
            yy <- c(i-0.5, i-0.5, i+0.5, i+0.5)
            polygon(xx, yy, col="grey90", border=NA)
        }
        abline(v=bins, col="grey", lty="dotted")
        axis(1, at=odd.bins)
        ## polygons for odd bins (go down from series line)
        for(i in seq(from=1, to=nseries)){
            y.deviation <- i-0.25
            ## whole segs
            xx <- c(segs.df[i, ], recursive=TRUE)
            xx <- c(xx, rev(xx))
            yy <- c(i, i, y.deviation, y.deviation)
            polygon(xx, yy, col=col.pal[3], border=NA)
            ## complete segs
            xx <- c(com.segs.df[i, ], recursive=TRUE)
            xx <- c(xx, rev(xx))
            polygon(xx, yy, col=col.pal[2], border=NA)
            ## flags
            xx <- c(flag.segs.df[i, ], recursive=TRUE)
            xx <- c(xx, rev(xx))
            polygon(xx, yy, col=col.pal[1], border=NA)
            ## guides
            guides.x <- c(odd.bins, recursive=T)
            guides.x <- sort(guides.x[!duplicated(guides.x)])
            guides.x <- guides.x[guides.x >= segs.df[i, 1]]
            guides.x <- guides.x[guides.x <= segs.df[i, 2]]
            segments(guides.x, i, guides.x, y.deviation, col="white")
        }

################################################################################
        ## even segs
        this.seq <- seq(from=2, to=nrow(bins), by=2)
        even.bins <- bins[this.seq, , drop=FALSE]
        even.p.val <- p.val[neworder, this.seq, drop=FALSE]
        neven.bins <- nrow(even.bins)
###    com.segs[,]=1
###    flag.segs[,]=NA
        com.segs <- as.data.frame(matrix(1, ncol=nseries, nrow=nyrs))
        rownames(com.segs) <- rnames
        colnames(com.segs) <- cnames
        flag.segs <- as.data.frame(matrix(NA, ncol=nseries, nrow=nyrs))
        rownames(flag.segs) <- rnames
        colnames(flag.segs) <- cnames

        ## loop through even.bins
        tmp <- even.p.val > pcrit
        for(i in 1:nseries){
            for(j in 1:neven.bins){
                mask <- yrs%in%seq(from=even.bins[j, 1], to=even.bins[j, 2]-1)
                mask2 <- yrs%in%seq(from=even.bins[j, 1], to=even.bins[j, 2])
                ## note lack of complete overlap
                if(any(is.na(segs[mask, i]))) com.segs[mask, i] <- NA
                if(is.na(tmp[i, j])) com.segs[mask, i] <- NA
                if(!is.na(tmp[i, j]) && tmp[i, j]) flag.segs[mask2, i] <- 1
            }
            ## make sure any segment past year range is NA for
            ## com.segs this is a fix because there was an errant 1 in
            ## floating series that was not getting tagged because of
            ## the mask indexing: odd.bins[j, 2]-1 that deals with the
            ## bin overlap issue. This should fix it.
            series.yrs <- yr.range(segs[, i])
            mask3 <- yrs%in%seq(series.yrs[1], series.yrs[2])
            com.segs[!mask3, i] <- NA
        }
        idx.small <- yrs < min(even.bins)
        idx.large <- yrs > max(even.bins)
        com.segs[idx.small, ] <- NA
        com.segs[idx.large, ] <- NA
        flag.segs[idx.small, ] <- NA
        flag.segs[idx.large, ] <- NA

        com.segs.df <- data.frame(t(apply(com.segs, 2, yr.range)))
        names(com.segs.df) <- c("first.yr", "last.yr")
        flag.segs.df <- data.frame(t(apply(flag.segs, 2, yr.range)))
        names(flag.segs.df) <- c("first.yr", "last.yr")

        ## guidelines
        axis(3, at=even.bins)
        ## polygons for even bins (go up from series line)
        for(i in seq(from=1, to=nseries)){
            y.deviation <- i+0.25
            ## whole segs
            xx <- c(segs.df[i, ], recursive=TRUE)
            xx <- c(xx, rev(xx))
            yy <- c(i, i, y.deviation, y.deviation)
            polygon(xx, yy, col=col.pal[3], border=NA)
            ## complete segs
            xx <- c(com.segs.df[i, ], recursive=TRUE)
            xx <- c(xx, rev(xx))
            polygon(xx, yy, col=col.pal[2], border=NA)
            ## flags
            xx <- c(flag.segs.df[i, ], recursive=TRUE)
            xx <- c(xx, rev(xx))
            polygon(xx, yy, col=col.pal[1], border=NA)
            ## guides
            guides.x <- c(even.bins, recursive=T)
            guides.x <- sort(guides.x[!duplicated(guides.x)])
            guides.x <- guides.x[guides.x >= segs.df[i, 1]]
            guides.x <- guides.x[guides.x <= segs.df[i, 2]]
            segments(guides.x, i, guides.x, y.deviation, col="white")

        }

        ## finish up plotting
        nsegs <- ncol(segs)
        odd.seq <- seq(from=1, to=nsegs, by=2)
        even.seq <- seq(from=2, to=nsegs, by=2)
        cnames.segs <- colnames(segs)
        axis(2, at=odd.seq,
             labels=cnames.segs[odd.seq], srt=45,
             tick=FALSE, las=2, cex.axis=label.cex)
        axis(4, at=even.seq,
             labels=cnames.segs[even.seq], srt=45,
             tick=FALSE, las=2, cex.axis=label.cex)
        abline(h=1:nseries, col="white")
        box()
        par(op)
    }
    res <- list(res.cor, res.pval, overall.cor, segavg.cor, seg.flags, bins)
    names(res) <- c("spearman.rho", "p.val", "overall", "avg.seg.rho",
         "flags", "bins")
    res
}
