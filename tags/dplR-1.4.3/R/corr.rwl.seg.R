corr.rwl.seg <- function(rwl, seg.length=50, bin.floor=100, n=NULL,
                         prewhiten = TRUE, pcrit=0.05, biweight=TRUE,
                         make.plot = TRUE, label.cex=1, ...){

    ## helper function
    yr.range <- function(x, yr.vec=as.numeric(names(x))) {
        if(any(mask <- !is.na(x))) range(yr.vec[mask])
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
    seg.lag <- seg.length / 2
    nseries <- ncol(rwl)
    if(nseries < 2) stop("At least 2 series are needed in 'rwl'")
    yrs <- as.numeric(rnames)
    nyrs <- length(yrs)
    min.yr <- min(yrs)
    max.yr <- max(yrs)
    if(is.null(bin.floor) || bin.floor == 0) min.bin <- min.yr
    else min.bin <- ceiling(min.yr/bin.floor) * bin.floor
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
        tmp <- cor.test(series, master,
                        method = "spearman", alternative = "g")
        overall.cor[i, 1] <- tmp$estimate
        overall.cor[i, 2] <- tmp$p.val
    }
    ## avg seg correlation
    segavg.cor <- colMeans(res.cor, na.rm=T)

    ## make a list of problem segments
    seg.flags <- rep(NA, nseries)
    names(seg.flags) <- cnames
    flag.logical <- res.pval >= pcrit
    flag.logical[is.na(flag.logical)] <- FALSE
    for(i in 1:length(seg.flags))
        seg.flags[i] <- paste(names(flag.logical[i, flag.logical[i, ]]),
                              collapse = ", ")
    seg.flags <- seg.flags[seg.flags != ""]

    ## plot
    if(make.plot){
        p.val <- res.pval
        segs <- rwi
        extreme.year <- apply(segs, 2, yr.range, yr.vec=yrs)
        first.year <- extreme.year[1, ]
        rsult <- sort.int(first.year, decreasing=FALSE, index.return=TRUE)
        neworder <- rsult$ix
        segs <- segs[, neworder, drop=FALSE]
        segs.mat <- t(extreme.year[, neworder])
        first.year <- first.year[neworder]
        last.year <- extreme.year[2, neworder]

        nsegs <- ncol(segs)
        op <- par(no.readonly=TRUE)
        on.exit(par(op), add=TRUE)
        col.pal <- c("#E41A1C", "#377EB8", "#4DAF4A")
        par(mar=c(4, 5, 4, 5) + 0.1, mgp=c(1.25, 0.25, 0), tcl=0.25)
        plot(yrs, segs[, 1], type="n", ylim=c(0, nsegs),
             axes=FALSE, ylab="", xlab=gettext("Year"),
             sub=gettextf("Segments: length=%d,lag=%d", seg.length, seg.lag,
             domain="R-dplR"),
             ...)
        ## bounding poly for even series
        xx <- c(min.yr-100, max.yr+100)
        xx <- c(xx, rev(xx))
        for(i in seq(from=1, to=nseries, by=2)){
            yy <- c(i-0.5, i-0.5, i+0.5, i+0.5)
            polygon(xx, yy, col="grey90", border=NA)
        }
        abline(v=bins, col="grey", lty="dotted")

        ## First odd segs, then even segs
        y.offset <- c(-0.25, 0.25)
        ax <- c(1, 3)
        for(odd.even in 1:2){
            this.seq <- seq(from=odd.even, to=nbins, by=2)
            these.bins <- bins[this.seq, , drop=FALSE]
            com.segs <- matrix(1, ncol=nseries, nrow=nyrs)
            flag.segs <- matrix(NA, ncol=nseries, nrow=nyrs)
            ## loop through these.bins
            tmp <- p.val[neworder, this.seq, drop=FALSE] > pcrit
            for(i in 1:nseries){
                for(j in 1:nrow(these.bins)){
                    ## minus 1 deals with edge in segment graphing
                    mask <- yrs%in%seq(from = these.bins[j, 1],
                                       to = these.bins[j, 2]-1)
                    ## note lack of complete overlap
                    if(any(is.na(segs[mask, i])))
                        com.segs[mask, i] <- NA
                    if(is.na(tmp[i, j]))
                        com.segs[mask, i] <- NA
                    else if(tmp[i, j]){
                        mask2 <- yrs%in%seq(from = these.bins[j, 1],
                                            to = these.bins[j, 2])
                        flag.segs[mask2, i] <- 1
                    }
                }
                ## make sure any segment past year range is NA for
                ## com.segs this is a fix because there was an errant
                ## 1 in floating series that was not getting tagged
                ## because of the mask indexing: these.bins[j, 2]-1
                ## that deals with the bin overlap issue. This should
                ## fix it.
                mask3 <- yrs%in%seq(first.year[i], last.year[i])
                com.segs[!mask3, i] <- NA
            }
            idx.small.large <-
                yrs < min(these.bins) | yrs > max(these.bins)
            com.segs[idx.small.large, ] <- NA
            flag.segs[idx.small.large, ] <- NA

            com.segs.mat <-
                t(apply(com.segs, 2, yr.range, yr.vec=yrs))
            flag.segs.mat <-
                t(apply(flag.segs, 2, yr.range, yr.vec=yrs))

            axis(ax[odd.even], at=these.bins)
            ## polygons for these bins (go down or up from series line)
            y.deviation <- y.offset[odd.even]
            guides.x.base <- c(these.bins, recursive=T)
            guides.x.base <- sort(guides.x.base[!duplicated(guides.x.base)])
            for(i in seq(from=1, to=nseries)){
                y.deviation <- y.deviation + 1
                ## whole segs
                xx <- segs.mat[i, ]
                xx <- c(xx, rev(xx))
                yy <- c(i, i, y.deviation, y.deviation)
                polygon(xx, yy, col=col.pal[3], border=NA)
                ## complete segs
                xx <- com.segs.mat[i, ]
                xx <- c(xx, rev(xx))
                polygon(xx, yy, col=col.pal[2], border=NA)
                ## flags
                xx <- flag.segs.mat[i, ]
                xx <- c(xx, rev(xx))
                polygon(xx, yy, col=col.pal[1], border=NA)
                ## guides
                guides.x <- guides.x.base[guides.x.base >= segs.mat[i, 1]]
                guides.x <- guides.x[guides.x <= segs.mat[i, 2]]
                segments(guides.x, i, guides.x, y.deviation, col="white")
            }
        }

        ## finish up plotting
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
    }

    list(spearman.rho = res.cor, p.val = res.pval, overall = overall.cor,
         avg.seg.rho = segavg.cor, flags = seg.flags, bins = bins)
}
