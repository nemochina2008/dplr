`detrend.series` <-
    function(y, y.name = "", make.plot = TRUE,
             method = c("Spline", "ModNegExp", "Mean"),
             nyrs = NULL, f = 0.5, pos.slope = FALSE,
             zero.is.missing = TRUE)
{
    known.methods <- c("Spline", "ModNegExp", "Mean")
    method2 <- match.arg(arg = method,
                         choices = known.methods,
                         several.ok = TRUE)
    stopifnot(is.logical(zero.is.missing), length(zero.is.missing) == 1)
    ## Remove NA from the data (they will be reinserted later)
    good.y <- which(!is.na(y))
    n.good <- length(good.y)
    if (zero.is.missing && n.good > 0) {
        good.y <- good.y[y[good.y] != 0]
        n.good <- length(good.y)
    }
    if(n.good == 0) {
        stop("all values are 'NA'")
    } else if(!zero.is.missing && any(diff(good.y) != 1)) {
        stop("'NA's are not allowed in the middle of the series")
    }
    y2 <- y[seq(from = good.y[1], to = good.y[n.good], by = 1)]
    if (zero.is.missing) {
        good.y2 <- good.y - (good.y[1] - 1)
        good.flag <- rep(FALSE, length(y2))
        good.flag[good.y2] <- TRUE
        y2[!good.flag] <- NA_real_
    } else {
        ## Recode any zero values to 0.001
        good.y2 <- seq_along(y2)
        good.flag <- rep(TRUE, length(y2))
        y2[y2 == 0] <- 0.001
    }

    resids <- list()

    if ("ModNegExp" %in% method2) {
        ## Nec or lm
        nec.func <- function(Y, flag) {
            n <- length(Y)
            first.part <- rep(FALSE, n)
            first.part[seq_len(floor(n * 0.1))] <- TRUE
            a <- mean(Y[first.part & flag])
            b <- -0.01
            last.part <- rep(FALSE, n)
            last.part[floor(n * 0.9):n] <- TRUE
            k <- mean(Y[last.part & flag])
            nec <- nls(formula = Y[flag] ~ a * exp(b * which(flag)) + k,
                       start = list(a=a, b=b, k=k))
            if(coef(nec)[2] >= 0) stop()
            fits <- predict(nec)
            n.fits <- length(fits)
            if(fits[1] < fits[n.fits]) stop()
            if(fits[n.fits] < 0) stop()
            fits
        }
        ModNegExp <- try(nec.func(y2, good.flag), silent=TRUE)
        if (class(ModNegExp)=="try-error") {
            ## Straight line via linear regression
            lm1 <- lm(y2[good.flag] ~ which(good.flag))
            ModNegExp <- predict(lm1)
            if (coef(lm1)[2] > 0 && !pos.slope) {
                ModNegExp <- rep(mean(y2[good.flag]), n.good)
            }
        }
        resids$ModNegExp <- y2[good.flag] / ModNegExp
        do.mne <- TRUE
    } else {
        do.mne <- FALSE
    }

    if ("Spline" %in% method2) {
        ## Smoothing spline
        ## "n-year spline" as the spline whose frequency response is
        ## 50%, or 0.50, at a wavelength of 67%n years if nyrs and f
        ## are NULL
        if (is.null(nyrs)) {
            nyrs2 <- floor((good.y2[n.good] - good.y2[1] + 1) * 0.67)
        } else {
            nyrs2 <- nyrs
        }
        Spline <- ffcsaps(y=y2[good.flag], x=good.y2, nyrs=nyrs2, f=f)
        resids$Spline <- y2[good.flag] / Spline
        do.spline <- TRUE
    } else {
        do.spline <- FALSE
    }

    if ("Mean" %in% method2) {
        ## Fit a horiz line
        Mean <- rep(mean(y2[good.flag]), n.good)
        resids$Mean <- y2[good.flag] / Mean
        do.mean <- TRUE
    } else {
        do.mean <- FALSE
    }

    resids <- data.frame(resids)

    if(make.plot){
        op <- par(no.readonly=TRUE)
        on.exit(par(op))
        par(mar=c(2.5, 2.5, 2.5, 0.5) + 0.1, mgp=c(1.5, 0.5, 0))
        n.rows <- 1 + ncol(resids)
        mat <- matrix(seq_len(n.rows), n.rows, 1)
        layout(mat,
               widths=rep(0.5, ncol(mat)),
               heights=rep(1, nrow(mat)))

        plot(y2, type="l", ylab="mm",
             xlab=gettext("Age (Yrs)", domain="R-dplR"),
             main=gettextf("Raw Series %s", y.name, domain="R-dplR"))
        if (do.spline) {
            lines(good.y2, Spline, col="green", lwd=2)
        }
        if (do.mne) {
            lines(good.y2, ModNegExp, col="red", lwd=2)
        }
        if (do.mean) {
            lines(good.y2, Mean, col="blue", lwd=2)
        }

        if(do.spline){
            plot.resSpline <- rep(NA_real_, length(y2))
            plot.resSpline[good.flag] <- resids$Spline
            plot(plot.resSpline, type="l", col="green",
                 main=gettext("Spline", domain="R-dplR"),
                 xlab=gettext("Age (Yrs)", domain="R-dplR"),
                 ylab=gettext("RWI", domain="R-dplR"))
            abline(h=1)
        }

        if(do.mne){
            plot.resModnegexp <- rep(NA_real_, length(y2))
            plot.resModnegexp[good.flag] <- resids$ModNegExp
            plot(plot.resModnegexp, type="l", col="red",
                 main=gettext("Neg. Exp. Curve or Straight Line",
                 domain="R-dplR"),
                 xlab=gettext("Age (Yrs)", domain="R-dplR"),
                 ylab=gettext("RWI", domain="R-dplR"))
            abline(h=1)
        }

        if(do.mean){
            plot.resMean <- rep(NA_real_, length(y2))
            plot.resMean[good.flag] <- resids$Mean
            plot(plot.resMean, type="l", col="blue",
                 main=gettext("Horizontal Line (Mean)", domain="R-dplR"),
                 xlab=gettext("Age (Yrs)", domain="R-dplR"),
                 ylab=gettext("RWI", domain="R-dplR"))
            abline(h=1)
        }
    }

    resids2 <- matrix(NA, ncol=length(resids), nrow=length(y))
    resids2 <- data.frame(resids2)
    names(resids2) <- names(resids)
    if(!is.null(names(y))) row.names(resids2) <- names(y)
    resids2[good.y, ] <- resids

    ## Reorder columns of output to match the order of the argument
    ## "method".
    resids2 <- resids2[, method2]
    ## Make sure names (years) are included if there is only one method
    if(!is.data.frame(resids2)) names(resids2) <- names(y)

    resids2
}
