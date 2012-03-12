`detrend.series` <-
    function(y, y.name = "", make.plot = TRUE,
             method = c("Spline", "ModNegExp", "Mean"),
             nyrs = NULL, f = 0.5, pos.slope = FALSE)
{
    known.methods <- c("Spline", "ModNegExp", "Mean")
    method2 <- match.arg(arg = method,
                         choices = known.methods,
                         several.ok = TRUE)
    check.flags(pos.slope, make.plot)
    ## Remove NA from the data (they will be reinserted later)
    good.y <- which(!is.na(y))
    n.good <- length(good.y)
    if(n.good == 0) {
        stop("all values are 'NA'")
    }
    y2 <- y[seq(from = good.y[1], to = good.y[n.good], by = 1)]
    n.y2 <- length(y2)
    good.y2 <- good.y - (good.y[1] - 1)
    good.flag <- rep(FALSE, n.y2)
    good.flag[good.y2] <- TRUE

    ## ## Recode any zero values to 0.001
    ## good.y2 <- seq_along(y2)
    ## good.flag <- rep(TRUE, n.y2)
    ## y2[y2 == 0] <- 0.001

    y3 <- y2[good.y2]

    resids <- list()

    if ("ModNegExp" %in% method2) {
        ## Nec or lm
        nec.func <- function(Y, not.na.flag, not.na.idx) {
            n <- length(Y)
            first.part <- rep(FALSE, n)
            first.part[seq_len(floor(n * 0.1))] <- TRUE
            a <- mean(Y[first.part & not.na.flag])
            b <- -0.01
            last.part <- rep(FALSE, n)
            last.part[floor(n * 0.9):n] <- TRUE
            k <- mean(Y[last.part & not.na.flag])
            nec <- nls(formula = Y[not.na.idx] ~ a * exp(b * not.na.idx) + k,
                       start = list(a=a, b=b, k=k))
            if (coef(nec)[2] >= 0) stop()
            fits <- predict(nec)
            last.fit <- fits[length(fits)]
            if (fits[1] < last.fit) stop()
            if (last.fit < 0) stop()
            fits
        }
        ModNegExp <- try(nec.func(y2, good.flag, good.y2), silent=TRUE)
        if (class(ModNegExp)=="try-error") {
            ## Straight line via linear regression
            lm1 <- lm(y3 ~ good.y2)
            ModNegExp <- predict(lm1)
            if (coef(lm1)[2] > 0 && !pos.slope) {
                ModNegExp <- rep(mean(y3), n.good)
            }
        }
        resids$ModNegExp <- y3 / ModNegExp
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
        Spline <- ffcsaps(y=y3, x=good.y2, nyrs=nyrs2, f=f)
        resids$Spline <- y3 / Spline
        do.spline <- TRUE
    } else {
        do.spline <- FALSE
    }

    if ("Mean" %in% method2) {
        ## Fit a horiz line
        Mean <- rep(mean(y3), n.good)
        resids$Mean <- y3 / Mean
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
            plot.resSpline <- rep(NA_real_, n.y2)
            plot.resSpline[good.flag] <- resids$Spline
            plot(plot.resSpline, type="l", col="green",
                 main=gettext("Spline", domain="R-dplR"),
                 xlab=gettext("Age (Yrs)", domain="R-dplR"),
                 ylab=gettext("RWI", domain="R-dplR"))
            abline(h=1)
        }

        if(do.mne){
            plot.resModnegexp <- rep(NA_real_, n.y2)
            plot.resModnegexp[good.flag] <- resids$ModNegExp
            plot(plot.resModnegexp, type="l", col="red",
                 main=gettext("Neg. Exp. Curve or Straight Line",
                 domain="R-dplR"),
                 xlab=gettext("Age (Yrs)", domain="R-dplR"),
                 ylab=gettext("RWI", domain="R-dplR"))
            abline(h=1)
        }

        if(do.mean){
            plot.resMean <- rep(NA_real_, n.y2)
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
