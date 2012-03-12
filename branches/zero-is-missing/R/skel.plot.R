`skel.plot` <-
    function(rw.vec, yr.vec = NULL, sname = "", filt.weight = 9,
             dat.out = FALSE, master=FALSE, plot=TRUE, metadata.out = FALSE)
{
    if (nchar(sname) > 7) {
        stop("'sname' must be a character string less than 8 characters long")
    }
    check.flags(metadata.out, plot, master, dat.out)

    na.mask <- is.na(rw.vec)
    notna.idx <- which(!na.mask)
    n.notna <- length(notna.idx)
    if (n.notna > 0) {
        rw.vec2 <- rw.vec[notna.idx[1]:notna.idx[n.notna]]
    } else {
        rw.vec2 <- numeric(0)
    }

    n.val <- length(rw.vec2)
    if (n.val > 840) {
        cat(gettextf("input series has length of %d\n", n.val))
        stop("long series (> 840) must be split into multiple plots")
    }
    if (n.val < filt.weight) {
        cat(gettextf("input series has length of %d", n.val),
            gettextf("'filt.weight' is %f\n", filt.weight), sep=", ")
        stop("'filt.weight' must not be larger than length of input series")
    }

    ## should wrap this into a function called skel.calc that returns the
    ## dates and skel

    ## if no yr then....
    if (is.null(yr.vec)) {
        yr.vec2 <- 0:(n.val - 1)
    } else {
        yr.vec2 <- yr.vec[notna.idx[1]:notna.idx[n.notna]]
    }
    original.na <- is.na(rw.vec2)       # NOTE! HUOM! metadata
    ## pad down to the nearest 10 if not already there
    min.yr <- min(yr.vec2)
    pad0 <- floor(min.yr / 10) * 10
    if (pad0 != min.yr) {
        pad.length <- min.yr - pad0
        rw.vec2 <- c(rep(NA_real_, pad.length), rw.vec2)
        yr.vec2 <- c(seq(from=pad0, by=1, length.out=pad.length), yr.vec2)
    } else {
        pad.length <- 0
    }

    ## detrend and pad
    rw.dt <- hanning(rw.vec2, filt.weight)
    skel <- rep(NA, length(rw.vec2))
    ## calc rel growth
    n.diff <- length(rw.vec2) - 1
    idx <- 2:n.diff
    temp.diff <- diff(rw.vec2)
    skel[idx] <- rowMeans(cbind(temp.diff[-n.diff],
                                -temp.diff[-1]), na.rm = FALSE) / rw.dt[idx]
    processing.na <- is.na(skel)        # NOTE! HUOM! metadata
    processing.na[original.na] <- FALSE
    skel[skel > 0] <- NA
    ## rescale from 0 to 10
    na.flag <- is.na(skel)
    if (all(na.flag)) {
        skel.range <- c(NA, NA)
    } else {
        skel.range <- range(skel[!na.flag])
    }
    newrange <- c(10, 1)
    mult.scalar <-
        (newrange[2] - newrange[1]) / (skel.range[2] - skel.range[1])
    skel <- newrange[1] + (skel - skel.range[1]) * mult.scalar
    skel[skel < 3] <- NA
    skel <- ceiling(skel)
    threshold.na <- is.na(skel)         # NOTE! HUOM! metadata
    threshold.na[original.na] <- FALSE
    threshold.na[processing.na] <- FALSE

    ## Variables for plotting
    ## page width
    pw <- 254
    ## page height
    ph <- 178
    ## row height
    rh <- 22
    ## row width
    rw <- 240
    ## spacer for text and dashed cutting lines
    spcr <- 5

    ## break series into sections of 120 years with an index
    yrs.col <- rw / 2 # n years per row
    n <- length(skel)
    n.rows <- ceiling(n / yrs.col)
    m <- seq_len(n.rows)
    row.index <- rep(m, each = yrs.col)[seq_len(n)]
    skel.df <- data.frame(yr=yr.vec2, skel)
    if (metadata.out) {
        skel.df <- cbind(skel.df, data.frame(original.NA, processing.NA,
                                             threshold.NA))
    }
    if (plot) {
        ## master page
        grid.newpage()
        vps <- list()
        y <- ph
        for (i in seq_len(min(n.rows, 7))) {
            y <- y - (rh + spcr)
            vps[[i]] <-
                viewport(x=unit(3, "mm"),
                         y=unit(y, "mm"),
                         width=unit(246, "mm"), height=unit(rh, "mm"),
                         just=c("left", "bottom"), name=LETTERS[i])
        }
        tree <-
            vpTree(viewport(width=unit(pw, "mm"), height=unit(ph, "mm"),
                            name="page"),
                   do.call(vpList, vps))

        ## set up page with the right number of rows
        pushViewport(tree)
        ## seq for 0 to plot width by 2mm
        tmp.1 <- seq(from=0, to=rw, by=2)
        tmp.2 <- seq(from=0, to=rh, by=2)
        tmp.3 <- seq(from=0, to=rw, by=20)
        for (i in m) {

            seekViewport(LETTERS[i])
            ## working code goes here - e.g., skelplot!
            grid.segments(x0=unit(tmp.1, "mm"), y0=unit(0, "mm"),
                          x1=unit(tmp.1, "mm"), y1=unit(rh, "mm"),
                          gp = gpar(col="green", lineend = "square", linejoin = "round"))
            grid.segments(x0=unit(0, "mm"), y0=unit(tmp.2, "mm"),
                          x1=unit(rw, "mm"), y1=unit(tmp.2, "mm"),
                          gp = gpar(col="green", lineend = "square", linejoin = "round"))

            ## decadal lines
            grid.segments(x0=unit(tmp.3, "mm"), y0=unit(0, "mm"),
                          x1=unit(tmp.3, "mm"), y1=unit(rh, "mm"),
                          gp = gpar(col = "black", lwd = 1.5, lty = "dashed",
                          lineend = "square", linejoin = "round"))

            ## lines on top and bottom of plot
            grid.lines(x=unit(c(0, rw), "mm"),
                       y=unit(c(rh, rh), "mm"),
                       gp=gpar(lwd = 2, lineend = "square", linejoin = "round"))
            grid.lines(x=unit(c(0, rw), "mm"),
                       y=unit(c(0, 0), "mm"),
                       gp=gpar(lwd = 2, lineend = "square", linejoin = "round"))
            ## plot x axis
            ## get this row's data
            skel.sub <- skel.df[row.index == i, ]
            end.yr <- length(skel.sub$yr)
            ticks <- seq(from=0, to=rw / 2, by=10)
            init.lab <- min(skel.sub$yr)
            x.labs <- seq(from=init.lab, length.out = length(ticks), by=10)
            if (master) {
                textY <- rh - 22.5
                textJust <- c("center", "top")
            } else {
                textY <- rh + 0.5
                textJust <- c("center", "bottom")
            }
            grid.text(label = x.labs,
                      x=unit(ticks * 2, "mm"),
                      y=unit(textY, "mm"),
                      just = textJust,
                      gp = gpar(fontsize=10))
            ## plot data
            notna.flag <- !is.na(skel.sub$skel)
            yr.m1 <- skel.sub$yr[notna.flag] - 1
            skel.temp <- skel.sub$skel[notna.flag]
            X <- (which(notna.flag) - 1) * 2
            n.X <- length(X)
            Y <- rep(0, n.X)
            X <- rep(X, each=2)
            Y <- c(Y, skel.temp * 2)
            Y <- Y[as.numeric(rbind(seq_len(n.X), seq_len(n.X)+n:X))]
            if (master) {
                Y <- 22 - Y
            }
            grid.lines(x=unit(X, "mm"),
                       y=unit(Y, "mm"),
                       gp = gpar(col = "black", lwd = 2,
                       lineend = "square", linejoin = "round"))
        }
        ## end arrow
        end.mm <- (end.yr - 1) * 2
        grid.lines(x=unit(c(end.mm, end.mm), "mm"),
                   y=unit(c(rh, 0), "mm"),
                   gp = gpar(lwd = 2, lineend = "square", linejoin = "round"))
        if (master) {
            Y <- c(rh, 16, 16)
        } else {
            Y <- c(0, 6, 6)
        }
        grid.polygon(x=unit(c(end.mm, end.mm, end.mm + 2), "mm"),
                     y=unit(Y, "mm"),
                     gp=gpar(fill = "black", lineend = "square", linejoin = "round"))
        ## start arrow and sample id
        start.mm <- pad.length * 2
        grid.lines(x=unit(c(start.mm, start.mm), "mm"),
                   y=unit(c(rh, 0), "mm"),
                   gp = gpar(lwd = 2, lineend = "square", linejoin = "round"))
        fontsize.sname <- ifelse(nchar(sname) > 6, 9, 10)
        if (master) {
            textY <- 1
            polyY <- c(rh, 16, 16)
            textJust <- c("left", "bottom")
        } else {
            textY <- rh - 1
            polyY <- c(0, 6, 6)
            textJust <- c("right", "bottom")
        }
        grid.polygon(x=unit(c(start.mm, start.mm, start.mm - 2), "mm"),
                     y=unit(polyY, "mm"),
                     gp=gpar(fill = "black", lineend = "square", linejoin = "round"))
        grid.text(label = sname,
                  x=unit(start.mm - 1, "mm"),
                  y=unit(textY, "mm"),
                  just = textJust,
                  rot = 90,
                  gp = gpar(fontsize=fontsize.sname))
    }
    if (dat.out) {
        return(skel.df)
    }
}
