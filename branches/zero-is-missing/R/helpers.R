### Creates a rwl data.frame with consecutive years
complete.rwl.df <- function(rwl, as.mat=FALSE) {
    cnames <- names(rwl)
    rwl.mat <- as.matrix(rwl)
    if (nrow(rwl.mat) == 0) {
        if (as.mat) {
            return(rwl.mat)
        } else {
            return(rwl)
        }
    }
    rnames <- rownames(rwl.mat)
    if (is.null(rnames)) {
        stop("automatic names found in 'rwl', must be explicit years")
    }
    yrs <- as.numeric(rnames)
    if (any(is.na(yrs))) {
        stop("non-numeric row name found in 'rwl'")
    }
    if (any(yrs != round(yrs))) {
        stop("row names of 'rwl' must represent integral-valued numbers")
    }
    min.yr <- min(yrs)
    max.yr <- max(yrs)
    yrs2 <- min.yr : max.yr
    rwl2 <- matrix(NA_real_,
                   nrow = max.yr - min.yr + 1,
                   ncol = length(rwl),
                   dimnames = list(as.character(yrs2), cnames))
    for (rname in rnames) {
        rwl2[rname, ] <- rwl.mat[rname, ]
    }
    if (as.mat) {
        rwl2
    } else {
        as.data.frame(rwl2)
    }
}

### Creates a series (vector) with consecutive years.
### Names must be present.
complete.series <- function(series) {
    names.series <- names(series)
    yrs <- as.numeric(names.series)
    min.yr <- min(yrs)
    max.yr <- max(yrs)
    series2 <- rep(NA_real_, max.yr - min.yr + 1)
    names(series2) <- as.character(min.yr : max.yr)
    series2[names.series] <- series
    series2
}

### Checks that all arguments are TRUE or FALSE
check.flags <- function(...) {
    flag.bad <- vapply(list(...),
                       function(x) { !(identical(x, TRUE) ||
                                       identical(x, FALSE)) },
                       TRUE,
                       USE.NAMES = FALSE)
    if (any(flag.bad)) {
        offending <- vapply(match.call(expand.dots=TRUE)[c(FALSE, flag.bad)],
                            deparse, "")
        stop(gettextf("must be TRUE or FALSE: %s",
                      paste(sQuote(offending), collapse=", "),
                      domain="R-dplR"),
             domain = NA)
    }
}

### Function to check if x is equivalent to its integer
### representation. Note: Returns FALSE for values that fall outside
### the range of the integer type. The result has the same shape as x;
### at least vector and array x are supported.
is.int <- function(x) {
    suppressWarnings(y <- x == as.integer(x))
    y[is.na(y)] <- FALSE
    y
}

### Converts from "year and suffix" presentation to dplR internal
### years, where year 0 (e.g. as a row name) is actually year 1 BC
dplr.year <- function(year, suffix) {
    switch(toupper(suffix),
           AD = ifelse(year > 0, year, as.numeric(NA)),
           BC = ifelse(year > 0, 1-year, as.numeric(NA)),
           BP = ifelse(year > 0, 1950-year, as.numeric(NA)),
           as.numeric(NA))
}

### Prints the contents of a matrix row together with column labels.
### By default, doesn't print NA values.  If show.all.na is TRUE,
### reports all-NA rows as one NA.
row.print <- function(x, drop.na=TRUE, show.all.na=TRUE, collapse=", ") {
    if (drop.na) {
        not.na <- !is.na(x)
        if (any(not.na)) {
            paste(colnames(x)[not.na], x[not.na], sep=": ", collapse=collapse)
        } else if (show.all.na) {
            as.character(NA)
        }
    } else {
        paste(colnames(x), x, sep=": ", collapse=collapse)
    }
}

### Returns indices of rows in matrix X that match with pattern.
row.match <- function(X, pattern) {
    which(apply(X, 1,
                function(x) {
                    all(is.na(x) == is.na(pattern)) &&
                    all(x == pattern, na.rm = TRUE)
                }))
}

### Increasing sequence.
### The equivalent of the C loop 'for(i=from;i<=to;i++){}'
### can be achieved by writing 'for(i in inc(from,to)){}'.
### Note that for(i in from:to) fails to do the same if to < from.
inc <- function(from, to) {
    if (is.numeric(to) && is.numeric(from) && to >= from) {
        seq(from=from, to=to)
    } else {
        integer(length=0)
    }
}

### Decreasing sequence. See inc.
dec <- function(from, to) {
    if (is.numeric(to) && is.numeric(from) && to <= from) {
        seq(from=from, to=to)
    } else {
        integer(length=0)
    }
}

### Levinson-Durbin algorithm for esmating an AR process from
### an autocovariance sequence.
### The same algorithm is used in stats::ar.yw().
levDurb <- function(acov) {
    order.max <- length(acov) - 1
    b <- vector(mode = "list", length = order.max)
    partialacf <- numeric(order.max)
    v <- rep(NA_real_, order.max)
    b[[1]] <- partialacf[1] <- b.temp <- acov[2] / acov[1]
    v[1] <- acov[1] * (1 - b.temp * b.temp)
    for (k in seq(from = 2, by = 1, length.out = order.max - 1)) {
        b.kk <- partialacf[k] <-
            (acov[k+1] - sum(b.temp * acov[k:2])) / v[k-1]
        b[[k]] <- b.temp <- c(b.temp - b.kk * rev(b.temp), b.kk)
        v[k] <- v[k-1] * (1 - b.kk * b.kk)
    }
    list(coeffs = b, var = v, partialacf = partialacf)
}

### AR function for chron, normalize1, normalize.xdate, ...
### The results should agree with those obtained by using standard ar.yw()
### in a non-NA case.
ar.func <- function(x, pass.na=TRUE, lag.max=NULL) {
    if (identical(pass.na, TRUE)) {
        na.act <- na.pass
    } else {
        na.act <- na.fail
    }
    idx.good <- which(!is.na(x))
    n.good <- length(idx.good)

    if (n.good > 0) {
        idx.seq <- idx.good[1]:idx.good[n.good]
        x.used <- x[idx.seq]
        n.used <- length(x.used)
        x.mean <- mean(x.used, na.rm = TRUE)
        x.demeaned <- x.used - x.mean
        ## Autocovariance sequence
        ## NOTE: if pass.na and NAs present, sequence may not be valid
        ACs <- acf(x.demeaned, type = "covariance",
                   na.action = na.act, lag.max = lag.max,
                   demean = FALSE, plot = FALSE)
        acov <- drop(ACs$acf)
        lags <- drop(ACs$lag)

        ## Check which lags are present in the autocovariance sequence
        notna <- !is.na(acov)
        acov <- acov[notna]
        lags <- lags[notna]
        lags.required <- seq(from = 0, by = 1,
                             length.out = max(1, length(lags)))
        lags.present <- lags.required %in% lags
        if (all(lags.present)) {
            acov <- acov[match(lags.required, lags)]
            order.max <- length(acov) - 1
        } else {
            first.missing <- lags.required[min(which(!lags.present))]
            order.max <- first.missing - 1
            acov <- acov[match(seq_len(first.missing) - 1, lags)]
        }

        ## Check if acov is OK, try to fix if it isn't
        if (order.max >= 1 && n.good < n.used) {
            acov.mat <- toeplitz(acov)
            eigval <- eigen(acov.mat, symmetric=TRUE)$values
            if (any(eigval < 0)) {
                acov.new <- try(nearPD(acov.mat, corr = FALSE,
                                       keepDiag = TRUE,
                                       do2eigen = TRUE,
                                       ensureSymmetry = FALSE),
                                silent = FALSE)
                if (inherits(acov.new, "try-error")) {
                    order.max <- 0
                } else {
                    n.rc <- order.max + 1
                    acov.mean <- numeric(n.rc)
                    col.idx <- rep(1:n.rc, each=n.rc)
                    row.idx <- rep(1:n.rc, n.rc)
                    new.mat <- acov.new$mat
                    for (k in 0:order.max) {
                        acov.mean[k+1] <-
                            mean(new.mat[abs(row.idx - col.idx) == k])
                    }
                    eigval2 <- eigen(toeplitz(acov.mean),
                                     symmetric=TRUE)$values
                    if (any(eigval2 < 0)) {
                        order.max <- 0
                    } else {
                        acov <- acov.mean
                    }
                }
            }
        }

        if (order.max >= 1) {
            ## AR process
            ARproc <- levDurb(acov)
            var.temp <- c(acov[1], ARproc$var)
            if (any(is.na(var.temp) | var.temp <= 0)) {
                order.max <- 0
            }
        }
        if (order.max >= 1) {
            aic <- n.good * log(var.temp) + 2 * seq_len(1 + order.max)
                partialacf <- ARproc$partialacf
            aic <- aic - min(aic)
            ## Select order by minimum AIC
            chosen.order <- which.min(aic) - 1
        } else if (order.max == 0) {
            var.temp <- acov[1]
            aic <- n.good * log(var.temp) + 2
            partialacf <- numeric(0)
            chosen.order <- order.max
        } else {
            aic <- numeric(0)
            partialacf <- numeric(0)
            chosen.order <- order.max
        }
        names(aic) <- as.character(seq(from = 0, by = 1,
                                       length.out = order.max + 1))

        if (chosen.order > 0) {
            ar <- ARproc$coeffs[[chosen.order]]
            y <- rep(NA_real_, length(x))
            resid.plus.mean <-
                x.used - c(rep(NA_real_, chosen.order),
                           drop(embed(x.demeaned[-n.used],
                                      chosen.order) %*% ar))
            ## Check that variance really decreased.
            ## Don't know if the contrary is possible in a case with NAs.
            if (n.good < n.used &&
                var(resid.plus.mean, na.rm=TRUE) >= var(x.used, na.rm=TRUE)) {
                chosen.order <- 0
                aic[-1] <- NA # unreliable values
            }
        }
        if (chosen.order > 0) {
            y[idx.seq] <- resid.plus.mean
            var.pred <- var.temp[chosen.order + 1] * n.good /
                (sum(!is.na(resid.plus.mean)) - 1)
        } else {
            ar <- numeric(0)
            var.pred <- var(x.used, na.rm=TRUE)
            y <- x
        }
        if (order.max < 0) {
            chosen.order <- NA_real_
            order.max <- NA_real_
        }
    } else {
        y <- x
        chosen.order <- NA_real_
        order.max <- NA_real_
        var.pred <- NA_real_
        aic <- ar <- partialacf <- numeric(0)
    }
    list(y = y, order = chosen.order, aic = aic, ar = ar,
         order.max = order.max, partialacf = partialacf,
         var.pred = var.pred)
}

### Range of years. Used in cms, rcs, rwl.stats, seg.plot, spag.plot, ...
yr.range <- function(x, yr.vec = as.numeric(names(x))) {
    range(yr.vec[!is.na(x)])
}

### Multiple ranges of years.
yr.ranges <- function(x, yr.vec = as.numeric(names(x))) {
    na.flag <- is.na(x)
    idx.good <- which(!na.flag)
    idx.bad <- which(na.flag)
    n <- length(x)
    res <- matrix(nrow=ceiling(n / 2), ncol=2)
    k <- 0
    while (length(idx.good) > 0) {
        first.good <- idx.good[1]
        idx.bad <- idx.bad[idx.bad > first.good]
        if (length(idx.bad) > 0) {
            first.bad <- idx.bad[1]
        } else {
            first.bad <- n + 1
        }
        idx.good <- idx.good[idx.good > first.bad]
        res[k <- k + 1, ] <- yr.vec[c(first.good, first.bad - 1)]
    }
    res[seq_len(k), , drop=FALSE]
}

### Used in cms, rcs, ...
sortByIndex <- function(x) {
    lowerBound <- which.min(is.na(x))
    c(x[lowerBound:length(x)], rep(NA, lowerBound - 1))
}

### Increment the given number (vector) x by one in the given base.
### Well, kind of: we count up to and including base (not base-1), and
### the smallest digit is one. Basically, we have a shift of one because
### of array indices starting from 1 instead of 0.  In case another
### digit is needed in the front, the result vector y grows.
count.base <- function(x, base) {
    n.x <- length(x)
    pos <- n.x
    y <- x
    y[pos] <- y[pos] + 1
    while (y[pos] == base + 1) {
        y[pos] <- 1
        if (pos == 1) {
            temp <- vector(mode="integer", length=n.x+1)
            temp[-1] <- y
            pos <- 2
            y <- temp
        }
        pos <- pos - 1
        y[pos] <- y[pos] + 1
    }
    y
}

### Compose a new name by attaching a suffix, which may partially
### replace the original name depending on the limit imposed on the
### length of names.
compose.name <- function(orig.name, alphabet, idx, limit) {
    idx.length <- length(idx)
    if (!is.null(limit) && idx.length > limit) {
        new.name <- ""
    } else {
        last.part <- paste(alphabet[idx], collapse="")
        if (is.null(limit)) {
            new.name <- paste0(orig.name, last.part)
        } else {
            new.name <- paste0(strtrim(orig.name, limit - idx.length),
                               last.part)
        }
    }
    new.name
}

### Fix names so that they are unique and no longer than the given
### length.  A reasonable effort will be done in the search for a set of
### unique names, although some stones will be left unturned. The
### approach should be good enough for all but the most pathological
### cases. The output vector keeps the names of the input vector.
fix.names <- function(x, limit=NULL, mapping.fname="", mapping.append=FALSE,
                      basic.charset=TRUE) {
    write.map <- FALSE
    n.x <- length(x)
    x.cut <- x
    rename.flag <- rep(FALSE, n.x)
    if (basic.charset) {
        bad.chars <- paste(c("[^",LETTERS,letters,0:9,"]"),collapse="")
        idx.bad <- grep(bad.chars, x.cut, perl=TRUE)
        if (length(idx.bad) > 0) {
            warning("characters outside a-z, A-Z, 0-9 present: renaming series")
            if (nchar(mapping.fname) > 0) {
                write.map <- TRUE
            }
            rename.flag[idx.bad] <- TRUE
            ## Remove inappropriate characters (replace with nothing)
            x.cut[idx.bad] <- gsub(bad.chars, "", x.cut[idx.bad])
        }
    }
    if (!is.null(limit)) {
        over.limit <- nchar(x.cut) > limit
        if (any(over.limit)) {
            warning("some names are too long: renaming series")
            if (nchar(mapping.fname) > 0) {
                write.map <- TRUE
            }
            rename.flag[over.limit] <- TRUE
            x.cut[over.limit] <- strtrim(x.cut[over.limit], limit)
        }
    }
    unique.cut <- unique(x.cut)
    n.unique <- length(unique.cut)
    ## Check if there are duplicate names after truncation and removal
    ## of inappropriate characters.  No duplicates => nothing to do
    ## beyond this point, except return the result.
    if (n.unique == n.x) {
        y <- x.cut
    } else {
        warning("duplicate names present: renaming series")
        if (nchar(mapping.fname) > 0) {
            write.map <- TRUE
        }

        y <- character(length=n.x)
        names(y) <- names(x)
        alphanumeric <- c(0:9, LETTERS, letters)
        n.an <- length(alphanumeric)
        ## First pass: Keep already unique names
        for (i in 1:n.unique) {
            idx.this <- which(x.cut %in% unique.cut[i])
            n.this <- length(idx.this)
            if (n.this == 1) {
                y[idx.this] <- x.cut[idx.this]
            }
        }

        if (!is.null(limit)) {
            x.cut <- strtrim(x.cut, limit-1)
        }
        x.cut[y != ""] <- NA
        unique.cut <- unique(x.cut) # may contain NA
        n.unique <- length(unique.cut)
        ## Second pass (exclude names that were set in the first pass):
        ## Make rest of the names unique
        for (i in 1:n.unique) {
            this.substr <- unique.cut[i]
            if (is.na(this.substr)) {# skip NA
                next
            }
            idx.this <- which(x.cut %in% this.substr)
            n.this <- length(idx.this)
            suffix.count <- 0
            for (j in 1:n.this){
                still.looking <- TRUE
                while (still.looking) {
                    suffix.count <- count.base(suffix.count, n.an)
                    proposed <-
                        compose.name(unique.cut[i],alphanumeric,suffix.count,limit)
                    if (nchar(proposed) == 0) {
                        warning("could not remap a name: some series will be missing")
                        still.looking <- FALSE
                        ## F for Fail...
                        proposed <- paste0(unique.cut[i], "F")
                    } else if (!any(y %in% proposed)) {
                        still.looking <- FALSE
                    }
                }
                this.idx <- idx.this[j]
                y[this.idx] <- proposed
                rename.flag[this.idx] <- TRUE
            }
        }
    }
    if (write.map) {
        if (mapping.append && file.exists(mapping.fname)) {
            map.file <- file(mapping.fname, "a")
        } else {
            map.file <- file(mapping.fname, "w")
        }
        for (i in which(rename.flag)) {
            if (x[i] != y[i]) {
                cat(x[i], "\t", y[i], "\n", file=map.file, sep = "")
            }
        }
        close(map.file)
    }
    y
}
