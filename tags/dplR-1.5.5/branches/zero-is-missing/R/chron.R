`chron` <-
    function(x, prefix="xxx", biweight=TRUE, prewhiten=FALSE, ids=NULL,
             x.out = FALSE)
{
    if(!is.data.frame(x)) {
        stop("'x' must be a data.frame")
    }
    check.flags(prewhiten, biweight, x.out)
    prefix.str <- as.character(prefix)
    if (length(prefix.str) != 1 || nchar(prefix.str) > 3) {
        stop("'prefix' must be a character string with less than 4 characters")
    }
    x2 <- complete.rwl.df(x, TRUE)
    rnames <- rownames(x2)
    cnames <- colnames(x2)
    n.cores <- length(x)
    if (is.null(ids)) {
        ids2 <- data.frame(tree=seq_along(x), core=rep(1, n.cores))
    } else {
        if (!is.data.frame(ids) || !all(c("tree", "core") %in% names(ids))) {
            stop("'ids' must be a data.frame with columns 'tree' and 'core'")
        }
        if (!all(vapply(ids, is.numeric, TRUE))) {
            stop("'ids' must have numeric columns")
        }
        ## If all column names in 'x' are present in the set of row
        ## names in 'ids', arrange 'ids' to matching order
        rownames.ids <- row.names(ids)
        if (!is.null(rownames.ids) && all(cnames %in% rownames.ids)) {
            ids2 <- ids[cnames, c("tree", "core")]
        } else if (nrow(ids) == n.cores) {
            ids2 <- ids[c("tree", "core")]
        } else {
            stop("dimension problem: ", "'ncol(x)' != 'nrow(ids)'")
        }
        row.names(ids2) <- NULL
        unique.ids <- unique(ids2)
        n.unique <- nrow(unique.ids)
        if (n.unique < n.cores) {
            ## If more than one columns of 'x' share a tree/core ID pair,
            ## the columns are averaged (robustly) and treated as one core
            x.temp <- matrix(NA_real_, nrow(x2), n.unique)
            names.temp <- character(n.unique)
            for (i in seq_len(n.unique)) {
                these.cols <- row.match(ids2, unique.ids[i, ])
                if (biweight) {
                    x.temp[, i] <-
                        apply(x2[, these.cols, drop=FALSE], 1, tbrm, C=9)
                } else {
                    x.temp[, i] <-
                        rowMeans(x2[, these.cols, drop=FALSE], na.rm=TRUE)
                }
                names.temp[i] <-
                    paste0(cnames[these.cols], collapse=";")
            }
            ids2 <- unique.ids
            cnames <- names.temp
            x2 <- x.temp
            rownames(x2) <- rnames
            colnames(x2) <- cnames
            message("Series with matching tree/core IDs have been averaged")
            message("Each unique tree/core ID pair adds 1 to sampling depth")
        }
    }
    tree.freq <- table(ids2$tree)
    samps <- rowSums(!is.na(x2))
    orig.na <- apply(is.na(x2), 2, sum)
    n.orig <- nrow(x2) - orig.na
    if (prewhiten) {
        ar.tmp <- apply(x2, 2, ar.func)
        x.ar <- vapply(ar.tmp, function(x) x$y, numeric(nrow(x2)))
        rownames(x.ar) <- rnames
        colnames(x.ar) <- cnames
        ar.order <- vapply(ar.tmp, function(x) x$order, 0)
        isna.xar <- is.na(x.ar)
        ar.na <- apply(isna.xar, 2, sum)
        n.ar <- nrow(x2) - ar.na
        samps.res <- rowSums(!isna.xar)
        names(ar.order) <- cnames
    }

    ## Series prior to within-tree averaging,
    ## except that series with shared IDs have been averaged
    x2.pre <- x2
    if (prewhiten) {
        x.ar.pre <- x.ar
    }

    if (any(tree.freq != tree.freq[1]) || (biweight && any(tree.freq != 1))) {
        unique.trees <- as.numeric(names(tree.freq))
        n.unique <- length(unique.trees)
        x.temp <- matrix(NA_real_, nrow(x2), n.unique)
        if (prewhiten) {
            x.ar.temp <- matrix(NA_real_, nrow(x.ar), n.unique)
        }
        one.flag <- tree.freq == 1
        multi.flag <- !one.flag
        multi.idx <- which(multi.flag)
        foo <- ids2$tree %in% unique.trees[one.flag]
        x.temp[, one.flag] <- x2[, foo, drop=FALSE]
        if (prewhiten) {
            x.ar.temp[, one.flag] <- x.ar[, foo, drop=FALSE]
        }
        if (biweight) {
            for (i in multi.idx) {
                foo <- ids2$tree == unique.trees[i]
                x.temp[, i] <- apply(x2[, foo, drop=FALSE], 1, tbrm)
                if (prewhiten) {
                    x.ar.temp[, i] <- apply(x.ar[, foo, drop=FALSE], 1, tbrm)
                }
            }
        } else {
            for (i in multi.idx) {
                foo <- ids2$tree == unique.trees[i]
                x.temp[, i] <- rowMeans(x2[, foo, drop=FALSE], na.rm=TRUE)
                if (prewhiten) {
                    x.ar.temp[, i] <-
                        rowMeans(x.ar[, foo, drop=FALSE], na.rm=TRUE)
                }
            }
        }
        x2 <- x.temp
        rownames(x2) <- rnames
        if (prewhiten) {
            x.ar <- x.ar.temp
            rownames(x.ar) <- rnames
        }
    }
    ## (Robust) mean of means
    if (biweight) {
        std <- apply(x2, 1, tbrm, C=9)
        if (prewhiten) {
            res <- apply(x.ar, 1, tbrm, C=9)
            res[is.nan(res)] <- NA_real_
        }
    } else {
        std <- rowMeans(x2, na.rm=TRUE)
        if (prewhiten) {
            res <- rowMeans(x.ar, na.rm=TRUE)
            res[is.nan(res)] <- NA_real_
        }
    }
    if (prewhiten) {
        out <- data.frame(std, res, samps, samps.res)
        names(out) <- c(paste0(prefix.str, "std"), paste0(prefix.str, "res"),
                        "samp.depth", "samp.res")
        attr(out, "stats") <-
            data.frame(order=ar.order, n.std=n.orig, n.res=n.ar)
        row.names(out) <- rnames
        if (x.out) {
            list(chron=out, x.std=x2.pre, x.res=x.ar.pre)
        } else {
            out
        }
    } else {
        out <- data.frame(std, samps)
        names(out) <- c(paste0(prefix.str, "std"), "samp.depth")
        attr(out, "stats") <- data.frame(n.std=n.orig)
        row.names(out) <- rnames
        if (x.out) {
            list(chron=out, x.std=x2.pre)
        } else {
            out
        }
    }
}
