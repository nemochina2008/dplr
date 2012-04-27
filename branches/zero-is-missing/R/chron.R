`chron` <-
    function(x, prefix="xxx", biweight=TRUE, prewhiten=FALSE, ids=NULL)
{
    if(!is.data.frame(x)) {
        stop("'x' must be a data.frame")
    }
    check.flags(prewhiten, biweight)
    prefix.str <- as.character(prefix)
    if (length(prefix.str) != 1 || nchar(prefix.str) > 3) {
        stop("'prefix' must be a character string with less than 4 characters")
    }
    x2 <- complete.rwl.df(x, TRUE)
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
        colnames.x <- names(x)
        ## If all column names in 'x' are present in the set of row
        ## names in 'ids', arrange 'ids' to matching order
        rownames.ids <- row.names(ids)
        if (!is.null(rownames.ids) && all(colnames.x %in% rownames.ids)) {
            ids2 <- ids[colnames.x, c("tree", "core")]
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
            for (i in seq_len(n.unique)) {
                these.cols <- row.match(ids2, unique.ids[i, ])
                if (biweight) {
                    x.temp[, i] <-
                        apply(x2[, these.cols, drop=FALSE], 1, tbrm, C=9)
                } else {
                    x.temp[, i] <-
                        rowMeans(x2[, these.cols, drop=FALSE], na.rm=TRUE)
                }
            }
            ids2 <- unique.ids
            x2 <- x.temp
            message("Series with matching tree/core IDs have been averaged")
            message("Each unique tree/core ID pair adds 1 to sampling depth")
        }
    }
    tree.freq <- table(ids2$tree)
    samps <- rowSums(!is.na(x2))
    if (prewhiten) {
        x.ar <- apply(x2, 2, function(x) ar.func(x)$y)
        samps.res <- rowSums(!is.na(x.ar))
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
        if (prewhiten) {
            x.ar <- x.ar.temp
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
        ## TODO: Return order of AR model
        names(out) <- c(paste0(prefix.str, "std"), paste0(prefix.str, "res"),
                        "samp.depth", "samp.res")
    } else {
        out <- data.frame(std, samps)
        names(out) <- c(paste0(prefix.str, "std"), "samp.depth")
    }
    row.names(out) <- rownames(x2)
    out
}
