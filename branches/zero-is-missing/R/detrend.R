`detrend` <-
    function(rwl, y.name = names(rwl), make.plot = FALSE,
             method = c("Spline", "ModNegExp", "Mean"),
             ...)
{
    known.methods <- c("Spline", "ModNegExp", "Mean")
    method2 <- match.arg(arg = method,
                         choices = known.methods,
                         several.ok = TRUE)
    if (!is.data.frame(rwl)) {
        stop("'rwl' must be a data.frame")
    }
    rn <- row.names(rwl)

    args <- list(...)
    args[["method"]] <- method2
    if (!make.plot &&
        ("Spline" %in% method2 || "ModNegExp" %in% method2) &&
        !inherits(try(suppressWarnings(req.it <-
                                       require(iterators, quietly=TRUE)),
                      silent = TRUE),
                  "try-error") && req.it &&
        !inherits(try(suppressWarnings(req.fe <-
                                       require(foreach, quietly=TRUE)),
                      silent = TRUE),
                  "try-error") && req.fe) {
        it.rwl <- iter(rwl, by = "col")
        ## a way to get rid of "no visible binding" NOTE in R CMD check
        rwl.i <- NULL
        args[["make.plot"]] <- FALSE
        out <- foreach(rwl.i=it.rwl, .packages="dplR") %dopar% {
            args[["y"]] <- rwl.i
            fits <- do.call("detrend.series", args)
            if (is.data.frame(fits)) {
                row.names(fits) <- rn
            }
            fits
        }
    } else {
        out <- list()
        args[["make.plot"]] <- make.plot
        for (i in seq_len(ncol(rwl))) {
            args[["y"]] <- rwl[[i]]
            if (make.plot) {
                args[["y.name"]] <- y.name[i]
            }
            fits <- do.call("detrend.series", args)
            if (is.data.frame(fits)) {
                row.names(fits) <- rn
            }
            out[[i]] <- fits
        }
    }
    names(out) <- names(rwl)
    if (length(method2) == 1) {
        out <- data.frame(out, row.names = rn)
        names(out) <- y.name
    }
    out
}
