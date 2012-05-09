combine.rwl <- function(x, y = NULL) {
    combinator <- function(x, y) {
        dim.x2 <- ncol(x)
        dim.y2 <- ncol(y)
        if (dim.x2 > 0 && dim.y2 > 0) {
            dim.x1 <- nrow(x)
            dim.y1 <- nrow(y)
            dim2 <- dim.x2 + dim.y2
            if (dim.x1 > 0 && dim.y1 > 0) {
                years.x <- rownames(x)
                years.y <- rownames(y)
                min.x <- as.numeric(years.x[1])
                min.y <- as.numeric(years.y[1])
                max.x <- as.numeric(years.x[length(years.x)])
                max.y <- as.numeric(years.y[length(years.y)])
                min.year <- min(min.x, min.y)
                years <- min.year:max(max.x, max.y)
                new <- matrix(NA, nrow = length(years), ncol = dim2)
                new[(min.x-min.year+1):(max.x-min.year+1),seq_len(dim.x2)] <- x
                new[(min.y-min.year+1):(max.y-min.year+1),(dim.x2+1):dim2] <- y
                rownames(new) <- years
                colnames(new) <- c(colnames(x), colnames(y))
            } else if (dim.x1 > 0) {
                years <- rownames(x)
                new <- matrix(NA, nrow = dim.x1, ncol = dim2)
                new[, seq_len(dim.x2)] <- x
                rownames(new) <- years
                colnames(new) <- c(colnames(x), colnames(y))
            } else if (dim.y1 > 0) {
                years <- rownames(y)
                new <- matrix(NA, nrow = dim.y1, ncol = dim2)
                new[, seq_len(dim.x2)] <- y
                rownames(new) <- years
                colnames(new) <- c(colnames(x), colnames(y))
            } else {
                new <- matrix(NA, nrow = 0, ncol = dim2)
                colnames(new) <- c(colnames(x), colnames(y))
            }
            new
        } else if (dim.y2 == 0) {
            x
        } else {
            y
        }
    }
    ## check, if x is a list with data.frames inside it
    ## if yes: forget about y, and loop through all items of x and apply
    ## combinator() one by one
    ## if no: just use combinator() with x and y
    if (is.list(x)) {
        n <- length(x)
        if (n > 0 && all(vapply(x, is.data.frame, TRUE))) {
            new.mat <- complete.rwl.df(x[[1]], TRUE)
            for (i in inc(2, n)) {
                new.mat <- combinator(new.mat, complete.rwl.df(x[[i]], TRUE))
            }
        } else if (is.data.frame(x) && is.data.frame(y)) {
            new.mat <- combinator(complete.rwl.df(x, TRUE),
                                  complete.rwl.df(y, TRUE))
        } else {
            stop("Nothing to combine here. Please supply data.frames formatted according to the data standards in dplR.")
        }
    } else {
        stop("Nothing to combine here. Please supply data.frames formatted according to the data standards in dplR.")
    }
    as.data.frame(new.mat)
}
