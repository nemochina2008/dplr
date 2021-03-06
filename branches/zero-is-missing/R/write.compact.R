write.compact <- function(rwl.df, fname, append=FALSE, prec=0.01,
                          mapping.fname="", mapping.append=FALSE, ...) {
    line.term <- "\x0D\x0A" # CR+LF, ASCII carriage return and line feed
    if (!is.data.frame(rwl.df)) {
        stop("'rwl.df' must be a data.frame")
    }
    if (!(prec == 0.01 || prec == 0.001)) {
        stop("'prec' must equal 0.01 or 0.001")
    }
    if (append && !file.exists(fname)) {
        stop(gettextf("file %s does not exist, cannot append", fname))
    }
    check.flags(mapping.append, append)

    ## Loop through series and write each one
    nseries <- ncol(rwl.df)
    yrs.all <- row.names(rwl.df)

    line.width <- 80 # max line width
    prec.rproc <- ifelse(prec == 0.01, 100, 1000) # reciprocal of precision
    max.field.width.width <-
        nchar(nchar(round(max(rwl.df, na.rm=TRUE) * prec.rproc)))
    max.n.width <- nchar(nrow(rwl.df)) # conservative
    max.i.width <- max(nchar(yrs.all)) # conservative
    ## Conservative length limit for the name of each series
    name.width <-
        line.width - max.field.width.width - max.n.width - max.i.width - 17

    col.names <- fix.names(x=names(rwl.df), limit=name.width,
                           mapping.fname=mapping.fname,
                           mapping.append=mapping.append, basic.charset=TRUE)

    ## Sort years using increasing order, reorder rwl.df accordingly
    yrs.all <- as.numeric(yrs.all)
    yrs.order <- sort.list(yrs.all)
    yrs.all <- yrs.all[yrs.order]
    rwl.df2 <- rwl.df[yrs.order, , drop=FALSE]

    if (append) {
        rwl.out <- file(fname, "a")
    } else {
        rwl.out <- file(fname, "w")
    }
    on.exit(close(rwl.out))
    missing.str <- 0

    for (l in seq_len(nseries)) {
        series <- rwl.df2[[l]]
        idx <- !is.na(series)
        yrs <- yrs.all[idx]
        series <- round(prec.rproc * series[idx])

        min.year <- min(yrs)
        max.year <- max(yrs)
        nyrs <- max.year - min.year + 1

        rwl.df.name <- col.names[l]

        ## Find missing data.
        missing.years <- setdiff(min.year:max.year, yrs)
        ## Mark missing data.
        if (length(missing.years) > 0) {
            yrs <- c(yrs, missing.years)
            series <- c(series, rep(missing.str, times=length(missing.years)))
            series.order <- sort.list(yrs)
            yrs <- yrs[series.order]
            series <- series[series.order]
        }
        ## Find negative values and mark as missing data
        series[series < 0] <- missing.str

        series <- as.character(series)
        field.width <- max(nchar(series))
        n.fields <- floor(line.width / field.width)
        n.lines <- floor(nyrs / n.fields)
        remainder <- nyrs - n.lines * n.fields

        ## Write header
        head1 <- paste0(nyrs, "=N", " ", min.year, "=I", " ")
        head2 <- paste0(ifelse(prec == 0.01, -2, -3), "(", n.fields, "F",
                        field.width, ".0)~")
        n.space <- line.width - nchar(head1) - nchar(head2) - nchar(rwl.df.name)
        if (n.space < 1) {
            ## since names are cut to length, this should not happen
            stop(gettextf("series %s: header line would be too long",
                          rwl.df.name))
        }
        cat(head1, rwl.df.name, rep(" ", n.space), head2, line.term,
            file=rwl.out, sep="")

        ## Write full lines
        for (i in seq_len(n.lines)) {
            end.idx <- i * n.fields
            line.str <- sprintf("%*s", field.width,
                                series[(end.idx - n.fields + 1) : end.idx])
            cat(line.str, line.term, file=rwl.out, sep="")
        }
        ## Write possibly remaining shorter line
        if (remainder > 0) {
            end.idx <- length(series)
            line.str <- sprintf("%*s", field.width,
                                series[(end.idx - remainder + 1) : end.idx])
            cat(line.str, line.term, file=rwl.out, sep="")
        }
    }
}
