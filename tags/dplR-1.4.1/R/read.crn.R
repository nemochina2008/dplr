`read.crn` <- function(fname, header=NULL, encoding = getOption("encoding"))
{
    ## Open the data file for reading
    con <- file(fname, encoding = encoding)
    on.exit(close(con))
    if(is.null(header)){
        ## Try to determine if the file has a header. This is failable.
        ## Find out if an ITRDB header (3 lines) in file
        hdr1 <- readLines(con, n=1)
        if(length(hdr1) == 0)
            stop("File is empty")
        if(nchar(hdr1) < 10)
            stop("First line in the crn file ends before col 10")
        yrcheck <- suppressWarnings(as.numeric(substr(hdr1, 7, 10)))
        if(is.null(yrcheck) || length(yrcheck)!=1 || is.na(yrcheck) |
           yrcheck < -1e04 || yrcheck > 1e04) {
            cat("There appears to be a header in the crn file\n")
            is.head <- TRUE
        }
        else {
            cat("There does not appear to be a header in the crn file\n")
            is.head <- FALSE # No header lines
        }
    } else if(!is.logical(header)){
        stop("Header must be NULL or logical")
    } else{
        is.head <- header
    }
    if(is.head){
        ## Read 4th line - should be first data line
        dat1 <- readLines(con, n=4)
        if(length(dat1) < 4)
            stop("File has under 4 lines")
        dat1 <- dat1[4]
    } else{
        dat1 <- readLines(con, n=1)
        if(length(dat1) == 0)
            stop("File is empty")
    }
    if(nchar(dat1) < 10)
        stop("First data line ends before col 10")
    yrcheck <- as.numeric(substr(dat1, 7, 10))
    if(is.null(yrcheck) || length(yrcheck)!=1 || is.na(yrcheck) ||
       yrcheck < -1e04 || yrcheck > 1e04)
        stop("Cols 7-10 of first data line not a year")
    ## Look at last line to determine if Chronology Statistics are present
    ## if nchar <=63 then there is a stats line
    nlines <- length(readLines(con, n=-1))
    ## Read file
    skip.lines <- ifelse(is.head, 3, 0)
    ## Do nothing. read.fwf closes (and destroys ?!?) the file connection
    on.exit()
    ## Get chron stats if needed
    chron.stats <- read.fwf(con, c(6, 4, 6, 6, 6, 7, 9, 9, 10),
                            skip=nlines-1, strip.white=TRUE)
    ## Unintuitively, the connection object seems to have been destroyed
    ## by the previous read.fwf.  We need to create a new one.
    con <- file(fname, encoding = encoding)
    ## Really read file
    dat <- read.fwf(con, c(6, 4, rep(c(4, 3), 10)),
                    skip=skip.lines, strip.white=TRUE)
    ## If columns 3 in chron.stats is an integer then there is no
    ## statistics line
    if(is.numeric(chron.stats[, 3]) &&
       !is.int(as.numeric(chron.stats[, 3]))){
        names(chron.stats) <-
            c("SiteID", "nYears", "AC[1]", "StdDev", "MeanSens",
              "MeanRWI", "IndicesSum", "IndicesSS", "MaxSeries")
        cat("Embedded chronology statistics\n")
        print(chron.stats)
        ## Chop off last row of dat
        dat <- dat[-nrow(dat), , drop=FALSE]
    }

    series <- dat[, 1]
    series.ids <- unique(series)
    nseries <- length(series.ids)
    cat("There are ", nseries, " series\n", sep="")
    series.index <- match(series, series.ids)
    min.year <- (min(dat[, 2]) %/% 10) * 10
    max.year <- ((max(dat[, 2])+10) %/% 10) * 10
    span <- max.year - min.year + 1
    ncol.crn.mat <- nseries + 1
    crn.mat <- matrix(NA, ncol=ncol.crn.mat, nrow=span)
    colnames(crn.mat) <- c(as.character(series.ids), "samp.depth")
    rownames(crn.mat) <- min.year:max.year
    for(i in 1:nseries){
        decade.yr <- dat[series.index==i, 2]
        ## RWI
        x <- dat[series.index == i, -c(1, 2, seq(from=4, to=22, by=2)),
                 drop=FALSE]
        ## All sample depths
        y <- dat[series.index == i, -c(1, 2, seq(from=3, to=21, by=2)),
                 drop=FALSE]
        for(j in 1:nrow(x)) {
            yr <- decade.yr[j]
            if(j == 1) yr <- min.year
            for(k in 1:ncol(x)){
                if(is.na(x[j, k])) break
                crn.mat[as.character(yr), i] <- x[j, k]
                ## If i is one then make samp depth
                if(i == 1)
                    crn.mat[as.character(yr), ncol.crn.mat] <- y[j, k]
                yr <- yr + 1
            }
        }
    }
    ## Clean up NAs
    crn.mat[which(crn.mat[, -ncol.crn.mat] == 9990)] <- NA # column-major order
    crn.mat <-
        crn.mat[!apply(is.na(crn.mat[, -ncol.crn.mat, drop=FALSE]), 1, all),
                ,
                drop=FALSE]
    ## If samp depth is all 1 then dump it
    sd.one <- all(crn.mat[, ncol.crn.mat] == 1)
    if(sd.one) {
        save.names <- colnames(crn.mat)[-ncol.crn.mat]
        crn.mat <- crn.mat[, -ncol.crn.mat, drop=FALSE]
        crn.mat <- crn.mat/1000
        crn.df <- as.data.frame(crn.mat)
        colnames(crn.df) <- save.names
        cat("All embedded sample depths are one...Dumping from matrix\n")
    }
    else {
        crn.mat[, 1:nseries] <- crn.mat[, 1:nseries]/1000
        crn.df <- as.data.frame(crn.mat)
    }
    crn.df
}
