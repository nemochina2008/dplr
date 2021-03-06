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
row.print <- function(x, drop.na=TRUE, show.all.na=TRUE, collapse=", "){
    if(drop.na){
        not.na <- !is.na(x)
        if(any(not.na))
            paste(colnames(x)[not.na], x[not.na], sep=": ", collapse=collapse)
        else if(show.all.na)
            as.character(NA)
    } else{
        paste(colnames(x), x, sep=": ", collapse=collapse)
    }
}

### Returns indices of rows in matrix X that match with pattern.
row.match <- function(X, pattern){
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
inc <- function(from, to){
    if(is.numeric(to) && is.numeric(from) && to >= from)
        seq(from=from, to=to)
    else
        integer(length=0)
}

### Decreasing sequence. See inc.
dec <- function(from, to){
    if(is.numeric(to) && is.numeric(from) && to <= from)
        seq(from=from, to=to)
    else
        integer(length=0)
}

### AR function for chron, normalize1, normalize.xdate, ...
ar.func <- function(x){
    y <- x
    idx.goody <- !is.na(y)
    ar1 <- ar(y[idx.goody])
    y[idx.goody] <- ar1$resid+ar1$x.mean
    y
}

### Range of years. Used in cms, rcs, rwl.stats, seg.plot, spag.plot, ...
yr.range <- function(x, yr.vec = as.numeric(names(x))){
    range(yr.vec[!is.na(x)])
}

### Used in cms, rcs, ...
sortByIndex <- function(x){
    lowerBound <- which.min(is.na(x))
    c(x[lowerBound:length(x)], rep(NA, lowerBound-1))
}

### Increment the given number (vector) x by one in the given base.
### Well, kind of: we count up to and including base (not base-1), and
### the smallest digit is one. Basically, we have a shift of one because
### of array indices starting from 1 instead of 0.  In case another
### digit is needed in the front, the result vector y grows.
count.base <- function(x, base){
    n.x <- length(x)
    pos <- n.x
    y <- x
    y[pos] <- y[pos] + 1
    while (y[pos] == base + 1){
        y[pos] <- 1
        if (pos == 1){
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
compose.name <- function(orig.name, alphabet, idx, limit){
    idx.length <- length(idx)
    if (!is.null(limit) && idx.length > limit){
        new.name <- ""
    } else{
        last.part <- paste(alphabet[idx], collapse="")
        if(is.null(limit))
            new.name <- paste(orig.name, last.part, sep="")
        else
            new.name <- paste(strtrim(orig.name,limit-idx.length),
                              last.part, sep="")
    }
    new.name
}

### Fix names so that they are unique and no longer than the given
### length.  A reasonable effort will be done in the search for a set of
### unique names, although some stones will be left unturned. The
### approach should be good enough for all but the most pathological
### cases. The output vector keeps the names of the input vector.
fix.names <- function(x, limit=NULL, mapping.fname="", mapping.append=FALSE,
                      basic.charset=TRUE){
    write.map <- FALSE
    n.x <- length(x)
    x.cut <- x
    rename.flag <- rep(FALSE, n.x)
    if(basic.charset){
        bad.chars <- paste(c("[^",LETTERS,letters,0:9,"]"),collapse="")
        idx.bad <- grep(bad.chars, x.cut, perl=TRUE)
        if(length(idx.bad) > 0){
            warning("characters outside a-z, A-Z, 0-9 present: renaming series")
            if(nchar(mapping.fname)>0)
                write.map <- TRUE
            rename.flag[idx.bad] <- TRUE
            ## Remove inappropriate characters (replace with nothing)
            x.cut[idx.bad] <- gsub(bad.chars, "", x.cut[idx.bad])
        }
    }
    if(!is.null(limit)){
        over.limit <- nchar(x.cut) > limit
        if(any(over.limit)){
            warning("some names are too long: renaming series")
            if(nchar(mapping.fname)>0)
                write.map <- TRUE
            rename.flag[over.limit] <- TRUE
            x.cut[over.limit] <- strtrim(x.cut[over.limit], limit)
        }
    }
    unique.cut <- unique(x.cut)
    n.unique <- length(unique.cut)
    ## Check if there are duplicate names after truncation and removal
    ## of inappropriate characters.  No duplicates => nothing to do
    ## beyond this point, except return the result.
    if (n.unique == n.x){
        y <- x.cut
    } else {
        warning("duplicate names present: renaming series")
        if(nchar(mapping.fname)>0)
            write.map <- TRUE

        y <- character(length=n.x)
        names(y) <- names(x)
        alphanumeric <- c(0:9, LETTERS, letters)
        n.an <- length(alphanumeric)
        ## First pass: Keep already unique names
        for (i in 1:n.unique){
            idx.this <- which(x.cut %in% unique.cut[i])
            n.this <- length(idx.this)
            if (n.this == 1)
                y[idx.this] <- x.cut[idx.this]
        }

        if(!is.null(limit))
            x.cut <- strtrim(x.cut, limit-1)
        x.cut[y != ""] <- NA
        unique.cut <- unique(x.cut) # may contain NA
        n.unique <- length(unique.cut)
        ## Second pass (exclude names that were set in the first pass):
        ## Make rest of the names unique
        for (i in 1:n.unique){
            this.substr <- unique.cut[i]
            if (is.na(this.substr)) # skip NA
                next
            idx.this <- which(x.cut %in% this.substr)
            n.this <- length(idx.this)
            suffix.count <- 0
            for (j in 1:n.this){
                still.looking <- TRUE
                while (still.looking){
                    suffix.count <- count.base(suffix.count, n.an)
                    proposed <-
                        compose.name(unique.cut[i],alphanumeric,suffix.count,limit)
                    if (nchar(proposed)==0){
                        warning("could not remap a name: some series will be missing")
                        still.looking <- FALSE
                        ## F for Fail...
                        proposed <- paste(unique.cut[i], "F", sep="")
                    } else if (!any(y %in% proposed)){
                        still.looking <- FALSE
                    }
                }
                this.idx <- idx.this[j]
                y[this.idx] <- proposed
                rename.flag[this.idx] <- TRUE
            }
        }
    }
    if(write.map){
        if(mapping.append && file.exists(mapping.fname))
            map.file <- file(mapping.fname, "a")
        else
            map.file <- file(mapping.fname, "w")
        for(i in which(rename.flag)){
            if(x[i] != y[i])
                cat(x[i], "\t", y[i], "\n", file=map.file, sep = "")
        }
        close(map.file)
    }
    y
}
