read.compact <- function(fname, zero.as.na = TRUE)
{
    check.flags(zero.as.na)
    res <- .Call(dplR.rcompact, path.expand(fname))
    min.year <- res[[1]]
    max.year <- res[[2]]
    series.ids <- res[[3]]
    series.min <- res[[4]]
    series.max <- res[[5]]
    series.mplier <- res[[6]]
    rw.mat <- res[[7]]
    project.comments <- res[[8]]
    rownames(rw.mat) <- min.year:max.year
    nseries <- ncol(rw.mat)
    ## Negative values are undoubtedly best treated as NA
    rw.mat[rw.mat < 0] <- NA
    if (zero.as.na) {
        rw.mat[rw.mat == 0] <- NA
    }

    cat(sprintf(ngettext(nseries,
                         "There is %d series\n",
                         "There are %d series\n",
                         domain="R-dplR"),
                nseries))
    cat(paste0(seq_len(nseries), "\t",
               series.ids, "\t",
               series.min, "\t",
               series.max, "\t",
               series.mplier, "\n"), sep="")
    if (length(project.comments) > 0) {
        cat(gettext("Comments:", domain="R-dplR"),
            paste(project.comments, collapse="\n"), sep="\n")
    }

    rw.df <- as.data.frame(rw.mat)
    names(rw.df) <- series.ids
    rw.df
}
