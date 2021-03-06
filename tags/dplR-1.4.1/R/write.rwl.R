write.rwl <-
    function(rwl.df, fname, format=c("tucson", "compact", "tridas"), ...)
{
    format <- match.arg(format)
    switch(format,
           tucson = write.tucson(rwl.df, fname, ...),
           compact = write.compact(rwl.df, fname, ...),
           tridas = write.tridas(rwl.df, fname, ...))
}
