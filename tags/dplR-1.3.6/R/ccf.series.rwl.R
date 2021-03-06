ccf.series.rwl <- function(rwl,series, series.yrs=as.numeric(names(series)),
  seg.length=50,bin.floor=100,n=NULL, prewhiten = TRUE, biweight=TRUE,
  pcrit=0.05, lag.max=5, make.plot = TRUE,...){

  #run error checks
  qa.xdate(rwl,seg.length,n,bin.floor)
  if(lag.max > seg.length) stop("lag.max > seg.length. bad idea.")
  seg.lag=seg.length/2

  # Normalize.
  tmp = normalize.xdate(rwl,series,n,prewhiten,biweight)
  master = tmp$master
   # trim master so there are no NaN like dividing when
  # only one series for instance.
  idx.good = !is.nan(master)
  master = master[idx.good]
  yrs = as.numeric(names(master))

  series = tmp$series
  # trim series in case it was submitted stright from the rwl
  idx.good = !is.na(series)
  series.yrs = series.yrs[idx.good]
  series = series[idx.good]

  # clip series to master dimensions
  series = series[series.yrs %in% yrs]
  series.yrs = as.numeric(names(series))
  # clip master to series dimensions
  master = master[yrs %in% series.yrs]
  yrs = as.numeric(names(master))
  nyrs = length(series.yrs)

  if(is.null(bin.floor) || bin.floor == 0) min.bin = min(series.yrs)
  else min.bin = ceiling(min(series.yrs)/bin.floor)*bin.floor
  to = max(series.yrs)-seg.length-seg.length
  if(min.bin > to){
      cat("maximum year in (filtered) series:",max(series.yrs),"\n")
      cat("first bin begins: ",min.bin,"\n")
      cat("Cannot fit two segments (not enough years in the series).\n")
      stop("Shorten segment length or adjust the bin floor.")
  }
  bins = seq(from=min.bin, to=to+seg.length, by=seg.lag)
  bins = cbind(bins, bins+seg.length)
  nbins = nrow(bins)
  bin.names = paste(bins[,1],".", bins[,2],sep="")

  # structures for results
  lag.vec = seq(from=-lag.max, to=lag.max, by=1)
  res.cor = matrix(NA,length(lag.vec),nbins)
  rownames(res.cor)=paste('lag',lag.vec,sep='.')
  colnames(res.cor)=bin.names

  #loop through bins
  for(j in 1:nbins){
    mask = yrs%in%seq(from=bins[j,1], to=bins[j,2])
    # cor is NA if there is not complete overlap
    if(any(is.na(series[mask])) |
       any(is.na(master[mask])) |
       !any(mask) | table(mask)[2] < seg.length
       ){
        bin.ccf = NA
    }
    else {
      tmp = ccf(master[mask], series[mask], lag.max=lag.max,plot=FALSE)
      bin.ccf = c(tmp$acf)
    }
    res.cor[,j] = bin.ccf
  }
  # plot
  if(make.plot){
    ccf.df = data.frame(r=c(res.cor,recursive=T),
      bin=rep(colnames(res.cor),each=length(lag.vec)),
      lag=rep(lag.vec,nbins))
    # reorder bins so that lattice definitely keeps them
    # in ascending order (i.e., no factor order funnies with
    # long series)
    num.bins=bins[,1]
    ord.num=order(num.bins)
    char.bins=as.character(bins[,1])
    ord.char=order(char.bins)
    foo=data.frame(num.bins,ord.num,char.bins,ord.char)
    ccf.df$bin=factor(ccf.df$bin,levels(ccf.df$bin)[order(foo$ord.char)])

    sig = qnorm((1 + 1 - pcrit)/2)/sqrt(seg.length)
    sig = c(-sig,sig)
    ccf.plot = xyplot(r ~ lag | bin, data = ccf.df,
      ylim=range(ccf.df$r,sig,na.rm=T)*1.1,
      xlab='Lag', ylab='Correlation',col.line = NA,
      panel = function(x, y, col, ...) {
        panel.abline(h=seq(from=-1, to=1, by=0.1), lty='solid', col='gray')
        panel.abline(v=lag.vec,lty='solid',col='gray')
        panel.abline(h=0,v=0,lwd=2)
        panel.abline(h=sig,lwd=2,lty='dashed')
        col = ifelse(y > 0, '#E41A1C','#377EB8')
        # segments, dots for all r
        panel.segments(x1=x, y1=0, x2=x, y2=y, col = col, lwd= 2)
        panel.dotplot(x, y, col = col, cex = 1.25, ...)
      },...)
    trellis.par.set(strip.background = list(col = c('transparent')))
    print(ccf.plot)
  }
  res = list(res.cor,bins)
  names(res) = c('ccf', 'bins')
  res
}
