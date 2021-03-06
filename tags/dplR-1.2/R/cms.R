cms <- function(rwl,po) {
  # support funcs
  yr.range = function(x) {
    yr.vec = as.numeric(names(x))
    mask = !is.na(x)
    range(yr.vec[mask])
  }

  sortByIndex<-function(x){
    n<-length(x)
    lowerBound<-which.min(is.na(x))
    nonNACount<-sum(!is.na(x))
    c(x[lowerBound:n], rep(NA, lowerBound-1))
  }

  biologicalTrend<-function(theDat){
    tt<-theDat[,1]
    n<-NROW(theDat[,1])
    err1<-array(0,n)
    err2<-array(0,n)
    err3<-array(0,n)
    err4<-array(0,n)
    err6<-array(0,n)

    for (i in c(1:n)){
      err1[i]<-(theDat[i:i,2])^4
      err2[i]<--1*(2*((theDat[i:i,2])^2)*(2*theDat[i:i,1]+1))
      ans<-polyroot(c(err1[i],err2[i],1))
      err3[i]<-ans[1]
      err4[i]<-ans[2]
    }
    err5<-Re(err4)
    med<-median(err5)
    for (i in c(1:n)){
      err6[i]<-sqrt(med*(theDat[i:i,1]+1))-sqrt(med*theDat[i:i,1])
    }
    indicies<-cbind(tt,err6)
    indicies
  }
 #main func
  if(ncol(rwl) != nrow(po)) { stop('dimension problem: ncol(rw) != nrow(po)') }
  if(!all(po[,1] %in% colnames(rwl))) { stop('Series ids in po and rwl do not match') }
  series.yrs = apply(rwl, 2, yr.range)
  rownames(series.yrs) <- c('first','last')

  rwl.ord<-apply(rwl, 2, sortByIndex)
  rwca<-data.frame(matrix(NA, ncol=ncol(rwl.ord), nrow=sum(nrow(rwl.ord) + max(po[,2]))))
  colnames(rwca) <- colnames(rwl)
  for (i in 1:ncol(rwl.ord)){
    series = colnames(rwl.ord)[i]
    yrs2pith = po[po[,1] %in% series,2]
    rwca[(yrs2pith):(yrs2pith + nrow(rwl.ord)-1),i]<-rwl.ord[,i]
  }

  # divide each series by c curve and restore to cal years
  rwi = rwl
  yrs = as.numeric(rownames(rwi))
  for(i in 1:ncol(rwca)){
   index = cbind(which(!is.na(rwca[,i])),na.omit(rwca[,i]))
   index = biologicalTrend(index)
    y = na.omit(rwca[,i])/index[,2]
    first = series.yrs[1,i]
    last = series.yrs[2,i]
    rwi[yrs %in% first:last,i] = y
  }
  rwi
}