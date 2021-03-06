`read.rwl` <-
function(fname, header=NULL)
{
  # Open the data file for reading
  dat=file(fname,"r")
  if(is.null(header)){
    # Try to determine if the file has a header. This is failable.
    # Find out if an ITRDB header  (3 lines) in file
    hdr1=readLines(dat,n=1)
    if(nchar(hdr1) < 12) stop("First line in .rwl file ends before col 12")
    yrcheck=suppressWarnings(as.numeric(substr(hdr1,7,10)))
    if(is.null(yrcheck) | length(yrcheck)!=1 | is.na(yrcheck) |
       yrcheck < -1e04 | yrcheck > 1e04) {
      cat("There appears to be a header in the rwl file\n")
      ihead=TRUE
    }
    else {
     	ihead=FALSE # No header lines
      cat("There does not appear to be a header in the rwl file\n")
    }
    close(dat)
    dat=file(fname,"r")
  }
  else ihead = header
  if(ihead){
    # Read 4th line - should be first data line
    dat1=readLines(dat,n=4)[-c(1:3)]
  }
  else dat1=readLines(dat,n=1)
  yrcheck=as.numeric(substr(dat1,9,12))
  if(is.null(yrcheck) | length(yrcheck)!=1) {
    stop("Cols 9-12 of first data line not a year")
  }
  close(dat)

  skip.lines=ifelse(ihead,3,0)
  dat=read.fwf(fname,c(8,4,rep(6,10)),skip=skip.lines,strip.white=TRUE,
                  blank.lines.skip=TRUE)
  # Remove any blank lines at the end of the file, for instance
  dat=dat[!apply(is.na(dat),1,all),]

  series=dat[,1]
  series.ids=unique(series)
  nseries=length(series.ids)

  cat("There are ",nseries," series\n",sep="")

  series.index=match(series,series.ids)
  min.year=(min(dat[,2]) %/% 10) * 10
  max.year=((max(dat[,2])+10) %/% 10) * 10
  span=max.year - min.year + 1
  rw.mat=matrix(NA,ncol=nseries,nrow=span)
  colnames(rw.mat)=as.character(series.ids)
  rownames(rw.mat)=min.year:max.year

  for(i in 1:nseries){
    id=unique(dat[series.index==i,1])
    decade.yr=dat[series.index==i,2]
    first.yr=dat[series.index==i,2][1]
    x=dat[series.index==i,-c(1,2)]
    for(j in 1:nrow(x)) {
      yr=decade.yr[j]
      for(k in 1:ncol(x)){
        if(is.na(x[j,k])) break
        rw.mat[as.character(yr),i]=x[j,k]
        yr=yr + 1
      }
    }
  }
  # Clean up NAs as either 999 or -9999
  if(any(rw.mat==-999,na.rm=TRUE) | any(rw.mat==999,na.rm=TRUE)) {
    prec=0.01
    rw.mat[rw.mat==999]=NA
    rw.mat[rw.mat==-999]=NA
  }
  if(any(rw.mat==-9999,na.rm=TRUE) | any(rw.mat==9999,na.rm=TRUE)){
    prec=0.001
    rw.mat[rw.mat==-9999]=NA
    rw.mat[rw.mat==9999]=NA
  }
  # Convert to mm
  rw.mat=rw.mat * prec
  # Fix internal NAs. These are coded as 0 in the DPL programs
  fix.internal.na=function(x){
    for(i in 2:length(x)){
      if(!is.na(x[i-1])&is.na(x[i])&!is.na(x[i+1])) x[i]=0
    }
    x
  }
  rw.mat=apply(rw.mat,2,fix.internal.na)
  # Remove any blank lines at the end of the file,for instance
  rw.mat=rw.mat[!apply(is.na(rw.mat),1,all),]
  rw.df=as.data.frame(rw.mat)
  rw.df
}

