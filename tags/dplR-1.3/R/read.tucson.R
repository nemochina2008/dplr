`read.tucson` <-
function(fname, header=NULL,long=FALSE)
{
  ## Open the data file for reading
  dat=file(fname,"r")
  if(is.null(header)){
    ## Try to determine if the file has a header. This is failable.
    ## 3 lines in file
    hdr1=readLines(dat,n=1)
    if(length(hdr1)==0){
      close(dat)
      stop("File is empty")
    }
    if(nchar(hdr1) < 12){
      close(dat)
      stop("First line in .rwl file ends before col 12")
    }
    yrcheck=suppressWarnings(as.numeric(substr(hdr1,9,12)))
    if(is.null(yrcheck) || length(yrcheck)!=1 || is.na(yrcheck) ||
       yrcheck < -1e04 || yrcheck > 1e04) {
      cat("There appears to be a header in the rwl file\n")
      is.head=TRUE
    }
    else {
      is.head=FALSE # No header lines
      cat("There does not appear to be a header in the rwl file\n")
    }
    seek(dat, where=0, origin="start")
  } else if(!is.logical(header)){
    close(dat)
    stop("header must be NULL or logical")
  } else{
    is.head = header
  }
  if(is.head){
    ## Read 4th line - should be first data line
    dat1=readLines(dat,n=4)
    if(length(dat1)<4){
      close(dat)
      stop("File has under 4 lines")
    }
    dat1=dat1[4]
  } else{
    dat1=readLines(dat,n=1)
    if(length(dat1)==0){
      close(dat)
      stop("File is empty")
    }
  }
  close(dat)
  yrcheck=as.numeric(substr(dat1,9,12))
  if(is.null(yrcheck) || length(yrcheck)!=1) {
    stop("Cols 9-12 of first data line not a year")
  }

  skip.lines=ifelse(is.head,3,0)
  if(long){
    dat=read.fwf(fname,c(7,5,rep(6,10)),skip=skip.lines,strip.white=TRUE,
      blank.lines.skip=TRUE)
  }
  else{
    dat=read.fwf(fname,c(8,5,rep(6,10)),skip=skip.lines,strip.white=TRUE,
      blank.lines.skip=TRUE)
  }
  ## Remove any blank lines at the end of the file, for instance
  dat=dat[!apply(is.na(dat),1,all),]

  series=dat[,1]
  series.ids=unique(series)
  nseries=length(series.ids)

  cat("There are ",nseries," series\n",sep="")

  series.index=match(series,series.ids)
  min.year=(min(dat[,2]) %/% 10) * 10
  max.year=((max(dat[,2])+10) %/% 10) * 10
  span=max.year - min.year + 1

  rw.vec=NA*vector(mode="numeric", length=nseries*span)
  series.min=rep.int(as.integer(max.year+1), nseries)
  series.max=rep.int(as.integer(min.year-1), nseries)
  prec.rproc=rep.int(as.integer(1), nseries)
  x=as.matrix(dat[,-c(1,2)])
  decade.yr=dat[,2]
  .C(rwl.readloop,
     series.index,
     decade.yr,
     as.vector(x, mode="integer"),
     nrow(x),
     ncol(x),
     as.integer(min.year),
     rw.vec,
     as.integer(span),
     as.integer(nseries),
     series.min,
     series.max,
     prec.rproc,
     NAOK=TRUE, DUP=FALSE)
  cat(paste(1:nseries, "\t",
            series.ids, "\t",
            series.min, "\t",
            series.max, "\t",
            1/prec.rproc, "\n", sep=""), sep="")
  rw.mat=matrix(rw.vec,ncol=nseries,nrow=span)
  rownames(rw.mat)=min.year:max.year
  
  ## Convert values < 0 to NA (either precision)
  rw.mat[rw.mat<0]=NA
  ## The operations in the loop depend on the precision of each series.
  ## It's not exactly clear whether the Tucson format allows mixed
  ## precisions in the same file, but we can support that in any case.
  for(i in 1:nseries){
    this.prec.rproc = prec.rproc[i]
    if(this.prec.rproc == 100){
      ## Convert stop marker (and any other) 999 to NA (precision 0.01)
      rw.mat[rw.mat[,i]==999,i]=NA
    } else if(this.prec.rproc != 1000){
      stop(paste("Precision unknown in series", series.ids[i]))
    }
    ## Convert to mm
    rw.mat[,i] = rw.mat[,i] / this.prec.rproc
  }

  ## trim the front and back of the output file to remove blank rows. The
  ## ham-handed approach below avoids removing areas with internal gaps such as
  ## those with floating, but dated segments.
  ##rw.mat=rw.mat[!apply(is.na(rw.mat),1,all),]
  ## subset first 11 years to trim out leading NAs
  foo <- as.matrix(rw.mat[1:11,])
  foo.yrs <- as.numeric(rownames(foo))
  min.year0 <- min(foo.yrs[!apply(is.na(foo),1,all)])
  ## subset last 11 years to trim out ending NAs
  foo <- as.matrix(rw.mat[(nrow(rw.mat)-11):nrow(rw.mat),])
  foo.yrs <- as.numeric(rownames(foo))
  max.year0 <- max(foo.yrs[!apply(is.na(foo),1,all)])
  ## trim
  yrs <- min.year0:max.year0
  rw.mat=as.matrix(rw.mat[as.numeric(rownames(rw.mat)) %in% yrs,])
  ## Fix internal NAs. These are coded as 0 in the DPL programs
  fix.internal.na=function(x){
    for(i in 2:(length(x)-1)){
      if(!is.na(x[i-1])&is.na(x[i])&!is.na(x[i+1])) x[i]=0
    }
    x
  }
  rw.mat=apply(rw.mat,2,fix.internal.na)
  rw.df=as.data.frame(rw.mat)
  colnames(rw.df)=as.character(series.ids)
  rw.df
}
