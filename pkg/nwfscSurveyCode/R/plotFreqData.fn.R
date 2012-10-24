plotFreqData.fn <-
function(dat,inch=0.15,ylab="Bins",xlab="Year",zero2NAs=T,...) {
    #This function plots frequency data as bubble plots
    #You may want to change all zeros to NA's so that those observations are not plotted.
    #   If you don't then set zero2NAs=F
    x <- as.numeric(as.character(dat$year))
    gender <- dat$gender[1]
    dat <- dat[,-c(1:6)]
    if(gender==0) {
        #dat <- dat[,1:(ncol(dat)/2)]
        dat <- dat[,-match("U0.1",names(dat))]
        #print(names(dat))
        dat <- dat[,-match("U0",names(dat))]
        #print(names(dat))
    }
    if(gender==3) {
        dat <- dat[,-match("F0",names(dat))]
        dat <- dat[,-match("M0",names(dat))]
    }
    numLens <- ncol(dat)/2
    y <- as.numeric(substring(names(dat),2))
    y <- y[1:numLens]
    if(zero2NAs) {dat[dat==0] <- NA}
    
    if(gender==0) {
        z <- c(unlist(dat[,1:numLens]),max(dat))
        symbols(c(rep(x,length(y)),0),c(rep(y,each=length(x)),0),circles=z,main="Unsexed+Males+Females",inches=inch,xlab=xlab,ylab=ylab,xlim=range(x),...)
    }
    if(gender==3) {
        z <- c(unlist(dat[,1:numLens]),max(dat))
        symbols(c(rep(x,length(y)),0),c(rep(y,each=length(x)),0),circles=z,main="Female",inches=inch,xlab=xlab,ylab=ylab,xlim=range(x),...)
        z <- c(unlist(dat[,(numLens+1):ncol(dat)]),max(dat))
        symbols(c(rep(x,length(y)),0),c(rep(y,each=length(x)),0),circles=z,main="Male",inches=inch,xlab=xlab,ylab=ylab,xlim=range(x),...)
    }
}
