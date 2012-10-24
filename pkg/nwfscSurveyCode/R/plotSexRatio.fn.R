plotSexRatio.fn <-
function(len,fn=median,circleSize=0.1) {
    ratioF <- len$NumF/(len$NumF+len$NumM)
    yF <- lapply(split(ratioF,floor(len$Length)),fn)
    x <- names(split(ratioF,floor(len$Length)))
    nobs <- unlist(lapply(split(ratioF,floor(len$Length)),length))
    plot(x,yF,type="l",col="red")
    symbols(x,yF,circles=nobs,inches=circleSize,fg="red",bg=rgb(1,0,0,alpha=0.5),add=T)
}
