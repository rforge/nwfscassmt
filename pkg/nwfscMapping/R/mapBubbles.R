mapBubbles <-
function(dat,depths=NULL,plotMap=T,
            maxValue=max(dat$Z),maxSize=0.3,circle.fg="black",circle.bg=rgb(0,0,0,0.3),type="surface",
            zeroPch=13,zeroCex=0.8,zeroCol=gray(0.5),
            legendPos="bottomright",legendBreaks=NULL,legendType="vert",legendTitle=NULL,legendCex=0.8,legendTitleCex=legendCex,...) 
{
#type=surface is sqrt
    if(plotMap) {
        plotMap(westCoastLL,tck=c(-0.02),...)
        if(exists("westCoastLL")) {
            addLines(westCoastLL)
        } else {
            print("Run data(westCoastLL) before calling this function to plot West Coast shoreline")
        }
        if(exists("WCstatesInlandPBS")) {
            addLines(WCstatesInlandPBS)
        } else {
            print("Run data(WCstatesInlandPBS) before calling this function to plot state boundaries")
        }
        if(!is.null(depths)) {
            for(i in 1:length(depths)) {
                addLines(depths[[i]],col=gray(0.5))
            }
        }
    }
    ind0 <- dat$Z<=0
    ind <- dat$Z>0
    points(dat[ind0,"Longitude"],dat[ind0,"Latitude"],pch=zeroPch,cex=zeroCex,col=zeroCol)
    #symbols(dat[ind0,"Longitude"],dat[ind0,"Latitude"],squares=rep(1,sum(ind0)),add=T,fg=fgZero,bg=fgZero,inches=zeroSize)
    tmp <- dat[ind,]
    event <- as.EventData(data.frame(EID=1:nrow(tmp),X=tmp$Longitude,Y=tmp$Latitude,Z=tmp$Z))
    .addBubbles.ach(event,type=type,z.max=maxValue,max.size=maxSize,symbol.fg=circle.fg,symbol.bg=circle.bg,symbol.zero=F,
                legend.pos=legendPos,legend.breaks=legendBreaks,legend.type=legendType,legend.title=legendTitle,legend.cex=legendCex,legend.title.cex=legendTitleCex)
    #symbols(c(-1000,dat[ind,]$Longitude),c(-1000,dat[ind,]$Latitude),circles=sqrt(c(maxValue,dat[ind,]$Z)),inches=maxSize,add=T,fg=circle.fg,bg=circle.bg)
    
    invisible()
}
