multiMap <-
function(nrows=1,ncols=1,lonLims,latLims,adjLonRange=TRUE,fnc,...) {
    #lonLims and latLims are matrices with each row is a panel and each column is min and max
    numMaps <- nrows*ncols
    if(ncol(latLims)!=2) stop("latLims has to have two columns for min latitude and max latitude")
    if(ncol(lonLims)!=2) stop("lonLims has to have two columns for min longitude and max longitude")
    if(nrow(latLims)!=numMaps) stop("Need 1 row in latLims for each map panel")
    if(nrow(lonLims)!=numMaps) stop("Need 1 row in lonLims for each map panel")

    #find a common range
    midLat <- apply(latLims,1,mean)
    midLon <- apply(lonLims,1,mean)
    dis <- rep(NA,numMaps)
    if(adjLonRange) {
        for(i in 1:numMaps) {
            dis[i] <- distance(lonLims[i,1],midLat[i],lonLims[i,2],midLat[i])
        }
        for(i in 1:numMaps) {
            diffDis <- max(dis) - dis[i]
            lonLims[i,1] <- uniroot(function(x,lon2,lat){distance(x,lat,lon2,lat)-diffDis/2},c(-180,lonLims[i,1]),lon2=lonLims[i,1],lat=midLat[i])$root
            lonLims[i,2] <- uniroot(function(x,lon2,lat){distance(x,lat,lon2,lat)-diffDis/2},c(0,lonLims[i,2]),lon2=lonLims[i,2],lat=midLat[i])$root
        }
    }
    if(!adjLonRange) {
        stop("adjust latitude not implemented yet\n")
        #d <- apply(tmp,1,diff)
        #ind <- which(d!=max(d))
        #for(i in ind) {  #adjust the lims that are smaller than max
        #    tmp[i,1] <- tmp[i,1] - (max(d)-d[i])/2
        #    tmp[i,2] <- tmp[i,2] + (max(d)-d[i])/2
        #}
        #if(adjLonRange) lonLims <- tmp
        #if(!adjLonRange) latLims <- tmp
    }
    
    par(mfrow=c(nrows,ncols))
    x <- apply(latLims,1,diff)
    x <- which(x==max(x))
    fnc(xlim=lonLims[x[1],],ylim=latLims[x[1],],...)
    plt1 <- par()$plt
    par(mfrow=c(nrows,ncols))
    for(i in 1:numMaps){
        fnc(xlim=lonLims[i,],ylim=latLims[i,],plt=NULL,...)  #use plt=NULL or plt=par()$plt so that same as first plot
        par(plt=plt1)
    }
    invisible()
}
