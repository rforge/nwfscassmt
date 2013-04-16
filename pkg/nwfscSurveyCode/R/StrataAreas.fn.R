StrataAreas.fn <- function(strat.df,SA3,convertFactor=0.01) {
    #calculates the area of your strata using John Wallace's SA3 file
    #this code is stolen from within John Wallace's 2011 GLMM code
    #a convertFactor of 0.01 convert hectares to km2
    S <- strat.df
    S$area <- NA

    for ( i in 1:nrow(S)) {
        maxLat <- max(c(S$START_LATITUDE.1[i],S$START_LATITUDE.2[i])) >= SA3$MAX_LAT_DD
        minLat <- min(c(S$START_LATITUDE.1[i],S$START_LATITUDE.2[i])) <= SA3$MIN_LAT_DD
        if(sum(maxLat) == 0 | sum(minLat) == 0) stop("A latitude in your strata is not available in SA3.\nEither use ana available latitude or supply your own area.") 
        maxDep <- max(c(S$BOTTOM_DEPTH.1[i],S$BOTTOM_DEPTH.2[i])) >= SA3$MAX_DEPTH_M
        minDep <- min(c(S$BOTTOM_DEPTH.1[i],S$BOTTOM_DEPTH.2[i])) <= SA3$MIN_DEPTH_M
        if(sum(maxDep) == 0 | sum(minDep) == 0) stop("A depth in your strata is not available in SA3.\nEither use an available depth or supply your own area.") 
        R <- SA3[maxLat & minLat & maxDep & minDep,]
        S$area[i] <- sum(R$AREA_HECTARES)*convertFactor
    }
    S
}
