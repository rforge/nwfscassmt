##########
# Plot tows on a map
##########
MapData = function(Data, SA3, strata.limits, FileName, Folder){
  
  # Distilled from line 426-539 of "Survey.Biomass.GlmmBUGS.ver.3.00.R" from John Wallace's code
  jpeg(file=paste(Folder,FileName,"TowMap.jpg",sep=""),width=8,height=8,res=200,units="in")
  
  # Draw box
  plot(c(-55, -1280), c(32, 50.5), xlab = "Depth (m)", ylab = "Latitude", xlim=c(-1280, -55), ylim=c(32, 49), type = "n")
  abline(v= -unique(c(SA3$MIN_DEPTH_M, SA3$MAX_DEPTH_M)), h=unique(c(SA3$MIN_LAT_DD, SA3$MAX_LAT_DD)), col="grey78")
  abline(h=34.5, v=-c(30, 100, 300, 700)*1.8288, col='red')
  
  # Draw centroid of bins
  avelat <- apply(cbind(SA3$MIN_LAT_DD, SA3$MAX_LAT_DD), 1, mean)
  avedep <- apply(cbind(SA3$MIN_DEPTH_M, SA3$MAX_DEPTH_M), 1, mean)
  points(-avedep, avelat, cex=0.5)
  
  # Label areas
  text(-1235, 34.7, "Starting in 2004 NWFSC survey sampling density changes at Pt. Conception (34.5)", adj=0, col='red')
  text(-1235, 49.25, "INPFC Areas", adj=0, col='blue')
  abline(h=c(32, 36, 40.5, 43, 47.5), col='blue')
  text(-1235, 33.25, "Conception", adj=0, col='blue')
  text(-1235, 38.25, "Monterey", adj=0, col='blue')
  text(-1235, 41.75, "Eureka", adj=0, col='blue')
  text(-1235, 44.5, "Columbia", adj=0, col='blue')
  text(-1235, 48.25, "Vancouver", adj=0, col='blue')
  
  # Plot strata
  S = strata.limits
  for (i in 1:nrow(S)) {
    polygon(-c(S$MinDepth[i], S$MaxDepth[i] , S$MaxDepth[i], S$MinDepth[i]), c(S$SLat[i], S$SLat[i], S$NLat[i], S$NLat[i]), col = rainbow(nrow(S), alpha=0.3)[i])
    text(-mean(c(S$MinDepth[i], S$MaxDepth[i])), mean(c(S$SLat[i], S$NLat[i])), S$STRATA[i], cex=1.2)
  }
  
  # Plot Absence tows
  points(-Data$BEST_DEPTH_M[Data$HAUL_WT_KG==0], Data$BEST_LAT_DD[Data$HAUL_WT_KG==0], pch=16, cex=0.5, col=rgb(red=1,0,0,alpha=0.2))
  
  # Plot presence tows by year
  DataPos <- Data[Data$HAUL_WT_KG > 0 & !is.na(Data$HAUL_WT_KG),] # Temp Sp.pos - redefined below
  CU <- c("black", "green", "blue", "cyan", "purple", "grey", "orange", "hotpink", "brown", "darkolivegreen2" ,"darkslategrey",   "deepskyblue1"  , "violet", "cyan", "magenta", "lightsalmon", "gold")
  for(i in 1:length(unique(DataPos$PROJECT_CYCLE))) {
    Which = which(DataPos$PROJECT_CYCLE==unique(DataPos$PROJECT_CYCLE)[i])
    points(-DataPos$BEST_DEPTH_M[Which], DataPos$BEST_LAT_DD[Which], pch=16, cex=0.5, col= CU[i])
  }
  dev.off()
}