##############################################################################################
### Plot latitudinal and depth distributions of positive tows, along with percent positive ###
##############################################################################################
#Inputs
#Survey.dat: Survey data used in the GLMM function. This can contain multiple species.
#strata.limits.spp: The strata used in the GLMM function. This can either be one entry or a list of entries for multiple species.
#spp.strata.vec: This vector indicates which of the species included in the Survey.dat file are also included int he strata.limits.spp file. 0 = not included, 1 = included.
#spp.col.start: the column in which species biomass begins
#sep.yr: Plot each year separate (=T) or all on one plot (=F)

Survey.CPUE.strata<-function(Survey.dat,strata.limits.spp,spp.strata.vec=c(0,0,0,0,1,1,0,1,1,1),spp.col.start=8,sep.yr=F)
{
  sep.yr=F
  survey.years<-sort(unique(Survey.dat$YEAR))
  col.bubs<-rainbow(length(survey.years))
  Survey.dat.binom<-Survey.dat
  j<-k<-1
  for(ii in spp.col.start:ncol(Survey.dat))
  {
    par(mar=c(5,2,4,5))
    Survey.dat.binom[Survey.dat[,ii]>0,ii]<-1
    if(sep.yr==F)
    {
      bubble3xy(-subset(Survey.dat,YEAR==survey.years[1])$BEST_DEPTH_M,subset(Survey.dat,YEAR==survey.years[1])$BEST_LAT_DD,subset(Survey.dat,YEAR==survey.years[1])[,ii]/subset(Survey.dat,YEAR==survey.years[1])$AREA_SWEPT_MSQ,xlab="Depth (m)",ylab="",xlim=c(-800,0),ylim=c(32,50),col=c(col.bubs[1],col.bubs[1]),axis1=FALSE,axis2=FALSE,main=paste(unlist(strsplit(colnames(Survey.dat)[ii],"_")),collapse=" "))
      for(i in 2:length(survey.years))
      {bubble3xy(-subset(Survey.dat,YEAR==survey.years[i])$BEST_DEPTH_M,subset(Survey.dat,YEAR==survey.years[i])$BEST_LAT_DD,subset(Survey.dat,YEAR==survey.years[i])[,ii]/subset(Survey.dat,YEAR==survey.years[i])$AREA_SWEPT_MSQ,col=c(col.bubs[i],col.bubs[i]),add=TRUE)}
      axis(1,seq(-800,0,50),rev(seq(0,800,50)))
      axis(4,seq(32,50,2),seq(32,50,2),las=2)
      mtext("Latitude (degrees)",4,line=3,las=0)
      smartlegend("left","top",survey.years,pch=1,col=col.bubs,bty="n",)  
      text(-750,32,paste("%+ = ",round(sum(Survey.dat.binom[,ii])/length(Survey.dat.binom[,ii]),2)*100,"%",sep=""),cex=1.25)
      if(spp.strata.vec[ii-(spp.col.start-1)]==1)
      {
        if(is.data.frame(strata.limits.spp)==TRUE){for(jj in 1:nrow(strata.limits.spp)){rect(-strata.limits.spp[jj,5],strata.limits.spp[jj,3],-strata.limits.spp[jj,4],strata.limits.spp[jj,2],lwd=2)}}
        if(is.data.frame(strata.limits.spp)==FALSE)
        {
          for(jj in 1:nrow(strata.limits.spp[[j]])){rect(-strata.limits.spp[[j]][jj,5],strata.limits.spp[[j]][jj,3],-strata.limits.spp[[j]][jj,4],strata.limits.spp[[j]][jj,2],lwd=2)}
          j<-j+1
        }
      }  
    }
    
    if(sep.yr==T)
    {
      for(i in 1:length(survey.years))
      {
        bubble3xy(-subset(Survey.dat,YEAR==survey.years[i])$BEST_DEPTH_M,subset(Survey.dat,YEAR==survey.years[i])$BEST_LAT_DD,subset(Survey.dat,YEAR==survey.years[i])[,ii]/subset(Survey.dat,YEAR==survey.years[i])$AREA_SWEPT_MSQ,xlab="Depth (m)",ylab="",xlim=c(-800,0),ylim=c(32,50),col=c(col.bubs[i],col.bubs[i]),axis1=FALSE,axis2=FALSE,main=paste(paste(unlist(strsplit(colnames(Survey.dat)[ii],"_")),collapse=" "),": ",survey.years[i],sep=""))
        axis(1,seq(-800,0,50),rev(seq(0,800,50)))
        axis(4,seq(32,50,2),seq(32,50,2),las=2)
        mtext("Latitude (degrees)",4,line=3,las=0)
        text(-750,32,paste("%+ = ",round(sum(subset(Survey.dat.binom,YEAR==survey.years[i])[,ii])/length(subset(Survey.dat.binom,YEAR==survey.years[i])[,ii]),2)*100,"%",sep=""),cex=1.25)
        if(spp.strata.vec[ii-(spp.col.start-1)]==0){strata<-matrix(-999,2,5)}
        if(spp.strata.vec[ii-(spp.col.start-1)]==1){strata<-strata.limits.spp[[k]]}
        for(jj in 1:nrow(strata)){rect(-strata[jj,5],strata[jj,3],-strata[jj,4],strata[jj,2],lwd=2)}      
      }
      if(spp.strata.vec[ii-(spp.col.start-1)]==1){k<-k+1}  
    }
  }  
}

#Example
Survey.CPUE.strata(subset(Dmod.tri.shelf.dat,YEAR>1977 & YEAR<1995),strata.limits.spp.earlytri,spp.strata.vec=c(0,0,0,0,1,1,0,1,1,1),spp.col.start=8,sep.yr=F)
