PosteriorPredictive = function(Data, Model, maxDims=6, FileName, Folder=NA){
  if(is.na(Folder)) Folder = paste(getwd(),"/",sep="")  
  
  # Attach stuff
  attach(Model$BUGSoutput$sims.list)
  #attach(Data)
  Dist = Model$likelihood
  nonZeros = which(Data[,'isNonZeroTrawl']==TRUE)
  u.nz = matrix(NA, nrow=nrow(B.pos), ncol=length(nonZeros))
  
  # Warnings about mismatch between data and model
  if(nlevels(strata)!=ncol(cMx(Sdev)) | nlevels(year)!=ncol(Ydev) | nlevels(strataYear)!=ncol(SYdev)| nlevels(vesselYear)!=ncol(VYdev)){
    stop("Model and data do not match")
  }
  
  if(Dist=="lognormal"){
    sigma = sqrt(log(1+(CV[,1]^2)))
    for(i in 1:length(nonZeros)){
      u.nz[,i] <- exp( cMx(Sdev)[,strata[nonZeros[i]]] + Ydev[,year[nonZeros[i]]] + VYdev[,vesselYear[nonZeros[i]]] + SYdev[,strataYear[nonZeros[i]]] + B.pos[,1]*logeffort[nonZeros[i]] + B.pos[,2]*logeffort2[nonZeros[i]] )
    }
  }
  if(Dist=="gamma"){
    gamma.a = oneOverCV2 = 1/(CV[,1]^2)
    for(i in 1:length(nonZeros)){
      Temp = cMx(Sdev)[,strata[nonZeros[i]]] + Ydev[,year[nonZeros[i]]] + VYdev[,vesselYear[nonZeros[i]]] + SYdev[,strataYear[nonZeros[i]]] + B.pos[,1]*logeffort[nonZeros[i]] + B.pos[,2]*logeffort2[nonZeros[i]]
      u.nz[,i] <- exp( ifelse(Temp>100,100,Temp) ) # Eric included this in the JAGS code
    }
  }
  if(Dist=="invGaussian"){
    oneOverCV2 = 1/(CV[,1]^2)
    for(i in 1:length(nonZeros)){
      Temp = cMx(Sdev)[,strata[nonZeros[i]]] + Ydev[,year[nonZeros[i]]] + VYdev[,vesselYear[nonZeros[i]]] + SYdev[,strataYear[nonZeros[i]]] + B.pos[,1]*logeffort[nonZeros[i]] + B.pos[,2]*logeffort2[nonZeros[i]]
      u.nz[,i] <- exp( ifelse(Temp>100,100,Temp) ) # Eric included this in the JAGS code
    }
  }
  if(Dist=="lognormalECE"){
    u.nz2 = array(NA, dim=dim(u.nz))
    sigma = sqrt(log(1+(CV[,1]^2)))
    sigma2 = sqrt(log(1+(CV[,2]^2)))
    for(i in 1:length(nonZeros)){
      Temp = cMx(Sdev)[,strata[nonZeros[i]]] + Ydev[,year[nonZeros[i]]] + VYdev[,vesselYear[nonZeros[i]]] + SYdev[,strataYear[nonZeros[i]]] + B.pos[,1]*logeffort[nonZeros[i]] + B.pos[,2]*logeffort2[nonZeros[i]]
      u.nz[,i] <- exp( ifelse(Temp>100,100,Temp) ) # Eric included this in the JAGS code
    }
    for(i in 1:length(nonZeros)){
      u.nz2[,i] <- u.nz[,i] * ratio
    }
  }
  if(Dist=="gammaECE"){
    u.nz2 = array(NA, dim=dim(u.nz))
    gamma.a = 1/(CV[,1]^2)
    gamma.a2 = 1/(CV[,2]^2)    
    for(i in 1:length(nonZeros)){
      Temp = cMx(Sdev)[,strata[nonZeros[i]]] + Ydev[,year[nonZeros[i]]] + VYdev[,vesselYear[nonZeros[i]]] + SYdev[,strataYear[nonZeros[i]]] + B.pos[,1]*logeffort[nonZeros[i]] + B.pos[,2]*logeffort2[nonZeros[i]]
      u.nz[,i] <- exp( ifelse(Temp>100,100,Temp) ) # Eric included this in the JAGS code
    }
    for(i in 1:length(nonZeros)){
      u.nz2[,i] <- u.nz[,i] * ratio
    }
  }
  
  Q = rep(NA, nrow(Data)) # vector to track quantiles for each observation
  Nstrat = length(unique(strataYear[nonZeros]))
  Ncol=ceiling(sqrt(Nstrat)); Nrow = ceiling(Nstrat/Ncol)
  Nstrat = length(unique(strataYear[nonZeros]))
  Nplots = ceiling( Nstrat / maxDims^2 )
  for(PlotI in 1:Nplots){
    ParSet = (PlotI-1)*maxDims^2 + 1:min(Nstrat-(PlotI-1)*maxDims^2,maxDims^2)
    Ncol=ceiling(sqrt(length(ParSet)))
    Nrow = ceiling(length(ParSet)/Ncol)
    jpeg(paste(Folder,"/",FileName,"Posterior_Predictive_",PlotI,".jpg",sep=""),width=Ncol*1.5,height=Nrow*1.5,res=200,units="in")
    par(mfrow=c(Nrow,Ncol), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
    for(StrataYearI in ParSet){
      Which = which(strataYear[nonZeros]==unique(strataYear[nonZeros])[StrataYearI])
      plot(Data[nonZeros[Which],'HAUL_WT_KG'],ylab="",xlab="",log="y",main=unique(strataYear[nonZeros])[StrataYearI],col="blue", ylim=range(Data[nonZeros,'HAUL_WT_KG']))
      # mean(u.nz[,2])
      for(ObsI in 1:length(Which)){
        if(Dist=="lognormal"){     
          y = rlnorm(n=1000,meanlog=log(u.nz[,Which[ObsI]]),sdlog=sigma)   # Plotting in log-space
        }
        if(Dist=="gamma"){     
          b = gamma.a / u.nz[,Which[ObsI]];    
          y = rgamma(n=1000,shape=gamma.a,rate=b)
        }
        if(Dist=="invGaussian"){     
          lambda = u.nz[,Which[ObsI]]*oneOverCV2
          y = rinvgauss(n=1000,mu=u.nz[,Which[ObsI]],lambda=lambda)
        }
        if(Dist=="lognormalECE"){     
          ECE = rbinom(n=nrow(u.nz)*10, size=1, prob=p.ece[,2])
          y = rlnorm(n=1000, meanlog=log(u.nz[,Which[ObsI]])*(1-ECE)+log(u.nz2[,Which[ObsI]])*ECE, sdlog=sigma*(1-ECE)+sigma2*ECE)
        }
        if(Dist=="gammaECE"){     
          b = gamma.a / u.nz[,Which[ObsI]];    
          b2 = gamma.a2 / u.nz2[,Which[ObsI]];    
          ECE = rbinom(n=nrow(u.nz)*10, size=1, prob=p.ece[,2])
          y = rgamma(n=1000, shape=gamma.a*(1-ECE)+gamma.a2*ECE, rate=b*(1-ECE)+b2*ECE)
        }
        Q[nonZeros[Which[ObsI]]] = mean(y>Data[nonZeros[Which[ObsI]],'HAUL_WT_KG'])
        Quantiles = quantile(y,prob=c(0.025,0.25,0.75,0.975))
        lines(x=c(ObsI,ObsI), y=Quantiles[2:3], lwd=2)
        lines(x=c(ObsI,ObsI), y=Quantiles[c(1,4)], lwd=1,lty="dotted")
        if(Data[nonZeros[Which[ObsI]],'HAUL_WT_KG']>max(Quantiles) | Data[nonZeros[Which[ObsI]],'HAUL_WT_KG']<min(Quantiles)){
          points(x=ObsI,y=Data[nonZeros[Which[ObsI]],'HAUL_WT_KG'],pch=4,col="red",cex=2)
        }
      }
    }
    dev.off()
  }
  
  # Q-Q plot
  jpeg(paste(Folder,"/",FileName,"Q-Q_plot.jpg",sep=""),width=4,height=4,res=200,units="in")
  par(mfrow=c(1,1), mar=c(2,2,2,0), mgp=c(1.25,0.25,0), tck=-0.02)
  Qtemp = na.omit(Q)
  Order = order(Qtemp)
  plot(x=seq(0,1,length=length(nonZeros)), y=Qtemp[Order], main="Q-Q plot", xlab="Uniform", ylab="Empirical")
  abline(a=0,b=1)
  dev.off()
  
  # Detach stuff
  #detach(Data)
  detach(Model$BUGSoutput$sims.list)
}