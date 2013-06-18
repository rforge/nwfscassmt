################################################################################
# This is an attempt to code George Sugihara's simplex method of non linear
#time series analysis in R.
#
#             DRH       5/4/12
#_______________________________________________________________________________
E.simplex=function(
    dat=c()
    ,X
    ,Y
    ,E
    ,plot=T
    ,tau1=1
    ,cross.val=F
    ,diag1=F
    ,file.name="E.simplex"
    ,path=getwd()
    ) {

  # All the data, (library and prediction sets) should be stored in dat.
  #The library vector is in X and the prediction vector is in Y.  E is the
  #vector of embedding dimensions you want to test.  Plot = T is the default
  #setting that will make some plots of output.  tau1 is the lag step for
  #embedding (default = 1). Tau is also the prediction window - there may be a
  #way to separate these later on so that tau and the prediction window are
  #not necessarily equivalent.  diag1 is a switch that turns on various diagnostic
  #plots.  file.name will store the prediction vector and other results - this
  #should have a path name if you want this to go somewhere other than the
  #working directory.  The working directory can be specified by path as well
  #other wise it will default to the current working directory.
#_______________________________________________________________________________

p1 <- try(setwd(path))
if (class(p1)=="try-error") setwd(getwd())
print(paste("Your working directory is: ",getwd(),sep=""))

#move into E[i] dimensional state space and for each Y find the E[i}+1
#nearest neighbors in the various lagged states of X.  The nearest neighbors
#depend critically on the dimension of the state space.  First we need to
#make the set of vectors that describes each point in state space.  The
#number of points in state space will decrease as E[i] increases
new.plot<-1
save.plot<-3

  for(i in 1:length(E)){

    #every three iterations we need a new plot window
    if(i==new.plot) {
      set.plot.window()
      new.plot<-new.plot+3
      }

    #A switch is required here - if the prediction set doesn't exist
    #then the cross validation procedure should be switched on.
    if(length(Y)==0) {
      cross.val<-T
      #Also need a test to determine if the dimensionality requested by the user
      #exceeds the constraints imposed by the length of the time series
      gappers<-((E[i]-1)+(E[i]-2)*(tau1-1))
      incompletes<-((E[i]-1)*tau1)+1
      check1<-(length(X)-1)
      check1<-check1 - gappers - incompletes  #vectors that cross the time gap
      #given a cross validation procedure
      check2<-(E[i]+tau1)
      if( check1 <= check2 ) {
        print("$$$$$ ERROR - Insufficient sample size for required dimensionality $$$$$")
        print(paste("n = ",length(X),", E = ",E[i],", tau = ",tau1,sep=""))
        return("EXITING FUNCTION E. SIMPLEX")
      }
    }

    #Going to add an overide switch if the user wants to run cross validation
    #despite having a viable prediction set...
    if(cross.val==T) {
      Pred<-cross.val.pred(X,E[i],tau1,diag1)
      #the cross validation procedure can use the same architecture as the
      #traditional library/prediction sets, it just need the library and
      #prediction sets to be built iteratively and then carried forward
      } else {    #this is the traditional set up...

        #First check the dimensionality of the problem - can we proceed?
        if( length(Y) <= E[i]*tau1 ) {
          print("$$$$$ ERROR - Insufficient sample size for required dimensionality $$$$$$")
          print(paste("Prediction vector n = ",length(Y),
            ", E = ",E[i],", tau = ",tau1,sep=""))
          return("EXITING FUNCTION E. SIMPLEX")
        }

        #Here we create the library and prediction sets
        LSS<-c()
        PSS<-c()
        LSS<-embed.library(X,LSS,E[i],tau1)
        PSS<-embed.predictor(Y,PSS,E[i],tau1)


         #now we have the library and prediction vectors filled and named
         #______________________________________________________________________

         #The next step is to find the nearest neighbors - R has a built in
         #Euclidian distance measure - dist().  We just need to iterate through
         #the prediction data frame
         nearest<-matrix(nrow=dim(PSS)[1],ncol=(E[i]+1))
         for (j in 1:dim(PSS)[1]) {
            nearest[j,]<-find.mins(PSS[j,],LSS,(E[i]+1),tau1)
         }
         #nearest now holds the row numbers for the E[i] closest neighbor in LSS
         #of each row in PSS.
         #______________________________________________________________________

         #Call the function that makes the predictions...
         Pred<-predictions(E[i],LSS,PSS,nearest,tau1,diag1)
         #______________________________________________________________________
      } #end of cross validation boolean switch

   #Now we can compare the Prediction to the prediction vector and see how we
   #did - first make plot variables
   plot.x<-c()
   plot.y<-c()
   if(cross.val==F) {
      plot.x<-dat
      #filling the prediction vector takes some doing - we must stop predicting
      #at the last spot - (E[i]-1).  Otherwise you run out of room to lag...
      plot.y[1:length(X)+(E[i]-1)]<-NA
      plot.y[(length(X)+E[i]):(length(X)+(E[i]-1)+length(Pred))]<-
            Pred[1:(length(Pred))]
    }

   if(cross.val==T) {
      plot.x<-X        #There is no possible comparison for the tau1 first
      #observed values, since you can't make a prediction for them, but we will
      #always plot it anyway (a check on the vector position assignments).
      plot.y<-Pred    #Likewise there is possible comparison for the last
      #prediction becasue there is no matching observation - till next year!
    }

   #plot the model fits and the observed values on the time series scale
   if(cross.val==F) plot.series(plot.y,plot.x,E[i],tau1)
   if(cross.val==T) plot.series.cross(plot.y,plot.x,E[i],tau1)
   #plot the correlation between observed and predicted
   plot.cor(plot.x,plot.y,E[i],tau1)

   #calculate some stats for performance of the predictor
   rho<-cor(plot.y,plot.x,use="pairwise.complete.obs")
   mae<-mean(na.omit(plot.y-plot.x))
   rmse<-sqrt(mean( (na.omit(plot.y-plot.x))^2) )

   line1<-"__________________________________________________________________"
   print(line1)
   print(paste("n = ",length(X),"; prediction vector n = ",length(Y),"; E = "
      ,E[i],"; tau = ",tau1,sep=""))
   print(paste("rho = ",rho,"; MAE = ",mae,"; RMSE = ",rmse,sep=""))
   rez<-data.frame(plot.y,plot.x)
   names(rez)<-c("Predicted","Observed")
   print(rez) #some output that should go into an output file...

   #store some results in the output file
   app<-FALSE
   if(i>1) app<-TRUE
   sampl<-data.frame(length(X),length(Y),E[i],tau1)
   names(sampl)<-c("obs. vec. n","pred. vec. n","E","tau")
   write.table(sampl
      ,file=paste(file.name,"csv",sep=".")
      ,sep=","
      ,row.names=F
      ,append=app
    )
   output<-data.frame(rho,mae,rmse)
   write.table(output
      ,file=paste(file.name,"csv",sep=".")
      ,sep=","
      ,row.names=F
      ,append=TRUE
    )
   write.table(rez
      ,file=paste(file.name,"csv",sep=".")
      ,sep=","
      ,row.names=F
      ,append=TRUE
    )
   write.table(line1
      ,file=paste(file.name,"csv",sep=".")
      ,row.names=F
      ,col.names=F
      ,append=TRUE
    )


    #every three iterations we need to save the plot window
    if(i==save.plot) {
      savePlot(filename=paste(file.name,(save.plot-2),"pdf",sep="."),type="pdf")
      save.plot<-save.plot+3
      }
    if(i==length(E)){
        savePlot(filename=paste(file.name,(save.plot-2),"pdf",sep="."),type="pdf")
      }

  } #E loop


}  #close function E.simplex

################################################################################

################################################################################
embed.library=function(x,LS,E1,ta1) {
#Here we create the library sets, by making lagged vectors as
#needed. x holds the original library set
#LS will be the matrix of embedded library vectors
#E is the embedding dimension, ta1 is the times step
  LS<-c()
  LS<-data.frame(x[E1:length(x)])
  names(LS)[1]<-"Xt"
  for(j in 1:(E1-1)) {
    LS<-data.frame(LS,x[(E1-j*ta1):(length(x)-j)])
    names(LS)[j+1]<-paste("Xt-",j,sep="")
  }
  return(LS)
}
################################################################################

################################################################################
embed.predictor=function(y,PS,E1,ta1) {
#Here we create the prediction sets, by making lagged vectors as
#needed. y holds the original prediction set
#PS will be the embedded prediction vectors,
#E is the embedding dimension, ta1 is the times step
  PS<-c()
  PS<-data.frame(y[E1:length(y)])
  names(PS)[1]<-"Xt"
  for(j in 1:(E1-1)) {
    PS<-data.frame(PS,y[(E1-j*ta1):(length(y)-j)])
    names(PS)[j+1]<-paste("Xt-",j,sep="")
  }
  return(PS)
}
################################################################################

################################################################################
find.mins=function(spot,mat,mins,tau1){
# This function will find the E minimum Euclidian distances between spot and
# the rows of mat.  It will return the row numbers of these nearest neighbors
# The prediction based on these distances will itself be based on the
# exponential weighted distance at time t and applied to time t+1. Therefore
# you cannot use the last tau1 observations in the library set, becasue time
# t+tau1 doesn't exist!  Best to take the last tau1 rows out of the competition.
  mat<-mat[-(((dim(mat)[1]+1)-tau1):dim(mat)[1]),]
  names(spot)<-names(mat)
  dvec<-dist(rbind(spot,mat)) #distances calculated here
  min1<-order(as.numeric(as.matrix(dvec)[-1,1]))  #order them
  m<-c()
  m<-min1[1:mins]
  return(m)  #return the pointers...
} #close function find.mins
#_______________________________________________________________________________

################################################################################
predictions=function(E1,LS,PS,nearest1,tau1,diag1) {
   #We have to estimate what the prediction of the library set
   #is and compare it to what actually ocurred. The prediction is based on the
   #exponential weighted distance at time t and applied to time t+1

   Pred<-c()
   for (j in 1:(dim(PS)[1])) {                #E1:(dim(PS)[1]-1)) {
      wt<-c() #initialize the wt vector
      for(k in 1:(E1+1)) {   #first determine the exponential weights.  These are
          #scaled to the distance between the nearest of the neighbors
          numer<-as.numeric(as.matrix(dist(rbind(LS[nearest1[j,k],],PS[j,])))[2])
          denom<-as.numeric(as.matrix(dist(rbind(LS[nearest1[j,1],],PS[j,])))[2])
          wt[k]<-exp(-(numer/denom) )
      }
      #now get the sum
      wt.sum<-sum(wt)
      #The prediction is now the weighted average of the library vectors at
      #time t+1
      wt.avg<-0
      for(k in 1:(E1+1)) {
          wt.avg<-wt.avg+(wt[k]*LS[(nearest1[j,k]+tau1),1])   #apply the
          #weighting to the point one time step in the future
        }
    Pred[j]<-wt.avg/wt.sum

    #graphical exploration of the prediction routine
    #print(PS)

    if(diag1==T) plot.prediction.graphics(PS,LS[nearest1[j,],],LS[(nearest1[j,]+tau1),],LS)


    #error trap
      if(is.finite(Pred[j])==F) {
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        print(c(" $$$ERROR - non finite prediction $$$$ ",wt.avg,wt.sum))
        print(c("n = ",dim(LS)[1]))
        print(c("weighting vector = ",wt))
        print(c("nearest neghbor pointers are: ",(nearest1[j,]+tau1)))
        print(c("nearest neighbor values: ",LS[(nearest1[j,]+tau1),1]))
        print(c("prediction vector is: ",as.data.frame(PS)))
        print(c("observation matrix is: ",as.data.frame(LS)))
        print("$$$END ERROR REPORT source: predictions() $$$")
        print("______________________________________________________")
      }
    }
  return(Pred)
  } # close function predictions
#_______________________________________________________________________________


################################################################################
cross.val.pred=function(X,E1,tau1,diag1,...) {
  #This function iteratively assigns prediction and library vectors in a cross
  #validation procedure.  It carries the algorithm through to a final prediction
  #vector based on simplex projection and calling other functions defined above.
  Predict<-c()
  #Remember that holding the first value out of sample doesn't work.  There
  #is no X(t-tau1) to base the predictions on - so the removals will start at
  #E1 + tau1 (minimally 3 in a two dimensional simplex).  Note that it is
  #impossible to create a vector of length E1 from observations before the
  #E1+tau-1 position in the prediction set (where the lag = tau) -
  #without spanning the gap, which is a no-no.
  off.set<-((E1-1)*tau1)
  for(i in (off.set+1):length(X)) {
    LSS<-c()
    PSS<-c()
    X2<-c() #we want to make sure to preserve the original data so make a temp..
    X2<-data.frame(X[1:length(X)]) #build the library set
    PSS<-data.frame(X2[i,1]) #pick an observation far enough into the library set
    #that you can make a full orediction vector, while accounting for the
    #current lag and put it into the prediction set
    LSS<-as.data.frame(X2[-i,1]) #then remove the prediction point
    names(LSS)[1]<-"Xt"
    names(PSS)[1]<-"Xt"


    #now create the lagged columns
    for(j in 1:(E1-1)) {
      #initialize the lag column
      LSS.lag<-c()
      LSS.lag<-LSS[1:(dim(LSS)[1]-(tau1)),j] #fill it with lagged observation
      LSS<-LSS[-(1:tau1),1:j] #drop the first rows to match lagged vector
      LSS<-data.frame(LSS,LSS.lag)
      PSS<-data.frame(PSS,X[(i-(j*tau1))]) #lagged prediction set
    }

    #fix dataframe the names up again...
    names(LSS)[1]<-"Xt"
    names(PSS)[1]<-"Xt"
    for(j in 2:E1) {
        names(LSS)[j]<-paste("Xt-",j-1,sep="")
        names(PSS)[j]<-paste("Xt-",j-1,sep="")
      }

    #after creating the lagged data frame..
    gap1<-(E1-1)+((E1-2)*(tau1-1))
    LSS<-LSS[-((i-off.set):(i-off.set+gap1)),] #drop the rows that span the gap in the time
    #series - the sequence must be adjusted down one row due to the removals
    #above...

    #______________________________________________________________________

    #Find the nearest neigbors for each coordinate vector
    nearest<-matrix(nrow=dim(PSS)[1],ncol=(E1+1))
    for (j in 1:dim(PSS)[1]) {
        nearest[j,]<-find.mins(PSS[j,],LSS,(E1+1),tau1)
      }
    #______________________________________________________________________

    #predict X(t+1) for all t>1
    Predict[i]<-predictions(E1,LSS,PSS,nearest,tau1,diag1)
    #______________________________________________________________________
  }
  return(Predict)
} #close function cross.val.pred
################################################################################

################################################################################
set.plot.window=function(){
  #start with some graphical book keeping - this just lays out the graphical
  #window in a pattern useful for looking at simplex prediction results...
   plot.new()
   resize.win(Width=20, Height=15, xpinch=30, ypinch=30)
   par(mar=c(4,5,5,2)+0.1,cex.axis=1.5,cex.lab=2,cex.main=2)
   l.out<-matrix(c(1:6),2,3, byrow=F)
   nf <- layout(l.out,respect=F,heights=c(4,4,5,5))
}
################################################################################

################################################################################
#adjust window size
resize.win <- function(Width=10, Height=10, xpinch=20, ypinch=30)
{
        # works for windows
    #dev.new(width=6, height=6)
    dev.off();
    windows(record=TRUE, width=Width, height=Height)
}
################################################################################

################################################################################
plot.series.cross=function(pred.x,obs.x,E1,tau1) {
# This one is just a plot of the predictions vs. the observed time series
#_______________________________________________________________________________

   l.b<-min(c(min(pred.x,na.rm = T),min(obs.x,na.rm = T)))
   u.b<-max(c(max(pred.x,na.rm = T),max(obs.x,na.rm = T)))

   mlab<-paste("E = ",E1," Tau = ",tau1,sep="")
   plot(obs.x~c(1:length(obs.x))
    ,type="b"
    ,xlab="Time"
    ,ylab="Index"
    #,main=paste("Model Fit in",E1,"Dimensions",sep=" ")
    ,main=mlab
    ,xlim=c(0,(length(obs.x)+tau1))
    ,ylim=c((l.b-1),(u.b+1))   #standardized variables should have room...
    ,data=na.omit(data.frame(obs.x))
    )

   p.x<-na.omit(pred.x)

   #where the predixtion vector starts given the lag (tau1) is a bit tricky
   st.point<-(((E1-1)*tau1)+1)
   index1<-c((st.point+1):(length(p.x)+st.point))

      #print(c(length(p.x),length(index1)))

   if(tau1>1) index1[length(index1)]<-index1[length(index1)]+(tau1-1)
   lines(p.x~index1
      ,col="red"
      ,type="b"
      #,data=na.omit(data.frame(pred.x))
    )
   c1<-cor(obs.x,pred.x,use="pairwise.complete.obs")
   c1<-round(c1,4)
   #sub1<-bquote(rho == .(c1))
   rmse<-round(sqrt(mean(na.omit(obs.x-pred.x)^2)),4)
   sub1<-paste("RMSE = ",rmse,sep="")
   mtext(sub1,line= .21) #put the prediction skill above the
} #end function plot.series.cross
################################################################################

################################################################################
plot.cor=function(plot.x,plot.y,E1,tau1) {
# A simple plot of the observed vs. predicted vectors and their correlation
   #test the correlation between predicted and observed
   c1<-cor(plot.y,plot.x,use="pairwise.complete.obs")
   rng<-ifelse(max(abs(na.omit(plot.x)))>max(abs(na.omit(plot.y))),
      max(abs(plot.x)),max(abs(plot.y)))
   rng<-round(rng+1,0)
   xlm<-c(-rng,rng) #three standard deviations in a standardized data set
   ylm<-xlm
   plot(plot.y~plot.x
       ,ylab="Predicted"
       ,xlab="Observed"
       ,main=paste("Simplex Prediction in ",E1," Dimensions",sep="")
       ,ylim=ylm
       ,xlim=xlm
    )
   mtext(bquote(rho == .(c1)), line= .21) #put the prediction skill above the
   #plot
   if(length(na.omit(plot.y))>1) {
     abline(coef(lm(plot.y~plot.x)) #plot the simple linear regression line for
     #observed vs. predicted
        ,lty=2
      )
   }
   abline(0,1) #plot the one to one line

}  #end fuction plot.cor
################################################################################

################################################################################
plot.prediction.graphics=function(PS,neib1,neib2,LS){
#This is a grahical check on the prediction routine
  print("$$$$$$$$$$$$$$$$$$$$$$")
  print((PS))
  print((neib1))
  print((neib2))
  print(dim(LS))
  print("$$$$$$$$$$$$$$$$$$$$$$")
  plot(LS[,1]~LS[,2]
      ,ylab="Library time t"
      ,xlab="Library time t-tau"
    )

  points(PS[,1]~PS[,2]
    ,col="red"
    ,pch=7
  )

  points(neib1[,1]~neib1[,2]
    ,col="blue"
    ,pch=3
  )
}
################################################################################

################################################################################
plot.series=function(pred.x,obs.x,E1,tau1) {
# This one is just a plot of the predictions vs. the observed time series
#_______________________________________________________________________________

   l.b<-min(c(min(pred.x,na.rm = T),min(obs.x,na.rm = T)))
   u.b<-max(c(max(pred.x,na.rm = T),max(obs.x,na.rm = T)))

   mlab<-paste("E = ",E1," Tau = ",tau1,sep="")
   plot(obs.x~c(1:length(obs.x))
    ,type="b"
    ,xlab="Time"
    ,ylab="Index"
    #,main=paste("Model Fit in",E1,"Dimensions",sep=" ")
    ,main=mlab
    ,xlim=c(0,(length(obs.x)+tau1))
    ,ylim=c((l.b-1),(u.b+1))   #standardized variables should have room...
    ,data=na.omit(data.frame(obs.x))
    )

   p.x<-na.omit(pred.x)

   #find where the predixtion vector starts
   st.point<-(length(obs.x)-length(p.x))
   index1<-c((st.point+1):(length(p.x)+st.point))

   if(tau1>1) index1[length(index1)]<-index1[length(index1)]+(tau1-1)
   lines(p.x~index1
      ,col="red"
      ,type="b"
      #,data=na.omit(data.frame(pred.x))
    )
   c1<-cor(obs.x,pred.x,use="pairwise.complete.obs")
   c1<-round(c1,4)
   #sub1<-bquote(rho == .(c1))
   rmse<-round(sqrt(mean(na.omit(obs.x-pred.x)^2)),4)
   sub1<-paste("RMSE = ",rmse,sep="")
   mtext(sub1,line= .21) #put the prediction skill above the
} #end function plot.series.cross
################################################################################