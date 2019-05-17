### SCRIPT OF A FUNCTION TO MAKE RASTERS WITH DS AND NS TASK

#In the function you'll have to define the dataset with the data of each rat during each session (if you used the code I sent you, this should be an object called "alldata"), the dataframe used as an index (csacqidx), the name of the rat and the session. Also the lower and upper limit of your window around cue onset, which I set as -15 and 15 by default.

### This is the function itself:

makeDSandNSraster<- function (dataset, indexframe, rat, day, lowerLimit=-15, upperLimit=15){
  
  #load data and name it
  alldata <- get(load(dataset))
  csacqidx <- get(load(indexframe))
  
  #Select data
  v <- 1:length(csacqidx$subject) 
  i= v[(csacqidx$subject==rat) & (csacqidx$session==day)]

  toPlot <- alldata[[i]]

  
  #Define the windows around each cue
  cue <- c()
  minWin <- c()
  maxWin <- c()
  kindCue <- c()
  DSend <- c()
  reward <- c()
  NSentry <- c()
  
  for(i in 1:length(toPlot$allCues)){
    minWin <- c(minWin, toPlot$allCues[i]+lowerLimit)
    cue <- c(cue, toPlot$allCues[i])
    maxWin <- c(maxWin, toPlot$allCues[i]+upperLimit)
    kindCue <- c(kindCue, toPlot$orderCues[i])
    DSend <- c(DSend, toPlot$CSplusEnd[i])
    reward <- c(reward, toPlot$rewarddelivery[i])
    NSentry <- c(NSentry, toPlot$CSminusresponse[i])
  }
  
  
  #Data.frame with data per trial
  FirstEnt <- c()
  FirstEntRaster <- c()
  for(i in 1:length(cue)){
    if(any(toPlot$receptacleentries>cue[i] & toPlot$receptacleentries<(cue[i]+10))==T){a <-c(toPlot$receptacleentries[toPlot$receptacleentries>cue[i] & toPlot$receptacleentries<(cue[i]+10)])}else{a <- NA}
    FirstEnt <- c(FirstEnt, min(a))
    FirstEntRaster <- c(FirstEntRaster, FirstEnt[i]-cue[i])
  }
  
  trials<- data.frame(
    trial = 1:length(toPlot$allCues),
    cue = cue,
    kindCue = kindCue,
    entry=FirstEntRaster
  )
  
  #Number of Trials
  noTrials <- 1:length(toPlot$allCues)
  
  #Identify to which row (window) each entry and exit belong if any
  inWinEnt <- c()
  inWinExit <- c()
  for(i in 1:length(toPlot$receptacleentries)){
    
    anyWinEnt <- (toPlot$receptacleentries[i]>minWin) == (toPlot$receptacleentries[i]<maxWin)
    whichWinEnt <- noTrials[anyWinEnt] 
    if(length(whichWinEnt)==1){a <- whichWinEnt}
    if(length(whichWinEnt)==0){a <- 0}
    
    anyWinExit <- (toPlot$receptacleexits[i]>minWin) == (toPlot$receptacleexits[i]<maxWin)
    whichWinExit <- noTrials[anyWinExit]
    if(length(whichWinExit)==1){b <- whichWinExit}
    if(length(whichWinExit)==0){b <- 0}
    
    #Case: two consecuetive cues with overlapping windows.
    if(length(c(whichWinEnt,whichWinExit))==3){a=median(c(whichWinEnt,whichWinExit)); b=median(c(whichWinEnt, whichWinExit))} #A. Entry or exit fall in only one window, the other one falls in two. TAke the median because that is the value that is repeated twice and will be in the middle.
    if(length(c(whichWinEnt,whichWinExit))==4){a=min(whichWinEnt); b=min(whichWinExit)} #B. Both the entry and exit fall in both windows. Choose only the smallest one
    
    inWinEnt <- c(inWinEnt,a)
    inWinExit <- c(inWinExit,b)
    
  }
  
  
  #Fix cases in which the last exit occurs during the last cue and is not recorded. Make the last exit be the same timestamp as the last entry.
  entry=toPlot$receptacleentries
  exit=toPlot$receptacleexits
  if(length(entry)>length(exit)){exit <- c(exit, max(entry))}
  
  
  #Make a (long, one element per entry) vector with the info about cue, minWin, maxWin and kindofCue per entry. Also, the relative time of entries and exit with respect to the cue in their window
  minWinLong <- c()
  cueLong <- c()
  maxWinLong <- c()
  kindOfCueLong <- c()
  entryRaster <- c()
  exitRaster <- c()
  DSendLong <- c()
  rewardLong <- c()
  NSentryLong <- c()
   
  for(i in 1:length(inWinEnt)){
    index<- inWinEnt[i]
    
    minWinSel <- minWin[index]
    if(length(minWinSel)==0) {minWinSel <- NA}
    minWinLong <- c(minWinLong, minWinSel)
    
    cueSel <- cue[index]
    if(length(cueSel)==0) {cueSel <- NA}
    cueLong <- c(cueLong, cueSel)
    
    maxWinSel <- maxWin[index]
    if(length(maxWinSel)==0) {maxWinSel <- NA}
    maxWinLong <- c(maxWinLong, maxWinSel)
    
    kindOfCueSel <- kindCue[index]
    if(length(kindOfCueSel)==0) {kindOfCueSel <- NA}
    kindOfCueLong <- c(kindOfCueLong, kindOfCueSel)
    
    kindOfCueSel <- kindCue[index]
    if(length(kindOfCueSel)==0) {kindOfCueSel <- NA}
    kindOfCueLong <- c(kindOfCueLong, kindOfCueSel)
    
    
    #Calculate ent and exit relative to cue. 
    a <- entry[i]-cueLong[i]
    z <- exit[i]-cueLong[i]
    
    #Fix cases where one of them falls outside the window. Two cases: when 
    if(inWinEnt[i]<inWinExit[i] && inWinEnt[i]==0){a <- lowerLimit}
    if(inWinEnt[i]>inWinExit[i] && inWinExit[i]==0){z<- upperLimit}
    
    entryRaster <- c(entryRaster, a)
    exitRaster <- c(exitRaster, z)
    
  }
  
  #Data frame for plotting the beh rasters
  rasterData <- data.frame (entry=entry, exit=exit, minWin=minWinLong, cue=cueLong, maxWin=maxWinLong, inWinEnt=inWinEnt, inWinExit=inWinExit, kindCue=kindOfCueLong, entryRaster=entryRaster, exitRaster=exitRaster)
  
  
  
  
  #Define plot dimensions
  plot.new()
  plot.window (xlim=c(lowerLimit, upperLimit), ylim=c(length(toPlot$orderCues), 1))
  
  
  ###DS or NS segments
  
  #10s after cue
  cueLength=10
  abline(v=cueLength, col="darkred")
  
  #Gray bars behind DS trials
  for(i in 1:length(toPlot$orderCues)){
    if(toPlot$orderCues[i]==1){
      segments(lowerLimit, i, upperLimit, i, col="gray90", lwd=5)
    }
  }
  
  #Segments for entries
  for(i in 1:length(rasterData$entry)){
    segments(entryRaster[i], inWinEnt[i], exitRaster[i], inWinEnt[i])
    }
  
  #Points for entries
  for(i in 1:length(trials$trial)){
    if(trials$kindCue[i]==1){points(trials$entry[i], trials$trial[i], col="blue", pch=16)} #rewarded trials
    if(trials$kindCue[i]==1 && is.na(trials$entry[i])){points(cueLength, trials$trial[i], col="orange", pch=16)} #Missed DS trials
  }
  
  for(i in 1:length(trials$trial)){
    if(trials$kindCue[i]==2){points(trials$entry[i],trials$trial[i], col="darkorange3", pch=1)}
  } #Responded to NS trials
  
  
  #Cue
  abline(v=0, col="red", lwd=2)
  
  
  
  #Define axes
  axis(1, at=seq(lowerLimit, upperLimit, by=5), labels=seq(lowerLimit, upperLimit, by=5), cex.axis = 1, las=1)
  
  axis(2, at=c(1, seq(5, length(toPlot$orderCues), 5)), labels=c(1, seq(5, length(toPlot$orderCues), 5)) , cex.axis = 1, las = 2, pos=-15)
  
  mtext("Time from cue onset (s)", side=1, line=2.5, font=3, col="gray50")
  mtext("Trial", side=2, line=2.5, font=3, col="gray50")
  
  title(paste(rat," " , "Session:", day))

}



## Once you run this function, you can just run this line changing the parameters and the raster will appear in your "Plots" window, for example:

makeDSandNSraster(dataset="C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 3/Beh pilots/Badge 11a/Training/Data for R/alldata.rdat", indexframe="C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/Experiment 3/Beh pilots/Badge 11a/Training/Data for R/csacqidx.ridx", rat="MV153", day="1", lowerLimit=-15, upperLimit=15)

