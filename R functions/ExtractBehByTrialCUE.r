ExtractBehByTrialCUE <- function (object=allNeuronsDS, cueLength=10, startt=startt, binw=binw, psthmin=psthmin, psthmax=psthmax, event=c("DS_onset", "NS_onset"), side="both", first="left"){
        
        
        ###############################
        # Select all trials and units #
        ###############################
        
        #Build list of all summary (nexdata) files (one per NEX file)
        nexdata <- lapply(seq(1, length(allNeuronsDS$nexdata)), function(x){allNeuronsDS$nexdata[[x]]})
        
        #Timestamps of all behavioral events
        DFbeh <- do.call("rbind", lapply(seq(1, length(nexdata)), function(x){
                events <- event #Which events are called for in the function
                eventsIdx <- match(events, names(nexdata[[x]])) #Index of those events in nexdata
                allEvents <- sort(unlist(nexdata[[x]][eventsIdx])) #Data array w all events sorted
                
                #Build data frame w timestamps of the cues and session info
                eventsDF <- data.frame(rat=nexdata[[x]]$ratname, 
                                       expt=nexdata[[x]]$expt, 
                                       evName=gsub('[[:digit:]]+', '', names(allEvents)),#Name of the event, 
                                       evDetail=names(allEvents),
                                       evTime=as.numeric(allEvents),
                                       trial=1:length(allEvents))
                
                #If my data frame is built around CUE ONSET data, I'll want to know latency, etc. of entries and exits
                #if(sum(grepl("onset", event))>=1){
                events2 <- c("Entry") #Identify the second event of interest, in this case "ENTRY"
                events2Idx <- match(events2, names(nexdata[[x]])) 
                entries <- as.numeric(sort(unlist(nexdata[[x]][events2Idx]))) #All entries
                entriesTrial <- findInterval(entries, eventsDF$evTime) #Which entries happened in which trial
                events3 <- c("Exit") #Identify third event of interest, in this case "EXIT"
                events3Idx <- match(events3, names(nexdata[[x]])) 
                exits <- as.numeric(sort(unlist(nexdata[[x]][events3Idx]))) #All exits
                exitsTrial <- findInterval(exits, eventsDF$evTime) #Which exits happened in which trial
                #For each trial in this session
                entryData <- c()
                entryData <- do.call("rbind", lapply(seq(1, length(eventsDF$trial)), function(y){
                        trlIdx <- eventsDF$trial[y]  #Identify the trial
                        masterEntry <- entries[entriesTrial==trlIdx] #Entries on that trial
                        FirstEnt <- masterEntry[1]; SecondEnt <- masterEntry[2]  #First and second entry on that trial
                        if(is.na(FirstEnt)){FirstEnt <- 0}; if(is.na(SecondEnt)){SecondEnt <- 0} #If no entry was made, make it 0
                        resp <- NA #Did the animal respond before the end of the cue?
                        if(FirstEnt>eventsDF$evTime[y] && FirstEnt<=(eventsDF$evTime[y]+cueLength)){resp=TRUE} else {resp=FALSE}
                        latency <- if(resp==T){FirstEnt-eventsDF$evTime[y]} else {latency=NA} #Latency of the cued response if any
                        masterExit <- exits[exitsTrial==trlIdx] #Exits on that trial
                        FirstExt <- masterExit[1]; SecondExt <- masterExit[2] #First and second exits on that trial
                        if(is.na(FirstExt)){FirstExt <- 0}; if(is.na(SecondExt)){SecondExt <- 0} #Assign 0 if no exit was made on that trial
                        FirstEntDur <- NA; SecondEntDur <- NA; IEI1 <- NA #1st and 2nd entry duration and inter entry interval (btwn 1st and 2nd entry)
                        if(FirstExt>FirstEnt){FirstEntDur <- FirstExt-FirstEnt} else {FirstEntDur <- NA}
                        if(SecondEnt>FirstExt){IEI1 <- SecondEnt-FirstExt} else {IEI1 <- NA}
                        if(SecondExt>SecondEnt){SecondEntDur <- SecondExt-SecondEnt} else {SecondEntDur <- NA}
                        
                        dat <- c(trlIdx, resp, latency, FirstEnt, FirstExt, SecondEnt, SecondExt, FirstEntDur, IEI1, SecondEntDur)
                        
                }))
                
                nms <- c("trial", "resp", "latency", "FirstEnt", "FirstExt", "SecondEnt", "SecondExt", "FirtEntDur", "IEI1", "SecondEntDur")
                entryData <- as.data.frame(entryData)
                names(entryData) <- nms
                
                # Calculate session-wise data by kind of trial (DSRR, NSRR, etc.) and add it to eventsDF
                DSsubset <- subset(entryData, grepl("DS", eventsDF$evName)) #Select DS trials in entryData
                DSRR <- sum(DSsubset$resp)/length(DSsubset$resp)
                DSlat <- mean(DSsubset$latency, na.rm=T)
                DSlatSE <- sd(DSsubset$latency, na.rm=T)/sqrt(sum(!is.na(DSsubset$latency)))
                
                NSsubset <- subset(entryData, grepl("NS", eventsDF$evName)) #Select NS trials in entryData
                NSRR <- sum(NSsubset$resp)/length(NSsubset$resp)
                NSlat <- mean(NSsubset$latency, na.rm=T)
                NSlatSE <- sd(NSsubset$latency, na.rm=T)/sqrt(sum(!is.na(NSsubset$latency)))
                
                TaskAccProb <- DSRR-NSRR #Prob of discrimination DS vs NS cues
                TaskAccLat <- NSlat-DSlat #How much faster were they to respond to DS vs. NS
                
                eventsDF$DSRR <- DSRR
                eventsDF$NSRR <- NSRR
                eventsDF$DSlatSess <- DSlat
                eventsDF$NSlatSess <- NSlat
                eventsDF$DSlatSE <- DSlatSE
                eventsDF$NSlatSE <- NSlatSE
                eventsDF$TaskAccProb <- TaskAccProb
                eventsDF$TaskAccLat <- TaskAccLat
                
                
                merge(eventsDF, entryData, by="trial")
                
                #Timestamps of all units
                #DFneu <- c()
        }))
}