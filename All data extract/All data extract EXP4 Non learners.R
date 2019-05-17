
#############################################################
### EXPERIMENT 4: UNILATERAL INFUSIONS, NON-LEARNERS      ###
#############################################################

### LOAD IMPORTANT LIBRARIES
install.packages("matrixStats")
install.packages('ez')
install.packages('dplyr')

library(matrixStats)
library(ez)
library(dplyr)

setwd("E:/Dropbox/NMDA/")

#Define colors
colindx <- c("#2171b5", "#cb181d") #Strong blue and red
colindxB <- c("#bdd7e7", "#fcae91") #Less strong blue and red
colindxC <- c("#eff3ff", "#fb6a4a") #Even less strong blue and red
colindxD <- c("#6baed6", "#fee5d9") #Lightest blue and red


### DEFINE ALL THE IMPORTANT FOLDERS

funcdirect <- paste(getwd(), "/R functions/", sep="")
datafolder <- paste(getwd(), "/EXP4_Unilateral AP5 Non learners/MedPC files/", sep="")
dataForRdir <- paste(getwd(), "/EXP4_Unilateral AP5 Non learners/Data for R/", sep="")
dataForRCumulative <- paste(getwd(), "/EXP4_Unilateral AP5 Non learners/Data for R cumulative/", sep="")
behGraphFolder <- paste(getwd(), "/EXP4_Unilateral AP5 Non learners/Graphs/Behavior/", sep="")
CPfuncFolder <- paste(getwd(), '/R functions/Change_Point-master/', sep="")
CPGraphFolder <- paste(getwd(), "/EXP4_Unilateral AP5 Non learners/Graphs/Behavior/Change point/", sep="")
MixedGraphFolder <- paste(getwd(), "/EXP4_Unilateral AP5 Non learners/Graphs/Mixed/", sep="")
NEXfiles <- paste(getwd(), "/EXP4_Unilateral AP5 Non learners/NEX files/", sep="")
IndcumPerfGraphFolder <- paste(getwd(), "/EXP4_Unilateral AP5 Non learners/Graphs/Behavior/Cumulative perf/", sep="")
PerfRelToCPFolder <- paste(getwd(), "/EXP4_Unilateral AP5 Non learners/Graphs/Behavior/Beh rel to CP/", sep="")
BySessFolder <- paste(MixedGraphFolder, "By Sess from CP/", sep="")
ScatterplotFolder <- paste(MixedGraphFolder, "Scatterplot/", sep="")


### Load necessary functions
load(file=paste(funcdirect, "MedPCextract.r", sep=""))
load(file=paste(funcdirect, "mpcextract_blockSingle.Rfunc", sep=""))
load(file=paste(funcdirect, "CPextract.R", sep=""))
load(file=paste(funcdirect, "cumulativeIndGraphs.R", sep=""))
load(file=paste(funcdirect, "PerformanceFromCP.R", sep=""))
load(file=paste(funcdirect, "PrePostCP_Perf.R", sep=""))
load(file=paste(funcdirect, "avgPerfByBin.R", sep=""))
load(file=paste(funcdirect, "CPextractMultipleCrit.R", sep=""))
load(file=paste(funcdirect, "neuralhist.r", sep=""))
load(file=paste(funcdirect, "FRbyNEURONbyBINcue.r", sep=""))
load(file=paste(funcdirect, "plotFRandCP.r", sep=""))
load(file=paste(funcdirect, "errBars.r", sep=""))
load(file=paste(funcdirect, "errCloud.r", sep=""))
load(file=paste(funcdirect, "masterDFsummary.r", sep=""))
load(file=paste(funcdirect, "cp_wrapper3.r", sep=""))
load(file=paste(funcdirect, "PLOT_perRatTrial_FRandBEH.r", sep=""))
load(file=paste(funcdirect, "RealCP.r", sep=""))
load(file=paste(funcdirect, "plotFRandCPhistogram.R", sep=""))
load(file=paste(funcdirect, "megaplot.r", sep=""))
load(file=paste(funcdirect, "sessFromCPdata.r", sep=""))
load(file=paste(funcdirect, "giveStars.R", sep=""))


#### GENERATE IMPORTANT DATA OBJECTS
MedPCextract(MovAvg="Impinged only", funcdirect = funcdirect, datafolder = datafolder, dataForRdir = dataForRdir, dataForRCumulative=dataForRCumulative, cuelength=10)

# Load important behavior-related objects
files <- paste(dataForRdir, list.files(dataForRdir), sep="")
filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep="")
for(i in 1:length(files)){load(files[[i]])}
for(i in 1:length(filesCum)){load(filesCum[[i]])}


#### BEHAVIORAL GRAPHS

# RASTERS
load(file=paste(funcdirect, "MAKERASTER.r", sep=""))
MAKERASTER(i=x, data=alldata, idxdata=csacqidx)
MAKERASTER(i=4, data=alldata, idxdata=csacqidx)


###
# I'M USING THE 10S PRE-CUE PERIOD AS BASELINE
# In 9% of the trials, the animal's head was inside the receptacle at the time the ITI window started. In those trials, I'll correct the ITI latency by assigning the MEAN ITI latency in the +/- 5 trials around that trial
###

#AVERAGE PERFORMANCE S+ VS S-. 5 trial bins.
avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(DStaskAcc, NStaskAcc), 
             cues=c("S+", "S-"), index="Specificity", legendLocation="topleft", y_axis_label="Performance index", 
             behGraphFolder=behGraphFolder, plot=T)
        
avgPerfByBin(binsize=5, colors=c(colindx[1], "darkblue", "gray40"), data=list(DSrespAll, NSrespAll, ITIrespRatio), 
             cues=c("S+", "S-"), index="Response ratio", legendLocation="topleft", y_axis_label="Response ratio", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(DSlatency, NSlatency), 
             cues=c("S+", "S-"), index="Latency", legendLocation="bottomright", y_axis_label="Latency (s)", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(ITIlatency), 
             cues=c("ITI"), index="Latency", legendLocation=-1, y_axis_label="Latency (s)", 
             behGraphFolder=behGraphFolder, plot=T)



### CHANGE POINT: S+ SPECIFICITY
# Change point analysis using the method in Gallistel et al., 2004. I'll use S+ SPECIFICITY as the main performance index. But I also need to look at the results that other performance indexes give me.

CPdata <- CPextract(GallCrit=1.3, plot=T, CPfuncFolder=CPfuncFolder, idx=idx, CPGraphFolder=CPGraphFolder, dataForRdir=dataForRdir, dataForRCumulative=dataForRCumulative)



#CHANGE POINT: Other performance indexes (response ratio, latency, ITI latency)
#This next line takes a long time. It's a calculation of change point using the other indexes (not the S+ specificity index). The method "intersect" calculates CPs using different criteria (from least to most conervative) and finds the earliest point that intersects
CPdataAllCrit <- CPextractMultipleCrit(GallCrit=c(1.3, 2, 4, 5, 6), idx, CPfuncFolder, CPGraphFolder, dataForRdir, dataForRCumulative)
RealCPbyIndex <- RealCP(data=CPdataAllCrit, method="intersect", 
       index_names=c("DSlatency", "DSrespAll", "DStaskAcc", "DStimeToSpare", "ITIlatency", 
                     "NSlatency", "NSrespAll", "NStaskAcc", "NStimeToSpare"), rats=rats)
        
#Find the CPs of the indexes I'm interested in:
targetIndexes <- c("DSrespAll", "DSlatency", "ITIlatency")
sel <- match(targetIndexes, rownames(RealCPbyIndex))
DSspecificity <- CPdata$CP #This one is the only one that can adopt negative values. So I had to calculate it differently.
CPframeExp4 <- rbind(DSspecificity, RealCPbyIndex[sel,])


#CUMULATIVE INDIVIDUAL PERFORMANCE
#All rats
cumulativeIndGraphs(numSess=6, numCues=35, sessLines=F, smartRat = "NA", dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=IndcumPerfGraphFolder, imgFormat="pdf", limit=210)



#####################################3333
## LET'S ASSIGN 1 AS THE CP and NS.PreCP TRIAL (SINCE THERE IS NO CP FOR THESE ANIMALS, WE'LL ALIGN EVERYTHING TO THE FIRST TRIAL)
CPdata$CP <- rep(1, nrow(CPdata))
CPdata$CPsess <- rep(1, nrow(CPdata))

CPdata$NS.PreCP <- rep(1, nrow(CPdata))

#MEAN PERFORMANCE AS A FUNCTION OF DISTANCE FROM CP
PerformanceFromCP(relTrialMin=0, relTrialMax=210, trialBinSize=5, typegraph="both", imgFormat="pdf", 
                  graphFolder=PerfRelToCPFolder)


#THIS WILL BE IMPORTANT LATER: SESSION IN WHICH THE CP TOOK PLACE
#I'll use this later to figure out in which session the CP took place
cueKind=1 #cueKind=1 is S+, cueKind=2 is S-
ratsTrialIdxPerSess <- lapply(seq(1, length(idx)), function(l){
        trlPerSess <- sapply(seq(1, length(idx[[l]])), function(m){
                if(cueKind==1){trlPerSess <- length(alldata[[idx[[l]][m]]]$CSpluscue)}
                if(cueKind==2){trlPerSess <- length(alldata[[idx[[l]][m]]]$CSminuscue)}
                return(trlPerSess)
        })
        cumsum(trlPerSess)
})

#Session in which the CP took place per rat (each item is a rat, the order is given by the object "rats")
CPvector <- CPdata$CP
sessionCPperRat <- sapply(seq(1, length(CPvector)), function(l){findInterval(CPvector[l], ratsTrialIdxPerSess[[l]])+1})


################################################################################
############ NEURONAL DATA #####################################################
################################################################################

#### RUN THIS FOR NEURONS IN VEHICLE SIDE
funcdirect <- paste(getwd(), "/R functions/", sep="")
allNeuronsDS_NL <-neuralhist(funcdirect=funcdirect, path=NEXfiles, event=1, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="both")
allNeuronsDS_VEH_NL = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=1, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsNS_VEH_NL = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=2, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsDSresponded_VEH_NL = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=5, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsDSmissed_VEH_NL = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=6, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsNSresponded_VEH_NL = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=7, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsNSmissed_VEH_NL = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=8, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsEntryDS_VEH_NL = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsEntryNS_VEH_NL = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsEntryITI_VEH_NL = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")


#### RUN THIS FOR NEURONS IN AP5 SIDE
allNeuronsDS_AP5_NL = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=1, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsNS_AP5_NL = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=2, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsDSresponded_AP5_NL = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=5, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsDSmissed_AP5_NL = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=6, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsNSresponded_AP5_NL = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=7, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsNSmissed_AP5_NL = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=8, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsEntryDS_AP5_NL = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsEntryNS_AP5_NL = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsEntryITI_AP5_NL = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")


#sAVE THESE OBJECTS
save(allNeuronsDS_VEH_NL, file=paste(dataForRdir, 'allNeuronsDS_VEH_NL.rdat', sep=""))
save(allNeuronsNS_VEH_NL, file=paste(dataForRdir, 'allNeuronsNS_VEH_NL.rdat', sep=""))
save(allNeuronsDSresponded_VEH_NL, file=paste(dataForRdir, 'allNeuronsDSresponded_VEH_NL.rdat', sep=""))
save(allNeuronsDSmissed_VEH_NL, file=paste(dataForRdir, 'allNeuronsDSmissed_VEH_NL.rdat', sep=""))
save(allNeuronsNSresponded_VEH_NL, file=paste(dataForRdir, 'allNeuronsNSresponded_VEH_NL.rdat', sep=""))
save(allNeuronsNSmissed_VEH_NL, file=paste(dataForRdir, 'allNeuronsNSmissed_VEH_NL.rdat', sep=""))
save(allNeuronsEntryDS_VEH_NL, file=paste(dataForRdir, 'allNeuronsEntryDS_VEH_NL.rdat', sep=""))
save(allNeuronsEntryNS_VEH_NL, file=paste(dataForRdir, 'allNeuronsEntryNS_VEH_NL.rdat', sep=""))
save(allNeuronsEntryITI_VEH_NL, file=paste(dataForRdir, 'allNeuronsEntryITI_VEH_NL.rdat', sep=""))
save(allNeuronsDS_AP5_NL, file=paste(dataForRdir, 'allNeuronsDS_AP5_NL.rdat', sep=""))
save(allNeuronsNS_AP5_NL, file=paste(dataForRdir, 'allNeuronsNS_AP5_NL.rdat', sep=""))
save(allNeuronsDSresponded_AP5_NL, file=paste(dataForRdir, 'allNeuronsDSresponded_AP5_NL.rdat', sep=""))
save(allNeuronsDSmissed_AP5_NL, file=paste(dataForRdir, 'allNeuronsDSmissed_AP5_NL.rdat', sep=""))
save(allNeuronsNSresponded_AP5_NL, file=paste(dataForRdir, 'allNeuronsNSresponded_AP5_NL.rdat', sep=""))
save(allNeuronsNSmissed_AP5_NL, file=paste(dataForRdir, 'allNeuronsNSmissed_AP5_NL.rdat', sep=""))
save(allNeuronsEntryDS_AP5_NL, file=paste(dataForRdir, 'allNeuronsEntryDS_AP5_NL.rdat', sep=""))
save(allNeuronsEntryNS_AP5_NL, file=paste(dataForRdir, 'allNeuronsEntryNS_AP5_NL.rdat', sep=""))
save(allNeuronsEntryITI_AP5_NL, file=paste(dataForRdir, 'allNeuronsEntryITI_AP5_NL.rdat', sep=""))

#Load files again      
#files <- paste(dataForRdir, list.files(dataForRdir), sep=""); for(i in 1:length(files)){load(files[[i]])}

#BUILD A DATAFRAME WITH FIRING RATE OF EACH NEURON ON EACH TRIAL WITH BEH INFO ABOUT EACH TRIAL
     
#VEH side
masterDF_DS_VEH_NL <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS_VEH_NL, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, 
                                      BLneudata=allNeuronsDS_VEH_NL)

#By modality
masterDF_DS_VEH_Tone_NL <- filter(masterDF_DS_VEH_NL, masterDF_DS_VEH_NL$modality==98)
masterDF_DS_VEH_Light_NL <- filter(masterDF_DS_VEH_NL, masterDF_DS_VEH_NL$modality==99)



masterDF_DSresponded_VEH_NL <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSresponded_VEH_NL, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDSresponded_VEH_NL)
masterDF_DSmissed_VEH_NL <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSmissed_VEH_NL, 
                                               CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                               BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDSmissed_VEH_NL)
masterDF_NS_VEH_NL <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_VEH_NL, 
                                            CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                            BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_VEH_NL)
masterDF_NSresponded_VEH_NL <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_VEH_NL, 
                                      CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_VEH_NL)
masterDF_NSmissed_VEH_NL <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_VEH_NL, 
                                      CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_VEH_NL)


#The "entry-related" masterDF don't give me FR info of every trial, only of those trials in which the animal responded
masterDF_DSEntry_VEH_NL <- FRbyNEURONbyBINcue(eventType="entry", cue="S+", neudata=allNeuronsEntryDS_VEH_NL, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDS_VEH_NL)
masterDF_NSEntry_VEH_NL <- FRbyNEURONbyBINcue(eventType="entry", cue="S+", neudata=allNeuronsEntryNS_VEH_NL, 
                                           CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                           BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_VEH_NL)

#AP5 side
masterDF_DS_AP5_NL <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS_AP5_NL, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDS_AP5_NL)


#By modality
masterDF_DS_AP5_Tone_NL <- filter(masterDF_DS_AP5_NL, masterDF_DS_AP5_NL$modality==98)
masterDF_DS_AP5_Light_NL <- filter(masterDF_DS_AP5_NL, masterDF_DS_AP5_NL$modality==99)


masterDF_DSresponded_AP5_NL <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSresponded_AP5_NL, 
                                               CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                               BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDSresponded_AP5_NL)
masterDF_DSmissed_AP5_NL <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSmissed_AP5_NL, 
                                            CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                            BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDSmissed_AP5_NL)
masterDF_NS_AP5_NL <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_AP5_NL, 
                                      CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_AP5_NL)
masterDF_NSresponded_AP5_NL <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_AP5_NL, 
                                               CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                               BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_AP5_NL)
masterDF_NSmissed_AP5_NL <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_AP5_NL, 
                                            CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                            BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_AP5_NL)

masterDF_DSEntry_AP5_NL <- FRbyNEURONbyBINcue(eventType="entry", cue="S+", neudata=allNeuronsEntryDS_AP5_NL, 
                                           CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                           BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDS_AP5_NL)
masterDF_NSEntry_AP5_NL <- FRbyNEURONbyBINcue(eventType="entry", cue="S+", neudata=allNeuronsEntryNS_AP5_NL, 
                                           CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                           BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_AP5_NL)



save(masterDF_DS_VEH_NL, file=paste(dataForRdir, "masterDF_DS_VEH_NL.rdat", sep=""))
save(masterDF_DSresponded_VEH_NL, file=paste(dataForRdir, "masterDF_DSresponded_VEH_NL.rdat", sep=""))
save(masterDF_DSmissed_VEH_NL, file=paste(dataForRdir, "masterDF_DSmissed_VEH_NL.rdat", sep=""))
save(masterDF_NS_VEH_NL, file=paste(dataForRdir, "masterDF_NS_VEH_NL.rdat", sep=""))
save(masterDF_NSresponded_VEH_NL, file=paste(dataForRdir, "masterDF_NSresponded_VEH_NL.rdat", sep=""))
save(masterDF_NSmissed_VEH_NL, file=paste(dataForRdir, "masterDF_NSmissed_VEH_NL.rdat", sep=""))
save(masterDF_DSEntry_VEH_NL, file=paste(dataForRdir, "masterDF_DSEntry_VEH_NL.rdat", sep=""))
save(masterDF_NSEntry_VEH_NL, file=paste(dataForRdir, "masterDF_NSEntry_VEH_NL.rdat", sep=""))

save(masterDF_DS_AP5_NL, file=paste(dataForRdir, "masterDF_DS_AP5_NL.rdat", sep=""))
save(masterDF_DSresponded_AP5_NL, file=paste(dataForRdir, "masterDF_DSresponded_AP5_NL.rdat", sep=""))
save(masterDF_DSmissed_AP5_NL, file=paste(dataForRdir, "masterDF_DSmissed_AP5_NL.rdat", sep=""))
save(masterDF_NS_AP5_NL, file=paste(dataForRdir, "masterDF_NS_AP5_NL.rdat", sep=""))
save(masterDF_NSresponded_AP5_NL, file=paste(dataForRdir, "masterDF_NSresponded_AP5_NL.rdat", sep=""))
save(masterDF_NSmissed_AP5_NL, file=paste(dataForRdir, "masterDF_NSmissed_AP5_NL.rdat", sep=""))
save(masterDF_DSEntry_AP5_NL, file=paste(dataForRdir, "masterDF_DSEntry_AP5_NL.rdat", sep=""))
save(masterDF_NSEntry_AP5_NL, file=paste(dataForRdir, "masterDF_NSEntry_AP5_NL.rdat", sep=""))

#masterDFsumm <- masterDFsummary(masterDF=masterDF_DS, neudata=allNeuronsDS, format="Zscores", postCueWdw=c(0, 0.4), postSpikeWdw=c(0.4, 5), orEntry=T, cuelength=10)


### PLOT FR POST-CUE (400MS WINDOW) as a function of distance to change point

#Redefine colors
colindx <- c("#2171b5", "#cb181d") #Strong blue and red

#CUE
plotFRandCP(experiment="Exp 4 Non Learners VEH side", cue=c("S+", "S-"), masterDF=list(masterDF_DS_VEH_NL, masterDF_NS_VEH_NL), 
            graphFolder=MixedGraphFolder, 
            trialBinSize=5, WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
            colindx =c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), capped=T, capValue = c(0, 210), 
            yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_VEH_NL)

plotFRandCP(experiment="Exp 4 Non Learners AP5 side", cue=c("S+", "S-"), masterDF=list(masterDF_DS_AP5_NL, masterDF_NS_AP5_NL), 
            graphFolder=MixedGraphFolder, 
            trialBinSize=5, WdwStart=100, WdwEnd=400, dataProcess="Zscores", legLabels=c("S+", "S-"),
            correctOnly=FALSE, colindx=c(colindx[2], "darkred"), capped=T, capValue = c(0, 210), 
            yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_AP5_NL)

plotFRandCP(experiment="Exp 4 Non Learners BOTH SIDES", cue=c("S+", "S-", "S+", "S-"), 
            masterDF=list(masterDF_DS_VEH_NL, masterDF_NS_VEH_NL, masterDF_DS_AP5_NL, masterDF_NS_AP5_NL), 
            graphFolder=MixedGraphFolder, trialBinSize=5, WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
            colindx =c("#2171b5", "darkblue", "#cb181d", "darkred"), legLabels=c("S+ VEH_NL", "S- VEH_NL", "S+ AP5_NL", "S- AP5_NL"), 
            capped=T, capValue = c(0, 210), cueExcOnly = F,
            yAxMinZ = -1, yAxMaxZ = 10, yAxMaxRaw = 10, neudata=allNeuronsDS_VEH_NL)



#% Exc
PercCueExc_VEH <-  plotFRandCP(experiment="Exp 4 Non learners VEH side", cue=c("S+"), masterDF=list(masterDF_DS_VEH_NL), 
                               graphFolder=MixedGraphFolder, 
            trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, 
            colindx =c(colindx[1]), legLabels=c("S+"), capped=T, capValue = c(0, 210), 
            yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_VEH_NL)

save(PercCueExc_VEH, file=paste(dataForRdir, "PercCueExc_VEH.rdat", sep=""))

# [[1]]
#   bin trialBins CueEx notCueExc CueInh notCueInh Total N
# 1   1         0     4        52     15        41      56
# 2   2        35     6        80     25        61      86
# 3   3        70     6        84     23        67      90
# 4   4       105     2        44     13        33      46
# 5   5       140     2        29      7        24      31
# 6   6       175     9        25      3        31      34


PercCueExc_AP5 <- plotFRandCP(experiment="Exp 4 Non learners AP5 side", cue=c("S+"), masterDF=list(masterDF_DS_AP5_NL), 
                              graphFolder=MixedGraphFolder, 
            trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, 
            colindx =c(colindx[2]), legLabels=c("S+"), capped=T, capValue = c(0, 210), 
            yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_AP5_NL)

save(PercCueExc_AP5, file=paste(dataForRdir, "PercCueExc_AP5.rdat", sep=""))

# [[1]]
#   bin trialBins CueEx notCueExc CueInh notCueInh Total N
# 1   1         0     3        38      9        32      41
# 2   2        35     2        46      6        42      48
# 3   3        70     2        58      5        55      60
# 4   4       105     3        48     10        41      51
# 5   5       140     3        45     12        36      48
# 6   6       175     5        32      8        29      37



#Is the proportion of cue-exc neurons modulated by a) the drug; b) the amount of training? (run a different chi.sq test for bins 1, 3 and 5 and for bins 2, 4 and 6, that way the units on each cell will be independent for sure)
#For Chi-sq analysis. I need to make sure that each cell has a completely different population of neurons. Bc each session has 40 trials of each kind, by only comparing bins that are separated by 40 trials I can guarantee that.


### VEH SIDE
contingency_table_comp1<- PercCueExc_VEH[[1]][is.even(PercCueExc_VEH[[1]]$bin)==FALSE, c(2, 3, 4)]
contingency_table_comp2 <-  PercCueExc_VEH[[1]][is.even(PercCueExc_VEH[[1]]$bin)==TRUE, c(2, 3, 4)]

inhCont_table_comp1 <- PercCueExc_VEH[[1]][is.even(PercCueExc_VEH[[1]]$bin)==FALSE, c(2, 5, 6)]
inhCont_table_comp2 <- PercCueExc_VEH[[1]][is.even(PercCueExc_VEH[[1]]$bin)==TRUE, c(2, 5, 6)]

#Chi-sq analysis excitations:
#Bins 1, 3 and 5
critData <- contingency_table_comp1[, -1]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 0.018785, df = 2, p-value = 0.9907

fisher.test(critData) #p=1

#Bins 2, 4 and 6
critData <- contingency_table_comp2[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 12.477, df = 2, p-value = 0.001952
# 

fisher.test(critData) # p=0.004394

#Post-hoc comparisons for % cue-excited units
# Is the distribution of exc/not exc neurons different from the first (-120 to -80) and the third (-40 to 0) bin. 
contTableExc_VEH <- PercCueExc_VEH[[1]][ , c(3, 4)]

#Bin 1 vs. 3
fisher.test(contTableExc_VEH[c(1, 3), ], alternative="less") #95% CI: 0.000000 3.884636 Odds ratio=1.076372 p=0.6785
#Bin 3 vs. 5
fisher.test(contTableExc_VEH[c(3, 5), ], alternative="less") #95% CI:0.0000000 7.453629 odds ratio=1.035443 p= 0.6564

#Bin 2 vs. 4
fisher.test((contTableExc_VEH[c(2, 4), ]), alternative="less") #95% CI: 0.0000000 11.70066 odds ratio= 1.644225  p= 0.837
#Bin 4 vs. 6
fisher.test((contTableExc_VEH[c(4, 6), ]), alternative="less") #95% CI:  0.000000 0.5664629 odds ratio=0.1295608 p= 0.005784

#Correct p value all comparisons
p.adjust(p=c(0.6785, 0.6564, 0.837, 0.005784), method="holm")
# 1.000000 1.000000 1.000000 0.023136

#Chi-sq for inhibitions:
#Bins 1, 3 and 5
critData <- inhCont_table_comp1[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 0.18777, df = 2, p-value = 0.9463

fisher.test(critData) # p=0.00439

#Bins 2, 4 and 6
critData <- inhCont_table_comp2[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 5.8045, df = 2, p-value = 0.0549

fisher.test(critData) # p=0.04434



contTableInh_VEH <- PercCueExc_VEH[[1]][ , c(5, 6)]

#Bin 1 vs. 3
fisher.test(contTableInh[c(1, 3), ], alternative="less") #95% CI: 0.000000 2.144489 Odds ratio=1.065281 p=0.6423
#Bin 3 vs. 5
fisher.test(contTableInh[c(3, 5), ], alternative="less") #95% CI:0.0000000  3.076922 odds ratio=1.175382 p= 0.7112

#Bin 2 vs. 4
fisher.test(contTableInh[c(2, 4), ], alternative="less") #95% CI: 0.000000 1.163983 Odds ratio=0.4780756 p=0.09493
#Bin 4 vs. 6
fisher.test(contTableInh[c(4, 6), ], alternative="less") #95% CI:0.00000 42.70053  odds ratio=6.688493 p= 0.9995




### AP5 SIDE
contingency_table_comp1<- PercCueExc_AP5[[1]][is.even(PercCueExc_AP5[[1]]$bin)==FALSE, c(2, 3, 4)]
contingency_table_comp2 <-  PercCueExc_AP5[[1]][is.even(PercCueExc_AP5[[1]]$bin)==TRUE, c(2, 3, 4)]

inhCont_table_comp1 <- PercCueExc_AP5[[1]][is.even(PercCueExc_AP5[[1]]$bin)==FALSE, c(2, 5, 6)]
inhCont_table_comp2 <- PercCueExc_AP5[[1]][is.even(PercCueExc_AP5[[1]]$bin)==TRUE, c(2, 5, 6)]

#Chi-sq analysis excitations:
#Bins 1, 3 and 5
critData <- contingency_table_comp1[, -1]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 0.86892, df = 2, p-value = 0.6476


#Bins 2, 4, 6, and 8
critData <- contingency_table_comp2[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 2.9386, df = 2, p-value = 0.2301

#Post-hoc comparisons for % cue-excited units
# Is the distribution of exc/not exc neurons different from the first (-120 to -80) and the third (-40 to 0) bin. 
contTableExc_AP5 <- PercCueExc_AP5[[1]][ , c(3, 4)]

#Bin 1 vs. 3
fisher.test(contTableExc_AP5[c(1, 3), ], alternative="less") #95% CI: 0.000000 19.05124 Odds ratio=2.270241 p= 0.9138
#Bin 3 vs. 5
fisher.test(contTableExc_AP5[c(3, 5), ], alternative="less") #95% CI:0.0000000 3.480343 odds ratio=0.5204374 p= 0.3947

#Bin 2 vs. 4
fisher.test((contTableExc_AP5[c(2, 4), ]), alternative="less") #95% CI: 0.0000000 4.684328 odds ratio= 0.6981769  p=  0.529
#Bin 4 vs. 6
fisher.test((contTableExc_AP5[c(4, 6), ]), alternative="less") #95% CI:  0.000000 1.774594 odds ratio=0.4042873 p=0.1961

#Correct p value all comparisons
p.adjust(p=c( 0.9138, 0.3947, 0.529, 0.1961), method="holm")
# 1.0000 1.0000 1.0000 0.7844


#Chi-sq for inhibitions:
#Bins 1, 3 and 5
critData <- inhCont_table_comp1[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 5.9381, df = 2, p-value = 0.05135

#Bins 2, 4 and 6
critData <- inhCont_table_comp2[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 1.4121, df = 2, p-value = 0.4936

contTableInh_AP5 <- PercCueExc_AP5[[1]][ , c(5, 6)]
#########

### VEH VS AP5
#Bin 1: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH[1, ], contTableExc_AP5[1, ])) #95% CI: 0.1549583 7.0445630  Odds ratio=0.9746103 p= 1

#Bin 2: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH[2, ], contTableExc_AP5[2, ])) #95% CI: 0.291676 18.100732 Odds ratio=1.71855 p= 0.7109

#Bin 3: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH[3, ], contTableExc_AP5[3, ])) #95% CI: 0.3527817 21.5874711 Odds ratio=2.062275 p= 0.4769

#Bin 4: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH[4, ], contTableExc_AP5[4, ])) #95% CI: 0.05841177 6.68175896 Odds ratio=0.7296343  p= 1

#Bin 5: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH[5, ], contTableExc_AP5[5, ])) #95% CI: 0.08174756 9.60940018 Odds ratio=1.03406 p= 1

#Bin 6: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH[6, ], contTableExc_AP5[6, ])) #95% CI: 0.5954174 9.8009064 Odds ratio=2.276913 p= 0.2351

p.adjust(p=c(1, 0.7109, 0.4769, 1, 1, 0.2351), method="holm") #1 1 1 1 1 1

### VEH VS AP5
#Bin 1: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH[1, ], contTableInh_AP5[1, ])) #p= 0.6402

#Bin 2: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH[2, ], contTableInh_AP5[2, ])) #p=0.03342

#Bin 3: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH[3, ], contTableInh_AP5[3, ])) #p=0.009662

#Bin 4: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH[4, ], contTableInh_AP5[4, ])) #p=0.3478

#Bin 5: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH[5, ], contTableInh_AP5[5, ])) #p= 1

#Bin 6: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH[6, ], contTableInh_AP5[6, ])) # p= 0.1932

p.adjust(p=c(0.6402, 0.03342, 0.00966, 0.3478, 1, 0.1932), method="holm")
#1.00000 0.16710 0.05796 1.00000 1.00000 0.77280


##### ACTUALLY, THE SESSIONS WERE 35 TRIALS LONG SO I CAN TREAT TRIAL BINS AS SESSIONS (INDEPENDENT DATA ON EACH BIN, NO CHANCE FOR NEURONS OVERLAPPING BINS)

### VEH SIDE
contingency_table_comp<- PercCueExc_VEH[[1]][, c(2, 3, 4)]
inhCont_table_comp <- PercCueExc_VEH[[1]][, c(2, 5, 6)]

#
fisher.test(contingency_table_comp[, -1]) #p = 0.03005
fisher.test(inhCont_table_comp[, -1]) #p = 0.2596


#EXCITATIONS, CONSECUTIVE BIN COMPARISONS
#Bin 1 vs. 2:
fisher.test(contingency_table_comp[c(1, 2), -1]) #p=1
#Bin 2 vs. 3:
fisher.test(contingency_table_comp[c(2, 3), -1]) #p=1
#Bin 3 vs. 4:
fisher.test(contingency_table_comp[c(3, 4), -1]) #p=0.7165
#Bin 4 vs. 5:
fisher.test(contingency_table_comp[c(4, 5), -1]) #p=1
#Bin 5 vs. 6:
fisher.test(contingency_table_comp[c(5, 6), -1]) #p=0.04653

p.adjust(p=c(1, 1, 0.7165, 1, 0.0465), method="holm")
#1.0000 1.0000 1.0000 1.0000 0.2325

#INHIBITIONS, CONSECUTIVE BIN COMPARISONS
#Bin 1 vs. 2:
fisher.test(inhCont_table_comp[c(1, 2), -1]) #p=0.8496
#Bin 2 vs. 3:
fisher.test(inhCont_table_comp[c(2, 3), -1]) #p=0.6159
#Bin 3 vs. 4:
fisher.test(inhCont_table_comp[c(3, 4), -1]) #p=0.8376
#Bin 4 vs. 5:
fisher.test(inhCont_table_comp[c(4, 5), -1]) #p=0.6089
#Bin 5 vs. 6:
fisher.test(inhCont_table_comp[c(5, 6), -1]) #p=0.1737

p.adjust(p=c(0.8496, 0.6159, 0.8376, 0.6089, 0.1737), method="holm")
#1.0000 1.0000 1.0000 1.0000 0.8685





### AP5 SIDE
contingency_table_comp<- PercCueExc_AP5[[1]][, c(2, 3, 4)]
inhCont_table_comp <- PercCueExc_AP5[[1]][, c(2, 5, 6)]

#
fisher.test(contingency_table_comp[, -1]) #p = 0.5149
fisher.test(inhCont_table_comp[, -1]) #p = 0.164


#EXCITATIONS, CONSECUTIVE BIN COMPARISONS
#Bin 1 vs. 2:
fisher.test(contingency_table_comp[c(1, 2), -1]) #p=0.6583
#Bin 2 vs. 3:
fisher.test(contingency_table_comp[c(2, 3), -1]) #p=1
#Bin 3 vs. 4:
fisher.test(contingency_table_comp[c(3, 4), -1]) #p=0.6596
#Bin 4 vs. 5:
fisher.test(contingency_table_comp[c(4, 5), -1]) #p=1
#Bin 5 vs. 6:
fisher.test(contingency_table_comp[c(5, 6), -1]) #p=0.2867

p.adjust(p=c(0.6583, 1, 0.6596, 1, 0.2867), method="holm")
#1.0000 1.0000 1.0000 1.0000 1.0000

#INHIBITIONS, CONSECUTIVE BIN COMPARISONS
#Bin 1 vs. 2:
fisher.test(inhCont_table_comp[c(1, 2), -1]) #p=0.2673
#Bin 2 vs. 3:
fisher.test(inhCont_table_comp[c(2, 3), -1]) #p=0.5337
#Bin 3 vs. 4:
fisher.test(inhCont_table_comp[c(3, 4), -1]) #p=0.1001
#Bin 4 vs. 5:
fisher.test(inhCont_table_comp[c(4, 5), -1]) #p=0.6303
#Bin 5 vs. 6:
fisher.test(inhCont_table_comp[c(5, 6), -1]) #p=0.7998

p.adjust(p=c(0.2673, 0.5337, 0.1001, 0.6303, 0.7998), method="holm")
#1.0000 1.0000 0.5005 1.0000 1.0000



###########################
#BEFORE ENTRY
plotFRandCP_Entry(experiment="Exp 4 Non learners VEH side", cue=c("S+", "S-"), 
                  masterDF=list(masterDF_DSEntry_VEH_NL, masterDF_NSEntry_VEH_NL), graphFolder=MixedGraphFolder, 
            trialBinSize=35, WdwStart=-2000, WdwEnd=0, dataProcess="Zscores", correctOnly=FALSE, 
            colindx =c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), capped=T, capValue = c(0, 210), 
            yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS_VEH_NL, wdwLabel="Pre S+ Entry")

plotFRandCP_Entry(experiment="Exp 4 Non learners AP5 side", cue=c("S+", "S-"), 
                  masterDF=list(masterDF_DSEntry_AP5_NL, masterDF_NSEntry_AP5_NL), graphFolder=MixedGraphFolder, 
                  trialBinSize=35, WdwStart=-2000, WdwEnd=0, dataProcess="Zscores", correctOnly=FALSE, 
                  colindx =c(colindx[2], "darkred"), legLabels=c("S+", "S-"), capped=T, capValue = c(0, 210), 
                  yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS_AP5, wdwLabel="Pre S+ Entry")

plotFRandCP_Entry(experiment="Exp 4 Non learners VEH vs AP5 side", cue=c("S+", "S+"), 
                  masterDF=list(masterDF_DSEntry_VEH_NL, masterDF_DSEntry_AP5_NL), graphFolder=MixedGraphFolder, 
                  trialBinSize=35, WdwStart=-2000, WdwEnd=0, dataProcess="Zscores", correctOnly=FALSE, 
                  colindx = colindx, legLabels=c("S+", "S-"), capped=T, capValue = c(0, 210), 
                   yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS_AP5, wdwLabel="Pre S+ Entry")



#AFTER ENTRY
plotFRandCP_Entry(experiment="Exp 4 Non learners VEH side", cue=c("S+", "S-"), 
                  masterDF=list(masterDF_DSEntry_VEH_NL, masterDF_NSEntry_VEH_NL), graphFolder=MixedGraphFolder, 
                  trialBinSize=35, WdwStart=0, WdwEnd=1000, dataProcess="Zscores", correctOnly=FALSE, 
                  colindx =c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), capped=T, capValue = c(0, 210), 
                  yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS_VEH_NL, wdwLabel="Post S+ Entry")

plotFRandCP_Entry(experiment="Exp 4 Non learners AP5 side", cue=c("S+", "S-"), 
                  masterDF=list(masterDF_DSEntry_AP5_NL, masterDF_NSEntry_AP5_NL), graphFolder=MixedGraphFolder, 
                  trialBinSize=35, WdwStart=0, WdwEnd=1000, dataProcess="Zscores", correctOnly=FALSE, 
                  colindx =c(colindx[2], "darkred"), legLabels=c("S+", "S-"), capped=T, capValue = c(0, 210), 
                  yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS_AP5_NL, wdwLabel="Post S+ Entry")

plotFRandCP_Entry(experiment="Exp 4 Non learners VEH vs AP5 side", cue=c("S+", "S+"), 
                  masterDF=list(masterDF_DSEntry_VEH_NL, masterDF_DSEntry_AP5_NL), graphFolder=MixedGraphFolder, 
                  trialBinSize=35, WdwStart=0, WdwEnd=1000, dataProcess="Zscores", correctOnly=FALSE, 
                  colindx = colindx, legLabels=c("S+", "S-"), capped=T, capValue = c(0, 210), 
                  yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS_AP5_NL, wdwLabel="Post S+ Entry")


### CALCULATE DIFFERENCE IN FR AROUND S+ VS ROUND AS A FUNCTION OF DISTANCE TO CHANGE POINT

giveStars <- function(p.vals){
        sapply(seq(1, length(p.vals)), function(x){
                if(p.vals[x]>0.05){a <- " "}
                if(p.vals[x]<=0.05){a <- "*"}
                if(p.vals[x]<=0.01){a <- "**"}
                if(p.vals[x]<=0.001){a <- "***"}
                return(a)
        })
}

#Post cue VEH
VEH_35bin <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 4 Non learners VEH", 
                                masterDF=list(masterDF_DS_VEH_NL, masterDF_NS_VEH_NL), 
                                   graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="Zscores", 
                                   correctOnly=FALSE, points=F, lines=F, color=c(colindx[1], colindxB[1]), legLabels=c("S+", "S-"), 
                                   yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, capValue=c(0, 210), cueExcOnly=F, 
                                neudata=allNeuronsDS_NL)

VEH_35bin_CueExc <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 4 Non learners VEH Cue Exc Only", 
                                masterDF=list(masterDF_DS_VEH_NL, masterDF_NS_VEH_NL), 
                                graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="Zscores", 
                                correctOnly=FALSE, points=F, lines=F, color=c(colindx[1], colindxB[1]), legLabels=c("S+", "S-"), 
                                yAxMinZ=-2, yAxMaxZ=20, yAxMaxRaw=7, capped=T, capValue=c(0, 210), cueExcOnly=T, 
                                neudata=allNeuronsDS_NL)


uniqueBin <- unique(VEH_35bin[[1]]$bin)

Diff_DSNS_VEH_35bins <- do.call("rbind", lapply(seq(1, length(uniqueBin)), function(i){
        DSbinSel <- VEH_35bin[[1]][VEH_35bin[[1]]$bin==uniqueBin[i], ]
        NSbinSel <- VEH_35bin[[2]][VEH_35bin[[2]]$bin==uniqueBin[i], ]
        
        matchUnitIdx <- match(DSbinSel$unitIdx, NSbinSel$unitIdx)
        NSfiring <- NSbinSel$byUnitFR[matchUnitIdx]
        
        DSbinSel$byUnitFR_NS <- NSfiring
        DSbinSel$DiffFR <- DSbinSel$byUnitFR-NSfiring
        
        DSbinSel
})
)


#This one compares S+ firing rate among bins only (and make a boxplot of it)
compareBins(data=Diff_DSNS_VEH_35bins, cueExcOnly=FALSE, color=colindx[1], ymin=-2, ymax=8, 
            graphFolder=MixedGraphFolder, experiment="Exp 4 Non learners VEH", points=F, 
            comparisons=list(c(1, 3), c(3, 5)), 
            cue="S+")


# [[1]]
# idx     bin         n                      comp        W       p.val      p.adj sig
# W     1 1 vs. 3 56 vs. 90     0 to 35 vs. 70 to 105 2108 0.081674700 0.24502410    
# W2    2 2 vs. 4 86 vs. 46   35 to 70 vs. 105 to 140 1906 0.407546610 0.40754661    
# W1    3 3 vs. 5 90 vs. 31  70 to 105 vs. 140 to 175 1158 0.092586343 0.24502410    
# W11   4 4 vs. 6 46 vs. 34 105 to 140 vs. 175 to 210  514 0.004350405 0.01740162   *
#         
#         [[2]]
# bins         n               comparison          W      p.val     p.adj sig
# W  1 vs. 3 56 vs. 90    0 to 35 vs. 70 to 105 2108 0.08167470 0.1633494    
# W1 3 vs. 5 90 vs. 31 70 to 105 vs. 140 to 175 1158 0.09258634 0.1633494    



PostCueFR_DSvsNS_fromCP_VEH_NL <- compareDSvsNSfromCP(masterDF=list(masterDF_DS_VEH_NL, masterDF_NS_VEH_NL), 
                                                   trialBinSize=35, event="cue", correctOnly=F, paired=T, 
                                                   WdwStart=100, WdwEnd=400, capped=T, capValue=c(0, 210), 
                                                   dataProcess="Zscores")

PostCueFR_DSvsNS_fromCP_VEH_NL$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP_VEH_NL$p, method="holm") 
PostCueFR_DSvsNS_fromCP_VEH_NL$sig <- giveStars(PostCueFR_DSvsNS_fromCP_VEH_NL$p.adj)
# bin           V                p  n      p.adj        sig
# 0 to 35       853     0.244710435 55  0.97884174    
# 35 to 70      1794    0.735950204 81  1.00000000    
# 70 to 105     2538    0.002038770 86  0.01223262      *
# 105 to 140    471     0.927807761 38  1.00000000    
# 140 to 175    80      0.415863037 18  1.00000000    
# 175 to 210    274     0.005522653 26  0.02761327      *


#Post cue AP5
AP5_35bin <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 4 AP5 Non learners", 
                                masterDF=list(masterDF_DS_AP5_NL, masterDF_NS_AP5_NL), 
                                graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="Zscores", 
                                correctOnly=FALSE, points=F, lines=F, color=c(colindx[2], colindxB[2]), legLabels=c("S+", "S-"), 
                                yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, capValue=c(0, 210), cueExcOnly=F, neudata=allNeuronsDS_NL)

AP5_35bin_CueExc <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 4 AP5 Non learners Cue Exc Only", 
                                masterDF=list(masterDF_DS_AP5_NL, masterDF_NS_AP5_NL), 
                                graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="Zscores", 
                                correctOnly=FALSE, points=F, lines=F, color=c(colindx[2], colindxB[2]), legLabels=c("S+", "S-"), 
                                yAxMinZ=-3, yAxMaxZ=20, yAxMaxRaw=7, capped=T, capValue=c(0, 210), cueExcOnly=T, neudata=allNeuronsDS_NL)

uniqueBin <- unique(AP5_35bin[[1]]$bin)

Diff_DSNS_AP5_35bins <- do.call("rbind", lapply(seq(1, length(uniqueBin)), function(i){
        DSbinSel <- AP5_35bin[[1]][AP5_35bin[[1]]$bin==uniqueBin[i], ]
        NSbinSel <- AP5_35bin[[2]][AP5_35bin[[2]]$bin==uniqueBin[i], ]
        
        matchUnitIdx <- match(DSbinSel$unitIdx, NSbinSel$unitIdx)
        NSfiring <- NSbinSel$byUnitFR[matchUnitIdx]
        
        DSbinSel$byUnitFR_NS <- NSfiring
        DSbinSel$DiffFR <- DSbinSel$byUnitFR-NSfiring
        
        DSbinSel
})
)

compareBins(data=Diff_DSNS_AP5_35bins, cueExcOnly=F, color=colindx[2], ymin=-2, ymax=8, 
            graphFolder=MixedGraphFolder, experiment="Exp 4 Non learners AP5", points=F, 
            comparisons=list(c(1, 3), c(3, 5)), 
            cue="S+")

# bin           n                      comp             W       p.val      p.adj        sig
# 1 vs. 3       41 vs. 60     0 to 35 vs. 70 to 105     988 0.280388493 0.84116548    
# 2 vs. 4       48 vs. 51   35 to 70 vs. 105 to 140     1032 0.570176025 0.91428764    
# 3 vs. 5       60 vs. 48  70 to 105 vs. 140 to 175     1216 0.457143819 0.91428764    
# 4 vs. 6       51 vs. 37 105 to 140 vs. 175 to 210     546 0.004573757 0.01829503      *
#         
# [[2]]
# bins          n               comparison                      W       p.val           p.adj sig
# 1 vs. 3       41 vs. 60       0 to 35 vs. 70 to 105           988     0.2803885       0.560777    
# 3 vs. 5       60 vs. 48       70 to 105 vs. 140 to 175        1216    0.4571438       0.560777    




DiffFR_ByBin(data=Diff_DSNS_AP5_35bins, cueExcOnly=FALSE, color=colindx[2], ymin=-2, ymax=18,
             graphFolder=MixedGraphFolder, experiment="Exp 4 Non learners AP5", points=TRUE, 
             comparisons=list(c(1, 3), c(3, 5)))

# [[1]]
# idx     bin                      comp   W       p.val       p.adj sig
# 1 1 vs. 3     0 to 35 vs. 70 to 105 768 0.818530823 1.000000000    
# 2 2 vs. 4   35 to 70 vs. 105 to 140 385 0.956329010 1.000000000    
# 3 3 vs. 5  70 to 105 vs. 140 to 175 377 0.616617768 1.000000000    
# 4 4 vs. 6 105 to 140 vs. 175 to 210  81 0.001169272 0.004677089  **
#         
# [[2]]
# bins               comparison   W     p.val p.adj sig
# 1 vs. 3    0 to 35 vs. 70 to 105 768 0.8185308     1    
# 3 vs. 5 70 to 105 vs. 140 to 175 377 0.6166178     1    



PostCueFR_DSvsNS_fromCP_AP5_NL <- compareDSvsNSfromCP(masterDF=list(masterDF_DS_AP5_NL, masterDF_NS_AP5_NL), cueExcOnly = FALSE,
                                                   trialBinSize=35, event="cue", correctOnly=F, paired=T, 
                                                   WdwStart=100, WdwEnd=400, capped=T, capValue=c(0, 210), 
                                                   dataProcess="Zscores")

PostCueFR_DSvsNS_fromCP_AP5_NL$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP_AP5_NL$p, method="holm") 
PostCueFR_DSvsNS_fromCP_AP5_NL$sig <- giveStars(PostCueFR_DSvsNS_fromCP_AP5_NL$p.adj)

#       bin     V          p    n     p.adj sig
# V     0 to 35 403 0.68202325 38 0.6820233    
# V1   35 to 70 178 0.21946526 24 0.6775757    
# V2  70 to 105 283 0.22067851 36 0.6775757    
# V3 105 to 140 126 0.16939393 25 0.6775757    
# V4 140 to 175  71 0.10808372 20 0.5404186    
# V5 175 to 210  96 0.02062988 15 0.1237793    



#Post cue VEH vs AP5
PostCueFR_DS_fromCP_VEHvsAP5 <- compareVEHvsAP5fromCP(masterDF=list(masterDF_DS_VEH_NL, masterDF_DS_AP5_NL), trialBinSize=35, 
                                                      event="cue", correctOnly=F, cueExcOnly = F, paired=F, WdwStart=100, WdwEnd=400, 
                                                      capped=T, capValue=c(0, 210), dataProcess="Zscores")
PostCueFR_DS_fromCP_VEHvsAP5$p.adj <- p.adjust(PostCueFR_DS_fromCP_VEHvsAP5$p, method="holm") 
PostCueFR_DS_fromCP_VEHvsAP5$sig <- giveStars(PostCueFR_DS_fromCP_VEHvsAP5$p.adj)

# bin                   n    W          p     p.adj sig
# 0 to 35 56 vs.        41 1131 0.25199649 1.0000000    
# 35 to 70 86 vs.       48 1829 0.49798739 1.0000000    
# 70 to 105 90 vs.      60 2835 0.08213281 0.4106640    
# 105 to 140 46 vs.     51 1164 0.26393244 1.0000000    
# 140 to 175 31 vs.     48  823 0.06560493 0.3936296    
# 175 to 210 34 vs.     37  569 0.38052143 1.0000000    



#Pre entry
PreEntryFR_DSvsNS_fromCP_VEH_NL <- compareDSvsNSfromCP(masterDF=list(masterDF_DSEntry_VEH_NL, masterDF_NSEntry_VEH_NL), 
                                                    trialBinSize=35, event="entry", correctOnly=F, WdwStart=-2000, WdwEnd=0, 
                                                    capped=T, capValue=c(0, 210), dataProcess="Zscores")

PreEntryFR_DSvsNS_fromCP_VEH_NL$p.adj <- p.adjust(PreEntryFR_DSvsNS_fromCP_VEH_NL$p, method="holm") 
PreEntryFR_DSvsNS_fromCP_VEH_NL$sig <- giveStars(PreEntryFR_DSvsNS_fromCP_VEH_NL$p.adj)
save(PreEntryFR_DSvsNS_fromCP_VEH_NL, file=paste(dataForRdir, "PreEntryFR_DSvsNS_fromCP_VEH_NL.rdat"))
#       bin        V           p  n       p.adj sig
# V     0 to 35  749 0.431813013 55 1.000000000    
# V1   35 to 70  607 0.386061761 50 1.000000000    
# V2  70 to 105 1015 0.692445548 61 1.000000000    
# V3 105 to 140  133 0.155897141 20 0.779485703    
# V4 140 to 175   37 0.382324219 11 1.000000000    
# V5 175 to 210  233 0.001363158 23 0.008178949  **
        
PreEntryFR_DSvsNS_fromCP_AP5_NL <- compareDSvsNSfromCP(masterDF=list(masterDF_DSEntry_AP5_NL, masterDF_NSEntry_AP5_NL), 
                                                    trialBinSize=35, event="entry", correctOnly=F, WdwStart=-2000, 
                                                    WdwEnd=0, capped=T, capValue=c(0, 210), dataProcess="Zscores")

PreEntryFR_DSvsNS_fromCP_AP5_NL$p.adj <- p.adjust(PreEntryFR_DSvsNS_fromCP_AP5_NL$p, method="holm") 
PreEntryFR_DSvsNS_fromCP_AP5_NL$sig <- giveStars(PreEntryFR_DSvsNS_fromCP_AP5_NL$p.adj)
save(PreEntryFR_DSvsNS_fromCP_AP5_NL, file=paste(dataForRdir, "PreEntryFR_DSvsNS_fromCP_AP5_NL.rdat"))
# bin           V           p   n     p.adj sig
# V     0 to 35 502 0.028419286 38 0.1420964    
# V1   35 to 70 141 0.327819347 22 0.9669091    
# V2  70 to 105 344 0.322303033 35 0.9669091    
# V3 105 to 140 123 0.227937162 24 0.9117486    
# V4 140 to 175 115 0.364253044 20 0.9669091    
# V5 175 to 210 106 0.003356934 15 0.0201416   *

#Post entry
PostEntryFR_DSvsNS_fromCP_VEH_NL <- compareDSvsNSfromCP(masterDF=list(masterDF_DSEntry_VEH_NL, masterDF_NSEntry_VEH_NL), 
                                                        trialBinSize=35, event="entry", correctOnly=F, WdwStart=0, WdwEnd=1000, 
                                                        capped=T, capValue=c(0, 210), dataProcess="Zscores")
PostEntryFR_DSvsNS_fromCP_VEH_NL$p.adj <- p.adjust(PostEntryFR_DSvsNS_fromCP_VEH_NL$p, method="holm") 
PostEntryFR_DSvsNS_fromCP_VEH_NL$sig <- giveStars(PostEntryFR_DSvsNS_fromCP_VEH_NL$p.adj)
save(PostEntryFR_DSvsNS_fromCP_VEH_NL, file=paste(dataForRdir, "PostEntryFR_DSvsNS_fromCP_VEH_NL.rdat"))
# bin           V               p       n       p.adj sig
# 0 to 35       1057    0.0081873729    55 0.040936865   *
# 35 to 70      762     0.1156519998    50 0.231304000    
# 70 to 105     1413    0.0003977494    61 0.002386496  **
# 105 to 140    168     0.0085906982    20 0.040936865   *
# 140 to 175    40      0.2885742187    11 0.288574219    
# 175 to 210    208     0.0163348913    23 0.049004674   *


PostEntryFR_DSvsNS_fromCP_AP5_NL <- compareDSvsNSfromCP(masterDF=list(masterDF_DSEntry_AP5_NL, masterDF_NSEntry_AP5_NL), 
                                                        trialBinSize=35, event="entry", correctOnly=F, 
                                                        WdwStart=0, WdwEnd=1000, capped=T, capValue=c(0, 210), dataProcess="Zscores")
PostEntryFR_DSvsNS_fromCP_AP5_NL$p.adj <- p.adjust(PostEntryFR_DSvsNS_fromCP_AP5_NL$p, method="holm") 
PostEntryFR_DSvsNS_fromCP_AP5_NL$sig <- giveStars(PostEntryFR_DSvsNS_fromCP_AP5_NL$p.adj)
save(PostEntryFR_DSvsNS_fromCP_AP5_NL, file=paste(dataForRdir, "PostEntryFR_DSvsNS_fromCP_AP5_NL.rdat"))
# bin            V            p  n        p.adj sig
# V     0 to 35 381 5.626400e-01 38 0.5626399890    
# V1   35 to 70 223 4.706383e-04 22 0.0018825531  **
# V2  70 to 105 536 6.962306e-05 35 0.0003481153 ***
# V3 105 to 140 224 1.699257e-02 24 0.0339851379   *
# V4 140 to 175 183 1.162529e-03 20 0.0034875870  **
# V5 175 to 210 120 3.051758e-05 15 0.0001831055 ***



### PLOT FR IN THE POSTCUE WINDOW AND TAIL PERIOD BY RAT BY TRIAL BY NEURON
PLOT_perRatTrial_FRandBEH(masterDFsumm=masterDFsumm, wdw="PostCue", threshold=2, MixedGraphFolder=MixedGraphFolder)
PLOT_perRatTrial_FRandBEH(masterDFsumm=masterDFsumm, wdw="PostSpike", threshold=2, MixedGraphFolder=MixedGraphFolder)


### PLOT FR AROUND THE CUE (PTSH) as a function of distance to change point
plotFRandCPhistogram(experiment="Exp 4 Non learners", masterDF=list(masterDF_DS_VEH_NL, masterDF_DS_AP5_NL), 
                     graphFolder=MixedGraphFolder, 
                     trialBinSize=35, dataProcess="Zscores", correctOnly=FALSE, color=colindx,
                     capped=T, capValue = c(0, 210), yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 3, 
                     psthmin=0.5, psthmax=2, imgFormat="pdf", neudata=allNeuronsDS_VEH_NL)


#Same thing but by session from CP instead of trial
plotPSTHfromSessCP(experiment="Exp 4 Non learners", masterDF=list(masterDF_DS_VEH_NL, masterDF_DS_AP5_NL), 
                          graphFolder=MixedGraphFolder, dataProcess="Zscores", comp=c("VEH", "AP5"),
                          correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(0, 5), 
                          yAxMinZ = -1, yAxMaxZ = 4, yAxMaxRaw = 10, psthmin=-0.5, psthmax=2, 
                          imgFormat="pdf", neudata=allNeuronsDS_VEH_NL)

### PLOT FR AROUND ENTRIES (PTSH) as a function of distance to change point
plotFRandCPhistogram_Entry(experiment="Exp 4", masterDF=list(masterDF_DSEntry_VEH, masterDF_DSEntry_AP5), graphFolder=MixedGraphFolder, 
                     trialBinSize=30, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
                     capped=T, capValue = c(-90, 90), yAxMinZ = -1, yAxMaxZ = 2, yAxMaxRaw = 3, 
                     WdwStart=-2000, WdwEnd=2000, imgFormat="pdf", neudata=allNeuronsEntryDS)



### PLOT BOXPLOTS AROUND CUE BY SESS FROM CP AND DRUG GROUP
plotBoxplotfromSessCP(experiment="Exp 4 Non learners", masterDF=list(masterDF_DS_VEH_NL, masterDF_DS_AP5_NL), 
                comp=c("VEH", "AP5"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                correctOnly=FALSE, cueExcOnly=F, color=colindx, sessFromCP = c(0, 5), 
                yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS_NL, morethanIQR=T)







#COMPARE WITH OTHER GROUPS
# compareCPs(data=list(CPdataExp4, CPdataExp10tr, CPdataHYBunilAP5, CPdataHYBbilAP5) , imgFormat="png", expNames=c("Task1", "Task2", "UnilAP5", "BilAP5"), colindx=c("gray30", "gray45", "gray60", "gray75"), graphFolder=behGraphFolder, minSess=5, graphs=c("Change point", "Non learners", "Dynamic Interval", "Asymptote"))


# FOCUS ON THE UNITS RECORDED ON THE DAY BEHAVIOR CHANGED
UnitHeatMap(data=masterDF, sessFromCP=4, FRparameters=allNeuronsDS$parameters, 
            folder=BySessFolder, winmin=0, winmax=400, BLmin=-2000, BLmax=0)



# Megaplot with firing in the 400ms window after 4 events: S+ responded to, S+ missed, S- responded to, S- missed
# This function makes the plot and also saves a data frame that will be useful for the next graph
megaplot(data=list(allNeuronsDSresponded_VEH, allNeuronsDSmissed_VEH, allNeuronsNSresponded_VEH, allNeuronsNSmissed_VEH), 
         CPdata=CPdata, csacqidx=csacqidx, BLwdw=2, winmin=0, winmax=400,
         colpalette="Rainbow", minFR=-1, maxFR=3, graphFolder=MixedGraphFolder, dataForRdir=dataForRdir, 
         ZcolLabels=c("ZDSresp", "ZDSmissed", "ZNSresp", "ZNSmissed"))
        





## SCATTERPLOT WHERE EACH SESSION IS A DOT AND THE AXES ARE --> X axis: performance index; Y axis: mean FR. 

### VEH SIDE
#Run megaplot for VEH, then load toplot, then run the third line line:
megaplot(data=list(allNeuronsDS_VEH, allNeuronsNS_VEH), CPdata=CPdata, csacqidx=csacqidx, BLwdw=2, winmin=0, winmax=400,
         colpalette="Rainbow", minFR=-2, maxFR=6, graphFolder=MixedGraphFolder, dataForRdir=dataForRdir,
         ZcolLabels=c("ZDS", "ZNS"), arrangeBy="ZDS")
load(file=paste(dataForRdir, "toplot.rdat", sep=""))
toplotVEH <- toplot

### Prepare objects to make scatterplots
toplotVEH$uniqSession <- paste(toplotVEH$rat, toplotVEH$expt)

### If that's what you want, subset the CUE-EXCITED neurons only
toplotVEH <- toplotVEH[toplotVEH$DSExc==TRUE, ]

toplotSumm_VEH <- do.call("rbind", lapply(seq(1, length(unique(toplotVEH$uniqSession))), function(x){
        selSess <- toplotVEH[toplotVEH$uniqSession==unique(toplotVEH$uniqSession)[x], ]
        meanZDS <- mean(selSess$ZDS, na.rm=T)
        meanZNS <- mean(selSess$ZNS, na.rm=T)
        sessSumm <- selSess[1, ]
        sessSumm$ZDS <- meanZDS
        sessSumm$ZNS <- meanZNS
        sessSumm$units <- nrow(selSess)
        return(sessSumm)
})
)

toplotSumm_VEH <- as.data.frame(toplotSumm_VEH)
CPsessIndex <- toplotSumm_VEH$BeforeCP==0

#Performance index
FRandPerf_Scatterplot(data=toplotSumm_VEH, graphFolder=ScatterplotFolder, CPtoo=TRUE, 
                      CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="PerfIndex")

#Latency
FRandPerf_Scatterplot(data=toplotSumm_VEH, graphFolder=ScatterplotFolder, CPtoo=TRUE, 
                      CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="Latency")

#Response ratio
FRandPerf_Scatterplot(data=toplotSumm_VEH, graphFolder=ScatterplotFolder, CPtoo=TRUE, 
                      CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="RespRatio")


#Run megaplot for AP5, then load toplot, then run this line:
megaplot(data=list(allNeuronsDS_AP5, allNeuronsNS_AP5), CPdata=CPdata, csacqidx=csacqidx, BLwdw=2, winmin=0, winmax=400,
         colpalette="Rainbow", minFR=-2, maxFR=6, graphFolder=MixedGraphFolder, dataForRdir=dataForRdir,
         ZcolLabels=c("ZDS", "ZNS"), arrangeBy="ZDS")
load(file=paste(dataForRdir, "toplot.rdat", sep=""))
toplotAP5 <- toplot


### Prepare objects to make scatterplots
toplotAP5$uniqSession <- paste(toplotAP5$rat, toplotAP5$expt)

### If that's what you want, subset the CUE-EXCITED neurons only
toplotAP5 <- toplotAP5[toplotAP5$DSExc==TRUE, ]

toplotSumm_AP5 <- do.call("rbind", lapply(seq(1, length(unique(toplotAP5$uniqSession))), function(x){
        selSess <- toplotAP5[toplotAP5$uniqSession==unique(toplotAP5$uniqSession)[x], ]
        meanZDS <- mean(selSess$ZDS, na.rm=T)
        meanZNS <- mean(selSess$ZNS, na.rm=T)
        sessSumm <- selSess[1, ]
        sessSumm$ZDS <- meanZDS
        sessSumm$ZNS <- meanZNS
        sessSumm$units <- nrow(selSess)
        return(sessSumm)
})
)

toplotSumm_AP5 <- as.data.frame(toplotSumm_AP5)
CPsessIndex <- toplotSumm_AP5$BeforeCP==0

#Performance index
FRandPerf_Scatterplot(data=toplotSumm_AP5, graphFolder=ScatterplotFolder, CPtoo=TRUE, 
                      CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="PerfIndex")

#Latency
FRandPerf_Scatterplot(data=toplotSumm_AP5, graphFolder=ScatterplotFolder, CPtoo=TRUE, 
                      CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="Latency")

#Response ratio
FRandPerf_Scatterplot(data=toplotSumm_AP5, graphFolder=ScatterplotFolder, CPtoo=TRUE, 
                      CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="RespRatio")


# % OF BINS EXCITED/INHIBITED

#Matrix in which rows are 50ms bins after the cue, columns are individual neurons and the values indicate if the neuron was EXCITED (ExcBins) or INHIBITED (InhBins) on that bin
ExcBins <- KC.sigbins(path=NEXfiles, startt=0, endt=10000, event=1, BLwdw=2, PostEvent_wdw=3, pbin=0.05, funcdirect=funcdirect)
InhBins <- KC.inhib.sigbins(path=NEXfiles, startt=0, endt=10000, event=1, BLwdw=2, PostEvent_wdw=3, pbin=0.05, funcdirect=funcdirect)
        
## 


### SESS FROM CP DATA

sessMinus1and1 <- sessFromCPdata(data=masterDF_DS_VEH, sessFromCP=c(-1, 1), FRparameters=allNeuronsDS_VEH$parameters, 
               winmin=0, winmax=400, BLmin=-2000, BLmax=0)

allSess_DS_VEH <- sessFromCPdata(data=masterDF_DS_VEH, sessFromCP=c(-4:4), FRparameters=allNeuronsDS_VEH$parameters, 
                                 winmin=0, winmax=400, BLmin=-2000, BLmax=0)

allSess_DS_AP5 <- sessFromCPdata(data=masterDF_DS_AP5, sessFromCP=c(-4:4), FRparameters=allNeuronsDS_AP5$parameters, 
                                 winmin=0, winmax=400, BLmin=-2000, BLmax=0)


allSess_NS_VEH <- sessFromCPdata(data=masterDF_NS, sessFromCP=c(-4:4), FRparameters=allNeuronsNS_VEH$parameters, 
                             winmin=0, winmax=400, BLmin=-2000, BLmax=0)

#
byRat_VEHvsAP5bySess <- summBySess(data=list(allSess_DS_VEH, allSess_DS_AP5), compLabels=c("VEH", "AP5"),
                             CPonly=F, cueExcOnly=F, groupPrePost=F, results="byRat")

byUnit_VEHvsAP5bySess <- summBySess(data=list(allSess_DS_VEH, allSess_DS_AP5), compLabels=c("VEH", "AP5"),
                                   CPonly=F, cueExcOnly=F, groupPrePost=F, results="byUnit")


lapply(seq(1, length(byRat_VEHvsAP5bySess[[1]])), function(d){
        rbind(byRat_VEHvsAP5bySess[[1]][[d]], byRat_VEHvsAP5bySess[[2]][[d]])
})  

lapply(seq(1, length(byUnit_VEHvsAP5bySess[[1]])), function(d){
        rbind(byUnit_VEHvsAP5bySess[[1]][[d]], byUnit_VEHvsAP5bySess[[2]][[d]])
})  

lapply(seq(1, length(byUnit_VEHvsAP5bySess)), function(c){
        
        PreCPsess <- 1:4
        CPsess <- 5
        PostCPsess <- 6:9
        
        plot.new()
        plot.window(xlim=c(0, 2), ylim=c(-4, 8))
        
        lapply(seq(1, length(PreCPsess)), function(d){
                data <- byUnit_VEHvsAP5bySess[[c]][[PreCPsess[d]]]
                data <- data[data$CueExcSess==T, ]
                points(x=rep(c, nrow(data)), y=data$FR_Zsc, pch=19)
                summ <- summary(data$FR_Zsc)
                rect(xleft = c-0.05, xright = c+0.05, ybottom = summ[2], ytop = summ[5])
                
        })
        
        lapply(seq(1, length(PostCPsess)), function(d){
                data <- byUnit_VEHvsAP5bySess[[c]][[PostCPsess[d]]]
                data <- data[data$CueExcSess==T, ]
                points(x=rep(c, nrow(data)), y=data$FR_Zsc, pch=19)
                summ <- summary(data$FR_Zsc)
                rect(xleft = c-0.05, xright = c+0.05, ybottom = summ[2], ytop = summ[5])
                
        })
})
        

#Index of how many units animals contribute to each session
unitsPerSessfromCP <- list(
        t(do.call("rbind", lapply(seq(1, 9), function(x){
                table(allSess_DS_VEH[[x]]$rat, allSess_DS_VEH[[x]]$TrialNumber)[,1]}))), 
        
        t(do.call("rbind", lapply(seq(1, 9), function(x){
                table(allSess_DS_AP5[[x]]$rat, allSess_DS_AP5[[x]]$TrialNumber)[,1]})))
        )

#I'm going to calculate the avg FR on each trial (pool data from different units and animals together)
DS_meanFRperTrial_fromCPsess <- lapply(seq(1, length(allSess_DS)), function(m){
        
        sessFromCP <- unique(allSess_DS[[m]]$sessFromCP)
        
        trialIdx <- unique(allSess_DS[[m]]$TrialNumber)
        
        dat <- sapply(seq(1, length(trialIdx)), function(n){
                selTrial <- allSess_DS[[m]][allSess_DS[[m]]$TrialNumber==n, ]
                meanFRperTrial <- mean(selTrial$Fr_Zsc, na.rm=T)
                
        })
        
        data.frame(trial=trialIdx, sessFromCP=sessFromCP, FR_Zsc=dat)
})

DS_meanFRperTrial_fromCPsess <- do.call("rbind", DS_meanFRperTrial_fromCPsess)
DS_meanFRperTrial_fromCPsess$idx <- 1:nrow(DS_meanFRperTrial_fromCPsess)


NS_meanFRperTrial_fromCPsess <- lapply(seq(1, length(allSess_NS)), function(m){
        
        sessFromCP <- unique(allSess_NS[[m]]$sessFromCP)
        
        trialIdx <- unique(allSess_NS[[m]]$TrialNumber)
        
        dat <- sapply(seq(1, length(trialIdx)), function(n){
                selTrial <- allSess_NS[[m]][allSess_NS[[m]]$TrialNumber==n, ]
                meanFRperTrial <- mean(selTrial$Fr_Zsc, na.rm=T)
                
        })
        
        data.frame(trial=trialIdx, sessFromCP=sessFromCP, FR_Zsc=dat)
})

NS_meanFRperTrial_fromCPsess <- do.call("rbind", NS_meanFRperTrial_fromCPsess)
NS_meanFRperTrial_fromCPsess$idx <- 1:nrow(NS_meanFRperTrial_fromCPsess)


meanFRperTrial_fromCPsess <- cbind(DS_meanFRperTrial_fromCPsess, NS_meanFRperTrial_fromCPsess$FR_Zsc)
colnames(meanFRperTrial_fromCPsess) <- c("trial", "sessFromCP", "DS_FR_Zsc", "idx", "NS_FR_Zsc")

#Make long format data frame with this data
LF_meanFRperTrial_fromCPsess <- melt(meanFRperTrial_fromCPsess, id.vars=c("sessFromCP", "trial", "idx"))    
        
ezANOVA(data=LF_meanFRperTrial_fromCPsess, dv=value, between=sessFromCP, within=variable, wid=idx, type=3)







####

        