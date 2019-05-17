
setwd("E:/Dropbox/NMDA/")

#################################################################################
### COMBINING ALL DATA OF RECORDINGS + NO INFUSIONS (DOWN AND NOT DOWN RATS)  ###
#################################################################################

### LOAD IMPORTANT LIBRARIES
install.packages("matrixStats")
install.packages('ez')

library(matrixStats)
library(ez)

setwd("E:/Dropbox (EinsteinMed)/NMDA/")

#Define colors
colindx <- c("#2171b5", "#cb181d") #Strong blue and red
colindxB <- c("#bdd7e7", "#fcae91") #Less strong blue and red
colindxC <- c("#eff3ff", "#fb6a4a") #Even less strong blue and red
colindxD <- c("#6baed6", "#fee5d9") #Lightest blue and red


### DEFINE ALL THE IMPORTANT FOLDERS

funcdirect <- paste(getwd(), "/R functions/", sep="")
datafolder <- paste(getwd(), "/EXP1B_Not down/MedPC files/", sep="")
dataForRdir <- paste(getwd(), "/EXP1B_Not down/Data for R/", sep="")
dataForRCumulative <- paste(getwd(), "/EXP1B_Not down/Data for R cumulative/", sep="")
behGraphFolder <- paste(getwd(), "/EXP1B_Not down/Graphs/Behavior/", sep="")
CPfuncFolder <- paste(funcdirect, 'Change_Point-master/', sep="")
CPGraphFolder <- paste(getwd(), "/EXP1B_Not down/Graphs/Behavior/Change point/", sep="")
MixedGraphFolder <- paste(getwd(), "/EXP1B_Not down/Graphs/Mixed/", sep="")
NEXfiles <- paste(getwd(), "/EXP1B_Not down/NEX files/", sep="")
IndcumPerfGraphFolder <- paste(getwd(), "/EXP1B_Not down/Graphs/Behavior/Cumulative perf/", sep="")
PerfRelToCPFolder <- paste(getwd(), "/EXP1B_Not down/Graphs/Behavior/Beh rel to CP/", sep="")
BySessFolder <- paste(MixedGraphFolder, "By Sess from CP/", sep="")
ScatterplotFolder <- paste(MixedGraphFolder, "Scatterplot/", sep="")
otherGraphFolder <- paste(getwd(), "/EXP1B_Not down/Graphs/Others/", sep="")

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
load(file=paste(funcdirect, "plotFRandCPhistogram.r", sep=""))
load(file=paste(funcdirect, "megaplot.r", sep=""))
load(file=paste(funcdirect, "sessFromCPdata.r", sep=""))
load(file=paste(funcdirect, "compareDSvsNSfromCP.R", sep=""))
load(file=paste(funcdirect, "plotFRBoxPlotandCP.R", sep = ""))
load(file=paste(funcdirect, "DiffFR_ByBin.R", sep=""))
load(file=paste(funcdirect, "cumulative_CP.R", sep=""))

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
MAKERASTER(i=x, data=alldata, idxdata=csacqidx) #"i" indicates the number of the session in "alldata". To know what day of what animal it is, that info is in "csacqidx" (or in "idx")
MAKERASTER(i=4, data=alldata, idxdata=csacqidx)


###
# I'M USING THE 10S PRE-CUE PERIOD AS BASELINE
# In 9% of the trials, the animal's head was inside the receptacle at the time the ITI window started. In those trials, I'll correct the ITI latency by assigning the MEAN ITI latency in the +/- 5 trials around that trial
###

#AVERAGE PERFORMANCE S+ VS S-. 5 trial bins.
avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(DStaskAcc, NStaskAcc), 
             cues=c("S+", "S-"), index="Specificity", legendLocation="topleft", y_axis_label="Performance index", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(DSrespAll, NSrespAll), 
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


#CPdata with 10s precue ITI + moving average
# rat  CP CPsess   slopePre slopePost FirstCPOver80pc DynamicInterv nSess NS.PreCP
# 1  MV42 156      5  0.1423823  4.638940             156             0     6      151
# 2  MV41 120      4 -0.2282295  2.877726             120             0     6      123
# 3  MV43 164      5  0.3007178  3.800131             164             0     6      152
# 4  MV44 166      5  0.3516210  1.433214             166             0     6      161
# 5  MV45  92      3  0.2927708  3.051936             107            15     6       87
# 6  MV46  61      2  0.1178718  4.431484              73            12     6       60
# 7  MV48 157      5 -0.7796931  1.658166             157             0     6      153
# 8  MV50 124      5  1.1235425  1.494282             143            19     6      120
# 9  MV53 174      5  0.4674086  2.168733             174             0     6      171
# 10 MV56  40      2 -0.2912341  3.914329              75            35     6       37
# 11 MV57  83      3 -0.2497722  3.711857             116            33     5       77
# 12 MV60  58      2 -0.3351285  2.735201              82            24     6       54
# 13 MV59 123      4  0.2740007  2.981928             123             0     6      121

CPdataND <- CPdata
alldataND <- alldata
csacqidxND <- csacqidx
ratsND <- rats
idxND <- idx

save(csacqidxND, file=paste(dataForRdir, "csacqidxND.ridx", sep=""))
save(alldataND, file=paste(dataForRdir, "alldataND.rdat", sep=""))
save(ratsND, file=paste(dataForRdir, "ratsND.rdat", sep=""))
save(idxND, file=paste(dataForRdir, "idxND.rdat", sep=""))


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
CPframeND <- rbind(DSspecificity, RealCPbyIndex[sel,])

save(CPframeND, file=paste(dataForRdir, "CPframeND.rdat", sep="")) #CP values for different behavioral indexes
save(CPdataND, file=paste(dataForRdir, "CPdataND.rdat", sep="")) #CP values and other details for PERFORMANCE INDEX (S+ specificity)


#CUMULATIVE INDIVIDUAL PERFORMANCE
cumulativeIndGraphs(numSess=6, numCues=35, sessLines=F, smartRat="NA", limit=210, dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=IndcumPerfGraphFolder, imgFormat="pdf")

#MEAN PERFORMANCE AS A FUNCTION OF DISTANCE FROM CP
PerformanceFromCP(relTrialMin=-100, relTrialMax=100, trialBinSize=5, typegraph="both", imgFormat="pdf", 
                  dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=PerfRelToCPFolder, CPdata=CPdata, 
                  idx=idx, rats = rats, alldata = alldata)


#PRE CP VS. POST CP PERFORMANCE

# S+ SPECIFICITY
PrevsPostFolder <- paste(PerfRelToCPFolder, "Average Pre vs Post/", sep="")
PrePostCP_Perf(data=DSrespAll, CP=CPdata$CP, y_axis_label="S+ response ratio", graphFolder=PrevsPostFolder, plot=T) #V = 0, p-value = 0.0002441
PrePostCP_Perf(data=NSrespAll, CP=CPdata$CP, y_axis_label="S- response ratio", graphFolder=PrevsPostFolder, plot=T) #t = 1.4904, df = 12, p-value = 0.1619
PrePostCP_Perf(data=DSlatency, CP=CPdata$CP, y_axis_label="S+ latency", graphFolder=PrevsPostFolder, plot=T) #t = 7.4541, df = 12, p-value = 7.693e-06
PrePostCP_Perf(data=NSlatency, CP=CPdata$CP, y_axis_label="S- latency", graphFolder=PrevsPostFolder, plot=T) #t = -2.6775, df = 12, p-value = 0.02013
PrePostCP_Perf(data=DStaskAcc, CP=CPdata$CP, y_axis_label="S+ specificity", graphFolder=PrevsPostFolder, plot=T) #t = -7.9401, df = 12, p-value = 4.061e-06
PrePostCP_Perf(data=NStaskAcc, CP=CPdata$CP, y_axis_label="S- specificity", graphFolder=PrevsPostFolder, plot=T) #t = 0.85764, df = 12, p-value = 0.4079
PrePostCP_Perf(data=ITIlatency, CP=CPdata$CP, y_axis_label="ITI latency", graphFolder=PrevsPostFolder, plot=T) #t = -1.8576, df = 12, p-value = 0.08793

p.adjust(c(0.0002441, 7.693e-06, 4.061e-06, 0.08793), method="holm")
#4.8820e-04 2.3079e-05 1.6244e-05 8.7930e-02



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

save(ratsTrialIdxPerSess, file=paste(dataForRdir, "ratsTrialIdxPerSess.rdat", sep=""))
save(CPvector, file=paste(dataForRdir, "CPvector.rdat", sep=""))
save(sessionCPperRat, file=paste(dataForRdir, "sessionCPperRat.rdat", sep=""))


### PREDICTORS OF CHANGE POINT
hist(CPdata$CP, breaks = seq(0, 250, 10), col="black", border="white")


CPdata$PreCPRewTrials <- sapply(seq(1, length(rats)), function(x){
        if(!is.na(CPdata$CP[x])){
                rewTrialsPreCP <- sum(DSrespAll[[x]][1:CPdata$CP[x]])
        } else {rewTrialsPreCP <- NA}
        
        rewTrialsPreCP
}) 


CPdata$TotalRewards <- sapply(seq(1, length(rats)), function(x){
        totalrewards <- sum(DSrespAll[[x]], na.rm=T)
        totalrewards
}) 

###TOTAL REWARDS UNDER 2s latency before CP
CPdata$ShortLatRewPreCP <- sapply(seq(1, length(rats)), function(x){
        if(!is.na(CPdata$CP[x])){
                under4strials <- DSlatency[[x]]<=2
                shortLatTr <- sum(under4strials[1:CPdata$CP[x]])
                
        } else {shortLatTr <- NA}
        
        shortLatTr
}) 

CPdata$TotalRewardsShortLat <- sapply(seq(1, length(rats)), function(x){
        totalrewardsSL <- sum(DSlatency[[x]]<=2, na.rm=T)
        totalrewardsSL
})


#TOTAL TRIALS TO CP VS. NUMBER OF REWARDS
plot.new()
plot.window(xlim=c(0, 2), ylim=c(0, 216))

rect(xleft=0, xright=0.5, ybottom=0, ytop=mean(CPdata$CP, na.rm=T), col="gray80")
#errBars(x=0.25, y=mean(CPdata$CP, na.rm=T), err=sd(CPdata$CP, na.rm=T)/(sqrt(sum(!is.na(CPdata$CP)))))
points(x=rep(0.25, nrow(CPdata)), y=CPdata$CP, pch=19, cex=1)

rect(xleft=0.5, xright=1, ybottom=0, ytop=mean(CPdata$PreCPRewTrials, na.rm=T), col="gray30")
#errBars(x=0.8, y=mean(CPdata$PreCPRewTrials, na.rm=T), err=sd(CPdata$PreCPRewTrials, na.rm=T)/(sqrt(sum(!is.na(CPdata$PreCPRewTrials)))))
points(x=rep(0.75, nrow(CPdata)), y=CPdata$PreCPRewTrials, pch=19, cex=1)

rect(xleft=1, xright=1.5, ybottom=0, ytop=mean(CPdata$ShortLatRewPreCP, na.rm=T), col="white")
#errBars(x=1.3, y=mean(CPdata$ShortLatRewPreCP, na.rm=T), err=sd(CPdata$ShortLatRewPreCP, na.rm=T)/(sqrt(sum(!is.na(CPdata$ShortLatRewPreCP)))))
points(x=rep(1.25, nrow(CPdata)), y=CPdata$ShortLatRewPreCP, pch=19, cex=1)


for(i in 1:nrow(CPdata)){
        lines(x=c(0.25, 0.75, 1.25), y=c(CPdata$CP[i], CPdata$PreCPRewTrials[i], CPdata$ShortLatRewPreCP[i]))
}

axis(side=2, las=2)
axis(side=1, at=c(0.25, 0.75, 1.25), labels=c("Total # trials to CP", "# Of rewarded trials to CP", "# rewards with < 2 latency"), font=2, cex.axis=1.2)

mtext(text="Trials", side=2, line=2.5, font=2, cex=1.2)
title("Total trials vs. Rewarded trials to Change point")


#####
# CUMULATIVE PROPORTION OF ANIMALS WITH CHANGE POINT BY SESSION

#Index of "Not Down" files in the order they are in "allNeuronsDS"
allNeuronsDS_NDidx <- sapply(seq(1, length(allNeuronsDS$nexdata)), function(x){allNeuronsDS$nexdata[[x]]$down}) == "NotDown"

#Index of "Not Down" rats
rats_NDidx <- unique(sapply(seq(1, length(allNeuronsDS$nexdata)), function(x){allNeuronsDS$nexdata[[x]]$rat})[allNeuronsDS_NDidx])


cumulative_CP(Exp="Not Down", CPdata=CPdata[CPdata$rat %in% rats_NDidx, ], numSess=6, byAnimal=TRUE, byNeuron=FALSE, 
              nexdata=allNeuronsDS$nexdata[c(allNeuronsDS_NDidx)], graphFolder=otherGraphFolder)

cumulative_CP(Exp="Not Down", CPdata=CPdata[CPdata$rat %in% rats_NDidx, ], numSess=6, byAnimal=F, byNeuron=T, 
              nexdata=allNeuronsDS$nexdata[c(allNeuronsDS_NDidx)], graphFolder=otherGraphFolder)






################################################################################
############ NEURONAL DATA #####################################################
################################################################################

allNeuronsDS = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNS = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=T, side="both")

allNeuronsDSresponded = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=5, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="both")
allNeuronsNSresponded = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=7, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="both")

allNeuronsDSmissed = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=6, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="both")
allNeuronsNSmissed = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=8, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side ="both")

### ENTRIES
allNeuronsEntryDS = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T)
allNeuronsEntryNS = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T)
allNeuronsEntryITI = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T)


#sAVE THESE OBJECTS
save(allNeuronsDS, file=paste(dataForRdir, 'allNeuronsDS.rdat', sep=""))
save(allNeuronsNS, file=paste(dataForRdir, 'allNeuronsNS.rdat', sep=""))
save(allNeuronsDSresponded, file=paste(dataForRdir, 'allNeuronsDSresponded.rdat', sep=""))
save(allNeuronsNSresponded, file=paste(dataForRdir, 'allNeuronsNSresponded.rdat', sep=""))
save(allNeuronsDSmissed, file=paste(dataForRdir, 'allNeuronsDSmissed.rdat', sep=""))
save(allNeuronsNSmissed, file=paste(dataForRdir, 'allNeuronsNSmissed.rdat', sep=""))
save(allNeuronsEntryDS, file=paste(dataForRdir, 'allNeuronsEntryDS.rdat', sep=""))
save(allNeuronsEntryNS, file=paste(dataForRdir, 'allNeuronsEntryNS.rdat', sep=""))
save(allNeuronsEntryITI, file=paste(dataForRdir, 'allNeuronsEntryITI.rdat', sep=""))

#Load files again      
#files <- paste(dataForRdir, list.files(dataForRdir), sep=""); for(i in 1:length(files)){load(files[[i]])}

####################33
### CHECK THAT THE CRITERIA FOR CUE EXCITED UNITS MAKES SENSE
#Time course of excitation
toplotTC <- megaplotTC(data=list(allNeuronsDS), graphmin=500, graphmax=5000, cuewinmin=100, cuewinmax=400,
                       colpalette="Rainbow", minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, 
                       ZcolLabels=c("ZDS"), arrangeBy=c("ZDS"))

#Histograms with percentage of cue excited units per bin
CueExcCheck(experiment="Exp Combo no infusions", data=toplotTC, graphFolder=MixedGraphFolder, 
            cexEx=1, cexInh=0.5, autom.FRrange=T, FRbin=0.5)

#BUILD A DATAFRAME WITH FIRING RATE OF EACH NEURON ON EACH TRIAL WITH BEH INFO ABOUT EACH TRIAL
masterDF_DS <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS, CPvector=CPdata$CP, sessionCPperRat = sessionCPperRat, funcdirect=funcdirect, dataForRdir=dataForRdir, BLduration=2, BLneudata=allNeuronsDS)
masterDF_DScorrect <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSresponded, CPvector=CPdata$CP, sessionCPperRat = sessionCPperRat, funcdirect=funcdirect, dataForRdir=dataForRdir, BLduration=2, BLneudata=allNeuronsDS)
masterDF_DSmissed <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSmissed, CPvector=CPdata$CP, sessionCPperRat = sessionCPperRat, funcdirect=funcdirect, dataForRdir=dataForRdir, BLduration=2, BLneudata=allNeuronsDS)
masterDF_NS <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsNS, CPvector=CPdata$NS.PreCP, sessionCPperRat = sessionCPperRat, funcdirect=funcdirect, dataForRdir=dataForRdir, BLduration=2, BLneudata=allNeuronsNS)
masterDF_NSresponded <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNSresponded, 
                                           CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                           BLduration=2, sessionCPperRat=sessionCPperRat, BLneudata=allNeuronsNS)
masterDF_NSmissed<- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNSmissed, 
                                       CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                       BLduration=2, sessionCPperRat=sessionCPperRat, BLneudata=allNeuronsNS)


#The "entry-related" masterDF don't give me FR info of every trial, only of those trials in which the animal responded
masterDF_DSEntry <- FRbyNEURONbyBINcue(eventType="entry", cue="S+", neudata=allNeuronsEntryDS, 
                                       CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                       BLduration=2, sessionCPperRat=sessionCPperRat, BLneudata=allNeuronsDS)
masterDF_NSEntry<- FRbyNEURONbyBINcue(eventType="entry", cue="S+", neudata=allNeuronsEntryNS, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=sessionCPperRat, BLneudata=allNeuronsNS)


#Save masterDF objects
save(masterDF_DS, file=paste(dataForRdir, "masterDF_DS.rdat", sep=""))
save(masterDF_DScorrect, file=paste(dataForRdir, "masterDF_DScorrect.rdat", sep=""))
save(masterDF_DSmissed, file=paste(dataForRdir, "masterDF_DSmissed.rdat", sep=""))
save(masterDF_NS, file=paste(dataForRdir, "masterDF_NS.rdat", sep=""))
save(masterDF_DSEntry, file=paste(dataForRdir, "masterDF_DSEntry.rdat", sep=""))
save(masterDF_NSEntry, file=paste(dataForRdir, "masterDF_NSEntry.rdat", sep=""))
save(masterDF_DS_NotDown, file=paste(dataForRdir, "masterDF_DS_NotDown.rdat", sep=""))
save(masterDF_NS_NotDown, file=paste(dataForRdir, "masterDF_NS_NotDown.rdat", sep=""))


# 
# ### BY MODALITY & DOWN/NOTDOWN
# masterDF_DS_Tone_Down <- filter(masterDF_DS, masterDF_DS$modality==98, masterDF_DS$down=="Down")
# masterDF_DS_Light_Down <- filter(masterDF_DS, masterDF_DS$modality==99, masterDF_DS$down=="Down")
# masterDF_DS_Tone_NotDown <- filter(masterDF_DS, masterDF_DS$modality==98, masterDF_DS$down=="NotDown")
# masterDF_DS_Light_NotDown <- filter(masterDF_DS, masterDF_DS$modality==99, masterDF_DS$down=="NotDown")
# masterDF_NS_Light_NotDown <- filter(masterDF_NS, masterDF_NS$modality==98, masterDF_NS$down=="NotDown")
# masterDF_NS_Tone_NotDown <- filter(masterDF_NS, masterDF_NS$modality==99, masterDF_NS$down=="NotDown")
# 
# 
# ### BY DOWN/NOTDOWN
# masterDF_DS_Down <- filter(masterDF_DS, masterDF_DS$down=="Down")
# masterDF_DS_NotDown <- filter(masterDF_DS, masterDF_DS$down=="NotDown")
# masterDF_NS_Down <- filter(masterDF_NS, masterDF_NS$down=="Down")
# masterDF_NS_NotDown <- filter(masterDF_NS, masterDF_NS$down=="NotDown")
# 
# ### BY TONE/LIGHT
# masterDF_DS_Light <- filter(masterDF_DS, masterDF_DS$modality==99)
# masterDF_DS_Tone <- filter(masterDF_DS, masterDF_DS$modality==98)
# 
# masterDF_NS_Light <- filter(masterDF_NS, masterDF_NS$modality==98)
# masterDF_NS_Tone <- filter(masterDF_NS, masterDF_NS$modality==99)
# 




#masterDFsumm <- masterDFsummary(masterDF=masterDF_DS, neudata=allNeuronsDS, format="Zscores", postCueWdw=c(0, 0.4), postSpikeWdw=c(0.4, 5), orEntry=T, cuelength=10)


### PLOT FR POST-CUE (400MS WINDOW) as a function of distance to change point


#CUE 100-400
plotFRandCP(experiment="Exp 3", masterDF=list(masterDF_DS), graphFolder=MixedGraphFolder, 
            trialBinSize=5, WdwStart=100, WdwEnd=400, dataProcess="raw", correctOnly=FALSE, 
            colindx =colindx[1], capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 10, 
            yAxMaxRaw = 10, neudata=allNeuronsDS)

plotFRandCP(experiment="Exp 3", masterDF=list(masterDF_DS), graphFolder=MixedGraphFolder, 
            trialBinSize=5, WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
            colindx=colindx[1], capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 10, 
            yAxMaxRaw = 10, neudata=allNeuronsDS)

plotFRandCP(experiment="Exp 3", masterDF=list(masterDF_DS), graphFolder=MixedGraphFolder, 
            trialBinSize=5, WdwStart=100, WdwEnd=400, dataProcess="raw", correctOnly=TRUE, 
            colindx=colindx[1], capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 10, 
            yAxMaxRaw = 10, neudata=allNeuronsDS)

plotFRandCP(experiment="Exp 3", masterDF=list(masterDF_DS), graphFolder=MixedGraphFolder, trialBinSize=5, 
            WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=TRUE, colindx=colindx[1], capped=T, 
            capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 10, yAxMaxRaw = 10, neudata=allNeuronsDS)

plotFRandCP(cue=c("S+", "S-"), experiment="Exp 1 S+ vs S-", masterDF=list(masterDF_DS, masterDF_NS), 
            legLabels=c("S+", "S-"), graphFolder=MixedGraphFolder, trialBinSize=5, 
            WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-140, 140), yAxMinZ = -1, yAxMaxZ = 10, yAxMaxRaw = 10, 
            neudata=c(allNeuronsDS, allNeuronsNS), cueExcOnly = F)


plotFRandCP(cue=c("S+", "S-"), experiment="Exp 1 NOT DOWN S+ vs S-", 
            masterDF=list(masterDF_DS_NotDown, masterDF_NS_NotDown), 
            legLabels=c("S+", "S-"), graphFolder=MixedGraphFolder, trialBinSize=5, 
            WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-140, 140), yAxMinZ = -1, yAxMaxZ = 6, yAxMaxRaw = 10, 
            neudata=c(allNeuronsDS, allNeuronsNS), cueExcOnly = F)

PercResults <- plotFRandCP(cue=c("S+"), experiment="Exp 1 NOT DOWN S+", 
                           masterDF=list(masterDF_DS_NotDown), 
                           legLabels=c("S+"), graphFolder=MixedGraphFolder, trialBinSize=35, 
                           WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), capped=T, 
                           capValue = c(-140, 140), yAxMinZ = -1, yAxMaxZ = 6, yAxMaxRaw = 10, 
                           neudata=c(allNeuronsDS), cueExcOnly = F)

save(PercResults, file=paste(dataForRdir, "PercResults.R", sep=""))

# [[1]]
#     bin     trialBins       CueEx   notCueExc       CueInh          notCueInh
# 1   1         -140            20        66            39              47
# 2   2         -105            31        70            16              85
# 3   3         -70             71       107            26              152
# 4   4         -35             83       100            43              140
# 5   5          0              75        68            39              104
# 6   6         35              54        46            23              77
# 7   7         70              23        20            12              31
# 8   8         105             14        11            8               17


plotFRandCP(cue=c("S-"), experiment="Exp 1 NOT DOWN S-", 
            masterDF=list(masterDF_NS_NotDown), 
            legLabels=c("S+"), graphFolder=MixedGraphFolder, trialBinSize=35, 
            WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-140, 140), yAxMinZ = -1, yAxMaxZ = 6, yAxMaxRaw = 10, 
            neudata=c(allNeuronsDS), cueExcOnly = F)
# [[1]]
#     bin trialBins CueEx notCueExc CueInh notCueInh
# 1   1      -140     7        79     42        44
# 2   2      -105     3        98     27        74
# 3   3       -70     2       176     28       150
# 4   4       -35     9       174     36       147
# 5   5         0    10       119     40        89
# 6   6        35     4        86     34        56
# 7   7        70     2        28     10        20
# 8   8       105     1        13      8         6


plotFRandCP(cue=c("S+", "S-"), experiment="Exp 1 NOT DOWN S+ vs S- CUE EXC only", 
            masterDF=list(masterDF_DS_NotDown, masterDF_NS_NotDown), 
            legLabels=c("S+", "S-"), graphFolder=MixedGraphFolder, trialBinSize=5, 
            WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 6, yAxMaxRaw = 10, 
            neudata=c(allNeuronsDS, allNeuronsNS), cueExcOnly = T)


plotFRandCP(cue=c("S+", "S+"), experiment="Exp 3 Tone vs Light S+", 
            masterDF=list(masterDF_DS_Tone, masterDF_DS_Light), 
            legLabels=c("Tone S+", "Light S+"), graphFolder=MixedGraphFolder, trialBinSize=5, 
            WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-100, 100), yAxMinZ = -5, yAxMaxZ = 10, yAxMaxRaw = 10, 
            neudata=allNeuronsDS, cueExcOnly = F)

plotFRandCP(cue=c("S+"), experiment="Exp 3 S+ PERCEXC Tone S+", 
            masterDF=list(masterDF_DS_Tone), 
            legLabels=c("Tone S+"), graphFolder=MixedGraphFolder, trialBinSize=15, 
            WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, 
            colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-90, 90), yAxMinZ = -5, yAxMaxZ = 6, yAxMaxRaw = 10, 
            neudata=allNeuronsDS, cueExcOnly = F)

plotFRandCP(cue=c("S+"), experiment="Exp 3 S+ PERCEXC Light S+", 
            masterDF=list(masterDF_DS_Light), 
            legLabels=c("Tone S+"), graphFolder=MixedGraphFolder, trialBinSize=15, 
            WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, 
            colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-90, 90), yAxMinZ = -5, yAxMaxZ = 6, yAxMaxRaw = 10, 
            neudata=allNeuronsDS, cueExcOnly = F)



plotFRandCP(cue=c("S+", "S-"), experiment="Exp 1 S+ Tone vs S- Light", 
            masterDF=list(masterDF_DS_Tone_NotDown, masterDF_NS_Light_NotDown), 
            legLabels=c("Tone S+", "Light S-"), graphFolder=MixedGraphFolder, trialBinSize=5, 
            WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-120, 120), yAxMinZ = -5, yAxMaxZ = 5, yAxMaxRaw = 10, 
            neudata=allNeuronsDS, cueExcOnly = F)

plotFRandCP(cue=c("S+"), experiment="Exp 1 S+ Tone", 
            masterDF=list(masterDF_DS_Tone_NotDown), 
            legLabels=c("Tone S+"), graphFolder=MixedGraphFolder, trialBinSize=40, 
            WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, colindx=c(colindx[1]), capped=T, 
            capValue = c(-120, 120), yAxMinZ = -5, yAxMaxZ = 5, yAxMaxRaw = 10, 
            neudata=allNeuronsDS, cueExcOnly = F)

plotFRandCP(cue=c("S+"), experiment="Exp 1 S+ Light", 
            masterDF=list(masterDF_DS_Light_NotDown), 
            legLabels=c("Light S+"), graphFolder=MixedGraphFolder, trialBinSize=40, 
            WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, colindx=c(colindx[1]), capped=T, 
            capValue = c(-120, 120), yAxMinZ = -5, yAxMaxZ = 5, yAxMaxRaw = 10, 
            neudata=allNeuronsDS, cueExcOnly = F)




#########
### BOXPLOTS vs. CHANGEPOINT BY BLOCKS OF TRIALS
####

byBinbyUnit <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 NOTDOWN S+ vs S-", 
                                  masterDF=list(masterDF_DS_NotDown, masterDF_NS_NotDown), points=F,
                                  graphFolder=MixedGraphFolder, trialBinSize=35, 
                                  WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                  color=c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), yAxMinZ=-2, 
                                  yAxMaxZ=6, yAxMaxRaw=7, 
                                  capped=T, capValue=c(-140, 140), cueExcOnly=F, neudata=allNeuronsDS)

nUnits <- length(unique(byBinbyUnit[[1]]$unitIdx))



### BY MODALITY
byBinbyUnit_Tone <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 NOTDOWN S+ Tone vs S- Light", 
                                       masterDF=list(masterDF_DS_Tone_NotDown, masterDF_NS_Light_NotDown), points=F,
                                       graphFolder=MixedGraphFolder, trialBinSize=20, 
                                       WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                       color=c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), yAxMinZ=-2, 
                                       yAxMaxZ=6, yAxMaxRaw=7, 
                                       capped=T, capValue=c(-120, 120), cueExcOnly=T, neudata=allNeuronsDS)


byBinbyUnit_Light <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 NOTDOWN S+ Light vs S- Tone", 
                                        masterDF=list(masterDF_DS_Light_NotDown, masterDF_NS_Tone_NotDown), points=F,
                                        graphFolder=MixedGraphFolder, trialBinSize=20, 
                                        WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                        color=c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), yAxMinZ=-2, 
                                        yAxMaxZ=6, yAxMaxRaw=7, 
                                        capped=T, capValue=c(-120, 120), cueExcOnly=T, neudata=allNeuronsDS)





#### PLOT FR IN ORDER from trial 1-n
CP1 <- CPdata
CP1$CP <- 1
CP1$CPsess <- 1

masterDF_DS_From1 <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS, 
                                        CPvector=CP1$CP, sessionCPperRat = CP1$CPsess, 
                                        funcdirect=funcdirect, dataForRdir=dataForRdir, BLduration=2, 
                                        BLneudata=allNeuronsDS)

masterDF_NS_From1 <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS, 
                                        CPvector=CP1$CP, sessionCPperRat = CP1$CPsess, 
                                        funcdirect=funcdirect, dataForRdir=dataForRdir, BLduration=2, 
                                        BLneudata=allNeuronsNS)

realCPsessDS <- sapply(seq(1, nrow(masterDF_DS_From1)), function(n){
        rat <- as.character(masterDF_DS_From1$rat[n])
        realCPsess <- CPdata[CPdata$rat==rat, ]$CPsess
        realCPsess
})

realCPsessNS <- sapply(seq(1, nrow(masterDF_NS_From1)), function(n){
        rat <- as.character(masterDF_NS_From1$rat[n])
        realCPsess <- CPdata[CPdata$rat==rat, ]$CPsess
        realCPsess
})

masterDF_DS_From1$realCPsess <- realCPsessDS
masterDF_NS_From1$realCPsess <- realCPsessNS        

save(masterDF_DS_From1, file=paste(dataForRdir, "masterDF_DS_From1.rdat", sep=""))
save(masterDF_NS_From1, file=paste(dataForRdir, "masterDF_NS_From1.rdat", sep=""))

#All rats (down and not down,by modality)
masterDF_DS_Tone <- masterDF_DS_From1[masterDF_DS_From1$modality=="98", ]
masterDF_DS_Light <- masterDF_DS_From1[masterDF_DS_From1$modality=="99", ]
masterDF_NS_Light <- masterDF_NS_From1[masterDF_NS_From1$modality=="98", ]
masterDF_NS_Tone <- masterDF_NS_From1[masterDF_NS_From1$modality=="99", ]

plotBoxplotfromSessCP(experiment="All no inf Sess 1-6", masterDF=list(masterDF_DS_Tone, masterDF_NS_Light,
                                                                      masterDF_DS_Light, masterDF_NS_Tone), 
                      comp=c("S+", "S-", "S+", "S-"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=c(colindx[1], colindxB[1], colindx[2], colindxB[2]), 
                      sessFromCP = c(0, 5), lines=TRUE, yAxMinZ=-2, yAxMaxZ = 6, yAxMaxRaw = 10, WdwStart=100, 
                      WdwEnd=400, analysis = TRUE, removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=F,
                      factors = list(Cue=c("S+", "S-"), Modality=c("Tone S+", "Light S+")))


plotBoxplotfromSessCP(experiment="All no inf Sess 1-6 S+ tone", masterDF=list(masterDF_DS_Tone, masterDF_NS_Light), 
                      comp=c("S+", "S-"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=c(colindx[1], colindxB[1]), 
                      sessFromCP = c(0, 5), lines=TRUE, yAxMinZ=-2, yAxMaxZ = 6, yAxMaxRaw = 10, WdwStart=100, 
                      WdwEnd=400, analysis = TRUE, removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=F,
                      factors = list(Cue=c("S+", "S-")))

plotBoxplotfromSessCP(experiment="All no inf Sess 1-6 S+ light", masterDF=list(masterDF_DS_Light, masterDF_NS_Tone), 
                      comp=c("S+", "S-"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=c(colindx[1], colindxB[1]), 
                      sessFromCP = c(0, 5), lines=TRUE, yAxMinZ=-2, yAxMaxZ = 6, yAxMaxRaw = 10, WdwStart=100, 
                      WdwEnd=400, analysis = TRUE, removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=F,
                      factors = list(Cue=c("S+", "S-")))


corr.FRandBeh.PerUnit(experiment= "All no inf Sess 1-6", 
                      masterDF=masterDF_DS, dataProcess="Spikes",
                      WdwStart=100, WdwEnd=400, neudata=allNeuronsDS, 
                      fromCP=TRUE, bySess=TRUE, sessRange = c(-5, 5), collapsePrePost=TRUE, 
                      byTrial=FALSE, capped=TRUE, capValue=c(-90, 90), trialBinSize=15,
                      cueExcOnly=TRUE, correctOnly=FALSE,
                      graphFolder=MixedGraphFolder)


#Not down only
masterDF_DS_From1_NotDown <- masterDF_DS_From1[masterDF_DS_From1$down=="NotDown", ]
masterDF_NS_From1_NotDown <- masterDF_NS_From1[masterDF_NS_From1$down=="NotDown", ]


### See mean FR by session in order
plotBoxplotfromSessCP(experiment="Exp 1_Sess1-6", masterDF=list(masterDF_DS_From1_NotDown, masterDF_NS_From1_NotDown), 
                      comp=c("S+", "S-"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(0, 5), lines=TRUE,
                      yAxMinZ=-2, yAxMaxZ = 6, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, analysis = TRUE,
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=F,
                      factors = list(Cue=c("S+", "S-")))




#By modality
masterDF_DS_From1_Tone_NotDown <- masterDF_DS_From1[masterDF_DS_From1$down=="NotDown"  & masterDF_DS_From1$modality=="98", ]
masterDF_NS_From1_Light_NotDown <- masterDF_NS_From1[masterDF_NS_From1$down=="NotDown"  & masterDF_NS_From1$modality=="98", ]
masterDF_DS_From1_Light_NotDown <- masterDF_DS_From1[masterDF_DS_From1$down=="NotDown"  & masterDF_DS_From1$modality=="99", ]
masterDF_NS_From1_Tone_NotDown <- masterDF_NS_From1[masterDF_NS_From1$down=="NotDown"  & masterDF_NS_From1$modality=="99", ]


byBinbyUnit_Sess1_ByModality <- plotFRBoxPlotandCP(cue=c("S+", "S-", "S+", "S-"), 
                                                   experiment="Exp 1 By session BOXPLOTS NOTDOWN S+ vs S-", 
                                                   masterDF=list(masterDF_DS_From1_Tone_NotDown, masterDF_NS_From1_Light_NotDown,
                                                                 masterDF_DS_From1_Light_NotDown, masterDF_NS_From1_Tone_NotDown), 
                                                   points=F, lines=T, graphFolder=MixedGraphFolder, trialBinSize=35, 
                                                   WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                                   color=c(colindx[1], colindxB[1], colindx[2], colindxB[2]), 
                                                   legLabels=c("S+ Tone ND", "S- Light ND", "S+ Light ND", "S- Tone ND"), 
                                                   yAxMinZ=-2, yAxMaxZ=6, yAxMaxRaw=7, capped=T, capValue=c(1, 210), cueExcOnly=F, 
                                                   neudata=allNeuronsDS)


#Sess1
masterDF_DS_From1_1_Tone_NotDown <- masterDF_DS_From1[masterDF_DS_From1$down=="NotDown" & masterDF_DS_From1$sessfromCPsess==0 & masterDF_DS_From1$modality=="98", ]
masterDF_NS_From1_1_Light_NotDown <- masterDF_NS_From1[masterDF_NS_From1$down=="NotDown" & masterDF_NS_From1$sessfromCPsess==0 & masterDF_NS_From1$modality=="98", ]
masterDF_DS_From1_1_Light_NotDown <- masterDF_DS_From1[masterDF_DS_From1$down=="NotDown" & masterDF_DS_From1$sessfromCPsess==0 & masterDF_DS_From1$modality=="99", ]
masterDF_NS_From1_1_Tone_NotDown <- masterDF_NS_From1[masterDF_NS_From1$down=="NotDown" & masterDF_NS_From1$sessfromCPsess==0 & masterDF_NS_From1$modality=="99", ]


byBinbyUnit_Sess1_ByModality <- plotFRBoxPlotandCP(cue=c("S+", "S-", "S+", "S-"), 
                                                   experiment="Exp 1 Sess 1 Trial 1-35 BOXPLOTS NOTDOWN S+ vs S-", 
                                                   masterDF=list(masterDF_DS_From1_1_Tone_NotDown, masterDF_NS_From1_1_Light_NotDown,
                                                                 masterDF_DS_From1_1_Light_NotDown, masterDF_NS_From1_1_Tone_NotDown), 
                                                   points=F, lines=T, graphFolder=MixedGraphFolder, trialBinSize=35, 
                                                   WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                                   color=c(colindx[1], colindxB[1], colindx[2], colindxB[2]), 
                                                   legLabels=c("S+ Tone ND", "S- Light ND", "S+ Light ND", "S- Tone ND"), 
                                                   yAxMinZ=-2, yAxMaxZ=6, yAxMaxRaw=7, capped=T, capValue=c(1, 35), cueExcOnly=F, 
                                                   neudata=allNeuronsDS)


save(byBinbyUnit_Sess1_ByModality, file=paste(dataForRdir, "byBinbyUnit_Sess1_ByModality.rdat", sep=""))

#Sess2 (Excluding animals that had a CP on Sess 2)
masterDF_DS_From1_2_Tone_NotDown <- masterDF_DS_From1[masterDF_DS_From1$down=="NotDown" & masterDF_DS_From1$sessfromCPsess==1 & masterDF_DS_From1$modality=="98", ]
masterDF_NS_From1_2_Light_NotDown <- masterDF_NS_From1[masterDF_NS_From1$down=="NotDown" & masterDF_NS_From1$sessfromCPsess==1 & masterDF_NS_From1$modality=="98", ]
masterDF_DS_From1_2_Light_NotDown <- masterDF_DS_From1[masterDF_DS_From1$down=="NotDown" & masterDF_DS_From1$sessfromCPsess==1 & masterDF_DS_From1$modality=="99", ]
masterDF_NS_From1_2_Tone_NotDown <- masterDF_NS_From1[masterDF_NS_From1$down=="NotDown" & masterDF_NS_From1$sessfromCPsess==1 & masterDF_NS_From1$modality=="99", ]


byBinbyUnit_Sess2_ByModality <- plotFRBoxPlotandCP(cue=c("S+", "S-", "S+", "S-"), 
                                                   experiment="Exp 1 Sess 2 Trial 1-35 BOXPLOTS NOTDOWN S+ vs S-", 
                                                   masterDF=list(masterDF_DS_From1_2_Tone_NotDown, masterDF_NS_From1_2_Light_NotDown,
                                                                 masterDF_DS_From1_2_Light_NotDown, masterDF_NS_From1_2_Tone_NotDown), 
                                                   points=F, lines=T, graphFolder=MixedGraphFolder, trialBinSize=35, 
                                                   WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                                   color=c(colindx[1], colindxB[1], colindx[2], colindxB[2]), 
                                                   legLabels=c("S+ Tone ND", "S- Light ND", "S+ Light ND", "S- Tone ND"), 
                                                   yAxMinZ=-2, yAxMaxZ=6, yAxMaxRaw=7, capped=T, capValue=c(1, 35), cueExcOnly=F, 
                                                   neudata=allNeuronsDS)





byBinbyUnit <- plotFRBoxPlotandCP(cue=c("S+", "S-", "S+", "S-"), 
                                  experiment="Exp 1 Sess 2 Trial 1-35 BOXPLOTS NOTDOWN S+ vs S-", 
                                  masterDF=list(masterDF_DS_From1_2_Tone_NotDown, masterDF_NS_From1_2_Light_NotDown,
                                                masterDF_DS_From1_2_Light_NotDown, masterDF_NS_From1_2_Tone_NotDown), 
                                  points=F, lines=T, graphFolder=MixedGraphFolder, trialBinSize=35, 
                                  WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                  color=c(colindx[1], colindxB[1], colindx[2], colindxB[2]), 
                                  legLabels=c("S+ Tone ND", "S- Light ND", "S+ Light ND", "S- Tone ND"), 
                                  yAxMinZ=-2, yAxMaxZ=6, yAxMaxRaw=7, capped=T, capValue=c(1, 35), cueExcOnly=F, 
                                  neudata=allNeuronsDS)







names(byBinbyUnit) <- c("S+", "S-")

forAnalysis <- do.call("rbind", lapply(seq(1, length(byBinbyUnit)), function(x){
        byBinbyUnit[[x]]$cue <- names(byBinbyUnit)[x]
        colsForAnalysis <- c("unitIdx", "rat", "byUnitFR", "cue", "bin")
        colIdx <- colnames(byBinbyUnit[[x]]) %in% colsForAnalysis
        selForAnalysis <- byBinbyUnit[[x]][, colIdx]
        return(selForAnalysis)
})
)

ezANOVA(data=forAnalysis, wid=unitIdx, within = c("bin", "cue"), dv=byUnitFR)

wilcox.test(x=forAnalysis$byUnitFR[forAnalysis$cue=="S+"],
            y=forAnalysis$byUnitFR[forAnalysis$cue=="S-"],
            paired=T)
# data:  first 5 trials of every animal, FR on S+ trials vs. FR on S- trials 
# W = 2739, p-value = 6.932e-06
# alternative hypothesis: true location shift is not equal to 0





#ENTRY
plotFRandCP_Entry(experiment="Exp 3", cue="S+", wdwLabel="Post S+ Entry", masterDF=list(masterDF_EntryDS), graphFolder=MixedGraphFolder, 
                  trialBinSize=30, WdwStart =0, WdwEnd=1000, dataProcess="Zscores", colindx =colindx[1], 
                  capped=T, capValue = c(-90, 90), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS)

plotFRandCP_Entry(experiment="Exp 3", cue="S-", wdwLabel="Post S- Entry", masterDF=list(masterDF_EntryNS), graphFolder=MixedGraphFolder, 
                  trialBinSize=30, WdwStart =0, WdwEnd=1000, dataProcess="Zscores", colindx =colindx[1], 
                  capped=T, capValue = c(-90, 90), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryNS)


plotFRandCP_Entry(experiment="Exp 3", cue=c("S+", "S-"), wdwLabel="Post S+ and S- Entry", masterDF=list(masterDF_EntryDS, masterDF_EntryNS), graphFolder=MixedGraphFolder, 
                  trialBinSize=30, WdwStart =0, WdwEnd=1000, dataProcess="Zscores", colindx =c(colindx[1], "darkblue"), 
                  capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 3, legLabels=c("S+", "S-"),
                  yAxMaxRaw = 10, neudata=allNeuronsEntryDS)

plotFRandCP_Entry(experiment="Exp 3", cue=c("S+", "S-"), wdwLabel="Pre S+ and S- Entry", masterDF=list(masterDF_EntryDS, masterDF_EntryNS), graphFolder=MixedGraphFolder, 
                  trialBinSize=30, WdwStart =-2000, WdwEnd=0, dataProcess="Zscores", colindx =c(colindx[1], "darkblue"), 
                  capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 1.5, legLabels=c("S+", "S-"),
                  yAxMaxRaw = 10, neudata=allNeuronsEntryDS)

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

#Post cue
PostCueFR_DSvsNS_fromCP <- compareDSvsNSfromCP(masterDF=list(masterDF_DS, masterDF_NS), 
                                               trialBinSize=15, event="cue", correctOnly=F, 
                                               WdwStart=100, WdwEnd=400, capped=T, capValue=c(-90, 90), 
                                               dataProcess="Zscores")
PostCueFR_DSvsNS_fromCP$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP$p, method="holm") 
PostCueFR_DSvsNS_fromCP$sig <- giveStars(PostCueFR_DSvsNS_fromCP$p.adj)

save(PostCueFR_DSvsNS_fromCP, file=paste(dataForRdir, "PostCueFR_DSvsNS_fromCP.rdat"))

#Post Cue Not Down only
PostCueFR_DSvsNS_fromCP_NotDown <- compareDSvsNSfromCP(masterDF=list(masterDF_DS_NotDown, masterDF_NS_NotDown), 
                                                       trialBinSize=35, event="cue", correctOnly=F, WdwStart=100, 
                                                       WdwEnd=400, capped=T, capValue=c(-140, 140), 
                                                       dataProcess="Zscores")

PostCueFR_DSvsNS_fromCP_NotDown$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP_NotDown$p, method="holm") 
PostCueFR_DSvsNS_fromCP_NotDown$sig <- giveStars(PostCueFR_DSvsNS_fromCP_NotDown$p.adj)
save(PostCueFR_DSvsNS_fromCP_NotDown, file=paste(dataForRdir, "PostCueFR_DSvsNS_fromCP_NotDown.rdat"))
# bin     V            p   n        p.adj sig
# V   -120 to -105   469 2.191700e-04  33 6.575100e-04 ***
#         V1   -105 to -90  2072 3.689263e-05  73 1.844632e-04 ***
#         V2    -90 to -75  3888 4.931851e-12  92 6.411407e-11 ***
#         V3    -75 to -60  2271 1.037938e-10  69 1.141732e-09 ***
#         V4    -60 to -45 12043 8.428291e-17 166 1.348526e-15 ***
#         V5    -45 to -30  5468 1.928989e-12 111 2.700584e-11 ***
#         V6    -30 to -15  8542 8.552519e-13 142 1.282878e-11 ***
#         V7      -15 to 0  5661 9.388451e-10 117 9.388451e-09 ***
#         V8       0 to 15  4317 8.427414e-09 102 6.741931e-08 ***
#         V9      15 to 30  4318 1.010848e-11  98 1.213018e-10 ***
#         V10     30 to 45  2072 6.127815e-07  70 4.289471e-06 ***
#         V11     45 to 60  2332 5.652191e-09  72 5.086972e-08 ***
#         V12     60 to 75   864 2.074738e-06  44 1.244843e-05 ***
#         V13     75 to 90   392 3.042044e-04  30 6.575100e-04 ***
#         V14    90 to 105   715 5.785305e-05  41 2.314122e-04 ***
#         V15   105 to 120    88 6.103516e-04  13 6.575100e-04 ***


#In 35 trial bins:
# bin     V             p               n        p.adj sig
# -140 to -105  2974    1.019527e-06    86      3.058582e-06 ***
# -105 to -70   4774    4.809712e-14    101     2.404856e-13 ***
# -70 to -35    14665   1.115045e-22    178     8.920357e-22 ***
# -35 to 0      14313   1.064584e-16    183     7.452088e-16 ***
# 0 to 35       8954    8.680977e-15    143     5.208586e-14 ***
# 35 to 70      4453    1.708589e-11    100     6.834354e-11 ***
# 70 to 105     784     4.299764e-05    43      8.599528e-05 ***
# 105 to 140    297     4.401803e-05    25      8.599528e-05 ***


#Pre entry
PreEntryFR_DSvsNS_fromCP <- compareDSvsNSfromCP(masterDF=list(masterDF_EntryDS, masterDF_EntryNS), trialBinSize=30, event="entry", correctOnly=F, WdwStart=-2000, WdwEnd=0, capped=T, capValue=c(-120, 120), dataProcess="Zscores")
PreEntryFR_DSvsNS_fromCP$p.adj <- p.adjust(PreEntryFR_DSvsNS_fromCP$p, method="holm") 
PreEntryFR_DSvsNS_fromCP$sig <- giveStars(PreEntryFR_DSvsNS_fromCP$p.adj)
save(PreEntryFR_DSvsNS_fromCP, file=paste(dataForRdir, "PreEntryFR_DSvsNS_fromCP.rdat"))


#Post entry
PostEntryFR_DSvsNS_fromCP <- compareDSvsNSfromCP(masterDF=list(masterDF_EntryDS, masterDF_EntryNS), trialBinSize=30, event="entry", correctOnly=F, WdwStart=0, WdwEnd=1000, capped=T, capValue=c(-120, 120), dataProcess="Zscores")
PostEntryFR_DSvsNS_fromCP$p.adj <- p.adjust(PostEntryFR_DSvsNS_fromCP$p, method="holm") 
PostEntryFR_DSvsNS_fromCP$sig <- giveStars(PostEntryFR_DSvsNS_fromCP$p.adj)
save(PostEntryFR_DSvsNS_fromCP, file=paste(dataForRdir, "PostEntryFR_DSvsNS_fromCP.rdat"))





### PLOT FR IN THE POSTCUE WINDOW AND TAIL PERIOD BY RAT BY TRIAL BY NEURON
PLOT_perRatTrial_FRandBEH(masterDFsumm=masterDFsumm, wdw="PostCue", threshold=2, MixedGraphFolder=MixedGraphFolder)
PLOT_perRatTrial_FRandBEH(masterDFsumm=masterDFsumm, wdw="PostSpike", threshold=2, MixedGraphFolder=MixedGraphFolder)


### PLOT FR AROUND THE CUE (PTSH) as a function of distance to change point
plotFRandCPhistogram(experiment="Exp 3", masterDF=list(masterDF_DS, masterDF_NS), graphFolder=MixedGraphFolder, 
                     trialBinSize=15, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
                     capped=T, capValue = c(-90, 90), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 3, 
                     psthmin=0.5, psthmax=2, imgFormat="pdf", neudata=allNeuronsDS)


### PLOT FR AROUND ENTRIES (PTSH) as a function of distance to change point
plotFRandCPhistogram_Entry(experiment="Exp 3", masterDF=list(masterDF_EntryDS, masterDF_EntryNS), graphFolder=MixedGraphFolder, 
                           trialBinSize=30, dataProcess="Zscores", correctOnly=FALSE, color=c(colindx[1], "darkblue"), 
                           capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 2, yAxMaxRaw = 3, 
                           WdwStart=-1500, WdwEnd=1500, imgFormat="pdf", neudata=allNeuronsEntryDS)


#Same thing but by session from CP instead of trial
plotPSTHfromSessCP(experiment="Exp 3", masterDF=list(masterDF_DS_Tone, masterDF_DS_Light), 
                   graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                   correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-2, 2), 
                   yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, psthmin=-0.5, psthmax=2, 
                   imgFormat="pdf", neudata=allNeuronsDS)



### PLOT BOXPLOTS AROUND CUE BY SESS FROM CP AND GROUP

### THIS ONE IS THE SAME BUT GROUPING THE SESSIONS BEFORE AND AFTER THE CP

#Light vs. tone in DOWN animals
plotBoxplotPrePostCP(experiment="EXP 3 DOWN", masterDF=list(masterDF_DS_Tone_Down, 
                                                            masterDF_DS_Light_Down), 
                     factors=list(Modality=c("Tone", "Light")), removeCPsess = T,
                     comp=c("Tone S+ Down", "Light S+ Down"), ANOVA = FALSE,
                     graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, cueExcOnly=FALSE, 
                     color=colindx, sessFromCP = c(-5, 5), yAxMinZ = -2, yAxMaxZ = 6, yAxMaxRaw = 10, 
                     WdwStart=100, WdwEnd=400, removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)

#                       Comparison    W        p.val        p.adj       sig
# W     Pre CP Tone vs. Pre CP Light 1082 1.359494e-01 2.718987e-01    
# W1  Post CP Tone vs. Post CP Light  140 5.278148e-01 5.278148e-01    
# W2    Pre CP Tone vs. Post CP Tone  199 1.730899e-07 6.923595e-07     ***
# W11 Pre CP Light vs. Post CP Light   85 8.237366e-03 2.471210e-02       *



#Light vs. tone in NOTDOWN ANIMALS
plotBoxplotPrePostCP(experiment="EXP 1 NOTDOWN", masterDF=list(masterDF_DS_Tone_NotDown, 
                                                               masterDF_DS_Light_NotDown), 
                     factors=list(Modality=c("Tone", "Light")), removeCPsess = T,
                     comp=c("Tone S+ Down", "Light S+ Down"), ANOVA=FALSE,
                     graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, cueExcOnly=FALSE, 
                     color=colindx, sessFromCP = c(-5, 5), yAxMinZ = -2, yAxMaxZ = 6, yAxMaxRaw = 10, 
                     WdwStart=100, WdwEnd=400, removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)
# Comparison    W        p.val       p.adj sig
# W     Pre CP Tone vs. Pre CP Light 3106 0.0031974634 0.009592390  **
# W1  Post CP Tone vs. Post CP Light 1357 0.4988188096 0.498818810    
# W2    Pre CP Tone vs. Post CP Tone  340 0.0061353974 0.012270795   *
# W11 Pre CP Light vs. Post CP Light 8392 0.0004210872 0.001684349  **



### THIS ONE IS THE SAME BUT GROUPING THE SESSIONS BEFORE AND AFTER THE CP
plotBoxplotPrePostCP(experiment="No infus. by Mod and Down", masterDF=list(masterDF_DS_Tone_Down, 
                                                                           masterDF_DS_Light_Down, masterDF_DS_Tone_NotDown, masterDF_DS_Light_NotDown), 
                     factors=list(Modality=c("Tone", "Light"), Down=c("Down", "NotDown")), removeCPsess = T,
                     comp=c("Tone S+ Down", "Light S+ Down", "Tone S+ Not Down", "Light S+ Not Down"), 
                     graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, cueExcOnly=FALSE, 
                     color=c(colindx, colindxB), sessFromCP = c(-4, 5), yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, 
                     WdwStart=100, WdwEnd=400, removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)


plotBoxplotPrePostCP(experiment="No infus. Down vs Not Down TAILshortsc", masterDF=list(masterDF_DS_Down, masterDF_DS_NotDown), 
                     factors=list(Down=c("Down", "NotDown")), removeCPsess = TRUE,
                     comp=c("Down", "Not Down"), graphFolder=MixedGraphFolder, 
                     dataProcess="Zscores", correctOnly=FALSE, cueExcOnly=TRUE, 
                     color=c(colindx[1], colindxB[1]), sessFromCP = c(-5, 5), yAxMinZ = -2, 
                     yAxMaxZ = 4, yAxMaxRaw = 10, WdwStart=750, WdwEnd=2000, removeOutliers=F, 
                     imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)

### 100-400ms after S+, cue excited units only
#                                   Comparison          W        p.val              p.adj sig
#W              Pre CP Down vs. Pre CP NotDown          1955    6.595299e-02 1.319060e-01    
#W1             Post CP Down vs. Post CP NotDown        1841    1.326034e-01 1.326034e-01    
#W2             Pre CP Down vs. Post CP Down            189     4.789009e-10 1.436703e-09 ***
#W11            Pre CP NotDown vs. Post CP NotDown      2437    1.873476e-10 7.493906e-10 ***


### 
#                             Comparison    W       p.val       p.adj sig
# W       Pre CP Down vs. Pre CP NotDown 1967 0.072199149 0.216597448    
# W1    Post CP Down vs. Post CP NotDown 1885 0.087508248 0.216597448    
# W2        Pre CP Down vs. Post CP Down  459 0.001314881 0.005259522  **
# W11 Pre CP NotDown vs. Post CP NotDown 4674 0.173321691 0.216597448    

#################################################333
### BOXPLOTS BY CHUNK OF TRIALS FROM CHANGE POINT
###################################################

NotDown35bin <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 Not Down", masterDF=list(masterDF_DS_NotDown, masterDF_NS_NotDown), 
                                   graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="Zscores", 
                                   correctOnly=FALSE, points=F, lines=F, color=c(colindx[1], colindxB[1]), legLabels=c("S+", "S-"), 
                                   yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, capValue=c(-140, 140), cueExcOnly=F, neudata=allNeuronsDS)

uniqueBin <- unique(NotDown35bin[[1]]$bin)

Diff_DSNS_NotDown_35bins <- do.call("rbind", lapply(seq(1, length(uniqueBin)), function(i){
        DSbinSel <- NotDown35bin[[1]][NotDown35bin[[1]]$bin==uniqueBin[i], ]
        NSbinSel <- NotDown35bin[[2]][NotDown35bin[[2]]$bin==uniqueBin[i], ]
        
        matchUnitIdx <- match(DSbinSel$unitIdx, NSbinSel$unitIdx)
        NSfiring <- NSbinSel$byUnitFR[matchUnitIdx]
        
        DSbinSel$byUnitFR_NS <- NSfiring
        DSbinSel$DiffFR <- DSbinSel$byUnitFR-NSfiring
        
        DSbinSel
})
)


DiffFR_ByBin(data=Diff_DSNS_NotDown_35bins, cueExcOnly=TRUE, color=colindx[1], ymin=-2, ymax=18,
             graphFolder=MixedGraphFolder, experiment="Exp 1 Not Down", points=TRUE, 
             comparisons=list(c(2, 4), c(4, 6)))
# [[1]]
# idx     bin                        comp    W        p.val        p.adj sig
# W     1 1 vs. 3 -140 to -105 vs. -70 to -35  455 7.360720e-03 2.944288e-02   *
#         W3    2 2 vs. 4    -105 to -70 vs. -35 to 0 1135 1.681229e-01 5.043687e-01    
# W1    3 3 vs. 5      -70 to -35 vs. 0 to 35 1523 8.767221e-05 4.383610e-04 ***
#         W11   4 4 vs. 6       -35 to 0 vs. 35 to 70 1071 2.933992e-06 1.760395e-05 ***
#         W2    5 5 vs. 7       0 to 35 vs. 70 to 105  567 6.055316e-01 9.035712e-01    
# W21   6 6 vs. 8     35 to 70 vs. 105 to 140  166 4.517856e-01 9.035712e-01    
# 
# [[2]]
# bins               comparison    W        p.val        p.adj sig
# W  2 vs. 4 -105 to -70 vs. -35 to 0 1135 1.681229e-01 1.681229e-01    
# W1 4 vs. 6    -35 to 0 vs. 35 to 70 1071 2.933992e-06 5.867984e-06 ***       





#This one compares S+ firing rate among bins only (and make a boxplot of it)
compareBins(data=Diff_DSNS_NotDown_35bins, cueExcOnly=TRUE, color=colindx[1], ymin=-2, ymax=8, 
            graphFolder=MixedGraphFolder, experiment="Exp 1 Not Down", points=F, 
            comparisons=list(c(1, 4), c(4, 8)), 
            cue="S+")
# [[1]]
# bin           n                        comp           W        p.val        p.adj sig
# 1 vs. 3       20 vs. 71 -140 to -105 vs. -70 to -35  429 3.590332e-03 0.0143613260   *
# 2 vs. 4       31 vs. 83    -105 to -70 vs. -35 to 0 1063 7.778560e-02 0.2333567972    
# 3 vs. 5       71 vs. 75      -70 to -35 vs. 0 to 35 1864 8.906571e-04 0.0044532857  **
# 4 4 vs. 6     83 vs. 54       -35 to 0 vs. 35 to 70 1391 9.135982e-05 0.0005481589 ***
# 5 5 vs. 7     75 vs. 23       0 to 35 vs. 70 to 105  858 4.866258e-01 0.6328163766    
# 6 6 vs. 8     54 vs. 14     35 to 70 vs. 105 to 140  346 3.164082e-01 0.6328163766    
# 
# [[2]]
# bins         n                comparison      W       p.val      p.adj sig
# 1 vs. 4 20 vs. 83 -140 to -105 vs. -35 to 0   545 0.008847681 0.01243276   *
# 4 vs. 8 83 vs. 14   -35 to 0 vs. 105 to 140   337 0.006216380 0.01243276   *




############################
### SCATTERPLOTS BY SESSION
############################

#This function generates a heatmap with firing rate per unit per session in the events designed by "data" (adjust ZcolLabels to designate those events)
#It also generates a data frame ("toplot") with info per unit per session

toplot <- megaplot(data=list(allNeuronsDS, allNeuronsNS), CPdata=CPdata, csacqidx=csacqidx, BLwdw=2, winmin=100, winmax=400,
                   colpalette="Rainbow", minFR=-2, maxFR=6, graphFolder=MixedGraphFolder, dataForRdir=dataForRdir,
                   ZcolLabels=c("ZDS", "ZNS"), arrangeBy="ZDS")


## SCATTERPLOT WHERE EACH SESSION IS A DOT AND THE AXES ARE --> X axis: performance index; Y axis: mean FR. 

##I'm going to choose the CUE EXCITED UNITS ONLY before summarizing the data by session
toplot <- toplot[toplot$DSExc==TRUE, ]

toplot$uniqSession <- paste(toplot$rat, toplot$expt)

toplotSumm <- do.call("rbind", lapply(seq(1, length(unique(toplot$uniqSession))), function(x){
        selSess <- toplot[toplot$uniqSession==unique(toplot$uniqSession)[x], ]
        meanZDS <- mean(selSess$ZDS, na.rm=T)
        meanZNS <- mean(selSess$ZNS, na.rm=T)
        sessSumm <- selSess[1, ]
        sessSumm$ZDS <- meanZDS
        sessSumm$ZNS <- meanZNS
        sessSumm$units <- nrow(selSess)
        return(sessSumm)
})
)

toplotSumm <- as.data.frame(toplotSumm)
CPsessIndex <- toplotSumm$BeforeCP==0


#Divide by condition down x version
toplotSumm_DownTone <- toplotSumm[toplotSumm$down=="Down" & toplotSumm$version=="98", ]
toplotSumm_DownLight <- toplotSumm[toplotSumm$down=="Down" & toplotSumm$version=="99", ]
toplotSumm_NotDownTone <- toplotSumm[toplotSumm$down=="NotDown" & toplotSumm$version=="98", ]
toplotSumm_NotDownLight <- toplotSumm[toplotSumm$down=="NotDown" & toplotSumm$version=="99", ]

#Divide by condition (down vs not down)
toplotSumm_Down <- toplotSumm[toplotSumm$down=="Down", ]
toplotSumm_NotDown <- toplotSumm[toplotSumm$down=="NotDown", ]

#Divide by condition (S+ light vs S+ tone)
toplotSumm_Tone <- toplotSumm[toplotSumm$version=="98", ]
toplotSumm_Light <- toplotSumm[toplotSumm$version=="99", ]


#DOWN: FR vs Perf Index, Latency and RR (by session)
FRandPerf_Scatterplot(condition="Down", data=toplotSumm_Down, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="PerfIndex")


FRandPerf_Scatterplot(condition="Down_Tone", data=toplotSumm_DownTone, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="PerfIndex")
FRandPerf_Scatterplot(condition="Down_Light", data=toplotSumm_DownLight, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="PerfIndex")


#NOT DOWN: FR vs Performance index, latency and response ratio (by session)

######## BEFORE CP, cue-exc only
FRandPerf_Scatterplot(condition="NOT DOWN Cue exc only Pre CP", data=toplotSumm_NotDown, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="PerfIndex", chunk=c(-1))
FRandPerf_Scatterplot(condition="NOT DOWN Cue exc only Pre CP", data=toplotSumm_NotDown, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="Latency", chunk=c(-1))
FRandPerf_Scatterplot(condition="NOT DOWN Cue exc only Pre CP", data=toplotSumm_NotDown, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="RespRatio", chunk=c(-1))


######### AFTER CP, cue-exc only
FRandPerf_Scatterplot(condition="NOT DOWN Cue exc only Post CP", data=toplotSumm_NotDown, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="PerfIndex", chunk=c(1))
FRandPerf_Scatterplot(condition="NOT DOWN Cue exc only Post CP", data=toplotSumm_NotDown, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="Latency", chunk=c(1))
FRandPerf_Scatterplot(condition="NOT DOWN Cue exc only Post CP", data=toplotSumm_NotDown, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="RespRatio", chunk=c(1))


######### NOT DOWN, ALL SESSIONS, cue-exc only
FRandPerf_Scatterplot(condition="NOT DOWN Cue exc", data=toplotSumm_NotDown, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="PerfIndex")
FRandPerf_Scatterplot(condition="NOT DOWN Cue exc", data=toplotSumm_NotDown, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="Latency")
FRandPerf_Scatterplot(condition="NOT DOWN Cue exc", data=toplotSumm_NotDown, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="RespRatio")


######### DOWN, ALL SESSIONS, cue-exc only
FRandPerf_Scatterplot(condition="DOWN Cue exc", data=toplotSumm_Down, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="PerfIndex")
FRandPerf_Scatterplot(condition="DOWN Cue exc", data=toplotSumm_Down, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="Latency")
FRandPerf_Scatterplot(condition="DOWN Cue exc", data=toplotSumm_Down, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="RespRatio")



FRandPerf_Scatterplot(condition="Tone", data=toplotSumm_NotDown, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="PerfIndex")

FRandPerf_Scatterplot(condition="Light", data=toplotSumm_Light, graphFolder=ScatterplotFolder, 
                      CPtoo=TRUE, CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="PerfIndex")


#############################

# FOCUS ON THE UNITS RECORDED ON THE DAY BEHAVIOR CHANGED
UnitHeatMap(data=masterDF, sessFromCP=4, FRparameters=allNeuronsDS$parameters, 
            folder=BySessFolder, winmin=100, winmax=300, BLmin=-2000, BLmax=0)



# Megaplot with firing in the 400ms window after 4 events: S+ responded to, S+ missed, S- responded to, S- missed
# This function makes the plot and also saves a data frame that will be useful for the next graph
megaplot(data=list(allNeuronsDSresponded, allNeuronsDSmissed, allNeuronsNSresponded, allNeuronsNSmissed), CPdata=CPdata, csacqidx=csacqidx, BLwdw=2, winmin=0, winmax=400,
         colpalette="Rainbow", minFR=-1, maxFR=3, graphFolder=MixedGraphFolder, dataForRdir=dataForRdir, 
         ZcolLabels=c("ZDSresp", "ZDSmissed", "ZNSresp", "ZNSmissed") )

megaplot(data=list(allNeuronsDS, allNeuronsNS), CPdata=CPdata, csacqidx=csacqidx, BLwdw=2, winmin=0, winmax=400,
         colpalette="Rainbow", minFR=-2, maxFR=6, graphFolder=MixedGraphFolder, dataForRdir=dataForRdir,
         ZcolLabels=c("ZDS", "ZNS"), arrangeBy="ZDS")





#toplotSummNoCPsess <- toplotSumm[CPsessIndex==F, ]



plot.new()
plot.window(xlim=c(-3, 6), ylim=c(-2, 10))

sapply(c(-1, 0, 1), function(x){
        dataSel <- toplotSumm[toplotSumm$BeforeCP==x,]
        if(x==1){pch=19; col=colindx[1]}
        if(x==0){pch=19; col="purple"}
        if(x==-1){pch=19; col=colindx[2]}
        perf <- as.numeric(as.character(dataSel$CSplusTA))
        FR <- as.numeric(as.character(dataSel$ZDS))
        
        points(x=perf, y=FR, pch=pch, cex=2, col=col)
        
        axis(side=1, at=seq(-3, 6, 1), cex.axis=1.4)
        axis(side=2, at=seq(-2, 10, 1), cex.axis=1.4, las=2)
        
        abline(h=0, lty=2)
        abline(v=0)
        
        mtext(side=2, line=2.5, text="Firing rate (Z sc.)", cex=1.5, font=2)
        mtext(side=1, line=2.5, text="Performance index", cex=1.5, font=2)
        
})

perfAll <- as.numeric(as.character(toplotSumm$CSplusTA))
FRAll <- as.numeric(as.character(toplotSumm$ZDS))

fit <- lm(FRAll ~ perfAll)
Rsq <- summary(fit)$r.squared
coeff <- summary(fit)$coefficients[2,1]
p.val <- summary(fit)$coefficients[2,4]

text(x=-2.5, y=6, labels = "With outliers", font=2, col="gray50")
text(x=-2.5, y=5.5, labels = paste("Coeff. = ", round(coeff, 2), sep=""))
text(x=-2.5, y=5, labels = paste("p = ", round(p.val, 5), sep=""))
text(x=-2.5, y=4.5, labels = paste("R.sq = ", round(Rsq, 5), sep=""))

abline(a=summary(fit)$coefficients[1, 1], b=summary(fit)$coefficients[2,1], lwd=2, col="gray50")

#Find outliers
cooksd <- cooks.distance(fit)

#Plot outlier test
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line. The cutoff point is 4 times the mean of cook´s distances
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels

outliers <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))])
points(x=perfAll[outliers]+0.2, y=FRAll[outliers]+0.1, pch=8)

fit <- lm(FRAll[-outliers] ~ perfAll[-outliers])
Rsq <- summary(fit)$r.squared
coeff <- summary(fit)$coefficients[2,1]
p.val <- summary(fit)$coefficients[2,4]

text(x=-0.5, y=6, labels = "Without outliers", font=2, col="black")
text(x=-0.5, y=5.5, labels = paste("Coeff. = ", round(coeff, 2), sep=""))
text(x=-0.5, y=5, labels = paste("p = ", round(p.val, 5), sep=""))
text(x=-0.5, y=4.5, labels = paste("R.sq = ", round(Rsq, 5), sep=""))


abline(a=summary(fit)$coefficients[1, 1], b=summary(fit)$coefficients[2,1], lwd=2)

legend(x=4.5, y=0.5, col=colindx, pch=19, legend=c("After CP", "Before CP"))



# Plot scatterplot: x axis:Perf Idx, y axis=FR, dots=per rat per trial mean FR (pool all neurons from that session) before CP vs after CP

# DSdat and NSdat are created with the function compareDSvsNSfromCP and saved in dataForRDir. If I need any other structure of data for DSdat and NSdat (e.g. entry-related data, run compareDSvsNSfromCP with the appropriate parameters)
load(file=paste(dataForRdir, "DSdat.rdat", sep=""))
load(file=paste(dataForRdir, "NSdat.rdat", sep=""))

#Function to divide sessions into 3 chunks: before CP, CP, after CP

DFbyChunk <- ByChunkByRatByTrial(data=masterDF_DS, WdwStart=0, WdwEnd=400, eventType="cue", cue="S+", parameters=allNeuronsDS$parameters)


makeScatterplot <- function(data=DFbyChunk[[1]]){
        
        plot.new()
        plot.window(xlim=c(-10, 10), ylim=c(-5, 15))
        
        points(x=data$CueSpecif, y=data$meanFR_Zsc, pch=19, cex=0.5, col="blue")
        
        axis(side=1, cex.axis=1.4)
        axis(side=2, cex.axis=1.4, las=2)
        
        fit <- summary(lm(data$meanFR_Zsc ~ data$CueSpecif))
        
        abline(a=fit$coefficients[1, 1], b=fit$coefficients[2, 1], col="darkblue")
}




# % OF BINS EXCITED/INHIBITED

#Matrix in which rows are 50ms bins after the cue, columns are individual neurons and the values indicate if the neuron was EXCITED (ExcBins) or INHIBITED (InhBins) on that bin
ExcBins <- KC.sigbins(path=NEXfiles, startt=0, endt=10000, event=1, BLwdw=2, PostEvent_wdw=3, pbin=0.05, funcdirect=funcdirect)
InhBins <- KC.inhib.sigbins(path=NEXfiles, startt=0, endt=10000, event=1, BLwdw=2, PostEvent_wdw=3, pbin=0.05, funcdirect=funcdirect)

## 

#########################################
### % EXCITED BY BIN IN NOT-DOWN RATS ###
#########################################

#Units are considered cue excited if their FR on 2 consecutive 50ms windows exceeds the 99.9% upper limit of a Poisson distr. calculated using 2s baseline
PercResults <- plotFRandCP(cue=c("S+"), experiment="Exp 1 Not Down S+ PERCEXC", masterDF=list(masterDF_DS_NotDown), 
                           legLabels=c("S+"), graphFolder=MixedGraphFolder, trialBinSize=35, 
                           WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), 
                           capped=T, capValue = c(-140, 140), yAxMinZ = -1, yAxMaxZ = 6, yAxMaxRaw = 10, 
                           neudata=c(allNeuronsDS), cueExcOnly = F)
# bin trialBins CueEx notCueExc CueInh notCueInh
# 1   1      -140    10        76      5        81
# 2   2      -105     8        93      3        98
# 3   3       -70    28       150      1       177
# 4   4       -35    41       142      7       176
# 5   5         0    41       102      6       137
# 6   6        35    38        62      3        97
# 7   7        70    17        26      0        43
# 8   8       105     9        16      0        25


#Is the proportion of cue-exc neurons modulated by the amount of training? (run a different chi.sq test for bins 1, 3 and 5 and for bins 2, 4 and 6, that way the units on each cell will be independent for sure)
#For Chi-sq analysis. I need to make sure that each cell has a completely different population of neurons. Bc each session has 35trials of each kind (in the NOT DOWN group), by only comparing bins that are separated by 35 trials I can guarantee that.
contingency_table_comp1 <- t(PercResults[[1]][is.even(1:nrow(PercResults[[1]])), c(2, 3, 4)])
contingency_table_comp2 <- t(PercResults[[1]][is.even(1:nrow(PercResults[[1]]))==FALSE, c(2, 3, 4)])

inhCont_table_comp1 <- t(PercResults[[1]][is.even(1:nrow(PercResults[[1]])), c(2, 5, 6)])
inhCont_table_comp2 <- t(PercResults[[1]][is.even(1:nrow(PercResults[[1]]))==FALSE, c(2, 5, 6)])

#Chi-sq analysis excitations:
#Bins -105 to -70, -35 to 0, 35 to 70 and 105 to 140
critData <- contingency_table_comp1[-1, ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 27.649, df = 3, p-value = 4.304e-06

#Bins -140 to -105, -70 to -35, 0 to 35 and 70 to 105
critData <- contingency_table_comp2[-1, ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 21.234, df = 3, p-value = 9.415e-05
# 

p.adjust(p=c(4.304e-06, 9.415e-05), method="holm")
# [1] 8.608e-06 9.415e-05

#Post-hoc comparisons for % cue-excited units
contTableExc <- PercResults[[1]][ , c(3, 4)]

#Bin -140 to -105 vs. -70- -35
fisher.test(contTableExc[c(1, 3), ], alternative="less") #95% CI:  0.000000 1.421061 Odds ratio=0.7057797 p=0.2438
#Bin -70 to -35 vs. 0-35
fisher.test(contTableExc[c(3, 5), ], alternative="less") #95% CI:0.0000000 0.7584871 odds ratio=0.4655249 p=0.003864
#Bin 0-35 vs. 70-105
fisher.test(contTableExc[c(5, 7), ], alternative="less") #95% CI:0.0000000 1.197249 odds ratio=0.6164848 p=0.1236


#Bin -105 to -70 vs. -35 to 0
fisher.test((contTableExc[c(2, 4), ]), alternative="less") #95% CI: 0.0000000 0.61091 odds ratio=0.2990564 p= 0.001131
#Bin -35 to 0 vs. 35 to 70
fisher.test((contTableExc[c(4, 6), ]), alternative="less") #95% CI:  0.000000  0.7647095 odds ratio=0.472416 p=0.004243
#Bin 35 to 70 vs. 105 to 140
fisher.test((contTableExc[c(6, 8), ]), alternative="less") #95% CI:  0.000000  2.639951 odds ratio=1.088891 p=0.6572



p.adjust(p=c(0.2438, 0.003864, 0.1236, 0.001131, 0.004243, 0.6572), method="holm")
# 0.487600 0.019320 0.370800 0.006786 0.019320 0.657200

#Chi-sq for inhibitions:
#Bins -105 to -70, -35 to 0, 35 to 70 and 105 to 140
critData <- inhCont_table_comp1[-1, ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 1.0939, df = 3, p-value = 0.7785

#Bins -140 to -105, -70 to -35, 0 to 35 and 70 to 105
critData <- inhCont_table_comp2[-1, ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 8.7867, df = 3, p-value = 0.03227
#########

p.adjust(p=c(0.7785, 0.03227), method="holm")
#[1] 0.77850 0.06454

##########################
## BOXPLOTS BY SESSION
##########################

plotBoxplotfromSessCP(experiment="Exp 1_Sess1-6 by Modality w Lines", 
                      masterDF=list(masterDF_DS_From1_Tone_NotDown, masterDF_NS_From1_Light_NotDown,
                                    masterDF_DS_From1_Light_NotDown, masterDF_NS_From1_Tone_NotDown), 
                      comp=c("S+", "S-", "S+", "S-"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=c(colindx[1], colindxB[1], colindx[2], colindxB[2]), 
                      sessFromCP = c(0, 5), lines=T, analysis=F,
                      yAxMinZ=-2, yAxMaxZ = 6, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=F,
                      factors=list(Modality=c("Tone S+", "Light S+"), Cue=c("S+", "S-")))



### SESS FROM CP DATA

sessMinus1and1 <- sessFromCPdata(data=masterDF_DS, sessFromCP=c(-1, 1), FRparameters=allNeuronsDS$parameters, 
                                 winmin=0, winmax=400, BLmin=-2000, BLmax=0)

allSess_DS <- sessFromCPdata(data=masterDF_DS, sessFromCP=c(-4:4), FRparameters=allNeuronsDS$parameters, 
                             winmin=0, winmax=400, BLmin=-2000, BLmax=0)

allSess_NS <- sessFromCPdata(data=masterDF_NS, sessFromCP=c(-4:4), FRparameters=allNeuronsDS$parameters, 
                             winmin=0, winmax=400, BLmin=-2000, BLmax=0)

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




