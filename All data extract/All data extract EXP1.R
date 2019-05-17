
#############################################################
### EXPERIMENT 1: NO INFUSIONS, SWITCH SHORT TO LONG ITI  ###
#############################################################

### LOAD IMPORTANT LIBRARIES
install.packages("matrixStats")
install.packages('ez')
intall.packages('R.utils')

library(matrixStats)
library(ez)
library(R.utils)

setwd("E:/Dropbox (EinsteinMed)/NMDA/")

#Define colors
colindx <- c("#2171b5", "#cb181d") #Strong blue and red
colindxB <- c("#bdd7e7", "#fcae91") #Less strong blue and red
colindxC <- c("#eff3ff", "#fb6a4a") #Even less strong blue and red
colindxD <- c("#6baed6", "#fee5d9") #Lightest blue and red


### DEFINE ALL THE IMPORTANT FOLDERS

funcdirect <- paste(getwd(), "/R functions/", sep="")
datafolder <- paste(getwd(), "/EXP1_NAc FR acquisition/MedPC files/", sep="")
dataForRdir <- paste(getwd(), "/EXP1_NAc FR acquisition/Data for R/", sep="")
dataForRCumulative <- paste(getwd(), "/EXP1_NAc FR acquisition/Data for R cumulative/", sep="")
behGraphFolder <- paste(getwd(), "/EXP1_NAc FR acquisition/Graphs/Behavior/", sep="")
CPfuncFolder <- paste(getwd(), '/Change_Point-master/', sep="")
CPGraphFolder <- paste(getwd(), "/EXP1_NAc FR acquisition/Graphs/Behavior/Change point/", sep="")
rasterGraphFolder <- paste(behGraphFolder, "/Rasters/", sep="")
MixedGraphFolder <- paste(getwd(), "/EXP1_NAc FR acquisition/Graphs/Mixed/", sep="")
neuGraphFolder <- paste(getwd(), "/EXP1_NAc FR acquisition/Graphs/Neuronal/", sep="")
otherGraphFolder <- paste(getwd(), "/EXP1_NAc FR acquisition/Graphs/Others/", sep="")
NEXfiles <- paste(getwd(), "/EXP1_NAc FR acquisition/NEX files/", sep="")
IndcumPerfGraphFolder <- paste(getwd(), "/EXP1_NAc FR acquisition/Graphs/Behavior/Cumulative perf/", sep="")
PerfRelToCPFolder <- paste(getwd(), "/EXP1_NAc FR acquisition/Graphs/Behavior/Beh rel to CP/", sep="")
BySessFolder <- paste(MixedGraphFolder, "By Sess from CP/", sep="")
ScatterplotFolder <- paste(MixedGraphFolder, "100-400/Scatterplot/", sep="")
PrevsPostFolder <- paste(PerfRelToCPFolder, "Average Pre vs Post/", sep="")


### Load necessary functions
sourceDirectory(funcdirect, modifiedOnly=FALSE) #It doesn't work for some reason, just load functions separately while I fix it

load(file=paste(funcdirect, "MedPCextract.r", sep=""))
load(file=paste(funcdirect, "mpcextract_blockSingle.Rfunc", sep=""))
load(file=paste(funcdirect, "CPextract.R", sep=""))
load(file=paste(funcdirect, "cumulativeIndGraphs.R", sep=""))
load(file=paste(funcdirect, "PerformanceFromCP.R", sep=""))
load(file=paste(funcdirect, "PrePostCP_Perf.R", sep=""))
load(file=paste(funcdirect, "avgPerfByBin.R", sep=""))
load(file=paste(funcdirect, "behRASTERS.R", sep=""))
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
load(file=paste(funcdirect, "is.even.R", sep=""))
load(file=paste(funcdirect, "FRandPerf_Scatterplot.R", sep=""))
load(file=paste(funcdirect, "giveStars.R", sep=""))
load(file=paste(funcdirect, "plotBoxplotfromSessCP.R", sep=""))
load(file=paste(funcdirect, "plotFRBoxPlotandCP.R", sep=""))
load(file=paste(funcdirect, "DiffFR_ByBin.R", sep=""))
load(file=paste(funcdirect, "FRandPerf_Scatterplot.R", sep=""))
load(file=paste(funcdirect, "cumulative_CP.R", sep=""))
load(file=paste(funcdirect, "plotBoxplotPrePostCP.R", sep=""))
load(file=paste(funcdirect, "compareBins.R", sep=""))
load(file=paste(funcdirect, "compareDSvsNSfromCP.R", sep=""))
load(file=paste(funcdirect, "megaplotTC.R", sep=""))
load(file=paste(funcdirect, "orderCuesPerRat.R", sep=""))
load(file=paste(funcdirect, "PerfFromCPperRat.R", sep=""))
load(file=paste(funcdirect, "plotFRandCPhistogram_Entry.R", sep=""))


#### GENERATE IMPORTANT DATA OBJECTS
MedPCextract(MovAvg="Impinged only", funcdirect = funcdirect, datafolder = datafolder, 
             dataForRdir = dataForRdir, dataForRCumulative=dataForRCumulative, cuelength=10,
             DifferentTaskVersion=TRUE)

# Load important behavior-related objects
files <- paste(dataForRdir, list.files(dataForRdir), sep="")
filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep="")
for(i in 1:length(files)){load(files[[i]])}
for(i in 1:length(filesCum)){load(filesCum[[i]])}


#### BEHAVIORAL GRAPHS

# RASTERS
#This function makes all of the possible combination of rasters for all rats and all days
sapply(seq(1, length(rats)), function(x){
        sapply(seq(1, 6), function(y){
                sapply(seq(1, 2), function(z){
                        if(z==1){comboPick="FALSE"}
                        if(z==2){comboPick="TRUE"}
                        behRASTERS(subject = rats[x], day = y, lowerLimit=20, upperLimit=20, combined=comboPick, graphFolder=rasterGraphFolder)
                        
                })
        })
})



###
# I'M USING THE 10S PRE-CUE PERIOD AS BASELINE
# In 9% of the trials, the animal's head was inside the receptacle at the time the ITI window started. In those trials, I'll correct the ITI latency by assigning the MEAN ITI latency in the +/- 5 trials around that trial
###

#AVERAGE PERFORMANCE S+ VS S-. 5 trial bins.
avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(DStaskAcc, NStaskAcc), 
             cues=c("S+", "S-"), index="Specificity", legendLocation="topleft", y_axis_label="Performance index", 
             behGraphFolder=behGraphFolder, plot=T)
        
RR <- avgPerfByBin(binsize=5, colors=c(colindx[1], "darkblue"), data=list(DSrespAll, NSrespAll), 
             cues=c("S+", "S-"), index="Response ratio", legendLocation="topleft", y_axis_label="Response ratio", 
             behGraphFolder=behGraphFolder, plot=T)

Lat <- avgPerfByBin(binsize=5, colors=c(colindx[1], "darkblue", "gray30"), data=list(DSlatency, NSlatency, ITIlatency), 
             cues=c("S+", "S-", "ITI"), index="Latency", legendLocation="bottomright", y_axis_label="Latency (s)", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(ITIlatency), 
             cues=c("ITI"), index="Latency", legendLocation=-1, y_axis_label="Latency (s)", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c("gray50"), data=list(ITIrespRatio), 
             cues=c("ITI"), index="Response ratio", legendLocation=-1, y_axis_label="Response ratio", 
             behGraphFolder=behGraphFolder, plot=T, limit=240)



### CHANGE POINT: S+ SPECIFICITY
# Change point analysis using the method in Gallistel et al., 2004. I'll use S+ SPECIFICITY as the main performance index. But I also need to look at the results that other performance indexes give me.

CPdata <- CPextract(GallCrit=1.3, plot=T, CPfuncFolder=CPfuncFolder, idx=idx, CPGraphFolder=CPGraphFolder, dataForRdir=dataForRdir, dataForRCumulative=dataForRCumulative)



#Create a uniquely named object for the important objects of this experiment (in case I want to have thse objects for several experiments in my environment at the same time)
CPdataExp1 <- CPdata
alldataExp1 <- alldata
csacqidxExp1 <- csacqidx
ratsExp1 <- rats
idxExp1 <- idx

save(csacqidxExp1, file=paste(dataForRdir, "csacqidxExp1.ridx", sep=""))
save(alldataExp1, file=paste(dataForRdir, "alldataExp1.rdat", sep=""))
save(ratsExp1, file=paste(dataForRdir, "ratsExp1.rdat", sep=""))
save(idxExp1, file=paste(dataForRdir, "idxExp1.rdat", sep=""))

#### When in the session is the CP more likely to occur? ###
CPdata$CP - (CPdata$CPsess-1)*40
# [1] 37 17 26 17 13 19

t.test(x=CPdata$CP - (CPdata$CPsess-1)*40, y=rep(20, 6))
Welch Two Sample t-test

# data:  CPdata$CP - (CPdata$CPsess - 1) * 40 and rep(20, 6)
# t = 0.42174, df = 5, p-value = 0.6907
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#         -7.642747 10.642747
# sample estimates:
#         mean of x mean of y 
# 21.5      20.0 

####
# Performance from CP by animal. I'll make two graphs per rat: one with prob of entry and one with latency


orderCues <- orderCuesPerRat(alldata, rats, idx) #Order of cues. Useful for dealing with ITI cumulative data.
save(orderCues, file=paste(dataForRdir, "ordercues.rdat", sep=""))

#Response ratio
PerfFromCPperRat(behIndices=list(DSrespAll, NSrespAll, ITIrespRatio), beh="Response Ratio",
                 XLIM=c(-20, 20), YLIM=c(0, 1), binsize=5, numSess=7, numCues=40,
                 CPdata=CPdata, graphFolder=behGraphFolder, orderCues=orderCues,
                 color=c(colindx[1], "darkblue", "gray40"))

#Latency        
PerfFromCPperRat(behIndices=list(DSlatency, NSlatency, ITIlatency), beh="Latency",
                 XLIM=c(-20, 20), YLIM=c(0, 10), binsize=5, numSess=7, numCues=40,
                 CPdata=CPdata, graphFolder=behGraphFolder, orderCues=orderCues,
                 color=c(colindx[1], "darkblue", "gray40"))

#####
# CUMULATIVE PROPORTION OF ANIMALS WITH CHANGE POINT BY SESSION

cumulative_CP(Exp="Exp1", CPdata=CPdata, numSess=6, byAnimal=TRUE, byNeuron=FALSE, 
              nexdata=allNeuronsDS$nexdata, graphFolder=otherGraphFolder)

cumulative_CP(Exp="Exp1", CPdata=CPdata, numSess=6, byAnimal=F, byNeuron=T, 
              nexdata=allNeuronsDS$nexdata, graphFolder=otherGraphFolder)


#CHANGE POINT: Other performance indexes (response ratio, latency, ITI latency)
#This next line takes a long time. It's a calculation of change point using the other indexes (not the S+ specificity index). The method "intersect" calculates CPs using different criteria (from least to most conervative) and finds the earliest point that intersects
CPdataAllCrit <- CPextractMultipleCrit(GallCrit=c(1.3, 2, 4, 5, 6), idx, CPfuncFolder, CPGraphFolder, dataForRdir, dataForRCumulative)
RealCPbyIndex <- RealCP(data=CPdataAllCrit, method="intersect", 
       index_names=c("DSlatency", "DSrespAll", "DStaskAcc", "DStimeToSpare", "ITIlatency", 
                     "ITIrespRatio", "NSlatency", "NSrespAll", "NStaskAcc", "NStimeToSpare"), rats=rats)
        
#Find the CPs of the indexes I'm interested in:
targetIndexes <- c("DSrespAll", "DSlatency", "ITIlatency")
sel <- match(targetIndexes, rownames(RealCPbyIndex))
DSspecificity <- CPdata$CP #This one is the only one that can adopt negative values. So I had to calculate it differently.
CPframeExp1 <- rbind(DSspecificity, RealCPbyIndex[sel,])

save(CPframeExp1, file=paste(dataForRdir, "CPframeExp1.rdat", sep="")) #CP values for different behavioral indexes
save(CPdataExp1, file=paste(dataForRdir, "CPdataExp1.rdat", sep="")) #CP values and other details for PERFORMANCE INDEX (S+ specificity)


#CUMULATIVE INDIVIDUAL PERFORMANCE
cumulativeIndGraphs(numSess=7, numCues=40, sessLines=F, smartRat="NA", limit=280, 
                    dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, 
                    graphFolder=IndcumPerfGraphFolder, imgFormat="pdf")

#MEAN PERFORMANCE AS A FUNCTION OF DISTANCE FROM CP
PerformanceFromCP(relTrialMin=-100, relTrialMax=100, trialBinSize=5, typegraph="both", imgFormat="pdf", 
                  dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=PerfRelToCPFolder, CPdata=CPdata, 
                  idx=idx, rats = rats, alldata = alldata)


#PRE CP VS. POST CP PERFORMANCE
PrePostCP_Perf(data=DSrespAll, CP=CPdata$CP, y_axis_label="S+ response ratio", graphFolder=PrevsPostFolder, plot=T) #t = -6.0698, df = 5, p-value = 0.001753. No correction
PrePostCP_Perf(data=NSrespAll, CP=CPdata$CP, y_axis_label="S- response ratio", graphFolder=PrevsPostFolder, plot=T) #t = -0.018975, df = 5, p-value = 0.9856
PrePostCP_Perf(data=DSlatency, CP=CPdata$CP, y_axis_label="S+ latency", graphFolder=PrevsPostFolder, plot=T) #t = 6.8489, df = 5, p-value = 0.001013
PrePostCP_Perf(data=NSlatency, CP=CPdata$CP, y_axis_label="S- latency", graphFolder=PrevsPostFolder, plot=T) #t = -0.42179, df = 5, p-value = 0.6907
PrePostCP_Perf(data=DStaskAcc, CP=CPdata$CP, y_axis_label="S+ specificity", graphFolder=PrevsPostFolder, plot=T) #t = -11.968, df = 5, p-value = 7.181e-05
PrePostCP_Perf(data=NStaskAcc, CP=CPdata$CP, y_axis_label="S- specificity", graphFolder=PrevsPostFolder, plot=T) #t = 0.94549, df = 5, p-value = 0.3878
PrePostCP_Perf(data=ITIlatency, CP=CPdata$CP, y_axis_label="ITI latency", graphFolder=PrevsPostFolder, plot=T) #t = 0.18549, df = 5, p-value = 0.8601
PrePostCP_Perf(data=ITIrespRatio, CP=CPdata$CP, y_axis_label="ITI resp ratio", graphFolder=PrevsPostFolder, plot=T) #t = -0.41793, df = 5, p-value = 0.6933

p.adjust(c(0.001753, 7.18e-05, 0.001013, 0.8601), method="holm")
#0.0035060 0.0002872 0.0030390 0.8601000

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


### PREDICTORS OF CHANGE POINT
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

save(CPdata, file=paste(dataForRdir, "CPdata.rdat", sep=""))


#TOTAL TRIALS TO CP VS. NUMBER OF REWARDS
hist(CPdata$CP, breaks = seq(0, 250, 10), col="black", border="white")


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
allNeuronsEntryDS = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, BLwdw = 5, targetwdw = 1.5)
allNeuronsEntryNS = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, BLwdw = 5, targetwdw = 1.5)
allNeuronsEntryITI = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, BLwdw = 5, targetwdw = 1.5)


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
toplotTC <- megaplotTC(experiment="Exp1", data=list(allNeuronsDS), graphmin=500, graphmax=2000, cuewinmin=100, cuewinmax=400,
           colpalette="Rainbow", minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, title="DS onset",
           ZcolLabels=c("ZDS"), arrangeBy=c("ZDS"))

toplotTC <- megaplotTC(experiment="Exp1", data=list(allNeuronsNS), graphmin=500, graphmax=2000, cuewinmin=100, cuewinmax=400,
                       colpalette="Rainbow", minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, title="NS onset",
                       ZcolLabels=c("ZDS"), arrangeBy=c("ZDS"))

toplotTC <- megaplotTC(experiment="Exp1", data=list(allNeuronsEntryDS), graphmin=500, graphmax=5000, cuewinmin=100, cuewinmax=400,
                       colpalette="Rainbow", minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, title="DS entry",
                       ZcolLabels=c("ZDS"), arrangeBy=c("ZDS"))

toplotTC <- megaplotTC(experiment="Exp1", data=list(allNeuronsEntryNS), graphmin=500, graphmax=2000, cuewinmin=100, cuewinmax=400,
                       colpalette="Rainbow", minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, title="NS entry",
                       ZcolLabels=c("ZDS"), arrangeBy=c("ZDS"))

toplotTC <- megaplotTC(experiment="Exp1", data=list(allNeuronsEntryITI), graphmin=500, graphmax=2000, cuewinmin=100, cuewinmax=400,
                       colpalette="Rainbow", minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, title="ITI entry",
                       ZcolLabels=c("ZDS"), arrangeBy=c("ZDS"))

toplotTC <- megaplotTC(experiment="Exp1 Cue and Entry", data=list(allNeuronsDS, allNeuronsEntryDS, allNeuronsEntryNS, allNeuronsEntryITI), graphmin=500, 
                       graphmax=2000, cuewinmin=100, cuewinmax=400,
                       colpalette="Rainbow", minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, title="Entry",
                       ZcolLabels=c("ZDS", "ZDS entry", "Z NS entry", "Z ITI entry"), arrangeBy=c("ZDS"))


#Histograms with percentage of cue excited units per bin
CueExcCheck(experiment="Exp 3", data=toplotTC, graphFolder=MixedGraphFolder, 
            autom.FRrange=T, FRbin=1)
        

#BUILD A DATAFRAME WITH FIRING RATE OF EACH NEURON ON EACH TRIAL WITH BEH INFO ABOUT EACH TRIAL
masterDF_DS <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS, CPvector=CPdata$CP, sessionCPperRat = sessionCPperRat, funcdirect=funcdirect, dataForRdir=dataForRdir, BLduration=2, BLneudata=allNeuronsDS)
masterDF_DScorrect <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSresponded, CPvector=CPdata$CP, sessionCPperRat = sessionCPperRat, funcdirect=funcdirect, dataForRdir=dataForRdir, BLduration=2, BLneudata=allNeuronsDS)
masterDF_DSmissed <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSmissed, CPvector=CPdata$CP, sessionCPperRat = sessionCPperRat, funcdirect=funcdirect, dataForRdir=dataForRdir, BLduration=2, BLneudata=allNeuronsDS)
masterDF_NS <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsNS, CPvector=CPdata$CP, sessionCPperRat = sessionCPperRat, funcdirect=funcdirect, dataForRdir=dataForRdir, BLduration=2, BLneudata=allNeuronsNS)
masterDF_NSresponded <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNSresponded, 
                                           CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                           BLduration=2, sessionCPperRat=sessionCPperRat, BLneudata=allNeuronsNS)
masterDF_NSmissed<- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNSmissed, 
                                       CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                       BLduration=2, sessionCPperRat=sessionCPperRat, BLneudata=allNeuronsNS)


#The "entry-related" masterDF don't give me FR info of every trial, only of those trials in which the animal responded
masterDF_DSEntry <- FRbyNEURONbyBINcue(eventType="entry", cue="S+", neudata=allNeuronsEntryDS, 
                                       CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                       BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDS)

masterDF_NSEntry<- FRbyNEURONbyBINcue(eventType="entry", cue="S-", neudata=allNeuronsEntryNS, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS)

masterDF_ITIEntry<- FRbyNEURONbyBINcue(eventType="entry", cue="S+", neudata=allNeuronsEntryITI, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDS)


### BY MODALITY
masterDF_DS_Tone <- filter(masterDF_DS, masterDF_DS$modality==98)
masterDF_DS_Light <- filter(masterDF_DS, masterDF_DS$modality==99)
masterDF_NS_Light <- filter(masterDF_NS, masterDF_NS$modality==98)
masterDF_NS_Tone <- filter(masterDF_NS, masterDF_NS$modality==99)


save(masterDF_DS, file=paste(dataForRdir, "masterDF_DS.rdat", sep=""))
save(masterDF_DScorrect, file=paste(dataForRdir, "masterDF_DScorrect.rdat", sep=""))
save(masterDF_DSmissed, file=paste(dataForRdir, "masterDF_DSmissed.rdat", sep=""))
save(masterDF_NS, file=paste(dataForRdir, "masterDF_NS.rdat", sep=""))
save(masterDF_DSEntry, file=paste(dataForRdir, "masterDF_DSEntry.rdat", sep=""))
save(masterDF_NSEntry, file=paste(dataForRdir, "masterDF_NSEntry.rdat", sep=""))
save(masterDF_ITIEntry, file=paste(dataForRdir, "masterDF_ITIEntry.rdat", sep=""))

#masterDFsumm <- masterDFsummary(masterDF=masterDF_DS, neudata=allNeuronsDS, format="Zscores", postCueWdw=c(0, 0.4), postSpikeWdw=c(0.4, 5), orEntry=T, cuelength=10)


### PLOT FR POST-CUE (400MS WINDOW) as a function of distance to change point


#POST CUE (100-400ms)
plotFRandCP(cue=c("S+", "S-"), experiment="Exp 1 S+ vs S-", masterDF=list(masterDF_DS, masterDF_NS), 
            legLabels=c("S+", "S-"), graphFolder=MixedGraphFolder, trialBinSize=40, 
            WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-160, 160), yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, 
            neudata=list(allNeuronsDS, allNeuronsNS), cueExcOnly = T)


#I want to subset the cue-excited units. But I need to calculate that by bin (bc a unit might become excited towards the end of the session, for ex., but not be cue-excited at the beginning)
cueExcDataByBin <- cueExcByBin(masterDF=masterDF_DS, neudata=allNeuronsDS, capValue=c(-160, 160),
        trialBinSize=40, WdwStart=100, WdwEnd=400, BLforPoisson=5000)
        
plotFRandCP(cue=c("S+", "S-"), experiment="Exp 1 S+ Tone vs S- Light", masterDF=list(masterDF_DS_Tone, masterDF_NS_Light), 
            legLabels=c("S+", "S-"), graphFolder=MixedGraphFolder, trialBinSize=5, 
            WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-120, 120), yAxMinZ = -2, yAxMaxZ = 8, yAxMaxRaw = 10, 
            neudata=c(allNeuronsDS, allNeuronsNS), cueExcOnly = F)

plotFRandCP(cue=c("S+", "S-"), experiment="Exp 1 S+ Light vs S- Tone", masterDF=list(masterDF_DS_Light, masterDF_NS_Tone), 
            legLabels=c("S+", "S-"), graphFolder=MixedGraphFolder, trialBinSize=5, 
            WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-120, 120), yAxMinZ = -2, yAxMaxZ = 8, yAxMaxRaw = 10, 
            neudata=c(allNeuronsDS, allNeuronsNS), cueExcOnly = F)



#Units are considered cue excited if their FR on 3 consecutive 50ms windows exceeds the 99.9% upper limit of a Poisson distr. calculated using 5s baseline
PercResults <- plotFRandCP(cue=c("S+"), experiment="Exp 1 S+ PERCEXC", masterDF=list(masterDF_DS), 
            legLabels=c("S+"), graphFolder=MixedGraphFolder, trialBinSize=40, 
            WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, colindx=c(colindx[1], "gray30"), 
            capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 6, yAxMaxRaw = 10, 
            neudata=c(allNeuronsDS), cueExcOnly = F)
save(PercResults, file=paste(dataForRdir, "PercResults.rdat", sep = ""))



PercResults_Tone <- plotFRandCP(cue=c("S+"), experiment="Exp 1 S+ Tone PERCEXC", masterDF=list(masterDF_DS_Tone), 
                           legLabels=c("S+"), graphFolder=MixedGraphFolder, trialBinSize=40, 
                           WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), 
                           capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 6, yAxMaxRaw = 10, 
                           neudata=c(allNeuronsDS), cueExcOnly = F)

PercResults_Light <- plotFRandCP(cue=c("S+"), experiment="Exp 1 S+ Light PERCEXC", masterDF=list(masterDF_DS_Light), 
                                legLabels=c("S+"), graphFolder=MixedGraphFolder, trialBinSize=40, 
                                WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), 
                                capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 6, yAxMaxRaw = 10, 
                                neudata=c(allNeuronsDS), cueExcOnly = F)



#Is the proportion of cue-exc neurons modulated by the amount of training? (run a different chi.sq test for bins 1, 3 and 5 and for bins 2, 4 and 6, that way the units on each cell will be independent for sure)
#For Chi-sq analysis. I need to make sure that each cell has a completely different population of neurons. Bc each session has 40 trials of each kind, by only comparing bins that are separated by 40 trials I can guarantee that.
contingency_table_comp1<- t(PercResults[[1]][is.even(1:nrow(PercResults[[1]]))==FALSE, c(2, 3, 4)])
contingency_table_comp2 <- t(PercResults[[1]][is.even(1:nrow(PercResults[[1]])), c(2, 3, 4)])

inhCont_table_comp1 <- t(PercResults[[1]][is.even(1:nrow(PercResults[[1]]))==FALSE, c(2, 5, 6)])
inhCont_table_comp2 <- t(PercResults[[1]][is.even(1:nrow(PercResults[[1]])), c(2, 5, 6)])

#Chi-sq analysis excitations:
#Bins 1, 3, 5 and 7
critData <- contingency_table_comp1[-1, ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# data:  critData
# X-squared = 15.565, df = 2, p-value = 0.000417


#Bins 2, 4, 6, and 8
critData <- contingency_table_comp2[-1, ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# data:  critData
# X-squared = 18.895, df = 2, p-value = 7.887e-05


#Fisher's Exact test for count data. Actually the one I should use because I have cells with less than 5 observations
#Bins 1, 3, 5
critData <- contingency_table_comp1[-1, ]
fisher.test(critData)
#p-value = 0.0003168
#alternative hypothesis: two.sided

#Bins 2, 4, 6
critData <- contingency_table_comp2[-1, ]
fisher.test(critData)
# p-value = 8.643e-05
# alternative hypothesis: two.sided


#Post-hoc comparisons for % cue-excited units
# Is the distribution of exc/not exc neurons different from the first (-120 to -80) and the third (-40 to 0) bin. 
contTableExc <- PercResults[[1]][ , c(3, 4)]

#Bin 1 vs. 3
fisher.test(contTableExc[c(1, 3), ], alternative="less") #95% CI: 0.000000 0.9250037 Odds ratio=0.4547229 p=0.03221
#Bin 3 vs. 5
fisher.test(contTableExc[c(3, 5), ], alternative="less") #95% CI:0.0000000 0.7590285 odds ratio=0.317144 p=0.01145


#Bin 2 vs. 4
fisher.test((contTableExc[c(2, 4), ]), alternative="less") #95% CI: 0.0000000 0.9250037 odds ratio=0.4547229  p= 0.03221
#Bin 4 vs. 6
fisher.test((contTableExc[c(4, 6), ]), alternative="less") #95% CI:  0.000000 0.7590285 odds ratio=0.317144 p=0.01145
#Bin 6 vs. 8
fisher.test((contTableExc[c(6, 8), ]), alternative="less") #95% CI:  0.000000 4.598937 odds ratio=1.069993 p=0.6923

#Correct p value of 2 vs. 4 and 4 vs. 6 for multiple comparisons
p.adjust(p=c(0.03221, 0.01145), method="holm")
# 0.03221 0.02290

#Chi-sq for inhibitions:
#Bins 1, 3 and 5
critData <- inhCont_table_comp1[-1, ]
chisq.test(x=critData)

# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 8.4253, df = 2, p-value = 0.01481

#Bins 2, 4, 6 and 8
critData <- inhCont_table_comp2[-1, ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 7.2874, df = 2, p-value = 0.02615



#Fisher exact test for inhibitions:
#Bins 1, 3 and 5
critData <- inhCont_table_comp1[-1, ]
fisher.test(critData)
#p-value = 0.01341
#alternative hypothesis: two.sided

#Bins 2, 4 and 6
critData <- inhCont_table_comp2[-1, ]
fisher.test(critData)
# p-value = 0.03079
# alternative hypothesis: two.sided


#Post-hoc comparisons for % cue-inihibited units
# E.g. Is the distribution of inh/not inh neurons different from the first (-120 to -80) and the third (-40 to 0) bin. 
contTableInh <- PercResults[[1]][ , c(5, 6)]

#Bin 1 vs. 3
fisher.test((contTableInh[c(1, 3), ])) #95% CI:  0.153341 1.342072 odds ratio=0.4780756  p= 0.1712
#Bin 3 vs. 5
fisher.test((contTableInh[c(3, 5), ])) #95% CI:    1.459279 62.994158 odds ratio=6.688493 p=0.005111

p.adjust(p=c(0.1712, 0.005111), method="holm") # 0.171200 0.010222

#########

##### REPEAT THESE GRAPHS WITH 20 TRIAL BINS
#########
### BOXPLOTS vs. CHANGEPOINT BY BLOCKS OF TRIALS
####

#POST CUE (100-400ms)
PercResults20bin <- plotFRandCP(cue=c("S+"), experiment="Exp 1 S+ PERCEXC", masterDF=list(masterDF_DS), 
                                legLabels=c("S+"), graphFolder=MixedGraphFolder, trialBinSize=20, 
                                WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, colindx=c(colindx[1], "gray30"), 
                                capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 6, yAxMaxRaw = 10, 
                                neudata=c(allNeuronsDS), cueExcOnly = F)


save(PercResults20bin, file=paste(dataForRdir, "PercResults20bin.rdat", sep = ""))

# I'll refer to bins by numbers. The numbers refer to the order of the bin in a sequence of -120 to 120 trials by 20 bins. (On the graphs I may have less bins than this, so it could be confusing)
#Is the proportion of cue-exc neurons modulated by the amount of training? (run a different chi.sq test for bins 1, 3 and 5 and for bins 2, 4 and 6, that way the units on each cell will be independent for sure)
#For Chi-sq analysis. I need to make sure that each cell has a completely different population of neurons. Bc each session has 40 trials of each kind, by only comparing bins that are separated by 40 trials I can guarantee that.
contingency_table_comp1<- t(PercResults[[1]][is.even(1:nrow(PercResults[[1]]))==FALSE, c(2, 3, 4)])
contingency_table_comp2 <- t(PercResults[[1]][is.even(1:nrow(PercResults[[1]])), c(2, 3, 4)])

inhCont_table_comp1 <- t(PercResults[[1]][is.even(1:nrow(PercResults[[1]]))==FALSE, c(2, 5, 6)])
inhCont_table_comp2 <- t(PercResults[[1]][is.even(1:nrow(PercResults[[1]])), c(2, 5, 6)])

#Chi-sq analysis excitations:
#Bins 1, 3, 5, 7 and 9
critData <- contingency_table_comp1[-1, ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# data:  critData
# X-squared = 26.068, df = 4, p-value = 3.066e-05


#Bins 2, 4, 6, 8 and 10
critData <- contingency_table_comp2[-1, ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# data:  critData
# X-squared = 27.014, df = 4, p-value = 1.975e-05


#Fisher's Exact test for count data. Actually the one I should use because I have cells with less than 5 observations
#Bins 1, 3, 5, 7 and 9
critData <- contingency_table_comp1[-1, ]
fisher.test(critData)
#p-value = 2.549e-05
#alternative hypothesis: two.sided

#Bins 2, 4, 6, 8 and 10
critData <- contingency_table_comp2[-1, ]
fisher.test(critData)
# p-value = 1.547e-05
# alternative hypothesis: two.sided


#Post-hoc comparisons for % cue-excited units
# Is the distribution of exc/not exc neurons different from the first (-120 to -80) and the third (-40 to 0) bin. 
contTableExc <- PercResults[[1]][ , c(3, 4)]

#Bin 2 vs. 5
fisher.test((contTableExc[c(2, 5), ]), alternative="less") #95% CI: 0.0000000 0.4855377 odds ratio=0.2156722  p= 0.0003773
#Bin 5 vs. 8
fisher.test((contTableExc[c(5, 8), ]), alternative="less") #95% CI: 0.0000000 1.429954 odds ratio=0.4894083  p= 0.1645


#Correct p value of 2 vs. 5 and 5 vs. 8 for multiple comparisons
p.adjust(p=c(0.0003773, 0.1645), method="holm")
# 0.0007546 0.1645000

#Chi-sq for inhibitions:
#Bins 1, 3, 5, 7 and 9
critData <- inhCont_table_comp1[-1, ]
chisq.test(x=critData)

# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 15.531, df = 4, p-value = 0.003717

#Bins 2, 4, 6, 8 and 10
critData <- inhCont_table_comp2[-1, ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 15.036, df = 4, p-value = 0.004628



#Fisher exact test for inhibitions:
#Bins 1, 3, 5, 7 and 9
critData <- inhCont_table_comp1[-1, ]
fisher.test(critData)
#p-value = 0.006429
#alternative hypothesis: two.sided

#Bins 2, 4, 6, 8 and 10
critData <- inhCont_table_comp2[-1, ]
fisher.test(critData)
# p-value = 0.003442
# alternative hypothesis: two.sided


#Post-hoc comparisons for % cue-inihibited units
# E.g. Is the distribution of inh/not inh neurons different from the first (-120 to -80) and the third (-40 to 0) bin. 
contTableInh <- PercResults[[1]][ , c(5, 6)]

# #Bin 2 vs. 4
# fisher.test((contTableInh[c(2, 4), ])) #95% CI:  0.153341 1.342072 odds ratio=0.4780756  p= 0.1712
# #Bin 4 vs. 6
# fisher.test((contTableInh[c(4, 6), ])) #95% CI:  1.459279 62.994158 odds ratio=6.688493 p=0.005111
# #Bin 6 vs. 8
# fisher.test((contTableInh[c(6, 8), ])) #95% CI:   0.03885865 50.52994188 odds ratio=0.8035521 p=1


# #Bin 2 vs. 5
fisher.test((contTableInh[c(2, 5), ])) #95% CI:    0.830413 7.506103 odds ratio=2.41374 p=0.09297
# #Bin 5 vs. 8
fisher.test((contTableInh[c(5, 8), ])) #95% CI:   0.03885865 245.3073816 odds ratio=5.241182 p=0.1422




Down20bin <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 S+ - S- FR per bin All units", masterDF=list(masterDF_DS, masterDF_NS), 
                                   graphFolder=MixedGraphFolder, trialBinSize=20, WdwStart=100, WdwEnd=400, dataProcess="Zscores", 
                                   correctOnly=FALSE, points=F, lines=F, color=c(colindx[1], colindxB[1]), legLabels=c("S+", "S-"), 
                                   yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, capValue=c(-100, 100), cueExcOnly=F, neudata=allNeuronsDS)

uniqueBin <- unique(Down20bin[[1]]$bin)

Down20bin_Diff <- do.call("rbind", lapply(seq(1, length(uniqueBin)), function(i){
        DSbinSel <- Down20bin[[1]][Down20bin[[1]]$bin==uniqueBin[i], ]
        NSbinSel <- Down20bin[[2]][Down20bin[[2]]$bin==uniqueBin[i], ]
        
        matchUnitIdx <- match(DSbinSel$unitIdx, NSbinSel$unitIdx)
        NSfiring <- NSbinSel$byUnitFR[matchUnitIdx]
        
        DSbinSel$byUnitFR_NS <- NSfiring
        DSbinSel$DiffFR <- DSbinSel$byUnitFR-NSfiring
        
        DSbinSel
})
)


# This is S+ minus S- firing rate
DiffFR_ByBin(data=Down20bin_Diff, cueExcOnly=FALSE, color=colindx[1], ymin=-2, ymax=8,
             graphFolder=MixedGraphFolder, experiment="Exp 1 Down", points=F, 
             comparisons=list(c(3, 6), c(6, 9)))

# [[1]]
# idx      bin                       comp    W        p.val       p.adj sig
# W     1  1 vs. 3 -100 to -80 vs. -60 to -40  497 7.947317e-01 0.819523546    
# W4    2  2 vs. 4  -80 to -60 vs. -40 to -20 1134 6.409683e-02 0.320484158    
# W1    3  3 vs. 5    -60 to -40 vs. -20 to 0  574 1.556904e-02 0.093414259    
# W11   4  4 vs. 6     -40 to -20 vs. 0 to 20 1072 2.627417e-04 0.001839192  **
# W2    5  5 vs. 7      -20 to 0 vs. 20 to 40  382 4.060375e-05 0.000324830 ***
# W21   6  6 vs. 8       0 to 20 vs. 40 to 60  619 9.228628e-02 0.369145106    
# W3    7  7 vs. 9      20 to 40 vs. 60 to 80  267 2.822191e-01 0.819523546    
# W31   8 8 vs. 10     40 to 60 vs. 80 to 100  129 2.731745e-01 0.819523546    
# 
# [[2]]
# bins             comparison   W        p.val        p.adj sig
# W  3 vs. 6 -60 to -40 vs. 0 to 20 495 7.262019e-05 0.0001452404 ***
# W1 6 vs. 9   0 to 20 vs. 60 to 80 367 4.274001e-02 0.0427400065   *
        


#all units
compareBins(data=Down20bin_Diff, cueExcOnly=F, color=colindx[1], ymin=-2, ymax=8, 
            graphFolder=MixedGraphFolder, experiment="Exp 1 Down", points=TRUE, 
            comparisons=list(c(2, 5), c(5, 8)), cue="S+")
# [[1]]
# idx      bin         n                       comp     W        p.val        p.adj sig
# W     1  1 vs. 3 26 vs. 34 -100 to -80 vs. -60 to -40 542.0 9.327052e-01 0.9327051585    
# W4    2  2 vs. 4 45 vs. 61  -80 to -60 vs. -40 to -20 906.0 1.447801e-03 0.0086868046  **
# W1    3  3 vs. 5 34 vs. 47    -60 to -40 vs. -20 to 0 555.0 9.601479e-03 0.0480073958   *
# W11   4  4 vs. 6 61 vs. 56     -40 to -20 vs. 0 to 20 940.5 1.426351e-05 0.0001141080 ***
# W2    5  5 vs. 7 47 vs. 33      -20 to 0 vs. 20 to 40 388.5 7.923350e-05 0.0005546345 ***
# W21   6  6 vs. 8 56 vs. 27       0 to 20 vs. 40 to 60 718.0 3.577396e-01 0.8867818440    
# W3    7  7 vs. 9 33 vs. 18      20 to 40 vs. 60 to 80 269.0 2.955939e-01 0.8867818440    
# W31   8 8 vs. 10 27 vs. 11     40 to 60 vs. 80 to 100 114.0 1.391865e-01 0.5567460278    
# 
# [[2]]
# bins         n              comparison   W        p.val        p.adj sig
# W  2 vs. 5 45 vs. 47 -80 to -60 vs. -20 to 0 564 4.193377e-05 8.386754e-05 ***
# W1 5 vs. 8 47 vs. 27   -20 to 0 vs. 40 to 60 463 2.722932e-02 2.722932e-02   *
###


#Cue-excited units only
compareBins(data=Down20bin_Diff, cueExcOnly=T, color=colindx[1], ymin=-2, ymax=8, 
            graphFolder=MixedGraphFolder, experiment="Exp 1 Down Cue-Exc only 20 trial", points=FALSE, 
            comparisons=list(c(2, 5), c(5, 8)), cue="S+")
# [[1]]
# idx      bin         n                       comp     W        p.val       p.adj sig
# W     1  1 vs. 3  9 vs. 14 -100 to -80 vs. -60 to -40  74.0 0.7616466183 0.924436641    
# W4    2  2 vs. 4 14 vs. 32  -80 to -60 vs. -40 to -20 140.0 0.0226524589 0.113570551    
# W1    3  3 vs. 5 14 vs. 32    -60 to -40 vs. -20 to 0 137.0 0.0189284252 0.113570551    
# W11   4  4 vs. 6 32 vs. 40     -40 to -20 vs. 0 to 20 338.5 0.0003235092 0.002588073  **
# W2    5  5 vs. 7 32 vs. 28      -20 to 0 vs. 20 to 40 229.5 0.0006184408 0.004329085  **
# W21   6  6 vs. 8 40 vs. 22       0 to 20 vs. 40 to 60 433.0 0.4622183207 0.924436641    
# W3    7  7 vs. 9 28 vs. 15      20 to 40 vs. 60 to 80 179.0 0.2207384276 0.662215283    
# W31   8 8 vs. 10  22 vs. 9     40 to 60 vs. 80 to 100  75.0 0.1562933670 0.625173468    
# 
# [[2]]
# bins         n              comparison   W        p.val        p.adj sig
# W  2 vs. 5 14 vs. 32 -80 to -60 vs. -20 to 0  89 0.0004408038 0.0008816076 ***
# W1 5 vs. 8 32 vs. 22   -20 to 0 vs. 40 to 60 247 0.0327505716 0.0327505716   *
#         

#S+ vs. S-
PostCueFR_DSvsNS_fromCP_20trials <- compareDSvsNSfromCP(masterDF=list(masterDF_DS, masterDF_NS), cueExcOnly=F,
                                               trialBinSize=20, event="cue", correctOnly=F, WdwStart=100, 
                                               WdwEnd=400, capped=T, capValue=c(-100, 100), dataProcess="Zscores")
PostCueFR_DSvsNS_fromCP_20trials$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP_20trials$p, method="holm") 
PostCueFR_DSvsNS_fromCP_20trials$sig <- giveStars(PostCueFR_DSvsNS_fromCP_20trials$p.adj)

save(PostCueFR_DSvsNS_fromCP_20trials, file=paste(dataForRdir, "PostCueFR_DSvsNS_fromCP_20trials.rdat"))
# bin    V            p  n        p.adj sig
# V  -100 to -80  252 2.625954e-02 26 1.050382e-01    
# V1  -80 to -60  434 1.761315e-01 45 3.043375e-01    
# V2  -60 to -40  339 1.521687e-01 33 3.043375e-01    
# V3  -40 to -20 1167 5.621178e-02 61 1.686354e-01    
# V4    -20 to 0  857 7.666897e-04 47 4.600138e-03  **
# V5     0 to 20 1455 4.274842e-08 56 3.847358e-07 ***
# V6    20 to 40  555 1.629815e-09 33 1.629815e-08 ***
# V7    40 to 60  357 3.330410e-06 27 2.664328e-05 ***
# V8    60 to 80  156 5.226135e-04 18 3.658295e-03  **
# V9   80 to 100   64 1.464844e-03 11 7.324219e-03  **

### post-spike window (750-2000)

### POST SPIKE FIRING RATE ALL UNITS
Down20binPostSpike <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 POST SPIKE S+ - S- FR per bin ALL units", 
                                         masterDF=list(masterDF_DS, masterDF_NS), graphFolder=MixedGraphFolder, 
                                         trialBinSize=20, WdwStart=750, WdwEnd=2000, dataProcess="Zscores", 
                                         correctOnly=FALSE, points=F, lines=F, color=c(colindx[1], colindxB[1]), 
                                         legLabels=c("S+", "S-"), yAxMinZ=-1, yAxMaxZ=4, yAxMaxRaw=7, capped=T, 
                                         capValue=c(-100, 100), cueExcOnly=F, neudata=allNeuronsDS)


uniqueBin <- unique(Down20binPostSpike[[1]]$bin)

Down20binPostSpike_Diff <- do.call("rbind", lapply(seq(1, length(uniqueBin)), function(i){
        DSbinSel <- Down20binPostSpike[[1]][Down20binPostSpike[[1]]$bin==uniqueBin[i], ]
        NSbinSel <- Down20binPostSpike[[2]][Down20binPostSpike[[2]]$bin==uniqueBin[i], ]
        
        matchUnitIdx <- match(DSbinSel$unitIdx, NSbinSel$unitIdx)
        NSfiring <- NSbinSel$byUnitFR[matchUnitIdx]
        
        DSbinSel$byUnitFR_NS <- NSfiring
        DSbinSel$DiffFR <- DSbinSel$byUnitFR-NSfiring
        
        DSbinSel
})
)


### This plots boxplots of FR from CP in the bins designated when creating Down20bin_Diff and it also gives me analyses (compares the bins I indicate in the "comparisons" parameter)
compareBins(data=Down20binPostSpike_Diff, cueExcOnly=FALSE, color=colindx[1], ymin=-1, ymax=4, 
            graphFolder=MixedGraphFolder, experiment="Exp 1 Down POST SPIKE All units", points=FALSE, 
            comparisons=list(c(2, 5), c(5, 8)), 
            cue="S+")

# [[1]]
# idx      bin         n                       comp      W       p.val      p.adj sig
# W     1  1 vs. 3 26 vs. 34 -100 to -80 vs. -60 to -40  341.0 0.067294841 0.26917937    
# W4    2  2 vs. 4 45 vs. 61  -80 to -60 vs. -40 to -20  967.0 0.004816695 0.02890017   *
# W1    3  3 vs. 5 34 vs. 47    -60 to -40 vs. -20 to 0  666.0 0.103019787 0.30905936    
# W11   4  4 vs. 6 61 vs. 56     -40 to -20 vs. 0 to 20 1170.0 0.001680087 0.01344069   *
# W2    5  5 vs. 7 47 vs. 33      -20 to 0 vs. 20 to 40  658.0 0.127218906 0.30905936    
# W21   6  6 vs. 8 56 vs. 27       0 to 20 vs. 40 to 60  860.0 0.845128305 0.84512831    
# W3    7  7 vs. 9 33 vs. 18      20 to 40 vs. 60 to 80  207.5 0.039692430 0.19846215    
# W31   8 8 vs. 10 27 vs. 11     40 to 60 vs. 80 to 100   67.0 0.003838393 0.02686875   *
#         
#         [[2]]
# bins         n              comparison   W        p.val        p.adj sig
# W  2 vs. 5 45 vs. 47 -80 to -60 vs. -20 to 0 601 0.0001455178 0.0002910357 ***
# W1 5 vs. 8 47 vs. 27   -20 to 0 vs. 40 to 60 658 0.6055382061 0.6055382061    


### POST SPIKE FIRING RATE cue excited only
Down20binPostSpike <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 POST SPIKE S+ - S- FR per bin CUE EXC", 
                                         masterDF=list(masterDF_DS, masterDF_NS), graphFolder=MixedGraphFolder, 
                                         trialBinSize=20, WdwStart=750, WdwEnd=2000, dataProcess="Zscores", 
                                         correctOnly=FALSE, points=F, lines=F, color=c(colindx[1], colindxB[1]), 
                                         legLabels=c("S+", "S-"), yAxMinZ=-1, yAxMaxZ=4, yAxMaxRaw=7, capped=T, 
                                         capValue=c(-100, 100), cueExcOnly=T, neudata=allNeuronsDS)


uniqueBin <- unique(Down20binPostSpike[[1]]$bin)

Down20binPostSpike_Diff <- do.call("rbind", lapply(seq(1, length(uniqueBin)), function(i){
        DSbinSel <- Down20binPostSpike[[1]][Down20binPostSpike[[1]]$bin==uniqueBin[i], ]
        NSbinSel <- Down20binPostSpike[[2]][Down20binPostSpike[[2]]$bin==uniqueBin[i], ]
        
        matchUnitIdx <- match(DSbinSel$unitIdx, NSbinSel$unitIdx)
        NSfiring <- NSbinSel$byUnitFR[matchUnitIdx]
        
        DSbinSel$byUnitFR_NS <- NSfiring
        DSbinSel$DiffFR <- DSbinSel$byUnitFR-NSfiring
        
        DSbinSel
})
)


### This plots boxplots of FR from CP in the bins designated when creating Down20bin_Diff and it also gives me analyses (compares the bins I indicate in the "comparisons" parameter)
compareBins(data=Down20binPostSpike_Diff, cueExcOnly=TRUE, color=colindx[1], ymin=-1, ymax=4, 
            graphFolder=MixedGraphFolder, experiment="Exp 1 Down POST SPIKE", points=FALSE, 
            comparisons=list(c(2, 5), c(5, 8)), 
            cue="S+")
# [[1]]
# idx      bin         n                       comp     W       p.val      p.adj sig
# W     1  1 vs. 3  9 vs. 14 -100 to -80 vs. -60 to -40  41.0 0.089680490 0.27766273    
# W4    2  2 vs. 4 14 vs. 32  -80 to -60 vs. -40 to -20 161.5 0.069415683 0.27766273    
# W1    3  3 vs. 5 14 vs. 32    -60 to -40 vs. -20 to 0 128.0 0.010664665 0.06398799    
# W11   4  4 vs. 6 32 vs. 40     -40 to -20 vs. 0 to 20 466.0 0.024423477 0.12211739    
# W2    5  5 vs. 7 32 vs. 28      -20 to 0 vs. 20 to 40 472.0 0.640719673 1.00000000    
# W21   6  6 vs. 8 40 vs. 22       0 to 20 vs. 40 to 60 506.0 0.834683505 1.00000000    
# W3    7  7 vs. 9 28 vs. 15      20 to 40 vs. 60 to 80 116.5 0.008895272 0.06226690    
# W31   8 8 vs. 10  22 vs. 9     40 to 60 vs. 80 to 100  36.0 0.002459415 0.01967532   *
#         
#         [[2]]
# bins         n              comparison   W       p.val       p.adj sig
# W  2 vs. 5 14 vs. 32 -80 to -60 vs. -20 to 0  96 0.000858768 0.001717536  **
# W1 5 vs. 8 32 vs. 22   -20 to 0 vs. 40 to 60 412 0.855238031 0.855238031    



#S+ vs. S-
PostSpikeFR_DSvsNS_fromCP_20 <- compareDSvsNSfromCP(masterDF=list(masterDF_DS, masterDF_NS), cueExcOnly=F,
                                                 trialBinSize=20, event="cue", correctOnly=F, WdwStart=750, 
                                                 WdwEnd=2000, capped=T, capValue=c(-100, 100), dataProcess="Zscores")
PostSpikeFR_DSvsNS_fromCP_20$p.adj <- p.adjust(PostSpikeFR_DSvsNS_fromCP_20$p, method="holm") 
PostSpikeFR_DSvsNS_fromCP_20$sig <- giveStars(PostSpikeFR_DSvsNS_fromCP_20$p.adj)

save(PostSpikeFR_DSvsNS_fromCP_20, file=paste(dataForRdir, "PostSpikeFR_DSvsNS_fromCP_20.rdat"))
# bin    V            p  n        p.adj sig
# V  -100 to -80  133 1.456881e-01 26 2.913763e-01    
# V1  -80 to -60  500 4.792712e-01 44 4.792712e-01    
# V2  -60 to -40  389 2.640985e-02 33 7.922954e-02    
# V3  -40 to -20 1566 4.227612e-06 61 3.804851e-05 ***
# V4    -20 to 0  955 6.272229e-06 47 4.390561e-05 ***
# V5     0 to 20 1366 1.835889e-06 56 1.835889e-05 ***
# V6    20 to 40  509 4.541478e-06 33 3.804851e-05 ***
# V7    40 to 60  289 7.532187e-03 27 3.012875e-02   *
# V8    60 to 80  152 1.167297e-03 18 5.836487e-03  **
# V9   80 to 100   66 4.882812e-04 11 2.929687e-03  **
#         


############ 40 TRIAL BINS ################
###
Down40bin <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 S+ - S- FR per bin ALL UNITS", 
                                masterDF=list(masterDF_DS, masterDF_NS), 
                                graphFolder=MixedGraphFolder, trialBinSize=40, 
                                WdwStart=100, WdwEnd=400, dataProcess="Zscores", 
                                correctOnly=FALSE, points=F, lines=F, color=c(colindx[1], colindxB[1]), 
                                legLabels=c("S+", "S-"), 
                                yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, 
                                capValue=c(-120, 120), cueExcOnly=F, neudata=allNeuronsDS)

uniqueBin <- unique(Down40bin[[1]]$bin)

Down40bin_Diff <- do.call("rbind", lapply(seq(1, length(uniqueBin)), function(i){
        DSbinSel <- Down40bin[[1]][Down40bin[[1]]$bin==uniqueBin[i], ]
        NSbinSel <- Down40bin[[2]][Down40bin[[2]]$bin==uniqueBin[i], ]
        
        matchUnitIdx <- match(DSbinSel$unitIdx, NSbinSel$unitIdx)
        NSfiring <- NSbinSel$byUnitFR[matchUnitIdx]
        
        DSbinSel$byUnitFR_NS <- NSfiring
        DSbinSel$DiffFR <- DSbinSel$byUnitFR-NSfiring
        
        DSbinSel
})
)


DiffFR_ByBin(data=Down40bin_Diff, cueExcOnly=T, color=colindx[1], ymin=-2, ymax=18,
             graphFolder=MixedGraphFolder, experiment="Exp 3 Down", points=TRUE, 
             comparisons=list(c(1, 3), c(3, 5)))
# [[1]]
# bin                        comp     W        p.val        p.adj sig
# 1 vs. 3 -160 to -120 vs. -80 to -40 155 9.432529e-01 9.432529e-01    
# 2 vs. 4    -120 to -80 vs. -40 to 0 240 5.267141e-02 2.106857e-01    
# 3 vs. 5      -80 to -40 vs. 0 to 40 164 2.881846e-06 1.729107e-05 ***
# 4 vs. 6       -40 to 0 vs. 40 to 80 396 1.090979e-02 5.454893e-02    
# 5 vs. 7       0 to 40 vs. 80 to 120 269 1.814523e-01 5.443569e-01    
# 6 vs. 8     40 to 80 vs. 120 to 160 147 1.852654e-01 5.443569e-01    

# [[2]]
# bins                comparison   W      p.val      p.adj sig
# 1 vs. 4 -160 to -120 vs. -40 to 0 183 0.23663546 0.23663546    
# 4 vs. 6     -40 to 0 vs. 40 to 80 396 0.01090979 0.02181957   *
        
        
        
### This function plots boxplots of FR from CP in the bins designated when creating Down40bin_Diff and it also gives me analyses (compares the bins I indicate in the "comparisons" parameter)

#all units
compareBins(data=Down40bin_Diff, cueExcOnly=F, color=colindx[1], ymin=-2, ymax=8, 
            graphFolder=MixedGraphFolder, experiment="Exp 3 Down", points=TRUE, 
            comparisons=list(c(1, 3), c(3, 5)), cue="S+")
# [[1]]
# idx     bin         n                     comp    W        p.val        p.adj sig
# W     1 1 vs. 3 45 vs. 68 -120 to -80 vs. -40 to 0 1438 2.957515e-01 2.957515e-01    
# W2    2 2 vs. 4 55 vs. 63   -80 to -40 vs. 0 to 40  622 1.061581e-09 4.246326e-09 ***
# W1    3 3 vs. 5 68 vs. 37    -40 to 0 vs. 40 to 80  800 1.074478e-03 3.223434e-03  **
# W11   4 4 vs. 6 63 vs. 17    0 to 40 vs. 80 to 120  426 9.992554e-02 1.998511e-01    
# 
# [[2]]
# bins         n               comparison    W       p.val       p.adj sig
# W  1 vs. 3 45 vs. 68 -120 to -80 vs. -40 to 0 1438 0.295751522 0.295751522    
# W1 3 vs. 5 68 vs. 37    -40 to 0 vs. 40 to 80  800 0.001074478 0.002148956  **
        



# cue-excited only
compareBins(data=Down40bin_Diff, cueExcOnly=TRUE, color=colindx[1], ymin=-2, ymax=9, 
            graphFolder=MixedGraphFolder, experiment="Exp 3 Down Cue Exc", points=FALSE, 
            comparisons=list(c(1, 3), c(3, 5)), 
            cue="S+")

# [[1]]
#   bin            n                        comp   W        p.val        p.adj sig
# 1 1 vs. 3 11 vs. 21 -160 to -120 vs. -80 to -40  80 9.049446e-01 0.9049445865    
# 2 2 vs. 4 17 vs. 39    -120 to -80 vs. -40 to 0 103 2.719981e-02 0.1182131791    
# 3 3 vs. 5 21 vs. 46      -80 to -40 vs. 0 to 40  93 2.460918e-05 0.0001476551 ***
# 4 4 vs. 6 39 vs. 30       -40 to 0 vs. 40 to 80 279 2.364264e-02 0.1182131791    
# 5 5 vs. 7 46 vs. 14       0 to 40 vs. 80 to 120 141 4.373402e-02 0.1312020497    
# 6 6 vs. 8 30 vs. 12     40 to 80 vs. 120 to 160 126 4.516285e-01 0.9032570276    
# 
# [[2]]
# bins           n               comparison   W      p.val      p.adj sig
# 2 vs. 4 17 vs. 39 -120 to -80 vs. -40 to 0 103 0.02719981 0.04728527   *
# 4 vs. 6 39 vs. 30    -40 to 0 vs. 40 to 80 279 0.02364264 0.04728527   *


#### ANOVA WITH NON-CONSECUTIVE BINS
Down40binCueExc <- Down40bin_Diff[Down40bin_Diff$CueExcited==TRUE, ]
OddBins <- c(-120, -40, 40) #Bins 1, 3 and 5
EvenBins <- c(-80, 0, 80) #Bins 2, 4, and 6

OddBinData <- Down40binCueExc[Down40binCueExc$bin %in% OddBins, ]
EvenBinData <- Down40binCueExc[Down40binCueExc$bin %in% EvenBins, ]

ezANOVA(data=OddBinData, dv=byUnitFR, between = bin, wid=unitIdx)
# $`ANOVA`
#       Effect DFn DFd      F            p  p<.05       ges
# 1       bin   1  84  9.905037 0.002281307     * 0.1054793

ezANOVA(data=EvenBinData, dv=byUnitFR, between = bin, wid=unitIdx)
# $`ANOVA`
# Effect DFn DFd        F          p            p<.05        ges
#  bin   1   79         20.72529 1.896121e-05     * 0.2078238





####

### POST SPIKE FIRING RATE ALL UNITS
Down40binPostSpike <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 3 POST SPIKE S+ - S- FR per bin ALL units", 
                                         masterDF=list(masterDF_DS, masterDF_NS), graphFolder=MixedGraphFolder, 
                                         trialBinSize=40, WdwStart=750, WdwEnd=2000, dataProcess="Zscores", 
                                         correctOnly=FALSE, points=F, lines=F, color=c(colindx[1], colindxB[1]), 
                                         legLabels=c("S+", "S-"), yAxMinZ=-1, yAxMaxZ=4, yAxMaxRaw=7, capped=T, 
                                         capValue=c(-160, 160), cueExcOnly=F, neudata=allNeuronsDS)


uniqueBin <- unique(Down40binPostSpike[[1]]$bin)

Down40binPostSpike_Diff <- do.call("rbind", lapply(seq(1, length(uniqueBin)), function(i){
        DSbinSel <- Down40binPostSpike[[1]][Down40binPostSpike[[1]]$bin==uniqueBin[i], ]
        NSbinSel <- Down40binPostSpike[[2]][Down40binPostSpike[[2]]$bin==uniqueBin[i], ]
        
        matchUnitIdx <- match(DSbinSel$unitIdx, NSbinSel$unitIdx)
        NSfiring <- NSbinSel$byUnitFR[matchUnitIdx]
        
        DSbinSel$byUnitFR_NS <- NSfiring
        DSbinSel$DiffFR <- DSbinSel$byUnitFR-NSfiring
        
        DSbinSel
})
)


### This plots boxplots of FR from CP in the bins designated when creating Down40bin_Diff and it also gives me analyses (compares the bins I indicate in the "comparisons" parameter)
compareBins(data=Down40binPostSpike_Diff, cueExcOnly=FALSE, color=colindx[1], ymin=-1, ymax=4, 
            graphFolder=MixedGraphFolder, experiment="Exp 3 Down POST SPIKE All units", points=FALSE, 
            comparisons=list(c(1, 4), c(4, 6)), 
            cue="S+")

# [[1]]
# idx     bin         n                        comp    W        p.val        p.adj   sig
# 1     1 vs. 3 33 vs. 55 -160 to -120 vs. -80 to -40 1014 8.217967e-01 8.217967e-01    
# 2     2 vs. 4 45 vs. 68    -120 to -80 vs. -40 to 0 1325 1.151826e-01 3.455478e-01    
# 3     3 vs. 5 55 vs. 63      -80 to -40 vs. 0 to 40  860 1.274572e-06 7.647430e-06 ***
# 4     4 vs. 6 68 vs. 37       -40 to 0 vs. 40 to 80 1002 4.327797e-02 1.731119e-01    
# 5     5 vs. 7 63 vs. 17       0 to 40 vs. 80 to 120  336 9.629274e-03 4.814637e-02   *
# 6     6 vs. 8 37 vs. 15     40 to 80 vs. 120 to 160  224 1.439113e-01 3.455478e-01    
# 
# [[2]]
# bins         n               comparison         W      p.val      p.adj sig
# W  2 vs. 4 45 vs. 68 -120 to -80 vs. -40 to 0 1325 0.11518260 0.11518260    
# W1 4 vs. 6 68 vs. 37    -40 to 0 vs. 40 to 80 1002 0.04327797 0.08655594    





### POST SPIKE FIRING RATE cue excited only
Down40binPostSpike <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 3 POST SPIKE S+ - S- FR per bin CUE EXC", 
                                         masterDF=list(masterDF_DS, masterDF_NS), graphFolder=MixedGraphFolder, 
                                         trialBinSize=40, WdwStart=750, WdwEnd=2000, dataProcess="Zscores", 
                                         correctOnly=FALSE, points=F, lines=F, color=c(colindx[1], colindxB[1]), 
                                         legLabels=c("S+", "S-"), yAxMinZ=-1, yAxMaxZ=4, yAxMaxRaw=7, capped=T, 
                                         capValue=c(-160, 160), cueExcOnly=T, neudata=allNeuronsDS)


uniqueBin <- unique(Down40binPostSpike[[1]]$bin)

Down40binPostSpike_Diff <- do.call("rbind", lapply(seq(1, length(uniqueBin)), function(i){
        DSbinSel <- Down40binPostSpike[[1]][Down40binPostSpike[[1]]$bin==uniqueBin[i], ]
        NSbinSel <- Down40binPostSpike[[2]][Down40binPostSpike[[2]]$bin==uniqueBin[i], ]
        
        matchUnitIdx <- match(DSbinSel$unitIdx, NSbinSel$unitIdx)
        NSfiring <- NSbinSel$byUnitFR[matchUnitIdx]
        
        DSbinSel$byUnitFR_NS <- NSfiring
        DSbinSel$DiffFR <- DSbinSel$byUnitFR-NSfiring
        
        DSbinSel
})
)


### This plots boxplots of FR from CP in the bins designated when creating Down40bin_Diff and it also gives me analyses (compares the bins I indicate in the "comparisons" parameter)
compareBins(data=Down40binPostSpike_Diff, cueExcOnly=TRUE, color=colindx[1], ymin=-1, ymax=4, 
            graphFolder=MixedGraphFolder, experiment="Exp 3 Down POST SPIKE", points=FALSE, 
            comparisons=list(c(2, 4), c(4, 6)), 
            cue="S+")

# [[1]]
# bin         n                        comp      W        p.val       p.adj sig
# 1 vs. 3 11 vs. 21 -160 to -120 vs. -80 to -40 105 0.3480414104 0.696082821    
# 2 vs. 4 17 vs. 39    -120 to -80 vs. -40 to 0 199 0.0088059517 0.035223807   *
# 3 vs. 5 21 vs. 46      -80 to -40 vs. 0 to 40 233 0.0002699168 0.001619501  **
# 4 vs. 6 39 vs. 30       -40 to 0 vs. 40 to 80 487 0.1198184086 0.359455226    
# 5 vs. 7 46 vs. 14       0 to 40 vs. 80 to 120 159 0.0018476242 0.009238121  **
# 6 6 vs. 8 30 vs. 12     40 to 80 vs. 120 to 160 170 0.3970997605 0.696082821    
# 
# [[2]]
# bins         n               comparison   W       p.val     p.adj sig
# 2 vs. 4 17 vs. 39 -120 to -80 vs. -40 to 0 199 0.008805952 0.0176119   *
# 4 vs. 6 39 vs. 30    -40 to 0 vs. 40 to 80 487 0.119818409 0.1198184    



####
# BY MODALITY

byBinbyUnit_Tone <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 3 S+ Tone vs S- Light", 
                    masterDF=list(masterDF_DS_Tone, masterDF_NS_Light), 
                    graphFolder=MixedGraphFolder, trialBinSize=20, points = F,
                    WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                    color=c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), yAxMinZ=-2, yAxMaxZ=10, yAxMaxRaw=7, 
                    capped=T, capValue=c(-120, 120), cueExcOnly=T, neudata=allNeuronsDS)

byBinbyUnit_Light <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 3 S+ Light vs S- Tone", 
                                       masterDF=list(masterDF_DS_Light, masterDF_NS_Tone), 
                                       graphFolder=MixedGraphFolder, trialBinSize=20, points = F,
                                       WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                       color=c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), yAxMinZ=-2, yAxMaxZ=10, yAxMaxRaw=7, 
                                       capped=T, capValue=c(-120, 120), cueExcOnly=T, neudata=allNeuronsDS)



### ANALYZE NON-CONSECUTIVE BINS
allBins <- unique(byBinbyUnit[[1]]$bin)
allBinsIdx <- 1:length(allBins)
binGroups <- list(c(allBins[is.even(allBinsIdx)]), 
                  c(allBins[is.even(allBinsIdx)==FALSE]))

comp <- c("S+", "S-")

        
forAnalysis <- lapply(seq(1, length(binGroups)), function(d){
                
        do.call("rbind", lapply(seq(1, length(byBinbyUnit)), function(c){
                selBinGroup <- binGroups[[d]]
                
                selDat <- byBinbyUnit[[c]][byBinbyUnit[[c]]$bin %in% selBinGroup, ]
                
                idx <- 1:nrow(selDat)
                
                data.frame(idx, cue=comp[c], binGroup=d, bin=selDat$bin, byUnitFR=selDat$byUnitFR)
                
        })
        )
        
})

names(forAnalysis) <- binGroups

# GROUP #1 OF NON-CONSECUTIVE BINS: c(-80 to -40, 0 to 40 and 80 to 120)
ezANOVA(data=forAnalysis[[1]], between = c(cue, bin), dv=byUnitFR, wid=idx)
# $`ANOVA`
#    Effect DFn DFd        F            p p<.05       ges
# 1     cue   1 266 50.86604 9.354200e-12     * 0.1605285
# 2     bin   1 266 41.42174 5.689220e-10     * 0.1347391
# 3 cue:bin   1 266 32.10929 3.782278e-08     * 0.1077098

# GROUP #2 OF NON-CONSECUTIVE BINS: c(-120 to -80, -40 to 0 and 40 to 80)
ezANOVA(data=forAnalysis[[2]], between = c(cue, bin), dv=byUnitFR, wid=idx)
# $`ANOVA`
# Effect DFn DFd        F            p p<.05        ges
# 1     cue   1 296 29.79268 1.018027e-07     * 0.09144675
# 2     bin   1 296 12.44113 4.866125e-04     * 0.04033551
# 3 cue:bin   1 296 15.86105 8.585727e-05     * 0.05085936

# -120 to 80 vs. -40 to 0
bin1 <- forAnalysis[[2]][forAnalysis[[2]]$bin==-120, ]; bin2 <- forAnalysis[[2]][forAnalysis[[2]]$bin==-40, ]
wilcox.test(bin1$byUnitFR, bin2$byUnitFR, paired=F, alternative="less") #W = 5827, p-value = 0.2716

# -40 to 0 vs. 40 to 80
bin1 <- forAnalysis[[2]][forAnalysis[[2]]$bin==-40, ]; bin2 <- forAnalysis[[2]][forAnalysis[[2]]$bin==40, ]
wilcox.test(bin1$byUnitFR, bin2$byUnitFR, paired=F, alternative="less") #W = 4362, p-value = 0.05575

# -80 to -40 vs. 0 to 40
bin1 <- forAnalysis[[1]][forAnalysis[[1]]$bin==-80, ]; bin2 <- forAnalysis[[1]][forAnalysis[[1]]$bin==-0, ]
wilcox.test(bin1$byUnitFR, bin2$byUnitFR, paired=F, alternative="less") #W = 4238, p-value = 1.342e-07

# 0 to 40 vs. 80 to 120
bin1 <- forAnalysis[[1]][forAnalysis[[1]]$bin==0, ]; bin2 <- forAnalysis[[1]][forAnalysis[[1]]$bin==80, ]
wilcox.test(bin1$byUnitFR, bin2$byUnitFR, paired=F, alternative="less") #W = 1999, p-value = 0.2761






#ENTRY
plotFRandCP_Entry(experiment="Exp 3", cue="S+", wdwLabel="Post S+ Entry", masterDF=list(masterDF_EntryDS), graphFolder=MixedGraphFolder, 
            trialBinSize=30, WdwStart =0, WdwEnd=1000, dataProcess="Zscores", colindx =colindx[1], 
            capped=T, capValue = c(-90, 90), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS)

plotFRandCP_Entry(experiment="Exp 3", cue="S-", wdwLabel="Post S- Entry", masterDF=list(masterDF_EntryNS), graphFolder=MixedGraphFolder, 
                  trialBinSize=30, WdwStart =0, WdwEnd=1000, dataProcess="Zscores", colindx =colindx[1], 
                  capped=T, capValue = c(-90, 90), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryNS)


plotFRandCP_Entry(experiment="Exp 3", cue=c("S+", "S-"), wdwLabel="Post S+ and S- Entry", masterDF=list(masterDF_DSEntry, masterDF_NSEntry), graphFolder=MixedGraphFolder, 
                  trialBinSize=40, WdwStart =0, WdwEnd=2000, dataProcess="Zscores", colindx =c(colindx[1], "darkblue"), 
                  capped=T, capValue = c(-160, 160), yAxMinZ = -1, yAxMaxZ = 2, legLabels=c("S+", "S-"),
                  yAxMaxRaw = 10, neudata=allNeuronsEntryDS)

plotFRandCP_Entry(experiment="Exp 3", cue=c("S+", "S-"), wdwLabel="Pre S+ and S- Entry", masterDF=list(masterDF_DSEntry, masterDF_NSEntry), graphFolder=MixedGraphFolder, 
                  trialBinSize=40, WdwStart =-2000, WdwEnd=0, dataProcess="Zscores", colindx =c(colindx[1], "darkblue"), 
                  capped=T, capValue = c(-160, 160), yAxMinZ = -1, yAxMaxZ = 2, legLabels=c("S+", "S-"),
                  yAxMaxRaw = 10, neudata=allNeuronsEntryDS)


################################################################################################
### CALCULATE DIFFERENCE IN FR AROUND S+ VS ROUND AS A FUNCTION OF DISTANCE TO CHANGE POINT
###############################################################################################

#Post cue
PostCueFR_DSvsNS_fromCP <- compareDSvsNSfromCP(masterDF=list(masterDF_DS, masterDF_NS), cueExcOnly=F,
                                               trialBinSize=40, event="cue", correctOnly=F, WdwStart=100, 
                                               WdwEnd=400, capped=T, capValue=c(-120, 120), dataProcess="Zscores")
PostCueFR_DSvsNS_fromCP$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP$p, method="holm") 
PostCueFR_DSvsNS_fromCP$sig <- giveStars(PostCueFR_DSvsNS_fromCP$p.adj)

save(PostCueFR_DSvsNS_fromCP, file=paste(dataForRdir, "PostCueFR_DSvsNS_fromCP.rdat"))

#            bin    V            p  n        p.adj sig
# V  -120 to -80  779 1.314045e-03 45 3.942136e-03  **
# V1  -80 to -40  781 4.649484e-01 55 4.649484e-01    
# V2    -40 to 0 1633 2.494777e-03 68 4.989553e-03  **
# V3     0 to 40 1865 2.263487e-09 63 1.358092e-08 ***
# V4    40 to 80  647 4.608228e-07 37 2.304114e-06 ***
# V5   80 to 120  148 7.629395e-05 17 3.051758e-04 ***
        


#Post spike
PostSpikeFR_DSvsNS_fromCP <- compareDSvsNSfromCP(masterDF=list(masterDF_DS, masterDF_NS), cueExcOnly=F,
                                               trialBinSize=40, event="cue", correctOnly=F, WdwStart=750, 
                                               WdwEnd=2000, capped=T, capValue=c(-160, 160), dataProcess="Zscores")
PostSpikeFR_DSvsNS_fromCP$p.adj <- p.adjust(PostSpikeFR_DSvsNS_fromCP$p, method="holm") 
PostSpikeFR_DSvsNS_fromCP$sig <- giveStars(PostSpikeFR_DSvsNS_fromCP$p.adj)

save(PostSpikeFR_DSvsNS_fromCP, file=paste(dataForRdir, "PostSpikeFR_DSvsNS_fromCP.rdat"))
# bin           V            p          n        p.adj  sig
# -160 to -120  440 1.771138e-03        33 5.313414e-03  **
# -120 to -80   421 1.407504e-01        45 2.815008e-01    
# -80 to -40    667 2.592153e-01        54 2.815008e-01    
# -40 to 0      1969 5.846024e-07       68 4.092217e-06 ***
# 0 to 40       1813 1.817857e-08       63 1.454285e-07 ***
# 40 to 80      576 2.167893e-04        37 1.281738e-03  **
# 80 to 120     139 8.392334e-04        17 3.356934e-03  **
# 120 to 160    116 2.136230e-04        15 1.281738e-03  **


#Post cue by modality
PostCueFR_DSvsNS_fromCP_Tone <- compareDSvsNSfromCP(masterDF=list(masterDF_DS_Tone, masterDF_NS_Light), 
                                               trialBinSize=20, event="cue", correctOnly=F, WdwStart=100, 
                                               WdwEnd=400, capped=T, capValue=c(-120, 120), dataProcess="Zscores")
PostCueFR_DSvsNS_fromCP_Tone$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP_Tone$p, method="holm") 
PostCueFR_DSvsNS_fromCP_Tone$sig <- giveStars(PostCueFR_DSvsNS_fromCP_Tone$p.adj)

save(PostCueFR_DSvsNS_fromCP_Tone, file=paste(dataForRdir, "PostCueFR_DSvsNS_fromCP_Tone.rdat"))

#
PostCueFR_DSvsNS_fromCP_Light <- compareDSvsNSfromCP(masterDF=list(masterDF_DS_Light, masterDF_NS_Tone), 
                                                    trialBinSize=20, event="cue", correctOnly=F, WdwStart=100, 
                                                    WdwEnd=400, capped=T, capValue=c(-120, 120), dataProcess="Zscores")
PostCueFR_DSvsNS_fromCP_Light$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP_Light$p, method="holm") 
PostCueFR_DSvsNS_fromCP_Light$sig <- giveStars(PostCueFR_DSvsNS_fromCP_Light$p.adj)

save(PostCueFR_DSvsNS_fromCP_Light, file=paste(dataForRdir, "PostCueFR_DSvsNS_fromCP_Light.rdat"))

#Number of units per bin
plot.new()
nUnits <- PostCueFR_DSvsNS_fromCP$n
bp <- barplot(nUnits, col="black", names.arg = seq(-90, 75, 15), axes=T, axisnames=T)
text(x=bp, y=3, labels = nUnits, col="white", font=3)
mtext(side=1,line=2.5, font=2, cex=1.4, text="Bin")
mtext(side=2,line=2.5, font=2, cex=1.4, text="Number of units")


#Pre entry
PreEntryFR_DSvsNS_fromCP <- compareDSvsNSfromCP(masterDF=list(masterDF_DSEntry, masterDF_NSEntry), 
                                                trialBinSize=40, event="entry", correctOnly=F, 
                                                WdwStart=-2000, WdwEnd=0, capped=T, 
                                                capValue=c(-160, 160), dataProcess="Zscores")
PreEntryFR_DSvsNS_fromCP$p.adj <- p.adjust(PreEntryFR_DSvsNS_fromCP$p, method="holm") 
PreEntryFR_DSvsNS_fromCP$sig <- giveStars(PreEntryFR_DSvsNS_fromCP$p.adj)
save(PreEntryFR_DSvsNS_fromCP, file=paste(dataForRdir, "PreEntryFR_DSvsNS_fromCP.rdat"))
# bin           V            p          n        p.adj sig
# -160 to -120  195 2.142429e-03        21 1.071215e-02   *
# -120 to -80   224 1.699257e-02        24 3.398514e-02   *
# -80 to -40    363 3.249185e-02        32 3.398514e-02   *
# -40 to 0      808 6.931978e-07        42 4.852384e-06 ***
# 0 to 40       477 1.429580e-07        31 1.143664e-06 ***
# 40 to 80      130 2.136230e-04        16 1.281738e-03  **
# 80 to 120     120 2.578735e-03        16 1.071215e-02   *
# 120 to 160    44 3.906250e-03         9 1.171875e-02   *


#Post entry
PostEntryFR_DSvsNS_fromCP <- compareDSvsNSfromCP(masterDF=list(masterDF_DSEntry, masterDF_NSEntry), 
                                                 trialBinSize=40, event="entry", correctOnly=F, 
                                                 WdwStart=0, WdwEnd=500, capped=T, 
                                                 capValue=c(-120, 120), dataProcess="Zscores")
PostEntryFR_DSvsNS_fromCP$p.adj <- p.adjust(PostEntryFR_DSvsNS_fromCP$p, method="holm") 
PostEntryFR_DSvsNS_fromCP$sig <- giveStars(PostEntryFR_DSvsNS_fromCP$p.adj)
save(PostEntryFR_DSvsNS_fromCP, file=paste(dataForRdir, "PostEntryFR_DSvsNS_fromCP.rdat"))

#           bin   V            p  n        p.adj sig
# V  -120 to -80 257 6.998777e-04 24 3.499389e-03  **
# V1  -80 to -40 346 6.438768e-02 32 1.931630e-01    
# V2    -40 to 0 785 4.372096e-06 42 2.623258e-05 ***
# V3     0 to 40 235 6.028901e-01 31 6.028901e-01    
# V4    40 to 80  33 3.695679e-02 16 1.478271e-01    
# V5   80 to 120  53 2.318726e-01 16 4.637451e-01    


 


# ### PLOT FR IN THE POSTCUE WINDOW AND TAIL PERIOD BY RAT BY TRIAL BY NEURON
# PLOT_perRatTrial_FRandBEH(masterDFsumm=masterDFsumm, wdw="PostCue", threshold=2, MixedGraphFolder=MixedGraphFolder)
# PLOT_perRatTrial_FRandBEH(masterDFsumm=masterDFsumm, wdw="PostSpike", threshold=2, MixedGraphFolder=MixedGraphFolder)


### PLOT FR AROUND THE CUE (PTSH) as a function of distance to change point
plotFRandCPhistogram(experiment="Exp 1 Cue Exc only", masterDF=list(masterDF_DS, masterDF_NS), graphFolder=MixedGraphFolder, 
                     trialBinSize=40, dataProcess="Zscores", correctOnly=FALSE, color=c(colindx[1], "darkblue"), 
                     capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 3, cueExcOnly=F,
                     psthmin=0.5, psthmax=2, imgFormat="pdf", neudata=allNeuronsDS)
        

### PLOT FR AROUND ENTRIES (PTSH) as a function of distance to change point
plotFRandCPhistogram_Entry(experiment="Exp 1 Entry Hist over 2s lat", masterDF=list(masterDF_EntryDS, masterDF_EntryNS), 
                           graphFolder=MixedGraphFolder, cueExcOnly=F, trialBinSize=40, dataProcess="Zscores", 
                           correctOnly=FALSE, color=c(colindx[1], "darkblue", "gray40"), 
                     capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 2, yAxMaxRaw = 3, 
                     WdwStart=-2000, WdwEnd=5000, imgFormat="pdf", neudata=allNeuronsEntryDS)

plotFRandCPhistogram_Entry(experiment="Exp 1 Entry Hist over 5s lat", masterDF=list(masterDF_EntryDS, masterDF_EntryNS), 
                           graphFolder=MixedGraphFolder, cueExcOnly=F, trialBinSize=40, dataProcess="Zscores", 
                           correctOnly=FALSE, color=c(colindx[1], "darkblue", "gray40"), 
                           capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 2, yAxMaxRaw = 3, 
                           WdwStart=-5000, WdwEnd=5000, imgFormat="pdf", neudata=allNeuronsEntryDS)

plotFRandCPhistogram_Entry(experiment="Exp 1 Entry Hist over 7s lat", masterDF=list(masterDF_EntryDS, masterDF_EntryNS), 
                           graphFolder=MixedGraphFolder, cueExcOnly=F, trialBinSize=40, dataProcess="Zscores", 
                           correctOnly=FALSE, color=c(colindx[1], "darkblue", "gray40"), 
                           capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 2, yAxMaxRaw = 3, 
                           WdwStart=-7000, WdwEnd=5000, imgFormat="pdf", neudata=allNeuronsEntryDS)



### PLOT FR AROUND ENTRIES (CP) as a function of distance to change point

plotFRandCP_Entry(cue=c("S+", "S-"), wdwLabel="Post S+", experiment="Exp 1 FR 2s preentry when latency is more than 5s", masterDF=list(masterDF_EntryDS, masterDF_EntryNS), 
                  graphFolder=MixedGraphFolder, trialBinSize=40, WdwStart=-2000, WdwEnd=0, dataProcess="Zscores", minLat=5000,
                  correctOnly=F, colindx=c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), yAxMinZ=-1, yAxMaxZ=2, yAxMaxRaw=7, 
                  capped=T, capValue=c(-120, 120), cueExcOnly=F, neudata=allNeuronsDS)

plotFRandCP_Entry(cue=c("S+", "S-"), wdwLabel="Post S+", experiment="Exp 1 FR 2s preentry when latency is more than 7s", masterDF=list(masterDF_EntryDS, masterDF_EntryNS), 
                  graphFolder=MixedGraphFolder, trialBinSize=40, WdwStart=-1000, WdwEnd=0, dataProcess="Zscores", minLat=6000,
                  correctOnly=F, colindx=c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), yAxMinZ=-1, yAxMaxZ=3, yAxMaxRaw=7, 
                  capped=T, capValue=c(-120, 120), cueExcOnly=F, neudata=allNeuronsDS)



#Same thing but by session from CP instead of trial
plotPSTHfromSessCP(experiment="Exp 1", masterDF=list(masterDF_DS_Tone, masterDF_DS_Light), 
                   graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                   correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-2, 2), 
                   yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, psthmin=-0.5, psthmax=2, 
                   imgFormat="pdf", neudata=allNeuronsDS)

### PLOT BOXPLOTS AROUND CUE BY SESS FROM CP AND DRUG GROUP
plotBoxplotfromSessCP(experiment="Exp 1", masterDF=list(masterDF_DS_Tone, masterDF_DS_Light), 
                      comp=c("Tone S+", "Light S+"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-5, 5), 
                      yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)




###################################################################33
### I want to analyze SESSION #1 ONLY. Bc all of my analyses are organized around time of CPsess, I need to change my CPdata so that CPsess is 1 for all of my subjects and then analyze sess 1 alone
CPdata$CP <- 1
CPdata$CPsess <- 1 #Remember to reset to the real CPdata object if you want to do real CP-related analyses after this

#Session in which the CP took place per rat (each item is a rat, the order is given by the object "rats")
CPvector <- CPdata$CP
sessionCPperRat <- sapply(seq(1, length(CPvector)), function(l){findInterval(CPvector[l], ratsTrialIdxPerSess[[l]])+1})

masterDF_DS_From1 <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS, 
                                        CPvector=CPdata$CP, sessionCPperRat = sessionCPperRat, 
                                        funcdirect=funcdirect, dataForRdir=dataForRdir, BLduration=2, 
                                        BLneudata=allNeuronsDS)

### BY MODALITY
masterDF_DS_Tone_From1 <- filter(masterDF_DS_From1, masterDF_DS$modality==98)
masterDF_DS_Light_From1 <- filter(masterDF_DS_From1, masterDF_DS$modality==99)
masterDF_NS_Light_From1 <- filter(masterDF_NS_From1, masterDF_NS$modality==98)
masterDF_NS_Tone_From1 <- filter(masterDF_NS_From1, masterDF_NS$modality==99)


save(masterDF_DS_From1, file=paste(dataForRdir, "masterDF_DS_From1.rdat"))
save(masterDF_DS_Tone_From1, file=paste(dataForRdir, "masterDF_DS_Tone_From1.rdat"))
save(masterDF_DS_Light_From1, file=paste(dataForRdir, "masterDF_DS_Light_From1.rdat"))
save(masterDF_NS_Tone_From1, file=paste(dataForRdir, "masterDF_NS_Tone_From1.rdat"))
save(masterDF_NS_Light_From1, file=paste(dataForRdir, "masterDF_NS_Light_From1.rdat"))

plotFRandCP(cue = "S+", experiment = "Exp 1", trialBinSize = 40, WdwStart=100, WdwEnd=400,
            dataProcess="PercCueExc", capped=T, capValue=c(0, 240), masterDF=list(masterDF_DS_From1))



plotBoxplotfromSessCP(experiment="Exp 1_Sess1-6 by Modality w Lines", 
                      masterDF=list(masterDF_DS_Tone_From1, masterDF_NS_Light_From1,
                                    masterDF_DS_Light_From1, masterDF_NS_Tone_From1), 
                      comp=c("S+", "S-", "S+", "S-"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=TRUE, color=c(colindx[1], colindxB[1], colindx[2], colindxB[2]), 
                      sessFromCP = c(0, 6), lines=T, analysis=F,
                      yAxMinZ=-2, yAxMaxZ = 6, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=F,
                      factors=list(Modality=c("Tone S+", "Light S+"), Cue=c("S+", "S-")))


plotBoxplotfromSessCP(experiment="Exp 1_Sess1-6", masterDF=list(masterDF_DS_From1, masterDF_NS_From1), 
                      comp=c("S+", "S-"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(0, 5), 
                      yAxMinZ=-2, yAxMaxZ = 6, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=F,
                      factors = list(cue=c("S+", "S-")))

# [[1]]
# [[1]]$`ANOVA`
#       Effect DFn DFd         F            p p<.05        ges
# 2      sess   5 358  4.870363 2.526248e-04     * 0.06368955
# 3      comp   1 358 70.641734 1.024361e-15     * 0.16480368
# 4 sess:comp   5 358  3.889315 1.909275e-03     * 0.05152140
# 
# [[1]]$`Levene's Test for Homogeneity of Variance`
# DFn DFd      SSn     SSd        F            p p<.05
# 1  11 358 236.6918 1530.46 5.033286 2.595365e-07     *
#         
#         
#[[2]]
#       W     p.val                 p.adj  sig
# W     892     4.231460e-04 1.269438e-03  **
# W1    693     7.324301e-03 1.464860e-02   *
# W2    531     4.025788e-01 4.025788e-01    
# W3    1141    4.981296e-05 1.992519e-04 ***
# W4    417     6.790527e-06 3.395264e-05 ***
# W5    595     2.260649e-10 1.356389e-09 ***

plotPSTHfromSessCP(experiment="Exp 3_Sess1-6", masterDF=list(masterDF_DS_From1, masterDF_NS_From1), 
                   graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                   correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(0, 5), 
                   yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, psthmin=-0.5, psthmax=2, 
                   imgFormat="pdf", neudata=allNeuronsDS)



### BY MODALITY

plotBoxplotfromSessCP(experiment="Exp 1_Sess1-6_By modality S+ tone", masterDF=list(masterDF_DS_Tone_From1, masterDF_NS_Light_From1), 
                      comp=c("S+", "S-"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(0, 5), 
                      yAxMinZ=-2, yAxMaxZ = 6, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=F,
                      factors = list(cue=c("S+", "S-")))

plotBoxplotfromSessCP(experiment="Exp 1_Sess1-6_By modality S+ light", masterDF=list(masterDF_DS_Light_From1, masterDF_NS_Tone_From1), 
                      comp=c("S+", "S-"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(0, 5), 
                      yAxMinZ=-2, yAxMaxZ = 6, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=F,
                      factors = list(cue=c("S+", "S-")))



#####################################################################

### THIS ONE IS THE SAME BUT GROUPING THE SESSIONS BEFORE AND AFTER THE CP
plotBoxplotPrePostCP(experiment="Exp 3", masterDF=list(masterDF_DS_Tone, masterDF_DS_Light), factors=list(mod=c("Tone", "Light")),
                     comp=c("Tone S+", "Light S+"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                     correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-5, 5), 
                     yAxMinZ = -2, yAxMaxZ = 8, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                     removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)



# $`ANOVA`
# Effect          DFn   DFd           F            p p<.05          ges
# Pre vs Post     1     145 29.35540883 2.448763e-07     * 1.683653e-01
# Tone vs. Light  1     145  0.00562708 9.403071e-01       3.880594e-05
# Interaction     1     145  0.07788280 7.805850e-01       5.368344e-04
# 
# $`Levene's Test for Homogeneity of Variance`
# DFn DFd      SSn      SSd        F            p p<.05
# 1   3 145 160.8292 877.9472 8.854075 1.989338e-05     *


#               Comparison              W        p.val        p.adj sig
# W     Pre CP Tone vs. Pre CP Light 1082 1.359494e-01 2.718987e-01    
# W1  Post CP Tone vs. Post CP Light  140 5.278148e-01 5.278148e-01    
# W2    Pre CP Tone vs. Post CP Tone  199 1.730899e-07 6.923595e-07 ***
# W11 Pre CP Light vs. Post CP Light   85 8.237366e-03 2.471210e-02   *




# #COMPARE WITH OTHER GROUPS
# compareCPs(data=list(CPdataExp4, CPdataExp10tr, CPdataHYBunilAP5, CPdataHYBbilAP5) , imgFormat="png", expNames=c("Task1", "Task2", "UnilAP5", "BilAP5"), colindx=c("gray30", "gray45", "gray60", "gray75"), graphFolder=behGraphFolder, minSess=5, graphs=c("Change point", "Non learners", "Dynamic Interval", "Asymptote"))


# FOCUS ON THE UNITS RECORDED ON THE DAY BEHAVIOR CHANGED
UnitHeatMap(data=masterDF, sessFromCP=4, FRparameters=allNeuronsDS$parameters, 
            folder=BySessFolder, winmin=100, winmax=300, BLmin=-2000, BLmax=0)



# Megaplot with firing in the 400ms window after 4 events: S+ responded to, S+ missed, S- responded to, S- missed
# This function makes the plot and also saves a data frame that will be useful for the next graph
megaplot(data=list(allNeuronsDSresponded, allNeuronsDSmissed, allNeuronsNSresponded, allNeuronsNSmissed), 
         CPdata=CPdata, csacqidx=csacqidx, BLwdw=2, winmin=100, winmax=400, colpalette="Rainbow", minFR=-3, 
         maxFR=3, graphFolder=MixedGraphFolder, dataForRdir=dataForRdir, 
         ZcolLabels=c("ZDSresp", "ZDSmissed", "ZNSresp", "ZNSmissed"))
        
megaplot(data=list(allNeuronsDS, allNeuronsNS), CPdata=CPdata, csacqidx=csacqidx, BLwdw=2, winmin=50, winmax=400,
         colpalette="Rainbow", minFR=-2, maxFR=6, graphFolder=MixedGraphFolder, dataForRdir=dataForRdir,
         ZcolLabels=c("ZDS", "ZNS"), arrangeBy="ZDS")



## SCATTERPLOT WHERE EACH SESSION IS A DOT AND THE AXES ARE --> X axis: performance index; Y axis: mean FR. 
#Load the version of toplot that I need for this scatterplot. "toplot" is made and save by the function "megaplot". So run the megaplot line that will give me the desired "toplot"
load(file=paste(dataForRdir, "toplot.rdat", sep=""))


toplot$uniqSession <- paste(toplot$rat, toplot$expt)
### FILTER NEURONS THAT ARE s+ EXCITED ONLY
toplot <- toplot[toplot$DSExc==TRUE, ]

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

#Performance index
FRandPerf_Scatterplot(data=toplotSumm, graphFolder=ScatterplotFolder, CPtoo=TRUE, 
                      CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="PerfIndex")

#Latency
FRandPerf_Scatterplot(data=toplotSumm, graphFolder=ScatterplotFolder, CPtoo=TRUE, 
                      CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="Latency")

#Response ratio
FRandPerf_Scatterplot(data=toplotSumm, graphFolder=ScatterplotFolder, CPtoo=TRUE, 
                      CPmerged=FALSE, color=c("darkblue", colindx[1], colindxB[1]), 
                      behIndex="RespRatio")


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

### BY DAY
masterDF_DS_Day1_D <- masterDF_DS_From1[masterDF_DS_From1$sessfromCPsess==0, ]
masterDF_NS_Day1_D <- masterDF_NS_From1[masterDF_NS_From1$sessfromCPsess==0, ]

masterDF_DS_Day2_D <- masterDF_DS_From1[masterDF_DS_From1$sessfromCPsess==1, ]
masterDF_NS_Day2_D <- masterDF_NS_From1[masterDF_NS_From1$sessfromCPsess==1, ]

masterDF_DS_Day3_D <- masterDF_DS_From1[masterDF_DS_From1$sessfromCPsess==2, ]
masterDF_NS_Day3_D <- masterDF_NS_From1[masterDF_NS_From1$sessfromCPsess==2, ]

masterDF_DS_Day4_D <- masterDF_DS_From1[masterDF_DS_From1$sessfromCPsess==3, ]
masterDF_NS_Day4_D <- masterDF_NS_From1[masterDF_NS_From1$sessfromCPsess==3, ]

masterDF_DS_Day5_D <- masterDF_DS_From1[masterDF_DS_From1$sessfromCPsess==4, ]
masterDF_NS_Day5_D <- masterDF_NS_From1[masterDF_NS_From1$sessfromCPsess==4, ]

masterDF_DS_Day6_D <- masterDF_DS_From1[masterDF_DS_From1$sessfromCPsess==5, ]
masterDF_NS_Day6_D <- masterDF_NS_From1[masterDF_NS_From1$sessfromCPsess==5, ]


## DAY 1 DIFFERENT CONDITIONS
masterDF_DS_Day1_Tone_D <- masterDF_DS_From1[masterDF_DS_From1$sessfromCPsess==0 & masterDF_DS_From1$modality=="98", ]
masterDF_DS_Day1_Light_D <- masterDF_DS_From1[masterDF_DS_From1$sessfromCPsess==0 & masterDF_DS_From1$modality=="99", ]
masterDF_NS_Day1_Light_D <- masterDF_NS_From1[masterDF_NS_From1$sessfromCPsess==0 & masterDF_NS_From1$modality=="98", ]
masterDF_NS_Day1_Tone_D <- masterDF_NS_From1[masterDF_NS_From1$sessfromCPsess==0 & masterDF_NS_From1$modality=="99", ]



byBinbyUnit_Trial1to15_ByModality <- plotFRBoxPlotandCP(cue=c("S+", "S+", "S-", "S-"), experiment="Exp 1 Sess 1 Trial 1-5 BOXPLOTS DOWN S+ vs S-", 
                                                       masterDF=list(masterDF_DS_Day1_Tone_D, masterDF_DS_Day1_Light_D, 
                                                                     masterDF_NS_Day1_Light_D, masterDF_NS_Day1_Tone_D), 
                                                       points=T, graphFolder=MixedGraphFolder, trialBinSize=15, 
                                                       WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                                       color=c(colindx, colindxB), legLabels=c("S+ Tone D", "S+ Light D", "S- Light D", "S- Tone D"), 
                                                       yAxMinZ=-3, yAxMaxZ=10, yAxMaxRaw=7, capped=T, capValue=c(1, 15), cueExcOnly=F, 
                                                       neudata=allNeuronsDS)

byBinbyUnit_Trial20to40_ByModality <- plotFRBoxPlotandCP(cue=c("S+", "S+", "S-", "S-"), experiment="Exp 1 Sess 1 Trial 15-30 BOXPLOTS DOWN S+ vs S-", 
                                                        masterDF=list(masterDF_DS_Day1_Tone_D, masterDF_DS_Day1_Light_D, 
                                                                      masterDF_NS_Day1_Light_D, masterDF_NS_Day1_Tone_D), 
                                                        points=T, graphFolder=MixedGraphFolder, trialBinSize=20, 
                                                        WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                                        color=c(colindx, colindxB), legLabels=c("S+ Tone D", "S+ Light D", "S- Light D", "S- Tone D"), 
                                                        yAxMinZ=-3, yAxMaxZ=10, yAxMaxRaw=7, capped=T, capValue=c(20, 40), cueExcOnly=F, 
                                                        neudata=allNeuronsDS)

## both modalities together, by 10 trial bin: S+ vs. S-. For analysis.
byBinbyUnit_Sess1_Trial1to10 <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 Sess 1 1-10 BOXPLOTS DOWN S+ vs S-", 
                                                        masterDF=list(masterDF_DS_Day1_D, masterDF_NS_Day1_D), 
                                                        points=F, lines=T, graphFolder=MixedGraphFolder, trialBinSize=10, 
                                                        WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                                        color=c(colindx, colindxB), legLabels=c("S+ D","S- D"), 
                                                        yAxMinZ=-3, yAxMaxZ=10, yAxMaxRaw=7, capped=T, capValue=c(1, 10), cueExcOnly=F, 
                                                        neudata=allNeuronsDS)

byBinbyUnit_Sess1_Trial11to20 <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 Sess 1 11-20 BOXPLOTS DOWN S+ vs S-", 
                                                   masterDF=list(masterDF_DS_Day1_D, masterDF_NS_Day1_D), 
                                                   points=F, lines=T, graphFolder=MixedGraphFolder, trialBinSize=10, 
                                                   WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                                   color=c(colindx, colindxB), legLabels=c("S+ D","S- D"), 
                                                   yAxMinZ=-3, yAxMaxZ=10, yAxMaxRaw=7, capped=T, capValue=c(11, 20), cueExcOnly=F, 
                                                   neudata=allNeuronsDS)

byBinbyUnit_Sess1_Trial21to30 <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 Sess 1 11-20 BOXPLOTS DOWN S+ vs S-", 
                                                    masterDF=list(masterDF_DS_Day1_D, masterDF_NS_Day1_D), 
                                                    points=F, lines=T, graphFolder=MixedGraphFolder, trialBinSize=10, 
                                                    WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                                    color=c(colindx, colindxB), legLabels=c("S+ D","S- D"), 
                                                    yAxMinZ=-3, yAxMaxZ=10, yAxMaxRaw=7, capped=T, capValue=c(21, 30), cueExcOnly=F, 
                                                    neudata=allNeuronsDS)

byBinbyUnit_Sess1_Trial31to40 <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 Sess 1 11-20 BOXPLOTS DOWN S+ vs S-", 
                                                    masterDF=list(masterDF_DS_Day1_D, masterDF_NS_Day1_D), 
                                                    points=F, lines=T, graphFolder=MixedGraphFolder, trialBinSize=10, 
                                                    WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                                    color=c(colindx, colindxB), legLabels=c("S+ D","S- D"), 
                                                    yAxMinZ=-3, yAxMaxZ=10, yAxMaxRaw=7, capped=T, capValue=c(31, 40), cueExcOnly=F, 
                                                    neudata=allNeuronsDS)




nUnits <- length(byBinbyUnit_Sess1_Trial1to10[[1]]$unitIdx)

save(byBinbyUnit_Sess1_Trial1to10, file=paste(dataForRdir, "byBinbyUnit_Trial1to10.rdat", sep=""))


forAnalysis <- do.call("rbind", lapply(seq(1, length(byBinbyUnit_Sess1_Trial31to40)), function(x){
        colsForAnalysis <- c("unitIdx", "rat", "byUnitFR", "cue", "bin")
        colIdx <- colnames(byBinbyUnit_Sess1_Trial31to40[[x]]) %in% colsForAnalysis
        selForAnalysis <- byBinbyUnit_Sess1_Trial31to40[[x]][, colIdx]
        return(selForAnalysis)
})
)

ezANOVA(data=forAnalysis, wid=unitIdx, within = c("bin", "cue"), dv=byUnitFR)

wilcox.test(x=forAnalysis$byUnitFR[forAnalysis$cue=="S+"],
            y=forAnalysis$byUnitFR[forAnalysis$cue=="S-"],
            paired=F)
# Wilcoxon signed rank test
# 
# data:  forAnalysis$byUnitFR[forAnalysis$cue == "S+"] and forAnalysis$byUnitFR[forAnalysis$cue == "S-"]
# V = 537, p-value = 0.0001283
# alternative hypothesis: true location shift is not equal to 0

# Bin  1:10: W = 696, p-value = 0.3318
# Bin 11:20: W = 810, p-value = 0.02
# Bin 21:30: W = 873, p-value = 0.001941
# Bin 31: 40: W = 865, p-value = 0.002703

p.adjust(p=c(0.332, 0.02, 0.0019, 0.0027), method='holm') #0.3320 0.0400 0.0076 0.0081

##################################
### WHAT'S GOING ON ON DAY 1?  ###
##################################

#This will give me a summary by trial for S+ and for S- on day 1
Day1_TrialUnit <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 Sess 1 1-40 BOXPLOTS DOWN S+ vs S-", 
                                                   masterDF=list(masterDF_DS_Day1_D, masterDF_NS_Day1_D), 
                                                   points=F, lines=F, graphFolder=MixedGraphFolder, trialBinSize=5, 
                                                   WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                                   color=c(colindx, colindxB), legLabels=c("S+ D","S- D"), 
                                                   yAxMinZ=-3, yAxMaxZ=10, yAxMaxRaw=7, capped=T, capValue=c(0, 40), cueExcOnly=F, 
                                                   neudata=allNeuronsDS)

Day2_TrialUnit <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 Sess 2 40-80 BOXPLOTS DOWN S+ vs S-", 
                                     masterDF=list(masterDF_DS_Day2_D, masterDF_NS_Day2_D), 
                                     points=F, lines=F, graphFolder=MixedGraphFolder, trialBinSize=1, 
                                     WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                     color=c(colindx, colindxB), legLabels=c("S+ D","S- D"), 
                                     yAxMinZ=-3, yAxMaxZ=10, yAxMaxRaw=7, capped=T, capValue=c(41, 80), cueExcOnly=F, 
                                     neudata=allNeuronsDS)

Day3_TrialUnit <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 Sess 3 81-120 BOXPLOTS DOWN S+ vs S-", 
                                     masterDF=list(masterDF_DS_Day3_D, masterDF_NS_Day3_D), 
                                     points=F, lines=F, graphFolder=MixedGraphFolder, trialBinSize=1, 
                                     WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                     color=c(colindx, colindxB), legLabels=c("S+ D","S- D"), 
                                     yAxMinZ=-3, yAxMaxZ=10, yAxMaxRaw=7, capped=T, capValue=c(81, 120), cueExcOnly=F, 
                                     neudata=allNeuronsDS)

Day4_TrialUnit <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 Sess 4 121-160 BOXPLOTS DOWN S+ vs S-", 
                                     masterDF=list(masterDF_DS_Day4_D, masterDF_NS_Day4_D), 
                                     points=F, lines=F, graphFolder=MixedGraphFolder, trialBinSize=1, 
                                     WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                     color=c(colindx, colindxB), legLabels=c("S+ D","S- D"), 
                                     yAxMinZ=-3, yAxMaxZ=10, yAxMaxRaw=7, capped=T, capValue=c(121, 160), cueExcOnly=F, 
                                     neudata=allNeuronsDS)


Day5_TrialUnit <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 Sess 5 161-200 BOXPLOTS DOWN S+ vs S-", 
                                     masterDF=list(masterDF_DS_Day5_D, masterDF_NS_Day5_D), 
                                     points=F, lines=F, graphFolder=MixedGraphFolder, trialBinSize=1, 
                                     WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                     color=c(colindx, colindxB), legLabels=c("S+ D","S- D"), 
                                     yAxMinZ=-3, yAxMaxZ=10, yAxMaxRaw=7, capped=T, capValue=c(161, 200), cueExcOnly=F, 
                                     neudata=allNeuronsDS)

Day6_TrialUnit <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 1 Sess 6 201-240 BOXPLOTS DOWN S+ vs S-", 
                                     masterDF=list(masterDF_DS_Day6_D, masterDF_NS_Day6_D), 
                                     points=F, lines=F, graphFolder=MixedGraphFolder, trialBinSize=1, 
                                     WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
                                     color=c(colindx, colindxB), legLabels=c("S+ D","S- D"), 
                                     yAxMinZ=-3, yAxMaxZ=10, yAxMaxRaw=7, capped=T, capValue=c(201, 240), cueExcOnly=F, 
                                     neudata=allNeuronsDS)


save(Day1_TrialUnit, file=paste(dataForRdir, "Day1_TrialUnit.rdat", sep=""))
save(Day2_TrialUnit, file=paste(dataForRdir, "Day2_TrialUnit.rdat", sep=""))
save(Day3_TrialUnit, file=paste(dataForRdir, "Day3_TrialUnit.rdat", sep=""))
save(Day4_TrialUnit, file=paste(dataForRdir, "Day4_TrialUnit.rdat", sep=""))
save(Day5_TrialUnit, file=paste(dataForRdir, "Day5_TrialUnit.rdat", sep=""))
save(Day6_TrialUnit, file=paste(dataForRdir, "Day6_TrialUnit.rdat", sep=""))


#I'm going to make a pdf with a different animal per page and, for each animal, all of the recorded units and what they did on every trial

Session_Polygraph(data=Day1_TrialUnit, cues=c("S+", "S-"), day=1, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)
Session_Polygraph(data=Day1_TrialUnit, cues=c("S+"), day=1, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)
Session_Polygraph(data=Day1_TrialUnit, cues=c("S-"), day=1, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)

 
Session_Polygraph(data=Day2_TrialUnit, cues=c("S+", "S-"), day=2, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)
Session_Polygraph(data=Day2_TrialUnit, cues=c("S+"), day=2, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)
Session_Polygraph(data=Day2_TrialUnit, cues=c("S-"), day=2, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)

Session_Polygraph(data=Day3_TrialUnit, cues=c("S+", "S-"), day=3, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)
Session_Polygraph(data=Day3_TrialUnit, cues=c("S+"), day=3, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)
Session_Polygraph(data=Day3_TrialUnit, cues=c("S-"), day=3, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)

Session_Polygraph(data=Day4_TrialUnit, cues=c("S+", "S-"), day=4, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)
Session_Polygraph(data=Day4_TrialUnit, cues=c("S+"), day=4, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)
Session_Polygraph(data=Day4_TrialUnit, cues=c("S-"), day=4, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)

Session_Polygraph(data=Day5_TrialUnit, cues=c("S+", "S-"), day=5, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)
Session_Polygraph(data=Day5_TrialUnit, cues=c("S+"), day=5, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)
Session_Polygraph(data=Day5_TrialUnit, cues=c("S-"), day=5, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)

Session_Polygraph(data=Day6_TrialUnit, cues=c("S+", "S-"), day=6, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)
Session_Polygraph(data=Day6_TrialUnit, cues=c("S+"), day=6, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)
Session_Polygraph(data=Day6_TrialUnit, cues=c("S-"), day=6, FRmax=10, FRmin=-2, plotFolder=neuGraphFolder)



############################################################################################
### RELATIONSHIP BETWEEN FIRING RATE AND BEHAVIOR PER UNIT                               ###
############################################################################################

corr.FRandBeh.PerUnit(experiment= "Exp 3 pre and post CP", 
                masterDF=masterDF_DS, 
                WdwStart=100, WdwEnd=400, neudata=allNeuronsDS, 
                fromCP=TRUE, bySess=FALSE, sessRange = c(-5, 5), collapsePrePost=FALSE, 
                byTrial=FALSE, capped=TRUE, capValue=c(-90, 90), trialBinSize=15,
                cueExcOnly=FALSE, correctOnly=FALSE,
                graphFolder=MixedGraphFolder)



###################################################################################
### RELATIONSHIP BETWEEN CUE-evoked AND ENTRY-evoked FIRING RATE                ###
###################################################################################

# Data frame indicating which units were recorded before, during or after the change point session
unitsDF <- unitsVSchangepoint(FRdata=allNeuronsDS, CPdata=CPdata, collapseCP=T)

#############

toplotTC <- megaplotTC(experiment="Exp1 Cue and Entry", data=list(allNeuronsDS, allNeuronsEntryDS, allNeuronsEntryNS, allNeuronsEntryITI), 
                       graphmin=500, graphmax=2000, cuewinmin=100, cuewinmax=400,
                       colpalette="Rainbow", minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, title="Entry",
                       ZcolLabels=c("ZDS", "ZDS entry", "Z NS entry", "Z ITI entry"), arrangeBy=c("ZDS"))

#Use the object I just created to compare FR after S+ vs. after entry. If I change graphmin and graphmax in the previous function,
# I need to change that inside the next function as well
CUEvsENTRYfiring(cuedata=toplotTC[[1]], entrydata=toplotTC[[2]], graphFolder=MixedGraphFolder,
         title="S+ onset vs. S+ entry", cueExcOnly=F, cuemin=100, cuemax=400, 
         entrymin=100, entrymax=1500)
        
CUEvsENTRYfiring(cuedata=toplotTC[[1]], entrydata=toplotTC[[3]], unitVScp=unitsDF,graphFolder=MixedGraphFolder, 
                 title="S+ onset vs. S- entry", cueExcOnly=F, cuemin=100, cuemax=400, entrymin=100, entrymax=1500)

CUEvsENTRYfiring(cuedata=toplotTC[[1]], entrydata=toplotTC[[4]], unitVScp=unitsDF,graphFolder=MixedGraphFolder,
                 title="S+ onset vs. ITI entry", cueExcOnly=F, cuemin=100, cuemax=400, entrymin=100, entrymax=1500)




########################################################################################
### BOXPLOT WITH FR RIGHT BEFORE ENTRY IN LONG-LATENCY ENTRIES
########################################################################################


PreEntry5sLat <- plotFRBoxPlotandCP_Entry(cue=c("S+", "S-"), experiment="Exp 1 2s Pre entry FR when latency over 5s from CP", 
                masterDF=list(masterDF_EntryDS, masterDF_EntryNS), graphFolder=MixedGraphFolder, 
                trialBinSize=40, WdwStart=-2000, WdwEnd=0, dataProcess="Zscores", minLat=5000,
                correctOnly=FALSE, points=F, lines=F, color=c(colindx[1], "darkblue"), 
                legLabels=c("S+", "S-"), yAxMinZ=-1, yAxMaxZ=4, yAxMaxRaw=7, capped=T, 
                capValue=c(-120, 120), cueExcOnly=F, BLneudata=allNeuronsDS)





#S+ vs. S-
DSvsNS_PreEntry <- do.call("rbind", 
        lapply(seq(1, length(unique(PreEntry5sLat[[1]]$correspBins))), function(i){
        bin <- unique(PreEntry5sLat[[1]]$correspBins)[i]
        
        datSel1 <- PreEntry5sLat[[1]][PreEntry5sLat[[1]]$correspBins==i, ]
        datSel2 <- PreEntry5sLat[[2]][PreEntry5sLat[[2]]$correspBins==i, ]
        
        comp <- wilcox.test(x=datSel1$byUnitFR, y=datSel2$byUnitFR, alternative = "greater")
        
        data.frame(W=comp$statistic, p=comp$p.value)
})
)

p.adjust(p=c(DSvsNS_PreEntry$p), method="holm")

# Comparing FR in the 2s before entry in S+ vs. S- trials in which latency to enter was 5 or more seconds.
# W           p.adj
# W   379.0 0.998981995
# W1  934.5 0.614210552
# W2 1933.0 0.028515820
# W3 2416.0 0.017735803
# W4  176.0 0.029007794
# W5  146.0 0.006298355


#across bin comp

comp1 <- PreEntry5sLat[[1]][PreEntry5sLat[[1]]$correspBins==1, ]$byUnitFR
comp2 <- PreEntry5sLat[[1]][PreEntry5sLat[[1]]$correspBins==3, ]$byUnitFR
comp3 <- PreEntry5sLat[[1]][PreEntry5sLat[[1]]$correspBins==5, ]$byUnitFR

wilcox.test(comp1, comp2, alternative="less") #W = 952, p-value = 0.001719
wilcox.test(comp2, comp3, alternative="less") #W = 256, p-value = 0.02061


#Number of units per condition (different bc less >5s entries in S+ than S-). That's why I treated these as different samples in the Wilcoxon test.
nunitsDS <- sapply(seq(1, 6), function(i){
        length(PreEntry5sLat[[1]][PreEntry5sLat[[1]]$correspBins==i, ]$byUnitFR)
})
#[1] 42 38 68 63 12 11

nunitsNS <- sapply(seq(1, 6), function(i){
        length(PreEntry5sLat[[2]][PreEntry5sLat[[2]]$correspBins==i, ]$byUnitFR)
})
#[1] 31 51 47 63 22 17



#### KAPLAN-MEIER GRAPH

