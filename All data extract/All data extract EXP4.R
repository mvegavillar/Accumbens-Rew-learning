
#############################################################
### EXPERIMENT 4: UNILATERAL INFUSIONS, SHORT ITI         ###
#############################################################

### LOAD IMPORTANT LIBRARIES
install.packages("matrixStats")
install.packages('ez')
install.packages('dplyr')

library(matrixStats)
library(ez)
library(dplyr)

setwd("E:/Dropbox (EinsteinMed)/NMDA/")

#Define colors
colindx <- c("#2171b5", "#cb181d") #Strong blue and red
colindxB <- c("#bdd7e7", "#fcae91") #Less strong blue and red
colindxC <- c("#eff3ff", "#fb6a4a") #Even less strong blue and red
colindxD <- c("#6baed6", "#fee5d9") #Lightest blue and red


### DEFINE ALL THE IMPORTANT FOLDERS

funcdirect <- paste(getwd(), "/R functions/", sep="")
datafolder <- paste(getwd(), "/EXP4_Unilateral AP5/MedPC files/", sep="")
dataForRdir <- paste(getwd(), "/EXP4_Unilateral AP5/Data for R/", sep="")
dataForRCumulative <- paste(getwd(), "/EXP4_Unilateral AP5/Data for R cumulative/", sep="")
behGraphFolder <- paste(getwd(), "/EXP4_Unilateral AP5/Graphs/Behavior/", sep="")
CPfuncFolder <- paste(getwd(), '/R functions/Change_Point-master/', sep="")
CPGraphFolder <- paste(getwd(), "/EXP4_Unilateral AP5/Graphs/Behavior/Change point/", sep="")
MixedGraphFolder <- paste(getwd(), "/EXP4_Unilateral AP5/Graphs/Mixed/", sep="")
NEXfiles <- paste(getwd(), "/EXP4_Unilateral AP5/NEX files/", sep="")
IndcumPerfGraphFolder <- paste(getwd(), "/EXP4_Unilateral AP5/Graphs/Behavior/Cumulative perf/", sep="")
PerfRelToCPFolder <- paste(getwd(), "/EXP4_Unilateral AP5/Graphs/Behavior/Beh rel to CP/", sep="")
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
load(file=paste(funcdirect, "megaplotTC.R", sep=""))
load(file=paste(funcdirect, "compareDSvsNSfromCP.R", sep=""))
load(file=paste(funcdirect, "plotFRBoxPlotandCP.R", sep=""))
load(file=paste(funcdirect, "compareVEHvsAP5fromCP.R", sep=""))
load(file=paste(getwd(), "/EXP4_Unilateral AP5/R functions/plotFRandCPhistogram_Entry.R", sep=""))
load(file=paste(getwd(), "/EXP4_Unilateral AP5/R functions/plotFRandCP_Entry.R", sep=""))


#### GENERATE IMPORTANT DATA OBJECTS
MedPCextract(MovAvg="Impinged only", funcdirect = funcdirect, datafolder = datafolder, dataForRdir = dataForRdir, dataForRCumulative=dataForRCumulative, cuelength=10)

# Load important behavior-related objects
files <- paste(dataForRdir, list.files(dataForRdir), sep="")
filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep="")
for(i in 1:length(files)){load(files[[i]])}
for(i in 1:length(filesCum)){load(filesCum[[i]])}


#
#### When in the session is the CP more likely to occur? ###
CPtrial <- CPdata$CP - (CPdata$CPsess-1)*35
sort(CPtrial[!is.na(CPtrial)])
# [1] 37 17 26 17 13 19

# data:  CPtrial[!is.na(CPtrial)] and c(20, 20, 20, 20, 35, 35, 35, 35, 35, 35, 35)
# t = -0.4539, df = 19.24, p-value = 0.655
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#         -9.175550  5.902823
# sample estimates:
#         mean of x mean of y 
# 27.90909  29.54545 

#### BEHAVIORAL GRAPHS

# RASTERS
idx[rats=="MV29"]
load(file=paste(funcdirect, "MAKERASTER.r", sep=""))
MAKERASTER(i=55, data=alldata, idxdata=csacqidx)
MAKERASTER(i=61, data=alldata, idxdata=csacqidx)
MAKERASTER(i=64, data=alldata, idxdata=csacqidx)
MAKERASTER(i=70, data=alldata, idxdata=csacqidx)
MAKERASTER(i=74, data=alldata, idxdata=csacqidx)
MAKERASTER(i=79, data=alldata, idxdata=csacqidx)

###
# I'M USING THE 10S PRE-CUE PERIOD AS BASELINE
# In 9% of the trials, the animal's head was inside the receptacle at the time the ITI window started. In those trials, I'll correct the ITI latency by assigning the MEAN ITI latency in the +/- 5 trials around that trial
###

#AVERAGE PERFORMANCE S+ VS S-. 5 trial bins.
#All rats
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



#LEARNERS ONLY.
L_Index <- !is.na(CPdata$CP)

avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(DStaskAcc[c(L_Index)], NStaskAcc[c(L_Index)]), 
             cues=c("S+", "S-"), index="Specificity", legendLocation="topleft", y_axis_label="Performance index", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(DSrespAll[c(L_Index)], NSrespAll[c(L_Index)]), 
             cues=c("S+", "S-"), index="Response ratio", legendLocation="topleft", y_axis_label="Response ratio", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(DSlatency[c(L_Index)], NSlatency[c(L_Index)]), 
             cues=c("S+", "S-"), index="Latency", legendLocation="bottomright", y_axis_label="Latency (s)", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(ITIlatency[c(L_Index)]), 
             cues=c("ITI"), index="Latency", legendLocation=-1, y_axis_label="Latency (s)", 
             behGraphFolder=behGraphFolder, plot=T)


#NON LEARNERS ONLY.
NL_Index <- is.na(CPdata$CP)

avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(DStaskAcc[c(NL_Index)], NStaskAcc[c(NL_Index)]), 
             cues=c("S+", "S-"), index="Specificity", legendLocation="topleft", y_axis_label="Performance index", 
             behGraphFolder=behGraphFolder, plot=T, YMIN=-2, YMAX=5)

avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(DSrespAll[c(NL_Index)], NSrespAll[c(NL_Index)]), 
             cues=c("S+", "S-"), index="Response ratio", legendLocation="topleft", y_axis_label="Response ratio", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(DSlatency[c(NL_Index)], NSlatency[c(NL_Index)]), 
             cues=c("S+", "S-"), index="Latency", legendLocation="bottomright", y_axis_label="Latency (s)", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), data=list(ITIlatency[c(NL_Index)]), 
             cues=c("ITI"), index="Latency", legendLocation=-1, y_axis_label="Latency (s)", 
             behGraphFolder=behGraphFolder, plot=T)



### CHANGE POINT: S+ SPECIFICITY
# Change point analysis using the method in Gallistel et al., 2004. I'll use S+ SPECIFICITY as the main performance index. But I also need to look at the results that other performance indexes give me.

CPdata <- CPextract(GallCrit=1.3, plot=T, CPfuncFolder=CPfuncFolder, idx=idx, CPGraphFolder=CPGraphFolder, dataForRdir=dataForRdir, dataForRCumulative=dataForRCumulative)


CPdataExp4 <- CPdata
alldataExp4 <- alldata
csacqidxExp4 <- csacqidx
ratsExp4 <- rats
idxExp4 <- idx

save(csacqidxExp4, file=paste(dataForRdir, "csacqidxExp4.ridx", sep=""))
save(alldataExp4, file=paste(dataForRdir, "alldataExp4.rdat", sep=""))
save(ratsExp4, file=paste(dataForRdir, "ratsExp4.rdat", sep=""))
save(idxExp4, file=paste(dataForRdir, "idxExp4.rdat", sep=""))

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

###NON LEARNERS:TOTAL REWARDS
NL_CPdata <- CPdata[is.na(CPdata$CP), ]
plot.new()
plot.window(xlim=c(0, 2), ylim=c(0, 210))
rect(xleft=0, xright=1, ybottom=0, ytop = mean(NL_CPdata$TotalRewards), col="gray30")
points(x=rep(0.5, nrow(NL_CPdata)), y=NL_CPdata$TotalRewards, pch=19, col="black")
rect(xleft=1, xright=2, ybottom=0, ytop=mean(NL_CPdata$TotalRewardsShortLat, col="white"))
points(x=rep(1.5, nrow(NL_CPdata)), y=NL_CPdata$TotalRewardsShortLat, pch=19, col="black")

for(i in 1:nrow(NL_CPdata)){
        lines(x=c(0.5, 1.5), y=c(NL_CPdata$TotalRewards[i], NL_CPdata$TotalRewardsShortLat[i]))
}


axis(side=2, las=2)
mtext(side=2, text="Total # rewards", line=2.5, font=2, cex=1.4)
axis(side=1, at=c(0.5, 1.5), labels=c("Total rewards", "# rewards with < 2 latency"), font=2, cex.axis=1.2)


#LEARNERS: TOTAL TRIALS TO CP VS. NUMBER OF REWARDS
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

save(CPframeExp4, file=paste(dataForRdir, "CPframeExp3.rdat", sep="")) #CP values for different behavioral indexes
save(CPdataExp4, file=paste(dataForRdir, "CPdataExp3.rdat", sep="")) #CP values and other details for PERFORMANCE INDEX (S+ specificity)


#CUMULATIVE INDIVIDUAL PERFORMANCE
#Smart rats
cumulativeIndGraphs(numSess=6, numCues=35, sessLines=F, smartRat = T, dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=IndcumPerfGraphFolder, imgFormat="pdf", limit=210)

#Dumb rats
cumulativeIndGraphs(numSess=6, numCues=35, sessLines=F, smartRat = F, dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=IndcumPerfGraphFolder, imgFormat="pdf", limit=210)

#All rats
cumulativeIndGraphs(numSess=6, numCues=35, sessLines=F, smartRat = "NA", dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=IndcumPerfGraphFolder, imgFormat="pdf", limit=210)


#MEAN PERFORMANCE AS A FUNCTION OF DISTANCE FROM CP
PerformanceFromCP(relTrialMin=-100, relTrialMax=100, trialBinSize=5, typegraph="both", imgFormat="pdf", 
                  dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=PerfRelToCPFolder, CPdata=CPdata, 
                  idx=idx, rats = rats, alldata = alldata)

PerformanceFromCP(relTrialMin=-100, relTrialMax=100, trialBinSize=5, typegraph="both", imgFormat="pdf", 
                  dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=PerfRelToCPFolder, CPdata=CPdata, 
                  idx=idx, rats = rats, alldata = alldata)


#PRE CP VS. POST CP PERFORMANCE
PrevsPostFolder <- paste(PerfRelToCPFolder, "Average Pre vs Post/", sep="")
PrePostCP_Perf(data=DSrespAll, CP=CPdata$CP, y_axis_label="S+ response ratio", graphFolder=PrevsPostFolder, plot=T) #t = -5.0607, df = 10, p-value = 0.0004914 No correction
PrePostCP_Perf(data=NSrespAll, CP=CPdata$CP, y_axis_label="S- response ratio", graphFolder=PrevsPostFolder, plot=T) #t = 0.82127, df = 10, p-value =  0.4306
PrePostCP_Perf(data=DSlatency, CP=CPdata$CP, y_axis_label="S+ latency", graphFolder=PrevsPostFolder, plot=T) #t = 5.938, df = 10, p-value = 0.0001435
PrePostCP_Perf(data=NSlatency, CP=CPdata$CP, y_axis_label="S- latency", graphFolder=PrevsPostFolder, plot=T) #t = -0.77407, df = 10, p-value = 0.4568
PrePostCP_Perf(data=DStaskAcc, CP=CPdata$CP, y_axis_label="S+ specificity", graphFolder=PrevsPostFolder, plot=T) #t = -10.21, df = 10, p-value = 1.313e-06
PrePostCP_Perf(data=NStaskAcc, CP=CPdata$CP, y_axis_label="S- specificity", graphFolder=PrevsPostFolder, plot=T) #t = 0.42843,, df = 10, p-value =  0.6774
PrePostCP_Perf(data=ITIlatency, CP=CPdata$CP, y_axis_label="ITI latency", graphFolder=PrevsPostFolder, plot=T) #t = 1.2465, df = 10, p-value = 0.241



p.adjust(c(0.0004914, 0.0001435, 1.313e-06,  0.241), method="holm") #9.828e-04 4.305e-04 5.252e-06 2.410e-01

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
allNeuronsDS <-neuralhist(funcdirect=funcdirect, path=NEXfiles, event=1, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="both")
allNeuronsDS_VEH = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=1, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsNS_VEH = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=2, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsDSresponded_VEH = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=5, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsDSmissed_VEH = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=6, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsNSresponded_VEH = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=7, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsNSmissed_VEH = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=8, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsEntryDS_VEH = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle", BLwdw = 5, targetwdw = 1.5)
allNeuronsEntryNS_VEH = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle", BLwdw = 5, targetwdw = 1.5)
allNeuronsEntryITI_VEH = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle", BLwdw = 5, targetwdw = 1.5)


#### RUN THIS FOR NEURONS IN AP5 SIDE
allNeuronsDS_AP5 = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=1, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsNS_AP5 = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=2, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsDSresponded_AP5 = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=5, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsDSmissed_AP5 = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=6, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsNSresponded_AP5 = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=7, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsNSmissed_AP5 = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=8, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsEntryDS_AP5 = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug", BLwdw = 5, targetwdw = 1.5)
allNeuronsEntryNS_AP5 = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug", BLwdw = 5, targetwdw = 1.5)
allNeuronsEntryITI_AP5 = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug", BLwdw = 5, targetwdw = 1.5)


#sAVE THESE OBJECTS
save(allNeuronsDS_VEH, file=paste(dataForRdir, 'allNeuronsDS_VEH.rdat', sep=""))
save(allNeuronsNS_VEH, file=paste(dataForRdir, 'allNeuronsNS_VEH.rdat', sep=""))
save(allNeuronsDSresponded_VEH, file=paste(dataForRdir, 'allNeuronsDSresponded_VEH.rdat', sep=""))
save(allNeuronsDSmissed_VEH, file=paste(dataForRdir, 'allNeuronsDSmissed_VEH.rdat', sep=""))
save(allNeuronsNSresponded_VEH, file=paste(dataForRdir, 'allNeuronsNSresponded_VEH.rdat', sep=""))
save(allNeuronsNSmissed_VEH, file=paste(dataForRdir, 'allNeuronsNSmissed_VEH.rdat', sep=""))
save(allNeuronsEntryDS_VEH, file=paste(dataForRdir, 'allNeuronsEntryDS_VEH.rdat', sep=""))
save(allNeuronsEntryNS_VEH, file=paste(dataForRdir, 'allNeuronsEntryNS_VEH.rdat', sep=""))
save(allNeuronsEntryITI_VEH, file=paste(dataForRdir, 'allNeuronsEntryITI_VEH.rdat', sep=""))
save(allNeuronsDS_AP5, file=paste(dataForRdir, 'allNeuronsDS_AP5.rdat', sep=""))
save(allNeuronsNS_AP5, file=paste(dataForRdir, 'allNeuronsNS_AP5.rdat', sep=""))
save(allNeuronsDSresponded_AP5, file=paste(dataForRdir, 'allNeuronsDSresponded_AP5.rdat', sep=""))
save(allNeuronsDSmissed_AP5, file=paste(dataForRdir, 'allNeuronsDSmissed_AP5.rdat', sep=""))
save(allNeuronsNSresponded_AP5, file=paste(dataForRdir, 'allNeuronsNSresponded_AP5.rdat', sep=""))
save(allNeuronsNSmissed_AP5, file=paste(dataForRdir, 'allNeuronsNSmissed_AP5.rdat', sep=""))
save(allNeuronsEntryDS_AP5, file=paste(dataForRdir, 'allNeuronsEntryDS_AP5.rdat', sep=""))
save(allNeuronsEntryNS_AP5, file=paste(dataForRdir, 'allNeuronsEntryNS_AP5.rdat', sep=""))
save(allNeuronsEntryITI_AP5, file=paste(dataForRdir, 'allNeuronsEntryITI_AP5.rdat', sep=""))

#Load files again      
#files <- paste(dataForRdir, list.files(dataForRdir), sep=""); for(i in 1:length(files)){load(files[[i]])}

#BUILD A DATAFRAME WITH FIRING RATE OF EACH NEURON ON EACH TRIAL WITH BEH INFO ABOUT EACH TRIAL
     
#VEH side
masterDF_DS_VEH <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS_VEH, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, 
                                      BLneudata=allNeuronsDS_VEH)

#By modality
masterDF_DS_VEH_Tone <- filter(masterDF_DS_VEH, masterDF_DS_VEH$modality==98)
masterDF_DS_VEH_Light <- filter(masterDF_DS_VEH, masterDF_DS_VEH$modality==99)



masterDF_DSresponded_VEH <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSresponded_VEH, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDSresponded_VEH)
masterDF_DSmissed_VEH <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSmissed_VEH, 
                                               CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                               BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDSmissed_VEH)
masterDF_NS_VEH <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_VEH, 
                                            CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                            BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_VEH)
masterDF_NSresponded_VEH <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_VEH, 
                                      CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_VEH)
masterDF_NSmissed_VEH <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_VEH, 
                                      CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_VEH)


#The "entry-related" masterDF don't give me FR info of every trial, only of those trials in which the animal responded
masterDF_DSEntry_VEH <- FRbyNEURONbyBINcue(eventType="entry", cue="S+", neudata=allNeuronsEntryDS_VEH, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDS_VEH)
masterDF_NSEntry_VEH <- FRbyNEURONbyBINcue(eventType="entry", cue="S+", neudata=allNeuronsEntryNS_VEH, 
                                           CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                           BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_VEH)
masterDF_ITIEntry_VEH <- FRbyNEURONbyBINcue(eventType="entry", cue="S+", neudata=allNeuronsEntryITI_VEH, 
                                           CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                           BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDS_VEH)


#AP5 side
masterDF_DS_AP5 <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS_AP5, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDS_AP5)


#By modality
masterDF_DS_AP5_Tone <- filter(masterDF_DS_AP5, masterDF_DS_AP5$modality==98)
masterDF_DS_AP5_Light <- filter(masterDF_DS_AP5, masterDF_DS_AP5$modality==99)


masterDF_DSresponded_AP5 <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSresponded_AP5, 
                                               CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                               BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDSresponded_AP5)
masterDF_DSmissed_AP5 <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSmissed_AP5, 
                                            CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                            BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDSmissed_AP5)
masterDF_NS_AP5 <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_AP5, 
                                      CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_AP5)
masterDF_NSresponded_AP5 <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_AP5, 
                                               CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                               BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_AP5)
masterDF_NSmissed_AP5 <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_AP5, 
                                            CPvector=CPdata$NS.PreCP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                            BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_AP5)

masterDF_DSEntry_AP5 <- FRbyNEURONbyBINcue(eventType="entry", cue="S+", neudata=allNeuronsEntryDS_AP5, 
                                           CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                           BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDS_AP5)
masterDF_NSEntry_AP5 <- FRbyNEURONbyBINcue(eventType="entry", cue="S+", neudata=allNeuronsEntryNS_AP5, 
                                           CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                           BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsNS_AP5)
masterDF_ITIEntry_AP5 <- FRbyNEURONbyBINcue(eventType="entry", cue="S+", neudata=allNeuronsEntryITI_AP5, 
                                            CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                            BLduration=2, sessionCPperRat=CPdata$CPsess, BLneudata=allNeuronsDS_AP5)



save(masterDF_DS_VEH, file=paste(dataForRdir, "masterDF_DS_VEH.rdat", sep=""))
save(masterDF_DSresponded_VEH, file=paste(dataForRdir, "masterDF_DSresponded_VEH.rdat", sep=""))
save(masterDF_DSmissed_VEH, file=paste(dataForRdir, "masterDF_DSmissed_VEH.rdat", sep=""))
save(masterDF_NS_VEH, file=paste(dataForRdir, "masterDF_NS_VEH.rdat", sep=""))
save(masterDF_NSresponded_VEH, file=paste(dataForRdir, "masterDF_NSresponded_VEH.rdat", sep=""))
save(masterDF_NSmissed_VEH, file=paste(dataForRdir, "masterDF_NSmissed_VEH.rdat", sep=""))
save(masterDF_DSEntry_VEH, file=paste(dataForRdir, "masterDF_DSEntry_VEH.rdat", sep=""))
save(masterDF_NSEntry_VEH, file=paste(dataForRdir, "masterDF_NSEntry_VEH.rdat", sep=""))
save(masterDF_ITIEntry_VEH, file=paste(dataForRdir, "masterDF_ITIEntry_VEH.rdat", sep=""))

save(masterDF_DS_AP5, file=paste(dataForRdir, "masterDF_DS_AP5.rdat", sep=""))
save(masterDF_DSresponded_AP5, file=paste(dataForRdir, "masterDF_DSresponded_AP5.rdat", sep=""))
save(masterDF_DSmissed_AP5, file=paste(dataForRdir, "masterDF_DSmissed_AP5.rdat", sep=""))
save(masterDF_NS_AP5, file=paste(dataForRdir, "masterDF_NS_AP5.rdat", sep=""))
save(masterDF_NSresponded_AP5, file=paste(dataForRdir, "masterDF_NSresponded_AP5.rdat", sep=""))
save(masterDF_NSmissed_AP5, file=paste(dataForRdir, "masterDF_NSmissed_AP5.rdat", sep=""))
save(masterDF_DSEntry_AP5, file=paste(dataForRdir, "masterDF_DSEntry_AP5.rdat", sep=""))
save(masterDF_NSEntry_AP5, file=paste(dataForRdir, "masterDF_NSEntry_AP5.rdat", sep=""))
save(masterDF_ITIEntry_AP5, file=paste(dataForRdir, "masterDF_ITIEntry_AP5.rdat", sep=""))

#masterDFsumm <- masterDFsummary(masterDF=masterDF_DS, neudata=allNeuronsDS, format="Zscores", postCueWdw=c(0, 0.4), postSpikeWdw=c(0.4, 5), orEntry=T, cuelength=10)




####################33
### CHECK THAT THE CRITERIA FOR CUE EXCITED UNITS MAKES SENSE
#Time course of excitation
toplotTC_VEH <- megaplotTC(experiment="Exp 4 VEH 7z", data=list(allNeuronsDS_VEH), graphmin=500, 
                           graphmax=3000, cuewinmin=100, cuewinmax=400, colpalette="Rainbow", 
                           minFR=-2, maxFR=7, graphFolder=MixedGraphFolder, 
                           ZcolLabels=c("ZDS"), arrangeBy=c("ZDS"))

toplotTC_AP5 <- megaplotTC(experiment="Exp 4 AP5 7z", data=list(allNeuronsDS_AP5), graphmin=500, 
                           graphmax=3000, cuewinmin=100, cuewinmax=400,
                           colpalette="Rainbow", minFR=-2, maxFR=7, graphFolder=MixedGraphFolder, 
                           ZcolLabels=c("ZDS"), arrangeBy=c("ZDS"))

toplotTC_VEH_DSEntry <- megaplotTC(experiment="Exp 4 VEH DS Entry", data=list(allNeuronsEntryDS_VEH), graphmin=500, 
                           graphmax=4000, cuewinmin=100, cuewinmax=400, colpalette="Rainbow", 
                           minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, 
                           ZcolLabels=c("ZDS"), arrangeBy=c("ZDS"))

toplotTC_VEH_NSEntry <- megaplotTC(experiment="Exp 4 VEH NS Entry", data=list(allNeuronsEntryNS_VEH), graphmin=500, 
                                   graphmax=4000, cuewinmin=100, cuewinmax=400, colpalette="Rainbow", 
                                   minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, 
                                   ZcolLabels=c("ZDS"), arrangeBy=c("ZDS"))

toplotTC_VEH_ITIEntry <- megaplotTC(experiment="Exp 4 VEH ITI Entry", data=list(allNeuronsEntryITI_VEH), graphmin=500, 
                                   graphmax=4000, cuewinmin=100, cuewinmax=400, colpalette="Rainbow", 
                                   minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, 
                                   ZcolLabels=c("ZDS"), arrangeBy=c("ZDS"))

toplotTC_AP5_DSEntry <- megaplotTC(experiment="Exp 4 AP5 DS Entry", data=list(allNeuronsEntryDS_AP5), graphmin=500, 
                                   graphmax=4000, cuewinmin=100, cuewinmax=400, colpalette="Rainbow", 
                                   minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, 
                                   ZcolLabels=c("ZDS"), arrangeBy=c("ZDS"))

toplotTC_AP5_NSEntry <- megaplotTC(experiment="Exp 4 AP5 NS Entry", data=list(allNeuronsEntryNS_AP5), graphmin=500, 
                                   graphmax=4000, cuewinmin=100, cuewinmax=400, colpalette="Rainbow", 
                                   minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, 
                                   ZcolLabels=c("ZDS"), arrangeBy=c("ZDS"))

toplotTC_AP5_ITIEntry <- megaplotTC(experiment="Exp 4 AP5 ITI Entry", data=list(allNeuronsEntryITI_AP5), graphmin=500, 
                                    graphmax=4000, cuewinmin=100, cuewinmax=400, colpalette="Rainbow", 
                                    minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, 
                                    ZcolLabels=c("ZDS"), arrangeBy=c("ZDS"))



#Histograms with percentage of cue excited units per bin
CueExcCheck(experiment="Exp 4 VEH", data=toplotTC_VEH, graphFolder=MixedGraphFolder, 
            autom.FRrange=T, FRbin=0.3, cexEx=0.8, cexInh =0.35)

CueExcCheck(experiment="Exp 4 AP5", data=toplotTC_AP5, graphFolder=MixedGraphFolder, 
            autom.FRrange=T, FRbin=0.3, cexEx=0.8, cexInh =0.35)

### PLOT FR POST-CUE (400MS WINDOW) as a function of distance to change point

#Redefine colors
colindx <- c("#2171b5", "#cb181d") #Strong blue and red

#CUE
plotFRandCP(experiment="Exp 4 VEH side", cue=c("S+", "S-"), masterDF=list(masterDF_DS_VEH, masterDF_NS_VEH), graphFolder=MixedGraphFolder, 
            trialBinSize=5, WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
            colindx =c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), capped=T, capValue = c(-105, 105), 
            yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_VEH)

plotFRandCP(experiment="Exp 4 AP5 side", cue=c("S+", "S-"), masterDF=list(masterDF_DS_AP5, masterDF_NS_AP5), graphFolder=MixedGraphFolder, 
            trialBinSize=5, WdwStart=100, WdwEnd=400, dataProcess="Zscores", legLabels=c("S+", "S-"),
            correctOnly=FALSE, colindx=c(colindx[2], "darkred"), capped=T, capValue = c(-105, 105), 
            yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_AP5)

plotFRandCP(experiment="Exp 4 BOTH SIDES Cue exc only", cue=c("S+", "S-", "S+", "S-"), masterDF=list(masterDF_DS_VEH, masterDF_NS_VEH, masterDF_DS_AP5, masterDF_NS_AP5), 
            graphFolder=MixedGraphFolder, trialBinSize=5, WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
            colindx =c("#2171b5", "darkblue", "#cb181d", "darkred"), legLabels=c("S+ VEH", "S- VEH", "S+ AP5", "S- AP5"), 
            capped=T, capValue = c(-105, 105), cueExcOnly = T,
            yAxMinZ = -1, yAxMaxZ = 10, yAxMaxRaw = 10, neudata=allNeuronsDS_VEH)

plotFRandCP(experiment="Exp 4 both sides", cue=c("S+", "S+"), masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
            graphFolder=MixedGraphFolder, trialBinSize=5, WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
            colindx =c(colindx[1],colindx[2]), legLabels=c("S+ VEH", "S+ AP5"), capped=T, capValue = c(-90, 90), 
            yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_VEH)




#### 35 TRIAL BINS: PERCENTAGE OF CUE-EVOKED EXCITATIONS
#% Exc. The function doesn't work well, I have to run the parameters first and then go to the function and run its contents. Fix this.
PercCueExc_VEH <-  plotFRandCP(experiment="Exp 4 VEH side", cue=c("S+"), masterDF=list(masterDF_DS_VEH), 
                               graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=100, WdwEnd=400, 
                               dataProcess="PercCueExc", correctOnly=FALSE, colindx =c(colindx[1]), 
                               legLabels=c("S+"), capped=T, capValue = c(-105, 105), yAxMinZ = -1, yAxMaxZ = 8, 
                               yAxMaxRaw = 10, neudata=allNeuronsDS_VEH)

#bin   trialBins CueEx notCueExc CueInh notCueInh   N
#1   1      -105    19        56     29        46  75
#2   2       -70    45        92     36       101 137
#3   3       -35    63        66     22       107 129
#4   4         0    87        44     21       110 131
#5   5        35    68        24     13        79  92
#6   6        70    23        13      5        31  36


PercCueExc_AP5 <- plotFRandCP(experiment="Exp 4 AP5 side2", cue=c("S+"), masterDF=list(masterDF_DS_AP5), graphFolder=MixedGraphFolder, 
            trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, 
            colindx =c(colindx[2]), legLabels=c("S+"), capped=T, capValue = c(-105, 105), 
            yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_AP5)
# bin   trialBins CueEx notCueExc CueInh notCueInh   N
# 1   1      -105    14        52     23        43  66
# 2   2       -70    24        79     31        72 103
# 3   3       -35    31        65     22        74  96
# 4   4         0    43        59     24        78 102
# 5   5        35    37        48     19        66  85
# 6   6        70    36        39     17        58  75
# 

PercCueExc_VEH <- PercCueExc_VEH[[1]]
PercCueExc_AP5 <- PercCueExc_AP5[[1]]

save(PercCueExc_VEH, file=paste(dataForRdir, "PercCueExc_VEH.rdat", sep=""))
save(PercCueExc_AP5, file=paste(dataForRdir, "PercCueExc_AP5.rdat", sep=""))

#Is the proportion of cue-exc neurons modulated by a) the drug; b) the amount of training? (run a different chi.sq test for bins 1, 3 and 5 and for bins 2, 4 and 6, that way the units on each cell will be independent for sure)
#For Chi-sq analysis. I need to make sure that each cell has a completely different population of neurons. Bc each session has 40 trials of each kind, by only comparing bins that are separated by 40 trials I can guarantee that.


### VEH SIDE
contingency_table_comp1<- PercCueExc_VEH[is.even(PercCueExc_VEH$bin)==FALSE, c(2, 3, 4)]
contingency_table_comp2 <-  PercCueExc_VEH[is.even(PercCueExc_VEH$bin)==TRUE, c(2, 3, 4)]

inhCont_table_comp1 <- PercCueExc_VEH[is.even(PercCueExc_VEH$bin)==FALSE, c(2, 5, 6)]
inhCont_table_comp2 <- PercCueExc_VEH[is.even(PercCueExc_VEH$bin)==TRUE, c(2, 5, 6)]

#Chi-sq analysis excitations:
#Bins 1, 3 and 5
critData <- contingency_table_comp1[, -1]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 39.32, df = 2, p-value = 2.896e-09

fisher.test(critData)
# data:  critData
# p-value = 1.4e-09
# alternative hypothesis: two.sided



#Bins 2, 4, 6, and 8
critData <- contingency_table_comp2[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 32.911, df = 2, p-value = 7.137e-08

fisher.test(critData)
# data:  critData
# p-value = 5.577e-08
# alternative hypothesis: two.sided



#Post-hoc comparisons for % cue-excited units
# Is the distribution of exc/not exc neurons different from the first (-120 to -80) and the third (-40 to 0) bin. 
contTableExc_VEH <- PercCueExc_VEH[ , c(3, 4)]

#Bin 1 vs. 3
fisher.test(contTableExc_VEH[c(1, 3), ], alternative="less") #95% CI: 0.000000  0.6268223 Odds ratio=0.3572457 p=0.0006973
#Bin 3 vs. 5
fisher.test(contTableExc_VEH[c(3, 5), ], alternative="less") #95% CI:0.0000000 0.5697446 odds ratio=0.3385928 p= 0.0001357

#Bin 2 vs. 4
fisher.test((contTableExc_VEH[c(2, 4), ]), alternative="less") #95% CI: 0.0000000 2.270799 odds ratio= 0.2487772  p= 2.979e-08
#Bin 4 vs. 6
fisher.test((contTableExc_VEH[c(4, 6), ]), alternative="less") #95% CI:  0.000000 1.736326 odds ratio= 1.116834 p=0.6879

#Correct p value all comparisons
p.adjust(p=c(0.0006973, 0.0001357, 2.979e-08, 0.6879), method="holm")
# 1.3946e-03 4.0710e-04 1.1916e-07 6.8790e-01

#Chi-sq for inhibitions:
#Bins 1, 3 and 5
critData <- inhCont_table_comp1[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 17.492, df = 2, p-value = 0.0001591

fisher.test(critData) #p-value = 0.000287

#Bins 2, 4 and 6
critData <- inhCont_table_comp2[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 5.3955, df = 2, p-value = 0.06736
# 


fisher.test(x=critData) #p-value = 0.07836



contTableInh_VEH <- PercCueExc_VEH[ , c(5, 6)]

#Bin 1 vs. 3
fisher.test(contTableInh_VEH[c(1, 3), ], alternative="greater") #95% CI:1.679992      Inf Odds ratio=3.048006 p=0.0006165
#Bin 3 vs. 5
fisher.test(contTableInh_VEH[c(3, 5), ], alternative="greater") #95% CI:0.6303721       Inf odds ratio= 1.248218  p-value = 0.3471



### AP5 SIDE
contingency_table_comp1<- PercCueExc_AP5[is.even(PercCueExc_AP5$bin)==FALSE, c(2, 3, 4)]
contingency_table_comp2 <-  PercCueExc_AP5[is.even(PercCueExc_AP5$bin)==TRUE, c(2, 3, 4)]

inhCont_table_comp1 <- PercCueExc_AP5[is.even(PercCueExc_AP5$bin)==FALSE, c(2, 5, 6)]
inhCont_table_comp2 <- PercCueExc_AP5[is.even(PercCueExc_AP5$bin)==TRUE, c(2, 5, 6)]

#Chi-sq analysis excitations:
#Bins 1, 3 and 5
critData <- contingency_table_comp1[, -1]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 8.402, df = 2, p-value = 0.01498

fisher.test(critData) #p-value = 0.01457

#Bins 2, 4, 6, and 8
critData <- contingency_table_comp2[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 13.376, df = 2, p-value = 0.001246

fisher.test(x=critData) #p-value = 0.001003

p.adjust(p=c(0.01498, 0.001246), method="holm") # 0.014980 0.002492


#Post-hoc comparisons for % cue-excited units
# Is the distribution of exc/not exc neurons different from the first (-120 to -80) and the third (-40 to 0) bin. 
contTableExc_AP5 <- PercCueExc_AP5[ , c(3, 4)]

#Bin 1 vs. 3
fisher.test(contTableExc_AP5[c(1, 3), ], alternative="less") #95% CI: 0.000000 1.101231 Odds ratio=0.5664783 p= 0.08469
#Bin 3 vs. 5
fisher.test(contTableExc_AP5[c(3, 5), ], alternative="less") #95% CI:0.0000000 1.076245 odds ratio=0.6203911 p= 0.08012

#Bin 2 vs. 4
fisher.test((contTableExc_AP5[c(2, 4), ]), alternative="less") #95% CI: 0.0000000 0.7217133 odds ratio= 0.4186572  p= 0.003075
#Bin 4 vs. 6
fisher.test((contTableExc_AP5[c(4, 6), ]), alternative="less") #95% CI:  0.000000 1.36554 odds ratio=0.7906157 p=0.2676

#Correct p value all comparisons
p.adjust(p=c( 0.08469, 0.08012, 0.003075, 0.2676), method="holm")
# 0.24036 0.24036 0.01230 0.26760





#Chi-sq for inhibitions:
#Bins 1, 3 and 5
critData <- inhCont_table_comp1[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 3.7551, df = 2, p-value = 0.153

fisher.test(critData) #p-value = 0.1561


#Bins 2, 4 and 6
critData <- inhCont_table_comp2[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 1.6553, df = 2, p-value =  0.4371
# 

fisher.test(x=critData) #p-value = 0.4518

contTableInh_AP5 <- PercCueExc_AP5[ , c(5, 6)]


#Bin 1 vs. 3
fisher.test(contTableInh_AP5[c(1, 3), ], alternative="greater") #95% CI: 0.9447521       Inf Odds ratio=1.792487 p= 0.06897
#Bin 3 vs. 5
fisher.test(contTableInh_AP5[c(3, 5), ], alternative="greater") #95% CI:0.5425087       Inf odds ratio=1.032542 p= 0.5356



#########

### % CUE EXCITED
### VEH VS AP5
#Bin 1: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH[1, ], contTableExc_AP5[1, ]), alternative="less") #95% CI: 0.000000 2.644727  Odds ratio=1.258141 p= 0.7807

#Bin 2: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH[2, ], contTableExc_AP5[2, ]), alternative="greater") #95% CI: 0.000000 inf Odds ratio=1.606878 p= 0.06973

#Bin 3: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH[3, ], contTableExc_AP5[3, ]), alternative="greater") #95% CI: 0.000000 inf Odds ratio=1.995273 p= 0.009084

#Bin 4: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH[4, ], contTableExc_AP5[4, ]), alternative="greater") #95% CI: 0.000000 inf Odds ratio=2.700761  p= 0.0001753

#Bin 5: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH[5, ], contTableExc_AP5[5, ]), alternative="greater") #95% CI: 0.000000 inf Odds ratio=3.646766 p=  3.379e-05

#Bin 6: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH[6, ], contTableExc_AP5[6, ]), alternative="greater") #95% CI: 0.000000 inf Odds ratio=2.629861 p= 0.01534

giveStars(p.adjust(p=c(0.7807, 0.06973, 0.009084, 0.0001753, 3.379e-05, 0.01534), method="holm"))
# 0.78070000 0.13946  0.03633600 0.00087650 0.00020274 0.04602000


### % CUE INHIBITED
### VEH VS AP5
#Bin 1: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH[1, ], contTableInh_AP5[1, ]), alternative="greater") #95% CI: 0.000000 2.644727  Odds ratio=1.258141 p= 0.3848

#Bin 2: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH[2, ], contTableInh_AP5[2, ]), alternative="less") #95% CI: 0.000000 1.389301 Odds ratio=0.8285122 p= 0.3052

#Bin 3: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH[3, ], contTableInh_AP5[3, ]), alternative="less") #95% CI: 0.000000 1.273573 Odds ratio=0.6927576 p= 0.1769

#Bin 4: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH[4, ], contTableInh_AP5[4, ]), alternative="less") #95% CI: 0.000000 0.1021 Odds ratio=0.6217663  p= 0.1021

#Bin 5: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH[5, ], contTableInh_AP5[5, ]), alternative="less") #95% CI: 0.000000 1.176476 Odds ratio=0.5734464 p=  0.1103

#Bin 6: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH[6, ], contTableInh_AP5[6, ]), alternative="less") #95% CI: 0.000000 1.50749 Odds ratio=0.553058 p= 0.2048

p.adjust(p=c(0.3848, 0.3052, 0.1769, 0.1021, 0.1103, 0.2048), method="holm")
#0.7076 0.7076 0.7076 0.6126 0.6126 0.7076




#Correct p value of ALL comparisons:
VEHcomps.pval <- c(0.0006973, 0.0001357, 2.979e-08, 0.6879)
AP5comps.pval <-  c(0.08469, 0.08012, 0.003075, 0.2676)
VEHvsAP5compsEXC.pval <- c(0.7807, 0.06973, 0.009084, 0.0001753, 3.379e-05, 0.01534)
VEHvsAP5compsINH.pval <- c(0.3848, 0.3052, 0.1769, 0.1021, 0.1103, 0.2048)

p.adjust(p=c(VEHcomps.pval, AP5comps.pval, VEHvsAP5compsEXC.pval, VEHvsAP5compsINH.pval), method="holm")

VEHcomps.pval.adj <- giveStars(c( 1.11568e-02, 2.44260e-03, 5.95800e-07, 1.00000e+00))
AP5comps.pval.adj <- giveStars(c( 8.81320e-01, 8.81320e-01, 4.61250e-02, 1.00000e+00))
VEHvsAP5compsEXC.pval.adj <- giveStars(c(1.00000e+00, 8.36760e-01, 1.27176e-01, 2.98010e-03, 6.42010e-04, 1.99420e-01))
VEHvsAP5compsINH.pval.adj <- giveStars(c(1.00000e+00, 1.00000e+00, 1.00000e+00, 9.18900e-01, 9.18900e-01, 1.00000e+00))







###########################
#### 20 TRIAL BINS: PERCENTAGE OF CUE-EVOKED EXCITATIONS
#% Exc. The function doesn't work well, I have to run the parameters first and then go to the function and run its contents. Fix this.
PercCueExc_VEH_20 <-  plotFRandCP(experiment="Exp 4 VEH side", cue=c("S+"), masterDF=list(masterDF_DS_VEH), 
                               graphFolder=MixedGraphFolder, trialBinSize=20, WdwStart=100, WdwEnd=400, 
                               dataProcess="PercCueExc", correctOnly=FALSE, colindx =c(colindx[1]), 
                               legLabels=c("S+"), capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 8, 
                               yAxMaxRaw = 10, neudata=allNeuronsDS_VEH)
# [[1]]
# bin trialBins CueEx notCueExc CueInh notCueInh
# 1    1      -100    18        54     27        45
# 2    2       -80    41        87     35        93
# 3    3       -60    29        43     11        61
# 4    4       -40    63        63     21       105
# 5    5       -20    38        31     13        56
# 6    6         0    86        44     20       110
# 7    7        20    53        21     11        63
# 8    8        40    56        17      9        64
# 9    9        60    23        12      6        29
# 10  10        80     8         8      2        14


PercCueExc_AP5_20 <- plotFRandCP(experiment="Exp 4 AP5 side 20", cue=c("S+"), masterDF=list(masterDF_DS_AP5), graphFolder=MixedGraphFolder, 
                              trialBinSize=20, WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, 
                              colindx =c(colindx[2]), legLabels=c("S+"), capped=T, capValue = c(-100, 100), 
                               yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_AP5)
# [[1]]
# bin trialBins CueEx notCueExc CueInh notCueInh
# 1    1      -100    13        43     16        40
# 2    2       -80    23        74     30        67
# 3    3       -60    16        49     20        45
# 4    4       -40    30        69     26        73
# 5    5       -20    20        32     11        41
# 6    6         0    43        58     24        77
# 7    7        20    23        32     13        42
# 8    8        40    34        46     18        62
# 9    9        60    30        39     17        52
# 10  10        80    22        23     11        34
# # 

PercCueExc_VEH_20 <- PercCueExc_VEH_20[[1]]
PercCueExc_AP5_20 <- PercCueExc_AP5_20[[1]]

save(PercCueExc_VEH_20, file=paste(dataForRdir, "PercCueExc_VEH.rdat", sep=""))
save(PercCueExc_AP5_20, file=paste(dataForRdir, "PercCueExc_AP5.rdat", sep=""))

#Is the proportion of cue-exc neurons modulated by a) the drug; b) the amount of training? (run a different chi.sq test for bins 1, 3 and 5 and for bins 2, 4 and 6, that way the units on each cell will be independent for sure)
#For Chi-sq analysis. I need to make sure that each cell has a completely different population of neurons. Bc each session has 40 trials of each kind, by only comparing bins that are separated by 40 trials I can guarantee that.


### VEH SIDE
contingency_table_comp1<- PercCueExc_VEH_20[is.even(PercCueExc_VEH_20$bin)==FALSE, c(2, 3, 4)]
contingency_table_comp2 <-  PercCueExc_VEH_20[is.even(PercCueExc_VEH_20$bin)==TRUE, c(2, 3, 4)]

inhCont_table_comp1 <- PercCueExc_VEH_20[is.even(PercCueExc_VEH_20$bin)==FALSE, c(2, 5, 6)]
inhCont_table_comp2 <- PercCueExc_VEH_20[is.even(PercCueExc_VEH_20$bin)==TRUE, c(2, 5, 6)]

#Chi-sq analysis excitations:
#Bins 1, 3, 5, 7 and 9
critData <- contingency_table_comp1[, -1]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 38.727, df = 4, p-value = 7.931e-08

fisher.test(critData)
# data:  critData
# p-value = 4.798e-08
# alternative hypothesis: two.sided



#Bins 2, 4, 6, 8 and 10
critData <- contingency_table_comp2[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 48.612, df = 4, p-value = 7.034e-10

fisher.test(critData)
# data:  critData
# p-value = 5.577e-08
# alternative hypothesis: two.sided



#Post-hoc comparisons for % cue-excited units
# Is the distribution of exc/not exc neurons different from the first and the third bin (for example). 
contTableExc_VEH_20 <- PercCueExc_VEH_20[ , c(3, 4)]

#Bin 2 vs. 5
fisher.test(contTableExc_VEH_20[c(2, 5), ], alternative="less") #95% CI: 0.000000  0.6672858 Odds ratio=0.3864282 p=0.00141
 #Bin 5 vs. 8
fisher.test(contTableExc_VEH_20[c(5, 8), ], alternative="less") #95% CI:0.0000000 0.7243426 odds ratio=0.374803 p= 0.005289

#Correct p value all comparisons
p.adjust(p=c(0.00141, 0.005289), method="holm")
#  0.002820 0.005289


contTableInh_VEH_20 <- PercCueExc_VEH_20[ , c(5, 6)]

#Bin 2 vs. 5
fisher.test(contTableInh_VEH_20[c(2, 5), ], alternative="greater") #95% CI:0.8437931      Inf Odds ratio=1.617342 p=0.1239
#Bin 5 vs. 8
fisher.test(contTableInh_VEH_20[c(5, 8), ], alternative="greater") #95% CI:0.6915767 Inf odds ratio= 1.64495  p-val = 0.2006



### AP5 SIDE
contingency_table_comp1<- PercCueExc_AP5_20[is.even(PercCueExc_AP5_20$bin)==FALSE, c(2, 3, 4)]
contingency_table_comp2 <-  PercCueExc_AP5_20[is.even(PercCueExc_AP5_20$bin)==TRUE, c(2, 3, 4)]

inhCont_table_comp1 <- PercCueExc_AP5_20[is.even(PercCueExc_AP5_20$bin)==FALSE, c(2, 5, 6)]
inhCont_table_comp2 <- PercCueExc_AP5_20[is.even(PercCueExc_AP5_20$bin)==TRUE, c(2, 5, 6)]

#Chi-sq analysis excitations:
#Bins 1, 3 and 5
critData <- contingency_table_comp1[, -1]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 10.111, df = 4, p-value = 0.03859

fisher.test(critData) #p-value = 0.03694

#Bins 2, 4, 6, and 8
critData <- contingency_table_comp2[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 14.355, df = 4, p-value = 0.006244

fisher.test(x=critData) #p-value = 0.005615


#Post-hoc comparisons for % cue-excited units
# Is the distribution of exc/not exc neurons different from the first (-120 to -80) and the third (-40 to 0) bin. 
contTableExc_AP5_20 <- PercCueExc_AP5_20[ , c(3, 4)]

#Bin 2 vs. 5
fisher.test(contTableExc_AP5_20[c(2, 5), ], alternative="less") #95% CI: 0.000000 0.9821509 Odds ratio=0.499756 p= 0.0452
#Bin 5 vs. 8
fisher.test(contTableExc_AP5_20[c(5, 8), ], alternative="less") #95% CI:0.0000000 1.631982 odds ratio=0.8466628 p= 0.3907

#Correct p value all comparisons
p.adjust(p=c( 0.0452, 0.3907), method="holm")
# 0.0904 0.39070



#Chi-sq for inhibitions:
#Bins 1, 3 and 5
critData <- inhCont_table_comp1[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 1.8244, df = 4, p-value = 0.768

fisher.test(critData) #p-value = 0.7779


#Bins 2, 4 and 6
critData <- inhCont_table_comp2[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 2.059, df = 4, p-value = 0.7249
# 

fisher.test(x=critData) #p-value = 0.7397

contTableInh_AP5_20 <- PercCueExc_AP5_20[ , c(5, 6)]


#Bin 2 vs. 5
fisher.test(contTableInh_AP5_20[c(2, 5), ], alternative="greater") #95% CI: 0.8063394       Inf Odds ratio=1.663381 p= 0.1395
#Bin 5 vs. 8
fisher.test(contTableInh_AP5_20[c(5, 8), ], alternative="greater") #95% CI:  0.4104399  Inf odds ratio=0.9246439 p= 0.6516



#########

### % CUE EXCITED
### VEH VS AP5
#Bin 1: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH_20[1, ], contTableExc_AP5_20[1, ]), alternative="less") #95% CI: 0.000000 2.393145  Odds ratio=1.101717 p= 0.6692

#Bin 2: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH_20[2, ], contTableExc_AP5_20[2, ]), alternative="greater") #95% CI: 0.8820046 inf Odds ratio=1.513451 p= 0.1107

#Bin 3: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH_20[3, ], contTableExc_AP5_20[3, ]), alternative="greater") #95% CI: 1.047456 inf Odds ratio=2.054446 p= 0.03813

#Bin 4: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH_20[4, ], contTableExc_AP5_20[4, ]), alternative="greater") #95% CI: 1.393957 inf Odds ratio=2.291379  p= 0.002132

#Bin 5: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH_20[5, ], contTableExc_AP5_20[5, ]), alternative="greater") #95% CI: 0.99437 inf Odds ratio=1.950252 p=  0.05162

#Bin 6: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH_20[6, ], contTableExc_AP5_20[6, ]), alternative="greater") #95% CI: 1.770777 inf Odds ratio=3.474471 p= 0.0002758

#Bin 7: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH_20[7, ], contTableExc_AP5_20[7, ]), alternative="greater") #95% CI: 1.770777 inf Odds ratio=3.474471  p= 0.0006186

#Bin 8: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH_20[8, ], contTableExc_AP5_20[8, ]), alternative="greater") #95% CI: 2.336163 inf Odds ratio=4.409919  p=  1.446e-05

#Bin 9: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH_20[9, ], contTableExc_AP5_20[9, ]), alternative="greater") #95% CI: 1.130981 inf Odds ratio=2.469474 p= 0.02598

#Bin 10: VEH vs. AP5
fisher.test(rbind(contTableExc_VEH_20[10, ], contTableExc_AP5_20[10, ]), alternative="greater") #95% CI: 0.3435645 inf Odds ratio=1.044715 p= 0.5845


p.adjust(p=c(0.6692, 0.1107, 0.03813, 0.002132, 0.05162, 0.0002758, 0.0006186, 1.446e-05, 0.02598, 0.5845))
# 1.0000000 0.3321000 0.1906500 0.0149240 0.2064800 0.0024822 0.0049488 0.0001446 0.1558800 1.0000000

giveStars(p.adjust(p=c(0.6692, 0.1107, 0.03813, 0.002132, 0.05162, 0.0002758, 0.0006186, 1.446e-05, 0.02598, 0.5845), method="holm"))



### % CUE INHIBITED
### VEH VS AP5
#Bin 1: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH_20[1, ], contTableInh_AP5_20[1, ]), alternative="greater") 
#95% CI: 0.7490693 Inf  Odds ratio=1.495257 p= 0.1918

#Bin 2: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH_20[2, ], contTableInh_AP5_20[2, ]), alternative="less") 
#95% CI: 0.00000 1.42798 Odds ratio=0.8411593 p= 0.3296

#Bin 3: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH_20[3, ], contTableInh_AP5_20[3, ]), alternative="less") 
#95% CI: 0.0000000 0.8789947 Odds ratio=0.4084578 p= 0.02481

#Bin 4: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH_20[4, ], contTableInh_AP5_20[4, ]), alternative="less") 
#95% CI: 0.00000 1.01832 Odds ratio=0.4084578  p= 0.02481

#Bin 5: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH_20[5, ], contTableInh_AP5_20[5, ]), alternative="less") 
#95% CI: 0.000000 2.043556 Odds ratio=0.8663043 p=  0.4633

#Bin 6: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH_20[6, ], contTableInh_AP5_20[6, ]), alternative="less") 
#95% CI: 0.000000 1.071944 Odds ratio=0.5847401 p= 0.07538

#Bin 7: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH_20[7, ], contTableInh_AP5_20[7, ]), alternative="less") 
#95% CI: 0.000000 1.314831 Odds ratio=0.5666588  p=  0.1499

#Bin 8: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH_20[8, ], contTableInh_AP5_20[8, ]), alternative="less") 
#95% CI: 0.000000 1.089567 Odds ratio=0.4866311 p=  0.07472

#Bin 9: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH_20[9, ], contTableInh_AP5_20[9, ]), alternative="less") 
#95% CI: 0.000000 1.659497 Odds ratio=0.635533 p= 0.2712

#Bin 10: VEH vs. AP5
fisher.test(rbind(contTableInh_VEH_20[10, ], contTableInh_AP5_20[10, ]), alternative="less") 
#95% CI: 0.000000 2.005612 Odds ratio=0.4468194 p= 0.2667



# 
# #Correct p value of ALL comparisons:
# VEHcomps.pval <- c(0.0006973, 0.0001357, 2.979e-08, 0.6879)
# AP5comps.pval <-  c(0.08469, 0.08012, 0.003075, 0.2676)
# VEHvsAP5compsEXC.pval <- c(0.7807, 0.06973, 0.009084, 0.0001753, 3.379e-05, 0.01534)
# VEHvsAP5compsINH.pval <- c(0.3848, 0.3052, 0.1769, 0.1021, 0.1103, 0.2048)
# 
# p.adjust(p=c(VEHcomps.pval, AP5comps.pval, VEHvsAP5compsEXC.pval, VEHvsAP5compsINH.pval), method="holm")
# 
# VEHcomps.pval.adj <- giveStars(c( 1.11568e-02, 2.44260e-03, 5.95800e-07, 1.00000e+00))
# AP5comps.pval.adj <- giveStars(c( 8.81320e-01, 8.81320e-01, 4.61250e-02, 1.00000e+00))
# VEHvsAP5compsEXC.pval.adj <- giveStars(c(1.00000e+00, 8.36760e-01, 1.27176e-01, 2.98010e-03, 6.42010e-04, 1.99420e-01))
# VEHvsAP5compsINH.pval.adj <- giveStars(c(1.00000e+00, 1.00000e+00, 1.00000e+00, 9.18900e-01, 9.18900e-01, 1.00000e+00))
# 







###########################
#BEFORE ENTRY
plotFRandCP_Entry(experiment="Exp 4 VEH side", cue=c("S+", "S-"), masterDF=list(masterDF_DSEntry_VEH, masterDF_NSEntry_VEH), graphFolder=MixedGraphFolder, 
            trialBinSize=35, WdwStart=-2000, WdwEnd=0, dataProcess="Zscores", correctOnly=FALSE, 
            colindx =c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), capped=T, capValue = c(-105, 105), 
            yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS_VEH, wdwLabel="Pre S+ Entry")

plotFRandCP_Entry(experiment="Exp 4 AP5 side", cue=c("S+", "S-"), masterDF=list(masterDF_DSEntry_AP5, masterDF_NSEntry_AP5), graphFolder=MixedGraphFolder, 
                  trialBinSize=35, WdwStart=-2000, WdwEnd=0, dataProcess="Zscores", correctOnly=FALSE, 
                  colindx =c(colindx[2], "darkred"), legLabels=c("S+", "S-"), capped=T, capValue = c(-105, 105), 
                  yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS_AP5, wdwLabel="Pre S+ Entry")

plotFRandCP_Entry(experiment="Exp 4 VEH vs AP5 side", cue=c("S+", "S+"), masterDF=list(masterDF_DSEntry_VEH, masterDF_DSEntry_AP5), graphFolder=MixedGraphFolder, 
                  trialBinSize=35, WdwStart=-2000, WdwEnd=0, dataProcess="Zscores", correctOnly=FALSE, 
                  colindx = colindx, legLabels=c("S+", "S-"), capped=T, capValue = c(-105, 105), 
                   yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS_AP5, wdwLabel="Pre S+ Entry")



#AFTER ENTRY
plotFRandCP_Entry(experiment="Exp 4 VEH side", cue=c("S+", "S-"), masterDF=list(masterDF_DSEntry_VEH, masterDF_NSEntry_VEH), graphFolder=MixedGraphFolder, 
                  trialBinSize=35, WdwStart=0, WdwEnd=1000, dataProcess="Zscores", correctOnly=FALSE, 
                  colindx =c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), capped=T, capValue = c(-105, 105), 
                  yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS_VEH, wdwLabel="Post S+ Entry")

plotFRandCP_Entry(experiment="Exp 4 AP5 side", cue=c("S+", "S-"), masterDF=list(masterDF_DSEntry_AP5, masterDF_NSEntry_AP5), graphFolder=MixedGraphFolder, 
                  trialBinSize=35, WdwStart=0, WdwEnd=1000, dataProcess="Zscores", correctOnly=FALSE, 
                  colindx =c(colindx[2], "darkred"), legLabels=c("S+", "S-"), capped=T, capValue = c(-105, 105), 
                  yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS_AP5, wdwLabel="Post S+ Entry")

plotFRandCP_Entry(experiment="Exp 4 VEH vs AP5 side", cue=c("S+", "S+"), masterDF=list(masterDF_DSEntry_VEH, masterDF_DSEntry_AP5), graphFolder=MixedGraphFolder, 
                  trialBinSize=35, WdwStart=0, WdwEnd=1000, dataProcess="Zscores", correctOnly=FALSE, 
                  colindx = colindx, legLabels=c("S+", "S-"), capped=T, capValue = c(-105, 105), 
                  yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS_AP5, wdwLabel="Post S+ Entry")


### PLOT FR AROUND ENTRIES (PTSH) as a function of distance to change point
plotFRandCPhistogram_Entry(experiment="Exp 4 Entry Hist over 2s lat", masterDF=list(masterDF_DSEntry_VEH, masterDF_DSEntry_AP5), 
                           graphFolder=MixedGraphFolder, cueExcOnly=F, trialBinSize=35, dataProcess="Zscores", 
                           correctOnly=FALSE, color=colindx, 
                           capped=T, capValue = c(-105, 105), yAxMinZ = -1, yAxMaxZ = 2, yAxMaxRaw = 3, 
                           WdwStart=-2000, WdwEnd=5000, imgFormat="pdf", neudata=allNeuronsEntryDS)

plotFRandCPhistogram_Entry(experiment="Exp 4 Entry Hist over 5s lat", masterDF=list(masterDF_DSEntry_VEH, masterDF_DSEntry_AP5), 
                           graphFolder=MixedGraphFolder, cueExcOnly=F, trialBinSize=35, dataProcess="Zscores", 
                           correctOnly=FALSE, color=colindx,
                           capped=T, capValue = c(-105, 105), yAxMinZ = -1, yAxMaxZ = 2, yAxMaxRaw = 3, 
                           WdwStart=-5000, WdwEnd=5000, imgFormat="pdf", neudata=allNeuronsEntryDS)

plotFRandCPhistogram_Entry(experiment="Exp 4 Entry Hist over 8s lat", masterDF=list(masterDF_DSEntry_VEH, masterDF_DSEntry_AP5), 
                           graphFolder=MixedGraphFolder, cueExcOnly=F, trialBinSize=35, dataProcess="Zscores", 
                           correctOnly=FALSE, color=colindx, 
                           capped=T, capValue = c(-105, 105), yAxMinZ = -1, yAxMaxZ = 2, yAxMaxRaw = 3, 
                           WdwStart=-8000, WdwEnd=5000, imgFormat="pdf", neudata=allNeuronsEntryDS)



### PLOT FR AROUND ENTRIES (CP) as a function of distance to change point

plotFRandCP_Entry(cue=c("S+", "S+"), wdwLabel="Post S+", experiment="Exp 4 FR 2s preentry when latency is more than 5s", 
                  masterDF=list(masterDF_DSEntry_VEH, masterDF_DSEntry_AP5), 
                  graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=-2000, WdwEnd=0, dataProcess="Zscores", minLat=5000,
                  correctOnly=F, colindx=colindx, legLabels=c("VEH", "AP5"), yAxMinZ=-1, yAxMaxZ=2, yAxMaxRaw=7, 
                  capped=T, capValue=c(-105, 105), cueExcOnly=F, neudata=allNeuronsDS)

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


#PRE cue VEH vs. AP5 (BASELINE)
BL_35bin <- plotFRBoxPlotandCP(cue=c("S+", "S+"), experiment="Exp 4 BASELINE firing from CP VEH vs. AP5", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                                graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=-2000, WdwEnd=0, dataProcess="raw", 
                                correctOnly=FALSE, points=F, lines=F, color=colindx, legLabels=c("S+", "S+"), 
                                yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, capValue=c(-105, 105), cueExcOnly=F, neudata=allNeuronsDS)

BL_All_trials <- plotFRBoxPlotandCP(cue=c("S+", "S+"), experiment="Exp 4 BASELINE VEH vs. AP5", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                               graphFolder=MixedGraphFolder, trialBinSize=420, WdwStart=-2000, WdwEnd=0, dataProcess="raw", 
                               correctOnly=FALSE, points=T, lines=F, color=colindx, legLabels=c("S+", "S+"), 
                               yAxMinZ=-1, yAxMaxZ=15, yAxMaxRaw=7, capped=T, capValue=c(-210, 210), cueExcOnly=F, neudata=allNeuronsDS)


#Post cue VEH
VEH_35bin <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 4 VEH", masterDF=list(masterDF_DS_VEH, masterDF_NS_VEH), 
                                   graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="Zscores", 
                                   correctOnly=FALSE, points=F, lines=F, color=c(colindx[1], colindxB[1]), legLabels=c("S+", "S-"), 
                                   yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, capValue=c(-105, 105), cueExcOnly=F, neudata=allNeuronsDS)


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
compareBins(data=Diff_DSNS_VEH_35bins, cueExcOnly=F, color=colindx[1], ymin=-2, ymax=8, 
            graphFolder=MixedGraphFolder, experiment="Exp 4 VEH", points=F, 
            comparisons=list(c(1, 3), c(3, 5)), 
            cue="S+")
# [[1]]
# idx     bin         n                     comp    W        p.val        p.adj sig
# 1 vs. 3 19 vs. 63 -105 to -70 vs. -35 to 0  192 4.059723e-06 1.623889e-05 ***
# 2 vs. 4 45 vs. 87   -70 to -35 vs. 0 to 35 1221 2.052271e-04 6.156813e-04 ***
# 3 vs. 5 63 vs. 68    -35 to 0 vs. 35 to 70 1688 1.835001e-02 1.835001e-02   *
# 4 vs. 6 87 vs. 23    0 to 35 vs. 70 to 105  678 8.971326e-03 1.794265e-02   *
#         
#         [[2]]
# bins         n               comparison    W        p.val        p.adj sig
# 1 vs. 3 19 vs. 63 -105 to -70 vs. -35 to 0  192 4.059723e-06 8.119445e-06 ***
# 3 vs. 5 63 vs. 68    -35 to 0 vs. 35 to 70 1688 1.835001e-02 1.835001e-02   *
#         


# [[1]] VEH ALL UNITS
# bin           n                     comp    W        p.val        p.adj sig
# 1 vs. 3  75 vs. 129 -105 to -70 vs. -35 to 0 2732 9.843124e-07 2.952937e-06 ***
# 2 vs. 4 137 vs. 131   -70 to -35 vs. 0 to 35 5117 1.359566e-08 5.438264e-08 ***
# 3 vs. 5  129 vs. 92    -35 to 0 vs. 35 to 70 3950 1.341498e-04 2.682996e-04 ***
# 4 vs. 6  131 vs. 36    0 to 35 vs. 70 to 105 1707 5.730776e-02 5.730776e-02    
# 
# [[2]]
# bins          n               comparison    W        p.val        p.adj sig
# 1 vs. 3 75 vs. 129 -105 to -70 vs. -35 to 0 2732 9.843124e-07 1.968625e-06 ***
# 3 vs. 5 129 vs. 92    -35 to 0 vs. 35 to 70 3950 1.341498e-04 1.341498e-04 ***


#This one is like the previous one but using firing rate after S+ minus after S- as a dependent variable (y axis)
DiffFR_ByBin(data=Diff_DSNS_VEH_35bins, cueExcOnly=TRUE, color=colindx[1], ymin=-2, ymax=18,
             graphFolder=MixedGraphFolder, experiment="Exp 4 VEH", points=TRUE, 
             comparisons=list(c(1, 3), c(3, 5)))
# [[1]] cue excited only (S+ fr minus S- firing rate per unit)
# bin                     comp          W        p.val       p.adj sig
# 1 vs. 3 -105 to -70 vs. -35 to 0  104 0.0019505572 0.005851672  **
# 2 vs. 4   -70 to -35 vs. 0 to 35  740 0.0024963796 0.005851672  **
# 3 vs. 5    -35 to 0 vs. 35 to 70 1057 0.0006830055 0.002732022  **
# 4 vs. 6    0 to 35 vs. 70 to 105  460 0.1601646509 0.160164651    

# bins               comparison    W        p.val       p.adj sig
# 1 vs. 3 -105 to -70 vs. -35 to 0  104 0.0019505572 0.001950557  **
# 3 vs. 5    -35 to 0 vs. 35 to 70 1057 0.0006830055 0.001366011  **

PostCueFR_DSvsNS_fromCP_VEH <- compareDSvsNSfromCP(masterDF=list(masterDF_DS_VEH, masterDF_NS_VEH), 
                                                   trialBinSize=35, event="cue", correctOnly=F, paired=T, 
                                                   WdwStart=100, WdwEnd=400, capped=T, capValue=c(-105, 105), 
                                                   dataProcess="Zscores")

PostCueFR_DSvsNS_fromCP_VEH$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP_VEH$p, method="holm") 
PostCueFR_DSvsNS_fromCP_VEH$sig <- giveStars(PostCueFR_DSvsNS_fromCP_VEH$p.adj)

# bin           V            p          n        p.adj          sig
# -105 to -70   1273 3.360339e-01       73      3.360339e-01    
# -70 to -35    6257 5.884856e-05       134     1.176971e-04 ***
# -35 to 0      6238 1.064615e-10       122     4.258458e-10 ***
# 0 to 35       6955 5.280913e-14       125     2.882187e-13 ***
# 35 to 70      3823 4.803645e-14       89      2.882187e-13 ***
# 70 to 105     224 9.059906e-06        21      2.717972e-05 ***

#Post cue AP5
AP5_35bin <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 4 AP5", masterDF=list(masterDF_DS_AP5, masterDF_NS_AP5), 
                                graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="Zscores", 
                                correctOnly=FALSE, points=F, lines=F, color=c(colindx[2], colindxB[2]), legLabels=c("S+", "S-"), 
                                yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, capValue=c(-105, 105), cueExcOnly=F, neudata=allNeuronsDS)


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

compareBins(data=Diff_DSNS_AP5_35bins, cueExcOnly=T, color=colindx[2], ymin=-2, ymax=8, 
            graphFolder=MixedGraphFolder, experiment="Exp 4 AP5", points=F, 
            comparisons=list(c(1, 3), c(3, 5)), 
            cue="S+")

# [[1]] ALL UNITS
# bin           n                     comp      W       p.val     p.adj sig
# 1 vs. 3   66 vs. 96 -105 to -70 vs. -35 to 0 2058 0.003546225 0.0141849   *
# 2 vs. 4 103 vs. 102   -70 to -35 vs. 0 to 35 3907 0.034510170 0.1035305    
# 3 vs. 5   96 vs. 85    -35 to 0 vs. 35 to 70 3745 0.763052411 0.7630524    
# 4 vs. 6  102 vs. 75    0 to 35 vs. 70 to 105 3141 0.225433970 0.4508679    
# 
# [[2]]
# bins         n               comparison    W       p.val       p.adj sig
# 1 vs. 3 66 vs. 96 -105 to -70 vs. -35 to 0 2058 0.003546225 0.007092449  **
# 3 vs. 5 96 vs. 85    -35 to 0 vs. 35 to 70 3745 0.763052411 0.763052411    


# CUE EXCITED ONLY
# [[1]]
# idx     bin         n                     comp   W       p.val      p.adj sig
# 1 vs. 3 14 vs. 31 -105 to -70 vs. -35 to 0 122 0.009603846 0.03841538   *
# 2 vs. 4 24 vs. 43   -70 to -35 vs. 0 to 35 346 0.012950995 0.03885299   *
# 3 vs. 5 31 vs. 37    -35 to 0 vs. 35 to 70 534 0.316596296 0.42452383    
# 4 vs. 6 43 vs. 36    0 to 35 vs. 70 to 105 692 0.212261916 0.42452383    
# 
# [[2]]
# bins         n               comparison   W       p.val      p.adj sig
# 1 vs. 3 14 vs. 31 -105 to -70 vs. -35 to 0 122 0.009603846 0.01920769   *
# 3 vs. 5 31 vs. 37    -35 to 0 vs. 35 to 70 534 0.316596296 0.31659630    



DiffFR_ByBin(data=Diff_DSNS_AP5_35bins, cueExcOnly=TRUE, color=colindx[2], ymin=-2, ymax=18,
             graphFolder=MixedGraphFolder, experiment="Exp 4 AP5", points=TRUE, 
             comparisons=list(c(1, 3), c(3, 5)))

# [[1]] ALL UNITS
# idx     bin                        comp   W      p.val     p.adj sig
# W     1 1 vs. 3 -140 to -105 vs. -70 to -35   4 0.05263158 0.2175159    
# W3    2 2 vs. 4    -105 to -70 vs. -35 to 0  56 0.03348018 0.2008811    
# W1    3 3 vs. 5      -70 to -35 vs. 0 to 35 190 0.04350318 0.2175159    
# W11   4 4 vs. 6       -35 to 0 vs. 35 to 70 230 0.12095389 0.3628617    
# W2    5 5 vs. 7       0 to 35 vs. 70 to 105 264 0.14973434 0.3628617    
# W21   6 6 vs. 8     35 to 70 vs. 105 to 140  67 0.46866448 0.4686645    
# 
# [[2]]
# bins               comparison      W      p.val      p.adj sig
# 2 vs. 4 -105 to -70 vs. -35 to 0  56 0.03348018 0.06696035    
# 4 vs. 6    -35 to 0 vs. 35 to 70 230 0.12095389 0.12095389  


# [[1]] CUE-EXCITED UNITS ONLY
# idx           comp               W       p.val      p.adj     sig
# 1 vs. 3 -105 to -70 vs. -35 to 0 103 0.005200363 0.02080145   *
# 2 vs. 4   -70 to -35 vs. 0 to 35 413 0.090738083 0.21107145    
# 3 vs. 5    -35 to 0 vs. 35 to 70 513 0.430114794 0.43011479    
# 4 vs. 6    0 to 35 vs. 70 to 105 366 0.070357149 0.21107145    
# 
# [[2]]
# bins           comparison             W       p.val      p.adj        sig
# 1 vs. 3 -105 to -70 vs. -35 to 0      103 0.005200363 0.01040073      *
# 3 vs. 5    -35 to 0 vs. 35 to 70      513 0.430114794 0.43011479    



PostCueFR_DSvsNS_fromCP_AP5 <- compareDSvsNSfromCP(masterDF=list(masterDF_DS_AP5, masterDF_NS_AP5), cueExcOnly = F,
                                                   trialBinSize=35, event="cue", correctOnly=F, paired=T, 
                                                   WdwStart=100, WdwEnd=400, capped=T, capValue=c(-105, 105), 
                                                   dataProcess="Zscores")

PostCueFR_DSvsNS_fromCP_AP5$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP_AP5$p, method="holm") 
PostCueFR_DSvsNS_fromCP_AP5$sig <- giveStars(PostCueFR_DSvsNS_fromCP_AP5$p.adj)

# bin           V            p          n        p.adj sig
# -105 to -70   1331 3.860282e-04       59 3.860282e-04 ***
# -70 to -35    3600 5.396597e-06       97 1.618979e-05 ***
# -35 to 0      3302 1.134632e-08       88 6.807792e-08 ***
# 0 to 35       3658 1.584687e-07       95 7.923434e-07 ***
# 35 to 70      2447 8.004642e-07       77 3.201857e-06 ***
# 70 to 105     770 1.281043e-05        42 2.562085e-05 ***


padjall <- p.adjust(c(PostCueFR_DSvsNS_fromCP_VEH$p, PostCueFR_DSvsNS_fromCP_AP5$p), method="holm") 

PostCueFR_DSvsNS_fromCP_VEH$p.adj <- padjall[1:nrow(PostCueFR_DSvsNS_fromCP_VEH)]
PostCueFR_DSvsNS_fromCP_AP5$p.adj <- padjall[(nrow(PostCueFR_DSvsNS_fromCP_VEH)+1):length(padjall)]



#Post cue VEH vs AP5
VEHAP5_35bin <- plotFRBoxPlotandCP(cue=c("S+", "S+"), experiment="Exp 4 VEH vs AP5 Spike", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                                graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="Zscores", 
                                correctOnly=FALSE, points=F, lines=F, color=colindx, legLabels=c("S+", "S+"), 
                                yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, capValue=c(-105, 105), cueExcOnly=F, neudata=allNeuronsDS_VEH)

VEHAP5_35bin_Tail <- plotFRBoxPlotandCP(cue=c("S+", "S+"), experiment="Exp 4 VEH vs AP5 Tail", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                                   graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=750, WdwEnd=2000, dataProcess="Zscores", 
                                   correctOnly=FALSE, points=F, lines=F, color=colindx, legLabels=c("S+", "S+"), 
                                   yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, capValue=c(-105, 105), cueExcOnly=F, neudata=allNeuronsDS_VEH)


VEHAP5_35bin_CueExc <- plotFRBoxPlotandCP(cue=c("S+", "S+"), experiment="Exp 4 VEH vs AP5 Spike Cue Exc", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                                   graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="Zscores", 
                                   correctOnly=FALSE, points=F, lines=F, color=colindx, legLabels=c("S+", "S+"), 
                                   yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, capValue=c(-105, 105), cueExcOnly=T, neudata=allNeuronsDS_VEH)

VEHAP5_35bin_CueExc_Tail <- plotFRBoxPlotandCP(cue=c("S+", "S+"), experiment="Exp 4 VEH vs AP5 Tail Cue Exc", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                                          graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=750, WdwEnd=2000, dataProcess="Zscores", 
                                          correctOnly=FALSE, points=F, lines=F, color=colindx, legLabels=c("S+", "S+"), 
                                          yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, capValue=c(-105, 105), cueExcOnly=T, neudata=allNeuronsDS_VEH)


VEHAP5_35bin_DSEntry <- plotFRBoxPlotandCP(cue=c("S+", "S+"), experiment="Exp 4 VEH vs AP5 DS Entry", masterDF=list(masterDF_DSEntry_VEH, masterDF_DSEntry_AP5), 
                                               graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=0, WdwEnd=1000, dataProcess="Zscores", 
                                               correctOnly=FALSE, points=F, lines=F, color=colindx, legLabels=c("S+", "S+"), 
                                               yAxMinZ=-1, yAxMaxZ=3, yAxMaxRaw=7, capped=T, capValue=c(-105, 105), cueExcOnly=T, neudata=allNeuronsDS_VEH)

VEHAP5_35bin_NSEntry <- plotFRBoxPlotandCP(cue=c("S+", "S+"), experiment="Exp 4 VEH vs AP5 NS Entry", masterDF=list(masterDF_NSEntry_VEH, masterDF_NSEntry_AP5), 
                                           graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=0, WdwEnd=1000, dataProcess="Zscores", 
                                           correctOnly=FALSE, points=F, lines=F, color=colindx, legLabels=c("S+", "S+"), 
                                           yAxMinZ=-1, yAxMaxZ=3, yAxMaxRaw=7, capped=T, capValue=c(-105, 105), cueExcOnly=T, neudata=allNeuronsDS_VEH)

VEHAP5_35bin_ITIEntry <- plotFRBoxPlotandCP(cue=c("S+", "S+"), experiment="Exp 4 VEH vs AP5 ITI Entry", masterDF=list(masterDF_ITIEntry, masterDF_ITIEntry_AP5), 
                                           graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=0, WdwEnd=1000, dataProcess="Zscores", 
                                           correctOnly=FALSE, points=F, lines=F, color=colindx, legLabels=c("S+", "S+"),
                                           yAxMinZ=-1, yAxMaxZ=3, yAxMaxRaw=7, capped=T, capValue=c(-105, 105), cueExcOnly=T, neudata=allNeuronsDS_VEH)


# All units VEH vs AP5 side by bin of trials with respect to CP: 100-400ms after S+
PostCueFR_DS_fromCP_VEHvsAP5_Spike <- compareVEHvsAP5fromCP(masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), trialBinSize=35, 
                                                      event="cue", correctOnly=F, cueExcOnly = F, paired=F, WdwStart=100, WdwEnd=400, 
                                                      capped=T, capValue=c(-105, 105), dataProcess="Zscores")
PostCueFR_DS_fromCP_VEHvsAP5_Spike$p.adj <- p.adjust(PostCueFR_DS_fromCP_VEHvsAP5_Spike$p, method="holm") 
PostCueFR_DS_fromCP_VEHvsAP5_Spike$sig <- giveStars(PostCueFR_DS_fromCP_VEHvsAP5_Spike$p.adj)
# bin                   n               W               p        p.adj sig
# -105 to -70           75 vs.  66      1888 0.0279973334 0.0839920003    
# -70 to -35            137 vs. 103     6370 0.3627862732 0.5941488691    
# -35 to 0              129 vs. 96      5691 0.2970744345 0.5941488691    
# 0 to 35               131 vs. 102     7249 0.0036266233 0.0181331167   *
# 35 to 70              92 vs.  85      4768 0.0001337584 0.0008025505 ***
# 70 to 105             36 vs.  75      1510 0.0091221475 0.0364885901   *



# All units VEH vs AP5 side by bin of trials with respect to CP: 750-2000ms after S+ (excluding trials in which the latency was less than 2s):
PostCueFR_DS_fromCP_VEHvsAP5_Tail <- compareVEHvsAP5fromCP(masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), trialBinSize=35, 
                                                      event="cue", correctOnly=F, cueExcOnly = F, paired=F, WdwStart=750, WdwEnd=2000, 
                                                      capped=T, capValue=c(-105, 105), dataProcess="Zscores")
PostCueFR_DS_fromCP_VEHvsAP5_Tail$p.adj <- p.adjust(PostCueFR_DS_fromCP_VEHvsAP5_Tail$p, method="holm") 
PostCueFR_DS_fromCP_VEHvsAP5_Tail$sig <- giveStars(PostCueFR_DS_fromCP_VEHvsAP5_Tail$p.adj)
# bin               n           W            p        p.adj     sig
# W  -105 to -70   75 vs. 66    2183 7.393506e-01 9.541374e-01    
# W1  -70 to -35  137 vs. 103   6577 4.770687e-01 9.541374e-01    
# W2    -35 to 0  129 vs. 96    6186 4.870135e-02 1.461041e-01    
# W3     0 to 35  131 vs. 102   8559 2.265593e-08 1.359356e-07 ***
# W4    35 to 70   74 vs. 83    3830 6.086496e-05 2.434598e-04 ***
# W5   70 to 105   36 vs. 75    1821 2.912077e-06 1.456039e-05 ***



# All units VEH vs AP5 side by bin of trials with respect to CP: BASELINE PERIOD (2s Window right before the cue):
PostCueFR_DS_fromCP_VEHvsAP5_BL <- compareVEHvsAP5fromCP(masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), trialBinSize=420, 
                                                           event="cue", correctOnly=F, cueExcOnly = F, paired=F, WdwStart=-2000, WdwEnd=0, 
                                                           capped=T, capValue=c(-210, 210), dataProcess="Zscores")
PostCueFR_DS_fromCP_VEHvsAP5_BL$p.adj <- p.adjust(PostCueFR_DS_fromCP_VEHvsAP5_BL$p, method="holm") 
PostCueFR_DS_fromCP_VEHvsAP5_BL$sig <- giveStars(PostCueFR_DS_fromCP_VEHvsAP5_BL$p.adj)

# In 35 trials from CP
#           bin           n    W          p  p.adj sig
#W  -105 to -70   75 vs. 66 2480 0.74075130 1.0000    
#W1  -70 to -35 137 vs. 103 5882 0.09361666 0.5617    
#W2    -35 to 0  129 vs. 96 4933 0.11756841 0.5617    
#W3     0 to 35 131 vs. 102 5396 0.10554280 0.5617    
#W4    35 to 70   92 vs. 85 3681 0.40078444 1.0000    
#W5   70 to 105   36 vs. 75 1121 0.36346789 1.0000  


# All trials together:
# bin           n            W         p    p.adj       sig
# W -210 to 210 337 vs. 305 45592 0.3495017 0.3495017    


# CUE-EXCITED units only VEH vs AP5 side by bin of trials with respect to CP: 100-400ms after S+
PostCueFR_DS_fromCP_VEHvsAP5_Spike_CueExc <- compareVEHvsAP5fromCP(masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), trialBinSize=35, 
                                                            event="cue", correctOnly=F, cueExcOnly = T, paired=F, WdwStart=100, WdwEnd=400, 
                                                            capped=T, capValue=c(-105, 105), dataProcess="Zscores")
PostCueFR_DS_fromCP_VEHvsAP5_Spike_CueExc$p.adj <- p.adjust(PostCueFR_DS_fromCP_VEHvsAP5_Spike_CueExc$p, method="holm") 
PostCueFR_DS_fromCP_VEHvsAP5_Spike_CueExc$sig <- giveStars(PostCueFR_DS_fromCP_VEHvsAP5_Spike_CueExc$p.adj)
### cue excited only
#       bin         n        W          p     p.adj sig
#W  -105 to -70 19 vs. 14   94 0.08146491 0.4073245    
#W1  -70 to -35 45 vs. 24  644 0.09685300 0.4073245    
#W2    -35 to 0 63 vs. 31  997 0.43610772 0.4835020    
#W3     0 to 35 87 vs. 43 2071 0.16116733 0.4835020    
#W4    35 to 70 68 vs. 37 1399 0.17298094 0.4835020    
#W5   70 to 105 23 vs. 36  530 0.03623141 0.2173885    



# CUE-EXCITED units only VEH vs AP5 side by bin of trials with respect to CP: 750-2000ms after S+
PostCueFR_DS_fromCP_VEHvsAP5_Tail_CueExc <- compareVEHvsAP5fromCP(masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), trialBinSize=35, 
                                                                   event="cue", correctOnly=F, cueExcOnly = T, paired=F, WdwStart=750, WdwEnd=2000, 
                                                                   capped=T, capValue=c(-105, 105), dataProcess="Zscores")
PostCueFR_DS_fromCP_VEHvsAP5_Tail_CueExc$p.adj <- p.adjust(PostCueFR_DS_fromCP_VEHvsAP5_Tail_CueExc$p, method="holm") 
PostCueFR_DS_fromCP_VEHvsAP5_Tail_CueExc$sig <- giveStars(PostCueFR_DS_fromCP_VEHvsAP5_Tail_CueExc$p.adj)
# bin         n    W            p        p.adj sig
# W  -105 to -70 19 vs. 14  115 2.644946e-01 0.5289892810    
# W1  -70 to -35 45 vs. 24  624 1.475508e-01 0.4426525163    
# W2    -35 to 0 63 vs. 31 1054 2.678747e-01 0.5289892810    
# W3     0 to 35 87 vs. 43 2406 4.056008e-03 0.0162240325   *
# W4    35 to 70 56 vs. 37 1431 9.789170e-04 0.0048945849  **
# W5   70 to 105 23 vs. 36  659 4.061035e-05 0.0002436621 ***





#Pre entry
PreEntryFR_DSvsNS_fromCP_VEH <- compareDSvsNSfromCP(masterDF=list(masterDF_DSEntry_VEH, masterDF_NSEntry_VEH), 
                                                    trialBinSize=35, event="entry", correctOnly=F, WdwStart=-2000, WdwEnd=0, 
                                                    capped=T, capValue=c(-105, 105), dataProcess="Zscores")

PreEntryFR_DSvsNS_fromCP_VEH$p.adj <- p.adjust(PreEntryFR_DSvsNS_fromCP_VEH$p, method="holm") 
PreEntryFR_DSvsNS_fromCP_VEH$sig <- giveStars(PreEntryFR_DSvsNS_fromCP_VEH$p.adj)
save(PreEntryFR_DSvsNS_fromCP_VEH, file=paste(dataForRdir, "PreEntryFR_DSvsNS_fromCP_VEH.rdat"))
#            bin    V            p  n        p.adj sig
# V  -105 to -70 1130 5.236606e-01 67 5.236606e-01    
# V1  -70 to -35 2041 6.240843e-06 71 2.496337e-05 **
# V2    -35 to 0 2040 2.337739e-05 72 7.013217e-05 ***
# V3     0 to 35 3345 1.097489e-10 86 6.584934e-10 ***
# V4    35 to 70  252 4.768372e-07 22 2.384186e-06 ***
# V5   70 to 105  101 4.272461e-04 14 8.544922e-04 ***
        
PreEntryFR_DSvsNS_fromCP_AP5 <- compareDSvsNSfromCP(masterDF=list(masterDF_DSEntry_AP5, masterDF_NSEntry_AP5), 
                                                    trialBinSize=35, event="entry", correctOnly=F, WdwStart=-2000, 
                                                    WdwEnd=0, capped=T, capValue=c(-105, 105), dataProcess="Zscores")

PreEntryFR_DSvsNS_fromCP_AP5$p.adj <- p.adjust(PreEntryFR_DSvsNS_fromCP_AP5$p, method="holm") 
PreEntryFR_DSvsNS_fromCP_AP5$sig <- giveStars(PreEntryFR_DSvsNS_fromCP_AP5$p.adj)
save(PreEntryFR_DSvsNS_fromCP_AP5, file=paste(dataForRdir, "PreEntryFR_DSvsNS_fromCP_AP5.rdat"))
# bin    V            p  n        p.adj sig
# -105 to -70     590 2.483864e-01 51 0.2483863750    
# -70 to -35      1616 5.938311e-05 64 0.0003562987 **
# -35 to 0        1750 1.039769e-01 77 0.2079538365    
# 0 to 35         950 3.621088e-03 51 0.0144843505   *
# 35 to 70        373 6.599286e-03 31 0.0197978583   *
# 70 to 105       743 7.218808e-05 42 0.0003609404 ***


#Post entry
PostEntryFR_DSvsNS_fromCP_VEH <- compareDSvsNSfromCP(masterDF=list(masterDF_DSEntry_VEH, masterDF_NSEntry_VEH), trialBinSize=30, event="entry", correctOnly=F, WdwStart=0, WdwEnd=1000, capped=T, capValue=c(-90, 90), dataProcess="Zscores")
PostEntryFR_DSvsNS_fromCP_VEH$p.adj <- p.adjust(PostEntryFR_DSvsNS_fromCP_VEH$p, method="holm") 
PostEntryFR_DSvsNS_fromCP_VEH$sig <- giveStars(PostEntryFR_DSvsNS_fromCP_VEH$p.adj)
save(PostEntryFR_DSvsNS_fromCP_VEH, file=paste(dataForRdir, "PostEntryFR_DSvsNS_fromCP_VEH.rdat"))

PostEntryFR_DSvsNS_fromCP_AP5 <- compareDSvsNSfromCP(masterDF=list(masterDF_DSEntry_AP5, masterDF_NSEntry_AP5), trialBinSize=30, event="entry", correctOnly=F, WdwStart=0, WdwEnd=1000, capped=T, capValue=c(-90, 90), dataProcess="Zscores")
PostEntryFR_DSvsNS_fromCP_AP5$p.adj <- p.adjust(PostEntryFR_DSvsNS_fromCP_AP5$p, method="holm") 
PostEntryFR_DSvsNS_fromCP_AP5$sig <- giveStars(PostEntryFR_DSvsNS_fromCP_AP5$p.adj)
save(PostEntryFR_DSvsNS_fromCP_AP5, file=paste(dataForRdir, "PostEntryFR_DSvsNS_fromCP_AP5.rdat"))




### PLOT FR IN THE POSTCUE WINDOW AND TAIL PERIOD BY RAT BY TRIAL BY NEURON
PLOT_perRatTrial_FRandBEH(masterDFsumm=masterDFsumm, wdw="PostCue", threshold=2, MixedGraphFolder=MixedGraphFolder)
PLOT_perRatTrial_FRandBEH(masterDFsumm=masterDFsumm, wdw="PostSpike", threshold=2, MixedGraphFolder=MixedGraphFolder)


### PLOT FR AROUND THE CUE (PTSH) as a function of distance to change point
plotFRandCPhistogram(experiment="Exp 4 Cue Exc", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), graphFolder=MixedGraphFolder, 
                     trialBinSize=35, dataProcess="Zscores", correctOnly=F, color=colindx, cueExcOnly = T, events=c("S+","S+"),
                     capped=T, capValue = c(-105, 105), yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 3, 
                     psthmin=0.5, psthmax=2, imgFormat="pdf", neudata=allNeuronsDS_VEH)


#Same thing but by session from CP instead of trial
plotPSTHfromSessCP(experiment="Exp 4", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                          graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                          correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-3, 3), 
                          yAxMinZ = -1, yAxMaxZ = 4, yAxMaxRaw = 10, psthmin=-0.5, psthmax=2, 
                          imgFormat="pdf", neudata=allNeuronsDS_VEH)

plotPSTHfromSessCP(experiment="Exp 4 Tone vs Light VEH", masterDF=list(masterDF_DS_VEH_Tone, masterDF_DS_VEH_Light), 
                   graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                   correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-2, 2), 
                   yAxMinZ = -1, yAxMaxZ = 4, yAxMaxRaw = 10, psthmin=-0.5, psthmax=2, 
                   imgFormat="pdf", neudata=allNeuronsDS_VEH)

plotPSTHfromSessCP(experiment="Exp 4 Tone vs Light AP5", masterDF=list(masterDF_DS_AP5_Tone, masterDF_DS_AP5_Light), 
                   graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                   correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-2, 2), 
                   yAxMinZ = -1, yAxMaxZ = 4, yAxMaxRaw = 10, psthmin=-0.5, psthmax=2, 
                   imgFormat="pdf", neudata=allNeuronsDS_AP5)


### PLOT FR AROUND ENTRIES (PTSH) as a function of distance to change point
plotFRandCPhistogram_Entry(experiment="Exp 4", masterDF=list(masterDF_DSEntry_VEH, masterDF_DSEntry_AP5), graphFolder=MixedGraphFolder, 
                     trialBinSize=30, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
                     capped=T, capValue = c(-90, 90), yAxMinZ = -1, yAxMaxZ = 2, yAxMaxRaw = 3, 
                     WdwStart=-2000, WdwEnd=2000, imgFormat="pdf", neudata=allNeuronsEntryDS)



### PLOT BOXPLOTS AROUND CUE BY SESS FROM CP AND DRUG GROUP
plotBoxplotfromSessCP(experiment="Exp 4", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                comp=c("VEH", "AP5"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                correctOnly=FALSE, cueExcOnly=F, color=colindx, sessFromCP = c(-3, 3), 
                yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)


plotBoxplotfromSessCP(experiment="Exp 4 VEH", masterDF=list(masterDF_DS_VEH_Tone, masterDF_DS_VEH_Light), 
                      comp=c("Tone S+", "Light S+"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-2, 2), 
                      yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)

plotBoxplotfromSessCP(experiment="Exp 4 AP5", masterDF=list(masterDF_DS_AP5_Tone, masterDF_DS_AP5_Light), 
                      comp=c("Tone S+", "Light S+"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-2, 2), 
                      yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)

### PLOT BOXPLOTS AROUND CUE BY SESS FROM CP AND KIND OF CUE AND DRUG GROUP
plotBoxplotfromSessCP(experiment="Exp 4", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                      comp=c("VEH", "AP5"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=T, color=colindx, sessFromCP = c(-3, 3), 
                      yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, WdwStart=0, WdwEnd=400, 
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS_VEH, morethanIQR=T)


### THIS ONE IS THE SAME BUT GROUPING THE SESSIONS BEFORE AND AFTER THE CP
#All units: VEH vs. AP5 before and after CP
plotBoxplotPrePostCP(experiment="Exp 4", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                comp=c("VEH", "AP5"), factors=list(drug=c("VEH", "AP5")),  ANOVA=T,
                graphFolder=MixedGraphFolder, dataProcess="Zscores", removeCPsess = T,
                correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-5, 5), 
                yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS_VEH, morethanIQR=T)


### ALL UNITS
# $`ANOVA`
# Effect        DFn     DFd         F            p p<.05        ges
# sess          1       505 73.112748 1.466446e-16     * 0.12646797
# drug          1       505  6.271276 1.258527e-02     * 0.01226604
# sess:drug     1       505 10.173639 1.513124e-03     * 0.01974798
# 
# $`Levene's Test for Homogeneity of Variance`
# DFn DFd      SSn      SSd        F           p p<.05
# 1   3 505 264.2481 1546.145 28.76946 3.46182e-17     *
#         
# wilcTest
# Comparison                            W        p.val        p.adj sig
# W     Pre CP VEH vs. Pre CP AP5       8730 6.946541e-02 6.946541e-02    
# W1    Post CP VEH vs. Post CP AP5     7964 5.078548e-04 1.523564e-03  **
# W2    Pre CP VEH vs. Post CP VEH      4362 3.111524e-12 1.244610e-11 ***
# W11   Pre CP AP5 vs. Post CP AP5      5903 1.038793e-02 2.077585e-02   *
# 


#All units: VEH vs AP5 before and after CP but only the tail of the excitation (750-2000ms after S+)
plotBoxplotPrePostCP(experiment="Exp 4 TAIL", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                     comp=c("VEH", "AP5"), factors=list(drug=c("VEH", "AP5")),  ANOVA=T,
                     graphFolder=MixedGraphFolder, dataProcess="Zscores", removeCPsess = T,
                     correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-5, 5), 
                     yAxMinZ = -2, yAxMaxZ = 6, yAxMaxRaw = 10, WdwStart=750, WdwEnd=2000, 
                     removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS_VEH, morethanIQR=T)
# [[1]]
# [[1]]$`ANOVA`
# Effect DFn DFd         F            p p<.05        ges
# 2      sess   1 505  6.283754 1.249806e-02     * 0.01229015
# 3      drug   1 505 47.307268 1.800650e-11     * 0.08565389
# 4 sess:drug   1 505 43.392151 1.127202e-10     * 0.07912613
# 
# [[1]]$`Levene's Test for Homogeneity of Variance`
# DFn DFd      SSn      SSd        F            p p<.05
# 1   3 505 13.74217 91.96831 25.15285 3.512871e-15     *
#         
#         
# [[2]]
#       Comparison                      W    p.val        p.adj         sig
# W     Pre CP VEH vs. Pre CP AP5       9878 4.168633e-01 4.168633e-01    
# W1    Post CP VEH vs. Post CP AP5     9692 4.848299e-12 1.939320e-11  ***
# W2    Pre CP VEH vs. Post CP VEH      5753 1.654502e-06 3.309004e-06  ***
# W11   Pre CP AP5 vs. Post CP AP5      9748 5.267440e-07 1.580232e-06  *** #Tail of exc. is significantly SMALLER after CP in the AP5 units. When testing whether it's LARGER: W1 Pre CP AP5 vs. Post CP AP5 W=9748 p=9.999995e-01


#Correcting p values for all these comparisons (100-400ms window and 750-2000ms window):
p.adjust(p=c(6.946541e-02, 5.078548e-04, 3.111524e-12, 1.038793e-02, 
             4.168633e-01, 4.848299e-12, 1.654502e-06, 5.267440e-07), method="holm")
#1.389308e-01 2.031419e-03 2.489219e-11 3.116379e-02 4.168633e-01 3.393809e-11 8.272510e-06 3.160464e-06




#Cue-exc units only
plotBoxplotPrePostCP(experiment="Exp 4 Cue exc only", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                     comp=c("VEH", "AP5"), factors=list(drug=c("VEH", "AP5")), removeCPsess = T,
                     graphFolder=MixedGraphFolder, dataProcess="Zscores", ANOVA=T, 
                     correctOnly=FALSE, cueExcOnly=TRUE, color=colindx, sessFromCP = c(-5, 5), 
                     yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                     removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS_VEH, morethanIQR=T)


### CUE-EXCITED UNITS ONLY
# [[1]]
# [[1]]$`ANOVA`
# Effect DFn DFd           F            p p<.05         ges
# 2      sess   1 208 31.49301308 6.343454e-08     * 0.131498672
# 3      drug   1 208  0.08869364 7.661418e-01       0.000426230
# 4 sess:drug   1 208  0.29549048 5.873055e-01       0.001418612
# 
# [[1]]$`Levene's Test for Homogeneity of Variance`
# DFn DFd      SSn      SSd        F           p p<.05
# 1   3 208 51.99183 576.0997 6.257193 0.000436711     *
#         
#         
#         [[2]]
# Comparison    W        p.val        p.adj sig
# W     Pre CP VEH vs. Pre CP AP5  715 6.672550e-01 6.672550e-01    
# W1  Post CP VEH vs. Post CP AP5 2513 1.150834e-01 2.301668e-01    
# W2   Pre CP VEH vs. Post CP VEH  826 4.324871e-08 1.729948e-07 ***
# W11  Pre CP AP5 vs. Post CP AP5  492 2.336408e-03 7.009225e-03  **

#Cue-exc units only, TAIL (750-2000ms)
plotBoxplotPrePostCP(experiment="Exp 4 Cue exc only 750-2000ms", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                     comp=c("VEH", "AP5"), factors=list(drug=c("VEH", "AP5")), removeCPsess = T,
                     graphFolder=MixedGraphFolder, dataProcess="Zscores", ANOVA=T,
                     correctOnly=FALSE, cueExcOnly=TRUE, color=colindx, sessFromCP = c(-5, 5), 
                     yAxMinZ = -2, yAxMaxZ = 3, yAxMaxRaw = 10, WdwStart=750, WdwEnd=2000, 
                     removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS_VEH, morethanIQR=T)

# [[1]]
# [[1]]$`ANOVA`
# Effect DFn DFd         F            p         p<.05        ges
# 2      sess   1 208  2.123191 1.465920e-01            0.01010450
# 3      drug   1 208 16.031468 8.667462e-05     *      0.07155900
# 4 sess:drug   1 208  5.910995 1.589365e-02     *      0.02763296
# 
# [[1]]$`Levene's Test for Homogeneity of Variance`
# DFn DFd      SSn      SSd        F           p p<.05
# 1   3 208 6.436823 45.00502 9.916369 3.86575e-06     *
#         
#         
#         [[2]]
# Comparison    W        p.val        p.adj sig
# W     Pre CP VEH vs. Pre CP AP5  799 9.368682e-02 9.368682e-02    
# W1  Post CP VEH vs. Post CP AP5 3306 1.192913e-06 4.771653e-06 ***
# W2   Pre CP VEH vs. Post CP VEH 1352 3.165603e-03 9.496809e-03  **
# W11  Pre CP AP5 vs. Post CP AP5  994 3.339285e-02 6.678571e-02    




plotBoxplotPrePostCP(experiment="Exp 4 VEH", masterDF=list(masterDF_DS_VEH_Tone, masterDF_DS_VEH_Light), 
                     comp=c("Tone S+", "Light S+"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                     correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-5, 5), 
                     yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                     removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)

plotBoxplotPrePostCP(experiment="Exp 4 AP5", masterDF=list(masterDF_DS_AP5_Tone, masterDF_DS_AP5_Light), 
                     comp=c("Tone S+", "Light S+"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                     correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-5, 5), 
                     yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                     removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)


### ENTRY DATA

#All units: VEH vs. AP5 before and after CP
plotBoxplotPrePostCP(experiment="Exp 4 DS entry", masterDF=list(masterDF_DSEntry_VEH, masterDF_DSEntry_AP5), 
                     comp=c("VEH", "AP5"), factors=list(drug=c("VEH", "AP5")),  ANOVA=T,
                     graphFolder=MixedGraphFolder, dataProcess="Zscores", removeCPsess = T,
                     correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-5, 5), 
                     yAxMinZ = -2, yAxMaxZ = 2, yAxMaxRaw = 10, WdwStart=0, WdwEnd=500, 
                     removeOutliers=F, imgFormat="pdf", neudata=allNeuronsEntryDS_VEH, morethanIQR=T)

# [[1]]$`ANOVA`
# Effect DFn DFd         F           p p<.05         ges
# 2      sess   1 505 8.2616703 0.004219814     * 0.016096410
# 3      comp   1 505 0.1330268 0.715466876       0.000263350
# 4 sess:comp   1 505 1.3891336 0.239106502       0.002743214
# 
# [[1]]$`Levene's Test for Homogeneity of Variance`
# DFn DFd       SSn      SSd         F         p p<.05
# 1   3 505 0.1478786 51.47684 0.4835748 0.6938394      
# 
# 
# [[2]]
# Comparison    W       p.val      p.adj sig
# W     Pre CP VEH vs. Pre CP AP5 8828 0.090870261 0.27261078    
# W1  Post CP VEH vs. Post CP AP5 6602 0.305757024 0.44126789    
# W2   Pre CP VEH vs. Post CP VEH 6975 0.014998186 0.05999274   
# W11  Pre CP AP5 vs. Post CP AP5 6727 0.220633945 0.44126789    


plotBoxplotPrePostCP(experiment="Exp 4 NS entry", masterDF=list(masterDF_NSEntry_VEH, masterDF_NSEntry_AP5), 
                     comp=c("VEH", "AP5"), factors=list(drug=c("VEH", "AP5")),  ANOVA=T,
                     graphFolder=MixedGraphFolder, dataProcess="Zscores", removeCPsess = T,
                     correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-5, 5), 
                     yAxMinZ = -2, yAxMaxZ = 2, yAxMaxRaw = 10, WdwStart=0, WdwEnd=500, 
                     removeOutliers=F, imgFormat="pdf", neudata=allNeuronsEntryDS_VEH, morethanIQR=T)

# [[1]]$`ANOVA`
#       Effect DFn DFd           F          p p<.05          ges
# 2      sess   1 504  0.01933025 0.889479401       0.0000383522
# 3      comp   1 504 10.06956404 0.001599413     * 0.0195879405
# 4 sess:comp   1 504  0.69535974 0.404742329       0.0013777811
# 
# [[1]]$`Levene's Test for Homogeneity of Variance`
# DFn DFd        SSn     SSd          F        p p<.05
# 1   3 504 0.01868257 51.1336 0.06138179 0.980084      
# 
# 
# [[2]]
# Comparison    W       p.val       p.adj sig
# W     Pre CP VEH vs. Pre CP AP5 7630 0.01531513  0.06126052  
# W1  Post CP VEH vs. Post CP AP5 5557 0.052439010 0.157317031    
# W2   Pre CP VEH vs. Post CP VEH 8399 0.367614845 0.414021413    
# W11  Pre CP AP5 vs. Post CP AP5 7576 0.207010706 0.414021413    

plotBoxplotPrePostCP(experiment="Exp 4 ITI entry", masterDF=list(masterDF_ITIEntry_VEH, masterDF_ITIEntry_AP5), 
                     comp=c("VEH", "AP5"), factors=list(drug=c("VEH", "AP5")),  ANOVA=T,
                     graphFolder=MixedGraphFolder, dataProcess="Zscores", removeCPsess = T,
                     correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-5, 5), 
                     yAxMinZ = -2, yAxMaxZ = 2, yAxMaxRaw = 10, WdwStart=0, WdwEnd=500, 
                     removeOutliers=F, imgFormat="pdf", neudata=allNeuronsEntryDS_VEH, morethanIQR=T)

# [[1]]$`ANOVA`
# Effect DFn DFd          F         p p<.05          ges
# 2      sess   1 505 1.93014254 0.1653551       3.807512e-03
# 3      comp   1 505 1.89263014 0.1695154       3.733789e-03
# 4 sess:comp   1 505 0.01000739 0.9203547       1.981623e-05
# 
# [[1]]$`Levene's Test for Homogeneity of Variance`
# DFn DFd       SSn      SSd         F         p p<.05
# 1   3 505 0.6290748 147.7049 0.7169312 0.5421935      

#       Comparison                      W         p.val     p.adj sig
#W      Pre CP VEH vs. Pre CP AP5       9043 0.15417585 0.4625276    
#W1     Post CP VEH vs. Post CP AP5     5704 0.09311757 0.3724703    
#W2     Pre CP VEH vs. Post CP VEH      9277 0.16314670 0.4625276    
#W11    Pre CP AP5 vs. Post CP AP5      7439 0.28758131 0.4625276    

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






#################################

#Same thing but by session from CP instead of trial in order. To do that, I'm going to change the CP data to 1
LearnerIdx <- !is.na(CPdata$CP)
CPdata$CP[LearnerIdx] <- 1
CPdata$CPsess[LearnerIdx] <- 1

#VEH side
masterDF_DS_VEH_From1 <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS_VEH, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, 
                                      BLneudata=allNeuronsDS_VEH)

masterDF_NS_VEH_From1 <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsNS_VEH, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, 
                                      BLneudata=allNeuronsNS_VEH)


plotBoxplotfromSessCP(experiment="Exp 4 VEH by session", masterDF=list(masterDF_DS_VEH_From1, masterDF_NS_VEH_From1), 
                      comp=c("S+", "S-"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(0, 5), 
                      yAxMinZ = -2, yAxMaxZ = 6, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)

PostCueFR_DSvsNS_fromCP_VEH_BySess <- compareDSvsNSfromCP(masterDF=list(masterDF_DS_VEH_From1, masterDF_NS_VEH_From1), cueExcOnly = F,
                                                   trialBinSize=35, event="cue", correctOnly=F, paired=F, 
                                                   WdwStart=100, WdwEnd=400, capped=T, capValue=c(-105, 105), 
                                                   dataProcess="Zscores")

PostCueFR_DSvsNS_fromCP_VEH_BySess$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP_VEH_BySess$p, method="holm") 
PostCueFR_DSvsNS_fromCP_VEH_BySess$sig <- giveStars(PostCueFR_DSvsNS_fromCP_VEH_BySess$p.adj)

# bin           V            p          n        p.adj sig
# -105 to -70   1331 3.860282e-04       59 3.860282e-04 ***
# -70 to -35    3600 5.396597e-06       97 1.618979e-05 ***
# -35 to 0      3302 1.134632e-08       88 6.807792e-08 ***
# 0 to 35       3658 1.584687e-07       95 7.923434e-07 ***
# 35 to 70      2447 8.004642e-07       77 3.201857e-06 ***
# 70 to 105     770 1.281043e-05        42 2.562085e-05 ***


cumulative_CP(Exp="Exp4 VEH", CPdata=CPdata, numSess=6, byAnimal=F, byNeuron=T, 
              nexdata=allNeuronsDS_VEH$nexdata, graphFolder=MixedGraphFolder)


#AP5 side
masterDF_DS_AP5_From1 <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS_AP5, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, 
                                      BLneudata=allNeuronsDS_AP5)

masterDF_NS_AP5_From1 <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsNS_AP5, 
                                      CPvector=CPdata$CP, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                      BLduration=2, sessionCPperRat=CPdata$CPsess, 
                                      BLneudata=allNeuronsNS_AP5)


plotBoxplotfromSessCP(experiment="Exp 4 AP5 by session", masterDF=list(masterDF_DS_AP5_From1, masterDF_NS_AP5_From1), 
                      comp=c("S+", "S-"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(0, 5), 
                      yAxMinZ = -2, yAxMaxZ = 6, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)


PostCueFR_DSvsNS_fromCP_AP5_BySess <- compareDSvsNSfromCP(masterDF=list(masterDF_DS_AP5_From1, masterDF_NS_AP5_From1), cueExcOnly = F,
                                                          trialBinSize=35, event="cue", correctOnly=F, paired=T, 
                                                          WdwStart=100, WdwEnd=400, capped=T, capValue=c(-105, 105), 
                                                          dataProcess="Zscores")

PostCueFR_DSvsNS_fromCP_AP5_BySess$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP_AP5_BySess$p, method="holm") 
PostCueFR_DSvsNS_fromCP_AP5_BySess$sig <- giveStars(PostCueFR_DSvsNS_fromCP_AP5_BySess$p.adj)


 
cumulative_CP(Exp="Exp4 AP5", CPdata=CPdata, numSess=6, byAnimal=F, byNeuron=T, 
              nexdata=allNeuronsDS_AP5$nexdata, graphFolder=MixedGraphFolder)


save(masterDF_DS_VEH_From1, file=paste(dataForRdir, "masterDF_DS_VEH_From1.rdat", sep=""))
save(masterDF_NS_VEH_From1, file=paste(dataForRdir, "masterDF_NS_VEH_From1.rdat", sep=""))
save(masterDF_DS_AP5_From1, file=paste(dataForRdir, "masterDF_DS_AP5_From1.rdat", sep=""))
save(masterDF_NS_AP5_From1, file=paste(dataForRdir, "masterDF_NS_AP5_From1.rdat", sep=""))


### All units VEH vs AP5 side by bin of trials with respect to CP: BASELINE PERIOD (2s Window right before the cue). not by bin
PostCueFR_DS_From1_VEHvsAP5_BL <- compareVEHvsAP5fromCP(masterDF=list(masterDF_DS_VEH_From1, masterDF_DS_AP5_From1), trialBinSize=210, 
                                                         event="cue", correctOnly=F, cueExcOnly = F, paired=F, WdwStart=-2000, WdwEnd=0, 
                                                         capped=T, capValue=c(1, 210), dataProcess="Zscores")
PostCueFR_DS_From1_VEHvsAP5_BL$p.adj <- p.adjust(PostCueFR_DS_From1_VEHvsAP5_BL$p, method="holm") 
PostCueFR_DS_From1_VEHvsAP5_BL$sig <- giveStars(PostCueFR_DS_From1_VEHvsAP5_BL$p.adj)

#         bin           n           W      p        p.adj sig
# W       1 to NA       337 vs. 305 45027 0.2590486 0.2590486    


##########################################################################
### 20 TRIAL BINS
##########################################################################

#Post cue VEH
VEH_20bin <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 4 VEH", 
                                masterDF=list(masterDF_DS_VEH, masterDF_NS_VEH), 
                                graphFolder=MixedGraphFolder, trialBinSize=20, WdwStart=100, WdwEnd=400, 
                                dataProcess="Zscores", correctOnly=FALSE, points=F, lines=F, 
                                color=c(colindx[1], colindxB[1]), legLabels=c("S+", "S-"), 
                                yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, 
                                capValue=c(-120, 120), cueExcOnly=F, neudata=allNeuronsDS_VEH)


uniqueBin <- unique(VEH_20bin[[1]]$bin)

Diff_DSNS_VEH_20bins <- do.call("rbind", lapply(seq(1, length(uniqueBin)), function(i){
        DSbinSel <- VEH_20bin[[1]][VEH_20bin[[1]]$bin==uniqueBin[i], ]
        NSbinSel <- VEH_20bin[[2]][VEH_20bin[[2]]$bin==uniqueBin[i], ]
        
        matchUnitIdx <- match(DSbinSel$unitIdx, NSbinSel$unitIdx)
        NSfiring <- NSbinSel$byUnitFR[matchUnitIdx]
        
        DSbinSel$byUnitFR_NS <- NSfiring
        DSbinSel$DiffFR <- DSbinSel$byUnitFR-NSfiring
        
        DSbinSel
})
)

#This one compares S+ firing rate among bins only (and make a boxplot of it)
compareBins(data=Diff_DSNS_VEH_20bins, cueExcOnly=F, color=colindx[1], ymin=-2, ymax=8, 
            graphFolder=MixedGraphFolder, experiment="Exp 4 VEH", points=F, 
            comparisons=list(c(3, 6), c(6, 9)), 
            cue="S+")
# [[1]]
# idx       bin           n                        comp      W        p.val        p.adj sig
# W     1   1 vs. 3  63 vs. 128 -120 to -100 vs. -80 to -60 3283.5 3.809391e-02 1.904696e-01    
# W5    2   2 vs. 4   72 vs. 72  -100 to -80 vs. -60 to -40 1768.0 1.076802e-03 9.691219e-03  **
# W1    3   3 vs. 5 128 vs. 126   -80 to -60 vs. -40 to -20 5079.0 1.964275e-06 1.964275e-05 ***
# W11   4   4 vs. 6   72 vs. 69     -60 to -40 vs. -20 to 0 1763.0 8.894694e-03 7.115755e-02    
# W2    5   5 vs. 7 126 vs. 130      -40 to -20 vs. 0 to 20 6366.0 1.248942e-02 7.493653e-02    
# W21   6   6 vs. 8   69 vs. 74       -20 to 0 vs. 20 to 40 2165.0 1.886320e-01 5.658960e-01    
# W3    7   7 vs. 9  130 vs. 73        0 to 20 vs. 40 to 60 3540.0 9.386944e-03 7.115755e-02    
# W31   8  8 vs. 10   74 vs. 35       20 to 40 vs. 60 to 80  949.0 4.093962e-02 1.904696e-01    
# W4    9  9 vs. 11   73 vs. 16      40 to 60 vs. 80 to 100  541.0 7.009865e-01 1.000000e+00    
# W41  10 10 vs. 12   35 vs. 16     60 to 80 vs. 100 to 120  281.0 8.783626e-01 1.000000e+00    
# 
# [[2]]
# bins          n              comparison    W        p.val        p.adj sig
# W  3 vs. 6 128 vs. 69 -80 to -60 vs. -20 to 0 2246 1.641497e-07 3.282995e-07 ***
# W1 6 vs. 9  69 vs. 73   -20 to 0 vs. 40 to 60 1709 4.590662e-03 4.590662e-03  **



compareBins(data=Diff_DSNS_VEH_20bins, cueExcOnly=T, color=colindx[1], ymin=-2, ymax=8, 
            graphFolder=MixedGraphFolder, experiment="Exp 4 VEH", points=F, 
            comparisons=list(c(3, 6), c(6, 10)), 
            cue="S+")


# [[1]]
# idx       bin         n                        comp      W        p.val        p.adj sig
# W     1   1 vs. 3 15 vs. 41 -120 to -100 vs. -80 to -60  250.5 1.479252e-01 0.5917007147    
# W5    2   2 vs. 4 18 vs. 29  -100 to -80 vs. -60 to -40  101.0 1.449676e-04 0.0013047085  **
# W1    3   3 vs. 5 41 vs. 63   -80 to -60 vs. -40 to -20  673.0 1.971785e-05 0.0001971785 ***
# W11   4   4 vs. 6 29 vs. 38     -60 to -40 vs. -20 to 0  421.0 5.075893e-02 0.3553124798    
# W2    5   5 vs. 7 63 vs. 86      -40 to -20 vs. 0 to 20 2308.0 6.190640e-02 0.3714383717    
# W21   6   6 vs. 8 38 vs. 53       -20 to 0 vs. 20 to 40  966.0 3.722393e-01 0.7444786444    
# W3    7   7 vs. 9 86 vs. 56        0 to 20 vs. 40 to 60 2053.5 6.974550e-02 0.3714383717    
# W31   8  8 vs. 10 53 vs. 23       20 to 40 vs. 60 to 80  378.0 4.502121e-03 0.0360169684   *
# W4    9  9 vs. 11  56 vs. 8      40 to 60 vs. 80 to 100  181.0 1.941381e-01 0.5917007147    
# W41  10 10 vs. 12  23 vs. 8     60 to 80 vs. 100 to 120  109.0 7.813806e-01 0.7813806414    
# 
# [[2]]
# bins         n              comparison   W        p.val        p.adj sig
# W  3 vs. 6 41 vs. 38 -80 to -60 vs. -20 to 0 286 2.081033e-07 4.162066e-07 ***
# W1 6 vs. 9 38 vs. 56   -20 to 0 vs. 40 to 60 874 7.214533e-02 7.214533e-02    


#[[2]]
#bins         n              comparison   W        p.val        p.adj sig
#W   3 vs. 6 41 vs. 38 -80 to -60 vs. -20 to 0 286 2.081033e-07 4.162066e-07 ***
#W1 6 vs. 10 38 vs. 23   -20 to 0 vs. 60 to 80 263 4.537571e-03 4.537571e-03  **





#Post cue AP5
AP5_20bin <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 4 AP5", 
                                masterDF=list(masterDF_DS_AP5, masterDF_NS_AP5), 
                                graphFolder=MixedGraphFolder, trialBinSize=20, WdwStart=100, WdwEnd=400, 
                                dataProcess="Zscores", correctOnly=FALSE, points=F, lines=F, 
                                color=c(colindx[1], colindxB[1]), legLabels=c("S+", "S-"), 
                                yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, 
                                capValue=c(-120, 120), cueExcOnly=F, neudata=allNeuronsDS_AP5)


uniqueBin <- unique(AP5_20bin[[1]]$bin)

Diff_DSNS_AP5_20bins <- do.call("rbind", lapply(seq(1, length(uniqueBin)), function(i){
        DSbinSel <- AP5_20bin[[1]][AP5_20bin[[1]]$bin==uniqueBin[i], ]
        NSbinSel <- AP5_20bin[[2]][AP5_20bin[[2]]$bin==uniqueBin[i], ]
        
        matchUnitIdx <- match(DSbinSel$unitIdx, NSbinSel$unitIdx)
        NSfiring <- NSbinSel$byUnitFR[matchUnitIdx]
        
        DSbinSel$byUnitFR_NS <- NSfiring
        DSbinSel$DiffFR <- DSbinSel$byUnitFR-NSfiring
        
        DSbinSel
})
)

#This one compares S+ firing rate among bins only (and make a boxplot of it)
compareBins(data=Diff_DSNS_AP5_20bins, cueExcOnly=F, color=colindx[1], ymin=-2, ymax=8, 
            graphFolder=MixedGraphFolder, experiment="Exp 4 AP5", points=F, 
            comparisons=list(c(3, 6), c(6, 10)), 
            cue="S+")
# [[1]]
# idx       bin          n                        comp      W        p.val       p.adj sig
# W     1   1 vs. 3  29 vs. 97 -120 to -100 vs. -80 to -60 1262.0 0.6017925701 1.000000000    
# W5    2   2 vs. 4  56 vs. 65  -100 to -80 vs. -60 to -40 1408.0 0.1186651063 0.818625211    
# W1    3   3 vs. 5  97 vs. 99   -80 to -60 vs. -40 to -20 2978.5 0.0000844924 0.000844924 ***
# W11   4   4 vs. 6  65 vs. 52     -60 to -40 vs. -20 to 0 1124.0 0.0564624495 0.471436576    
# W2    5   5 vs. 7 99 vs. 101      -40 to -20 vs. 0 to 20 4875.5 0.8904874902 1.000000000    
# W21   6   6 vs. 8  52 vs. 55       -20 to 0 vs. 20 to 40 1167.0 0.6242933245 1.000000000    
# W3    7   7 vs. 9 101 vs. 80        0 to 20 vs. 40 to 60 3716.0 0.6855602334 1.000000000    
# W31   8  8 vs. 10  55 vs. 69       20 to 40 vs. 60 to 80 1739.0 0.7408569105 1.000000000    
# W4    9  9 vs. 11  80 vs. 45      40 to 60 vs. 80 to 100 1289.0 0.0523818417 0.471436576    
# W41  10 10 vs. 12  69 vs. 44     60 to 80 vs. 100 to 120 1148.5 0.1169464587 0.818625211    
# 
# [[2]]
# bins         n              comparison    W       p.val      p.adj sig
# W   3 vs. 6 97 vs. 52 -80 to -60 vs. -20 to 0 1550 0.005525344 0.01105069   *
# W1 6 vs. 10 52 vs. 69   -20 to 0 vs. 60 to 80 1582 0.767163441 0.76716344    



compareBins(data=Diff_DSNS_AP5_20bins, cueExcOnly=T, color=colindx[1], ymin=-2, ymax=8, 
            graphFolder=MixedGraphFolder, experiment="Exp 4 AP5", points=F, 
            comparisons=list(c(3, 6), c(6, 9)), 
            cue="S+")

# [[1]]
# idx       bin         n                        comp     W        p.val       p.adj sig
# W     1   1 vs. 3  6 vs. 23 -120 to -100 vs. -80 to -60  29.0 0.0153593533 0.122874826    
# W5    2   2 vs. 4 13 vs. 16  -100 to -80 vs. -60 to -40 102.0 0.4741943638 1.000000000    
# W1    3   3 vs. 5 23 vs. 30   -80 to -60 vs. -40 to -20 162.0 0.0003906602 0.003906602  **
# W11   4   4 vs. 6 16 vs. 20     -60 to -40 vs. -20 to 0  84.0 0.0074406581 0.066965923    
# W2    5   5 vs. 7 30 vs. 43      -40 to -20 vs. 0 to 20 657.0 0.5553725187 1.000000000    
# W21   6   6 vs. 8 20 vs. 23       -20 to 0 vs. 20 to 40 226.0 0.4664020778 1.000000000    
# W3    7   7 vs. 9 43 vs. 34        0 to 20 vs. 40 to 60 719.0 0.4533248820 1.000000000    
# W31   8  8 vs. 10 23 vs. 30       20 to 40 vs. 60 to 80 335.0 0.4328957112 1.000000000    
# W4    9  9 vs. 11 34 vs. 22      40 to 60 vs. 80 to 100 315.0 0.1646259271 1.000000000    
# W41  10 10 vs. 12 30 vs. 22     60 to 80 vs. 100 to 120 295.5 0.2644287350 1.000000000    
# 
# [[2]]
# bins         n              comparison   W       p.val       p.adj sig
# W  3 vs. 6 23 vs. 20 -80 to -60 vs. -20 to 0 124 0.004595229 0.009190459  **
# W1 6 vs. 9 20 vs. 34   -20 to 0 vs. 40 to 60 338 0.489380943 0.489380943    


# [[2]]
# bins         n              comparison   W       p.val       p.adj sig
# W   3 vs. 6 23 vs. 20 -80 to -60 vs. -20 to 0 124 0.004595229 0.009190459  **
# W1 6 vs. 10 20 vs. 30   -20 to 0 vs. 60 to 80 265 0.249000033 0.249000033    
# 


#### HISTOGRAM
plotFRandCPhistogram(experiment="Exp 4 All units", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), graphFolder=MixedGraphFolder, 
                     trialBinSize=20, dataProcess="Zscores", correctOnly=F, color=c(colindx[1], "darkblue"), 
                     cueExcOnly = F, events=c("S+","S+"),
                     capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 3, 
                     psthmin=0.5, psthmax=2, imgFormat="pdf", neudata=allNeuronsDS_VEH)

plotFRandCPhistogram(experiment="Exp 4 Cue excited", masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), graphFolder=MixedGraphFolder, 
                     trialBinSize=20, dataProcess="Zscores", correctOnly=F, color=c(colindx[1], "darkblue"), 
                     cueExcOnly = T, events=c("S+","S+"),
                     capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 3, 
                     psthmin=0.5, psthmax=2, imgFormat="pdf", neudata=allNeuronsDS_VEH)




##### S+ vs S-

# VEH
PostCueFR_DSvsNS_fromCP_VEH <- compareDSvsNSfromCP(masterDF=list(masterDF_DS_VEH, masterDF_NS_VEH), cueExcOnly = F,
                                                   trialBinSize=20, event="cue", correctOnly=F, paired=T, 
                                                   WdwStart=100, WdwEnd=400, capped=T, capValue=c(-120, 120), 
                                                   dataProcess="Zscores")

PostCueFR_DSvsNS_fromCP_VEH$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP_VEH$p, method="holm") 
PostCueFR_DSvsNS_fromCP_VEH$sig <- giveStars(PostCueFR_DSvsNS_fromCP_VEH$p.adj)

# bin    V            p   n        p.adj sig
# V   -120 to -100  136 3.987975e-02  29 1.196393e-01    
# V1   -100 to -80  860 6.024391e-02  52 1.204878e-01    
# V2    -80 to -60 3278 1.136270e-01 107 1.204878e-01    
# V3    -60 to -40 1925 3.287149e-05  70 1.972290e-04 ***
# V4    -40 to -20 6165 7.806722e-11 121 7.806722e-10 ***
# V5      -20 to 0 1564 1.517169e-07  59 1.213735e-06 ***
# v6       0 to 20 6697 9.898376e-13 124 1.187805e-11 ***
# V7      20 to 40 2437 1.496644e-10  72 1.346979e-09 ***
# V8      40 to 60  977 5.002221e-12  44 5.502443e-11 ***
# V9      60 to 80  489 1.756009e-06  32 1.229206e-05 ***
# V10    80 to 100   85 1.708984e-03  13 8.544922e-03  **
# V11   100 to 120   96 2.014160e-03  14 8.544922e-03  **


# AP5
PostCueFR_DSvsNS_fromCP_AP5 <- compareDSvsNSfromCP(masterDF=list(masterDF_DS_AP5, masterDF_NS_AP5), cueExcOnly = F,
                                                   trialBinSize=20, event="cue", correctOnly=F, paired=T, 
                                                   WdwStart=100, WdwEnd=400, capped=T, capValue=c(-120, 120), 
                                                   dataProcess="Zscores")

PostCueFR_DSvsNS_fromCP_AP5$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP_AP5$p, method="holm") 
PostCueFR_DSvsNS_fromCP_AP5$sig <- giveStars(PostCueFR_DSvsNS_fromCP_AP5$p.adj)

# bin    V            p  n        p.adj sig
# V   -120 to -100  172 2.728348e-01 24 2.728348e-01    
# V1   -100 to -80  908 4.459079e-02 53 8.918157e-02    
# V2    -80 to -60 3188 6.170370e-05 93 3.085185e-04 ***
# V3    -60 to -40 1637 3.465306e-07 61 3.465306e-06 ***
# V4    -40 to -20 3845 1.031251e-10 93 1.237501e-09 ***
# V5      -20 to 0  926 2.889173e-07 45 3.178090e-06 ***
# V6       0 to 20 3428 1.025071e-05 95 8.200571e-05 ***
# V7      20 to 40 1108 2.852451e-06 50 2.567206e-05 ***
# V8      40 to 60  408 5.777781e-04 31 2.311112e-03  **
# V9      60 to 80 1559 7.466593e-04 65 2.311112e-03  **
# V10    80 to 100  761 2.343283e-05 42 1.640298e-04 ***
# V11   100 to 120  717 5.093772e-05 41 3.056263e-04 ***


padjall <- p.adjust(c(PostCueFR_DSvsNS_fromCP_VEH$p, PostCueFR_DSvsNS_fromCP_AP5$p), method="holm") 

PostCueFR_DSvsNS_fromCP_VEH$p.adj <- padjall[1:nrow(PostCueFR_DSvsNS_fromCP_VEH)]
PostCueFR_DSvsNS_fromCP_AP5$p.adj <- padjall[(nrow(PostCueFR_DSvsNS_fromCP_VEH)+1):length(padjall)]



# All units VEH vs AP5 side by bin of trials with respect to CP: 100-400ms after S+
PostCueFR_DS_fromCP_VEHvsAP5_Spike <- compareVEHvsAP5fromCP(masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                                                            trialBinSize=20,  event="cue", correctOnly=F, 
                                                            cueExcOnly = F, paired=F, WdwStart=100, WdwEnd=400, 
                                                            capped=T, capValue=c(-80, 80), dataProcess="Zscores")


PostCueFR_DS_fromCP_VEHvsAP5_Spike$p.adj <- p.adjust(PostCueFR_DS_fromCP_VEHvsAP5_Spike$p, method="holm") 
PostCueFR_DS_fromCP_VEHvsAP5_Spike$sig <- giveStars(PostCueFR_DS_fromCP_VEHvsAP5_Spike$p.adj)
# bin           n    W            p        p.adj sig
# W  -80 to -60  128 vs. 97 5071 3.430035e-02 1.372014e-01    
# W1 -60 to -40   72 vs. 65 2289 2.872609e-01 5.745219e-01    
# W2 -40 to -20  126 vs. 99 5596 5.680816e-01 5.745219e-01    
# W3   -20 to 0   69 vs. 52 1650 1.277980e-01 3.833939e-01    
# W4    0 to 20 130 vs. 101 7210 3.264350e-03 1.958610e-02   *
# W5   20 to 40   74 vs. 55 2233 1.794074e-02 8.970369e-02    
# W6   40 to 60   73 vs. 80 3754 9.701067e-06 7.760854e-05 ***
# W7   60 to 80   35 vs. 69 1490 8.601614e-04 6.021130e-03  **


PostCueFR_DS_fromCP_VEHvsAP5_Spike <- compareVEHvsAP5fromCP(masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                                                            trialBinSize=20,  event="cue", correctOnly=F, 
                                                            cueExcOnly = T, paired=F, WdwStart=100, WdwEnd=400, 
                                                            capped=T, capValue=c(-80, 80), dataProcess="Zscores")


PostCueFR_DS_fromCP_VEHvsAP5_Spike$p.adj <- p.adjust(PostCueFR_DS_fromCP_VEHvsAP5_Spike$p, method="holm") 
PostCueFR_DS_fromCP_VEHvsAP5_Spike$sig <- giveStars(PostCueFR_DS_fromCP_VEHvsAP5_Spike$p.adj)
# bin         n    W          p     p.adj sig
# W  -80 to -60 41 vs. 23  381 0.10482776 0.5241388    
# W1 -60 to -40 29 vs. 16  324 0.01446056 0.1081896    
# W2 -40 to -20 63 vs. 30  912 0.39469484 0.8141336    
# W3   -20 to 0 38 vs. 20  418 0.27137786 0.8141336    
# W4    0 to 20 86 vs. 43 2052 0.15583587 0.6233435    
# W5   20 to 40 53 vs. 23  640 0.63702343 0.8141336    
# W6   40 to 60 56 vs. 34 1127 0.07322047 0.4393228    
# W7   60 to 80 23 vs. 30  468 0.01352371 0.1081896  



PostCueFR_DS_fromCP_VEHvsAP5_Tail <- compareVEHvsAP5fromCP(masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                                                            trialBinSize=20,  event="cue", correctOnly=F, 
                                                            cueExcOnly = F, paired=F, WdwStart=750, WdwEnd=2000, 
                                                            capped=T, capValue=c(-80, 80), dataProcess="Zscores")


PostCueFR_DS_fromCP_VEHvsAP5_Tail$p.adj <- p.adjust(PostCueFR_DS_fromCP_VEHvsAP5_Tail$p, method="holm") 
PostCueFR_DS_fromCP_VEHvsAP5_Tail$sig <- giveStars(PostCueFR_DS_fromCP_VEHvsAP5_Tail$p.adj)
# bin           n    W            p        p.adj sig
# W  -80 to -60  128 vs. 97 6214 2.662670e-01 7.988009e-01    
# W1 -60 to -40   72 vs. 65 2156 4.836167e-01 9.227237e-01    
# W2 -40 to -20  125 vs. 98 5568 4.613619e-01 9.227237e-01    
# W3   -20 to 0   69 vs. 52 2023 3.314181e-04 1.325672e-03  **
# W4    0 to 20 130 vs. 101 7932 1.004631e-05 6.027789e-05 ***
# W5   20 to 40   74 vs. 55 2648 1.141973e-05 6.027789e-05 ***
# W6   40 to 60   58 vs. 38 1453 3.575525e-06 2.502868e-05 ***
# W7   60 to 80   35 vs. 69 1831 6.058191e-09 4.846553e-08 ***



PostCueFR_DS_fromCP_VEHvsAP5_Tail <- compareVEHvsAP5fromCP(masterDF=list(masterDF_DS_VEH, masterDF_DS_AP5), 
                                                           trialBinSize=20,  event="cue", correctOnly=F, 
                                                           cueExcOnly = T, paired=F, WdwStart=750, WdwEnd=2000, 
                                                           capped=T, capValue=c(-80, 80), dataProcess="Zscores")


PostCueFR_DS_fromCP_VEHvsAP5_Tail$p.adj <- p.adjust(PostCueFR_DS_fromCP_VEHvsAP5_Tail$p, method="holm") 
PostCueFR_DS_fromCP_VEHvsAP5_Tail$sig <- giveStars(PostCueFR_DS_fromCP_VEHvsAP5_Tail$p.adj)
# bin         n    W            p        p.adj sig
# W  -80 to -60 41 vs. 23  499 3.538534e-01 7.077068e-01    
# W1 -60 to -40 29 vs. 16  237 5.512983e-01 7.077068e-01    
# W2 -40 to -20 63 vs. 30  823 1.590051e-01 4.770153e-01    
# W3   -20 to 0 38 vs. 20  485 4.371546e-02 2.185773e-01    
# W4    0 to 20 86 vs. 43 2149 6.728163e-02 2.691265e-01    
# W5   20 to 40 53 vs. 23  855 2.801184e-03 1.960829e-02   *
# W6   40 to 60 46 vs. 17  542 9.382288e-03 5.629373e-02    
# W7   60 to 80 23 vs. 30  592 1.280599e-06 1.024479e-05 ***





###################################################################################
### RELATIONSHIP BETWEEN CUE-evoked AND ENTRY-evoked FIRING RATE                ###
###################################################################################

### VEH SIDE
VEH_toplotTC <- megaplotTC(experiment="Exp4 VEH Cue and Entry", 
                           data=list(allNeuronsDS_VEH, allNeuronsEntryDS_VEH, 
                                     allNeuronsEntryNS_VEH, allNeuronsEntryITI_VEH), 
                       graphmin=500, graphmax=2000, cuewinmin=100, cuewinmax=400,
                       colpalette="Rainbow", minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, 
                       title="Entry", ZcolLabels=c("ZDS", "ZDS entry", "Z NS entry", "Z ITI entry"), arrangeBy=c("ZDS"))

save(VEH_toplotTC, file=paste(dataForRdir, "VEH_toplotTC.rdat", sep=""))

#Use the object I just created to compare FR after S+ vs. after entry. If I change graphmin and graphmax in the previous function,
# I need to change that inside the next function as well
CUEvsENTRYfiring(cuedata=VEH_toplotTC[[1]], entrydata=VEH_toplotTC[[2]], graphFolder=MixedGraphFolder,
                 title="VEH S+ onset vs. S+ entry", cueExcOnly=F, cuemin=100, cuemax=400, collapseCP=TRUE,
                 entrymin=0, entrymax=1500)

CUEvsENTRYfiring(cuedata=VEH_toplotTC[[1]], entrydata=VEH_toplotTC[[3]], unitVScp=unitsDF,graphFolder=MixedGraphFolder, 
                 title="VEH S+ onset vs. S- entry", cueExcOnly=F, cuemin=100, cuemax=400, collapseCP=TRUE,
                 entrymin=0, entrymax=1500)

CUEvsENTRYfiring(cuedata=VEH_toplotTC[[1]], entrydata=VEH_toplotTC[[4]], unitVScp=unitsDF,graphFolder=MixedGraphFolder,
                 title="VEH S+ onset vs. ITI entry", cueExcOnly=F, cuemin=100, cuemax=400, collapseCP=TRUE,
                 entrymin=0, entrymax=1500)



### AP5 SIDE
AP5_toplotTC <- megaplotTC(experiment="Exp4 AP5 Cue and Entry", 
                           data=list(allNeuronsDS_AP5, allNeuronsEntryDS_AP5, 
                                     allNeuronsEntryNS_AP5, allNeuronsEntryITI_AP5), 
                           graphmin=500, graphmax=2000, cuewinmin=100, cuewinmax=400,
                           colpalette="Rainbow", minFR=-3, maxFR=3, graphFolder=MixedGraphFolder, 
                           title="Entry", ZcolLabels=c("ZDS", "ZDS entry", "Z NS entry", "Z ITI entry"), 
                           arrangeBy=c("ZDS"))

save(AP5_toplotTC, file=paste(dataForRdir, "AP5_toplotTC.rdat", sep=""))

#Use the object I just created to compare FR after S+ vs. after entry. If I change graphmin and graphmax in the previous function,
# I need to change that inside the next function as well
CUEvsENTRYfiring(cuedata=AP5_toplotTC[[1]], entrydata=AP5_toplotTC[[2]], graphFolder=MixedGraphFolder,
                 title="AP5 S+ onset vs. S+ entry", cueExcOnly=F, cuemin=100, cuemax=400, collapseCP=TRUE,
                 entrymin=0, entrymax=1500)

CUEvsENTRYfiring(cuedata=AP5_toplotTC[[1]], entrydata=AP5_toplotTC[[3]], unitVScp=unitsDF,graphFolder=MixedGraphFolder, 
                 title="AP5 S+ onset vs. S- entry", cueExcOnly=F, cuemin=100, cuemax=400, collapseCP=TRUE,
                 entrymin=100, entrymax=1500)

CUEvsENTRYfiring(cuedata=AP5_toplotTC[[1]], entrydata=AP5_toplotTC[[4]], unitVScp=unitsDF,graphFolder=MixedGraphFolder,
                 title="AP5 S+ onset vs. ITI entry", cueExcOnly=F, cuemin=100, cuemax=400, collapseCP=TRUE,
                 entrymin=100, entrymax=1500)


# Percentage of S+, S- and ITI entry-excited and inhibited neurons before and after CP and comparison of FR among these populations of units
# I'm going to use VEH_toplotTC [[2:4]], which represent FR aligned to S+, S- and ITI entry. 
# The % of exc/inh are also adjusted to each one of those events. Criteria in neurohist.

DSentry_VEHvsAP5 <- entry.perc.and.firing(data=list(VEH_toplotTC[[2]], AP5_toplotTC[[2]]), collapseCP=T, 
                     groupLabels=c("VEH", "AP5"), color=colindx, eventName="S+ Entry", 
                     funcdirect=funcdirect, graphFolder=MixedGraphFolder)


save(DSentry_VEHvsAP5, file=paste(dataForRdir, "DSentry_VEHvsAP5.rdat", sep=""))

#Adjust p values. 12 tests. It goes in order from to the topleft corner of the graph down the first column, then second then third column down.
p.vals <- unlist(sapply(seq(1, length(DSentry_VEHvsAP5)), function(x){
        unlist(sapply(seq(1, length(DSentry_VEHvsAP5[[x]])), function(y){
                DSentry_VEHvsAP5[[x]][[y]]$p.value
        }))
}))

p.adjust(p=p.vals, method="holm")
#[1] 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00 1.000000e+00
# [9] 3.618792e-06 3.000154e-01 5.311355e-01 1.000000e+00

NSentry_VEHvsAP5 <- entry.perc.and.firing(data=list(VEH_toplotTC[[3]], AP5_toplotTC[[3]]), collapseCP=T, 
                                          groupLabels=c("VEH", "AP5"), color=colindx, eventName="S- Entry", 
                                          funcdirect=funcdirect, graphFolder=MixedGraphFolder)

#Adjust p values. 12 tests. It goes in order from to the topleft corner of the graph down the first column, then second then third column down.
p.vals <- unlist(sapply(seq(1, length(NSentry_VEHvsAP5)), function(x){
        unlist(sapply(seq(1, length(NSentry_VEHvsAP5[[x]])), function(y){
                p <- NSentry_VEHvsAP5[[x]][[y]]
                if(is.na(p)){pval <- NA} else {p <- p$p.value}
                p
        }))
}))

p.adjust(p=p.vals, method="holm")
#[1] 1.000000e+00 1.403094e-01 1.000000e+00 1.000000e+00 6.650364e-01 1.000000e+00           NA 2.525046e-01
#[9] 2.599162e-05 1.000000e+00           NA 4.086509e-01

ITIentry_VEHvsAP5 <- entry.perc.and.firing(data=list(VEH_toplotTC[[4]], AP5_toplotTC[[4]]), collapseCP=T, 
                                          groupLabels=c("VEH", "AP5"), color=colindx, eventName="ITI Entry", 
                                          funcdirect=funcdirect, graphFolder=MixedGraphFolder)

p.vals <- unlist(sapply(seq(1, length(ITIentry_VEHvsAP5)), function(x){
        unlist(sapply(seq(1, length(ITIentry_VEHvsAP5[[x]])), function(y){
                p <- ITIentry_VEHvsAP5[[x]][[y]]
                if(is.na(p)){pval <- NA} else {p <- p$p.value}
                p
        }))
})) 

p.adjust(p=p.vals, method="holm")
#[1] 0.53699759 0.53699759 0.53699759         NA 1.00000000 1.00000000 0.50242466         NA 0.49566112
#[10] 0.01306523 0.01607137 0.53699759






########################################################################################
### BOXPLOT WITH FR RIGHT BEFORE ENTRY IN LONG-LATENCY ENTRIES
########################################################################################


PreEntry5sLat <- plotFRBoxPlotandCP_Entry(cue=c("S+", "S+"), experiment="Exp 4 2s Pre entry FR when latency over 5s from CP CUE EXC", 
                                          masterDF=list(masterDF_DSEntry_VEH, masterDF_DSEntry_AP5), graphFolder=MixedGraphFolder, 
                                          trialBinSize=35, WdwStart=-2000, WdwEnd=0, dataProcess="Zscores", minLat=5000,
                                          correctOnly=FALSE, points=F, lines=F, color=colindx, 
                                          legLabels=c("VEH", "AP5"), yAxMinZ=-1, yAxMaxZ=4, yAxMaxRaw=7, capped=T, 
                                          capValue=c(-105, 105), cueExcOnly=T, BLneudata=allNeuronsDS_VEH)





#VEH vs. AP5
VEHvsAP5_PreEntry <- do.call("rbind", 
                           lapply(seq(1, length(unique(PreEntry5sLat[[1]]$correspBins))), function(i){
                                   bin <- unique(PreEntry5sLat[[1]]$correspBins)[i]
                                   
                                   datSel1 <- PreEntry5sLat[[1]][PreEntry5sLat[[1]]$correspBins==i, ]
                                   datSel2 <- PreEntry5sLat[[2]][PreEntry5sLat[[2]]$correspBins==i, ]
                                   
                                   comp <- wilcox.test(x=datSel1$byUnitFR, y=datSel2$byUnitFR, alternative = "greater")
                                   
                                   data.frame(W=comp$statistic, p=comp$p.value)
                           })
)

# Comparing FR in the 2s before entry VEH vs AP5 side in S+ trials in which latency to enter was 5 or more seconds.
# W            p.adj
# W  1622 0.6817585471
# W1 2583 0.2246616741
# W2 3449 0.0328868795
# W3 2740 0.0404253619
# W4 1477 0.0009510481
# W5  709 0.0001469833

p.adjust(p=c(VEHvsAP5_PreEntry$p), method="holm")
# [1] 0.6817585471 0.2246616741 0.0328868795 0.0404253619
# [5] 0.0009510481 0.0001469833


#across bin comp VEH SIDE
comp1 <- PreEntry5sLat[[1]][PreEntry5sLat[[1]]$correspBins==1, ]$byUnitFR
comp2 <- PreEntry5sLat[[1]][PreEntry5sLat[[1]]$correspBins==3, ]$byUnitFR
comp3 <- PreEntry5sLat[[1]][PreEntry5sLat[[1]]$correspBins==5, ]$byUnitFR

wilcox.test(comp1, comp2, alternative="less") #W = 219, p-value = 0.1732
wilcox.test(comp2, comp3, alternative="less") #W = 103, p-value = 0.008814

#across bin comp
comp1 <- PreEntry5sLat[[2]][PreEntry5sLat[[2]]$correspBins==1, ]$byUnitFR
comp2 <- PreEntry5sLat[[2]][PreEntry5sLat[[2]]$correspBins==3, ]$byUnitFR
comp3 <- PreEntry5sLat[[2]][PreEntry5sLat[[2]]$correspBins==5, ]$byUnitFR

wilcox.test(comp1, comp2, alternative="less") #W = 129, p-value = 0.7568
wilcox.test(comp2, comp3, alternative="less") #W = 88, p-value = 0.7694


#Number of units per condition (different bc less >5s entries in S+ than S-). That's why I treated these as different samples in the Wilcoxon test.
nunitsVEH <- sapply(seq(1, 6), function(i){
        length(PreEntry5sLat[[1]][PreEntry5sLat[[1]]$correspBins==i, ]$byUnitFR)
})
#[1] 21 22 25 28 15  2

nunitsAP5 <- sapply(seq(1, 6), function(i){
        length(PreEntry5sLat[[2]][PreEntry5sLat[[2]]$correspBins==i, ]$byUnitFR)
})
#[1] 15 20 15 11 10  9