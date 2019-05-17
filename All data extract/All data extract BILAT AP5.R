
setwd("E:/Dropbox/NMDA/")

#############################################################
### BILATERAL AP5 INFUSIONS + RECORDINGS                  ###
#############################################################

### LOAD IMPORTANT LIBRARIES
install.packages("matrixStats")
install.packages('ez')

library(matrixStats)
library(ez)

setwd("E:/Dropbox/NMDA/")

#Define colors
colindx <- c("#2171b5", "#cb181d") #Strong blue and red
colindxB <- c("#bdd7e7", "#fcae91") #Less strong blue and red
colindxC <- c("#eff3ff", "#fb6a4a") #Even less strong blue and red
colindxD <- c("#6baed6", "#fee5d9") #Lightest blue and red


### DEFINE ALL THE IMPORTANT FOLDERS
funcdirect <- paste(getwd(), "/R functions/", sep="")
datafolder <- paste(getwd(), "/EXP3_Bilateral AP5/MedPC files/", sep="")
dataForRdir <- paste(getwd(), "/EXP3_Bilateral AP5/Data for R/", sep="")
dataForRCumulative <- paste(getwd(), "/EXP3_Bilateral AP5/Data for R cumulative/", sep="")
behGraphFolder <- paste(getwd(), "/EXP3_Bilateral AP5/Graphs/Behavior/", sep="")
CPfuncFolder <- paste(getwd(), '/R functions/Change_Point-master/', sep="")
CPGraphFolder <- paste(getwd(), "/EXP3_Bilateral AP5/Graphs/Behavior/Change point/", sep="")
MixedGraphFolder <- paste(getwd(), "/EXP3_Bilateral AP5/Graphs/Mixed/", sep="")
NEXfiles <- paste(getwd(), "/EXP3_Bilateral AP5/NEX files/", sep="")
IndcumPerfGraphFolder <- paste(getwd(), "/EXP3_Bilateral AP5/Graphs/Behavior/Cumulative perf/", sep="")
PerfRelToCPFolder <- paste(getwd(), "/EXP3_Bilateral AP5/Graphs/Behavior/Beh rel to CP/", sep="")
BySessFolder <- paste(MixedGraphFolder, "By Sess from CP/", sep="")


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
avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), YMIN=-2, YMAX=5, data=list(DStaskAcc, NStaskAcc), 
             cues=c("S+", "S-"), index="Specificity", legendLocation="topleft", y_axis_label="Performance index", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c(colindx[1], "darkblue", "gray50"), data=list(DSrespAll, NSrespAll, ITIrespRatio), 
             cues=c("S+", "S-"), index="Response ratio", legendLocation="topleft", y_axis_label="Response ratio", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), YMIN=2, YMAX=10, data=list(DSlatency, NSlatency), 
             cues=c("S+", "S-"), index="Latency", legendLocation="bottomright", y_axis_label="Latency (s)", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c("gray20", "gray50"), YMIN=2, YMAX=10, data=list(ITIlatency), 
             cues=c("ITI"), index="Latency", legendLocation=-1, y_axis_label="Latency (s)", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c("#0d2955", "#2678fa"), YMIN=-3, YMAX=6, data=list(DStaskAcc, NStaskAcc), 
             cues=c("S+", "S-"), index="Specificity", legendLocation="topleft", y_axis_label="Performance index", 
             behGraphFolder=behGraphFolder, plot=T)

avgPerfByBin(binsize=5, colors=c("black", "gray30"), YMIN=0, YMAX=10, data=list(DSlatency, ITIlatency), 
             cues=c("S+", "ITI"), index="Latency", legendLocation="topleft", y_axis_label="Latency", 
             behGraphFolder=behGraphFolder, plot=T)


### CHANGE POINT: S+ SPECIFICITY
# Change point analysis using the method in Gallistel et al., 2004. I'll use S+ SPECIFICITY as the main performance index. But I also need to look at the results that other performance indexes give me.

CPdata <- CPextract(GallCrit=1.3, plot=T, CPfuncFolder=CPfuncFolder, idx=idx, CPGraphFolder=CPGraphFolder, dataForRdir=dataForRdir, dataForRCumulative=dataForRCumulative)


CPdataExp5 <- CPdata
alldataExp5 <- alldata
csacqidxExp5 <- csacqidx
ratsExp5 <- rats
idxExp5 <- idx

save(csacqidxExp5, file=paste(dataForRdir, "csacqidxExp5.ridx", sep=""))
save(alldataExp5, file=paste(dataForRdir, "alldataExp5.rdat", sep=""))
save(ratsExp5, file=paste(dataForRdir, "ratsExp5.rdat", sep=""))
save(idxExp5, file=paste(dataForRdir, "idxExp5.rdat", sep=""))


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
CPframeExp5 <- rbind(DSspecificity, RealCPbyIndex[sel,])

save(CPframeExp5, file=paste(dataForRdir, "CPframeExp5.rdat", sep="")) #CP values for different behavioral indexes
save(CPdataExp5, file=paste(dataForRdir, "CPdataExp5.rdat", sep="")) #CP values and other details for PERFORMANCE INDEX (S+ specificity)


#CUMULATIVE INDIVIDUAL PERFORMANCE
# #Smart rats
# cumulativeIndGraphs(numSess=6, numCues=35, sessLines=F, smartRat = T, dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=IndcumPerfGraphFolder, imgFormat="pdf", limit=210)
# 
# #Dumb rats
# cumulativeIndGraphs(numSess=6, numCues=35, sessLines=F, smartRat = F, dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=IndcumPerfGraphFolder, imgFormat="pdf", limit=210)

#All rats
cumulativeIndGraphs(numSess=6, numCues=35, sessLines=F, smartRat = "NA", dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=IndcumPerfGraphFolder, imgFormat="pdf", limit=210)


#MEAN PERFORMANCE AS A FUNCTION OF DISTANCE FROM CP
CPdata$CP <- rep(1, nrow(CPdata)); CPdata$CPsess <- rep(1, nrow(CPdata))
PerformanceFromCP(relTrialMin=0, relTrialMax=210, trialBinSize=5, typegraph="both", imgFormat="pdf", 
                  dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=PerfRelToCPFolder, CPdata=CPdata, 
                  idx=idx, rats = rats, alldata = alldata)


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

###TOTAL REWARDS

plot.new()
plot.window(xlim=c(0, 2), ylim=c(0, 210))
rect(xleft=0, xright=1, ybottom=0, ytop = mean(CPdata$TotalRewards), col="gray30")
points(x=rep(0.5, nrow(CPdata)), y=CPdata$TotalRewards, pch=19, col="black")
rect(xleft=1, xright=2, ybottom=0, ytop=mean(CPdata$TotalRewardsShortLat, col="white"))
points(x=rep(1.5, nrow(CPdata)), y=CPdata$TotalRewardsShortLat, pch=19, col="black")

for(i in 1:nrow(CPdata)){
        lines(x=c(0.5, 1.5), y=c(CPdata$TotalRewards[i], CPdata$TotalRewardsShortLat[i]))
}


axis(side=2, las=2)
mtext(side=2, text="Total # rewards", line=2.5, font=2, cex=1.4)
axis(side=1, at=c(0.5, 1.5), labels=c("Total rewards", "# rewards with < 2 latency"), font=2, cex.axis=1.2)


#lEARNERS: TOTAL TRIALS TO CP VS. NUMBER OF REWARDS
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

#### RUN THIS FOR NEURONS IN VEHICLE SIDE
funcdirect <- paste(getwd(), "/R functions/", sep="")
allNeuronsDS <-neuralhist(funcdirect=funcdirect, path=NEXfiles, event=1, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="both")
allNeuronsNS = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=2, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="both")
allNeuronsDSresponded = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=5, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="both")
allNeuronsDSmissed = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=6, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="both")
allNeuronsNSresponded = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=7, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="both")
allNeuronsNSmissed = neuralhist(funcdirect=funcdirect, path=NEXfiles, event=8, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="both")
allNeuronsEntryDS = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="both")
allNeuronsEntryNS = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="both")
allNeuronsEntryITI = neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="both")


#sAVE THESE OBJECTS
save(allNeuronsDS, file=paste(dataForRdir, 'allNeuronsDS.rdat', sep=""))
save(allNeuronsNS, file=paste(dataForRdir, 'allNeuronsNS.rdat', sep=""))
save(allNeuronsDSresponded, file=paste(dataForRdir, 'allNeuronsDSresponded.rdat', sep=""))
save(allNeuronsDSmissed, file=paste(dataForRdir, 'allNeuronsDSmissed.rdat', sep=""))
save(allNeuronsNSresponded, file=paste(dataForRdir, 'allNeuronsNSresponded.rdat', sep=""))
save(allNeuronsNSmissed, file=paste(dataForRdir, 'allNeuronsNSmissed.rdat', sep=""))
save(allNeuronsEntryDS, file=paste(dataForRdir, 'allNeuronsEntryDS.rdat', sep=""))
save(allNeuronsEntryNS, file=paste(dataForRdir, 'allNeuronsEntryNS.rdat', sep=""))
save(allNeuronsEntryITI, file=paste(dataForRdir, 'allNeuronsEntryITI.rdat', sep=""))


##############################
### FR from trial 1 
##############################
#Only 2 out of 5 rats showed a CP and it looks like an artifact (a short interval of >0 performance index is picked up by the algorithm, but is not a real CP.
#So, for further analyses, I'll look at the data from trial 1, and all my code is defined around CP trial, so change CP trial to 1 so that I don't have to rewrite all the code

CPdata$CP <- rep(1, nrow(CPdata))
CPdata$CPsess <- rep(1, nrow(CPdata))
CPvector <- CPdata$CP
sessionCPperRat <- CPdata$CPsess


masterDF_DS_Bil <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS, 
                                   CPvector=CPdata$CP, sessionCPperRat=sessionCPperRat, 
                                   funcdirect=funcdirect, dataForRdir=dataForRdir, BLduration=2, 
                                   BLneudata = allNeuronsDS)

masterDF_NS_Bil <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS, 
                                      CPvector=CPdata$CP, sessionCPperRat=sessionCPperRat, 
                                      funcdirect=funcdirect, dataForRdir=dataForRdir, BLduration=2, 
                                      BLneudata = allNeuronsNS)


save(masterDF_DS_Bil, file=paste(dataForRdir, "masterDF_DS_Bil.rdat", sep=""))
save(masterDF_NS_Bil, file=paste(dataForRdir, "masterDF_NS_Bil.rdat", sep=""))

#By modality
masterDF_DS_Bil_Tone <- filter(masterDF_DS_Bil, masterDF_DS_Bil$modality==98)
masterDF_DS_Bil_Light <- filter(masterDF_DS_Bil, masterDF_DS_Bil$modality==99)

masterDF_NS_Bil_Tone <- filter(masterDF_NS_Bil, masterDF_NS_Bil$modality==98)
masterDF_NS_Bil_Light <- filter(masterDF_NS_Bil, masterDF_NS_Bil$modality==99)


#Plot firing rate after the cue from CP in 5 trial bins
plotFRandCP(experiment="Bilateral AP5", masterDF=list(masterDF_DS_Bil), 
            cue="S+", graphFolder=MixedGraphFolder, trialBinSize=5, WdwStart=100, 
            WdwEnd=400, dataProcess="raw", correctOnly=FALSE, colindx="darkred", 
            capped=T, capValue = c(0, 200), yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 10, 
            neudata=allNeuronsDS)


plotFRandCP(experiment="Bilateral AP5", masterDF=list(masterDF_DS_Bil), 
            cue="S+", graphFolder=MixedGraphFolder, trialBinSize=5, WdwStart=100, 
            WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, colindx="darkred", 
            capped=T, capValue = c(0, 200), yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 10, 
            neudata=allNeuronsDS)


PercCueExc_BILAP5 <- plotFRandCP(experiment="Bilateral AP5", masterDF=list(masterDF_DS_Bil), 
            cue="S+", graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=100, 
            WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, colindx="darkred", 
            capped=T, capValue = c(0, 210), yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 10, 
            neudata=allNeuronsDS)

save(PercCueExc_BILAP5, file=paste(dataForRdir, "PercCueExc_BILAP5.rdat", sep=""))


### COMPARISON OF PROPORTIONS OF CUE EXCITED UNITS
contingency_table_comp1<- PercCueExc_BILAP5[[1]][is.even(PercCueExc_BILAP5[[1]]$bin)==FALSE, c(2, 3, 4)]
contingency_table_comp2 <-  PercCueExc_BILAP5[[1]][is.even(PercCueExc_BILAP5[[1]]$bin)==TRUE, c(2, 3, 4)]

inhCont_table_comp1 <- PercCueExc_BILAP5[[1]][is.even(PercCueExc_BILAP5[[1]]$bin)==FALSE, c(2, 5, 6)]
inhCont_table_comp2 <- PercCueExc_BILAP5[[1]][is.even(PercCueExc_BILAP5[[1]]$bin)==TRUE, c(2, 5, 6)]

contingency_table <- PercCueExc_BILAP5[[1]][, c(2, 3, 4)]
inhCont_table <- PercCueExc_BILAP5[[1]][,c(2, 5, 6)]

#Chi-sq analysis excitations:
#Bins 1, 3 and 5
critData <- contingency_table_comp1[, -1]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 2.4955, df = 2, p-value = 0.2872


#Bins 2, 4 and 6
critData <- contingency_table_comp2[,-1 ]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 4.427, df = 2, p-value = 0.1093

#Chi-sq analysis excitations:
#All bins
critData <- contingency_table[, -1]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 7.0984, df = 5, p-value = 0.2134


#Fisher exact test for count data:
#All bins

#Excitations
critData <- contingency_table[, -1]
fisher.test(critData)
# data:  critData
# p-value = 0.2718
# alternative hypothesis: two.sided


#Excitations
critData <- inhCont_table[ , -1]
fisher.test(critData)
# data:  critData
# p-value = 0.9478
# alternative hypothesis: two.sided




#Post-hoc comparisons for % cue-excited units
# Is the distribution of exc/not exc neurons different from the first (-120 to -80) and the third (-40 to 0) bin. 
contTableExc_VEH <- PercCueExc_BILAP5[[1]][ , c(3, 4)]

#Bin 1 vs. 3
fisher.test(contTableExc_VEH[c(1, 3), ], alternative="less") #95% CI: 0.000000 0.1995 Odds ratio=0.5222663 p=0.1995
#Bin 3 vs. 5
fisher.test(contTableExc_VEH[c(3, 5), ], alternative="less") #95% CI:0.0000000 99.96224 odds ratio=4.157763 p= 0.9734

#Bin 2 vs. 4
fisher.test((contTableExc_VEH[c(2, 4), ]), alternative="less") #95% CI: 0.0000000 1.353198 odds ratio= 0.4018825  p= 0.1206
#Bin 4 vs. 6
fisher.test((contTableExc_VEH[c(4, 6), ]), alternative="less") #95% CI:  0.000000 Inf odds ratio=Inf p= 1

#Correct p value all comparisons
p.adjust(p=c(0.1995, 0.9734, 0.1206, 1), method="holm")
# 0.5985 1.0000 0.4824 1.0000



#Chi-sq for inhibitions:
#Bins 1, 3 and 5
critData <- inhCont_table_comp1[,-1 ]
chisq.test(x=critData) #X-squared = 0.9071, df = 2, p-value = 0.6354

#Bins 2, 4 and 6
critData <- inhCont_table_comp2[,-1 ]
chisq.test(x=critData) #X-squared = 0.39636, df = 2, p-value = 0.8202

#Chi-sq analysis inhibitions:
#All bins
critData <- inhCont_table[, -1]
chisq.test(x=critData)
# Pearson's Chi-squared test
# 
# data:  critData
# X-squared = 1.3041, df = 5, p-value = 0.9345

plotFRandCP(experiment="Bilateral AP5 DS Tone vs light", masterDF=list(masterDF_DS_Bil_Tone, masterDF_DS_Bil_Light), 
            cue=c("S+", "S+"), graphFolder=MixedGraphFolder, trialBinSize=5, WdwStart=100, legLabels = c("Tone S+", "Light S+"),
            WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, colindx=colindx, 
            capped=T, capValue = c(0, 180), yAxMinZ = -1, yAxMaxZ = 10, yAxMaxRaw = 10, 
            neudata=allNeuronsDS)

plotFRandCP(experiment="Bilateral AP5 NS Tone vs light", masterDF=list(masterDF_NS_Bil_Light, masterDF_NS_Bil_Tone), 
            cue=c("S-", "S-"), graphFolder=MixedGraphFolder, trialBinSize=5, WdwStart=100, 
            WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, colindx=colindx, legLabels=c("Tone S-", "Light S-"),
            capped=T, capValue = c(0, 180), yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, 
            neudata=allNeuronsNS)


##### BOXPLOTS AROUND CHANGEPOIN IN BINS OF TRIALS
BIL_35bin <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 5 Bil AP5", 
                                masterDF=list(masterDF_DS_Bil, masterDF_NS_Bil), 
                                graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="Zscores", 
                                correctOnly=FALSE, points=F, lines=F, color=c(colindx[2], "darkred"), legLabels=c("S+", "S-"), 
                                yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, capValue=c(0, 210), cueExcOnly=F, 
                                neudata=allNeuronsDS)



BIL_35bin_CueExc <- plotFRBoxPlotandCP(cue=c("S+", "S-"), experiment="Exp 5 bil AP5 Cue Exc Only", 
                                       masterDF=list(masterDF_DS_Bil, masterDF_NS_Bil), 
                                       graphFolder=MixedGraphFolder, trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="Zscores", 
                                       correctOnly=FALSE, points=F, lines=F, color=c(colindx[2], "darkred"), legLabels=c("S+", "S-"), 
                                       yAxMinZ=-1, yAxMaxZ=8, yAxMaxRaw=7, capped=T, capValue=c(0, 210), cueExcOnly=T, 
                                       neudata=allNeuronsDS)


uniqueBin <- unique(BIL_35bin[[1]]$bin)

Diff_DSNS_BIL_35bin <- do.call("rbind", lapply(seq(1, length(uniqueBin)), function(i){
        DSbinSel <- BIL_35bin[[1]][BIL_35bin[[1]]$bin==uniqueBin[i], ]
        NSbinSel <- BIL_35bin[[2]][BIL_35bin[[2]]$bin==uniqueBin[i], ]
        
        matchUnitIdx <- match(DSbinSel$unitIdx, NSbinSel$unitIdx)
        NSfiring <- NSbinSel$byUnitFR[matchUnitIdx]
        
        DSbinSel$byUnitFR_NS <- NSfiring
        DSbinSel$DiffFR <- DSbinSel$byUnitFR-NSfiring
        
        DSbinSel
})
)


#This one compares S+ firing rate among bins only (and make a boxplot of it)
compareBins(data=Diff_DSNS_BIL_35bin, cueExcOnly=FALSE, color=colindx[1], ymin=-2, ymax=8, 
            graphFolder=MixedGraphFolder, experiment="Exp 5 Bil AP5", points=F, 
            comparisons=list(c(1, 3), c(3, 5)), 
            cue="S+")
# [[1]]
# idx     bin         n                      comp   W      p.val     p.adj sig
# W     1 1 vs. 3 61 vs. 35     0 to 35 vs. 70 to 105 896 0.09651408 0.3860563    
# W2    2 2 vs. 4 57 vs. 23   35 to 70 vs. 105 to 140 565 0.16935184 0.5080555    
# W1    3 3 vs. 5 35 vs. 18  70 to 105 vs. 140 to 175 339 0.67577103 1.0000000    
# W11   4 4 vs. 6 23 vs. 10 105 to 140 vs. 175 to 210 136 0.91516978 1.0000000    
# 
# [[2]]
# bins         n               comparison   W      p.val     p.adj sig
# W  1 vs. 3 61 vs. 35    0 to 35 vs. 70 to 105 896 0.09651408 0.1930282    
# W1 3 vs. 5 35 vs. 18 70 to 105 vs. 140 to 175 339 0.67577103 0.6757710    


PostCueFR_DSvsNS_BILAP5 <- compareDSvsNSfromCP(masterDF=list(masterDF_DS_Bil, masterDF_NS_Bil), 
                                                      trialBinSize=35, event="cue", correctOnly=F, paired=T, 
                                                      WdwStart=100, WdwEnd=400, capped=T, capValue=c(0, 210), 
                                                      dataProcess="Zscores")

PostCueFR_DSvsNS_BILAP5$p.adj <- p.adjust(PostCueFR_DSvsNS_BILAP5$p, method="holm") 
PostCueFR_DSvsNS_BILAP5$sig <- giveStars(PostCueFR_DSvsNS_BILAP5$p.adj)

# bin   V           p  n      p.adj sig
# V     0 to 35 483 0.000452605 61 0.00271563  **
# V1   35 to 70 768 0.322463133 57 0.95800781    
# V2  70 to 105 250 0.147355526 35 0.73677763    
# V3 105 to 140  72 0.186798096 19 0.74719238    
# V4 140 to 175  51 0.319335938 15 0.95800781    
# V5 175 to 210  14 0.320312500  8 0.95800781    


### PLOT FR AROUND THE CUE (PTSH) as a function of distance to change point
plotFRandCPhistogram(experiment="Exp 5", masterDF=list(masterDF_DS_Bil, masterDF_NS_Bil), graphFolder=MixedGraphFolder, 
                     trialBinSize=35, dataProcess="Zscores", correctOnly=FALSE, color=c("darkred", "gray30"), 
                     capped=T, capValue = c(0, 210), yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 3, 
                     psthmin=0.5, psthmax=2, imgFormat="pdf", neudata=allNeuronsDS)



### PLOT FR AROUND THE CUE (PSTH) as a function of distance to SESSION in which the CP took place (by session)
plotPSTHfromSessCP(experiment="Exp 5 S+ Tone vs Light", masterDF=list(masterDF_DS_Bil_Tone, masterDF_DS_Bil_Light), 
                   graphFolder=MixedGraphFolder, dataProcess="Zscores", comp=c("Tone S+", "Light S+"),
                   correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(0, 5), 
                   yAxMinZ = -1, yAxMaxZ = 4, yAxMaxRaw = 10, psthmin=-0.5, psthmax=2, 
                   imgFormat="pdf", neudata=allNeuronsDS)


plotPSTHfromSessCP(experiment="Exp 5 S- Tone vs Light", masterDF=list(masterDF_NS_Bil_Light, masterDF_NS_Bil_Tone), 
                   graphFolder=MixedGraphFolder, dataProcess="Zscores", comp=c("Tone S-", "Light S-"),
                   correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(0, 5), 
                   yAxMinZ = -1, yAxMaxZ = 4, yAxMaxRaw = 10, psthmin=-0.5, psthmax=2, 
                   imgFormat="pdf", neudata=allNeuronsDS)





### PLOT BOXPLOTS AROUND CUE BY SESS FROM CP AND MODALITY
# I only have 1 unit on session 1 in the Light S+ modality, so these will generate the boxplots but the statistical analyses won't work
plotBoxplotfromSessCP(experiment="Exp 5 S+ Tone vs Light", masterDF=list(masterDF_DS_Bil_Tone, masterDF_DS_Bil_Light), 
                      comp=c("Tone S+", "Light S+"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(0, 5), 
                      yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)

plotBoxplotfromSessCP(experiment="Exp 5 S+ Tone vs Light", masterDF=list(masterDF_DS_Bil_Tone, masterDF_DS_Bil_Light), 
                      comp=c("Tone S+", "Light S+"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(0, 5), 
                      yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, WdwStart=100, WdwEnd=400, 
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)



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


