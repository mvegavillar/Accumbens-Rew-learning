
setwd("E:/Dropbox/NMDA/")

######################################################################
### EXPERIMENT 1 TRAINING: NO INFUSIONS, SWITCH SHORT TO LONG ITI  ###
######################################################################

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
#funcdirect <- paste(getwd(), "/EXP4_Unilateral AP5/R functions/", sep="")
datafolder <- paste(getwd(), "/EXP1_Performance/Training/MedPC files/", sep="")
dataForRdir <- paste(getwd(), "/EXP1_Performance/Training/Data for R/", sep="")
dataForRCumulative <- paste(getwd(), "/EXP1_Performance/Training/Data for R cumulative/", sep="")
behGraphFolder <- paste(getwd(), "/EXP1_Performance/Training/Graphs/Behavior/", sep="")
CPfuncFolder <- paste(funcdirect, 'Change_Point-master/', sep="")
CPGraphFolder <- paste(getwd(), "/EXP1_Performance/Training/Graphs/Behavior/Change point/", sep="")
MixedGraphFolder <- paste(getwd(), "/EXP1_Performance/Training/Graphs/Mixed/", sep="")
NEXfiles <- paste(getwd(), "/EXP1_Performance/Training/NEX files/", sep="")
IndcumPerfGraphFolder <- paste(getwd(), "/EXP1_Performance/Training/Graphs/Behavior/Cumulative perf/", sep="")
PerfRelToCPFolder <- paste(getwd(), "/EXP1_Performance/Training/Graphs/Behavior/Beh rel to CP/", sep="")
BySessFolder <- paste(MixedGraphFolder, "By Sess from CP/", sep="")


# ### Load necessary functions. All of the files in that folder are functions I want and then there is a folder with CP functions. Load all the functions but that
# functionNames <- paste(funcdirect, list.files(funcdirect), sep="")
# functionNames <- functionNames[-grep("Change_Point-master", functionNames)]
# for(i in 1:length(functionNames)){source(functionNames[[i]])}


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

#     rat  CP CPsess   slopePre slopePost FirstCPOver80pc DynamicInterv nSess NS.PreCP
# 1  MV41 120      4 -0.2282295  2.877726             120             0     6      123
# 2  MV42 156      5  0.1423823  4.638940             156             0     6      151
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


CPdataExp1 <- CPdata
alldataExp1 <- alldata
csacqidxExp1 <- csacqidx
ratsExp1 <- rats
idxExp1 <- idx

save(csacqidxExp1, file=paste(dataForRdir, "csacqidxExp1.ridx", sep=""))
save(alldataExp1, file=paste(dataForRdir, "alldataExp1.rdat", sep=""))
save(ratsExp1, file=paste(dataForRdir, "ratsExp1.rdat", sep=""))
save(idxExp1, file=paste(dataForRdir, "idxExp1.rdat", sep=""))


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
CPframeExp1 <- rbind(DSspecificity, RealCPbyIndex[sel,])

save(CPframeExp1, file=paste(dataForRdir, "CPframeExp1.rdat", sep="")) #CP values for different behavioral indexes
save(CPdataExp1, file=paste(dataForRdir, "CPdataExp1.rdat", sep="")) #CP values and other details for PERFORMANCE INDEX (S+ specificity)


#CUMULATIVE INDIVIDUAL PERFORMANCE
cumulativeIndGraphs(numSess=7, numCues=35, sessLines=F, smartRat="NA", limit=210, dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=IndcumPerfGraphFolder, imgFormat="pdf")

#MEAN PERFORMANCE AS A FUNCTION OF DISTANCE FROM CP
PerformanceFromCP(relTrialMin=-100, relTrialMax=100, trialBinSize=5, typegraph="both", imgFormat="pdf", 
                  dataForRCumulative=dataForRCumulative, dataForRdir=dataForRdir, graphFolder=PerfRelToCPFolder, CPdata=CPdata, 
                  idx=idx, rats = rats, alldata = alldata)


#PRE CP VS. POST CP PERFORMANCE

# S+ SPECIFICITY
plot.new()
plot.window(xlim=c(0, 0.5), ylim=c(-1.5, 4))

rect(xleft=0, xright=0.25, ybottom=0, ytop=mean(CPdata$slopePre, na.rm=T), border = "white", col="gray60")
rect(xleft=0.25, xright=0.5, ybottom=0, ytop=mean(CPdata$slopePost, na.rm=T), border = "white", col="gray60")

sapply(seq(1, nrow(CPdata)), function(x){lines(x=c(0.125, 0.375), y=c(CPdata$slopePre[x], CPdata$slopePost[x]), lwd=2, col="gray20")})

abline(h=0, lty=3)
axis(side=1, at=c(0.125, 0.375), labels=c("Before CP", "After CP"), cex.axis=1.4, tick = FALSE)
axis(side=2, at=seq(-1, 4, 1), las=2, cex.axis=1.4)
mtext(side=2, line=2.5, text="S+ specificity", cex=1.5, font=2)

t.test(x=CPdata$slopePre, y=CPdata$slopePost, paired=T) #t = -13.05, df = 5, p-value = 4.713e-05. No correction

PrevsPostFolder <- paste(PerfRelToCPFolder, "Average Pre vs Post/", sep="")
PrePostCP_Perf(data=DSrespAll, CP=CPdata$CP, y_axis_label="S+ response ratio", graphFolder=PrevsPostFolder, plot=T) #V=0, P=0.000241
PrePostCP_Perf(data=NSrespAll, CP=CPdata$CP, y_axis_label="S- response ratio", graphFolder=PrevsPostFolder, plot=T) #t=1.4904, df=12, p=0.1619
PrePostCP_Perf(data=DSlatency, CP=CPdata$CP, y_axis_label="S+ latency", graphFolder=PrevsPostFolder, plot=T) #t = 7.4541, df=12, p=0.000007693
PrePostCP_Perf(data=NSlatency, CP=CPdata$CP, y_axis_label="S- latency", graphFolder=PrevsPostFolder, plot=T) #t=-2.6775, df=12, p=0.02013
PrePostCP_Perf(data=DStaskAcc, CP=CPdata$CP, y_axis_label="S+ specificity", graphFolder=PrevsPostFolder, plot=T) #t = -7.9401, df=12, p=0.00000406
PrePostCP_Perf(data=NStaskAcc, CP=CPdata$CP, y_axis_label="S- specificity", graphFolder=PrevsPostFolder, plot=T) #t = 0.85764, df = 12, p-value = 0.4079
PrePostCP_Perf(data=ITIlatency, CP=CPdata$CP, y_axis_label="ITI latency", graphFolder=PrevsPostFolder, plot=T) #t = -1.8576, df = 12, p-value = 0.08793

p.adjust(c(0.000241, 0.000007693, 0.00000406, 0.08793), method="holm")
#4.8200e-04 2.3079e-05 1.6240e-05 8.7930e-02

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
files <- paste(dataForRdir, list.files(dataForRdir), sep=""); for(i in 1:length(files)){load(files[[i]])}

#BUILD A DATAFRAME WITH FIRING RATE OF EACH NEURON ON EACH TRIAL WITH BEH INFO ABOUT EACH TRIAL
masterDF_DS <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS, CPvector=CPdata$CP, 
                                  sessionCPperRat = sessionCPperRat, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                  BLduration=2, BLneudata=allNeuronsDS)

masterDF_DScorrect <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSresponded, CPvector=CPdata$CP, 
                                         sessionCPperRat = sessionCPperRat, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                         BLduration=2, BLneudata=allNeuronsDS)

masterDF_DSmissed <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSmissed, CPvector=CPdata$CP, 
                                        sessionCPperRat = sessionCPperRat, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                        BLduration=2, BLneudata=allNeuronsDS)

masterDF_NS <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsNS, CPvector=CPdata$CP, 
                                  sessionCPperRat = sessionCPperRat, funcdirect=funcdirect, dataForRdir=dataForRdir, 
                                  BLduration=2, BLneudata=allNeuronsNS)

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

#By modality
masterDF_DS_Tone <- filter(masterDF_DS, masterDF_DS$modality==98)
masterDF_DS_Light <- filter(masterDF_DS, masterDF_DS$modality==99)



save(masterDF_DS, file=paste(dataForRdir, "masterDF_DS.rdat", sep=""))
save(masterDF_DScorrect, file=paste(dataForRdir, "masterDF_DScorrect.rdat", sep=""))
save(masterDF_DSmissed, file=paste(dataForRdir, "masterDF_DSmissed.rdat", sep=""))
save(masterDF_NS, file=paste(dataForRdir, "masterDF_NS.rdat", sep=""))
save(masterDF_DSEntry, file=paste(dataForRdir, "masterDF_DSEntry.rdat", sep=""))
save(masterDF_NSEntry, file=paste(dataForRdir, "masterDF_NSEntry.rdat", sep=""))

#masterDFsumm <- masterDFsummary(masterDF=masterDF_DS, neudata=allNeuronsDS, format="Zscores", postCueWdw=c(0, 0.4), postSpikeWdw=c(0.4, 5), orEntry=T, cuelength=10)


### PLOT FR POST-CUE (400MS WINDOW) as a function of distance to change point


#CUE
plotFRandCP(experiment="Exp 1", masterDF=list(masterDF_DS), graphFolder=MixedGraphFolder, 
            trialBinSize=5, WdwStart=100, WdwEnd=300, dataProcess="raw", correctOnly=FALSE, 
            colindx =colindx[1], capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 5, 
            yAxMaxRaw = 10, neudata=allNeuronsDS)

plotFRandCP(experiment="Exp 1", masterDF=list(masterDF_DS), graphFolder=MixedGraphFolder, 
            trialBinSize=5, WdwStart=100, WdwEnd=300, dataProcess="Zscores", correctOnly=FALSE, 
            colindx=colindx[1], capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 5, 
            yAxMaxRaw = 10, neudata=allNeuronsDS)

plotFRandCP(experiment="Exp 1", masterDF=list(masterDF_DS), graphFolder=MixedGraphFolder, 
            trialBinSize=5, WdwStart=100, WdwEnd=300, dataProcess="raw", correctOnly=TRUE, 
            colindx=colindx[1], capped=T, capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 5, 
            yAxMaxRaw = 10, neudata=allNeuronsDS)

plotFRandCP(experiment="Exp 1", masterDF=list(masterDF_DS), graphFolder=MixedGraphFolder, trialBinSize=5, 
            WdwStart=100, WdwEnd=300, dataProcess="Zscores", correctOnly=TRUE, colindx=colindx[1], capped=T, 
            capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 4, yAxMaxRaw = 10, neudata=allNeuronsDS)

plotFRandCP(cue=c("S+", "S+"), experiment="Exp 1 correct vs missed", masterDF=list(masterDF_DScorrect, masterDF_DSmissed), 
            legLabels=c("Correct S+ trials", "Missed S+ trials"), graphFolder=MixedGraphFolder, trialBinSize=5, 
            WdwStart=100, WdwEnd=300, dataProcess="Zscores", correctOnly=FALSE, colindx=colindx, capped=T, 
            capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 10, 
            neudata=c(allNeuronsDSresponded, allNeuronsDSmissed))

plotFRandCP(cue=c("S+", "S-"), experiment="Exp 1 S+ vs S-", masterDF=list(masterDF_DS, masterDF_NS), 
            legLabels=c("S+", "S-"), graphFolder=MixedGraphFolder, trialBinSize=5, 
            WdwStart=100, WdwEnd=300, dataProcess="Zscores", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-100, 100), yAxMinZ = -1, yAxMaxZ = 6, yAxMaxRaw = 10, 
            neudata=c(allNeuronsDS, allNeuronsNS), cueExcOnly = F)


plotFRandCP(cue=c("S+"), experiment="Exp 1 S+ PERCEXC", masterDF=list(masterDF_DS), 
            legLabels=c("S+"), graphFolder=MixedGraphFolder, trialBinSize=15, 
            WdwStart=100, WdwEnd=300, dataProcess="PercCueExc", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-90, 90), yAxMinZ = -1, yAxMaxZ = 6, yAxMaxRaw = 10, 
            neudata=c(allNeuronsDS, allNeuronsNS), cueExcOnly = F)


### BY MODALITY
plotFRandCP(cue=c("S+", "S+"), experiment="Exp 1 Tone vs Light S+", 
            masterDF=list(masterDF_DS_Tone, masterDF_DS_Light), 
            legLabels=c("Tone S+", "Light S+"), graphFolder=MixedGraphFolder, trialBinSize=5, 
            WdwStart=100, WdwEnd=300, dataProcess="Zscores", correctOnly=FALSE, colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-90, 90), yAxMinZ = -5, yAxMaxZ = 6, yAxMaxRaw = 10, 
            neudata=allNeuronsDS, cueExcOnly = F)

plotFRandCP(cue=c("S+"), experiment="Exp 1 S+ PERCEXC Tone S+", 
            masterDF=list(masterDF_DS_Tone), 
            legLabels=c("Tone S+"), graphFolder=MixedGraphFolder, trialBinSize=15, 
            WdwStart=100, WdwEnd=300, dataProcess="PercCueExc", correctOnly=FALSE, 
            colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-90, 90), yAxMinZ = -5, yAxMaxZ = 6, yAxMaxRaw = 10, 
            neudata=allNeuronsDS, cueExcOnly = F)

plotFRandCP(cue=c("S+"), experiment="Exp 1 S+ PERCEXC Light S+", 
            masterDF=list(masterDF_DS_Light), 
            legLabels=c("Tone S+"), graphFolder=MixedGraphFolder, trialBinSize=15, 
            WdwStart=100, WdwEnd=300, dataProcess="PercCueExc", correctOnly=FALSE, 
            colindx=c(colindx[1], "darkblue"), capped=T, 
            capValue = c(-90, 90), yAxMinZ = -5, yAxMaxZ = 6, yAxMaxRaw = 10, 
            neudata=allNeuronsDS, cueExcOnly = F)


#ENTRY
plotFRandCP_Entry(experiment="Exp 1", cue="S+", wdwLabel="Post S+ Entry", masterDF=list(masterDF_DSEntry), 
                  graphFolder=MixedGraphFolder, 
                  trialBinSize=30, WdwStart =0, WdwEnd=1000, dataProcess="Zscores", colindx =colindx[1], 
                  capped=T, capValue = c(-90, 90), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryDS)

plotFRandCP_Entry(experiment="Exp 1", cue="S-", wdwLabel="Post S- Entry", masterDF=list(masterDF_NSEntry), graphFolder=MixedGraphFolder, 
                  trialBinSize=30, WdwStart =0, WdwEnd=1000, dataProcess="Zscores", colindx =colindx[1], 
                  capped=T, capValue = c(-90, 90), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, neudata=allNeuronsEntryNS)


plotFRandCP_Entry(experiment="Exp 1", cue=c("S+", "S-"), wdwLabel="Post S+ and S- Entry", masterDF=list(masterDF_DSEntry, masterDF_NSEntry), graphFolder=MixedGraphFolder, 
                  trialBinSize=30, WdwStart =0, WdwEnd=1000, dataProcess="Zscores", colindx =c(colindx[1], "darkblue"), 
                  capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 3, legLabels=c("S+", "S-"),
                  yAxMaxRaw = 10, neudata=allNeuronsEntryDS)

plotFRandCP_Entry(experiment="Exp 1", cue=c("S+", "S-"), wdwLabel="Pre S+ and S- Entry", masterDF=list(masterDF_DSEntry, masterDF_NSEntry), graphFolder=MixedGraphFolder, 
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
PostCueFR_DSvsNS_fromCP <- compareDSvsNSfromCP(masterDF=list(masterDF_DS, masterDF_NS), trialBinSize=15, 
                                               event="cue", correctOnly=F, WdwStart=100, WdwEnd=300, capped=T, 
                                               capValue=c(-90, 90), dataProcess="Zscores")

PostCueFR_DSvsNS_fromCP$p.adj <- p.adjust(PostCueFR_DSvsNS_fromCP$p, method="holm") 
PostCueFR_DSvsNS_fromCP$sig <- giveStars(PostCueFR_DSvsNS_fromCP$p.adj)

save(PostCueFR_DSvsNS_fromCP, file=paste(dataForRdir, "PostCueFR_DSvsNS_fromCP.rdat"))

#Pre entry
PreEntryFR_DSvsNS_fromCP <- compareDSvsNSfromCP(masterDF=list(masterDF_DSEntry, masterDF_NSEntry), trialBinSize=30, 
                                                event="entry", correctOnly=F, WdwStart=-2000, WdwEnd=0, capped=T, 
                                                capValue=c(-90, 90), dataProcess="Zscores")

PreEntryFR_DSvsNS_fromCP$p.adj <- p.adjust(PreEntryFR_DSvsNS_fromCP$p, method="holm") 
PreEntryFR_DSvsNS_fromCP$sig <- giveStars(PreEntryFR_DSvsNS_fromCP$p.adj)
save(PreEntryFR_DSvsNS_fromCP, file=paste(dataForRdir, "PreEntryFR_DSvsNS_fromCP.rdat"))


#Post entry
PostEntryFR_DSvsNS_fromCP <- compareDSvsNSfromCP(masterDF=list(masterDF_DSEntry, masterDF_NSEntry), trialBinSize=30, 
                                                 event="entry", correctOnly=F, WdwStart=0, WdwEnd=1000, capped=T, 
                                                 capValue=c(-90, 90), dataProcess="Zscores")
PostEntryFR_DSvsNS_fromCP$p.adj <- p.adjust(PostEntryFR_DSvsNS_fromCP$p, method="holm") 
PostEntryFR_DSvsNS_fromCP$sig <- giveStars(PostEntryFR_DSvsNS_fromCP$p.adj)
save(PostEntryFR_DSvsNS_fromCP, file=paste(dataForRdir, "PostEntryFR_DSvsNS_fromCP.rdat"))



####################################
#Same thing but by session from CP instead of trial
plotPSTHfromSessCP(experiment="Exp 1", masterDF=list(masterDF_DS_Tone, masterDF_DS_Light), 
                   graphFolder=MixedGraphFolder, dataProcess="Zscores", comp=c("Tone S+", "Light S+"),
                   correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-2, 2), 
                   yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 10, psthmin=-0.5, psthmax=2, 
                   imgFormat="pdf", neudata=allNeuronsDS)


### PLOT BOXPLOTS AROUND CUE BY SESS FROM CP AND DRUG GROUP
plotBoxplotfromSessCP(experiment="Exp 1", masterDF=list(masterDF_DS_Tone, masterDF_DS_Light), 
                      comp=c("Tone S+", "Light S+"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                      correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-2, 2), 
                      yAxMinZ = -2, yAxMaxZ = 12, yAxMaxRaw = 10, WdwStart=100, WdwEnd=300, 
                      removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)


### THIS ONE IS THE SAME BUT GROUPING THE SESSIONS BEFORE AND AFTER THE CP
plotBoxplotPrePostCP(experiment="Exp 1", masterDF=list(masterDF_DS_Tone, masterDF_DS_Light), 
                     comp=c("Tone S+", "Light S+"), graphFolder=MixedGraphFolder, dataProcess="Zscores", 
                     correctOnly=FALSE, cueExcOnly=FALSE, color=colindx, sessFromCP = c(-5, 5), 
                     yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, WdwStart=100, WdwEnd=300, 
                     removeOutliers=F, imgFormat="pdf", neudata=allNeuronsDS, morethanIQR=T)





### PLOT FR IN THE POSTCUE WINDOW AND TAIL PERIOD BY RAT BY TRIAL BY NEURON
#PLOT_perRatTrial_FRandBEH(masterDFsumm=masterDFsumm, wdw="PostCue", threshold=2, MixedGraphFolder=MixedGraphFolder)
#PLOT_perRatTrial_FRandBEH(masterDFsumm=masterDFsumm, wdw="PostSpike", threshold=2, MixedGraphFolder=MixedGraphFolder)


### PLOT FR AROUND THE CUE (PTSH) as a function of distance to change point
plotFRandCPhistogram(experiment="Exp 1", masterDF=list(masterDF_DS, masterDF_NS), graphFolder=MixedGraphFolder, 
                     trialBinSize=15, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
                     capped=T, capValue = c(-90, 90), yAxMinZ = -1, yAxMaxZ = 3, yAxMaxRaw = 3, 
                     psthmin=0.5, psthmax=2, imgFormat="pdf", neudata=allNeuronsDS)


### PLOT FR AROUND ENTRIES (PTSH) as a function of distance to change point
plotFRandCPhistogram_Entry(experiment="Exp 1", masterDF=list(masterDF_DSEntry, masterDF_NSEntry), graphFolder=MixedGraphFolder, 
                           trialBinSize=30, dataProcess="Zscores", correctOnly=FALSE, color=c(colindx[1], "darkblue"), 
                           capped=T, capValue = c(-120, 120), yAxMinZ = -1, yAxMaxZ = 2, yAxMaxRaw = 3, 
                           WdwStart=-2000, WdwEnd=2000, imgFormat="pdf", neudata=allNeuronsEntryDS)



# #COMPARE WITH OTHER GROUPS
# compareCPs(data=list(CPdataExp4, CPdataExp10tr, CPdataHYBunilAP5, CPdataHYBbilAP5) , imgFormat="png", expNames=c("Task1", "Task2", "UnilAP5", "BilAP5"), colindx=c("gray30", "gray45", "gray60", "gray75"), graphFolder=behGraphFolder, minSess=5, graphs=c("Change point", "Non learners", "Dynamic Interval", "Asymptote"))


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




## SCATTERPLOT WHERE EACH SESSION IS A DOT AND THE AXES ARE --> X axis: performance index; Y axis: mean FR. 
#Load the version of toplot that I need for this scatterplot. "toplot" is made and save by the function "megaplot". So run the megaplot line that will give me the desired "toplot"
load(file=paste(dataForRdir, "toplot.rdat", sep=""))

toplotSumm <- toplot %>%
        group_by(rat, CPfromSess) %>%
        summarise(BeforeCP=unique(BeforeCP),
                  CSplusTA=unique(CSplusTA),
                  CSminusTA=unique(CSminusTA),
                  CSplusLat=unique(CSplusLat),
                  CSminusLat=unique(CSminusLat),
                  DSRR=unique(DSRR),
                  NSRR=unique(NSRR),
                  ZDS=mean(as.numeric(ZDS, na.rm=T)),
                  ZNS=mean(as.numeric(ZNS, na.rm=T))
        )

toplotSumm <- as.data.frame(toplotSumm)
CPsessIndex <- toplotSumm$BeforeCP==0

toplotSummNoCPsess <- toplotSumm[CPsessIndex==F, ]

plot.new()
plot.window(xlim=c(-3, 6), ylim=c(-1, 6.5))
sapply(c(-1, 1), function(x){
        dataSel <- toplotSummNoCPsess[toplotSummNoCPsess$BeforeCP==x,]
        if(x==1){pch=19; col=colindx[1]}
        if(x==-1){pch=19; col=colindx[2]}
        perf <- as.numeric(as.character(dataSel$CSplusTA))
        FR <- as.numeric(as.character(dataSel$ZDS))
        
        points(x=perf, y=FR, pch=pch, cex=2, col=col)
        
        axis(side=1, at=seq(-3, 6, 1), cex.axis=1.4)
        axis(side=2, at=seq(-1, 6.5, 1), cex.axis=1.4, las=2)
        
        abline(h=0, lty=2)
        abline(v=0)
        
        mtext(side=2, line=2.5, text="Firing rate (Z sc.)", cex=1.5, font=2)
        mtext(side=1, line=2.5, text="Performance index", cex=1.5, font=2)
        
})

perfAll <- as.numeric(as.character(toplotSummNoCPsess$CSplusTA))
FRAll <- as.numeric(as.character(toplotSummNoCPsess$ZDS))

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


