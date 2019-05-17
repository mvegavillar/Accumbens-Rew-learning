
######################################################################################################
### EXPERIMENT 4: UNILATERAL AP5 INFUSIONS: EXTINCTION TEST (MV )                                  ###
######################################################################################################


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
datafolder_L <- paste(getwd(), "/EXP4_Unilateral AP5 Extinction/MedPC files/Learners/", sep="") #Rats that learned the task before the test
datafolder_NL <- paste(getwd(), "/EXP4_Unilateral AP5 Extinction/MedPC files/Nonlearners/", sep="") #Rats that didn't learn the task before the test
dataForRdir_L <- paste(getwd(), "/EXP4_Unilateral AP5 Extinction/Data for R/Learners/", sep="")
dataForRdir_NL <- paste(getwd(), "/EXP4_Unilateral AP5 Extinction/Data for R/Nonlearners/", sep="")
dataForRCumulative_L <- paste(getwd(), "/EXP4_Unilateral AP5 Extinction/Data for R cumulative/Learners/", sep="")
dataForRCumulative_NL <- paste(getwd(), "/EXP4_Unilateral AP5 Extinction/Data for R cumulative/Nonlearners/", sep="")
behGraphFolder <- paste(getwd(), "/EXP4_Unilateral AP5 Extinction/Graphs/Behavior/", sep="")
CPfuncFolder <- paste(getwd(), '/R functions/Change_Point-master/', sep="")
CPGraphFolder <- paste(getwd(), "/EXP4_Unilateral AP5 Extinction/Graphs/Behavior/Change point/", sep="")
MixedGraphFolder <- paste(getwd(), "/EXP4_Unilateral AP5 Extinction/Graphs/Mixed/", sep="")
NEXfiles_L <- paste(getwd(), "/EXP4_Unilateral AP5 Extinction/NEX files/Learners/", sep="")
NEXfiles_NL <- paste(getwd(), "/EXP4_Unilateral AP5 Extinction/NEX files/Nonlearners/", sep="")
IndcumPerfGraphFolder <- paste(getwd(), "/EXP4_Unilateral AP5 Extinction/Graphs/Behavior/Cumulative/", sep="")

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
load(file=paste(funcdirect, "SessionBOXPLOT.R", sep=""))
load(file=paste(funcdirect, "SessionPSTH.R", sep=""))


#############################################################################################################
########### LEARNERS: MV27, 30, 31, 51, 52, 55 ##############################################################
#############################################################################################################

### EXTRACT BASIC BEHAVIOR OBJECTS
MedPCextract(MovAvg="Impinged only", funcdirect = funcdirect, datafolder = datafolder_L, 
             dataForRdir = dataForRdir_L, dataForRCumulative=dataForRCumulative_L, cuelength=10)


# Load important behavior-related objects
files <- paste(dataForRdir_L, list.files(dataForRdir_L), sep="")
filesCum <- paste(dataForRCumulative_L, list.files(dataForRCumulative_L), sep="")
for(i in 1:length(files)){load(files[[i]])}
for(i in 1:length(filesCum)){load(filesCum[[i]])}

csacqidx$session <- rep(1, nrow(csacqidx))

#### BEHAVIORAL GRAPHS

#Individual performance
DStaskAccbyTrial <- sapply(seq(1, 35), function(x){
        sapply(seq(1, length(DStaskAcc)), function(j){
                DStaskAcc[[j]][x]
        })
})
DStaskAccMeanTrial <- colMeans(DStaskAccbyTrial, na.rm=T)
DStaskAccSEMTrial <- colSds(DStaskAccbyTrial, na.rm=T)


NStaskAccbyTrial <- sapply(seq(1, 35), function(x){
        sapply(seq(1, length(NStaskAcc)), function(j){
                NStaskAcc[[j]][x]
        })
})
NStaskAccMeanTrial <- colMeans(NStaskAccbyTrial, na.rm=T)
NStaskAccSEMTrial <- colSds(NStaskAccbyTrial, na.rm=T)

#trial by trial performance on the extinction test
plot.new()
plot.window(xlim=c(0, 35), ylim=c(0, 5.5))
lines(DStaskAccMeanTrial, col="black", lwd=2)

#5 trial bin performance on the extinction test

binCuts <- seq(1, length(DStaskAccMeanTrial), by=5)
binIndex <- 1:length(DStaskAccMeanTrial)

binAssignment <- findInterval(binIndex, binCuts)


#DS SPECIFICITY
DStaskAcc5trialBins <- sapply(seq(1, length(binCuts)), function(m){
        mean(DStaskAccMeanTrial[binAssignment==m], na.rm=T)
})

DStaskAcc5trialBinsSEM <- sapply(seq(1, length(binCuts)), function(m){
        sd(DStaskAccMeanTrial[binAssignment==m], na.rm=T)/sqrt(sum(binAssignment==m))
})


#NS SPECIFICITY
NStaskAcc5trialBins <- sapply(seq(1, length(binCuts)), function(m){
        mean(NStaskAccMeanTrial[binAssignment==m], na.rm=T)
})

NStaskAcc5trialBinsSEM <- sapply(seq(1, length(binCuts)), function(m){
        sd(NStaskAccMeanTrial[binAssignment==m], na.rm=T)/sqrt(sum(binAssignment==m))
})



plot.new()
plot.window(xlim=c(0, length(binCuts)), ylim=c(-2, 6))

points(DStaskAcc5trialBins, pch=19, cex=2)
lines(DStaskAcc5trialBins, lwd=2)
errBars(x=seq(1, length(binCuts)), y=DStaskAcc5trialBins, err=DStaskAcc5trialBinsSEM, hatLength = 0.05)

points(NStaskAcc5trialBins, pch=19, cex=2, col="gray40")
lines(NStaskAcc5trialBins, lwd=2, col="gray40")
errBars(x=seq(1, length(binCuts)), y=NStaskAcc5trialBins, err=NStaskAcc5trialBinsSEM, hatLength = 0.05)

labVals <- c("1-5", "6-10", "11-15", "16-20", "21-25", "26-30", "31-35")
axis(side=1, at=seq(1, 7, by=1), labels=labVals, cex.axis=2, las=2)
axis(side=2, at=seq(-2, 6, by=1), cex.axis=2, las=2, pos=0.5)
abline(h=0, lty=3)

legend("topright", legend=c("S+", "S-"), col=c("black", "gray40"), lwd=2)

### For analysis:
DS_PerfIndex_Longform <- do.call("rbind", lapply(seq(1, nrow(DStaskAccbyTrial)), function(x){
        ratByBin <- sapply(unique(binAssignment), function(y){
                mean(DStaskAccbyTrial[x, binAssignment==y], na.rm=T)
        })
        data.frame(cue="S+", rat=rats[x], PerfIndex=ratByBin, bin=unique(binAssignment))
})
)

NS_PerfIndex_Longform <- do.call("rbind", lapply(seq(1, nrow(NStaskAccbyTrial)), function(x){
        ratByBin <- sapply(unique(binAssignment), function(y){
                mean(NStaskAccbyTrial[x, binAssignment==y], na.rm=T)
        })
        data.frame(cue="S-", rat=rats[x], PerfIndex=ratByBin, bin=unique(binAssignment))
})
)

PerfIndex_ForANOVA <- rbind(DS_PerfIndex_Longform, NS_PerfIndex_Longform)

ezANOVA(data=PerfIndex_ForANOVA, wid=rat, within = c(cue, bin), dv=PerfIndex) #Check this for questions about Mauchly's test and GG correction. https://www.r-exercises.com/2016/11/29/repeated-measures-anova-in-r-exercises/
# $`ANOVA`
#    Effect DFn DFd          F            p p<.05       ges
# 1     cue   1   5 119.926109 0.0001104032     * 0.8477629
# 2     bin   1   5   2.381756 0.1834075990       0.1937092
# 3 cue:bin   1   5  10.526758 0.0228349341     * 0.3567954

# $`Mauchly's Test for Sphericity`: it tests the null hypothesis that variance across each level of a particular within-subjects factor is equal 
# Effect            W            p p<.05
# 3    bin 7.068002e-17 1.118896e-11     *
# Assumption of sphericity is violated for the "bin" factor, which increases the likelihood of Type II error (failing to reject a null hypothesis that is actually wrong aka less power)         
#Greenhouse-Geiser corrects for that and recalculates p values of the effects that rely on "bin" (main effect of bin and interaction cue*bin). The results are the same.

#$`Sphericity Corrections`
#    Effect       GGe      p[GG] p[GG]<.05       HFe       p[HF] p[HF]<.05
# 3     bin 0.4565120 0.18282855           1.0640495 0.117037877          
# 4 cue:bin 0.4143815 0.02049974         * 0.8564852 0.002310229         *

#With "aov" instead of "ezANOVA"
PerfIndex_ForANOVA$cue <- factor(PerfIndex_ForANOVA$cue)
PerfIndex_ForANOVA$bin <- factor(PerfIndex_ForANOVA$bin)

options(contrasts = c("contr.sum","contr.poly"))
aovtest_Perf_L <- summary(with(PerfIndex_ForANOVA, 
                aov(PerfIndex ~ cue*bin + Error(rat/(cue*bin)))
))

qqnorm(PerfIndex_ForANOVA$PerfIndex)
hist(PerfIndex_ForANOVA$PerfIndex)

save(aovtest_Perf_L, file=paste(dataForRdir_L, "aovtest_Perf_L.rdat", sep=""))


# RASTERS
load(file=paste(funcdirect, "MAKERASTER.r", sep=""))
MAKERASTER(i=1, data=alldata, idxdata=csacqidx); title(main=paste(rats[[1]], "Extinction test"))
MAKERASTER(i=2, data=alldata, idxdata=csacqidx); title(main=paste(rats[[2]], "Extinction test"))
MAKERASTER(i=3, data=alldata, idxdata=csacqidx); title(main=paste(rats[[3]], "Extinction test"))
MAKERASTER(i=4, data=alldata, idxdata=csacqidx); title(main=paste(rats[[4]], "Extinction test"))
MAKERASTER(i=6, data=alldata, idxdata=csacqidx); title(main=paste(rats[[6]], "Extinction test"))


#### NEURONAL DATA_Extract firing rate data related to cues and entries on BOTH SIDES (VEH and AP5)
allNeuronsDS_L_VEH = neuralhist (funcdirect=funcdirect, path=NEXfiles_L, event=1, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=T, side="vehicle")
allNeuronsNS_L_VEH= neuralhist (funcdirect=funcdirect, path=NEXfiles_L, event=2, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=T, side="vehicle")
allNeuronsDSresponded_L_VEH= neuralhist (funcdirect=funcdirect, path=NEXfiles_L, event=5, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsNSresponded_L_VEH= neuralhist (funcdirect=funcdirect, path=NEXfiles_L, event=7, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsDSmissed_L_VEH= neuralhist (funcdirect=funcdirect, path=NEXfiles_L, event=6, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsNSmissed_L_VEH= neuralhist (funcdirect=funcdirect, path=NEXfiles_L, event=8, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")

allNeuronsDS_L_AP5 = neuralhist (funcdirect=funcdirect, path=NEXfiles_L, event=1, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=T, side="drug")
allNeuronsNS_L_AP5= neuralhist (funcdirect=funcdirect, path=NEXfiles_L, event=2, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=T, side="drug")
allNeuronsDSresponded_L_AP5= neuralhist (funcdirect=funcdirect, path=NEXfiles_L, event=5, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsNSresponded_L_AP5= neuralhist (funcdirect=funcdirect, path=NEXfiles_L, event=7, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsDSmissed_L_AP5= neuralhist (funcdirect=funcdirect, path=NEXfiles_L, event=6, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsNSmissed_L_AP5= neuralhist (funcdirect=funcdirect, path=NEXfiles_L, event=8, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")

#sAVE THESE OBJECTS
save(allNeuronsDS_L_VEH, file=paste(dataForRdir, 'allNeuronsDS_L_VEH.rdat', sep=""))
save(allNeuronsNS_L_VEH, file=paste(dataForRdir, 'allNeuronsNS_L_VEH.rdat', sep=""))
save(allNeuronsDSresponded_L_VEH, file=paste(dataForRdir, 'allNeuronsDSresponded_L_VEH.rdat', sep=""))
save(allNeuronsNSresponded_L_VEH, file=paste(dataForRdir, 'allNeuronsNSresponded_L_VEH.rdat', sep=""))
save(allNeuronsDSmissed_L_VEH, file=paste(dataForRdir, 'allNeuronsDSmissed_L_VEH.rdat', sep=""))
save(allNeuronsNSmissed_L_VEH, file=paste(dataForRdir, 'allNeuronsNSmissed_L_VEH.rdat', sep=""))
save(allNeuronsDS_L_AP5, file=paste(dataForRdir, 'allNeuronsDS_L_AP5.rdat', sep=""))
save(allNeuronsNS_L_AP5, file=paste(dataForRdir, 'allNeuronsNS_L_AP5.rdat', sep=""))
save(allNeuronsDSresponded_L_AP5, file=paste(dataForRdir, 'allNeuronsDSresponded_L_AP5.rdat', sep=""))
save(allNeuronsNSresponded_L_AP5, file=paste(dataForRdir, 'allNeuronsNSresponded_L_AP5.rdat', sep=""))
save(allNeuronsDSmissed_L_AP5, file=paste(dataForRdir, 'allNeuronsDSmissed_L_AP5.rdat', sep=""))
save(allNeuronsNSmissed_L_AP5, file=paste(dataForRdir, 'allNeuronsNSmissed_L_AP5.rdat', sep=""))


#BUILD A DATAFRAME WITH FIRING RATE OF EACH NEURON ON EACH TRIAL WITH BEH INFO ABOUT EACH TRIAL

masterDF_DS_L_VEH <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS_L_VEH, 
                                      CPvector=rep(1, length(rats)), funcdirect=funcdirect, 
                                      dataForRdir=dataForRdir_L, 
                                      BLduration=2, sessionCPperRat=rep(1, length(rats)), 
                                      BLneudata=allNeuronsDS_L_VEH)
masterDF_DS_L_AP5 <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS_L_AP5, 
                                        CPvector=rep(1, length(rats)), funcdirect=funcdirect, 
                                        dataForRdir=dataForRdir_L, 
                                        BLduration=2, sessionCPperRat=rep(1, length(rats)), 
                                        BLneudata=allNeuronsDS_L_AP5)
masterDF_NS_L_VEH <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_L_VEH, 
                                        CPvector=rep(1, length(rats)), funcdirect=funcdirect, 
                                        dataForRdir=dataForRdir_L, 
                                        BLduration=2, sessionCPperRat=rep(1, length(rats)), 
                                        BLneudata=allNeuronsNS_L_VEH)
masterDF_NS_L_AP5 <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_L_AP5, 
                                        CPvector=rep(1, length(rats)), funcdirect=funcdirect, 
                                        dataForRdir=dataForRdir_L, 
                                        BLduration=2, sessionCPperRat=rep(1, length(rats)), 
                                        BLneudata=allNeuronsNS_L_AP5)
masterDF_DSMissed_L_VEH <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSmissed_L_VEH, 
                                        CPvector=rep(1, length(rats)), funcdirect=funcdirect, 
                                        dataForRdir=dataForRdir_L, 
                                        BLduration=2, sessionCPperRat=rep(1, length(rats)), 
                                        BLneudata=allNeuronsDSmissed_L_VEH)
masterDF_DSMissed_L_AP5 <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSmissed_L_AP5, 
                                        CPvector=rep(1, length(rats)), funcdirect=funcdirect, 
                                        dataForRdir=dataForRdir_L, 
                                        BLduration=2, sessionCPperRat=rep(1, length(rats)), 
                                        BLneudata=allNeuronsDSmissed_L_AP5)
masterDF_DSResp_L_VEH <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSresponded_L_VEH, 
                                              CPvector=rep(1, length(rats)), funcdirect=funcdirect, 
                                              dataForRdir=dataForRdir_L, 
                                              BLduration=2, sessionCPperRat=rep(1, length(rats)), 
                                              BLneudata=allNeuronsDSresponded_L_VEH)
masterDF_DSResp_L_AP5 <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDSresponded_L_AP5, 
                                              CPvector=rep(1, length(rats)), funcdirect=funcdirect, 
                                              dataForRdir=dataForRdir_L, 
                                              BLduration=2, sessionCPperRat=rep(1, length(rats)), 
                                              BLneudata=allNeuronsDSresponded_L_AP5)


save(masterDF_DS_L_VEH, file=paste(dataForRdir_L, "masterDF_DS_L_VEH.rdat", sep=""))
save(masterDF_DS_L_AP5, file=paste(dataForRdir_L, "masterDF_DS_L_AP5.rdat", sep=""))
save(masterDF_NS_L_VEH, file=paste(dataForRdir_L, "masterDF_NS_L_VEH.rdat", sep=""))
save(masterDF_NS_L_AP5, file=paste(dataForRdir_L, "masterDF_NS_L_AP5.rdat", sep=""))
save(masterDF_DSMissed_L_VEH, file=(paste(dataForRdir_L, "masterDF_DSMissed_L_VEH.rdat", sep="")))
save(masterDF_DSMissed_L_AP5, file=(paste(dataForRdir_L, "masterDF_DSMissed_L_AP5.rdat", sep="")))
save(masterDF_DSResp_L_VEH, file=(paste(dataForRdir_L, "masterDF_DSResp_L_VEH.rdat", sep="")))
save(masterDF_DSResp_L_AP5, file=(paste(dataForRdir_L, "masterDF_DSResp_L_AP5.rdat", sep="")))



### PLOT FR POST-CUE (100-400MS WINDOW) as a function of distance from 1st trial
#CUE
plotFRandCP(experiment="Exp 4 EXT Learners VEH side", cue=c("S+", "S-"), 
            masterDF=list(masterDF_DS_L_VEH, masterDF_NS_L_VEH), graphFolder=MixedGraphFolder, 
            trialBinSize=5, WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
            colindx =c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), capped=T, capValue = c(1, 35), 
            yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_L_VEH)
plotFRandCP(experiment="Exp 4 EXT Learners VEH side CUE EXC", cue=c("S+", "S-"), 
            masterDF=list(masterDF_DS_L_VEH, masterDF_NS_L_VEH), graphFolder=MixedGraphFolder, cueExcOnly = TRUE,
            trialBinSize=5, WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
            colindx =c(colindx[1], "darkblue"), legLabels=c("S+", "S-"), capped=T, capValue = c(1, 35), 
            yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_L_VEH)

plotFRandCP(experiment="Exp 4 EXT Learners VEH side PERCEXC", cue=c("S+"), 
            masterDF=list(masterDF_DS_L_VEH), graphFolder=MixedGraphFolder, cueExcOnly = FALSE,
            trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, 
            colindx =c(colindx[1], "darkblue"), legLabels=c("S+"), capped=T, capValue = c(1, 35), 
            yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_L_VEH)


plotFRandCP(experiment="Exp 4 EXT Learners AP5 side", cue=c("S+", "S-"), 
            masterDF=list(masterDF_DS_L_AP5, masterDF_NS_L_AP5), graphFolder=MixedGraphFolder, 
            trialBinSize=5, WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
            colindx =c(colindx[2], "darkred"), legLabels=c("S+", "S-"), capped=T, capValue = c(1, 35), 
            yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_L_AP5)
plotFRandCP(experiment="Exp 4 EXT Learners AP5 side CUE EXC", cue=c("S+", "S-"), 
            masterDF=list(masterDF_DS_L_AP5, masterDF_NS_L_AP5), graphFolder=MixedGraphFolder, cueExcOnly = TRUE,
            trialBinSize=5, WdwStart=100, WdwEnd=400, dataProcess="Zscores", correctOnly=FALSE, 
            colindx =c(colindx[2], "darkred"), legLabels=c("S+", "S-"), capped=T, capValue = c(1, 35), 
            yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_L_AP5)
plotFRandCP(experiment="Exp 4 EXT Learners AP5 side PercCueExc", cue=c("S+"), 
            masterDF=list(masterDF_DS_L_AP5), graphFolder=MixedGraphFolder, 
            trialBinSize=35, WdwStart=100, WdwEnd=400, dataProcess="PercCueExc", correctOnly=FALSE, 
            colindx =c(colindx[2], "darkred"), legLabels=c("S+"), capped=T, capValue = c(1, 35), 
            yAxMinZ = -1, yAxMaxZ = 8, yAxMaxRaw = 10, neudata=allNeuronsDS_L_AP5)



### PLOT FR AROUND THE CUE (PTSH) 
#All units
SessionPSTH(experiment="Exp 4 Unil AP5 Learners EXT", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
        graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
        yAxMinZ = -2, yAxMaxZ = 8, yAxMaxRaw = 10, psthmin=0.5, psthmax=2, imgFormat="pdf", 
        neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=F, 
        legendLabels=c("S+ VEH side L", "S+ AP5 side L"))

#Cue excited only
SessionPSTH(experiment="Exp 4 Unil AP5 Learners EXT Cue exc only", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
            graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
            yAxMinZ = -2, yAxMaxZ = 16, yAxMaxRaw = 10, psthmin=0.5, psthmax=2, imgFormat="pdf", capped = F,
            neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=T, 
            legendLabels=c("S+ VEH side L", "S+ AP5 side L"))


#Make data frame with info about proportion of cue-exc neurons (this I know from calculations inside the previous function, you can check the % cue exc. printed on the graph to calculate this)
ncueExcVEH <- 12
ncueExcAP5 <- 8
forChiSq <- data.frame(side=c("VEH", "AP5"), cueexc=c(ncueExcVEH, ncueExcAP5), noncueexc=c(38-ncueExcVEH, 39-ncueExcAP5))
chisq.test(forChiSq[, -1]) 
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  forChiSq[, -1]
# X-squared = 0.71784, df = 1, p-value = 0.3969

fisher.test(forChiSq[, -1]) #0.3073 



barplot(forChiSq$cueexc/(forChiSq$cueexc+forChiSq$noncueexc), col=colindx, border = NA, ylim=c(0, 1))


#First few trials:
SessionPSTH(experiment="Exp 4 Unil AP5 Learners EXT trials 1 to 10", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
            graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
            yAxMinZ = -2, yAxMaxZ = 8, yAxMaxRaw = 10, psthmin=0.5, psthmax=2, imgFormat="pdf", capped=T, capValue = c(1, 10),
            neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=F, legendLabels=c("S+ VEH side L", "S+ AP5 side L"))

SessionPSTH(experiment="Exp 4 Unil AP5 Learners EXT trials 1 to 10 Cue Exc only", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
            graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
            yAxMinZ = -2, yAxMaxZ = 15, yAxMaxRaw = 10, psthmin=0.5, psthmax=2, imgFormat="pdf", capped=T, capValue = c(1, 10),
            neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=T, legendLabels=c("S+ VEH side L", "S+ AP5 side L"))

#Middle trials:
SessionPSTH(experiment="Exp 4 Unil AP5 Learners EXT trials 11 to 20", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
            graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
            yAxMinZ = -2, yAxMaxZ = 8, yAxMaxRaw = 10, psthmin=0.5, psthmax=2, imgFormat="pdf", capped=T, capValue = c(11, 20),
            neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=F, legendLabels=c("S+ VEH side L", "S+ AP5 side L"))

SessionPSTH(experiment="Exp 4 Unil AP5 Learners EXT trials 11 to 20 Cue Exc only", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
            graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
            yAxMinZ = -2, yAxMaxZ = 8, yAxMaxRaw = 10, psthmin=0.5, psthmax=2, imgFormat="pdf", capped=T, capValue = c(11, 20),
            neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=T, legendLabels=c("S+ VEH side L", "S+ AP5 side L"))

#Last trials:
SessionPSTH(experiment="Exp 4 Unil AP5 Learners EXT trials 21 to 30", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
            graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
            yAxMinZ = -2, yAxMaxZ = 8, yAxMaxRaw = 10, psthmin=0.5, psthmax=2, imgFormat="pdf", capped=T, capValue = c(21, 30),
            neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=F, legendLabels=c("S+ VEH side L", "S+ AP5 side L"))

SessionPSTH(experiment="Exp 4 Unil AP5 Learners EXT trials 21 to 30 Cue Exc only", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
            graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
            yAxMinZ = -2, yAxMaxZ = 8, yAxMaxRaw = 10, psthmin=0.5, psthmax=2, imgFormat="pdf", capped=T, capValue = c(21, 30),
            neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=T, legendLabels=c("S+ VEH side L", "S+ AP5 side L"))



#Responded and missed S+ trials separately (all units)
SessionPSTH(experiment="Exp 4 Unil AP5 Learners EXT Responded trials", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
            graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=T, color=colindx,  
            yAxMinZ = -2, yAxMaxZ = 8, yAxMaxRaw = 10, psthmin=0.5, psthmax=2, imgFormat="pdf", 
            neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=F, legendLabels=c("S+ resp VEH side L", "S+ resp AP5 side L"))


SessionPSTH(experiment="Exp 4 Unil AP5 EXT Learners Missed trials", masterDF=list(masterDF_DSMissed_L_VEH, masterDF_DSMissed_L_AP5), 
            graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
            yAxMinZ = -2, yAxMaxZ = 8, yAxMaxRaw = 10, psthmin=0.5, psthmax=2, imgFormat="pdf", 
            neudata=list(allNeuronsDSmissed_L_VEH, allNeuronsDSmissed_L_AP5), cueExcOnly=F, legendLabels=c("S+ missed VEH side L", "S+ missed AP5 side L"))


# Boxplot for a specific window and analysis
# 100-400ms. All units. 
SessionBOXPLOT(experiment="Exp 4 Unil AP5 Learners EXT", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
            graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
            yAxMinZ = -2, yAxMaxZ = 6, yAxMaxRaw = 10, imgFormat="pdf", WdwStart = 100, WdwEnd = 400,
            neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=F, 
            legendLabels=c("S+ VEH side L", "S+ AP5 side L"))
# [[2]]
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  first and second
# W = 1034400, p-value = 2.2e-16
# alternative hypothesis: true location shift is greater than 0


# 750-2000ms. All units. 
SessionBOXPLOT(experiment="Exp 4 Unil AP5 Learners EXT TAIL", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
               graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
               yAxMinZ = -2, yAxMaxZ = 6, yAxMaxRaw = 10, imgFormat="pdf", WdwStart = 750, WdwEnd = 2000,
               neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=F, 
               legendLabels=c("S+ VEH side L", "S+ AP5 side L"))

# [[2]]
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  first and second
# W = 906790, p-value = 7.753e-07
# alternative hypothesis: true location shift is greater than 0


p.adjust(p=c(2.2e-16, 7.753e-07), method = "holm") #4.400e-16 7.753e-07



# -2000ms to 0. All units. 
SessionBOXPLOT(experiment="Exp 4 Unil AP5 Learners BASELINE", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
               graphFolder=MixedGraphFolder, dataProcess="raw", correctOnly=FALSE, color=colindx, 
               yAxMinZ = -2, yAxMaxZ = 6, yAxMaxRaw = 10, imgFormat="pdf", WdwStart = -2000, WdwEnd = 0,
               neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=F, 
               legendLabels=c("S+ VEH side L", "S+ AP5 side L"))

#[[2]]

#Wilcoxon rank sum test with continuity correction

#data:  first and second
#W = 816250, p-value = 0.5185
#alternative hypothesis: true location shift is greater than 0


# 100-400ms. Raw scoresAll units. 
SessionBOXPLOT(experiment="Exp 4 Unil AP5 Learners PEAK CUE EXC RAW", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
               graphFolder=MixedGraphFolder, dataProcess="raw", correctOnly=FALSE, color=colindx, 
               yAxMinZ = -2, yAxMaxZ = 6, yAxMaxRaw = 10, imgFormat="pdf", WdwStart = 100, WdwEnd = 400,
               neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=F, 
               legendLabels=c("S+ VEH side L", "S+ AP5 side L"))

#[[2]]

#Wilcoxon rank sum test with continuity correction

#data:  first and second
#W = 1070300, p-value = 2.841e-16
#alternative hypothesis: true location shift is not equal to 0




#### CUE-EXCITED UNITS ALONE
#100-400ms. Cue-excited units alone
SessionBOXPLOT(experiment="Exp 4 Unil AP5 Learners EXT 100 t0 400", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
               graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
               yAxMinZ = -2, yAxMaxZ = 12, yAxMaxRaw = 10, imgFormat="pdf", WdwStart = 100, WdwEnd = 400,
               neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=T, capped=T,
               legendLabels=c("S+ VEH side L", "S+ AP5 side L"))
# [[2]]
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  first and second
# W = 14223, p-value < 2.2e-16
# alternative hypothesis: true location shift is greater than 0


#750-2000ms Cue-excited units alone
SessionBOXPLOT(experiment="Exp 4 Unil AP5 Learners EXT 750 to 2000", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
               graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
               yAxMinZ = -2, yAxMaxZ = 12, yAxMaxRaw = 10, imgFormat="pdf", WdwStart = 750, WdwEnd = 2000,
               neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), cueExcOnly=T, capped=T,
               legendLabels=c("S+ VEH side L", "S+ AP5 side L"))
# [[2]]
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  first and second
# W = 11088, p-value = 0.0001746
# alternative hypothesis: true location shift is greater than 0




########### NON LEARNERS: MV28, 29 ##############################################################

### EXTRACT BASIC BEHAVIOR OBJECTS
MedPCextract(MovAvg="Impinged only", funcdirect = funcdirect, datafolder = datafolder_NL, 
             dataForRdir = dataForRdir_NL, dataForRCumulative=dataForRCumulative_NL, cuelength=10)


# Load important behavior-related objects
files <- paste(dataForRdir_NL, list.files(dataForRdir_NL), sep="")
filesCum <- paste(dataForRCumulative_NL, list.files(dataForRCumulative_NL), sep="")
for(i in 1:length(files)){load(files[[i]])}
for(i in 1:length(filesCum)){load(filesCum[[i]])}

csacqidx$session <- rep(1, nrow(csacqidx))

#### BEHAVIORAL GRAPHS

#Individual performance
DStaskAccbyTrial <- sapply(seq(1, 35), function(x){
        sapply(seq(1, length(DStaskAcc)), function(j){
                DStaskAcc[[j]][x]
        })
})
DStaskAccMeanTrial <- colMeans(DStaskAccbyTrial, na.rm=T)
DStaskAccSEMTrial <- colSds(DStaskAccbyTrial, na.rm=T)


NStaskAccbyTrial <- sapply(seq(1, 35), function(x){
        sapply(seq(1, length(NStaskAcc)), function(j){
                NStaskAcc[[j]][x]
        })
})
NStaskAccMeanTrial <- colMeans(NStaskAccbyTrial, na.rm=T)
NStaskAccSEMTrial <- colSds(NStaskAccbyTrial, na.rm=T)

#5 trial bin performance on the extinction test

binCuts <- seq(1, length(DStaskAccMeanTrial), by=5)
binIndex <- 1:length(DStaskAccMeanTrial)

binAssignment <- findInterval(binIndex, binCuts)


#DS SPECIFICITY
DStaskAcc5trialBins <- sapply(seq(1, length(binCuts)), function(m){
        mean(DStaskAccMeanTrial[binAssignment==m], na.rm=T)
})

DStaskAcc5trialBinsSEM <- sapply(seq(1, length(binCuts)), function(m){
        sd(DStaskAccMeanTrial[binAssignment==m], na.rm=T)/sqrt(sum(binAssignment==m))
})


#NS SPECIFICITY
NStaskAcc5trialBins <- sapply(seq(1, length(binCuts)), function(m){
        mean(NStaskAccMeanTrial[binAssignment==m], na.rm=T)
})

NStaskAcc5trialBinsSEM <- sapply(seq(1, length(binCuts)), function(m){
        sd(NStaskAccMeanTrial[binAssignment==m], na.rm=T)/sqrt(sum(binAssignment==m))
})



plot.new()
plot.window(xlim=c(0, length(binCuts)), ylim=c(-4, 6))

points(DStaskAcc5trialBins, pch=19, cex=2)
lines(DStaskAcc5trialBins, lwd=2)
errBars(x=seq(1, length(binCuts)), y=DStaskAcc5trialBins, err=DStaskAcc5trialBinsSEM, hatLength = 0.05)

points(NStaskAcc5trialBins, pch=19, cex=2, col="gray40")
lines(NStaskAcc5trialBins, lwd=2, col="gray40")
errBars(x=seq(1, length(binCuts)), y=NStaskAcc5trialBins, err=NStaskAcc5trialBinsSEM, hatLength = 0.05)

labVals <- c("1-5", "6-10", "11-15", "16-20", "21-25", "26-30", "31-35")
axis(side=1, at=seq(1, 7, by=1), labels=labVals, cex.axis=2, las=2)
axis(side=2, at=seq(-4, 6, by=1), cex.axis=2, las=2, pos=0.5)
abline(h=0, lty=3)

legend("topright", legend=c("S+", "S-"), col=c("black", "gray40"), lwd=2)

### For analysis:
DS_PerfIndex_Longform <- do.call("rbind", lapply(seq(1, nrow(DStaskAccbyTrial)), function(x){
        ratByBin <- sapply(unique(binAssignment), function(y){
                mean(DStaskAccbyTrial[x, binAssignment==y], na.rm=T)
        })
        data.frame(cue="S+", rat=rats[x], PerfIndex=ratByBin, bin=unique(binAssignment))
})
)

NS_PerfIndex_Longform <- do.call("rbind", lapply(seq(1, nrow(NStaskAccbyTrial)), function(x){
        ratByBin <- sapply(unique(binAssignment), function(y){
                mean(NStaskAccbyTrial[x, binAssignment==y], na.rm=T)
        })
        data.frame(cue="S-", rat=rats[x], PerfIndex=ratByBin, bin=unique(binAssignment))
})
)

PerfIndex_ForANOVA <- rbind(DS_PerfIndex_Longform, NS_PerfIndex_Longform)

ezANOVA(data=PerfIndex_ForANOVA, wid=rat, within = c(cue, bin), dv=PerfIndex) #Check this for questions about Mauchly's test and GG correction. https://www.r-exercises.com/2016/11/29/repeated-measures-anova-in-r-exercises/
# $`ANOVA`
# Effect DFn DFd            F            p p<.05       ges
# 1     cue   1   1 1.869764e+00 0.4019855461       0.4669608
# 2     bin   1   1 5.750234e+00 0.2515236555       0.7534574
# 3 cue:bin   1   1 2.221958e+06 0.0004270829     * 0.2145717

#With "aov" instead of "ezANOVA"
PerfIndex_ForANOVA$cue <- factor(PerfIndex_ForANOVA$cue)
PerfIndex_ForANOVA$bin <- factor(PerfIndex_ForANOVA$bin)

options(contrasts = c("contr.sum","contr.poly"))
aovtest_Perf_NL <- summary(with(PerfIndex_ForANOVA, 
                       aov(PerfIndex ~ cue*bin + Error(rat/(cue*bin)))
))

qqnorm(PerfIndex_ForANOVA$PerfIndex)
hist(PerfIndex_ForANOVA$PerfIndex)

save(aovtest_Perf_NL, file=paste(dataForRdir_L, "aovtest_Perf_NL.rdat", sep=""))






########## FR NON LEARNERS
#### NEURONAL DATA_Extract firing rate data related to cues and entries on BOTH SIDES (VEH and AP5)
allNeuronsDS_NL_VEH = neuralhist (funcdirect=funcdirect, path=NEXfiles_NL, event=1, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=T, side="vehicle")
allNeuronsNS_NL_VEH= neuralhist (funcdirect=funcdirect, path=NEXfiles_NL, event=2, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=T, side="vehicle")
allNeuronsDSresponded_NL_VEH= neuralhist (funcdirect=funcdirect, path=NEXfiles_NL, event=5, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsNSresponded_NL_VEH= neuralhist (funcdirect=funcdirect, path=NEXfiles_NL, event=7, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsDSmissed_NL_VEH= neuralhist (funcdirect=funcdirect, path=NEXfiles_NL, event=6, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")
allNeuronsNSmissed_NL_VEH= neuralhist (funcdirect=funcdirect, path=NEXfiles_NL, event=8, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="vehicle")

allNeuronsDS_NL_AP5 = neuralhist (funcdirect=funcdirect, path=NEXfiles_NL, event=1, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=T, side="drug")
allNeuronsNS_NL_AP5= neuralhist (funcdirect=funcdirect, path=NEXfiles_NL, event=2, startt=0, endt=10000, binw=50, psthmin=10, psthmax=10, cueexonly=F, allResults=T, side="drug")
allNeuronsDSresponded_NL_AP5= neuralhist (funcdirect=funcdirect, path=NEXfiles_NL, event=5, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsNSresponded_NL_AP5= neuralhist (funcdirect=funcdirect, path=NEXfiles_NL, event=7, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsDSmissed_NL_AP5= neuralhist (funcdirect=funcdirect, path=NEXfiles_NL, event=6, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")
allNeuronsNSmissed_NL_AP5= neuralhist (funcdirect=funcdirect, path=NEXfiles_NL, event=8, startt=0,endt=10000, binw=50, psthmin=10, psthmax=10,cueexonly=F, allResults=T, side="drug")

#sAVE THESE OBJECTS
save(allNeuronsDS_NL_VEH, file=paste(dataForRdir_NL, 'allNeuronsDS_NL_VEH.rdat', sep=""))
save(allNeuronsNS_NL_VEH, file=paste(dataForRdir_NL, 'allNeuronsNS_NL_VEH.rdat', sep=""))
save(allNeuronsDSresponded_NL_VEH, file=paste(dataForRdir_NL, 'allNeuronsDSresponded_NL_VEH.rdat', sep=""))
save(allNeuronsNSresponded_NL_VEH, file=paste(dataForRdir_NL, 'allNeuronsNSresponded_NL_VEH.rdat', sep=""))
save(allNeuronsDSmissed_NL_VEH, file=paste(dataForRdir_NL, 'allNeuronsDSmissed_NL_VEH.rdat', sep=""))
save(allNeuronsNSmissed_NL_VEH, file=paste(dataForRdir_NL, 'allNeuronsNSmissed_NL_VEH.rdat', sep=""))
save(allNeuronsDS_NL_AP5, file=paste(dataForRdir_NL, 'allNeuronsDS_NL_AP5.rdat', sep=""))
save(allNeuronsNS_NL_AP5, file=paste(dataForRdir_NL, 'allNeuronsNS_NL_AP5.rdat', sep=""))
save(allNeuronsDSresponded_NL_AP5, file=paste(dataForRdir_NL, 'allNeuronsDSresponded_NL_AP5.rdat', sep=""))
save(allNeuronsNSresponded_NL_AP5, file=paste(dataForRdir_NL, 'allNeuronsNSresponded_NL_AP5.rdat', sep=""))
save(allNeuronsDSmissed_NL_AP5, file=paste(dataForRdir_NL, 'allNeuronsDSmissed_NL_AP5.rdat', sep=""))
save(allNeuronsNSmissed_NL_AP5, file=paste(dataForRdir_NL, 'allNeuronsNSmissed_NL_AP5.rdat', sep=""))

masterDF_DS_NL_VEH <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS_NL_VEH, 
                                        CPvector=rep(1, length(rats)), funcdirect=funcdirect, 
                                        dataForRdir=dataForRdir_NL, 
                                        BLduration=2, sessionCPperRat=rep(1, length(rats)), 
                                        BLneudata=allNeuronsDS_NL_VEH)
masterDF_DS_NL_AP5 <- FRbyNEURONbyBINcue(eventType="cue", cue="S+", neudata=allNeuronsDS_NL_AP5, 
                                        CPvector=rep(1, length(rats)), funcdirect=funcdirect, 
                                        dataForRdir=dataForRdir_NL, 
                                        BLduration=2, sessionCPperRat=rep(1, length(rats)), 
                                        BLneudata=allNeuronsDS_NL_AP5)
masterDF_NS_NL_VEH <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_NL_VEH, 
                                        CPvector=rep(1, length(rats)), funcdirect=funcdirect, 
                                        dataForRdir=dataForRdir_NL, 
                                        BLduration=2, sessionCPperRat=rep(1, length(rats)), 
                                        BLneudata=allNeuronsNS_NL_VEH)
masterDF_NS_NL_AP5 <- FRbyNEURONbyBINcue(eventType="cue", cue="S-", neudata=allNeuronsNS_NL_AP5, 
                                        CPvector=rep(1, length(rats)), funcdirect=funcdirect, 
                                        dataForRdir=dataForRdir_NL, 
                                        BLduration=2, sessionCPperRat=rep(1, length(rats)), 
                                        BLneudata=allNeuronsNS_NL_AP5)


save(masterDF_DS_NL_VEH, file=paste(dataForRdir, "masterDF_DS_NL_VEH.rdat", sep=""))
save(masterDF_DS_NL_AP5, file=paste(dataForRdir, "masterDF_DS_NL_AP5.rdat", sep=""))
save(masterDF_NS_NL_VEH, file=paste(dataForRdir, "masterDF_NS_NL_VEH.rdat", sep=""))
save(masterDF_NS_NL_AP5, file=paste(dataForRdir, "masterDF_NS_NL_AP5.rdat", sep=""))


### PLOT FR POST-CUE (100-400MS WINDOW) as a function of distance from 1st trial
### PLOT FR AROUND THE CUE (PTSH) 
#All units
SessionPSTH(experiment="Exp 4 Unil AP5 Non Learners EXT", 
            masterDF=list(masterDF_DS_NL_VEH, masterDF_DS_NL_AP5), 
            graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
            yAxMinZ = -2, yAxMaxZ = 8, yAxMaxRaw = 10, psthmin=0.5, psthmax=2, imgFormat="pdf", 
            neudata=list(allNeuronsDS_NL_VEH, allNeuronsDS_NL_AP5), cueExcOnly=F, 
            legendLabels=c("S+ VEH side NL", "S+ AP5 side NL"))

#Make data frame with info about proportion of cue-exc neurons (this I know from calculations inside the previous function, you can check the % cue exc. printed on the graph to calculate this)
ncueExcVEH <- 1
ncueExcAP5 <- 1
forChiSq <- data.frame(side=c("VEH", "AP5"), cueexc=c(ncueExcVEH, ncueExcAP5), noncueexc=c(26-ncueExcVEH, 15-ncueExcAP5))
chisq.test(forChiSq[, -1]) 
# Pearson's Chi-squared test with Yates' continuity correction
# data:  forChiSq[, -1]
# X-squared = 1.0839e-30, df = 1, p-value = 1

fisher.test(forChiSq[, -1]) # p = 1 


barplot(forChiSq$cueexc/(forChiSq$cueexc+forChiSq$noncueexc), col=colindx, border = NA, ylim=c(0, 1))

#Cue excited only
SessionPSTH(experiment="Exp 4 Unil AP5 Non Learners EXT Cue exc only", 
            masterDF=list(masterDF_DS_NL_VEH, masterDF_DS_NL_AP5), 
            graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
            yAxMinZ = -2, yAxMaxZ = 8, yAxMaxRaw = 10, psthmin=0.5, psthmax=2, imgFormat="pdf", 
            neudata=list(allNeuronsDS_NL_VEH, allNeuronsDS_NL_AP5), cueExcOnly=T, 
            legendLabels=c("S+ VEH side NL", "S+ AP5 side NL"))



#Boxplot for a specific window and analysis
SessionBOXPLOT(experiment="Exp 4 Unil AP5 Non learners EXT", masterDF=list(masterDF_DS_NL_VEH, masterDF_DS_NL_AP5), 
               graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
               yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, imgFormat="pdf", WdwStart = 100, WdwEnd = 400,
               neudata=list(allNeuronsDS_VEH_NL, allNeuronsDS_AP5_NL), cueExcOnly=F, 
               legendLabels=c("S+ VEH side NL", "S+ AP5 side NL"))
# Wilcoxon rank sum test with continuity correction
# 
# data:  first and second
# W = 239040, p-value = 0.1008
# alternative hypothesis: true location shift is greater than 0



SessionBOXPLOT(experiment="Exp 4 Unil AP5 Non learners EXT tail", masterDF=list(masterDF_DS_NL_VEH, masterDF_DS_NL_AP5), 
               graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
               yAxMinZ = -2, yAxMaxZ = 10, yAxMaxRaw = 10, imgFormat="pdf", WdwStart = 750, WdwEnd = 2000,
               neudata=list(allNeuronsDS_VEH_NL, allNeuronsDS_AP5_NL), cueExcOnly=F, 
               legendLabels=c("S+ VEH side NL", "S+ AP5 side NL"))

# [[2]]
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  first and second
# W = 233180, p-value = 0.3166
# alternative hypothesis: true location shift is greater than 0

