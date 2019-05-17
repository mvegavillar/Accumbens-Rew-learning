

#############################################################
### EXPERIMENT 5A: EARLY AP5 VS VEH TEST                  ###
#############################################################

### LOAD IMPORTANT LIBRARIES
install.packages("matrixStats")
install.packages('ez')

library(matrixStats)
library(ez)

Exp1folder <- "E:/Dropbox/NMDA/EXP5_Performance/"

##########################
##########################
### LOAD FUNCTIONS     ###
##########################
##########################

funcdirect <- "E:/Dropbox/NMDA/R functions/"
CPfuncFolder <- paste(funcdirect, 'Change_Point-master/', sep="")

#Load functions
load(file=paste(funcdirect, "MedPCextract.r", sep=""))
load(file=paste(funcdirect, "mpcextract_blockSingle.Rfunc", sep=""))
load(file=paste(funcdirect, "CPextract.r", sep=""))
load(file=paste(funcdirect, "neuralhist.r", sep=""))
load(file=paste(funcdirect, "FRbyNEURONbyBINcue.r", sep=""))
load(file=paste(funcdirect, "errBars.r", sep=""))
load(file=paste(funcdirect, "errCloud.r", sep=""))
load(file=paste(funcdirect, "psthInf.r", sep=""))
load(file=paste(funcdirect, "BinIndexCalculator.R", sep=""))
load(file=paste(funcdirect, "dotplot.r", sep=""))
load(file=paste(funcdirect, "KC.sigbins.R", sep=""))
load(file=paste(funcdirect, "KC.inhib.sigbins.R", sep=""))
load(file=paste(funcdirect, "prePostInf_FR.r", sep=""))

###########################
###########################
### DEFINE FOLDERS      ###
###########################
###########################


# Define folders for one group *OR* the other before running the rest of the code (not both because then you'll just rewrite the folders you defined for the first group)

### EARLY VEH #################################################################
subTestFolder <- paste(Exp1folder, "Early VEH/", sep="")
datafolder <- paste(subTestFolder, "MedPC files/", sep="")
dataForRdir <- paste(subTestFolder, "Data for R/", sep="")
dataForRCumulative <- paste(subTestFolder, "Data for R cumulative/", sep="")
dataForRCumulativeEarlyVEH <- dataForRCumulative
behGraphFolder <- paste(subTestFolder, "Graphs/Behavior/", sep="")
neuGraphFolder <- paste(Exp1folder, "Graphs/Neuronal/", sep="")
MixedGraphFolder <- paste(subTestFolder, "Graphs/Mixed/", sep="")
CPGraphFolder <- paste(subTestFolder, "Graphs/Behavior/Change point/", sep="")
NEXfiles <- paste(subTestFolder, "NEX files/", sep="")
preVsPostFRFolder <- "E:/Dropbox/NMDA/EXP5_Performance/Graphs/Neuronal/FR pre vs post scatterplot/"


### EARLY AP5 #################################################################
subTestFolder <- paste(Exp1folder, "Early AP5/", sep="")
datafolder <- paste(subTestFolder, "MedPC files/", sep="")
dataForRdir <- paste(subTestFolder, "Data for R/", sep="")
dataForRCumulative <- paste(subTestFolder, "Data for R cumulative/", sep="")
dataForRCumulativeEarlyAP5 <- dataForRCumulative
behGraphFolder <- paste(subTestFolder, "Graphs/Behavior/", sep="")
MixedGraphFolder <- paste(subTestFolder, "Graphs/Mixed/", sep="")
CPGraphFolder <- paste(subTestFolder, "Graphs/Behavior/Change point/", sep="")
NEXfiles <- paste(subTestFolder, "NEX files/", sep="")


################################
################################
### CREATE IMPORTANT OBJECTS ###
################################
################################

# Extract behavioral data from MedPC files. This function saves the generated objects in the "dataForRdir". You have to load them (see next line of code) to bring them to your environment.
# This will give you a few error messages if, in any file, the first cue comes on after 5s of session onset. Ignore it, it just assigns NA to that trial, which is what you want.
# The parameter 'consumeRewWdw' is just the segment of the ITI that we discard (for ITI latency calculations) bc we assume that, if the animal got a reward on the previous trial, he might still be consuming the reward.
MedPCextract(MovAvg="Impinged only", cuelength=10, funcdirect = funcdirect, datafolder = datafolder, dataForRdir = dataForRdir, dataForRCumulative=dataForRCumulative)

# Load the behavior-related objects that you generated with the previous function. The main objects that we loaded are 'alldata' (detailed data by session by animal) and 'csacqidx' (an index of all the files). Name all sessions on csacqidx the same (i.e. '1')
files <- paste(dataForRdir, list.files(dataForRdir), sep=""); for(i in 1:length(files)){load(files[[i]])}
filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep=""); for(i in 1:length(filesCum)){load(filesCum[[i]])}

#This function will create the following objects: DSbinIdx, NSbinIdx and AllCueBinIdx. These are indexes indicating, for each rat and each kind of event, to what bin the events belong. The were generated and saved in the dataForRCumulative folder, so I need to load them
binsize <- 600
BinIndexCalculator(data=alldata, binsize=binsize, sessLength = 9000); filesCum <- paste(dataForRCumulative, list.files(dataForRCumulative), sep=""); for(i in 1:length(filesCum)){load(filesCum[[i]])}

# Create an object with data per bin for each one of our behavioral parameters
# Response ratio:
minBinNo <- min(sapply(DSbinIdx, max))

DSrespRatioByBin <- lapply(seq(1, length(DSrespAll)), function(x){
        sapply(seq(1, minBinNo), function(y){
                DSinBin <- DSrespAll[[x]][DSbinIdx[[x]]==y]
                DSrespRatio <- sum(DSinBin)/length(DSinBin)
        })
})

NSrespRatioByBin <- lapply(seq(1, length(NSrespAll)), function(x){
        sapply(seq(1, minBinNo), function(y){
                NSinBin <- NSrespAll[[x]][NSbinIdx[[x]]==y]
                NSrespRatio <- sum(NSinBin)/length(NSinBin)
        })
})


# Latency:
DSlatencyByBin <- lapply(seq(1, length(DSlatency)), function(x){
        sapply(seq(1, minBinNo), function(y){
                DSinBin <- DSlatency[[x]][DSbinIdx[[x]]==y]
                DSlatencyByBin <- mean(DSinBin, na.rm=T)
        })
})

NSlatencyByBin <- lapply(seq(1, length(NSlatency)), function(x){
        sapply(seq(1, minBinNo), function(y){
                NSinBin <- NSlatency[[x]][NSbinIdx[[x]]==y]
                NSlatencyByBin <- mean(NSinBin, na.rm=T)
        })
})

# Task Accuracy
DStaskAccByBin <- lapply(seq(1, length(DStaskAcc)), function(x){
        sapply(seq(1, minBinNo), function(y){
                DSinBin <- DStaskAcc[[x]][DSbinIdx[[x]]==y]
                DStaskAccByBin <- mean(DSinBin, na.rm=T)
        })
})

NStaskAccByBin <- lapply(seq(1, length(NStaskAcc)), function(x){
        sapply(seq(1, minBinNo), function(y){
                NSinBin <- NStaskAcc[[x]][NSbinIdx[[x]]==y]
                NStaskAccByBin <- mean(NSinBin, na.rm=T)
        })
})

# ITI latency
ITIlatByBin <- lapply(seq(1, length(ITIlatency)), function(x){
        sapply(seq(1, minBinNo), function(y){
                ITIlatInBin <- ITIlatency[[x]][AllCueBinIdx[[x]]==y]
                ITIlatByBin <- mean(ITIlatInBin, na.rm=T)
        })
})

# Make an object with the DStaskAccByBin for later statistical analyses
# Change "drug="VEH"" when I run this for AP5
DStaskAccByBin_LongFormat <- do.call("rbind", lapply(seq(1, length(DStaskAccByBin)), function(k){
        rat=as.character(rats[[k]])
        a <- 1:length(DStaskAccByBin[[k]])
        data.frame(rat=rat, bin=a, drug="VEH", perf=DStaskAccByBin[[k]], index="DS")
}))


NStaskAccByBin_LongFormat <- do.call("rbind", lapply(seq(1, length(NStaskAccByBin)), function(k){
        rat=as.character(rats[[k]])
        a <- 1:length(NStaskAccByBin[[k]])
        data.frame(rat=rat, bin=a, drug="VEH", perf=NStaskAccByBin[[k]], index="NS")
}))



### Make a long-format object with all these data for statistical analyses

#Run all of the above lines FIRST for VEH rats and then this line:
# Let's create objects to help us select the bins of interest for the pre and the post
# The infusion took place after 30min and it lasted 12min. I'm going to use the 30min before the infusion as baseline and the 30min after the infusion as the post.
PreInfLength <- 30*60            #In sec, baseline period (in my MedPC code it's always 30min)
PostInfStart <- (30*60)+12*60    #In sec, time of infusion end (when my postinfusion period starts)
PostInfEnd <- PostInfStart+30*60 #In sec, end of the window of interest after infusion (I made it 30min to match BL)

BLbinIndex <- (1:minBinNo)[1:(PreInfLength/binsize)] #Bins that correspond with the baseline period (the 30min before the infusion). 
PostInfBinIndex <- (1:minBinNo)[ceiling(PostInfStart/binsize):(PostInfEnd/binsize)] #Bins that correspond with the postinfusion period I want to study. In this case, the 30min after infusion.

byBinDataEarlyVEH <- list(DSrespRatioByBin, NSrespRatioByBin, DSlatencyByBin, NSlatencyByBin, DStaskAccByBin, NStaskAccByBin, ITIlatByBin)
IndexLabel <- c("S+.RR", "S-.RR", "S+.Latency", "S-.Latency", "S+.Spec.", "S-.Spec.", "ITI.Latency.")

EarlyVEH_LongFormat <- do.call("rbind", lapply(seq(1, length(byBinDataEarlyVEH)), function(x){ #For each index
        mat <- do.call("rbind", byBinDataEarlyVEH[[x]])
        if(length(BLbinIndex)>1){
                BLmean <- rowMeans(mat[,BLbinIndex], na.rm=T) #Mean by subject PRE infusion  
        } else {
                BLmean <- mean(mat[,BLbinIndex], na.rm=T) #Mean by subject PRE infusion
        }
        
        PostMean <- rowMeans(mat[,PostInfBinIndex], na.rm=T) #Mean by subject POST infusion
        ratnames <- paste("VEH", 1:nrow(mat), sep="_")
        return(data.frame(Drug="VEH", Rat=ratnames, Index=IndexLabel[x], Infusion=c(rep("Pre", nrow(mat)), rep("Post", nrow(mat))), Performance=c(BLmean, PostMean)))
}))

EarlyVEH_DStaskAccByBin_LongFormat <- DStaskAccByBin_LongFormat
EarlyVEH_NStaskAccByBin_LongFormat <- NStaskAccByBin_LongFormat

save(EarlyVEH_LongFormat, file=paste(dataForRdir, "EarlyVEH_LongFormat.rdat", sep=""))
save(EarlyVEH_DStaskAccByBin_LongFormat, file=paste(dataForRdir, "EarlyVEH_DStaskAccByBin_LongFormat.rdat", sep=""))
save(EarlyVEH_NStaskAccByBin_LongFormat, file=paste(dataForRdir, "EarlyVEH_NStaskAccByBin_LongFormat.rdat", sep=""))


#Then repeat for AP5 rats and run these lines
byBinDataEarlyAP5 <- list(DSrespRatioByBin, NSrespRatioByBin, DSlatencyByBin, NSlatencyByBin, DStaskAccByBin, NStaskAccByBin, ITIlatByBin)
IndexLabel <- c("S+.RR", "S-.RR", "S+.Latency", "S-.Latency", "S+.Spec.", "S-.Spec.", "ITI.Latency.")

EarlyAP5_LongFormat <- do.call("rbind", lapply(seq(1, length(byBinDataEarlyAP5)), function(x){ #For each index
        mat <- do.call("rbind", byBinDataEarlyAP5[[x]])
        BLmean <- rowMeans(mat[,BLbinIndex], na.rm=T) #Mean by subject PRE infusion
        PostMean <- rowMeans(mat[,PostInfBinIndex], na.rm=T) #Mean by subject POST infusion
        ratnames <- paste("AP5", 1:nrow(mat), sep="_")
        return(data.frame(Drug="AP5", Rat=ratnames, Index=IndexLabel[x], Infusion=c(rep("Pre", nrow(mat)), rep("Post", nrow(mat))), Performance=c(BLmean, PostMean)))
}))

EarlyAP5_DStaskAccByBin_LongFormat <- DStaskAccByBin_LongFormat
EarlyAP5_NStaskAccByBin_LongFormat <- NStaskAccByBin_LongFormat


save(EarlyAP5_LongFormat, file=paste(dataForRdir, "EarlyAP5_LongFormat.rdat", sep=""))
save(EarlyAP5_DStaskAccByBin_LongFormat, file=paste(dataForRdir, "EarlyAP5_DStaskAccByBin_LongFormat.rdat", sep=""))
save(EarlyAP5_NStaskAccByBin_LongFormat, file=paste(dataForRdir, "EarlyAP5_NStaskAccByBin_LongFormat.rdat", sep=""))

###

Early_LongFormat <- rbind(EarlyVEH_LongFormat, EarlyAP5_LongFormat)
Early_LongFormatByBin <- rbind(EarlyVEH_DStaskAccByBin_LongFormat, EarlyAP5_DStaskAccByBin_LongFormat)
Early_LongFormatByBin_DSandNS <- rbind(EarlyVEH_DStaskAccByBin_LongFormat, EarlyVEH_NStaskAccByBin_LongFormat,
                                       EarlyAP5_DStaskAccByBin_LongFormat, EarlyAP5_NStaskAccByBin_LongFormat)

save(Early_LongFormat, file=paste(dataForRdir, "Early_LongFormat.rdat", sep=""))
save(Early_LongFormatByBin, file=paste(dataForRdir, "Early_LongFormatByBin.rdat", sep=""))
save(Early_LongFormatByBin_DSandNS, file=paste(dataForRdir, "Early_LongFormatByBin_DSandNS.rdat", sep=""))



# Extract neuronal data from NEX files. 

#VEH test: data aligned to DS and NS onset BEFORE and AFTER the infusion.
postInfTargetWdw <- 1800+12*60+30*60 #For the post infusion window, I'll choose the period between the end of the infusion +30'.
allNeuronsDSEarlyVEHPreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsDSEarlyVEHPostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSEarlyVEHPreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSEarlyVEHPostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryDSEarlyVEHPreInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryDSEarlyVEHPostInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryNSEarlyVEHPreInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryNSEarlyVEHPostInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryITIEarlyVEHPreInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryITIEarlyVEHPostInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")


#AP5 test: data aligned to DS and NS onset BEFORE and AFTER the infusion.
postInfTargetWdw <- 1800+12*60+30*60 #For the post infusion window, I'll choose the period between the end of the infusion +30'.
allNeuronsDSEarlyAP5PreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsDSEarlyAP5PostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=1, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSEarlyAP5PreInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsNSEarlyAP5PostInf <-  neuralhist (funcdirect=funcdirect, path=NEXfiles, event=2, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryDSEarlyAP5PreInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryDSEarlyAP5PostInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=9, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryNSEarlyAP5PreInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryNSEarlyAP5PostInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=14, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryITIEarlyAP5PreInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=0, endt=1800, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")
allNeuronsEntryITIEarlyAP5PostInf <- neuralhist (funcdirect=funcdirect, path=NEXfiles, event=10, startt=2520, endt=postInfTargetWdw, binw=50, psthmin=2, psthmax=10, cueexonly=F, allResults=T, side="both")


### GIVE THESE OBJECTS A UNIQUE NAME

## VEH SIDE
# csacqidxEarlyVEH <- csacqidx
# alldataEarlyVEH <- alldata
# ratsEarlyVEH <- rats
# idxEarlyVEH <- idx
# cumDataEarlyVEH <- list(DSrespAll, DStaskAcc, DStimeToSpare, NSrespAll, NStaskAcc, NStimeToSpare)
# 
# ## AP5 SIDE
# csacqidxEarlyAP5 <- csacqidx
# alldataEarlyAP5 <- alldata
# ratsEarlyAP5 <- rats
# idxEarlyAP5 <- idx
# cumDataEarlyAP5 <- list(DSrespAll, DStaskAcc, DStimeToSpare, NSrespAll, NStaskAcc, NStimeToSpare)



######################################################
######################################################
### PLOT GRAPHS                                    ###
######################################################
######################################################


###################
### 1. BEHAVIOR ###
###################

# Let's create objects to help us select the bins of interest for the pre and the post
# The infusion took place after 30min and it lasted 12min. I'm going to use the 30min before the infusion as baseline and the 30min after the infusion as the post.
PreInfLength <- 30*60            #In sec
PostInfStart <- (30*60)+12*60    #In sec
PostInfEnd <- PostInfStart+30*60 #In sec

BLbinIndex <- (1:minBinNo)[1:(PreInfLength/binsize)]
PostInfBinIndex <- (1:minBinNo)[ceiling(PostInfStart/binsize):(PostInfEnd/binsize)]

#Function for plotting lines more easily. I just need to adjust the data I feed the function, the color and the points
plotPrePostLines <- function(data, color, pch, scores, jitter=0){
        mat <- do.call("rbind", data) #Create matrix in which rows are different rats and columns are bins
        
        if(scores=="absolute"){
                BLmean <- rowMeans(mat[,BLbinIndex], na.rm=T) #Mean by subject PRE infusion
                PostMean <- rowMeans(mat[,PostInfBinIndex], na.rm=T) #Mean by subject POST infusion
        }
        
        if(scores=="percentBL"){
                BLmeanAll <- mean(rowMeans(mat[,BLbinIndex], na.rm=T), na.rm=T) #Mean all subjects PRE infusion
                PostMeanAll <- mean(rowMeans(mat[,PostInfBinIndex], na.rm=T), na.rm=T) #Mean all subjects POST infusion
                
                BLmeanEach <- rowMeans(mat[,BLbinIndex], na.rm=T) #Mean by subject PRE infusion
                PostMeanEach <- rowMeans(mat[,PostInfBinIndex], na.rm=T) #Mean by subject POST infusion
                
                BLmean <- (BLmeanEach/BLmeanEach)*100 #Mean by subject PRE infusion in terms of percentage of BL performance of that same subject (it has to be 100%)
                PostMean <- (PostMeanEach/BLmeanEach)*100 #Mean by subject POST infusion in terms of percentage of BL performance 
        }
        
        lines(x=c(0, 1), y=c(mean(BLmean), mean(PostMean)), col=color, cex=2)
        errBars(x=c(0, 1), y=c(mean(BLmean), mean(PostMean)), err=c(sd(BLmean)/sqrt(length(BLmean)), sd(PostMean)/sqrt(length(PostMean))), color=color, jitter=jitter)
        points(x=c(0, 1), y=c(mean(BLmean), mean(PostMean)), pch=pch, col=color, cex=2)
        if(pch==22){points(x=c(0, 1), y=c(mean(BLmean), mean(PostMean)), pch=pch, col=color, cex=2, bg="white")}
        
}


#Function for plotting bars more easily. I just need to adjust the data I feed the function, the color and the points
plotPrePostBars <- function(data, color, xmiddle, barwidth, labelY, colLabel){
        mat <- do.call("rbind", data) #Create matrix in which rows are different rats and columns are bins
        BLmean <- rowMeans(mat[,BLbinIndex], na.rm=T) #Mean by subject PRE infusion
        PostMean <- rowMeans(mat[,PostInfBinIndex], na.rm=T) #Mean by subject POST infusion
        
        #Pre
        rect(xleft=xmiddle-barwidth, xright = xmiddle, ybottom=0, ytop=mean(BLmean), col=color, border="white")
        text(x=xmiddle-barwidth/2, y=labelY, labels = "Pre", col=colLabel, font=2)
        #Post
        rect(xleft=xmiddle, xright=xmiddle+barwidth, ybottom=0, ytop=mean(PostMean), col=color, border="white")
        text(x=xmiddle+barwidth/2, y=labelY, labels = "Post", col=colLabel, font=2)
        #Individual lines
        for(i in 1:length(data)){lines(x=c(xmiddle-barwidth/2, xmiddle+barwidth/2), y=c(BLmean[i], PostMean[i]))}
}

#Define colors
colindx <- c("#2171b5", "#cb181d") #Strong blue and red
colindxB <- c("#bdd7e7", "#fcae91") #Less strong blue and red
colindxC <- c("#eff3ff", "#fb6a4a") #Even less strong blue and red
colindxD <- c("#6baed6", "#fee5d9") #Lightest blue and red
     


#### 1.1. RESPONSE RATIO

### 1.1.1. Response ratio: S+ and S- responding pre vs. post infusion in AP5 vs. VEH
# In the objects 'byBinDataEarlyVEH' and 'byBinDataEarlyAP5', the first and second items are DSrespratio and NSrespratio by subject by bin

## 1.1.1.1. Absolute scores
plot.new()
par(oma=c(2,2,2,2))
plot.window(xlim=c(0, 1), ylim=c(0, 1))

plotPrePostLines(data=byBinDataEarlyVEH[[1]], color=colindx[1], pch=15, scores="absolute") #VEH group, S+
plotPrePostLines(data=byBinDataEarlyAP5[[1]], color=colindx[2], pch=15, scores="absolute") #AP5 group, S+

plotPrePostLines(data=byBinDataEarlyVEH[[2]], color=colindx[1], pch=22, scores="absolute", jitter=0.015) #VEH group, S-
plotPrePostLines(data=byBinDataEarlyAP5[[2]], color=colindx[2], pch=22, scores="absolute", jitter=0.015) #AP5 group, S-

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 1, by=0.2, labels=seq(0, 1, 0.2)), font=2, las=2, pos=-0.1)
mtext(side=2, line=4, text="Proportion", cex=1.4, font=2)

legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n")



## 1.1.1.2. Percentage of BL 
plot.new()
par(oma=c(2,2,2,2))
plot.window(xlim=c(0, 1), ylim=c(0, 120))

plotPrePostLines(data=byBinDataEarlyVEH[[1]], color=colindx[1], pch=15, scores="percentBL") #VEH group, S+
plotPrePostLines(data=byBinDataEarlyAP5[[1]], color=colindx[2], pch=15, scores="percentBL") #AP5 group, S+

#plotPrePostLines(data=byBinDataEarlyVEH[[2]], color=colindx[1], pch=22, scores="percentBL", jitter=0.015) #VEH group, S-. It's confusing so I'm not plotting it
#plotPrePostLines(data=byBinDataEarlyAP5[[2]], color=colindx[2], pch=22, scores="percentBL", jitter=0.015) #AP5 group, S-

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 120, by=20), labels=seq(0, 120, 20), font=2, las=2, pos=-0.1)
mtext(side=2, line=4, text="% of BL response ratio", font=2, cex.axis=1.5)

#legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n")


### 1.1.2. Response ratio: S+ and S- responding by bin on test day in AP5 vs. VEH

## 1.1.2.1. Absolute scores
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 1))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=0, ytop=1.2, col="gray95", border="white")

#lapply(seq(1, length(ratsEarlyVEH)), function(x) {lines(byBinDataEarlyVEH[[1]][[x]], col=colindx[1])})
#lapply(seq(1, length(ratsEarlyAP5)), function(x) {lines(byBinDataEarlyAP5[[1]][[x]], col=colindx[2])})

matVEH <- do.call("rbind", byBinDataEarlyVEH[[1]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataEarlyAP5[[1]])

lines(colMeans(matVEH), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEH), err=colSds(matVEH)/sqrt(nrow(matVEH)), color=colindx[1])
points(colMeans(matVEH), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5), err=colSds(matAP5)/sqrt(nrow(matAP5)), color=colindx[2])
points(colMeans(matAP5), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(0, 1, by=0.2), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.4)
mtext(side=2, line=2, text="Proportion", font=2, cex=1.4)

legend("bottomright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)



## 1.1.2.2. Percentage of BL
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 120))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=0, ytop=120, col="gray95", border="white")

#Get data ready
matVEH <- do.call("rbind", byBinDataEarlyVEH[[1]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataEarlyAP5[[1]])

VEHbl <- colMeans(matVEH[,BLbinIndex], na.rm=T); AP5bl <- colMeans(matAP5[,BLbinIndex], na.rm=T)

matVEHperc <- (matVEH/VEHbl)*100
matAP5perc <- (matAP5/AP5bl)*100

#Plot
lines(colMeans(matVEHperc), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEHperc), err=colSds(matVEHperc)/sqrt(nrow(matVEHperc)), color=colindx[1])
points(colMeans(matVEHperc), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5perc), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5perc), err=colSds(matAP5perc)/sqrt(nrow(matAP5perc)), color=colindx[2])
points(colMeans(matAP5perc), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(0, 120, by=20), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.4)
mtext(side=2, line=2, text="% of baseline", font=2, cex=1.4)

legend("bottomright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)


### 1.1.3. Response ratio: barplots of pre and post infusion, S+ vs S- and VEH vs AP5
plot.new()
plot.window(xlim=c(0, 6), ylim=c(0, 1))

#S+ both groups pre and post infusion
plotPrePostBars(data=byBinDataEarlyVEH[[1]], color=colindx[1], xmiddle=1, barwidth=0.5, colLabel = "white", labelY = 0.05)
plotPrePostBars(data=byBinDataEarlyAP5[[1]], color=colindx[2], xmiddle=2, barwidth=0.5, colLabel = "white", labelY = 0.05)

#S- both groups pre and post infusion
plotPrePostBars(data=byBinDataEarlyVEH[[2]], color=colindxB[1], xmiddle=3.25, barwidth=0.5, colLabel = "black", labelY = 0.05)
plotPrePostBars(data=byBinDataEarlyAP5[[2]], color=colindxB[2], xmiddle=4.25, barwidth=0.5, colLabel = "black", labelY = 0.05)

#Axis and labels
axis(side=1, tick = F, at=c(1.5, 3.75), labels=c("S+", "S-"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 1, 0.2), cex.axis=1, font=2, las=2)
mtext(side=2, line=2.5, text="Response ratio", cex=1.4, font=2)
rect(xleft=3, xright=3.5, ybottom=0.95, ytop=1, col=colindx[1], border="white")
rect(xleft=3.5, xright=4, ybottom=0.95, ytop=1, col=colindxB[1], border="white")
rect(xleft=3, xright=3.5, ybottom=0.85, ytop=0.9, col=colindx[2], border="white")
rect(xleft=3.5, xright=4, ybottom=0.85, ytop=0.9, col=colindxB[2], border="white")
text(x=4.5, y=0.98, labels="VEH", cex=1.5)
text(x=4.5, y=0.88, labels="AP5", cex=1.5)



#### 1.2. CUED LATENCY
# In the objects 'byBinDataEarlyVEH' and 'byBinDataEarlyAP5', the 3rd and 4th items are DSlatency and NSlatency by subject by bin

### 1.2.1. Cued latency: S+ and S- latency pre vs. post infusion in AP5 vs. VEH

## 1.2.1.1. Absolute scores
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 10))

plotPrePostLines(data=byBinDataEarlyVEH[[3]], color=colindx[1], pch=15, scores="absolute") #VEH group, S+
plotPrePostLines(data=byBinDataEarlyAP5[[3]], color=colindx[2], pch=15, scores="absolute") #AP5 group, S+

plotPrePostLines(data=byBinDataEarlyVEH[[4]], color=colindx[1], pch=22, scores="absolute", jitter=0.015) #VEH group, S-
plotPrePostLines(data=byBinDataEarlyAP5[[4]], color=colindx[2], pch=22, scores="absolute", jitter=0.015) #AP5 group, S-

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 10, by=2, labels=seq(0, 10, 2)), font=2, las=2, pos=-0.05)

mtext(side=2, line=2.5, text="Latency (s)", cex=1.4, font=2)

legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n")


## 1.2.1.2. Percentage of BL 
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 400))

plotPrePostLines(data=byBinDataEarlyVEH[[3]], color=colindx[1], pch=15, scores="percentBL") #VEH group, S+
plotPrePostLines(data=byBinDataEarlyAP5[[3]], color=colindx[2], pch=15, scores="percentBL") #AP5 group, S+

#plotPrePostLines(data=byBinDataEarlyVEH[[2]], color=colindx[1], pch=22, scores="percentBL", jitter=0.015) #VEH group, S-. It's confusing so I'm not plotting it
#plotPrePostLines(data=byBinDataEarlyAP5[[2]], color=colindx[2], pch=22, scores="percentBL", jitter=0.015) #AP5 group, S-

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 400, by=50), labels=seq(0, 400, 50), font=2, las=2, pos=-0.05)
mtext(side=2, line=3, text="% of BL latency", font=2, cex=1.4)

#legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n", cex=1.5)


### 1.2.2. Cued latency: S+ and S- responding by bin on test day in AP5 vs. VEH

## 1.2.2.1. Absolute scores
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 10))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=0, ytop=10, col="gray95", border="white")

#lapply(seq(1, length(ratsEarlyVEH)), function(x) {lines(byBinDataEarlyVEH[[1]][[x]], col=colindx[1])})
#lapply(seq(1, length(ratsEarlyAP5)), function(x) {lines(byBinDataEarlyAP5[[1]][[x]], col=colindx[2])})

matVEH <- do.call("rbind", byBinDataEarlyVEH[[3]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataEarlyAP5[[3]])

lines(colMeans(matVEH), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEH), err=colSds(matVEH)/sqrt(nrow(matVEH)), color=colindx[1])
points(colMeans(matVEH), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5), err=colSds(matAP5)/sqrt(nrow(matAP5)), color=colindx[2])
points(colMeans(matAP5), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(0, 10, by=2), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.4)
mtext(side=2, line=1, text="Latency (s)", font=2, cex=1.4)

legend("topright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)

## 1.2.2.2. Percentage of BL
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 300))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=0, ytop=300, col="gray95", border="white")

#Get data ready
matVEH <- do.call("rbind", byBinDataEarlyVEH[[3]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataEarlyAP5[[3]])

VEHbl <- colMeans(matVEH[,BLbinIndex], na.rm=T); AP5bl <- colMeans(matAP5[,BLbinIndex], na.rm=T)

matVEHperc <- (matVEH/VEHbl)*100
matAP5perc <- (matAP5/AP5bl)*100

#Plot
lines(colMeans(matVEHperc), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEHperc), err=colSds(matVEHperc)/sqrt(nrow(matVEHperc)), color=colindx[1])
points(colMeans(matVEHperc), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5perc), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5perc), err=colSds(matAP5perc)/sqrt(nrow(matAP5perc)), color=colindx[2])
points(colMeans(matAP5perc), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(0, 300, by=50), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.2)
mtext(side=2, line=1, text="% of baseline latency", font=2, cex=1.2)

legend("topright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)

### 1.2.3. Cued latency: barplots of pre and post infusion, S+ vs S- and VEH vs AP5
plot.new()
plot.window(xlim=c(0, 6), ylim=c(0, 10))

#S+ both groups pre and post infusion
plotPrePostBars(data=byBinDataEarlyVEH[[3]], color=colindx[1], xmiddle=1, barwidth=0.5, colLabel = "white", labelY = 0.5)
plotPrePostBars(data=byBinDataEarlyAP5[[3]], color=colindx[2], xmiddle=2, barwidth=0.5, colLabel = "white", labelY = 0.5)

#S- both groups pre and post infusion
plotPrePostBars(data=byBinDataEarlyVEH[[4]], color=colindxB[1], xmiddle=3.25, barwidth=0.5, colLabel = "black", labelY = 0.5)
plotPrePostBars(data=byBinDataEarlyAP5[[4]], color=colindxB[2], xmiddle=4.25, barwidth=0.5, colLabel = "black", labelY = 0.5)

#Axis and labels
axis(side=1, tick = F, at=c(1.5, 3.75), labels=c("S+", "S-"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 10, 2), cex.axis=1, font=2, las=2)
mtext(side=2, line=2.5, text="Latency (s)", cex=1.4, font=2)
rect(xleft=0, xright=0.5, ybottom=9.5, ytop=10, col=colindx[1], border="white")
rect(xleft=0.5, xright=1, ybottom=9.5, ytop=10, col=colindxB[1], border="white")
rect(xleft=0, xright=0.5, ybottom=8.5, ytop=9, col=colindx[2], border="white")
rect(xleft=0.5, xright=1, ybottom=8.5, ytop=9, col=colindxB[2], border="white")
text(x=1.5, y=9.8, labels="VEH", cex=1.5)
text(x=1.5, y=8.8, labels="AP5", cex=1.5)



#### 1.3. ITI latency

### 1.3.1. ITI latency: ITI latency pre vs. post infusion in AP5 vs. VEH

## 1.3.1.1. Absolute scores
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 10))

plotPrePostLines(data=byBinDataEarlyVEH[[7]], color=colindx[1], pch=15, scores="absolute") #VEH group, S+
plotPrePostLines(data=byBinDataEarlyAP5[[7]], color=colindx[2], pch=15, scores="absolute") #AP5 group, S+

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 10, by=2), font=2, las=2, pos=-0.05)

mtext(side=2, line=2.5, text="ITI latency (s)", cex=1.4, font = 2)

legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n")

## 1.3.1.2. Percentage of BL 
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 150))

plotPrePostLines(data=byBinDataEarlyVEH[[7]], color=colindx[1], pch=15, scores="percentBL") #VEH group, S+
plotPrePostLines(data=byBinDataEarlyAP5[[7]], color=colindx[2], pch=15, scores="percentBL") #AP5 group, S+

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1, font=2)
axis(side=2, at=seq(0, 200, by=50), font=2, las=2, pos=-0.05)
mtext(side=2, line=3, text="% of BL ITI latency", font=2, cex.axis=1.5)

#legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomright", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n", cex=1.5)

### 1.3.2. ITI latency: ITI latency by bin on test day in AP5 vs. VEH

## 1.3.2.1. Absolute scores
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 10))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=0, ytop=10, col="gray95", border="white")

#lapply(seq(1, length(ratsEarlyVEH)), function(x) {lines(byBinDataEarlyVEH[[1]][[x]], col=colindx[1])})
#lapply(seq(1, length(ratsEarlyAP5)), function(x) {lines(byBinDataEarlyAP5[[1]][[x]], col=colindx[2])})

matVEH <- do.call("rbind", byBinDataEarlyVEH[[7]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataEarlyAP5[[7]])

lines(colMeans(matVEH), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEH), err=colSds(matVEH)/sqrt(nrow(matVEH)), color=colindx[1])
points(colMeans(matVEH), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5), err=colSds(matAP5)/sqrt(nrow(matAP5)), color=colindx[2])
points(colMeans(matAP5), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(0, 10, by=2), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.2)
mtext(side=2, line=1, text="ITI latency (s)", font=2, cex=1.2)

legend("bottomright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)


## 1.3.2.2. Percentage of BL
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 150))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=0, ytop=150, col="gray95", border="white")

#Get data ready
matVEH <- do.call("rbind", byBinDataEarlyVEH[[7]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataEarlyAP5[[7]])

VEHbl <- colMeans(matVEH[,BLbinIndex], na.rm=T); AP5bl <- colMeans(matAP5[,BLbinIndex], na.rm=T)

matVEHperc <- (matVEH/VEHbl)*100
matAP5perc <- (matAP5/AP5bl)*100

#Plot
lines(colMeans(matVEHperc), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEHperc), err=colSds(matVEHperc)/sqrt(nrow(matVEHperc)), color=colindx[1])
points(colMeans(matVEHperc), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5perc), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5perc), err=colSds(matAP5perc)/sqrt(nrow(matAP5perc)), color=colindx[2])
points(colMeans(matAP5perc), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(0, 150, by=50), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.2)
mtext(side=2, line=1, text="% of baseline ITI latency", font=2, cex=1.2)

legend("bottomright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)


### 1.3.3. ITI latency: barplots of pre and post infusion, S+ vs S- and VEH vs AP5
plot.new()
plot.window(xlim=c(0, 6), ylim=c(0, 10))

#S+ both groups pre and post infusion
plotPrePostBars(data=byBinDataEarlyVEH[[7]], color=colindx[1], xmiddle=1, barwidth=0.5, colLabel = "white", labelY = 0.5)
plotPrePostBars(data=byBinDataEarlyAP5[[7]], color=colindx[2], xmiddle=2, barwidth=0.5, colLabel = "white", labelY = 0.5)

#Axis and labels
axis(side=2, at=seq(0, 10, 2), cex.axis=1, font=2, las=2)
mtext(side=2, line=2.5, text="ITI latency (s)", cex=1.4, font=2)
rect(xleft=0, xright=0.5, ybottom=9.5, ytop=10, col=colindx[1], border="white")
rect(xleft=0, xright=0.5, ybottom=8.5, ytop=9, col=colindx[2], border="white")
text(x=1, y=9.8, labels="VEH", cex=1.5)
text(x=1, y=8.8, labels="AP5", cex=1.5)



#### 1.4. CUED SPECIFICITY

### 1.4.1. Cue specificity: S+ and S- specificity pre vs. post infusion in AP5 vs. VEH
# In the objects 'byBinDataEarlyVEH' and 'byBinDataEarlyAP5', the 5th and 6th items are DStaskAccuracy and NStaskAccuracy by subject by bin

## 1.4.1.1. Absolute scores
plot.new()
plot.window(xlim=c(0, 1), ylim=c(-2, 6))

abline(h=0, lty=3)

plotPrePostLines(data=byBinDataEarlyVEH[[5]], color=colindx[1], pch=15, scores="absolute") #VEH group, S+
plotPrePostLines(data=byBinDataEarlyAP5[[5]], color=colindx[2], pch=15, scores="absolute") #AP5 group, S+

plotPrePostLines(data=byBinDataEarlyVEH[[6]], color=colindx[1], pch=22, scores="absolute", jitter=0.015) #VEH group, S-
plotPrePostLines(data=byBinDataEarlyAP5[[6]], color=colindx[2], pch=22, scores="absolute", jitter=0.015) #AP5 group, S-

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(-2, 6, by=1), font=2, las=2, pos=-0.04)
mtext(side=2, line=2, text="S+ Specificity (s)", font=2, cex=1.4)

legend("topright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("topleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n")

## 1.4.1.2. Percentage of BL 
plot.new()
plot.window(xlim=c(0, 1), ylim=c(0, 140))

plotPrePostLines(data=byBinDataEarlyVEH[[5]], color=colindx[1], pch=15, scores="percentBL") #VEH group, S+
plotPrePostLines(data=byBinDataEarlyAP5[[5]], color=colindx[2], pch=15, scores="percentBL") #AP5 group, S+

#plotPrePostLines(data=byBinDataEarlyVEH[[2]], color=colindx[1], pch=22, scores="percentBL", jitter=0.015) #VEH group, S-. It's confusing so I'm not plotting it
#plotPrePostLines(data=byBinDataEarlyAP5[[2]], color=colindx[2], pch=22, scores="percentBL", jitter=0.015) #AP5 group, S-

axis(side=1, at=c(0, 1), labels=c("Preinfusion", "Postinfusion"), cex.axis=1.4, font=2)
axis(side=2, at=seq(0, 140, by=20), font=2, las=2, pos=-0.1)
mtext(side=2, line=4, text="% of BL S+ specificity", font=2, cex=1.4)

#legend("bottomright", legend = c("S+", "S-"), pch = c(15, 22), bty = "n" )
legend("bottomleft", legend=c("VEH", "AP5"), lty=1, col=colindx, bty="n", cex=1.2)


### 1.4.2.  Cue specificity: S+ and S- specificity by bin on test day in AP5 vs. VEH

## 1.4.2.1. Absolute scores
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(-2, 7))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=-2, ytop=6, col="gray95", border="white")
abline(h=0, lty=3)

#lapply(seq(1, length(ratsEarlyVEH)), function(x) {lines(byBinDataEarlyVEH[[1]][[x]], col=colindx[1])})
#lapply(seq(1, length(ratsEarlyAP5)), function(x) {lines(byBinDataEarlyAP5[[1]][[x]], col=colindx[2])})

matVEH <- do.call("rbind", byBinDataEarlyVEH[[5]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataEarlyAP5[[5]])

lines(colMeans(matVEH), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEH), err=colSds(matVEH)/sqrt(nrow(matVEH)), color=colindx[1])
points(colMeans(matVEH), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5), err=colSds(matAP5)/sqrt(nrow(matAP5)), color=colindx[2])
points(colMeans(matAP5), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(-2, 6, by=2), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.2)
mtext(side=2, line=1, text="S+ specificity (s)", font=2, cex=1.2)

legend("bottomright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)



## 1.4.2.1.B. s- Cued specificity. Absolute scores
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(-2, 7))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=-2, ytop=6, col="gray95", border="white")
abline(h=0, lty=3)

#lapply(seq(1, length(ratsEarlyVEH)), function(x) {lines(byBinDataEarlyVEH[[1]][[x]], col=colindx[1])})
#lapply(seq(1, length(ratsEarlyAP5)), function(x) {lines(byBinDataEarlyAP5[[1]][[x]], col=colindx[2])})

matVEH <- do.call("rbind", byBinDataEarlyVEH[[6]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataEarlyAP5[[6]])

lines(colMeans(matVEH), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEH), err=colSds(matVEH)/sqrt(nrow(matVEH)), color=colindx[1])
points(colMeans(matVEH), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5), err=colSds(matAP5)/sqrt(nrow(matAP5)), color=colindx[2])
points(colMeans(matAP5), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(-2, 6, by=2), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.2)
mtext(side=2, line=1, text="S+ specificity (s)", font=2, cex=1.2)

legend("bottomright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.5)





## 1.4.2.2. Percentage of BL
plot.new()
plot.window(xlim=c(0, minBinNo), ylim=c(0, 150))

#Mark infusion period
screenPerSec <- minBinNo/(12*binsize) #Length of one second in the X axis
infusionStart <- 1800; infusionEnd <- 1800+12*60
infusionStartScreen <- infusionStart*screenPerSec; infusionEndScreen <- infusionEnd*screenPerSec
rect(xleft=infusionStartScreen, xright=infusionEndScreen, ybottom=0, ytop=140, col="gray95", border="white")
#abline(h=100, lty=3)

#Get data ready
matVEH <- do.call("rbind", byBinDataEarlyVEH[[5]]) #Create matrix in which rows are different rats and columns are bins
matAP5 <- do.call("rbind", byBinDataEarlyAP5[[5]])

VEHbl <- colMeans(matVEH[,BLbinIndex], na.rm=T); AP5bl <- colMeans(matAP5[,BLbinIndex], na.rm=T)

matVEHperc <- (matVEH/VEHbl)*100
matAP5perc <- (matAP5/AP5bl)*100

#Plot
lines(colMeans(matVEHperc), col=colindx[1], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matVEHperc), err=colSds(matVEHperc)/sqrt(nrow(matVEHperc)), color=colindx[1])
points(colMeans(matVEHperc), col=colindx[1], pch=15, cex=1.5)


lines(colMeans(matAP5perc), col=colindx[2], lwd=2)
errBars(x=seq(1, minBinNo), y=colMeans(matAP5perc), err=colSds(matAP5perc)/sqrt(nrow(matAP5perc)), color=colindx[2])
points(colMeans(matAP5perc), col=colindx[2], pch=15, cex=1.5)

axis(side=1, at=seq(1, minBinNo, by=1), labels=seq(binsize/60, (minBinNo*binsize)/60, by=binsize/60), font=2)
axis(side=2, at=seq(0, 140, by=20), font=2, las=2, pos=0.5)
mtext(side=1, line=2.5, text = "Time (min)", font=2, cex=1.2)
mtext(side=2, line=1, text="% of baseline S+ specificity", font=2, cex=1.2)

legend("bottomright", legend=c("VEH", "AP5"), lty=1, lwd=2, col=colindx, bty="n", cex=1.2)


### 1.4.3. Cued specificity: barplots of pre and post infusion, S+ vs S- and VEH vs AP5
plot.new()
plot.window(xlim=c(0, 6), ylim=c(-2, 6))

#S+ both groups pre and post infusion
plotPrePostBars(data=byBinDataEarlyVEH[[5]], color=colindx[1], xmiddle=1, barwidth=0.5, colLabel = "white", labelY = 0.5)
plotPrePostBars(data=byBinDataEarlyAP5[[5]], color=colindx[2], xmiddle=2, barwidth=0.5, colLabel = "white", labelY = 0.5)

#S- both groups pre and post infusion
plotPrePostBars(data=byBinDataEarlyVEH[[6]], color=colindxB[1], xmiddle=3.25, barwidth=0.5, colLabel = "black", labelY = 0.5)
plotPrePostBars(data=byBinDataEarlyAP5[[6]], color=colindxB[2], xmiddle=4.25, barwidth=0.5, colLabel = "black", labelY = 0.5)

#Axis and labels
axis(side=1, tick = F, at=c(1.5, 3.75), labels=c("S+", "S-"), cex.axis=1.4, font=2)
axis(side=2, at=seq(-2, 6, 2), cex.axis=1, font=2, las=2)
mtext(side=2, line=2.5, text="Cued specificity", cex=1.4, font=2)
rect(xleft=3, xright=3.5, ybottom=5.5, ytop=6, col=colindx[1], border="white")
rect(xleft=3.5, xright=4, ybottom=5.5, ytop=6, col=colindxB[1], border="white")
rect(xleft=3, xright=3.5, ybottom=4.5, ytop=5, col=colindx[2], border="white")
rect(xleft=3.5, xright=4, ybottom=4.5, ytop=5, col=colindxB[2], border="white")
text(x=4.5, y=5.8, labels="VEH", cex=1.5)
text(x=4.5, y=4.8, labels="AP5", cex=1.5)







######################################################
######################################################
### STATISTICAL ANALYSES                           ###
######################################################
######################################################

Early_LongFormat #This is our object of reference

indexes <- unique(Early_LongFormat$Index)


## S+ Response ratio
DSRR <- subset(x=Early_LongFormat, Early_LongFormat$Index==indexes[1])
vehap5prepost.test <- ezANOVA(data=DSRR, dv=Performance, within=Infusion, between=Drug, wid=Rat, type=3)
# $`ANOVA`
#          Effect DFn DFd        F            p p<.05       ges
# 2          Drug   1   9 50.99187 5.419987e-05     * 0.6694572
# 3      Infusion   1   9 14.72451 3.983238e-03     * 0.5124852
# 4 Drug:Infusion   1   9 15.22492 3.608893e-03     * 0.5208306

#The interaction was significant. As post-hoc test, I'll split the dataset into the groups and, within each group, 
# use a paired t-test for the pre vs. post
DSRR_VEH <- subset(DSRR, Drug=="VEH")
DSRR_AP5 <- subset(DSRR, Drug=="AP5")

vehtest <- t.test(x=DSRR_VEH$Performance[DSRR_VEH$Infusion=="Pre"], 
                  y=DSRR_VEH$Performance[DSRR_VEH$Infusion=="Post"], 
                  paired=T, alternative="greater") #t(5)= -0.44473, p=0.66244712
ap5test <- t.test(x=DSRR_AP5$Performance[DSRR_AP5$Infusion=="Pre"], 
                  y=DSRR_AP5$Performance[DSRR_AP5$Infusion=="Post"], 
                  paired=T, alternative="greater") #t(4)=3.5043, p=0.02479956

p.adjust(p=c(vehtest$p.value, ap5test$p.value), method="holm") # 0.66244712 0.02479956


## S- Response ratio
NSRR <- subset(x=Early_LongFormat, Early_LongFormat$Index==indexes[2])
ezANOVA(data=NSRR, dv=Performance, within=Infusion, between=Drug, wid=Rat, type=3) #Nothing was significant
#         Effect DFn DFd         F          p p<.05        ges
#2          Drug   1   9 0.5167905 0.49045629       0.02197836
#3      Infusion   1   9 4.2795119 0.06850836       0.22445136
#4 Drug:Infusion   1   9 2.8361728 0.12645010       0.16093399


### DS RR AND NS RR
DSRR_NSRR <- rbind(DSRR, NSRR)
vehap5prepost.DSRRNSRR.test <- ezANOVA(data=DSRR_NSRR, dv=Performance, within=c(Infusion, Index), between=Drug, wid=Rat, type=3)
# $`ANOVA`
# Effect DFn DFd           F            p p<.05        ges
# 2                Drug   1   9  18.7354306 1.909441e-03     * 0.33331346
# 3            Infusion   1   9  15.7074285 3.288008e-03     * 0.35874208
# 5               Index   1   9 136.2991734 9.727380e-07     * 0.67663028
# 4       Drug:Infusion   1   9  13.7593887 4.849745e-03     * 0.32888322
# 6          Drug:Index   1   9  19.0464483 1.812275e-03     * 0.22624393
# 7      Infusion:Index   1   9   0.6816821 4.303494e-01       0.02229936
# 8 Drug:Infusion:Index   1   9   1.3698530 2.718946e-01       0.04382443



## S+ latency
DSlat <- subset(x=Early_LongFormat, Early_LongFormat$Index==indexes[3])
vehAP5.PrePost.test.Lat <- ezANOVA(data=DSlat, dv=Performance, within=Infusion, between=Drug, wid=Rat, type=3)

#The interaction was significant. As post-hoc test, I'll split the dataset into the groups and, within each group, use a paired t-test for the pre vs. post
DSlat_VEH <- subset(DSlat, Drug=="VEH")
DSlat_AP5 <- subset(DSlat, Drug=="AP5")

vehtest <- t.test(x=DSlat_VEH$Performance[DSlat_VEH$Infusion=="Pre"], 
                  y=DSlat_VEH$Performance[DSlat_VEH$Infusion=="Post"], 
                  paired=T, alternative="less") #t(5)= 0.70908, p=0.74502111
ap5test <- t.test(x=DSlat_AP5$Performance[DSlat_AP5$Infusion=="Pre"], 
                  y=DSlat_AP5$Performance[DSlat_AP5$Infusion=="Post"], 
                  paired=T, alternative="less") #t(4)=-3.0849, p=0.03675593

p.adjust(p=c(vehtest$p.value, ap5test$p.value), method="holm") #0.74502111 0.03675593

# $ANOVA
# Effect DFn DFd        F              p p<.05       ges
# 2          Drug   1   9 88.27426 0.000005998055     * 0.7005632
# 3      Infusion   1   9 11.00157 0.008985138980     * 0.4820834
# 4 Drug:Infusion   1   9 12.03800 0.007053418470     * 0.5045832


## S- latency
NSlat <- subset(x=Early_LongFormat, Early_LongFormat$Index==indexes[4])
ezANOVA(data=NSlat, dv=Performance, within=Infusion, between=Drug, wid=Rat, type=3)

#Nothing was significant:
# $ANOVA
         # Effect DFn DFd          F          p p<.05        ges
# 2          Drug   1   9 0.05383936 0.82170236       0.00198363
# 3      Infusion   1   9 4.40878841 0.06514878       0.24648147
# 4 Drug:Infusion   1   9 2.69415562 0.13513460       0.16659113

###DS AND NS LATENCY
DSNSlat <- rbind(DSlat, NSlat)
ezANOVA(data=DSNSlat, dv=Performance, within=c(Infusion, Index), between=Drug, wid=Rat, type=3)
# $`ANOVA`
#                Effect DFn DFd          F            p p<.05        ges
# 2                Drug   1   9  33.656896 2.590016e-04     * 0.42003301
# 3            Infusion   1   9  13.306906 5.335225e-03     * 0.38005858
# 5               Index   1   9 410.898873 8.064236e-09     * 0.79177201
# 4       Drug:Infusion   1   9  12.228779 6.755671e-03     * 0.36036261
# 6          Drug:Index   1   9  71.014208 1.457603e-05     * 0.39655808
# 7      Infusion:Index   1   9   2.051263 1.858816e-01       0.06567702
# 8 Drug:Infusion:Index   1   9   3.488967 9.461821e-02       0.10679325


## ITI latency
ITIlat <- subset(x=Early_LongFormat, Early_LongFormat$Index==indexes[7])
ITIlat.aov_test <- ezANOVA(data=ITIlat, dv=Performance, within=Infusion, between=Drug, wid=Rat, type=3)
# $ANOVA
#          Effect DFn DFd        F          p p<.05       ges
# 2          Drug   1   9 6.241261 0.03395910     * 0.3262542
# 3      Infusion   1   9 5.825659 0.03901543     * 0.1633909
# 4 Drug:Infusion   1   9 9.284387 0.01385936     * 0.2373706

#The interaction was significant. As post-hoc test, I'll split the dataset into the groups and, within each group, use a paired t-test for the pre vs. post
ITIlat_VEH <- subset(ITIlat, Drug=="VEH")
ITIlat_AP5 <- subset(ITIlat, Drug=="AP5")

vehtest <- t.test(x=ITIlat_VEH$Performance[ITIlat_VEH$Infusion=="Pre"], 
                  y=ITIlat_VEH$Performance[ITIlat_VEH$Infusion=="Post"], 
                  paired=T, alternative="less") #t(5)= -0.29979, p=0.38819996
ap5test <- t.test(x=ITIlat_AP5$Performance[ITIlat_AP5$Infusion=="Pre"], 
                  y=ITIlat_AP5$Performance[ITIlat_AP5$Infusion=="Post"], 
                  paired=T, alternative="less") #t(4)=-2.9156, p=0.04343353

p.adjust(p=c(vehtest$p.value, ap5test$p.value), method="holm") #0.73076380 0.04378754



## S+ specificity
DSspec <- subset(x=Early_LongFormat, Early_LongFormat$Index==indexes[5])
spec.test <- ezANOVA(data=DSspec, dv=Performance, within=Infusion, between=Drug, wid=Rat, type=3)
# $ANOVA
#          Effect DFn DFd         F          p p<.05       ges
# 2          Drug   1   9 12.857249 0.00587814     * 0.3688511
# 3      Infusion   1   9  4.437363 0.06443437       0.2256135
# 4 Drug:Infusion   1   9  4.277391 0.06856527       0.2192633

#The interaction was significant. As post-hoc test, I'll split the dataset into the groups and, within each group, use a paired t-test for the pre vs. post
DSspec_VEH <- subset(DSspec, Drug=="VEH")
DSspec_AP5 <- subset(DSspec, Drug=="AP5")

vehtest <- t.test(x=DSspec_VEH$Performance[DSspec_VEH$Infusion=="Pre"], y=DSspec_VEH$Performance[DSspec_VEH$Infusion=="Post"], paired=T, alternative="greater") #t = 0.061249, df = 5, p-value = 0.4768
ap5test <- t.test(x=DSspec_AP5$Performance[DSspec_AP5$Infusion=="Pre"], y=DSspec_AP5$Performance[DSspec_AP5$Infusion=="Post"], paired=T, alternative="greater") #t = 2.0081, df = 4, p-value = 0.04752
p.adjust(p=c(vehtest$p.value, ap5test$p.value), method="holm") #0.4767669 0.1150489


## S- specificity
NSspec <- subset(x=Early_LongFormat, Early_LongFormat$Index==indexes[6])
ezANOVA(data=NSspec, dv=Performance, within=Infusion, between=Drug, wid=Rat, type=3)
# $ANOVA
#          Effect DFn DFd         F          p p<.05        ges
# 2          Drug   1   9 3.6204787 0.08949286       0.18594300
# 3      Infusion   1   9 0.3845323 0.55056370       0.01813095
# 4 Drug:Infusion   1   9 0.5718089 0.46887081       0.02672518


###S+ and S- specificity together
DSNS.spec <- rbind(DSspec, NSspec)
ezANOVA(data=DSNS.spec, dv=Performance, within=c(Index, Infusion), between=Drug, wid=Rat, type=3)
 # $`ANOVA`
#                Effect DFn DFd           F            p p<.05        ges
# 2                Drug   1   9   0.8882808 3.705530e-01       0.03748976
# 3               Index   1   9 386.6251147 1.054812e-08     * 0.77438776
# 5            Infusion   1   9   8.1118638 1.914740e-02     * 0.11151720
# 4          Drug:Index   1   9  44.9048800 8.842506e-05     * 0.28502838
# 6       Drug:Infusion   1   9   2.9049243 1.225049e-01       0.04301426
# 7      Index:Infusion   1   9   1.2421203 2.939386e-01       0.05060416
# 8 Drug:Index:Infusion   1   9   3.0622127 1.140584e-01       0.11614290




### PERFORMANCE INDEX BY BIN
Early_LongFormatByBin$bin <- as.character(Early_LongFormatByBin$bin)

bins <- unique(Early_LongFormatByBin$bin)
#newBinsIndex <- c(1, 1, 1, 0, 2, 2, 2, 3, 3, 3, 4, 4) #I want to create 30min bins instead of 10min bins because to compare so many bins reduces my p values a lot when adjusting
newBinsIndex <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4) #I want to create 30min bins instead of 10min bins because to compare so many bins reduces my p values a lot when adjusting

newBinsVals <- sapply(seq(1, nrow(Early_LongFormatByBin)), function(l){
        sel <- as.numeric(Early_LongFormatByBin$bin[l])
        newBinsIndex[sel]
})

Early_LongFormatByBin$Bigbins <- newBinsVals

smallbins.aov <- summary(aov(perf ~ drug * bins + Error(rat/bins), data=Early_LongFormatByBin))
bigbins.aov <- summary(aov(perf ~ drug * Bigbins + Error(rat/(Bigbins)), data=Early_LongFormatByBin))

#Including kind of cue (DS vs NS) as a within factor too
Early_LongFormatByBin_DSandNS$bin <- as.character(Early_LongFormatByBin_DSandNS$bin)
bins <- unique(Early_LongFormatByBin_DSandNS$bin)
newBinsIndex <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4) #I want to create 30min bins instead of 10min bins because to compare so many bins reduces my p values a lot when adjusting
newBinsVals <- sapply(seq(1, nrow(Early_LongFormatByBin_DSandNS)), function(l){
        sel <- as.numeric(Early_LongFormatByBin_DSandNS$bin[l])
        newBinsIndex[sel]
})

Early_LongFormatByBin_DSandNS$Bigbins <- newBinsVals

ezANOVA(data=Early_LongFormatByBin_DSandNS, dv=perf, within=c(index, Bigbins), between=drug, wid=rat, type=3)


### Results of the Mixed-effects (1 within, 1 btwn-subject factor) ANOVA with the original 10min bins
# Error: rat
#            Df Sum Sq Mean Sq F value  Pr(>F)   
# drug       1  232.0  232.01   12.12 0.00692 **

#         Residuals  9  172.3   19.14                   
# ---
#         Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Error: rat:bin
#           Df Sum Sq Mean Sq F value Pr(>F)
# bin       11  48.71   4.429   1.433  0.170
# drug:bin  11  54.56   4.960   1.605  0.109
# Residuals 99 306.02   3.091 



### Results of the Mixed-effects (1 within, 1 btwn-subject factor) ANOVA with the BIG BINS (30 min)
# Error: rat
#           Df Sum Sq Mean Sq F value  Pr(>F)   
# drug       1  232.0  232.01   12.12 0.00692 **
#         Residuals  9  172.3   19.14                   
# ---
#         Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Error: rat:Bigbins
#              Df Sum Sq Mean Sq F value Pr(>F)  
# Bigbins       1   9.40    9.40   1.511 0.2502  
# drug:Bigbins  1  31.79   31.79   5.111 0.0501 .
# Residuals     9  55.98    6.22                 
# ---
#         Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Error: Within
#            Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 110  312.1   2.837 

#Just to double-check
ezANOVA(data=Early_LongFormatByBin, dv=perf, within=Bigbins, between=drug, wid=rat, type=1)
# $`ANOVA`
#         Effect DFn DFd         F           p p<.05        ges
# 1         drug   1   9 12.119157 0.006924822     * 0.50405442
# 2      Bigbins   1   9  1.510524 0.250223403       0.03953167
# 3 drug:Bigbins   1   9  5.111312 0.050109313       0.12224727

ttestPerBin <- do.call("rbind", lapply(seq(1, length(unique(Early_LongFormatByBin$bin))), function(m){
       bindex <- unique(Early_LongFormatByBin$bin)[m]
        tst <- t.test(x=EarlyVEH_DStaskAccByBin_LongFormat$perf[EarlyVEH_DStaskAccByBin_LongFormat$bin==bindex], 
                      y=EarlyAP5_DStaskAccByBin_LongFormat$perf[EarlyAP5_DStaskAccByBin_LongFormat$bin==bindex], paired=F, alternative="greater") 
        
        data.frame(bin=bindex, t=tst$statistic, df=tst$parameter, p=tst$p.value)
})
)

#Adjust the t test p values and also the p values of the ANOVA (using the small bins)
padjusted <- p.adjust(p=c(0.00692, 0.17, 0.109, ttestPerBin$p), method="holm")


ttestPerBin$p.adjusted <- padjusted[-c(1:3)]

# bin           t       df               p p.adjusted
# t     1  1.83145811 5.304251 0.061588185 0.23241359
# t1    2 -0.09012978 8.012100 0.534801740 0.69922324
# t2    3  0.40517671 6.078269 0.349611621 0.69922324
# t3    4  2.12717518 6.890756 0.035793356 0.18215673
# t4    5  2.92773629 7.727084 0.009899460 0.09899460
# t5    6  2.65123274 8.242288 0.014228686 0.11500863
# t6    7  3.59036615 8.723828 0.003069769 0.03683723
# t7    8  2.79555498 4.178221 0.023336776 0.16335743
# t8    9  1.83100539 6.076459 0.058103398 0.23241359
# t9   10  2.82825212 6.977097 0.012778737 0.11500863
# t10  11  3.00188759 8.884798 0.007562351 0.08318586
# t11  12  2.21600537 7.294561 0.030359454 0.18215673

newBinsVals <- sapply(seq(1, nrow(EarlyVEH_DStaskAccByBin_LongFormat)), function(l){
        sel <- EarlyVEH_DStaskAccByBin_LongFormat$bin[l]
        newBinsIndex[sel]
})

EarlyVEH_DStaskAccByBin_LongFormat$Bigbins <- newBinsVals

newBinsVals <- sapply(seq(1, nrow(EarlyAP5_DStaskAccByBin_LongFormat)), function(l){
        sel <- EarlyAP5_DStaskAccByBin_LongFormat$bin[l]
        newBinsIndex[sel]
})

EarlyAP5_DStaskAccByBin_LongFormat$Bigbins <- newBinsVals


ttestPerBigBin <- do.call("rbind", lapply(seq(1, length(unique(Early_LongFormatByBin$Bigbins))), function(m){
        bindex <- unique(Early_LongFormatByBin$Bigbins)[m]
        tst <- t.test(x=EarlyVEH_DStaskAccByBin_LongFormat$perf[EarlyVEH_DStaskAccByBin_LongFormat$Bigbins==bindex], 
                      y=EarlyAP5_DStaskAccByBin_LongFormat$perf[EarlyAP5_DStaskAccByBin_LongFormat$Bigbins==bindex], paired=F, alternative="greater") 
        
        data.frame(Bigbins=bindex, t=tst$statistic, df=tst$parameter, p=tst$p.value)
})
)

#Adjust the t test p values (using the big bins)
padjusted <- p.adjust(p=c(ttestPerBigBin$p), method="holm")


ttestPerBigBin$p.adjusted <- p.adjust(p=ttestPerBigBin$p, method="holm")

###Making 31-40 its own separate bin
#    Bigbins         t        df              p    p.adjusted
# t        1 0.8568663 23.645328 0.200060020152 0.34000000000 #1-30 min. I'll use this as the PRE window
# t1       0 2.1271752  6.890756 0.035793355919 0.14317342368 #31-40 min. I discard this because it's the time at which the infusion is taking place
# t2       2 5.4051789 28.468158 0.000004356194 0.00003484955 #41-70 min. I'll use this as the POST window
# t3       3 4.4557853 21.256654 0.000106486237 0.00074540366 #71-100 min
# t4       4 3.7628519 18.653759 0.000676119814 0.00405671888 #101-120 min


### including 31-40min
# Bigbins         t       df            p   p.adjusted
# t        1 0.3039309 22.70217 0.3819743836 0.3819743836
# t1       2 4.2915647 25.93471 0.0001093981 0.0004375925
# t2       3 4.0209585 25.39606 0.0002292323 0.0006876970
# t3       4 3.5534131 27.75646 0.0006919061 0.0013838123






#################################################################
####### NEURONAL FIRING                                ##########
#################################################################

###############################################################
#1. PSTH pre and post infusion
###############################################################

######
#VEH #
######

# S+ Onset
psthInf(formatDat="Zscores", group="VEH", event="S+", comp=c("Pre VEH injection", "Post VEH injection"), 
        expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, 
        col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=2, binw=50, neudata=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf), 
        stimulus="cue", imgFormat="pdf", BLNeuData=0)  

psthInf(formatDat="raw", group="VEH", event="S+", comp=c("Pre VEH injection", "Post VEH injection"), 
        expName = "Early", errShade=T, ymax=26, graphFolder=neuGraphFolder, 
        col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf), 
        stimulus="cue", imgFormat="pdf", BLNeuData=0)       

#S- Onset
psthInf(formatDat="Zscores", group="VEH", event="S-", comp=c("Pre VEH injection", "Post VEH injection"), 
        expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), 
        infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsNSEarlyVEHPreInf, allNeuronsNSEarlyVEHPostInf),
        stimulus="cue", imgFormat="pdf", BLNeuData=0)       

psthInf(formatDat="raw", group="VEH", event="S-", comp=c("Pre VEH injection", "Post VEH injection"), 
        expName = "Early", errShade=T, ymax=26, graphFolder=neuGraphFolder, col=c("black", colindx[1]), 
        infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsNSEarlyVEHPreInf, allNeuronsNSEarlyVEHPostInf), 
        stimulus="cue", imgFormat="pdf", BLNeuData=0)       


#S+ Entry
psthInf(formatDat="Zscores", group="VEH", event="S+ Entry", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryDSEarlyVEHPreInf, allNeuronsEntryDSEarlyVEHPostInf), stimulus="entry",  BLNeuData=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf), imgFormat="pdf")       
psthInf(formatDat="raw", group="VEH", event="S+ Entry", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryDSEarlyVEHPreInf, allNeuronsEntryDSEarlyVEHPostInf), stimulus="entry", BLNeuData=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf), imgFormat="pdf")       

#S- Entry
psthInf(formatDat="Zscores", group="VEH", event="S- Entry", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryNSEarlyVEHPreInf, allNeuronsEntryNSEarlyVEHPostInf), stimulus="entry", BLNeuData=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf), imgFormat="pdf")       
psthInf(formatDat="raw", group="VEH", event="S- Entry", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryNSEarlyVEHPreInf, allNeuronsEntryNSEarlyVEHPostInf), stimulus="entry", BLNeuData=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf), imgFormat="pdf")       

#ITI Entry
psthInf(formatDat="Zscores", group="VEH", event="ITI Entry", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryITIEarlyVEHPreInf, allNeuronsEntryITIEarlyVEHPostInf), stimulus="entry", BLNeuData=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf), imgFormat="pdf")       
psthInf(formatDat="raw", group="VEH", event="ITI Entry", comp=c("Pre VEH injection", "Post VEH injection"), expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[1]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryITIEarlyVEHPreInf, allNeuronsEntryITIEarlyVEHPostInf), stimulus="entry", BLNeuData=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf), imgFormat="pdf")       


########
#AP5
########

#S+ Onset
psthInf(formatDat="Zscores", group="AP5", event="S+", comp=c("Pre AP5 injection", "Post AP5 injection"), 
        expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), 
        infTime=1800, infDur=12*60, xmin=0.5, xmax=2, binw=50, 
        neudata=list(allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf), stimulus="cue", imgFormat="pdf", BLNeuData=0)

psthInf(formatDat="raw", group="AP5", event="S+", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Early", errShade=T, ymax=26, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf), stimulus="cue", imgFormat="pdf", BLNeuData=0)       

#S- Onset
psthInf(formatDat="Zscores", group="AP5", event="S-", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsNSEarlyAP5PreInf, allNeuronsNSEarlyAP5PostInf), stimulus="cue", imgFormat="pdf", BLNeuData=0)       
psthInf(formatDat="raw", group="AP5", event="S-", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Early", errShade=T, ymax=26, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsNSEarlyAP5PreInf, allNeuronsNSEarlyAP5PostInf), stimulus="cue", imgFormat="pdf", BLNeuData=0)       


#S+ Entry
psthInf(formatDat="Zscores", group="AP5", event="S+ Entry", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryDSEarlyAP5PreInf, allNeuronsEntryDSEarlyAP5PostInf), stimulus="entry", BLNeuData=list(allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf), imgFormat="pdf")       
psthInf(formatDat="raw", group="AP5", event="S+ Entry", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryDSEarlyAP5PreInf, allNeuronsEntryDSEarlyAP5PostInf), stimulus="entry", BLNeuData=list(allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf), imgFormat="pdf")       


#S- Entry
psthInf(formatDat="Zscores", group="AP5", event="S- Entry", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryNSEarlyAP5PreInf, allNeuronsEntryNSEarlyAP5PostInf), stimulus="entry", BLNeuData=list(allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf), imgFormat="pdf")       
psthInf(formatDat="raw", group="AP5", event="S- Entry", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryNSEarlyAP5PreInf, allNeuronsEntryNSEarlyAP5PostInf), stimulus="entry", BLNeuData=list(allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf), imgFormat="pdf")       


#ITI Entry
psthInf(formatDat="Zscores", group="AP5", event="ITI Entry", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryITIEarlyAP5PreInf, allNeuronsEntryITIEarlyAP5PostInf), stimulus="entry", BLNeuData=list(allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf), imgFormat="pdf")       
psthInf(formatDat="raw", group="AP5", event="ITI Entry", comp=c("Pre AP5 injection", "Post AP5 injection"), expName = "Early", errShade=T, ymax=14, graphFolder=neuGraphFolder, col=c("black", colindx[2]), infTime=1800, infDur=12*60, 
        xmin=0.5, xmax=1.5, binw=50, neudata=list(allNeuronsEntryITIEarlyAP5PreInf, allNeuronsEntryITIEarlyAP5PostInf), stimulus="entry", BLNeuData=list(allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf), imgFormat="pdf")       




###############################################################
#2. POINTS pre and post infusion around time of cue
###############################################################

dotplot(neudata=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf, allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf),
        expName="Early", dot="Medians", Lines=T, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400)
dotplot(neudata=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf, allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf),
        expName="Early", dot="Means", Lines=T, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400)
dotplot(neudata=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf, allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf),
        expName="Early", dot="Medians", Lines=F, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400)
dotplot(neudata=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf, allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf),
        expName="Early", dot="Means", Lines=F, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400)

#Same but with boxplot instead of dotplot

#100-400ms
dotPlotByGroupDS <- dotplot(boxplot=T, neudata=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf, allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf),
        expName="Early 100-400", Lines=T, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400, ytop=12)
dotPlotByGroupNS <- dotplot(boxplot=T, neudata=list(allNeuronsNSEarlyVEHPreInf, allNeuronsNSEarlyVEHPostInf, allNeuronsNSEarlyAP5PreInf, allNeuronsNSEarlyAP5PostInf),
                            expName="Early 100-400", Lines=T, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400, ytop=12)
save(dotPlotByGroupDS, file=paste(dataForRdir, "dotPlotByGroupDS.rdat", sep=""))
save(dotPlotByGroupNS, file=paste(dataForRdir, "dotPlotByGroupNS.rdat", sep=""))


#750-2000
dotPlotByGroupDS_TAIL <- dotplot(boxplot=T, neudata=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf, allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf),
        expName="Early 750-2000", Lines=T, col=colindx, plotWidth=0.3, event="S+", winmin=750, winmax=2000, ytop=12)

dotPlotByGroupNS_TAIL <- dotplot(boxplot=T, neudata=list(allNeuronsNSEarlyVEHPreInf, allNeuronsNSEarlyVEHPostInf, allNeuronsNSEarlyAP5PreInf, allNeuronsNSEarlyAP5PostInf),
        expName="Early 750-2000", Lines=T, ytop=12, ybottom=-2, col=colindx, plotWidth=0.3, event="S-", winmin=750, winmax=2000)

save(dotPlotByGroupDS_TAIL, file=paste(dataForRdir, "dotPlotByGroupDS_TAIL.rdat", sep=""))
save(dotPlotByGroupNS_TAIL, file=paste(dataForRdir, "dotPlotByGroupNS_TAIL.rdat", sep=""))



#Raw
#100-400ms after the S+
# dotplot(boxplot=T, neudata=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf, allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf),
#         expName="Early 100-400", Lines=T, formatDat="Raw", ytop=20, ybottom=0, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400, comp=c("VEH", "AP5"))
# dotplot(boxplot=T, neudata=list(allNeuronsNSEarlyVEHPreInf, allNeuronsNSEarlyVEHPostInf, allNeuronsNSEarlyAP5PreInf, allNeuronsNSEarlyAP5PostInf),
#         expName="Early 100-400", Lines=T, formatDat="Raw", ytop=20, ybottom=0, col=colindx, plotWidth=0.3, event="S-", winmin=100, winmax=400, comp=c("VEH", "AP5"))
# dotplot(boxplot=T, neudata=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf, allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf),
#         expName="Early 100-400", Lines=F, formatDat="Raw", ytop=20, ybottom=0, col=colindx, plotWidth=0.3, event="S+", winmin=100, winmax=400, comp=c("VEH", "AP5"))
# dotplot(boxplot=T, neudata=list(allNeuronsNSEarlyVEHPreInf, allNeuronsNSEarlyVEHPostInf, allNeuronsNSEarlyAP5PreInf, allNeuronsNSEarlyAP5PostInf),
#         expName="Early 100-400", Lines=F, formatDat="Raw", ytop=20, ybottom=0, col=colindx, plotWidth=0.3, event="S-", winmin=100, winmax=400, comp=c("VEH", "AP5"))



### STATISTICAL TEST
#Comparison of cue-evoked firing rate (100ms-400ms window) pre vs post injection for both groups. I'll use a Wilcoxon paired test 
# separately for each group 

load(file=paste(dataForRdir, "dotPlotByGroup.rdat", sep=""))
dotPlotByGroupEarly <- dotPlotByGroupDS
dotplotDataVEH_Early <- dotPlotByGroupEarly$VEH
dotplotDataAP5_Early <- dotPlotByGroupEarly$AP5

dotPlotByGroupEarlyNS <- dotPlotByGroupNS
dotplotDataVEH_Early_NS <- dotPlotByGroupEarlyNS$VEH
dotplotDataAP5_Early_NS <- dotPlotByGroupEarlyNS$AP5

#100-400ms after S+
# DS_VEH_Early pre vs. post
wilcox.test(x=dotplotDataVEH_Early[,1], y=dotplotDataVEH_Early[,2], paired=T, alternative = "greater") #V = 22, p=0.3203; pcorrected=6.4060e-01 ; #Raw scores: V=24, p=0.4609
# DS_AP5_Early pre vs. post
wilcox.test(x=dotplotDataAP5_Early[,1], y=dotplotDataAP5_Early[,2], paired=T, alternative = "greater") #V = 424, p=9.122e-06; pcorrected= 5.4732e-05; Raw scores: V=80, p=0.0011
# NS_VEH_Early pre vs. post
wilcox.test(x=dotplotDataVEH_Early_NS[,1], y=dotplotDataVEH_Early_NS[,2], paired=T, alternative = "greater") #V = 22, p=8.2020e-02; pcorrected=1
# NS_AP5_Early pre vs. post
wilcox.test(x=dotplotDataAP5_Early_NS[,1], y=dotplotDataAP5_Early_NS[,2], paired=T, alternative = "greater") #V = 89, p=9.9890e-01; pcorrected= 1.0944e-04
# DS vs NS: Early VEH PRE
wilcox.test(x=dotplotDataVEH_Early[,1], y=dotplotDataVEH_Early_NS[,1], paired=T, alternative = "greater") #V = 36, p=1.9530e-02; pcorrected=3.1200e-02
# DS vs. NS: Early AP5 PRE
wilcox.test(x=dotplotDataAP5_Early[,1], y=dotplotDataAP5_Early_NS[,1], paired=T, alternative = "greater") #V = 442, p=4.7680e-06 ; pcorrected=9.5360e-06
# DS vs. NS: Early VEH POST
wilcox.test(x=dotplotDataVEH_Early[,2], y=dotplotDataVEH_Early_NS[,2], paired=T, alternative = "greater") #V = 35, p=3.1252e-02; pcorrected=3.1689e-02 
# DS vs. NS: Early AP5 POST
wilcox.test(x=dotplotDataAP5_Early[,2], y=dotplotDataAP5_Early_NS[,2], paired=T, alternative = "greater") #V = 430, p=2.7979e-05; pcorrected=5.5958e-05

p.adjust(p=c(0.3203, 9.122e-06, 0.02734, 0.9989, 0.003906, 5.96e-07, 0.007813, 3.997e-06))
#  6.4060e-01 5.4732e-05 8.2020e-02 9.9890e-01 1.9530e-02 4.7680e-06 3.1252e-02 2.7979e-05


#750-2000ms after S+
load(file=paste(dataForRdir, "dotPlotByGroupDS_TAIL.rdat", sep=""))
load(file=paste(dataForRdir, "dotPlotByGroupNS_TAIL.rdat", sep=""))

# Create objects that will be useful for the comparisons
dotPlotByGroupEarly <- dotPlotByGroupDS_TAIL
dotplotDataVEH_Early <- dotPlotByGroupEarly$VEH
dotplotDataAP5_Early <- dotPlotByGroupEarly$AP5

dotPlotByGroupEarlyNS <- dotPlotByGroupNS_TAIL
dotplotDataVEH_Early_NS <- dotPlotByGroupEarlyNS$VEH
dotplotDataAP5_Early_NS <- dotPlotByGroupEarlyNS$AP5

#750-2000ms after S+
# DS_VEH_Early pre vs. post
wilcox.test(x=dotplotDataVEH_Early[,1], y=dotplotDataVEH_Early[,2], paired=T, alternative = "greater") #V = 6, p=0.9609; pcorrected=6.4060e-01 ; #Raw scores: V=24, p=0.4609
# DS_AP5_Early pre vs. post
wilcox.test(x=dotplotDataAP5_Early[,1], y=dotplotDataAP5_Early[,2], paired=T, alternative = "greater") #V = 150, p=0.956; pcorrected= 5.4732e-05; Raw scores: V=80, p=0.0011
# NS_VEH_Early pre vs. post
wilcox.test(x=dotplotDataVEH_Early_NS[,1], y=dotplotDataVEH_Early_NS[,2], paired=T, alternative = "greater") #V = 22, p=0.3203; pcorrected=1
# NS_AP5_Early pre vs. post
wilcox.test(x=dotplotDataAP5_Early_NS[,1], y=dotplotDataAP5_Early_NS[,2], paired=T, alternative = "greater") #V = 162, p=0.927; pcorrected= 1.0944e-04
# DS vs NS: Early VEH PRE
wilcox.test(x=dotplotDataVEH_Early[,1], y=dotplotDataVEH_Early_NS[,1], paired=T, alternative = "greater") #V = 31, p=0.03906; pcorrected=3.1200e-02
# DS vs. NS: Early AP5 PRE
wilcox.test(x=dotplotDataAP5_Early[,1], y=dotplotDataAP5_Early_NS[,1], paired=T, alternative = "greater") #V = 260, p=0.2919 ; pcorrected=9.5360e-06
# DS vs. NS: Early VEH POST
wilcox.test(x=dotplotDataVEH_Early[,2], y=dotplotDataVEH_Early_NS[,2], paired=T, alternative = "greater") #V = 30, p=0.05469; pcorrected=3.1689e-02 
# DS vs. NS: Early AP5 POST
wilcox.test(x=dotplotDataAP5_Early[,2], y=dotplotDataAP5_Early_NS[,2], paired=T, alternative = "greater") #V = 292, p=0.1143; pcorrected=5.5958e-05







#Correct p values taking into account the other 2 wilcoxon tests from the "Late" test
p.adjust(p=c(0.9609, 0.956, 0.3203, 0.927, 0.03906, 0.2919, 0.05469, 0.1143), method="holm")
# 1.00000 1.00000 1.00000 1.00000 0.31248 1.00000 0.38283 0.68580



###############################################################
#3. SCATTERPLOT pre and post infusion around time of cue
###############################################################

prePostInf_FR(data=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf), 
             dataformat="Raw", BLwdw=2, winmin=100, winmax=400, col_labels="purple",
             comparison="Early VEH Pre vs. Post", graphfolder=preVsPostFRFolder, 
             xmin=0, ymin=0)
        
prePostInf_FR(data=list(allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf), 
              dataformat="Raw", BLwdw=2, winmin=100, winmax=400, col_labels="purple",
              comparison="Early AP5 Pre vs. Post", graphfolder=preVsPostFRFolder)


#Baseline
prePostInf_FR(data=list(allNeuronsDSEarlyVEHPreInf, allNeuronsDSEarlyVEHPostInf), 
              dataformat="Raw", BLwdw=2, winmin=-2000, winmax=0,
              comparison="Baseline_Early VEH Pre vs. Post", graphfolder=preVsPostFRFolder)

prePostInf_FR(data=list(allNeuronsDSEarlyAP5PreInf, allNeuronsDSEarlyAP5PostInf), 
              dataformat="Raw", BLwdw=2, winmin=-2000, winmax=0,
              comparison="Baseline_Early AP5 Pre vs. Post", graphfolder=preVsPostFRFolder)






################################################
### EXCITATION AND INHIBITION BY BIN
################################################
postInfTargetWdw <- 1800+12*60+30*60 #For the post infusion window, I'll choose the period between the end of the infusion +30'.

#Matrix in which rows are 50ms bins after the cue, columns are individual neurons and the values indicate if the neuron was EXCITED (ExcBins) or INHIBITED (InhBins) on that bin
NEXfiles <- "E:/Dropbox/NMDA/EXP1_Performance/Early VEH/NEX files/"
EarlyVEHPreInf_ExcBins <- KC.sigbins(path=NEXfiles, startt=0, endt=1800, event=1, BLwdw=5, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)
EarlyVEHPostInf_ExcBins <- KC.sigbins(path=NEXfiles, startt=2520, endt=postInfTargetWdw, event=1, BLwdw=5, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)
EarlyVEHPreInf_InhBins <- KC.inhib.sigbins(path=NEXfiles, startt=0, endt=1800, event=1, BLwdw=5, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)
EarlyVEHPostInf_InhBins <- KC.inhib.sigbins(path=NEXfiles, startt=2520, endt=postInfTargetWdw, event=1, BLwdw=5, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)

# In neuralhist, I flagged neurons as CUE-EXCITED if they were excited (>99.9% confidence interval of a Poisson distribution given by BL firing) for 3 consecutive 10ms bins in the 500ms window after the cue. I used the 2s precue window as baseline to define my Poisson distribution.
#EarlyVEH_ExcUnits <- unlist(allNeuronsDSEarlyVEHPreInf$cueexidx) #Index of cue-excited units

#Redefine NEXfiles now so that it sends the function to the AP5 files and repeat
NEXfiles <- "E:/Dropbox/NMDA/EXP1_Performance/Early AP5/NEX files/"
EarlyAP5PreInf_ExcBins <- KC.sigbins(path=NEXfiles, startt=0, endt=1800, event=1, BLwdw=5, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)
EarlyAP5PostInf_ExcBins <- KC.sigbins(path=NEXfiles, startt=2520, endt=postInfTargetWdw, event=1, BLwdw=5, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)
EarlyAP5PreInf_InhBins <- KC.inhib.sigbins(path=NEXfiles, startt=0, endt=1800, event=1, BLwdw=5, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)
EarlyAP5PostInf_InhBins <- KC.inhib.sigbins(path=NEXfiles, startt=2520, endt=postInfTargetWdw, event=1, BLwdw=5, PostEvent_wdw=1, pbin=0.05, funcdirect=funcdirect)

#Save these files
save(EarlyVEHPreInf_ExcBins, file=paste(dataForRdir, "EarlyVEHPreInf_ExcBins.rdat", sep=""))
save(EarlyVEHPostInf_ExcBins, file=paste(dataForRdir, "EarlyVEHPostInf_ExcBins.rdat", sep=""))
save(EarlyVEHPreInf_InhBins, file=paste(dataForRdir, "EarlyVEHPreInf_InhBins.rdat", sep=""))
save(EarlyVEHPostInf_InhBins, file=paste(dataForRdir, "EarlyVEHPostInf_InhBins.rdat", sep=""))

save(EarlyAP5PreInf_ExcBins, file=paste(dataForRdir, "EarlyAP5PreInf_ExcBins.rdat", sep=""))
save(EarlyAP5PostInf_ExcBins, file=paste(dataForRdir, "EarlyAP5PostInf_ExcBins.rdat", sep=""))
save(EarlyAP5PreInf_InhBins, file=paste(dataForRdir, "EarlyAP5PreInf_InhBins.rdat", sep=""))
save(EarlyAP5PostInf_InhBins, file=paste(dataForRdir, "EarlyAP5PostInf_InhBins.rdat", sep=""))


### DEFINE WHAT UNITS QUALIFY AS CUE-EXCITED BASED ON THE BINS THEY'RE SIGNIFICANTLY CUE-EXCITED (criterion: 3 consecutive 50ms bins)

###############################################
### PROPORTION OF CUE EXCITED NEURONS
###############################################
# % of CUE-EXCITED UNITS 
#This function tells me, based on the "ExcBins" matrix, which units qualify as cue-excited based on my criterion
CueExcIndex <- function(excbybin, threhold=3){
        
        sapply(seq(1, ncol(excbybin)), function(x){
                exc <- as.numeric(excbybin[1:10, x])
                oneruns <- rle(exc)$lengths[rle(exc)$values==1]
                cueexc <- FALSE
                if(length(oneruns)>0 & sum(oneruns)>=threshold){cueexc <- TRUE}
                cueexc
        })
}

EarlyVEHPreInf_ExcUnits <-  CueExcIndex(EarlyVEHPreInf_ExcBins) 
EarlyVEHPostInf_ExcUnits <- CueExcIndex(EarlyVEHPostInf_ExcBins)
EarlyAP5PreInf_ExcUnits <- CueExcIndex(EarlyAP5PreInf_ExcBins)
EarlyAP5PostInf_ExcUnits <- CueExcIndex(EarlyAP5PostInf_ExcBins)


contTable_EarlyVEH <- t(data.frame(Pre=as.matrix(table(EarlyVEHPreInf_ExcUnits)), Post=as.matrix(table(EarlyVEHPostInf_ExcUnits))))
contTable_EarlyAP5 <- t(data.frame(Pre=as.matrix(table(EarlyAP5PreInf_ExcUnits)), Post=as.matrix(table(EarlyAP5PostInf_ExcUnits))))

chisq.test(contTable_EarlyVEH) #X-squared = 0, df = 1, p-value = 1
chisq.test(contTable_EarlyAP5) #X-squared = 0.28202, df = 1, p-value = 0.5954

fisher.test(contTable_EarlyVEH) #CI: 0.1557766 18.7666075, odds ratio= 1.613716,  p-value = 1
fisher.test(contTable_EarlyAP5) #CI:  0.2007157 2.1061495, odds ratio=0.6585366,  p-value = 0.5959




#Early VEH
plot.new()
par(mar=c(2, 6, 2, 2))
plot.window(xlim=c(0, 2), ylim=c(0, 1))
rect(xleft=0, xright=1, ybottom=0, ytop=sum(EarlyVEHPreInf_ExcUnits)/length(EarlyVEHPreInf_ExcUnits), col="gray30", border = F)
rect(xleft=1, xright=2, ybottom=0, ytop=sum(EarlyVEHPostInf_ExcUnits)/length(EarlyVEHPostInf_ExcUnits), col=colindx[1], border= F)
axis(side=1, at=c(0.5, 1.5), tick = F, labels=c("Pre", "Post"), cex.axis=1.5, font=2)
axis(side=2, at=seq(0, 1, 0.25), cex.axis=1.4, las=2)
mtext(side=2, line=4, text="Proportion", cex=1.5, font=2)


#Early AP5
plot.new()
par(mar=c(2, 6, 2, 2))
plot.window(xlim=c(0, 2), ylim=c(0, 1))
rect(xleft=0, xright=1, ybottom=0, ytop=sum(EarlyAP5PreInf_ExcUnits)/length(EarlyAP5PreInf_ExcUnits), col="gray30", border = F)
rect(xleft=1, xright=2, ybottom=0, ytop=sum(EarlyAP5PostInf_ExcUnits)/length(EarlyAP5PostInf_ExcUnits), col=colindx[2], border= F)
axis(side=1, at=c(0.5, 1.5), tick = F, labels=c("Pre", "Post"), cex.axis=1.5, font=2)
axis(side=2, at=seq(0, 1, 0.25), cex.axis=1.4, las=2)
mtext(side=2, line=4, text="Proportion", cex=1.5, font=2)



############################################################################
#Plot % bins excited/inhibited before and after infusion of VEH or AP5
############################################################################

#Function to calculate the percentage of units exc/inh to apply on the objects that I created with KC.sigbins.R and KC.inhib.sigbins.R
PercBins <- function(sigBinData){
        sapply(seq(1, nrow(sigBinData)), function (x){
                sum(sigBinData[x,])/ncol(sigBinData)
        })
}


#Early VEH All units
plot.new()
plot.window(xlim = c(0, nrow(EarlyVEHPreInf_ExcBins)), ylim=c(0, 1))
abline(h=seq(-1, 1, by=0.25), col="gray90")

lines(x=seq(1, nrow(EarlyVEHPreInf_ExcBins)), y=PercBins(EarlyVEHPreInf_ExcBins), col="gray30", lwd=2)
lines(x=seq(1, nrow(EarlyVEHPostInf_ExcBins)), y=PercBins(EarlyVEHPostInf_ExcBins), col="blue", lwd=2)

lines(x=seq(1, nrow(EarlyVEHPreInf_InhBins)), y=-PercBins(EarlyVEHPreInf_InhBins), col="gray30", lwd=2)
lines(x=seq(1, nrow(EarlyVEHPostInf_InhBins)), y=-PercBins(EarlyVEHPostInf_InhBins), col="blue", lwd=2)

axis(side=1, at=seq(0, nrow(EarlyVEHPreInf_ExcBins), by=10), labels=seq(0, 1, by=0.5), cex.axis=1.4)
axis(side=2, las=2, at=seq(0, 1, by=0.5), labels=seq(0, 100, 50), cex.axis=1.4)
axis(side=2, las=2, at=seq(0, -1, by=-0.5), labels=seq(0, 100, 50), cex.axis=1.4)
mtext(side=1, text="Time from S+ onset (s)", font=2, cex=1.5, line=2.5)
mtext(side=2, text="% Excited", at=0.5, font=2, cex = 1.5, line=2.5)
mtext(side=2, text="% Inhibited", at=-0.5, font=2, cex = 1.5, line=2.5)

#Early VEH Cue-excited units only
plot.new()
plot.window(xlim = c(0, nrow(EarlyVEHPreInf_ExcBins[,EarlyVEHPreInf_ExcUnits])), ylim=c(0, 1))
abline(h=seq(-1, 1, by=0.25), col="gray90")

lines(x=seq(1, nrow(EarlyVEHPreInf_ExcBins)), y=PercBins(EarlyVEHPreInf_ExcBins[,EarlyVEHPreInf_ExcUnits]), col="gray30", lwd=2)
lines(x=seq(1, nrow(EarlyVEHPostInf_ExcBins)), y=PercBins(EarlyVEHPostInf_ExcBins[,EarlyVEHPreInf_ExcUnits]), col="blue", lwd=2)

axis(side=1, at=seq(0, nrow(EarlyVEHPreInf_ExcBins), by=10), labels=seq(0, 1, by=0.5), cex.axis=1.4)
axis(side=2, las=2, at=seq(0, 1, by=0.5), labels=seq(0, 100, 50), cex.axis=1.4)
axis(side=2, las=2, at=seq(0, -1, by=-0.5), labels=seq(0, 100, 50), cex.axis=1.4)
mtext(side=1, text="Time from S+ onset (s)", font=2, cex=1.5, line=2.5)
mtext(side=2, text="% Excited", at=0.5, font=2, cex = 1.5, line=2.5)



#Early AP5 All units
plot.new()
plot.window(xlim = c(0, nrow(EarlyAP5PreInf_ExcBins)), ylim=c(0, 1))
abline(h=seq(-1, 1, by=0.25), col="gray90")

lines(x=seq(1, nrow(EarlyAP5PreInf_ExcBins)), y=PercBins(EarlyAP5PreInf_ExcBins), col="gray30", lwd=2)
lines(x=seq(1, nrow(EarlyAP5PostInf_ExcBins)), y=PercBins(EarlyAP5PostInf_ExcBins), col="red", lwd=2)

lines(x=seq(1, nrow(EarlyAP5PreInf_InhBins)), y=-PercBins(EarlyAP5PreInf_InhBins), col="gray30", lwd=2)
lines(x=seq(1, nrow(EarlyAP5PostInf_InhBins)), y=-PercBins(EarlyAP5PostInf_InhBins), col="red", lwd=2)

axis(side=1, at=seq(0, nrow(EarlyVEHPreInf_ExcBins), by=10), labels=seq(0, 1, by=0.5), cex.axis=1.4)
axis(side=2, las=2, at=seq(0, 1, by=0.5), labels=seq(0, 100, 50), cex.axis=1.4)
axis(side=2, las=2, at=seq(0, -1, by=-0.5), labels=seq(0, 100, 50), cex.axis=1.4)
mtext(side=1, text="Time from S+ onset (s)", font=2, cex=1.5, line=2.5)
mtext(side=2, text="% Excited", at=0.5, font=2, cex = 1.5, line=2.5)
mtext(side=2, text="% Inhibited", at=-0.5, font=2, cex = 1.5, line=2.5)


#Early AP5 Cue-excited units only
plot.new()
plot.window(xlim = c(0, nrow(EarlyAP5PreInf_ExcBins[,EarlyAP5PreInf_ExcUnits])), ylim=c(0, 1))
abline(h=seq(-1, 1, by=0.25), col="gray90")

lines(x=seq(1, nrow(EarlyAP5PreInf_ExcBins)), y=PercBins(EarlyAP5PreInf_ExcBins[,EarlyAP5PreInf_ExcUnits]), col="gray30", lwd=2)
lines(x=seq(1, nrow(EarlyAP5PostInf_ExcBins)), y=PercBins(EarlyAP5PostInf_ExcBins[,EarlyAP5PreInf_ExcUnits]), col=colindx[2], lwd=2)

axis(side=1, at=seq(0, nrow(EarlyVEHPreInf_ExcBins), by=10), labels=seq(0, 1, by=0.5), cex.axis=1.4)
axis(side=2, las=2, at=seq(0, 1, by=0.5), labels=seq(0, 100, 50), cex.axis=1.4)
axis(side=2, las=2, at=seq(0, -1, by=-0.5), labels=seq(0, 100, 50), cex.axis=1.4)
mtext(side=1, text="Time from S+ onset (s)", font=2, cex=1.5, line=2.5)
mtext(side=2, text="% Excited", at=0.5, font=2, cex = 1.5, line=2.5)


########
#Dot plot (or boxplot) that says: of the cue-excited neurons, during what % of bins were those units excited before and after injection 

EarlyVEHPre_cueExcOnly_ExcPerBin <- EarlyVEHPreInf_ExcBins[,EarlyVEH_ExcUnits]
EarlyVEHPost_cueExcOnly_ExcPerBin <- EarlyVEHPostInf_ExcBins[,EarlyVEH_ExcUnits]

EarlyVEHPre_ExcDotPlot <- colSums(EarlyVEHPre_cueExcOnly_ExcPerBin)/nrow(EarlyVEHPre_cueExcOnly_ExcPerBin)
EarlyVEHPost_ExcDotPlot <- colSums(EarlyVEHPost_cueExcOnly_ExcPerBin)/nrow(EarlyVEHPost_cueExcOnly_ExcPerBin)


EarlyAP5Pre_cueExcOnly_ExcPerBin <- EarlyAP5PreInf_ExcBins[,EarlyAP5_ExcUnits]
EarlyAP5Post_cueExcOnly_ExcPerBin <- EarlyAP5PostInf_ExcBins[,EarlyAP5_ExcUnits]

EarlyAP5Pre_ExcDotPlot <- colSums(EarlyAP5Pre_cueExcOnly_ExcPerBin)/nrow(EarlyAP5Pre_cueExcOnly_ExcPerBin)
EarlyAP5Post_ExcDotPlot <- colSums(EarlyAP5Post_cueExcOnly_ExcPerBin)/nrow(EarlyAP5Post_cueExcOnly_ExcPerBin)


#Make dotplot
makeBoxPlot <- function(data, xmin, xmax, color){
        rect(xleft=xmin, xright=xmax, ybottom=summary(data)[2], ytop=summary(data)[5], border=color, lwd=2) #IQR
        segments(x0=xmin, x1=xmax, y0=summary(data)[3], col=color, lwd=2) #Median
        segments(x0=xmin, x1=xmax, y0=summary(data)[4], col="black", lwd=2) #Mean
}

plot.new()
plot.window(xlim=c(0, 2), ylim=c(0, 0.5))

makeBoxPlot(data=EarlyVEHPre_ExcDotPlot, xmin=0, xmax=0.3, color=colindx[1])
makeBoxPlot(data=EarlyVEHPost_ExcDotPlot, xmin=0.4, xmax=0.7, color=colindx[1])

makeBoxPlot(data=EarlyAP5Pre_ExcDotPlot, xmin=1, xmax=1.3, color=colindx[2])
makeBoxPlot(data=EarlyAP5Post_ExcDotPlot, xmin=1.4, xmax=1.7, color=colindx[2])


axis(side=1, at=c(0.15, 0.55), labels =c("Pre", "Post"), cex.axis=1.5, font=2)
axis(side=1, at=c(1.15, 1.55), labels =c("Pre", "Post"), cex.axis=1.5, font=2)
axis(side=2, las=2, cex.axis=1.4)
mtext(side=2, line=3, text="% excited bins", cex=1.5, font=2)

wilcox.test(EarlyVEHPre_ExcDotPlot, EarlyVEHPost_ExcDotPlot, paired=T) #V = 6.5, p-value = 0.8918
wilcox.test(EarlyAP5Pre_ExcDotPlot, EarlyAP5Post_ExcDotPlot, paired=T) #V =  194.5, p-value = 0.0008807


