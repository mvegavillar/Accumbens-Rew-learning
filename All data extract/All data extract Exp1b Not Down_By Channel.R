
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

  

############### ANALYSIS AVERAGING DATA PER CHANNEL (AS PER SALEEM'S EMAIL 3/26/19)

#Channel index

ch_idx <- do.call("rbind",
                  
                  lapply(seq(1, length(allNeuronsDS$nexdata)), function(x){
                          
                          nexdataSel <- allNeuronsDS$nexdata[[x]]
                          
                          subject <- nexdataSel$ratname
                          day <- nexdataSel$expt
                          
                          channelsIdx <- names(nexdataSel)[grep("sig", names(nexdataSel))]
                          
                          channelNumber <- as.numeric(gsub("([^0-9.]+)", "", channelsIdx))
                          
                          data.frame(subject, day, channelNumber)
                          
                  })
)

ch_idx$allUnitIdx <- 1:nrow(ch_idx)

#toplot has average firing rate of each unit, in Z scores, 100-400 ms after DS and after NS (it was calculated using megaplot())

#toplot <- megaplot(data=list(allNeuronsDS, allNeuronsNS), CPdata=CPdata, csacqidx=csacqidx, BLwdw=2, winmin=100, winmax=400,
#                   colpalette="Rainbow", minFR=-2, maxFR=6, graphFolder=MixedGraphFolder, dataForRdir=dataForRdir,
#                   ZcolLabels=c("ZDS", "ZNS"), arrangeBy="ZDS")

load(file=paste(dataForRdir, "toplot.rdat", sep=""))

#Data frame with channel info, firing info and behavior info

ch_df <- merge(ch_idx, toplot)

ch_df$CHANNEL <- paste(ch_df$subject, ch_df$channelNumber)

#Average activity per channel

CH_DF <- do.call('rbind', lapply(seq(1, length(unique(ch_df$CHANNEL))), function(x){
        
        ch <- unique(ch_df$CHANNEL)[x]
        
        CHsel <- ch_df[ch_df$CHANNEL==ch,]
        
        do.call("rbind", lapply(seq(1, length(unique(CHsel$day))), function(y){
                
                day <- unique(CHsel$day)[y]
                
                dayCHsel <- CHsel[CHsel$day==day, ]
                
                toprint <- dayCHsel[1,]
                
                toprint$ZDS <- mean(dayCHsel$ZDS, na.rm=T)
                toprint$ZNS <- mean(dayCHsel$ZNS, na.rm=T)
                
                toprint$DSminusNS <- mean(dayCHsel$ZDS-dayCHsel$ZNS, na.rm = T)
                toprint$DSplusNS <- mean(dayCHsel$ZDS+dayCHsel$ZNS, na.rm = T)
                
                toprint$specFR <- mean(dayCHsel$ZDS-dayCHsel$ZNS, na.rm = T)/mean(dayCHsel$ZDS+dayCHsel$ZNS, na.rm = T)
                
                toprint
        })
        )
})
)



missingData <- table(CH_DF$CHANNEL, CH_DF$CPfromSess)[,c(4:1, 5:9)]

save(missingData, file=paste(dataForRdir, "missingData.rdat", sep=""))
write.csv(missingData, file=paste(dataForRdir, "missingData.csv", sep=""))

#Subset the channels that had data on day -1, 0 AND 1 from CP.

goodChannels <- do.call("rbind", lapply(seq(1, length(unique(CH_DF$CHANNEL))), function(x){
        selCH <- CH_DF[CH_DF$CHANNEL==unique(CH_DF$CHANNEL)[x], ]
        
        targetDays <- c(-1, 0, 1)
        sessPerCH <- selCH$CPfromSess
        
        sessIdx <- sessPerCH %in% targetDays
        
        if(sum(sessIdx)==3){
                selCH[sessIdx, ]
        }
})
)

library(ez)

ezANOVA(goodChannels, dv=ZDS, within=CPfromSess, wid=CHANNEL)
# $`ANOVA`
# Effect      DFn DFd        F            p p<.05       ges
# CPfromSess   2  28  9.178532 0.0008601803     * 0.1313192
# 
# $`Mauchly's Test for Sphericity`
# Effect         W          p p<.05
# CPfromSess 0.5584255 0.02266074     *
#         
# $`Sphericity Corrections`
# Effect      GGe       p[GG]       p[GG]<.05       HFe       p[HF] p[HF]<.05
# CPfromSess 0.693686   0.003599983         * 0.7457042 0.002818656         *


ezANOVA(goodChannels, dv=ZDS, within=CPfromSess, wid=CHANNEL)
#$`ANOVA`
#      Effect DFn DFd        F           p p<.05       ges
#2 CPfromSess   2  28 7.211332 0.002977429     * 0.1223831
#
#$`Mauchly's Test for Sphericity`
#      Effect         W        p p<.05
#2 CPfromSess 0.6847527 0.085304      
#
#$`Sphericity Corrections`
#  Effect       GGe       p[GG]     p[GG]<.05       HFe       p[HF] p[HF]<.05
#2 CPfromSess 0.7603133 0.006969813         * 0.8337517 0.005364974         *

#Plot:
plot.new()
plot.window(xlim=c(0, 4), ylim=c(-4, 12))

datapoints <- do.call("cbind", lapply(seq(1, length(unique(goodChannels$CPfromSess))), function(x){
        
        selCH <- goodChannels[goodChannels$CPfromSess==unique(goodChannels$CPfromSess)[x], ]
        
        summ <- summary(selCH$DSminusNS)
        
        median <- summ[3]; Q1 <- summ[2]; Q3 <- summ[5]
        
        rect(xleft=x-0.3, xright=x+0.3,
             ybottom=Q1, ytop=Q3, col="gray50", border="white")
        
        segments(x0=x-0.3, x1=x+0.3, y0=median, y1=median, col="white", lwd=3)
        
        dots <- selCH$DSminusNS
        
        points(x=rep(x, length(dots)), y=dots, pch=16, col="black")
        
        dots
        
})
)

axis(side=1, at=seq(1:3), labels=c(-1, 0, 1)); mtext(side=1, line=2.5, text="Session from CP", font=2)
axis(side=2, at=seq(-4, 12, 2), las=2); mtext(side=2, line=2.5, text="Post-S+ - Post-S- firing rate (Z sc.)", font=2)

abline(h=0, lty=3)


### We can actually treat DS vs. NS firing rate as another within-subject factor. Let's do that

#install.packages("tidyr")
library(tidyr)

#Make a long format version of the data frame in which "Cue" is a variable (DS vs. NS) and FR is the firing rate in Z scores for that cue
data_long <- gather(goodChannels, Cue, FR, ZDS:ZNS, factor_key=TRUE)


ezANOVA(data_long, dv=FR, within=c(CPfromSess, Cue), wid=CHANNEL)
# $`ANOVA`
#           Effect DFn DFd         F            p p<.05        ges
# 2     CPfromSess   2  28  5.092506 0.0129930602     * 0.05524460
# 3            Cue   1  14 19.631766 0.0005700603     * 0.35606325
# 4 CPfromSess:Cue   2  28  9.178532 0.0008601803     * 0.07190359
# 
# $`Mauchly's Test for Sphericity`
#           Effect         W          p p<.05
# 2     CPfromSess 0.8456711 0.33636410      
# 4 CPfromSess:Cue 0.5584255 0.02266074     *
#         
# $`Sphericity Corrections`
#       Effect           GGe       p[GG] p[GG]<.05       HFe       p[HF] p[HF]<.05
# 2 CPfromSess     0.8663042 0.017580392         * 0.9777599 0.013661473         *
# 4 CPfromSess:Cue 0.6936860 0.003599983         * 0.7457042 0.002818656         *


#Plot

plot.new()
plot.window(xlim=c(0, 4), ylim=c(-2, 12))

colors <- c(colindx[1], "darkblue")

datapoints <- lapply(seq(1, length(unique(data_long$CPfromSess))), function(x){
        
        #Select data from one session (with respect to CP)
        selCH <- data_long[data_long$CPfromSess==unique(data_long$CPfromSess)[x], ]
        
        #For each kind of cue (DS vs. NS)
        do.call("cbind", lapply(seq(1, length(unique(selCH$Cue))), function(y){
                selCue <- selCH[selCH$Cue==unique(selCH$Cue)[y], ]
                
                summ <- summary(selCue$FR)
                
                median <- summ[3]; Q1 <- summ[2]; Q3 <- summ[5]
                
                boxwidth <- 0.3
                XLEFT <- x+(-2+y)*0.3
                XRIGHT <- XLEFT + 0.3
                
                rect(xleft=XLEFT, xright=XRIGHT,
                     ybottom=Q1, ytop=Q3, col=colors[y], border="white")
                
                segments(x0=XLEFT, x1=XRIGHT, y0=median, y1=median, col="white", lwd=3)
                
                dots <- selCue$FR
                # 
                # points(x=rep(XLEFT+(XRIGHT-XLEFT)/2, length(dots)), y=dots, pch=16, col="black")
                # 
                dots
                
        })
        )
})

targetSessions <- c(-1, 0, 1)
axis(side=1, at=seq(1:3), labels=targetSessions); mtext(side=1, line=2.5, text="Session from CP", cex=1.4, font=2)
axis(side=2, at=seq(-2, 12, by=2), las=2); mtext(side=2, line=2.5, text="Firing 100-400 ms after the cue (Z sc.)", cex=1.4, font=2)

abline(h=0, lty=3)
legend("topright", fill=colors, legend=c("S+", "S-"), cex=1.4)
text(x=3, y=-1.5, labels = paste("n = ", nrow(datapoints[[1]]), " channels", sep=""))

#Post-hoc comparisons (S+ vs. S- firing rate) across sessions from CP
PostHoc_Comp <- do.call("rbind", lapply(seq(1, length(datapoints)), function(x){
        comp <- t.test(x=datapoints[[x]][,1], y=datapoints[[x]][,2], paired=T, alternative="greater")
        
        data.frame(sessFromCP=targetSessions[x], t=comp$statistic, df=comp$parameter, p=comp$p.value)
})
)

PostHoc_Comp$p.adj <- p.adjust(p=PostHoc_Comp$p, method="holm")
PostHoc_Comp$sig <- giveStars(PostHoc_Comp$p.adj)
#    sessFromCP        t df            p        p.adj sig
# t          -1 2.798498 14 0.0071111483 0.0071111483  **
# t1          0 3.504997 14 0.0017501370 0.0035002739  **
# t2          1 5.097081 14 0.0000812829 0.0002438487 ***

sapply(seq(1, nrow(PostHoc_Comp)), function(x){
        text(x=x, y=8, labels = PostHoc_Comp$sig[x], cex=1.4)
})


########## ACTUALLY, I SHOULD USE NON-PARAMETRIC COMPARISONS! ####################
#Across session comparisons:
everyAcrossSessComp <- do.call("rbind", lapply(seq(1, 2), function(x){
        comp=wilcox.test(x=datapoints[[x]][,1], y=datapoints[[x+1]][,1], paired=T, alternative="less")
        data.frame(sess=paste(targetSessions[x], " vs. ", targetSessions[x+1], sep=""),
                   statistic=comp$statistic, p=comp$p.value, alt=comp$alternative)
})
)
#       sess statistic          p  alt
#V  -1 vs. 0        27 0.03186035 less
#V1  0 vs. 1        29 0.04162598 less

#Only -1 vs. 1 comparison
wilcox.test(x=datapoints[[1]][,1], y=datapoints[[3]][,1], paired=T, alternative="less")
# Wilcoxon signed rank test

# data:  datapoints[[1]][, 1] and datapoints[[3]][, 1]
# V = 3, p-value = 0.0001526
# alternative hypothesis: true location shift is less than 0


#Post-hoc comparisons (S+ vs. S- firing rate) across sessions from CP
PostHoc_Comp <- do.call("rbind", lapply(seq(1, length(datapoints)), function(x){
        comp <- wilcox.test(x=datapoints[[x]][,1], y=datapoints[[x]][,2], paired=T, alternative="greater")
        
        data.frame(sessFromCP=targetSessions[x], W=comp$statistic, p=comp$p.value, alt=comp$alternative)
})
)

p.adjust(p=c(0.0001526, PostHoc_Comp$p), method="holm")
#[1] 0.0006104000 0.0041809082 0.0015258789 0.0006408691