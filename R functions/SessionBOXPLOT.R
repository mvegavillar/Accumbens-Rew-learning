SessionBOXPLOT <- function(experiment="Exp 4 Unil AP5 EXT", masterDF=list(masterDF_DS_L_VEH, masterDF_DS_L_AP5), 
         graphFolder=MixedGraphFolder, dataProcess="Zscores", correctOnly=FALSE, color=colindx, 
         yAxMinZ = -1, yAxMaxZ = 5, yAxMaxRaw = 10, imgFormat="pdf", WdwStart=100, WdwEnd=400,
         capped=F, capValue=c(1, 15), neudata=list(allNeuronsDS_L_VEH, allNeuronsDS_L_AP5), 
         cueExcOnly=FALSE, legendLabels=c("S+ VEH side L", "S+ AP5 side L")){
        
        
        # Install and call necessary packages
        if(!require(dplyr)){install.packages("dplyr")}
        library(dplyr)
        if(!require(matrixStats)){install.packages("matrixStats")}
        library(matrixStats)
        
        
        if(correctOnly==TRUE){trialSel="correctOnly"} else {trialSel="all trials"}
        if(cueExcOnly==TRUE){unitSel="CueExcOnly"} else {unitSel="all units"}
        
        if(imgFormat=="pdf"){pdf(file=paste(graphFolder, experiment, "BOXPLOT", dataProcess, trialSel, unitSel, ".pdf", sep="_"))}
        if(imgFormat=="png"){pnf(filename=paste(graphFolder, experiment, "BOXPLOT", dataProcess, trialSel, unitSel, ".pdf", sep="_"))}
        
        binw <- neudata[[1]]$parameters$binw
        if(is.null(binw)){binw <- neudata$parameters$binw} #Preventive troubleshooting: In case I define neudata not as a list (as it should), but as a single object
        
        minBin <- WdwStart/binw
        maxBin <- WdwEnd/binw
        
        selBins <- minBin:maxBin
        
        
        plot.new()
        
        byUnit <- do.call("list", lapply(seq(1, length(masterDF)), function(c){
                
                FRcols <- (1:ncol(masterDF[[c]]))[is.element(colnames(masterDF[[c]]), 1:ncol(masterDF[[c]]))]
                
                subFRcolNames <- (unique(masterDF[[c]]$CueBin)+minBin):(unique(masterDF[[c]]$CueBin)+maxBin) #Select bins to be plotted
                subFRcols <- colnames(masterDF[[c]]) %in% as.character(subFRcolNames)
                
                dataSel <- masterDF[[c]]
                
                #Cue excited units taking into account the whole session
                cueexidx <- neudata[[c]]$cueexidx
                propcueex <- sum(unlist(cueexidx))/length(unlist(cueexidx))
                
                if(capped==TRUE){
                        dataSel <- dataSel[dataSel$trialIdx>=capValue[1] & dataSel$trialIdx<=capValue[2], ]
                        
                        masterlist <- neudata[[c]]$masterlist
                        
                        #Cue exc units taking into account those specific trials only
                        cueexidx = lapply(seq(1, length(masterlist)), function(x) { #For each rat
                                lapply(seq(1, length(masterlist[[x]])), function(y) { #For each unit
                                        
                                        trialSel <- capValue[1]:capValue[2]
                                        pbin = .05
                                        threshold = 3
                                        BLwdw=2
                                        
                                        allvals = unlist(masterlist[[x]][[y]][trialSel])
                                        hcounts = hist(allvals[which(allvals >= -BLwdw & allvals <= .5)], breaks = seq(-BLwdw, .5, pbin), plot = F)$counts 
                                        
                                        freq = (hcounts/pbin)/length(masterlist[[x]][[y]])  #just for comparison's sake, to see how counts compare to actual frequency
                                        baseline = mean(freq[1:(BLwdw/pbin)])   
                                        
                                        baselinecounts=hcounts[1:(BLwdw/pbin)]
                                        critwin = hcounts[(BLwdw/pbin + 1):((BLwdw+.5)/pbin)]         #this is the window (500 ms) in which we look for excitation
                                        
                                        critval = poisson.test(x = sum(baselinecounts), T = length(baselinecounts), conf.level = 0.999)$conf.int[2]   #computes the upper limit of the confindence interval based on the baseline
                                        
                                        diffs = diff(which(critwin > critval))       #computes the differences in indices for bins exceeding the critical value (to check whether excited bins are consecutive or not)
                                        cueex = F
                                        if(length(which(rle(diffs)$values == 1 & rle(diffs)$lengths >= threshold))>0) (cueex = T)   #looks for consecutive bins (diff equal to 1) that are at least 2 bins long -or whatever #of bins specified by the argument "threshold" (rle length of at least 2)
                                        
                                        return(cueex)  
                                })
                        })
                        
                        
                        # propcueex <- lapply(seq(1, length(cueexidx)), function(x){
                        #         length(which(cueexidx[[x]]==T))/length(cueexidx[[x]])
                        # })
                        
                        propcueex <- sum(unlist(cueexidx))/length(unlist(cueexidx))
                }
                
                #Subset the trials (correct only or all) and the units (cue excited or not) to be used for the PSTH
                if(correctOnly==TRUE){
                        dataSel <- dataSel[!is.na(dataSel$CueResponse), ] 
                }
                
                
                if(cueExcOnly==TRUE){
                        
                        cueExc <- unlist(cueexidx)[dataSel$allUnitIdx]
                        dataSel$CueExcited <- cueExc
                        
                        dataSel <- dataSel[dataSel$CueExcited==TRUE, ]
                        
                }
                
                ZscoreCalc <- function(x, avg, sd){(x-avg)/sd}
                
                #Calculate the mean average rate per bin around the time of the event of interest (in Z scores or raw scores)
                if(dataProcess=="Zscores"){
                        plot.window(xlim=c(0, length(masterDF)+0.5), ylim=c(yAxMinZ, yAxMaxZ))
                        numericDF <- apply(dataSel[, subFRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                        BLaverage <- as.numeric(format(dataSel$BLavg, digits=2)) #Baseline info is in integer format. If I just say numeric, it'll remove the decimal point and do sth weird. So I have to recur to this roundabout way.
                        BLsd <- as.numeric(format(dataSel$BLsd, digits=2))
                        numericDFZsc <- sapply(seq(1, nrow(numericDF)), function(w){
                                sapply(seq(1, length(numericDF[w,])), function(v){
                                        ZscoreCalc(x=numericDF[w,][v], avg=BLaverage[w], sd=BLsd[w])})
                        })
                        numericDFZsc <- t(numericDFZsc)
                        
                        MeanByUnit <- rowMeans(numericDFZsc, na.rm=T)
                        SEMbyUnit <- rowSds(numericDFZsc, na.rm=T)/sqrt(ncol(numericDFZsc))
                        yAxMax=yAxMaxZ
                        
                        
                } else {
                        plot.window(xlim=c(0, length(masterDF)+0.5), ylim=c(0, yAxMaxRaw))
                        numericDF <- apply(dataSel[,subFRcols], MARGIN=2, as.numeric) #Convert selected FR columns into numeric
                        MeanByUnit <- rowMeans(numericDF, na.rm=T)
                        SEMbyUnit <- rowSds(numericDF, na.rm=T)/sqrt(ncol(numericDF))
                        yAxMax=yAxMaxRaw
                        labelLeg=paste(yAxMax, "(Hz)")
                }
                
                
                
                #Calculate the indices of central tendency and dispersion for the boxplot
                summFR <- summary(MeanByUnit)
                
                Q1=summFR[2]; Q3=summFR[5]; mean=summFR[4]; median=summFR[3]
                
                #Plot the boxplots and mean (black) and median (white)
                rect(xleft=c-0.2, xright=c+0.2, ybottom = Q1, ytop=Q3, col = color[c], border = NA)
                segments(x0=c-0.2, x1=c+0.2, y0=mean, y1=mean, col = "black", lwd=2)
                segments(x0=c-0.2, x1=c+0.2, y0=median, y1=median, col = "white", lwd=2)
                
                results <- list(MeanByUnit)
                names(results) <- legendLabels[c]
                
                nUnits <- length(unique(dataSel$allUnitIdx))
                text(x=10, y=5-c, labels=paste("n=", nUnits, "; ", "Prop Cue Exc=", round(propcueex, 2)), col=colindx[c])
                
                return(results)              
        })
        )
        
        abline(h=0, col="gray30", lty=3)
        axis(side=1,at=seq(1, length(legendLabels)), labels=legendLabels)
        axis(side=2, las=2)
        mtext(side=2, "Firing rate", line=2.5, font=2, cex=1.2)   
        
        forPrint <- byUnit
        
        if(length(byUnit)==2){
                first <- unlist(byUnit[[1]])
                second <- unlist(byUnit[[2]])
                
                # if(mean(first, na.rm=T) < mean(second, na.rm=T)){alt.pick <- "less"}
                # if(mean(first, na.rm=T) > mean(second, na.rm=T)){alt.pick <- "greater"}
                # if(mean(first, na.rm=T) == mean(second, na.rm=T)){alt.pick <- "two.sided"}
                
                alt.pick <- "two.sided"
                
                wtest <- wilcox.test(x=first, y=second, paired=F, alternative=alt.pick)
                
                if(wtest$p.value<0.05){text(labels=giveStars(wtest$p.value), x=1.5, y=yAxMax)}
                
                forPrint <- list(byUnit, wtest)
        }
        
        
        dev.off()
        
        print(forPrint)
        
        
}