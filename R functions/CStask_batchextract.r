funcdirect="C:/Users/Cindy/Desktop/R scripts/Functions/"

load (paste(funcdirect,"CStask.Rfunc",sep=""))
load (paste(funcdirect,"mpcextract_blockSingle.Rfunc",sep=""))
load (paste(funcdirect,"mpcextractSingle.Rfunc",sep=""))

mpcdatadir = "C:/Users/Cindy/Desktop/Test/"

allfiles = list.files(mpcdatadir)

alldata = list()
for (i in(1:length(allfiles))) {                                          
  rawmpc = mpcextractSingle.Rfunc(paste(mpcdatadir, allfiles[i], sep = ""))
  toanalyze = CStask.Rfunc(rawmpc)

  EXPTindx=grep("Experiment",rawmpc[[1]])
  EXPT=substr (rawmpc[[1]][EXPTindx],13,nchar(rawmpc[[1]][EXPTindx]))
  testname=unlist (strsplit(rawmpc$info[1],"Subject"))
  ratname=substr (testname[2],2,nchar(testname[2]))
  idate=unlist (strsplit(rawmpc$info[1],"!"))
  testdate=substr(idate[2],1,10)
  
  CSpluscue = toanalyze$CSpluscue
  CSminuscue = toanalyze$CSminuscue
  receptacleentries = toanalyze$receptacleentries
  rewarddelivery = toanalyze$rewarddelivery 
  
  CSplusresponse = vector()
  for(j in 1:length(CSpluscue)) {
    CSplusentries = receptacleentries[receptacleentries > CSpluscue[j] & receptacleentries < CSpluscue[j]+5]
    if (length(CSplusentries) > 0)(CSplusresponse[j] = min(CSplusentries))
    else(CSplusresponse[j] = NA)
  }
  CSplusresponse = CSplusresponse[-which(is.na(CSplusresponse))]
  
  CSminusresponse = vector()
  for(k in 1:length(CSminuscue)) {
    CSminusentries = receptacleentries[receptacleentries > CSminuscue[k] & receptacleentries < CSminuscue[k]+5]
    if (length(CSminusentries) > 0)(CSminusresponse[k] = min(CSminusentries))
    else(CSminusresponse[j] = NA)
  }
  
  CSminusresponse = CSminusresponse[-which(is.na(CSminusresponse))]
  
  
  alldata[[i]] = list(ratname = ratname, testdate = testdate, expt = EXPT, receptacleentries = receptacleentries, receptacleexits = toanalyze$receptacleexits,
    CSminuscue = CSminuscue, CSpluscue = CSpluscue, CSminusresponse = CSminusresponse, CSplusresponse = CSplusresponse, laseron = toanalyze$laseron,  rewarddelivery = toanalyze$rewarddelivery)
  
  if(sum(alldata[[i]]$CSplusresponse == 0)) (alldata[[i]]$CSplusresponse = 0) 
  if(sum(alldata[[i]]$CSminusresponse == 0)) (alldata[[i]]$CSminusresponse = 0)
     
  }

allrats=unlist (lapply(alldata,FUN=function(x)({
toget=x$ratname
return(toget)
})))

allexpts=unlist (lapply(alldata,FUN=function(x)({
toget2=x$expt
return(toget2)
})))

alldates=unlist (lapply(alldata,FUN=function(x)({
toget3=x$testdate
return(toget3)
})))                         

