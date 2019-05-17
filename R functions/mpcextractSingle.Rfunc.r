mpcextractSingle.Rfunc=function(mpcfile, outputunknowns = F, writetofile = F){
outputunknowns = F
writetofile = F
#load ("C:/Users/Kevin Caref/Desktop/RScripts/Functions/Cued_Licking_3Solutions.Rfunc")
rawin = scan(mpcfile, what = "char", sep = "\n")
DATE =  grep("Start Date:", rawin)
BOX = grep("Box: ", rawin)
DATASTART = grep("MSN:", rawin)
blocklist = mpcextract_blockSingle.Rfunc(rawin)
bigoutlist = c(blocklist)
return(bigoutlist)
}

save(mpcextractSingle.Rfunc,file="C:/Users/Cindy/Desktop/R scripts/Functions/mpcextractSingle.Rfunc")