#function for parsing raw Med PC files

mpcextract_blockSingle.Rfunc=function(rawblock){
singlevars = grep("^[A-Z]{1}:[[:space:]]{0,}[0-9]", rawblock)
arrayvars = grep("^[A-Z]{1}:$", rawblock)
listout = list(info = as.vector(rawblock[c(1:10)]))
if( length(singlevars) > 0) {
  svarlist = as.numeric(substr(rawblock[singlevars], start= 3, stop = 100))
  svarlist = as.list(svarlist)
  names(svarlist) = tolower(substr(rawblock[singlevars], start= 1, stop = 1))
  }

  if( length(singlevars) == 0){svarlist = list(singlevars = "No single variables found")}
if(length(arrayvars) > 0){
  for(j in 1:length(arrayvars)){
    if(j < length(arrayvars)) {iAB = (arrayvars[j]+1):(arrayvars[j+1] - 1)}
    if(j == length(arrayvars)) {iAB = (arrayvars[j]+1):(length(rawblock))}
    arrayblock = rawblock[iAB]
    avars = substr(arrayblock, start = 8, stop = 100)
    avars = unlist(strsplit(avars, split = "[[:space:]]{1,}"))
    avars = as.numeric(avars)
    avars = avars[!is.na(avars)]
    if(j == 1 ){avarlist = list(avars)}
      else {avarlist[[j]] = avars}
    } #End for(j . . .
    names(avarlist) = tolower(substr(rawblock[arrayvars], start= 1, stop = 1))
  }# End if(length(arrayvars . . .
    if(length(arrayvars) == 0){avarlist = list(arrayvars = "No array variables found")}
listout = c(listout, svarlist, avarlist)
return(listout)
}



save(mpcextract_blockSingle.Rfunc,file="C:/Users/Cindy/Desktop/R scripts/Functions/mpcextract_blockSingle.Rfunc")