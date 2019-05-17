CStask.Rfunc=function(bigoutelement){
out = list(info = bigoutelement$info)

out$receptacleentries = bigoutelement$w
out$receptacleexits = bigoutelement$x
out$CSminuscue = bigoutelement$t
out$CSpluscue = bigoutelement$s
out$laseron = bigoutelement$y
out$rewarddelivery = bigoutelement$u

return(out)
}


save(CStask.Rfunc,file="C:/Users/Mercedes/Desktop/OneDrive - Cuny GradCenter/R functions/CStask.Rfunc")
