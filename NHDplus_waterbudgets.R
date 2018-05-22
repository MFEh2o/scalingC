# water budgets for NHLDplus lakes using lakeCat and our estimates of evap
# SEJ
rm(list=ls())

setwd("~/Documents/Research/LTER_EAGER/scalingC/")

# load lakeCAT tables
setwd("lakeCAT_downloads/")
prism=read.csv("PRISM_1981_2010.csv",header=TRUE,stringsAsFactors=FALSE)
prism=prism[order(prism$COMID),]
runoff=read.csv("Runoff.csv",header=TRUE,stringsAsFactors=FALSE)
runoff=runoff[order(runoff$COMID),]

# load and clean evap tables
setwd("../NHDplusOutput/")
evap=read.csv("NHDplusWaterbody_Evaporation.csv",header=TRUE,stringsAsFactors=FALSE)
evap=evap[evap$COMID%in%runoff$COMID,]
# 6 duplicated COMIDs, lakeCAT found 8...
dups=evap$COMID[duplicated(evap$COMID)]
temp=evap[evap$COMID%in%dups,]
temp=temp[order(temp$COMID),] # information is exactly duplicated so simply removing dups
evap=evap[!duplicated(evap$COMID),]
evap=evap[order(evap$COMID),]

sum(prism$COMID==runoff$COMID)
sum(prism$COMID==evap$COMID)

# create water budgets
# NOTE: lakeCAT convention is that catchment (Cat) is the local land area directly contributing to the lake
#       watershed (Ws) is "upslope" land areas that contribute to other lakes, but then flow into the focal lake
#       it is my understanding that the Cat area is included in the Ws area
#       "Likewise, we used the term watershed to denote a set of nested catchments that represent the full contributing area to a downslope lake."

# For now, I'm calculating budgets like we did for NHLD that ignores Ws and Cat distinction
# I'm just using the Ws area and flows; interesting to compare results when considering lake linkages...

# the long-term average water mass balance is: Qin + Pin = E + Qout
#   Qin = RunoffWs/(365*1000)*WsAreaSqKm*1e6 [m3 day-1]
#   Pin = Precip8110Ws/(365*1000)*AREASQKM*1e6 [m3 day-1]
#   E = EVAPpenman/1000*AREASQKM*1e6 [m3 day-1]
#   Qout = max((Qin + Pin - E), 0)
NHDplus_waterbudgets=data.frame(COMID=prism$COMID,LONG=evap$LONG,LAT=evap$LAT,AREASQKM=evap$AREASQKM,WsAreaSqKm=runoff$WsAreaSqKm,inStreamCat=runoff$inStreamCat,Qin=NA,Pin=NA,E=NA,Qout=NA,flag=NA,stringsAsFactors=FALSE)

NHDplus_waterbudgets$Qin=runoff$RunoffWs/(365*1000)*runoff$WsAreaSqKm*1e6 # m3 day-1
# 321 lakes are missing runoff quantities for the watershed and catchment; use median watershed yield and precip instead...
# could do this a lot better geographically...
NHDplus_waterbudgets$flag[is.na(NHDplus_waterbudgets$Qin)]=1   # NA for runoff in lakeCAT
plot(NHDplus_waterbudgets$LONG[is.na(NHDplus_waterbudgets$Qin)],NHDplus_waterbudgets$LAT[is.na(NHDplus_waterbudgets$Qin)],pch=16,cex=0.5)
library(fields)
US(add=TRUE)
wsYield=runoff$RunoffWs/prism$Precip8110Cat
range(wsYield,na.rm=TRUE)
wsYieldMED=median(wsYield,na.rm=TRUE)  #0.2598
NHDplus_waterbudgets$Qin[is.na(NHDplus_waterbudgets$Qin)]=prism$Precip8110Ws[is.na(NHDplus_waterbudgets$Qin)]*wsYieldMED/(365*1000)*runoff$WsAreaSqKm[is.na(NHDplus_waterbudgets$Qin)]*1e6  # m3 day-1

NHDplus_waterbudgets$Pin=prism$Precip8110Ws/(365*1000)*evap$AREASQKM*1e6  # m3 day-1
NHDplus_waterbudgets$E=evap$EVAPpenman/1000*evap$AREASQKM*1e6             # m3 day-1
putativeQout=NHDplus_waterbudgets$Qin+NHDplus_waterbudgets$Pin-NHDplus_waterbudgets$E
NHDplus_waterbudgets$Qout=ifelse(putativeQout>0,putativeQout,0)           # m3 day-1

#write.csv(NHDplus_waterbudgets,"NHDplus_waterbudgets.csv",row.names=FALSE)

# look at balancing of budgets
range(NHDplus_waterbudgets$Qin)
range(NHDplus_waterbudgets$Pin)
range(NHDplus_waterbudgets$E)
range(NHDplus_waterbudgets$Qout)

sum(putativeQout<0)  # 5398
sum(NHDplus_waterbudgets$Qout==0)  # 5398

sum(NHDplus_waterbudgets$Qin==0)   # 1666

sum((NHDplus_waterbudgets$Qin+NHDplus_waterbudgets$Pin)==0)   # 0
sum((NHDplus_waterbudgets$Pin-NHDplus_waterbudgets$E)<=0)   # 25207
sum((NHDplus_waterbudgets$Qin+NHDplus_waterbudgets$Pin-NHDplus_waterbudgets$E)<=0) # 5398

# estimate volumes from areas
# from Cael et al. 2017 even though it isn't for individual lakes... all other recent models include buffer slope or something
NHDplus_waterbudgets$VOL=10^(1.204*log10(NHDplus_waterbudgets$AREASQKM*1e6)-0.629)     # m3
range(NHDplus_waterbudgets$VOL)
sum(NHDplus_waterbudgets$VOL==0)  # 7 have a zero volume, but that is becaue their area is 0
NHDplus_waterbudgets$flag[NHDplus_waterbudgets$AREASQKM==0]=2  # lake area of 0
NHDplus_waterbudgets$HRT=NHDplus_waterbudgets$VOL/(NHDplus_waterbudgets$E+NHDplus_waterbudgets$Qout) # days

#look at HRT predictions --- not too bad!
range(NHDplus_waterbudgets$HRT)
summary(NHDplus_waterbudgets$HRT)
hist(log10(NHDplus_waterbudgets$HRT),main=NULL,xlab="log10(residence time) [days]")

#look at importance of evap, etc.
fracEvap=NHDplus_waterbudgets$E/(NHDplus_waterbudgets$E+NHDplus_waterbudgets$Qout)
hist(fracEvap,main=NULL,xlab="fraction of water export via evaporation")

plot(fracEvap,log10(NHDplus_waterbudgets$HRT),xlab="fraction of water export via evaporation",ylab="log10(residence time) [days]")

#****** compare to Brooks et al. (2014) isotope-based observations
# load brooks data from ratesVfates repo
brooks=read.csv("../NLAobservations/Brooksetal.csv",header=TRUE,stringsAsFactors=FALSE)
brooks=brooks[brooks$VISIT_NO==1,]
brooks=brooks[order(brooks$SITE_ID),]

# load nla 2007 site info
nla2007=read.csv("../NLAobservations/nla2007_sampledlakeinformation_20091113.csv",header=TRUE,stringsAsFactors=FALSE)
nla2007=nla2007[nla2007$VISIT_NO==1,]
nla2007=nla2007[order(nla2007$SITE_ID),]

sum(nla2007$SITE_ID==brooks$SITE_ID)
dim(nla2007)
dim(brooks)

brooks$COM_ID=nla2007$COM_ID

# try matching with COMIDs
NHDplus4comp=NHDplus_waterbudgets[NHDplus_waterbudgets$COMID%in%brooks$COM_ID,]
dim(NHDplus4comp)  # only 1055 here; should be 1157...
missingBrooksLakes=brooks[!(brooks$COM_ID%in%NHDplus_waterbudgets$COMID),]
sum(duplicated(brooks$COM_ID))
brooks$COM_ID[duplicated(brooks$COM_ID)]
# seems like a few COMIDs are duplicated in brooks, there are also a bunch (~100) that aren't in NHDplus for some reason...

#use minimum distance between lat/long instead...
deltaLat=abs(NHDplus_waterbudgets$LAT-brooks$LAT_DD[1])
deltaLon=abs(NHDplus_waterbudgets$LONG-brooks$LON_DD[1])
dist=sqrt(deltaLat^2+deltaLon^2)

minDist=min(dist)
NHDplus4compDist=NHDplus_waterbudgets[dist==min(dist),]
for(i in 2:nrow(brooks)){
  deltaLat=abs(NHDplus_waterbudgets$LAT-brooks$LAT_DD[i])
  deltaLon=abs(NHDplus_waterbudgets$LONG-brooks$LON_DD[i])
  dist=sqrt(deltaLat^2+deltaLon^2)
  NHDplus4compDist=rbind(NHDplus4compDist,NHDplus_waterbudgets[dist==min(dist),])
  minDist=c(minDist,min(dist))
}

range(minDist)
hist(log10(minDist))  # some of these are pretty far apart...
sum(brooks$COM_ID==NHDplus4compDist$COMID)  # 1009 out of 1157...
range(minDist[brooks$COM_ID==NHDplus4compDist$COMID])
range(minDist[brooks$COM_ID!=NHDplus4compDist$COMID]) # some of these seem right and others don't
matched=hist(log10(minDist[brooks$COM_ID==NHDplus4compDist$COMID]),breaks=seq(-5,0,0.1))
notmatched=hist(log10(minDist[brooks$COM_ID!=NHDplus4compDist$COMID]),breaks=seq(-5,0,0.1)) # some of these seem right and others don't
barplot(rbind(matched$counts,notmatched$counts),beside=TRUE,names.arg=matched$breaks[-1])

# color codes for comparisons
matchCol=rep('black',nrow(brooks))
matchCol[brooks$COM_ID!=NHDplus4compDist$COMID]='red'

plot(log10(brooks$RT*365),log10(NHDplus4compDist$HRT),col=matchCol,xlab="Brooks et al. 2014 HRT (days)",ylab="lakeCAT modeled HRT (days)")
abline(a=0,b=1)

brooksRTcdf=ecdf(log10(brooks$RT*365+1))
modeledRTcdf=ecdf(log10(NHDplus4compDist$HRT))
plot(brooksRTcdf)
plot(modeledRTcdf,add=TRUE,col='red')


plot(brooks$E_I,NHDplus4compDist$E/NHDplus4compDist$Qin,col=matchCol,xlab="Brooks et al. 2014 Evap:Inflow",ylab="lakeCAT modeled Evap:Inflow")
abline(a=0,b=1)
plot(brooks$E_I,NHDplus4compDist$E/NHDplus4compDist$Qin,col=matchCol,xlab="Brooks et al. 2014 Evap:Inflow",ylab="lakeCAT modeled Evap:Inflow",xlim=c(0,1.5),ylim=c(0,1.5))
abline(a=0,b=1)
sum(NHDplus4compDist$E/NHDplus4compDist$Qin>1.5) # 59 out of 1157 not shown on plot because very high

brooksEIcdf=ecdf(brooks$E_I)
modelEIcdf=ecdf(NHDplus4compDist$E/NHDplus4compDist$Qin)
plot(brooksEIcdf)
plot(modelEIcdf,add=TRUE,col='red')
