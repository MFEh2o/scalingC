# PRISM monthlys data download and calculation of PET using Hamon equation
# SEJ

# running this on jones macpro to get done in a somewhat timely fashion by parallelization
# copied required files (Runoff.csv and NHDplusWaterbody_Evaporation.csv) 
# from directories in local version of github repo to
# ~/Documents/jones/HamonEvap_May2018

rm(list=ls())

setwd("~/Documents/jones/HamonEvap_May2018/")

# load lakeCAT and evap table to get lake lat-lon
runoff=read.csv("Runoff.csv",header=TRUE,stringsAsFactors=FALSE)
runoff=runoff[order(runoff$COMID),]
evap=read.csv("NHDplusWaterbody_Evaporation.csv",header=TRUE,stringsAsFactors=FALSE)
evap=evap[evap$COMID%in%runoff$COMID,]
# 6 duplicated COMIDs, lakeCAT found 8...
dups=evap$COMID[duplicated(evap$COMID)]
temp=evap[evap$COMID%in%dups,]
temp=temp[order(temp$COMID),] # information is exactly duplicated so simply removing dups
evap=evap[!duplicated(evap$COMID),]
evap=evap[order(evap$COMID),]

# work with prism data
library(prism)
library(raster)
options(prism.path="./prismDownloads")

# pull monthly normals for mean air temperature from prism database
for(i in 1:12){
  get_prism_normals(type="tmean",resolution="4km",mon=i,keepZip=FALSE)
}

# generate prism raster stack to extract data from
TMEANstack=prism_stack(ls_prism_data()[,1])
coordList=as.list(as.data.frame(t(evap[,2:3])))

# define function for pulling monthly norms for each NHDplus lake
pullMonthlyTMEAN<-function(location){
  unlist(extract(TMEANstack,matrix(location,nrow=1),buffer=10))
}

library(parallel)
out=mclapply(coordList,pullMonthlyTMEAN,mc.cores=8)

monthlyTMEANlong=data.frame(COMID=rep(evap$COMID,each=12),month=rep(1:12,nrow(evap)),TMEAN=unlist(out),stringsAsFactors=FALSE)

monthlyTMEANwide=reshape(monthlyTMEANlong,idvar="COMID",timevar="month",direction="wide")

write.csv(monthlyTMEANwide,"monthlyTMEANnorms_prism.txt",row.names=FALSE)

# put a copy of this output in /scalingC/NHDplusOutput/