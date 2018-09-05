# PRISM monthlys data download and calculation of PET using Hamon equation
# SEJ
rm(list=ls())

# load lakeCAT and evap table to get lake lat-lon
setwd("~/Documents/Research/LTER_EAGER/scalingC/lakeCAT_downloads/")
runoff=read.csv("Runoff.csv",header=TRUE,stringsAsFactors=FALSE)
runoff=runoff[order(runoff$COMID),]
setwd("../NHDplusOutput/")
evap=read.csv("NHDplusWaterbody_Evaporation.csv",header=TRUE,stringsAsFactors=FALSE)
evap=evap[evap$COMID%in%runoff$COMID,]
# 6 duplicated COMIDs, lakeCAT found 8...
dups=evap$COMID[duplicated(evap$COMID)]
temp=evap[evap$COMID%in%dups,]
temp=temp[order(temp$COMID),] # information is exactly duplicated so simply removing dups
evap=evap[!duplicated(evap$COMID),]
evap=evap[order(evap$COMID),]

# load monthly normals extracted from PRISM using hamonEvapFromPrismParallel.R
setwd("~/Documents/Research/LTER_EAGER/scalingC/NHDplusOutput/")
tmean=read.csv("monthlyTMEANnorms_prism.txt",header=TRUE)

sum(tmean$COMID==evap$COMID)

# function for calculating evap using Hamon's Equation from Rosenberry et al. 2007 and Zwart
evap_func <- function(airT,jDay,lat){
  #saturated Vapor Density
  svd<-5.018+(.32321*airT)+(0.0081847*airT^2)+(0.00031243*airT^3)
  
  #Daylength
  degToRad<-2*pi/360
  radToDeg<-180/pi
  
  #day angle gamma (radians)
  dayAngle<-2*pi*(jDay-1)/365
  
  #declination of the sun 'delta' (radians)
  dec<-0.006918-0.399912*cos(dayAngle)+0.070257*sin(dayAngle)-
    0.006758*cos(2*dayAngle)+0.000907*sin(2*dayAngle)-0.002697*
    cos(3*dayAngle)+0.00148*sin(3*dayAngle)
  
  # sunrise hour angle 'omega' (degrees)
  latRad<-lat*degToRad
  sunriseHourAngle<-acos(-tan(latRad)*tan(dec))*radToDeg
  
  #sunrise and sunset times (decimal hours, relative to solar time)
  sunrise<-12-sunriseHourAngle/15
  sunset<-12+sunriseHourAngle/15
  dayLength<-sunset-sunrise #day length in hours
  
  evap = 0.55*((dayLength/12)^2)*(svd/100)*25.4 #calculates evaporation for each jDay (units are mm/day)
  return(evap)
}


# day length constants
degToRad<-2*pi/360
radToDeg<-180/pi

#day angle gamma (radians)
dayAngle<-2*pi*((1:365)-1)/365

#mean monthly declination of the sun 'delta' (radians)
dec<-0.006918-0.399912*cos(dayAngle)+0.070257*sin(dayAngle)-
  0.006758*cos(2*dayAngle)+0.000907*sin(2*dayAngle)-0.002697*
  cos(3*dayAngle)+0.00148*sin(3*dayAngle)

dec12=c(mean(dec[1:31]),mean(dec[32:59]),mean(dec[60:90]),mean(dec[91:120]),mean(dec[121:151]),mean(dec[152:181]),mean(dec[182:212]),mean(dec[213:243]),mean(dec[244:273]),mean(dec[274:304]),mean(dec[305:334]),mean(dec[335:365]))

evap_func <- function(airT,lat){
  #saturated Vapor Density
  svd<-5.018+(.32321*airT)+(0.0081847*airT^2)+(0.00031243*airT^3)
  
  # sunrise hour angle 'omega' (degrees)
  latRad<-lat*degToRad
  sunriseHourAngle<-acos(-tan(latRad)*tan(dec12))*radToDeg
  
  #sunrise and sunset times (decimal hours, relative to solar time)
  sunrise<-12-sunriseHourAngle/15
  sunset<-12+sunriseHourAngle/15
  dayLength<-sunset-sunrise #day length in hours
  
  evap12 = 0.55*((dayLength/12)^2)*(svd/100)*25.4 #calculates evaporation for each jDay (units are mm/day)
  annualEvap=mean(evap12)
  return(annualEvap)
}

hamonEvap=numeric(nrow(evap))
for(i in 1:nrow(evap)){
  print(i)
  hamonEvap[i]=evap_func(as.numeric(tmean[i,2:13]),evap$LAT[i])
}

# THINK ABOUT ICE COVER???

# output results
evap$EVAPhamon=hamonEvap

write.csv(evap,"NHDplusWaterbody_EvaporationIncludingHamon.csv",row.names=FALSE)
