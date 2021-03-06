# NLDAS_VIC evaporation calculations
# SEJ
# Using output from VIC run asssociated with the North American Land Data Assimilation System Phase 2 (NLDAS-2) project 
# to estimate evaporation according to the Penman equation

#******************
# function for calculating evap from VIC output from NLDAS-2
#******************
rm(list=ls())

evap_calc<-function(pressure,tair,shortwave,longwave,SH){
#evap_calc<-function(elev,tair,shortwave,longwave,SH){
# REQUIRED INPUTS FOR THE LOCATION
# elevation [elev; m]
# air temperature [tair; C]
# incoming shortwave radiation [shortwave; ????]
# incoming longwave radiation [longwave]
# specific humidity [SH; Kg Kg-1]

# OUTPUT 
# evaporation [evap; mm day-1]

# PARAMETERS
LAPSE_PM = -0.006  #environmental lapse rate [C m-1]
PS_PM = 101300  #sea level air pressure [Pa]
CP_PM = 1013  #specific heat of moist air at constant pressure [J kg-1 C-1]
Z0_Lower = 0.001  #roughness
d_Lower = 0.0054  #displacement
von_K = 0.4  #Von Karman constant for evapotranspiration
K2 = von_K*von_K
C = 2.16679 #constant for specific humidity to vapor pressure conversion [g K J-1] 
A_SVP = 0.61078  #A term in saturated vapor pressure calculation
B_SVP = 17.269  #B term in saturated vapor pressure calculation
C_SVP = 237.3  #C term in saturated vapor pressure calculation
SEC_PER_DAY = 86400  #seconds per day
H2O_SURF_ALBEDO = 0.08  #albedo of water surface
STEFAN_B = 5.6696e-8  #stefan-boltzmann constant [W/m^2/K^4]
# height of wind measuremtn is 10 m (not sure if this is needed)

# INTERMEDIATE EQUATIONS
#h=287/9.81*((tair+273.15)+0.5*elev*LAPSE_PM)  #scale height in the atmosphere [m?]

#pz=PS_PM*exp(-elevation/h)  #surface air pressure [Pa]
pz=pressure  #surface air pressure [Pa]

lv=2501000-2361*tair  #latent heat of vaporization [J Kg-1]

gamma=1628.6*pz/lv  #psychrometric constant [Pa C-1]

r_air=0.003486*pz/(275+tair)  #air density [Kg m-3]

rs=0  #minimal stomatal resistance [s m-1]
rc=0  #

ra=log((2+(1/0.63-1)*d_Lower)/Z0_Lower)*log((2+(1/0.63-1)*d_Lower)/(0.1*Z0_Lower))/K2 #aerodynamic resistance [s m-1]
  
rarc=0  #architectural resistance [s m-1]

svp=A_SVP*exp((B_SVP*tair)/(C_SVP+tair))  #saturated vapor pressure [Pa]

vp=SH*r_air*1000*(tair+273.15)/C  #vapor pressure [Pa]

vpd=svp-vp  #vapor pressure deficit [Pa]

slope=(B_SVP*C_SVP)/((C_SVP+tair)*(C_SVP+tair))*svp  #slope of saturated vapor pressure curve [Pa K-1]

net_short = (1-H2O_SURF_ALBEDO)*shortwave  # [W m-2]

# From VIC func_surf_energy_bal.c
#Tmp = Ts + KELVIN; Ts is soil temperature or surface temperature and KELVIN=273.15
#LongBareOut = STEFAN_B * Tmp * Tmp * Tmp * Tmp; 

#a lake study (Binyamin et al. 2006; Int. J. Climatol. 26: 2261-2273) used E*sigma*T^4 as outgoing and E*longwave as incoming, where
# E is emissivity of the water surface (0.97), sigma is the Stefan-Boltzman constant, and T is the 
# surface water temperature

longBareOut=STEFAN_B * (tair+273.15)^4  # assuming air (or maybe water) temperature is somewhere close to soil temperature... [W m-2]

net_long = longwave-longBareOut  #[W m-2]

rad=net_short+net_long  #[W m-2]

evap=(slope*rad+r_air*CP_PM*vpd/ra)/(lv*(slope+gamma*(1+(rc+rarc)/ra)))*SEC_PER_DAY  #evaporation [mm day-1]

# The units on this seem to work out except slope is per degree K and CP_PM is per degree C
# Similar problem with slope (per degree K) and gamma (per degree C) in the denominator...

# Presumably VIC folks have this fixed and the units are just mislabeled???

return(evap)
}

# loading netcdf file created from grib file from NLDAS-2 output
setwd("~/Documents/Research/LTER_EAGER/scalingC/NLDAS_forcing_files/")

library(ncdf4)

ncFiles=list.files(pattern="nc")

# set up arrays to store forcings
i=1
nc=nc_open(ncFiles[i])
summary(nc$var)
#var11-Temperature [K]
#var51-Specific humidity [kg kg-1]
#var1-Pressure [Pa]
#var33-Zonal wind speed (U) [m s-1]
#var34-Meridional wind speed (V) [m s-1]
#var205-Downward long wave radiation [W m-2]
#var153-Cloud mixing ratio [kg kg-1]
#var157-Convective Available Potential Energy [J kg-1]
#var228-Potential evaporation [kg m-2]
#var61-monthly total precipitation [kg m-2] or [mm]
#var204-downward shortwave radiation [W m-2]

shortwave=ncvar_get(nc,"var204")
nc_close(nc)

store_shortwave=array(NA,dim=c(nrow(shortwave),ncol(shortwave),length(ncFiles)))
store_longwave=store_shortwave
store_precip=store_shortwave
store_tair=store_shortwave
store_sh=store_shortwave
store_pressure=store_shortwave

store_evap=store_shortwave

for(i in 1:length(ncFiles)){
  nc=nc_open(ncFiles[i])
  
  store_shortwave[,,i]=ncvar_get(nc,"var204")
  store_longwave[,,i]=ncvar_get(nc,"var205")
  store_precip[,,i]=ncvar_get(nc,"var61")
  store_tair[,,i]=ncvar_get(nc,"var11")
  store_sh[,,i]=ncvar_get(nc,"var51")
  store_pressure[,,i]=ncvar_get(nc,"var1")
  
  nc_close(nc)
  
  # calculate monthly evap and then aggregate
  store_evap[,,i]=evap_calc(store_pressure[,,i],store_tair[,,i]-273.15,store_shortwave[,,i],store_longwave[,,i],store_sh[,,i])
}

sumPrecip=apply(store_precip,c(1,2),FUN=sum)  #kg m-2 year-1 or mm year-1
sumPrecipDaily=sumPrecip/365  # mm day-1
filled.contour(1:nrow(sumPrecip),1:ncol(sumPrecip),log10(sumPrecip))

meanEvap=-1*apply(store_evap,c(1,2),FUN=mean)  #mm day-1
filled.contour(1:nrow(meanEvap),1:ncol(meanEvap),meanEvap)

# make spatial dataframe with long, lat, and meanEvap numbers
library(sp)

long=seq(-124.9375,-67.0625,by=0.125)
lat=seq(25.0625,52.9375,by=0.125)
meanEvapDF=data.frame(long=rep(long,length(lat)),lat=rep(lat,each=length(long)),meanEvap=as.numeric(meanEvap))
coordinates(meanEvapDF) <- ~long + lat
proj4string(meanEvapDF) <- CRS("+init=epsg:4326")

# ****
# run code in usingDOENREL_evap.R, then run following code
#*****

meanEvapCONUS=meanEvapDF[!is.na(over(meanEvapDF,usaPoly)),]


colRange=colorRamp(c('blue','red'))
colorsPenman=rgb(colRange(meanEvapCONUS@data[,1]/max(c(meanEvapCONUS@data[,1],evap.krigedCONUS@data[,1])))/255)
colorsNREL=rgb(colRange(evap.krigedCONUS@data[,1]/max(c(meanEvapCONUS@data[,1],evap.krigedCONUS@data[,1])))/255)

colorsPenman_only=rgb(colRange(meanEvapCONUS@data[,1]/max(meanEvapCONUS@data[,1]))/255)
colorsNREL_only=rgb(colRange(evap.krigedCONUS@data[,1]/max(evap.krigedCONUS@data[,1]))/255)


plot(meanEvapCONUS@coords[,1],meanEvapCONUS@coords[,2],pch=15,col=colorsPenman,xlab="long. (degrees)",ylab="lat. (degrees)")
plot(evap.krigedCONUS@coords[,1],evap.krigedCONUS@coords[,2],pch=15,col=colorsNREL,xlab="long. (degrees)",ylab="lat. (degrees)")
plot(meanEvapCONUS@coords[,1],meanEvapCONUS@coords[,2],pch=15,col=colorsPenman_only,xlab="long. (degrees)",ylab="lat. (degrees)")
plot(evap.krigedCONUS@coords[,1],evap.krigedCONUS@coords[,2],pch=15,col=colorsNREL_only,xlab="long. (degrees)",ylab="lat. (degrees)")

#********
# grid designations and coordinates from nldas
#********

# the NLDAS domain extends from 25 to 53 North and -125 to -67 West
# latitude/longitude values represent center of 1/8th degree grid boxes
# Position      GridColumn   GridRow   Longitude   Latitude
# Lower left    1            1         -124.9375   25.0625
# Lower right   464          1         -67.0625    25.0625
# Upper right   464          224       -67.0625    52.9375
# Upper left    1            224       -124.9375   52.9375

setwd("..")

#nldas land versus water
nldasLvW=read.table("NLDASmask_UMDunified.asc",header=FALSE,stringsAsFactors=FALSE)
colnames(nldasLvW)=c("colNum",'rowNum','lat','long','waterVland')
# column 1: grid column number
# column 2: grid row number
# column 3: Latitude (center of 1/8th-degree grid boxes)
# column 4: Longitude (center of 1/8th-degree grid boxes)
# column 5: Mask Value (0=water, 1=land)

#nldasCONUS=read.table("NLDASmask_CONUS.asc",header=FALSE,stringsAsFactors=FALSE)
#colnames(nldasCONUS)=c("colNum",'rowNum','lat','long','CONUS')
# column 1: grid column number
# column 2: grid row number
# column 3: Latitude (center of 1/8th-degree grid boxes)
# column 4: Longitude (center of 1/8th-degree grid boxes)
# column 5: Mask Value (0=outside of CONUS, 1=CONUS)

# load NHDplus lake coordinates to get Evap for all
setwd("NHDplusOutput")
NHDlakes=read.csv("NHDplusWaterbody_LatLong.csv",header=TRUE,stringsAsFactors = FALSE)

NHDevap=data.frame(COMID=NHDlakes$COMID,LONG=NHDlakes$long,LAT=NHDlakes$lat,AREASQKM=NHDlakes$AREASQKM,FTYPE=NHDlakes$FTYPE,FCODE=NHDlakes$FCODE,EVAPpenman=NA,EVAPnrel=NA,stringsAsFactors=FALSE)

meanEvap_NLDASform=t(meanEvap)
meanEvapNREL_NLDASform=t(NRELevapCONUS)

# clean up R environment to deal with big matrices
toRemove=ls()
toRemove=toRemove[!(toRemove%in%c("meanEvap_NLDASform","meanEvapNREL_NLDASform","nldasLvW","NHDevap"))]
rm(list=toRemove)

# loop through grids of Evap
long=seq(-125,-67.125,by=0.125)
lat=seq(25,52.875,by=0.125)
lengthLong=length(long)
lengthLat=length(lat)
N=lengthLong*lengthLat
for(a in 1:length(long)){
  for(b in 1:length(lat)){
    print(((a-1)*lengthLat+b)/N*100)
    if((!is.na(meanEvap_NLDASform[b,a]) & (!is.na(meanEvapNREL_NLDASform[b,a])))){
      minlong=long[a]
      minlat=lat[b]
      maxlong=minlong+0.125
      maxlat=minlat+0.125
    
      matchLakes=which(((NHDevap$LONG>=minlong & NHDevap$LONG<maxlong) & (NHDevap$LAT>=minlat & NHDevap$LAT<maxlat)))
      if(length(matchLakes)>0){
        NHDevap$EVAPpenman[matchLakes]=meanEvap_NLDASform[b,a]
        NHDevap$EVAPnrel[matchLakes]=meanEvapNREL_NLDASform[b,a]
      }
    }
  }
}

#write.csv(NHDevap,"NHDplusWaterbody_Evaporation.csv",row.names=FALSE)






#*****
# old code for NLA2007 lakes
#*****
setwd("NLA2007")
lakes=read.csv("NLA2007_SampledLakeInformation_20091113.csv",header=TRUE,fill=TRUE,stringsAsFactors=FALSE)
lakes=lakes[lakes$VISIT_NO==1,]


library(fields)

lat=nc$dim[[2]]$vals
long=nc$dim[[1]]$vals

# interpolating daily NLDAS precip for NLA 2007 lakes
obs=list(x=long,y=lat,z=sumPrecipDaily)
loc=as.matrix(lakes[,9:10])

NLDASprecip=cbind(lakes$SITE_ID,interp.surface(obs,loc))

colnames(NLDASprecip)=c("SITE_ID","dailyNLDASprecip_mm")

#write.table(NLDASprecip,"NLA2007_nldas2000s_precip.txt",row.names=FALSE,sep="\t")


# interpolating daily NLDAS water surface evap for NLA 2007 lakes
obs=list(x=long,y=lat,z=meanEvap)
loc=as.matrix(lakes[,9:10])

NLDASevap=cbind(lakes$SITE_ID,interp.surface(obs,loc))

colnames(NLDASevap)=c("SITE_ID","dailyNLDASevap_mm")

#write.table(NLDASevap,"NLA2007_nldas2000s_evap.txt",row.names=FALSE,sep="\t")