## Using Carly's GRFP Model (rates v. fates + Diehl pelagic GPP, but has sediment P state variable) as a starting place
## biggest issue currently is that it assumes well-mixed water column
## also ignores labile v. recalcitrant C
rm(list=ls())

library(deSolve)

# define the model

#state variables
# C - terrestrial dissolved organic carbon [g C]
# A - phytoplankton phosphorus [g P]
# P - pelagic bioavailable phosphorus [g P]
# S - sediment bioavailable phosphorus [g P]

# parameters
#kbg - contribution of water to light attenuation coefficient [?]; -0.05 this needs to be checked with Patrick
#kc - DOC contribution to light attenuation coefficient [?]; 0.00000042
#ka - phytoplankton contribution to light attenuation coefficient [?]; 0.000015
#I0 - incoming light [uE]; 300
#Cprecip - concentration of DOC in precipitation [g C m-3]; 1
#Pprecip - concentration of phosphorus in precipitation [g P m-3]; 0.005
#d - decay rate of DOC [day-1]; 0.01
#pa - max growth rate of phytoplankton [day-1]; 1
#ha - half saturation constant for light for phytoplankton [uE]; 100
#ma - half saturation constant for bioavailable phosphorus [g P m-3]; 0.003
#la - phytoplankton mortality rate [day-1]; 0.1
#r - settling rate of phytoplankton [m day-1]; 0.1
#ap - sediment phosphorus release rate [m day-1]; 0.05
#b - rate of permenant burial of sediment phosphorus [day-1]; 0.001
#zsed - depth of the active layer of sediments [m]; 0.01

#Cin - concentration of DOC in stream [g C m-3]; 10 - need to make a function of catchment landcover
#Pin - concentration of phospohrus in stream [g P m-3]; 0.040 - need to make a function of catchment landcover
#Al - lake area [m2]; from NLA
#Ac - catchment area [m2]; from NLA
#zmax - lake max and mean depth (assuming cylindrical lake) [m]; from NLA, but need to make f(Al)
#map - mean annual precipitation for the lake catchment [m yr-1]; from NLDAS for 2000s
#ec - evapotranspiration rate of the lake catchment [m day-1]; from cida for 2000s
#el - evaporation rate from the lake; [m day-1]; from NLDAS for 2000s

tstep=function(t, y, parms) {
  with(as.list(c(y,parms)),{
  
    ###### Intermediate calculations
    V=Al*zmax  # lake volume [m3]
    Vsed=Al*zsed  # volume of active lake sediments [m3]
    precip=map/365  # daily precipitation [m d-1]
    Qin=(precip - ec)*Ac  # watershed hydrologic input [m3 d-1]
    Qprecip=precip*Al  # water input from direct precipitation [m3 d-1]
    Qout=Qin + (precip-el)*Al  # outflow through lake outlet [m3 d-1]
    Qevap=el*Al  # water loss via evaporation [m3 d-1]
    
    kd=kbg + (kc*C) + (ka*A)  # light attenuation coefficient of lake water [m-1] 
    Izmax=I0*exp(-kd*zmax)  # light available at bottom of the lake [uE]
  
    ###### Differential equations
    dC=(Qin*Cin) - ((C/V)*Qout) + (Qprecip*Cprecip) - d*C  # terrestrial dissolved organic carbon pool [g C]
    dA=(((pa/(kd*zmax))*log((ha+I0)/(ha+Izmax))*((P/V)/((P/V)+ma))) - la - (r/zmax) - (Qout/V))*A  # phytoplankton pool [g P]
    dP=(Qin*Pin) - ((P/V)*Qout) + (Qprecip*Pprecip) - ((pa/(kd*zmax))*log((ha+I0)/(ha+Izmax))*((P/V)/((P/V)+ma)))*A + (la*A) + ((ap/zmax)*((S/Vsed)-(P/V)))*V  # dissolved bioavailable phosphorus pool [g P]
    dS=(r/zmax)*A - ((ap/zmax)*((S/Vsed)-(P/V)))*V - b*S  # bioavailable sediment phosphorus pool [g P]
    
    return(list(c(dC=dC, dA=dA, dP=dP, dS=dS)))
  })
}

###### load forcings for EPA National Lake Assessment 2007 lakes --> IGNORING LANDUSE VARIATION FOR NOW (i.e. Cin and Pin are constant across lakes)
# load NLA 2007 files
lakes=read.csv("NLA2007_SampledLakeInformation_20091113.csv",header=TRUE,fill=TRUE,stringsAsFactors=FALSE)
lakes=lakes[lakes$VISIT_NO==1,]
basins=read.csv("nla2007_basin_landuse_metrics_20061022.csv",header=TRUE,fill=TRUE,stringsAsFactors=FALSE)

# sort NLA files to ensure rows match
basins=basins[order(basins[,1]),]
lakes=lakes[order(lakes[,1]),]
dim(basins)
dim(lakes)
sum(basins[,1]==lakes[,1])

# pull columns wanted for simulations
toSim=cbind(lakes[,c("SITE_ID","LON_DD","LAT_DD","LAKE_ORIGIN","HUC_8","DEPTHMAX","LAKEAREA")],basins[,c("BASINAREA_KM2")])
colnames(toSim)[6:8]=c("DEPTHMAX_m","LAKEAREA_KM2","BASINAREA_KM2")

# load forcings we've generated
cidaET=read.table("NLA2007_cida2000s_ET.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
cidaET=cidaET[order(cidaET[,1]),]
PRISM=read.table("NLA2007_prism2000s_MAP.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
PRISM=PRISM[order(PRISM[,1]),]
NLDASprecip=read.table("NLA2007_nldas2000s_precip.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
NLDASprecip=NLDASprecip[order(NLDASprecip[,1]),]
NLDASevap=read.table("NLA2007_nldas2000s_evap.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
NLDASevap=NLDASevap[order(NLDASevap[,1]),]

# quick compare PRISM and NLDAS precip --> r=0.97
plot(PRISM$MAP_mm,NLDASprecip$dailyNLDASprecip_mm*365,xlab="PRISM mean annual precipitation (mm)",ylab="NLDAS mean annual precipitation (mm)")
abline(a=0,b=1,lwd=3,col='red')
cor(PRISM$MAP_mm[!is.na(NLDASprecip$dailyNLDASprecip_mm)],NLDASprecip$dailyNLDASprecip_mm[!is.na(NLDASprecip$dailyNLDASprecip_mm)],)

# go with NLDAS precip for consistency with evap.
toSim=cbind(toSim,NLDASprecip$dailyNLDASprecip_mm,cidaET$MAET_mm/365,NLDASevap$dailyNLDASevap_mm)
colnames(toSim)[9:11]=c('dailyNLDASprecip_mm','dailyCIDAet_mm','dailyNLDASevap_mm')

# remove the lake with NA for precip and evap --> FIX THIS IN NLDAS PULL!!!
toSim=toSim[-is.na(toSim$dailyNLDASprecip_mm),]


###### simulate for EPA NLA 2007 lakes
NLA2007_modelEquil=matrix(NA,nrow(toSim),4)
times=1:1000

for(i in 1:nrow(NLA2007_modelEquil)){
  print(i)
  # defining initial states, parameters, and time over which to solve
  x0=c(C=1*10000*2, A=0.005*10000*2, P=0.005*10000*2, S=1)
  # hacky fix to convert depth at sampling to mean --> dividing by 5???
  parms=c(kbg=-0.05, kc=0.00000042, ka=0.000015, I0=300, zmax=toSim[i,6]/5, Cin=10, Cprecip=1, d=0.01, pa=1, ha=100,
          ma=0.003, la=0.1, r=0.1, Pin=0.040, Pprecip=0.005, ap=0.05, b=0.001, map=toSim[i,9]*365/1000, 
          ec=toSim[i,10]/1000, Al=toSim[i,7]*1e6, el=toSim[i,11]/1000, zsed=0.01, Ac=toSim[i,8]*1e6)
  
  # simulate model
  run=ode(y = x0, times = times, func = tstep, parms = parms, method = "lsoda")

  # store equilibrium values
  NLA2007_modelEquil[i,]=run[nrow(run),2:ncol(run)]
}

# Getting DLSODA warnings for a number of lakes - check for hydrology budget working out for lakes...

