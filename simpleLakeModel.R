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

#Al - lake area [m2]; 1e4
#Ac - catchment area [m2]; 1e6
#zmax - lake max and mean depth (assuming cylindrical lake) [m]; 2
#zsed - depth of the active layer of sediments [m]; 0.01
#Cin - concentration of DOC in stream [g C m-3]; 10
#Pin - concentration of phospohrus in stream [g P m-3]; 0.040
#map - mean annual precipitation for the lake catchment [m yr-1]; 0.8
#ec - evapotranspiration rate of the lake catchment [m day-1]; 0.0014
#el - evaporation rate from the lake; [m day-1]; 0.0025

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

###### defining initial states, parameters, and time over which to solve
x0=c(C=1*10000*2, A=0.005*10000*2, P=0.005*10000*2, S=1)
parms=c(kbg=-0.05, kc=0.00000042, ka=0.000015, I0=300, zmax=2, Cin=10, Cprecip=1, d=0.01, pa=1, ha=100,
        ma=0.003, la=0.1, r=0.1, Pin=0.040, Pprecip=0.005, ap=0.05, b=0.001, map=0.8, 
        ec=0.0014, Al=10000, el=0.0025, zsed=0.01, Ac=1000000)
times=1:250

###### simulate and plot results
run=ode(y = x0, times = times, func = tstep, parms = parms, method = "lsoda")

dev.new()
par(mfrow=c(2,2))
plot(run[,1], run[,3], type = 'l', lwd=3,col='darkgreen',xlab='Time (days)', ylab='phytoplankton biomass (g P)')
plot(run[,1],run[,2],type='l',lwd=3,col='brown',xlab='Time (days)',ylab='terrestrial DOC (g C)')
plot(run[,1],run[,4],type='l',lwd=3,col='blue',xlab='Time (days)',ylab='pelagic phosphorus (g P)')
plot(run[,1],run[,5],type='l',lwd=3,col='black',xlab='Time (days)',ylab='sediment phosphorus (g P)')