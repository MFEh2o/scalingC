# NLDAS_VIC evaporation calculations
# SEJ
# Using output from VIC run asssociated with the North American Land Data Assimilation System Phase 2 (NLDAS-2) project 
# to estimate evaporation according to the Penman equation

# REQUIRED INPUTS FOR THE LOCATION
# elevation [elev; m]
# air temperature [tair; C]
# incoming shortwave radiation [shortwave; ????]
# incoming longwave radiation [longwave]
# specific humidit [SH; Kg Kg-1]

# OUTPUT 
# evaporation [evap; mm day-1]

# PARAMETERS
LAPSE_PM = -0.006  #environmental lapse rate [C m-1]
PS_PM = 101300  #sea level air pressure [Pa]
Z0_Lower = 0.0045  #roughness
d_Lower = 0  #displacement
von_K = 0.4  #Von Karman constant for evapotranspiration
K2 = von_K*von_K
C = 2.16679 #constant for specific humidity to vapor pressure conversion [g K J-1] 
SVP_A = 0.61078  #A term in saturated vapor pressure calculation
SVP_B = 17.269  #B term in saturated vapor pressure calculation
SVP_C = 237.3  #C term in saturated vapor pressure calculation
SEC_PER_DAY = 86400  #seconds per day
H2O_SURF_ALBEDO = 0.08  #albedo of water surface
STEFAN_B = 5.6696e-8  #stefan-boltzmann constant [W/m^2/K^4]
# height of wind measuremtn is 10 m (not sure if this is needed)

# INTERMEDIATE EQUATIONS
h=287/9.81*((tair+273.15)+0.5*elev*LAPSE_PM)  #scale height in the atmosphere [????]

pz=PS_PM*exp(-elevation/h)  #surface air pressure [???]

lv=2501000-2361*tair  #latent heat of vaporization [J Kg-1]

gamma=1628.6*pz/lv  #psychrometric constant [Pa C-1]

r_air=0.003486*pz/(275+tair)  #air density [Kg m-3]

rs=0  #minimal stomatal resistance [s m-1]
rc=0  #

ra=log((2+(1/0.63-1)*d_Lower)/Z0_Lower)*log((2+(1/0.63-1)*d_Lower)/(0.1*Z0_Lower))/K2 #aerodynamic resistance [???]
  
rarc=0  #architectural resistance [s m-1]

svp=A_SVP*exp((B_SVP*tair)/(C_SVP+tair))  #saturated vapor pressure [Pa]

vp=SH*r_air*1000*(tair+273.15)/C  #vapor pressure [Pa]

vpd=svp-vp  #vapor pressure deficit [Pa]

slope=(B_SVP*C_SVP)/((C_SVP+tair)*(C_SVP+tair))*svp  #slope of saturated vapor pressure curve [Pa K-1]

net_short = (1-H2O_SURF_ALBEDO)*shortwave

# From VIC func_surf_energy_bal.c
#Tmp = Ts + KELVIN; Ts is soil temperature or surface temperature and KELVIN=273.15
#LongBareOut = STEFAN_B * Tmp * Tmp * Tmp * Tmp; 

#a lake study (Binyamin et al. 2006; Int. J. Climatol. 26: 2261-2273) used E*sigma*T^4 as outgoing and E*longwave as incoming, where
# E is emissivity of the water surface (0.97), sigma is the Stefan-Boltzman constant, and T is the 
# surface water temperature

longBareOut=STEFAN_B * (tair+273.15)^4  # assuming air (or maybe water) temperature is somewhere close to soil temperature...

net_long = longwave-longBareOut

rad=net_short+net_long

evap=(slope*rad+r_air*CP_PM*vpd/ra)/(lv*(slope+gamma*(1+(rc+rarc)/ra)))*SEC_PER_DAY  #evaporation [mm day-1]
