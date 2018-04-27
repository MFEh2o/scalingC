# Using NHD hydrography snapshots to get lat-long for lake centroids and associate with COMID
rm(list=ls())
setwd("~/Documents/Research/LTER_EAGER/scalingC/NHDplusDownloads/")

library(rgdal)
library(rgeos)
#library(raster)
#library(geosphere)

regions=list.files()
regions=regions[-1]

i=1
setwd(paste(regions[i],"/NHDSnapshot/Hydrography",sep=""))
shape<-readOGR(dsn=".",layer="NHDWaterbody")
shapeDF=as.data.frame(shape)
shapeCent=gCentroid(shape,byid=TRUE,id=shapeDF$COMID)
out=data.frame(COMID=rownames(shapeCent@coords),long=shapeCent@coords[,1],lat=shapeCent@coords[,2],AREASQKM=shapeDF$AREASQKM,FTYPE=shapeDF$FTYPE,FCODE=shapeDF$FCODE)
setwd("~/Documents/Research/LTER_EAGER/scalingC/NHDplusDownloads/")

for(i in 2:length(regions)){
  setwd(paste(regions[i],"/NHDSnapshot/Hydrography",sep=""))

  shape<-readOGR(dsn=".",layer="NHDWaterbody")
  shapeDF=as.data.frame(shape)
  colnames(shapeDF)=toupper(colnames(shapeDF))
  shapeCent=gCentroid(shape,byid=TRUE,id=shapeDF$COMID)
  #dim(shapeCent@coords)
  #head(shapeCent@coords)
  app=data.frame(COMID=rownames(shapeCent@coords),long=shapeCent@coords[,1],lat=shapeCent@coords[,2],AREASQKM=shapeDF$AREASQKM,FTYPE=shapeDF$FTYPE,FCODE=shapeDF$FCODE)
  out=rbind(out,app)
  setwd("~/Documents/Research/LTER_EAGER/scalingC/NHDplusDownloads/")
}

setwd("~/Documents/Research/LTER_EAGER/scalingC")
#write.csv(out,"NHDplusWaterbody_LatLong.csv",row.names=FALSE)

#### FCODE definitions
# Feature Type FCode Description
# ESTUARY 49300 feature type only: no attributes
# ICE MASS 37800 feature type only: no attributes
# LAKE/POND 39000 feature type only: no attributes
# LAKE/POND 39001 Hydrographic Category|intermittent
# LAKE/POND 39004 Hydrographic Category|perennial
# LAKE/POND 39005 Hydrographic Category|intermittent; Stage|high water elevation
# LAKE/POND 39006 Hydrographic Category|intermittent; Stage|date of photography
# LAKE/POND 39009 Hydrographic Category|perennial; Stage|average water elevation
# LAKE/POND 39010 Hydrographic Category|perennial; Stage|normal pool
# LAKE/POND 39011 Hydrographic Category|perennial; Stage|date of photography
# LAKE/POND 39012 Hydrographic Category|perennial; Stage|spillway elevation
# PLAYA 36100 feature type only: no attributes
# RESERVOIR 43600 feature type only: no attributes
# RESERVOIR 43601 Reservoir Type|aquaculture
# RESERVOIR 43603 Reservoir Type|decorative pool
# RESERVOIR 43604 Reservoir Type|disposal-tailings pond; Construction Material|earthen
# RESERVOIR 43605 Reservoir Type|disposal-tailings pond
# RESERVOIR 43606 Reservoir Type|disposal-unspecified
# RESERVOIR 43607 Reservoir Type|evaporator
# RESERVOIR 43608 Reservoir Type|swimming pool
# RESERVOIR 43609 Reservoir Type|treatment-cooling pond
# RESERVOIR 43610 Reservoir Type|treatment-filtration pond
# RESERVOIR 43611 Reservoir Type|treatment-settling pond
# RESERVOIR 43612 Reservoir Type|treatment-sewage treatment pond
# RESERVOIR 43613 Reservoir Type|water storage; Construction Material|nonearthen
# RESERVOIR 43614 Reservoir Type|water storage; Construction Material|earthen; Hydrographic Category|intermittent
# RESERVOIR 43615 Reservoir Type|water storage; Construction Material|earthen; Hydrographic Category|perennial
# RESERVOIR 43617 Reservoir Type|water storage
# RESERVOIR 43618 Reservoir Type|unspecified; Construction Material|earthen
# RESERVOIR 43619 Reservoir Type|unspecified; Construction Material|nonearthen
# RESERVOIR 43621 Reservoir Type|water storage; Hydrographic category|perennial
# RESERVOIR 43623 Reservoir Type|evaporator; Construction Material|earthen
# RESERVOIR 43624 Reservoir Type|treatment
# RESERVOIR 43625 Reservoir Type|disposal
# RESERVOIR 43626 Reservoir Type|disposal; Construction Material|nonearthen
# SWAMP/MARSH 46600 feature type only: no attributes
# SWAMP/MARSH 46601 Hydrographic category|intermittent
#SWAMP/MARSH 46602 Hydrographic category|perennial