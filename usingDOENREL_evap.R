# Using evap calculations from meteorological (and pan evap?) sites
# data coallated into "typical meteorological years" by DOE/NREL
# found in Cavusoglu et al. (2017) Nature Communications; DOI: 10.1038/s41467-017-00581-w
rm(list=ls())

setwd("~/Documents/Research/LTER_EAGER/scalingC/DOENRELdata/")
d=read.csv("dBmTMY3.csv",header=TRUE,stringsAsFactors=FALSE)

# extract site info
sites=d[1:934,2:6]

# remove puerto rico, virgin islands, and guam
sites=sites[!(sites$STATE%in%c("GU","PR","VI")),]
d=d[!(d$STATE%in%c("GU","PR","VI")),]

# plot sites
plot(sites$LONG,sites$LAT)
library(fields)
US(add=TRUE)

# generating annual average open-water evaporation
TMYevap=sites
TMYevap$evap=NA
for(i in 1:nrow(TMYevap)){
  cur=d[d$USAF==TMYevap$USAF[i],]
  TMYevap$evap[i]=sum(cur$EV0)/365
}
range(TMYevap$evap)

# kriging TMY mean annual evap
# help from: https://rpubs.com/nabilabd/118172 & https://stackoverflow.com/questions/43436466/create-grid-in-r-for-kriging-in-gstat
library(sp)
library(rgdal)
library(raster)
library(gstat)
library(rgeos)

# pull just lat-long and evap info
evap2krig=TMYevap[,4:6]

# convert to spatialpointsdataframe
coordinates(evap2krig) <- ~LONG + LAT
# set projection then convert to utm for kriging
proj4string(evap2krig) <- CRS("+init=epsg:4326")
evap2krigUTM=spTransform(evap2krig,CRS("+proj=utm"))

# define grid - using NLDAS grid
# set the origin
ori <- SpatialPoints(cbind(-124.9375,25.0625),proj4string=CRS("+init=epsg:4326"))

# define how many cells for x and y axis
x_cell <- 464
y_cell <- 224
# define the resolution to be 1/8th degree
cell_size <- 1/8
# create the extent
ext <- extent(coordinates(ori)[1,1],coordinates(ori)[1,1]+(x_cell*cell_size),coordinates(ori)[1,2],coordinates(ori)[1,2]+(y_cell*cell_size))
# convert extent to raster
ras<-raster(ext)
# set resolution of raster and set values to 0
res(ras)<-c(cell_size,cell_size)
ras[]<-0
# project raster
projection(ras)<-CRS("+init=epsg:4326")

# convert raster to spatial pixel
st_grid<-rasterToPoints(ras,spatial=TRUE)
gridded(st_grid)<-TRUE

# convert to utm for kriging
st_gridUTM=spTransform(st_grid,CRS("+proj=utm"))

# fit variogram
evap.vgm<-variogram(evap~1,evap2krig)
#show.vgms()
evap.fit<-fit.variogram(evap.vgm,model=vgm(1,"Lin",0))
plot(evap.vgm,evap.fit)  # not great at short distance, but ok

evap.kriged2<-krige(evap ~ 1, evap2krigUTM,st_gridUTM,model=evap.fit)

# put back into long-lat for use
evap.kriged2proj=spTransform(evap.kriged2,CRS("+init=epsg:4326"))

#spplot(evap.kriged2proj,"var1.pred")

# alternative method for kriging
#library(automap)
#evap.kriged<-autoKrige(evap~1,evap2krigUTM,st_gridUTM)

# make matrix of evap for map plot - need to round coords because conversion from UTM creates small differences???
evap.kriged2proj@coords[,1]=round(evap.kriged2proj@coords[,1],3)
evap.kriged2proj@coords[,2]=round(evap.kriged2proj@coords[,2],3)
NRELevap=matrix(NA,x_cell,y_cell)
longs=sort(unique(evap.kriged2proj@coords[,1]))
lats=sort(unique(evap.kriged2proj@coords[,2]))
for(i in length(lats):1){
  cur=evap.kriged2proj[evap.kriged2proj@coords[,2]==lats[i],]
  NRELevap[,i]=cur@data$var1.pred[order(cur@coords[,1],decreasing=FALSE)]
}

filled.contour(longs,lats,NRELevap)

#**********
# mask evap.kriged2proj with CONUS polygon
#**********
# get polygon for masking grid
library(ggplot2)

usa=map_data("usa")
usa$id=usa$group

usaSPDF=usa
coordinates(usaSPDF) <- ~long + lat

# function to convert points to a polygon
points2polygons <- function(df,data) {
  get.grpPoly <- function(group,ID,df) {
    Polygon(coordinates(df[df$id==ID & df$group==group,]))
  }
  get.spPoly  <- function(ID,df) {
    Polygons(lapply(unique(df[df$id==ID,]$group),get.grpPoly,ID,df),ID)
  }
  spPolygons  <- SpatialPolygons(lapply(unique(df$id),get.spPoly,df))
  SpatialPolygonsDataFrame(spPolygons,match.ID=T,data=data)
}

usaPoly=points2polygons(usaSPDF,data.frame(data=rep(0,length(unique(usaSPDF$id)))))
proj4string(usaPoly) <- CRS("+init=epsg:4326")
plot(usaPoly)


# actual filtering of points in evap.kriged2proj
class(usaPoly) <- "SpatialPolygons"
evap.krigedCONUS=evap.kriged2proj[!is.na(over(evap.kriged2proj,usaPoly)),]

# make matrix of evapCONUS for map plot - need to round coords because conversion from UTM creates small differences???
evap.krigedCONUS@coords[,1]=round(evap.krigedCONUS@coords[,1],3)
evap.krigedCONUS@coords[,2]=round(evap.krigedCONUS@coords[,2],3)
NRELevapCONUS=matrix(NA,x_cell,y_cell)

for(i in length(lats):1){
  cur=evap.krigedCONUS[evap.krigedCONUS@coords[,2]==lats[i],]
  if(nrow(cur)>0){
    NRELevapCONUS[match(cur@coords[,1],longs),i]=cur@data$var1.pred
  }
}

filled.contour(longs,lats,NRELevapCONUS)


# plot
meanEvapDF=data.frame(long=evap.krigedCONUS@coords[,1],lat=evap.krigedCONUS@coords[,2],meanEvap=evap.krigedCONUS@data[,1])

colRange=colorRamp(c('blue','red'))
colors=rgb(colRange(meanEvapDF$meanEvap/max(meanEvapDF$meanEvap))/255)
plot(meanEvapDF$long,meanEvapDF$lat,pch=15,col=colors)
