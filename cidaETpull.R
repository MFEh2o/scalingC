# accessing CIDA evapotranspiration data from netcdf format
library(ncdf4)
etnc=nc_open("cida_yearly_aET.nc")

et=ncvar_get(etnc)
years=2000:(2000+dim(et)[3]-1)
lat=etnc$dim[[2]]$vals
long=etnc$dim[[3]]$vals

# average across 2000's
mean_et=apply(et[,,1:10],c(1,2),mean)
cv_et=apply(et[,,1:10],c(1,2),function(x) sd(x)/mean(x))

# load list of NLA2007 lakes
lakes=read.csv("NLA2007_SampledLakeInformation_20091113.csv",header=TRUE,fill=TRUE,stringsAsFactors=FALSE)
lakes=lakes[lakes$VISIT_NO==1,]

# interpolating average ET over the 2000's for NLA 2007 lakes
library(fields)

obsMean=list(x=long,y=lat,z=mean_et)
obsCV=list(x=long,y=lat,z=cv_et)
loc=as.matrix(lakes[,9:10])

ET=cbind(lakes$SITE_ID,interp.surface(obsMean,loc),interp.surface(obsCV,loc))

colnames(ET)=c("SITE_ID","MAET_mm","MAET_CV")

#write.table(ET,"NLA2007_cida2000s_ET.txt",row.names=FALSE,sep="\t")