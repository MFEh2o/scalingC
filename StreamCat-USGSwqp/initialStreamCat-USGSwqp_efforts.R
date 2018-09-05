rm(list=ls())

library(dataRetrieval)

library(hydrolinks)



# coords in vector longitude, latitude
# buffer in decimal degrees (1 degree ~100 km on a side)
# parameter code from USGS database: 00666 -> TP, 00681 -> DOC

sites=try(whatNWISsites(bBox=boxCoords,parameterCd=param),silent=TRUE)


WIdischarge_sites=whatNWISsites(stateCd="WI",parameterCd="00060")
WItp_sites=whatNWISsites(stateCd="WI",parameterCd="00666")
WIdoc_sites=whatNWISsites(stateCd="WI",parameterCd="00681")

WI=list(q=WIdischarge_sites,tp=WItp_sites,doc=WIdoc_sites)

sum(WI$q$site_no%in%WI$tp$site_no)
sum(WI$q$site_no%in%WI$doc$site_no)
sum(WI$tp$site_no%in%WI$doc$site_no)

WIoverlap=WI[[which(unlist(lapply(WI,nrow))==max(unlist(lapply(WI,nrow))))]]$site_no
WIoverlap=WIoverlap[WIoverlap%in%WI[[1]]$site_no]
WIoverlap=WIoverlap[WIoverlap%in%WI[[2]]$site_no]
WIoverlap=WIoverlap[WIoverlap%in%WI[[3]]$site_no]

WIoverlapSites=WI[[which(unlist(lapply(WI,nrow))==max(unlist(lapply(WI,nrow))))]]
WIoverlapSites=WIoverlapSites[WIoverlapSites$site_no%in%WIoverlap,]


temp=link_to_flowlines(WIoverlapSites$dec_lat_va[1],WIoverlapSites$dec_long_va[1],WIoverlapSites$site_no[1])


# pulling StreamCat data
setwd("~/Documents/Research/LTER_EAGER/scalingC/StreamCat-USGSwqp/StreamCat_Data/")

library(curl)
url="ftp://newftp.epa.gov/EPADataCommons/ORD/NHDPlusLandscapeAttributes/StreamCat/HydroRegions/"

h=new_handle(dirlistonly=TRUE)
con=curl(url,"r",h)
tbl=read.table(con,stringsAsFactors=FALSE,fill=TRUE)
close(con)
head(tbl)

urls = paste0(url,grep("Region05",tbl[,1],value=TRUE))
fls=basename(urls)
for(i in 1:4){
  curl_fetch_disk(urls[i],fls[i])
}

