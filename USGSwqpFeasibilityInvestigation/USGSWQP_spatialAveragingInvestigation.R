### evaluating spatial averaging of USGS waterchem data for 2007 NLA sites
### 12-9-16
### SEJ

##***************** TAKEHOME****************##
#### seems like this isn't great coverage for the NLA lakes...
### something like 80% of lakes have USGS measurements within ~50 km
### not terrible, but not super close either...
### more like 25-30% within a HUC8

### will start with landcover models first, as these are simpler anyway



### tabulating how many stream sites are available for each lake within the surrounding HUC8 and boxes of varied size

rm(list=ls())
library(dataRetrieval)

# load NLA 2007 lakes
lakes=read.csv("../NLA2007_SampledLakeInformation_20091113.csv",stringsAsFactors=FALSE)
lakes=lakes[lakes$VISIT_NO==1,]
basin=read.csv("../nla2007_basin_landuse_metrics_20061022.csv",stringsAsFactors=FALSE)

# functions for counting sites in box around lat,lon with given dimensions with provided parameter code
countSites<-function(coords,buffer,param){
	# coords in vector longitude, latitude
	# buffer in decimal degrees (1 degree ~100 km on a side)
	# parameter code from USGS database: 00666 -> TP, 00681 -> DOC
	
	boxCoords=round(c(coords[1]-buffer,coords[2]-buffer,coords[1]+buffer,coords[2]+buffer),6)
	
	sites=try(whatNWISsites(bBox=boxCoords,parameterCd=param),silent=TRUE)
	
	if(class(sites)=="try-error"){
		return(0)
	}else{
		return(sum(sites$site_tp_cd=="ST"))
	}
}


# pre-define matrix for storing site counts
TPsiteCounts=matrix(NA,nrow(lakes),6)
DOCsiteCounts=matrix(NA,nrow(lakes),6)

# loop through lakes and count sites
for(i in 1:nrow(lakes)){
	print(i)
	curCoords=as.numeric(lakes[i,9:10])
	curHUC8=lakes$HUC_8[i]

	TPsiteCounts[i,2]=countSites(curCoords,0.125,"00666")
	TPsiteCounts[i,3]=countSites(curCoords,0.25,"00666")
	TPsiteCounts[i,4]=countSites(curCoords,0.5,"00666")
	TPsiteCounts[i,5]=countSites(curCoords,1,"00666")
	hucTP=try(whatNWISsites(huc=curHUC8,parameterCd="00666"),silent=TRUE)
	if(class(hucTP)=="try-error"){
		TPsiteCounts[i,6]=0
	}else{
		TPsiteCounts[i,6]=sum(hucTP$site_tp_cd=="ST")
	}
	TPsiteCounts[i,1]=lakes[i,1]
	
	DOCsiteCounts[i,2]=countSites(curCoords,0.125,"00681")
	DOCsiteCounts[i,3]=countSites(curCoords,0.25,"00681")
	DOCsiteCounts[i,4]=countSites(curCoords,0.5,"00681")
	DOCsiteCounts[i,5]=countSites(curCoords,1,"00681")
	hucDOC=try(whatNWISsites(huc=curHUC8,parameterCd="00681"),silent=TRUE)
	if(class(hucDOC)=="try-error"){
		DOCsiteCounts[i,6]=0
	}else{
		DOCsiteCounts[i,6]=sum(hucDOC$site_tp_cd=="ST")
	}
	DOCsiteCounts[i,1]=lakes[i,1]
}

colnames(TPsiteCounts)=c('NLAid','squareEighth','squareQuarter','squareHalf','squareOne','HUC8')
colnames(DOCsiteCounts)=colnames(TPsiteCounts)

#write.table(TPsiteCounts,"TPsiteCounts.txt",row.names=FALSE,sep="\t")
#write.table(DOCsiteCounts,"DOCsiteCounts.txt",row.names=FALSE,sep="\t")

# plot site count results
TPsiteCounts=read.table("TPsiteCounts.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
DOCsiteCounts=read.table("DOCsiteCounts.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)

dev.new()
par(mfrow=c(3,2))
hist(TPsiteCounts[,2],breaks=seq(0,max(TPsiteCounts[,-1])),xlim=c(0,10))
hist(TPsiteCounts[,3],breaks=seq(0,max(TPsiteCounts[,-1])),xlim=c(0,10))
hist(TPsiteCounts[,4],breaks=seq(0,max(TPsiteCounts[,-1])),xlim=c(0,10))
hist(TPsiteCounts[,5],breaks=seq(0,max(TPsiteCounts[,-1])),xlim=c(0,10))
hist(TPsiteCounts[,6],breaks=seq(0,max(TPsiteCounts[,-1])),xlim=c(0,10))

colSums((TPsiteCounts[,-1]>0)*1)/nrow(TPsiteCounts)

dev.new()
par(mfrow=c(3,2))
hist(DOCsiteCounts[,2],breaks=seq(0,max(DOCsiteCounts[,-1])),xlim=c(0,10))
hist(DOCsiteCounts[,3],breaks=seq(0,max(DOCsiteCounts[,-1])),xlim=c(0,10))
hist(DOCsiteCounts[,4],breaks=seq(0,max(DOCsiteCounts[,-1])),xlim=c(0,10))
hist(DOCsiteCounts[,5],breaks=seq(0,max(DOCsiteCounts[,-1])),xlim=c(0,10))
hist(DOCsiteCounts[,6],breaks=seq(0,max(DOCsiteCounts[,-1])),xlim=c(0,10))

colSums((DOCsiteCounts[,-1]>0)*1)/nrow(DOCsiteCounts)
