# PRISM data download and analysis for NE GLEON sites
# SEJ
rm(list=ls())

library(prism)

options(prism.path="./prismDownloads")

# pull annual precipitation from prism database
get_prism_annual(type="ppt",years=2000:2009,keepZip=FALSE)

# look at the files we have
ls_prism_data()

# generate list of precipitation data files to pull numbers based on lat-long
PPTto_slice=ls_prism_data()[,1]

# load list of NLA2007 lakes
lakes=read.csv("NLA2007_SampledLakeInformation_20091113.csv",header=TRUE,fill=TRUE,stringsAsFactors=FALSE)
lakes=lakes[lakes$VISIT_NO==1,]

annualPPT=matrix(NA,length(PPTto_slice),nrow(lakes))
for(i in 1:ncol(annualPPT)){
  print(i)
	curCoords=as.numeric(lakes[i,9:10])
	annualPPT[,i]=prism_slice(curCoords,ls_prism_data()[,1])$data[,1]
}

MAP=cbind(lakes$SITE_ID,colMeans(annualPPT),apply(annualPPT,2,function(x) sd(x)/mean(x)))
colnames(MAP)=c("SITE_ID","MAP_cm","MAP_CV")

#write.table(MAP,"NLA2007_prism2000s_MAP.txt",row.names=FALSE,sep="\t")