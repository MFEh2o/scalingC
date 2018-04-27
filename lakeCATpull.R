library(curl)
url="ftp://newftp.epa.gov/EPADataCommons/ORD/NHDPlusLandscapeAttributes/LakeCat/FinalTables/"

h = new_handle(dirlistonly=TRUE)
con = curl(url, "r", h)
tbl = read.table(con, stringsAsFactors=TRUE, fill=TRUE)
close(con)
head(tbl)

urls <- paste0(url, tbl[,1])
fls = basename(urls)

setwd("~/Documents/Research/LTER_EAGER/scalingC/lakeCAT_downloads/")

for(i in 2:nrow(tbl)){
  print(paste(i,'---',fls[i]))
  curl_fetch_disk(urls[i], fls[i])
}

# careful this one is quite big
curl_fetch_disk("ftp://newftp.epa.gov/EPADataCommons/ORD/NHDPlusLandscapeAttributes/LakeCat/Framework/LkCat_Frame_min.zip","LkCat_Frame_min.zip")
