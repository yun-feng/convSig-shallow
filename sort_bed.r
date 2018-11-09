wkdata<-read.table("/data/ted/multi/Ubiquitous_s.bed");
wkdata<-wkdata[order(wkdata[,3]),]
wkdata<-wkdata[order(wkdata[,2]),]
wkdata<-wkdata[order(wkdata[,1]),]
write.table(wkdata,"/data/ted/multi/Ubiquitous-sorted.bed",quote=F,row.names=F,sep="\t");
