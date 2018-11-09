setwd("/data/ted/COAD")
wkdata<-read.table("colorectal.txt",sep="\t",header=F)
wkdata<-wkdata[order(wkdata[,4]),c(1,2,3,5,6)]
wkdata<-wkdata[order(wkdata[,3]),]
wkdata<-wkdata[order(wkdata[,2]),]
write.table(wkdata,"colorectal-sorted.txt",sep="\t",quote=F,row.names = F,col.names =F)
