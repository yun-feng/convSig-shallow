setwd("C:\\MutationSignature\\TCGA\\COAD\\exp")

mydata<-read.table("p.txt",sep="\t",header=F)
#mydata<-t(mydata)

d <- dist(mydata, method = "euclidean") #
fit <- hclust(d, method="ward") 
plot(fit) 
groups <- cutree(fit, k=3) 
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=3, border="red")

library(ggfortify)
autoplot(prcomp(mydata))
autoplot(kmeans(mydata, 2), data = mydata)

setwd("C:\\MutationSignature\\comp\\base3_simu_sup_CRC")
theta<-read.table("signature.txt",sep="\t",header=T,row.names=1)
theta<-as.matrix(theta)#[,c(1,7,10,11,16,19,20)])
Nsig=dim(theta)[2]
iter=20;
K=10;
mutation<-read.table("mutation.txt",sep="\t",header=F)[,1]
theta<-theta[rank(mutation),]
setwd("C:\\MutationSignature\\TCGA\\eso")

err<-c()
theta<-as.matrix(read.table("eso_relu_conv.txt",sep="\t",header=F))
#theta<-as.matrix(read.table("exp\\conv0.txt",sep="\t",header=F))
Nsig=8
for (j in 0:19){
conv<-as.matrix(read.table(paste("exp\\conv",j,".txt",sep=""),sep="\t",header=F))
for(i in 1:ncol(conv)){
  #BayesNMF
  #conv[,i]<-unlist(conv[,i]*(opp[1,]/1000))
  theta[,i]<-theta[,i]/sum(theta[,i])
conv[,i]<-conv[,i]/sum(conv[,i])


}
temp_err<-c()
for(i in 1:(Nsig)){
    #temp_err<-c(temp_err,min(apply(abs(theta[i,]-conv),2,sum))/96)
    temp_err<-rbind(temp_err,max(apply(theta[,i]*conv,2,sum)/sqrt(apply(conv^2,2,sum)*sum(theta[,i]^2))))
  }
err<-cbind(err,temp_err)
}


library(ggplot2)
err<-c(err)
sig<-rep(1:8,20)
wkdata<-data.frame(err=err,sig=factor(sig))
wkdata_relu<-wkdata
p_relu<-ggplot(data=wkdata_relu)+geom_boxplot(aes(x=sig,y=err,fill=sig))+
  geom_point(aes(x=sig,y=err))+
  labs(x="Signature",y="Cosine Similarity")
p_relu

library(corrplot)
colnames(err)<-paste("Signature",1:8,sep=" ")
row.names(err)<-paste("Signature",1:8,sep=" ")
corrplot(err, method="number")

color<-rep("black",30)
color[c(1,6,10,18)]<-"red"
color[c(5,17,20,26)]<-"green"
color[c(9,11,12,14,15,19,21,28)]<-"blue"
plot(temp_err,col=color,ylim=c(0,1))

  lego96<-as.matrix(read.table("COAD.mut.txt",sep="\t",header=F))
  num=0;
  while(num<(Nsig-3)){
    res <- BayesNMF.L1KL(as.matrix(lego96),100000,a0,tol,Kcol,Kcol,1.0)
    num=length(which(apply(res[[1]],2,sum)>0))
  }
  feature=res[[2]][which(apply(res[[1]],2,sum)>0),]
  
  temp=matrix(0*96*num,nrow=num,ncol=96)
  temp[paste("V",1:96,sep="") %in% colnames(feature)]=feature
  
  p=res[[1]][,which(apply(res[[1]],2,sum)>0)]
  write.table(p,"BayesNMF\\p.txt",sep="\t",quote=F,row.names = F,col.names =F)
  write.table(temp,"BayesNMF\\theta.txt",sep="\t",quote=F,row.names = F,col.names =F)

  err<-c()
  nmf_num<-c()
    conv<-as.matrix(read.table("BayesNMF\\theta.txt",sep="\t",header=F))
    conv<-t(conv)
    
    
    
    
    conv=conv/apply(conv,1,sum)
    temp_err<-c()
    nmf_num<-c(nmf_num,nrow(conv))
    for(i in 1:Nsig){
      #temp_err<-c(temp_err,min(apply(abs(theta[i,]-conv),1,sum))/96)
      temp_err<-rbind(temp_err,apply(theta[,i]*t(conv),2,sum)/sqrt(apply(conv^2,1,sum)*sum(theta[,i]^2)))
    }
    err<-rbind(err,temp_err)
    
    for(i in 1:Nsig){
      #temp_err<-c(temp_err,min(apply(abs(theta[i,]-conv),1,sum))/96)
      temp_err<-rbind(temp_err,apply(theta[i,]*t(theta),2,sum)/sqrt(apply(theta^2,1,sum)*sum(theta[i,]^2)))
    }
    
    
wkdata<-as.matrix(read.table("COAD.mut.txt",sep="\t",header=F))
label<-read.table("signer_sig.txt",sep="\t",header=F)
label<-paste(label[,4],paste(label[,1],label[,2],label[,3],sep=""),sep=":")
wkdata<-rbind(label,wkdata)
sample<-paste("sample",1:nrow(wkdata),sep="")
write.table(wkdata,"signer_input.txt",sep="\t",quote=F,row.names = sample,col.names =label)

wkdata<-as.matrix(read.table("C:\\MutationSignature\\comp\\base3_simu_sup_CRC\\Background.txt",sep="\t",header=F)[[1]])
wkdata<-rep(wkdata,each=3)
wkdata<-t(matrix(rep(wkdata,367),nrow=96,ncol=367))
write.table(wkdata,"signer_opp.txt",sep="\t",quote=F,row.names = F,col.names =F)


setwd("C:\\MutationSignature\\TCGA\\COAD")
library(signeR)
mut <- read.table("signer_input.txt", header=TRUE, check.names=FALSE)
opp <- read.table("signer_opp.txt")
signatures <- signeR(M=mut, Opport=opp)
signatures <- signeR(M=mut, Opport=opp, nlim=c(3,11))
BICboxplot(signatures)
SignPlot(signatures$SignExposures)


mut2 <- read.table(system.file("extdata","21_breast_cancers.mutations.txt",
                              package="signeR"), header=TRUE, check.names=FALSE)
mut <- read.table("signer_input.txt", header=TRUE, check.names=FALSE)
mut2<-mut2[,rank(colnames(mut))]
wkdata<-as.matrix(t(mut2))[,1:5]

setwd("C:\\MutationSignature\\TCGA\\COAD")
wkdata<-as.matrix(read.table("C:\\MutationSignature\\TCGA\\21bc\\exp_conv.txt",sep="\t",header=F))
label<-read.table("signer_sig.txt",sep="\t",header=F)[,4]
mut <- read.table("signer_input.txt", header=TRUE, check.names=FALSE)
row.names(wkdata)<-colnames(mut)
label<-label[order(colnames(mut))]
wkdata<-wkdata[order(colnames(mut)),]
for(i in 1:ncol(wkdata)){
  wkdata[,i]<-wkdata[,i]/sum(wkdata[,i])
}
theta<-c(wkdata)

type<-rep(row.names(wkdata),ncol(wkdata))

mut<-rep(label,ncol(wkdata))

sig<-paste("sig",rep(1:ncol(wkdata),each=96),sep="")

#relu
sig<-paste("sig",rep(c(4,5,1,2,3),each=96),sep="")
#exp
sig<-paste("sig",rep(c(4,3,2,1,5),each=96),sep="")

wkdata_relu<-data.frame(theta=theta,type=type,mut=factor(mut),sig=sig)
library(ggplot2)
p_relu<-ggplot(data=wkdata_relu)+geom_bar(aes(x=type,y=theta,fill=mut),stat="identity")+
  facet_wrap(~sig,ncol=1)+
  labs(y=NULL,x="Mutation Type")+
  ggtitle("ConvSig(ReLU)")+
  theme(axis.text.x=element_blank(),plot.title = element_text(hjust=0.5))
p_relu

wkdata_exp<-data.frame(theta=theta,type=type,mut=factor(mut),sig=sig)
p_exp<-ggplot(data=wkdata_exp)+geom_bar(aes(x=type,y=theta,fill=mut),stat="identity")+
  facet_wrap(~sig,ncol=1)+
  labs(y=NULL,x="Mutation Type")+
  ggtitle("ConvSig(exp)")+
  theme(axis.text.x=element_blank(),plot.title = element_text(hjust=0.5))
p_exp

library(cowplot)
theme_set(theme_gray())
p<-plot_grid(p_exp,p_relu,labels="AUTO",align="hv",axis="l",ncol=2,label_size = 25)
save_plot("C:\\MutationSignature\\TCGA\\21bc\\21bc.png",p,base_width=25.2,base_height = 14)

setwd("C:\\MutationSignature\\TCGA\\eso")
p<-as.matrix(read.table("relu\\P_tran.txt",sep="\t",header=F))
mut<-as.matrix(read.table("eso.mut.txt",sep="\t",header=F))
mut<-mut[20:119,]
conv<-as.matrix(read.table("relu\\conv_tran.txt",sep="\t",header=F))
syn<-p%*%t(conv)
back<-as.matrix(read.table("eso.wt.txt",sep="\t",header=F))
for(i in 1:100){
syn[i,]<-syn[i,]*back
}

temp_err<-c()
for(i in 1:100){
  #BayesNMF
  #conv[,i]<-unlist(conv[,i]*(opp[1,]/1000))
  #mut[i,]<-mut[i,]/sum(mut[i,])
  #syn[i,]<-syn[i,]/sum(syn[i,])
  temp_err<-rbind(temp_err,sum(mut[i,]*syn[i,])/sqrt(sum(mut[i,]^2)*sum(syn[i,]^2)))
  
}







