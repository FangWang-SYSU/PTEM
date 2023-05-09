workpath="./PTEM"
setwd(workpath)
library(diptest)
library(LaplacesDemon)
source("model.function.R")
model <- c("Neutral","Branching","Linear","Punctuated")

#####Figure 1B. The probability density distribution of CNA frequency corresponding to each evolution pattern
model <- c("Neutral","Branching","Linear","Punctuated")
CNVfre=list()
k=1
for (m in model){
  CNVdata = read.csv(paste("./simulation/CNVdensity/",m,".CNV.Fre.txt",sep=""),sep="\t")
  fre <- CNVdata$freq
  CNVfre[[k]]=fre
  k=k+1
}
pdf(file = "./simulation/CNVdensity/modality.distribution.pdf",width=6,height = 6)
par(mfrow=c(2,2))
plot(density(CNVfre[[1]],width=0.3),main=model[1],lwd=2,xlim=c(0,1),xlab="Frequency of CNAs")
plot(density(CNVfre[[2]],width = 0.2),main=model[2],lwd=2,xlim=c(0,1),xlab="Frequency of CNAs")
plot(density(CNVfre[[3]]),main=model[3],lwd=2,xlim=c(0,1),xlab="Frequency of CNAs")
plot(density(CNVfre[[4]]),main=model[4],lwd=2,xlim=c(0,1),xlab="Frequency of CNAs")
dev.off()


#######Figure2A. The heatmap of simulation single cell profile
for (m in model){
  CNVs <- read.table(paste("./simulation/CNVprofile/",m,".CNVs.N0.txt",sep=""),sep="\t",header = T)
  breaks1 <- seq(0,2,by = 0.2)
  breaks2 <- seq(2,max(CNVs),by=0.2)
  breaks <- c(breaks1,breaks2)
  gradient1 <- colorpanel(length(breaks1),"darkblue","blue","white")
  gradient2 <- colorpanel(length(breaks2),"white","red","darkred")
  hm.col <- c(gradient1,gradient2[2:length(gradient2)])
  pdf(file=paste("./simulation/CNVprofile/",m,".heatmap.N0.pdf",sep=""),width = 5,height = 5,useDingbats = F)
  pheatmap(CNVs,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,col=hm.col)
  dev.off()
}



#######Estimating performance for distinguish neutral and selective model
####Figure 2B and C: power on trainning set
library(ROCR)
modalityres <- c()
for (m in model){
  CNVdata = read.csv(paste("./simulation/trainningset/N0_",m,"_Freq.txt",sep=""),sep="\t")
  dataCount <- unique(CNVdata$Nth_simu)
  alpha <- sapply(dataCount, function(i,CNVdata){
    subdata <- CNVdata[CNVdata$Nth_simu==i&CNVdata$freq>0.05,]
    fre <- subdata$freq
    modalityEst <- dip.test(fre)
    pvalue <- modalityEst$p.value
    lable <- is.unimodal(fre)
    return(c(pvalue,lable))
  },CNVdata)
  modalityres <- rbind(modalityres,data.frame(pvalue = alpha[1,],lable = alpha[2,],model = m))
}

modalityres$evolution <- modalityres$model
modalityres$evolution[modalityres$evolution!="Neutral"]="Selection"
pred <- prediction(1-modalityres$pvalue, modalityres$evolution)
perfB <- performance(pred,"tpr","fpr")
aucB  <- performance(pred, 'auc')
aucB  <- unlist(slot(aucB,"y.values"));
pdf(file = "./simulation/trainningset/neutra.ROC.pdf",width=8,height = 5)
par(mfrow=c(1,2))
boxplot(modalityres$pvalue[modalityres$evolution=="Neutral"],modalityres$pvalue[modalityres$evolution=="Selection"],ylab="p.value",axes=F)
axis(side=1,at=c(1,2),labels=c("Neutral","Selection"))
axis(side=2)
abline(h=0.05,lty=2)
box()
plot(perfB,lwd=2,main = paste("Neutral AUC = ",round(aucB,3),sep=""))  
dev.off()
####Figure 2D：The effect of noise
noise <- c(0:10)
res <- c()
for (m in model){
  CNVdata = read.csv(paste("./simulation/noise/",m,".CNV.txt",sep=""),sep="\t")
  for (n in noise){
    subdata=CNVdata[CNVdata$noise==n,]
    dataCount <- unique(subdata$Nth_simu)
    alpha <- sapply(dataCount, function(i,subdata){
      subdata1 <- subdata[subdata$Nth_simu==i&subdata$freq>0.05,]
      fre <- subdata1$freq
      modalityEst <- dip.test(fre)
      pvalue <- modalityEst$p.value
      lable <- is.unimodal(fre)
      return(c(pvalue,lable))
    },subdata)
    res <- rbind(res,data.frame(pvalue = alpha[1,],lable = alpha[2,],noise = n,model = m)) 
  }
}
res$evolution <- res$model
res$evolution[res$evolution!="Neutral"]="Selection"
AUC <- c()
for (n in noise){
  subres=res[res$noise == n,]
  pred <- prediction(1-subres$pvalue, subres$evolution)
  perfB <- performance(pred,"tpr","fpr")
  aucB  <- performance(pred, 'auc')
  aucB  <- unlist(slot(aucB,"y.values"));
  AUC <- rbind(AUC,data.frame(noise = n, AUC = aucB))
}
pdf(file = "./simulation/noise/Neutral.nosise.AUC.pdf",width = 5,height = 5)
plot(AUC$AUC,type="o",pch=20,cex=2,ylim=c(0,1),lwd=2,col=1,xlim=c(0,11),xlab="Noise",ylab="AUC")
dev.off()

#########Prediction of selective evolution pattern and assessment
#####Figure 2E: density plot of corrected diversity score 
res <- c()
for (m in model){
  CNVdata = read.csv(paste("./simulation/trainningset/N0_",m,"_Freq.txt",sep=""),sep="\t")
  dataCount <- unique(CNVdata$Nth_simu)
  cds <- sapply(dataCount, function(i,CNVdata){
    cds <- selectModel(CNVdata$freq[CNVdata$Nth_simu==i])
    return(cds)
  },CNVdata)
  res <- rbind(res,data.frame(cds = cds,model = m))
}
branchcds = res$cds[res$model=='Branching']
linearcds = res$cds[res$model=='Linear']
punctuatedcds = res$cds[res$model=='Punctuated']
pdf(file="./simulation/trainningset/CDS.density.pdf",width = 4,height = 4)
plot(density(branchcds),xlim=c(0,1.2),ylim=c(0,7),lwd=2,col=1,xlab="diversity",main="")
par(new=T)
plot(density(linearcds),xlim=c(0,1.2),ylim=c(0,7),lwd=2,col=2,xlab="",ylab="",main="")
par(new=T)
plot(density(punctuatedcds),xlim=c(0,1.2),ylim=c(0,7),lwd=2,col=3,xlab="",ylab="",main="")
legend("topleft",legend = c("Branching","Linear","Punctuated"),lwd=2,col=c(1,2,3),bty="n")
dev.off()

######trainning the distribution parameter for each selective pattern corresponding to Table 1
#Branching
logBranch <- function(theta){
  mu = theta[1]
  sigma = theta[2]
  logBranch <- -sum(dnorm(branchcds,mu,sigma,log = T))
  return(logBranch)
}
optB <- optim(c(mean(branchcds),sd(branchcds)),logBranch)
#Linear
logLinear <- function(theta){
  mu = theta[1]
  sigma = theta[2]
  logLinear <- -sum(dnorm(linearcds,mu,sigma,log = T))
  return(logLinear)
}
optL <- optim(c(mean(linearcds),sd(linearcds)),logLinear)
#Punctuated
logPunctuated <- function(theta){
  mu = theta[1]
  sigma = theta[2]
  logPunctuated <- -sum(dnorm(punctuatedcds,mu,sigma,log = T))
  return(logPunctuated)
}
optP <- optim(c(mean(punctuatedcds),sd(punctuatedcds)),logPunctuated)


sigma <- c(optB$par[2],optL$par[2],optP$par[2])
mu <- c(optB$par[1],optL$par[1],optP$par[1])

approachPar <- data.frame(model = c("Branching","Linear","Punctuated"),mu = round(mu,3),sigma = round(sigma,3))
write.table(approachPar,"./simulation/trainningset/selectiveModel.predict.para.txt",sep="\t",col.names = T,row.names = F,quote = F)

###Figure 2F: accuracy of trainning set
res1 = res[res$model!="Neutral",]
selecEst <- ModelPredict(res1$cds)
res1 <- cbind(res1,selecEst)
res1$branchLabel <- "F"
res1$branchLabel[res1$model=="Branching"]="T"
res1$linearLabel <- "F"
res1$linearLabel[res1$model=="Linear"] <- "T"
res1$punctuatedLabel <- "F"
res1$punctuatedLabel[res1$model=="Punctuated"] <- "T"

pred <- prediction(res1$branch, res1$branchLabel)
perfB <- performance(pred,"tpr","fpr")
aucB  <- performance(pred, 'auc')
aucB  <- unlist(slot(aucB,"y.values"));

pred <- prediction(res1$linear, res1$linearLabel)
perfL <- performance(pred,"tpr","fpr")
aucL  <- performance(pred, 'auc')
aucL  <- unlist(slot(aucL,"y.values"));

pred <- prediction(res1$punctuated, res1$punctuatedLabel)
perfP <- performance(pred,"tpr","fpr")
aucP  <- performance(pred, 'auc')
aucP  <- unlist(slot(aucP,"y.values"));
pdf(file = "./simulation/trainningset/selectiveModel.ROC.pdf",width=6,height = 6)
plot(perfB,lwd=2,col=1)  
plot(perfL,lwd=2,add=T,col=2)  
plot(perfP,lwd=2,add=T,col=3) 
legend("bottomright",legend = paste(c("Branching AUC = ","Linear AUC = ","Punctuated AUC = "),c(round(aucB,2),round(aucL,2),round(aucP,2)),sep=""),col=c(1:3),lwd=2,,bty="n")
dev.off()

######Noise effect on selective model
###Figure 2H
library(ggplot2)
res <- c()
for (m in model){
  CNVdata = read.csv(paste("./simulation/noise/",m,".CNV.txt",sep=""),sep="\t")
  for (n in noise){
    subdata=CNVdata[CNVdata$noise==n,]
    dataCount <- unique(subdata$Nth_simu)
    cds <- sapply(dataCount, function(i,subdata){
      cds <- selectModel(subdata$freq[subdata$Nth_simu==i])
      return(cds)
    },subdata)
    res <- rbind(res,data.frame(cds = cds,noise = n,model = m)) 
  }
}
res1 <- res[res$model!="Neutral",] 
res1$noise <- as.factor(res1$noise)
p <- ggplot(res1, aes(x=model, y=cds,color = noise)) + geom_boxplot()
p <- p + geom_hline(yintercept=0.88, linetype="dashed", color = 1)
p <- p + geom_hline(yintercept=0.428+0.127, linetype="dashed", color = 1)
pdf(file="./simulation/noise/CDS.nosie.pdf",width = 5,height = 6)
p
dev.off()

###Figure 2I
selecEst <- ModelPredict(res1$cds)
res1 <- cbind(res1,selecEst)
res1$branchLabel <- "F"
res1$branchLabel[res1$model=="Branching"]="T"
res1$linearLabel <- "F"
res1$linearLabel[res1$model=="Linear"] <- "T"
res1$punctuatedLabel <- "F"
res1$punctuatedLabel[res1$model=="Punctuated"] <- "T"

AUC <- c()
for (n in noise){
  subres=res1[res1$noise == n,]
  pred <- prediction(subres$branch, subres$branchLabel)
  perfB <- performance(pred,"tpr","fpr")
  aucB  <- performance(pred, 'auc')
  aucB  <- unlist(slot(aucB,"y.values"));
  pred <- prediction(subres$linear, subres$linearLabel)
  perfL <- performance(pred,"tpr","fpr")
  aucL  <- performance(pred, 'auc')
  aucL  <- unlist(slot(aucL,"y.values"));
  pred <- prediction(subres$punctuated, subres$punctuatedLabel)
  perfP <- performance(pred,"tpr","fpr")
  aucP  <- performance(pred, 'auc')
  aucP  <- unlist(slot(aucP,"y.values"));
  AUC <- rbind(AUC,data.frame(noise = n, AUCB = aucB, AUCL = aucL, AUCP = aucP))
}
pdf(file = "./simulation/noise/selectiveModel.nosise.AUC.pdf",width = 5,height = 5)
plot(AUC$AUCB,type="o",pch=20,cex=2,ylim=c(0,1),lwd=2,col=1,xlim=c(0,11),xlab="Noise",ylab="AUC")
par(new=T)
plot(AUC$AUCL,type="o",pch=20,cex=2,ylim=c(0,1),lwd=2,col=2,xlim=c(0,11),xlab="",ylab="")
par(new=T)
plot(AUC$AUCP,type="o",pch=20,cex=2,ylim=c(0,1),lwd=2,col=3,xlim=c(0,11),xlab="",ylab="")
legend("bottomleft",legend = c("Branching","Linear","Punctuated"),lwd=2,col=c(1,2,3),bty="n")
dev.off()






