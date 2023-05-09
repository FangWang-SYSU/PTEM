library(diptest)
library(LaplacesDemon)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
workpath="./PTEM"
setwd(workpath)
source("model.function.R")

##input CNV frequence inferred by inferCNV
RNAdata <- read.csv("./scRNA/inferCNV.txt",sep="\t")
cancerID <- unique(RNAdata$cancer)
modelres <- c()
for (cancer in cancerID){
  CNVevent <- RNAdata[RNAdata$cancer==cancer,]
  sampleID <- unique(CNVevent$sampleID)
  if (cancer == "Breastcancer"|cancer == "CRC"){
    threshold = 0.02
  }else{
    threshold = 0.05
  }
  for (ID in sampleID){
    subdata=CNVevent[CNVevent$sampleID==ID,]
    modalityEst <- dip.test(subdata$fre[subdata$fre>threshold])
    modelEst <- data.frame(ID = ID,cancer=cancer,modality.P = modalityEst$p.value)
    cds <- selectModel(subdata$fre,cutoff1 = 0.01,cutoff2 = 1)
    modelEst <- cbind(modelEst,ModelPredict(cds))
    modelEst$Predict[modelEst$modality.P>0.05]="Neutral"
    modelres <- rbind(modelres,modelEst)
  }
}
modelres

write.table(modelres,"./scRNA/scRNA.data.eolutionmodel.txt",sep = "\t",col.names = T,row.names = F,quote = F)


####Figure 4E
library(RColorBrewer)
RNAmodel <- read.csv("./scRNA/scRNA.data.eolutionmodel.txt",sep="\t")
modelFrac <- table(RNAmodel$Predict[RNAmodel$cancer=="ccRCC"])
modelFrac <- rbind(modelFrac,table(RNAmodel$Predict[RNAmodel$cancer=="CRC"]))
modelFrac <- rbind(modelFrac,table(RNAmodel$Predict[RNAmodel$cancer=="Lungcancer"]))
modelFrac <- rbind(modelFrac,table(RNAmodel$Predict[RNAmodel$cancer=="Breastcancer"]))
row.names(modelFrac)= c("ccRCC","Colorectal cancer","Lung cancer","Breast cancer")
modelFrac1 <- apply(modelFrac, 1, function(x){return(x/sum(x))})
mycol=brewer.pal(3, "Dark2")
pdf(file = "./scRNA/RNAmodel.frac.pdf",width = 6,height = 6)
par(mar=c(4,8,4,1))
barplot(modelFrac1,horiz = T,las=2,col=mycol,xlab="Fraction")
legend("top",legend = c("Branching","Neutral","Punctuated"),col=mycol,pch = 15,cex=2)
dev.off()

#####The hallmark and metaprogram of tumor
#####Figure 4F and G 
##meta program
MPdata <- read.table("./scRNA/tumor.MP.score.txt",header = T,sep="\t")
sampleAnno <- MPdata[,1:2]
scaleData <- MPdata[,3:dim(MPdata)[2]]
col_gain = colorRamp2(c(-2, 0,2), c("blue", "white","red"))
sampleOrder <- row.names(sampleAnno)[order(sampleAnno$model)]
modelcolor=mycol
cancercolor=brewer.pal(4, "Set3")
mat_colors <- list(model=modelcolor,cancer=cancercolor)
names(mat_colors$model)=c("Branching","Neutral","Punctuated")
names(mat_colors$cancer)=c("Breast cancer","CRC","ccRCC","Lung cancer")
pdf(file="./scRNA/MP.heatmap.pdf",width = 7,height = 6)
pheatmap(t(scaleData[sampleOrder,]),color = col_gain,show_colnames = F,annotation_col = sampleAnno,cluster_cols = F,annotation_colors = mat_colors,row_km=3)
dev.off()


MPdata <- read.table("./scRNA/tumor.Hallmark.score.txt",header = T,sep="\t")
sampleAnno <- MPdata[,1:2]
scaleData <- MPdata[,3:dim(MPdata)[2]]
col_gain = colorRamp2(c(-2, 0,2), c("blue", "white","red"))
sampleOrder <- row.names(sampleAnno)[order(sampleAnno$model)]
modelcolor=mycol
cancercolor=brewer.pal(4, "Set3")
mat_colors <- list(model=modelcolor,cancer=cancercolor)
names(mat_colors$model)=c("Branching","Neutral","Punctuated")
names(mat_colors$cancer)=c("Breast cancer","CRC","ccRCC","Lung cancer")

pdf(file="./scRNA/HALLMARK.heatmap.pdf",width = 7,height = 6)
pheatmap(t(scaleData[sampleOrder,]),color = col_gain,show_colnames = F,annotation_col = sampleAnno,cluster_cols = F,annotation_colors = mat_colors,row_km=3)
dev.off()


###The tumor microenvironment composition
###Figure 4H to K
cellpro <- read.csv("./scRNA/Cell_Composition_total.txt",sep="\t")
cellpro$celltype[cellpro$celltype=="Nephron"]="Epi"

celltypenames <- unique(cellpro$celltype)
pdf(file = "./scRNA/tumorcell.proportion.pdf",width = 4,height = 4)
cellpro1 <- cellpro[cellpro$celltype=="Epi",]
p <- ggboxplot(cellpro1, x = "Predict", y = "prop_in_total",
               color = "Predict", palette = "jco",
               add = "jitter")
p+ stat_compare_means()
dev.off()
pdf(file = "./scRNA/Tcell.proportion.pdf",width = 4,height = 4)
cellpro1 <- cellpro[cellpro$celltype=="T",]
p <- ggboxplot(cellpro1, x = "Predict", y = "prop_in_total",
               color = "Predict", palette = "jco",
               add = "jitter")
p+ stat_compare_means()
dev.off()


cellpro$Predict <- factor(cellpro$Predict,levels = c("Punctuated","Branching","Neutral"))
cellpro1 <- cellpro[cellpro$celltype=="Epi",]
pdf(file="./scRNA/tumor.pro.cancer.pdf",width=4.5,height = 3)
p <- ggboxplot(cellpro1, x = "CancerTpye", y = "prop_in_total",
               color = "Predict", palette = "jco",
               add = "jitter")
p
dev.off()
pdf(file="./scRNA/T.pro.cancer.pdf",width=4.5,height = 3)
cellpro1 <- cellpro[cellpro$celltype=="T",]
p <- ggboxplot(cellpro1, x = "CancerTpye", y = "prop_in_total",
               color = "Predict", palette = "jco",
               add = "jitter")
p
dev.off()









