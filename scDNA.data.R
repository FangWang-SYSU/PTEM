####application on real data
workpath="./PTEM"
setwd(workpath)
source("model.function.R")

###scDNA-seq
library(dplyr)
#####evolution pattern prediction on TNBC scDNA-seq data
###Figure 3A
TNBCdata <- read.csv("./scDNA/TNBC.CNV.fre.txt",sep="\t")
sampleID <- unique(TNBCdata$sampleID)
modelres <- c()
for (ID in sampleID){
  CNVevent <- TNBCdata[TNBCdata$sampleID==ID,]
  modalityEst <- dip.test(CNVevent$fre[CNVevent$fre>0.05])
  modelEst <- data.frame(ID = ID,modality.P = modalityEst$p.value)
  cds <- selectModel(CNVevent$fre,cutoff1 = 0.01,cutoff2 = 0.76)
  modelEst <- cbind(modelEst,ModelPredict(cds))
  modelEst$Predict[modelEst$modality.P>0.05]="Neutral"
  modelres <- rbind(modelres,modelEst)
}
modelres

mdel_tab <- table(modelres$Predict)
df=data.frame(model=names(mdel_tab),value=as.numeric(mdel_tab))
bp<- ggplot(df, aes(x="", y=value, fill=model))+
  geom_bar(width = 1, stat = "identity")
pie <- bp + coord_polar("y", start=0)
pdf(file = "./scDNA/TNBC.scDNA.modelpie.pdf",width = 5,height = 5)
pie + scale_fill_brewer(palette="Dark2")
dev.off()

####Figure 3B
TNBCdata <- read.csv("./scDNA/TNBC.posttreatment.CNV.txt",sep="\t")
sampleID <- unique(TNBCdata$sampleID)
treatres <- c()
for (ID in sampleID){
  CNVevent <- TNBCdata[TNBCdata$sampleID==ID,]
  modalityEst <- dip.test(CNVevent$fre[CNVevent$fre>0.05])
  modelEst <- data.frame(ID = ID,modality.P = modalityEst$p.value)
  cds <- selectModel(CNVevent$fre,cutoff1 = 0.01,cutoff2 = 0.76)
  modelEst <- cbind(modelEst,ModelPredict(cds))
  modelEst$Predict[modelEst$modality.P>0.05]="Neutral"
  modelEst$treat = "Post"
  treatres <- rbind(treatres,modelEst)
}

treatModel <- treatres
treatModel$postModel <- treatModel$Predict
treatModel$preMode <- modelres$Predict[match(treatModel$ID,modelres$ID)]
NR <- treatModel
NRdata <- data.frame(source = c("Punctuated","Punctuated"),target = c("Branching","Punctuated"),value = c(1/4,3/4))
NRdata <- NRdata %>% mutate_if(is.character, as.factor) # convert to factor
NRdata$target <- paste(NRdata$target, " ", sep="")
nodes <- data.frame(name=c(as.character(NRdata$source), as.character(NRdata$target)) %>% unique())#制作nodes
NRdata$IDsource=match(NRdata$source, nodes$name)-1 
NRdata$IDtarget=match(NRdata$target, nodes$name)-1
ColourScal='d3.scaleOrdinal() .range(["#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1" ,"#6BAED6", "#4292C6" ,"#2171B5" ,"#084594"])'

library(webshot)
library(networkD3)

sankeyplot <- sankeyNetwork(Links = NRdata, Nodes = nodes,
                            Source = "IDsource", Target = "IDtarget",
                            Value = "value", NodeID = "name", 
                            sinksRight=FALSE, colourScale=ColourScal, nodeWidth=50, fontSize=16, nodePadding=10)
saveNetwork(sankeyplot, file = "./scDNA/TNBC.treat.sankeyplot.html") #save to html
webshot::webshot("./scDNA/TNBC.treat.sankeyplot.html", file =  "./scDNA/TNBC.treat.sankeyplot.pdf", vwidth = 500, vheight = 500)

#####Figure 3C
library(ggplot2)
library(ggpubr)
modelres <- read.csv("./scDNA/TNBC.model.CNVscore.txt",sep="\t")
pdf(file = "./scDNA/CNVscore.boxplot.pdf",width = 3,height = 4)
p <- ggboxplot(modelres, x = "Predict", y = "CNVscore",
               color = "Predict", palette = "jco",
               add = "jitter")
p+ stat_compare_means(method = "t.test")
dev.off()






