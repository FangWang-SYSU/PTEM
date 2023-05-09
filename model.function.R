library(diptest)


selectModel <- function(CNVfre,cutoff1 = 0.05,cutoff2 = 1){
  subfre <- CNVfre[CNVfre>cutoff1&CNVfre<=cutoff2]
  fre <- round(subfre,1)
  fre_tab <- table(fre)
  refer <- data.frame(group = seq(0,1,by= 0.1),fre = rep(0,length=11))
  refer$fre[!is.na(match(refer$group,names(fre_tab)))]=fre_tab
  fre_tab <- fre_tab/sum(fre_tab)
  CNVdiversity = -sum(as.numeric(names(fre_tab))*fre_tab*log(fre_tab))
  return(CNVdiversity)
}

ModelPredict <- function(CNVdiversity){
  model <- c("Branching","Linear","Punctuated")
  proB <- dnorm(CNVdiversity,mean=0.688,sd=0.082)
  proL <- dnorm(CNVdiversity,mean=0.930,sd=0.078)
  proP <- dnorm(CNVdiversity,mean=0.428,sd=0.127)
  pro <- data.frame(branch = proB,linear = proL, punctuated = proP)
  prosum <- apply(pro,1,sum)
  postPro <- apply(pro, 2, function(x,prosum){return(x/prosum)},prosum=prosum)
  if(length(CNVdiversity)>1){
    postPro <- as.data.frame(postPro) 
  }else{
    postPro <- data.frame(branch=postPro[1],linear=postPro[2],punctuated=postPro[3]) 
    row.names(postPro) <- 1
  }
  predictM <- apply(postPro,1,function(x,model){
    k = which.max(x)
    return(model[k])
  },model)
  postPro$Predict <- predictM
  return(postPro)
}
