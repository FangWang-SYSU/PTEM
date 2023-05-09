# PTEM
A two-step approach to Predict Tumor Evolution Model

## Description
This is an apporch to predict tumor evolution model using single cell copy number profile. It will estiamted:
* It the tumo cell population follows neutral or selective evolution model.
* If the tumor cell popoulation follows selective eovlution model, it is linear, branching or punctuated eovlution model.

## System requirements and dependency
This approach runs on R (version > 4.0) and has dependency on the R packages: diptest.

## Usage
Please download and copy the distribution to your specific location.
```
  setwd("./PTEM")
  source("model.function.R")
```
Input a CNV frequncey vector calculated from a single cell coy number profile. Each element in the vector corresponds to the frequency of genome region/gene gain or loss in the tumo cell population.
  >*Generate a CNV frequency vector from neutral evolution model i.e. unimodal distribution* 
>```
>  CNVfre1 <- rnorm(500,mean=0.3,sd=0.02)
>```
  >*Generate a CNV frequency vector from selective evolution model i.e. multimodal distribution*
>```
>  CNVfre2 <- c(rnorm(100,mean=0.1,sd=0.02),rnorm(200,mean=0.4,sd=0.02),rnorm(200,mean=0.6,sd=0.02))
>```
Estimating if the CNV frequency vector follows neutral evolution
```
  modalityEst1 <- dip.test(CNVfre1)
  modalityEst2 <- dip.test(CNVfre2)
```
>*For scDNA-seq, keep CNV frequency > 0.05; For scRNA-seq, keep CNV frequenct > 0.02*

**Neutral evolution model**
>```
>modalityEst1
>
>Hartigans' dip test for unimodality / multimodality
>
>data:  CNVfre1
>D = 0.013069, p-value = 0.8705
>alternative hypothesis: non-unimodal, i.e., at least bimodal
>```

**Selective evolution model**
>```
>modalityEst2
>
>Hartigans' dip test for unimodality / multimodality
>
>data:  CNVfre2
>D = 0.11944, p-value < 2.2e-16
>alternative hypothesis: non-unimodal, i.e., at least bimodal
>```

For selectiv evolution model, further distinguish linea, branching and punctuated evolution model
```
cds <- selectModel(CNVfre2,cutoff1=0.01,cutoff2=0.7)
modelEst <- ModelPredict(cds)
```
>```
>modelEst
>
>   branch       linear punctuated    Predict
>1 0.01829491 5.971049e-09  0.9817051 Punctuated
>```
>*The first to the third column are the posterior probabilities of corresponding selective eovlution model. The Predict column returns the model with the largest posterior probability.*

  
