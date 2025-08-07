## ----knitr_init, echo=FALSE, cache=FALSE---------------------------------------------------
library(knitr)
library(kableExtra)
## Global options
options(max.print="100")
opts_chunk$set(echo=TRUE,
               cache=TRUE,
               prompt=TRUE,
               collapse=FALSE,
               comment=NA,
               strip.white=TRUE,
               message=FALSE,
               warning=FALSE,
               width=65,
               tidy.opts=list(width.cutoff=65, tidy=TRUE))


## ----setup, include=FALSE, echo=FALSE------------------------------------------------------
  require("knitr")
  opts_knit$set(root.dir = "~/Documents/Research/Workshops/IRRI_2025/Module5")


## ----warning=FALSE, message=FALSE----------------------------------------------------------
  library(AGHmatrix)
  library(BGLR)
  library(lme4)
  library(ggplot2)
  library(sommer)


## ------------------------------------------------------------------------------------------
  rm(list=ls()) # remove History
# Read the phenotypic data
  BLUEs.all<-read.csv(file="./Data/BLUES.ALL.csv")
  BLUEs.all<-subset(BLUEs.all, Environment=="ENV1")


## ------------------------------------------------------------------------------------------
  geno<-readRDS("./Data/GBS_datav2.rds")
  dim(geno)
# Match genotype with Phenotype
  Ids<-unique(BLUEs.all$Genotype)
  length(Ids)
# Now subet the genotype Data based on IDs
  geno<-geno[row.names(geno)%in%Ids,]
  dim(geno)


## ------------------------------------------------------------------------------------------
  GM<- Gmatrix(SNPmatrix=geno, missingValue=NA, 
                            maf=0.05, method="VanRaden")
  dim(GM)


## ------------------------------------------------------------------------------------------
# Let us take position of markers based on columns
   qtl <- c(1120, 1126, 1128, 1129)

# Create a ETA list
    ETA <- list(
      list(X = geno[, qtl], model = 'FIXED', probIn = 0.10),
      list(K = GM, X = geno[, -qtl], model = 'RKHS', probIn = 0.10)
    )

# Fit the model
  model_fix<- BGLR(y=BLUEs.all$BLUEs, ETA=ETA, burnIn=500, nIter=2000, 
                        verbose=FALSE)
# Extract the GBVs
  GEBVs_fixed<-data.frame(GEBVs=model_fix$yHat)

