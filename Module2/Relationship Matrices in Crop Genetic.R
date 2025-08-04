###############################################################################
# Course: Fundamentals of Genomic Prediction-2025
# Module 1: Relationship Matrices in R. 
# Authors: Waseem Hussain and  Mahender Anumalla 

###############################################################################


  #install.packages("AGHmatrix")
  library ("AGHmatrix")
  library(DT)


## -----------------------------------------------------------------------------------
# Read the Pedigree Data
  ped<-read.csv(file="./Data/ped.finalv2.csv", header = TRUE)
  head(ped)
# Use Amatrix function to build pedigree matirx
  ped.matrix<-Amatrix(ped, ploidy=2)
  dim(ped.matrix)
  ped.matrix[1:10,1:5]


## ----message=FALSE, warning=FALSE---------------------------------------------------
# Read SNP data
  geno<-read.csv(file="./Data/geno.data.csv", header = TRUE, row.names = 1)
  geno<-as.matrix(geno) # Convert as matrix
# Build the VanRaden 2008 G matrix
  G_additive <- Gmatrix(SNPmatrix=geno, missingValue=NA, 
                          maf=0.05, method="VanRaden")


## -----------------------------------------------------------------------------------
  heatmap(G_additive)


## ----message=FALSE, warning=FALSE---------------------------------------------------
  G_Dominace <- Gmatrix(SNPmatrix=geno, missingValue=NA, 
                    maf=0.05, method="Su")


## -----------------------------------------------------------------------------------
  heatmap(G_Dominace)

