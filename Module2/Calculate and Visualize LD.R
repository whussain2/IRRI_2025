###############################################################################
# Module 2: Linkage Disequillibrum in R. 
# Authors: Waseem Hussain and  Mahender Anumalla 

###############################################################################

# Set working Directory
  setwd("~/Documents/Research/Workshops/IRRI-IIRR_2025/Module2")

  rm(list=ls()) # remove the previous history
# Install
  #install.packages("BGLR")
  #install.packages("genetics")
# Installing snpStats package from Bioconuctor
  #if (!requireNamespace("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")
#BiocManager::install("snpStats")
# Load the packages
  library(BGLR)
  library(genetics)
  library(pheatmap)
  #library(LDheatmap)

# Read the mice marker data
  #mice.X<-read.csv(file="mice.X.csv", header = TRUE)
# Load the mice data
  data(mice)
# Subset the mice data, first 20 markers
  mice.20<- mice.X[, 1:20] # use the first 10 markers
# Visualize first 5 rows and columns
  mice.20[1:5, 1:5] # Data is coded 0, 1 and 2
# Make genotypes
  mice.20.G<- makeGenotypes(mice.20, convert=c(colnames(mice.20)), 
                            method=as.genotype.allele.count)
# Visualize first 5 rows and columns
  mice.20.G[1:5, 1:5] # Data is coded 0, 1 and 2
# Now calculate the LD
  LD.20<- LD(mice.20.G) # This will return the list
  names (LD.20)
# Extract r2 ( Hill and Robertson (1968)
  r2<-LD.20$`R^2`
# Copy upper part of matrix to lower for visualizations
  lowerTriangle(r2) <- upperTriangle(r2)
# Convert Diagonal to 1
  diag(r2)<-1

ld.map<-heatmap(r2)
ld.map
#########################END##################################################