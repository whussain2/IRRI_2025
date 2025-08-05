###############################################################################
# Module 3: Trainingset Optimization in R. 
# Authors: Waseem Hussain and  Mahender Anumalla 
# Date: August 6, 2024
###############################################################################

# # Load the Required Libraries
#   rm(list=ls()) # Remove previous work
#   library(dplyr)
#   library(readxl)
#   library(TrainSel)
#   library(factoextra)
#   library(FactoMineR)
#   library(rrBLUP)


## ----message=FALSE, echo=FALSE, eval=FALSE----------------------------------
#   geno <- read.delim("./data/oyt.drt_num.txt", row.names=1)
#   dim(geno)
#   geno <- as.matrix(geno)
# # Impute missing Data
#   imputation <- A.mat(geno, impute.method = "EM", return.imputed = T, min.MAF = 0.05)
#   geno <- imputation$imputed


## ----message=FALSE, echo=TRUE, eval=FALSE-----------------------------------
# ###########################################
# # Generating a training population
# ###########################################
# 
# #Objective: Training set optimization consists of choosing a set of training individuals that will better predict un-phenotyped germplasm in a TS
# 
#   library(TrainSel)
#   control=TrainSelControl()
# 
# # When using a mixed model based CDmin statistic we need to prepare the data using the ’MakeTrainSel-Data function
# # We will use marker data to prepare data for mixed model based criterion.
#   TSData<-MakeTrainSelData(M=geno)
# ## We use ’TrainSel’ function to select an Unordered Sample ("UOS") of size 20 from the first 100 individuals in the dataset
# #as a default it picks, corresponding to the first 100 rows of data. The default settings do not need to be modified for small to medium-sized optimization problems. We select settypes as "OS"power if to perform multi-objective optimization, else we will use "UOS".
#   controlOuter<-TrainSelControl()
#   controlOuter$niterations<-5
#   T.set<-TrainSel(Data=TSData, Candidates=list(1:100), setsizes=c(20),settypes="UOS", control=
#                  control)
# 
# # Subset the best complete genotypic data of the best picked 10 genotypes from the genotype file.
#   T.set<- T.set$BestSol_int
# # Subset the best complete genotypic data of the best picked 100 genotypes from the genotype file
#   TS.final <- geno[T.set,]


## ----echo=TRUE, include=TRUE, warning=FALSE, message=FALSE, eval=FALSE------
# # First Get genomic matrix
# #Xs <- scale(geno.oyt, center = TRUE)
# # Construct G matrix
#   GM <- geno %*% t(geno)/ncol(geno)
#   dim(GM)


## ----echo=TRUE, include=TRUE, warning=FALSE, message=FALSE, eval=FALSE------
# #  Get lines from training set
#   line.names<-data.frame(Genotype=ifelse(row.names(geno)%in%row.names(TS.final), "Un-selected", "Selected"))
# # Perfom PCA on GM matrix
#   pca<- PCA(GM, graph = FALSE)
# # Get the Biplot
#   biplot.gm<-fviz_pca_biplot(pca, palette = "jco", axes = c(1, 2), geom="point",pointsize = 3,
#                 addEllipses = FALSE,col.ind = line.names$Genotype,
#                 geom.var= "none", label = "none", legend.title = "Group")+
#   theme_minimal()+
#   # add and modify the title to plot
#   theme (
#     plot.title = element_blank(),
#     # add and modify title to x axis
#     axis.title.x = element_text(color="black", size=12),
#     # add and modify title to y axis
#     axis.title.y = element_text(color="black", size=12)) +
#   # modify the axis text
#   theme(axis.text= element_text(color = "black", size = 12))
#   biplot.gm

