###############################################################################
# Module 3: Training Set Optimization in R. 
# Authors: Waseem Hussain and  Mahender Anumalla 
# Date: August 6, 2024
###############################################################################

# Load the Required Libraries
  rm(list=ls()) # Remove previous work
  library(AGHmatrix)
  library(STPGA) # Package for Training Optimization
  library(FactoMineR)
  library(factoextra)


## ----warning = FALSE, message = FALSE---------------------------------------
# Genotype data
  geno<-readRDS("./Data/geno.oyt.filtered.rds") # Upload marker data
  dim(geno) # Check dimensions
  geno[is.na(geno)] <-0
# Recode the matrix in 2,1 and 0 numeric format
  geno[geno==1]<-2
  geno[geno==-1]<-0
  geno[geno==0]<-0
  kable(head(geno)[1:6,1:2])


## ----warning = FALSE, message = FALSE---------------------------------------
# Get subset of first 200 genotypes
  geno<-geno[1:200, ] # Subset 
  geno<-as.matrix(geno) # Convert as matrix
# Build the VanRaden 2008 G matrix
  G_additive <- Gmatrix(SNPmatrix=geno, missingValue=NA, 
                          maf=0.05, method="VanRaden")


## ----message=FALSE, warning=FALSE-------------------------------------------
# Define the number of individuals to select for the training set
  nSel <- 60  # Can be adjusted based on population size

# Optimize the training set using CDmean in the argument errorstat

    result <- GenAlgForSubsetSelectionNoTest(
      P =  G_additive,  # marker matrix or relationship matrix
      ntoselect = nSel, # Number of individuals in the training set
      npop = 100,  # population size for the genetic algorithm
      nelite = 5, # population size for the genetic algorithm
      mutprob = 0.5, # probability of mutation for each generated solution
      mutintensity = 1,
      niterations = 50, # number of iterations.
      minitbefstop = 50, # number of iterations before stopping
      tabu = FALSE,
      plotiters = FALSE, # plot the convergence
      lambda = 1e-6, # scalar shrinkage parameter
      errorstat = "CDMEAN"
    )


## ----message=FALSE, warning=FALSE-------------------------------------------
  selected_ids <- data.frame(TrainingGenotypes=result$`Solution with rank 1`)
  kable(head(selected_ids$TrainingGenotypes))


## ----message=FALSE, warning=FALSE-------------------------------------------
#  Get lines from training set
  line.names<-data.frame(Genotype=ifelse(row.names(G_additive)%in%selected_ids$TrainingGenotypes, "Selected", "Un-selected"))
  row.names( line.names)<-row.names(G_additive)

# Perfom PCA on GM matrix
  pca<- PCA(G_additive, graph = FALSE)
# Get the Biplot
  biplot.gm<-fviz_pca_biplot(pca, palette = "jco", axes = c(1, 2), geom="point",pointsize = 3,
                           addEllipses = FALSE,col.ind = line.names$Genotype,
                           geom.var= "none", label = "none", legend.title = "Group")+
  theme_minimal()+
  # add and modify the title to plot
  theme (
    plot.title = element_blank(),
    # add and modify title to x axis
    axis.title.x = element_text(color="black", size=12), 
    # add and modify title to y axis
    axis.title.y = element_text(color="black", size=12)) +
  # modify the axis text
  theme(axis.text= element_text(color = "black", size = 12)) 
biplot.gm
#' ###################End########################################


