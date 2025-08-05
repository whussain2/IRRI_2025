###############################################################################
# Module 3: Genomic Prediction Models in R. 
# Authors: Waseem Hussain and  Mahender Anumalla 
# Date: August 6, 2024
###############################################################################


# Load the Required Libraries
  rm(list=ls()) # Remove previous work
  library(rrBLUP)
  library(BGLR)
  library(AGHmatrix)
  library(ggplot2)
  library(DT)
  library(cvTools)


## ----warning = FALSE, message = FALSE---------------------------------------
# Read the phenotypic data set
  pheno<-readRDS("./Data/pheno.filtered.rds")
  pheno<-pheno[,1:5] # Removing missing columns
# View data as table
  head(pheno)


## ----warning = FALSE, message = FALSE---------------------------------------
# Genotype data
  geno<-readRDS("./Data/geno.oyt.filtered.rds")
  dim(geno) # Check dimensions
  kable(head(geno)[1:6,1:2])
  geno[is.na(geno)] <-0


## ---------------------------------------------------------------------------
# Mixed model in rrBLUP
  model.rrblup<- mixed.solve(y =pheno[,2], Z = geno) # RR-BLUP


## ---------------------------------------------------------------------------
# Get Marker effects
  kable((model.rrblup$u)[1:5])
# Get GEBVs (Multiple marker effects by Matrix)
  GEBVs = data.frame(EBVs=geno %*% model.rrblup$u)
  kable(head(GEBVs))


## ----echo=TRUE, eval=TRUE---------------------------------------------------
# Set the ETA
###----------- Linear predictor 
ETA = list(A = list(X = geno,   # X (Markers) in ETA
                    model = 'BRR'))  # Specify the model
                                 
# Phenotypic data
Y = pheno[,2]

# Run Model for trait 1
model.rrblupB = BGLR(y=Y,         # Phenotypic data
               ETA=ETA,          # Model ETA
               nIter=1200,       # Number of iterations for the model
               burnIn=120,       # Burnin iterations
               thin=5,           # Sampling throughout iterations
               verbose = FALSE)


## ----echo=TRUE,eval=TRUE----------------------------------------------------
# Mean estimated
  model.rrblupB$mu
#Variance component for SNP  marker effects
  model.rrblupB$ETA[[1]]$varB
#Variance component for the residual effect
  model.rrblupB$varE
# SNP effects
  length(model.rrblupB$ETA[[1]]$b) # markers effects
#Get GEBVs
  GEBVs.B =data.frame(EBVs=geno%*% model.rrblupB$ETA[[1]]$b)
   kable(head(GEBVs.B))


## ----warning=FALSE, message=FALSE-------------------------------------------
# Marker data as matrix
  geno1<-as.matrix(geno)
  geno1[geno1 ==1] <- 2 # Conversion
  geno1[geno1 ==-1] <- 1 # Conversion
  geno1[geno1 ==0] <- NA # Conversion
# Build the G matrix
  G_matrix = Gmatrix(SNPmatrix=geno1, method = "VanRaden", 
                   ploidy = 2,integer = FALSE, missingValue =NA)


## ----echo=TRUE, eval=TRUE, warning = FALSE, message = FALSE-----------------
# Run the GBLUP model
  model.gblup<- mixed.solve(y = pheno[,2], Z = G_matrix) # G-BLUP
  gebvs.gblup<-data.frame(Gebvs=model.gblup$u)
  # Add intercept
  intercept<-model.gblup$beta
  gebvs.gblup$Gebvs<-gebvs.gblup$Gebvs+intercept
  kable(head(gebvs.gblup))


## ----echo=TRUE, eval=TRUE---------------------------------------------------
# Pheno
  Y = as.matrix(pheno[,2])
# Model
  model.BayesA<- BGLR(y=Y,# Phenotypic data blues non-stress
             ETA=list(list(X=geno, model='BayesA')),          # Model Bayes A
             nIter=1200,       # Number of iterations for the model
             burnIn=120,       # Burnin iterations
             thin=5,
             verbose = FALSE)           # Sampling throughout iterations
# Mean estimated
  model.BayesA$mu
# Variance component for the residual effect
  model.BayesA$varE
# SNP effects
  length(model.BayesA$ETA[[1]]$b) #'b' slot
  
# Assignment: How to get Estimated Breeding Values


## ----echo=TRUE, eval=TRUE---------------------------------------------------
# Pheno
  Y = as.matrix(pheno[,2])
#Model
  model.BayesB<- BGLR(y=Y,# Phenotypic data blues non-stress
             ETA=list(list(X=geno, model='BayesB')),          # Model Bayes A
             nIter=1200,       # Number of iterations for the model
             burnIn=120,       # Burnin iterations
             thin=5,
             verbose = FALSE)           # Sampling throughout iterations

# Mean estimated
  model.BayesB$mu
# Variance component for the residual effect
  model.BayesB$varE
# SNP effects
  length(model.BayesB$ETA[[1]]$b) #'b' slot


## ----echo=TRUE, eval=TRUE---------------------------------------------------
# Pheno
  Y = as.matrix(pheno[,2])
#Model
  model.BayesC<- BGLR(y=Y,# Phenotypic data blues non-stress
             ETA=list(list(X=geno, model='BayesC')),          # Model Bayes A
             nIter=1200,       # Number of iterations for the model
             burnIn=120,       # Burnin iterations
             thin=5,
             verbose = FALSE)           # Sampling throughout iterations

# Mean estimated
  model.BayesC$mu
# Variance component for the residual effect
  model.BayesC$varE
# SNP effects
  length(model.BayesC$ETA[[1]]$b) #'b' slot


## ---------------------------------------------------------------------------

# Plotting the effects
  BayesB_rrBLUP = data.frame("Marker_BayesB" = model.BayesB$ETA[[1]]$b,
                           "Marker_rrBLUP" = model.rrblupB$ETA[[1]]$b)

  ggplot(BayesB_rrBLUP, aes(x = Marker_BayesB, y = Marker_rrBLUP))+
  geom_point()+
  theme_classic()



## ---------------------------------------------------------------------------
# Get the number of folds and replications
  nfolds = 4 # Four folds
  nrep = 2
#  Order Pheno and geno files
  pheno = pheno[order(pheno$Designation, decreasing = FALSE),] 
  geno= geno[order(row.names(geno)),]
  all(rownames(geno) == pheno$Designation)
  pheno$Designation<-as.factor(pheno$Designation)
# Create a factor OF idS
  pheno$ID = factor(1:nrow(pheno))
  set.seed(100)
# Create fold assignment using forloop
  sort<- list()
  for(a in 1:nrep){
    for(j in 1:nfolds){
      folds <- cvFolds(nlevels(pheno$Designation),type = "random", K = 4, R = 1)
      Sample <- cbind(folds$which,folds$subsets)
      cv <- split(Sample[,2], f=Sample[,1])
    }
    sort[[a]] <- cv  
  }
  rm(a, folds, j, cv, Sample)
  
# Save outputs
  fold.list = sort
  results = list()
  Out = list()
# Get the model and run it
  ETA = list(list(X = geno, model = "BayesB"))
  for (z in 1:length(fold.list)) {
    for (i in 1:nfolds){
      
# Training set
      train_data <- pheno
# Validation set
      train_data[train_data$ID %in% fold.list[[z]][[i]], "blues_ns"] <- NA 
# Fitting model 
      BB_M <- BGLR(y = train_data$blues_ns, ETA = ETA, nIter = 10000, 
                   burnIn = 5000, thin = 5, verbose = F)
# GEBV
      Pred <- data.frame(Yhat = BB_M$yHat, G = pheno$ID)
      rownames(Pred) <- rownames(geno)
# Predicted GEBV
      results[[i]] <- Pred[Pred[, "G"] %in% fold.list[[z]][[i]], ] 
# Remove
      rm(BB_M, Pred, train_data)
    }
    GEBV <- do.call(rbind, results)
    GEBV <- GEBV[order(GEBV$G), ]
# Log
    log <- all(GEBV$G == pheno$ID)
# GET Results
    Out[[z]] <- data.frame(
      Rep = z,
      Log = log,
      Ac = round(cor(GEBV$Yhat, pheno$blues_ns, use = "na.or.complete"), 3),
      MSPE = round(mean(((GEBV$Yhat - pheno$blues_ns)^2), na.rm = T), 3)
    )
  }

