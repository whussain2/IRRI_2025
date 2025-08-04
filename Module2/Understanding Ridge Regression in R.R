###############################################################################
# Course: Fundamentals of Genomic Prediction-2025
# Module 2: Understanding Ridge Regression. 
# Authors: Waseem Hussain and  Mahender Anumalla 

###############################################################################


## -----------------------------------------------------------------------------------
# Matrix, 3 rows and 5 columns
  set.seed(1)
  n <- 3
  m <- 5
  X <- matrix(rbinom(n = n * m, size = 2, prob = 0.5), nrow = n, ncol = m)
  X


## -----------------------------------------------------------------------------------
  det(t(X) %*% X)


## -----------------------------------------------------------------------------------
# determinant
  det(t(X) %*% X + diag(1, m))


## -----------------------------------------------------------------------------------
# Load Package 
  library(BGLR)
  rm(list=ls()) # Remove previous history
# Read marker file in .ped format
  Geno<-read_ped("./Data/sativas413.ped")
# Set dimenions
  p=Geno$p
  n=Geno$n
  Geno=Geno$x
# Now load .fam file having genotype/acession names
  FAM <- read.table("./Data/sativas413.fam")
# Now read .mpa file containing map information 
  MAP <- read.table("./Data/sativas413.map")
# Let us now recode the marker data in .ped file
  Geno[Geno == 2] <- NA  # Converting missing data to NA
  #Geno[Geno == 0] <- 0   # Converting 0 data to 0
  #Geno[Geno == 1] <- 1   # Converting 1 to 1
  Geno[Geno == 3] <- 2   # Converting 3 to 2
# Now convert the marker data into matrix and transponse and check dimensions
  Geno <- matrix(Geno, nrow=p, ncol=n, byrow=TRUE)
  Geno <- t(Geno)
  dim(Geno)
  Geno[1:10, 1:10]
# Convert the missing, impute as mean
  for (j in 1:ncol(Geno)) {
    Geno[, j] <- ifelse(is.na(Geno[, j]), mean(Geno[, j], na.rm = TRUE), 
                        Geno[,j])
    }
  Geno[1:5, 1:4]
# Now assign the row and column names to marker file
  colnames(Geno)<-MAP$V2
# Adding line names stored in column second and pasted NSFTV_ID_ to each line name.
  row.names(Geno)<-paste0("NSFTV_",FAM$V2)
  #saveRDS(Geno,"./Data/Geno.coverted.rds")


## -----------------------------------------------------------------------------------
# Phenotypic data 
  pheno<-read.csv(file="./Data/rice.csv",
                   header=TRUE)
# Convert the missing data into mean
  for (j in 1:ncol(pheno)) {
    pheno[, j] <- ifelse(is.na(pheno[, j]), mean(pheno[, j], na.rm = TRUE), 
                         pheno[, j])
  }


## -----------------------------------------------------------------------------------
# Create a phenotype vector
  pheno1 <- as.vector(pheno$Flowering.time.at.Arkansas) # phenotype vector 
# Create intercept
  intercept <- rep(1, length(pheno1)) # intercept vector 
# Create an X matrix of SNPs (first 20) and add intercept
  X <- cbind(intercept, Geno[,1:20]) # the intercept and the SNP matrix including the first 100 SNPs
# Check dimesnions
  length(pheno1) 
  dim(X) 
# Fit OLS through equation
  ols1<- solve(t(X) %*% X) %*% t(X) %*% pheno1
  head(ols1) # estimates of the intercept and the first five SNP effects


## -----------------------------------------------------------------------------------
# Create subset of markers 1000
  X2 <- cbind(intercept, Geno[, 1:1000]) # the intercept and the whole SNP matrix 
  dim(X2) 
# Fit all Markers
  # use the solve() function
  #ols2<- solve(t(X2) %*% X2) %*% t(X2) %*% pheno1
  # Check the error
  # use the lm() function
  summary(lm(pheno1 ~ -1 + X2)) # check the warning 


## -----------------------------------------------------------------------------------
# Same 1000 markers
  X2 <- cbind(intercept, Geno[, 1:1000]) # the intercept and the whole SNP matrix 
  dim(X2) 
  # use the solve() function
  ols2<- solve(t(X2) %*% X2+diag(1, 1001)) %*% t(X2) %*% pheno1
# Check the error

