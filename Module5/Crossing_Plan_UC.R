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


## ------------------------------------------------------------------------------------------
# Install Library
  #library(devtools)
  #install_github('Resende-Lab/SimpleMating')
  library(SimpleMating)
  library(AGHmatrix)


## ------------------------------------------------------------------------------------------
# 1.  Marker Data
  geno<-read.csv(file="./Data/Geno_Marker.csv", header=TRUE,
              row.names = 1)
  geno<-as.matrix(geno)
  dim(geno)
#kable(head(geno))

# 2. Markers effects for target trait
  MarEff<-read.csv(file="./Data/Marker_effects.csv", header=TRUE,
                 row.names = 1)
  #kable(head(MarEff))
# 3. Genetic Map or LD Matrix
  Map_info<-read.csv(file="./Data/Map_info.csv", header = TRUE,
                   row.names = 1)
  #kable(head(Map_info))


## ------------------------------------------------------------------------------------------
#GRM<- (lines_Geno %*% t(lines_Geno)) / ncol(lines_Geno)
  GRM<- Gmatrix(geno, missingValue=NA, 
             maf=0.05, method="VanRaden")
  dim(GRM)
  heatmap(GRM)


## ------------------------------------------------------------------------------------------
# Select parents from geno
Parents <- rownames(geno)


## ------------------------------------------------------------------------------------------
# Create a Half sib plan
Cross_plan <- planCross(TargetPop = Parents,
                        MateDesign = "half") # see other plans in arguments


## ------------------------------------------------------------------------------------------
# Usefulness of trait number 1
  UC_crosses<- getUsefA(MatePlan = Cross_plan,
                     Markers = geno,
                     addEff = MarEff[, 1], # Trait with Marker effects
                     Map.In = Map_info,
                     K = GRM,
                     propSel = 0.05, # Proportion selected
                     Type = "RIL", # Type DH or RILs
                     Generation = 6) # Which generation
kable(head(UC_crosses[[1]], 10))


## ------------------------------------------------------------------------------------------
# Use Usefulness Data above
  maxGainPlan <- selectCrosses(data =  UC_crosses[[2]], 
                n.cross = 20, # Set Target crosses
                max.cross = 2, # Maximum crosses parent can be part
                min.cross = 1, # Minimum crosses parent can be part
                culling.pairwise.k = 1) # covariance between cross, restriction
# Mating Plan
maxGainPlan[[2]]

