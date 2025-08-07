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


## ------------------------------------------------------------------------------------------
#png(file="~/boxplot.png", width=10, height =6, units = 'in',res=300) 
        boxplot<-ggplot(BLUEs.all, aes(x=Environment, y=BLUEs))+ 
          geom_boxplot(aes(fill=Environment))+ # fill by timepoint to give different color
          #scale_fill_manual(values = c("", ""))+
          #scale_color_manual(values = c("", ""))+
          theme_classic()+ #choose the theme for background
          labs(title="",x="Environment", y = "Grain Yield (Kg/ha)")+#Add the labels to the plots
          theme (plot.title = element_text(color="black", size=14, face="bold",hjust=0.5), # add and modify the title to plot
                 axis.title.x = element_text(color="black", size=12, face="bold"), # add and modify title to x axis
                 axis.title.y = element_text(color="black", size=12, face="bold")) + # add and modify title to y axis
          theme(axis.text= element_text(face = "bold", color = "black", size = 10))+ # modify the axis text
          theme(legend.position="none") # remove the theme from the plot
          #aes(x = fct_inorder(timepoint))+ # order the levels
          #dev.off()
         # ggplotly(p)
    boxplot


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


## ----message=FALSE, warning=FALSE----------------------------------------------------------

# Fit Model
  gs.model1<- mmes(BLUEs~Environment, # Environment Fixed
                  random= ~ vsm(ism(Genotype), Gu=GM), #vsm is covariance function to assign matrices
                  rcov= ~ units, # Residuals have no-covariance
                  data=BLUEs.all, verbose = FALSE)
# Get summary
  summary(gs.model1)

# Extract the GEBVs (random effects)
  estimated.all<-data.frame(GEBVs= gs.model1$u) # stored in u
  estimated.all$GEBVs<-estimated.all$GEBVs+ gs.model1$b[1] # adding intercept
  gs.model1$AIC # Check AIC and BIC values
  gs.model1$BIC
  kable(head(estimated.all)) # View as table


## ----message=FALSE, warning=FALSE----------------------------------------------------------
  gs.model3<- mmes(BLUEs~Environment,
                random= ~ vsm(dsm(Environment),ism(Genotype), Gu=GM),
                rcov= ~ units,
                data=BLUEs.all, verbose = FALSE)
  summary(gs.model3)


## ----eval=FALSE----------------------------------------------------------------------------
#   gs.model6<- mmes(BLUEs~Environment,
#                     random= ~ vsm(usm(Environment),ism(Genotype), Gu=GM),
#                     rcov= ~ units,
#                     data=BLUEs.all, verbose = FALSE)
#   summary(gs.model6)
# 


## ----eval=FALSE----------------------------------------------------------------------------
#     E <- diag(length(unique(BLUEs.all$Environment)))
#     rownames(E) <- colnames(E) <- unique(BLUEs.all$Environment)
#     Ei <- solve(E)
#     Gi <- solve(GM)
#     EGi <- kronecker(Ei,Gi, make.dimnames = TRUE)
#     Ei <- as(as(as( Ei,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
#     Gi <- as(as(as( Gi,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
#     EGi <- as(as(as( EGi,  "dMatrix"), "generalMatrix"), "CsparseMatrix")
#     attr(Gi, "inverse")=TRUE
#     attr(EGi, "inverse")=TRUE
#     model5<- mmes(BLUEs~Environment,
#                   random= ~ vsm(ism(Genotype), Gu=Gi) + vsm(ism(Environment:Genotype), Gu=EGi),
#                   rcov= ~ units,
#                   data=BLUEs.all, verbose = FALSE)
#     summary(gs.model5)


## ------------------------------------------------------------------------------------------
# Convert Environment to design matrix
  X <- model.matrix(~ Environment, data = BLUEs.all)
# Convert Genotype to design matrix
  Z <- model.matrix(~ Genotype - 1, data = BLUEs.all)
# Genotype kernel
  KG <- tcrossprod(Z %*%  GM)
# GxE interaction kernel
  KGE <- tcrossprod(X) * KG  # Element-wise multiplication

# ETA list for BGLR Package
    ETA <- list(
      ENV = list(X = X, model = "FIXED"),
      G   = list(K = KG, model = "RKHS"),
      GxE = list(K = KGE, model = "RKHS")
    )
# Fit the Model
modelGxE<- BGLR(y = BLUEs.all$BLUEs,
           ETA = ETA,
           nIter = 1000,
           burnIn = 100,
           thin = 2,
           verbose = FALSE)


## ------------------------------------------------------------------------------------------
# Variance components
  modelGxE$ETA$G$varU      # Genotype variance
  modelGxE$ETA$GxE$varU    # GxE variance
  modelGxE$varE            # Residual variance
# Breeding values (genotype effects)
  GEBVs <- data.frame(GEBVs_All=modelGxE$ETA$G$u)
  #BGLR directly outputs the genotype predictions as yHat
  GEBVs2<- data.frame(modelGxE$yHat)

