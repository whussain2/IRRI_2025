###############################################################################
# Module 3: G x E in R. 
# Authors: Waseem Hussain and  Mahender Anumalla 
# Date: August 6, 2024
###############################################################################

# Load the Required Libraries
  rm(list=ls()) # Remove previous work
  library(lme4)
  library(sommer)
  library(AGHmatrix)
  library(ggplot2)
  library(DT)
  library(statgenSTA)
  library(statgenGxE)
  library(metan)
  library(readxl)
  library(arm)


## ---------------------------------------------------------------------------
# Read the data
  demo.data<- read_excel("./Data/demo.data.xlsx", sheet ="Filtered")
# Convert variables into appropriate data types
  demo.data$Genotype<-as.factor(demo.data$Genotype) # Genotypes as factor
  demo.data$Block<-as.factor(demo.data$Block) # Block as factor
  demo.data$Row<-as.factor(demo.data$Row) # Row as factor
  demo.data$Rep<-as.factor(demo.data$Rep) # Replication as factor
  demo.data$Column<-as.factor(demo.data$Column) # Column as factor
  demo.data$Environment<-as.factor(demo.data$Environment)  # Env. as factor
  demo.data$Yield<-as.numeric(demo.data$Yield)  # Env. as factor
# Order the factor levels
  demo.data$Environment<- factor(demo.data$Environment, 
                        levels = c("Env1","Env2", "Env3",  "Env4", "Env5","Env6",
                                   "Env7",  "Env8",  "Env9",  "Env10", "Env11"))
# Visualize as table
  kable(head(demo.data))
  str(demo.data)


## ----message=FALSE, warning=FALSE-------------------------------------------
  ggplot(demo.data, aes(x=Environment, y=Yield))+ 
  geom_boxplot(color="black", fill="pink")+
  # fill by timepoint to give different color
  #scale_fill_manual(values = c("", ""))+
  #scale_color_manual(values = c("", ""))+
  theme_classic()+ #choose the theme for background
  labs(title="",x="Environments", y = "Grain Yield")+#Add the labels to the plots
  theme (plot.title = element_text(color="black", size=14, 
        face="bold",hjust=0.5), # add and modify the title to plot
  axis.title.x = element_text(color="black", size=12),
  axis.title.y = element_text(color="black", size=12, 
      face="bold")) + # add and modify title to y axis
  theme(axis.text= element_text(color = "black", size = 10))+ # modify the axis text
  theme(legend.position="none") # remove the theme from the plot
  


## ----warning=FALSE, message=FALSE-------------------------------------------
# Run linear model using lm function of R
model1<-lm(formula = 
Yield~Genotype+Environment+Genotype*Environment+Environment:Rep+ Environment:Rep:Block,
          data=demo.data)
# Get ANNOVA using anova function
anova(model1)
#test<-data.frame(model1$coefficients)


## ---------------------------------------------------------------------------
## Create a TD object from dropsPheno.
  dropsss<- statgenSTA::createTD(data =demo.data, genotype = "Genotype", 
                                 trial = "Environment")
## Perform a Finlay-Wilkinson analysis for all trials.
  dropsFW <- gxeFw(TD = dropsss, trait = "Yield")
## Create line plot for Finlay Wilkinson analysis.
  plot(dropsFW, plotType = "line")


## ---------------------------------------------------------------------------
# Developm AMMI Model
Model_ammi<-performs_ammi(demo.data,
              env = Environment,
              gen = Genotype,
              rep = Rep,
              block=Block,
              resp = Yield)
# ANNOVA
  Model_ammi$Yield$ANOVA
# PCA Components
  Model_ammi$Yield$PCA

# AMMI  Biplot
  plot_scores(Model_ammi, type = 1)
# Stability indices
  stab_indexes <- AMMI_indexes(Model_ammi)
  head (stab_indexes$Yield)[1:5,1:7]


## ---------------------------------------------------------------------------
gge_model <- gge(demo.data,
                 env = Environment,
                 gen = Genotype,
                 resp = Yield,
                 centering = 2,
                 svp = "symmetrical", 
                 scaling = 0)

# Mean performance vs. stability 
  plot(gge_model, type = 2, size.text.gen = 2)
# Discriminativeness vs. representativeness
  plot(gge_model, type = 4, size.text.gen = 3)



## ---------------------------------------------------------------------------
demo.data$Environment<-as.factor(demo.data$Environment)
Model3.lme4<-lmer(Yield~Rep+(1|Genotype)+(1|Environment:Genotype)+
          +(1|Environment:Rep:Block), data=demo.data)



## ---------------------------------------------------------------------------
  summary(Model3.lme4)


## ---------------------------------------------------------------------------
  plot(Model3.lme4)


## ---------------------------------------------------------------------------
  Ve<- data.frame (VarCorr(Model3.lme4))
  Ve


## ---------------------------------------------------------------------------
  std.err<-se.ranef(Model3.lme4)$Genotype
  v_BLUP<- mean(std.err)
# Heritability/Reliability 
  h2<- (1-((v_BLUP)^2/(Ve[2,4]*2)))*100
  h2


## ---------------------------------------------------------------------------
  BLUPs<-data.frame(BLUps.GY=coef(Model3.lme4)$Genotype[,1])
  head(BLUPs)


## ---------------------------------------------------------------------------
library(statgenSTA)
TD_STA <- createTD(demo.data, 
                   trial = "Environment",
                   genotype = "Genotype", 
                   rowCoord = "Row", 
                   colCoord = "Column",
                   repId = "Rep",
                   subBlock = "Block",
                   trDesign = "res.ibd")

## Create layout plot with variety labels
#plot(TD_STA, 
       # plotType = "layout",
        #showGeno = TRUE)
## Model specification (using engine = "SpATS")
sta_model_SpATS <- fitTD(TD = TD_STA, 
                    traits = "Yield",
                    design = "res.ibd",
                    what = "fixed",
                    spatial = TRUE,
                    engine = "SpATS")

#plot(sta_model_SpATS, 
   #  plotType = "spatial", 
    # traits = 'Yield')

## Extract all available statistics from the fitted model.
  extr <- extractSTA(sta_model_SpATS)
## Extract only the BLUEs from the fitted model.
  BLUEs <- extractSTA(sta_model_SpATS,
                    what = "BLUEs")

