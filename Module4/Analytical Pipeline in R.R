###############################################################################
# Module 4: Analytical Pipeline in R
# Authors: Waseem Hussain and  Mahender Anumalla 
# Date: August 7, 2025
###############################################################################

# Load the Required Libraries
  rm(list=ls()) # Remove previous work
  library(rrBLUP)
  library(BGLR)
  library(AGHmatrix)
  library(ggplot2)
  library(DT)
  #library(cvTools)
  library(dplyr)
  library(lme4)
  library(arm)
  library(statgenSTA)


## --------------------------------------------------------------
demo.data.filtered<-read.csv(file="./Data/demo.data.filtered.csv",
                             header = TRUE)
# factor conversion if below are not in factors 
columns<-c("Environment", "Genotype", "Rep", "Block", "Row", "Column", "Line.type")
demo.data.filtered[, columns]<-lapply(columns, function(x) as.factor(demo.data.filtered[[x]]))
demo.data.filtered$Yield<-as.numeric(demo.data.filtered$Yield)
demo.data.filtered$HT<-as.numeric(demo.data.filtered$HT)
demo.data.filtered$DTF<-as.numeric(demo.data.filtered$DTF)

# Subset the required columns
demo.data.filtered<-demo.data.filtered[, c("Environment", "Genotype", "Rep", 
                                           "Block", "Row", "Column", "Line.type",
                                           "Yield", "HT", "DTF")]
# First we will arrange the rows and columns for spatial analysis.
# Now we will subset the environments and Yields for analysis
demo.data.filtered<-data.frame(demo.data.filtered%>% group_by(Environment)%>%arrange(Row, Column)) # arrange by row and column
demo.data.filtered<-data.frame(demo.data.filtered%>% arrange(Environment)) # Arrange by environment

#demo.data.filtered<-demo.data.filtered[!demo.data.filtered$Environment %in% c("Env2", "Env5","Env8", "Env9"), ]
# View as table in file
head(demo.data.filtered)


## --------------------------------------------------------------
ggplot(data = demo.data.filtered, aes(x = Environment, y = Yield, fill = Rep))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +# fill by timepoint to give different color
  #scale_fill_manual(values = c("", ""))+
  #scale_color_manual(values = c("", ""))
  theme (plot.title = element_text(color="black", size=12,hjust=0.5, face = "bold"), # add and modify the title to plot
         axis.title.x = element_text(color="black", size=12, face = "bold"), # add and modify title to x axis
         axis.title.y = element_text(color="black", size=12, face="bold")) + # add and modify title to y axis
  #scale_y_continuous(limits=c(0,15000), breaks=seq(0,15000,1000), expand = c(0, 0))+
  theme(axis.text= element_text(color = "black", size = 10)) # modify the axis text



## --------------------------------------------------------------
  ggplot(data = demo.data.filtered, aes(x = Environment, y = Yield, fill = Rep))+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +# fill by timepoint to give different color
  #scale_fill_manual(values = c("", ""))+
  #scale_color_manual(values = c("", ""))
  theme (plot.title = element_text(color="black", size=12,hjust=0.5, face = "bold"), # add and modify the title to plot
  axis.title.x = element_text(color="black", size=12, face = "bold"), # add and modify title to x axis
  axis.title.y = element_text(color="black", size=12, face="bold")) + # add and modify title to y axis
#scale_y_continuous(limits=c(0,15000), breaks=seq(0,15000,1000), expand = c(0, 0))+
 theme(axis.text= element_text(color = "black", size = 10)) # modify the axis text
  


## --------------------------------------------------------------
# Subset the environment 1
  sub.data<-subset(demo.data.filtered, Environment=="Env1")
  sub.data<-droplevels.data.frame(sub.data)


## --------------------------------------------------------------
# Now apply model
  model1<-lmer(Yield~Rep+(1|Genotype)+ (1|Rep:Block), data =sub.data)


## --------------------------------------------------------------
# Summarise the results
  summary(model1)


## --------------------------------------------------------------
  Ve<- VarCorr(model1)
  Ve


## --------------------------------------------------------------
# Plot the residual plot
  plot(fitted(model1), resid(model1), type="pearson")
  abline(0,0, col="blue")
# Plot QQ plot
  qqnorm(resid(model1))
# Residual plot
  plot(residuals(model1,type="pearson"), main='Model residuals', 
  ylab='Pearson residual value')


## --------------------------------------------------------------
# ANOVA
  anova(model1)


## --------------------------------------------------------------
  BLUEs<-fixef(model1)
  BLUEs


## --------------------------------------------------------------
# Extract the Random effects
  BLUPs<-data.frame(Blups.yield=ranef(model1)$Genotype)
  GV<-data.frame(BLUps.GY=coef(model1)$Genotype[,1]) #Genotype values (Blups +Intercept)


## --------------------------------------------------------------
# Extract the variance components
  Ve<- data.frame (VarCorr(model1))
  Ve
# Now calculate heritability using variance components
  genotype.var=Ve[1,4]
  error.var=Ve[2,4]
# Now heritability
  h2=genotype.var/(genotype.var+error.var)*100
  h2
# Reliability
  std.err<-se.ranef(model1)$Genotype
  v_BLUP<- mean(std.err)
# Heritability/Reliability 
  h2<- (1-((v_BLUP)^2/(Ve[1,4]*2)))*100
  h2


## --------------------------------------------------------------
# Run the analysis and check reliability
  # For Non-Stress Data using DTF as co-variate
  demo.data.filtered$Environment<- as.character(demo.data.filtered$Environment)
  un.exp<- unique(demo.data.filtered$Environment)
  for(i in 1:length(un.exp)){
    sub<- droplevels.data.frame(demo.data.filtered[which(demo.data.filtered$Environment==un.exp[i]),])

      model<-lmer(Yield~Rep+(1|Genotype)+ (1|Rep:Block), data =sub)
      #BLUPs<-data.frame(Blups.yield=ranef(model)$Genotype, Environment=un.exp[i])
        BLUPs<-data.frame(BLUps.GY=coef(model)$Genotype[,1], Environment=un.exp[i])
    if(i>1){
      BLUPs.all<-rbind(BLUPs.all, BLUPs)
    }
    else{
      BLUPs.all<- BLUPs
    }
  }
# Save the BLUES out put file for Genomic Predictions
 #estimates.all$Genotype<-gsub("^.{8}", "",  estimates.all$Genotype)


## --------------------------------------------------------------
# Linear model to get ANOVA
  demo.data.filtered$Environment<-as.factor(demo.data.filtered$Environment)
  model.anova<-lm(formula = Yield~Genotype+Environment+Genotype*Environment+Environment:Rep+ Environment:Rep:Block,
    data=demo.data.filtered)
# Get ANOVA
  anova(model.anova)


## --------------------------------------------------------------
# 
model2<- lmer(Yield~Rep+(1|Genotype)+(1|Environment)+
          (1|Environment:Rep)+(1|Environment:Rep:Block),
          data=demo.data.filtered)

#plot residuals 
plot(residuals(model2,type="pearson"), main='Model residuals', 
ylab='Pearson residual value')
#var.test(Yield~Environment,data=demo.data.filtered)


## --------------------------------------------------------------
demo.data.filtered$Environment<-as.factor(demo.data.filtered$Environment)
Model3.lme4<-lmer(Yield~Genotype+(1|Rep)+(1|Environment:Genotype)+
          +(1|Environment:Rep:Block), data=demo.data.filtered)



## --------------------------------------------------------------
  summary(Model3.lme4)


## --------------------------------------------------------------
plot(Model3.lme4)


## --------------------------------------------------------------
Ve<- data.frame (VarCorr(Model3.lme4))
Ve


## --------------------------------------------------------------
#std.err<-se.ranef(Model3.lme4)$Genotype
#v_BLUP<- mean(std.err)
# Heritability/Reliability 
#h2<- (1-((v_BLUP)^2/(Ve[2,4]*2)))*100
#h2


## --------------------------------------------------------------
# BLUEs
  BLUEs.all<-data.frame(BLUEs.Yield=fixef(Model3.lme4))
  #BLUPs<-data.frame(BLUps.GY=coef(Model3.lme4)$Genotype[,1])
  head(BLUEs.all)


## --------------------------------------------------------------
library(statgenSTA)
TD_STA <- createTD(sub.data, 
                   trial = "Environment",
                   genotype = "Genotype", 
                   rowCoord = "Row", 
                   colCoord = "Column",
                   repId = "Rep",
                   subBlock = "Block",
                   trDesign = "res.ibd")

## Create layout plot with variety labels
plot(TD_STA, 
        plotType = "layout",
        showGeno = TRUE)


## Model specification (using engine = "SpATS")
sta_model_SpATS <- fitTD(TD = TD_STA, 
                    traits = "Yield",
                    design = "res.ibd",
                    what = "fixed",
                    spatial = TRUE,
                    engine = "SpATS")

plot(sta_model_SpATS, 
     plotType = "spatial", 
     traits = 'Yield')

## Extract all available statistics from the fitted model.
  extr <- extractSTA(sta_model_SpATS)
## Extract only the BLUEs from the fitted model.
  BLUEs <- extractSTA(sta_model_SpATS,
                    what = "BLUEs")

