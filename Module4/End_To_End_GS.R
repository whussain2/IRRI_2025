###############################################################################
# Module 4: End_to_End_GS Pipeline in R
# Authors: Waseem Hussain and  Mahender Anumalla 
# Date: August 7, 2025
###############################################################################

# Load the Required Libraries
  rm(list=ls()) # Remove previous work
  library(rrBLUP)
  library(BGLR)
  library(AGHmatrix)
  #library(cvTools)
  library(dplyr)
  library(arm)
  library(ggplot2)
  library(DT)
  #library(cvTools)
  library(lme4)
  library(reshape2)
  library(data.table)
  library(sommer)


## --------------------------------------------------------------
# Read the phenotypic data
  pheno<-read.csv(file="./Data/RAW.DATA1.csv", header=TRUE) 
  str(pheno)
  pheno$Environment<-as.factor(pheno$Environment) # Convert environment as factor
  pheno$Season<-as.factor(pheno$Season) # Convert Season as factor
  pheno$Genotype<-as.factor(pheno$Genotype) # Convert Genotype as factor
  pheno$Block<-as.factor(pheno$Block) # Convert Block as factor
  pheno$Replication<-as.factor(pheno$Replication) # Convert Replcation as factor
  #table(pheno$Environment)


## --------------------------------------------------------------
 str(pheno)


## --------------------------------------------------------------
# Check the missing data in each location
   Data.missing<-data.frame(pheno %>%group_by(Environment) %>%
                        summarise_each(funs(sum(is.na(.))/length(.))))
  # Extract the three variables
  Data.missing<-Data.missing[, c("Environment", "Yield", "Replication", "Block")]
  Data.missing<-melt(setDT(Data.missing), id.vars = c("Environment"), 
                     variable.name = "Trait")
    
# Plot the missing plot for Grain Yield
   #png(file = "./Outputs/Plots/Missing.data.png", width =12, 
    # height =6, units = "in", res = 600)
    ggplot(Data.missing, aes(x=Environment, y=value))+ 
     geom_point(size=3) + 
     geom_segment(aes(x=Environment, 
                             xend=Environment, 
                             y=0, 
                            yend=value)) + 
      labs(title="", y="Proportion of Missing Data", x="Environment" )+
     theme_classic()+
     theme(axis.text.x = element_text(angle=90, vjust=0.6))+
     #gghighlight(max(value) > .05, label_key =Environment)+
     facet_wrap(~Trait , ncol = 3,nrow=1,scales = "free")+
     theme (plot.title = element_text(color="black", size=14, hjust=0.5),
                  axis.title.x = element_text(color="black", size=24),
                  axis.title.y = element_text(color="black", size=24))+
     theme(axis.text= element_text(color = "black", size = 10))


## ----warning=FALSE---------------------------------------------
# Get the Box plot
  ggplot(pheno, aes(x=Environment, y=Yield))+ 
    geom_boxplot(aes(fill=Environment))+
    theme_classic()+
    geom_boxplot()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +# fill by timepoint to give different color
     #scale_fill_manual(values = c("", ""))+
    #scale_color_manual(values = c("", ""))
    theme (plot.title = element_text(color="black", size=12,hjust=0.5, face = "bold"), # add and modify the title to plot
          axis.title.x = element_text(color="black", size=12, face = "bold"), # add and modify title to x axis
          axis.title.y = element_text(color="black", size=12, face="bold"))  # add and modify title to y axis
     #scale_y_continuous(limits=c(0,15000), breaks=seq(0,15000,1000), expand = c(0, 0)) theme(legend.title = element_text(colour="black", size=16), legend.position = "none",




## --------------------------------------------------------------
  summary.Yield<-data.frame(pheno %>% 
                group_by(Environment)%>% 
                summarize(Mean = mean(Yield, na.rm=TRUE),
                Median= median(Yield, na.rm=TRUE),
                SD =sd(Yield, na.rm=TRUE),
                Min.=min(Yield, na.rm=TRUE),
                Max.=max(Yield, na.rm=TRUE),
                CV=sd(Yield, na.rm=TRUE)/mean(Yield, na.rm=TRUE)*100,
                St.err= sd(Yield, na.rm=TRUE)/sqrt(length(Yield))
                ))
  summary.Yield


## --------------------------------------------------------------
# Use lme4 R package
  pheno$Environment<- as.character(pheno$Environment) # Environment as character for use for loop
  un.exp<- unique(pheno$Environment) # unique environments
  for(i in 1:length(un.exp)){ # use for loop from 1 to n environments
    sub<- droplevels.data.frame(pheno[which(pheno$Environment==un.exp[i]),]) # Run per environment
    if (length(levels(sub$Season))>1){ # If season is more than one fit below model
      
      model<-lmer(Yield ~ Genotype+ (1|Season)+ (1|Block), data=sub)

    } else { # if season is just one fit then this model
      model<-lmer(Yield ~ Genotype + (1|Replication), data=sub)
    }
    estimates<-data.frame(BLUEs=fixef(model)[-1], Environment=un.exp[i]) # Extract BLUEs
    estimates$BLUEs<-estimates$BLUEs+fixef(model)[1] # Add intercept
    estimates$Genotype<-row.names(estimates) # Add names
    if(i>1){
      BLUEs.all<-rbind(BLUEs.all, estimates) # combine all across environments
    }
    else{
      BLUEs.all<- estimates
    }
  }
# Save the BLUES out put file for Genomic Predictions
  BLUEs.all$Genotype<-gsub("^.{8}", "",  BLUEs.all$Genotype)
 # Save the file in csv formate
 #write.csv( BLUEs.all, file="BLUES.ALL.csv")


## --------------------------------------------------------------
  head(BLUEs.all)


## --------------------------------------------------------------
  BLUEs.all<-read.csv(file="./Data/BLUES.ALL.csv")


## --------------------------------------------------------------
  geno<-readRDS("./Data/GBS_datav2.rds")
  dim(geno)
# Match genotype with Phenotype
  Ids<-unique(BLUEs.all$Genotype)
  length(Ids)
# Now subet the genotype Data based on IDs
  geno<-geno[row.names(geno)%in%Ids,]
  dim(geno)


## --------------------------------------------------------------
  GM<- Gmatrix(SNPmatrix=geno, missingValue=NA, 
                            maf=0.05, method="VanRaden")
  dim(GM)


## --------------------------------------------------------------
  #heatmap(GM)


## --------------------------------------------------------------
# Write the model
  g_blup<- mmer(BLUEs~1, # BLUES as response variable name of column in data
              random=~vsr(Genotype,Gu=GM)+Environment, # vsr function take covariance structure
              rcov=~units, nIters=3,data=BLUEs.all,verbose = FALSE) 
  summary(g_blup)
  #BLUPs<-g_blup$U$`u:Genotype`$BLUEs
 # g_blup$Beta[1,3]
  estimated<-data.frame(GEBVs= g_blup$U$`u:Genotype`$BLUEs) # Extract the Random effects
  estimated$GEBVs<-estimated$GEBVs+ g_blup$Beta[1,3] # Add intercept (mean)
  head( estimated) # Show in Table


## ----warning=FALSE---------------------------------------------
# Get the Box plot
  ggplot(estimated, aes(y = GEBVs)) +
      geom_boxplot(fill="skyblue")+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +# fill by timepoint to give different color
     #scale_fill_manual(values = c("", ""))+
    #scale_color_manual(values = c("", ""))
    theme (plot.title = element_text(color="black", size=12,hjust=0.5, face = "bold"), # add and modify the title to plot
          axis.title.x = element_text(color="black", size=12, face = "bold"), # add and modify title to x axis
          axis.title.y = element_text(color="black", size=12, face="bold"))  # add and modify title to y axis
     #scale_y_continuous(limits=c(0,15000), breaks=seq(0,15000,1000), expand = c(0, 0))



## --------------------------------------------------------------
# Read the same marker data
  geno_all<-readRDS("./Data/GBS_datav2.rds")
  dim(geno_all)


## --------------------------------------------------------------
  GM_all<- Gmatrix(SNPmatrix=geno_all, missingValue=NA, 
                            maf=0.05, method="VanRaden")
  dim(GM_all)
  heatmap(GM_all)


## --------------------------------------------------------------
  gs.model1<- mmer(BLUEs~1,
              random=~vsr(Genotype,Gu=GM_all)+Environment,
              rcov=~units, nIters=3,data=BLUEs.all,verbose = FALSE) 
  summary(gs.model1)
  #gs.model1$U$`u:Genotype`$BLUEs
  #gs.model1$Beta[1,3]
  estimated.all<-data.frame(GEBVs= gs.model1$U$`u:Genotype`$BLUEs)
  estimated.all$GEBVs<-estimated.all$GEBVs+ gs.model1$Beta[1,3]
   kable(head(estimated.all)) # Show in Table


## --------------------------------------------------------------
ggplot(data=estimated, aes(GEBVs))+
          #geom_density(alpha = 0.5)+
          geom_histogram(fill="pink", color="black")+
          #theme_few()+ #use white theme
          labs(title="",x="Value", y = "Count")



## --------------------------------------------------------------
#var comps
  sm <- summary( gs.model1)$varcomp
  sm


## --------------------------------------------------------------
  library(sommer)
  vg <- sm[grepl('Genotype', row.names(sm)), 1]
  ve <- sm[grepl('units', row.names(sm)), 1]
  #hertability <- vpredict(gs.model1, h2 ~ V1 / (V1 + V2))
  #hertability
  #h2 <- vpredict(gs.model1,dam.prop~V1/ (V1 + V2+V3))
  #h2
  #h2<-vg/(vg+ve)
  #h2


## --------------------------------------------------------------
# Get prediction Error Variance 
  pev <- diag(gs.model1$PevU$`u:Genotype`$BLUEs)
# Get Reliability
  reliability<- data.frame(r2=1 - pev / vg) # Recall the formulla of reliability
  head(reliability)


## --------------------------------------------------------------
  estimated.all$Genotype<-rownames(estimated.all) # Assign column to genotypes
  reliability$Genotype<-rownames(reliability) # Assign names to column
  GEBVs.all<-merge(estimated.all, reliability, by="Genotype")
# Arrange the BLUPs in decreasing order
  GEBVs.all<-GEBVs.all%>%arrange(desc(GEBVs))
# Now select Top 40
  GEBVs.top50<-data.frame(GEBVs.all[1:50, ])
  kable( GEBVs.top50)
  #write.csv(GEBVs.all, file="./Data/GEBVs.all.csv")


## --------------------------------------------------------------
# Miantain order from top to bottom
  GEBVs.top50$Genotype <- factor(GEBVs.top50$Genotype, 
                                 levels=unique(GEBVs.top50$Genotype))
# Visualize as bar plot
  bar.plot<-ggplot(data=GEBVs.top50, aes(x=Genotype, y=GEBVs)) +
  geom_bar(stat="identity", width=0.5, fill="blue")+
  theme_classic()+
  labs(title="GEBVs of Top Ranked Genotypes",x="Genotype", y = "Breeding Value")+
  #scale_y_continuous(limits = c(0, 6000), breaks = seq(0, 6000, by = 500))+
  theme (plot.title = element_text(color="black", size=1, face="bold", hjust=0),
           axis.title.x = element_text(color="black", size=10, face="bold"),
           axis.title.y = element_text(color="black", size=10, face="bold")) +
   theme(axis.text= element_text(color = "black", size = 8))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
 bar.plot

