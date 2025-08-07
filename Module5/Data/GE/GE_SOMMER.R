


data(DT_example)
DT <- BLUEs.all
A <- GM


# Using subset()
subset_envs <- c("Env1", "Env2", "Env3", "Env4")
BLUES.subset <- subset(BLUEs.all, Environment %in% subset_envs)

# Single Environment Model

gs.model1<- mmes(BLUEs~1,
                  random= ~ vsm(ism(Genotype), Gu=GM),
                  rcov= ~ units,
                  data=BLUEs.all, verbose = FALSE)
summary(gs.model1)

estimated.all<-data.frame(GEBVs= gs.model1$u)
estimated.all$GEBVs<-estimated.all$GEBVs+ gs.model1$b
kable(head(estimated.all)) # S


# 2 MET: main effect model

gs.model2<- mmes(BLUEs~Environment,
                random= ~ vsm(ism(Genotype), Gu=GM),
                rcov= ~ units,
                data=BLUEs.all, verbose = FALSE)

summary(gs.model2)


estimated.all<-data.frame(GEBVs= gs.model1$u)
estimated.all$GEBVs<-estimated.all$GEBVs+ gs.model1$b
gs.model1$AIC
gs.model1$BIC
kable(head(estimated.all)) # S


