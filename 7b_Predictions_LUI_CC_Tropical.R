
# Here, values for certain fixed effect combinations are predicted and 
# expressed as percentage change for use in the text of the paper. 


# directories
predictsDir <- "C:/Users/Kyra/Documents/GLITRS/Code/5_RunLUIClimateModels/"
modDir <- "C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/"
outDir <- "C:/Users/Kyra/Documents/GLITRS/Code/7_Predictions/"

if(!dir.exists(outDir)) dir.create(outDir)

# load libraries
library(StatisticalModels)
library(predictsFunctions)
library(webshot)
library(gt)
source('C:/Users/Kyra/Documents/GLITRS/Data/0_Functions.R')


# read in the predicts data
predictsSites <- readRDS(file = paste0(predictsDir,"PREDICTSSitesClimate_Data.rds"))
trop <- readRDS(file = paste0(modDir,"trop.rds"))
nontrop <- readRDS(file = paste0(modDir,"nontrop.rds"))


#### Hyp 1: land use effect only ####

# see LUI_plots.R

#### Hyp 2: Land use and climate anomaly interaction ####

# first one, looking at SCA of 1

# load in models
load(file = paste0(moddir, "/MeanAnomalyModelAbund_trop.rdata"))
load(file = paste0(moddir, "/MeanAnomalyModelAbund_nontrop.rdata"))
load(file = paste0(moddir, "/MeanAnomalyModelRich_trop.rdata"))
load(file = paste0(moddir, "/MeanAnomalyModelRich_nontrop.rdata"))
load(file = paste0(moddir, "/MaxAnomalyModelAbund_trop.rdata"))
load(file = paste0(moddir, "/MaxAnomalyModelAbund_nontrop.rdata"))
load(file = paste0(moddir, "/MaxAnomalyModelRich_trop.rdata"))
load(file = paste0(moddir, "/MaxAnomalyModelRich_nontrop.rdata"))

## Mean Anomaly, Trop ##
#create matrix for predictions 
# SCA = 1
# abun and richness = 0

# what is the rescaled value of SCA of 1
BackTransformCentreredPredictor(transformedX = 0.27, originalX = trop$StdTmeanAnomaly) # 0.27 gives about 1 

# what is the rescaled value of SCA of 0
BackTransformCentreredPredictor(transformedX = -1.58, originalX = trop$StdTmeanAnomaly) # -1.58 gives about 0 

# these values might not be close enough to 1 or 0, and are leading to mismatches in the values in Predictions and the values on the plot

# reference is primary with 0 climate change so have 0 for that row

data_tab <- expand.grid(LUI = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High","Primary vegetation","Secondary vegetation", "Agriculture_Low", "Agriculture_High"),
                        Order = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera"),
                        LogAbund = 0,
                        Species_richness = 0)

# add column with SCA, values repeating 6 times
StdTmeanAnomalyRS = rep(c(-1.58,-1.58,-1.58,-1.58,0.27,0.27,0.27,0.27),times=6)

# add SCA to the data_tab
data_tab<-cbind(data_tab,StdTmeanAnomalyRS)

# factor the LU info
data_tab$LUI <- factor(data_tab$LUI, levels = levels(predictsSites$LUI))

## Abundance, Tropical ##

# predict results
result.ab.trop <- PredictGLMER(model = MeanAnomalyModelAbund_trop$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.ab.trop <- exp(result.ab.trop)-0.01

# add in the LU info
result.ab.trop$LUI <- data_tab$LUI

# add in the Order info
result.ab.trop$Order <- data_tab$Order

# break into Orders
Order<- paste0("",result.ab.trop$Order)
list.result.ab.trop <- split(result.ab.trop,Order)
list2env(list.result.ab.trop,globalenv())

# express as a percentage of primary
Coleoptera$perc <- ((Coleoptera$y/Coleoptera$y[1]) * 100) - 100
Diptera$perc <- ((Diptera$y/Diptera$y[1]) * 100) - 100
Hemiptera$perc <- ((Hemiptera$y/Hemiptera$y[1]) * 100) - 100
Hymenoptera$perc <- ((Hymenoptera$y/Hymenoptera$y[1]) * 100) - 100
Lepidoptera$perc <- ((Lepidoptera$y/Lepidoptera$y[1]) * 100) - 100
Orthoptera$perc <- ((Orthoptera$y/Orthoptera$y[1]) * 100) - 100

# put it back together
result.ab.trop <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Orthoptera)

# add in SCA vals
result.ab.trop$SCA <- rep(c(0,0,0,0,1, 1, 1, 1),times=6)

# add zone factor
result.ab.trop$Zone <- as.factor("Tropical")

## Species Richness, Tropical ##

# predict the results
result.sr.trop <- PredictGLMER(model = MeanAnomalyModelRich_trop$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.sr.trop <- exp(result.sr.trop)

# add in the LU info
result.sr.trop$LUI <- data_tab$LUI

# add in the Order info
result.sr.trop$Order <- data_tab$Order

# break into Orders
Order<- paste0("",result.sr.trop$Order)
list.result.sr.trop <- split(result.sr.trop,Order)
list2env(list.result.sr.trop,globalenv())

# express as a percentage of primary
Coleoptera$perc <- ((Coleoptera$y/Coleoptera$y[1]) * 100) - 100
Diptera$perc <- ((Diptera$y/Diptera$y[1]) * 100) - 100
Hemiptera$perc <- ((Hemiptera$y/Hemiptera$y[1]) * 100) - 100
Hymenoptera$perc <- ((Hymenoptera$y/Hymenoptera$y[1]) * 100) - 100
Lepidoptera$perc <- ((Lepidoptera$y/Lepidoptera$y[1]) * 100) - 100
Orthoptera$perc <- ((Orthoptera$y/Orthoptera$y[1]) * 100) - 100

# put it back together
result.sr.trop <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Orthoptera)

# add in SCA vals
result.sr.trop$SCA <- rep(c(0,0,0,0,1, 1, 1, 1),times=6)

# add zone factor
result.sr.trop$Zone <- as.factor("Tropical")

## Mean Anomaly, NonTrop ##
#create matrix for predictions 
# SCA = 1
# abun and richness = 0

# what is the rescaled value of SCA of 1
BackTransformCentreredPredictor(transformedX = 2.2, originalX = nontrop$StdTmeanAnomaly) # -1.77 gives about 1 

# what is the rescaled value of SCA of 0
BackTransformCentreredPredictor(transformedX = -1.77, originalX = nontrop$StdTmeanAnomaly) # -1.77 gives about 0 

# these values might not be close enough to 1 or 0, and are leading to mismatches in the values in Predictions and the values on the plot

# reference is primary with 0 climate change so have 0 for that row

data_tab <- expand.grid(LUI = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High","Primary vegetation","Secondary vegetation", "Agriculture_Low", "Agriculture_High"),
                        Order = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera"),
                        LogAbund = 0,
                        Species_richness = 0)

# add column with SCA, values repeating 6 times
StdTmeanAnomalyRS = rep(c(-1.77,-1.77,-1.77,-1.77,2.2,2.2,2.2,2.2),times=6)

# add SCA to the data_tab
data_tab<-cbind(data_tab,StdTmeanAnomalyRS)

# factor the LU info
data_tab$LUI <- factor(data_tab$LUI, levels = levels(predictsSites$LUI))


## Abundance, NonTropical ##

# predict results
result.ab.nontrop <- PredictGLMER(model = MeanAnomalyModelAbund_nontrop$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.ab.nontrop <- exp(result.ab.nontrop)-0.01

# add in the LU info
result.ab.nontrop$LUI <- data_tab$LUI

# add in the Order info
result.ab.nontrop$Order <- data_tab$Order

# break into Orders
Order<- paste0("",result.ab.nontrop$Order)
list.result.ab.nontrop <- split(result.ab.nontrop,Order)
list2env(list.result.ab.nontrop,globalenv())

# express as a percentage of primary
Coleoptera$perc <- ((Coleoptera$y/Coleoptera$y[1]) * 100) - 100
Diptera$perc <- ((Diptera$y/Diptera$y[1]) * 100) - 100
Hemiptera$perc <- ((Hemiptera$y/Hemiptera$y[1]) * 100) - 100
Hymenoptera$perc <- ((Hymenoptera$y/Hymenoptera$y[1]) * 100) - 100
Lepidoptera$perc <- ((Lepidoptera$y/Lepidoptera$y[1]) * 100) - 100
Orthoptera$perc <- ((Orthoptera$y/Orthoptera$y[1]) * 100) - 100

# put it back together
result.ab.nontrop <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Orthoptera)

# add in SCA vals
result.ab.nontrop$SCA <- rep(c(0,0,0,0,1, 1, 1, 1),times=6)

# add zone factor
result.ab.nontrop$Zone <- as.factor("NonTropical")

## Species Richness, NonTropical ##

# predict the results
result.sr.nontrop <- PredictGLMER(model = MeanAnomalyModelRich_nontrop$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.sr.nontrop <- exp(result.sr.nontrop)

# add in the LU info
result.sr.nontrop$LUI <- data_tab$LUI

# add in the Order info
result.sr.nontrop$Order <- data_tab$Order

# break into Orders
Order<- paste0("",result.sr.nontrop$Order)
list.result.sr.nontrop <- split(result.sr.nontrop,Order)
list2env(list.result.sr.nontrop,globalenv())

# express as a percentage of primary
Coleoptera$perc <- ((Coleoptera$y/Coleoptera$y[1]) * 100) - 100
Diptera$perc <- ((Diptera$y/Diptera$y[1]) * 100) - 100
Hemiptera$perc <- ((Hemiptera$y/Hemiptera$y[1]) * 100) - 100
Hymenoptera$perc <- ((Hymenoptera$y/Hymenoptera$y[1]) * 100) - 100
Lepidoptera$perc <- ((Lepidoptera$y/Lepidoptera$y[1]) * 100) - 100
Orthoptera$perc <- ((Orthoptera$y/Orthoptera$y[1]) * 100) - 100

# put it back together
result.sr.nontrop <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Orthoptera)

# add in SCA vals
result.sr.nontrop$SCA <- rep(c(0,0,0,0,1, 1, 1, 1),times=6)

# add zone factor
result.sr.nontrop$Zone <- as.factor("NonTropical")

# combine results into a table for saving
all_res <- rbind(result.ab.nontrop, result.ab.trop, result.sr.nontrop, result.sr.trop)

all_res$measure <- c(rep("ab", 96), rep("sr", 96))

# save as png
percentage_change_LUI_CC <- all_res %>% gt()
gtsave(percentage_change_LUI_CC,"C:/Users/Kyra/Documents/GLITRS/Code/7_Predictions/MeanAnom_PercentageChange_LUI_CC_Tropical.png")


# save table as csv
write.csv(all_res, file = paste0(outDir, "/MeanAnom_PercentageChange_LUI_CC_Tropical.csv"))

## Max Anomaly , Tropical ##
#create matrix for predictions 
# SCA = 1
# abun and richness = 0

# what is the rescaled value of SCA of 1
BackTransformCentreredPredictor(transformedX = -0.425, originalX = trop$StdTmaxAnomaly) # -0.425 gives about 1 

# what is the rescaled value of SCA of 0
BackTransformCentreredPredictor(transformedX = -1.035, originalX = trop$StdTmaxAnomaly) # -1.035 gives about 0

# reference is primary with 0 climate change so have 0 for that row

data_tab <- expand.grid(LUI = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High","Primary vegetation","Secondary vegetation", "Agriculture_Low", "Agriculture_High"),
                        Order = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera"),
                        LogAbund = 0,
                        Species_richness = 0)

# add column with SCA, values repeating 6 times
StdTmaxAnomalyRS = rep(c(-1.035,-1.035,-1.035,-1.035,-0.425,-0.425,-0.425,-0.425),times=6)

# add SCA to the data_tab
data_tab<-cbind(data_tab,StdTmaxAnomalyRS)

# factor the LU info
data_tab$LUI <- factor(data_tab$LUI, levels = levels(predictsSites$LUI))

## Abundance, Tropical ##

# predict results
result.ab.trop <- PredictGLMER(model = MaxAnomalyModelAbund_trop$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.ab.trop <- exp(result.ab.trop)-0.01

# add in the LU info
result.ab.trop$LUI <- data_tab$LUI

# add in the Order info
result.ab.trop$Order <- data_tab$Order

# break into Orders
Order<- paste0("",result.ab.trop$Order)
list.result.ab.trop <- split(result.ab.trop,Order)
list2env(list.result.ab.trop,globalenv())

# express as a percentage of primary
Coleoptera$perc <- ((Coleoptera$y/Coleoptera$y[1]) * 100) - 100
Diptera$perc <- ((Diptera$y/Diptera$y[1]) * 100) - 100
Hemiptera$perc <- ((Hemiptera$y/Hemiptera$y[1]) * 100) - 100
Hymenoptera$perc <- ((Hymenoptera$y/Hymenoptera$y[1]) * 100) - 100
Lepidoptera$perc <- ((Lepidoptera$y/Lepidoptera$y[1]) * 100) - 100
Orthoptera$perc <- ((Orthoptera$y/Orthoptera$y[1]) * 100) - 100

# put it back together
result.ab.trop <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Orthoptera)

# add in SCA vals
result.ab.trop$SCA <- rep(c(0,0,0,0,1, 1, 1, 1),times=6)

# add zone factor
result.ab.trop$Zone <- as.factor("Tropical")

## Species Richness, Tropical ##

# predict the results
result.sr.trop <- PredictGLMER(model = MaxAnomalyModelRich_trop$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.sr.trop <- exp(result.sr.trop)

# add in the LU info
result.sr.trop$LUI <- data_tab$LUI

# add in the Order info
result.sr.trop$Order <- data_tab$Order

# break into Orders
Order<- paste0("",result.sr.trop$Order)
list.result.sr.trop <- split(result.sr.trop,Order)
list2env(list.result.sr.trop,globalenv())

# express as a percentage of primary
Coleoptera$perc <- ((Coleoptera$y/Coleoptera$y[1]) * 100) - 100
Diptera$perc <- ((Diptera$y/Diptera$y[1]) * 100) - 100
Hemiptera$perc <- ((Hemiptera$y/Hemiptera$y[1]) * 100) - 100
Hymenoptera$perc <- ((Hymenoptera$y/Hymenoptera$y[1]) * 100) - 100
Lepidoptera$perc <- ((Lepidoptera$y/Lepidoptera$y[1]) * 100) - 100
Orthoptera$perc <- ((Orthoptera$y/Orthoptera$y[1]) * 100) - 100

# put it back together
result.sr.trop <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Orthoptera)

# add in SCA vals
result.sr.trop$SCA <- rep(c(0,0,0,0,1, 1, 1, 1),times=6)

# add zone factor
result.sr.trop$Zone <- as.factor("Tropical")

## Max Anomaly , NonTropical ##
#create matrix for predictions 
# SCA = 1
# abun and richness = 0

# what is the rescaled value of SCA of 1
BackTransformCentreredPredictor(transformedX = 0.31, originalX = nontrop$StdTmaxAnomaly) # 0.31 gives about 1 

# what is the rescaled value of SCA of 0
BackTransformCentreredPredictor(transformedX = -1.02, originalX = nontrop$StdTmaxAnomaly) # -1.02 gives about 0

# reference is primary with 0 climate change so have 0 for that row

data_tab <- expand.grid(LUI = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High","Primary vegetation","Secondary vegetation", "Agriculture_Low", "Agriculture_High"),
                        Order = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera"),
                        LogAbund = 0,
                        Species_richness = 0)

# add column with SCA, values repeating 6 times
StdTmaxAnomalyRS = rep(c(-1.02,-1.02,-1.02,-1.02,0.31,0.31,0.31,0.31),times=6)

# add SCA to the data_tab
data_tab<-cbind(data_tab,StdTmaxAnomalyRS)

# factor the LU info
data_tab$LUI <- factor(data_tab$LUI, levels = levels(predictsSites$LUI))

## Abundance, NonTropical ##

# predict results
result.ab.nontrop <- PredictGLMER(model = MaxAnomalyModelAbund_nontrop$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.ab.nontrop <- exp(result.ab.nontrop)-0.01

# add in the LU info
result.ab.nontrop$LUI <- data_tab$LUI

# add in the Order info
result.ab.nontrop$Order <- data_tab$Order

# break into Orders
Order<- paste0("",result.ab.nontrop$Order)
list.result.ab.nontrop <- split(result.ab.nontrop,Order)
list2env(list.result.ab.nontrop,globalenv())

# express as a percentage of primary
Coleoptera$perc <- ((Coleoptera$y/Coleoptera$y[1]) * 100) - 100
Diptera$perc <- ((Diptera$y/Diptera$y[1]) * 100) - 100
Hemiptera$perc <- ((Hemiptera$y/Hemiptera$y[1]) * 100) - 100
Hymenoptera$perc <- ((Hymenoptera$y/Hymenoptera$y[1]) * 100) - 100
Lepidoptera$perc <- ((Lepidoptera$y/Lepidoptera$y[1]) * 100) - 100
Orthoptera$perc <- ((Orthoptera$y/Orthoptera$y[1]) * 100) - 10

# put it back together
result.ab.nontrop <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Orthoptera)

# add in SCA vals
result.ab.nontrop$SCA <- rep(c(0,0,0,0,1, 1, 1, 1),times=6)

# add zone factor
result.ab.nontrop$Zone <- as.factor("NonTropical")


## Species Richness, NonTropical ##

# predict the results
result.sr.nontrop <- PredictGLMER(model = MaxAnomalyModelRich_nontrop$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.sr.nontrop <- exp(result.sr.nontrop)

# add in the LU info
result.sr.nontrop$LUI <- data_tab$LUI

# add in the Order info
result.sr.nontrop$Order <- data_tab$Order

# break into Orders
Order<- paste0("",result.sr.nontrop$Order)
list.result.sr.nontrop <- split(result.sr.nontrop,Order)
list2env(list.result.sr.nontrop,globalenv())

# express as a percentage of primary
Coleoptera$perc <- ((Coleoptera$y/Coleoptera$y[1]) * 100) - 100
Diptera$perc <- ((Diptera$y/Diptera$y[1]) * 100) - 100
Hemiptera$perc <- ((Hemiptera$y/Hemiptera$y[1]) * 100) - 100
Hymenoptera$perc <- ((Hymenoptera$y/Hymenoptera$y[1]) * 100) - 100
Lepidoptera$perc <- ((Lepidoptera$y/Lepidoptera$y[1]) * 100) - 100
Orthoptera$perc <- ((Orthoptera$y/Orthoptera$y[1]) * 100) - 100

# put it back together
result.sr.nontrop <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Orthoptera)

# add in SCA vals
result.sr.nontrop$SCA <- rep(c(0,0,0,0,1, 1, 1, 1),times=6)

# add zone factor
result.sr.nontrop$Zone <- as.factor("NonTropical")

# combine results into a table for saving
all_res <- rbind(result.ab.nontrop, result.ab.trop, result.sr.nontrop, result.sr.trop)

all_res$measure <- c(rep("ab", 96), rep("sr", 96))

# save as png
percentage_change_LUI_CC <- all_res %>% gt()
gtsave(percentage_change_LUI_CC,"C:/Users/Kyra/Documents/GLITRS/Code/7_Predictions/MaxAnom_PercentageChange_LUI_CC_Tropical.png")


# save table as csv
write.csv(all_res, file = paste0(outDir, "/MaxAnom_PercentageChange_LUI_CC_Tropical.csv"))
