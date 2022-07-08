##%######################################################%##
#                                                          #
####        Predicting values from model outputs        ####
#                                                          #
##%######################################################%##

# Here, values for certain fixed effect combinations are predicted and 
# expressed as percentage change for use in the text of the paper. 


# directories
inDir <- "C:/Users/Kyra/Documents/GLITRS/Code/5_RunLUIClimateModels/"
outDir <- "C:/Users/Kyra/Documents/GLITRS/Code/7_Predictions/"

if(!dir.exists(outDir)) dir.create(outDir)


# sink(paste0(outDir,"log.txt"))
# 
# t.start <- Sys.time()
# 
# print(t.start)

# load libraries
library(StatisticalModels)
library(predictsFunctions)
library(webshot)
library(gt)
source('C:/Users/Kyra/Documents/GLITRS/Data/0_Functions.R')


# read in the predicts data
predictsSites <- readRDS(paste0(inDir,"PREDICTSSitesClimate_Data.rds"))


#### Hyp 1: land use effect only ####

# see LUI_plots.R

#### Hyp 2: Land use and climate anomaly interaction ####

# looking at STA of 1

# load in models
load(file = paste0(inDir, "MeanAnomalyModelAbund.rdata"))
load(file = paste0(inDir, "MeanAnomalyModelRich.rdata"))
load(file = paste0(inDir, "MaxAnomalyModelAbund.rdata"))
load(file = paste0(inDir, "MaxAnomalyModelRich.rdata"))

## Mean Anomaly ##
#create matrix for predictions 
# STA = 1
# abun and richness = 0

# what is the rescaled value of STA of 1
BackTransformCentreredPredictor(transformedX = 0.97, originalX = predictsSites$StdTmeanAnomaly) # 0.97 gives about 1 

# what is the rescaled value of STA of 0
BackTransformCentreredPredictor(transformedX = -1.32, originalX = predictsSites$StdTmeanAnomaly) # -1.32 gives about 0 

# reference is primary with 0 climate change so have 0 for that row

data_tab <- expand.grid(LUI = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High","Primary vegetation","Secondary vegetation", "Agriculture_Low", "Agriculture_High"),
                       Order = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera"),
                       LogAbund = 0,
                       Species_richness = 0)

# add column with STA, values repeating 6 times
StdTmeanAnomalyRS = rep(c(-1.32,-1.32,-1.32,-1.32,0.97,0.97,0.97,0.97),times=6)

# add STA to the data_tab
data_tab<-cbind(data_tab,StdTmeanAnomalyRS)

# factor the LUI info
data_tab$LUI <- factor(data_tab$LUI, levels = levels(predictsSites$LUI))

### Abundance ###

# predict results
result.ab <- PredictGLMER(model = MeanAnomalyModelAbund$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.ab <- exp(result.ab)-0.01

# add in the LU info
result.ab$LUI <- data_tab$LUI

# add in the Order info
result.ab$Order <- data_tab$Order

# break into Orders
Order<- paste0("",result.ab$Order)
list.result.ab <- split(result.ab,Order)
list2env(list.result.ab,globalenv())

# express as a percentage of primary
Coleoptera$perc <- ((Coleoptera$y/Coleoptera$y[1]) * 100) - 100
Diptera$perc <- ((Diptera$y/Diptera$y[1]) * 100) - 100
Hemiptera$perc <- ((Hemiptera$y/Hemiptera$y[1]) * 100) - 100
Hymenoptera$perc <- ((Hymenoptera$y/Hymenoptera$y[1]) * 100) - 100
Lepidoptera$perc <- ((Lepidoptera$y/Lepidoptera$y[1]) * 100) - 100
Orthoptera$perc <- ((Orthoptera$y/Orthoptera$y[1]) * 100) - 100

# put it back together
result.ab <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Orthoptera)

# add in STA vals
result.ab$STA <- rep(c(0,0,0,0,1, 1, 1, 1),times=6)

## Species Richness ##

# predict the results
result.sr <- PredictGLMER(model = MeanAnomalyModelRich$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.sr <- exp(result.sr)

# add in the LU info
result.sr$LUI <- data_tab$LUI

# add in the Order info
result.sr$Order <- data_tab$Order

# break into Orders
Order<- paste0("",result.sr$Order)
list.result.sr <- split(result.sr,Order)
list2env(list.result.sr,globalenv())

# express as a percentage of primary
Coleoptera$perc <- ((Coleoptera$y/Coleoptera$y[1]) * 100) - 100
Diptera$perc <- ((Diptera$y/Diptera$y[1]) * 100) - 100
Hemiptera$perc <- ((Hemiptera$y/Hemiptera$y[1]) * 100) - 100
Hymenoptera$perc <- ((Hymenoptera$y/Hymenoptera$y[1]) * 100) - 100
Lepidoptera$perc <- ((Lepidoptera$y/Lepidoptera$y[1]) * 100) - 100
Orthoptera$perc <- ((Orthoptera$y/Orthoptera$y[1]) * 100) - 100

# put it back together
result.sr <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Orthoptera)

# add in STA vals
result.sr$STA <- rep(c(0,0,0,0,1, 1, 1, 1),times=6)

# combine results into a table for saving
all_res <- rbind(result.ab, result.sr)

all_res$measure <- c(rep("ab", 48), rep("sr", 48))

# save as png
percentage_change_LUI_CC <- all_res %>% gt()
gtsave(percentage_change_LUI_CC,"C:/Users/Kyra/Documents/GLITRS/Code/7_Predictions/MeanAnom_PercentageChange_LUI_CC.png")

# save table
write.csv(all_res, file = paste0(outDir, "/MeanAnom_PercentageChange_LU_CC.csv"))

## Max Anomaly ##
#create matrix for predictions 
# STA = 1
# abun and richness = 0

# what is the rescaled value of STA of 1
BackTransformCentreredPredictor(transformedX = -0.14, originalX = predictsSites$StdTmaxAnomaly) # -0.14 gives about 1 

# what is the rescaled value of STA of 0
BackTransformCentreredPredictor(transformedX = -0.89, originalX = predictsSites$StdTmaxAnomaly) # -0.89 gives about 0

# reference is primary with 0 climate change so have 0 for that row

data_tab <- expand.grid(LUI = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High","Primary vegetation","Secondary vegetation", "Agriculture_Low", "Agriculture_High"),
                        Order = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera"),
                        LogAbund = 0,
                        Species_richness = 0)

# add column with STA, values repeating 6 times
StdTmaxAnomalyRS = rep(c(-0.89,-0.89,-0.89,-0.89,-0.14,-0.14,-0.14,-0.14),times=6)

# add STA to the data_tab
data_tab<-cbind(data_tab,StdTmaxAnomalyRS)

# factor the LU info
data_tab$LUI <- factor(data_tab$LUI, levels = levels(predictsSites$LUI))

## Species Richness ##

# predict the results
result.sr <- PredictGLMER(model = MaxAnomalyModelRich$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.sr <- exp(result.sr)

# add in the LU info
result.sr$LUI <- data_tab$LUI

# add in the Order info
result.sr$Order <- data_tab$Order

# break into Orders
Order<- paste0("",result.sr$Order)
list.result.sr <- split(result.sr,Order)
list2env(list.result.sr,globalenv())

# express as a percentage of primary
Coleoptera$perc <- ((Coleoptera$y/Coleoptera$y[1]) * 100) - 100
Diptera$perc <- ((Diptera$y/Diptera$y[1]) * 100) - 100
Hemiptera$perc <- ((Hemiptera$y/Hemiptera$y[1]) * 100) - 100
Hymenoptera$perc <- ((Hymenoptera$y/Hymenoptera$y[1]) * 100) - 100
Lepidoptera$perc <- ((Lepidoptera$y/Lepidoptera$y[1]) * 100) - 100
Orthoptera$perc <- ((Orthoptera$y/Orthoptera$y[1]) * 100) - 100

# put it back together
result.sr <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Orthoptera)

# add in STA vals
result.sr$STA <- rep(c(0,0,0,0,1, 1, 1, 1),times=6)

### Abundance ###

# predict results
result.ab <- PredictGLMER(model = MaxAnomalyModelAbund$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.ab <- exp(result.ab)-0.01

# add in the LU info
result.ab$LUI <- data_tab$LUI

# add in the Order info
result.ab$Order <- data_tab$Order

# break into Orders
Order<- paste0("",result.sr$Order)
list.result.ab <- split(result.ab,Order)
list2env(list.result.ab,globalenv())

# express as a percentage of primary
Coleoptera$perc <- ((Coleoptera$y/Coleoptera$y[1]) * 100) - 100
Diptera$perc <- ((Diptera$y/Diptera$y[1]) * 100) - 100
Hemiptera$perc <- ((Hemiptera$y/Hemiptera$y[1]) * 100) - 100
Hymenoptera$perc <- ((Hymenoptera$y/Hymenoptera$y[1]) * 100) - 100
Lepidoptera$perc <- ((Lepidoptera$y/Lepidoptera$y[1]) * 100) - 100
Orthoptera$perc <- ((Orthoptera$y/Orthoptera$y[1]) * 100) - 100

# put it back together
result.ab <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Orthoptera)

# add in STA vals
result.ab$STA <- rep(c(0,0,0,0,1, 1, 1, 1),times=6)

# combine results into a table for saving
all_res <- rbind(result.ab, result.sr)

all_res$measure <- c(rep("ab", 48), rep("sr", 48))

# save as png
percentage_change_LUI_CC <- all_res %>% gt()
gtsave(percentage_change_LUI_CC,"C:/Users/Kyra/Documents/GLITRS/Code/7_Predictions/MaxAnom_PercentageChange_LUI_CC.png")


# save table
write.csv(all_res, file = paste0(outDir, "/MaxAnom_PercentageChange_LUI_CC.csv"))


# t.end <- Sys.time()
# 
# print(round(t.end - t.start,0))
# 
# sink()
