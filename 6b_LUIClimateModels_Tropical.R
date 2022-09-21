###Create Models for all insect orders in predicts for standardised climate anomaly and Land interactions ###

## get set up ##

# directories 
predictsDir<- "C:/Users/Kyra/Documents/GLITRS/Code/5_RunLUIClimateModels/"
outDir <- "C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/"
tabDir <- "C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Tables/"
if(!dir.exists(outDir)) dir.create(outDir)
if(!dir.exists(tabDir)) dir.create(tabDir)

# sink(paste0(outDir,"log_LUI_ClimateModels.txt"))
# 
# t.start <- Sys.time()
# 
# print(t.start)

# load libraries
library(devtools)
library(StatisticalModels)
library(predictsFunctions)

# source additional functions
source("C:/Users/Kyra/Documents/GLITRS/Data/0_Functions.R")

# read in the predicts data from the LUI Climate Models
predictsSites <- readRDS(paste0(predictsDir,"PREDICTSSitesClimate_Data.rds"))

# split into two data frames

trop <- predictsSites[predictsSites$Realm == "Tropical", ]
nontrop <- predictsSites[predictsSites$Realm == "NonTropical", ]

# take a look at possible correlations between variables
# Tropical
cor(trop$avg_temp, trop$TmeanAnomaly)

# 0.07376417

cor(trop$avg_temp, trop$StdTmeanAnomaly)

#0.2118244

cor(trop$TmeanAnomaly, trop$StdTmeanAnomaly)

# -0.1501622

# NonTropical
cor(nontrop$avg_temp, nontrop$TmeanAnomaly)

# -0.246491

cor(nontrop$avg_temp, nontrop$StdTmeanAnomaly)

#0.2118244

cor(nontrop$TmeanAnomaly, nontrop$StdTmeanAnomaly)

# 0.2639987

# save the datasets
saveRDS(object = trop,file = paste0(outDir,"trop.rds"))
saveRDS(object = nontrop,file = paste0(outDir,"nontrop.rds"))


## Model Selection - Tropical ##

# 1. Abundance, Mean Anomaly, Tropical

model_data <- trop[!is.na(trop$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MeanAnomalyModelAbund_trop <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                                fixedStruct = "UI2 * StdTmeanAnomalyRS * Order",
                                randomStruct = "(1|SS)+(1|SSB)",
                                saveVars = c("SSBS"))

# get summary
summary(MeanAnomalyModelAbund_trop$model)


# save the model output
save(MeanAnomalyModelAbund_trop, file = paste0(outDir, "MeanAnomalyModelAbund_trop.rdata"))



# 2. Richness, Mean Anomaly, Tropical

model_data <- trop[!is.na(trop$StdTmeanAnomalyRS), ]

MeanAnomalyModelRich_trop <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmeanAnomalyRS * Order",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                               saveVars = c("SSBS"))

summary(MeanAnomalyModelRich_trop$model)

save(MeanAnomalyModelRich_trop, file = paste0(outDir, "MeanAnomalyModelRich_trop.rdata"))

# 3. Abundance, Max Anomaly, Trop

model_data <- trop[!is.na(trop$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MaxAnomalyModelAbund_trop <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                               fixedStruct = "UI2 * StdTmaxAnomalyRS * Order",
                               randomStruct = "(1|SS)+(1|SSB)",
                               saveVars = c("SSBS"))

summary(MaxAnomalyModelAbund_trop$model)

save(MaxAnomalyModelAbund_trop, file = paste0(outDir, "MaxAnomalyModelAbund_trop.rdata"))

# 4. Richness, Max Anomaly, Tropical

model_data <- trop[!is.na(trop$StdTmeanAnomalyRS),]

MaxAnomalyModelRich_trop <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                              fixedStruct = "UI2 * StdTmaxAnomalyRS * Order",
                              randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                              saveVars = c("SSBS"))

# save model output
save(MaxAnomalyModelRich_trop, file = paste0(outDir, "MaxAnomalyModelRich_trop.rdata"))


## Model Selection - NonTropical ##

# 1. Abundance, Mean Anomaly, Nontropical

model_data <- nontrop[!is.na(nontrop$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MeanAnomalyModelAbund_nontrop <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                                    fixedStruct = "UI2 * StdTmeanAnomalyRS * Order",
                                    randomStruct = "(1|SS)+(1|SSB)",
                                    saveVars = c("SSBS"))

# get summary
summary(MeanAnomalyModelAbund_nontrop$model)


# save the model output
save(MeanAnomalyModelAbund_nontrop, file = paste0(outDir, "MeanAnomalyModelAbund_nontrop.rdata"))



# 2. Richness, Mean Anomaly, Nontropical

model_data <- nontrop[!is.na(nontrop$StdTmeanAnomalyRS), ]

MeanAnomalyModelRich_nontrop <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                                   fixedStruct = "UI2 * StdTmeanAnomalyRS * Order",
                                   randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                   saveVars = c("SSBS"))

summary(MeanAnomalyModelRich_nontrop$model)

save(MeanAnomalyModelRich_nontrop, file = paste0(outDir, "MeanAnomalyModelRich_nontrop.rdata"))

# 3. Abundance, Max Anomaly, Nontrop

model_data <- nontrop[!is.na(nontrops$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MaxAnomalyModelAbund_nontrop <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                                   fixedStruct = "UI2 * StdTmaxAnomalyRS * Order",
                                   randomStruct = "(1|SS)+(1|SSB)",
                                   saveVars = c("SSBS"))

summary(MaxAnomalyModelAbund_nontrop$model)

save(MaxAnomalyModelAbund_nontrop, file = paste0(outDir, "MaxAnomalyModelAbund_nontrop.rdata"))

# 4. Richness, Max Anomaly, Nontropical

model_data <- nontrop[!is.na(nontrop$StdTmeanAnomalyRS),]

MaxAnomalyModelRich_nontrop <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                                  fixedStruct = "UI2 * StdTmaxAnomalyRS * Order",
                                  randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                  saveVars = c("SSBS"))

# save model output
save(MaxAnomalyModelRich_nontrop, file = paste0(outDir, "MaxAnomalyModelRich_nontrop.rdata"))


#### save model output tables for use in supplementary information ####

# tropical

tab_model(MeanAnomalyModelAbund_trop$model, transform = NULL, file = paste0(tabDir,"AbunMeanAnom_output_trop_table.html"))
summary(MeanAnomalyModelAbund_trop$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelAbund_trop$model) # check the R2 values 

tab_model(MeanAnomalyModelRich_trop$model, transform = NULL, file = paste0(tabDir, "RichMeanAnom_output_trop_table.html"))
summary(MeanAnomalyModelRich_trop$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelRich_trop$model) # check the R2 values

tab_model(MaxAnomalyModelAbund_trop$model, transform = NULL, file = paste0(tabDir,"AbunMaxAnom_output_trop_table.html"))
summary(MaxAnomalyModelAbund_trop$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelAbund_trop$model) # check the R2 values 

tab_model(MaxAnomalyModelRich_trop$model, transform = NULL, file = paste0(tabDir,"RichMaxAnom_output_trop_table.html"))
summary(MaxAnomalyModelRich_trop$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelRich_trop$model) # check the R2 values 

# nontropical
 
tab_model(MeanAnomalyModelAbund_nontrop$model, transform = NULL, file = paste0(tabDir,"AbunMeanAnom_output_nontrop_table.html"))
summary(MeanAnomalyModelAbund_nontrop$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelAbund_nontrop$model) # check the R2 values 

tab_model(MeanAnomalyModelRich_nontrop$model, transform = NULL, file = paste0(tabDir,"RichMeanAnom_output_nontrop_table.html"))
summary(MeanAnomalyModelRich_nontrop$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelRich_nontrop$model) # check the R2 values

tab_model(MaxAnomalyModelAbund_nontrop$model, transform = NULL, file = paste0(tabDir,"AbunMaxAnom_output_nontrop_table.html"))
summary(MaxAnomalyModelAbund_nontrop$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelAbund_nontrop$model) # check the R2 values 

tab_model(MaxAnomalyModelRich_nontrop$model, transform = NULL, file = paste0(tabDir,"RichMaxAnom_output_nontrop_table.html"))
summary(MaxAnomalyModelRich_nontrop$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelRich_nontrop$model) # check the R2 values 

