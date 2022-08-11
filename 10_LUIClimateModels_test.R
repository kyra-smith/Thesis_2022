### Test models with and without Order as a fixed effect ###
## Get set up ##

# directories 
inDir <- "C:/Users/Kyra/Documents/GLITRS/Code/5_RunLUIClimateModels/"
outDir <- "C:/Users/Kyra/Documents/GLITRS/Code/10_Additional_Tests/"
if(!dir.exists(inDir)) dir.create(inDir)
if(!dir.exists(outDir)) dir.create(outDir)


sink(paste0(outDir,"log_LUI_ClimateModels_test.txt"))

t.start <- Sys.time()

print(t.start)

# load libraries
packages_model <- c("devtools","StatisticalModels", "predictsFunctions", "ggplot2", "cowplot", "sjPlot","dplyr")
suppressWarnings(suppressMessages(lapply(packages_model, require, character.only = TRUE)))

packages_plot <- c("patchwork", "dplyr", "yarg", "lme4", "gt", "broom.mixed", "MASS")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# source in additional functions
source("C:/Users/Kyra/Documents/GLITRS/Data/0_Functions.R")

# read in the predicts climate data
predictsSites <- readRDS(paste0(inDir,"PREDICTSSitesClimate_Data.rds"))

## Model Selection ##

# 1. Abundance, mean anomaly

model_data <- predictsSites[!is.na(predictsSites$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MeanAnomalyModelAbund <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                               fixedStruct = "Order * LUI * StdTmeanAnomalyRS",
                               randomStruct = "(1|SS)+(1|SSB)",
                               saveVars = c("SSBS"))

MeanAnomalyModelAbund_test <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                               fixedStruct = "LUI * StdTmeanAnomalyRS",
                               randomStruct = "(1|SS)+(1|SSB)",
                               saveVars = c("SSBS"))

# take a look at the AICs
AIC_MeanAbund<-print(AIC(MeanAnomalyModelAbund_test$model,MeanAnomalyModelAbund$model))

# 2. Richness, mean anomaly

model_data <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ]

MeanAnomalyModelRich <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                              fixedStruct = "Order * LUI * StdTmeanAnomalyRS",
                              randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                              saveVars = c("SSBS"))

MeanAnomalyModelRich_test <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                              fixedStruct = "LUI * StdTmeanAnomalyRS",
                              randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                              saveVars = c("SSBS"))

# take a look at the AICs
AIC_MeanRich<-print(AIC(MeanAnomalyModelRich_test$model,MeanAnomalyModelRich$model))

# 3. Abundance, max anomaly

model_data <- predictsSites[!is.na(predictsSites$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MaxAnomalyModelAbund <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                              fixedStruct = "Order * LUI * StdTmaxAnomalyRS",
                              randomStruct = "(1|SS)+(1|SSB)",
                              saveVars = c("SSBS"))

MaxAnomalyModelAbund_test <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                              fixedStruct = "LUI * StdTmaxAnomalyRS",
                              randomStruct = "(1|SS)+(1|SSB)",
                              saveVars = c("SSBS"))

# take a look at the AICs
AIC_MaxAbund<-print(AIC(MaxAnomalyModelAbund_test$model,MaxAnomalyModelAbund$model))

# 4. Richness, max anomaly

model_data <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS),]

MaxAnomalyModelRich <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                             fixedStruct = "Order * LUI * StdTmaxAnomalyRS",
                             randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                             saveVars = c("SSBS"))

MaxAnomalyModelRich_test <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                             fixedStruct = "LUI * StdTmaxAnomalyRS",
                             randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                             saveVars = c("SSBS"))

# take a look at the AICs
AIC_MaxRich<-print(AIC(MaxAnomalyModelRich_test$model, MaxAnomalyModelRich$model))

# put together and save .csv file
AICs_test <- rbind(AIC_MeanAbund,AIC_MaxAbund,AIC_MeanRich,AIC_MaxRich)

write.csv(AICs_test, file = paste0(outDir,"AICs_test.csv"))

#### save model output tables for use in supplementary information ####
# use function from sjPlot library to save neat versions of model output table
# conditional: the conditional R2 value, i.e. the variance explained by fixed and random effects 
# marginal: the marginal R2 value, i.e. the variance explained by the fixed effects

tab_model(MeanAnomalyModelAbund_test$model, transform = NULL, file = paste0(outDir,"Tables/AbunMeanAnom_test_output_table.html"))
summary(MeanAnomalyModelAbund_test$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelAbund_test$model) # check the R2 values 

tab_model(MeanAnomalyModelRich_test$model, transform = NULL, file = paste0(outDir,"Tables/RichMeanAnom_test_output_table.html"))
summary(MeanAnomalyModelRich_test$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelRich_test$model) # check the R2 values

tab_model(MaxAnomalyModelAbund_test$model, transform = NULL, file = paste0(outDir,"Tables/AbunMaxAnom_test_output_table.html"))
summary(MaxAnomalyModelAbund_test$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelAbund_test$model) # check the R2 values 

tab_model(MaxAnomalyModelRich_test$model, transform = NULL, file = paste0(outDir,"Tables/RichMaxAnom_test_output_table.html"))
summary(MaxAnomalyModelRich_test$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelRich_test$model) # check the R2 values 

# table of AICs
# species richness and abundance together

# load AIC.csv file
AICs <- read.csv("C:/Users/Kyra/Documents/GLITRS/Code/10_Additional_Tests/AICs_test.csv", header=TRUE, stringsAsFactors=FALSE)
# AIC values pulled from CSV file

# selection_table <- data.frame("Response" = c(rep("Species richness", 4),
#                                              rep("Total abundance", 4)),
#                               "Climate Anomaly" = c("Mean","Mean","Maximum","Maximum",
#                                                     "Mean","Mean","Maximum","Maximum"),
#                               "Model" = c("Species_richness ~ LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Species_richness ~ Order * LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Species_richness ~ LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Species_richness ~ Order * LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Total_abundance ~ LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB)",
#                                           "Total_abundance ~ Order * LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB)",
#                                           "Total_abundance ~ LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB)",
#                                           "Total_abundance ~ Order * LUI * StdTmeanAnomalyRS + (1|SS) + (1|SSB)"),
#                               "AIC" = c(AIC(MeanAnomalyModelRich_test$model), AIC(MeanAnomalyModelRich$model), AIC(MaxAnomalyModelRich_test$model), AIC(MaxAnomalyModelRich$model),  
#                                         AIC(MeanAnomalyModelAbund_test$model), AIC(MeanAnomalyModelAbund$model), AIC(MaxAnomalyModelAbund_test$model), AIC(MaxAnomalyModelAbund$model))) %>%
#   group_by(Response) %>%                              
#   mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
#   ungroup() %>%
#   gt()

# add model descriptions
AICs$Model <- c("Abundance ~ LUI * STA + (1|SS) + (1|SSB)",
                "Abundance ~ Order * LUI * STA + (1|SS) + (1|SSB)",
                "Abundance ~ LUI * SMTA + (1|SS) + (1|SSB)",
                "Abundance ~ Order * LUI * SMTA + (1|SS) + (1|SSB)",
                "Species richness ~ LUI * STA + (1|SS) + (1|SSB) + (1|SSBS)",
                "Species richness ~ Order * LUI * STA + (1|SS) + (1|SSB) + (1|SSBS)",
                "Species richness ~ LUI * SMTA + (1|SS) + (1|SSB) + (1|SSBS)",
                "Species richness ~ Order * LUI * SMTA + (1|SS) + (1|SSB) + (1|SSBS)")


# add climate anomaly
AICs$Anomaly <- c(rep("Mean",2),rep("Maximum",2),rep("Mean",2),rep("Maximum",2))

# add climate anomaly
AICs$Response <- c(rep("Abundance",4),rep("Species_richness",4))

AICs <- as.data.frame(AICs)

# drop columns 'X' and 'df' and re-order columns
AICs <- subset (AICs,select = c(Response,Anomaly,Model,AIC))

# make gt table of model selection
# AIC_select <- AICs %>% 
#   group_by(Response,Anomaly) %>%                              
#   mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
#   ungroup() %>%
#   select(Anomaly,Model,AIC,deltaAIC) %>%
#   gt(rowname_col = "Model") %>%
#   tab_row_group(
#     label = "Species Richness",
#     rows = starts_with("Species richness")
#   ) %>%
#   tab_row_group(
#     label = "Abundance",
#     rows = starts_with("Abundance")
#   ) %>%
#   tab_stubhead(label = "Models")

AIC_select <- AICs %>% 
  group_by(Response,Anomaly) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  select(Model,Anomaly,AIC,deltaAIC) %>%
  gt(rowname_col = "Model") %>%
  tab_row_group(
    label = "Species Richness",
    rows = starts_with("Species richness")
  ) %>%
  tab_row_group(
    label = "Abundance",
    rows = starts_with("Abundance")
  ) %>% 
  cols_align(
    align = "center",
    columns = c(Model, Anomaly, AIC, deltaAIC)
  )%>%
  tab_stubhead(label = "Models") %>%
  cols_label(
    Model = md ("Models"),
    Anomaly = md("Temperature Anomaly"),
    AIC = md("AIC"),
    deltaAIC = md("deltaAIC")
  ) 

# save

gtsave(AIC_select,"C:/Users/Kyra/Documents/GLITRS/Code/10_Additional_Tests/LUI_CC_ModelSelection.png")

t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()
