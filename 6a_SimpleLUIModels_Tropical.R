## Get set up ##

inDir <- "C:/Users/Kyra/Documents/GLITRS/Code/1_CheckPrepareData/"
outDir <- "C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/"
predsDir <- "C:/Users/Kyra/Documents/GLITRS/Code/7_Predictions/"
if(!dir.exists(outDir)) dir.create(outDir)

# sink(paste0(outDir,"log_SimpleLUIModels_Tropical.txt"))

# t.start <- Sys.time()
# 
# print(t.start)

# load libraries
packages_model <- c("StatisticalModels", "predictsFunctions", "ggplot2", "cowplot", "sjPlot","dplyr")
suppressWarnings(suppressMessages(lapply(packages_model, require, character.only = TRUE)))

# source in additional functions
source("C:/Users/Kyra/Documents/GLITRS/Data/0_Functions.R")

packages_plot <- c("patchwork", "dplyr", "yarg", "lme4", "gt", "broom.mixed", "MASS")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# read in the PREDICTS Site data
sites <- readRDS(file = paste0(inDir,"PREDICTSSiteData.rds"))

# abundance data subset
table(sites[!is.na(sites$LogAbund), 'Realm'])
# NonTropical    Tropical 
#        6005        3014 

# species richness data subset
table(sites[!is.na(sites$Species_richness), 'Realm'])
# NonTropical    Tropical 
#        6288        3167 

# split by land use classes
table(sites$LUI, sites$Realm)
#                       NonTropical Tropical
# Agriculture_High            1820      822
# Agriculture_Low             1455      576
# Primary vegetation          1412      810
# Secondary vegetation        1601      959

# create separate data sets for Tropical and NonTropical sites
trop <- sites[sites$Realm == "Tropical", ]
nontrop <- sites[sites$Realm == "NonTropical", ]

## Species richness models - Tropical ##

# remove NAs in the specified columns
model_data_sr_trop <- na.omit(trop[,c('Species_richness','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order')])

# order data
model_data_sr_trop$LUI <- factor(model_data_sr_trop$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_sr_trop$Order <- factor(model_data_sr_trop$Order, levels = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera"))

# relevel
model_data_sr_trop$LUI <- relevel(model_data_sr_trop$LUI, ref = "Primary vegetation")

# summaries
length(unique(model_data_sr_trop$SS)) # 96
length(unique(model_data_sr_trop$SSBS)) # 1719

# look at the spread of land use/use intensity categories
print(table(model_data_sr_trop$LUI)) 

# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#                810                  959                  576                  822

# Run species richness models using GLMER function from StatisticalModels

# effect of land use (the null (intercept-only) model)
sm0_trop <-GLMER (modelData = model_data_sr_trop,responseVar = "Species_richness",fitFamily = "poisson",
                  fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# effect of LUI (LUI)
sm3_trop <- GLMER(modelData = model_data_sr_trop,responseVar = "Species_richness",fitFamily = "poisson",
                  fixedStruct = "LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# effect of order only
sm0.2_trop <- GLMER(modelData = model_data_sr_trop,responseVar = "Species_richness",fitFamily = "poisson",
                    fixedStruct = "Order",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effects of order and LUI
sm3.2_trop <- GLMER(modelData = model_data_sr_trop,responseVar = "Species_richness",fitFamily = "poisson",
                    fixedStruct = "Order+LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effect of order and LUI
sm3.3_trop <- GLMER(modelData = model_data_sr_trop,responseVar = "Species_richness",fitFamily = "poisson",
                    fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

## Species richness models - NonTropical ##

# remove NAs in the specified columns
model_data_sr_nontrop <- na.omit(nontrop[,c('Species_richness','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order')])

# order data
# exclude Neuroptera and Thysanoptera
model_data_sr_nontrop$LUI <- factor(model_data_sr_nontrop$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_sr_nontrop$Order <- factor(model_data_sr_nontrop$Order, levels = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera"))

# relevel
model_data_sr_nontrop$LUI <- relevel(model_data_sr_nontrop$LUI, ref = "Primary vegetation")

# summaries
length(unique(model_data_sr_nontrop$SS)) # 162
length(unique(model_data_sr_nontrop$SSBS)) # 4353

# look at the spread of land use/use intensity categories
print(table(model_data_sr_nontrop$LUI))
# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#               1412                 1601                 1455                 1820 

# effect of land use (the null (intercept-only) model)
sm0_nontrop <-GLMER (modelData = model_data_sr_nontrop,responseVar = "Species_richness",fitFamily = "poisson",
                     fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# effect of LUI (LUI)
sm3_nontrop <- GLMER(modelData = model_data_sr_nontrop,responseVar = "Species_richness",fitFamily = "poisson",
                     fixedStruct = "LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# effect of order only
sm0.2_nontrop <- GLMER(modelData = model_data_sr_nontrop,responseVar = "Species_richness",fitFamily = "poisson",
                       fixedStruct = "Order",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effects of order and LUI
sm3.2_nontrop <- GLMER(modelData = model_data_sr_nontrop,responseVar = "Species_richness",fitFamily = "poisson",
                       fixedStruct = "Order+LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effect of order and LUI
sm3.3_nontrop <- GLMER(modelData = model_data_sr_nontrop,responseVar = "Species_richness",fitFamily = "poisson",
                       fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# take a look at the AICs
AIC_sm_realm<-print(AIC(sm0_trop$model,sm3_trop$model,sm0.2_trop$model,sm3.2_trop$model,sm3.3_trop$model,
                        sm0_nontrop$model,sm3_nontrop$model,sm0.2_nontrop$model,sm3.2_nontrop$model,sm3.3_nontrop$model))

# Warning message:
#   In AIC.default(sm0_trop$model, sm3_trop$model, sm0.2_trop$model,  :
#                    models are not all fitted to the same number of observations

write.csv(AIC_sm_realm,"C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/AIC_sm_realm.csv", row.names = TRUE)

## Abundance models - Tropical ##

# remove NAs in the specified columns
model_data_ab_trop <- na.omit(trop[,c('LogAbund','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order')])

# order data
model_data_ab_trop$LUI <- factor(model_data_ab_trop$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_ab_trop$Order <- factor(model_data_ab_trop$Order, levels = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera"))

# relevel
model_data_ab_trop$LUI <- relevel(model_data_ab_trop$LUI, ref = "Primary vegetation")

# summaries
length(unique(model_data_ab_trop$SS)) # 85
length(unique(model_data_ab_trop$SSBS)) #1566

# look at the spread of land use/use intensity categories
print(table(model_data_ab_trop$LUI))
# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#                749                  908                  569                  788

# Run abundance models using 'GLMER' function from StatisticalModels

am0_trop <- GLMER(modelData = model_data_ab_trop,responseVar = "LogAbund",fitFamily = "gaussian",
                  fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3_trop <- GLMER(modelData = model_data_ab_trop,responseVar = "LogAbund",fitFamily = "gaussian",
                  fixedStruct = "LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am0.2_trop <- GLMER(modelData = model_data_ab_trop,responseVar = "LogAbund",fitFamily = "gaussian",
                    fixedStruct = "Order",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3.2_trop <- GLMER(modelData = model_data_ab_trop,responseVar = "LogAbund",fitFamily = "gaussian",
                    fixedStruct = "Order+LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3.3_trop <- GLMER(modelData = model_data_ab_trop,responseVar = "LogAbund",fitFamily = "gaussian",
                    fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

## Abundance models - NonTropical ##

# remove NAs in the specified columns
model_data_ab_nontrop <- na.omit(nontrop[,c('LogAbund','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order')])

# order data
model_data_ab_nontrop$LUI <- factor(model_data_ab_nontrop$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_ab_nontrop$Order <- factor(model_data_ab_nontrop$Order, levels = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera"))

# relevel
model_data_ab_nontrop$LUI <- relevel(model_data_ab_nontrop$LUI, ref = "Primary vegetation")

# summaries
length(unique(model_data_ab_nontrop$SS)) # 153
length(unique(model_data_ab_nontrop$SSBS)) # 4170

# look at the spread of land use/use intensity categories
print(table(model_data_ab_nontrop$LUI))
# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#               1367                 1507                 1439                 1692

# Run abundance models using 'GLMER' function from StatisticalModels

am0_nontrop <- GLMER(modelData = model_data_ab_nontrop,responseVar = "LogAbund",fitFamily = "gaussian",
                     fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3_nontrop <- GLMER(modelData = model_data_ab_nontrop,responseVar = "LogAbund",fitFamily = "gaussian",
                     fixedStruct = "LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am0.2_nontrop <- GLMER(modelData = model_data_ab_nontrop,responseVar = "LogAbund",fitFamily = "gaussian",
                       fixedStruct = "Order",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3.2_nontrop <- GLMER(modelData = model_data_ab_nontrop,responseVar = "LogAbund",fitFamily = "gaussian",
                       fixedStruct = "Order+LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3.3_nontrop <- GLMER(modelData = model_data_ab_nontrop,responseVar = "LogAbund",fitFamily = "gaussian",
                       fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)


# take a look at the AICs
AIC_am_realm<-print(AIC(am0_trop$model,am3_trop$model,am0.2_trop$model,am3.2_trop$model,am3.3_trop$model,
                        am0_nontrop$model,am3_nontrop$model,am0.2_nontrop$model,am3.2_nontrop$model,am3.3_nontrop$model))

# Warning message:
#   In AIC.default(am0_trop$model, am3_trop$model, am0.2_trop$model,  :
#                    models are not all fitted to the same number of observations

write.csv(AIC_am_realm,"C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/AIC_am_realm.csv", row.names = TRUE)


saveRDS(object = sm0_trop ,file = paste0(outDir,"sm0_trop.rds"))
saveRDS(object = sm3_trop ,file = paste0(outDir,"sm3_trop.rds"))
saveRDS(object = sm0.2_trop ,file = paste0(outDir,"sm0.2_trop.rds"))
saveRDS(object = sm3.2_trop ,file = paste0(outDir,"sm3.2_trop.rds"))
saveRDS(object = sm3.3_trop ,file = paste0(outDir,"sm3.3_trop.rds"))
saveRDS(object = sm0_nontrop ,file = paste0(outDir,"sm0_nontrop.rds"))
saveRDS(object = sm3_nontrop ,file = paste0(outDir,"sm3_nontrop.rds"))
saveRDS(object = sm0.2_nontrop ,file = paste0(outDir,"sm0.2_nontrop.rds"))
saveRDS(object = sm3.2_nontrop ,file = paste0(outDir,"sm3.2_nontrop.rds"))
saveRDS(object = sm3.3_nontrop ,file = paste0(outDir,"sm3.3_nontrop.rds"))
saveRDS(object = am0_trop ,file = paste0(outDir,"am0_trop.rds"))
saveRDS(object = am3_trop ,file = paste0(outDir,"am3_trop.rds"))
saveRDS(object = am0.2_trop ,file = paste0(outDir,"am0.2_trop.rds"))
saveRDS(object = am3.2_trop ,file = paste0(outDir,"am3.2_trop.rds"))
saveRDS(object = am3.3_trop ,file = paste0(outDir,"am3.3_trop.rds"))
saveRDS(object = am0_nontrop ,file = paste0(outDir,"am0_nontrop.rds"))
saveRDS(object = am3_nontrop ,file = paste0(outDir,"am3_nontrop.rds"))
saveRDS(object = am0.2_nontrop ,file = paste0(outDir,"am0.2_nontrop.rds"))
saveRDS(object = am3.2_nontrop ,file = paste0(outDir,"am3.2_nontrop.rds"))
saveRDS(object = am3.3_nontrop ,file = paste0(outDir,"am3.3_nontrop.rds"))
saveRDS(object = model_data_sr_trop ,file = paste0(outDir,"model_data_sr_trop.rds"))
saveRDS(object = model_data_sr_nontrop ,file = paste0(outDir,"model_data_sr_nontrop.rds"))
saveRDS(object = model_data_ab_trop ,file = paste0(outDir,"model_data_ab_trop.rds"))
saveRDS(object = model_data_ab_nontrop ,file = paste0(outDir,"model_data_ab_nontrop.rds"))

## Plot Figures ##
# read in model data
# sm0_trop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/sm0_trop.rds"))
# sm3_trop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/sm3_trop.rds"))
# sm0.2_trop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/sm0.2_trop.rds"))
# sm3.2_trop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/sm3.2_trop.rds"))
sm3.3_trop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/sm3.3_trop.rds"))
# am0_trop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/am0_trop.rds"))
# am3_trop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/am3_trop.rds"))
# am0.2_trop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/am0.2_trop.rds"))
# am3.2_trop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/am3.2_trop.rds"))
am3.3_trop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/am3.3_trop.rds"))
# sm0_nontrop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/sm0_nontrop.rds"))
# sm3_nontrop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/sm3_nontrop.rds"))
# sm0.2_nontrop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/sm0.2_nontrop.rds"))
# sm3.2_nontrop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/sm3.2_nontrop.rds"))
sm3.3_nontrop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/sm3.3_nontrop.rds"))
# am0_nontrop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/am0_nontrop.rds"))
# am3_nontrop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/am3_nontrop.rds"))
# am0.2_nontrop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/am0.2_nontrop.rds"))
# am3.2_nontrop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/am3.2_nontrop.rds"))
am3.3_nontrop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/am3.3_nontrop.rds"))
model_data_sr_trop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/model_data_sr_trop.rds"))
model_data_ab_trop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/model_data_ab_trop.rds"))
model_data_sr_nontrop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/model_data_sr_nontrop.rds"))
model_data_ab_nontrop <- readRDS(file = paste0("C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/model_data_ab_nontrop.rds"))

# select model_data for only three orders (Coleoptera, Hymenoptera, and Lepidoptera)
model_data_sr_trop <- filter(model_data_sr_trop, Order %in% c('Coleoptera', 'Hymenoptera', 'Lepidoptera'))
model_data_ab_trop <- filter(model_data_ab_trop, Order %in% c('Coleoptera', 'Hymenoptera', 'Lepidoptera'))
model_data_sr_nontrop <- filter(model_data_sr_nontrop, Order %in% c('Coleoptera', 'Hymenoptera', 'Lepidoptera'))
model_data_ab_nontrop <- filter(model_data_sr_nontrop, Order %in% c('Coleoptera', 'Hymenoptera', 'Lepidoptera'))

# table of AICs
selection_table_trop <- data.frame("Zone" = c(rep("Tropical", 10)),
                                   "Response" = c(rep("Species richness", 5),
                                                  rep("Total abundance", 5)),
                                   "Model" = c("Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS)",
                                               "Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                               "Species_richness ~ Order + (1|SS) + (1|SSB) + (1|SSBS)",
                                               "Species_richness ~ Order + LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                               "Species_richness ~ Order * LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                               "Total_abundance ~ 1 + (1|SS) + (1|SSB)",
                                               "Total_abundance ~ LUI + (1|SS) + (1|SSB)",
                                               "Total_abundance ~ Order + (1|SS) + (1|SSB)",
                                               "Total_abundance ~ Order + LUI + (1|SS) + (1|SSB)",
                                               "Total_abundance ~ Order * LUI + (1|SS) + (1|SSB)"),
                                   "AIC" = c(AIC(sm0_trop$model), AIC(sm3_trop$model), AIC(sm0.2_trop$model), AIC(sm3.2_trop$model), AIC(sm3.3_trop$model),
                                             AIC(am0_trop$model), AIC(am3_trop$model), AIC(am0.2_trop$model), AIC(am3.2_trop$model),AIC(am3.3_trop$model))) %>%
  group_by(Response) %>%
  arrange(AIC) %>%
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  gt()

selection_table_nontrop <- data.frame("Zone" = c(rep("NonTropical", 10)),
                                      "Response" = c(rep("Species richness", 5),
                                                     rep("Total abundance", 5)),
                                      "Model" = c("Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS)",
                                                  "Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                                  "Species_richness ~ Order + (1|SS) + (1|SSB) + (1|SSBS)",
                                                  "Species_richness ~ Order + LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                                  "Species_richness ~ Order * LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                                  "Total_abundance ~ 1 + (1|SS) + (1|SSB)",
                                                  "Total_abundance ~ LUI + (1|SS) + (1|SSB)",
                                                  "Total_abundance ~ Order + (1|SS) + (1|SSB)",
                                                  "Total_abundance ~ Order + LUI + (1|SS) + (1|SSB)",
                                                  "Total_abundance ~ Order * LUI + (1|SS) + (1|SSB)"),
                                      "AIC" = c(AIC(sm0_nontrop$model), AIC(sm3_nontrop$model), AIC(sm0.2_nontrop$model), AIC(sm3.2_nontrop$model), AIC(sm3.3_nontrop$model),  
                                                AIC(am0_nontrop$model), AIC(am3_nontrop$model), AIC(am0.2_nontrop$model), AIC(am3.2_nontrop$model),AIC(am3.3_nontrop$model))) %>%
  group_by(Response) %>%
  arrange(AIC) %>%
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  gt()

# save the tables
gtsave(selection_table_trop,"C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/output/LUIModels_Selection_trop.png")
gtsave(selection_table_nontrop,"C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/output/LUIModels_Selection_nontrop.png")

## Species Richness, NonTropical ##

richness_metric_nontrop <- predict_effects(iterations = 1000,
                                           model = sm3.3_nontrop$model,
                                           model_data = model_data_sr_nontrop,
                                           response_variable = "Species_richness",
                                           fixed_number = 2,
                                           fixed_column = c("Order", "LUI"),
                                           factor_number_1 = 3,
                                           factor_number_2 = 4,
                                           neg_binom = FALSE)
# predictions
result.sr.nontrop <- fin_conf
result.sr.nontrop <- dplyr::select(result.sr.nontrop,-c(Species_richness))
richness_metric_nontrop

model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_sr_nontrop$LUI), 3)[1:12]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# richness data
model_data(richness_metric_nontrop)

## Species Richness, Tropical ##

richness_metric_trop <- predict_effects(iterations = 1000,
                                        model = sm3.3_trop$model,
                                        model_data = model_data_sr_trop,
                                        response_variable = "Species_richness",
                                        fixed_number = 2,
                                        fixed_column = c("Order", "LUI"),
                                        factor_number_1 = 3,
                                        factor_number_2 = 4,
                                        neg_binom = FALSE)
result.sr.trop <- fin_conf
result.sr.trop <- dplyr::select(result.sr.trop,-c(Species_richness))
richness_metric_trop

model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_sr_trop$LUI), 3)[1:12]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# richness data
model_data(richness_metric_trop)


## Abundance, NonTropical ##

abundance_metric_nontrop <- predict_effects(iterations = 1000,
                                            model = am3.3_nontrop$model,
                                            model_data = model_data_ab_nontrop,
                                            response_variable = "LogAbund",
                                            fixed_number = 2,
                                            fixed_column = c("Order", "LUI"),
                                            factor_number_1 = 3,
                                            factor_number_2 = 4,
                                            neg_binom = FALSE)
result.ab.nontrop <- fin_conf
result.ab.nontrop <- dplyr::select(result.ab.nontrop,-c(LogAbund))
abundance_metric_nontrop

model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_ab_nontrop$LUI), 3)[1:12]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# abundance data
model_data(abundance_metric_nontrop)

## Abundance, Tropical ##

abundance_metric_trop <- predict_effects(iterations = 1000,
                                         model = am3.3_trop$model,
                                         model_data = model_data_ab_trop,
                                         response_variable = "LogAbund",
                                         fixed_number = 2,
                                         fixed_column = c("Order", "LUI"),
                                         factor_number_1 = 3,
                                         factor_number_2 = 4,
                                         neg_binom = FALSE)
result.ab.trop <- fin_conf
result.ab.trop <- dplyr::select(result.ab.trop,-c(LogAbund))
abundance_metric_trop

model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_ab_trop$LUI), 3)[1:12]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# abundance data
model_data(abundance_metric_trop)

# save predictions
# combine results into a table for saving
all_res <- rbind(result.ab.nontrop, result.ab.trop, result.sr.nontrop, result.sr.trop)
all_res$measure <- c(rep("ab",24 ), rep("sr", 24))
all_res$Zone <- c(rep("nontrop",12), rep("trop", 12),rep("nontrop",12), rep("trop",12))

# save as gt table
percentage_change_LUI_tropical <- all_res %>% gt()
# gtsave() doesn't like using the "predsDir" filepath, have to use the whole path
gtsave(percentage_change_LUI_tropical,"C:/Users/Kyra/Documents/GLITRS/Paper/Code/percentage_change_LUI_tropical.png")

# save table
write.csv(all_res, file = paste0(predsDir,"percentage_change_LUI_tropical.csv"))


# plot realms separately

(richness_metric_trop + xlab(NULL) + guides(colour = FALSE) + ggtitle("A") + scale_y_continuous("Species richness diff. (%)") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())) + 
  (abundance_metric_trop + xlab(NULL) + ggtitle("B") + scale_y_continuous("Total abundance diff. (%)") + theme(axis.text.x = element_text(size = 10,angle=45,margin=margin(t=20)), axis.ticks = element_blank(), legend.position = "bottom",legend.text = element_text(size = 10), legend.title = element_text(size = 11)) + guides(colour = guide_legend("Land-use intensity")))+ plot_layout(ncol = 1) 

ggsave("LUI_predictions_tropical.jpeg", device ="jpeg", path = outDir, width=25, height=15, units="cm", dpi = 350)

(richness_metric_nontrop + xlab(NULL) + guides(colour = FALSE) + ggtitle("A") + scale_y_continuous("Species richness diff. (%)") + theme(axis.text.x = element_blank(), axis.ticks = element_blank())) + 
  (abundance_metric_nontrop + xlab(NULL) + ggtitle("B") + scale_y_continuous("Total abundance diff. (%)") + theme(axis.text.x = element_text(size = 10,angle=45,margin=margin(t=20)), axis.ticks = element_blank(), legend.position = "bottom",legend.text = element_text(size = 10), legend.title = element_text(size = 11)) + guides(colour = guide_legend("Land-use intensity")))+ plot_layout(ncol = 1) 

ggsave("LUI_predictions_nontropical.jpeg", device ="jpeg", path = outDir, width=25, height=15, units="cm", dpi = 350)

# plot realms together

title_nontrop <- ggdraw() + 
  draw_label(
    "NonTropical",
    fontface = 'bold',
    size = 10)

title_trop <- ggdraw() + 
  draw_label(
    "Tropical",
    fontface = 'bold',
    size = 10)

richness_metric_nontrop <- richness_metric_nontrop + 
  xlab(NULL) + 
  ylab("Species richness diff. (%)")+
  ggtitle("A") + 
  scale_y_continuous(breaks = c(-100,-50, 0, 50, 100, 150, 200), limits = c(-100, 200)) + 
  theme(axis.text.x = element_blank(),
        axis.title.y = element_text(angle=90),
        legend.position = "none")

richness_metric_trop <- richness_metric_trop + 
  xlab(NULL) + 
  ylab(NULL)+
  guides(colour = "none") + 
  ggtitle("B") + 
  scale_y_continuous(breaks = c(-100,-50, 0, 50, 100, 150, 200), limits = c(-100, 200)) + 
  theme(axis.text.x = element_blank())

abundance_metric_nontrop <- abundance_metric_nontrop + 
  xlab(NULL) +
  ylab("Total abundance diff. (%)")+
  ggtitle("C") + 
  scale_y_continuous(breaks = c(-100,-50, 0, 50, 100, 150, 200), limits = c(-100, 200)) + 
  theme(axis.text.x = element_text(size = 10,
                                   angle=45,
                                   margin=margin(t=20)), 
        legend.position = "none") + 
  guides(colour = guide_legend("Land-use intensity"))

abundance_metric_trop <- abundance_metric_trop + 
  xlab(NULL) + 
  ylab(NULL)+
  ggtitle("D") + 
  scale_y_continuous(breaks = c(-100,-50, 0, 50, 100, 150, 200), limits = c(-100, 200)) +
  theme(axis.text.x = element_text(size = 10,
                                   angle=45,
                                   margin=margin(t=20)), 
        legend.position = "none")

legend <- get_legend(
  richness_metric_nontrop +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(),
          legend.text = element_text(size = 10),
          legend.title = element_blank())
)

LUI_predictions_realm<-cowplot::plot_grid(title_nontrop, title_trop, richness_metric_nontrop, richness_metric_trop,abundance_metric_nontrop,abundance_metric_trop,
                                          ncol=2,
                                          rel_heights = c(0.1,1,1.4))

LUI_predictions_realm <- cowplot::plot_grid(LUI_predictions_realm, legend, 
                                            ncol = 1,
                                            rel_heights = c(1,0.1))

ggsave("LUI_predictions_realms.jpeg", device ="jpeg", path = outDir, width=40, height=24, units="cm", dpi = 350)


# t.end <- Sys.time()
# 
# print(round(t.end - t.start,0))
# 
# sink()
