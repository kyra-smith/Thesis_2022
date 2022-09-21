## Get set up ##

inDir <- "C:/Users/Kyra/Documents/GLITRS/Code/1_CheckPrepareData/"
outDir <- "C:/Users/Kyra/Documents/GLITRS/Code/2_RunSimpleLUIModel/Output/"
predsDir <- "C:/Users/Kyra/Documents/GLITRS/Code/7_Predictions/"
# if(!dir.exists(outDir)) dir.create(outDir)

# sink(paste0(outDir,"log_SimpleLUIModels.txt"))
# 
# t.start <- Sys.time()
# 
# print(t.start)

# load libraries
packages_model <- c("StatisticalModels", "predictsFunctions", "ggplot2", "cowplot", "sjPlot","dplyr")
suppressWarnings(suppressMessages(lapply(packages_model, require, character.only = TRUE)))

packages_plot <- c("patchwork", "dplyr", "yarg", "lme4", "gt", "broom.mixed", "MASS","webshot")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# source in additional functions
source("C:/Users/Kyra/Documents/GLITRS/Data/0_Functions.R")

# read in the Site data
sites <- readRDS(file = paste0(inDir,"PREDICTSSiteData.rds"))

## Species Richness Model ##

# remove NAs in the specified columns
model_data_sr <- na.omit(sites[,c('Species_richness','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order')])

# order data
model_data_sr$LUI <- factor(model_data_sr$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_sr$Order <- factor(model_data_sr$Order, levels = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera"))

# relevel
model_data_sr$LUI <- relevel(model_data_sr$LUI, ref = "Primary vegetation")

#save model_data_sr
saveRDS(object = model_data_sr ,file = paste0(outDir,"model_data_sr.rds"))

# summaries
length(unique(model_data_sr$SS)) # 258
length(unique(model_data_sr$SSBS)) # 6072

# look at the spread of land use/use intensity categories
print(table(model_data_sr$LUI))

# Primary vegetation Secondary vegetation      Agriculture_Low     Agriculture_High 
#               2222                 2560                 2031                 2642 

# Run species richness models using GLMER function from StatisticalModels

# null (intercept-only) model
sm0 <-GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# effect of Land Use
sm1 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LandUse",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effects of land use and use intensity
sm2 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# effect of LUI
sm3 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# interactive effects of land use and use intensity
sm4 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# effect of order only
sm0.2 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effects of order and land use
sm1.2 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order+LandUse",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effects of order and land use and use intensity
sm2.2 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order+LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effects of order and LUI
sm3.2 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order+LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effects of order and interactive effects of land use and use intensity
sm4.2 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order+LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# interactive effect of order and land use
sm1.3 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order*LandUse",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# interactive effects of order and additive effects of land use and use intensity
sm2.3 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order*LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# additive effect of order and LUI
sm3.3 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# interactive effects of order, land use, and use intensity
sm4.3 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
               fixedStruct = "Order*LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

# take a look at the AICs
AIC_sm<-print(AIC(sm0$model,sm1$model,sm2$model,sm3$model,sm4$model,
                  sm0.2$model,sm1.2$model,sm2.2$model,sm3.2$model,sm4.2$model,
                  sm1.3$model,sm2.3$model,sm3.3$model,sm4.3$model))

write.csv(AIC_sm, file = paste0(outDir,"AIC_sm.csv"))

# Run richness models using 'glmr' function from lme4
# need these to run allFit()
# need to run allFit() to test the models

g_sm0 <- glmer(Species_richness~1+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm1 <- glmer(Species_richness~LandUse+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm2 <- glmer(Species_richness~LandUse+Use_intensity+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm3 <- glmer(Species_richness~LUI+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm4 <- glmer(Species_richness~LandUse*Use_intensity+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm0.2 <- glmer(Species_richness~Order+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm1.2 <- glmer(Species_richness~Order+LandUse+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm2.2 <- glmer(Species_richness~Order+LandUse+Use_intensity+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm3.2 <- glmer(Species_richness~Order+LUI+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm4.2 <- glmer(Species_richness~Order+LandUse*Use_intensity+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm1.3 <- glmer(Species_richness~Order*LandUse+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm2.3 <- glmer(Species_richness~Order*LandUse+Use_intensity+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm3.3 <- glmer(Species_richness~Order*LUI+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)
g_sm4.3 <- glmer(Species_richness~Order*LandUse*Use_intensity+(1|SS)+(1|SSB)+(1|SSBS), data = model_data_sr, family = poisson)

# take a look at the AICs
AIC_g_sm<-print(AIC(g_sm0,g_sm1,g_sm2,g_sm3,g_sm4,
                    g_sm0.2,g_sm1.2,g_sm2.2,g_sm3.2,g_sm4.2,
                    g_sm1.3,g_sm2.3,g_sm3.3,g_sm4.3))

write.csv(AIC_g_sm, file = paste0(outDir,"AIC_g_sm.csv"))

# Testing richness models with allFit()

g_sm0.2_all <-allFit(g_sm0.2)
g_sm1.2_all <-allFit(g_sm1.2)
g_sm2.2_all <-allFit(g_sm2.2)
g_sm3.2_all <-allFit(g_sm3.2)
g_sm4.2_all <-allFit(g_sm4.2)
g_sm1.3_all <-allFit(g_sm1.3)
g_sm2.3_all <-allFit(g_sm2.3)
g_sm3.3_all <-allFit(g_sm3.3)
g_sm4.3_all <-allFit(g_sm4.3)

## Abundance Models  ##

# remove NAs in the specified columns
model_data_ab <- na.omit(sites[,c('LogAbund','LandUse','Use_intensity','LUI','SS','SSB','SSBS','Order')])

# order data
model_data_ab$LUI <- factor(model_data_ab$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
model_data_ab$Order <- factor(model_data_ab$Order, levels = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera"))

# relevel
model_data_ab$LUI <- relevel(model_data_ab$LUI, ref = "Primary vegetation")

#save model_data_ab
saveRDS(object = model_data_ab ,file = paste0(outDir,"model_data_ab.rds"))

# summaries
length(unique(model_data_ab$SS))
length(unique(model_data_ab$SSBS)) 

# look at the spread of land use/use intensity categories
print(table(model_data_ab$LUI))

# Run abundance models using 'GLMER' function from StatisticalModels

am0 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am1 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LandUse",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am4 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am0.2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am1.2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order+LandUse",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am2.2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order+LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3.2<- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
              fixedStruct = "Order+LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am4.2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order+LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am1.3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order*LandUse",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am2.3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order*LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3.3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order*LUI",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am4.3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
               fixedStruct = "Order*LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

# take a look at the AICs
AIC_am <-print(AIC(am0$model,am1$model,am2$model,am3$model,am4$model,
                   am0.2$model,am1.2$model,am2.2$model,am3.2$model,am4.2$model,
                   am1.3$model,am2.3$model,am3.3$model,am4.3$model))

write.csv(AIC_am, file = paste0(outDir,"AIC_am.csv"))

# save 
saveRDS(object = sm0 ,file = paste0(outDir,"sm0.rds"))
saveRDS(object = sm3 ,file = paste0(outDir,"sm3.rds"))
saveRDS(object = sm0.2 ,file = paste0(outDir,"sm0.2.rds"))
saveRDS(object = sm3.2 ,file = paste0(outDir,"sm3.2.rds"))
saveRDS(object = sm3.3 ,file = paste0(outDir,"sm3.3.rds"))
saveRDS(object = am0 ,file = paste0(outDir,"am0.rds"))
saveRDS(object = am3 ,file = paste0(outDir,"am3.rds"))
saveRDS(object = am0.2 ,file = paste0(outDir,"am0.2.rds"))
saveRDS(object = am3.2 ,file = paste0(outDir,"am3.2.rds"))
saveRDS(object = am3.3 ,file = paste0(outDir,"am3.3.rds"))

# Plot Figures #

# read in model data
sm0 <- readRDS(file = paste0(outDir,"sm0.rds"))
sm3 <- readRDS(file = paste0(outDir,"sm3.rds"))
sm0.2 <- readRDS(file = paste0(outDir,"sm0.2.rds"))
sm3.2 <- readRDS(file = paste0(outDir,"sm3.2.rds"))
sm3.3 <- readRDS(file = paste0(outDir,"sm3.3.rds"))
am0 <- readRDS(file = paste0(outDir,"am0.rds"))
am3 <- readRDS(file = paste0(outDir,"am3.rds"))
am0.2 <- readRDS(file = paste0(outDir,"am0.2.rds"))
am3.2 <- readRDS(file = paste0(outDir,"am3.2.rds"))
am3.3 <- readRDS(file = paste0(outDir,"am3.3.rds"))
model_data_sr <- readRDS(file = paste0(outDir,"model_data_sr.rds"))
model_data_ab <- readRDS(file = paste0(outDir,"model_data_ab.rds"))

# table of AICs
# species richness and abundance together
selection_table <- data.frame("Response" = c(rep("Species richness", 5),
                                             rep("Abundance", 5)),
                              "Model" = c("Species richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS)",
                                          "Species richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                          "Species richness ~ Order + (1|SS) + (1|SSB) + (1|SSBS)",
                                          "Species richness ~ Order + LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                          "Species richness ~ Order * LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                          "Abundance ~ 1 + (1|SS) + (1|SSB)",
                                          "Abundance ~ LUI + (1|SS) + (1|SSB)",
                                          "Abundance ~ Order + (1|SS) + (1|SSB)",
                                          "Abundance ~ Order + LUI + (1|SS) + (1|SSB)",
                                          "Abundance ~ Order * LUI + (1|SS) + (1|SSB)"),
                              "AIC" = c(AIC(sm0$model), AIC(sm3$model), AIC(sm0.2$model), AIC(sm3.2$model), AIC(sm3.3$model),  
                                        AIC(am0$model), AIC(am3$model), AIC(am0.2$model), AIC(am3.2$model),AIC(am3.3$model))) %>%
  group_by(Response) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  select(Model,AIC,deltaAIC) %>%
  gt(rowname_col = "Model") %>%
       tab_row_group(
         label = "Species Richness",
         rows = starts_with("Species richness")
       ) %>%
       tab_row_group(
         label = "Abundance",
         rows = starts_with("Abundance")
       )%>% 
  cols_align(
    align = "center",
    columns = c(Model, AIC, deltaAIC)
  )%>%
       tab_stubhead(label = "Models")

gtsave(selection_table,"C:/Users/Kyra/Documents/GLITRS/Code/2_RunSimpleLUIModel/Output/LUIModels_Selection1.png")

selection_table_sr <- data.frame("Response" = c(rep("Species richness", 5)),
                                 "Model" = c("Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS)",
                                             "Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                             "Species_richness ~ Order + (1|SS) + (1|SSB) + (1|SSBS)",
                                             "Species_richness ~ Order + LUI + (1|SS) + (1|SSB) + (1|SSBS)",
                                             "Species_richness ~ Order * LUI + (1|SS) + (1|SSB) + (1|SSBS)"),
                                 "AIC" = c(AIC(sm0$model), AIC(sm3$model), AIC(sm0.2$model), AIC(sm3.2$model), AIC(sm3.3$model))) %>%
  group_by(Response) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  gt()

gtsave(selection_table_sr,"C:/Users/Kyra/Documents/GLITRS/Code/2_RunSimpleLUIModel/Output/LUIModels_Selection_Rich.png")

selection_table_ab <- data.frame("Response" = c(rep("Total abundance", 5)),
                                 "Model" = c("Total_abundance ~ 1 + (1|SS) + (1|SSB)",
                                             "Total_abundance ~ LUI + (1|SS) + (1|SSB)",
                                             "Total_abundance ~ Order + (1|SS) + (1|SSB)",
                                             "Total_abundance ~ Order + LUI + (1|SS) + (1|SSB)",
                                             "Total_abundance ~ Order * LUI + (1|SS) + (1|SSB)"),
                                 "AIC" = c(AIC(am0$model), AIC(am3$model), AIC(am0.2$model), AIC(am3.2$model),AIC(am3.3$model))) %>%
  group_by(Response) %>%                              
  mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
  ungroup() %>%
  gt()


gtsave(selection_table_ab,"C:/Users/Kyra/Documents/GLITRS/Code/2_RunSimpleLUIModel/Output/LUIModels_Selection_Abund.png")

# save model output tables for use in supplementary information 
# use function from sjPlot library to save neat versions of model output table
# conditional: the conditional R2 value, i.e. the variance explained by fixed and random effects 
# marginal: the marginal R2 value, i.e. the variance explained by the fixed effects
tab_model(am3.3$model, transform = NULL, file = paste0(outDir,"Abun_output_table.html"))
summary(am3.3$model) # check the table against the outputs
R2GLMER(am3.3$model) # check the R2 values 
# $conditional
# [1] 0.7728262
# 
# $marginal
# [1] 0.09685826

tab_model(sm3.3$model, transform = NULL, file = paste0(outDir,"Rich_output_table.html"))
summary(sm3.3$model) # check the table against the outputs
R2GLMER(sm3.3$model) # check the R2 values 
# $conditional
# [1] 0.7315692
# 
# $marginal
# [1] 0.1149284

## Species Richness Plot ##
richness_metric <- predict_effects(iterations = 1000,
                                   model = sm3.3$model,
                                   model_data = model_data_sr,
                                   response_variable = "Species_richness",
                                   fixed_number = 2,
                                   fixed_column = c("Order", "LUI"),
                                   factor_number_1 = 6,
                                   factor_number_2 = 4,
                                   neg_binom = FALSE)

# rename prediction data frame and drop "Species_richness" column
result.sr <- fin_conf
result.sr <- select(result.sr,-c(Species_richness))

richness_metric

model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_sr$LUI), 6)[1:24]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# richness data
model_data(richness_metric)

## Abundance Plot ##

abundance_metric <- predict_effects(iterations = 1000,
                                    model = am3.3$model,
                                    model_data = model_data_ab,
                                    response_variable = "LogAbund",
                                    fixed_number = 2,
                                    fixed_column = c("Order", "LUI"),
                                    factor_number_1 = 6,
                                    factor_number_2 = 4,
                                    neg_binom = FALSE)

# rename prediction data frame and drop "Abundance" column
result.ab <- fin_conf
result.ab <- select(result.ab,-c(LogAbund))


abundance_metric

model_data <- function(model_plot){
  ggplot_build(model_plot)$data[[2]] %>%
    dplyr::select(y) %>%
    cbind(ggplot_build(model_plot)$data[[3]] %>% dplyr::select(ymin, ymax)) %>%
    mutate(LUI = rep(levels(model_data_ab$LUI), 6)[1:24]) %>%
    dplyr::select(LUI, y, ymin, ymax)
}

# abundance data
model_data(abundance_metric)

# combine results into a table for saving
all_res <- rbind(result.ab, result.sr)
all_res$measure <- c(rep("ab", 24), rep("sr", 24))

# save as table
percentage_change_LUI <- all_res %>% gt()
gtsave(percentage_change_LUI,"C:/Users/Kyra/Documents/GLITRS/Code/7_Predictions/percentage_change_LUI.png")

# save as .csv
write.csv(all_res, file = paste0(predsDir,"percentage_change_LUI.csv"))


# plot
richness <- richness_metric + 
  xlab(NULL) + 
  ylab("Species richness diff. (%)") + 
  guides(scale = "none") + ggtitle("A") + 
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75), limits = c(-100, 75)) + 
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175,200,225,250), limits = c(-100, 250)) + 
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(),
        legend.position = "none")


abundance <- abundance_metric + xlab(NULL) + 
  ylab("Total abundance diff. (%)") + ggtitle("B") + 
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175,200,225,250), limits = c(-100, 250)) + 
  theme(axis.text.x = element_text(size = 9,angle=45,margin=margin(t=20)), 
        axis.ticks = element_blank(), 
        legend.position = "none")


# get the legend
legend <- get_legend(
  abundance +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(), 
          legend.text = element_text(size = 11), 
          legend.title = element_blank())
)

LUI_predictions <- cowplot::plot_grid(richness, abundance, legend,ncol=1, rel_heights = c(1,1.25,0.2))

# save plots
#ggsave("LUI_predictions.jpeg", device ="jpeg", path = outDir, width=25, height=20, units="cm", dpi = 350)
ggsave("LUI_predictions_diffaxes.jpeg", device ="jpeg", path = outDir, width=25, height=20, units="cm", dpi = 350)

# t.end <- Sys.time()
# 
# print(round(t.end - t.start,0))
# 
# sink()
