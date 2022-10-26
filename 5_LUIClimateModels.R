## Get set up ##

# directories 
inDir<- "C:/Users/Kyra/Documents/GLITRS/Code/4_PREDICTSMatchClimateIndex/"
outDir <- "C:/Users/Kyra/Documents/GLITRS/Code/5_RunLUIClimateModels/"
plotDir <- "C:/Users/Kyra/Documents/GLITRS/Code/5_RunLUIClimateModels/Plots/"
if(!dir.exists(outDir)) dir.create(outDir)
if(!dir.exists(plotDir)) dir.create(plotDir)

# sink(paste0(outDir,"log_LUI_ClimateModels.txt"))
# 
# t.start <- Sys.time()
# 
# print(t.start)

# load libraries
packages_model <- c("devtools","StatisticalModels", "predictsFunctions", "ggplot2", "cowplot", "sjPlot","dplyr")
suppressWarnings(suppressMessages(lapply(packages_model, require, character.only = TRUE)))

packages_plot <- c("patchwork", "dplyr", "yarg", "lme4", "gt", "broom.mixed", "MASS")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# source in additional functions
source("C:/Users/Kyra/Documents/GLITRS/Data/0_Functions.R")

###Create Models for all insect orders in predicts for standardised climate anomaly and Land interactions ###

# read in the predicts data
predictsSites <- readRDS(paste0(inDir,"PREDICTSSites_Climate.rds"))
predictsSites <- predictsSites@data

# remove groups Blattodea, Neuroptera, Other, Thysanoptera, Trichoptera 
# only necessary if these were not already removed when preparing the dataset
predictsSites <- predictsSites %>% filter(Order %in% c("Hymenoptera", "Coleoptera", "Lepidoptera", "Diptera", "Orthoptera", "Hemiptera")) %>% droplevels()

# set LUI as factor and set reference level
predictsSites$LUI <- factor(predictsSites$LUI)
predictsSites$LUI <- relevel(predictsSites$LUI,ref="Primary vegetation")

# rescale the climate variables
predictsSites$StdTmeanAnomalyRS <- StdCenterPredictor(predictsSites$StdTmeanAnomaly)
predictsSites$StdTmaxAnomalyRS <- StdCenterPredictor(predictsSites$StdTmaxAnomaly)

# rescaling abundance and log values
predictsSites <- RescaleAbundance(predictsSites)

# Charlie added this line as later bits were throwing errors
predictsSites <- droplevels(predictsSites)

# some of the climate values are NA since they do not meet the thresholds
predictsSites <- predictsSites[!is.na(predictsSites$avg_temp), ]

# take a look at possible correlations between variables
cor(predictsSites$avg_temp, predictsSites$TmeanAnomaly)

# -0.246491

cor(predictsSites$avg_temp, predictsSites$StdTmeanAnomaly)

#0.2118244

cor(predictsSites$TmeanAnomaly, predictsSites$StdTmeanAnomaly)

# 0.2639987

# save the dataset
saveRDS(object = predictsSites,file = paste0(outDir,"PREDICTSSitesClimate_Data.rds"))

# predictsSites <- readRDS(file = paste0(outDir,"PREDICTSSitesClimate_Data.rds"))


## Model Selection ##

# 1. Abundance, mean anomaly

model_data <- predictsSites[!is.na(predictsSites$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MeanAnomalyModelAbund <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                               fixedStruct = "LUI * StdTmeanAnomalyRS * Order",
                               randomStruct = "(1|SS)+(1|SSB)",
                               saveVars = c("SSBS"))

# get summary
summary(MeanAnomalyModelAbund$model)


# save the model output
save(MeanAnomalyModelAbund, file = paste0(outDir, "MeanAnomalyModelAbund.rdata"))


# 2. Richness, mean anomaly

model_data <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ]

MeanAnomalyModelRich <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                              fixedStruct = "LUI * StdTmeanAnomalyRS * Order",
                              randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                              saveVars = c("SSBS"))

summary(MeanAnomalyModelRich$model)

save(MeanAnomalyModelRich, file = paste0(outDir, "MeanAnomalyModelRich.rdata"))

# 3. Abundance, max anomaly

model_data <- predictsSites[!is.na(predictsSites$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]


MaxAnomalyModelAbund <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
                              fixedStruct = "LUI * StdTmaxAnomalyRS * Order",
                              randomStruct = "(1|SS)+(1|SSB)",
                              saveVars = c("SSBS"))

summary(MaxAnomalyModelAbund$model)

save(MaxAnomalyModelAbund, file = paste0(outDir, "MaxAnomalyModelAbund.rdata"))

# 4. Richness, max anomaly

model_data <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS),]

MaxAnomalyModelRich <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                             fixedStruct = "LUI * StdTmaxAnomalyRS * Order",
                             randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                             saveVars = c("SSBS"))

# save model output
save(MaxAnomalyModelRich, file = paste0(outDir, "MaxAnomalyModelRich.rdata"))


#### save model output tables for use in supplementary information ####
# use function from sjPlot library to save neat versions of model output table
# conditional: the conditional R2 value, i.e. the variance explained by fixed and random effects 
# marginal: the marginal R2 value, i.e. the variance explained by the fixed effects

tab_model(MeanAnomalyModelAbund$model, transform = NULL, file = paste0(outDir,"Tables/AbunMeanAnom_output_table.html"))
summary(MeanAnomalyModelAbund$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelAbund$model) # check the R2 values 

tab_model(MeanAnomalyModelRich$model, transform = NULL, file = paste0(outDir,"Tables/RichMeanAnom_output_table.html"))
summary(MeanAnomalyModelRich$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelRich$model) # check the R2 values

tab_model(MaxAnomalyModelAbund$model, transform = NULL, file = paste0(outDir,"Tables/AbunMaxAnom_output_table.html"))
summary(MaxAnomalyModelAbund$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelAbund$model) # check the R2 values 

tab_model(MaxAnomalyModelRich$model, transform = NULL, file = paste0(outDir,"Tables/RichMaxAnom_output_table.html"))
summary(MaxAnomalyModelRich$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelRich$model) # check the R2 values 

## Plot results ##

# load models
predictsSites <- readRDS(file = paste0(outDir,"PREDICTSSitesClimate_Data.rds"))
load(paste0(outDir, "MeanAnomalyModelAbund.rdata"))
load(paste0(outDir, "MeanAnomalyModelRich.rdata"))
load(paste0(outDir, "MaxAnomalyModelAbund.rdata"))
load(paste0(outDir, "MaxAnomalyModelRich.rdata"))


# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

## Abundance, Mean Anomaly ##

nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
                        length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")))

# have to maintain Orthoptera, b/c it's part of the model and the predictions won't run without it
# remove later, when plotting

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# RECORD WHICH ROW IS REFERENCE FOR LATER (see 'Values')
refRow <- which((nd$LUI=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))
# 4th row, every 400 rows 

# set quantiles
QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$LUI=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# create list of matrices
number_of_chunks = 6
list_a.preds.tmean <- lapply(seq(1, NROW(a.preds.tmean), ceiling(NROW(a.preds.tmean)/number_of_chunks)),
       function(i) a.preds.tmean[i:min(i + ceiling(NROW(a.preds.tmean)/number_of_chunks) - 1, NROW(a.preds.tmean)),])

# name them
names(list_a.preds.tmean) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")
                             
# keep all but Orthoptera
list_a.preds.tmean <- list_a.preds.tmean[names(list_a.preds.tmean) %in% c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")]
                      
# tim's suggestion
list_a.preds.tmean <- lapply(list_a.preds.tmean,FUN=function(x){
 sweep (x=x, MARGIN = 2, STATS=x[4,],FUN="/") 
})

# print to global environment
list2env(list_a.preds.tmean,globalenv())
                             
# can now remove the extra orders from nd
nd <- filter(nd, Order %in% c('Coleoptera', 'Diptera', 'Hemiptera','Hymenoptera', 'Lepidoptera'))

# split nd by order
Order<- paste0("nd_",nd$Order)
# create a list of data frames
by_Order <- split(nd,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles

Coleoptera[which((nd_Coleoptera$LUI=="Primary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Primary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Secondary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Secondary vegetation") & (nd_Coleoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_Low") & (nd_Coleoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_Low") & (nd_Coleoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_High") & (nd_Coleoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd_Coleoptera$LUI=="Agriculture_High") & (nd_Coleoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd_Diptera$LUI=="Primary vegetation") & (nd_Diptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd_Diptera$LUI=="Primary vegetation") & (nd_Diptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Secondary vegetation") & (nd_Diptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Secondary vegetation") & (nd_Diptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_Low") & (nd_Diptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_Low") & (nd_Diptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_High") & (nd_Diptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd_Diptera$LUI=="Agriculture_High") & (nd_Diptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd_Hemiptera$LUI=="Primary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Primary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Secondary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Secondary vegetation") & (nd_Hemiptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_Low") & (nd_Hemiptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_Low") & (nd_Hemiptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_High") & (nd_Hemiptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd_Hemiptera$LUI=="Agriculture_High") & (nd_Hemiptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd_Hymenoptera$LUI=="Primary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Primary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Secondary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Secondary vegetation") & (nd_Hymenoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_Low") & (nd_Hymenoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_Low") & (nd_Hymenoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_High") & (nd_Hymenoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd_Hymenoptera$LUI=="Agriculture_High") & (nd_Hymenoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd_Lepidoptera$LUI=="Primary vegetation") & (nd_Lepidoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Primary vegetation") & (nd_Lepidoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Secondary vegetation") & (nd_Lepidoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Secondary vegetation") & (nd_Lepidoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Agriculture_Low") & (nd_Lepidoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Agriculture_Low") & (nd_Lepidoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Agriculture_High") & (nd_Lepidoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd_Lepidoptera$LUI=="Agriculture_High") & (nd_Lepidoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                    FUN = median,na.rm=TRUE))*100)-100
nd_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                 FUN = median,na.rm=TRUE))*100)-100
nd_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                   FUN = median,na.rm=TRUE))*100)-100
nd_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                  FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                  FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
nd_Coleoptera$LUI <- factor(nd_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Diptera$LUI <- factor(nd_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Hemiptera$LUI <- factor(nd_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Hymenoptera$LUI <- factor(nd_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_Lepidoptera$LUI <- factor(nd_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot

p_coleoptera <- ggplot(data = nd_Coleoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Coleoptera$PredLower, ymax = nd_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Coleoptera")

p_diptera <- ggplot(data = nd_Diptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Diptera$PredLower, ymax = nd_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Diptera")

p_hemiptera <- ggplot(data = nd_Hemiptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Hemiptera$PredLower, ymax = nd_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Hemiptera")

p_hymenoptera <- ggplot(data = nd_Hymenoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Hymenoptera$PredLower, ymax = nd_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Hymenoptera")

p_lepidoptera <- ggplot(data = nd_Lepidoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Lepidoptera$PredLower, ymax = nd_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Lepidoptera")

# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(), 
          legend.text = element_text(size = 11), 
          legend.title = element_blank())
)


# put them all together to save them
MeanAnomAbund <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera)
MeanAnomAbund <- cowplot::plot_grid(MeanAnomAbund,legend,ncol=1, rel_heights = c(1,0.1))

# save the ggplots
ggsave(filename = paste0(plotDir, "MeanAnomAbund.pdf"), plot = MeanAnomAbund, width = 200, height = 150, units = "mm", dpi = 300)
# ggsave(filename = paste0(plotDir, "MeanAnomAbund_extended yaxis.pdf"), plot = MeanAnomAbund, width = 200, height = 150, units = "mm", dpi = 300)

## Richness, Mean Anomaly ##

nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
                        length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelRich$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")))

# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# RECORD WHICH ROW IS REFERENCE FOR LATER (see 'Values')
# reference row is 4th row, every 400 rows
refRow <- which((nd2$LUI=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))))

# set quantiles
QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$LUI=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
sr.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model,data = nd2)

# back transform the abundance values
sr.preds.tmean <- exp(sr.preds.tmean)-0.01

# create list of matrices
number_of_chunks = 6
list_sr.preds.tmean <- lapply(seq(1, NROW(sr.preds.tmean), ceiling(NROW(sr.preds.tmean)/number_of_chunks)),
                              function(i) sr.preds.tmean[i:min(i + ceiling(NROW(sr.preds.tmean)/number_of_chunks) - 1, NROW(sr.preds.tmean)),])
# name them
names(list_sr.preds.tmean) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")

# keep all but Orthoptera
list_sr.preds.tmean <- list_sr.preds.tmean[names(list_sr.preds.tmean) %in% c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")]

# tim's suggestion
list_sr.preds.tmean <- lapply(list_sr.preds.tmean,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[4,],FUN="/") 
})
                              
# can now remove the extra orders from nd2
nd2 <- filter(nd2, Order %in% c('Coleoptera', 'Diptera', 'Hemiptera','Hymenoptera', 'Lepidoptera'))                             

list2env(list_sr.preds.tmean,globalenv())

# split nd by order
Order<- paste0("nd2_",nd2$Order)
# create a list of data frames
by_Order <- split(nd2,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles
Coleoptera[which((nd2_Coleoptera$LUI=="Primary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Primary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Secondary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Secondary vegetation") & (nd2_Coleoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_Low") & (nd2_Coleoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_Low") & (nd2_Coleoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_High") & (nd2_Coleoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd2_Coleoptera$LUI=="Agriculture_High") & (nd2_Coleoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd2_Diptera$LUI=="Primary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Primary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Secondary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Secondary vegetation") & (nd2_Diptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Agriculture_Low") & (nd2_Diptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Agriculture_Low") & (nd2_Diptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Agriculture_High") & (nd2_Diptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd2_Diptera$LUI=="Agriculture_High") & (nd2_Diptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd2_Hemiptera$LUI=="Primary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Primary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Secondary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Secondary vegetation") & (nd2_Hemiptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_Low") & (nd2_Hemiptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_Low") & (nd2_Hemiptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_High") & (nd2_Hemiptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd2_Hemiptera$LUI=="Agriculture_High") & (nd2_Hemiptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd2_Hymenoptera$LUI=="Primary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Primary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Secondary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Secondary vegetation") & (nd2_Hymenoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_Low") & (nd2_Hymenoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_Low") & (nd2_Hymenoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_High") & (nd2_Hymenoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd2_Hymenoptera$LUI=="Agriculture_High") & (nd2_Hymenoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd2_Lepidoptera$LUI=="Primary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Primary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Secondary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Secondary vegetation") & (nd2_Lepidoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_Low") & (nd2_Lepidoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_Low") & (nd2_Lepidoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_High") & (nd2_Lepidoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd2_Lepidoptera$LUI=="Agriculture_High") & (nd2_Lepidoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot
nd2_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd2_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                  FUN = median,na.rm=TRUE))*100)-100
nd2_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                 FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                 FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                    FUN = median,na.rm=TRUE))*100)-100
nd2_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                      FUN = median,na.rm=TRUE))*100)-100
nd2_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                      FUN = median,na.rm=TRUE))*100)-100
nd2_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
nd2_Coleoptera$LUI <- factor(nd2_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_Diptera$LUI <- factor(nd2_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_Hemiptera$LUI <- factor(nd2_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_Hymenoptera$LUI <- factor(nd2_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_Lepidoptera$LUI <- factor(nd2_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot
p_coleoptera <- ggplot(data = nd2_Coleoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Coleoptera$PredLower, ymax = nd2_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Coleoptera")

p_diptera <- ggplot(data = nd2_Diptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Diptera$PredLower, ymax = nd2_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Diptera")

p_hemiptera <- ggplot(data = nd2_Hemiptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Hemiptera$PredLower, ymax = nd2_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Hemiptera")

p_hymenoptera <- ggplot(data = nd2_Hymenoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Hymenoptera$PredLower, ymax = nd2_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Hymenoptera")

p_lepidoptera <- ggplot(data = nd2_Lepidoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Lepidoptera$PredLower, ymax = nd2_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Lepidoptera")


# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(), 
          legend.text = element_text(size = 11), 
          legend.title = element_blank())
)


# put them all together to save them
MeanAnomRich <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera)
MeanAnomRich <- cowplot::plot_grid(MeanAnomRich,legend,ncol=1, rel_heights = c(1,0.1))

# save the ggplots
ggsave(filename = paste0(plotDir, "MeanAnomRich.pdf"), plot = MeanAnomRich, width = 200, height = 150, units = "mm", dpi = 300)
# ggsave(filename = paste0(plotDir, "MeanAnomRich_extended yaxis.pdf"), plot = MeanAnomRich, width = 200, height = 150, units = "mm", dpi = 300)

## Abundance, Max Anomaly ##

nd3 <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
                       to = max(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
                       length.out = 200),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelAbund$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")))


nd3$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

nd3$LogAbund <- 0
nd3$Species_richness <- 0

# Check the reference row manually and record it for later
refRow <- which((nd3$LUI=="Primary vegetation") & (nd3$StdTmaxAnomaly==min(abs(nd3$StdTmaxAnomaly))))


QPV <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund$data$LUI=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund$data$LUI=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund$data$LUI=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund$data$LUI=="Agriculture_High"],
  probs = exclQuantiles)

# if(!is.null(MaxAnomalyModelAbund$model)){

a.preds.tmax <- PredictGLMERRandIter(model = MaxAnomalyModelAbund$model,data = nd3, nIters = 10000)
a.preds.tmax <- exp(a.preds.tmax)-0.01

# split into 6 groups
number_of_chunks = 6
list_a.preds.tmax <- lapply(seq(1, NROW(a.preds.tmax), ceiling(NROW(a.preds.tmax)/number_of_chunks)),
                            function(i) a.preds.tmax[i:min(i + ceiling(NROW(a.preds.tmax)/number_of_chunks) - 1, NROW(a.preds.tmax)),])
# name the matrices
names(list_a.preds.tmax) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")
                            
# keep all but Orthoptera
list_a.preds.tmax <- list_a.preds.tmax[names(list_a.preds.tmax) %in% c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")]                            

# unload into global environment
list2env(list_a.preds.tmax,globalenv())

# sweep out according to the 114th row
list_a.preds.tmax <- lapply(list_a.preds.tmax,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[114,],FUN="/") 
})

list2env(list_a.preds.tmax,globalenv())
                            
 # can now remove the extra orders from nd3
nd3 <- filter(nd3, Order %in% c('Coleoptera', 'Diptera', 'Hemiptera','Hymenoptera', 'Lepidoptera'))

# split nd by order
Order<- paste0("nd3_",nd3$Order)
# create a list of data frames
by_Order <- split(nd3,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles
Coleoptera[which((nd3_Coleoptera$LUI=="Primary vegetation") & (nd3_Coleoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd3_Coleoptera$LUI=="Primary vegetation") & (nd3_Coleoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd3_Coleoptera$LUI=="Secondary vegetation") & (nd3_Coleoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd3_Coleoptera$LUI=="Secondary vegetation") & (nd3_Coleoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd3_Coleoptera$LUI=="Agriculture_Low") & (nd3_Coleoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd3_Coleoptera$LUI=="Agriculture_Low") & (nd3_Coleoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd3_Coleoptera$LUI=="Agriculture_High") & (nd3_Coleoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd3_Coleoptera$LUI=="Agriculture_High") & (nd3_Coleoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd3_Diptera$LUI=="Primary vegetation") & (nd3_Diptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd3_Diptera$LUI=="Primary vegetation") & (nd3_Diptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd3_Diptera$LUI=="Secondary vegetation") & (nd3_Diptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd3_Diptera$LUI=="Secondary vegetation") & (nd3_Diptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd3_Diptera$LUI=="Agriculture_Low") & (nd3_Diptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd3_Diptera$LUI=="Agriculture_Low") & (nd3_Diptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd3_Diptera$LUI=="Agriculture_High") & (nd3_Diptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd3_Diptera$LUI=="Agriculture_High") & (nd3_Diptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd3_Hemiptera$LUI=="Primary vegetation") & (nd3_Hemiptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd3_Hemiptera$LUI=="Primary vegetation") & (nd3_Hemiptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd3_Hemiptera$LUI=="Secondary vegetation") & (nd3_Hemiptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd3_Hemiptera$LUI=="Secondary vegetation") & (nd3_Hemiptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd3_Hemiptera$LUI=="Agriculture_Low") & (nd3_Hemiptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd3_Hemiptera$LUI=="Agriculture_Low") & (nd3_Hemiptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd3_Hemiptera$LUI=="Agriculture_High") & (nd3_Hemiptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd3_Hemiptera$LUI=="Agriculture_High") & (nd3_Hemiptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd3_Hymenoptera$LUI=="Primary vegetation") & (nd3_Hymenoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd3_Hymenoptera$LUI=="Primary vegetation") & (nd3_Hymenoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd3_Hymenoptera$LUI=="Secondary vegetation") & (nd3_Hymenoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd3_Hymenoptera$LUI=="Secondary vegetation") & (nd3_Hymenoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd3_Hymenoptera$LUI=="Agriculture_Low") & (nd3_Hymenoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd3_Hymenoptera$LUI=="Agriculture_Low") & (nd3_Hymenoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd3_Hymenoptera$LUI=="Agriculture_High") & (nd3_Hymenoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd3_Hymenoptera$LUI=="Agriculture_High") & (nd3_Hymenoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd3_Lepidoptera$LUI=="Primary vegetation") & (nd3_Lepidoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd3_Lepidoptera$LUI=="Primary vegetation") & (nd3_Lepidoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd3_Lepidoptera$LUI=="Secondary vegetation") & (nd3_Lepidoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd3_Lepidoptera$LUI=="Secondary vegetation") & (nd3_Lepidoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd3_Lepidoptera$LUI=="Agriculture_Low") & (nd3_Lepidoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd3_Lepidoptera$LUI=="Agriculture_Low") & (nd3_Lepidoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd3_Lepidoptera$LUI=="Agriculture_High") & (nd3_Lepidoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd3_Lepidoptera$LUI=="Agriculture_High") & (nd3_Lepidoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot
nd3_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd3_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                  FUN = median,na.rm=TRUE))*100)-100
nd3_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                 FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                 FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                    FUN = median,na.rm=TRUE))*100)-100
nd3_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                      FUN = median,na.rm=TRUE))*100)-100
nd3_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                      FUN = median,na.rm=TRUE))*100)-100
nd3_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# set factor levels
# nd$LUI <- factor(nd$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_Coleoptera$LUI <- factor(nd3_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_Diptera$LUI <- factor(nd3_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_Hemiptera$LUI <- factor(nd3_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_Hymenoptera$LUI <- factor(nd3_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_Lepidoptera$LUI <- factor(nd3_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot
p_coleoptera <- ggplot(data = nd3_Coleoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Coleoptera$PredLower, ymax = nd3_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150), limits = c(-100, 150)) +
  #scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150,175, 200), limits = c(-100, 200)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Coleoptera")

p_diptera <- ggplot(data = nd3_Diptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Diptera$PredLower, ymax = nd3_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150), limits = c(-100, 150)) +
  #scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150,175, 200), limits = c(-100, 200)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Diptera")

p_hemiptera <- ggplot(data = nd3_Hemiptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Hemiptera$PredLower, ymax = nd3_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150), limits = c(-100, 150)) +
  #scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150,175, 200), limits = c(-100, 200)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Hemiptera")

p_hymenoptera <- ggplot(data = nd3_Hymenoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Hymenoptera$PredLower, ymax = nd3_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150), limits = c(-100, 150)) +
  #scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150,175, 200), limits = c(-100, 200)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Hymenoptera")

p_lepidoptera <- ggplot(data = nd3_Lepidoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Lepidoptera$PredLower, ymax = nd3_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150), limits = c(-100, 150)) +
  #scale_y_continuous(breaks = c(-100, -75, -50, -25, 0, 25, 50, 75, 100, 125, 150,175, 200), limits = c(-100, 200)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Lepidoptera")

# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(), 
          legend.text = element_text(size = 11), 
          legend.title = element_blank())
)


# put them all together to save them
MaxAnomAbund <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera)
MaxAnomAbund <- cowplot::plot_grid(MaxAnomAbund,legend,ncol=1, rel_heights = c(1,0.1))

# save the ggplots
ggsave(filename = paste0(plotDir, "MaxAnomAbund.pdf"), plot = MaxAnomAbund, width = 200, height = 150, units = "mm", dpi = 300)

## Richness, Max Anomaly ##

nd4 <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelRich$data$StdTmaxAnomalyRS),
                       to = max(MaxAnomalyModelRich$data$StdTmaxAnomalyRS),
                       length.out = 200),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelRich$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")))


nd4$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

nd4$LogAbund <- 0
nd4$Species_richness <- 0

# Check the reference row manually and record it for later
refRow <- which((nd4$LUI=="Primary vegetation") & (nd4$StdTmaxAnomaly==min(abs(nd4$StdTmaxAnomaly))))
# 114th row

QPV <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich$data$LUI=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich$data$LUI=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich$data$LUI=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich$data$LUI=="Agriculture_High"],
  probs = exclQuantiles)

sr.preds.tmax <- PredictGLMERRandIter(model = MaxAnomalyModelRich$model,data = nd4, nIters = 10000)
sr.preds.tmax <- exp(sr.preds.tmax)-0.01

# split into 6 groups
number_of_chunks = 6
list_sr.preds.tmax <- lapply(seq(1, NROW(sr.preds.tmax), ceiling(NROW(sr.preds.tmax)/number_of_chunks)),
                             function(i) sr.preds.tmax[i:min(i + ceiling(NROW(sr.preds.tmax)/number_of_chunks) - 1, NROW(sr.preds.tmax)),])
# name the matrices
names(list_sr.preds.tmax) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")
                             
# keep all but Orthoptera
list_sr.preds.tmax <- list_sr.preds.tmax[names(list_sr.preds.tmax) %in% c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera")]                             

# unload into global environment
list2env(list_sr.preds.tmax,globalenv())

# sweep out according to the nth row (Check manually)
list_sr.preds.tmax <- lapply(list_sr.preds.tmax,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[114,],FUN="/") 
})

list2env(list_sr.preds.tmax,globalenv())
                             
# can now remove the extra orders from nd4
nd4 <- filter(nd4, Order %in% c('Coleoptera', 'Diptera', 'Hemiptera','Hymenoptera', 'Lepidoptera'))                             

# split nd by order
Order<- paste0("nd4_",nd4$Order)
# create a list of data frames
by_Order <- split(nd4,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles
Coleoptera[which((nd4_Coleoptera$LUI=="Primary vegetation") & (nd4_Coleoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd4_Coleoptera$LUI=="Primary vegetation") & (nd4_Coleoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd4_Coleoptera$LUI=="Secondary vegetation") & (nd4_Coleoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd4_Coleoptera$LUI=="Secondary vegetation") & (nd4_Coleoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd4_Coleoptera$LUI=="Agriculture_Low") & (nd4_Coleoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd4_Coleoptera$LUI=="Agriculture_Low") & (nd4_Coleoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd4_Coleoptera$LUI=="Agriculture_High") & (nd4_Coleoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd4_Coleoptera$LUI=="Agriculture_High") & (nd4_Coleoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd4_Diptera$LUI=="Primary vegetation") & (nd4_Diptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd4_Diptera$LUI=="Primary vegetation") & (nd4_Diptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd4_Diptera$LUI=="Secondary vegetation") & (nd4_Diptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd4_Diptera$LUI=="Secondary vegetation") & (nd4_Diptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd4_Diptera$LUI=="Agriculture_Low") & (nd4_Diptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd4_Diptera$LUI=="Agriculture_Low") & (nd4_Diptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd4_Diptera$LUI=="Agriculture_High") & (nd4_Diptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd4_Diptera$LUI=="Agriculture_High") & (nd4_Diptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd4_Hemiptera$LUI=="Primary vegetation") & (nd4_Hemiptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd4_Hemiptera$LUI=="Primary vegetation") & (nd4_Hemiptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd4_Hemiptera$LUI=="Secondary vegetation") & (nd4_Hemiptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd4_Hemiptera$LUI=="Secondary vegetation") & (nd4_Hemiptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd4_Hemiptera$LUI=="Agriculture_Low") & (nd4_Hemiptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd4_Hemiptera$LUI=="Agriculture_Low") & (nd4_Hemiptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd4_Hemiptera$LUI=="Agriculture_High") & (nd4_Hemiptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd4_Hemiptera$LUI=="Agriculture_High") & (nd4_Hemiptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd4_Hymenoptera$LUI=="Primary vegetation") & (nd4_Hymenoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd4_Hymenoptera$LUI=="Primary vegetation") & (nd4_Hymenoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd4_Hymenoptera$LUI=="Secondary vegetation") & (nd4_Hymenoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd4_Hymenoptera$LUI=="Secondary vegetation") & (nd4_Hymenoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd4_Hymenoptera$LUI=="Agriculture_Low") & (nd4_Hymenoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd4_Hymenoptera$LUI=="Agriculture_Low") & (nd4_Hymenoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd4_Hymenoptera$LUI=="Agriculture_High") & (nd4_Hymenoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd4_Hymenoptera$LUI=="Agriculture_High") & (nd4_Hymenoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd4_Lepidoptera$LUI=="Primary vegetation") & (nd4_Lepidoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd4_Lepidoptera$LUI=="Primary vegetation") & (nd4_Lepidoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd4_Lepidoptera$LUI=="Secondary vegetation") & (nd4_Lepidoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd4_Lepidoptera$LUI=="Secondary vegetation") & (nd4_Lepidoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd4_Lepidoptera$LUI=="Agriculture_Low") & (nd4_Lepidoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd4_Lepidoptera$LUI=="Agriculture_Low") & (nd4_Lepidoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd4_Lepidoptera$LUI=="Agriculture_High") & (nd4_Lepidoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd4_Lepidoptera$LUI=="Agriculture_High") & (nd4_Lepidoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd4_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd4_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                  FUN = median,na.rm=TRUE))*100)-100
nd4_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                 FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                 FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                    FUN = median,na.rm=TRUE))*100)-100
nd4_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                      FUN = median,na.rm=TRUE))*100)-100
nd4_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                      FUN = median,na.rm=TRUE))*100)-100
nd4_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

#set factor levels
nd4_Coleoptera$LUI <- factor(nd4_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_Diptera$LUI <- factor(nd4_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_Hemiptera$LUI <- factor(nd4_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_Hymenoptera$LUI <- factor(nd4_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_Lepidoptera$LUI <- factor(nd4_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot
p_coleoptera <- ggplot(data = nd4_Coleoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Coleoptera$PredLower, ymax = nd4_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250), limits = c(-100, 200)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Coleoptera")

p_diptera <- ggplot(data = nd4_Diptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Diptera$PredLower, ymax = nd4_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250), limits = c(-100, 200)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Diptera")

p_hemiptera <- ggplot(data = nd4_Hemiptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Hemiptera$PredLower, ymax = nd4_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250), limits = c(-100, 200)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Hemiptera")

p_hymenoptera <- ggplot(data = nd4_Hymenoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Hymenoptera$PredLower, ymax = nd4_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250), limits = c(-100, 200)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Hymenoptera")

p_lepidoptera <- ggplot(data = nd4_Lepidoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Lepidoptera$PredLower, ymax = nd4_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250), limits = c(-100, 200)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        #legend.position = c(0.2, 0.8),
        #legend.background = element_blank(), 
        #legend.text = element_text(size = 6), 
        #legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Lepidoptera")

# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(), 
          legend.text = element_text(size = 11), 
          legend.title = element_blank())
)


# put them all together to save them
MaxAnomRich <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera)
MaxAnomRich <- cowplot::plot_grid(MaxAnomRich,legend,ncol=1, rel_heights = c(1,0.1))

# save the ggplots
ggsave(filename = paste0(plotDir, "MaxAnomRich.pdf"), plot = MaxAnomRich, width = 200, height = 150, units = "mm", dpi = 300)
#ggsave(filename = paste0(plotDir, "MaxAnomRich_extended yaxis.pdf"), plot = MaxAnomRich, width = 200, height = 150, units = "mm", dpi = 300)


# t.end <- Sys.time()
# 
# print(round(t.end - t.start,0))
# 
# sink()
