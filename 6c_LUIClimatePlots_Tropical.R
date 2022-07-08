
## Plot results - separate plots ##

# get set up
# directories 
predictsDir<- "C:/Users/Kyra/Documents/GLITRS/Code/5_RunLUIClimateModels/"
inDir<- "C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Output/"
outDir <- "C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Plots/"
if(!dir.exists(outDir)) dir.create(outDir)

# load libraries
library(devtools)
library(StatisticalModels)
library(predictsFunctions)
library(sjPlot)
library(cowplot)

packages_plot <- c("devtools","StatisticalModels", "predictsFunctions","cowplot", "sjPlot")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# source additional functions
source("C:/Users/Kyra/Documents/GLITRS/Data/0_Functions.R")

# load data sets and models
predictsSites <- readRDS(file = paste0(predictsDir,"PREDICTSSitesClimate_Data.rds"))
trop <- readRDS(file = paste0(inDir,"trop.rds"))
nontrop <- readRDS(file = paste0(inDir,"nontrop.rds"))
load(paste0(inDir, "MeanAnomalyModelAbund_trop.rdata"))
load(paste0(inDir, "MeanAnomalyModelRich_trop.rdata"))
load(paste0(inDir, "MaxAnomalyModelAbund_trop.rdata"))
load(paste0(inDir, "MaxAnomalyModelRich_trop.rdata"))
load(paste0(inDir, "MeanAnomalyModelAbund_nontrop.rdata"))
load(paste0(inDir, "MeanAnomalyModelRich_nontrop.rdata"))
load(paste0(inDir, "MaxAnomalyModelAbund_nontrop.rdata"))
load(paste0(inDir, "MaxAnomalyModelRich_nontrop.rdata"))

# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

## Abundance, Mean Anomaly ##
## Tropical ##

nd_trop <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund_trop$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund_trop$data$StdTmeanAnomalyRS),
                        length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund_trop$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")))

# back transform the predictors
nd_trop$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd_trop$StdTmeanAnomalyRS,
  originalX = trop$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd_trop$LogAbund <- 0
nd_trop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# Record refRow for later use (see 'Values')
refRow <- which((nd_trop$LUI=="Primary vegetation") & (nd_trop$StdTmeanAnomaly==min(abs(nd_trop$StdTmeanAnomaly))))
# 1st row of each order

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelAbund_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_trop$data$LUI=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelAbund_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_trop$data$LUI=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelAbund_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_trop$data$LUI=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelAbund_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_trop$data$LUI=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean.trop <- PredictGLMERRandIter(model = MeanAnomalyModelAbund_trop$model,data = nd_trop)

# back transform the abundance values
a.preds.tmean.trop <- exp(a.preds.tmean.trop)-0.01


# split by order into matrices, then name them
number_of_chunks = 6
list_a.preds.tmean.trop <- lapply(seq(1, NROW(a.preds.tmean.trop), ceiling(NROW(a.preds.tmean.trop)/number_of_chunks)),
                                  function(i) a.preds.tmean.trop[i:min(i + ceiling(NROW(a.preds.tmean.trop)/number_of_chunks) - 1, NROW(a.preds.tmean.trop)),])

names(list_a.preds.tmean.trop) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")
list2env(list_a.preds.tmean.trop,globalenv())

# sweep out refRow
list_a.preds.tmean.trop <- lapply(list_a.preds.tmean.trop,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[1,],FUN="/") 
})

list2env(list_a.preds.tmean.trop,globalenv())

# split nd_trop by order
Order<- paste0("nd_trop_",nd_trop$Order)
# create a list of data frames
by_Order <- split(nd_trop,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles
Coleoptera[which((nd_trop_Coleoptera$LUI=="Primary vegetation") & (nd_trop_Coleoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd_trop_Coleoptera$LUI=="Primary vegetation") & (nd_trop_Coleoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd_trop_Coleoptera$LUI=="Second_tropary vegetation") & (nd_trop_Coleoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd_trop_Coleoptera$LUI=="Second_tropary vegetation") & (nd_trop_Coleoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd_trop_Coleoptera$LUI=="Agriculture_Low") & (nd_trop_Coleoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd_trop_Coleoptera$LUI=="Agriculture_Low") & (nd_trop_Coleoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd_trop_Coleoptera$LUI=="Agriculture_High") & (nd_trop_Coleoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd_trop_Coleoptera$LUI=="Agriculture_High") & (nd_trop_Coleoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd_trop_Diptera$LUI=="Primary vegetation") & (nd_trop_Diptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd_trop_Diptera$LUI=="Primary vegetation") & (nd_trop_Diptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd_trop_Diptera$LUI=="Second_tropary vegetation") & (nd_trop_Diptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd_trop_Diptera$LUI=="Second_tropary vegetation") & (nd_trop_Diptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd_trop_Diptera$LUI=="Agriculture_Low") & (nd_trop_Diptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd_trop_Diptera$LUI=="Agriculture_Low") & (nd_trop_Diptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd_trop_Diptera$LUI=="Agriculture_High") & (nd_trop_Diptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd_trop_Diptera$LUI=="Agriculture_High") & (nd_trop_Diptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd_trop_Hemiptera$LUI=="Primary vegetation") & (nd_trop_Hemiptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd_trop_Hemiptera$LUI=="Primary vegetation") & (nd_trop_Hemiptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd_trop_Hemiptera$LUI=="Second_tropary vegetation") & (nd_trop_Hemiptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd_trop_Hemiptera$LUI=="Second_tropary vegetation") & (nd_trop_Hemiptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd_trop_Hemiptera$LUI=="Agriculture_Low") & (nd_trop_Hemiptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd_trop_Hemiptera$LUI=="Agriculture_Low") & (nd_trop_Hemiptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd_trop_Hemiptera$LUI=="Agriculture_High") & (nd_trop_Hemiptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd_trop_Hemiptera$LUI=="Agriculture_High") & (nd_trop_Hemiptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd_trop_Hymenoptera$LUI=="Primary vegetation") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd_trop_Hymenoptera$LUI=="Primary vegetation") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd_trop_Hymenoptera$LUI=="Second_tropary vegetation") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd_trop_Hymenoptera$LUI=="Second_tropary vegetation") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd_trop_Hymenoptera$LUI=="Agriculture_Low") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd_trop_Hymenoptera$LUI=="Agriculture_Low") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd_trop_Hymenoptera$LUI=="Agriculture_High") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd_trop_Hymenoptera$LUI=="Agriculture_High") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd_trop_Lepidoptera$LUI=="Primary vegetation") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd_trop_Lepidoptera$LUI=="Primary vegetation") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd_trop_Lepidoptera$LUI=="Second_tropary vegetation") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd_trop_Lepidoptera$LUI=="Second_tropary vegetation") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd_trop_Lepidoptera$LUI=="Agriculture_Low") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd_trop_Lepidoptera$LUI=="Agriculture_Low") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd_trop_Lepidoptera$LUI=="Agriculture_High") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd_trop_Lepidoptera$LUI=="Agriculture_High") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd_trop_Orthoptera$LUI=="Primary vegetation") & (nd_trop_Orthoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd_trop_Orthoptera$LUI=="Primary vegetation") & (nd_trop_Orthoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd_trop_Orthoptera$LUI=="Second_tropary vegetation") & (nd_trop_Orthoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd_trop_Orthoptera$LUI=="Second_tropary vegetation") & (nd_trop_Orthoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd_trop_Orthoptera$LUI=="Agriculture_Low") & (nd_trop_Orthoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd_trop_Orthoptera$LUI=="Agriculture_Low") & (nd_trop_Orthoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd_trop_Orthoptera$LUI=="Agriculture_High") & (nd_trop_Orthoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd_trop_Orthoptera$LUI=="Agriculture_High") & (nd_trop_Orthoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot

nd_trop_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                         FUN = median,na.rm=TRUE))*100)-100
nd_trop_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_trop_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_trop_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                      FUN = median,na.rm=TRUE))*100)-100
nd_trop_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_trop_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_trop_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                        FUN = median,na.rm=TRUE))*100)-100
nd_trop_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                       FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_trop_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                       FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_trop_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd_trop_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_trop_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_trop_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd_trop_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_trop_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_trop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                         FUN = median,na.rm=TRUE))*100)-100
nd_trop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_trop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
nd_trop_Coleoptera$LUI <- factor(nd_trop_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Diptera$LUI <- factor(nd_trop_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Hemiptera$LUI <- factor(nd_trop_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Hymenoptera$LUI <- factor(nd_trop_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Lepidoptera$LUI <- factor(nd_trop_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Orthoptera$LUI <- factor(nd_trop_Orthoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot

p_coleoptera <- ggplot(data = nd_trop_Coleoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_trop_Coleoptera$PredLower, ymax = nd_trop_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_diptera <- ggplot(data = nd_trop_Diptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_trop_Diptera$PredLower, ymax = nd_trop_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_hemiptera <- ggplot(data = nd_trop_Hemiptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_trop_Hemiptera$PredLower, ymax = nd_trop_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_hymenoptera <- ggplot(data = nd_trop_Hymenoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_trop_Hymenoptera$PredLower, ymax = nd_trop_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_lepidoptera <- ggplot(data = nd_trop_Lepidoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_trop_Lepidoptera$PredLower, ymax = nd_trop_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_orthoptera <- ggplot(data = nd_trop_Orthoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_trop_Orthoptera$PredLower, ymax = nd_trop_Orthoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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
  ggtitle("Orthoptera")

# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(), 
          legend.text = element_text(size = 6), 
          legend.title = element_blank())
)


# put them all together to save them
MeanAnomAbund_trop <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,p_orthoptera)
#MeanAnomAbund_trop <- cowplot::plot_grid(MeanAnomAbund_trop,legend,ncol=1, rel_heights = c(1,0.1))

# save the ggplots
#ggsave(filename = paste0(outDir, "MeanAnomAbund_trop.pdf"), plot = MeanAnomAbund_trop, width = 200, height = 150, units = "mm", dpi = 300)
#ggsave(filename = paste0(outDir, "MeanAnomAbund_trop_extended yaxis.pdf"), plot = MeanAnomAbund_trop, width = 200, height = 150, units = "mm", dpi = 300)

## Abundance, Mean Anomaly ##
## NonTropical ##

nd_nontrop <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund_nontrop$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund_nontrop$data$StdTmeanAnomalyRS),
                        length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund_nontrop$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")))

# back transform the predictors
nd_nontrop$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd_nontrop$StdTmeanAnomalyRS,
  originalX = nontrop$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd_nontrop$LogAbund <- 0
nd_nontrop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# Record refRow for later use (see 'Values')
refRow <- which((nd_nontrop$LUI=="Primary vegetation") & (nd_nontrop$StdTmeanAnomaly==min(abs(nd_nontrop$StdTmeanAnomaly))))
# first row

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelAbund_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_nontrop$data$LUI=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelAbund_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_nontrop$data$LUI=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelAbund_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_nontrop$data$LUI=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelAbund_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_nontrop$data$LUI=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean.nontrop <- PredictGLMERRandIter(model = MeanAnomalyModelAbund_nontrop$model,data = nd_nontrop)

# back transform the abundance values
a.preds.tmean.nontrop <- exp(a.preds.tmean.nontrop)-0.01


# split by order into matrices, then name them
number_of_chunks = 6
list_a.preds.tmean.nontrop <- lapply(seq(1, NROW(a.preds.tmean.nontrop), ceiling(NROW(a.preds.tmean.nontrop)/number_of_chunks)),
                                     function(i) a.preds.tmean.nontrop[i:min(i + ceiling(NROW(a.preds.tmean.nontrop)/number_of_chunks) - 1, NROW(a.preds.tmean.nontrop)),])

names(list_a.preds.tmean.nontrop) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")
list2env(list_a.preds.tmean.nontrop,globalenv())

# sweep out refRow
list_a.preds.tmean.nontrop <- lapply(list_a.preds.tmean.nontrop,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[1,],FUN="/") 
})

list2env(list_a.preds.tmean.nontrop,globalenv())

# split nd_nontrop by order
Order<- paste0("nd_nontrop_",nd_nontrop$Order)
# create a list of data frames
by_Order <- split(nd_nontrop,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles

Coleoptera[which((nd_nontrop_Coleoptera$LUI=="Primary vegetation") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd_nontrop_Coleoptera$LUI=="Primary vegetation") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd_nontrop_Coleoptera$LUI=="Second_nontropary vegetation") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd_nontrop_Coleoptera$LUI=="Second_nontropary vegetation") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd_nontrop_Coleoptera$LUI=="Agriculture_Low") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd_nontrop_Coleoptera$LUI=="Agriculture_Low") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd_nontrop_Coleoptera$LUI=="Agriculture_High") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd_nontrop_Coleoptera$LUI=="Agriculture_High") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd_nontrop_Diptera$LUI=="Primary vegetation") & (nd_nontrop_Diptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd_nontrop_Diptera$LUI=="Primary vegetation") & (nd_nontrop_Diptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd_nontrop_Diptera$LUI=="Second_nontropary vegetation") & (nd_nontrop_Diptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd_nontrop_Diptera$LUI=="Second_nontropary vegetation") & (nd_nontrop_Diptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd_nontrop_Diptera$LUI=="Agriculture_Low") & (nd_nontrop_Diptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd_nontrop_Diptera$LUI=="Agriculture_Low") & (nd_nontrop_Diptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd_nontrop_Diptera$LUI=="Agriculture_High") & (nd_nontrop_Diptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd_nontrop_Diptera$LUI=="Agriculture_High") & (nd_nontrop_Diptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd_nontrop_Hemiptera$LUI=="Primary vegetation") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd_nontrop_Hemiptera$LUI=="Primary vegetation") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd_nontrop_Hemiptera$LUI=="Second_nontropary vegetation") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd_nontrop_Hemiptera$LUI=="Second_nontropary vegetation") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd_nontrop_Hemiptera$LUI=="Agriculture_Low") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd_nontrop_Hemiptera$LUI=="Agriculture_Low") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd_nontrop_Hemiptera$LUI=="Agriculture_High") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd_nontrop_Hemiptera$LUI=="Agriculture_High") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd_nontrop_Hymenoptera$LUI=="Primary vegetation") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd_nontrop_Hymenoptera$LUI=="Primary vegetation") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd_nontrop_Hymenoptera$LUI=="Second_nontropary vegetation") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd_nontrop_Hymenoptera$LUI=="Second_nontropary vegetation") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd_nontrop_Hymenoptera$LUI=="Agriculture_Low") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd_nontrop_Hymenoptera$LUI=="Agriculture_Low") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd_nontrop_Hymenoptera$LUI=="Agriculture_High") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd_nontrop_Hymenoptera$LUI=="Agriculture_High") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd_nontrop_Lepidoptera$LUI=="Primary vegetation") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd_nontrop_Lepidoptera$LUI=="Primary vegetation") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd_nontrop_Lepidoptera$LUI=="Second_nontropary vegetation") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd_nontrop_Lepidoptera$LUI=="Second_nontropary vegetation") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd_nontrop_Lepidoptera$LUI=="Agriculture_Low") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd_nontrop_Lepidoptera$LUI=="Agriculture_Low") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd_nontrop_Lepidoptera$LUI=="Agriculture_High") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd_nontrop_Lepidoptera$LUI=="Agriculture_High") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd_nontrop_Orthoptera$LUI=="Primary vegetation") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd_nontrop_Orthoptera$LUI=="Primary vegetation") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd_nontrop_Orthoptera$LUI=="Second_nontropary vegetation") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd_nontrop_Orthoptera$LUI=="Second_nontropary vegetation") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd_nontrop_Orthoptera$LUI=="Agriculture_Low") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd_nontrop_Orthoptera$LUI=="Agriculture_Low") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd_nontrop_Orthoptera$LUI=="Agriculture_High") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd_nontrop_Orthoptera$LUI=="Agriculture_High") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd_nontrop_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                            FUN = median,na.rm=TRUE))*100)-100
nd_nontrop_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_nontrop_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_nontrop_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                         FUN = median,na.rm=TRUE))*100)-100
nd_nontrop_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_nontrop_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_nontrop_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                           FUN = median,na.rm=TRUE))*100)-100
nd_nontrop_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_nontrop_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_nontrop_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                             FUN = median,na.rm=TRUE))*100)-100
nd_nontrop_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_nontrop_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_nontrop_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                             FUN = median,na.rm=TRUE))*100)-100
nd_nontrop_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_nontrop_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_nontrop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                            FUN = median,na.rm=TRUE))*100)-100
nd_nontrop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_nontrop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# set factor levels
# nd$LUI <- factor(nd$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Coleoptera$LUI <- factor(nd_nontrop_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Diptera$LUI <- factor(nd_nontrop_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Hemiptera$LUI <- factor(nd_nontrop_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Hymenoptera$LUI <- factor(nd_nontrop_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Lepidoptera$LUI <- factor(nd_nontrop_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Orthoptera$LUI <- factor(nd_nontrop_Orthoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot

p_coleoptera <- ggplot(data = nd_nontrop_Coleoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_nontrop_Coleoptera$PredLower, ymax = nd_nontrop_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700), limits = c(-100, 700)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) + 
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

p_diptera <- ggplot(data = nd_nontrop_Diptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_nontrop_Diptera$PredLower, ymax = nd_nontrop_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700), limits = c(-100, 700)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
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

p_hemiptera <- ggplot(data = nd_nontrop_Hemiptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_nontrop_Hemiptera$PredLower, ymax = nd_nontrop_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700), limits = c(-100, 700)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
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

p_hymenoptera <- ggplot(data = nd_nontrop_Hymenoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_nontrop_Hymenoptera$PredLower, ymax = nd_nontrop_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700), limits = c(-100, 700)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
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

p_lepidoptera <- ggplot(data = nd_nontrop_Lepidoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_nontrop_Lepidoptera$PredLower, ymax = nd_nontrop_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700), limits = c(-100, 700)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
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

p_orthoptera <- ggplot(data = nd_nontrop_Orthoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd_nontrop_Orthoptera$PredLower, ymax = nd_nontrop_Orthoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700), limits = c(-100, 700)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
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
  ggtitle("Orthoptera")


# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(), 
          legend.text = element_text(size = 6), 
          legend.title = element_blank())
)


# put them all together to save them
MeanAnomAbund_nontrop <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,p_orthoptera)
#MeanAnomAbund_nontrop <- cowplot::plot_grid(MeanAnomAbund_nontrop,legend,ncol=1, rel_heights = c(1,0.1))

# save the ggplot
#ggsave(filename = paste0(outDir, "MeanAnomAbund_nontrop.pdf"), plot = MeanAnomAbund_nontrop, width = 200, height = 150, units = "mm", dpi = 300)
#ggsave(filename = paste0(outDir, "MeanAnomAbund_nontrop_extended yaxis.pdf"), plot = MeanAnomAbund_nontrop, width = 200, height = 150, units = "mm", dpi = 300)

# plot Tropical and NonTropical plots together
# to do this, comment the lines wherein cowplot is used to add a legend to each plot and where ggsave is used to save individual plots (405,408,409,779,782,783)
# re-run from the start, to this point
# comment lines according to which style of plot is desired (stacked or side-by-side)

#MeanAnomAbundRealms <-cowplot::plot_grid(MeanAnomAbund_nontrop, NULL, MeanAnomAbund_trop, ncol=1, labels = c('A',"",'B'),rel_heights = c(1,0.1,1)) # stacked
MeanAnomAbundRealms <-cowplot::plot_grid(MeanAnomAbund_nontrop, NULL, MeanAnomAbund_trop, ncol=3, labels = c('A',"",'B'),rel_widths = c(1,0.1,1)) # side-by-side

MeanAnomAbundRealms <-cowplot::plot_grid(MeanAnomAbundRealms, legend, ncol=1, rel_heights = c(1,0.1))

#ggsave(filename = paste0(outDir, "MeanAnomAbundRealms_stacked.pdf"), plot = MeanAnomAbundRealms, width = 200, height = 300, units = "mm", dpi = 300)
ggsave(filename = paste0(outDir, "MeanAnomAbundRealms_sidebyside.pdf"), plot = MeanAnomAbundRealms, width = 400, height = 150, units = "mm", dpi = 300)

## Richness, Mean Anomaly ##
## Tropical ##

nd2_trop <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelRich_trop$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelRich_trop$data$StdTmeanAnomalyRS),
                        length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelRich_trop$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")))

# back transform the predictors
nd2_trop$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2_trop$StdTmeanAnomalyRS,
  originalX = trop$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd2_trop$LogAbund <- 0
nd2_trop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# Record refRow for later use
refRow <- which((nd2_trop$LUI=="Primary vegetation") & (nd2_trop$StdTmeanAnomaly==min(abs(nd2_trop$StdTmeanAnomaly))))
# the first row

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelRich_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_trop$data$LUI=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelRich_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_trop$data$LUI=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelRich_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_trop$data$LUI=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelRich_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_trop$data$LUI=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
sr.preds.tmean.trop <- PredictGLMERRandIter(model = MeanAnomalyModelRich_trop$model,data = nd2_trop)

# back transform the abundance values
sr.preds.tmean.trop <- exp(sr.preds.tmean.trop)-0.01

# split by order into matrices, then name them
number_of_chunks = 6
list_sr.preds.tmean.trop <- lapply(seq(1, NROW(sr.preds.tmean.trop), ceiling(NROW(sr.preds.tmean.trop)/number_of_chunks)),
                                   function(i) sr.preds.tmean.trop[i:min(i + ceiling(NROW(sr.preds.tmean.trop)/number_of_chunks) - 1, NROW(sr.preds.tmean.trop)),])
names(list_sr.preds.tmean.trop) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")
list2env(list_sr.preds.tmean.trop,globalenv())


# sweep out refRow
list_sr.preds.tmean.trop <- lapply(list_sr.preds.tmean.trop,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[1,],FUN="/") 
})

list2env(list_sr.preds.tmean.trop,globalenv())

# split nd2_trop by order
Order<- paste0("nd2_trop_",nd2_trop$Order)
# create a list of data frames
by_Order <- split(nd2_trop,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles

Coleoptera[which((nd2_trop_Coleoptera$LUI=="Primary vegetation") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd2_trop_Coleoptera$LUI=="Primary vegetation") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd2_trop_Coleoptera$LUI=="Second2_tropary vegetation") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd2_trop_Coleoptera$LUI=="Second2_tropary vegetation") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd2_trop_Coleoptera$LUI=="Agriculture_Low") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd2_trop_Coleoptera$LUI=="Agriculture_Low") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd2_trop_Coleoptera$LUI=="Agriculture_High") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd2_trop_Coleoptera$LUI=="Agriculture_High") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd2_trop_Diptera$LUI=="Primary vegetation") & (nd2_trop_Diptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd2_trop_Diptera$LUI=="Primary vegetation") & (nd2_trop_Diptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd2_trop_Diptera$LUI=="Second2_tropary vegetation") & (nd2_trop_Diptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd2_trop_Diptera$LUI=="Second2_tropary vegetation") & (nd2_trop_Diptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd2_trop_Diptera$LUI=="Agriculture_Low") & (nd2_trop_Diptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd2_trop_Diptera$LUI=="Agriculture_Low") & (nd2_trop_Diptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd2_trop_Diptera$LUI=="Agriculture_High") & (nd2_trop_Diptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd2_trop_Diptera$LUI=="Agriculture_High") & (nd2_trop_Diptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd2_trop_Hemiptera$LUI=="Primary vegetation") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd2_trop_Hemiptera$LUI=="Primary vegetation") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd2_trop_Hemiptera$LUI=="Second2_tropary vegetation") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd2_trop_Hemiptera$LUI=="Second2_tropary vegetation") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd2_trop_Hemiptera$LUI=="Agriculture_Low") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd2_trop_Hemiptera$LUI=="Agriculture_Low") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd2_trop_Hemiptera$LUI=="Agriculture_High") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd2_trop_Hemiptera$LUI=="Agriculture_High") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd2_trop_Hymenoptera$LUI=="Primary vegetation") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd2_trop_Hymenoptera$LUI=="Primary vegetation") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd2_trop_Hymenoptera$LUI=="Second2_tropary vegetation") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd2_trop_Hymenoptera$LUI=="Second2_tropary vegetation") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd2_trop_Hymenoptera$LUI=="Agriculture_Low") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd2_trop_Hymenoptera$LUI=="Agriculture_Low") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd2_trop_Hymenoptera$LUI=="Agriculture_High") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd2_trop_Hymenoptera$LUI=="Agriculture_High") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd2_trop_Lepidoptera$LUI=="Primary vegetation") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd2_trop_Lepidoptera$LUI=="Primary vegetation") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd2_trop_Lepidoptera$LUI=="Second2_tropary vegetation") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd2_trop_Lepidoptera$LUI=="Second2_tropary vegetation") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd2_trop_Lepidoptera$LUI=="Agriculture_Low") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd2_trop_Lepidoptera$LUI=="Agriculture_Low") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd2_trop_Lepidoptera$LUI=="Agriculture_High") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd2_trop_Lepidoptera$LUI=="Agriculture_High") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd2_trop_Orthoptera$LUI=="Primary vegetation") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd2_trop_Orthoptera$LUI=="Primary vegetation") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd2_trop_Orthoptera$LUI=="Second2_tropary vegetation") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd2_trop_Orthoptera$LUI=="Second2_tropary vegetation") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd2_trop_Orthoptera$LUI=="Agriculture_Low") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd2_trop_Orthoptera$LUI=="Agriculture_Low") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd2_trop_Orthoptera$LUI=="Agriculture_High") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd2_trop_Orthoptera$LUI=="Agriculture_High") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot

nd2_trop_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd2_trop_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_trop_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_trop_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                       FUN = median,na.rm=TRUE))*100)-100
nd2_trop_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                      FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_trop_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                      FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_trop_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                         FUN = median,na.rm=TRUE))*100)-100
nd2_trop_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_trop_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_trop_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                           FUN = median,na.rm=TRUE))*100)-100
nd2_trop_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_trop_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_trop_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                           FUN = median,na.rm=TRUE))*100)-100
nd2_trop_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_trop_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_trop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd2_trop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_trop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
nd2_trop_Coleoptera$LUI <- factor(nd2_trop_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Diptera$LUI <- factor(nd2_trop_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Hemiptera$LUI <- factor(nd2_trop_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Hymenoptera$LUI <- factor(nd2_trop_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Lepidoptera$LUI <- factor(nd2_trop_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Orthoptera$LUI <- factor(nd2_trop_Orthoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot 

p_coleoptera <- ggplot(data = nd2_trop_Coleoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_trop_Coleoptera$PredLower, ymax = nd2_trop_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_diptera <- ggplot(data = nd2_trop_Diptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_trop_Diptera$PredLower, ymax = nd2_trop_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_hemiptera <- ggplot(data = nd2_trop_Hemiptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_trop_Hemiptera$PredLower, ymax = nd2_trop_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_hymenoptera <- ggplot(data = nd2_trop_Hymenoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_trop_Hymenoptera$PredLower, ymax = nd2_trop_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_lepidoptera <- ggplot(data = nd2_trop_Lepidoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_trop_Lepidoptera$PredLower, ymax = nd2_trop_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_orthoptera <- ggplot(data = nd2_trop_Orthoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_trop_Orthoptera$PredLower, ymax = nd2_trop_Orthoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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
  ggtitle("Orthoptera")


# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(), 
          legend.text = element_text(size = 6), 
          legend.title = element_blank())
)


# put them all together to save them
MeanAnomRich_trop <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,p_orthoptera)
#MeanAnomRich_trop <- cowplot::plot_grid(MeanAnomRich_trop,legend,ncol=1, rel_heights = c(1,0.1))

# save the ggplots
#ggsave(filename = paste0(outDir, "MeanAnomRich_trop.pdf"), plot = MeanAnomRich_trop, width = 200, height = 150, units = "mm", dpi = 300)
#ggsave(filename = paste0(outDir, "MeanAnomRich_trop_extended yaxis.pdf"), plot = MeanAnomRich_trop, width = 200, height = 150, units = "mm", dpi = 300)


## Richness, Mean Anomaly ##
## NonTropical ##

nd2_nontrop <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelRich_nontrop$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelRich_nontrop$data$StdTmeanAnomalyRS),
                        length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelRich_nontrop$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")))

# back transform the predictors
nd2_nontrop$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2_nontrop$StdTmeanAnomalyRS,
  originalX = nontrop$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd2_nontrop$LogAbund <- 0
nd2_nontrop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# does this for each Order
# Record refRow for later use (see 'Values')
refRow <- which((nd2_nontrop$LUI=="Primary vegetation") & (nd2_nontrop$StdTmeanAnomaly==min(abs(nd2_nontrop$StdTmeanAnomaly))))
# the first row, every 400 rows

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelRich_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_nontrop$data$LUI=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelRich_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_nontrop$data$LUI=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelRich_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_nontrop$data$LUI=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelRich_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_nontrop$data$LUI=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
sr.preds.tmean.nontrop <- PredictGLMERRandIter(model = MeanAnomalyModelRich_nontrop$model,data = nd2_nontrop)

# back transform the abundance values
sr.preds.tmean.nontrop <- exp(sr.preds.tmean.nontrop)-0.01


# split by order into matrices, then name them
number_of_chunks = 6
list_sr.preds.tmean.nontrop <- lapply(seq(1, NROW(sr.preds.tmean.nontrop), ceiling(NROW(sr.preds.tmean.nontrop)/number_of_chunks)),
                                      function(i) sr.preds.tmean.nontrop[i:min(i + ceiling(NROW(sr.preds.tmean.nontrop)/number_of_chunks) - 1, NROW(sr.preds.tmean.nontrop)),])

names(list_sr.preds.tmean.nontrop) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")
list2env(list_sr.preds.tmean.nontrop,globalenv())

# sweep out refRow
list_sr.preds.tmean.nontrop <- lapply(list_sr.preds.tmean.nontrop,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[1,],FUN="/") 
})

list2env(list_sr.preds.tmean.nontrop,globalenv())

# split nd2_nontrop by order
Order<- paste0("nd2_nontrop_",nd2_nontrop$Order)
# create a list of data frames
by_Order <- split(nd2_nontrop,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles

Coleoptera[which((nd2_nontrop_Coleoptera$LUI=="Primary vegetation") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd2_nontrop_Coleoptera$LUI=="Primary vegetation") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd2_nontrop_Coleoptera$LUI=="Secondary vegetation") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd2_nontrop_Coleoptera$LUI=="Secondary vegetation") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd2_nontrop_Coleoptera$LUI=="Agriculture_Low") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd2_nontrop_Coleoptera$LUI=="Agriculture_Low") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd2_nontrop_Coleoptera$LUI=="Agriculture_High") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd2_nontrop_Coleoptera$LUI=="Agriculture_High") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd2_nontrop_Diptera$LUI=="Primary vegetation") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd2_nontrop_Diptera$LUI=="Primary vegetation") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd2_nontrop_Diptera$LUI=="Secondary vegetation") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd2_nontrop_Diptera$LUI=="Secondary vegetation") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd2_nontrop_Diptera$LUI=="Agriculture_Low") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd2_nontrop_Diptera$LUI=="Agriculture_Low") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd2_nontrop_Diptera$LUI=="Agriculture_High") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd2_nontrop_Diptera$LUI=="Agriculture_High") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd2_nontrop_Hemiptera$LUI=="Primary vegetation") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd2_nontrop_Hemiptera$LUI=="Primary vegetation") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd2_nontrop_Hemiptera$LUI=="Secondary vegetation") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd2_nontrop_Hemiptera$LUI=="Secondary vegetation") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd2_nontrop_Hemiptera$LUI=="Agriculture_Low") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd2_nontrop_Hemiptera$LUI=="Agriculture_Low") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd2_nontrop_Hemiptera$LUI=="Agriculture_High") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd2_nontrop_Hemiptera$LUI=="Agriculture_High") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd2_nontrop_Hymenoptera$LUI=="Primary vegetation") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd2_nontrop_Hymenoptera$LUI=="Primary vegetation") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd2_nontrop_Hymenoptera$LUI=="Secondary vegetation") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd2_nontrop_Hymenoptera$LUI=="Secondary vegetation") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd2_nontrop_Hymenoptera$LUI=="Agriculture_Low") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd2_nontrop_Hymenoptera$LUI=="Agriculture_Low") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd2_nontrop_Hymenoptera$LUI=="Agriculture_High") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd2_nontrop_Hymenoptera$LUI=="Agriculture_High") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd2_nontrop_Lepidoptera$LUI=="Primary vegetation") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd2_nontrop_Lepidoptera$LUI=="Primary vegetation") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd2_nontrop_Lepidoptera$LUI=="Secondary vegetation") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd2_nontrop_Lepidoptera$LUI=="Secondary vegetation") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd2_nontrop_Lepidoptera$LUI=="Agriculture_Low") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd2_nontrop_Lepidoptera$LUI=="Agriculture_Low") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd2_nontrop_Lepidoptera$LUI=="Agriculture_High") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd2_nontrop_Lepidoptera$LUI=="Agriculture_High") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd2_nontrop_Orthoptera$LUI=="Primary vegetation") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd2_nontrop_Orthoptera$LUI=="Primary vegetation") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd2_nontrop_Orthoptera$LUI=="Secondary vegetation") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd2_nontrop_Orthoptera$LUI=="Secondary vegetation") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd2_nontrop_Orthoptera$LUI=="Agriculture_Low") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd2_nontrop_Orthoptera$LUI=="Agriculture_Low") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd2_nontrop_Orthoptera$LUI=="Agriculture_High") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd2_nontrop_Orthoptera$LUI=="Agriculture_High") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd2_nontrop_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                             FUN = median,na.rm=TRUE))*100)-100
nd2_nontrop_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_nontrop_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_nontrop_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd2_nontrop_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_nontrop_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_nontrop_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                            FUN = median,na.rm=TRUE))*100)-100
nd2_nontrop_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_nontrop_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_nontrop_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                              FUN = median,na.rm=TRUE))*100)-100
nd2_nontrop_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_nontrop_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_nontrop_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                              FUN = median,na.rm=TRUE))*100)-100
nd2_nontrop_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_nontrop_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_nontrop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                             FUN = median,na.rm=TRUE))*100)-100
nd2_nontrop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_nontrop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
nd2_nontrop_Coleoptera$LUI <- factor(nd2_nontrop_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Diptera$LUI <- factor(nd2_nontrop_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Hemiptera$LUI <- factor(nd2_nontrop_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Hymenoptera$LUI <- factor(nd2_nontrop_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Lepidoptera$LUI <- factor(nd2_nontrop_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Orthoptera$LUI <- factor(nd2_nontrop_Orthoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot

p_coleoptera <- ggplot(data = nd2_nontrop_Coleoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_nontrop_Coleoptera$PredLower, ymax = nd2_nontrop_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500), limits = c(-100, 500)) +
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

p_diptera <- ggplot(data = nd2_nontrop_Diptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_nontrop_Diptera$PredLower, ymax = nd2_nontrop_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500), limits = c(-100, 500)) +
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

p_hemiptera <- ggplot(data = nd2_nontrop_Hemiptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_nontrop_Hemiptera$PredLower, ymax = nd2_nontrop_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500), limits = c(-100, 500)) +
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

p_hymenoptera <- ggplot(data = nd2_nontrop_Hymenoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_nontrop_Hymenoptera$PredLower, ymax = nd2_nontrop_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500), limits = c(-100, 500)) +
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

p_lepidoptera <- ggplot(data = nd2_nontrop_Lepidoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_nontrop_Lepidoptera$PredLower, ymax = nd2_nontrop_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500), limits = c(-100, 500)) +
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

p_orthoptera <- ggplot(data = nd2_nontrop_Orthoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_nontrop_Orthoptera$PredLower, ymax = nd2_nontrop_Orthoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.25, 0.5, 0.75, 1), limits = c(0, 1)) +
  #scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000), limits = c(-100, 1000)) +
  scale_y_continuous(breaks = c(-100, 0, 100, 200, 300, 400, 500), limits = c(-100, 500)) +
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
  ggtitle("Orthoptera")

# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(), 
          legend.text = element_text(size = 6), 
          legend.title = element_blank())
)


# put them all together to save them
MeanAnomRich_nontrop <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,p_orthoptera)
#MeanAnomRich_nontrop <- cowplot::plot_grid(MeanAnomRich_nontrop,legend,ncol=1, rel_heights = c(1,0.1))

# save the ggplots
#ggsave(filename = paste0(outDir, "MeanAnomRich_nontrop.pdf"), plot = MeanAnomRich_nontrop, width = 200, height = 150, units = "mm", dpi = 300)
#ggsave(filename = paste0(outDir, "MeanAnomRich_nontrop_extended yaxis.pdf"), plot = MeanAnomRich_nontrop, width = 200, height = 150, units = "mm", dpi = 300)

# plot Tropical and NonTropical plots together
# to do this, comment the lines wherein cowplot is used to add a legend to each plot and where ggsave is used to save individual plots (405,408,409,1539,1542,1543 )
# re-run from the start, to this point
# comment lines according to which style of plot is desired (stacked or side-by-side)

#MeanAnomRichRealms <-cowplot::plot_grid(MeanAnomRich_nontrop, NULL, MeanAnomRich_trop, ncol=1, labels = c('A',"",'B'),rel_heights = c(1,0.1,1)) # stacked
MeanAnomRichRealms <-cowplot::plot_grid(MeanAnomRich_nontrop, NULL, MeanAnomRich_trop, ncol=3, labels = c('A',"",'B'),rel_widths = c(1,0.1,1)) # side-by-side

MeanAnomRichRealms <-cowplot::plot_grid(MeanAnomRichRealms, legend, ncol=1, rel_heights = c(1,0.1))

#ggsave(filename = paste0(outDir, "MeanAnomRichRealms_stacked.pdf"), plot = MeanAnomRichRealms, width = 200, height = 300, units = "mm", dpi = 300)
ggsave(filename = paste0(outDir, "MeanAnomRichRealms_sidebyside.pdf"), plot = MeanAnomRichRealms, width = 400, height = 150, units = "mm", dpi = 300)


## Abundance, Max Anomaly ##
## Tropical ##

nd3_trop <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelAbund_trop$data$StdTmaxAnomalyRS),
                       to = max(MaxAnomalyModelAbund_trop$data$StdTmaxAnomalyRS),
                       length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelAbund_trop$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")))

# back transform the predictors
nd3_trop$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3_trop$StdTmaxAnomalyRS,
  originalX = trop$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd3_trop$LogAbund <- 0
nd3_trop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# Record for later use (see 'Values' or check the data set itself for the row in Coleoptera that satisfies the conditions)
refRow <- which((nd3_trop$LUI=="Primary vegetation") & (nd3_trop$StdTmaxAnomaly==min(abs(nd3_trop$StdTmaxAnomaly))))
# 56th row

QPV <- quantile(x = MaxAnomalyModelAbund_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_trop$data$LUI=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAnomalyModelAbund_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_trop$data$LUI=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAnomalyModelAbund_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_trop$data$LUI=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAnomalyModelAbund_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_trop$data$LUI=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmax.trop <- PredictGLMERRandIter(model = MaxAnomalyModelAbund_trop$model,data = nd3_trop)

# back transform the abundance values
a.preds.tmax.trop <- exp(a.preds.tmax.trop)-0.01


# split by order into matrices, then name them
number_of_chunks = 6
list_a.preds.tmax.trop <- lapply(seq(1, NROW(a.preds.tmax.trop), ceiling(NROW(a.preds.tmax.trop)/number_of_chunks)),
                                 function(i) a.preds.tmax.trop[i:min(i + ceiling(NROW(a.preds.tmax.trop)/number_of_chunks) - 1, NROW(a.preds.tmax.trop)),])
# success!
# creates list of matrices
# name them
names(list_a.preds.tmax.trop) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")
list2env(list_a.preds.tmax.trop,globalenv())


# sweep out refRow
list_a.preds.tmax.trop <- lapply(list_a.preds.tmax.trop,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[56,],FUN="/") 
})

list2env(list_a.preds.tmax.trop,globalenv())

# split nd3_trop by order
Order<- paste0("nd3_trop_",nd3_trop$Order)
# create a list of data frames
by_Order <- split(nd3_trop,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles

Coleoptera[which((nd3_trop_Coleoptera$LUI=="Primary vegetation") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd3_trop_Coleoptera$LUI=="Primary vegetation") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd3_trop_Coleoptera$LUI=="Secondary vegetation") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd3_trop_Coleoptera$LUI=="Secondary vegetation") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd3_trop_Coleoptera$LUI=="Agriculture_Low") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd3_trop_Coleoptera$LUI=="Agriculture_Low") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd3_trop_Coleoptera$LUI=="Agriculture_High") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd3_trop_Coleoptera$LUI=="Agriculture_High") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd3_trop_Diptera$LUI=="Primary vegetation") & (nd3_trop_Diptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd3_trop_Diptera$LUI=="Primary vegetation") & (nd3_trop_Diptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd3_trop_Diptera$LUI=="Secondary vegetation") & (nd3_trop_Diptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd3_trop_Diptera$LUI=="Secondary vegetation") & (nd3_trop_Diptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd3_trop_Diptera$LUI=="Agriculture_Low") & (nd3_trop_Diptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd3_trop_Diptera$LUI=="Agriculture_Low") & (nd3_trop_Diptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd3_trop_Diptera$LUI=="Agriculture_High") & (nd3_trop_Diptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd3_trop_Diptera$LUI=="Agriculture_High") & (nd3_trop_Diptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd3_trop_Hemiptera$LUI=="Primary vegetation") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd3_trop_Hemiptera$LUI=="Primary vegetation") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd3_trop_Hemiptera$LUI=="Secondary vegetation") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd3_trop_Hemiptera$LUI=="Secondary vegetation") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd3_trop_Hemiptera$LUI=="Agriculture_Low") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd3_trop_Hemiptera$LUI=="Agriculture_Low") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd3_trop_Hemiptera$LUI=="Agriculture_High") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd3_trop_Hemiptera$LUI=="Agriculture_High") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd3_trop_Hymenoptera$LUI=="Primary vegetation") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd3_trop_Hymenoptera$LUI=="Primary vegetation") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd3_trop_Hymenoptera$LUI=="Secondary vegetation") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd3_trop_Hymenoptera$LUI=="Secondary vegetation") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd3_trop_Hymenoptera$LUI=="Agriculture_Low") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd3_trop_Hymenoptera$LUI=="Agriculture_Low") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd3_trop_Hymenoptera$LUI=="Agriculture_High") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd3_trop_Hymenoptera$LUI=="Agriculture_High") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd3_trop_Lepidoptera$LUI=="Primary vegetation") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd3_trop_Lepidoptera$LUI=="Primary vegetation") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd3_trop_Lepidoptera$LUI=="Secondary vegetation") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd3_trop_Lepidoptera$LUI=="Secondary vegetation") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd3_trop_Lepidoptera$LUI=="Agriculture_Low") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd3_trop_Lepidoptera$LUI=="Agriculture_Low") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd3_trop_Lepidoptera$LUI=="Agriculture_High") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd3_trop_Lepidoptera$LUI=="Agriculture_High") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd3_trop_Orthoptera$LUI=="Primary vegetation") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd3_trop_Orthoptera$LUI=="Primary vegetation") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd3_trop_Orthoptera$LUI=="Secondary vegetation") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd3_trop_Orthoptera$LUI=="Secondary vegetation") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd3_trop_Orthoptera$LUI=="Agriculture_Low") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd3_trop_Orthoptera$LUI=="Agriculture_Low") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd3_trop_Orthoptera$LUI=="Agriculture_High") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd3_trop_Orthoptera$LUI=="Agriculture_High") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd3_trop_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd3_trop_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_trop_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_trop_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                       FUN = median,na.rm=TRUE))*100)-100
nd3_trop_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                      FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_trop_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                      FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_trop_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                         FUN = median,na.rm=TRUE))*100)-100
nd3_trop_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_trop_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_trop_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                           FUN = median,na.rm=TRUE))*100)-100
nd3_trop_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_trop_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_trop_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                           FUN = median,na.rm=TRUE))*100)-100
nd3_trop_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_trop_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_trop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd3_trop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_trop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
nd3_trop_Coleoptera$LUI <- factor(nd3_trop_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Diptera$LUI <- factor(nd3_trop_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Hemiptera$LUI <- factor(nd3_trop_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Hymenoptera$LUI <- factor(nd3_trop_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Lepidoptera$LUI <- factor(nd3_trop_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Orthoptera$LUI <- factor(nd3_trop_Orthoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot

p_coleoptera <- ggplot(data = nd3_trop_Coleoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_trop_Coleoptera$PredLower, ymax = nd3_trop_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_diptera <- ggplot(data = nd3_trop_Diptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_trop_Diptera$PredLower, ymax = nd3_trop_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_hemiptera <- ggplot(data = nd3_trop_Hemiptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_trop_Hemiptera$PredLower, ymax = nd3_trop_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_hymenoptera <- ggplot(data = nd3_trop_Hymenoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_trop_Hymenoptera$PredLower, ymax = nd3_trop_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_lepidoptera <- ggplot(data = nd3_trop_Lepidoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_trop_Lepidoptera$PredLower, ymax = nd3_trop_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_orthoptera <- ggplot(data = nd3_trop_Orthoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_trop_Orthoptera$PredLower, ymax = nd3_trop_Orthoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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
  ggtitle("Orthoptera")

# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(), 
          legend.text = element_text(size = 6), 
          legend.title = element_blank())
)


# put them all together to save them
MaxAnomAbund_trop <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,p_orthoptera)
#MaxAnomAbund_trop <- cowplot::plot_grid(MaxAnomAbund_trop,legend,ncol=1,rel_heights = c(1,0.1))

# save them
#ggsave(filename = paste0(outDir, "MaxAnomAbund_trop.pdf"), plot = MaxAnomAbund_trop, width = 200, height = 150, units = "mm", dpi = 300)
#ggsave(filename = paste0(outDir, "MaxAnomAbund_trop_extended yaxis.pdf"), plot = MaxAnomAbund_trop, width = 200, height = 150, units = "mm", dpi = 300)

## Abundance, Max Anomaly ##
## NonTropical ##

nd3_nontrop <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelAbund_nontrop$data$StdTmaxAnomalyRS),
                       to = max(MaxAnomalyModelAbund_nontrop$data$StdTmaxAnomalyRS),
                       length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelAbund_nontrop$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")))

# back transform the predictors
nd3_nontrop$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3_nontrop$StdTmaxAnomalyRS,
  originalX = nontrop$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd3_nontrop$LogAbund <- 0
nd3_nontrop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# Record for later use (see 'Values' or check the data set itself for the row in Coleoptera that satisfies the conditions)
refRow <- which((nd3_nontrop$LUI=="Primary vegetation") & (nd3_nontrop$StdTmaxAnomaly==min(abs(nd3_nontrop$StdTmaxAnomaly))))
# 15th row, every 400 rows (had to check myself)

# adjust plot 1: max anomaly and abundance

QPV <- quantile(x = MaxAnomalyModelAbund_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_nontrop$data$LUI=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAnomalyModelAbund_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_nontrop$data$LUI=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAnomalyModelAbund_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_nontrop$data$LUI=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAnomalyModelAbund_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_nontrop$data$LUI=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmax.nontrop <- PredictGLMERRandIter(model = MaxAnomalyModelAbund_nontrop$model,data = nd3_nontrop)

# back transform the abundance values
a.preds.tmax.nontrop <- exp(a.preds.tmax.nontrop)-0.01


# split by order into matrices, then name them
number_of_chunks = 6
list_a.preds.tmax.nontrop <- lapply(seq(1, NROW(a.preds.tmax.nontrop), ceiling(NROW(a.preds.tmax.nontrop)/number_of_chunks)),
                                    function(i) a.preds.tmax.nontrop[i:min(i + ceiling(NROW(a.preds.tmax.nontrop)/number_of_chunks) - 1, NROW(a.preds.tmax.nontrop)),])

names(list_a.preds.tmax.nontrop) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")
list2env(list_a.preds.tmax.nontrop,globalenv())

# sweep out refRow
list_a.preds.tmax.nontrop <- lapply(list_a.preds.tmax.nontrop,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[15,],FUN="/") 
})

list2env(list_a.preds.tmax.nontrop,globalenv())

# split nd3_nontrop by order
Order<- paste0("nd3_nontrop_",nd3_nontrop$Order)
# create a list of data frames
by_Order <- split(nd3_nontrop,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles
Coleoptera[which((nd3_nontrop_Coleoptera$LUI=="Primary vegetation") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd3_nontrop_Coleoptera$LUI=="Primary vegetation") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd3_nontrop_Coleoptera$LUI=="Secondary vegetation") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd3_nontrop_Coleoptera$LUI=="Secondary vegetation") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd3_nontrop_Coleoptera$LUI=="Agriculture_Low") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd3_nontrop_Coleoptera$LUI=="Agriculture_Low") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd3_nontrop_Coleoptera$LUI=="Agriculture_High") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd3_nontrop_Coleoptera$LUI=="Agriculture_High") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd3_nontrop_Diptera$LUI=="Primary vegetation") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd3_nontrop_Diptera$LUI=="Primary vegetation") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd3_nontrop_Diptera$LUI=="Secondary vegetation") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd3_nontrop_Diptera$LUI=="Secondary vegetation") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd3_nontrop_Diptera$LUI=="Agriculture_Low") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd3_nontrop_Diptera$LUI=="Agriculture_Low") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd3_nontrop_Diptera$LUI=="Agriculture_High") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd3_nontrop_Diptera$LUI=="Agriculture_High") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd3_nontrop_Hemiptera$LUI=="Primary vegetation") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd3_nontrop_Hemiptera$LUI=="Primary vegetation") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd3_nontrop_Hemiptera$LUI=="Secondary vegetation") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd3_nontrop_Hemiptera$LUI=="Secondary vegetation") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd3_nontrop_Hemiptera$LUI=="Agriculture_Low") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd3_nontrop_Hemiptera$LUI=="Agriculture_Low") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd3_nontrop_Hemiptera$LUI=="Agriculture_High") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd3_nontrop_Hemiptera$LUI=="Agriculture_High") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd3_nontrop_Hymenoptera$LUI=="Primary vegetation") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd3_nontrop_Hymenoptera$LUI=="Primary vegetation") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd3_nontrop_Hymenoptera$LUI=="Secondary vegetation") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd3_nontrop_Hymenoptera$LUI=="Secondary vegetation") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd3_nontrop_Hymenoptera$LUI=="Agriculture_Low") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd3_nontrop_Hymenoptera$LUI=="Agriculture_Low") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd3_nontrop_Hymenoptera$LUI=="Agriculture_High") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd3_nontrop_Hymenoptera$LUI=="Agriculture_High") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd3_nontrop_Lepidoptera$LUI=="Primary vegetation") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd3_nontrop_Lepidoptera$LUI=="Primary vegetation") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd3_nontrop_Lepidoptera$LUI=="Secondary vegetation") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd3_nontrop_Lepidoptera$LUI=="Secondary vegetation") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd3_nontrop_Lepidoptera$LUI=="Agriculture_Low") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd3_nontrop_Lepidoptera$LUI=="Agriculture_Low") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd3_nontrop_Lepidoptera$LUI=="Agriculture_High") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd3_nontrop_Lepidoptera$LUI=="Agriculture_High") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd3_nontrop_Orthoptera$LUI=="Primary vegetation") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd3_nontrop_Orthoptera$LUI=="Primary vegetation") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd3_nontrop_Orthoptera$LUI=="Secondary vegetation") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd3_nontrop_Orthoptera$LUI=="Secondary vegetation") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd3_nontrop_Orthoptera$LUI=="Agriculture_Low") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd3_nontrop_Orthoptera$LUI=="Agriculture_Low") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd3_nontrop_Orthoptera$LUI=="Agriculture_High") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd3_nontrop_Orthoptera$LUI=="Agriculture_High") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd3_nontrop_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                             FUN = median,na.rm=TRUE))*100)-100
nd3_nontrop_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_nontrop_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_nontrop_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd3_nontrop_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_nontrop_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_nontrop_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                            FUN = median,na.rm=TRUE))*100)-100
nd3_nontrop_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_nontrop_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_nontrop_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                              FUN = median,na.rm=TRUE))*100)-100
nd3_nontrop_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_nontrop_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_nontrop_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                              FUN = median,na.rm=TRUE))*100)-100
nd3_nontrop_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_nontrop_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_nontrop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                             FUN = median,na.rm=TRUE))*100)-100
nd3_nontrop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_nontrop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
nd3_nontrop_Coleoptera$LUI <- factor(nd3_nontrop_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Diptera$LUI <- factor(nd3_nontrop_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Hemiptera$LUI <- factor(nd3_nontrop_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Hymenoptera$LUI <- factor(nd3_nontrop_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Lepidoptera$LUI <- factor(nd3_nontrop_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Orthoptera$LUI <- factor(nd3_nontrop_Orthoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot

p_coleoptera <- ggplot(data = nd3_nontrop_Coleoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_nontrop_Coleoptera$PredLower, ymax = nd3_nontrop_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_diptera <- ggplot(data = nd3_nontrop_Diptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_nontrop_Diptera$PredLower, ymax = nd3_nontrop_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_hemiptera <- ggplot(data = nd3_nontrop_Hemiptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_nontrop_Hemiptera$PredLower, ymax = nd3_nontrop_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_hymenoptera <- ggplot(data = nd3_nontrop_Hymenoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_nontrop_Hymenoptera$PredLower, ymax = nd3_nontrop_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_lepidoptera <- ggplot(data = nd3_nontrop_Lepidoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_nontrop_Lepidoptera$PredLower, ymax = nd3_nontrop_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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


p_orthoptera <- ggplot(data = nd3_nontrop_Orthoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_nontrop_Orthoptera$PredLower, ymax = nd3_nontrop_Orthoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
  #scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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
  ggtitle("Orthoptera")

# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(), 
          legend.text = element_text(size = 6), 
          legend.title = element_blank())
)


# put them all together to save them
MaxAnomAbund_nontrop <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,p_orthoptera)
#MaxAnomAbund_nontrop <- cowplot::plot_grid(MaxAnomAbund_nontrop,legend,ncol=1,rel_heights = c(1,0.1))

# save them
#ggsave(filename = paste0(outDir, "MaxAnomAbund_nontrop.pdf"), plot = MaxAnomAbund_nontrop, width = 200, height = 150, units = "mm", dpi = 300)
#ggsave(filename = paste0(outDir, "MaxAnomAbund_nontrop_extended yaxis.pdf"), plot = MaxAnomAbund_nontrop, width = 200, height = 150, units = "mm", dpi = 300)

# plot Tropical and NonTropical plots together
# to do this, comment the lines wherein cowplot is used to add a legend to each plot and where ggsave is used to save individual plots
# re-run from the start, to this point
# comment lines according to which style of plot is desired (stacked or side-by-side)

#MaxAnomAbundRealms <-cowplot::plot_grid(MaxAnomAbund_nontrop, NULL, MaxAnomAbund_trop, ncol=1, labels = c('A',"",'B'),rel_heights = c(1,0.1,1)) # stacked
MaxAnomAbundRealms <-cowplot::plot_grid(MaxAnomAbund_nontrop, NULL, MaxAnomAbund_trop, ncol=3, labels = c('A',"",'B'),rel_widths = c(1,0.1,1)) # side-by-side

MaxAnomAbundRealms <-cowplot::plot_grid(MaxAnomAbundRealms, legend, ncol=1, rel_heights = c(1,0.1))

#ggsave(filename = paste0(outDir, "MaxAnomAbundRealms_stacked.pdf"), plot = MaxAnomAbundRealms, width = 200, height = 300, units = "mm", dpi = 300)
ggsave(filename = paste0(outDir, "MaxAnomAbundRealms_sidebyside.pdf"), plot = MaxAnomAbundRealms, width = 400, height = 150, units = "mm", dpi = 300)


## Richness, Max Anomaly ##
## Tropical ##

nd4_trop <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelRich_trop$data$StdTmaxAnomalyRS),
                       to = max(MaxAnomalyModelRich_trop$data$StdTmaxAnomalyRS),
                       length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelRich_trop$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")))

# back transform the predictors
nd4_trop$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4_trop$StdTmaxAnomalyRS,
  originalX = trop$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd4_trop$LogRich <- 0
nd4_trop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# # Record for later use (see 'Values' or check the data set itself for the row in Coleoptera that satisfies the conditions)
refRow <- which((nd4_trop$LUI=="Primary vegetation") & (nd4_trop$StdTmaxAnomaly==min(abs(nd4_trop$StdTmaxAnomaly))))
# 56th row

QPV <- quantile(x = MaxAnomalyModelRich_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_trop$data$LUI=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAnomalyModelRich_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_trop$data$LUI=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAnomalyModelRich_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_trop$data$LUI=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAnomalyModelRich_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_trop$data$LUI=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
sr.preds.tmax.trop <- PredictGLMERRandIter(model = MaxAnomalyModelRich_trop$model,data = nd4_trop)

# back transform the abundance values
sr.preds.tmax.trop <- exp(sr.preds.tmax.trop)-0.01


# split by order into matrices, then name them
number_of_chunks = 6
list_sr.preds.tmax.trop <- lapply(seq(1, NROW(sr.preds.tmax.trop), ceiling(NROW(sr.preds.tmax.trop)/number_of_chunks)),
                                  function(i) sr.preds.tmax.trop[i:min(i + ceiling(NROW(sr.preds.tmax.trop)/number_of_chunks) - 1, NROW(sr.preds.tmax.trop)),])

names(list_sr.preds.tmax.trop) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")
list2env(list_sr.preds.tmax.trop,globalenv())


# sweep out refRow
list_sr.preds.tmax.trop <- lapply(list_sr.preds.tmax.trop,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[56,],FUN="/") 
})

list2env(list_sr.preds.tmax.trop,globalenv())

# split nd4_trop by order
Order<- paste0("nd4_trop_",nd4_trop$Order)
# create a list of data frames
by_Order <- split(nd4_trop,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles
Coleoptera[which((nd4_trop_Coleoptera$LUI=="Primary vegetation") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd4_trop_Coleoptera$LUI=="Primary vegetation") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd4_trop_Coleoptera$LUI=="Secondary vegetation") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd4_trop_Coleoptera$LUI=="Secondary vegetation") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd4_trop_Coleoptera$LUI=="Agriculture_Low") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd4_trop_Coleoptera$LUI=="Agriculture_Low") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd4_trop_Coleoptera$LUI=="Agriculture_High") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd4_trop_Coleoptera$LUI=="Agriculture_High") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd4_trop_Diptera$LUI=="Primary vegetation") & (nd4_trop_Diptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd4_trop_Diptera$LUI=="Primary vegetation") & (nd4_trop_Diptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd4_trop_Diptera$LUI=="Secondary vegetation") & (nd4_trop_Diptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd4_trop_Diptera$LUI=="Secondary vegetation") & (nd4_trop_Diptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd4_trop_Diptera$LUI=="Agriculture_Low") & (nd4_trop_Diptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd4_trop_Diptera$LUI=="Agriculture_Low") & (nd4_trop_Diptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd4_trop_Diptera$LUI=="Agriculture_High") & (nd4_trop_Diptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd4_trop_Diptera$LUI=="Agriculture_High") & (nd4_trop_Diptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd4_trop_Hemiptera$LUI=="Primary vegetation") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd4_trop_Hemiptera$LUI=="Primary vegetation") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd4_trop_Hemiptera$LUI=="Secondary vegetation") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd4_trop_Hemiptera$LUI=="Secondary vegetation") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd4_trop_Hemiptera$LUI=="Agriculture_Low") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd4_trop_Hemiptera$LUI=="Agriculture_Low") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd4_trop_Hemiptera$LUI=="Agriculture_High") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd4_trop_Hemiptera$LUI=="Agriculture_High") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd4_trop_Hymenoptera$LUI=="Primary vegetation") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd4_trop_Hymenoptera$LUI=="Primary vegetation") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd4_trop_Hymenoptera$LUI=="Secondary vegetation") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd4_trop_Hymenoptera$LUI=="Secondary vegetation") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd4_trop_Hymenoptera$LUI=="Agriculture_Low") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd4_trop_Hymenoptera$LUI=="Agriculture_Low") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd4_trop_Hymenoptera$LUI=="Agriculture_High") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd4_trop_Hymenoptera$LUI=="Agriculture_High") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd4_trop_Lepidoptera$LUI=="Primary vegetation") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd4_trop_Lepidoptera$LUI=="Primary vegetation") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd4_trop_Lepidoptera$LUI=="Secondary vegetation") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd4_trop_Lepidoptera$LUI=="Secondary vegetation") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd4_trop_Lepidoptera$LUI=="Agriculture_Low") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd4_trop_Lepidoptera$LUI=="Agriculture_Low") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd4_trop_Lepidoptera$LUI=="Agriculture_High") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd4_trop_Lepidoptera$LUI=="Agriculture_High") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd4_trop_Orthoptera$LUI=="Primary vegetation") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd4_trop_Orthoptera$LUI=="Primary vegetation") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd4_trop_Orthoptera$LUI=="Secondary vegetation") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd4_trop_Orthoptera$LUI=="Secondary vegetation") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd4_trop_Orthoptera$LUI=="Agriculture_Low") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd4_trop_Orthoptera$LUI=="Agriculture_Low") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd4_trop_Orthoptera$LUI=="Agriculture_High") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd4_trop_Orthoptera$LUI=="Agriculture_High") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd4_trop_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd4_trop_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_trop_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_trop_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                       FUN = median,na.rm=TRUE))*100)-100
nd4_trop_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                      FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_trop_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                      FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_trop_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                         FUN = median,na.rm=TRUE))*100)-100
nd4_trop_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_trop_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_trop_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                           FUN = median,na.rm=TRUE))*100)-100
nd4_trop_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_trop_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_trop_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                           FUN = median,na.rm=TRUE))*100)-100
nd4_trop_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_trop_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_trop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd4_trop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_trop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# set factor levels
nd4_trop_Coleoptera$LUI <- factor(nd4_trop_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Diptera$LUI <- factor(nd4_trop_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Hemiptera$LUI <- factor(nd4_trop_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Hymenoptera$LUI <- factor(nd4_trop_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Lepidoptera$LUI <- factor(nd4_trop_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Orthoptera$LUI <- factor(nd4_trop_Orthoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot

p_coleoptera <- ggplot(data = nd4_trop_Coleoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_trop_Coleoptera$PredLower, ymax = nd4_trop_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
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

p_diptera <- ggplot(data = nd4_trop_Diptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_trop_Diptera$PredLower, ymax = nd4_trop_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
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

p_hemiptera <- ggplot(data = nd4_trop_Hemiptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_trop_Hemiptera$PredLower, ymax = nd4_trop_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
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

p_hymenoptera <- ggplot(data = nd4_trop_Hymenoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_trop_Hymenoptera$PredLower, ymax = nd4_trop_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
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

p_lepidoptera <- ggplot(data = nd4_trop_Lepidoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_trop_Lepidoptera$PredLower, ymax = nd4_trop_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
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

p_orthoptera <- ggplot(data = nd4_trop_Orthoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_trop_Orthoptera$PredLower, ymax = nd4_trop_Orthoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
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
  ggtitle("Orthoptera")


# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(), 
          legend.text = element_text(size = 6), 
          legend.title = element_blank())
)


# put them all together to save them
MaxAnomRich_trop <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,p_orthoptera)
#MaxAnomRich_trop <- cowplot::plot_grid(MaxAnomRich_trop,legend,ncol=1,rel_heights = c(1,0.1))

# save them
#ggsave(filename = paste0(outDir, "MaxAnomRich_trop.pdf"), plot = MaxAnomRich_trop, width = 200, height = 150, units = "mm", dpi = 300)

## NonTropical ##

nd4_nontrop <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelRich_nontrop$data$StdTmaxAnomalyRS),
                       to = max(MaxAnomalyModelRich_nontrop$data$StdTmaxAnomalyRS),
                       length.out = 100),
  LUI=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelRich_nontrop$data$LUI)),
  Order=factor(c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")))

# back transform the predictors
nd4_nontrop$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4_nontrop$StdTmaxAnomalyRS,
  originalX = nontrop$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd4_nontrop$LogRich <- 0
nd4_nontrop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# Record for later use (see 'Values' or check the data set itself for the row in Coleoptera that satisfies the conditions)
refRow <- which((nd4_nontrop$LUI=="Primary vegetation") & (nd4_nontrop$StdTmaxAnomaly==min(abs(nd4_nontrop$StdTmaxAnomaly))))
# 15th row

# adjust plot 1: max anomaly and abundance

QPV <- quantile(x = MaxAnomalyModelRich_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_nontrop$data$LUI=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAnomalyModelRich_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_nontrop$data$LUI=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAnomalyModelRich_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_nontrop$data$LUI=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAnomalyModelRich_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_nontrop$data$LUI=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
sr.preds.tmax.nontrop <- PredictGLMERRandIter(model = MaxAnomalyModelRich_nontrop$model,data = nd4_nontrop)

# back transform the abundance values
sr.preds.tmax.nontrop <- exp(sr.preds.tmax.nontrop)-0.01


# split by order into matrices, then name them
number_of_chunks = 6
list_sr.preds.tmax.nontrop <- lapply(seq(1, NROW(sr.preds.tmax.nontrop), ceiling(NROW(sr.preds.tmax.nontrop)/number_of_chunks)),
                                     function(i) sr.preds.tmax.nontrop[i:min(i + ceiling(NROW(sr.preds.tmax.nontrop)/number_of_chunks) - 1, NROW(sr.preds.tmax.nontrop)),])

names(list_sr.preds.tmax.nontrop) <- c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")
list2env(list_sr.preds.tmax.nontrop,globalenv())

# sweep out refRow
list_sr.preds.tmax.nontrop <- lapply(list_sr.preds.tmax.nontrop,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[15,],FUN="/") 
})

list2env(list_sr.preds.tmax.nontrop,globalenv())

# split nd4_nontrop by order
Order<- paste0("nd4_nontrop_",nd4_nontrop$Order)
# create a list of data frames
by_Order <- split(nd4_nontrop,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles
Coleoptera[which((nd4_nontrop_Coleoptera$LUI=="Primary vegetation") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd4_nontrop_Coleoptera$LUI=="Primary vegetation") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd4_nontrop_Coleoptera$LUI=="Secondary vegetation") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd4_nontrop_Coleoptera$LUI=="Secondary vegetation") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd4_nontrop_Coleoptera$LUI=="Agriculture_Low") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd4_nontrop_Coleoptera$LUI=="Agriculture_Low") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd4_nontrop_Coleoptera$LUI=="Agriculture_High") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd4_nontrop_Coleoptera$LUI=="Agriculture_High") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd4_nontrop_Diptera$LUI=="Primary vegetation") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd4_nontrop_Diptera$LUI=="Primary vegetation") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd4_nontrop_Diptera$LUI=="Secondary vegetation") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd4_nontrop_Diptera$LUI=="Secondary vegetation") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd4_nontrop_Diptera$LUI=="Agriculture_Low") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd4_nontrop_Diptera$LUI=="Agriculture_Low") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd4_nontrop_Diptera$LUI=="Agriculture_High") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd4_nontrop_Diptera$LUI=="Agriculture_High") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd4_nontrop_Hemiptera$LUI=="Primary vegetation") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd4_nontrop_Hemiptera$LUI=="Primary vegetation") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd4_nontrop_Hemiptera$LUI=="Secondary vegetation") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd4_nontrop_Hemiptera$LUI=="Secondary vegetation") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd4_nontrop_Hemiptera$LUI=="Agriculture_Low") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd4_nontrop_Hemiptera$LUI=="Agriculture_Low") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd4_nontrop_Hemiptera$LUI=="Agriculture_High") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd4_nontrop_Hemiptera$LUI=="Agriculture_High") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd4_nontrop_Hymenoptera$LUI=="Primary vegetation") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd4_nontrop_Hymenoptera$LUI=="Primary vegetation") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd4_nontrop_Hymenoptera$LUI=="Secondary vegetation") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd4_nontrop_Hymenoptera$LUI=="Secondary vegetation") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd4_nontrop_Hymenoptera$LUI=="Agriculture_Low") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd4_nontrop_Hymenoptera$LUI=="Agriculture_Low") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd4_nontrop_Hymenoptera$LUI=="Agriculture_High") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd4_nontrop_Hymenoptera$LUI=="Agriculture_High") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd4_nontrop_Lepidoptera$LUI=="Primary vegetation") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd4_nontrop_Lepidoptera$LUI=="Primary vegetation") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd4_nontrop_Lepidoptera$LUI=="Secondary vegetation") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd4_nontrop_Lepidoptera$LUI=="Secondary vegetation") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd4_nontrop_Lepidoptera$LUI=="Agriculture_Low") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd4_nontrop_Lepidoptera$LUI=="Agriculture_Low") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd4_nontrop_Lepidoptera$LUI=="Agriculture_High") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd4_nontrop_Lepidoptera$LUI=="Agriculture_High") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd4_nontrop_Orthoptera$LUI=="Primary vegetation") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd4_nontrop_Orthoptera$LUI=="Primary vegetation") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd4_nontrop_Orthoptera$LUI=="Secondary vegetation") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd4_nontrop_Orthoptera$LUI=="Secondary vegetation") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd4_nontrop_Orthoptera$LUI=="Agriculture_Low") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd4_nontrop_Orthoptera$LUI=="Agriculture_Low") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd4_nontrop_Orthoptera$LUI=="Agriculture_High") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd4_nontrop_Orthoptera$LUI=="Agriculture_High") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd4_nontrop_Coleoptera$PredMedian <- ((apply(X = Coleoptera,MARGIN = 1,
                                             FUN = median,na.rm=TRUE))*100)-100
nd4_nontrop_Coleoptera$PredUpper <- ((apply(X = Coleoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_nontrop_Coleoptera$PredLower <- ((apply(X = Coleoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_nontrop_Diptera$PredMedian <- ((apply(X = Diptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd4_nontrop_Diptera$PredUpper <- ((apply(X = Diptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_nontrop_Diptera$PredLower <- ((apply(X = Diptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_nontrop_Hemiptera$PredMedian <- ((apply(X = Hemiptera,MARGIN = 1,
                                            FUN = median,na.rm=TRUE))*100)-100
nd4_nontrop_Hemiptera$PredUpper <- ((apply(X = Hemiptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_nontrop_Hemiptera$PredLower <- ((apply(X = Hemiptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_nontrop_Hymenoptera$PredMedian <- ((apply(X = Hymenoptera,MARGIN = 1,
                                              FUN = median,na.rm=TRUE))*100)-100
nd4_nontrop_Hymenoptera$PredUpper <- ((apply(X = Hymenoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_nontrop_Hymenoptera$PredLower <- ((apply(X = Hymenoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_nontrop_Lepidoptera$PredMedian <- ((apply(X = Lepidoptera,MARGIN = 1,
                                              FUN = median,na.rm=TRUE))*100)-100
nd4_nontrop_Lepidoptera$PredUpper <- ((apply(X = Lepidoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_nontrop_Lepidoptera$PredLower <- ((apply(X = Lepidoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_nontrop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                             FUN = median,na.rm=TRUE))*100)-100
nd4_nontrop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_nontrop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
nd4_nontrop_Coleoptera$LUI <- factor(nd4_nontrop_Coleoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Diptera$LUI <- factor(nd4_nontrop_Diptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Hemiptera$LUI <- factor(nd4_nontrop_Hemiptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Hymenoptera$LUI <- factor(nd4_nontrop_Hymenoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Lepidoptera$LUI <- factor(nd4_nontrop_Lepidoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Orthoptera$LUI <- factor(nd4_nontrop_Orthoptera$LUI, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# plot

p_coleoptera <- ggplot(data = nd4_nontrop_Coleoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_nontrop_Coleoptera$PredLower, ymax = nd4_nontrop_Coleoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
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

p_diptera <- ggplot(data = nd4_nontrop_Diptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_nontrop_Diptera$PredLower, ymax = nd4_nontrop_Diptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
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

p_hemiptera <- ggplot(data = nd4_nontrop_Hemiptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_nontrop_Hemiptera$PredLower, ymax = nd4_nontrop_Hemiptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
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

p_hymenoptera <- ggplot(data = nd4_nontrop_Hymenoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_nontrop_Hymenoptera$PredLower, ymax = nd4_nontrop_Hymenoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
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

p_lepidoptera <- ggplot(data = nd4_nontrop_Lepidoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_nontrop_Lepidoptera$PredLower, ymax = nd4_nontrop_Lepidoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
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

p_orthoptera <- ggplot(data = nd4_nontrop_Orthoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = LUI), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_nontrop_Orthoptera$PredLower, ymax = nd4_nontrop_Orthoptera$PredUpper, fill = LUI), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300), limits = c(-100, 300)) +
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
  ggtitle("Orthoptera")

# get the legend
legend <- get_legend(
  p_coleoptera +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.background = element_blank(), 
          legend.text = element_text(size = 6), 
          legend.title = element_blank())
)


# put them all together to save them
MaxAnomRich_nontrop <- cowplot::plot_grid(p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,p_orthoptera)
#MaxAnomRich_nontrop <- cowplot::plot_grid(MaxAnomRich_nontrop,legend,ncol=1,rel_heights = c(1,0.1))

# save them
#ggsave(filename = paste0(outDir, "MaxAnomRich_nontrop.pdf"), plot = MaxAnomRich_nontrop, width = 200, height = 150, units = "mm", dpi = 300)
#ggsave(filename = paste0(outDir, "MaxAnomRich_nontrop_extended yaxis.pdf"), plot = MaxAnomRich_nontrop, width = 200, height = 150, units = "mm", dpi = 300)

# plot Tropical and NonTropical plots together
# to do this, comment the lines wherein cowplot is used to add a legend to each plot and where ggsave is used to save individual plots
# re-run from the start, to this point
# comment lines according to which style of plot is desired (stacked or side-by-side)

MaxAnomRichRealms <-cowplot::plot_grid(MaxAnomRich_nontrop, NULL, MaxAnomRich_trop, ncol=1, labels = c('A',"",'B'),rel_heights = c(1,0.1,1)) # stacked
#MaxAnomRichRealms <-cowplot::plot_grid(MaxAnomRich_nontrop, NULL, MaxAnomRich_trop, ncol=3, labels = c('A',"",'B'),rel_widths = c(1,0.1,1)) # side-by-side

MaxAnomRichRealms <-cowplot::plot_grid(MaxAnomRichRealms, legend, ncol=1, rel_heights = c(1,0.1))

ggsave(filename = paste0(outDir, "MaxAnomRichRealms_stacked.pdf"), plot = MaxAnomRichRealms, width = 200, height = 300, units = "mm", dpi = 300)
#ggsave(filename = paste0(outDir, "MaxAnomRichRealms_sidebyside.pdf"), plot = MaxAnomRichRealms, width = 400, height = 150, units = "mm", dpi = 300)

