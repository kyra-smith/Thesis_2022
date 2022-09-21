
## Plot results ##

# get set up
# directories 
inDir<- "C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/"
outDir <- "C:/Users/Kyra/Documents/GLITRS/Code/6_TropicalModels/Plots/"
if(!dir.exists(outDir)) dir.create(outDir)

# load libraries
library(devtools)
library(StatisticalModels)
library(predictsFunctions)
source("C:/Users/Kyra/Documents/GLITRS/Data/0_Functions.R")
library(sjPlot)
library(cowplot)

predictsSites <- readRDS(file = paste0(inDir,"PREDICTSSiteData.rds"))
trop <- readRDS(file = paste0(inDir,"trop.rds"))
nontrop <- readRDS(file = paste0(inDir,"nontrop.rds"))
load(paste0(inDir, "/MeanAnomalyModelAbund_trop.rdata"))
load(paste0(inDir, "/MeanAnomalyModelRich_trop.rdata"))
load(paste0(inDir, "/MaxAnomalyModelAbund_trop.rdata"))
load(paste0(inDir, "/MaxAnomalyModelRich_trop.rdata"))
load(paste0(inDir, "/MeanAnomalyModelAbund_nontrop.rdata"))
load(paste0(inDir, "/MeanAnomalyModelRich_nontrop.rdata"))
load(paste0(inDir, "/MaxAnomalyModelAbund_nontrop.rdata"))
load(paste0(inDir, "/MaxAnomalyModelRich_nontrop.rdata"))

# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

## Abundance, Mean Anomaly ##
  ## Tropical ##

nd_trop <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund_trop$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund_trop$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund_trop$data$UI2)),
  Order=factor(c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")))

# back transform the predictors
nd_trop$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd_trop$StdTmeanAnomalyRS,
  originalX = trop$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd_trop$LogAbund <- 0
nd_trop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# does this for each Order
# need to record it to use later on
refRow <- which((nd_trop$UI2=="Primary vegetation") & (nd_trop$StdTmeanAnomaly==min(abs(nd_trop$StdTmeanAnomaly))))
  # the first row, every 400 rows

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelAbund_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_trop$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelAbund_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_trop$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelAbund_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_trop$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelAbund_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_trop$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean.trop <- PredictGLMERRandIter(model = MeanAnomalyModelAbund_trop$model,data = nd_trop)

# back transform the abundance values
a.preds.tmean.trop <- exp(a.preds.tmean.trop)-0.01


# split up by order
number_of_chunks = 11
list_a.preds.tmean.trop <- lapply(seq(1, NROW(a.preds.tmean.trop), ceiling(NROW(a.preds.tmean.trop)/number_of_chunks)),
                             function(i) a.preds.tmean.trop[i:min(i + ceiling(NROW(a.preds.tmean.trop)/number_of_chunks) - 1, NROW(a.preds.tmean.trop)),])
# success!
# creates list of matrices
# name them
names(list_a.preds.tmean.trop) <- c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")
list2env(list_a.preds.tmean.trop,globalenv())


# tim's suggestion
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

Blattodea[which((nd_trop_Blattodea$UI2=="Primary vegetation") & (nd_trop_Blattodea$StdTmeanAnomalyRS > QPV[2])),] <- NA
Blattodea[which((nd_trop_Blattodea$UI2=="Primary vegetation") & (nd_trop_Blattodea$StdTmeanAnomalyRS < QPV[1])),] <- NA
Blattodea[which((nd_trop_Blattodea$UI2=="Second_tropary vegetation") & (nd_trop_Blattodea$StdTmeanAnomalyRS < QSV[1])),] <- NA
Blattodea[which((nd_trop_Blattodea$UI2=="Second_tropary vegetation") & (nd_trop_Blattodea$StdTmeanAnomalyRS > QSV[2])),] <- NA
Blattodea[which((nd_trop_Blattodea$UI2=="Agriculture_Low") & (nd_trop_Blattodea$StdTmeanAnomalyRS < QAL[1])),] <- NA
Blattodea[which((nd_trop_Blattodea$UI2=="Agriculture_Low") & (nd_trop_Blattodea$StdTmeanAnomalyRS > QAL[2])),] <- NA
Blattodea[which((nd_trop_Blattodea$UI2=="Agriculture_High") & (nd_trop_Blattodea$StdTmeanAnomalyRS < QAH[1])),] <- NA
Blattodea[which((nd_trop_Blattodea$UI2=="Agriculture_High") & (nd_trop_Blattodea$StdTmeanAnomalyRS > QAH[2])),] <- NA

Coleoptera[which((nd_trop_Coleoptera$UI2=="Primary vegetation") & (nd_trop_Coleoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd_trop_Coleoptera$UI2=="Primary vegetation") & (nd_trop_Coleoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd_trop_Coleoptera$UI2=="Second_tropary vegetation") & (nd_trop_Coleoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd_trop_Coleoptera$UI2=="Second_tropary vegetation") & (nd_trop_Coleoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd_trop_Coleoptera$UI2=="Agriculture_Low") & (nd_trop_Coleoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd_trop_Coleoptera$UI2=="Agriculture_Low") & (nd_trop_Coleoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd_trop_Coleoptera$UI2=="Agriculture_High") & (nd_trop_Coleoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd_trop_Coleoptera$UI2=="Agriculture_High") & (nd_trop_Coleoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd_trop_Diptera$UI2=="Primary vegetation") & (nd_trop_Diptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd_trop_Diptera$UI2=="Primary vegetation") & (nd_trop_Diptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd_trop_Diptera$UI2=="Second_tropary vegetation") & (nd_trop_Diptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd_trop_Diptera$UI2=="Second_tropary vegetation") & (nd_trop_Diptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd_trop_Diptera$UI2=="Agriculture_Low") & (nd_trop_Diptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd_trop_Diptera$UI2=="Agriculture_Low") & (nd_trop_Diptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd_trop_Diptera$UI2=="Agriculture_High") & (nd_trop_Diptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd_trop_Diptera$UI2=="Agriculture_High") & (nd_trop_Diptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd_trop_Hemiptera$UI2=="Primary vegetation") & (nd_trop_Hemiptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd_trop_Hemiptera$UI2=="Primary vegetation") & (nd_trop_Hemiptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd_trop_Hemiptera$UI2=="Second_tropary vegetation") & (nd_trop_Hemiptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd_trop_Hemiptera$UI2=="Second_tropary vegetation") & (nd_trop_Hemiptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd_trop_Hemiptera$UI2=="Agriculture_Low") & (nd_trop_Hemiptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd_trop_Hemiptera$UI2=="Agriculture_Low") & (nd_trop_Hemiptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd_trop_Hemiptera$UI2=="Agriculture_High") & (nd_trop_Hemiptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd_trop_Hemiptera$UI2=="Agriculture_High") & (nd_trop_Hemiptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd_trop_Hymenoptera$UI2=="Primary vegetation") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd_trop_Hymenoptera$UI2=="Primary vegetation") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd_trop_Hymenoptera$UI2=="Second_tropary vegetation") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd_trop_Hymenoptera$UI2=="Second_tropary vegetation") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd_trop_Hymenoptera$UI2=="Agriculture_Low") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd_trop_Hymenoptera$UI2=="Agriculture_Low") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd_trop_Hymenoptera$UI2=="Agriculture_High") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd_trop_Hymenoptera$UI2=="Agriculture_High") & (nd_trop_Hymenoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd_trop_Lepidoptera$UI2=="Primary vegetation") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd_trop_Lepidoptera$UI2=="Primary vegetation") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd_trop_Lepidoptera$UI2=="Second_tropary vegetation") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd_trop_Lepidoptera$UI2=="Second_tropary vegetation") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd_trop_Lepidoptera$UI2=="Agriculture_Low") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd_trop_Lepidoptera$UI2=="Agriculture_Low") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd_trop_Lepidoptera$UI2=="Agriculture_High") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd_trop_Lepidoptera$UI2=="Agriculture_High") & (nd_trop_Lepidoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Neuroptera[which((nd_trop_Neuroptera$UI2=="Primary vegetation") & (nd_trop_Neuroptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Neuroptera[which((nd_trop_Neuroptera$UI2=="Primary vegetation") & (nd_trop_Neuroptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Neuroptera[which((nd_trop_Neuroptera$UI2=="Second_tropary vegetation") & (nd_trop_Neuroptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Neuroptera[which((nd_trop_Neuroptera$UI2=="Second_tropary vegetation") & (nd_trop_Neuroptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Neuroptera[which((nd_trop_Neuroptera$UI2=="Agriculture_Low") & (nd_trop_Neuroptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Neuroptera[which((nd_trop_Neuroptera$UI2=="Agriculture_Low") & (nd_trop_Neuroptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Neuroptera[which((nd_trop_Neuroptera$UI2=="Agriculture_High") & (nd_trop_Neuroptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Neuroptera[which((nd_trop_Neuroptera$UI2=="Agriculture_High") & (nd_trop_Neuroptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd_trop_Orthoptera$UI2=="Primary vegetation") & (nd_trop_Orthoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd_trop_Orthoptera$UI2=="Primary vegetation") & (nd_trop_Orthoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd_trop_Orthoptera$UI2=="Second_tropary vegetation") & (nd_trop_Orthoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd_trop_Orthoptera$UI2=="Second_tropary vegetation") & (nd_trop_Orthoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd_trop_Orthoptera$UI2=="Agriculture_Low") & (nd_trop_Orthoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd_trop_Orthoptera$UI2=="Agriculture_Low") & (nd_trop_Orthoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd_trop_Orthoptera$UI2=="Agriculture_High") & (nd_trop_Orthoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd_trop_Orthoptera$UI2=="Agriculture_High") & (nd_trop_Orthoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Other[which((nd_trop_Other$UI2=="Primary vegetation") & (nd_trop_Other$StdTmeanAnomalyRS > QPV[2])),] <- NA
Other[which((nd_trop_Other$UI2=="Primary vegetation") & (nd_trop_Other$StdTmeanAnomalyRS < QPV[1])),] <- NA
Other[which((nd_trop_Other$UI2=="Second_tropary vegetation") & (nd_trop_Other$StdTmeanAnomalyRS < QSV[1])),] <- NA
Other[which((nd_trop_Other$UI2=="Second_tropary vegetation") & (nd_trop_Other$StdTmeanAnomalyRS > QSV[2])),] <- NA
Other[which((nd_trop_Other$UI2=="Agriculture_Low") & (nd_trop_Other$StdTmeanAnomalyRS < QAL[1])),] <- NA
Other[which((nd_trop_Other$UI2=="Agriculture_Low") & (nd_trop_Other$StdTmeanAnomalyRS > QAL[2])),] <- NA
Other[which((nd_trop_Other$UI2=="Agriculture_High") & (nd_trop_Other$StdTmeanAnomalyRS < QAH[1])),] <- NA
Other[which((nd_trop_Other$UI2=="Agriculture_High") & (nd_trop_Other$StdTmeanAnomalyRS > QAH[2])),] <- NA

Thysanoptera[which((nd_trop_Thysanoptera$UI2=="Primary vegetation") & (nd_trop_Thysanoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Thysanoptera[which((nd_trop_Thysanoptera$UI2=="Primary vegetation") & (nd_trop_Thysanoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Thysanoptera[which((nd_trop_Thysanoptera$UI2=="Second_tropary vegetation") & (nd_trop_Thysanoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Thysanoptera[which((nd_trop_Thysanoptera$UI2=="Second_tropary vegetation") & (nd_trop_Thysanoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Thysanoptera[which((nd_trop_Thysanoptera$UI2=="Agriculture_Low") & (nd_trop_Thysanoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Thysanoptera[which((nd_trop_Thysanoptera$UI2=="Agriculture_Low") & (nd_trop_Thysanoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Thysanoptera[which((nd_trop_Thysanoptera$UI2=="Agriculture_High") & (nd_trop_Thysanoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Thysanoptera[which((nd_trop_Thysanoptera$UI2=="Agriculture_High") & (nd_trop_Thysanoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Trichoptera[which((nd_trop_Trichoptera$UI2=="Primary vegetation") & (nd_trop_Trichoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Trichoptera[which((nd_trop_Trichoptera$UI2=="Primary vegetation") & (nd_trop_Trichoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Trichoptera[which((nd_trop_Trichoptera$UI2=="Second_tropary vegetation") & (nd_trop_Trichoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Trichoptera[which((nd_trop_Trichoptera$UI2=="Second_tropary vegetation") & (nd_trop_Trichoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Trichoptera[which((nd_trop_Trichoptera$UI2=="Agriculture_Low") & (nd_trop_Trichoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Trichoptera[which((nd_trop_Trichoptera$UI2=="Agriculture_Low") & (nd_trop_Trichoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Trichoptera[which((nd_trop_Trichoptera$UI2=="Agriculture_High") & (nd_trop_Trichoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Trichoptera[which((nd_trop_Trichoptera$UI2=="Agriculture_High") & (nd_trop_Trichoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd_trop_Blattodea$PredMedian <- ((apply(X = Blattodea,MARGIN = 1,
                                   FUN = median,na.rm=TRUE))*100)-100
nd_trop_Blattodea$PredUpper <- ((apply(X = Blattodea,MARGIN = 1,
                                  FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_trop_Blattodea$PredLower <- ((apply(X = Blattodea,MARGIN = 1,
                                  FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

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

nd_trop_Neuroptera$PredMedian <- ((apply(X = Neuroptera,MARGIN = 1,
                                    FUN = median,na.rm=TRUE))*100)-100
nd_trop_Neuroptera$PredUpper <- ((apply(X = Neuroptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_trop_Neuroptera$PredLower <- ((apply(X = Neuroptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_trop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                    FUN = median,na.rm=TRUE))*100)-100
nd_trop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_trop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                   FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_trop_Other$PredMedian <- ((apply(X = Other,MARGIN = 1,
                               FUN = median,na.rm=TRUE))*100)-100
nd_trop_Other$PredUpper <- ((apply(X = Other,MARGIN = 1,
                              FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_trop_Other$PredLower <- ((apply(X = Other,MARGIN = 1,
                              FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_trop_Thysanoptera$PredMedian <- ((apply(X = Thysanoptera,MARGIN = 1,
                                      FUN = median,na.rm=TRUE))*100)-100
nd_trop_Thysanoptera$PredUpper <- ((apply(X = Thysanoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_trop_Thysanoptera$PredLower <- ((apply(X = Thysanoptera,MARGIN = 1,
                                     FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_trop_Trichoptera$PredMedian <- ((apply(X = Trichoptera,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd_trop_Trichoptera$PredUpper <- ((apply(X = Trichoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_trop_Trichoptera$PredLower <- ((apply(X = Trichoptera,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
# nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Blattodea$UI2 <- factor(nd_trop_Blattodea$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Coleoptera$UI2 <- factor(nd_trop_Coleoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Diptera$UI2 <- factor(nd_trop_Diptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Hemiptera$UI2 <- factor(nd_trop_Hemiptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Hymenoptera$UI2 <- factor(nd_trop_Hymenoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Lepidoptera$UI2 <- factor(nd_trop_Lepidoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Neuroptera$UI2 <- factor(nd_trop_Neuroptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Orthoptera$UI2 <- factor(nd_trop_Orthoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Other$UI2 <- factor(nd_trop_Other$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Thysanoptera$UI2 <- factor(nd_trop_Thysanoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_trop_Trichoptera$UI2 <- factor(nd_trop_Trichoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# add zone factor
nd_trop_Blattodea$Zone <- as.factor("Tropical")
nd_trop_Coleoptera$Zone <- as.factor("Tropical")
nd_trop_Diptera$Zone <- as.factor("Tropical")
nd_trop_Hemiptera$Zone <- as.factor("Tropical")
nd_trop_Hymenoptera$Zone <- as.factor("Tropical")
nd_trop_Lepidoptera$Zone <- as.factor("Tropical")
nd_trop_Neuroptera$Zone <- as.factor("Tropical")
nd_trop_Orthoptera$Zone <- as.factor("Tropical")
nd_trop_Other$Zone <- as.factor("Tropical")
nd_trop_Thysanoptera$Zone <- as.factor("Tropical")
nd_trop_Trichoptera$Zone <- as.factor("Tropical")

  ## NonTropical ##

nd_nontrop <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund_nontrop$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund_nontrop$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund_nontrop$data$UI2)),
  Order=factor(c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")))

# back transform the predictors
nd_nontrop$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd_nontrop$StdTmeanAnomalyRS,
  originalX = nontrop$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd_nontrop$LogAbund <- 0
nd_nontrop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# does this for each Order
# reference row is 4th row, every 400 rows (see 'Values')
refRow <- which((nd_nontrop$UI2=="Primary vegetation") & (nd_nontrop$StdTmeanAnomaly==min(abs(nd_nontrop$StdTmeanAnomaly))))
# the first row, every 400 rows

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelAbund_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_nontrop$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelAbund_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_nontrop$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelAbund_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_nontrop$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelAbund_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_nontrop$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean.nontrop <- PredictGLMERRandIter(model = MeanAnomalyModelAbund_nontrop$model,data = nd_nontrop)

# back transform the abundance values
a.preds.tmean.nontrop <- exp(a.preds.tmean.nontrop)-0.01


# another try!
number_of_chunks = 11
list_a.preds.tmean.nontrop <- lapply(seq(1, NROW(a.preds.tmean.nontrop), ceiling(NROW(a.preds.tmean.nontrop)/number_of_chunks)),
                                  function(i) a.preds.tmean.nontrop[i:min(i + ceiling(NROW(a.preds.tmean.nontrop)/number_of_chunks) - 1, NROW(a.preds.tmean.nontrop)),])
# success!
# creates list of matrices
# name them
names(list_a.preds.tmean.nontrop) <- c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")
list2env(list_a.preds.tmean.nontrop,globalenv())

# tim's suggestion
list_a.preds.tmean.nontrop <- lapply(list_a.preds.tmean.nontrop,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[1,],FUN="/") 
})

list2env(list_a.preds.tmean.nontrop,globalenv())

# works!

# Charlie's input: if there are a few, use facet.wrap, if there are more, use ggplot and then cowplot

# split nd_nontrop by order
Order<- paste0("nd_nontrop_",nd_nontrop$Order)
# create a list of data frames
by_Order <- split(nd_nontrop,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles

Blattodea[which((nd_nontrop_Blattodea$UI2=="Primary vegetation") & (nd_nontrop_Blattodea$StdTmeanAnomalyRS > QPV[2])),] <- NA
Blattodea[which((nd_nontrop_Blattodea$UI2=="Primary vegetation") & (nd_nontrop_Blattodea$StdTmeanAnomalyRS < QPV[1])),] <- NA
Blattodea[which((nd_nontrop_Blattodea$UI2=="Second_nontropary vegetation") & (nd_nontrop_Blattodea$StdTmeanAnomalyRS < QSV[1])),] <- NA
Blattodea[which((nd_nontrop_Blattodea$UI2=="Second_nontropary vegetation") & (nd_nontrop_Blattodea$StdTmeanAnomalyRS > QSV[2])),] <- NA
Blattodea[which((nd_nontrop_Blattodea$UI2=="Agriculture_Low") & (nd_nontrop_Blattodea$StdTmeanAnomalyRS < QAL[1])),] <- NA
Blattodea[which((nd_nontrop_Blattodea$UI2=="Agriculture_Low") & (nd_nontrop_Blattodea$StdTmeanAnomalyRS > QAL[2])),] <- NA
Blattodea[which((nd_nontrop_Blattodea$UI2=="Agriculture_High") & (nd_nontrop_Blattodea$StdTmeanAnomalyRS < QAH[1])),] <- NA
Blattodea[which((nd_nontrop_Blattodea$UI2=="Agriculture_High") & (nd_nontrop_Blattodea$StdTmeanAnomalyRS > QAH[2])),] <- NA

Coleoptera[which((nd_nontrop_Coleoptera$UI2=="Primary vegetation") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd_nontrop_Coleoptera$UI2=="Primary vegetation") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd_nontrop_Coleoptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd_nontrop_Coleoptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd_nontrop_Coleoptera$UI2=="Agriculture_Low") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd_nontrop_Coleoptera$UI2=="Agriculture_Low") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd_nontrop_Coleoptera$UI2=="Agriculture_High") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd_nontrop_Coleoptera$UI2=="Agriculture_High") & (nd_nontrop_Coleoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd_nontrop_Diptera$UI2=="Primary vegetation") & (nd_nontrop_Diptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd_nontrop_Diptera$UI2=="Primary vegetation") & (nd_nontrop_Diptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd_nontrop_Diptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Diptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd_nontrop_Diptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Diptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd_nontrop_Diptera$UI2=="Agriculture_Low") & (nd_nontrop_Diptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd_nontrop_Diptera$UI2=="Agriculture_Low") & (nd_nontrop_Diptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd_nontrop_Diptera$UI2=="Agriculture_High") & (nd_nontrop_Diptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd_nontrop_Diptera$UI2=="Agriculture_High") & (nd_nontrop_Diptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd_nontrop_Hemiptera$UI2=="Primary vegetation") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd_nontrop_Hemiptera$UI2=="Primary vegetation") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd_nontrop_Hemiptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd_nontrop_Hemiptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd_nontrop_Hemiptera$UI2=="Agriculture_Low") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd_nontrop_Hemiptera$UI2=="Agriculture_Low") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd_nontrop_Hemiptera$UI2=="Agriculture_High") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd_nontrop_Hemiptera$UI2=="Agriculture_High") & (nd_nontrop_Hemiptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd_nontrop_Hymenoptera$UI2=="Primary vegetation") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd_nontrop_Hymenoptera$UI2=="Primary vegetation") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd_nontrop_Hymenoptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd_nontrop_Hymenoptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd_nontrop_Hymenoptera$UI2=="Agriculture_Low") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd_nontrop_Hymenoptera$UI2=="Agriculture_Low") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd_nontrop_Hymenoptera$UI2=="Agriculture_High") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd_nontrop_Hymenoptera$UI2=="Agriculture_High") & (nd_nontrop_Hymenoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd_nontrop_Lepidoptera$UI2=="Primary vegetation") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd_nontrop_Lepidoptera$UI2=="Primary vegetation") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd_nontrop_Lepidoptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd_nontrop_Lepidoptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd_nontrop_Lepidoptera$UI2=="Agriculture_Low") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd_nontrop_Lepidoptera$UI2=="Agriculture_Low") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd_nontrop_Lepidoptera$UI2=="Agriculture_High") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd_nontrop_Lepidoptera$UI2=="Agriculture_High") & (nd_nontrop_Lepidoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Neuroptera[which((nd_nontrop_Neuroptera$UI2=="Primary vegetation") & (nd_nontrop_Neuroptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Neuroptera[which((nd_nontrop_Neuroptera$UI2=="Primary vegetation") & (nd_nontrop_Neuroptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Neuroptera[which((nd_nontrop_Neuroptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Neuroptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Neuroptera[which((nd_nontrop_Neuroptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Neuroptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Neuroptera[which((nd_nontrop_Neuroptera$UI2=="Agriculture_Low") & (nd_nontrop_Neuroptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Neuroptera[which((nd_nontrop_Neuroptera$UI2=="Agriculture_Low") & (nd_nontrop_Neuroptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Neuroptera[which((nd_nontrop_Neuroptera$UI2=="Agriculture_High") & (nd_nontrop_Neuroptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Neuroptera[which((nd_nontrop_Neuroptera$UI2=="Agriculture_High") & (nd_nontrop_Neuroptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd_nontrop_Orthoptera$UI2=="Primary vegetation") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd_nontrop_Orthoptera$UI2=="Primary vegetation") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd_nontrop_Orthoptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd_nontrop_Orthoptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd_nontrop_Orthoptera$UI2=="Agriculture_Low") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd_nontrop_Orthoptera$UI2=="Agriculture_Low") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd_nontrop_Orthoptera$UI2=="Agriculture_High") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd_nontrop_Orthoptera$UI2=="Agriculture_High") & (nd_nontrop_Orthoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Other[which((nd_nontrop_Other$UI2=="Primary vegetation") & (nd_nontrop_Other$StdTmeanAnomalyRS > QPV[2])),] <- NA
Other[which((nd_nontrop_Other$UI2=="Primary vegetation") & (nd_nontrop_Other$StdTmeanAnomalyRS < QPV[1])),] <- NA
Other[which((nd_nontrop_Other$UI2=="Second_nontropary vegetation") & (nd_nontrop_Other$StdTmeanAnomalyRS < QSV[1])),] <- NA
Other[which((nd_nontrop_Other$UI2=="Second_nontropary vegetation") & (nd_nontrop_Other$StdTmeanAnomalyRS > QSV[2])),] <- NA
Other[which((nd_nontrop_Other$UI2=="Agriculture_Low") & (nd_nontrop_Other$StdTmeanAnomalyRS < QAL[1])),] <- NA
Other[which((nd_nontrop_Other$UI2=="Agriculture_Low") & (nd_nontrop_Other$StdTmeanAnomalyRS > QAL[2])),] <- NA
Other[which((nd_nontrop_Other$UI2=="Agriculture_High") & (nd_nontrop_Other$StdTmeanAnomalyRS < QAH[1])),] <- NA
Other[which((nd_nontrop_Other$UI2=="Agriculture_High") & (nd_nontrop_Other$StdTmeanAnomalyRS > QAH[2])),] <- NA

Thysanoptera[which((nd_nontrop_Thysanoptera$UI2=="Primary vegetation") & (nd_nontrop_Thysanoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Thysanoptera[which((nd_nontrop_Thysanoptera$UI2=="Primary vegetation") & (nd_nontrop_Thysanoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Thysanoptera[which((nd_nontrop_Thysanoptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Thysanoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Thysanoptera[which((nd_nontrop_Thysanoptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Thysanoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Thysanoptera[which((nd_nontrop_Thysanoptera$UI2=="Agriculture_Low") & (nd_nontrop_Thysanoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Thysanoptera[which((nd_nontrop_Thysanoptera$UI2=="Agriculture_Low") & (nd_nontrop_Thysanoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Thysanoptera[which((nd_nontrop_Thysanoptera$UI2=="Agriculture_High") & (nd_nontrop_Thysanoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Thysanoptera[which((nd_nontrop_Thysanoptera$UI2=="Agriculture_High") & (nd_nontrop_Thysanoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Trichoptera[which((nd_nontrop_Trichoptera$UI2=="Primary vegetation") & (nd_nontrop_Trichoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Trichoptera[which((nd_nontrop_Trichoptera$UI2=="Primary vegetation") & (nd_nontrop_Trichoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Trichoptera[which((nd_nontrop_Trichoptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Trichoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Trichoptera[which((nd_nontrop_Trichoptera$UI2=="Second_nontropary vegetation") & (nd_nontrop_Trichoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Trichoptera[which((nd_nontrop_Trichoptera$UI2=="Agriculture_Low") & (nd_nontrop_Trichoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Trichoptera[which((nd_nontrop_Trichoptera$UI2=="Agriculture_Low") & (nd_nontrop_Trichoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Trichoptera[which((nd_nontrop_Trichoptera$UI2=="Agriculture_High") & (nd_nontrop_Trichoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Trichoptera[which((nd_nontrop_Trichoptera$UI2=="Agriculture_High") & (nd_nontrop_Trichoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd_nontrop_Blattodea$PredMedian <- ((apply(X = Blattodea,MARGIN = 1,
                                        FUN = median,na.rm=TRUE))*100)-100
nd_nontrop_Blattodea$PredUpper <- ((apply(X = Blattodea,MARGIN = 1,
                                       FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_nontrop_Blattodea$PredLower <- ((apply(X = Blattodea,MARGIN = 1,
                                       FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

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

nd_nontrop_Neuroptera$PredMedian <- ((apply(X = Neuroptera,MARGIN = 1,
                                         FUN = median,na.rm=TRUE))*100)-100
nd_nontrop_Neuroptera$PredUpper <- ((apply(X = Neuroptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_nontrop_Neuroptera$PredLower <- ((apply(X = Neuroptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_nontrop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                         FUN = median,na.rm=TRUE))*100)-100
nd_nontrop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_nontrop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_nontrop_Other$PredMedian <- ((apply(X = Other,MARGIN = 1,
                                    FUN = median,na.rm=TRUE))*100)-100
nd_nontrop_Other$PredUpper <- ((apply(X = Other,MARGIN = 1,
                                   FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_nontrop_Other$PredLower <- ((apply(X = Other,MARGIN = 1,
                                   FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_nontrop_Thysanoptera$PredMedian <- ((apply(X = Thysanoptera,MARGIN = 1,
                                           FUN = median,na.rm=TRUE))*100)-100
nd_nontrop_Thysanoptera$PredUpper <- ((apply(X = Thysanoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_nontrop_Thysanoptera$PredLower <- ((apply(X = Thysanoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd_nontrop_Trichoptera$PredMedian <- ((apply(X = Trichoptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd_nontrop_Trichoptera$PredUpper <- ((apply(X = Trichoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd_nontrop_Trichoptera$PredLower <- ((apply(X = Trichoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
# nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Blattodea$UI2 <- factor(nd_nontrop_Blattodea$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Coleoptera$UI2 <- factor(nd_nontrop_Coleoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Diptera$UI2 <- factor(nd_nontrop_Diptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Hemiptera$UI2 <- factor(nd_nontrop_Hemiptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Hymenoptera$UI2 <- factor(nd_nontrop_Hymenoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Lepidoptera$UI2 <- factor(nd_nontrop_Lepidoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Neuroptera$UI2 <- factor(nd_nontrop_Neuroptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Orthoptera$UI2 <- factor(nd_nontrop_Orthoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Other$UI2 <- factor(nd_nontrop_Other$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Thysanoptera$UI2 <- factor(nd_nontrop_Thysanoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd_nontrop_Trichoptera$UI2 <- factor(nd_nontrop_Trichoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# add zone factor
nd_nontrop_Blattodea$Zone <- as.factor("NonTropical")
nd_nontrop_Coleoptera$Zone <- as.factor("NonTropical")
nd_nontrop_Diptera$Zone <- as.factor("NonTropical")
nd_nontrop_Hemiptera$Zone <- as.factor("NonTropical")
nd_nontrop_Hymenoptera$Zone <- as.factor("NonTropical")
nd_nontrop_Lepidoptera$Zone <- as.factor("NonTropical")
nd_nontrop_Neuroptera$Zone <- as.factor("NonTropical")
nd_nontrop_Orthoptera$Zone <- as.factor("NonTropical")
nd_nontrop_Other$Zone <- as.factor("NonTropical")
nd_nontrop_Thysanoptera$Zone <- as.factor("NonTropical")
nd_nontrop_Trichoptera$Zone <- as.factor("NonTropical")

# put use rbind to add nd_nontrop to nd_trop to make one data table for plotting

nd_Blattodea <- rbind(nd_trop_Blattodea,nd_nontrop_Blattodea)
nd_Coleoptera <- rbind(nd_trop_Coleoptera,nd_nontrop_Coleoptera)
nd_Diptera <- rbind(nd_trop_Diptera,nd_nontrop_Diptera)
nd_Hemiptera <- rbind(nd_trop_Hemiptera,nd_nontrop_Hemiptera)
nd_Hymenoptera <- rbind(nd_trop_Hymenoptera,nd_nontrop_Hymenoptera)
nd_Lepidoptera <- rbind(nd_trop_Lepidoptera,nd_nontrop_Lepidoptera)
nd_Neuroptera <- rbind(nd_trop_Neuroptera,nd_nontrop_Neuroptera)
nd_Orthoptera <- rbind(nd_trop_Orthoptera,nd_nontrop_Orthoptera)
nd_Other <- rbind(nd_trop_Other,nd_nontrop_Other)
nd_Thysanoptera <- rbind(nd_trop_Thysanoptera,nd_nontrop_Thysanoptera)
nd_Trichoptera <- rbind(nd_trop_Trichoptera,nd_nontrop_Trichoptera)

# plot

p_blattodea <- ggplot(data = nd_Blattodea, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Blattodea$PredLower, ymax = nd_Blattodea$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) + 
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) + # extended y axis for orders with high larger confidence intervals
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        # legend.background = element_blank(), 
        # legend.text = element_text(size = 6), 
        # legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Blattodea")

p_coleoptera <- ggplot(data = nd_Coleoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Coleoptera$PredLower, ymax = nd_Coleoptera$PredUpper, fill = UI2), alpha = 0.2) +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Diptera$PredLower, ymax = nd_Diptera$PredUpper, fill = UI2), alpha = 0.2) +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Hemiptera$PredLower, ymax = nd_Hemiptera$PredUpper, fill = UI2), alpha = 0.2) +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Hymenoptera$PredLower, ymax = nd_Hymenoptera$PredUpper, fill = UI2), alpha = 0.2) +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Lepidoptera$PredLower, ymax = nd_Lepidoptera$PredUpper, fill = UI2), alpha = 0.2) +
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

p_neuroptera <- ggplot(data = nd_Neuroptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Neuroptera$PredLower, ymax = nd_Neuroptera$PredUpper, fill = UI2), alpha = 0.2) +
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
  ggtitle("Neuroptera")

p_other <- ggplot(data = nd_Other, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Other$PredLower, ymax = nd_Other$PredUpper, fill = UI2), alpha = 0.2) +
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
  ggtitle("Other")

p_orthoptera <- ggplot(data = nd_Orthoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Orthoptera$PredLower, ymax = nd_Orthoptera$PredUpper, fill = UI2), alpha = 0.2) +
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
  ggtitle("Orthoptera")

p_thysanoptera <- ggplot(data = nd_Thysanoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Thysanoptera$PredLower, ymax = nd_Thysanoptera$PredUpper, fill = UI2), alpha = 0.2) +
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
  ggtitle("Thysanoptera")

p_trichoptera <- ggplot(data = nd_Trichoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd_Trichoptera$PredLower, ymax = nd_Trichoptera$PredUpper, fill = UI2), alpha = 0.2) +
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
  ggtitle("Trichoptera")

# get the legend
legend <- get_legend(
  p_blattodea +
    guides(color = guide_legend(nrow = 1),
           linetype = guide_legend (nrow=1)) +
    theme(legend.position = "right",
          legend.background = element_blank(), 
          legend.text = element_text(size = 6), 
          legend.title = element_blank())
)


# put them all together to save them
MeanAnomAbund <- cowplot::plot_grid(p_blattodea, p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,p_neuroptera,p_orthoptera,p_other,p_thysanoptera,p_trichoptera,legend)

# save the ggplots
ggsave(filename = paste0(outDir, "MeanAnomAbund.pdf"), plot = MeanAnomAbund, width = 200, height = 150, units = "mm", dpi = 300)
#ggsave(filename = paste0(outDir, "MeanAnomAbund_extended yaxis.pdf"), plot = MeanAnomRich, width = 200, height = 150, units = "mm", dpi = 300)

## Richness, Mean Anomaly ##
## Tropical ##

nd2_trop <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelRich_trop$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelRich_trop$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelRich_trop$data$UI2)),
  Order=factor(c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")))

# back transform the predictors
nd2_trop$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2_trop$StdTmeanAnomalyRS,
  originalX = trop$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd2_trop$LogAbund <- 0
nd2_trop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# does this for each Order
# need to record it to use later on
refRow <- which((nd2_trop$UI2=="Primary vegetation") & (nd2_trop$StdTmeanAnomaly==min(abs(nd2_trop$StdTmeanAnomaly))))
# the first row, every 400 rows

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelRich_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_trop$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelRich_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_trop$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelRich_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_trop$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelRich_trop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_trop$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
sr.preds.tmean.trop <- PredictGLMERRandIter(model = MeanAnomalyModelRich_trop$model,data = nd2_trop)

# back transform the abundance values
sr.preds.tmean.trop <- exp(sr.preds.tmean.trop)-0.01


# split up by order
number_of_chunks = 11
list_sr.preds.tmean.trop <- lapply(seq(1, NROW(sr.preds.tmean.trop), ceiling(NROW(sr.preds.tmean.trop)/number_of_chunks)),
                                   function(i) sr.preds.tmean.trop[i:min(i + ceiling(NROW(sr.preds.tmean.trop)/number_of_chunks) - 1, NROW(sr.preds.tmean.trop)),])
# success!
# creates list of matrices
# name them
names(list_sr.preds.tmean.trop) <- c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")
list2env(list_sr.preds.tmean.trop,globalenv())


# tim's suggestion
list_sr.preds.tmean.trop <- lapply(list_sr.preds.tmean.trop,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[1,],FUN="/") 
})

list2env(list_sr.preds.tmean.trop,globalenv())

# works!

# Charlie's input: if there are a few, use facet.wrap, if there are more, use ggplot and then cowplot

# split nd2_trop by order
Order<- paste0("nd2_trop_",nd2_trop$Order)
# create a list of data frames
by_Order <- split(nd2_trop,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles

Blattodea[which((nd2_trop_Blattodea$UI2=="Primary vegetation") & (nd2_trop_Blattodea$StdTmeanAnomalyRS > QPV[2])),] <- NA
Blattodea[which((nd2_trop_Blattodea$UI2=="Primary vegetation") & (nd2_trop_Blattodea$StdTmeanAnomalyRS < QPV[1])),] <- NA
Blattodea[which((nd2_trop_Blattodea$UI2=="Second2_tropary vegetation") & (nd2_trop_Blattodea$StdTmeanAnomalyRS < QSV[1])),] <- NA
Blattodea[which((nd2_trop_Blattodea$UI2=="Second2_tropary vegetation") & (nd2_trop_Blattodea$StdTmeanAnomalyRS > QSV[2])),] <- NA
Blattodea[which((nd2_trop_Blattodea$UI2=="Agriculture_Low") & (nd2_trop_Blattodea$StdTmeanAnomalyRS < QAL[1])),] <- NA
Blattodea[which((nd2_trop_Blattodea$UI2=="Agriculture_Low") & (nd2_trop_Blattodea$StdTmeanAnomalyRS > QAL[2])),] <- NA
Blattodea[which((nd2_trop_Blattodea$UI2=="Agriculture_High") & (nd2_trop_Blattodea$StdTmeanAnomalyRS < QAH[1])),] <- NA
Blattodea[which((nd2_trop_Blattodea$UI2=="Agriculture_High") & (nd2_trop_Blattodea$StdTmeanAnomalyRS > QAH[2])),] <- NA

Coleoptera[which((nd2_trop_Coleoptera$UI2=="Primary vegetation") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd2_trop_Coleoptera$UI2=="Primary vegetation") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd2_trop_Coleoptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd2_trop_Coleoptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd2_trop_Coleoptera$UI2=="Agriculture_Low") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd2_trop_Coleoptera$UI2=="Agriculture_Low") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd2_trop_Coleoptera$UI2=="Agriculture_High") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd2_trop_Coleoptera$UI2=="Agriculture_High") & (nd2_trop_Coleoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd2_trop_Diptera$UI2=="Primary vegetation") & (nd2_trop_Diptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd2_trop_Diptera$UI2=="Primary vegetation") & (nd2_trop_Diptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd2_trop_Diptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Diptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd2_trop_Diptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Diptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd2_trop_Diptera$UI2=="Agriculture_Low") & (nd2_trop_Diptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd2_trop_Diptera$UI2=="Agriculture_Low") & (nd2_trop_Diptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd2_trop_Diptera$UI2=="Agriculture_High") & (nd2_trop_Diptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd2_trop_Diptera$UI2=="Agriculture_High") & (nd2_trop_Diptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd2_trop_Hemiptera$UI2=="Primary vegetation") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd2_trop_Hemiptera$UI2=="Primary vegetation") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd2_trop_Hemiptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd2_trop_Hemiptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd2_trop_Hemiptera$UI2=="Agriculture_Low") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd2_trop_Hemiptera$UI2=="Agriculture_Low") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd2_trop_Hemiptera$UI2=="Agriculture_High") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd2_trop_Hemiptera$UI2=="Agriculture_High") & (nd2_trop_Hemiptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd2_trop_Hymenoptera$UI2=="Primary vegetation") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd2_trop_Hymenoptera$UI2=="Primary vegetation") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd2_trop_Hymenoptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd2_trop_Hymenoptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd2_trop_Hymenoptera$UI2=="Agriculture_Low") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd2_trop_Hymenoptera$UI2=="Agriculture_Low") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd2_trop_Hymenoptera$UI2=="Agriculture_High") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd2_trop_Hymenoptera$UI2=="Agriculture_High") & (nd2_trop_Hymenoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd2_trop_Lepidoptera$UI2=="Primary vegetation") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd2_trop_Lepidoptera$UI2=="Primary vegetation") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd2_trop_Lepidoptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd2_trop_Lepidoptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd2_trop_Lepidoptera$UI2=="Agriculture_Low") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd2_trop_Lepidoptera$UI2=="Agriculture_Low") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd2_trop_Lepidoptera$UI2=="Agriculture_High") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd2_trop_Lepidoptera$UI2=="Agriculture_High") & (nd2_trop_Lepidoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Neuroptera[which((nd2_trop_Neuroptera$UI2=="Primary vegetation") & (nd2_trop_Neuroptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Neuroptera[which((nd2_trop_Neuroptera$UI2=="Primary vegetation") & (nd2_trop_Neuroptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Neuroptera[which((nd2_trop_Neuroptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Neuroptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Neuroptera[which((nd2_trop_Neuroptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Neuroptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Neuroptera[which((nd2_trop_Neuroptera$UI2=="Agriculture_Low") & (nd2_trop_Neuroptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Neuroptera[which((nd2_trop_Neuroptera$UI2=="Agriculture_Low") & (nd2_trop_Neuroptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Neuroptera[which((nd2_trop_Neuroptera$UI2=="Agriculture_High") & (nd2_trop_Neuroptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Neuroptera[which((nd2_trop_Neuroptera$UI2=="Agriculture_High") & (nd2_trop_Neuroptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd2_trop_Orthoptera$UI2=="Primary vegetation") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd2_trop_Orthoptera$UI2=="Primary vegetation") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd2_trop_Orthoptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd2_trop_Orthoptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd2_trop_Orthoptera$UI2=="Agriculture_Low") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd2_trop_Orthoptera$UI2=="Agriculture_Low") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd2_trop_Orthoptera$UI2=="Agriculture_High") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd2_trop_Orthoptera$UI2=="Agriculture_High") & (nd2_trop_Orthoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Other[which((nd2_trop_Other$UI2=="Primary vegetation") & (nd2_trop_Other$StdTmeanAnomalyRS > QPV[2])),] <- NA
Other[which((nd2_trop_Other$UI2=="Primary vegetation") & (nd2_trop_Other$StdTmeanAnomalyRS < QPV[1])),] <- NA
Other[which((nd2_trop_Other$UI2=="Second2_tropary vegetation") & (nd2_trop_Other$StdTmeanAnomalyRS < QSV[1])),] <- NA
Other[which((nd2_trop_Other$UI2=="Second2_tropary vegetation") & (nd2_trop_Other$StdTmeanAnomalyRS > QSV[2])),] <- NA
Other[which((nd2_trop_Other$UI2=="Agriculture_Low") & (nd2_trop_Other$StdTmeanAnomalyRS < QAL[1])),] <- NA
Other[which((nd2_trop_Other$UI2=="Agriculture_Low") & (nd2_trop_Other$StdTmeanAnomalyRS > QAL[2])),] <- NA
Other[which((nd2_trop_Other$UI2=="Agriculture_High") & (nd2_trop_Other$StdTmeanAnomalyRS < QAH[1])),] <- NA
Other[which((nd2_trop_Other$UI2=="Agriculture_High") & (nd2_trop_Other$StdTmeanAnomalyRS > QAH[2])),] <- NA

Thysanoptera[which((nd2_trop_Thysanoptera$UI2=="Primary vegetation") & (nd2_trop_Thysanoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Thysanoptera[which((nd2_trop_Thysanoptera$UI2=="Primary vegetation") & (nd2_trop_Thysanoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Thysanoptera[which((nd2_trop_Thysanoptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Thysanoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Thysanoptera[which((nd2_trop_Thysanoptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Thysanoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Thysanoptera[which((nd2_trop_Thysanoptera$UI2=="Agriculture_Low") & (nd2_trop_Thysanoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Thysanoptera[which((nd2_trop_Thysanoptera$UI2=="Agriculture_Low") & (nd2_trop_Thysanoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Thysanoptera[which((nd2_trop_Thysanoptera$UI2=="Agriculture_High") & (nd2_trop_Thysanoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Thysanoptera[which((nd2_trop_Thysanoptera$UI2=="Agriculture_High") & (nd2_trop_Thysanoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Trichoptera[which((nd2_trop_Trichoptera$UI2=="Primary vegetation") & (nd2_trop_Trichoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Trichoptera[which((nd2_trop_Trichoptera$UI2=="Primary vegetation") & (nd2_trop_Trichoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Trichoptera[which((nd2_trop_Trichoptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Trichoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Trichoptera[which((nd2_trop_Trichoptera$UI2=="Second2_tropary vegetation") & (nd2_trop_Trichoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Trichoptera[which((nd2_trop_Trichoptera$UI2=="Agriculture_Low") & (nd2_trop_Trichoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Trichoptera[which((nd2_trop_Trichoptera$UI2=="Agriculture_Low") & (nd2_trop_Trichoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Trichoptera[which((nd2_trop_Trichoptera$UI2=="Agriculture_High") & (nd2_trop_Trichoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Trichoptera[which((nd2_trop_Trichoptera$UI2=="Agriculture_High") & (nd2_trop_Trichoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd2_trop_Blattodea$PredMedian <- ((apply(X = Blattodea,MARGIN = 1,
                                         FUN = median,na.rm=TRUE))*100)-100
nd2_trop_Blattodea$PredUpper <- ((apply(X = Blattodea,MARGIN = 1,
                                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_trop_Blattodea$PredLower <- ((apply(X = Blattodea,MARGIN = 1,
                                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

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

nd2_trop_Neuroptera$PredMedian <- ((apply(X = Neuroptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd2_trop_Neuroptera$PredUpper <- ((apply(X = Neuroptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_trop_Neuroptera$PredLower <- ((apply(X = Neuroptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_trop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd2_trop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_trop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_trop_Other$PredMedian <- ((apply(X = Other,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd2_trop_Other$PredUpper <- ((apply(X = Other,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_trop_Other$PredLower <- ((apply(X = Other,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_trop_Thysanoptera$PredMedian <- ((apply(X = Thysanoptera,MARGIN = 1,
                                            FUN = median,na.rm=TRUE))*100)-100
nd2_trop_Thysanoptera$PredUpper <- ((apply(X = Thysanoptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_trop_Thysanoptera$PredLower <- ((apply(X = Thysanoptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_trop_Trichoptera$PredMedian <- ((apply(X = Trichoptera,MARGIN = 1,
                                           FUN = median,na.rm=TRUE))*100)-100
nd2_trop_Trichoptera$PredUpper <- ((apply(X = Trichoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_trop_Trichoptera$PredLower <- ((apply(X = Trichoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
# nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Blattodea$UI2 <- factor(nd2_trop_Blattodea$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Coleoptera$UI2 <- factor(nd2_trop_Coleoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Diptera$UI2 <- factor(nd2_trop_Diptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Hemiptera$UI2 <- factor(nd2_trop_Hemiptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Hymenoptera$UI2 <- factor(nd2_trop_Hymenoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Lepidoptera$UI2 <- factor(nd2_trop_Lepidoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Neuroptera$UI2 <- factor(nd2_trop_Neuroptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Orthoptera$UI2 <- factor(nd2_trop_Orthoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Other$UI2 <- factor(nd2_trop_Other$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Thysanoptera$UI2 <- factor(nd2_trop_Thysanoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_trop_Trichoptera$UI2 <- factor(nd2_trop_Trichoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# add zone factor
nd2_trop_Blattodea$Zone <- as.factor("Tropical")
nd2_trop_Coleoptera$Zone <- as.factor("Tropical")
nd2_trop_Diptera$Zone <- as.factor("Tropical")
nd2_trop_Hemiptera$Zone <- as.factor("Tropical")
nd2_trop_Hymenoptera$Zone <- as.factor("Tropical")
nd2_trop_Lepidoptera$Zone <- as.factor("Tropical")
nd2_trop_Neuroptera$Zone <- as.factor("Tropical")
nd2_trop_Orthoptera$Zone <- as.factor("Tropical")
nd2_trop_Other$Zone <- as.factor("Tropical")
nd2_trop_Thysanoptera$Zone <- as.factor("Tropical")
nd2_trop_Trichoptera$Zone <- as.factor("Tropical")

## NonTropical ##

nd2_nontrop <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelRich_nontrop$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelRich_nontrop$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelRich_nontrop$data$UI2)),
  Order=factor(c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")))

# back transform the predictors
nd2_nontrop$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2_nontrop$StdTmeanAnomalyRS,
  originalX = nontrop$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd2_nontrop$LogAbund <- 0
nd2_nontrop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# does this for each Order
# reference row is 4th row, every 400 rows (see 'Values')
refRow <- which((nd2_nontrop$UI2=="Primary vegetation") & (nd2_nontrop$StdTmeanAnomaly==min(abs(nd2_nontrop$StdTmeanAnomaly))))
# the first row, every 400 rows

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelRich_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_nontrop$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelRich_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_nontrop$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelRich_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_nontrop$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelRich_nontrop$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich_nontrop$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
sr.preds.tmean.nontrop <- PredictGLMERRandIter(model = MeanAnomalyModelRich_nontrop$model,data = nd2_nontrop)

# back transform the abundance values
sr.preds.tmean.nontrop <- exp(sr.preds.tmean.nontrop)-0.01


# another try!
number_of_chunks = 11
list_sr.preds.tmean.nontrop <- lapply(seq(1, NROW(sr.preds.tmean.nontrop), ceiling(NROW(sr.preds.tmean.nontrop)/number_of_chunks)),
                                      function(i) sr.preds.tmean.nontrop[i:min(i + ceiling(NROW(sr.preds.tmean.nontrop)/number_of_chunks) - 1, NROW(sr.preds.tmean.nontrop)),])
# success!
# creates list of matrices
# name them
names(list_sr.preds.tmean.nontrop) <- c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")
list2env(list_sr.preds.tmean.nontrop,globalenv())

# tim's suggestion
list_sr.preds.tmean.nontrop <- lapply(list_sr.preds.tmean.nontrop,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[1,],FUN="/") 
})

list2env(list_sr.preds.tmean.nontrop,globalenv())

# works!

# Charlie's input: if there are a few, use facet.wrap, if there are more, use ggplot and then cowplot

# split nd2_nontrop by order
Order<- paste0("nd2_nontrop_",nd2_nontrop$Order)
# create a list of data frames
by_Order <- split(nd2_nontrop,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles

Blattodea[which((nd2_nontrop_Blattodea$UI2=="Primary vegetation") & (nd2_nontrop_Blattodea$StdTmeanAnomalyRS > QPV[2])),] <- NA
Blattodea[which((nd2_nontrop_Blattodea$UI2=="Primary vegetation") & (nd2_nontrop_Blattodea$StdTmeanAnomalyRS < QPV[1])),] <- NA
Blattodea[which((nd2_nontrop_Blattodea$UI2=="Secondary vegetation") & (nd2_nontrop_Blattodea$StdTmeanAnomalyRS < QSV[1])),] <- NA
Blattodea[which((nd2_nontrop_Blattodea$UI2=="Secondary vegetation") & (nd2_nontrop_Blattodea$StdTmeanAnomalyRS > QSV[2])),] <- NA
Blattodea[which((nd2_nontrop_Blattodea$UI2=="Agriculture_Low") & (nd2_nontrop_Blattodea$StdTmeanAnomalyRS < QAL[1])),] <- NA
Blattodea[which((nd2_nontrop_Blattodea$UI2=="Agriculture_Low") & (nd2_nontrop_Blattodea$StdTmeanAnomalyRS > QAL[2])),] <- NA
Blattodea[which((nd2_nontrop_Blattodea$UI2=="Agriculture_High") & (nd2_nontrop_Blattodea$StdTmeanAnomalyRS < QAH[1])),] <- NA
Blattodea[which((nd2_nontrop_Blattodea$UI2=="Agriculture_High") & (nd2_nontrop_Blattodea$StdTmeanAnomalyRS > QAH[2])),] <- NA

Coleoptera[which((nd2_nontrop_Coleoptera$UI2=="Primary vegetation") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd2_nontrop_Coleoptera$UI2=="Primary vegetation") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd2_nontrop_Coleoptera$UI2=="Secondary vegetation") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd2_nontrop_Coleoptera$UI2=="Secondary vegetation") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd2_nontrop_Coleoptera$UI2=="Agriculture_Low") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd2_nontrop_Coleoptera$UI2=="Agriculture_Low") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd2_nontrop_Coleoptera$UI2=="Agriculture_High") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd2_nontrop_Coleoptera$UI2=="Agriculture_High") & (nd2_nontrop_Coleoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd2_nontrop_Diptera$UI2=="Primary vegetation") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd2_nontrop_Diptera$UI2=="Primary vegetation") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd2_nontrop_Diptera$UI2=="Secondary vegetation") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd2_nontrop_Diptera$UI2=="Secondary vegetation") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd2_nontrop_Diptera$UI2=="Agriculture_Low") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd2_nontrop_Diptera$UI2=="Agriculture_Low") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd2_nontrop_Diptera$UI2=="Agriculture_High") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd2_nontrop_Diptera$UI2=="Agriculture_High") & (nd2_nontrop_Diptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd2_nontrop_Hemiptera$UI2=="Primary vegetation") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd2_nontrop_Hemiptera$UI2=="Primary vegetation") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd2_nontrop_Hemiptera$UI2=="Secondary vegetation") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd2_nontrop_Hemiptera$UI2=="Secondary vegetation") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd2_nontrop_Hemiptera$UI2=="Agriculture_Low") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd2_nontrop_Hemiptera$UI2=="Agriculture_Low") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd2_nontrop_Hemiptera$UI2=="Agriculture_High") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd2_nontrop_Hemiptera$UI2=="Agriculture_High") & (nd2_nontrop_Hemiptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd2_nontrop_Hymenoptera$UI2=="Primary vegetation") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd2_nontrop_Hymenoptera$UI2=="Primary vegetation") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd2_nontrop_Hymenoptera$UI2=="Secondary vegetation") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd2_nontrop_Hymenoptera$UI2=="Secondary vegetation") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd2_nontrop_Hymenoptera$UI2=="Agriculture_Low") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd2_nontrop_Hymenoptera$UI2=="Agriculture_Low") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd2_nontrop_Hymenoptera$UI2=="Agriculture_High") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd2_nontrop_Hymenoptera$UI2=="Agriculture_High") & (nd2_nontrop_Hymenoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd2_nontrop_Lepidoptera$UI2=="Primary vegetation") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd2_nontrop_Lepidoptera$UI2=="Primary vegetation") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd2_nontrop_Lepidoptera$UI2=="Secondary vegetation") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd2_nontrop_Lepidoptera$UI2=="Secondary vegetation") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd2_nontrop_Lepidoptera$UI2=="Agriculture_Low") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd2_nontrop_Lepidoptera$UI2=="Agriculture_Low") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd2_nontrop_Lepidoptera$UI2=="Agriculture_High") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd2_nontrop_Lepidoptera$UI2=="Agriculture_High") & (nd2_nontrop_Lepidoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Neuroptera[which((nd2_nontrop_Neuroptera$UI2=="Primary vegetation") & (nd2_nontrop_Neuroptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Neuroptera[which((nd2_nontrop_Neuroptera$UI2=="Primary vegetation") & (nd2_nontrop_Neuroptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Neuroptera[which((nd2_nontrop_Neuroptera$UI2=="Secondary vegetation") & (nd2_nontrop_Neuroptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Neuroptera[which((nd2_nontrop_Neuroptera$UI2=="Secondary vegetation") & (nd2_nontrop_Neuroptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Neuroptera[which((nd2_nontrop_Neuroptera$UI2=="Agriculture_Low") & (nd2_nontrop_Neuroptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Neuroptera[which((nd2_nontrop_Neuroptera$UI2=="Agriculture_Low") & (nd2_nontrop_Neuroptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Neuroptera[which((nd2_nontrop_Neuroptera$UI2=="Agriculture_High") & (nd2_nontrop_Neuroptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Neuroptera[which((nd2_nontrop_Neuroptera$UI2=="Agriculture_High") & (nd2_nontrop_Neuroptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd2_nontrop_Orthoptera$UI2=="Primary vegetation") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd2_nontrop_Orthoptera$UI2=="Primary vegetation") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd2_nontrop_Orthoptera$UI2=="Secondary vegetation") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd2_nontrop_Orthoptera$UI2=="Secondary vegetation") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd2_nontrop_Orthoptera$UI2=="Agriculture_Low") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd2_nontrop_Orthoptera$UI2=="Agriculture_Low") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd2_nontrop_Orthoptera$UI2=="Agriculture_High") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd2_nontrop_Orthoptera$UI2=="Agriculture_High") & (nd2_nontrop_Orthoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Other[which((nd2_nontrop_Other$UI2=="Primary vegetation") & (nd2_nontrop_Other$StdTmeanAnomalyRS > QPV[2])),] <- NA
Other[which((nd2_nontrop_Other$UI2=="Primary vegetation") & (nd2_nontrop_Other$StdTmeanAnomalyRS < QPV[1])),] <- NA
Other[which((nd2_nontrop_Other$UI2=="Secondary vegetation") & (nd2_nontrop_Other$StdTmeanAnomalyRS < QSV[1])),] <- NA
Other[which((nd2_nontrop_Other$UI2=="Secondary vegetation") & (nd2_nontrop_Other$StdTmeanAnomalyRS > QSV[2])),] <- NA
Other[which((nd2_nontrop_Other$UI2=="Agriculture_Low") & (nd2_nontrop_Other$StdTmeanAnomalyRS < QAL[1])),] <- NA
Other[which((nd2_nontrop_Other$UI2=="Agriculture_Low") & (nd2_nontrop_Other$StdTmeanAnomalyRS > QAL[2])),] <- NA
Other[which((nd2_nontrop_Other$UI2=="Agriculture_High") & (nd2_nontrop_Other$StdTmeanAnomalyRS < QAH[1])),] <- NA
Other[which((nd2_nontrop_Other$UI2=="Agriculture_High") & (nd2_nontrop_Other$StdTmeanAnomalyRS > QAH[2])),] <- NA

Thysanoptera[which((nd2_nontrop_Thysanoptera$UI2=="Primary vegetation") & (nd2_nontrop_Thysanoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Thysanoptera[which((nd2_nontrop_Thysanoptera$UI2=="Primary vegetation") & (nd2_nontrop_Thysanoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Thysanoptera[which((nd2_nontrop_Thysanoptera$UI2=="Secondary vegetation") & (nd2_nontrop_Thysanoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Thysanoptera[which((nd2_nontrop_Thysanoptera$UI2=="Secondary vegetation") & (nd2_nontrop_Thysanoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Thysanoptera[which((nd2_nontrop_Thysanoptera$UI2=="Agriculture_Low") & (nd2_nontrop_Thysanoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Thysanoptera[which((nd2_nontrop_Thysanoptera$UI2=="Agriculture_Low") & (nd2_nontrop_Thysanoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Thysanoptera[which((nd2_nontrop_Thysanoptera$UI2=="Agriculture_High") & (nd2_nontrop_Thysanoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Thysanoptera[which((nd2_nontrop_Thysanoptera$UI2=="Agriculture_High") & (nd2_nontrop_Thysanoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA

Trichoptera[which((nd2_nontrop_Trichoptera$UI2=="Primary vegetation") & (nd2_nontrop_Trichoptera$StdTmeanAnomalyRS > QPV[2])),] <- NA
Trichoptera[which((nd2_nontrop_Trichoptera$UI2=="Primary vegetation") & (nd2_nontrop_Trichoptera$StdTmeanAnomalyRS < QPV[1])),] <- NA
Trichoptera[which((nd2_nontrop_Trichoptera$UI2=="Secondary vegetation") & (nd2_nontrop_Trichoptera$StdTmeanAnomalyRS < QSV[1])),] <- NA
Trichoptera[which((nd2_nontrop_Trichoptera$UI2=="Secondary vegetation") & (nd2_nontrop_Trichoptera$StdTmeanAnomalyRS > QSV[2])),] <- NA
Trichoptera[which((nd2_nontrop_Trichoptera$UI2=="Agriculture_Low") & (nd2_nontrop_Trichoptera$StdTmeanAnomalyRS < QAL[1])),] <- NA
Trichoptera[which((nd2_nontrop_Trichoptera$UI2=="Agriculture_Low") & (nd2_nontrop_Trichoptera$StdTmeanAnomalyRS > QAL[2])),] <- NA
Trichoptera[which((nd2_nontrop_Trichoptera$UI2=="Agriculture_High") & (nd2_nontrop_Trichoptera$StdTmeanAnomalyRS < QAH[1])),] <- NA
Trichoptera[which((nd2_nontrop_Trichoptera$UI2=="Agriculture_High") & (nd2_nontrop_Trichoptera$StdTmeanAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd2_nontrop_Blattodea$PredMedian <- ((apply(X = Blattodea,MARGIN = 1,
                                            FUN = median,na.rm=TRUE))*100)-100
nd2_nontrop_Blattodea$PredUpper <- ((apply(X = Blattodea,MARGIN = 1,
                                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_nontrop_Blattodea$PredLower <- ((apply(X = Blattodea,MARGIN = 1,
                                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

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

nd2_nontrop_Neuroptera$PredMedian <- ((apply(X = Neuroptera,MARGIN = 1,
                                             FUN = median,na.rm=TRUE))*100)-100
nd2_nontrop_Neuroptera$PredUpper <- ((apply(X = Neuroptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_nontrop_Neuroptera$PredLower <- ((apply(X = Neuroptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_nontrop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                             FUN = median,na.rm=TRUE))*100)-100
nd2_nontrop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_nontrop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_nontrop_Other$PredMedian <- ((apply(X = Other,MARGIN = 1,
                                        FUN = median,na.rm=TRUE))*100)-100
nd2_nontrop_Other$PredUpper <- ((apply(X = Other,MARGIN = 1,
                                       FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_nontrop_Other$PredLower <- ((apply(X = Other,MARGIN = 1,
                                       FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_nontrop_Thysanoptera$PredMedian <- ((apply(X = Thysanoptera,MARGIN = 1,
                                               FUN = median,na.rm=TRUE))*100)-100
nd2_nontrop_Thysanoptera$PredUpper <- ((apply(X = Thysanoptera,MARGIN = 1,
                                              FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_nontrop_Thysanoptera$PredLower <- ((apply(X = Thysanoptera,MARGIN = 1,
                                              FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2_nontrop_Trichoptera$PredMedian <- ((apply(X = Trichoptera,MARGIN = 1,
                                              FUN = median,na.rm=TRUE))*100)-100
nd2_nontrop_Trichoptera$PredUpper <- ((apply(X = Trichoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2_nontrop_Trichoptera$PredLower <- ((apply(X = Trichoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
# nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Blattodea$UI2 <- factor(nd2_nontrop_Blattodea$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Coleoptera$UI2 <- factor(nd2_nontrop_Coleoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Diptera$UI2 <- factor(nd2_nontrop_Diptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Hemiptera$UI2 <- factor(nd2_nontrop_Hemiptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Hymenoptera$UI2 <- factor(nd2_nontrop_Hymenoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Lepidoptera$UI2 <- factor(nd2_nontrop_Lepidoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Neuroptera$UI2 <- factor(nd2_nontrop_Neuroptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Orthoptera$UI2 <- factor(nd2_nontrop_Orthoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Other$UI2 <- factor(nd2_nontrop_Other$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Thysanoptera$UI2 <- factor(nd2_nontrop_Thysanoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd2_nontrop_Trichoptera$UI2 <- factor(nd2_nontrop_Trichoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# add zone factor
nd2_nontrop_Blattodea$Zone <- as.factor("NonTropical")
nd2_nontrop_Coleoptera$Zone <- as.factor("NonTropical")
nd2_nontrop_Diptera$Zone <- as.factor("NonTropical")
nd2_nontrop_Hemiptera$Zone <- as.factor("NonTropical")
nd2_nontrop_Hymenoptera$Zone <- as.factor("NonTropical")
nd2_nontrop_Lepidoptera$Zone <- as.factor("NonTropical")
nd2_nontrop_Neuroptera$Zone <- as.factor("NonTropical")
nd2_nontrop_Orthoptera$Zone <- as.factor("NonTropical")
nd2_nontrop_Other$Zone <- as.factor("NonTropical")
nd2_nontrop_Thysanoptera$Zone <- as.factor("NonTropical")
nd2_nontrop_Trichoptera$Zone <- as.factor("NonTropical")

# put use rbind to add nd2_nontrop to nd2_trop to make one data table for plotting

nd2_Blattodea <- rbind(nd2_trop_Blattodea,nd2_nontrop_Blattodea)
nd2_Coleoptera <- rbind(nd2_trop_Coleoptera,nd2_nontrop_Coleoptera)
nd2_Diptera <- rbind(nd2_trop_Diptera,nd2_nontrop_Diptera)
nd2_Hemiptera <- rbind(nd2_trop_Hemiptera,nd2_nontrop_Hemiptera)
nd2_Hymenoptera <- rbind(nd2_trop_Hymenoptera,nd2_nontrop_Hymenoptera)
nd2_Lepidoptera <- rbind(nd2_trop_Lepidoptera,nd2_nontrop_Lepidoptera)
nd2_Neuroptera <- rbind(nd2_trop_Neuroptera,nd2_nontrop_Neuroptera)
nd2_Orthoptera <- rbind(nd2_trop_Orthoptera,nd2_nontrop_Orthoptera)
nd2_Other <- rbind(nd2_trop_Other,nd2_nontrop_Other)
nd2_Thysanoptera <- rbind(nd2_trop_Thysanoptera,nd2_nontrop_Thysanoptera)
nd2_Trichoptera <- rbind(nd2_trop_Trichoptera,nd2_nontrop_Trichoptera)

# plot

p_blattodea <- ggplot(data = nd2_Blattodea, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Blattodea$PredLower, ymax = nd2_Blattodea$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) + 
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) + # extended y axis for orders with high larger confidence intervals
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        # legend.background = element_blank(), 
        # legend.text = element_text(size = 6), 
        # legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Blattodea")

p_coleoptera <- ggplot(data = nd2_Coleoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Coleoptera$PredLower, ymax = nd2_Coleoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Diptera$PredLower, ymax = nd2_Diptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Hemiptera$PredLower, ymax = nd2_Hemiptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Hymenoptera$PredLower, ymax = nd2_Hymenoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Lepidoptera$PredLower, ymax = nd2_Lepidoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_neuroptera <- ggplot(data = nd2_Neuroptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Neuroptera$PredLower, ymax = nd2_Neuroptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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
  ggtitle("Neuroptera")

p_other <- ggplot(data = nd2_Other, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Other$PredLower, ymax = nd2_Other$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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
  ggtitle("Other")

p_orthoptera <- ggplot(data = nd2_Orthoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Orthoptera$PredLower, ymax = nd2_Orthoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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

p_thysanoptera <- ggplot(data = nd2_Thysanoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Thysanoptera$PredLower, ymax = nd2_Thysanoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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
  ggtitle("Thysanoptera")

p_trichoptera <- ggplot(data = nd2_Trichoptera, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd2_Trichoptera$PredLower, ymax = nd2_Trichoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
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
  ggtitle("Trichoptera")

# get the legend
legend <- get_legend(
  p_blattodea +
    guides(color = guide_legend(nrow = 1),
           linetype = guide_legend (nrow=1)) +
    theme(legend.position = "right",
          legend.background = element_blank(), 
          legend.text = element_text(size = 6), 
          legend.title = element_blank())
)


# put them all together to save them
MeanAnomRich <- cowplot::plot_grid(p_blattodea, p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,p_neuroptera,p_orthoptera,p_other,p_thysanoptera,p_trichoptera,legend)

# save them
ggsave(filename = paste0(outDir, "MeanAnomRich.pdf"), plot = MeanAnomRich, width = 200, height = 150, units = "mm", dpi = 300)
#ggsave(filename = paste0(outDir, "MeanAnomRich_extended yaxis.pdf"), plot = MeanAnomRich, width = 200, height = 150, units = "mm", dpi = 300)

## Abundance, Max Anomaly ##
## Tropical ##

nd3_trop <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelAbund_trop$data$StdTmaxAnomalyRS),
                       to = max(MaxAnomalyModelAbund_trop$data$StdTmaxAnomalyRS),
                       length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelAbund_trop$data$UI2)),
  Order=factor(c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")))

# back transform the predictors
nd3_trop$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3_trop$StdTmaxAnomalyRS,
  originalX = trop$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd3_trop$LogAbund <- 0
nd3_trop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# does this for each Order
# need to record it to use later on
refRow <- which((nd3_trop$UI2=="Primary vegetation") & (nd3_trop$StdTmaxAnomaly==min(abs(nd3_trop$StdTmaxAnomaly))))
# don't know why this isn't woking - will work it out later, but the refRow is 55th row every 400 rows (I checked)

QPV <- quantile(x = MaxAnomalyModelAbund_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_trop$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAnomalyModelAbund_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_trop$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAnomalyModelAbund_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_trop$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAnomalyModelAbund_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_trop$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmax.trop <- PredictGLMERRandIter(model = MaxAnomalyModelAbund_trop$model,data = nd3_trop)

# back transform the abundance values
a.preds.tmax.trop <- exp(a.preds.tmax.trop)-0.01


# split up by order
number_of_chunks = 11
list_a.preds.tmax.trop <- lapply(seq(1, NROW(a.preds.tmax.trop), ceiling(NROW(a.preds.tmax.trop)/number_of_chunks)),
                                 function(i) a.preds.tmax.trop[i:min(i + ceiling(NROW(a.preds.tmax.trop)/number_of_chunks) - 1, NROW(a.preds.tmax.trop)),])
# success!
# creates list of matrices
# name them
names(list_a.preds.tmax.trop) <- c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")
list2env(list_a.preds.tmax.trop,globalenv())


# tim's suggestion
list_a.preds.tmax.trop <- lapply(list_a.preds.tmax.trop,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[55,],FUN="/") 
})

list2env(list_a.preds.tmax.trop,globalenv())

# works!

# Charlie's input: if there are a few, use facet.wrap, if there are more, use ggplot and then cowplot

# split nd3_trop by order
Order<- paste0("nd3_trop_",nd3_trop$Order)
# create a list of data frames
by_Order <- split(nd3_trop,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles

Blattodea[which((nd3_trop_Blattodea$UI2=="Primary vegetation") & (nd3_trop_Blattodea$StdTmaxAnomalyRS > QPV[2])),] <- NA
Blattodea[which((nd3_trop_Blattodea$UI2=="Primary vegetation") & (nd3_trop_Blattodea$StdTmaxAnomalyRS < QPV[1])),] <- NA
Blattodea[which((nd3_trop_Blattodea$UI2=="Secondary vegetation") & (nd3_trop_Blattodea$StdTmaxAnomalyRS < QSV[1])),] <- NA
Blattodea[which((nd3_trop_Blattodea$UI2=="Secondary vegetation") & (nd3_trop_Blattodea$StdTmaxAnomalyRS > QSV[2])),] <- NA
Blattodea[which((nd3_trop_Blattodea$UI2=="Agriculture_Low") & (nd3_trop_Blattodea$StdTmaxAnomalyRS < QAL[1])),] <- NA
Blattodea[which((nd3_trop_Blattodea$UI2=="Agriculture_Low") & (nd3_trop_Blattodea$StdTmaxAnomalyRS > QAL[2])),] <- NA
Blattodea[which((nd3_trop_Blattodea$UI2=="Agriculture_High") & (nd3_trop_Blattodea$StdTmaxAnomalyRS < QAH[1])),] <- NA
Blattodea[which((nd3_trop_Blattodea$UI2=="Agriculture_High") & (nd3_trop_Blattodea$StdTmaxAnomalyRS > QAH[2])),] <- NA

Coleoptera[which((nd3_trop_Coleoptera$UI2=="Primary vegetation") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd3_trop_Coleoptera$UI2=="Primary vegetation") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd3_trop_Coleoptera$UI2=="Secondary vegetation") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd3_trop_Coleoptera$UI2=="Secondary vegetation") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd3_trop_Coleoptera$UI2=="Agriculture_Low") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd3_trop_Coleoptera$UI2=="Agriculture_Low") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd3_trop_Coleoptera$UI2=="Agriculture_High") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd3_trop_Coleoptera$UI2=="Agriculture_High") & (nd3_trop_Coleoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd3_trop_Diptera$UI2=="Primary vegetation") & (nd3_trop_Diptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd3_trop_Diptera$UI2=="Primary vegetation") & (nd3_trop_Diptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd3_trop_Diptera$UI2=="Secondary vegetation") & (nd3_trop_Diptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd3_trop_Diptera$UI2=="Secondary vegetation") & (nd3_trop_Diptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd3_trop_Diptera$UI2=="Agriculture_Low") & (nd3_trop_Diptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd3_trop_Diptera$UI2=="Agriculture_Low") & (nd3_trop_Diptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd3_trop_Diptera$UI2=="Agriculture_High") & (nd3_trop_Diptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd3_trop_Diptera$UI2=="Agriculture_High") & (nd3_trop_Diptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd3_trop_Hemiptera$UI2=="Primary vegetation") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd3_trop_Hemiptera$UI2=="Primary vegetation") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd3_trop_Hemiptera$UI2=="Secondary vegetation") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd3_trop_Hemiptera$UI2=="Secondary vegetation") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd3_trop_Hemiptera$UI2=="Agriculture_Low") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd3_trop_Hemiptera$UI2=="Agriculture_Low") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd3_trop_Hemiptera$UI2=="Agriculture_High") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd3_trop_Hemiptera$UI2=="Agriculture_High") & (nd3_trop_Hemiptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd3_trop_Hymenoptera$UI2=="Primary vegetation") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd3_trop_Hymenoptera$UI2=="Primary vegetation") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd3_trop_Hymenoptera$UI2=="Secondary vegetation") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd3_trop_Hymenoptera$UI2=="Secondary vegetation") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd3_trop_Hymenoptera$UI2=="Agriculture_Low") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd3_trop_Hymenoptera$UI2=="Agriculture_Low") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd3_trop_Hymenoptera$UI2=="Agriculture_High") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd3_trop_Hymenoptera$UI2=="Agriculture_High") & (nd3_trop_Hymenoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd3_trop_Lepidoptera$UI2=="Primary vegetation") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd3_trop_Lepidoptera$UI2=="Primary vegetation") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd3_trop_Lepidoptera$UI2=="Secondary vegetation") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd3_trop_Lepidoptera$UI2=="Secondary vegetation") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd3_trop_Lepidoptera$UI2=="Agriculture_Low") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd3_trop_Lepidoptera$UI2=="Agriculture_Low") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd3_trop_Lepidoptera$UI2=="Agriculture_High") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd3_trop_Lepidoptera$UI2=="Agriculture_High") & (nd3_trop_Lepidoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Neuroptera[which((nd3_trop_Neuroptera$UI2=="Primary vegetation") & (nd3_trop_Neuroptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Neuroptera[which((nd3_trop_Neuroptera$UI2=="Primary vegetation") & (nd3_trop_Neuroptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Neuroptera[which((nd3_trop_Neuroptera$UI2=="Secondary vegetation") & (nd3_trop_Neuroptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Neuroptera[which((nd3_trop_Neuroptera$UI2=="Secondary vegetation") & (nd3_trop_Neuroptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Neuroptera[which((nd3_trop_Neuroptera$UI2=="Agriculture_Low") & (nd3_trop_Neuroptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Neuroptera[which((nd3_trop_Neuroptera$UI2=="Agriculture_Low") & (nd3_trop_Neuroptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Neuroptera[which((nd3_trop_Neuroptera$UI2=="Agriculture_High") & (nd3_trop_Neuroptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Neuroptera[which((nd3_trop_Neuroptera$UI2=="Agriculture_High") & (nd3_trop_Neuroptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd3_trop_Orthoptera$UI2=="Primary vegetation") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd3_trop_Orthoptera$UI2=="Primary vegetation") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd3_trop_Orthoptera$UI2=="Secondary vegetation") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd3_trop_Orthoptera$UI2=="Secondary vegetation") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd3_trop_Orthoptera$UI2=="Agriculture_Low") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd3_trop_Orthoptera$UI2=="Agriculture_Low") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd3_trop_Orthoptera$UI2=="Agriculture_High") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd3_trop_Orthoptera$UI2=="Agriculture_High") & (nd3_trop_Orthoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Other[which((nd3_trop_Other$UI2=="Primary vegetation") & (nd3_trop_Other$StdTmaxAnomalyRS > QPV[2])),] <- NA
Other[which((nd3_trop_Other$UI2=="Primary vegetation") & (nd3_trop_Other$StdTmaxAnomalyRS < QPV[1])),] <- NA
Other[which((nd3_trop_Other$UI2=="Secondary vegetation") & (nd3_trop_Other$StdTmaxAnomalyRS < QSV[1])),] <- NA
Other[which((nd3_trop_Other$UI2=="Secondary vegetation") & (nd3_trop_Other$StdTmaxAnomalyRS > QSV[2])),] <- NA
Other[which((nd3_trop_Other$UI2=="Agriculture_Low") & (nd3_trop_Other$StdTmaxAnomalyRS < QAL[1])),] <- NA
Other[which((nd3_trop_Other$UI2=="Agriculture_Low") & (nd3_trop_Other$StdTmaxAnomalyRS > QAL[2])),] <- NA
Other[which((nd3_trop_Other$UI2=="Agriculture_High") & (nd3_trop_Other$StdTmaxAnomalyRS < QAH[1])),] <- NA
Other[which((nd3_trop_Other$UI2=="Agriculture_High") & (nd3_trop_Other$StdTmaxAnomalyRS > QAH[2])),] <- NA

Thysanoptera[which((nd3_trop_Thysanoptera$UI2=="Primary vegetation") & (nd3_trop_Thysanoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Thysanoptera[which((nd3_trop_Thysanoptera$UI2=="Primary vegetation") & (nd3_trop_Thysanoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Thysanoptera[which((nd3_trop_Thysanoptera$UI2=="Secondary vegetation") & (nd3_trop_Thysanoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Thysanoptera[which((nd3_trop_Thysanoptera$UI2=="Secondary vegetation") & (nd3_trop_Thysanoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Thysanoptera[which((nd3_trop_Thysanoptera$UI2=="Agriculture_Low") & (nd3_trop_Thysanoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Thysanoptera[which((nd3_trop_Thysanoptera$UI2=="Agriculture_Low") & (nd3_trop_Thysanoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Thysanoptera[which((nd3_trop_Thysanoptera$UI2=="Agriculture_High") & (nd3_trop_Thysanoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Thysanoptera[which((nd3_trop_Thysanoptera$UI2=="Agriculture_High") & (nd3_trop_Thysanoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Trichoptera[which((nd3_trop_Trichoptera$UI2=="Primary vegetation") & (nd3_trop_Trichoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Trichoptera[which((nd3_trop_Trichoptera$UI2=="Primary vegetation") & (nd3_trop_Trichoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Trichoptera[which((nd3_trop_Trichoptera$UI2=="Secondary vegetation") & (nd3_trop_Trichoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Trichoptera[which((nd3_trop_Trichoptera$UI2=="Secondary vegetation") & (nd3_trop_Trichoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Trichoptera[which((nd3_trop_Trichoptera$UI2=="Agriculture_Low") & (nd3_trop_Trichoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Trichoptera[which((nd3_trop_Trichoptera$UI2=="Agriculture_Low") & (nd3_trop_Trichoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Trichoptera[which((nd3_trop_Trichoptera$UI2=="Agriculture_High") & (nd3_trop_Trichoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Trichoptera[which((nd3_trop_Trichoptera$UI2=="Agriculture_High") & (nd3_trop_Trichoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd3_trop_Blattodea$PredMedian <- ((apply(X = Blattodea,MARGIN = 1,
                                         FUN = median,na.rm=TRUE))*100)-100
nd3_trop_Blattodea$PredUpper <- ((apply(X = Blattodea,MARGIN = 1,
                                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_trop_Blattodea$PredLower <- ((apply(X = Blattodea,MARGIN = 1,
                                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

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

nd3_trop_Neuroptera$PredMedian <- ((apply(X = Neuroptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd3_trop_Neuroptera$PredUpper <- ((apply(X = Neuroptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_trop_Neuroptera$PredLower <- ((apply(X = Neuroptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_trop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd3_trop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_trop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_trop_Other$PredMedian <- ((apply(X = Other,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd3_trop_Other$PredUpper <- ((apply(X = Other,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_trop_Other$PredLower <- ((apply(X = Other,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_trop_Thysanoptera$PredMedian <- ((apply(X = Thysanoptera,MARGIN = 1,
                                            FUN = median,na.rm=TRUE))*100)-100
nd3_trop_Thysanoptera$PredUpper <- ((apply(X = Thysanoptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_trop_Thysanoptera$PredLower <- ((apply(X = Thysanoptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_trop_Trichoptera$PredMedian <- ((apply(X = Trichoptera,MARGIN = 1,
                                           FUN = median,na.rm=TRUE))*100)-100
nd3_trop_Trichoptera$PredUpper <- ((apply(X = Trichoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_trop_Trichoptera$PredLower <- ((apply(X = Trichoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
# nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Blattodea$UI2 <- factor(nd3_trop_Blattodea$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Coleoptera$UI2 <- factor(nd3_trop_Coleoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Diptera$UI2 <- factor(nd3_trop_Diptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Hemiptera$UI2 <- factor(nd3_trop_Hemiptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Hymenoptera$UI2 <- factor(nd3_trop_Hymenoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Lepidoptera$UI2 <- factor(nd3_trop_Lepidoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Neuroptera$UI2 <- factor(nd3_trop_Neuroptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Orthoptera$UI2 <- factor(nd3_trop_Orthoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Other$UI2 <- factor(nd3_trop_Other$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Thysanoptera$UI2 <- factor(nd3_trop_Thysanoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_trop_Trichoptera$UI2 <- factor(nd3_trop_Trichoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# add zone factor
nd3_trop_Blattodea$Zone <- as.factor("Tropical")
nd3_trop_Coleoptera$Zone <- as.factor("Tropical")
nd3_trop_Diptera$Zone <- as.factor("Tropical")
nd3_trop_Hemiptera$Zone <- as.factor("Tropical")
nd3_trop_Hymenoptera$Zone <- as.factor("Tropical")
nd3_trop_Lepidoptera$Zone <- as.factor("Tropical")
nd3_trop_Neuroptera$Zone <- as.factor("Tropical")
nd3_trop_Orthoptera$Zone <- as.factor("Tropical")
nd3_trop_Other$Zone <- as.factor("Tropical")
nd3_trop_Thysanoptera$Zone <- as.factor("Tropical")
nd3_trop_Trichoptera$Zone <- as.factor("Tropical")

## NonTropical ##

nd3_nontrop <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelAbund_nontrop$data$StdTmaxAnomalyRS),
                       to = max(MaxAnomalyModelAbund_nontrop$data$StdTmaxAnomalyRS),
                       length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelAbund_nontrop$data$UI2)),
  Order=factor(c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")))

# back transform the predictors
nd3_nontrop$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3_nontrop$StdTmaxAnomalyRS,
  originalX = nontrop$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd3_nontrop$LogAbund <- 0
nd3_nontrop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# does this for each Order
# reference row is 4th row, every 400 rows (see 'Values')
refRow <- which((nd3_nontrop$UI2=="Primary vegetation") & (nd3_nontrop$StdTmaxAnomaly==min(abs(nd3_nontrop$StdTmaxAnomaly))))
# 15th row, every 400 rows

# adjust plot 1: max anomaly and abundance

QPV <- quantile(x = MaxAnomalyModelAbund_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_nontrop$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAnomalyModelAbund_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_nontrop$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAnomalyModelAbund_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_nontrop$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAnomalyModelAbund_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund_nontrop$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmax.nontrop <- PredictGLMERRandIter(model = MaxAnomalyModelAbund_nontrop$model,data = nd3_nontrop)

# back transform the abundance values
a.preds.tmax.nontrop <- exp(a.preds.tmax.nontrop)-0.01


# another try!
number_of_chunks = 11
list_a.preds.tmax.nontrop <- lapply(seq(1, NROW(a.preds.tmax.nontrop), ceiling(NROW(a.preds.tmax.nontrop)/number_of_chunks)),
                                    function(i) a.preds.tmax.nontrop[i:min(i + ceiling(NROW(a.preds.tmax.nontrop)/number_of_chunks) - 1, NROW(a.preds.tmax.nontrop)),])
# success!
# creates list of matrices
# name them
names(list_a.preds.tmax.nontrop) <- c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")
list2env(list_a.preds.tmax.nontrop,globalenv())

# tim's suggestion
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

Blattodea[which((nd3_nontrop_Blattodea$UI2=="Primary vegetation") & (nd3_nontrop_Blattodea$StdTmaxAnomalyRS > QPV[2])),] <- NA
Blattodea[which((nd3_nontrop_Blattodea$UI2=="Primary vegetation") & (nd3_nontrop_Blattodea$StdTmaxAnomalyRS < QPV[1])),] <- NA
Blattodea[which((nd3_nontrop_Blattodea$UI2=="Secondary vegetation") & (nd3_nontrop_Blattodea$StdTmaxAnomalyRS < QSV[1])),] <- NA
Blattodea[which((nd3_nontrop_Blattodea$UI2=="Secondary vegetation") & (nd3_nontrop_Blattodea$StdTmaxAnomalyRS > QSV[2])),] <- NA
Blattodea[which((nd3_nontrop_Blattodea$UI2=="Agriculture_Low") & (nd3_nontrop_Blattodea$StdTmaxAnomalyRS < QAL[1])),] <- NA
Blattodea[which((nd3_nontrop_Blattodea$UI2=="Agriculture_Low") & (nd3_nontrop_Blattodea$StdTmaxAnomalyRS > QAL[2])),] <- NA
Blattodea[which((nd3_nontrop_Blattodea$UI2=="Agriculture_High") & (nd3_nontrop_Blattodea$StdTmaxAnomalyRS < QAH[1])),] <- NA
Blattodea[which((nd3_nontrop_Blattodea$UI2=="Agriculture_High") & (nd3_nontrop_Blattodea$StdTmaxAnomalyRS > QAH[2])),] <- NA

Coleoptera[which((nd3_nontrop_Coleoptera$UI2=="Primary vegetation") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd3_nontrop_Coleoptera$UI2=="Primary vegetation") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd3_nontrop_Coleoptera$UI2=="Secondary vegetation") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd3_nontrop_Coleoptera$UI2=="Secondary vegetation") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd3_nontrop_Coleoptera$UI2=="Agriculture_Low") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd3_nontrop_Coleoptera$UI2=="Agriculture_Low") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd3_nontrop_Coleoptera$UI2=="Agriculture_High") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd3_nontrop_Coleoptera$UI2=="Agriculture_High") & (nd3_nontrop_Coleoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd3_nontrop_Diptera$UI2=="Primary vegetation") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd3_nontrop_Diptera$UI2=="Primary vegetation") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd3_nontrop_Diptera$UI2=="Secondary vegetation") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd3_nontrop_Diptera$UI2=="Secondary vegetation") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd3_nontrop_Diptera$UI2=="Agriculture_Low") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd3_nontrop_Diptera$UI2=="Agriculture_Low") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd3_nontrop_Diptera$UI2=="Agriculture_High") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd3_nontrop_Diptera$UI2=="Agriculture_High") & (nd3_nontrop_Diptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd3_nontrop_Hemiptera$UI2=="Primary vegetation") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd3_nontrop_Hemiptera$UI2=="Primary vegetation") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd3_nontrop_Hemiptera$UI2=="Secondary vegetation") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd3_nontrop_Hemiptera$UI2=="Secondary vegetation") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd3_nontrop_Hemiptera$UI2=="Agriculture_Low") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd3_nontrop_Hemiptera$UI2=="Agriculture_Low") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd3_nontrop_Hemiptera$UI2=="Agriculture_High") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd3_nontrop_Hemiptera$UI2=="Agriculture_High") & (nd3_nontrop_Hemiptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd3_nontrop_Hymenoptera$UI2=="Primary vegetation") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd3_nontrop_Hymenoptera$UI2=="Primary vegetation") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd3_nontrop_Hymenoptera$UI2=="Secondary vegetation") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd3_nontrop_Hymenoptera$UI2=="Secondary vegetation") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd3_nontrop_Hymenoptera$UI2=="Agriculture_Low") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd3_nontrop_Hymenoptera$UI2=="Agriculture_Low") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd3_nontrop_Hymenoptera$UI2=="Agriculture_High") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd3_nontrop_Hymenoptera$UI2=="Agriculture_High") & (nd3_nontrop_Hymenoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd3_nontrop_Lepidoptera$UI2=="Primary vegetation") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd3_nontrop_Lepidoptera$UI2=="Primary vegetation") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd3_nontrop_Lepidoptera$UI2=="Secondary vegetation") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd3_nontrop_Lepidoptera$UI2=="Secondary vegetation") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd3_nontrop_Lepidoptera$UI2=="Agriculture_Low") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd3_nontrop_Lepidoptera$UI2=="Agriculture_Low") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd3_nontrop_Lepidoptera$UI2=="Agriculture_High") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd3_nontrop_Lepidoptera$UI2=="Agriculture_High") & (nd3_nontrop_Lepidoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Neuroptera[which((nd3_nontrop_Neuroptera$UI2=="Primary vegetation") & (nd3_nontrop_Neuroptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Neuroptera[which((nd3_nontrop_Neuroptera$UI2=="Primary vegetation") & (nd3_nontrop_Neuroptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Neuroptera[which((nd3_nontrop_Neuroptera$UI2=="Secondary vegetation") & (nd3_nontrop_Neuroptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Neuroptera[which((nd3_nontrop_Neuroptera$UI2=="Secondary vegetation") & (nd3_nontrop_Neuroptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Neuroptera[which((nd3_nontrop_Neuroptera$UI2=="Agriculture_Low") & (nd3_nontrop_Neuroptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Neuroptera[which((nd3_nontrop_Neuroptera$UI2=="Agriculture_Low") & (nd3_nontrop_Neuroptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Neuroptera[which((nd3_nontrop_Neuroptera$UI2=="Agriculture_High") & (nd3_nontrop_Neuroptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Neuroptera[which((nd3_nontrop_Neuroptera$UI2=="Agriculture_High") & (nd3_nontrop_Neuroptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd3_nontrop_Orthoptera$UI2=="Primary vegetation") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd3_nontrop_Orthoptera$UI2=="Primary vegetation") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd3_nontrop_Orthoptera$UI2=="Secondary vegetation") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd3_nontrop_Orthoptera$UI2=="Secondary vegetation") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd3_nontrop_Orthoptera$UI2=="Agriculture_Low") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd3_nontrop_Orthoptera$UI2=="Agriculture_Low") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd3_nontrop_Orthoptera$UI2=="Agriculture_High") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd3_nontrop_Orthoptera$UI2=="Agriculture_High") & (nd3_nontrop_Orthoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Other[which((nd3_nontrop_Other$UI2=="Primary vegetation") & (nd3_nontrop_Other$StdTmaxAnomalyRS > QPV[2])),] <- NA
Other[which((nd3_nontrop_Other$UI2=="Primary vegetation") & (nd3_nontrop_Other$StdTmaxAnomalyRS < QPV[1])),] <- NA
Other[which((nd3_nontrop_Other$UI2=="Secondary vegetation") & (nd3_nontrop_Other$StdTmaxAnomalyRS < QSV[1])),] <- NA
Other[which((nd3_nontrop_Other$UI2=="Secondary vegetation") & (nd3_nontrop_Other$StdTmaxAnomalyRS > QSV[2])),] <- NA
Other[which((nd3_nontrop_Other$UI2=="Agriculture_Low") & (nd3_nontrop_Other$StdTmaxAnomalyRS < QAL[1])),] <- NA
Other[which((nd3_nontrop_Other$UI2=="Agriculture_Low") & (nd3_nontrop_Other$StdTmaxAnomalyRS > QAL[2])),] <- NA
Other[which((nd3_nontrop_Other$UI2=="Agriculture_High") & (nd3_nontrop_Other$StdTmaxAnomalyRS < QAH[1])),] <- NA
Other[which((nd3_nontrop_Other$UI2=="Agriculture_High") & (nd3_nontrop_Other$StdTmaxAnomalyRS > QAH[2])),] <- NA

Thysanoptera[which((nd3_nontrop_Thysanoptera$UI2=="Primary vegetation") & (nd3_nontrop_Thysanoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Thysanoptera[which((nd3_nontrop_Thysanoptera$UI2=="Primary vegetation") & (nd3_nontrop_Thysanoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Thysanoptera[which((nd3_nontrop_Thysanoptera$UI2=="Secondary vegetation") & (nd3_nontrop_Thysanoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Thysanoptera[which((nd3_nontrop_Thysanoptera$UI2=="Secondary vegetation") & (nd3_nontrop_Thysanoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Thysanoptera[which((nd3_nontrop_Thysanoptera$UI2=="Agriculture_Low") & (nd3_nontrop_Thysanoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Thysanoptera[which((nd3_nontrop_Thysanoptera$UI2=="Agriculture_Low") & (nd3_nontrop_Thysanoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Thysanoptera[which((nd3_nontrop_Thysanoptera$UI2=="Agriculture_High") & (nd3_nontrop_Thysanoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Thysanoptera[which((nd3_nontrop_Thysanoptera$UI2=="Agriculture_High") & (nd3_nontrop_Thysanoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Trichoptera[which((nd3_nontrop_Trichoptera$UI2=="Primary vegetation") & (nd3_nontrop_Trichoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Trichoptera[which((nd3_nontrop_Trichoptera$UI2=="Primary vegetation") & (nd3_nontrop_Trichoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Trichoptera[which((nd3_nontrop_Trichoptera$UI2=="Secondary vegetation") & (nd3_nontrop_Trichoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Trichoptera[which((nd3_nontrop_Trichoptera$UI2=="Secondary vegetation") & (nd3_nontrop_Trichoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Trichoptera[which((nd3_nontrop_Trichoptera$UI2=="Agriculture_Low") & (nd3_nontrop_Trichoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Trichoptera[which((nd3_nontrop_Trichoptera$UI2=="Agriculture_Low") & (nd3_nontrop_Trichoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Trichoptera[which((nd3_nontrop_Trichoptera$UI2=="Agriculture_High") & (nd3_nontrop_Trichoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Trichoptera[which((nd3_nontrop_Trichoptera$UI2=="Agriculture_High") & (nd3_nontrop_Trichoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd3_nontrop_Blattodea$PredMedian <- ((apply(X = Blattodea,MARGIN = 1,
                                            FUN = median,na.rm=TRUE))*100)-100
nd3_nontrop_Blattodea$PredUpper <- ((apply(X = Blattodea,MARGIN = 1,
                                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_nontrop_Blattodea$PredLower <- ((apply(X = Blattodea,MARGIN = 1,
                                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

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

nd3_nontrop_Neuroptera$PredMedian <- ((apply(X = Neuroptera,MARGIN = 1,
                                             FUN = median,na.rm=TRUE))*100)-100
nd3_nontrop_Neuroptera$PredUpper <- ((apply(X = Neuroptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_nontrop_Neuroptera$PredLower <- ((apply(X = Neuroptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_nontrop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                             FUN = median,na.rm=TRUE))*100)-100
nd3_nontrop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_nontrop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_nontrop_Other$PredMedian <- ((apply(X = Other,MARGIN = 1,
                                        FUN = median,na.rm=TRUE))*100)-100
nd3_nontrop_Other$PredUpper <- ((apply(X = Other,MARGIN = 1,
                                       FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_nontrop_Other$PredLower <- ((apply(X = Other,MARGIN = 1,
                                       FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_nontrop_Thysanoptera$PredMedian <- ((apply(X = Thysanoptera,MARGIN = 1,
                                               FUN = median,na.rm=TRUE))*100)-100
nd3_nontrop_Thysanoptera$PredUpper <- ((apply(X = Thysanoptera,MARGIN = 1,
                                              FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_nontrop_Thysanoptera$PredLower <- ((apply(X = Thysanoptera,MARGIN = 1,
                                              FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3_nontrop_Trichoptera$PredMedian <- ((apply(X = Trichoptera,MARGIN = 1,
                                              FUN = median,na.rm=TRUE))*100)-100
nd3_nontrop_Trichoptera$PredUpper <- ((apply(X = Trichoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3_nontrop_Trichoptera$PredLower <- ((apply(X = Trichoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
# nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Blattodea$UI2 <- factor(nd3_nontrop_Blattodea$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Coleoptera$UI2 <- factor(nd3_nontrop_Coleoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Diptera$UI2 <- factor(nd3_nontrop_Diptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Hemiptera$UI2 <- factor(nd3_nontrop_Hemiptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Hymenoptera$UI2 <- factor(nd3_nontrop_Hymenoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Lepidoptera$UI2 <- factor(nd3_nontrop_Lepidoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Neuroptera$UI2 <- factor(nd3_nontrop_Neuroptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Orthoptera$UI2 <- factor(nd3_nontrop_Orthoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Other$UI2 <- factor(nd3_nontrop_Other$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Thysanoptera$UI2 <- factor(nd3_nontrop_Thysanoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd3_nontrop_Trichoptera$UI2 <- factor(nd3_nontrop_Trichoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# add zone factor
nd3_nontrop_Blattodea$Zone <- as.factor("NonTropical")
nd3_nontrop_Coleoptera$Zone <- as.factor("NonTropical")
nd3_nontrop_Diptera$Zone <- as.factor("NonTropical")
nd3_nontrop_Hemiptera$Zone <- as.factor("NonTropical")
nd3_nontrop_Hymenoptera$Zone <- as.factor("NonTropical")
nd3_nontrop_Lepidoptera$Zone <- as.factor("NonTropical")
nd3_nontrop_Neuroptera$Zone <- as.factor("NonTropical")
nd3_nontrop_Orthoptera$Zone <- as.factor("NonTropical")
nd3_nontrop_Other$Zone <- as.factor("NonTropical")
nd3_nontrop_Thysanoptera$Zone <- as.factor("NonTropical")
nd3_nontrop_Trichoptera$Zone <- as.factor("NonTropical")

# put use rbind to add nd3_nontrop to nd3_trop to make one data table for plotting

nd3_Blattodea <- rbind(nd3_trop_Blattodea,nd3_nontrop_Blattodea)
nd3_Coleoptera <- rbind(nd3_trop_Coleoptera,nd3_nontrop_Coleoptera)
nd3_Diptera <- rbind(nd3_trop_Diptera,nd3_nontrop_Diptera)
nd3_Hemiptera <- rbind(nd3_trop_Hemiptera,nd3_nontrop_Hemiptera)
nd3_Hymenoptera <- rbind(nd3_trop_Hymenoptera,nd3_nontrop_Hymenoptera)
nd3_Lepidoptera <- rbind(nd3_trop_Lepidoptera,nd3_nontrop_Lepidoptera)
nd3_Neuroptera <- rbind(nd3_trop_Neuroptera,nd3_nontrop_Neuroptera)
nd3_Orthoptera <- rbind(nd3_trop_Orthoptera,nd3_nontrop_Orthoptera)
nd3_Other <- rbind(nd3_trop_Other,nd3_nontrop_Other)
nd3_Thysanoptera <- rbind(nd3_trop_Thysanoptera,nd3_nontrop_Thysanoptera)
nd3_Trichoptera <- rbind(nd3_trop_Trichoptera,nd3_nontrop_Trichoptera)

# plot

p_blattodea <- ggplot(data = nd3_Blattodea, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Blattodea$PredLower, ymax = nd3_Blattodea$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) + 
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) + # extended y axis for orders with high larger confidence intervals
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        # legend.background = element_blank(), 
        # legend.text = element_text(size = 6), 
        # legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Blattodea")

p_coleoptera <- ggplot(data = nd3_Coleoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Coleoptera$PredLower, ymax = nd3_Coleoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Diptera$PredLower, ymax = nd3_Diptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Hemiptera$PredLower, ymax = nd3_Hemiptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Hymenoptera$PredLower, ymax = nd3_Hymenoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Lepidoptera$PredLower, ymax = nd3_Lepidoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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

p_neuroptera <- ggplot(data = nd3_Neuroptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Neuroptera$PredLower, ymax = nd3_Neuroptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  ggtitle("Neuroptera")

p_other <- ggplot(data = nd3_Other, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Other$PredLower, ymax = nd3_Other$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  ggtitle("Other")

p_orthoptera <- ggplot(data = nd3_Orthoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Orthoptera$PredLower, ymax = nd3_Orthoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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

p_thysanoptera <- ggplot(data = nd3_Thysanoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Thysanoptera$PredLower, ymax = nd3_Thysanoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  ggtitle("Thysanoptera")

p_trichoptera <- ggplot(data = nd3_Trichoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd3_Trichoptera$PredLower, ymax = nd3_Trichoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  ggtitle("Trichoptera")

# get the legend
legend <- get_legend(
  p_blattodea +
    guides(color = guide_legend(nrow = 1),
           linetype = guide_legend (nrow=1)) +
    theme(legend.position = "right",
          legend.background = element_blank(), 
          legend.text = element_text(size = 6), 
          legend.title = element_blank())
)


# put them all together to save them
MaxAnomAbund <- cowplot::plot_grid(p_blattodea, p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,p_neuroptera,p_orthoptera,p_other,p_thysanoptera,p_trichoptera,legend)

# save them
ggsave(filename = paste0(outDir, "MaxAnomAbund.pdf"), plot = MaxAnomAbund, width = 200, height = 150, units = "mm", dpi = 300)
#ggsave(filename = paste0(outDir, "MaxAnomAbund_extended yaxis.pdf"), plot = MaxAnomAbund, width = 200, height = 150, units = "mm", dpi = 300)

## Richness, Max Anomaly ##
## Tropical ##

nd4_trop <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelRich_trop$data$StdTmaxAnomalyRS),
                       to = max(MaxAnomalyModelRich_trop$data$StdTmaxAnomalyRS),
                       length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelRich_trop$data$UI2)),
  Order=factor(c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")))

# back transform the predictors
nd4_trop$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4_trop$StdTmaxAnomalyRS,
  originalX = trop$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd4_trop$LogRich <- 0
nd4_trop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# does this for each Order
# need to record it to use later on
refRow <- which((nd4_trop$UI2=="Primary vegetation") & (nd4_trop$StdTmaxAnomaly==min(abs(nd4_trop$StdTmaxAnomaly))))
# don't know why this isn't woking - will work it out later, but the refRow is 55th row every 400 rows (I checked)

QPV <- quantile(x = MaxAnomalyModelRich_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_trop$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAnomalyModelRich_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_trop$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAnomalyModelRich_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_trop$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAnomalyModelRich_trop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_trop$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
sr.preds.tmax.trop <- PredictGLMERRandIter(model = MaxAnomalyModelRich_trop$model,data = nd4_trop)

# back transform the abundance values
sr.preds.tmax.trop <- exp(sr.preds.tmax.trop)-0.01


# split up by order
number_of_chunks = 11
list_sr.preds.tmax.trop <- lapply(seq(1, NROW(sr.preds.tmax.trop), ceiling(NROW(sr.preds.tmax.trop)/number_of_chunks)),
                                  function(i) sr.preds.tmax.trop[i:min(i + ceiling(NROW(sr.preds.tmax.trop)/number_of_chunks) - 1, NROW(sr.preds.tmax.trop)),])
# success!
# creates list of matrices
# name them
names(list_sr.preds.tmax.trop) <- c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")
list2env(list_sr.preds.tmax.trop,globalenv())


# tim's suggestion
list_sr.preds.tmax.trop <- lapply(list_sr.preds.tmax.trop,FUN=function(x){
  sweep (x=x, MARGIN = 2, STATS=x[55,],FUN="/") 
})

list2env(list_sr.preds.tmax.trop,globalenv())

# split nd4_trop by order
Order<- paste0("nd4_trop_",nd4_trop$Order)
# create a list of data frames
by_Order <- split(nd4_trop,Order)
list2env(by_Order,globalenv())

# remove anything above and below the quantiles

Blattodea[which((nd4_trop_Blattodea$UI2=="Primary vegetation") & (nd4_trop_Blattodea$StdTmaxAnomalyRS > QPV[2])),] <- NA
Blattodea[which((nd4_trop_Blattodea$UI2=="Primary vegetation") & (nd4_trop_Blattodea$StdTmaxAnomalyRS < QPV[1])),] <- NA
Blattodea[which((nd4_trop_Blattodea$UI2=="Secondary vegetation") & (nd4_trop_Blattodea$StdTmaxAnomalyRS < QSV[1])),] <- NA
Blattodea[which((nd4_trop_Blattodea$UI2=="Secondary vegetation") & (nd4_trop_Blattodea$StdTmaxAnomalyRS > QSV[2])),] <- NA
Blattodea[which((nd4_trop_Blattodea$UI2=="Agriculture_Low") & (nd4_trop_Blattodea$StdTmaxAnomalyRS < QAL[1])),] <- NA
Blattodea[which((nd4_trop_Blattodea$UI2=="Agriculture_Low") & (nd4_trop_Blattodea$StdTmaxAnomalyRS > QAL[2])),] <- NA
Blattodea[which((nd4_trop_Blattodea$UI2=="Agriculture_High") & (nd4_trop_Blattodea$StdTmaxAnomalyRS < QAH[1])),] <- NA
Blattodea[which((nd4_trop_Blattodea$UI2=="Agriculture_High") & (nd4_trop_Blattodea$StdTmaxAnomalyRS > QAH[2])),] <- NA

Coleoptera[which((nd4_trop_Coleoptera$UI2=="Primary vegetation") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd4_trop_Coleoptera$UI2=="Primary vegetation") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd4_trop_Coleoptera$UI2=="Secondary vegetation") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd4_trop_Coleoptera$UI2=="Secondary vegetation") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd4_trop_Coleoptera$UI2=="Agriculture_Low") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd4_trop_Coleoptera$UI2=="Agriculture_Low") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd4_trop_Coleoptera$UI2=="Agriculture_High") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd4_trop_Coleoptera$UI2=="Agriculture_High") & (nd4_trop_Coleoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd4_trop_Diptera$UI2=="Primary vegetation") & (nd4_trop_Diptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd4_trop_Diptera$UI2=="Primary vegetation") & (nd4_trop_Diptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd4_trop_Diptera$UI2=="Secondary vegetation") & (nd4_trop_Diptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd4_trop_Diptera$UI2=="Secondary vegetation") & (nd4_trop_Diptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd4_trop_Diptera$UI2=="Agriculture_Low") & (nd4_trop_Diptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd4_trop_Diptera$UI2=="Agriculture_Low") & (nd4_trop_Diptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd4_trop_Diptera$UI2=="Agriculture_High") & (nd4_trop_Diptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd4_trop_Diptera$UI2=="Agriculture_High") & (nd4_trop_Diptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd4_trop_Hemiptera$UI2=="Primary vegetation") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd4_trop_Hemiptera$UI2=="Primary vegetation") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd4_trop_Hemiptera$UI2=="Secondary vegetation") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd4_trop_Hemiptera$UI2=="Secondary vegetation") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd4_trop_Hemiptera$UI2=="Agriculture_Low") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd4_trop_Hemiptera$UI2=="Agriculture_Low") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd4_trop_Hemiptera$UI2=="Agriculture_High") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd4_trop_Hemiptera$UI2=="Agriculture_High") & (nd4_trop_Hemiptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd4_trop_Hymenoptera$UI2=="Primary vegetation") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd4_trop_Hymenoptera$UI2=="Primary vegetation") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd4_trop_Hymenoptera$UI2=="Secondary vegetation") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd4_trop_Hymenoptera$UI2=="Secondary vegetation") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd4_trop_Hymenoptera$UI2=="Agriculture_Low") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd4_trop_Hymenoptera$UI2=="Agriculture_Low") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd4_trop_Hymenoptera$UI2=="Agriculture_High") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd4_trop_Hymenoptera$UI2=="Agriculture_High") & (nd4_trop_Hymenoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd4_trop_Lepidoptera$UI2=="Primary vegetation") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd4_trop_Lepidoptera$UI2=="Primary vegetation") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd4_trop_Lepidoptera$UI2=="Secondary vegetation") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd4_trop_Lepidoptera$UI2=="Secondary vegetation") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd4_trop_Lepidoptera$UI2=="Agriculture_Low") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd4_trop_Lepidoptera$UI2=="Agriculture_Low") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd4_trop_Lepidoptera$UI2=="Agriculture_High") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd4_trop_Lepidoptera$UI2=="Agriculture_High") & (nd4_trop_Lepidoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Neuroptera[which((nd4_trop_Neuroptera$UI2=="Primary vegetation") & (nd4_trop_Neuroptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Neuroptera[which((nd4_trop_Neuroptera$UI2=="Primary vegetation") & (nd4_trop_Neuroptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Neuroptera[which((nd4_trop_Neuroptera$UI2=="Secondary vegetation") & (nd4_trop_Neuroptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Neuroptera[which((nd4_trop_Neuroptera$UI2=="Secondary vegetation") & (nd4_trop_Neuroptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Neuroptera[which((nd4_trop_Neuroptera$UI2=="Agriculture_Low") & (nd4_trop_Neuroptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Neuroptera[which((nd4_trop_Neuroptera$UI2=="Agriculture_Low") & (nd4_trop_Neuroptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Neuroptera[which((nd4_trop_Neuroptera$UI2=="Agriculture_High") & (nd4_trop_Neuroptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Neuroptera[which((nd4_trop_Neuroptera$UI2=="Agriculture_High") & (nd4_trop_Neuroptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd4_trop_Orthoptera$UI2=="Primary vegetation") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd4_trop_Orthoptera$UI2=="Primary vegetation") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd4_trop_Orthoptera$UI2=="Secondary vegetation") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd4_trop_Orthoptera$UI2=="Secondary vegetation") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd4_trop_Orthoptera$UI2=="Agriculture_Low") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd4_trop_Orthoptera$UI2=="Agriculture_Low") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd4_trop_Orthoptera$UI2=="Agriculture_High") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd4_trop_Orthoptera$UI2=="Agriculture_High") & (nd4_trop_Orthoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Other[which((nd4_trop_Other$UI2=="Primary vegetation") & (nd4_trop_Other$StdTmaxAnomalyRS > QPV[2])),] <- NA
Other[which((nd4_trop_Other$UI2=="Primary vegetation") & (nd4_trop_Other$StdTmaxAnomalyRS < QPV[1])),] <- NA
Other[which((nd4_trop_Other$UI2=="Secondary vegetation") & (nd4_trop_Other$StdTmaxAnomalyRS < QSV[1])),] <- NA
Other[which((nd4_trop_Other$UI2=="Secondary vegetation") & (nd4_trop_Other$StdTmaxAnomalyRS > QSV[2])),] <- NA
Other[which((nd4_trop_Other$UI2=="Agriculture_Low") & (nd4_trop_Other$StdTmaxAnomalyRS < QAL[1])),] <- NA
Other[which((nd4_trop_Other$UI2=="Agriculture_Low") & (nd4_trop_Other$StdTmaxAnomalyRS > QAL[2])),] <- NA
Other[which((nd4_trop_Other$UI2=="Agriculture_High") & (nd4_trop_Other$StdTmaxAnomalyRS < QAH[1])),] <- NA
Other[which((nd4_trop_Other$UI2=="Agriculture_High") & (nd4_trop_Other$StdTmaxAnomalyRS > QAH[2])),] <- NA

Thysanoptera[which((nd4_trop_Thysanoptera$UI2=="Primary vegetation") & (nd4_trop_Thysanoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Thysanoptera[which((nd4_trop_Thysanoptera$UI2=="Primary vegetation") & (nd4_trop_Thysanoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Thysanoptera[which((nd4_trop_Thysanoptera$UI2=="Secondary vegetation") & (nd4_trop_Thysanoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Thysanoptera[which((nd4_trop_Thysanoptera$UI2=="Secondary vegetation") & (nd4_trop_Thysanoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Thysanoptera[which((nd4_trop_Thysanoptera$UI2=="Agriculture_Low") & (nd4_trop_Thysanoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Thysanoptera[which((nd4_trop_Thysanoptera$UI2=="Agriculture_Low") & (nd4_trop_Thysanoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Thysanoptera[which((nd4_trop_Thysanoptera$UI2=="Agriculture_High") & (nd4_trop_Thysanoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Thysanoptera[which((nd4_trop_Thysanoptera$UI2=="Agriculture_High") & (nd4_trop_Thysanoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Trichoptera[which((nd4_trop_Trichoptera$UI2=="Primary vegetation") & (nd4_trop_Trichoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Trichoptera[which((nd4_trop_Trichoptera$UI2=="Primary vegetation") & (nd4_trop_Trichoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Trichoptera[which((nd4_trop_Trichoptera$UI2=="Secondary vegetation") & (nd4_trop_Trichoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Trichoptera[which((nd4_trop_Trichoptera$UI2=="Secondary vegetation") & (nd4_trop_Trichoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Trichoptera[which((nd4_trop_Trichoptera$UI2=="Agriculture_Low") & (nd4_trop_Trichoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Trichoptera[which((nd4_trop_Trichoptera$UI2=="Agriculture_Low") & (nd4_trop_Trichoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Trichoptera[which((nd4_trop_Trichoptera$UI2=="Agriculture_High") & (nd4_trop_Trichoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Trichoptera[which((nd4_trop_Trichoptera$UI2=="Agriculture_High") & (nd4_trop_Trichoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd4_trop_Blattodea$PredMedian <- ((apply(X = Blattodea,MARGIN = 1,
                                         FUN = median,na.rm=TRUE))*100)-100
nd4_trop_Blattodea$PredUpper <- ((apply(X = Blattodea,MARGIN = 1,
                                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_trop_Blattodea$PredLower <- ((apply(X = Blattodea,MARGIN = 1,
                                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

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

nd4_trop_Neuroptera$PredMedian <- ((apply(X = Neuroptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd4_trop_Neuroptera$PredUpper <- ((apply(X = Neuroptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_trop_Neuroptera$PredLower <- ((apply(X = Neuroptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_trop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                          FUN = median,na.rm=TRUE))*100)-100
nd4_trop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_trop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_trop_Other$PredMedian <- ((apply(X = Other,MARGIN = 1,
                                     FUN = median,na.rm=TRUE))*100)-100
nd4_trop_Other$PredUpper <- ((apply(X = Other,MARGIN = 1,
                                    FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_trop_Other$PredLower <- ((apply(X = Other,MARGIN = 1,
                                    FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_trop_Thysanoptera$PredMedian <- ((apply(X = Thysanoptera,MARGIN = 1,
                                            FUN = median,na.rm=TRUE))*100)-100
nd4_trop_Thysanoptera$PredUpper <- ((apply(X = Thysanoptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_trop_Thysanoptera$PredLower <- ((apply(X = Thysanoptera,MARGIN = 1,
                                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_trop_Trichoptera$PredMedian <- ((apply(X = Trichoptera,MARGIN = 1,
                                           FUN = median,na.rm=TRUE))*100)-100
nd4_trop_Trichoptera$PredUpper <- ((apply(X = Trichoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_trop_Trichoptera$PredLower <- ((apply(X = Trichoptera,MARGIN = 1,
                                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
# nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Blattodea$UI2 <- factor(nd4_trop_Blattodea$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Coleoptera$UI2 <- factor(nd4_trop_Coleoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Diptera$UI2 <- factor(nd4_trop_Diptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Hemiptera$UI2 <- factor(nd4_trop_Hemiptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Hymenoptera$UI2 <- factor(nd4_trop_Hymenoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Lepidoptera$UI2 <- factor(nd4_trop_Lepidoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Neuroptera$UI2 <- factor(nd4_trop_Neuroptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Orthoptera$UI2 <- factor(nd4_trop_Orthoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Other$UI2 <- factor(nd4_trop_Other$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Thysanoptera$UI2 <- factor(nd4_trop_Thysanoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_trop_Trichoptera$UI2 <- factor(nd4_trop_Trichoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# add zone factor
nd4_trop_Blattodea$Zone <- as.factor("Tropical")
nd4_trop_Coleoptera$Zone <- as.factor("Tropical")
nd4_trop_Diptera$Zone <- as.factor("Tropical")
nd4_trop_Hemiptera$Zone <- as.factor("Tropical")
nd4_trop_Hymenoptera$Zone <- as.factor("Tropical")
nd4_trop_Lepidoptera$Zone <- as.factor("Tropical")
nd4_trop_Neuroptera$Zone <- as.factor("Tropical")
nd4_trop_Orthoptera$Zone <- as.factor("Tropical")
nd4_trop_Other$Zone <- as.factor("Tropical")
nd4_trop_Thysanoptera$Zone <- as.factor("Tropical")
nd4_trop_Trichoptera$Zone <- as.factor("Tropical")

## NonTropical ##

nd4_nontrop <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelRich_nontrop$data$StdTmaxAnomalyRS),
                       to = max(MaxAnomalyModelRich_nontrop$data$StdTmaxAnomalyRS),
                       length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelRich_nontrop$data$UI2)),
  Order=factor(c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")))

# back transform the predictors
nd4_nontrop$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4_nontrop$StdTmaxAnomalyRS,
  originalX = nontrop$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd4_nontrop$LogRich <- 0
nd4_nontrop$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
# does this for each Order
# reference row is 4th row, every 400 rows (see 'Values')
refRow <- which((nd4_nontrop$UI2=="Primary vegetation") & (nd4_nontrop$StdTmaxAnomaly==min(abs(nd4_nontrop$StdTmaxAnomaly))))
# 15th row, every 400 rows

# adjust plot 1: max anomaly and abundance

QPV <- quantile(x = MaxAnomalyModelRich_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_nontrop$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAnomalyModelRich_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_nontrop$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAnomalyModelRich_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_nontrop$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAnomalyModelRich_nontrop$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich_nontrop$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
sr.preds.tmax.nontrop <- PredictGLMERRandIter(model = MaxAnomalyModelRich_nontrop$model,data = nd4_nontrop)

# back transform the abundance values
sr.preds.tmax.nontrop <- exp(sr.preds.tmax.nontrop)-0.01


# another try!
number_of_chunks = 11
list_sr.preds.tmax.nontrop <- lapply(seq(1, NROW(sr.preds.tmax.nontrop), ceiling(NROW(sr.preds.tmax.nontrop)/number_of_chunks)),
                                     function(i) sr.preds.tmax.nontrop[i:min(i + ceiling(NROW(sr.preds.tmax.nontrop)/number_of_chunks) - 1, NROW(sr.preds.tmax.nontrop)),])
# success!
# creates list of matrices
# name them
names(list_sr.preds.tmax.nontrop) <- c("Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Neuroptera","Orthoptera","Other","Thysanoptera","Trichoptera")
list2env(list_sr.preds.tmax.nontrop,globalenv())

# tim's suggestion
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

Blattodea[which((nd4_nontrop_Blattodea$UI2=="Primary vegetation") & (nd4_nontrop_Blattodea$StdTmaxAnomalyRS > QPV[2])),] <- NA
Blattodea[which((nd4_nontrop_Blattodea$UI2=="Primary vegetation") & (nd4_nontrop_Blattodea$StdTmaxAnomalyRS < QPV[1])),] <- NA
Blattodea[which((nd4_nontrop_Blattodea$UI2=="Secondary vegetation") & (nd4_nontrop_Blattodea$StdTmaxAnomalyRS < QSV[1])),] <- NA
Blattodea[which((nd4_nontrop_Blattodea$UI2=="Secondary vegetation") & (nd4_nontrop_Blattodea$StdTmaxAnomalyRS > QSV[2])),] <- NA
Blattodea[which((nd4_nontrop_Blattodea$UI2=="Agriculture_Low") & (nd4_nontrop_Blattodea$StdTmaxAnomalyRS < QAL[1])),] <- NA
Blattodea[which((nd4_nontrop_Blattodea$UI2=="Agriculture_Low") & (nd4_nontrop_Blattodea$StdTmaxAnomalyRS > QAL[2])),] <- NA
Blattodea[which((nd4_nontrop_Blattodea$UI2=="Agriculture_High") & (nd4_nontrop_Blattodea$StdTmaxAnomalyRS < QAH[1])),] <- NA
Blattodea[which((nd4_nontrop_Blattodea$UI2=="Agriculture_High") & (nd4_nontrop_Blattodea$StdTmaxAnomalyRS > QAH[2])),] <- NA

Coleoptera[which((nd4_nontrop_Coleoptera$UI2=="Primary vegetation") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Coleoptera[which((nd4_nontrop_Coleoptera$UI2=="Primary vegetation") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Coleoptera[which((nd4_nontrop_Coleoptera$UI2=="Secondary vegetation") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Coleoptera[which((nd4_nontrop_Coleoptera$UI2=="Secondary vegetation") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Coleoptera[which((nd4_nontrop_Coleoptera$UI2=="Agriculture_Low") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Coleoptera[which((nd4_nontrop_Coleoptera$UI2=="Agriculture_Low") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Coleoptera[which((nd4_nontrop_Coleoptera$UI2=="Agriculture_High") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Coleoptera[which((nd4_nontrop_Coleoptera$UI2=="Agriculture_High") & (nd4_nontrop_Coleoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Diptera[which((nd4_nontrop_Diptera$UI2=="Primary vegetation") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Diptera[which((nd4_nontrop_Diptera$UI2=="Primary vegetation") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Diptera[which((nd4_nontrop_Diptera$UI2=="Secondary vegetation") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Diptera[which((nd4_nontrop_Diptera$UI2=="Secondary vegetation") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Diptera[which((nd4_nontrop_Diptera$UI2=="Agriculture_Low") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Diptera[which((nd4_nontrop_Diptera$UI2=="Agriculture_Low") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Diptera[which((nd4_nontrop_Diptera$UI2=="Agriculture_High") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Diptera[which((nd4_nontrop_Diptera$UI2=="Agriculture_High") & (nd4_nontrop_Diptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hemiptera[which((nd4_nontrop_Hemiptera$UI2=="Primary vegetation") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hemiptera[which((nd4_nontrop_Hemiptera$UI2=="Primary vegetation") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hemiptera[which((nd4_nontrop_Hemiptera$UI2=="Secondary vegetation") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hemiptera[which((nd4_nontrop_Hemiptera$UI2=="Secondary vegetation") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hemiptera[which((nd4_nontrop_Hemiptera$UI2=="Agriculture_Low") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hemiptera[which((nd4_nontrop_Hemiptera$UI2=="Agriculture_Low") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hemiptera[which((nd4_nontrop_Hemiptera$UI2=="Agriculture_High") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hemiptera[which((nd4_nontrop_Hemiptera$UI2=="Agriculture_High") & (nd4_nontrop_Hemiptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Hymenoptera[which((nd4_nontrop_Hymenoptera$UI2=="Primary vegetation") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Hymenoptera[which((nd4_nontrop_Hymenoptera$UI2=="Primary vegetation") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Hymenoptera[which((nd4_nontrop_Hymenoptera$UI2=="Secondary vegetation") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Hymenoptera[which((nd4_nontrop_Hymenoptera$UI2=="Secondary vegetation") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Hymenoptera[which((nd4_nontrop_Hymenoptera$UI2=="Agriculture_Low") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Hymenoptera[which((nd4_nontrop_Hymenoptera$UI2=="Agriculture_Low") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Hymenoptera[which((nd4_nontrop_Hymenoptera$UI2=="Agriculture_High") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Hymenoptera[which((nd4_nontrop_Hymenoptera$UI2=="Agriculture_High") & (nd4_nontrop_Hymenoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Lepidoptera[which((nd4_nontrop_Lepidoptera$UI2=="Primary vegetation") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Lepidoptera[which((nd4_nontrop_Lepidoptera$UI2=="Primary vegetation") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Lepidoptera[which((nd4_nontrop_Lepidoptera$UI2=="Secondary vegetation") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Lepidoptera[which((nd4_nontrop_Lepidoptera$UI2=="Secondary vegetation") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Lepidoptera[which((nd4_nontrop_Lepidoptera$UI2=="Agriculture_Low") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Lepidoptera[which((nd4_nontrop_Lepidoptera$UI2=="Agriculture_Low") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Lepidoptera[which((nd4_nontrop_Lepidoptera$UI2=="Agriculture_High") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Lepidoptera[which((nd4_nontrop_Lepidoptera$UI2=="Agriculture_High") & (nd4_nontrop_Lepidoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Neuroptera[which((nd4_nontrop_Neuroptera$UI2=="Primary vegetation") & (nd4_nontrop_Neuroptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Neuroptera[which((nd4_nontrop_Neuroptera$UI2=="Primary vegetation") & (nd4_nontrop_Neuroptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Neuroptera[which((nd4_nontrop_Neuroptera$UI2=="Secondary vegetation") & (nd4_nontrop_Neuroptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Neuroptera[which((nd4_nontrop_Neuroptera$UI2=="Secondary vegetation") & (nd4_nontrop_Neuroptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Neuroptera[which((nd4_nontrop_Neuroptera$UI2=="Agriculture_Low") & (nd4_nontrop_Neuroptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Neuroptera[which((nd4_nontrop_Neuroptera$UI2=="Agriculture_Low") & (nd4_nontrop_Neuroptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Neuroptera[which((nd4_nontrop_Neuroptera$UI2=="Agriculture_High") & (nd4_nontrop_Neuroptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Neuroptera[which((nd4_nontrop_Neuroptera$UI2=="Agriculture_High") & (nd4_nontrop_Neuroptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Orthoptera[which((nd4_nontrop_Orthoptera$UI2=="Primary vegetation") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Orthoptera[which((nd4_nontrop_Orthoptera$UI2=="Primary vegetation") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Orthoptera[which((nd4_nontrop_Orthoptera$UI2=="Secondary vegetation") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Orthoptera[which((nd4_nontrop_Orthoptera$UI2=="Secondary vegetation") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Orthoptera[which((nd4_nontrop_Orthoptera$UI2=="Agriculture_Low") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Orthoptera[which((nd4_nontrop_Orthoptera$UI2=="Agriculture_Low") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Orthoptera[which((nd4_nontrop_Orthoptera$UI2=="Agriculture_High") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Orthoptera[which((nd4_nontrop_Orthoptera$UI2=="Agriculture_High") & (nd4_nontrop_Orthoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Other[which((nd4_nontrop_Other$UI2=="Primary vegetation") & (nd4_nontrop_Other$StdTmaxAnomalyRS > QPV[2])),] <- NA
Other[which((nd4_nontrop_Other$UI2=="Primary vegetation") & (nd4_nontrop_Other$StdTmaxAnomalyRS < QPV[1])),] <- NA
Other[which((nd4_nontrop_Other$UI2=="Secondary vegetation") & (nd4_nontrop_Other$StdTmaxAnomalyRS < QSV[1])),] <- NA
Other[which((nd4_nontrop_Other$UI2=="Secondary vegetation") & (nd4_nontrop_Other$StdTmaxAnomalyRS > QSV[2])),] <- NA
Other[which((nd4_nontrop_Other$UI2=="Agriculture_Low") & (nd4_nontrop_Other$StdTmaxAnomalyRS < QAL[1])),] <- NA
Other[which((nd4_nontrop_Other$UI2=="Agriculture_Low") & (nd4_nontrop_Other$StdTmaxAnomalyRS > QAL[2])),] <- NA
Other[which((nd4_nontrop_Other$UI2=="Agriculture_High") & (nd4_nontrop_Other$StdTmaxAnomalyRS < QAH[1])),] <- NA
Other[which((nd4_nontrop_Other$UI2=="Agriculture_High") & (nd4_nontrop_Other$StdTmaxAnomalyRS > QAH[2])),] <- NA

Thysanoptera[which((nd4_nontrop_Thysanoptera$UI2=="Primary vegetation") & (nd4_nontrop_Thysanoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Thysanoptera[which((nd4_nontrop_Thysanoptera$UI2=="Primary vegetation") & (nd4_nontrop_Thysanoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Thysanoptera[which((nd4_nontrop_Thysanoptera$UI2=="Secondary vegetation") & (nd4_nontrop_Thysanoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Thysanoptera[which((nd4_nontrop_Thysanoptera$UI2=="Secondary vegetation") & (nd4_nontrop_Thysanoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Thysanoptera[which((nd4_nontrop_Thysanoptera$UI2=="Agriculture_Low") & (nd4_nontrop_Thysanoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Thysanoptera[which((nd4_nontrop_Thysanoptera$UI2=="Agriculture_Low") & (nd4_nontrop_Thysanoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Thysanoptera[which((nd4_nontrop_Thysanoptera$UI2=="Agriculture_High") & (nd4_nontrop_Thysanoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Thysanoptera[which((nd4_nontrop_Thysanoptera$UI2=="Agriculture_High") & (nd4_nontrop_Thysanoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA

Trichoptera[which((nd4_nontrop_Trichoptera$UI2=="Primary vegetation") & (nd4_nontrop_Trichoptera$StdTmaxAnomalyRS > QPV[2])),] <- NA
Trichoptera[which((nd4_nontrop_Trichoptera$UI2=="Primary vegetation") & (nd4_nontrop_Trichoptera$StdTmaxAnomalyRS < QPV[1])),] <- NA
Trichoptera[which((nd4_nontrop_Trichoptera$UI2=="Secondary vegetation") & (nd4_nontrop_Trichoptera$StdTmaxAnomalyRS < QSV[1])),] <- NA
Trichoptera[which((nd4_nontrop_Trichoptera$UI2=="Secondary vegetation") & (nd4_nontrop_Trichoptera$StdTmaxAnomalyRS > QSV[2])),] <- NA
Trichoptera[which((nd4_nontrop_Trichoptera$UI2=="Agriculture_Low") & (nd4_nontrop_Trichoptera$StdTmaxAnomalyRS < QAL[1])),] <- NA
Trichoptera[which((nd4_nontrop_Trichoptera$UI2=="Agriculture_Low") & (nd4_nontrop_Trichoptera$StdTmaxAnomalyRS > QAL[2])),] <- NA
Trichoptera[which((nd4_nontrop_Trichoptera$UI2=="Agriculture_High") & (nd4_nontrop_Trichoptera$StdTmaxAnomalyRS < QAH[1])),] <- NA
Trichoptera[which((nd4_nontrop_Trichoptera$UI2=="Agriculture_High") & (nd4_nontrop_Trichoptera$StdTmaxAnomalyRS > QAH[2])),] <- NA


# Get the median, upper and lower quants for the plot

nd4_nontrop_Blattodea$PredMedian <- ((apply(X = Blattodea,MARGIN = 1,
                                            FUN = median,na.rm=TRUE))*100)-100
nd4_nontrop_Blattodea$PredUpper <- ((apply(X = Blattodea,MARGIN = 1,
                                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_nontrop_Blattodea$PredLower <- ((apply(X = Blattodea,MARGIN = 1,
                                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

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

nd4_nontrop_Neuroptera$PredMedian <- ((apply(X = Neuroptera,MARGIN = 1,
                                             FUN = median,na.rm=TRUE))*100)-100
nd4_nontrop_Neuroptera$PredUpper <- ((apply(X = Neuroptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_nontrop_Neuroptera$PredLower <- ((apply(X = Neuroptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_nontrop_Orthoptera$PredMedian <- ((apply(X = Orthoptera,MARGIN = 1,
                                             FUN = median,na.rm=TRUE))*100)-100
nd4_nontrop_Orthoptera$PredUpper <- ((apply(X = Orthoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_nontrop_Orthoptera$PredLower <- ((apply(X = Orthoptera,MARGIN = 1,
                                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_nontrop_Other$PredMedian <- ((apply(X = Other,MARGIN = 1,
                                        FUN = median,na.rm=TRUE))*100)-100
nd4_nontrop_Other$PredUpper <- ((apply(X = Other,MARGIN = 1,
                                       FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_nontrop_Other$PredLower <- ((apply(X = Other,MARGIN = 1,
                                       FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_nontrop_Thysanoptera$PredMedian <- ((apply(X = Thysanoptera,MARGIN = 1,
                                               FUN = median,na.rm=TRUE))*100)-100
nd4_nontrop_Thysanoptera$PredUpper <- ((apply(X = Thysanoptera,MARGIN = 1,
                                              FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_nontrop_Thysanoptera$PredLower <- ((apply(X = Thysanoptera,MARGIN = 1,
                                              FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4_nontrop_Trichoptera$PredMedian <- ((apply(X = Trichoptera,MARGIN = 1,
                                              FUN = median,na.rm=TRUE))*100)-100
nd4_nontrop_Trichoptera$PredUpper <- ((apply(X = Trichoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4_nontrop_Trichoptera$PredLower <- ((apply(X = Trichoptera,MARGIN = 1,
                                             FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


# set factor levels
# nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Blattodea$UI2 <- factor(nd4_nontrop_Blattodea$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Coleoptera$UI2 <- factor(nd4_nontrop_Coleoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Diptera$UI2 <- factor(nd4_nontrop_Diptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Hemiptera$UI2 <- factor(nd4_nontrop_Hemiptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Hymenoptera$UI2 <- factor(nd4_nontrop_Hymenoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Lepidoptera$UI2 <- factor(nd4_nontrop_Lepidoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Neuroptera$UI2 <- factor(nd4_nontrop_Neuroptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Orthoptera$UI2 <- factor(nd4_nontrop_Orthoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Other$UI2 <- factor(nd4_nontrop_Other$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Thysanoptera$UI2 <- factor(nd4_nontrop_Thysanoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
nd4_nontrop_Trichoptera$UI2 <- factor(nd4_nontrop_Trichoptera$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

# add zone factor
nd4_nontrop_Blattodea$Zone <- as.factor("NonTropical")
nd4_nontrop_Coleoptera$Zone <- as.factor("NonTropical")
nd4_nontrop_Diptera$Zone <- as.factor("NonTropical")
nd4_nontrop_Hemiptera$Zone <- as.factor("NonTropical")
nd4_nontrop_Hymenoptera$Zone <- as.factor("NonTropical")
nd4_nontrop_Lepidoptera$Zone <- as.factor("NonTropical")
nd4_nontrop_Neuroptera$Zone <- as.factor("NonTropical")
nd4_nontrop_Orthoptera$Zone <- as.factor("NonTropical")
nd4_nontrop_Other$Zone <- as.factor("NonTropical")
nd4_nontrop_Thysanoptera$Zone <- as.factor("NonTropical")
nd4_nontrop_Trichoptera$Zone <- as.factor("NonTropical")

# put use rbind to add nd4_nontrop to nd4_trop to make one data table for plotting

nd4_Blattodea <- rbind(nd4_trop_Blattodea,nd4_nontrop_Blattodea)
nd4_Coleoptera <- rbind(nd4_trop_Coleoptera,nd4_nontrop_Coleoptera)
nd4_Diptera <- rbind(nd4_trop_Diptera,nd4_nontrop_Diptera)
nd4_Hemiptera <- rbind(nd4_trop_Hemiptera,nd4_nontrop_Hemiptera)
nd4_Hymenoptera <- rbind(nd4_trop_Hymenoptera,nd4_nontrop_Hymenoptera)
nd4_Lepidoptera <- rbind(nd4_trop_Lepidoptera,nd4_nontrop_Lepidoptera)
nd4_Neuroptera <- rbind(nd4_trop_Neuroptera,nd4_nontrop_Neuroptera)
nd4_Orthoptera <- rbind(nd4_trop_Orthoptera,nd4_nontrop_Orthoptera)
nd4_Other <- rbind(nd4_trop_Other,nd4_nontrop_Other)
nd4_Thysanoptera <- rbind(nd4_trop_Thysanoptera,nd4_nontrop_Thysanoptera)
nd4_Trichoptera <- rbind(nd4_trop_Trichoptera,nd4_nontrop_Trichoptera)

# plot

p_blattodea <- ggplot(data = nd4_Blattodea, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Blattodea$PredLower, ymax = nd4_Blattodea$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) + 
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) + # extended y axis for orders with high larger confidence intervals
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        # legend.background = element_blank(), 
        # legend.text = element_text(size = 6), 
        # legend.title = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2)) + 
  ggtitle("Blattodea")

p_coleoptera <- ggplot(data = nd4_Coleoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Coleoptera$PredLower, ymax = nd4_Coleoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Diptera$PredLower, ymax = nd4_Diptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Hemiptera$PredLower, ymax = nd4_Hemiptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Hymenoptera$PredLower, ymax = nd4_Hymenoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Lepidoptera$PredLower, ymax = nd4_Lepidoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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

p_neuroptera <- ggplot(data = nd4_Neuroptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Neuroptera$PredLower, ymax = nd4_Neuroptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  ggtitle("Neuroptera")

p_other <- ggplot(data = nd4_Other, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Other$PredLower, ymax = nd4_Other$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  ggtitle("Other")

p_orthoptera <- ggplot(data = nd4_Orthoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Orthoptera$PredLower, ymax = nd4_Orthoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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

p_thysanoptera <- ggplot(data = nd4_Thysanoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Thysanoptera$PredLower, ymax = nd4_Thysanoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  ggtitle("Thysanoptera")

p_trichoptera <- ggplot(data = nd4_Trichoptera, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(linetype = Zone, col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd4_Trichoptera$PredLower, ymax = nd4_Trichoptera$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0,0.5, 1, 1.5, 2), limits = c(0, 2)) +
  #scale_y_continuous(breaks = c(-100,-75, -50, -25, 0, 25, 50, 75, 100, 125,150,175), limits = c(-100, 175)) +
  scale_y_continuous(breaks = c(-100, -50,  0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500), limits = c(-100, 500)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Maximum n/Temperature Anomaly") +
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
  ggtitle("Trichoptera")

# get the legend
legend <- get_legend(
  p_blattodea +
    guides(color = guide_legend(nrow = 1),
           linetype = guide_legend (nrow=1)) +
    theme(legend.position = "right",
          legend.background = element_blank(), 
          legend.text = element_text(size = 6), 
          legend.title = element_blank())
)


# put them all together to save them
MaxAnomRich <- cowplot::plot_grid(p_blattodea, p_coleoptera,p_diptera,p_hemiptera,p_hymenoptera,p_lepidoptera,p_neuroptera,p_orthoptera,p_other,p_thysanoptera,p_trichoptera,legend)

# save them
ggsave(filename = paste0(outDir, "MaxAnomRich.pdf"), plot = MaxAnomRich, width = 200, height = 150, units = "mm", dpi = 300)
#ggsave(filename = paste0(outDir, "MaxAnomRich_extended yaxis.pdf"), plot = MaxAnomRich, width = 200, height = 150, units = "mm", dpi = 300)
