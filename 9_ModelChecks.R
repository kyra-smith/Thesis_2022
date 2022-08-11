##%######################################################%##
#                                                          #
####                    Model Checks                    ####
#                                                          #
##%######################################################%##

rm(list = ls())

# directories
LUIDir <- "C:/Users/Kyra/Documents/GLITRS/Code/2_RunSimpleLUIModel/Output/"
LUICCDir <- "C:/Users/Kyra/Documents/GLITRS/Code/5_RunLUIClimateModels/"
predictsDataDir <- "C:/Users/Kyra/Documents/GLITRS/Code/5_RunLUIClimateModels/"
outdir <- "C:/Users/Kyra/Documents/GLITRS/Code/9_ModelChecks/"
if(!dir.exists(outdir)) dir.create(outdir)


#sink(paste0(outdir,"log.txt"))

#t.start <- Sys.time()

#print(t.start)

# load libraries
library(cowplot)
library(ggplot2)
library(StatisticalModels)
library(gridGraphics)

# read in model files
LUIAbund <- readRDS(file = paste0(LUIDir,"am3.3.rds"))
LUIRich <- readRDS(file = paste0(LUIDir,"sm3.3.rds"))

load(paste0(LUICCDir, "MeanAnomalyModelAbund.rdata"))
load(paste0(LUICCDir, "MeanAnomalyModelRich.rdata"))
load(paste0(LUICCDir, "MaxAnomalyModelAbund.rdata"))
load(paste0(LUICCDir, "MaxAnomalyModelRich.rdata"))


# try to write a function to run these checks on each model output and save a pdf of all figs.

mod_list <- c("LUIAbund", "LUIRich","MeanAnomalyModelAbund", "MeanAnomalyModelRich", "MaxAnomalyModelAbund", "MaxAnomalyModelRich")
mod_list1 <- c("LUIAbund", "LUIRich")
mod_list2 <- c("MeanAnomalyModelAbund", "MeanAnomalyModelRich", "MaxAnomalyModelAbund", "MaxAnomalyModelRich")

# load dataset
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSitesClimate_Data.rds"))
predictsSites <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ]

#x <- mod_list[1]

# loop through models and generate plots for checking assumptions etc
for(x in mod_list){
  
  # only do this one if not SR model
  # grepl() searches for matches to "Abun" in "x", the list of models
  if(grepl("Abun", x) ==1) {
    
    
    ## 1. Checking the fitted vs residuals relationship
    p1 <- plot(get(x)$model)
  }
  
  ## 2. Normality of Residuals
  pdf(NULL)
  dev.control(displaylist="enable")
  qqnorm(resid(get(x)$model), main = "")
  qqline(resid(get(x)$model))
  p2 <- recordPlot()
  invisible(dev.off())
  
  
  if(grepl("Abun", x) == 1) {
    
    
    predData <- predictsSites[!is.na(predictsSites$LogAbund), ]
    
    
    # 3. plot of observed vs fitted values
    pdf(NULL)
    dev.control(displaylist="enable")
    plot(predData$LogAbund,fitted(get(x)), 
         xlab = "Observed values", ylab = "Fitted values") 
    abline(a = 0, b = 1, col = "red", lwd = 2)
    p3 <- recordPlot()
    invisible(dev.off())
    
    
    cowplot::plot_grid(p1,p2,p3,
                       labels = c("A.", "B.", "C."))
    
    ggsave(file = paste0(outdir, x, "_model_checks.pdf"), height = 10, width = 10)
    
    rm(p1, p2, p3)
    
    
  }else{
    
    
    # 4. plot of observed vs fitted values
    
    pdf(NULL)
    dev.control(displaylist="enable")
    plot(predictsSites$Species_richness,fitted(get(x)), 
         xlab = "Observed values", ylab = "Fitted values") 
    abline(a = 0, b = 1, col = "red", lwd = 2)
    p3 <- recordPlot()
    invisible(dev.off())
    
    cowplot::plot_grid(p2,p3,
                       labels = c("A.", "B."))#
    
    ggsave(file = paste0(outdir, x, "_model_checks.pdf")) 
    
    rm(p2, p3)
    
  }
  
}


#t.end <- Sys.time()

#print(round(t.end - t.start,0))

#sink()