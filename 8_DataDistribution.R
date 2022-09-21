# Get Set Up

inDir <- "C:/Users/Kyra/Documents/GLITRS/Code/4_PREDICTSMatchClimateIndex/"
outDir <- "C:/Users/Kyra/Documents/GLITRS/Code/8_DataDistribution/"
if(!dir.exists(outDir)) dir.create(outDir)


## check spread of climate data within groups ##
library(ggplot2)
library(dplyr)
library(cowplot)

# load PREDICTSSites_ClimateData.rds
sites <- readRDS(paste0(inDir,"PREDICTSSites_Climate.rds"))

# change PREDICTSSites_Climate.rds (LargeSpatialPointsDataFrame) to sites, a data frame
sites <- as.data.frame(sites)

# check for missing values
na <- sites[rowSums(is.na(sites)) > 0, ] 
# 3239 obs. of 36 variables

# which columns contain these values
names(which(colSums(is.na(na)) > 0))
# [11] "Total_abundance"   "ChaoR"             "Richness_rarefied" "LogAbund"          "avg_temp"          "TmeanAnomaly"      "StdTmeanAnomaly"  
# [8] "max_temp"          "TmaxAnomaly"       "StdTmaxAnomaly"  

# remove diversity metric columns and just look at the climate data, order, longitude, and latitude columns
na <- sites %>% select(Order,coords.x1, coords.x2, avg_temp,TmeanAnomaly, StdTmeanAnomaly, max_temp, TmaxAnomaly, StdTmaxAnomaly)
# 11865 obs. of 9 variables
na <- na[rowSums(is.na(na)) > 0, ] 
# 26 obs. of 9 variables
# now let's investigate
# diptera has 1 missing value, at long,lat -48.475833, -27.59083 (This is a mistake - the coordinates are in the Atlantic Ocean)
# hymenoptera has 25 missing values
  # 7 at roughly the same coordinates (long,lat, ~7.6, ~ 46.3) (This is in Switzerland)
  # 18 at roughly the same coordinates (long,lat, ~ -7.5, ~ 57.6) (In the North Atlantic, off the coast of North Scotland)

# remove the rows with missing values for avg_temp (which means the other parameters are also NA)
sites <- sites %>% filter(!is.na(avg_temp))
# 11839 obs. of 36 variables

# number of active months
p_n_months<-ggplot(sites)+geom_boxplot(aes(x=Order,y=n_months,fill=Order)) +
  labs(title = "Number of Active Months",
       x = "Order",
       y = "Number of Active Months / Year")+
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5))

# save it 
# ggsave(filename = paste0(outDir, "Active Months.pdf"), height = 4, width = 8)
ggsave("ActiveMonths.jpeg", device ="jpeg", path = outDir, width=20, height=10, units="cm", dpi = 350)

# latitude
p_latitude<-ggplot(sites)+geom_boxplot(aes(x=Order,y=coords.x2,fill=Order)) +
  labs(title = "Latitudinal Distribution",
       x = "Order",
       y = "Latitude")+
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5))

# save it 
# ggsave(filename = paste0(outDir, "Latitude.pdf"), height = 4, width = 8)
ggsave("Latitude.jpeg", device ="jpeg", path = outDir, width=20, height=10, units="cm", dpi = 350)

# avg_temp
p_mean_temp<-ggplot(sites)+geom_boxplot(aes(x=Order,y=avg_temp,fill=Order)) +
  labs(y = "Mean Site Temperature (°C)")+
  theme_gray() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.y=element_text(size=6),
        legend.position="none")

# TmeanAnomaly
p_TmeanAnomaly<-ggplot(sites)+geom_boxplot(aes(x=Order,y=TmeanAnomaly,fill=Order)) +
  labs(y = "Absolute Temperature Anomaly (°C)")+
  theme_gray() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.y=element_text(size=6),
        legend.position="none")

# StdTmeanAnomaly
p_StdTmeanAnomaly<-ggplot(sites)+geom_boxplot(aes(x=Order,y=StdTmeanAnomaly,fill=Order)) +
  labs(x = "Order",
       y = "Standardized Temperature Anomaly")+
  theme_gray() +
  theme(axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=6),
        axis.title.y=element_text(size=8),
        axis.text.y=element_text(size=6),
        legend.position="none")

# title for mean plots
title_mean <- ggdraw() + 
  draw_label(
    "Mean Site Temperatures and Anomalies",
    fontface = 'bold',
    size = 10)

# put together the plots of mean temp
mean_plots <- cowplot::plot_grid(title_mean,p_mean_temp,p_TmeanAnomaly,p_StdTmeanAnomaly,
                                ncol=1,
                                rel_heights = c(0.1, 1,1,1))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

# save it
ggsave("MeanTemperatureDistribution.jpeg", device ="jpeg", path = outDir, width=20, height=20, units="cm", dpi = 350)

# max_temp
p_max_temp<-ggplot(sites)+geom_boxplot(aes(x=Order,y=max_temp,fill=Order)) +
  labs(y = "Maximum Site Temperature (°C)")+
  theme_gray()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=6),
        legend.position="none")

# TmaxAnomaly
p_TmaxAnomaly<-ggplot(sites)+geom_boxplot(aes(x=Order,y=TmaxAnomaly,fill=Order)) +
  labs(y = "Absolute Maximum Temperature Anomaly (°C)")+
  theme_gray() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=6),
        legend.position="none")

# StdTmaxAnomaly
p_StdTmaxAnomaly<-ggplot(sites)+geom_boxplot(aes(x=Order,y=StdTmaxAnomaly,fill=Order)) +
  labs(x = "Order",
       y = "Standardized Maximum Temperature Anomaly")+
  theme_gray() +
  theme(axis.title.x=element_text(size=10),
        axis.text.x=element_text(size=6),
        axis.title.y=element_text(size=6),
        axis.text.y=element_text(size=6),
        legend.position="none")

# title for max plots
title_max <- ggdraw() + 
  draw_label(
    "Maximum Site Temperatures and Anomalies",
    fontface = 'bold',
    size = 10)

# put together the plots of max temp
max_plots <- cowplot::plot_grid(title_max,p_max_temp,p_TmaxAnomaly,p_StdTmaxAnomaly,
                        ncol=1,
                        rel_heights = c(0.1,1,1,1))+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))

# save it 
ggsave("MaxTemperatureDistribution.jpeg", device ="jpeg", path = outDir, width=20, height=20, units="cm", dpi = 350)

# summary stats
# mean temperature
summary_mean <- sites %>%
  group_by(Order) %>%
  summarise(Observations = length(Order),
            mean_Mean_Temp = mean(avg_temp),
            mean_Temp_Anomaly = mean(TmeanAnomaly),
            mean_Std_Temp_Anomaly = mean(StdTmeanAnomaly))
summary_mean

# max temperature
summary_max <- sites %>%
  group_by(Order) %>%
  summarise(Observations = length(Order),
            mean_Max_Temp = mean(max_temp),
            mean_Temp_Anomaly = mean(TmaxAnomaly),
            mean_Std_Temp_Anomaly = mean(StdTmaxAnomaly))
summary_max

# percentiles
summary_meanQ <- sites %>%
  group_by(Order) %>%
  summarise(quantile = scales::percent(c(0.05, 0.25, 0.5, 0.75, 0.95)),
            Qavg_temp = quantile(avg_temp, c(0.05, 0.25, 0.5, 0.75, 0.95)),
            QTmeanAnomaly = quantile(TmeanAnomaly, c(0.05, 0.25, 0.5, 0.75, 0.95)),
            QStdTmeanAnomaly = quantile(StdTmeanAnomaly, c(0.05, 0.25, 0.5, 0.75, 0.95)))

# percentiles
summary_maxQ <- sites %>%
  group_by(Order) %>%
  summarise(quantile = scales::percent(c(0.05, 0.25, 0.5, 0.75, 0.95)),
            Qmax_temp = quantile(max_temp, c(0.05, 0.25, 0.5, 0.75, 0.95)),
            QTmaxAnomaly = quantile(TmaxAnomaly, c(0.05, 0.25, 0.5, 0.75, 0.95)),
            QStdTmaxAnomaly = quantile(StdTmaxAnomaly, c(0.05, 0.25, 0.5, 0.75, 0.95)))

# I want to look at the distribution in climate anomalies between Tropical and Non-Tropical sites

# sites$Realm <- ifelse(sites$Latitude >= -23.5 & sites$Latitude <= 23.5, "Tropical", "NonTropical")
# sites$Realm <- factor(sites$Realm, levels = c("NonTropical", "Tropical"))
# sites <- sites %>%
#   filter(!is.na(Realm))

# now for plots!
# number of active months
p_n_months_realm<-ggplot(sites)+geom_boxplot(aes(x=Realm,y=n_months,fill=Realm)) +
  labs(title = "Number of Active Months",
       x = "Realm",
       y = "Number of Active Months / Year")+
  scale_fill_manual(values=c("#C1FFC1","#CAE1FF"))+
  theme_gray() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=12),
        legend.position="none")

# save it 
ggsave(filename = paste0(outDir, "Active Months by Realm.pdf"), height = 6, width = 8)

# just looking at Standardized Climate Anomaly
# StdTmeanAnomaly
p_StdTmeanAnomaly_realm<-ggplot(sites)+geom_boxplot(aes(x=Realm,y=StdTmeanAnomaly,fill=Realm)) +
  labs(x = "Realm",
       y = "Standardized Climate Anomaly (Mean Temp)")+
  scale_fill_manual(values=c("#C1FFC1","#CAE1FF"))+
  theme_gray() +
  theme(axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=12),
        legend.position="none")

ggsave(filename = paste0(outDir, "StdMeanAnom_Tropical.pdf"), height = 6, width = 8)


# StdTmaxAnomaly
p_StdTmaxAnomaly_realm<-ggplot(sites)+geom_boxplot(aes(x=Realm,y=StdTmaxAnomaly,fill=Realm)) +
  labs(x = "Realm",
       y = "Standardized Climate Anomaly (Max Temp)")+
  scale_fill_manual(values=c("#C1FFC1","#CAE1FF"))+
  theme_gray() +
  theme(axis.title.x=element_text(size=14),
        axis.text.x=element_text(size=12),
        axis.title.y=element_text(size=14),
        axis.text.y=element_text(size=12),
        legend.position="none")

ggsave(filename = paste0(outDir, "StdMaxAnom_Tropical.pdf"), height = 6, width = 8)

# selection_table <- data.frame("Order" = c("Coleoptera","Diptera","Hemiptera","Hymenoptera","Lepidoptera","Orthoptera")),
#                               "Count" = c("Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Species_richness ~ LUI + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Species_richness ~ Order + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Species_richness ~ Order + LUI + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Species_richness ~ Order * LUI + (1|SS) + (1|SSB) + (1|SSBS)",
#                                           "Total_abundance ~ 1 + (1|SS) + (1|SSB)",
#                                           "Total_abundance ~ LUI + (1|SS) + (1|SSB)",
#                                           "Total_abundance ~ Order + (1|SS) + (1|SSB)",
#                                           "Total_abundance ~ Order + LUI + (1|SS) + (1|SSB)",
#                                           "Total_abundance ~ Order * LUI + (1|SS) + (1|SSB)"),
#                               "AIC" = c(AIC(sm0$model), AIC(sm3$model), AIC(sm0.2$model), AIC(sm3.2$model), AIC(sm3.3$model),  
#                                         AIC(am0$model), AIC(am3$model), AIC(am0.2$model), AIC(am3.2$model),AIC(am3.3$model))) %>%
#   group_by(Response) %>%                              
#   mutate(deltaAIC = cumsum(c(0, diff(AIC)))) %>%
#   ungroup() %>%
#   gt()
# 
# gtsave(selection_table,"C:/Users/Kyra/Documents/GLITRS/Code/2_RunSimpleLUIModel/Output/LUIModels_Selection.png")

# Summary Table
summary_LUI <- sites %>%
  group_by(Order) %>%
  summarise(Sites = length(Order),
            Primary_vegetation = length(LUI[LUI == "Primary vegetation"]),
            Secondary_vegetation = length(LUI[LUI == "Secondary vegetation"]),
            Agriculture_Low = length(LUI[LUI == "Agriculture_Low"]),
            Agriculture_High = length(LUI[LUI == "Agriculture_High"]))
summary_LUI

summary_lat <- sites %>%
  group_by(Order) %>%
  summarise(Sites = length(Order),
            Tropical = length(Realm[Realm == "Tropical"]),
            NonTropical = length(Realm[Realm == "NonTropical"]),
            Percent_NonTropical = (NonTropical/Sites)*100,
            Percent_North = (length(Latitude[Latitude >= 0])/Sites)*100)
summary_lat

summary_climate <- sites %>%
  group_by(Order) %>%
  summarise(Sites = length(Order),
            AvgTemp = mean(avg_temp), 
            RangeofMean = max(avg_temp)-min(avg_temp),
            AvgSTA = mean(StdTmeanAnomaly),
            MaxTemp = max(max_temp),
            RangeofMax = max(max_temp)-min(max_temp),
            AvgSMTA = mean(StdTmaxAnomaly))
summary_climate

# save
write.csv(summary_LUI,"C:/Users/Kyra/Documents/GLITRS/Code/8_DataDistribution/summary_LUI.csv", row.names = TRUE)
write.csv(summary_lat,"C:/Users/Kyra/Documents/GLITRS/Code/8_DataDistribution/summary_lat.csv", row.names = TRUE)
write.csv(summary_climate,"C:/Users/Kyra/Documents/GLITRS/Code/8_DataDistribution/summary_climate.csv", row.names = TRUE)
