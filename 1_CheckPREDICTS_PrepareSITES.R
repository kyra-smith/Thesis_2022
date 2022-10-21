# get set up
  
# set up directories
dataDir <- "C:/Users/Kyra/Documents/GLITRS/Data/"
outDir <- "C:/Users/Kyra/Documents/GLITRS/Code/1_CheckPrepareData/"
if(!dir.exists(outDir)) dir.create(outDir)

# Load required libraries
library(predictsFunctions)
library(dplyr)
library(ggplot2)
source("C:/Users/Kyra/Documents/GLITRS/Data/0_Functions.R")

packages_plot <- c("patchwork", "dplyr", "yarg", "lme4", "gt", "broom.mixed", "MASS","webshot")
suppressWarnings(suppressMessages(lapply(packages_plot, require, character.only = TRUE)))

# Set the path to copy of the database
predicts.path <- paste0(dataDir,"database.rds")

# Read in the PREDICTS data
predicts <- ReadPREDICTS(predicts.path)
# 3250404 obs. of 67 variables

# Select only data for insects
predicts <- predicts[(predicts$Class=="Insecta"),]
# 935078 obs. of 67 variables

# Correct effort-sensitive abundance measures (assumes linear relationship between effort and recorded abundance)
predicts <- CorrectSamplingEffort(diversity = predicts)
# Correcting 0 missing sampling effort values
# Re-scaling sampling effort
# Correcting 870378 values for sensitivity to sampling effort # matches Outhwaite et al.

# check diversity metrics
table(predicts$Diversity_metric)

# insects should not have diversity metric "percent cover", this is a mistake in the database
# remove those entries that are the problem
predicts <- predicts[!predicts$Diversity_metric == "percent cover", ]
# 934845 obs. of 67 variables


# MergeSites

# Merge sites that have the same coordinates, from the same study and same taxonomic family (e.g. multiple traps on a single transect)
predicts <- MergeSites(diversity = predicts)

nrow(predicts)
# 826292 obs. of 67 variables

# # predicts.complete
# 
# # remove entries where land use or use intensity info is missing
# predicts.complete <- droplevels(predicts[(predicts$Predominant_land_use!="Cannot decide"),])
# predicts.complete <- droplevels(predicts.complete[(predicts.complete$Use_intensity!="Cannot decide"),])
# 
# nrow(predicts.complete)
# # 779,912 obs. of 67 variables # matches Outhwaite et al
# 
# # remove entries where Order is missing
# predicts.complete <- droplevels(predicts.complete[(predicts.complete$Order!=""),])
# 
# nrow(predicts.complete)
# # 779,670 obs. of 67 variables
# 
# species <- unique(predicts.complete[,c('Order','Taxon_name_entered')])
# 
# order.count <- tapply(X = species$Taxon_name_entered,
#                       INDEX = species$Order,
#                       FUN = function(sp) length(unique(sp)))
# order.count


# prepare PREDICTS

# remove entries without Order
predicts <- droplevels(predicts[(predicts$Order!=""),])
# 826016 obs. of 67 variables

# check
table(predicts$Order)

# recode Isoptera to Blattodea and Phthiraptera to Psocodea
# Isoptera = infraorder within Blattodea; Phtiraptera = parvorder within Psocodea (Catalogue of Life)

predicts$Order <- recode_factor(predicts$Order,
                                Isoptera = "Blattodea",
                                Phthiraptera = "Psocodea")
# 826016 obs. of 67 variables

# summarize predicts statistics by order
predicts_summary_all <- predicts %>%
  group_by(Order) %>%
  summarise(Count_Sites=length(Order),
            Unique_Sites = n_distinct(SSBS),
            Unique_Species = n_distinct(Taxon_name_entered)) %>% arrange(desc(Unique_Sites))
predicts_summary_all


# # optional histogram
# hist<-ggplot(summary_orders, aes(x=CountRecords))+  
#   geom_histogram()+
#   labs(x=bquote(Count~Records), y=bquote(Orders))
# hist

#remove orders with less than 600 unique sites (SSBS)
# these are the orders: Psocodea, Dermaptera, Archaeognatha Mantodea, Odonata, Embioptera, Ephemeroptera, Mecoptera, Zoraptera, Siphonaptera, Zygentoma, Phasmida, Thysanoptera, Trichoptera, Blattodea, Neuroptera
# remaining orders are: Hymenoptera, Coleoptera, Lepidoptera, Diptera, Orthoptera, Hemiptera

predicts <- predicts %>% filter(Order %in% c("Hymenoptera", "Coleoptera", "Lepidoptera", "Diptera", "Orthoptera", "Hemiptera")) %>% droplevels()

# 814738 obs. of 67 variables

# summarize predicts statistics by order
predicts_summary <- predicts %>%
  group_by(Order) %>%
  summarise(Count_Sites=length(Order),
            Unique_Sites = n_distinct(SSBS),
            Unique_Species = n_distinct(Taxon_name_entered))%>% arrange(desc(Unique_Sites))
predicts_summary

# convert Order to a "factor"
predicts$Order <- as.factor(predicts$Order)

# save as csv
write.csv(predicts_summary,"C:/Users/Kyra/Documents/GLITRS/Code/1_CheckPrepareData/predicts_summary.csv", row.names = TRUE)


# Split predicts into separate data frames according to insect Order 

# use split function to split the predicts data frame into 6 data frames (1/Order)
OrderName <- paste0("",predicts$Order)

by_Order <- split(predicts,OrderName)

# extract data frames from list into global environment
list2env(by_Order,globalenv())

# Calculate site metrics of diversity for each order, include extra columns:
  # Predominant_land_use
  # SSB
  # SSBS
  # Biome
  # Order

# droplevels() drops unused factor levels. This is particularly useful if we want to drop factor levels that are no longer used due to subsetting a vector or a data frame (as we did with split()). 

Coleoptera <- droplevels(Coleoptera)
Coleoptera <- SiteMetrics(diversity = Coleoptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order"))
Diptera <- droplevels(Diptera)
Diptera <- SiteMetrics(diversity = Diptera,
                       extra.cols = c("Predominant_land_use",
                                      "SSB","SSBS", "Biome","Order"))
Hemiptera <- droplevels(Hemiptera)
Hemiptera <- SiteMetrics(diversity = Hemiptera,
                         extra.cols = c("Predominant_land_use",
                                        "SSB","SSBS", "Biome","Order"))
Hymenoptera <- droplevels(Hymenoptera)
Hymenoptera <- SiteMetrics(diversity = Hymenoptera,
                           extra.cols = c("Predominant_land_use",
                                          "SSB","SSBS", "Biome","Order"))
Lepidoptera <- droplevels(Lepidoptera)
Lepidoptera <- SiteMetrics(diversity = Lepidoptera,
                           extra.cols = c("Predominant_land_use",
                                          "SSB","SSBS", "Biome","Order"))
Orthoptera <- droplevels(Orthoptera)
Orthoptera <- SiteMetrics(diversity = Orthoptera,
                          extra.cols = c("Predominant_land_use",
                                         "SSB","SSBS", "Biome","Order"))


# merge all sites_Order data frames into one called "sites"
# merge using rbind()
sites <- rbind(Coleoptera,Diptera,Hemiptera,Hymenoptera,Lepidoptera,Orthoptera)


# First, we will rearrange the land-use classification a bit
# rename "Predominant_land_use" to "LandUse"
sites$LandUse <- paste(sites$Predominant_land_use)

# Drop classification where land use could not be identified
sites$LandUse[(sites$LandUse=="Cannot decide")] <- NA

# drop classification where use intensity could not be identified
sites$Use_intensity[sites$Use_intensity=="Cannot decide"] <- NA

# Now make the variable a factor, and set the reference level to primary vegetation
sites$LandUse <- factor(sites$LandUse)
sites$LandUse <- relevel(sites$LandUse,ref="Primary vegetation")

# combine LandUse (LU) and Use Intensity (UI) into new variable Land Use Intensity (LUI)
sites$LUI <- paste0(sites$LandUse,'_',sites$Use_intensity)
sites$LUI[grep("NA",sites$LUI)] <- NA # where "NA" appears in the UI field, drop classification

# recode according to land use and use intensity combinations
sites$LUI <- dplyr::recode(sites$LUI,
                           'Primary vegetation_Minimal use' = 'Primary vegetation',
                           'Cropland_Light use' = 'Agriculture_High',
                           'Secondary vegetation (indeterminate age)_Minimal use' = 'Secondary vegetation',
                           'Urban_Light use' = 'Urban',
                           'Secondary vegetation (indeterminate age)_Light use' = 'Secondary vegetation',
                           'Cropland_Intense use' = 'Agriculture_High',
                           'Cropland_Minimal use' = 'Agriculture_Low',
                           'Pasture_Light use' = 'Agriculture_Low',
                           'Pasture_Minimal use' = 'Agriculture_Low',
                           'Intermediate secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Mature secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Secondary vegetation (indeterminate age)_Intense use' = 'Secondary vegetation',
                           'Pasture_Intense use' = 'Agriculture_High',
                           'Urban_Minimal use' = 'Urban',
                           'Primary vegetation_Light use' = 'Primary vegetation',
                           'Young secondary vegetation_Light use' = 'Secondary vegetation',
                           'Urban_Intense use' = 'Urban',
                           'Primary vegetation_Intense use' = 'Primary vegetation',
                           'Young secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Mature secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Plantation forest_Minimal use' = 'Agriculture_Low',
                           'Plantation forest_Intense use' = 'Agriculture_High',
                           'Young secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Plantation forest_Light use' = 'Agriculture_High',
                           'Mature secondary vegetation_Light use' = 'Secondary vegetation',
                           'Intermediate secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Intermediate secondary vegetation_Light use' = 'Secondary vegetation')

# 11410 obs. of 25 variables

# ...not sure, but it looks like we are re-classifying all secondary vegetation as "Light use"
# but why after we have already re-coded the land use type and intensity at all sites?
sites$Use_intensity[((sites$LandUse=="Mature secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"
sites$Use_intensity[((sites$LandUse=="Intermediate secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"
sites$Use_intensity[((sites$LandUse=="Young secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"

# remove the urban sites and sites that are NA in LUI
sites <- sites[!sites$LUI == "Urban", ]
sites <- sites[!is.na(sites$LUI), ]

# 9461 obs. of 25 variables

sites <- droplevels(sites)

# transform abundance values 
sites$LogAbund <- log(sites$Total_abundance+1)

# Remove sites without coordinates
sites <- sites[!is.na(sites$Latitude), ]

# sites: 9455 obs. of 26 variables

# create a new variable designating sites as Tropical or Non-tropical
# assign new variable for tropical/temperate, convert to factor, and filter out NA
sites$Realm <- ifelse(sites$Latitude >= -23.5 & sites$Latitude <= 23.5, "Tropical", "NonTropical")
sites$Realm <- factor(sites$Realm, levels = c("NonTropical", "Tropical"))
sites <- sites %>%
  filter(!is.na(Realm))

# save as csv
write.csv(sites_summary,"C:/Users/Kyra/Documents/GLITRS/Code/1_CheckPrepareData/sites_summary.csv", row.names = TRUE)

# summarize sites by order
# summarize sites_all
sites_summary <- sites %>%
  mutate(LogAbund = ifelse(is.na(LogAbund), 0, LogAbund)) %>% # replace NA values with 0 for ease of summarising
  group_by(Order) %>%
  summarise(Unique_Sites = length(Order),
            Primary_vegetation = length(LUI[LUI == "Primary vegetation"]),
            Secondary_vegetation = length(LUI[LUI == "Secondary vegetation"]),
            Agriculture_Low = length(LUI[LUI == "Agriculture_Low"]),
            Agriculture_High = length(LUI[LUI == "Agriculture_High"]),
            Tropical = length(Realm[Realm == "Tropical"]),
            NonTropical= length(Realm[Realm == "NonTropical"]),
            Abundance = sum(LogAbund>0) ,
            SpeciesRichness = sum(Species_richness>0)) %>%
ungroup() %>%
gt() %>%
  tab_spanner(
    label = "Land-use-intensity",
    columns = c(Primary_vegetation, Secondary_vegetation, Agriculture_Low, Agriculture_High)
  ) %>%
  tab_spanner(
    label = "Latitudinal Realm",
    columns = c(Tropical,NonTropical)
  ) %>%
  tab_spanner(
    label = "Diversity Metric",
    columns = c(Abundance,SpeciesRichness)
  )  %>% 
  cols_align(
    align = "center",
    columns = c(Order, Unique_Sites, Primary_vegetation, Secondary_vegetation, Agriculture_Low,Agriculture_High,Tropical,NonTropical,Abundance,SpeciesRichness)
    )%>%
  cols_label(
    Order = md("Order"),
    Unique_Sites = md("Sites"),
    Primary_vegetation = md("Primary vegetation"),
    Secondary_vegetation = md("Secondary vegetation"),
    Agriculture_Low = md("Low-intensity agriculture"),
    Agriculture_High = md("High-intensity agriculture"),
    Tropical = md("Tropical"),
    NonTropical = md("Non-tropical"),
    Abundance = md("Abundance"),
    SpeciesRichness = md("Species richness")
  )

sites_summary

gtsave(sites_summary,"C:/Users/Kyra/Documents/GLITRS/Code/1_CheckPrepareData/sites_summary.png")

# save the prepared dataset
saveRDS(object = sites,file = paste0(outDir,"PREDICTSSiteData.rds")) 

## plot points of sites ##

# plot the raster in ggplot
map.world <- map_data('world')

# map of sites, coloured by order
p_points_colour <-ggplot() +
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= "grey", colour="grey", size=0.2) +
  geom_point(data = sites, aes(x = Longitude, y = Latitude, colour = factor(Order)), shape = 20) +
  theme(axis.title = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())


# save plot
ggsave(filename = paste0(outDir, "/PREDICTS_points_map_colour.pdf"), height = 4, width = 8)
