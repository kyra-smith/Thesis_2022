# Shouldn't need this bit

# # build base map for fertiliser plot
# get_basemap <- function(){
#   
#   # download full basemap
#   base_map <- getMap(resolution = "high")
#   
#   # convert to correction projection
#   proj4string(base_map) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
#   
#   # return basemap
#   return(base_map)

###Functions
library(StatisticalModels)
library(predictsFunctions)
library(dplyr)


##Function for standardizing and centering predictor variables on very different scales, e.g. percNH
StdCenterPredictor <- function(x) {
  variable <- x
  sd <- sd(na.omit(variable))
  mean <- mean(na.omit(variable))
  variable.s <- (variable - mean)/sd
  return(variable.s)
}

BackTransformCentreredPredictor <- function(transformedX,originalX){
  
  sd <- sd(na.omit(originalX))
  mean <- mean(na.omit(originalX))
  backTX <- (transformedX * sd) + mean
  
  return(backTX)
  
}

RescaleAbundance <- function(sites){
  
  sites <- droplevels(sites)
  
  StudyMaxAbund <- suppressWarnings(tapply(
    X = sites$Total_abundance,INDEX = sites$SS,
    FUN = max,na.rm=TRUE))
  
  StudyMaxAbund[StudyMaxAbund == -Inf] <- NA
  
  AllMaxAbund <- StudyMaxAbund[match(sites$SS,names(StudyMaxAbund))]
  
  sites$Total_abundance_RS <- sites$Total_abundance/AllMaxAbund
  
  sites$LogAbund <- log(sites$Total_abundance_RS+0.01)
  
  return(sites)
}

##For a set of merged sites in predicts:
##This functions finds climate data, calculates scaled abundance values, logs abundance values,
##organizes the LandUse class, calculate standardise climate anomalies, and finds the percNH data
##In order to work, you must run the climate_percNH.R script and use the csv from that as the
## predicts_sites_filename arguement

organize_land_use_find_climate <- function(aggregated_sites) {
  
  sites <- aggregated_sites
  predicts_sites <- readRDS("Outputs/predicts_sites_info.rds")
  
  print("Organising Land Use")
  #####Organizing Land use
  sites$LandUse <- paste(sites$Predominant_land_use)
  sites$LandUse[(sites$LandUse=="Cannot decide")] <- NA
  sites$LandUse[(sites$LandUse=="Secondary vegetation (indeterminate age)")] <- NA
  sites$LandUse[(sites$LandUse=="Urban")] <- NA
  sites$LandUse <- factor(sites$LandUse)
  sites$LandUse <- relevel(sites$LandUse,ref="Primary vegetation")
  sites$Use_intensity[sites$Use_intensity=="Cannot decide"] <- NA
  sites$LandUse <- relevel(sites$LandUse, ref="Primary vegetation")
  sites$Use_intensity <- relevel(sites$Use_intensity, ref="Minimal use")
  
  ###Create agriculture intensity land class
  
  sites$Land <- paste(sites$LandUse)
  sites$Use_intensity <- as.character(sites$Use_intensity)
  sites$Use_intensity[is.na(sites$Use_intensity)] <- "Unknown"
  sites$Land <- interaction(sites$Land, sites$Use_intensity)
  sites$Land <- as.character(sites$Land)
  sites$Land[sites$Land=="Cropland.Intense use"] <- "High Agriculture"
  sites$Land[sites$Land=="Cropland.Light use"] <- "Low Agriculture"
  sites$Land[sites$Land=="Cropland.Minimal use"] <- "Low Agriculture"
  sites$Land[sites$Land=="Cropland.Unknown"] <- NA
  sites$Land[sites$Land=="Plantation forest.Light use"] <- "High Agriculture"
  sites$Land[sites$Land=="Plantation forest.Intense use"] <- "High Agriculture"
  sites$Land[sites$Land=="Plantation forest.Minimal use"] <- "Low Agriculture"
  sites$Land[sites$Land=="Plantation forest.Unknown"] <- NA
  sites$Land[sites$Land=="Pasture.Light use"] <- "Low Agriculture"
  sites$Land[sites$Land=="Pasture.Minimal use"] <- "Low Agriculture"
  sites$Land[sites$Land=="Pasture.Intense use"] <- "High Agriculture"
  sites$Land[sites$Land=="Pasture.Unknown"] <- NA
  sites$Land[sites$Land=="Primary vegetation.Intense use"] <- "PV"
  sites$Land[sites$Land=="Primary vegetation.Light use"] <- "PV"
  sites$Land[sites$Land=="Primary vegetation.Minimal use"] <- "PV"
  sites$Land[sites$Land=="Primary vegetation.Unknown"] <- "PV"
  sites$Land[sites$Land=="Intermediate secondary vegetation.Intense use"] <- "SV"
  sites$Land[sites$Land=="Intermediate secondary vegetation.Light use"] <- "SV"
  sites$Land[sites$Land=="Intermediate secondary vegetation.Minimal use"] <- "SV"
  sites$Land[sites$Land=="Intermediate secondary vegetation.Unknown"] <- "SV"
  sites$Land[sites$Land=="Mature secondary vegetation.Intense use"] <- "SV"
  sites$Land[sites$Land=="Mature secondary vegetation.Light use"] <- "SV"
  sites$Land[sites$Land=="Mature secondary vegetation.Minimal use"] <- "SV"
  sites$Land[sites$Land=="Mature secondary vegetation.Unknown"] <- "SV"
  sites$Land[sites$Land=="Young secondary vegetation.Intense use"] <- "SV"
  sites$Land[sites$Land=="Young secondary vegetation.Unknown"] <- "SV"
  sites$Land[sites$Land=="Young secondary vegetation.Minimal use"] <- "SV"
  sites$Land[sites$Land=="Young secondary vegetation.Light use"] <- "SV"
  sites$Land <- factor(sites$Land, levels = c("PV", "SV", "Low Agriculture", "High Agriculture"))
  
  print(table(sites$Land))
  
  ##Calculate log abundance
  sites$log_abundance <- log(sites$Total_abundance+1)
  
  
  ##Calculate Scaled abundance values
  n_sites <- length(unique(sites$SSBS))
  n_studies <- length(unique(sites$SS))
  print(paste("Calculating Scaled Abundance for",n_sites , "sites across",n_studies , "studies"))
  sites$Scaled_Abundance <- rep(3,1)
  for(i in 1:length(unique(sites$SS))) {
    s <- unique(sites$SS)[i]
    g <- sites[sites$SS==s, ]
    max <- max(g$Total_abundance)
    l <- g$Total_abundance / max
    
    sites[sites$SS==s, ]$Scaled_Abundance <- l
  }
  not_abundance <- which(sites$Scaled_Abundance==3)
  if(length(not_abundance)>0) sites[not_abundance, ]$Scaled_Abundance <- NA
  sites$log_scaled_abundance <- log(sites$Scaled_Abundance+0.01)
  
  ###Finding climate data for sites
  
  ##One Needs to be a character vector to match a factor
  predicts_sites$SSBS <- as.character(predicts_sites$SSBS)
  print("Finding Climate Data")
  site_info <- dplyr::select(predicts_sites, SSBS,climate_anomaly:SV_10000)
  sites <- left_join(sites, site_info)
  
  ##Combine PV and SV
  sites$NH_1000 <- sites$PV_1000 + sites$SV_1000
  sites$NH_3000 <- sites$PV_3000 + sites$SV_3000
  sites$NH_5000 <- sites$PV_5000 + sites$SV_5000
  sites$NH_10000 <- sites$PV_10000 + sites$SV_10000
  
  print("Standardising Climate Data")
  
  
  sites$std_climate_anomaly <- sites$climate_anomaly / sites$historic_sd
  sites$std_tmax_anomaly <- sites$tmax_anomaly / sites$historic_sd_tmax
  sites$std_tmax_quarter_anomaly <- sites$tmax_quarter_anomaly / sites$historic_sd_tmax
  
  ##Standardising and centering variables on very different scales
  sites$NH_1000.s <- StdCenterPredictor(sites$NH_1000)
  sites$NH_3000.s <- StdCenterPredictor(sites$NH_3000)
  sites$NH_5000.s <- StdCenterPredictor(sites$NH_5000)
  sites$NH_10000.s <- StdCenterPredictor(sites$NH_10000)
  
  print("Done")
  
  return(sites)
}
max_quarter_fast <- function(x) {
  
  mean(Rfast::Sort(x, TRUE)[1:3])
  
}



vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

# edited version of StatisticalModels function for use with glmmTMB model

PredictGLMERRandIter_TMB <- function(model,data,nIters=1000){
  
  # stopifnot((class(model)[1] == "lmerMod") | (class(model)[1] == "glmerMod"))
  
  preds <- sapply(X = 1:nIters,FUN = function(i){
    
    mm<-model.matrix(terms(model),data)
    
    if(ncol(mm)>length(lme4::fixef(model)[1]$cond)){
      mm <- mm[,-which(!(names(mm[1,]) %in% names(
        lme4::fixef(model))))]
    }
    
    fe.draw <- mvrnorm(n = 1,mu = fixef(model)[1]$cond,Sigma = vcov(model)[1]$cond)
    
    y <- mm %*% fe.draw
    
  })
  
  return(preds)
}

### function for extracting covariance matrix and calculating uncertainties
iterate_covar <- function(i, model, prediction_data, factor_no_1, factor_no_2, fixed_effect_no, neg_binom){
  
  # if the family is negative binomial zero inflated, need to index for the fixed effects
  if(neg_binom == TRUE){
    
    # extract fixed effects from covariance matrix 
    coefs <- mvrnorm(n = 1, mu = fixef(object = model)[1]$cond, Sigma = vcov(object = model)[1]$cond)
  }
  
  # if not negative binomial, can just take the fixed effects straight from the model
  else{
    coefs <- mvrnorm(n = 1, mu = fixef(object = model), Sigma = vcov(object = model))
  }
  
  # then extract the terms from the model matrix
  mm <- model.matrix(terms(model), prediction_data)
  
  # set up for loop for adjusting for seperate taxonomic class
  counter <- 1
  subset_predictions <- c()
  
  if(fixed_effect_no == 2){
    if (ncol(mm) > length(lme4::fixef(model))) {
      mm <- mm[, -which(!(names(mm[1, ]) %in% names(lme4::fixef(model))))]
    }
    y <- mm %*% coefs
    
    # run for loop to take exponential, and adjust as rescaled percentage for each sample
    for(i in 1:(factor_no_1)){
      if(i < 6){
        subset_predictions <- c(subset_predictions, (exp(y[counter:(counter+(factor_no_2-1))]) / exp(y)[counter] * 100))
      }
      
      else{
        subset_predictions <- c(subset_predictions, (exp(y[counter:(counter+(factor_no_2-1))]) / exp(y)[counter] * 100))
      }
      # step up the counter to move along the vector
      counter <- counter + factor_no_2
    }
  }
  
  # if fixed effect is 1, calculate effect sizes relative to the first value (primary vegetation)
  else{
    y <- mm %*% coefs
    subset_predictions <- (exp(y)/exp(y)[1])*100
  }
  
  # return the vector of adjusted values for that sample
  return(subset_predictions)
}


## Predictions ##

 predict_effects <- function(model, 
                             iterations, 
                             model_data, 
                             response_variable, 
                             factor_number_1,
                             factor_number_2,
                             fixed_number,
                             fixed_column,
                             neg_binom){

# print the fixed effects for the model as a reference for predictions
print(fixef(model))

# commented this bit because I know how many fixed effects I have ##

# if number of fixed effects = 1, create a dataframe from the model_data with one fixed and response
if(fixed_number == 1){
  unique_data <- data.frame(factor(as.character(levels(model_data[,fixed_column])),
                                   levels = levels(model_data[,fixed_column])), 0) %>%
    unique()
  
  # assigned the correct column names according to fixed effect and response arguments
  colnames(unique_data) <- c(fixed_column, response_variable)
}

# # if fixed effect is not 1 (2), build dataframe for the 2 fixed effects
else{
  
  # set up the dataframe of unique factors and reorder by the model data
  unique_data <- data.frame(factor(as.character(model_data[, fixed_column[1]]),
                                   levels = levels(model_data[, fixed_column[1]])),
                            factor(as.character(model_data[, fixed_column[2]]),
                                   levels = levels(model_data[, fixed_column[2]])), 0) %>%
    unique()
  
  # reorder the data frame by factor order
  colnames(unique_data) <- c(fixed_column[1], fixed_column[2], response_variable)
  col_1 <- unique_data %>% pull(fixed_column[1])
  col_2 <- unique_data %>% pull(fixed_column[2])
  unique_data <- unique_data[with(unique_data, order(col_1, col_2)), ]
  
  # run iterations of extraction of coefficients
  preds.emp <- sapply(X = 1:iterations, iterate_covar, model, prediction_data = unique_data, factor_no_1 = factor_number_1, factor_no_2 = factor_number_2, fixed_effect_no = fixed_number, neg_binom)
  
  # extract the median, upper interval, and lower interval for samples
  preds.emp.summ <- data.frame(Median = apply(X = preds.emp, MARGIN = 1, FUN = median),
                               Upper = apply(X = preds.emp, MARGIN = 1, FUN = quantile, probs = 0.975),
                               Lower = apply(X = preds.emp, MARGIN = 1, FUN = quantile, probs = 0.025))
  
  # bind prediction data back onto the median, upper, and lower intervals, and adjust as percentage change - note factor order
  fin_conf <- cbind(unique_data, preds.emp.summ) %>%
    mutate(Median = Median - 100) %>%
    mutate(Upper = Upper - 100) %>%
    mutate(Lower = Lower - 100)
  
  if(fixed_number == 2){
    # convert the reference factor (primary vegetation) to NA so doesn't plot error bar
    fin_conf$Upper[fin_conf$LUI == "Primary vegetation"] <- NA
    fin_conf$Lower[fin_conf$LUI == "Primary vegetation"] <- NA
    
    # print final dataframe before plotting
    print(fin_conf)
    .GlobalEnv$fin_conf <- fin_conf
    
    # construct the plot for that biodiversity metric - note the y axis label needs to be specific for that metric
    output_plot <- ggplot(fin_conf) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(aes(x = fin_conf[,1], y = Median, colour = fin_conf[,2]), position = position_dodge(0.5), size = 2) +
      #geom_errorbar(aes(x = fin_conf[,1], ymin = Lower, ymax = Upper, colour = fin_conf[,2]), position = position_dodge(0.5), size = 1, width = 1/(factor_number_2)) +
      geom_errorbar(aes(x = fin_conf[,1], ymin = Lower, ymax = Upper, colour = fin_conf[,2]), position = position_dodge(0.5), size = 0.5, width = 0.5) +
      scale_colour_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00"), name = "Land Use Intensity") +
      #scale_y_continuous(name = gsub(x = paste(response_variable, "difference (%)", sep = " "), pattern = "_", replacement = " ")) +
      theme_bw() +
      theme(panel.grid = element_blank())
  }
  
  # if number of fixed effects =  1 (and single combined LUI) convert LUI into two columns for ggplot 
  else{
    # set up vectors for intensity and land-use types - can be moved as an argument?
    land_use_intensity <- c("MU", "LU", "IU")
    land_use_type <- c("PV", "MSV", "ISV", "YSV", "PF", "P", "C", "U")
    fin_conf$veg_type <- fin_conf[,fixed_column]
    
    # create empty vectors, check for each intensity in veg type, and assign to a new vector
    assign_intensity <- c()
    intensity <- c()
    for(j in 1:length(fin_conf$veg_type)){
      for(i in 1:length(land_use_intensity)){
        assign_intensity <- grep(land_use_intensity[i], fin_conf$veg_type[j])
        if(length(assign_intensity > 0)){
          intensity <-  c(intensity, land_use_intensity[i])
        }
      }
    }
    
    # assign the correctly sorted intensity to a column of dataframe 
    fin_conf$intensity <- intensity 
    
    # amend columns - veg_type
    fin_conf$veg_type <- gsub("MU", "", fin_conf$veg_type)
    fin_conf$veg_type <- gsub("LU", "", fin_conf$veg_type)
    fin_conf$veg_type <- gsub("IU", "", fin_conf$veg_type)
    
    # order factors for intensity and land use type
    fin_conf$intensity <- factor(fin_conf$intensity, levels = land_use_intensity,
                                 labels = c("Minimal use", "Light use", "Intense use"))
    fin_conf$veg_type <- factor(fin_conf$veg_type, levels = land_use_type, 
                                labels = c("Primary", "MSV", "ISV", "YSV", "Plantation", "Pasture", "Cropland", "Urban"))
    
    # print final dataframe before plotting
    print(fin_conf)
    
    # construct the plot for that biodiversity metric - note the y axis label needs to be specific for that metric
    output_plot <- ggplot(fin_conf) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(aes(x = veg_type, y = Median, shape = intensity, colour = veg_type), position = position_dodge(0.75), size = 2.5) +
      geom_errorbar(aes(x = veg_type, ymin = Lower, shape = intensity, ymax = Upper, colour = veg_type), position = position_dodge(0.75), size = 1, width = 0.25) +
      scale_colour_manual(values = c("#E69F00", "#009E73", "#F0E442","#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")) +
      scale_y_continuous(name = gsub(x = paste(response_variable, "difference (%)", sep = " "), pattern = "_", replacement = " ")) +
      xlab(NULL) +
      labs(shape = NULL) +
      guides(colour = FALSE) +
      theme_bw() +
      theme(legend.position = c(0.15, 0.2), legend.background = element_rect(colour = "black"), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5, size = 13, colour = c("#E69F00", "#009E73", "#F0E442","#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")))
  }
  
  #   # return the plot for that biodiversity metric
  return(output_plot)
  
}
 }
 
 
# function for predicting continuous variable
predict_continuous <- function(model,
                               model_data, 
                               response_variable,
                               categorical_variable,
                               continuous_variable,
                               continuous_transformation,
                               random_variable,
                               colour_palette){
  
  # set up the prediction dataframe
  prediction_data <- model_data[, c(response_variable, 
                                    random_variable[1], 
                                    random_variable[2], 
                                    random_variable[3], 
                                    categorical_variable[1], 
                                    continuous_variable[1])]
  
  # remove any incomplete rows (NAs) from the prediction data
  prediction_data <- prediction_data[complete.cases(prediction_data),]
  
  # predict the values for the model
  y_value <- c(StatisticalModels::PredictGLMER(model, data = prediction_data, se.fit = TRUE, seMultiplier = 1.96))[[1]]
  y_value_plus <- c(StatisticalModels::PredictGLMER(model, data = prediction_data, se.fit = TRUE, seMultiplier = 1.96))[[2]]
  y_value_minus <- c(StatisticalModels::PredictGLMER(model, data = prediction_data, se.fit = TRUE, seMultiplier = 1.96))[[3]]
  
  # bind the predicted values to the prediction data
  bound_values <- data.frame(cbind(prediction_data,
                                   y_value, 
                                   y_value_plus, 
                                   y_value_minus, 
                                   metric = response_variable,
                                   continuous_transformation(prediction_data[, continuous_variable[1]])))
  
  # rename last column after transformation
  colnames(bound_values)[ncol(bound_values)] <- paste(continuous_variable[1], "transform", sep = "_")
  
  # rename the response variable column "response_variable" for later plot function
  bound_values <- bound_values %>%
    rename("response_variable" = all_of(response_variable))
  
  # print the final dataframe of predicted values
  print(bound_values)
  
}

# function for plotting the output from a continuous variable glmer - predict_continuous()
# need to amend to read in categorical variable rather than zone/order
plot_fert_response <- function(data_set, categorical_variable){
  plot_obj <- ggplot(data_set) +
    geom_line(aes_string(x = "fert_transform", y = "y_value", colour = categorical_variable, linetype = "significant"), size = 0.7) +
    geom_ribbon(aes_string(x = "fert_transform", ymin = "y_value_minus", ymax = "y_value_plus", fill = categorical_variable), alpha = 0.4) +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 2.39794, 2.69897, 3, 3.39794), labels = c(0.1, 1, 10, 100, 250, 500, 1000, 2500)) +
    scale_linetype_manual("", values = c("solid", "dashed"), labels = c("Non-significant", "Significant")) +
    theme_bw() +
    theme(panel.grid = element_blank())
  return(plot_obj)
}