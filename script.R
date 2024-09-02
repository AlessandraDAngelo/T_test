# R script for calculating the T-test statistic for experimental marine methane concentration and isotope ratio.
# D'Angelo A. and Loose B. (2020) - University of Rhode Island, Graduate School of Oceanography

# Load dataframe "data"

library(dplyr) 

# Splitting file and creating character names to subset the df to calculate the standard error of linear regression for each sample
Bag_num_split <- split(data, data$bag.Num, drop = FALSE)
new_names <- as.character(unique(MOSAiC.methox.results$bag.Num))

for (i in 1:length(Bag_num_split)) {
  assign(new_names[i], Bag_num_split[[i]])
  filename <- paste(new_names[i], ".csv")  
  write.csv(Bag_num_split[[i]], filename, row.names = FALSE)  
}

my_files <- list.files(path="", pattern=".csv", full.names=TRUE, recursive=FALSE)

bag <- as.list(rep("", length(my_files)))

# T-test for methane concentrations
for(i in 1:length(my_files))
{
  # Defining bag[i]
  bag[[i]] <- read.csv(my_files[i])
  # Create a column with the intercept value of methane moles (tot_mass) versus dates
  bag[[i]]$intercept_T <- lm(bag[[i]]$tot_mass~bag[[i]]$valueDate)$coef[1]
  # Create a column with the slope value of methane moles (tot_mass) versus dates
  bag[[i]]$slope_T <- lm(bag[[i]]$tot_mass~bag[[i]]$valueDate)$coef[2]
  # Use the absolute value for the statistic test T
  bag[[i]]$slope_T_abs <- abs(lm(bag[[i]]$tot_mass~bag[[i]]$valueDate)$coef[2])
  # Create a column with expected values. Use + instead of - to have subtraction.
  bag[[i]]$yhat <- (as.numeric(bag[[i]]$valueDate * bag[[i]]$slope_T)+(as.numeric(bag[[i]]$intercept_T)))
  # Create a column for (y-yhat)^2
  bag[[i]]$DF <- (sum(table(as.numeric(bag[[i]]$tot_mass)))-2)
  bag[[i]]$y_yhat <- as.numeric((bag[[i]]$tot_mass - bag[[i]]$yhat)^2)/bag[[i]]$DF
  # Create a column for (x-xbar)^2
  bag[[i]]$mean_valueDate <- mean(bag[[i]]$valueDate)
  bag[[i]]$x_xbar <- (bag[[i]]$valueDate - bag[[i]]$mean_valueDate)^2
  # Standard error
  a <- sqrt(sum(na.omit(bag[[i]]$y_yhat)))
  b <- sqrt(sum(na.omit(bag[[i]]$x_xbar)))
  bag[[i]]$SE <- a/b
  # Confidence interval 0.95
  P <- 0.95
  # Critical value 
  bag[[i]]$t <- qt(0.95, bag[[i]]$DF)
  # Statistic test T
  bag[[i]]$T <- bag[[i]]$slope_T_abs/bag[[i]]$SE
  bag[[i]]$Null.hypothesis <- ifelse(bag[[i]]$T > bag[[i]]$t, "Reject", "Accept")
  write.csv(bag[[i]], paste("Bag", bag[[i]]$bag.Num[1] , ".csv"), row.names = FALSE) 
}

# T-test for methane isotopic ratio

for(i in 1:length(my_files))
{
  # Defining bag[i]
  bag[[i]] <- read.csv(my_files[i])
  # Create a column with the intercept value of methane isotope ratio (iso.data.calibrated) versus dates
  bag[[i]]$intercept_T <- lm(bag[[i]]$iso.data.calibrated ~bag[[i]]$valueDate)$coef[1]
  # Create a column with the slope value of methane isotope ratio (iso.data.calibrated) versus dates
  bag[[i]]$slope_T <- lm(bag[[i]]$iso.data.calibrated~bag[[i]]$valueDate)$coef[2]
  # Use the absolute value for the statistic test T
  bag[[i]]$slope_T_abs <- abs(lm(bag[[i]]$iso.data.calibrated~bag[[i]]$valueDate)$coef[2])
  # Create a column with expected values. Use + instead of - to have subtraction. In Excel, we had subtraction using -.
  bag[[i]]$yhat <- (as.numeric(bag[[i]]$valueDate * bag[[i]]$slope_T)+(as.numeric(bag[[i]]$intercept_T)))
  # Create a column for (y-yhat)^2
  bag[[i]]$DF <- (sum(table(as.numeric(bag[[i]]$iso.data.calibrated)))-2)
  bag[[i]]$y_yhat <- as.numeric((bag[[i]]$iso.data.calibrated - bag[[i]]$yhat)^2)/bag[[i]]$DF
  # Create a column for (x-xbar)^2
  bag[[i]]$mean_valueDate <- mean(bag[[i]]$valueDate)
  bag[[i]]$x_xbar <- (bag[[i]]$valueDate - bag[[i]]$mean_valueDate)^2
  # Standard error
  a <- sqrt(sum(na.omit(bag[[i]]$y_yhat)))
  b <- sqrt(sum(na.omit(bag[[i]]$x_xbar)))
  bag[[i]]$SE <- a/b
  # Confidence interval 0.95
  P <- 0.95
  # Critical value 
  bag[[i]]$t <- qt(P, bag[[i]]$DF)
  # Statistic test T
  bag[[i]]$T <- bag[[i]]$slope_T_abs/bag[[i]]$SE
  bag[[i]]$Null.hypothesis <- ifelse(bag[[i]]$T > bag[[i]]$t, "Reject", "Accept")
  write.csv(bag[[i]], paste("Bag", bag[[i]]$bag.Num[1] , ".csv"), row.names = FALSE) 
}
