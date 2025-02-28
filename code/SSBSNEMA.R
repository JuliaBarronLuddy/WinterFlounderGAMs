library(ecodata)
library(dplyr)
library(here)
library(tidyverse)
library(mgcv)
library(ggplot2)
################################################################################

#SNEMA age data wrangling
snema_NAA <- snema_age[which(snema_age$AGE > 3),] #make ages 4 plus
snema_NAA <- snema_NAA[!snema_NAA$YEAR  < 1982,] #make years 1982 plus to match WAA data
snema_NAA <- subset(snema_NAA, select = c(SURVEY, SEASON, YEAR, AGE, NO_AT_AGE)) #select the columns we want to use
snema_NAA <- snema_NAA %>%
  mutate(Age_Group = ifelse(AGE >= 7, 7, AGE)) %>% #Create a new column to combine ages 7+
  group_by(YEAR, Age_Group) %>% #Group by Year and Age_Group
  dplyr::summarise(MEAN = mean(NO_AT_AGE, na.rm = TRUE), .groups = 'drop') %>% # Compute the average
  pivot_wider(
    names_from = Age_Group,
    values_from = MEAN,
    names_prefix = "Age"
  ) %>%
  mutate(MEAN = replace_na(rowMeans(across(starts_with("Age")), na.rm = TRUE), 0)) # Calculate mean of Age4, Age5, Age6, and Age7
names(snema_NAA)[names(snema_NAA) == '4'] <- 'Age4' #match column names of AGE to WAA data
names(snema_NAA)[names(snema_NAA) == '5'] <- 'Age5' 
names(snema_NAA)[names(snema_NAA) == '6'] <- 'Age6' 
names(snema_NAA)[names(snema_NAA) == '7'] <- 'Age7' #again, I know that this is 7+, not 7

#Troubleshooting NAA data wrangling
# SNEMA age data wrangling
snema_NAA <- snema_age[which(snema_age$AGE > 3),] # Keep ages 4+
snema_NAA <- snema_NAA[!snema_NAA$YEAR < 1982,] # Keep years 1982+
snema_NAA <- subset(snema_NAA, select = c(SURVEY, SEASON, YEAR, AGE, NO_AT_AGE)) # Select relevant columns

# Create Age_Group and calculate mean NO_AT_AGE by YEAR and Age_Group
snema_NAA <- snema_NAA %>%
  mutate(Age_Group = ifelse(AGE >= 7, 7, AGE)) %>%  # Combine ages 7+
  group_by(YEAR, Age_Group) %>%  
  summarise(MEAN = mean(NO_AT_AGE, na.rm = TRUE), .groups = 'drop') %>%  # Compute the average NO_AT_AGE by Year and Age_Group
  ungroup() %>%  # Ensure data is ungrouped before pivoting
  pivot_wider(  # Pivot to get separate columns for each Age_Group
    names_from = Age_Group,
    values_from = MEAN,
    names_prefix = "Age"
  ) %>%
  mutate(
    MEAN = rowMeans(across(starts_with("Age")), na.rm = TRUE) # Calculate mean of Age4, Age5, Age6, Age7, ensuring NA values are handled
  )

# Check the result
head(snema_NAA)

# Check and replace any rows where all age columns are NA (to prevent NaN in MEAN)
snema_NAA <- snema_NAA %>%
  mutate(
    MEAN = ifelse(is.na(MEAN), 0, MEAN)  # Replace NaN MEAN with 0 (or another value, as needed)
  )

# Rename columns to match WAA data
names(snema_NAA)[names(snema_NAA) == '4'] <- 'Age4' # Match AGE to WAA data
names(snema_NAA)[names(snema_NAA) == '5'] <- 'Age5' 
names(snema_NAA)[names(snema_NAA) == '6'] <- 'Age6' 
names(snema_NAA)[names(snema_NAA) == '7'] <- 'Age7' # Again, this is 7+, not 7


#WAA SNEMA data wrangling
waa_snema <- read.csv(here("data/WAA_snema.csv"))
names(waa_snema)[names(waa_snema) == 'Year'] <- 'YEAR' #match case of YEAR to NAA data
names(waa_snema)[names(waa_snema) == 'Age7+'] <- 'Age7'
names(waa_snema)[names(waa_snema) == 'Age7.'] <- 'Age7'
waa_snema <- subset(waa_snema, select = c(YEAR, Age4, Age5, Age6, Age7)) #get rid of premature fish ages (<4 years)
waa_snema <- waa_snema %>%
  mutate(MEAN = rowMeans(across(c(Age4, Age5, Age6, Age7)), na.rm = TRUE)) #find mean of all ages and include as separate column
waa_snema <- as_tibble(waa_snema)




##############################################################
get_ssb <- function(snema_NAA, waa_snema) {
  
  # Check column names
  print("Columns in snema_NAA:")
  print(names(snema_NAA))
  
  print("Columns in waa_snema:")
  print(names(waa_snema))
  
  # Select relevant columns from snema_NAA
  if (!all(c("YEAR", "Age4", "Age5", "Age6", "Age7", "MEAN") %in% names(snema_NAA))) {
    stop("Missing required columns in snema_NAA.")
  }
  SSB_stock <- snema_NAA[, c("YEAR", "Age4", "Age5", "Age6", "Age7", "MEAN")]
  
  # Debug: Check initial SSB_stock structure
  print("Initial SSB_stock:")
  print(head(SSB_stock))
  
  # Merge with waa_snema
  if (!"YEAR" %in% names(waa_snema)) {
    stop("waa_snema must have a 'YEAR' column.")
  }
  SSB_stock <- merge(SSB_stock, waa_snema, by = "YEAR")
  
  # Debug: Check merged data
  print("Merged SSB_stock:")
  print(head(SSB_stock))
  print(names(SSB_stock))  # Print column names
  
  # Check for NAs
  if (any(is.na(SSB_stock))) {
    print("NAs present in SSB_stock:")
    print(SSB_stock[rowSums(is.na(SSB_stock)) > 0, ])
  }
  
  # Ensure correct columns are present after merge
  age_columns <- c("Age4.x", "Age5.x", "Age6.x", "Age7.x", "MEAN.x")
  if (!all(age_columns %in% names(SSB_stock))) {
    stop(paste("Missing columns in SSB_stock:", paste(age_columns[!age_columns %in% names(SSB_stock)], collapse = ", ")))
  }
  
  # Calculate SSB using the correct column names
  SSB_stock$SSB <- rowSums(
    SSB_stock[, c("Age4.x", "Age5.x", "Age6.x", "Age7.x")] * SSB_stock$MEAN.x, 
    na.rm = TRUE
  )
  
  # Debug: Check SSB calculation
  print("Calculated SSB values:")
  print(SSB_stock$SSB)  # Print the SSB values
  
  # Keep only relevant columns
  SSB_stock <- SSB_stock[, c("YEAR", "SSB")]
  
  # Debug: Check final SSB stock
  print("Final SSB_stock:")
  print(head(SSB_stock))  # Print first few rows before aggregation
  
  # Aggregate SSB by YEAR
  SSB_aggregated <- aggregate(SSB ~ YEAR, data = SSB_stock, FUN = sum, na.action = na.omit)
  
  return(SSB_aggregated)
}
######################################################################

#call function
ssb_snema <- get_ssb(snema_NAA, waa_snema)

#Run the function and store the result
ssb_results_snema <- get_ssb(snema_NAA, waa_snema)

#Print the results
print(ssb_results_snema)

#Check the structure of the resulting data frame
str(ssb_results_snema)

#Check for any NA values
if (any(is.na(ssb_results_snema))) {
  print("There are NA values in the results.")
} else {
  print("No NA values in the results.")
}

#Print the first few rows of the result
print(head(ssb_results_snema))

#Get summary statistics
print(summary(ssb_results_snema))

#################################
#plotting code
ggplot(data = ssb_results_snema, aes(x = YEAR, y = SSB)) +
  geom_line(size = 1) +        
  geom_point(color = "red", size = 2) +         
  labs(title = "SNE/MA Spawning Stock Biomass (SSB)",
       x = "Year",
       y = "SSB") +
  theme_minimal() +                           
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  
    axis.title = element_text(size = 12),              
    axis.text = element_text(size = 10)                
  )

###########debugging!!!
# Check column names
print("Columns in snema_NAA:")
print(names(snema_NAA))

print("Columns in waa_snema:")
print(names(waa_snema))

# Select relevant columns from snema_NAA
if (!all(c("YEAR", "Age4", "Age5", "Age6", "Age7", "MEAN") %in% names(snema_NAA))) {
  stop("Missing required columns in snema_NAA.")
}

# Merge with waa_snema
if (!"YEAR" %in% names(waa_snema)) {
  stop("waa_snema must have a 'YEAR' column.")
}
SSB_stock <- merge(snema_NAA, waa_snema, by = "YEAR", all.x = TRUE)

# Debug: Check merged data
print("Merged SSB_stock:")
print(head(SSB_stock))
print(names(SSB_stock))  # Print column names

# Check for NAs
if (any(is.na(SSB_stock))) {
  print("NAs present in SSB_stock:")
  print(SSB_stock[rowSums(is.na(SSB_stock)) > 0, ])
}

# Ensure correct columns are present after merge
age_columns <- c("Age4.x", "Age5.x", "Age6.x", "Age7.x", "MEAN.y")
if (!all(age_columns %in% names(SSB_stock))) {
  stop(paste("Missing columns in SSB_stock:", paste(age_columns[!age_columns %in% names(SSB_stock)], collapse = ", ")))
}

# Debug: Check Age and MEAN values before multiplication
print("Age values before multiplication:")
print(head(SSB_stock[, c("Age4.x", "Age5.x", "Age6.x", "Age7.x")]))

print("MEAN.y values before multiplication:")
print(head(SSB_stock$MEAN.y))

# Calculate SSB
SSB_stock$SSB <- rowSums(
  SSB_stock[, c("Age4.x", "Age5.x", "Age6.x", "Age7.x")] * SSB_stock$MEAN.y, 
  na.rm = TRUE
)

# Debug: Check SSB calculation
print("Calculated SSB values:")
print(head(SSB_stock$SSB))  # Print first few SSB values

# Keep only relevant columns
SSB_stock <- SSB_stock[, c("YEAR", "SSB")]

# Debug: Check final SSB_stock
print("Final SSB_stock:")
print(head(SSB_stock))

# Aggregate SSB by YEAR
SSB_aggregated <- aggregate(SSB ~ YEAR, data = SSB_stock, FUN = sum, na.action = na.omit)

# Print final result
print("Final Aggregated SSB:")
print(head(SSB_aggregated))

# Check structure
str(SSB_aggregated)

# Check for NA values
if (any(is.na(SSB_aggregated))) {
  print("There are NA values in the results.")
} else {
  print("No NA values in the results.")
}

# Get summary statistics
print(summary(SSB_aggregated))
ssb_results_snema <- SSB_stock
ssb_snema <- ssb_results_snema







