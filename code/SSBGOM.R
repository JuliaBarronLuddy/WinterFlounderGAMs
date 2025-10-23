library(ecodata)
library(dplyr)
library(here)
library(tidyverse)
library(mgcv)
library(ggplot2)
################################################################################

#GOM age data wrangling
gom_age <- read.csv(here("data/gom_age.csv"))
gom_NAA <- gom_age[which(gom_age$AGE > 3),] #make ages 4 plus
gom_NAA <- gom_NAA[!gom_NAA$YEAR  < 1982,] #make years 1982 plus to match WAA data
gom_NAA <- subset(gom_NAA, select = c(SURVEY, SEASON, YEAR, AGE, NO_AT_AGE)) #select the columns we want to use
gom_NAA <- gom_NAA %>%
  mutate(Age_Group = ifelse(AGE >= 7, 7, AGE)) %>% #Create a new column to combine ages 7+
  group_by(YEAR, Age_Group) %>% #Group by Year and Age_Group
  dplyr::summarise(MEAN = mean(NO_AT_AGE, na.rm = TRUE), .groups = 'drop') %>% # Compute the average
  pivot_wider(
    names_from = Age_Group,
    values_from = MEAN,
    names_prefix = "Age"
  ) %>%
  mutate(MEAN = replace_na(rowMeans(across(starts_with("Age")), na.rm = TRUE), 0)) # Calculate mean of Age4, Age5, Age6, and Age7
names(gom_NAA)[names(gom_NAA) == '4'] <- 'Age4' #match column names of AGE to WAA data
names(gom_NAA)[names(gom_NAA) == '5'] <- 'Age5' 
names(gom_NAA)[names(gom_NAA) == '6'] <- 'Age6' 
names(gom_NAA)[names(gom_NAA) == '7'] <- 'Age7' #again, I know that this is 7+, not 7

#WAA GOM data wrangling
waa_gom <- read.csv(here("data/WAA_gom.csv"))
waa_gom <- subset(waa_gom, select = c(SURVEY, YEAR, AGE, MEAN)) #select the columns we want to use)) 
names(waa_gom)[names(waa_gom) == 'MEAN'] <- 'Weight'
waa_gom <- waa_gom[which(waa_gom$AGE > 3),]

waa_gom <- waa_gom %>%
  group_by(YEAR, AGE) %>%
  summarise(Weight = mean(Weight), .groups = "drop")

waa_gom <- waa_gom %>%
mutate(Age_Group = ifelse(AGE >= 7, 7, AGE)) %>% #Create a new column to combine ages 7+
  group_by(YEAR, Age_Group) %>% #Group by Year and Age_Group
  dplyr::summarise(MEAN = mean(Weight, na.rm = TRUE), .groups = 'drop') %>% # Compute the average
  pivot_wider(
    names_from = Age_Group,
    values_from = MEAN,
    names_prefix = "Age"
  ) %>%
  mutate(MEAN = replace_na(rowMeans(across(starts_with("Age")), na.rm = TRUE), 0))


print(waa_gom)


########### GET SSB ###########

get_ssb <- function(gom_NAA, waa_gom) {
  
  # Check column names
  print("Columns in gom_NAA:")
  print(names(gom_NAA))
  
  print("Columns in waa_gom:")
  print(names(waa_gom))
  
  # Select relevant columns from gom_NAA
  SSB_stock <- gom_NAA[, c("YEAR", "Age4", "Age5", "Age6", "Age7", "MEAN")]
  
  # Debug: Check initial SSB_stock structure
  print("Initial SSB_stock:")
  print(head(SSB_stock))
  
  # Merge with waa_gom
  SSB_stock <- merge(SSB_stock, waa_gom, by = "YEAR", all = TRUE)
  
  # Debug: Check merged data
  print("Merged SSB_stock:")
  print(head(SSB_stock))
  print(names(SSB_stock))  # Print column names
  
  # Check for NAs
  if (any(is.na(SSB_stock))) {
    print("NAs present in SSB_stock:")
    print(SSB_stock[is.na(SSB_stock), ])
  }
  
  # Ensure correct columns are present
  age_columns <- c("Age4.x", "Age5.x", "Age6.x", "Age7.x", "MEAN.x")
  if (!all(age_columns %in% names(SSB_stock))) {
    stop(paste("Missing columns in SSB_stock:", paste(age_columns[!age_columns %in% names(SSB_stock)], collapse = ", ")))
  }
  
  # Calculate SSB using the correct column names
  SSB_stock$SSB <- (
    (SSB_stock$Age4.x * SSB_stock$MEAN.x * (SSB_stock$Age4.x > 0)) +
      (SSB_stock$Age5.x * SSB_stock$MEAN.x * (SSB_stock$Age5.x > 0)) +
      (SSB_stock$Age6.x * SSB_stock$MEAN.x * (SSB_stock$Age6.x > 0)) +
      (SSB_stock$Age7.x * SSB_stock$MEAN.x * (SSB_stock$Age7.x > 0))
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
  SSB_gom <- aggregate(SSB ~ YEAR, data = SSB_stock, FUN = sum)
  
  return(SSB_gom)
}

########

#Run the function and store the result
ssb_gom <- get_ssb(gom_NAA, waa_gom)

#Print the results
print(ssb_gom)

#Check the structure of the resulting data frame
str(ssb_gom)

#Check for any NA values
if (any(is.na(ssb_gom))) {
  print("There are NA values in the results.")
} else {
  print("No NA values in the results.")
}

#Print the first few rows of the result
print(head(ssb_gom))

#Get summary statistics
print(summary(ssb_gom))

#################################
ggplot(data = ssb_gom, aes(x = YEAR, y = SSB)) +
  geom_line(size = 1) +         # Line for SSB
  geom_point(color = "red", size = 2) +         # Points for each year
  labs(title = "Gulf of Maine Spawning Stock Biomass (SSB)",
       x = "Year",
       y = "SSB") +
  theme_minimal() +                             # Minimal theme for a clean look
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Center and size title
    axis.title = element_text(size = 12),               # Size of axis titles
    axis.text = element_text(size = 10)                 # Size of axis text
  )










