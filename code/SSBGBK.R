library(ecodata)
library(dplyr)
library(here)
library(tidyverse)
library(mgcv)
library(ggplot2)
################################################################################

#GBK age data wrangling
gbk_age <- read.csv(here("data/gbk_age.csv"))
gbk_NAA <- gbk_age[which(gbk_age$AGE > 3),] #make ages 4 plus
gbk_NAA <- gbk_NAA[!gbk_NAA$YEAR  < 1982,] #make years 1982 plus to match WAA data
gbk_NAA <- subset(gbk_NAA, select = c(SURVEY, SEASON, YEAR, AGE, NO_AT_AGE)) #select the columns we want to use
gbk_NAA <- gbk_NAA %>%
  mutate(Age_Group = ifelse(AGE >= 7, 7, AGE)) %>% #Create a new column to combine ages 7+
  group_by(YEAR, Age_Group) %>%                #Group by Year and Age_Group
  dplyr::summarise(MEAN = mean(NO_AT_AGE, na.rm = TRUE), .groups = 'drop') # Compute the average
#I know in my heart that age 7 represents 7+
#pivot so data can match columns in WAA data
gbk_NAA <- gbk_NAA %>%
  pivot_wider(
    names_from = Age_Group,
    values_from = MEAN
  )
names(gbk_NAA)[names(gbk_NAA) == '4'] <- 'Age4' #match column names of AGE to WAA data
names(gbk_NAA)[names(gbk_NAA) == '5'] <- 'Age5' 
names(gbk_NAA)[names(gbk_NAA) == '6'] <- 'Age6' 
names(gbk_NAA)[names(gbk_NAA) == '7'] <- 'Age7' #again, I know that this is 7+, not 7
gbk_NAA <- gbk_NAA %>%
  mutate(MEAN = rowMeans(across(starts_with("Age")), na.rm = TRUE)) #find mean of all ages and include as separate column


#WAA gbk data wrangling
waa_gbk <- read.csv(here("data/WAA_gbk.csv"))
names(waa_gbk)[names(waa_gbk) == 'Year'] <- 'YEAR' #match case of YEAR to NAA data
waa_gbk <- subset(waa_gbk, select = c(YEAR, Age4, Age5, Age6, Age7)) #get rid of premature fish ages (<4 years)
waa_gbk <- waa_gbk %>%
  mutate(MEAN = rowMeans(across(starts_with("Age")), na.rm = TRUE)) #find mean of all ages and include as separate column
waa_gbk <- as_tibble(waa_gbk)


#now it's time to get SSB, using Jamie's function
#et_ssb<- function(gbk_NAA, waa_gbk){
  #select relevant columns from gbk_NAA
  #SSB_stock <- gbk_NAA[,c(1:5)]
  #names(SSB_stock)<-c("YEAR", "Age4","Age5","Age6","Age7")
  #merge with waa_gbk
  #SSB_stock<-merge(SSB_stock,waa_gbk,by=c("YEAR"),all=TRUE)
  #calculate ssb
  #SSB_stock$SSB <- (SSB_stock$Age4 * SSB_stock$MEAN * (SSB_stock$Age4 == 4))+
                    #(SSB_stock$Age5 * SSB_stock$MEAN * (SSB_stock$Age5 == 5))+
                    #(SSB_stock$Age6 * SSB_stock$MEAN * (SSB_stock$Age6 == 6))+
                    #(SSB_stock$Age7 * SSB_stock$MEAN * (SSB_stock$Age7 == 7))+ 
  #keep only relevant columns
  #SSB_stock<-SSB_stock[,c("YEAR", "SSB")]
  #aggregate SSB by year
  #SSB_aggregated<-aggregate(SSB~YEAR,SSB_stock,FUN=sum)
  #return(SSB_aggregated)
#}

##############################################################
#Chatgpt test
get_ssb <- function(gbk_NAA, waa_gbk) {
  
  # Check column names
  print("Columns in gbk_NAA:")
  print(names(gbk_NAA))
  
  print("Columns in waa_gbk:")
  print(names(waa_gbk))
  
  # Select relevant columns from gbk_NAA
  SSB_stock <- gbk_NAA[, c("YEAR", "Age4", "Age5", "Age6", "Age7", "MEAN")]
  
  # Debug: Check initial SSB_stock structure
  print("Initial SSB_stock:")
  print(head(SSB_stock))
  
  # Merge with waa_gbk
  SSB_stock <- merge(SSB_stock, waa_gbk, by = "YEAR", all = TRUE)
  
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
  SSB_aggregated <- aggregate(SSB ~ YEAR, data = SSB_stock, FUN = sum)
  
  return(SSB_aggregated)
}
######################################################################



#Run the function and store the result
ssb_gbk <- get_ssb(gbk_NAA, waa_gbk)

#Print the results
print(ssb_gbk)

#Check the structure of the resulting data frame
str(ssb_gbk)

#Check for any NA values
if (any(is.na(ssb_gbk))) {
  print("There are NA values in the results.")
} else {
  print("No NA values in the results.")
}

#Print the first few rows of the result
print(head(ssb_gbk))

#Get summary statistics
print(summary(ssb_gbk))

#################################
ggplot(data = ssb_gbk, aes(x = YEAR, y = SSB)) +
  geom_line(size = 1) +         # Line for SSB
  geom_point(color = "red", size = 2) +         # Points for each year
  labs(title = "Georges Bank Spawning Stock Biomass (SSB)",
       x = "Year",
       y = "SSB") +
  theme_minimal() +                             # Minimal theme for a clean look
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Center and size title
    axis.title = element_text(size = 12),               # Size of axis titles
    axis.text = element_text(size = 10)                 # Size of axis text
  )





























