library(ecodata)
library(dplyr)
library(here)
library(tidyverse)
library(mgcv)
################################################################################
#Georges Bank

#Data load
gbk_NMFS_age1 <- gbk_age1[which(gom_age1$SURVEY=='NMFS spring BTS'),]
gbk_NMFS_age1 <- gbk_NMFS_age1 %>% select(c('YEAR', 'NO_AT_AGE'))
names(gbk_NMFS_age1)[names(gbk_NMFS_age1) == 'YEAR'] <- 'Year'

#Bottom temperature
bottomtemp <- ecodata::bottom_temp_comp
bottomtemp <- bottomtemp[which(bottomtemp$EPU=='GOM'),]
bottomtemp <- bottomtemp[which(bottomtemp$Var=='Spring_Bottom Temp Anomaly'),]
bottomtemp <- bottomtemp[-which(bottomtemp$Source=='PSY'),]
bottomtemp <- subset(bottomtemp, select = -c(EPU, Source, Var, Units))
names(bottomtemp)[names(bottomtemp) == 'Value'] <- 'BTAnom'
names(bottomtemp)[names(bottomtemp) == 'Time'] <- 'Year'

#Cold pool index
coldpool <- ecodata::cold_pool
coldpool <- coldpool[which(coldpool$Var=='cold_pool_index'),]
coldpool <- coldpool[-which(coldpool$source=='PSY'),]
names(coldpool)[names(coldpool) == 'Value'] <- 'CPI'
names(coldpool)[names(coldpool) == 'Time'] <- 'Year'

#GSI
gsi <- ecodata::gsi
gsi <- gsi[which(gsi$Var=='gulf stream index'),]
gsi <- gsi %>% 
  separate(Time, into = c("Year", "Month"), sep = "\\.", convert = TRUE)
gsiyear <- gsi %>%
  group_by(Year) %>%
  summarize(Mean_Value = mean(Value))
names(gsiyear)[names(gsiyear) == 'Mean_Value'] <- 'GSI'

#Warm Core Ring census
warmcore <- read.csv(here("data/wcr_census.csv"))










