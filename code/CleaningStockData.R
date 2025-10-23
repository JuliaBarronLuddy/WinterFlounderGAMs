#AGE 1 PROCESSING FOR CLEAN DATA
# Load required packages
library(dplyr)
library(readr)
library(here)

# Load raw age datasets
gom_age <- read_csv(here("data/gom_age.csv"))
gbk_age <- read_csv(here("data/gbk_age.csv"))
snema_age <- read_csv(here("data/snema_age.csv"))

# Function to process each stock's raw data
process_age_data <- function(df, stock_name) {
  df %>%
    filter(AGE == 1) %>%                                  # Keep only age-1
    dplyr::select(Year = YEAR, SURVEY, NO_AT_AGE) %>%            # Select relevant columns
    mutate(stock = stock_name) %>%                        # Add stock identifier
    filter(!is.na(NO_AT_AGE))                             # Drop rows with missing values
}

# Apply function to each stock
gom_age1   <- process_age_data(gom_age, "GOM")
gbk_age1   <- process_age_data(gbk_age, "GB")
snema_age1 <- process_age_data(snema_age, "SNEMA")

# Combine all into one age1df
age1df <- bind_rows(gom_age1, gbk_age1, snema_age1)

# Peek at the result
head(age1df)


#WAA CLEANING
# Clean GOM data: remove MEAN column and add stock label
waa_gom_clean <- waa_gom %>%
  dplyr::select(-MEAN) %>%
  mutate(stock = "GOM")

# Clean GBK data: add stock label
waa_gbk_clean <- waa_gbk %>%
  mutate(stock = "GBK")

# Clean SNEMA data: add stock label
waa_snema_clean <- waa_snema %>%
  mutate(stock = "SNEMA")

# Combine all three
waadf <- bind_rows(waa_gom_clean, waa_gbk_clean, waa_snema_clean) %>%
  relocate(stock, .after = Year) %>%  # Optional: move stock next to Year
  arrange(stock, Year)


#COND CLEANING
conddf <- conddf %>%
  rename(Year = year) %>%
  pivot_longer(
    cols = starts_with("condition_"),
    names_to = "stock",
    values_to = "Condition"
  ) %>%
  mutate(
    stock = case_when(
      stock == "condition_GOM" ~ "GOM",
      stock == "condition_GBK" ~ "GBK",
      stock == "condition_MAB" ~ "SNEMA",  # Rename proxy
      TRUE ~ stock
    )
  ) %>%
  arrange(stock, Year)

#COG CLEANING
cogdf <- cogdf %>%
  rename(Year = year) %>%
  dplyr::select(
    Year,
    fall_cog_lat, fall_cog_depth,
    spring_cog_lat, spring_cog_depth
  ) %>%
  arrange(Year)

#RSSB CLEANING
library(dplyr)
library(tidyr)
library(stringr)

# Reshape RSSB data from wide to long
rssb_clean <- rssbdf %>%
  pivot_longer(
    cols = -Year,
    names_to = "Stock_Season",
    values_to = "rssb"
  ) %>%
  # Separate column into stock and season
  separate(Stock_Season, into = c("stock", "season", "metric"), sep = "_") %>%
  filter(metric == "rssb") %>%  # Only keep rssb metric
  dplyr::select(-metric) %>%
  pivot_wider(
    names_from = season,
    values_from = rssb,
    names_prefix = "rssb_"
  ) %>%
  mutate(stock = recode(stock,
                        "SNEMA" = "SNEMA",
                        "GBK" = "GBK",
                        "GOM" = "GOM")) %>%
  dplyr::select(Year, stock, rssb_Fall = rssb_Fall, rssb_Spring = rssb_Spring) %>%
  arrange(stock, Year)


