library(dplyr)
library(tidyr)
library(lubridate)
library(here)

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
### Environmental Covariates
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------

### ---------------------------
### 1. Load and Prep Datasets
### ---------------------------

# Stocks
stocks <- c("GOM", "GB", "SNEMA")

# Create survey table for each stock
survey_years <- 1980:2023
survey_table <- expand.grid(Year = survey_years, Season = c("Spring", "Fall"), stock = stocks) %>%
  mutate(Survey_Start = case_when(
    Season == "Spring" ~ ymd(paste(Year, "02", "01", sep = "-")),
    Season == "Fall"   ~ ymd(paste(Year, "09", "01", sep = "-"))
  ))

# --- WCR
wcrdata <- read.csv(here("data/wcr_census.csv")) %>%
  mutate(DOB = as.Date(DOB),
         year = year(DOB),
         month = month(DOB)) %>%
  group_by(year, month) %>%
  summarise(wcr_count = n(), .groups = 'drop') %>%
  mutate(date = ymd(paste(year, month, "01", sep = "-"))) %>%
  complete(date = seq(min(date), max(date), by = "month")) %>%
  mutate(wcr_count = replace_na(wcr_count, 0))

# --- Bottom Temperature (by stock)
btdata <- read.csv(here("data/GLORYS_monthly_BottomT_winter_flounder_1993_2024.csv")) %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  filter(var.name == "BottomT") %>%
  mutate(
    date = ymd(paste(year, month, "01", sep = "-")),
    bt_mean = mean
  ) %>%
  dplyr::select(date, stock, bt_mean)

# --- Surface Temperature (by stock)
sstdata <- read.csv(here("data/monthly_surfaceT_by_stock.csv")) %>%
  pivot_wider(names_from = statistic, values_from = value) %>%
  filter(var.name == "SurfaceT") %>%
  mutate(
    date = ymd(paste(year, month, "01", sep = "-")),
    sst_mean = mean
  ) %>%
  dplyr::select(date, stock, sst_mean) %>%
  mutate(stock = recode(stock,
                              "GBK" = "GB"))

# --- GSI
gsi <- ecodata::gsi %>%
  filter(Var == "gulf stream index") %>%
  separate(Time, into = c("Year", "Month"), sep = "\\.", convert = TRUE) %>%
  mutate(date = ymd(paste(Year, Month, "01", sep = "-"))) %>%
  dplyr::select(date, gsi = Value)

# --- AMO
AMO_data <- read.table(here("data/amo.txt"), fileEncoding = "UTF-8")
colnames(AMO_data) <- c("Year", as.character(1:12))
amo_long <- pivot_longer(AMO_data, cols = -Year, names_to = "month", values_to = "amo") %>%
  mutate(month = as.integer(month),
         date = ymd(paste(Year, month, "01", sep = "-"))) %>%
  dplyr::select(date, amo)

# --- NAO
NAO_data <- read.table(here("data/nao.txt"), fileEncoding = "UTF-8")
colnames(NAO_data) <- c("Year", as.character(1:12))
nao_long <- pivot_longer(NAO_data, cols = -Year, names_to = "month", values_to = "nao") %>%
  mutate(month = as.integer(month),
         date = ymd(paste(Year, month, "01", sep = "-"))) %>%
  dplyr::select(date, nao)

# --- SSB (by stock, yearly-will be lagged) – 
ssb_gbk <- ssb_gbk %>%
  mutate(stock = "GB")

ssb_gom <- ssb_gom %>%
  mutate(stock = "GOM")

ssb_snema <- ssb_snema %>%
  mutate(stock = "SNEMA") %>%
  rename(year = Year)

ssb_by_stock <- bind_rows(ssb_gbk, ssb_gom, ssb_snema)
ssb_by_stock <- ssb_by_stock %>%
  dplyr::select(year, stock, SSB) %>%
  rename(Year = year)
  
write.csv(ssb_by_stock, here("data/ssb_by_stock.csv"), row.names = FALSE)


### ---------------------------
### 2. Lag Function
### ---------------------------

get_lagged_mean_rowwise <- function(env_df, months_back, var_name, stock_col = NULL) {
  this_row <- cur_data()
  start_date <- as.Date(this_row$Survey_Start)
  lag_start <- start_date %m-% months(months_back - 1)
  lag_end   <- start_date
  
  df <- env_df
  if (!is.null(stock_col) && "stock" %in% names(env_df)) {
    df <- df %>% filter(stock == this_row$stock)
  }
  
  df %>%
    filter(date >= lag_start & date <= lag_end) %>%
    summarise(mean_val = mean(.data[[var_name]], na.rm = TRUE)) %>%
    pull(mean_val)
}

### ---------------------------
### 3. Compute All Lagged Variables
### ---------------------------

lagged_env <- survey_table %>%
  rowwise() %>%
  mutate(
    # WCR – not stock-specific
    WCR_lag_12mo = get_lagged_mean_rowwise(wcrdata, 12, "wcr_count"),
    WCR_lag_6mo  = get_lagged_mean_rowwise(wcrdata, 6,  "wcr_count"),
    
    # BT & SST – stock-specific
    BT_lag_12mo  = get_lagged_mean_rowwise(btdata, 12, "bt_mean", stock_col = "stock"),
    BT_lag_6mo   = get_lagged_mean_rowwise(btdata, 6,  "bt_mean", stock_col = "stock"),
    SST_lag_12mo = get_lagged_mean_rowwise(sstdata, 12, "sst_mean", stock_col = "stock"),
    SST_lag_6mo  = get_lagged_mean_rowwise(sstdata, 6,  "sst_mean", stock_col = "stock"),
    
    # All others – not stock-specific
    GSI_lag_12mo = get_lagged_mean_rowwise(gsi, 12, "gsi"),
    GSI_lag_6mo  = get_lagged_mean_rowwise(gsi, 6,  "gsi"),
    AMO_lag_12mo = get_lagged_mean_rowwise(amo_long, 12, "amo"),
    AMO_lag_6mo  = get_lagged_mean_rowwise(amo_long, 6,  "amo"),
    NAO_lag_12mo = get_lagged_mean_rowwise(nao_long, 12, "nao"),
    NAO_lag_6mo  = get_lagged_mean_rowwise(nao_long, 6,  "nao")
  ) %>%
  ungroup()


### ---------------------------
### 4. Merge SSB 
### ---------------------------

final_env_data <- lagged_env %>%
  left_join(ssb_by_stock, by = c("Year", "stock"))

### ---------------------------
### 5. Ready for GAMs
### ---------------------------

head(final_env_data)
# Filter final dataset to start from 1995 to ensure full lag coverage
final_env_data<- final_env_data %>%
  filter(Year >= 1995)

###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
### Stock Dynamic Indices
###-----------------------------------------------------------------------------
###-----------------------------------------------------------------------------
#make sure you run through cleaningstockdata.r before doing this

library(dplyr)
library(tidyr)

# Process GOM

age1_wide_gom <- age1df %>%
  filter(stock == "GOM") %>%
  dplyr::select(Year, Survey = SURVEY, Age1 = NO_AT_AGE) %>%
  arrange(Year) %>%
  pivot_wider(names_from = Survey, values_from = Age1, names_prefix = "Age1_") %>%
  rename_with(~ gsub(" ", "_", .x)) %>%
  mutate(Year = as.integer(Year))

rssb_stock_gom <- rssb_clean %>%
  filter(stock == "GOM") %>%
  dplyr::select(Year, rssb_Fall, rssb_Spring)

waa_stock_gom <- waadf %>%
  filter(stock == "GOM") %>%
  mutate(Year = as.integer(Year)) %>%
  arrange(Year)

cond_stock_gom <- conddf %>%
  filter(stock == "GOM") %>%
  dplyr::select(Year, Condition) %>%
  mutate(Year = as.integer(Year)) %>%
  arrange(Year)

cog_stock_gom <- cogdf %>%
  mutate(Year = as.integer(Year)) %>%
  dplyr::select(Year, fall_cog_lat, fall_cog_depth, spring_cog_lat, spring_cog_depth)

stock_indices_gom <- age1_wide_gom %>%
  full_join(rssb_stock_gom, by = "Year") %>%
  full_join(waa_stock_gom, by = "Year") %>%
  full_join(cond_stock_gom, by = "Year") %>%
  full_join(cog_stock_gom, by = "Year") %>%
  mutate(stock = "GOM") %>%
  arrange(Year)

# Process GBK


age1_wide_gbk <- age1df %>%
  filter(stock == "GB") %>%
  dplyr::select(Year, Survey = SURVEY, Age1 = NO_AT_AGE) %>%
  arrange(Year) %>%
  pivot_wider(names_from = Survey, values_from = Age1, names_prefix = "Age1_") %>%
  rename_with(~ gsub(" ", "_", .x)) %>%
  mutate(Year = as.integer(Year))

rssb_stock_gbk <- rssb_clean %>%
  filter(stock == "GBK") %>%
  dplyr::select(Year, rssb_Fall, rssb_Spring)

waa_stock_gbk <- waadf %>%
  filter(stock == "GBK") %>%
  mutate(Year = as.integer(Year)) %>%
  arrange(Year)

cond_stock_gbk <- conddf %>%
  filter(stock == "GBK") %>%
  dplyr::select(Year, Condition) %>%
  mutate(Year = as.integer(Year)) %>%
  arrange(Year)

cog_stock_gbk <- cogdf %>%
  mutate(Year = as.integer(Year)) %>%
  dplyr::select(Year, fall_cog_lat, fall_cog_depth, spring_cog_lat, spring_cog_depth)

stock_indices_gbk <- age1_wide_gbk %>%
  full_join(rssb_stock_gbk, by = "Year") %>%
  full_join(waa_stock_gbk, by = "Year") %>%
  full_join(cond_stock_gbk, by = "Year") %>%
  full_join(cog_stock_gbk, by = "Year") %>%
  mutate(stock = "GBK") %>%
  arrange(Year)

# Process SNEMA


age1_wide_snema <- age1df %>%
  filter(stock == "SNEMA") %>%
  dplyr::select(Year, Survey = SURVEY, Age1 = NO_AT_AGE) %>%
  arrange(Year) %>%
  pivot_wider(names_from = Survey, values_from = Age1, names_prefix = "Age1_") %>%
  rename_with(~ gsub(" ", "_", .x)) %>%
  mutate(Year = as.integer(Year))

rssb_stock_snema <- rssb_clean %>%
  filter(stock == "SNEMA") %>%
  dplyr::select(Year, rssb_Fall, rssb_Spring)

waa_stock_snema <- waadf %>%
  filter(stock == "SNEMA") %>%
  mutate(Year = as.integer(Year)) %>%
  arrange(Year)

cond_stock_snema <- conddf %>%
  filter(stock == "SNEMA") %>%
  dplyr::select(Year, Condition) %>%
  mutate(Year = as.integer(Year)) %>%
  arrange(Year)

cog_stock_snema <- cogdf %>%
  mutate(Year = as.integer(Year)) %>%
  dplyr::select(Year, fall_cog_lat, fall_cog_depth, spring_cog_lat, spring_cog_depth)

stock_indices_snema <- age1_wide_snema %>%
  full_join(rssb_stock_snema, by = "Year") %>%
  full_join(waa_stock_snema, by = "Year") %>%
  full_join(cond_stock_snema, by = "Year") %>%
  full_join(cog_stock_snema, by = "Year") %>%
  mutate(stock = "SNEMA") %>%
  arrange(Year)

# Combine all stocks into one dataset
combined_stock_indices <- bind_rows(stock_indices_gom, stock_indices_gbk, stock_indices_snema)

# View the first few rows
print(head(combined_stock_indices))
final_stock_indices<- combined_stock_indices %>%
  filter(Year >= 1995)
