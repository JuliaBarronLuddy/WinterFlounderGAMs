library(tidyverse)
library(dplyr)
library(corrplot)
library(here)
library(reshape2)
library(ggplot2)
library(ggcorrplot)
library(vegan)
library(pheatmap)

#SPEARMAN RANK CORRELATION OF ALL CHANGEPOINTS

EnvYears <- data.frame(
  GBK = c("", "", "1993", "1993, 2014", "2009", "2010", "2014", "2015", "", "", "2001, 2012"),
  GOM = c("", "", "2003", "", "1985, 2001", "1985, 2012", "1983", "1983", "", "", "1999"),
  SNEMA = c("1983, 1992", "1983", "2015", "", "1998", "1993, 1999", "1993", "1991, 2012", "1991, 2012", "1991", ""),
  row.names = c("R_EnvYearFall", "R_EnvYearSpring", "RSSB_EnvYearFall", "RSSB_EnvYearSpring",
                "WAA_EnvYear1", "WAA_EnvYear2", "WAA_EnvYear3", "WAA_EnvYear4", "WAA_EnvYear5", 
                "WAA_EnvYear6", "COND_EnvYear"),
  stringsAsFactors = FALSE
)

# Convert to long format and extract changepoint years
env_years_long <- EnvYears %>%
  rownames_to_column(var = "EnvVar") %>%
  pivot_longer(cols = GBK:SNEMA, names_to = "Stock", values_to = "Years") %>%
  mutate(YearList = str_split(Years, ",\\s*")) %>%
  unnest(YearList) %>%
  filter(YearList != "") %>%
  mutate(Year = as.integer(YearList))

# Extract dynamic type from EnvVar
env_years_long <- env_years_long %>%
  mutate(Dynamic = case_when(
    str_starts(EnvVar, "RSSB") ~ "RSSB",
    str_starts(EnvVar, "R")    ~ "R",
    str_starts(EnvVar, "WAA")  ~ "WAA",
    str_starts(EnvVar, "COND") ~ "COND",
    TRUE ~ "Other"
  )) %>%
  mutate(StockDynamic = paste(Stock, Dynamic, sep = "_"))


presence_matrix_sd <- env_years_long %>%
  dplyr::select(StockDynamic, Year) %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = Year,
    values_from = present,
    values_fill = list(present = 0)  
  ) %>%
  column_to_rownames("StockDynamic")










# View the matrix
head(presence_matrix_sd)

# Spearman correlation
spearman_sd <- cor(t(presence_matrix_sd), method = "spearman")

# Round and view
round(spearman_sd, 2)














env_years_long

#Create a unique list of Stock_Year
env_years_long <- env_years_long %>%
  mutate(Stock_Year = paste0(Stock, "_", YearList))

#Create presence/absence matrix: EnvYear × Stock_Year
presence_matrix <- env_years_long %>%
  dplyr::select(EnvYear, Stock_Year) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Stock_Year, values_from = present, values_fill = 0) %>%
  column_to_rownames(var = "EnvYear")

presence_matrix

# Convert to binary presence (0/1)
presence_binary <- presence_matrix
presence_binary[presence_binary > 0] <- 1

# Assuming presence_matrix rows are named with EnvYear, e.g., "WAA_EnvYear1"
env_year_types <- data.frame(
  EnvYear = rownames(presence_binary)
) %>%
  mutate(Type = case_when(
    str_detect(EnvYear, "^WAA") ~ "WAA",
    str_detect(EnvYear, "^RSSB") ~ "RSSB",
    str_detect(EnvYear, "^R") ~ "R",
    str_detect(EnvYear, "^COND") ~ "COND",
    TRUE ~ "Other"
  ))


# Add Type column
presence_df <- presence_binary %>%
  as.data.frame() %>%
  rownames_to_column("EnvYear") %>%
  left_join(env_year_types, by = "EnvYear")

# Summarise by type: any 1s in each Stock_Year for that type
presence_by_type_stockyear <- presence_df %>%
  group_by(Type) %>%
  summarise(across(where(is.numeric), ~ as.integer(any(.x == 1)))) %>%
  column_to_rownames("Type")

# Transpose: each type is a row
type_matrix <- as.matrix(presence_by_type)

png(here("Figures/Raw_data_trends/Changepoint/Heatmap_Type_StockYear.png"),
    width = 2400, height = 1600, units = "px", res = 300)

pheatmap(presence_by_type_stockyear,
         main = "Changepoint Types Across Stocks and Years",
         color = colorRampPalette(c("white", "steelblue", "darkblue"))(100),
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         fontsize_number = 10,
         fontsize = 12)

dev.off()









# Jaccard distance (between types)
jaccard_dist <- vegdist(type_matrix, method = "jaccard", binary = TRUE)

# Convert distance to similarity
jaccard_sim <- 1 - as.matrix(jaccard_dist)

#Inspect
print(round(jaccard_sim, 2))

png(here("Figures/Raw_data_trends/Changepoint/HeatMap.png"), width = 1800, height = 1600, units = "px", res = 300)

pheatmap(jaccard_sim,
         main = "Jaccard Similarity Between Changepoint Types",
         color = colorRampPalette(c("white", "lightblue", "darkblue"))(100),
         display_numbers = TRUE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         fontsize_number = 10)
dev.off()


# Save corrplot version
png(here("Figures/Raw_data_trends/Changepoint/CorrPlot_Jaccard.png"),
    width = 1800, height = 1600, units = "px", res = 300)

corrplot(jaccard_sim,
         method = "color",
         type = "upper",
         order = "hclust",
         addCoef.col = "black",
         tl.col = "black",
         tl.srt = 45,
         cl.pos = "r",
         is.corr = FALSE)  # Important: Jaccard isn't correlation!

dev.off()

# Transpose to get types as rows, stock-years as columns
spearman_matrix <- cor(t(presence_by_type), method = "spearman")

# View the result
round(spearman_matrix, 2)

png(here("Figures/Raw_data_trends/Changepoint/Spearman.png"),
    width = 1800, height = 1600, units = "px", res = 300)
corrplot(spearman_matrix,
         method = "color",
         type = "upper",
         order = "hclust",
         addCoef.col = "black",
         tl.col = "black",
         tl.srt = 45,
         diag = FALSE)
dev.off()


