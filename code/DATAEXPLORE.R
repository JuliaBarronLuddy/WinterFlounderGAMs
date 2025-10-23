library(ecodata)
library(dplyr)
library(here)
library(tidyverse)
library(mgcv)

##GLORYS Bottom Temp in shape of WFL stock areas
#Data load
btdata <- read.csv(here("data/GLORYS_monthly_BottomT_winter_flounder_1993_2024.csv"))
btdata <- btdata %>%
  pivot_wider(names_from = statistic, values_from = value) 
btdata <- btdata %>%
  mutate(date = ymd(paste(year, month, "01", sep ="-")))
#make mean yearly values into a sep. dataset
bt_yearly <- btdata %>%
  filter(var.name == "BottomT") %>%
  group_by(year, stock) %>%
  summarize(mean_yearly_temp = mean(mean, na.rm = TRUE), .groups = "drop")

#plotting mean yearly values over time by stock
ggplot(bt_yearly, aes(x = year, y = mean_yearly_temp, color = stock)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Mean Yearly Bottom Temperature by Stock",
    x = "Year",
    y = "Mean Bottom Temperature (°C)",
    color = "Stock"
  ) +
  theme_minimal()

#linear model
lm_bt <- lm(mean_yearly_temp ~ year * stock, data = bt_yearly)
shapiro.test(resid(lm_bt)) #p value is 0.4463 which means data is normally distributed
summary(lm_bt)
#summary shows not a significant increase over time even though all stock bt is increasing

#plot with lm trend lines
ggplot(bt_yearly, aes(x = year, y = mean_yearly_temp, color = stock)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    title = "Yearly Mean Bottom Temperature with Linear Trends",
    x = "Year",
    y = "Mean Bottom Temp (°C)"
  ) +
  theme_minimal()

