library(ggplot2)
library(corrplot)

#must have all the dfs in environment ready to go for this to work
print(age1df) 
print (rssbdf)
print(conddf)
print(waadf)


#----------------- AGE 1 RECRUITMENT ------------------
#Fall Plot
with(age1df, {
  all_values <- c(gom_fall_age1_nmfs, snema_fall_age1_nmfs, gbk_fall_age1_nmfs)
  all_values <- all_values[is.finite(all_values)]
  
  col_gom   <- "#1f77b4"  
  col_snema <- "#d62728"  
  col_gbk   <- "#2ca02c"  
  

  plot(year, gom_fall_age1_nmfs, type = "l", col = col_gom, lwd = 2,
       ylim = range(all_values), ylab = "NAA 1", xlab = "Years", main = "Fall Age 1 Recruitment")
  
  lines(year, snema_fall_age1_nmfs, col = col_snema, lwd = 2)
  lines(year, gbk_fall_age1_nmfs, col = col_gbk, lwd = 2)
})

legend("topleft",
       legend = c("GOM", "SNEMA", "GBK"),
       col = c("#1f77b4", "#d62728", "#2ca02c"),
       lty = 1,
       lwd = 2)

#Spring Plot
with(age1df, {
  all_values <- c(gom_spring_age1_nmfs, snema_spring_age1_nmfs, gbk_spring_age1_nmfs)
  all_values <- all_values[is.finite(all_values)]
  
  col_gom   <- "#1f77b4"  
  col_snema <- "#d62728"  
  col_gbk   <- "#2ca02c"  
  
  
  plot(year, gom_spring_age1_nmfs, type = "l", col = col_gom, lwd = 2,
       ylim = range(all_values), ylab = "NAA 1", xlab = "Years", main = "Spring Age 1 Recruitment")
  
  lines(year, snema_spring_age1_nmfs, col = col_snema, lwd = 2)
  lines(year, gbk_spring_age1_nmfs, col = col_gbk, lwd = 2)
})

legend("topleft",
       legend = c("GOM", "SNEMA", "GBK"),
       col = c("#1f77b4", "#d62728", "#2ca02c"),
       lty = 1,
       lwd = 2)

#let's add in madmf indices to this corrplot too
snema_NMFS_age1 <- snema_age1[which(snema_age1$SURVEY=='NMFS spring BTS'),]
snema_NMFS_age1 <- snema_NMFS_age1 %>% dplyr::select(YEAR, NO_AT_AGE)
names(snema_NMFS_age1)[names(snema_NMFS_age1) == 'YEAR'] <- 'Year'

snema_MADMF_age1 <- snema_age1[which(snema_age1$SURVEY=='MADMF spring BTS'),]
snema_MADMF_age1 <- snema_MADMF_age1 %>% dplyr::select(YEAR, NO_AT_AGE)
names(snema_MADMF_age1)[names(snema_MADMF_age1) == 'YEAR'] <- 'Year'
snema_MADMF_age1 <- na.omit(snema_MADMF_age1)

# Rename NO_AT_AGE columns to differentiate after merge
snema_NMFS_age1 <- snema_NMFS_age1 %>%
  rename(NO_AT_AGE1_NMFS = NO_AT_AGE)

snema_MADMF_age1 <- snema_MADMF_age1 %>%
  rename(NO_AT_AGE1_MADMF = NO_AT_AGE)

# Merge datasets by Year (inner join to keep only years in both)
merged_snema_age1 <- inner_join(snema_NMFS_age1, snema_MADMF_age1, by = "Year")

print(merged_snema_age1)

#All of age 1 correlation test with Pearson correlation coefficient
age1mergedcor <- cor(merged_snema_age1[, c("NO_AT_AGE1_NMFS", "NO_AT_AGE1_MADMF")],
               use = "complete.obs")  # handles any NA values

corrplot(age1mergedcor, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, # text label color and angle
         addCoef.col = "black",        # show correlation values
         number.cex = 0.8,             # size of numbers
         col = colorRampPalette(c("red", "white", "blue"))(200))  # color scale

#0.06 pearson correlation between both surveys, so some similarity but mostly different

#All of age 1 correlation test with Pearson correlation coefficient
age1cor <- cor(age1df[, c("gom_fall_age1_nmfs", "snema_fall_age1_nmfs", "gbk_fall_age1_nmfs", "gom_spring_age1_nmfs", "snema_spring_age1_nmfs", "gbk_spring_age1_nmfs")],
    use = "complete.obs")  # handles any NA values

corrplot(age1cor, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, # text label color and angle
         addCoef.col = "black",        # show correlation values
         number.cex = 0.8,             # size of numbers
         col = colorRampPalette(c("red", "white", "blue"))(200))  # color scale

#----------------- RECRUITS PER SPAWNER ------------------
#Fall Plot
with(rssbdf, {
  all_values <- c(GOM_Fall_rssb, SNEMA_Fall_rssb, GBK_Fall_rssb)
  all_values <- all_values[is.finite(all_values)]
  
  col_gom   <- "#1f77b4"  
  col_snema <- "#d62728"  
  col_gbk   <- "#2ca02c"  
  
  
  plot(Year, GOM_Fall_rssb, type = "l", col = col_gom, lwd = 2,
       ylim = range(all_values), ylab = "R/SSB", xlab = "Years", main = "Fall Recruits per Spawner")
  
  lines(Year, SNEMA_Fall_rssb, col = col_snema, lwd = 2)
  lines(Year, GBK_Fall_rssb, col = col_gbk, lwd = 2)
})

legend("topleft",
       legend = c("GOM", "SNEMA", "GBK"),
       col = c("#1f77b4", "#d62728", "#2ca02c"),
       lty = 1,
       lwd = 2)

#Spring Plot
with(rssbdf, {
  all_values <- c(GOM_Spring_rssb, SNEMA_Spring_rssb, GBK_Spring_rssb)
  all_values <- all_values[is.finite(all_values)]
  
  col_gom   <- "#1f77b4"  
  col_snema <- "#d62728"  
  col_gbk   <- "#2ca02c"  
  
  
  plot(Year, GOM_Spring_rssb, type = "l", col = col_gom, lwd = 2,
       ylim = range(all_values), ylab = "R/SSB", xlab = "Years", main = "Spring Recruits per Spawner")
  
  lines(Year, SNEMA_Spring_rssb, col = col_snema, lwd = 2)
  lines(Year, GBK_Spring_rssb, col = col_gbk, lwd = 2)
})

legend("topleft",
       legend = c("GOM", "SNEMA", "GBK"),
       col = c("#1f77b4", "#d62728", "#2ca02c"),
       lty = 1,
       lwd = 2)

#All of R/SSB correlation test with Pearson correlation coefficient
rssbcor <- cor(rssbdf[, c("GOM_Fall_rssb", "SNEMA_Fall_rssb", "GBK_Fall_rssb", "GOM_Spring_rssb", "SNEMA_Spring_rssb", "GBK_Spring_rssb")],
               use = "complete.obs")  # handles any NA values

corrplot(rssbcor, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, # text label color and angle
         addCoef.col = "black",        # show correlation values
         number.cex = 0.8,             # size of numbers
         col = colorRampPalette(c("red", "white", "blue"))(200))  # color scale

#----------------- CONDITION ------------------
#Plot
with(conddf, {
  all_values <- c(condition_GOM, condition_MAB, condition_GBK)
  all_values <- all_values[is.finite(all_values)]
  
  col_gom   <- "#1f77b4"  
  col_mab <- "#d62728"  
  col_gbk   <- "#2ca02c"  
  
  
  plot(year, condition_GOM, type = "l", col = col_gom, lwd = 2,
       ylim = range(all_values), ylab = "Kn=W/W′", xlab = "Years", main = "Winter Flounder Relative Condition")
  
  lines(year, condition_MAB, col = col_mab, lwd = 2)
  lines(year, condition_GBK, col = col_gbk, lwd = 2)
})

legend("topleft",
       legend = c("GOM", "MAB", "GBK"),
       col = c("#1f77b4", "#d62728", "#2ca02c"),
       lty = 1,
       lwd = 2)

#Correlation test with Pearson correlation coefficient
condcor <- cor(conddf[, c("condition_GOM", "condition_MAB", "condition_GBK")],
               use = "complete.obs")  # handles any NA values

corrplot(condcor, method = "color", type = "upper",
         tl.col = "black", tl.srt = 45, # text label color and angle
         addCoef.col = "black",        # show correlation values
         number.cex = 0.8,             # size of numbers
         col = colorRampPalette(c("red", "white", "blue"))(200))  # color scale

#----------------- WAA 1 ------------------
 
  #Correlation test with Pearson correlation coefficient
  waacor <- cor(waadf[, c("age1GBK", "age2GBK", "age3GBK", "age4GBK", "age5GBK", "age6GBK", "age7GBK", 
                          "age1SNEMA", "age2SNEMA", "age3SNEMA", "age4SNEMA", "age5SNEMA", "age6SNEMA", "age7SNEMA",
                          "age1GOM", "age2GOM", "age3GOM", "age4GOM", "age5GOM", "age6GOM", "age7GOM")],
                 use = "complete.obs")  # handles any NA values
  
corrplot(waacor, method = "color", type = "upper", 
         tl.col = "black", tl.srt = 45, 
         addCoef.col = "black",  
         number.cex = 0.8, 
         col = colorRampPalette(c("red", "white", "blue"))(200)) 















