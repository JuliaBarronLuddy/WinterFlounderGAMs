library(ecodata)
library(dplyr)
library(here)
library(tidyverse)
library(mgcv)
################################################################################
#Southern New England Mid Atlantic

##Stock data load
#Recruitment
snema_NMFS_age1 <- snema_age1[which(snema_age1$SURVEY=='NMFS spring BTS'),]
snema_NMFS_age1 <- snema_NMFS_age1 %>% dplyr::select(YEAR, NO_AT_AGE)
names(snema_NMFS_age1)[names(snema_NMFS_age1) == 'YEAR'] <- 'Year'

snema_MADMF_age1 <- snema_age1[which(snema_age1$SURVEY=='MADMF spring BTS'),]
snema_MADMF_age1 <- snema_MADMF_age1 %>% dplyr::select(YEAR, NO_AT_AGE)
names(snema_MADMF_age1)[names(snema_MADMF_age1) == 'YEAR'] <- 'Year'
snema_MADMF_age1 <- na.omit(snema_MADMF_age1)

#pull SSB from SSB calculation datafile SSBSNEMA.R (and r/ssb)
ssb_snema <- ssb_snema
print (rssbdf)

#WAA
waa_snema <- read.csv(here("data/WAA_snema.csv"))
colnames(waa_snema)[colnames(waa_snema) == "Age7."] <- "Age7" #getting rid of random period in column
snema_waa1 <- waa_snema[!is.na(waa_snema$Age1), c("Year", "Age1")]
snema_waa2 <- waa_snema[!is.na(waa_snema$Age2), c("Year", "Age2")]
snema_waa3 <- waa_snema[!is.na(waa_snema$Age3), c("Year", "Age3")]
snema_waa4 <- waa_snema[!is.na(waa_snema$Age4), c("Year", "Age4")]
snema_waa5 <- waa_snema[!is.na(waa_snema$Age5), c("Year", "Age5")]
snema_waa6 <- waa_snema[!is.na(waa_snema$Age6), c("Year", "Age6")]
snema_waa7 <- waa_snema[!is.na(waa_snema$Age7.), c("Year", "Age7")]

#CONDITION
condition <- ecodata::condition #we want Var -- Winter flounder
condition <- condition[which(condition$Var=='Winter flounder'),]
colnames(condition)[colnames(condition) == "Time"] <- "Year"
cond_sne <- condition[which(condition$EPU=='MAB'),]

#DISTRIBUTION IS BELOW
##########################

##Ecovars
#BT
btdata <- read.csv(here("data/GLORYS_monthly_BottomT_winter_flounder_1993_2024.csv"))
btdata <- btdata %>%
  pivot_wider(names_from = statistic, values_from = value) 
btdata <- btdata %>%
  mutate(date = ymd(paste(year, month, "01", sep ="-")))
#make mean yearly values into a sep. dataset
bt_yearly_snema <- btdata %>%
  filter(var.name == "BottomT", stock == "SNEMA") %>%
  group_by(year, stock) %>%
  summarize(mean_yearly_temp = mean(mean, na.rm = TRUE), .groups = "drop") %>%
  rename(
    BottomT = mean_yearly_temp,
    Year = year
  ) %>%
  dplyr::select(Year, BottomT)

#SST
sstdata <- read.csv(here("data/monthly_surfaceT_by_stock.csv"))
sstdata <- sstdata %>%
  pivot_wider(names_from = statistic, values_from = value) 
sstdata <- sstdata %>%
  mutate(date = ymd(paste(year, month, "01", sep ="-")))
#make mean yearly values into a sep. dataset
sst_yearly_snema <- sstdata %>%
  filter(var.name == "SurfaceT", stock == "SNEMA") %>%
  group_by(year, stock) %>%
  summarize(mean_yearly_temp = mean(mean, na.rm = TRUE), .groups = "drop") %>%
  rename(
    SurfaceT = mean_yearly_temp,
    Year = year
  ) %>%
  dplyr::select(Year, SurfaceT)

#GSI
gsi <- ecodata::gsi
gsi <- gsi[which(gsi$Var=='gulf stream index'),]
gsi <- gsi %>% 
  separate(Time, into = c("Year", "Month"), sep = "\\.", convert = TRUE)
gsiyear <- gsi %>%
  group_by(Year) %>%
  summarize(Mean_Value = mean(Value))
names(gsiyear)[names(gsiyear) == 'Mean_Value'] <- 'GSI'

#AMO
AMO_data <- read.table(here("data/amo.txt"), fileEncoding = "UTF-8")
colnames(AMO_data) <- c("Year", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
amoyear <- AMO_data %>%
  rowwise() %>%
  mutate(Mean = mean(c_across(2:13), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(Year, Mean)
names(amoyear)[names(amoyear) == 'Mean'] <- 'AMO'
amoyear <- amoyear[amoyear$Year !=2023, ] #got rid of last column because it was not accurate

#NAO
NAO_data <- read.table(here("data/nao.txt"), fileEncoding = "UTF-8")
colnames(NAO_data) <- c("Year", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
naoyear <- NAO_data %>%
  rowwise() %>%
  mutate(Mean = mean(c_across(all_of(as.character(2:12))), na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::select(Year, Mean)
names(naoyear)[names(naoyear) == 'Mean'] <- 'NAO'


#COMBINE ALL ECOVARS
ecovars_snema <- bt_yearly_snema %>%
  inner_join(sst_yearly_snema, by = "Year") %>%
  inner_join(gsiyear, by = "Year") %>%
  inner_join(amoyear, by = "Year") %>%
  inner_join(naoyear, by = "Year") 

head(ecovars_snema)
str(ecovars_snema)
summary(ecovars_snema)

plot(BottomT ~Year, data=ecovars_snema, main="Exploratory Visual", xlab="Year", ylab="Degrees Lat Change & Temp", type="l", lwd=3,col="#00608A", cex.lab=1.4,cex.axis=1.1)
lines(GSI ~ Year, data=ecovars_snema, xlab="Year", type="l",col="#EA4F12",lwd=3)
lines(AMO ~ Year, data=ecovars_snema, xlab="Year", type="l",col="#00736D",lwd=3)
lines(NAO ~ Year, data=ecovars_snema, xlab="Year", type="l",col="#08229D",lwd=3)

summary(ecovars_snema[, c("GSI", "AMO", "NAO")])

#test to scale
# Normalize selected columns
ecovars_snema_scaled <- ecovars_snema
ecovars_snema_scaled$BottomT <- scale(ecovars_snema$BottomT)
ecovars_snema_scaled$GSI <- scale(ecovars_snema$GSI)
ecovars_snema_scaled$AMO <- scale(ecovars_snema$AMO)
ecovars_snema_scaled$NAO <- scale(ecovars_snema$NAO)

#plot scaled data
plot(BottomT ~ Year, data=ecovars_snema_scaled, type="l", lwd=3, col="#00608A",
     main="Exploratory Visual (Scaled)", xlab="Year", ylab="Standardized Value",
     cex.lab=1.4, cex.axis=1.1)
lines(GSI ~ Year, data=ecovars_snema_scaled, col="#EA4F12", lwd=3)
lines(AMO ~ Year, data=ecovars_snema_scaled, col="#00736D", lwd=3) 
lines(NAO ~ Year, data=ecovars_snema_scaled, col="#08229D", lwd=3)

#visual test for outliers 
dotchart(ecovars_snema$BottomT, main = "Bottom Temp")
dotchart(ecovars_snema$GSI, main = "GSI")
dotchart(ecovars_snema$AMO, main = "AMO")
dotchart(ecovars_snema$NAO, main = "NAO")
boxplot(ecovars_snema$BottomT, main = "Bottom Temp")
boxplot(ecovars_snema$GSI, main = "GSI")
boxplot(ecovars_snema$AMO, main = "AMO")
boxplot(ecovars_snema$NAO, main = "NAO")
boxplot(snema_NMFS_age1$NO_AT_AGE, main = "NMFS Recrutiment")
boxplot(snema_MADMF_age1$NO_AT_AGE, main = "MADMF Recruitment")
#check distribution
hist(ecovars_snema$BottomT, main = "Bottom Temp")
hist(ecovars_snema$GSI, main = "GSI") 
hist(ecovars_snema$AMO, main = "AMO") 
hist(ecovars_snema$NAO, main = "NAO") 
hist(snema_NMFS_age1$NO_AT_AGE, main = "NMFS Recrutiment")
hist(snema_MADMF_age1$NO_AT_AGE, main = "MADMF Recruitment")

#SUMMER FLOUNDER LAT (PREDATION/COMPETITION)
summerfl <- read.csv(here("data/summerflounderlatSPRING.csv")) #using coglat as main variable
names(summerfl)[names(summerfl) == 'YEAR'] <- 'Year'
names(summerfl)[names(summerfl) == 'COGLAT'] <- 'SFLatSpring'
summerfl <- summerfl[, c("Year", "SFLatSpring")]

shapiro_fun(summerfl) #lat is not normal, trying gaussian first
hist(summerfl$SFLatSpring)

snema_comp1_data <- snema_NMFS_age1 %>%
  inner_join(summerfl, by = "Year") %>%
  inner_join(ecovars_snema, by = "Year")
snema_comp2_data <- snema_MADMF_age1 %>%
  inner_join(summerfl, by = "Year") %>%
  inner_join(ecovars_snema, by = "Year")

#checking for co-linearity
Mypairs(snema_comp1_data) 
Mypairs(snema_comp2_data) 
corvif(snema_comp1_data) 
corvif(snema_comp2_data) 
corvif(snema_comp1_data)
corvif(snema_comp1_data[c(2:3)])
corvif(snema_comp1_data[c(2:4)])
corvif(snema_comp2_data)
corvif(snema_comp2_data[c(2:3)])
corvif(snema_comp2_data[c(2:4)])

#Shapiro Test (function)
shapiro_fun <- function(data) {
  for (i in seq_along(data)) {
    p_value <- shapiro.test(data[[i]])$p.value
    cat(names(data)[i], ": ", if (p_value > 0.05) "Normal" else "Not Normal", "\n")
  }
}

shapiro_fun(ecovars_snema) #all are normal
shapiro.test(ecovars_snema$BottomT) 
shapiro.test(ecovars_snema$GSI)
shapiro.test(ecovars_snema$AMO)
shapiro.test(ecovars_snema$NAO)

#Run functions as is, NO EDIT###############
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

panel.cor <- function(x, y, digits=1, prefix="", cex.cor = 6){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1=cor(x,y,use="pairwise.complete.obs")
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r1, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) { cex <- 0.9/strwidth(txt) } else {
    cex = cex.cor}
  text(0.5, 0.5, txt, cex = cex * r)
}
Mypairs <- function(Z) {
  MyVarx <- colnames(Z)
  pairs(Z, labels = MyVarx,
        cex.labels =  2,
        lower.panel = function(x, y, digits=2, prefix="", cex.cor = 7) {
          panel.cor(x, y, digits, prefix, cex.cor)},
        upper.panel =  function(x, y) points(x, y,
                                             pch = 16, cex = 0.8,
                                             col = gray(0.1)))
}


######run functions with data
Mypairs(ecovars_snema) 
corvif(ecovars_snema) 
corvif(ecovars_snema[c(2:4)])
corvif(ecovars_snema[c(3:4)])
corvif(ecovars_snema[c(1,4)])

#GAMS RUN

snema_gam1data <- snema_NMFS_age1 %>%
  inner_join(ecovars_snema, by = "Year")

snema_gam1data <- snema_gam1data[order(snema_gam1data$Year),]

snema_gam2data <- snema_MADMF_age1 %>%
  inner_join(ecovars_snema, by = "Year")

snema_gam2data <- snema_gam2data[order(snema_gam2data$Year),]

#s(BottomT, k=5) + s(SurfaceT, k=5) + s(GSI, k=5) + s(AMO, k=5) + s(NAO, k=5)


#NMFS GAUSSIAN
snema_gam1 <- gam(log(NO_AT_AGE)~s(SurfaceT, k=5), family=gaussian(), method="REML", data=snema_gam1data)
#REML: Restricted maximum likelihood approach to smoothing
summary(snema_gam1)
gam.check(snema_gam1) 
concurvity(snema_gam1)
AIC(snema_gam1) 
#RUN 1, AIC 81.96, P= s(BottomT) 0.504, s(SurfaceT) 0.144, s(GSI) 0.358, s(AMO) 0.164, s(NAO) 0.356 => DELETING BOTTOM T
#RUN 2, AIC 80.48, P= s(SurfaceT) 0.144, s(GSI) 0.353, s(AMO) 0.157, s(NAO) 0.416 => DELETING NAO
#RUN 3, AIC 79.47, P= s(SurfaceT) 0.168, s(GSI) 0.420, s(AMO) 0.234 => DELETING GSI
#RUN 4, AIC 78.44, P= s(SurfaceT) 0.0974 ., s(AMO) 0.3628  => DELETING AMO
#RUN 5, AIC 77.70, P= s(SurfaceT) 0.229 => NOTHING SIGNIFICANT 

#NMFS GAMMA
snema_gamma1 <- gam(NO_AT_AGE~s(AMO, k=5), family=Gamma(), method="REML", data=snema_gam1data)
#REML: Restricted maximum likelihood approach to smoothing
summary(snema_gamma1)
gam.check(snema_gamma1) 
concurvity(snema_gamma1)
AIC(snema_gamma1) 
#RUN 1, AIC -45.79, P= s(BottomT) 0.3831, s(SurfaceT) 0.0713 ., s(GSI) 0.3747,s(AMO) 0.0626 ., s(NAO) 0.4863 => DELETING NAO
#RUN 2, AIC -47.42, P= s(BottomT) 0.4482, s(SurfaceT) 0.1423, s(GSI) 0.2221, s(AMO) 0.0872 . => DELETING BOTTOM T
#RUN 3, AIC -48.87, P= s(SurfaceT) 0.1907,s(GSI) 0.2249, s(AMO) 0.0957 . => DELETING GSI
#RUN 4, AIC -49.23, P= s(SurfaceT) 0.425, s(AMO) 0.249 => DELETING SURFACE T
#RUN 5, AIC -50.52, P= s(AMO) 0.36 => NOTHING SIGNIFICANT

#NMFS TWEEDIE
snema_tw1 <- gam(NO_AT_AGE~s(SurfaceT, k=5), family=tw(), method="REML", data=snema_gam1data)
#REML: Restricted maximum likelihood approach to smoothing
summary(snema_tw1)
gam.check(snema_tw1) 
concurvity(snema_tw1)
AIC(snema_tw1) 
#RUN 1, AIC -44.07, P= s(BottomT) 0.4735, s(SurfaceT) 0.0826 ., s(GSI) 0.2405, s(AMO) 0.0607 ., s(NAO) 0.3185  => DELETING BOTTOM T
#RUN 2, AIC -45.59, P= s(SurfaceT) 0.0624 ., s(GSI) 0.2415, s(AMO) 0.0601 ., s(NAO) 0.3901  => DELETING NAO
#RUN 3, AIC -47.81, P= s(SurfaceT) 0.0521 ., s(GSI) 0.1326, s(AMO) 0.1019  => DELETING GSI
#RUN 4, AIC -47.82, P= s(SurfaceT) 0.208, s(AMO) 0.282 => DELETING AMO
#RUN 5, AIC -48.21, P= s(SurfaceT) 0.435, => NOTHING SIGNIFICANT

#MADMF GAUSSIAN 
snema_gam2 <- gam(log(NO_AT_AGE)~s(SurfaceT, k=5), family=gaussian(), method="REML", data=snema_gam2data)
#REML: Restricted maximum likelihood approach to smoothing
summary(snema_gam2)
gam.check(snema_gam2)
concurvity(snema_gam2)
AIC(snema_gam2) 
#RUN 1, AIC 77.54, P= s(BottomT) 0.822, s(SurfaceT) 0.253, s(GSI) 0.518, s(AMO) 0.598, s(NAO) 0.595 => DELETING BOTTOM T
#RUN 2, AIC 75.69, P= s(SurfaceT) 0.0591 ., s(GSI) 0.5074, s(AMO) 0.5446, s(NAO) 0.5446  => DELETING NAO AND AMO
#RUN 3, AIC 70.57, P= s(SurfaceT) 0.0494 *, s(GSI) 0.3250 => DELETING GSI
#RUN 4, AIC 71.68, P= s(SurfaceT) 0.0138 * => SURFACE T SIGNIFICANT

#MADMF GAMMA
snema_gamma2 <- gam(NO_AT_AGE~s(SurfaceT, k=5) + s(GSI, k=5) + s(NAO, k=5), family=Gamma(), method="REML", data=snema_gam2data)
#REML: Restricted maximum likelihood approach to smoothing
summary(snema_gamma2)
gam.check(snema_gamma2)
concurvity(snema_gamma2)
AIC(snema_gamma2)
#RUN 1, AIC 149.22, P= s(BottomT) 0.8215, s(SurfaceT) 0.0686 ., s(GSI) 0.0712 ., s(AMO) 0.6541,  s(NAO) 0.2495  => DELETING BOTTOM T
#RUN 2, AIC 146.89, P= s(SurfaceT) 0.0173 *, s(GSI) 0.0626 ., s(AMO) 0.6725, s(NAO) 0.2361  => DELETING AMO
#RUN 3, AIC 147.65, P= s(SurfaceT) 0.0160 *, s(GSI) 0.0565 ., s(NAO) 0.0688 . => KEEP ALL?

#MADMF TWEEDIE
snema_tw2 <- gam(NO_AT_AGE~s(SurfaceT, k=5), family=tw(), method="REML", data=snema_gam2data)
#REML: Restricted maximum likelihood approach to smoothing
summary(snema_tw2)
gam.check(snema_tw2) 
concurvity(snema_tw2)
AIC(snema_tw2) 
#RUN 1, AIC 151.98, P= s(BottomT) 0.9686, s(SurfaceT) 0.0649 ., s(GSI) 0.1337, s(AMO) 0.5107, s(NAO) 0.7087 => DELETING BOTTOM T
#RUN 2, AIC 149.74, P= s(SurfaceT) 0.00966 **, s(GSI) 0.12503, s(AMO) 0.48116, s(NAO) 0.72435  => DELETING NAO
#RUN 3, AIC 147.61, P= s(SurfaceT) 0.00762 **, s(GSI) 0.13509, s(AMO) 0.21357 => DELETING AMO
#RUN 4, AIC 148.20, P= s(SurfaceT) 0.0114 *, s(GSI) 0.2458 => DELETING GSI
#RUN 4, AIC 147.61, P= s(SurfaceT) 0.00881 ** => SURFACE T SIGNIFICANT

#PLOT WINNER GAUSSIAN MADMF
png(here("Figures/Raw_data_trends/GAMS/SNEMA/BTRecruitmentMADMF.png"),
    width = 2400, height = 2000, units = "px", res = 300)
plot.gam(
  snema_gam2,
  xlab = "Bottom Temperature (°C)",
  ylab = "Partial Effect",
  select = 1,
  cex.lab = 1.5,
  cex.axis = 1.4,
  rug = TRUE,
  shade = TRUE,
  col = "black",
  shade.col = "#E9E9E9",
  lwd = 2,
  main = "Partial Effect of Bottom Temperature on SNEMA Age-1 Recruitment (MADMF)"
)
dev.off()

#---------------------------------------
#Let's see if Summer Flounder spring latitude has an effect on WF recruitment in addition to BottomT

summerfl <- read.csv(here("data/summerflounderlatSPRING.csv")) #using coglat as main variable
names(summerfl)[names(summerfl) == 'YEAR'] <- 'Year'
names(summerfl)[names(summerfl) == 'COGLAT'] <- 'SFLatSpring'
summerfl <- summerfl[, c("Year", "SFLatSpring")]

shapiro_fun(summerfl) #lat is not normal, trying gaussian first
hist(summerfl$SFLatSpring)

snema_comp1_data <- snema_NMFS_age1 %>%
  inner_join(summerfl, by = "Year") %>%
  inner_join(ecovars_snema, by = "Year")
snema_comp2_data <- snema_MADMF_age1 %>%
  inner_join(summerfl, by = "Year") %>%
  inner_join(ecovars_snema, by = "Year")

#checking for co-linearity
Mypairs(snema_comp1_data) 
Mypairs(snema_comp2_data) 
corvif(snema_comp1_data) 
corvif(snema_comp2_data) 
corvif(snema_comp1_data)
corvif(snema_comp1_data[c(2:3)])
corvif(snema_comp1_data[c(2:4)])
corvif(snema_comp2_data)
corvif(snema_comp2_data[c(2:3)])
corvif(snema_comp2_data[c(2:4)])

#s(BottomT, k=5) + s(SurfaceT, k=5) +s(SFLatSpring, k=5)+ s(GSI, k=5) + s(AMO, k=5) + s(NAO, k=5)

#NMFS GAUSSIAN
snema_compgam1 <- gam(log(NO_AT_AGE)~s(BottomT, k=5)+s(SurfaceT, k=5)+ s(SFLatSpring, k=5)+s(GSI, k=5)+s(NAO, k=5)+s(AMO, k=5), family=gaussian(), method="REML", data=snema_comp1_data)
summary(snema_compgam1)
gam.check(snema_compgam1)
concurvity(snema_compgam1)
AIC(snema_compgam1) 
#RUN 1, AIC , P= 
#RUN 2, AIC , P= 
#RUN 3, AIC , P= 
#RUN 4, AIC , P= 

#NMFS GAMMA
snema_compgamma1 <- gam(NO_AT_AGE~s(BottomT, k=5)+s(SurfaceT, k=5)+s(SFLatSpring, k=5)+s(GSI, k=5)+s(NAO, k=5)+s(AMO, k=5), family=Gamma(), method="REML", data=snema_comp1_data)
summary(snema_compgamma1)
gam.check(snema_compgamma1)
concurvity(snema_compgamma1)
AIC(snema_compgamma1)
#RUN 1, AIC , P= 
#RUN 2, AIC , P= 
#RUN 3, AIC , P= 
#RUN 4, AIC , P= 

#NMFS TWEEDIE 
snema_comptw1 <- gam(NO_AT_AGE~s(BottomT, k=5)+s(SurfaceT, k=5)+s(SFLatSpring, k=5)+s(GSI, k=5)+s(NAO, k=5)+s(AMO, k=5), family=tw(), method="REML", data=snema_comp1_data)
summary(snema_comptw1)
gam.check(snema_comptw1)
concurvity(snema_comptw1)
AIC(snema_comptw1)
#RUN 1, AIC , P= 
#RUN 2, AIC , P= 
#RUN 3, AIC , P= 
#RUN 4, AIC , P= 

#MADMF GAUSSIAN
snema_compgam2 <- gam(log(NO_AT_AGE)~s(BottomT, k=5)+s(SurfaceT, k=5)+s(SFLatSpring, k=5)+s(GSI, k=5)+s(NAO, k=5)+s(AMO, k=5), family=gaussian(), method="REML", data=snema_comp2_data)
summary(snema_compgam2)
gam.check(snema_compgam2)
concurvity(snema_compgam2)
AIC(snema_compgam2)
#RUN 1, AIC , P= 
#RUN 2, AIC , P= 
#RUN 3, AIC , P= 
#RUN 4, AIC , P= 

#MADMF GAMMA
snema_compgamma2 <- gam(NO_AT_AGE~s(BottomT, k=5)+s(SurfaceT, k=5)+s(SFLatSpring, k=5)+s(GSI, k=5)+s(NAO, k=5)+s(AMO, k=5), family=Gamma(), method="REML", data=snema_comp2_data)
summary(snema_compgamma2)
gam.check(snema_compgamma2)
concurvity(snema_compgamma2)
AIC(snema_compgamma2)
#RUN 1, AIC , P= 
#RUN 2, AIC , P= 
#RUN 3, AIC , P= 
#RUN 4, AIC , P= 

#MADMF TWEEDIE
snema_comptw2 <- gam(NO_AT_AGE~s(BottomT, k=5)+s(SurfaceT, k=5)+s(SFLatSpring, k=5)+s(GSI, k=5)+s(NAO, k=5)+s(AMO, k=5), family=tw(), method="REML", data=snema_comp2_data)
summary(snema_comptw2)
gam.check(snema_comptw2)
concurvity(snema_comptw2)
AIC(snema_comptw2)
#RUN 1, AIC , P= 
#RUN 2, AIC , P= 
#RUN 3, AIC , P= 
#RUN 4, AIC , P= 


png(here("Figures/Raw_data_trends/GAMS/SNEMA/SFlatonAge1MADMF.png"),
    width = 2600, height = 2000, units = "px", res = 300)
plot.gam(
  snema_compgam2,
  xlab = "Latitude",
  ylab = "Partial Effect",
  select = 1,
  cex.lab = 1.5,
  cex.axis = 1.4,
  rug = TRUE,
  shade = TRUE,
  col = "black",
  shade.col = "#E9E9E9",
  lwd = 2,
  main = "Partial Effect of Summer Flounder Latitude on SNEMA WF Age-1 Recruitment (MADMF)" 
)
dev.off()



#-------------------------------------------------------------------------------
#RECRUITS PER SPAWNER

snema_rssb <- rssbdf[, c("Year", "SNEMA_Fall_rssb", "SNEMA_Spring_rssb")]

snema_rssb_gamdata1 <- snema_rssb %>%
  inner_join(ecovars_snema, by = "Year")

snema_rssb_gamdata1 <- snema_rssb_gamdata1[order(snema_rssb_gamdata1$Year),]

#s(BottomT, k=5) + s(SurfaceT, k=5) + s(GSI, k=5) + s(AMO, k=5) + s(NAO, k=5)

#FALL GAM
snema_fallrssb_gam <- gam(log(SNEMA_Fall_rssb)~s(AMO, k=5), family=gaussian(), method="REML", data=snema_rssb_gamdata1)
#REML: Restricted maximum likelihood approach to smoothing
summary(snema_fallrssb_gam)
gam.check(snema_fallrssb_gam) 
concurvity(snema_fallrssb_gam)
AIC(snema_fallrssb_gam) 
#RUN 1, AIC= 88.80, P= s(BottomT) 0.880, s(SurfaceT) 0.674, s(GSI) 0.551, s(AMO) 0.303, s(NAO) 0.982 => DELETING NAO 
#RUN 2, AIC= 86.83, P= s(BottomT) 0.867, s(SurfaceT) 0.651, s(GSI) 0.542, s(AMO) 0.236 => DELETING BOTTOM T
#RUN 3, AIC= 84.68, P= s(SurfaceT) 0.673, s(GSI) 0.540, s(AMO) 0.223 => DELETING SURFACE T
#RUN 4, AIC= 82.89, P= s(GSI) 0.6413, s(AMO) 0.0819 . => DELETING GSI
#RUN 5, AIC= 81.14, P= s(AMO) 0.0746 . => AMO BARELY SIGNIFICANT


#SPRING GAM
snema_springrssb_gam <- gam(log(SNEMA_Spring_rssb)~s(SurfaceT, k=5), family=gaussian(), method="REML", data=snema_rssb_gamdata1)
#REML: Restricted maximum likelihood approach to smoothing
summary(snema_springrssb_gam)
gam.check(snema_springrssb_gam) 
concurvity(snema_springrssb_gam)
AIC(snema_springrssb_gam) 
#RUN 1, AIC= 84.68, P= s(BottomT) 0.666, s(SurfaceT) 0.380, s(GSI) 0.870, s(AMO) 0.251, s(NAO) 0.948 => DELETING NAO
#RUN 2, AIC= 82.72, P= s(BottomT) 0.644, s(SurfaceT) 0.347, s(GSI) 0.868, s(AMO) 0.183 => DELETING GSI
#RUN 3, AIC= 80.75, P= s(BottomT) 0.634, s(SurfaceT) 0.182, s(AMO) 0.168 => DELETING BOTTOM T
#RUN 4, AIC= 79.02, P= s(SurfaceT) 0.0751 ., s(AMO) 0.1483 => DELETING AMO
#RUN 5, AIC= 80.22, P= s(SurfaceT) 0.0365 * => SURFACE T SIGNIFICANT

#png(here("Figures/Raw_data_trends/GAMS/SNEMA/AMO_GSI_onRSSB_Spring.png"),
    #width = 2600, height = 2000, units = "px", res = 300)
#plot.gam(
  #snema_springrssb_gam,
  #xlab = "GSI",
  #ylab = "Partial Effect",
  #select = 2,
  #cex.lab = 1.5,
  #cex.axis = 1.4,
  #rug = TRUE,
  #shade = TRUE,
  #col = "black",
  #shade.col = "#E9E9E9",
  #lwd = 2,
  #main = "Partial Effect of GSI on SNEMA Spring Recruits per Spawner" 
#)
#dev.off()

#-------------------------------------------------------------------------------
#CONDITION 
#Keep in mind that the COND data set uses the Mid Atlantic Bight NOT snema! 

print(cond_sne)
#mean condition = units
cond_sne <- cond_sne[, c("Year", "Value")]
cond_sne <- cond_sne %>%
  rename(MeanCondition = Value)

snema_cond_gamdata <- cond_sne %>%
  inner_join(ecovars_snema, by = "Year")

snema_cond_gamdata <- snema_cond_gamdata[order(snema_cond_gamdata$Year),]

#s(BottomT, k=5) + s(SurfaceT, k=5) + s(GSI, k=5) + s(AMO, k=5) + s(NAO, k=5)
snema_cond_gam <- gam(MeanCondition~s(BottomT, k=5) + s(SurfaceT, k=5), family=gaussian(), method="REML", data=snema_cond_gamdata)
summary(snema_cond_gam)
gam.check(snema_cond_gam) 
concurvity(snema_cond_gam)
AIC(snema_cond_gam) 
#RUN 1 AIC=-117.40, p= s(BottomT) 0.0333 *, s(SurfaceT) 0.0989 ., s(GSI) 0.5138, s(AMO) 0.5591, s(NAO) 0.5902 => DELETING NAO
#RUN 2 AIC=-119.02, p= s(BottomT) 0.0198 *, s(SurfaceT) 0.0561 ., s(GSI) 0.5179, s(AMO) 0.3152 => DELETING GSI
#RUN 3 AIC=-120.50, p= s(BottomT) 0.01594 *, s(SurfaceT) 0.00498 **, s(AMO) 0.14085 => DELETING AMO
#RUN 4 AIC=-119.92, p= s(BottomT) 0.0196 *, s(SurfaceT) 0.0106 * => BOTTOM T AND SURFACE T SIGNIFICANT

png(here("Figures/Raw_data_trends/GAMS/SNEMA/BottomT_GSI_onMABCondition.png"),
    width = 2600, height = 2000, units = "px", res = 300)
plot.gam(
  snema_cond_gam,
  xlab = "GSI",
  ylab = "Partial Effect",
  select = 2,
  cex.lab = 1.5,
  cex.axis = 1.4,
  rug = TRUE,
  shade = TRUE,
  col = "black",
  shade.col = "#E9E9E9",
  lwd = 2,
  main = "Partial Effect of GSI on Mid-Atlantic Bight (SNEMA) Mean Condition" 
)
dev.off()

#-------------------------------------------------------------------------------
#DISTRIBUTION
#make sure to include ssb from ssb_snema file
#distribution too: cogdf

#get ready to join
names(ssb_snema)[names(ssb_snema) == 'year'] <- 'Year'
names(cogdf)[names(cogdf) == 'year'] <- 'Year'
cogdf <- cogdf[, c("Year", "fall_cog_lat", "fall_cog_depth", "spring_cog_lat", "spring_cog_depth")]

snema_distribution_gamdata <- ssb_snema %>%
  inner_join(ecovars_snema, by = "Year") %>%
  inner_join(summerfl, by = "Year") %>%
  inner_join(cogdf, by = "Year") 

#checking for co-linearity
Mypairs(snema_distribution_gamdata) 

hist(snema_distribution_gamdata$fall_cog_lat) #right skewed gaussian logged or gamma
hist(snema_distribution_gamdata$spring_cog_lat) #gaussian or gamma
hist(snema_distribution_gamdata$fall_cog_depth) #gaussian or gamma
hist(snema_distribution_gamdata$spring_cog_depth) #gaussian or gamma

#GAUSSIAN AUTO GAM
#------------------------------------------------------------------------------
library(mgcv)
library(ggplot2)

auto_gam_gaussian <- function(response, predictors, data, family = gaussian(), k = 10,
                              p_thresh = 0.07, max_steps = 20, verbose = TRUE,
                              plot_path = "figures/Raw_data_trends/GAMS/SNEMA/") {
  
  remaining_predictors <- predictors
  results <- list()
  step <- 1
  
  while (length(remaining_predictors) > 0 && step <= max_steps) {
    
    # Safety check
    if (length(remaining_predictors) == 0) {
      if (verbose) message("No predictors left to model. Exiting loop.")
      break
    }
    
    # Build the formula with smooth terms
    smooth_terms <- paste0("s(", remaining_predictors, ", k=", k, ")", collapse = " + ")
    formula_text <- paste(response, "~", smooth_terms)
    model_formula <- as.formula(formula_text)
    
    start_time <- Sys.time()
    
    # Fit the GAM
    model <- gam(model_formula, data = data, family = family, method = "REML")
    model_summary <- summary(model)
    aic_value <- AIC(model)
    
    # Extract p-values for smooth terms (full precision for logic)
    pvals <- model_summary$s.table[, "p-value"]
    pvals_rounded <- round(pvals, 4)
    
    # Get k-index diagnostics
    kcheck <- suppressWarnings(gam.check(model, rep = 100, verbose = FALSE))
    k_index_table <- kcheck$k.check  # matrix with k-index and p-values
    
    end_time <- Sys.time()
    elapsed <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
    
    # Store results
    results[[step]] <- list(
      step = step,
      formula = formula_text,
      AIC = aic_value,
      p_values = pvals,
      k_index = k_index_table,
      model = model
    )
    
    # ---- Output ----
    if (verbose) {
      message("---- Run ", step, " ----")
      message("Formula: ", formula_text)
      message("AIC: ", round(aic_value, 2))
      message("P-values (rounded):")
      print(pvals_rounded)
      
      # Print k-index warnings
      message("Checking k-index diagnostics:")
      for (i in seq_len(nrow(k_index_table))) {
        k_row <- k_index_table[i, ]
        if (k_row["k-index"] < 1 && k_row["p-value"] < 0.05) {
          message("⚠️  ", rownames(k_index_table)[i], ": k-index = ", round(k_row["k-index"], 3),
                  " with p = ", signif(k_row["p-value"], 4), " → consider increasing k")
        }
      }
      message("Step runtime: ", elapsed, " seconds\n")
    }
    
    # ---- Drop term logic ----
    if (all(pvals < p_thresh)) {
      if (verbose) message("✅ All terms significant. Stopping iteration.")
      break
    }
    
    max_pval <- max(pvals)
    worst_term <- names(pvals)[which.max(pvals)]
    # Extract predictor name from smooth term (s(variable,k=...))
    worst_var <- sub("s\\(([^,]+),.*", "\\1", worst_term)
    
    if (max_pval >= p_thresh) {
      if (worst_var %in% remaining_predictors) {
        if (verbose) message("Dropping term '", worst_var, "' with p = ", round(max_pval, 4), "\n")
        remaining_predictors <- setdiff(remaining_predictors, worst_var)
      } else {
        warning("⚠️ Variable '", worst_var, "' not found in predictor list. Stopping to avoid infinite loop.")
        break
      }
    } else {
      if (verbose) message("No non-significant term above threshold. Stopping iteration.")
      break
    }
    
    step <- step + 1
  }
  
  # ---- Plot significant smooths from final model ----
  final_model <- results[[length(results)]]$model
  final_pvals <- results[[length(results)]]$p_values
  
  sig_terms <- names(final_pvals)[final_pvals < p_thresh]
  
  if (length(sig_terms) > 0) {
    if (!dir.exists(plot_path)) dir.create(plot_path, recursive = TRUE)
    
    for (term in sig_terms) {
      # Clean term name for filename
      clean_term <- gsub("s\\(|\\)", "", term)
      plot_filename <- paste0(plot_path, response, "_step", step, "_", clean_term, "_GAMplot.png")
      
      png(plot_filename, width = 800, height = 600)
      plot(final_model, select = which(names(final_pvals) == term),
           shade = TRUE, shade.col = "lightblue",
           main = paste0("Smooth for ", clean_term))
      dev.off()
    }
    
    if (verbose) message("📊 Significant smooth plots saved to: ", plot_path)
  }
  
  return(results)
}


#AUTO GAMMA GAM
#------------------------------------------------------------------------------

auto_gam_gamma <- function(response, predictors, data, family = Gamma(), k = 5,
                              p_thresh = 0.07, max_steps = 20, verbose = TRUE) {
  
  remaining_predictors <- predictors
  results <- list()
  step <- 1
  
  while (length(remaining_predictors) > 0 && step <= max_steps) {
    
    # Build the GAM formula dynamically
    smooth_terms <- paste0("s(", remaining_predictors, ", k=", k, ")", collapse = " + ")
    formula_text <- paste(response, "~", smooth_terms)
    model_formula <- as.formula(formula_text)
    
    # Fit the model
    model <- gam(model_formula, data = data, family = family, method = "REML")
    model_summary <- summary(model)
    aic_value <- AIC(model)
    
    # Extract p-values for smooth terms
    pvals <- model_summary$s.table[, "p-value"]
    
    # Save results
    results[[step]] <- list(
      step = step,
      formula = formula_text,
      AIC = aic_value,
      p_values = round(pvals, 4),
      model = model
    )
    
    # Print info
    if (verbose) {
      cat("---- Run", step, "----\n")
      cat("Formula:", formula_text, "\n")
      cat("AIC:", round(aic_value, 2), "\n")
      print(round(pvals, 4))
      cat("\n")
    }
    
    # Check if all terms are below threshold
    if (all(pvals < p_thresh)) {
      if (verbose) cat("All terms significant. Stopping.\n")
      break
    }
    
    # Identify the term with the highest p-value
    max_pval <- max(pvals)
    worst_term <- names(pvals)[which.max(pvals)]  # e.g. "s(SurfaceT)"
    worst_var <- sub("^s\\(([^,]+)\\).*", "\\1", worst_term)
    
    
    if (max_pval >= p_thresh) {
      if (verbose) cat("Dropping:", worst_var, "with p =", round(max_pval, 4), "\n\n")
      # Remove the worst predictor from the list using logical indexing
      if (worst_var %in% remaining_predictors) {
        remaining_predictors <- remaining_predictors[remaining_predictors != worst_var]
      } else {
        if (verbose) cat("Warning: Variable", worst_var, "not found in predictors. Stopping to avoid infinite loop.\n")
        break
      }
    } else {
      if (verbose) cat("No term above threshold. Stopping.\n")
      break
    }
    
    # Stop if no predictors left
    if (length(remaining_predictors) == 0) {
      if (verbose) cat("No predictors left. Stopping.\n")
      break
    }
    
    step <- step + 1
  }
  
  return(results)
}
#-------------------------------------------------------------------------------

#s(BottomT, k=5) + s(SurfaceT, k=5) s(SFLatSpring, k=5) + s(GSI, k=5) + s(AMO, k=5) + s(NAO, k=5) + s(SSB, k=5)
#SPRING LAT GAUSSIAN
snema_distgamlat <- auto_gam_gaussian(
  response = "spring_cog_lat",
  predictors = c("BottomT", "SurfaceT", "SFLatSpring", "GSI", "AMO", "NAO", "SSB"),
  data = snema_distribution_gamdata,
  p_thresh = 0.07
)
#Formula: spring_cog_lat ~ s(BottomT, k=5) + s(SSB, k=5) 
#AIC: -21.53 
#s(BottomT)     s(SSB) 
#0.0453     0.0004

#SPRING LAT GAMMA
snema_distgammalat <- auto_gam_gamma(
  response = "spring_cog_lat",
  predictors = c("BottomT", "SurfaceT", "SFLatSpring", "GSI", "AMO", "NAO", "SSB"),
  data = snema_distribution_gamdata,
  p_thresh = 0.07
)
#Formula: spring_cog_lat ~ s(BottomT, k=5) + s(SSB, k=5) 
#AIC: -21.1 
#s(BottomT)     s(SSB) 
#0.0451     0.0004 

png(here("Figures/Raw_data_trends/GAMS/SNEMA/SSB_SpringLatDist.png"),
    width = 2600, height = 2000, units = "px", res = 300)
plot.gam(
  snema_distgammalat,
  xlab = "SSB",
  ylab = "Partial Effect",
  select = 2,
  cex.lab = 1.5,
  cex.axis = 1.4,
  rug = TRUE,
  shade = TRUE,
  col = "black",
  shade.col = "#E9E9E9",
  lwd = 2,
  main = "Partial Effect of SNEMA SSB on Spring COG Latitude " 
)
dev.off()

png(here("Figures/Raw_data_trends/GAMS/SNEMA/BT_SpringLatDist.png"),
    width = 2600, height = 2000, units = "px", res = 300)
plot.gam(
  snema_distgammalat,
  xlab = "Bottom Temperature (°C)",
  ylab = "Partial Effect",
  select = 1,
  cex.lab = 1.5,
  cex.axis = 1.4,
  rug = TRUE,
  shade = TRUE,
  col = "black",
  shade.col = "#E9E9E9",
  lwd = 2,
  main = "Partial Effect of SNEMA Bottom Temp on Spring COG Latitude " 
)
dev.off()

#SPRING DEPTH GAUSSIAN
snema_distgamdepth <- auto_gam_gaussian(
  response = "spring_cog_depth",
  predictors = c("BottomT", "SurfaceT", "GSI", "AMO", "NAO", "SSB"),
  data = snema_distribution_gamdata,
  p_thresh = 0.07
)
#Formula: spring_cog_depth ~ s(BottomT, k=5) + s(SSB, k=5) 
#AIC: 167.5 
#s(BottomT)     s(SSB) 
#0.0151     0.0605 


#SPRING DEPTH GAMMA
snema_distgammadepth <- gam(abs(spring_cog_depth)~s(BottomT, k=5), family=Gamma(), method="REML", data=snema_distribution_gamdata)
summary(snema_distgammadepth)
gam.check(snema_distgammadepth)
concurvity(snema_distgammadepth)
AIC(snema_distgammadepth)
#RUN 1, AIC 173.94, P= s(BottomT) 0.0276 *, s(SSB) 0.1829, s(GSI) 0.3936, s(NAO) 0.2574, s(AMO) 0.8883 => REMOVING AMO
#RUN 2, AIC 171.94, P= s(BottomT) 0.0171 *, s(SSB) 0.1446, s(GSI) 0.4099, s(NAO) 0.2456  => REMOVING GSI
#RUN 3, AIC 169.65, P= s(BottomT) 0.0124 *, s(SSB) 0.0458 *, s(NAO) 0.3247 => REMOVING NAO 
#RUN 4, AIC 168.80, P= s(BottomT) 0.0155 *, s(SSB) 0.0596 . => REMOVING SSB
#RUN 5, AIC 172.32, P= s(BottomT) 0.23 => NOTHING SIGNIFICANT

#FALL LAT GAUSSIAN
snema_distgamlat1 <- gam(log(fall_cog_lat)~s(SSB, k=5)+s(GSI, k=5), family=gaussian(), method="REML", data=snema_distribution_gamdata)
summary(snema_distgamlat1)
gam.check(snema_distgamlat1)
concurvity(snema_distgamlat1)
AIC(snema_distgamlat1) 
#RUN 1, AIC -259.15, P= s(BottomT) 0.085 ., s(SSB) 0.770, s(GSI) 0.270, s(NAO) 0.790, s(AMO) 0.850 => REMOVING AMO
#RUN 2, AIC -258.41, P= s(BottomT) 0.7676, s(SSB) 4.48e-06 ***, s(GSI) 0.0903 ., s(NAO) 0.6948  => REMOVING BOTTOM T
#RUN 3, AIC -260.30, P= s(SSB) 1.29e-06 ***, s(GSI) 0.0585 ., s(NAO) 0.6141 => REMOVING NAO  
#RUN 4, AIC -261.85, P= s(SSB) 1.13e-06 ***, s(GSI) 0.05 . => SSB AND GSI SIGNIFICANT 

#FALL LAT GAMMA
snema_distgammalat1 <- gam(fall_cog_lat~s(SSB, k=5)+s(GSI, k=5), family=Gamma(), method="REML", data=snema_distribution_gamdata)
summary(snema_distgammalat1)
gam.check(snema_distgammalat1)
concurvity(snema_distgammalat1)
AIC(snema_distgammalat1)
#RUN 1, AIC -48.83, P= s(BottomT) 0.095 ., s(SSB) 0.725, s(GSI) 0.250, s(NAO) 0.705, s(AMO) 0.905 => REMOVING AMO
#RUN 2, AIC -48.27, P= s(BottomT) 0.7684, s(SSB) 4.51e-06 ***, s(GSI) 0.0884 ., s(NAO) 0.6955 => REMOVING BOTTOM T 
#RUN 3, AIC -50.47, P= s(SSB) 1.3e-06 ***, s(GSI) 0.0571 ., s(NAO) 0.6147 => REMOVING NAO
#RUN 4, AIC -52.56, P= s(SSB) 1.14e-06 ***, s(GSI) 0.0487 * => SSB AND GSI SIGNIFICANT

png(here("Figures/Raw_data_trends/GAMS/SNEMA/SSB_FallLatDist.png"),
    width = 2600, height = 2000, units = "px", res = 300)
plot.gam(
  snema_distgammalat1,
  xlab = "SSB",
  ylab = "Partial Effect",
  select = 1,
  cex.lab = 1.5,
  cex.axis = 1.4,
  rug = TRUE,
  shade = TRUE,
  col = "black",
  shade.col = "#E9E9E9",
  lwd = 2,
  main = "Partial Effect of SNEMA SSB on Fall COG Latitude " 
)
dev.off()

png(here("Figures/Raw_data_trends/GAMS/SNEMA/GSI_FallLatDist.png"),
    width = 2600, height = 2000, units = "px", res = 300)
plot.gam(
  snema_distgammalat1,
  xlab = "GSI",
  ylab = "Partial Effect",
  select = 2,
  cex.lab = 1.5,
  cex.axis = 1.4,
  rug = TRUE,
  shade = TRUE,
  col = "black",
  shade.col = "#E9E9E9",
  lwd = 2,
  main = "Partial Effect of GSI on Fall COG Latitude " 
)
dev.off()

#FALL DEPTH GAUSSIAN
snema_distgamdepth1 <- gam(fall_cog_depth~s(SSB, k=5), family=gaussian(), method="REML", data=snema_distribution_gamdata)
summary(snema_distgamdepth1)
gam.check(snema_distgamdepth1)
concurvity(snema_distgamdepth1)
AIC(snema_distgamdepth1) 
#RUN 1, AIC 166.04, P= s(BottomT) 0.91703, s(SSB) 0.00626 **, s(GSI) 0.79724, s(NAO) 0.46555, s(AMO) 0.11778 => REMOVING BT
#RUN 2, AIC 163.92, P= s(SSB) 0.00449 **, s(GSI) 0.81550, s(NAO) 0.41136, s(AMO) 0.07572 . => REMOVING GSI
#RUN 3, AIC 161.82, P= s(SSB) 0.000242 ***, s(NAO) 0.356372, s(AMO) 0.063714 . => REMOVING NAO  
#RUN 4, AIC 162.68, P= s(SSB) 0.000184 ***, s(AMO) 0.056811 . => REMOVING AMO
#RUN 5, AIC 166.68, P= s(SSB) 2.29e-05 *** => SSB SIGNIFICANT



png(here("Figures/Raw_data_trends/GAMS/SNEMA/SSB_FallDepthDist.png"),
    width = 2600, height = 2000, units = "px", res = 300)
plot.gam(
  snema_distgamdepth1,
  xlab = "SSB",
  ylab = "Partial Effect",
  select = 1,
  cex.lab = 1.5,
  cex.axis = 1.4,
  rug = TRUE,
  shade = TRUE,
  col = "black",
  shade.col = "#E9E9E9",
  lwd = 2,
  main = "Partial Effect of SSB on Fall COG Depth " 
)
dev.off()

#FALL DEPTH GAMMA
snema_distgammadepth1 <- gam(abs(fall_cog_depth)~s(SSB, k=5)+s(AMO, k=5), family=Gamma(), method="REML", data=snema_distribution_gamdata)
summary(snema_distgammadepth1)
gam.check(snema_distgammadepth1)
concurvity(snema_distgammadepth1)
AIC(snema_distgammadepth1)
#RUN 1, AIC 163.58, P= s(BottomT) 0.96563, s(SSB) 0.00603 **, s(GSI) 0.79705, s(NAO) 0.33771, s(AMO) 0.04929 * => REMOVING BOTTOM T
#RUN 2, AIC 161.01, P= s(SSB) 0.00438 **, s(GSI) 0.77450, s(NAO) 0.28492, s(AMO) 0.02803 * => REMOVING GSI
#RUN 3, AIC 158.73, P= s(SSB) 0.000138 ***, s(NAO) 0.246134, s(AMO) 0.022391 * => REMOVING NAO  
#RUN 4, AIC 159.83, P= s(SSB) 0.000153 ***, s(AMO) 0.026809 * => SSB AND AMO SIGNIFICANT 

png(here("Figures/Raw_data_trends/GAMS/SNEMA/SSB_FallDepthDist.png"),
    width = 2600, height = 2000, units = "px", res = 300)
plot.gam(
  snema_distgammadepth1,
  xlab = "SSB",
  ylab = "Partial Effect",
  select = 1,
  cex.lab = 1.5,
  cex.axis = 1.4,
  rug = TRUE,
  shade = TRUE,
  col = "black",
  shade.col = "#E9E9E9",
  lwd = 2,
  main = "Partial Effect of SNEMA SSB on Fall COG abs(Depth) " 
)
dev.off()

png(here("Figures/Raw_data_trends/GAMS/SNEMA/AMO_FallDepthDist.png"),
    width = 2600, height = 2000, units = "px", res = 300)
plot.gam(
  snema_distgammadepth1,
  xlab = "AMO",
  ylab = "Partial Effect",
  select = 2,
  cex.lab = 1.5,
  cex.axis = 1.4,
  rug = TRUE,
  shade = TRUE,
  col = "black",
  shade.col = "#E9E9E9",
  lwd = 2,
  main = "Partial Effect of AMO on Fall COG abs(Depth) " 
)
dev.off()

# Find minimum depth (should be negative)
min_depth <- min(snema_distribution_gamdata$fall_cog_depth, na.rm = TRUE)

# Add a constant to make all values positive (e.g., shift up by |min_depth| + small buffer)
snema_distribution_gamdata$fall_cog_depth_shifted <- snema_distribution_gamdata$fall_cog_depth + abs(min_depth) + 1

snema_distgammadepth_shifted <- gam(fall_cog_depth_shifted ~ s(SSB, k=5), 
                                    family = Gamma(), method = "REML", 
                                    data = snema_distribution_gamdata)
summary(snema_distgammadepth_shifted)

png(here("Figures/Raw_data_trends/GAMS/SNEMA/SSB_FallDepthDist_shifted.png"),
    width = 2600, height = 2000, units = "px", res = 300)
plot.gam(
  snema_distgammadepth_shifted,
  xlab = "SSB",
  ylab = "Partial Effect",
  select = 1,
  cex.lab = 1.5,
  cex.axis = 1.4,
  rug = TRUE,
  shade = TRUE,
  col = "black",
  shade.col = "#E9E9E9",
  lwd = 2,
  main = "Partial Effect of SSB on Fall COG Depth"
)
dev.off()

AIC(snema_distgammadepth_shifted)

