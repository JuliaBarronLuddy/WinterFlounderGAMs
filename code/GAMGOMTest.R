library(ecodata)
library(dplyr)
library(here)
library(tidyverse)
library(mgcv)
################################################################################
##Data Exploration

#Data load
gom_age <- read.csv(here("data/gom_age.csv"))
gom_age1 <- gom_age[which(gom_age$AGE=='1'),]
gom_NMFS_age1 <- gom_age1[which(gom_age1$SURVEY=='NMFS spring BTS'),]
gom_NMFS_age1 <- gom_NMFS_age1 %>% select(c('YEAR', 'NO_AT_AGE'))
names(gom_NMFS_age1)[names(gom_NMFS_age1) == 'YEAR'] <- 'Year'

#Bottom temperature
bottomtemp <- ecodata::bottom_temp_comp
bottomtemp <- bottomtemp[which(bottomtemp$EPU=='GOM'),]
bottomtemp <- bottomtemp[which(bottomtemp$Var=='Spring_Bottom Temp Anomaly'),]
bottomtemp <- bottomtemp[-which(bottomtemp$Source=='PSY'),]
bottomtemp <- subset(bottomtemp, select = -c(EPU, Source, Var, Units))
names(bottomtemp)[names(bottomtemp) == 'Value'] <- 'BTAnom'
names(bottomtemp)[names(bottomtemp) == 'Time'] <- 'Year'

#Surface temperature
surfacetemp <- ecodata::seasonal_oisst_anom
surfacetemp <- surfacetemp[which(surfacetemp$EPU=='GOM'),]
surfacetemp <- surfacetemp[which(surfacetemp$Var=='Spring'),]
names(surfacetemp)[names(surfacetemp) == 'Value'] <- 'SSTAnom'
names(surfacetemp)[names(surfacetemp) == 'Time'] <- 'Year'
surfacetemp <- subset(surfacetemp, select = -c(EPU, Var))

#GSI
gsi <- ecodata::gsi
gsi <- gsi[which(gsi$Var=='gulf stream index'),]
gsi <- gsi %>% 
  separate(Time, into = c("Year", "Month"), sep = "\\.", convert = TRUE)
gsiyear <- gsi %>%
  group_by(Year) %>%
  summarize(Mean_Value = mean(Value))
names(gsiyear)[names(gsiyear) == 'Mean_Value'] <- 'GSI'

#Combining 3 ecovs
ecovars <- bottomtemp %>%
  full_join(surfacetemp, by = "Year") %>%
  full_join(gsiyear, by = "Year")

ecovars <- ecovars[order(ecovars$Year),]
ecovars <- subset(ecovars, Year > 1976)

head(ecovars)
str(ecovars)
summary(ecovars)

plot(GSI ~Year, data=ecovars, main="Exploratory Visual", xlab="Year", ylab="Degrees Lat Change & Temp", type="l", lwd=3,col="#00608A", cex.lab=1.4,cex.axis=1.1)
lines(SSTAnom ~Year, data=ecovars, xlab="Year", type="l",col="#EA4F12",lwd=3)
lines(BTAnom ~Year, data=ecovars, xlab="Year", type="l",col="#00736D",lwd=3)


#quick visual checks for outliers in data
dotchart(ecovars$GSI, main = "GSI")
dotchart(ecovars$SSTAnom, main = "SST Anomalies")
dotchart(ecovars$BTAnom, main = "BT Anomalies")
boxplot(ecovars$GSI, main = "GSI")
boxplot(ecovars$SSTAnom, main = "SST Anomalies")
boxplot(ecovars$BTAnom, main = "BT Anomalies")
boxplot(gom_NMFS_age1$NO_AT_AGE, main = "No. at age 1")
#check distribution
hist(ecovars$GSI, main = "GSI") 
hist(ecovars$SSTAnom, main = "SST Anomalies") 
hist(ecovars$BTAnom, main = "BT Anomalies") 
hist(gom_NMFS_age1$NO_AT_AGE, main = "No. at age 1")

#Shapiro Test with Jamie's function
shapiro_fun <- function(data) {
  for (i in seq_along(data)) {
    p_value <- shapiro.test(data[[i]])$p.value
    cat(names(data)[i], ": ", if (p_value > 0.05) "Normal" else "Not Normal", "\n")
  }
}

shapiro_fun(ecovars) #only BTAnom is not normal, rest are normal
shapiro.test(ecovars$GSI) 
shapiro.test(ecovars$SSTAnom)
shapiro.test(ecovars$BTAnom)

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
Mypairs(ecovars) 
corvif(ecovars) 
corvif(ecovars[c(2:4)])
corvif(ecovars[c(3:4)])
corvif(ecovars[c(1,4)])

########################################################################

#Run GAM

gam1data <- gom_NMFS_age1 %>%
  full_join(ecovars, by = "Year")

gam1data <- gam1data[order(gam1data$Year),]

gam1 <- gam(NO_AT_AGE~s(BTAnom, k=5)+s(SSTAnom, k=5), family=tw(), method="REML", data=gam1data)
#REML: Restricted maximum likelihood approach to smoothing
summary(gam1)
gam.check(gam1)
concurvity(gam1)
AIC(gam1) #-46 gamma, -43 tw

plot.gam(gam1,xlab="BT Anomalies",ylab="Partial Effect of BT Anomalies on Age 1 Recruitment", select=1, cex.lab=1.5,cex.axis=1.4,rug=TRUE,shade = TRUE,col = "Black",shade.col="#E9E9E9",lwd = 2)
#rug(gam1$BTAnom, ticksize=0.1, side=1, lwd=1.75,col="Black")
plot.gam(gam1,xlab="SST Anomalies",ylab="Partial Effect of SST Anomalies on Age 1 Recruitment", select=2, cex.lab=1.5,cex.axis=1.4,rug=TRUE,shade = TRUE,col = "Black",shade.col="#E9E9E9",lwd = 2)
plot.gam(gam1,xlab="GSI",ylab="Partial Effect of GSI on Age 1 Recruitment", select=3, cex.lab=1.5,cex.axis=1.4,rug=TRUE,shade = TRUE,col = "Black",shade.col="#E9E9E9",lwd = 2)










