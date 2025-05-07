#changepoint analysis using EnvCpt package
library(EnvCpt)
library(changepoint)
library(here)
#to make plot figures:
library(png)
library(grid)
library(gridExtra)
library(plotrix) #to make cropped abline

#calling all age datasets for each stock
gom_age <- read.csv(here("data/gom_age.csv"))
gbk_age <- read.csv(here("data/gbk_age.csv"))
snema_age <- read.csv(here("data/snema_age.csv"))
#filtering by recruitment
gom_age1 <- gom_age[which(gom_age$AGE=='1'),]
gbk_age1 <- gbk_age[which(gbk_age$AGE=='1'),]
snema_age1 <- snema_age[which(snema_age$AGE=='1'),]
#filtering recruits by season
gom1spring <- gom_age1[which(gom_age1$SEASON=='SPRING'),]
gom1fall <- gom_age1[which(gom_age1$SEASON=='FALL'),]
gbk1spring <- gbk_age1[which(gbk_age1$SEASON=='SPRING'),]
gbk1fall <- gbk_age1[which(gbk_age1$SEASON=='FALL'),]
snema1spring <- snema_age1[which(snema_age1$SEASON=='SPRING'),]
snema1fall <- snema_age1[which(snema_age1$SEASON=='FALL'),]
#filtering recruits by season by survey
gom1springNMFS <- gom1spring[which(gom1spring$SURVEY=='NMFS spring BTS'),]
gom1springMADMF <- gom1spring[which(gom1spring$SEASON=='MADMF spring BTS'),]
gom1fallNMFS <- gom1fall[which(gom1fall$SURVEY=='NMFS fall BTS'),]
gom1fallMADMF <- gom1fall[which(gom1fall$SURVEY=='MADMF fall BTS'),]
gbk1springNMFS <- gbk1spring
gbk1fallNMFS <- gbk1fall
snema1springNMFS <- snema1spring[which(snema1spring$SURVEY=='NMFS spring BTS'),]
snema1springMADMF <- snema1spring[which(snema1spring$SEASON=='MADMF spring BTS'),]
snema1fallNMFS <- snema1fall[which(snema1fall$SURVEY=='NMFS fall BTS'),]
snema1fallMADMF <- snema1fall[which(snema1fall$SURVEY=='MADMF fall BTS'),]


#### ENVCPT FUNCTION ####
trackyear<-data.frame(matrix(ncol = 1, nrow = 0))
colnames(trackyear) <- c('cpyear')
envcpt_fun<- function(data,dfnames,foldertype,yax,ylabel,datatype){
  for (k in 2:length(data)){
    clean_data <- na.omit(data[c(1:46), c(1, k)])  #keep only non-NA rows
    fit_envcpt = envcpt(na.omit(data[c(1:46),k]),models="meancpt",minseglen=5)  #fit mean model
    ints = unlist((fit_envcpt$meancpt@param.est[1])) #get mean cpt values
    cp = fit_envcpt$meancpt@cpts #get changepoint year numbers
    cp = cp[!cp %in% c(1,nrow(na.omit(data[k])))]
    cpyear<-na.omit(as.data.frame(data[c(1,k)]))[cp,1] #get changepoint year
    trackyear2<-as.data.frame(cpyear)
    trackyear<- rbind(trackyear,trackyear2)
    
    #get first & last valid years
    startyear <- dplyr::first(clean_data[, 1])
    finalyear <- dplyr::last(clean_data[, 1])
    
    #make plot
    png(here(paste0("figures/Raw_data_trends/Changepoint/",foldertype,"/EnvCPT/",dfnames[k-1],".png")),width = 449, height = 374.5, units = "px",res=90)
    par(mar=c(4,4.5,2,0.5))
    plot(data[,1],data[,k],main=paste0(c(dfnames[k-1],yax),collapse=" "),ylab=ylabel, xlab="Year",type="l",lwd=3,col="#00608A",cex.lab=1.5,cex.axis=1.2,xaxt='n')
    axis(side=1,at=c(startyear,cpyear[1],cpyear[2],finalyear),labels=c(startyear,cpyear[1],cpyear[2],finalyear))
    abline(v=c(cpyear[1],cpyear[2],cpyear[3]),col="red",lwd=3,lty=2)
    if (length(cp)==0){abline(h=ints[1],col="red",lwd=1,lty=1)}else{
      if (length(cp)==1){ablineclip(h=ints[1],col="red",lwd=1,lty=1,x1=startyear-1,x2=cpyear[1])
        ablineclip(h=ints[2],col="red",lwd=1,lty=1,x1=cpyear[1],x2=finalyear)}else{
          if (length(cp)==2){ablineclip(h=ints[1],col="red",lwd=1,lty=1,x1=startyear-1,x2=cpyear[1])
            ablineclip(h=ints[2],col="red",lwd=1,lty=1,x1=cpyear[1],x2=cpyear[2])
            ablineclip(h=ints[3],col="red",lwd=1,lty=1,x1=cpyear[2],x2=finalyear)}else{
              if (length(cp)>=3){ablineclip(h=ints[1],col="red",lwd=1,lty=1,x1=startyear-1,x2=cpyear[1])
                ablineclip(h=ints[2],col="red",lwd=1,lty=1,x1=cpyear[1],x2=cpyear[2])
                ablineclip(h=ints[3],col="red",lwd=1,lty=1,x1=cpyear[2],x2=cpyear[3])
                ablineclip(h=ints[4],col="red",lwd=1,lty=1,x1=cpyear[3],x2=finalyear)}
              
            }}}
    legend("topright",inset=c(0.03,0.03), legend=c(datatype, "Segmented Mean","Changepoint"), col=c("#00608A", "red", "red"), lty=c(1,1,2),lwd=c(3,1,3), cex=1.0)
    dev.off()
  }
  list2env(trackyear, envir = .GlobalEnv)
}

####### call function for Age 1 recruitment (NMFS survey) ########
#extract data for each year and create individual data frames
gom_fall <- data.frame(year = gom1fallNMFS[1:42, 7], gom_fall_age1_nmfs = gom1fallNMFS[1:42, 9])
gom_spring <- data.frame(year = gom1springNMFS[1:43, 7], gom_spring_age1_nmfs = gom1springNMFS[1:43, 9])
gbk_fall <- data.frame(year = gbk1fallNMFS[1:42, 7], gbk_fall_age1_nmfs = gbk1fallNMFS[1:42, 9])
gbk_spring <- data.frame(year = gbk1springNMFS[1:36, 7], gbk_spring_age1_nmfs = gbk1springNMFS[1:36, 9])
snema_fall <- data.frame(year = snema1fallNMFS[1:43, 7], snema_fall_age1_nmfs = snema1fallNMFS[1:43, 9])
snema_spring <- data.frame(year = snema1springNMFS[1:44, 7], snema_spring_age1_nmfs = snema1springNMFS[1:44, 9])

#perform full joins by year using Reduce and merge
age1df <- Reduce(function(x, y) merge(x, y, by = "year", all = TRUE), 
                 list(gom_fall, gom_spring, gbk_fall, gbk_spring, snema_fall, snema_spring))

#check result
colnames(age1df)<-c('year',
                    'gom_fall_age1_nmfs','gom_spring_age1_nmfs',
                    'gbk_fall_age1_nmfs','gbk_spring_age1_nmfs',
                    'snema_fall_age1_nmfs','snema_spring_age1_nmfs')
head(age1df)
age1names<-colnames(age1df[2:7])
#run function
envcpt_fun(age1df,age1names,foldertype="recruitment",yax="Abundance",ylabel="Recruitment",datatype="Recruitment Data")
R_Envyears<-as.data.frame(cpyear)


##### Changepoint package #####
####Changepoint Function for Changepoint package####
trackyear<-data.frame(matrix(ncol = 1, nrow = 0))
colnames(trackyear) <- c('cpyear')
changepoint_fun<- function(data,dfnames,foldertype,yax,ylabel,datatype){
  for (k in 2:length(data)){
    fit_cpt = changepoint::cpt.mean(na.omit(data[c(1:46),k]),method = "AMOC",penalty="AIC",minseglen = 5)  # Fit mean model
    ints = param.est(fit_cpt)$mean #get mean cpt values
    cp = cpts(fit_cpt) #get changepoint year numbers
    cpyear<-na.omit(as.data.frame(data[c(1,k)]))[cp,1] #get changepoint year
    trackyear2<-as.data.frame(cpyear)
    trackyear<- rbind(trackyear,trackyear2)
    
    startyear<-dplyr::first(na.omit(data[data[,k]==(dplyr::first(na.omit(data[c(1:46),k]))),1])) #get first non NA year of data for every column
    finalyear<-dplyr::last(na.omit(data[data[,k]==(dplyr::last(na.omit(data[c(1:46),k]))),1])) #get final year of data for every column
    #make plot
    png(here(paste0("Figures/Raw_data_trends/Changepoint/",foldertype,"/changepoint/",dfnames[k-1],".png")),width = 449, height = 374.5, units = "px",res=90)
    par(mar=c(4,4.5,2,0.5))
    plot(data[,1],data[,k],main=paste0(c(dfnames[k-1],yax),collapse=" "),ylab=ylabel, xlab="Year",type="l",lwd=3,col="#00608A",cex.lab=1.5,cex.axis=1.2,xaxt='n')
    
    abline(v=c(cpyear[1],cpyear[2],cpyear[3]),col="red",lwd=3,lty=2)
    if (length(cp)==0){abline(h=ints[1],col="red",lwd=1,lty=1)}else{
      if (length(cp)==1){ablineclip(h=ints[1],col="red",lwd=1,lty=1,x1=startyear-1,x2=cpyear[1])
        ablineclip(h=ints[2],col="red",lwd=1,lty=1,x1=cpyear[1],x2=finalyear)}else{
          if (length(cp)==2){ablineclip(h=ints[1],col="red",lwd=1,lty=1,x1=startyear-1,x2=cpyear[1])
            ablineclip(h=ints[2],col="red",lwd=1,lty=1,x1=cpyear[1],x2=cpyear[2])
            ablineclip(h=ints[3],col="red",lwd=1,lty=1,x1=cpyear[2],x2=finalyear)}else{
              if (length(cp)>=3){ablineclip(h=ints[1],col="red",lwd=1,lty=1,x1=startyear-1,x2=cpyear[1])
                ablineclip(h=ints[2],col="red",lwd=1,lty=1,x1=cpyear[1],x2=cpyear[2])
                ablineclip(h=ints[3],col="red",lwd=1,lty=1,x1=cpyear[2],x2=cpyear[3])
                ablineclip(h=ints[4],col="red",lwd=1,lty=1,x1=cpyear[3],x2=finalyear)}
              
            }}}
    
    axis(side=1,at=c(startyear,cpyear,finalyear),labels=c(startyear,cpyear,finalyear))
    #    abline(v=cpyear,col="red",lwd=3,lty=2)
    #    ablineclip(h=ints[1],col="red",lwd=1,lty=1,x1=startyear-1,x2=if(length(cp)==0){finalyear}else{cpyear})
    #    ablineclip(h=ints[2],col="red",lwd=1,lty=1,x1=if(length(cp)==0){startyear-1}else{cpyear},x2=finalyear)
    legend("topright",inset=c(0.03,0.03), legend=c(datatype, "Segmented Mean","Changepoint"), col=c("#00608A", "red", "red"), lty=c(1,1,2),lwd=c(3,1,3), cex=1.0)
    dev.off()
  }
  list2env(trackyear, envir = .GlobalEnv)
}
#####################change point analysis for Age 1 recruitment data for winter flounder####
nm <- list.files(path =here("data/Final_Data_for_Modelers/Recruitment"), pattern = ".csv", full.names = TRUE)
nm2 <- list.files(path =here("data/Final_Data_for_Modelers/Recruitment"), pattern = ".csv", full.names =FALSE)
list2env(lapply(setNames(nm, make.names(gsub("*.csv$", "",nm2))),read.csv),envir=.GlobalEnv)
rm(nm,nm2)
#### Call function for age1 data####
changepoint_fun(age1df[-c(8:9)],age1names,foldertype="recruitment",yax="Abundance",ylabel="Recruitment",datatype="Recruitment Data")
R_cpyears<-as.data.frame(cpyear)


###########R/SSB Calculations for SNE/MA##############
#Merge recruitment dataset with ssb dataset 

#Check heads
print(head(ssb_snema))
print(head(age1df))

#change name of year column to match
names(ssb_snema)[names(ssb_snema) == 'YEAR'] <- 'year'

#select only the columns I want out of the recrutiment df 
age1snema <- age1df %>%
  dplyr::select(year, snema_fall_age1_nmfs, snema_spring_age1_nmfs)

#merge to create new dataset!
recruitment_snema <- merge(ssb_snema, age1snema, by = "year")

#create R/SSB!
recruitment_snema$rssb_spring<-recruitment_snema$snema_fall_age1_nmfs/lag(recruitment_snema[,"SSB"])
recruitment_snema$rssb_fall<-recruitment_snema$snema_spring_age1_nmfs/lag(recruitment_snema[,"SSB"])

###########R/SSB Calculations for GBK##############
#Merge recruitment dataset with ssb dataset 

#Check heads
print(head(ssb_gbk))
print(head(age1df))

#change name of year column to match
names(ssb_gbk)[names(ssb_gbk) == 'YEAR'] <- 'year'

#select only the columns I want out of the recrutiment df 
age1gbk <- age1df %>%
  dplyr::select(year, gbk_fall_age1_nmfs, gbk_spring_age1_nmfs)

#merge to create new dataset!
recruitment_gbk <- merge(ssb_gbk, age1gbk, by = "year")

#create R/SSB!
recruitment_gbk$rssb_spring<-recruitment_gbk$gbk_fall_age1_nmfs/lag(recruitment_gbk[,"SSB"])
recruitment_gbk$rssb_fall<-recruitment_gbk$gbk_spring_age1_nmfs/lag(recruitment_gbk[,"SSB"])

#########NOW COMBINE ALL#########
#make dataframe with everything combined (sans GOM because no WAA table right now)

rssbdf<-data.frame(recruitment_snema[c(2:40),1],
                   recruitment_snema[c(2:40),5],recruitment_snema[c(2:40),6],
                   recruitment_gbk[c(2:40),5],recruitment_gbk[c(2:40),6])
colnames(rssbdf)<-c('Year',
                    'SNEMA_Fall_rssb','SNEMA_Spring_rssb',
                    'GBK_Fall_rssb','GBK_Spring_rssb')

rssbnames<-colnames(rssbdf[2:4])
envcpt_fun(rssbdf[-c(8:9)],rssbnames,foldertype="recruitment",yax="",ylabel="Winter Flounder R/SSB",datatype="Recruitment Data")
RSSB_Envyears<-as.data.frame(cpyear)

########CHATGPT##########

# Create the combined dataset with properly aligned rows
rssbdf <- data.frame(
  recruitment_snema[1:40, 1],  
  recruitment_snema[1:40, 5],  
  recruitment_snema[1:40, 6],  
  recruitment_gbk[1:40, 5],  
  recruitment_gbk[1:40, 6]
)

# Assign correct column names
colnames(rssbdf) <- c('Year', 'SNEMA_Fall_rssb', 'SNEMA_Spring_rssb', 'GBK_Fall_rssb', 'GBK_Spring_rssb')

# Ensure correct indexing of column names
rssbnames <- colnames(rssbdf)[2:5]

# Run the function
envcpt_fun(rssbdf, rssbnames, foldertype="recruitment", yax="", ylabel="Winter Flounder R/SSB", datatype="Recruitment Data")

# Store results
RSSB_Envyears <- as.data.frame(cpyear)


#### ENVCPT FUNCTION ####
trackyear <- data.frame(matrix(ncol = 1, nrow = 0))
colnames(trackyear) <- c('cpyear')

envcpt_fun <- function(data, dfnames, foldertype, yax, ylabel, datatype) {
  for (k in 2:ncol(data)) {
    
    # Ensure we are not removing key data with na.omit()
    fit_envcpt <- envcpt(data[c(1:46), k], models="meancpt", minseglen=5)
    
    # Extract mean change points
    ints <- unlist((fit_envcpt$meancpt@param.est[1])) 
    cp <- fit_envcpt$meancpt@cpts  
    cp <- cp[!cp %in% c(1, nrow(na.omit(data[k])))]  # Remove first and last row if selected
    
    # Ensure cpyear is not empty
    cpyear <- if (length(cp) > 0) na.omit(as.data.frame(data[c(1, k)]))[cp, 1] else NA
    trackyear2 <- as.data.frame(cpyear)
    trackyear <- rbind(trackyear, trackyear2)
    
    # Get start and end years
    startyear <- dplyr::first(na.omit(data[data[, k] == dplyr::first(na.omit(data[c(1:46), k])), 1]))
    finalyear <- dplyr::last(na.omit(data[data[, k] == dplyr::last(na.omit(data[c(1:46), k])), 1]))
    
    # Create and save the plot
    png(here(paste0("figures/Raw_data_trends/Changepoint/", foldertype, "/EnvCPT/", dfnames[k - 1], ".png")),
        width = 449, height = 374.5, units = "px", res = 90)
    par(mar = c(4, 4.5, 2, 0.5))
    plot(data[, 1], data[, k], main = paste0(dfnames[k - 1], yax), 
         ylab = ylabel, xlab = "Year", type = "l", lwd = 3, col = "#00608A", 
         cex.lab = 1.5, cex.axis = 1.2, xaxt = 'n')
    
    # Ensure cpyear is valid before plotting
    if (!is.na(cpyear[1])) {
      axis(side = 1, at = c(startyear, cpyear, finalyear), labels = c(startyear, cpyear, finalyear))
      abline(v = cpyear, col = "red", lwd = 3, lty = 2)
    }
    
    # Plot segmented mean lines safely
    if (length(cp) == 0) {
      abline(h = ints[1], col = "red", lwd = 1, lty = 1)
    } else {
      ablineclip(h = ints[1], col = "red", lwd = 1, lty = 1, x1 = startyear - 1, x2 = cpyear[1])
      if (length(cp) >= 1) ablineclip(h = ints[2], col = "red", lwd = 1, lty = 1, x1 = cpyear[1], x2 = finalyear)
      if (length(cp) >= 2) ablineclip(h = ints[3], col = "red", lwd = 1, lty = 1, x1 = cpyear[2], x2 = finalyear)
      if (length(cp) >= 3) ablineclip(h = ints[4], col = "red", lwd = 1, lty = 1, x1 = cpyear[3], x2 = finalyear)
    }
    
    # Add legend
    legend("topright", inset = c(0.03, 0.03), 
           legend = c(datatype, "Segmented Mean", "Changepoint"), 
           col = c("#00608A", "red", "red"), 
           lty = c(1, 1, 2), lwd = c(3, 1, 3), cex = 1.0)
    
    dev.off()
  }
  list2env(trackyear, envir = .GlobalEnv)
}

####Changepoint Function for Changepoint package####
trackyear<-data.frame(matrix(ncol = 1, nrow = 0))
colnames(trackyear) <- c('cpyear')
changepoint_fun<- function(data,dfnames,foldertype,yax,ylabel,datatype){
  for (k in 2:length(data)){
    fit_cpt = changepoint::cpt.mean(na.omit(data[c(1:46),k]),method = "AMOC",penalty="AIC",minseglen = 5)  # Fit mean model
    ints = param.est(fit_cpt)$mean #get mean cpt values
    cp = cpts(fit_cpt) #get changepoint year numbers
    cpyear<-na.omit(as.data.frame(data[c(1,k)]))[cp,1] #get changepoint year
    trackyear2<-as.data.frame(cpyear)
    trackyear<- rbind(trackyear,trackyear2)
    
    startyear<-dplyr::first(na.omit(data[data[,k]==(dplyr::first(na.omit(data[c(1:46),k]))),1])) #get first non NA year of data for every column
    finalyear<-dplyr::last(na.omit(data[data[,k]==(dplyr::last(na.omit(data[c(1:46),k]))),1])) #get final year of data for every column
    #make plot
    png(here(paste0("Figures/Raw_data_trends/Changepoint/",foldertype,"/changepoint/",dfnames[k-1],".png")),width = 449, height = 374.5, units = "px",res=90)
    par(mar=c(4,4.5,2,0.5))
    plot(data[,1],data[,k],main=paste0(c(dfnames[k-1],yax),collapse=" "),ylab=ylabel, xlab="Year",type="l",lwd=3,col="#00608A",cex.lab=1.5,cex.axis=1.2,xaxt='n')
    
    abline(v=c(cpyear[1],cpyear[2],cpyear[3]),col="red",lwd=3,lty=2)
    if (length(cp)==0){abline(h=ints[1],col="red",lwd=1,lty=1)}else{
      if (length(cp)==1){ablineclip(h=ints[1],col="red",lwd=1,lty=1,x1=startyear-1,x2=cpyear[1])
        ablineclip(h=ints[2],col="red",lwd=1,lty=1,x1=cpyear[1],x2=finalyear)}else{
          if (length(cp)==2){ablineclip(h=ints[1],col="red",lwd=1,lty=1,x1=startyear-1,x2=cpyear[1])
            ablineclip(h=ints[2],col="red",lwd=1,lty=1,x1=cpyear[1],x2=cpyear[2])
            ablineclip(h=ints[3],col="red",lwd=1,lty=1,x1=cpyear[2],x2=finalyear)}else{
              if (length(cp)>=3){ablineclip(h=ints[1],col="red",lwd=1,lty=1,x1=startyear-1,x2=cpyear[1])
                ablineclip(h=ints[2],col="red",lwd=1,lty=1,x1=cpyear[1],x2=cpyear[2])
                ablineclip(h=ints[3],col="red",lwd=1,lty=1,x1=cpyear[2],x2=cpyear[3])
                ablineclip(h=ints[4],col="red",lwd=1,lty=1,x1=cpyear[3],x2=finalyear)}
              
            }}}
    
    axis(side=1,at=c(startyear,cpyear,finalyear),labels=c(startyear,cpyear,finalyear))
    #    abline(v=cpyear,col="red",lwd=3,lty=2)
    #    ablineclip(h=ints[1],col="red",lwd=1,lty=1,x1=startyear-1,x2=if(length(cp)==0){finalyear}else{cpyear})
    #    ablineclip(h=ints[2],col="red",lwd=1,lty=1,x1=if(length(cp)==0){startyear-1}else{cpyear},x2=finalyear)
    legend("topright",inset=c(0.03,0.03), legend=c(datatype, "Segmented Mean","Changepoint"), col=c("#00608A", "red", "red"), lty=c(1,1,2),lwd=c(3,1,3), cex=1.0)
    dev.off()
  }
  list2env(trackyear, envir = .GlobalEnv)
}
#####################change point analysis for r/ssb data for winter flounder#######################
#### Call function for age1 data####
changepoint_fun(age1df[-c(8:9)],age1names,foldertype="recruitment",yax="Abundance",ylabel="Recruitment",datatype="Recruitment Data")
envcpt_fun(rssbdf, rssbnames, foldertype="recruitment", yax="", ylabel="Winter Flounder R/SSB", datatype="Recruitment Data")






