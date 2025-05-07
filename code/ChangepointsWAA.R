#changepoint analysis using EnvCpt package
library(EnvCpt)
library(changepoint)
library(here)
#to make plot figures:
library(png)
library(grid)
library(gridExtra)
library(plotrix) #to make cropped abline
library(dplyr) #renaming columns

#load in age data for each stock
gom_age <- read.csv(here("data/gom_age.csv"))
gbk_age <- read.csv(here("data/gbk_age.csv"))
snema_age <- read.csv(here("data/snema_age.csv"))

#load in waa data for each stock
#waa_gom <- read.csv(here("data/WAA_gom.csv"))
waa_gbk <- read.csv(here("data/WAA_gbk.csv"))
waa_snema <- read.csv(here("data/WAA_snema.csv"))
colnames(waa_snema)[colnames(waa_snema) == "Age7."] <- "Age7" #getting rid of random period in column

#NAA data filtering by AGE
gom_age1 <- gom_age[which(gom_age$AGE=='1'),]
gbk_age1 <- gbk_age[which(gbk_age$AGE=='1'),]
snema_age1 <- snema_age[which(snema_age$AGE=='1'),]
gom_age2 <- gom_age[which(gom_age$AGE=='2'),]
gbk_age2 <- gbk_age[which(gbk_age$AGE=='2'),]
snema_age2 <- snema_age[which(snema_age$AGE=='2'),]
gom_age3 <- gom_age[which(gom_age$AGE=='3'),]
gbk_age3 <- gbk_age[which(gbk_age$AGE=='3'),]
snema_age3 <- snema_age[which(snema_age$AGE=='3'),]
gom_age4 <- gom_age[which(gom_age$AGE=='4'),]
gbk_age4 <- gbk_age[which(gbk_age$AGE=='4'),]
snema_age4 <- snema_age[which(snema_age$AGE=='4'),]
gom_age5 <- gom_age[which(gom_age$AGE=='5'),]
gbk_age5 <- gbk_age[which(gbk_age$AGE=='5'),]
snema_age5 <- snema_age[which(snema_age$AGE=='5'),]
gom_age6 <- gom_age[which(gom_age$AGE=='6'),]
gbk_age6 <- gbk_age[which(gbk_age$AGE=='6'),]
snema_age6 <- snema_age[which(snema_age$AGE=='6'),]
gom_age7 <- gom_age[which(gom_age$AGE=='7'),]
gbk_age7 <- gbk_age[which(gbk_age$AGE=='7'),]
snema_age7 <- snema_age[which(snema_age$AGE=='7'),]

#waa data filtering by AGE and filtering out rows with NA values.
gbk_waa1 <- waa_gbk[!is.na(waa_gbk$Age1), c("Year", "Age1")]
gbk_waa2 <- waa_gbk[!is.na(waa_gbk$Age2), c("Year", "Age2")]
gbk_waa3 <- waa_gbk[!is.na(waa_gbk$Age3), c("Year", "Age3")] 
gbk_waa4 <- waa_gbk[!is.na(waa_gbk$Age4), c("Year", "Age4")]
gbk_waa5 <- waa_gbk[!is.na(waa_gbk$Age5), c("Year", "Age5")]
gbk_waa6 <- waa_gbk[!is.na(waa_gbk$Age6), c("Year", "Age6")]
gbk_waa7 <- waa_gbk[!is.na(waa_gbk$Age7), c("Year", "Age7")]
snema_waa1 <- waa_snema[!is.na(waa_snema$Age1), c("Year", "Age1")]
snema_waa2 <- waa_snema[!is.na(waa_snema$Age2), c("Year", "Age2")]
snema_waa3 <- waa_snema[!is.na(waa_snema$Age3), c("Year", "Age3")]
snema_waa4 <- waa_snema[!is.na(waa_snema$Age4), c("Year", "Age4")]
snema_waa5 <- waa_snema[!is.na(waa_snema$Age5), c("Year", "Age5")]
snema_waa6 <- waa_snema[!is.na(waa_snema$Age6), c("Year", "Age6")]
snema_waa7 <- waa_snema[!is.na(waa_snema$Age7.), c("Year", "Age7.")]

#### ENVCPT FUNCTION ####
trackyear<-data.frame(matrix(ncol = 1, nrow = 0))
colnames(trackyear) <- c('cpyear')
envcpt_fun<- function(data,dfnames,foldertype,yax,ylabel,datatype){
  for (k in 2:length(data)){
    clean_data <- na.omit(data[c(1:48), c(1, k)])  #keep only non-NA rows
    fit_envcpt = envcpt(na.omit(data[c(1:48),k]),models="meancpt",minseglen=5)  #fit mean model
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

####### call function ########
#perform full joins by year using Reduce and merge
waadf <- Reduce(function(x, y) merge(x, y, by = "Year", all = TRUE), 
                 list(waa_gbk, waa_snema))
waadf <- waadf[-c(1), ] #get rid of first row because there are NAs for GBK

#check result
colnames(waadf)<-c('year',
                    'age1GBK','age2GBK', 'age3GBK', 'age4GBK', 'age5GBK', 'age6GBK', 'age7GBK',
                   'age1SNEMA', 'age2SNEMA', 'age3SNEMA', 'age4SNEMA', 'age5SNEMA', 'age6SNEMA', 'age7SNEMA')
head(waadf)
waanames<-colnames(waadf[2:15])
#run function
envcpt_fun(waadf,waanames,foldertype="waa",yax="Weight at Age",ylabel="Weight",datatype="WAA Data")


