#changepoint analysis using EnvCpt package
library(EnvCpt)
library(changepoint)
library(here)
#to make plot figures:
library(png)
library(grid)
library(gridExtra)
library(plotrix) #to make cropped abline

#dismap data is species-wide, can I possibly separate by stock using lat/long?
#load spring and fall datasets from dismap
cogfall <- read.csv(here("data/distributionmetricsFall.csv"))
cogspring <- read.csv(here("data/distributionmetricsSpring.csv"))

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

#merge spring and fall datasets by year
cogdf <- Reduce(function(x, y) merge(x, y, by = "YEAR", all = TRUE), 
                 list(cogfall,cogspring))

#check result and rename columns
colnames(cogdf)<-c('year',
                    'fall_cog_lat','fall_cog_depth','fall_min_lat','fall_max_lat',
                    'spring_cog_lat','spring_cog_depth','spring_min_lat','spring_max_lat')

########DEPTH###########
cogdepthnames<-c('fall_cog_depth','spring_cog_depth')
envcpt_fun(cogdf[, c(1,3,7)],cogdepthnames,foldertype="cog",yax="COG",ylabel="Depth",datatype="COG Data")

#########LATITUDE#######
coglatnames<-c('fall_cog_lat','spring_cog_lat')
envcpt_fun(cogdf[, c(1,2,6)],coglatnames,foldertype="cog",yax="COG",ylabel="Latitude",datatype="COG Data")















