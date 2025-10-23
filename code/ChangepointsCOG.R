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



#######ORIGINAL Function that drops changepoint if it exceeds .25 sd from both segement means#####
trackyear <- data.frame(matrix(ncol = 1, nrow = 0))
colnames(trackyear) <- c('cpyear')

envcpt_fun <- function(data, dfnames, foldertype, yax, ylabel, datatype) {
  for (k in 2:length(data)) { #loop through each column (variable) in the data
    clean_data <- na.omit(data[1:46, c(1, k)])  # keep only non-NA rows, isolate time series
    data_column <- na.omit(data[1:46, k])
    threshold <- max(0.25 * sd(data_column, na.rm = TRUE), 0.2) #set threshold for valid mean shifts -> adapts to the natural variability in the data (0.25 × SD). Ensures the threshold doesn't get too small with low-variance series (minimum of 0.2)
    fit_envcpt <- envcpt(data_column, models = "meancpt", minseglen = 5,
                         penalty = "Manual", pen.value = 4 * log(length(data_column))) #detecting changepoints in mean only while applying manual penalty to control sensitivity and (4*log(n)) to reduce overfitting 
    ints <- unlist(fit_envcpt$meancpt@param.est[[1]])
    cp <- fit_envcpt$meancpt@cpts
    cp <- cp[!cp %in% c(1, length(data_column))]  # remove changepoints at endpoints
    
    # Remove changepoints caused by singular extreme datapoints
    if (length(cp) > 0) {
      new_cp <- c()
      for (i in seq_along(cp)) {
        idx <- cp[i]
        window <- data_column[(idx - 1):(idx + 1)]
        if (length(window) < 3) next
        pre_mean <- mean(data_column[1:idx], na.rm = TRUE)
        post_mean <- mean(data_column[(idx + 1):length(data_column)], na.rm = TRUE)
        sd_val <- sd(data_column, na.rm = TRUE)
        outlier_count <- sum(abs(window - pre_mean) > 0.5 * sd_val & abs(window - post_mean) > 0.5 * sd_val, na.rm = TRUE)
        if (outlier_count <= 1) next  # skip this changepoint
        new_cp <- c(new_cp, idx)
      } # chunk above is to filter out changepoints caused by singular extreme values that can be false positives. I define a small 3-point window around the changepoint and keep it if at least two of those points differ significantly (>0.5xSD) from both pre- and post-mean
      cp <- new_cp
      if (length(cp) == 0) {
        ints <- ints[1]
      } else {
        ints <- unlist(fit_envcpt$meancpt@param.est[[1]])[c(1, which(cp %in% fit_envcpt$meancpt@cpts) + 1)]
      } #chunk above recalculates segment means if changepoints were removed due to above filtering
    }
    
    # Filter changepoints by mean difference threshold
    if (length(ints) > 1) {
      mean_diffs <- abs(diff(ints))
      keep_idx <- which(mean_diffs >= threshold) #filter changepoints by meaningful segemnt different that lead to a real shift in mean, this filters out noise and tiny shifts 
      if (length(keep_idx) > 0) {
        cp <- cp[keep_idx]
        ints <- ints[c(1, keep_idx + 1)]
      } else {
        cp <- numeric(0)
        ints <- ints[1]
      }
    }
    
    cpyear <- na.omit(as.data.frame(data[1:46, c(1, k)]))[cp, 1] #Extract changepoint years
    trackyear2 <- as.data.frame(cpyear)
    trackyear <- rbind(trackyear, trackyear2)
    
    # Get first & last valid years
    startyear <- dplyr::first(clean_data[, 1])
    finalyear <- dplyr::last(clean_data[, 1])
    
    # Make plot
    png(here(paste0("figures/Raw_data_trends/Changepoint/", foldertype, "/EnvCPT/", dfnames[k - 1], ".png")),
        width = 449, height = 374.5, units = "px", res = 90)
    par(mar = c(4, 4.5, 2, 0.5))
    plot(data[, 1], data[, k], main = paste0(c(dfnames[k - 1], yax), collapse = " "),
         ylab = ylabel, xlab = "Year", type = "l", lwd = 3, col = "#00608A",
         cex.lab = 1.5, cex.axis = 1.2, xaxt = 'n')
    
    axis_labels <- unique(na.omit(c(startyear, cpyear, finalyear)))
    axis(side = 1, at = axis_labels, labels = axis_labels)
    
    if (length(cpyear) > 0) abline(v = cpyear, col = "red", lwd = 3, lty = 2)
    
    # Segmented means
    if (length(cp) == 0) {
      abline(h = ints[1], col = "red", lwd = 1, lty = 1)
    } else if (length(cp) == 1) {
      ablineclip(h = ints[1], col = "red", lwd = 1, lty = 1, x1 = startyear - 1, x2 = cpyear[1])
      ablineclip(h = ints[2], col = "red", lwd = 1, lty = 1, x1 = cpyear[1], x2 = finalyear)
    } else if (length(cp) == 2) {
      ablineclip(h = ints[1], col = "red", lwd = 1, lty = 1, x1 = startyear - 1, x2 = cpyear[1])
      ablineclip(h = ints[2], col = "red", lwd = 1, lty = 1, x1 = cpyear[1], x2 = cpyear[2])
      ablineclip(h = ints[3], col = "red", lwd = 1, lty = 1, x1 = cpyear[2], x2 = finalyear)
    } else if (length(cp) >= 3) {
      ablineclip(h = ints[1], col = "red", lwd = 1, lty = 1, x1 = startyear - 1, x2 = cpyear[1])
      ablineclip(h = ints[2], col = "red", lwd = 1, lty = 1, x1 = cpyear[1], x2 = cpyear[2])
      ablineclip(h = ints[3], col = "red", lwd = 1, lty = 1, x1 = cpyear[2], x2 = cpyear[3])
      ablineclip(h = ints[4], col = "red", lwd = 1, lty = 1, x1 = cpyear[3], x2 = finalyear)
    }
    
    legend("topright", inset = c(0.03, 0.03),
           legend = c(datatype, "Segmented Mean", "Changepoint"),
           col = c("#00608A", "red", "red"),
           lty = c(1, 1, 2), lwd = c(3, 1, 3), cex = 1.0)
    dev.off()
  }
  list2env(trackyear, envir = .GlobalEnv) #store changepoint years to global environment
}












