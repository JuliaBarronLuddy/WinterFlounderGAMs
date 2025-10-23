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

library(ecodata)
condition <- ecodata::condition #we want Var -- Winter flounder
condition <- condition[which(condition$Var=='Winter flounder'),]
colnames(condition)[colnames(condition) == "Time"] <- "Year"

cond_gom <- condition[which(condition$EPU=='GOM'),]
cond_gbk <- condition[which(condition$EPU=='GB'),]
cond_sne <- condition[which(condition$EPU=='MAB'),]

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

#perform full joins by year using Reduce and merge
conddf <- Reduce(function(x, y) merge(x, y, by = "Year", all = TRUE), 
                list(cond_gom, cond_gbk, cond_sne))
conddf <- conddf[-c(26), ] #get rid of 2017 because SSE didnt have a value for that year

#get rid of columns that aren't needed for analysis
conddf <- conddf[, -c(2, 3, 5, 6, 7, 9, 10, 11, 13)]



#check result and rename columns
colnames(conddf)<-c('year',
                   'condition_GOM','condition_GBK', 'condition_MAB')
condnames <- colnames(conddf[2:4])

#run function
envcpt_fun(conddf,condnames,foldertype="condition",yax="Relative Condition",ylabel="Relative Condition",datatype="EcoData Condition")
COND_Envyears <- as.data.frame(cpyear)











