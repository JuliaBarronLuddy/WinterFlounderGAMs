### Estimating cod SSB outside stock assessment models
### Approach is to estimate SSB based on abundance of age 4+ fish
library(here)

### will use numbers at age data from trawl survey to do this
Cod_NAA<-read.csv(here("data/cod_NAA.csv"))
### WAA data from Charles
WAA_EGOM<-read.csv(here("data/WAA/EGOM_mean_weight_at_age.csv"))
WAA_WGOM<-read.csv(here("data/WAA/WGOM_mean_weight_at_age.csv"))
WAA_GBK<-read.csv(here("data/WAA/GBK_mean_weight_at_age.csv"))
WAA_SNEMA<-read.csv(here("data/WAA/SNEMA_mean_weight_at_age.csv"))


##### Stock specific SSB estimates#######
#NAA data wrangling
Cod_NAA<-Cod_NAA[,c(2:3,10:15,17:18)]
Cod_NAA = Cod_NAA[Cod_NAA$SURVEY == "NEFSC_BTS",]
Cod_NAA = Cod_NAA[!Cod_NAA$YEAR  < 1982,]
Cod_NAA = Cod_NAA[!Cod_NAA$YEAR  > 2021,]
NAA_EGOM = Cod_NAA[Cod_NAA$STOCK == "EGOM",]
NAA_WGOM = Cod_NAA[Cod_NAA$STOCK == "WGOM",]
NAA_GBK = Cod_NAA[Cod_NAA$STOCK == "GBK",]
NAA_SNE = Cod_NAA[Cod_NAA$STOCK == "SNE",]

#WAA data wrnagling
WAA_EGOM<-WAA_EGOM[,c(6:8,10)]
WAA_WGOM<-WAA_WGOM[,c(6:8,10)]
WAA_GBK<-WAA_GBK[,c(6:8,10)]
WAA_SNEMA<-WAA_SNEMA[,c(6:8,10)]

###clip to years I want (1982-2019)###
WAA_EGOM = WAA_EGOM[!WAA_EGOM$YEAR  < 1982,]
WAA_EGOM = WAA_EGOM[!WAA_EGOM$AGE  < 4,]
WAA_EGOM = WAA_EGOM[!WAA_EGOM$AGE  > 9,]

WAA_WGOM = WAA_WGOM[!WAA_WGOM$YEAR  < 1982,]
WAA_WGOM = WAA_WGOM[!WAA_WGOM$AGE  < 4,]
WAA_WGOM = WAA_WGOM[!WAA_WGOM$AGE  > 9,]

WAA_GBK = WAA_GBK[!WAA_GBK$YEAR  < 1982,]
WAA_GBK = WAA_GBK[!WAA_GBK$AGE  < 4,]
WAA_GBK = WAA_GBK[!WAA_GBK$AGE  > 9,]

WAA_SNEMA = WAA_SNEMA[!WAA_SNEMA$YEAR  < 1982,]
WAA_SNEMA = WAA_SNEMA[!WAA_SNEMA$AGE  < 4,]
WAA_SNEMA = WAA_EGOM[!WAA_SNEMA$AGE  > 9,]

## Create SSB dataframe
### GET SSB ESTIMATES FOR EACH STOCK #####
#function
get_ssb<- function(NAA_data, WAA_data,SSB_stock){
  
  SSB_stock<-NAA_data[,c(1:9)]
  names(SSB_stock)<-c("YEAR", "SEASON","Age4NAA","Age5NAA","Age6NAA","Age7NAA","Age8NAA","Age9NAA","STOCK")
  SSB_stock<-merge(SSB_stock,WAA_data,by=c("YEAR","SEASON"),all=TRUE)
  SSB_stock$SSB<-((SSB_stock$Age4NAA*SSB_stock$MEAN*(SSB_stock$AGE == 4))+
                    (SSB_stock$Age5NAA*SSB_stock$MEAN*(SSB_stock$AGE == 5))+
                    (SSB_stock$Age6NAA*SSB_stock$MEAN*(SSB_stock$AGE == 6))+
                   (SSB_stock$Age7NAA*SSB_stock$MEAN*(SSB_stock$AGE == 7))+
                   (SSB_stock$Age8NAA*SSB_stock$MEAN*(SSB_stock$AGE == 8))+
                   (SSB_stock$Age9NAA*SSB_stock$MEAN*(SSB_stock$AGE == 9))) #9+ in Jamie's heart
  SSB_stock<-SSB_stock[,c(1,2,12)]
  SSB_stock<-aggregate(SSB~YEAR+SEASON,SSB_stock,FUN=sum)
  
}
 
SSB_WGOM<-NAA_WGOM[,c(1:9)]
names(SSB_WGOM)<-c("YEAR", "SEASON","Age4NAA","Age5NAA","Age6NAA","Age7NAA","Age8NAA","Age9NAA","STOCK")
SSB_WGOM<-merge(SSB_WGOM,WAA_WGOM,by=c("YEAR","SEASON"),all=TRUE)
SSB_WGOM$SSB<-((SSB_WGOM$Age4NAA*SSB_WGOM$MEAN*(SSB_WGOM$AGE == 4))+
                 (SSB_WGOM$Age5NAA*SSB_WGOM$MEAN*(SSB_WGOM$AGE == 5))+
                 (SSB_WGOM$Age6NAA*SSB_WGOM$MEAN*(SSB_WGOM$AGE == 6))+
                  (SSB_WGOM$Age7NAA*SSB_WGOM$MEAN*(SSB_WGOM$AGE == 7))+
                  (SSB_WGOM$Age8NAA*SSB_WGOM$MEAN*(SSB_WGOM$AGE == 8))+
                  (SSB_WGOM$Age9NAA*SSB_WGOM$MEAN*(SSB_WGOM$AGE == 9)))
SSB_WGOM<-SSB_WGOM[,c(1,2,12)]
SSB_WGOM<-aggregate(SSB~YEAR+SEASON,SSB_WGOM,FUN=sum)


SSB_WGOM<-get_ssb(NAA_WGOM,WAA_WGOM,SSB_WGOM)
SSB_EGOM<-get_ssb(NAA_EGOM,WAA_EGOM,SSB_EGOM)
SSB_GBK<-get_ssb(NAA_GBK,WAA_GBK,SSB_GBK)
SSB_SNE<-get_ssb(NAA_SNE,WAA_SNEMA,SSB_SNE)

##### GET SSB ESTIMATES FOR ALL STOCKS COMBINED #####

SSB_allstocks<-list(SSB_EGOM,SSB_WGOM,SSB_GBK,SSB_SNE)
SSB_allstocks<-Reduce(function(x, y) merge(x, y, all=TRUE), SSB_allstocks)
SSB_allstocks<-aggregate(SSB~YEAR+SEASON,SSB_allstocks,FUN=sum)

##### combine seasons and EGOM+WGOM to compare to most recent stock assesment ssb plots###
SSB_GOM<-list(SSB_EGOM,SSB_WGOM)
SSB_GOM<-Reduce(function(x, y) merge(x, y, all=TRUE), SSB_GOM)
SSB_GOM<-aggregate(SSB~YEAR,SSB_GOM,FUN=sum)
#### plot estimates ####

df<-data.frame((SSB_allstocks[SSB_allstocks$SEASON =="FALL",]["YEAR"]),(SSB_allstocks[SSB_allstocks$SEASON =="FALL",]["SSB"]))
plot(df, type="l")

df<-data.frame((SSB_allstocks[SSB_allstocks$SEASON =="SPRING",]["YEAR"]),(SSB_allstocks[SSB_allstocks$SEASON =="SPRING",]["SSB"]))
plot(df, type="l")


#plot GOM estimates

plot(SSB_GOM$YEAR,SSB_GOM$SSB, type="l",main="NEFSC Trawl Data SSB Estimates ages 4+")

#### Make Final dfs and save #####
SSB_Fall_all<-SSB_allstocks[SSB_allstocks$SEASON =="FALL",]
SSB_Spring_all<-SSB_allstocks[SSB_allstocks$SEASON =="SPRING",]

SSB_Fall_EGOM<-SSB_EGOM[SSB_EGOM$SEASON =="FALL",]
SSB_Spring_EGOM<-SSB_EGOM[SSB_EGOM$SEASON =="SPRING",]

SSB_Fall_WGOM<-SSB_WGOM[SSB_WGOM$SEASON =="FALL",]
SSB_Spring_WGOM<-SSB_WGOM[SSB_WGOM$SEASON =="SPRING",]

SSB_Fall_GBK<-SSB_GBK[SSB_GBK$SEASON =="FALL",]
SSB_Spring_GBK<-SSB_GBK[SSB_GBK$SEASON =="SPRING",]

SSB_Fall_SNE<-SSB_SNE[SSB_SNE$SEASON =="FALL",]
SSB_Spring_SNE<-SSB_SNE[SSB_SNE$SEASON =="SPRING",]

#write.csv(SSB_Fall_all,here("data/SSB_estimates/SSB_Fall_all.csv"), row.names = FALSE)
#write.csv(SSB_Spring_all,here("data/SSB_estimates/SSB_Spring_all.csv"), row.names = FALSE)

#write.csv(SSB_Fall_EGOM,here("data/SSB_estimates/SSB_Fall_EGOM.csv"), row.names = FALSE)
#write.csv(SSB_Spring_EGOM,here("data/SSB_estimates/SSB_Spring_EGOM.csv"), row.names = FALSE)

#write.csv(SSB_Fall_WGOM,here("data/SSB_estimates/SSB_Fall_WGOM.csv"), row.names = FALSE)
#write.csv(SSB_Spring_WGOM,here("data/SSB_estimates/SSB_Spring_WGOM.csv"), row.names = FALSE)

#write.csv(SSB_Fall_GBK,here("data/SSB_estimates/SSB_Fall_GBK.csv"), row.names = FALSE)
#write.csv(SSB_Spring_GBK,here("data/SSB_estimates/SSB_Spring_GBK.csv"), row.names = FALSE)

#write.csv(SSB_Fall_SNE,here("data/SSB_estimates/SSB_Fall_SNE.csv"), row.names = FALSE)
#write.csv(SSB_Spring_SNE,here("data/SSB_estimates/SSB_Spring_SNE.csv"), row.names = FALSE)

#plot ssb estimates by stock area and region
par(mar=c(2,4,2,1), mfrow=c(2,4))
plot(SSB_Fall_EGOM$YEAR,SSB_Fall_EGOM$SSB, type="l",ylim=c(0,20),xlim=c(1982,2019),ylab="SSB",lwd=3,col="#00608A",main="Fall EGOM")
plot(SSB_Fall_WGOM$YEAR,SSB_Fall_WGOM$SSB, type="l",ylim=c(0,20),xlim=c(1982,2019),ylab="SSB",lwd=3,col="#00608A",main="Fall WGOM")
plot(SSB_Fall_GBK$YEAR,SSB_Fall_GBK$SSB, type="l",ylim=c(0,20),xlim=c(1982,2019),ylab="SSB",lwd=3,col="#00608A",main="Fall GBK")
plot(SSB_Fall_SNE$YEAR,SSB_Fall_SNE$SSB, type="l",ylim=c(0,20),xlim=c(1982,2019),ylab="SSB",lwd=3,col="#00608A",main="Fall SNE")
plot(SSB_Spring_EGOM$YEAR,SSB_Spring_EGOM$SSB, type="l",ylim=c(0,20),xlim=c(1982,2019),ylab="SSB",lwd=3,col="#00608A",main="Spring EGOM")
plot(SSB_Spring_WGOM$YEAR,SSB_Spring_WGOM$SSB, type="l",ylim=c(0,20),xlim=c(1982,2019),ylab="SSB",lwd=3,col="#00608A",main="Spring WGOM")
plot(SSB_Spring_GBK$YEAR,SSB_Spring_GBK$SSB, type="l",ylim=c(0,85),xlim=c(1982,2019),ylab="SSB",lwd=3,col="#00608A",main="Spring GBK")
plot(SSB_Spring_SNE$YEAR,SSB_Spring_SNE$SSB, type="l",ylim=c(0,20),xlim=c(1982,2019),ylab="SSB",lwd=3,col="#00608A",main="Spring SNE")
