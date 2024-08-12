rm(list=ls())


library(ggplot2)
library(plyr)
library(dplyr)


homewd = "/Users/carabrook/Developer/cambodia-dengue-national"

#load the province-level transmission
beta.df <- read.csv(file = paste0(homewd, "/data/beta_TSIR_fit_national.csv"), header = T, stringsAsFactors = F)
head(beta.df)


#also load and prep the climate data at the province level
katwd = "/Users/carabrook/Developer/cambodia-dengue-province"
temp.dat <- read.csv(file = paste0(katwd, "/Climate_data/Cambodia_province_climate/data_final/biweek_temp.csv"), header = T, stringsAsFactors = F )
precip.dat <- read.csv(file = paste0(katwd, "/Climate_data/Cambodia_province_climate/data_final/biweek_ppt.csv"), header = T, stringsAsFactors = F )

names(temp.dat)[names(temp.dat)=="adm1_name"] <- "provname"
names(precip.dat)[names(precip.dat)=="adm1_name"] <- "provname"

#temp.dat = subset(temp.dat, provname!="Administrative unit not available")
#precip.dat = subset(precip.dat, provname!="Administrative unit not available")
head(temp.dat) # year and biweek
head(precip.dat) # year and biweek

#need to correct precip.dat for 2020 and 2021 (even though will be dropped in merge)
precip.dat$year <- as.numeric(as.character(precip.dat$year))
precip.dat$precip_mm[precip.dat$year>=2020] <- precip.dat$precip_mm[precip.dat$year>=2020]/1000
precip.dat$precip_m[precip.dat$year>=2020] <- precip.dat$precip_m[precip.dat$year>=2020]/1000

precip.dat$year <-as.factor(precip.dat$year)
ggplot(precip.dat) + geom_line(aes(x=biwk, y=precip_mm, color=year)) + facet_wrap(provname~.)

#still some funky stuff in Kep, Pailin, Phnom Penh, Preah Sihanouk

#sum at national level
temp.sum <- ddply(temp.dat, .(year, biwk), summarise, temp_C = mean(temp_C))
precip.sum <- ddply(precip.dat, .(year, biwk), summarise, precip_mm = sum(precip_mm))

# and merge with case data
names(temp.sum)[names(temp.sum)=="biwk"] <- "biweek"
names(precip.sum)[names(precip.sum)=="biwk"] <- "biweek"

head(beta.df)
unique(beta.df$year)
head(temp.sum)
head(precip.sum)

climate.merge <- merge(beta.df, temp.sum, by = c("year", "biweek"), all.x = T)
head(climate.merge)
climate.merge <- merge(climate.merge, precip.sum, by = c("year", "biweek"), all.x = T)
head(climate.merge)

climate.merge$year <- as.factor(climate.merge$year)
ggplot(climate.merge) + geom_line(aes(x=biweek, y=precip_mm, color=year)) + facet_wrap(provname~.)
ggplot(climate.merge) + geom_line(aes(x=biweek, y=temp_C, color=year)) + facet_wrap(provname~.)

#save this

write.csv(climate.merge, file = paste0(homewd, "/data/climate_beta_nat.csv"), row.names = F)

# and test for lags now 

#split by epiyear
climate.split <- dlply(climate.merge, .(epiyr))

# Look for cross-correlations between 
find.lags <- function(dat){
  
  #first, take out the actual epiyear that was not involved in the caculation
  dat = subset(dat, year!=unique(dat$epiyr))
  
  #precip
  dat.lag <- cbind.data.frame(lag = print(ccf(dat$precip_mm, dat$beta))$lag, acf=print(ccf(dat$precip_mm, dat$beta))$acf)
  dat.lag$variable <- "precip_mm"
  dat.lag$lag[dat.lag$acf==max(dat.lag$acf[dat.lag$lag<0 & dat.lag$lag>-27])] 
  
  
  
  dat2 = cbind.data.frame(lag = print(ccf(dat$temp_C, dat$beta))$lag, acf=print(ccf(dat$temp_C, dat$beta))$acf)
  dat2$variable <- "temp_C"
  dat2$lag[dat2$acf==max(dat2$acf[dat2$lag<0 & dat2$lag>-27])]
  dat.lag <- rbind(dat.lag, dat2)
  
  max.lag <- dlply(dat.lag, .(variable))
  get.lag <- function(df){
    lag = df$lag[df$acf==max(df$acf[df$lag<0 &df$lag>-27])]
    if(length(lag)>1){
      lag = max(lag[lag<0])
    }
    df.out = cbind.data.frame(variable=unique(df$variable), lag=lag)
    return(df.out)
  }
  max.lag <- data.table::rbindlist(lapply(max.lag, get.lag))
  
  max.lag$epiyr <- unique(dat$epiyr)
  max.lag$provname <- unique(dat$provname)
  
  return(max.lag)
  
  
}

lag.climate.out <- lapply(climate.split, find.lags)

lag.climate.df <- data.table::rbindlist(lag.climate.out)

#and save lag data
write.csv(lag.climate.df, file=paste0(homewd,"/data/lags_climate_nat.csv"), row.names = F)

#now make shifted dataset based on lags and save

ddply(lag.climate.df, .(variable), summarise, mean_lag=mean(lag), median_lag = median(lag)) 
#median lag for precip at national level is 2 biweeks (4 weeks) and for temp is 5 biweeks (10 weeks)
#for province, this was 2 and 7

#this is consistent with wagner. use it!

#make shifted dataset and save
head(climate.merge)

climate.shift <- climate.merge[8:nrow(climate.merge),] #starts in year 8 after 7 timestep shift
climate.shift$temp_C_lag <- climate.merge$temp_C[1:(length(climate.merge$temp_C)-7)]
head(climate.shift)

#and the shifted precip.
climate.shift$precip_mm_lag <- climate.merge$precip_mm[6:(length(climate.merge$precip_mm)-2)]

#write data and do lagged regression in another script
write.csv(climate.shift, file = paste0(homewd, "/data/lagged-nat-clim-beta.csv"), row.names = F)



