rm(list=ls())


library(ggplot2)
library(plyr)
library(dplyr)
library(mgcv)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national-public"

#load and prep the climate data at the province level
temp.dat <- read.csv(file = paste0(homewd, "/data/biweek_temp.csv"), header = T, stringsAsFactors = F )
precip.dat <- read.csv(file = paste0(homewd, "/data/biweek_ppt.csv"), header = T, stringsAsFactors = F )
head(temp.dat) # year and biweek
head(precip.dat) # year and biweek


temp.dat = subset(temp.dat, provname!="Administrative unit not available")
precip.dat = subset(precip.dat, provname!="Administrative unit not available")


#setdiff(unique(temp.dat$provname), unique(beta.df$provname)) #"Mondul Kiri"    "Oddar Meanchey" "Ratanak Kiri"   "Siemreap"       "Tboung Khmum"  
temp.dat$provname[temp.dat$provname=="Siemreap"] <- "Siem Reap"
temp.dat$provname[temp.dat$provname=="Oddar Meanchey"] <- "Otdar Meanchey"
precip.dat$provname[precip.dat$provname=="Siemreap"] <- "Siem Reap"
precip.dat$provname[precip.dat$provname=="Oddar Meanchey"] <- "Otdar Meanchey"

#and plot to check
temp.dat$year <-as.factor(temp.dat$year)

ggplot(temp.dat) + geom_line(aes(x=biwk, y=temp_C, color=year)) + facet_wrap(provname~.)

precip.dat$year <-as.factor(precip.dat$year)

ggplot(precip.dat) + geom_line(aes(x=biwk, y=precip_mm, color=year)) + facet_wrap(provname~.)

head(precip.dat)
head(temp.dat)


# and test whether precip and temp are increasing through time
precip.dat$provname <- as.factor(precip.dat$provname)
temp.dat$provname <- as.factor(temp.dat$provname)

#remove Tboung Khmum becuase it was created too recently
# precip.dat = subset(precip.dat, provname!="Tboung Khmum")
# temp.dat = subset(temp.dat, provname!="Tboung Khmum")


# Here, we look at precipitation through time, allowing for 
# a different intercept for each province

temp.dat$year <- as.numeric(as.character(temp.dat$year))
precip.dat$year <- as.numeric(as.character(precip.dat$year))

gam1 <- gam(precip_mm~ year:provname + #slope specific by province
                       s(biwk, bs="cc", k=7) + #controls for internal annual cycles
                       s(provname, bs="re"), data=precip.dat) #y-intercept specific by province too
summary(gam1) 
#no change in precip over time 

#send to tables
tableS1b <- cbind.data.frame(GAM =rep("interannual_precip", length(summary(gam1)$p.coeff)), provname=sapply(strsplit(names(summary(gam1)$p.coeff), "provname"), "[", 2), slope = round(summary(gam1)$p.coeff, 3), CI= paste0("[",round(summary(gam1)$p.coeff-1.96*summary(gam1)$se[1:length(summary(gam1)$p.coeff)], 3), "-", round(summary(gam1)$p.coeff+1.96*summary(gam1)$se[1:length(summary(gam1)$p.coeff)], 3), "]"),  stat=round(summary(gam1)$p.t,2), p_val=round(summary(gam1)$p.pv,3))
tableS1b$dev_explained <- round(summary(gam1)$dev.expl,3)*100

#...are there any years that are significant deviations?
precip.dat$year <- as.factor(precip.dat$year)
gam1b <- gam(precip_mm~ s(year, bs="re") +
                        s(biwk, bs="cc", k=7) + #controls for internal annual cycles
                        s(provname, bs="re"), data=precip.dat) #y-intercept specific by province too
summary(gam1b)
AIC(gam1, gam1b)

tableS1D <- cbind.data.frame(GAM = rep("anomaly_precip",3), smooth = c("year", "biweek", "province"), deg_freedom = round(summary(gam1b)$s.table[,1], 2), stat = round(summary(gam1b)$s.table[,3], 1), p_val = round(summary(gam1b)$s.table[,4], 3))



source(paste0(homewd, "/figure-development/Fig1/mollentze-streicker-2020-functions.R"))

year.df <- get_partial_effects(gam1b, var="year")
plot.partial(df=year.df, var="year", response_var = "precip_mm")

precip.dat$year <- as.numeric(as.character(precip.dat$year))
precip.predict <- cbind.data.frame(provname = rep(unique(precip.dat$provname), each=18), year = rep(unique(precip.dat$year), 25))
precip.predict$biwk = 10
precip.predict$precip_mm <- predict.gam(gam1, newdata = precip.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$fit
precip.predict$precip_mm_lci <- predict.gam(gam1, newdata = precip.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$fit - 1.96*predict.gam(gam1, newdata = precip.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$se
precip.predict$precip_mm_uci <- predict.gam(gam1, newdata = precip.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$fit + 1.96*predict.gam(gam1, newdata = precip.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$se


# Here, we look at temperature through time, allowing for
# a different intercept for each province
gam2 <- gam(temp_C ~ year:provname + 
              s(biwk, bs="cc", k=7) + 
              s(provname, bs="re"),
              data=temp.dat)

summary(gam2) 
# temp is increasing slightly for all provinces
tableS1 <- cbind.data.frame(GAM =rep("interannual_temp", length(summary(gam2)$p.coeff)), provname=sapply(strsplit(names(summary(gam2)$p.coeff), "provname"), "[", 2), slope = round(summary(gam2)$p.coeff, 3), CI= paste0("[",round(summary(gam2)$p.coeff-1.96*summary(gam2)$se[1:length(summary(gam2)$p.coeff)], 3), "-", round(summary(gam2)$p.coeff+1.96*summary(gam2)$se[1:length(summary(gam2)$p.coeff)], 3), "]"),  stat=round(summary(gam2)$p.t,2), p_val=round(summary(gam2)$p.pv,3))
tableS1$dev_explained <- round(summary(gam2)$dev.expl,3)*100




#are there any temperature years that are temperature anomalies?
temp.dat$year <- as.factor(temp.dat$year)
gam2b <- gam(temp_C ~ s(year, bs="re") +
              s(biwk, bs="cc", k=7) + 
              s(provname, bs="re"),
            data=temp.dat)

summary(gam2b) 
year.df.temp <- get_partial_effects(gam2b, var="year")
plot.partial(df=year.df.temp, var="year", response_var = "temp_C") #2019 and 2012 were hot years but 2007 was not
temp.dat$year <- as.numeric(as.character(temp.dat$year))

tableS1C <- cbind.data.frame(GAM = rep("anomaly_temp",3), smooth = c("year", "biweek", "province"), deg_freedom = round(summary(gam2b)$s.table[,1], 2), stat = round(summary(gam2b)$s.table[,3], 1), p_val = round(summary(gam2b)$s.table[,4], 3))

tableS1A <- rbind(tableS1, tableS1b)
tableS1B <- rbind(tableS1C,tableS1D) 

write.csv(tableS1A, file = paste0(homewd, "/data/tableS1_part1.csv"), row.names = F)
write.csv(tableS1B, file = paste0(homewd, "/data/tableS1_part2.csv"), row.names = F)
# here are predictions for best fit model
temp.predict <- cbind.data.frame(provname = rep(unique(temp.dat$provname), each=18), year = rep(unique(temp.dat$year), 25))

temp.predict$biwk = 10
temp.predict$temp_C <- predict.gam(gam2, newdata = temp.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$fit
temp.predict$temp_C_lci <- predict.gam(gam2, newdata = temp.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$fit - 1.96*predict.gam(gam2, newdata = temp.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$se
temp.predict$temp_C_uci <- predict.gam(gam2, newdata = temp.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$fit + 1.96*predict.gam(gam2, newdata = temp.predict, exclude = c("s(biwk)"), type="response", se.fit = TRUE)$se


# and merge with case data
temp.dat <- dplyr::select(temp.dat, -(temp_K))
names(temp.dat)[names(temp.dat)=="biwk"] <- "biweek"
names(precip.dat)[names(precip.dat)=="biwk"] <- "biweek"


# now load the transmission data to calculate lags for tsir
beta.df <- read.csv(file = paste0(homewd, "/data/beta_TSIR_fit_province.csv"), header = T, stringsAsFactors = F)
head(beta.df)


unique(beta.df$year)

climate.merge <- merge(beta.df, temp.dat, by = c("provname", "year", "biweek"), all.x = T)
head(climate.merge)
climate.merge <- merge(climate.merge, precip.dat, by = c("provname", "year", "biweek"), all.x = T)
head(climate.merge)

climate.merge$year <- as.factor(climate.merge$year)

#change provname for plotting
climate.merge$provname[climate.merge$provname=="Otdar Meanchey"] <- "Oddar Meanchey"


#drop the sus reconstruction info
climate.merge <- dplyr::select(climate.merge, -(rsquared), -(sus_reconstruction))

#save this for calculating regression
write.csv(climate.merge, file = paste0(homewd, "/data/climate_beta_prov.csv"), row.names = F)


#now merge climate data with cases and population data for wavelet analysis
head(temp.dat)
head(precip.dat)

temp.dat$year <- as.numeric(as.character(temp.dat$year))
precip.dat$year <- as.numeric(as.character(precip.dat$year))

#load tsir data
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_dat_province.csv"), header = T, stringsAsFactors = F)
head(tsir.dat)
unique(tsir.dat$provname)
tsir.dat$provname[tsir.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"
precip.dat$provname <- as.character(precip.dat$provname)
temp.dat$provname <- as.character(temp.dat$provname)
precip.dat$provname[precip.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"
temp.dat$provname[temp.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"

precip.merge <- dplyr::select(precip.dat, -(precip_m))


#and merge on almost everything
wave.merge <- merge(tsir.dat, precip.merge, by=c("provname", "year", "biweek"))
head(wave.merge)

wave.merge <- merge(wave.merge, temp.dat, by=c("provname", "year", "biweek"))
head(wave.merge)

wave.merge <- arrange(wave.merge, provname, time)

#and save
write.csv(wave.merge, file = paste0(homewd, "/data/wavelet_input_dat_prov.csv"), row.names = F)


wave.merge$year <- as.factor(wave.merge$year)
pS2 <- ggplot(wave.merge) + geom_line (aes(x=biweek, y=precip_mm, color=year)) + facet_wrap(provname~.) + 
        theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
        ylab("biweekly sum precipitation (mm)") + xlab("biweek of year")

ggsave(file = paste0(homewd,"/final-figures/FigS2.png"),
       plot = pS2,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)


pS1 <- ggplot(wave.merge) + geom_line(aes(x=biweek, y=temp_C, color=year)) + facet_wrap(provname~.) + 
        theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
        ylab(bquote('biweekly mean temperature ('^0~'C)')) + xlab("biweek of year")

ggsave(file = paste0(homewd,"/final-figures/FigS1.png"),
       plot = pS1,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)


# and plot the time series foward
head(climate.merge)
head(precip.predict)
precip.predict$time <- precip.predict$year +.5
temp.predict$time <- temp.predict$year +.5


precip.predict$provname <- as.character(precip.predict$provname)
temp.predict$provname <- as.character(temp.predict$provname)
precip.predict$provname[precip.predict$provname=="Otdar Meanchey"] <- "Oddar Meanchey"
temp.predict$provname[temp.predict$provname=="Otdar Meanchey"] <- "Oddar Meanchey"

precip.predict$provname <- as.factor(precip.predict$provname)
temp.predict$provname <- as.factor(temp.predict$provname)


pS4 <- ggplot(wave.merge) +
       geom_line(aes(x=time, y=precip_mm, color=provname), show.legend = F) + facet_wrap(provname~.) + 
       theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"), axis.title.x = element_blank())+
       geom_ribbon(data=precip.predict, aes(time, ymin=precip_mm_lci, ymax=precip_mm_uci),alpha=.4) +
       geom_line(data=precip.predict, aes(time, precip_mm), size=.8) + ylab("biweekly sum precipitation (mm)")
  
ggsave(file = paste0(homewd,"/final-figures/FigS4.png"),
       plot = pS4,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)


pS3 <- ggplot(wave.merge) +
    geom_line(aes(x=time, y=temp_C, color=provname), show.legend = F) + facet_wrap(provname~.) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"), axis.title.x = element_blank())+
    geom_ribbon(data=temp.predict, aes(time, ymin=temp_C_lci, ymax=temp_C_uci),alpha=.4) +
    geom_line(data=temp.predict, aes(time, temp_C), size=.8) + ylab(bquote('biweekly mean temperature ('^0~'C)'))

ggsave(file = paste0(homewd,"/final-figures/FigS3.png"),
       plot = pS3,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)

#climate predictors by province - just to see it
p2 <- ggplot(temp.predict) +
  theme_bw() + theme(panel.grid = element_blank(), axis.title.x = element_blank(), legend.title = element_blank(),
                     axis.title.y = element_text(size=16), axis.text = element_text(size=13),
                     plot.margin = unit(c(.2,.3,.1,.3), "cm"), legend.position = "bottom")+
  geom_ribbon(data=temp.predict, aes(time, ymin=temp_C_lci, ymax=temp_C_uci, fill=provname),alpha=.4) +
  geom_line(data=temp.predict, aes(time, temp_C, color=provname), size=.8) + 
  ylab(bquote('annual mean temperature prediction ('^0~'C)')) #+ 
  #guides(fill=guide_legend(ncol=1), color= guide_legend(ncol=1))

print(p2)


p3 <- ggplot(precip.predict) +
  theme_bw() + theme(panel.grid = element_blank(), axis.title.x = element_blank(), legend.title = element_blank(),
                     axis.title.y = element_text(size=16), axis.text = element_text(size=13),
                     plot.margin = unit(c(.2,.3,.1,.3), "cm"), legend.position = "bottom")+
  geom_ribbon(data=precip.predict, aes(time, ymin=precip_mm_lci, ymax=precip_mm_uci, fill=provname),alpha=.4) +
  geom_line(data=precip.predict, aes(time, precip_mm, color=provname), size=.8) + 
  ylab("sum biweekly precipitation (mm) prediction by year") #+ 
#guides(fill=guide_legend(ncol=1), color= guide_legend(ncol=1))

print(p3)

# now, pair with the population data and drop the beta info - send off to the wavelet analysis

#now plot temp and precip as Zscores in next script




# and test for lags now - split by province first

climate.split <- dlply(climate.merge, .(provname, epiyr))

# Look for cross-correlations between 
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

ddply(lag.climate.df, .(variable), summarise, mean_lag=mean(lag), median_lag = median(lag)) 
#mean precip = -5, temp = -6
#median precip = -2, temp = -7

#and save lag data
write.csv(lag.climate.df, file=paste0(homewd,"/data/lags_climate_prov.csv"), row.names = F)

#and manipulate to make table
library(reshape2)
melt(lag.climate.df)
df.list <- dlply(lag.climate.df, .(provname))
make.table <- function(df){
  df2 <- cbind.data.frame(provname=unique(df$provname), temp_2007 = df$lag[df$variable=="temp_C"  & df$epiyr==2007], precip_2007 = df$lag[df$variable=="precip_mm"  & df$epiyr==2007], temp_2012 = df$lag[df$variable=="temp_C"  & df$epiyr==2012], precip_2012 = df$lag[df$variable=="precip_mm"  & df$epiyr==2012], temp_2019 = df$lag[df$variable=="temp_C"  & df$epiyr==2019], precip_2019 = df$lag[df$variable=="precip_mm"  & df$epiyr==2019])
  return(df2)
}
tableS3 <- data.table::rbindlist(lapply(df.list, make.table))

write.csv(tableS3, file = paste0(homewd, "/data/tableS3.csv"), row.names = F)
#make lagged data based on optimal lags by province (here, different from Wagner et al. 2020)

#make shifted dataset and save
head(climate.merge)
head(lag.climate.df)
climate.merge <- arrange(climate.merge, provname, time)
lag.climate.df <- arrange(lag.climate.df, provname, epiyr)

#and plot
climate.merge.split <- dlply(climate.merge,.(provname, epiyr))
lag.climate.split <- dlply(lag.climate.df, .(provname, epiyr)) #same length

make.lag.dat <- function(df, lag.df){
  #first, convert all (+) values to 0
  lag.df$lag[lag.df$lag>0] <- 0
  lag.df$lag <- abs(lag.df$lag)
  max.lag <- max(lag.df$lag)
  
  #and shift
  shift.dat <- df[(max.lag+1):nrow(df),]
  if(max.lag==lag.df$lag[lag.df$variable=="temp_C"]){
    shift.dat$temp_C_lag <- df$temp_C[1:(length(df$temp_C)-max.lag)]
    precip.start = (((length(df$precip_mm)-lag.df$lag[lag.df$variable=="precip_mm"]))-length(shift.dat$temp_C_lag) +1)
    shift.dat$precip_mm_lag <- df$precip_mm[precip.start:((length(df$precip_mm)-lag.df$lag[lag.df$variable=="precip_mm"]))]
  }else if(max.lag==lag.df$lag[lag.df$variable=="precip_mm"]){
    shift.dat$precip_mm_lag <- df$precip_mm[1:(length(df$precip_mm)-max.lag)]
    temp.start = (((length(df$temp_C)-lag.df$lag[lag.df$variable=="temp_C"]))-length(shift.dat$precip_mm) +1)
    shift.dat$temp_C_lag <- df$temp_C[temp.start:((length(df$temp_C)-lag.df$lag[lag.df$variable=="temp_C"]))]
  }
  shift.dat <- dplyr::select(shift.dat, provname, year, biweek, epiyr, time, cases, beta, betalow, betahigh, temp_C, precip_m, precip_mm, temp_C_lag, precip_mm_lag)
  return(shift.dat)
}

#and apply
climate.lag.dat.split <- mapply(make.lag.dat, df = climate.merge.split, lag.df=lag.climate.split, SIMPLIFY = F)

climate.lag.df <- data.table::rbindlist(climate.lag.dat.split)
climate.lag.df <- arrange(climate.lag.df, provname, time)
# climate.shift <- climate.merge[8:nrow(climate.merge),] #starts in year 8 after 7 timestep shift
# climate.shift$temp_C_lag <- climate.merge$temp_C[1:(length(climate.merge$temp_C)-7)]
# head(climate.shift)
# 
# #and the shifted precip.
# climate.shift$precip_mm_lag <- climate.merge$precip_mm[6:(length(climate.merge$precip_mm)-2)]

#write data and do lagged regression in another script
write.csv(climate.lag.df, file = paste0(homewd, "/data/lagged-prov-clim-beta.csv"), row.names = F)

