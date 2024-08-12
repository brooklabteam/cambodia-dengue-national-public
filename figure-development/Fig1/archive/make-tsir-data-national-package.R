rm(list=ls())


library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(epitools)
library(reshape2)
library(tsiR)

#make national population and birth data

homewd= "/Users/carabrook/Developer/cambodia-dengue-national/"
setwd(homewd)

#feed to tsir data
# time = time series of timesteps by 1 week
# cases = series of weekly cases (same length as time)
# births = annual time series of births
# population = annual time series of population

#load data 
#and get national pop data
pop.dat <- read.csv(file=paste0(homewd, "data/world_bank_cambodia.csv"), header = T, stringsAsFactors = F)
head(pop.dat)

#get population vector
pop.vec <- pop.dat[2,5:ncol(pop.dat)]
names(pop.vec) <- seq(1960,2020,1) #assume this is census at the end of each year,
#so one timestep before gives you the population at the beginning of the year
pop.vec <- c(unlist(pop.vec[which(names(pop.vec)=="2001"):which(names(pop.vec)=="2020")]))
#pop.vec[length(pop.vec)] <- pop.vec[length(pop.vec)-1]

#do the same for births - these are births per 1000 people
#get total births
birth.vec <- pop.dat[1,5:ncol(pop.dat)]
names(birth.vec) <- seq(1960,2020,1) #assume this is census at the end of each year,
birth.vec <- birth.vec[which(names(birth.vec)=="2002"):which(names(birth.vec)=="2020")]
birth.vec['2020'] <- birth.vec['2019'] #assume this is the same as prior year

#now scale up by population size to get total births per year
birth.vec = c(unlist(birth.vec*(pop.vec/1000)))
#time = c(2002:2021)

#and cases

dat <- read.csv(file = paste0(homewd, "/data/DENV-Dist-Diag.csv"), header=T, stringsAsFactors = F)
head(dat) 

#plot time series of each type by province by year
unique(dat$diagnostic) #df, dhf, dss
dat$date <- as.Date(dat$date, format = "%m/%d/%y")

sort(unique(dat$week_report)) # 1:52
sort(unique(dat$year)) # 2002: 2020

dat$provname <- "national"


#and sum by week of year and year and province
dat.sum <- ddply(dat,.(year, week_report), summarise, cases=sum(case))
head(dat.sum)


dat.sum$year <- as.factor(dat.sum$year)

dat.sum$provname = "national"
#and plot 
ggplot(dat.sum) + geom_line(aes(x=week_report, y=cases, color=year)) + 
  facet_wrap(~provname, scales = "free_y")

#now assign a "time" for each of the biweek values
time.df <- cbind.data.frame(week_report = 1:52, biweek =rep(1:26,each =2))
time.df <- arrange(time.df, week_report, biweek)
head(time.df)

#and add this to each
dat.sum <- merge(dat.sum, time.df, by="week_report", all.x = T)
#dat.sum$time <- as.numeric(dat.sum$time + as.numeric(as.character(dat.sum$year)) )
head(dat.sum)
#unique(dat.sum$time - trunc(dat.sum$time)) #26 possible timesteps per year
dat.sum$week_report <- as.numeric(as.character(dat.sum$week_report))

head(dat.sum)

dat.sum <- arrange(dat.sum, provname, year, biweek, week_report)

head(dat.sum)

ggplot(dat.sum) + geom_line(aes(x=week_report, y=cases, color=year)) + 
  facet_wrap(~provname, scales = "free_y")


# now write a function that multiplies pop and births by province proportion,
# creates a time column, and generates tsir data by province

load(paste0(homewd, "/data/cambodia_province_proportions.Rdata"))
head(prop.prov)
prop.prov <- rbind(prop.prov, c("national", 1))
prop.prov$pop_prop <- as.numeric(prop.prov$pop_prop)


build.prov.tsir <- function(prov.case.dat, nat.births, nat.pop, prop.prov){
  #first, scale the pop and birth vectors by prov
  prov.proportion = prop.prov$pop_prop[prop.prov==unique(prov.case.dat$provname)]
  prop.pop <- nat.pop*prov.proportion
  prop.births <- nat.births*prov.proportion
  
  #add time column based on week and biweek
  #time.df<- cbind.data.frame(week_report=1:52, time = rep((0:25)/26, each=2), biweek = rep(1:26, each=2))
  
  time.df <- cbind.data.frame(year = rep(2002:2020, each=52), week_report=rep(1:52, 19), time = rep(rep((0:25)/26, each=2), 19) , biweek = rep(rep(1:26, each=2),19))
  
  
  prov.case.dat <- merge(time.df, prov.case.dat, by=c("week_report", "biweek", "year"), all.x = T)
  
  #prov.case.dat <- merge(prov.case.dat, time.df, by=c("week_report", "biweek"), all.x = T)
  prov.case.dat$time <- as.numeric(as.character(prov.case.dat$year)) + prov.case.dat$time
  prov.case.dat$time <- as.numeric(as.character(prov.case.dat$time))
  
  prov.case.dat <- arrange(prov.case.dat, time)
  
  #head(prov.case.dat)
  
  prov.case.dat$cases[is.na(prov.case.dat$cases)] <- 1
  prov.case.dat$provname[is.na(prov.case.dat$provname)] <- unique(prov.case.dat$provname[!is.na(prov.case.dat$provname)])
  
  prov.sum = ddply(prov.case.dat, .(year, biweek), summarise, cases=sum(cases), time=unique(time))
  #head(prov.sum)
  
  prov.sum <- arrange(prov.sum, time)
  #prov.case.dat$time <- year(prov.case.dat$epiwk) + yday(prov.case.dat$epiwk)/365
  
  intbirths <- approx(prop.births, n = length(prov.sum$time))$y/(26)
  intpop <- approx(prop.pop, n = length(prov.sum$time))$y
  
  
  out.prov <- cbind.data.frame(time=prov.sum$time, cases = prov.sum$cases, births = intbirths, pop=intpop, year = prov.sum$year, biweek = prov.sum$biweek)
  out.prov$provname <- unique(prov.case.dat$provname)
  
  return(out.prov)
  
  
}


# and build tSIR dataset by province
tsir.prov <- build.prov.tsir(prov.case.dat = dat.sum, nat.births=birth.vec, nat.pop=pop.vec, prop.prov=prop.prov)


#and save - no climate, but this is fine
head(tsir.prov)

sort(unique(tsir.prov$time-trunc(tsir.prov$time)))  #26 timesteps!, 

#within a year and a province only 26
sort(unique(tsir.prov$time[trunc(tsir.prov$time)==2014]-trunc(tsir.prov$time[trunc(tsir.prov$time)==2014]))) #26

ggplot(tsir.prov) + geom_line(aes(x=time, y=cases))
tsir.prov$year <- as.factor(tsir.prov$year)
ggplot(tsir.prov) + geom_line(aes(x=biweek, y=cases, color=year)) + facet_wrap(~provname)
tsir.prov$year <- as.numeric(as.character(tsir.prov$year))

write.csv(tsir.prov, file = paste0(homewd, "/data/tsir_dat_national.csv"), row.names = F)



