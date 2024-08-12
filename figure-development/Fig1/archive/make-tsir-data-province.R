


rm(list=ls())


library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(epitools)
library(reshape2)

#make national population and birth data

homewd= "/Users/carabrook/Developer/cambodia-dengue-national/"
setwd(homewd)

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
time = c(2002:2021)

#now expand into biweeks

make.biweekly.dat <- function(pop, births){
  pop.add <- pop[2:length(pop)] - pop[1:length(pop)-1]
  pop <- pop[2:length(pop)]
  
  births <- rep(birth.vec, each=26)/26
  
  pop.list <- as.list(pop)
  pop.add.list <- as.list(pop.add)
  out <- mapply(expand.pop, N=pop.list, N.add=pop.add.list, SIMPLIFY = F)
  pop.out <-c(unlist(out))
  
  #and for biweeks too
  time = rep(as.numeric(names(pop)[1]):as.numeric(names(pop)[length(names(pop))]), each=26)
  lts = as.numeric(names(pop)[length(names(pop))]) - as.numeric(names(pop)[1])
  biwk.timestep = cumsum(rep(1/26, 26))
  biwk.timestep <- biwk.timestep[1:length(biwk.timestep)-1]
  biwk.timestep <- round(c(0, biwk.timestep),2)
  biwk.timestep <- rep_len(biwk.timestep, length.out = length(time))
  time <- time + biwk.timestep
  
  #and join
  dat.out <- cbind.data.frame(time=time, births=births, pop=pop.out)
  #
  return(dat.out)
}
expand.pop <- function(N, N.add){
  
  N.new <- rep(N, each=26)#this is pop by biweek with each biweek getting the total at the end of each year
  
  start.N = N-N.add
  N.add.biweek <- cumsum(rep(N.add/26, 26))
  
  N.new <- rep(start.N, 26)
  N.new <- N.new+ N.add.biweek
  
  return(N.new)
  
}

dat.population <- make.biweekly.dat(pop=pop.vec, births = birth.vec)

head(dat.population)
rownames(dat.population) <- c()

# now scale by province based on population size
dat.population$provname <- "national"

load(paste0(homewd, "/data/cambodia_province_proportions.Rdata"))
head(prop.prov)


prop.list=dlply(prop.prov,.(provname))
make.prov.data <- function(prop.df, nat.dat){
  new.df = nat.dat
  new.df$provname <- unique(prop.df$provname)
  new.df$pop <- new.df$pop*prop.df$pop_prop
  new.df$births <- new.df$births*prop.df$pop_prop
  return(new.df)
  
}

prop.dat <- data.table::rbindlist(lapply(prop.list, make.prov.data, nat.dat=dat.population))
head(prop.dat)

prop.dat$year <- trunc(prop.dat$time)
prop.dat <- arrange(prop.dat, provname, time)
prop.dat = subset(prop.dat, provname!="Tboung Khmum")
unique(prop.dat$provname)
length(unique(prop.dat$year))
nrow(prop.dat)
length(rep(rep(1:26, 19), 24))
prop.dat$biwk <- rep(rep(1:26, 19), 24)

# now add cases and climate
case.dat <-read.csv(file = paste0(homewd, "/data/climate_cases_prov.csv"), header = T, stringsAsFactors = F)
head(case.dat)
case.dat$dates <- as.Date(case.dat$dates, format = "%m/%d/%y")

#and merge

merge.dat <- merge(prop.dat, case.dat, by=c("provname", "year", "biwk"), all.x = T)

head(merge.dat) # this dataset has everything...

write.csv(merge.dat, file = paste0(homewd, "/data/tsir_dat_province.csv"), row.names = F)



# case.dat <- read.csv(paste0(homewd, "data/case_data.csv"), header = T, stringsAsFactors = F)
# dat = case.dat
# 
# get.data.tsir <- function(dat, pop.dat){
#   
#   #make a vector of weeks, years, months
#   week <- week(dat$date)
#   week[week==53] <- 52 #now have a vector of 52 weeks in a year
#   year <- year(dat$date) #here we do the same for years
#   month <- month(dat$date) #and here the same for months
#   month.index <- rep(NA,52) #here make a vector of NA 52 times
#   for (j in 1:52) {month.index[j] <- mean(month[week==j]) }#here fill in that vector to assign the 52 weeks in the year to a given month (1-12)
#   month.index[1] <- 1; month.index[52] <- 12
#   
#   
#   u.yr <- unique(year) ; u.yr <- u.yr[order(u.yr)]  #here make a vector of the unique years in increasing order over which the time series spans
#   yr.index <- rep(u.yr,each=52) #here make the same vector that has a year assigned to each of those 52 entries. you now have 2 vectors, for the weeks and years the dataset spans
#   #make vector of number of positive cases by biweek 
#   #across the entire  time series for dengue (bwn)
#   
#   wks <- 1:52
#   biwks <- rep(1:26, each=2)
#   
#   dat$biwk <- dat$week
#   for (i in 1:52){
#     dat$biwk[dat$biwk==wks[i]] <- biwks[i]
#   }
#   
#   #summarize case data by biweek
#   date.sum <- ddply(dat, .(year, biwk), summarise, cases = length(year))
#   #head(date.sum)
#   index.dat <- cbind.data.frame(biwk=rep(1:26, length(u.yr)),year=rep(u.yr, each=26))
#   index.dat <- merge(index.dat, date.sum, by=c("year", "biwk"), all.x = T, sort = F)
#   index.dat$cases[is.na(index.dat$cases)] <- 0
#   
#   #vector of cases by biweek across the timeseries
#   bwn <- index.dat$cases
#   yrs <- rep(u.yr,each=26) #make a vector that is the same length and assigns a year to each of these 2 week periods
#   weeks <- rep(1:26,length(u.yr)) #do the same for weeks
#   months <- round(rep(0.5*(month.index[seq(1,52,by=2)]+month.index[seq(1,52,by=2)+1]), length(u.yr))) #do the same for months
#   
#   
#   #now join with the pop.dat
#   
#   out.dat <- cbind.data.frame(yr=yrs, month=months, biweek=weeks, cases=bwn)
#   
#   
#   #and get a time to correspond to each 
#   time.vec = round(cumsum(rep(1/26, 26)),2)
#   time.vec <- time.vec[1:length(time.vec)-1]
#   time.vec <- c(0, time.vec)
#   time.vec <- rep_len(time.vec, length.out = length(out.dat$yr))
#   out.dat$time <- out.dat$yr + time.vec
#   
#   #and merge
#   new.dat <- merge(out.dat, pop.dat, by="time", all.x = T, sort = F)
#   #and slim to what is needed
#   new.dat <- dplyr::select(new.dat, time, cases, births, pop, yr, biweek)
#   
#   
#   ## return things needed for TSIR
#   return(new.dat)
#   
#   
# }
# tsir.dat <- get.data.tsir(dat=case.dat, pop.dat=dat.population)

#and save this 
write.csv(tsir.dat, file = paste0(homewd, "/data/tsir_data.csv"), row.names = F)


