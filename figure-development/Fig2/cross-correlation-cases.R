rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(geosphere)

#calculate pearson's cross correlation coefficient 
#between time series of two provinces. 
#First, do it annually with the annual case data
# And, annually, with the reconstructed annual cycle data
# Then do it in overlapping 5-year windows for the multi-annual cycle data

#first, load the time series data
homewd= "/Users/carabrook/Developer/cambodia-dengue-national-public"

#load tsir data
corr.dat <- read.csv(paste0(homewd, "/data/synchrony_case_data_oct15.csv"), header = T, stringsAsFactors = F)
head(corr.dat)
tail(corr.dat) #goes from beginning of 2002 to end of 2020

#link to centorid data so we can also calculate geographic
#distance between provinces


#load and attach centroid of each province
centroid.prov <- read.csv(file = paste0(homewd, "/data/centroid_provinces.csv"), header = T, stringsAsFactors = F)
head(centroid.prov)
centroid.merge <- dplyr::select(centroid.prov, -(mean_elevation_m))

setdiff(unique(corr.dat$provname), unique(centroid.merge$provname))
setdiff(unique(centroid.merge$provname), unique(corr.dat$provname))
corr.dat$provname[corr.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"

#merge
corr.dat <- merge(corr.dat, centroid.merge, by = "provname")




cross.corr <- function(df){
  out.df <- cbind.data.frame(year=unique(df$year), corr=cor(df$cases_this_prov, df$cases_other_prov))
  return(out.df)
  
}
cross.corr.cycle <- function(df){
  out.df <- cbind.data.frame(year=unique(df$year), corr=cor(df$annual_cycle_this_prov, df$annual_cycle_other_prov))
  return(out.df)
  
}
cross.corr.multi <- function(df){
  out.df <- cbind.data.frame(year_range=paste0(min(df$year),"-" , max(df$year)), mid_year=(min(df$year) +((max(df$year)-min(df$year))/2)), corr=cor(df$multi_this_prov, df$multi_other_prov))
  return(out.df)
  
}
assess.corr.annual <- function(df2, df1){
  #set aside prior data
  df1.hold = df1
  df2.hold = df2
  
  #first, join the two series for the full length that they match up
  df1 <- dplyr::select(df1, time, year, cases_per_100000)
  df2 <- dplyr::select(df2, time, year, cases_per_100000)
  names(df1) <- c("time", "year", "cases_this_prov")
  names(df2) <- c("time", "year", "cases_other_prov")
  df.join <- merge(df1,df2, by = c("time", "year"))
  
  #now get the full time series correlation
  full.corr <- cor(df.join$cases_this_prov, df.join$cases_other_prov)
  
  #then, split both by year and look within
  df.year.split <- dlply(df.join, .(year))
  
  df.year.out <- lapply(df.year.split, cross.corr)
  df.year.out <- data.table::rbindlist(df.year.out)
  #and assess cross corr within each year
  
  # now, add in the identifying details
  df.year.out$provname = unique(df1.hold$provname)
  df.year.out$comp_prov <- unique(df2.hold$provname)
  df.year.out$full_ts_corr <- full.corr
  #and add in the phylogenetic distance between these two
  df.year.out$dist_m <- distm(c(unique(df1.hold$longitude), unique(df1.hold$latitude)), c(unique(df2.hold$longitude), unique(df2.hold$latitude)), fun = distHaversine)
  df.year.out$dist_km <- df.year.out$dist_m/1000
  df.year.out$dist_from_PP_m <- distm(c(unique(df1.hold$longitude), unique(df1.hold$latitude)), c(104.8397, 11.58035), fun = distHaversine)
  df.year.out$dist_from_PP_km <- df.year.out$dist_from_PP_m/1000
  
  return(df.year.out)
}
assess.corr.annual.cycle <- function(df2, df1){
  #set aside prior data
  df1.hold = df1
  df2.hold = df2
  
  #first, join the two series for the full length that they match up
  df1 <- dplyr::select(df1, time, year, reconstructed_annual_period)
  df2 <- dplyr::select(df2, time, year, reconstructed_annual_period)
  names(df1) <- c("time", "year", "annual_cycle_this_prov")
  names(df2) <- c("time", "year", "annual_cycle_other_prov")
  df.join <- merge(df1,df2, by = c("time", "year"))
  
  #now get the full time series correlation
  full.corr <- cor(df.join$annual_cycle_this_prov, df.join$annual_cycle_other_prov)
  
  #then, split both by year and look within
  df.year.split <- dlply(df.join, .(year))
  
  df.year.out <- lapply(df.year.split, cross.corr.cycle)
  df.year.out <- data.table::rbindlist(df.year.out)
  #and assess cross corr within each year
  
  # now, add in the identifying details
  df.year.out$provname = unique(df1.hold$provname)
  df.year.out$comp_prov <- unique(df2.hold$provname)
  df.year.out$full_ts_corr <- full.corr
  #and add in the phylogenetic distance between these two
  df.year.out$dist_m <- distm(c(unique(df1.hold$longitude), unique(df1.hold$latitude)), c(unique(df2.hold$longitude), unique(df2.hold$latitude)), fun = distHaversine)
  df.year.out$dist_km <- df.year.out$dist_m/1000
  df.year.out$dist_from_PP_m <- distm(c(unique(df1.hold$longitude), unique(df1.hold$latitude)), c(104.8397, 11.58035), fun = distHaversine)
  df.year.out$dist_from_PP_km <- df.year.out$dist_from_PP_m/1000
  
  return(df.year.out)
}
assess.corr.5year <- function(df2, df1){
  #set aside prior data
  df1.hold = df1
  df2.hold = df2
  
  #first, join the two series for the full length that they match up
  df1 <- dplyr::select(df1, time, year, reconstructed_multi_period)
  df2 <- dplyr::select(df2, time, year, reconstructed_multi_period)
  names(df1) <- c("time", "year", "multi_this_prov")
  names(df2) <- c("time", "year", "multi_other_prov")
  df.join <- merge(df1,df2, by = c("time", "year"))
  
  #now get the full time series correlation
  full.corr <- cor(df.join$multi_this_prov, df.join$multi_other_prov)
  
  #then, split both by 5 year moving time windows and look within each one
  df.year.split <- list()
  for(i in 1:(length(unique(df.join$year))-4)){
    index.end = i + 4
    year.vect = unique(df.join$year)
    sub.dat = subset(df.join, year>=year.vect[i] & year<= year.vect[index.end])
    df.year.split[[i]] <- sub.dat
  }
  
  
  df.year.out <- lapply(df.year.split, cross.corr.multi)
  df.year.out <- data.table::rbindlist(df.year.out)
  #and assess cross corr within each year
  
  # now, add in the identifying details
  df.year.out$provname = unique(df1.hold$provname)
  df.year.out$comp_prov <- unique(df2.hold$provname)
  df.year.out$full_ts_corr <- full.corr
  #and add in the phylogenetic distance between these two
  df.year.out$dist_m <- distm(c(unique(df1.hold$longitude), unique(df1.hold$latitude)), c(unique(df2.hold$longitude), unique(df2.hold$latitude)), fun = distHaversine)
  df.year.out$dist_km <- df.year.out$dist_m/1000
  df.year.out$dist_from_PP_m <- distm(c(unique(df1.hold$longitude), unique(df1.hold$latitude)), c(104.8397, 11.58035), fun = distHaversine)
  df.year.out$dist_from_PP_km <- df.year.out$dist_from_PP_m/1000
  
  return(df.year.out)
}
assess.cross.corr <- function(provname1, df.all, cycle_type){
  df.other <- subset(df.all, provname!=provname1)
  
  #now split other by province 
  df.other <- arrange(df.other, provname, time)
  df.split <- dlply(df.other, .(provname))
  
  df.now = subset(df.all, provname==provname1)
  df.now <- arrange(df.now, provname, time)
  
  #and apply other function across this list of provinces to compare
  if(cycle_type=="annual_case"){
    out.provs.list <- lapply(df.split, assess.corr.annual, df1 = df.now)  
  }else if(cycle_type=="annual_cycle"){
    out.provs.list <- lapply(df.split, assess.corr.annual.cycle, df1 = df.now)  
  }else if(cycle_type=="multi"){
    out.provs.list <- lapply(df.split, assess.corr.5year, df1 = df.now)  
  }
  
  out.provs.df <- data.table::rbindlist(out.provs.list)
  
  #ggplot(out.provs.df) + geom_line(aes(x=year, y = corr, color = comp_prov)) + geom_line(aes(x=year, y=full_ts_corr, color=comp_prov),size=1) + facet_wrap(~comp_prov)
   return(out.provs.df) 
}

#now, split the provinces and run for all
provname.list <- as.list(unique(corr.dat$provname))

out.pearsons <- lapply(provname.list, assess.cross.corr, df.all=corr.dat, cycle_type = "annual_case")
pearsons.df <- data.table::rbindlist(out.pearsons)
head(pearsons.df)

#save these data
write.csv(pearsons.df, file = paste0(homewd, "/data/pearsons_correlations_provinces.csv"),row.names = F)

#and do the multi version
out.pearsons.multi <- lapply(provname.list, assess.cross.corr, df.all=corr.dat, cycle_type = "multi")
pearsons.df.multi <- data.table::rbindlist(out.pearsons.multi)
head(pearsons.df.multi)

#save these data
write.csv(pearsons.df.multi, file = paste0(homewd, "/data/pearsons_correlations_provinces_multi.csv"),row.names = F)


#and annual cycle data

out.pearsons.annual.cycle <- lapply(provname.list, assess.cross.corr, df.all=corr.dat, cycle_type = "annual_cycle")
pearsons.df.annual.cycle <- data.table::rbindlist(out.pearsons.annual.cycle)
head(pearsons.df.annual.cycle)

#save these data
write.csv(pearsons.df.annual.cycle, file = paste0(homewd, "/data/pearsons_correlations_provinces_annual_cycle.csv"),row.names = F)


