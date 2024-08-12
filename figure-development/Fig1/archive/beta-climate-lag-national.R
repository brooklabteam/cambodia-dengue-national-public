rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)

#set homewd
homewd = "/Users/carabrook/Developer/cambodia-dengue-national"

# load the national beta estimates from tsir
beta.df <- read.csv(file = paste0(homewd, "/data/beta_TSIR_fit_national.csv"), header = T, stringsAsFactors = F)

head(beta.df)
# pair with the case data

tsir.dat <- read.csv(paste0(homewd, "/data/tsir_dat_national.csv"), header = T, stringsAsFactors = F)

head(tsir.dat)
unique(tsir.dat$time- trunc(tsir.dat$time))
tsir.dat$biwk


#load the national-level climate data
