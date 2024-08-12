

rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(lubridate)
library(matrixcalc)


# at the province level, compare each run with and without age data
# and with and without waning immunity

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig3/4-sero-fits/ferg-fit-wane/"))
load("sigma-fit-all.Rdata")
sigma.fit # there is a value for sigma!
1/sigma.fit$sigma #927 years duration immunity (basically 0 immune waning!)

