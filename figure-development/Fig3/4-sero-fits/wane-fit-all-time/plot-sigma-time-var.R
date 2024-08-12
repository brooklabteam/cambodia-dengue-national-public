

rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(lubridate)
library(matrixcalc)


# at the province level, compare each run with and without age data
# and with and without waning immunity

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig3/4-sero-fits/wane-fit-all-time/"))
load("sigma-fit-all-time-var.Rdata")
#load("sigma-profile.Rdata")

#sigma.fit <- dplyr::select(sigma.fit, year, sigma, neg_llik, convergence)
#sigma.fit <- merge(sigma.fit, out.sigma.CI, by="year", all.x = T)
#head(sigma.fit)

p1 <- ggplot(sigma.fit) + geom_line(aes(x=year, y=sigma)) + 
                          geom_ribbon(aes(x=year, ymin=lci_sigma, ymax=uci_sigma), alpha=.3) + 
                          geom_point(aes(x=year, y=sigma)) #+ scale_y_log10()
p1

#and save 
write.csv(sigma.fit, file = paste0(homewd, "/data/sigma-fit.csv"), row.names = F)
