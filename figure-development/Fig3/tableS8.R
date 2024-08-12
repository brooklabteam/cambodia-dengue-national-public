rm(list=ls())

library(plyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(mgcv)
library(reshape2)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national-public"
setwd(homewd)

fit.dat <- read.csv(file = paste0(homewd,"/data/prov-fits-FOI.csv"), header = T,stringsAsFactors = F)
head(fit.dat)
fit.dat$lambda <- signif(fit.dat$lambda,4)
fit.dat$CI <- paste0("[", signif(fit.dat$lci,3), "-", signif(fit.dat$uci,3), "]")
#fit.dat$llik <- round(fit.dat$llik,2)
TableS8 <- dplyr::select(fit.dat, provname, year, lambda, CI, llik)
head(TableS8)
names(TableS8)[names(TableS8)=="llik"] <- "llik_FOI"
#and load the one from comparison
com.dat <- read.csv(file=paste0(homewd,"/data/TableS8_part2.csv"), header = T, stringsAsFactors = F)
head(com.dat)
com.dat <- dplyr::select(com.dat, -(llik_FOI))
TableS8 <- merge(TableS8, com.dat, by=c("provname"))
TableS8$llik_FOI <- round(TableS8$llik_FOI, 0)
TableS8$llik_age <- round(TableS8$llik_age, 0)
TableS8$llik_sigma <- round(TableS8$llik_sigma, 0)
TableS8$llik_sigma_age <- round(TableS8$llik_sigma_age, 0)


TableS8$AIC_FOI <- round(TableS8$AIC_FOI, 0)
TableS8$AIC_age <- round(TableS8$AIC_age, 0)
TableS8$AIC_sigma <- round(TableS8$AIC_sigma, 0)
TableS8$AIC_sigma_age <- round(TableS8$AIC_sigma_age, 0)


#save 
write.csv(TableS8, file = paste0(homewd,"/data/TableS8_final.csv"), row.names = F)
