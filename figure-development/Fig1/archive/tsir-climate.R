


rm(list=ls())
## tsiR examples

## load the package and dependencies
## kernlab is for the gaussian process 
## the rest are for plotting 
require(tsiR)
require(kernlab)
require(ggplot2)
require(reshape2)
require(grid)
require(plyr)
library(lubridate)

#epidemic years 2007, 2012, 2019

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
#homewd = "/home/rstudio/"
setwd(paste0(homewd))

#load lagged national climate data and try to fit TSIR... as a whole piece???
# 
# dat <- read.csv(file = paste0(homewd, "/data/lagged-nat-clim.csv"), header = T, stringsAsFactors = F)
# 
# head(dat)
# 
# dat <- arrange(dat, dates)


#load the functions
#here, fit the model in an "endemic" period to predict an epidemic year

#fits with an increased S
sim.com.epi.clim.nat.lag <- function(time.start1, time.start2, time.start3, dat, epiyr1,epiyr2,epiyr3, family){
  #
  
  if(family=="poisson"){
    
    
    fittedpars1 <- estpars(data=subset(dat, time >= time.start1 & time<epiyr1),
                           IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                           regtype='lm',family='poisson',link='log')
    
    fittedpars2 <- estpars(data=subset(dat, time >= time.start2 & time<epiyr2),
                           IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                           regtype='lm',family='poisson',link='log')
    
    fittedpars3 <- estpars(data=subset(dat, time >= time.start3 & time<epiyr3),
                           IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                           regtype='lm',family='poisson',link='log')
    
    ## The type of regression used in susceptible reconstruction. 
    # Options are 'gaussian', 'lm' (linear model), 'spline' (smooth.spline with 2.5 degrees freedom), 'lowess' (with f = 2/3, iter = 1), 'loess' (degree 1), and 'user' which is just a user inputed vector. Defaults to 'gaussian' and if that fails then defaults to loess.
    
    ## The family in the GLM regression. 
    # One can use any of the GLM ones, but the options are essentially 'poisson' (with link='log'), 'gaussian' (with link='log' or 'identity'), or 'quasipoisson' (with link='log'). Default is 'gaussian'.
    
    # simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
    #                           IP = 2,
    #                           parms=fittedpars,
    #                           #epidemics='break', threshold=3,
    #                           method='pois', nsim=100)
    
  }else if(family=="gaussian"){
    
    
    fittedpars1 <- estpars(data=subset(dat, time >= time.start1 & time<epiyr1),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='gaussian')
    
    fittedpars2 <- estpars(data=subset(dat, time >= time.start2 & time<epiyr2),
                           IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                           regtype='lm',family='gaussian')
    
    fittedpars3 <- estpars(data=subset(dat, time >= time.start3 & time<epiyr3),
                           IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                           regtype='lm',family='gaussian')
    
  }
    #now, get regression based on climate for the transmission rate
    #take time series of transmission rate (constant across the years but repeating each biweek)
    #and calculate cross-correlation with the climate variables
    
    #here, make time series of transmission rate
    beta.ts1 <- rep(fittedpars1$beta, (epiyr1-time.start1))
    beta.ts2 <- rep(fittedpars2$beta, (epiyr2-time.start2))
    beta.ts3 <- rep(fittedpars3$beta, (epiyr3-time.start3))
    
    dat.ts = subset(dat, year!=epiyr1 & year!=epiyr2 &year<epiyr3)
    
    
    #head(dat.ts)
    dat.ts$beta <- c(beta.ts1, beta.ts2, beta.ts3)
    
    
    #now get lags and shift dataset
    dat.lag <- cbind.data.frame(lag = print(ccf(dat.ts$precip_mm, dat.ts$beta, lag.max = 26))$lag, acf=print(ccf(dat.ts$precip_mm, dat.ts$beta, lag.max = 26))$acf)
    dat.lag$variable <- "precip_mm"
    
    #sub to only neg lags (climate precedes cases)
    dat.lag = subset(dat.lag, lag<0)
    
    dat.lag.precip <- min(dat.lag$lag[dat.lag$acf==max(dat.lag$acf)])#minimum lag which maximizes cross corr
    
    dat.lag.precip <- abs(dat.lag.precip)
    print(paste0("precip precedes cases by ", dat.lag.precip, " biweeks"))
    
    dat2 = cbind.data.frame(lag = print(ccf(dat.ts$temp_C, dat.ts$beta, lag.max = 26))$lag, acf=print(ccf(dat.ts$temp_C, dat.ts$beta, lag.max = 26))$acf)
    dat2$variable <- "temp_C"
    dat2 = subset(dat2, lag<0)
    dat.lag.temp <- min(dat2$lag[dat2$acf==max(dat2$acf)]) #minimum lag which maximizes cross corr
    dat.lag.temp <- abs(dat.lag.temp)
    
    print(paste0("temp precedes cases by ", dat.lag.temp, " biweeks"))
    
    lag.df <- cbind.data.frame(lag_var=c("precip_mm", "temp_C"), lag=c(dat.lag.precip, dat.lag.temp))
    
    return(lag.df)
  }
  

sim.com.epi.clim <- function(time.start, lag.dat, dat, epiyr, family){
  #
  
  if(family=="poisson"){
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='poisson',link='log')
    ## The type of regression used in susceptible reconstruction. 
    # Options are 'gaussian', 'lm' (linear model), 'spline' (smooth.spline with 2.5 degrees freedom), 'lowess' (with f = 2/3, iter = 1), 'loess' (degree 1), and 'user' which is just a user inputed vector. Defaults to 'gaussian' and if that fails then defaults to loess.
    
    ## The family in the GLM regression. 
    # One can use any of the GLM ones, but the options are essentially 'poisson' (with link='log'), 'gaussian' (with link='log' or 'identity'), or 'quasipoisson' (with link='log'). Default is 'gaussian'.
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  }else if(family=="gaussian"){
    
    
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='gaussian')
    
    
    #now, get regression based on climate for the transmission rate
    #take time series of transmission rate (constant across the years but repeating each biweek)
    #and calculate cross-correlation with the climate variables
    
    #here, make time series of transmission rate
    beta.ts <- rep(fittedpars$beta, (epiyr-time.start))
    dat.ts = subset(dat, time >= time.start & time<epiyr)
    
    #head(dat.ts)
    dat.ts$beta <- beta.ts
    
    
    dat.lag.precip = lag.dat$lag[lag.dat$lag_var=="precip_mm"]
    dat.lag.temp = lag.dat$lag[lag.dat$lag_var=="temp_C"]
    #dat.lag <- rbind(dat.lag, dat2)
    
    #now make the appropriately lagged dataset - 
    #make it for the entire time series to save yourself a headache later
    if(dat.lag.precip>dat.lag.temp){
      merge.shift <- dat[dat.lag.precip:length(dat$year),] 
      #head(merge.shift)
      merge.shift$precip_lag <- dat$precip_mm[1:(length(dat$precip_mm)-(dat.lag.precip-1))] 
      merge.shift$meantempLag <- dat$temp_C[(dat.lag.precip-dat.lag.temp+1):(length(dat$temp_C)-(dat.lag.temp-1))]
      
    }else if(dat.lag.precip<dat.lag.temp){
      merge.shift <- dat[dat.lag.temp:length(dat$year),] 
      #head(merge.shift)
      
      merge.shift$meantempLag <- dat$temp_C[1:(length(dat$temp_C)-(dat.lag.temp-1))]
      merge.shift$precip_lag <- dat$precip_mm[(dat.lag.temp-dat.lag.precip+1):(length(dat$precip_mm)-(dat.lag.precip-1))]
      
    }
    
    #head(merge.shift)  # here is your lagged dataset for regression / tsir
    
    
    
    #add "month" colum
    month.df <- cbind.data.frame(month=rep(1:12, each=2), biwk=1:24)
    month.df <- rbind(month.df, c(12,25))
    month.df$biwk <- month.df$biwk+1
    month.df <- rbind(c(1,1), month.df)
    
    beta.df <- dplyr::select(dat.ts, biwk, beta)
    beta.df <- beta.df[!duplicated(beta.df),]
    
    merge.shift <- merge(merge.shift, beta.df, by=c("biwk") , all.x = T)
    merge.shift <- merge(merge.shift,month.df, by=c("biwk") , all.x = T)
    
    merge.shift <- arrange(merge.shift, time)
    #head(merge.shift)  # here is your lagged dataset for regression / tsir
    
    merge.shift$beta_log = log(merge.shift$beta +1)
    merge.shift$meantempLag_log <- log(merge.shift$meantempLag)
    merge.shift$precip_lag_log<- log(merge.shift$precip_lag)
    
    #head(merge.shift)  # here is your lagged dataset for regression / tsir
    merge.fit <- subset(merge.shift, time >= time.start & time<epiyr)
    
    #this will be the panel regression once we feel good about the province-level data
    #will include interaction predictiors of province by year and provice by month
    
    m1 <- glm(beta_log~precip_lag+meantempLag + year + month, data=merge.fit)
    m2 <- glm(beta_log~precip_lag_log+meantempLag_log + year +month, data=merge.fit)
    
    #AIC(m1,m2)
    dat.pred= subset(merge.shift, time >= epiyr & time<=(epiyr+1))
    
    #overwrite beta
    dat.pred$beta_log <- predict(m2, newdata = dat.pred, type="response")
    
    dat.pred$beta <- exp(dat.pred$beta_log)-1
    
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              #method='pois',
                              nsim=100)
    
    
    
  } 
  #saveRDS(fittedpars, paste0(homewd,"/sim.com.epi_estpars_",epiyr,".RDS") )
  #saveRDS(simfitted, paste0(homewd,"/sim.com.epi_",epiyr,".RDS") )
  
  
  
  
  Its=simfitted$res$mean*fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time=simfitted$res$time, I=Its)
  Ifinal<-Its[length(Its)]
  
  
  
  #now try to predict the epidemic both with and without the increased S
  #compute the most likely model fit with the proper increase which gives the 
  #closest match of model to observed cases
  
  fracincvec<-seq(1,3,0.01)
  #fracincdist<-matrix(0,1,1)
  #ratioImax<-matrix(0,1,length(fracincvec))
  #ratioImaxdist<-matrix(0,1,1)
  #Send<-matrix(0,1,1)
  sm.sq <- rep(NA, length(fracincvec))
  sm.sq.clim <- rep(NA, length(fracincvec))
  for(i in 1:length(fracincvec)){
    
    Sbegin=simfitted$sbar+simfitted$Z
    Sepi<-Sbegin[length(Sbegin)]
    SepiInc<-fracincvec[i]*Sbegin[length(Sbegin)]
    
    #this is for no increased S
    dat.fit = subset(dat, time >= time.start & time<epiyr)
    
    predict_ts <- predicttsir(times=dat.fit$time,
                              births = dat.fit$births,
                              beta = fittedpars$beta,
                              alpha = fittedpars$alpha,
                              S0 = Sbegin[1],
                              I0 = dat.fit$cases[1],
                              nsim=100,
                              stochastic = T)
    
    
    
    #and the epi year-no increase with climate
    predict_epi_noInc_clim <- predicttsir(times=dat.pred$time,
                                     births = dat.pred$births,
                                     beta = dat.pred$beta, #this the climate-driven prediction for beta
                                     alpha = fittedpars$alpha,
                                     S0 = Sepi,
                                     I0 = Ifinal,
                                     nsim=100,
                                     stochastic = T)
    
    #saveRDS(predict_epi_noInc, paste0(homewd,"/sim.com.epi_predict_epi_noInc_",epiyr,".RDS") )
    
    #and including an increase with climare
    predict_epi_Inc_clim <- predicttsir(times=dat.pred$time,
                                   births = dat.pred$births,
                                   beta = dat.pred$beta, #this the climate-driven prediction for beta
                                   alpha = fittedpars$alpha,
                                   S0 = SepiInc,
                                   I0 = Ifinal,
                                   nsim=100,
                                   stochastic = T)
    
    
    #and no climate - no increase
    predict_epi_noInc <- predicttsir(times=dat.pred$time,
                                          births = dat.pred$births,
                                          beta = fittedpars$beta, #this without climate in beta
                                          alpha = fittedpars$alpha,
                                          S0 = Sepi,
                                          I0 = Ifinal,
                                          nsim=100,
                                          stochastic = T)
    
    #no climate - yes increase
    
    predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                        births = dat.pred$births,
                                        beta = fittedpars$beta, #this without climate in beta
                                        alpha = fittedpars$alpha,
                                        S0 = SepiInc,
                                        I0 = Ifinal,
                                        nsim=100,
                                        stochastic = T)

    

    #now look at data and model predictions for the epidemic year only
    IPredEpi = cbind.data.frame(time=predict_epi_Inc$I$time, mean_Inc=predict_epi_Inc$I$mean)
    names(IPredEpi) <- c("time", "model_predicted_absolute_cases")
    IPredEpi$model_predicted_reported_cases <- IPredEpi$model_predicted_absolute_cases/fittedpars$rho[length(fittedpars$rho)]
    
    
    IPredEpiClim = cbind.data.frame(time=predict_epi_Inc_clim$I$time, mean_Inc=predict_epi_Inc_clim$I$mean)
    names(IPredEpiClim) <- c("time", "model_predicted_absolute_cases_clim")
    IPredEpiClim$model_predicted_reported_cases_clim <- IPredEpiClim$model_predicted_absolute_cases_clim/fittedpars$rho[length(fittedpars$rho)]
    
    
    #Iall$reported_cases <- Iall$value/mean(fittedpars$rho)
    #head(Iall)
    #head(all.dat)
    
    #names(Iall) <- c("time", "model", "model_predicted_absolute_cases", "model_predicted_reported_cases")
    #and merge on time
    all.dat.merge <- merge(dat.pred, IPredEpi, by="time", all.x = T, sort=F)
    all.dat.merge.clim <- merge(dat.pred, IPredEpiClim, by="time", all.x = T, sort=F)
    
    
    #and get sum of sq differences
    all.dat.merge$sq_diff <- (all.dat.merge$cases - all.dat.merge$model_predicted_reported_cases)^2
    all.dat.merge.clim$sq_diff <- (all.dat.merge.clim$cases - all.dat.merge.clim$model_predicted_reported_cases)^2
    
    sm.sq[i] = sum(all.dat.merge$sq_diff)
    sm.sq.clim[i] = sum(all.dat.merge.clim$sq_diff)
    
    # ##saveRDS(predict_epi_Inc, paste0(homewd,"/sim.com.epi_predict_epi_Inc_",epiyr,".RDS") )
    # # 
    #  IPredEpi = melt(cbind.data.frame(time=predict_epi_Inc$I$time, mean_noInc=predict_epi_noInc$I$mean, mean_Inc=predict_epi_Inc$I$mean),id.vars = "time")
    #  ISimTS = cbind.data.frame(time=predict_ts$I$time, variable = rep("TSIR fit", length(predict_ts$I$time)), value=predict_ts$I$mean)
    # # 
    #  Iall <- rbind(IPredEpi, ISimTS)
    #  Iall$variable <- as.character(Iall$variable)
    #  Iall$variable[Iall$variable=="mean_Inc"] <- "prediction TSIR-increased S"
    #  Iall$variable[Iall$variable=="mean_noInc"] <- "prediction TSIR with climate"
    # # 
    #  all.dat <- subset(dat, time<=(epiyr+1))
    #  all.dat$cases <- mean(fittedpars$rho)*all.dat$cases
    #  ggplot(data=subset(Iall, variable!="prediction TSIR-increased S")) + 
    #    geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
    #    geom_line(aes(x=time, y=value, color=variable)) 
    # # 
    # # ggplot(data=Iall) + 
    # #   geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
    # #   geom_line(aes(x=time, y=value, color=variable)) 
    # # 
    # 
    # #capture the increased S that best recovers the peak
    # ratioImax[i]<-max(Iall$value[Iall$variable=="prediction TSIR-increased S"])/max(all.dat$cases)
  }
  #spot<-which(abs(ratioImax-1)==min(abs(ratioImax-1)))
  spot<-which(sm.sq==min(sm.sq))
  spot.clim<-which(sm.sq.clim==min(sm.sq.clim))
  fracincout <- fracincvec[spot]
  fracincout.clim <- fracincvec[spot.clim]
  #ratioImaxout<-ratioImax[spot]
  S_noinc<-Sepi
  
  dat.out <- cbind.data.frame(epidemic_year=epiyr, epiSnoic= S_noinc, frac_incS = fracincout, frac_incS_clim = fracincout.clim)
  
  return(dat.out)
}
#fits with an increased Beta
sim.com.beta <- function(time.start, dat, epiyr, family){
  #
  
  if(family=="poisson"){
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='poisson',link='log')
    ## The type of regression used in susceptible reconstruction. 
    # Options are 'gaussian', 'lm' (linear model), 'spline' (smooth.spline with 2.5 degrees freedom), 'lowess' (with f = 2/3, iter = 1), 'loess' (degree 1), and 'user' which is just a user inputed vector. Defaults to 'gaussian' and if that fails then defaults to loess.
    
    ## The family in the GLM regression. 
    # One can use any of the GLM ones, but the options are essentially 'poisson' (with link='log'), 'gaussian' (with link='log' or 'identity'), or 'quasipoisson' (with link='log'). Default is 'gaussian'.
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  }else if(family=="gaussian"){
    
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='gaussian')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              #method='pois',
                              nsim=100)
    
    
    
  } 
  
  
  Its=simfitted$res$mean*fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time=simfitted$res$time, I=Its)
  Ifinal<-Its[length(Its)]
  
  
  
  #now try to predict the epidemic both with and without the increased S
  #compute the most likely model fit with the proper increase which gives the 
  #closest match of model to observed cases
  
  fracincvec<-seq(1,5,0.01)
  #fracincdist<-matrix(0,1,1)
  #ratioImax<-matrix(0,1,length(fracincvec))
  #ratioImaxdist<-matrix(0,1,1)
  #Send<-matrix(0,1,1)
  sm.sq <- rep(NA, length(fracincvec))
  for(i in 1:length(fracincvec)){
    
    Sbegin=simfitted$sbar+simfitted$Z
    Sepi<-Sbegin[length(Sbegin)]
    #SepiInc<-fracincvec[i]*Sbegin[length(Sbegin)]
    
    
    #and look for the beta to multiply
    incBeta <- fracincvec[i]*fittedpars$beta
    #this is for no increased S
    dat.fit = subset(dat, time >= time.start & time<epiyr)
    
    predict_ts <- predicttsir(times=dat.fit$time,
                              births = dat.fit$births,
                              beta = fittedpars$beta,
                              alpha = fittedpars$alpha,
                              S0 = Sbegin[1],
                              I0 = dat.fit$cases[1],
                              nsim=100,
                              stochastic = T)
    
    
    
    #and the epi year-no increase
    dat.pred= subset(dat, time >= epiyr & time<=(epiyr+1))
    predict_epi_noInc <- predicttsir(times=dat.pred$time,
                                     births = dat.pred$births,
                                     beta = fittedpars$beta,
                                     alpha = fittedpars$alpha,
                                     S0 = Sepi,
                                     I0 = Ifinal,
                                     nsim=100,
                                     stochastic = T)
    
    
    
    #and including an increase on the beta
    predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                   births = dat.pred$births,
                                   beta = incBeta,
                                   alpha = fittedpars$alpha,
                                   S0 = Sepi,
                                   I0 = Ifinal,
                                   nsim=100,
                                   stochastic = T)
    
    
    #now look at data and model predictions for the epidemic year only
    IPredEpi = cbind.data.frame(time=predict_epi_Inc$I$time, mean_Inc=predict_epi_Inc$I$mean)
    names(IPredEpi) <- c("time", "model_predicted_absolute_cases")
    IPredEpi$model_predicted_reported_cases <- IPredEpi$model_predicted_absolute_cases/fittedpars$rho[length(fittedpars$rho)]

    #Iall$reported_cases <- Iall$value/mean(fittedpars$rho)
    #head(Iall)
    #head(all.dat)
    
    #names(Iall) <- c("time", "model", "model_predicted_absolute_cases", "model_predicted_reported_cases")
    #and merge on time
    all.dat.merge <- merge(dat.pred, IPredEpi, by="time", all.x = T, sort=F)
    
    
    #and get sum of sq differences
    all.dat.merge$sq_diff <- (all.dat.merge$cases - all.dat.merge$model_predicted_reported_cases)^2
    
    sm.sq[i] = sum(all.dat.merge$sq_diff)
    
  }
  spot<-which(sm.sq==min(sm.sq))
  fracincout <- fracincvec[spot]
  #ratioImaxout<-ratioImax[spot]
  S_noinc<-Sepi
  
  dat.out <- cbind.data.frame(epidemic_year=epiyr, epiSnoic= S_noinc, frac_incBeta = fracincout)
  
  return(dat.out)
}
#here, run fitted version with and without increased S
sim.with.increaseS <- function(fracincS,time.start, dat, epiyr, family){
  
  if(family=="poisson"){
    
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='poisson',link='log')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  }else if (family=="gaussian"){
    # 
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='gaussian')
    # 
    # #now simulate tsir
    
    # 
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              #method='pois', 
                              nsim=100)
  }
  #saveRDS(fittedpars, paste0(homewd,"/sim.with.increaseS_estpars_",epiyr,".RDS") )
  #saveRDS(simfitted, paste0(homewd,"/sim.with.increaseS_simfitted_",epiyr,".RDS") )
  #
  Its=simfitted$res$mean*fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time=simfitted$res$time, I=Its)
  Ifinal<-Its[length(Its)]
  
  #and sim with and without increased S
  
  
  #now try to predict the epidemic both with and without the increased S
  #compute the most likely model fit with the proper increase which gives the 
  #closest match of model to observed cases
  
  
  
  Sbegin=simfitted$sbar+simfitted$Z
  
  Sepi<-Sbegin[length(Sbegin)]
  SepiInc<-fracincS*Sepi
  
  #this is for no increased S
  dat.fit = subset(dat, time >= time.start & time<epiyr)
  
  predict_ts <- predicttsir(times=dat.fit$time,
                            births = dat.fit$births,
                            beta = fittedpars$beta,
                            alpha = fittedpars$alpha,
                            S0 = Sbegin[1],
                            I0 = dat.fit$cases[1],
                            nsim=100,
                            stochastic = T)
  
  
  
  #and the epi year-no increase
  dat.pred= subset(dat, time >= epiyr & time< (epiyr+1))
  
  predict_epi_noInc <- predicttsir(times=dat.pred$time,
                                   births = dat.pred$births,
                                   beta = fittedpars$beta,
                                   alpha = fittedpars$alpha,
                                   S0 = Sepi,
                                   I0 = Ifinal,
                                   nsim=100,
                                   stochastic = T)
  #saveRDS(predict_epi_noInc, paste0(homewd,"/sim.with.increaseS_predict_epi_noInc_",epiyr,".RDS") )
  
  #and including an increase
  predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                 births = dat.pred$births,
                                 beta = fittedpars$beta,
                                 alpha = fittedpars$alpha,
                                 S0 = SepiInc,
                                 I0 = Ifinal,
                                 nsim=100,
                                 stochastic = T)
  
  #saveRDS(predict_epi_Inc, paste0(homewd,"/sim.with.increaseS_predict_epi_Inc_",epiyr,".RDS") )
  
  
  IPredEpi = melt(cbind.data.frame(time=predict_epi_Inc$I$time, mean_noInc=predict_epi_noInc$I$mean, mean_Inc=predict_epi_Inc$I$mean),id.vars = "time")
  ISimTS = cbind.data.frame(time=predict_ts$I$time, variable = rep("TSIR fit", length(predict_ts$I$time)), value=predict_ts$I$mean)
  
  Iall <- rbind(IPredEpi, ISimTS)
  Iall$variable <- as.character(Iall$variable)
  Iall$variable[Iall$variable=="mean_Inc"] <- "prediction TSIR-increased S"
  Iall$variable[Iall$variable=="mean_noInc"] <- "prediction TSIR"
  
  all.dat <- subset(dat, time<=(epiyr+1))
  all.dat$cases <- mean(fittedpars$rho)*all.dat$cases
  # ggplot(data=subset(Iall, variable!="prediction TSIR-increased S")) + 
  #   geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
  #   geom_line(aes(x=time, y=value, color=variable)) 
  
  all.dat <- dplyr::select(all.dat, time, cases)
  all.dat$reported_cases <- all.dat$cases/mean(fittedpars$rho)
  all.dat$variable <- "data"
  
  names(Iall)[names(Iall)=="value"] <- "cases"
  Iall$reported_cases <- Iall$cases/mean(fittedpars$rho)
  Iall <- dplyr::select(Iall, names(all.dat))
  
  all.out <- rbind(Iall, all.dat)
  
  
  p2 <- ggplot(data=all.out) + 
    geom_line(aes(x=time, y=cases, color=variable)) 
  
  return(all.out)
  
}
sim.with.increaseBeta <- function(fracincBeta,time.start, dat, epiyr, family){
  
  if(family=="poisson"){
    
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='poisson',link='log')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  }else if (family=="gaussian"){
    # 
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='gaussian')
    # 
    # #now simulate tsir
    
    # 
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              #method='pois', 
                              nsim=100)
  }
  #saveRDS(fittedpars, paste0(homewd,"/sim.with.increaseS_estpars_",epiyr,".RDS") )
  #saveRDS(simfitted, paste0(homewd,"/sim.with.increaseS_simfitted_",epiyr,".RDS") )
  #
  Its=simfitted$res$mean*fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time=simfitted$res$time, I=Its)
  Ifinal<-Its[length(Its)]
  

  
  Sbegin=simfitted$sbar+simfitted$Z
  
  Sepi<-Sbegin[length(Sbegin)]
  #SepiInc<-fracincS*Sepi
  BetaInc <- fittedpars$beta*fracincBeta
  #this is for no increased S
  dat.fit = subset(dat, time >= time.start & time<epiyr)
  
  predict_ts <- predicttsir(times=dat.fit$time,
                            births = dat.fit$births,
                            beta = fittedpars$beta,
                            alpha = fittedpars$alpha,
                            S0 = Sbegin[1],
                            I0 = dat.fit$cases[1],
                            nsim=100,
                            stochastic = T)
  
  
  
  #and the epi year-no increase
  dat.pred= subset(dat, time >= epiyr & time< (epiyr+1))
  
  predict_epi_noInc <- predicttsir(times=dat.pred$time,
                                   births = dat.pred$births,
                                   beta = fittedpars$beta,
                                   alpha = fittedpars$alpha,
                                   S0 = Sepi,
                                   I0 = Ifinal,
                                   nsim=100,
                                   stochastic = T)
  #saveRDS(predict_epi_noInc, paste0(homewd,"/sim.with.increaseS_predict_epi_noInc_",epiyr,".RDS") )
  
  #and including an increase in beta
  predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                 births = dat.pred$births,
                                 beta = BetaInc,
                                 alpha = fittedpars$alpha,
                                 S0 = Sepi,
                                 I0 = Ifinal,
                                 nsim=100,
                                 stochastic = T)
  
  #saveRDS(predict_epi_Inc, paste0(homewd,"/sim.with.increaseS_predict_epi_Inc_",epiyr,".RDS") )
  
  
  IPredEpi = melt(cbind.data.frame(time=predict_epi_Inc$I$time, mean_noInc=predict_epi_noInc$I$mean, mean_Inc=predict_epi_Inc$I$mean),id.vars = "time")
  ISimTS = cbind.data.frame(time=predict_ts$I$time, variable = rep("TSIR fit", length(predict_ts$I$time)), value=predict_ts$I$mean)
  
  Iall <- rbind(IPredEpi, ISimTS)
  Iall$variable <- as.character(Iall$variable)
  Iall$variable[Iall$variable=="mean_Inc"] <- "prediction TSIR-increased Beta"
  Iall$variable[Iall$variable=="mean_noInc"] <- "prediction TSIR"
  
  all.dat <- subset(dat, time<=(epiyr+1))
  all.dat$cases <- mean(fittedpars$rho)*all.dat$cases
   # ggplot(data=subset(Iall, variable!="prediction TSIR-increased S")) + 
   #   geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
   #   geom_line(aes(x=time, y=value, color=variable)) 
   # 
  all.dat <- dplyr::select(all.dat, time, cases)
  all.dat$reported_cases <- all.dat$cases/mean(fittedpars$rho)
  all.dat$variable <- "data"
  
  names(Iall)[names(Iall)=="value"] <- "cases"
  Iall$reported_cases <- Iall$cases/mean(fittedpars$rho)
  Iall <- dplyr::select(Iall, names(all.dat))
  
  all.out <- rbind(Iall, all.dat)
  
  
  p2 <- ggplot(data=all.out) + 
    geom_line(aes(x=time, y=cases, color=variable)) 
  
  return(all.out)
  
}
sim.com.noS <- function(time.start, dat, epiyr, family){
  
  
  if(family=="poisson"){
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='poisson',link='log')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  }else if(family=="gaussian"){
    
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='gaussian')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              #method='pois',
                              nsim=100)
    
    
  }
  
  
  
  
  
  
  Its=simfitted$res$mean*fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time=simfitted$res$time, I=Its)
  Ifinal<-Its[length(Its)]
  
  
  
  Sbegin=simfitted$sbar+simfitted$Z
  Sepi<-Sbegin[length(Sbegin)]
  
  
  #this is for no increased S
  dat.fit = subset(dat, time >= time.start & time<epiyr)
  
  predict_ts <- predicttsir(times=dat.fit$time,
                            births = dat.fit$births,
                            beta = fittedpars$beta,
                            alpha = fittedpars$alpha,
                            S0 = Sbegin[1],
                            I0 = dat.fit$cases[1],
                            nsim=100,
                            stochastic = T)
  
  
  
  #and the epi year-no increase
  dat.pred= subset(dat, time >= epiyr & time<=(epiyr+1))
  predict_epi_noInc <- predicttsir(times=dat.pred$time,
                                   births = dat.pred$births,
                                   beta = fittedpars$beta,
                                   alpha = fittedpars$alpha,
                                   S0 = Sepi,
                                   I0 = Ifinal,
                                   nsim=100,
                                   stochastic = T)
  
  
  
  IPredEpi = melt(cbind.data.frame(time=predict_epi_noInc$I$time, mean_noInc=predict_epi_noInc$I$mean),id.vars = "time")
  ISimTS = cbind.data.frame(time=predict_ts$I$time, variable = rep("TSIR fit", length(predict_ts$I$time)), value=predict_ts$I$mean)
  
  Iall <- rbind(IPredEpi, ISimTS)
  Iall$variable <- as.character(Iall$variable)
  #Iall$variable[Iall$variable=="mean_Inc"] <- "prediction TSIR-increased S"
  Iall$variable[Iall$variable=="mean_noInc"] <- "prediction TSIR"
  
  all.dat <- subset(dat, time<=(epiyr+1))
  all.dat$cases <- mean(fittedpars$rho)*all.dat$cases
  # ggplot(data=subset(Iall, variable!="prediction TSIR-increased S")) + 
  #   geom_line(data=all.dat, aes(x=time, y=cases),linetype=2) +
  #   geom_line(aes(x=time, y=value, color=variable)) 
  
  
  all.dat <- dplyr::select(all.dat, time, cases)
  all.dat$reported_cases <- all.dat$cases/mean(fittedpars$rho)
  all.dat$variable <- "data"
  
  names(Iall)[names(Iall)=="value"] <- "cases"
  Iall$reported_cases <- Iall$cases/mean(fittedpars$rho)
  Iall <- dplyr::select(Iall, names(all.dat))
  
  all.out <- rbind(Iall, all.dat)
  
  
  p2 <- ggplot(data=all.out) + 
    geom_line(aes(x=time, y=cases, color=variable)) 
  
  return(all.out)
  
  
  
}
sim.return.S.incS <- function(fracincS,time.start, dat, epiyr, family){
  
  if(family=="poisson"){
    
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='poisson',link='log')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  }else if (family=="gaussian"){
    # 
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='gaussian')
    # 
    # #now simulate tsir
    
    # 
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              #method='pois', 
                              nsim=100)
  }
  # 
  Its=simfitted$res$mean*fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time=simfitted$res$time, I=Its)
  Ifinal<-Its[length(Its)]
  
  #and sim with and without increased S
  
  
  #now try to predict the epidemic both with and without the increased S
  #compute the most likely model fit with the proper increase which gives the 
  #closest match of model to observed cases
  
  
  Sbegin=simfitted$sbar+simfitted$Z
  Sepi<-Sbegin[length(Sbegin)]
  SepiInc<-fracincS*Sepi
  #betabegin=simfitted$beta
  #betaepi<-betabegin[length(betabegin)]
  #betaepiInc<-fracincS*betaepi
  
  #this is for no increased S
  dat.fit = subset(dat, time >= time.start & time<epiyr)
  
  predict_ts <- predicttsir(times=dat.fit$time,
                            births = dat.fit$births,
                            beta = fittedpars$beta,
                            alpha = fittedpars$alpha,
                            S0 = Sbegin[1],
                            I0 = dat.fit$cases[1],
                            nsim=100,
                            stochastic = T)
  
  
  #now, return the data frame of susceptibles vs births
  
  dat.S <- cbind.data.frame(time = dat.fit$time, births=dat.fit$births, sus_mean = predict_ts$S$mean, sus_low= predict_ts$S$low, sus_high=predict_ts$S$high)
  dat.S$variable <- "TSIR fit"
  #head(dat.S)
  
  
  #and the epi year-no increase
  dat.pred= subset(dat, time >= epiyr & time< (epiyr+1))
  predict_epi_noInc <- predicttsir(times=dat.pred$time,
                                   births = dat.pred$births,
                                   beta = fittedpars$beta,
                                   alpha = fittedpars$alpha,
                                   S0 = Sepi,
                                   I0 = Ifinal,
                                   nsim=100,
                                   stochastic = T)
  
  #and including an increase
  predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                 births = dat.pred$births,
                                 beta = fittedpars$beta,
                                 alpha = fittedpars$alpha,
                                 S0 = SepiInc,
                                 I0 = Ifinal,
                                 nsim=100,
                                 stochastic = T)
  
  #now combine the S data
  dat.S.noinc <- cbind.data.frame(time = dat.pred$time, births=dat.pred$births, sus_mean = predict_epi_noInc$S$mean, sus_low= predict_epi_noInc$S$low, sus_high=predict_epi_noInc$S$high)
  dat.S.noinc$variable <- "prediction TSIR"
  
  dat.S.inc <- cbind.data.frame(time = dat.pred$time, births=dat.pred$births, sus_mean = predict_epi_Inc$S$mean, sus_low= predict_epi_Inc$S$low, sus_high=predict_epi_Inc$S$high)
  dat.S.inc$variable <- "prediction TSIR-increased S"
  
  dat.S.all <- rbind(dat.S, dat.S.noinc, dat.S.inc)
  dat.S.all$year <- trunc(dat.S.all$time)
  dat.S.yr <- dlply(dat.S.all, .(year))
  
  take.min <- function(df){
    df1 <- subset(df, time == min(df$time))
    return(df1)
  }
  
  #then, take only the even year data
  dat.S.yr <- data.table::rbindlist(lapply(dat.S.yr, take.min))
  
  #and return this
  return(dat.S.yr)
}
sim.return.S.noinc <- function(fracincS,time.start, dat, epiyr, family){
  
  if(family=="poisson"){
    
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='poisson',link='log')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  }else if (family=="gaussian"){
    # 
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='gaussian')
    # 
    # #now simulate tsir
    
    # 
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              #method='pois', 
                              nsim=100)
  }
  # 
  Its=simfitted$res$mean*fittedpars$rho #account for underreporting
  dat.sim <- cbind.data.frame(time=simfitted$res$time, I=Its)
  Ifinal<-Its[length(Its)]
  
  #and sim with and without increased S
  
  
  #now try to predict the epidemic both with and without the increased S
  #compute the most likely model fit with the proper increase which gives the 
  #closest match of model to observed cases
  
  
  Sbegin=simfitted$sbar+simfitted$Z
  Sepi<-Sbegin[length(Sbegin)]
  #SepiInc<-fracincS*Sepi
  
  
  #this is for no increased S
  dat.fit = subset(dat, time >= time.start & time<epiyr)
  
  predict_ts <- predicttsir(times=dat.fit$time,
                            births = dat.fit$births,
                            beta = fittedpars$beta,
                            alpha = fittedpars$alpha,
                            S0 = Sbegin[1],
                            I0 = dat.fit$cases[1],
                            nsim=100,
                            stochastic = T)
  
  
  #now, return the data frame of susceptibles vs births
  
  dat.S <- cbind.data.frame(time = dat.fit$time, births=dat.fit$births, sus_mean = predict_ts$S$mean, sus_low= predict_ts$S$low, sus_high=predict_ts$S$high)
  dat.S$variable <- "TSIR fit"
  #head(dat.S)
  
  
  #and the epi year-no increase
  dat.pred= subset(dat, time >= epiyr & time< (epiyr+1))
  predict_epi_noInc <- predicttsir(times=dat.pred$time,
                                   births = dat.pred$births,
                                   beta = fittedpars$beta,
                                   alpha = fittedpars$alpha,
                                   S0 = Sepi,
                                   I0 = Ifinal,
                                   nsim=100,
                                   stochastic = T)
  # 
  # #and including an increase
  # predict_epi_Inc <- predicttsir(times=dat.pred$time,
  #                                births = dat.pred$births,
  #                                beta = fittedpars$beta,
  #                                alpha = fittedpars$alpha,
  #                                S0 = SepiInc,
  #                                I0 = Ifinal,
  #                                nsim=100,
  #                                stochastic = T)
  # 
  #now combine the S data
  dat.S.noinc <- cbind.data.frame(time = dat.pred$time, births=dat.pred$births, sus_mean = predict_epi_noInc$S$mean, sus_low= predict_epi_noInc$S$low, sus_high=predict_epi_noInc$S$high)
  dat.S.noinc$variable <- "prediction TSIR"
  # 
  # dat.S.inc <- cbind.data.frame(time = dat.pred$time, births=dat.pred$births, sus_mean = predict_epi_Inc$S$mean, sus_low= predict_epi_Inc$S$low, sus_high=predict_epi_Inc$S$high)
  # dat.S.inc$variable <- "prediction TSIR-increased S"
  
  dat.S.all <- rbind(dat.S, dat.S.noinc)#, dat.S.inc)
  dat.S.all$year <- trunc(dat.S.all$time)
  dat.S.yr <- dlply(dat.S.all, .(year))
  
  take.min <- function(df){
    df1 <- subset(df, time == min(df$time))
    return(df1)
  }
  
  #then, take only the even year data
  dat.S.yr <- data.table::rbindlist(lapply(dat.S.yr, take.min))
  
  #and return this
  return(dat.S.yr)
}
plot.comp <- function(dat, filename){
  dat$variable = factor(dat$variable, levels=c("data",
                                               "TSIR fit",
                                               "prediction TSIR",
                                               "prediction TSIR-increased S"))
  colz = c('data' = "black", 'TSIR fit' = "blue", 'prediction TSIR' = "red", 'prediction TSIR-increased S' = "seagreen")
  typez = c('data' = 2, 'TSIR fit' = 1, 'prediction TSIR' = 1, 'prediction TSIR-increased S' = 1)
  p1 <- ggplot(data=dat) + 
    geom_line(aes(x=time, y=reported_cases, color=variable, linetype=variable), size=1) +
    scale_color_manual(values=colz) +
    scale_linetype_manual(values=typez) +theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=18),
          legend.text = element_text(size=10),
          legend.position = c(.25,.82),
          axis.text = element_text(size=14), legend.title = element_blank()) +
    ylab("reported cases")
  
  print(p1)
  
  ggsave(file =filename,
         plot=p1,
         units="mm",  
         width=50, 
         height=40, 
         scale=3, 
         dpi=300)
  
  
}


###################################################################################
###################################################################################
###################################################################################
###################################################################################


#load data into tsir form
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_data.csv"), header = T, stringsAsFactors = F)

head(tsir.dat)

#also load the climate data
clim.dat <- read.csv(file = paste0(homewd, "/data/clim_biwk_nat.csv"), header = T, stringsAsFactors = F)
head(clim.dat)

#just as a gut check, plot mean temp and sum precip by year across the time series
clim.sum = ddply(clim.dat, .(year), summarise, mean_temp = mean(temp_C), sum_precip=sum(precip_mm))

ggplot(data=clim.sum) + geom_line(aes(x=year, y=mean_temp, color=year))
ggplot(data=clim.sum) + geom_line(aes(x=year, y=sum_precip, color=year))

names(tsir.dat)[names(tsir.dat)=="biweek"] <- "biwk"
names(tsir.dat)[names(tsir.dat)=="yr"] <- "year"

#merge climate data
tsir.dat <- merge(tsir.dat, clim.dat, by =c("year", "biwk") , all.x = T)
head(tsir.dat)

tsir.dat <- arrange(tsir.dat, time)

nat.lag <- sim.com.epi.clim.nat.lag(dat = tsir.dat,
                                    time.start1 =  2002,
                                    time.start2 =  2008,
                                    time.start3 =  2013,
                                    family="gaussian",
                                    epiyr1 = 2007,
                                    epiyr2 = 2012,
                                    epiyr3 = 2019)

#now feed these lags in to panel regression and prediction of epidemic years

#fit to pre epidemic period with increased S
#addig climate demands even more susceptibles! this year was a cold one!
sim.2007 <- sim.com.epi.clim(dat = tsir.dat,
                        time.start =  min(tsir.dat$time),
                        lag.dat =nat.lag,
                        family="gaussian",
                        epiyr = 2007)

#do all the pre-pandemic fits

#now for 2012 - fit with increased S
#adding climate again demands more susceptibles - 
#climatic conditions were suboptimal in this year
sim.2012 <- sim.com.epi.clim(dat = tsir.dat,
                             time.start =  min(tsir.dat$time[tsir.dat$time>=2008]),
                             lag.dat =nat.lag,
                             family="gaussian",
                             epiyr = 2012)

#and fit 2019 with increased S - 
#climate accounts for a bit of difference here. not so in the other epidemic years
sim.2019 <- sim.com.epi.clim(dat = tsir.dat,
                             time.start =  min(tsir.dat$time[tsir.dat$time>=2013]),
                             lag.dat =nat.lag,
                             family="gaussian",
                             epiyr = 2019)


#and fit with increased beta
sim.beta.2007 <- sim.com.beta(dat = tsir.dat,
                        time.start =  min(tsir.dat$time),
                        family="gaussian",
                        epiyr = 2007)

#now return and recover Sdat for increased S
Sdat2007 <- sim.return.S.incS(fracincS = sim.2007$frac_incS,
                              dat=tsir.dat,
                              time.start =  min(tsir.dat$time),
                              family="gaussian",
                              epiyr = 2007)


#and sim with increased S
comp.2007 <- sim.with.increaseS(
  fracincS = sim.2007$frac_incS,
  dat = tsir.dat,
  time.start =  min(tsir.dat$time),
  family = "gaussian",
  epiyr = 2007
)


#and sim with increased beta
comp.beta.2007 <- sim.with.increaseBeta(
  fracincBeta = sim.beta.2007$frac_incBeta,
  dat = tsir.dat,
  time.start =  min(tsir.dat$time),
  family = "gaussian",
  epiyr = 2007
)









#and with increased beta
sim.beta.2012 <- sim.com.beta(dat = tsir.dat,
                        time.start =  min(tsir.dat$time[tsir.dat$time>=2008]),
                        family="gaussian",
                        epiyr = 2012)


#now return S from increased S
Sdat2012 <- sim.return.S.incS(fracincS = sim.2012$frac_incS,
                              dat=subset(tsir.dat, time >=min(tsir.dat$time[tsir.dat$time>=2008]) & time<2013),
                              time.start =  min(tsir.dat$time[tsir.dat$time>=2008]),
                              family="gaussian",
                              epiyr = 2012)


#now sim with increased S
comp.2012 <- sim.with.increaseS(
  fracincS = sim.2012$frac_incS,
  dat = subset(tsir.dat, time >= min(tsir.dat$time[tsir.dat$time >=
                                                     2008]) & time < 2013),
  time.start =  min(tsir.dat$time[tsir.dat$time >=
                                    2008]),
  family = "gaussian",
  epiyr = 2012
)

#and sim with increased beta
comp.beta.2012  <- sim.with.increaseBeta(
  fracincBeta = sim.beta.2012$frac_incBeta,
  dat = subset(tsir.dat, time >= min(tsir.dat$time[tsir.dat$time >=
                                                     2008]) & time < 2013),
  time.start =  min(tsir.dat$time[tsir.dat$time >=
                                    2008]),
  family = "gaussian",
  epiyr = 2012
)








#and with increased beta
sim.beta.2019 <- sim.com.beta(dat = tsir.dat,
                        time.start =  min(tsir.dat$time[tsir.dat$time>=2013]),
                        family="gaussian",
                        epiyr = 2019)


#now return S from increased S
Sdat2019 <- sim.return.S.incS(fracincS = sim.2019$frac_incS,
                              dat=subset(tsir.dat, time >=min(tsir.dat$time[tsir.dat$time>=2013]) &time<2020),
                              time.start =  min(tsir.dat$time[tsir.dat$time>=2013]),
                              family="gaussian",
                              epiyr = 2019)


#now return and simulate with results
comp.2019 <- sim.with.increaseS(
  fracincS = sim.2019$frac_incS,
  dat = subset(tsir.dat, time >= min(tsir.dat$time[tsir.dat$time >=
                                                     2013]) & time < 2020),
  time.start =  min(tsir.dat$time[tsir.dat$time >=
                                    2013]),
  family = "gaussian",
  epiyr = 2019
)


#and sim with increased beta
comp.beta.2019  <- sim.with.increaseBeta(
  fracincBeta = sim.beta.2019$frac_incBeta,
  dat = subset(tsir.dat, time >= min(tsir.dat$time[tsir.dat$time >=
                                                     2013]) & time < 2020),
  time.start =  min(tsir.dat$time[tsir.dat$time >=
                                    2013]),
  family = "gaussian",
  epiyr = 2019
)





#and combine for the susceptible plot
SdatIncSim <- rbind(Sdat2007, Sdat2012, Sdat2019)

SdatIncSim$sim <- "TSIR-increased-S"
 
#and without the breaks
###and the whole time series
#and the results for 2019

sim.all <- sim.com.epi(dat = tsir.dat,
                       time.start =  min(tsir.dat$time[tsir.dat$time>=2002]),
                       family="gaussian",
                       epiyr = 2019)

Sdat2019fullsim <- sim.return.S.noinc(fracincS = sim.all$frac_incS,
                                      dat=subset(tsir.dat, time>=2002),
                                      time.start =  min(tsir.dat$time[tsir.dat$time>=2002]),
                                      family="gaussian",
                                      epiyr = 2019)

Sdat2019fullsim$sim <- "TSIR"

#and combine and plot

Sdatcombined <- rbind(SdatIncSim, Sdat2019fullsim)
Sdatcombined$epiyr <- 0
Sdatcombined$epiyr[Sdatcombined$year==2007|Sdatcombined$year==2012|Sdatcombined$year==2019] <- 1
Sdatcombined$epiyr <- factor(Sdatcombined$epiyr)

head(Sdatcombined)
unique(Sdatcombined$sim)

#and write data file
#write.csv(Sdatcombined, file = paste0(homewd, "/SdatCara.csv"), row.names = F)
saveRDS(Sdatcombined,file = paste0(homewd, "/data/Sdatcombined.rds"))



##################################################################
####################################################################
#####################################################################
######################################################################

#and save the sims of increased S and not for the joint dataset
all.comp <- rbind(comp.2007, comp.2012, comp.2019)

#and combine beta sims and add to dataset
comp.beta <- rbind(comp.beta.2007, comp.beta.2012, comp.beta.2019)
comp.beta = subset(comp.beta, variable == "prediction TSIR-increased Beta")

#and combine
all.comp <- rbind(all.comp, comp.beta)

#and save all comp as data
write.csv(all.comp, file = paste0(homewd, "/data/TSIR_fitted_timeseries_estimates.csv"), row.names = F)

my_blue<-"#2980B9"
my_orange<-"#D35400"
my_red<-"tomato"
my_green<-"seagreen"
my_yellow<-"#F1C40F"
my_purple <- "purple"



#now plot the whole thing
# all.IncBeta <- rbind(comp.2007, comp.2012, comp.2019)
# all.IncBeta <- subset(all.IncBeta, variable == "prediction TSIR-increased Beta")
# all.IncBeta$epiyr <- 1
# 
# #and join
# Sdatcombined <- rbind(Sdatcombined, all.IncBeta)
#and plot all
dat=all.comp
dat$variable = factor(dat$variable, levels=c("data",
                                             "TSIR fit",
                                             "prediction TSIR",
                                             "prediction TSIR-increased S",
                                             "prediction TSIR-increased Beta"))
colz = c('data' = "black", 'TSIR fit' = my_blue, 'prediction TSIR' = my_green, 'prediction TSIR-increased S' =my_red, 'prediction TSIR-increased Beta' = my_purple)
typez = c('data' = 2, 'TSIR fit' = 1, 'prediction TSIR' = 1, 'prediction TSIR-increased S' = 1, 'prediction TSIR-increased Beta' = 1)


pl_a <- ggplot(data=dat, aes(x=time, y=reported_cases, color=variable, linetype=variable)) + 
  geom_line(data=dplyr::filter(dat, variable=="data"), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR fit" & time<2007), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR fit" & time>=2008 & time<2012), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR fit" & time>=2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR" & time>=2007 & time<2008), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR" & time>=2012 & time<2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR" & time>=2019 & time<2020), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased S" & time>=2007 & time<2008), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased Beta" & time>=2007 & time<2008), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased S" & time>=2012 & time<2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased Beta" & time>=2012 & time<2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased S" & time>=2019 & time<2020), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction TSIR-increased Beta" & time>=2019 & time<2020), size=1) +
  scale_color_manual(values=colz) + # coord_cartesian(ylim = c(0,1500))+
  scale_linetype_manual(values=typez) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        legend.text = element_text(size=10),
        legend.position = c(.1,.82),
        axis.text = element_text(size=14), legend.title = element_blank()) +
  ylab(paste0("reported cases"))

plot(pl_a)











