rm(list=ls())

require(tsiR)
require(kernlab)
require(ggplot2)
require(reshape2)
require(grid)
require(plyr)



# as one script, split by province, pull in and reconstruct susceptibles by
# best method, fit TSIR, and run future year with climate projected beta
# that you have already saved elsewhere


homewd = "/Users/carabrook/Developer/cambodia-dengue-national-public"
setwd(paste0(homewd))

#load the functions
###################################################################################

#this function finds the optimal fraction of increased susceptibles both
#with and without the climate-driven beta
find.frac.incS <- function(simfitted1, fittedpars1, dat1, epiyr1, sbar1, clim.beta){
  
  
  
  Its=simfitted1$res$mean*fittedpars1$rho #account for underreporting
  Ifinal<-Its[length(Its)]
  
  dat.sim <- cbind.data.frame(time=simfitted1$res$time, I=Its)
  #now try to predict the epidemic both with and without the increased S
  #compute the most likely model fit with the proper increase which gives the 
  #closest match of model to observed cases
  
  #make a vector of the increases to tested
  
  fracincvec <- seq(1,3,0.01)
  fracincvec.list <- as.list(fracincvec)
  
  #and apply this function over that vector to get deviations from the actual data
  out.inc.S <- lapply(fracincvec.list, run.increased.S, data.series=dat1, epiyr1=epiyr1,  simfitted1=simfitted1, fittedpars1=fittedpars1, beta.clim=clim.beta$beta)
  
  #now, fit the increase that has the lowest sum of squares for both the standard beta and the climate-driven beta
  sm.sq.noclim <- c(unlist(sapply(out.inc.S, '[[', 1)))
  sm.sq.clim <- c(unlist(sapply(out.inc.S, '[[', 2)))
  
  spot.noclim<-which(sm.sq.noclim==min(sm.sq.noclim))
  spot.clim<-which(sm.sq.clim==min(sm.sq.clim))
  fracincout.noclim <- fracincvec[spot.noclim]
  fracincout.clim <- fracincvec[spot.clim]
  
  Sbegin=simfitted1$sbar+simfitted1$Z
  S_noinc<-Sbegin[length(Sbegin)]
  
  
  dat.out <- cbind.data.frame(epidemic_year=epiyr1, epiSnoic= S_noinc, frac_incS_noclim = fracincout.noclim, frac_incS_clim=fracincout.clim)
  
  return(dat.out)
  
  
  
  
}
run.increased.S <- function(fracincS, data.series, epiyr1, simfitted1, fittedpars1, beta.clim){
  
  Its=simfitted1$res$mean*fittedpars1$rho #account for under-reporting. 
  # This is a single number in the case of years where the best S reconstruction
  # was with a linear regression or a vector the length of the time series in the
  # case of years where the best S reconstruction was with a Gaussian regression
  
  Ifinal<-Its[length(Its)]
  
  time.start =  min(data.series$time)
  
  #the fitting data
  dat.fit = subset(data.series, time >= time.start & time<epiyr1)
  #the prediction data
  dat.pred= subset(data.series, time >= epiyr1 & time<=(epiyr1+1))
  
  
  
  Sbegin=simfitted1$sbar+simfitted1$Z
  Sepi<-Sbegin[length(Sbegin)]
  SepiInc<-fracincS*Sbegin[length(Sbegin)]
  
  #this is for no increased S
  
  
  # predict_ts <- predicttsir(times=dat.fit$time,
  #                           births = dat.fit$births,
  #                           beta = fittedpars1$beta,
  #                           alpha = fittedpars1$alpha,
  #                           S0 = Sbegin[1],
  #                           I0 = dat.fit$cases[1],
  #                           nsim=100,
  #                           stochastic = T)
  # 
  
  
  # and the epi prediction year-no increased S 
  # and using the same beta as the rest of the year
  
  
  
  # predict_epi_noInc <- predicttsir(times=dat.pred$time,
  #                                  births = dat.pred$births,
  #                                  beta = fittedpars1$beta,
  #                                  alpha = fittedpars1$alpha,
  #                                  S0 = Sepi,
  #                                  I0 = Ifinal,
  #                                  nsim=100,
  #                                  stochastic = T)
  # 
  
  # and epi-year prediction including increased S
  # and using the same beta as the rest of the year
  predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                 births = dat.pred$births,
                                 beta = fittedpars1$beta,
                                 alpha = fittedpars1$alpha,
                                 S0 = SepiInc,
                                 I0 = Ifinal,
                                 nsim=100,
                                 stochastic = T)
  
  
  # # and epi-year prediction with climate-driven beta but no increased S
  # predict_epi_noInc_clim <- predicttsir(times=dat.pred$time,
  #                                       births = dat.pred$births,
  #                                       beta = beta.clim,
  #                                       alpha = fittedpars1$alpha,
  #                                       S0 = Sepi,
  #                                       I0 = Ifinal,
  #                                       nsim=100,
  #                                       stochastic = T)
  # 
  
  
  #and the epiyear prediction with climate-driven beta AND increased S
  predict_epi_Inc_clim <- predicttsir(times=dat.pred$time,
                                      births = dat.pred$births,
                                      beta = beta.clim,
                                      alpha = fittedpars1$alpha,
                                      S0 = SepiInc,
                                      I0 = Ifinal,
                                      nsim=100,
                                      stochastic = T)
  
  
  
  
  #now look at data and model predictions for the epidemic year only across all the possibilities
  
  #here just TSIR
  
  #  IPredEpi1 = cbind.data.frame(time=predict_epi_noInc$I$time, mean_Inc=predict_epi_noInc$I$mean)
  #  names(IPredEpi1) <- c("time", "model_predicted_absolute_cases")
  #  IPredEpi1$model_predicted_reported_cases <- IPredEpi1$model_predicted_absolute_cases/fittedpars$rho[length(fittedpars$rho)]
  #  IPredEpi1$sim_type <- "TSIR fit"
  # # 
  # # here with increased S
  IPredEpi = cbind.data.frame(time=predict_epi_Inc$I$time, mean_Inc=predict_epi_Inc$I$mean)
  names(IPredEpi) <- c("time", "no_clim_model_predicted_absolute_cases")
  IPredEpi$no_clim_model_predicted_reported_cases <- IPredEpi$no_clim_model_predicted_absolute_cases/mean(fittedpars1$rho)
  
  
  #here with increased S and climate-driven beta
  IPredEpi2 = cbind.data.frame(time=predict_epi_Inc_clim$I$time, mean_Inc=predict_epi_Inc_clim$I$mean)
  names(IPredEpi2) <- c("time", "clim_model_predicted_absolute_cases")
  IPredEpi2$clim_model_predicted_reported_cases <- IPredEpi2$clim_model_predicted_absolute_cases/mean(fittedpars1$rho)
  
  
  
  
  #Iall$reported_cases <- Iall$value/mean(fittedpars$rho)
  #head(Iall)
  #head(all.dat)
  
  #names(Iall) <- c("time", "model", "model_predicted_absolute_cases", "model_predicted_reported_cases")
  #and merge on time
  all.dat.merge <- merge(dat.pred, IPredEpi, by="time", all.x = T, sort=F)
  all.dat.merge <- merge(all.dat.merge, IPredEpi2, by="time", all.x = T, sort=F)
  
  # ggplot(data=all.dat.merge) + 
  #   geom_line(aes(x=time, y=cases),linetype=2) +
  #   geom_line(aes(x=time, y=model_predicted_reported_cases), color="tomato") 
  # # 
  # 
  
  #and get sum of sq differences for increased S
  all.dat.merge$sq_diff_noclim <- (all.dat.merge$cases - all.dat.merge$no_clim_model_predicted_reported_cases)^2
  all.dat.merge$sq_diff_clim <- (all.dat.merge$cases - all.dat.merge$clim_model_predicted_reported_cases)^2
  
  #and the sum of squared differences for increased S from the climate driven beta
  
  sm.sq.noclim = sum(all.dat.merge$sq_diff_noclim)
  sm.sq.clim = sum(all.dat.merge$sq_diff_clim)
  
  return(list(sm.sq.noclim, sm.sq.clim))
  
}
wrap.pipeline.TSIR <- function(df, epiyr1, sbar1, sus.dat, epi.beta.df){
  print(unique(df$provname))
  #first, identify the right susceptible reconstruction
  prov.choose = subset(sus.dat, provname == unique(df$provname) & epiyr == epiyr1)
  
  #and pick the right epiyear beta
  beta.epi = subset(epi.beta.df, provname == unique(df$provname) & epiyr == epiyr1)
  beta.epi <- arrange(beta.epi, time)
  
  suffix = unique(df$provname)
  #print(suffix)
  suffix <- gsub(pattern=" ", replacement = "_", x=suffix)
  
  time.start =  min(df$time)
  
  dat.fit = subset(df, time >= time.start & time<epiyr1)
  
  
  if(prov.choose$sus_reconstruction=="lm"){
    fittedpars <- estpars(data=dat.fit,
                          IP=2, 
                          alpha=NULL, 
                          sbar=sbar1,
                          xreg = "cumcases",
                          regtype='lm',
                          family='poisson',
                          link='log')
    
    
    simfitted <- simulatetsir(data=dat.fit,
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  }else if(prov.choose$sus_reconstruction=="gaussian"){
    
    fittedpars <- estpars(data=dat.fit,
                          IP=2, 
                          alpha=NULL, 
                          sbar=sbar1, 
                          xreg = "cumcases",
                          regtype='gaussian',
                          family='poisson',
                          link='log')
    
    
    #p1 <- plotsbar(fittedpars2)
    
    simfitted <- simulatetsir(data=dat.fit,
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  } 
  
  
  #now, take above and find the optimal increased S both with and without climate-driven beta
  incS.df <- find.frac.incS(simfitted1 = simfitted, fittedpars1 = fittedpars, dat1 = df, epiyr1 = epiyr1, sbar1 = sbar1, clim.beta = beta.epi)
  incS.df$diff_S_clim <- incS.df$frac_incS_noclim - incS.df$frac_incS_clim
  incS.df$provname <- suffix
  
 
  # Run the fitted data series
  
  Its=simfitted$res$mean*fittedpars$rho #account for underreporting
  Ifinal<-Its[length(Its)]
  
  Sbegin=simfitted$sbar+simfitted$Z
  Sepi<-Sbegin[length(Sbegin)]
  SepiIncNoClim<-Sbegin[length(Sbegin)]*incS.df$frac_incS_noclim
  SepiIncClim<-Sbegin[length(Sbegin)]*incS.df$frac_incS_clim
  
  
  
   predict_ts <- predicttsir(times=dat.fit$time,
                             births = dat.fit$births,
                             beta = fittedpars$beta,
                             alpha = fittedpars$alpha,
                             S0 = Sbegin[1],
                             I0 = dat.fit$cases[1],
                             nsim=100,
                             stochastic = T)
   
   
   predict_ts_lci <- predicttsir(times=dat.fit$time,
                             births = dat.fit$births,
                             beta = fittedpars$contact$betalow,
                             alpha = fittedpars$alpha,
                             S0 = Sbegin[1],
                             I0 = dat.fit$cases[1],
                             nsim=100,
                             stochastic = T)
   
   predict_ts_uci <- predicttsir(times=dat.fit$time,
                                 births = dat.fit$births,
                                 beta = fittedpars$contact$betahigh,
                                 alpha = fittedpars$alpha,
                                 S0 = Sbegin[1],
                                 I0 = dat.fit$cases[1],
                                 nsim=100,
                                 stochastic = T)
   
  
   #now run  all 4 time series for the epidemic year
   
   #prediction data
   dat.pred= subset(df, time >= epiyr1 & time<=(epiyr1+1))
   
   #1. project TSIR fit (no increased S + normal beta)
   
   predict_epi_noInc <- predicttsir(times=dat.pred$time,
                                    births = dat.pred$births,
                                    beta = fittedpars$beta,
                                    alpha = fittedpars$alpha,
                                    S0 = Sepi,
                                    I0 = Ifinal,
                                    nsim=100,
                                    stochastic = T)
   
   
   
   predict_epi_noInc_lci <- predicttsir(times=dat.pred$time,
                                    births = dat.pred$births,
                                    beta = fittedpars$contact$betalow,
                                    alpha = fittedpars$alpha,
                                    S0 = Sepi,
                                    I0 = Ifinal,
                                    nsim=100,
                                    stochastic = T)
   
   predict_epi_noInc_uci <- predicttsir(times=dat.pred$time,
                                        births = dat.pred$births,
                                        beta = fittedpars$contact$betahigh,
                                        alpha = fittedpars$alpha,
                                        S0 = Sepi,
                                        I0 = Ifinal,
                                        nsim=100,
                                        stochastic = T)
  
   
  #2. project epi year with increased S and normal beta
   
  # and epi-year prediction including increased S
  # and using the same beta as the rest of the year
  predict_epi_Inc <- predicttsir(times=dat.pred$time,
                                 births = dat.pred$births,
                                 beta = fittedpars$beta,
                                 alpha = fittedpars$alpha,
                                 S0 = SepiIncNoClim,
                                 I0 = Ifinal,
                                 nsim=100,
                                 stochastic = T)
  
  predict_epi_Inc_lci <- predicttsir(times=dat.pred$time,
                                 births = dat.pred$births,
                                 beta = fittedpars$contact$betalow,
                                 alpha = fittedpars$alpha,
                                 S0 = SepiIncNoClim,
                                 I0 = Ifinal,
                                 nsim=100,
                                 stochastic = T)
  
  
  predict_epi_Inc_uci <- predicttsir(times=dat.pred$time,
                                     births = dat.pred$births,
                                     beta = fittedpars$contact$betahigh,
                                     alpha = fittedpars$alpha,
                                     S0 = SepiIncNoClim,
                                     I0 = Ifinal,
                                     nsim=100,
                                     stochastic = T)
  
  
  
  #3. project epi year with no increased S but climate-driven beta
  
  
  
   predict_epi_noInc_clim <- predicttsir(times=dat.pred$time,
                                         births = dat.pred$births,
                                         beta = beta.epi$beta,
                                         alpha = fittedpars$alpha,
                                         S0 = Sepi,
                                         I0 = Ifinal,
                                         nsim=100,
                                         stochastic = T)
   
   
   
   predict_epi_noInc_clim_lci <- predicttsir(times=dat.pred$time,
                                         births = dat.pred$births,
                                         beta = beta.epi$betalow,
                                         alpha = fittedpars$alpha,
                                         S0 = Sepi,
                                         I0 = Ifinal,
                                         nsim=100,
                                         stochastic = T)
   
   predict_epi_noInc_clim_uci <- predicttsir(times=dat.pred$time,
                                             births = dat.pred$births,
                                             beta = beta.epi$betahigh,
                                             alpha = fittedpars$alpha,
                                             S0 = Sepi,
                                             I0 = Ifinal,
                                             nsim=100,
                                             stochastic = T)
  
  
  #4. project epi year with increased S and climate-driven beta
   
  predict_epi_Inc_clim <- predicttsir(times=dat.pred$time,
                                      births = dat.pred$births,
                                      beta = beta.epi$beta,
                                      alpha = fittedpars$alpha,
                                      S0 = SepiIncClim,
                                      I0 = Ifinal,
                                      nsim=100,
                                      stochastic = T)
  
  predict_epi_Inc_clim_lci <- predicttsir(times=dat.pred$time,
                                      births = dat.pred$births,
                                      beta = beta.epi$betalow,
                                      alpha = fittedpars$alpha,
                                      S0 = SepiIncClim,
                                      I0 = Ifinal,
                                      nsim=100,
                                      stochastic = T)
  
  predict_epi_Inc_clim_uci <- predicttsir(times=dat.pred$time,
                                          births = dat.pred$births,
                                          beta = beta.epi$betahigh,
                                          alpha = fittedpars$alpha,
                                          S0 = SepiIncClim,
                                          I0 = Ifinal,
                                          nsim=100,
                                          stochastic = T)
  
  
  
  
  #now look at data and model predictions for the epidemic year only across all the possibilities
  
  
  #rho will depend on gaussian vs. linear
  
  #1. here just TSIR
  
  IPredEpi1 = cbind.data.frame(time=predict_epi_noInc$I$time, mean_Inc=predict_epi_noInc$I$mean, mean_Inc_lci = predict_epi_noInc_lci$I$mean, mean_Inc_uci = predict_epi_noInc_uci$I$mean)
  names(IPredEpi1) <- c("time", "model_predicted_absolute_cases", "model_predicted_absolute_cases_lci", "model_predicted_absolute_cases_uci")
  # For the epidemic year prediction, use the mean reporting rate across the TSIR fit. This will be 
  # the same value always in the case of S reconstruction by lm but will be a true mean in the case of 
  # years with S reconstruction by gaussian regression
  IPredEpi1$model_predicted_reported_cases <- IPredEpi1$model_predicted_absolute_cases/mean(fittedpars$rho)
  IPredEpi1$model_predicted_reported_cases_lci <- IPredEpi1$model_predicted_absolute_cases_lci/mean(fittedpars$rho)
  IPredEpi1$model_predicted_reported_cases_uci <- IPredEpi1$model_predicted_absolute_cases_uci/mean(fittedpars$rho)
  IPredEpi1$sim_type <- "TSIR-prediction"
   
  
  #2.  here with increased S and normal beta
  IPredEpi2 = cbind.data.frame(time=predict_epi_Inc$I$time, mean_Inc=predict_epi_Inc$I$mean, mean_Inc_lci=predict_epi_Inc_lci$I$mean, mean_Inc_uci=predict_epi_Inc_uci$I$mean)
  names(IPredEpi2) <- c("time", "model_predicted_absolute_cases", "model_predicted_absolute_cases_lci", "model_predicted_absolute_cases_uci")
  IPredEpi2$model_predicted_reported_cases <- IPredEpi2$model_predicted_absolute_cases/mean(fittedpars$rho)
  IPredEpi2$model_predicted_reported_cases_lci <- IPredEpi2$model_predicted_absolute_cases_lci/mean(fittedpars$rho)
  IPredEpi2$model_predicted_reported_cases_uci <- IPredEpi2$model_predicted_absolute_cases_uci/mean(fittedpars$rho)
  IPredEpi2$sim_type <- "increased-S-standard-beta"
  
  
  #3. here with no increase in S but a climate-driven beta
  IPredEpi3 = cbind.data.frame(time=predict_epi_noInc_clim$I$time, mean_Inc=predict_epi_noInc_clim$I$mean, mean_Inc_lci=predict_epi_noInc_clim_lci$I$mean, mean_Inc_uci=predict_epi_noInc_clim_uci$I$mean)
  names(IPredEpi3) <- c("time", "model_predicted_absolute_cases", "model_predicted_absolute_cases_lci", "model_predicted_absolute_cases_uci")
  IPredEpi3$model_predicted_reported_cases <- IPredEpi3$model_predicted_absolute_cases/mean(fittedpars$rho)
  IPredEpi3$model_predicted_reported_cases_lci <- IPredEpi3$model_predicted_absolute_cases_lci/mean(fittedpars$rho)
  IPredEpi3$model_predicted_reported_cases_uci <- IPredEpi3$model_predicted_absolute_cases_uci/mean(fittedpars$rho)
  IPredEpi3$sim_type <- "no-increase-S-climate-beta"
  
  
  #4. here with both increased S and climate driven beta
  
  IPredEpi4 = cbind.data.frame(time=predict_epi_Inc_clim$I$time, mean_Inc=predict_epi_Inc_clim$I$mean, mean_Inc_lci=predict_epi_Inc_clim_lci$I$mean, mean_Inc_uci=predict_epi_Inc_clim_uci$I$mean)
  names(IPredEpi4) <- c("time", "model_predicted_absolute_cases", "model_predicted_absolute_cases_lci", "model_predicted_absolute_cases_uci")
  IPredEpi4$model_predicted_reported_cases <- IPredEpi4$model_predicted_absolute_cases/mean(fittedpars$rho)
  IPredEpi4$model_predicted_reported_cases_lci <- IPredEpi4$model_predicted_absolute_cases_lci/mean(fittedpars$rho)
  IPredEpi4$model_predicted_reported_cases_uci <- IPredEpi4$model_predicted_absolute_cases_uci/mean(fittedpars$rho)
  IPredEpi4$sim_type <- "increased-S-climate-beta"
  
  
  # fitted TSIR
  IFit = cbind.data.frame(time=predict_ts$I$time, mean_Inc=predict_ts$I$mean, mean_Inc_lci=predict_ts_lci$I$mean, mean_Inc_uci=predict_ts_uci$I$mean)
  names(IFit) <- c("time", "model_predicted_absolute_cases", "model_predicted_absolute_cases_lci", "model_predicted_absolute_cases_uci")
  IFit$model_predicted_reported_cases <- IFit$model_predicted_absolute_cases/mean(fittedpars$rho)
  IFit$model_predicted_reported_cases_lci <- IFit$model_predicted_absolute_cases_lci/mean(fittedpars$rho)
  IFit$model_predicted_reported_cases_uci <- IFit$model_predicted_absolute_cases_uci/mean(fittedpars$rho)
  
  IFit$sim_type <- "TSIR-fit"
  
  #bind all the predictions
  
  IPred <- rbind(IPredEpi1, IPredEpi2, IPredEpi3, IPredEpi4, IFit)
  
  
  

  #ggplot(data = IPred) + geom_line(aes(x=time, y= model_predicted_reported_cases, color=sim_type))
  
  
  #and merge with data on time
  
  df.merge = subset(df, time<=(epiyr1+1))
  all.dat.merge <- merge(df.merge, IPred, by="time", all.x = T, sort=F)
  
   #  ggplot(data=all.dat.merge) + 
   #    geom_line(aes(x=time, y=cases),linetype=2) +
   #    geom_line(aes(x=time, y=model_predicted_reported_cases, color=sim_type)) 
   # # 
   # 
   #and return both sets of data
   
   return(list(incS.df, all.dat.merge))
}


###################################################################################

#load province data in tsir form - with associated beta info and climate (including climate projections)
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_dat_province.csv"), header = T, stringsAsFactors = F)
head(tsir.dat)
 
length(unique(tsir.dat$provname)) #25 unique provinces
 
 #remove province with uneven time series
tsir.dat <- subset(tsir.dat, provname!="Tboung Khmum" & provname!="Mondul Kiri" & provname!= "Ratanak Kiri")
tsir.dat <- arrange(tsir.dat, provname, time) # just verify it is in order
length(unique(tsir.dat$provname)) #22 used for climate projections

#load beta info for fitted years with optimal susceptible reconstruction
beta.df <- read.csv(file = paste0(homewd, "/data/beta_TSIR_fit_province.csv"), header = T, stringsAsFactors = F)
sus.merge <- ddply(beta.df, .(provname, epiyr), summarise, rsquared=unique(rsquared), sus_reconstruction=unique(sus_reconstruction))
head(sus.merge)
sus.merge <- arrange(sus.merge, provname, epiyr)
unique(sus.merge$provname)

#and load the epiyear beta predictions using climate (these are GAM predictions)
clim.dat <- read.csv(file= paste0(homewd, "/data/tsir_dat_beta_climate_province.csv"), header = T, stringsAsFactors = F)
unique(clim.dat$provname)
clim.dat$provname[clim.dat$provname=="Oddar Meanchey"] <- "Otdar Meanchey"
clim.dat <- arrange(clim.dat, provname, time)
head(clim.dat)

#split by province and pre-epidemic period
tsir.split.2007 <- dlply(tsir.dat,.(provname))

tsir.dat.2012 <- subset(tsir.dat, time > 2007.9999)
tsir.split.2012 <- dlply(tsir.dat.2012,.(provname))

tsir.dat.2019 <- subset(tsir.dat, time > 2012.9999)
tsir.split.2019 <- dlply(tsir.dat.2019,.(provname))


# and apply it across the board. save the data and plot it after

clim.project.tsir.2007 <- lapply(tsir.split.2007, wrap.pipeline.TSIR, epiyr1=2007, sbar1 = NULL, sus.dat = sus.merge, epi.beta.df = clim.dat)
clim.project.tsir.2012 <- lapply(tsir.split.2012, wrap.pipeline.TSIR, epiyr1=2012, sbar1 = NULL, sus.dat = sus.merge, epi.beta.df = clim.dat)
clim.project.tsir.2019 <- lapply(tsir.split.2019, wrap.pipeline.TSIR, epiyr1=2019, sbar1 = NULL, sus.dat = sus.merge, epi.beta.df = clim.dat)


# save the first element from each which is the ratio of susceptible increases 
# needed to recapture the epidemic data, both with and without the climate-driven beta
sus.df.2007 <-data.table::rbindlist(sapply(clim.project.tsir.2007,'[', 1))
sus.df.2012 <-data.table::rbindlist(sapply(clim.project.tsir.2012,'[', 1))
sus.df.2019 <-data.table::rbindlist(sapply(clim.project.tsir.2019,'[', 1))

sus.clim.df <- rbind.data.frame(sus.df.2007, sus.df.2012, sus.df.2019) 
#in 2007, - values for 5 provinces, + values for 11 provinces, 0 values for 6 provinces
#in 2012, - values for 14 provinces, + values for 3 provinces, 0 values for 5 provinces
#in 2019, - values for 12 provinces, + values for 7 provinces, 0 values for 3 provinces
write.csv(sus.clim.df, paste0(homewd,"/data/sus_increase_by_clim_by_prov.csv"), row.names = F)

#and compare by t test
t.test(sus.clim.df$frac_incS_noclim[sus.clim.df$epidemic_year==2007], sus.clim.df$frac_incS_clim[sus.clim.df$epidemic_year==2007], paired = T, alternative = "greater")
t.test(sus.clim.df$frac_incS_noclim[sus.clim.df$epidemic_year==2012], sus.clim.df$frac_incS_clim[sus.clim.df$epidemic_year==2012], paired = T, alternative = "greater")
t.test(sus.clim.df$frac_incS_noclim[sus.clim.df$epidemic_year==2019], sus.clim.df$frac_incS_clim[sus.clim.df$epidemic_year==2019], paired = T, alternative = "greater")

#and write to table s5
tableS5 <- dplyr::select(sus.clim.df, epidemic_year, frac_incS_noclim, frac_incS_clim, provname)
write.csv(tableS5, file=paste0(homewd,"/data/tableS5.csv"), row.names=F)

# and save the second element from each, which are the projection outputs for plotting
project.df.2007 <-data.table::rbindlist(sapply(clim.project.tsir.2007,'[', 2))
project.df.2012 <-data.table::rbindlist(sapply(clim.project.tsir.2012,'[', 2))
project.df.2019 <-data.table::rbindlist(sapply(clim.project.tsir.2019,'[', 2))

project.df.all <- rbind.data.frame(project.df.2007, project.df.2012, project.df.2019)

#save
write.csv(project.df.all, paste0(homewd,"/data/climate_projections_TSIR_by_prov.csv"), row.names = F)



project.df.all$sim_type <- factor(project.df.all$sim_type, levels=c("TSIR-fit", "TSIR-prediction", "increased-S-standard-beta", "no-increase-S-climate-beta", "increased-S-climate-beta"))

colz = c('TSIR-fit'="#00B0F6",'TSIR-prediction' = "#E76BF3", 'increased-S-standard-beta' = "#A3A500", 'no-increase-S-climate-beta'="#00BF7D",'increased-S-climate-beta' = "#F8766D" )

ggplot(subset(project.df.all, provname=="Battambang")) + facet_wrap(~provname, scales = "free_y") +
  geom_line(aes(x=time, y=cases), linetype=2) + scale_color_manual(values=colz) + scale_fill_manual(values=colz) + theme_bw() +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
  geom_ribbon(data = subset(project.df.all, provname=="Battambang" &time<2008),
              aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(project.df.all, provname=="Battambang" &time<2008),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) +
  geom_ribbon(data = subset(project.df.all, provname=="Battambang" & time>=2008 & time<2012  & sim_type=="TSIR-fit"),
              aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(project.df.all, provname=="Battambang" & time>=2008 & time<2012  & sim_type=="TSIR-fit"),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) +
  geom_ribbon(data = subset(project.df.all, provname=="Battambang" & time>=2012 & time<2013 & sim_type!="TSIR-fit"),
              aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(project.df.all, provname=="Battambang" & time>=2012 & time<2013 & sim_type!="TSIR-fit"),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) +
  geom_ribbon(data = subset(project.df.all, provname=="Battambang" & time>=2013 & time<2019  & sim_type=="TSIR-fit"),
              aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(project.df.all, provname=="Battambang" & time>=2013 & time<2019  & sim_type=="TSIR-fit"),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) +
  geom_ribbon(data = subset(project.df.all, provname=="Battambang" & time>=2019 & sim_type!="TSIR-fit"),
              aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(project.df.all, provname=="Battambang" & time>=2019 & sim_type!="TSIR-fit"),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) 
  
  
# geom_line(data = subset(project.df.all, provname=="Battambang" &time<2008),
#             aes(x=time, y=model_predicted_reported_cases, color=sim_type))
  



p1 <- ggplot(project.df.all) + facet_wrap(~provname, scales = "free_y") + ylab("reported cases") +
  geom_line(aes(x=time, y=cases), linetype=2) + scale_color_manual(values=colz, labels=scales::parse_format()) +
  scale_fill_manual(values=colz, labels=scales::parse_format()) + theme_bw() +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"), 
        legend.title = element_blank(), legend.position = c(.91,.08), axis.title.x = element_blank()) +
  geom_ribbon(data = subset(project.df.all, time<2008),
            aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(project.df.all, time<2008),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) +
  geom_ribbon(data = subset(project.df.all,time>=2008 & time<2012  & sim_type=="TSIR-fit"),
              aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(project.df.all,time>=2008 & time<2012  & sim_type=="TSIR-fit"),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) +
  geom_ribbon(data = subset(project.df.all, time>=2012 & time<2013 & sim_type!="TSIR-fit"),
              aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(project.df.all, time>=2012 & time<2013 & sim_type!="TSIR-fit"),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) +
  geom_ribbon(data = subset(project.df.all, time>=2013 & time<2019  & sim_type=="TSIR-fit"),
              aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(project.df.all, time>=2013 & time<2019  & sim_type=="TSIR-fit"),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) +
  geom_ribbon(data = subset(project.df.all, time>=2019 & sim_type!="TSIR-fit"),
              aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(project.df.all, time>=2019 & sim_type!="TSIR-fit"),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) 
  



ggsave(file = paste0(homewd, "/final-figures/TSIR-Prov-Climate-Projections.png"),
       plot = p1,
       units="mm",  
       width=100, 
       height=70, 
       scale=3, 
       dpi=300)

