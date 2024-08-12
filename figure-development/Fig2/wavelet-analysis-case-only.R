rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(WaveletComp)
library(mgcv)
library(reshape2)

#making wavelet reconstructions of the case data only
#first, load the time series data
homewd= "/Users/carabrook/Developer/cambodia-dengue-national-public"

#load tsir data
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_dat_province.csv"), header = T, stringsAsFactors = F)
head(tsir.dat)
tail(tsir.dat) #goes from beginning of 2002 to end of 2020

#incidence per 100,000 pop
tsir.dat$cases_per_100000 <- (tsir.dat$cases/tsir.dat$pop)*100000

tsir.dat <- arrange(tsir.dat, provname, time)



#first, for each time series per province, run a function to collect:
#(a) the reconstructed period for each timestep, both annual
#(b) and multi-annual
#(c) the average wavelet power per timestep, for annual
#(d) and multi-annual

prov.split <- dlply(tsir.dat, .(provname))

prov.rank.multi <- function(df2, df1, pval.cut){
  
  df1a <-  df1
  df2a <- df2
  df1 <- dplyr::select(df1, time, cases_per_1000)
  df2 <- dplyr::select(df2, time, cases_per_1000)
  
  names(df1)[names(df1)=="cases_per_1000"] <- "cases_this_prov"
  names(df2)[names(df2)=="cases_per_1000"] <- "cases_other_prov"
  
  df.merge <- merge(df1, df2, by="time", all.x = T)
  head(df.merge)
  
  df.merge <- df.merge[complete.cases(df.merge),]
  df.merge <- arrange(df.merge, time)
  
  corr.prov <- analyze.coherency(df.merge, my.pair = c("cases_this_prov","cases_other_prov"),
                                 loess.span = 0,
                                 dt = 1/26, dj = 1/100,
                                 window.type.t = 1, window.type.s = 1,
                                 #window.size.t = (1), #examine coherence year-by-year
                                 window.size.t = (26*5), #examine coherence year-by-year - here 5 years
                                 #window.size.t = (26*5), #examine coherence year-by-year
                                 window.size.s = (1/4), #periods on the order of 
                                 lowerPeriod = 2, #shortest possible period in years (multi)
                                 upperPeriod = 7, #largest possible period (in weeks; here, 7 years)
                                 make.pval = TRUE, n.sim = 100)
  
  
  # # # 
  #      wc.image(corr.prov, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #               periodlab= "period (in years)", spec.time.axis = list(at = seq(1,26*18, 26), labels = 2002:2019), which.arrow.sig = "wt", color.key = "interval")
  # # # # # # 
  # # # # 
  #     wc.image(corr.prov,which.image = "wc",  
  #              color.key = "interval", n.levels = 250,
  #              siglvl.contour = 0.1, siglvl.arrow = 0.01,
  #              legend.params = list(lab = "wavelet coherence levels"),
  #              spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #              periodlab= "period (in years)", spec.time.axis = list(at = seq(1,26*18, 26), labels = 2002:2019))
  # # # # # 
  # # # 
  # 
  #is there a statistically significant coherence through time?
  coherence.mat <- corr.prov$Coherence
  pval.mat <- corr.prov$Coherence.pval
  pval.mat <- pval.mat< pval.cut
  
  out.mat <- coherence.mat*pval.mat
  out.mat[out.mat==0] <- NA
  
  provname <- sub(" ", "_", unique(df2a$provname), fixed = T)
  
  #df.merge$coherence_annual <- colSums(pval.mat)
  df.merge$coherence_multi <- colMeans(out.mat, na.rm = T)
  
  
  df2add <- dplyr::select(df.merge, time, coherence_multi)
  
  df1a <- merge(df1a, df2add, by="time", all.x = T)
  
  df1a$coherence_prov <- provname
  
  df.out <- dplyr::select(df1a, time, coherence_multi, coherence_prov)
  df.out <- arrange(df.out, time)
  #ggplot(df.out) + geom_line(aes(x=time, y=coherence_multi))
  
  
  #add in cross wavelet power
  power.mat <- corr.prov$Power.xy
  pval.mat <- corr.prov$Power.xy.pval
  pval.mat <- pval.mat< pval.cut
  
  power.mat <- power.mat*pval.mat
  power.mat[power.mat==0] <- NA
  
  df.out$cross_wave_power_multi <- colMeans(power.mat, na.rm=T)
  
  
  #ggplot(df.out) + geom_line(aes(x=time, y=cross_wave_power_multi))
  return(df.out)
  
  
}
prov.rank.annual <- function(df2, df1, pval.cut){
  
  df1a <-  df1
  df2a <- df2
  df1 <- dplyr::select(df1, time, cases_per_1000)
  df2 <- dplyr::select(df2, time, cases_per_1000)
  
  names(df1)[names(df1)=="cases_per_1000"] <- "cases_this_prov"
  names(df2)[names(df2)=="cases_per_1000"] <- "cases_other_prov"
  
  df.merge <- merge(df1, df2, by="time", all.x = T)
  head(df.merge)
  
  df.merge <- df.merge[complete.cases(df.merge),]
  df.merge <- arrange(df.merge, time)
  
  corr.prov <- analyze.coherency(df.merge, my.pair = c("cases_this_prov","cases_other_prov"),
                                 loess.span = 0,
                                 dt = 1/26, dj = 1/100,
                                 window.type.t = 1, window.type.s = 1,
                                 window.size.t = 26, #examine coherence year-by-year
                                 window.size.s = (1/4), #periods on the order of 
                                 lowerPeriod = 13/26, #shortest possible period in years (annual)
                                 upperPeriod = 2, #largest possible period (2 years)
                                 make.pval = TRUE, n.sim = 100)
  
  
  # 
  #     wc.image(corr.prov, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #              periodlab= "period (in years)", spec.time.axis = list(at = seq(1,26*18, 26), labels = 2002:2019), which.arrow.sig = "wt", color.key = "interval")
  # # # # # 
  # # # # 
  #     wc.image(corr.prov,which.image = "wc",  
  #              color.key = "interval", n.levels = 250,
  #              siglvl.contour = 0.1, siglvl.arrow = 0.05,
  #              legend.params = list(lab = "wavelet coherence levels"),
  #              spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #              periodlab= "period (in years)", spec.time.axis = list(at = seq(1,26*18, 26), labels = 2002:2019))
  # # # # # 
  # # # 
  # 
  #is there a statistically significant coherence through time?
  coherence.mat <- corr.prov$Coherence
  pval.mat <- corr.prov$Coherence.pval
  pval.mat <- pval.mat< pval.cut
  
  out.mat <- coherence.mat*pval.mat
  out.mat[out.mat==0] <- NA
  
  provname <- sub(" ", "_", unique(df2a$provname), fixed = T)
  
  #df.merge$coherence_annual <- colSums(pval.mat)
  df.merge$coherence_annual <- colMeans(out.mat, na.rm = T)
  
  
  df2add <- dplyr::select(df.merge, time, coherence_annual)
  
  df1a <- merge(df1a, df2add, by="time", all.x = T)
  
  df1a$coherence_prov <- provname
  
  df.out <- dplyr::select(df1a, time, coherence_annual, coherence_prov)
  df.out <- arrange(df.out, time)
  #ggplot(df.out) + geom_line(aes(x=time, y=coherence_annual))
  
  
  #add in cross wavelet power
  power.mat <- corr.prov$Power.xy
  pval.mat <- corr.prov$Power.xy.pval
  pval.mat <- pval.mat< pval.cut
  
  power.mat <- power.mat*pval.mat
  power.mat[power.mat==0] <- NA
  
  df.out$cross_wave_power_annual <- colMeans(power.mat, na.rm=T)
  
  
  #ggplot(df.out) + geom_line(aes(x=time, y=cross_wave_power_annual))
  return(df.out)
  
  
}
get.mean.period <- function(col.df){
  
  #bind
  dat.period <- cbind.data.frame(period=col.df$period_length, power=col.df$value)
  
  mean_period = weighted.mean(x=col.df$period_length, w=col.df$value)
  
  return(mean_period)
}
get.max.period <- function(col.df){
  
  #bind
  #dat.period <- cbind.data.frame(period=col.df$period_length, power=col.df$value)
  
  max_period = col.df$period_length[col.df$value==max(col.df$value)]
  
  if(length(max_period)>1){
    max_period = min(max_period) #if more than one, take the shortest
  }
  
  return(max_period)
}
get.wavelet.case.dat <- function(dat, dat.all){
  
  dat.all = subset(dat.all, provname !=unique(dat$provname))
  #first, get reconstructed period for annual and multiannual using the biweekly case data
  
  #(a) annual
  anal.dat.annual <-  analyze.wavelet(dat,
                                     my.series = ncol(dat), #cases per 100000
                                     loess.span = 1, #detrending (e.g. smoothing) over annual timestep
                                     dt = 1/26,#this allows for annual timestep if biweekly timeseries
                                     dj = 1/100, #default =1/20. image gets clearer as number on the bottom gets bigger
                                     lowerPeriod = 1/26,#shortest possible period (one biweek)
                                     upperPeriod = 2, #largest possible period ( here, 2 years)
                                     make.pval = TRUE, n.sim = 100)
  
   
  #   wt.image(anal.dat.annual, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #            periodlab= "period (in years)", spec.time.axis = list(at = seq(1,26*18, 26), labels = 2002:2019))
  # # # #
  
  
  anal.dat = reconstruct(anal.dat.annual, only.sig = T)
  
  
  dat$reconstructed_annual_period <- anal.dat$series$cases_per_100000.r
  
  # (b) multi
  anal.dat.multi <-  analyze.wavelet(dat,
                                     my.series = "cases_per_100000", #cases per 1000
                                     loess.span = 1,
                                     dt = 1/26,#this allows for annual timestep
                                     dj = 1/100, #default =1/20. image gets clearer as number on the bottom gets bigger
                                     lowerPeriod = 2,#shortest possible period (2 years)
                                     upperPeriod = 7, #largest possible period (in weeks; here, 20 years)
                                     make.pval = TRUE, n.sim = 100)
  
  #  wt.image(anal.dat.multi, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #           periodlab= "period (in years)", spec.time.axis = list(at = seq(1,26*18, 26), labels = 2002:2019))
  # # 
  # 
  multi.dat = reconstruct(anal.dat.multi, only.sig = T)
  
  
  #reconstruct(anal.dat.multi, only.ridge = T)
  dat$reconstructed_multi_period <- multi.dat$series$cases_per_100000.r
  
  #and, as bonus, get the average length of the statistically significant multi-annual dengue cycles
  #col.split = as.list(as.data.frame(anal.dat.multi$Power))
  
  # 
  # power.mat <- anal.dat.multi$Power
  # pval.mat <- anal.dat.multi$Power.pval
  # pval.mat <- pval.mat<1
  # power.df <- power.mat*pval.mat
  # power.df[power.df==0] <- NA
  
  power.df = as.data.frame(anal.dat.multi$Power)
  names(power.df) <- 1:length(dat$reconstructed_multi_period)
  power.df$period <- as.numeric(rownames(power.df))
  
  
  power.df.melt <- melt(power.df, id.vars = "period")
  names(power.df.melt)[names(power.df.melt)=="variable"] <- "time"
  power.df.melt$time <- as.numeric(as.character(power.df.melt$time))
  
  period.df = cbind.data.frame(period=1:length(unique(power.df.melt$period)), period_length=anal.dat.multi$Period)
  
  power.df.melt <- merge(power.df.melt, period.df, by="period", all.x=T)
  #power.df.melt <-   power.df.melt[complete.cases(power.df.melt),]
  power.df.melt <- arrange(power.df.melt, time)
  
  col.split <- dlply(power.df.melt, .(time))
  
  #ggplot(power.df.melt) + geom_tile(aes(x=time, y=period, fill=value))
  
  
  mean.period.list <- lapply(col.split, get.mean.period)
  # period.df <- data.table::rbindlist(mean.period.list)
  # names(period.df) <- c("time_index","multi_period_mean")
  # 
  # dat$time_index <- 1:nrow(dat)
  # dat <- merge(dat, period.df, by="time_index", all.x = T, sort = F)
  #max.period.list <- lapply(col.split, get.max.period)
  
  
  
  dat$multi_period_mean <-  c(unlist(mean.period.list))
  
  #ggplot(dat) + geom_line(aes(x=time, y=multi_period_mean))
  
  #(c) then, average statistically significant wavelet power for annual
  power.mat <- anal.dat.annual$Power
  pval.mat <- anal.dat.annual$Power.pval
  pval.mat <- pval.mat<0.05
  
  avg.mat <- power.mat*pval.mat
  avg.mat[avg.mat==0] <- NA
  
  
  dat$avg_wave_power_annual <- (colMeans(avg.mat, na.rm = T))
  #ggplot(dat) + geom_line(aes(x=time, y=avg_wave_power_annual)) + geom_vline(aes(xintercept=2007), color="red") + geom_vline(aes(xintercept=2012, color="red")) + geom_vline(aes(xintercept=2019, color="red"))
  
  #(d) average stat sig wavelet power for multiannual
  power.mat <- anal.dat.multi$Power
  pval.mat <- anal.dat.multi$Power.pval
  pval.mat <- pval.mat<0.05
  
  avg.mat <- power.mat*pval.mat
  avg.mat[avg.mat==0] <- NA
  
  dat$avg_wave_power_multi <- (colMeans(avg.mat, na.rm = T))
  #ggplot(dat) + geom_line(aes(x=time, y=avg_wave_power_multi)) + geom_vline(aes(xintercept=2007), color="red") + geom_vline(aes(xintercept=2012, color="red")) + geom_vline(aes(xintercept=2019, color="red"))
  
  
  #head(dat)
  dat <- arrange(dat, time)
  
  return(dat)
  
}
get.wavelet.oni <- function(dat, dat.all, oni.df, pval.cut.coherence, pval.cut.power){
  
  dat.all = subset(dat.all, provname !=unique(dat$provname))
  #first, get reconstructed period for annual and multiannual using the biweekly case data
  
  #(e) average wavelet coherency with ONI (multi only and must be at the monthly scale)
  biwk.month.dat <- cbind.data.frame(biweek=1:26, month=c(1,rep(1:12, each=2),12)) 
  dat <- merge(dat, biwk.month.dat, by="biweek")
  dat <- arrange(dat, time)
  case.sum <- ddply(dat, .(provname, year, month), summarise, time = min(time), cases_per_100000 = sum(cases_per_1000))
  
  
  #add oni
  
  oni.merge <- dplyr::select(oni.dat, year, oni_anomaly, month)
  case.sum <- merge(case.sum, oni.merge, by=c("year", "month"))
  
  #and look for the cross correlation
  case.sum <- arrange(case.sum, provname, time)
  
  corr.oni <- analyze.coherency(case.sum, my.pair = c("oni_anomaly","cases_per_1000"),
                                loess.span = 0,
                                dt = 1/12, 
                                dj = 1/100,
                                window.type.t = 1, window.type.s = 1,
                                #window.size.t = 1, #examine coherence year-by-year
                                window.size.t = (12*5), #examine coherence on a 5-year window for ENSO
                                window.size.s = (1/4), #periods on the order of 
                                lowerPeriod = 2, #shortest possible period in years
                                upperPeriod = 7, #largest possible period (in weeks; here, 7 years)
                                make.pval = TRUE, n.sim = 100)
  # # # 
  #      wc.image(corr.oni, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),siglvl.contour = 0.01,
  #               periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019), which.arrow.sig = "wt", color.key = "interval")
  #               
  # # # # # # 
  # # # # # 
  #      wc.image(corr.oni,which.image = "wc",  
  #              color.key = "interval", n.levels = 250,
  #               siglvl.contour = 0.1, siglvl.arrow = 0.05, #contours around associations significant by pval =0.1. arrows at pval=0.05
  #               legend.params = list(lab = "wavelet coherence levels"),
  #               spec.period.axis = list(at = c(2:10), labels = c(2:10)),
  #               periodlab= "period (in years)", spec.time.axis = list(at = seq(1,12*18, 12), labels = 2002:2019))
  # # # # # # 
  # # # 
  #monthly average significant coherence
  #remove those not significant
  coherence.mat <- corr.oni$Coherence
  pval.mat <- corr.oni$Coherence.pval
  pval.mat <- pval.mat< pval.cut.coherence
  
  out.mat <- coherence.mat*pval.mat
  out.mat[out.mat==0] <- NA
  #plot(colMeans(out.mat))
  #wc.avg(corr.oni)
  
  
  
  case.sum$avg_wave_coherence_oni <- (colMeans(out.mat,na.rm = T))
  #ggplot(case.sum) + geom_line(aes(x=time, y=avg_wave_coherence_oni)) + geom_vline(aes(xintercept=2007), color="red") + geom_vline(aes(xintercept=2012, color="red")) + geom_vline(aes(xintercept=2019, color="red"))
  
  #with(dat, plot(month_date, avg_cross_power, type="b"))
  
  power.table <- corr.oni$Power.xy
  power.table.pval <- corr.oni$Power.xy.pval
  power.table.pval <- power.table.pval< pval.cut.power
  
  power.mat <- power.table*power.table.pval
  power.mat[power.mat==0] <- NA
  
  
  case.sum$avg_cross_power_oni <-  (colMeans(power.mat, na.rm = T))
  #ggplot(case.sum) + geom_line(aes(x=time, y=avg_cross_power_oni)) + geom_vline(aes(xintercept=2007), color="red") + geom_vline(aes(xintercept=2012, color="red")) + geom_vline(aes(xintercept=2019, color="red"))
  
  
  return(case.sum)
  
}


prov.split.out <- lapply(prov.split, get.wavelet.case.dat, dat.all=tsir.dat)
prov.split.df <- data.table::rbindlist(prov.split.out)
head(prov.split.df)
names(prov.split.df)


write.csv(prov.split.df, file = paste0(homewd, "/data/synchrony_case_data_oct15.csv"), row.names = F)
