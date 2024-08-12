# testing our framework using the simplest SIR model

rm(list=ls())

#homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
#setwd(paste0(homewd, "/figure-development/Fig5/fit-sim"))
#set wdset wd
#.libPaths("~/R/x86_64-pc-linux-gnu-library/3.4")

library(plyr)
library(dplyr)
library(matrixcalc)
library(Matrix)
library(ggplot2)
library(mvtnorm)
library(reshape2)



sum.yr.all.sim <- function(df){
  
  df.sum <- ddply(df, .(year, age), summarise, Nage = sum(count))
  df.sum$Nage <- round(df.sum$Nage, 0)
  
  df.out = cbind.data.frame(age=1:max(df$age))
  
  df.out <- merge(df.out, df.sum, by="age", all.x = T, sort = F)
  df.out$Nage[is.na(df.out$Nage)] <- 0
  df.out$year[is.na(df.out$year)] <- unique(df$year)
  #df.out <- rbind(c(0,0), df.out)
  #bind
  df.add <- cbind.data.frame(age=0, year = unique(df$year), Nage=0)
  df.out <- rbind(df.add, df.out)
  df.out <- arrange(df.out, age)
  df.out$cum_cases <- cumsum(df.out$Nage)
  df.out$n <- sum(df.sum$Nage)
  df.out$cum_prop_cases <- df.out$cum_cases/df.out$n
  
  return(df.out)
  
}
model.age.incidence.series <- function(par.dat, age_vect, year.start){
  
  # you will have one lambda for each year in the time series
  # and one N-sero for each year as well
  # and X ages with varying probabilities of infection within those years
  
  lambda.start = par.dat$lambda[1]
  
  
  
  lts = nrow(par.dat)
  
  #proportion exposed in any way - needs to long enough for each age within each year
  pexposed= rep(NA, (length(age_vect)-1))
  pexposed[1] <- 0
  pexposed = rep(list(pexposed),lts) 
  
  #naive
  pnaive= rep(NA, (length(age_vect)-1))
  pnaive[1] <- 1
  pnaive = rep(list(pnaive),lts) 
  
  
  #exposed to a single serotype only (this is one serotype's worth of primary infections)
  pprim= rep(NA, (length(age_vect)-1))
  pprim[1] <- 0
  pprim = rep(list(pprim),lts)
  
  #multitypic infections
  pmulti= rep(NA, (length(age_vect)-1))
  pmulti[1] <- 0
  pmulti = rep(list(pmulti),lts)
  
  age_tracker= rep(NA, (length(age_vect)-1))
  age_tracker[1] <- 0
  age_tracker = rep(list(age_tracker),lts)
  
  year_tracker= rep(NA, (length(age_vect)-1))
  year_tracker[1] <- year.start
  year_tracker = rep(list(year_tracker),lts)
  
  
  
  
  for (i in 1:lts){
    
    if(i < max(age_vect)){
      age_vect_tmp = age_vect[1:which(age_vect==i)]
    }else{
      age_vect_tmp = age_vect
    }
    
    for(a in 2:(length(age_vect_tmp))){ 
      # for-loops over all possible ages in our data range across all the years. 
      # all ages get rounded up, so babies <1 year are counted as 1 year old,
      # and experienced a hazard of infection in the year prior
      
      
      
      
      age = age_vect_tmp[a]
      
      year.now = par.dat$year[i]
      year.par = subset(par.dat, year<=year.now)
      N_sero = year.par$N_sero
      lambda =  year.par$lambda #one per year
      
      # make a vector of durations for the lambdas across the time series
      # in this model, it is easy because the duration is only 1 year...
      # you may eventually also want partial year durations in a different model
      #if (age>1){
      dur = rep(1, length(lambda))  
      #}else if(age<=1){
      # dur = rep(0, length(lambda))  
      #}
      
      #if you were born after the time series began
      #some of those ones may need to be replaced with 0s
      if(i>1){
        #diff.count <- i-age
        tot.class = ceiling(age)
        #diff.count <- sum(dur)-age_trunc
        
        if(tot.class<length(dur)){
          dur[1:(length(dur)-tot.class)]<- 0
        }
      }
      
      # #if 0, this kid spent 1 year in each age class and the dur vector does not need to be changed
      if(sum(dur) != age){
        warning(paste0("mismatch in age and integration at age=", age, " and i=", i))
      }
      
      
      
      # now, we integrate the hazard of becoming infected over all the years in each
      
      # First, integrand A deals with avoiding infection with all strains until now
      # Then, integrand B deals with getting exposed to that select strain.
      
      # first, integrand A  (avoiding infection with all other strains) 
      inte_A = sum(dur*lambda*N_sero)
      
      # then integrand B (getting exposed to a single strain)
      inte_B = sum(dur*lambda)
      
      # now fill in the probabilities based on these : 
      #this is the probability of being naive to all serotypes at this point in time
      pnaive[[i]][[a]] <- exp(-inte_A)
      
      #this is the probability of being exposed to any serotype at the current time
      pexposed[[i]][[a]] <- (1-exp(-inte_A))
      
      # this is the probability of being exposed to just one serotype (a primary infection)
      # here, this is expressed for just a single target strain 
      # (would need to multiply by the number of circulating strains if you wanted 
      # the probability of a primary infection with ANY strain)
      pprim[[i]][[a]] <- exp(-inte_A)*(exp(inte_B)-1)
      
      
      # if not primarily infected or naive, this should be a multitypic infection
      # but remember it could be a primary infection with any of the four serotypes 
      # so you need to account for that here by again multiplying by the 
      # number of circulating serotypes internally
      pmulti[[i]][[a]] = 1 - pnaive[[i]][[a]] - (N_sero[i]*(pprim[[i]][[a]]))
      age_tracker[[i]][[a]] = age
      year_tracker[[i]][[a]] =year.now
      
    }
  }
  
  #some of the age distributions will be shorter
  
  #and get the estimates of each
  #p.out = cbind.data.frame(year = rep(seq((year.start), (year.start+lts-1), 1), each=length(age_vect)),
  p.out <-  cbind.data.frame(year=c(unlist(year_tracker)),
                             age=c(c(unlist(age_tracker))), 
                             exposed=c(c(unlist(pexposed))),
                             naive=c(c(unlist(pnaive))),
                             one_prim=c(unlist(pprim)), 
                             multi = c(unlist(pmulti)))
  
  merge.par <- dplyr::select(par.dat, year, N_sero)
  p.out <- merge(p.out, merge.par, by="year", all.x = T, sort = F)
  
  
  p.out$all_prim <- p.out$one_prim*p.out$N_sero
  
  
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(cbind(p.out["exposed"], p.out["naive"])) #should be 1
  p.out$sum_naive_prim_multi <- rowSums(cbind(p.out["naive"], p.out["all_prim"], p.out["multi"])) #should be 1
  
  
  p.out <- arrange(p.out, year, age)
  
  # ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  + 
  # geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  # 
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  p.out$cum_prop_cases <- p.out$multi
  p.out <- p.out[complete.cases(p.out),]
  
  
  p.add <- cbind.data.frame(year = unique(p.out$year), age = rep(0, length(unique(p.out$year))), exposed=0, naive=1, one_prim=0, multi=0, N_sero = N_sero[1], all_prim=0, sum_exp_naive=1, sum_naive_prim_multi=1, cum_prop_cases =0)
  
  
  p.sum <- rbind(p.add, p.out)
  p.sum <- arrange(p.sum, year, age)
  
  p.sum <- p.sum[!duplicated(p.sum),]
  
  p1 <- ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + 
    geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  # print(p1)
  # 
  #ggplot(data=subset(p.sum, year==1985)) + geom_point(aes(x=age, y=multi))  + 
  #geom_line(aes(x=age, y=multi)) 
  
  return(p.sum) #returns prevalence by age for each year for fitting to the years for which we have data
}
log.lik.fit.all <- function(par, par.dat, dat, year.start){
  
  par.dat$lambda <- par
  
  
  age_vect=seq(0, max(dat$age), by=1)
  
  out.mod <- model.age.incidence.series(par.dat = par.dat, year.start = year.start, age_vect=age_vect)  
  
  #ggplot(data=out.mod) + geom_point(aes(x=age,y= cum_prop_cases)) + facet_wrap(~year)
  
  # now, select only the years for fitting for which there are data
  # (the model can produce data back further in time based on the first year FOI from the dataset)
  out.mod = subset(out.mod, year >=min(dat$year))
  
  
  # merge model with data
  out.mod <- arrange(out.mod, year, age)
  dat <- arrange(dat, year, age)
  dat.merge <- dplyr::select(dat, age, year, Nage, cum_cases, n)
  out.merge <- merge(out.mod, dat.merge, by= c("year", "age"))
  out.merge$cum_prop_cases_data <- out.merge$cum_cases/out.merge$n
  
  out.merge <- arrange(out.merge, year, age)
  
  # #plot model with data
  # ggplot(data=out.merge) + facet_wrap(~year) +  geom_point(aes(x=age, y=cum_prop_cases_data)) + geom_line(aes(x=age, y=cum_prop_cases_data)) + geom_line(aes(x=age,y= cum_prop_cases), color="tomato") 
  # #   
  
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(out.merge$age)){
    ll=ll+dbinom(out.merge$cum_cases[i],out.merge$n[i],prob=out.merge$cum_prop_cases[i],log=T)
    #print(ll)
  }
  
  # 
  
  par.dat$lambda_high <- 0 
  par.dat$lambda_high[par.dat$lambda>1] <- 1 #penalty if lambda is >1
  tot_over = sum(par.dat$lambda_high)
  
  if(tot_over>0){
    
    ll <- ll-(tot_over*1000000)
    
  }
  
  if(ll==-Inf){
    ll <- -1000000
  }
  
  return(-ll)
}
model.age.incidence.series.age <- function(par.dat, age_vect, year.start, age_mult){
  
  # you will have one lambda for each year in the time series
  # and one N-sero for each year as well
  # and X ages with varying probabilities of infection within those years
  
  lambda.start = par.dat$lambda[1]
  
  
  
  lts = nrow(par.dat)
  
  #proportion exposed in any way - needs to long enough for each age within each year
  pexposed= rep(NA, (length(age_vect)-1))
  pexposed[1] <- 0
  pexposed = rep(list(pexposed),lts) 
  
  #naive
  pnaive= rep(NA, (length(age_vect)-1))
  pnaive[1] <- 1
  pnaive = rep(list(pnaive),lts) 
  
  
  #exposed to a single serotype only (this is one serotype's worth of primary infections)
  pprim= rep(NA, (length(age_vect)-1))
  pprim[1] <- 0
  pprim = rep(list(pprim),lts)
  
  #multitypic infections
  pmulti= rep(NA, (length(age_vect)-1))
  pmulti[1] <- 0
  pmulti = rep(list(pmulti),lts)
  
  age_tracker= rep(NA, (length(age_vect)-1))
  age_tracker[1] <- 0
  age_tracker = rep(list(age_tracker),lts)
  
  year_tracker= rep(NA, (length(age_vect)-1))
  year_tracker[1] <- year.start
  year_tracker = rep(list(year_tracker),lts)
  
  year.start= min(par.dat$year)
  
  
  for (i in 1:lts){
    
    if(i < max(age_vect)){
      age_vect_tmp = age_vect[1:which(age_vect==i)]
    }else{
      age_vect_tmp = age_vect
    }
    
    for(a in 2:(length(age_vect_tmp))){ 
      # for-loops over all possible ages in our data range across all the years. 
      # all ages get rounded up, so babies <1 year are counted as 1 year old,
      # and experienced a hazard of infection in the year prior
      
      
      
      
      age = age_vect_tmp[a]
      
      year.now = par.dat$year[i]
      year.par = subset(par.dat, year<=year.now)
      N_sero = year.par$N_sero
      lambda =  year.par$lambda #one per year
      
      # make a vector of durations for the lambdas across the time series
      # in this model, it is easy because the duration is only 1 year...
      # you may eventually also want partial year durations in a different model
      #if (age>1){
      dur = rep(1, length(lambda))  
      #}else if(age<=1){
      # dur = rep(0, length(lambda))  
      #}
      
      #if you were born after the time series began
      #some of those ones may need to be replaced with 0s
      if(i>1){
        #diff.count <- i-age
        tot.class = ceiling(age)
        #diff.count <- sum(dur)-age_trunc
        
        if(tot.class<length(dur)){
          dur[1:(length(dur)-tot.class)]<- 0
        }
      }
      
      # #if 0, this kid spent 1 year in each age class and the dur vector does not need to be changed
      if(sum(dur) != age){
        warning(paste0("mismatch in age and integration at age=", age, " and i=", i))
      }
      
      
      #now, add in age-specific lambda multiplier for values of lambda for which dur is not 0
      # also add in 10 age-specific multiplicative factors for ages:
      
      #1-2,3-4,5-6,7-9,10-12,13-15,16-19,20-29,30-39,40+
      
      if(year.now <=2010){
        age.mult.df <- cbind.data.frame(mult=age_mult[1:8], 
                                        age_min=c(1,3,5,7,10,13,16,20),
                                        age_max=c(2,4,6,9,12,15,19,90))  
      }else if (year.now >2010){
        age.mult.df <- cbind.data.frame(mult=age_mult[9:19], 
                                        age_min=c(1,3,5,7,10,13,16,20,30,40,50),
                                        age_max=c(2,4,6,9,12,15,19,29,39,49,90))  
      }
      
      
      
      age.mult.df$dur = (age.mult.df$age_max-age.mult.df$age_min)+1
      
      #rep it as a multiplier by year
      age.mult.vect = c(unlist(mapply(rep, x=as.list(age.mult.df$mult), each=as.list(age.mult.df$dur), SIMPLIFY = F)))
      
      #now, make it the length of dur for which there are not 0s
      
      age.mult.vect = age.mult.vect[1:length(dur[dur>0])]
      
      age.mult.vect = c(dur[dur==0], age.mult.vect)
      
      if(length( age.mult.vect) != length(dur)){
        warning(paste0("age multiplier incorrect length at age=", age, " and i=", i))
      }
      
      
      #just to be safe, convert the 0s to 1s so it does not modulate anything
      age.mult.vect[age.mult.vect==0] <- 1
      
      # now, modulate lambda by multiplying with age.mult.vect
      lambda.age = lambda*age.mult.vect
      
      # now use lambda.age below
      
      # now, we integrate the hazard of becoming infected over all the years in each
      
      # First, integrand A deals with avoiding infection with all strains until now
      # Then, integrand B deals with getting exposed to that select strain.
      
      # first, integrand A  (avoiding infection with all other strains) 
      inte_A = sum(dur*lambda.age*N_sero)
      
      # them integrand B (getting exposed to a single strain)
      inte_B = sum(dur*lambda.age)
      
      # now fill in the probabilities based on these : 
      #this is the probability of being naive to all serotypes at this point in time
      pnaive[[i]][[a]] <- exp(-inte_A)
      
      #this is the probability of being exposed to any serotype at the current time
      pexposed[[i]][[a]] <- (1-exp(-inte_A))
      
      # this is the probability of being exposed to just one serotype (a primary infection)
      # here, this is expressed for just a single target strain 
      # (would need to multiply by the number of circulating strains if you wanted 
      # the probability of a primary infection with ANY strain)
      pprim[[i]][[a]] <- exp(-inte_A)*(exp(inte_B)-1)
      
      
      # if not primarily infected or naive, this should be a multitypic infection
      # but remember it could be a primary infection with any of the four serotypes 
      # so you need to account for that here by again multiplying by the 
      # number of circulating serotypes internally
      pmulti[[i]][[a]] = 1 - pnaive[[i]][[a]] - (N_sero[i]*(pprim[[i]][[a]]))
      age_tracker[[i]][[a]] = age
      year_tracker[[i]][[a]] =year.now
      
    }
  }
  
  #some of the age distributions will be shorter
  
  #and get the estimates of each
  #p.out = cbind.data.frame(year = rep(seq((year.start), (year.start+lts-1), 1), each=length(age_vect)),
  p.out <-  cbind.data.frame(year=c(unlist(year_tracker)),
                             age=c(c(unlist(age_tracker))), 
                             exposed=c(c(unlist(pexposed))),
                             naive=c(c(unlist(pnaive))),
                             one_prim=c(unlist(pprim)), 
                             multi = c(unlist(pmulti)))
  
  add.df <- dplyr::select(par.dat, year, N_sero)
  p.out <- merge(p.out, add.df, by=c("year"), all.x = T, sort=FALSE)
  
  p.out$all_prim <- p.out$one_prim*p.out$N_sero
  
  
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(cbind(p.out["exposed"], p.out["naive"])) #should be 1
  p.out$sum_naive_prim_multi <- rowSums(cbind(p.out["naive"], p.out["all_prim"], p.out["multi"])) #should be 1
  
  
  p.out <- arrange(p.out, year, age)
  
  #ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  +  geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  p.out$cum_prop_cases <- p.out$multi
  p.out <- p.out[complete.cases(p.out),]
  
  
  p.add <- cbind.data.frame(year = unique(p.out$year), age = rep(0, length(unique(p.out$year))), exposed=0, naive=1, one_prim=0, multi=0, N_sero = N_sero[1], all_prim=0, sum_exp_naive=1, sum_naive_prim_multi=1, cum_prop_cases =0)
  
  
  p.sum <- rbind(p.add, p.out)
  p.sum <- arrange(p.sum, year, age)
  p.sum$provname = unique(par.dat$provname)
  
  #p1 <- ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  # print(p1)
  # 
  #ggplot(data=subset(p.sum, year==1985)) + geom_point(aes(x=age, y=multi))  + 
  #geom_line(aes(x=age, y=multi)) 
  
  return(p.sum) #returns prevalence by age for each year for fitting to the years for which we have data
}
log.lik.age <- function(par, par.dat, dat){
  
  #par.dat$lambda does not change here
  age.mult <- par
  
  age_vect=seq(0, max(dat$age), by=1)
  
  
  #make year start split
  year.start <- min(par.dat$year)
  
  out.mod <- model.age.incidence.series.age(par.dat= par.dat, year.start=year.start, age_mult=age.mult, age_vect=age_vect)
  
  
  #ggplot(data=out.mod) + geom_point(aes(x=age,y= cum_prop_cases)) + facet_wrap(~year)
  
  # now, select only the years for fitting for which there are data
  # (the model can produce data back further in time based on the first year FOI from the dataset)
  out.mod = subset(out.mod, year >=min(dat$year))
  
  
  # merge model with data
  out.mod <- arrange(out.mod, year, age)
  dat <- arrange(dat,year, age)
  dat.merge <- dplyr::select(dat, age, year, Nage, cum_cases, n)
  out.merge <- merge(out.mod, dat.merge, by= c("year", "age"))
  out.merge$cum_prop_cases_data <- out.merge$cum_cases/out.merge$n
  
  out.merge <- arrange(out.merge, year, age)
  
  # #plot model with data
  # ggplot(data=out.merge)  + facet_wrap(~year) +  geom_point(aes(x=age, y=cum_prop_cases_data)) + geom_line(aes(x=age, y=cum_prop_cases_data))+geom_line(aes(x=age,y= cum_prop_cases), color="tomato")
  # #   
  
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(out.merge$age)){
    ll=ll+dbinom(out.merge$cum_cases[i],out.merge$n[i],prob=out.merge$cum_prop_cases[i],log=T)
    #print(ll)
  }
  
  # 
  
  par.dat$lambda_high <- 0 
  par.dat$lambda_high[par.dat$lambda>1] <- 1 #penalty if lambda is >1
  tot_over = sum(par.dat$lambda_high)
  
  if(tot_over>0){
    
    ll <- ll-(tot_over*1000000)
    
  }
  
  if(ll==-Inf){
    ll <- -1000000
  }
  
  return(-ll)
}
fit.all.LBFGSB.three <- function(dat.all,  burnin, lambda.guess, N.sero.fix, sim_type, fit.CI){
  
  
  #and only take those with values
  dat3 = subset(dat.all, count>0 & age < max(dat.all$age))
  
  dat = subset(dat3,  year >= (min(dat.all$year)+ burnin))
  
  
  
  #sum by age an
  
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  
  year.dat.sum <- lapply(year.dat, sum.yr.all.sim)
  df.out <- data.table::rbindlist( year.dat.sum)
  
  #head(df.out)
  #ggplot(df.out) + geom_line(aes(x=age,y=cum_prop_cases)) + facet_wrap(~year)
  
  #get dist back here for the oldest individual in the first year of the dataset
  dist.back =   max(df.out$age[df.out$year==min(df.out$year) & df.out$Nage>0]) 
  
  
  # # # # # # # # #head(df.out)
  #              ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
  #                 geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  # # # # # # # # # # # # # # # 
  # # # # #make your guess parameters
  # #lambda is takes data from the previous year and creates infections in this year
  
  #now, make the appropriate number of serotypes and lambdas
  
  if(length(N.sero.fix)==1 & length(lambda.guess)==1){ #here, number of serotypes is fixed across the time series
    par.dat <- cbind.data.frame(year= ((min(dat$year)-dist.back +1):max(dat$year)),
                                lambda = rep(lambda.guess, length((min(dat$year)-dist.back +1):max(dat$year))),
                                N_sero = rep(N.sero.fix, length((min(dat$year)-dist.back +1):max(dat$year))))
    
  }else if (length(N.sero.fix)>1 & length(lambda.guess)==1){ #here you can vary the sero-strains by providing your own vector
    par.dat <- cbind.data.frame(year= ((min(dat$year)-dist.back +1):max(dat$year)),
                                lambda = rep(lambda.guess, length((min(dat$year)-dist.back +1):max(dat$year))),
                                N_sero = N.sero.fix)
    
  }else if (length(N.sero.fix)>1 & length(lambda.guess)>1){ #here you have pre-prepped both the lambda vector and the serotype vector
    par.dat <- cbind.data.frame(year= ((min(dat$year)-dist.back +1):max(dat$year)),
                                lambda = lambda.guess,
                                N_sero = N.sero.fix)
    
    
  }else if (length(N.sero.fix)==1 & length(lambda.guess)>1){ # here you have pre-prepped only the lambda but are going to fix serotypes across the dataset
    par.dat <- cbind.data.frame(year= ((min(dat$year)-dist.back +1):max(dat$year)),
                                lambda = lambda.guess,
                                N_sero =rep(N.sero.fix, length((min(dat$year)-dist.back +1):max(dat$year))))
    
  }
  
  #now, for all lambda prior to 1970, assume a very, very, very low FOI
  #par.dat$lambda[par.dat$year<1970] <- 0.0000001
  #and fit it cumulatively across the entire time series
  
  
  #if this is the first year, you need
  #log.lambda.guess <- log(par.dat$lambda)
  lambda.guess <- par.dat$lambda
  
  
  out.NS <- optim(par = lambda.guess, 
                  fn=log.lik.fit.all, 
                  method = "L-BFGS-B",
                  lower = 0.0000001,
                  upper=1,
                  par.dat=par.dat, 
                  year.start = min(par.dat$year),#this is the year the model will start iterating with the first birth cohort
                  dat=df.out, 
                  control = list(maxit=1000),
                  hessian = T)
  
  
  # 
  # out.NS <- nlm(f=log.lik.fit.all, 
  #               p = log.lambda.guess, 
  #               par.dat=par.dat, 
  #               year.start = min(par.dat$year),#this is the year the model will start iterating with the first birth cohort
  #               dat=df.out, 
  #               iterlim=1000,
  #               hessian = T)
  
  
  # par.dat$lambda <- exp(out.NS$estimate)
  # par.dat$llik <- out.NS$minimum
  # par.dat$convergence <- out.NS$code
  
  
  
  par.dat$lambda <- out.NS$par
  par.dat$llik <- out.NS$value
  par.dat$convergence <- out.NS$convergence
  
  if(fit.CI==TRUE){
    
    if (is.positive.definite(out.NS$hessian)==TRUE){
      hess <- solve(out.NS$hessian)
      prop_sigma <-sqrt(diag(hess))
      upper<-out.NS$par+1.96*prop_sigma
      lower<-out.NS$par-1.96*prop_sigma
      CI <-data.frame(lower=lower, upper=upper)
      CI$lower[CI$lower<0] <- 0
      CI$upper[CI$upper<0] <- 0
      
      par.dat$lci <- CI$lower
      par.dat$uci <- CI$upper
      
    }else{
      par.dat$lci <- "not yet"
      par.dat$uci <- "not yet"
      
    }
    
  }
  
  #and return
  par.dat$sim_type <- sim_type
  
  return(par.dat)
  
}
fit.age.LBFGSB.all.three <- function(dat.all,par.dat, burnin, fit.CI, sim_type){
  
  #and only take those with values
  dat3 = subset(dat.all, count>0 & age < max(dat.all$age))
  
  dat = subset(dat3,  year >= (min(dat.all$year)+ burnin))
  
  
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  
  year.dat.sum <- lapply(year.dat, sum.yr.all.sim)
  df.out <- data.table::rbindlist( year.dat.sum)
  
  #head(df.out)
  #now take a guess at 10 age-specific multiplicative factors for 3 different eras:
  age.guess <- rep(1,19)
  
  
  out.NS <- optim(par = age.guess, 
                  fn=log.lik.age, 
                  method = "L-BFGS-B",
                  lower = 0.0000001, # slowest sigma
                  upper=1, # fastest sigma - not allowing for waning more rapid than a single year
                  par.dat=par.dat,
                  dat=df.out,
                  control = list(maxit=1000),
                  hessian = T)
  
  
  #lambda is already fixed
  
  #return sigma as its own data table
  
  
  age.mult.df <- cbind.data.frame(age_mult=c(out.NS$par), 
                                  age_min=c(1,3,5,7,10,13,16,20, 1,3,5,7,10,13,16,20,30,40,50),
                                  age_max=c(2,4,6,9,12,15,19,90, 2,4,6,9,12,15,19,29,39,49,90),
                                  year_range = c(rep("2002-2010",8), rep("2011-2020", 11)))
  
  
  
  age.mult.df$llik <- out.NS$value
  age.mult.df$convergence <- out.NS$convergence 
  
  # will need to rerun the subsets of the data with these universal
  # age modifiers to get likelihood values for those
  
  
  
  if(fit.CI==TRUE){
    
    if (is.positive.definite(out.NS$hessian)==TRUE){
      hess <- solve(out.NS$hessian)
      prop_sigma <-sqrt(diag(hess))
      upper<-out.NS$par+1.96*prop_sigma
      lower<-out.NS$par-1.96*prop_sigma
      CI <-data.frame(lower=lower, upper=upper)
      CI$lower[CI$lower<0] <- 0
      CI$upper[CI$upper<0] <- 0
      
      
      age.mult.df$lci_mult <- CI$lower
      age.mult.df$uci_mult <- CI$upper
      
    }else{
      age.mult.df$lci_mult <- "not yet"
      age.mult.df$uci_mult <- "not yet"
      
    }
    
  }
  
  
  #and return
  age.mult.df$sim_type <- sim_type
  return(age.mult.df)
  
}
log.lik.wane.fit.all <- function(par, par.dat, dat){
  
  #par.dat$lambda does not change here
  #sigma gets added to the rest of the parameters
  par.dat$sigma <- NA
  par.dat$sigma[par.dat$year<min(dat$year)] <- 0
  
  sigma.df = cbind.data.frame(year=min(dat$year):max(dat$year), sigma=exp(par))
  
  for(i in 1:nrow(sigma.df)){
    par.dat$sigma[par.dat$year==sigma.df$year[i]]<- sigma.df$sigma[i]
  }
  
  
  age_vect=seq(0, max(dat$age), by=1)
  
  
  
  out.mod  <-model.age.incidence.series.wane(par.dat= par.dat, year.start=min(par.dat$year), age_vect=age_vect)
  
  #ggplot(data=out.mod) + geom_line(aes(x=age,y= cum_prop_cases, color=provname)) + facet_wrap(~year)
  
  # now, select only the years for fitting for which there are data
  # (the model can produce data back further in time based on the first year FOI from the dataset)
  out.mod = subset(out.mod, year >=min(dat$year))
  
  # merge model with data
  out.mod <- arrange(out.mod, year, age)
  dat <- arrange(dat, year, age)
  dat.merge <- dplyr::select(dat,  age, year, Nage, cum_cases, n)
  out.merge <- merge(out.mod, dat.merge, by= c("year", "age"))
  out.merge$cum_prop_cases_data <- out.merge$cum_cases/out.merge$n
  
  out.merge <- arrange(out.merge, year, age)
  
  
  
  
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(out.merge$age)){
    ll=ll+dbinom(out.merge$cum_cases[i],out.merge$n[i],prob=out.merge$cum_prop_cases[i],log=T)
    #print(ll)
  }
  
  
  
  if(ll==-Inf){
    ll <- -1000000
  }
  
  return(-ll)
}
model.age.incidence.series.wane <- function(par.dat, age_vect, year.start){
  
  # you will have one lambda for each year in the time series
  # and one N-sero for each year as well
  # and X ages with varying probabilities of infection within those years
  
  lambda.start = par.dat$lambda[1]
  
  
  
  lts = nrow(par.dat)
  
  #proportion exposed in any way - needs to long enough for each age within each year
  pexposed= rep(NA, (length(age_vect)-1))
  pexposed[1] <- 0
  pexposed = rep(list(pexposed),lts) 
  
  #naive
  pnaive= rep(NA, (length(age_vect)-1))
  pnaive[1] <- 1
  pnaive = rep(list(pnaive),lts) 
  
  
  #exposed to a single serotype only (this is one serotype's worth of primary infections)
  pprim= rep(NA, (length(age_vect)-1))
  pprim[1] <- 0
  pprim = rep(list(pprim),lts)
  
  #multitypic infections
  pmulti= rep(NA, (length(age_vect)-1))
  pmulti[1] <- 0
  pmulti = rep(list(pmulti),lts)
  
  age_tracker= rep(NA, (length(age_vect)-1))
  age_tracker[1] <- 0
  age_tracker = rep(list(age_tracker),lts)
  
  year_tracker= rep(NA, (length(age_vect)-1))
  year_tracker = rep(list(year_tracker),lts)
  year_tracker[[1]][1] <- year.start
  
  year.start= min(par.dat$year)
  
  
  for (i in 1:lts){
    
    if(i < max(age_vect)){
      age_vect_tmp = age_vect[1:which(age_vect==i)]
    }else{
      age_vect_tmp = age_vect
    }
    
    for(a in 2:(length(age_vect_tmp))){ 
      # for-loops over all possible ages in our data range across all the years. 
      # all ages get rounded up, so babies <1 year are counted as 1 year old,
      # and experienced a hazard of infection in the year prior
      
      
      
      
      age = age_vect_tmp[a]
      
      year.now = par.dat$year[i]
      year.par = subset(par.dat, year<=year.now)
      N_sero = year.par$N_sero
      lambda =  year.par$lambda #one per year
      sigma = year.par$sigma
      # make a vector of durations for the lambdas across the time series
      # in this model, it is easy because the duration is only 1 year...
      # you may eventually also want partial year durations in a different model
      #if (age>1){
      dur = rep(1, length(lambda))  
      #}else if(age<=1){
      # dur = rep(0, length(lambda))  
      #}
      
      #if you were born after the time series began
      #some of those ones may need to be replaced with 0s
      if(i>1){
        #diff.count <- i-age
        tot.class = ceiling(age)
        #diff.count <- sum(dur)-age_trunc
        
        if(tot.class<length(dur)){
          dur[1:(length(dur)-tot.class)]<- 0
        }
      }
      
      # #if 0, this kid spent 1 year in each age class and the dur vector does not need to be changed
      if(sum(dur) != age){
        warning(paste0("mismatch in age and integration at age=", age, " and i=", i))
      }
      
      
      
      # now, we integrate the hazard of becoming infected over all the years in each
      
      
      
      #First, integrand A (hazard of exposure to all strains )
      inte_A = (sum(dur*lambda*(N_sero)))
      
      #Then, hazard of exposure to strain i only
      inte_B = (sum(dur*lambda))
      
      # Then, probability of exposure to strain i and then strain k (where k does not equal i)
      # Does not require avoiding exposure with any other strains
      inte_C = (sum(dur*lambda*(2)))
      
      
      #Then, cumulative waning immunity hazard for one strain
      inte_E = sum(sigma*dur) 
      
      #And cumulative waning hazard for all the other strains
      inte_D = (sum((N_sero-1)*sigma*dur))
      
      
      
      #could also write as sigma[1]*(a-1) (since a starts at class 2 we need to subtract 1 from a)
      
      
      # now fill in the probabilities based on these : 
      #this is the probability of being naive to all serotypes at this point in time
      pnaive[[i]][[a]] <- exp(-inte_A)
      
      #this is the probability of being exposed to any serotype at the current time
      pexposed[[i]][[a]] <- (1-exp(-inte_A))
      
      # this is the probability of being exposed to just one serotype (a primary infection)
      # here, this is expressed for just a single target strain 
      # (would need to multiply by the number of circulating strains if you wanted 
      # the probability of a primary infection with ANY strain)
      pprim[[i]][[a]] <- exp(-inte_A)*(exp(inte_B)-1) + (1-exp(-inte_C))*(exp(-inte_D))*(1-exp(-inte_E))
      
      
      # if not primarily infected or naive, this should be a multitypic infection
      # but remember it could be a primary infection with any of the four serotypes 
      # so you need to account for that here by again multiplying by the 
      # number of circulating serotypes internally
      pmulti[[i]][[a]] = 1 - pnaive[[i]][[a]] - (N_sero[i]*(pprim[[i]][[a]]))
      age_tracker[[i]][[a]] = age
      year_tracker[[i]][[a]] =year.now
      
    }
  }
  
  #some of the age distributions will be shorter
  
  #and get the estimates of each
  #p.out = cbind.data.frame(year = rep(seq((year.start), (year.start+lts-1), 1), each=length(age_vect)),
  p.out <-  cbind.data.frame(year=c(unlist(year_tracker)),
                             age=c(c(unlist(age_tracker))), 
                             exposed=c(c(unlist(pexposed))),
                             naive=c(c(unlist(pnaive))),
                             one_prim=c(unlist(pprim)), 
                             multi = c(unlist(pmulti)))
  
  
  # merge to add N_sero
  add.par <- dplyr::select(par.dat, year, N_sero)
  
  p.out <- merge(p.out, add.par, by="year", all.x = T, sort = F)
  
  
  p.out$all_prim <- p.out$one_prim*p.out$N_sero
  
  
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(cbind(p.out["exposed"], p.out["naive"])) #should be 1
  p.out$sum_naive_prim_multi <- rowSums(cbind(p.out["naive"], p.out["all_prim"], p.out["multi"])) #should be 1
  
  
  p.out <- arrange(p.out, year, age)
  
  # ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  +  geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  p.out$cum_prop_cases <- p.out$multi
  p.out <- p.out[complete.cases(p.out),]
  
  
  p.add <- cbind.data.frame(year = unique(p.out$year), age = rep(0, length(unique(p.out$year))), exposed=0, naive=1, one_prim=0, multi=0, N_sero = N_sero[1], all_prim=0, sum_exp_naive=1, sum_naive_prim_multi=1, cum_prop_cases =0)
  
  
  p.sum <- rbind(p.add, p.out)
  p.sum <- arrange(p.sum, year, age)
  p.sum$provname = unique(par.dat$provname)
  
  
  p.sum <- p.sum[!duplicated(p.sum),]
  
  p1 <- ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + 
    geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  # print(p1)
  # 
  #ggplot(data=subset(p.sum, year==1985)) + geom_point(aes(x=age, y=multi))  + 
  #geom_line(aes(x=age, y=multi)) 
  
  return(p.sum) #returns prevalence by age for each year for fitting to the years for which we have data
}
fit.wane.LBFGSB.all.three <- function(dat.all,par.dat, burnin, sigma.guess, fit.CI, sim_type){
  
  
  #and only take those with values
  dat3 = subset(dat.all, count>0 & age < max(dat.all$age))
  
  dat = subset(dat3,  year >= (min(dat.all$year)+ burnin))
  
  
  
  #
  
  # par.dat <- cbind.data.frame(year= ((min(dat$year)-dist.back +1):max(dat$year)),
  #                             lambda = rep(lambda.fix, length((min(dat$year)-dist.back +1):max(dat$year))),
  #                             N_sero = rep(N.sero.fix, length((min(dat$year)-dist.back +1):max(dat$year))))
  # 
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  
  year.dat.sum <- lapply(year.dat, sum.yr.all.sim)
  df.out <- data.table::rbindlist( year.dat.sum)
  #ggplot(df.out) +geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  
  out.NS <- optim(par = log(rep(sigma.guess, length((min(df.out$year):max(df.out$year))))), 
                  fn=log.lik.wane.fit.all, 
                  method = "L-BFGS-B",
                  lower = rep((-16.2), length((min(dat$year):max(dat$year)))), # slowest sigma
                  upper=rep((-0.69), length((min(dat$year):max(dat$year)))),
                  par.dat=par.dat,
                  dat=df.out,
                  control = list(maxit=1000),
                  hessian = T)
  
  
  #lambda is already fixed
  
  #return sigma as its own data table
  
  sigma.df <-cbind.data.frame(year = ((min(df.out$year):max(df.out$year))), sigma=out.NS$par, dur_immunity = (1/out.NS$par))
  sigma.df$neg_llik=out.NS$value
  sigma.df$convergence = out.NS$convergence
  
  if(fit.CI==TRUE){
    
    if (is.positive.definite(out.NS$hessian)==TRUE){
      hess <- solve(out.NS$hessian)
      prop_sigma <-sqrt(diag(hess))
      upper<-out.NS$par+1.96*prop_sigma
      lower<-out.NS$par-1.96*prop_sigma
      CI <-data.frame(lower=lower, upper=upper)
      CI$lower[CI$lower<0] <- 0
      CI$upper[CI$upper<0] <- 0
      
      
      sigma.df$lci_sigma <- CI$lower
      sigma.df$uci_sigma <- CI$upper
      
      sigma.df$uci_dur_imm <- 1/sigma.df$lci_sigma 
      sigma.df$lci_dur_imm <- 1/sigma.df$uci_sigma 
      
      
    }else{
      sigma.df$lci_sigma <- NA
      sigma.df$uci_sigma <- NA
      
      sigma.df$lci_dur_imm <- NA
      sigma.df$uci_dur_imm <- NA
      
    }
    
  }
  
  #and return
  sigma.df$sim_type <- sim_type
  return(sigma.df)
  
}



#load the foi fits for cambodia
# fit.dat <- read.csv(file = "prov-fits-FOI.csv", stringsAsFactors = F, header = T)
# nrow(fit.dat[fit.dat$provname=="National",]) #40 years - run for 60 before this

# #now fit and recover age -FOI fixed
# par.dat <- cbind.data.frame(year = 1921:2020,lambda=c(rep(.2,60),fit.dat$lambda[fit.dat$provname=="National"]))
# par.dat$N_sero=2
# par.dat = subset(par.dat, year>=2000)


load("hyp1-fit-lambda.Rdata")
load("comp-dat-sim.Rdata")

cam.sim = subset(comp.dat, hyp=="H1: Increasing Tertiary\nCase Detection")


hyp1.fit.wane <- fit.wane.LBFGSB.all.three(dat.all = cam.sim,
                                           burnin=20,
                                           par.dat=hyp1.fit.lambda,
                                           sigma.guess = 1/100,
                                           sim_type="hyp1_wane", 
                                           fit.CI=T)

save(hyp1.fit.wane, file = "hyp1-fit-wane.Rdata")
