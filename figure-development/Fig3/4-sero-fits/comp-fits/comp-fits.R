

rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(lubridate)
library(matrixcalc)

# at the province level, compare each run with and without age data
# and with and without waning immunity

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig3/comp-fits/"))

#functions
sum.yr.all.prov <- function(df){
  
  df.sum <- ddply(df, .(year, age), summarise, Nage = length(age))
  
  
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
  df.out$provname = unique(df$provname)
  return(df.out)
  
}
sum.yr.all <- function(df){
  
  df.sum <- ddply(df, .(year, age), summarise, Nage = length(age))
  
  
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
sum.yr.yr <- function(df, age_vect){
  
  df.sum <- ddply(df, .(year, age), summarise, Nage = length(age))
  
  
  df.out = cbind.data.frame(age=1:max(age_vect))
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
get.min.year <- function(df){
  out = min(df$year)
  return(out)
}
log.lik.wane.fit.one <- function(par.dat, dat){
  
  age_vect=seq(0, max(dat$age), by=1)
  
  out.mod <- model.age.incidence.series.wane(par.dat= par.dat, year.start=min(par.dat$year), age_vect=age_vect)
  
  #ggplot(data=out.mod) + geom_line(aes(x=age,y= cum_prop_cases, color=provname)) + facet_wrap(~year)
  
  # now, select only the years for fitting for which there are data
  # (the model can produce data back further in time based on the first year FOI from the dataset)
  out.mod = subset(out.mod, year >=min(dat$year))
  
  # merge model with data
  out.mod <- arrange(out.mod,  year, age)
  dat <- arrange(dat, year, age)
  dat.merge <- dplyr::select(dat, age, year, Nage, cum_cases, n)
  out.merge <- merge(out.mod, dat.merge, by= c("year", "age"))
  out.merge$cum_prop_cases_data <- out.merge$cum_cases/out.merge$n
  
  out.merge <- arrange(out.merge, year, age)
  
  # ggplot(data=out.merge) + 
  #   facet_wrap(~year) +  geom_point(aes(x=age, y=cum_prop_cases_data)) +
  #   geom_line(aes(x=age, y=cum_prop_cases_data), size=1) +
  #   geom_line(aes(x=age,y= cum_prop_cases), color="tomato") 
  # 
  # 
  
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(out.merge$age)){
    ll=ll+dbinom(out.merge$cum_cases[i],out.merge$n[i],prob=out.merge$cum_prop_cases[i],log=T)
    #print(ll)
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
      sigma = rep(unique(par.dat$sigma), length(lambda))
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
      
      #Then, cumulative waning immunity hazard for one strain
      inte_C = sum(sigma*dur) 
      
      #And cumulative waning hazard for all the others
      inte_D = (sum((N_sero-1)*sigma*dur))
      
      #Then, probability of exposure to any other strains besides i
      inte_E = (sum(dur*lambda*(N_sero-1)))
      
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
      pprim[[i]][[a]] <- exp(-inte_A)*(exp(inte_B)-1) + (exp(-inte_A)*(exp(inte_B)-1)*(1-exp(-inte_E))*(1-exp(-inte_C))*(exp(-inte_D)))
      
      
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
  
  
  
  p.sum <- p.sum[!duplicated(p.sum),]
  
  p1 <- ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + 
    geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  # print(p1)
  # 
  #ggplot(data=subset(p.sum, year==1985)) + geom_point(aes(x=age, y=multi))  + 
  #geom_line(aes(x=age, y=multi)) 
  
  return(p.sum) #returns prevalence by age for each year for fitting to the years for which we have data
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
  year_tracker = rep(list(year_tracker),lts)
  year_tracker[[1]][1] <- year.start
  
  
  
  
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
     # if(i>1){
        #diff.count <- i-age
        tot.class = ceiling(age)
        #diff.count <- sum(dur)-age_trunc
        
        if(tot.class<length(dur)){
          dur[1:(length(dur)-tot.class)]<- 0
        }
    #  }
      
      # #if 0, this kid spent 1 year in each age class and the dur vector does not need to be changed
      if(sum(dur) != age){
        warning(paste0("mismatch in age and integration at age=", age, " and i=", i))
      }
      
      
      #now, add in age-specific lambda multiplier for values of lambda for which dur is not 0
      # also add in 10 age-specific multiplicative factors for ages:
      
      #1-2,3-4,5-6,7-9,10-12,13-15,16-19,20-29,30-39,40+
      
      age.mult.df <- cbind.data.frame(mult=age_mult, 
                                      age_min=c(1,3,5,7,10,13,16,20,30,40,50,60,70,80),
                                      age_max=c(2,4,6,9,12,15,19,29,39,49,59,69,79,90))
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

  
  #p1 <- ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  # print(p1)
  # 
  #ggplot(data=subset(p.sum, year==1985)) + geom_point(aes(x=age, y=multi))  + 
  #geom_line(aes(x=age, y=multi)) 
  
  return(p.sum) #returns prevalence by age for each year for fitting to the years for which we have data
}
log.lik.wane.age.fit.all <- function(age.mult, par.dat, dat){
  
  
  
  age_vect=seq(0, max(dat$age), by=1)
  
  
  
  
  #make year start split
  year.start = min(par.dat$year)
  
  out.mod <-model.age.incidence.series.wane(par.dat= par.dat, year.start=year.start, age_mult=age.mult, age_vect=age_vect)
  #out.mod <- data.table::rbindlist(out.mod.list)
  #ggplot(data=out.mod) + geom_line(aes(x=age,y= cum_prop_cases, color=provname)) + facet_wrap(~year)
  
  # now, select only the years for fitting for which there are data
  # (the model can produce data back further in time based on the first year FOI from the dataset)
  out.mod = subset(out.mod, year >=min(dat$year))
  
  # merge model with data
  out.mod <- arrange(out.mod,  year, age)
  dat <- arrange(dat,  year, age)
  dat.merge <- dplyr::select(dat,age, year, Nage, cum_cases, n)
  out.merge <- merge(out.mod, dat.merge, by= c("year", "age"))
  out.merge$cum_prop_cases_data <- out.merge$cum_cases/out.merge$n
  
  out.merge <- arrange(out.merge, year, age)
  
  #  ggplot(data=out.merge) + 
  #    facet_wrap(~year) +  geom_point(aes(x=age, y=cum_prop_cases_data)) +
  #    geom_line(aes(x=age, y=cum_prop_cases_data), size=1) +
  #    geom_line(aes(x=age,y= cum_prop_cases), color="tomato") 
  # # 
  
  
  
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(out.merge$age)){
    ll=ll+dbinom(out.merge$cum_cases[i],out.merge$n[i],prob=out.merge$cum_prop_cases[i],log=T)
    #print(ll)
  }
  
  
  return(-ll)
}
model.age.incidence.series.wane <- function(par.dat, age_vect, year.start, age_mult){
  
  
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
      sigma = rep(unique(par.dat$sigma), length(lambda))
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
      
      
      #1-2,3-4,5-6,7-9,10-12,13-15,16-19,20-29,30-39,40+
      
      age.mult.df <- cbind.data.frame(mult=age_mult, 
                                      age_min=c(1,3,5,7,10,13,16,20,30,40,50,60,70,80),
                                      age_max=c(2,4,6,9,12,15,19,29,39,49,59,69,79,90))
      
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
      
      
      # now, we integrate the hazard of becoming infected over all the years in each
      
      #First, integrand A (hazard of exposure to all strains )
      inte_A = (sum(dur*lambda.age*(N_sero)))
      
      #Then, hazard of exposure to strain i only
      inte_B = (sum(dur*lambda.age))
      
      #Then, cumulative waning immunity hazard for one strain
      inte_C = sum(sigma*dur) 
      
      #And cumulative waning hazard for all the others
      inte_D = (sum((N_sero-1)*sigma*dur))
      
      #Then, probability of exposure to any other strains besides i
      inte_E = (sum(dur*lambda.age*(N_sero-1)))
      
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
      pprim[[i]][[a]] <- exp(-inte_A)*(exp(inte_B)-1) + (exp(-inte_A)*(exp(inte_B)-1)*(1-exp(-inte_E))*(1-exp(-inte_C))*(exp(-inte_D)))
      
      
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
log.lik.fit.age <- function(par.dat, age.mult, dat){
  
  
  age_vect=seq(0, max(dat$age), by=1)
  
  
  out.mod <- model.age.incidence.series.age(par.dat= par.dat, year.start=min(par.dat$year), age_mult=age.mult$age_mult, age_vect=age_vect)
  
  
  #ggplot(data=out.mod) + geom_point(aes(x=age,y= cum_prop_cases)) + facet_wrap(~year)
  
  # now, select only the years for fitting for which there are data
  # (the model can produce data back further in time based on the first year FOI from the dataset)
  out.mod = subset(out.mod, year >=min(dat$year))
  
  
  # merge model with data
  out.mod <- arrange(out.mod,  year, age)
  dat <- arrange(dat,  year, age)
  dat.merge <- dplyr::select(dat,  age, year, Nage, cum_cases, n)
  out.merge <- merge(out.mod, dat.merge, by= c( "year", "age"))
  out.merge$cum_prop_cases_data <- out.merge$cum_cases/out.merge$n
  
  out.merge <- arrange(out.merge,  year, age)
  
  # #plot model with data
  # ggplot(data=out.merge) + 
  #    facet_wrap(~year) +  geom_point(aes(x=age, y=cum_prop_cases_data)) +
  #    geom_line(aes(x=age, y=cum_prop_cases_data), size=1) +
  #    geom_line(aes(x=age,y= cum_prop_cases), color="tomato") 
  #    
  # #   
  
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(out.merge$age)){
    ll=ll+dbinom(out.merge$cum_cases[i],out.merge$n[i],prob=out.merge$cum_prop_cases[i],log=T)
    #print(ll)
  }
  
  
  
 
  return(-ll)
}
run.comp.models <- function(dat,par.dat, sigma, age.mult){
  
  #get dist back here for the oldest individual in the first year of the dataset
  #dist.back =   max(dat$age[dat$year==min(dat$year)]) 
  
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  
  year.dat.sum <- lapply(year.dat, sum.yr.all)
  df.out <- data.table::rbindlist( year.dat.sum)
  
  df.tmp = subset(df.out, year>2018)
  #you already have log-lik for the raw fit
  #now get it for age-fit
  par.dat$llik_age <- log.lik.fit.age(par.dat = par.dat, age.mult=age.mult, dat=df.out)
  
  
  #and sigma fit
  par.dat$sigma = sigma
  par.dat$llik_sigma = log.lik.wane.fit.one(par.dat = par.dat, dat=df.out)
  
  #and sigma combined with age
  par.dat$llik_sigma_age = log.lik.wane.age.fit.all(age.mult = age.mult$age_mult, par.dat = par.dat, dat=df.out)
  
  par.slim <- dplyr::select(par.dat, llik, llik_age, llik_sigma, provname)
  par.slim <- par.slim[!duplicated(par.slim),]
  
  #and calculate AIC for each
  #AIC = 2k - 2*log(lik)
  names(par.slim)[names(par.slim)=="llik"] <- "llik_FOI"
  par.slim$AIC_FOI <- (2*nrow(par.dat))+(2*par.slim$llik_FOI)
  par.slim$AIC_age <- (2*(nrow(par.dat)+10))+(2*par.slim$llik_age)
  par.slim$AIC_sigma <- (2*(nrow(par.dat)+1))+(2*par.slim$llik_sigma)
  
  
  return(par.slim)
  
}


#load the data, split by province, and run it to get llik with age added and with sigma added


dat <- read.csv(file = "DENV-Dist-Diag.csv", header=T, stringsAsFactors = F)
head(dat) 

#plot time series of each type by province by year
unique(dat$dianostic) #df, dhf, dss
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
head(dat)

#now attach biwks (2 biweeks per week)
wks = 1:52
biwks <- rep(1:26, each=2)
wk.biwk <- cbind.data.frame(week_report=wks, biwk=biwks)

sort(unique(dat$week_report))
unique(dat$provname)

dat <- merge(dat, wk.biwk, by="week_report")
#head(dat) #now plot as the lowest epiweek date per biweek

dat$age <- ceiling(dat$age)

#head(dat)
#min(dat$age)

dat <- arrange(dat, procode, date)

load("prov-fits-FOI.Rdata")
load("age-mult-refit.Rdata")
load("sigma-fit-all-time-var.Rdata")


dat.nat = dat
dat.nat$provname <- "National"

dat <- rbind(dat, dat.nat)
dat <- arrange(dat, provname, year)
fit.dat <- arrange(fit.dat, provname, year)

#split by prov
fit.split <- dlply(fit.dat, .(provname))
dat.split <- dlply(dat, .(provname))

#and apply comparison script
comp.split <- mapply(run.comp.models, dat=dat.split, par.dat=fit.split, MoreArgs = list(sigma = sigma.fit$sigma, age.mult = out), SIMPLIFY = F)

comp.df <- data.table::rbindlist(comp.split)
