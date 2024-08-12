#profile age modifiers

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

#load the age modifiers that need profiling


#and load the function to do so
# first, load the functions
model.age.incidence.series.age.mult.wane <- function(par.dat, age_vect, year.start, age_mult){
  
  
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
      
      
      #1-2,3-4,5-6,7-9,10-12,13-15,16-19,20-29,30-39,40+
      #1-2,3-4,5-6,7-9,10-12,13-15,16-19,20-29,30-39,40+
      if(year.now <=2010){
        age.mult.df <- cbind.data.frame(mult=age_mult[1:13], 
                                        age_min=c(1,3,5,7,10,13,16,20,30,40,50,60,70),
                                        age_max=c(2,4,6,9,12,15,19,29,39,49,59,69,90))  
      }else if (year.now >2010){
        age.mult.df <- cbind.data.frame(mult=age_mult[14:26], 
                                        age_min=c(1,3,5,7,10,13,16,20,30,40,50,60,70),
                                        age_max=c(2,4,6,9,12,15,19,29,39,49,59,69,90))  
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
get.min.year <- function(df){
  out = min(df$year)
  return(out)
}
run.model.data.all <- function(dat,par.dat, sigma.fit, age.mult.df){
  
  
  
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  year.start = min(par.dat$year)
  
  nat.dat <- dat
  nat.dat$provname  = "National"
  dat <- rbind(dat, nat.dat)
  
  #first, prep the data
  year.dat <- dlply(dat, .(provname, year))
  
  
  year.dat.sum <- lapply(year.dat, sum.yr.all.prov)
  
  df.out <- data.table::rbindlist( year.dat.sum)
  
  
  #and run the models across all the provinces
  
  age_vect=seq(0, max(dat$age), by=1)
  
  #merge sigma with par.dat 
  #and sigma fit only
  sigma.merge <- dplyr::select(sigma.fit, year, sigma)
  par.dat <- merge(par.dat, sigma.merge, by="year", all.x = T)
  par.dat$sigma[is.na(par.dat$sigma)] <- 0
  
  
  #now, split the par.dat by provinces and apply simultaneously
  par.dat.split <- dlply(par.dat, .(provname))
  
  #make year start split
  year.start.split <- lapply(par.dat.split, get.min.year)
  
  #apply function across all par and year start subsets
  # age mod and age vect are constant
  
  out.mod.list <- mapply(model.age.incidence.series.age.mult.wane, par.dat= par.dat.split, year.start=year.start.split, MoreArgs = list(age_mult = age.mult.df$age_mult, age_vect=age_vect), SIMPLIFY = FALSE)
  
  # and bind
  
  out.mod.df <- data.table::rbindlist(out.mod.list)
  
  #head(out.mod.df)
  out.mod.df <- arrange(out.mod.df, provname, year, age)
  
  df.out <- arrange(df.out, provname, year, age)
  
  min.df.prov <- ddply(df.out, .(provname), summarise, min_year = min(year))
  out.mod.df <- merge(out.mod.df, min.df.prov, by="provname", all.x = T, sort = F)
  
  # and remove modeling projections from years with no data
  out.mod.df$year_early <- 1
  out.mod.df$year_early[out.mod.df$year>=out.mod.df$min_year] <- 0
  out.mod.df <- subset(out.mod.df, year_early==0)
  mod.df <- dplyr::select(out.mod.df, year, age, cum_prop_cases, provname)
  df.merge <- dplyr::select(df.out, year, age, cum_prop_cases, provname)
  names(mod.df) <- c("year", "age", "cum_prop_cases", "provname")
  mod.df$type <- "model"
  names(df.merge) <- c("year", "age", "cum_prop_cases", "provname")
  df.merge$type <- "data"
  
  #and combine
  all.mod.dat <- rbind(df.merge, mod.df)
  return(all.mod.dat)
  
}
