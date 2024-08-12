

rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(lubridate)
library(matrixcalc)


# take the province-specific lambdas that have converged and 
# estimate an age specific multiplier for these data
# here, this is lambda across all provinces

#homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
#setwd(paste0(homewd, "/figure-development/Fig3/age-fit-alt-two-short-trim/"))

#compare N serotype hypothesis across all years to cumulative case data
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
model.age.incidence.series.age.mult <- function(par.dat, age_vect, year.start, age_mult){
  
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
  
  #year.start= min(par.dat$year)
  
  
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
get.min.year <- function(df){
  out = min(df$year)
  return(out)
}
log.lik.fit.all.prov.age <- function(par, par.dat, dat){
  
  #par.dat$lambda does not change here
  age.mult <- par
  
  age_vect=seq(0, max(dat$age), by=1)
  
  #now, split the par.dat by provinces and apply simultaneously
  par.dat.split <- dlply(par.dat, .(provname))
  
  #make year start split
  year.start.split <- lapply(par.dat.split, get.min.year)
  
  out.mod.list <- mapply(model.age.incidence.series.age.mult, par.dat= par.dat.split, year.start=year.start.split, MoreArgs = list(age_mult=age.mult, age_vect=age_vect), SIMPLIFY = FALSE)
  
  out.mod <- data.table::rbindlist(out.mod.list)
  #ggplot(data=out.mod) + geom_point(aes(x=age,y= cum_prop_cases)) + facet_wrap(provname~year)
  
  # now, select only the years for fitting for which there are data
  # (the model can produce data back further in time based on the first year FOI from the dataset)
  out.mod = subset(out.mod, year >=min(dat$year))
  
  
  # merge model with data
  out.mod <- arrange(out.mod, provname, year, age)
  dat <- arrange(dat, provname, year, age)
  dat.merge <- dplyr::select(dat, provname, age, year, Nage, cum_cases, n)
  out.merge <- merge(out.mod, dat.merge, by= c("provname", "year", "age"))
  out.merge$cum_prop_cases_data <- out.merge$cum_cases/out.merge$n
  
  out.merge <- arrange(out.merge, provname, year, age)
  
  # #plot model with data
  # ggplot(data=subset(out.merge, provname=="Kampot")) +  facet_wrap(~year) +  geom_point(aes(x=age, y=cum_prop_cases_data)) + geom_line(aes(x=age, y=cum_prop_cases_data)) +geom_line(aes(x=age,y= cum_prop_cases), color="tomato") 
  # #   
  
  #how likely are the data, given the model as truth?
  ll=0
  for (i in 1:length(out.merge$age)){
    ll=ll+dbinom(out.merge$cum_cases[i],out.merge$n[i],prob=out.merge$cum_prop_cases[i],log=T)
    #print(ll)
  }
  
  # 
  
  
  
  if(ll==-Inf){
    ll <- -1000000
  }
  
  return(-ll)
}
fit.age.mult.LBFGSB.all <- function(dat,par.dat, fit.CI){
  
  
  #get dist back here for the oldest individual in the first year of the dataset
  #dist.back =   max(dat$age[dat$year==min(dat$year)]) 
  
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  
  
  #first, prep the data - split by province and year
  year.dat <- dlply(dat, .(provname, year))
  
  
  year.dat.sum <- lapply(year.dat, sum.yr.all.prov)
  
  df.out <- data.table::rbindlist( year.dat.sum)
  
  
  #now take a guess at 10 age-specific multiplicative factors for 3 different eras:
  age.guess <- rep(1,19)
  
  #lambda.guess <- par.dat$lambda
  
  
  out.NS <- optim(par = age.guess, 
                  fn=log.lik.fit.all.prov.age, 
                  method = "L-BFGS-B",
                  lower = 0.0000001, # the range for the age guess
                  upper=10, # the range for the age guess
                  par.dat=par.dat,
                  dat=df.out,
                  #control = list(maxit=1000),
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
  
  #lambda is already fixed
  
  #return age modifier as it's own data table
  
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
  
  return(age.mult.df)
  
}


#load the data, split by province and fit across the range

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

age.fit <- fit.age.mult.LBFGSB.all(dat=dat,
                          par.dat=fit.dat,
                          fit.CI = TRUE)


save(age.fit, file = "fit-many-age-mult.Rdata")
