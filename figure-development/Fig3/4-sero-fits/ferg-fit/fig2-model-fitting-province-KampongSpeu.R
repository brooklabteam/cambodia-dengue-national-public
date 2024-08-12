

rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(lubridate)
library(matrixcalc)


#homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
#setwd(paste0(homewd, "/figure-development/Fig3/ferg-fit/"))

#compare N serotype hypothesis across all years to cumulative case data
#functions
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
sum.yr <- function(df, age_vect){
  
  df.sum <- ddply(df, .(age), summarise, Nage = length(age))
  
  
  df.out = cbind.data.frame(age=1:max(age_vect))
  df.out <- merge(df.out, df.sum, by="age", all.x = T, sort = F)
  df.out[is.na(df.out)]<- 0
  df.out <- rbind(c(0,0), df.out)
  df.out <- arrange(df.out, age)
  df.out$cum_cases <- cumsum(df.out$Nage)
  df.out$n <- sum(df.out$Nage)
  df.out$cum_prop_cases <- cumsum(df.out$Nage)/sum(df.out$Nage)
  df.out$year = unique(df$year)
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
      
      
      
      # now, we integrate the hazard of becoming infected over all the years in each
      
      # First, integrand A deals with avoiding infection with all strains until now
      # Then, integrand B deals with getting exposed to that select strain.
      
      # first, integrand A  (avoiding infection with all other strains) 
      inte_A = sum(dur*lambda*N_sero)
      
      # them integrand B (getting exposed to a single strain)
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
  
  p.out$N_sero <- N_sero
  
  p.out$all_prim <- p.out$one_prim*p.out$N_sero
  
  
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(cbind(p.out["exposed"], p.out["naive"])) #should be 1
  p.out$sum_naive_prim_multi <- rowSums(cbind(p.out["naive"], p.out["all_prim"], p.out["multi"])) #should be 1
  
  
  p.out <- arrange(p.out, year, age)
  
  # ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  + 
  # geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  p.out$cum_prop_cases <- p.out$multi
  p.out <- p.out[complete.cases(p.out),]
  
  
  p.add <- cbind.data.frame(year = unique(p.out$year), age = rep(0, length(unique(p.out$year))), exposed=0, naive=1, one_prim=0, multi=0, N_sero = N_sero[1], all_prim=0, sum_exp_naive=1, sum_naive_prim_multi=1, cum_prop_cases =0)
  
  
  p.sum <- rbind(p.add, p.out)
  p.sum <- arrange(p.sum, year, age)
  
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
  
  out.mod <- model.age.incidence.series(par.dat = par.dat, year.start,
                                        age_vect=age_vect)  
  
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
  # ggplot(data=out.merge) + geom_line(aes(x=age,y= cum_prop_cases), color="tomato") + facet_wrap(~year) +  geom_point(aes(x=age, y=cum_prop_cases_data)) + geom_line(aes(x=age, y=cum_prop_cases_data))
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
fit.all.yrs.seq.yr.LBFGSB <- function(dat,  lambda.guess, N.sero.fix,  fit.CI){
  
  
  #get dist back here for the oldest individual in the first year of the dataset
  dist.back =   max(dat$age[dat$year==min(dat$year)]) 
  
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  
  
  #first, prep the data
  year.dat <- dlply(dat, .(year))
  
  
  year.dat.sum <- lapply(year.dat, sum.yr.all)
  df.out <- data.table::rbindlist( year.dat.sum)
  
  
  
  # #head(df.out)
  #      ggplot(data=df.out) + geom_point(aes(x=age, y=cum_prop_cases)) +
  #         geom_line(aes(x=age, y=cum_prop_cases)) + facet_wrap(~year)
  # # # # # # # 
  # #make your guess parameters
  #lambda is takes data from the previous year and creates infections in this year
  
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
                  dat=df.out, hessian = T)
  
  
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
  
  return(par.dat)
  
}
get.lliks <- function(index.par, par.dat, df,  n.iterations, age_vect, year.start){
  
  #sample from 0 to 300, using 3000 numbers
  lambda.choose = par.dat$lambda[index.par]
  lambda.test = seq(par.dat$lambda_min[index.par], par.dat$lambda_max[index.par], length=n.iterations)
  
  #and
  store.llik <- list()
  for (i in 1:length(lambda.test)){
    
    #replace the parameter of interest with your sampled value and return the likelihood
    lambda.list <- par.dat$lambda
    lambda.list[index.par] <- lambda.test[i]
    
    out.lik <- log.lik.fit.all(par = log( lambda.list),
                               par.dat=par.dat, 
                               age_vect=age_vect,
                               year.start = year.start,
                               dat=df)
    
    store.llik[[i]] <- out.lik
    
    
  }
  
  #and bind the value and the likelihood
  llik.out <- cbind.data.frame(par = lambda.test,neg_llik=c(unlist(store.llik)))
  
  #with(llik.out, plot(par, neg_llik, type="l"))
  
  llik.out <- subset(llik.out, neg_llik<Inf)
  tmp=smooth.spline(llik.out$par, llik.out$neg_llik) 
  new=seq(par.dat$lambda_min[index.par],par.dat$lambda_max[index.par], by=.001) 
  interp = predict(tmp, new)$y
  #and then use the fact that the profile likelihood is χ2-distributed to erect 95% confidence intervals:
  mle1=new[which.min(interp)]
  tmp3=(predict(tmp, new)$y-min(predict(tmp, new)$y))-qchisq(0.95,1)
  CI <- range(new[tmp3<0])
  
  llik.out$par[llik.out$neg_llik==min(llik.out$neg_llik)]
  #add the info about what par
  llik.out$par_test_index = index.par
  llik.out$year = par.dat$year[index.par]
  
  #and return the CIs instead
  
  ret.llik <- par.dat[index.par,]
  ret.llik$lci <- CI[1]
  ret.llik$uci <- CI[2]
  
  #and return
  return( ret.llik )
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
dat = subset(dat, provname=="Kampong Speu")

out <- fit.all.yrs.seq.yr.LBFGSB(dat = dat,
                               lambda.guess=0.1,
                               N.sero.fix= 4,
                               fit.CI = TRUE)
out$provname = unique(dat$provname)

namehold = unique(dat$provname)
namehold = sub(" ", "-", namehold, fixed = T)

save(out, file = paste0("fit-prov-", namehold, ".Rdata"))
