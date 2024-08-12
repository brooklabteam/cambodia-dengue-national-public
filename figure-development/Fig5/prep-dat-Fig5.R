rm(list = ls())


library(ggplot2)
#library(bobfunctions2)
library(plyr)
library(dplyr)
library(stringr)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national-public"
setwd(paste0(homewd, "/figure-development/Fig5/"))

#load the output from the previous trials
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-no-wane.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-no-wane-lci.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-no-wane-uci.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-rep-2007.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-rep-2007-lci.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-rep-2007-uci.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-rep-2019.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-rep-2019-lci.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/sim-new/cam-sim-geno-rep-2019-uci.Rdata"))


summarise.age.dist.wane <- function(dat, year.start){
  
  dat1 = subset(dat,year>=year.start)
  
  #summarise into primary, secondary, tertiary, quaternary cases
  
  dat1$case_type <- (str_extract(dat1$class, "[aA-zZ]+"))
  dat1$case_sum <- nchar(dat1$class)
  
  dat1$serotype <- NA  
  dat1$serotype[dat1$case_type=="I"] <- str_sub(dat1$class[dat1$case_type=="I"], start=-1)
  
  dat1$state <- NA  
  dat1$state[dat1$case_type=="S"] <- "S"
  dat1$state[dat1$case_type=="P" | dat1$case_type=="Pms"] <- "Temporary-Heterotypic-Immunity"
  dat1$state[dat1$case_type=="Pm" ] <- "Pm"
  
  dat1$state[dat1$case_type=="I" & dat1$case_sum==2] <- "Primary-Infection"
  
  # you want to track how often you get "reinfections" within a serotype - 
  # these have to be serotype 4. and then, they need to have previously also experienced
  # infection with serotype 1
  
  dat1$state[dat1$class=="I14"] <- "Primary-Re-Infection"
  dat1$state[dat1$serotype=="4" & dat1$case_sum==4 & dat1$case_type=="I" & dat1$class!="I234" & dat1$class!="I324"] <- "Secondary-Re-Infection"
  dat1$state[dat1$serotype=="4" & dat1$case_sum==5 & dat1$case_type=="I"] <- "Tertiary-Re-Infection"
  
  #and the rest of the cases
  dat1$state[dat1$case_type=="I" & dat1$case_sum==3 & is.na(dat1$state)] <- "Secondary-Infection"
  dat1$state[dat1$class=="I41"] <- "Not-Possible"
  dat1$state[dat1$case_type=="I" & dat1$case_sum==4 & is.na(dat1$state)] <- "Tertiary-Infection"
  dat1$state[dat1$case_type=="I" & dat1$case_sum==5 & is.na(dat1$state) & dat1$serotype!=1] <- "Tertiary-Infection-After-Reinfection"
  
  #and there are a few that still need editing - can't have 1 after 4
  dat1$state[dat1$case_type=="I" & dat1$case_sum==5 & dat1$serotype==1] <- "Not-Possible"
  dat1$state[dat1$class=="I241" | dat1$class=="I341"| dat1$class=="I412"|dat1$class=="I413"|dat1$class=="I421" | dat1$class=="I431"] <- "Not-Possible"
  dat1$state[dat1$class=="I142" | dat1$class=="I143"] <- "Secondary-Infection-After-Reinfection"
  
  
  dat2 = subset(dat1, state!="Not-Possible")
  
  #tracking reinfections with serotype 2
  
  #and sum by year
  df.sum <- ddply(dat2,.(year, age, state), summarise, count = sum(count))
  
  
  
  df.sum$count<- round(df.sum$count,0)
  
  
  df.sum <- df.sum[complete.cases(df.sum),]
  
  #and just focus on infections
  #df.sum.1 <- df.sum
  
  #and return these
  df.sum.I = subset(df.sum, state!="Pm" & state!="Temporary-Heterotypic-Immunity" & state!="S")
  return(df.sum.I)
}

age.out.2007 = summarise.age.dist.wane(dat=out.cam.geno.rep.2007, year.start = min(out.cam.geno.rep.2007$year))
age.out.2007.lci = summarise.age.dist.wane(dat=out.cam.geno.rep.2007.lci, year.start = min(out.cam.geno.rep.2007$year))
age.out.2007.uci = summarise.age.dist.wane(dat=out.cam.geno.rep.2007.uci, year.start = min(out.cam.geno.rep.2007$year))

age.out.2019 = summarise.age.dist.wane(dat=out.cam.geno.rep.2019, year.start = min(out.cam.geno.rep.2019$year))
age.out.2019.lci = summarise.age.dist.wane(dat=out.cam.geno.rep.2019.lci, year.start = min(out.cam.geno.rep.2019$year))
age.out.2019.uci = summarise.age.dist.wane(dat=out.cam.geno.rep.2019.uci, year.start = min(out.cam.geno.rep.2019$year))

age.out.nointro = summarise.age.dist.wane(dat=out.cam.geno, year.start = min(out.cam.geno$year))
age.out.nointro.lci = summarise.age.dist.wane(dat=out.cam.geno.lci, year.start = min(out.cam.geno$year))
age.out.nointro.uci = summarise.age.dist.wane(dat=out.cam.geno.uci, year.start = min(out.cam.geno$year))


#and select what gets counted as symptomatic
select.symptom <- function(df, criteria){
  
  
  
  
  #first, select what you want.
  if(criteria=="Secondary+Reinfection"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Primary-Re-Infection" | state=="Secondary-Re-Infection" | state=="Tertiary-Re-Infection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
  }else if(criteria=="Secondary+Reinfection+ReSecondary"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Primary-Re-Infection" | state=="Secondary-Re-Infection" | state=="Tertiary-Re-Infection" | state == "Secondary-Infection-After-Reinfection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
  }else if(criteria == "Secondary-Only"){
    df1 = subset(df, state =="Secondary-Infection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
  }else if (criteria=="Secondary+Secondary-Reinfection"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Secondary-Re-Infection")
    df1$state[df1$state!="Secondary-Infection"] <- "Repeat-Infection"
    
  }else if (criteria=="Secondary-Extension"){
    df1 = subset(df, state =="Secondary-Infection" | state=="Secondary-Re-Infection" | state=="Secondary-Infection-After-Reinfection")
    df1$state[df1$state!="Secondary-Infection" & df1$state!="Secondary-Infection-After-Reinfection"] <- "Repeat-Infection"
    
  }else if (criteria=="Increasing-Tertiary"){
    df1.sec = subset(df, state =="Secondary-Infection" )
    df1.tert = subset(df, state=="Tertiary-Infection" )
    df.perc <- cbind.data.frame(year=unique(df1.tert$year), perc_obs=seq(0,.15, length.out = length(unique(df1.tert$year))))
    df1.tert <- merge(df1.tert, df.perc, by="year", all.x=T)
    df1.tert$count <- df1.tert$count*df1.tert$perc_obs
    df1.tert <- dplyr::select(df1.tert, -(perc_obs))
    
    df1 <- rbind(df1.sec, df1.tert)
    
    
    
  }
  
  
  
  df1$case_type = "symptomatic"
  
  return(df1)
  
  
}

age.out.2007 = subset(age.out.2007, year <2021)
age.sub.2007 = select.symptom(df=age.out.2007,criteria = "Secondary-Extension")
age.sub.2007$hyp = "H4: Genotype Replacement\n+ Waning Immunity (2007)"


age.out.2007.lci = subset(age.out.2007.lci, year <2021)
age.sub.2007.lci = select.symptom(df=age.out.2007.lci,criteria = "Secondary-Extension")
age.sub.2007.lci$hyp = "H4: Genotype Replacement\n+ Waning Immunity (2007)"


age.out.2007.uci = subset(age.out.2007.uci, year <2021)
age.sub.2007.uci = select.symptom(df=age.out.2007.uci,criteria = "Secondary-Extension")
age.sub.2007.uci$hyp = "H4: Genotype Replacement\n+ Waning Immunity (2007)"


age.out.2019 = subset(age.out.2019, year <2021)
age.sub.2019 = select.symptom(df=age.out.2019,criteria = "Secondary-Extension")
age.sub.2019$hyp = "H4: Genotype Replacement\n+ Waning Immunity (2019)"

age.out.2019.lci = subset(age.out.2019.lci, year <2021)
age.sub.2019.lci = select.symptom(df=age.out.2019.lci,criteria = "Secondary-Extension")
age.sub.2019.lci$hyp = "H4: Genotype Replacement\n+ Waning Immunity (2019)"

age.out.2019.uci = subset(age.out.2019.uci, year <2021)
age.sub.2019.uci = select.symptom(df=age.out.2019.uci,criteria = "Secondary-Extension")
age.sub.2019.uci$hyp = "H4: Genotype Replacement\n+ Waning Immunity (2019)"


age.sub.2019.b = select.symptom(df=age.out.2019,criteria = "Secondary-Only")
age.sub.2019.b$hyp = "H2: Genotype Intro\n+ Normal Immunity (2019)"

age.sub.2019.lci.b = select.symptom(df=age.out.2019.lci,criteria = "Secondary-Only")
age.sub.2019.lci.b$hyp = "H2: Genotype Intro\n+ Normal Immunity (2019)"

age.sub.2019.uci.b = select.symptom(df=age.out.2019.uci,criteria = "Secondary-Only")
age.sub.2019.uci.b$hyp = "H2: Genotype Intro\n+ Normal Immunity (2019)"


age.sub.2019.c = select.symptom(df=age.out.2019,criteria = "Increasing-Tertiary")
age.sub.2019.c$hyp = "H3: Genotype Intro + Increasing\nTertiary Case Detection"

age.sub.2019.lci.c = select.symptom(df=age.out.2019.lci,criteria = "Increasing-Tertiary")
age.sub.2019.lci.c$hyp = "H3: Genotype Intro + Increasing\nTertiary Case Detection"

age.sub.2019.uci.c = select.symptom(df=age.out.2019.uci,criteria = "Increasing-Tertiary")
age.sub.2019.uci.c$hyp = "H3: Genotype Intro + Increasing\nTertiary Case Detection"


age.out.nointro = subset(age.out.nointro, year<2021)
age.sub.tert = select.symptom(df=age.out.nointro,criteria = "Increasing-Tertiary")
age.sub.tert$hyp = "H1: Increasing Tertiary\nCase Detection"

age.out.nointro.lci = subset(age.out.nointro.lci, year<2021)
age.sub.tert.lci = select.symptom(df=age.out.nointro.lci,criteria = "Increasing-Tertiary")
age.sub.tert.lci$hyp = "H1: Increasing Tertiary\nCase Detection"

age.out.nointro.uci = subset(age.out.nointro.uci, year<2021)
age.sub.tert.uci = select.symptom(df=age.out.nointro.uci,criteria = "Increasing-Tertiary")
age.sub.tert.uci$hyp = "H1: Increasing Tertiary\nCase Detection"



age.sub.H0 = select.symptom(df=age.out.nointro,criteria = "Secondary-Only")
age.sub.H0$hyp = "H0: Normal Demographic\nSimulation"



age.sub.H0.lci = select.symptom(df=age.out.nointro.lci,criteria = "Secondary-Only")
age.sub.H0.lci$hyp = "H0: Normal Demographic\nSimulation"



age.sub.H0.uci = select.symptom(df=age.out.nointro.uci,criteria = "Secondary-Only")
age.sub.H0.uci$hyp = "H0: Normal Demographic\nSimulation"


#put all the data together
comp.dat <- rbind(age.sub.H0, age.sub.tert, age.sub.2019.b, age.sub.2019.c, age.sub.2019, age.sub.2007)
#comp.dat$hyp <- factor(comp.dat$hyp, levels = c("H0: Normal Demographic\nSimulation", "H1: Increasing Tertiary\nCase Detection", "H2: Genotype Replacement\n+ Waning Immunity (2019)", "H2: Genotype Replacement\n+ Waning Immunity (2007)"))
comp.dat$hyp <- factor(comp.dat$hyp, levels = c("H0: Normal Demographic\nSimulation", 
                                                "H1: Increasing Tertiary\nCase Detection",
                                                "H2: Genotype Intro\n+ Normal Immunity (2019)",
                                                "H3: Genotype Intro + Increasing\nTertiary Case Detection",
                                                "H4: Genotype Replacement\n+ Waning Immunity (2019)", 
                                                "H4: Genotype Replacement\n+ Waning Immunity (2007)"))

#and save for fitting
save(comp.dat, file = "comp-dat-sim.Rdata") 

comp.dat.lci <- rbind(age.sub.H0.lci, age.sub.tert.lci, age.sub.2019.lci.b, age.sub.2019.lci.c, age.sub.2019.lci, age.sub.2007.lci)
comp.dat.lci$hyp <- factor(comp.dat.lci$hyp, levels = c("H0: Normal Demographic\nSimulation", 
                                                        "H1: Increasing Tertiary\nCase Detection",
                                                        "H2: Genotype Intro\n+ Normal Immunity (2019)",
                                                        "H3: Genotype Intro + Increasing\nTertiary Case Detection",
                                                        "H4: Genotype Replacement\n+ Waning Immunity (2019)", 
                                                        "H4: Genotype Replacement\n+ Waning Immunity (2007)"))

save(comp.dat.lci, file = "comp-dat-sim-lci.Rdata") 
comp.dat.uci <- rbind(age.sub.H0.uci, age.sub.tert.uci, age.sub.2019.uci.b, age.sub.2019.uci.c, age.sub.2019.uci,  age.sub.2007.uci)
comp.dat.uci$hyp <- factor(comp.dat.uci$hyp, levels = c("H0: Normal Demographic\nSimulation", 
                                                        "H1: Increasing Tertiary\nCase Detection",
                                                        "H2: Genotype Intro\n+ Normal Immunity (2019)",
                                                        "H3: Genotype Intro + Increasing\nTertiary Case Detection",
                                                        "H4: Genotype Replacement\n+ Waning Immunity (2019)", 
                                                        "H4: Genotype Replacement\n+ Waning Immunity (2007)"))

save(comp.dat.uci, file = "comp-dat-sim-uci.Rdata") 






