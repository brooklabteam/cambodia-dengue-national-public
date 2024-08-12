rm(list = ls())


library(ggplot2)
#library(bobfunctions2)
library(plyr)
library(dplyr)
library(stringr)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national-public"
setwd(paste0(homewd, "/figure-development/Fig5/"))

#load the output from the previous trials
load(paste0(homewd,"/figure-development/Fig5/comp-dat-sim.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/comp-dat-sim-lci.Rdata"))
load(paste0(homewd,"/figure-development/Fig5/comp-dat-sim-uci.Rdata"))



#then, feed into plotting

#helper functions
mean.age <- function(df){
  df$mult <- df$age*df$count
  mean.age <- sum(df$mult)/sum(df$count)
  df2 <- cbind.data.frame(year=unique(df$year), mean_age=mean.age, hyp = unique(df$hyp))
  return(df2)
}
replicate.data.type <- function(df, slim.quant){
  #print(unique(df$year))
  if(df$count>0){
    new.dat = cbind.data.frame(age=rep(df$age,(df$count)), case=rep(1, (df$count)))  
    new.dat$year <- unique(df$year)
    new.dat$state <- unique(df$state)
    new.dat$hyp <- unique(df$hyp)
    
    #then, cut to 5% of the cases:
    #all should be the same, so just take the top 5% of rows
    n.row= round(nrow(new.dat)*(slim.quant),0)
    
    new.dat <- new.dat[1:n.row,]
    
    
    return(new.dat)
  }
  
  
}
cum.sum.year <- function(df){
  df.sum <- ddply(df,.(age), summarise, cases=sum(count))
  
  df.sum$cum_cases = cumsum(df.sum$cases)
  df.sum$cum_prop_cases <- df.sum$cum_cases/sum(df.sum$cases)
  df.sum$year <- unique(df$year)
  df.sum$hyp <- unique(df$hyp)
  return(df.sum)
}


#now, choose those that will go in the supplement vs. the main text
#supplement first
comp.dat$hyp <- as.character(comp.dat$hyp)
comp.dat.lci$hyp <- as.character(comp.dat.lci$hyp)
comp.dat.uci$hyp <- as.character(comp.dat.uci$hyp)

dat.supp = subset(comp.dat, hyp == "H2: Genotype Intro\n+ Normal Immunity (2019)" | hyp == "H3: Genotype Intro + Increasing\nTertiary Case Detection")
dat.supp.lci = subset(comp.dat.lci, hyp == "H2: Genotype Intro\n+ Normal Immunity (2019)" | hyp == "H3: Genotype Intro + Increasing\nTertiary Case Detection")
dat.supp.uci = subset(comp.dat.uci, hyp == "H2: Genotype Intro\n+ Normal Immunity (2019)" | hyp == "H3: Genotype Intro + Increasing\nTertiary Case Detection")

comp.dat = subset(comp.dat, hyp != "H2: Genotype Intro\n+ Normal Immunity (2019)" & hyp != "H3: Genotype Intro + Increasing\nTertiary Case Detection")
comp.dat.lci = subset(comp.dat.lci, hyp != "H2: Genotype Intro\n+ Normal Immunity (2019)" & hyp != "H3: Genotype Intro + Increasing\nTertiary Case Detection")
comp.dat.uci = subset(comp.dat.uci, hyp != "H2: Genotype Intro\n+ Normal Immunity (2019)" & hyp != "H3: Genotype Intro + Increasing\nTertiary Case Detection")

comp.dat$hyp[comp.dat$hyp=="H4: Genotype Replacement\n+ Waning Immunity (2019)"] <- "H2: Genotype Replacement\n+ Waning Immunity (2019)"
comp.dat$hyp[comp.dat$hyp=="H4: Genotype Replacement\n+ Waning Immunity (2007)"] <- "H2: Genotype Replacement\n+ Waning Immunity (2007)"

comp.dat.lci$hyp[comp.dat.lci$hyp=="H4: Genotype Replacement\n+ Waning Immunity (2019)"] <- "H2: Genotype Replacement\n+ Waning Immunity (2019)"
comp.dat.lci$hyp[comp.dat.lci$hyp=="H4: Genotype Replacement\n+ Waning Immunity (2007)"] <- "H2: Genotype Replacement\n+ Waning Immunity (2007)"

comp.dat.uci$hyp[comp.dat.uci$hyp=="H4: Genotype Replacement\n+ Waning Immunity (2019)"] <- "H2: Genotype Replacement\n+ Waning Immunity (2019)"
comp.dat.uci$hyp[comp.dat.uci$hyp=="H4: Genotype Replacement\n+ Waning Immunity (2007)"] <- "H2: Genotype Replacement\n+ Waning Immunity (2007)"

comp.dat$hyp <- factor(comp.dat$hyp, levels = c("H0: Normal Demographic\nSimulation", "H1: Increasing Tertiary\nCase Detection", "H2: Genotype Replacement\n+ Waning Immunity (2019)", "H2: Genotype Replacement\n+ Waning Immunity (2007)"))
comp.dat.lci$hyp <- factor(comp.dat.lci$hyp, levels = c("H0: Normal Demographic\nSimulation", "H1: Increasing Tertiary\nCase Detection", "H2: Genotype Replacement\n+ Waning Immunity (2019)", "H2: Genotype Replacement\n+ Waning Immunity (2007)"))
comp.dat.uci$hyp <- factor(comp.dat.uci$hyp, levels = c("H0: Normal Demographic\nSimulation", "H1: Increasing Tertiary\nCase Detection", "H2: Genotype Replacement\n+ Waning Immunity (2019)", "H2: Genotype Replacement\n+ Waning Immunity (2007)"))
#first, make the first column

column.1 <- function(dat, dat.lci, dat.uci, year.start){
  dat1 = subset(dat, year >= year.start) 
  dat1.lci = subset(dat.lci, year >= year.start) 
  dat1.uci = subset(dat.uci, year >= year.start) 
  
  denv.case = subset(dat1, case_type == "symptomatic")
  denv.case.lci = subset(dat1.lci, case_type == "symptomatic")
  denv.case.uci = subset(dat1.uci, case_type == "symptomatic")
  
  
  #dat.ts <- ddply(denv.case, .(year, time), summarise, count=sum(count))
  dat.ts <- ddply(denv.case, .(hyp, year), summarise, count=sum(count))
  dat.ts.lci <- ddply(denv.case.lci, .(hyp,year), summarise, count=sum(count))
  dat.ts.uci <- ddply(denv.case.uci, .(hyp, year), summarise, count=sum(count))
  #dat.ts$hyp <- factor(dat.ts$hyp)
  dat.ts$plot_type <- "total reported cases\nfrom deterministic simulation"
  # dat.ts.lci$hyp <- factor(dat.ts.lci$hyp)
  dat.ts.lci$plot_type <- "total reported cases\nfrom deterministic simulation"
  # dat.ts.uci$hyp <- factor(dat.ts.uci$hyp)
  dat.ts.uci$plot_type <- "total reported cases\nfrom deterministic simulation"
  # 
  dat.ts.lci.merge <- dplyr::select(dat.ts.lci, hyp, year, count)
  dat.ts.uci.merge <- dplyr::select(dat.ts.uci, hyp, year, count)
  names(dat.ts.lci.merge)[3] <- "lci"
  names(dat.ts.uci.merge)[3] <- "uci"
  # 
  dat.ts <- merge(dat.ts, dat.ts.lci.merge, by =c("hyp", "year"))
  dat.ts <- merge(dat.ts, dat.ts.uci.merge, by =c("hyp", "year"))
  dat.ts$lci_new <- NA
  dat.ts$uci_new <- NA
  dat.ts$lci_new[dat.ts$lci<dat.ts$uci] <- dat.ts$lci[dat.ts$lci<dat.ts$uci]
  
  #dat.ts$lci_new[dat.ts$lci>dat.ts$uci] <- dat.ts$uci[dat.ts$lci>dat.ts$uci]
  #dat.ts$uci_new[dat.ts$lci>dat.ts$uci] <- dat.ts$lci[dat.ts$lci>dat.ts$uci]
  dat.ts$uci_new[dat.ts$lci<dat.ts$uci] <- dat.ts$uci[dat.ts$lci<dat.ts$uci]
  # 
  
  #dat.ts$count_dss <- dat.ts$count*perc_dss
  #dat.ts$count_mort <- dat.ts$count_dss*perc_mort #some fraction of dss cases is the total mortality
  
  p1 <- ggplot(dat.ts) + theme_bw() + facet_grid(hyp~plot_type, switch = "y") +
    #geom_vline(aes(xintercept=2012), linetype=2)+
    #geom_vline(aes(xintercept=2019), linetype=2) +
    #geom_ribbon(aes(x=year, ymin=lci_new, ymax=uci_new),alpha=.3) +
    geom_ribbon(aes(x=year, ymin=lci, ymax=uci),alpha=.3) +
    geom_line(aes(x=year, y=count), size=.3) + 
    #geom_line(aes(x=year, y=count_dss), size=.3, color="navy") + 
    #geom_line(aes(x=year, y=count_mort), size=.3, color="green") + 
    theme(panel.grid = element_blank(), 
          strip.background.y  = element_rect(fill="white"),
          strip.background.x  = element_rect(fill="white"),
          axis.title = element_blank(), 
          plot.margin = unit(c(.2,.1,.2,.1), "cm"), 
          axis.text = element_text(size=12),
          strip.placement = "outside",
          strip.text.x.top  = element_text(size=12),
          strip.text.y.left  = element_text(size=12)) +
    #      coord_cartesian(xlim=c(2015,2020), ylim=c(0,13000)) + #scale_y_log10() +
    #scale_y_continuous(breaks=c(2500,7500,12500)) +
    scale_x_continuous(breaks=c(2000,2005,2010,2015, 2020)) +
    geom_vline(aes(xintercept=2007), linetype=2, color="red") +
    geom_vline(aes(xintercept=2012), linetype=2, color="red") +
    geom_vline(aes(xintercept=2019), linetype=2, color="red")
  
  #scale_x_continuous(breaks=c(2016,2018, 2020)) 
          
  return(p1)
}
col1 <- column.1(dat =comp.dat, dat.lci = comp.dat.lci, dat.uci = comp.dat.uci, year.start = 2002)
#column 2 age distribution of cases
column.2 <- function (dat, year.start){
  
  dat1 = subset(dat, year >= year.start) 
  
  denv.case = subset(dat1, case_type == "symptomatic")
  
  
  
  
  #denv.case$round <- factor(denv.case$round, levels=c("secondary", "tertiary"))
  #subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
   if(length(unique(denv.case$age))>1){
     denv.case = subset(denv.case, age<max(denv.case$age))  
   }
  # #denv.case$year <- trunc(denv.case$time)
  
  #first, split by year
  #then, make a "case" out of each
  #dat.split <- dlply(denv.case,.(year))
  
    df.sum = ddply(denv.case,.(hyp, year,age, state),summarise, count=sum(count))  
  
  
  
  
  #and build out
  #df.sum$count<-round(df.sum$count,0)
  
  df.sum$count<- round(df.sum$count,0)
  #df.sum$count<- ceiling(df.sum$count)
  
  
  
  
  
  
  df.year.only <- dlply(df.sum,.(hyp, year))
  
  #get mean age
  mean.df <- data.table::rbindlist(lapply(df.year.only, mean.age))
  mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  
  
  #split by age and year and state
  df.age.list <- dlply(df.sum,.(hyp, year, age, state))
  
  
  
  
  
  dat.age <- data.table::rbindlist(lapply(df.age.list, replicate.data.type, slim.quant=1))
  dat.age$state[dat.age$state!="Secondary-Infection"] <- "alternative"
  dat.age$state[dat.age$state=="Secondary-Infection"] <- "secondary"
  dat.age$plot_type <- "age distribution of\nreported cases by year"
  #dat.age$hyp = unique(dat.age$hyp)
  
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  
  
  colz = c('secondary'="black", "alternative"="royalblue3")
  #and plot
  p1 <- ggplot(dat.age) + facet_grid(hyp~plot_type) + theme_bw()+
    geom_jitter(aes(x=year, y=age, color=state), width=.09, height=.09, size=.09, alpha=.5, show.legend = F) +
    scale_color_manual(values=colz) +
    geom_violin(aes(x=year,y=age, group=year),  color="gray55", draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
    geom_line(data=mean.df, aes(x=year, y=mean_age), color="tomato", size=.8) +
    theme(panel.grid = element_blank(), 
          strip.background.y  = element_blank(),
          strip.background.x  = element_rect(fill="white"),
          axis.title = element_blank(), 
          #axis.title.y = element_text(size=12), 
          axis.text = element_text(size=10),
          strip.text.x.top = element_text(size=12),
          strip.text.y =element_blank()) + 
    scale_x_continuous(breaks=c(c(2000,2005,2010,2015,2020)))
    
  
  return(p1)
  
}
col2 <- column.2(dat =comp.dat, year.start = 2002)
#column 3 cumulative proportion of cases
column.3 <- function(dat, year.start){
  
  dat1 = subset(dat, year >= year.start) 
  
  denv.case = subset(dat1, case_type == "symptomatic")
  
  
  
  
  
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  
  
  #head(denv.case)
  denv.case = subset(denv.case, year>=year.start)
  
  denv.split <- dlply(denv.case, .(hyp, year))
  
  
  denv.cum <- lapply(denv.split, cum.sum.year)
  denv.dat <- data.table::rbindlist(denv.cum)
  denv.dat$year <- as.factor(denv.dat$year)
  denv.dat$plot_type <- "cumulative proportion\nof cases by age"
  
  #head(denv.dat)
  
  p1 <- ggplot(data=denv.dat) + facet_grid(hyp~plot_type) + theme_bw()+
        geom_line(aes(x=age, y=cum_prop_cases, color=year)) +
        scale_color_viridis_d(direction=-1,option="turbo") +
        theme(panel.grid = element_blank(), 
              strip.background.y  = element_blank(),
              strip.background.x  = element_rect(fill="white"),
          legend.title = element_blank(),
          axis.title = element_blank(), 
          legend.key.size = unit(c(.4), "cm"), 
          legend.position = c(.6,.91),
          legend.text = element_text(size=10),
          #axis.title.y = element_text(size=12), 
          axis.text = element_text(size=10),
          strip.text.y = element_blank(),
          strip.text.x.top = element_text(size=12)) + 
        guides(color=guide_legend(ncol=3)) 
  
  return(p1)
  
}
col3 <- column.3(dat =comp.dat,  year.start = 2002)


#and put them together

Fig5 <- cowplot::plot_grid(col1, col2, col3, nrow = 1, ncol=3, labels=c("A", "B", "C"),rel_widths = c(1.15,1,1),  label_size = 22, align = "v", label_x = c(0,-.03,0))

ggsave(filename = paste0(homewd, "/final-figures/Fig5.png"),
       plot = Fig5,
       units="mm",  
       width=110, 
       height=85, 
       scale=3, 
       dpi=300)

ggsave(filename = paste0(homewd, "/final-figures/Fig5.pdf"),
       plot = Fig5,
       units="mm",  
       width=110, 
       height=85, 
       scale=3, 
       dpi=300)

# and the supplementary figure (S21) just includes those simulations that did not make the main text cut
# same as plot above
dat.supp$hyp[dat.supp$hyp=="H2: Genotype Intro\n+ Normal Immunity (2019)"] <- "H3: Genotype Intro (2019)\n+ Normal Immunity"
dat.supp$hyp[dat.supp$hyp=="H3: Genotype Intro + Increasing\nTertiary Case Detection"] <- "H4: Genotype Intro (2019)\n+ Increasing Tertiary Case Detection"

dat.supp.lci$hyp[dat.supp.lci$hyp=="H2: Genotype Intro\n+ Normal Immunity (2019)"] <- "H3: Genotype Intro (2019)\n+ Normal Immunity"
dat.supp.lci$hyp[dat.supp.lci$hyp=="H3: Genotype Intro + Increasing\nTertiary Case Detection"] <- "H4: Genotype Intro (2019)\n+ Increasing Tertiary Case Detection"

dat.supp.uci$hyp[dat.supp.uci$hyp=="H2: Genotype Intro\n+ Normal Immunity (2019)"] <- "H3: Genotype Intro (2019)\n+ Normal Immunity"
dat.supp.uci$hyp[dat.supp.uci$hyp=="H3: Genotype Intro + Increasing\nTertiary Case Detection"] <- "H4: Genotype Intro (2019)\n+ Increasing Tertiary Case Detection"


column.1.supp <- function(dat, dat.lci, dat.uci, year.start){
  dat1 = subset(dat, year >= year.start) 
  dat1.lci = subset(dat.lci, year >= year.start) 
  dat1.uci = subset(dat.uci, year >= year.start) 
  
  denv.case = subset(dat1, case_type == "symptomatic")
  denv.case.lci = subset(dat1.lci, case_type == "symptomatic")
  denv.case.uci = subset(dat1.uci, case_type == "symptomatic")
  
  
  #dat.ts <- ddply(denv.case, .(year, time), summarise, count=sum(count))
  dat.ts <- ddply(denv.case, .(hyp, year), summarise, count=sum(count))
  dat.ts.lci <- ddply(denv.case.lci, .(hyp,year), summarise, count=sum(count))
  dat.ts.uci <- ddply(denv.case.uci, .(hyp, year), summarise, count=sum(count))
  #dat.ts$hyp <- factor(dat.ts$hyp)
  dat.ts$plot_type <- "total reported cases\nfrom deterministic simulation"
  # dat.ts.lci$hyp <- factor(dat.ts.lci$hyp)
  dat.ts.lci$plot_type <- "total reported cases\nfrom deterministic simulation"
  # dat.ts.uci$hyp <- factor(dat.ts.uci$hyp)
  dat.ts.uci$plot_type <- "total reported cases\nfrom deterministic simulation"
  # 
  dat.ts.lci.merge <- dplyr::select(dat.ts.lci, hyp, year, count)
  dat.ts.uci.merge <- dplyr::select(dat.ts.uci, hyp, year, count)
  names(dat.ts.lci.merge)[3] <- "lci"
  names(dat.ts.uci.merge)[3] <- "uci"
  # 
  dat.ts <- merge(dat.ts, dat.ts.lci.merge, by =c("hyp", "year"))
  dat.ts <- merge(dat.ts, dat.ts.uci.merge, by =c("hyp", "year"))
  dat.ts$lci_new <- NA
  dat.ts$uci_new <- NA
  dat.ts$lci_new[dat.ts$lci<dat.ts$uci] <- dat.ts$lci[dat.ts$lci<dat.ts$uci]
  dat.ts$lci_new[dat.ts$lci>dat.ts$uci] <- dat.ts$uci[dat.ts$lci>dat.ts$uci]
  dat.ts$uci_new[dat.ts$lci>dat.ts$uci] <- dat.ts$lci[dat.ts$lci>dat.ts$uci]
  dat.ts$uci_new[dat.ts$lci<dat.ts$uci] <- dat.ts$uci[dat.ts$lci<dat.ts$uci]
  # 
  
  #dat.ts$count_dss <- dat.ts$count*perc_dss
  #dat.ts$count_mort <- dat.ts$count_dss*perc_mort #some fraction of dss cases is the total mortality
  
  p1 <- ggplot(dat.ts) + theme_bw() + facet_grid(hyp~plot_type, switch = "y") +
    #geom_vline(aes(xintercept=2012), linetype=2)+
    #geom_vline(aes(xintercept=2019), linetype=2) +
    #geom_ribbon(aes(x=year, ymin=lci_new, ymax=uci_new),alpha=.3) +
    geom_ribbon(aes(x=year, ymin=lci, ymax=uci),alpha=.3) +
    geom_line(aes(x=year, y=count), size=.3) + 
    #geom_line(aes(x=year, y=count_dss), size=.3, color="navy") + 
    #geom_line(aes(x=year, y=count_mort), size=.3, color="green") + 
    theme(panel.grid = element_blank(), 
          strip.background.y  = element_rect(fill="white"),
          strip.background.x  = element_rect(fill="white"),
          axis.title = element_blank(), 
          plot.margin = unit(c(.2,.1,.2,.1), "cm"), 
          axis.text = element_text(size=12),
          strip.placement = "outside",
          strip.text.x.top  = element_text(size=12),
          strip.text.y.left  = element_text(size=12)) +
          coord_cartesian( ylim=c(0,100000)) + #scale_y_log10() +
    #scale_y_continuous(breaks=c(2500,7500,12500)) +
    scale_x_continuous(breaks=c(2000,2005,2010,2015, 2020)) +
    geom_vline(aes(xintercept=2007), linetype=2, color="red") +
    geom_vline(aes(xintercept=2012), linetype=2, color="red") +
    geom_vline(aes(xintercept=2019), linetype=2, color="red")
  
  #scale_x_continuous(breaks=c(2016,2018, 2020)) 
  
  return(p1)
}
column.2.supp <- function (dat, year.start){
  
  dat1 = subset(dat, year >= year.start) 
  
  denv.case = subset(dat1, case_type == "symptomatic")
  
  
  
  
  #denv.case$round <- factor(denv.case$round, levels=c("secondary", "tertiary"))
  #subset(dat1,  class == "I12" | class=="I13" | class=="I21"  | class=="I23"  | class=="I31"  | class=="I32")
  
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  # #denv.case$year <- trunc(denv.case$time)
  
  #first, split by year
  #then, make a "case" out of each
  #dat.split <- dlply(denv.case,.(year))
  
  df.sum = ddply(denv.case,.(hyp, year,age, state),summarise, count=sum(count))  
  
  
  
  
  #and build out
  #df.sum$count<-round(df.sum$count,0)
  
  df.sum$count<- round(df.sum$count,0)
  #df.sum$count<- ceiling(df.sum$count)
  
  
  
  
  
  
  df.year.only <- dlply(df.sum,.(hyp, year))
  
  #get mean age
  mean.df <- data.table::rbindlist(lapply(df.year.only, mean.age))
  mean.df$mean_age[is.na(mean.df$mean_age)] <- 0
  
  
  #split by age and year and state
  df.age.list <- dlply(df.sum,.(hyp, year, age, state))
  
  
  
  
  
  dat.age <- data.table::rbindlist(lapply(df.age.list, replicate.data.type, slim.quant=1))
  dat.age$state[dat.age$state!="Secondary-Infection"] <- "alternative"
  dat.age$state[dat.age$state=="Secondary-Infection"] <- "secondary"
  dat.age$plot_type <- "age distribution of\nreported cases by year"
  #dat.age$hyp = unique(dat.age$hyp)
  
  #then, evenly slim the dataset to 10% of its current values, evenly by age and year
  
  
  colz = c('secondary'="black", "alternative"="royalblue3")
  #and plot
  p1 <- ggplot(dat.age) + facet_grid(hyp~plot_type) + theme_bw()+
    geom_jitter(aes(x=year, y=age, color=state), width=.09, height=.09, size=.09, alpha=.5, show.legend = F) +
    scale_color_manual(values=colz) + coord_cartesian(ylim=c(0,100)) +
    geom_violin(aes(x=year,y=age, group=year),  color="gray55", draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA)+
    geom_line(data=mean.df, aes(x=year, y=mean_age), color="tomato", size=.8) +
    theme(panel.grid = element_blank(), 
          strip.background.y  = element_blank(),
          strip.background.x  = element_rect(fill="white"),
          axis.title = element_blank(), 
          #axis.title.y = element_text(size=12), 
          axis.text = element_text(size=10),
          strip.text.x.top = element_text(size=12),
          strip.text.y =element_blank()) + 
    scale_x_continuous(breaks=c(c(2000,2005,2010,2015,2020)))
  
  
  return(p1)
  
}
column.3.supp <- function(dat, year.start){
  
  dat1 = subset(dat, year >= year.start) 
  
  denv.case = subset(dat1, case_type == "symptomatic")
  
  
  
  
  
  if(length(unique(denv.case$age))>1){
    denv.case = subset(denv.case, age<max(denv.case$age))  
  }
  
  
  #head(denv.case)
  denv.case = subset(denv.case, year>=year.start)
  
  denv.split <- dlply(denv.case, .(hyp, year))
  
  
  denv.cum <- lapply(denv.split, cum.sum.year)
  denv.dat <- data.table::rbindlist(denv.cum)
  denv.dat$year <- as.factor(denv.dat$year)
  denv.dat$plot_type <- "cumulative proportion\nof cases by age"
  
  #head(denv.dat)
  
  p1 <- ggplot(data=denv.dat) + facet_grid(hyp~plot_type) + theme_bw()+
    geom_line(aes(x=age, y=cum_prop_cases, color=year)) +
    scale_color_viridis_d(direction=-1,option="turbo") +
    theme(panel.grid = element_blank(), 
          strip.background.y  = element_blank(),
          strip.background.x  = element_rect(fill="white"),
          legend.title = element_blank(),
          axis.title = element_blank(), 
          legend.key.size = unit(c(.4), "cm"), 
          legend.position = c(.6,.8),
          legend.text = element_text(size=10),
          #axis.title.y = element_text(size=12), 
          axis.text = element_text(size=10),
          strip.text.y = element_blank(),
          strip.text.x.top = element_text(size=12)) + 
    guides(color=guide_legend(ncol=3)) 
  
  return(p1)
  
}

col1Supp <- column.1.supp(dat=dat.supp, dat.lci = dat.supp.lci, dat.uci = dat.supp.uci, year.start = 2002)
col2Supp <- column.2.supp(dat=dat.supp,  year.start = 2002)
col3Supp <- column.3.supp(dat=dat.supp, year.start = 2002)




FigS21 <- cowplot::plot_grid(col1Supp, col2Supp, col3Supp,  nrow = 1, ncol=3, labels=c("A", "B", "C"),rel_widths = c(1.15,1.05,1),  label_size = 22, align = "v", label_x = c(0,0,-.03))

ggsave(filename = paste0(homewd, "/final-figures/FigS21.png"),
       plot = FigS21,
       units="mm",  
       width=110, 
       height=60, 
       scale=3, 
       dpi=300)



