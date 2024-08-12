rm(list=ls())


library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(epitools)
library(reshape2)

#make national population and birth data

homewd= "/Users/carabrook/Developer/cambodia-dengue-national-public/"
setwd(homewd)

#load data 
#and get national pop data
pop.dat <- read.csv(file=paste0(homewd, "data/world_bank_cambodia.csv"), header = T, stringsAsFactors = F)
head(pop.dat)

#get population vector
pop.vec <- pop.dat[2,5:ncol(pop.dat)]
names(pop.vec) <- seq(1960,2020,1) #assume this is census at the end of each year,
#so one timestep before gives you the population at the beginning of the year
pop.vec <- c(unlist(pop.vec[which(names(pop.vec)=="1960"):which(names(pop.vec)=="2020")]))
#pop.vec[length(pop.vec)] <- pop.vec[length(pop.vec)-1]
#plot(x=2001:2020, y=pop.vec) pop is increasing

#do the same for births - these are births per 1000 people
#get total births
birth.vec <- pop.dat[1,5:ncol(pop.dat)]
names(birth.vec) <- seq(1960,2020,1) #assume this is census at the end of each year,
birth.vec <- birth.vec[which(names(birth.vec)=="1960"):which(names(birth.vec)=="2020")]
birth.vec['2020'] <- birth.vec['2019'] #assume this is the same as prior year

#plot(x=2002:2020, y=birth.vec[1,]) birth rate is declining

#now scale up by population size to get total births per year
#births per 1000 ppl. births. pop/1000 * births gives you 
#total population birth rate. so this is big because the 
#population goes up even as the births go down!

#save old births just because
births.per.1000 <- birth.vec

#and get the annual death rate
death.vec <- pop.dat[3,5:ncol(pop.dat)]
names(death.vec) <- seq(1960,2020,1) #assume this is census at the end of each year,
death.vec <- death.vec[which(names(death.vec)=="1960"):which(names(death.vec)=="2020")]
death.vec['2020'] <- death.vec['2019'] #assume this is the same as prior year

#and link
#pop.vec = pop.vec[2:length(pop.vec)]
dat.all <- cbind.data.frame(year=rep(names(pop.vec),3),value=c(unlist(c(pop.vec,birth.vec, death.vec))), metric = rep(c("total\npopulation", "births per\n1000 ppl", "deaths per\n1000 ppl"), each=length(death.vec)))
head(dat.all)
dat.all$year <- as.numeric(dat.all$year)

dat.all$source <- "pre-dataset"
dat.all$source[dat.all$year>=2002] <- "with-dataset"


#dat.all = subset(dat.all, year>=2002)
#and plot
FigS14A <- ggplot(data=dat.all) + 
  facet_grid(metric~source, scales = "free", switch = "y") +
  geom_line(aes(x=year, y=value, color=metric), show.legend = F, size=1) + 
  theme_bw() + scale_color_manual(values=c( "darkcyan", "darkorchid3", "navy")) +
  scale_x_continuous(breaks=c(1960, 1980, 2000,  2010, 2020)) +
  theme(panel.grid = element_blank(), 
        panel.spacing.x = unit(c(0), "cm"),
        strip.background = element_rect(fill="white"),
        axis.title = element_blank(), 
        plot.margin = unit(c(.2,.3,.2,.2), "cm"),
        axis.text = element_text(size=14),
        strip.text = element_text(size=16),
        strip.placement = "outside")

#and pair with the duration of the multi-annual cycles

get.mean.period <- function(col.df){
  
  #bind
  dat.period <- cbind.data.frame(period=col.df$period_length, power=col.df$value)
  
  mean_period = weighted.mean(x=col.df$period_length, w=col.df$value)
  
  return(mean_period)
}
get.national.period <- function(nat.dat, popdat){
  
  
  popdat <- dplyr::select(popdat, year, pop,provname)
  popnat <- ddply(popdat, .(year), summarise, pop=sum(pop))
  head(popnat)
  
  nat.dat$year <- trunc(nat.dat$time)
  nat.dat <- merge(nat.dat, popnat, by=c("year"), all.x = T)
  head(nat.dat)
  nat.dat$cases_per_1000 <- (nat.dat$cases/nat.dat$pop)*1000
  
  nat.dat.multi <-  analyze.wavelet(nat.dat,
                                    my.series = "cases_per_1000", #cases per 1000
                                    #loess.span = 0,
                                    dt = 1/26,#this allows for annual timestep
                                    dj = 1/100, #default =1/20. image gets clearer as number on the bottom gets bigger
                                    lowerPeriod = 2,#shortest possible period (2)
                                    upperPeriod = 20, #largest possible period (in weeks; here, 10 years)
                                    make.pval = TRUE, n.sim = 100)
  
  # wt.image(anal.dat.multi, n.levels = 250, legend.params = list(lab = "cross-wavelet power levels"), spec.period.axis = list(at = c(2:10,12,14,16), labels = c(2:10,12,14,16)),
  #                    periodlab= "period (in years)", spec.time.axis = list(at = seq(1,26*18, 26), labels = 2002:2019))
  #          #
  
  multi.dat = WaveletComp::reconstruct(nat.dat.multi, only.sig=T)
  
  nat.dat$reconstructed_multi_period <- multi.dat$series$cases_per_1000.r
  
  power.df = as.data.frame(nat.dat.multi$Power)
  names(power.df) <- 1:length(nat.dat$reconstructed_multi_period)
  power.df$period <- as.numeric(rownames(power.df))
  
  power.df.melt <- melt(power.df, id.vars = "period")
  names(power.df.melt)[names(power.df.melt)=="variable"] <- "time"
  power.df.melt$time <- as.numeric(as.character(power.df.melt$time))
  
  period.df = cbind.data.frame(period=1:length(unique(power.df.melt$period)), period_length=nat.dat.multi$Period)
  
  power.df.melt <- merge(power.df.melt, period.df, by="period", all.x=T)
  
  col.split <- dlply(power.df.melt, .(time))
  
  #ggplot(power.df.melt) + geom_tile(aes(x=time, y=period, fill=value))
  
  mean.period.list <- lapply(col.split, get.mean.period)
  #max.period.list <- lapply(col.split, get.max.period)
  
  nat.dat$multi_period_mean <-  c(unlist(mean.period.list))
  #nat.dat$multi_period_max <-  c(unlist(max.period.list))
  
  return(nat.dat)
}

wave.dat <- read.csv(file = paste0(homewd, "/data/synchrony_case_data_oct15.csv"), header=T, stringsAsFactors = F)
head(wave.dat)

nat.dat <- ddply(wave.dat, .(time), summarise, cases=sum(cases))
head(nat.dat)

#join with pop data
popdat <- read.csv(file = paste0(homewd, "/data/cambodia_pop_dat.csv"), header=T, stringsAsFactors = F)
head(popdat)

nat.dat <- get.national.period(nat.dat = nat.dat,popdat = popdat)


colz = scales::hue_pal()(length(unique((wave.dat$provname)))) #25
name.list <- sort(unique(as.character(wave.dat$provname))) #alphabetical
name.list <- c(name.list[name.list!="Ratanak Kiri" & name.list!="Mondul Kiri" & name.list!= "Tboung Khmum"], "Mondul Kiri", "Ratanak Kiri", "Tboung Khmum") #these three get moved to the end
names(colz) <- name.list


FigS14B <- ggplot(wave.dat) + theme_bw() + ylab("mean period duration for multi-annual cycles\nin dengue incidence per 1000 ppl (years)") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=16),
        axis.text = element_text(size=14), legend.title = element_blank(),
        plot.margin = unit(c(.2,.2,.2,.2), "cm")) +
  scale_color_manual(values=colz) +
  geom_line(aes(x=time,y=multi_period_mean, color=provname)) +
  geom_line(data=nat.dat, aes(x=time, y= multi_period_mean), linewidth=1.5) +
  guides(color=guide_legend(ncol=1))

FigS14B


FigS14 <- cowplot::plot_grid(FigS14A, FigS14B, ncol=2, nrow = 1, labels = c("A", "B"), label_size = 22, rel_widths = c(1,1.1))



ggsave(file = paste0(homewd, "/final-figures/FigS14.png"),
       plot= FigS14,
       units="mm",  
       width=120, 
       height=55, 
       scale=3, 
       dpi=300)



library(lme4)
library(lmerTest)
wave.dat$provname <- as.factor(wave.dat$provname)

m1 <- lmer(multi_period_mean~time:provname + (1|provname), data=wave.dat)
summary(m1)  
# most provinces have a significantly decreasing period length
# for a few, it is increasing: Kampot, Kampong Speu, 
# Kandal, Kep, Prey Veng, Pursat, Svay Rieng, Takeo
m2 <- lmer(multi_period_mean~year:provname + (1|provname), data=wave.dat)
summary(m2)   
#same patterns here. same ones increasing vs. decreasing
#no sig direction for Kampong Cham and Kampong Chhnang


