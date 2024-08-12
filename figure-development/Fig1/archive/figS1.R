rm(list=ls())


library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(epitools)
library(reshape2)

#make national population and birth data

homewd= "/Users/carabrook/Developer/cambodia-dengue-national/"
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

#save to data
write.csv(dat.all, file = paste0(homewd, "data/pop_data_full.csv"), row.names = F)

dat.all = subset(dat.all, year>=2002)
#and plot
p <- ggplot(data=dat.all) + facet_grid(metric~., scales = "free_y", switch = "y") +
    geom_line(aes(x=year, y=value, color=metric), show.legend = F, size=1) + 
    theme_bw() + scale_color_manual(values=c( "darkcyan", "darkorchid3", "navy")) +
    theme(panel.grid = element_blank(), 
          strip.background = element_rect(fill="white"),
          axis.title = element_blank(), 
          axis.text = element_text(size=14),
          strip.text = element_text(size=16),
          strip.placement = "outside")


#and save
ggsave(file =paste0(homewd, "final-figures/figS1.png"),
       plot=p,
       units="mm",  
       width=50, 
       height=70, 
       scale=3, 
       dpi=300)
