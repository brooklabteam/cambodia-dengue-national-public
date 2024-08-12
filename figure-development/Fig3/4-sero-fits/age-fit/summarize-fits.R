rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig3/4-sero-fits/age-fit/"))

#load and combine the fits and plot them together (except not the ones with Inf)
load("fit-many-age-mult.Rdata")
age.fit #not converged -- need to rerun
load("refit-many-age-mult.Rdata")
age.refit #converged
#age.mult = out

#and reshape data for plotting with age min and age max as a single "age" column

age.long <- melt(age.refit, id.vars = c("year_range", "age_mult", "lci_mult", "uci_mult"), measure.vars = c("age_min", "age_max"))

head(age.long)
names(age.long)[names(age.long)=="value"] <- "age"

p1 <- ggplot(age.long) + theme_bw() +
  geom_line(aes(x=age, y=age_mult), size=1) + facet_grid(~year_range) +
  geom_ribbon(aes(x=age, ymin=lci_mult, ymax=uci_mult), alpha=.3) +
  theme(panel.grid = element_blank(), axis.title = element_text(size=13), axis.text = element_text(size=10)) + ylab("age multiplier") +
  geom_hline(aes(yintercept=1), linetype=2)
print(p1)

#save this plot

ggsave(filename = "age-multiplier.png",
       plot = p1,
       units="mm",  
       width=55, 
       height=40, 
       scale=3, 
       dpi=300)


## and including waning 

load("fit-many-age-mult-wane.Rdata")
#age.mult = out

#and reshape data for plotting with age min and age max as a single "age" column

age.long.wane <- melt(age.fit.wane, id.vars = c("year_range", "age_mult", "lci_mult", "uci_mult"), measure.vars = c("age_min", "age_max"))

head(age.long.wane)
names(age.long.wane)[names(age.long.wane)=="value"] <- "age"

p2 <- ggplot(age.long.wane) + theme_bw() +
  geom_line(aes(x=age, y=age_mult), size=1) + facet_grid(~year_range) +
  #geom_ribbon(aes(x=age, ymin=lci_mult, ymax=uci_mult), alpha=.3) +
  theme(panel.grid = element_blank(), axis.title = element_text(size=13), axis.text = element_text(size=10)) + ylab("age multiplier") 
print(p2)

#and all together
age.long.wane$type <- "wane"
age.long$type <- "foi"

age.long <- rbind(age.long, age.long.wane)


p3 <- ggplot(age.long) + theme_bw() + 
  geom_line(aes(x=age, y=age_mult, color=type), size=1) + facet_grid(~year_range) +
  #geom_ribbon(aes(x=age, ymin=lci_mult, ymax=uci_mult), alpha=.3) +
  theme(panel.grid = element_blank(), axis.title = element_text(size=13), axis.text = element_text(size=10)) + ylab("age multiplier") 
print(p3)
