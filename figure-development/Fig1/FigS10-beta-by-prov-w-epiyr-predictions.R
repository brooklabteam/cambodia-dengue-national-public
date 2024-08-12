rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)

#set wd
homewd = "/Users/carabrook/Developer/cambodia-dengue-national-public"

#load the data for the betas by province

beta.dat <- read.csv(file = paste0(homewd, "/data/beta_TSIR_fit_province.csv"), header=T, stringsAsFactors = F)
head(beta.dat)
unique(beta.dat$provname)
unique(beta.dat$year)
beta.dat = subset(beta.dat, year!=2007 & year!=2012& year!=2019)
beta.dat$epiyr <- as.factor(beta.dat$epiyr)

#now, add in your epiyear predictions, from the climate regression
epi.dat <- read.csv(file = paste0(homewd, "/data/tsir_dat_beta_climate_province.csv"), header = T, stringsAsFactors = F)
head(epi.dat)
epi.dat$epiyr <- as.factor(epi.dat$epiyr)

#change name of province for plotting
epi.dat$provname[epi.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"
beta.dat$provname[beta.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"

legend.dummy <- cbind.data.frame(x=c(0,0), y=c(0,0), fit_type = c("'TSIR-fitted'~beta", "'climate-projected'~beta"), provname="Prey Veng")
legend.dummy$fit_type <- factor(legend.dummy$fit_type, levels=c("'TSIR-fitted'~beta","'climate-projected'~beta"))


FigS10 <- ggplot(data=subset(beta.dat, provname!="Mondul Kiri" & provname!="Ratanak Kiri")) + theme_bw() +
  geom_line(data=legend.dummy, aes(x=x,y=y,linetype=fit_type)) +
  scale_linetype_manual(name = "transmission type",  labels = scales::parse_format(), values = c(1,6))+
  geom_ribbon(aes(x=biweek, ymin= betalow, ymax=betahigh, fill=epiyr), alpha=.2) + 
  geom_ribbon(data = subset(epi.dat,provname!="Mondul Kiri" & provname!="Ratanak Kiri"), aes(x=biweek, ymin= betalow, ymax=betahigh, fill=epiyr), alpha=.2) + 
  geom_line(data = subset(epi.dat, provname!="Mondul Kiri" & provname!="Ratanak Kiri"), aes(x=biweek, y= beta, color=epiyr), size=.8, linetype="twodash") +
  geom_line(aes(x=biweek, y= beta, color=epiyr), size=.8) + scale_fill_manual(name="epidemic year", values=c("tomato", "cornflowerblue", "seagreen")) +
  scale_color_manual(name="epidemic year", values=c("tomato", "cornflowerblue", "seagreen")) +
  facet_wrap(provname~., scales = "free_y")+ylab(bquote(beta~', biweekly transmission')) +
  xlab("biweek of year") +
  theme(panel.grid = element_blank(), legend.position = c(.86,.06),
        strip.background = element_rect(fill="white"), legend.box = "horizontal",
        axis.title = element_text(size=18), axis.text = element_text(size=13))


ggsave(file = paste0(homewd, "/final-figures/FigS10.png"),
       plot = FigS10,
       units="mm",  
       width=100, 
       height=70, 
       scale=3, 
       dpi=300)


legend.dummy$provname<- "Battambang"
FigS10sub <- ggplot(data=subset(beta.dat, provname=="Battambang")) + theme_bw() +
  geom_line(data=legend.dummy, aes(x=x,y=y,linetype=fit_type)) +
  scale_linetype_manual(name = "transmission type",  labels = scales::parse_format(), values = c(1,6))+
  geom_ribbon(aes(x=biweek, ymin= betalow, ymax=betahigh, fill=epiyr), alpha=.2) + 
  geom_ribbon(data = subset(epi.dat,provname=="Battambang"), aes(x=biweek, ymin= betalow, ymax=betahigh, fill=epiyr), alpha=.2) + 
  geom_line(data = subset(epi.dat, provname=="Battambang"), aes(x=biweek, y= beta, color=epiyr), size=.8, linetype="twodash") +
  geom_line(aes(x=biweek, y= beta, color=epiyr), size=.8) + scale_fill_manual(name="epidemic year", values=c("tomato", "cornflowerblue", "seagreen")) +
  scale_color_manual(name="epidemic year", values=c("tomato", "cornflowerblue", "seagreen")) +
  facet_wrap(provname~., scales = "free_y")+ylab(bquote(beta~', biweekly transmission')) +
  xlab("biweek of year") +
  theme(panel.grid = element_blank(), #legend.position = c(.86,.06),
        strip.background = element_rect(fill="white"), 
        legend.box = "vertical", strip.text = element_text(size=18),
        axis.title = element_text(size=18), 
        axis.text = element_text(size=13))


ggsave(file = paste0(homewd, "/final-figures/poster/Beta-Fig.png"),
       plot = FigS10sub,
       units="mm",  
       width=60, 
       height=40, 
       scale=3, 
       dpi=300)


#plot together by province to see
p1 <- ggplot(data = epi.dat) + theme_bw() + 
  geom_ribbon(aes(x=biweek, ymin= betalow, ymax=betahigh, fill=provname), alpha=.2) + 
  geom_line(aes(x=biweek, y= beta, color=provname), size=.8) +
  facet_wrap(epiyr~., scales = "free_y")+ylab(bquote(beta~', biweekly transmission')) +
  xlab("biweek of year") +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        strip.background = element_rect(fill="white"), legend.title = element_blank(),
        axis.title = element_text(size=18), axis.text = element_text(size=13)) +
  guides(fill=guide_legend(ncol=8),color=guide_legend(ncol=8))
print(p1)

