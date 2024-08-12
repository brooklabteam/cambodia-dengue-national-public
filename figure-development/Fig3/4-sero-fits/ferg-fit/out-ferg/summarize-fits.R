rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/figure-development/Fig3/4-sero-fits/ferg-fit/out-ferg/"))

#load and combine the fits and plot them together (except not the ones with Inf)

load("fit-prov-Banteay-Meanchey.Rdata")
fit.dat <- out

load("fit-prov-Battambang.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kampong-Cham.Rdata")
fit.dat <- rbind(fit.dat, out)


load("fit-prov-Kampong-Chhnang.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kampong-Speu.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kampong-Thom.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kampot.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kandal.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Kep.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Koh-Kong.Rdata")
fit.dat <- rbind(fit.dat, out)


load("fit-prov-Kratie.Rdata")
fit.dat <- rbind(fit.dat, out)


load("fit-prov-Mondul-Kiri.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Otdar-Meanchey.Rdata")
fit.dat <- rbind(fit.dat, out)


load("fit-prov-Pailin.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Phnom-Penh.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Preah-Sihanouk.Rdata")
fit.dat <- rbind(fit.dat, out)


load("fit-prov-Preah-Vihear.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Prey-Veng.Rdata")
fit.dat <- rbind(fit.dat, out)


load("fit-prov-Pursat.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Ratanak-Kiri.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Siem-Reap.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Stung-Treng.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Svay-Rieng.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Takeo.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-prov-Tboung-Khmum.Rdata")
fit.dat <- rbind(fit.dat, out)

load("fit-national.Rdata")
fit.dat <- rbind(fit.dat, out)


head(fit.dat)
fit.dat$provname[fit.dat$provname=="national"] <- "National"
unique(fit.dat$convergence) #all converged!
unique(fit.dat$provname[fit.dat$convergence==1])# 9 provinces..."Kampong Cham" "Kampong Speu" "Kampong Thom" "Kampot"       "Kandal"       "Phnom Penh"  
#fit.dat <- subset(fit.dat, lambda<10)

colz = scales::hue_pal()(length(unique(fit.dat$provname)))
names(colz) <- unique(fit.dat$provname)

colz[names(colz)=="National"] <- "black"
colz[names(colz)=="Phnom Penh"] <- "red"

fit.dat$lci[fit.dat$lci=="not yet"] <- NA
fit.dat$uci[fit.dat$uci=="not yet"] <- NA
fit.dat$lci <-as.numeric(fit.dat$lci)
fit.dat$uci <-as.numeric(fit.dat$uci)

#khmer.rouge = cbind.data.frame(year=c(1975:1979), ymin=0, ymax=1)


p1 <- ggplot(subset(fit.dat, provname!="Tboung Khmum")) + 
      #geom_ribbon(data=khmer.rouge, aes(x=year, ymin=ymin, ymax=ymax), fill="cornflowerblue", alpha=.3)+
      geom_point(aes(x=year, y=lambda, color=provname)) +
      geom_ribbon(aes(x=year, ymin=lci, ymax=uci, fill=provname), alpha=.3) +    
      geom_line(aes(x=year, y=lambda, color=provname)) +theme_bw() + 
      geom_point(data = subset(fit.dat, provname=="Phnom Penh"), aes(x=year, y=lambda), color="red") +
      geom_ribbon(data = subset(fit.dat, provname=="Phnom Penh"),aes(x=year, ymin=lci, ymax=uci), alpha=.3, fill="red") +    
      geom_line(data = subset(fit.dat, provname=="Phnom Penh"),aes(x=year, y=lambda), size =.8, color="red")+
      geom_point(data = subset(fit.dat, provname=="National"), aes(x=year, y=lambda), color="black") +
      geom_ribbon(data = subset(fit.dat, provname=="National"),aes(x=year, ymin=lci, ymax=uci), alpha=.3, fill="black") +    
      geom_line(data = subset(fit.dat, provname=="National"),aes(x=year, y=lambda), size =1, color="black") +
      theme(panel.grid = element_blank()) + #coord_cartesian(ylim = c(0,1)) +
      geom_vline(xintercept = 2007, linetype=2) + scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
      geom_vline(xintercept = 2012, linetype=2) +
      geom_vline(xintercept = 2019, linetype=2) 
      

print(p1)

ggsave(filename = "FOI-by-province.png",
       plot = p1,
       units="mm",  
       width=85, 
       height=50, 
       scale=3, 
       dpi=300)


#save data 
save(fit.dat, file = "prov-fits-FOI.Rdata")

p2 <- ggplot(fit.dat) + 
  geom_point(aes(x=year, y=lambda, color=provname), show.legend = F) +
  geom_ribbon(aes(x=year, ymin=lci, ymax=uci, fill=provname), alpha=.3,show.legend = F) +    
  geom_line(aes(x=year, y=lambda, color=provname), show.legend = F) +theme_bw() + 
  geom_point(data = subset(fit.dat, provname=="Phnom Penh"), aes(x=year, y=lambda), color="red") +
  geom_ribbon(data = subset(fit.dat, provname=="Phnom Penh"),aes(x=year, ymin=lci, ymax=uci), alpha=.3, fill="red") +    
  geom_line(data = subset(fit.dat, provname=="Phnom Penh"),aes(x=year, y=lambda), size =.8, color="red")+
  geom_point(data = subset(fit.dat, provname=="National"), aes(x=year, y=lambda), color="black") +
  geom_ribbon(data = subset(fit.dat, provname=="National"),aes(x=year, ymin=lci, ymax=uci), alpha=.3, fill="black") +    
  geom_line(data = subset(fit.dat, provname=="National"),aes(x=year, y=lambda), size =1, color="black") +
  theme(panel.grid = element_blank(), strip.background =  element_rect(fill="white")) + #coord_cartesian(ylim = c(0,1)) +
  geom_vline(xintercept = 2007, linetype=2) + scale_color_manual(values=colz) +
  geom_vline(xintercept = 2012, linetype=2) +
  geom_vline(xintercept = 2019, linetype=2) +facet_wrap(~provname, ncol=4)

print(p2)

ggsave(filename = "FOI-by-province-facet.png",
       plot = p2,
       units="mm",  
       width=100, 
       height=90, 
       scale=3, 
       dpi=300)


min.all <- ddply(fit.dat, .(provname), summarise, min_year = min(year))
sort(unique(min.all$min_year)) #1966 forTboung Khmum and 1981 after that
