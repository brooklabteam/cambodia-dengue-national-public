rm(list=ls())


library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(epitools)
library(reshape2)
library(ggh4x)

#this one is just the plot of coherence with ONI (or lack thereof)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national-public/"
setwd(homewd)

oni.dat <- read.csv(file= paste0(homewd, "/data/oni_coherence.csv"), header = T, stringsAsFactors = F)
clim.dat <- read.csv(file= paste0(homewd, "/data/synchrony_climate_data_oct15.csv"), header = T, stringsAsFactors = F)


#load and attach centroid of each province
centroid.prov <- read.csv(file = paste0(homewd, "/data/centroid_provinces.csv"), header = T, stringsAsFactors = F)
head(centroid.prov)

setdiff(unique(oni.dat$provname), unique(centroid.prov$provname))
setdiff(unique(centroid.prov$provname), unique(oni.dat$provname))
setdiff(unique(clim.dat$provname), unique(centroid.prov$provname))
setdiff(unique(centroid.prov$provname), unique(clim.dat$provname))
oni.dat$provname[oni.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"

oni.dat <- merge(oni.dat, centroid.prov, by="provname")
clim.dat <- merge(clim.dat, centroid.prov, by="provname")

#and plot by latitude
oni.dat <- arrange(oni.dat, latitude, time)
clim.dat <- arrange(clim.dat, latitude, time)

oni.dat$provname <- factor(oni.dat$provname, levels=unique(oni.dat$provname))
clim.dat$provname <- factor(clim.dat$provname, levels=unique(clim.dat$provname))


#now make color bar
colz = scales::hue_pal()(length(unique((oni.dat$provname)))) #25
name.list <- sort(unique(as.character(oni.dat$provname))) #alphabetical
name.list <- c(name.list[name.list!="Ratanak Kiri" & name.list!="Mondul Kiri" & name.list!= "Tboung Khmum"], "Mondul Kiri", "Ratanak Kiri", "Tboung Khmum") #these three get moved to the end
names(colz) <- name.list

#now, reorder it in the same order that the provinces are plotted.
centroid.prov$color= NA
for(i in 1:length(colz)){
  centroid.prov$color[centroid.prov$provname==names(colz)[i]] <- colz[i]
}
centroid.prov <- arrange(centroid.prov, latitude)
colz = centroid.prov$color
names(colz) <- centroid.prov$provname

# Only colour strips in x-direction
strip <- strip_themed(background_y = elem_list_rect(fill = colz))

#and coherence
vert.df <- cbind.data.frame(xint = c(2007, 2008, 2012, 2013, 2019,2020))

#ggplot(subset(oni.dat, provname=="Battambang")) + geom_line(aes(x=time, y=avg_wave_coherence_oni))

#first, cross-wave power, annual cases and temperature
FigS13Ab <- ggplot(data=clim.dat) +  coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time, y=provname, fill=avg_cross_power_temp, color=avg_cross_power_temp), show.legend = F) + 
  scale_fill_viridis_c( option="inferno", name="Average cross-wavelet\npower (annual cases):\ndengue IR and temperature", trans="sqrt", limits=c(0,20)) +
  scale_color_viridis_c( option="inferno", name="Average cross-wavelet\npower (annual cases):\ndengue IR and temperature", trans="sqrt", limits=c(0,20)) +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(.1,.5,.1,1), "cm"),
                     strip.text = element_text(size=8),# legend.position = "bottom",
                     panel.spacing = unit(c(0), "cm"),
                     #legend.title = element_text(size=9),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.border = element_rect(linewidth=0),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) 


#and the distribution by month
dist.dat.clim <- ddply(clim.dat, .(time, year), summarise, median_temp_cross = quantile(avg_cross_power_temp, na.rm=T)["50%"],  min_temp_cross = quantile(avg_cross_power_temp, na.rm=T)["25%"], max_temp_cross = quantile(avg_cross_power_temp, na.rm=T)["75%"],  median_precip_cross = quantile(avg_cross_power_precip, na.rm=T)["50%"],  min_precip_cross = quantile(avg_cross_power_precip, na.rm=T)["25%"], max_precip_cross = quantile(avg_cross_power_precip, na.rm=T)["75%"],  median_temp_cross_multi_recons = quantile(avg_cross_power_temp_multi_recons, na.rm=T)["50%"],  min_temp_cross_multi_recons = quantile(avg_cross_power_temp_multi_recons, na.rm=T)["25%"], max_temp_cross_multi = quantile(avg_cross_power_temp_multi_recons, na.rm=T)["75%"], median_precip_cross_multi_recons = quantile(avg_cross_power_precip_multi_recons, na.rm=T)["50%"],  min_precip_cross_multi_recons = quantile(avg_cross_power_precip_multi_recons, na.rm=T)["25%"], max_precip_cross_multi = quantile(avg_cross_power_precip_multi_recons, na.rm=T)["75%"])
head(dist.dat.clim)
tail(dist.dat.clim)

#and plot
FigS13Aa <- ggplot(data=dist.dat.clim) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.3,.5,.1,1), "cm")) +
  geom_ribbon(aes(x=time, ymin=min_temp_cross, ymax=max_temp_cross), alpha=.3) + ylab("distribution biweekly\nannual dengue IR\ncross-wavelet power\nwith temperature") +
  geom_line(aes(x=time, y=median_temp_cross), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)

FigS13A <- cowplot::plot_grid(FigS13Aa, FigS13Ab, ncol=1, nrow=2, rel_heights = c(.25,1))


#and b is annual cases and precip
FigS13Bb <- ggplot(data=clim.dat) +  coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time, y=provname, fill=avg_cross_power_precip, color=avg_cross_power_precip)) + 
  scale_fill_viridis_c( option="inferno", name="Average cross-wavelet power:\ndengue annual IR and\ntemperature (A) or\nprecipitation (B)", trans="sqrt", limits=c(0,20)) +
  scale_color_viridis_c( option="inferno", name="Average cross-wavelet power:\ndengue annual IR and\ntemperature (A) or\nprecipitation (B)", trans="sqrt", limits=c(0,20))+
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(.1,.5,.1,1), "cm"),
                     strip.text = element_text(size=8), legend.position = "bottom",
                     panel.spacing = unit(c(0), "cm"),
                     legend.title = element_text(size=9),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.border = element_rect(linewidth=0),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) 


#and the distribution by month
dist.dat.clim <- ddply(clim.dat, .(time, year), summarise, median_temp_cross = quantile(avg_cross_power_temp, na.rm=T)["50%"],  min_temp_cross = quantile(avg_cross_power_temp, na.rm=T)["25%"], max_temp_cross = quantile(avg_cross_power_temp, na.rm=T)["75%"],  median_precip_cross = quantile(avg_cross_power_precip, na.rm=T)["50%"],  min_precip_cross = quantile(avg_cross_power_precip, na.rm=T)["25%"], max_precip_cross = quantile(avg_cross_power_precip, na.rm=T)["75%"],  median_temp_cross_multi_recons = quantile(avg_cross_power_temp_multi_recons, na.rm=T)["50%"],  min_temp_cross_multi_recons = quantile(avg_cross_power_temp_multi_recons, na.rm=T)["25%"], max_temp_cross_multi = quantile(avg_cross_power_temp_multi_recons, na.rm=T)["75%"], median_precip_cross_multi_recons = quantile(avg_cross_power_precip_multi_recons, na.rm=T)["50%"],  min_precip_cross_multi_recons = quantile(avg_cross_power_precip_multi_recons, na.rm=T)["25%"], max_precip_cross_multi = quantile(avg_cross_power_precip_multi_recons, na.rm=T)["75%"])
head(dist.dat.clim)
tail(dist.dat.clim)

#and plot
FigS13Ba <- ggplot(data=dist.dat.clim) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.3,.5,.1,1), "cm")) +
  geom_ribbon(aes(x=time, ymin=min_precip_cross, ymax=max_precip_cross), alpha=.3) + ylab("distribution biweekly\nannual dengue IR\ncross-wavelet power\nwith precipitation") +
  geom_line(aes(x=time, y=median_precip_cross), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)

FigS13B <- cowplot::plot_grid(FigS13Ba, FigS13Bb, ncol=1, nrow=2, rel_heights = c(.25,1))


#and cross wave power of multi-annual cycles and temperature


#first, cross-wave power, annual cases and temperature
FigS13Cb <- ggplot(data=clim.dat) +  coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time, y=provname, fill=avg_cross_power_temp_multi_recons, color=avg_cross_power_temp_multi_recons), show.legend = F) + 
  scale_fill_viridis_c( option="inferno", name="Average cross-wavelet\npower (multiannual):\ndengue cycles and temperature", trans="sqrt",limits=c(0,5)) +
  scale_color_viridis_c( option="inferno", name="Average cross-wavelet\npower (multiannual):\ndengue cycles and temperature", trans="sqrt",limits=c(0,5)) +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(.1,.5,.1,1), "cm"),
                     strip.text = element_text(size=8), legend.position = "bottom",
                     panel.spacing = unit(c(0), "cm"),
                     legend.title = element_text(size=9),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.border = element_rect(linewidth=0),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) 



#and plot
FigS13Ca <- ggplot(data=dist.dat.clim) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.3,.5,.1,1), "cm")) +
  geom_ribbon(aes(x=time, ymin=min_temp_cross_multi_recons, ymax=max_temp_cross_multi), alpha=.3) + ylab("distribution biweekly\nmulti-annual dengue\ncycles cross-wavelet\npower with temperature") +
  geom_line(aes(x=time, y=median_temp_cross_multi_recons), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)

FigS13C <- cowplot::plot_grid(FigS13Ca, FigS13Cb, ncol=1, nrow=2, rel_heights = c(.25,1))

#now, cross power multi-annual cycles with precip
FigS13Db <- ggplot(data=clim.dat) +  coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time, y=provname, fill=avg_cross_power_precip_multi_recons, color=avg_cross_power_precip_multi_recons)) + 
  scale_fill_viridis_c( option="inferno", name="Average cross-wavelet power:\nmultiannual dengue cycles\nand temperature (C) or\nprecipitation (B)", trans="sqrt",limits=c(0,5)) +
  scale_color_viridis_c( option="inferno", name="Average cross-wavelet power:\nmultiannual dengue cycles\nand temperature (C) or\nprecipitation (B)", trans="sqrt",limits=c(0,5)) +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(.1,.5,.1,1), "cm"),
                     strip.text = element_text(size=8), legend.position = "bottom",
                     panel.spacing = unit(c(0), "cm"),
                     legend.title = element_text(size=9),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.border = element_rect(linewidth=0),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) 



#and plot
FigS13Da <- ggplot(data=dist.dat.clim) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.3,.5,.1,1), "cm")) +
  geom_ribbon(aes(x=time, ymin=min_precip_cross_multi_recons, ymax=max_precip_cross_multi), alpha=.3) + ylab("distribution biweekly\nmulti-annual dengue\ncycles cross-wavelet\npower with precipitation") +
  geom_line(aes(x=time, y=median_precip_cross_multi_recons), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)

FigS13D <- cowplot::plot_grid(FigS13Da, FigS13Db, ncol=1, nrow=2, rel_heights = c(.25,1))



#cross wavelet power - ONI
FigS13Eb <- ggplot(data=oni.dat) +  coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time, y=provname, fill=avg_cross_power_oni, color=avg_cross_power_oni)) + 
  scale_fill_viridis_c( option="inferno", name="Average cross-wavelet\npower (multiannual):\ndengue IR and ONI", trans="sqrt") +# limits=c(0.5,1)) +
  scale_color_viridis_c( option="inferno", name="Average cross-wavelet\npower (multiannual):\ndengue IR and ONI", trans="sqrt")+#, limits=c(0.5,1)) +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(.1,.5,.1,1), "cm"),
                     strip.text = element_text(size=8), legend.position = "bottom",
                     panel.spacing = unit(c(0), "cm"),
                     legend.title = element_text(size=9),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.border = element_rect(linewidth=0),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) 


#and the distribution by month - oni
dist.dat <- ddply(oni.dat, .(time, year, month), summarise, median_oni_cross = quantile(avg_cross_power_oni, na.rm=T)["50%"],  min_oni_cross = quantile(avg_cross_power_oni, na.rm=T)["25%"], max_oni_cross = quantile(avg_cross_power_oni, na.rm=T)["75%"], median_oni_coherence = quantile(avg_wave_coherence_oni, na.rm=T)["50%"],  min_oni_coherence = quantile(avg_wave_coherence_oni, na.rm=T)["25%"], max_oni_coherence = quantile(avg_wave_coherence_oni, na.rm=T)["75%"])
head(dist.dat)
tail(dist.dat)

#and plot
FigS13Ea <- ggplot(data=dist.dat) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.3,.5,.1,1), "cm")) +
  geom_ribbon(aes(x=time, ymin=min_oni_cross, ymax=max_oni_cross), alpha=.3) + ylab("distribution monthly\nmulti-annual dengue\ncycle cross-wavelet\npower with ONI") +
  geom_line(aes(x=time, y=median_oni_cross), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)

FigS13E <- cowplot::plot_grid(FigS13Ea, FigS13Eb, ncol=1, nrow=2, rel_heights = c(.25,1))

FigS13AB <- cowplot::plot_grid(FigS13A, FigS13B, nrow=2, ncol=1, labels = c("A", "B"), rel_heights = c(1,1.2), label_size = 22)
FigS13CD <- cowplot::plot_grid(FigS13C, FigS13D, nrow=2, ncol=1, labels = c("C", "D"), rel_heights = c(1,1.2), label_size = 22)
FigS13EF <- cowplot::plot_grid(FigS13E, nrow=2, ncol=1, labels = c("E", ""), rel_heights = c(1.06,1), label_size = 22)

FigS13 <- cowplot::plot_grid(FigS13AB, FigS13CD, FigS13EF, ncol=3, nrow=1) + theme(plot.background = element_rect(fill="white"))


ggsave(file = paste0(homewd, "/final-figures/FigS13.png"),
       plot= FigS13,
       units="mm",  
       width=150, 
       height=110, 
       scale=3, 
       dpi=300)
