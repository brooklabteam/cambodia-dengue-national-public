rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(sf)
library(ggh4x)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national-public"

setwd(homewd)

wave.dat <- read.csv(file = paste0(homewd, "/data/wavelet_input_dat_prov.csv"), header = T, stringsAsFactors = F)
head(wave.dat)

#load and attach centroid of each province
centroid.prov <- read.csv(file = paste0(homewd, "/data/centroid_provinces.csv"), header = T, stringsAsFactors = F)
head(centroid.prov)

setdiff(unique(wave.dat$provname), unique(centroid.prov$provname))
setdiff(unique(centroid.prov$provname), unique(wave.dat$provname))
wave.dat <- merge(wave.dat, centroid.prov, by="provname")

#z score script is (x-mu)/sigma

#get z score for the entire data set
mean.temp = mean(wave.dat$temp_C)
sd.temp = sd(wave.dat$temp_C)

mean.precip = mean(wave.dat$precip_mm)
sd.precip = sd(wave.dat$precip_mm)

#and by province
mean.temp.prov = ddply(wave.dat, .(provname), summarise, mean_temp = mean(temp_C), sd_temp = sd(temp_C))
mean.precip.prov = ddply(wave.dat, .(provname), summarise, mean_precip = mean(precip_mm), sd_precip = sd(precip_mm))


#whole dataset
wave.dat$z_score_temp <- (wave.dat$temp_C-mean.temp)/sd.temp
wave.dat$z_score_precip <- (wave.dat$precip_mm-mean.precip)/sd.precip

#and controlling by province
wave.dat <- merge(wave.dat, mean.temp.prov, by="provname")
wave.dat <- merge(wave.dat, mean.precip.prov, by="provname")

#and new Z
wave.dat$z_score_temp_prov <- (wave.dat$temp_C - wave.dat$mean_temp)/wave.dat$sd_temp
wave.dat$z_score_precip_prov <- (wave.dat$precip_mm- wave.dat$mean_precip)/wave.dat$sd_precip

#and plot it out as a heat map, ranked by latitude of the centroid of each province

jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))


wave.dat <- arrange(wave.dat, latitude)
wave.dat$provname <- factor(wave.dat$provname, levels=unique(wave.dat$provname))

vert.df <- cbind.data.frame(xint = c(2007, 2008, 2012, 2013, 2019,2020))

#arranged North at top, South at bottom, by latitude of province
# FigS5B <- ggplot(data=wave.dat) + 
#       geom_tile(aes(x=time, y=provname, fill=z_score_temp)) + 
#       scale_fill_gradientn(colors=jet.colors(7), name="z-score of\nbiweekly\ntemperature") +
#       theme_bw() + theme(panel.grid = element_blank(), axis.title = element_blank()) +
#       geom_vline(data=vert.df, aes(xintercept=xint), color="black", size=1)
#       
#       
# FigS5B


#and controlling within province - better

#get a color bar that matches the one used in figure 1 (there were only 22 provinces there)
# "Tboung Khmum" "Mondul Kiri"  "Ratanak Kiri" need to be moved to the end

colz = scales::hue_pal()(length(unique((wave.dat$provname)))) #25
name.list <- sort(unique(as.character(wave.dat$provname))) #alphabetical
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

FigS5Cb <- ggplot(data=wave.dat) + coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time, y=provname, fill=z_score_temp_prov)) + 
  scale_fill_gradientn(colors=jet.colors(7), name="z-score of\nbiweekly\ntemperature") +
  theme_bw() + theme(panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(0,.2,.2,1.1), "cm"),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.spacing = unit(c(0), "cm"),
                     panel.border = element_rect(linewidth=0)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="black", size=1) 




#and precip
wave.dat$provname <- as.character(wave.dat$provname)
wave.dat <- arrange(wave.dat, latitude)
wave.dat$provname <- factor(wave.dat$provname, levels=unique(wave.dat$provname))

#arranged North at top, South at bottom, by latitude of province
# #using national z score
# FigS6Ab <- ggplot(data=wave.dat) + 
#   geom_tile(aes(x=time, y=provname, fill=z_score_precip)) + 
#   scale_fill_gradientn(colors=jet.colors(7), name="z-score of\nbiweekly\nprecipitation") +
#   theme_bw() + theme(panel.grid = element_blank(), axis.title = element_blank()) + 
#   geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1)
# 
# FigS6Ab
# 

#and controlling within a province - better
FigS6Cb <- ggplot(data=wave.dat) + coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time, y=provname, fill=z_score_precip_prov)) + 
  scale_fill_gradientn(colors=jet.colors(7), name="z-score of\nbiweekly\nprecipitation") +
  theme_bw() + theme(panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(0,.2,.2,1.1), "cm"),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.spacing = unit(c(0), "cm"),
                     panel.border = element_rect(linewidth=0)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) 


#and pair with biweekly distribution of the z-scores median and IQR
dist.dat <- ddply(wave.dat, .(time, year, biweek), summarise, median_temp_z = quantile(z_score_temp_prov)["50%"],  min_temp_z = quantile(z_score_temp_prov)["25%"], max_temp_z = quantile(z_score_temp_prov)["75%"], median_precip_z = quantile(z_score_precip_prov)["50%"], min_precip_z = quantile(z_score_precip_prov)["25%"], max_precip_z = quantile(z_score_precip_prov)["75%"])
head(dist.dat)

#and plot
FigS5Ca <- ggplot(data=dist.dat) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
               axis.title.y = element_text(size=13), axis.text.y = element_text(size=10),
               axis.text.x = element_blank(), axis.ticks.x = element_blank(),
               plot.margin = unit(c(.2,3,0,1.1), "cm")) +
               geom_ribbon(aes(x=time, ymin=min_temp_z, ymax=max_temp_z), alpha=.3) + ylab("biweekly temperature\nz-score distribution") +
               geom_line(aes(x=time, y=median_temp_z), size=1) +
               geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)


FigS6Ca <- ggplot(data=dist.dat) +theme_bw() +  coord_cartesian(xlim=c(2001.9, 2020.1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=13), axis.text.y = element_text(size=10),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.2,3,0,1.1), "cm")) +
  geom_ribbon(aes(x=time, ymin=min_precip_z, ymax=max_precip_z), alpha=.3) + ylab("biweekly precipitation\nz-score distribution") +
  geom_line(aes(x=time, y=median_precip_z), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)



FigS5C <- cowplot::plot_grid(FigS5Ca, FigS5Cb, ncol = 1, nrow = 2, rel_heights = c(.3,1))
FigS6C <- cowplot::plot_grid(FigS6Ca, FigS6Cb, ncol = 1, nrow = 2, rel_heights = c(.3,1))

#and the GAM plots for each
#which years represent significant temperature anomalies?
wave.dat <- arrange(wave.dat, provname, time)
wave.dat$provname <- as.factor(wave.dat$provname)
wave.dat$year <- as.factor(wave.dat$year)
#temp 
gam1 <- gam(temp_C~ s(year, bs="re") +
                    s(biweek, bs="cc", k=7) + #controls for internal annual cycles
                    s(provname, bs="re"), #y-intercept specific by province 
                    data=wave.dat) 
summary(gam1)


#precip
gam2 <- gam(precip_mm~ s(year, bs="re") +
            s(biweek, bs="cc", k=7) + #controls for internal annual cycles
            s(provname, bs="re"), #y-intercept specific by province 
            data=wave.dat) 
summary(gam2)

#and visualize those deviations
source(paste0(homewd, "/figure-development/Fig1/mollentze-streicker-2020-functions.R"))

temp.df <- get_partial_effects(fit=gam1, var="year")
precip.df <- get_partial_effects(fit=gam2, var="year")

plot.partial(temp.df, var="year", response_var = "temp (*C)")
plot.partial(precip.df, var="year", response_var = "precip (mm)")


#and add to the dataset

temp.df = temp.df[[1]]
head(temp.df)
temp.df$color <- "NotSig"
temp.df$color[temp.df$IsSignificant=="Yes" & temp.df$y>0] <- "Pos"
temp.df$color[temp.df$IsSignificant=="Yes" & temp.df$y<0] <- "Neg"
colz2 <- c('NotSig' = "gray60", 'Pos'="red", 'Neg'="blue")

FigS5B <- ggplot(data=temp.df) + theme_bw() +
          theme(panel.grid = element_blank(), axis.title.x = element_blank()) +
          ylab("partial effect of year on temperature")+
          geom_hline(aes(yintercept=0))+
          geom_point(aes(x=year, y=y, color=color), show.legend = F, size=3) +
          geom_linerange(aes(x=year, ymin=ylower, ymax=yupper, color=color), show.legend = F) +
          scale_color_manual(values=colz2) 



precip.df = precip.df[[1]]
head(precip.df)
precip.df$color <- "NotSig"
precip.df$color[precip.df$IsSignificant=="Yes" & precip.df$y>0] <- "Pos"
precip.df$color[precip.df$IsSignificant=="Yes" & precip.df$y<0] <- "Neg"
colz2 <- c('NotSig' = "gray60", 'Pos'="red", 'Neg'="blue")

FigS6B <- ggplot(data=precip.df) + theme_bw() +
  theme(panel.grid = element_blank(), axis.title.x = element_blank()) +
  ylab("partial effect of year on precipitation")+
  geom_hline(aes(yintercept=0))+
  geom_point(aes(x=year, y=y, color=color), show.legend = F, size=3) +
  geom_linerange(aes(x=year, ymin=ylower, ymax=yupper, color=color), show.legend = F) +
  scale_color_manual(values=colz2) 


#and put these all together with map
cam = sf::st_read(paste0(homewd, "/data/province-shape/khm_admbnda_adm1_gov_20181004.shp"))
head(cam)
unique(cam$ADM1_EN)
unique(cam$ADM1_EN)
cam$ADM1_EN[cam$ADM1_EN=="Siemreap"] <- "Siem Reap"


FigCam <- ggplot(cam) +
  geom_sf(aes(fill=ADM1_EN), linewidth=.6, show.legend = F) +
  geom_label(data=centroid.prov, aes(x=longitude, y=latitude, label=provname), size=2, label.size = 0, 
             label.padding =unit(c(.1), "lines")) +
  theme_bw() +  
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_manual(values= colz)



FigS5AB <- cowplot::plot_grid(FigS5B, FigCam, ncol=1, nrow=2, labels = c("A", "B"), label_size = 22, label_y = c(1,.98), label_x = .02) 
FigS6AB <- cowplot::plot_grid(FigS6B, FigCam, ncol=1, nrow=2, labels = c("A", "B"), label_size = 22,  label_y = c(1,.98), label_x = .02)
#FigS5 <- cowplot::plot_grid(FigS5A, FigS5BC, ncol=2, nrow=1, rel_widths = c(1,.8), labels = c("A",""),label_size = 22)

FigS5 <- cowplot::plot_grid(FigS5AB, FigS5C, ncol=2, nrow=1, rel_widths = c(.7,1), labels = c("","C"),label_size = 22,  label_y = c(1,1.005)) + theme(plot.background = element_rect(fill="white"))
FigS6 <- cowplot::plot_grid(FigS6AB, FigS6C, ncol=2, nrow=1, rel_widths = c(.7,1), labels = c("","C"),label_size = 22, label_y = c(1,1.005))+ theme(plot.background = element_rect(fill="white"))


### and save

ggsave(file = paste0(homewd, "/final-figures/FigS5.png"),
       plot= FigS5,
       units="mm",  
       width=120, 
       height=85, 
       scale=3, 
       dpi=300)



ggsave(file = paste0(homewd, "/final-figures/FigS6.png"),
       plot= FigS6,
       units="mm",  
       width=120, 
       height=85, 
       scale=3, 
       dpi=300)


#archive - different orders for the provinces
###################

#first, temperature
   #and arranged by longitude - east (top) to west (bottom)
   wave.dat$provname <- as.character(wave.dat$provname)
   wave.dat <- arrange(wave.dat, longitude)
   wave.dat$provname <- factor(wave.dat$provname, levels=unique(wave.dat$provname))
   p2 <- ggplot(data=wave.dat) + 
     geom_tile(aes(x=time, y=provname, fill=z_score_temp)) + 
     scale_fill_gradientn(colors=jet.colors(7), name="z-score of\nbiweekly\ntemperature") +
     theme_bw() + theme(panel.grid = element_blank(), axis.title = element_blank()) + 
     geom_vline(aes(xintercept=c(2007)), color="red")+
     geom_vline(aes(xintercept=c(2008)), color="red") +
     geom_vline(aes(xintercept=c(2012)), color="red") +
     geom_vline(aes(xintercept=c(2013)), color="red") +
     geom_vline(aes(xintercept=c(2019)), color="red") +
     geom_vline(aes(xintercept=c(2020)), color="red")
   #scale_fill_viridis_c(option="inferno", name="z-score of\nbiweekly temperature") 
   p2
   
   #and by elevation
   wave.dat$provname <- as.character(wave.dat$provname)
   wave.dat <- arrange(wave.dat, mean_elevation_m)
   wave.dat$provname <- factor(wave.dat$provname, levels=unique(wave.dat$provname))
   p3 <- ggplot(data=wave.dat) + 
     geom_tile(aes(x=time, y=provname, fill=z_score_temp)) + 
     scale_fill_gradientn(colors=jet.colors(7), name="z-score of\nbiweekly\ntemperature") +
     theme_bw() + theme(panel.grid = element_blank(), axis.title = element_blank()) + 
     geom_vline(aes(xintercept=c(2007)), color="red")+
     geom_vline(aes(xintercept=c(2008)), color="red") +
     geom_vline(aes(xintercept=c(2012)), color="red") +
     geom_vline(aes(xintercept=c(2013)), color="red") +
     geom_vline(aes(xintercept=c(2019)), color="red") +
     geom_vline(aes(xintercept=c(2020)), color="red")
   #scale_fill_viridis_c(option="inferno", name="z-score of\nbiweekly temperature") 
   p3
   
   
   
#######
# and precip 
   
   #and arranged by longitude -  east (top) to west (bottom)
   wave.dat$provname <- as.character(wave.dat$provname)
   wave.dat <- arrange(wave.dat, longitude)
   wave.dat$provname <- factor(wave.dat$provname, levels=unique(wave.dat$provname))
   p5 <- ggplot(data=wave.dat) + 
     geom_tile(aes(x=time, y=provname, fill=z_score_precip)) + 
     scale_fill_gradientn(colors=jet.colors(7), name="z-score of\nbiweekly\nprecipitation") +
     theme_bw() + theme(panel.grid = element_blank(), axis.title = element_blank()) + 
     geom_vline(aes(xintercept=c(2007)), color="red")+
     geom_vline(aes(xintercept=c(2008)), color="red") +
     geom_vline(aes(xintercept=c(2012)), color="red") +
     geom_vline(aes(xintercept=c(2013)), color="red") +
     geom_vline(aes(xintercept=c(2019)), color="red") +
     geom_vline(aes(xintercept=c(2020)), color="red")
   
   p5
   
   #and by elevation
   wave.dat$provname <- as.character(wave.dat$provname)
   wave.dat <- arrange(wave.dat, mean_elevation_m)
   wave.dat$provname <- factor(wave.dat$provname, levels=unique(wave.dat$provname))
   p6 <- ggplot(data=wave.dat) + 
     geom_tile(aes(x=time, y=provname, fill=z_score_temp)) + 
     scale_fill_gradientn(colors=jet.colors(7), name="z-score of\nbiweekly\nprecipitation") +
     theme_bw() + theme(panel.grid = element_blank(), axis.title = element_blank()) + 
     geom_vline(aes(xintercept=c(2007)), color="red")+
     geom_vline(aes(xintercept=c(2008)), color="red") +
     geom_vline(aes(xintercept=c(2012)), color="red") +
     geom_vline(aes(xintercept=c(2013)), color="red") +
     geom_vline(aes(xintercept=c(2019)), color="red") +
     geom_vline(aes(xintercept=c(2020)), color="red")
   
   p6
   