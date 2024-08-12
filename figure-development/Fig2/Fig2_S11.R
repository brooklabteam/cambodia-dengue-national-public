rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(reshape2)
library(patchwork)
library(RColorBrewer)
library(ggh4x)


homewd = "/Users/carabrook/Developer/cambodia-dengue-national-public"
setwd(homewd)


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
  max.period.list <- lapply(col.split, get.max.period)
  
  nat.dat$multi_period_mean <-  c(unlist(mean.period.list))
  nat.dat$multi_period_max <-  c(unlist(max.period.list))
  
  return(nat.dat)
}


dat <- read.csv(file = paste0(homewd, "/data/synchrony_case_data_oct15.csv"), header=T, stringsAsFactors = F)
head(dat)
names(dat)
dat$provname[dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"


#load and attach centroid of each province
centroid.prov <- read.csv(file = paste0(homewd, "/data/centroid_provinces.csv"), header = T, stringsAsFactors = F)
head(centroid.prov)

setdiff(unique(dat$provname), unique(centroid.prov$provname))
setdiff(unique(centroid.prov$provname), unique(dat$provname))

dat <- merge(dat, centroid.prov, by="provname")

#and plot by latitude
dat <- arrange(dat, latitude, time)

dat$provname <- factor(dat$provname, levels=unique(dat$provname))

#display.brewer.all()


# define jet colormap
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# 
# Fig2Ab <- ggplot(data=dat) + theme_bw() +
#   theme(panel.grid = element_blank(), axis.title = element_blank(),
#         panel.background = element_rect(fill="gray30"),
#         axis.text = element_text(size=14), legend.position = "bottom") +
#   scale_fill_gradientn(colors = jet.colors(7), name="Reconstructed\nannual cycles", limits=c(-2,4)) +
#   scale_color_gradientn(colors = jet.colors(7), name="Reconstructed\nannual cycles", limits=c(-2,4)) +
#   geom_tile(aes(x=time, y=provname, fill=reconstructed_annual_period, color=reconstructed_annual_period))#, color=epiyear))
# 


colz = scales::hue_pal()(length(unique((dat$provname)))) #25
name.list <- sort(unique(as.character(dat$provname))) #alphabetical
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
vert.df <- cbind.data.frame(xint = c(2007, 2008, 2012, 2013, 2019,2020))
# 
# Fig2Ab <- ggplot(data=dat) +  coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
#   facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
#   geom_tile(aes(x=time,  y=provname, fill=scale_annual, color=scale_annual), show.legend = FALSE) + 
#   scale_fill_gradientn(colors = jet.colors, name="Reconstructed\nannual cycles", trans="log10")+#limits=c(-100,280)) +
#   scale_color_gradientn(colors = jet.colors, name="Reconstructed\nannual cycles", trans="log10")+# limits=c(-100,280)) +
#   theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
#                      plot.margin = unit(c(.1,.5,.1,1), "cm"),
#                      strip.text = element_text(size=8), legend.position = "bottom",
#                      panel.spacing = unit(c(0), "cm"),
#                      legend.title = element_text(size=9),
#                      strip.text.y.left = element_text(angle=0, size=6),
#                      panel.border = element_rect(linewidth=0),
#                      panel.grid = element_blank(), axis.title = element_blank(),
#                      axis.text = element_text(size=14)) +
#   geom_vline(data=vert.df, aes(xintercept=xint), color="black", size=1) 
# 
# 
# 


colorz<-c("#00007F","#0000F0", "#0062FF", "#00D4FF", "#FFF800", "#FFC200","#FF5600", "#FF2100","#C80000", "#7F0000")
Fig2Ab <- ggplot(data=dat) +  coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time,  y=provname, fill=reconstructed_annual_period, color=reconstructed_annual_period)) + 
  scale_fill_stepsn(colours = colorz, name="Reconstructed\nannual cycles", breaks=c(-80, -50, -1, 0,1,3,10,30,70,300), values=rescale(x=c(-80, -50, -1, 0,1,3,10,30,70,300), to=c(0,1)))+
  scale_color_stepsn(colours = colorz, name="Reconstructed\nannual cycles", breaks=c(-80, -50, -1,0,1,3,10,30,70,300), values=rescale(x=c(-80, -50, -1, 0,1,3,10,30,70,300), to=c(0,1)))+
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(.1,.5,.1,1), "cm"),
                     strip.text = element_text(size=8), legend.position = "bottom",
                     panel.spacing = unit(c(0), "cm"),
                     legend.title = element_text(size=9),
                     legend.text  = element_text(size=7),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.border = element_rect(linewidth=0),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="black", size=1) 

Fig2Ableg <- cowplot::get_legend(Fig2Ab)

Fig2Abnoleg <- ggplot(data=dat) +  coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time,  y=provname, fill=reconstructed_annual_period, color=reconstructed_annual_period), show.legend = F) + 
  scale_fill_stepsn(colours = colorz, name="Reconstructed\nannual cycles", breaks=c(-80, -50, -1, 0,1,3,10,30,70,300), values=rescale(x=c(-80, -50, -1, 0,1,3,10,30,70,300), to=c(0,1)))+
  scale_color_stepsn(colours = colorz, name="Reconstructed\nannual cycles", breaks=c(-80, -50, -1,0,1,3,10,30,70,300), values=rescale(x=c(-80, -50, -1, 0,1,3,10,30,70,300), to=c(0,1)))+
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(.1,.5,.1,1), "cm"),
                     strip.text = element_text(size=8), legend.position = "bottom",
                     panel.spacing = unit(c(0), "cm"),
                     legend.title = element_text(size=9),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.border = element_rect(linewidth=0),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="black", size=1) 


colorz2<-c("#00007F","#0000F0", "#0062FF", "#7FFF7F", "#FFF800", "#FFC200","#FF5600", "#FF2100","#C80000", "#7F0000")
Fig2Bb <- ggplot(data=dat) + coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time,  y=provname, fill=reconstructed_multi_period, color=reconstructed_multi_period)) + 
  scale_fill_stepsn(colours = colorz2, name="Reconstructed\nmultiannual cycles", breaks=c(-85, -50,-5, 0,3,10,30,50,70,100), values=rescale(x=c(-85, -50,-5,0,3,10,30,50,70,100), to=c(0,1)))+
  scale_color_stepsn(colours = colorz2, name="Reconstructed\nmultiannual cycles", breaks=c(-85, -50,-5,0,3,10,30,50,70,100), values=rescale(x=c(-85, -50,-5,0,3,10,30,50,70,100), to=c(0,1)))+
  #scale_fill_gradientn(colors = jet.colors(7), name="Reconstructed cycles\n(top: annual; bottom: multi-annual)",limits=c(-100,280)) +
  #scale_color_gradientn(colors = jet.colors(7), name="Reconstructed cycles\n(top: annual; bottom: multi-annual)",limits=c(-100,280)) +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(.1,.5,.1,1), "cm"),
                     strip.text = element_text(size=8),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.spacing = unit(c(0), "cm"),
                     legend.title = element_text(size=9),
                     legend.text  = element_text(size=6.5),
                     legend.position = "bottom",
                     panel.border = element_rect(linewidth=0),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     #panel.background = element_rect(fill="gray30"),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="black", size=1) 

Fig2Bbleg <- cowplot::get_legend(Fig2Bb )

Fig2Bbnoleg <- ggplot(data=dat) + coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time,  y=provname, fill=reconstructed_multi_period, color=reconstructed_multi_period), show.legend = F) + 
  scale_fill_stepsn(colours = colorz2, name="Reconstructed\nmultiannual cycles", breaks=c(-85, -50,-5, 0,3,10,30,50,70,100), values=rescale(x=c(-85, -50,-5,0,3,10,30,50,70,100), to=c(0,1)))+
  scale_color_stepsn(colours = colorz2, name="Reconstructed\nmultiannual cycles", breaks=c(-85, -50,-5,0,3,10,30,50,70,100), values=rescale(x=c(-85, -50,-5,0,3,10,30,50,70,100), to=c(0,1)))+
  #scale_fill_gradientn(colors = jet.colors(7), name="Reconstructed cycles\n(top: annual; bottom: multi-annual)",limits=c(-100,280)) +
  #scale_color_gradientn(colors = jet.colors(7), name="Reconstructed cycles\n(top: annual; bottom: multi-annual)",limits=c(-100,280)) +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(.1,.5,.1,1), "cm"),
                     strip.text = element_text(size=8),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.spacing = unit(c(0), "cm"),
                     legend.title = element_text(size=9),
                     legend.position = "bottom",
                     panel.border = element_rect(linewidth=0),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     #panel.background = element_rect(fill="gray30"),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="black", size=1) 






#and get the monthly distribution of annual and biannual cases

dist.dat <- ddply(dat, .(time, year, biweek), summarise, median_annual = quantile(reconstructed_annual_period)["50%"],  min_annual = quantile(reconstructed_annual_period)["25%"], max_annual = quantile(reconstructed_annual_period)["75%"], median_multi = quantile(reconstructed_multi_period)["50%"], min_multi = quantile(reconstructed_multi_period)["25%"], max_multi = quantile(reconstructed_multi_period)["75%"])
head(dist.dat)



#and plot
Fig2Aa <- ggplot(data=dist.dat) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2021.1), ylim=c(-15,40), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.3,.5,0,.9), "cm")) +
  geom_ribbon(aes(x=time, ymin=min_annual, ymax=max_annual), alpha=.3) + ylab("biweekly annual\nreconstructed\ncycle distribution") +
  geom_line(aes(x=time, y=median_annual), size=1) +
  #scale_y_continuous(breaks=c(-.1,0,.2,.4)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)


Fig2Ba <- ggplot(data=dist.dat) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.2,.5,0,.8), "cm")) +
  geom_ribbon(aes(x=time, ymin=min_multi, ymax=max_multi), alpha=.3) + ylab("biweekly multi-annual\nreconstructed\ncycle distribution") +
  geom_line(aes(x=time, y=median_multi), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)


#and construct
Fig2A <- cowplot::plot_grid(Fig2Aa, Fig2Abnoleg, ncol=1, nrow=2, rel_heights = c(.25,1))
Fig2B <- cowplot::plot_grid(Fig2Ba, Fig2Bbnoleg, ncol=1, nrow=2, rel_heights = c(.25,1))


#add legends
Fig2ABlegends <- cowplot::plot_grid(Fig2Ableg, Fig2Bbleg, ncol=2, nrow=1) + theme(plot.margin = unit(c(0,0,0,0), "cm")) 



#and plot as A/B next to coherence
Fig2A2B <- cowplot::plot_grid(Fig2A, Fig2B, Fig2ABlegends, ncol=1, nrow=3, labels = c("A", "B", ""), label_size = 22, rel_heights = c(1,1,.1), label_x = -0.01)


# Now add proportion of provinces with significant annual correlation - by pearsons correlation
# load data
pearsons.df <- read.csv(file = paste0(homewd, "/data/pearsons_correlations_provinces.csv"), header = T, stringsAsFactors = F)


#then, look at avg synchrony per province per year
pearsons.avg.year <- ddply(pearsons.df, .(provname, year), summarise, mean_corr=mean(corr, na.rm=T))
pearsons.avg.year <- merge(pearsons.avg.year, centroid.prov, by="provname")

pearsons.avg.year <- arrange(pearsons.avg.year, latitude, year)
pearsons.avg.year$provname <- factor(pearsons.avg.year$provname, levels=unique(pearsons.avg.year$provname))
unique(pearsons.avg.year$provname)
vert.df2 <- cbind.data.frame(xint = c(2006.5, 2007.5, 2011.5, 2012.5,  2018.5, 2019.5))

# Only colour strips in x-direction
strip <- strip_themed(background_y = elem_list_rect(fill = colz))

#plot as heatmap
pearsons.avg.year$year_plot=pearsons.avg.year$year + .5

Fig2Cb <- ggplot(data=pearsons.avg.year) +  coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=year_plot, y=provname, fill=mean_corr, color=mean_corr), show.legend = F) +
  scale_fill_viridis_c( option="inferno", name="Average annual pairwise province Pearson's correlation coefficient", na.value = "black", limits=c(0,1)) +
  scale_color_viridis_c( option="inferno", name="Average annual pairwise province Pearson's correlation coefficient", na.value = "black", limits=c(0,1)) +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(.1,.5,.1,1), "cm"),
                     strip.text = element_text(size=8), legend.position = "bottom",
                     panel.spacing = unit(c(0), "cm"),
                     legend.title = element_text(size=9),
                     legend.title.align = 0,
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.border = element_rect(linewidth=0),
                     panel.background = element_rect(fill="gray50"),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=.8) 

#and the annual distributions
pearsons.df$year_plot <- pearsons.df$year + 0.5
dist.dat.D <- ddply(pearsons.df, .(year, year_plot), summarise, median_corr = quantile(corr, na.rm=T)["50%"],  min_corr = quantile(corr, na.rm=T)["25%"], max_corr = quantile(corr, na.rm=T)["75%"])
head(dist.dat.D)
dist.dat.D$year <- as.numeric(as.character(dist.dat.D$year))
dist.dat.D$year_plot <- as.numeric(as.character(dist.dat.D$year_plot))
dist.dat.D <- rbind(dist.dat.D, dist.dat.D[nrow(dist.dat.D),])
dist.dat.D$year[length(dist.dat.D$year)] <- 2021
dist.dat.D$year_plot[length(dist.dat.D$year_plot)] <- 2021
#and plot
Fig2Ca <- ggplot(data=dist.dat.D) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.3,.5,.1,1), "cm")) +
  geom_ribbon(aes(x=year_plot, ymin=min_corr, ymax=max_corr), alpha=.3) + ylab("annual distribution\npairwise province\nPearson's correlation ") +
  geom_line(aes(x=year_plot, y=median_corr), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)


Fig2C <- cowplot::plot_grid(Fig2Ca, Fig2Cb, ncol=1, nrow=2, rel_heights = c(.25,1))


# and now pearsons for 5-year moving multi-annual cycles
pearsons.df.multi <- read.csv(file = paste0(homewd, "/data/pearsons_correlations_provinces_multi.csv"), header = T, stringsAsFactors = F)
head(pearsons.df.multi)
tail(pearsons.df.multi)

#then, look at avg synchrony per province per mid_year
pearsons.avg.year.multi <- ddply(pearsons.df.multi, .(provname, mid_year, year_range), summarise, mean_corr=mean(corr, na.rm=T))
pearsons.avg.year.multi <- merge(pearsons.avg.year.multi, centroid.prov, by="provname")

pearsons.avg.year.multi <- arrange(pearsons.avg.year.multi, latitude, mid_year)
pearsons.avg.year.multi$provname <- factor(pearsons.avg.year.multi$provname, levels=unique(pearsons.avg.year.multi$provname))

vert.df2 <- cbind.data.frame(xint = c(2006.5, 2007.5, 2011.5, 2012.5,  2018.5, 2019.5))

# Only colour strips in x-direction
strip <- strip_themed(background_y = elem_list_rect(fill = colz))

#plot as heatmap with overlapping years
pearsons.avg.year.multi$min_year <- as.numeric(sapply(strsplit(pearsons.avg.year.multi$year_range, "-"), '[', 1))
pearsons.avg.year.multi$max_year <- as.numeric(sapply(strsplit(pearsons.avg.year.multi$year_range, "-"), '[', 2))
pearsons.avg.year.multi$max_year <- pearsons.avg.year.multi$max_year +1
#and get the two other years in the series
pearsons.avg.year.multi$year_2 <- pearsons.avg.year.multi$min_year +1
pearsons.avg.year.multi$year_4 <- pearsons.avg.year.multi$mid_year +1


Fig2Db <- ggplot(data=pearsons.avg.year.multi) +  coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=min_year, y=provname, fill=mean_corr, color=mean_corr), alpha=.5) +
  geom_tile(aes(x=year_2, y=provname, fill=mean_corr, color=mean_corr), alpha=.5) +
  geom_tile(aes(x=mid_year, y=provname, fill=mean_corr, color=mean_corr), alpha=.5) +
  geom_tile(aes(x=year_4, y=provname, fill=mean_corr, color=mean_corr), alpha=.5) +
  geom_tile(aes(x=max_year, y=provname, fill=mean_corr, color=mean_corr), alpha=.5) +
  scale_fill_viridis_c( option="inferno", name="Average pairwise province Pearson's correlation\ncoefficient (top: annual; bottom: multi-annual)", na.value = "black", limits=c(0,1)) +
  scale_color_viridis_c( option="inferno", name="Average pairwise province Pearson's correlation\ncoefficient (top: annual; bottom: multi-annual)", na.value = "black", limits=c(0,1)) +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(.1,.5,.1,1), "cm"),
                     strip.text = element_text(size=8), legend.position = "bottom",
                     panel.spacing = unit(c(0), "cm"),
                     legend.title = element_text(size=9),
                     legend.title.align = 0,
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.border = element_rect(linewidth=0),
                     #panel.background = element_rect(fill="gray50"),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=.8) 

#and the multiannual distributions
# to rep every year in the dataset, get a separate value for each year of the 5 year period

pearsons.df.multi.split <- dlply(pearsons.df.multi, .(mid_year))

expand.df <- function(df){
  
  df$min_year <- as.numeric(sapply(strsplit(df$year_range, "-"), '[', 1))
  df$max_year <- as.numeric(sapply(strsplit(df$year_range, "-"), '[', 2))
  df$max_year <- df$max_year +1
  
  year_range = unique(df$min_year):unique(df$max_year)
  
  df.1 <- df
  df.2 <- df
  df.3 <- df
  df.4 <- df
  df.5 <- df
  
  df.1$plot_year <- year_range[1]
  df.2$plot_year <- year_range[2]
  df.3$plot_year <- year_range[3]
  df.4$plot_year <- year_range[4]
  df.5$plot_year <- year_range[5]
  
  df.out <- rbind(df.1, df.2, df.3, df.4, df.5)
  return(df.out)  
}

pearsons.ex <- data.table::rbindlist(lapply(pearsons.df.multi.split, expand.df))
head(pearsons.ex)
dist.dat.multi <- ddply(pearsons.ex, .(plot_year), summarise, median_corr = quantile(corr, na.rm=T)["50%"],  min_corr = quantile(corr, na.rm=T)["25%"], max_corr = quantile(corr, na.rm=T)["75%"])
head(dist.dat.multi)
tail(dist.dat.multi)
dist.dat.multi <- rbind(dist.dat.multi, dist.dat.multi[nrow(dist.dat.multi),])
dist.dat.multi$plot_year[length(dist.dat.multi$plot_year)] <- 2021

#and plot
Fig2Da <- ggplot(data=dist.dat.multi) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.3,.5,.1,1), "cm")) +
  geom_ribbon(aes(x=plot_year, ymin=min_corr, ymax=max_corr), alpha=.3) + ylab("multi-annual distribution\npairwise province\nPearson's correlation ") +
  geom_line(aes(x=plot_year, y=median_corr), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)

Fig2D <- cowplot::plot_grid(Fig2Da, Fig2Db, ncol=1, nrow=2, rel_heights = c(.25,1))

Fig2C2D <- cowplot::plot_grid(Fig2C, Fig2D, ncol=1, nrow=2, labels = c("C", "D"), label_size = 22, rel_heights = c(1,1.1), label_x = -0.01)

#and put them all together 
Fig2 <- cowplot::plot_grid(Fig2A2B, Fig2C2D, ncol=2, nrow = 1) + theme(panel.background = element_rect(fill="white"), plot.background = element_rect(fill="white"))


ggsave(file = paste0(homewd, "/final-figures/Fig2.png"),
       plot= Fig2,
       units="mm",  
       width=120, 
       height=110, 
       scale=3, 
       dpi=300)

ggsave(file = paste0(homewd, "/final-figures/Fig2.pdf"),
       plot= Fig2,
       units="mm",  
       width=120, 
       height=110, 
       scale=3, 
       dpi=300)


#and now for the supplementary figure S11!



###################################################################### 
###################################################################### 
#and Fig S11 with wavelet power for annual and multiannual cycles 
#(another way to support the results from above)

FigS11Ab <- ggplot(data=dat) +  coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time,  y=provname, fill=avg_wave_power_annual, color=avg_wave_power_annual)) + 
  scale_fill_viridis_c(trans="sqrt", option="inferno", name="Average power annual\ncycles from dengue incidence") +
  scale_color_viridis_c(trans="sqrt", option="inferno", name="Average power annual\ncycles from dengue incidence") +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(0,.5,.1,.7), "cm"),
                     strip.text = element_text(size=8), legend.position = "bottom",
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.spacing = unit(c(0), "cm"),
                     panel.border = element_rect(linewidth=0),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) 

FigS11Bb <- ggplot(data=dat) +  coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time,  y=provname, fill=avg_wave_power_multi, color=avg_wave_power_multi)) + 
  scale_fill_viridis_c(trans="sqrt", option="inferno", name="Average power multi-annual\ncycles from dengue incidence",  breaks=c(3,6,9,12)) +
  scale_color_viridis_c(trans="sqrt", option="inferno", name="Average power multi-annual\ncycles from dengue incidence", breaks=c(3,6,9,12)) +
  theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                     plot.margin = unit(c(0,.4,.1,.7), "cm"),
                     strip.text = element_text(size=8), legend.position = "bottom",
                     panel.spacing = unit(c(0), "cm"),
                     strip.text.y.left = element_text(angle=0, size=6),
                     panel.border = element_rect(linewidth=0),
                     panel.grid = element_blank(), axis.title = element_blank(),
                     axis.text = element_text(size=14)) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) 



#and get the monthly distribution of annual and biannual cases

dist.power.dat <- ddply(dat, .(time, year, biweek), summarise, median_annual = quantile(avg_wave_power_annual, na.rm=T)["50%"],  min_annual = quantile(avg_wave_power_annual, na.rm=T)["25%"], max_annual = quantile(avg_wave_power_annual, na.rm=T)["75%"], median_multi = quantile(avg_wave_power_multi, na.rm=T)["50%"], min_multi = quantile(avg_wave_power_multi, na.rm=T)["25%"], max_multi = quantile(avg_wave_power_multi, na.rm=T)["75%"])
head(dist.power.dat)
tail(dist.power.dat)



#and plot
FigS11Aa <- ggplot(data=dist.power.dat) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=9), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.5,.6,0,.8), "cm")) +
  geom_ribbon(aes(x=time, ymin=min_annual, ymax=max_annual), alpha=.3) + ylab("biweekly annual\nwavelet power\ndistribution") +
  geom_line(aes(x=time, y=median_annual), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)


FigS11Ba <- ggplot(data=dist.power.dat) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=9), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.5,.5,0,.9), "cm")) +
  geom_ribbon(aes(x=time, ymin=min_multi, ymax=max_multi), alpha=.3) + ylab("biweekly multi-\nannual wavelet\npower distribution") +
  geom_line(aes(x=time, y=median_multi), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)


#and construct
FigS11A <- cowplot::plot_grid(FigS11Aa, FigS11Ab, ncol=1, nrow=2, rel_heights = c(.25,1))
FigS11B <- cowplot::plot_grid(FigS11Ba, FigS11Bb, ncol=1, nrow=2, rel_heights = c(.25,1))



# Now and add in the proportion of other provinces with significant annual cross wavelet power per timestep
new.dat <-  read.csv(file = paste0(homewd, "/data/annual_multi_coherence_oct15.csv"), header=T, stringsAsFactors = F)
head(new.dat)

#first, test that annual is the same as above - this was done with a threshold p.value of 0.01
new.dat$is_sig_coherence <- 0
new.dat$is_sig_coherence[!is.na(new.dat$coherence)]<- 1
new.dat$is_sig_power <- 0
new.dat$is_sig_power[!is.na(new.dat$cross_wave_power)]<- 1

#and sum by time and province
new.sum <- ddply(new.dat, .(cycle_type, provname, time), summarise, avgCoherence=mean(coherence, na.rm=T), NposCoherence=sum(is_sig_coherence, na.rm=T), NtotCoherence=length(is_sig_coherence), avgPower = mean(cross_wave_power, na.rm=T), NPosPower = sum(is_sig_power, na.rm = T), NtotPower=length(is_sig_power))
head(new.sum)
new.sum$proportion_coherent <- new.sum$NposCoherence/new.sum$NtotCoherence
new.sum$proportion_coherent_power <- new.sum$NPosPower/new.sum$NtotPower
unique(new.sum$cycle_type)

FigS11Cb <- ggplot(data=subset(new.sum, cycle_type=="annual_cases")) +  coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time, y=provname, fill=proportion_coherent_power, color=proportion_coherent_power))+ 
  scale_fill_viridis_c( option="inferno", name="Proportion of provinces\nwith coherent annual\ncross wavelet power", limits=c(0,1), labels=c("0",".25",".5",".75","1")) +
  scale_color_viridis_c( option="inferno", name="Proportion of provinces\nwith coherent annual\ncross wavelet power", limits=c(0,1), labels=c("0",".25",".5",".75","1")) +
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

dist.dat.annual.multi <-  ddply(new.sum, .(cycle_type, time), summarise, median_prop = quantile(proportion_coherent_power)["50%"],  min_prop = quantile(proportion_coherent_power)["25%"], max_prop = quantile(proportion_coherent_power)["75%"])
unique(dist.dat.annual.multi$cycle_type)

FigS11Ca <- ggplot(data=subset(dist.dat.annual.multi, cycle_type=="annual_cases")) +theme_bw() + coord_cartesian(xlim=c(2001.9, 2021.1), ylim = c(0,1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.3,.5,0,.8), "cm")) +
  geom_ribbon(aes(x=time, ymin=min_prop, ymax=max_prop), alpha=.3) + ylab("biweekly distribution %\nprovinces w/coherent annual\ncross wavelet power") +
  geom_line(aes(x=time, y=median_prop), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)


FigS11C <- cowplot::plot_grid(FigS11Ca, FigS11Cb, ncol=1, nrow=2, rel_heights = c(.25,1))

#and multi

FigS11Db <- ggplot(data=subset(new.sum, cycle_type=="multi_cycles")) +  coord_cartesian(xlim=c(2001.9, 2021.1), expand = F)+
  facet_nested(provname~., scales = "free_y", space="free_y", switch = "y", strip = strip, labeller = label_wrap_gen(width=6)) +
  geom_tile(aes(x=time, y=provname, fill=proportion_coherent_power, color=proportion_coherent_power))+ 
  scale_fill_viridis_c( option="inferno", name="Proportion of provinces with\ncoherent multi-annual\ncross-wavelet power", limits=c(0,1), labels=c("0",".25",".5",".75","1")) +
  scale_color_viridis_c( option="inferno", name="Proportion of provinces with\ncoherent multi-annual\ncross-wavelet power", limits=c(0,1), labels=c("0",".25",".5",".75","1")) +
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


FigS11Da <- ggplot(data=subset(dist.dat.annual.multi, cycle_type=="multi_cycles")) +theme_bw() + 
  coord_cartesian(xlim=c(2001.9, 2021.1), ylim = c(0,1), expand = F)+
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=8), axis.text.y = element_text(size=8),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = unit(c(.3,.5,0,.8), "cm")) +
  geom_ribbon(aes(x=time, ymin=min_prop, ymax=max_prop), alpha=.3) + 
  ylab("biweekly distribution %\nprovinces w/coherent multi-\nannual cross-wavelet power") +
  geom_line(aes(x=time, y=median_prop), size=1) +
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1) +
  geom_hline(aes(yintercept=0), linetype=2)


FigS11D <- cowplot::plot_grid(FigS11Da, FigS11Db, ncol=1, nrow=2, rel_heights = c(.25,1))


FigS11AB <- cowplot::plot_grid(FigS11A, FigS11B, ncol=1, nrow=2, labels = c("A", "B"), label_size = 22, rel_heights = c(1,1.1), label_x = -0.01)
FigS11CD <- cowplot::plot_grid(FigS11C, FigS11D, ncol=1, nrow=2, labels = c("C", "D"), label_size = 22,rel_heights = c(1,1.1), label_x = -0.01)

#and put them all together 
FigS11 <- cowplot::plot_grid(FigS11AB, FigS11CD, ncol=2, nrow = 1)


ggsave(file = paste0(homewd, "/final-figures/FigS11.png"),
       plot= FigS11,
       units="mm",  
       width=120, 
       height=110, 
       scale=3, 
       dpi=300)

