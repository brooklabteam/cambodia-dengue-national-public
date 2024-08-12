rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(WaveletComp)
library(mgcv)
library(reshape2)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national-public"
setwd(homewd)

get.mean.period <- function(col.df){
  
  #bind
  dat.period <- cbind.data.frame(period=col.df$period_length, power=col.df$value)
  
  mean_period = weighted.mean(x=col.df$period_length, w=col.df$value)
  
  return(mean_period)
}
get.max.period <- function(col.df){
  
  #bind
  #dat.period <- cbind.data.frame(period=col.df$period_length, power=col.df$value)
  
  max_period = col.df$period_length[col.df$value==max(col.df$value)]
  
  if(length(max_period)>1){
    max_period = min(max_period) #if more than one, take the shortest
  }
  
  return(max_period)
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
  max.period.list <- lapply(col.split, get.max.period)
  
  nat.dat$multi_period_mean <-  c(unlist(mean.period.list))
  nat.dat$multi_period_max <-  c(unlist(max.period.list))
  
  return(nat.dat)
}



dat <- read.csv(file = paste0(homewd, "/data/synchrony_data_aug12.csv"), header=T, stringsAsFactors = F)
head(dat)
names(dat)


#dat$month_date <- as.Date(dat$month_date)#, format = "%m/%d/%y")


#and plot
dat <- arrange(dat, desc(pop), time)

dat$provname <- factor(dat$provname, levels=unique(dat$provname))
dat$epiyear <- "no"
#dat$year <- year(dat$month_date)
dat$epiyear[dat$year==2007| dat$year==2012| dat$year==2019] <- "yes"
dat$epiyear <- as.factor(dat$epiyear)

# epicolz = c('no'=NA, 'yes'="red3") 
# 
# colorz<- c('-1' = "navy", '0'="blue", '1'="cyan", '2'="green", '3'="yellow", '4'="orange", '5'="red")
  #

library(patchwork)
library(RColorBrewer)
display.brewer.all()


# define jet colormap
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))


p1 <- ggplot(data=dat) + theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        panel.background = element_rect(fill="gray30"),
        axis.text = element_text(size=14), legend.position = "bottom") +
  scale_fill_gradientn(colors = jet.colors(7), name="Reconstructed\nannual cycles", limits=c(-2,4)) +
  scale_color_gradientn(colors = jet.colors(7), name="Reconstructed\nannual cycles", limits=c(-2,4)) +
  #scale_fill_viridis_c(option="turbo", name="Reconstructed\nannual cycles", limits=c(-2,4)) + 
  #scale_color_viridis_c(option="turbo", name="Reconstructed\nannual cycles", limits=c(-2,4))  + 
  geom_tile(aes(x=time, y=provname, fill=reconstructed_annual_period, color=reconstructed_annual_period))#, color=epiyear))


p2 <- ggplot(data=dat) + theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        panel.background = element_rect(fill="gray30"),
        axis.text = element_text(size=14), legend.position = "bottom") +
  scale_fill_gradientn(colors = jet.colors(7), name="Reconstructed\nmulti-annual cycles", limits=c(-2,4)) +
  scale_color_gradientn(colors = jet.colors(7), name="Reconstructed\nmulti-annual cycles", limits=c(-2,4)) +
  #scale_fill_viridis_c(option="turbo", name="Reconstructed\nmulti-annual cycles", limits=c(-2,4)) + 
  #scale_color_viridis_c(option="turbo", name="Reconstructed\nmulti-annual cycles", limits=c(-2,4))  +
  geom_tile(aes(x=time, y=provname, fill=reconstructed_multi_period), color=NA)#, color=epiyear))


p1 | p2


####### period length
nat.dat <- ddply(dat, .(time), summarise, cases=sum(cases))
head(nat.dat)

#join with pop data
popdat <- read.csv(file = paste0(homewd, "/data/cambodia_pop_dat.csv"), header=T, stringsAsFactors = F)
head(popdat)

nat.dat <- get.national.period(nat.dat = nat.dat,popdat = popdat)


p3 <- ggplot(dat) + theme_bw() + ylab("mean period duration,\ndengue incidence per 1000 ppl") +
  theme(panel.grid = element_blank(), legend.position = "bottom", axis.title.x = element_blank(), axis.title.y = element_text(size=16),
        axis.text = element_text(size=14)) +
  geom_line(aes(x=time,y=multi_period_mean, color=provname)) +
  geom_line(data=nat.dat, aes(x=time, y= multi_period_mean), linewidth=1.5)

p3

# 
# p3b <- ggplot(dat) + theme_bw() + ylab("dominant period duration,\ndengue incidence peri 1000 ppl") +
#   theme(panel.grid = element_blank(), legend.position = "bottom", axis.title.x = element_blank(), axis.title.y = element_text(size=16),
#         axis.text = element_text(size=14)) +
#   geom_line(aes(x=time,y=multi_period_max, color=provname)) +
#   geom_line(data=nat.dat, aes(x=time, y= multi_period_max), linewidth=1.5)
# 
# p3b


p3c <- ggplot(dat) + theme_bw() + ylab("mean period duration,\ndengue incidence per 1000 ppl") +
  theme(panel.grid = element_blank(), legend.position = "bottom", axis.title.x = element_blank(), axis.title.y = element_text(size=16),
        axis.text = element_text(size=14)) +
  geom_line(aes(x=time,y=multi_period_mean, color=provname)) +
  #geom_line(data=nat.dat, aes(x=month_date, y= multi_period), linewidth=1.5) + 
  facet_wrap(~provname, ncol=4)

p3c


library(lme4)
library(lmerTest)
dat$provname
m1 <- lmer(multi_period_mean~time:provname + (1|provname), data=dat)
summary(m1)   
m2 <- lmer(multi_period_mean~year:provname + (1|provname), data=dat)
summary(m2)   


#less here - power
p4 <- ggplot(data=dat) + theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        panel.background = element_rect(fill="gray30"),
        axis.text = element_text(size=14), legend.position = "bottom") +
  scale_fill_viridis_c(trans="sqrt", option="inferno", name="Average power\nannual cycles") +# scale_color_manual(values=epicolz) +
  scale_color_viridis_c(trans="sqrt", option="inferno", name="Average power\nannual cycles") +# scale_color_manual(values=epicolz) +
  geom_tile(aes(x=time, y=provname, fill=avg_wave_power_annual))#, color=epiyear))

p5 <- ggplot(data=dat) + theme_bw() +
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(),
        panel.background = element_rect(fill="gray30"),
        panel.spacing.x = unit(0, "mm"),
        axis.text = element_text(size=14), legend.position = "bottom") +
  scale_fill_viridis_c(trans="sqrt", option="inferno", name="Average power\nmultiannual cycles") +# scale_color_manual(values=epicolz) +
  scale_color_viridis_c(trans="sqrt", option="inferno", name="Average power\nmultiannual cycles") +# scale_color_manual(values=epicolz) +
  geom_tile(aes(x=time, y=provname, fill=avg_wave_power_multi))#, color=epiyear))


p4 | p5

#and coherence
vert.df <- cbind.data.frame(xint = c(2007, 2008, 2012, 2013, 2019,2020))
 p6 <- ggplot(data=dat) + theme_bw() +
   theme(panel.grid = element_blank(), axis.title = element_blank(),
         panel.background = element_rect(fill="gray30"),
         axis.text = element_text(size=14),legend.position = "bottom") +
   scale_fill_viridis_c( option="inferno", name="Proportion of provinces\nwith coherent annual cycles") +# scale_color_manual(values=epicolz) +
   scale_color_viridis_c( option="inferno", name="Proportion of provinces\nwith coherent annual cycles") +# scale_color_manual(values=epicolz) 
   geom_tile(aes(x=time, y=provname, fill=proportion_coherence_annual))+ 
   geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1)
p6
 # 
# # p7 <- ggplot(data=dat) + theme_bw() +
#   theme(panel.grid = element_blank(), axis.title = element_blank(),
#         panel.background = element_rect(fill="gray30"),
#         axis.text = element_text(size=14),legend.position = "bottom") +
#   scale_fill_viridis_c( option="inferno", name="Proportion of provinces\nwith coherent multi-annual cycles") +# scale_color_manual(values=epicolz) +
#   scale_color_viridis_c( option="inferno", name="Proportion of provinces\nwith coherent multi-annual cycles") +# scale_color_manual(values=epicolz) 
#   geom_tile(aes(x=time, y=provname, fill=proportion_coherence_multi))#, color=epiyear))

#p6# | p7



# #temp
# p10 <- ggplot(data=dat) + theme_bw() +
#   theme(panel.grid = element_blank(), axis.title = element_blank(),
#         panel.background = element_rect(fill="gray30"),
#         axis.text = element_text(size=14),legend.position = "bottom") +
#   scale_fill_viridis_c( option="inferno", name="Coherence (multi):\n dengueIR and\ntemperature") +# scale_color_manual(values=epicolz) +
#   scale_color_viridis_c( option="inferno", name="Coherence (multi):\n dengue IR and\ntemperature") +# scale_color_manual(values=epicolz) 
#   geom_tile(aes(x=time, y=provname, fill=avg_wave_coherency_temp_multi))#, color=epiyear))
# 

#precip

# p11 <- ggplot(data=dat) + theme_bw() +
#   theme(panel.grid = element_blank(), axis.title = element_blank(),
#         panel.background = element_rect(fill="gray30"),
#         axis.text = element_text(size=14),legend.position = "bottom") +
#   scale_fill_viridis_c( option="inferno", name="Coherence (multi):\ndengue IR and\nprecipitation") +# scale_color_manual(values=epicolz) +
#   scale_color_viridis_c( option="inferno", name="Coherence (multi):\ndengue IR and\nprecipitation") +# scale_color_manual(values=epicolz) 
#   geom_tile(aes(x=time, y=provname, fill=avg_wave_coherency_precip_multi))#, color=epiyear))

#p10 | p11


#wave power
#temp
 #p12 <- ggplot(data=dat) + theme_bw() +
#   theme(panel.grid = element_blank(), axis.title = element_blank(),
#         panel.background = element_rect(fill="gray30"),
#         axis.text = element_text(size=14),legend.position = "bottom") +
#   scale_fill_viridis_c( option="inferno", name="Average cross-wavelet\npower (multi): dengue IR\nand temperature") +# scale_color_manual(values=epicolz) +
#   scale_color_viridis_c( option="inferno", name="Average cross-wavelet\npower (multi): dengue IR\nand temperature") +# scale_color_manual(values=epicolz) 
#   geom_tile(aes(x=month_date, y=provname, fill=avg_cross_power_temp_multi))#, color=epiyear))
# 
# 
# #precip
# 
# p13 <- ggplot(data=dat) + theme_bw() +
#   theme(panel.grid = element_blank(), axis.title = element_blank(),
#         panel.background = element_rect(fill="gray30"),
#         axis.text = element_text(size=14),legend.position = "bottom") +
#   scale_fill_viridis_c( option="inferno", name="Average cross-wavelet\npower (multi): dengue IR\nand precipitation") +# scale_color_manual(values=epicolz) +
#   scale_color_viridis_c( option="inferno", name="Average cross-wavelet\npower (multi): dengue IR\nand precipitation") +# scale_color_manual(values=epicolz) 
#   geom_tile(aes(x=month_date, y=provname, fill=avg_cross_power_precip_multi))#, color=epiyear))
# 
# p12 | p13
# 

#and the annual
#temp
 p14 <- ggplot(data=dat) + theme_bw() +
   theme(panel.grid = element_blank(), axis.title = element_blank(),
         panel.background = element_rect(fill="gray30"),
         axis.text = element_text(size=14),legend.position = "bottom") +
   scale_fill_viridis_c( option="inferno", name="Coherence (annual):\ndengue IR and\ntemperature") +# scale_color_manual(values=epicolz) +
   scale_color_viridis_c( option="inferno", name="Coherence (annual):\ndengue IR and\ntemperature") +# scale_color_manual(values=epicolz) 
   geom_tile(aes(x=time, y=provname, fill=avg_wave_coherency_temp))+
   geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1)
# 
# 
# #precip
# 
 p15 <- ggplot(data=dat) + theme_bw() +
   theme(panel.grid = element_blank(), axis.title = element_blank(),
         panel.background = element_rect(fill="gray30"),
         axis.text = element_text(size=14),legend.position = "bottom") +
   scale_fill_viridis_c( option="inferno", name="Coherence (annual):\ndengue IR and\nprecipitation") +# scale_color_manual(values=epicolz) +
   scale_color_viridis_c( option="inferno", name="Coherence (annual):\ndengue IR and\nprecipitation") +# scale_color_manual(values=epicolz) 
   geom_tile(aes(x=time, y=provname, fill=avg_wave_coherency_precip))+
   geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1)
# 
 p14 | p15


#wave power
#temp
p16 <- ggplot(data=dat) + theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        panel.background = element_rect(fill="gray30"),
        axis.text = element_text(size=14),legend.position = "bottom") +
  scale_fill_viridis_c( option="inferno", name="Average cross-wavelet\npower (annual): dengue IR\nand temperature", trans="sqrt") +# scale_color_manual(values=epicolz) +
  scale_color_viridis_c( option="inferno", name="Average cross-wavelet\npower (annual): dengue IR\nand temperature", trans="sqrt") +# scale_color_manual(values=epicolz) 
  geom_tile(aes(x=time, y=provname, fill=avg_cross_power_temp))+
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1)


#precip

p17 <- ggplot(data=dat) + theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        panel.background = element_rect(fill="gray30"),
        axis.text = element_text(size=14),legend.position = "bottom") +
  scale_fill_viridis_c( option="inferno", name="Average cross-wavelet\npower (annual): dengue IR\nand precipitation", trans="sqrt") +# scale_color_manual(values=epicolz) +
  scale_color_viridis_c( option="inferno", name="Average cross-wavelet\npower (annual): dengue IR\nand precipitation", trans="sqrt") +# scale_color_manual(values=epicolz) 
  geom_tile(aes(x=time, y=provname, fill=avg_cross_power_precip))+
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1)

p16 | p17



#annual cycles are synchronized across provinces and with temperature and precip in epidemic years

p4 | p16  | p17

#and coherence with oni

oni.dat <- read.csv(file = paste0(homewd, "/data/oni_coherence.csv"), header=T, stringsAsFactors = F)
head(oni.dat)

#and plot
p18 <- ggplot(data=oni.dat) + theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        panel.background = element_rect(fill="gray30"),
        axis.text = element_text(size=14),legend.position = "bottom") +
  scale_fill_viridis_c( option="inferno", name="Average cross-wavelet\npower (multiannual):\ndengue IR and ONI", trans="sqrt") +# scale_color_manual(values=epicolz) +
  scale_color_viridis_c( option="inferno", name="Average cross-wavelet\npower (multiannual):\ndengue IR and ONI", trans="sqrt") +# scale_color_manual(values=epicolz) 
  geom_tile(aes(x=time, y=provname, fill=avg_cross_power_oni))+
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1)

p18

#and coherence

p19 <- ggplot(data=oni.dat) + theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        panel.background = element_rect(fill="gray30"),
        axis.text = element_text(size=14),legend.position = "bottom") +
  scale_fill_viridis_c( option="inferno", name="Average wavelet\ncoherence (multiannual): dengue IR\nand ONI", trans="sqrt") +# scale_color_manual(values=epicolz) +
  scale_color_viridis_c( option="inferno", name="Average wavelet\ncoherence (multiannual): dengue IR\nand ONI", trans="sqrt") +# scale_color_manual(values=epicolz) 
  geom_tile(aes(x=time, y=provname, fill=avg_wave_coherence_oni))+
  geom_vline(data=vert.df, aes(xintercept=xint), color="red", size=1)

p19


