rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(sf)
library(ggspatial)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national-public"

setwd(homewd)

# first plot the provinces by mean temp
# then attach the plots


cam = sf::st_read(paste0(homewd, "/data/province-shape/khm_admbnda_adm1_gov_20181004.shp"))
head(cam)
unique(cam$ADM1_EN)

# now load the climate data

climdat <- read.csv(file =  paste0(homewd, "/data/wavelet_input_dat_prov.csv"), header = T, stringsAsFactors = F)
head(climdat)

# get mean temp across the study period
clim.avg <- ddply(climdat, .(provname), summarise, temp_C=mean(temp_C))

#and sum precip by year and mean of that
precip.sum <- ddply(climdat, .(provname, year), summarise, sum_ppt_mm=sum(precip_mm))
precip.avg <- ddply(precip.sum,.(provname), summarise, mean_precip_annual_mm = mean(sum_ppt_mm))

#and join
clim.avg <- merge(clim.avg, precip.avg, by="provname")

setdiff(unique(clim.avg$provname), unique(cam$ADM1_EN))
setdiff(unique(cam$ADM1_EN), unique(clim.avg$provname)) #"Mondul Kiri"  "Ratanak Kiri" "Siemreap"     "Tboung Khmum"

cam$ADM1_EN[cam$ADM1_EN=="Siemreap"] <- "Siem Reap"

names(clim.avg)[names(clim.avg)=="provname"] <- "ADM1_EN"
cam_merge <- merge(cam,clim.avg, by = "ADM1_EN", all.x=T, sort=F)


FigTempLeg <- ggplot(cam_merge) + geom_sf(aes(fill=temp_C), size =.4) + 
  theme_bw() + coord_sf(xlim=c(100,140), ylim=c(9,15), expand = T) +
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.key.height = unit(c(.4),"cm"), 
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) +
  scale_fill_gradient(low = 'lightgoldenrod1',  high = 'firebrick', 
                      na.value = "grey50",
                      labels=scales::parse_format(),
                      name = bquote(atop("mean temperature,", "2002-2019 ("^0~"C)"))) 
  


LegTemp <- cowplot::get_legend(FigTempLeg)


FigTemp <- ggplot(cam_merge) +
  geom_sf(aes(fill=temp_C, color=ADM1_EN), linewidth=.6, show.legend = F) +
  theme_bw() +  
  coord_sf(xlim=c(96,114), ylim=c(4,22.6), expand = T) + 
  #scale_linewidth_continuous(range=c(10,10)) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        #legend.position = c(.12,.9),
        #legend.key.height = unit(c(.6),"cm"), 
        #legend.position = element_blank(),
        #strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) +
        scale_fill_gradient(low = 'lightgoldenrod1',  high = 'firebrick', 
                      na.value = "grey50",
                      labels=scales::parse_format()) +
  annotation_scale(pad_x = unit(16, "cm"), pad_y = unit(12, "cm"),
                   width_hint=.05, text_cex = .9)



# and load your TSIR projections
TSIR.ts <- read.csv(file = paste0(homewd,"/data/climate_projections_TSIR_by_prov.csv"), header = T, stringsAsFactors = F)



TSIR.ts$sim_type <- factor(TSIR.ts$sim_type, levels=c("TSIR-fit", "TSIR-prediction", "increased-S-standard-beta", "no-increase-S-climate-beta", "increased-S-climate-beta"))

colz = c('TSIR-fit'="#00B0F6",'TSIR-prediction' = "#E76BF3", 'increased-S-standard-beta' = "#A3A500", 'no-increase-S-climate-beta'="#00BF7D",'increased-S-climate-beta' = "#F8766D" )

#save legend from this facet plot
p1 <- ggplot(TSIR.ts) + facet_wrap(~provname, scales = "free_y") + ylab("reported cases") +
  geom_line(aes(x=time, y=cases), linetype=2) + scale_color_manual(values=colz, labels=scales::parse_format()) +
  scale_fill_manual(values=colz, labels=scales::parse_format()) + theme_bw() +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"), legend.text.align = 0,
        legend.title = element_blank(), legend.position = c(.91,.08), axis.title.x = element_blank()) +
  geom_ribbon(data = subset(TSIR.ts, time<2008),
              aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(TSIR.ts, time<2008),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) +
  geom_ribbon(data = subset(TSIR.ts,time>=2008 & time<2012  & sim_type=="TSIR-fit"),
              aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(TSIR.ts,time>=2008 & time<2012  & sim_type=="TSIR-fit"),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) +
  geom_ribbon(data = subset(TSIR.ts, time>=2012 & time<2013 & sim_type!="TSIR-fit"),
              aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(TSIR.ts, time>=2012 & time<2013 & sim_type!="TSIR-fit"),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) +
  geom_ribbon(data = subset(TSIR.ts, time>=2013 & time<2019  & sim_type=="TSIR-fit"),
              aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(TSIR.ts, time>=2013 & time<2019  & sim_type=="TSIR-fit"),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) +
  geom_ribbon(data = subset(TSIR.ts, time>=2019 & sim_type!="TSIR-fit"),
              aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2) +
  geom_line(data = subset(TSIR.ts, time>=2019 & sim_type!="TSIR-fit"),
            aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8) 


TSIR_legend <- cowplot::get_legend(p1)


#get a list of plots by province, then drop them on the map
make.plot.list <- function(df){
  
  df <- as.data.frame(df)
  
  df$sim_type <- factor(df$sim_type, levels=c("TSIR-fit", "TSIR-prediction", "increased-S-standard-beta", "no-increase-S-climate-beta", "increased-S-climate-beta"))
  
  colz = c('TSIR-fit'="#00B0F6",'TSIR-prediction' = "#E76BF3", 'increased-S-standard-beta' = "#A3A500", 'no-increase-S-climate-beta'="#00BF7D",'increased-S-climate-beta' = "#F8766D" )
  
  
  
  p1 <- ggplot(df) + facet_wrap(~provname, scales = "free_y") + ylab("reported cases") +
    geom_line(aes(x=time, y=cases), linetype=2) + scale_color_manual(values=colz, labels=scales::parse_format()) +
    scale_fill_manual(values=colz, labels=scales::parse_format()) + theme_bw() +
    theme(panel.grid = element_blank(), strip.background = element_rect(fill="white"), strip.text = element_text(size=11),
          legend.title = element_blank(), axis.title = element_blank()) +
    geom_ribbon(data = subset(df, time<2008),
                aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2, show.legend = F) +
    geom_line(data = subset(df, time<2008),
              aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8, show.legend = F) +
    geom_ribbon(data = subset(df,time>=2008 & time<2012  & sim_type=="TSIR-fit"),
                aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2, show.legend = F) +
    geom_line(data = subset(df,time>=2008 & time<2012  & sim_type=="TSIR-fit"),
              aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8, show.legend = F) +
    geom_ribbon(data = subset(df, time>=2012 & time<2013 & sim_type!="TSIR-fit"),
                aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2, show.legend = F) +
    geom_line(data = subset(df, time>=2012 & time<2013 & sim_type!="TSIR-fit"),
              aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8, show.legend = F) +
    geom_ribbon(data = subset(df, time>=2013 & time<2019  & sim_type=="TSIR-fit"),
                aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2, show.legend = F) +
    geom_line(data = subset(df, time>=2013 & time<2019  & sim_type=="TSIR-fit"),
              aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8, show.legend = F) +
    geom_ribbon(data = subset(df, time>=2019 & sim_type!="TSIR-fit"),
                aes(x=time, ymin=model_predicted_reported_cases_lci, ymax=model_predicted_reported_cases_uci, fill=sim_type), alpha=.2, show.legend = F) +
    geom_line(data = subset(df, time>=2019 & sim_type!="TSIR-fit"),
              aes(x=time, y=model_predicted_reported_cases, color=sim_type), size=.8, show.legend = F) 
  
  return(p1)
  
}

TSIR.ts$provname[TSIR.ts$provname=="Otdar Meanchey"] <- "Oddar Meanchey"

data.list = dlply(TSIR.ts,.(provname))

plot.list <- lapply(data.list, make.plot.list)

#now turn them into grobs and give them x,y max and min coordinates
grobfun = function(min_x, max_x, min_y, max_y, subplots) {
  p1 <- annotation_custom(ggplotGrob(subplots),
                    xmin = min_x, ymin = min_y,
                    xmax = max_x, ymax = max_y)
  return(p1)
}
join.dat <- function(pts.cam, plot.list){
  names_prov = names(plot.list)
  plot.sub <- list()
  
  for(i in 1:length(pts.cam$provname)){
    tmp = pts.cam$provname[i]
    plot.sub[[i]] <- plot.list[[tmp]]
  }
  
  out.grob.list <- mapply(grobfun, min_x = as.list(pts.cam$xmin), max_x=as.list(pts.cam$xmax), min_y=as.list(pts.cam$ymin), max_y = as.list(pts.cam$ymax), subplots=plot.sub, SIMPLIFY = F)
  
  return(out.grob.list)
}


#load the centroid data 
pts.cam <- read.csv(file = paste0(homewd, "/data/cambodia_province_lat_long.csv"), header = T, stringsAsFactors = F)
pts.cam <- pts.cam[complete.cases(pts.cam),]


out.list <- join.dat(pts.cam = pts.cam, plot.list = plot.list)

Fig1Temp <- FigTemp + 
            geom_segment(data=pts.cam, aes(x=longitude, xend=xcorner, y=latitude, yend=ycorner, colour=provname), size=.8, show.legend = F)


Fig1Final <- Fig1Temp + out.list[[1]] + out.list[[2]] + out.list[[3]] + out.list[[4]] + out.list[[5]] +
            out.list[[6]] + out.list[[7]] + out.list[[8]] + out.list[[9]] + out.list[[10]]  + out.list[[11]] +
            out.list[[12]]+ out.list[[13]] + out.list[[14]] + out.list[[15]] + out.list[[16]] +
            out.list[[17]]+ out.list[[18]] + out.list[[19]] + out.list[[20]] + out.list[[21]] +  out.list[[22]] +
            annotation_custom(TSIR_legend, xmin = 111, ymin = 20.8, xmax = 	114.8, ymax = 23.6) +
            annotation_custom(LegTemp, xmin = 95.2, ymin = 21.2, xmax = 	98.7, ymax = 23.2) 



# 
# ggsave(file = paste0(homewd, "/final-figures/Fig1temp.png"),
#        plot= Fig1Temp,
#        units="mm",  
#        width=200, 
#        height=100, 
#        scale=3, 
#        dpi=500)



ggsave(file = paste0(homewd, "/final-figures/Fig1.png"),
       plot= Fig1Final,
       units="mm",  
       width=110, 
       height=100, 
       scale=3, 
       dpi=500)


ggsave(file = paste0(homewd, "/final-figures/Fig1.pdf"),
       plot= Fig1Final,
       units="mm",  
       width=110, 
       height=100, 
       scale=3, 
       dpi=500)

#and save

############################################################################
############################################################################

# Below: versions of the country by color block and precip

FigProv <- ggplot(cam_merge) + geom_sf(aes(fill=ADM1_EN), size =.4) + 
  theme_bw() +
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) 


FigPrecip <- ggplot(cam_merge) + geom_sf(aes(fill=mean_precip_annual_mm), size =.4) + 
  theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom",
        axis.title = element_blank(),
        plot.margin = unit(c(.3,.2,0,3),"lines")) +
  scale_fill_gradient(low = 'lightblue',  high = 'navy', 
                      na.value = "grey50",
                      labels=scales::parse_format(),
                      name="mean annual precipitation, 2002-2019 (mm)") 