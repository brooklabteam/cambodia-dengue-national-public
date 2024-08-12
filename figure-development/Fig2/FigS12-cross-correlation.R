rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(mgcv)
library(ggh4x)


homewd= "/Users/carabrook/Developer/cambodia-dengue-national-public"

#load data
pearsons.df <- read.csv(file = paste0(homewd, "/data/pearsons_correlations_provinces.csv"), header = T, stringsAsFactors = F)
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_dat_province.csv"), header = T, stringsAsFactors = F)
head(tsir.dat)
centroid.prov <- read.csv(file = paste0(homewd, "/data/centroid_provinces.csv"), header = T, stringsAsFactors = F)
head(centroid.prov)

#load and prep the climate data at the province level
temp.dat <- read.csv(file = paste0(homewd, "/data/biweek_temp.csv"), header = T, stringsAsFactors = F )
precip.dat <- read.csv(file = paste0(homewd, "/data/biweek_ppt.csv"), header = T, stringsAsFactors = F )
head(temp.dat) # year and biweek
head(precip.dat) # year and biweek


temp.dat = subset(temp.dat, provname!="Administrative unit not available")
precip.dat = subset(precip.dat, provname!="Administrative unit not available")


#setdiff(unique(temp.dat$provname), unique(beta.df$provname)) #"Mondul Kiri"    "Oddar Meanchey" "Ratanak Kiri"   "Siemreap"       "Tboung Khmum"  
temp.dat$provname[temp.dat$provname=="Siemreap"] <- "Siem Reap"
precip.dat$provname[precip.dat$provname=="Siemreap"] <- "Siem Reap"
temp.dat$provname[temp.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"
precip.dat$provname[precip.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"




pearsons.sum <- ddply(pearsons.df, .(provname, comp_prov), summarise, full_ts_corr=unique(full_ts_corr), dist_m = unique(dist_m), dist_km = unique(dist_km))
head(pearsons.sum)


colz = scales::hue_pal()(length(unique((tsir.dat$provname)))) #25
name.list <- sort(unique(as.character(tsir.dat$provname))) #alphabetical
name.list[name.list =="Otdar Meanchey"] <- "Oddar Meanchey"
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

p1 <- ggplot(pearsons.sum) + theme_bw() + ylab("full time series pearson's correlation between provinces") +
  xlab("distance between province centroids (km)") + coord_cartesian(ylim=c(0,1))+
  geom_point(aes(x=dist_km, y=full_ts_corr, color=comp_prov), size=3) +
  scale_color_manual(values=colz) +theme(panel.grid = element_blank(), legend.title = element_blank(),
                                         legend.position = "bottom") +facet_wrap(~provname)
p1

#then, get the average synchrony per province
pearsons.avg <- ddply(pearsons.sum, .(provname), summarise, mean_corr=mean(full_ts_corr ))

#then, look at avg synchrony per province per year
pearsons.avg.year <- ddply(pearsons.df, .(provname, year), summarise, mean_corr=mean(corr, na.rm=T))

pearsons.avg.year <- merge(pearsons.avg.year, centroid.prov, by="provname")

pearsons.avg.year <- arrange(pearsons.avg.year, latitude, year)
pearsons.avg.year$provname <- factor(pearsons.avg.year$provname, levels=unique(pearsons.avg.year$provname))


#and plot 
#vert.df <- cbind.data.frame(xint = c(2007, 2012,  2019))


vert.df <- cbind.data.frame(xint = c(2006.5, 2007.5, 2011.5, 2012.5,  2018.5, 2019.5))

# Only colour strips in x-direction
strip <- strip_themed(background_y = elem_list_rect(fill = colz))

# plot as heatmap - fig 2

#now add some GAM plots for FigS14

# add pop size, mean temp, mean precip, and elevation as covariates in this
# and calculate predictors of the correlation coefficient
centroid.add <- dplyr::select(centroid.prov, provname, mean_elevation_m)
pearsons.df <- merge(pearsons.df, centroid.add, by="provname")
tsir.dat$provname[tsir.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"

tsir.add <- ddply(tsir.dat, .(provname, year), summarise, pop=mean(pop))
#tsir.add <- ddply(tsir.dat, .(provname), summarise, pop=mean(pop))
pearsons.df <- merge(pearsons.df, tsir.add, by= c("provname", "year"))

#pearsons.df <- merge(pearsons.df, tsir.add, by= c("provname"))

#and temp and precip
mean.temp <- ddply(temp.dat, .(provname, year), summarise, mean_temp_C = mean(temp_C))

#mean.temp <- ddply(temp.dat, .(provname), summarise, mean_temp_C = mean(temp_C))
setdiff(unique(mean.temp$year), unique(pearsons.df$year))
setdiff(unique(pearsons.df$year), unique(mean.temp$year)) #no 2020 in climate data
setdiff(unique(mean.temp$provname), unique(pearsons.df$provname))

pearsons.df <- merge(pearsons.df, mean.temp, by=c("provname", "year"))

#pearsons.df <- merge(pearsons.df, mean.temp, by=c("provname"))

unique(precip.dat$provname)
precip.sum <- ddply(precip.dat, .(provname, year), summarise, sum_precip = sum(precip_mm))

#precip.sum <- ddply(precip.dat, .(provname), summarise, sum_precip = sum(precip_mm))

pearsons.df <- merge(pearsons.df, precip.sum, by=c("provname", "year"))

#pearsons.df <- merge(pearsons.df, precip.sum, by=c("provname"))

head(pearsons.df)
pearsons.df$log10_elevation_m <- log10(pearsons.df$mean_elevation_m)
pearsons.df$log10_temp_C <- log10(pearsons.df$mean_temp_C)
pearsons.df$log10_pop <- log10(pearsons.df$pop)
pearsons.df$log10_precip <- log10(pearsons.df$sum_precip)
pearsons.df$provname <- as.factor(pearsons.df$provname)
pearsons.df$year <- as.factor(pearsons.df$year)
#pearsons.df$year <- as.numeric(as.character(pearsons.df$year))

gam1 <- gam(corr ~ dist_m:provname +
                  s(log10_temp_C, bs="tp") +
                  s(log10_pop, bs="tp") +
                  s(log10_precip, bs="tp") +
                  s(year, bs="re") +
                  s(provname, bs="re"), #different y-intercept for each province
                  data = pearsons.df)

summary(gam1) 

#and write to table S6
TableS6A <- cbind.data.frame(province=c("intercept", unique(as.character(pearsons.df$provname))), slope = signif(summary(gam1)$p.table[,1],3), CI= paste0("[", signif((summary(gam1)$p.table[,1]-1.96*summary(gam1)$p.table[,2]),3),"-", signif((summary(gam1)$p.table[,1]+1.96*summary(gam1)$p.table[,2]),3), "]"), statistic = signif(summary(gam1)$p.table[,3],2), p_val = signif(summary(gam1)$p.table[,4],3))
TableS6A <- TableS6A[2:nrow(TableS6A),]
rownames(TableS6A) <- c()
write.csv(TableS6A, file = paste0(homewd,"/data/TableS6_part1.csv"), row.names = F)

#and the smoother table
TableS6B <- cbind.data.frame(smoother=c("log10_temp", "log10_pop", "log10_precip", "year", "province_random"), degrees_of_freedom = signif(summary(gam1)$s.table[,1],2), statistic = signif(summary(gam1)$s.table[,3],3), p_value = signif(summary(gam1)$s.table[,4],3))
rownames(TableS6B) <- c()
write.csv(TableS6B, file = paste0(homewd,"/data/TableS6_part2.csv"), row.names = F)

gam1b <- gam(corr ~ 
              s(dist_m, bs="tp")+
              s(log10_temp_C, bs="tp") +
              s(log10_pop, bs="tp") +
              s(log10_precip, bs="tp") +
              s(year, bs="re") +
              s(provname, bs="re"), #different y-intercept for each province
            data = pearsons.df)

summary(gam1b) 


#strong negative correlation with increasing distance and correlation, sig for most provinces
source(paste0(homewd, "/figure-development/Fig1/mollentze-streicker-2020-functions.R"))

year.df <- get_partial_effects(fit=gam1, var="year")
year.df <- year.df[[1]]
prov.df <- get_partial_effects(fit=gam1, var="provname") 
prov.df <- prov.df[[1]]

out.temp <- plot(gam1)[[1]]
temp.df <- cbind.data.frame(log10_temp_C=out.temp$x, y=out.temp$fit, ylower = (out.temp$fit - 1.96*out.temp$se),  yupper = (out.temp$fit + 1.96*out.temp$se))
temp.df$IsSig <- "Pos"
temp.df$temp_C <- 10^(temp.df$log10_temp_C)

out.precip <- plot(gam1)[[3]]
precip.df <- cbind.data.frame(log10_precip_mm=out.precip$x, y=out.precip$fit, ylower = (out.precip$fit - 1.96*out.precip$se),  yupper = (out.precip$fit + 1.96*out.precip$se))
precip.df$IsSig <- "Pos"
precip.df$precip_mm <- 10^(precip.df$log10_precip_mm)

out.pop <- plot(gam1)[[2]]
pop.df <- cbind.data.frame(log10_pop=out.pop$x, y=out.pop$fit, ylower = (out.pop$fit - 1.96*out.pop$se),  yupper = (out.pop$fit + 1.96*out.pop$se))
pop.df$IsSig <- "Pos"
pop.df$pop <- 10^(pop.df$log10_pop)


#and the fixed effects
out <- summary(gam1)
fixed.df <- cbind.data.frame(provname = c("intercept", sort(as.character(unique(pearsons.df$provname)))), y = out$p.coeff, ylower=(out$p.coeff -1.96*out$p.table[,2]), yupper=(out$p.coeff +1.96*out$p.table[,2]), pval =out$p.pv)
rownames(fixed.df) <- c()
fixed.df$IsSig <- "NotSig"
fixed.df$IsSig[fixed.df$y>0 & fixed.df$pval<0.01] <- "Pos"
fixed.df$IsSig[fixed.df$y<0 & fixed.df$pval<0.01] <- "Neg"

year.df$IsSig <- "NotSig"
year.df$IsSig[year.df$IsSignificant=="Yes" & year.df$y>0] <- "Pos"
year.df$IsSig[year.df$IsSignificant=="Yes" & year.df$y<0] <- "Neg"

#and plot
colz = c('Pos' = "red", 'Neg'="blue", 'NotSig' = "gray50" )

fixed.df = subset(fixed.df, provname!="intercept")
FigS12A <- ggplot(data=fixed.df) + theme_bw() + coord_flip() +
           theme(panel.grid = element_blank(), axis.title.y = element_blank()) +
           ylab(bquote(atop("fixed effect of province and geographic distance", "interaction on pairwise correlation coefficient,"~rho)))+
           geom_hline(aes(yintercept=0))+ 
           geom_point(aes(x=provname, y=y, color=IsSig), show.legend = F, size=3) +
           geom_linerange(aes(x=provname, ymin=ylower, ymax=yupper, color=IsSig), show.legend = F) +
           scale_color_manual(values=colz) 
# 
# FigS14B <- ggplot(data=year.df) + theme_bw() + coord_flip() +
#            theme(panel.grid = element_blank(), axis.title.y = element_blank()) +
#            ylab(bquote(atop("partial effect of year on", "pairwise correlation coefficient,"~rho)))+
#            geom_hline(aes(yintercept=0))+ 
#            geom_point(aes(x=year, y=y, color=IsSig), show.legend = F, size=3) +
#            geom_linerange(aes(x=year, ymin=ylower, ymax=yupper, color=IsSig), show.legend = F) +
#            scale_color_manual(values=colz) 


FigS12B <- ggplot(data=year.df) + theme_bw() + 
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(angle=300, vjust = 0)) +
  ylab(bquote(atop("partial effect of year on", "pairwise correlation coefficient,"~rho)))+
  geom_hline(aes(yintercept=0))+ 
  geom_point(aes(x=year, y=y, color=IsSig), show.legend = F, size=3) +
  geom_linerange(aes(x=year, ymin=ylower, ymax=yupper, color=IsSig), show.legend = F) +
  scale_color_manual(values=colz) 


FigS12C <- ggplot(data=temp.df) + theme_bw() + scale_x_log10() +
          theme(panel.grid = element_blank()) +
          ylab(bquote(atop("partial effect of temperature on", "pairwise correlation coefficient,"~rho)))+
          xlab(bquote("mean biweekly temperature ("^0~"C)"))+
          geom_hline(aes(yintercept=0))+ 
          geom_line(aes(x=temp_C, y=y, color=IsSig), show.legend = F, size=3) +
          geom_ribbon(aes(x=temp_C, ymin=ylower, ymax=yupper, fill=IsSig),alpha=.3,  show.legend = F) +
          scale_color_manual(values=colz) +scale_fill_manual(values=colz) 

FigS12D <- ggplot(data=precip.df) + theme_bw() + scale_x_log10() +
  theme(panel.grid = element_blank()) +
  ylab(bquote(atop("partial effect of precipitation on", "pairwise correlation coefficient,"~rho)))+
  xlab("mean total annual precipitation (mm)")+
  geom_hline(aes(yintercept=0))+ 
  geom_line(aes(x=precip_mm, y=y, color=IsSig), show.legend = F, size=3) +
  geom_ribbon(aes(x=precip_mm, ymin=ylower, ymax=yupper, fill=IsSig),alpha=.3,  show.legend = F) +
  scale_color_manual(values=colz) +scale_fill_manual(values=colz) 

    

FigS12E <- ggplot(data=pop.df) + theme_bw() + scale_x_log10() +
           theme(panel.grid = element_blank()) +
           ylab(bquote(atop("partial effect of province population on", "pairwise correlation coefficient,"~rho)))+
           xlab("mean population size")+
           geom_hline(aes(yintercept=0))+ 
           geom_line(aes(x=pop, y=y, color=IsSig), show.legend = F, size=3) +
           geom_ribbon(aes(x=pop, ymin=ylower, ymax=yupper, fill=IsSig),alpha=.3,  show.legend = F) +
           scale_color_manual(values=colz) +scale_fill_manual(values=colz)       
 

FigS12BCDE  <- cowplot::plot_grid(FigS12B, FigS12C, FigS12D, FigS12E, ncol=2, nrow = 2, align = "hv", labels = c("B", "C", "D", "E"), label_size = 22)

FigS12 <- cowplot::plot_grid(FigS12A, FigS12BCDE, ncol=2, nrow=1, labels = c("A", ""), rel_widths = c(1.05,2), label_size = 22, label_x = -.02)

ggsave(file = paste0(homewd, "/final-figures/FigS12.png"),
       plot= FigS12,
       units="mm",  
       width=110, 
       height=70, 
       scale=3, 
       dpi=300)
