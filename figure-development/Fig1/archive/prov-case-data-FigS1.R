rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(mgcv)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)

#plot time series and time trend by province

dat <- read.csv(file = paste0(homewd, "/data/DENV-Dist-Diag.csv"), header=T, stringsAsFactors = F)
head(dat)
#change Oddar Meanchey
dat$provname[dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"
# hist(dat$age)
# sum.dat <- ddply(dat, .(year, age), summarise, N=length(age))
# subset(sum.dat, age>50)
# 
# # p1 <- ggplot(data=dat) + geom_violin(aes(x=year, y=age, group=year), 
# #             draw_quantiles = c(0,.25,.5,.75)) +
# #   geom_jitter(aes(x=year, y=age), width=.1, size=.1, alpha=.2) + facet_wrap(~provname)
# #   
# # p1
# 
# #and the max age of dengue infection
# age.max <- ddply(dat, .(year), summarise, max_age = max(age))
# ggplot(age.max) + geom_point(aes(x=year, y = max_age)) + geom_line(aes(x=year, y=max_age))

popdat <- read.csv(file = paste0(homewd, "/data/cambodia_pop_dat.csv"), header=T, stringsAsFactors = F)

#sum by epidemic week, type, province
dat.prov <- ddply(dat, .(year, epimonth, provname, diagnostic), summarise, cases = sum(case))
#dat.nat <- ddply(dat, .(year, epimonth, diagnostic), summarise, cases = sum(case))
#dat.nat$epimonth <- as.Date(dat.nat$epimonth, format = "%m/%d/%y")
dat.prov$epimonth <- as.Date(dat.prov$epimonth, format = "%m/%d/%y")


#add pop
popmerge <- dplyr::select(popdat, year, pop, provname)
dat.prov <- merge(dat.prov, popdat, by=c("year", "provname"), all.x = T)
head(dat.prov)

library(lubridate)

dat.prov$month <- month(dat.prov$epimonth)


#need to eliminate one province with short time series
dat.prov = subset(dat.prov, provname!="Tboung Khmum")
dat.prov$provname <- as.factor(dat.prov$provname)
dat.prov$cases_per_1000 <- (dat.prov$cases/dat.prov$pop)*1000


m1 <- gam(cases_per_1000 ~ 
            year:provname + #slope specific by province
            s(month, k=7, bs="cc") + #controls for internal annual cycles
            s(provname, bs="re"), #y-intercept specific by province too
           dat = subset(dat.prov, diagnostic=="df"))
summary(m1) #all increasing


out=summary(m1)

prov.df <- cbind.data.frame(provname = unique(dat.prov$provname), 
                            slope= out$p.coeff[2:length(out$p.coeff)],
                            se = out$se[2:length(out$p.coeff)],
                            p_val = out$p.pv[2:length(out$p.coeff)])

prov.df$slope_uci <- prov.df$slope + 1.96*prov.df$se
prov.df$slope_lci <- prov.df$slope - 1.96*prov.df$se
prov.df$sig <- "sig"
prov.df$sig[prov.df$p_val>0.05] <- "not_sig"
rownames(prov.df) <- c()
prov.df$diagnostic = "df"

#and predict
predict.dat <- cbind.data.frame(provname=rep(unique(dat.prov$provname), each = length(unique(dat.prov$year))), year=rep(unique(dat.prov$year), length(unique(dat.prov$provname))))
predict.dat$month <- 1

predict.dat$predict_cases <- predict.gam(m1,newdata = predict.dat,  exclude = "s(month)", type = "response")
predict.dat$predict_cases_lci  <- predict.dat$predict_cases - 1.96*predict.gam(m1, newdata = predict.dat,  exclude = "s(month)", type = "response", se.fit = T)$se
predict.dat$predict_cases_uci  <- predict.dat$predict_cases + 1.96*predict.gam(m1, newdata = predict.dat,  exclude = "s(month)", type = "response", se.fit = T)$se


predict.dat$diagnostic = "df"



#and repeat for dhf and dss

m2 <- gam(cases_per_1000 ~ year:provname + 
            #s(year, by=provname, bs="tp", k=3) + #slope specific by province
            s(month, k=7, bs="cc") + #controls for internal annual cycles
            s(provname, bs="re"), #y-intercept specific by province too
          dat = subset(dat.prov, diagnostic=="dhf"))
summary(m2) #all sig. all positive


out=summary(m2)

prov.df.2 <- cbind.data.frame(provname = unique(dat.prov$provname), 
                            slope= out$p.coeff[2:length(out$p.coeff)],
                            se = out$se[2:length(out$p.coeff)],
                            p_val = out$p.pv[2:length(out$p.coeff)])

prov.df.2$slope_uci <- prov.df.2$slope + 1.96*prov.df.2$se
prov.df.2$slope_lci <- prov.df.2$slope - 1.96*prov.df.2$se
prov.df.2$sig <- "sig"
prov.df.2$sig[prov.df.2$p_val>0.05] <- "not_sig"
rownames(prov.df.2) <- c()

prov.df.2$diagnostic <- "dhf"

#and predict
predict.dat.2 <- cbind.data.frame(provname=rep(unique(dat.prov$provname), each = length(unique(dat.prov$year))), year=rep(unique(dat.prov$year), length(unique(dat.prov$provname))))
predict.dat.2$month <- 1

predict.dat.2$predict_cases <- predict.gam(m2,newdata = predict.dat.2,  exclude = "s(month)", type = "response")
predict.dat.2$predict_cases_lci  <- predict.dat.2$predict_cases - 1.96*predict.gam(m2, newdata = predict.dat.2,  exclude = "s(month)", type = "response", se.fit = T)$se
predict.dat.2$predict_cases_uci  <- predict.dat.2$predict_cases + 1.96*predict.gam(m2, newdata = predict.dat.2,  exclude = "s(month)", type = "response", se.fit = T)$se


predict.dat.2$diagnostic = "dhf"


#and dss

m3 <- gam(cases_per_1000 ~ year:provname + 
            #s(year, by=provname, bs="tp", k=3) + #slope specific by province
            s(month, k=7, bs="cc") + #controls for internal annual cycles
            s(provname, bs="re"), #y-intercept specific by province too
          dat = subset(dat.prov, diagnostic=="dss"))
summary(m3) #all but one significant, all neg


out=summary(m3)

prov.df.3 <- cbind.data.frame(provname = unique(dat.prov$provname), 
                              slope= out$p.coeff[2:length(out$p.coeff)],
                              se = out$se[2:length(out$p.coeff)],
                              p_val = out$p.pv[2:length(out$p.coeff)])

prov.df.3$slope_uci <- prov.df.3$slope + 1.96*prov.df.3$se
prov.df.3$slope_lci <- prov.df.3$slope - 1.96*prov.df.3$se
prov.df.3$sig <- "sig"
prov.df.3$sig[prov.df.3$p_val>0.05] <- "not_sig"
rownames(prov.df.3) <- c()

prov.df.3$diagnostic <- "dss"

#and predict
predict.dat.3 <- cbind.data.frame(provname=rep(unique(dat.prov$provname), each = length(unique(dat.prov$year))), year=rep(unique(dat.prov$year), length(unique(dat.prov$provname))))
predict.dat.3$month <- 1

predict.dat.3$predict_cases <- predict.gam(m3,newdata = predict.dat.3,  exclude = "s(month)", type = "response")
predict.dat.3$predict_cases_lci  <- predict.dat.3$predict_cases - 1.96*predict.gam(m3, newdata = predict.dat.3,  exclude = "s(month)", type = "response", se.fit = T)$se
predict.dat.3$predict_cases_uci  <- predict.dat.3$predict_cases + 1.96*predict.gam(m3, newdata = predict.dat.3,  exclude = "s(month)", type = "response", se.fit = T)$se

predict.dat.3$diagnostic = "dss"


#and combine
prov.df <- rbind(prov.df,prov.df.2, prov.df.3)
predict.dat <- rbind(predict.dat, predict.dat.2, predict.dat.3)

predict.dat$epiyear <- as.Date(paste0(predict.dat$year, "-01-01"))

colz <- c("df" = "navy", dhf="maroon", dss="darkseagreen4")
#and supplemental figure
FigS1A <- ggplot(predict.dat) + scale_color_manual(values=colz) +scale_fill_manual(values=colz) +
          geom_line(data= subset(dat.prov, diagnostic=="df"), aes(x=epimonth, y=cases_per_1000), color="gray70") +
          geom_line(aes(x=epiyear, y=predict_cases, color=diagnostic)) + theme_bw() +
          geom_ribbon(aes(x=epiyear, ymin=predict_cases_lci, ymax=predict_cases_uci, fill=diagnostic), alpha=.3) + 
          facet_wrap(~provname, ncol=4) + ylab("monthly predicted annual cases per 1000 ppl") +
          theme(axis.title.x = element_blank(), panel.grid = element_blank(), 
                strip.background = element_rect(fill="white"),
                plot.margin = unit(c(.3,.2,.2,.2),"lines"),
          legend.title = element_blank(), legend.position = "bottom",
          axis.text = element_text(size=14), axis.title.y = element_text(size=16),
          legend.text = element_text(size=12)) + coord_cartesian(ylim=c(0,.5))
      

# ggsave(paste0(homewd, "/fig-new/FigS1.png"),
#        plot = FigS1,
#        bg='white',
#        width=90, 
#        height=80, 
#        scale=3.2, 
#        dpi=300,
#        units = "mm")



### and then FigS1B is the map with the slopes by province
library(sf)
cam = sf::st_read(paste0(homewd, "/data/province-shape/khm_admbnda_adm1_gov_20181004.shp"))
head(cam)
unique(cam$ADM1_EN)

prov.df$provname <- as.character(prov.df$provname)

setdiff(unique(prov.df$provname), unique(cam$ADM1_EN))
setdiff(unique(cam$ADM1_EN), unique(prov.df$provname)) #"Tboung Khmum"

cam$ADM1_EN[cam$ADM1_EN=="Siemreap"] <- "Siem Reap"
#cam$ADM1_EN[cam$ADM1_EN=="Oddar Meanchey"] <- "Otdar Meanchey"

#and add the blanks
head(prov.df)
prov.add <- subset(prov.df, provname=="Kandal")
prov.add$provname <- "Tboung Khmum"
prov.add$sig <- "not_sig"
prov.add$slope <- prov.add$se <- prov.add$p_val <- prov.add$slope_lci <- prov.add$slope_uci <- NA

prov.df <- rbind(prov.df, prov.add)
#and merge with the slopes
#dat.merge <- dplyr::select(case.dat, y, IsSignificant, provname, diagnostic)

names(prov.df)[names(prov.df)=="provname"] <- "ADM1_EN"
cam_merge <- merge(cam, prov.df, by = "ADM1_EN", all.x=T, sort=F)

cam_sub = subset(cam_merge, is.na(slope))
head(cam_merge)

unique(prov.df$sig)
colorz <- c('sig' = "black", 'not_sig'="grey80")

#library(patchwork)

FigS1B <- ggplot(cam_merge) + geom_sf(aes(fill=slope, color=sig), size =.4) + 
  theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(fill="white"),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(.3,.2,0,3),"lines")) +
  facet_grid(~diagnostic) +
  scale_fill_gradient2(low = 'navy', mid = 'white', high = 'firebrick', 
                       na.value = "grey50",
                       name="slope of cases per 1000\npeople through time") +
  scale_color_manual(values = colorz) +
  guides(color="none")



FigS1 <- cowplot::plot_grid(FigS1A, FigS1B, ncol=1, nrow=2, rel_heights  = c(1.8,1), labels = c("A", "B"), label_size = 22, align = "hv", axis="tb", label_x = -.005, label_y = 1.005)


ggsave(paste0(homewd, "/final-figures/FigS1.png"),
       plot = FigS1,
       bg='white',
       width=110, 
       height=110, 
       scale=3, 
       dpi=300,
       units = "mm")
