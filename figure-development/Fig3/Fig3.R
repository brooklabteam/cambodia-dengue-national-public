rm(list=ls())

library(plyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(mgcv)
library(reshape2)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national-public"
setwd(homewd)

# plot age distribution by province,
# mean age by province 
# and FOI by province, as multipanels

#then, decide how to arrange all together

dat = read.csv(file = paste0(homewd, "/data/DENV-Dist-Diag.csv"), header=T, stringsAsFactors = F)
unique(dat$dianostic) #df, dhf, dss
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
head(dat)
subset(dat, age >70)

#now attach biwks and doy (2 biweeks per month)
wks = 1:52
biwks <- rep(1:26, each=2)

#and attach a doy that corresponds to each biweek
doy.biwk <- seq(1,365,length.out=26)

wk.biwk <- cbind.data.frame(week_report=wks, biwk=biwks, doy=doy.biwk)


sort(unique(dat$week_report))
unique(dat$provname)

dat <- merge(dat, wk.biwk, by="week_report")
dat <- arrange(dat, procode, date)

#okay... plot ages by biweek by year
head(dat) 
dat$time <- (dat$doy/365) + dat$year

dat$provname[dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"
colz = scales::hue_pal()((length(unique(dat$provname))))
names(colz) <- unique(dat$provname)



p1 <- ggplot(dat) + geom_jitter(aes(x=year, y=age), width=.1, size=.1, alpha=.3) + 
      facet_wrap(~provname) + theme_bw() + 
      theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.background = element_rect(fill="white")) +
      geom_violin(aes(x=year,y=age, group=year, color = provname),  
                  draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) 

# add in the mean age of infection here too - as a linear regression
dat$provname <- as.factor(dat$provname)
m1 <- gam(age~year:provname +
              s(provname, bs="re"), 
              data=dat)
summary(m1)

#and save this for Table S7
TableS7 <- cbind.data.frame(province=c("intercept", unique(as.character(dat$provname))), slope = signif(summary(m1)$p.table[,1],3), CI= paste0("[", signif((summary(m1)$p.table[,1]-1.96*summary(m1)$p.table[,2]),3),"-", signif((summary(m1)$p.table[,1]+1.96*summary(m1)$p.table[,2]),3), "]"), statistic = signif(summary(m1)$p.table[,3],2), p_val = signif(summary(m1)$p.table[,4],3))
TableS7  <- TableS7 [2:nrow(TableS7),]
rownames(TableS7 ) <- c()
write.csv(TableS7 , file = paste0(homewd,"/data/TableS7.csv"), row.names = F)





#and get a prediction for each province
pred.df <- cbind.data.frame(provname=rep(unique(dat$provname), each = length(unique(dat$year))) , year=rep(sort(unique(dat$year)), length(unique(dat$provname))))
#eliminate Tboung Khmum before data starts
pred.df$keep = 1
pred.df$keep[pred.df$provname=="Tboung Khmum" & pred.df$year< min(dat$year[dat$provname=="Tboung Khmum"])] <- 0

pred.df = subset(pred.df, keep==1)
pred.df$age_pred <- predict(m1, newdata = pred.df, type="response")
pred.df$age_pred_lci <- predict(m1, newdata = pred.df, type="response") - 1.96*(predict(m1, newdata = pred.df, type="response", se.fit = T)$se.fit)
pred.df$age_pred_uci <- predict(m1, newdata = pred.df, type="response") + 1.96*(predict(m1, newdata = pred.df, type="response", se.fit = T)$se.fit)




# add to plot

FigS15 <- ggplot(dat) + 
  geom_vline(xintercept = 2007, linetype=2, size=.3, color="gray65") + 
  geom_vline(xintercept = 2012, linetype=2, size=.3, color="gray65") + 
  geom_vline(xintercept = 2019, linetype=2, size=.3, color="gray65") +
  geom_jitter(aes(x=year, y=age), width=.1, size=.1, alpha=.3) + 
  facet_wrap(~provname) + theme_bw() + ylab("age of dengue infection") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), strip.background = element_rect(fill="white")) +
  geom_violin(aes(x=year,y=age, group=year, color = provname),  
              draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) +
  geom_ribbon(data=pred.df, aes(x=year, ymin=age_pred_lci, ymax=age_pred_uci), alpha=.3) +
  geom_line(data=pred.df, aes(x=year, y=age_pred))  




FigS15

ggsave(filename = paste0(homewd, "/final-figures/FigS15.png"),
       plot = FigS15,
       units="mm",  
       width=120, 
       height=85, 
       scale=3, 
       dpi=300)


#and a cumulative one for inclusion in the actual paper
#refit with no provname factored into slope

m2 <- gam(age ~ year + 
            s(provname, bs = "re"),
          data=dat)

summary(m2)
AIC(m1,m2)


pred.df.national <- cbind.data.frame(provname=rep("Battambang", length(unique(dat$year))), year=sort(unique(dat$year)))
pred.df.national$age_pred <- predict(m2, newdata = pred.df.national, type="response", exclude = "s(provname)")
pred.df.national$age_pred_lci <- predict(m2, newdata = pred.df.national, type="response", exclude = "s(provname)") - 1.96*(predict(m2, newdata = pred.df.national, exclude = "s(provname)", type="response", se.fit = T)$se.fit)
pred.df.national$age_pred_uci <- predict(m2, newdata = pred.df.national, type="response", exclude = "s(provname)") + 1.96*(predict(m2, newdata = pred.df.national, exclude = "s(provname)", type="response", se.fit = T)$se.fit)



Fig3B <- ggplot(data=dat) + 
  geom_vline(xintercept = 2007, linetype=2, size=.3, color="lightblue") + 
  geom_vline(xintercept = 2012, linetype=2, size=.3, color="lightblue") + 
  geom_vline(xintercept = 2019, linetype=2, size=.3, color="lightblue") +
  geom_jitter(aes(x=year, y=age), width=.09, size=.09, alpha=.2, show.legend = F) + 
  theme_bw() + ylab("age of dengue infection") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(),
        axis.text=element_text(size=13), axis.title.y = element_text(size=16),
        plot.margin = unit(c(.2,.2,.2,.2), "cm"))+
  geom_violin(aes(x=year,y=age, group=year),  color="gray55",
              draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) +
  geom_ribbon(data=pred.df.national, aes(x=year, ymin=age_pred_lci, ymax=age_pred_uci), fill="tomato", alpha=.3) +
  geom_line(data=pred.df.national, aes(x=year, y=age_pred), color="tomato", size=1) 


# and add in a panel with the mean age of infection
# or the slope of the increasing age from the GAM on a map

#load map

#and put these all together with map
cam = sf::st_read(paste0(homewd, "/data/province-shape/khm_admbnda_adm1_gov_20181004.shp"))
head(cam)
unique(cam$ADM1_EN)
unique(cam$ADM1_EN)
cam$ADM1_EN[cam$ADM1_EN=="Siemreap"] <- "Siem Reap"

# join with mean age and slope age
age.df <- cbind.data.frame(provname = unique(pred.df$provname), slope_age = summary(m1)$p.coeff[2:length(summary(m1)$p.coeff)])
mean.df <- ddply(dat, .(provname), summarise, mean_age=mean(age))
mean.year.df <- ddply(dat, .(provname,year), summarise, mean_age=mean(age))
mean.year.df = subset(mean.year.df, year==2020)
mean.year.df <- dplyr::select(mean.year.df, -(year))
names(mean.year.df)[names(mean.year.df)=="mean_age"] <- "mean_age_2020"
age.df <- merge(age.df, mean.df, by="provname")
names(age.df)[names(age.df)=="provname"] <- "ADM1_EN"
names(mean.year.df)[names(mean.year.df)=="provname"] <- "ADM1_EN"

cam.merge <- merge(cam, age.df, by="ADM1_EN")
cam.merge <- merge(cam.merge, mean.year.df, by="ADM1_EN")

Fig3A <- ggplot(cam.merge) +
  geom_sf(aes(fill=mean_age_2020, color=ADM1_EN), linewidth=.6) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm")) +
  scale_fill_gradient(low = 'lightgoldenrod1',  high = 'firebrick', 
                      na.value = "grey50",
                      trans="log10",
                      labels=scales::parse_format(),
                      name = "mean age of dengue\ninfection, 2020") +
  guides(color="none")


Fig3AB <- cowplot::plot_grid(Fig3A, Fig3B, nrow=1, ncol=2, labels = c("A", "B"), label_size = 22)


# Now lambda and sigma
fit.dat <- read.csv(file = paste0(homewd,"/data/prov-fits-FOI.csv"), header = T, stringsAsFactors = F)
age.fit <- read.csv(file = paste0(homewd,"/data/age-mult-profile.csv"), header = T, stringsAsFactors = F)
sigma.fit <- read.csv(file = paste0(homewd,"/data/sigma-fit.csv"), header = T, stringsAsFactors = F)

fit.dat$provname[fit.dat$provname=="Otdar Meanchey"] <- "Oddar Meanchey"
colz = scales::hue_pal()((length(unique(fit.dat$provname))-1))
names(colz) <- sort(unique(fit.dat$provname[fit.dat$provname!="National"]))
colz <- c(colz, "black")
names(colz)[length(colz)] <- "National"

#mult by 4?
# fit.dat$all_lambda <- fit.dat$N_sero*fit.dat$lambda
# fit.dat$all_lci <- fit.dat$N_sero*fit.dat$lci
# fit.dat$all_uci <- fit.dat$N_sero*fit.dat$uci

Fig3Cc <- ggplot(subset(fit.dat, provname!="Tboung Khmum")) + 
  geom_point(aes(x=year, y=lambda, color=provname), show.legend = F) +
  ylab(bquote("annual force of infection, per serotype, per capita,"~lambda)) +
  geom_ribbon(aes(x=year, ymin=lci, ymax=uci, fill=provname), alpha=.3, show.legend = F) +    
  geom_line(aes(x=year, y=lambda, color=provname), show.legend = F) +theme_bw() + 
  geom_point(data = subset(fit.dat, provname=="National"), aes(x=year, y=lambda), color="black") +
  geom_ribbon(data = subset(fit.dat, provname=="National"),aes(x=year, ymin=lci, ymax=uci), alpha=.3, fill="black") +    
  geom_line(data = subset(fit.dat, provname=="National"),aes(x=year, y=lambda), size =1, color="black") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=16),
        axis.text = element_text(size=13), plot.margin = unit(c(0,.2,.6,.4), "cm")) + #coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
  geom_vline(xintercept = 2007, linetype=2) + 
  geom_vline(xintercept = 2012, linetype=2) +
  geom_vline(xintercept = 2019, linetype=2) 

# Fig3Cc <- ggplot(subset(fit.dat, provname!="Tboung Khmum")) + 
#   geom_point(aes(x=year, y=lambda, color=provname), show.legend = F) +
#   ylab(bquote("annual force of infection, per capita,"~lambda)) +
#   geom_ribbon(aes(x=year, ymin=all_lci, ymax=all_uci, fill=provname), alpha=.3, show.legend = F) +    
#   geom_line(aes(x=year, y=all_lambda, color=provname), show.legend = F) +theme_bw() + 
#   geom_point(data = subset(fit.dat, provname=="National"), aes(x=year, y=all_lambda), color="black") +
#   geom_ribbon(data = subset(fit.dat, provname=="National"),aes(x=year, ymin=all_lci, ymax=all_uci), alpha=.3, fill="black") +    
#   geom_line(data = subset(fit.dat, provname=="National"),aes(x=year, y=all_lambda), size =1, color="black") +
#   theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=16),
#         axis.text = element_text(size=13), plot.margin = unit(c(0,.2,.6,.4), "cm")) + #coord_cartesian(ylim = c(0,1)) +
#   scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
#   geom_vline(xintercept = 2007, linetype=2) + 
#   geom_vline(xintercept = 2012, linetype=2) +
#   geom_vline(xintercept = 2019, linetype=2) 


# now add the birth rate panel to it at the top
pop.dat <- read.csv(file=paste0(homewd, "/data/world_bank_cambodia.csv"), header = T, stringsAsFactors = F)
head(pop.dat)

#get total births  - these are births per 1000 people
birth.vec <- pop.dat[1,5:ncol(pop.dat)]
death.vec <- pop.dat[3,5:ncol(pop.dat)]
names(birth.vec) <- seq(1960,2020,1) #assume this is census at the end of each year,
names(death.vec) <- seq(1960,2020,1) #assume this is census at the end of each year,
birth.vec <- birth.vec[which(names(birth.vec)=="1960"):which(names(birth.vec)=="2020")]
death.vec <- death.vec[which(names(death.vec)=="1960"):which(names(death.vec)=="2020")]
birth.vec['2020'] <- birth.vec['2019'] #assume this is the same as prior year
death.vec['2020'] <- death.vec['2019'] #assume this is the same as prior year

dat.births <- cbind.data.frame(year=names(birth.vec),births=c(unlist(birth.vec)))
head(dat.births)
dat.births$year <- as.numeric(dat.births$year)
dat.births = subset(dat.births, year >= min(fit.dat$year[fit.dat$provname!="Tboung Khmum"]))

dat.deaths <- cbind.data.frame(year=names(death.vec),deaths=c(unlist(death.vec)))
head(dat.deaths)
dat.deaths$year <- as.numeric(dat.deaths$year)
dat.deaths = subset(dat.deaths, year >= min(fit.dat$year[fit.dat$provname!="Tboung Khmum"]))


Fig3Ca <- ggplot(data=dat.births) + 
  geom_line(aes(x=year, y=births), show.legend = F, size=1) + 
  theme_bw() + ylab("births per\n1000 ppl") +
  theme(panel.grid = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(.2,.2,0,.31), "cm"),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=16)) +
  geom_vline(xintercept = 2007, linetype=2) + 
  geom_vline(xintercept = 2012, linetype=2) +
  geom_vline(xintercept = 2019, linetype=2) 


Fig3Cb <- ggplot(data=dat.deaths) + 
  geom_line(aes(x=year, y=deaths), show.legend = F, size=1) + 
  theme_bw() + ylab("deaths per\n1000 ppl") +
  theme(panel.grid = element_blank(), 
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(.2,.2,0,.31), "cm"),
        axis.text.y = element_text(size=13),
        axis.title.y = element_text(size=16)) +
  geom_vline(xintercept = 2007, linetype=2) + 
  geom_vline(xintercept = 2012, linetype=2) +
  geom_vline(xintercept = 2019, linetype=2) 


Fig3C <- cowplot::plot_grid(Fig3Ca, Fig3Cb, Fig3Cc, nrow=3, ncol = 1, rel_heights = c(.24,.24,1), labels = c("C", "", ""), label_size = 22, label_y = 1.1)

#store the supp figure 


FigS16 <- ggplot(subset(fit.dat, provname!="Tboung Khmum")) + facet_wrap(~provname) +
  geom_point(aes(x=year, y=lambda, color=provname), show.legend = F) +
  ylab(bquote("annual force of infection, per capita,"~lambda)) +
  geom_ribbon(aes(x=year, ymin=lci, ymax=uci, fill=provname), alpha=.3, show.legend = F) +    
  geom_line(aes(x=year, y=lambda, color=provname), show.legend = F) +theme_bw() + 
  geom_point(data = subset(fit.dat, provname=="National"), aes(x=year, y=lambda), color="black") +
  geom_ribbon(data = subset(fit.dat, provname=="National"),aes(x=year, ymin=lci, ymax=uci), alpha=.3, fill="black") +    
  geom_line(data = subset(fit.dat, provname=="National"),aes(x=year, y=lambda), size =1, color="black") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size=14), strip.text = element_text(size=14), 
        axis.text = element_text(size=12), plot.margin = unit(c(.2,.2,.2,.2), "cm"), strip.background = element_rect(fill="white")) + #coord_cartesian(ylim = c(0,1)) +
  geom_vline(xintercept = 2007, linetype=2) + scale_color_manual(values=colz) + scale_fill_manual(values=colz) +
  geom_vline(xintercept = 2012, linetype=2) + scale_x_continuous(breaks = c(1990, 2000, 2010, 2020)) +
  geom_vline(xintercept = 2019, linetype=2) 

print(FigS16)


ggsave(filename = paste0(homewd, "/final-figures/FigS16.png"),
       plot = FigS16,
       units="mm",  
       width=120, 
       height=85, 
       scale=3, 
       dpi=300)


# #now, the age multiplier becomes panel D
# #age.mult.profile <- read.csv(file = paste0(homewd, "/data/age-mult-profile.csv"), header=T,stringsAsFactors = F)
# 
# 
# #now, the age multiplier becomes panel D
# load(paste0(homewd, "/figure-development/Fig3/age-fit-alt-two-short-trim/refit-many-age-mult.Rdata"))

# embed the age multiplier
#age.long <- melt(age.mult.profile, id.vars = c("age_mult", "lci_mult", "uci_mult", "year_range"), measure.vars = c("age_min", "age_max"))
age.long <- melt(age.fit, id.vars = c("age_mult", "lci_mult", "uci_mult", "year_range"), measure.vars = c("age_min", "age_max"))
head(age.long)
names(age.long)[names(age.long)=="value"] <- "age"

Fig3D <- ggplot(age.long) + theme_bw() + facet_grid(~year_range) + theme_bw()+
  geom_hline(aes(yintercept=1), linetype=2, color="gray85") + coord_cartesian(ylim = c(0,3.5))+
  geom_line(aes(x=age, y=age_mult), size=1)  + ylab(bquote("age multiplier for"~lambda))+ 
  geom_ribbon(aes(x=age, ymin=lci_mult, ymax=uci_mult), alpha=.2) +
  theme(panel.grid = element_blank(), axis.title = element_text(size=16), 
        plot.margin = unit(c(.2,.2,.2,.7), "cm"), 
        strip.background = element_rect(fill="white"), strip.text = element_text(size=16),
        axis.text = element_text(size=12)) 
print(Fig3D)


# then, do the cumulative age vs. model at the national level


# first, load the functions
model.age.incidence.series.age.mult.wane <- function(par.dat, age_vect, year.start, age_mult){
  
  
  # you will have one lambda for each year in the time series
  # and one N-sero for each year as well
  # and X ages with varying probabilities of infection within those years
  
  lambda.start = par.dat$lambda[1]
  
  
  
  lts = nrow(par.dat)
  
  #proportion exposed in any way - needs to long enough for each age within each year
  pexposed= rep(NA, (length(age_vect)-1))
  pexposed[1] <- 0
  pexposed = rep(list(pexposed),lts) 
  
  #naive
  pnaive= rep(NA, (length(age_vect)-1))
  pnaive[1] <- 1
  pnaive = rep(list(pnaive),lts) 
  
  
  #exposed to a single serotype only (this is one serotype's worth of primary infections)
  pprim= rep(NA, (length(age_vect)-1))
  pprim[1] <- 0
  pprim = rep(list(pprim),lts)
  
  #multitypic infections
  pmulti= rep(NA, (length(age_vect)-1))
  pmulti[1] <- 0
  pmulti = rep(list(pmulti),lts)
  
  age_tracker= rep(NA, (length(age_vect)-1))
  age_tracker[1] <- 0
  age_tracker = rep(list(age_tracker),lts)
  
  year_tracker= rep(NA, (length(age_vect)-1))
  year_tracker = rep(list(year_tracker),lts)
  year_tracker[[1]][1] <- year.start
  
  year.start= min(par.dat$year)
  
  
  for (i in 1:lts){
    
    if(i < max(age_vect)){
      age_vect_tmp = age_vect[1:which(age_vect==i)]
    }else{
      age_vect_tmp = age_vect
    }
    
    for(a in 2:(length(age_vect_tmp))){ 
      # for-loops over all possible ages in our data range across all the years. 
      # all ages get rounded up, so babies <1 year are counted as 1 year old,
      # and experienced a hazard of infection in the year prior
      
      
      
      
      age = age_vect_tmp[a]
      
      year.now = par.dat$year[i]
      year.par = subset(par.dat, year<=year.now)
      N_sero = year.par$N_sero
      lambda =  year.par$lambda #one per year
      sigma = year.par$sigma
      # make a vector of durations for the lambdas across the time series
      # in this model, it is easy because the duration is only 1 year...
      # you may eventually also want partial year durations in a different model
      #if (age>1){
      dur = rep(1, length(lambda))  
      #}else if(age<=1){
      # dur = rep(0, length(lambda))  
      #}
      
      #if you were born after the time series began
      #some of those ones may need to be replaced with 0s
      if(i>1){
        #diff.count <- i-age
        tot.class = ceiling(age)
        #diff.count <- sum(dur)-age_trunc
        
        if(tot.class<length(dur)){
          dur[1:(length(dur)-tot.class)]<- 0
        }
      }
      
      # #if 0, this kid spent 1 year in each age class and the dur vector does not need to be changed
      if(sum(dur) != age){
        warning(paste0("mismatch in age and integration at age=", age, " and i=", i))
      }
      
      
      #1-2,3-4,5-6,7-9,10-12,13-15,16-19,20-29,30-39,40+
      #1-2,3-4,5-6,7-9,10-12,13-15,16-19,20-29,30-39,40+
      # if(year.now <=2010){
      #   age.mult.df <- cbind.data.frame(mult=age_mult[1:13], 
      #                                   age_min=c(1,3,5,7,10,13,16,20,30,40,50,60,70),
      #                                   age_max=c(2,4,6,9,12,15,19,29,39,49,59,69,90))  
      # }else if (year.now >2010){
      #   age.mult.df <- cbind.data.frame(mult=age_mult[14:26], 
      #                                   age_min=c(1,3,5,7,10,13,16,20,30,40,50,60,70),
      #                                   age_max=c(2,4,6,9,12,15,19,29,39,49,59,69,90))  
      # }
      # 
      
      
      if(year.now <=2010){
        age.mult.df <- cbind.data.frame(mult=age_mult[1:8], 
                                        age_min=c(1,3,5,7,10,13,16,20),
                                        age_max=c(2,4,6,9,12,15,19,90))  
      }else if (year.now >2010){
        age.mult.df <- cbind.data.frame(mult=age_mult[9:19], 
                                        age_min=c(1,3,5,7,10,13,16,20,30,40,50),
                                        age_max=c(2,4,6,9,12,15,19,29,39,49,90))  
      }
      
      age.mult.df$dur = (age.mult.df$age_max-age.mult.df$age_min)+1
      
      #rep it as a multiplier by year
      age.mult.vect = c(unlist(mapply(rep, x=as.list(age.mult.df$mult), each=as.list(age.mult.df$dur), SIMPLIFY = F)))
      
      #now, make it the length of dur for which there are not 0s
      
      age.mult.vect = age.mult.vect[1:length(dur[dur>0])]
      
      age.mult.vect = c(dur[dur==0], age.mult.vect)
      
      if(length( age.mult.vect) != length(dur)){
        warning(paste0("age multiplier incorrect length at age=", age, " and i=", i))
      }
      
      
      #just to be safe, convert the 0s to 1s so it does not modulate anything
      age.mult.vect[age.mult.vect==0] <- 1
      
      # now, modulate lambda by multiplying with age.mult.vect
      lambda.age = lambda*age.mult.vect
      
      
      # now, we integrate the hazard of becoming infected over all the years in each
      
      
      #First, integrand A (hazard of exposure to all strains )
      inte_A = (sum(dur*lambda*(N_sero)))
      
      #Then, hazard of exposure to strain i only
      inte_B = (sum(dur*lambda))
      
      #Then, probability of exposure to strain i and then strain k (where k does not equal i)
      inte_C = (sum(dur*lambda*(2)))
      
      
      #Then, cumulative waning immunity hazard for one strain
      inte_E = sum(sigma*dur) 
      
      #And cumulative waning hazard for all the other strains
      inte_D = (sum((N_sero-1)*sigma*dur))
      
      
      
      #could also write as sigma[1]*(a-1) (since a starts at class 2 we need to subtract 1 from a)
      
      
      # now fill in the probabilities based on these : 
      #this is the probability of being naive to all serotypes at this point in time
      pnaive[[i]][[a]] <- exp(-inte_A)
      
      #this is the probability of being exposed to any serotype at the current time
      pexposed[[i]][[a]] <- (1-exp(-inte_A))
      
      # this is the probability of being exposed to just one serotype (a primary infection)
      # here, this is expressed for just a single target strain 
      # (would need to multiply by the number of circulating strains if you wanted 
      # the probability of a primary infection with ANY strain)
      pprim[[i]][[a]] <- exp(-inte_A)*(exp(inte_B)-1) + (1-exp(-inte_C))*(exp(-inte_D))*(1-exp(-inte_E))
      
      
      # if not primarily infected or naive, this should be a multitypic infection
      # but remember it could be a primary infection with any of the four serotypes 
      # so you need to account for that here by again multiplying by the 
      # number of circulating serotypes internally
      pmulti[[i]][[a]] = 1 - pnaive[[i]][[a]] - (N_sero[i]*(pprim[[i]][[a]]))
      age_tracker[[i]][[a]] = age
      year_tracker[[i]][[a]] =year.now
      
    }
  }
  
  #some of the age distributions will be shorter
  
  #and get the estimates of each
  #p.out = cbind.data.frame(year = rep(seq((year.start), (year.start+lts-1), 1), each=length(age_vect)),
  p.out <-  cbind.data.frame(year=c(unlist(year_tracker)),
                             age=c(c(unlist(age_tracker))), 
                             exposed=c(c(unlist(pexposed))),
                             naive=c(c(unlist(pnaive))),
                             one_prim=c(unlist(pprim)), 
                             multi = c(unlist(pmulti)))
  
  
  # merge to add N_sero
  add.par <- dplyr::select(par.dat, year, N_sero)
  
  p.out <- merge(p.out, add.par, by="year", all.x = T, sort = F)
  
  
  p.out$all_prim <- p.out$one_prim*p.out$N_sero
  
  
  
  p.out[p.out<0] <- 0
  
  #check the proportions
  p.out$sum_exp_naive <- rowSums(cbind(p.out["exposed"], p.out["naive"])) #should be 1
  p.out$sum_naive_prim_multi <- rowSums(cbind(p.out["naive"], p.out["all_prim"], p.out["multi"])) #should be 1
  
  
  p.out <- arrange(p.out, year, age)
  
  # ggplot(data=p.out) + geom_point(aes(x=age, y=multi))  +  geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  #and the multitypic distribution corresponds to the cumulative proportion of cases
  p.out$cum_prop_cases <- p.out$multi
  p.out <- p.out[complete.cases(p.out),]
  
  
  p.add <- cbind.data.frame(year = unique(p.out$year), age = rep(0, length(unique(p.out$year))), exposed=0, naive=1, one_prim=0, multi=0, N_sero = N_sero[1], all_prim=0, sum_exp_naive=1, sum_naive_prim_multi=1, cum_prop_cases =0)
  
  
  p.sum <- rbind(p.add, p.out)
  p.sum <- arrange(p.sum, year, age)
  p.sum$provname = unique(par.dat$provname)
  
  
  p.sum <- p.sum[!duplicated(p.sum),]
  
  p1 <- ggplot(data=p.sum) + geom_point(aes(x=age, y=multi))  + 
    geom_line(aes(x=age, y=multi)) +facet_wrap(~year)
  
  
  # print(p1)
  # 
  #ggplot(data=subset(p.sum, year==1985)) + geom_point(aes(x=age, y=multi))  + 
  #geom_line(aes(x=age, y=multi)) 
  
  return(p.sum) #returns prevalence by age for each year for fitting to the years for which we have data
}
sum.yr.all.prov <- function(df){
  
  df.sum <- ddply(df, .(year, age), summarise, Nage = length(age))
  
  
  df.out = cbind.data.frame(age=1:max(df$age))
  
  df.out <- merge(df.out, df.sum, by="age", all.x = T, sort = F)
  df.out$Nage[is.na(df.out$Nage)] <- 0
  df.out$year[is.na(df.out$year)] <- unique(df$year)
  #df.out <- rbind(c(0,0), df.out)
  #bind
  df.add <- cbind.data.frame(age=0, year = unique(df$year), Nage=0)
  df.out <- rbind(df.add, df.out)
  df.out <- arrange(df.out, age)
  df.out$cum_cases <- cumsum(df.out$Nage)
  df.out$n <- sum(df.sum$Nage)
  df.out$cum_prop_cases <- df.out$cum_cases/df.out$n
  df.out$provname = unique(df$provname)
  return(df.out)
  
}
get.min.year <- function(df){
  out = min(df$year)
  return(out)
}
run.model.data.all <- function(dat,par.dat, sigma.fit, age.mult.df){
  
  
  
  
  #for the first year in the dataset, 
  #estimate foi just from the cross-sectional data
  year.start = min(par.dat$year)
  
  nat.dat <- dat
  nat.dat$provname  = "National"
  dat <- rbind(dat, nat.dat)
  
  #first, prep the data
  year.dat <- dlply(dat, .(provname, year))
  
  
  year.dat.sum <- lapply(year.dat, sum.yr.all.prov)
  
  df.out <- data.table::rbindlist( year.dat.sum)
  
  
  #and run the models across all the provinces
  
  age_vect=seq(0, max(dat$age), by=1)
  
  #merge sigma with par.dat 
  #and sigma fit only
  sigma.merge <- dplyr::select(sigma.fit, year, sigma)
  par.dat <- merge(par.dat, sigma.merge, by="year", all.x = T)
  par.dat$sigma[is.na(par.dat$sigma)] <- 0
  
  
  #now, split the par.dat by provinces and apply simultaneously
  par.dat.split <- dlply(par.dat, .(provname))
  
  #make year start split
  year.start.split <- lapply(par.dat.split, get.min.year)
  
  #apply function across all par and year start subsets
  # age mod and age vect are constant
  
  out.mod.list <- mapply(model.age.incidence.series.age.mult.wane, par.dat= par.dat.split, year.start=year.start.split, MoreArgs = list(age_mult = age.mult.df$age_mult, age_vect=age_vect), SIMPLIFY = FALSE)
  
  # and bind
  
  out.mod.df <- data.table::rbindlist(out.mod.list)
  #ggplot(out.mod.df) + geom_line(aes(x=age, y=cum_prop_cases, color=year)) + facet_wrap(~provname)
  #head(out.mod.df)
  out.mod.df <- arrange(out.mod.df, provname, year, age)
  
  df.out <- arrange(df.out, provname, year, age)
  
  min.df.prov <- ddply(df.out, .(provname), summarise, min_year = min(year))
  out.mod.df <- merge(out.mod.df, min.df.prov, by="provname", all.x = T, sort = F)
  
  # and remove modeling projections from years with no data
  out.mod.df$year_early <- 1
  out.mod.df$year_early[out.mod.df$year>=out.mod.df$min_year] <- 0
  out.mod.df <- subset(out.mod.df, year_early==0)
  mod.df <- dplyr::select(out.mod.df, year, age, cum_prop_cases, provname)
  df.merge <- dplyr::select(df.out, year, age, cum_prop_cases, provname)
  names(mod.df) <- c("year", "age", "cum_prop_cases", "provname")
  mod.df$type <- "model"
  names(df.merge) <- c("year", "age", "cum_prop_cases", "provname")
  df.merge$type <- "data"
  
  #and combine
  all.mod.dat <- rbind(df.merge, mod.df)
  return(all.mod.dat)
  
}

#run the model and arrange the data
#out.plot <- run.model.data.all(par.dat = fit.dat, dat = dat, age.mult.df = age.mult.profile, sigma.fit = sigma.fit)
out.plot <- run.model.data.all(par.dat = fit.dat, dat = dat, age.mult.df = age.fit, sigma.fit = sigma.fit)
head(out.plot)
out.plot <- arrange(out.plot, provname, year)
out.plot$provname <- factor(out.plot$provname, levels = unique(out.plot$provname))
out.plot$year <- as.factor(out.plot$year)
linez = c('data'=3 , 'model'=1)

#first the supplement with all the provinces

FigS17 <- ggplot(subset(out.plot, provname!="Tboung Khmum")) + facet_wrap(~provname, ncol=5) +
  geom_line(aes(x=age, y=cum_prop_cases, linetype=type, color=year), alpha=.8, size=.8) +
  scale_linetype_manual(values=linez) +
  theme_bw() + ylab("cumulative proportion cases") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(size=14), strip.text = element_text(size=14), 
        axis.text = element_text(size=12), plot.margin = unit(c(.2,.2,.2,.2), "cm"), 
        legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_blank(), strip.background = element_rect(fill="white")) +
  guides(linetype="none", color=guide_legend(ncol=8)) + 
  scale_color_viridis_d(direction = -1)

ggsave(filename = paste0(homewd, "/final-figures/FigS17.png"),
       plot = FigS17,
       units="mm",  
       width=120, 
       height=85, 
       scale=3, 
       dpi=300)

# now the one with the national scale
Fig3Ea <- ggplot(subset(out.plot, provname=="National")) +
  geom_line(aes(x=age, y=cum_prop_cases, linetype=type, color=year), alpha=.8, size=.8) +
  scale_linetype_manual(values=linez) +
  theme_bw() + ylab("cumulative proportion cases (national)") +
  theme(panel.grid = element_blank(),
        legend.text = element_text(size=12),
        axis.title = element_text(size=16), strip.text = element_text(size=14), 
        axis.text = element_text(size=13), plot.margin = unit(c(.2,.2,.05,.2), "cm"), 
        legend.position = c(.85,.46), legend.box = "horizontal",legend.key.width = unit(c(.4), "cm"), 
        legend.title = element_blank(), strip.background = element_rect(fill="white")) +
  guides( color=guide_legend(ncol=1)) + 
  scale_color_viridis_d(direction = -1, option = "turbo")


#embed sigma
head(sigma.fit)

Fig3Eb <- ggplot(sigma.fit) + theme_bw() + 
  geom_line(aes(x=year, y=sigma), size=1)  + ylab(bquote(atop("annual rate of waning","multitypic immunity,"~sigma)))+ 
  geom_ribbon(aes(x=year, ymin=lci_sigma, ymax=uci_sigma), alpha=.3) +
  theme(panel.grid = element_blank(), axis.title = element_text(size=16), axis.title.x = element_blank(),
        axis.text = element_text(size=11))
print(Fig3Eb)


#and add it in as an inset
Fig3E <- Fig3Ea + annotation_custom(ggplotGrob(Fig3Eb), 
                                      xmin = 19, xmax = 78, ymin = .01, ymax = .54)


print(Fig3E)


# and pair with Fig3D
Fig3DE <- cowplot::plot_grid(Fig3D, Fig3E, ncol=1, nrow=2, labels = c("D", "E"), label_size = 22, rel_heights = c(.45,1), label_y = c(1.05))


#and bottom 2/3 of plot
Fig3low <- cowplot::plot_grid(Fig3C, Fig3DE, ncol=2, nrow=1)

# and all together
Fig3 <- cowplot::plot_grid(Fig3AB, Fig3low, ncol=1, nrow = 2, rel_heights = c(.5,1)) + theme(plot.background = element_rect(fill="white"))


ggsave(filename = paste0(homewd, "/final-figures/Fig3.png"),
       plot = Fig3,
       units="mm",  
       width=120, 
       height=110, 
       scale=3, 
       dpi=300)

# 
# 
 ggsave(filename = paste0(homewd, "/final-figures/Fig3.pdf"),
        plot = Fig3,
        units="mm",  
        width=120, 
        height=110, 
       scale=3, 
        dpi=300)
# 






