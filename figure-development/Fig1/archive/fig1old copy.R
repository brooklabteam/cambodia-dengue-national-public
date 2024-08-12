

library(plyr)
library(dplyr)
library(ggplot2)

rm(list=ls())
# This file is for figure 1. 

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"




#see the script "generate-combined-tsir-increased-S-increased-beta-dataset.R" 
# to build the datasets that are used here


# set up size and color
xy_size<-13
my_red="tomato"
my_blue<-"#2980B9"
my_orange<-"#D35400"
my_green<-"seagreen"
my_yellow<-"#F1C40F"
color_2007<-"purple"
color_2012<-"turquoise3"
color_2019<-"tomato"


#epidemic years 2007, 2012, 2019

my_blue<-"#2980B9"
my_orange<-"#D35400"
my_red<-"tomato"
my_green<-"seagreen"
my_yellow<-"#F1C40F"
my_purple <- "purple"


#load data 
dat <- read.csv(file =paste0(homewd, "/data/TSIR_fitted_timeseries_estimates.csv"),header =T, stringsAsFactors = F)

#and plot all
#first, convert to something that can be parsed
dat$variable[dat$variable=="TSIR fit"] <- "TSIR~fit"
dat$variable[dat$variable=="prediction TSIR"] <- "prediction~TSIR"
dat$variable[dat$variable=="prediction TSIR-increased S"] <- "prediction~TSIR-increased~S"
dat$variable[dat$variable=="prediction TSIR-increased Beta"] <- "prediction~TSIR-increased~beta"
# 
# dat$variable = factor(dat$variable, levels=c("data",
#                                              "TSIR fit",
#                                              "prediction TSIR",
#                                              "prediction TSIR-increased S",
#                                              "prediction TSIR-increased Beta"))

dat$variable = factor(dat$variable, levels=c("data",
                                             "TSIR~fit",
                                             "prediction~TSIR",
                                             "prediction~TSIR-increased~S",
                                             "prediction~TSIR-increased~beta"))

#colz = c('data' = "black", 'TSIR fit' = my_blue, 'prediction TSIR' = my_green, 'prediction TSIR-increased S' =my_purple, 'prediction TSIR-increased Beta' = my_red)
#typez = c('data' = 2, 'TSIR fit' = 1, 'prediction TSIR' = 1, 'prediction TSIR-increased S' = 1, 'prediction TSIR-increased Beta' = 1)


colz = c('data' = "black", 'TSIR~fit' = my_blue, 'prediction~TSIR' = my_green, 'prediction~TSIR-increased~S' =my_purple, 'prediction~TSIR-increased~beta' = my_red)
typez = c('data' = 2, 'TSIR~fit' = 1, 'prediction~TSIR' = 1, 'prediction~TSIR-increased~S' = 1, 'prediction~TSIR-increased~beta' = 1)


library(scales)
pl_a <- ggplot(data=dat, aes(x=time, y=reported_cases, color=variable, linetype=variable)) + 
  geom_line(data=dplyr::filter(dat, variable=="data"), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR~fit" & time<2007), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR~fit" & time>=2008 & time<2012), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="TSIR~fit" & time>=2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction~TSIR" & time>=2007 & time<2008), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction~TSIR" & time>=2012 & time<2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction~TSIR" & time>=2019 & time<2020), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction~TSIR-increased~S" & time>=2007 & time<2008), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction~TSIR-increased~S" & time>=2012 & time<2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction~TSIR-increased~S" & time>=2019 & time<2020), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction~TSIR-increased~beta" & time>=2007 & time<2008), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction~TSIR-increased~beta" & time>=2012 & time<2013), size=1) +
  geom_line(data=dplyr::filter(dat, variable=="prediction~TSIR-increased~beta" & time>=2019 & time<2020), size=1) +
  scale_color_manual(values=colz, labels = label_parse()) + # coord_cartesian(ylim = c(0,1500))+
  scale_linetype_manual(values=typez, labels = label_parse()) +theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=18),
        legend.text = element_text(size=12),
        legend.position = c(.25,.82),
        axis.text = element_text(size=14), legend.title = element_blank()) +
  ylab(paste0("reported cases"))


plot(pl_a)








######################################################################################
#### Fig B | log ######################################
######################################################################################

library(tsiR)

sim.fit.return.beta <- function(time.start, dat, epiyr, family){
  
  
  if(family=="poisson"){
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                          IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                          regtype='lm',family='poisson',link='log')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                              IP = 2,
                              parms=fittedpars,
                              #epidemics='break', threshold=3,
                              method='pois', nsim=100)
    
  }else if(family=="gaussian"){
    
    fittedpars <- estpars(data=subset(dat, time >= time.start & time<epiyr),
                                                 IP=2, alpha=NULL, sbar=NULL, xreg = "cumcases",
                                                 regtype='lm',family='gaussian')
    
    simfitted <- simulatetsir(data=subset(dat, time >= time.start & time<epiyr),
                                                         IP = 2,
                                                         parms=fittedpars,
                                                         #epidemics='break', threshold=3,
                                                         #method='pois',
                                                         nsim=100)
                           
    
  }
  
   
   
   df.out = fittedpars$contact
   df.out$epiyr <- epiyr
   
  
  
  return(df.out)
}

#load data into tsir form
# tsir_data.csv is generated by make-tsir-data.R
tsir.dat <- read.csv(paste0(homewd, "/data/tsir_data.csv"), header = T, stringsAsFactors = F)

head(tsir.dat)

#first try to predict 2007-epidemic
#here we fit to 2002:2007

## @Yimei, I know you made some changes to the fitting script for TSIR
## Below assumes regtype = "lm" and family = "gaussian".
## Please make sure it matches

sim.2007 <- sim.fit.return.beta(dat = tsir.dat,
                        time.start =  min(tsir.dat$time),
                        family="gaussian",
                        epiyr = 2007)


#now for 2012
sim.2012 <- sim.fit.return.beta(dat = tsir.dat,
                        time.start =  min(tsir.dat$time[tsir.dat$time>=2008]),
                        family="gaussian",
                        epiyr = 2012)


#and the results for 2019
sim.2019 <- sim.fit.return.beta(dat = tsir.dat,
                        time.start =  min(tsir.dat$time[tsir.dat$time>=2013]),
                        family="gaussian",
                        epiyr = 2019)

#and combine
fit.all <- rbind(sim.2007, sim.2012, sim.2019)
fit.all$epiyr[fit.all$epiyr==2007] <- "2002-2006"
fit.all$epiyr[fit.all$epiyr==2012] <- "2008-2011"
fit.all$epiyr[fit.all$epiyr==2019] <- "2013-2018"
fit.all$epiyr <- as.factor(fit.all$epiyr)

fit.all$beta_per_1000 <- fit.all$beta*1000
fit.all$beta_low_per_1000 <- fit.all$betalow*1000
fit.all$beta_hi_per_1000 <- fit.all$betahigh*1000

p1 <- ggplot(data=fit.all) +#facet_grid(epiyr~., scales = "free_y") +
      geom_line(aes(x=time, y=beta, color=epiyr)) +
      geom_ribbon(aes(x=time, ymin=betalow, ymax=betahigh, fill=epiyr), alpha=.3)




#or as lines and points
pl_b <- ggplot(data=fit.all) + #facet_grid(epiyr~., scales = "free_y") +
  geom_point(aes(x=time, y=beta, color=epiyr)) +
  geom_linerange(aes(x=time, ymin=betalow, ymax=betahigh, color=epiyr), alpha=.3)+
  geom_line(aes(x=time, y=beta, color=epiyr)) +
  #geom_ribbon(aes(x=time, ymin=betalow, ymax=betahigh, fill=epiyr), alpha=.3)+
  theme_bw() + scale_y_log10()+
  scale_color_manual(values=c( "darkcyan", "darkorchid3", "navy"), 
                     name="fitting period") +
  scale_fill_manual(values=c( "darkcyan", "darkorchid3", "navy"), 
                       name="fitting period") + 
  ylab(bquote(log[10](beta)~',transmission rate')) +
  #xlab("week of year")+
  scale_x_continuous(breaks=c(1*2,3*2,5*2,7*2,9*2,11*2), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
  theme(panel.grid = element_blank(), 
        legend.position = c(.128,.885),
        legend.text = element_text(size=14),
        strip.background = element_rect(fill="white"),
        axis.title.y = element_text(size=16), 
        axis.title.x = element_blank(),
        axis.text = element_text(size=14), legend.title = element_blank(),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))

#with per 1000
#or as lines and points
pl_b <- ggplot(data=fit.all) + #facet_grid(epiyr~., scales = "free_y") +
  geom_point(aes(x=time, y=beta_per_1000, color=epiyr)) +
  geom_linerange(aes(x=time, ymin=beta_low_per_1000, ymax=beta_hi_per_1000, color=epiyr), alpha=.3)+
  geom_line(aes(x=time, y=beta_per_1000, color=epiyr)) +
  #geom_ribbon(aes(x=time, ymin=betalow, ymax=betahigh, fill=epiyr), alpha=.3)+
  theme_bw() + scale_y_log10()+
  scale_color_manual(values=c( "darkcyan", "darkorchid3", "navy"), 
                     name="fitting period") +
  scale_fill_manual(values=c( "darkcyan", "darkorchid3", "navy"), 
                    name="fitting period") + 
  #ylab(bquote(log[10](beta)~',transmission rate')) +
  ylab(bquote(atop(beta~', biweekly transmission rate', '(/1000 ppl)'))) +
  #ylab(bquote(beta~', biweekly transmission rate (/1000 ppl)')) +
  #ylab("transmission rate, per 1000 ppl") +
  #xlab("week of year")+
  scale_x_continuous(breaks=c(1*2,3*2,5*2,7*2,9*2,11*2), labels = c("Jan", "Mar", "May", "Jul", "Sep", "Nov"))+
  theme(panel.grid = element_blank(), 
        legend.position = c(.128,.885),
        legend.text = element_text(size=14),
        strip.background = element_rect(fill="white"),
        axis.title.y = element_text(size=16), 
        axis.title.x = element_blank(),
        axis.text = element_text(size=14), legend.title = element_blank(),legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))








######################################################################################
#### Fig C | time & sus & color birth rate ######################################
######################################################################################



# Sdatcombined is generated by "generate-combined-tsir-increased-S-increased-beta-dataset.R"
Sdatcombined <- readRDS(paste0(homewd, "/data/Sdatcombined.rds"))
colz = c("other years"="black", "epidemic years"="red")
Sdatcombined



Sdatcombined$epiyr <- "other years"
for (eachrow in seq(1,dim(Sdatcombined)[1])){
  arow <- Sdatcombined[eachrow,]
  if (arow$time == 2007 | arow$time == 2012 | arow$time == 2019 ){
    Sdatcombined[eachrow,]$epiyr <- "epidemic years"
  }
    
}
# Sdatcombined




################################################################
# Fig 3C, time vs sus, birth change to birth/1000 ppl
################################################################


# this is from the original birth data with added  births_per_1000 column for each year. (The birth rate is actually declining while the population size is increasing.)
find_births<-read.csv(paste0(homewd,"/data/tsir_data_birth_updated.csv"))


# choose a number every 26 rows
births_update<- cbind.data.frame(seq(2002,2020,1), sapply(split(find_births$births_per_1000, rep(1:(nrow(find_births)/26), each=26)), mean))
colnames(births_update)<-c("time","births_per1000")
births_update


Sdatcombined_update<-merge(x = Sdatcombined, y =births_update, by = "time", all.x = TRUE)

Sdatcombined = subset(Sdatcombined_update, select = -c(births) )
names(Sdatcombined)[names(Sdatcombined) == 'births_per1000'] <- 'births per 1000 ppl'
Sdatcombined



colz = c("other years"="black", "epidemic years"="red")

Sdatcombined$epiyr <- "other years"
for (eachrow in seq(1,dim(Sdatcombined)[1])){
  arow <- Sdatcombined[eachrow,]
  if (arow$time == 2007 | arow$time == 2012 | arow$time == 2019 ){
    Sdatcombined[eachrow,]$epiyr <- "epidemic years"
  }
    
}

# remove TSIR and only show TSIR-increased S

Sdatcombined<-Sdatcombined[!(Sdatcombined$sim %in% "TSIR"),]



Sdatcombined$variable<-sub("-", "- \n ",Sdatcombined$variable)  


a<-c("black","red","red","black","red")

my_blue<-"#2980B9"
my_green<-"green3"
colz = c("prediction TSIR"=my_green, "prediction TSIR- \n increased S"="tomato","TSIR fit"=my_blue)
options(scipen=10000) # we don't want to plot with 1e+00


pl_c<-ggplot(data=Sdatcombined) + 
  geom_point(aes(x=year, y=sus_mean, 
                 color=variable, shape=variable), size=5, stroke = 1) + labs(shape=NULL, colour=NULL) +# only remove two legend
  scale_color_manual(values=colz) +

  ggnewscale::new_scale_color() +
  geom_point(aes(x=year, y=sus_mean, 
                 color=`births per 1000 ppl`, shape=variable), size=3) + 
  
  theme_bw() +
  theme(panel.grid = element_blank(),strip.background  =element_rect(fill="white")) +
  labs(x ="Year", y = "Susceptible Count")+
  theme(
  axis.title.x = element_text(color="black", size=14),
  axis.title.y = element_text(color="black", size=14))+
  scale_y_log10() + scale_color_viridis_c()+ geom_vline(xintercept = c(2007,2012,2019) , linetype=2, 
                                                                                        color = "red", size=0.5)+
  scale_x_continuous(limits=c(2002, 2020), breaks=c(2002,2007,2012,2015,2019))+
  theme(axis.text.x = element_text( hjust = 0.5, colour = a),strip.text = element_text(colour = 'black',size = 14))  +
  
  theme(
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size=16),
        legend.text = element_text(size=12),legend.title=element_text(size=14),
        legend.position = c(.847,.723),legend.box.background = element_rect(colour = "black"),legend.background = element_blank(),
        axis.text = element_text(size=14))+# rremove("xy.title")+# theme(legend.title = element_blank()) + 
  theme(legend.spacing.y = unit(0, "cm"))+ # space between legend rows
  theme(legend.title.align = 1.08)




######################################################################################
#### Fig D | reported case & predicted case TSIR & increasedS ########################
######################################################################################

dat_true<-dat[dat$variable=="data",]

dat_pred_update<-merge(x = dat, y =dat_true, by = "time", all.x = TRUE)

dat_df<-cbind.data.frame(dat_pred_update$time, dat_pred_update$reported_cases.x, dat_pred_update$reported_cases.y, dat_pred_update$variable.x,dat_pred_update$variable.y)

colnames(dat_df)<-c("time","reported_cases","true_cases","pred_type","data")




# make every year different shape


dat_pred_3<-dat_df
dat_pred_3$year<-floor(dat_pred_3$time)

non_pei_years<- as.data.frame(dplyr::filter(dat_pred_3,  (year>=2002 & year <2007) | (year>=2008 & year <2012) | (year>=2013 & year <2019) | (year>2019) ) )


dat_pred_3[dat_pred_3$time %in% non_pei_years$time,]$year<-"other_year"




dat <- dat_pred_3

colz = c(
  "prediction TSIR" = "seagreen",
  "prediction TSIR-increased S" = "tomato",
  "TSIR fit" = my_blue,
  "data" = "black"
)

colz = c(
  "prediction~TSIR" = "seagreen",
  "prediction~TSIR-increased~S" = "purple",
  "prediction~TSIR-increased~beta" = "tomato",
  "TSIR~fit" = my_blue,
  "data" = "black"
)


circle_size <- 3
pl_d <-
  ggplot(data = dat, aes(x = true_cases, y = reported_cases)) + geom_point(
    aes(
      x = true_cases,
      y = reported_cases,
      color = pred_type,
      shape = year
    ),
    size = circle_size,
    stroke = 2
  ) +
  labs(shape = NULL, colour = NULL) +  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # aspect.ratio = 1,
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 14)
  ) +theme_bw() +
  theme(panel.grid = element_blank(),strip.background  =element_rect(fill="white"))+ 
# +geom_point(aes(x=true_cases, y=reported_cases,
#                fill=year), size=2, stroke = 0.05)
  ylim(0, 9000) + xlim(0, 9000) + scale_color_manual(values = colz, labels=label_parse()) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.position = c(.218, .745),
    legend.box.background = element_rect(colour = "black"),
    legend.background = element_blank(),
    axis.text = element_text(size = 14)
  ) + labs(x = "Reported Cases", y = "Predicted Cases")+theme(legend.spacing.y = unit(-0.1, "cm"))# + # space between legend rows
  # theme(legend.title.align = 1.08)














######################################################################################
######################################################################################
######################################################################################
############################# subplots together ######################################
######################################################################################
######################################################################################



pl<-plot_grid(pl_a, NULL, pl_b, pl_c, NULL, pl_d, ncol = 3, nrow = 2,align = "vh",labels = c("A","", "B","C","","D"), label_size =xy_size+20, hjust=-0.5, vjust = c(1.3,1.3,1.4,1.4),rel_widths = c(1, 0.08,1))+ theme(plot.margin = unit(c(1,1,2,1), "cm")) 




ggsave(paste0(homewd, "/final-figures/fig1.png"),
  plot = pl,
  device = NULL,
  path = NULL,
  scale = 1,
  width = 17,
  height = 12,
  bg='white',
  units = c("in"))









