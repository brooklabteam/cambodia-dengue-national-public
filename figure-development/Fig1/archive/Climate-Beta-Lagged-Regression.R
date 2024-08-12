rm(list=ls())


library(ggplot2)
library(plyr)
library(dplyr)
library(lme4)
library(lmerTest)
library(sjPlot)
library(mgcv)


homewd = "/Users/carabrook/Developer/cambodia-dengue-national"


# compare non-parametric binned regression (using GAMs) with a linear regression model
# of how transmission (beta) responds to climate (temp + precipitation)


# first load the province data
dat <- read.csv(file = paste0(homewd, "/data/lagged-prov-clim-beta.csv"), header = T, stringsAsFactors = F)
head(dat)

#sort again by time
dat <- arrange(dat, provname, time)


#remove the epidemic years with the betas that were not fitted
dat.epi = subset(dat, year==2007 | year==2012 | year==2019) #set these aside
dat.epi$beta <- NA
dat.epi$betalow <- NA
dat.epi$betahigh <- NA

dat = subset(dat, year!=2007 & year!=2012 & year!=2019)

# #and national data
# nat.dat <- read.csv(file = paste0(homewd, "/data/lagged-nat-clim-beta.csv"), header = T, stringsAsFactors = F)
# head(nat.dat)

# log beta for the response
dat$log_beta <- log((dat$beta))



# and make random effects into factors for regression
dat$provname <- as.factor(dat$provname)
dat$biweek <- as.factor(dat$biweek)
dat$year <- as.factor(dat$year)
dat$epiyr <- as.factor(dat$epiyr)

briere_out <- function(temp_c, c_val, minT, maxT){
  briere_out = c_val*temp_c*(temp_c-minT)*((maxT-temp_c)^(1/2))
  
  return(briere_out)
}

dat$temp_C_lag_briere <- briere_out(temp_c = dat$temp_C_lag, c_val = 4.99*10^-4, minT = 25.27, maxT = 38.63)



m1 <- lmer(log_beta~temp_C_lag + precip_mm_lag + (1 | provname / biweek) + (1 | provname/year) + (1 | epiyr), 
           data = dat, control = lmerControl(optimizer ="Nelder_Mead"))


m1b <- lmer(log_beta~temp_C_lag + precip_mm_lag + (1 | provname / biweek)  + (1 | epiyr), 
           data = dat, control = lmerControl(optimizer ="Nelder_Mead"))


summary(m1)


# Now, plot coefficients on beta (for this model they are constant across x values)
est.df <- get_model_data(m1, type="est")
est.df$term <- as.character(est.df$term)
names(est.df)[names(est.df)=="term"] <- "predictor"
est.df$predictor[est.df$predictor=="temp_C_lag"] <- "'lagged temp ('^0~'C)'"
est.df$predictor[est.df$predictor=="precip_mm_lag"] <- "'lagged precip (mm)'"

est.df$posneg <- "pos"
est.df$posneg[est.df$estimate<0 & est.df$conf.high<0] <- "neg"
est.df$posneg[est.df$estimate<0 & est.df$conf.high>0] <- "notsig"
est.df$posneg[est.df$estimate>0 & est.df$conf.low<0] <- "notsig"
est.df$p.stars[est.df$p.stars==""] <- "."

colz <- c('pos' = "red3", 'neg'="cornflowerblue", 'notsig'="gray50")


pA <- ggplot(data=est.df) + coord_cartesian(ylim=c(-.007,0.016)) + scale_color_manual(values=colz) +
  geom_label(aes(x=predictor, y=0.0158, label=p.stars), label.size = 0, size=8) + 
  facet_grid(~predictor, labeller = label_parsed,  switch = "x", scales = "free") +
  geom_point(aes(x=predictor, y = estimate, color=posneg), show.legend = F, size=3) +
  ylab(bquote('coefficients on'~beta)) + geom_hline(aes(yintercept=0), linetype=2)+ theme_bw() +
  geom_linerange(aes(x=predictor, ymin=conf.low, ymax=conf.high, color=posneg), show.legend = F, size=1) +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 16), axis.ticks.x=element_blank(),
        strip.background = element_blank(), axis.text = element_text(siz=16), strip.text = element_text(size=16),
        axis.title.x = element_blank(), axis.text.y = element_text(size = 13), axis.text.x = element_blank(),
        plot.margin = unit(c(.2,1.2,.2,.2), "cm"), panel.spacing = unit(c(0), "cm")) + 
  scale_y_continuous(breaks=c(-0.005,0, 0.005, 0.01))

pA
# then, also do a linear regression for precipitation paired with a Briere response for temperature


m2 <- lmer(log_beta~temp_C_lag_briere + 
             precip_mm_lag + (1 | provname / biweek) + (1 | provname/year) + (1 | epiyr), 
           data = dat, control = lmerControl(optimizer ="Nelder_Mead"))# 'optimx', optCtrl=list(method='nlminb')))


summary(m2)

plot_model(m2, type="est")

est.df <- get_model_data(m2, type="est")
est.df$term <- as.character(est.df$term)
names(est.df)[names(est.df)=="term"] <- "predictor"
est.df$predictor[est.df$predictor=="temp_C_lag_briere"] <- "'lagged temp ('^0~'C)'"
est.df$predictor[est.df$predictor=="precip_mm_lag"] <- "atop(' ','lagged precip (mm)')"

est.df$posneg <- "pos"
est.df$posneg[est.df$estimate<0 & est.df$conf.high<0] <- "neg"
est.df$posneg[est.df$estimate<0 & est.df$conf.high>0] <- "notsig"
est.df$posneg[est.df$estimate>0 & est.df$conf.low<0] <- "notsig"
est.df$p.stars[est.df$p.stars==""] <- "."

# now, multiply the coefficient by the briere output to get the actual value
# so 


ex.dat = cbind.data.frame(temp_C=20:40, Briere_funct=briere_out(20:40, c_val = 4.99 * 10^-4, minT = 25.27, maxT = 38.63))
ex.dat$Briere_funct[is.na(ex.dat$Briere_funct)]<- 0
#ex.dat$Briere_funct[ex.dat$Briere_funct<0]<- 0
ex.dat$beta_coeff <- est.df$estimate[est.df$predictor=="'lagged temp ('^0~'C)'"]*ex.dat$Briere_funct
ex.dat$beta_coeff_low <- est.df$conf.low[est.df$predictor=="'lagged temp ('^0~'C)'"]*ex.dat$Briere_funct
ex.dat$beta_coeff_high <- est.df$conf.high[est.df$predictor=="'lagged temp ('^0~'C)'"]*ex.dat$Briere_funct

#plot

p2a <- ggplot(data=subset(est.df, predictor=="atop(' ','lagged precip (mm)')")) +  
  scale_color_manual(values=colz) + coord_cartesian(ylim=c(-.007,0.016)) +
  geom_label(aes(x=predictor, y=0.0155, label=p.stars), label.size = 0, size=8) + 
  facet_grid(~predictor, labeller = label_parsed,  switch = "x", scales = "free") +
  geom_point(aes(x=predictor, y = estimate, color=posneg), show.legend = F, size=3) +
  ylab(bquote('coefficients on'~beta)) + geom_hline(aes(yintercept=0), linetype=2)+ theme_bw() +
  geom_linerange(aes(x=predictor, ymin=conf.low, ymax=conf.high, color=posneg), show.legend = F, size=1) +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 16), axis.ticks.x=element_blank(),
        strip.background = element_blank(),  strip.text = element_text(size=16),
        axis.title.x = element_blank(), axis.text.y = element_text(size = 13), axis.text.x = element_blank(),
        plot.margin = unit(c(.2,0,.4,.6), "cm"), panel.spacing = unit(c(0), "cm")) 


#and the Briere
ex.dat$predictor <- "'lagged temp ('^0~'C)'"
p2b <- ggplot(ex.dat) + 
        geom_line(aes(x=temp_C, y=beta_coeff), size=1, color="red3") + 
        geom_label(aes(x=30, y=0.128, label="***"), label.size = 0, size=8) + 
        coord_cartesian(ylim = c(0,.13)) + #ylab(bquote('coefficients on'~beta)) +
        facet_grid(~predictor, labeller = label_parsed,  switch = "x", scales = "free") +
        geom_ribbon(aes(x=temp_C, ymin=beta_coeff_low, ymax=beta_coeff_high), alpha=.3, fill="red3") + 
        geom_hline(aes(yintercept=0), linetype=2)+ theme_bw() + scale_y_continuous(position="right")+
        theme(panel.grid = element_blank(), axis.title.y = element_text(size = 16), 
        strip.background = element_blank(), strip.text = element_text(size=16),
        strip.placement = "outside", axis.title.y.right = element_blank(),
        axis.title = element_blank(), axis.text = element_text(size = 13), 
        plot.margin = unit(c(.2,.2,.1,0), "cm"), panel.spacing = unit(c(0), "cm")) 


pB <- cowplot::plot_grid(p2a,p2b, nrow=1, ncol=2, rel_widths = c(1.1,1))


# Now try a GAM (no need for Biere since this is already captured in the flexible model approach)

library(mgcv)
gam1 <- gam(log_beta~s(temp_C_lag, k=7, bs='tp') + 
                    s(precip_mm_lag, k=7, bs='tp') +
                    s(provname, biweek, bs="re") +
                    s(provname, year, bs="re") +
                    s(epiyr, bs="re"), 
                    data = dat)

#linear response GAM
summary(gam1)# R-sq.(adj) =  0.994.  Deviance explained = 99.5%

test.df = cbind.data.frame(temp_C_lag=20:34)
test.df$precip_mm_lag <- 20
test.df$biweek <- 10
test.df$year <- 2005
test.df$provname = "Kandal"
test.df$epiyr = 2007

temp.df <- cbind.data.frame(x=20:34, beta_coeff=c(predict.gam(gam1, newdata = test.df, exclude = c("s(precip_mm_lag)", "s(provname,biweek)", "s(provname,year)", "s(epiyr)"), type = "terms")))
temp.df$beta_coeff_low <- c(predict.gam(gam1, newdata = test.df, exclude = c("s(precip_mm_lag)", "s(provname,biweek)", "s(provname,year)",  "s(epiyr)"), type = "terms")-1.96*predict.gam(gam1, newdata = test.df, exclude = c("s(precip_mm_lag)", "s(provname,biweek)", "s(provname,year)", "s(epiyr)"), type = "terms", se.fit = TRUE)$se.fit)
temp.df$beta_coeff_high <-c(predict.gam(gam1, newdata = test.df, exclude = c("s(precip_mm_lag)", "s(provname,biweek)", "s(provname,year)",  "s(epiyr)"), type = "terms")+1.96*predict.gam(gam1, newdata = test.df, exclude = c("s(precip_mm_lag)", "s(provname,biweek)", "s(provname,year)",  "s(epiyr)"), type = "terms", se.fit = TRUE)$se.fit)
temp.df$predictor <- "'lagged temp ('^0~'C)'"



test.df = cbind.data.frame(precip_mm_lag = 0:30)
test.df$temp_C_lag = 25
test.df$biweek <- 10
test.df$year <- 2005
test.df$provname = "Kandal"
test.df$epiyr = 2007


precip.df <- cbind.data.frame(x=0:30, beta_coeff=c(predict.gam(gam1, newdata = test.df, exclude = c("s(temp_C_lag)", "s(provname,biweek)", "s(provname,year)",  "s(epiyr)"), type = "terms")))
precip.df$beta_coeff_low <- c(predict.gam(gam1, newdata = test.df, exclude = c("s(temp_C_lag)", "s(provname,biweek)", "s(provname,year)",  "s(epiyr)"), type = "terms")-1.96*predict.gam(gam1, newdata = test.df, exclude = c("s(temp_C_lag)", "s(provname,biweek)", "s(provname,year)",  "s(epiyr)"), type = "terms", se.fit = TRUE)$se.fit)
precip.df$beta_coeff_high <-c(predict.gam(gam1, newdata = test.df, exclude = c("s(temp_C_lag)", "s(provname,biweek)", "s(provname,year)",  "s(epiyr)"), type = "terms")+1.96*predict.gam(gam1, newdata = test.df, exclude = c("s(temp_C_lag)", "s(provname,biweek)", "s(provname,year)",  "s(epiyr)"), type = "terms", se.fit = TRUE)$se.fit)
precip.df$predictor <- "'lagged precip (mm)'"


gam.df.linear <- rbind(temp.df, precip.df)

star.df <- cbind.data.frame(x=c(15,27), y=c(0.145,0.145), label=c("***", "**"), predictor = c("'lagged precip (mm)'", "'lagged temp ('^0~'C)'"))


pC <- ggplot(data = gam.df.linear) + 
  geom_line(aes(x=x, y=beta_coeff), size=1, color="red3") + ylab(bquote('coefficients on'~beta)) +
  geom_label(data=star.df, aes(x=x,y=y,label=label), label.size = 0, size=8) +
  geom_ribbon(aes(x=x, ymin=beta_coeff_low, ymax=beta_coeff_high), alpha=.3, fill="red3") + theme_bw() + 
  facet_grid(~predictor, scales = "free_x", labeller = label_parsed, switch = "x") + 
  geom_hline(aes(yintercept=0), linetype=2)+ 
  theme(panel.grid = element_blank(), strip.text = element_text(size=16), 
        panel.spacing = unit(c(0), "cm"),
        plot.margin = unit(c(.2,1.2,.2,.5), "cm"),
        strip.placement = "outside", strip.background = element_blank(),
        axis.title.x = element_blank(), axis.text = element_text(size=13), 
        axis.title.y = element_text(size=16))


FigSupp <- cowplot::plot_grid(pA, pB, pC, nrow=3, ncol=1, labels=c("A", "B", "C"), label_size=22, rel_heights = c(1,1.1,1.1))


ggsave(file = paste0(homewd, "/final-figures/SuppFig_BetaClimResponse.png"),
       plot = FigSupp,
       units="mm",  
       width=70, 
       height=100, 
       scale=3, 
       dpi=300)



# Now also do the prediction of the beta 
# for the epidemic years (based on their climate values)
# from this fitted model GAM

head(dat.epi)
dat.epi$provname <- as.factor(dat.epi$provname)
dat.epi$biweek <- as.factor(dat.epi$biweek)
dat.epi$year <- as.factor(dat.epi$year)
dat.epi$epiyr <- as.factor(dat.epi$epiyr)

dat.epi$beta <- c(exp(predict.gam(gam1, newdata = dat.epi, type="response")))
dat.epi$betalow <- c(exp(predict.gam(gam1, newdata = dat.epi, type="response") - 1.96*(predict.gam(gam1, newdata = dat.epi, type="response", se.fit = T)$se.fit)))
dat.epi$betahigh <- c(exp(predict.gam(gam1, newdata = dat.epi, type="response") - 1.96*(predict.gam(gam1, newdata = dat.epi, type="response", se.fit = T)$se.fit)))


#also try with the linear regression
dat.epi2 = dat.epi
dat.epi2$beta <- c(exp(predict(m1b, newdata = dat.epi2, type="response"))) #can't predict CIs from mixed effects model
dat.epi2$betalow <-NA
dat.epi2$betahigh <-NA

# attach to the main dataset
head(dat)
dat <- dplyr::select(dat, -(temp_C_lag_briere), -(log_beta))

dat$year <- as.numeric(as.character(dat$year))
dat$biweek <- as.numeric(as.character(dat$biweek))

dat.all <- rbind(dat, dat.epi)
dat.all2 <- rbind(dat, dat.epi2)


dat.all <- arrange(dat.all, provname, year, biweek)
dat.all2 <- arrange(dat.all2, provname, year, biweek)

# Now pull in the pop and births, merge, and send to TSIR
tsir.prov <- read.csv(file=paste0(homewd, "/data/tsir_dat_province.csv"), header = T, stringsAsFactors = F)
head(tsir.prov)

dat.all <- dplyr::select(dat.all, -(time), -(cases))
dat.all2 <- dplyr::select(dat.all2, -(time), -(cases))

tsir.clim <- merge(dat.all, tsir.prov, by=c("provname", "year", "biweek"), all.x = T)
head(tsir.clim)

write.csv(tsir.clim, paste0(homewd, "/data/tsir_dat_beta_climate_province.csv"), row.names = F)

tsir.clim2 <- merge(dat.all2, tsir.prov, by=c("provname", "year", "biweek"), all.x = T)
head(tsir.clim2)

write.csv(tsir.clim2, paste0(homewd, "/data/tsir_dat_beta_climate_province_linear.csv"), row.names = F)

############################################################
############################################################

rm(list=ls())

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"


# Also repeat at the national level
dat <- read.csv(file = paste0(homewd, "/data/lagged-nat-clim-beta.csv"), header = T, stringsAsFactors = F)
head(dat)

#remove the epidemic years that were not fitted
dat.epi = subset(dat, year==2007 | year==2012 | year==2019) #set these aside
dat.epi$beta <- NA
dat.epi$betalow <- NA
dat.epi$betahigh <- NA

dat = subset(dat, year!=2007 & year!=2012 & year!=2019)

# log beta for the response
dat$log_beta <- log((dat$beta))

# and make random effects into factors for regressiion
unique(dat$provname) #just national
dat$biweek <- as.factor(dat$biweek)
dat$year <- as.factor(dat$year)

briere_out <- function(temp_c, c_val, minT, maxT){
  briere_out = c_val*temp_c*(temp_c-minT)*((maxT-temp_c)^(1/2))
  
  return(briere_out)
}

dat$temp_C_lag_briere <- briere_out(temp_c = dat$temp_C_lag, c_val = 4.99*10^-4, minT = 25.27, maxT = 38.63)



m1 <- lmer(log_beta~temp_C_lag + precip_mm_lag + (1 |  biweek) + (1 | year), 
           data = dat)# control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))

summary(m1) # not very sig


# Now, plot coefficients on beta (for this model they are constant across x values)
est.df <- get_model_data(m1, type="est")
est.df$term <- as.character(est.df$term)
names(est.df)[names(est.df)=="term"] <- "predictor"
est.df$predictor[est.df$predictor=="temp_C_lag"] <- "'lagged temp ('^0~'C)'"
est.df$predictor[est.df$predictor=="precip_mm_lag"] <- "'lagged precip (mm)'"

est.df$posneg <- "pos"
est.df$posneg[est.df$estimate<0 & est.df$conf.high<0] <- "neg"
est.df$posneg[est.df$estimate<0 & est.df$conf.high>0] <- "notsig"
est.df$posneg[est.df$estimate>0 & est.df$conf.low<0] <- "notsig"
est.df$p.stars[est.df$p.stars==""] <- "."

colz <- c('pos' = "red3", 'neg'="cornflowerblue", 'notsig'="gray50")


pA <- ggplot(data=est.df) + scale_color_manual(values=colz) +
  geom_label(aes(x=predictor, y=0.0158, label=p.stars), label.size = 0, size=8) + 
  facet_grid(~predictor, labeller = label_parsed,  switch = "x", scales = "free") +
  geom_point(aes(x=predictor, y = estimate, color=posneg), show.legend = F, size=3) +
  ylab(bquote('coefficients on'~beta)) + geom_hline(aes(yintercept=0), linetype=2)+ theme_bw() +
  geom_linerange(aes(x=predictor, ymin=conf.low, ymax=conf.high, color=posneg), show.legend = F, size=1) +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 16), axis.ticks.x=element_blank(),
        strip.background = element_blank(), axis.text = element_text(siz=16), strip.text = element_text(size=16),
        axis.title.x = element_blank(), axis.text.y = element_text(size = 13), axis.text.x = element_blank(),
        plot.margin = unit(c(.2,1.2,.2,.2), "cm"), panel.spacing = unit(c(0), "cm")) + 
  scale_y_continuous(breaks=c(-0.01,-0.005,0, 0.005, 0.01, 0.015))

pA
# then, also do a linear regression for precipitation paired with a Briere response for temperature


m2 <- lmer(log_beta~temp_C_lag_briere + 
             precip_mm_lag + (1 |  biweek) + (1 | year), 
           data = dat)#, control = lmerControl(optimizer ="Nelder_Mead"))# 'optimx', optCtrl=list(method='nlminb')))


summary(m2)

plot_model(m2, type="est")

est.df <- get_model_data(m2, type="est")
est.df$term <- as.character(est.df$term)
names(est.df)[names(est.df)=="term"] <- "predictor"
est.df$predictor[est.df$predictor=="temp_C_lag_briere"] <- "'lagged temp ('^0~'C)'"
est.df$predictor[est.df$predictor=="precip_mm_lag"] <- "atop(' ','lagged precip (mm)')"

est.df$posneg <- "pos"
est.df$posneg[est.df$estimate<0 & est.df$conf.high<0] <- "neg"
est.df$posneg[est.df$estimate<0 & est.df$conf.high>0] <- "notsig"
est.df$posneg[est.df$estimate>0 & est.df$conf.low<0] <- "notsig"
est.df$p.stars[est.df$p.stars==""] <- "."

# now, multiply the coefficient by the briere output to get the actual value
# so 


ex.dat = cbind.data.frame(temp_C=20:40, Briere_funct=briere_out(20:40, c_val = 4.99 * 10^-4, minT = 25.27, maxT = 38.63))
ex.dat$Briere_funct[is.na(ex.dat$Briere_funct)]<- 0
#ex.dat$Briere_funct[ex.dat$Briere_funct<0]<- 0
ex.dat$beta_coeff <- est.df$estimate[est.df$predictor=="'lagged temp ('^0~'C)'"]*ex.dat$Briere_funct
ex.dat$beta_coeff_low <- est.df$conf.low[est.df$predictor=="'lagged temp ('^0~'C)'"]*ex.dat$Briere_funct
ex.dat$beta_coeff_high <- est.df$conf.high[est.df$predictor=="'lagged temp ('^0~'C)'"]*ex.dat$Briere_funct

#plot

p2a <- ggplot(data=subset(est.df, predictor=="atop(' ','lagged precip (mm)')")) +  
  scale_color_manual(values=colz) + coord_cartesian(ylim=c(-.00025,0.0005)) +
  geom_label(aes(x=predictor, y=0.0155, label=p.stars), label.size = 0, size=8) + 
  facet_grid(~predictor, labeller = label_parsed,  switch = "x", scales = "free") +
  geom_point(aes(x=predictor, y = estimate, color=posneg), show.legend = F, size=3) +
  ylab(bquote('coefficients on'~beta)) + geom_hline(aes(yintercept=0), linetype=2)+ theme_bw() +
  geom_linerange(aes(x=predictor, ymin=conf.low, ymax=conf.high, color=posneg), show.legend = F, size=1) +
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 16), axis.ticks.x=element_blank(),
        strip.background = element_blank(),  strip.text = element_text(size=16),
        axis.title.x = element_blank(), axis.text.y = element_text(size = 13), axis.text.x = element_blank(),
        plot.margin = unit(c(.2,0,.4,.6), "cm"), panel.spacing = unit(c(0), "cm")) 


#and the Briere
ex.dat$predictor <- "'lagged temp ('^0~'C)'"
p2b <- ggplot(ex.dat) + 
  geom_line(aes(x=temp_C, y=beta_coeff), size=1) + 
  #geom_label(aes(x=30, y=0.128, label="***"), label.size = 0, size=8) + 
  coord_cartesian(ylim = c(0,.13)) + #ylab(bquote('coefficients on'~beta)) +
  facet_grid(~predictor, labeller = label_parsed,  switch = "x", scales = "free") +
  geom_ribbon(aes(x=temp_C, ymin=beta_coeff_low, ymax=beta_coeff_high), alpha=.3) + 
  geom_hline(aes(yintercept=0), linetype=2)+ theme_bw() + scale_y_continuous(position="right")+
  theme(panel.grid = element_blank(), axis.title.y = element_text(size = 16), 
        strip.background = element_blank(), strip.text = element_text(size=16),
        strip.placement = "outside", axis.title.y.right = element_blank(),
        axis.title = element_blank(), axis.text = element_text(size = 13), 
        plot.margin = unit(c(.2,.2,.1,0), "cm"), panel.spacing = unit(c(0), "cm")) 


pB <- cowplot::plot_grid(p2a,p2b, nrow=1, ncol=2, rel_widths = c(1.1,1))


# Now try a GAM (no need for Biere since this is already captured in the flexible model approach)

library(mgcv)
gam1 <- gam(log_beta~s(temp_C_lag, k=7, bs='tp') + 
              s(precip_mm_lag, k=7, bs='tp') +
              s( biweek, bs="re") +
              s( year, bs="re"), 
            data = dat)# control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))

#linear response GAM
summary(gam1)# R-sq.(adj) =  0.994.  Deviance explained = 99.5%

test.df = cbind.data.frame(temp_C_lag=20:34)
test.df$precip_mm_lag <- 20
test.df$biweek <- 10
test.df$year <- 2005
test.df$provname = "national"

temp.df <- cbind.data.frame(x=20:34, beta_coeff=c(predict.gam(gam1, newdata = test.df, exclude = c("s(precip_mm_lag)", "s(biweek)", "s(year)"), type = "terms")))
temp.df$beta_coeff_low <- c(predict.gam(gam1, newdata = test.df, exclude = c("s(precip_mm_lag)", "s(biweek)", "s(year)"), type = "terms")-1.96*predict.gam(gam1, newdata = test.df, exclude = c("s(precip_mm_lag)", "s(biweek)", "s(year)"), type = "terms", se.fit = TRUE)$se.fit)
temp.df$beta_coeff_high <-c(predict.gam(gam1, newdata = test.df, exclude = c("s(precip_mm_lag)", "s(biweek)", "s(year)"), type = "terms")+1.96*predict.gam(gam1, newdata = test.df, exclude = c("s(precip_mm_lag)", "s(biweek)", "s(year)"), type = "terms", se.fit = TRUE)$se.fit)
temp.df$predictor <- "'lagged temp ('^0~'C)'"



test.df = cbind.data.frame(precip_mm_lag = 0:30)
test.df$temp_C_lag = 25
test.df$biweek <- 10
test.df$year <- 2005
test.df$provname = "national"


precip.df <- cbind.data.frame(x=0:30, beta_coeff=c(predict.gam(gam1, newdata = test.df, exclude = c("s(temp_C_lag)", "s(biweek)", "s(year)"), type = "terms")))
precip.df$beta_coeff_low <- c(predict.gam(gam1, newdata = test.df, exclude = c("s(temp_C_lag)", "s(biweek)", "s(year)"), type = "terms")-1.96*predict.gam(gam1, newdata = test.df, exclude = c("s(temp_C_lag)", "s(biweek)", "s(year)"), type = "terms", se.fit = TRUE)$se.fit)
precip.df$beta_coeff_high <-c(predict.gam(gam1, newdata = test.df, exclude = c("s(temp_C_lag)", "s(biweek)", "s(year)"), type = "terms")+1.96*predict.gam(gam1, newdata = test.df, exclude = c("s(temp_C_lag)", "s(biweek)", "s(year)"), type = "terms", se.fit = TRUE)$se.fit)
precip.df$predictor <- "'lagged precip (mm)'"


gam.df.linear <- rbind(temp.df, precip.df)

#star.df <- cbind.data.frame(x=c(15,27), y=c(0.145,0.145), label=c("***", "**"), predictor = c("'lagged precip (mm)'", "'lagged temp ('^0~'C)'"))


pC <- ggplot(data = gam.df.linear) + 
  geom_line(aes(x=x, y=beta_coeff), size=1) + ylab(bquote('coefficients on'~beta)) +
  #geom_label(data=star.df, aes(x=x,y=y,label=label), label.size = 0, size=8) +
  geom_ribbon(aes(x=x, ymin=beta_coeff_low, ymax=beta_coeff_high), alpha=.3) + theme_bw() + 
  facet_grid(~predictor, scales = "free_x", labeller = label_parsed, switch = "x") + 
  geom_hline(aes(yintercept=0), linetype=2)+ 
  theme(panel.grid = element_blank(), strip.text = element_text(size=16), 
        panel.spacing = unit(c(0), "cm"),
        plot.margin = unit(c(.2,1.2,.2,.5), "cm"),
        strip.placement = "outside", strip.background = element_blank(),
        axis.title.x = element_blank(), axis.text = element_text(size=13), 
        axis.title.y = element_text(size=16))


FigSuppNat <- cowplot::plot_grid(pA, pB, pC, nrow=3, ncol=1, labels=c("A", "B", "C"), label_size=22, rel_heights = c(1,1.1,1.1))


print(FigSuppNat)

#none of these are significant, so let's stop here
