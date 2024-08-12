rm(list=ls())

library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)


homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(homewd)

# plot age distribution by province,
# mean age by province 
# and FOI by province, as multipanels

#then, decide how to arrange all together

dat = read.csv(file = paste0(homewd, "/data/mort-dat.csv"), header=T, stringsAsFactors = F)
head(dat)


dat.melt <- melt(dat, id.vars = c("province", "year", "type"))
names(dat.melt)[4:5] <- c("month", "count")
unique(dat.melt$month)
dat.melt$month <- as.numeric(factor(dat.melt$month))
dat.melt$month<- as.character(dat.melt$month)
dat.melt$month[dat.melt$month!="12" & dat.melt$month!="11" & dat.melt$month!="10"] <- paste0("0", dat.melt$month[dat.melt$month!="12" & dat.melt$month!="11" & dat.melt$month!="10"])
dat.melt$year <- gsub(pattern=",", replacement = "", x= dat.melt$year, fixed = T)
dat.melt$year <- as.numeric(dat.melt$year)
dat.melt$date <- as.Date(paste0(dat.melt$year, "-", dat.melt$month, "-01"))
dat.melt$count <- as.numeric(dat.melt$count)
dat.melt$count[is.na(dat.melt$count)] <- 0

#total deaths
ggplot(subset(dat.melt, type=="deaths")) + geom_line(aes(x=date, y=count, color=province)) + facet_wrap(~province)

#and cfr
dat.cases = subset(dat.melt, type=="cases")
dat.deaths = subset(dat.melt, type=="deaths")
dat.deaths$cases <- dat.cases$count
dat.deaths$cfr <- dat.deaths$count/dat.cases$count
(sum(dat.deaths$count)/sum(dat.cases$count))*100 #0.5% CFR overall

ggplot(dat.deaths) + geom_line(aes(x=date, y=cfr, color=province)) + facet_wrap(~province) + coord_cartesian(ylim=c(0,.2))

#join cfr with dss data below

#quantify declines in the cfr through time by year
dat.year <- ddply(dat.deaths, .(province, year), summarise, tot_cases = sum(cases), tot_deaths = sum(count))
dat.year$cfr <- dat.year$tot_deaths/dat.year$tot_cases

head(dat.year)
dat.year$cfr <- as.numeric(dat.year$cfr)

unique(dat.year$province)

ggplot(dat.year) + geom_line(aes(x=year, y=cfr, color=province)) + facet_wrap(~province)



library(mgcv)

dat.year$province <- as.factor(dat.year$province)
dat.year$year <- as.numeric(dat.year$year)
dat.year$cfr <- as.numeric(dat.year$cfr)

dat.year <- dat.year[complete.cases(dat.year),]

head(dat.year)

m1 <- gam(tot_deaths~year + offset(log(tot_cases)) + s(province, bs="re"), data=dat.year, family = poisson(link="log"))
summary(m1) #significant declining cfr through time

#and also look at the proportion of diagnoses of different types for clinical cases
dat = read.csv(file = paste0(homewd, "/data/DENV-Dist-Diag.csv"), header=T, stringsAsFactors = F)
unique(dat$dianostic) #df, dhf, dss
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
head(dat)


p1 <- ggplot(subset(dat, provname=="Battambang")) + geom_jitter(aes(x=year, y=age, color=diagnostic), width=.1, size=.1, alpha=.3) + 
  facet_wrap(~provname) + theme_bw() + facet_grid(~diagnostic) +
  geom_vline(xintercept = 2007, linetype=2, size=.3, color="gray65") + 
  geom_vline(xintercept = 2012, linetype=2, size=.3, color="gray65") + 
  geom_vline(xintercept = 2019, linetype=2, size=.3, color="gray65") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        strip.background = element_rect(fill="white")) +
  geom_violin(aes(x=year,y=age, group=year, color = diagnostic),  position = position_dodge(),
              draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) 
p1



p2 <- ggplot(subset(dat, provname=="Phnom Penh")) + geom_jitter(aes(x=year, y=age, color=diagnostic), width=.1, size=.1, alpha=.3) + 
  facet_wrap(~provname) + theme_bw() + facet_grid(~diagnostic) +
  geom_vline(xintercept = 2007, linetype=2, size=.3, color="gray65") + 
  geom_vline(xintercept = 2012, linetype=2, size=.3, color="gray65") + 
  geom_vline(xintercept = 2019, linetype=2, size=.3, color="gray65") +
  theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
        strip.background = element_rect(fill="white")) +
  geom_violin(aes(x=year,y=age, group=year, color = diagnostic),  position = position_dodge(),
              draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) 
p2


dat$epimonth <- as.Date(dat$epimonth, format = "%m/%d/%y")

rate.df <- dlply(dat,.(epimonth, provname))

tot.sum <- function(df){
  tot.df = sum(df$case[df$diagnostic=="df"])
  tot.df[is.na(tot.df)] <- 0
  tot.dhf = sum(df$case[df$diagnostic=="dhf"])
  tot.dhf[is.na(tot.dhf)] <- 0
  tot.dss = sum(df$case[df$diagnostic=="dss"])
  tot.dss[is.na(tot.dss)] <- 0
  dhf.prop = tot.dhf/(tot.df + tot.dss +  tot.dhf)
  dss.prop = tot.dss/(tot.df + tot.dss +  tot.dhf)
  df.prop = tot.df/(tot.df + tot.dss +  tot.dhf)
  
  out.df = cbind.data.frame(provname=unique(df$provname),  year=unique(df$year),  epimonth=unique(df$epimonth), prop_dhf=dhf.prop, prop_dss=dss.prop, prop_df =  df.prop)
  return(out.df)
}

out.df <- data.table::rbindlist(lapply(rate.df, tot.sum))

head(out.df)

out.melt <- melt(out.df, id.vars = c("provname", "year", "epimonth"))
head(out.melt)
names(out.melt)[names(out.melt)=="variable"] <- "diagnostic"
out.melt$diagnostic  <- as.character(out.melt$diagnostic )
out.melt$diagnostic <- gsub(pattern = "prop_", replacement = "", x=out.melt$diagnostic, fixed = T)


ggplot(out.melt) + geom_line(aes(x=epimonth, y=value, color=diagnostic)) + 
                  facet_wrap(~provname) + scale_y_log10()


#and by year
tot.sum.year <- function(df){
  tot.df = sum(df$case[df$diagnostic=="df"])
  tot.df[is.na(tot.df)] <- 0
  tot.dhf = sum(df$case[df$diagnostic=="dhf"])
  tot.dhf[is.na(tot.dhf)] <- 0
  tot.dss = sum(df$case[df$diagnostic=="dss"])
  tot.dss[is.na(tot.dss)] <- 0
  dhf.prop = tot.dhf/(tot.df + tot.dss +  tot.dhf)
  dss.prop = tot.dss/(tot.df + tot.dss +  tot.dhf)
  df.prop = tot.df/(tot.df + tot.dss +  tot.dhf)
  tot_case = tot.df + tot.dss +  tot.dhf
  
  out.df = cbind.data.frame(provname=unique(df$provname),  year=unique(df$year),   prop_dhf=dhf.prop, prop_dss=dss.prop, prop_df =  df.prop, tot_case=tot_case)
  return(out.df)
}
rate.df.year <- dlply(dat,.(year, provname))

out.df.year <- data.table::rbindlist(lapply(rate.df.year, tot.sum.year))

head(out.df.year)

out.melt.year <- melt(out.df.year, id.vars = c("provname", "year", "tot_case"))
head(out.melt.year)
names(out.melt.year)[names(out.melt.year)=="variable"] <- "diagnostic"
out.melt.year$diagnostic  <- as.character(out.melt.year$diagnostic )
out.melt.year$diagnostic <- gsub(pattern = "prop_", replacement = "", x=out.melt.year$diagnostic, fixed = T)

#this one is the keeper:

head(out.melt.year)
dat.year$province <- as.character(dat.year$province)
out.dss = subset(out.melt.year, diagnostic=="dss")
setdiff(unique(out.dss$provname), unique(dat.year$province))
setdiff(unique(dat.year$province), unique(out.dss$provname))

dat.year$province[dat.year$province=="B.Meanchey"] <- "Banteay Meanchey"
dat.year$province[dat.year$province=="Kg.Chhnang"] <- "Kampong Chhnang"
dat.year$province[dat.year$province=="KampongSpeu"] <- "Kampong Speu"
dat.year$province[dat.year$province=="KampongThom"] <- "Kampong Thom"
dat.year$province[dat.year$province=="MondulKiri"] <- "Mondul Kiri"
dat.year$province[dat.year$province=="PreahVihear"] <- "Preah Vihear"
dat.year$province[dat.year$province=="PreyVeng"] <- "Prey Veng"
dat.year$province[dat.year$province=="RattanaKiri"] <- "Ratanak Kiri"
dat.year$province[dat.year$province=="K.Pr.Sihaknouk"] <- "Preah Sihanouk"
dat.year$province[dat.year$province=="StungTreng"] <- "Stung Treng"
dat.year$province[dat.year$province=="SvayRieng"] <- "Svay Rieng"
dat.year$province[dat.year$province=="O.Meanchey"] <- "Otdar Meanchey"
dat.year$province[dat.year$province=="Paillin"] <- "Pailin"

setdiff(unique(out.dss$provname), unique(dat.year$province))
setdiff(unique(dat.year$province), unique(out.dss$provname))

merge.cfr <- dplyr::select(dat.year, province, year, cfr )

names(out.dss)[names(out.dss)=="provname"] <- "province"
out.dss <- merge(out.dss, merge.cfr, by = c("province", "year"), all.x=T)
head(out.dss)
out.dss1 = dplyr::select(out.dss, province, year, value)
out.dss2 = dplyr::select(out.dss, province, year, cfr)
out.dss1$type <- "proportion-dss"
out.dss2$type <- "CFR"

names(out.dss1) [names(out.dss1)=="value"] <- "proportion"
names(out.dss2) [names(out.dss2)=="cfr"] <- "proportion"

out.dss <- rbind(out.dss1, out.dss2)

FigS21 <-  ggplot(out.dss) + geom_line(aes(x=year, y=proportion, color=type), size=1) + 
          ylab("proportion") + theme_bw() +
          theme(panel.grid = element_blank(), axis.title.x = element_blank(), 
                legend.title = element_blank(),
                legend.position = c(.87,.1),
                axis.title.y = element_text(size =16), axis.text = element_text(size=12),
                strip.background = element_rect(fill="white"))+
                facet_wrap(~province) #+ scale_y_log10()


ggsave(file = paste0(homewd, "/final-figures/FigS21.png"),
       plot= FigS21,
       units="mm",  
       width=110, 
       height=90, 
       scale=3, 
       dpi=300)


#include

#and do stats on year

p3 <- ggplot(subset(out.melt.year, diagnostic=="dss")) + geom_line(aes(x=year, y=value)) + ylab("proportion of cases diagnosed as DSS") +
      facet_wrap(~provname) + theme(axis.title.x = element_blank())#+ scale_y_log10()

p3

out.melt.year$provname <- as.factor(out.melt.year$provname)

#m2 <- lmer(value~year + (1|provname), data = subset(out.melt.year, diagnostic=="dss"))

#summary(m2)

#or, as gam
out.dss <- subset(out.melt.year, diagnostic=="dss")
out.dss$tot_dss <- round(out.dss$value*out.dss$tot_case,0)

m2 <- gam(tot_dss~year + s(provname, bs="re") + offset(log(tot_case)), data = out.dss, family = "poisson")
summary(m2) #sig neg decline


#What about predictors of the status?

head(dat)
dat$diag_dhf_dss <- 0
dat$diag_dhf_dss[dat$diagnostic=="dhf" | dat$diagnostic=="dss"] <- 1
dat$diag_dss <- 0
dat$diag_dss[dat$diagnostic=="dss"] <- 1

dat$provname <- as.factor(dat$provname)
dat$date <- as.Date(dat$date)
dat$month <- as.numeric(dat$month)
dat$year <- as.numeric(dat$year)
dat$sex[dat$sex=="m"] <- "M"
dat$sex <- as.factor(dat$sex)
dat$age <- as.numeric(dat$age)

m3 <- gam(diag_dhf_dss~s(year, k=7, bs="tp") +
               s(month, k=7, bs="cc") +
               s(provname, bs="re") +
               s(age, k=7, bs="tp") +
              s(sex, bs="re"), dat = dat, family = "binomial")
summary(m3)

source(paste0(homewd, "/figure-development/Fig1/mollentze-streicker-2020-functions.R"))

year.df <- get_partial_effects_continuous(gamFit=m3, var="year")
month.df <- get_partial_effects_continuous(gamFit=m3, var="month")
age.df <- get_partial_effects_continuous(gamFit=m3, var="age")
prov.df <- get_partial_effects(fit=m3, var="provname")
sex.df <- get_partial_effects(fit=m3, var="sex")

plot.partial(prov.df, var="provname", response_var = "severe dengue (dhf/dss)")
plot.partial(sex.df, var="sex", response_var = "severe dengue (dhf/dss)") #females more likely to be dhf/dss vs. df
plot.partial.cont(age.df, log = F, var="age", response_var = "severe dengue (dhf/dss)", alt_var = "age")
plot.partial.cont(month.df, log = F, var="month", response_var = "severe dengue (dhf/dss)", alt_var = "month")
plot.partial.cont(year.df, log = F, var="year", response_var = "severe dengue (dhf/dss)", alt_var = "year")

#and try just dss

m4 <- gam(diag_dss~s(year, k=7, bs="tp") +
            s(month, k=7, bs="cc") +
            s(provname, bs="re") +
            s(age, k=7, bs="tp") +
            s(sex, bs="re"), dat = dat, family = "binomial")
summary(m4)



year.df <- get_partial_effects_continuous(gamFit=m4, var="year")
month.df <- get_partial_effects_continuous(gamFit=m4, var="month")
age.df <- get_partial_effects_continuous(gamFit=m4, var="age")
prov.df <- get_partial_effects(fit=m4, var="provname")
sex.df <- get_partial_effects(fit=m4, var="sex")

plot.partial(prov.df, var="provname", response_var = "dengue shock syndrome")
plot.partial(sex.df, var="sex", response_var = "dengue shock syndrome") #females more likely to be dss vs. df
plot.partial.cont(age.df, log = F, var="age", response_var = "dengue shock syndrome", alt_var = "age") #wow! - way low in later ages
plot.partial.cont(month.df, log = F, var="month", response_var = "dengue shock syndrome", alt_var = "month")
plot.partial.cont(year.df, log = F, var="year", response_var = "dengue shock syndrome", alt_var = "year") # decreasing through time



#age distribution?

ggplot(dat) + geom_violin(aes(x=diagnostic,y=age)) + 
              geom_jitter(aes(x=diagnostic,y=age, color=year), width = .1, alpha=.3) +
              scale_color_continuous(trans="log10")
#and add to the dataset

# 
# head(out.melt.year)
# 
# ggplot(out.df) + geom_line(aes(x=year, y=prop_dhf, color=provname)) + facet_wrap(~provname)
# 
# ggplot(out.df) + geom_line(aes(x=year, y=prop_dss, color=provname)) + facet_wrap(~provname)
# 
# ggplot(out.df) + geom_line(aes(x=year, y=prop_df, color=provname)) + facet_wrap(~provname)
# 
# 
# out.melt <- melt(out.df, id.vars = c("provname", "year"))
# head(out.melt)
# names(out.melt)[3] <- "diagnostic"
# 
# ggplot(out.melt) + geom_line(aes(x=year, y=value, color=diagnostic)) + facet_wrap(~provname) + scale_y_log10()
# 
