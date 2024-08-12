

rm(list=ls())


library(ggplot2)
library(ggtree)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(rnaturalearth)
library(rnaturalearthdata)
library(treeio)
library(ape)
library(sp)
library(dismo)
library(geosphere)
library(rgeos)
library(sf)




#make three timetrees and one map all together
#read in tree from beast
#first, make map

#and load the metadata
#homewd= "/home/rstudio"
homewd = "/Users/carabrook/Developer/cambodia-dengue-national-public"
setwd(homewd)


#load the trees
tree1 <- read.beast(file = paste0(homewd, "/BEAST-tree/denv1-out-final/DENV1avg.tree"))
tree2 <- read.beast(file = paste0(homewd, "/BEAST-tree/denv2-out-final/DENV2avg.tree"))



tree1dat <- cbind.data.frame(tip_name = tree1@phylo$tip.label)
#tree1dat$beast_name <-tree1dat$tip_name
tree2dat <- cbind.data.frame(tip_name = tree2@phylo$tip.label)
#tree2dat$beast_name <-tree2dat$tip_name

head(tree1dat)
head(tree2dat)

#add serotype
tree1dat$DENV.serotype <- "DENV-1"
tree2dat$DENV.serotype <- "DENV-2"

#join
treedat <- rbind(tree1dat, tree2dat)

#now get the accession number and the date
treedat$date <- as.Date(sapply(strsplit(treedat$tip_name, "_"), function(x) x[[2]]))
treedat$accession_num <- sapply(strsplit(treedat$tip_name, "_"), function(x) x[[1]])

#and load the old data to get your country
genbank.dat <- read.csv(file = paste0(homewd, "/gather-sequences/All_Seq_SE_Asia.csv"), header = T, stringsAsFactors = F)
head(genbank.dat)
genbank.dat <- dplyr::select(genbank.dat, Accession, Country)
names(genbank.dat) <- c("accession_num", "country")

#and join
treedat <- merge(treedat, genbank.dat, by="accession_num")
head(treedat)
unique(treedat$country)

#and, until the new trees finish, swap the dates for the cambodia sequences for the 
# #precise dates here
# cam.dat <- read.csv(file = paste0(homewd, "/gather-sequences/cambodia-seq-details.csv"), header = T, stringsAsFactors = F)
# head(cam.dat)
# cam.dat <- dplyr::select(cam.dat, accession_num, date)
# cam.dat$date <- as.Date(cam.dat$date, format="%m/%d/%y")
# names(cam.dat) <- c("accession_num", "precise_date")
# 
# treedat <- merge(treedat, cam.dat, by="accession_num", all.x=T)
# head(treedat)
# 
# treedat$date[!is.na(treedat$precise_date)] <- treedat$precise_date[!is.na(treedat$precise_date)]

#and drop column 
# treedat <- dplyr::select(treedat, -(precise_date))
treedat$date

#now, for all the cambodia sequences, add in lat-long info and the DENV genotype
gis.dat <- read.csv(file = paste0(homewd, "/data/lat-long-cambodia.csv"), header = T, stringsAsFactors = F)
head(gis.dat)
gis.dat <- dplyr::select(gis.dat, -(DENV.serotype), -(country), -(date), -(year))
length(unique(gis.dat$accession_num)) #260

gis.dat[duplicated(gis.dat$accession_num),]

treedat <- merge(treedat, gis.dat, by = "accession_num", all.x= T)
head(treedat)

#sub.dat <- subset(treedat, country == "Cambodia" & is.na(lat))
#write.csv(sub.dat, file = "fill_cambodia.csv", row.names = F)
#and, finally, for the Kampong-Speu provinces, add in the age and sex of the host

#load host dat
host.dat <- read.csv(file = paste0(homewd, "/data/IDseq_PAGODAS_ALL_metadata_through_2020_CLEAN.csv"), header = T, stringsAsFactors = F)
head(host.dat)
host.dat <- dplyr::select(host.dat, NIH.ID, age, sex, date)
host.dat$date <- as.Date(host.dat$date, format = "%m/%d/%y")

#treemerge <- merge(treedat, host.dat, by = c("NIH.ID"), all.x=T)
treemerge <- merge(treedat, host.dat, by = c("NIH.ID", "date"), all.x=T)
head(treemerge)

write.csv(treemerge, file = paste0(homewd, "/data/beasttree_metadata.csv"), row.names = F)
#above is 95% correct - the "tip_name" column will change in the new version

#now go make transmission trees

#and build figures 4 and 5

