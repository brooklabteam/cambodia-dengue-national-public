rm(list=ls())

#time to make Fig3A

library(ggplot2)
library(ggtree)
library(ape)
library(ggnewscale)

homewd= "/Users/carabrook/Developer/cambodia-dengue-national/"

setwd(paste0(homewd, "/figure-development/FigS18"))

#load the DENV1 tree
treeD1 <-  read.tree(file = paste0(homewd, "figure-development/FigS18/raxml-out/denv1ML/T3.raxml.supportFBP"))

#and DENV2
treeD2 <-  read.tree(file = paste0(homewd, "figure-development/FigS18/raxml-out/denv2ML/T3.raxml.supportFBP"))
#root it


rooted.D1 <- root(treeD1, which(treeD1$tip.label == "NC_002640_DENV4"))
#take a quick look in base R
plot(rooted.D1)

#and DENV2
rooted.D2 <- root(treeD2, which(treeD2$tip.label == "NC_002640_DENV4"))
#take a quick look in base R
plot(rooted.D2)

#load tree data prepared from elsewhere
dat <-read.csv(file = paste0(homewd, "figure-development/FigS18/ML-Sequences.csv"))
head(dat)


#check subgroup names
unique(dat$Subclade[dat$Serotype=="DENV-1"])
unique(dat$Subclade[dat$Serotype=="DENV-2"])

#add denv4
dat.add <- c("NC_002640", NA, NA, "DENV-4", "outgroup")
dat <- rbind(dat, dat.add)

dat.denv1 <- subset(dat, Serotype!="DENV-2")

tree.denv1 <- cbind.data.frame(tip_label = rooted.D1$tip.label)
tree.denv1$Accession <- sapply(strsplit(tree.denv1$tip_label, split = "_"), function(x) x[[1]])

dat.all.denv1 <- merge(tree.denv1, dat.denv1, by ="Accession", all.x = T, sort=F)
dat.all.denv1$Accession[dat.all.denv1$tip_label=="NC_002640_DENV4"] <- "NC_002640"
dat.all.denv1$Serotype[dat.all.denv1$tip_label=="NC_002640_DENV4"] <- "DENV-4"

dat.all.denv1 <- dplyr::select(dat.all.denv1, tip_label, Accession, Locality, Serotype, Subclade, Collection_Year)
# unique(dat.all.denv1$Subclade)
# subset(dat.all.denv1, is.na(Subclade))
dat.all.denv1$Subclade[dat.all.denv1$Serotype=="DENV-4"] <- "Outgroup"

colz1 = c('Genotype-I' = "mediumseagreen",
          'Genotype-II' = "royalblue", 
          'Genotype-III' ="darkgoldenrod1",
          'Old-Cambodia'= "palevioletred",
          'New-Cambodia' = "mediumpurple4", 
          'Outgroup' = "black")



#take a glance
p1 <- ggtree(rooted.D1) %<+% dat.all.denv1 + geom_tippoint(aes(fill=Subclade), shape=21 ) +
  geom_tiplab(size=1) + geom_nodelab(size=1) +
  scale_fill_manual(values=colz1) + theme(legend.position = c(.2,.85), 
                                          legend.title = element_blank())

p1 #looks great - change to clade label and collapse the outgroup

#now shrink the long branch length for the outgroup
sort(rooted.D1$edge.length) #.03 is the next longest
rooted.D1$edge.length[rooted.D1$edge.length==max(rooted.D1$edge.length)] <- .15

#and add cladebars
#here for genotype II
gen2 <- MRCA(rooted.D1, which(rooted.D1$tip.label == "DQ285561_Seychelles_2004" ),which(rooted.D1$tip.label == "AB204803_Micronesia_2004"))
gen3 <- MRCA(rooted.D1, which(rooted.D1$tip.label == "EU081258_Singapore_2005" ),which(rooted.D1$tip.label == "EF122231_French_Guiana_2006"))
gen1 <- MRCA(rooted.D1, which(rooted.D1$tip.label == "AY726550_Myanmar_2001" ),which(rooted.D1$tip.label == "OK159947_Cambodia_2019"))


dat.all.denv1$Subclade
#dat.all.denv1$Subclade[is.na(dat.all.denv1$Subclade)] <- "Outgroup"
dat.all.denv1$Subclade <- factor(dat.all.denv1$Subclade, 
                                 levels=c("Genotype-I", "Genotype-II", "Genotype-III",
                                          "Old-Cambodia", "New-Cambodia", "Outgroup")) 


# #and make new tip labels
# names(dat.all.denv1)[names(dat.all.denv1)=="tip_label"] <- "old_tip_label"
# dat.all.denv1$year <- sapply(strsplit(dat.all.denv1$old_tip_label, "_"),  function(x) x[[length(x)]] )
# dat.all.denv1$New_Accession[is.na(dat.all.denv1$New_Accession)] <- dat.all.denv1$Accession[is.na(dat.all.denv1$New_Accession)]
# dat.all.denv1$tip_label <- paste0(dat.all.denv1$New_Accession,"_",  dat.all.denv1$Locality, "_", dat.all.denv1$year)

# # #and overwrite
# # names(tree.denv1)[names(tree.denv1)=="tip_label"] <- "old_tip_label"
# # new.dat <- merge(tree.denv1, dat.all.denv1, by = "old_tip_label", all.x=T, sort=F)
# # 
# # rooted.D1$tip.label <- new.dat$tip_label
# 
# dat.all.denv1 <- dplyr::select(dat.all.denv1, tip_label, Accession, Locality, Serotype, Subclade)



#take a glance
pA <- ggtree(rooted.D1) %<+% dat.all.denv1 + geom_tippoint(aes(fill=Subclade), shape=21 ) +
  geom_tiplab(size=1, hjust=-.1) + #geom_nodelab(size=1, hjust=2) +
  scale_fill_manual(values=colz1, name=NULL) + theme(legend.position = c(.3,.83), 
                                          legend.background = element_rect(color="black")) +
  geom_cladelabel(node=gen2, label="Genotype-II", color="royalblue", offset = .0335) +
  geom_cladelabel(node=gen3, label="Genotype-III", color="darkgoldenrod1",  offset = .082, extend = c(0,5)) +
  geom_cladelabel(node=gen1, label="Genotype-I", color="mediumseagreen",  offset = .018) +
  xlim(c(0,.17)) +geom_treescale(fontsize=4, x=.005,y=110, linesize = .5) 

pA

A.dat <- pA$data

A.dat$node_fill <- NA
A.dat$node_fill[(length(dat.all.denv1$tip_label)+1):length(A.dat$label)] <- as.numeric(A.dat$label[(length(dat.all.denv1$tip_label)+1): length(A.dat$label)])#fill with label


pA2 <- pA  %<+% A.dat + 
        ggnewscale::new_scale_fill() + 
        geom_nodepoint(aes(fill=node_fill), shape=21, color="black", size=1, stroke=.1) + 
        scale_fill_continuous(low="yellow", high="red", limits=c(0,100), name="Bootstrap\nSupport") + theme(legend.box = "horizontal", legend.title = element_text(size=10)) +
        ggnewscale::new_scale_fill() + geom_tippoint(aes(fill=Subclade), shape=21 ) + scale_fill_manual(values=colz1, name=NULL)
#and do DENV2

dat.denv2 <- subset(dat, Serotype!="DENV-1")

tree.denv2 <- cbind.data.frame(tip_label = rooted.D2$tip.label)
tree.denv2$Accession <- NA
tree.denv2$Accession <- sapply(strsplit(tree.denv2$tip_label, split = "_Cambodia"), function(x) x[[1]])
tree.denv2$Accession[tree.denv2$Accession==tree.denv2$tip_label] <- NA
tree.denv2$Accession[is.na(tree.denv2$Accession)] <- sapply(strsplit(tree.denv2$tip_label[is.na(tree.denv2$Accession)], split = "_"), function(x) x[[1]])

dat.all.denv2 <- merge(tree.denv2, dat.denv2, by ="Accession", all.x = T, sort=F)
dat.all.denv2$Accession[dat.all.denv2$tip_label=="NC_002640_DENV4"] <- "NC_002640"
dat.all.denv2$Serotype[dat.all.denv2$tip_label=="NC_002640_DENV4"] <- "DENV-4"
dat.all.denv2$Subclade[dat.all.denv2$tip_label=="NC_002640_DENV4"] <- "Outgroup"

dat.all.denv2 <- dplyr::select(dat.all.denv2, tip_label, Accession, Locality, Serotype, Subclade, Collection_Year)


#take a glance
p2 <- ggtree(rooted.D2) %<+% dat.all.denv2 + geom_tippoint(aes(fill=Subclade), shape=21 ) +
  geom_tiplab(size=1) + geom_nodelab(size=1) + #scale_fill_manual(values=colz1) + 
  theme(legend.position = c(.2,.85), 
  legend.title = element_blank())

p2 #looks great - change to clade label and collapse the outgroup

unique(dat.all.denv2$Subclade)

#now shrink the long branch length for the outgroup
sort(rooted.D2$edge.length) #.07 is the next longest
rooted.D2$edge.length[rooted.D2$edge.length==max(rooted.D2$edge.length)] <- .14



dat.all.denv2$Subclade[dat.all.denv2$Subclade=="American"] <- "Genotype-I"
dat.all.denv2$Subclade[dat.all.denv2$Subclade=="Asian-American"] <- "Genotype-III"
dat.all.denv2$Subclade[dat.all.denv2$Subclade=="Asian-I"] <- "Genotype-V"
dat.all.denv2$Subclade[dat.all.denv2$Subclade=="Asian-II"] <- "Genotype-IV"


dat.all.denv2$Subclade[dat.all.denv2$Subclade=="Cosmopolitan-I-C-1B" |dat.all.denv2$Subclade=="Cosmopolitan-I-C-1A"] <- "Cosmopolitan-I"

dat.all.denv2$Subclade <- factor(dat.all.denv2$Subclade, 
                                 levels=c("Genotype-I", 
                                          "Cosmopolitan-I",
                                          "Cosmopolitan-II", 
                                          "Cosmopolitan-III",
                                          "Genotype-III", 
                                          "Genotype-IV", 
                                          "Genotype-V", 
                                          "Old-Cambodia",
                                          "New-Cambodia",
                                          "Outgroup")) 


colz2 = c('Genotype-I'= "cornflowerblue", 
          'Cosmopolitan-I' = "tomato",
          'Cosmopolitan-II' = "mediumseagreen", 
          'Cosmopolitan-III' = "royalblue",
          'Genotype-III' = "darkorange1", 
          'Genotype-IV' = "darkgoldenrod1", 
          'Genotype-V' = "forestgreen", 
          'Old-Cambodia'="palevioletred",
          'New-Cambodia'="mediumpurple4",
          'Outgroup' = "black")

cosmoIII <- MRCA(rooted.D2, which(rooted.D2$tip.label == "EU179857_Brunei_2005" ),which(rooted.D2$tip.label == "OL414753_Cambodia_2020"))
cosmoII <- MRCA(rooted.D2, which(rooted.D2$tip.label == "GQ398259_Indonesia_1976" ),which(rooted.D2$tip.label == "GQ398258_Indonesia_1975"))
cosmoI <- MRCA(rooted.D2, which(rooted.D2$tip.label == "JX475906_India_2009" ),which(rooted.D2$tip.label == "KF041234_Pakistan_2011"))
American <- MRCA(rooted.D2, which(rooted.D2$tip.label == "GQ868592_Columbia_1986" ),which(rooted.D2$tip.label == "HM582099_Fiji_1971"))
AsAm <- MRCA(rooted.D2, which(rooted.D2$tip.label == "HQ999999_Guatemala_2009" ),which(rooted.D2$tip.label == "FJ639700_Cambodia_2002"))
Asian1 <- MRCA(rooted.D2, which(rooted.D2$tip.label == "GQ868591_Thailand_1964" ),which(rooted.D2$tip.label == "OL414721_Cambodia_2019"))
Asian2 <- MRCA(rooted.D2, which(rooted.D2$tip.label == "HQ891024_Taiwan_2008" ),which(rooted.D2$tip.label == "JF730050_Puerto_Rico_2007"))

# #rooted.D2$node.label[as.numeric(rooted.D2$node.label)<80] <- ""
# 
# 
# #and make new tip labels
# names(dat.all.denv2)[names(dat.all.denv2)=="tip_label"] <- "old_tip_label"
# dat.all.denv2$year <- sapply(strsplit(dat.all.denv2$old_tip_label, "_"),  function(x) x[[length(x)]] )
# dat.all.denv2$New_Accession[is.na(dat.all.denv2$New_Accession)] <- dat.all.denv2$Accession[is.na(dat.all.denv2$New_Accession)]
# dat.all.denv2$tip_label <- paste0(dat.all.denv2$New_Accession,"_",  dat.all.denv2$Locality, "_", dat.all.denv2$year)
# 
# 
# #and overwrite
# names(tree.denv2)[names(tree.denv2)=="tip_label"] <- "old_tip_label"
# new.dat <- merge(tree.denv2, dat.all.denv2, by = "old_tip_label", all.x=T, sort=F)
# 
# rooted.D2$tip.label <- new.dat$tip_label
# 
# dat.all.denv2 <- dplyr::select(dat.all.denv2, tip_label, Accession, Locality, Serotype, Subclade)



#take a glance
pB <- ggtree(rooted.D2) %<+% dat.all.denv2 + geom_tippoint(aes(fill=Subclade), shape=21 ) +
  geom_tiplab(size=1, hjust=-.1) + #geom_nodelab(size=1, hjust=2) +
  scale_fill_manual(values=colz2, name=NULL) + theme(legend.position = c(.13,.76), 
                                          legend.background = element_rect(color="black"),
                                          legend.title = element_blank()) +
  geom_cladelabel(node=cosmoIII, label="Cosmopolitan-III", color="royalblue", offset = .038, extend = c(0,8), fontsize = 3.5) +
  geom_cladelabel(node=cosmoII, label="Cosmopolitan-II", color="mediumseagreen",  offset = .081, extend = c(1,1), fontsize = 3.5) +
  geom_cladelabel(node=cosmoI, label="Cosmopolitan-I", color="tomato",  offset = .0475, fontsize = 3.5) +
  geom_cladelabel(node=American, label="Genotype-I", color="cornflowerblue",  offset = .014, fontsize = 3.5) +
  geom_cladelabel(node=AsAm , label="Genotype-III", color="darkorange1",  offset = .058, fontsize = 3.5) +
  geom_cladelabel(node=Asian1 , label="Genotype-V", color="forestgreen",  offset = .056, fontsize = 3.5) +
  geom_cladelabel(node=Asian2 , label="Genotype-IV", color="darkgoldenrod1",  offset = .13, extend = c(1,1), fontsize = 3.5) +
  xlim(c(0,.16)) +geom_treescale( x=.002,y=110, linesize = .5, fontsize = 3.5) 

pB


B.dat <- pB$data

B.dat$node_fill <- NA
B.dat$node_fill[(length(dat.all.denv2$tip_label)+1):length(B.dat$label)] <- as.numeric(B.dat$label[(length(dat.all.denv2$tip_label)+1):length(B.dat$label)])#fill with label


pB2 <- pB  %<+% B.dat + 
  ggnewscale::new_scale_fill() + 
  geom_nodepoint(aes(fill=node_fill), shape=21, color="black", size=1, stroke=.1, show.legend = F) + 
  scale_fill_continuous(low="yellow", high="red", limits=c(0,100)) +
  ggnewscale::new_scale_fill() + 
  geom_tippoint(aes(fill=Subclade), shape=21 ) +
  scale_fill_manual(values=colz2, name=NULL) 
#and do DENV2

#and together

FigS3 <- cowplot::plot_grid(pA2, pB2, ncol = 2, nrow = 1, labels = c("A", "B"),label_size = 18)


ggsave(file = paste0(homewd, "/final-figures/FigS18.png"),
       plot=FigS3,
       units="mm",  
       width=100, 
       height=60, 
       scale=2.8)#, 



# now quantify by year!
cam.dat.denv1 = subset(dat.all.denv1, Locality == "Cambodia")
cam.dat.denv2 = subset(dat.all.denv2, Locality == "Cambodia")
cam.dat.denv1$year <- sapply(strsplit(cam.dat.denv1$tip_label, "Cambodia_"), '[', 2)
cam.dat.denv2$year <- sapply(strsplit(cam.dat.denv2$tip_label, "Cambodia_"), '[', 2)


cam.dat.denv1 = subset(cam.dat.denv1, year >=2019)
cam.dat.denv2 = subset(cam.dat.denv2, year >=2019)

tot <- rbind(cam.dat.denv1, cam.dat.denv2)
tot$Subclade <- as.character(tot$Subclade)
unique(tot$Subclade)
tot$Subclade[tot$Serotype=="DENV-1"] <- "Genotype-I"
tot$Subclade[tot$Serotype=="DENV-2"] <- "Cosmopolitan-III"
tot$Subclade[tot$Accession=="OL414721"] <- "Genotype-V"
tot$Subclade[tot$Accession=="OL414720"] <- "Genotype-V"
tot$Subclade[tot$Accession=="OL414744"] <- "Genotype-V"
tot$Subclade[tot$Accession=="OL414719"] <- "Genotype-V"
tot$Subclade[tot$Accession=="OL414743"] <- "Genotype-V"
tot$Subclade[tot$Accession=="OL414735"] <- "Genotype-V"
tot$Subclade[tot$Accession=="OL414734"] <- "Genotype-V"
tot$Subclade[tot$Accession=="109-1990S_L1"] <- "Genotype-V"
tot$Subclade[tot$Accession=="OQ678012"] <- "Genotype-V"
tot$Subclade[tot$Accession=="OL414761"] <- "Genotype-V"
tot$Subclade[tot$Accession=="OL414729"] <- "Genotype-V"
tot$Subclade[tot$Accession=="OL414728"] <- "Genotype-V"

# and summarize

# by year
dat.sum <- ddply(tot,.(year, Serotype, Subclade), summarise, N=length(Subclade))
# DENV-1 vs. DENV-2
# 2019 : 54 DENV-1 (61%) and 35 DENV-2 (39%)
# 2020 : 6 DENV-1 (20.7%) and 35 DENV-2 (79%)
# 2021 : 2 DENV-1 (25%) and 9 DENV-2 (75%)
# 2022 : 0 DENV-1 (0%) and 54 DENV-2 (100%)

# Cosmopolitan vs. Asian 1
# 2019 : 9 Genotype-V (25.7%) and 26 Cosmo (74.3%)
# 2020 : 2 Genotype-V (8.7%) and 21 Cosmo (91.3%)
# 2021 : 0 Genotype-V (0%) and 9 Cosmo (100%)
# 2022 : 1 Genotype-V (1.9%) and 53 Cosmo (98.1%)


#load beast data with ages
beast.dat <- read.csv(file = paste0(homewd, "data/beasttree_metadata.csv"))
head(beast.dat)
names(beast.dat)[names(beast.dat)=="accession_num"] <- "Accession"
beast.dat <- dplyr::select(beast.dat, Accession, lat, long, notes, age, sex)

tot.merge <- merge(tot, beast.dat, by="Accession", all.x=T)
head(tot.merge)
tail((tot.merge))


ggplot(data=tot.merge) + geom_boxplot(aes(x=Subclade, y=age, color=year))

# the one 2022 that is not DENV-2 cosmopolitan is from Phnom Penh
