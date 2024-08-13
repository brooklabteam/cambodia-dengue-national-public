


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
#library(rgeos)
library(sf)
library(lubridate)
library(ggspatial)


#make three timetrees and one map all together
#read in tree from beast
#first, make map

#and load the metadata
#homewd= "/home/rstudio"
homewd = "/Users/carabrook/Developer/cambodia-dengue-national-public"
setwd(homewd)


###################
###################

#first panel, map of cambodia with all the sequences

#and get submap of kampong speu witht the points of cambodia sequences
cam = sf::st_read(paste0(homewd, "/data/province-shape/khm_admbnda_adm1_gov_20181004.shp"))
#sub = subset(cam, name=="Kampong Speu")


dat <- read.csv(file = paste0(homewd, "/data/beasttree_metadata.csv"), header = T, stringsAsFactors = F)

dat$date <- as.Date(dat$date, format = "%m/%d/%y")
dat$year <- year(dat$date)
#dat = subset(dat, DENV.serotype=="DENV-2")
#names(dat)
#dat$date <- as.Date(dat$date, format = "%m/%d/%y")
#dat.plot = subset(dat, !is.na(lat) & DENV.serotype=="DENV-1" & year >2018) #49 sequences
dat.plot = dat
#make new labels
#dat.plot$date <- as.Date(dat.plot$date)

head(dat.plot)
dat.plot$date <- as.character(as.Date(dat.plot$date))
#dat.plot$date <- as.character(as.Date(dat.plot$date, format = "%m/%d/%y"))
dat.plot$new_label = paste(dat.plot$NIH.ID, dat.plot$date, sep = " ")
dat.plot$new_label


#dat.plot  <- dplyr::select(dat.plot, new_label, unique.ID, NIH.ID, CZB.ID, beast_name, date, year, lat, long, accession_num, DENV.serotype)
head(dat.plot)
get.centroid.bb <- function(x){
  N <- length(x)  # Number of polygons
  # Initialise data.frame
  Centroids.bb <- data.frame(matrix(NA, N, 2, dimnames = list(NULL, c("long", "lat"))))
  for(i in 1:N){
    # Bounding box of polygon
    bb <- bbox(x@polygons[[i]])
    # Compute centroid
    Centroids.bb[i,] <- c(
      0.5 * (bb[1,1] + bb[1,2]),
      0.5 * (bb[2,1] + bb[2,2]))
  }
  return(Centroids.bb)
}
jitter.dup <- function(pf,perc_jitter){
  if(nrow(pf)>1){
    nseq = nrow(pf)
    pf$lat[2:nrow(pf)] <- rnorm((nseq-1), mean=pf$lat[1], sd = perc_jitter)
    pf$long[2:nrow(pf)] <- rnorm((nseq-1), mean=pf$long[1], sd = perc_jitter)
  }
  return(pf)
}
get.spatial.object <- function(dat1, dist.thresh, denv.serotype, perc_jitter){
  
  #add in some local points and cluster them 
  xy <- SpatialPointsDataFrame(
    matrix(c(dat1$long, dat1$lat), ncol=2), data.frame(ID=seq(1:length(dat1$lat))),
    proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  
  # use the distm function to generate a geodesic distance matrix in meters
  mdist <- distm(xy)
  # cluster all points using a hierarchical clustering approach
  hc <- hclust(as.dist(mdist), method="complete")
  # define the distance threshold in meters. 
  #d=2000 #this is the typical
  d=dist.thresh #this is within 10km
  #d=500
  
  
  # define clusters based on a tree "height" cutoff "d" and add them to the SpDataFrame
  xy$clust <- cutree(hc, h=d)
  
  # expand the extent of plotting frame
  xy@bbox[] <- as.matrix(extend(extent(xy),0.001))
  
  
  # get the centroid coords for each cluster
  
  
  cent <- matrix(ncol=2, nrow=max(xy$clust))
  for (i in 1:max(xy$clust)){
    
    #convert to polygon
    ### Convert the SpatialPointsDataFrame to SpatialPolygons
    Sr1 = Polygon(subset(xy, clust==i))
    Srs1 = Polygons(list(Sr1), "s1")
    SpP = SpatialPolygons(list(Srs1), 1:1, proj4string= CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
    #plot(SpP, col = 3:3, pbg="white", add=T) 
    
    
    ### Convert the SpatialPolygons to SpatialPolygonsDataFrame
    shape_pol <- SpatialPolygonsDataFrame(SpP, match.ID=F, data= data.frame(x=xy[1:1,1], y=xy[1:1,2]))
    #shape_pol ### can be write as shapefile
    #plot(shape_pol, col = 4, add=T)
    
    cent[i,1:2] <- c(unlist(get.centroid.bb(shape_pol)))
  }
  
  
  # compute circles around the centroid coords using a radius of d meters
  # from the dismo package
  ci <- circles(cent, d=d, lonlat=T)
  
  dat1$cluster_ID <- NA
  dat1$cluster_ID[!is.na(dat1$lat)] <- xy$clust
  dat1$cluster_ID <- as.factor(dat1$cluster_ID)
  head(dat1)
  
  slim.gis <- dplyr::select(dat1, lat, long, cluster_ID, DENV.serotype)
  slim.gis  <- slim.gis[complete.cases(slim.gis),]
  #coord.dat <- st_as_sf(slim.gis, coords = c( "long", "lat"), crs = 4326)
  #circle.dat <- st_as_sf(ci@polygons, coords=c("X1", "X2"), crs=4326)
  
  #or the centroid dat
  if(denv.serotype=="all"){
    cent.gis <- ddply(slim.gis, .(cluster_ID), summarize, N=length(lat), DENV.serotype="all")
  }else if(denv.serotype!="all"){
    cent.gis <- ddply(slim.gis, .(DENV.serotype, cluster_ID), summarize, N=length(lat))
  }
  
  cent.merge <- cbind.data.frame(lat = cent[,2], long=cent[,1], cluster_ID = 1:nrow(cent))
  #cent.merge$cluster
  #cent.gis$long <- cent[,1]
  #cent.gis$lat <- cent[,2]
  cent.gis <- merge(cent.gis, cent.merge, by="cluster_ID", all_x=T)
  
  #for any values that are duplicates, go ahead and jitter
  cent.list <- dlply(cent.gis, .(lat))
  
  cent.gis <- data.table::rbindlist(lapply(cent.list, jitter.dup, perc_jitter=perc_jitter))
  
  circlpts.dat <- st_as_sf(cent.gis, coords = c( "long", "lat"), crs = 4326)
  
  
  
  
  return(list(circlpts.dat, dat1)) 
}


#first load the data
df1 <- read.csv(file = paste0(homewd,"/data/Cambodia-Serotype-Genotype.csv"), header = T, stringsAsFactors = F)
head(df1)
df1$date <- as.Date(df1$date, format = "%m/%d/%y")
df1$year <- year(df1$date)
names(df1)[names(df1)=="Genotype"] <- "DENV.serotype"
df1$DENV.serotype[df1$DENV.serotype=="Cosmopolitan"] <- "DENV-2-Cosmopolitan"
df1$DENV.serotype[df1$DENV.serotype=="Asian-1"] <- "DENV-2-Genotype-V"
#df1.hold = df1
df1 = subset(df1, !is.na(lat))
df1 = subset(df1, !is.na(long))


#total
ddply(df1,.(study, province), summarise, N_tot=length(sample_name))
ddply(df1,.( province), summarise, N_tot=length(sample_name))
#IDSeq or PAGODAs: 193+13+45 = 251 with 238 from KS
#NDCP = 21 with 7 from KS

#27 sequences not from KS
#21 sequences from NDCP 


#and make your spatial object
all.denv.partial <- get.spatial.object(dat1=df1, dist.thresh = 50000, denv.serotype = "merge", perc_jitter = 0.1)
denv.map.partial <- all.denv.partial[[1]]
dat.plot.cluster.partial <- all.denv.partial[[2]]
dat.plot.cluster.partial <- dplyr::select(dat.plot.cluster.partial,sample_name, cluster_ID)


#and plot
sero.cols = c("DENV-1"="forestgreen", "DENV-2-Genotype-V"="dodgerblue", "DENV-2-Cosmopolitan" = "navy", "DENV-4"="tomato")

pCamSummary <- ggplot(cam) + 
  geom_sf(fill="#9590FF", color="black", alpha=.8, size =.4) + #"#9590FF"
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        legend.box.background = element_blank(),
        #legend.position = c(.7,.05),
        legend.position = c(.6,.87),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        #legend.direction = "horizontal",
        legend.direction = "vertical",
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) +
  geom_sf(data=denv.map.partial, aes(fill=DENV.serotype, size=N),
          color="black",# "gray90",
          shape=21, stroke=.8) +
  coord_sf(xlim=c(102.3,107.63), ylim = c(10.26,15.8)) +
  scale_size_continuous(trans = "log10", name="number of\ngenotyped\ncases") +
  scale_fill_manual(values=sero.cols) + #+name="DENV serotype/\ngenotype", 
  guides(fill="none") +
  annotation_scale(pad_x = unit(8, "cm"), pad_y = unit(2.7, "cm"),
                   width_hint=.18, text_cex = .9)
#       fill=guide_legend(override.aes = list(color="black")))


#and get the time series of all cases
head(df1)



df.sum <- ddply(df1,.(year, DENV.serotype), summarise, N=length(sample_name))
df.sum.year <- ddply(df1,.(year), summarise, Ntot=length(sample_name))

df.sum <- merge(df.sum, df.sum.year, by="year", all.x=T)
head(df.sum)
df.sum$proportion <- df.sum$N/df.sum$Ntot
df.sum$DENV.serotype <- factor(df.sum$DENV.serotype, levels=c("DENV-1", "DENV-4",
                                                              "DENV-2-Genotype-V",
                                                              "DENV-2-Cosmopolitan"))


sero.cols2 = c("DENV-1"="forestgreen", "DENV-2-Genotype-V"="dodgerblue", "DENV-2-Cosmopolitan" = "navy", "DENV-4"="tomato")

pBar <- ggplot(data=df.sum) + scale_fill_manual(values=sero.cols2) +
  geom_bar(aes(x=year, y=N, fill=DENV.serotype), 
           position = "stack", stat = "identity", show.legend = F) +
  theme_bw() + ylab("number of genotyped cases") + 
  #coord_cartesian(ylim=c(0,150)) +
  scale_y_continuous(breaks= c(30,60,90,120)) +
  theme(#panel.background = element_blank(),
    panel.grid = element_blank(),
    #legend.position = c(.27,.87),
    #legend.text = element_text(size=8),
    #legend.direction = "horizontal",
    #legend.title = element_text(size=10),
    #legend.box.background = element_rect(color="black"),
    axis.title.x = element_blank(), 
    plot.margin =  unit(c(.1,.3,.7,.1), units="cm"),
    axis.title.y = element_text(size=16),
    axis.text = element_text(size=14)) #+
#guides(fill=guide_legend(ncol=2))

#and ages
head(df1)
#ggplot(data=df1) + geom_jitter(aes(x=DENV.serotype, y=age), width=.1, size=.1, alpha=.3) + 
 # geom_violin(aes(x=DENV.serotype, y=age,  color=DENV.serotype), draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) #+ 
#facet_grid(~year)


m1 <- lm(log10(age)~DENV.serotype, data = df1)
summary(m1)

# histogram
hist(m1$residuals)#pretty normal


# QQ-plot
library(car)
qqPlot(m1$residuals,
       id = FALSE # id = FALSE to remove point identification
)

library(sjPlot)
plot_model(m1, type="est") #denv cosmopolitand and DENV-4 are both older in age than DENV-1
summary(m1)

df1$DENV.serotype[df1$DENV.serotype=="DENV-2-Genotype-V"] <- "DENV-2\n(Geno-V)"
df1$DENV.serotype[df1$DENV.serotype=="DENV-2-Cosmopolitan"] <- "DENV-2\n(Cosmo)"
df1$DENV.serotype <- factor(df1$DENV.serotype, levels=c("DENV-1", "DENV-2\n(Geno-V)","DENV-2\n(Cosmo)", "DENV-4" ))

label.dat = cbind.data.frame(DENV.serotype=c("DENV-1", "DENV-2\n(Geno-V)","DENV-2\n(Cosmo)", "DENV-4" ), text = c("", "", "***", "***"))
label.dat$DENV.serotype <- factor(label.dat$DENV.serotype, levels=c("DENV-1", "DENV-2\n(Geno-V)","DENV-2\n(Cosmo)", "DENV-4" ))

sero.cols = c("DENV-1"="forestgreen", "DENV-2\n(Geno-V)"="dodgerblue", "DENV-2\n(Cosmo)" = "navy", "DENV-4"="tomato")

##nicer plot
pViolin <- ggplot(data=df1) + ylab("age of genotyped cases") + 
  facet_grid(~DENV.serotype, scales = "free_x") +
  geom_jitter(aes(x=DENV.serotype, y=age), width=.1, size=.1, alpha=.3) + theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x  = element_blank(), 
        strip.background = element_rect(fill="white", color="white"),
        axis.ticks = element_blank(),
        axis.text.y  = element_text(size=8), 
        plot.margin =  unit(c(0,0,0,0), units="cm"),
        panel.spacing = unit(0, units="cm"),
        strip.text = element_text(size=8),
        axis.title = element_text(size=8), 
        axis.title.x = element_blank())+ coord_cartesian(ylim=c(0,55)) +
  geom_violin(aes(x=DENV.serotype, y=age,  color=DENV.serotype), draw_quantiles = c(0,.25,.5,.75), show.legend=F, fill=NA) +
  geom_label(data = label.dat, aes(x=DENV.serotype, y= 53, label=text), label.size = 0, size=5) +
  scale_y_continuous(breaks = seq(0,50, by=10)) + 
  scale_x_discrete(position="top") +
  scale_color_manual(values = sero.cols)



#embed violin in map 102.3338 ymin: 10.26162 xmax: 107.6277 ymax: 14.69025

pA <-  pCamSummary + annotation_custom(ggplotGrob(pViolin), xmin = 101.5, xmax = 105.1,  ymin = 13.8, ymax = 16) #add map as embedded subplot

#pB <-  pBar + annotation_custom(ggplotGrob(pViolin), xmin = 2018, xmax = 2021.5,  ymin = 97, ymax = 150) #add map as embedded subplot

# #remake
# pBC <- cowplot::plot_grid(pViolin, pBar, ncol=1, nrow =2,  labels = c("B", "C"), label_size = 22, label_x = c(-.03,-.03) )


#and put together with map
pAB <- cowplot::plot_grid(pA, pBar, ncol = 1, nrow = 2, labels=c("A", "B"), label_size = 22, rel_heights = c(1,1), label_x = c(-.01,-.01))
#pAB <- cowplot::plot_grid(pCamSummary, pB, ncol = 1, nrow = 2, labels=c("A", "B"), label_size = 22, rel_heights = c(1,1))
#pS3_KPS<- cowplot::plot_grid(pKPS,  nrow=1, ncol=1, labels = c("A"), label_size = 22)


#first, the colors
colorz = scales::hue_pal()(10)
#and the names
names(colorz) <- c("Singapore", "Brunei", "Vietnam","Malaysia", "Philippines", "Thailand", "Myanmar", "Cambodia", "Laos", "Indonesia")

#world map
world <- ne_countries(scale = "medium", returnclass = "sf")
#subset
SEA = subset(world, sovereignt=="Singapore"| sovereignt=="Brunei"| sovereignt=="Vietnam" | sovereignt==  "Cambodia"| sovereignt== "Philippines"| sovereignt== "Thailand"| sovereignt== "Myanmar"| sovereignt== "Malaysia"| sovereignt== "Laos"| sovereignt== "Indonesia")


#plot in color
# gene world map

pSEA <- ggplot(data = SEA) +
  geom_sf(aes(fill=sovereignt), show.legend = F, color="black", size=.2) +
  #labs( x = "Longitude", y = "Latitude") +
  scale_fill_manual(values = colorz) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines"),
        #plot.margin = unit(c(10,3,10,3),"lines"), 
        #panel.background  = element_rect(size=1, fill = NULL, color = "black"),
        panel.background  = element_rect(fill="white"),
        panel.border  = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank()) 
coord_sf(xlim=c(95,170), expand = T)
psub<-pSEA+
  theme(plot.margin = unit(c(0,0,0,-1.5), "cm"))

#now make into a grob, so it can get embedded as a subplot in the phylogeny



#now, make th submap



#now add the trees
tree1 <- read.beast(file = paste0(homewd, "/BEAST-tree/denv-beast/beast-out/denv1/denv1Avg.tree"))
tree2 <- read.beast(file = paste0(homewd, "/BEAST-tree/denv-beast/beast-out/denv2/denv2Avg.tree"))


#replace these (did not have accession numbers yet when first made tree)
tree1@phylo$tip.label[tree1@phylo$tip.label=="DENS-OB-067_2021-05-29"] <- "PP470671_2021-05-29" 
tree1@phylo$tip.label[tree1@phylo$tip.label=="DENS-OB-068_2021-05-29"] <- "PP470672_2021-05-29" 
tree1@phylo$tip.label[tree1@phylo$tip.label=="DENS-OB-070_2021-05-30"] <- "PP470673_2021-05-30" 


tree1dat <- cbind.data.frame(tip_name = tree1@phylo$tip.label)

tree2dat <- cbind.data.frame(tip_name = tree2@phylo$tip.label)


head(tree1dat)
head(tree2dat)

#and load the metadata

 

dat <- read.csv(file = paste0(homewd, "/data/beasttree_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)


setdiff(tree1dat$tip_name, dat$tip_name) #these get edited down below
setdiff(tree2dat$tip_name, dat$tip_name)



#check the format
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
#dat$date <- as.Date(dat$date, format = "%m/%d/%y")



mrsd.denv1 <- max(dat$date[dat$DENV.serotype=="DENV-1"]) #"2021-05-29"
mrsd.denv2 <- max(dat$date[dat$DENV.serotype=="DENV-2"])#"2022-10-10"

node.tree1 <- MRCA(tree1, which(tree1@phylo$tip.label== "OL412678_2019-07-25" ),which(tree1@phylo$tip.label == "ON046271_2004-11-20"))


pA1 <- ggtree(tree1, mrsd=mrsd.denv1, color="forestgreen")  + 
  geom_cladelab(node=node.tree1, label="Genotype-I", textcolor="seagreen", barcolor="seagreen", fontsize=6,
                offset =-37, angle=270, offset.text = -12, vjust=2, hjust=.5)  +
  theme_tree2() + coord_cartesian(xlim=c(1930,2030), ylim=c(0,400)) + 
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1, show.legend = F) +
  scale_fill_continuous(low="yellow", high="red", limits=c(0,1)) +
  scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  #scale_fill_continuous(low="yellow", high="red", limits=c(0,1))+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) +
  annotation_custom(ggplotGrob(psub), xmin = 1935, xmax = 1980,  ymin = 140, ymax = 385) #add map as embedded subplot



#node.tree2.1 <- MRCA(tree2, which(tree2@phylo$tip.label== "OL414741_2019-07-15" ),which(tree2@phylo$tip.label == "KU509277_2010-07-31"))
node.tree2 <- MRCA(tree2, which(tree2@phylo$tip.label== "PP411229_2022-06-01" ),which(tree2@phylo$tip.label == "MK506264_2007-06-15"))
#node.tree2.2 <- MRCA(tree2, which(tree2@phylo$tip.label== "OL414721_2019-07-15"),which(tree2@phylo$tip.label == "KF744400_2000-07-31"))
node.tree2.1 <- MRCA(tree2, which(tree2@phylo$tip.label== "OQ426897_2019-03-06" ),which(tree2@phylo$tip.label == "KF921930_2002-07-31"))

pB2 <- ggtree(tree2, mrsd=mrsd.denv2, color="navy")  + theme_tree2() + 
  coord_cartesian(xlim=c(1930,2030),  ylim=c(0,400))+
  geom_cladelab(node=node.tree2, label="Cosmopolitan-III", textcolor="tomato",barcolor="tomato",
              offset =-37, angle=270, offset.text = -12, fontsize=6, vjust=2, hjust=.5)  +
  geom_cladelab(node=node.tree2.1, label="Genotype-V", textcolor="navy", barcolor="navy", fontsize=6,vjust=2, hjust=.5,
               offset =-37, angle=270, offset.text = -12)  +
  #geom_tiplab(size=1)+
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1) +
  scale_fill_continuous(low="yellow", high="red", limits=c(0,1)) +
  scale_x_continuous(breaks=c(1950, 1975, 2000, 2020))+
  theme(legend.position = c(.15,.85), 
        legend.key.size = unit(.4, units="cm"),
        legend.text = element_text(size=9),
        legend.title = element_text(size=10))



tree1merge <- merge(x=tree1dat, y=dat, by="tip_name", all.x = T, sort = F)
tree2merge <- merge(x=tree2dat, y=dat, by="tip_name", all.x = T, sort = F)
head(tree1merge)
tail(tree1merge)
head(tree2merge)



#and attach clusterID - for suppplot
#dat.clust.save <- dplyr::select(dat.clust.save, -(new_label))
#tree1merge <- merge(x=tree1merge, y=dat.clust.save, by="beast_name", all.x = T, sort = F)
#tree2merge <- merge(x=tree2merge, y=dat.clust.save, by="beast_name", all.x = T, sort = F)

tree1merge$new_label = sapply(strsplit(tree1merge$tip_name, "_"), function(x) x[[1]])
tree1merge$new_label <- paste0(tree1merge$new_label, " ", as.character(tree1merge$date))

tree1merge$new_seq = "no"
tree1merge$new_seq[tree1merge$country=="Cambodia" & !is.na(tree1merge$NIH.ID)] <- "yes"
tree1merge$new_seq <- as.factor(tree1merge$new_seq)

tree1merge$CambodiaSeq <- "no"
tree1merge$CambodiaSeq[tree1merge$country=="Cambodia"] <- "yes"



#unique(tree2merge$country)
#tree2merge$new_label = sapply(strsplit(tree2merge$tip_name, "_"), function(x) x[[1]])
#tree2merge$new_label <- paste0(tree2merge$new_label, " ", as.character(tree2merge$date))

tree2merge$new_seq = "no"
tree2merge$new_seq[tree2merge$country=="Cambodia" &  !is.na(tree2merge$NIH.ID)] <- "yes"
tree2merge$new_seq <- as.factor(tree2merge$new_seq)

tree2merge$CambodiaSeq <- "no"
tree2merge$CambodiaSeq[tree2merge$country=="Cambodia"] <- "yes"

shapez = c("yes"=24, "no"=21)


pA2 <- pA1 %<+% tree1merge +
  ggnewscale::new_scale_fill() +
  #geom_label(aes(label=new_label), label.size = NA, size=4, hjust=-.06) +
  #geom_tippoint(aes(shape=CambodiaSeq), color="black", fill="black", size = 4)+#,show.legend = F) +
  geom_tippoint(aes(fill=country, shape=new_seq), 
                size = 2, show.legend = F, stroke=.1, color="black") +
  theme(legend.position = c(.08,.8), 
        legend.key.size = unit(.4, units="cm"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=11),
        #plot.background  = element_rect(size=3, fill = NULL, color = "black"),
        axis.text = element_text(size=18)) +
  scale_fill_manual(values=colorz) +
  scale_shape_manual(values=shapez) 




tree1.tiplabel<-tree1@phylo[["tip.label"]]




######################################################################
#################### add bar to subplot B tree #######################
######################################################################



pB <-pB2 %<+% tree2merge + 
  ggnewscale::new_scale_fill()  +
  geom_tippoint(aes(fill=country, shape=new_seq), size = 2, show.legend = F, stroke=.1, color="black") +
  theme(axis.text = element_text(size=18)) + 
  #plot.background  = element_rect(size=3, fill = NULL, color = "black"),
  #legend.position = "none",
  #legend.title = element_blank()) + 
  #coord_cartesian(c(2000,2021)) +
  scale_shape_manual(values=shapez) +scale_fill_manual(values=colorz)


pCD <-cowplot::plot_grid(pA2,pB,nrow=2,ncol=1,labels=c("C", "D"),label_size=22)+
  theme(plot.margin = unit(c(0,0,.5,0), "cm"))




######################################################################
#### now, combine with transmission trees plot (D) ###################
######################################################################

#load the transmission tree data
#here for both serotypes
all.denv <- read.csv(file=paste0(homewd, "/data/AllDENVtransTreeDat.csv"), header = T, stringsAsFactors = F)
all.denv$distance <- all.denv$distance/1000 #convert to km

#ggplot(dat) + geom_point(aes(x=evol_time,y=distance)) + facet_grid(~DENV.subtype) 

#are they in the same season?
all.denv$season <- "yes"
all.denv$season[all.denv$pairtime2-all.denv$pairtime1>.5] <- "no"

#now, only look at those within a season
all.denv <- subset(all.denv, season=="yes")
unique(all.denv$paired)


#dat=subset(all.denv,DENV.serotype=="DENV-1")
mrca_thresh=.5
#and write over to get the transmission trees
#geothresh <- list(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) 
geothresh <- list(.2,.4,.6,.8,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) 
#here's the threshold list of distances over which to compute

#here, decides whether each pair is in a transmission chain, based on the mrca_thresh (in years)
is.trans.chain <- function(dat, mrca_thresh){
  dat$trans_chain <- 0
  dat$trans_chain[dat$tMRCA <= mrca_thresh] <- 1
  return(dat)
  
}
get.prop.chain.discrete<- function(thresh, dat, thresh.bin){
  if(which(thresh.bin==thresh)==1){
    min.dist = 0
  }else{
    min.dist = thresh.bin[[(which(thresh.bin==thresh)-1)]]  
  }
  tot.chain <- sum(subset(dat, distance<=thresh & distance>min.dist)$trans_chain)
  N.pairs <- length(subset(dat, distance<=thresh& distance>min.dist)$trans_chain)
  if(N.pairs==0){
    prop.chain <-NA
    dat.new <- cbind.data.frame(distance = thresh, prop=NA, prop_lci = NA, prop_uci= NA, tot_pairs_in_chain=NA, tot_pairs=NA)
  }else{
    prop.chain <- tot.chain/N.pairs
    CI <- binom.test(x=tot.chain, n=N.pairs, alternative = "two.sided", conf.level = .95)
    dat.new <- cbind.data.frame(distance = thresh, prop=prop.chain, prop_lci = CI$conf.int[1], prop_uci= CI$conf.int[2], tot_pairs_in_chain=tot.chain, tot_pairs=N.pairs)
    
  }
  
  
  
  
  
  #and get CIs
  
  
  return(dat.new)
}
get.prop.chain.max <- function(thresh, dat){
  
  tot.chain <- sum(subset(dat, distance<=thresh)$trans_chain)
  N.pairs <- length(subset(dat, distance<=thresh)$trans_chain)
  prop.chain <- tot.chain/N.pairs
  
  CI <- binom.test(x=tot.chain, n=N.pairs, alternative = "two.sided", conf.level = .95)
  
  
  #and get CIs
  
  dat.new <- cbind.data.frame(distance = thresh, prop=prop.chain, prop_lci = CI$conf.int[1], prop_uci= CI$conf.int[2], tot_pairs_in_chain=tot.chain, tot_pairs=N.pairs)
  return(dat.new)
}
make.trans.chains <- function(dat, geothresh, mrca_thresh, character){
  
  dat1 = is.trans.chain(dat=dat, mrca_thresh = mrca_thresh)
  
  #and proportion of transmission chains
  dat.prop = lapply(geothresh, get.prop.chain.max, dat=dat1)
  #dat.prop = lapply(geothresh, get.prop.chain.discrete, dat=dat1, thresh.bin=geothresh)
  prop.df = data.table::rbindlist( dat.prop)
  prop.df$character = character
  prop.df$DENV.serotype=unique(dat$DENV.serotype)
  prop.df$transchain_threshold = mrca_thresh
  prop.df$character <- character
  return(prop.df)
  
}

combine.chain.prop<-function(mrca.thresh){
  out.prop1 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-1"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  #out.prop3 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2" & paired=="DENV-2-Cosmopolitan/DENV-2-Cosmopolitan"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  #out.prop5 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2" & paired=="DENV-2-Asian-1/DENV-2-Asian-1"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  out.prop7 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  #head(out.prop1)
  out.prop1$DENV.subtype <- "DENV-1"
  #out.prop3$DENV.subtype <- "DENV-2-Cosmopolitan"
  out.prop7$DENV.subtype <- "DENV-2-All"
  #out.prop5$DENV.subtype <- "DENV-2-Asian-1"
  #out.all <-rbind(out.prop1, out.prop3, out.prop5, out.prop7)
  out.all <-rbind(out.prop1, out.prop7)
  return(out.all)
  
}
make.chain.diff <- function(mrca.thresh){
  out.prop1 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-1"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  out.prop3 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2" & paired=="DENV-2-Cosmopolitan/DENV-2-Cosmopolitan"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  out.prop5 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2" & paired=="DENV-2-Asian-1/DENV-2-Asian-1"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  out.prop7 <- make.trans.chains(dat=subset(all.denv,DENV.serotype=="DENV-2"), geothresh = geothresh, mrca_thresh = mrca.thresh, character = "whole-dataset")
  head(out.prop1)
  
  out.prop <- out.prop3 #Cosmopolitan is dominant and we subtract from it
  out.prop$prop_diff <- out.prop3$prop - out.prop1$prop
  out.prop$prop_diff_lci <-  out.prop$prop_lci-out.prop1$prop
  out.prop$prop_diff_uci <-  out.prop$prop_uci-out.prop1$prop
  out.prop$comp <- "DENV-1"
  
  
  out.prop2 <- out.prop3 #Cosmopolitan is dominant and we subtract from it
  out.prop2$prop_diff <- out.prop2$prop - out.prop7$prop
  out.prop2$prop_diff_lci <-  out.prop2$prop_lci-out.prop7$prop
  out.prop2$prop_diff_uci <-  out.prop2$prop_uci-out.prop7$prop
  #out.prop2$comp <- "DENV-2-Asian-1"#
  out.prop2$comp <- "All-DENV-2"
  
  out.prop <- rbind(out.prop, out.prop2)
  
  out.prop <- dplyr::select(out.prop, distance, prop_diff, prop_diff_lci, prop_diff_uci, transchain_threshold, comp)
  return (out.prop)
}
mean.trans.chains <- function(dat, geothresh, mrca_thresh){
  
  dat1 = is.trans.chain(dat=dat, mrca_thresh = mrca_thresh)
  
  #get the probability that two cases are from the same chain
  prob <- sum(dat1$trans_chain)/length(dat1$trans_chain)
  prob.uci <- prop.test(x=sum(dat1$trans_chain), n=length(dat1$trans_chain), alternative = "t", conf.level = .95)$conf.int[1]
  prob.lci <- prop.test(x=sum(dat1$trans_chain), n=length(dat1$trans_chain), alternative = "t", conf.level = .95)$conf.int[2]
  
  eff.chains <- 1/prob #22.264 chains circulating in the region (Denv1)
  eff.chains.lci <- 1/prob.lci #17.019 chains circulating in the region (Denv1)
  eff.chains.uci <- 1/prob.uci #29.298 chains circulating in the region (Denv1)
  
  #and return
  chain.df <- cbind.data.frame(DENV.serotype=unique(dat1$DENV.serotype), DENV.subtype=unique(dat1$DENV.subtype), N_chains=eff.chains, N_chains_lci=eff.chains.lci, N_chains_uci = eff.chains.uci)
  
  return(chain.df)
  #pop size KP: 877,523... this is a little higher than the prediction for rural Thailand and lower
  #than the prediction for Bangkok
  #and get the mean number of transmission chains for this region that are in circulation
  
  
}


#out.pt5 = make.chain.diff(.5)
out.pt5 = combine.chain.prop(mrca.thresh=.5)
head(out.pt5)
tail(out.pt5)
#out.prop = rbind(out.pt3, out.pt5, out.1)#, out.3)
#head(out.prop)
#out.pt5$DENV.subtype<- factor(out.pt5$DENV.subtype, levels = c("DENV-1", "DENV-2-All", "DENV-2-Cosmopolitan"))
#out.pt5$DENV.subtype<- factor(out.pt5$DENV.subtype, levels = c("DENV-1", "DENV-2-All", "DENV-2-Asian-1", "DENV-2-Cosmopolitan"))
out.pt5$DENV.serotype<- factor(out.pt5$DENV.serotype, levels = c("DENV-1", "DENV-2"))


#colzB=c("DENV-1"="mediumseagreen", "DENV-2-All"="navy", "DENV-2-Cosmopolitan"="dodgerblue")
#colzB=c("DENV-1"="mediumseagreen", "DENV-2-All"="navy", "DENV-2-Asian-1"="cyan",  "DENV-2-Cosmopolitan"="dodgerblue")
colzB=c("DENV-1"="mediumseagreen", "DENV-2"="navy")


pD <- ggplot(data=out.pt5) + theme_bw()+
  #facet_grid(dummy_label~.) +
  #geom_line(aes(x=distance, y=prop, color=sex),show.legend = F) +
  geom_ribbon(aes(x=distance, ymin=prop_lci, ymax=prop_uci, fill=DENV.serotype, group=DENV.serotype), alpha=.3) +
  geom_line(aes(x=distance, y=prop, color=DENV.serotype, group=DENV.serotype)) +
  theme(panel.grid = element_blank(), axis.title = element_text(size=16), 
        #axis.title.x = element_blank(), axis.text.x = element_blank(), 
        #axis.ticks.x = element_blank(),
        plot.margin = unit(c(.2,.3,.1,.5), "cm"),
        legend.title = element_blank(),
        strip.background = element_rect(fill="white"), strip.text = element_text(size = 18),
        axis.text = element_text(size=14), legend.text = element_text(size=12),
        legend.position = c(.73,.84)) + coord_cartesian(ylim=c(0,1), xlim=c(0,5), expand = F)+
  scale_color_manual(values=colzB) + 
  scale_fill_manual(values=colzB) + 
  #scale_color_manual(values=colz, name = "transmission chain\nthreshold (yrs)") + 
  #scale_fill_manual(values=colz, name = "transmission chain\nthreshold (yrs)") + 
  ylab("Proportion same transmission chain")  +
  xlab("Max distance between cases (km)")


######################################################################
################## and combine with data from Salje (C) ##############
######################################################################

#and plot with the salje data
salje.dat <- read.csv(file=paste0(homewd, "/data/salje_chains.csv"), header = T, stringsAsFactors = F)
head(salje.dat)

#first, take just the KS data
head(all.denv)
dat2 <- read.csv(file=paste0(homewd, "/data/Cambodia-Serotype-Genotype.csv"), header = T, stringsAsFactors = F)
head(dat2)
dat2 = subset(dat2, sequence_completeness=="complete" & province=="Kampong Speu")
dat2 <- dplyr::select(dat2, genbank_accession)
dat2$flag <- 1
names(dat2)[names(dat2)=="genbank_accession"] <- "accession_early"
all.denv <- merge(all.denv, dat2, by="accession_early", all.x=T)
all.denv = subset(all.denv, flag==1)
all.denv <- dplyr::select(all.denv, -(flag))

#and get mean chains
denv.1 = subset(all.denv, DENV.serotype=="DENV-1")
denv.1$DENV.subtype <- "DENV-1"
denv.1.mean = mean.trans.chains(dat=denv.1, geothresh, mrca_thresh)

denv.2 = subset(all.denv, DENV.serotype=="DENV-2")
denv.2$DENV.subtype <- "DENV-2"
denv.2.mean = mean.trans.chains(dat=denv.2, geothresh, mrca_thresh)


# denv.cosmo.2 = subset(all.denv, DENV.serotype=="DENV-2" & DENV.subtype=="DENV-2-Cosmopolitan")
# denv.cosmo.2$DENV.subtype <- "DENV-2-Cosmopolitan"
# denv.cosmo.2.mean = mean.trans.chains(dat=denv.cosmo.2, geothresh, mrca_thresh)

#all.denv.mean <- rbind(denv.1.mean, denv.2.mean, denv.cosmo.2.mean)
all.denv.mean <- rbind(denv.1.mean, denv.2.mean)
#fewer circulating chains for denv-2 vs. 1 and even fewer for co


#salje.dat$locale[salje.dat$locale=="Bangkok"] <- "Salje et al. Bangkok"
salje.dat$locale[salje.dat$locale=="Thai_Countryside"] <- "Rural Thailand"# "Salje et al.\nThai rural"
all.denv.mean$DENV.subtype[all.denv.mean$DENV.subtype=="DENV-1"] <- "Cambodia DENV-1"
all.denv.mean$DENV.subtype[all.denv.mean$DENV.subtype=="DENV-2"] <- "Cambodia DENV-2"
#all.denv.mean$DENV.subtype[all.denv.mean$DENV.subtype=="DENV-2-Cosmopolitan"] <- "Cambodia\nDENV-2-Cosmopolitan"

salje.dat$study <- "Salje et al. 2017"
all.denv.mean$study <- "Kampong Speu 2019-2022"

#colznew <- c('Bangkok' = "black", 'Rural Thailand' = "gray60", 'Cambodia DENV-1' = "forestgreen", 'Cambodia All-DENV-2' = "navy")
shapeznew <- c('Salje et al. 2017' = 21, "Kampong Speu 2019-2022" = 24)
colznew <- c('Bangkok' = "black", 'Rural Thailand' = "gray60", 'Cambodia DENV-1' = "forestgreen", 'Cambodia DENV-2' = "navy", 'Salje et al. 2017' = "black", "Kampong Speu 2019-2022" = "red")

#and plot
pC <- ggplot(data=salje.dat) + 
  geom_point(aes(x=pop_size, y=eff_chains, fill=locale, 
                 shape=study, color=locale), size=4, color="black", show.legend = F) +
  geom_errorbar(aes(x=pop_size, ymin=lci, ymax=uci, color=locale)) +
  scale_y_log10() + 
  scale_x_log10(breaks=c(1e+03, 1e+04, 1e+05, 1e+06, 1e+07), 
                labels=c(1,10,100,1000, 10000)) + theme_bw() +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.title = element_text(size=18),
        axis.text = element_text(size=14),
        plot.margin = unit(c(.3,.2,.1,.5), "cm"),
        #legend.background = element_rect(color="black"),
        legend.position = c(.35,.88)) +
  ylab("Effective # Chains") +
  xlab("Population Size (x1000)") +
  geom_point(data=all.denv.mean,aes(x=877523., y=N_chains, fill=DENV.subtype, shape=study, color=study), size=5,  stroke=2) +
  scale_color_manual(values=colznew, guide=NULL) + scale_fill_manual(values=colznew, guide=NULL) +
  scale_shape_manual(values = shapeznew) 


#put C and D together


pEF <-cowplot::plot_grid(pD,pC,nrow=2,ncol=1,labels=c("E", "F"),label_size=22)+
  theme(plot.margin = unit(c(0,0,0,0), "cm"))


#and all together

Fig4 <- cowplot::plot_grid(pAB, pCD, pEF, nrow=1, ncol = 3) + theme(plot.background = element_rect(fill ="white"))+ 
  theme_classic()+  theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),
                          axis.ticks=element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank())

ggsave(file = paste0(homewd, "/final-figures/Fig4.png"),
       plot= Fig4,
       units="mm",  
       width=140, 
       height=85, 
       scale=3, 
       dpi=300)

ggsave(file = paste0(homewd, "/final-figures/Fig4.pdf"),
       plot= Fig4,
       units="mm",  
       width=140, 
       height=85, 
       scale=3, 
       dpi=300)

############################
############################
############################

#Fig S20

# replot at several thresholds and store as a supplementary figure



#out.pt5 = make.chain.diff(.5)
out.1 = combine.chain.prop(mrca.thresh=1)
out.pt08 = combine.chain.prop(mrca.thresh=(1/12))
out.pt3 = combine.chain.prop(mrca.thresh=.25)
out.pt9 = combine.chain.prop(mrca.thresh=.75)
head(out.1)
head(out.pt3)
head(out.pt9)

out.prop = rbind(out.pt08, out.pt3, out.pt5, out.1)#, out.3)
#head(out.prop)
#out.pt5$DENV.subtype<- factor(out.pt5$DENV.subtype, levels = c("DENV-1", "DENV-2-All", "DENV-2-Cosmopolitan"))
#out.pt5$DENV.subtype<- factor(out.pt5$DENV.subtype, levels = c("DENV-1", "DENV-2-All", "DENV-2-Asian-1", "DENV-2-Cosmopolitan"))
out.prop$DENV.serotype<- factor(out.prop$DENV.serotype, levels = c("DENV-1", "DENV-2"))
out.prop$differentiate <- NA
out.prop$differentiate[out.prop$DENV.serotype=="DENV-1" & out.prop$transchain_threshold==(1/12)] <- "DENV-1: 1 mon"
out.prop$differentiate[out.prop$DENV.serotype=="DENV-1" & out.prop$transchain_threshold==0.25] <- "DENV-1: 3 mon"
out.prop$differentiate[out.prop$DENV.serotype=="DENV-1" & out.prop$transchain_threshold==0.5] <- "DENV-1: 6 mon"
#out.prop$differentiate[out.prop$DENV.serotype=="DENV-1" & out.prop$transchain_threshold==0.75] <- "DENV-1: 9 months"
out.prop$differentiate[out.prop$DENV.serotype=="DENV-1" & out.prop$transchain_threshold==1] <- "DENV-1: 12 mon"

out.prop$differentiate[out.prop$DENV.serotype=="DENV-2" & out.prop$transchain_threshold==(1/12)] <- "DENV-2: 1 mon"
out.prop$differentiate[out.prop$DENV.serotype=="DENV-2" & out.prop$transchain_threshold==0.25] <- "DENV-2: 3 mon"
out.prop$differentiate[out.prop$DENV.serotype=="DENV-2" & out.prop$transchain_threshold==0.5] <- "DENV-2: 6 mon"
#out.prop$differentiate[out.prop$DENV.serotype=="DENV-2" & out.prop$transchain_threshold==0.75] <- "DENV-2: 9 months"
out.prop$differentiate[out.prop$DENV.serotype=="DENV-2" & out.prop$transchain_threshold==1] <- "DENV-2: 12 mon"

out.prop$differentiate <- factor(out.prop$differentiate, levels = c("DENV-1: 1 mon", "DENV-1: 3 mon", "DENV-1: 6 mon", "DENV-1: 12 mon",
                                                                    "DENV-2: 1 mon", "DENV-2: 3 mon", "DENV-2: 6 mon", "DENV-2: 12 mon"))


#colzB=c("DENV-1"="mediumseagreen", "DENV-2-All"="navy", "DENV-2-Cosmopolitan"="dodgerblue")
#colzB=c("DENV-1"="mediumseagreen", "DENV-2-All"="navy", "DENV-2-Asian-1"="cyan",  "DENV-2-Cosmopolitan"="dodgerblue")
#colzB=c("DENV-1"="mediumseagreen", "DENV-2"="navy")
# colzC=c("DENV-1: 3 mon. MRCA" = "lightgreen" , "DENV-1: 6 mon. MRCA" = "green2", "DENV-1: 9 mon. MRCA"="green3", "DENV-1: 12 mon. MRCA"="forestgreen",
#         "DENV-2: 3 mon. MRCA"="lightblue1", "DENV-2: 6 mon. MRCA"="cornflowerblue", "DENV-2: 9 mon. MRCA"="blue", "DENV-2: 12 mon. MRCA"="navy")

colzC=c("DENV-1: 1 mon" = "lightgreen" , "DENV-1: 3 mon" = "green2", "DENV-1: 6 mon"="green3", "DENV-1: 12 mon"="forestgreen",
        "DENV-2: 1 mon"="lightblue1", "DENV-2: 3 mon"="cornflowerblue", "DENV-2: 6 mon"="blue", "DENV-2: 12 mon"="navy")


pS20 <- ggplot(data=out.prop) + theme_bw()+ #facet_wrap(~transchain_threshold) +
  facet_grid(~DENV.serotype) +  
  #geom_line(aes(x=distance, y=prop, color=sex),show.legend = F) +
  geom_ribbon(aes(x=distance, ymin=prop_lci, ymax=prop_uci, fill=differentiate, group=differentiate), alpha=.3) +
  geom_line(aes(x=distance, y=prop, color=differentiate, group=differentiate)) +
  theme(panel.grid = element_blank(), axis.title = element_text(size=16), 
        #axis.title.x = element_blank(), axis.text.x = element_blank(), 
        #axis.ticks.x = element_blank(),
        plot.margin = unit(c(.2,.3,.1,.5), "cm"),
        #panel.spacing = unit(c(.3), "cm"),
        #legend.title = element_blank(),
        strip.background = element_rect(fill="white"), strip.text = element_text(size = 18),
        axis.text = element_text(size=14), legend.text = element_text(size=10),
        legend.title = element_text(size=11),
        legend.position = c(.19,.8)) + coord_cartesian(ylim=c(0,1), xlim=c(0,5), expand = F)+
  scale_color_manual(values=colzC, name="MRCA threshold") + scale_x_continuous(breaks=c(1,2,3,4,5)) +
  scale_fill_manual(values=colzC, name="MRCA threshold") + 
  #scale_color_manual(values=colz, name = "transmission chain\nthreshold (yrs)") + 
  #scale_fill_manual(values=colz, name = "transmission chain\nthreshold (yrs)") + 
  ylab("Proportion same transmission chain")  +
  xlab("Max distance between cases (km)") + guides(fill=guide_legend(ncol=2), color=guide_legend(ncol=2))


ggsave(file = paste0(homewd, "/final-figures/FigS20.png"),
       plot= pS20,
       units="mm",  
       width=80, 
       height=40, 
       scale=3, 
       dpi=300)


