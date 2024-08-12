
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


#make three timetrees and one map all together
#read in tree from beast
#first, make map

#and load the metadata
homewd = "/Users/carabrook/Developer/cambodia-dengue-national-public"
#homewd= "/home/rstudio"
setwd(homewd)


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
        axis.title = element_blank()) +
  geom_rect(xmin = 103.8, xmax = 104.8, ymin = 11, ymax = 12.2, 
            fill = NA, colour = "black", size = 1)  +
  coord_sf(xlim=c(95,170), expand = T)




#and get submap of kampong speu witht the points of cambodia sequences
cam = sf::st_read(paste0(homewd, "/data/province-shape/khm_admbnda_adm1_gov_20181004.shp"))
#cam = sf::st_read(paste0(homewd, "/data/shapefile/provinces.shp"))
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

#now select those used for geospatial analysis
dat.plot = subset(dat.plot, !is.na(lat) & year >2018)#185 sequences
length(dat.plot$DENV.serotype[dat.plot$DENV.serotype=="DENV-1"]) #60 - because there are 3 DENV-1 sequences with no lat/long
length(dat.plot$DENV.serotype[dat.plot$DENV.serotype=="DENV-2"]) #120

dat.plot.1 = subset(dat.plot, DENV.serotype=="DENV-1")
dat.plot.2 = subset(dat.plot, DENV.serotype=="DENV-2")

# all.denv <- get.spatial.object(dat1=dat.plot, dist.thresh = 5000, denv.serotype = "all")
# denv.map <- all.denv[[1]]
# dat.plot.cluster <- all.denv[[2]]
# dat.plot.cluster <- dplyr::select(dat.plot.cluster,new_label, cluster_ID)

all.denv1 <- get.spatial.object(dat1=dat.plot.1, dist.thresh = 10000, denv.serotype = "all", perc_jitter = 0.05)
denv.map1 <- all.denv1[[1]]
dat.plot.cluster1 <- all.denv1[[2]]
dat.plot.cluster1 <- dplyr::select(dat.plot.cluster1,new_label, cluster_ID)

all.denv2 <- get.spatial.object(dat1=dat.plot.2, dist.thresh = 10000, denv.serotype = "all", perc_jitter = 0.05)
denv.map2 <- all.denv2[[1]]
dat.plot.cluster2 <- all.denv2[[2]]
dat.plot.cluster2 <- dplyr::select(dat.plot.cluster2,new_label, cluster_ID)

# #here, add cluster ID to the dataset
# dat.plot <- merge(dat.plot, dat.plot.cluster, by="new_label", all.x=TRUE)

#here, add cluster ID to the dataset
dat.plot.1 <- merge(dat.plot.1, dat.plot.cluster1, by="new_label", all.x=TRUE)
dat.plot.2 <- merge(dat.plot.2, dat.plot.cluster2, by="new_label", all.x=TRUE)


#head(dat.plot)

# create your own color palette 
colorz1 = sample(rainbow(length(unique(dat.plot.1$cluster_ID))), size=length(unique(dat.plot.1$cluster_ID)), replace=F)
names(colorz1) <- sort(unique(dat.plot.1$cluster_ID))

colorz2 = sample(rainbow(length(unique(dat.plot.2$cluster_ID))), size=length(unique(dat.plot.2$cluster_ID)), replace=F)
names(colorz2) <- sort(unique(dat.plot.2$cluster_ID))


# get map
pKPS.1 <- ggplot(cam) + geom_sf(fill="forestgreen", color="black", alpha=.4, size =.4) + 
  # theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) +
  geom_sf(data=denv.map1, aes(fill=cluster_ID), color="black",
          shape=21, size=2.5, stroke=1, show.legend = F) +
  scale_size_continuous(range=c(3,20)) +
  guides(fill="none") + scale_fill_manual(values=colorz1)


pKPS.2 <- ggplot(cam) + geom_sf(fill="navy", color="black",alpha=.4,  size =.4) + 
  # theme_bw() + 
  theme(panel.background = element_blank()) + 
  theme(panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        # panel.background = element_rect(fill=NULL, color = "white", size=3),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0),"lines")) +
  geom_sf(data=denv.map2, aes(fill=cluster_ID), color="black",
          shape=21, size=2.5, stroke=1, show.legend = F) +
  scale_size_continuous(range=c(3,20)) +
  guides(fill="none") + scale_fill_manual(values=colorz2)



#now we add in the trees

dat.clust.save.1 <- dplyr::select(dat.plot.1, tip_name, new_label, cluster_ID)
dat.clust.save.2 <- dplyr::select(dat.plot.2, tip_name, new_label, cluster_ID)

dat <- read.csv(file = paste0(homewd, "/data/beasttree_metadata.csv"), header = T, stringsAsFactors = F)

#dat = subset(dat, DENV.serotype=="DENV-2")
names(dat)
dat$date <- as.Date(dat$date, format = "%m/%d/%y")
#dat.plot = subset(dat, !is.na(lat) & DENV.serotype=="DENV-1" & year >2018) #49 sequences
dat.plot = dat
#make new labels
dat.plot$date <- as.Date(dat.plot$date)
dat.plot$accession_num <- c(unlist(sapply(strsplit(dat.plot$tip_name, "_"), function(x) x[[1]])))
head(dat.plot)
dat.plot$new_label = paste(dat.plot$accession_num, dat.plot$date, sep = " ")


#dat.plot  <- dplyr::select(dat.plot, new_label, unique.ID, NIH.ID, CZB.ID, beast_name, date, year, lat, long, accession_num, DENV.serotype)
head(dat.plot)

#now build the same trees as before, except select out the lineages of interest and color based on cluster instead
tree1 <- read.beast(file = paste0(homewd, "/BEAST-tree/denv-beast/beast-out/denv1/denv1Avg.tree"))
tree2 <- read.beast(file = paste0(homewd, "/BEAST-tree/denv-beast/beast-out/denv2/denv2Avg.tree"))

tree1@phylo$tip.label[tree1@phylo$tip.label=="DENS-OB-067_2021-05-29"] <- "PP470671_2021-05-29"
tree1@phylo$tip.label[tree1@phylo$tip.label=="DENS-OB-068_2021-05-29"] <- "PP470672_2021-05-29"
tree1@phylo$tip.label[tree1@phylo$tip.label=="DENS-OB-070_2021-05-30"] <- "PP470673_2021-05-30"
  
tree1dat <- cbind.data.frame(tip_name = tree1@phylo$tip.label)
#replace

tree2dat <- cbind.data.frame(tip_name = tree2@phylo$tip.label)


head(tree1dat)
head(tree2dat)

#and load the metadata
dat <- read.csv(file = paste0(homewd, "/data/beasttree_metadata.csv"), header = T, stringsAsFactors = F)
head(dat)


setdiff(tree1dat$tip_name, dat$tip_name) 
setdiff(tree2dat$tip_name, dat$tip_name) 
#replace

#check the format
dat$date <- as.Date(dat$date, format = "%m/%d/%y")


mrsd.denv1 <- max(dat$date[dat$DENV.serotype=="DENV-1"]) #"2021-05-30"
mrsd.denv2 <- max(dat$date[dat$DENV.serotype=="DENV-2"])#"2022-12-01"



node.tree1 <- MRCA(tree1, which(tree1@phylo$tip.label== "OL412678_2019-07-25" ),which(tree1@phylo$tip.label == "ON046271_2004-11-20"))

tree1@data$posterior[tree1@data$posterior>=.9] <- 1
tree1@data$posterior[tree1@data$posterior<.9] <- 0
tree1@data$posterior <- as.factor(tree1@data$posterior)

postz = c("black", "white")
names(postz) <- as.numeric(c("1", "0"))



#this is the large tree in panel A, Figure 4
pB1 <- ggtree(tree1, mrsd=mrsd.denv1, color="forestgreen")  + 
  #geom_cladelab(node=node.tree1, label="Genotype I", textcolor="seagreen", barcolor="seagreen", fontsize=6,
   #             offset =-37, angle=270, offset.text = -12, vjust=2, hjust=.5)  +
  theme_tree2() + coord_cartesian(xlim=c(1930,2030), ylim=c(0,400)) + 
  #geom_range(range='length_0.95_HPD', color='red', alpha=.6, size=2) +
  #geom_tiplab(size=1) +
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1, show.legend = F) +
  scale_x_continuous(breaks=c(2018, 2020,  2022))+
  scale_fill_manual(values=postz)+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8)) 


#and here is panel B, Figure 4


node.tree2 <- MRCA(tree2, which(tree2@phylo$tip.label== "PP411229_2022-06-01" ),which(tree2@phylo$tip.label == "MK506264_2007-06-15"))
node.tree2.1 <- MRCA(tree2, which(tree2@phylo$tip.label== "OQ426897_2019-03-06" ),which(tree2@phylo$tip.label == "KF921930_2002-07-31"))



tree2@data$posterior[tree2@data$posterior>=.9] <- 1
tree2@data$posterior[tree2@data$posterior<.9] <- 0
tree2@data$posterior <- as.factor(tree2@data$posterior)




pC2 <- ggtree(tree2, mrsd=mrsd.denv2, color="navy")  + theme_tree2() + 
  coord_cartesian(xlim=c(1930,2030),  ylim=c(0,450))+
  #geom_cladelab(node=node.tree2, label="Cosmopolitan I", textcolor="tomato",barcolor="tomato",
   #             offset =-37, angle=270, offset.text = -12, fontsize=6, vjust=2, hjust=.5)  +
  #geom_cladelab(node=node.tree2.1, label="Asian I", textcolor="navy", barcolor="navy", fontsize=6,vjust=2, hjust=.5,
   #             offset =-37, angle=270, offset.text = -12)  +
  #geom_tiplab(size=1)+
  geom_nodepoint(aes(fill=posterior), shape=21, color="black", size=1, stroke=.1, show.legend = F) +
  scale_fill_manual(values=postz)+
  scale_x_continuous(breaks=c(2018, 2020,  2022))+
  theme(legend.position = c(.2,.8), 
        legend.key.size = unit(.3, units="cm"),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))


tree1merge <- merge(x=tree1dat, y=dat, by="tip_name", all.x = T, sort = F)
tree2merge <- merge(x=tree2dat, y=dat, by="tip_name", all.x = T, sort = F)
head(tree1merge)
head(tree2merge)


tree1merge$new_label = sapply(strsplit(tree1merge$tip_name, "_"), function(x) x[[1]])
tree1merge$new_label <- paste0(tree1merge$new_label, " ", as.character(tree1merge$date))

tree1merge$new_seq = "no"
#tree1merge$new_seq[tree1merge$country=="Cambodia" & !is.na(tree1merge$sex)] <- "yes"
length(tree1merge$new_seq[tree1merge$new_seq=="yes"])#63
tree1merge$new_seq[tree1merge$country=="Cambodia" & !is.na(tree1merge$NIH.ID)] <- "yes"

tree1merge$new_seq <- as.factor(tree1merge$new_seq)

tree1merge$CambodiaSeq <- "no"
tree1merge$CambodiaSeq[tree1merge$country=="Cambodia"] <- "yes"


unique(tree2merge$country)
tree2merge$new_label = sapply(strsplit(tree2merge$tip_name, "_"), function(x) x[[1]])
tree2merge$new_label <- paste0(tree2merge$new_label, " ", as.character(tree2merge$date))

tree2merge$new_seq = "no"
tree2merge$new_seq[tree2merge$country=="Cambodia" & !is.na(tree2merge$NIH.ID)] <- "yes"
length(tree2merge$new_seq[tree2merge$new_seq=="yes"])#120
tree2merge$new_seq <- as.factor(tree2merge$new_seq)

tree2merge$CambodiaSeq <- "no"
tree2merge$CambodiaSeq[tree2merge$country=="Cambodia"] <- "yes"

shapez = c("yes"=24, "no"=21)



#and merge the data with the cluster ID
dat.clust.save.1 <- dplyr::select(dat.clust.save.1, -(new_label))
dat.clust.save.2 <- dplyr::select(dat.clust.save.2, -(new_label))
tree1merge <- merge(tree1merge, dat.clust.save.1, by = "tip_name", all.x = T)
tree2merge <- merge(tree2merge, dat.clust.save.2, by = "tip_name", all.x = T)

#
colorz2.1 = c("no"="white","yes"="black")
alphaz = c("no"=0,"yes"=1)


tree1.tiplabel<-tree1@phylo[["tip.label"]]



pB <- pB1 %<+% tree1merge +
  ggnewscale::new_scale_fill() +
  #geom_tiplab(size=.8, aes(color=new_seq,alpha=new_seq), nudge_x=-9.5, show.legend=F) +
  coord_cartesian(xlim=c(2017.5,2023), ylim=c(185,375)) +
  #geom_tiplab(size=1, aes(alpha=new_seq, color=new_seq), hjust=-3) +
  #geom_label(aes(label=new_label), label.size = NA, size=4, hjust=-.06) +
  #geom_tippoint(aes(shape=CambodiaSeq), color="black", fill="black", size = 4)+#,show.legend = F) +
  geom_tippoint(aes(fill=cluster_ID, 
                    shape=new_seq, 
                    alpha=new_seq, 
                    color=new_seq), 
                size = 2, show.legend = F, stroke=.5) +
  theme(legend.position = c(.1,.75), 
        legend.key.size = unit(.4, units="cm"),
        legend.text = element_text(size=8),
        legend.title = element_text(size=10),
        plot.margin =  unit(c(.3,.3,.3,1.8), units="cm"),
        #plot.background  = element_rect(size=3, fill = NULL, color = "black"),
        axis.text = element_text(size=18)) +
  scale_color_manual(values=colorz2.1) +
  scale_fill_manual(values=colorz1) +
  scale_shape_manual(values=shapez) +
  scale_alpha_manual(values=alphaz) 
#scale_color_manual(values=colorz)#+





#Fig 4C
pC <-pC2 %<+% tree2merge + 
  ggnewscale::new_scale_fill()  + 
  coord_cartesian(xlim=c(2017.5,2023), ylim=c(125,400)) +
  #geom_tiplab(size=.8, aes(color=new_seq, alpha=new_seq), nudge_x=-3, show.legend=F) +
  #geom_tiplab(size=1) +
  geom_tippoint(aes(fill=cluster_ID, shape=new_seq, alpha=new_seq, color=new_seq), 
                size = 2, show.legend = F, stroke=.5) +
  theme(axis.text = element_text(size=18), 
        #plot.background  = element_rect(size=3, fill = NULL, color = "black"),
        legend.position = "none",
        plot.margin =  unit(c(.3,.3,.3,1.8), units="cm"),
        legend.title = element_blank()) + 
  #coord_cartesian(c(2000,2021)) +
  scale_shape_manual(values=shapez) +scale_fill_manual(values=colorz2) + 
  scale_alpha_manual(values=alphaz) + scale_color_manual(values=colorz2.1) 




#and put it all together


pTop<- cowplot::plot_grid(pKPS.1, pB,  pC, labels=c("A", "C"), nrow=1, ncol=2, label_size = 22)

pBottom <- cowplot::plot_grid(pKPS.2, pC,  pC, labels=c("B", "D"), nrow=1, ncol=2, label_size = 22)




pALL <- cowplot::plot_grid(pTop, pBottom,  nrow=2, ncol = 1, rel_heights = c(1,1)) + theme_bw() + theme(panel.background = element_rect(fill="white"), plot.background = element_rect(fill="white")) 


ggsave(file = paste0(homewd, "/final-figures/FigS19.png"),
       plot= pALL,
       units="mm",  
       width=100, 
       height=90, 
       scale=3, 
       dpi=300)





####################################################################
####################################################################
####################################################################

#finally now we subset the tree to just the recent Cambodia clades for analysis
node2_sub <-MRCA(tree2, which(tree2@phylo$tip.label== "PP411228_2022-06-01" ),which(tree2@phylo$tip.label == "OL414746_2020-06-16"))

#take the upper DENV-1 clade (there is also a lower clade)
node1_sub <-MRCA(tree1, which(tree1@phylo$tip.label== "OL412703_2019-08-09"),which(tree1@phylo$tip.label == "GU131923_2005-07-31"))



#get
library(ape)
library(treeio)
library(ggplot2)
library(ggtree)
library(tidytree)

tree1sub <- tree_subset(tree=tree1, node =  node1_sub, levels_back=0)
tree2sub <- tree_subset(tree=tree2, node =  node2_sub, levels_back=0)


p1 <- ggtree(tree1sub) %<+% tree1merge + geom_tiplab() + 
  geom_tippoint(aes(color=country))

#get mrca for these 
nodedenv1 <- MRCA(tree1, which(tree1@phylo$tip.label== "OL412703_2019-08-09"),which(tree1@phylo$tip.label == "MW265679_2015-10-06"))
node.dat.denv1 <- pB1$data #from above
subset(node.dat.denv1, node==nodedenv1) #node height is 8.32  years
mrsd.denv1 - (8.32  *365) # "2013-02-04"
node.dat.denv1$height_0.95_HPD[[nodedenv1]] #7.754796 8.890605

mrsd.denv1 - (node.dat.denv1$height_0.95_HPD[[node2cosmomrca]][1]*365) # "2016-02-26"
mrsd.denv1 - (node.dat.denv1$height_0.95_HPD[[node2cosmomrca]][2]*365) # "2014-11-01"


#and the cosmo 2 sequences

#and mrca for cambodia with previously reported cambodian sequences
node2cosmomrca <-MRCA(tree2, which(tree2@phylo$tip.label== "PP411228_2022-06-01" ),which(tree2@phylo$tip.label == "FJ639698_2002-07-31"))

node.dat.denv2 <- pC$data #from above
subset(node.dat.denv2, node==node2cosmomrca) #node height is 77.1  years
mrsd.denv2 - (77.0  *365) # "1945-12-20"
node.dat.denv2$height_0.95_HPD[[node2cosmomrca]] #73.51945 80.88328

mrsd.denv2 - (node.dat.denv2$height_0.95_HPD[[node2cosmomrca]][1]*365) # "1949-06-12"
mrsd.denv2 - (node.dat.denv2$height_0.95_HPD[[node2cosmomrca]][2]*365) # "1942-02-01"


#and look at age distribution of cases 
head(df.sum)

#get the sum of cases per year
#2019: 94 genotypes. 59.6% (56/94) are DENV-1 and 37.2% (35/94) are DENV-2. of those DENV-2, 74.3% (26/35) are Cosmo
#2020: 32 genotypes. 18.8% (6/32) are DENV-1 and 78.1% (25/32) are DENV-2. of those DENV-2, 92% (23/25) are Cosmo
#2021: 28 genotypes. 10.7% (3/28) are DENV-1 and 89.3% (25/28) are DENV-2. of those DENV-2, 100% are Cosmo
#2022: 118 genotypes. 0.85% (1/118) are DENV-1 and 91.5% (108/118) are DENV-2. of those DENV-2, 99.1%% (107/108) are Cosmo


