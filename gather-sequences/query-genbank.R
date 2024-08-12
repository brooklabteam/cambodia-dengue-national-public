rm(list=ls())

library(dplyr)
library(plyr)
library(lubridate)

#set to your home
homewd = "/Users/carabrook/Developer/cambodia-dengue-national-public"
dat <- read.csv(file=paste0(homewd, "/gather-sequences/allSEasiadengue.csv"), header = T, stringsAsFactors = F)
head(dat)
names(dat)
names(dat)[names(dat)=="Organism_Name"] <- "Serotype"
names(dat)[names(dat)=="Geo_Location"] <- "Country"
dat$Serotype[dat$Serotype=="dengue virus type 1" |dat$Serotype=="dengue virus type I" | dat$Serotype=="dengue virus 1" ] <- "DENV-1"
dat$Serotype[dat$Serotype=="Dengue virus type 2" | dat$Serotype=="dengue virus type 2" |dat$Serotype=="dengue virus 2" ] <- "DENV-2"
unique(dat$Serotype)
#choose
dat <- dplyr::select(dat, Accession, Serotype, Country, Collection_Date, Host)
head(dat)
unique(dat$Host)
#unique(dat$Geo_Location)
#dat$Country <- sapply(strsplit(dat$Geo_Location, split=":"), "[",1)
unique(dat$Country)
#dat$Country[dat$Country=="Viet Nam"] <- "Vietnam"

#and edit dates
unique(dat$Collection_Date)
dat$Collection_Date[dat$Collection_Date==""] <- NA
dat <- dat[!is.na(dat$Collection_Date),]
unique(dat$Collection_Date)
dat$nchar_date <- nchar(dat$Collection_Date)
unique(dat$Collection_Date[dat$nchar_date==4])
dat$Collection_Date[dat$nchar_date==4] <- paste0(dat$Collection_Date[dat$nchar_date==4], "-07-31")
unique(dat$Collection_Date[dat$nchar_date==6] )
dat$Collection_Date[dat$nchar_date==6] <- as.character(as.Date(dat$Collection_Date[dat$nchar_date==6], format = "%m/%d/%y"))
unique(dat$nchar_date)
dat$has_dash <- lapply(strsplit(dat$Collection_Date, split="/"), length)
dat$Collection_Date[dat$nchar_date==7 & dat$has_dash==3] <- as.character(as.Date(dat$Collection_Date[dat$nchar_date==7 & dat$has_dash==3], format = "%m/%d/%y"))
dat$Collection_Date[dat$nchar_date==7 & dat$has_dash==1] <- paste0(dat$Collection_Date[dat$nchar_date==7 & dat$has_dash==1], "-15")
dat$Collection_Date[dat$nchar_date==8] <- as.character(as.Date(dat$Collection_Date[dat$nchar_date==8], format = "%m/%d/%y"))
unique(dat$Collection_Date)
dat$Collection_Date <- as.Date(dat$Collection_Date )

#now only past 2000
dat = subset(dat, Collection_Date >="2000-01-01")
head(dat)
tail(dat)

#then, split and collect three genotypes per year at random from each country
dat$year <- year(dat$Collection_Date)

dat.cambodia <- subset(dat, Country=="Cambodia")
dat <- subset(dat, Country!="Cambodia")
unique(dat$Serotype)
dat <- subset(dat, !is.na(dat$Serotype))
dat.split <- dlply(dat, .(Serotype, Country, year))

select.genomes <- function(df){
  
  if(nrow(df)>3){
   subseq <- sample(df$Accession, size=3, replace = F)
   
   df$flag <- 0
   for(i in 1:3){
     df$flag[df$Accession==subseq[i]] <- 1
   }
   df.out = subset(df, flag==1)
  df.out <- dplyr::select(df.out, -(flag))
  }else{
    df.out = df
  }
  return(df.out)
  
}

dat.out <- lapply(dat.split, select.genomes)
dat.sub <- data.table::rbindlist(dat.out)
head(dat.sub)


length(dat.cambodia$Accession[dat.cambodia$Serotype=="DENV-1"])#197
length(dat.cambodia$Accession[dat.cambodia$Serotype=="DENV-2"])#179
dat.sub <- rbind(dat.sub, dat.cambodia)
names(dat.sub)[names(dat.sub)=="year"] <- "Year"

dat.sub <- dplyr::select(dat.sub, Accession, Serotype, Collection_Date, Year, Country)
dat.sub <- arrange(dat.sub, Serotype, Country, Year)

#dat.sub.add = subset(dat.sub, Country=="Indonesia")
head(dat.sub)

#and remove the three duplicate cambodia sequences (DENV2)
dat.sub$flag <- 0
dat.sub$flag[dat.sub$Accession=="OQ000263"] <- 1
dat.sub$flag[dat.sub$Accession=="OL414731"] <- 1
dat.sub$flag[dat.sub$Accession=="OP999339"] <- 1

dat.sub = subset(dat.sub, flag==0)

dat.sub <- dplyr::select(dat.sub, -(flag))

write.csv(dat.sub, file = paste0(homewd, "/gather-sequences/All_Seq_SE_Asia.csv"), row.names = F)

#and get the piece to call from genbank
#separate for dengue 1 and 2
denv1 = subset(dat.sub, Serotype=="DENV-1")
denv2 = subset(dat.sub, Serotype=="DENV-2")

fileConn<-file(paste0(homewd,"/gather-sequences/DENV1-NCBI.txt"))
writeLines(paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=",paste(denv1$Accession, collapse = ",")), fileConn)
close(fileConn)

#and dengue 2
fileConn<-file(paste0(homewd,"/gather-sequences/DENV2-NCBI.txt"))
writeLines(paste0("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&retmode=text&id=",paste(denv2$Accession, collapse = ",")), fileConn)
close(fileConn)

#then 


