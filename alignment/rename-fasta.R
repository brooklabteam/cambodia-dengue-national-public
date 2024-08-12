rm(list=ls())

library(plyr)
library(dplyr)
library(ggplot2)
library(seqinr)

homewd = "/Users/carabrook/Developer/cambodia-dengue-national"
setwd(paste0(homewd, "/alignment/"))


#load prename
denv1.pre <- read.fasta(file = "GenBankDENV1_prename.fasta", as.string = T, forceDNAtolower = F)
names(denv1.pre)
denv2.pre <- read.fasta(file = "GenBankDENV2_prename.fasta", as.string = T, forceDNAtolower = F)
names(denv2.pre)

#load the metadata for the sequence files
all.dat <- read.csv(file = paste0(homewd, "/gather-sequences/All_Seq_SE_Asia.csv"), header = T, stringsAsFactors = F)
head(all.dat)

all.dat$Collection_Date <- as.Date(all.dat$Collection_Date, format = "%m/%d/%y")

#replace the dates on our  cambodia sequences with the more accurate (more precise) ones from our database
new.dat <- read.csv(file = paste0(homewd, "/data/beasttree_metadata.csv"))

head(new.dat)
new.dat <- dplyr::select(new.dat, Accession, date)
new.dat$date <- as.Date(new.dat$date, format = "%m/%d/%y")
names(new.dat) <- c("Accession", "precise_date")

#and merge
all.dat <- merge(all.dat, new.dat, by = "Accession", all.x = T)
head(all.dat)
all.dat$Collection_Date[!is.na(all.dat$precise_date)] <- all.dat$precise_date[!is.na(all.dat$precise_date)] 


all.dat <- dplyr::select(all.dat, -(precise_date))
head(all.dat)
all.dat$Collection_Date[all.dat$Country=="Cambodia"]

#generate a beast appropriate name for each sequence
all.dat$beast_name <- paste0(all.dat$Accession, "_", all.dat$Collection_Date)
all.denv1 = subset(all.dat, Serotype=="DENV-1")
all.denv2 = subset(all.dat, Serotype=="DENV-2")

#now collect names from the sequences as is
name.dat.denv1 <- cbind.data.frame(name = names(denv1.pre))
name.dat.denv2 <- cbind.data.frame(name = names(denv2.pre))

#and lose the points
name.dat.denv1$Accession <- sapply(strsplit(gsub("\\.", "_", name.dat.denv1$name),split = "_"), function(x) x[[1]])
name.dat.denv2$Accession <- sapply(strsplit(gsub("\\.", "_", name.dat.denv2$name),split = "_"), function(x) x[[1]])
#name.dat.denv2$Accession <- NA
#name.dat.denv2$Accession[1:(length(name.dat.denv2$name)-63)] <- sapply(strsplit(gsub("\\.", "_", name.dat.denv2$name[1:(length(name.dat.denv2$name)-63)]),split = "_"), function(x) x[[1]])
#name.dat.denv2$Accession[is.na(name.dat.denv2$Accession)] <- name.dat.denv2$name[is.na(name.dat.denv2$Accession)] 

#and merge on accession
merge.denv1 <- merge(name.dat.denv1, all.denv1, by = "Accession", all.x = T, sort = F)
merge.denv2 <- merge(name.dat.denv2, all.denv2, by = "Accession", all.x = T, sort = F)

head(merge.denv1)
head(merge.denv2)

# dat.1 <- read.fasta(file = paste0(homewd, "/alignment/allDENV1_beast.fasta"), as.string = T, forceDNAtolower = F)
# names(dat.1)
# 
# dat.2 <- read.fasta(file = paste0(homewd, "/alignment/allDENV2_beast.fasta"), as.string = T, forceDNAtolower = F)
# names(dat.2)

#and write with beast name
write.fasta(denv1.pre, names=merge.denv1$beast_name,  file = paste0(homewd, "/alignment/allDENV1_beast.fasta"), as.string = T)
write.fasta(denv2.pre, names=merge.denv2$beast_name,  file = paste0(homewd, "/alignment/allDENV2_beast.fasta"), as.string = T)
