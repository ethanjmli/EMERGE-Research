library(apLCMS)
library(readr)
library(readxl)
library(dplyr)
library(xMSanalyzer)
library(stringr)
setwd("/Users/ethan_j_li/Documents/Emory/Research")

hilic_results <- read_delim("/Users/ethan_j_li/Documents/Emory/Research/2017/ComBat_mzcalibrated_untargeted_averaged_featuretable_hilic.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
c18_results <- read_delim("/Users/ethan_j_li/Documents/Emory/Research/2017/ComBat_mzcalibrated_untargeted_averaged_featuretable_c18.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)


#####HILIC 2017 S1 Annotations
######
HRef<- read.csv("/Users/ethan_j_li/Documents/Emory/Research/Hilic total_confirmed.csv")

setwd("/Users/ethan_j_li/Documents/Emory/Research/2017/HILIC 2017 PreComp Regressions")

list.files(pattern="h.*s1_2017.csv$") 
list.filenames <- list.files(pattern="h.*s1_2017.csv$")
list.filenames

annotation_list <- list()
# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  annotation_list[[i]] <- read.csv(list.filenames[i])
}
# add the names of these files to the list
names(annotation_list) <- list.filenames
str(annotation_list)

setwd("/Users/ethan_j_li/Documents/Emory/Research/2017/HILIC 2017 Annotations")

for(i in 1:length(annotation_list)){
  MWASHilicS1 <- hilic_results[,1:2]
  hilic_2017_s1<-annotation_list[[i]]
  MWASHilicS1$hilic_s1_unadjusted_p <- hilic_2017_s1[,5]
  MWASHilicS1$hilic_s1_pfdr <- hilic_2017_s1[,6] 
  Sig<-MWASHilicS1[,1:2] ##Column # for mz and retention time
  LVL1<-HRef[,5:6] ##Column # for mz and retention time
  masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
  LVL1.sig<-slice(MWASHilicS1,masteroverlap$index.A)
  LVL1.match<-slice(HRef,masteroverlap$index.B)
  LVL1total<-cbind(LVL1.match,LVL1.sig)
  #LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
  LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
  #write.csv(MWASh20s1p,"MWASs1p.csv")
  write.csv(LVL1total,paste0(str_sub(list.filenames[i],1,-5),"_MWASs1p_lvl1_2017.csv"))
  print(list.filenames[i])
}


#####HILIC 2020 S2 Annotations
######

setwd("/Users/ethan_j_li/Documents/Emory/Research/2017/HILIC 2017 PreComp Regressions")

list.files(pattern="h.*s2_2017.csv$") 
list.filenames <- list.files(pattern="h.*s2_2017.csv$")
list.filenames

annotation_list <- list()
# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  annotation_list[[i]] <- read.csv(list.filenames[i])
}
# add the names of these files to the list
names(annotation_list) <- list.filenames
str(annotation_list)

setwd("/Users/ethan_j_li/Documents/Emory/Research/2017/HILIC 2017 Annotations")

for(i in 1:length(annotation_list)){
  MWASHilicS2 <- hilic_results[,1:2]
  hilic_2017_s2<-annotation_list[[i]]
  MWASHilicS2$hilic_s2_unadjusted_p <- hilic_2017_s2[,5]
  MWASHilicS2$hilic_s2_pfdr <- hilic_2017_s2[,6] 
  Sig<-MWASHilicS2[,1:2] ##Column # for mz and retention time
  LVL1<-HRef[,5:6] ##Column # for mz and retention time
  masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
  LVL1.sig<-slice(MWASHilicS2,masteroverlap$index.A)
  LVL1.match<-slice(HRef,masteroverlap$index.B)
  LVL1total<-cbind(LVL1.match,LVL1.sig)
  #LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
  LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
  #write.csv(MWASh20s1p,"MWASs1p.csv")
  write.csv(LVL1total,paste0(str_sub(list.filenames[i],1,-5),"_MWASs2p_lvl1_2017.csv"))
  print(list.filenames[i])
}

####C18 Annotations
####
CRef <- read.csv("/Users/ethan_j_li/Documents/Emory/Research/c18 total_confirmed.csv")

#####C18 2020 S1 Annotations
######

setwd("/Users/ethan_j_li/Documents/Emory/Research/2017/C18 2017 PreComp Regressions")

list.files(pattern="c.*1_2017.csv$") 
list.filenames <- list.files(pattern="c.*1_2017.csv$")
list.filenames

annotation_list <- list()
# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  annotation_list[[i]] <- read.csv(list.filenames[i])
}
# add the names of these files to the list
names(annotation_list) <- list.filenames
str(annotation_list)

setwd("/Users/ethan_j_li/Documents/Emory/Research/2017/C18 2017 Annotations")

for(i in 1:length(annotation_list)){
  MWASC18S1 <- c18_results[,1:2]
  c18_2017_s1<-annotation_list[[i]]
  MWASC18S1$c18_s1_unadjusted_p <- c18_2017_s1[,5]
  MWASC18S1$c18_s1_pfdr <- c18_2017_s1[,6] 
  Sig<-MWASC18S1[,1:2] ##Column # for mz and retention time
  LVL1<-CRef[,5:6] ##Column # for mz and retention time
  masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
  LVL1.sig<-slice(MWASC18S1,masteroverlap$index.A)
  LVL1.match<-slice(CRef,masteroverlap$index.B)
  LVL1total<-cbind(LVL1.match,LVL1.sig)
  #LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
  LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
  #write.csv(MWASh20s1p,"MWASs1p.csv")
  write.csv(LVL1total,paste0(str_sub(list.filenames[i],1,-5),"_MWASs1p_lvl1_2017.csv"))
  print(list.filenames[i])
}

#####C18 2020 S2 Annotations
######

setwd("/Users/ethan_j_li/Documents/Emory/Research/2017/C18 2017 PreComp Regressions")

list.files(pattern="c.*2_2017.csv$") 
list.filenames <- list.files(pattern="c.*2_2017.csv$")
list.filenames

annotation_list <- list()
# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  annotation_list[[i]] <- read.csv(list.filenames[i])
}
# add the names of these files to the list
names(annotation_list) <- list.filenames
str(annotation_list)

setwd("/Users/ethan_j_li/Documents/Emory/Research/2017/C18 2017 Annotations")

for(i in 1:length(annotation_list)){
  MWASC18S2 <- c18_results[,1:2]
  c18_2017_s2<-annotation_list[[i]]
  MWASC18S2$c18_s2_unadjusted_p <- c18_2017_s2[,5]
  MWASC18S2$c18_s2_pfdr <- c18_2017_s2[,6]
  Sig<-MWASC18S2[,1:2] ##Column # for mz and retention time
  LVL1<-CRef[,5:6] ##Column # for mz and retention time
  masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
  LVL1.sig<-slice(MWASC18S2,masteroverlap$index.A)
  LVL1.match<-slice(CRef,masteroverlap$index.B)
  LVL1total<-cbind(LVL1.match,LVL1.sig)
  #LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
  LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
  #write.csv(MWASh20s1p,"MWASs1p.csv")
  write.csv(LVL1total,paste0(str_sub(list.filenames[i],1,-5),"_MWASs2p_lvl1_2017.csv"))
  print(list.filenames[i])
}



