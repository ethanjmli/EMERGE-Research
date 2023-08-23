library(apLCMS)
library(readr)
library(readxl)
library(dplyr)
library(xMSanalyzer)
library(stringr)
setwd("/Users/ethan_j_li/Documents/Emory/Research")

hilic_results <- read_delim("/Users/ethan_j_li/Documents/Emory/Research/2020/ComBat_mzcalibrated_untargeted_mediansummarized_featuretable_hilic.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
c18_results <- read_delim("/Users/ethan_j_li/Documents/Emory/Research/2020/ComBat_mzcalibrated_untargeted_mediansummarized_featuretable_c18.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)
newcolnames <- c("unadjusted p", "fdr adjusted p")
####HILIC Annotations####

HRef<- read.csv("/Users/ethan_j_li/Documents/Emory/Research/Hilic total_confirmed.csv") #Read in HILIC reference

#####HILIC 2020 S1 Annotations######


setwd("/Users/ethan_j_li/Documents/Emory/Research/2020/HILIC PreComp Regressions 2020")
setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/HILIC 2020 Final Regressions")

list.files(pattern="h.*s1_Final.csv$") 
list.filenames <- list.files(pattern="h.*s1_Final.csv$")
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

setwd("/Users/ethan_j_li/Documents/Emory/Research/2020/HILIC 2020 Annotations")
setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/HILIC 2020 Final Annotations")

for(i in 1:length(annotation_list)){
  MWASHilicS1 <- hilic_results[,1:2] #get mz and rt for LCMS features
  hilic_2020_s1<-annotation_list[[i]] #load MWAS results
  MWASHilicS1$hilic_s1_unadjusted_p <- hilic_2020_s1[,5] #get unadjusted p value from MWAS
  MWASHilicS1$hilic_s1_pfdr <- hilic_2020_s1[,6] #get adjusted p value from MWAS
  names(MWASHilicS1)[3:4]<-newcolnames #prepare for rowbinding
  Sig<-MWASHilicS1[,1:2] ##Column # for mz and retention time from LCMS
  LVL1<-HRef[,5:6] ##Column # for mz and retention time from reference
  masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS') #Find overlap between 
                                                                                            #reference and LCMS
  LVL1.sig<-slice(MWASHilicS1,masteroverlap$index.A) #Extract LCMS features that overlapped 
  LVL1.match<-slice(HRef,masteroverlap$index.B) #Extract reference features that overlapped
  names(LVL1.match)[5]<-'ref_mz' #Distinguish LCMS mz and Reference mz
  LVL1total<-cbind(LVL1.match,LVL1.sig) 
  #LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
  LVL1total<-LVL1total[!duplicated(LVL1total[,2]),] #Remove duplicate detected features
  #write.csv(MWASh20s1p,"MWASs1p.csv")
  write.csv(LVL1total,paste0(str_sub(list.filenames[i],1,-5),"_MWASs1p_lvl1.csv"))
  print(list.filenames[i]) #Verify loop is running
}


#####HILIC 2020 S2 Annotations######

setwd("/Users/ethan_j_li/Documents/Emory/Research/2020/HILIC PreComp Regressions 2020")
setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT/HILIC 2020 HFT Regression")

list.files(pattern="h.*s2_Final.csv$") 
list.filenames <- list.files(pattern="h.*s2_Final.csv$")
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

setwd("/Users/ethan_j_li/Documents/Emory/Research/2020/HILIC 2020 Annotations")
setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/HILIC 2020 Final Annotations")
list.files()

for(i in 1:length(annotation_list)){
  MWASHilicS2 <- hilic_results[,1:2]
  hilic_2020_s2<-annotation_list[[i]]
  MWASHilicS2$hilic_s2_unadjusted_p <- hilic_2020_s2[,5]
  MWASHilicS2$hilic_s2_pfdr <- hilic_2020_s2[,6]
  names(MWASHilicS2)[3:4]<-newcolnames
  Sig<-MWASHilicS2[,1:2] ##Column # for mz and retention time
  LVL1<-HRef[,5:6] ##Column # for mz and retention time
  masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
  LVL1.sig<-slice(MWASHilicS2,masteroverlap$index.A)
  LVL1.match<-slice(HRef,masteroverlap$index.B)
  names(LVL1.match)[5]<-'ref_mz'
  LVL1total<-cbind(LVL1.match,LVL1.sig)
  #LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
  LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
  #write.csv(MWASh20s1p,"MWASs1p.csv")
  write.csv(LVL1total,paste0(str_sub(list.filenames[i],1,-5),"_MWASs2p_lvl1.csv"))
  print(list.filenames[i])
}


####C18 Annotations####
####

CRef <- read.csv("/Users/ethan_j_li/Documents/Emory/Research/c18 total_confirmed.csv") #Read in C18 Reference

#####C18 2020 S1 Annotations####


setwd("/Users/ethan_j_li/Documents/Emory/Research/2020/C18 PreComp Regressions 2020")
setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/C18 2020 Final Regressions")

list.files(pattern="c.*s1_Final.csv$") 
list.filenames <- list.files(pattern="c.*s1_Final.csv$")
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

setwd("/Users/ethan_j_li/Documents/Emory/Research/2020/C18 2020 Annotations")
setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/C18 2020 Final Annotations")


for(i in 1:length(annotation_list)){
  MWASC18S1 <- c18_results[,1:2]
  c18_2020_s1<-annotation_list[[i]]
  MWASC18S1$c18_s1_unadjusted_p <- c18_2020_s1[,5]
  MWASC18S1$c18_s1_pfdr <- c18_2020_s1[,6]
  names(MWASC18S1)[3:4]<-newcolnames
  Sig<-MWASC18S1[,1:2] ##Column # for mz and retention time
  LVL1<-CRef[,5:6] ##Column # for mz and retention time
  masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
  LVL1.sig<-slice(MWASC18S1,masteroverlap$index.A)
  LVL1.match<-slice(CRef,masteroverlap$index.B)
  names(LVL1.match)[5]<-'ref_mz'
  LVL1total<-cbind(LVL1.match,LVL1.sig)
  #LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
  LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
  #write.csv(MWASh20s1p,"MWASs1p.csv")
  write.csv(LVL1total,paste0(str_sub(list.filenames[i],1,-5),"_MWASs1p_lvl1.csv"))
  print(list.filenames[i])
}

#####C18 2020 S2 Annotations####


setwd("/Users/ethan_j_li/Documents/Emory/Research/2020/C18 PreComp Regressions 2020")
setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT/C18 2020 HFT Regression")

list.files(pattern="c.*s2_Final.csv$") 
list.filenames <- list.files(pattern="c.*s2_Final.csv$")
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

setwd("/Users/ethan_j_li/Documents/Emory/Research/2020/C18 2020 Annotations")
setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/C18 2020 Final Annotations")

for(i in 1:length(annotation_list)){
  MWASC18S2 <- c18_results[,1:2]
  c18_2020_s2<-annotation_list[[i]]
  MWASC18S2$c18_s2_unadjusted_p <- c18_2020_s2[,5]
  MWASC18S2$c18_s2_pfdr <- c18_2020_s2[,6]
  names(MWASC18S2)[3:4]<-newcolnames
  Sig<-MWASC18S2[,1:2] ##Column # for mz and retention time
  LVL1<-CRef[,5:6] ##Column # for mz and retention time
  masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
  LVL1.sig<-slice(MWASC18S2,masteroverlap$index.A)
  LVL1.match<-slice(CRef,masteroverlap$index.B)
  names(LVL1.match)[5]<-'ref_mz'
  LVL1total<-cbind(LVL1.match,LVL1.sig)
  #LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
  LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
  #write.csv(MWASh20s1p,"MWASs1p.csv")
  write.csv(LVL1total,paste0(str_sub(list.filenames[i],1,-5),"_MWASs2p_lvl1.csv"))
  print(list.filenames[i])
}

####Making combined HILIC s1 s2 files and filtering on unadjusted p <0.01####

setwd("/Users/ethan_j_li/Documents/Emory/Research/2020/HILIC 2020 Annotations")
setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/HILIC 2020 Final Annotations")
setwd("C:/Users/EJLI4/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/HILIC 2020 Final Annotations")

list.files()
list.files(pattern="h.*lvl1.csv$")
list.filenames2<-list.files(pattern="h.*lvl1.csv$")
list.filenames2

combined_list <- list()
for (i in 1:length(list.filenames2))
{
  combined_list[[i]] <- read.csv(list.filenames2[i])
}

names(combined_list) <- list.filenames2
str(combined_list)

for (i in 1:5){
  a<-(i+(i-1))
  print(a)
  b<-(i+i)
  print(b)
  print(list.filenames2[a])
  print(list.filenames2[b])
  combined_list[[a]]$origin <- str_sub(list.filenames2[a],1,-27)
  combined_list[[b]]$origin <- str_sub(list.filenames2[b],1,-27)
  test <- rbind(combined_list[[a]],combined_list[[b]])
  sigtest <-test[test[,10]<0.01 & !is.na(test[,10]),]
  #test <- test[order(test[,3],test[,10],decreasing=FALSE),]
  sigtest<-sigtest[!duplicated(sigtest[,2]),] 
  
  write.csv(sigtest,paste0(str_sub(list.filenames2[a],1,-27),"_combined_annotation.csv"))
}


list.files(pattern="combined.*.csv")
list.filenames <- list.files(pattern="combined.*.csv")

rbind_list <- list()
# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  rbind_list[[i]] <- read.csv(list.filenames[i])
}
# add the names of these files to the list
names(rbind_list) <- list.filenames
str(rbind_list)

hilic_combined_annotation <- rbind(rbind_list[[1]],rbind_list[[2]],rbind_list[[3]],rbind_list[[4]],rbind_list[[5]])
write.csv(hilic_combined_annotation,"HILIC Final combo annotation.csv")

####Making combined C18 s1 s2 files and filtering on unadjusted p <0.01####

setwd("/Users/ethan_j_li/Documents/Emory/Research/2020/C18 2020 Annotations")
setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/C18 2020 Final Annotations")

list.files(pattern="c.*lvl1.csv$")
list.filenames2<-list.files(pattern="c.*lvl1.csv$")
list.filenames2

combined_list <- list()
for (i in 1:length(list.filenames2))
{
  combined_list[[i]] <- read.csv(list.filenames2[i])
}

names(combined_list) <- list.filenames2
str(combined_list)

for (i in 1:5){
  a<-(i+(i-1))
  print(a)
  b<-(i+i)
  print(b)
  print(list.filenames2[a])
  print(list.filenames2[b])
  test <- rbind(combined_list[[a]],combined_list[[b]])
  test$origin <- str_sub(list.filenames2[a],1,-27)
  sigtest <-test[test[,10]<0.01 & !is.na(test[,10]),]
  sigtest<-sigtest[!duplicated(sigtest[,2]),] 
  write.csv(sigtest,paste0(str_sub(list.filenames2[a],1,-27),"_combined_annotation.csv"))
}

list.files(pattern="combined.*.csv")
list.filenames <- list.files(pattern="combined.*.csv")

rbind_list <- list()
# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  rbind_list[[i]] <- read.csv(list.filenames[i])
}
# add the names of these files to the list
names(rbind_list) <- list.filenames
str(rbind_list)

C18_combined_annotation <- rbind(rbind_list[[1]],rbind_list[[2]],rbind_list[[3]],rbind_list[[4]],rbind_list[[5]])
write.csv(C18_combined_annotation,"C18 Final combo annotation.csv")

