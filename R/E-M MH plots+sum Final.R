library(stringr)

###MH Plots 
####Positive Mode####
setwd("/Users/ethan_j_li/Documents/Emory/Research/2020")
posMO<- read.delim(file = "ComBat_mzcalibrated_untargeted_mediansummarized_featuretable_hilic.txt")

setwd("/Users/ethan_j_li/Documents/Emory/Research/2020/HILIC PreComp Regressions 2020")
setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/HILIC 2020 Final Regressions")
list.files(pattern="h.*.csv$") 
list.filenames <- list.files(pattern="h.*.csv$")

lmfit_crude <-list()

for (i in 1:length(list.filenames))
{
  
  lmfit_crude[[i]]<- cbind(posMO[,1:2],read.csv(list.filenames[i]))
  
}
names(lmfit_crude) <- list.filenames
str(lmfit_crude)


## hilic p-values 0.05
setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/HILIC 2020 Final E-M Plots")

for (j in 1:length(lmfit_crude)) {
  
  png(file = paste0(substring(list.filenames[j],1,nchar(list.filenames[j])-4),"_p0.05.png"),width=550,height=400)
  par(family = "sans",mar = c(4,4,2,2))
  ## retention time
  xvec = lmfit_crude[[j]][,"time"]
  yvec = log10(lmfit_crude[[j]][,"P"])*(-1)
  efficients = lmfit_crude[[j]][,"Est"]
  ythresh = -log10(0.05)
  colorvec = c("red","black","blue","darkgreen")
  pchvec = c(21,24)
  xincrement = 50
  yincrement = 1
  background.points.col="gray50"
  
  min_lim<-min(0,xvec)
  max_val<-max(xvec)
  plot(xvec,yvec,xlab="Retention Time (s)",ylab="-log10(pvalue)",xaxt="n",yaxt="n",cex=0.4,cex.main=1.5,xlim = c(0,300))
  axis(1, at=seq(0,300, by=xincrement) , las=2,cex.axis = 1.5)
  axis(2, at=seq(0 , (max(yvec,na.rm = TRUE)+2), by=yincrement) , las=2,cex.axis = 1.5)
  points(xvec,yvec,col=background.points.col,bg=background.points.col,cex=0.4,pch=21)
  
  goodip <- which(yvec>=ythresh)
  
  posefficients <- which(efficients >0)
  negefficients <- which(efficients <0)
  
  for(i in intersect(goodip,posefficients)){
    points(xvec[i],yvec[i],col=colorvec[1],cex=1,pch=pchvec[1],bg=colorvec[1])
  }
  
  for(i in intersect(goodip,negefficients)){
    points(xvec[i],yvec[i],col=colorvec[3],cex=0.8,pch=pchvec[2],bg=colorvec[3])
  }
  
  if(length(goodip) > 0){
    abline(h = ythresh,col = "green",lty=2,lwd=0.8)
    text(50, ythresh - 0.4, "p-value < 0.05", col = "black",cex = 1.5)
  }
  
  dev.off()
  
}


## hilic p-values 0.01
setwd("/Users/ethan_j_li/Documents/Emory/Research/2020/HILIC 2020 E-M Plots")

for (j in 1:length(lmfit_crude)) {
  
  png(file = paste0(substring(list.filenames[j],1,nchar(list.filenames[j])-4),"_p0.01.png"),width=550,height=400)
  par(family = "sans",mar = c(4,4,2,2))
  ## retention time
  xvec = lmfit_crude[[j]][,"time"]
  yvec = log10(lmfit_crude[[j]][,"P"])*(-1)
  efficients = lmfit_crude[[j]][,"Est"]
  ythresh = -log10(0.01)
  colorvec = c("red","black","blue","darkgreen")
  pchvec = c(21,24)
  xincrement = 50
  yincrement = 1
  background.points.col="gray50"
  
  min_lim<-min(0,xvec)
  max_val<-max(xvec)
  plot(xvec,yvec,xlab="Retention Time (s)",ylab="-log10(pvalue)",xaxt="n",yaxt="n",cex=0.4,cex.main=1.5,xlim = c(0,300))
  axis(1, at=seq(0,300, by=xincrement) , las=2,cex.axis = 1.5)
  axis(2, at=seq(0 , (max(yvec,na.rm = TRUE)+2), by=yincrement) , las=2,cex.axis = 1.5)
  points(xvec,yvec,col=background.points.col,bg=background.points.col,cex=0.4,pch=21)
  
  goodip <- which(yvec>=ythresh)
  
  posefficients <- which(efficients >0)
  negefficients <- which(efficients <0)
  
  for(i in intersect(goodip,posefficients)){
    points(xvec[i],yvec[i],col=colorvec[1],cex=1,pch=pchvec[1],bg=colorvec[1])
  }
  
  for(i in intersect(goodip,negefficients)){
    points(xvec[i],yvec[i],col=colorvec[3],cex=0.8,pch=pchvec[2],bg=colorvec[3])
  }
  
  if(length(goodip) > 0){
    abline(h = ythresh,col = "green",lty=2,lwd=0.8)
    text(50, ythresh - 0.4, "p-value < 0.05", col = "black",cex = 1.5)
  }
  
  dev.off()
  
}


###MH Plots 
####Negative Mode####
setwd("/Users/ethan_j_li/Documents/Emory/Research/2020")
NegMO<- read.delim(file = "ComBat_mzcalibrated_untargeted_mediansummarized_featuretable_c18.txt")

setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/C18 2020 Final Regressions")
list.files(pattern="c.*.csv$") 
list.filenames <- list.files(pattern="c.*.csv$")

lmfit_crude <-list()

for (i in 1:length(list.filenames))
{
  
  lmfit_crude[[i]]<- cbind(NegMO[,1:2],read.csv(list.filenames[i]))
  
}
names(lmfit_crude) <- list.filenames
str(lmfit_crude)


## c18 p-values 0.05
setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/C18 2020 Final E-M Plots")

for (j in 1:length(lmfit_crude)) {
  
  png(file = paste0(substring(list.filenames[j],1,nchar(list.filenames[j])-4),"_p0.05.png"),width=550,height=400)
  par(family = "sans",mar = c(4,4,2,2))
  ## retention time
  xvec = lmfit_crude[[j]][,"time"]
  yvec = log10(lmfit_crude[[j]][,"P"])*(-1)
  efficients = lmfit_crude[[j]][,"Est"]
  ythresh = -log10(0.05)
  colorvec = c("red","black","blue","darkgreen")
  pchvec = c(21,24)
  xincrement = 50
  yincrement = 1
  background.points.col="gray50"
  
  min_lim<-min(0,xvec)
  max_val<-max(xvec)
  plot(xvec,yvec,xlab="Retention Time (s)",ylab="-log10(pvalue)",xaxt="n",yaxt="n",cex=0.4,cex.main=1.5,xlim = c(0,300))
  axis(1, at=seq(0,300, by=xincrement) , las=2,cex.axis = 1.5)
  axis(2, at=seq(0 , (max(yvec,na.rm = TRUE)+2), by=yincrement) , las=2,cex.axis = 1.5)
  points(xvec,yvec,col=background.points.col,bg=background.points.col,cex=0.4,pch=21)
  
  goodip <- which(yvec>=ythresh)
  
  posefficients <- which(efficients >0)
  negefficients <- which(efficients <0)
  
  for(i in intersect(goodip,posefficients)){
    points(xvec[i],yvec[i],col=colorvec[1],cex=1,pch=pchvec[1],bg=colorvec[1])
  }
  
  for(i in intersect(goodip,negefficients)){
    points(xvec[i],yvec[i],col=colorvec[3],cex=0.8,pch=pchvec[2],bg=colorvec[3])
  }
  
  if(length(goodip) > 0){
    abline(h = ythresh,col = "green",lty=2,lwd=0.8)
    text(50, ythresh - 0.4, "p-value < 0.05", col = "black",cex = 1.5)
  }
  
  dev.off()
  
}



## c18 p-values 0.01
setwd("/Users/tanyouran/Desktop/Liang/Thesis Materials/E-M c18 plots")

for (j in 1:length(lmfit_crude)) {
  
  png(file = paste0(substring(list.filenames[j],1,nchar(list.filenames[j])-4),"_p0.01.png"),width=550,height=400)
  par(family = "sans",mar = c(4,4,2,2))
  ## retention time
  xvec = lmfit_crude[[j]][,"time"]
  yvec = log10(lmfit_crude[[j]][,"P"])*(-1)
  efficients = lmfit_crude[[j]][,"Est"]
  ythresh = -log10(0.01)
  colorvec = c("red","black","blue","darkgreen")
  pchvec = c(21,24)
  xincrement = 50
  yincrement = 1
  background.points.col="gray50"
  
  min_lim<-min(0,xvec)
  max_val<-max(xvec)
  plot(xvec,yvec,xlab="Retention Time (s)",ylab="-log10(pvalue)",xaxt="n",yaxt="n",cex=0.4,cex.main=1.5,xlim = c(0,300))
  axis(1, at=seq(0,300, by=xincrement) , las=2,cex.axis = 1.5)
  axis(2, at=seq(0 , (max(yvec,na.rm = TRUE)+2), by=yincrement) , las=2,cex.axis = 1.5)
  points(xvec,yvec,col=background.points.col,bg=background.points.col,cex=0.4,pch=21)
  
  goodip <- which(yvec>=ythresh)
  
  posefficients <- which(efficients >0)
  negefficients <- which(efficients <0)
  
  for(i in intersect(goodip,posefficients)){
    points(xvec[i],yvec[i],col=colorvec[1],cex=1,pch=pchvec[1],bg=colorvec[1])
  }
  
  for(i in intersect(goodip,negefficients)){
    points(xvec[i],yvec[i],col=colorvec[3],cex=0.8,pch=pchvec[2],bg=colorvec[3])
  }
  
  if(length(goodip) > 0){
    abline(h = ythresh,col = "green",lty=2,lwd=0.8)
    text(50, ythresh - 0.4, "p-value < 0.05", col = "black",cex = 1.5)
  }
  
  dev.off()
  
}
 

###Making TXT files for Mummichog Analyses
###Positive Mode

## import multiple csv files to a list
setwd("/Users/tanyouran/Desktop/Liang/Thesis Materials/E-M hillic")
list.files(pattern="E-M.*.csv$") 
list.filenames <- list.files(pattern="E-M.*.csv$")
list.filenames
# create an empty list that will serve as a container to receive the incoming files
lmfit_crude <- list()
# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  lmfit_crude[[i]] <- read.csv(list.filenames[i])
}
# add the names of these files to the list
names(lmfit_crude) <- list.filenames
str(lmfit_crude)

library(tidyverse)
polname<-str_remove_all(list.filenames, ".csv")


for (i in 1:8)
{
  print(i)
  print(polname[i])
  tryCatch(
    {
      ACE_ref<-posMO[1:2]
      ACE_ref$p.value<-lmfit_crude[[i]][4]
      ACE_ref$t.value<-lmfit_crude[[i]][7]
      ACE_ref<-ACE_ref[complete.cases(ACE_ref[ , 3]),]
      ACE_ref2<-matrix(as.numeric(unlist(ACE_ref)),nrow=nrow(ACE_ref))
      colnames(ACE_ref2)<-c("mz","time","p.value","t.value")
      
      write.table(ACE_ref2, paste0("/Users/tanyouran/Desktop/Liang/Thesis Materials/Mummichog_E-M/",polname[i],".txt"),sep="\t",row.names=FALSE)
    },error=function(e){})
}


###Making TXT files for Mummichog Analyses
###Negative Mode

## import multiple csv files to a list
setwd("/Users/tanyouran/Desktop/Liang/Thesis Materials/E-M c18")
list.files(pattern="E-M.*.csv$") 
list.filenames <- list.files(pattern="E-M.*.csv$")
list.filenames

# create an empty list that will serve as a container to receive the incoming files
lmfit_crude <- list()
# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  lmfit_crude[[i]] <- read.csv(list.filenames[i])
}
# add the names of these files to the list
names(lmfit_crude) <- list.filenames
str(lmfit_crude)


polname<-str_remove_all(list.filenames, ".csv")

for (i in 1:8)
{
  print(i)
  print(polname[i])
  tryCatch(
    {
      ACE_ref<-NegMO[1:2]
      ACE_ref$p.value<-lmfit_crude[[i]][4]
      ACE_ref$t.value<-lmfit_crude[[i]][7]
      ACE_ref<-ACE_ref[complete.cases(ACE_ref[ , 3]),]
      ACE_ref2<-matrix(as.numeric(unlist(ACE_ref)),nrow=nrow(ACE_ref))
      colnames(ACE_ref2)<-c("mz","time","p.value","t.value")
      
      write.table(ACE_ref2, paste0("/Users/tanyouran/Desktop/Liang/Thesis Materials/Mummichog_E-M/",polname[i],".txt"),sep="\t",row.names=FALSE)
    },error=function(e){})
}




library("stringr")
library("readxl")

####Summary of total&overlaping features####
# M-O Int 0.005
#list.filenames <- list.files(path = "/Users/tanyouran/Desktop/Liang/Thesis Materials/Mummichog_M-O_int(0.005)/Sum_Path", pattern="mcg_pathwayanalysis.*.xlsx$",full.names = TRUE, recursive = FALSE) 
list.filenames <- list.files(path = "/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/Final Heatmapping", pattern="mcg_pathwayanalysis.*.xlsx$",full.names = TRUE, recursive = FALSE) 
res<-substring(list.filenames,92,nchar(list.filenames)-5)
res

myData1 <- read_excel('/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/Final Heatmapping/mcg_pathwayanalysis_h_g_dm_s1_Final_2020.xlsx')[1:119,]
Mum_Pathway <- myData1[order(myData1$pathway),] [1]
Mum_Pathway_V1<-myData1[order(myData1$pathway),] [1]
Mum_Pathway_V2<-myData1[order(myData1$pathway),] [1]

# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  Path_Tem<- read_excel(list.filenames[i])[1:119,]
  Path_Tem<-Path_Tem[order(Path_Tem$pathway),]
  Mum_Pathway[1+i]<-Path_Tem$`p-value`
  Mum_Pathway_V1[1+i]<-Path_Tem$`overlap_size`  
  Mum_Pathway_V2[1+i]<-Path_Tem$`pathway_size`
  colnames(Mum_Pathway)[1+i]<-res[i]
  colnames(Mum_Pathway_V1)[1+i]<-res[i]
  colnames(Mum_Pathway_V2)[1+i]<-res[i]
  print(colnames(Mum_Pathway)[1+i])
  print(colnames(Mum_Pathway_V1[1+i]))
  print(colnames(Mum_Pathway_V2[1+i]))
}

# number of pathway with adjusted p.value<0.05
coln<-dim(Mum_Pathway)[2]

Mum_Pathway$Total0.05<-NA
Mum_Pathway$Total0.10<-NA
for(i in 1:119)
{
  Mum_Pathway$Total0.05[i]<-sum(Mum_Pathway[i,2:coln]<=0.05)
  Mum_Pathway$Total0.10[i]<-sum(Mum_Pathway[i,2:coln]<=0.1)
}

#average overlap features 
Mum_Pathway_V1$overlap_average<-0    
Mum_Pathway_V1 <- as.data.frame(Mum_Pathway_V1)

for(i in 1:nrow(Mum_Pathway_V1)){
  Mum_Pathway_V1$overlap_average[i] <- ceiling(mean(as.numeric(Mum_Pathway_V1[i,2:coln]), na.rm = TRUE))
}

#average total features
Mum_Pathway_V2$total_average<-0
Mum_Pathway_V2<-as.data.frame(Mum_Pathway_V2)
for(i in 1:nrow(Mum_Pathway_V2)){
  Mum_Pathway_V2$total_average[i] <-ceiling(mean(as.numeric(Mum_Pathway_V2[i,2:coln]),na.rm=TRUE))
}

#final dataset
Mum_Pathway<-cbind(Mum_Pathway,Mum_Pathway_V1[,"overlap_average"],Mum_Pathway_V2[,"total_average"])
Mum_Pathway_V3<-Mum_Pathway[order(-Mum_Pathway$Total0.05,-Mum_Pathway$Total0.1),]
names(Mum_Pathway_V3)[24:25]<-paste(c("number of overlapping features","number of metabolites in pathway"))
write.csv(Mum_Pathway_V3,"/Users/ethan_j_li/OneDrive - Emory University/Research/2020 HFT+Sub+Gast/Final Heatmapping/Pathway Summary Final PreComp 2020.csv")





# Summary of total&overlaping features
# M-O NoInt 0.005
list.filenames <- list.files(path = "/Users/tanyouran/Desktop/Liang/Thesis Materials/Mummichog_M-O_noint(0.005)/Sum_Path", pattern="mcg_pathwayanalysis.*.xlsx$",full.names = TRUE, recursive = FALSE) 
res<-substring(list.filenames,105,nchar(list.filenames)-5)
res

myData1 <- read_excel('/Users/tanyouran/Desktop/Liang/Thesis Materials/Mummichog_M-O_noint(0.005)/Sum_Path/mcg_pathwayanalysis_cbirthganoint_stg1.xlsx')[1:119,]
Mum_Pathway <- myData1[order(myData1$pathway),] [1]
Mum_Pathway_V1<-myData1[order(myData1$pathway),] [1]
Mum_Pathway_V2<-myData1[order(myData1$pathway),] [1]

# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  Path_Tem<- read_excel(list.filenames[i])[1:119,]
  Path_Tem<-Path_Tem[order(Path_Tem$pathway),]
  Mum_Pathway[1+i]<-Path_Tem$`p-value`
  Mum_Pathway_V1[1+i]<-Path_Tem$`overlap_size`  
  Mum_Pathway_V2[1+i]<-Path_Tem$`pathway_size`
  colnames(Mum_Pathway)[1+i]<-res[i]
  colnames(Mum_Pathway_V1)[1+i]<-res[i]
  colnames(Mum_Pathway_V2)[1+i]<-res[i]
  print(colnames(Mum_Pathway)[1+i])
  print(colnames(Mum_Pathway_V1[1+i]))
  print(colnames(Mum_Pathway_V2[1+i]))
}

# number of pathway with adjusted p.value<0.05
coln<-dim(Mum_Pathway)[2]

Mum_Pathway$Total0.05<-NA
Mum_Pathway$Total0.10<-NA
for(i in 1:119)
{
  Mum_Pathway$Total0.05[i]<-sum(Mum_Pathway[i,2:coln]<=0.05)
  Mum_Pathway$Total0.10[i]<-sum(Mum_Pathway[i,2:coln]<=0.1)
}

#average overlap features 
Mum_Pathway_V1$overlap_average<-0    
Mum_Pathway_V1 <- as.data.frame(Mum_Pathway_V1)

for(i in 1:nrow(Mum_Pathway_V1)){
  Mum_Pathway_V1$overlap_average[i] <- ceiling(mean(as.numeric(Mum_Pathway_V1[i,2:coln]), na.rm = TRUE))
}

#average total features
Mum_Pathway_V2$total_average<-0
Mum_Pathway_V2<-as.data.frame(Mum_Pathway_V2)
for(i in 1:nrow(Mum_Pathway_V2)){
  Mum_Pathway_V2$total_average[i] <-ceiling(mean(as.numeric(Mum_Pathway_V2[i,2:coln]),na.rm=TRUE))
}

#final dataset
Mum_Pathway<-cbind(Mum_Pathway,Mum_Pathway_V1[,"overlap_average"],Mum_Pathway_V2[,"total_average"])
Mum_Pathway_V3<-Mum_Pathway[order(-Mum_Pathway$Total0.05,-Mum_Pathway$Total0.1),]
names(Mum_Pathway_V3)[16:17]<-paste(c("number of overlapping features","number of metabolites in pathway"))
write.csv(Mum_Pathway_V3,"/Users/tanyouran/Desktop/Liang/Thesis Materials/Mummichog_M-O_noint(0.005)/Sum_Path_M-O_noint.csv")






###Positive Mode
setwd("/Users/tanyouran/Desktop/Liang/Thesis Materials/DATA")
posMO<- read.table(file = "ComBat_mzcalibrated_untargeted_averaged_featuretable_hilic.txt",header = T)

## import multiple csv files to a list
setwd("/Users/tanyouran/Desktop/Liang/Thesis Materials/E-M hillic")
list.files(pattern="E-M.*.csv$") 
list.filenames <- list.files(pattern="E-M.*.csv$")


# create an empty list that will serve as a container to receive the incoming files
lmfit_crude <- list()
# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  lmfit_crude[[i]] <- read.csv(list.filenames[i])
}
# add the names of these files to the list
names(lmfit_crude) <- list.filenames
str(lmfit_crude)


###Check out the sig num of the features in positive mode

colnames(lmfit_crude[[8]][5])

MWAStats <- as.data.frame(matrix(nrow = 8, ncol = 7))
colnames(MWAStats)<-c("E-Int-stg","FDR0.05","FDR0.20","Raw0.0005", "Raw0.005","Raw0.01","Raw0.05")
for (i in 1:8)
{
  tryCatch(
    {
      MWAStats[i,1]<-names(lmfit_crude)[i]
      MWAStats[i,2]<-sum(lmfit_crude[[i]][5]<0.05,na.rm = TRUE)
      MWAStats[i,3]<-sum(lmfit_crude[[i]][5]<0.20,na.rm = TRUE)
      MWAStats[i,4]<-sum(lmfit_crude[[i]][4]<0.0005,na.rm = TRUE)
      MWAStats[i,5]<-sum(lmfit_crude[[i]][4]<0.005,na.rm = TRUE)
      MWAStats[i,6]<-sum(lmfit_crude[[i]][4]<0.01,na.rm = TRUE)
      MWAStats[i,7]<-sum(lmfit_crude[[i]][4]<0.05,na.rm = TRUE)
    },error=function(e){})
}

write.csv(MWAStats,"/Users/tanyouran/Desktop/Liang/Thesis Materials/E-M hillic/TotalStats_pos.csv")





################
###Negative Mode
setwd("/Users/tanyouran/Desktop/Liang/Thesis Materials/DATA")
NegMO<- read.table(file = "ComBat_mzcalibrated_untargeted_averaged_featuretable_c18.txt",header = T)

setwd("/Users/tanyouran/Desktop/Liang/Thesis Materials/E-M c18")
list.files(pattern="E-M.*.csv$") 
list.filenames <- list.files(pattern="E-M.*.csv$")


# create an empty list that will serve as a container to receive the incoming files
lmfit_crude <- list()
# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  lmfit_crude[[i]] <- read.csv(list.filenames[i])
}
# add the names of these files to the list
names(lmfit_crude) <- list.filenames
str(lmfit_crude)


###Check out the sig num of the features in positive mode

colnames(lmfit_crude[[8]][5])

MWAStats <- as.data.frame(matrix(nrow = 8, ncol = 7))
colnames(MWAStats)<-c("E-Int-stg","FDR0.05","FDR0.20","Raw0.0005", "Raw0.005","Raw0.01","Raw0.05")
for (i in 1:8)
{
  tryCatch(
    {
      MWAStats[i,1]<-names(lmfit_crude)[i]
      MWAStats[i,2]<-sum(lmfit_crude[[i]][5]<0.05,na.rm = TRUE)
      MWAStats[i,3]<-sum(lmfit_crude[[i]][5]<0.20,na.rm = TRUE)
      MWAStats[i,4]<-sum(lmfit_crude[[i]][4]<0.0005,na.rm = TRUE)
      MWAStats[i,5]<-sum(lmfit_crude[[i]][4]<0.005,na.rm = TRUE)
      MWAStats[i,6]<-sum(lmfit_crude[[i]][4]<0.01,na.rm = TRUE)
      MWAStats[i,7]<-sum(lmfit_crude[[i]][4]<0.05,na.rm = TRUE)
    },error=function(e){})
}

write.csv(MWAStats,"/Users/tanyouran/Desktop/Liang/Thesis Materials/E-M c18/TotalStats_neg.csv")







###Making Scripts for Mummichog

setwd("/Users/tanyouran/Desktop/Liang/Thesis Materials/Mummichog_E-M")
list.files(pattern=".txt$") 
list.filenames <- list.files(pattern=".txt$")
list.filenames
list.filenames[1]

output<-unlist(substring(list.filenames[1],5,nchar(list.filenames)-4))[1]
test<-paste0("python /Users/tanyouran/Downloads/mummichog-1.0.9/mummichog/main.py -f /Users/tanyouran/Desktop/Liang/",list.filenames[1]," -o ",output, 
             "_0.005 -c 0.005 -m negative -z TRUE")
test

earthmum<-as.data.frame(matrix(nrow =8*2, ncol = 1))
###Neg @ 0.005

for (i in 1:8)
{
  output<-unlist(substring(list.filenames[i],5,nchar(list.filenames[i])-4))[1]
  earthmum[i,1]<-paste0("python /Users/tanyouran/Downloads/mummichog-1.0.9/mummichog/main.py -f /Users/tanyouran/Desktop/Liang/",list.filenames[i]," -o ",output, 
                        "_0.005 -c 0.005 -m negative -z TRUE")
}


###Pos @ 0.005
output<-unlist(substring(list.filenames[9],5,nchar(list.filenames[9])-4))[1]
test<-paste0("python /Users/tanyouran/Downloads/mummichog-1.0.9/mummichog/main.py -f /Users/tanyouran/Desktop/Liang/",list.filenames[9]," -o ",output, 
             "_0.005 -c 0.005 -m positive -z TRUE")
test

for (i in 9:16)
{
  output<-unlist(substring(list.filenames[i],5,nchar(list.filenames[i])-4))[1]
  earthmum[i,1]<-paste0("python /Users/tanyouran/Downloads/mummichog-1.0.9/mummichog/main.py -f /Users/tanyouran/Desktop/Liang/",list.filenames[i]," -o ",output, 
                        "_0.005 -c 0.005 -m positive -z TRUE")
}



write.table(earthmum, "earthmum.txt",sep="\t",row.names=FALSE)





rm(list = ls())

