############################################################################
############# EH590R Lecture 11 Metabolomics Part 1  #######################
############# Donghai Liang, 04/05/2021             #######################
###########################################################################

# Clean up memory
rm(list=ls())

# If you have not installed package lmerTest before
install.packages("lmerTest")

# below are the packages required for this session
library(lme4)
library(lmerTest)   ###These two packages are needed to execute linear mixed effect model
library(readr)

################## Step 1 Metabolomics Analytic Database Assembly ###########################

################## Note that here we skip an important data preprocessing step (aka, metabolomic feature extraction)
################## The feature extraction step is very time consuming, which will take days/weeks to complete
################## Here are the papers on the feature extraction packages we applied on this dataset
################## PMID: 19414529; PMID: 23323971; and PMID: 32807888
################## If you are interested in learning more about this, feel free to contact me afterwards

#######Step 1.1 Import and Clean the Metabolomic Feature Table###########
setwd("C:/Users/SSBUCKE/OneDrive - Emory University/MWAS Materials for Susan")

##Import the feature table
mft<- read_csv("ComBat_corrected_p1_U_p2_feature_table_DRIVE_hilic_plasma.csv")

##Check the dimension of the dataset
dim(mft) #checks how many variables
View(mft) 
###In this dataset, each row is a metabolomic feature-- hence, there are 13,419 features in this dataset
###Starting from column #10, each column represents the signal intensity (relative concentrations)for 
###each metabolic features on a biological sample

#check the first 10 columns
summary(mft[,1:10])

##Check data quality
##NumPres.All.Samples/NumPres.Biological.Samples/median_CV/Q.score
hist(mft$NumPres.All.Samples)
sum(mft$NumPres.All.Samples<29)  ##570*5%
#Checking how many samples are present in at least 15-20%.  ~4K features are present in less than 29 participants

hist(mft$median_CV)
sum(mft$median_CV>30)

plot(mft$Qscore,mft$median_CV)

##QAQC Filtering on metabolomics feature table
mdat<-subset(mft,mft$median_CV<=30 & mft$NumPres.All.Samples>29)
###data cleaning - rm data with median_CV greater than 30 and number of samples less than 29)

##check data quality on the new datset again
hist(mdat$NumPres.All.Samples)
sum(mdat$NumPres.All.Samples<29)  

hist(mdat$median_CV)
sum(mdat$median_CV>30)

#######Step 1.2 Calculate Mean Intensity for each metabolic feature across triplicates in Metabolomic Feature Table###########
View(mdat[,1:20])
###Import the sample sequence ID file
infofile <- plasma_sequencefileHILIC #used import dataset in Environment 
View(infofile)

###Merge Info file with feature table sample column names
colnames(mdat)[10:dim(mdat)[2]]
## create a vector of all of the column names 
colnames.mdat<-colnames(mdat)[10:dim(mdat)[2]]
View(colnames.mdat)
File.n <- as.data.frame(substr(colnames.mdat,1,nchar(colnames.mdat)-5)) #need to remove the extra char. to link
#substr function is very helpful with character variables and manipulating words 
View(File.n)
colnames(File.n)<-'File_Name' # Cleaning up the column title - want to use the same name as the file you are merging to
View(File.n)

##Merge two datasets and see if the sample info matches with the column names in the feature table
#need to do this to make sure the data is in the correct order
infofile.mdat <- merge.data.frame(File.n,infofile,by ='File_Name',all.x = T,sort = F)
View(infofile.mdat)
##If there is discrepancy, there will be NA in the Sample ID column
sum(is.na(infofile.mdat[,2]))

##Replacing the colnames in the feature table *aka, sequence ID) with sample ID
intensity <- as.matrix(mdat[,10:dim(mdat)[2]])
View(intensity)
colnames(intensity) <- infofile.mdat$Sample_ID #assigning sample ID names to the column names in intensity 
View(intensity)
###Replace missing values (0) with NA
intensity[intensity == 0] <- NA
View(intensity)
mdat1 <- cbind(mdat[,1:9],intensity)
View(mdat1)           ##Upon this step, you finish replacing the colnames in the feature table with sample ID, and replacing all missing values (0) with NA

##Start calculating means across triplicates
##Method-1 use for loops-- a bit time consuming
means_1<-data.frame(matrix(NA, nrow = dim(intensity)[1], ncol = dim(intensity)[2]/3 ))
dim(means_1)
for (i in 1:dim(intensity)[1])
{
  for (j in 1:(dim(intensity)[2]/3))
  {
    means_1[i,j]<-mean(c(intensity[i,3*j-2],intensity[i,3*j-1],intensity[i,3*j]),na.rm=T)
  }
}
View(means_1)

##Method-2 use lapply function-- way more time efficient  -- https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/lapply
colnum <- 1:ncol(intensity)
ind <- as.data.frame(matrix(colnum,byrow = F,nrow = 3))
means <- as.data.frame(lapply(ind, function(i) rowMeans(intensity[,i],na.rm = T)))
View(means)

##Generate new feature tables with mean intensities across triplicates
samp_id <- substr(as.vector(colnames(intensity)),1,nchar(colnames(intensity))-2)[c(T,rep(F,2))]  ##c(T,rep(F,2)) means 'TRUE FALSE FALSE'
samp_id
colnames(means) <- samp_id
View(means)
metabo <- cbind(mdat1[,1:9],means)
View(metabo)

##Log transformation
dim(metabo)
###check metabolomics feature intensity distribution -- check the fifth metabolic feature
hist(as.numeric(metabo[5,10:196]))
hist(as.numeric(log2(metabo[5,10:196])))
###Apply universal log2/log10 transformation to the feature table
metabo.log<-log2(metabo[,10:196])
View(metabo.log)
metlog<-cbind(metabo[,1:9],metabo.log)
View(metlog)
###Convert all NaN values to NA
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
metlog[is.nan(metlog)]<- NA
View(metlog)

##Upon this step, you have: 
## 1) filtered the feature table by removing low-presence and signals with bad data quality
## 2) replaced the colnames in the feature table with sample ID
## 3) replaced all missing values (0) with NA
## 4) calculated average across triplicates for each bio-sample for each metabolic feature
## 5) log-transformed the data 

##The feature table is clean and now ready for MWAS analysis
##save the dataset
write.csv('Met_mean_triplicate_log2_HILIC.csv')

#######Step 1.3 Merge the clean Metabolomic Feature Table with demographic table and pollutant table###########
###Import and clean the demographic Table
demo <- read_csv("demographic.csv")
View(demo)
dim(demo)
summary(demo)

##Create new variable to assign correct move in days for each stage
demo$movindays<-NA
for(i in 1:175)
{
  if(demo$stage[i]==1){demo$movindays[i]<-demo$days_1[i]}
  if(demo$stage[i]==2){demo$movindays[i]<-demo$days_2[i]}
  if(demo$stage[i]==3){demo$movindays[i]<-demo$days_3[i]}
  if(demo$stage[i]==4){demo$movindays[i]<-demo$days_4[i]}
}
summary(demo$movindays)

###Import the traffic related air pollution exposure assessment data
drweekMet <- read_csv("drweekMet.csv")
View(drweekMet)
##Create Pollution Variables in the demographic table with traffic related air pollution exposure assessment data
demo$BCout<-NA
demo$COout<-NA
demo$NOout<-NA
demo$NO2out<-NA
demo$NOxout<-NA
demo$PM25out<-NA

demo$BCin<-NA
demo$COin<-NA
demo$NOin<-NA
demo$NO2in<-NA
demo$NOxin<-NA
demo$PM25in<-NA



for(i in 1:175)
{
  if(demo$stage[i]==1)
  {
    if(demo$dorm_group[i]=="P")
    {
      demo$BCout[i]<-drweekMet$BCNDO[2]
      demo$BCin[i]<-drweekMet$BCNDI[2]
      demo$COout[i]<-drweekMet$CONDO[2]
      demo$COin[i]<-drweekMet$CONDI[2]
      demo$NOout[i]<-drweekMet$NONDO[2]
      demo$NOin[i]<-drweekMet$NONDI[2]
      demo$NO2out[i]<-drweekMet$NO2NDO[2]
      demo$NO2in[i]<-drweekMet$NO2NDI[2]
      demo$NOxout[i]<-drweekMet$NOxNDO[2]
      demo$NOxin[i]<-drweekMet$NOxNDI[2]
      demo$PM25out[i]<-drweekMet$PM25NDO[2]
      demo$PM25in[i]<-drweekMet$PM25NDI[2]
    }
    else if(demo$dorm_group[i]=="W")
    {
      demo$BCout[i]<-drweekMet$BCFDO[2]
      demo$BCin[i]<-drweekMet$BCFDI[2]
      demo$COout[i]<-drweekMet$COFDO[2]
      demo$COin[i]<-drweekMet$COFDI[2]
      demo$NOout[i]<-drweekMet$NOFDO[2]
      demo$NOin[i]<-drweekMet$NOFDI[2]
      demo$NO2out[i]<-drweekMet$NO2FDO[2]
      demo$NO2in[i]<-drweekMet$NO2FDI[2]
      demo$NOxout[i]<-drweekMet$NOxFDO[2]
      demo$NOxin[i]<-drweekMet$NOxFDI[2]
      demo$PM25out[i]<-drweekMet$PM25FDO[2]
      demo$PM25in[i]<-drweekMet$PM25FDI[2]
    }
    
  }
  
  else if(demo$stage[i]==2)
  {
    if(demo$dorm_group[i]=="P")
    {
      demo$BCout[i]<-drweekMet$BCNDO[7]
      demo$BCin[i]<-drweekMet$BCNDI[7]
      demo$COout[i]<-drweekMet$CONDO[7]
      demo$COin[i]<-drweekMet$CONDI[7]
      demo$NOout[i]<-drweekMet$NONDO[7]
      demo$NOin[i]<-drweekMet$NONDI[7]
      demo$NO2out[i]<-drweekMet$NO2NDO[7]
      demo$NO2in[i]<-drweekMet$NO2NDI[7]
      demo$NOxout[i]<-drweekMet$NOxNDO[7]
      demo$NOxin[i]<-drweekMet$NOxNDI[7]
      demo$PM25out[i]<-drweekMet$PM25NDO[7]
      demo$PM25in[i]<-drweekMet$PM25NDI[7]
    }
    else if(demo$dorm_group[i]=="W")
    {
      demo$BCout[i]<-drweekMet$BCFDO[7]
      demo$BCin[i]<-drweekMet$BCFDI[7]
      demo$COout[i]<-drweekMet$COFDO[7]
      demo$COin[i]<-drweekMet$COFDI[7]
      demo$NOout[i]<-drweekMet$NOFDO[7]
      demo$NOin[i]<-drweekMet$NOFDI[7]
      demo$NO2out[i]<-drweekMet$NO2FDO[7]
      demo$NO2in[i]<-drweekMet$NO2FDI[7]
      demo$NOxout[i]<-drweekMet$NOxFDO[7]
      demo$NOxin[i]<-drweekMet$NOxFDI[7]
      demo$PM25out[i]<-drweekMet$PM25FDO[7]
      demo$PM25in[i]<-drweekMet$PM25FDI[7]
    }
    
  } 
  
  else if(demo$stage[i]==3)
  {
    if(demo$dorm_group[i]=="P")
    {
      demo$BCout[i]<-drweekMet$BCNDO[10]
      demo$BCin[i]<-drweekMet$BCNDI[10]
      demo$COout[i]<-drweekMet$CONDO[10]
      demo$COin[i]<-drweekMet$CONDI[10]
      demo$NOout[i]<-drweekMet$NONDO[10]
      demo$NOin[i]<-drweekMet$NONDI[10]
      demo$NO2out[i]<-drweekMet$NO2NDO[10]
      demo$NO2in[i]<-drweekMet$NO2NDI[10]
      demo$NOxout[i]<-drweekMet$NOxNDO[10]
      demo$NOxin[i]<-drweekMet$NOxNDI[10]
      demo$PM25out[i]<-drweekMet$PM25NDO[10]
      demo$PM25in[i]<-drweekMet$PM25NDI[10]
    }
    else if(demo$dorm_group[i]=="W")
    {
      demo$BCout[i]<-drweekMet$BCFDO[10]
      demo$BCin[i]<-drweekMet$BCFDI[10]
      demo$COout[i]<-drweekMet$COFDO[10]
      demo$COin[i]<-drweekMet$COFDI[10]
      demo$NOout[i]<-drweekMet$NOFDO[10]
      demo$NOin[i]<-drweekMet$NOFDI[10]
      demo$NO2out[i]<-drweekMet$NO2FDO[10]
      demo$NO2in[i]<-drweekMet$NO2FDI[10]
      demo$NOxout[i]<-drweekMet$NOxFDO[10]
      demo$NOxin[i]<-drweekMet$NOxFDI[10]
      demo$PM25out[i]<-drweekMet$PM25FDO[10]
      demo$PM25in[i]<-drweekMet$PM25FDI[10]
    }
    
  } 
  else if(demo$stage[i]==4)
  {
    if(demo$dorm_group[i]=="P")
    {
      demo$BCout[i]<-drweekMet$BCNDO[13]
      demo$BCin[i]<-drweekMet$BCNDI[13]
      demo$COout[i]<-drweekMet$CONDO[13]
      demo$COin[i]<-drweekMet$CONDI[13]
      demo$NOout[i]<-drweekMet$NONDO[13]
      demo$NOin[i]<-drweekMet$NONDI[13]
      demo$NO2out[i]<-drweekMet$NO2NDO[13]
      demo$NO2in[i]<-drweekMet$NO2NDI[13]
      demo$NOxout[i]<-drweekMet$NOxNDO[13]
      demo$NOxin[i]<-drweekMet$NOxNDI[13]
      demo$PM25out[i]<-drweekMet$PM25NDO[13]
      demo$PM25in[i]<-drweekMet$PM25NDI[13]
    }
    else if(demo$dorm_group[i]=="W")
    {
      demo$BCout[i]<-drweekMet$BCFDO[13]
      demo$BCin[i]<-drweekMet$BCFDI[13]
      demo$COout[i]<-drweekMet$COFDO[13]
      demo$COin[i]<-drweekMet$COFDI[13]
      demo$NOout[i]<-drweekMet$NOFDO[13]
      demo$NOin[i]<-drweekMet$NOFDI[13]
      demo$NO2out[i]<-drweekMet$NO2FDO[13]
      demo$NO2in[i]<-drweekMet$NO2FDI[13]
      demo$NOxout[i]<-drweekMet$NOxFDO[13]
      demo$NOxin[i]<-drweekMet$NOxFDI[13]
      demo$PM25out[i]<-drweekMet$PM25FDO[13]
      demo$PM25in[i]<-drweekMet$PM25FDI[13]
    }
    
  } 
}

summary(demo)
summary(demo$COin)

##Merge the demographic table with the metabolomcis feature table
View(demo)
View(metlog)
###To merge, need to transform the metabolomic feature table
metdat<-as.data.frame(t(metlog[,10:196]))
metdat$sampleid<-rownames(metdat)

finaldt<-merge(demo,metdat,by="sampleid",all.x = T)

###Now, your feature table has been linked with the demographic/pollutant data. The final analytic database has been assembled 
###You are now ready for the MWAS data analysis

################## Step 2 Metabolomics-wide Association Study ###########################
###2.1 MWAS in Cross-sectional Study###
###create a cross-sectional dataset
finaldt_s1<-subset(finaldt,finaldt$stage==1)
summary(finaldt_s1[,1:40])
###use generalized linear regression/or linear regression as the MWAS model --glm()/lm
###https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/glm
##Take the 5th meatabolic feature for example
View(finaldt_s1)
hist(finaldt_s1$V5)
test_model1<-glm(finaldt_s1$V5~finaldt_s1$COout,family=gaussian())
summary(test_model1)
test_model2<-lm(finaldt_s1$V5~finaldt_s1$COout)
summary(test_model2)
##Need to control for other covariates
test_model3<-glm(finaldt_s1$V5~finaldt_s1$COout+finaldt_s1$age+finaldt_s1$gender+finaldt_s1$race+finaldt_s1$bmi
                 +finaldt_s1$college_yr+finaldt_s1$movindays,family=gaussian())
summary(test_model3)

###Now, how shall we do this for all these 6,100 features?
###Use for loop
###Need to create a data frame to store the model output
COoutHILIC<-data.frame(matrix(nrow = 6100, ncol = 4))
colnames(COoutHILIC)<-c("b","std","t","p")

for (i in 1:6100) 
{
  
      glmfit<-glm(finaldt_s1[,i+31]~finaldt_s1$COout+finaldt_s1$age+finaldt_s1$gender+finaldt_s1$race+finaldt_s1$bmi
                  +finaldt_s1$college_yr+finaldt_s1$movindays,family=gaussian())
      COoutHILIC[i,1]<-summary(glmfit)$coefficients[2,1]
      COoutHILIC[i,2]<-summary(glmfit)$coefficients[2,2]
      COoutHILIC[i,3]<-summary(glmfit)$coefficients[2,3]
      COoutHILIC[i,4]<-summary(glmfit)$coefficients[2,4]
      
      if(i%%10 ==0) {print(i)}
}

View(COoutHILIC)

###What happened? --check error
summary(finaldt_s1$V62)

###How to solve this?
###Use try catch function
for (i in 1:6100) 
{
  tryCatch(
    {
  glmfit<-glm(finaldt_s1[,i+31]~finaldt_s1$COout+finaldt_s1$age+finaldt_s1$gender+finaldt_s1$race+finaldt_s1$bmi
              +finaldt_s1$college_yr+finaldt_s1$movindays,family=gaussian())
  COoutHILIC[i,1]<-summary(glmfit)$coefficients[2,1]
  COoutHILIC[i,2]<-summary(glmfit)$coefficients[2,2]
  COoutHILIC[i,3]<-summary(glmfit)$coefficients[2,3]
  COoutHILIC[i,4]<-summary(glmfit)$coefficients[2,4]
  
  if(i%%100 ==0) {print(i)}
  }, error=function(e){})
}
View(COoutHILIC)
##Merge the MWAS output file with the feature table
COout_hp_s1<-cbind(metlog[,1:9],COoutHILIC)
View(COout_hp_s1)
##Address multiple comparison testing-- false positive discovery rate correction
hist(COout_hp_s1$p)
sum(COout_hp_s1$p<0.05,na.rm = T)
COout_hp_s1$pfdr<-p.adjust(COout_hp_s1$p,method="BH")
hist(COout_hp_s1$pfdr)
sum(COout_hp_s1$pfdr<0.05,na.rm = T)

##Save MWAS Result
write.csv(COout_hp_s1,paste0(dir_out,'HILIC_CO_out_s1.csv'))

###2.2 MWAS in Longitudinal Study with repeated measurement###
###use linear mixed effect model as the MWAS model --lmer()
###https://www.rdocumentation.org/packages/lme4/versions/1.1-26/topics/lmer
###Introduction on Mixed effect model https://stats.idre.ucla.edu/other/mult-pkg/introduction-to-linear-mixed-models/
library(lme4)
library(lmerTest)
##Take the 5th meatabolic feature for example
View(finaldt)
hist(finaldt$V5)
test_model1<-lmer(finaldt$V5~COout+(1|drid),data=finaldt)
summary(test_model1)
##Need to control for other covariates
test_model2<-lmer(finaldt$V5~COout+age+gender+race+bmi+college_yr+movindays+(1|drid),data=finaldt)
summary(test_model2)

###Now, how shall we do this for all these 6,100 features?
###Use for loop
###Need to create a data frame to store the model output
COoutHILICmix<-data.frame(matrix(nrow = 6100, ncol = 4))
colnames(COoutHILICmix)<-c("b","std","t","p")

###Use try catch function
for (i in 1:6100) 
{
  tryCatch(
    {
      lmerfit<-lmer(finaldt[,i+31]~COout+age+gender+race+bmi+college_yr+movindays+(1|drid),data=finaldt)
      COoutHILICmix[i,1]<-summary(lmerfit)$coefficients[2,1]
      COoutHILICmix[i,2]<-summary(lmerfit)$coefficients[2,2]
      COoutHILICmix[i,3]<-summary(lmerfit)$coefficients[2,4]
      COoutHILICmix[i,4]<-summary(lmerfit)$coefficients[2,5]
      
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
View(COoutHILICmix)
##Merge the MWAS output file with the feature table
COout_hp_mix<-cbind(metlog[,1:9],COoutHILICmix)
View(COout_hp_mix)
##Address multiple comparison testing-- false positive discovery rate correction
hist(COout_hp_mix$p)
sum(COout_hp_mix$p<0.05,na.rm = T)
COout_hp_mix$pfdr<-p.adjust(COout_hp_mix$p,method="BH")
hist(COout_hp_mix$pfdr)
sum(COout_hp_mix$pfdr<0.05,na.rm = T)

###Save Dataset
write.csv(COout_hp_mix,paste0(dir_out,'HILIC_CO_out_mix.csv'))


###Step 3 MWAS Study Interpretation ###
###Many steps are involved in this stage, including MWAS model results visualization, Pathway Enrichment Analysis, Chemical Annotation and Confirmation ###
###Numerous statistical and bioinformatics tools are involved in this process, including, but not limited to R, Python, Matlab, CytoScape, ect ###
###Below is an example on generating figures to visualize MWAS model results
###Visualize the MWAS model results

###Manhatthan plots
par(family = "sans",mar = c(4,2.5,2,1))
## retention time
xvec = COout_hp_s1$time
yvec = log10(COout_hp_s1$pfdr)*(-1)
coefficients = COout_hp_s1$b
ythresh = -log10(0.05)
colorvec = c("red","black","blue","darkgreen")
pchvec = c(21,24)
xincrement = 50
yincrement = 1
background.points.col="gray50"

plot(xvec,yvec,xlab="Retention Time (s)",ylab="-log10(pfdr)",xaxt="n",yaxt="n",cex=0.4,cex.main=1.5,xlim = c(0,300),main=paste0("HILIC ","CO-outdoor ","Stage 1"))
axis(1, at=seq(0,300, by=xincrement) , las=2,cex.axis = 1.5)
axis(2, at=seq(0 , (max(yvec,na.rm=TRUE)+2), by=yincrement) , las=2,cex.axis = 1.5)
points(xvec,yvec,col=background.points.col,bg=background.points.col,cex=0.4,pch=21)

goodip <- which(yvec>=ythresh)

posefficients <- which(coefficients >0)
negefficients <- which(coefficients <0)

for(i in intersect(goodip,posefficients)){
  points(xvec[i],yvec[i],col=colorvec[1],cex=1,pch=pchvec[1],bg=colorvec[1])
}

for(i in intersect(goodip,negefficients)){
  points(xvec[i],yvec[i],col=colorvec[3],cex=0.8,pch=pchvec[2],bg=colorvec[3])
}

if(length(goodip) > 0){
  abline(h = ythresh,col = "firebrick",lty=2,lwd=0.8)
  text(500, ythresh - 0.4, "p-value < 0.05", col = "black",cex = 1.5)
}

###Volcano plots
par(family = "sans",mar = c(4,2.5,2,1))
## retention time
xvec = COout_hp_s1$b
yvec = log10(COout_hp_s1$pfdr)*(-1)
coefficients = COout_hp_s1$b
ythresh = -log10(0.05)
colorvec = c("red","black","blue","darkgreen")
pchvec = c(21,24)
xincrement = 0.005
yincrement = 1
background.points.col="gray50"


plot(xvec,yvec,xlab="Beta Coefficients",ylab="-log10(pfdr)",xaxt="n",yaxt="n",cex=0.4,cex.main=1.5,xlim = c(-0.05,0.05),main=paste0("HILIC ","CO-outdoor ","Stage 1"))
axis(1, at=seq(-0.05,0.05, by=xincrement) , las=2,cex.axis = 1.5)
axis(2, at=seq(0 , (max(yvec,na.rm=TRUE)+2), by=yincrement) , las=2,cex.axis = 1.5)
points(xvec,yvec,col=background.points.col,bg=background.points.col,cex=0.4,pch=21)

goodip <- which(yvec>=ythresh)

posefficients <- which(coefficients >0)
negefficients <- which(coefficients <0)

for(i in intersect(goodip,posefficients)){
  points(xvec[i],yvec[i],col=colorvec[1],cex=1,pch=pchvec[1],bg=colorvec[1])
}

for(i in intersect(goodip,negefficients)){
  points(xvec[i],yvec[i],col=colorvec[3],cex=0.8,pch=pchvec[2],bg=colorvec[3])
}

if(length(goodip) > 0){
  abline(h = ythresh,col = "firebrick",lty=2,lwd=0.8)
  text(500, ythresh - 0.4, "p-value < 0.05", col = "black",cex = 1.5)
}

###We can also generate csv file for pathway enrichment analysis via mummichog (a bioinformatic tool)
###https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003123
View(COout_hp_s1)
COout_hp_s1_mum<-COout_hp_s1[,1:2]
COout_hp_s1_mum$p.value<-COout_hp_s1$pfdr
COout_hp_s1_mum$t.value<-COout_hp_s1$t
COout_hp_s1_mum<-na.omit(COout_hp_s1_mum)
View(COout_hp_s1_mum)

write.table(COout_hp_s1_mum,paste0(dir_out,"COout_hp_s1_mum.txt"),sep ="\t",row.names = FALSE, col.names = T)

###Then use this table and go to http://mummichog-2.appspot.com/ for pathway enrichment analysis 
###You can also use Python