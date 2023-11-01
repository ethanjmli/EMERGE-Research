###Read in the dataset
#Initiate necessary R packages
#install.packages("apLCMS")
#library(apLCMS)
install.packages("remotes")
remotes::install_github("yufree/apLCMS")
library(apLCMS)

#install.packages("xMSanalyzer")
#library(xMSanalyzer)

library(dplyr)
library(readr)
library(readxl)
setwd("/Users/ethan_j_li/Documents/Emory/Research/2017")
##2020
#?read_delim
#Read a delimited file into the object, in this case HILIC column results into object"hilic_results"
#Question: what is "\t" doing?
#h20 <- read_delim("/Users/tanyouran/Desktop/Liang/Thesis Materials/DATA_2 update/ComBat_mzcalibrated_untargeted_mediansummarized_featuretable_hilic.txt", 
#   "\t", escape_double = FALSE, trim_ws = TRUE)
#My version
hilic_results <- read_delim("/Users/ethan_j_li/Documents/Emory/Research/2017/ComBat_mzcalibrated_untargeted_averaged_featuretable_hilic.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
#Read a delimited file into the object, in this case the C18 column results into object "c18_results"
#c20 <- read_delim("/Users/tanyouran/Desktop/Liang/Thesis Materials/DATA_2 update/ComBat_mzcalibrated_untargeted_mediansummarized_featuretable_c18.txt", 
#                 "\t", escape_double = FALSE, trim_ws = TRUE)
c18_results <- read_delim("/Users/ethan_j_li/Documents/Emory/Research/2017/ComBat_mzcalibrated_untargeted_averaged_featuretable_c18.txt", 
                          "\t", escape_double = FALSE, trim_ws = TRUE)

###Link the demographic to the metabolomics
###Study ID Sample ID files
#Read delimited files into more objects, in this case the HILIC and C18 sample ids 
#into objects hilic_ids and c18_ids

#HPid20 <- read.delim("/Users/tanyouran/Desktop/Liang/Thesis Materials/DATA_2 update/R24_mapping_hilicpos.txt")
hilic_ids <- read.delim("/Users/ethan_j_li/Documents/Emory/Research/2017/hilicpos_sample_id_mapfile_combined.txt")

#CPid20 <- read.delim("/Users/tanyouran/Desktop/Liang/Thesis Materials/DATA_2 update/R24_mapping_c18neg.txt")
c18_ids <- read.delim("/Users/ethan_j_li/Documents/Emory/Research/2017/c18neg_sample_id_mapfile_combined.txt")

###Study ID Demongraphic 

#Read demographic info, e.g. dates of visit, birth outcome, by StudyID into object "demographicbystudyid"
#Sid<-read_excel("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Newborn/Demo/MPTB_merged2021_12032021_corrected.xlsm",2)
#Sid<-read.csv("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Newborn/Demo/MPTB_merged2021_05062021.csv")
demographicsbystudyid<-read.csv("/Users/ethan_j_li/Documents/Emory/Research/MPTB_AA Cohort Data_N547_addZscore.csv")

###2020
##HP
#Assign HPid20 to HPid for manipulation
#Assign hilic_ids to hilicID for manipulation
hilicID<-hilic_ids

#Create variable seq in hilicID, value of seq = row number
hilicID$seq<-(1:nrow(hilicID))

#Create variable subjectid in hilicID, values of which = substrings of variable Sample.ID 
#consisting of characters 1 through 5. Creates common variable w/demographicbystudyid for merging
#Ex: nist1_1 --> nist1, GOO73-2-Red2.1_003 --> G0073
#OG = Sample.ID, new = SampleID
hilicID$subjectid<-substr(hilicID$SampleID,1,5)

#Create variable stage in hilicID, values of which = character 7 of Sample.ID
#Corresponds to the stage of the visit of each observation, i.e. did this obs occur 
#in stage 1 (early) or stage 2 (late) of pregnancy
#Ex: G0232-1-Red2.1_001 --> 1
#OG = Sample.ID, new = SampleID
hilicID$stage<-substr(hilicID$SampleID,7,7)

#Order object demographicsbystudyid by variable subjectid
demographicsbystudyid<-demographicsbystudyid[order(demographicsbystudyid$subjectid),]

###Merge into a complete SID PID Demongraphic Table
#Merge HPid, Sid by variable subjectid into object "total"
#Merge hilicID and demographicsbystudyid by variable subjectid into object "total"
#Combine demographic information and HILIC column info by subjectid
total <- merge(hilicID,demographicsbystudyid,by=c("subjectid"))
# totals1<-subset(total,total$stage==1)
# totals2<-subset(total,total$stage==2)

###Transform the Metabolomic Dataset
#Convert transpose of h20, t(h20), into a data frame with as.data.frame(t(h20)) and assign to object "HPt"
#Convert transpose of hilic_results, t(hilic_results), into a data frame with 
#as.data.frame(t(hilic_results)) and assign to object "HRt", prepare for merging 
#Has metabolic features as columns and variables as rows
HRt<-as.data.frame(t(hilic_results))

#Retrieve the rownames of HPt, checking for successful transposition
#Rownames of transposed (HPt) should be column names of pre-transposed (h20)
#Rownames of (HRt) should be column names of pre-transposed (hilic_results)
rownames(HRt)

#Assign rownames of HPt into list HPt$File.Name
#Create variable File.Name in HRt, values of which = rownames of HRt = column 
#names of hilic_results = diff variables observed in HILIC
HRt$FileName<-rownames(HRt) 

#Check for successful assignment by HRt$File.Name<-rownames(HRt)
HRt$FileName 
table(HRt$FileName %in% rownames(HRt)) #TRUE 
                                        # 674 

##Merge total, the merge of hilicID and demographicsbystudyid, with HRt, the transposition of hilic_results, 
#by the common variable File.name, into object "Metotal"
#Metotal = a table containing, for every sample/file name, basic HILIC column info,
#demographics, and HILIC column results
#Question: why do we have many more observations in total (HILIC info, demographics) than HRt(column results)
Metotal<-merge(total,HRt, by = "FileName")
#Metotal: 190 variables in total + 16482 variables in HRt - 1 overlap (File.name) = 16671 variables in Metotal
#Metabolic features start at Metotal col# 191 = metabolic feature 1 (colname V1)
# Metotals<-Metotal[,201:300]

#Analogous to PROC FREQ, breakdown the frequency of each value of the variable "stage" in Metotal. 
#Total = 650, = # of obs in Metotal = # of samples run through HILIC (check, question)
table(Metotal$stage)    #   1   2   3 
                        # 316 256   1

###Metotal is the Final Dataset, contains basic column info, demographics, HILIC column results
###Create Gestational Age Variable

#Initialize variable gast in Metotal, initial value = NA (not available) for each observation
Metotal$gast<-NA
for (i in 1:nrow(Metotal))
{
  if (Metotal$stage[i]=='1')
  {
    Metotal$gast[i]<-Metotal$gasamp1_wk.10[i]
  }
  else if (Metotal$stage[i]=='2')
  {
    Metotal$gast[i]<-Metotal$gasamp2_wk.10[i]
  }
}
summary(Metotal$gasamp1_wk.10)
summary(Metotal$gasamp2_wk.10)
summary(Metotal$gast)


###Create nulliparious
summary(Metotal$parity)
#Initialize variable nullp, initial value = NA for all observations
Metotal$nullp<-NA
for (i in 1:nrow(Metotal))
{
  if (Metotal$parity[i]==0)
  {
    Metotal$nullp[i]<-0
  }
  else if (Metotal$parity[i]>0)
  {
    Metotal$nullp[i]<-1
  }
}
table(Metotal$nullp) #0   1 
                    #259 314    

###Log-transformation
#Get names of columns 284-289 in Metotal
#colnames(Metotal[,284:289])    ##284:289 #V94, V95, V96, V97, V98, V99
colnames(Metotal[,185:195])     ##185:195 metabolic feature data begins at 191
#Get names of columns 13897:13902 in Metotal
#colnames(Metotal[,13897:13902])  ##metabolic feature data ends at 16671; 16581 metabolic features (V1 - V16481)
colnames(Metotal[,16663:16673])

#Create metotal containing columns 1-284 of Metotal, gast, and nullp
#AKA only demographic info, gestational age, nulliparity
metotal<-Metotal[,1:190]
metotal$gast<-Metotal$gast
metotal$nullp<-Metotal$nullp

# metotal[,286:13901]<-log2(Metotal[,283:13898])   get demographic table first 
# metotal[metotal == -Inf] <- NA


####################################
table(metotal$stage) # stage1: 330 / stage2: 270 / -:1 / P: 49 
#Question: What do "-" and "P" represent?
# subset to subjects attending stage1 and stage2
#Generate metotalmis, all observations in metotal with stage ==1 or stage ==2
#AKA all observations recorded in stage 1 or stage 2
metotalmis<-subset(metotal,stage==1| stage==2) #600 observations, exclude 1 - and 49 P

###Demographic
#Checking for # of observations in each birth category
table(metotalmis$birthoutcome,useNA="always")
table(metotalmis$Comp_SpPTBvsALLFullTerm,useNA="always")   # used in PTB model  NA 168 NewNA = 122
table(metotalmis$Comp_SpPTBvsHealthyFullTerm,useNA="always")  # used in sPTB model  NA 217 NewNA = 215

table(metotalmis$Comp_AllPTBvsHealthyFullTerm,useNA="always")  # NA 192 NewNA = 190
table(metotalmis$Comp_SpPTBSpETBvsHealthyFullTerm,useNA="always")  #NA 145 NewNA = 143


table(metotalmis[which(metotalmis$Comp_SpPTBvsHealthyFullTerm!=" "),]$Comp_SpPTBvsALLFullTerm,useNA="always")
# NA 33 NewNA = 0

table(metotalmis[which(metotalmis$Comp_SpPTBvsHealthyFullTerm!=" "),]$Comp_AllPTBvsHealthyFullTerm,useNA="always")
# NA 0

table(metotalmis[which(metotalmis$Comp_SpPTBvsHealthyFullTerm!=" "),]$Comp_SpPTBSpETBvsHealthyFullTerm,useNA="always")
# NA 0

table(metotalmis[which(metotalmis$Comp_SpPTBSpETBvsHealthyFullTerm!=" "),]$Comp_AllPTBvsHealthyFullTerm,useNA="always")
# NA 72

table(metotalmis[which(metotalmis$Comp_AllPTBvsHealthyFullTerm!=" "),]$Comp_SpPTBvsHealthyFullTerm,useNA="always")
# NA 25

table(metotalmis[which(metotalmis$Comp_SpPTBSpETBvsHealthyFullTerm!=" "),]$Comp_SpPTBvsHealthyFullTerm,useNA="always")
# NA 72

###############
# subset to exclude all PTB and sPTB non missing
# metotalnmis<-subset(metotalmis,Comp_SpPTBvsALLFullTerm !=" " 
#                     | Comp_SpPTBSpETBvsHealthyFullTerm !=" " 
#                     | Comp_SpPTBvsHealthyFullTerm !=" "
#                     | Comp_AllPTBvsHealthyFullTerm !=" ")  # get 562 obs in total

###### We should ONLY include women with PTB ETB FTB and should not include spontaneous/elective abortions
table(metotalmis$birthoutcome, useNA = "always")
#create metotalnmis, subset of metotalmis only including observations with values of 
#birthoutcome == EarlyTerm, FullTerm, Preterm
#Essentially, metotalnmis includes observations who attended 1 or 2 visits, with 
#birth outcomes of Early Term, Full Term, or Preterm birth, excluding 
metotalnmis<-subset(metotalmis,birthoutcome=="EarlyTerm" | birthoutcome=="FullTerm" |
                      birthoutcome=="Preterm")    # 546 obs in total, excluding 1 Preterm-STILLBIRTH, presumably aborted


table(metotalnmis$stage) # 290 stage1 and 256 stage 2  
                         #metotalnmis for later modeling

#################
# need to figure out how many unique subjects in total 
table(table(metotalnmis$subjectid)) # subjectid E0586 repeated twice on first visit
#remove this repeated obs
#metotalnmis<-metotalnmis[!metotalnmis$Sample.ID=="E0586-1-Red2.1_007",] #Now 598 total obs after eliminating 1 repeated obs
#table(table(metotalnmis$subjectid))
#252 appear twice, 42 appear once

#Assign to rn1, the row names of the subset of metotalnmis whose subjectid only appears once 
#AKA the subset of metotalnmis who appear twice
rn1<-rownames(subset(table(metotalnmis$subjectid),table(metotalnmis$subjectid)==1))  
#Assign to rn2, the row names of the subset of metotalnmis with subject id does not only appear once
#AKA the subset of metotalnmis who appear once
rn2<-rownames(subset(table(metotalnmis$subjectid),table(metotalnmis$subjectid)!=1))

#Assign to metotalv1, the observations of metotalnmis whose subjectids are in rn1
metotalv1<-metotalnmis[which(metotalnmis$subjectid %in% rn1),]           ######single stage 
table(metotalv1$subjectid)    # 42 obs- 42 usbjects

#Assign to metotalv2, the observations of metotalnmis whose subjectids are in rn2
metotalv2<-metotalnmis[which(metotalnmis$subjectid %in% rn2),]     ###### both stages
table(metotalv2$subjectid)   # 504 obs = 252 subjects total, 2 visits each: 252+42 = 294 subjects
# 328 subjects visit1-270 visit2
#252 subjects each had 2 observations = 504 observations in metotalv2, + 42 subjects w/ only 1 observation in metotalv1 = 294 total subjects
metotalv2s<-subset(metotalv2,stage==1)
#Create a, table consisting only of unique subjects from metotalv1 and metotalv2. 331 total unique subjects
#metotalv2 contains all subjects who appear twice, first visit stage ==1, second visit stage ==2
#By only taking the observations whose stage == 1, we take only one instance of each, get unique subjects
a<-rbind(metotalv1,metotalv2s) 

#Get a list of variables in a
dput(names(a))
#Question: why only these variables?
#List of continuous variables in a
contVars<-c("age_enrollment","birthga","FirstPrenatalBMI")
#List of categorical variables in  a
catVars<-c("Education_4.level","income","PregInsurance_2.level",
           #my addition#
           "PresInsurance_3.level", 
           #my addition end#
           "nullp","priorabortion","Site_Prenatal","Depression_Preg", "MarriedCohab_Not","priorpreterm",
           "Sex","g_htn","iugr","g_dm","mode","birthoutcome","TobaccoUse_MR","AlcoholUse_MR")

#Create Table One for continuous variables from a listed in contVars
library(tableone)
tab0<-CreateTableOne(vars=contVars, data = a)
print(tab0, showAllLevels = TRUE, test = FALSE)

#Create Table One for categorical variables from a listed in catVars
tab2<-CreateCatTable(vars=catVars, data= a, includeNA = TRUE)
print(tab2, catDigits = 3,pDigits = 3,showAllLevels = TRUE, test = FALSE)

#What we did here is create "a", a table with only unique subjects, then extract the 
#baseline demographic data (at visit 1, all observations have stage ==1), categorical and 
#continuous, into Table 1

# tab2<-CreateCatTable(vars=catVars, data= a, includeNA = TRUE,strata="birthoutcome")
# print(tab2, catDigits = 3,pDigits = 3,showAllLevels = TRUE, test = FALSE)


########## For gestational weeks
table(metotalnmis$stage)
#Assign to "b", a subset of metotalnmis with values of stage == 1
#metotalnmis = subset of metotalmis only including observations with values of birthoutcome == EarlyTerm, FullTerm, Preterm
#metotalmis = all observations in metotal with stage ==1 or stage ==2
b<-subset(metotalnmis,stage==1) #3 unique subjects only attended a stage == 2 visit, only 328 obs in b 
#Only appeared once, but stage == 2 only, why b != a 

dput(names(b))
contVars<-c("gast")

#Create table of mean gestational age for all observations in metotalnmis whose stage ==1
tab0<-CreateTableOne(vars=contVars, data = b)
print(tab0, showAllLevels = TRUE, test = FALSE)


c<-subset(metotalnmis,stage==2)
dput(names(c))

#Create table of mean gestational age for all observations in metotalnmis whose stage == 2
tab0<-CreateTableOne(vars=contVars, data = c)
print(tab0, showAllLevels = TRUE, test = FALSE)

#############Look into outcomes of interest in our analysis
dput(names(a))
#Create new list of categorical variables of interest consisting of birth outcomes of interest
catVars<-c("Comp_SpPTBvsALLFullTerm","Comp_AllPTBvsHealthyFullTerm","Comp_InducedPTBvsHealthyFullTerm", "Comp_SpPTBvsHealthyFullTerm",
           "Comp_AllETBvsHealthyFullTerm","Comp_SpETBvsHealthyFullTerm", "Comp_InducedETBvsHealthyFullTerm","Comp_SpPTBSpETBvsHealthyFullTerm","g_htn","iugr","g_dm")

#Generate tables for categorical outcome variables of interest, stratified on one variable at a time
tab2<-CreateCatTable(vars=catVars, data= a, includeNA = TRUE)
print(tab2, catDigits = 3,pDigits = 3,showAllLevels = TRUE, test = FALSE)

tab2<-CreateCatTable(vars=catVars, data= a, includeNA = TRUE,strata = "Comp_SpPTBvsHealthyFullTerm")
print(tab2, catDigits = 3,pDigits = 3,showAllLevels = TRUE, test = FALSE)

tab2<-CreateCatTable(vars=catVars, data= a, includeNA = TRUE,strata = "Comp_SpETBvsHealthyFullTerm")
print(tab2, catDigits = 3,pDigits = 3,showAllLevels = TRUE, test = FALSE)

tab2<-CreateCatTable(vars=catVars, data= a, includeNA = TRUE,strata = "Comp_InducedPTBvsHealthyFullTerm")
print(tab2, catDigits = 3,pDigits = 3,showAllLevels = TRUE, test = FALSE)

tab2<-CreateCatTable(vars=catVars, data= a, includeNA = TRUE,strata = "Comp_InducedETBvsHealthyFullTerm")
print(tab2, catDigits = 3,pDigits = 3,showAllLevels = TRUE, test = FALSE)

######################
####create prenatal complication indicator variable---g_htn,iugr,g_dm,Depression_Preg
# metotalnmis$pre_complication<-ifelse(metotalnmis$g_htn!=0 | metotalnmis$g_dm==1 | 
#                                       metotalnmis$iugr==1 | metotalnmis$Depression_Preg==1,1,0)
# table(metotalnmis$pre_complication)

### create prenatal complication factor variable ---0,1,2 only include ghtn gdm
#Question: logic behind this code?
metotalnmis$pre_complication<-ifelse(metotalnmis$g_htn==1 & metotalnmis$g_dm==1,2, 
                                     ifelse(metotalnmis$g_htn==2 & metotalnmis$g_dm==1,2, 
                                            ifelse(metotalnmis$g_htn==0 & metotalnmis$g_dm==0,0,1)))
table(metotalnmis$pre_complication)

#########log transform the dataset--unify the row observations
#Recall metotal = demographic + column info + metabolic features
Metotalnmis<-subset(Metotal,stage==1| stage==2)
Metotalnmis<-subset(Metotalnmis,birthoutcome=="EarlyTerm" | birthoutcome=="FullTerm" |
                      birthoutcome=="Preterm")
#Metotalnmis<-Metotalnmis[!Metotalnmis$Sample.ID=="E0586-1-Red2.1_007",] Not needed, not a duplicate obs in 2017

#metotalnmis[,288:13903]<-log2(Metotalnmis[,285:13808])
#log transform intensities of metabolic features
metotalnmis[,194:16674]<-log2(Metotalnmis[,191:16671])
#Replace all values of -Infinity resulting from the log2 transform above with NA values
#AKA all values in Metotalnmis which were very small
metotalnmis[metotalnmis == -Inf] <- NA
#save(metotalnmis,file="metotalnmisH.RData")

setwd("/Users/ethan_j_li/Documents/Emory/Research/2017/HILIC 2017 PreComp Regressions")

####### Need to create combined group--PTB+ETB=PETB
metotalnmis$Comp_PTBETBvsHealthyFullTerm<-ifelse(metotalnmis$Comp_AllETBvsHealthyFullTerm==1 | metotalnmis$Comp_AllPTBvsHealthyFullTerm==1,1,
                                                 ifelse(metotalnmis$Comp_AllETBvsHealthyFullTerm==0,0,NA))

####### Need to create combined group--miPTB+miETB=miPETB
metotalnmis$Comp_InducedPTBETBvsHealthyFullTerm<-ifelse(metotalnmis$Comp_InducedPTBvsHealthyFullTerm==1 | metotalnmis$Comp_InducedETBvsHealthyFullTerm==1,1,
                                                        ifelse(metotalnmis$Comp_InducedPTBvsHealthyFullTerm==0,0,NA))

###Stage 1
metotals1<-subset(metotalnmis,metotalnmis$stage==1)
#Why do we only analyze stage 1 visits?

## sPTB vs All full term Done
#Old = h20, new = hilic_results
h_SpPTBvsALLFullTerm_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_SpPTBvsALLFullTerm_s1)<-c("Est","Std","t-value","P")

#Question: regress intensity of metabolic features recorded at stage 1/visit 1, 
#in metotals1, against Spontaneous Preterm birth, Sex, age, nulliparity, education, 
#first prenatal bmi, as recorded in metotals1
#Iterate through each metabolic feature present in hilic_results (Ethan Li: 13616 present)
for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      #Regress every metabolic feature against birth outcome, sex, enrollment age, nulliparity, education, first prenatal bmi recorded at stage 1
      lmfit<-lm( metotals1[,i+193] ~Comp_SpPTBvsALLFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      #Write regression coefficients
      h_SpPTBvsALLFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
#Adjust p-values w/Benjamin-Hochberg correction for multiple hypothesis testing
h_SpPTBvsALLFullTerm_s1$pfdr<-p.adjust(h_SpPTBvsALLFullTerm_s1$P,method="BH")
#Output regression results to csv file
write.csv(h_SpPTBvsALLFullTerm_s1,"h_SpPTBvsALLFullTerm_s1.csv")


## sPTB sETB vs Healthy full term Done
h_SpPTBSpETBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_SpPTBSpETBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+193] ~Comp_SpPTBSpETBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      h_SpPTBSpETBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_SpPTBSpETBvsHealthyFullTerm_s1$pfdr<-p.adjust(h_SpPTBSpETBvsHealthyFullTerm_s1$P,method="BH")
write.csv(h_SpPTBSpETBvsHealthyFullTerm_s1,"h_SpPTBSpETBvsHealthyFullTerm_s1.csv")

## PTB ETB vs Healthy full term Done
h_PTBETBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_PTBETBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+193] ~Comp_PTBETBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      h_PTBETBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_PTBETBvsHealthyFullTerm_s1$pfdr<-p.adjust(h_PTBETBvsHealthyFullTerm_s1$P,method="BH")
#h_PTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_PTBETBvsHealthyFullTerm_s1)
write.csv(h_PTBETBvsHealthyFullTerm_s1,"h_PTBETBvsHealthyFullTerm_s1.csv")


## sPTB  vs Healthy full term Done
h_SpPTBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_SpPTBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+193] ~Comp_SpPTBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      h_SpPTBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_SpPTBvsHealthyFullTerm_s1$pfdr<-p.adjust(h_SpPTBvsHealthyFullTerm_s1$P,method="BH")
write.csv(h_SpPTBvsHealthyFullTerm_s1,"h_SpPTBvsHealthyFullTerm_s1.csv")


## AllPTB  vs Healthy full term Done
h_AllPTBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_AllPTBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+193] ~Comp_AllPTBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      h_AllPTBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_AllPTBvsHealthyFullTerm_s1$pfdr<-p.adjust(h_AllPTBvsHealthyFullTerm_s1$P,method="BH")
write.csv(h_AllPTBvsHealthyFullTerm_s1,"h_AllPTBvsHealthyFullTerm_s1.csv")


## Induced PTB  vs Healthy full term Done
h_InducedPTBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_InducedPTBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+193] ~Comp_InducedPTBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      h_InducedPTBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_InducedPTBvsHealthyFullTerm_s1$pfdr<-p.adjust(h_InducedPTBvsHealthyFullTerm_s1$P,method="BH")
write.csv(h_InducedPTBvsHealthyFullTerm_s1,"h_InducedPTBvsHealthyFullTerm_s1.csv")


## All ETB  vs Healthy full term Done
h_AllETBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_AllETBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+193] ~Comp_AllETBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      h_AllETBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_AllETBvsHealthyFullTerm_s1$pfdr<-p.adjust(h_AllETBvsHealthyFullTerm_s1$P,method="BH")
write.csv(h_AllETBvsHealthyFullTerm_s1,"h_AllETBvsHealthyFullTerm_s1.csv")

## spETB  vs Healthy full term Done
h_SpETBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_SpETBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+193] ~Comp_SpETBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      h_SpETBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_SpETBvsHealthyFullTerm_s1$pfdr<-p.adjust(h_SpETBvsHealthyFullTerm_s1$P,method="BH")
write.csv(h_SpETBvsHealthyFullTerm_s1,"h_SpETBvsHealthyFullTerm_s1.csv")

## Induced ETB  vs Healthy full term Done
h_InducedETBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_InducedETBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+193] ~Comp_InducedETBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      h_InducedETBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_InducedETBvsHealthyFullTerm_s1$pfdr<-p.adjust(h_InducedETBvsHealthyFullTerm_s1$P,method="BH")
write.csv(h_InducedETBvsHealthyFullTerm_s1,"h_InducedETBvsHealthyFullTerm_s1.csv")

## Induced PTB ETB  vs Healthy full term Done
h_InducedPTBETBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_InducedPTBETBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+193] ~Comp_InducedPTBETBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      h_InducedPTBETBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_InducedPTBETBvsHealthyFullTerm_s1$pfdr<-p.adjust(h_InducedPTBETBvsHealthyFullTerm_s1$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(h_InducedPTBETBvsHealthyFullTerm_s1,"h_InducedPTBETBvsHealthyFullTerm_s1.csv")

#g_htn S1 HILIC regression Done
h_g_htn_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_g_htn_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+193] ~g_htn+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      h_g_htn_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_g_htn_s1$pfdr<-p.adjust(h_g_htn_s1$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(h_g_htn_s1,"h_g_htn_s1_2017.csv")

#prenatal complications S1 HILIC Done
h_pre_complication_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_pre_complication_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+193] ~as.factor(pre_complication)+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      h_pre_complication_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_pre_complication_s1$pfdr<-p.adjust(h_pre_complication_s1$P,method="BH")
write.csv(h_pre_complication_s1,"h_pre_complication_s1_2017.csv")

#g_dm S1 HILIC regression Done
h_g_dm_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_g_dm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+193] ~g_dm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      h_g_dm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_g_dm_s1$pfdr<-p.adjust(h_g_dm_s1$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(h_g_dm_s1,"h_g_dm_s1_2017.csv")

#PROM S1 HILIC regression Done
h_PROM_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_PROM_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+193] ~PROM+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      h_PROM_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_PROM_s1$pfdr<-p.adjust(h_PROM_s1$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(h_PROM_s1,"h_PROM_s1_2017.csv")

#IUGR S1 HILIC regression Done
h_IUGR_s1<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_IUGR_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+193] ~Comp_IUGRvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      h_IUGR_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_IUGR_s1$pfdr<-p.adjust(h_IUGR_s1$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(h_IUGR_s1,"h_IUGR_s1_2017.csv")

###Stage 2 HILIC 
metotals2<-subset(metotalnmis,metotalnmis$stage==2)

## sPTB vs All full term Done
h_SpPTBvsALLFullTerm_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_SpPTBvsALLFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~Comp_SpPTBvsALLFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_SpPTBvsALLFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_SpPTBvsALLFullTerm_s2$pfdr<-p.adjust(h_SpPTBvsALLFullTerm_s2$P,method="BH")
write.csv(h_SpPTBvsALLFullTerm_s2,"h_SpPTBvsALLFullTerm_s2.csv")


## sPTB sETB vs Healthy full term Done
h_SpPTBSpETBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_SpPTBSpETBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~Comp_SpPTBSpETBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_SpPTBSpETBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_SpPTBSpETBvsHealthyFullTerm_s2$pfdr<-p.adjust(h_SpPTBSpETBvsHealthyFullTerm_s2$P,method="BH")
write.csv(h_SpPTBSpETBvsHealthyFullTerm_s2,"h_SpPTBSpETBvsHealthyFullTerm_s2.csv")


## sPTB  vs Healthy full term Done
h_SpPTBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_SpPTBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~Comp_SpPTBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_SpPTBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_SpPTBvsHealthyFullTerm_s2$pfdr<-p.adjust(h_SpPTBvsHealthyFullTerm_s2$P,method="BH")
write.csv(h_SpPTBvsHealthyFullTerm_s2,"h_SpPTBvsHealthyFullTerm_s2.csv")


## AllPTB  vs Healthy full term Done
h_AllPTBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_AllPTBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~Comp_AllPTBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_AllPTBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_AllPTBvsHealthyFullTerm_s2$pfdr<-p.adjust(h_AllPTBvsHealthyFullTerm_s2$P,method="BH")
write.csv(h_AllPTBvsHealthyFullTerm_s2,"h_AllPTBvsHealthyFullTerm_s2.csv")

## Induced PTB  vs Healthy full term Done
h_InducedPTBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_InducedPTBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~Comp_InducedPTBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_InducedPTBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_InducedPTBvsHealthyFullTerm_s2$pfdr<-p.adjust(h_InducedPTBvsHealthyFullTerm_s2$P,method="BH")
write.csv(h_InducedPTBvsHealthyFullTerm_s2,"h_InducedPTBvsHealthyFullTerm_s2.csv")

## AllETB  vs Healthy full term Done
h_AllETBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_AllETBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~Comp_AllETBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_AllETBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_AllETBvsHealthyFullTerm_s2$pfdr<-p.adjust(h_AllETBvsHealthyFullTerm_s2$P,method="BH")
write.csv(h_AllETBvsHealthyFullTerm_s2,"h_AllETBvsHealthyFullTerm_s2.csv")

## spETB  vs Healthy full term Done
h_SpETBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_SpETBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~Comp_SpETBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_SpETBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_SpETBvsHealthyFullTerm_s2$pfdr<-p.adjust(h_SpETBvsHealthyFullTerm_s2$P,method="BH")
write.csv(h_SpETBvsHealthyFullTerm_s2,"h_SpETBvsHealthyFullTerm_s2.csv")

## Induced ETB  vs Healthy full term Done
h_InducedETBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_InducedETBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~Comp_InducedETBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_InducedETBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_InducedETBvsHealthyFullTerm_s2$pfdr<-p.adjust(h_InducedETBvsHealthyFullTerm_s2$P,method="BH")
write.csv(h_InducedETBvsHealthyFullTerm_s2,"h_InducedETBvsHealthyFullTerm_s2.csv")

## PTB ETB vs Healthy full term Done
h_PTBETBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_PTBETBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~Comp_PTBETBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_PTBETBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_PTBETBvsHealthyFullTerm_s2$pfdr<-p.adjust(h_PTBETBvsHealthyFullTerm_s2$P,method="BH")
#h_PTBETBvsHealthyFullTerm_s2<-cbind(h20[,c(1:2)],h_PTBETBvsHealthyFullTerm_s2)
write.csv(h_PTBETBvsHealthyFullTerm_s2,"h_PTBETBvsHealthyFullTerm_s2.csv")

## Induced PTB ETB  vs Healthy full term Done
h_InducedPTBETBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_InducedPTBETBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~Comp_InducedPTBETBvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_InducedPTBETBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_InducedPTBETBvsHealthyFullTerm_s2$pfdr<-p.adjust(h_InducedPTBETBvsHealthyFullTerm_s2$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s2<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s2)
write.csv(h_InducedPTBETBvsHealthyFullTerm_s2,"h_InducedPTBETBvsHealthyFullTerm_s2.csv")

#pre_complication stage 2
h_pre_complication_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_pre_complication_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~as.factor(pre_complication)+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_pre_complication_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_pre_complication_s2$pfdr<-p.adjust(h_pre_complication_s2$P,method="BH")
write.csv(h_pre_complication_s2,"h_pre_complication_s2_2017.csv")

#g_htn S2 HILIC regression Done
h_g_htn_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_g_htn_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~g_htn+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_g_htn_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_g_htn_s2$pfdr<-p.adjust(h_g_htn_s2$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(h_g_htn_s2,"h_g_htn_s2_2017.csv")

#g_dm S2 HILIC regression Done
h_g_dm_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_g_dm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~g_dm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_g_dm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_g_dm_s2$pfdr<-p.adjust(h_g_dm_s2$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(h_g_dm_s2,"h_g_dm_s2_2017.csv")

#PROM S2 HILIC regression Done
h_PROM_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_PROM_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~PROM+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_PROM_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_PROM_s2$pfdr<-p.adjust(h_PROM_s2$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(h_PROM_s2,"h_PROM_s2_2017.csv")

#IUGR S2 HILIC regression Current
h_IUGR_s2<-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(h_IUGR_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(hilic_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+193] ~Comp_IUGRvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      h_IUGR_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
h_IUGR_s2$pfdr<-p.adjust(h_IUGR_s2$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(h_IUGR_s2,"h_IUGR_s2_2017.csv")
#######################################
##########CP
#CPid<-CPid20
c18ID <- c18_ids
#CPid$seq<-(1:nrow(CPid))
c18ID$seq<-(1:nrow(c18ID))
#CPid$subjectid<-substr(CPid$Sample.ID,1,5)
c18ID$subjectid<-substr(c18ID$SampleID,1,5)
#CPid$stage<-substr(CPid$Sample.ID,7,7)
c18ID$stage <- substr(c18ID$SampleID,7,7)

#Sid<-Sid[order(Sid$subjectid),]
demographicsbystudyid<-demographicsbystudyid[order(demographicsbystudyid$subjectid),]

###Merge into a complete SID PID Demongraphic Table 
#Demographics + C18 column basic info
total2 <- merge(c18ID,demographicsbystudyid,by=c("subjectid"))

###Transform the Metabolomic Dataset
#CPt<-as.data.frame(t(c20))
CRt <-as.data.frame(t(c18_results))
rownames(CRt)
CRt$FileName<-rownames(CRt)
CRt$FileName
###Merge Two together
#Demographic + C18 column basic info + C18 metabolic features
Metotal2<-merge(total2,CRt, by="FileName")
table(Metotal2$stage)

###Create Gestational Age Variable
Metotal2$gast<-NA
for (i in 1:nrow(Metotal2))
{
  if (Metotal2$stage[i]=='1')
  {
    Metotal2$gast[i]<-Metotal2$gasamp1_wk.10[i]
  }
  else if (Metotal2$stage[i]=='2')
  {
    Metotal2$gast[i]<-Metotal2$gasamp2_wk.10[i]
  }
}
summary(Metotal2$gasamp1_wk.10)
summary(Metotal2$gasamp2_wk.10)
summary(Metotal2$gast)


###Create nulliparious
summary(Metotal2$parity)
Metotal2$nullp<-NA
for (i in 1:nrow(Metotal2))
{
  if (Metotal2$parity[i]==0)
  {
    Metotal2$nullp[i]<-0
  }
  else if (Metotal2$parity[i]>0)
  {
    Metotal2$nullp[i]<-1
  }
}
table(Metotal2$nullp)

metotal2<-Metotal2[,1:190] #Demographics+ C18 basic column info + gestational age + nulliparity
metotal2$gast<-Metotal2$gast
metotal2$nullp<-Metotal2$nullp
#table(metotal2$stage) 
#Demographics+ C18 basic column info + gestational age + nulliparity only inc stage 1 + stage 2
metotalmis2<-subset(metotal2,stage==1| stage==2)
# metotalnmis<-subset(metotalmis,Comp_SpPTBvsALLFullTerm !=" " 
#                     | Comp_SpPTBSpETBvsHealthyFullTerm !=" " 
#                     | Comp_SpPTBvsHealthyFullTerm !=" "
#                     | Comp_AllPTBvsHealthyFullTerm !=" ") 
#table(metotalmis2$birthoutcome)
#Demographics+ C18 basic column info + gestational age + nulliparity only inc stage 1 + stage 2 w/o abortion
metotalnmis2<-subset(metotalmis2,birthoutcome=="EarlyTerm" | birthoutcome=="FullTerm" |
                       birthoutcome=="Preterm") 

table(table(metotalnmis2$subjectid)) # why here different between pos and neg
# incorrect coding after I compare the metotal for neg/pos
# E0616 in neg but not in pos; remove it
#metotalnmis2<-metotalnmis2[!metotalnmis2$Sample.ID=="E0616-1-Red2.1_001",]
table(table(metotalnmis2$subjectid)) 
# 1   2 
#42 252  252 appear twice, 42 appear once, matches HILIC data

#creat prenatal complication indicator variable---g_htn,iugr,g_dm,Depression_Preg
# metotalnmis$pre_complication<-ifelse(metotalnmis$g_htn!=0 | metotalnmis$g_dm==1 | 
#                                        metotalnmis$iugr==1 | metotalnmis$Depression_Preg==1,1,0)
# table(metotalnmis$pre_complication)

### create prenatal complication factor variable ---0,1,2 only include ghtn gdm
metotalnmis2$pre_complication<-ifelse(metotalnmis2$g_htn==1 & metotalnmis2$g_dm==1,2,
                                      ifelse(metotalnmis2$g_htn==2 & metotalnmis2$g_dm==1,2,
                                             ifelse(metotalnmis2$g_htn==0 & metotalnmis2$g_dm==0,0,1)))

#########log transform the dataset--unify the row observations
#Metotal2 = Demographics+ C18 basic column info + gestational age + nulliparity + C18 metabolic features
#Make Metotalnmis2 match metotalnmis2 in # of obs by applying same exclusions for stage, birth outcome, duplicate obs
Metotalnmis2<-subset(Metotal2,stage==1| stage==2)
Metotalnmis2<-subset(Metotalnmis2,birthoutcome=="EarlyTerm" | birthoutcome=="FullTerm" |
                       birthoutcome=="Preterm")
#Metotalnmis2<-Metotalnmis2[!Metotalnmis2$Sample.ID=="E0616-1-Red2.1_001",]

#metotalnmis2[,288:12187]<-log2(Metotalnmis2[,285:12184]) 
metotalnmis2[,194:13236]<-log2(Metotalnmis2[,191:13233])   
metotalnmis2[metotalnmis2 == -Inf] <- NA

#setwd("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/NEG/")
setwd("/Users/ethan_j_li/Documents/Emory/Research/2017/C18 2017 PreComp Regressions")

####### Need to create combined group--PTB+ETB=PETB
metotalnmis2$Comp_PTBETBvsHealthyFullTerm<-ifelse(metotalnmis2$Comp_AllETBvsHealthyFullTerm==1 | metotalnmis2$Comp_AllPTBvsHealthyFullTerm==1,1,
                                                  ifelse(metotalnmis2$Comp_AllETBvsHealthyFullTerm==0,0,NA))

####### Need to create combined group--miPTB+miETB=miPETB
metotalnmis2$Comp_InducedPTBETBvsHealthyFullTerm<-ifelse(metotalnmis2$Comp_InducedPTBvsHealthyFullTerm==1 | metotalnmis2$Comp_InducedETBvsHealthyFullTerm==1,1,
                                                         ifelse(metotalnmis2$Comp_InducedPTBvsHealthyFullTerm==0,0,NA))

###Stage 1
#Question: why separate analyses for stage 1 and stage 2 data?
metotals1_2<-subset(metotalnmis2,metotalnmis2$stage==1) #267+64-3

## sPTB vs All full term
c_SpPTBvsALLFullTerm_s1<-data.frame(matrix(nrow = nrow(c18_results), ncol = 4))
colnames(c_SpPTBvsALLFullTerm_s1)<-c("Est","Std","t-value","P")
#Metabolite data starts at col 194 = 1 + 193
#Change age --> age_enrollment for regression
#Regression fails on certain metabolic feats, not sure why? No intensity value from C18 LCMS analysis?
for (i in 1:nrow(c18_results) ) 
{
  tryCatch(
    {
      # lmfit<-lm( metotals1_2[,i+193] ~Comp_SpPTBvsALLFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1_2)
      lmfit<-lm( metotals1_2[,i+193] ~Comp_SpPTBvsALLFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1_2)
      c_SpPTBvsALLFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_SpPTBvsALLFullTerm_s1$pfdr<-p.adjust(c_SpPTBvsALLFullTerm_s1$P,method="BH")
write.csv(c_SpPTBvsALLFullTerm_s1,"c_SpPTBvsALLFullTerm_s1.csv")


## sPTB sETB vs Healthy full term
c_SpPTBSpETBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_SpPTBSpETBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+287] ~Comp_SpPTBSpETBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      c_SpPTBSpETBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_SpPTBSpETBvsHealthyFullTerm_s1$pfdr<-p.adjust(c_SpPTBSpETBvsHealthyFullTerm_s1$P,method="BH")
write.csv(c_SpPTBSpETBvsHealthyFullTerm_s1,"c_SpPTBSpETBvsHealthyFullTerm_s1.csv")


## sPTB  vs Healthy full term
c_SpPTBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_SpPTBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+287] ~Comp_SpPTBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      c_SpPTBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_SpPTBvsHealthyFullTerm_s1$pfdr<-p.adjust(c_SpPTBvsHealthyFullTerm_s1$P,method="BH")
write.csv(c_SpPTBvsHealthyFullTerm_s1,"c_SpPTBvsHealthyFullTerm_s1.csv")


## AllPTB  vs Healthy full term
c_AllPTBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_AllPTBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+287] ~Comp_AllPTBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      c_AllPTBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_AllPTBvsHealthyFullTerm_s1$pfdr<-p.adjust(c_AllPTBvsHealthyFullTerm_s1$P,method="BH")
write.csv(c_AllPTBvsHealthyFullTerm_s1,"c_AllPTBvsHealthyFullTerm_s1.csv")

## Induced PTB  vs Healthy full term
c_InducedPTBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_InducedPTBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+287] ~Comp_InducedPTBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      c_InducedPTBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_InducedPTBvsHealthyFullTerm_s1$pfdr<-p.adjust(c_InducedPTBvsHealthyFullTerm_s1$P,method="BH")
write.csv(c_InducedPTBvsHealthyFullTerm_s1,"c_InducedPTBvsHealthyFullTerm_s1.csv")

## All ETB  vs Healthy full term
c_AllETBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_AllETBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+287] ~Comp_AllETBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      c_AllETBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_AllETBvsHealthyFullTerm_s1$pfdr<-p.adjust(c_AllETBvsHealthyFullTerm_s1$P,method="BH")
write.csv(c_AllETBvsHealthyFullTerm_s1,"c_AllETBvsHealthyFullTerm_s1.csv")

## sETB  vs Healthy full term
c_SpETBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_SpETBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+287] ~Comp_SpETBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      c_SpETBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_SpETBvsHealthyFullTerm_s1$pfdr<-p.adjust(c_SpETBvsHealthyFullTerm_s1$P,method="BH")
write.csv(c_SpETBvsHealthyFullTerm_s1,"c_SpETBvsHealthyFullTerm_s1.csv")

## Induced ETB  vs Healthy full term
c_InducedETBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_InducedETBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+287] ~Comp_InducedETBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      c_InducedETBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_InducedETBvsHealthyFullTerm_s1$pfdr<-p.adjust(c_InducedETBvsHealthyFullTerm_s1$P,method="BH")
write.csv(c_InducedETBvsHealthyFullTerm_s1,"c_InducedETBvsHealthyFullTerm_s1.csv")

## PTB ETB vs Healthy full term
c_PTBETBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_PTBETBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+287] ~Comp_PTBETBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      c_PTBETBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_PTBETBvsHealthyFullTerm_s1$pfdr<-p.adjust(c_PTBETBvsHealthyFullTerm_s1$P,method="BH")
#c_PTBETBvsHealthyFullTerm_s1<-cbind(c20[,c(1:2)],c_PTBETBvsHealthyFullTerm_s1)
write.csv(c_PTBETBvsHealthyFullTerm_s1,"c_PTBETBvsHealthyFullTerm_s1.csv")

## Induced PTB ETB  vs Healthy full term
c_InducedPTBETBvsHealthyFullTerm_s1<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_InducedPTBETBvsHealthyFullTerm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1[,i+287] ~Comp_InducedPTBETBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1)
      c_InducedPTBETBvsHealthyFullTerm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_InducedPTBETBvsHealthyFullTerm_s1$pfdr<-p.adjust(c_InducedPTBETBvsHealthyFullTerm_s1$P,method="BH")
#c_InducedPTBETBvsHealthyFullTerm_s1<-cbind(c20[,c(1:2)],c_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(c_InducedPTBETBvsHealthyFullTerm_s1,"c_InducedPTBETBvsHealthyFullTerm_s1.csv")


#preg complications regression C18 Stage 1
c_pre_complication_s1<-data.frame(matrix(nrow = nrow(c18_results), ncol = 4))
colnames(c_pre_complication_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c18_results) ) 
{
  tryCatch(
    {
      # lmfit<-lm( metotals1_2[,i+193] ~Comp_SpPTBvsALLFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1_2)
      lmfit<-lm( metotals1_2[,i+193] ~pre_complication+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1_2)
      c_pre_complication_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_pre_complication_s1$pfdr<-p.adjust(c_pre_complication_s1$P,method="BH")
write.csv(c_pre_complication_s1,"c_pre_complication_s1_2017.csv")


#g_htn S1 C18 regression Done
c_g_htn_s1<-data.frame(matrix(nrow = nrow(c18_results), ncol = 4))
colnames(c_g_htn_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c18_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1_2[,i+193] ~g_htn+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1_2)
      c_g_htn_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_g_htn_s1$pfdr<-p.adjust(c_g_htn_s1$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(c_g_htn_s1,"c_g_htn_s1_2017.csv")

#g_dm S1 C18 regression Done
c_g_dm_s1<-data.frame(matrix(nrow = nrow(c18_results), ncol = 4))
colnames(c_g_dm_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c18_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1_2[,i+193] ~g_dm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1_2)
      c_g_dm_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_g_dm_s1$pfdr<-p.adjust(c_g_dm_s1$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(c_g_dm_s1,"c_g_dm_s1_2017.csv")

#PROM S1 C18 regression Done
c_PROM_s1<-data.frame(matrix(nrow = nrow(c18_results), ncol = 4))
colnames(c_PROM_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c18_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1_2[,i+193] ~PROM+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1_2)
      c_PROM_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_PROM_s1$pfdr<-p.adjust(c_PROM_s1$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(c_PROM_s1,"c_PROM_s1_2017.csv")

#IUGR S1 C18 regression Done
c_IUGR_s1<-data.frame(matrix(nrow = nrow(c18_results), ncol = 4))
colnames(c_IUGR_s1)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c18_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals1_2[,i+193] ~Comp_IUGRvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1_2)
      c_IUGR_s1[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_IUGR_s1$pfdr<-p.adjust(c_IUGR_s1$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(c_IUGR_s1,"c_IUGR_s1_2017.csv")

###Stage 2
metotals2_2<-subset(metotalnmis2,metotalnmis2$stage==2)

## sPTB vs All full term
c_SpPTBvsALLFullTerm_s2<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_SpPTBvsALLFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+287] ~Comp_SpPTBvsALLFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      c_SpPTBvsALLFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_SpPTBvsALLFullTerm_s2$pfdr<-p.adjust(c_SpPTBvsALLFullTerm_s2$P,method="BH")
write.csv(c_SpPTBvsALLFullTerm_s2,"c_SpPTBvsALLFullTerm_s2.csv")


## sPTB sETB vs Healthy full term
c_SpPTBSpETBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_SpPTBSpETBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+287] ~Comp_SpPTBSpETBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      c_SpPTBSpETBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_SpPTBSpETBvsHealthyFullTerm_s2$pfdr<-p.adjust(c_SpPTBSpETBvsHealthyFullTerm_s2$P,method="BH")
write.csv(c_SpPTBSpETBvsHealthyFullTerm_s2,"c_SpPTBSpETBvsHealthyFullTerm_s2.csv")


## sPTB  vs Healthy full term
c_SpPTBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_SpPTBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+287] ~Comp_SpPTBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      c_SpPTBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_SpPTBvsHealthyFullTerm_s2$pfdr<-p.adjust(c_SpPTBvsHealthyFullTerm_s2$P,method="BH")
write.csv(c_SpPTBvsHealthyFullTerm_s2,"c_SpPTBvsHealthyFullTerm_s2.csv")


## AllPTB  vs Healthy full term
c_AllPTBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_AllPTBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+287] ~Comp_AllPTBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      c_AllPTBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_AllPTBvsHealthyFullTerm_s2$pfdr<-p.adjust(c_AllPTBvsHealthyFullTerm_s2$P,method="BH")
write.csv(c_AllPTBvsHealthyFullTerm_s2,"c_AllPTBvsHealthyFullTerm_s2.csv")

## Induced PTB  vs Healthy full term
c_InducedPTBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_InducedPTBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+287] ~Comp_InducedPTBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      c_InducedPTBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_InducedPTBvsHealthyFullTerm_s2$pfdr<-p.adjust(c_InducedPTBvsHealthyFullTerm_s2$P,method="BH")
write.csv(c_InducedPTBvsHealthyFullTerm_s2,"c_InducedPTBvsHealthyFullTerm_s2.csv")

## ALL ETB vs Healthy full term
c_AllETBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_AllETBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+287] ~Comp_AllETBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      c_AllETBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_AllETBvsHealthyFullTerm_s2$pfdr<-p.adjust(c_AllETBvsHealthyFullTerm_s2$P,method="BH")
write.csv(c_AllETBvsHealthyFullTerm_s2,"c_AllETBvsHealthyFullTerm_s2.csv")

## sETB vs Healthy full term
c_SpETBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_SpETBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+287] ~Comp_SpETBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      c_SpETBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_SpETBvsHealthyFullTerm_s2$pfdr<-p.adjust(c_SpETBvsHealthyFullTerm_s2$P,method="BH")
write.csv(c_SpETBvsHealthyFullTerm_s2,"c_SpETBvsHealthyFullTerm_s2.csv")

## Induced ETB vs Healthy full term
c_InducedETBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_InducedETBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+287] ~Comp_InducedETBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      c_InducedETBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_InducedETBvsHealthyFullTerm_s2$pfdr<-p.adjust(c_InducedETBvsHealthyFullTerm_s2$P,method="BH")
write.csv(c_InducedETBvsHealthyFullTerm_s2,"c_InducedETBvsHealthyFullTerm_s2.csv")

## PTB ETB vs Healthy full term
c_PTBETBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_PTBETBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+287] ~Comp_PTBETBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      c_PTBETBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_PTBETBvsHealthyFullTerm_s2$pfdr<-p.adjust(c_PTBETBvsHealthyFullTerm_s2$P,method="BH")
#c_PTBETBvsHealthyFullTerm_s2<-cbind(c20[,c(1:2)],c_PTBETBvsHealthyFullTerm_s2)
write.csv(c_PTBETBvsHealthyFullTerm_s2,"c_PTBETBvsHealthyFullTerm_s2.csv")

## Induced PTB ETB  vs Healthy full term
c_InducedPTBETBvsHealthyFullTerm_s2<-data.frame(matrix(nrow = nrow(c20), ncol = 4))
colnames(c_InducedPTBETBvsHealthyFullTerm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c20) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2[,i+287] ~Comp_InducedPTBETBvsHealthyFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2)
      c_InducedPTBETBvsHealthyFullTerm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_InducedPTBETBvsHealthyFullTerm_s2$pfdr<-p.adjust(c_InducedPTBETBvsHealthyFullTerm_s2$P,method="BH")
c_InducedPTBETBvsHealthyFullTerm_s2<-cbind(c20[,c(1:2)],c_InducedPTBETBvsHealthyFullTerm_s2)
write.csv(c_InducedPTBETBvsHealthyFullTerm_s2,"c_InducedPTBETBvsHealthyFullTerm_s2.csv")

#g_htn S2 C18 regression Done
c_g_htn_s2<-data.frame(matrix(nrow = nrow(c18_results), ncol = 4))
colnames(c_g_htn_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c18_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2_2[,i+193] ~g_htn+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2_2)
      c_g_htn_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_g_htn_s2$pfdr<-p.adjust(c_g_htn_s2$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(c_g_htn_s2,"c_g_htn_s2_2017.csv")

#g_dm S2 C18 regression Done
c_g_dm_s2<-data.frame(matrix(nrow = nrow(c18_results), ncol = 4))
colnames(c_g_dm_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c18_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2_2[,i+193] ~g_dm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2_2)
      c_g_dm_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_g_dm_s2$pfdr<-p.adjust(c_g_dm_s2$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(c_g_dm_s2,"c_g_dm_s2_2017.csv")

#PROM S2 C18 regression Done
c_PROM_s2<-data.frame(matrix(nrow = nrow(c18_results), ncol = 4))
colnames(c_PROM_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c18_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2_2[,i+193] ~PROM+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2_2)
      c_PROM_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_PROM_s2$pfdr<-p.adjust(c_PROM_s2$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(c_PROM_s2,"c_PROM_s2_2017.csv")

#IUGR S2 C18 regression Done
c_IUGR_s2<-data.frame(matrix(nrow = nrow(c18_results), ncol = 4))
colnames(c_IUGR_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c18_results) ) 
{
  tryCatch(
    {
      lmfit<-lm( metotals2_2[,i+193] ~Comp_IUGRvsHealthyFullTerm+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2_2)
      c_IUGR_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_IUGR_s2$pfdr<-p.adjust(c_IUGR_s2$P,method="BH")
#h_InducedPTBETBvsHealthyFullTerm_s1<-cbind(h20[,c(1:2)],h_InducedPTBETBvsHealthyFullTerm_s1)
write.csv(c_IUGR_s2,"c_IUGR_s2_2017.csv")

#pre_complication S2 C18 Regression Current
c_pre_complication_s2<-data.frame(matrix(nrow = nrow(c18_results), ncol = 4))
colnames(c_pre_complication_s2)<-c("Est","Std","t-value","P")

for (i in 1:nrow(c18_results) ) 
{
  tryCatch(
    {
      # lmfit<-lm( metotals1_2[,i+193] ~Comp_SpPTBvsALLFullTerm+as.factor(Sex)+age+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals1_2)
      lmfit<-lm( metotals2_2[,i+193] ~pre_complication+as.factor(Sex)+age_enrollment+as.factor(nullp)+as.factor(Education_4.level)+FirstPrenatalBMI, data=metotals2_2)
      c_pre_complication_s2[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}
c_pre_complication_s2$pfdr<-p.adjust(c_pre_complication_s2$P,method="BH")
write.csv(c_pre_complication_s2,"c_pre_complication_s2_2017.csv")

################ Get statistics under different cut-off p-values
####POS
## import multiple csv files to a list
setwd("/Users/ethan_j_li/Documents/Emory/Research/2017/HILIC 2017 Regressions/")
list.files(pattern="h.*.csv$") 
list.filenames <- list.files(pattern="h.*.csv$")

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

colnames(lmfit_crude[[1]][6])

MWAStatsHILICPreComp <- as.data.frame(matrix(nrow = 10, ncol = 7))
colnames(MWAStatsHILICPreComp)<-c("Outcome","FDR0.05","FDR0.20","Raw0.0005", "Raw0.005","Raw0.01","Raw0.05")
for (i in 1:10)
{
  tryCatch(
    {
      MWAStatsHILICPreComp[i,1]<-names(lmfit_crude)[i]
      MWAStatsHILICPreComp[i,2]<-sum(lmfit_crude[[i]][6]<0.05,na.rm = TRUE)
      MWAStatsHILICPreComp[i,3]<-sum(lmfit_crude[[i]][6]<0.20,na.rm = TRUE)
      MWAStatsHILICPreComp[i,4]<-sum(lmfit_crude[[i]][5]<0.0005,na.rm = TRUE)
      MWAStatsHILICPreComp[i,5]<-sum(lmfit_crude[[i]][5]<0.005,na.rm = TRUE)
      MWAStatsHILICPreComp[i,6]<-sum(lmfit_crude[[i]][5]<0.01,na.rm = TRUE)
      MWAStatsHILICPreComp[i,7]<-sum(lmfit_crude[[i]][5]<0.05,na.rm = TRUE)
    },error=function(e){})
}

write.csv(MWAStatsHILICPreComp,"MWAStatsHILICPreComp_2017.csv")

#HILIC PreComp S1-S2 overlaps
overlap <- as.data.frame(matrix(nrow=5,ncol=7))
names(overlap) <- c("Comparison", "Overlap P<0.05", "Overlap P<0.01", "Overlap P<0.005", "Overlap P<0.0005","Overlap BHP<0.2","Overlap BHP<0.05")
a <- 1
for (i in seq(from=1,to=9,by=2))
{
  df1 <- lmfit_crude[[i]]
  
  
  df2 <- lmfit_crude[[i+1]]
  
  overlap[a,1] <- substring(list.filenames[i],1,nchar(list.filenames[i])-7)
  overlap[a,2] <- length(df1$X[df1$P<0.05 & is.na(df1$P)==F & df2$P<0.05 & is.na(df2$P)==F])
  overlap[a,3] <- length(df1$X[df1$P<0.01 & is.na(df1$P)==F & df2$P<0.01 & is.na(df2$P)==F])
  overlap[a,4] <- length(df1$X[df1$P<0.005 & is.na(df1$P)==F & df2$P<0.005 & is.na(df2$P)==F])
  overlap[a,5] <- length(df1$X[df1$P<0.0005 & is.na(df1$P)==F & df2$P<0.0005 & is.na(df2$P)==F])
  overlap[a,6] <- length(df1$X[df1$pdfr<0.2 & is.na(df1$pdfr)==F & df2$pdfr<0.2 & is.na(df2$pdfr)==F])
  overlap[a,7] <- length(df1$X[df1$pdfr<0.05 & is.na(df1$pdfr)==F & df2$pdfr<0.05 & is.na(df2$pdfr)==F])
  a <- a+1
}

write.csv(overlap,"HILIC_Precomp_Overlap_2017.csv")

#HILIC Cross-PreComp Overlaps


####NEG
## import multiple csv files to a list
setwd("/Users/ethan_j_li/Documents/Emory/Research/2017/C18 2017 Regressions/")
list.files(pattern="c.*.csv$") 
list.filenames <- list.files(pattern="c.*.csv$")


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


###Check out the sig num of the features in negative mode Question: positive/negative = HILIC/C18, referring to MS?

colnames(lmfit_crude[[1]][6])

MWAStatsC18PreComp <- as.data.frame(matrix(nrow = 10, ncol = 7))
colnames(MWAStatsC18PreComp)<-c("Outcome","FDR0.05","FDR0.20","Raw0.0005", "Raw0.005","Raw0.01","Raw0.05")
for (i in 1:10)
{
  tryCatch(
    {
      MWAStatsC18PreComp[i,1]<-names(lmfit_crude)[i]
      MWAStatsC18PreComp[i,2]<-sum(lmfit_crude[[i]][6]<0.05,na.rm = TRUE)
      MWAStatsC18PreComp[i,3]<-sum(lmfit_crude[[i]][6]<0.20,na.rm = TRUE)
      MWAStatsC18PreComp[i,4]<-sum(lmfit_crude[[i]][5]<0.0005,na.rm = TRUE)
      MWAStatsC18PreComp[i,5]<-sum(lmfit_crude[[i]][5]<0.005,na.rm = TRUE)
      MWAStatsC18PreComp[i,6]<-sum(lmfit_crude[[i]][5]<0.01,na.rm = TRUE)
      MWAStatsC18PreComp[i,7]<-sum(lmfit_crude[[i]][5]<0.05,na.rm = TRUE)
    },error=function(e){})
}

write.csv(MWAStatsC18PreComp,"MWAStatsC18PreComp_2017.csv")

#C18 PreComp S1-S2 Overlap
overlap <- as.data.frame(matrix(nrow=5,ncol=7))
names(overlap) <- c("Comparison", "Overlap P<0.05", "Overlap P<0.01", "Overlap P<0.005", "Overlap P<0.0005","Overlap BHP<0.2","Overlap BHP<0.05")
a <- 1
for (i in seq(from=1,to=9,by=2))
{
  df1 <- lmfit_crude[[i]]
  
  
  df2 <- lmfit_crude[[i+1]]
  
  overlap[a,1] <- substring(list.filenames[i],1,nchar(list.filenames[i])-7)
  overlap[a,2] <- length(df1$X[df1$P<0.05 & is.na(df1$P)==F & df2$P<0.05 & is.na(df2$P)==F])
  overlap[a,3] <- length(df1$X[df1$P<0.01 & is.na(df1$P)==F & df2$P<0.01 & is.na(df2$P)==F])
  overlap[a,4] <- length(df1$X[df1$P<0.005 & is.na(df1$P)==F & df2$P<0.005 & is.na(df2$P)==F])
  overlap[a,5] <- length(df1$X[df1$P<0.0005 & is.na(df1$P)==F & df2$P<0.0005 & is.na(df2$P)==F])
  overlap[a,6] <- length(df1$X[df1$pdfr<0.2 & is.na(df1$pdfr)==F & df2$pdfr<0.2 & is.na(df2$pdfr)==F])
  overlap[a,7] <- length(df1$X[df1$pdfr<0.05 & is.na(df1$pdfr)==F & df2$pdfr<0.05 & is.na(df2$pdfr)==F])
  a <- a+1
}

write.csv(overlap,"C18_Precomp_Overlap_2017.csv")

###Making TXT files for Mummichog Analyses-----raw P=0.05
###Positive Mode

## import multiple csv files to a list
setwd("/Users/ethan_j_li/Documents/Emory/Research/2017/HILIC 2017 PreComp Regressions")
list.files(pattern="h.*.csv$") 
list.filenames <- list.files(pattern="h.*.csv$")
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
polname

for (i in 1:10)
{
  print(i)
  print(polname[i])
  tryCatch(
    {
      #ACE_ref<-h20[1:2]   # Use after cleaning metabolites No.non-zero
      ACE_ref <- hilic_results[1:2]
      ACE_ref$p.value<-lmfit_crude[[i]][5]
      ACE_ref$t.value<-lmfit_crude[[i]][4]
      ACE_ref2<-matrix(as.numeric(unlist(ACE_ref)),nrow=nrow(ACE_ref))
      ACE_ref2<-ACE_ref2[complete.cases(ACE_ref2[ ,3]),]
      colnames(ACE_ref2)<-c("mz","time","p.value","t.value")
      
      write.table(ACE_ref2, paste0("/Users/ethan_j_li/Documents/Emory/Research/2017/",polname[i],".txt"),sep="\t",row.names=FALSE)
    },error=function(e){})
}


###Making TXT files for Mummichog Analyses----raw P=0.05
###Negative Mode

## import multiple csv files to a list
setwd("/Users/ethan_j_li/Documents/Emory/Research/2017/C18 2017 PreComp Regressions/")
list.files(pattern="c.*.csv$") 
list.filenames <- list.files(pattern="c.*.csv$")
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
polname

for (i in 1:10)
{
  print(i)
  print(polname[i])
  tryCatch(
    {
      #ACE_ref<-c20[1:2]
      ACE_ref<-c18_results[1:2]
      ACE_ref$p.value<-lmfit_crude[[i]][5]
      ACE_ref$t.value<-lmfit_crude[[i]][4]
      ACE_ref2<-matrix(as.numeric(unlist(ACE_ref)),nrow=nrow(ACE_ref))
      ACE_ref2<-ACE_ref2[complete.cases(ACE_ref2[ ,3]),]   # remove NA
      colnames(ACE_ref2)<-c("mz","time","p.value","t.value")
      
      write.table(ACE_ref2, paste0("/Users/ethan_j_li/Documents/Emory/Research/2017/",polname[i],".txt"),sep="\t",row.names=FALSE)
    },error=function(e){})
}


###Making Scripts for Mummichog---make sure no space and only underline!!

setwd("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/Mummichog_input")
list.files(pattern="*.txt$") 
list.filenames <- list.files(pattern=".txt$")
list.filenames

output<-unlist(substring(list.filenames[1],1,nchar(list.filenames)-4))[1]
test<-paste0("python /Users/tanyouran/Downloads/mummichog-1.0.9/mummichog/main.py -f /Users/tanyouran/Desktop/Liang/mummichog_results/Mummichog_input/",list.filenames[1]," -o ",output, 
             "_p0.05 -c 0.05 -m negative -z TRUE")
test

earthmum<-as.data.frame(matrix(nrow =20*2, ncol = 1))

###Neg @ 0.05
for (i in 1:20)
{
  output<-unlist(substring(list.filenames[i],1,nchar(list.filenames[i])-4))[1]
  earthmum[i,1]<-paste0("python /Users/tanyouran/Downloads/mummichog-1.0.9/mummichog/main.py -f /Users/tanyouran/Desktop/Liang/mummichog_results/Mummichog_input/",list.filenames[i]," -o ",output, 
                        "_p0.05 -c 0.05 -m negative -z TRUE")
}

###Pos @ 0.05
for (i in 21:40)
{
  output<-unlist(substring(list.filenames[i],1,nchar(list.filenames[i])-4))[1]
  earthmum[i,1]<-paste0("python /Users/tanyouran/Downloads/mummichog-1.0.9/mummichog/main.py -f /Users/tanyouran/Desktop/Liang/mummichog_results/Mummichog_input/",list.filenames[i]," -o ",output, 
                        "_p0.05 -c 0.05 -m positive -z TRUE")
}

setwd("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/Mummichog_output0.05/")
write.table(earthmum, "earthmum.txt",sep="\t",row.names=FALSE)


########## summary of significant features--generating heatmap--p0.05--ETB
library("stringr")

list.filenames <- list.files(path = "/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/Mummichog_output0.05/Summary_ETB/", pattern="mcg_pathwayanalysis.*.xlsx$",full.names = TRUE, recursive = FALSE) 
res<-substring(list.filenames,149,nchar(list.filenames)-11)
res

myData1 <- read_excel('/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/Mummichog_output0.05/Summary_ETB/mcg_pathwayanalysis_c_AllETBvsHealthyFullTerm_s1_p0.05.xlsx')[1:119,]
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
Mum_Pathway<-Mum_Pathway[order(-Mum_Pathway$Total0.05,-Mum_Pathway$Total0.10),]
write.csv(Mum_Pathway,"/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/Mummichog_output0.05/sETBsum_raw0.05.csv")

########## summary of significant features--generating heatmap--p0.05--PTB
list.filenames <- list.files(path = "/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/Mummichog_output0.05/Summary_PTB/", pattern="mcg_pathwayanalysis.*.xlsx$",full.names = TRUE, recursive = FALSE) 
res<-substring(list.filenames,149,nchar(list.filenames)-11)
res

myData1 <- read_excel('/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/Mummichog_output0.05/Summary_PTB/mcg_pathwayanalysis_c_AllPTBvsHealthyFullTerm_s1_p0.05.xlsx')[1:119,]
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
Mum_Pathway<-Mum_Pathway[order(-Mum_Pathway$Total0.05,-Mum_Pathway$Total0.10),]
write.csv(Mum_Pathway,"/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/Mummichog_output0.05/sPTBsum_raw0.05.csv")


################
##add two columns--mz and retention time to the MWAS output for further annotation
#POS
setwd("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/POS/")
list.files(pattern="h.*.csv$") 
list.filenames <- list.files(pattern="h.*.csv$")
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

for (i in 1:16)
{
  print(i)
  tryCatch(
    {
      ACE_com<-h20[1:2]
      ACE_com$Est<-lmfit_crude[[i]][2]
      ACE_com$Std<-lmfit_crude[[i]][3]
      ACE_com$t.value<-lmfit_crude[[i]][4]
      ACE_com$P<-lmfit_crude[[i]][5]
      ACE_com$pfdr<-lmfit_crude[[i]][6]
      ACE_com2<-matrix(as.numeric(unlist(ACE_com)),nrow=nrow(ACE_com))
      colnames(ACE_com2)<-c("mz","time","Est","Std","t.value","P","pfdr")
      write.csv(ACE_com2, paste0("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/POS/",polname[i],".csv"),na = "NA", row.names = F)
    },error=function(e){})
}

for (i in 1:length(list.filenames)) assign(list.filenames[i], read.csv(list.filenames[i]))

#NEG
setwd("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/NEG")
list.files(pattern="c.*.csv$") 
list.filenames <- list.files(pattern="c.*.csv$")
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

polname<-str_remove_all(list.filenames, ".csv")

for (i in 1:16)
{
  print(i)
  tryCatch(
    {
      ACE_com<-c20[1:2]
      ACE_com$Est<-lmfit_crude[[i]][2]
      ACE_com$Std<-lmfit_crude[[i]][3]
      ACE_com$t.value<-lmfit_crude[[i]][4]
      ACE_com$P<-lmfit_crude[[i]][5]
      ACE_com$pfdr<-lmfit_crude[[i]][6]
      ACE_com2<-matrix(as.numeric(unlist(ACE_com)),nrow=nrow(ACE_com))
      colnames(ACE_com2)<-c("mz","time","Est","Std","t.value","P","pfdr")
      write.csv(ACE_com2, paste0("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/NEG/",polname[i],".csv"),na = "NA", row.names = F)
    },error=function(e){})
}


for (i in 1:length(list.filenames)) assign(list.filenames[i], read.csv(list.filenames[i]))


###Annotation and Summary
HMSMS<- read.csv("/Users/tanyouran/Desktop/Liang/ECHO_AA_Cohort/Newborn/Annotation/Hilic total_confirmed.csv") ##See my file attachments--this is a list of authentic reference standard at HILIC column
CMSMS <- read.csv("/Users/tanyouran/Desktop/Liang/ECHO_AA_Cohort/Newborn/Annotation/c18 total_confirmed.csv") ##See my file attachments--this is a list of authentic reference standard at C18 column

##read in dataset
setwd("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/POS")
list.filenames <- list.files(pattern="h.*.csv$")
list.filenames
for (i in 1:length(list.filenames)) assign(list.filenames[i], read.csv(list.filenames[i]))


##h20
##s1
##p
MWASh20s1p<-h20[,1:2]
MWASh20s1p$h20_AllPTBvsHealthyFullTerm_s1<-h_AllPTBvsHealthyFullTerm_s1.csv[,6]
#MWASh20s1p$h20_SpPTBvsALLFullTerm_s1<-h_SpPTBvsALLFullTerm_s1.csv[,6]
MWASh20s1p$h20_SpPTBvsHealthyFullTerm_s1<-h_SpPTBvsHealthyFullTerm_s1.csv[,6]
MWASh20s1p$h20_InducedPTBvsHealthyFullTerm_s1<-h_InducedPTBvsHealthyFullTerm_s1.csv[,6]
MWASh20s1p$h20_AllETBvsHealthyFullTerm_s1<-h_AllETBvsHealthyFullTerm_s1.csv[,6]
MWASh20s1p$h20_SpETBvsHealthyFullTerm_s1<-h_SpETBvsHealthyFullTerm_s1.csv[,6]
MWASh20s1p$h20_InducedETBvsHealthyFullTerm_s1<-h_InducedETBvsHealthyFullTerm_s1.csv[,6]
MWASh20s1p$h_PTBETBvsHealthyFullTerm_s1<-h_PTBETBvsHealthyFullTerm_s1.csv[,6]
MWASh20s1p$h20_SpPTBSpETBvsHealthyFullTerm_s1<-h_SpPTBSpETBvsHealthyFullTerm_s1.csv[,6]
MWASh20s1p$h_InducedPTBETBvsHealthyFullTerm_s1<-h_InducedPTBETBvsHealthyFullTerm_s1.csv[,6]

MWASh20s1p$totalsig<-NA
MWASh20s1p$mediator<-0
for (i in 1:nrow(h20))
{
  tryCatch(
    {
      
      MWASh20s1p$totalsig[i]<-sum(MWASh20s1p[i,3:11]<=0.05)
      if(sum(MWASh20s1p[i,3:11]<=0.05)>0)
      {
        MWASh20s1p$mediator[i]<-1
      }}, error=function(e){})
}

#est
MWASh20s1p$h20_AllPTBvsHealthyFullTerm_s1_beta<-h_AllPTBvsHealthyFullTerm_s1.csv[,3]
MWASh20s1p$h20_SpPTBvsHealthyFullTerm_s1_beta<-h_SpPTBvsHealthyFullTerm_s1.csv[,3]
MWASh20s1p$h20_InducedPTBvsHealthyFullTerm_s1_beta<-h_InducedPTBvsHealthyFullTerm_s1.csv[,3]
#MWASh20s1p$h20_SpPTBvsALLFullTerm_s1_beta<-h_SpPTBvsALLFullTerm_s1.csv[,3]
MWASh20s1p$h20_AllETBvsHealthyFullTerm_s1_beta<-h_AllETBvsHealthyFullTerm_s1.csv[,3]
MWASh20s1p$h20_SpETBvsHealthyFullTerm_s1_beta<-h_SpETBvsHealthyFullTerm_s1.csv[,3]
MWASh20s1p$h20_InducedETBvsHealthyFullTerm_s1_beta<-h_InducedETBvsHealthyFullTerm_s1.csv[,3]
MWASh20s1p$h_PTBETBvsHealthyFullTerm_s1_beta<-h_PTBETBvsHealthyFullTerm_s1.csv[,3]
MWASh20s1p$h20_SpPTBSpETBvsHealthyFullTerm_s1_beta<-h_SpPTBSpETBvsHealthyFullTerm_s1.csv[,3]
MWASh20s1p$h_InducedPTBETBvsHealthyFullTerm_s1_beta<-h_InducedPTBETBvsHealthyFullTerm_s1.csv[,3]

#std
MWASh20s1p$h20_AllPTBvsHealthyFullTerm_s1_std<-h_AllPTBvsHealthyFullTerm_s1.csv[,4]
MWASh20s1p$h20_SpPTBvsHealthyFullTerm_s1_std<-h_SpPTBvsHealthyFullTerm_s1.csv[,4]
MWASh20s1p$h20_InducedPTBvsHealthyFullTerm_s1_std<-h_InducedPTBvsHealthyFullTerm_s1.csv[,4]
#MWASh20s1p$h20_SpPTBvsALLFullTerm_s1_std<-h_SpPTBvsALLFullTerm_s1.csv[,4]
MWASh20s1p$h20_AllETBvsHealthyFullTerm_s1_std<-h_AllETBvsHealthyFullTerm_s1.csv[,4]
MWASh20s1p$h20_SpETBvsHealthyFullTerm_s1_std<-h_SpETBvsHealthyFullTerm_s1.csv[,4]
MWASh20s1p$h20_InducedETBvsHealthyFullTerm_s1_std<-h_InducedETBvsHealthyFullTerm_s1.csv[,4]
MWASh20s1p$h_PTBETBvsHealthyFullTerm_s1_std<-h_PTBETBvsHealthyFullTerm_s1.csv[,4]
MWASh20s1p$h20_SpPTBSpETBvsHealthyFullTerm_s1_std<-h_SpPTBSpETBvsHealthyFullTerm_s1.csv[,4]
MWASh20s1p$h_InducedPTBETBvsHealthyFullTerm_s1_std<-h_InducedPTBETBvsHealthyFullTerm_s1.csv[,4]

Sig<-MWASh20s1p[,1:2] ##Column # for mz and retention time
LVL1<-HMSMS[,5:6] ##Column # for mz and retention time
masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
LVL1.sig<-slice(MWASh20s1p,masteroverlap$index.A)
LVL1.match<-slice(HMSMS,masteroverlap$index.B)
LVL1total<-cbind(LVL1.match,LVL1.sig)
LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
write.csv(MWASh20s1p,"MWASs1p.csv")
write.csv(LVL1total,"MWASs1p_lvl1.csv")


##h20
##s2
##p
MWASh20s2p<-h20[,1:2]
MWASh20s2p$h20_AllPTBvsHealthyFullTerm_s2<-h_AllPTBvsHealthyFullTerm_s2.csv[,6]
#MWASh20s2p$h20_SpPTBvsALLFullTerm_s2<-h_SpPTBvsALLFullTerm_s2.csv[,6]
MWASh20s2p$h20_SpPTBvsHealthyFullTerm_s2<-h_SpPTBvsHealthyFullTerm_s2.csv[,6]
MWASh20s2p$h20_InducedPTBvsHealthyFullTerm_s2<-h_InducedPTBvsHealthyFullTerm_s2.csv[,6]
MWASh20s2p$h20_AllETBvsHealthyFullTerm_s2<-h_AllETBvsHealthyFullTerm_s2.csv[,6]
MWASh20s2p$h20_SpETBvsHealthyFullTerm_s2<-h_SpETBvsHealthyFullTerm_s2.csv[,6]
MWASh20s2p$h20_InducedETBvsHealthyFullTerm_s2<-h_InducedETBvsHealthyFullTerm_s2.csv[,6]
MWASh20s2p$h_PTBETBvsHealthyFullTerm_s2<-h_PTBETBvsHealthyFullTerm_s2.csv[,6]
MWASh20s2p$h20_SpPTBSpETBvsHealthyFullTerm_s2<-h_SpPTBSpETBvsHealthyFullTerm_s2.csv[,6]
MWASh20s2p$h_InducedPTBETBvsHealthyFullTerm_s2<-h_InducedPTBETBvsHealthyFullTerm_s2.csv[,6]

MWASh20s2p$totalsig<-NA
MWASh20s2p$mediator<-0
for (i in 1:nrow(h20))
{
  tryCatch(
    {
      
      MWASh20s2p$totalsig[i]<-sum(MWASh20s2p[i,3:11]<=0.05)
      if(sum(MWASh20s2p[i,3:11]<=0.05)>0)
      {
        MWASh20s2p$mediator[i]<-1
      }}, error=function(e){})
}

#est
MWASh20s2p$h20_AllPTBvsHealthyFullTerm_s2_beta<-h_AllPTBvsHealthyFullTerm_s2.csv[,3]
MWASh20s2p$h20_SpPTBvsHealthyFullTerm_s2_beta<-h_SpPTBvsHealthyFullTerm_s2.csv[,3]
MWASh20s2p$h20_InducedPTBvsHealthyFullTerm_s2_beta<-h_InducedPTBvsHealthyFullTerm_s2.csv[,3]
#MWASh20s2p$h20_SpPTBvsALLFullTerm_s2_beta<-h_SpPTBvsALLFullTerm_s2.csv[,3]
MWASh20s2p$h20_AllETBvsHealthyFullTerm_s2_beta<-h_AllETBvsHealthyFullTerm_s2.csv[,3]
MWASh20s2p$h20_SpETBvsHealthyFullTerm_s2_beta<-h_SpETBvsHealthyFullTerm_s2.csv[,3]
MWASh20s2p$h20_InducedETBvsHealthyFullTerm_s2_beta<-h_InducedETBvsHealthyFullTerm_s2.csv[,3]
MWASh20s2p$h_PTBETBvsHealthyFullTerm_s2_beta<-h_PTBETBvsHealthyFullTerm_s2.csv[,3]
MWASh20s2p$h20_SpPTBSpETBvsHealthyFullTerm_s2_beta<-h_SpPTBSpETBvsHealthyFullTerm_s2.csv[,3]
MWASh20s2p$h_InducedPTBETBvsHealthyFullTerm_s2_beta<-h_InducedPTBETBvsHealthyFullTerm_s2.csv[,3]

#std
MWASh20s2p$h20_AllPTBvsHealthyFullTerm_s2_std<-h_AllPTBvsHealthyFullTerm_s2.csv[,4]
MWASh20s2p$h20_SpPTBvsHealthyFullTerm_s2_std<-h_SpPTBvsHealthyFullTerm_s2.csv[,4]
MWASh20s2p$h20_InducedPTBvsHealthyFullTerm_s2_std<-h_InducedPTBvsHealthyFullTerm_s2.csv[,4]
#MWASh20s2p$h20_SpPTBvsALLFullTerm_s2_std<-h_SpPTBvsALLFullTerm_s2.csv[,4]
MWASh20s2p$h20_AllETBvsHealthyFullTerm_s2_std<-h_AllETBvsHealthyFullTerm_s2.csv[,4]
MWASh20s2p$h20_SpETBvsHealthyFullTerm_s2_std<-h_SpETBvsHealthyFullTerm_s2.csv[,4]
MWASh20s2p$h20_InducedETBvsHealthyFullTerm_s2_std<-h_InducedETBvsHealthyFullTerm_s2.csv[,4]
MWASh20s2p$h_PTBETBvsHealthyFullTerm_s2_std<-h_PTBETBvsHealthyFullTerm_s2.csv[,4]
MWASh20s2p$h20_SpPTBSpETBvsHealthyFullTerm_s2_std<-h_SpPTBSpETBvsHealthyFullTerm_s2.csv[,4]
MWASh20s2p$h_InducedPTBETBvsHealthyFullTerm_s2_std<-h_InducedPTBETBvsHealthyFullTerm_s2.csv[,4]

Sig<-MWASh20s2p[,1:2] ##Column # for mz and retention time
LVL1<-HMSMS[,5:6] ##Column # for mz and retention time
masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
LVL1.sig<-slice(MWASh20s2p,masteroverlap$index.A)
LVL1.match<-slice(HMSMS,masteroverlap$index.B)
LVL1total<-cbind(LVL1.match,LVL1.sig)
LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
write.csv(MWASh20s2p,"MWASs2p.csv")
write.csv(LVL1total,"MWASs2p_lvl1.csv")


##c20
##read in dataset
setwd("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/NEG/")
list.filenames <- list.files(pattern="c.*.csv$")
list.filenames
for (i in 1:length(list.filenames)) assign(list.filenames[i], read.csv(list.filenames[i]))

##s1
##p
MWASc20s1p<-c20[,1:2]
MWASc20s1p$c20_AllPTBvsHealthyFullTerm_s1<-c_AllPTBvsHealthyFullTerm_s1.csv[,6]
MWASc20s1p$c20_SpPTBvsHealthyFullTerm_s1<-c_SpPTBvsHealthyFullTerm_s1.csv[,6]
MWASc20s1p$c20_InducedPTBvsHealthyFullTerm_s1<-c_InducedPTBvsHealthyFullTerm_s1.csv[,6]
#MWASc20s1p$c20_SpPTBvsALLFullTerm_s1<-c_SpPTBvsALLFullTerm_s1.csv[,6]
MWASc20s1p$c20_AllETBvsHealthyFullTerm_s1<-c_AllETBvsHealthyFullTerm_s1.csv[,6]
MWASc20s1p$c20_SpETBvsHealthyFullTerm_s1<-c_SpETBvsHealthyFullTerm_s1.csv[,6]
MWASc20s1p$c20_InducedETBvsHealthyFullTerm_s1<-c_InducedETBvsHealthyFullTerm_s1.csv[,6]
MWASc20s1p$c_PTBETBvsHealthyFullTerm_s1<-c_PTBETBvsHealthyFullTerm_s1.csv[,6]
MWASc20s1p$c20_SpPTBSpETBvsHealthyFullTerm_s1<-c_SpPTBSpETBvsHealthyFullTerm_s1.csv[,6]
MWASc20s1p$c_InducedPTBETBvsHealthyFullTerm_s1<-c_InducedPTBETBvsHealthyFullTerm_s1.csv[,6]

MWASc20s1p$totalsig<-NA
MWASc20s1p$mediator<-0
for (i in 1:nrow(c20))
{
  tryCatch(
    {
      
      MWASc20s1p$totalsig[i]<-sum(MWASc20s1p[i,3:11]<=0.05)
      if(sum(MWASc20s1p[i,3:11]<=0.05)>0)
      {
        MWASc20s1p$mediator[i]<-1
      }}, error=function(e){})
}

#est
MWASc20s1p$c20_AllPTBvsHealthyFullTerm_s1_beta<-c_AllPTBvsHealthyFullTerm_s1.csv[,3]
#MWASc20s1p$c20_SpPTBvsALLFullTerm_s1_beta<-c_SpPTBvsALLFullTerm_s1.csv[,3]
MWASc20s1p$c20_SpPTBvsHealthyFullTerm_s1_beta<-c_SpPTBvsHealthyFullTerm_s1.csv[,3]
MWASc20s1p$c20_InducedPTBvsHealthyFullTerm_s1_beta<-c_InducedPTBvsHealthyFullTerm_s1.csv[,3]
MWASc20s1p$c20_AllETBvsHealthyFullTerm_s1_beta<-c_AllETBvsHealthyFullTerm_s1.csv[,3]
MWASc20s1p$c20_SpETBvsHealthyFullTerm_s1_beta<-c_SpETBvsHealthyFullTerm_s1.csv[,3]
MWASc20s1p$c20_InducedETBvsHealthyFullTerm_s1_beta<-c_InducedETBvsHealthyFullTerm_s1.csv[,3]
MWASc20s1p$c_PTBETBvsHealthyFullTerm_s1_beta<-c_PTBETBvsHealthyFullTerm_s1.csv[,3]
MWASc20s1p$c20_SpPTBSpETBvsHealthyFullTerm_s1_beta<-c_SpPTBSpETBvsHealthyFullTerm_s1.csv[,3]
MWASc20s1p$c_InducedPTBETBvsHealthyFullTerm_s1_beta<-c_InducedPTBETBvsHealthyFullTerm_s1.csv[,3]

#std
MWASc20s1p$c20_AllPTBvsHealthyFullTerm_s1_std<-c_AllPTBvsHealthyFullTerm_s1.csv[,4]
#MWASc20s1p$c20_SpPTBvsALLFullTerm_s1_std<-c_SpPTBvsALLFullTerm_s1.csv[,4]
MWASc20s1p$c20_SpPTBvsHealthyFullTerm_s1_std<-c_SpPTBvsHealthyFullTerm_s1.csv[,4]
MWASc20s1p$c20_InducedPTBvsHealthyFullTerm_s1_std<-c_InducedPTBvsHealthyFullTerm_s1.csv[,4]
MWASc20s1p$c20_AllETBvsHealthyFullTerm_s1_std<-c_AllETBvsHealthyFullTerm_s1.csv[,4]
MWASc20s1p$c20_SpETBvsHealthyFullTerm_s1_std<-c_SpETBvsHealthyFullTerm_s1.csv[,4]
MWASc20s1p$c20_InducedETBvsHealthyFullTerm_s1_std<-c_InducedETBvsHealthyFullTerm_s1.csv[,4]
MWASc20s1p$c_PTBETBvsHealthyFullTerm_s1_std<-c_PTBETBvsHealthyFullTerm_s1.csv[,4]
MWASc20s1p$c20_SpPTBSpETBvsHealthyFullTerm_s1_std<-c_SpPTBSpETBvsHealthyFullTerm_s1.csv[,4]
MWASc20s1p$c_InducedPTBETBvsHealthyFullTerm_s1_std<-c_InducedPTBETBvsHealthyFullTerm_s1.csv[,4]

Sig<-MWASc20s1p[,1:2] ##Column # for mz and retention time
LVL1<-CMSMS[,5:6] ##Column # for mz and retention time
masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
LVL1.sig<-slice(MWASc20s1p,masteroverlap$index.A)
LVL1.match<-slice(CMSMS,masteroverlap$index.B)
LVL1total<-cbind(LVL1.match,LVL1.sig)
LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
write.csv(MWASc20s1p,"MWASs1n.csv")
write.csv(LVL1total,"MWASs1n_lvl1.csv")


##c20
##s2
##p
MWASc20s2p<-c20[,1:2]
MWASc20s2p$c20_AllPTBvsHealthyFullTerm_s2<-c_AllPTBvsHealthyFullTerm_s2.csv[,6]
MWASc20s2p$c20_SpPTBvsHealthyFullTerm_s2<-c_SpPTBvsHealthyFullTerm_s2.csv[,6]
MWASc20s2p$c20_InducedPTBvsHealthyFullTerm_s2<-c_InducedPTBvsHealthyFullTerm_s2.csv[,6]
#MWASc20s2p$c20_SpPTBvsALLFullTerm_s2<-c_SpPTBvsALLFullTerm_s2.csv[,6]
MWASc20s2p$c20_AllETBvsHealthyFullTerm_s2<-c_AllETBvsHealthyFullTerm_s2.csv[,6]
MWASc20s2p$c20_SpETBvsHealthyFullTerm_s2<-c_SpETBvsHealthyFullTerm_s2.csv[,6]
MWASc20s2p$c20_InducedETBvsHealthyFullTerm_s2<-c_InducedETBvsHealthyFullTerm_s2.csv[,6]
MWASc20s2p$c_PTBETBvsHealthyFullTerm_s2<-c_PTBETBvsHealthyFullTerm_s2.csv[,6]
MWASc20s2p$c20_SpPTBSpETBvsHealthyFullTerm_s2<-c_SpPTBSpETBvsHealthyFullTerm_s2.csv[,6]
MWASc20s2p$c_InducedPTBETBvsHealthyFullTerm_s2<-c_InducedPTBETBvsHealthyFullTerm_s2.csv[,6]

MWASc20s2p$totalsig<-NA
MWASc20s2p$mediator<-0
for (i in 1:nrow(c20))
{
  tryCatch(
    {
      
      MWASc20s2p$totalsig[i]<-sum(MWASc20s2p[i,3:11]<=0.05)
      if(sum(MWASc20s2p[i,3:11]<=0.05)>0)
      {
        MWASc20s2p$mediator[i]<-1
      }}, error=function(e){})
}

#est
MWASc20s2p$c20_AllPTBvsHealthyFullTerm_s2_beta<-c_AllPTBvsHealthyFullTerm_s2.csv[,3]
#MWASc20s2p$c20_SpPTBvsALLFullTerm_s2_beta<-c_SpPTBvsALLFullTerm_s2.csv[,3]
MWASc20s2p$c20_SpPTBvsHealthyFullTerm_s2_beta<-c_SpPTBvsHealthyFullTerm_s2.csv[,3]
MWASc20s2p$c20_InducedPTBvsHealthyFullTerm_s2_beta<-c_InducedPTBvsHealthyFullTerm_s2.csv[,3]
MWASc20s2p$c20_AllETBvsHealthyFullTerm_s2_beta<-c_AllETBvsHealthyFullTerm_s2.csv[,3]
MWASc20s2p$c20_SpETBvsHealthyFullTerm_s2_beta<-c_SpETBvsHealthyFullTerm_s2.csv[,3]
MWASc20s2p$c20_InducedETBvsHealthyFullTerm_s2_beta<-c_InducedETBvsHealthyFullTerm_s2.csv[,3]
MWASc20s2p$c_PTBETBvsHealthyFullTerm_s2_beta<-c_PTBETBvsHealthyFullTerm_s2.csv[,3]
MWASc20s2p$c20_SpPTBSpETBvsHealthyFullTerm_s2_beta<-c_SpPTBSpETBvsHealthyFullTerm_s2.csv[,3]
MWASc20s2p$c_InducedPTBETBvsHealthyFullTerm_s2_beta<-c_InducedPTBETBvsHealthyFullTerm_s2.csv[,3]

#std
MWASc20s2p$c20_AllPTBvsHealthyFullTerm_s2_std<-c_AllPTBvsHealthyFullTerm_s2.csv[,4]
#MWASc20s2p$c20_SpPTBvsALLFullTerm_s2_std<-c_SpPTBvsALLFullTerm_s2.csv[,4]
MWASc20s2p$c20_SpPTBvsHealthyFullTerm_s2_std<-c_SpPTBvsHealthyFullTerm_s2.csv[,4]
MWASc20s2p$c20_InducedPTBvsHealthyFullTerm_s2_std<-c_InducedPTBvsHealthyFullTerm_s2.csv[,4]
MWASc20s2p$c20_AllETBvsHealthyFullTerm_s2_std<-c_AllETBvsHealthyFullTerm_s2.csv[,4]
MWASc20s2p$c20_SpETBvsHealthyFullTerm_s2_std<-c_SpETBvsHealthyFullTerm_s2.csv[,4]
MWASc20s2p$c20_InducedETBvsHealthyFullTerm_s2_std<-c_InducedETBvsHealthyFullTerm_s2.csv[,4]
MWASc20s2p$c_PTBETBvsHealthyFullTerm_s2_std<-c_PTBETBvsHealthyFullTerm_s2.csv[,4]
MWASc20s2p$c20_SpPTBSpETBvsHealthyFullTerm_s2_std<-c_SpPTBSpETBvsHealthyFullTerm_s2.csv[,4]
MWASc20s2p$c_InducedPTBETBvsHealthyFullTerm_s2_std<-c_InducedPTBETBvsHealthyFullTerm_s2.csv[,4]

Sig<-MWASc20s2p[,1:2] ##Column # for mz and retention time
LVL1<-CMSMS[,5:6] ##Column # for mz and retention time
masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
LVL1.sig<-slice(MWASc20s2p,masteroverlap$index.A)
LVL1.match<-slice(CMSMS,masteroverlap$index.B)
LVL1total<-cbind(LVL1.match,LVL1.sig)
LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
write.csv(MWASc20s2p,"MWASs2n.csv")
write.csv(LVL1total,"MWASs2n_lvl1.csv")


#################################
######################## Two newly added models PTBETB, miPTBETB results are summarized separately
###################
###Annotation and Summary
HMSMS<- read.csv("/Users/tanyouran/Desktop/Liang/ECHO_AA_Cohort/Newborn/Annotation/Hilic total_confirmed.csv") ##See my file attachments--this is a list of authentic reference standard at HILIC column
CMSMS <- read.csv("/Users/tanyouran/Desktop/Liang/ECHO_AA_Cohort/Newborn/Annotation/c18 total_confirmed.csv") ##See my file attachments--this is a list of authentic reference standard at C18 column

##read in dataset
setwd("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/POS")
list.filenames <- list.files(pattern="h.*.csv$")
list.filenames
for (i in 1:length(list.filenames)) assign(list.filenames[i], read.csv(list.filenames[i]))

##h20
##s1
##p
MWASh20s1p<-h20[,1:2]
MWASh20s1p$h_PTBETBvsHealthyFullTerm_s1<-h_PTBETBvsHealthyFullTerm_s1.csv[,6]
MWASh20s1p$h_InducedPTBETBvsHealthyFullTerm_s1<-h_InducedPTBETBvsHealthyFullTerm_s1.csv[,6]

MWASh20s1p$totalsig<-NA
for (i in 1:nrow(h20))
{
  tryCatch(
    {
      
      MWASh20s1p$totalsig[i]<-sum(MWASh20s1p[i,3:4]<=0.05)
    }, error=function(e){})
}

#est
MWASh20s1p$h_PTBETBvsHealthyFullTerm_s1_beta<-h_PTBETBvsHealthyFullTerm_s1.csv[,3]
MWASh20s1p$h_InducedPTBETBvsHealthyFullTerm_s1_beta<-h_InducedPTBETBvsHealthyFullTerm_s1.csv[,3]

#std
MWASh20s1p$h_PTBETBvsHealthyFullTerm_s1_std<-h_PTBETBvsHealthyFullTerm_s1.csv[,4]
MWASh20s1p$h_InducedPTBETBvsHealthyFullTerm_s1_std<-h_InducedPTBETBvsHealthyFullTerm_s1.csv[,4]

Sig<-MWASh20s1p[,1:2] ##Column # for mz and retention time
LVL1<-HMSMS[,5:6] ##Column # for mz and retention time
masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
LVL1.sig<-slice(MWASh20s1p,masteroverlap$index.A)
LVL1.match<-slice(HMSMS,masteroverlap$index.B)
LVL1total<-cbind(LVL1.match,LVL1.sig)
LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
write.csv(LVL1total,"MWASs1p_lvl1.csv")

##h20
##s2
##p
MWASh20s2p<-h20[,1:2]
MWASh20s2p$h_PTBETBvsHealthyFullTerm_s2<-h_PTBETBvsHealthyFullTerm_s2.csv[,6]
MWASh20s2p$h_InducedPTBETBvsHealthyFullTerm_s2<-h_InducedPTBETBvsHealthyFullTerm_s2.csv[,6]

MWASh20s2p$totalsig<-NA
for (i in 1:nrow(h20))
{
  tryCatch(
    {
      
      MWASh20s2p$totalsig[i]<-sum(MWASh20s2p[i,3:4]<=0.05)
    }, error=function(e){})
}

#est
MWASh20s2p$h_PTBETBvsHealthyFullTerm_s2_beta<-h_PTBETBvsHealthyFullTerm_s2.csv[,3]
MWASh20s2p$h_InducedPTBETBvsHealthyFullTerm_s2_beta<-h_InducedPTBETBvsHealthyFullTerm_s2.csv[,3]

#std
MWASh20s2p$h_PTBETBvsHealthyFullTerm_s2_std<-h_PTBETBvsHealthyFullTerm_s2.csv[,4]
MWASh20s2p$h_InducedPTBETBvsHealthyFullTerm_s2_std<-h_InducedPTBETBvsHealthyFullTerm_s2.csv[,4]

Sig<-MWASh20s2p[,1:2] ##Column # for mz and retention time
LVL1<-HMSMS[,5:6] ##Column # for mz and retention time
masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
LVL1.sig<-slice(MWASh20s2p,masteroverlap$index.A)
LVL1.match<-slice(HMSMS,masteroverlap$index.B)
LVL1total<-cbind(LVL1.match,LVL1.sig)
LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
write.csv(LVL1total,"MWASs2p_lvl1.csv")

##c20
##read in dataset
setwd("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/NEG/")
list.filenames <- list.files(pattern="c.*.csv$")
list.filenames
for (i in 1:length(list.filenames)) assign(list.filenames[i], read.csv(list.filenames[i]))

##s1
##p
MWASc20s1p<-c20[,1:2]
MWASc20s1p$c_PTBETBvsHealthyFullTerm_s1<-c_PTBETBvsHealthyFullTerm_s1.csv[,6]
MWASc20s1p$c_InducedPTBETBvsHealthyFullTerm_s1<-c_InducedPTBETBvsHealthyFullTerm_s1.csv[,6]

MWASc20s1p$totalsig<-NA

for (i in 1:nrow(c20))
{
  tryCatch(
    {
      
      MWASc20s1p$totalsig[i]<-sum(MWASc20s1p[i,3:4]<=0.05)
    }, error=function(e){})
}

#est
MWASc20s1p$c_PTBETBvsHealthyFullTerm_s1_beta<-c_PTBETBvsHealthyFullTerm_s1.csv[,3]
MWASc20s1p$c_InducedPTBETBvsHealthyFullTerm_s1_beta<-c_InducedPTBETBvsHealthyFullTerm_s1.csv[,3]

#std
MWASc20s1p$c_PTBETBvsHealthyFullTerm_s1_std<-c_PTBETBvsHealthyFullTerm_s1.csv[,4]
MWASc20s1p$c_InducedPTBETBvsHealthyFullTerm_s1_std<-c_InducedPTBETBvsHealthyFullTerm_s1.csv[,4]

Sig<-MWASc20s1p[,1:2] ##Column # for mz and retention time
LVL1<-CMSMS[,5:6] ##Column # for mz and retention time
masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
LVL1.sig<-slice(MWASc20s1p,masteroverlap$index.A)
LVL1.match<-slice(CMSMS,masteroverlap$index.B)
LVL1total<-cbind(LVL1.match,LVL1.sig)
LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
write.csv(LVL1total,"MWASs1n_lvl1.csv")


##s2
##p
MWASc20s2p<-c20[,1:2]
MWASc20s2p$c_PTBETBvsHealthyFullTerm_s2<-c_PTBETBvsHealthyFullTerm_s2.csv[,6]
MWASc20s2p$c_InducedPTBETBvsHealthyFullTerm_s2<-c_InducedPTBETBvsHealthyFullTerm_s2.csv[,6]

MWASc20s2p$totalsig<-NA

for (i in 1:nrow(c20))
{
  tryCatch(
    {
      
      MWASc20s2p$totalsig[i]<-sum(MWASc20s2p[i,3:4]<=0.05)
    }, error=function(e){})
}

#est
MWASc20s2p$c_PTBETBvsHealthyFullTerm_s2_beta<-c_PTBETBvsHealthyFullTerm_s2.csv[,3]
MWASc20s2p$c_InducedPTBETBvsHealthyFullTerm_s2_beta<-c_InducedPTBETBvsHealthyFullTerm_s2.csv[,3]

#std
MWASc20s2p$c_PTBETBvsHealthyFullTerm_s2_std<-c_PTBETBvsHealthyFullTerm_s2.csv[,4]
MWASc20s2p$c_InducedPTBETBvsHealthyFullTerm_s2_std<-c_InducedPTBETBvsHealthyFullTerm_s2.csv[,4]

Sig<-MWASc20s2p[,1:2] ##Column # for mz and retention time
LVL1<-CMSMS[,5:6] ##Column # for mz and retention time
masteroverlap <- find.Overlapping.mzs(Sig,LVL1, mz.thresh = 10,alignment.tool = 'apLCMS')
LVL1.sig<-slice(MWASc20s2p,masteroverlap$index.A)
LVL1.match<-slice(CMSMS,masteroverlap$index.B)
LVL1total<-cbind(LVL1.match,LVL1.sig)
LVL1total<-subset(LVL1total,LVL1total$totalsig>0)
LVL1total<-LVL1total[!duplicated(LVL1total[,2]),]  
write.csv(LVL1total,"MWASs2n_lvl1.csv")

##########################
#####I need to rerun previous annotation code from start to get a four large file containing P, est, std
### then restrict the mz to the following to make the chemical annotatio table
## now I have MWASh20s1p, MWASh20s2p, MWASc20s1p,MWASc20s2p
## select specific mz
c1<-c("101.0246","116.0718","128.0355","129.0195","129.0558","130.0874","135.0304","153.0193",
      "174.0885","175.0248","180.0665","187.1338","191.0196","197.0434","355.2853")
c2<-c("128.0355","129.0195","135.0304","143.1077","164.0716","167.0208","174.0885","191.0196",
      "267.0732")
h1<-c("104.1072","126.022","132.0656","150.0583","153.0396","165.0546","182.0812","200.1282",
      "241.0312","284.2947","315.2331","347.2215","357.3003","361.201","391.2842")
h2<-c("134.0448","149.0232","176.0707","200.1282","260.1855","315.2331","316.248","344.2792",
      "347.2215","357.3003","385.3466")
MWASh20s1p<-MWASh20s1p[MWASh20s1p$mz %in% h1,]
MWASh20s2p<-MWASh20s2p[MWASh20s2p$mz %in% h2,]
MWASc20s1p<-MWASc20s1p[MWASc20s1p$mz %in% c1,]
MWASc20s2p<-MWASc20s2p[MWASc20s2p$mz %in% c2,]
write.csv(MWASh20s1p,"MWASh20s1p_chemtable.csv",row.names = F)
write.csv(MWASh20s2p,"MWASh20s2p_chemtable.csv",row.names = F)
write.csv(MWASc20s1p,"MWASc20s1p_chemtable.csv",row.names = F)
write.csv(MWASc20s2p,"MWASc20s2p_chemtable.csv",row.names = F)

#############Calculate overlapped feature in Visit1&2
Pos<-cbind(MWASs1p[,], MWASs2p[,])
nrow(subset(Pos,Pos$h20_AllPTBvsHealthyFullTerm_s1<0.0005 & Pos$h20_AllPTBvsHealthyFullTerm_s2<0.005))
nrow(subset(Pos,Pos$h20_AllPTBvsHealthyFullTerm_s1<0.005 & Pos$h20_AllPTBvsHealthyFullTerm_s2<0.005))
nrow(subset(Pos,Pos$h20_AllPTBvsHealthyFullTerm_s1<0.01 & Pos$h20_AllPTBvsHealthyFullTerm_s2<0.01))
nrow(subset(Pos,Pos$h20_AllPTBvsHealthyFullTerm_s1<0.05 & Pos$h20_AllPTBvsHealthyFullTerm_s2<0.05))

nrow(subset(Pos,Pos$h20_SpPTBSpETBvsHealthyFullTerm_s1<0.0005 & Pos$h20_SpPTBSpETBvsHealthyFullTerm_s2<0.005))
nrow(subset(Pos,Pos$h20_SpPTBSpETBvsHealthyFullTerm_s1<0.005 & Pos$h20_SpPTBSpETBvsHealthyFullTerm_s2<0.005))
nrow(subset(Pos,Pos$h20_SpPTBSpETBvsHealthyFullTerm_s1<0.01 & Pos$h20_SpPTBSpETBvsHealthyFullTerm_s2<0.01))
nrow(subset(Pos,Pos$h20_SpPTBSpETBvsHealthyFullTerm_s1<0.05 & Pos$h20_SpPTBSpETBvsHealthyFullTerm_s2<0.05))

nrow(subset(Pos,Pos$h20_SpPTBvsALLFullTerm_s1<0.0005 & Pos$h20_SpPTBvsALLFullTerm_s2<0.005))
nrow(subset(Pos,Pos$h20_SpPTBvsALLFullTerm_s1<0.005 & Pos$h20_SpPTBvsALLFullTerm_s2<0.005))
nrow(subset(Pos,Pos$h20_SpPTBvsALLFullTerm_s1<0.01 & Pos$h20_SpPTBvsALLFullTerm_s2<0.01))
nrow(subset(Pos,Pos$h20_SpPTBvsALLFullTerm_s1<0.05 & Pos$h20_SpPTBvsALLFullTerm_s2<0.05))

nrow(subset(Pos,Pos$h20_SpPTBvsHealthyFullTerm_s1<0.0005 & Pos$h20_SpPTBvsHealthyFullTerm_s2<0.005))
nrow(subset(Pos,Pos$h20_SpPTBvsHealthyFullTerm_s1<0.005 & Pos$h20_SpPTBvsHealthyFullTerm_s2<0.005))
nrow(subset(Pos,Pos$h20_SpPTBvsHealthyFullTerm_s1<0.01 & Pos$h20_SpPTBvsHealthyFullTerm_s2<0.01))
nrow(subset(Pos,Pos$h20_SpPTBvsHealthyFullTerm_s1<0.05 & Pos$h20_SpPTBvsHealthyFullTerm_s2<0.05))

Neg<-cbind(MWASs1n[,], MWASs2n[,])
nrow(subset(Neg,Neg$c20_AllPTBvsHealthyFullTerm_s1<0.0005 & Neg$c20_AllPTBvsHealthyFullTerm_s2<0.005))
nrow(subset(Neg,Neg$c20_AllPTBvsHealthyFullTerm_s1<0.005 & Neg$c20_AllPTBvsHealthyFullTerm_s2<0.005))
nrow(subset(Neg,Neg$c20_AllPTBvsHealthyFullTerm_s1<0.01 & Neg$c20_AllPTBvsHealthyFullTerm_s2<0.01))
nrow(subset(Neg,Neg$c20_AllPTBvsHealthyFullTerm_s1<0.05 & Neg$c20_AllPTBvsHealthyFullTerm_s2<0.05))

nrow(subset(Neg,Neg$c20_SpPTBSpETBvsHealthyFullTerm_s1<0.0005 & Neg$c20_SpPTBSpETBvsHealthyFullTerm_s2<0.005))
nrow(subset(Neg,Neg$c20_SpPTBSpETBvsHealthyFullTerm_s1<0.005 & Neg$c20_SpPTBSpETBvsHealthyFullTerm_s2<0.005))
nrow(subset(Neg,Neg$c20_SpPTBSpETBvsHealthyFullTerm_s1<0.01 & Neg$c20_SpPTBSpETBvsHealthyFullTerm_s2<0.01))
nrow(subset(Neg,Neg$c20_SpPTBSpETBvsHealthyFullTerm_s1<0.05 & Neg$c20_SpPTBSpETBvsHealthyFullTerm_s2<0.05))

nrow(subset(Neg,Neg$c20_SpPTBvsALLFullTerm_s1<0.0005 & Neg$c20_SpPTBvsALLFullTerm_s2<0.005))
nrow(subset(Neg,Neg$c20_SpPTBvsALLFullTerm_s1<0.005 & Neg$c20_SpPTBvsALLFullTerm_s2<0.005))
nrow(subset(Neg,Neg$c20_SpPTBvsALLFullTerm_s1<0.01 & Neg$c20_SpPTBvsALLFullTerm_s2<0.01))
nrow(subset(Neg,Neg$c20_SpPTBvsALLFullTerm_s1<0.05 & Neg$c20_SpPTBvsALLFullTerm_s2<0.05))

nrow(subset(Neg,Neg$c20_SpPTBvsHealthyFullTerm_s1<0.0005 & Neg$c20_SpPTBvsHealthyFullTerm_s2<0.005))
nrow(subset(Neg,Neg$c20_SpPTBvsHealthyFullTerm_s1<0.005 & Neg$c20_SpPTBvsHealthyFullTerm_s2<0.005))
nrow(subset(Neg,Neg$c20_SpPTBvsHealthyFullTerm_s1<0.01 & Neg$c20_SpPTBvsHealthyFullTerm_s2<0.01))
nrow(subset(Neg,Neg$c20_SpPTBvsHealthyFullTerm_s1<0.05 & Neg$c20_SpPTBvsHealthyFullTerm_s2<0.05))



library(stringr)
setwd("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/NEG/")
list.files(pattern="c.*.csv$") 
list.filenames <- list.files(pattern="c.*.csv$")
polname<-str_remove_all(list.filenames, ".csv")

for (i in 1:length(polname)) 
  assign(polname[i], read.csv(list.filenames[i]))
Neg<-as.data.frame(cbind(c_AllETBvsHealthyFullTerm_s1$P,c_AllETBvsHealthyFullTerm_s2$P,c_AllPTBvsHealthyFullTerm_s1$P,
                         c_AllPTBvsHealthyFullTerm_s2$P,c_InducedETBvsHealthyFullTerm_s1$P,c_InducedETBvsHealthyFullTerm_s2$P,
                         c_InducedPTBvsHealthyFullTerm_s1$P,c_InducedPTBvsHealthyFullTerm_s2$P,c_SpETBvsHealthyFullTerm_s1$P,
                         c_SpETBvsHealthyFullTerm_s2$P,c_SpPTBSpETBvsHealthyFullTerm_s1$P,c_SpPTBSpETBvsHealthyFullTerm_s2$P,
                         c_SpPTBvsALLFullTerm_s1$P,c_SpPTBvsALLFullTerm_s2$P,c_SpPTBvsHealthyFullTerm_s1$P,
                         c_SpPTBvsHealthyFullTerm_s2$P))
colnames(Neg)<-c("ETBH1","ETBH2","PTBH1","PTBH2","miETBH1","miETBH2","miPTBH1","miPTBH2","sETBH1","sETBH2",
                 "sPETBH1","sPETBH2","sPTBA1","sPTBA2","sPTBH1","sPTBH2")

nrow(subset(Neg,Neg$ETBH1<0.0005 & Neg$ETBH2<0.0005))
nrow(subset(Neg,Neg$ETBH1<0.005 & Neg$ETBH2<0.005))
nrow(subset(Neg,Neg$ETBH1<0.01 & Neg$ETBH2<0.01))
nrow(subset(Neg,Neg$ETBH1<0.05 & Neg$ETBH2<0.05))

nrow(subset(Neg,Neg$sETBH1<0.0005 & Neg$sETBH2<0.0005))
nrow(subset(Neg,Neg$sETBH1<0.005 & Neg$sETBH2<0.005))
nrow(subset(Neg,Neg$sETBH1<0.01 & Neg$sETBH2<0.01))
nrow(subset(Neg,Neg$sETBH1<0.05 & Neg$sETBH2<0.05))

nrow(subset(Neg,Neg$miETBH1<0.0005 & Neg$miETBH2<0.0005))
nrow(subset(Neg,Neg$miETBH1<0.005 & Neg$miETBH2<0.005))
nrow(subset(Neg,Neg$miETBH1<0.01 & Neg$miETBH2<0.01))
nrow(subset(Neg,Neg$miETBH1<0.05 & Neg$miETBH2<0.05))

nrow(subset(Neg,Neg$PTBH1<0.0005 & Neg$PTBH2<0.0005))
nrow(subset(Neg,Neg$PTBH1<0.005 & Neg$PTBH2<0.005))
nrow(subset(Neg,Neg$PTBH1<0.01 & Neg$PTBH2<0.01))
nrow(subset(Neg,Neg$PTBH1<0.05 & Neg$PTBH2<0.05))

nrow(subset(Neg,Neg$sPTBH1<0.0005 & Neg$sPTBH2<0.0005))
nrow(subset(Neg,Neg$sPTBH1<0.005 & Neg$sPTBH2<0.005))
nrow(subset(Neg,Neg$sPTBH1<0.01 & Neg$sPTBH2<0.01))
nrow(subset(Neg,Neg$sPTBH1<0.05 & Neg$sPTBH2<0.05))

nrow(subset(Neg,Neg$miPTBH1<0.0005 & Neg$miPTBH2<0.0005))
nrow(subset(Neg,Neg$miPTBH1<0.005 & Neg$miPTBH2<0.005))
nrow(subset(Neg,Neg$miPTBH1<0.01 & Neg$miPTBH2<0.01))
nrow(subset(Neg,Neg$miPTBH1<0.05 & Neg$miPTBH2<0.05))

nrow(subset(Neg,Neg$sPETBH1<0.0005 & Neg$sPETBH2<0.0005))
nrow(subset(Neg,Neg$sPETBH1<0.005 & Neg$sPETBH2<0.005))
nrow(subset(Neg,Neg$sPETBH1<0.01 & Neg$sPETBH2<0.01))
nrow(subset(Neg,Neg$sPETBH1<0.05 & Neg$sPETBH2<0.05))

nrow(subset(Neg,Neg$sPTBA1<0.0005 & Neg$sPTBA2<0.0005))
nrow(subset(Neg,Neg$sPTBA1<0.005 & Neg$sPTBA2<0.005))
nrow(subset(Neg,Neg$sPTBA1<0.01 & Neg$sPTBA2<0.01))
nrow(subset(Neg,Neg$sPTBA1<0.05 & Neg$sPTBA2<0.05))

nrow(subset(Neg,Neg$ETBH1<0.0005 & Neg$sETBH1<0.0005 & miETBH1<0.0005 & sPETBH1<0.0005))
nrow(subset(Neg,Neg$ETBH1<0.005 & Neg$sETBH1<0.005 & miETBH1<0.005 & sPETBH1<0.005))
nrow(subset(Neg,Neg$ETBH1<0.01 & Neg$sETBH1<0.01 & miETBH1<0.01 & sPETBH1<0.01))
nrow(subset(Neg,Neg$ETBH1<0.05 & Neg$sETBH1<0.05 & miETBH1<0.05 & sPETBH1<0.05 ))

nrow(subset(Neg,Neg$ETBH2<0.0005 & Neg$sETBH2<0.0005 & miETBH2<0.0005 & sPETBH2<0.0005))
nrow(subset(Neg,Neg$ETBH2<0.005 & Neg$sETBH2<0.005 & miETBH2<0.005 & sPETBH2<0.005))
nrow(subset(Neg,Neg$ETBH2<0.01 & Neg$sETBH2<0.01 & miETBH2<0.01 & sPETBH2<0.01))
nrow(subset(Neg,Neg$ETBH2<0.05 & Neg$sETBH2<0.05 & miETBH2<0.05 & sPETBH2<0.05))
nrow(subset(Neg,Neg$ETBH1<0.05 & Neg$sETBH1<0.05 & miETBH1<0.05 & sPETBH1<0.05 & Neg$ETBH2<0.05 & Neg$sETBH2<0.05 & miETBH2<0.05 & sPETBH2<0.05))


nrow(subset(Neg,Neg$PTBH1<0.0005 & Neg$sPTBH1<0.0005 & miPTBH1<0.0005 & sPETBH1<0.0005 & sPTBA1<0.0005))
nrow(subset(Neg,Neg$PTBH1<0.005 & Neg$sPTBH1<0.005 & miPTBH1<0.005 & sPETBH1<0.005 & sPTBA1<0.005))
nrow(subset(Neg,Neg$PTBH1<0.01 & Neg$sPTBH1<0.01 & miPTBH1<0.01 & sPETBH1<0.01 & sPTBA1<0.01))
nrow(subset(Neg,Neg$PTBH1<0.05 & Neg$sPTBH1<0.05 & miPTBH1<0.05 & sPETBH1<0.05 & sPTBA1<0.05))

nrow(subset(Neg,Neg$PTBH2<0.0005 & Neg$sPTBH2<0.0005 & miPTBH2<0.0005 & sPETBH2<0.0005 & sPTBA2<0.0005))
nrow(subset(Neg,Neg$PTBH2<0.005 & Neg$sPTBH2<0.005 & miPTBH2<0.005 & sPETBH2<0.005 & sPTBA2<0.005))
nrow(subset(Neg,Neg$PTBH2<0.01 & Neg$sPTBH2<0.01 & miPTBH2<0.01 & sPETBH2<0.01 & sPTBA2<0.01))
nrow(subset(Neg,Neg$PTBH2<0.05 & Neg$sPTBH2<0.05 & miPTBH2<0.05 & sPETBH2<0.05 & sPTBA2<0.05))
nrow(subset(Neg,Neg$PTBH2<0.05 & Neg$sPTBH2<0.05 & miPTBH2<0.05 & sPETBH2<0.05 & sPTBA2<0.05 & Neg$PTBH1<0.05 & Neg$sPTBH1<0.05 & Neg$miPTBH1<0.05 & Neg$sPETBH1<0.05 & Neg$sPTBA1<0.05))


#####read in csv
setwd("/Users/tanyouran/OneDrive - Emory University/ECHO_AA_Cohort/Preterm_race/Updates_results_sPTB/POS/")
list.files(pattern="h.*.csv$") 
list.filenames <- list.files(pattern="h.*.csv$")
polname<-str_remove_all(list.filenames, ".csv")

for (i in 1:length(polname)) 
  assign(polname[i], read.csv(list.filenames[i]))

Pos<-as.data.frame(cbind(h_AllETBvsHealthyFullTerm_s1$P,h_AllETBvsHealthyFullTerm_s2$P,h_AllPTBvsHealthyFullTerm_s1$P,
                         h_AllPTBvsHealthyFullTerm_s2$P,h_InducedETBvsHealthyFullTerm_s1$P,h_InducedETBvsHealthyFullTerm_s2$P,
                         h_InducedPTBvsHealthyFullTerm_s1$P,h_InducedPTBvsHealthyFullTerm_s2$P,h_SpETBvsHealthyFullTerm_s1$P,
                         h_SpETBvsHealthyFullTerm_s2$P,h_SpPTBSpETBvsHealthyFullTerm_s1$P,h_SpPTBSpETBvsHealthyFullTerm_s2$P,
                         h_SpPTBvsALLFullTerm_s1$P,h_SpPTBvsALLFullTerm_s2$P,h_SpPTBvsHealthyFullTerm_s1$P,
                         h_SpPTBvsHealthyFullTerm_s2$P))
colnames(Pos)<-c("ETBH1","ETBH2","PTBH1","PTBH2","miETBH1","miETBH2","miPTBH1","miPTBH2","sETBH1","sETBH2",
                 "sPETBH1","sPETBH2","sPTBA1","sPTBA2","sPTBH1","sPTBH2")

nrow(subset(Pos,Pos$ETBH1<0.0005 & Pos$ETBH2<0.0005))
nrow(subset(Pos,Pos$ETBH1<0.005 & Pos$ETBH2<0.005))
nrow(subset(Pos,Pos$ETBH1<0.01 & Pos$ETBH2<0.01))
nrow(subset(Pos,Pos$ETBH1<0.05 & Pos$ETBH2<0.05))

nrow(subset(Pos,Pos$sETBH1<0.0005 & Pos$sETBH2<0.0005))
nrow(subset(Pos,Pos$sETBH1<0.005 & Pos$sETBH2<0.005))
nrow(subset(Pos,Pos$sETBH1<0.01 & Pos$sETBH2<0.01))
nrow(subset(Pos,Pos$sETBH1<0.05 & Pos$sETBH2<0.05))

nrow(subset(Pos,Pos$miETBH1<0.0005 & Pos$miETBH2<0.0005))
nrow(subset(Pos,Pos$miETBH1<0.005 & Pos$miETBH2<0.005))
nrow(subset(Pos,Pos$miETBH1<0.01 & Pos$miETBH2<0.01))
nrow(subset(Pos,Pos$miETBH1<0.05 & Pos$miETBH2<0.05))

nrow(subset(Pos,Pos$PTBH1<0.0005 & Pos$PTBH2<0.0005))
nrow(subset(Pos,Pos$PTBH1<0.005 & Pos$PTBH2<0.005))
nrow(subset(Pos,Pos$PTBH1<0.01 & Pos$PTBH2<0.01))
nrow(subset(Pos,Pos$PTBH1<0.05 & Pos$PTBH2<0.05))

nrow(subset(Pos,Pos$sPTBH1<0.0005 & Pos$sPTBH2<0.0005))
nrow(subset(Pos,Pos$sPTBH1<0.005 & Pos$sPTBH2<0.005))
nrow(subset(Pos,Pos$sPTBH1<0.01 & Pos$sPTBH2<0.01))
nrow(subset(Pos,Pos$sPTBH1<0.05 & Pos$sPTBH2<0.05))

nrow(subset(Pos,Pos$miPTBH1<0.0005 & Pos$miPTBH2<0.0005))
nrow(subset(Pos,Pos$miPTBH1<0.005 & Pos$miPTBH2<0.005))
nrow(subset(Pos,Pos$miPTBH1<0.01 & Pos$miPTBH2<0.01))
nrow(subset(Pos,Pos$miPTBH1<0.05 & Pos$miPTBH2<0.05))

nrow(subset(Pos,Pos$sPTBA1<0.0005 & Pos$sPTBA2<0.0005))
nrow(subset(Pos,Pos$sPTBA1<0.005 & Pos$sPTBA2<0.005))
nrow(subset(Pos,Pos$sPTBA1<0.01 & Pos$sPTBA2<0.01))
nrow(subset(Pos,Pos$sPTBA1<0.05 & Pos$sPTBA2<0.05))

nrow(subset(Pos,Pos$sPETBH1<0.0005 & Pos$sPETBH2<0.0005))
nrow(subset(Pos,Pos$sPETBH1<0.005 & Pos$sPETBH2<0.005))
nrow(subset(Pos,Pos$sPETBH1<0.01 & Pos$sPETBH2<0.01))
nrow(subset(Pos,Pos$sPETBH1<0.05 & Pos$sPETBH2<0.05))

nrow(subset(Pos,Pos$ETBH1<0.0005 & Pos$sETBH1<0.0005 & miETBH1<0.0005 & sPETBH1<0.0005))
nrow(subset(Pos,Pos$ETBH1<0.005 & Pos$sETBH1<0.005 & miETBH1<0.005 & sPETBH1<0.005))
nrow(subset(Pos,Pos$ETBH1<0.01 & Pos$sETBH1<0.01 & miETBH1<0.01 & sPETBH1<0.01))
nrow(subset(Pos,Pos$ETBH1<0.05 & Pos$sETBH1<0.05 & miETBH1<0.05 & sPETBH1<0.05 ))

nrow(subset(Pos,Pos$ETBH2<0.0005 & Pos$sETBH2<0.0005 & miETBH2<0.0005 & sPETBH2<0.0005))
nrow(subset(Pos,Pos$ETBH2<0.005 & Pos$sETBH2<0.005 & miETBH2<0.005 & sPETBH2<0.005))
nrow(subset(Pos,Pos$ETBH2<0.01 & Pos$sETBH2<0.01 & miETBH2<0.01 & sPETBH2<0.01))
nrow(subset(Pos,Pos$ETBH2<0.05 & Pos$sETBH2<0.05 & miETBH2<0.05 & sPETBH2<0.05))
nrow(subset(Pos,Pos$ETBH1<0.05 & Pos$sETBH1<0.05 & miETBH1<0.05 & sPETBH1<0.05 & Pos$ETBH2<0.05 & Pos$sETBH2<0.05 & miETBH2<0.05 & sPETBH2<0.05))


nrow(subset(Pos,Pos$PTBH1<0.0005 & Pos$sPTBH1<0.0005 & miPTBH1<0.0005 & sPETBH1<0.0005 & sPTBA1<0.0005))
nrow(subset(Pos,Pos$PTBH1<0.005 & Pos$sPTBH1<0.005 & miPTBH1<0.005 & sPETBH1<0.005 & sPTBA1<0.005))
nrow(subset(Pos,Pos$PTBH1<0.01 & Pos$sPTBH1<0.01 & miPTBH1<0.01 & sPETBH1<0.01 & sPTBA1<0.01))
nrow(subset(Pos,Pos$PTBH1<0.05 & Pos$sPTBH1<0.05 & miPTBH1<0.05 & sPETBH1<0.05 & sPTBA1<0.05))

nrow(subset(Pos,Pos$PTBH2<0.0005 & Pos$sPTBH2<0.0005 & miPTBH2<0.0005 & sPETBH2<0.0005 & sPTBA2<0.0005))
nrow(subset(Pos,Pos$PTBH2<0.005 & Pos$sPTBH2<0.005 & miPTBH2<0.005 & sPETBH2<0.005 & sPTBA2<0.005))
nrow(subset(Pos,Pos$PTBH2<0.01 & Pos$sPTBH2<0.01 & miPTBH2<0.01 & sPETBH2<0.01 & sPTBA2<0.01))
nrow(subset(Pos,Pos$PTBH2<0.05 & Pos$sPTBH2<0.05 & miPTBH2<0.05 & sPETBH2<0.05 & sPTBA2<0.05 & Pos$PTBH1<0.05 & Pos$sPTBH1<0.05 & miPTBH1<0.05 & sPETBH1<0.05 & sPTBA1<0.05))


Neg<-as.data.frame(cbind(c_AllETBvsHealthyFullTerm_s1$pfdr,c_AllETBvsHealthyFullTerm_s2$pfdr,c_AllPTBvsHealthyFullTerm_s1$pfdr,
                         c_AllPTBvsHealthyFullTerm_s2$pfdr,c_InducedETBvsHealthyFullTerm_s1$pfdr,c_InducedETBvsHealthyFullTerm_s2$pfdr,
                         c_InducedPTBvsHealthyFullTerm_s1$pfdr,c_InducedPTBvsHealthyFullTerm_s2$pfdr,c_SpETBvsHealthyFullTerm_s1$pfdr,
                         c_SpETBvsHealthyFullTerm_s2$pfdr,c_SpPTBSpETBvsHealthyFullTerm_s1$pfdr,c_SpPTBSpETBvsHealthyFullTerm_s2$pfdr,
                         c_SpPTBvsALLFullTerm_s1$pfdr,c_SpPTBvsALLFullTerm_s2$pfdr,c_SpPTBvsHealthyFullTerm_s1$pfdr,
                         c_SpPTBvsHealthyFullTerm_s2$pfdr))
colnames(Neg)<-c("ETBH1","ETBH2","PTBH1","PTBH2","miETBH1","miETBH2","miPTBH1","miPTBH2","sETBH1","sETBH2",
                 "sPETBH1","sPETBH2","sPTBA1","sPTBA2","sPTBH1","sPTBH2")

nrow(subset(Neg,Neg$ETBH1<0.2 & Neg$ETBH2<0.2))
nrow(subset(Neg,Neg$sETBH1<0.2 & Neg$sETBH2<0.2))
nrow(subset(Neg,Neg$miETBH1<0.2 & Neg$miETBH2<0.2))
nrow(subset(Neg,Neg$PTBH1<0.2 & Neg$PTBH2<0.2))
nrow(subset(Neg,Neg$sPTBH1<0.2 & Neg$sPTBH2<0.2))
nrow(subset(Neg,Neg$miPTBH1<0.2 & Neg$miPTBH2<0.2))
nrow(subset(Neg,Neg$ETBH1<0.2 & Neg$sETBH1<0.2 & miETBH1<0.2 & sPETBH1<0.2))
nrow(subset(Neg,Neg$ETBH2<0.2 & Neg$sETBH2<0.2 & miETBH2<0.2 & sPETBH2<0.2))


Pos<-as.data.frame(cbind(h_AllETBvsHealthyFullTerm_s1$pfdr,h_AllETBvsHealthyFullTerm_s2$pfdr,h_AllPTBvsHealthyFullTerm_s1$pfdr,
                         h_AllPTBvsHealthyFullTerm_s2$pfdr,h_InducedETBvsHealthyFullTerm_s1$pfdr,h_InducedETBvsHealthyFullTerm_s2$pfdr,
                         h_InducedPTBvsHealthyFullTerm_s1$pfdr,h_InducedPTBvsHealthyFullTerm_s2$pfdr,h_SpETBvsHealthyFullTerm_s1$pfdr,
                         h_SpETBvsHealthyFullTerm_s2$pfdr,h_SpPTBSpETBvsHealthyFullTerm_s1$pfdr,h_SpPTBSpETBvsHealthyFullTerm_s2$pfdr,
                         h_SpPTBvsALLFullTerm_s1$pfdr,h_SpPTBvsALLFullTerm_s2$pfdr,h_SpPTBvsHealthyFullTerm_s1$pfdr,
                         h_SpPTBvsHealthyFullTerm_s2$pfdr))
colnames(Pos)<-c("ETBH1","ETBH2","PTBH1","PTBH2","miETBH1","miETBH2","miPTBH1","miPTBH2","sETBH1","sETBH2",
                 "sPETBH1","sPETBH2","sPTBA1","sPTBA2","sPTBH1","sPTBH2")

nrow(subset(Pos,Pos$ETBH1<0.2 & Pos$ETBH2<0.2))
nrow(subset(Pos,Pos$sETBH1<0.2 & Pos$sETBH2<0.2))
nrow(subset(Pos,Pos$miETBH1<0.2 & Pos$miETBH2<0.2))
nrow(subset(Pos,Pos$miPTBH1<0.2 & Pos$miPTBH2<0.2))
nrow(subset(Pos,Pos$PTBH1<0.2 & Pos$PTBH2<0.2))
nrow(subset(Pos,Pos$sPTBH1<0.2 & Pos$sPTBH2<0.2))
nrow(subset(Pos,Pos$PETBH1<0.2 & Pos$PETBH2<0.2))
nrow(subset(Pos,Pos$ETBH2<0.2 & Pos$sETBH2<0.2 & miETBH2<0.2 & sPETBH2<0.2))





