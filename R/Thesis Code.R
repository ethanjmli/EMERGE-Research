####Packages####
#remotes::install_github("yufree/xMSanalyzer")
pacman::p_load(remotes,
               rgl,
               apLCMS,
               GO.db,
               xcms,
               dplyr,
               readr,
               readxl,
               here,
               tableone,
               gam,
               hydroTSM)
library(xMSanalyzer)

####Get Detected Features####
#####Get HILIC Features####
hilic_results <- read_delim(here::here("Data","2020","ComBat_mzcalibrated_untargeted_mediansummarized_featuretable_hilic.txt"), 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

#####Get C18 Features####
c18_results <- read_delim(here::here("Data","2020","ComBat_mzcalibrated_untargeted_mediansummarized_featuretable_c18.txt"), 
                          "\t", escape_double = FALSE, trim_ws = TRUE)


####Get Sample IDs####
#####Get HILIC sample ids####
hilic_ids <- read.delim(here::here("Data","2020","R24_mapping_hilicpos.txt"))

#####Get C18 Sample IDs####
c18_ids <- read.delim(here::here("Data","2020","R24_mapping_c18neg.txt"))


####Get Study Demographics w/Temp####
demo <- read.csv(here::here("Data","MPTB_AA Cohort Data_N547_addZscore.csv"))
demographicsbystudyid<-read.csv(here::here("Data","AA Cohort Demographics + Temperature + Humidity.csv"))


colnames(demographicsbystudyid)
demographicsbystudyid <- demographicsbystudyid%>%
  select(-full_address)

####2020 HILIC Analysis####
#####Merging HILIC sample IDs and demographics with Feature table#####
hilic_ids$seq<-(1:nrow(hilic_ids))

#Extract subject ids from hilic sample IDs
hilic_ids$subjectid<-substr(hilic_ids$Sample.ID,1,5)

#Extract study stage (visit #) from hilic sample IDs
hilic_ids$stage<-substr(hilic_ids$Sample.ID,7,7)

#Merge hilic sample IDs and demographics on subject id
demographicsbystudyid<-demographicsbystudyid[order(demographicsbystudyid$subjectid),]

demographics.IDs.total <- hilic_ids %>%
  merge(demographicsbystudyid,by="subjectid")

#Merge demographics with HILIC results on filename (filename comes from HILIC sample IDs)
HRt<-as.data.frame(t(hilic_results))
rownames(HRt)
HRt$File.Name<-rownames(HRt) 
HRt$File.Name 
table(HRt$File.Name %in% rownames(HRt))
table(HRt$File.Name %in% colnames(hilic_results))

table(demographics.IDs.total$File.Name %in% HRt$File.Name)

Metotal<- demographics.IDs.total %>%
  merge(HRt, by="File.Name")

table(Metotal$stage,useNA="always")   
# -    1    2    P <NA> 
# 1  330  270   49    0

sum(is.na(Metotal$temp_exposed_status)) #0

#####Create Gestational Age Variable#####
Metotal$gast<-NA
for (i in 1:nrow(Metotal)){
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


#####Create nulliparity variable#####
summary(Metotal$parity)
#Initialize variable nullp, initial value = NA for all observations
#1 = not nulliparous, 0 = nulliparous
Metotal$nullp<-NA
for (i in 1:nrow(Metotal)){
  if (Metotal$parity[i]==0)
  {
    Metotal$nullp[i]<-0
  }
  else if (Metotal$parity[i]>0)
  {
    Metotal$nullp[i]<-1
  }
}

Metotal <- Metotal %>%
  mutate(nullp = factor(nullp, levels = c(0,1)))

table(Metotal$nullp) #Yes Nulliparous Not Nulliparous 
                     #283             367 

#####Clean Covariate Data#####
colnames(Metotal[185:195]) 
colnames(Metotal[196:206]) 
colnames(Metotal[201:205]) #Metabolites begin at Column 201
colnames(Metotal[13815:13818]) #Metabolites end at Column 13816
colnames(Metotal[13816:13818])

#Extract covariates (everything except metabolite intensity)
Metotal_no_metabolites <-Metotal %>%  #[,c(1:190,"gast","nullp")] #metotaldemo
  dplyr::select(1:200,gast,nullp)


#Clean covariates
table(Metotal_no_metabolites$stage, useNA= "always")    #-   1   2   P <NA>
                                                        #1 330 270  49    0

Metotal_no_metabolites <- Metotal_no_metabolites %>%
  filter(stage ==1 | stage == 2)

table(Metotal_no_metabolites$stage, useNA = "always") #   1    2 <NA> 
                                                      # 330  270    0
nrow(Metotal_no_metabolites) #600
# Metotal_no_metabolites_stage_1 <-subset(Metotal_no_metabolites,stage==1)   
# nrow(Metotal_no_metabolites_stage_1) #330

#Check # of observations for each complication
##GHTN
table(Metotal_no_metabolites$Comp_gHTNvsHealthyFullTerm, useNA= "always") 
# 0    1 <NA> 
# 302   77  221 

##GDM
table(Metotal_no_metabolites$Comp_gDMvsHealthyFullTerm, useNA= "always")
# 0    1 <NA> 
# 302   23  275  

##IUGR
table(Metotal_no_metabolites$Comp_IUGRvsHealthyFullTerm, useNA= "always")
#   0    1 <NA> 
# 302   37  261

##PROM
table(Metotal_no_metabolites$Comp_PROMvsHealthyFullTerm, useNA= "always")
#   0    1 <NA> 
# 302   24  274


#Exclude spontaneous/elective (stillbirth) abortion outcomes
table(Metotal_no_metabolites$birthoutcome, useNA= "always") # EarlyTerm    FullTerm     Preterm Preterm-STILLBIRTH      <NA> 
                                                            #        97         395         107                  1         0

Metotal_no_metabolites_noSTB <-Metotal_no_metabolites %>%
  dplyr::filter(!(birthoutcome=="Preterm-STILLBIRTH" ))

nrow(Metotal_no_metabolites_noSTB)
#[1] 599

table(Metotal_no_metabolites_noSTB$birthoutcome)  
#EarlyTerm  FullTerm   Preterm 
#97       395       107 

table(Metotal_no_metabolites_noSTB$stage)
#   1   2 
# 329 270

#Reexamine complications 
##GHTN
table(Metotal_no_metabolites_noSTB$Comp_gHTNvsHealthyFullTerm, useNA= "always") 
#0    1 <NA> 
#  302   77  220

##GDM
table(Metotal_no_metabolites_noSTB$Comp_gDMvsHealthyFullTerm, useNA= "always")
#0    1 <NA> 
#  302   23  274 

##IUGR
table(Metotal_no_metabolites_noSTB$Comp_IUGRvsHealthyFullTerm, useNA= "always")
#0    1 <NA> 
#  302   37  260 

##PROM
table(Metotal_no_metabolites_noSTB$Comp_PROMvsHealthyFullTerm, useNA= "always")
#0    1 <NA> 
#  302   24  273


#####Removing Duplicates#####
table(Metotal_no_metabolites_noSTB$subjectid)
table(table(Metotal_no_metabolites_noSTB$subjectid))   #1   2   3 
                                                       #64 266   1 
                                                       #Only 2 vists max, no subject id should appear >2 times

which(table(Metotal_no_metabolites_noSTB$subjectid)>2) #subject id E0586 appears more than 2 times (3 times)

subjectids <- data.frame(table(Metotal_no_metabolites_noSTB$subjectid)) %>% #confirm freq of E0586 = 3
  rename(subjectid = Var1)                                    
subjectids[which(subjectids$Freq>2),] #           subjectid Freq
                                      #       50      E0586    3

Metotal_no_metabolites_noSTB[Metotal_no_metabolites_noSTB$subjectid=="E0586",c("File.Name","Sample.ID","seq")]
Metotal_no_metabolites_noSTB[388,] == Metotal_no_metabolites_noSTB[389,]
Metotal_no_metabolites_noSTB[388,] == Metotal_no_metabolites_noSTB[402,]
Metotal_no_metabolites_noSTB_nodupes <- Metotal_no_metabolites_noSTB %>%
  filter(!(File.Name == "F069_200209_M429_175"))

sum(duplicated(Metotal_no_metabolites_noSTB[,-c(1,3,6)])) #1
sum(duplicated(Metotal_no_metabolites_noSTB_nodupes[,-c(1,3,6)])) #0

table(Metotal_no_metabolites_noSTB_nodupes$subjectid)
which(table(Metotal_no_metabolites_noSTB_nodupes$subjectid)>2) #None appear > 2 times

subjectids_nodupes <- data.frame(table(Metotal_no_metabolites_noSTB_nodupes$subjectid)) %>% 
  rename(subjectid = Var1)                                    
subjectids[which(subjectids_nodupes$Freq>2),] #0

#####Check Number of Unique Subjects#####
subjectids_1_time <-rownames(subset(table(Metotal_no_metabolites_noSTB_nodupes$subjectid),
                                    table(Metotal_no_metabolites_noSTB_nodupes$subjectid)==1))  
subjectids_2_times <- rownames(subset(table(Metotal_no_metabolites_noSTB_nodupes$subjectid),
                                     table(Metotal_no_metabolites_noSTB_nodupes$subjectid)!=1))

#Assign to metotalv1, the observations of metotalnmis whose subjectids are in rn1
metotalv1<-Metotal_no_metabolites_noSTB_nodupes[which(Metotal_no_metabolites_noSTB_nodupes$subjectid %in% subjectids_1_time),] ######single stage 
length(table(metotalv1$subjectid))    # 64 obs- 64 subjectids appear only once

#Assign to metotalv2, the observations of metotalnmis whose subjectids are in rn2
metotalv2<-Metotal_no_metabolites_noSTB_nodupes[which(Metotal_no_metabolites_noSTB_nodupes$subjectid %in% subjectids_2_times),] ###### both stages
length(table(metotalv2$subjectid))   #267 subjectids appear twice
# 534 obs- 267 subjects total, 2 visits each: 267+64=331 subjects
# 328 subjects visit1-270 visit2
#267 subjects each had 2 observations = 534 observations in metotalv2, + 64 subjects w/ only 1 observation in metotalv1 331 total subjects

table(metotalv1$stage)
metotalv1[which(metotalv1$stage == 2),"subjectid"] 
#331 total subjects; however, 3 only appeared once, in stage 2; therefore, 328 total stage 1 subjects/observations
#E0539, E0545, E0637
length(unique(Metotal_no_metabolites_noSTB_nodupes$subjectid))
#Confirm, 331 total

unique_subjects <- rbind(metotalv1, metotalv2[metotalv2$stage==1,]) #331 subjects
unique_subjects_analytic <- unique_subjects %>% #328 subjects; exclude 3 subjects w/only visit 2
  filter(stage == 1)

unique_subjects_analytic$sub_use <- NA
unique_subjects_analytic$sub_use <- ifelse(unique_subjects_analytic$TobaccoUse_MR ==1 | 
                                             unique_subjects_analytic$TobaccoUse_MRorSR ==1 | 
                                             unique_subjects_analytic$tobacco.AnyInPreg ==1 |
                                             unique_subjects_analytic$AlcoholUse_MR ==1 | 
                                             unique_subjects_analytic$AlcoholUse_MRorSR ==1 | 
                                             unique_subjects_analytic$alcohol.AnyInPreg ==1  |
                                             unique_subjects_analytic$MarijuanaUse_MR ==1 | 
                                             unique_subjects_analytic$MarijuanaUse_MRorSR ==1 | 
                                             unique_subjects_analytic$marijuana.AnyInPreg, 1, 0)
table(unique_subjects_analytic$sub_use)  #  0   1 
                                         #169 152

#####Create Table 1######
names(unique_subjects_analytic)

#Question: why only these variables?
#List of continuous variables in a
contVars<-c("age_enrollment","birthga","FirstPrenatalBMI","temperature","humidity")
#List of categorical variables in  a
catVars<-c("temp_exposed_status","Education_4.level","income",
           #my addition#
           "PresInsurance_3.level", 
           #my addition end#
           "nullp","priorabortion","Site_Prenatal","Depression_Preg", "MarriedCohab_Not","priorpreterm",
           "Sex","g_htn","g_dm","PROM","Comp_gHTNvsHealthyFullTerm","Comp_gDMvsHealthyFullTerm","Comp_PROMvsHealthyFullTerm",
           "Comp_IUGRvsHealthyFullTerm","mode","birthoutcome","TobaccoUse_MR","AlcoholUse_MR","MarijuanaUse_MR","sub_use")

##Create Table One for continuous variables from a listed in contVars
tab_Continuous<-CreateTableOne(vars=contVars, data = unique_subjects_analytic)
print(tab_Continuous, showAllLevels = TRUE, test = FALSE)

#                             level  Overall     
# n                                    328       
# age_enrollment (mean (SD))         24.93 (4.75)
# birthga (mean (SD))                38.20 (3.01)
# FirstPrenatalBMI (mean (SD))       28.60 (7.65)
# temperature (mean (SD))            26.04 (7.05)
# humidity (mean (SD))                0.01 (0.00)

#Create Table One for categorical variables from a listed in catVars
tab_Categorical<-CreateCatTable(vars=catVars, data= unique_subjects_analytic, includeNA = TRUE)
print(tab_Categorical, catDigits = 3,pDigits = 3,showAllLevels = TRUE, test = FALSE)

#                                level           Overall    
# n                                              328        
# temp_exposed_status (%)        0               293 (89.3) 
#                                1                35 (10.7) 
# Education_4.level (%)          1                54 (16.5) 
#                                2               123 (37.5) 
#                                3               100 (30.5) 
#                                4                51 (15.5) 
# income (%)                     0               152 (46.3) 
#                                1                41 (12.5) 
#                                2                32 ( 9.8) 
#                                3                38 (11.6) 
#                                4                21 ( 6.4) 
#                                5                19 ( 5.8) 
#                                6                25 ( 7.6) 
# PresInsurance_3.level (%)      1               116 (35.4) 
#                                2               145 (44.2) 
#                                3                67 (20.4) 
# nullp (%)                      Yes Nulliparous 147 (44.8) 
#                                Not Nulliparous 181 (55.2) 
# priorabortion (%)              0                 3 ( 0.9) 
#                                1               325 (99.1) 
# Site_Prenatal (%)              Emory           127 (38.7) 
#                                Grady           201 (61.3) 
# Depression_Preg (%)            0               314 (95.7) 
#                                1                12 ( 3.7) 
#                                1.1               2 ( 0.6) 
# MarriedCohab_Not (%)           0               157 (47.9) 
#                                1               171 (52.1) 
# priorpreterm (%)               0               288 (87.8) 
#                                1                29 ( 8.8) 
#                                2                 7 ( 2.1) 
#                                3                 4 ( 1.2) 
# Sex (%)                        1               163 (49.7) 
#                                2               165 (50.3) 
# g_htn (%)                      0               260 (79.3) 
#                                1                42 (12.8) 
#                                2                26 ( 7.9) 
# g_dm (%)                       0               315 (96.0) 
#                                1                13 ( 4.0) 
# PROM (%)                       0               314 (95.7) 
#                                1                14 ( 4.3) 
# Comp_gHTNvsHealthyFullTerm (%) 0               159 (48.5) 
#                                1                42 (12.8) 
#                                <NA>            127 (38.7) 
# Comp_gDMvsHealthyFullTerm (%)  0               159 (48.5) 
#                                1                13 ( 4.0) 
#                                <NA>            156 (47.6) 
# Comp_PROMvsHealthyFullTerm (%) 0               159 (48.5) 
#                                1                14 ( 4.3) 
#                                <NA>            155 (47.3) 
# Comp_IUGRvsHealthyFullTerm (%) 0               159 (48.5) 
#                                1                22 ( 6.7) 
#                                <NA>            147 (44.8) 
# mode (%)                       C-section        75 (22.9) 
#                                Vaginal         253 (77.1) 
# birthoutcome (%)               EarlyTerm        54 (16.5) 
#                                FullTerm        210 (64.0) 
#                                Preterm          64 (19.5) 
# TobaccoUse_MR (%)              0               278 (84.8) 
#                                1                50 (15.2) 
# AlcoholUse_MR (%)              0               302 (92.1) 
#                                1                26 ( 7.9) 
# MarijuanaUse_MR (%)            0               212 (64.6) 
#                                1               116 (35.4) 
# sub_use (%)                    0               169 (51.5) 
#                                1               152 (46.3) 
#                                <NA>              7 ( 2.1)


unique_subjects_analytic$subjectid %in% Metotal_no_metabolites_noSTB_nodupes[which(Metotal_no_metabolites_noSTB_nodupes$stage==1),]$subjectid
Metotal_no_metabolites_noSTB_nodupes[which(Metotal_no_metabolites_noSTB_nodupes$stage==1),]$subjectid %in% unique_subjects_analytic$subjectid

#####Subset to stage 1 and recheck Table 1#####
Metotal_no_metabolites_noSTB_nodupes_s1 <- Metotal_no_metabolites_noSTB_nodupes %>%
  filter(stage == 1)

contVars <- c(contVars, "gast")

##Create Table One for continuous variables from a listed in contVars
tab_Continuous<-CreateTableOne(vars=contVars, data = Metotal_no_metabolites_noSTB_nodupes_s1)
print(tab_Continuous, showAllLevels = TRUE, test = FALSE)
#gast (mean (SD))                   11.45 (2.20)

#Create Table One for categorical variables from a listed in catVars
tab_Categorical<-CreateCatTable(vars=catVars, data= Metotal_no_metabolites_noSTB_nodupes_s1, includeNA = TRUE)
print(tab_Categorical, catDigits = 3,pDigits = 3,showAllLevels = TRUE, test = FALSE)


#####Log transform metabolite intensities####
Metotal_metabolites <-Metotal %>%
  filter(stage == 1,
         birthoutcome=="EarlyTerm" | birthoutcome=="FullTerm" |birthoutcome=="Preterm",
         !(File.Name == "F069_200209_M429_175"))
Metotal_no_metabolites_noSTB_nodupes_s1$subjectid %in% Metotal_metabolites$subjectid
Metotal_metabolites$subjectid %in% Metotal_no_metabolites_noSTB_nodupes_s1$subjectid

Metotal_clean_transformed <- Metotal_no_metabolites_noSTB_nodupes_s1 %>%
  left_join(Metotal_metabolites[,c(2,201:13816)], by = join_by("subjectid"))

Metotal_clean_transformed[203:13818] <- log2(Metotal_clean_transformed[203:13818])
Metotal_clean_transformed[Metotal_clean_transformed == -Inf] <- NA
saveRDS(Metotal_clean_transformed, here::here("Data","Thesis","HILIC 2020 Analytic Dataset.rds"))


#####2020 HILIC Temperature vs Metabolite Regressions####

Metotal_clean_transformed <- readRDS(here::here("Data","Thesis","HILIC 2020 Analytic Dataset.rds")) %>%
  mutate(temp_exposed_status = factor(temp_exposed_status, levels = c(0,1),labels = c("Unexposed", "Exposed")),
         Sex = factor(Sex, levels = c(1,2), labels = c("Male","Female")),
         PresInsurance_3.level = factor(PresInsurance_3.level, levels = c(1,2,3),labels = c("Low-Income Medicaid","Right from the Start Medicaid","Private")),
         nullp = factor(nullp,levels = c(0,1), labels = c("Yes Nulliparous","Not Nulliparous")),
         AlcoholUse_MRorSR = factor(AlcoholUse_MRorSR, levels = c(0,1), labels = c("No","Yes")),
         TobaccoUse_MRorSR = factor(TobaccoUse_MRorSR, levels = c(0,1),labels = c("No","Yes")),
         conception_season = factor(conception_season)
         )

exposure_threshold_0.95 <- quantile(Metotal_clean_transformed$temperature,probs = (0.95))
exposure_threshold_0.75 <- quantile(Metotal_clean_transformed$temperature,probs = (0.75))

Metotal_clean_transformed<-Metotal_clean_transformed%>%
  mutate(temp_exposed_status_0.95 = case_when(temperature >= exposure_threshold_0.95 ~ 1,
                                              temperature < exposure_threshold_0.95 ~ 0),
         temp_exposed_status_0.75 = case_when(temperature >= exposure_threshold_0.75~1,
                                              temperature < exposure_threshold_0.75 ~0),
         Alc_and_Tobacco = case_when(AlcoholUse_MRorSR == "Yes" ~ 1,
                                     TobaccoUse_MRorSR == "Yes" ~ 1,
                                     TRUE ~ 0))

table(Metotal_clean_transformed$Alc_and_Tobacco, Metotal_clean_transformed$AlcoholUse_MRorSR)
table(Metotal_clean_transformed$Alc_and_Tobacco, Metotal_clean_transformed$TobaccoUse_MRorSR)

table(Metotal_clean_transformed$temp_exposed_status)
plot(Metotal_clean_transformed$gast, Metotal_clean_transformed$temperature)
plot(Metotal_clean_transformed$temperature,Metotal_clean_transformed$V120)
plot(Metotal_clean_transformed$temperature, Metotal_clean_transformed$Comp_gHTNvsHealthyFullTerm)
plot(Metotal_clean_transformed$temperature, Metotal_clean_transformed$Comp_gDMvsHealthyFullTerm)
plot(Metotal_clean_transformed$temperature, Metotal_clean_transformed$Comp_IUGRvsHealthyFullTerm)
plot(Metotal_clean_transformed$temperature, Metotal_clean_transformed$Comp_PROMvsHealthyFullTerm)
gam1<-gam(temperature~s(gast,df=327),data = Metotal_clean_transformed)
plot(gam1)
abline(lm(V120 ~ temperature, data = Metotal_clean_transformed))
summary(lm(V12242 ~ temperature, data = Metotal_clean_transformed))

is.factor(Metotal_clean_transformed$conception_season)

colnames(Metotal_clean_transformed[203:207]) #Metabolites begin at Column 203
colnames(Metotal_clean_transformed[13815:13818]) #Metabolites end at Column 13818

#####Categorical Temperature####
HILIC_cat_temp_vs_metabolite <-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(HILIC_cat_temp_vs_metabolite)<-c("Est","Std","t-value","P")

i<-1
for (i in 1:nrow(hilic_results) ) {
  tryCatch(
    {
      lmfit<-lm(Metotal_clean_transformed[,i+202] ~temp_exposed_status+age_enrollment+FirstPrenatalBMI+Sex+
                  PresInsurance_3.level+nullp + AlcoholUse_MRorSR + TobaccoUse_MRorSR + gast + conception_season +
                  humidity, 
                data=Metotal_clean_transformed)
      HILIC_cat_temp_vs_metabolite[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}

HILIC_cat_temp_vs_metabolite$pfdr<-p.adjust(HILIC_cat_temp_vs_metabolite$P,method="BH")
write.csv(HILIC_cat_temp_vs_metabolite,here::here("Data","Thesis","HILIC Categorical Temperature Regressions.csv"))


#####Continuous Temperature####
HILIC_cont_temp_vs_metabolite <-data.frame(matrix(nrow = nrow(hilic_results), ncol = 4))
colnames(HILIC_cont_temp_vs_metabolite)<-c("Est","Std","t-value","P")

i<-1
for (i in 1:nrow(hilic_results) ) {
  tryCatch(
    {
      lmfit<-lm(Metotal_clean_transformed[,i+202] ~temperature+age_enrollment+FirstPrenatalBMI+Sex+
                  PresInsurance_3.level+nullp + AlcoholUse_MRorSR + TobaccoUse_MRorSR + gast + conception_season +
                  humidity, 
                data=Metotal_clean_transformed)
      HILIC_cont_temp_vs_metabolite[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}

HILIC_cont_temp_vs_metabolite$pfdr<-p.adjust(HILIC_cont_temp_vs_metabolite$P,method="BH")
write.csv(HILIC_cont_temp_vs_metabolite,here::here("Data","Thesis","HILIC Continuous Temperature Regressions.csv"))

####2020 C18 Analysis####
#####Merging C18 sample IDs and demographics with Feature table#####
c18_ids$seq<-(1:nrow(c18_ids))

#Extract subject ids from c18 sample IDs
c18_ids$subjectid<-substr(c18_ids$Sample.ID,1,5)

#Extract study stage (visit #) from c18 sample IDs
c18_ids$stage<-substr(c18_ids$Sample.ID,7,7)

#Merge c18 sample IDs and demographics on subject id
demographicsbystudyid<-demographicsbystudyid[order(demographicsbystudyid$subjectid),]

demographics.IDs.total_2 <- c18_ids %>%
  merge(demographicsbystudyid,by="subjectid")

#Merge demographics with C18 results on filename (filename comes from C18 sample IDs)
CRt<-as.data.frame(t(c18_results))
rownames(CRt)
CRt$File.Name<-rownames(CRt) 
CRt$File.Name 
table(CRt$File.Name %in% rownames(CRt))
table(CRt$File.Name %in% colnames(c18_results))

table(demographics.IDs.total_2$File.Name %in% HRt$File.Name)

Metotal_2<- demographics.IDs.total_2 %>%
  merge(CRt, by="File.Name")

table(Metotal_2$stage)   

sum(is.na(Metotal_2$temp_exposed_status))

#####Create Gestational Age Variable#####
Metotal_2$gast<-NA
for (i in 1:nrow(Metotal_2)){
  if (Metotal_2$stage[i]=='1')
  {
    Metotal_2$gast[i]<-Metotal_2$gasamp1_wk.10[i]
  }
  else if (Metotal_2$stage[i]=='2')
  {
    Metotal_2$gast[i]<-Metotal_2$gasamp2_wk.10[i]
  }
}
summary(Metotal_2$gasamp1_wk.10)
summary(Metotal_2$gasamp2_wk.10)
summary(Metotal_2$gast)


#####Create nulliparity variable#####
summary(Metotal_2$parity)
#Initialize variable nullp, initial value = NA for all observations
#1 = not nulliparous, 0 = nulliparous
Metotal_2$nullp<-NA
for (i in 1:nrow(Metotal_2)){
  if (Metotal_2$parity[i]==0)
  {
    Metotal_2$nullp[i]<-0
  }
  else if (Metotal_2$parity[i]>0)
  {
    Metotal_2$nullp[i]<-1
  }
}

Metotal_2 <- Metotal_2 %>%
  mutate(nullp = factor(nullp, levels = c(0,1)))

table(Metotal_2$nullp) #0   1 
#284 366
#Discrepancy on row 389

#
Metotal$nullp == Metotal_2$nullp
Metotal$subjectid == Metotal_2$subjectid
Metotal$File.Name == Metotal_2$File.Name

#####Clean Covariate Data#####
colnames(Metotal_2[185:195]) 
colnames(Metotal_2[196:206]) 
colnames(Metotal_2[201:205]) #Metabolites begin at Column 201
colnames(Metotal_2[12100:12102]) #Metabolites end at Column 12100

#Extract covariates (everything except metabolite intensity)
Metotal_no_metabolites_2 <-Metotal_2 %>%  #[,c(1:190,"gast","nullp")] #metotaldemo
  dplyr::select(1:200,gast,nullp)


#Clean covariates
table(Metotal_no_metabolites_2$stage, useNA= "always")    #-   1   2   P 
#1 330 270  49 

Metotal_no_metabolites_2 <- Metotal_no_metabolites_2 %>%
  filter(stage ==1 | stage == 2)

table(Metotal_no_metabolites_2$stage, useNA = "always") #   1    2 <NA> 
                                                        # 330  270    0
nrow(Metotal_no_metabolites_2) #600
# Metotal_no_metabolites_stage_1 <-subset(Metotal_no_metabolites,stage==1)   
# nrow(Metotal_no_metabolites_stage_1) #330

#Check # of observations for each complication
##GHTN
table(Metotal_no_metabolites_2$Comp_gHTNvsHealthyFullTerm, useNA= "always") 
# 0    1 <NA> 
# 302   77  221 

##GDM
table(Metotal_no_metabolites_2$Comp_gDMvsHealthyFullTerm, useNA= "always")
# 0    1 <NA> 
# 302   23  275  

##IUGR
table(Metotal_no_metabolites_2$Comp_IUGRvsHealthyFullTerm, useNA= "always")
#   0    1 <NA> 
# 302   37  261

##PROM
table(Metotal_no_metabolites_2$Comp_PROMvsHealthyFullTerm, useNA= "always")
#   0    1 <NA> 
# 302   24  274


#Exclude spontaneous/elective (stillbirth) abortion outcomes
table(Metotal_no_metabolites_2$birthoutcome, useNA= "always")

Metotal_no_metabolites_noSTB_2 <-Metotal_no_metabolites_2 %>%
  dplyr::filter(!(birthoutcome=="Preterm-STILLBIRTH" ))

nrow(Metotal_no_metabolites_noSTB_2)
#[1] 599

table(Metotal_no_metabolites_noSTB_2$birthoutcome)  
#EarlyTerm  FullTerm   Preterm 
#97       395       107 

table(Metotal_no_metabolites_noSTB_2$stage)
#   1   2 
# 329 270

#Reexamine complications 
##GHTN
table(Metotal_no_metabolites_noSTB_2$Comp_gHTNvsHealthyFullTerm, useNA= "always") 
#0    1 <NA> 
#  302   77  220

##GDM
table(Metotal_no_metabolites_noSTB_2$Comp_gDMvsHealthyFullTerm, useNA= "always")
#0    1 <NA> 
#  302   23  274 

##IUGR
table(Metotal_no_metabolites_noSTB_2$Comp_IUGRvsHealthyFullTerm, useNA= "always")
#0    1 <NA> 
#  302   37  260 

##PROM
table(Metotal_no_metabolites_noSTB_2$Comp_PROMvsHealthyFullTerm, useNA= "always")
#0    1 <NA> 
#  302   24  273


#####Removing Duplicates#####
table(Metotal_no_metabolites_noSTB_2$subjectid)
table(table(Metotal_no_metabolites_noSTB_2$subjectid))   #1   2  
                                                         #65 267 
                                                    #Only 2 vists max, no subject id should appear >2 times

which(table(Metotal_no_metabolites_noSTB_2$subjectid)>2) #None

subjectids_2 <- data.frame(table(Metotal_no_metabolites_noSTB_2$subjectid)) %>% #confirm Freq of none = 3
  rename(subjectid = Var1)                                    
subjectids_2[which(subjectids_2$Freq>2),] #0


Metotal_no_metabolites_noSTB_2[Metotal_no_metabolites_noSTB_2$subjectid=="E0616",c("File.Name","Sample.ID","seq")]
Metotal_no_metabolites_noSTB_2[388,] == Metotal_no_metabolites_noSTB_2[389,]
Metotal_no_metabolites_noSTB_2[388,] == Metotal_no_metabolites_noSTB_2[402,]
Metotal_no_metabolites_noSTB_nodupes_2 <- Metotal_no_metabolites_noSTB_2 %>%
  filter(!(File.Name == "F069_200209_M429_175"))

sum(duplicated(Metotal_no_metabolites_noSTB_2[,-c(1,3,6)])) #0
sum(duplicated(Metotal_no_metabolites_noSTB_nodupes_2[,-c(1,3)])) #0

table(Metotal_no_metabolites_noSTB_nodupes_2$subjectid)
which(table(Metotal_no_metabolites_noSTB_nodupes_2$subjectid)>2) #None appear > 2 times

subjectids_nodupes_2 <- data.frame(table(Metotal_no_metabolites_noSTB_nodupes_2$subjectid)) %>% 
  rename(subjectid = Var1)                                    
subjectids_nodupes_2[which(subjectids_nodupes_2$Freq>2),] #0

#####Check Number of Unique Subjects#####
subjectids_1_time_2 <-rownames(subset(table(Metotal_no_metabolites_noSTB_nodupes_2$subjectid),
                                    table(Metotal_no_metabolites_noSTB_nodupes_2$subjectid)==1))  
subjectids_2_times_2 <- rownames(subset(table(Metotal_no_metabolites_noSTB_nodupes_2$subjectid),
                                      table(Metotal_no_metabolites_noSTB_nodupes_2$subjectid)!=1))

#Assign to metotalv1, the observations of metotalnmis whose subjectids are in rn1
metotalv1_2<-Metotal_no_metabolites_noSTB_nodupes_2[which(Metotal_no_metabolites_noSTB_nodupes_2$subjectid %in% subjectids_1_time_2),] ######single stage 
length(table(metotalv1_2$subjectid))    # 65 obs- 65 subjectids appear only once

#Assign to metotalv2, the observations of metotalnmis whose subjectids are in rn2
metotalv2_2<-Metotal_no_metabolites_noSTB_nodupes_2[which(Metotal_no_metabolites_noSTB_nodupes_2$subjectid %in% subjectids_2_times_2),] ###### both stages
length(table(metotalv2_2$subjectid))   #267 subjectids appear twice
# 534 obs- 267 subjects total, 2 visits each: 267+65=332 subjects
# 328 subjects visit1-270 visit2
#267 subjects each had 2 observations = 534 observations in metotalv2, + 64 subjects w/ only 1 observation in metotalv1 331 total subjects

table(metotalv1_2$stage)
metotalv1_2[which(metotalv1_2$stage == 2),] 
#332 total subjects; however, 3 only appeared once, in stage 2; therefore, 328 total stage 1 subjects/observations
#E0539, E0545, E0637
length(unique(Metotal_no_metabolites_noSTB_nodupes_2$subjectid))
#Confirm, 332 total

unique_subjects_2 <- rbind(metotalv1_2, metotalv2_2[metotalv2_2$stage==1,]) #332 subjects
unique_subjects_analytic_2 <- unique_subjects_2 %>% #329 subjects; exclude 3 subjects w/only visit 2
  filter(stage == 1)

unique_subjects_analytic_2$sub_use <- NA
unique_subjects_analytic_2$sub_use <- ifelse(unique_subjects_analytic_2$TobaccoUse_MR ==1 | 
                                               unique_subjects_analytic_2$TobaccoUse_MRorSR ==1 | 
                                               unique_subjects_analytic_2$tobacco.AnyInPreg ==1 |
                                               unique_subjects_analytic_2$AlcoholUse_MR ==1 | 
                                               unique_subjects_analytic_2$AlcoholUse_MRorSR ==1 | 
                                               unique_subjects_analytic_2$alcohol.AnyInPreg ==1  |
                                               unique_subjects_analytic_2$MarijuanaUse_MR ==1 | 
                                               unique_subjects_analytic_2$MarijuanaUse_MRorSR ==1 | 
                                               unique_subjects_analytic_2$marijuana.AnyInPreg, 1, 0)
table(unique_subjects_analytic_2$sub_use)

#####Create Table 1######
names(unique_subjects_analytic_2)

#Question: why only these variables?
#List of continuous variables in a
contVars<-c("age_enrollment","birthga","FirstPrenatalBMI")
#List of categorical variables in  a
catVars<-c("Education_4.level","income","PregInsurance_2.level",
           #my addition#
           "PresInsurance_3.level", 
           #my addition end#
           "nullp","priorabortion","Site_Prenatal","Depression_Preg", "MarriedCohab_Not","priorpreterm",
           "Sex","g_htn","g_dm","PROM","Comp_gHTNvsHealthyFullTerm","Comp_gDMvsHealthyFullTerm","Comp_PROMvsHealthyFullTerm",
           "Comp_IUGRvsHealthyFullTerm","mode","birthoutcome","TobaccoUse_MR","AlcoholUse_MR","MarijuanaUse_MR","sub_use")

##Create Table One for continuous variables from a listed in contVars
tab_Continuous_2<-CreateTableOne(vars=contVars, data = unique_subjects_analytic_2)
print(tab_Continuous_2, showAllLevels = TRUE, test = FALSE)

#Create Table One for categorical variables from a listed in catVars
tab_Categorical_2<-CreateCatTable(vars=catVars, data= unique_subjects_analytic_2, includeNA = TRUE)
print(tab_Categorical_2, catDigits = 3,pDigits = 3,showAllLevels = TRUE, test = FALSE)

unique_subjects_analytic_2$subjectid %in% Metotal_no_metabolites_noSTB_nodupes_2[which(Metotal_no_metabolites_noSTB_nodupes_2$stage==1),]$subjectid
Metotal_no_metabolites_noSTB_nodupes_2[which(Metotal_no_metabolites_noSTB_nodupes_2$stage==1),]$subjectid %in% unique_subjects_analytic_2$subjectid

#####Subset to stage 1 and recheck Table 1#####
Metotal_no_metabolites_noSTB_nodupes_s1_2 <- Metotal_no_metabolites_noSTB_nodupes_2 %>%
  filter(stage == 1)

contVars <- c(contVars, "gast")

##Create Table One for continuous variables from a listed in contVars
tab_Continuous_2<-CreateTableOne(vars=contVars, data = Metotal_no_metabolites_noSTB_nodupes_s1_2)
print(tab_Continuous_2, showAllLevels = TRUE, test = FALSE)

#Create Table One for categorical variables from a listed in catVars
tab_Categorical_2<-CreateCatTable(vars=catVars, data= Metotal_no_metabolites_noSTB_nodupes_s1_2, includeNA = TRUE)
print(tab_Categorical_2, catDigits = 3,pDigits = 3,showAllLevels = TRUE, test = FALSE)


#####Log transform metabolite intensities####
Metotal_metabolites_2 <-Metotal_2 %>%
  filter(stage == 1,
         birthoutcome=="EarlyTerm" | birthoutcome=="FullTerm" |birthoutcome=="Preterm",
         !(File.Name == "F069_200209_M429_175"))
Metotal_no_metabolites_noSTB_nodupes_s1_2$subjectid %in% Metotal_metabolites_2$subjectid
Metotal_metabolites_2$subjectid %in% Metotal_no_metabolites_noSTB_nodupes_s1_2$subjectid

#colnames(Metotal_metabolites_2[c(2,202,12101)])

Metotal_clean_transformed_2 <- Metotal_no_metabolites_noSTB_nodupes_s1_2 %>%
  left_join(Metotal_metabolites_2[,c(2,201:12100)], by = join_by("subjectid"))

colnames(Metotal_clean_transformed_2[c(203,12102)])

Metotal_clean_transformed_2[203:12102] <- log2(Metotal_clean_transformed_2[203:12102])
Metotal_clean_transformed_2[Metotal_clean_transformed_2 == -Inf] <- NA
saveRDS(Metotal_clean_transformed_2, here::here("Data","Thesis","C18 2020 Analytic Dataset.rds"))


#####2020 C18 Temperature vs Metabolite Regressions####

Metotal_clean_transformed_2 <- readRDS(here::here("Data","Thesis","C18 2020 Analytic Dataset.rds")) %>%
  mutate(temp_exposed_status = factor(temp_exposed_status, levels = c(0,1),labels = c("Unexposed", "Exposed")),
         Sex = factor(Sex, levels = c(1,2), labels = c("Male","Female")),
         PresInsurance_3.level = factor(PresInsurance_3.level, levels = c(1,2,3),labels = c("Low-Income Medicaid","Right from the Start Medicaid","Private")),
         nullp = factor(nullp,levels = c(0,1), labels = c("Yes Nulliparous","Not Nulliparous")),
         AlcoholUse_MRorSR = factor(AlcoholUse_MRorSR, levels = c(0,1), labels = c("No","Yes")),
         TobaccoUse_MRorSR = factor(TobaccoUse_MRorSR, levels = c(0,1),labels = c("No","Yes")),
         conception_season = factor(conception_season)
         )

exposure_threshold_0.95 <- quantile(Metotal_clean_transformed_2$temperature,probs = (0.95))
exposure_threshold_0.75 <- quantile(Metotal_clean_transformed_2$temperature,probs = (0.75))

Metotal_clean_transformed_2<-Metotal_clean_transformed_2%>%
  mutate(temp_exposed_status_0.95 = case_when(temperature >= exposure_threshold_0.95 ~ 1,
                                              temperature < exposure_threshold_0.95 ~ 0),
         temp_exposed_status_0.75 = case_when(temperature >= exposure_threshold_0.75~1,
                                              temperature < exposure_threshold_0.75 ~0),
         Alc_and_Tobacco = case_when(AlcoholUse_MRorSR == "Yes" ~ 1,
                                     TobaccoUse_MRorSR == "Yes" ~ 1,
                                     TRUE ~ 0))


hist(Metotal_clean_transformed_2$temperature)
plot(Metotal_clean_transformed_2$gast, Metotal_clean_transformed_2$temperature)
plot(Metotal_clean_transformed_2$temperature, Metotal_clean_transformed_2$Comp_gHTNvsHealthyFullTerm)
plot(Metotal_clean_transformed_2$temperature, Metotal_clean_transformed_2$Comp_gDMvsHealthyFullTerm)
plot(Metotal_clean_transformed_2$temperature, Metotal_clean_transformed_2$Comp_IUGRvsHealthyFullTerm)
plot(Metotal_clean_transformed_2$temperature, Metotal_clean_transformed_2$Comp_PROMvsHealthyFullTerm)

is.factor(Metotal_clean_transformed_2$conception_season)

colnames(Metotal_clean_transformed_2[203:207]) #Metabolites begin at Column 203
colnames(Metotal_clean_transformed_2[12099:12102]) #Metabolites end at Column 12102

#####Categorical Temperature####
C18_cat_temp_vs_metabolite <-data.frame(matrix(nrow = nrow(c18_results), ncol = 4))
colnames(C18_cat_temp_vs_metabolite)<-c("Est","Std","t-value","P")

i<-1
for (i in 1:nrow(c18_results) ) {
  tryCatch(
    {
      lmfit<-lm(Metotal_clean_transformed_2[,i+202] ~temp_exposed_status+age_enrollment+FirstPrenatalBMI+Sex+
                  PresInsurance_3.level+nullp + AlcoholUse_MRorSR + TobaccoUse_MRorSR + gast + conception_season +
                  humidity, 
                data=Metotal_clean_transformed_2)
      C18_cat_temp_vs_metabolite[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}

C18_cat_temp_vs_metabolite$pfdr<-p.adjust(C18_cat_temp_vs_metabolite$P,method="BH")
write.csv(C18_cat_temp_vs_metabolite,here::here("Data","Thesis","C18 Categorical Temperature Regressions.csv"))


#####Continuous Temperature####
C18_cont_temp_vs_metabolite <-data.frame(matrix(nrow = nrow(c18_results), ncol = 4))
colnames(C18_cont_temp_vs_metabolite)<-c("Est","Std","t-value","P")

i<-1
for (i in 1:nrow(c18_results) ) {
  tryCatch(
    {
      lmfit<-lm(Metotal_clean_transformed_2[,i+202] ~temperature+age_enrollment+FirstPrenatalBMI+Sex+
                  PresInsurance_3.level+nullp + AlcoholUse_MRorSR + TobaccoUse_MRorSR + gast + conception_season +
                  humidity, 
                data=Metotal_clean_transformed_2)
      C18_cont_temp_vs_metabolite[i,]<-summary(lmfit)$coefficients[2,]
      if(i%%100 ==0) {print(i)}
    }, error=function(e){})
}

C18_cont_temp_vs_metabolite$pfdr<-p.adjust(C18_cont_temp_vs_metabolite$P,method="BH")
write.csv(C18_cont_temp_vs_metabolite,here::here("Data","Thesis","C18 Continuous Temperature Regressions.csv"))




















