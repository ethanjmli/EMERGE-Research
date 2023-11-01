#R Coding for Metabolomics Analysis#

#Clean up memory by setting list to an empty list

rm(list=ls())

install.packages("lmerTest")

#required packages)
library(lme4)
library(lmerTest) #below 2 needed to execute lienar mixed effect mode
library(readr)

#Skipping feature extraction

#Step 1.1 Import and clean the Metabolomic Feature Table#
dir_in='/Users/ethan_j_li/Documents/Emory/'
dir_out = 'Users/ethan_j_li/Emory/Research/Output'

#Import feature table
mft <- read_csv(paste0(dir_in,'BIOS 500 HW 1.xlsx'))

#Check dimensions of the dataset
dim(mft)
View(mft)

#Variables of interest: m/z ratio (mz), column retention time (time) identify#
#unique metabolic features#

#Check summary of first ten columns
summary(mft[,1:10])

#Check data quality
hist(mft$NumPres.All.Samples) #Distribution of frequency of metabolite appearance
sum(mft$NumPres.All.Samples<29) #how many metabolic features are present in <29 samples

hist(mft$median_CV) #Distribution of median of coefficient of varaince
sum(mft$median_CV>30) #How many metabolites have CV>30

plot(mft$Qscore,mft$median_CV)

#Quality Assurance/Control Filtering on metabolics feature table

#Create new dataset mdat which is a subset of mft, consisting of all 
#metabolites in mft with CV <=30 and present in >29 samples
mdat <- subset(mft,mft$median_CV<=30&mft$NumPres.All.Samples>29)
dim(mdat)

#Check data quality on new dataset
hist(mdat$NumPres.All.Samples)
sum(mdat$NumPres.All.Samples < 29)

hist(mdat$median_CV)
sum(mdat$median_CV>30)

#Import sample sequence ID file
#Match blinded file names to sample ids
infofile <- read.table(file = past0(dir_in, "plasma_sequencefileHILIC.txt"), header = T)
view(infofile)

#Merge info file with feature table sample column names
#Replace randomly assigned sample column names with meaningful sample ids

#Isolate sample column names from mdat, columns 10--> end
colnames(mdat)[10:dim(mdat)[2]]
colnames.mdat<-colnames(mdat)[10:dim(mdat)[2]]
#Substring all sample column names, keep all characters until last 5 characters, 
#save substrings in new file
File.n <- as.data.frame(substr(colnames.mdat,1,nchar(colnames.mdat)-5))

#Merge two datasets and see if sample info matches with column names in feature table
infofile.mdat <-merge.data.frame(File.n,infofile, by='File_Name', all.x = T, sort = F)
View(infofile.mdat)
#If there is a discrepancy, there will be NA value in the sample ID column
#Check is there NA value?
sum(is.na(infofile.dat[,2]))

#Replacing colnames in the feature table, aka Sequence ID, with sample ID
intensity <- as.matrix(mdat[,10:dim(mdat[2])])
View(intensity)
colnames(intensity) <- infofile.mdat$Sample_ID
View(intensity)
#Replace missing values (0) with NA
intensity[intensity == 0] <- NA
mdat1 <- cbind(mdat[,1:9],intensity)
View(mdat1)

#Start calculating means across triplicates
#Method-1: use for loops
means_1 <- data.frame(matrix(NA, nrow = dim(intensity[1],ncol = dim(intensity[2]/3))))
for (i in 1:dim(intensity[1]))
{
  for (j in 1:(dim(intensity[2]/3)))
  {
    means_1[i,j]<-mean(c(intensity[i,3*j-2],intensity[i,3*j-1],intensity[i,3*j]),na.rm=T)
  }
}

#Method-2 use lapply function, more time efficient than for loops
#lapply uses a list the same length as X, each element of which is the result of 
#applying'FUN' to the corresponding element 'X'
colnum <- 1:ncol(intensity) #calculate # of columns
ind <- as.data.frame(matrix(colnum,byrow = F, nrow = 3)) 
means <- as.data.frame(lapply(ind, function(i) rowMeans(intensity[,i],na.rm=T)))

#generate new feature tables with mean intensities across triplicates
samp_id <- substr(as.vector(colnames(intensity)), 1,nhcar(colnames(intensity))-2)[c(T,rep(F,2))]
colnames(means) <- samp_id
#combine first nine columns from mdat with newly calculated means
metabo <-cbind(mdat1[,1:9],means)

#