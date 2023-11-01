#HILIC 2020 Stage 1 Pairwise Overlaps

setwd("/Users/ethan_j_li/Documents/Emory/Research/HILIC PreComp Regressions 2020")
list.files(pattern="h.*s1*.csv$") 
list.filenames <- list.files(pattern="h.*s1*.csv$")

# create an empty list that will serve as a container to receive the incoming files
lmfit_crude <- list()

# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  lmfit_crude[[i]] <- read.csv(list.filenames[i])
}

rm(i)

# add the names of these files to the list
names(lmfit_crude) <- list.filenames
str(lmfit_crude)

#Create dataframe to hold overlap counts
overlap <- as.data.frame(matrix(nrow=10,ncol=7))
names(overlap) <- c("Comparison", "Overlap P<0.05", "Overlap P<0.01", "Overlap P<0.005", "Overlap P<0.0005","Overlap BHP<0.2","Overlap BHP<0.05")

#Pairwise Comparisons
a<-2
overlaprow <- 1
for (firstdata in 1:(length(lmfit_crude)-1))
{
  df1 <- lmfit_crude[[firstdata]]
  
  for(seconddata in a:length(lmfit_crude))
  {
    df2 <- lmfit_crude[[seconddata]]
    
    
    overlap[overlaprow,1] <- paste("s1_",substring(list.filenames[firstdata],1,nchar(list.filenames[firstdata])-7), 
                                    substring(list.filenames[seconddata],1,nchar(list.filenames[seconddata])-7), 
                                    sep = "+")
    overlap[overlaprow,2] <- length(df1$X[df1$P <0.05 & df2$P < 0.05])
    overlap[overlaprow,3] <- length(df1$X[df1$P < 0.01 & df2$P < 0.01])
    overlap[overlaprow,4] <- length(df1$X[df1$P<0.005 & df2$P<0.005])
    overlap[overlaprow,5] <- length(df1$X[df1$P<0.0005 & df2$P<0.0005])
    overlap[overlaprow,6] <- length(df1$X[df1$pdfr<0.2 & df2$pdfr<0.2])
    overlap[overlaprow,7] <- length(df1$X[df1$pdfr<0.05 & df2$pdfr<0.05])
    overlaprow <- overlaprow+1
    print(paste(firstdata, seconddata,sep = " X "))
  }
  a<-a+1
}

write.csv(overlap,"HILIC Stage 1 Pairwise Overlaps.csv")

#############
#############

#HILIC 2020 Stage 2 Pairwise Overlaps 

setwd("/Users/ethan_j_li/Documents/Emory/Research/HILIC PreComp Regressions")
list.files(pattern="h.*s2*.csv$") 
list.filenames <- list.files(pattern="h.*s2*.csv$")

# create an empty list that will serve as a container to receive the incoming files
lmfit_crude <- list()

# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  lmfit_crude[[i]] <- read.csv(list.filenames[i])
}

rm(i)

# add the names of these files to the list
names(lmfit_crude) <- list.filenames
str(lmfit_crude)

#Create dataframe to hold overlap counts
overlap <- as.data.frame(matrix(nrow=10,ncol=7))
names(overlap) <- c("Comparison", "Overlap P<0.05", "Overlap P<0.01", "Overlap P<0.005", "Overlap P<0.0005","Overlap BHP<0.2","Overlap BHP<0.05")

#Pairwise Comparisons
a<-2
overlaprow <- 1
for (firstdata in 1:(length(lmfit_crude)-1))
{
  df1 <- lmfit_crude[[firstdata]]
  
  for(seconddata in a:length(lmfit_crude))
  {
    df2 <- lmfit_crude[[seconddata]]
    
    
    overlap[overlaprow,1] <- paste("s2_",substring(list.filenames[firstdata],1,nchar(list.filenames[firstdata])-7), 
                                   substring(list.filenames[seconddata],1,nchar(list.filenames[seconddata])-7), 
                                   sep = "+")
    overlap[overlaprow,2] <- length(df1$X[df1$P <0.05 & df2$P < 0.05])
    overlap[overlaprow,3] <- length(df1$X[df1$P < 0.01 & df2$P < 0.01])
    overlap[overlaprow,4] <- length(df1$X[df1$P<0.005 & df2$P<0.005])
    overlap[overlaprow,5] <- length(df1$X[df1$P<0.0005 & df2$P<0.0005])
    overlap[overlaprow,6] <- length(df1$X[df1$pdfr<0.2 & df2$pdfr<0.2])
    overlap[overlaprow,7] <- length(df1$X[df1$pdfr<0.05 & df2$pdfr<0.05])
    overlaprow <- overlaprow+1
    print(paste(firstdata, seconddata,sep = " X "))
  }
  a<-a+1
}

write.csv(overlap,"HILIC Stage 2 Pairwise Overlaps.csv")

#############
#############

#C18 2020 Stage 1 Pairwise Overlaps

setwd("/Users/ethan_j_li/Documents/Emory/Research/C18 PreComp Regressions")
list.files(pattern="c.*s1*.csv$") 
list.filenames <- list.files(pattern="c.*s1*.csv$")

# create an empty list that will serve as a container to receive the incoming files
lmfit_crude <- list()

# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  lmfit_crude[[i]] <- read.csv(list.filenames[i])
}

rm(i)

# add the names of these files to the list
names(lmfit_crude) <- list.filenames
str(lmfit_crude)

#Create dataframe to hold overlap counts
overlap2 <- as.data.frame(matrix(nrow=10,ncol=7))
names(overlap2) <- c("Comparison", "Overlap P<0.05", "Overlap P<0.01", "Overlap P<0.005", "Overlap P<0.0005","Overlap BHP<0.2","Overlap BHP<0.05")

#Pairwise Comparisons
a<-2
overlaprow <- 1
for (firstdata in 1:(length(lmfit_crude)-1))
{
  df1 <- lmfit_crude[[firstdata]]
  
  for(seconddata in a:length(lmfit_crude))
  {
    df2 <- lmfit_crude[[seconddata]]
    
    
    overlap2[overlaprow,1] <- paste("s1_",substring(list.filenames[firstdata],1,nchar(list.filenames[firstdata])-7), 
                                   substring(list.filenames[seconddata],1,nchar(list.filenames[seconddata])-7), 
                                   sep = "+")
    overlap2[overlaprow,2] <- length(df1$X[df1$P <0.05 & df2$P < 0.05 & is.na(df1$P) == F])
    overlap2[overlaprow,3] <- length(df1$X[df1$P < 0.01 & df2$P < 0.01 & is.na(df1$P) == F])
    overlap2[overlaprow,4] <- length(df1$X[df1$P<0.005 & df2$P<0.005 & is.na(df1$P) == F])
    overlap2[overlaprow,5] <- length(df1$X[df1$P<0.0005 & df2$P<0.0005 & is.na(df1$P) == F])
    overlap2[overlaprow,6] <- length(df1$X[df1$pdfr<0.2 & df2$pdfr<0.2 & is.na(df1$P) == F])
    overlap2[overlaprow,7] <- length(df1$X[df1$pdfr<0.05 & df2$pdfr<0.05 & is.na(df1$P) == F])
    overlaprow <- overlaprow+1
    print(paste(firstdata, seconddata,sep = " X "))
  }
  a<-a+1
}

write.csv(overlap2,"C18 Stage 1 Pairwise Overlaps.csv")

#############
#############

#C18 2020 Stage 2 Pairwise Overlaps

setwd("/Users/ethan_j_li/Documents/Emory/Research/C18 PreComp Regressions")
list.files(pattern="c.*s2*.csv$") 
list.filenames <- list.files(pattern="c.*s2*.csv$")

# create an empty list that will serve as a container to receive the incoming files
lmfit_crude <- list()

# create a loop to read in the list
for (i in 1:length(list.filenames))
{
  lmfit_crude[[i]] <- read.csv(list.filenames[i])
}

rm(i)

# add the names of these files to the list
names(lmfit_crude) <- list.filenames
str(lmfit_crude)

#Create dataframe to hold overlap counts
overlap <- as.data.frame(matrix(nrow=10,ncol=7))
names(overlap) <- c("Comparison", "Overlap P<0.05", "Overlap P<0.01", "Overlap P<0.005", "Overlap P<0.0005","Overlap BHP<0.2","Overlap BHP<0.05")

#Pairwise Comparisons
a<-2
overlaprow <- 1
for (firstdata in 1:(length(lmfit_crude)-1))
{
  df1 <- lmfit_crude[[firstdata]]
  
  for(seconddata in a:length(lmfit_crude))
  {
    df2 <- lmfit_crude[[seconddata]]
    
    
    overlap[overlaprow,1] <- paste("s2_",substring(list.filenames[firstdata],1,nchar(list.filenames[firstdata])-7), 
                                   substring(list.filenames[seconddata],1,nchar(list.filenames[seconddata])-7), 
                                   sep = "+")
    overlap[overlaprow,2] <- length(df1$X[df1$P <0.05 & df2$P < 0.05])
    overlap[overlaprow,3] <- length(df1$X[df1$P < 0.01 & df2$P < 0.01])
    overlap[overlaprow,4] <- length(df1$X[df1$P<0.005 & df2$P<0.005])
    overlap[overlaprow,5] <- length(df1$X[df1$P<0.0005 & df2$P<0.0005])
    overlap[overlaprow,6] <- length(df1$X[df1$pdfr<0.2 & df2$pdfr<0.2])
    overlap[overlaprow,7] <- length(df1$X[df1$pdfr<0.05 & df2$pdfr<0.05])
    overlaprow <- overlaprow+1
    print(paste(firstdata, seconddata,sep = " X "))
  }
  a<-a+1
}

write.csv(overlap,"C18 Stage 2 Pairwise Overlaps.csv")
