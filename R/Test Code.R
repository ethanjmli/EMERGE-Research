precomp <- read.csv("h_g_htn_2017_s2.csv")
precomp2 <- read.csv("h_pre_complication_2017_s2.csv")

precompT <- precomp$X[precomp$P < 0.0005]
precompT2 <- precomp2$X[precomp2$P < 0.0005]

for(i in 1:length(lmfit_crude))
{
  print(table(is.na(lmfit_crude[[i]][4])))
}

for(i in 1:length(precompT2))
{
  ifelse(precompT2[i] %in% precompT,print(precompT2[i]),NA)
}

table(precompT %in% precompT2)

precompT
precompT2

precompT3 <- (precompT %in% precompT2)

length(precomp$X[precomp$P <0.05 & precomp2$P < 0.05])

#HILIC Stage 1 Pairwise Overlaps

setwd("/Users/ethan_j_li/OneDrive - Emory University/Research/2020/HILIC PreComp Regressions 2020")
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
    
    #length(precomp$X[precomp$P <0.05 & precomp2$P < 0.05])
    overlap2[overlaprow,1] <- paste("s1_",substring(list.filenames[firstdata],1,nchar(list.filenames[firstdata])-7), 
                                   substring(list.filenames[seconddata],1,nchar(list.filenames[seconddata])-7), 
                                   sep = "+")
    overlap2[overlaprow,2] <- length(df1$X[df1$P <0.05 & df2$P < 0.05])
    overlap2[overlaprow,3] <- length(df1$X[df1$P < 0.01 & df2$P < 0.01])
    overlap2[overlaprow,4] <- length(df1$X[df1$P<0.005 & df2$P<0.005])
    overlap2[overlaprow,5] <- length(df1$X[df1$P<0.0005 & df2$P<0.0005])
    overlap2[overlaprow,6] <- length(df1$X[df1$pdfr<0.2 & df2$pdfr<0.2])
    overlap2[overlaprow,7] <- length(df1$X[df1$pdfr<0.05 & df2$pdfr<0.05])
    overlaprow <- overlaprow+1
  }
  a<-a+1
}
