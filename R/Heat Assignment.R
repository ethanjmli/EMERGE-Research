pacman::p_load(here,
               stringr,
               sf,
               tidyverse,
               hydroTSM)


####
temp_file_list <- list.files(here::here("Data","Temperature Data"))
file_dates <- as.Date(str_sub(temp_file_list,1,-26))
diagnostic_temp <- NULL

for(i in 1:length(temp_file_list)){
  print(file_dates[i])
  diagnostic_file <- readRDS(here::here("Data","Temperature Data",temp_file_list[i])) 
  diagnostic_temp <- c(diagnostic_temp, mean(diagnostic_file$t_max_Celsius))
}

plot(file_dates,diagnostic_temp)

####

subject_locations <- read.csv(here::here("Data","Atlanta AAMC Geocoded Visit 1 Prenatal Addresses.csv"))%>%
  st_as_sf(coords = c("longitude_arcgis","latitude_arcgis"), crs = 4269)

demographicsbystudyid_1<-read.csv(here::here("Data","MPTB_AA Cohort Data_N547_addZscore.csv")) %>%
  mutate(datevisit1 = as.Date(datevisit1, format = "%m/%d/%Y"),
         date_conception = datevisit1 - gasamp1_wk.10*7) %>%
  left_join(subject_locations, by = join_by(subjectid == SubjectID)) %>%
  filter(datevisit1 <= as.Date("2018-12-31"))

# a <- demographicsbystudyid %>%
#   filter(datevisit1 > as.Date("2018-12-31")) 49 subjects 
# [1] "E0737" "E0743" "E0744" "E0745" "E0746" "E0748" "E0749" "E0750" "E0751" "E0752" "E0753" "E0754" "E0755" "E0756" "E0757" "E0758" "E0759"
# [18] "E0760" "E0761" "E0763" "E0764" "E0766" "G0327" "G0358" "G0363" "G0364" "G0365" "G0368" "G0370" "G0371" "G0373" "G0375" "G0376" "G0377"
# [35] "G0378" "G0379" "G0380" "G0381" "G0382" "G0386" "G0387" "G0388" "G0389" "G0390" "G0394" "G0395" "G0397" "G0401" "G0403"



demographicsbystudyid_1$temperature <- NA
demographicsbystudyid_1$humidity <- NA


temp_file_list <- list.files(here::here("Data","Temperature Data"))
file_dates <- as.Date(str_sub(temp_file_list,1,-26))



subject <- 1

for(subject in 1:nrow(demographicsbystudyid_1)){
  
  temperature <- NULL
  humidity <- NULL
  
  date_sequence <- seq.Date(from = demographicsbystudyid_1$date_conception[subject], 
                            to =demographicsbystudyid_1$datevisit1[subject],
                            by="day")

  
  for(i in which(file_dates %in% date_sequence)){
    
    temperature_file <- readRDS(here::here("Data","Temperature Data",temp_file_list[i])) %>%
                                  st_as_sf(coords = c("Longitude","Latitude"),crs = 4269)
    
    nearest_point <- st_nearest_feature(demographicsbystudyid_1$geometry[subject], temperature_file)
    
    temperature <- c(temperature, temperature_file$t_max_Celsius[nearest_point])
    
    humidity <- c(humidity, temperature_file$humidity[nearest_point])
    print(i)
    
  }
  
  avg_temp <- mean(temperature, na.rm=T)
  avg_humidity <- mean(humidity, na.rm=T)
  demographicsbystudyid_1[subject,"temperature"] <- avg_temp
  demographicsbystudyid_1[subject,"humidity"] <- avg_humidity
  print(subject)
  
}


plot(demographicsbystudyid_1$date_conception,demographicsbystudyid_1$temperature)
plot(demographicsbystudyid_1$datevisit1,demographicsbystudyid_1$temperature)
plot(demographicsbystudyid_1$X.x,demographicsbystudyid_1$temperature)
plot(demographicsbystudyid_1$gasamp1_wk.10,demographicsbystudyid_1$temperature)
testplot <- lm(temperature ~ gasamp1_wk.10, data = demographicsbystudyid_1)
summary(testplot)
abline(testplot)


class(demographicsbystudyid_1$X.x)
class(demographicsbystudyid_1$X.y)
demographicsbystudyid_1$X.x == demographicsbystudyid_1$X.y

exposure_threshold <- quantile(demographicsbystudyid_1$temperature, probs = 0.9, na.rm=T)

demographicsbystudyid_1 <- demographicsbystudyid_1 %>%
  mutate(temp_exposed_status = ifelse(temperature >= exposure_threshold, 1, 0),
         temp_exposed_status = factor(temp_exposed_status, levels = c(0,1)),
         conception_season = hydroTSM::time2season(date_conception, out.fmt = "seasons"),
         conception_season = factor(conception_season,ordered=FALSE)) %>%
  select(-c(X.x,X.y,geometry)) 

is.factor(demographicsbystudyid_1$conception_season)

table(demographicsbystudyid_1$temp_exposed_status,demographicsbystudyid_1$conception_season)
#   spring  summer fall winter
# 0    110    130    0    105
# 1     18     32    0      0

write.csv(demographicsbystudyid_1, here::here("Data","AA Cohort Demographics + Temperature + Humidity.csv"))
