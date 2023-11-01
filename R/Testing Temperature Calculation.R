pacman::p_load(here,
               stringr,
               sf)
here()

date_sequence <- seq.Date(from = as.Date("2012-01-01"), to =as.Date("2012-04-05"),by="day")

temp_file_list <- list.files(here::here("Data","Temperature Data"))

file_dates <- as.Date(str_sub(temp_file_list,1,-26))

subject_locations <- read.csv(here::here("Data","Atlanta AAMC Geocoded Visit 1 Prenatal Addresses.csv"))%>%
                                st_as_sf(coords = c("longitude_arcgis","latitude_arcgis"), crs = 4269)

temperature <- NULL
humidity <- NULL
subject_locations$temperature <- NA
subject_locations$humidity <- NA

for(subject in 1:nrow(subject_locations)){


  for(i in which(file_dates %in% date_sequence)){
  
    temperature_file <- readRDS(here::here("Data","Temperature Data",temp_file_list[i])) %>%
                                  st_as_sf(coords = c("Longitude","Latitude"),crs = 4269)
    nearest_point <- st_nearest_feature(subject_locations$geometry[subject], temperature_file)
    
    temperature <- c(temperature, temperature_file$t_max_Celsius[nearest_point])
    
    humidity <- c(humidity, temperature_file$humidity[nearest_point])
    print(i)
  
       }
  
  avg_temp <- mean(temperature, na.rm=T)
  avg_humidity <- mean(humidity, na.rm=T)
  subject_locations[subject,"temperature"] <- avg_temp
  subject_locations[subject,"humidity"] <- avg_humidity
  print(subject)
}








