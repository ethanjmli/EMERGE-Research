pacman::p_load(ncdf4,
               tidyverse,
               sf,
               here,
               tidygeocoder,
               tmap,
               raster)
tmap_options(check.and.fix = TRUE)

####Geocoding####
######Google API: AIzaSyBZExW7V-cId5JNsDf9BX558t-PC9gxuxk####
setwd(here())
addresses <- read.csv(here::here("Data","AtlantaAAMC_Prenatal Address Data_9.22.2023, Visit1 and 2.csv"))
summary(addresses)
bbox <- addresses %>%
  filter(Visit == 1)%>%
  mutate(full_address = paste0(StreetAddress, ", ",City, ", ", State)) %>%
  geocode(address = full_address, method = "arcgis",lat = latitude_arcgis , long = longitude_arcgis)
  
sum(is.na(bbox$latitude_arcgis))
sum(is.na(bbox$longitude_arcgis))

write.csv(bbox, here::here("Data","Atlanta AAMC Geocoded Visit 1 Prenatal Addresses.csv"))

######
bbox <- read.csv(here::here("Data","Atlanta AAMC Geocoded Visit 1 Prenatal Addresses.csv"))

bbox.sf <- bbox %>%
  st_as_sf(coords = c("longitude_arcgis","latitude_arcgis"),crs=4269)

five_c <- st_read(here::here("Data","five_counties","five_counties.shp")) %>%
  st_transform(crs = 4269)

plot(st_geometry(bbox.sf))

tmap_options(check.and.fix = TRUE)
tm_shape(bbox.sf)+
  tm_dots()+
tm_shape(five_c)+
  tm_borders()

st_bbox(bbox.sf)

#      xmin      ymin      xmax      ymax 
# -84.80931  33.29111 -83.62932  34.12155 
rm(five_c)

####Test NCDF####
test_nc_data <- ncdf4::nc_open(file.path("H:/HRLDAS Raw ncdf/2010-2018 HRLDAS Raw ncdf/uhi_epi_201007_13_georgia.nc4"))
print(test_nc_data)

attributes(test_nc_data$var)
attributes(test_nc_data$dim)


lat <- ncvar_get(test_nc_data, "lat")
lat_fill <- ncatt_get(test_nc_data, "lat","_FillValue")
lat[lat == lat_fill]<-NA
lat

lon <- ncvar_get(test_nc_data, "lon")
lon_fill <- ncatt_get(test_nc_data, "lon","_FillValue")
lon[lon == lon_fill] <- NA
lon = lon-360

#print(c(nlon, nlat))
# 
# lat_df <- as.data.frame(lat)
# lon_df <- as.data.frame(lon)

time <- ncvar_get(test_nc_data, "time")
as.Date.numeric(time,origin = "1981-01-01")
date <- as.Date.numeric(time, origin = "1981-01-01")


t2_max_array <- ncvar_get(test_nc_data, "t2_max_cor") 
t2_fillvalue <- ncatt_get(test_nc_data, "t2_max_cor", "_FillValue")
t2_max_array[t2_max_array == t2_fillvalue$value] <- NA
t2_max_array_Celsius <- t2_max_array - 273 #Convert Kelvin to Celsius

humid <- ncvar_get(test_nc_data, "q2")
humid_fillvalue <- ncatt_get(test_nc_data, "q2","_FillValue")
humid[humid == humid_fillvalue$value] <- NA
#humidunits <- ncatt_get(test_nc_data, "q2","units")

c(lon)
c(lat)
c(t2_max_array_Celsius)
c(humid)

lonlat <- as.data.frame(cbind(c(lon),c(lat),c(t2_max_array_Celsius),c(humid))) %>%
  rename(Longitude = V1,
         Latitude = V2,
         t_max_Celsius = V3,
         humidity = V4)

lonlat.sf <- lonlat %>%
  st_as_sf(coords = c("Longitude","Latitude"),crs=4269)

nearestpoints <- st_nearest_feature(bbox.sf,lonlat.sf)

tm_shape(bbox.sf)+
  tm_dots(col = "red",size = 0.1)+
tm_shape(lonlat.sf)+
  tm_dots(col = "t_max_Celsius")

sum(is.na(lonlat.sf$t_max_Celsius))
sum(is.na(lonlat.sf$humidity))


bbox.sf$temp_max <- lonlat.sf$t_max_Celsius[nearestpoints]
bbox.sf$humidity <- lonlat.sf$humidity[nearestpoints]

here()
####Full Process NCDF####

setwd("H:/HRLDAS Raw ncdf/2010-2018 HRLDAS Raw ncdf/")
raw_nc_list <- list.files()
i<-1


for(i in 1:length(raw_nc_list)){
  raw_nc_file <- ncdf4::nc_open(file.path(raw_nc_list[i]))
  
  lat <- ncvar_get(raw_nc_file, "lat")
  lat_fill <- ncatt_get(raw_nc_file, "lat","_FillValue")
  lat[lat == lat_fill]<-NA
  
  lon <- ncvar_get(raw_nc_file, "lon")
  lon_fill <- ncatt_get(raw_nc_file, "lon","_FillValue")
  lon[lon == lon_fill] <- NA
  lon = lon-360
  
  
  time <- ncvar_get(raw_nc_file, "time")
  date <- as.Date.numeric(time, origin = "1981-01-01")
  
  t2_max_array <- ncvar_get(raw_nc_file, "t2_max_cor") 
  t2_fillvalue <- ncatt_get(raw_nc_file, "t2_max_cor", "_FillValue")
  t2_max_array[t2_max_array == t2_fillvalue$value] <- NA
  t2_max_array <- t2_max_array - 273 #Convert Kelvin to Celsius
  
  humid <- ncvar_get(raw_nc_file, "q2")
  humid_fillvalue <- ncatt_get(raw_nc_file, "q2","_FillValue")
  humid[humid == humid_fillvalue$value] <- NA
  
  lonlat.sf <- as.data.frame(cbind(c(lon),c(lat),c(t2_max_array),c(humid))) %>%
    rename(Longitude = V1,
           Latitude = V2,
           t_max_Celsius = V3,
           humidity = V4)
    # st_as_sf(coords = c("Longitude","Latitude"),crs=4269) %>%
    # cbind(st_coordinates(lonlat.sf))%>%
    # rename(Longitude = X,
    #        Latitude = Y)
  
  lonlat_writing <- lonlat.sf%>%
    # st_drop_geometry() %>%
    filter(row_number() %in% nearestpoints)
    
  saveRDS(lonlat_writing, here::here("Data","Temperature Data",paste(date, "Maximum Temperatures.rds")))  
  print(i)
}


a<- readRDS(here::here("Data","Temperature Data",paste(date, "Maximum Temperatures.rds")))
b <- readRDS(here::here("Data","Temperature Data",paste(date-1, "Maximum Temperatures.rds")))
c <- readRDS(here::here("Data","Temperature Data", "2010-08-15 Maximum Temperatures.rds"))
d <- readRDS(here::here("Data","Temperature Data", "2010-12-15 Maximum Temperatures.rds"))
e <- readRDS(here::here("Data","Temperature Data", "2012-01-15 Maximum Temperatures.rds"))
