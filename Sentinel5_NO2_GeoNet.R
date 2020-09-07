############################################################################################################################
############################### Sentinel-5P: Processing and handling NO2 measured values in R ##############################
############################################################################################################################

############################################################################################################################
# Author: Magdalena Halbgewachs (ESRI Germany, JMU Wuerzburg)
# Date: 27.07.2020

# RStudio version: 1.2.5033
# R version: 4.0.2

# Overview:
## - 1.) Download Sentinel-5P NetCDF files for specific region and time period
## - 2.) Processing Sentinel-5P data
## - 2.1) Save NO2 concentrations per day as TIFF file (+ application of a cloud filter)
## - 2.2) Rolling Window: Save average values over a period of 5 days as TIFF
## - 3.) Deleting the NetCDFs after successfully creating the TIFF files


# install required packages
#install.packages("devtools")                                   # for downloading an R package from a GitHub account, the devtools are required
#devtools::install_github("16EAGLE/getSpatialData", force=TRUE)
#install.packages("raster")
#install.packages("sf")
#install.packages("sp")
#install.packages("dismo")
#install.packages("geosphere")
#install.packages("rgeos")
#devtools::install_github("MBalthasar/S5Processor")
#install.packages("ncdf4")
#install.packages("ggplot2")
#install.packages("maptools")
#install.packages("rgdal")
#install.packages("lubridate")
#install.packages("cat")


# import necessary packages
library(getSpatialData)
library(raster)
library(sf)
library(sp)
library(dismo)
library(geosphere)
library(rgeos)
library(S5Processor)
library(ncdf4)
library(ggplot2)
library(maptools)
library(rgdal)
library(lubridate)
library(cat)




############################################################################################################################
#########################################
#### 1.) Download Sentinel-5P images ####
#########################################

# Define AOI (investigation area) -> Load polygon (shape layer)
aoi_data <- raster::shapefile("C:\\data\\AOI_bbox.shp")                 # adjust manually

# Set AOI using Shapefile
set_aoi(aoi_data)

# Display AOI
#view_aoi()

#--------------------------------------------------------------------------------------------------------------------------
# Define a time range and a specific satellite

#time_range <-  c("2020-04-06", "2020-04-08")                           # time range in the past, after 1.) move on to 2.1b.)
time_range <-  c(as.character(Sys.Date()-1), as.character(Sys.Date()))  # date yesterday, after 1.) move on to 2.1a.)

platform <- "Sentinel-5P"                                               # only images from Sentinel-5P are requested

#--------------------------------------------------------------------------------------------------------------------------
# Log in to the Copernicus Hub account
login_CopHub(username = "my_username", password="my_password")          # add your own user name and password, adjust manually

#--------------------------------------------------------------------------------------------------------------------------
# Query in the Sentinel Hub
records <- getSentinel_records(time_range = time_range, platform = platform)
#colnames(records)                                                      # prints column names

# Filter the resulting scenes by product level and product type
records_filtered <- records[which(records$product_type == "L2__NO2___"),]# filtered by "NO2"


#--------------------------------------------------------------------------------------------------------------------------
# Display filtered scenes

#View(records)                      # by time period and satellite
#View(records_filtered)             # by product type

#--------------------------------------------------------------------------------------------------------------------------
# Download of filtered scenes

# Define the output path where the downloaded scenes are saved
set_archive("C:\\data\\Sentinel5\\raw")                                 # adjust manually

datasets <- function() {
  getSentinel_data(records = records_filtered)                          # filtered scenes as input variable
}

r <- NULL
attempt <- 1
while(is.null(r) && attempt <= 10){                                     # tries to download the scenes 10x before script stops
  attempt <- attempt+1
  try(
    r <- datasets()
  )
}


#-------------------------------------------------------------------------------------------------------------------------
# write log file

# Specify the path in which the log file is created
log_con <- file("C:\\data\\Sentinel5\\logs\\Sentinel5_R_log.txt", open="a") # adjust manually

#writes date and time to log file after successful download
cat("Download successful on", as.character(Sys.time()), "\n", file = log_con)


############################################################################################################################
################################################################
#### 2.) Create TIFF file from all scenes of one day ###########
################################################################

#--------------------------------------------------------------------------------------------------------------------------
# 2.1a) Current recordings (from yesterday)
#--------------------------------------------------------------------------------------------------------------------------
# This paragraph must be executed if NetCDF files were downloaded with yesterdays time stamp and need to be processed, otherwise continue with 2.1a)

# Lists all downloaded images
files_folder <- "C:\\Sentinel5\\raw\\datasets\\Sentinel-5P"		              # adjust manually


# lists all pictures from yesterday
files_today <- list.files(files_folder, pattern= as.character(Sys.Date()-1,"%Y%m%d"), full.names=T, all.files=F, include.dirs = F)# listet nur gestrige Szenen

# saves retrieved scenes in vector
one_day <- c(files_today)


# load polygon (shapefile) of the final investigation area
aoi <- raster::shapefile("E:\\data\\Sentinel5\\boundaries.shp")             # adjust manually
aoi <- spTransform(aoi, CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")) # define coordinate system
aoi_buff <- gBuffer(aoi, width = 5000)                                                                        # Buffering of the area so that no pixel is missing at the borders ("width" in meters)

# lists all variables stored in downloaded NetCDF
S5P_1 <- S5P_process(input = one_day)

#-----------------------------------------------------------NO2---------------------------------------------------------
# Merge all scenes of a day, retrieve a specific variable (product) and convert the unit to molecules/cm2
S5P_mask_unit <- S5P_process(input = one_day, my_res = 10000,              # resolution 10000m
                             product = 6, my_aoi = aoi_buff,               # product 6 = nitrogendioxide_tropospheric_column (troposphere mole content of NO2)
                             extent_only = FALSE,
                             apply_scale_factor = T)                       # conversion of mol/m2 to molecules/cm2


#-------------------------------------------------------CLOUD FILTER----------------------------------------------------
# Merge all scenes of a day, retrieve a specific variable (product)
S5P_mask_unit_cloud <- S5P_process(input = one_day, my_res = 10000,
                             product = 5, my_aoi = aoi_buff,               #product 5 = qa_value (< 0.5 is cloud)
                             extent_only = FALSE)

#-----------------------------------------------------EXTRACT CLOUDS--------------------------------------------------
# Remove cloud cover from scene
S5P_NO2_final <- S5P_mask_unit
S5P_NO2_final[S5P_mask_unit_cloud < 0.5] <- -999                          # -999 as NoData value
#plot(S5P_NO2_final)                                                      # plot the scene with cloud filter

#-------------------------------------------------------------------------------------------------------------------------
# Crop the layer to the area of the AOI 
aoi_clip <- spTransform(aoi, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))  # Assign a new coordinate system

S5P_NO2_final_crop <- crop(S5P_NO2_final, aoi_clip)                       # rough cropping
S5P_NO2_final_mask <- mask(S5P_NO2_final_crop, aoi_clip)                  # exact cropping

#plot(S5P_NO2_final_mask)                                                 # Plot final layer
#plot(aoi_clip, add=T)                                                    # Plot with AOI

#-------------------------------------------------------------------------------------------------------------------------
# save final layer as TIFF

# Path where the layer is saved
out_folder <- "C:\\data\\Sentinel5\\2020\\"                               # adjust manually

# Naming is adjusted to the respective date
name <- as.character(Sys.Date()-1, "%Y%m%d")

# save raster as TIFF
writeRaster(S5P_NO2_final_mask, filename = paste0(out_folder, name, ".tif"), format="GTiff",overwrite=T)

#--------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------
# 2.1b) Images of earlier dates
#--------------------------------------------------------------------------------------------------------------------------
# This paragraph must be executed if NetCDF files were downloaded from earlier times and need to be processed
# -- Date must be entered individually and cannot be derived automatically from the system time --

# Lists all downloaded images
#files_folder <- "C:\\Sentinel5\\raw\\datasets\\Sentinel-5P"		          # adjust manually


# lists all images of a certain date ("NO2____20190328")
#files_past <- list.files(files_folder, pattern= "NO2____20190328", full.names=T, all.files=F, include.dirs = F)# adjust date manually

# saves retrieved scenes in vector
#one_day <- c(files_past)


# load polygon (shapefile) of the final investigation area
#aoi <- raster::shapefile("E:\\data\\Sentinel5\\boundaries.shp")          # adjust manually
#aoi <- spTransform(aoi, CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")) # define coordinate system
#aoi_buff <- gBuffer(aoi, width = 5000) 

# lists all variables stored in downloaded NetCDF
#S5P_1 <- S5P_process(input = one_day)

#-----------------------------------------------------------NO2---------------------------------------------------------
# Merge all scenes of a day, retrieve a specific variable (product) and convert the unit to molecules/cm2
#S5P_mask_unit <- S5P_process(input = one_day, my_res = 10000,            # resolution 10000m
                             #product = 6, my_aoi = aoi_buff,             # product 6 = nitrogendioxide_tropospheric_column (troposphere mole content of NO2)
                             #extent_only = FALSE,
                             #apply_scale_factor = T)                     # conversion of mol/m2 to molecules/cm2

#-------------------------------------------------------CLOUD FILTER-------------------------------------------------------
# Merge all scenes of a day, retrieve a specific variable (product)
#S5P_mask_unit_cloud <- S5P_process(input = one_day, my_res = 10000,
                                   #product = 5, my_aoi = aoi_buff,       #product 5 = qa_value (< 0.5 is cloud)
                                   #extent_only = FALSE)

#-----------------------------------------------------EXTRACT CLOUDS--------------------------------------------------
# Remove cloud cover from scene
#S5P_NO2_final <- S5P_mask_unit
#S5P_NO2_final[S5P_mask_unit_cloud < 0.5] <- -999                         # -999 as NoData value
#plot(S5P_NO2_final)                                                      # plot the scene with cloud filter

#-------------------------------------------------------------------------------------------------------------------------
# Zuschneiden des Layers auf die Fl?che der AOI (L?ndergrenzen)
# Crop the layer to the area of the AOI 
#aoi_clip <- spTransform(aoi, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))  # Assign a new coordinate system

#S5P_NO2_final_crop <- crop(S5P_NO2_final, aoi_clip)                      # rough cropping
#S5P_NO2_final_mask <- mask(S5P_NO2_final_crop, aoi_clip)                 # exact cropping

#plot(S5P_NO2_final_mask)                                                 # Plot final layer
#plot(aoi_clip, add=T)                                                    # Plot with AOI

#-------------------------------------------------------------------------------------------------------------------------
# save final layer as TIFF

# Path where the layer is saved
#out_folder <- "C:\\data\\Sentinel5\\2020\\"                              # adjust manually

# Naming is adjusted to the respective date (directly from filename)
#name <- substr(files_past[1], 100, 107)                                  # 100 and 107 results from position in path -> adjustment necessary to extract date from path

# Raster als TIFF speichern
#writeRaster(S5P_NO2_final_mask, filename = paste0(out_folder, name, ".tif"), format="GTiff",overwrite=T)


#-------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------
# Add log file

# successful creation and saving of the layer of the respective date
cat("NO2 raster successfully saved on", as.character(Sys.time()), "\n", file = log_con)


#--------------------------------------------------------------------------------------------------------------------------
# 2.2) Rolling Window
#--------------------------------------------------------------------------------------------------------------------------
# The Rolling Window displays the average value of the measured NO2 concentrations over 5 consecutive days and stores them in a TIFF

#--------------------------------------------------------------------------------------------------------------------------
# List of all stored NO2 TIFFS
date_folder <- "C:\\data\\Sentinel5\\2020"                              # adjust manually

# Listing only TIFF files in this path (pattern)
date_list <- list.files(date_folder, pattern= "\\.tif$", full.names=T, all.files=F, include.dirs = F)

#--------------------------------------------------------------------------------------------------------------------------
# Extracting the date from the TIFF name
date <- lubridate::ymd(basename(date_list))

#--------------------------------------------------------------------------------------------------------------------------
# For loop, which always lists five consecutive TIFFs (by date) and calculates and stores the average NO2 concentration as TIFF
for (i in 5:length(date_list)){                                       # starts with day 5 (fifth first day)
  window_dates <- date[i:(i-4)]                                       # includes the first 5 days
  window_files <- lapply(window_dates, function(x){
    date_str <- strftime(x, format="%Y%m%d")
    list.files(date_folder, pattern = date_str, full.names = T)[1]
  })
  
  in_stack <- stack()                                                 # Create an empty stack
  for (j in 1:length(window_files)){                                  # stacks 5 listed scenes
    in_stack <- raster::stack(window_files)
  }
  mean <- stackApply(in_stack, indices =  rep(1,nlayers(in_stack)), fun = "mean", na.rm = T)# Calculation of the mean value
  plot(mean)                                                                                # Plot the new grid with the mean values
  name <- as.character(c(window_dates[1]))                                                  # Naming the layer after the fifth date
  out_folder <- "C:\\data\\Sentinel5\\2020_rolling_window_5\\"                               # Path where the average TIFFs are stored (Adjust manually)
  filename = paste0(out_folder, name, ".tif")
  if (file.exists(filename)){                                                               # Skip the file if it already exists
    next
  }
  writeRaster(mean, filename = paste0(out_folder, name, ".tif"), format="GTiff",overwrite=F)# Save the raster as TIFF
}


#-------------------------------------------------------------------------------------------------------------------------
# add log file

# successful creation and saving of the Rolling Window layer of the respective date
cat("Rolling Window raster successfully saved on", as.character(Sys.time()), "\n", "\n", file = log_con)
 

############################################################################################################################
##################################################################################
############ 3.) Delete all NetDCFs after creating the TIFF files ############
##################################################################################

# Define the folder where all NetCDFs are located
CDF_folder <- "C:\\data\\Sentinel5\\raw\\datasets\\Sentinel-5P"         # adjust manually

# Specify list of all NetCDFs to be deleted
files_del <- list.files(CDF_folder, full.names=T, all.files=F, include.dirs = F)

# deletes all NetCDFs in specified folder
do.call(file.remove, list(files_del))


############################################################################################################################
############################################################################################################################
