# ENM Template for 'mega-species' Fraxinus
# Author: Lauren Jenkins
# 25 January 2022
# Last updated: 8 June 2022

rm(list = ls())

## load packages ##

# for visualization # 
library(sp)
library(sf)
library(ggplot2)

# maps #
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires) 

# for handling rasters # 
library(raster)
library(rgdal)
library(dismo)
# library(rgeos) # not used, but may be helpful when debugging

# additional tools #
library(tools)
library(units)
library(dplyr)
library(tidyr)

# enm tools #
library(enmSdm)

setwd('/Volumes/lj_mac_22/MOBOT/by_genus/fraxinus')

# set constants #
nad27_crs <- getCRS('nad27')
wgs84_crs <- getCRS('wgs84')

ll <- c('longitude', 'latitude')

yrs <- seq(22 , 1, by = -1) # timestep indicies helper
climYear <- 0 # only modeling present-day in this script

pc <- 5 # number of principal components used from environmental data

gcmList <- c('hadley', 'ccsm', 'ecbilt') # general circulation models for env data

# reference raster from ccsm env data for cell size, used in thinning
lorenzRast <- raster::raster('/Volumes/lj_mac_22/MOBOT/PVMvsENM/data_and_analyses/env_data/Lorenz/V2/ccsm_21-0k_all_tifs_LJ/0BP/an_avg_ETR.tif')

sink('./cleaning_script_genus.txt')
cat(paste0('GENUS: FRAXINUS\n'))
sp <- 'Fraxinus all'
speciesAb_ <- 'Frax_All'

load('/Volumes/lj_mac_22/MOBOT/PVMvsENM/workspaces/01 - Modeling Workspace - Fraxinus Range Maps (BIEN + Little)')
range <- rgeos::gUnion(rgeos::gUnion(rgeos::gUnion(littleRange_FraxAmer, littleRange_FraxCaro), 
                                     rgeos::gUnion(littleRange_FraxCusp, littleRange_FraxGreg)),
                       rgeos::gUnion(rgeos::gUnion(littleRange_FraxNigr, littleRange_FraxPenn), 
                                     rgeos::gUnion(littleRange_FraxProf, littleRange_FraxQuad)))

if (file.exists(paste0('./species_records/03_', gsub(' ', '_', tolower(sp)), '_final_records.rData'))) { # if already cleaned, load cleaned data set
  print('Already cleaned!\n')
} else { # otherwise, clean the data! #
  
  recordsFileName <- paste0('./species_records/01_', 
                            gsub(' ', '_', tolower(sp)), 
                            '_retained_records.rData')
  
  if (file.exists(recordsFileName)) { # load only retained occurrences
    occs <- load(recordsFileName) 
  } else { 
    ## standard cleaning procedures ##
    # remove records without coordinates or dates,
    # remove records from before 1900,
    # date_collected column is recorded twice, remove one of them,
    # remove imprecise coordinates
    
    # concatenate all (raw) occurrences from all Fraxinus species into one dataframe
    # for a "mega-Fraxinus" species
    # occsRaw #
    speciesList <- paste0('Fraxinus ', 
                          c('americana', 'caroliniana', 'cuspidata', 'greggii', 
                            'nigra', 'pennsylvanica', 'profunda', 'quadrangulata'))
    
    occsRaw <- data.frame()
    for(sp in speciesList) {
      species <- gsub(tolower(sp), pattern=' ', replacement='_')
      speciesFileName <- paste0('./data_and_analyses/species_records/01_', 
                                gsub(tolower(species), pattern = ' ', 
                                     replacement = '_'), '_retained_records.rds')
      
      occs_current <- readRDS(speciesFileName)
      occsRaw <- rbind(occsRaw, occs_current)
    }
    
    # remove records without coordinates or dates to ensure geographic reliability
    occs <- occsRaw[!(is.na(occsRaw$longitude) 
                      | is.na(occsRaw$latitude) 
                      | is.na(occsRaw$date_collected)), ]
    
    # remove records before 1900
    occs$year <- substr(occs$date_collected, 1, 4)
    occs <- occs[occs$year >= 1900, ]
    
    # if date_collected column is duplicated, remove one of them
    occs <- occs[!duplicated(as.list(occs))]
    
    # remove imprecise coordinates
    llOcc <- cbind(occs$longitude, occs$latitude)
    coordPrecis <- coordPrecision(llOcc) # treat lat/long as decimal degrees
    coordPrecisDMS <- coordPrecision(llOcc, dms = TRUE) # treat lat/long as DMS
    precision <- data.frame(llOcc, coordPrecis, coordPrecisDMS)
    
    # visualize the distribution of precision among occurrences
    par(mfrow=c(2,1))
    hist(precision$coordPrecis,
         breaks = 100,
         xlab='Coordinate uncertainty (m)',
         ylab='Number of records',
         main='decimal degrees',
         col='red')
    
    hist(precision$coordPrecisDMS,
         breaks=100,
         xlab='Coordinate uncertainty (m)',
         ylab='Number of records',
         main='dms',
         col='red' )
    
    # find the max between two precision measurements (dms vs decimal degrees)...
    occs$precision <- pmax(coordPrecis, coordPrecisDMS)
    
    # how many are deemed "imprecise"?
    cat(paste0("Number of occurrences with precision > 1,000 m = ", 
               length(which(occs$precision > 1000)), '\n'))
    
    # ...and remove occurrences where precision > 1,000 m
    occs <- occs[which(occs$precision <= 1000), ]
    
    ## visualize occurrences in each stage ##
    # convert to spatial object
    occsSp <- SpatialPointsDataFrame(occs[, ll], data = occs, proj4string = CRS(wgs84_crs))
    
    # plot cleaned occurrences
    par(mfrow = c(1,1))
    plot(range, col = scales::alpha('blue', 0.4), border = 'blue', 
         main = paste0("retained BIEN records w/ Little range map\n", sp))
    plot(occsSp, pch = 16, cex = 0.3, col = 'red', add = TRUE)
    map("state", add = TRUE)
    map("world", add = TRUE)
    
    # save retained records
    save(occs, file = recordsFileName)
  }
  
  ## thin data ##
  thinnedFileName <- paste0('./species_records/02_', 
                            gsub(' ', '_', tolower(sp)), 
                            '_thinned_records.rData')
  
  if (file.exists(thinnedFileName)) { # if already thinned, load that file
    load(thinnedFileName)
  } else { # otherwise, thin the occurrences!
    # eliminate duplicates within each cell, 
    # using Lorenz environmental raster for reference
    thinned <- elimCellDups(occs, lorenzRast, longLat = ll)
    save(thinned, file = thinnedFileName)
  }
  
  # create sf objects for visualization
  speciesSf_thinned <- st_as_sf(x = thinned,
                                coords = c(x = 'longitude',
                                           y = 'latitude'),
                                crs = wgs84_crs)
  speciesSf_occs <- st_as_sf(x = occs,
                             coords = c(x = 'longitude',
                                        y = 'latitude'),
                             crs = wgs84_crs)
  
  # load country/world basemap data for sf plotting
  world <- ne_countries(returnclass = 'sf')
  states <- ne_states(country = 'united states of america', returnclass = 'sf')
  canada <- ne_states(country = 'canada', returnclass = 'sf')
  
  ggplot(data = world) +
    theme_bw() +
    geom_sf(fill = "white") +
    geom_sf(data = states, fill = NA) + 
    geom_sf(data = canada, fill = NA) + 
    geom_sf(data = speciesSf_thinned, col = scales::alpha('red'), cex = 0.5) +
    coord_sf(xlim = c(min(occs$longitude), max(occs$longitude)), 
             ylim = c(min(occs$latitude), max(occs$latitude)), expand = TRUE) +
    xlab("Longitude") + 
    ylab("Latitude") +
    ggtitle(paste0(sp, ' occurences'), subtitle = "(Cleaned and cell duplicates eliminated)") +
    theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                          size = 0.5), panel.background = element_rect(fill = "lavender"))
  
  # create a buffer surrounding the Little range map to exclude occurrence data that
  # is unrealistically out of range
  bufferFileName <- paste0('./species_records/buffer/', 
                           gsub(' ', '_', tolower(sp)), 
                           '_buffer.rData')
  
  if(!file.exists(bufferFileName)) { # if buffer hasn't been calculated, calculate it
    
    # create a buffer that will include 97% (0.03 * # of occurrences) of the occurrences
    # this threshold value can change based on how much you want to filter "out of bounds" data
    threshold <- ceiling(0.03 * nrow(speciesSf_thinned))
    
    # convert maps to albers NA for buffer calculation
    rangeMapAlb <- st_transform(st_as_sf(x = range), getCRS('albersNA'))
    speciesSf_thinnedAlb <- st_transform(speciesSf_thinned, getCRS('albersNA'))
    speciesSf_occsAlb <- st_transform(speciesSf_occs, getCRS('albersNA'))
    
    # function for calculating the proper buffer distance to include desired % of occurrences
    calculate_buffer <- function(t) {
      while (TRUE) {
        t <- t + 10 # increase threshold value by 10 each iteration
        buffer_distance_temp <- as_units(t, "km") # convert to km units
        range_buffer <- st_buffer(rangeMapAlb, dist = buffer_distance_temp)
        
        # only keeped occurrences within buffered distance
        speciesSf_thinned_buffered <- speciesSf_thinnedAlb %>%
          mutate(within_range = lengths(st_within(x = speciesSf_thinnedAlb,
                                                  y = rangeMapAlb)),
                 within_buffer = lengths(st_within(x = speciesSf_thinnedAlb,
                                                   y = range_buffer)))
        
        ## sanity checking ... ##
        ## uncomment the following lines to see how the # of occurrences outside
        ## the buffer decreases as the buffer distance increases
        
        cat(paste0("Buffer distance = ", t, " km\n"))
        cat(paste0("Number of occurrences outside buffer = ",
                   length(which(speciesSf_thinned_buffered$within_buffer == 0)), '\n'))
        
        # if buffer has reached the thresholded value, stop increasing the buffer
        # & save it here
        if(length(which(speciesSf_thinned_buffered$within_buffer == 0)) < threshold) {
          if (t < 200) { # to accommodate for species with very few records,
            # if the calculated buffer is < 200 km,
            # double it to ensure enough records are included in model
            
            cat(paste0('Yes, buffer distance is < 200 (', t, ' km), so doubling buffer distance...\n'))
            buffer_distance <- as_units((t*2), 'km') # double buffer distance
            
            cat(paste0('New buffer distance = ', buffer_distance, ' km\n'))
            range_buffer <- st_buffer(rangeMapAlb, dist = buffer_distance)
            
            # recalculate buffers with new buffer distance
            speciesSf_thinned_buffered <<- speciesSf_thinnedAlb %>%
              mutate(within_range = lengths(st_within(x = speciesSf_thinnedAlb,
                                                      y = rangeMapAlb)),
                     within_buffer = lengths(st_within(x = speciesSf_thinnedAlb,
                                                       y = range_buffer)))
            speciesSf_occs_buffered <<- speciesSf_occsAlb %>% 
              mutate(within_range = lengths(st_within(x = speciesSf_occsAlb,
                                                      y = rangeMapAlb)),
                     within_buffer = lengths(st_within(x = speciesSf_occsAlb,
                                                       y = range_buffer)))
            
            # save buffer data to use later
            save(buffer_distance, range_buffer, speciesSf_thinned_buffered, 
                 speciesSf_occs_buffered, file = bufferFileName)
          } else { # if the calculated buffer is larger than 200 km, use that distance 
            buffer_distance <<- buffer_distance_temp
            speciesSf_occs_buffered <<- speciesSf_occsAlb %>% 
              mutate(within_range = lengths(st_within(x = speciesSf_occsAlb,
                                                      y = rangeMapAlb)),
                     within_buffer = lengths(st_within(x = speciesSf_occsAlb,
                                                       y = range_buffer)))
            
            # save buffer data to use later
            save(buffer_distance, range_buffer, speciesSf_thinned_buffered, 
                 speciesSf_occs_buffered, file = bufferFileName)
          }
          break
        }
      }
    }
    
    t <- 10 # start buffer distance at 10 km
    
    # calculate buffer
    calculate_buffer(t)
  }
  
  load(bufferFileName) # load in buffer
  
  #  final dataset with occurrences outside buffer removed
  speciesSf_final <- speciesSf_thinned %>%
    mutate(within_range = lengths(st_within(x = speciesSf_thinnedAlb,
                                            y = rangeMapAlb)),
           within_buffer = lengths(st_within(x = speciesSf_thinnedAlb,
                                             y = range_buffer))) %>%
    filter(within_buffer > 0)
  
  # buffers create circles around the edges of the range map
  # since the range map isn't one succinct polygon, there are serveral overlapping circles
  # so we need to unionize them to create one all-encompassing buffer
  range_buffer_final <- st_union(range_buffer)
  
  finalRecordsFileName <- paste0('./species_records/03_', 
                                 gsub(' ', '_', tolower(sp)), 
                                 '_final_records.rData')
  
  save(speciesSf_final, range_buffer_final, file = finalRecordsFileName)
  
  ## visualize data after buffering ##
  
  # cleaned & thinned occurrences with calculated buffer
  # occurrences that fall outside the range are colored red
  ggplot(data = world) +
    theme_bw() +
    geom_sf(fill = "white") +
    geom_sf(data = states, fill = NA) + 
    geom_sf(data = canada, fill = NA) +
    geom_sf(data = rangeMapAlb, 
            color = 'darkgreen',
            fill = 'green', 
            alpha = 0.4, 
            inherit.aes = FALSE) +
    geom_sf(data = range_buffer_final,
            color = 'yellow',
            fill = NA,
            inherit.aes = FALSE) +
    geom_sf(data = speciesSf_thinned_buffered,
            size = 0.75, 
            aes(color = within_buffer > 0),
            inherit.aes = FALSE) +
    scale_colour_manual(values = setNames(c('black','red'),
                                          c(T, F)), 
                        guide = "none") +
    guides(fill = "none") +
    coord_sf(xlim = c(min(thinned$longitude) - 5, 
                      max(thinned$longitude) + 5), 
             ylim = c(min(thinned$latitude), 
                      max(thinned$latitude) + 5), 
             expand = TRUE) +
    xlab("Longitude") + 
    ylab("Latitude") +
    ggtitle(paste0(sp, ' occurences'), 
            subtitle = paste0("(Cleaned & thinned, buffer = ",
                              buffer_distance, 
                              " km)")) +
    theme(panel.grid.major = element_line(color = gray(0.5), 
                                          linetype = "dashed", 
                                          size = 0.5), 
          panel.background = element_rect(fill = "lavender"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # cleaned & thinned occurrences with calculated buffer
  # records that fell outside buffer are removed
  # red points are occurrences that fall outside of the Little range map
  ggplot(data = world) +
    theme_bw() +
    geom_sf(fill = "white") +
    geom_sf(data = states, fill = NA) + 
    geom_sf(data = canada, fill = NA) +
    geom_sf(data = rangeMapAlb, 
            color = 'darkgreen',
            fill = 'green', 
            alpha = 0.4, 
            inherit.aes = FALSE) +
    geom_sf(data = range_buffer_final,
            color = 'yellow',
            fill = NA,
            inherit.aes = FALSE) +
    geom_sf(data = speciesSf_final,
            size = 0.5,
            aes(color = within_range > 0),
            inherit.aes = FALSE) +
    scale_colour_manual(values = setNames(c('black','red'),
                                          c(T, F)), 
                        guide = "none") +
    guides(fill = "none") +
    coord_sf(xlim = c(min(thinned$longitude) - 5, 
                      max(thinned$longitude) + 5), 
             ylim = c(min(thinned$latitude), 
                      max(thinned$latitude) + 5), 
             expand = TRUE) +
    xlab("Longitude") +
    ylab("Latitude") +
    ggtitle(paste0(sp, ' occurences'), 
            subtitle = paste0("(Cleaned, thinned, and duplicates eliminated, with buffer = ",
                              buffer_distance, 
                              " km)")) +
    theme(panel.grid.major = element_line(color = gray(0.5), 
                                          linetype = "dashed", 
                                          size = 0.5), 
          panel.background = element_rect(fill = "lavender"),
          plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  # How many occurrences are in the final dataset?
  cat(paste0('Final number of occurrences = ', nrow(speciesSf_final), '\n'))
  
  #   How many occurrences fell outside of the buffer?
  cat(paste0("Number of occurrences outside buffer & excluded = ", 
             length(which(speciesSf_thinned_buffered$within_buffer == 0)), '\n'))
  
}
sink()
sink('./ENM_script_genus.txt')
for (gcm in gcmList) {
  cat(paste0('\nGCM = ', gcm, ', Species: ', sp, '\n'))
  
  # identify study region
  studyRegionFileName <- '/Volumes/lj_mac_22/pollen/predictions-FRAXINUS_meanpred_iceMask.tif'
  studyRegionRasts <- brick(studyRegionFileName)
  
  # load climate for given gcm 
  load(paste0('/Volumes/lj_mac_22/MOBOT/by_genus/env_data/', gcm, 
              '/PCA_climate_rasters_pc', pc, '.Rdata'))
  load(paste0('/Volumes/lj_mac_22/MOBOT/by_genus/env_data/', gcm, 
              '/PCA_prcomp_pc', pc, '.Rdata'))
  
  # set constants for given gcm
  fileName <- paste0('/Volumes/lj_mac_22/MOBOT/by_genus/env_data/', gcm, 
                     '/clipped_data/envDataClipped_', climYear, 'YBP_pc', pc, '.tif')
  vars <- names(clim[[1]])
  
  # set constants for Lorenz gcm's 
  if (gcm == 'ccsm' | gcm == 'ecbilt') {
    clim <- lorenz
    workingFolder <- paste0('/Volumes/lj_mac_22/MOBOT/by_genus/env_data/', gcm, 
                            '/tifs/', climYear, 'BP')
    vars <- names(clim[[1]])
  }
  
  ## load environment data from PCA ##
  if (file.exists(fileName)) {  # if already clipped to respective year, load that 
    envData <- brick(fileName)
    names(envData) <- paste0('pca', 1:pc) # rename raster layers to match pc
  } else { # otherwise, clip data to correct climate year
    pcPrediction <- list()
    
    for (i in 1:length(clim)) { # label each env layer by variable name and year
      names(clim[[i]]) <- vars
      pcPrediction[i] <- raster::predict(clim[[i]], pca, index = 1:pc)
      names(pcPrediction[[i]]) <- paste0("pc", 1:pc, "_", (i-1)*1000, "KYBP")
    }
    
    envDataPca <- stack(pcPrediction)
    
    envYr <- pcPrediction[[(climYear/1000) + 1]] # keep only the PCA rasters for given climate year
    names(envYr) <- paste0('pca', 1:pc) # label layers by pc
    
    ## sanity check for ensuring the rasters are in the correct projection
    # print("Ensure that the projection of these rasters is WGS84:")
    # print(paste0("Projection of envYr = ", projection(envYr)))
    
    envDataClipped <- list()
    
    for (n in 1:nlayers(envYr)) { # clip PCAs to study extent for given species
      x <- envYr[[n]]
      x <- crop(x, extent(studyRegionRasts[[n]]))
      projection(x) <- getCRS("WGS84")
      envDataClipped[[n]] <- x
    }
    
    envData <- stack(envDataClipped)
    
    # plot(envData) # plot clipped environmental pca rasters
    writeRaster(envData, fileName, format = 'GTiff', overwrite = T)
  } 
  
  recordsFileName <- paste0('./species_records/03_', 
                            gsub(' ', '_', tolower(sp)), 
                            '_final_records.rData')
  load(recordsFileName) # load records for given species
  
  ## prepare occurrence data for maxent ##
  records <- data.frame(speciesSf_final)
  records$geometry <- gsub("[c()]", "", records$geometry) # clean ll format in geometry column
  
  # create separate columns for lat & long
  records <- separate(data = records, 
                      col = 'geometry', 
                      into = ll, 
                      sep = "\\,") 
  
  records$longitude <- as.double(records$longitude)
  records$latitude <- as.double(records$latitude)
  
  # extract environmental data at each occurrence point
  occsEnv <- raster::extract(envData, 
                             cbind(records$longitude, 
                                   records$latitude))
  occsEnvDf <- as.data.frame(occsEnv) # convert to dataframe
  records <- cbind(records, occsEnvDf) # add to records dataframe
  
  ## remove any records that fall in the water ##
  if (exists('water')) rm(water) # remove 'water' from previous species
  if (any(is.na(rowSums(occsEnvDf)))) { # define points in the water
    water <- records[which(is.na(rowSums(occsEnvDf))), ] 
    water <- SpatialPointsDataFrame(water[,ll], data = water,
                                    proj4 = getCRS('wgs84', TRUE))
  }
  
  if (any(is.na(rowSums(occsEnvDf)))) records <- records[-which(is.na(rowSums(occsEnvDf))), ] # remove records in water
  
  # convert to sp object for visualization
  recordsSp <- SpatialPointsDataFrame(records[, ll], data = records,
                                      proj4 = getCRS('wgs84', TRUE))
  
  # visualize points that fall in the water (colored in blue)
  plot(recordsSp, pch = 16, cex = 0.5, col = "red", 
       main = paste0(sp, ' occurrences (BIEN) thinned'))
  if (exists("water")) {
    plot(water, col = 'blue', add = TRUE)
  } 
  map("state", add = TRUE)
  map("world", add = TRUE)
  
  # save.image(paste0('./workspaces/04 - Modeling Workspace - Clipping ', 
  #                   sp, '_PC_', pc, '_GCM_', gcm))
  
  bufferFileName <- paste0('./species_records/buffer/', 
                           gsub(' ', '_', tolower(sp)), 
                           '_buffer.rData')
  load(bufferFileName)
  
  ## calculate calibration region at 320-km to extract bg sites from ##
  # draw from all of NA #
  calibBuffer <- st_buffer(st_transform(st_as_sf(x = recordsSp), getCRS('albersNA')),
                           dist = as_units(320, 'km'))
  calibBuffer <- st_union(calibBuffer) # unionize
  
  # convert to different crs objects 
  calibRegionSpAlb <- sp::spTransform(as(calibBuffer, 'Spatial'), getCRS('albersNA', TRUE))
  calibRegionSpWgs <- sp::spTransform(calibRegionSpAlb, getCRS('wgs84', TRUE))
  
  # set constants for retrieving background sites #
  bgFileName <- './background_sites/Random Background Sites across Study Region.Rdata'
  
  # load bg sites in calibration region if they have already been defined (bgTestSp, bgCalib, bgEnv, bg)
  if (file.exists(bgFileName)) load(bgFileName) else { 
    # otherwise, get 20,000 random background sites from calibration region
    bgTestSpAlb <- suppressWarnings(sp::spsample(calibRegionSpAlb, n=20000, 
                                                 type='random', iter = 10))
    
    bgTestSp <- sp::spTransform(bgTestSpAlb, wgs84_crs) # transform to wgs84
    bgCalib <- as.data.frame(coordinates(bgTestSp)) # lat/long of background sites
    names(bgCalib) <- ll
    
    save(bgTestSp, bgCalib, file = bgFileName, compress = T, overwrite = T)
  }
  
  # plot the bg sites to verify
  plot(bgTestSp, pch = 16, cex = 0.5, col = "red", 
       main = paste0(sp, ' background sites'))
  plot(calibRegionSpWgs, add = TRUE, border = 'blue')
  map("state", add = TRUE)
  map("world", add = TRUE)
  
  climate <- envData
  bgEnv <- raster::extract(climate, bgCalib) # extract environment at random background sites
  bgEnv <- as.data.frame(bgEnv) # convert to dataframe
  
  # remove any sites with NA for at least one variable #
  isNa <- is.na(rowSums(bgEnv))
  if (any(isNa)) {
    bgCalib <- bgCalib[-which(isNa), ]
    bgEnv <- bgEnv[-which(isNa), ]
  }
  
  bg <- cbind(bgCalib, bgEnv) # combine with coordinates
  names(bg)[1:2] <- ll # rename lat/long columns, respectively
  
  presBg <- c(rep(1, nrow(records)), rep(0, nrow(bg))) # identify presences
  occsEnv <- occsEnv[complete.cases(occsEnv), ] # remove NA values
  
  ## prepare env data frame for maxent ##
  env <- rbind(occsEnv, bgEnv)
  env <- cbind(presBg, env)
  env <- as.data.frame(env)
  
  env <- env[complete.cases(env), ] # remove NA values
  
  ## run maxent for species ##
  # model tuning for easy fine-tuning later
  envModel_tune <- enmSdm::trainMaxNet(data = env, resp = 'presBg', 
                                       classes = 'lpq', out = c('models', 'tuning'))
  envModel <- envModel_tune$models[[1]] # select best fitted model
  
  predictors <- c(paste0('pca', 1:pc))
  
  print('Line 587 = .XML ??')
  # prediction for given year
  envMap <- predict(
    climate[[predictors]], 
    envModel,
    filename = paste0('./models/predictions/', speciesAb_, '/GCM_', gcm,
                      '_PC', pc, '_', climYear, 'ybp'),
    clamp = F,
    format='GTiff',
    overwrite = T,
    type='cloglog')
  
  envMapSp <- rasterToPolygons(envMap) # convert to spatial object for plotting
  
  plot(range, border = 'blue', main = paste0('Maxent output, ', sp,
                                                ' occurrences'))
  plot(envMap, add = TRUE)
  plot(range, border = 'blue', add = TRUE)
  map("state", add = TRUE)
  map("world", add = TRUE)
  points(records$longitude, records$latitude, pch = 16, cex = 0.6, col = 'red')
  
  plot(envMap, main = paste0('Maxent output, ', 
                             sp,
                             ' occurrences'))
  plot(range, border = 'blue', add = TRUE)
  
  modelFileName <- paste0('./models/', speciesAb_, '_Maxent_PC', 
                          pc, '_GCM_', gcm, '.Rdata')
  save(envModel, file = modelFileName, compress = T, overwrite = T) # save model
  
  outputFileName <- paste0('./models/predictions/', speciesAb_, 
                           '/GCM_', gcm, '_PC', pc, '.rData')
  save(range, envMap, envModel, records, file = outputFileName, overwrite = T)
}
sink()

