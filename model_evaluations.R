# model evaluations
# author: Lauren Jenkins
# 12 April 2022
# last updated: 6 July 2022

rm(list = ls())

library(xlsx)

library(enmSdm)
library(geosphere)
library(sp)
library(dismo)

cluster <- F

if (cluster == T) { # constants for running on cluster
  args <- commandArgs(TRUE)
  
  evalType <- args[1]
  evalType <- as.character(evalType)
  
  gcmList <- args[2]
  gcmList <- unlist(gcmList)
  
  genus <- args[3]
  genus <- as.character(genus)
  
  speciesList <- args[4]
  speciesList <- strsplit(speciesList, split = ', ')
  speciesList <- unlist(speciesList)
  
  baseFolder <- '/mnt/research/TIMBER/PVMvsENM/'
  setwd('/mnt/research/TIMBER/PVMvsENM')
} else {
  evalType <- 'random'
  
  ## genus constants ##
  genus <- 'fraxinus'
  speciesList <- paste0('Fraxinus ', 
                        c('americana', 'caroliniana','cuspidata',
                          'greggii', 'nigra', 'pennsylvanica', 
                          'profunda', 'quadrangulata'))
  
  gcmList <- c('hadley', 'ccsm', 'ecbilt') # general circulation models for env data
  
  baseFolder <- '/Volumes/lj_mac_22/MOBOT/by_genus/'
  setwd(paste0(baseFolder, genus))
}

ll <- c('longitude', 'latitude')
pc <- 5
predictors <- c(paste0('pca', 1:pc))

if(evalType == 'random') {
  for (gcm in gcmList) {
    print(paste0("GCM = ", gcm))
    
    a <- data.frame(c(seq(1:5)))
    c <- data.frame(c(seq(1:5)))
    colnames(a)[1] <- colnames(c)[1] <- 'fold #'
    
    for(sp in speciesList) {
      print(paste0("Species = ", sp))
      
      speciesAb_ <- sub("(.{4})(.*)", "\\1_\\2", 
                        paste0(substr(sp,1,4), toupper(substr(sub("^\\S+\\s+", '', sp),1,1)), 
                               substr(sub("^\\S+\\s+", '', sp),2,4)))
      
      # set constants for retrieving background sites #
      bgFileName <- paste0(baseFolder, 
                           'background_sites/Random Background Sites across Study Region.Rdata')
      load(bgFileName) # load bg sites in calibration region
      
      modelFileName <- paste0('./models/predictions/', speciesAb_, 
                              '/GCM_', gcm, '_PC', pc, '.rData')
      load(modelFileName) # load model object, bg, and records for given species
      
      evalFolderName <- paste0('./models/model_evaluations/', speciesAb_, 
                               '/', evalType, '_k_folds/')
      if(!dir.exists(evalFolderName)) dir.create(evalFolderName, recursive = TRUE, 
                                                 showWarnings = FALSE)
      
      evalFolderName <- paste0(evalFolderName, gcm, '/')
      if(!dir.exists(evalFolderName)) dir.create(evalFolderName, recursive = TRUE, 
                                                 showWarnings = FALSE)
      
      # variable to store auc & cbi output
      auc <- cbi <- rep(NA, 5)
      
      if (evalType == 'random') {
        kPres <- kfold(records, k = 5) # k-folds for presences
        kBg <- kfold(bg, k = 5) # k-folds for backgrounds
        
        ### code for visualizing folds ###
        # plot(range, main = paste0(sp, ', k-fold #1'))
        # points(records$longitude, records$latitude)
        # points(records$longitude[kPres==1],
        #        records$latitude[kPres==1],
        #        bg='red',
        #        pch=21
        # )
        
        # legend('topright',
        #        legend=c('Training presence', 'Test presence'),
        #        pch=c(1, 16),
        #        col=c('black', 'red'),
        #        bg='white',
        #        cex=0.8
        # )
        
        for(i in 1:5) { # for each k-fold
          print(paste0('K-fold ', i, ':'))
          
          # create training data, with presences/absences vector of 0/1 with all points 
          # EXCEPT the ones in the fold
          envData <- rbind(records[kPres != i, predictors], bg[kBg != i, predictors])
          presBg <- c(rep(1, sum(kPres != i)), rep(0, sum(kBg != i)))
          trainData <- cbind(presBg, envData)
          
          model_tune <- enmSdm::trainMaxNet(data = trainData, resp = 'presBg', 
                                            classes = 'lpq', out = c('models', 'tuning'))
          model <- model_tune$models[[1]]
          
          # predict presences & background sites
          predPres <- raster::predict(model, 
                                      newdata = records[kPres == i,],
                                      clamp = F,
                                      type = 'cloglog')
          predBg <- raster::predict(model, 
                                    newdata = bg[kPres == i,],
                                    clamp = F,
                                    type = 'cloglog')
          
          save(model, model_tune, predPres, predBg, kPres, kBg,
               file = paste0(evalFolderName, '/model_', i, '.Rdata'), overwrite = T)
          
        }
      } else if (evalType == 'geo') {
        # create g-folds
        gPres <- geoFold(x = records, k = 5, minIn = 5, minOut = 10, longLat = ll)
        
        # now, we have our folds, but we want to divide the bg sites into 
        # folds based on where they are in relation to records
        
        # initialize vectors to store g-fold assignments
        gTestBg <- rep(NA, nrow(bg))
        
        # convert records to sp object for gDistance function
        sp.records <- records
        coordinates(sp.records) <- ~longitude + latitude
        sp.randomBg <- bg
        coordinates(sp.randomBg) <- ~longitude + latitude
        
        # divide bg sites between training & test
        nearest <- apply(gDistance(sp.records, sp.randomBg, byid = T), 1, which.min)
        
        for (k in 1:nrow(bg)) {
          gTestBg[k] <- gPres[nearest[k]]
        }
        
        for (m in 1:5) { # make training data frame with predictors 
          # and vector of 1/0 for presence/background
          
          envData <- rbind(
            records[gPres!=m, predictors],
            bg[gTestBg!=m, predictors]
          )
          
          presBg <- c(rep(1, sum(gPres!=m)), rep(0, sum(gTestBg!=m)))
          trainData <- cbind(presBg, envData)
          
          # maxent model
          model_tune <- enmSdm::trainMaxNet(data = trainData, resp = 'presBg', 
                                            classes = 'lpq', out = c('models', 'tuning'))
          model <- model_tune$models[[1]]
          
          # predict presences & background sites
          predPres <- raster::predict(model, 
                                      newdata = records[gPres == m,],
                                      clamp = F,
                                      type = 'cloglog')
          predBg <- raster::predict(model, 
                                    newdata = bg[gTestBg == m,],
                                    clamp = F,
                                    type = 'cloglog')
          
          save(model, predPres, predBg, gPres, gTestBg, model_tune, 
               file = paste0(evalFolderName, '/model ', m, '.Rdata'), compress = T)
          
        }
      }
      
      # evaluate
      thisEval <- evaluate(p = as.vector(predPres), a = as.vector(predBg))
      thisAuc <- thisEval@auc
      thisCbi <- contBoyce(pres = predPres, bg = predBg)
      
      # print(paste('AUC = ', round(thisAuc, 2), ' | CBI = ', round(thisCbi, 2)))
      
      auc[m] <- thisAuc
      cbi[m] <- thisCbi
    }
    
    save(auc, cbi, file = paste0(evalFolderName, '/auc_cbi_vals.Rdata'))
    
    a <- cbind(a, auc)
    c <- cbind(c, cbi)
    n <- ncol(a)
    colnames(a)[n] <- colnames(c)[n] <- sp
    
  }
  
  write.xlsx(a, file = paste0('./models/model_evaluations/', evalType, '_evals.xlsx'), 
             sheetName = paste0(gcm, '_auc'), append = T, row.names = F)
  write.xlsx(c, file = paste0('./models/model_evaluations/', evalType, '_evals.xlsx'), 
             sheetName = paste0(gcm, '_cbi'), append = T, row.names = F)
  
  save(a, c, file = paste0('./models/model_evaluations/', gcm, '_evals.Rdata'))
  
}