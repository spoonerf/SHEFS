
##########################################################################
# SHEFS Species distribution model
##########################################################################

library(dismo)
library(rgdal)  
library(raster) 
library(sp)
library(randomForest)
library(kernlab)
library(rJava)
library(doParallel)
library(here)



memory.limit(size=50000000)
#-----------------------------------------------------------------------------
# user section to set directories and chose input data
#-----------------------------------------------------------------------------

# data source
path<-here("Part1/")
func_group<- c('bushmeat','medical','pollinator')

# Input/output directories
predictors_dir <- paste(path,"/Part1_predictors/",func_group[1],"_predictors",sep='')
predictors_LU_dir <- paste(path,"/Part1_predictors/MODIS/",sep='')
occurrence_dir <- paste(path,"/Part1_occurrence/",func_group[1],"/",func_group[1],"_occurrence_in/",sep='')
background_file <- paste(occurrence_dir,func_group[1],"_background/new_",func_group[1],"_background.csv",sep='')

if(!dir.exists(paste(path,"/Part1_occurrence/",func_group[1],"/",func_group[1],"_output/",func_group[1],"_output_SDM/",sep=''))){
  dir.create(paste(path,"/Part1_occurrence/",func_group[1],"/",func_group[1],"_output/",func_group[1],"_output_SDM/",sep=''), recursive = T)
}
outdir_SDM <- paste(path,"/Part1_occurrence/",func_group[1],"/",func_group[1],"_output/",func_group[1],"_output_SDM/",sep='')

if(!dir.exists(paste(path,"/Part1_occurrence/",func_group[1],"/",func_group[1],"_output/",func_group[1],"_output_LU/",sep=''))){
  dir.create(paste(path,"/Part1_occurrence/",func_group[1],"/",func_group[1],"_output/",func_group[1],"_output_LU/",sep=''), recursive = T)
}
outdir_LU <- paste(path,"/Part1_occurrence/",func_group[1],"/",func_group[1],"_output/",func_group[1],"_output_LU/",sep='')

#-----------------------------------------------------------------------------
# data preparation
#-----------------------------------------------------------------------------

# load predictors and occurrence records
predictors<- list.files(predictors_dir, pattern = 'tif')
predictors_in<-paste(predictors_dir, predictors, sep="/")

predictors_r <- lapply(predictors_in, raster)
predictors <- stack(predictors_r)
print(paste('predictors loaded'))

species <- as.character(list.files(occurrence_dir, pattern = 'csv'))


#for (s in 1:length(species)){}
s=15
  occurrence_in <- read.csv(paste(occurrence_dir, species[s], sep=''))
  species_name <- substr(species[s],14,nchar(species[s])-4)   #this could be more reproducible although need to think of a way
  coordinates(occurrence_in) <- ~lon+lat
  crs(occurrence_in) <- "+proj=longlat +datum=WGS84 +no_defs"
  
  ext <- extent(occurrence_in)*1.2
  print(paste(species_name,' occurrence records loaded'))

  # create data frame with presence/background data and predictor values
  # occurence training and test set
  occurrence_lonlat=cbind.data.frame(occurrence_in$lon,occurrence_in$lat) 
  colnames(occurrence_lonlat) = c('lon', 'lat')
  
  #removing points which have missing environmental data
  envpres <- data.frame(extract(predictors, occurrence_lonlat))
  occurrence_lonlat<-occurrence_lonlat[-which(!complete.cases(envpres)),]
  envpres<-envpres[complete.cases(envpres),]
  
  set.seed(0)
  k<-5
  group_pres <- kfold(occurrence_lonlat, k=k)
  # pres_train <- occurrence_lonlat[group_pres != 1, ] 
  # pres_test <- occurrence_lonlat[group_pres == 1, ] 
  print(paste(species_name,'occurrence training and test set done'))

  # Background training and test set (bias correction)
  background_in <- read.csv(paste(background_file)) 
  coordinates(background_in) <- ~lon+lat
  crs(background_in) <- "+proj=longlat +datum=WGS84 +no_defs"
  background_raster <- rasterize(background_in, predictors[[1]], 'X', fun = min)
  background_raster_sp<-crop(background_raster, ext) ##added this so all background points are within the same extent as the models
  set.seed(0)
  backgr <- randomPoints(background_raster_sp, 1000) #removed !is.na as the function automatically excludes NA points and '!is.na' was cancelling it out somehow!
  colnames(backgr) = c('lon', 'lat')
  set.seed(0)
  group_bg <- kfold(backgr, k=k)
  # backg_train <- backgr[group_bg != 1, ]
  # backg_test <- backgr[group_bg == 1, ]
  print(paste(species_name,' background training and test set done'))

  # extract predictors training
  #train <- rbind(pres_train, backg_train)
  #set_frame_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
  envbackg <- data.frame(extract(predictors, backgr))
  #envtrainbackg <- data.frame(cbind(pa=set_frame_train, envtrain))
  print(paste(species_name,' predictors for presence and background extracted'))

  # # extract predictors test
  # #test <- rbind(pres_test, backg_test)
  # envtestpres <- data.frame( extract(predictors, pres_test) )
  # envtestbackg <- data.frame( extract(predictors, backg_test) )
  # print(paste(species_name,' test predictors extracted'))
  # 
  
  envpres_pa<-cbind(1, envpres)
  envbackg_pa<-cbind(0, envbackg)
  colnames(envpres_pa)[1]<-"pa"
  colnames(envbackg_pa)[1]<-"pa"
  env_all<-rbind(envpres_pa, envbackg_pa)
  
      #env_all<-rbind(envtrain, envtest)
    
    #-----------------------------------------------------------------------------
    # SDM Model ensemble
    #-----------------------------------------------------------------------------
    
      # BIOCLIM
      print(paste(' '))
      print(paste('+++++++++++++++'))
      print(paste('+++ Bioclim +++'))
      print(paste('+++++++++++++++'))
      
      
      evl_bc<- list()
      
      cl <- makeCluster(detectCores(), type='PSOCK')
      registerDoParallel(cl)
      clusterExport(cl,c('bioclim', 'evaluate', 'evl_bc', 'threshold', 'predict'))
  
        eval_bc<-foreach(i = 1:k, .combine = 'c') %dopar% {
              pres_train<-occurrence_lonlat[group_pres!=i,]
              pres_test<-occurrence_lonlat[group_pres==i,]
              backg_test<-backgr[group_bg==i,]
              bc <- bioclim(predictors,pres_train)
              evl_bc[[i]] <- dismo:::evaluate(pres_test, backg_test, bc,predictors,type="response")
              
              saveRDS(evl_bc, file = paste(outdir_SDM,species_name,"_eval_bc_",i,".ascii",sep=''),ascii=TRUE)
        
              tr_bc <- threshold(evl_bc[[i]], 'spec_sens')
              saveRDS(tr_bc, file = paste(outdir_SDM,species_name,"_tr_bc_",i, ".ascii",sep=''),ascii=TRUE)
        
              predict_bioclim <- predict(predictors, bc, ext=ext, progress='')
              saveRDS( predict_bioclim, file = paste(outdir_SDM,species_name,"_predict_bioclim_raw_",i,".ascii",sep=''),ascii=TRUE)
        
              predict_bioclim_pa <- predict_bioclim > tr_bc
              saveRDS( predict_bioclim_pa, file = paste(outdir_SDM,species_name,"_predict_bioclim_pa_",i,".ascii",sep=''),ascii=TRUE)
              
              print(evl_bc[[i]])
         }
      stopCluster(cl)
      
      
      auc_bc <- sapply(eval_bc, function(x){slot(x, "auc")} )
      print(auc_bc)
      bc_auc<-mean(auc_bc)
      
      # GENERALIZED LINEAR MODEL
      print(paste(' '))
      print(paste('+++++++++++'))
      print(paste('+++ GLM +++'))
      print(paste('+++++++++++'))
      
      model_glm1=NULL
      vars = names(predictors)
      
      for(i in 1:length(vars)){
        xx = combn(vars,i)
        if(is.null(dim(xx))){
          fla = paste("pa ~", paste(xx, collapse="+"))
          model_glm1[[length(model_glm1)+1]]=glm(as.formula(fla),family=binomial(link = "logit"), data=env_all)
        } else {
          for(j in 1:dim(xx)[2]){
            fla = paste("pa ~", paste(xx[1:dim(xx)[1],j], collapse="+"))
            model_glm1[[length(model_glm1)+1]]=glm(as.formula(fla),family=binomial(link = "logit"), data=env_all) 
          }
        }
      }
    
        AICs = NULL
        for(i in 1:length(model_glm1)){
          AICs[i] = AIC(model_glm1[[i]])
        }
    
        model_glm <- which(AICs==min(AICs))   
        
        
        
        evl_glm<- list()
        
        cl <- makeCluster((detectCores()-2), type='PSOCK')
        registerDoParallel(cl)
        clusterExport(cl,c('bioclim', 'evaluate', 'threshold', 'predict'))
        
        eval_glm<-foreach(i = 1:k, .combine = 'c') %dopar% {
          pres_train<-envpres_pa[group_pres!=i ,]
          pres_test<-envpres_pa[(group_pres==i) ,]
          back_test<-envbackg_pa[(group_bg==i),]
          back_train<-envbackg_pa[(group_bg!=i),]
          envtrain<-rbind(pres_train, back_train)
          
          yy<-names(coef(model_glm1[[model_glm]]))[-1]
          fla<-paste("pa ~", paste(yy, collapse="+"))
          model_glm_out<-glm(as.formula(fla),family=binomial(link = "logit"), data=envtrain)
          evl_glm[[i]] <- evaluate(pres_test, back_test,model= model_glm_out,type="response")
          saveRDS(evl_glm[[i]], file = paste(outdir_SDM,species_name,"_eval_glm",i,".ascii",sep=''),ascii=TRUE)
          
          tr_glm <- threshold(evl_glm[[i]], 'spec_sens')
          saveRDS(tr_glm, file = paste(outdir_SDM,species_name,"_tr_glm",i,".ascii",sep=''),ascii=TRUE)
          
          predict_glm <- predict(predictors, model_glm1[[model_glm]], ext=ext)   #adding in conversion back from logit space
          predict_glm <-raster:::calc(predict_glm, fun=function(x){ exp(x)/(1+exp(x))}) 
          saveRDS( predict_glm, file = paste(outdir_SDM,species_name,"_predict_glm_raw",i,".ascii",sep=''),ascii=TRUE)
          
          predict_glm_pa <- predict_glm > tr_glm
          saveRDS(predict_glm_pa, file = paste(outdir_SDM,species_name,"_predict_glm_pa",i,".ascii",sep=''),ascii=TRUE)
          
          print(evl_glm[[i]])
          }
        
        stopCluster(cl)
        
        auc_gam <- sapply(eval_glm, function(x){slot(x, "auc")} )
        print(auc_gam)
        
        gam_auc<-mean(auc_gam)
        
        
         print(paste('GLM done'))
      
      # MAXENT
      print(paste(' '))
      print(paste('++++++++++++++'))
      print(paste('+++ MAXENT +++'))
      print(paste('++++++++++++++'))
      jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
    
      
      eval_ma<- list()
      for (i in 1:k){
        pres_train<-occurrence_lonlat[group_pres!=i ,]
        pres_test<-occurrence_lonlat[(group_pres==i) ,]
        back_test<-backgr[(group_bg==i),]
        back_train<-backgr[(group_bg!=i),]
        model_ma <- maxent(predictors, pres_train)#,l1_regularizer=0.7)
        saveRDS(model_ma, file = paste(outdir_SDM,species_name,"_model_ma_",i,".ascii",sep=''),ascii=TRUE)
        
        eval_ma[[i]] <- evaluate(pres_test, back_test,model= model_ma,x = predictors)
        
        saveRDS(eval_ma[[i]], file = paste(outdir_SDM,species_name,"_eval_ma",i,".ascii",sep=''),ascii=TRUE)
        tr_ma <- threshold(eval_ma[[i]], 'spec_sens')
        saveRDS(tr_ma, file = paste(outdir_SDM,species_name,"_tr_ma",i,".ascii",sep=''),ascii=TRUE)
        predict_maxent <- predict(model_ma, predictors,ext=ext)  
        saveRDS(predict_maxent, file = paste(outdir_SDM,species_name,"_predict_maxent_raw_",i,".ascii",sep=''),ascii=TRUE)
        plot(predict_maxent)
        predict_maxent_pa <- predict_maxent >tr_ma
        plot(predict_maxent_pa)
        saveRDS(predict_maxent_pa, file = paste(outdir_SDM,species_name,"_predict_maxent_pa",i,".ascii",sep=''),ascii=TRUE)
        print(i)
      }
      
      auc_ma <- sapply( eval_ma, function(x){slot(x, "auc")} )
      print(auc_ma)
      
      ma_auc<-mean(auc_ma)
      
      # 
      # RANDOM FOREST
      print(paste(' '))
      print(paste('+++++++++++++++++++++'))
      print(paste('+++ Random forest +++'))
      print(paste('+++++++++++++++++++++'))
      
      model_rf1 <- paste("pa ~", paste(names(predictors), collapse=" + "))
     
      
      eval_rf<- list()
      for (i in 1:k){
        pres_train<-envpres_pa[group_pres!=i ,]
        pres_test<-envpres_pa[(group_pres==i) ,]
        back_test<-envbackg_pa[(group_bg==i),]
        back_train<-envbackg_pa[(group_bg!=i),]
        envtrain<-rbind(pres_train, back_train)
        model_rf <- randomForest(as.formula(model_rf1), data=envtrain, na.action=na.roughfix) 
        saveRDS(model_rf, file = paste(outdir_SDM,species_name,"_model_rf.ascii",sep=''),ascii=TRUE)
        
        eval_rf[[i]] <- evaluate(pres_test, back_test,model= model_rf,type="response")
        
        saveRDS(eval_rf[[i]], file = paste(outdir_SDM,species_name,"_eval_rf",i,".ascii",sep=''),ascii=TRUE)
        tr_rf <- threshold(eval_rf[[i]], 'spec_sens')
        saveRDS(tr_rf, file = paste(outdir_SDM,species_name,"_tr_rf",i,".ascii",sep=''),ascii=TRUE)
        predict_rf <- predict(model_rf, as.data.frame(predictors),ext=ext)  
        saveRDS(predict_rf, file = paste(outdir_SDM,species_name,"_predict_rf_raw_",i,".ascii",sep=''),ascii=TRUE)
        plot(predict_rf)
        predict_rf_pa <- predict_rf >tr_rf
        plot(predict_rf_pa)
        saveRDS(predict_rf_pa, file = paste(outdir_SDM,species_name,"_predict_rf_pa",i,".ascii",sep=''),ascii=TRUE)
        
        print(i)
      }
      
      auc_rf <- sapply(eval_rf, function(x){slot(x, "auc")} )
      print(auc_rf)
      
      rf_auc<-mean(auc_rf)
      
      
      model_rf <- randomForest(as.formula(model_rf1), data=envtrain, na.action=na.roughfix) 
      saveRDS(model_rf, file = paste(outdir_SDM,species_name,"_model_rf.ascii",sep=''),ascii=TRUE)
  
      eval_rf <- evaluate(envtestpres, envtestbackg, model_rf)
      saveRDS(eval_rf, file = paste(outdir_SDM,species_name,"_eval_rf.ascii",sep=''),ascii=TRUE)
      print(eval_rf)
    
      tr_rf <- threshold(eval_rf, 'spec_sens')
      saveRDS(tr_rf, file = paste(outdir_SDM,species_name,"_tr_rf.ascii",sep=''),ascii=TRUE)
      print(tr_rf)
    
      predict_rf <- predict(predictors, model_rf, ext=ext)  
      saveRDS(predict_rf, file = paste(outdir_SDM,species_name,"_predict_rf_raw.ascii",sep=''),ascii=TRUE)
      plot(predict_rf)
      print(predict_rf)
    
      predict_rf_pa <- predict_rf > tr_rf
      plot(predict_rf_pa)
      saveRDS(predict_rf_pa, file = paste(outdir_SDM,species_name,"_predict_rf_pa.ascii",sep=''),ascii=TRUE)

      plot(predict_rf)
      points(occurrence_in$lon, occurrence_in$lat)
      print(paste('Random forest done'))
      
      # MODEL ENSEMBLE
      print(paste(' '))
      print(paste('++++++++++++++++++++++'))
      print(paste('+++ Model ensemble +++'))
      print(paste('++++++++++++++++++++++'))

      ##need to run the models again on the whole data set - rather than just the training 
      
      model_bc_all<-bioclim(predictors, occurrence_lonlat)
      
      yy<-names(coef(model_glm1[[model_glm]]))[-1]
      fla<-paste("pa ~", paste(yy, collapse="+"))
      model_glm_all<-glm(as.formula(fla),family=binomial(link = "logit"), data=env_all)
      
      model_ma_all <- maxent(predictors, occurrence_lonlat)
      model_rf_all <- randomForest(as.formula(model_rf1), data=env_all, na.action=na.roughfix) 
      
      predict_bioclim_all <- predict(predictors, model_bc_all,ext=ext, progress='') 
      predict_glm_all<- predict(predictors, model_glm_all, ext = ext)
      predict_glm_all <-raster:::calc(predict_glm_all, fun=function(x){ exp(x)/(1+exp(x))})
      predict_maxent_all<-predict(model_ma_all, predictors,ext=ext)  
      predict_rf_all<-predict(model_rf_all, as.data.frame(predictors),ext=ext)  
      
      
      all_models <- stack(predict_bioclim_all, predict_glm_all, predict_maxent_all, predict_rf_all)  
      names(all_models) <- c("BIOCLIM", "GLM", "MAXENT","RANDOM FOREST")
      auc <- sapply(list(eval_bc, eval_glm, eval_ma, eval_rf), function(x) x@auc)
      saveRDS(auc, file = paste(outdir_SDM,species_name,"_all_models_auc.ascii",sep=''),ascii=TRUE)
      w <- (auc-0.5)^2
      ensemble_raster <- weighted.mean( all_models, w)  
      plot(ensemble_raster)
      saveRDS(ensemble_raster, file = paste(outdir_SDM,species_name,"_ensemble_raster_raw.ascii",sep=''),ascii=TRUE)

      print(paste('Ensemble mean done'))
      
      saveRDS(ensemble_raster, file = paste(outdir_SDM,species_name,"_ensemble_eval.ascii",sep=''),ascii=TRUE)
      print(paste('Ensemble evaluation done'))
      
      
      #evaluate ensemble
      
      source("ensemble_evaluate.R")
    
      ens<-ensemble_evaluate(occurrence_lonlat,backgr , ensemble_raster)
      tr_ensemble<-threshold(ens, stat="spec_sens")
      
      print(tr_ensemble)
            
      ensemble_raster_pa <- ensemble_raster > tr_ensemble
      plot(ensemble_raster_pa)
      points(occurrence_in$lon, occurrence_in$lat)
      saveRDS(ensemble_raster_pa, file = paste(outdir_SDM,species_name,"_ensemble_raster_pa.ascii",sep=''),ascii=TRUE)
      
      library(maptools)
      data(wrld_simpl)
      
      tiff(paste(outdir_SDM,species_name,'_ensemble_pa.tiff',sep=''))
      plot(ensemble_raster,legend = FALSE, col = rev(terrain.colors(2)), main=paste('Weighted ensemble mean - ', species_name),xlab="longitue", ylab="latitude")
      plot(wrld_simpl, add=TRUE)
      legend("bottomleft", legend = c("Absence", "Presence"),box.col = "white", fill = rev(terrain.colors(2)))
      box()
      points(occurrence_in$lon,occurrence_in$lat)  
      dev.off()
    
      print(paste(species_name, 'DONE', sep=' '))
    #}
      
      
      ####################################################################################
      ####################################################################################
      #  
      # Land use model
      #  
      ####################################################################################
      ####################################################################################
      
      # load predictor MODIS file
      setwd(predictors_LU_dir)
      f <- list.files(pattern='tif')
      
      cl <- makeCluster(detectCores(), type='PSOCK')
      registerDoParallel(cl)
        predictors_LU_in <- lapply(f, raster)
        predictors_LU_s <- stack(predictors_LU_in)
        predictors_LU <- crop(predictors_LU_s,extent(ensemble_raster))
      stopCluster(cl)
      print(paste('predictors LU loaded'))
      
      # mask land use with occurrence based on climate
      
      cl <- makeCluster(detectCores(), type='PSOCK')
      registerDoParallel(cl)
      ensemble_raster[ensemble_raster ==0 ] <- NA
      predictor_LU_masked <- mask(predictors_LU, ensemble_raster)  
      stopCluster(cl)
    
      set.seed(0)
      backgr_LU <- randomPoints(!is.na(predictor_LU_masked), 500) 
      points(backgr_LU)
      group <- kfold(backgr_LU, 5)
      backg_train <- backgr_LU[group != 1, ]
      backg_test <- backgr_LU[group == 1, ]
      colnames(backg_train) = c('lon', 'lat') 
      colnames(backg_test) = c('lon', 'lat')
      print(paste('background training and test set LU done'))
      
      # extract predictors training LU
      train <- rbind(pres_train, backg_train)
      set_frame_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
      envtrain <- extract(predictor_LU_masked, train)
      envtrain <- data.frame(cbind(pa=set_frame_train, envtrain))
      print(paste('training LU predictors extracted'))
      
      # extract predictors test LU
      test <- rbind(pres_test, backg_test)
      envtestpres <- data.frame( extract(predictor_LU_masked, pres_test) )
      envtestbackg <- data.frame( extract(predictor_LU_masked, backg_test) )
      print(paste('test predictors extracted'))
      
      #-----------------------------------------------------------------------------
      # LU Model ensemble
      #-----------------------------------------------------------------------------
      
      # BIOCLIM
      print(paste(' '))
      print(paste('+++++++++++++++'))
      print(paste('+++ Bioclim +++'))
      print(paste('+++++++++++++++'))
      
        model_bc <- bioclim(predictor_LU_masked, pres_train)  
        saveRDS(model_bc, file = paste(outdir_LU,species_name,"_model_bc_LU.ascii",sep=''),ascii=TRUE)
      
        eval_bc <- evaluate(pres_test, backg_test, model_bc, predictor_LU_masked) 
        saveRDS(eval_bc, file = paste(outdir_LU,species_name,"_eval_bc_LU.ascii",sep=''),ascii=TRUE)
        print(eval_bc)
      
        tr_bc <- threshold(eval_bc, 'spec_sens') 
        saveRDS(tr_bc, file = paste(outdir_LU,species_name,"_tr_bc_LU.ascii",sep=''),ascii=TRUE)
        print(tr_bc)
      
        predict_bioclim <- predict(predictor_LU_masked, model_bc, ext=ext, progress='') 
        saveRDS( predict_bioclim, file = paste(outdir_LU,species_name,"_predict_bioclim_raw_LU.ascii",sep=''),ascii=TRUE)
        print(predict_bioclim)
      
        predict_bioclim <- predict_bioclim > tr_bc
        saveRDS( predict_bioclim, file = paste(outdir_LU,species_name,"_predict_bioclim_pa_LU.ascii",sep=''),ascii=TRUE)

      plot(predict_bioclim)
      points(occurrence_in$lon, occurrence_in$lat)
      print(paste('Bioclim done'))
       
      # GENERALIZED LINEAR MODEL
      print(paste(' '))
      print(paste('+++++++++++'))
      print(paste('+++ GLM +++'))
      print(paste('+++++++++++'))
      
      cl <- makeCluster(detectCores(), type='PSOCK')
      registerDoParallel(cl)
      
        model_glm <- glm(pa ~ MODIS_LCType,
                       family = binomial(link = "logit"), data=envtrain)
      
        saveRDS( model_glm, file = paste(outdir_LU,species_name,"_model_glm_LU.ascii",sep=''),ascii=TRUE)
    
        eval_glm <- evaluate(envtestpres, envtestbackg, model_glm)  
        saveRDS( eval_glm, file = paste(outdir_LU,species_name,"_eval_glm_LU.ascii",sep=''),ascii=TRUE)
        print(eval_glm)
      
        tr_glm <- threshold(eval_glm, 'spec_sens')
        saveRDS( tr_glm, file = paste(outdir_LU,species_name,"_tr_glm_LU.ascii",sep=''),ascii=TRUE)
        print(tr_glm)
      
        predict_glm <- predict(predictor_LU_masked, model_glm, ext=ext)  
        saveRDS( predict_glm, file = paste(outdir_LU,species_name,"_predict_glm_raw_LU.ascii",sep=''),ascii=TRUE)
        print(predict_glm)
      
        predict_glm <- predict_glm >tr_glm
      stopCluster(cl)
      
      plot(predict_glm)
      points(occurrence_in$lon, occurrence_in$lat)
      saveRDS( predict_glm, file = paste(outdir_LU,species_name,"_predict_glm_pa_LU.ascii",sep=''),ascii=TRUE)
      print(paste('GLM done'))
      
      # MAXENT
      print(paste(' '))
      print(paste('++++++++++++++'))
      print(paste('+++ MAXENT +++'))
      print(paste('++++++++++++++'))
      jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

        model_ma <- maxent(predictor_LU_masked, pres_train)#,l1_regularizer=0.7)  
        saveRDS( model_ma, file = paste(outdir_LU,species_name,"_model_ma_LU.ascii",sep=''),ascii=TRUE)
     
        eval_ma <- evaluate(model_ma, p=pres_test, a=backgr_LU, x=predictor_LU_masked) 
        saveRDS(eval_ma, file = paste(outdir_LU,species_name,"_eval_ma_LU.ascii",sep=''),ascii=TRUE)
        print(eval_ma)
      
        tr_ma <- threshold(eval_ma, 'spec_sens')
        saveRDS(tr_ma, file = paste(outdir_LU,species_name,"_tr_ma_LU.ascii",sep=''),ascii=TRUE)
        print(tr_ma)
      
        predict_maxent <- predict(model_ma, predictor_LU_masked,ext=ext) 
        saveRDS(predict_maxent, file = paste(outdir_LU,species_name,"_predict_maxent_raw_LU.ascii",sep=''),ascii=TRUE)
        print(predict_maxent)
      
        predict_maxent <- predict_maxent > tr_ma
        saveRDS(predict_maxent, file = paste(outdir_LU,species_name,"_predict_maxent_pa_LU.ascii",sep=''),ascii=TRUE)

      plot(predict_maxent)
      points(occurrence_in$lon, occurrence_in$lat)
      print(paste('MAXENT done'))
      
      # RANDOM FOREST
      print(paste(' '))
      print(paste('+++++++++++++++++++++'))
      print(paste('+++ Random forest +++'))
      print(paste('+++++++++++++++++++++'))

      
        model_rf1 <- paste("pa ~", paste(names(predictor_LU_masked), collapse=" + "))
      
        model_rf <- randomForest(as.formula(model_rf1), data=envtrain, na.action=na.roughfix)  
        saveRDS(model_rf, file = paste(outdir_LU,species_name,"_model_rf_LU.ascii",sep=''),ascii=TRUE)
      
        eval_rf <- evaluate(envtestpres, envtestbackg, model_rf)  
        saveRDS(eval_rf, file = paste(outdir_LU,species_name,"_eval_rf_LU.ascii",sep=''),ascii=TRUE)
        print(eval_rf)
      
        tr_rf <- threshold(eval_rf, 'spec_sens')
        saveRDS(tr_rf, file = paste(outdir_LU,species_name,"_tr_rf_LU.ascii",sep=''),ascii=TRUE)
        print(tr_rf)
      
        predict_rf <- predict(predictor_LU_masked, model_rf, ext=ext)
        saveRDS(predict_rf, file = paste(outdir_LU,species_name,"_predict_rf_raw_LU.ascii",sep=''),ascii=TRUE)
        print(predict_rf)
      
        predict_rf <- predict_rf > tr_rf
        saveRDS(predict_rf, file = paste(outdir_LU,species_name,"_predict_rf_pa_LU.ascii",sep=''),ascii=TRUE)

      plot(predict_rf)
      points(occurrence_in$lon, occurrence_in$lat)
      print(paste('Random forest done'))
      
      # MODEL ENSEMBLE
      print(paste(' '))
      print(paste('++++++++++++++++++++++'))
      print(paste('+++ Model ensemble +++'))
      print(paste('++++++++++++++++++++++'))

        all_models <- stack(predict_bioclim, predict_glm, predict_maxent, predict_rf)
        names(all_models) <- c("BIOCLIM", "GLM", "MAXENT","RANDOM FOREST")
        auc <- sapply(list(eval_bc, eval_glm, eval_ma, eval_rf), function(x) x@auc)
        saveRDS(auc, file = paste(outdir_LU,species_name,"_all_models_auc_LU.ascii",sep=''),ascii=TRUE)
        w <- (auc-0.5)^2
        ensemble_raster_LU <- weighted.mean( all_models, w)
        saveRDS(ensemble_raster_LU, file = paste(outdir_LU,species_name,"_ensemble_raster_raw_LU.ascii",sep=''),ascii=TRUE)

      print(paste('Ensemble mean done'))
      
      #evaluate ensemble
      p <- extract(ensemble_raster, occurrence_lonlat)
      a <- extract(ensemble_raster, backgr)
      
      p <- stats::na.omit(p)
      a <- stats::na.omit(a)
      np <- length(p)
      na <- length(a)
      
      tr <- p
      tr <- sort(unique(round(tr, 8)))
      tr <- c(tr - 1e-04, tr[length(tr)] + c(0, 1e-04))
      
      N <- na + np
      xc <- new("ModelEvaluation")
      xc@presence <- p
      xc@absence <- a
      R <- sum(rank(c(p, a))[1:np]) - (np * (np + 1) / 2)
      xc@auc <- R / (as.numeric(na) * as.numeric(np))
      cr <- try(cor.test(c(p, a), c(rep(1, length(p)), rep(0, length(a)))),
                silent = TRUE)
      if (class(cr) != "try-error") {
        xc@cor <- cr$estimate
        xc@pcor <- cr$p.value
      }
      res <- matrix(ncol = 4, nrow = length(tr))
      colnames(res) <- c("tp", "fp", "fn", "tn")
      xc@t <- tr
      for (i in 1:length(tr)) {
        res[i, 1] <- length(p[p >= tr[i]])
        res[i, 2] <- length(a[a >= tr[i]])
        res[i, 3] <- length(p[p < tr[i]])
        res[i, 4] <- length(a[a < tr[i]])
      }
      xc@confusion <- res
      a <- res[, 1]
      b <- res[, 2]
      c <- res[, 3]
      d <- res[, 4]
      xc@np <- as.integer(np)
      xc@na <- as.integer(na)
      xc@prevalence <- (a + c) / N
      xc@ODP <- (b + d) / N
      xc@CCR <- (a + d) / N
      xc@TPR <- a / (a + c)
      xc@TNR <- d / (b + d)
      xc@FPR <- b / (b + d)
      xc@FNR <- c / (a + c)
      xc@PPP <- a / (a + b)
      xc@NPP <- d / (c + d)
      xc@MCR <- (b + c) / N
      xc@OR <- (a * d) / (c * b)
      prA <- (a + d) / N
      prY <- (a + b) / N * (a + c) / N
      prN <- (c + d) / N * (b + d) / N
      prE <- prY + prN
      xc@kappa <- (prA - prE) / (1 - prE)
      
      xc
      tr_ensemble <- threshold(xc,'spec_sens')
      
      saveRDS(xc, file = paste(outdir_LU,species_name,"_ensemble_eval_LU.ascii",sep=''),ascii=TRUE)
      print(paste('Ensemble evaluation LU done'))
      
      ensemble_raster_LU <- ensemble_raster_LU > tr_ensemble
      saveRDS(ensemble_raster_LU, file = paste(outdir_LU,species_name,"_ensemble_raster_pa_LU.ascii",sep=''),ascii=TRUE)
      
      tiff(paste(outdir_LU,species_name,'_ensemble_pa_LU.tiff',sep=''))
      plot(ensemble_raster_LU,legend = FALSE, col = rev(terrain.colors(2)), main=paste('Weighted ensemble mean SDM-LU - ', species_name),xlab="longitue", ylab="latitude")
      plot(wrld_simpl, add=TRUE)
      legend("bottomleft", legend = c("Absence", "Presence"),box.col = "white", fill = rev(terrain.colors(2)))
      box()
      points(occurrence_in$lon,occurrence_in$lat)  
      dev.off()
      
      print(paste(species_name, 'DONE', sep=' '))
    #  }
	