# gPipeline and wrappers : main pipeline for Microvelia longipes image analysis
#---------------------------------------------------------------------------------------------------

#' Waterstrider image analysis pipeline
#'
#' Finds individual's size, sex, presence of wing, tibia and femur length
#'
#' @param img_path image path
#' @param write_output logical, to save metrics as csv and print image analysis diagnostic plots
#' @param return_df logical, to return metrics in R
#' @param single_write logical, writing method for gWritePipeline()
#' @param thresh numerical value in range 0:1, global gray-value threshold
#' for the initial binarization of the input image
#' @param winged_names vector, two names for not-winged / winged individuals
#'
#' @export
gPipeline <- function(img_path, write_output=T, return_df=T,
                      single_write=T, thresh=0.8, winged_names=c(0,1)){
  message("1 - Loading and segmenting image")
  base_img <- imager::load.image(img_path)
  #Loading and Binarizing image
  img_bin <- binaryLoad(img_path, thresh) #load and binarize image
  #Segmenting image
  body_lab_points <- gFastSeg(img_bin) #segmentation
  #Body base metrics
  body_length_pix <- unlist(lapply(body_lab_points, bodyLength)) #length in pixels of bodies
  body_centroids <- lapply(body_lab_points, function(i) { #body centroids
    apply(i, 2, mean)
  })
  message("| \n2 - Computing scale")
  #Scale
  scale <- redRulerScale(base_img, msg=F) #get scale
  message("| \n3 - Computing body")
  body_length <- body_length_pix/scale[1] #conversion in mm
  error_margin <- body_length*scale[2]/scale[1] #scale error margin (ignoring pixel error)
  #clean bodies
  clean_body_points <- cleanBodyShape(body_lab_points)
  #Remove bodies
  im_nobodies <- gNobody(base_img = base_img, body_lab_points = clean_body_points, viz=F)
  #Size based ind crop
  l_crop <- gCrop(im_nobodies, body_centroids, body_length_pix, viz=F)
  #Fit GMM to image values
  message("| \n4 - Individual thresholding")
  #***Formats
  body_l_crops <- recropBodies(body_coords = clean_body_points, #Body to cropfull reference
                               base_img = base_img,
                               crop_coords = l_crop$crop_coords)
  body <- body_l_crops %>% imgAsCoords #body as coordinates, same reference as dilcont
  cen <- lapply(body, function(x){ apply(x,2,mean) }) #new centroids
  #GMM threshold
  l_cropbin <- gGMMThresh(l_crop$img, msg=F)
  #Clean data
  full <- gCleanBin(l_cropbin, cen, as_coords=F)
  
  message("| \n5 - Orienting individuals")
  #1 - Dilate body contour on same coords as leg crops
  dilbody <- dilBodies(body_img = body_l_crops,
                       body_length = body_length_pix) #contour of the dilation as the intersection line
  #2 - Dilated body contour as coords
  dilcont <- simpleCont(dilbody)
  #3 - Intersections (as boolean of rows in dilcont)
  inter <- overlapPoints(dilcont, lapply(full,imgAsCoords))
  #*** - Formats
  inter_coords <- lapply(inter, function(x) x$coords)
  inter_idx <- lapply(inter, function(x) x$index)
  #4 - Orientation (PCA+intersection|density)
  ang <- gOrientation(body, inter_coords)
  
  message("| \n6 - Leg segmentation")
  #5 - Hind legs insertions
  inser <- gLegInsertion(ori_angle=ang, dil_contour=dilcont, inter_index=inter_idx)
  #6 - Hind legs segmentation
  legs <- gLegSeg(gerris = full,
                  dilated_body = dilbody,
                  intersection_coords = inter_coords,
                  insertions = inser)
  message("| \n7 - Leg landmarking & measurement")
  #7 - Get Points of interest from leg
  leg_lm0 <- gLegLandmarks(legs, inser, viz=F, msg=F)
  #8 - Correct insertion landmark by connection leg to body
  leg_lm <- gConnectLeg(body, leg_lm0)
  #9 - Distances
  leg_size <- gMeasureLeg(leg_lm, scale) #leg segment sizes in microns
  #Sex and wing prediction using body contour
  message("| \n8 - Sex and wing prediction")
  #elliptic fourier harmonics reduced with PCA
  hPCA <- scoresEFA(body_l_crops, ang, nb_h=7, viz=F)
  #loading LDA models from package
  LDA_sex <- loadLDA("LDAsex")
  LDA_wingF <- loadLDA("LDAwingF")
  LDA_wingM <- loadLDA("LDAwingM")
  #predict sex using LDA model
  sex_prediction <- gPredictLDA(LDA_sex, hPCA) #LDA trained on 2PC of 7EF
  #Predict presence of wings using LDA model
  wing_prediction <- gPredictWing(PCA_scores = hPCA, sex_class = sex_prediction$class,
                                  LDA_f = LDA_wingF, LDA_m = LDA_wingM, wing_names = winged_names)
  
  #OUTPUTS
  clean_base_path <- sub("\\.(jpe?g|tif|png|bmp|gif|jpg)$", "", basename(img_path), ignore.case = TRUE)
  leg_res <- gSplitLegMeasures(leg_size)
  df_out <- data.frame(image =  rep(clean_base_path, length(body_lab_points)),
                       i = seq_along(body_lab_points),
                       L_body = body_length,
                       error_margin = error_margin,
                       left_tibia = leg_res$right_tibia,#inverting left/right due to error
                       right_tibia = leg_res$left_tibia,
                       left_femur = leg_res$right_femur,
                       right_femur = leg_res$left_femur,
                       sex = sex_prediction$class,
                       F_proba = sex_prediction$posterior[,'F'],
                       winged = wing_prediction$class,
                       w_proba = wing_prediction$posterior[,winged_names[2]])
  
  if(write_output){
    # Detection plot of the image
    x_cen <- sapply(body_centroids, function(x) x[1])
    y_cen <- sapply(body_centroids, function(x) x[2])
    detection_plot <- function(){ #save arguments to be called later
      gDetectionPlot(base_img=base_img, scale=scale[1], x=x_cen, y=y_cen)
    }
    # Individuals metrics plot
    i_plots <- lapply(seq_along(body_l_crops), function(i) {
      function(){ #save arguments to be called later
        gGerrisPlot(i, full, body, cen, dilcont, ang, legs,
                    leg_lm, leg_size, inser, clean_base_path)
      }
    })
    # Write plots and dataframe
    gWritePipeline(img_path, i_plots, detection_plot, df_out, dim_img=dim(base_img), single_write)
  }

  if(return_df){
    return(df_out)
  }
}


#' Write waterstriders pipeline outputs in local directory
#' 
#' Plot and dataframe Writing function for gPipeline, an alternative printing method is
#' implemented for gMultiPipeline by passing single_write=False
#' 
#' @param img_path valid path to processed image, used to define writing directory 
#' @param i_plots list of plots of individuals as created by gPipeline with write_output=T 
#' @param detection_plot whole image detection plot as created by gPipeline with write_output=T 
#' @param df_out corresponding dataframe created by gPipeline
#' @param dim_img dimensions of base image
#' @param single_write writing method
#'  True : a directory is created for the current image, containing plots and dataframe
#'  False : a single 'out' directory with all plots is created for the current directory
#'  ,gMultiPipeline is supposed to take in charge the creation of the dataframe
gWritePipeline <- function(img_path, i_plots, detection_plot, df_out, dim_img, single_write=T){
  clean_base_path <- sub("\\.(jpe?g|tif|png|bmp|gif|jpg)$", "",
                         basename(img_path), ignore.case = TRUE)
  if(single_write){
    local_dir <- paste0(dirname(img_path),"/out_", clean_base_path) #solo directory
    iplot_dir <- paste0(local_dir,"/individuals_",clean_base_path)
    create_dir(local_dir, msg=F)
    create_dir(iplot_dir, msg=F) #individuals dir
    dplot_path <- paste0(local_dir,"/detection_",clean_base_path,".png")
    df_path <- paste0(local_dir,"/data_",clean_base_path,".csv")
    write.table(df_out, file = df_path, row.names = FALSE, sep=";") 
  } else { #multi image processing, usually with gMultiPipeline
    out_dir <- paste0(dirname(img_path),"/out_",basename(dirname(img_path))) #output directory name
    create_dir(out_dir, msg=F) #create if not already done
    dplot_dir <- paste0(out_dir,"/detection")
    create_dir(dplot_dir, msg=F)
    iplot_dir <- paste0(out_dir,"/individuals_",clean_base_path)
    create_dir(iplot_dir, msg=F) #individuals dir
    dplot_path <- paste0(dplot_dir,"/detection_",clean_base_path,".png")
    #global dataframe created in gMultiPipeline
  }
  #individual's plot
  lapply(seq_along(i_plots), function(i){
    iplot_path <- paste0(iplot_dir,"/",clean_base_path,"_",i,".png")
    png(iplot_path, width = 480, height = 480)
    i_plots[[i]]() #call saved plot function
    dev.off()
  })
  #detection plot
  png(dplot_path, width = dim_img[1], height = dim_img[2])  # specify filename and dimensions
  detection_plot() #call saved plot function
  dev.off()
}

#' Pipeline iteration for multiple images 
#' 
#' gPipeline iteration for multiple images including specific writing option,
#' and appropriate console messages
#'
#' @param img_path_list a vector or list or valid image paths
#' @param return_df logical, to return the concatenated gPipeline dataframe
#' @export
gMultiPipeline <- function(img_path_list, return_df){
  img_path_length <- length(img_path_list)
  df_out_list <- list()
  message(paste0("Running pipeline for ", img_path_length, " images \n_____________________________"))
  #Run pipeline for each image
  for(img_path_i in 1:img_path_length){
    clean_base_path <- sub("\\.(jpe?g|tif|png|bmp|gif|jpg)$", "", basename(img_path_list[img_path_i]), ignore.case = TRUE)
    message(paste0("\n  ",img_path_i,"/",img_path_length, " - Processing image \"",clean_base_path,"\"\n  |"))
    df_out_list[[clean_base_path]] <- gPipeline(img_path_list[img_path_i], #save output dataframe in list
              write_output=T,
              return_df=T, #return dataframe of each image analysis
              single_write=F) #writing format for whole directory
  }
  #create and save dataframe
  multi_df_out <- do.call(rbind, df_out_list)
  dir_name <- dirname(img_path_list)[1] #single element for all the list
  main_dir <- basename(dir_name) #name to use for output directories = original directory
  df_out_path <- paste0(dir_name, "/out_", main_dir, "/data_", main_dir, ".csv") #output directory name
  write.table(multi_df_out, file = df_out_path, row.names = FALSE, sep=";") 
  if(return_df){
    return(multi_df_out)
  }
}

#' Main method to run gPipeline for one or multiple images in the same directory
#'
#' Checks if an input can be valid for either gPipline or gMultiPipeline,
#' and use the appropriate function, or return a specific error message
#'
#' @param img_path image path or vector of image paths to analyze 
#' @param return_df logical to return the output dataframe
#' @export
gRunPipeline <- function(img_path, return_df=T){
  #1 check if path is an image
  img_pattern <- "\\.(jpe?g|tif|png|bmp|gif|jpg)$"
  if( grepl(img_pattern, img_path, ignore.case = TRUE) ){ #if provided path is a single img path
    if( !file.exists(img_path) ){
      message(paste0("File ",img_path," appears to be an image but was not found in ", getwd()))
    } else {
      gPipeline(img_path, return_df=return_df)
    }
  } else { #path is not an image, assuming a folder containing images
    if( !file.exists(img_path) ){ #invalid path
      message(paste0("File ",img_path," was not found in ", getwd()))
    } else {
      if( !file.info(img_path)$isdir ) { #path valid but not a directory
        message(paste0(img_path, " is not an image or a directory"))
      } else { #path is a valid dir
        lf <- list.files(img_path, pattern=img_pattern, ignore.case=T, full.names=T) #img paths of the dir
        llf=length(lf)
        if(llf==0){
          message(paste0(img_path, " is a valid directory but does not contain any image"))
        } else if(llf==1){
          message("Running pipeline for the single image found in directory")
          return(gPipeline(lf, return_df=return_df))
        } else { #if multiple images are found in dir
          gMultiPipeline(lf, return_df=return_df) 
        }
      }
    }
  }
}