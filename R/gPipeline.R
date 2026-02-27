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
#' @param auto_scale a boolean. If 'TRUE' : red ruler based scale
#' @param predict_sex_wing a boolean. To predict sex and presence of wings
#' @param k.bin_thresh numerical in \code{[0, 1]}. Global gray-value threshold
#' for the initial binarization of the input image
#' @param k.clean_kernel numerical. Kernel size for morphological opening and closing 
#' of body points to get a cleaner shape.
#' @param k.crop_size_factor numerical. Individual's crop size relative to body length
#' @param k.gmm_slope_value numerical. y-value of individual's plot model
#' (derivative of GMM used to fit gray values) to use to define threshold between individual's
#' gray values and those of the background.
#' @param k.body_dilation_ratio numerical in \code{[0, 1]}. Dilation kernel size as body length ratio for creation
#' of dilated body.
#' @param k.leglm_search_w numerical. Distance (in number of points in contour) to use to compute leg contour angles
#' @param k.leglm_splines_df numerical. Number of splines to fit leg contour angle variation
#' @param k.body_size_px_range numerical vector of length 2. range in pixels of expected body size. Check imager::label(image_binary) with boxplots beforehand
#' @param f.red_thresh numerical in \code{[0, 1]}. Difference in intensity between red component
#' and other components to detect the red ruler. Only relevant if auto_scale == TRUE
#' @param f.clean_small_spots numerical. Minimal size (in number of pixels) allowed for a patch
#' of with pixels (label) found after individual thresholding.
#' @param f.leg_lim_ratio numerical vector of length 2. Values between 0 and 1 to define range where
#' legs insertions are expected to be found on the dilated body contour as ratio of body elongation
#' starting from the back (thus excluding head appendages and hind body part)
#' @param f.leg_inser_thresh numerical in \code{[0, 1]}. As ratio of leg half contour length.
#' Threshold to remove angular points too close to the leg insertion point when starting to look for
#' knee points in leg landmarks.
#' @param f.leg_inser_knee_thresh numerical in \code{[0, 1]}. Threshold for knee-insertion points dissimilarity
#' Tolerated ratio between the distances from the insertion point to each detected knee
#' A value close to 0 will ensure knee-insertion distances are similar
#' @param f.leg_ankle_knee_thresh numerical in \code{[0, 1]}. Threshold for knee-ankle points dissimilarity
#' As ratio of leg half contour length. When detecting ankles point using knee points, threshold
#' to detect if the difference between the two knee-ankle distances is realistically acceptable 
#' @param f.leg_n_inflexions numerical vector of length 2. Tolerated range of inflexion points
#' on leg contour. Usually 1 or 2 for extreme ends of the limb, and 2 per expected joint in the
#' middle of it.
#' @param f.segment_lengths_range numerical vector of length 2 n \code{[0, 1]}.
#' Range of tolerated leg segments sizes relative to leg half contour length.
#' @param return_everything boolean. to return key process step data in R.
#' 
#' @export
gPipeline <- function(# INPUT
                      img_path,
                      # OUTPUT
                      write_output=TRUE,
                      return_df=TRUE,
                      single_write=TRUE,
                      auto_scale=TRUE,
                      predict_sex_wing=TRUE,
                      return_everything=FALSE,
                      # In-pipeline parameters
                      # a. Key variables
                      k.bin_thresh = 0.8, #body
                      k.body_size_px_range = c(300,1500),
                      k.clean_kernel = 3,
                      k.crop_size_factor = 3.5,
                      k.gmm_slope_value = -25,
                      k.body_dilation_ratio = 0.3,
                      k.leglm_search_w = 6,
                      k.leglm_splines_df = 30,
                      # b. Filters
                      f.red_thresh = 0.05,
                      f.clean_small_spots = 25,
                      f.leg_lim_ratio = c(.25,.6),
                      f.leg_inser_thresh = .1,
                      f.leg_inser_knee_thresh = .2,
                      f.leg_ankle_knee_thresh = .12,
                      f.leg_n_inflexions = c(5, 7),
                      f.segment_lengths_range = c(.15, .6)
                      ){
  message("1 - Loading and segmenting image")
  base_img <- imager::load.image(img_path)
  #Loading and Binarizing image
  img_bin <- binaryLoad(img_path,
                        threshold = k.bin_thresh) #load and binarize image
  #Segmenting image
  body_lab_points <- gFastSeg(img_bin, px_range=k.body_size_px_range, viz=FALSE) #segmentation
  body_centroids <- lapply(body_lab_points, function(i) apply(i, 2, mean)) #body centroids
    
  message("| \n2 - Computing scale")
  #Scale
  scale <- pipelineScale(base_img, auto_scale,
                         red_thresh = f.red_thresh) #get scale
  message("| \n3 - Computing body")
  #Clean bodies
  clean_body_points <- cleanBodyShape(body_lab_points,
                                      kernel_size = k.clean_kernel)
  #Body length
  bl_output <- bodyLength(clean_body_points, return_ext=TRUE)
  body_length_pix <- bl_output$len
  body_length <- body_length_pix/scale[1] #conversion in mm
  if(auto_scale){
    error_margin <- body_length*scale[2]/scale[1] #scale error margin (ignoring pixel error)
  } else error_margin <- NA
  #Remove bodies
  im_nobodies <- gNobody(base_img = base_img, body_lab_points = clean_body_points, viz=FALSE)
  #Size based ind crop
  l_crop <- gCrop(im_nobodies, body_centroids, body_length_pix, viz = FALSE,
                  factor = k.crop_size_factor)
  #Fit GMM to image values
  message("| \n4 - Individual thresholding")
  #***Formats
  body_l_crops <- recropBodies(body_coords = clean_body_points, #Body to cropfull reference
                               base_img = base_img,
                               crop_coords = l_crop$crop_coords)
  body <- body_l_crops %>% imgAsCoords #body as coordinates, same reference as dilcont
  cen <- lapply(body, function(x){ apply(x,2,mean) }) #new centroids
  #GMM threshold
  l_cropbin <- gGMMThresh(l_crop$img, msg=FALSE,
                          y_root = k.gmm_slope_value)
  #Clean data
  full <- gCleanBin(l_cropbin, cen, as_coords=FALSE,
                    px_filter = f.clean_small_spots)
  
  message("| \n5 - Orienting individuals")
  #1 - Dilate body contour on same coords as leg crops
  dilbody <- dilBodies(body_img = body_l_crops,
                       body_length = body_length_pix, #contour of the dilation as the intersection line
                       dilation_ratio = k.body_dilation_ratio)
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
  inser <- gLegInsertion(ori_angle=ang, dil_contour=dilcont, inter_index=inter_idx,
                         leg_lim_ratio = f.leg_lim_ratio)
  #6 - Hind legs segmentation
  legs <- gLegSeg(gerris = full,
                  dilated_body = dilbody,
                  intersection_coords = inter_coords,
                  insertions = inser)
  message("| \n7 - Leg landmarking & measurement")
  #7 - Get Points of interest from leg
  leg_lm0 <- gLegLandmarksLoop(leg_coords = legs, insertion = inser, viz=FALSE, msg=FALSE,
                           search_w = k.leglm_search_w,
                           n_splines = k.leglm_splines_df,
                           inser_thresh = f.leg_inser_thresh,
                           knee_diff_thresh = f.leg_inser_knee_thresh,
                           tresh_ankle = f.leg_ankle_knee_thresh,
                           inflexion_pts_range = f.leg_n_inflexions,
                           segment_length_range = f.segment_lengths_range
                           )
  #8 - Correct insertion landmark by connection leg to body
  leg_lm <- gConnectLeg(body, leg_lm0)
  #9 - Distances
  leg_size <- gMeasureLeg(leg_lm, scale) #leg segment sizes in microns
  #Sex and wing prediction using body contour

  if(predict_sex_wing){
    message("| \n8 - Sex and wing prediction")
    #prediction pipelines: body img > contour > EFA > PCA > LDA, using models from package
    sex_prediction <- gBodyPredict(img_data = body_l_crops, #NA handling built-in
                                   angle = ang,
                                   what = "sex")
    wing_prediction <- gBodyPredict(img_data = body_l_crops,
                                    angle = ang,
                                    what = "wing")
    }
  
  #OUTPUTS
  clean_base_path <- sub("\\.(jpe?g|tif|png|bmp|gif|jpg)$", "", basename(img_path), ignore.case = TRUE)
  leg_res <- gSplitLegMeasures(leg_size)
  df_out <- data.frame(image =  rep(clean_base_path, length(body_lab_points)),
                       i = seq_along(body_lab_points),
                       L_body = body_length %>% round,
                       error_margin = error_margin,
                       left_tibia = leg_res$right_tibia %>% round,#inverting left/right due to error
                       right_tibia = leg_res$left_tibia %>% round,
                       left_femur = leg_res$right_femur %>% round,
                       right_femur = leg_res$left_femur %>% round)
  legText(df_out)
  if(predict_sex_wing){ #add columns for sex/wing predictions
    df_sw <- data.frame(
      sex = c("F","M")[sex_prediction$class],
      F_prob = sex_prediction$posterior[,"F"],
      M_prob = sex_prediction$posterior[,"M"],
      wing = (wing_prediction$class %>% as.numeric)-1, #scale 2-1 (model levels) to 1-0
      winged_prob = wing_prediction$posterior[,"2"],
      wingless_prob = wing_prediction$posterior[,"1"]
    )
   df_out <- cbind(df_out, df_sw)
 }
  
  if(write_output){
    # Detection plot of the image
    x_cen <- sapply(body_centroids, function(x) x[1])
    y_cen <- sapply(body_centroids, function(x) x[2])
    detection_plot <- function(){ #save arguments to be called later
      gDetectionPlot(base_img=base_img, scale=scale, x=x_cen, y=y_cen, auto_scale=auto_scale)
    }
    # Individuals metrics plot
    body_L_pts <- bl_output$body_L_pts
    i_plots <- lapply(seq_along(body_l_crops), function(i) {
      function(){ #save arguments to be called later
        gGerrisPlot(i, full, body, cen, dilcont, ang, legs,
                    leg_lm, leg_size, inser, clean_base_path, body_L_pts, body_length)
      }
    })
    # Write plots and dataframe
    gWritePipeline(img_path, i_plots, detection_plot, df_out, dim_img=dim(base_img), single_write)
  }
  
  output <- list()
  if(return_df){
    output[["df"]] <- df_out
  }
  if(return_everything){ #return mid-pipeline data (for model training or diagnostics for example)
    output[["process_data"]] <- list(body_img = body_l_crops, angle = ang, legs = legs,
                                     dilcont = dilcont, inter_idx = inter_idx,
                                     body_centroids = body_centroids, scale = scale, full = full,
                                     inser = inser, crop_base = l_crop)
  }
  return(output)
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
gWritePipeline <- function(img_path, i_plots, detection_plot, df_out, dim_img, single_write=TRUE){
  clean_base_path <- sub("\\.(jpe?g|tif|png|bmp|gif|jpg)$", "",
                         basename(img_path), ignore.case = TRUE)
  if(single_write){
    local_dir <- paste0(dirname(img_path),"/out_", clean_base_path) #solo directory
    iplot_dir <- paste0(local_dir,"/individuals_",clean_base_path)
    create_dir(local_dir, msg=FALSE)
    create_dir(iplot_dir, msg=FALSE) #individuals dir
    dplot_path <- paste0(local_dir,"/detection_",clean_base_path,".png")
    df_path <- paste0(local_dir,"/data_",clean_base_path,".csv")
    write.table(df_out, file = df_path, row.names = FALSE, sep=";") 
  } else { #multi image processing, usually with gMultiPipeline
    out_dir <- paste0(dirname(img_path),"/out_",basename(dirname(img_path))) #output directory name
    create_dir(out_dir, msg=FALSE) #create if not already done
    dplot_dir <- paste0(out_dir,"/detection")
    create_dir(dplot_dir, msg=FALSE)
    iplot_dir <- paste0(out_dir,"/individuals_",clean_base_path)
    create_dir(iplot_dir, msg=FALSE) #individuals dir
    dplot_path <- paste0(dplot_dir,"/detection_",clean_base_path,".png")
    #global dataframe will be created in gMultiPipeline
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
#' @param ... Additional arguments passed to `gPipeline()`
#' 
#' @return A list of results returned by `gPipeline()`.
#' @seealso [gPipeline()]
#' @export
gMultiPipeline <- function(img_path_list, return_df, ...){
  img_path_length <- length(img_path_list)
  df_out_list <- list()
  message(paste0("Running pipeline for ", img_path_length, " images \n_____________________________"))
  #Run pipeline for each image
  for(img_path_i in 1:img_path_length){
    clean_base_path <- sub("\\.(jpe?g|tif|png|bmp|gif|jpg)$", "",
                           basename(img_path_list[img_path_i]), ignore.case = TRUE)
    message(paste0("\n  ",img_path_i,"/",img_path_length,
                   " - Processing image \"",clean_base_path,"\"\n  |"))
    tryCatch({
      pipe_result <- gPipeline(img_path_list[img_path_i], #save output dataframe in list
                               write_output=TRUE,
                               return_df=TRUE, #return dataframe of each image analysis
                               single_write=FALSE, #writing format for whole directory
                               ...) #gPipeline arguments
      if(return_df) df_out_list[[clean_base_path]] <- pipe_result[["df"]]
    }, error = function(e) { #gPipeline full crash handling, probably needs debugging
      warning(paste("gPipeline Failure on",clean_base_path,
                    " no result will be returned for this image.\n",e$message))
    })
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
#' @param ... Additional arguments passed to `gPipeline()`
#' 
#' @return A list of results returned by `gPipeline()`.
#' @seealso [gPipeline()]
#' @export
gRunPipeline <- function(img_path, return_df=TRUE, ...){
  t0 <- Sys.time()
  #1 check if path is an image
  img_pattern <- "\\.(jpe?g|tif|png|bmp|gif|jpg)$"
  if( grepl(img_pattern, img_path, ignore.case = TRUE) ){ #if provided path is a single img path
    if( !file.exists(img_path) ){
      message(paste0("File ",img_path," appears to be an image but was not found in ", getwd()))
    } else {
      gPipeline(img_path, return_df=return_df, ...)
    }
  } else { #path is not an image, assuming a folder containing images
    if( !file.exists(img_path) ){ #invalid path
      message(paste0("File ",img_path," was not found in ", getwd()))
    } else {
      if( !file.info(img_path)$isdir ) { #path valid but not a directory
        message(paste0(img_path, " is not an image or a directory"))
      } else { #path is a valid dir
        lf <- list.files(img_path, pattern=img_pattern,
                         ignore.case=TRUE, full.names=TRUE) #img paths of the dir
        llf <- length(lf)
        if(llf==0){
          message(paste0(img_path, " is a valid directory but does not contain any image"))
        } else if(llf==1){
          message("Running pipeline for the single image found in directory")
          return(gPipeline(lf, return_df=return_df, ...))
        } else { #if multiple images are found in dir
          gMultiPipeline(lf, return_df=return_df, ...) 
        }
      }
    }
  }
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "min"))
  message(sprintf("**********\nJob done in %.2f min", elapsed))
}

#' Keep only one tibia and one femur column
#'
#' Transform result dataframe to only have one tibia and one femur columns
#' If both legs are present, return mean values. Removes NA lines.
#'
#' @param df_result result dataframe as outout of gPipeline or gRunPipeline
#' @export
gOneLeg <- function(df_result){
  
  names_v <- c("left_femur","left_tibia","right_femur","right_tibia")
  lleg <- list(
    L=df_result[,names_v[1:2]],
    R=df_result[,names_v[3:4]]
  )
  r_na <- !is.na(lleg$R$right_tibia)
  l_na <- !is.na(lleg$L$left_tibia)
  
  one_leg <- xor(l_na, r_na)
  left_leg <- one_leg & l_na #get left only
  right_leg <- one_leg & r_na #get right only
  both_leg <- l_na & r_na #get both
  df_leg_choice <- data.frame(both_leg,left_leg,right_leg)
  
  new_leg_l <- lapply(seq_len(nrow(df_result)), function(i){
    if(df_leg_choice[i,"both_leg"]){
      matrix_data <- cbind(as.numeric(lleg$R[i,]), as.numeric(lleg$L[i,]))
      as.numeric(rowMeans(matrix_data))
    } else if(df_leg_choice[i,"left_leg"]) {
      as.numeric(lleg$L[i,])
    } else if(df_leg_choice[i,"right_leg"]) {
      as.numeric(lleg$R[i,])
    } else { #none
      c(NA,NA)
    }
  })
  new_leg <- do.call(rbind, new_leg_l)
  df_result2 <- df_result[,!(names(df_result)) %in% names_v]
  df_result2[,c("femur","tibia")] <- new_leg
  df_result3 <- df_result2[!is.na(df_result2$femur),]
  
  return(df_result3)
}

#' Quick leg segmentation summary
#'
#' Femur segmentation summary for quick diagnostic in the console
#'
#' @param df A dataframe with df_out format from gPipeline
legText <- function(df){
  lf <- df$left_tibia
  rf <- df$right_tibia
  n <- nrow(df)
  nof <- sum(is.na(lf)&is.na(rf))
  ratiof <- round((sum(!is.na(lf)) + sum(!is.na(rf)))/(n*.02),1)
  message(
    paste0("| \n|_ Tibia length measured on at least one leg in ", n-nof, "/",
           n, " individuals.\n|_ Measured tibias ratio: ", ratiof, "%")
  )
}
