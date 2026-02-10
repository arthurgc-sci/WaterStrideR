# Gerroidea specific functions, used in gPipeline
#---------------------------------------------------------------------------------------------------

#' Fast segmentation with size filter
#'
#' Given a binary image as pixset, returns list of coordinates of every label with pixels in the min-max range
#' @param img_bin a pixset
#' @param px_range numerical vector of length 2. range in pixels of expected body size. Check imager::label(image_binary) with boxplots beforehand
#' @param viz boolean for visualization
#' @export
gFastSeg <- function(img_bin, px_range=c(300,1500), viz=TRUE){
  min <- px_range[1]
  max <- px_range[2]
  # 1 - global size threshold (keep bodies)
  lab <- imager::label(img_bin) #connected components detection on binary image
  lab_tab <- table(lab[lab>0]) #number of pixels of each white label
  valid_body_lab <- which(lab_tab>min & lab_tab<max)
  valid_lab_tab <- lab_tab[valid_body_lab]
  # 2 - 60% median threshold (remove eventual connected bodies)
  med_body_L <- lab_tab[valid_body_lab] %>% median #median body size (pixel surface)
  valid_lab_tab <- lab_tab[valid_body_lab] #label size of size filtered labels
  single_body_lab <- valid_body_lab[valid_lab_tab < med_body_L*1.6] #remove labels 60% bigger than the median
  body_lab_points <- lapply(single_body_lab, function(x) { #body points as xy coordinates
    which(lab == x, arr.ind = TRUE)[, 1:2]
  })
  names(body_lab_points) <- seq_along(body_lab_points) #names for indexation
  if(viz){
    plot(img_bin)
    L_seg <- length(body_lab_points)
    cols <- grDevices::rainbow(L_seg)
    for(i in 1:L_seg) points(body_lab_points[[i]],col=cols[i],cex=.01)
  }
  return(body_lab_points)
}

#' Plot gerris position and id with 1cm scale on original image
#'
#' @param base_img c.img or pixset of base image
#' @param scale mean pixel amount for 1 micron
#' @param x,y ordered position of individuals
#' @param auto_scale boolean, TRUE if scale was produced by redRulerScale, FALSE if it was produced by getScale 
#' @export
gDetectionPlot <- function(base_img, scale, x, y, auto_scale){
  if(auto_scale){
    scale_factor <- scale[1] * 10000
    scale_label <- "1 cm"
  } else {
    scale_factor <- scale[1] * scale[2]
    scale_label <- "unit"
  }
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot(base_img)
  x_start <- 100
  y_start <- ncol(base_img) - 100
  graphics::segments(x_start, y_start, x_start + scale_factor, 
                     y_start, col = "brown4", lwd = 4)
  text(x = (2 * x_start + scale_factor)/2, y = y_start - 
         50, labels = scale_label, col = "brown4", cex = 2)
  points(x, y, pch = 4, col = "white", cex = 2.5)
  text(x, y + 40, labels = seq_along(x), col = "cyan", cex = 2.5)
}

#' Cropping points list
#'
#' Accessory function for gCrop,
#'
#' @param list_points A list of numeric vectors, each representing a point (typically of length 2: x and y).
#' @param list_size A numeric vector or list of same length as \code{list_points}, representing the diameter or window size around each point.
#' @param scaling_factor A numeric scalar used to scale the size of each bounding box (default is 1).
#'
#' @return A list of numeric vectors of length 4 (xmin, xmax, ymin, ymax) representing the cropped bounding boxes.
#'
cropPoints <- function(list_points, list_size, scaling_factor=1){
  out <- mapply(function(pt, dia){
    crop_points <- sapply(pt, function(x) rep(x,2)) + c(-0.5,0.5,-0.5,0.5)*dia
  }, list_points, list_size*scaling_factor, SIMPLIFY = FALSE)
  return(out)
}

#' Crop individual gerris
#'
#' Square crop around individual gerris centroids from image base_img with a diameter of factor*size
#' Base image must be imported with imager package, centroids and sizes must be lists of same length
#' Returns a list of inverted grayscale cropped images
#'
#' @param base_img c.img. Base image suitable for waterstrider pipeline
#' @param centroids list of xy vectors. Body centroids coordinates of individuals
#' @param sizes list or vector. Body length of individuals
#' @param viz logical. visualization option
#' @param factor numerical. Numerical value for crop size relative to individual's body length
#'
#' @export
gCrop <- function(base_img, centroids, sizes, viz=FALSE, factor=3.5){
  #Error handling
  if(!imager::is.cimg(base_img)) stop("Error : base_img must be of type c.img")
  if(!is.list(centroids)) stop("Error : centroids must be a list")
  if(!is.list(sizes) && !is.vector(sizes)) stop("Error : sizes must be a vector or a list")
  if(length(sizes)!=length(centroids)) stop("Error : centroids and sizes lengths must be equal")
  #Crop
  crop_points0 <- cropPoints(centroids, sizes, factor)
  crop_points <- lapply(crop_points0, round) #'manual' round to avoid useless computation
  pb <- progress::progress_bar$new(total=length(centroids))
  res <- lapply(crop_points, function(crop_pt_df){
    xrange <- crop_pt_df[,1]
    yrange <- crop_pt_df[,2]
    x <- NA
    y <- NA #avoid warnings in devtools::check() and goodpractice::gp()
    crop <- imager::imsub(base_img, x %inr% xrange, y %inr% yrange)
    crop[which(crop==1)] <- max(crop[which(crop!=1)])
    out <- crop %>% imager::grayscale() %>% invert_grayscale #%>% correct_illumination #TODO handle ruler > R(img) negatively impact global segmentation performance
    pb$tick()
    return(out)
  })
  #Visuals
  if(viz){
    par(mfrow=c(1,1))
    plot(base_img) #plot base image
    cen_df <- do.call(rbind, centroids)
    text(cen_df, asp=1, pch=4, cex=0.5, col=2) #plot centroids
    lapply(crop_points, function(x){
      square_pts <- matrix(x[c(1,1,2,2,1,
                               3,4,4,3,3)], ncol=2) #format for square
      points(square_pts, type="l", col=3) #plot square
    }) %>% invisible #avoid NULL outputs of 'lapply' in the console caused by the use of lapply for graphics only
  }
  return(list(img=res, crop_coords=crop_points))
}

#'Remove gerris bodies
#'
#' Changes the values of given set of pixel coordinate body_lab_points to 0
#' in an image base_img
#'
#' @param base_img c.img. Base image suitable for waterstrider pipeline
#' @param body_lab_points a list of dataframes, a dataframe or a matrix
#' of coordinates of individual's body points 
#' @param viz logical.visualization option
#'
#' @export
gNobody <- function(base_img, body_lab_points, viz=FALSE){
  #Error handling
  if(!is.matrix(body_lab_points) & !is.list(body_lab_points) & !is.matrix(body_lab_points)) stop("Error : body_lab_points should be a list of dataframes, a dataframe or a matrix")
  if(!imager::is.cimg(base_img)) stop("Error : base_img must be of type c.img")
  if(!is.list(body_lab_points)){ #if it's a coordinates matrix or a dataframe
    allbodies <- body_lab_points
  } else {
    allbodies <- do.call(rbind, body_lab_points) #aggregate list of presumed matrices/dataframes
  }
  linear_index <- coordsToLinear(allbodies, base_img)#convert pixel coordinate format to linear index format for imager
  L_channel <- length(base_img[,,,1]) #number of pixels in the image (so per channel)
  im_rmb <- base_img
  im_rmb[c(linear_index, linear_index+L_channel, linear_index+L_channel*2)] <- 0 #set to 0 for all 3 channels based on linear index
  if(viz){
    plot(im_rmb)
  }
  return(im_rmb)
}

#' Clean individuals binary image
#'
#' Removes components touching the edges except for the body + removes small spots
#' @param l_img_bin An image as pixset or c.img or a list of images
#' @param centroid Numeric xy coordinates of the body centroid, or a list body centroids
#' @param as_coords A boolean
#' @param px_filter A numerical. Minimum amount pixel to keep label 
#' @export
gCleanBin <- function(l_img_bin, centroid, as_coords=TRUE, px_filter=25){
  if(!is.list(l_img_bin)){
    if(all(is.na(l_img_bin)) | all(is.na(centroid))) { return(NA) }
    if(!imager::is.cimg(l_img_bin) & !imager::is.pixset(l_img_bin)) {
      stop("l_img_bin must be a c.img or list of c.img") #error handling
    }
  } else { #if it is a list, check the content
    if(!imager::is.cimg(l_img_bin[[1]]) & !imager::is.pixset(l_img_bin[[1]])) {
      stop("l_img_bin must be a c.img or list of c.img") #error handling
    }
  }

  if(is.list(l_img_bin)){ #vectorization
    pb <- progress::progress_bar$new(total = length(l_img_bin))
    res <- mapply(function(img_bin, cen){
      pb$tick()
      return( gCleanBin(img_bin, cen, as_coords) )
    }, l_img_bin, centroid, SIMPLIFY=FALSE)
    return(res)
  }
  if(all(is.na(l_img_bin))) return(NA)
  #0 innit
  dimX <- nrow(l_img_bin); dimY <- ncol(l_img_bin)
  labs0 <- imager::label(im=l_img_bin) #label different connected white patches
  labs <- labs0
  labs[which(l_img_bin==0)] <- 0
  #1 Get body label
  central_lab_val <- labs[centroid[1], centroid[2], 1, 1] #central label = corresponding to the body of the individual of interest
  #2 Remove edge labels
  edge_labs <- c(labs[c(1, dimX),,1,1], labs[, c(1, dimY),1,1]) %>% unique #labels with contact to an edge
  edge_labs_cen <- edge_labs[edge_labs != central_lab_val] #keep body
  #3 Remove small spots
  unique_labs <- unique(labs)
  labs_val_r <- unique_labs[!(unique_labs %in% edge_labs_cen)]#remove edge contact
  labs_val_rsp <- labs_val_r[sapply(labs_val_r, function(x) sum(labs==x)) >= px_filter]
  img_res <- imager::as.cimg(labs %in% labs_val_rsp, dim=c(dimX, dimY, 1, 1))
  if(as_coords){ img_res <- img_res %>% imgAsCoords }

  return(img_res)
}

#' Crop individual's body using crop coordinates given by gCrop
#'
#' @param base_img c.img. Base image suitable for waterstrider pipeline
#' @param body_coords a list of body coordinates
#' @param crop_coords a list of 2x2 crop matrices
#'
#' @export
recropBodies <- function(body_coords, base_img, crop_coords){
  empty_base_img <- imager::imfill(dim = dim(base_img), val = 0) #create empty image with base image dimensions
  body_l_coords <- lapply(body_coords, coordsToLinear, img = base_img) #linear coordinates in the base image of body points
  body_l_crops <- mapply(function(body, crop){ #crop them with the same crop as for the legs to preserve coordinates
    empty_base_img[body] <- 1 #only fill the specific individuals point to avoid overlap related issues
    xrange <- crop[,1]; yrange <- crop[,2] #crop range
    x <- NA
    y <- NA #avoid warnings in devtools::check() and goodpractice::gp()
    body_crop <- imager::imsub(empty_base_img, x %inr% xrange, y %inr% yrange) #crop
    return(body_crop)
  }, body_l_coords, crop_coords, SIMPLIFY=FALSE)
  return(body_l_crops)
}

#' Body dilation proportionnal to body length
#'
#' @param body_img Binary pixset or c.img of waterstrider body
#' @param body_length A numerical for length of corresponding body
#' @param dilation_ratio Dilation kernel size as body length ratio
#' 
#' @export
dilBodies <- function(body_img, body_length, dilation_ratio=0.3){
  if(!imager::is.pixset(body_img) && !(imager::is.cimg(body_img))){ #check type 
    if(is.list(body_length) | is.vector(body_length)){
      res <- mapply(function(im, l){
        dilBodies(im, l, dilation_ratio=dilation_ratio)
      }
      , body_img, body_length, SIMPLIFY = FALSE) #vectorization
      return(res)
    } else {
      stop("Wrong formats for body_img or body_length in dilBodies()")
    }
  }
  size0 <- round(body_length * dilation_ratio + 1)  #adaptative dilation kernel based on body length
  size <- ifelse(size0 %% 2 == 0, size0 + 1, size0)  # Round to an odd number to avoid warnings in makeBrush
  kernel <- size %>% make_disc_kernel %>% imager::as.cimg() #create disc kernel suitable for imager based on size
  dil_body <- body_img %>% imager::dilate(mask = kernel) #dilation around body to ensure intersecting legs and not the body
  return(dil_body)
}

#' Gerris orientation angle
#'
#' Orientation based on body points PCA, then antero-posterior point density, or appendages repartition
#'
#' @param body_coords list or single dataframe or matrix of body points coordinates
#' @param intersection_coords list or single dataframe or matrix of 
#' full individual and dilated body contour intersection coordinates
#' @param diag a logical for output analysis diagnostic messages
#' 
#' @export
gOrientation <- function(body_coords, intersection_coords, diag=FALSE){
  if(is.list(body_coords) & is.list(intersection_coords)){
    return(mapply(gOrientation, body_coords, intersection_coords)) #list iteration
  }
  pca <- stats::prcomp(body_coords) #principal direction given by body PCA
  m_rot_PC1 <- pca$rotation[,1] #rotation matrix of PC1
  body_PC1 <- body_coords %*% m_rot_PC1 #body PC1 coordinates
  ang0 <- atan2(m_rot_PC1[2], m_rot_PC1[1]) #angle body elongation angle
  # DENSITY BASED 1D ORIENTATION (females)
  min_PC1 <- min(body_PC1)
  mid_PC1 <- (max(body_PC1)-min_PC1)/2 + min_PC1 #middle distance between rear and front of individual
  sum_right2 <- sum((body_PC1-mid_PC1) > 0) * 2 #number of points above mid_PC1 x2
  l_PC1 <- length(body_PC1) #total amount of points
  dif_PC1 <- abs(l_PC1-sum_right2) / l_PC1# % de différence entre droite et gauche de milieu de la range #diffetence in %
  if(dif_PC1 > 0.1){ #if density significant
    if(diag){message("Antero_posterior density based orientation (female)")}
    ang0 <- ang0 -pi
    if(sum_right2 / l_PC1 > 1){ #check head direction (least point density direction)
      angle <- ang0
    } else {
      angle <- ang0 + sign(ang0) * -pi #keep between -pi and pi
    }
    # INTERSECTION BASED 1D ORIENTATION
  } else {
    if(diag){message("Limbs intersections based orientation")}
    inter_PC1 <- intersection_coords %>% as.matrix %*% m_rot_PC1 #PC1 values of intersections
    ori_PC1 <- mean(inter_PC1 > mid_PC1) > 0.5 #if most intersection points are in the direction of PC1
    angle <- ang0 -pi * !ori_PC1 #angle of c
  }
  names(angle) <- NULL
  return(angle)
}

#' Gerris limbs intersection points
#'
#' For v a boolean of intersection points in the contour, assume limb intersections as continuous sequences
#' of intersection values in the contour. Compute first, last, and middle point of each sequence
#'
#' @param v a vector of logical values
#'
#' @export
gLimbInter <- function(v){
  seq_lim <- boolSeqLim(v, circular=TRUE)#Interval of all True/1 sequence in boolean vector
  f <- seq_lim$first #first point of each sequence
  l <- seq_lim$last #last
  mid <- ((f + l)/2) %>% round #middle point assuming no circularity and matching l & f
  l_t <- length(f) #number of sequences
  if(f[l_t]>l[l_t]){ #circularity issue detection > last before first
    l_cont <- length(v) #n points of contour
    mid[l_t] <- (f[l_t] + ( l[l_t]+l_cont+1-f[l_t] )/2) %% l_cont #circularity handling
  }
  seq_lim$mid <- mid #add to boolSeqLim list
  return(seq_lim)
}

#' Find rear leg insertions
#'
#' Returns rear leg insertions points as xy coordinates from the dilated body contour dilated_body
#' using head orientation angle ori_angle in radians and intersection between dilated body contour and full body
#' as boolean, like $index output of  <-
#'
#' @param ori_angle Numerical value or vector of numericals
#' for orientation angle of individuals in radians
#' @param dil_contour list or single dilated body contour as ordered xy coordinates matrix or dataframe
#' @param inter_index list or single vector of numericals for indexes of intersection points
#' in corresponding dilated contour
#' @param leg_lim_ratio numeric vector of length two. Values between 0 and 1 to define range where
#' legs are expected to be found as ratio of body elongation starting from the back (excluding head appendages and hind body part)
#' 
#' @export
gLegInsertion <- function(ori_angle, dil_contour, inter_index,
                          leg_lim_ratio=c(0.25, 0.6)){
  if(is.list(dil_contour) & is.list(inter_index)){
    res <- mapply(gLegInsertion, ori_angle, dil_contour, inter_index, SIMPLIFY=FALSE)
    return(res)
  }
  if(is.na(ori_angle) | length(unique(inter_index))==1){
    return(list(left=NA, right=NA)) #NA if no angle, no intersection or all intersection
  }
  bary <- apply(dil_contour, 2, mean) #centroid
  #insertion points as intersection sequences
  seq_inter <- gLimbInter(inter_index) #find intersection sequences = limb base points
  limbs_base <- dil_contour[seq_inter$mid,]
  #orient elements
  ang_rot <- -ori_angle #get rotation angle to orient individual
  rot_mat <- matrix(c(cos(ang_rot), -sin(ang_rot), sin(ang_rot), cos(ang_rot)), ncol = 2) #rotation matrix
  limbs_base_r <- limbs_base %>% as.matrix %*% rot_mat
  cont_r <- dil_contour %>% as.matrix %*% rot_mat
  bary_r <-  bary %*% rot_mat
  #remove limbs insertion on rear-most first 25% of elongation (likely not legs)
  cont_r_min <- min(cont_r[,1])
  int_elong_ratio <- (limbs_base_r[,1]-cont_r_min)/(max(cont_r[,1])-cont_r_min) #normalized position of insertion in elongation direction
  valid_limbs_id <- (int_elong_ratio > leg_lim_ratio[1]) & (int_elong_ratio < leg_lim_ratio[2]) #spacial constraint : no body and no antennas
  #find anterior insertion for both sides
  left_limbs_id <- limbs_base_r[,2] > bary_r[,2] #under y=0 on barycenter reference : right side
  inser <- list()
  valid_left <- left_limbs_id & valid_limbs_id
  valid_right <- !left_limbs_id & valid_limbs_id
  inser$left <- limbs_base[valid_left,][limbs_base_r[valid_left, 1] %>% which.min, ]
  inser$right <- limbs_base[valid_right,][limbs_base_r[valid_right,1] %>% which.min,]

  #if(is.null(inser)){ inser <- NA}
  return(inser)
}


#' Leg segmentation
#'
#' Leg segmentation using insertion points as output of gLegInsertion to detect matching limbs
#' first separates limbs by substracting dilated body to full individual
#'
#' @param gerris `cimg` or list of `cimg`
#' @param dilated_body `cimg` or list of `cimg`
#' @param intersection_coords `x,y` matrix or list of `x,y` matrices
#' @param insertions `x,y` matrix or list of `x,y` matrices
#' @export
gLegSeg <- function(gerris, dilated_body, intersection_coords, insertions){
  if(all(is.na(gerris))) return(list(right=NA, left=NA))
  if(is.list(gerris) & is.list(insertions)){
    pb <- progress::progress_bar$new(total=length(gerris))
    res <- mapply(function(a,b,c,d){
      pb$tick()
      return(gLegSeg(a,b,c,d))
    } , gerris, dilated_body, intersection_coords, insertions, SIMPLIFY = FALSE)
    return(res)
  }
    # find limbs
  limbs <- ((dilated_body - gerris) < 0) %>% imager::as.cimg() #full individual minus dilated body
  limbs[coordsToLinear(intersection_coords, limbs)] <- 1 #reattach intersections with body, required for detection
  split_limbs_img <- limbs %>% imager::split_connected() #images of separated limbs
  split_limbs <- lapply(split_limbs_img, imgAsCoords) #as xy coords 
    # get leg points
  leg_split_id <- lapply(insertions, function(inser){ #for each hind leg insertion point :
    tryCatch({
      lapply(split_limbs, function(x){ #for each limb as xy coordinates :
          overlapPoints(inser, x, coords = FALSE) == TRUE #check if limb is connected to hind leg insertion
        }) %>% unlist %>% which #return as index
      }, error = function(e) { NA } #error handling
    )   
  })
  legs <- lapply(leg_split_id, function(id) {
    if (any(is.na(id[1]))) return(NA)
    tryCatch({ split_limbs[id][[1]] }, #get limbs with index
             error = function(e) { NA })  #return NA if error occur (no matching limb)
  })

  return(legs)
}

#' Find knee and ankle of gerris leg
#'
#' Ankle and knee position using angle along the contour of the leg,
#' oriented with the body insertion point
#'
#' @param leg_coords numerical matrix or dataframe or list of those types.
#' leg points coordinates
#' @param insertion numerical matrix or dataframe or list of those types.
#' insertion points coordinates
#' @param inser_thresh Numerical value in range 0:1 for threshold to remove
#' angular points too close to the leg insertion point
#' (value as fraction of half contour length)
#' @param tresh_ankle Numerical value in range 0:1 for threshold to detect if
#' the difference between the two knee-ankle distances is realistically acceptable 
#' (value as fraction of half contour length)
#' @param viz logical. visualization of landmark positionning and peak along contour distance profile
#' @param viz_angle a boolean. For diagnostic plots of leg contour angle variation and peak fitting
#' @param msg logical. message option
#' @param inflexion_pts_range vector of two integers. Tolerated range of inflexion points on leg contour
#' @param knee_diff_thresh numerical in \code{[0, 1]}. 
#' Tolerated ratio between the distances from the insertion point to each detected knee
#' A value close to 0 will ensure knee-insertion distances are similar
#' @param segment_length_range vector of two numericals in \code{[0, 1]}. Range of tolerated leg segment size
#' as ratio of half leg contour distance
#' @param n_splines an integer. Number of splines to fit leg contour angle variation
#' @param search_w an integer. Distance (in number of points in contour) to use to compute leg contour angles
#' 
#' @export
gLegLandmarks <- function(leg_coords,
                          insertion,
                          inser_thresh=0.1,
                          tresh_ankle=0.12,
                          viz=FALSE,
                          viz_angle=FALSE,
                          msg=TRUE,
                          inflexion_pts_range = c(5,7),
                          knee_diff_thresh = 0.2,
                          segment_length_range = c(0.15, 0.6),
                          n_splines = 30,
                          search_w = 6) {
  
  if(is.null(insertion) | all(is.na(insertion)) |
     is.null(leg_coords) | all(is.na(leg_coords))){ # | length(leg_coords)==2 ???
    if(msg) message("FAILED TO RECOGNIZE ADEQUATE TYPES")
    return(NA)
  } #Na handling
  if(!("dim1" %in% names(insertion))){  #Vectorization assuming nesting "left" "right" structure
    if(msg) message("Attempting vectorization with $left and $right list structure")
    pb <- progress::progress_bar$new(total = length(leg_coords))
    res <- mapply(function(leg, ins){
      pb$tick()
      list(
        left = if( !(all(is.na(leg$left))) | !(all(is.na(ins$left))) ){
          gLegLandmarks(leg$left, ins$left, 
                        inser_thresh = inser_thresh,           # Explicit
                        tresh_ankle = tresh_ankle, 
                        viz = viz, 
                        msg = msg,
                        inflexion_pts_range = inflexion_pts_range,
                        knee_diff_thresh = knee_diff_thresh,
                        segment_length_range = segment_length_range,
                        n_splines = n_splines,
                        search_w = search_w)
        } else { NA },
        right = if( !(all(is.na(leg$right))) | !(all(is.na(ins$right))) ){
          gLegLandmarks(leg$right, ins$right, 
                        inser_thresh = inser_thresh,
                        tresh_ankle = tresh_ankle, 
                        viz = viz, 
                        msg = msg,
                        inflexion_pts_range = inflexion_pts_range,
                        knee_diff_thresh = knee_diff_thresh,
                        segment_length_range = segment_length_range,
                        n_splines = n_splines,
                        search_w = search_w)
        } else { NA }
      )
    }, leg_coords, insertion, SIMPLIFY = FALSE)

    return(res)
  }

  #0 - Contour and contour angles
  leg_img <- coordsAsImg(leg_coords, padding=2, return_offset=TRUE) #image conversion for contour algorithm
  leg_cont_offset <- simpleCont(leg_img$img) #contour of image with offset
  cont <- leg_cont_offset + matrix(leg_img$coord_offset[1,],
                                   nrow=nrow(leg_cont_offset), ncol=2, byrow=TRUE) #offset correction
  peaks <- tryCatch({contourTurns(cont, search_w = search_w,
                                  splines_df = n_splines, viz = viz_angle)}, #contour angles
           error = function(e) {
             if(msg) message("contourTurns() failed : ", e$message)
             return(NULL)
           })
  if(all(is.null(peaks))) {
    if(msg) message("no peak found")
    return(NA)
  }
  # 1 - Find insertion
  lcont <- nrow(cont) #n points of contour
  inser_id <- overlapPoints(cont, insertion, coords=FALSE) %>% which #insertion id in leg contour
  reord <- (inser_id:(inser_id + lcont-1)) %% lcont +1#contour index insertion points as origin
  cont_inser <- rbind(cont[reord,], cont[reord[1],]) #cicular contour
  peaks_inser <- sapply(peaks, function(x) which(reord %in% x)) %>% sort
  dist <- sqrt(diff(cont_inser[,1])^2 + diff(cont_inser[,2])^2) #distance between pairs of points
  dist_cont <- sum(dist)
  dist_peaks_str <- c(0,(dist %>% cumsum))[peaks_inser] #distance of each peak to the insertion point
  names(dist_peaks_str) <- peaks_inser #names to track index in contour
  dist_peaks_rev <- c(0, (dist_cont-dist_peaks_str)[-1] %>% rev) #same but counterclockwise
  names(dist_peaks_rev)[1] <- "1"
  dist_peaks <- list(dist_peaks_str, dist_peaks_rev)
  if(length(peaks_inser) >= min(inflexion_pts_range) &&
     length(peaks_inser) <= max(inflexion_pts_range)){ #incoherent number of inflexion points
    check_detection <- TRUE
    #2 Ignore inflexion points too close to the insertion, as they probably correspond to the inflexion
    cont_half_l <- dist_cont/2 #leg contour half length
    filt_id <- lapply(dist_peaks, function(x) (x >= cont_half_l*inser_thresh))#filtered id in dist_peaks & peaks_inser
    #3 Knees
    #TODO REPLACE knee_dist !!!
    filt_cont_pts <- lapply(filt_id, function(y) {
      ids <- y %>% names %>% as.numeric #numeric id
      filt_ids <- ids[y] #valid numeric id
      cont_inser[filt_ids,][1,] #1st corresponding contour point in each direction
    })
    knee_abs_dist <- sapply(filt_cont_pts, function(x){
      (insertion-x)^2 %>% sum %>% sqrt #distance
    })
    #old // mapply(function(x,y) {x[y][1]}, dist_peaks, filt_id)
    #old // if((diff(knee_dist / cont_half_l) %>% abs) < 0.2) {  #if the dif between the two knee-ankle distances is below 20% of half contour length
    knee_dist <- mapply(function(x,y) {x[y][1]}, dist_peaks, filt_id) #for later
    knee_abs_diff <- (knee_abs_dist %>% diff / sum(knee_abs_dist)) %>% abs
    if(knee_abs_diff < knee_diff_thresh) { 
      #4 Ankles
      ankle_dist <- mapply(function(x,y) {x[y][2]}, dist_peaks, filt_id)
      if( (diff(ankle_dist-knee_dist) %>% abs / cont_half_l) > tresh_ankle ){ #if the dif between the two knee-ankle distances is above tresh_ankle% of half contour length
        ankle_dist1 <- mapply(function(x,y,z) {x[y][z]}, dist_peaks, filt_id, c(2,3)) #distance between points 1-3(straight contour) and 1-2(reverse contour)
        ankle_dist2 <- mapply(function(x,y,z) {x[y][z]}, dist_peaks, filt_id, c(3,2)) #distance between points 1-2(straight contour) and 1-3(reverse contour)
        option1 <- (diff(ankle_dist1-knee_dist) %>% abs / cont_half_l) < tresh_ankle #knee-ankle distances is above tresh_ankle% of half contour length
        option2 <- (diff(ankle_dist2-knee_dist) %>% abs / cont_half_l) < tresh_ankle  #knee-ankle distances is above tresh_ankle% of half contour length
        if(is.na(option1)){ option1 <- FALSE }
        if(is.na(option2)){ option2 <- FALSE }
        if(option1) { ankle_dist <- ankle_dist1 }
        else if(option2) { ankle_dist <- ankle_dist2 }
        else {
          if(msg) message("incoherent knee-ankle distance")
          check_detection <- FALSE
        } #incoherent knee-ankle distance
      }
    } else {
      if(msg) message("incoherent knee distances")
      check_detection <- FALSE
    } #incoherent knee distances
  } else {
    if(msg) message("incoherent number of inflexion points")
    check_detection <- FALSE
  } #incoherent number of inflexion points
  if(check_detection){ #sucessful detection of components
    knee_pts <- do.call(rbind, filt_cont_pts)
    ankle_pts <- cont_inser[names(ankle_dist) %>% as.numeric, ]
    knee <- apply(knee_pts,2,mean)
    ankle <- apply(ankle_pts,2,mean)
    segment_lengths <- sapply(list(insertion, ankle), function(pt){  #for safety check
      ((pt -knee)^2 %>% sum %>% abs %>% sqrt)/cont_half_l
    })
    if( all(segment_lengths < min(segment_length_range) |
            segment_lengths > max(segment_length_range)) ){ #incoherent segment sizes relative to contour
      if(msg) message("incoherent segment sizes relative to contour")
      return(NA)
    }
    if(viz){
      par(mfrow=c(1,2))
      #P1 Inflexion points to insertion point distance + insertion threshold
      n_p <- length(dist_peaks_str)
      plot(x = dist_peaks_str, y=rep(1.1,n_p),pch=16,ylim=c(0.2,1.8),
           xlim = c(0,sum(dist)/2), main="Inflexion Profile")
      points(x = c(0,sum(dist)), y = c(1,1), col = "blue",pch=16)
      points(x = dist_peaks_rev, y = rep(0.9,n_p), pch = 16,
             ylim = c(0.5,1.5), xlim = c(0,sum(dist)))
      abline(v = sum(dist)/2*inser_thresh, col = "blue")
      legend("bottomright", legend = c(paste("thresh",inser_thresh),"mid leg"),
             col = c("blue",1), lty = 1, bty = "n")
      #P2 Contour reordered around insertion point + inflexion points
      plot(rbind(cont_inser,cont_inser[1,]), type="l", as=1,
           main="Leg Contour", lwd=4, col="gray40") #reordered cont
      points(cont_inser[peaks_inser,],pch=16,cex=1.5) #reordered peaks
      points(cont_inser[1,],pch=16,col="blue") #insertion as 1st point
      points(cont_inser[1,],pch=13,cex=2,col="blue") #insertion as 1st point
      points(knee_pts,cex=1,pch=16,col="forestgreen")
      lines(knee_pts, lwd=3, col="forestgreen")
      points(ankle_pts,cex=1,pch=16,col="darkorchid4")
      lines(ankle_pts, lwd=3, col="darkorchid4")
      lines(rbind(cont_inser[1,],knee,ankle), lwd=2, lty=4, col="red")
      points(rbind(knee,ankle),pch=16,cex=1,col="red") #insertion as 1st point
    }
    return(data.frame(insertion = cont_inser[1,] %>% unlist, knee = knee, ankle = ankle) %>% t)
  } else { return( NA )}
}

#' Connect leg landmarks with body
#'
#' As gLegLandmarks finds insertion point on dilated body, the insertion landmark is biased
#' Given body as xy coords and landmarks as output of gLegLandmarks, corrects the insertion point
#'
#' @param body dataframe or matrix of body points coordinates 
#' @param landmarks landmarks as output of gLegLandmarks()
#' @param viz logical.visualization option
#' 
#' @export
gConnectLeg <- function(body, landmarks, viz=FALSE){
  #Vectorization assuming left/right inner nesting for landmarks
  if(body %>% is.list){
    pb <- progress::progress_bar$new(total = length(body))
    res <- mapply(function(bod, lmk){
      pb$tick()
      list(
        left = gConnectLeg(bod, lmk$left, viz),
        right = gConnectLeg(bod, lmk$right, viz)
      )
    }, body, landmarks, SIMPLIFY = FALSE)
    return(res)
  }
  if(landmarks %>% is.na %>% all){ return(NA) }
  #Function
  kne <- landmarks["knee",] %>% as.numeric #knee landmark
  ins <- landmarks["insertion",] %>% as.numeric #dilated insertion landmark
  v_dir <- (ins-kne) #direction vector from knee to dilated insertion
  v_u <-  v_dir / (v_dir^2 %>% sum %>% sqrt) #unit direction vector
  new_ins <- ins #new insertion value
  new_ins_round <- new_ins %>% round #new rounded insertion value to check pixel correspondence
  iterations <- 0 #cap the number of iterations
  point_in_set <- function(pt, set) {
    any(set[, 1] == pt[1] & set[, 2] == pt[2])
  }
  while( !point_in_set(new_ins_round, body) && iterations<2000){ #if new point is not in body points
    new_ins <- new_ins + v_u #new point for insertion moving towards body
    new_ins_round <- new_ins %>% round #round for search in discrete pixel body
    iterations <- iterations+1
  }
  new_ins <- (new_ins - v_u) %>% round #avoid overlap
  #visualization
  if(viz){
    plot(rbind(body, landmarks), as=1, col="grey60", pch=15)
    points(landmarks, pch=16)
    lines(landmarks, col=1, lwd=4)
    points(new_ins[1], new_ins[2], pch=16, col="maroon", cex=2)
    lines(rbind(new_ins, landmarks[1,]), col="maroon", lwd=4)
    legend("bottomleft", legend = "Extended Leg", col = "maroon", lwd = 4)
  }
  #output
  if(iterations==2000){ #leg not connecting : most likely a landmarking error
    return( NA )
  } else {
    landmarks["insertion", ] <- new_ins
  }
  return(landmarks)
}

#' Measure tibia and femur
#'
#' Scale from pixel to micrometer, and landmarks as outputs of gConnectLeg or gLegLandmarks
#'
#' @param landmarks landmarks as output of gLegLandmarks()
#' @param scale a numerical value for conversion from pixels to micrometers
#' @export
gMeasureLeg <- function(landmarks, scale){
  if(landmarks[[1]] %>% is.list){
    res <- lapply(landmarks, function(lmk){
      list(
        left = tryCatch({gMeasureLeg(lmk$left, scale)}, error = function(e) { NA }),
        right = tryCatch({gMeasureLeg(lmk$right, scale)}, error = function(e) { NA })
      )})
    return(res)
  }
  if(landmarks %>% is.na %>% all){
    return( NA )
  }
  femur_pix <- eu_dist(landmarks["insertion",], landmarks["knee",])
  tibia_pix <- eu_dist(landmarks["knee",], landmarks["ankle",])
  femur <- femur_pix / scale[1]
  tibia <- tibia_pix / scale[1]
  #output
  res <- c(femur, tibia)
  names(res) <- c("femur", "tibia")
  return(res)
}

#' Individual plot for gPipeline
#'
#' Individual plot for gPipeline, designed to be used on lists of internal 
#' variables of gPipepline()
#'
#' @param i index of the individual to plot in the list-format variables of gPipeline()
#' @param full Similar to the internal variable of the same name in gPipeline()
#' @param body Similar to the internal variable of the same name in gPipeline()
#' @param cen Similar to the internal variable of the same name in gPipeline()
#' @param dilcont Similar to the internal variable of the same name in gPipeline()
#' @param ang Similar to the internal variable of the same name in gPipeline()
#' @param legs Similar to the internal variable of the same name in gPipeline()
#' @param leg_lm Similar to the internal variable of the same name in gPipeline()
#' @param leg_size Similar to the internal variable of the same name in gPipeline()
#' @param inser Similar to the internal variable of the same name in gPipeline()
#' @param clean_base_path Similar to the internal variable of the same name in gPipeline()
#' @export
gGerrisPlot <- function(i, full, body, cen, dilcont, ang, legs,
                        leg_lm, leg_size, inser, clean_base_path){
  legend_err <- c()
  full_coords <- full[[i]] %>% imgAsCoords
  body_coords <- body[[i]]
  par(mar=c(0,0,0,0))
  if(!anyNA(full[[i]])){
    plot(NA,  #empty plot
         xlim = range(full_coords[,1]) + c(-5,5),  # Added padding of 5 units
         ylim = (range(full_coords[,2]) + c(-5,5)) %>% rev,  #invert y to match image format
         asp = 1, xlab = "", ylab = "")
    symbols(full_coords[,1], full_coords[,2], #full binarization
            squares = rep(1, nrow(full_coords)),
            inches = FALSE, add = TRUE, fg = NA, bg = "grey50")
    symbols(body_coords[,1], body_coords[,2], #body
            squares = rep(1, nrow(body_coords)),
            inches = FALSE, add = TRUE, fg = NA, bg = "grey30")
  } else {
    plot(body[[i]], cex=0.1)
    text(x=body[[i]][,1]%>%mean,y=body[[i]][,2]%>%mean,labels='FULL BODY SEGMENTATION FAILED')
  }
  centroid <- cen[[i]]
  if(is.numeric(centroid) && length(centroid)==2){
    points(centroid[1], centroid[2], pch = 16, col = "magenta") #centroid
    points(dilcont[[i]], type = "l", lty = 2, col = "magenta3") #dilated contour
    #orientation
    if(!is.na(ang[[i]])){
      pt_ang <- c(cos(ang[[i]]), sin(ang[[i]])) * 20 + cen[[i]]
      lines(rbind(cen[[i]], pt_ang), lwd = 2, col = "magenta") #orientation line
    }
  }
  #legs
  lapply(legs[[i]], function(x){
    if(!anyNA(x) && length(x) != 2) {  # Only draw if valid
      symbols(x[,1], x[,2],
              squares = rep(1, nrow(x)),
              inches = FALSE, add = TRUE, fg = NA, bg = "gray15")
    }
  }) %>% invisible
  #leg landmarks
  mapply(function(x, size){ #for each side
    if(!anyNA(x) && nrow(x) > 0) {  # Only process if valid
      pal <- c("coral","coral3")
      pal0 <- c(pal[1], 0, pal[2])
      for(j in c(1,3)){ #i is already used in the function
        pt <- x[(j+1)%%3+1, ]
        if(!anyNA(pt)) {  # Check if point is valid
          points(pt[1], pt[2], pch = 16, col = pal0[j])
          lines(x[-j,], col = pal0[j], lwd = 2, lty = 1)
          ix <- j%%3+1
          if(ix==3){ix <- 2}
          size <- round(size)
          off <- 0.2
          ft <- 11
          colo <- "#000000E6" #E6 is for transparency
          posi <- 2
          text(x=pt[1]-off, y=pt[2]-off, font=ft, labels=size[ix], col=colo, pos=posi) #dark outline
          text(x=pt[1]+off, y=pt[2]+off, font=ft, labels=size[ix], col=colo, pos=posi)
          text(x=pt[1]-off, y=pt[2]+off, font=ft, labels=size[ix], col=colo, pos=posi)
          text(x=pt[1]+off, y=pt[2]-off, font=ft, labels=size[ix], col=colo, pos=posi)
          text(x=pt[1], y=pt[2], font=11, labels=size[ix],
               col=pal[ix], pos=posi)
        }
      }
    }
  }, leg_lm[[i]], leg_size[[i]]) %>% invisible
  lapply(inser[[i]], function(x){   # Insertions
    if(!anyNA(x)) {  # Only draw if valid
      points(x, pch = 13, col = "magenta3")
    }
  }) %>% invisible
  legend("bottomleft",
         legend = c(i, clean_base_path),
         bty = "n",
         text.col = c("black", "gray40"),
         cex = c(2, 1),
         text.font = c(11, 11),
         bg = NA)
}

#' Re-format output of gMeasureLeg for easy inclusion in dataframe
#'
#' @param leg_list output of gMeasureLeg()
#' @export
gSplitLegMeasures <- function(leg_list) {
  right_femur <- list()
  right_tibia <- list()
  left_femur <- list()
  left_tibia <- list()
  for (i in seq_along(leg_list)) {
    if (!is.na(leg_list[[i]]$right[1])) {
      right_femur[[i]] <- leg_list[[i]]$right["femur"]
      right_tibia[[i]] <- leg_list[[i]]$right["tibia"]
    } else {
      right_femur[[i]] <- NA
      right_tibia[[i]] <- NA
    }
    if (!is.na(leg_list[[i]]$left[1])) {
      left_femur[[i]] <- leg_list[[i]]$left["femur"]
      left_tibia[[i]] <- leg_list[[i]]$left["tibia"]
    } else {
      left_femur[[i]] <- NA
      left_tibia[[i]] <- NA
    }
  }
  l_out <- list(
    right_femur = right_femur,
    right_tibia = right_tibia,
    left_femur = left_femur,
    left_tibia = left_tibia
  )
  l_out2 <- lapply(l_out, unlist)
  l_out3 <- lapply(l_out2, function(x) stats::setNames(x, seq_along(x)) )
  return(l_out3)
}
