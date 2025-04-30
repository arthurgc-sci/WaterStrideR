# Gerroidea specific functions, used in gPipeline
#---------------------------------------------------------------------------------------------------

#' Fast segmentation with size filter
#'
#' Given a binary image as pixset, returns list of coordinates of every label with pixels in the min-max range
#' @param img_bin a pixset
#' @param min minimum number of pixels to keep label
#' @param max maximum number of pixels to keep label
#' @export
gFastSeg <- function(img_bin, min=300, max=1500){
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
  return(body_lab_points)
}

#' Plot gerris position and id with 1cm scale on original image
#'
#' @param base_img c.img or pixset of base image
#' @param scale mean pixel amount for 1 micron
#' @param x,y ordered position of individuals
#' @export
gDetectionPlot <- function(base_img, scale, x, y){
  par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
  plot(base_img)
  #scale
  x_start <- 100
  y_start <- ncol(base_img)-100
  graphics::segments(x_start, y_start, x_start+scale*10000, y_start, col = "brown4", lwd = 4)
  text(x = (2*x_start+scale[1]*10000)/2, y = y_start-50,
       labels = "1 cm", col = "brown4", cex = 2)
  #insects
  points(x, y, pch=4, col="white", cex=2.5)
  text(x, y+40,
       labels = seq_along(x), col = "cyan", cex = 2.5)
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
cropPoints = function(list_points, list_size, scaling_factor=1){
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
#' @export
gCrop <- function(base_img, centroids, sizes, viz=F, factor=3.5){
  #Error handling
  if(!is.cimg(base_img)) stop("Error : base_img must be of type c.img")
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
    crop <- imager::imsub(base_img, x %inr% xrange, y %inr% yrange)
    crop[which(crop==1)] <- max(crop[which(crop!=1)])
    out <- crop %>% grayscale %>% invert_grayscale #%>% correct_illumination #TODO handle ruler > R(img) negatively impact global segmentation performance
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
#' @export
gNobody <- function(base_img, body_lab_points, viz=F){
  #Error handling
  if(!is.matrix(body_lab_points) & !is.list(body_lab_points) & !is.matrix(body_lab_points)) stop("Error : body_lab_points should be a list of dataframes, a dataframe or a matrix")
  if(!is.cimg(base_img)) stop("Error : base_img must be of type c.img")
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

#' Binarize gerris using GMM derivative
#'
#' Binarize gerris by setting a fixed threshold per input image using GMM derivative
#'
#' @param gerris_crops c.img or list of c.img of cropped binary individuals
#' @param nG number of gaussians to fit in GMM
#' @param noise_n normal noise to add to histogram data, helps breaking plateaus
#' @param y_root see gMMThresh
#' @param msg logical for display of diagnostic messages
#' @export
gGMMThresh <- function(gerris_crops, nG=2, noise_n=0.005, y_root=-25, msg=T){
  if(is.list(gerris_crops)){ #Vectorization
    if(msg) message("GMM thresholding")
    pb <- progress::progress_bar$new(total=length(gerris_crops))
    res <- lapply(gerris_crops, function(x){
      pb$tick()
      gGMMThresh(x, nG, noise_n, y_root)
    })
    return(res)
  }
  g_val <- gerris_crops %>% as.numeric #convert cropped images to numeric
  g_val_clean <- g_val[g_val!=0 & g_val!=1] #removes eventual extreme values created by gNobody
  l_GMM <- GMM(img_val = g_val_clean, n = nG,
               keep_max = F, n_noise_sd = noise_n) #Fit two weighted gaussians to slighlty noised data
  thresh <- GMMdRoots(l_GMM$df, y = y_root, msg = F) #find roots of the GMM for y= -25 (a fixed slope)
  if(length(thresh) == 0) return(NA)
  binarized_gerris <- gerris_crops > max(thresh) #binarize using maximum solution of f'(x)
  return( binarized_gerris )
}

#' Clean individuals binary image
#'
#' Removes components touching the edges except for the body + removes small spots
#' @param l_img_bin An image as pixset or c.img or a list of images
#' @param centroid Numeric xy coordinates of the body centroid, or a list body centroids
#' @param as_coords A boolean
#' @export
gCleanBin <- function(l_img_bin, centroid, as_coords=T){
  if(!is.list(l_img_bin)){
    if(all(is.na(l_img_bin)) | all(is.na(centroid))) { return(NA) }
    if(!is.cimg(l_img_bin) & !is.pixset(l_img_bin)) {
      stop("l_img_bin must be a c.img or list of c.img") #error handling
    }
  } else { #if it is a list, check the content
    if(!is.cimg(l_img_bin[[1]]) & !is.pixset(l_img_bin[[1]])) {
      stop("l_img_bin must be a c.img or list of c.img") #error handling
    }
  }

  if(is.list(l_img_bin)){ #vectorization
    pb <- progress::progress_bar$new(total = length(l_img_bin))
    res <- mapply(function(img_bin, cen){
      pb$tick()
      return( gCleanBin(img_bin, cen, as_coords) )
    }, l_img_bin, centroid)
    return(res)
  }
  if(all(is.na(l_img_bin))) return(NA)
  #0 innit
  dimX <- nrow(l_img_bin); dimY <- ncol(l_img_bin)
  labs0 <- label(im=l_img_bin) #label different connected white patches
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
  labs_val_rsp <- labs_val_r[sapply(labs_val_r, function(x) sum(labs==x)) > 25]
  img_res <- imager::as.cimg(labs %in% labs_val_rsp, dim=c(dimX, dimY, 1, 1))
  if(as_coords){ img_res <- img_res %>% imgAsCoords }

  return(img_res)
}

#' Crop individual's body using crop coordinates given by gCrop
#'
#'base_img : a single c.img / body_coords : list of body coordinates / crop_coords : list of 2x2 crop matrices
#' @export
recropBodies <- function(body_coords, base_img, crop_coords){
  empty_base_img <- imfill(dim = dim(base_img), val = 0) #create empty image with base image dimensions
  body_l_coords <- lapply(body_coords, coordsToLinear, img = base_img) #linear coordinates in the base image of body points
  body_l_crops <- mapply(function(body, crop){ #crop them with the same crop as for the legs to preserve coordinates
    empty_base_img[body] <- 1 #only fill the specific individuals point to avoid overlap related issues
    xrange <- crop[,1]; yrange <- crop[,2] #crop range
    body_crop <- imager::imsub(empty_base_img, x %inr% xrange, y %inr% yrange) #crop
    return(body_crop)
  }, body_l_coords, crop_coords)
  return(body_l_crops)
}

#' Body dilation proportionnal to body length
#'
#' @export
dilBodies <- function(body_img, body_length){
  kernels <- lapply(body_length, function(x){ #adaptative dilation kernel based on body length
    size <- round(x * 1.3 - x + 1)
    size <- ifelse(size %% 2 == 0, size + 1, size)  # Round to an odd number to avoid warnings in makeBrush
    make_disc_kernel(size)
  })
  dil_body <- mapply(function(crop, kern){
    imager::dilate(crop, kern=kern) #TODO TO TEST, REPLACED EBI
  }, body_img, kernels) #dilation around body to ensure intersecting legs and not the body
  return(dil_body)
}

#' Gerris orientation angle
#'
#' Orientation based on body points PCA, then antero-posterior point density, or appendages repartition
#'
#' @export
gOrientation <- function(body_coords, intersection_coords, diag=F){
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


#' Accessory logicals for boolSeqLim
#'
bool_first <- function(x0,x1){
  if(x0==0 & x1==1){T} else {F}
}
bool_last <- function(x0,x1){
  if(x0==1 & x1==0){T} else {F}
}
#' Interval of all True/1 sequence in boolean vector
#'
#' Give matching index of first and last True values of sequences of True in v
#'
#' @export
boolSeqLim <- function(v, circular = T){ #circular :
  l_v <- length(v)
  if(circular){
    v0 <- c(v[l_v],v) # v
    v1 <- c(v,v[1]) # v+1
  } else {
    print('Non circular boolSeqLim non implemented yet')
    break()
  }
  first <- mapply(bool_first, v0, v1)[-(l_v+1)] #which 0 are followed by 1
  last <- mapply(bool_last, v0, v1)[-1] #which 1 are followed by 0
  i_first <- which(first) #as index
  i_last <- which(last)
  if(i_last[1] < i_first[1]){ #match both lists
    i_last <- c(i_last[-1], i_last[1]) #first id in last position
  }
  return(list("first" = i_first, "last" = i_last))
}

#' Gerris limbs intersection points
#'
#' For v a boolean of intersection points in the contour, assume limb intersections as continuous sequences
#' of intersection values in the contour. Compute first, last, and middle point of each sequence
#'
#' @export
gLimbInter <- function(v){
  seq_lim <- boolSeqLim(v, circular=T)#Interval of all True/1 sequence in boolean vector
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
#' @export
gLegInsertion <- function(ori_angle, dil_contour, inter_index){
  if(is.list(dil_contour) & is.list(inter_index)){
    res <- mapply(gLegInsertion, ori_angle, dil_contour, inter_index, SIMPLIFY=F)
    return(res)
  }
  if(is.na(ori_angle) | sum(inter_index)==0) return(list(left=NA, right=NA)) #NA if no angle or no intersection
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
  valid_limbs_id <- (int_elong_ratio > 0.25) & (int_elong_ratio < 0.6) #spacial constraint : no body and no antennas
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
#'
#' @export
gLegSeg <- function(gerris, dilated_body, intersection_coords, insertions){
  if(all(is.na(gerris))) return(list(right=NA, left=NA))
  if(is.list(gerris) & is.list(insertions)){
    pb <- progress::progress_bar$new(total=length(gerris))
    res <- mapply(function(a,b,c,d){
      pb$tick()
      return(gLegSeg(a,b,c,d))
    } , gerris, dilated_body, intersection_coords, insertions, SIMPLIFY = F)
    return(res)
  }
    # find limbs
  limbs <- ((dilated_body - gerris) < 0) %>% imager::as.cimg #full individual minus dilated body
  limbs[coordsToLinear(intersection_coords, limbs)] <- 1 #reattach intersections with body, required for detection
  split_limbs <- lapply((imager::split_connected(limbs)), imgAsCoords) #list xy coords of separated limbs
    # get leg points
  leg_split_id <- lapply(insertions, function(inser){ #for each hind leg insertion point :
    tryCatch({
      lapply(split_limbs, function(x){ #for each limb as xy coordinates :
        overlapPoints(inser, x, coords=F) == T #check if limb is connected to hind leg insertion
      })}, error = function(e) { NA } #error handling
    )  %>% unlist %>% which #return as index
  })
  legs <- lapply(leg_split_id, function(id) {
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
#' @export
gLegLandmarks <- function(leg_coords, insertion, inser_thresh=0.1, tresh_ankle=0.12, viz=F, msg=T){
  if(is.null(insertion) | all(is.na(insertion)) | is.null(leg_coords) | all(is.na(leg_coords))){
    return(NA)
  } #Na handling
  if(!("dim1" %in% names(insertion))){  #Vectorization assuming nesting "left" "right" structure
    pb <- progress::progress_bar$new(total = length(leg_coords))
    res <- mapply(function(leg, ins){
      pb$tick()
      list(
        left = if( !(all(is.na(leg$left))) | !(all(is.na(ins$left))) ){
          gLegLandmarks(leg$left, ins$left, inser_thresh, tresh_ankle, viz, msg)
        } else { NA } ,
        right = if( !(all(is.na(leg$right))) | !(all(is.na(ins$right))) ){
          gLegLandmarks(leg$right, ins$right, inser_thresh, tresh_ankle, viz, msg)
        } else { NA }
      )
    }, leg_coords, insertion, SIMPLIFY = FALSE)

    return(res)
  }

  #0 - Contour and contour angles
  leg_img <- coordsAsImg(leg_coords, padding=2, return_offset=T) #image conversion for contour algorithm
  leg_cont_offset <- simpleCont(leg_img$img) #contour of image with offset
  cont <- leg_cont_offset + matrix(leg_img$coord_offset[1,], nrow=nrow(leg_cont_offset), ncol=2, byrow=TRUE) #offset correction
  peaks <- tryCatch({contourTurns(cont, search_w = 6, splines_df = 30)}, #contour angles
           error = function(e) {
             if(msg) message("contourTurns() failed : ",e$message)
             return(NULL)
           })
  if(all(is.null(peaks))) { return(NA) }
  # 1 - Find insertion
  lcont <- nrow(cont) #n points of contour
  inser_id <- overlapPoints(cont, insertion, coords=F) %>% which #insertion id in leg contour
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
  if(length(peaks_inser) < 8 && length(peaks_inser) > 4){ #incoherent number of inflexion points
    check_detection <- T
    #2 Ignore inflexion points too close to the insertion, as they probably correspond to the inflexion
    cont_half_l <- dist_cont/2 #leg contour half length
    filt_id <- lapply(dist_peaks, function(x) !(x <= cont_half_l*inser_thresh))#filtered id in dist_peaks & peaks_inser
    #3 Knees
    #TODO REMPLACE knee_dist !!!
    filt_cont_pts <- lapply(filt_id, function(y) {
      ids <- y %>% names %>% as.numeric #numeric id
      filt_ids <- ids[y] #valid numeric id
      cont_inser[filt_ids,][1,] #corresponding contour point
    })
    knee_abs_dist <- sapply(filt_cont_pts, function(x){
      (insertion-x)^2 %>% sum %>% sqrt #distance
    })
    #old // mapply(function(x,y) {x[y][1]}, dist_peaks, filt_id)
    #old // if((diff(knee_dist / cont_half_l) %>% abs) < 0.2) { #if the dif between the two knee-ankle distances is below 20% of half contour length
    knee_dist <- mapply(function(x,y) {x[y][1]}, dist_peaks, filt_id) #for later
    knee_abs_diff <- (knee_abs_dist %>% diff / sum(knee_abs_dist)) %>% abs
    if(knee_abs_diff < 0.2) {
      #4 Ankles
      ankle_dist <- mapply(function(x,y) {x[y][2]}, dist_peaks, filt_id)
      if( (diff(ankle_dist-knee_dist) %>% abs / cont_half_l) > tresh_ankle ){ #if the dif between the two knee-ankle distances is above tresh_ankle% of half contour length
        ankle_dist1 <- mapply(function(x,y,z) {x[y][z]}, dist_peaks, filt_id, c(2,3)) #distance between points 1-3(straight contour) and 1-2(reverse contour)
        ankle_dist2 <- mapply(function(x,y,z) {x[y][z]}, dist_peaks, filt_id, c(3,2)) #distance between points 1-2(straight contour) and 1-3(reverse contour)
        option1 <- (diff(ankle_dist1-knee_dist) %>% abs / cont_half_l) < tresh_ankle #knee-ankle distances is above tresh_ankle% of half contour length
        option2 <- (diff(ankle_dist2-knee_dist) %>% abs / cont_half_l) < tresh_ankle  #knee-ankle distances is above tresh_ankle% of half contour length
        if(is.na(option1)){ option1 <- F }
        if(is.na(option2)){ option2 <- F }
        if(option1) { ankle_dist <- ankle_dist1 }
        else if(option2) { ankle_dist <- ankle_dist2 }
        else { check_detection <- F } #incoherent knee-ankle distance
      }
    } else { check_detection <- F } #incoherent knee distances
  } else { check_detection <- F } #incoherent number of inflexion points
  if(check_detection){ #sucessful detection of components
    knee_pts <- do.call(rbind, filt_cont_pts)
    ankle_pts <- cont_inser[names(ankle_dist) %>% as.numeric, ]
    knee <- apply(knee_pts,2,mean)
    ankle <- apply(ankle_pts,2,mean)
    segment_lengths <- sapply(list(insertion, ankle), function(pt){  #for safety check
      ((pt -knee)^2 %>% sum %>% abs %>% sqrt)/cont_half_l
    })
    if( all(segment_lengths < 0.15 | segment_lengths > 0.6) ){ #incoherent segment sizes relative to contour
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
      abline(v = sum(dist)/2*thresh, col = "blue")
      legend("bottomright", legend = c(paste("thresh",thresh),"mid leg"),
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
#' @export
gConnectLeg <- function(body, landmarks, viz=F){
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
  kne <- landmarks["knee",] #knee landmark
  ins <- landmarks["insertion",] #dilated insertion landmark
  v_dir <- (ins-kne) #direction vector from knee to dilated insertion
  v_u <-  v_dir / (v_dir^2 %>% sum %>% sqrt) #unit direction vector
  new_ins <- ins #new insertion value
  new_ins_round <- new_ins %>% round #new rounded insertion value to check pixel correspondence
  iterations <- 0 #cap the number of iterations
  while( !overlapPoints(new_ins_round, body, coords=F) && iterations<100){ #if new point is not in body points
    new_ins <- new_ins + v_u #new point for insertion moving towards body
    new_ins_round <- new_ins %>% round #round for search in discrete pixel body
    iterations <- iterations+1
  }
  new_ins <- new_ins - v_u  #avoid overlap
  #visualization
  if(viz){
    plot(rbind(body, landmarks), as=1, col="grey60", pch=15)
    points(landmarks, pch=16)
    lines(landmarks, col=1, lwd=4)
    points(new_ins[1], new_ins[2], pch=16, col="maroon", cex=2)
    lines(rbind(new_ins, landmarks[1,]), col="maroon", lwd=4)
    legend("bottomright", legend = "Extended Leg", col = "maroon", lwd = 4)
  }
  #output
  if(iterations==100){ #leg not connecting : most likely a landmarking error
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
         ylim = range(full_coords[,2]) + c(-5,5),  # Added padding of 5 units
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
    if(!anyNA(x) && nrow(x) > 0) {  # Only draw if valid
      symbols(x[,1], x[,2],
              squares = rep(1, nrow(x)),
              inches = FALSE, add = TRUE, fg = NA, bg = "gray15")
    }
  }) %>% invisible
  #leg landmarks
  mapply(function(x, size){ #for each side
    if(!anyNA(x) && nrow(x) > 0) {  # Only process if valid
      for(i in c(1,3)){
        pt <- x[(i+1)%%3+1, ]
        if(!anyNA(pt)) {  # Check if point is valid
          points(pt[1], pt[2], pch = 16, col = c("coral3",0,"coral")[i])
          lines(x[-i,], col = c("coral3",0,"coral")[i], lwd = 2, lty = 1)
          ix <- i
          if(ix==3){ix <- 2}
          size <- round(size, 1)
          off <- 0.2
          ft <- 11
          colo <- alpha("black", 0.9)
          posi <- 2
          text(x=pt[1]-off, y=pt[2]-off, font=ft, labels=size[ix], col=colo, pos=posi)
          text(x=pt[1]+off, y=pt[2]+off, font=ft, labels=size[ix], col=colo, pos=posi)
          text(x=pt[1]-off, y=pt[2]+off, font=ft, labels=size[ix], col=colo, pos=posi)
          text(x=pt[1]+off, y=pt[2]-off, font=ft, labels=size[ix], col=colo, pos=posi)
          text(x=pt[1], y=pt[2], font=11, labels=size[ix], col=c("coral3","coral")[ix], pos=posi)
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
  l_out3 <- lapply(l_out2, function(x) setNames(x, seq_along(x)) )
  return(l_out3)
}


#' PCA scores of EFA of oriented contours from binary image
#'
#' @param img_bin a list of binary image or pixset with a single white object
#' @param ori_angle a list of orientation angles of provided shapes
#' @param viz visualization option, additionally returns 3 plots
#' @export
scoresEFA <- function(img_bin, ori_angle=rep(0,length(img_bin)), nb_h=7 ,viz=F){
  na_ang <- ori_angle %>% is.na #valid angle boolean
  if(sum(na_ang) != 0){ #if there is na angles
    ori_angle_c <- ori_angle[!na_ang] #analysis on valid data only
    img_bin_c <- img_bin[!na_ang]
  } else {
    ori_angle_c <- ori_angle
    img_bin_c <- img_bin
  }
  body_cont <- img_bin_c %>% simpleCont #cleaned body contour
  body_cont_mat <- lapply(body_cont, function(x){ #format contours for Momocs
    x <- x[c(1:nrow(x),1),] #close contour
    x <- x %>% as.matrix #as matrix
    return(x)
  })
  ang_mat <- lapply(ori_angle_c, function(a){ #rotation matrices using detected orientation
    cos_a <- cos(a)
    sin_a <- sin(a)
    matrix(c(cos_a, sin_a, -sin_a, cos_a), ncol=2)
  })
  rot_cont <- mapply(function(c, m){ c %*% m }, body_cont_mat, ang_mat) #rotate
  cent_cont <- lapply(rot_cont ,centerContour) #center
  outs <- Momocs::Out(cent_cont) #outline Momocs object
  scale_cont <- outs %>% Momocs::coo_scale #scale
  efa <- scale_cont %>% Momocs::efourier(nb.h=nb_h, norm=T) #elliptic fourier analysis
  efa_pca <- prcomp(efa$coe) #principal components analysis
  if(viz){ #visualization
    par(mfrow=c(1,1))
    scale_cont %>% stack #outlines superposition
    scale_cont %>% panel #all outlines
    scale_cont %>% Momocs::calibrate_harmonicpower_efourier(nb.h=nb_h) #elliptic fourier calibration
    plot_PCA(Momocs::PCA(efa))
  }
  #re-adding missing lines
  scores <- efa_pca$x[,1:20]
  df_scores <- data.frame(matrix(NA, nrow = length(ori_angle)
                          , ncol = ncol(scores)), row.names = names(ori_angle))
  colnames(df_scores) <- colnames(scores)
  df_scores[!na_ang,] <- scores

  return(df_scores)
}

#' Predict factor using LDA model on harmonic fourier principal components
#'
#' @param LDA_model LDA model created with MASS, trained on PCA output
#' @param hPCA_scores harmonic fourier principal components as output of prcomp
#' @export
gPredictLDA <- function(LDA_model, hPCA_scores){
  if(class(LDA_model) != "lda") stop("LDA_model should be an object of class lda")

  valid_bool <- !(hPCA_scores[,1] %>% is.na) #non NA rows
  valid_hPCA <- hPCA_scores[valid_bool,] #na.omit(hPCA)
  n_PC <- LDA_model$means %>% ncol #number of PC expected by the model
  predictions <- LDA_model %>% predict(valid_hPCA[,1:n_PC]) #predictions for valid data
  #add NA rows for matching format between input and output
  p_pred0 <- predictions$posterior #extracting posterior probas and predicted class
  class_pred0 <- predictions$class
  l <- valid_bool %>% length #original length
  p_pred <- matrix(NA, nrow=l, ncol=ncol(p_pred0)) #empty objects with input size
  colnames(p_pred) <- colnames(p_pred0)
  class_pred <- NA %>% rep(l)
  p_pred[valid_bool,] <- p_pred0 #write new data on objects
  class_pred[valid_bool] <- class_pred0 %>% as.character #return classes names

  return(list(posterior = p_pred, class = class_pred %>% as.factor)) #same output formats as predict()
}


#' Predict presence of wing using LDAs
#'
#' @param PCA_scores matrix or dataframe with PC scores for individuals
#' @param sex_class factor with F & M levels of individuals
#' @param LDA_f MASS lda class object
#' @param LDA_m MASS lda class object
#' @param wing_names a vector of length 2 indicating how the presence of wings should be written
#' @export
gPredictWing <- function(PCA_scores, sex_class, LDA_f, LDA_m, wing_names = c(1,2)){
  if(length(sex_class) != nrow(PCA_scores)) stop("PCA_scores row number and sex_class length do not match")
  if(ncol(LDA_m$means) > ncol(PCA_scores)) stop("LDA model expects more variables than provided by PCA_scores")
  if(length(wing_names) != 2) stop("wing_names must be a vector of length 2 c('non-winged', 'winged')")
  #NA results for entire dataset
  l <- length(sex_class) #number of individuals
  post <- matrix(NA , nrow=l, ncol=2) #empty posterior probabilities
  class <- character(l) #empty predicted class
  #predictions for each group using their respective model
  female_id <-  (sex_class == "F") %>% which #female id
  if(length(female_id) > 0){ #if there is at least one female
    wing_f_prediction <- gPredictLDA(LDA_f, PCA_scores[female_id, , drop = FALSE]) #LDA trained on 10PC of 7EF
    post[female_id,] <- wing_f_prediction$posterior #post probas
    class[female_id] <- wing_f_prediction$class %>% as.character #predicted class
  }
  male_id <- (sex_class == "M") %>% which #male id
  if(length(male_id) > 0){ #if there is at least one female
    wing_m_prediction <- gPredictLDA(LDA_m, PCA_scores[male_id, , drop = FALSE]) #LDA trained on 10PC of 7EF
    post[male_id,] <- wing_m_prediction$posterior #post probas
    class[male_id] <- wing_m_prediction$class %>% as.character  #predicted class
  }
  if(!identical(wing_names, c(1,2))){ #update winged class names if needed
    colnames(post) <- wing_names
    class[class==1] <- wing_names[1]
    class[class==2] <- wing_names[2]
  }

  return(list(posterior = post, class = class %>% factor))
}

#' The main thing... for now
#'
#' FAIRE VERSION QUI TOURNE POUR UN DOSSIER D'IMAGES AVEC SORTIE UNIQUE
#'
#' @param img_path image path or list of image paths
#' @param write_df a boolean, to save metrics as csv
#' @param write_plots a boolean, to save image analysis diagnostic plots
#' @param df_output a boolean, to return metrics in R
#' @export
gPipeline <- function(img_path, write_df=T, write_plots=T, df_output=T){
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
  message("|"); message("2 - Computing scale")
  #Scale
  scale <- redRulerScale(base_img, msg=F) #get scale
  message("|"); message("3 - Computing body")
  body_length <- body_length_pix/scale[1] #conversion in mm
  error_margin <- body_length*scale[2]/scale[1] #scale error margin (ignoring pixel error)
  #clean bodies
  clean_body_points <- cleanBodyShape(body_lab_points)
  #Remove bodies
  im_nobodies <- gNobody(base_img = base_img, body_lab_points = clean_body_points, viz=F)
  #Size based ind crop
  l_crop <- gCrop(im_nobodies, body_centroids, body_length_pix, viz=F)
  #Fit GMM to image values
  message("|"); message("4 - Individual thresholding")
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

  message("|"); message("5 - Orienting individuals")
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

  message("|"); message("6 - Leg segmentation")
  #5 - Hind legs insertions
  inser <- gLegInsertion(ori_angle=ang, dil_contour=dilcont, inter_index=inter_idx)
  #6 - Hind legs segmentation
  legs <- gLegSeg(gerris = full,
                  dilated_body = dilbody,
                  intersection_coords = inter_coords,
                  insertions = inser)
  message("|"); message("7 - Leg landmarking & measurement")
  #7 - Get Points of interest from leg
  leg_lm0 <- gLegLandmarks(legs, inser, viz=F, msg=F)
  #8 - Correct insertion landmark by connection leg to body
  leg_lm <- gConnectLeg(body, leg_lm0)
  #9 - Distances
  leg_size <- gMeasureLeg(leg_lm, scale) #leg segment sizes in microns
  #Sex and wing prediction using body contour
  message("|");  message("8 - Sex and wing prediction")
  #elliptic fourier harmonics reduced with PCA
  hPCA <- scoresEFA(body_l_crops, ang, nb_h=7, viz=F)
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
  if(write_plots || write_df){
    setwd(dirname(img_path))
    create_dir("out", msg=F)
    local_dir <- paste0("out/", clean_base_path)
    create_dir(local_dir, msg=F)
  }
  if(write_plots){
    # Detection plot of the image
    png(paste0(local_dir,"/detection_",clean_base_path,".png"), width = nrow(base_img), height = ncol(base_img))  # specify filename and dimensions
    x_cen <- sapply(body_centroids, function(x) x[1])
    y_cen <- sapply(body_centroids, function(x) x[2])
    gDetectionPlot(base_img=base_img, scale=scale[1], x=x_cen, y=y_cen)
    dev.off()
    #Individuals metrics plot
    iplot_dir <- paste0(local_dir,"/individuals_",clean_base_path)
    create_dir(iplot_dir, msg=F)
    for(i in 1:length(body_l_crops)){
      grDevices::png(paste0(iplot_dir,"/",clean_base_path,"_",i,".png"), width=480, height=480)  # specify filename and dimensions
      gGerrisPlot(i, full, body, cen, dilcont, ang, legs,
                  leg_lm, leg_size, inser, clean_base_path)
      dev.off()
    }
  }
  if(write_df){
    write.table(df_out, paste0(local_dir,"/data_",clean_base_path,".csv"), row.names = FALSE, sep=";") #export
  }
  if(df_output){
    return(df_out)
  }
}
