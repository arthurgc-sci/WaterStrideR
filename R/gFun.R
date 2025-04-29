###INSTALL AND LOAD PACKAGES
list_packages = c('dplyr', 'ggplot2', 'magrittr', 'cluster', 'imager', 'mclust', 'purrr',
                  'Momocs', 'pracma', 'ade4', 'BiocManager', 'MASS', 'progress', 'rootSolve')
new_pack = list_packages[!(list_packages %in% installed.packages()[,"Package"])]
if(length(new_pack)){install.packages(new_pack)}
for(pkg in list_packages) {library(pkg, character.only = TRUE)}
if (!requireNamespace('EBImage', quietly = TRUE)) {
  BiocManager::install('EBImage')
}
library('EBImage')


#IMAGE#####################################################################################################
####
# Utility functions for imager::cimg objects.
####

#' DataFrame conversion
#'
#' Convert image to dataframe and expose color channel.
#'
#' @param img an imager::cimg
#'
#' @returns a data.frame
img_to_df <- function(img) {
  out <- img %>%
    as.data.frame() %>%
    mutate(channel = factor(cc, labels = c("R", "G", "B")))
  return(out)
}


#' Histogram equalization
#'
#' Flatten histogram by replacing the pixel value of an image by their rank.
#'
#' @param img an imager::cimg
#'
#' @returns an imager::cimg
hist_eq <- function(img) {
  return(as.cimg(ecdf(img)(img), dim = dim(img)))
}


#' Enhance contrast
#'
#' Enhance the contrasts of an image by running an histogram equalization
#' separately on each channel and combining the results.
#'
#' @param img an imager::cimg
#'
#' @returns an imager::cimg
enhance_contrast <- function(img) {
  out <- img %>%
    imsplit("cc") %>%
    map_il(hist_eq) %>%
    imappend("cc")
  return(out)
}


#' Reduce saturation
#'
#' Reduce the saturation of an image through HSV conversion.
#'
#' @param img an imager::cimg
#' @param ratio an integer, how much to divide the saturation by.
#'
#' @returns an imager::cimg
desaturation <- function(img, ratio = 2L) {
  out <- img %>%
    RGBtoHSV() %>%
    imsplit("cc") %>%
    modify_at(2L, ~ . / ratio) %>%
    imappend("cc") %>%
    HSVtoRGB()
  return(out)
}


#' Correct illumination
#'
#' Correct a gray-scaled image illumination by fitting a linear model and
#' removing the spatial trend.
#'
#' @param img an imager::cimg
#' @param nsamples an integer, pixel subsampling value.
#'
#' @returns an imager::cimg object
correct_illumination <- function(img, nsamples = 1e4L) {
  # convert to grayscale if needed
  if (rev(dim(img))[1L] > 1L) {
    img <- grayscale(img)
  }
  # linear regression trend
  trend <- img %>%
    as.data.frame() %>%
    sample_n(nsamples) %>%
    lm(value ~ x * y, data = .) %>%
    predict(img)
  out <- img - trend
  return(out)
}


#' Invert grayscale image
#'
#' @param img an imager::cimg
#'
#' @returns an imager::cimg
invert_grayscale <- function(img) {
  # convert to grayscale if needed
  if (rev(dim(img))[1L] > 1L) {
    img <- grayscale(img)
  }
  out <- max(img) - img
  return(out)
}


#' Binarize
#'
#' Binarize a grayscale image.
#'
#' @param img an imager::cimg
#' @param quantile a real, the quantile level used for thresholding.
#' @param G an integer, the number of mixture components.
#' @param ... mclust::densityMclust parameters
#'
#' @returns an imager::cimg
binarize <- function(img, quantile = 0.95, G = 1L, ...) {
  # convert to grayscale if needed and invert
  if (rev(dim(img))[1L] > 1L) {
    stop("A grayscale image is expected.")
  }
  # sample
  sample <- sample_histogram(img)
  # fit Gaussian mixture
  gm <- densityMclust(sample, G = G, plot = FALSE, ...)
  # threshold based on 95% quantile
  threshold <- qnorm(
    quantile, gm$parameters$mean[1L], sqrt(gm$parameters$variance$sigmasq[1L])
  )
  out <- img > threshold
  return(out)
}


#PIXSET#####################################################################################################
####
# Utility functions for imager::pixset objects.
####
#' Get centroid of a pixset.
#'
#' @param pixset an imager::pixset
#'
#' @returns an integer vector
get_centroid <- function(pixset) {
  centroid <- pixset %>%
    where() %>%
    colMeans() %>%
    as.integer()
  return(centroid)
}


#' Center a pixset.
#'
#' @param pixset an imager::pixset
#'
#' @returns an imager::pixset
px_center <- function(pixset) {
  centroid <- get_centroid(pixset)
  delta <- dim(pixset)[1L:2L] / 2L - centroid
  out <- pixset %>%
    as.cimg() %>%
    imshift(delta_x = delta[1L], delta_y = delta[2L]) %>%
    as.pixset()
  return(out)
}


#' Intersect borders
#'
#' Test whether a pixset intersect an image borders.
#'
#' @param pixset an imager::pixset
#' @param img an imager::cimg
#'
#' @returns a vector of boolean
intersect_borders <- function(pixset, img) {
  # extract the image borders (as pixels coordinates)
  borders <- img %>%
    px.borders() %>%
    where()
  # count the number of pixels intersecting the borders
  pix_on_borders <- suppressMessages(
    pixset %>%
      where() %>%
      inner_join(borders) %>%
      nrow()
  )
  return(pix_on_borders > 0L)
}


#' Combine pixset
#'
#' Combine pixsets into a reference pixset according to a shared neighbourhood.
#' The neighbourhood is defined as a squared stencil centered around each
#' candidate pixset.
#'
#' @param pixset_ref an imager::pixset, the reference.
#' @param pixset_list a list of imager::pixset, the list of pixsets to combine
#' with the reference.
#' @param padding an integer, the stencil length
#'
#' @returns a vector of boolean
combine <- function(pixset_ref, pixset_list, padding = 50L) {
  # convert candidat to list if there is only one
  if (!is.list(pixset_list) & is.pixset(pixset_list)) {
    pixset_list <- list(pixset_list)
  }
  # iterate over parts
  for (idx in seq_along(pixset_list)) {
    # get centroid
    centroid <- get_centroid(pixset_list[[idx]])
    # define stencil for neighborhood
    stencil <- square_stencil(centroid, padding, pixset_ref)
    # check overlap
    overlap <- pixset_ref %>%
      get.stencil(stencil, x = centroid[1L], y = centroid[2L]) %>%
      sum()
    # merge if overlap
    if (overlap > 0L) {
      pixset_ref <- parany(list(pixset_ref, pixset_list[[idx]]))
    }
  }
  return(pixset_ref)
}

combine_bis <- function(pixset_ref, pixset_list, padding = 25L) {
  # convert candidat to list if there is only one
  if (!is.list(pixset_list) & is.pixset(pixset_list)) {
    pixset_list <- list(pixset_list)
  }
  # get centroid
  centroid <- get_centroid(pixset_ref)
  # define stencil for neighborhood
  stencil <- square_stencil(centroid, padding, pixset_ref)
  # iterate over parts
  for (idx in seq_along(pixset_list)) {
    # check overlap
    overlap <- pixset_list[[idx]] %>%
      get.stencil(stencil, x = centroid[1L], y = centroid[2L]) %>%
      sum()
    # merge if overlap
    if (overlap > 0L) {
      pixset_ref <- parany(list(pixset_ref, pixset_list[[idx]]))
    }
  }
  return(pixset_ref)
}


#' Rotation angle
#'
#' Get the rotation angle to align the abdomen along the x-axis. The function
#' fits an ellipsoid hull around the pixset to derive the rotation angle.
#'
#' @param pixset an imager::pixset
#' @param size an integer, the morphological opening factor (optional)
#'
#' @returns a real
rotation_angle <- function(pixset, size = 6L) {
  # compute ellipsoid hull
  eh <- ehull(pixset, size)
  # get semi-axis lengths
  l_term <- (eh$cov[1L, 1L] + eh$cov[2L, 2L]) / 2L
  r_term <- sum(
    c((eh$cov[1L, 1L] - eh$cov[2L, 2L]) / 2L, eh$cov[1L, 2L]) ** 2L
  )
  lambda_1 <- l_term + sqrt(r_term)
  lambda_2 <- l_term - sqrt(r_term)
  axis_l <- c(sqrt(lambda_1), sqrt(lambda_2))
  # get rotation angle (from x-axis)
  if (eh$cov[1L, 2L] != 0L) {
    theta <- pi - atan2(lambda_1 - eh$cov[1L, 1L], eh$cov[1L, 2L])
  } else {
    theta <- ifelse(eh$cov[1L, 1L] >= eh$cov[2L, 2L], 0L, pi / 2L)
  }
  return(theta)
}



#DRAW#####################################################################################################
####
# Plotting functions.
####
#' Channel histogram
#'
#' Histogram ggplot of each color channel of an image.
#'
#' @param img an imager::cimg
#' @param ch a character vector, subset of ("R", "G", "B").
#'
#' @returns a ggplot2 graph
channel_hist <- function(img, ch = c("R", "G", "B")) {
  if (all(ch %in% c("R", "G", "B"))) {
    gg <- img %>%
      img_to_df() %>%
      filter(channel %in% ch) %>%
      ggplot(aes(value, col = channel)) +
      geom_histogram(bins = 30L) +
      facet_wrap(~channel)
  } else {
    stop("Invalid 'cc' argument.")
  }
  return(gg)
}


#' Raster channel
#'
#' Raster ggplot of an image channel.
#'
#' @param img an imager::cimg
#' @param ch character, specify the image channel, one of "R", "G", "B".
#'
#' @returns a ggplot2 graph
raster_channel <- function(img, ch = NULL) {
  if (ch %in% c("R", "G", "B")) {
    gg <- img %>%
      img_to_df() %>%
      filter(channel == ch) %>%
      ggplot(aes(x, y, fill = value)) +
      geom_raster() +
      scale_x_continuous(expand = c(0L, 0L)) +
      scale_y_continuous(expand = c(0L, 0L), trans = scales::reverse_trans()) +
      scale_fill_viridis_c(direction = 1L, option = "plasma")
  } else {
    stop("Invalid 'cc' argument.")
  }

  return(gg)
}



#UTILS#####################################################################################################
####
# Additional utility functions.
####
#' Sample pixel values from histogram
#'
#' @param obj an array-like object (img, pixset)
#' @param ratio a real, the sampling ratio.
#'
#' @returns an integer vector
sample_histogram <- function(obj, ratio = 0.02) {
  n_samples <- min(1e4L, prod(dim(obj)[1L:2L]))
  n_breaks <- as.integer(n_samples * ratio)
  hist <- hist(obj, breaks = n_breaks, plot = FALSE)
  bins <- with(
    hist, sample(length(mids), n_samples, p = density, replace = TRUE)
  )
  sample <- runif(
    length(bins), hist$breaks[bins], hist$breaks[bins + 1L]
  )
  return(sample)
}


#' Square stencil
#'
#' Define a square stencil within the boundaries of an image cropping points
#' that are out of bounds.
#'
#' @param centroid an integer vector, the stencil centroid.
#' @param pad an integer, the padding around the stencil centroid.
#' @param obj an imager::cimg or an imager::pixset
#'
#' @returns a data.frame
square_stencil <- function(centroid, pad, obj) {

  # get image dimensions
  if (is.cimg(obj) | is.pixset(obj)) {
    dims <- as.integer(dim(obj)[1L:2L])
  } else {
    stop("Invalid 'obj' argument.")
  }
  # check bounds
  bounds <- c(
    ifelse(centroid - pad < 1L, -centroid + 1L, -pad),
    ifelse(pad + centroid > dims, dims - centroid, pad)
  )
  # define stencil
  stencil <- expand.grid(
    dx = seq(bounds[1L], bounds[3L]),
    dy = seq(bounds[2L], bounds[4L])
  )
  return(stencil)
}


#' Centimeter to pixel
#'
#' Get the centimeter to pixel conversion factor from a ruler defined as
#' a pixset.
#'
#' @param ruler an imager::pixset
#' @param quantile a real, the quantile level used for thresholding.
#'
#' @returns a data.frame
cm_to_pixel <- function(ruler, quantile = 5e-3) {

  # get ruler contours size
  ct_size <- ruler %>%
    contours(nlevels = 1L) %>%
    map(~ length(.$x)) %>%
    unlist()
  # estimate the distribution by a Gaussian
  gm <- densityMclust(log10(ct_size), G = 1L, plot = FALSE)
  # threshold to discriminate square contours
  thresholds <- qnorm(
    c(quantile, 1L - quantile),
    gm$parameters$mean[1L],
    sqrt(gm$parameters$variance$sigmasq[1L])
  )
  # get the conversion factor cm to pixel
  conv_factor <- ct_size %>%
    keep(~ (. %inr% 10L ** thresholds)) %>%
    {. / 4L} %>%
    mean()
  return(conv_factor)
}


#' Ellipsoid hull
#'
#' Fit an ellipsoid hull for a pixset with optional morphological opening.
#'
#' @param pixset an imager::pixset
#' @param size an integer, the morphological opening factor.
#'
#' @returns a data.frame
ehull <- function(pixset, size = 6L) {
  out <- pixset %>%
    mopening_square(size) %>%
    as.pixset() %>%
    where() %>%
    as.matrix() %>%
    ellipsoidhull()
  return(out)
}

#mean value of white pixels
px_mean <- function(pixset) {
  out <- which(pixset == 1, arr.ind = TRUE)[,1:2] %>% colMeans
  return(out)
}

eu_dist <- function(v1, v2) {
  sqrt(sum((v1 - v2)^2))
}


#' Edge collision
#'
#' Check if any white pixel of a binary pixset is touching the edge of the image
#'
edge_touching <- function(pixset){
  image <- pixset[,,1,1] #extract values only
  out <- any(c(image[1, ] , image[nrow(image), ] , image[, 1] , image[, ncol(image)]))
  return <- out
}

#' Check white pixels
#'
#' Get white pixels coordinates of a pixset
#'
white_pix <- function(pixset){
  return(pixset > 0.99)
}
#' Get white pixels coordinates of a pixset for boolean files
#'
white_pix_bool <- function(pixset){
  return(which(pixset == T, arr.ind = TRUE))
}

#' Get crop coordinates
#'
#' Get coordinates of the upper left pixel of the bounding box used to crop a label
#'
get_0coords <- function(pixset) {
  white_pixels <- which(pixset==T, arr.ind = TRUE)
  point0 <- c(white_pixels[which.min(white_pixels[, 1]), 1],
              white_pixels[which.min(white_pixels[, 2]), 2])
  return(point0)
}


#' Remove small spots
#'
#' Remove small spots under a given size threshold on a binary image
#'
remove_small_spots <- function(pixset, size){
  seg_small <- split_connected(pixset)
  if (length(seg_small) == 1){
    return(list(pixset))
  }
  else {
    n_w_pix <- sapply(seg_small, function(x) nrow(white_pix_bool(x)))
    return(seg_small[ which(n_w_pix >= size) ])
  }
}


#' Centroid
#'
#' Calculate the centroid of a pixset = mean of the coordinates
#'
centroid <- function(pixset){
  white_pixels <- white_pix_bool(pixset)
  centroid_x <- mean(white_pixels[, 1])  #mean of x coordinates
  centroid_y <- mean(white_pixels[, 2])  #mean of y coordinates
  centroid_final <- c(centroid_x, centroid_y)
  return(centroid_final)
}

#' Convert pixset to binary image
#'
#' for jpg export, also it's the negative
#'
pixset_to_image <- function(pixset) {
  img <- as.cimg(pixset)  # Convert pixset to cimg object
  neg_img <- 1-img # Negative image for Momocs
  return(neg_img)
}

#' Create directory
#'
#' Create directory if doesn't already exists
#'
create_dir <- function(name, msg=T){
  if (!dir.exists(name)) {
    dir.create(name)
    if(msg) message("Created directory:", name)
  } else {
    if(msg) message("Directory already exists:", name)
  }
}

#' Center Contour
#'
#' Center baserd on outline centroid
#'
centerContour <- function(outline) {
  centroid <- colMeans(outline)  # Calculate centroid
  outline <- sweep(outline, 2, centroid)  # Subtract centroid from each point
  return(outline)
}

#' Major axis angle
#'
#' Find major axis angle using PCA on outline
#'
major_axis_angle <- function(outline) {
  pca <- prcomp(outline)   #PCA to find major axis
  angle <- atan2(pca$rotation[2,1], pca$rotation[1,1])  #calculate rotation angle
  return(angle)
}


#' Align major axis
#'
#' Do rotation based on major PCA axis
#'
align_major_axis <- function(outline) {
  angle <- major_axis_angle(outline)
  rot_mat <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), ncol = 2)
  aligned_outline <- t(rot_mat %*% t(outline))  # Apply rotation
  return(aligned_outline)
}


#' Mirror x or y
#'
#' Mirror set of points horizontally or vertically
#'
mirror_x=function(points){ #symétrie axiale X d'un objet [i,dim]
  meanX=mean(points[,1]) #inversion autour de la valeur abscisse moyenne
  points[,1]= -points[,1]+2*meanX
  return(points)
}
mirror_y=function(points){ #symétrie axiale Y d'un objet [i,dim]
  meanY=mean(points[,2]) #inversion autour de la valeur abscisse moyenne
  points[,2]= -points[,2]+2*meanY
  return(points)
}

#' Vector an  gle
#'
#' Angle of vector from point1 to point2
#'
v_angle <- function(point1, point2){
  dx <- point2[1] - point1[1]
  dy <- point2[2] - point1[2]
  theta <- atan2(dy, dx)
  return(theta)
}


#' Identify sequences of TRUE values in a cyclic boolean vector
#'
find_true_sequences <- function(boolean_vector) {
  true_sequences <- list()
  in_sequence <- FALSE
  start_idx <- NULL
  for (i in seq_along(boolean_vector)) {
    if (boolean_vector[i]) {
      if (!in_sequence) {  #start of a new sequence
        start_idx <- i
        in_sequence <- TRUE
      }
      if (i == length(boolean_vector)) {  #end of sequence at the end of the vector
        true_sequences[[length(true_sequences) + 1]] <- start_idx:i
      }
    } else if (in_sequence) { #end of current sequence
      true_sequences[[length(true_sequences) + 1]] <- start_idx:(i - 1)
      in_sequence <- FALSE
    }
  }
  return(true_sequences)
}


#' Function to calculate cyclic distance between two indices
#'
cyclic_distance <- function(index1, index2, length) {
  index1 <- index1 %% length
  index2 <- index2 %% length
  dist1 <- abs(index1 - index2)
  dist2 <- length - dist1
  return(min(dist1, dist2))
}

#' find the closest sequences in both directions
#'
find_closest_sequences <- function(rear_index, true_sequences, vector_length) {
  closest_forward <- NULL
  closest_backward <- NULL
  min_forward_dist <- Inf
  min_backward_dist <- Inf
  for (sequence in true_sequences) { #calculate forward distance to the start of the sequence
    forward_dist <- cyclic_distance(rear_index, sequence[1], vector_length)
    if (sequence[1] >= rear_index && forward_dist < min_forward_dist) {
      closest_forward <- sequence
      min_forward_dist <- forward_dist
    }
    backward_dist <- cyclic_distance(rear_index, sequence[length(sequence)], vector_length) #calculate backward distance to the end of the sequence
    if (sequence[length(sequence)] <= rear_index && backward_dist < min_backward_dist) {
      closest_backward <- sequence
      min_backward_dist <- backward_dist
    }
    #TODO MOCHE, REVOIR AVEC MODULO%%
    if (is.null(closest_backward)) {
      closest_backward <- true_sequences[[length(true_sequences)]]
    }
    #TODO end MOCHE
  }
  return(list(closest_forward = closest_forward, closest_backward = closest_backward))
}


#' Imager::contours but with simplified output format
#'
simple_imager_contour <- function(pixset){
  contour <- imager::contours(pixset)[[1]]
  s_contour <- cbind(contour[2]$x,contour[3]$y)
  colnames(s_contour) <- c("x", "y")
  return(s_contour)
}

#' Angle between two vectors
#'
normV=function(u){return(sqrt(scalar(u,u)))}
vv_angle <- function(u,v){
  return(acos( round(scalar(u,v)/(normV(u)*normV(v)),14) ))
  } #rad

#' Scalar product
#'
scalar <- function(u, v){
  if(length(u)==3){
    return(u[1]*v[1]+u[2]*v[2]+u[3]*v[3])
  }else if(length(u)==2){
    return(u[1]*v[1]+u[2]*v[2])
  }else{print('Incorrect vector size, must be 2 or 3')}
}

#' Contour angle for each point given a search window = index of the next value to check
#' @export
#'
contourAngles <- function(coords, search_w = 5){
  angcont=c()
  lcont=nrow(coords)
  for(i in 0:(lcont-1)){
    u=as.vector(unlist(coords[(i+search_w)%%lcont+1,]-coords[i%%lcont+1,]))
    v=as.vector(unlist(coords[(i-search_w)%%lcont+1,]-coords[i%%lcont+1,]))
    var=(-vv_angle(u,v)+pi)/pi
    angcont=c(angcont, var)
  }
  return(angcont)
}

#' Euclidian distance between two points
#'
#'TODO : delete > only use one eu_dist function
eu_dist2 <- function(point1, point2) {
  if (!all(c("x", "y") %in% names(point1)) || !all(c("x", "y") %in% names(point2))) { #ensure both points have "x" and "y" coordinates
    stop("Both points must have 'x' and 'y' coordinates.")
  }
  distance <- sqrt((point2["x"] - point1["x"])^2 + (point2["y"] - point1["y"])^2)
  return(distance)
}


#' Body length from body points
#'
bodyLength <- function(body_points){
  long_angle <- major_axis_angle(body_points) #elongation axis angle of the body
  proj_body <- body_points[,1]*cos(long_angle) + body_points[,2]*sin(long_angle) #coordinates on this axis
  return(max(proj_body)-min(proj_body)) #body length
}


#' Fast segmentation with size filter
#'
#' Given a binary image as pixset, returns list of coordinates of every label with pixels in the min-max range
#' @param img_bin a pixset
#' @param min minimum number of pixels to keep label
#' @param max maximum number of pixels to keep label
#' @export
#'
gFastSeg <- function(img_bin, min=300, max=1500){
  # 1 - global size threshold (keep bodies)
  lab <- label(img_bin) #connected components detection on binary image
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



#'Load and binarize image
#'
#'Load image from given path, negative filter, illumination correction filter, binarization based on threshold
#'
#' @export
binaryLoad <- function(img_path, threshold){
  img_pix <- load.image(img_path) #load image with imager
  img_gr <- grayscale(img_pix) #in b/w
  img_gr_inv <- img_gr %>%      #negative and illumination correction
    correct_illumination() %>%
    invert_grayscale()
  grayscale_values <- as.vector(img_gr_inv)
  img_bin <- img_gr_inv %>% threshold(threshold) #binarization based on threshold
  return(img_bin)
}

#' Accessory for redRulerScaler
#'
nZeros <- function(text, zeros=1){
  if(substr(text, zeros, zeros) == 0){
    return( nZeros(text = text,
                   zeros = zeros+1) )
  } else { return (zeros-1) }
}
#'Get scale using detection of red millimeter paper
#'
#' Detect ruler as biggest red label,
#'
#' @export
redRulerScale <- function(img, red_thresh=0.05, confidence_interval=0.95, viz=F, msg=T){
  # Get ruler as the biggest set of red pixels from the image
  if(msg) message("--- scaling - Extracting ruler...")
  IR <- correct_illumination(R(img)) #RGB components with linear correction of illumination
  GR <- correct_illumination(G(img))
  BR <- correct_illumination(B(img))
  red_mask <-  (IR>GR+red_thresh & IR>BR+red_thresh ) #red only component of the image (with few red and blue compared to red)
  lab0 <- label(red_mask)
  lab <- lab0[lab0>0] #filter segmentation of background
  lab_tab <- table(lab) #number of pixels of each label
  largest_label <- as.numeric(names(which.max(lab_tab)))
  ruler_img <- as.cimg(lab0 == largest_label) #ruler's binary image

  # Get indivual tiles from the ruler
  if(msg) message("--- scaling - Segmenting ruler's tiles...")
  ruler_crop <- crop.bbox(ruler_img, ruler_img==T) #cropped ruler
  ruler_neg <- ruler_crop %>% invert_grayscale #negative
  seg_ruler <- gFastSeg(ruler_neg, min=1, max=100000000) #negative ruler segmentation
  l_seg_r <- lapply(seg_ruler, nrow) %>% unlist #label's size
  median_threshold <- which(l_seg_r < 1.5 * median(l_seg_r) & l_seg_r > 0.5 * median(l_seg_r))
  seg_tiles <- seg_ruler[median_threshold] #removes label too big or too small
  if(length(seg_tiles) == 0) stop("Error : Ruler tiles not found")

  # Distance between tiles centroids
  if(msg) message("--- scaling - Computing tile's centroid distances...")
  mean_px_val <- lapply(seg_tiles, colMeans) #vector of xy mean position of the ruler's tiles
  tile_pos <- as.data.frame(do.call(rbind, mean_px_val)) #dataframe with xy mean position of the ruler's tiles
  n_tiles <- nrow(tile_pos) #number of tiles
  pb <- progress_bar$new(total = n_tiles-1)
  distances <- numeric(n_tiles * (n_tiles - 1) / 2) #number of possible combinations between points
  k <- 1
  for(i in 1:(n_tiles - 1)) {
    pb$tick()
    for(j in (i+1):n_tiles) {
      distances[k] <- eu_dist(tile_pos[i,], tile_pos[j,])
      k <- k+1
    }
  }

  # Get mean distance (pixels) between two adjacent tiles (mm)
  if(msg) message("--- scaling - Computing scale...")
  d_hist <- hist(distances, breaks = 500, plot=F) #inbetween tiles distances histogram
  counts_spiked <- spikePlateau(d_hist$counts) #add 'spikes' because pracma::findpeaks does not detect flat peaks, which can frequently happen in histograms
  h_peaks <- pracma::findpeaks(counts_spiked, threshold = 0.2 * max(counts_spiked)) #histogram peaks
  h_thresh <- d_hist$mids[h_peaks[1:2, 2]] %>% mean #mean value between the first two peaks, used as threshold to select 1st peaks surrounding values
  unit_dist <- distances[distances < h_thresh] / 1000
  unit_dist_mean <- unit_dist %>% mean #mean pixel amount for 1 micron
  se <- (unit_dist %>% sd)/(unit_dist %>% length %>% sqrt) #mean's Standard Error
  e_margin <- se*qnorm(1-(1-confidence_interval)/2) #Error margin for a given confidence interval

  #visualization to check correct execution
  if(viz){
    par(mfrow=c(1,2))
    plot(tile_pos, pch=16, cex=0.5, asp=1) #tiles centroids
    plot(d_hist, freq = T, xlim = c(0, max(distances)/5)) #distance histogram
    points(d_hist$mids[h_peaks[, 2]], h_peaks[, 1], col = "red", pch = 19, cex = 1.5) #plotted peaks
    par(mfrow=c(1,1))
  }

  #console print
  if(msg){
    exp_main <- unit_dist_mean %>% abs %>% log10 %>% floor
    e_same_exp <- (e_margin * 10^(-exp_main+1) ) %>% as.character #e_same_exp
    e_subp <- gsub("\\.", "", e_same_exp)
    n_zeros <- nZeros(e_subp)
    e_margin_format <- paste0("0.", substring(e_subp, 2, n_zeros+1), "e-0", abs(exp_main))
    cat('1 µm =', format(unit_dist_mean , scientific=T, digits=3),
        "\u00B1", e_margin_format,'px \n')
  }
  return(c("scale" = unit_dist_mean, "error" = e_margin))
}


#' Plot gerris position and id with 1cm scale on original image
#'
#' name : file name
#' base_img : base image
#' scale : mean pixel amount for 1 micron
#' x,y : ordered position of individuals
#'
#' @export
gDetectionPlot <- function(base_img, scale, x, y){
  par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
  plot(base_img)
  #scale
  x_start <- 100
  y_start <- ncol(base_img)-100
  segments(x_start, y_start, x_start+scale*10000, y_start, col = "brown4", lwd = 4)
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
  pb <- progress_bar$new(total=length(centroids))
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

#' GMM 'V' using Mclust
#'
#' returns the GMM model of a distribution as a numeric vector img_val and a dataframe of the parameters with possibility
#' of removing the one with the highest mean. n integer vector of components to test
#' if n length > 1 : the number of components will be decided using BIC
#'
#' @export
GMM <- function(img_val, n=1:9, keep_max=T, viz=F, n_noise_sd=0){
  if(is.cimg(img_val)){
    img_val <- img_val %>% as.numeric #conversion to numeric if the input is still an image
  }
  if(n_noise_sd != 0){
    img_val <- img_val + rnorm(length(img_val), mean=0, sd=n_noise_sd) #normal noise to cope with image discrete gray values
  }
  mod <- Mclust(img_val, modelNames="V", G=n, verbose=F) #fit model on cropped removed_body noised data
  df_mod <- data.frame(mod$parameters$mean,
                       sqrt(mod$parameters$variance$sigmasq), #relevant parameters as dataframe
                       mod$parameters$pro)
  names(df_mod) <- c("mu","sd","w") #set names
  if(!keep_max){
    df_mod <- df_mod[-which(df_mod[,"mu"]==max(df_mod[,"mu"])),] #remove rightmost component
    df_mod[,"w"] <- df_mod[,"w"]*1/sum(df_mod[,"w"]) #consequently adjust weights
  }
  if(viz){
    hist(img_val, breaks=500, freq=F, border=NA) #sum classes = 1
    x_vals <- seq(0, 1, by=10e-4)
    y_vals <- apply(df_mod, 1, function(x){
      dnorm(x_vals, mean = x["mu"], sd = x["sd"]) * x["w"]
    })
    ng <- ncol(y_vals)
    #print(str(y_vals))
    for(i in 1:ng){
      points(y_vals[,i]~x_vals, type='l', col=(1+i), lwd=2)
    }
    Y <- apply(y_vals, 1, sum)
    points(Y~x_vals, type='l', col=1, lty="dotted", lwd=2)
    legend("topright", legend=seq(1,ng), col=seq(1,ng)+1,lty=1, lwd=3)
  }
  return(list("df"=df_mod, "model"=mod)) # [[1]]: parameters dataframe [[2]]:Mclust global output
}

#' Quantile off GMM model dataframe gmm_model_df from GMM function
#'
#' @export
GMMquant <- function(gmm_model_df, q, viz=F){
  f <- function(x){
    apply(gmm_model_df, 1, function(df){
      pnorm(x, mean=df["mu"], sd=df["sd"]) * df["w"]
    }) %>% sum - q
  }
  x_quant <- uniroot(f, c(0,10))[1] %>% unlist #TODO INTERVAL???
  if(viz){ #visualization option
    x_vals <- seq(0, 1, by=10e-4)
    f2 <- function(x){
      apply(gmm_model_df, 1, function(df){
        dnorm(x, mean=df["mu"], sd=df["sd"]) * df["w"]
      }) %>% sum
    }
    y_vals <- sapply(x_vals , f2)
    lim <- c(0,1)
    if(q<0.99 & q>0.01){   #adjust scale for non extreme quantile values
      lim <- c(uniroot(function(x) f(x)+q-0.01, c(0,1))[1] %>% unlist,
               uniroot(function(x) f(x)+q-0.99, c(0,1))[1] %>% unlist)
    }
    plot(y_vals~x_vals, type="l", xlim=lim)
    abline(v=x_quant, col=2)
    points(x=x_quant, y=f2(x_quant), col=2, cex=1.3, pch=16)
  }
  return( x_quant %>% unlist )
}

#' Roots of the derivative of a GMM
#'
#' Give roots of y for the derivative of given GMM parameters as dataframe created with GMM function
#'
#' @export
GMMdRoots=function(GMM_df, y, viz=F, msg=T){
  if(!is.numeric(y)) stop("y must be a numeric value") #error handling
  if(!is.data.frame(GMM_df)) stop("GMM_df must be a dataframe")
  if(any(colnames(GMM_df) != c("mu","sd","w"))) stop("GMM_df must be a dataframe created with GMM")
  #FX
  roots <- uniroot.all(function(x) f_prime(x, GMM_df)-y , interval = c(0, 1), tol=10e-8) #solve with rootSolve package
  if(viz){
    roots_plot <- c(uniroot.all(function(x) f_prime(x, GMM_df)-5 , interval = c(0, 1), tol=10e-8), #auto adjust plot limits
                    uniroot.all(function(x) f_prime(x, GMM_df)+5 , interval = c(0, 1), tol=10e-8))
    plot(function(x) f_prime(x,GMM_df), xlim=c(min(roots_plot)-0.05, max(roots_plot)+0.05), #plot fx
         main = paste0("GMM derivative, y = ",y), ylab = "y", xlab = "x")
    points(x=roots, y=rep(y, length(roots)), col=2, pch=16) #plot found points
    abline(h=y,col=2)
  }
  if(length(roots)==0 && msg) message("Failed to find any solutions")
  return(roots)
}
#' Accesory function for GMMdRoots : creates derivative function of GMM
#'
f_prime <- function(x, GMM_df){
  sapply(x, function(xi){ # !!! VECTORIZATION !!!
    apply(GMM_df, 1, function(df){
      mu <- df["mu"]
      sd <- df["sd"]
      w <- df["w"]
      -((xi - mu) / sd^2) * dnorm(xi, mean = mu, sd = sd) * w
    }) %>% sum
  })
}

#' Binarize gerris using GMM derivative
#'
#' Binarize gerris by setting a fixed threshold per input image using GMM derivative
#'
#' @export
gGMMThresh <- function(gerris_crops, nG=2, noise_n=0.005, y_root=-25, msg=T){
  if(is.list(gerris_crops)){ #Vectorization
    if(msg) message("GMM thresholding")
    pb <- progress_bar$new(total=length(gerris_crops))
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
    pb <- progress_bar$new(total = length(l_img_bin))
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
  img_res <- as.cimg(labs %in% labs_val_rsp, dim=c(dimX, dimY, 1, 1))
  if(as_coords){ img_res <- img_res %>% imgAsCoords }

  return(img_res)
}

#' Convert xy coordinates dataset as square c.img with padding
#'
#' @export
coordsAsImg <- function(coords, padding=0, return_offset=F){
  coords <- coords %>% as.matrix #avoid issues related to dataframe format, considered as a type of list
  xy_range <- apply(coords, 2, range) #max and min values of coords
  max_range <- apply(xy_range, 2, diff) %>% max + 2*padding #longest dimensional span (x or y) + padding
  xy_mid <- apply(xy_range, 2, mean) #middle of coords to center the new image
  fxy_range <- sapply(xy_mid, function(x) x+c(-max_range, max_range)/2) #x and y range of the new image
  img_matrix <- matrix(0, nrow = max_range+1, ncol = max_range+1) #empty matrix as base for new img
  min_range <- apply(fxy_range, 2, function(x) x %>% min %>% round)
  coord_offset <- matrix(min_range, nrow=nrow(coords), ncol=2, byrow=T)
  xy_coords <- coords - coord_offset+1
  img_matrix[xy_coords] <- 1
  img <- as.cimg(img_matrix)
  if(!return_offset){
    return(img)
  } else {    #if return_offset=True, also returns the offset between old and new coordinates
    return(list(img=img, coord_offset=coord_offset-1))
  }
}

#' Image '1' values as xy coordinates
#'
#' @export
imgAsCoords <- function(img){
  if(is.list(img)){
    return(lapply(img, imgAsCoords))
  }
  if(all(is.na(img))) return(NA)
  return( which(img == 1, arr.ind = TRUE)[,1:2] )
}

#' Convert xy coordinates as linear index in given c.img
#'
#' @export
coordsToLinear <- function(coords, img){
  return((coords[,2]-1) * width(img) + coords[,1])
}

#' Clean body shape as xy coordinates using morphological closing following by opening
#'
#' Square closing then closing of size kernel_size on body coordinates body_lab_points
#' Returns either coordinates as_coord, or an image c.img
#'
#' @export
cleanBodyShape <- function(body_lab_points, kernel_size = 3, as_coords = T){
  if(!is.list(body_lab_points)){
    body_lab_points <- list(body_lab_points) #convert to list if theres a single element for correct use of lapply
  }
  body_img <- lapply(body_lab_points, coordsAsImg, return_offset = T, padding = kernel_size) #convert to img
  closed_body_img <- lapply(body_img, function(x) mclosing_square(x$img, size = kernel_size)) #morphological closing
  opened_body_img <- lapply(closed_body_img, mopening_square, size = kernel_size) #morphological opening
  res <- opened_body_img
  if(as_coords){
    res <- mapply(function(x, y, z){
      coords0 <- imgAsCoords(x) - (kernel_size + 2)
      offset <- z$coord_offset[1,] %>% round + kernel_size+2
      offset_matrix <- matrix(rep(offset, nrow(coords0)), ncol=2, byrow=T)
      coords <- coords0 + offset_matrix
      return(coords)
    }, opened_body_img, body_lab_points, body_img)
  }
  return(res)
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
    EBImage::makeBrush(size, shape='disc')
  })
  dil_body <- mapply(function(crop, kern){
    EBImage::dilate(crop, kern=kern)
  }, body_img, kernels) #dilation around body to ensure intersecting legs and not the body
  return(dil_body)
}

#' Simple contour coords of binary image
#'
#' @export
simpleCont <- function(img_bin){
  if(is.list(img_bin)){
    return(lapply(img_bin, simpleCont))
  }
  cont <- imager::contours(img_bin, nlevels=1) #contour of the dilation as the intersection line
  #note : this is suboptimal as imager::contours uses subpixel level alogrithm which is useless in this case
  df <- data.frame("dim1" = round(cont[[1]]$x), "dim2" = round(cont[[1]]$y)) #round
  cont_coords <-  df[!duplicated(df), ] #removes duplicates
  rownames(cont_coords) <- seq_along(cont_coords[,1]) #coherent rownames
  return(cont_coords)
}

#' With 2 lists of discrete 2D points : which points from set_a are in set_b
#'
#' @export
overlapPoints <- function(set_a, set_b, coords=T, index=T){
  set_a <- set_a %>% as.matrix
  set_b <- set_b %>% as.matrix
  ptInPts <- function(pt, coords){
    which(pt[1] == coords[,1] & pt[2] == coords[,2]) #2D vector both equal
  }
  if(is.list(set_a) & is.list(set_b)){ #vectorization
    res <- mapply(overlapPoints, set_a, set_b, SIMPLIFY = FALSE)
    return(res)
  }
  if(length(set_a)==2) { #formating required if there's a single point in either set
    set_a <- matrix(set_a, ncol = 2, byrow = TRUE)
  }
  if(length(set_b)==2) {
    set_b <- matrix(set_b, ncol = 2, byrow = TRUE)
  }
  interID <- apply(set_b, 1, ptInPts, coords=set_a) %>% unlist
  inter <- 1:nrow(set_a) %in% interID
  if(index && coords){
    return(list("coords"=set_a[inter,], "index"=inter))
  } else if(coords){
    return(set_a[inter,])
  } else if(index){
    return(inter)
  }
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
  pca <- prcomp(body_coords) #principal direction given by body PCA
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
    pb <- progress_bar$new(total=length(gerris))
    res <- mapply(function(a,b,c,d){
      pb$tick()
      return(gLegSeg(a,b,c,d))
    } , gerris, dilated_body, intersection_coords, insertions, SIMPLIFY = F)
    return(res)
  }
    # find limbs
  limbs <- ((dilated_body - gerris) < 0) %>% as.cimg #full individual minus dilated body
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



#' Norm of vector
#'
#' @export
normV <- function(u){ #2D vector norm
  sqrt(scalar(u,u))
}
#' Agle of vector
#'
#' @export
angleV <- function(u, v){
  round(scalar(u,v)/(normV(u)*normV(v)),14) %>% acos
}

#' Contour angle peaks, convex or concave
#'
#' Modified from fxLPJ.R / Peaks of angle contour for a givens search window, smoothed using splines
#'
#' @export
contourTurns <- function(coords, search_w=5, splines_df=30, angle_thresh=0.15, viz=0){
  if(is.list(coords) & !is.data.frame(coords)){ #Vectorization
    message("Contour angle peaks")
    pb <- progress_bar$new(total=length(coords))
    res <- lapply(coords, function(x){
      pb$tick()
      contourTurns(x, search_w, splines_df, angle_thresh)
      })
    return(res)
  }
  angcont <- c()
  lcont <- nrow(coords) #0:(lcont-1)
  for(i in 0:(lcont-1)){ #adjust modulo functions for R indexation starting with 1 making it unpractical
    u <- as.vector(unlist(coords[(i+search_w)%%lcont+1,]-coords[i%%lcont+1,]))
    v <- as.vector(unlist(coords[(i-search_w)%%lcont+1,]-coords[i%%lcont+1,]))
    var <- (-angleV(u,v)+pi)/pi #scale between 0 and 1
    angcont <- c(angcont, var)
  }
  index <- seq(1, lcont)
  m_index <- c(index-max(index), index, index+max(index)) #on met bout à bout 3 fois la séquence
  m_angcont <- rep(angcont, 3)
  #splines and max
  splines_df_3 <- splines_df*3
  if(splines_df_3 >= (length(unique(m_index))-1) ) {
    stop("Contour smaller than splines_df")
  }
  m_smoothed <- smooth.spline(m_index, m_angcont, df=splines_df_3) #splines
  m_pred <- predict(m_smoothed) #model prediction for initial values
  md_pred <- diff(m_pred$y)/diff(m_pred$x) #derivative
  cross0 <- which(diff(sign(md_pred))!=0)+1 #f'=0 -> plateau
  mid_sel <- c(max(index)+1, 2*max(index))#selection of middle of extended range
  ind_cross <- cross0[cross0 <= mid_sel[2] & cross0 >= mid_sel[1]] - lcont #f'=0 as index
  d_pred <- md_pred[mid_sel[1]:mid_sel[2]] #base range
  #get points (convex concave)
  first_convex <- if(d_pred[ind_cross[1]]<0){1} else{0} #if after the first point f' goes up or down
  convex <- ind_cross[seq(first_convex, length(ind_cross), by = 2)]
  concave <- ind_cross[seq(first_convex+1, length(ind_cross), by = 2)]
  yconc <- angcont[concave]; yconv <- angcont[convex]
  concave <- concave[yconc > angle_thresh] #delete inflection points found on small angles
  convex <- convex[yconv > angle_thresh]
  yconc <- yconc[yconc > angle_thresh]
  yconv <- yconv[yconv > angle_thresh]
  #viz
  colvcont <- sapply(angcont, function(x){ gray(1-x) }) #color values for plot
  if(viz==1){ #viz 1
    par(mfrow=c(1,2))
    plot(angcont, main='Splines and peak fitting', type='l', lwd=2)
    points(m_pred, type='l', col='green3', lwd=2)
    points(yconc~concave, col='blue', pch=16)
    points(yconv~convex, col='red', pch=16)
    plot(coords, col=colvcont, main='Contour angle variation', pch=16, as=1)
    points(coords[concave,], col='blue', cex=1, pch=13)
    points(coords[convex,], col='red', cex=1, pch=13)
    par(mfrow=c(1,1))
  }else if(viz==2){ #viz 2
    plot(coords, col=colvcont, main='Contour angle variation', pch=16, as=1)
    points(coords[concave,], col='blue', cex=2)
    points(coords[convex,], col='red', cex=2)
  }
  return(convex)
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
    knee_pts <- do.call(rbind,filt_cont_pts)
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
      plot(x=dist_peaks_str, y=rep(1.1,n_p),pch=16,ylim=c(0.2,1.8),
           xlim=c(0,sum(dist)/2), main="Inflexion Profile")
      points(x=c(0,sum(dist)),y=c(1,1),col="blue",pch=16)
      points(x=dist_peaks_rev, y=rep(0.9,n_p),pch=16,ylim=c(0.5,1.5),xlim=c(0,sum(dist)))
      abline(v=sum(dist)/2*thresh, col="blue")
      legend("bottomright", legend =c(paste("thresh",thresh),"mid leg"), col = c("blue",1), lty = 1, bty = "n")
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
  outs <- Out(cent_cont) #outline Momocs object
  scale_cont <- outs %>% coo_scale #scale
  efa <- scale_cont %>% efourier(nb.h=nb_h, norm=T) #elliptic fourier analysis
  efa_pca <- prcomp(efa$coe) #principal components analysis
  if(viz){ #visualization
    par(mfrow=c(1,1))
    scale_cont %>% stack #outlines superposition
    scale_cont %>% panel #all outlines
    scale_cont %>% calibrate_harmonicpower_efourier(nb.h=nb_h) #elliptic fourier calibration
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


#' Plot label output with colors
#'
#' @param lab_img c.img
#' @export
label2RGB <- function(lab_img) {
  num_labels <- max(lab_img)
  colors <- data.frame(col2rgb(rainbow(num_labels))) / 255
  dims <- dim(lab_img)
  r <- imfill(dims[1], dims[2])
  g <- imfill(dims[1], dims[2])
  b <- imfill(dims[1], dims[2])
  for (i in 1:num_labels) {
    mask <- lab_img == i
    r[mask] <- colors[1, i]
    g[mask] <- colors[2, i]
    b[mask] <- colors[3, i]
  }
  color_img <- imappend(list(r, g, b), "c")
  return(color_img)
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

#' Add pike to plateaus of a numeric vector
#'
#' @param v a numeric vector
#' @param increment increment coefficient "size of pyramid steps"
#' @export
spikePlateau <- function(v, increment = 1e-5){
  len <- length(v)
  same <- v == v[(1:len) %% len+1]
  boolseq <- boolSeqLim(same)
  addval <- rep(0, len)
  for(i in 1:length(boolseq$first)){
    first <- boolseq$first[i]
    ran <- boolseq$last[i] - first+1
    if(ran != 1){
      pyramid <- c(1:ceiling((ran+1)/2), (floor((ran+1)/2):1)[-1]) #test it for yourself with any int 'ran'
      add_i <- pyramid * increment
    } else {
      add_i <- increment #single increment for single value
    }
    addval[first:(first+ran-1)] <- add_i
  }

  return(v+addval)
}


#' The main thing... for now
#'
#' @param img_path image path or list of image paths
#' @param write_df a boolean, to save metrics as csv
#' @param write_plots a boolean, to save image analysis diagnostic plots
#' @param df_output a boolean, to return metrics in R
#' @export
gPipeline <- function(img_path, write_df=T, write_plots=T, df_output=T){
  message("1 - Loading and segmenting image")
  base_img <- load.image(img_path)
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
      png(paste0(iplot_dir,"/",clean_base_path,"_",i,".png"), width=480, height=480)  # specify filename and dimensions
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
