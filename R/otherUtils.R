#unused utils, mostly from Laurent Gilquin but also Arthur Gairin-Calvo
#---------------------------------------------------------------------------------------------------

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

#' Edge collision
#'
#' Check if any white pixel of a binary pixset is touching the edge of the image
#'
edge_touching <- function(pixset){
  image <- pixset[,,1,1] #extract values only
  out <- any(c(image[1, ] , image[nrow(image), ] , image[, 1] , image[, ncol(image)]))
  return <- out
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

#' Convert pixset to binary image
#'
#' for jpg export, also it's the negative
#'
pixset_to_image <- function(pixset) {
  img <- as.cimg(pixset)  # Convert pixset to cimg object
  neg_img <- 1-img # Negative image for Momocs
  return(neg_img)
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

#' Get white pixels coordinates of a pixset for boolean files
#'
white_pix_bool <- function(pixset){
  return(which(pixset == T, arr.ind = TRUE))
}

#' Check white pixels
#'
#' Get white pixels coordinates of a pixset
#'
white_pix <- function(pixset){
  return(pixset > 0.99)
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
v_angle_old <- function(point1, point2){
  dx <- point2[1] - point1[1]
  dy <- point2[2] - point1[2]
  theta <- atan2(dy, dx)
  return(theta)
}

#' Euclidian distance between two points
#'
eu_dist_old <- function(point1, point2) {
  if (!all(c("x", "y") %in% names(point1)) || !all(c("x", "y") %in% names(point2))) { #ensure both points have "x" and "y" coordinates
    stop("Both points must have 'x' and 'y' coordinates.")
  }
  distance <- sqrt((point2["x"] - point1["x"])^2 + (point2["y"] - point1["y"])^2)
  return(distance)
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