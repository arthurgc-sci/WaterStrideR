#utils related to insect leg analysis pipeline developped here mostly by Arthur Gairin-Calvo but also Laurent Gilquin
#---------------------------------------------------------------------------------------------------

#' Correct illumination
#'
#' Correct a gray-scaled image illumination by fitting a linear model and
#' removing the spatial trend.
#'
#' @importFrom magrittr %>%
#' @importFrom grDevices dev.off gray png
#' @importFrom graphics abline hist legend lines par points symbols text
#' @importFrom stats dnorm lm median pnorm prcomp predict qnorm rnorm sd uniroot na.omit
#' @importFrom utils stack write.table
#' @importFrom imager %inr%
#' @import MASS
#' @importFrom mclust Mclust mclustBIC
#'
#' @param img an imager::cimg
#' @param nsamples an integer, pixel subsampling value.
#'
#' @returns an imager::cimg object
#' @export
correct_illumination <- function(img, nsamples = 1e4L) {
  # convert to grayscale if needed
  if (rev(dim(img))[1L] > 1L) {
    img <- imager::grayscale(img)
  }
  # linear regression trend
  trend <- img %>%
    as.data.frame() %>%
    dplyr::sample_n(nsamples) %>%
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
#' @export
invert_grayscale <- function(img) {
  # convert to grayscale if needed
  if (rev(dim(img))[1L] > 1L) {
    img <- imager::grayscale(img)
  }
  out <- max(img) - img
  return(out)
}

#' Euclidian distance between two points
#'
#' @param v1 numerical vector of length 2
#' @param v2 numerical vector of length 2
#' @keywords internal
eu_dist <- function(v1, v2) {
  sqrt(sum((v1 - v2)^2))
}

#' Create directory
#'
#' Create directory if it doesn't already exists
#'
#' @param name a string. name and path of the directory to create
#' @param msg logical. console message option  
#' @keywords internal
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
#' @param outline matrix or dataframe of points coordinates
#' @keywords internal
centerContour <- function(outline) {
  centroid <- colMeans(outline)  # Calculate centroid
  outline <- sweep(outline, 2, centroid)  # Subtract centroid from each point
  return(outline)
}

#' Major axis angle
#'
#' Find major axis angle using PCA on outline
#'
#' @param outline matrix or dataframe of points coordinates
#' @keywords internal
major_axis_angle <- function(outline) {
  pca <- prcomp(outline)   #PCA to find major axis
  angle <- atan2(pca$rotation[2,1], pca$rotation[1,1])  #calculate rotation angle
  return(angle)
}

#' Identify sequences of TRUE values in a cyclic boolean vector
#'
#' @param boolean_vector a vector of boolean values
#' @export
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

#' Norm of vector
#'
#' @param u a numerical vector of length 2
#' @export
normV <- function(u){ #2D vector norm
  sqrt(scalar(u,u))
}

#' Angle of vector
#'
#' @param u a numerical vector of length 2
#' @param v a numerical vector of length 2
#' @export
angleV <- function(u, v){
  round(scalar(u,v)/(normV(u)*normV(v)),14) %>% acos
}

#' Scalar product
#'
#' @param u a numerical vector of length 2
#' @param v a numerical vector of length 2
#' @keywords internal
scalar <- function(u, v){
  if(length(u)==3){
    return(u[1]*v[1]+u[2]*v[2]+u[3]*v[3])
  }else if(length(u)==2){
    return(u[1]*v[1]+u[2]*v[2])
  }else{print('Incorrect vector size, must be 2 or 3')}
}

#' Contour angle for each point given a search window = index of the next value to check
#' 
#' @param coords dataframe or matrix of contour coordinates
#' @param search_w numerical integer. Index of the next value to check in the contour 
#' @export
contourAngles <- function(coords, search_w = 5){
  angcont=c()
  lcont=nrow(coords)
  for(i in 0:(lcont-1)){
    u=as.vector(unlist(coords[(i+search_w)%%lcont+1,]-coords[i%%lcont+1,]))
    v=as.vector(unlist(coords[(i-search_w)%%lcont+1,]-coords[i%%lcont+1,]))
    var=(-angleV(u,v)+pi)/pi
    angcont=c(angcont, var)
  }
  return(angcont)
}

#' Body length from body points
#'
#' @param body_points a dataframe or matrix of body points coordinates
#' @keywords internal
bodyLength <- function(body_points){
  long_angle <- major_axis_angle(body_points) #elongation axis angle of the body
  proj_body <- body_points[,1]*cos(long_angle) + body_points[,2]*sin(long_angle) #coordinates on this axis
  return(max(proj_body)-min(proj_body)) #body length
}


#'Load and binarize image
#'
#'Load image from given path, negative filter, illumination correction filter, binarization based on threshold
#'
#' @param img_path a string. path to the base image  
#' @param threshold a numerical in range 0:1 for binarization threshold 
#' @export 
binaryLoad <- function(img_path, threshold){
  img_pix <- imager::load.image(img_path) #load image with imager
  img_gr <- imager::grayscale(img_pix) #in b/w
  img_gr_inv <- img_gr %>%      #negative and illumination correction
    correct_illumination() %>%
    invert_grayscale()
  grayscale_values <- as.vector(img_gr_inv)
  img_bin <- img_gr_inv %>% imager::threshold(threshold) #binarization based on threshold
  return(img_bin)
}

#' Convert xy coordinates dataset as square c.img with padding
#'
#' @param coords matrix or dataframe of numerical coordinates
#' @param padding numerical integer of padding to add to the crop
#' @param return_offset logical. to also return offset value of recrop
#' @export
coordsAsImg <- function(coords, padding=0, return_offset=F){
  #safety
  coords <- coords %>% as.matrix #avoid issues related to dataframe format, considered as a type of list
  if(length(coords)==2){  #single points, no point in creating an image
    warning("coordsAsImg not implemented for single point, returning NA")
    return(NA)
  }
  #
  xy_range <- apply(coords, 2, range) #max and min values of coords
  max_range <- apply(xy_range, 2, diff) %>% max + 2*padding #longest dimensional span (x or y) + padding
  xy_mid <- apply(xy_range, 2, mean) #middle of coords to center the new image
  fxy_range <- sapply(xy_mid, function(x) x+c(-max_range, max_range)/2) #x and y range of the new image
  img_matrix <- matrix(0, nrow = max_range+1, ncol = max_range+1) #empty matrix as base for new img
  min_range <- apply(fxy_range, 2, function(x) x %>% min %>% round)
  coord_offset <- matrix(min_range, nrow=nrow(coords), ncol=2, byrow=T)
  xy_coords <- coords - coord_offset+1
  img_matrix[xy_coords] <- 1
  img <- imager::as.cimg(img_matrix)
  if(!return_offset){
    return(img)
  } else {    #if return_offset=True, also returns the offset between old and new coordinates
    return(list(img=img, coord_offset=coord_offset-1))
  }
}

#' Image '1' values as xy coordinates
#'
#' @param img a pixet or c.img image
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
#' @param img a c.img
#' @param coords matrix or dataframe of numerical coordinates
#' @export
coordsToLinear <- function(coords, img){
  return((coords[,2]-1) * imager::width(img) + coords[,1])
}


#' Clean body shape as xy coordinates using morphological closing following by opening
#'
#' Square closing then closing of size kernel_size on body coordinates body_lab_points
#' Returns either coordinates as_coord, or an image c.img
#'
#' @param body_lab_points body points as numerical matrix or dataframe
#' @param kernel_size numerical integer. closing and opening kernel size in pixels
#' @param as_coords logical. to return
#' @export
cleanBodyShape <- function(body_lab_points, kernel_size = 3, as_coords = T){
  if(!is.list(body_lab_points)){
    body_lab_points <- list(body_lab_points) #convert to list if theres a single element for correct use of lapply
  }
  body_img <- lapply(body_lab_points, coordsAsImg, return_offset = T, padding = kernel_size) #convert to img
  closed_body_img <- lapply(body_img, function(x) imager::mclosing_square(x$img, size = kernel_size)) #morphological closing
  opened_body_img <- lapply(closed_body_img, imager::mopening_square, size = kernel_size) #morphological opening
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

#' Simple contour coords of binary image
#'
#' @param img_bin a binary c.img
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
#' @param set_a a numerical matrix of 2D points
#' @param set_b a numerical matrix of 2D points
#' @param coords logical. to return numerical points coordinates
#' @param index logical. to return point's index
#' @export
overlapPoints <- function(set_a, set_b, coords=T, index=T){
  set_a <- set_a %>% as.matrix
  set_b <- set_b %>% as.matrix
  ptInPts <- function(pt, crds){
    which(pt[1] == crds[,1] & pt[2] == crds[,2]) #2D vector both equal
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
  interID <- apply(set_b, 1, ptInPts, crds=set_a) %>% unlist
  inter <- 1:nrow(set_a) %in% interID
  if(index && coords){
    return(list("coords"=set_a[inter,], "index"=inter))
  } else if(coords){
    return(set_a[inter,])
  } else if(index){
    return(inter)
  }
}

#' Contour angle peaks, convex or concave
#'
#' Modified from fxLPJ.R / Peaks of angle contour for a givens search window, smoothed using splines
#'
#' @param coords numerical dataframe or list of dataframe. Ordered contour points coordinates
#' @param search_w a numerical integer. Distance, in number of points of in the contour
#' to compute angles from
#' @param splines_df a numerical integer. spline's degrees of freedom
#' for contour angle histogram smoothing
#' @param angle_thresh a numerical. Angle value (in radians) under which 
#' inflection points will be dismissed 
#' @param viz logical. visualization option
#' @export
contourTurns <- function(coords, search_w=5, splines_df=30, angle_thresh=0.15, viz=0){
  if(is.list(coords) & !is.data.frame(coords)){ #Vectorization
    message("Contour angle peaks")
    pb <- progress::progress_bar$new(total=length(coords))
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
  m_index <- c(index-max(index), index, index+max(index)) #on met bout à bout 3 fois la sequence
  m_angcont <- rep(angcont, 3)
  #splines and max
  splines_df_3 <- splines_df*3
  if(splines_df_3 >= (length(unique(m_index))-1) ) {
    stop("Contour smaller than splines_df")
  }
  m_smoothed <- stats::smooth.spline(m_index, m_angcont, df=splines_df_3) #splines
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

#' Make Disc Kernel
#'
#' Make Filled Disc Kernel matrix for morphological dilations or erosions
#' 
#' @param size numeric value for disc size, will be rounded to the next odd number*
#' @keywords internal
make_disc_kernel <- function(size) {
  size <- ifelse(size %% 2 == 0, size + 1, size)
  radius <- floor(size / 2) + 0.5
  center <- (size + 1) / 2
  grid <- expand.grid(x = 1:size, y = 1:size)
  kernel <- with(grid, matrix(
    as.numeric(sqrt((x - center)^2 + (y - center)^2) < radius),
    nrow = size, ncol = size
  ))
  return(kernel)
}
