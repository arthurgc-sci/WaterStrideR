#utils related to insect leg analysis pipeline developped here mostly by Arthur Gairin-Calvo but also Laurent Gilquin
#---------------------------------------------------------------------------------------------------

#' Correct illumination
#'
#' Correct a gray-scaled image illumination by fitting a linear model and
#' removing the spatial trend.
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
eu_dist <- function(v1, v2) {
  sqrt(sum((v1 - v2)^2))
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

#' Norm of vector
#'
#' @export
normV <- function(u){ #2D vector norm
  sqrt(scalar(u,u))
}

#' Angle of vector
#'
#' @export
angleV <- function(u, v){
  round(scalar(u,v)/(normV(u)*normV(v)),14) %>% acos
}

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
    var=(-angleV(u,v)+pi)/pi
    angcont=c(angcont, var)
  }
  return(angcont)
}

#' Body length from body points
#'
bodyLength <- function(body_points){
  long_angle <- major_axis_angle(body_points) #elongation axis angle of the body
  proj_body <- body_points[,1]*cos(long_angle) + body_points[,2]*sin(long_angle) #coordinates on this axis
  return(max(proj_body)-min(proj_body)) #body length
}


#'Load and binarize image
#'
#'Load image from given path, negative filter, illumination correction filter, binarization based on threshold
#'
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
  lab0 <- imager::label(red_mask)
  lab <- lab0[lab0>0] #filter segmentation of background
  lab_tab <- table(lab) #number of pixels of each label
  largest_label <- as.numeric(names(which.max(lab_tab)))
  ruler_img <- imager::as.cimg(lab0 == largest_label) #ruler's binary image
  
  # Get indivual tiles from the ruler
  if(msg) message("--- scaling - Segmenting ruler's tiles...")
  ruler_crop <- imager::crop.bbox(ruler_img, ruler_img==T) #cropped ruler
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
  pb <- progress::progress_bar$new(total = n_tiles-1)
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
  mod <- mclust::Mclust(img_val, modelNames="V", G=n, verbose=F) #fit model on cropped removed_body noised data
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
  roots <- rootSolve::uniroot.all(function(x) f_prime(x, GMM_df)-y , interval = c(0, 1), tol=10e-8) #solve with rootSolve package
  if(viz){
    roots_plot <- c(rootSolve::uniroot.all(function(x) f_prime(x, GMM_df)-5 , interval = c(0, 1), tol=10e-8), #auto adjust plot limits
                    rootSolve::uniroot.all(function(x) f_prime(x, GMM_df)+5 , interval = c(0, 1), tol=10e-8))
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
  img <- imager::as.cimg(img_matrix)
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

#' Contour angle peaks, convex or concave
#'
#' Modified from fxLPJ.R / Peaks of angle contour for a givens search window, smoothed using splines
#'
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
#' @param size numeric value for disc size, will be rounded to the next odd number
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