#' Load LDA model for package
#'
#' @param model_name a string corresponding to the name of the model
#' @return A fitted LDA or PCA model
loadModel <- function(model_name) {
  model_name2 <- paste0(model_name, ".rds")
  model_path <- system.file("models", model_name2, package = "WaterStrideR")
  if (model_path == ""){
    stop("Model file not found in package's 'models' directory")
  }
  model <- readRDS(model_path) #load model, will stop if it does not exist
  if(!inherits(model, "lda") && !inherits(model, "prcomp")){
    stop("Loaded model must be of class 'lda' or 'prcomp'")
  }
  return(model)
}

#' Body contour and orientation 
#'
#' Orient body contours extracted from image using provided angles then define consistent
#' contour starting point as leftmost point
#'
#' @param img_bin c.img or list of c.img of waterstrider body image(s)
#' @param ori_angle angle, vector of angles or list of angles of the segmented bodies
#' @export
normBody <- function(img_bin, ori_angle){
  if(is.list(img_bin) | length(ori_angle)>1){ #VECTORIZATION
    #safety for batches
    na_ang <- ori_angle %>% is.na #valid angle boolean
    if(sum(na_ang) != 0){ #if there is na angles
      ori_angle_c <- ori_angle[!na_ang] #analysis on valid data only
      img_bin_c <- img_bin[!na_ang]
      if(length(img_bin_c)==0) stop("No non NA angle for any of the images provided")
    } else {
      ori_angle_c <- ori_angle
      img_bin_c <- img_bin
    }
    if(length(ori_angle_c)!=length(img_bin_c)) stop("mismatching lengths for img_bin and ori_angle")
    #vectorized FX
    res <- mapply(normBody, img_bin_c, ori_angle_c)
    return(res)
  }
  #safety for single angle-img pair
  if(!imager::is.cimg(img_bin)) stop("img_bin should be a c.img or list of c.img")
  if(!is.numeric(ori_angle)) stop("ori_angle should be a number o list of numbers")
  #FX
  body_cont <- img_bin %>% simpleCont #cleaned body contour
  body_cont_mat <- body_cont[c(seq_len(nrow(body_cont)), 1), ] # Close contour
  body_cont_mat <- body_cont_mat %>% as.matrix # As matrix, format for Momocs
  cos_a <- cos(ori_angle) # Rotation matrix using detected orientation
  sin_a <- sin(ori_angle)
  ang_mat <- matrix(c(cos_a, sin_a, -sin_a, cos_a), ncol = 2)
  rot_cont <- body_cont_mat %*% ang_mat # Rotate
  start_idx <- which.min(rot_cont[, 1]) # Consistent starting point
  aligned_cont <- rbind(rot_cont[start_idx:nrow(rot_cont), ], rot_cont[1:(start_idx - 1), ])   # Start from leftmost point (minimum X) to normalize contour insertion thus avoid 180° random rotations
  
  return(aligned_cont)
}

#' Momocs EFA from waterstrider body image
#'
#' @param ori_cont contour or list of contours (as 2D numeric dataframe of matrix) of waterstrider body image(s)
#' @param nb_h number of harmonics to compute
#' @param viz logical for plot of alignment, fourier space and PCA
#' @param resamp number of points to resample on each contour
#' @importFrom Momocs coo_scale coo_center coo_interpolate efourier panel calibrate_harmonicpower_efourier PCA plot_PCA
#' @export
scoresEFA <- function(ori_cont, nb_h=7, viz=FALSE, resamp=100){
  if(!is.list(ori_cont)){
    single <- TRUE
    #format check for single
    if(length(dim(ori_cont))==2){
      if(dim(ori_cont)[2]==2){
        if(!is.matrix(ori_cont)){
          ori_cont <- ori_cont %>% as.matrix
        }
        if(!is.numeric(ori_cont)) stop ("ori_cont must be numerical")
        ori_cont <- list(ori_cont) #format
      } else stop("ori_cont must have exactly two columns")
    } else stop("ori_cont must be a contour with 2 dimensions (rows and columns)")
  } else single <- FALSE
  outs <- Momocs::Out(ori_cont) #outline Momocs object
  outs_rs <- outs %>% coo_interpolate(n = resamp)  #resample
  efa <- outs_rs %>% efourier(nb.h=nb_h,  #elliptic fourier analysis
                              norm=TRUE) #center and scale to harmonic 1 unit size
  
  if(viz){ #visualization
    if(single){
      warning("No visualisation allowed for a single contour")
    } else {
      par(mfrow=c(1,1))
      outs_rs <- outs_rs %>% coo_center() %>% coo_scale()
      stack(outs_rs) #outlines superposition
      outs_rs %>% panel() #all outlines
      outs_rs %>% calibrate_harmonicpower_efourier(nb.h=nb_h) #elliptic fourier calibration
      plot_PCA(PCA(efa))
    }
  }

  return(efa)
}

#' Prediction pipeline for sex or wings
#'
#' Wrapper for gBodyPredict1, allowing for format control and vectorization
#' From segmented body binary images: contour, orientation, resampling, 
#' 
#' @param img_data c.img or list of c.img of segmented bodies
#' @param angle angle, vector of angles or list of angles of the segmented bodies
#' @param what "sex" or "wing". Data to be predicted.
#' @export
gBodyPredict <- function(img_data, angle, what){
  #Load required models and set corresponding parameters
  if(what=="sex"){
    rs <- 20 #group resampling
    pca <- loadModel("PCAsex") #group PCA
    lda <- loadModel("LDAsex") #group LDA
  } else if(what=="wing"){
    rs <- 40
    pca <- loadModel("PCAwing")
    lda <- loadModel("LDAwing")
  } else stop("'what' should be 'sex' or 'wing'")
  
  #Vectorization with safety
  if(is.list(img_data)){ #img_data is list
    ang_len <- length(angle)
    if(ang_len == length(img_data)){ #angle has more than 1 element
      if(ang_len == length(img_data)){ #matching length
        valid <- !( sapply(angle, function(x) all(is.na(x))) |
                    sapply(img_data, function(x) all(is.na(x))) ) #non na points, boolean
        if(sum(valid) >= 1){ #at least 1 valid point
          res <- vector(mode="list", length=ang_len)
          valid_res <- mapply(function(img, ang){ #make predictions
            gBodyPredict1(img=img, ang=ang, rs=rs, pca=pca, lda=lda)
          }, img_data[valid], angle[valid], SIMPLIFY=FALSE)
          res[valid] <- valid_res #index preservation !
          
          #output format
          res <- list(  #NA for non valid ang-img pair
            class=rep(NA,ang_len),
            posterior=data.frame(matrix(NA, nrow=ang_len, ncol=2)),
            x=rep(NA,ang_len)
          )
          colnames(res$posterior) <- colnames(valid_res[[1]]$posterior)
          res$class[valid] <- sapply(valid_res, \(x) x$class)
          res$posterior[valid, ] <- do.call(rbind, lapply(valid_res, \(x) x$posterior))
          res$x[valid] <- sapply(valid_res, \(i) i$x)
          
        } else {
          stop("no valid (non-NA) img-angle pair")
        }
      } else {
        stop("img_data and angle should have the same length")
      }
    } else { #img_data is list but only 1 angle is provided
      stop("img_data is list but do not have the same length as provided angle(s)") 
    }
  } else { #Predict single
    res <- gBodyPredict1(img=img_data, ang=angle, rs=rs, pca=pca, lda=lda)
  }
  
  return(res)
}

#' Prediction pipeline for sex or wings - single image, to be used in gBodyPredict
#'
#' @param img a c.img object. Segmented binary waterstrider body image.
#' @param ang a number. Orientation angle of that individual.
#' @param rs an integer. 
#' @param pca a prcomp object.
#' @param lda a MASS lda object.
#' @import MASS
#' @export
gBodyPredict1 <- function(img, ang, rs, pca, lda){
  #safety
  if(!imager::is.cimg(img)) stop("img should be a c.img")
  if(!is.numeric(ang)) stop("ang should be a number")
  
  #Pre-processing data
  ori_img <- normBody(img, ori_angle = ang) #orientation
  H <- pca$rotation %>% nrow /4 #number of EFA harmonics for the PCA
  PC <- lda$means %>% ncol #number of PC for the LDA
  efa <- scoresEFA(ori_img, nb_h=H, viz=FALSE, resamp=rs)
  pca_scores <- pca %>% predict(newdata = efa$coe)
  #predictions
  predictions <- lda %>% predict(newdata = pca_scores[,1:PC])
  return(predictions)
}
