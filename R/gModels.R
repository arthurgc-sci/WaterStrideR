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

#' PCA scores of EFA of oriented contours from binary image
#'
#' @param img_bin a list of binary image or pixset with a single white object
#' @param ori_angle a list of orientation angles of provided shapes
#' @param nb_h a numerical integer. number of elliptuc fourier harmonics to compute
#' @param viz visualization option, additionally returns 3 plots
#' @export
scoresEFA <- function(img_bin, ori_angle=rep(0,length(img_bin)), nb_h=7, viz=FALSE){
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
    x <- x[c(seq_len(nrow(x)),1),] #close contour
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
  scale_cont <- outs %>% Momocs::coo_scale() #scale
  efa <- scale_cont %>% Momocs::efourier(nb.h=nb_h, norm=TRUE) #elliptic fourier analysis
  efa_pca <- prcomp(efa$coe) #principal components analysis
  if(viz){ #visualization
    par(mfrow=c(1,1))
    scale_cont %>% stack #outlines superposition
    scale_cont %>% Momocs::panel() #all outlines
    scale_cont %>% Momocs::calibrate_harmonicpower_efourier(nb.h=nb_h) #elliptic fourier calibration
    Momocs::plot_PCA(Momocs::PCA(efa))
  }
  #re-adding missing lines
  scores <- efa_pca$x[,1:20]
  df_scores <- data.frame(matrix(NA, nrow = length(ori_angle)
                                 , ncol = ncol(scores)), row.names = names(ori_angle))
  colnames(df_scores) <- colnames(scores)
  df_scores[!na_ang,] <- scores
  
  return(df_scores)
}

# FULL MASS IMPORT REQUIRED SINCE predict.lda IS NOT EXPORTED FROM MASS
#' Predict factor using LDA model on harmonic fourier principal components
#'
#' @param LDA_model LDA model created with MASS, trained on PCA output
#' @param hPCA_scores harmonic fourier principal components as output of prcomp
#' @import MASS
#' @export
gPredictLDA <- function(LDA_model, hPCA_scores){
  if(!inherits(LDA_model, "lda")) stop("LDA_model should be an object of class lda")
  
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

#################NEW

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
  body_cont_mat <- body_cont[c(seq_along(body_cont), 1), ] # Close contour
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
scoresEFA2 <- function(ori_cont, nb_h=7, viz=FALSE,resamp=100){
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
#' @export
gBodyPredict1 <- function(img, ang, rs, pca, lda){
  #safety
  if(!imager::is.cimg(img)) stop("img should be a c.img")
  if(!is.numeric(ang)) stop("ang should be a number")
  
  #Pre-processing data
  ori_img <- normBody(img, ori_angle = ang) #orientation
  H <- pca$rotation %>% nrow /4 #number of EFA harmonics for the PCA
  PC <- lda$means %>% ncol #number of PC for the LDA
  efa <- scoresEFA2(ori_img, nb_h=H, viz=FALSE, resamp=rs)
  pca_scores <- pca %>% predict(newdata = efa$coe)
  #predictions
  predictions <- lda %>% predict(newdata = pca_scores[,1:PC])
  return(predictions)
}
