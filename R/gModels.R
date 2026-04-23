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
#' Orient body contours extracted from image using provided angles then define
#' consistent contour starting point as leftmost point
#'
#' @param img_bin c.img or list of c.img of waterstrider body image(s)
#' @param ori_angle angle, vector of angles or list of angles of the segmented
#'   bodies
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
  aligned_cont <- rbind(rot_cont[start_idx:nrow(rot_cont), ],
                        rot_cont[1:(start_idx - 1), ])   # Start from leftmost point (minimum X) to normalize contour insertion thus avoid 180° random rotations
  
  return(aligned_cont)
}

#' Momocs EFA from waterstrider body image
#'
#' @param ori_cont contour or list of contours (as 2D numeric dataframe of
#'   matrix) of waterstrider body image(s)
#' @param nb_h number of harmonics to compute
#' @param viz logical for plot of alignment, fourier space and PCA
#' @param resamp number of points to resample on each contour
#' @importFrom Momocs coo_scale coo_center coo_interpolate efourier panel
#'   calibrate_harmonicpower_efourier PCA plot_PCA
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
  outs_rs <- Momocs::coo_interpolate(outs, n = resamp)  #resample
  efa <- Momocs::efourier(outs_rs, nb.h=nb_h,  #elliptic fourier analysis
                              norm=TRUE) #center and scale to harmonic 1 unit size
  
  if(viz){ #visualization
    if(single){
      warning("No visualization implemented for a single contour")
    } else {
      par(mfrow=c(1,1))
      outs_rs <- outs_rs %>% coo_center() %>% coo_scale()
      stack(outs_rs) #outlines superposition
      outs_rs %>% panel() #all outlines
      outs_rs %>% calibrate_harmonicpower_efourier(nb.h=nb_h) #elliptic fourier calibration
      plot_PCA(Momocs::PCA(efa))
    }
  }

  return(efa)
}

#' Body contour normalization pipeline for LDA models
#' 
#' @param body c.img. Segmented body.
#' @param ang numerical. Body elongation angle.
#' @param resamp positive integer. Resampling size of body contour.
#' @export
gNorm <- function(body, ang, resamp){
  norm_cont <- normBody(img_bin = body, ori_angle = ang)
  if(!is.list(norm_cont)) norm_cont <- list(norm_cont) #handle single
  outs <- Momocs::Out(norm_cont) #outline Momocs object
  outs_rs <- outs %>%
    Momocs::coo_interpolate(n = resamp) %>%
    Momocs::coo_center() %>%
    Momocs::coo_scale()
  return(outs_rs)
}

#' Predict direction of head of Microvelia longipes
#'
#' Uses body elongation axis, limb-body intersections and body EFA coefficients
#' to predict direction of head as angle from centroid using trained LDA model.
#'
#' @param body c.img. Segmented body.
#' @param dil_cont numerical matrix of data frame. Dilated body contour.
#' @param inter_idx boolean vector matching dil_cont rows. Limb intersections 
#' with dilated body contour.
#' @param el_angle numerical in radians. Body elongation angle.
#' confidence threshold under which returned angle should be NA.
#' @param viz boolean. Visualization option.
#' @export
gDirection <- function(body, dil_cont, inter_idx, el_angle, viz = TRUE){
  #Intersection features
  inter_feat <- itxFeature(dil = dil_cont, int = inter_idx,
                           ang = el_angle)
  
  #Shape features
  outs_rs <- gNorm(body, el_angle, resamp = 100)
  pca_mod <- loadModel("PCAori_mini")
  n_harm <- nrow(pca_mod$rotation)/4 #number of EFA harmonics for the PCA
  efa <- scoresEFA(outs_rs, nb_h = n_harm, resamp = 100)
  pred_pca <- predict(pca_mod, newdata = efa$coe)
  
  #LDA Predictions
  lda_mod <- loadModel("LDAori")
  lda_names <- colnames(lda_mod$means)
  n_comp <- sum(grepl("PC", lda_names))
  feat_ori <- cbind(pred_pca[,1:n_comp, drop = FALSE], inter_feat)
  if(!all(colnames(feat_ori)==lda_names)){
    warning("Body direction LDA prediction feature build failed, returning NAs")
    return(rep(NA, length(el_angle)))
  }
  pred <- suppressWarnings(predict(lda_mod, newdata = as.data.frame(feat_ori)))
  to_redirect <- ifelse(pred$class=="R", 1, 0) #Correct orientation
  ang <- ((unlist(el_angle) + to_redirect*pi) %% (2*pi)) - pi
  conf_pred <- pred$posterior > 0.95 #Confidence filter
  conf_ok <- conf_pred[,"R", drop = FALSE] | conf_pred[,"L", drop = FALSE]
  ang[!conf_ok] <- NA
  
  if(viz){
    reorc <- lapply(1:length(outs_rs[[1]]), function(i){
      if(to_redirect[i]==0){
        return(outs_rs[[1]][[i]] %*% rotMat(pi))
      }
      else return(outs_rs[[1]][[i]])
    })
    outs_rsCC <- outs_rs
    outs_rsCC[[1]] <- reorc
    par(mfrow=c(1,2))
    stack(outs_rs,title="Before")
    stack(outs_rsCC,title="After") 
  }
  
  return(ang)
}

#' Predict sex and presence of wings of Microvelia longipes
#' 
#' Use body EFA coefficients to predict sex and presence of wings using
#' trained LDA model
#' 
#' @param body a c.img (or list of this type). Segmented body.
#' @param angle a numerical (or vector of this type). Orientation angle of 
#' corresponding body
#' @importFrom MASS lda
#' @export  
gPredict <- function(body, angle){
  if(length(body)==1 && is.list(body)) body <- body[[1]]
  if(is.list(angle)) angle <- unlist(angle)
  #Preprocessing
  outs_rs <- gNorm(body = body, ang = angle - pi, resamp = 100) #trained on reversed angles
  pca_mod <- loadModel("PCAbio_mini")
  n_harm <- nrow(pca_mod$rotation)/4 #number of EFA harmonics for the PCA
  efa <- scoresEFA(outs_rs, nb_h = n_harm, resamp = 100)
  pred_pca <- predict(pca_mod, newdata = efa$coe)
  #Prediction
  lda_mod <- loadModel("LDAbio")
  lda_names <- colnames(lda_mod$means)
  n_comp <- sum(grepl("PC", lda_names))
  pred <- predict(lda_mod, newdata = as.data.frame(pred_pca[,1:n_comp, drop = FALSE]))
  split_pred <- strsplit(as.character(pred$class), split="[.]") #format
  df_pred <- do.call(rbind, split_pred)
  #Confidence
  post <- pred$posterior
  pred_names <- colnames(post)
  f_id <- grep("F", pred_names) #index of female columns
  w_id <- grep("1", pred_names)
  f_prob <- rowSums(post[,f_id, drop = FALSE]) #cumprob for "F"
  w_prob <- rowSums(post[,w_id, drop = FALSE])
  sex_prob <- 0.5 + abs(f_prob - 0.5) #confidence in sex prediction
  wing_prob <- 0.5 + abs(w_prob - 0.5)
  #Output
  na_ang <- is.na(angle) #NA handling
  full_df <- as.data.frame(matrix(NA, nrow = length(na_ang), ncol = 4))
  colnames(full_df) <- c("sex", "wing", "prob_sex", "prob_wing")
  full_df[!na_ang, ] <- cbind(df_pred, sex_prob, wing_prob)
  return(full_df)
}

#' Body-limbs intersection feature builder for orientation prediction model.
#'
#' Creates feature vector with: Ratio of intersection points in 3 equal-sized
#' bins along dilated body elongation axis, and number of continuous
#' intersection sequences along dilated body contour.
#'
#' @param dil list of dilated contour
#' @param int list of limbs intersection index on dilcont
#' @param ang list of body elongation axis angle 
itxFeature <- function(dil, int, ang){
  if (is.data.frame(dil) || is.matrix(dil)) dil <- list(dil) #as list if single
  if (!is.list(int)) int <- list(int)
  if (length(ang)==1) ang <- list(ang)
  rdil <- mapply(function(d, r) as.matrix(d)%*%r, dil,
                 lapply(ang, rotMat), SIMPLIFY = FALSE)
  f1 <- mapply(function(dd,itx){
    d <- dd[,1]
    ran <- range(d)
    dr <- diff(ran)
    cuts <- cut(d[itx], breaks=ran[1]+c(-1, dr*0.333, dr*0.667, dr))
    return(table(cuts)/sum(itx))
  }, rdil, int, SIMPLIFY=FALSE)
  f2 <- sapply(int, \(x){
    if(length(unique(x))==1) 0
    else boolSeqLim(x)$first %>% length
  })
  dfout <- cbind(do.call(rbind,f1), f2)
  colnames(dfout) <- c("itxbin1","itxbin2","itxbin3","n_itx")
  return(dfout)
}