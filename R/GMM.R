#GMM related functions

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

#' GMM 'V' using Mclust
#'
#' returns the GMM model of a distribution as a numeric vector img_val and a dataframe of the parameters with possibility
#' of removing the one with the highest mean. n integer vector of components to test
#' if n length > 1 : the number of components will be decided using BIC
#'
#' @param img_val numerical gray values of an image
#' @param n numerical integer in range 1:9. number 'G' of components of the GMM
#' @param keep_max logical. to remove rightmost component and consequently
#' adjust weights
#' @param viz visualization option
#' @param n_noise_sd numerical. standard deviation of normal noise to add
#' to the data to help breaking plateaus in he data, not compatible with findPeaks() 
#' @export
GMM <- function(img_val, n=1:9, keep_max=T, viz=F, n_noise_sd=0){
  if(imager::is.cimg(img_val)){
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
#' @param gmm_model_df dataframe ($df) output of GMM()
#' @param q numerical in range 0:1. quantile
#' @param viz visualization option
#' @export
GMMquant <- function(gmm_model_df, q, viz=F){
  f <- function(x){
    apply(gmm_model_df, 1, function(df){
      pnorm(x, mean=df["mu"], sd=df["sd"]) * df["w"]
    }) %>% sum - q
  }
  x_quant <- stats::uniroot(f, c(0,10))[1] %>% unlist #TODO INTERVAL???
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
#' @param GMM_df dataframe ($df) output of GMM()
#' @param y numerical value to find the roots of
#' @param viz visualization option
#' @param msg onsole pessage option
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
#' @param x A numeric vector. The input values for which 
#' the derivative of the Gaussian Mixture Model is to be computed.
#' @param GMM_df dataframe ($df) output of GMM()
#' @keywords internal
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