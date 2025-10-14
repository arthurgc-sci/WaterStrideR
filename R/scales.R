#functions related to manual scaling

#' Open a Graphics Device Depending on OS
#'
#' Opens an interactive graphics device suitable for the current operating system.
#' - On **macOS**, uses `quartz()`
#' - On **Linux** and **Windows**, uses `X11()`
#'
#' @param title a string. Graphic window's name
#' @return
#' `TRUE` if the graphics device was successfully opened, `FALSE` otherwise.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   if (openGrDevice()) {
#'     plot(1:10)
#'   }
#' }
#' }
openGrDevice <- function(title = "R Graphics Window"){
  if (!interactive()) {
    stop("This function is interactive and must be used in an interactive session")
  }
  tryCatch({
    if(Sys.info()[["sysname"]] == "Darwin") { #MAC
      grDevices::quartz(title = title)
    } else { #windows or linux
      grDevices::X11(title = title)
    }
    return(TRUE)
  }, error = function(e){
    warning("⚠️ Failed to open graphic device: ", conditionMessage(e))
    return(FALSE)
  })
}

#' Manual scaling of an image
#' 
#' Scale imager image or pixset by setting scale value for a distance between
#' two pixels set manually on image
#' 
#' @param img a c.img or pixset
#' @param scale_value scale value. If NA: value will be asked in console
#' @param plot_col string or numerical for plot color
#' @param lwd numerical. line width on graphic windows between scale points.
#' @param cex_max numerical. black point size of scale points.
#' @return number of pixels per unit
#' @export
getScale <- function(img, scale_value = NA, plot_col="#ff9100",
                     lwd = 1, cex_max = 1.2){
  if(!imager::is.pixset(img) && !imager::is.cimg(img)){
    stop("img must be a pixset of c.img")
  }
  #1 - Get scale value
  if(is.na(scale_value)){ #manual scale value in R console
    repeat{
      scale_value <- readline(prompt = "📐 Enter scale value: ")
      num <- suppressWarnings(as.numeric(scale_value)) #NA if not a number
      if(!is.na(num)){ #if input is a number
        if(num > 0){ #if it's positive
          break
        }
      }
      cat("❌ Please enter a positive number") 
    }
  } else if(!is.numeric(scale_value) || length(scale_value) != 1){ #if not a number
    stop("❌ scale_value must be a single numeric value")
  }
  # 2 - Get scale pixel size
  if(!openGrDevice(title = "Click two points to set scale (this window can then be closed)")){ #opening interactive windows to bypass errors with locator() in Rstudio plot
    message("Using plot to set scale 
    If you are using Rstudio please set
    Tools>Global options>Appearance>Zoom = 100%") #using local plot if opening graphic device fails
  }
  cex_min=cex_max*.6 #colored point size (20% black outline)
  par(mar=rep(0,4))
  plot(img, ann=FALSE, asp=1)
  scale_pt1 <- graphics::locator(1) #click point on image
  points(scale_pt1, pch=16, col=1, cex=cex_max)
  points(scale_pt1, pch=16, col=plot_col, cex=cex_min)
  scale_pt2 <- graphics::locator(1)
  points(scale_pt2, pch=16, col=1, cex=cex_max)
  points(scale_pt2, pch=16, col=plot_col, cex=cex_min)
  graphics::segments(x0 = scale_pt1[[1]], y0 = scale_pt1[[2]], x1 = scale_pt2[[1]], y1 = scale_pt2[[2]],
           col = plot_col, lwd = 1)
  dist_px <- (unlist(scale_pt1)-unlist(scale_pt2))^2 %>% sum %>% sqrt #computing distance
  return(dist_px/as.numeric(scale_value)) #pixel per unit
}

#' Accessory for redRulerScaler
#'
#' @param text A string of characters. The text in which to count leading zeros.
#' @param zeros An integer. The number of leading zeros to check for in the `text`. Default is 1.
#' @return The number of leading zeros in the string `text`.
#' @keywords internal
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
#' @param img c.img or pixset. base image
#' @param red_thresh required difference in intensity between red component
#' and other components to detect the red ruler
#' @param confidence_interval a numerical value in range 0:1.
#' confidence interval value to build the error margin from
#' @param viz logical. visualization option
#' @param msg logical. console message option
#' @export
redRulerScale <- function(img, red_thresh=0.05, confidence_interval=0.95, viz=F, msg=T){
  # Get ruler as the biggest set of red pixels from the image
  if(msg) message("--- scaling - Extracting ruler...")
  IR <- correct_illumination(imager::R(img)) #RGB components with linear correction of illumination
  GR <- correct_illumination(imager::G(img))
  BR <- correct_illumination(imager::B(img))
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
    cat('1 \u00b5m =', format(unit_dist_mean , scientific=T, digits=3), #unicode for characters
        "\u00B1", e_margin_format,'px \n')
  }
  return(c("scale" = unit_dist_mean, "error" = e_margin))
}

#' Scale method selecter for gPipeline
#' 
#' @param img c.img or pixset
#' @param auto_scale boolean for redRulerScale or manual
#' @param red_thresh required difference in intensity between red component
#' and other components to detect the red ruler
pipelineScale <- function(img, auto_scale, red_thresh=0.05){
  if(auto_scale){
    scale <- redRulerScale(img, msg=F, red_thresh=red_thresh) #get scale
  } else {
    scale1 <- getScale(img)
    scale <- c(scale1, NA) #no error margin : not yet implemented here    
  }
}
