#functions related to manual scaling

#' Open a Graphics Device Depending on OS
#'
#' Opens an interactive graphics device suitable for the current operating
#' system.
#' - On **macOS**, uses `quartz()`
#' - On **Linux** and **Windows**, uses `X11()`
#'
#' @param title a string. Graphic window's name
#' @return `TRUE` if the graphics device was successfully opened, `FALSE`
#' otherwise.
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
    warning("Failed to open graphic device: ", conditionMessage(e))
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
      scale_value <- readline(prompt = "\U1F4D0 Enter scale value: ")
      num <- suppressWarnings(as.numeric(scale_value)) #NA if not a number
      if(!is.na(num)){ #if input is a number
        if(num > 0){ #if it's positive
          break
        }
      }
      cat("\u274C Please enter a positive number") 
    }
  } else if(!is.numeric(scale_value) || length(scale_value) != 1){ #if not a number
    stop("\u274C scale_value must be a single numeric value")
  }
  # 2 - Get scale pixel size
  if(!openGrDevice(title = "Click two points to set scale (this window can then be closed)")){ #opening interactive windows to bypass errors with locator() in Rstudio plot
    message("Using plot to set scale 
    If you are using Rstudio please set
    Tools>Global options>Appearance>Zoom = 100%") #using local plot if opening graphic device fails
  }
  cex_min <- cex_max*.6 #colored point size (20% black outline)
  par(mar=rep(0,4))
  plot(img, ann=FALSE, asp=1)
  scale_pt1 <- graphics::locator(1) #click point on image
  points(scale_pt1, pch=16, col=1, cex=cex_max)
  points(scale_pt1, pch=16, col=plot_col, cex=cex_min)
  scale_pt2 <- graphics::locator(1)
  points(scale_pt2, pch=16, col=1, cex=cex_max)
  points(scale_pt2, pch=16, col=plot_col, cex=cex_min)
  graphics::segments(x0 = scale_pt1[[1]], y0 = scale_pt1[[2]],
                     x1 = scale_pt2[[1]], y1 = scale_pt2[[2]],
           col = plot_col, lwd = 1)
  dist_px <- (unlist(scale_pt1)-unlist(scale_pt2))^2 %>% sum %>% sqrt #computing distance
  return(c(scale=dist_px/as.numeric(scale_value),
           scale_value=scale_value) %>% as.numeric) #pixel per unit
}

#' Red ruler tiles centroids
#'
#' From image with red ruler and no other big red object : returns dataframe of
#' centroids of each tiles of the scale
#'
#' @param img c.img including a red ruler and no other red object bigger than
#'   it.
#' @param red_thresh required difference in intensity between red component and
#'   other components to detect the red ruler
#' @param viz a boolean. visualization option
#' @param msg a boolean. Console output option
redRulerTiles <- function(img, red_thresh, viz=FALSE, msg=TRUE){
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
  ruler_crop <- imager::crop.bbox(ruler_img, ruler_img==TRUE) #cropped ruler
  ruler_neg <- ruler_crop %>% invert_grayscale #negative
  seg_ruler <- gFastSeg(ruler_neg, px_range=c(1,100000000), viz=viz) #negative ruler segmentation
  l_seg_r <- lapply(seg_ruler, nrow) %>% unlist #label's size
  median_threshold <- which(l_seg_r < 1.5 * median(l_seg_r) &
                              l_seg_r > 0.5 * median(l_seg_r))
  seg_tiles <- seg_ruler[median_threshold] #removes label too big or too small
  if(length(seg_tiles) == 0) stop("Error : Ruler tiles not found")
  
  # Distance between tiles centroids
  if(msg) message("--- scaling - Computing tile's centroid distances...")
  mean_px_val <- lapply(seg_tiles, colMeans) #vector of xy mean position of the ruler's tiles
  tile_pos <- as.data.frame(do.call(rbind, mean_px_val)) #dataframe with xy mean position of the ruler's tiles
  
  if(msg) message(paste0("--- scaling - Found ", nrow(tile_pos), " tiles."))
  if(viz){
    points(tile_pos, pch=16, col=0)
    points(tile_pos, pch=10, col=1)
  }
  
  return(tile_pos)
}

#' 4 Nearest neighbor distance for set of points
#'
#' @param tiles_centroid tiles centroid as output of redRulerTiles()
#' @param alpha confidence interval level
#' @param viz a boolean. visualization option
tilesUnitDist <- function(tiles_centroid, alpha=0.05, viz=FALSE){
  knn_res <- FNN::get.knn(tiles_centroid, k = 4) #4 nearest neighbod using KD-tree
  v4 <- as.vector(knn_res$nn.dist) #flatten
  med <- stats::quantile(v4, 0.3) #median for a first average estimate of 1mm with low leverage
  lo_f <- med * 0.5 #filter small outliers likely residues of segmentation error
  hi_f <- med * 1.207 #med+(sqrt(med)-med)/2 = mean value between 1mm and 1mm square diagonal
  v4f0 <- v4[v4 > lo_f & v4 < hi_f] 
  v4f <- v4f0[!duplicated(v4f0)]  #remove duplicates
  unit_dist_mean <- mean(v4f)
  #confidence interval
  n <- length(v4f)
  if(n < 120) warning(paste0("Only ", n, " measurements were found; at least 120 are required to keep the scale error margin below 1%. Try adjusting 'f.red_thresh' in gPipeline(), increasing the scale size, or scaling manually."))
  sem <- sd(v4f) / sqrt(n) #std error mean
  ci_alpha <- abs(qnorm(alpha/2) * sem) #confidence interval
  
  if(viz){
    h <- hist(v4, breaks = 50, plot = FALSE)
    bar_cols <- ifelse(h$mids >= lo_f & h$mids <= hi_f, "#4CAF50", "#aaaaaa")
    plot(h, col = bar_cols, border = NA, xlab = "Pixels",
         main = "Distance to 4 nearest neighbors of each tile")
    abline(v = c(lo_f, hi_f), col=2, lwd=1)
    legend("topright", bty="n", pch=15, col=c("#4CAF50", "#aaaaaa"),
           legend=c("Kept", "Filtered"))
    abline(v=unit_dist_mean, col="#235426", lwd=2)
  }
  
  return(c("scale" = unit_dist_mean, "error" = ci_alpha))
}

#'Get scale using detection of red millimeter paper
#'
#' Detect ruler as biggest red label
#' 
#' @param img c.img or pixset. base image
#' @param red_thresh required difference in intensity between red component
#' and other components to detect the red ruler
#' @param alpha_confidence a numerical value in range 0:1.
#' alpha level for confidence interval
#' @param viz logical. visualization option
#' @param msg logical. console message option
#' @export
redRulerScale <- function(img, red_thresh, alpha_confidence=0.05,
                          viz=FALSE, msg=TRUE){
  tiles_centroid <- redRulerTiles(img, red_thresh=red_thresh, viz=viz, msg=msg)
  sc <- tilesUnitDist(tiles_centroid, alpha=alpha_confidence, viz=viz)
  scale_um <- sc["scale"]/1000
  ic_um <- sc["error"]/1000
  if(msg){
    sc_txt <- format(scale_um , scientific=TRUE, digits=3) #px/micron
    ic_txt <- format(ic_um , scientific=TRUE, digits=3)
    cat('1 \u00b5m =', sc_txt, "\u00B1", ic_txt, "pixels")
  }
  return(c("scale" = scale_um, "error" = ic_um))
}

#' Scale method selecter for gPipeline
#' 
#' @param img c.img or pixset
#' @param auto_scale boolean for redRulerScale or manual
#' @param red_thresh required difference in intensity between red component
#' @param viz visualization option of redRulerScale
#' and other components to detect the red ruler
pipelineScale <- function(img, auto_scale, red_thresh, viz=FALSE){
  if(auto_scale){
    scale <- redRulerScale(img, msg=FALSE, red_thresh=red_thresh, viz=viz) #get scale
  } else {
    scale1 <- getScale(img)
    scale <- c(scale1, NA) #no error margin : not yet implemented here    
  }
}
