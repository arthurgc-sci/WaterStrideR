#gFastSeg
test_that("Connected components segmentation works - gFastSeg", {
  sq2 <- readRDS(test_path("fixtures", "img_2square.rds"))
  res <- gFastSeg(sq2, viz=FALSE, px_range=c(1,15))
  expect_equal(length(res),2)
  expect_equal(lapply(res, nrow) %>% unlist %>% as.numeric, c(9,8))
})

#gCleanBin
test_that("Small or edge-connex components removal from binary image", {
  sq3 <- readRDS(test_path("fixtures", "img_3square.rds"))
  res <- gCleanBin(l_img_bin=sq3, centroid=c(6,6), px_filter = 9)
  expect_equal(nrow(res),9)
})

#dilBodies
test_that("Size based body image dilation", {
  body_img <- readRDS(test_path("fixtures", "img_body.rds"))
  L <- bodyLength(imgAsCoords(body_img))
  expect_equal((dilBodies(body_img, L)==1) %>% which %>% length, 12953)
})

#gLegLandmarks + simpleCont + coordsAsImg + contourTurns + overlapPoints + 
test_that("Leg landmarks pipeline works", {
  pts_leg <- readRDS(test_path("fixtures", "img_leg.rds")) %>% imgAsCoords
  leg_lm <- gLegLandmarks(leg_coords=pts_leg,
                          insertion=c("dim1"=53, "dim2"=22),
                          inser_thresh=0.1,
                          tresh_ankle=0.12,
                          viz=FALSE,
                          msg=FALSE,
                          viz_angle=FALSE,
                          inflexion_pts_range = c(5,7),
                          knee_diff_thresh = 0.2,
                          segment_length_range = c(0.15, 0.6),
                          n_splines = 30,
                          search_w = 20)
  expect_equal(leg_lm %>% sum, 356.5)
})

#gConnectLeg
test_that("Reconnect insertion to non-dilated body", {
  pts_bod <- readRDS(test_path("fixtures", "pts_body.rds"))
  leg_lm <- data.frame(dim1=c(53,95,24),dim2=c(23,48,113.5))+50
  rownames(leg_lm) <- c("insertion","knee","ankle")
  res <- gConnectLeg(pts_bod,leg_lm)
  expect_equal(res["insertion",] %>% as.numeric, c(95,68))
})

#gMeasureLeg + gSplitLegMeasures
test_that("Leg measurment and output format", {
  rand6 <- c(-4,5,9,13,-12,-20)
  leg_lm1 <- list(
    right = data.frame(dim1=c(53,95,24) %>% as.numeric,
                       dim2=c(23,48,113.5) %>% as.numeric) + rand6,
    left = data.frame(dim1=-c(53,95,24) %>% as.numeric,
                      dim2=(c(23,48,113.5)) %>% as.numeric) + rand6
  )
  for(i in 1:2) rownames(leg_lm1[[i]]) <- c("insertion","knee","ankle")
  leg_lm2 <- list(
    right = data.frame(dim1=c(53,95,24) %>% as.numeric,
                       dim2=c(23,48,113.5) %>% as.numeric),
    left = data.frame(dim1=-c(53,95,24) %>% as.numeric,
                      dim2=(c(23,48,113.5)) %>% as.numeric)
  )
  for(i in 1:2) rownames(leg_lm2[[i]]) <- c("insertion","knee","ankle")
  leg_lm <- list(leg_lm1, leg_lm2)
  leglen <- gMeasureLeg(landmarks=leg_lm, scale=10)
  splitres <- gSplitLegMeasures(leglen)
  expect_equal(length(leglen)==2 &&
               length(leglen[[1]])==2, TRUE)
  expect_equal(leglen[[2]][[2]] %>% as.numeric %>% round(.,5),
               c(4.88774,9.65984))
  expect_equal(length(splitres)==4 &&
               length(splitres[[1]]==2), TRUE)
  expect_equal(splitres[[2]] %>% as.numeric %>% round(.,5), c(8.82907, 9.65984))
})
