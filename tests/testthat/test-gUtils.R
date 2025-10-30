#boolSeqLim
test_that("sequences of TRUE values are found in circular vector", {
  expect_equal(boolSeqLim(c(1,1,0,1,0,1,1),circular=TRUE),
               expected = list(first=c(4,6),last=c(4,2)))
  expect_equal(boolSeqLim(FALSE), list(first=NULL,last=NULL))
  expect_equal(boolSeqLim(TRUE), list(first=1,last=1))
})

#simpleCont
test_that("Simple binary contouring works",{
  img_leg <- readRDS(test_path("fixtures", "img_leg.rds"))
  cont_leg <- readRDS(test_path("fixtures", "cont_leg.rds"))
  expect_equal(simpleCont(img_leg), cont_leg)
  expect_equal(simpleCont(imager::as.pixset(img_leg)), cont_leg) #with format dif
})

#major_axis_angle
test_that("PCA and atan2 angle conversion",{
  cont_bod <- readRDS(test_path("fixtures", "cont_body.rds"))
  expect_equal(major_axis_angle(cont_bod),2.020087923)
})

#normV, scalar, angleV
test_that("basic angle operations work",{
  v1 <- c(4,3)
  v2 <- c(1,2)
  expect_equal(normV(v1),5)
  expect_equal(scalar(v1,v2),10)
  expect_equal(angleV(v1,v2),0.463647609)
})

#contourAngles + contourTurns
test_that("Angle computation along contour: contourAngles + contourTurns",{
  cont_bod <- readRDS(test_path("fixtures", "cont_body.rds"))
  contang <- contourAngles(cont_bod, search_w=10)
  expect_equal(mean(contang) %>% round(.,7),0.1189731)
  contturn <- contourTurns(cont_bod, search_w=5, splines_df=30, angle_thresh=0.15, viz=0)
  expect_equal(contturn, c(5, 37, 53, 72, 132, 175, 193, 213, 271))
})

#bodyLength + imgAsCoords
test_that("Images can be converted to points for measurments",{
  img_bod <- readRDS(test_path("fixtures", "img_body.rds"))
  pts_bod <- imgAsCoords(img_bod)
  expect_equal(bodyLength(pts_bod) %>% round(.,7), 136.841084)
})

#cleanBodyShape
test_that("Body shape filter with different output type",{
  pts_bod <- readRDS(test_path("fixtures", "pts_body.rds"))
  cbs1 <- cleanBodyShape(pts_bod, kernel_size=10, as_coords = TRUE)[[1]]
  cbs2 <- cleanBodyShape(pts_bod, kernel_size=10, as_coords = FALSE)[[1]]
  expect_equal(nrow(cbs1), sum(cbs2==1))
  expect_equal(nrow(cbs1), 5152)
})

#overlapPoints
test_that("overlapPoints flexibility",{
  set_a1 <- c(1,1)
  set_b1 <- c(1,1)
  set_b2 <- c(2,1)
  set_a2 <- c(7,2,2,4,1,3) %>% matrix(ncol=2)
  set_b3 <- c(7,2,2,0,0,4,1,3,0,0) %>% matrix(ncol=2)
  expect_equal(overlapPoints(set_a1, set_b1)$index,TRUE)
  expect_equal(overlapPoints(set_a1, set_b2)$index,FALSE)
  expect_equal(overlapPoints(set_a2, set_b2)$index,c(FALSE, TRUE, FALSE))
  expect_equal(overlapPoints(set_b3, set_a2)$index,c(TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(overlapPoints(set_a2, set_b3)$index,c(TRUE, TRUE, TRUE))
})

#make_disc_kernel
test_that("Kernel disc creation",{
  expect_equal(make_disc_kernel(4), matrix(c(0,1,1,1,0,rep(1,15),0,1,1,1,0),ncol=5) )
})
