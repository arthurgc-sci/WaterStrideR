#redRulerScale
test_that("Red ruler detection", {
  img_scale <- readRDS(test_path("fixtures", "img_scale.rds"))
  scale <- redRulerScale(img_scale, red_thresh=0.05, confidence_interval=0.95, viz=FALSE, msg=FALSE)
  expect_equal(scale %>% as.numeric, c(0.0374664710, 0.0004382986))
})
