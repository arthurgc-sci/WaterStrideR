#redRulerScale
test_that("Red ruler detection", {
  img_scale <- readRDS(test_path("fixtures", "img_scale.rds"))
  scale <- suppressWarnings(redRulerScale(img_scale, red_thresh=0.05, alpha_confidence=0.05, viz=TRUE, msg=FALSE))
  expect_equal(scale %>% as.numeric, c(0.040322547, 0.001407644))
})
