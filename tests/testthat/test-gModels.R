#loadModel
test_that("Models can be loaded", {
  lda <- loadModel("LDAsex")
  pca <- loadModel("PCAsex")
  expect_equal(inherits(pca, "prcomp"), TRUE)
  expect_equal(inherits(lda, "lda"), TRUE)
})

#normBody + scoresEFA
test_that("Elliptic harmonics extraction",{
  img_b <- readRDS(test_path("fixtures", "img_body.rds"))
  ori_cont <- normBody(img_bin=img_b, ori_angle=-1.07)
  efa <- scoresEFA(ori_cont, nb_h=3, viz=FALSE, resamp=20)
  expect_equal(inherits(efa, "OutCoe"), TRUE)
  expect_equal(inherits(efa, "Coe"), TRUE)
  expect_equal(efa$coe %>% length, 12)
  expect_equal(efa$coe[2] %>% round(., 7), -0.003956)
})

#everything in gModels involving sex prediction
test_that("correct feature prediction",{
  img_b <- readRDS(test_path("fixtures", "img_body.rds"))
  res <- gBodyPredict(img_data=img_b, angle=-1.07, "sex")
  expect_equal(res$posterior[,"F"] %>% as.numeric, 1)
})


