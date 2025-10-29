#GMM
test_that("Custom Gaussian Mixture Model fit",{
  img_gsolocut <- readRDS(test_path("fixtures", "img_gsolocut.rds"))
  g_val <- img_gsolocut %>% as.numeric #convert cropped images to numeric
  g_val_clean <- g_val[g_val!=0 & g_val!=1] #removes eventual extreme values created by gNobody
  res <- GMM(g_val_clean, n=2, keep_max=TRUE, viz=FALSE, n_noise_sd=0)
  expect_equal(res$df[1,] %>% as.numeric %>% round(.,2), c(0.14, 0.08, 0.49))
  expect_true(inherits(res$model, "Mclust"))
})

#GMMdRoots + f_prime
test_that("GMM Roots Solving",{
  gmm_demo <- readRDS(test_path("fixtures", "GMM_ex.rds"))
  res <- GMMdRoots(gmm_demo$df, 0.5, viz=FALSE, msg=FALSE)
  expect_equal(res %>% round(.,7), c(0.1397560, 0.3267837, 0.4974483))
})
