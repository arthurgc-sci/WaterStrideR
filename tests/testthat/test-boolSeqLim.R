#create test
test_that("sequences of TRUE values are found in circular vector", {
  expect_equal(boolSeqLim(c(1,1,0,1,0,1,1),circular=TRUE), list(first=c(4,6),last=c(4,2)))
})

           