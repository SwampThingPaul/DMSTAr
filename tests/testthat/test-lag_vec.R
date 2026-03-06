test_that("lag_vec shifts forward with zero padding", {
  x <- 1:5
  expect_equal(lag_vec(x, 0), x)
  expect_equal(lag_vec(x, 2), c(0,0,1,2,3))
})

test_that("lag_vec shifts backward with zero padding", {
  x <- 1:5
  expect_equal(lag_vec(x, -2), c(3,4,5,0,0))
})
