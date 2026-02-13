test_that("example dataset series loads and has expected columns", {
  data(series, package = "DMSTAr")  # adjust package name
  expect_true(is.data.frame(series))
  expect_true(all(c("Date","Flow","Conc","Rainfall","ET") %in% names(series)))
  expect_true(nrow(series) > 0)
})
