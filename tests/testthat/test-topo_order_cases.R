test_that("topo_order_cases returns upstream-to-downstream order", {
  routes <- data.frame(
    from_case = c("A","B"),
    stream    = c("outflow","outflow"),
    to_type   = c("CASE","CASE"),
    to_id     = c("B","C"),
    frac      = 1,
    lag_days  = 0L,
    stringsAsFactors = FALSE
  )
  ord <- topo_order_cases(routes, case_names = c("A","B","C"))
  expect_equal(ord, c("A","B","C"))
})

test_that("topo_order_cases errors on cycles", {
  routes <- data.frame(
    from_case = c("A","B"),
    stream    = c("outflow","outflow"),
    to_type   = c("CASE","CASE"),
    to_id     = c("B","A"),
    frac      = 1,
    lag_days  = 0L,
    stringsAsFactors = FALSE
  )
  expect_error(topo_order_cases(routes, case_names = c("A","B")), "cycle")
})
