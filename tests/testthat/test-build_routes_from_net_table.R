
test_that("build_routes_from_net_table parses CASE and OUTLET destinations", {
  net <- make_net_table_simple()
  routes <- build_routes_from_net_table(net, outlet_count = 5L, case_col = "CaseName")

  expect_s3_class(routes, "data.frame")
  expect_true(all(c("from_case","stream","to_type","to_id","frac","lag_days") %in% names(routes)))

  # Expect 2 routes: UP outflow->DN (CASE), DN outflow->1 (OUTLET)
  expect_equal(nrow(routes), 2)
  expect_equal(routes$from_case, c("UP","DN"))
  expect_equal(routes$stream,    c("outflow","outflow"))
  expect_equal(routes$to_type,   c("CASE","OUTLET"))
  expect_equal(as.character(routes$to_id), c("DN","1"))
})

test_that("build_routes_from_net_table rejects missing required columns", {
  net <- data.frame(CaseName = "A")
  expect_error(build_routes_from_net_table(net), "missing columns")
})

test_that("build_routes_from_net_table validates outlet range", {
  net <- make_net_table_simple()
  net$Outflow_to[2] <- "99"  # invalid
  expect_error(build_routes_from_net_table(net, outlet_count = 5L), "out of range")
})

test_that("build_routes_from_net_table ignores blanks and NA", {
  net <- make_net_table_simple()
  net$Outflow_to[1] <- "   "
  net$Seepage_to[1] <- NA
  routes <- build_routes_from_net_table(net, outlet_count = 5L)

  # Only DN->Outlet1 remains
  expect_equal(nrow(routes), 1)
  expect_equal(routes$from_case, "DN")
  expect_equal(routes$to_type, "OUTLET")
})
