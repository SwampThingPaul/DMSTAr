test_that("run_network_of_cases runs and returns expected structure", {
  cases <- make_cases_2(n = 10)
  net   <- make_net_table_simple()

  out <- run_network_of_cases(
    cases = cases,
    net_table = net,
    outlet_count = 5L,
    verbose = FALSE,
    # Qmethod = "RK4"
    # plus any dmsta_flowP_case args you want to standardize in tests
  )

  expect_true(all(c("order","case_results","routed_in","outlets","outlet_summary","ledger") %in% names(out)))

  # Order should place UP before DN
  expect_equal(out$order, c("UP","DN"))

  # Each case result should be a list-like dmsta_flowP_case return
  expect_true(is.list(out$case_results$UP))
  expect_true(is.list(out$case_results$DN))

  # routed_in contains Q/L vectors length = nday
  expect_equal(length(out$routed_in$UP$Q), 10)
  expect_equal(length(out$routed_in$DN$L), 10)

  # outlets is list length outlet_count
  expect_equal(length(out$outlets), 5)
  expect_equal(length(out$outlets[[1]]$Q), 10)

  # ledger should have one row per configured route
  expect_s3_class(out$ledger, "data.frame")
  expect_equal(nrow(out$ledger), 2)
})

test_that("run_network_of_cases errors on unknown from_case", {
  cases <- make_cases_2(n = 5)
  routes <- data.frame(
    from_case="NOPE", stream="outflow", to_type="OUTLET", to_id=1,
    frac=1, lag_days=0L, stringsAsFactors=FALSE
  )
  expect_error(run_network_of_cases(cases, routes=routes, verbose=FALSE), "Unknown from_case")
})

test_that("run_network_of_cases errors on unknown downstream CASE", {
  cases <- make_cases_2(n = 5)
  routes <- data.frame(
    from_case="UP", stream="outflow", to_type="CASE", to_id="NOPE",
    frac=1, lag_days=0L, stringsAsFactors=FALSE
  )
  expect_error(run_network_of_cases(cases, routes=routes, verbose=FALSE), "Unknown downstream CASE")
})
