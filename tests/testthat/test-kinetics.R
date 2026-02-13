test_that("compute_DMSTA_kvals returns PModel=1 for standard inputs", {
  K <- compute_DMSTA_kvals(C1000 = 1000, Cstar = 100, Ks = 1)
  expect_equal(K$PModel, 1L)
  expect_true(all(c("K1","K2","K3","PModel") %in% names(K)))
})

test_that("compute_DMSTA_kvals returns PModel=2 for transformed kinetics", {
  K <- compute_DMSTA_kvals(C1000 = 1000, Cstar = -0.2, Ks = 1)
  expect_equal(K$PModel, 2L)
})

test_that("build_P_kin_slots returns aligned K-length vectors", {
  # minimal params list; must include whatever your builders require
  p <- list(
    DutyCycle = 0.95,
    # STA inputs
    C1000 = 1000, Cstar = 100, Ks_per_yr = 1,
    Z1 = 10, Z2 = 20, Z3 = 30, K2Coef1 = 0.1, Chalf = 50, SeasonalFactor = 1,
    # PSTA inputs
    Ytrans = 1, Ysigma = 1, C1000_2 = 1000, ks_2 = 1, zh_2 = 10,
    # RES inputs
    k_depth_penalty = 0.5
  )

  out <- build_P_kin_slots(
    mods = c("STA","PSTA","RES"),
    pparams = p,
    Dpy = 365.25,
    DutyCycle = p$DutyCycle
  )

  expect_true(is.list(out))
  expect_equal(out$Kslots, 3L)
  validate_P_paramsK(out)  # should not error

  expect_length(out$K1, 3)
  expect_length(out$K2, 3)
  expect_length(out$K3, 3)
})
