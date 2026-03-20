test_that("dmsta_flow_series returns expected structure and class", {
  # Minimal series (3 days)
  series <- data.frame(
    Date = as.Date("2000-01-01") + 0:2,
    Qi = c(1, 1, 1),
    Rain = c(0, 0, 0),
    Et = c(0.005, 0.005, 0.005),
    Qr0 = 0,
    Qr1 = 0,
    Qr2 = 0,
    Zcontrol = c(1, 1, 1)
  )

  # Minimal params required by your hydro code paths
  params <- list(
    A_cell = 10,
    Zmin = 10,     # cm
    Q_a = 0, Q_b = 1,
    Zweir = 100, Q_zmin = 0,
    Width = 1,
    Bypass_elev = 0,
    Seepout_Rate = 0, Seepout_Elev = 0,
    Seepin_Rate = 0,  Seepin_Elev = 0,
    Qin_Frac = 1,
    IsaNode = NULL
  )

  out <- dmsta_flow_series(V_init = 10, series = series, params = params, Nsteps = 4L)

  expect_s3_class(out, "dmsta_hydro_result")
  expect_true(is.data.frame(out$results))
  expect_true(is.list(out$budgets))
  expect_true(is.list(out$meta))
  expect_true(all(c("Date","Qout","V_end") %in% names(out$results)))
})
