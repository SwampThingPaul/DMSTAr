test_that("latent storage is disabled when inflow is present", {

  params <- list(
    A_cell = 100,
    Zinit = 0,
    Zmin = 0,
    Zrelease = -50,
    dmsta_version = "2C2B",
    enable_latent_storage = TRUE,
    Bypass_elev = 9999,
    Q_zmin = 0,
    Zweir = 0,
    Seepin_Rate = 0,
    Seepout_Rate = 0,
    IsaNode = FALSE
  )

  inputs <- list(
    Date = as.Date("2000-01-01"),
    Qi = 1,   # inflow breaks latent condition
    Rain = 0,
    Et = 0,
    Zcontrol = 0,
    Qr0 = 0,
    Qr1 = 0,
    Qr2 = 0
  )

  out <- dmsta_flow_day(
    V = 0,
    inputs = inputs,
    params = params,
    Qmethod = "RK4"
  )

  expect_false(identical(out$meta$method, "LATENT"))
})
