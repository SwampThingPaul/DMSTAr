test_that("latent storage does not trigger when disabled", {

  params <- list(
    A_cell = 100,
    Zinit = 0,
    Zmin = 0,
    Zrelease = -50,
    dmsta_version = "2C2B",
    enable_latent_storage = FALSE,
    IsaNode = FALSE
  )

  inputs <- list(
    Date = as.Date("2000-01-01"),
    Qi = 0,
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
    params = params
  )

  expect_false(identical(out$meta$method, "LATENT"))
})
