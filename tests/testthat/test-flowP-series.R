test_that("dmsta_flowP_series runs and returns expected structure", {

  # --- Minimal input series (5 days) ---
  series <- data.frame(
    Date     = as.Date("2000-01-01") + 0:4,
    Qi       = c(1, 1.2, 0.8, 1, 1),     # volume/day
    Ci       = c(50, 60, 55, 50, 52),    # ppb
    Rain     = rep(0, 5),               # m/day
    Et       = rep(0.005, 5),           # m/day
    Zcontrol = rep(1, 5),               # m
    Qr0 = rep(0, 5),
    Qr1 = rep(0, 5),
    Qr2 = rep(0, 5),
    stringsAsFactors = FALSE
  )

  # --- Minimal params required by both hydrology + kinetics ---
  params <- list(
    # cell geometry and timing
    A_cell    = 10,       # area (km^2 in your DMSTA convention)
    DutyCycle = 0.95,
    Qimax = 0, Qomax = 0,

    # initial conditions for hydrology + P
    Zinit       = 100,    # cm
    C_init_ppb  = 50,
    Y_init_mgm2 = 0,

    # hydrology parameters (dmsta_deriv_flow expects these names)
    Zmin   = 10,          # cm -> converted to m internally
    Q_a    = 0,           # 0 => "volume control" style outflow when Z > Zcont
    Q_b    = 1,
    Zweir  = 100,         # cm
    Q_zmin = 0,           # cm
    Width  = 1,
    Zrelease = 0,
    Bypass_elev   = 0,    # cm
    Seepout_Rate  = 0,
    Seepout_Elev  = 0,    # cm
    Seepin_Rate   = 0,
    Seepin_Elev   = 0,    # cm
    Qin_Frac      = 1.0,  # used by dmsta_flow_day_steps; case code may override

    # --- Kinetics builder inputs (STA/PSTA/RES) ---
    # STA
    C1000       = 1000,
    Cstar       = 100,
    Ks_per_yr   = 1.0,
    Z1          = 10,     # cm
    Z2          = 20,     # cm
    Z3          = 30,     # cm
    K2Coef1     = 0.1,
    Chalf       = 50,
    SeasonalFactor = 0,

    # PSTA
    Ytrans      = 1,
    Ysigma      = 1,
    C1000_2     = 1000,
    ks_2        = 1.0,
    zh_2        = 10,     # cm

    # RES
    k_depth_penalty = 0.5,

    # constants builder inputs (optional, but makes behavior explicit)
    Cmax      = 2000,
    C_rain    = 0,
    DryDepo   = 0,
    seepin_conc = 0,
    seepage_c   = 0,
    fseep_recycle = 0,
    fseep_out     = 0
  )

  # ppar <- build_P_kin_slots(
  #   mods = c("STA", "PSTA", "RES"),
  #   registry = NULL,
  #   pparams = params,
  #   Dpy = 365.25,
  #   DutyCycle = params$DutyCycle
  # )
  # validate_P_paramsK(ppar)
  #
  # tanks <- dmsta_build_tanks(params$A_cell, 3)
  # init_P_state <- dmsta_p_init_state(
  #     tanks,
  #     Z_init_m = params$Zinit/100,
  #     C_init_ppb = params$C_init_ppb,
  #     Y_init_mgm2 = params$Y_init_mgm2
  #   )
  # Vint <- params$A_cell * params$Zinit/100
  #
  # constants <- list(
  #   Cmax = if (is.null(params$Cmax)) 2000 else params$Cmax,
  #   C_rain = if (is.null(params$C_rain)) 0 else params$C_rain,
  #   DryDepo = if (is.null(params$DryDepo)) 0 else (params$DryDepo / 365.25),
  #   seepin_conc = params$seepin_conc,
  #   seepout_conc_max = if (!is.null(params$seepage_c)) params$seepage_c else 0,
  #   fseep_recycle = if (!is.null(params$fseep_recycle)) params$fseep_recycle else 0,
  #   fseep_out = if (!is.null(params$fseep_out)) params$fseep_out else 0
  # )

  # out_day1_step <- DMSTAr:::dmsta_flowP_day_steps(
  #   V = Vint,
  #   inputs = series[1,],
  #   params = params,
  #   tanks = tanks,
  #   ppar = ppar,
  #   P_state = init_P_state,
  #   constants = constants
  # )
  # out_day1_step$accum
  # out_day1_step$hyd

  # out_day1 <- DMSTAr:::dmsta_flowP_day(
  #   V = Vint,
  #   inputs = series[1,],
  #   params = params,
  #   tanks = tanks,
  #   ppar = ppar,
  #   P_state = init_P_state,
  #   constants = constants
  # )
  # out_day1$results$P$loads
  # out_day1$results$P$flows
  # out_day1$results$P$conc
  #
  # out_day1$budgets$mass$storage

  # --- Run ---
  out <- dmsta_flowP_series(
    series   = series,
    params   = params,
    ttankS   = 3.0,
    Nsteps   = 4L,
    N_plant  = 3L,
    return_steps = FALSE
  )



  # --- Contract / structure ---
  expect_s3_class(out, "dmsta_result")
  expect_true(is.list(out))
  expect_true(is.data.frame(out$results))
  expect_true(is.list(out$budgets))
  expect_true(is.list(out$meta))

  # results rows match input days
  expect_equal(nrow(out$results), nrow(series))

  # key columns exist (adjust if you rename columns)
  must_have <- c(
    "Date",
    "V_end", "Z_end", "Z_avg", "V_cell_day",
    "Qin", "Qout", "Q_treated", "Q_rel1", "Q_rel2",
    "Cin", "C_out", "C_treated", "C_rel1", "C_rel2"
  )
  expect_true(all(must_have %in% names(out$results)))

  # budgets present and well-formed
  expect_true(is.data.frame(out$budgets$water))
  expect_true(is.data.frame(out$budgets$mass))
  expect_equal(nrow(out$budgets$water), nrow(series))
  expect_equal(nrow(out$budgets$mass), nrow(series))

  # meta should carry key components used by downstream code
  expect_true(is.list(out$meta$tanks))
  expect_true(is.list(out$meta$ppar))
  expect_true(is.list(out$meta$constants))

  # --- basic numeric sanity checks ---
  expect_true(all(is.finite(out$results$V_end)))
  expect_true(all(out$results$V_end >= 0))
  expect_true(all(is.finite(out$results$Qout)))

  # concentrations should be finite (may be 0 if flows are 0)
  expect_true(all(is.finite(out$results$C_out)))

})
