test_that("dmsta_flowP_case runs for a minimal 2-cell network and returns expected structure", {

  # If these helpers are not in the package namespace yet, skip gracefully.
  # (Remove these once everything is in the package and stable.)
  if (!exists("neighbors_zcontrol", mode = "function")) testthat::skip("neighbors_zcontrol() not available")
  if (!exists("fw", mode = "function"))                testthat::skip("fw() not available")
  if (!exists("julian_day", mode = "function"))        testthat::skip("julian_day() not available")

  # --- Boundary forcing series (5 days) ---
  series <- data.frame(
    Date     = as.Date("2000-01-01") + 0:4,
    Qi       = c(1, 1.2, 0.8, 1, 1),     # volume/day
    Ci       = c(50, 60, 55, 50, 52),    # ppb
    Rain     = rep(0, 5),                # m/day
    Et       = rep(0.005, 5),            # m/day
    Zcontrol = rep(1, 5),                # m
    stringsAsFactors = FALSE
  )

  # --- Minimal per-cell params (hydrology + kinetics + init state) ---
  base_params <- list(
    # geometry & timing
    A_cell    = 10,
    DutyCycle = 0.95,

    # initial conditions
    Zinit       = 100,   # cm
    C_init_ppb  = 50,
    Y_init_mgm2 = 0,

    # hydrology parameters used by dmsta_deriv_flow()
    Zmin   = 10,         # cm
    Q_a    = 0,
    Q_b    = 1,
    Zweir  = 100,        # cm
    Q_zmin = 0,          # cm
    Width  = 1,
    Bypass_elev   = 0,   # cm
    Seepout_Rate  = 0,
    Seepout_Elev  = 0,   # cm
    Seepin_Rate   = 0,
    Seepin_Elev   = 0,   # cm
    Qin_Frac      = 1.0, # NOTE: dmsta_flowP_case sets p_run$Qin_Frac = 1.0 internally

    # kinetics builder inputs for STA/PSTA/RES
    # STA
    C1000       = 1000,
    Cstar       = 100,
    Ks_per_yr   = 1.0,
    Z1          = 10,    # cm
    Z2          = 20,    # cm
    Z3          = 30,    # cm
    K2Coef1     = 0.1,
    Chalf       = 50,
    SeasonalFactor = 0,

    # PSTA
    Ytrans      = 1,
    Ysigma      = 1,
    C1000_2     = 1000,
    ks_2        = 1.0,
    zh_2        = 10,    # cm

    # RES
    k_depth_penalty = 0.5,

    # constants inputs (optional; makes defaults explicit)
    Cmax        = 2000,
    C_rain      = 0,
    DryDepo     = 0,
    seepin_conc = 0,
    seepage_c   = 0,
    fseep_recycle = 0,
    fseep_out     = 0
  )

  # Make cell2 slightly different to avoid “accidentally identical” edge behaviors
  params1 <- base_params
  params2 <- modifyList(base_params, list(A_cell = 8))

  # --- Build a minimal 2-cell network ---
  # Cell 1 routes treated outflow to cell 2.
  # Cell 2 is terminal (DownCell = 0 means out of system).
  cells <- list(
    dmsta_make_cell(
      label = "CELL1",
      params = params1,
      ttankS = 3.0,
      DownCell = 2L,
      Qin_Frac = 1.0,
      RecycleIndex = 1L
    ),
    dmsta_make_cell(
      label = "CELL2",
      params = params2,
      ttankS = 3.0,
      DownCell = 0L,
      Qin_Frac = 0.0,
      RecycleIndex = 2L
    )
  )

  # --- Run the case ---
  out <- dmsta_flowP_case(
    series = series,
    cells = cells,
    Nsteps = 4L,
    N_plant = 3L,
    max_iter = 1L,
    conv_tol = 0.01,
    return_cell_series = TRUE,
    keep_Q17 = TRUE
  )

  # --- Contract / structure checks ---
  expect_s3_class(out, "dmsta_network_result")
  expect_true(is.list(out))
  expect_true(is.list(out$results))
  expect_true(is.list(out$budgets))
  expect_true(is.list(out$meta))

  # Results structure: list(case=df, cells=list-of-df)
  expect_true(is.data.frame(out$results$case))
  expect_true(is.list(out$results$cells))
  expect_length(out$results$cells, length(cells))
  expect_true(all(vapply(out$results$cells, is.data.frame, logical(1))))

  # Rows match input days
  expect_equal(nrow(out$results$case), nrow(series))
  expect_true(all(vapply(out$results$cells, nrow, integer(1)) == nrow(series)))

  # Budgets structure
  expect_true(is.list(out$budgets$water))
  expect_true(is.list(out$budgets$mass))

  # Case budgets are data.frames
  expect_true(is.data.frame(out$budgets$water$case))
  expect_true(is.data.frame(out$budgets$mass$case))

  # Cell budgets are lists-of-data.frames (if returned)
  expect_true(is.list(out$budgets$water$cells))
  expect_true(is.list(out$budgets$mass$cells))
  expect_length(out$budgets$water$cells, length(cells))
  expect_length(out$budgets$mass$cells, length(cells))

  # --- Key columns exist in normalized case output ---
  # dmsta_case_components() creates these standardized columns:
  must_have_case <- c(
    "Date",
    "V_end", "Z_avg",
    "RainVol", "EtVol", "NetAtmo",
    "WB_in", "WB_out", "WB_err", "WB_rel",
    "Q_in_total", "L_in_total", "C_in_total",
    "Q_out_total", "L_out_total", "C_out_total"
  )
  expect_true(all(must_have_case %in% names(out$results$case)))

  # Optional Q17 standardized columns if keep_Q17 = TRUE and present
  if ("Q_seep_recycle" %in% names(out$results$case)) {
    expect_true(all(c("Q_seep_recycle","L_seep_recycle","C_seep_recycle") %in% names(out$results$case)))
  }

  # --- Numeric sanity ---
  expect_true(all(is.finite(out$results$case$Q_in_total)))
  expect_true(all(is.finite(out$results$case$Q_out_total)))
  expect_true(all(is.finite(out$results$case$L_in_total)))
  expect_true(all(is.finite(out$results$case$L_out_total)))
  expect_true(all(out$results$case$Q_out_total >= 0))
  expect_true(all(out$results$case$Q_in_total  >= 0))

  # --- Convergence metadata exists ---
  expect_true(is.list(out$meta$convergence))
  expect_true(all(c("converged","iterations_used","history") %in% names(out$meta$convergence)))
})
