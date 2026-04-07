#'
#' DMSTA hydrology engine (internal)
#'
#' Internal helpers implementing DMSTA hydrologic integration, including
#' within-day control depth interpolation, derivative evaluation, and RK4 stepping.
#'
#' @name internal_dmsta_hydro
#' @keywords internal
NULL

#' Evaluate one hydrologic derivative step (DMSTA DerivFlow equivalent)
#'
#' Computes instantaneous hydrologic fluxes and the rate of change of
#' storage volume for a single within-day sub-step, following the
#' DMSTA `DerivFlow` logic implemented in VBA (Module1).
#'
#' This function performs **no time integration**. It evaluates the
#' hydrologic balance at a specific fractional position within the day
#' (`StepFrac`) and returns the instantaneous derivative and
#' associated flux diagnostics. Numerical integration over a day is
#' handled by higher-level integrators (e.g., RK4, Euler).
#'
#' Control-depth interpolation, outflow constraints, bypass logic,
#' seepage, evapotranspiration shutdown, and release semantics follow
#' DMSTA conventions and are governed by the supplied parameters and
#' HydroIndex-style flags.
#'
#' @param V Numeric scalar. Current storage volume (hm^3).
#' @param A_cell Numeric scalar. Cell area (km^2).
#' @param Qi Numeric scalar. Inflow rate (hm^3/day).
#' @param Rain Numeric scalar. Rainfall rate (m/day).
#' @param Et Numeric scalar. Evapotranspiration rate (m/day).
#' @param Zcontrol Numeric scalar. Control depth for the current day (m).
#' @param Zcontrol1 Numeric scalar. Control depth for the previous day (m).
#'   Defaults to `Zcontrol`.
#' @param Zcontrol2 Numeric scalar. Control depth for the next day (m).
#'   Defaults to `Zcontrol`.
#' @param Step Integer. Sub-step index within the day (1-based).
#' @param StepFrac Numeric scalar. Fractional position within the current
#'   sub-step (0–1).
#' @param Nsteps Integer. Total number of sub-steps per day.
#' @param Ddt Numeric scalar. Sub-step duration (days).
#' @param params List of hydrologic parameters. Must include
#'   `interp_option` specifying control-depth interpolation mode.
#'
#' @param Q_zmin Numeric scalar. Minimum depth offset required for outflow (m).
#' @param Zweir Numeric scalar. Weir crest elevation (m).
#' @param Zmin Numeric scalar. Minimum allowable depth (m).
#' @param Vmin Numeric scalar. Minimum allowable volume (hm^3).
#'   If `NULL`, defaults to `Zmin * A_cell`.
#' @param Seepin_Rate Numeric scalar. Seep-in rate coefficient.
#' @param Seepin_Elev Numeric scalar. External groundwater elevation for seep-in (m).
#' @param Seepout_Rate Numeric scalar. Seep-out rate coefficient.
#' @param Seepout_Elev Numeric scalar. External groundwater elevation for seep-out (m).
#' @param RecycleQ Numeric scalar. Internal recycle inflow (hm^3/day).
#' @param Qrelease Numeric scalar. Explicit release flow (hm^3/day).
#' @param Bypass_elev Numeric scalar. Elevation above which full bypass occurs (m).
#' @param Qimax Numeric scalar. Maximum allowable inflow before bypass (hm^3/day).
#' @param Qomax Numeric scalar. Maximum allowable outflow (hm^3/day).
#'   Negative values indicate a minimum bypass constraint.
#' @param Q_a Numeric scalar. Hydraulic discharge coefficient.
#' @param Q_b Numeric scalar. Hydraulic discharge exponent.
#' @param Width Numeric scalar. Effective weir width (km).
#' @param ShutdownET Logical. If `TRUE`, evapotranspiration is reduced
#'   to prevent storage from falling below `Vmin`.
#'
#' @param has_outflow_constraint Logical. Indicates presence of an explicit
#'   outflow constraint (DMSTA HydroIndex(2) > 0).
#' @param has_depth_constraint Logical. Indicates presence of a depth
#'   constraint (DMSTA HydroIndex(1) > 0).
#' @param Qr_0 Numeric scalar. Fixed outflow rate used when
#'   `has_outflow_constraint = TRUE` and
#'   `has_depth_constraint = FALSE`.
#'
#' @return A named list containing instantaneous derivatives and
#' diagnostic quantities:
#' \describe{
#'   \item{Dvdt}{Rate of change of storage volume (hm^3/day).}
#'   \item{Qot}{Total outflow including releases (hm^3/day).}
#'   \item{Qo}{Overflow (depth-dependent) outflow only (hm^3/day).}
#'   \item{Qnet}{Net inflow excluding overflow (hm^3/day).}
#'   \item{Vnext}{Trial volume at end of sub-step before trimming (hm^3).}
#'   \item{Z}{Current depth (m).}
#'   \item{Zcont}{Interpolated and constrained control depth (m).}
#'   \item{Vcontrol}{Control volume corresponding to `Zcont` (hm^3).}
#'   \item{Etest}{Effective evapotranspiration rate after any shutdown (m/day).}
#'   \item{AtmoS}{Net atmospheric flux (rain – ET) as volume (hm^3/day).}
#'   \item{Seepin}{Seep-in flux (hm^3/day).}
#'   \item{Seepout}{Seep-out flux (hm^3/day).}
#'   \item{Bypass}{Bypass flow (hm^3/day).}
#'   \item{Delta}{Fractional position within the day used for interpolation.}
#' }
#'
#' @details
#' This function corresponds closely to the DMSTA VBA routine
#' `DerivFlow` and is intended for internal use by hydrologic
#' integrators. It should not be called directly by users.
#'
#' For strict DMSTA parity, control-depth interpolation behavior is governed
#' by `params$interp_option`, with `2` (mid-day values) as the
#' default.
#'
#'
#' @keywords internal
#' @rdname internal_dmsta_hydro

dmsta_DerivFlow <- function(
    V, A_cell,
    Qi, Rain, Et,
    Zcontrol, Zcontrol1 = Zcontrol, Zcontrol2 = Zcontrol,
    Step, StepFrac, Nsteps, Ddt,
    params,
    # "globals" / parameters that were module-level in VBA:
    Q_zmin = 0, Zweir = 0, Zmin = 0,
    Vmin = NULL,
    Seepin_Rate = 0, Seepin_Elev = 0,
    Seepout_Rate = 0, Seepout_Elev = 0,
    RecycleQ = 0,
    Qrelease = 0,
    Bypass_elev = 0,
    Qimax = 0,
    Qomax = 0,
    Q_a = 0,
    Q_b = 1,
    Width = 1,
    ShutdownET = TRUE,
    # HydroIndex semantics from VBA:
    has_outflow_constraint = FALSE,  # corresponds to HydroIndex(2) > 0
    has_depth_constraint  = FALSE,   # corresponds to HydroIndex(1) > 0
    Qr_0 = 0
) {
  # InterpOption (DMSTA control-depth interpolation mode)
  # InterpOption in VBA: 1,2,3
  InterpOption <- params$interp_option
  if (is.null(InterpOption) || !is.finite(InterpOption)) {
    InterpOption <- 2L  # DMSTA default: mid-day values
  }
  stopifnot(InterpOption %in% c(1L, 2L, 3L))# temporary sanity check

  # State at start of derivative evaluation
  VlasT <- V
  Z_prev <- if (A_cell > 0) VlasT / A_cell else 0
  Z <- if (A_cell > 0) V / A_cell else 0

  # control depth interpolation (VBA: Delta, Select Case InterpOption)
  Delta <- (Step - 1 + StepFrac) / Nsteps

  if (InterpOption == 1L) {
    # simple - assuming inputs are means
    Zcont <- Zcontrol
  } else if (InterpOption == 2L) {
    # interpolate between adjacent days - assuming inputs are mid-day values
    if (Delta <= 0.5) {
      Zcont <- Zcontrol1 + (Zcontrol - Zcontrol1) * (Delta + 0.5)
    } else {
      Zcont <- Zcontrol + (Zcontrol2 - Zcontrol) * (Delta - 0.5)
    }
  } else if (InterpOption == 3L) {
    # interpolate assuming inputs are end-of-day values
    Zcont <- Zcontrol1 * (1 - Delta) + Zcontrol * Delta
  } else {
    # defensive fallback to DMSTA default (2)
    if (Delta <= 0.5) {
      Zcont <- Zcontrol1 + (Zcontrol - Zcontrol1) * (Delta + 0.5)
    } else {
      Zcont <- Zcontrol + (Zcontrol2 - Zcontrol) * (Delta - 0.5)
    }
  }

  # minimum required depth for outflow (excluding releases)
  # Zcont <- max(Zcont + Q_zmin, Zweir, Zmin)
  Zcont <- max(
    (if (has_depth_constraint) Zcont else 0) + Q_zmin,
    Q_zmin,
    Zweir,
    Zmin
  )

  # control volume
  Vcontrol <- Zcont * A_cell

  # Vmin default (VBA had global Vmin; often equals Zmin*A_cell)
  if (is.null(Vmin)) Vmin <- Zmin * A_cell

  # seepage and atmosphere
  Seepin  <- Seepin_Rate  * max(0, Seepin_Elev  - Z) * A_cell
  Seepout <- Seepout_Rate * max(0, Z - Seepout_Elev) * A_cell

  Etest <- Et
  AtmoS <- (Rain - Etest) * A_cell

  # new bypass code (July 2010)
  if (Bypass_elev > 0 && Z_prev > Bypass_elev) {
    # depth constraint bypass
    Bypass <- Qi
  } else if (Qimax > 0 && Qi > Qimax) {
    # max inflow constraint
    Bypass <- Qi - Qimax
  } else if (Qomax < 0) {
    # min bypass constraint added Oct 2007
    Bypass <- min(-Qomax, Qi)
  } else {
    Bypass <- 0
  }

  # net inflow without overflow
  QneT <- Qi + AtmoS - Seepout + Seepin + RecycleQ - Qrelease - Bypass
  Qo <- 0 # overflow = depth-dependent outflow
  VnexT <- VlasT + Ddt * QneT # trial volume at end of time step with no overflow

  # constrain to minimum pool volume with no outflow
  if (VnexT < Vmin) {
    if (Seepout > 0) {
      # first try eliminating outflow seepage
      Seepout <- max(Seepout + (VnexT - Vmin) / Ddt, 0)
      QneT <- Qi + AtmoS - Seepout + Seepin + RecycleQ - Qrelease - Bypass
      VnexT <- VlasT + Ddt * QneT
    }
    # next try reducing ET
    if (VnexT < Vmin && isTRUE(ShutdownET)) {
      Etest <- max(Etest + (VnexT - Vmin) / Ddt / A_cell, 0)
      AtmoS <- (Rain - Etest) * A_cell
      QneT <- Qi + AtmoS - Seepout + Seepin + RecycleQ - Qrelease - Bypass
      VnexT <- VlasT + Ddt * QneT
    }

  } else if (isTRUE(has_outflow_constraint) && !isTRUE(has_depth_constraint)) {
    # VBA: HydroIndex(2) > 0 And HydroIndex(1) = 0  -> Qo = Qr_0
    Qo <- Qr_0

  } else if (Q_a < 0) {
    # constant volume system - not documented
    Qo <- QneT

  } else if (Z > Zcont) {
    # compute overflow
    if (Q_a == 0) {
      # forces outflow computed from control depth
      Qo <- (VlasT - Vcontrol) / Ddt
      if (Qo < 0) Qo <- 0
    } else {
      # hydraulic model
      Qo <- Q_a * Width * (Z - Zweir)^Q_b
    }

    # max outflow constraint
    if (Qomax > 0 && Qo > Qomax) Qo <- Qomax
  }

  # update with overflow and trim if dipping below control volume
  VnexT <- VlasT + (QneT - Qo) * Ddt
  if (VnexT < Vcontrol && Qo > 0) {
    Qo <- max(0, QneT - (Vcontrol - VlasT) / Ddt)
  }

  Dvdt <- QneT - Qo
  Qot  <- Qo + Qrelease

  # Return everything the VBA left in globals (handy for debugging/parity)
  list(
    Dvdt = Dvdt,
    Qot  = Qot,
    Qo   = Qo,
    Qnet = QneT,
    Vnext = VnexT,
    Z = Z,
    Zcont = Zcont,
    Vcontrol = Vcontrol,
    Etest = Etest,
    AtmoS = AtmoS,
    Seepin = Seepin,
    Seepout = Seepout,
    Bypass = Bypass,
    Delta = Delta
  )
}

#' Gate release flows based on storage relative to release depth
#'
#' Applies DMSTA-style release gating logic to determine which release
#' components are active based on current storage volume relative to a
#' specified release depth. When storage is at or below the release
#' elevation, discretionary release components are suppressed.
#'
#' This function implements the release gating behavior used in DMSTA
#' to prevent releases when pool elevation is below the allowable
#' release depth. It is a lightweight helper and does not perform
#' hydrologic integration or time stepping.
#'
#' @param Vo Numeric scalar. Current storage volume (hm^3).
#' @param A_cell Numeric scalar. Cell area (km^2).
#' @param Zrelease Numeric scalar. Release depth threshold (m). Releases
#'   are suppressed when `Vo <= Zrelease * A_cell`.
#' @param Qr_0 Numeric scalar. Fixed outflow or gate-controlled flow
#'   component (hm^3/day). Always passed through if finite.
#' @param Qr_1 Numeric scalar. First discretionary release component
#'   (hm^3/day).
#' @param Qr_2 Numeric scalar. Second discretionary release component
#'   (hm^3/day).
#'
#' @return A named list with elements:
#' \describe{
#'   \item{QrU_0}{Effective fixed/gated outflow component (hm^3/day).}
#'   \item{QrU_1}{Effective first discretionary release component
#'     after gating (hm^3/day).}
#'   \item{QrU_2}{Effective second discretionary release component
#'     after gating (hm^3/day).}
#'   \item{Qrelease}{Total discretionary release
#'     (`QrU_1 + QrU_2`, hm^3/day).}
#'   \item{Sspec}{Total specified outflow including all components
#'     (`QrU_0 + QrU_1 + QrU_2`, hm^3/day).}
#' }
#'
#' @details
#' If `Zrelease` is missing or non-finite, it is treated as zero,
#' allowing releases at all storage levels. When storage is below the
#' release depth, discretionary releases (`Qr_1` and `Qr_2`)
#' are set to zero, while the fixed component (`Qr_0`) remains
#' active if finite.
#'
#' This function corresponds to the release gating behavior embedded
#' in the DMSTA VBA hydrology routines (Module1) and is intended for
#' internal use by hydrology and flow–P orchestration functions.
#'
#'
#' @keywords internal

dmsta_gate_releases <- function(Vo, A_cell, Zrelease, Qr_0, Qr_1, Qr_2) {
  # Zrelease is a depth (same units as Z = V/A_cell). Often meters in your R workflow.
  if (!is.finite(Zrelease) || is.null(Zrelease)) Zrelease <- 0

  if (A_cell > 0 && Vo <= Zrelease * A_cell) {
    QrU_1 <- 0
    QrU_2 <- 0
  } else {
    QrU_1 <- if (is.finite(Qr_1)) Qr_1 else 0
    QrU_2 <- if (is.finite(Qr_2)) Qr_2 else 0
  }
  QrU_0 <- if (is.finite(Qr_0)) Qr_0 else 0
  list(
    QrU_0 = QrU_0,
    QrU_1 = QrU_1,
    QrU_2 = QrU_2,
    Qrelease = QrU_1 + QrU_2,
    Sspec = QrU_0 + QrU_1 + QrU_2
  )
}

#' Perform one RK4 hydrology sub-step (DMSTA DerivFlow kernel)
#'
#' Executes a single Runge–Kutta 4th-order (RK4) integration step for
#' DMSTA hydrology over a sub-day interval. This function evaluates
#' `dmsta_DerivFlow()` at multiple within-step locations and
#' returns the updated volume and diagnostic fluxes.
#'
#' This is a low-level integrator used internally by daily hydrology
#' drivers. It performs no looping over days.
#'
#' @param V Numeric scalar. Volume at the start of the sub-step (hm^3).
#' @param step_index Integer. Sub-step index within the day (1..Nsteps).
#' @param Nsteps Integer. Total number of sub-steps per day.
#' @param Dt Numeric scalar. Sub-step duration (days), typically 1/Nsteps.
#' @param inputs List of hydrologic forcings and control variables
#'   (Qi, Rain, Et, Zcontrol, Zcontrol_prev, Zcontrol_next, releases, recycle).
#' @param params List of hydrologic parameters passed to
#'   `dmsta_DerivFlow()`.
#'
#' @return A list containing updated volume and instantaneous flux
#'   diagnostics for the sub-step.
#'
#' @keywords internal

dmsta_rk4_hydro_step <- function(
    V,                      # start-of-substep volume (Vo in VBA)
    step_index,             # Step in 1..Nsteps
    Nsteps,                 # steps per day
    Dt,                     # 1/Nsteps
    inputs,                 # list: Qi, Rain, Et, Zcontrol, Zcontrol_prev, Zcontrol_next, RecycleQ, Qr0,Qr1,Qr2
    params                 # list: A_cell, Zrelease, plus all dmsta_DerivFlow params
) {
  # runge kutta integration
  A_cell <- params$A_cell
  Vo <- V

  # --- Gate releases using Vo (DMSTA semantics) ---
  gate <- dmsta_gate_releases(
    Vo = Vo, A_cell = A_cell,
    Zrelease = if (!is.null(params$Zrelease)) params$Zrelease else 0,
    Qr_0 = if (!is.null(inputs$Qr0)) inputs$Qr0 else 0,
    Qr_1 = if (!is.null(inputs$Qr1)) inputs$Qr1 else 0,
    Qr_2 = if (!is.null(inputs$Qr2)) inputs$Qr2 else 0
  )
  Qrelease <- gate$Qrelease

  # convenience: call derivative with shared args
  call_deriv <- function(V_stage, StepFrac, Ddt_stage) {
    dmsta_DerivFlow(
      V = V_stage, A_cell = A_cell,
      Qi = inputs$Qi,
      Rain = inputs$Rain,
      Et = inputs$Et,
      Zcontrol = inputs$Zcontrol,
      Zcontrol1 = inputs$Zcontrol_prev,
      Zcontrol2 = inputs$Zcontrol_next,
      Step = step_index,
      StepFrac = StepFrac,
      Nsteps = Nsteps,
      Ddt = Ddt_stage,
      params = params,
      # from inputs
      RecycleQ = if (!is.null(inputs$RecycleQ)) inputs$RecycleQ else 0,
      Qrelease = Qrelease,
      # from params (pass through; dmsta_DerivFlow defines defaults too)
      Q_zmin = params$Q_zmin,
      Zweir  = params$Zweir,
      Zmin   = params$Zmin,
      Vmin   = params$Vmin,
      Seepin_Rate  = params$Seepin_Rate,
      Seepin_Elev  = params$Seepin_Elev,
      Seepout_Rate = params$Seepout_Rate,
      Seepout_Elev = params$Seepout_Elev,
      Bypass_elev  = params$Bypass_elev,
      Qimax = params$Qimax,
      Qomax = params$Qomax,
      Q_a = params$Q_a,
      Q_b = params$Q_b,
      Width = params$Width,
      ShutdownET = params$ShutdownET,
      # HydroIndex-style flags
      has_outflow_constraint = isTRUE(inputs$has_outflow_constraint),
      has_depth_constraint  = isTRUE(inputs$has_depth_constraint),
      Qr_0 = if (!is.null(inputs$Qr0)) inputs$Qr0 else 0
    )
  }

  # --- RK4 stages (matches VBA) ---
  s1 <- call_deriv(Vo, StepFrac = 0.0, Ddt_stage = 0.5 * Dt)
  Dv1 <- s1$Dvdt * Dt
  V2  <- Vo + Dv1 / 2

  s2 <- call_deriv(V2, StepFrac = 0.5, Ddt_stage = 0.5 * Dt)
  Dv2 <- s2$Dvdt * Dt
  V3  <- Vo + Dv2 / 2

  s3 <- call_deriv(V3, StepFrac = 0.5, Ddt_stage = Dt)
  Dv3 <- s3$Dvdt * Dt
  V4  <- Vo + Dv3

  s4 <- call_deriv(V4, StepFrac = 1.0, Ddt_stage = Dt)
  Dv4 <- s4$Dvdt * Dt

  # --- Weighted averages (matches VBA) ---
  V_new <- Vo + (Dv1 + 2 * Dv2 + 2 * Dv3 + Dv4) / 6
  qout_total <- (s1$Qot + 2 * s2$Qot + 2 * s3$Qot + s4$Qot) / 6
  Etest <- (s1$Etest + 2 * s2$Etest + 2 * s3$Etest + s4$Etest) / 6
  Seepout <- (s1$Seepout + 2 * s2$Seepout + 2 * s3$Seepout + s4$Seepout) / 6
  Seepin  <- (s1$Seepin  + 2 * s2$Seepin  + 2 * s3$Seepin  + s4$Seepin ) / 6
  Bypass  <- (s1$Bypass  + 2 * s2$Bypass  + 2 * s3$Bypass  + s4$Bypass ) / 6

  # --- Split total outflow into treated vs releases (VBA Fr_0/Fr_1/Fr_2 logic) ---
  # VBA logic (paraphrased):
  # if qout_total > Sspec -> overflow; fractions use /qout_total
  # else if qout_total > 0 -> constrained releases dominate; fractions use /Sspec
  # else -> 0
  Sspec <- gate$Sspec
  if (qout_total > Sspec && qout_total > 0) {
    Fr_1 <- gate$QrU_1 / qout_total
    Fr_2 <- gate$QrU_2 / qout_total
  } else if (qout_total > 0 && Sspec > 0) {
    Fr_1 <- gate$QrU_1 / Sspec
    Fr_2 <- gate$QrU_2 / Sspec
  } else {
    Fr_1 <- 0
    Fr_2 <- 0
  }
  Fr_0 <- max(0, 1 - Fr_1 - Fr_2)

  Q_rel1   <- qout_total * Fr_1
  Q_rel2   <- qout_total * Fr_2
  Q_treated <- qout_total * Fr_0

  list(
    V_new = V_new,
    qout_total = qout_total,
    Q_treated = Q_treated,
    Q_rel1 = Q_rel1,
    Q_rel2 = Q_rel2,
    Etest = Etest,
    Seepout = Seepout,
    Seepin = Seepin,
    Bypass = Bypass,
    gate = gate
  )
}

#' Integrate DMSTA hydrology over one day using RK4
#'
#' Integrates storage volume and hydrologic fluxes over a single day
#' using fixed-step Runge–Kutta 4th-order integration. This function
#' loops over sub-steps and aggregates daily totals.
#'
#' This is the primary DMSTA-parity hydrology integrator.
#'
#' @param V Numeric scalar. Volume at the start of the day (hm^3).
#' @param inputs List of daily hydrologic forcings.
#' @param params List of hydrologic parameters.
#' @param Nsteps Integer. Number of RK4 sub-steps per day.
#'
#' @return A list containing end-of-day volume, daily flux totals,
#'   and optional per-substep diagnostics.
#'
#' @keywords internal

dmsta_rk4_hydro_day <- function(V, inputs, params, Nsteps = 4L) {
  Dt <- 1 / Nsteps
  V_start <- V

  step_out <- vector("list", Nsteps)

  # integrate within day
  for (k in seq_len(Nsteps)) {
    Vo <- V
    step_res <- dmsta_rk4_hydro_step(
      V = V,
      step_index = k,
      Nsteps = Nsteps,
      Dt = Dt,
      inputs = inputs,
      params = params
    )
    V <- step_res$V_new

    step_out[[k]] <- list(
      step = k, Vo = Vo, V = V, Dt = Dt,
      Qout = step_res$qout_total,
      Q_treated = step_res$Q_treated,
      Q_rel1 = step_res$Q_rel1,
      Q_rel2 = step_res$Q_rel2,
      SeepOut = step_res$Seepout,
      SeepIn  = step_res$Seepin,
      Bypass  = step_res$Bypass,
      Etest = step_res$Etest
    )
  }

  # daily totals by integrating substep-mean rates over Dt
  sum_dt <- sum(vapply(step_out, function(x) x$Dt, 0.0))
  Qout_day <- sum(vapply(step_out, function(x) x$Qout * x$Dt, 0.0))
  Qtrt_day <- sum(vapply(step_out, function(x) x$Q_treated * x$Dt, 0.0))
  Qr1_day  <- sum(vapply(step_out, function(x) x$Q_rel1 * x$Dt, 0.0))
  Qr2_day  <- sum(vapply(step_out, function(x) x$Q_rel2 * x$Dt, 0.0))
  SeepOut_day <- sum(vapply(step_out, function(x) x$SeepOut * x$Dt, 0.0))
  SeepIn_day  <- sum(vapply(step_out, function(x) x$SeepIn  * x$Dt, 0.0))
  Bypass_day  <- sum(vapply(step_out, function(x) x$Bypass  * x$Dt, 0.0))
  EtVol_day   <- sum(vapply(step_out, function(x) x$Etest * params$A_cell * x$Dt, 0.0))
  RainVol_day <- inputs$Rain * params$A_cell
  NetAtmo_day <- RainVol_day - EtVol_day

  # DMSTA-style daily-average depth and time-integrated average volume (optional diagnostics)
  Z_day <- sum(vapply(step_out, function(x) (x$V / params$A_cell) * x$Dt, 0.0)) / max(1e-12, sum_dt)
  V_day <- sum(vapply(step_out, function(x) ((x$Vo + x$V) / 2) * x$Dt, 0.0))

  list(
    V_end = V,          # end-of-day state
    V_day = V_day,      # time-integrated average volume (handy diagnostic)
    Z_day = Z_day,
    Qin = inputs$Qi,
    Qout = Qout_day,
    Q_treated = Qtrt_day,
    Q_rel1 = Qr1_day,
    Q_rel2 = Qr2_day,
    SeepOut = SeepOut_day,
    SeepIn  = SeepIn_day,
    Bypass  = Bypass_day,
    RainVol = RainVol_day,
    EtVol = EtVol_day,
    NetAtmo = NetAtmo_day,
    steps = step_out,
    V_start = V_start
  )
}

#' Integrate DMSTA hydrology over one day using Euler method
#'
#' Performs single-step Euler integration of DMSTA hydrology over
#' one day. This integrator is intended for diagnostic and debugging
#' purposes and does not provide DMSTA numerical parity.
#'
#' @param V Numeric scalar. Volume at the start of the day (hm^3).
#' @param inputs List of daily hydrologic forcings.
#' @param params List of hydrologic parameters.
#'
#' @return A list containing end-of-day volume and daily flux totals.
#'
#'
#' @keywords internal

dmsta_euler_hydro_day <- function(V, inputs, params) {
  # Euler
  A_cell <- params$A_cell

  # Gate releases using Vo for the day (Euler has a single step; Vo=V)
  gate <- dmsta_gate_releases(
    Vo = V, A_cell = A_cell,
    Zrelease = if (!is.null(params$Zrelease)) params$Zrelease else 0,
    Qr_0 = if (!is.null(inputs$Qr0)) inputs$Qr0 else 0,
    Qr_1 = if (!is.null(inputs$Qr1)) inputs$Qr1 else 0,
    Qr_2 = if (!is.null(inputs$Qr2)) inputs$Qr2 else 0
  )

  s <- dmsta_DerivFlow(
    V = V, A_cell = A_cell,
    Qi = inputs$Qi, Rain = inputs$Rain, Et = inputs$Et,
    Zcontrol = inputs$Zcontrol,
    Zcontrol1 = inputs$Zcontrol_prev,
    Zcontrol2 = inputs$Zcontrol_next,
    Step = 1, StepFrac = 0, Nsteps = 1, Ddt = 1,
    params = params,
    RecycleQ = if (!is.null(inputs$RecycleQ)) inputs$RecycleQ else 0,
    Qrelease = gate$Qrelease,
    Q_zmin = params$Q_zmin, Zweir = params$Zweir, Zmin = params$Zmin, Vmin = params$Vmin,
    Seepin_Rate = params$Seepin_Rate, Seepin_Elev = params$Seepin_Elev,
    Seepout_Rate = params$Seepout_Rate, Seepout_Elev = params$Seepout_Elev,
    Bypass_elev = params$Bypass_elev, Qimax = params$Qimax, Qomax = params$Qomax,
    Q_a = params$Q_a, Q_b = params$Q_b, Width = params$Width,
    ShutdownET = params$ShutdownET,
    has_outflow_constraint = isTRUE(inputs$has_outflow_constraint),
    has_depth_constraint  = isTRUE(inputs$has_depth_constraint),
    Qr_0 = if (!is.null(inputs$Qr0)) inputs$Qr0 else 0
  )

  V_new <- V + s$Dvdt * 1
  qout_total <- s$Qot

  # Use the same release fraction logic as RK4 output split
  Sspec <- gate$Sspec
  if (qout_total > Sspec && qout_total > 0) {
    Fr_1 <- gate$QrU_1 / qout_total
    Fr_2 <- gate$QrU_2 / qout_total
  } else if (qout_total > 0 && Sspec > 0) {
    Fr_1 <- gate$QrU_1 / Sspec
    Fr_2 <- gate$QrU_2 / Sspec
  } else {
    Fr_1 <- 0; Fr_2 <- 0
  }
  Fr_0 <- max(0, 1 - Fr_1 - Fr_2)

  list(
    V_end = V_new,
    Qout = qout_total,
    Q_treated = qout_total * Fr_0,
    Q_rel1 = qout_total * Fr_1,
    Q_rel2 = qout_total * Fr_2,
    SeepOut = s$Seepout,
    SeepIn  = s$Seepin,
    Bypass  = s$Bypass,
    EtVol = s$Etest * A_cell,    # m/d * km^2 -> volume/day in your unit system
    RainVol = inputs$Rain * A_cell,
    gate = gate,
    steps = list(list(Vo = V, V = V_new, Dt = 1, Qout = qout_total))
  )
}

# New RKF45 solver
# 0) Utility: trapezoid integration over irregular time grid

#' Trapezoidal integration helper
#'
#' Computes the trapezoidal integral of a time series.
#' Used internally for diagnostic aggregation of sub-step outputs.
#'
#' @param t Numeric vector. Time values.
#' @param y Numeric vector. Values to integrate.
#'
#' @return Numeric scalar. Trapezoidal integral.
#'
#' @keywords internal
.dmsta_trapz_integrate <- function(t, y) {
  n <- length(t)
  if (n < 2) return(0)
  sum(0.5 * (y[-1] + y[-n]) * diff(t))
}

# 1) Adaptive RKF45 (Fehlberg) integrator for y' = f(t, y)
# Internal helper: adaptive Runge–Kutta–Fehlberg (RKF45) integrator
# Solves y' = f(t, y) over [t0, t1]
# NOTE: Internal DMSTAr use only (hydrology ODE support)

#' Adaptive Runge–Kutta–Fehlberg (RKF45) ODE integrator integration diagnostics.
#'
#'
#' Integrates an ordinary differential equation using adaptive
#' Runge–Kutta–Fehlberg (4/5) stepping with local error control.
#'
#' This solver is intended for diagnostic and research use and
#' is not required for DMSTA parity.
#'
#' @param rhs Function computing dy/dt given (t, y).
#' @param y0 Numeric. Initial state.
#' @param t0 Numeric. Start time.
#' @param t1 Numeric. End time.
#' @param h0 Optional initial step size.
#' @param atol Absolute tolerance.
#' @param rtol Relative tolerance.
#' @param hmin Minimum step size.
#' @param hmax Maximum step size.
#' @param max_steps Integer. Maximum number of integration steps.
#'
#' @keywords internal

.dmsta_rkf45_integrate <- function(
    rhs,            # function(t, y) -> dy/dt
    y0,
    t0,
    t1,
    h0 = NULL,
    atol = 1e-8,
    rtol = 1e-6,
    hmin = 1e-10,
    hmax = NULL,
    max_steps = 200000L
) {
  if (!is.finite(t0) || !is.finite(t1)) stop("t0/t1 must be finite.")
  if (t1 <= t0) stop("t1 must be > t0.")

  if (is.null(hmax)) hmax <- (t1 - t0)
  if (is.null(h0))   h0   <- min(hmax, (t1 - t0) / 20)

  # Fehlberg coefficients
  a2 <- 1/4; a3 <- 3/8; a4 <- 12/13; a5 <- 1; a6 <- 1/2

  b21 <- 1/4
  b31 <- 3/32;    b32 <- 9/32
  b41 <- 1932/2197; b42 <- -7200/2197; b43 <- 7296/2197
  b51 <- 439/216;  b52 <- -8; b53 <- 3680/513; b54 <- -845/4104
  b61 <- -8/27;    b62 <- 2;  b63 <- -3544/2565; b64 <- 1859/4104; b65 <- -11/40

  # 4th and 5th order weights
  c4 <- c(25/216, 0, 1408/2565, 2197/4104, -1/5, 0)
  c5 <- c(16/135, 0, 6656/12825, 28561/56430, -9/50, 2/55)

  t <- t0
  y <- y0
  h <- h0
  steps <- 0L

  Tval <- t
  Yval <- y

  while (t < t1 && steps < max_steps) {
    steps <- steps + 1L

    h <- min(max(h, hmin), hmax)
    if (t + h > t1) h <- t1 - t

    k1 <- rhs(t, y)
    k2 <- rhs(t + a2*h, y + h*(b21*k1))
    k3 <- rhs(t + a3*h, y + h*(b31*k1 + b32*k2))
    k4 <- rhs(t + a4*h, y + h*(b41*k1 + b42*k2 + b43*k3))
    k5 <- rhs(t + a5*h, y + h*(b51*k1 + b52*k2 + b53*k3 + b54*k4))
    k6 <- rhs(t + a6*h, y + h*(b61*k1 + b62*k2 + b63*k3 + b64*k4 + b65*k5))

    y4 <- y + h*(c4[1]*k1 + c4[3]*k3 + c4[4]*k4 + c4[5]*k5)
    y5 <- y + h*(c5[1]*k1 + c5[3]*k3 + c5[4]*k4 + c5[5]*k5 + c5[6]*k6)

    err <- abs(y5 - y4)
    tol <- atol + rtol * max(abs(y), abs(y5))

    if (is.finite(err) && err <= tol) {
      t <- t + h
      y <- y5
      Tval <- c(Tval, t)
      Yval <- c(Yval, y)
    }

    if (!is.finite(err) || err == 0) {
      fac <- 2
    } else {
      fac <- 0.9 * (tol / err)^(1/5)
      fac <- max(0.2, min(5, fac))
    }
    h <- h * fac
  }

  if (steps >= max_steps) {
    warning("RKF45 reached max_steps; solution may be incomplete.")
  }

  list(time = Tval, y = Yval)
}

# 2) ODE substep integrator for one DMSTA within-day step k
# Integrates V over tau in [0,1] for the step.

#' Integrate DMSTA hydrology over one day using RKF45
#'
#' Integrates daily hydrology using adaptive Runge–Kutta–Fehlberg
#' stepping. Intended for numerical diagnostics and sensitivity
#' analysis rather than DMSTA parity runs.
#'
#' @param V Numeric scalar. Volume at start of day.
#' @param inputs List of daily hydrologic forcings.
#' @param params List of hydrologic parameters.
#' @param Nsteps Integer. Nominal steps per day.
#' @param atol Absolute tolerance.
#' @param rtol Relative tolerance.
#' @param TrackDeficit Logical. Track volume deficit diagnostics.
#' @param VDeficit0 Numeric. Initial deficit.
#' @param isa_node Logical. Whether cell is a node.
#' @param ... Additional options.
#'
#' @return Hydrology result object with adaptive diagnostics.
#'
#' @keywords internal
.dmsta_ode_one_step <- function(Vo, k, Nsteps, inputs, params,
                                atol = 1e-8, rtol = 1e-6,
                                max_steps = 200000L) {

  A_cell <- params$A_cell
  Dt <- 1 / Nsteps

  # --- Release gating ONCE per step using Vo (DMSTA semantics) ---
  gate <- dmsta_gate_releases(
    Vo = Vo, A_cell = A_cell,
    Zrelease = if (!is.null(params$Zrelease)) params$Zrelease else 0,
    Qr_0 = if (!is.null(inputs$Qr0)) inputs$Qr0 else 0,
    Qr_1 = if (!is.null(inputs$Qr1)) inputs$Qr1 else 0,
    Qr_2 = if (!is.null(inputs$Qr2)) inputs$Qr2 else 0
  )
  Qrelease_fixed <- gate$Qrelease

  # RHS: dV/dt (t in [0,1] is local within-step tau)
  f <- function(tau, V) {
    res <- dmsta_DerivFlow(
      V = V,
      A_cell = A_cell,
      Qi = inputs$Qi,
      Rain = inputs$Rain,
      Et = inputs$Et,
      Zcontrol = inputs$Zcontrol,
      Zcontrol1 = inputs$Zcontrol_prev,
      Zcontrol2 = inputs$Zcontrol_next,
      Step = k,
      StepFrac = tau,
      Nsteps = Nsteps,
      Ddt = Dt,
      params = params,
      RecycleQ = if (!is.null(inputs$RecycleQ)) inputs$RecycleQ else 0,
      Qrelease = Qrelease_fixed,
      Q_zmin = params$Q_zmin,
      Zweir  = params$Zweir,
      Zmin   = params$Zmin,
      Vmin   = params$Vmin,
      Seepin_Rate  = params$Seepin_Rate,
      Seepin_Elev  = params$Seepin_Elev,
      Seepout_Rate = params$Seepout_Rate,
      Seepout_Elev = params$Seepout_Elev,
      Bypass_elev  = params$Bypass_elev,
      Qimax = params$Qimax,
      Qomax = params$Qomax,
      Q_a = params$Q_a,
      Q_b = params$Q_b,
      Width = params$Width,
      ShutdownET = params$ShutdownET,
      has_outflow_constraint = isTRUE(inputs$has_outflow_constraint),
      has_depth_constraint  = isTRUE(inputs$has_depth_constraint),
      Qr_0 = if (!is.null(inputs$Qr0)) inputs$Qr0 else 0
    )
    res$Dvdt
  }

  # Integrate V over local tau in [0, 1]
  sol <- .dmsta_rkf45_integrate(
    rhs = f,
    y0  = Vo,
    t0  = 0,
    t1  = 1,
    atol = atol,
    rtol = rtol,
    max_steps = max_steps
  )

  tau <- sol$time
  Vtraj <- sol$y
  Vend <- tail(Vtraj, 1)

  # Re-evaluate dmsta_DerivFlow on the trajectory for flux averaging
  # (this is how we derive averaged Qot / seep / bypass / Etest for the step)
  eval_at <- function(tau_i, V_i) {
    dmsta_DerivFlow(
      V = V_i,
      A_cell = A_cell,
      Qi = inputs$Qi,
      Rain = inputs$Rain,
      Et = inputs$Et,
      Zcontrol = inputs$Zcontrol,
      Zcontrol1 = inputs$Zcontrol_prev,
      Zcontrol2 = inputs$Zcontrol_next,
      Step = k,
      StepFrac = tau_i,
      Nsteps = Nsteps,
      Ddt = Dt,
      params = params,
      RecycleQ = if (!is.null(inputs$RecycleQ)) inputs$RecycleQ else 0,
      Qrelease = Qrelease_fixed,
      Q_zmin = params$Q_zmin,
      Zweir  = params$Zweir,
      Zmin   = params$Zmin,
      Vmin   = params$Vmin,
      Seepin_Rate  = params$Seepin_Rate,
      Seepin_Elev  = params$Seepin_Elev,
      Seepout_Rate = params$Seepout_Rate,
      Seepout_Elev = params$Seepout_Elev,
      Bypass_elev  = params$Bypass_elev,
      Qimax = params$Qimax,
      Qomax = params$Qomax,
      Q_a = params$Q_a,
      Q_b = params$Q_b,
      Width = params$Width,
      ShutdownET = params$ShutdownET,
      has_outflow_constraint = isTRUE(inputs$has_outflow_constraint),
      has_depth_constraint  = isTRUE(inputs$has_depth_constraint),
      Qr_0 = if (!is.null(inputs$Qr0)) inputs$Qr0 else 0
    )
  }

  # Evaluate at all solver points
  ders <- Map(eval_at, tau, Vtraj)

  Qot_vec     <- vapply(ders, function(x) x$Qot, 0.0)
  Seepout_vec <- vapply(ders, function(x) x$Seepout, 0.0)
  Seepin_vec  <- vapply(ders, function(x) x$Seepin, 0.0)
  Bypass_vec  <- vapply(ders, function(x) x$Bypass, 0.0)
  Etest_vec   <- vapply(ders, function(x) x$Etest, 0.0)

  # Time-average over tau in [0,1]
  # average rate = integral(rate d tau) / 1
  qout_avg    <- .dmsta_trapz_integrate(tau, Qot_vec)
  seepout_avg <- .dmsta_trapz_integrate(tau, Seepout_vec)
  seepin_avg  <- .dmsta_trapz_integrate(tau, Seepin_vec)
  bypass_avg  <- .dmsta_trapz_integrate(tau, Bypass_vec)
  etest_avg   <- .dmsta_trapz_integrate(tau, Etest_vec)

  # Split total outflow into treated vs releases using DMSTA fraction logic based on average qout
  Sspec <- gate$Sspec
  if (qout_avg > Sspec && qout_avg > 0) {
    Fr_1 <- gate$QrU_1 / qout_avg
    Fr_2 <- gate$QrU_2 / qout_avg
  } else if (qout_avg > 0 && Sspec > 0) {
    Fr_1 <- gate$QrU_1 / Sspec
    Fr_2 <- gate$QrU_2 / Sspec
  } else {
    Fr_1 <- 0; Fr_2 <- 0
  }
  Fr_0 <- max(0, 1 - Fr_1 - Fr_2)

  list(
    Vo = Vo,
    V = Vend,
    Dt = Dt,
    Qout = qout_avg,
    Q_treated = qout_avg * Fr_0,
    Q_rel1 = qout_avg * Fr_1,
    Q_rel2 = qout_avg * Fr_2,
    SeepOut = seepout_avg,
    SeepIn  = seepin_avg,
    Bypass  = bypass_avg,
    Etest   = etest_avg,
    gate = gate,
    ode = list(tau = tau, V = Vtraj) # optional: trajectory for debugging
  )
}

# 3) Complete integrator_fun for dmsta_hydro_day(method="ODE")
# Signature matches the hook contract you defined earlier.
dmsta_rkf45_hydro_day <- function(V, inputs, params, Nsteps = 4L,
                                 atol = 1e-8, rtol = 1e-6,
                                 TrackDeficit = FALSE, VDeficit0 = 0,
                                 isa_node = FALSE,
                                 ...) {

  # Optional: Node shortcut (kept minimal; you can expand to match your node-step behavior)
  if (isTRUE(isa_node)) {
    # Mimic the node behavior used in your other drivers (rate-based, no ODE)
    qout <- inputs$Qi + if (!is.null(inputs$RecycleQ)) inputs$RecycleQ else 0
    return(list(
      V_end = V,
      Qin = inputs$Qi,
      Qout = qout,
      Q_treated = qout,
      Q_rel1 = 0, Q_rel2 = 0,
      SeepOut = 0, SeepIn = 0, Bypass = 0,
      RainVol = inputs$Rain * params$A_cell,
      EtVol = 0,
      NetAtmo = inputs$Rain * params$A_cell,
      steps = list(list(step = 1, Vo = V, V = V, Dt = 1, Qout = qout)),
      V_start = V
    ))
  }

  V_start <- V
  step_out <- vector("list", Nsteps)

  VDef <- VDeficit0
  Dt <- 1 / Nsteps

  for (k in seq_len(Nsteps)) {
    one <- .dmsta_ode_one_step(
      Vo = V, k = k, Nsteps = Nsteps,
      inputs = inputs, params = params,
      atol = atol, rtol = rtol
    )
    V <- one$V

    # Dryout handling (optional) — mirrors your step-level placement,
    # but applied using averaged step outflow qout (consistent with your RK4 driver placement)
    if (isTRUE(TrackDeficit)) {
      Vmin_eff <- if (!is.null(params$Vmin) && is.finite(params$Vmin)) params$Vmin else params$Zmin * params$A_cell
      qouT <- one$Qout

      Dm1 <- qouT
      VDef <- VDef + max(0, Vmin_eff - V)

      if (VDef > 0 && qouT > 0) {
        Delta <- qouT - VDef / Dt
        if (Delta < 0) {
          VDef <- VDef - qouT * Dt
          qouT <- 0
        } else {
          VDef <- 0
          qouT <- Delta
        }
        V <- V + (Dm1 - qouT) * Dt
      }

      # overwrite step outflow fields after dryout adjustment
      one$Qout <- qouT
      # recompute treated/release split based on new qout (keeping same gated QrU_*)
      Sspec <- one$gate$Sspec
      if (qouT > Sspec && qouT > 0) {
        Fr_1 <- one$gate$QrU_1 / qouT
        Fr_2 <- one$gate$QrU_2 / qouT
      } else if (qouT > 0 && Sspec > 0) {
        Fr_1 <- one$gate$QrU_1 / Sspec
        Fr_2 <- one$gate$QrU_2 / Sspec
      } else {
        Fr_1 <- 0; Fr_2 <- 0
      }
      Fr_0 <- max(0, 1 - Fr_1 - Fr_2)
      one$Q_treated <- qouT * Fr_0
      one$Q_rel1 <- qouT * Fr_1
      one$Q_rel2 <- qouT * Fr_2
    } else {
      # final dryout clamp (optional)
      Vmin_eff <- if (!is.null(params$Vmin) && is.finite(params$Vmin)) params$Vmin else params$Zmin * params$A_cell
      if (V < Vmin_eff) {
        qouT <- max(0, one$Qout - (Vmin_eff - V) / Dt)
        V <- Vmin_eff
        one$Qout <- qouT
        # update splits
        Sspec <- one$gate$Sspec
        if (qouT > Sspec && qouT > 0) {
          Fr_1 <- one$gate$QrU_1 / qouT
          Fr_2 <- one$gate$QrU_2 / qouT
        } else if (qouT > 0 && Sspec > 0) {
          Fr_1 <- one$gate$QrU_1 / Sspec
          Fr_2 <- one$gate$QrU_2 / Sspec
        } else {
          Fr_1 <- 0; Fr_2 <- 0
        }
        Fr_0 <- max(0, 1 - Fr_1 - Fr_2)
        one$Q_treated <- qouT * Fr_0
        one$Q_rel1 <- qouT * Fr_1
        one$Q_rel2 <- qouT * Fr_2
      }
    }

    step_out[[k]] <- c(list(step = k), one[c("Vo","V","Dt","Qout","Q_treated","Q_rel1","Q_rel2","SeepOut","SeepIn","Bypass","Etest")])
  }

  # Day totals (integrate substep-mean rates over Dt)
  Qout_day <- sum(vapply(step_out, function(x) x$Qout * x$Dt, 0.0))
  Qtrt_day <- sum(vapply(step_out, function(x) x$Q_treated * x$Dt, 0.0))
  Qr1_day  <- sum(vapply(step_out, function(x) x$Q_rel1 * x$Dt, 0.0))
  Qr2_day  <- sum(vapply(step_out, function(x) x$Q_rel2 * x$Dt, 0.0))
  SeepOut_day <- sum(vapply(step_out, function(x) x$SeepOut * x$Dt, 0.0))
  SeepIn_day  <- sum(vapply(step_out, function(x) x$SeepIn  * x$Dt, 0.0))
  Bypass_day  <- sum(vapply(step_out, function(x) x$Bypass  * x$Dt, 0.0))
  EtVol_day   <- sum(vapply(step_out, function(x) x$Etest * params$A_cell * x$Dt, 0.0))
  RainVol_day <- inputs$Rain * params$A_cell
  NetAtmo_day <- RainVol_day - EtVol_day

  # Optional diagnostics similar to your RK4 driver style
  Z_day <- sum(vapply(step_out, function(x) (x$V / params$A_cell) * x$Dt, 0.0)) / max(1e-12, sum(vapply(step_out, function(x) x$Dt, 0.0)))
  V_day <- sum(vapply(step_out, function(x) ((x$Vo + x$V) / 2) * x$Dt, 0.0))

  list(
    V_end = V,
    V_day = V_day,
    Z_day = Z_day,
    Qin = inputs$Qi,
    Qout = Qout_day,
    Q_treated = Qtrt_day,
    Q_rel1 = Qr1_day,
    Q_rel2 = Qr2_day,
    SeepOut = SeepOut_day,
    SeepIn  = SeepIn_day,
    Bypass  = Bypass_day,
    RainVol = RainVol_day,
    EtVol   = EtVol_day,
    NetAtmo = NetAtmo_day,
    steps = step_out,
    V_start = V_start,
    VDeficit_end = if (isTRUE(TrackDeficit)) VDef else NULL,
    method = "RKF45"
  )
}

#' Dispatch daily hydrology integration
#'
#' Selects and executes the requested hydrology integrator
#' (RK4, Euler, RKF45, or custom) for one simulation day.
#'
#' @param V Numeric scalar. Volume at start of day.
#' @param inputs List of daily hydrologic forcings.
#' @param params List of hydrologic parameters.
#' @param method Character. Hydrology integrator method.
#' @param Nsteps Integer. Sub-steps per day.
#' @param integrator_fun Optional custom integrator.
#' @param ... Additional arguments passed to integrator.
#'
#' @return Hydrology result object.
#'
#'
#' @keywords internal

dmsta_hydro_day <- function(
    V,
    inputs,
    params,
    method = c("RK4", "Euler", "RKF45", "custom"),
    Nsteps = 4L,
    integrator_fun = NULL,
    ...
) {
  method <- match.arg(method)

  registry <- list(
    RK4 = function(V, inputs, params, ...) {
      dmsta_rk4_hydro_day(V = V, inputs = inputs, params = params,
                          Nsteps = Nsteps, ...)
    },
    Euler = function(V, inputs, params, ...) {
      dmsta_euler_hydro_day(V = V, inputs = inputs, params = params, ...)
    },
    RKF45 = function(V, inputs, params, ...) {
      # Nsteps used only for DMSTA step semantics; solver is adaptive
      dmsta_rkf45_hydro_day(V = V, inputs = inputs, params = params,
                            Nsteps = Nsteps, ...)
    },
    custom = function(V, inputs, params, ...) {
      if (!is.function(integrator_fun)) {
        stop("method = 'custom' requires a function supplied via `integrator_fun`")
      }
      integrator_fun(V = V, inputs = inputs, params = params,
                     Nsteps = Nsteps, ...)
    }
  )
  ## future implementation consider of LSODA (Livermore Solver
  ### for Ordinary Differential Equations with Automatic method switching)

  registry[[method]](V, inputs, params, ...)
}

#' Determine Whether a Cell Should Be Treated as a Node
#'
#' Internal helper that replicates DMSTA VBA "IsaNode" behavior: if the cell
#' has non-positive area (or area is missing/non-finite), it is treated as a
#' routing node rather than a treatment cell. Node cells bypass storage-based
#' hydrology and treatment kinetics.
#'
#' @param A_cell Numeric scalar. Cell area (km^2). Non-positive values indicate a node.
#' @param IsaNode Logical scalar or NULL. Optional explicit override; if TRUE,
#'   the function returns TRUE regardless of `A_cell`.
#'
#' @return Logical scalar. TRUE if the cell should be handled as a node.
#'
#' @details
#' A cell is classified as a node when:
#' \itemize{
#'   \item `IsaNode` is explicitly `TRUE`, or
#'   \item `A_cell` is NULL, NA, non-finite, or `<= 0`.
#' }
#'
#' @keywords internal
#' @noRd
dmsta_is_node <- function(A_cell, IsaNode = NULL) {
  if (isTRUE(IsaNode)) return(TRUE)
  if (is.null(A_cell) || is.na(A_cell) || !is.finite(A_cell)) return(TRUE)
  A_cell <= 0
}


#' Node Routing Logic for DMSTA (Pass-through with Priority Rules)
#'
#' Internal helper implementing the DMSTA VBA "node" (IsaNode) routing logic.
#' When a cell is treated as a node (no area / no storage), total inflow is
#' passed through and optionally diverted according to DMSTA priority rules:
#' seepage out takes precedence, then full bypass, then low-flow bypass, then
#' outflow cap bypass.
#'
#' @param Qi Numeric scalar. Inflow rate to the node (hm^3/day).
#' @param RecycleQ Numeric scalar. Recycle inflow rate to the node (hm^3/day).
#'   Default is 0.
#' @param Seepout_Rate Numeric scalar. If `> 0`, the node routes *all*
#'   pass-through flow to seepage out (hm^3/day equivalent). Default is 0.
#' @param Qimax Numeric scalar. If `> 0` (and `Seepout_Rate <= 0`),
#'   the node routes *all* pass-through flow to bypass. Default is 0.
#' @param Qomax Numeric scalar. If negative, applies low-flow bypass
#'   (`min(Qi, -Qomax)`). If positive, caps outflow and bypasses remainder.
#'   Default is 0.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{Qout}{Numeric scalar. Routed outflow rate (hm^3/day).}
#'   \item{Etest}{Numeric scalar. Effective ET rate (m/day); always 0 for nodes.}
#'   \item{SeepOut}{Numeric scalar. Seepage outflow rate (hm^3/day).}
#'   \item{SeepIn}{Numeric scalar. Seepage inflow rate (hm^3/day); always 0 for nodes.}
#'   \item{Bypass}{Numeric scalar. Bypass flow rate (hm^3/day).}
#' }
#'
#' @details
#' The node begins with `qout = Qi + RecycleQ`. Then the following priority
#' routing rules are applied:
#' \enumerate{
#'   \item If `Seepout_Rate > 0`: route everything to seepage out.
#'   \item Else if `Qimax > 0`: bypass everything.
#'   \item Else if `Qomax < 0`: low-flow bypass `min(Qi, -Qomax)`.
#'   \item Else if `Qomax > 0`: cap outflow at `Qomax`, bypass the remainder.
#' }
#'
#' @keywords internal
#' @noRd
dmsta_node_route <- function(Qi, RecycleQ = 0,
                             Seepout_Rate = 0,
                             Qimax = 0,
                             Qomax = 0) {

  qout <- max(0, Qi + RecycleQ)   # VBA: qouT = Qi + RecycleQ
  Etest <- 0
  seepout <- 0
  seepin <- 0
  bypass <- 0

  # VBA node priority routing
  if (is.finite(Seepout_Rate) && Seepout_Rate > 0) {
    seepout <- qout
    qout <- 0
  } else if (is.finite(Qimax) && Qimax > 0) {
    bypass <- qout
    qout <- 0
  } else if (is.finite(Qomax) && Qomax < 0) {
    bypass <- min(Qi, -Qomax)
    bypass <- max(0, bypass)
    qout <- max(0, qout - bypass)
  } else if (is.finite(Qomax) && Qomax > 0) {
    bypass <- max(0, qout - Qomax)
    qout <- max(0, qout - bypass)
  }

  list(Qout = qout, Etest = Etest, SeepOut = seepout, SeepIn = seepin, Bypass = bypass)
}

#' Route a DMSTA "node" for one day (internal)
#'
#' Internal helper for DMSTA node-style routing (i.e., a link/node element with
#' no storage). This function calls `dmsta_node_route()` to compute outflow,
#' seepage, bypass, and ET diagnostics for the day, then returns a DMSTAr-style
#' structured result with `results`, `budgets`, and `meta`.
#'
#' Node semantics are consistent with the package parameter convention that
#' `IsaNode` may be derived when `A_cell <= 0`.
#'
#' @param Qi Numeric scalar. External inflow for the day (discharge units used
#'   by the hydrology engine; commonly hm^3/day in DMSTAr examples).
#' @param RecycleQ Numeric scalar. Optional recycle inflow for the day (same
#'   discharge units as `Qi`).
#' @param params Named list of parameters used by `dmsta_node_route()`. This
#'   function expects at least:
#'   \describe{
#'     \item{Seepout_Rate}{Numeric. Seepage outflow rate parameter (node routing).}
#'     \item{Qimax}{Numeric. Inflow capacity parameter (node routing).}
#'     \item{Qomax}{Numeric. Outflow capacity parameter (node routing).}
#'   }
#' @param Nsteps Integer. Number of substep records to emit in `meta$steps`.
#'   For nodes, fluxes are constant across the day; `Nsteps` mainly controls the
#'   structure/length of `meta$steps`.
#' @param Dt Numeric scalar. Substep duration in days. Default is `1.0`.
#'
#' @details
#' The returned `meta$steps` is a list of substep records. Each record follows
#' the internal convention of providing `Vo` (start volume), `V` (end volume),
#' and `Dt` (substep duration).
#'
#' Water budget diagnostics use the standard fields `WB_in`, `WB_out`, `WB_err`,
#' and `WB_rel` used throughout the hydrology engine.
#'
#' @return A structured list with components:
#' \describe{
#'   \item{results}{Named list of daily totals/diagnostics (see below).}
#'   \item{budgets}{Named list with `water` budget diagnostics and `mass = NULL`.}
#'   \item{meta}{Metadata including inputs, parameters used, and `steps`.}
#' }
#'
#' The `results` element contains:
#' \describe{
#'   \item{V_end}{End-of-day volume (node: typically 0).}
#'   \item{Qin}{Inflow discharge for the day (`Qi`).}
#'   \item{Qout}{Total outflow discharge for the day.}
#'   \item{Q_treated}{Outflow assigned to treated component (node: equals `Qout`).}
#'   \item{Q_rel1, Q_rel2}{Release components (node: 0 in this helper).}
#'   \item{SeepOut, SeepIn}{Seepage outflow/inflow for the day.}
#'   \item{Bypass}{Bypass discharge for the day.}
#'   \item{RainVol, EtVol, NetAtmo}{Atmospheric volume terms (node: set to 0 here).}
#' }
#'
#'
#' @keywords internal
#' @noRd
dmsta_node_step <- function(Qi, RecycleQ, params, Nsteps = 1L, Dt = 1.0) {
  nh <- dmsta_node_route(
    Qi = Qi,
    RecycleQ = RecycleQ,
    Seepout_Rate = params$Seepout_Rate,
    Qimax = params$Qimax,
    Qomax = params$Qomax
  )

  # Build a constant step_out (optional, but keeps meta$steps consistent)
  # was vector and loop
  step_out <- replicate(Nsteps, list(
    step = NA_integer_,
    Vo = 0, V = 0,
    Dt = Dt,
    Qi_eff = Qi,
    Qout = nh$Qout,
    Q_treated = nh$Qout,
    Q_rel1 = 0,
    Q_rel2 = 0,
    SeepOut = nh$SeepOut,
    SeepIn  = nh$SeepIn,
    Bypass  = nh$Bypass,
    Etest   = nh$Etest
  ), simplify = FALSE)

  # No storage in node mode
  V_end <- 0
  # No area-based rain/ET volumes in node mode
  RainVol_day <- 0
  EtVol_day <- 0
  NetAtmo_day <- 0

  # Water budget: include recycle on inflow like your normal budget does
  WB_in  <- Qi + RainVol_day + nh$SeepIn + RecycleQ
  WB_out <- nh$Qout + nh$SeepOut + EtVol_day + nh$Bypass
  WB_err <- (V_end - 0) - (WB_in - WB_out)
  WB_rel <- WB_err / max(1e-12, max(WB_in, WB_out))

  list(
    results = list(
      V_end = V_end,
      Qin = Qi,
      Qout = nh$Qout,
      Q_treated = nh$Qout,
      Q_rel1 = 0,
      Q_rel2 = 0,
      SeepOut = nh$SeepOut,
      SeepIn = nh$SeepIn,
      Bypass = nh$Bypass,
      RainVol = RainVol_day,
      EtVol = EtVol_day,
      NetAtmo = NetAtmo_day
    ),
    budgets = list(
      water = list(
        # Core diagnostics
        WB_in = WB_in,
        WB_out = WB_out,
        WB_err = WB_err,
        WB_rel = WB_rel,
        # Atmospheric
        RainVol = RainVol_day,
        EtVol   = EtVol_day,
        NetAtmo = NetAtmo_day,
        # Inflow components
        Qi_eff = Qi, RainVol = RainVol_day, SeepIn = nh$SeepIn, Recycle = RecycleQ,
        # Outflow components
        Qout = nh$Qout, SeepOut = nh$SeepOut, EtVol = EtVol_day, Bypass = nh$Bypass,
        # Storage change
        dV = V_end
      ),
      mass = NULL
    ),
    meta = list(
      V_start = 0,
      V_end   = V_end,
      Nsteps  = Nsteps,
      Dt      = Dt,
      Qi_eff  = Qi,
      params_used = params,
      inputs_used = NULL,  # optionally pass inputs if needed
      steps = step_out,
      IsaNode = TRUE
    )
  )
}

#' Simulate one day of DMSTA hydrology
#'
#' High-level daily hydrology orchestrator applying DMSTA
#' control-depth semantics, release gating, and integrator
#' selection.
#'
#' @param V Numeric scalar. Volume at start of day.
#' @param inputs List of daily inputs.
#' @param params List of hydrologic parameters.
#' @param Qmethod Character. Hydrology integrator.
#' @param Nsteps Integer. Sub-steps per day.
#' @param integrator_fun Optional custom integrator.
#' @param interp_option Control-depth interpolation mode.
#' @param ... Additional arguments.
#'
#' @return Daily hydrology results and diagnostics.
#'
#' @keywords internal
dmsta_flow_day <- function(
    V,
    inputs,
    params,
    Qmethod = c("RK4", "Euler", "RKF45", "custom"),
    Nsteps = 4L,
    integrator_fun = NULL,
    interp_option = 2L,
    ...
) {
  method <- match.arg(Qmethod)

  # Defensive copies
  inputs2 <- as.list(inputs)
  params2 <- as.list(params)

  V_start <- V
  Dt <- 1 / as.integer(Nsteps)

  # Defaults (DMSTAr conventions)
  if (is.null(interp_option)) interp_option <- 2L   # mid‑day
  params2$interp_option <- interp_option

  if (is.null(params2$Qin_Frac)) params2$Qin_Frac      <- 1.0
  if (is.null(inputs2$RecycleQ)) inputs2$RecycleQ      <- 0.0

  # Neighbors for Zcontrol (safe fallback)
  if (is.null(inputs2$Zcontrol_prev)) inputs2$Zcontrol_prev <- inputs2$Zcontrol
  if (is.null(inputs2$Zcontrol_next)) inputs2$Zcontrol_next <- inputs2$Zcontrol

  # Unit normalization (cm → m)
  params2$Zinit        <- cm_to_m(params2$Zinit)
  params2$Zmin         <- cm_to_m(params2$Zmin)
  params2$Zrelease     <- cm_to_m(params2$Zrelease)
  params2$Bypass_elev  <- cm_to_m(params2$Bypass_elev)
  params2$Q_zmin       <- cm_to_m(params2$Q_zmin)
  params2$Zweir        <- cm_to_m(params2$Zweir)
  params2$Seepin_Elev  <- cm_to_m(params2$Seepin_Elev)
  params2$Seepout_Elev <- cm_to_m(params2$Seepout_Elev)

  # Presence-based flags
  # (DMSTAr inputs: 0 = inactive)
  has_Qr0_series  <- any(is.finite(inputs2$Qr0) & inputs2$Qr0 != 0)
  has_Qr1_series  <- any(is.finite(inputs2$Qr1) & inputs2$Qr1 != 0)
  has_Qr2_series  <- any(is.finite(inputs2$Qr2) & inputs2$Qr2 != 0)
  has_depth_series <- any(is.finite(inputs2$Zcontrol) & inputs2$Zcontrol != 0)

  # Normalize inactive series to 0
  if (!has_Qr0_series) inputs2$Qr0 <- 0
  if (!has_Qr1_series) inputs2$Qr1 <- 0
  if (!has_Qr2_series) inputs2$Qr2 <- 0

  # Pass flags downstream
  inputs2$has_outflow_constraint <- isTRUE(has_Qr0_series)
  inputs2$has_depth_constraint   <- isTRUE(has_depth_series)
  params2$force_Q_out            <- isTRUE(has_Qr0_series)

  # Qin routing (Module1e semantics)
  Qi_eff <- params2$Qin_Frac * inputs2$Qi
  inputs2$Qi <- Qi_eff

  # Zrelease clamp to Zmin
  # ONLY when releases are present
  has_rel <- isTRUE(has_Qr1_series) || isTRUE(has_Qr2_series)
  if (isTRUE(has_rel)) {
    if (is.finite(params2$Zrelease) &&
        is.finite(params2$Zmin) &&
        params2$Zrelease < params2$Zmin) {
      params2$Zrelease <- params2$Zmin
    }
  }
  if (!is.finite(params2$Zrelease)) params2$Zrelease <- 0

  # NODE SHORT-CIRCUIT
  isa_node <- dmsta_is_node(
    A_cell  = params2$A_cell,
    IsaNode = params2$IsaNode
  )


  if (isTRUE(isa_node)) {
    return(
      dmsta_node_step(
        Qi       = Qi_eff,
        RecycleQ = inputs2$RecycleQ,
        params   = params2,
        Nsteps   = as.integer(Nsteps),
        Dt       = Dt
      )
    )
  }

  # Run hydrology integrator (storage cell)
  hyd <- dmsta_hydro_day(
    V = V,
    inputs = inputs2,
    params = params2,
    method = method,
    Nsteps = as.integer(Nsteps),
    integrator_fun = integrator_fun,
    ...
  )
  step_out <- hyd$steps

  # Daily totals (integrate substeps)
  sum_steps <- function(name, mult = 1) {
    if (!is.null(hyd[[name]]) && is.finite(hyd[[name]])) return(hyd[[name]])
    if (length(step_out) == 0) return(NA_real_)
    sum(vapply(step_out, function(s) s[[name]] * s$Dt * mult, 0.0))
  }

  Qout_day      <- sum_steps("Qout")
  Q_treated_day <- sum_steps("Q_treated")
  Q_rel1_day    <- sum_steps("Q_rel1")
  Q_rel2_day    <- sum_steps("Q_rel2")
  SeepOut_day   <- sum_steps("SeepOut")
  SeepIn_day    <- sum_steps("SeepIn")
  Bypass_day    <- sum_steps("Bypass")

  V_end <- if (!is.null(hyd$V_end) && is.finite(hyd$V_end)) {
    hyd$V_end
  } else if (length(step_out) == 0) {
    V
  } else if (!is.null(step_out[[length(step_out)]]$V) &&
             is.finite(step_out[[length(step_out)]]$V)) {
    step_out[[length(step_out)]]$V
  } else {
    V
  }

  A_cell <- params2$A_cell

  # Atmospheric volumes
  RainVol_day <- inputs2$Rain * A_cell
  EtVol_day <- if (!is.null(hyd$EtVol)) hyd$EtVol else
    sum(vapply(step_out, function(s) s$Etest * A_cell * s$Dt, 0.0))
  NetAtmo_day <- RainVol_day - EtVol_day

  # Mean depth diagnostic
  Z_day <- sum(vapply(step_out, function(s) (s$V / A_cell) * s$Dt, 0.0)) /
    max(1e-12, sum(vapply(step_out, function(s) s$Dt, 0.0)))

  # Water budget
  RecycleVol_day <- inputs2$RecycleQ

  WB_in  <- Qi_eff + RainVol_day + SeepIn_day + RecycleVol_day
  WB_out <- Qout_day + SeepOut_day + EtVol_day + Bypass_day
  dV <- V_end - V_start
  WB_err <- dV - (WB_in - WB_out)
  WB_rel <- WB_err / max(1e-12, max(WB_in, WB_out))

  # Outputs
  results <- list(
    V_end = V_end,
    Qin = Qi_eff,
    Qout = Qout_day,
    Q_treated = Q_treated_day,
    Q_rel1 = Q_rel1_day,
    Q_rel2 = Q_rel2_day,
    SeepOut = SeepOut_day,
    SeepIn = SeepIn_day,
    Bypass = Bypass_day,
    RainVol = RainVol_day,
    EtVol = EtVol_day,
    NetAtmo = NetAtmo_day,
    Z_end = Z_day
  )

  water_budget <- list(
    # Core diagnostics
    WB_in = WB_in,
    WB_out = WB_out,
    WB_err = WB_err,
    WB_rel = WB_rel,
    # Atmospheric
    RainVol = RainVol_day,
    EtVol   = EtVol_day,
    NetAtmo = NetAtmo_day,
    # Inflow components
    Qi_eff  = Qi_eff,
    RainVol = RainVol_day,
    SeepIn  = SeepIn_day,
    Recycle = RecycleVol_day,
    # Outflow components
    Qout    = Qout_day,
    SeepOut = SeepOut_day,
    EtVol   = EtVol_day,
    Bypass  = Bypass_day,
    # Storage change
    dV = V_end
  )

    meta <- list(
    V_start = V_start,
    V_end   = V_end,
    Nsteps  = as.integer(Nsteps),
    Dt      = Dt,
    Qi_eff  = Qi_eff,
    params_used = params2,
    inputs_used = inputs2,
    steps = step_out,
    method = method,
    flags = list(
      has_Qr0_series = has_Qr0_series,
      has_Qr1_series = has_Qr1_series,
      has_Qr2_series = has_Qr2_series,
      has_depth_series = has_depth_series
    )
  )

  out <- list(
    results = results,
    budgets = list(water = water_budget, mass = NULL),
    meta = meta
  )

  class(out) <- c("dmsta_flow_day_result", "list")
  out
}

#' Return daily hydrology results with normalized sub-steps
#'
#' Wrapper around `dmsta_flow_day()` that guarantees
#' consistent access to hydrology sub-step records for
#' coupling with constituent models.
#'
#' @param V Numeric scalar. Volume at start of day.
#' @param inputs List of daily inputs.
#' @param params List of hydrologic parameters.
#' @param Qmethod Character. Hydrology integrator.
#' @param Nsteps Integer. Sub-steps per day.
#' @param integrator_fun Optional custom integrator.
#' @param interp_option Control-depth interpolation mode.
#' @param ... Additional arguments.
#'
#' @return Hydrology object including `steps`.
#'
#' @keywords internal
dmsta_flow_day_steps <- function(
    V,
    inputs,
    params,
    Qmethod = c("RK4", "Euler", "RKF45", "custom"),
    Nsteps = 4L,
    integrator_fun = NULL,
    interp_option = 2L,
    ...
){
  method <- match.arg(Qmethod)
  # 1) Call the new hydrology engine (authoritative)
  day <- dmsta_flow_day(
    V = V,
    inputs = inputs,
    params = params,
    Qmethod = Qmethod,
    Nsteps = Nsteps,
    integrator_fun = integrator_fun,
    interp_option = interp_option,
    ...
  )


  # 2) Extract substeps safely (RK4 has them, Euler may not)
  steps <- NULL
  if (!is.null(day$meta) && !is.null(day$meta$steps)) {
    steps <- day$meta$steps
  } else if (!is.null(day$steps)) {
    # defensive: some node paths may expose steps here
    steps <- day$steps
  } else {
    steps <- vector("list", 0L)
  }

  # 3) Determine Qi_eff deterministically (do NOT rely on meta presence)
  qi_eff <- if (!is.null(day$meta$Qi_eff)) {
    day$meta$Qi_eff
  } else {
    params$Qin_Frac * inputs$Qi
  }

  # 4) If the integrator returned NO steps, fabricate ONE full-day step
  # This is intentional and correct for:
  #   - Euler
  #   - some RKF45 / custom integrators
  #   - node short-circuit cases
  if (length(steps) == 0L) {
    steps <- list(list(
      step = 1L,
      Dt   = 1.0,
      Vo   = V,
      V    = day$results$V_end,
      Qi_eff    = qi_eff,
      Qout      = day$results$Qout,
      Q_treated = day$results$Q_treated,
      Q_rel1    = day$results$Q_rel1,
      Q_rel2    = day$results$Q_rel2,
      SeepOut   = day$results$SeepOut,
      SeepIn    = day$results$SeepIn,
      Bypass    = day$results$Bypass
    ))
  }


  # 5) Normalize each step so downstream code sees a uniform contract

  steps <- lapply(seq_along(steps), function(k) {
    s <- steps[[k]]

    # step index
    if (is.null(s$step)) s$step <- k

    # timestep (required by P integrators)
    if (is.null(s$Dt)) {
      s$Dt <- 1.0 / as.numeric(Nsteps)
    }

    # normalize volume names
    if (is.null(s$Vo) && !is.null(s$V0))     s$Vo <- s$V0
    if (is.null(s$V)  && !is.null(s$V_end))  s$V  <- s$V_end

    # normalize inflow naming
    if (is.null(s$Qi_eff)) {
      if (!is.null(s$Qin)) s$Qi_eff <- s$Qin else s$Qi_eff <- qi_eff
    }

    # normalize treated / release naming
    if (is.null(s$Q_treated) && !is.null(s$Q_treated_ts))
      s$Q_treated <- s$Q_treated_ts
    if (is.null(s$Q_rel1) && !is.null(s$Q_rel1_ts))
      s$Q_rel1 <- s$Q_rel1_ts
    if (is.null(s$Q_rel2) && !is.null(s$Q_rel2_ts))
      s$Q_rel2 <- s$Q_rel2_ts

    s
  })


  # 6) Build output object WITHOUT mutating `day`
  res <- day$results
  wb  <- if (!is.null(day$budgets)) day$budgets$water else NULL

  out <- list(
    results = res,
    budgets = day$budgets,
    meta    = day$meta,
    steps   = steps
  )

  # 7) Legacy top-level aliases (required by dmsta_flowP_day)
  if (!is.null(res)) {
    out$V_end     <- res$V_end
    out$V_avg     <- res$V_avg
    out$Z_end     <- res$Z_end
    out$Z_avg     <- res$Z_avg

    out$Qin       <- res$Qin
    out$Qout      <- res$Qout
    out$Q_treated <- res$Q_treated
    out$Q_rel1    <- res$Q_rel1
    out$Q_rel2    <- res$Q_rel2
    out$SeepOut   <- res$SeepOut
    out$SeepIn    <- res$SeepIn
    out$Bypass    <- res$Bypass
    out$RainVol   <- res$RainVol
    out$EtVol     <- res$EtVol
    out$NetAtmo   <- res$NetAtmo
  }

  if (!is.null(wb)) {
    out$WB_in  <- wb$WB_in
    out$WB_out <- wb$WB_out
    out$WB_err <- wb$WB_err
    out$WB_rel <- wb$WB_rel

    # defensive legacy exposure
    if (is.null(out$RainVol)) out$RainVol <- wb$RainVol
    if (is.null(out$EtVol))   out$EtVol   <- wb$EtVol
    if (is.null(out$NetAtmo)) out$NetAtmo <- wb$NetAtmo
  }

  # 8) Enforce legacy contract early (fail fast if broken)
  stopifnot(
    is.list(out$steps),
    is.numeric(out$Qout),
    is.numeric(out$SeepOut),
    is.numeric(out$Bypass)
  )

  out
}


#' Run a DMSTA hydrology simulation over a time series
#'
#' Simulates hydrologic behavior for a single cell over
#' a multi-day time series, managing initialization,
#' rolling diagnostics, and result aggregation.
#'
#' @param series Data frame of daily inputs.
#' @param params List of hydrologic parameters.
#' @param V_init Optional initial volume.
#' @param Qmethod Character. Hydrology integrator.
#' @param Nsteps Integer. Sub-steps per day.
#' @param integrator_fun Optional custom integrator.
#' @param interp_option Control-depth interpolation mode.
#' @param ... Additional options.
#'
#' @return Object of class `"dmsta_result"`.
#'
#'
#' @examples
#' \dontrun{
#' # Read data (internal)
#' data(series)
#' series <- series[1:370,]; # for example, limit input file
#'
#' # Data formatting
#' series$Qi <- cfs_to_hm3d(series$Flow) # cfs to hm3/d
#' series$Rain <- in_to_m(series$Rainfall) # inches to meters per day
#' series$Et <- in_to_m(series$ET)
#' series$Zcontrol <- 0/100 # meters; setting to zero to see what happens
#' # If you have release series; otherwise set to 0
#' series$Qr0 <- 0   # constrained outflow (forced Q) if used
#' series$Qr1 <- 0   # release 1
#' series$Qr2 <- 0   # release 2
#' series$Ci <- series$Conc
#'
#' # input parameters
#' params <- list(
#'  A_cell = 2.19, # km2
#'  Zmin   = 2,              # cm
#'  Vmin   = 0,              # hm3
#'  Q_a = 1.0, # qcoef_a; discharge coef
#'  Q_b = 4.0, # qcoef_b; discharge exponent
#'  Zweir = 0, # cm; qcoef_offset; depth offset for outflow computation
#'  Q_zmin = 38,  # cm; qcoef_zmin
#'  Qomax = 0.0, # maximum discharge  hm3/day
#'  Qimax = 0, # maximum inflow  (hm3/day)
#'  Width = 1.55, # km
#'  Bypass_elev = 121.92, # z_byp; mean depth at which bypass begins (m) from input (cm)
#'  Seepout_Rate = 0.00789, # outflow seepage rate per unit head (m/day)/m from input cm/d/cm
#'  Seepout_Elev = 0.0, # elevation controling outflow seepage rate (input cm)
#'  Seepin_Rate  = 0.0,
#'  Seepin_Elev  = 0.0,
#'  ShutdownET = TRUE,
#'  force_Q_out = FALSE,
#'  wrap_interp = TRUE,
#'  Zinit = 40,  # cm; initial water column depth
#'  Qin_Frac = 0.22, # inflow_frac; fraction of basin flows going into this cell
#'  Zrelease = 0,   # cm; z_release; minimum depth for releases
#'  RecycleQ = 0,
#'  IsaNode = NULL,
#'  release_pause_days = 0L,
#'  enable_P_release = FALSE,
#'  K_release = 0,
#'  IsaNode = NULL
#')
#'
#'  V_init <- (cm_to_m(params$Zinit) * params$A_cell)
#'  out <- dmsta_flow_series(V_init, series, params, Nsteps = 4)
#'  out$results
#'  }
#'
#' @export

dmsta_flow_series <- function(
    series,
    params,
    V_init = NULL,
    Qmethod = c("RK4", "Euler", "RKF45", "custom"),
    Nsteps = 4L,
    integrator_fun = NULL,
    interp_option = 2L,
    ...
) {
  Qmethod <- match.arg(Qmethod)

  release_pause_days <- params$release_pause_days
  if (is.null(release_pause_days)) release_pause_days <- 0L

  n <- nrow(series)
  if (n < 1) stop("`series` must have at least one row.")

  # Initialize volume
  if (is.null(V_init)) {
    if (!is.null(params$Zinit) && !is.null(params$A_cell)) {
      V <- (cm_to_m(params$Zinit)) * params$A_cell
    } else {
      stop("Provide V_init or params$Zinit and params$A_cell.")
    }
  } else {
    V <- V_init
  }

  # Pre-allocate outputs
  results_list <- vector("list", n)
  water_budgets <- vector("list", n)

  for (i in seq_len(n)) {
    row <- series[i, , drop = FALSE]
    inputs <- as.list(row)

    # release‑pause control match CEPP PACR A2
    if (i <= release_pause_days) {
      inputs$Qr0 <- 0
      inputs$Qr1 <- 0
      inputs$Qr2 <- 0
    }

    # Neighbor control depths
    inputs$Zcontrol_prev <- if (i > 1) series$Zcontrol[i - 1] else series$Zcontrol[i]
    inputs$Zcontrol_next <- if (i < n) series$Zcontrol[i + 1] else series$Zcontrol[i]

    # Run one day
    day <- dmsta_flow_day(
      V = V,
      inputs = inputs,
      params = params,
      Qmethod = Qmethod,
      Nsteps = Nsteps,
      integrator_fun = integrator_fun,
      interp_option = interp_option,
      ...
    )

    # Save outputs
    results_list[[i]] <- c(list(Date = series$Date[i]), day$results)
    water_budgets[[i]] <- c(list(Date = series$Date[i]), day$budgets$water)

    # Advance state
    V <- day$results$V_end
  }

  results_df <- do.call(rbind.data.frame, results_list)
  water_budget_df <- do.call(rbind.data.frame, water_budgets)

  out <- list(
    results = results_df,
    budgets = list(
      water = water_budget_df,
      mass = NULL
    ),
    meta = list(
      V_init = if (is.null(V_init)) NA_real_ else V_init,
      V_end = V,
      Nsteps = as.integer(Nsteps),
      A_cell = params$A_cell,
      Qmethod = Qmethod,
      release_pause_days = release_pause_days
    )
  )

  class(out) <- c("dmsta_hydro_series_result", "list")
  out
}


