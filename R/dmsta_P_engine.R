
#' DMSTA phosphorus engine (internal)
#'
#' Internal helpers implementing DMSTA coupled phosphorus mass balance within
#' a tank series, including derivative evaluation and RK4 stepping.
#'
#' These functions are used by [dmsta_flowP_series()].
#'
#' @name internal_dmsta_phosphorus
#' @keywords internal
NULL

#' Derivative evaluation for phosphorus mass balance (single tank, single stage)
#'
#' Computes instantaneous derivatives `dMdt` and `dSdt` for one tank at
#' a single RK stage. This corresponds to DMSTA "DerivMass" logic and includes
#' uptake, recycle, sediment/burial, direct sedimentation, and advective transport
#' driven by tank outflow and seepage.
#'
#' @param state Named list containing tank state:
#'   \describe{
#'     \item{M}{Numeric scalar. Water-column mass/state for the tank (kg).}
#'     \item{S}{Numeric scalar. Sediment/areal store for the tank (kg).}
#'   }
#'
#' @param drivers Named list of stage-specific drivers and geometry. Expected fields include:
#'   \describe{
#'     \item{A_tank}{Tank area (same area units as `A_cell`).}
#'     \item{A_cell}{Cell area.}
#'     \item{V_tank_avg}{Average tank volume for the substep (e.g., hm^3).}
#'     \item{Vdel}{Change in tank volume over the substep (same units as `V_tank_avg`).}
#'     \item{StepFrac}{RK stage fraction in \{0, 0.5, 1\}.}
#'     \item{Ddt}{Effective duration used in limiter logic (days).}
#'     \item{Qo_tank}{Tank outflow rate (volume/day).}
#'     \item{Li_tank}{Tank inflow load rate (kg/day).}
#'     \item{Rain}{Rain rate (m/day).}
#'     \item{Seepout}{Cell seep-out flow rate (volume/day).}
#'     \item{Seepin}{Cell seep-in flow rate (volume/day).}
#'     \item{Z_plant}{Rolling/representative plant depth used for reservoir blending (m).}
#'     \item{julian}{Julian day (numeric) used for seasonal modulation.}
#'     \item{Mo_fix}{Reference M for the RK stage limiter (kg).}
#'     \item{So_fix}{Reference S for the RK stage limiter (kg).}
#'   }
#'
#' @param ppar Named list of phosphorus kinetic parameters, typically produced by
#'   `build_P_kin_slots()`. Expected fields include K-length vectors
#'   (`K1`, `K2`, `K3`, `Chalf`, `Z_1`, `Z_2`,
#'   `Z_3`, `K2Coef`) and global scalars (`PModel`, `Ytrans`,
#'   `Ysigma`, `SeasonalFactor`, `CZero`).
#'
#' @param constants Named list of constants and caps used during derivative evaluation.
#'   Expected fields include:
#'   \describe{
#'     \item{Cmax}{Concentration cap used in half-saturation terms.}
#'     \item{C_rain}{Concentration in rainfall (ppb).}
#'     \item{DryDepo}{Dry deposition term (mass/area/day).}
#'     \item{seepin_conc}{Seep-in concentration (ppb).}
#'   }
#'
#' @details
#' Internal diagnostics are computed from stage volume `vV`:
#' \itemize{
#'   \item `C = M / vV` (ppb-like concentration)
#'   \item `Y = S / A_tank` (areal store, mg/m^2 equivalent)
#'   \item `z = vV / A_tank` (m)
#' }
#'
#' Kinetics differ by `ppar$PModel`:
#' \itemize{
#'   \item `PModel == 2L`: transformed kinetics using slot-1 coefficients.
#'   \item otherwise: blended "transition" kinetics combining slots 1--2 and
#'     optionally slot-3 reservoir blending using `Z_plant`.
#' }
#'
#' @return A named list with elements:
#' \describe{
#'   \item{dMdt}{Numeric scalar. Derivative of `M` (kg/day).}
#'   \item{dSdt}{Numeric scalar. Derivative of `S` (kg/day).}
#'   \item{flux}{List of per-area flux terms (rates) used internally:
#'     `P_uptake`, `P_recycle`, `P_sed`, `P_direct`.}
#'   \item{diag}{List of diagnostics: `C`, `Y`, `z`.}
#' }
#'
#' @rdname internal_dmsta_phosphorus
#' @keywords internal

dmsta_DerivMass <- function(state, drivers,ppar,constants) {
  # state: list(M=kg, S=kg)
  # drivers: list with required fields described below
  Ddt <- drivers$Ddt
  stopifnot(is.finite(Ddt), Ddt > 0)


  M <- state$M
  S <- state$S

  A_tank <- drivers$A_tank
  A_cell <- drivers$A_cell
  stopifnot(is.finite(A_tank), A_tank > 0)

  # vV is stage volume inside the step: vV = V_Tank + (StepFrac - 0.5)*Vdel
  vV <- drivers$V_tank_avg + (drivers$StepFrac - 0.5) * drivers$Vdel
  if (!is.finite(vV) || vV <= 0) vV <- 1e-12

  # state-derived
  C <- M / vV          # ppb
  Y <- S / A_tank      # mg/m2 (DMSTA uses kg/km2 == mg/m2)
  zZ <- vV / A_tank    # m

  # depth and concentration modifiers for each module (1..3)
  Fz <- fC <- c(1, 1, 1)
  for (k in 1:3) {
    if (isTRUE(ppar$Z_1[k] > 0 && zZ < ppar$Z_1[k])){
      Fz[k] <- zZ / ppar$Z_1[k] # Fz is a list of 1 so no else Fz[k] <- 1
    }

    if (isTRUE(ppar$Chalf[k] > 0)) {
      fC[k] <- ppar$Chalf[k] / (ppar$Chalf[k] + min(C, constants$Cmax))
    } else {
      fC[k] <- 1
    }
    # CZero adjustment (resistant P)
    if (C < ppar$CZero) {
      fC[k] <- 0
    } else if (C > 0) {
      fC[k] <- fC[k] * (C - ppar$CZero) / C
    } else {
      fC[k] <- 0
    }
  }

  # kinetics
  P_uptake <- P_recycle <- P_sed <- P_direct <- 0

  if (isTRUE(ppar$PModel == 2L)) {
    # PModel 2
    k1A <- ppar$K1[1] * Fz[1] * fC[1]
    k2A <- ppar$K2[1]
    k3A <- ppar$K3[1]
    P_uptake  <- k1A * C
    P_recycle <- k2A * Y
    k3_used <- k3A
  } else {
    # PModel 1: NEWS transition (PSTA <-> SAV)
    Ftrans <- if (ppar$Ysigma <= 0) 1 else 1 / (1 + exp(-(Y - ppar$Ytrans) / ppar$Ysigma))

    k1A <- Ftrans * ppar$K1[1] * Fz[1] * fC[1] + (1 - Ftrans) * ppar$K1[2] * Fz[2] * fC[2]
    k2A <- Ftrans * ppar$K2[1] + (1 - Ftrans) * ppar$K2[2]
    k3A <- Ftrans * ppar$K3[1] + (1 - Ftrans) * ppar$K3[2]

    # reservoir depth penalty blend using Z_plant
    if (isTRUE(ppar$Z_3[1] > 0) &&
        is.finite(drivers$Z_plant) &&
        drivers$Z_plant > ppar$Z_2[1]){

      fres <- if (drivers$Z_plant >= ppar$Z_3[1]) 1
      else (drivers$Z_plant - ppar$Z_2[1]) / (ppar$Z_3[1] - ppar$Z_2[1])

      k1A <- fres * ppar$K1[3] * Fz[3] * fC[3] + (1 - fres) * k1A
      k2A <- fres * ppar$K2[3] + (1 - fres) * k2A
      k3A <- fres * ppar$K3[3] + (1 - fres) * k3A
    }

    P_uptake  <- k1A * C * Y
    P_recycle <- k2A * Y * Y
    k3_used <- k3A
  }

  # SeasonalFactor = [cstar_2] multiplier on uptake
  if (ppar$SeasonalFactor > 0) {
    xJ <- drivers$julian * 2 * pi / 365.25
    sF <- 1.03926 + ppar$SeasonalFactor * (0.02944 * cos(xJ) - 0.14851 * sin(xJ))
    P_uptake <- P_uptake * sF
  }

  # direct sedimentation (slot 1)
  P_direct <- ppar$K2Coef[1] * C * Fz[1]

  # biomass -> soil flux used as burial shortcut in VBA
  P_sed <- k3_used * Y

  # burial
  P_burial <- P_sed + P_direct

  # release (DMSTA for STAs, but in reservoirs with internal recycling this can matter later)
  P_release <- if (isTRUE(ppar$enable_P_release)) ppar$K_release else 0

  # mass fluxes per area in VBA form
  M_outflow <- (drivers$Qo_tank / A_tank + drivers$Seepout / A_cell) * C
  M_inflow  <- drivers$Li_tank / A_tank +
    (constants$DryDepo + drivers$Rain * constants$C_rain) +
    (drivers$Seepin * constants$seepin_conc / A_cell)

  # derivatives (kg/day)
  dMdt <- (M_inflow + P_recycle + P_release - M_outflow - P_uptake - P_direct) * A_tank
  dSdt <- (P_uptake - P_recycle - P_sed) * A_tank

  # limiter (mratiO = 0.01) to prevent collapse below a fraction of Mo/So
  mratiO <- 0.01
  Mo <- drivers$Mo_fix
  So <- drivers$So_fix

  # Safe Mass Limiter
  can_limit_M <-
    is.numeric(Mo) &&
    length(Mo) == 1L &&
    is.finite(Mo) &&
    is.finite(dMdt) &&
    is.finite(Ddt) &&
    Ddt > 0 &&
    Mo > 0

  if (isTRUE(can_limit_M) &&
      isTRUE((Mo + dMdt * Ddt) <= mratiO * Mo)) {

    dMdt <- -Mo * (1 - mratiO) / Ddt

    # S derivative unchanged (matches DMSTA)
    dSdt <- (P_uptake - P_recycle - P_sed) * A_tank
  }

  can_limit_S <-
    is.numeric(So) &&
    length(So) == 1L &&
    is.finite(So) &&
    is.finite(dSdt) &&
    is.finite(Ddt) &&
    Ddt > 0 &&
    So > 0

  if (isTRUE(can_limit_S) &&
      isTRUE((So + dSdt * Ddt) <= mratiO * So)) {

    dSdt <- -So * (1 - mratiO) / Ddt
  }

  list(
    dMdt = dMdt,
    dSdt = dSdt,
    flux = list(
      P_uptake = P_uptake,
      P_recycle = P_recycle,
      P_sed = P_sed,
      P_direct = P_direct,
      P_burial = P_burial,
      P_release = P_release
    ),
    diag = list(C = C, Y = Y, z = zZ)
  )
}

#' RK4 step for phosphorus mass balance in a single tank (internal)
#'
#' Advances tank states `M` and `S` forward by one substep of length
#' `Dt` (days) using a 4th-order Runge-Kutta integrator. Stage-specific
#' driver fields `StepFrac`, `Ddt`, `Mo_fix`, and `So_fix`
#' are injected before calling `dmsta_deriv_mass()`.
#'
#' @param M Numeric scalar. Current water-column state for the tank (kg).
#' @param S Numeric scalar. Current areal/sediment store for the tank (kg).
#' @param args_base Named list with element `drivers`, itself a named list
#'   containing all fields required by `dmsta_deriv_mass()` except
#'   `StepFrac`, `Ddt`, `Mo_fix`, `So_fix`.
#' @param Dt Numeric scalar. Step duration (days), typically `1/Nsteps`.
#' @param ppar Phosphorus kinetic parameters (see `dmsta_deriv_mass()`).
#' @param constants Constants list (see `dmsta_deriv_mass()`).
#' @param clamp Logical; if `TRUE`, clamp outputs to `Mmin` and `Smin`.
#' @param Mmin Numeric scalar. Minimum allowed `M` after clamping.
#' @param Smin Numeric scalar. Minimum allowed `S` after clamping.
#'
#' @return A named list containing updated states, mean derivatives, and
#' averaged diagnostics over the step:
#' \describe{
#'   \item{M_new, S_new}{Updated tank states (kg).}
#'   \item{dMdt_ts, dSdt_ts}{Step-mean derivatives (kg/day).}
#'   \item{P_uptake_ts, P_recycle_ts, P_sed_ts, P_direct_ts}{Step-mean flux terms (rates).}
#'   \item{C_ts, Y_ts, z_ts}{Step-mean diagnostics.}
#' }
#'
#' @rdname internal_dmsta_phosphorus
#' @keywords internal
dmsta_rk4_P_step <- function(M, S, args_base, Dt,ppar,constants,
                             clamp = TRUE, Mmin = 0, Smin = 0) {
  # M, S : current tank states (kg)
  # args_base : list with at least $drivers (a list)
  # Dt : step length (days), typically 1/Nsteps
  #
  # args_base$drivers must contain everything dmsta_deriv_mass expects
  # EXCEPT StepFrac, Ddt, Mo_fix, So_fix (these are set per stage here).

  if (is.null(args_base$drivers) || !is.list(args_base$drivers)) {
    stop("args_base$drivers must be a list.")
  }

  # Fixed reference masses for limiter (DMSTA semantics)
  Mo <- M
  So <- S

  stage_call <- function(Ms, Ss, step_frac, Ddt_use) {
    drv <- args_base$drivers
    drv$StepFrac <- step_frac
    drv$Ddt      <- Ddt_use
    drv$Mo_fix   <- Mo
    drv$So_fix   <- So
    dmsta_DerivMass(state = list(M = Ms, S = Ss),
                     drivers = drv,
                     ppar = ppar,
                     constants = constants
    )
  }

  # RK stage 1
  s1 <- stage_call(Mo, So, step_frac = 0.0, Ddt_use = 0.5 * Dt)
  M2 <- Mo + s1$dMdt * 0.5 * Dt
  S2 <- So + s1$dSdt * 0.5 * Dt

  # RK stage 2
  s2 <- stage_call(M2, S2, step_frac = 0.5, Ddt_use = 0.5 * Dt)
  M3 <- Mo + s2$dMdt * 0.5 * Dt
  S3 <- So + s2$dSdt * 0.5 * Dt

  # RK stage 3
  s3 <- stage_call(M3, S3, step_frac = 0.5, Ddt_use = Dt)
  M4 <- Mo + s3$dMdt * Dt
  S4 <- So + s3$dSdt * Dt

  #  RK stage 4
  s4 <- stage_call(M4, S4, step_frac = 1.0, Ddt_use = Dt)

  #  RK combine: averaged derivatives over step
  dMdt_ts <- (s1$dMdt + 2*s2$dMdt + 2*s3$dMdt + s4$dMdt) / 6
  dSdt_ts <- (s1$dSdt + 2*s2$dSdt + 2*s3$dSdt + s4$dSdt) / 6

  M_new <- Mo + dMdt_ts * Dt
  S_new <- So + dSdt_ts * Dt

  if (clamp) {
    if (isTRUE(!is.finite(M_new)) || isTRUE(M_new < Mmin)){M_new <- Mmin}
    if (isTRUE(!is.finite(S_new)) || isTRUE(S_new < Smin)){S_new <- Smin}
  }

  #  averaged flux diagnostics over step
  avg_flux <- function(name) {
    (s1$flux[[name]] + 2*s2$flux[[name]] + 2*s3$flux[[name]] + s4$flux[[name]]) / 6
  }
  #  averaged diag diagnostics over step (C, Y, z)
  avg_diag <- function(name) {
    (s1$diag[[name]] + 2*s2$diag[[name]] + 2*s3$diag[[name]] + s4$diag[[name]]) / 6
  }

  list(
    M_new = M_new,
    S_new = S_new,
    dMdt_ts = dMdt_ts,
    dSdt_ts = dSdt_ts,
    flux = list(
      P_uptake  = avg_flux("P_uptake"),
      P_recycle = avg_flux("P_recycle"),
      P_sed     = avg_flux("P_sed"),
      P_direct  = avg_flux("P_direct"),
      P_burial  = avg_flux("P_burial"),
      P_release = avg_flux("P_release")
    ),
    diag = list(
      C = avg_diag("C"),
      Y = avg_diag("Y"),
      z = avg_diag("z")
    )
  )

}

#' Integrate one sub-day phosphorus step using Euler method (DMSTA-compatible)
#'
#' Performs a single Euler (forward) integration step for DMSTA phosphorus
#' mass balances over a sub-day interval. This integrator evaluates
#' `dmsta_DerivMass()` once at mid-step (`StepFrac = 0.5`),
#' matching DMSTA conventions for within-day process evaluation.
#'
#' The Euler P integrator is intended primarily for diagnostic,
#' debugging, and unit-testing purposes. For production runs and
#' DMSTA parity comparisons, `dmsta_rk4_P_step()` should be used.
#'
#' @param M Numeric scalar. Current water-column phosphorus mass (kg).
#' @param S Numeric scalar. Current labile (biota/soil) phosphorus mass (kg).
#' @param args_base List containing at least `$drivers`, a list of
#'   hydrologic and geometric drivers required by `dmsta_DerivMass()`.
#'   These drivers are assumed constant over the sub-step.
#' @param Dt Numeric scalar. Sub-step duration (days).
#' @param ppar List of phosphorus model parameters (output of
#'   `build_P_kin_slots()`).
#' @param constants List of physical and chemical constants used by
#'   `dmsta_DerivMass()` (e.g., dry deposition, rain concentration,
#'   seepage concentrations).
#' @param clamp Logical. If `TRUE`, negative or non-finite updated
#'   masses are clamped to `Mmin` and `Smin`. Default is `TRUE`.
#' @param Mmin Numeric scalar. Minimum allowable water-column P mass (kg).
#'   Default is 0.
#' @param Smin Numeric scalar. Minimum allowable labile P mass (kg).
#'   Default is 0.
#'
#' @return A list with elements:
#' \describe{
#'   \item{M_new}{Updated water-column phosphorus mass (kg).}
#'   \item{S_new}{Updated labile phosphorus mass (kg).}
#'   \item{dMdt_ts}{Instantaneous rate of change of `M` (kg/day).}
#'   \item{dSdt_ts}{Instantaneous rate of change of `S` (kg/day).}
#'   \item{flux}{Named list of phosphorus flux rates (per area),
#'     including uptake, recycle, sedimentation, burial, and release.}
#'   \item{diag}{Named list of diagnostic variables (e.g., concentration,
#'     storage density, depth) evaluated at the derivative call.}
#' }
#'
#'
#' @rdname internal_dmsta_phosphorus
#' @keywords internal
dmsta_euler_P_step <- function(
    M, S,
    args_base,
    Dt,
    ppar,
    constants,
    clamp = TRUE,
    Mmin = 0,
    Smin = 0
) {

  stopifnot(is.list(args_base), is.list(args_base$drivers))

  # Fixed reference masses for limiter (DMSTA semantics)
  Mo <- M
  So <- S

  drv <- args_base$drivers
  drv$StepFrac <- 0.5        # mid-step (DMSTA-style)
  drv$Ddt      <- Dt
  drv$Mo_fix   <- Mo
  drv$So_fix   <- So

  # Single derivative evaluation
  s <- dmsta_DerivMass(
    state     = list(M = M, S = S),
    drivers   = drv,
    ppar      = ppar,
    constants = constants
  )

  # Euler update
  M_new <- M + s$dMdt * Dt
  S_new <- S + s$dSdt * Dt

  if (clamp) {
    if (isTRUE(!is.finite(M_new)) || isTRUE(M_new < Mmin)){M_new <- Mmin}
    if (isTRUE(!is.finite(S_new)) || isTRUE(S_new < Smin)){S_new <- Smin}
  }

  list(
    M_new = M_new,
    S_new = S_new,
    dMdt_ts = s$dMdt,
    dSdt_ts = s$dSdt,
    flux = s$flux,
    diag = s$diag
  )
}

#' Integrate one sub-day phosphorus step using a selected integrator adaptive RKF45 or custom solvers)
#' may be added in future development phases.
#'
#' @param method Character string specifying the P integrator to use.
#'   One of `"RK4"` or `"Euler"`.
#' @param M Numeric scalar. Current water-column phosphorus mass (kg).
#' @param S Numeric scalar. Current labile (biota/soil) phosphorus mass (kg).
#' @param args_base List containing at least `$drivers`, a list of
#'   hydrologic and geometric drivers required by `dmsta_DerivMass()`.
#' @param Dt Numeric scalar. Sub-step duration (days).
#' @param ppar List of phosphorus model parameters (output of
#'   `build_P_kin_slots()`).
#' @param constants List of physical and chemical constants used by
#'   `dmsta_DerivMass()`.
#' @param ... Additional arguments passed through to the selected
#'   integrator (e.g., clamping options).
#'
#' @return A list returned by the selected P integrator, containing
#'   updated phosphorus state variables, instantaneous rates, fluxes,
#'   and diagnostic values.
#'
#' @details
#' This function does not implement any numerical integration logic
#' itself. It simply dispatches to the appropriate integrator based on
#' `method`. For strict DMSTA parity and production simulations,
#' `method = "RK4"` is recommended.
#'
#' Dispatches phosphorus mass integration to the requested numerical
#' integrator. This function provides a unified interface for multiple
#' P integrators, mirroring the design of the DMSTAr hydrology engine.
#'
#' Currently supported methods are:
#' \itemize{
#'   \item `"RK4"`: Fourth-order Runge--Kutta integrator
#'     (default, DMSTA parity).
#'   \item `"Euler"`: Single-step Euler integrator
#'     (diagnostic/debugging use).
#' }
#'
#' @rdname internal_dmsta_phosphorus
#' @keywords internal

dmsta_P_step <- function(
    method = c("RK4", "Euler"),
    M, S,
    args_base,
    Dt,
    ppar,
    constants,
    ...
) {
  method <- match.arg(method)

  res <- switch(
    method,
    RK4 = dmsta_rk4_P_step(
      M = M, S = S,
      args_base = args_base,
      Dt = Dt,
      ppar = ppar,
      constants = constants,
      ...
    ),
    Euler = dmsta_euler_P_step(
      M = M, S = S,
      args_base = args_base,
      Dt = Dt,
      ppar = ppar,
      constants = constants,
      ...
    )
  )

  ## ENFORCE RETURN CONTRACT (CRITICAL)

  # If step did not produce updated masses, keep previous state
  if (is.null(res$M_new) || length(res$M_new) != 1L || !is.finite(res$M_new)) {
    res$M_new <- M
  }

  if (is.null(res$S_new) || length(res$S_new) != 1L || !is.finite(res$S_new)) {
    res$S_new <- S
  }

  res
}

#' Integrate hydrology–phosphorus dynamics over one day using sub-day steps
#'
#' Performs coupled hydrology and phosphorus integration over a single day
#' for a storage cell by looping over hydrologic sub-steps and integrating
#' phosphorus dynamics within each sub-step. Hydrology stepping is provided
#' by `dmsta_flow_day_steps()`, and phosphorus integration is delegated
#' to `dmsta_P_step()`.
#'
#' This function assumes a storage cell with one or more tanks in series.
#' Node behavior (cells without storage or P state) is handled upstream
#' in `dmsta_flowP_day()` and should not be routed through this function.
#'
#' @param V Numeric scalar. Volume at the start of the day (hm^3).
#' @param P_state List with elements `M` and `S`, numeric vectors
#'   of length `tanks$Ntanks` giving initial phosphorus state (kg)
#'   in each tank.
#' @param tanks List describing the tanks-in-series configuration
#'   (e.g., output of `dmsta_build_tanks()`).
#' @param inputs List of daily forcing inputs (Qi, Ci, Rain, Et, Zcontrol,
#'   releases, recycle flows, etc.), already prepared for the current day.
#' @param params List of hydrologic parameters for the cell.
#' @param ppar List of phosphorus kinetic parameters (output of
#'   `build_P_kin_slots()`).
#' @param constants List of physical and chemical constants used by
#'   `dmsta_DerivMass()`.
#' @param Qmethod Character string specifying the hydrology integrator
#'   (e.g., `"RK4"`, `"Euler"`, `"RKF45"`, `"custom"`).
#' @param Pmethod Character string specifying the phosphorus integrator
#'   (currently `"RK4"` or `"Euler"`).
#' @param Nsteps Integer. Number of sub-day hydrology steps per day.
#' @param Z_plant Numeric scalar. Rolling-mean depth (m) used for reservoir
#'   penalty blending in phosphorus kinetics.
#' @param integrator_fun Optional custom hydrology integrator function,
#'   used when `Qmethod = "custom"`.
#' @param interp_option Integer. Control-depth interpolation option
#'   (DMSTA semantics; default is mid-day).
#' @param ... Additional arguments passed to hydrology integrators.
#'
#' @return A list with elements:
#' \describe{
#'   \item{hyd}{Normalized hydrology object returned by
#'     `dmsta_flow_day_steps()`.}
#'   \item{P_state_end}{List with updated `M` and `S` vectors
#'     at end of day.}
#'   \item{accum}{Named list of daily accumulated flows, loads, and
#'     phosphorus mechanism totals (uptake, recycle, sedimentation, burial).}
#'   \item{storage}{List of starting and ending total phosphorus storage
#'     for the day.}
#' }
#'
#' @keywords internal

dmsta_flowP_day_steps <- function(
    V,
    P_state,
    tanks,
    inputs,
    params,
    ppar,
    constants,
    Qmethod = c("RK4", "Euler", "RKF45", "custom"),
    Pmethod = c("RK4", "Euler"),
    Nsteps = 4L,
    Z_plant = 0,
    integrator_fun = NULL,
    interp_option = 2L,
    ...
) {

  Qmethod <- match.arg(Qmethod)
  Pmethod <- match.arg(Pmethod)
  # 1) Hydrology with normalized steps (integrator-agnostic)
  hyd <- dmsta_flow_day_steps(
    V = V,
    inputs = inputs,
    params = params,
    Qmethod = Qmethod,
    Nsteps = Nsteps,
    integrator_fun = integrator_fun,
    interp_option = interp_option,
    ...
  )

  steps <- hyd$steps

  # 2) Validate and unpack P state
  Nt <- tanks$Ntanks
  if (length(P_state$M) != Nt || length(P_state$S) != Nt) {
    stop("P_state$M and P_state$S must have length = tanks$Ntanks")
  }

  M <- P_state$M
  S <- P_state$S

  # Storage starts
  M_start <- sum(M)
  S_start <- sum(S)

  # 3) Daily accumulators
  Q_treat <- L_treat <- 0
  Q_r1    <- L_r1    <- 0
  Q_r2    <- L_r2    <- 0
  Q_byp   <- L_byp   <- 0

  # optional seepage bookkeeping (matches DMSTA seepage conc cap pattern)
  Q_seep_rec <- L_seep_rec <- 0
  Q_seep_dis <- L_seep_dis <- 0

  # P budget accumulators for budget
  # treated inflow into tanks (includes recycle)
  Q_in_tanks <- L_in_tanks <- 0

  # split inflow into external flow part vs internal recycle transfer
  Q_in_flow    <- L_in_flow    <- 0
  Q_in_recycle <- L_in_recycle <- 0

  # internal mechanisms
  L_uptake  <- 0
  L_recycle <- 0
  L_sed     <- 0
  L_direct  <- 0
  L_burial  <- 0
  L_release <- 0
  # L_rain <- 0
  # L_drydep <- 0
  # L_seepin <- 0

  # convenience
  jul <- julian_day(inputs$Date)
  A_cell <- params$A_cell

  # 4) Substep loop
  for (hs in steps) {
    Dt <- hs$Dt
    Vo <- hs$Vo
    V1 <- hs$V
    Vavg <- 0.5 * (Vo + V1)

    # Hydrology rates
    Qi_eff   <- hs$Qi_eff
    qout_tot <- hs$Qout
    byp_rate <- hs$Bypass
    seepout  <- hs$SeepOut
    seepin   <- hs$SeepIn

    q_treated <- hs$Q_treated
    q_rel1    <- hs$Q_rel1
    q_rel2    <- hs$Q_rel2

    Ci <- inputs$Ci

    RecycleQ <- if (is.null(inputs$RecycleQ)) 0 else inputs$RecycleQ
    RecycleM <- if (is.null(inputs$RecycleM)) 0 else inputs$RecycleM

    Qin_total <- Qi_eff - byp_rate + RecycleQ
    Qin_total <- max(0, Qin_total)

    # Bypass bookkeeping
    Q_byp <- Q_byp + byp_rate * Dt
    L_byp <- L_byp + (byp_rate * Ci) * Dt

    # P budget: treated inflow into tanks (optional)
    q_in_flow_step <- max(0, Qi_eff - byp_rate)   # excludes recycle
    l_in_flow_step <- q_in_flow_step * Ci         # kg/day

    # totals integrated over day
    Q_in_flow <- Q_in_flow + q_in_flow_step * Dt
    L_in_flow <- L_in_flow + l_in_flow_step * Dt

    Q_in_recycle <- Q_in_recycle + RecycleQ * Dt
    L_in_recycle <- L_in_recycle + RecycleM * Dt

    Q_in_tanks <- Q_in_tanks + Qin_total * Dt
    L_in_tanks <- L_in_tanks + (l_in_flow_step + RecycleM) * Dt

    # Tank routing
    Qo_prev <- 0
    Lo_prev <- 0

    for (tk in seq_len(Nt)) {

      A_tk <- tanks$A_Tank[tk]
      Ftk  <- tanks$F_Tank[tk]
      Fcum <- tanks$Fcum[tk]

      V_tank_avg <- Vavg * Ftk
      Vdel <- (V1 - Vo) * Ftk

      if (tk == 1) {
        Qi_tank <- Qin_total
        Li_tank <- (Qi_eff - byp_rate) * Ci + RecycleM
      } else {
        Qi_tank <- Qo_prev
        Li_tank <- Lo_prev
      }

      Qo_tank <- Qin_total * (1 - Fcum) + qout_tot * Fcum
      Qo_tank <- max(0, Qo_tank)

      # Build drivers for P integrator
      args_base <- list(drivers = list(
        A_tank = A_tk,
        A_cell = A_cell,
        V_tank_avg = V_tank_avg,
        Vdel = Vdel,
        Qo_tank = Qo_tank,
        Li_tank = Li_tank,
        Rain = inputs$Rain,
        Seepout = seepout,
        Seepin = seepin,
        Z_plant = Z_plant,
        julian = jul
      ))

      # Integrate P for this tank and substep
      resP <- dmsta_P_step(
        method = Pmethod,
        M = M[tk],
        S = S[tk],
        args_base = args_base,
        Dt = Dt,
        ppar = ppar,
        constants = constants
      )

      M_old <- M[tk]
      M[tk] <- resP$M_new
      S[tk] <- resP$S_new

      # Mechanism accumulation
      L_uptake  <- L_uptake  + resP$flux$P_uptake  * A_tk * Dt
      L_recycle <- L_recycle + resP$flux$P_recycle * A_tk * Dt
      L_sed     <- L_sed     + resP$flux$P_sed     * A_tk * Dt
      L_direct  <- L_direct  + resP$flux$P_direct * A_tk * Dt
      L_burial  <- L_burial  + resP$flux$P_burial * A_tk * Dt
      L_release <- L_release + resP$flux$P_release * A_tk * Dt

      # Routing concentration
      # Mavg <- 0.5 * (resP$M_new + M[tk])
      Mavg <- 0.5 * (M_old + resP$M_new)
      Cavg <- Mavg / max(V_tank_avg, 1e-12)

      Lo_tank <- Qo_tank * Cavg
      Qo_prev <- Qo_tank
      Lo_prev <- Lo_tank

      if (tk == Nt) {
        Q_treat <- Q_treat + q_treated * Dt
        L_treat <- L_treat + q_treated * Cavg * Dt

        Q_r1 <- Q_r1 + q_rel1 * Dt
        L_r1 <- L_r1 + q_rel1 * Cavg * Dt

        Q_r2 <- Q_r2 + q_rel2 * Dt
        L_r2 <- L_r2 + q_rel2 * Cavg * Dt

        seepC <- Cavg
        if (!is.null(constants$seepout_conc_max) &&
            constants$seepout_conc_max > 0 &&
            seepC > constants$seepout_conc_max) {
          seepC <- constants$seepout_conc_max
        }

        f_rec <- if (is.null(constants$fseep_recycle)) 0 else constants$fseep_recycle
        f_out <- if (is.null(constants$fseep_out)) 0 else constants$fseep_out

        Qsr <- seepout * f_rec
        Qsd <- seepout * f_out

        Q_seep_rec <- Q_seep_rec + Qsr * Dt
        L_seep_rec <- L_seep_rec + Qsr * seepC * Dt

        Q_seep_dis <- Q_seep_dis + Qsd * Dt
        L_seep_dis <- L_seep_dis + Qsd * seepC * Dt
      }
    } # tank loop
  } # substep loop



  # 5) Return step-level outputs (dmsta_flowP_day will finish)
  list(
    hyd = hyd,
    P_state_end = list(M = M, S = S),
    accum = list(
      Q_out = Q_treat + Q_r1 + Q_r2 + Q_byp + Q_seep_dis,
      Q_treated = Q_treat,
      Q_rel1 = Q_r1,
      Q_rel2 = Q_r2,
      Q_bypass = Q_byp,
      Q_seep_recycle = Q_seep_rec,
      Q_seep_discharge = Q_seep_dis,
      Q_in_tanks = Q_in_tanks,
      Q_in_flow = Q_in_flow,
      Q_in_recycle = Q_in_recycle,
      L_out = L_treat + L_r1 + L_r2 + L_byp + L_seep_dis,
      L_treated = L_treat,
      L_rel1 = L_r1,
      L_rel2 = L_r2,
      L_bypass = L_byp,
      L_seep_recycle = L_seep_rec,
      L_seep_discharge = L_seep_dis,
      L_uptake = L_uptake,
      L_recycle = L_recycle,
      L_sed = L_sed,
      L_direct = L_direct,
      L_burial = L_burial,
      L_release = L_release,
      L_in_tanks = L_in_tanks,
      L_in_flow = L_in_flow,
      L_in_recycle = L_in_recycle
    ),
    storage = list(
      M_start = M_start,
      S_start = S_start,
      M_end = sum(M),
      S_end = sum(S)
    )
  )
}

#' Simulate one day of coupled hydrology and phosphorus dynamics
#'
#' Simulates a single day of DMSTA hydrology and phosphorus behavior for
#' either a storage cell or a node. Node behavior (no storage, no phosphorus
#' state) is handled algebraically, while storage cells are integrated using
#' sub-day hydrology steps coupled to phosphorus integration.
#'
#' This function acts as the top-level daily orchestrator for coupled
#' hydrology–phosphorus dynamics. Numerical integration details are
#' delegated to lower-level step functions.
#'
#' @param V Numeric scalar. Volume at the start of the day (hm^3).
#' @param P_state List with elements `M` and `S`, giving initial
#'   phosphorus state (kg). For node cells, this is ignored.
#' @param tanks List describing the tanks-in-series configuration.
#' @param inputs List of daily forcing inputs (Qi, Ci, Rain, Et, Zcontrol,
#'   releases, recycle flows, etc.).
#' @param params List of hydrologic parameters for the cell.
#' @param ppar List of phosphorus kinetic parameters.
#' @param constants List of physical and chemical constants used by
#'   phosphorus derivative calculations.
#' @param Qmethod Character string specifying the hydrology integrator
#'   (`"RK4"`, `"Euler"`, `"RKF45"`, or `"custom"`).
#' @param Pmethod Character string specifying the phosphorus integrator
#'   (`"RK4"` or `"Euler"`).
#' @param Nsteps Integer. Number of hydrology sub-steps per day.
#' @param Z_plant Numeric scalar. Rolling-mean depth (m) used for reservoir
#'   penalty blending.
#' @param integrator_fun Optional custom hydrology integrator function.
#' @param interp_option Integer. Control-depth interpolation option
#'   (DMSTA semantics).
#' @param ... Additional arguments passed to hydrology integrators.
#'
#' @return A list with elements:
#' \describe{
#'   \item{results}{List of daily end-of-day states and aggregated
#'     hydrology and phosphorus outputs.}
#'   \item{budgets}{List containing daily water and phosphorus mass budgets.}
#'   \item{meta}{Metadata describing the simulation day, methods used,
#'     and configuration flags.}
#' }
#'
#' @details
#' For strict DMSTA parity, use `Qmethod = "RK4"` and
#' `Pmethod = "RK4"` with no operational overrides. Optional operational
#' features (e.g., release pauses) are applied upstream in series-level
#' drivers and are not part of the DMSTA core formulation.
#'
#'
#' @keywords internal


dmsta_flowP_day <- function(
    V,
    P_state,
    tanks,
    inputs,
    params,
    ppar,
    constants,
    Qmethod = c("RK4", "Euler", "RKF45", "custom"),
    Pmethod = c("RK4", "Euler"),
    Nsteps = 4L,
    Z_plant = 0,
    integrator_fun = NULL,
    interp_option = 2L,
    ...
) {

  Qmethod <- match.arg(Qmethod)
  Pmethod <- match.arg(Pmethod)

  A_cell <- params$A_cell
  isa_node <- dmsta_is_node(A_cell, params$IsaNode)

  # NODE CASE (algebraic routing only, no P integration)
  if (isTRUE(isa_node)) {

    hyd <- dmsta_flow_day_steps(
      V = 0,
      inputs = inputs,
      params = params,
      Qmethod = Qmethod,
      Nsteps = Nsteps,
      integrator_fun = integrator_fun,
      interp_option = interp_option,
      ...
    )

    Qi_eff <- params$Qin_Frac * inputs$Qi
    Ci <- inputs$Ci

    RecycleQ <- if (is.null(inputs$RecycleQ)) 0 else inputs$RecycleQ
    RecycleM <- if (is.null(inputs$RecycleM)) 0 else inputs$RecycleM

    Qin_total <- Qi_eff - hyd$Bypass + RecycleQ
    Lin_total <- (Qi_eff - hyd$Bypass) * Ci + RecycleM
    Cout <- if (Qin_total > 0) Lin_total / Qin_total else 0

    Q_treat <- hyd$Qout
    L_treat <- Q_treat * Ci

    Q_byp <- hyd$Bypass
    L_byp <- Q_byp * Ci

    seepC <- Cout
    if (!is.null(constants$seepout_conc_max) &&
        constants$seepout_conc_max > 0) {
      seepC <- min(seepC, constants$seepout_conc_max)
    }

    f_rec <- if (is.null(constants$fseep_recycle)) 0 else constants$fseep_recycle
    f_out <- if (is.null(constants$fseep_out)) 0 else constants$fseep_out

    Q_seep_rec <- hyd$SeepOut * f_rec
    Q_seep_dis <- hyd$SeepOut * f_out
    L_seep_rec <- Q_seep_rec * seepC
    L_seep_dis <- Q_seep_dis * seepC

    results <- list(
      V_end = 0,
      Z_end = NA_real_,
      Z_avg = NA_real_,
      V_cell_day = 0,
      Qout = hyd$Qout,
      Q_treated = Q_treat,
      Q_rel1 = 0,
      Q_rel2 = 0,
      SeepOut = hyd$SeepOut,
      SeepIn = 0,
      Bypass = hyd$Bypass,
      P = list(
        flows = list(
          total = Q_treat + 0 + 0 + Q_byp + Q_seep_dis,
          treated = Q_treat, rel1 = 0, rel2 = 0, bypass = Q_byp,
          seep_recycle = Q_seep_rec, seep_discharge = Q_seep_dis
        ),
        loads = list(
          total = L_treat + 0 + 0 + L_byp + L_seep_dis,
          treated = L_treat, rel1 = 0, rel2 = 0, bypass = L_byp,
          seep_recycle = L_seep_rec, seep_discharge = L_seep_dis
        ),
        conc = list(
          # C_out = if(hyd$Qout>0) (L_treat + 0 + 0)/hyd$Qout else NA_real_,
          C_treated = if (Q_treat > 0) L_treat / Q_treat else NA_real_,
          C_rel1    = NA_real_,
          C_rel2    = NA_real_,
          C_out = if (hyd$Qout > 0) L_treat / hyd$Qout else NA_real_,
          C_bypass = if (Q_byp > 0) L_byp / Q_byp else NA_real_,
          C_seep_recycle = if (Q_seep_rec > 0) L_seep_rec / Q_seep_rec else NA_real_,
          C_seep_discharge = if (Q_seep_dis > 0) L_seep_dis / Q_seep_dis else NA_real_
        )
      ),
      P_state_end = list(M = numeric(0), S = numeric(0))
    )

    mass_budget <- list(
      storage = list(
        M_start = 0, S_start = 0,
        M_end = 0,   S_end = 0,
        dM = 0, dS = 0, dP = 0
      ),
      inflow_tanks = list(
        Q_in_tanks = 0,
        L_in_tanks = 0,
        C_in_tanks = NA_real_,
        Q_in_flow  = Qi_eff - hyd$Bypass,
        L_in_flow  = (Qi_eff - hyd$Bypass) * Ci,
        Q_in_recycle = RecycleQ,
        L_in_recycle = RecycleM
      ),
      inputs_external = list(
        L_rain = 0,
        L_drydep = 0,
        L_seepin = 0
      ),
      outputs_external = list(
        L_treated = L_treat,
        L_rel1 = 0,
        L_rel2 = 0,
        L_bypass = L_byp,
        L_seep_discharge = L_seep_dis
      ),
      transfers = list(
        L_seep_recycle_out = L_seep_rec,
        Q_seep_recycle_out = Q_seep_rec
      ),
      mechanisms = list(
        L_uptake = 0, L_recycle = 0, L_sed = 0,
        L_direct = 0, L_burial = 0, L_release = 0
      ),
      closure = list(
        Pin_external = (Qi_eff * Ci),
        Pout_external = (L_treat + L_byp + L_seep_dis),
        Perr_external = 0,
        Prel_external = 0,
        Pin_total = (Qi_eff * Ci),
        Pout_total = (L_treat + L_byp + L_seep_dis),
        Perr_total = 0,
        Prel_total = 0
      )
    )

    return(list(
      results = results,
      budgets = list(
        water = hyd$budgets$water,
        mass  = mass_budget
      ),
      meta = list(Date = inputs$Date, IsaNode = TRUE)
    ))

  }

  # STORAGE CELL CASE (delegated to step integrator)
  step_res <- dmsta_flowP_day_steps(
    V = V,
    P_state = P_state,
    tanks = tanks,
    inputs = inputs,
    params = params,
    ppar = ppar,
    constants = constants,
    Qmethod = Qmethod,
    Pmethod = Pmethod,
    Nsteps = Nsteps,
    Z_plant = Z_plant,
    integrator_fun = integrator_fun,
    interp_option = interp_option,
    ...
  )

  hyd <- step_res$hyd
  acc <- step_res$accum
  stor <- step_res$storage

  # Mass budget (compact but faithful)
  P_start <- stor$M_start + stor$S_start
  P_end <- stor$M_end + stor$S_end
  dM <- stor$M_end - stor$M_start
  dS <- stor$S_end - stor$S_start
  dP <- P_end - P_start

  # external (non-transfer) inputs to tanks are L_in_flow (treated inflow excluding recycle)
  # plus atmospheric + seep-in inputs (apply to the whole cell area)
  L_rain   <- (inputs$Rain * constants$C_rain) * A_cell
  L_drydep <- (constants$DryDepo) * A_cell
  L_seepin <- hyd$SeepIn * constants$seepin_conc

  # define outflows:
  L_out_external <- acc$L_treated + acc$L_rel1 + acc$L_rel2 + acc$L_byp + acc$L_seep_dis
  L_out_transfer <- acc$L_seep_rec  # internal transfer stream
  L_in_external  <- acc$L_in_flow + L_rain + L_drydep + L_seepin
  L_in_transfer  <- acc$L_in_recycle

  # closures:
  # total closure includes transfers
  Pin_total  <- L_in_external + L_in_transfer
  Pout_total <- L_out_external + L_out_transfer
  Perr_total <- dP - (Pin_total - Pout_total)
  Prel_total <- Perr_total / max(1e-12, max(Pin_total, Pout_total))

  # external-only closure (will not close if transfers exist and you don't track transit)
  Perr_external <- dP - (L_in_external - L_out_external)
  Prel_external <- Perr_external / max(1e-12, max(L_in_external, L_out_external))


  # Mean depth diagnostics
  A_cell <- params$A_cell
  V_cell_day <- sum(vapply(hyd$steps, function(s) {
    0.5 * (s$Vo + s$V) * s$Dt
  }, 0.0))

  Z_avg <- if (A_cell > 0) V_cell_day / A_cell else NA_real_
  Z_end <- if (A_cell > 0) hyd$V_end / A_cell else NA_real_

  # Results assembly
  results <- list(
    V_end      = hyd$V_end,
    Z_end      = Z_end,
    Z_avg      = Z_avg,
    V_cell_day = V_cell_day,
    Qin        = hyd$Qin,
    Qout       = hyd$Qout,
    Q_treated  = acc$Q_treated,
    Q_rel1     = acc$Q_rel1,
    Q_rel2     = acc$Q_rel2,
    SeepOut    = hyd$SeepOut,
    SeepIn     = hyd$SeepIn,
    Bypass     = hyd$Bypass,
    P = list(
      flows = list(
        treated = acc$Q_treated,
        rel1 = acc$Q_rel1,
        rel2 = acc$Q_rel2,
        bypass = acc$Q_bypass,
        seep_recycle = acc$Q_seep_recycle,
        seep_discharge = acc$Q_seep_discharge
      ),
      loads = list(
        treated = acc$L_treated,
        rel1 = acc$L_rel1,
        rel2 = acc$L_rel2,
        bypass = acc$L_bypass,
        seep_recycle = acc$L_seep_recycle,
        seep_discharge = acc$L_seep_discharge
      ),
      conc = list(
        C_treated = fw(acc$L_treated, acc$Q_treated),
        C_rel1    = fw(acc$L_rel1, acc$Q_rel1),
        C_rel2    = fw(acc$L_rel2, acc$Q_rel2),
        C_out     = fw(acc$L_out,acc$Q_out),
        C_bypass  = fw(acc$L_byp, acc$Q_byp),
        C_seep_recycle   = fw(acc$L_seep_rec, acc$Q_seep_rec),
        C_seep_discharge = fw(acc$L_seep_dis, acc$Q_seep_dis)
      )
    ),
    P_state_end = step_res$P_state_end
  )

  # Water budget passed through
  water_budget <- hyd$budgets$water

  mass_budget <- list(
    storage = c(stor, P_start = P_start, P_end = P_end,
                dM = dM,dS = dS,dP = dP),
    inflow_tanks = list(
      Q_in_tanks = acc$Q_in_tanks,
      L_in_tanks = acc$L_in_tanks,
      C_in_tanks = fw(acc$L_in_tanks, acc$Q_in_tanks),
      Q_in_flow  = acc$Q_in_flow,
      L_in_flow  = acc$L_in_flow,
      C_in_flow  = fw(acc$L_in_flow, acc$Q_in_flow),
      Q_in_recycle = acc$Q_in_recycle,
      L_in_recycle = acc$L_in_recycle
    ),
    inputs_external = list(
      L_rain = L_rain, L_drydep = L_drydep, L_seepin = L_seepin
    ),
    outputs_external = list(
      L_treated = acc$L_treated, L_rel1 = acc$L_r1, L_rel2 = acc$L_r2,
      L_bypass = acc$L_byp, L_seep_discharge = acc$L_seep_dis
    ),
    transfers = list(
      L_seep_recycle_out = acc$L_seep_rec,
      Q_seep_recycle_out = acc$Q_seep_rec
    ),
    mechanisms = list(
      L_uptake = acc$L_uptake,
      L_recycle = acc$L_recycle,
      L_sed = acc$L_sed,
      L_direct = acc$L_direct,
      L_burial = acc$L_burial,
      L_release = acc$L_release
    ),
    closure = list(
      Pin_total = Pin_total, Pout_total = Pout_total,
      Perr_total = Perr_total, Prel_total = Prel_total,
      Pin_external = L_in_external, Pout_external = L_out_external,
      Perr_external = Perr_external, Prel_external = Prel_external
    )
  )

  list(
    results = results,
    budgets = list(
      water = water_budget,
      mass  = mass_budget
    ),
    meta = list(
      Date = inputs$Date,
      V_start = V,
      V_end = hyd$V_end,
      A_cell = A_cell,
      Nsteps = Nsteps,
      Qmethod = Qmethod,
      Pmethod = Pmethod,
      Z_plant = Z_plant
    )
  )
}


#' Run a coupled hydrology–phosphorus simulation over a time series
#'
#' Simulates DMSTA hydrology and phosphorus dynamics over a multi-day
#' time series for a single cell. This function manages time-series
#' iteration, initialization, rolling diagnostics (e.g., `Z_plant`),
#' and aggregation of daily results, while delegating daily integration
#' to `dmsta_flowP_day()`.
#'
#' Series-level operational semantics (e.g., HydroIndex-style presence
#' flags, optional release warm-up periods) are applied here and passed
#' downstream as daily inputs.
#'
#' @param series Data frame containing daily forcing inputs. Must include
#'   at least `Date`, `Qi`, `Ci`, `Rain`, `Et`,
#'   and `Zcontrol`. Optional columns include release and recycle terms.
#' @param params List of hydrologic and phosphorus parameters.
#' @param pparams Optional list of phosphorus parameters to merge into
#'   `params`.
#' @param ttankS Numeric. Number of tanks in series (may be fractional).
#' @param Nsteps Integer. Number of hydrology sub-steps per day.
#' @param N_plant Integer. Window length (days) used to compute rolling
#'   mean depth for `Z_plant`.
#' @param Qmethod Character string specifying the hydrology integrator.
#' @param Pmethod Character string specifying the phosphorus integrator.
#' @param integrator_fun Optional custom hydrology integrator function.
#' @param interp_option Integer. Control-depth interpolation option.
#' @param ppar Optional precomputed phosphorus kinetic parameter list.
#' @param constants Optional list of constants used by phosphorus derivatives.
#' @param tanks Optional pre-built tanks-in-series configuration.
#' @param V_init Optional initial volume (hm^3). If `NULL`, derived
#'   from depth initialization rules.
#' @param init_P_state Optional initial phosphorus state. If `NULL`,
#'   initialized from depth and concentration parameters.
#' @param return_steps Logical. If `TRUE`, store daily hydrology
#'   sub-step outputs in the result metadata.
#' @param ... Additional arguments passed to daily hydrology integration.
#'
#' @return An object of class `"dmsta_result"` with elements:
#' \describe{
#'   \item{results}{Data frame of daily hydrology and phosphorus outputs.}
#'   \item{budgets}{List containing water and phosphorus mass budget
#'     data frames.}
#'   \item{meta}{Metadata describing model configuration, initialization,
#'     and methods used.}
#' }
#'
#'
#' @export

dmsta_flowP_series <- function(
    series,
    params,
    pparams = NULL,
    ttankS = 3.0,
    Nsteps = 4L,
    N_plant = 30L,
    Qmethod = c("RK4", "Euler", "RKF45", "custom"),
    Pmethod = c("RK4", "Euler"),
    integrator_fun = NULL,
    interp_option = 2L,
    ppar = NULL,
    constants = NULL,
    tanks = NULL,
    V_init = NULL,
    init_P_state = NULL,
    return_steps = FALSE,
    ...
) {

  Qmethod <- match.arg(Qmethod)
  Pmethod <- match.arg(Pmethod)

  # Basic checks
  if (!is.data.frame(series)) stop("'series' must be a data.frame.")
  req <- c("Date", "Qi", "Ci", "Rain", "Et", "Zcontrol")
  miss <- setdiff(req, names(series))
  if (length(miss) > 0) stop("series is missing columns: ", paste(miss, collapse = ", "))

  n <- nrow(series)
  if (n < 1) stop("series has zero rows.")

  # Merge params + pparams if provided
  # carry over from older dev versions
  if (!is.null(pparams)) params <- modifyList(params, pparams)

  # Build tanks once
  if (is.null(tanks)) {
    tanks <- dmsta_build_tanks(params$A_cell, ttankS)
  }

  # Build P kinetics slots once
  if (is.null(ppar)) {
    ppar <- build_P_kin_slots(
      mods = c("STA", "PSTA", "RES"),
      registry = NULL,
      pparams = params,
      Dpy = 365.25,
      DutyCycle = params$DutyCycle
    )
    validate_P_paramsK(ppar)
  }

  # Constants (as expected by dmsta_DerivMass)
  if (is.null(constants)) {

    seepin_conc_val <- 0
    if (!is.null(params$Seepin_Conc)) seepin_conc_val <- params$Seepin_Conc
    if (!is.null(params$seepin_conc)) seepin_conc_val <- params$seepin_conc

    constants <- list(
      Cmax = if (is.null(params$Cmax)) 2000 else params$Cmax,
      C_rain = if (is.null(params$C_rain)) 0 else params$C_rain,
      DryDepo = if (is.null(params$DryDepo)) 0 else (params$DryDepo / 365.25),
      seepin_conc = seepin_conc_val,
      seepout_conc_max = if (!is.null(params$seepage_c)) params$seepage_c else 0,
      fseep_recycle = if (!is.null(params$fseep_recycle)) params$fseep_recycle else 0,
      fseep_out = if (!is.null(params$fseep_out)) params$fseep_out else 0
    )
  }

  # SERIES-LEVEL HydroIndex-style flags (DMSTA semantics)
  has_Qr0_series <- if (!is.null(series$has_Qr0_series)) {
    isTRUE(series$has_Qr0_series)
  } else {
    any(is.finite(series$Qr0) & series$Qr0 != 0)
  }

  has_Qr1_series <- any(is.finite(series$Qr1) & series$Qr1 != 0)
  has_Qr2_series <- any(is.finite(series$Qr2) & series$Qr2 != 0)

  has_depth_constraint_series <- if (!is.null(series$has_depth_constraint)) {
    isTRUE(series$has_depth_constraint)
  } else {
    any(is.finite(series$Zcontrol) & series$Zcontrol != 0)
  }

  # Initial depth and volume (DMSTA-faithful)
  Zmin_m <- params$Zmin / 100

  if (isTRUE(has_depth_constraint_series)) {
    zc1 <- series$Zcontrol[1]
    if (!is.finite(zc1)) zc1 <- 0
    Z0_m <- max(zc1, Zmin_m)
  } else {
    Z0_m <- max(params$Zinit / 100, Zmin_m)
  }

  if (is.null(V_init)) {
    V_init <- params$A_cell * Z0_m
  } else {
    Z0_m <- max(V_init / params$A_cell, Zmin_m)
  }

  V <- V_init

  # Initial P state
  if (is.null(init_P_state)) {
    init_P_state <- dmsta_p_init_state(
      tanks,
      Z_init_m = Z0_m,
      C_init_ppb = params$C_init_ppb,
      Y_init_mgm2 = params$Y_init_mgm2
    )
  }
  P_state <- init_P_state

  # Rolling depth history for Z_plant
  Z_hist <- numeric(n)
  Z_hist[1] <- if (params$A_cell > 0) V / params$A_cell else NA_real_

  # Output containers
  results_df <- data.frame(
    Date = as.Date(series$Date),
    V_end = NA_real_, Z_end = NA_real_, Z_avg = NA_real_, V_cell_day = NA_real_,
    Qin = NA_real_,
    Qout = NA_real_, Q_treated = NA_real_, Q_rel1 = NA_real_, Q_rel2 = NA_real_,
    SeepOut = NA_real_, SeepIn = NA_real_, Bypass = NA_real_,
    Cin = series$Ci,
    C_out = NA_real_, C_treated = NA_real_, C_rel1 = NA_real_, C_rel2 = NA_real_, C_bypass = NA_real_,
    C_seep_recycle = NA_real_, C_seep_discharge = NA_real_,
    Lin = series$Qi * series$Ci,
    L_out = NA_real_, L_treated = NA_real_, L_rel1 = NA_real_, L_rel2 = NA_real_, L_bypass = NA_real_,
    L_seep_recycle = NA_real_, L_seep_discharge = NA_real_,
    stringsAsFactors = FALSE
  )

  water_df <- data.frame(
    Date = as.Date(series$Date),
    RainVol = NA_real_, EtVol = NA_real_, NetAtmo = NA_real_,
    WB_in = NA_real_, WB_out = NA_real_, WB_err = NA_real_, WB_rel = NA_real_,
    stringsAsFactors = FALSE
  )

  mass_df <- data.frame(
    Date = as.Date(series$Date),
    dP = NA_real_,
    Pin_total = NA_real_, Pout_total = NA_real_, Perr_total = NA_real_, Prel_total = NA_real_,
    Pin_external = NA_real_, Pout_external = NA_real_, Perr_external = NA_real_, Prel_external = NA_real_,
    Q_in_tanks = NA_real_, L_in_tanks = NA_real_, C_in_tanks = NA_real_,
    L_rain = NA_real_, L_drydep = NA_real_, L_seepin = NA_real_,
    L_treated = NA_real_, L_rel1 = NA_real_, L_rel2 = NA_real_, L_bypass = NA_real_, L_seep_discharge = NA_real_,
    L_seep_recycle_out = NA_real_,
    L_uptake = NA_real_, L_recycle = NA_real_, L_sed = NA_real_, L_direct = NA_real_,
    stringsAsFactors = FALSE
  )

  if (return_steps) {
    steps_store <- vector("list", n)
  } else {
    steps_store <- NULL
  }

  # MAIN DAY LOOP (delegates everything downstream)
  for (i in seq_len(n)) {
    day_inputs <- as.list(series[i, ])

    nz <- dmsta_zneighbors(i, series$Zcontrol)

    day_inputs$Zcontrol      <- nz$today
    day_inputs$Zcontrol_prev <- nz$prev_day
    day_inputs$Zcontrol_next <- nz$nxt

    day_inputs$has_depth_constraint <- has_depth_constraint_series
    day_inputs$has_Qr0_series <- has_Qr0_series
    day_inputs$has_Qr1_series <- has_Qr1_series
    day_inputs$has_Qr2_series <- has_Qr2_series

    day_inputs$Qr0 <- if (has_Qr0_series) day_inputs$Qr0 else 0
    day_inputs$Qr1 <- if (has_Qr1_series) day_inputs$Qr1 else 0
    day_inputs$Qr2 <- if (has_Qr2_series) day_inputs$Qr2 else 0

    if (is.null(day_inputs$RecycleQ)) day_inputs$RecycleQ <- 0
    if (is.null(day_inputs$RecycleM)) day_inputs$RecycleM <- 0

    if (i > 1 && params$A_cell > 0) Z_hist[i] <- V / params$A_cell
    i0 <- max(1L, i - N_plant + 1L)
    Z_plant <- mean(Z_hist[i0:i], na.rm = TRUE)

    res <- dmsta_flowP_day(
      V = V,
      P_state = P_state,
      tanks = tanks,
      inputs = day_inputs,
      params = params,
      ppar = ppar,
      constants = constants,
      Qmethod = Qmethod,
      Pmethod = Pmethod,
      Nsteps = Nsteps,
      Z_plant = Z_plant,
      integrator_fun = integrator_fun,
      interp_option = interp_option,
      ...
    )

    V <- res$results$V_end
    P_state <- res$results$P_state_end
    Z_hist[i] <- if (params$A_cell > 0) V / params$A_cell else NA_real_

    results_df$V_end[i] <- res$results$V_end
    results_df$Z_end[i] <- res$results$Z_end
    results_df$Z_avg[i] <- res$results$Z_avg
    results_df$V_cell_day[i] <- res$results$V_cell_day

    results_df$Qin[i] <- res$results$Qin
    results_df$Qout[i] <- res$results$Qout
    results_df$Q_treated[i] <- res$results$Q_treated
    results_df$Q_rel1[i] <- res$results$Q_rel1
    results_df$Q_rel2[i] <- res$results$Q_rel2
    results_df$SeepOut[i] <- res$results$SeepOut
    results_df$SeepIn[i] <- res$results$SeepIn
    results_df$Bypass[i] <- res$results$Bypass

    results_df$C_out[i] <- res$results$P$conc$C_out
    results_df$C_treated[i] <- res$results$P$conc$C_treated
    results_df$C_rel1[i] <- res$results$P$conc$C_rel1
    results_df$C_rel2[i] <- res$results$P$conc$C_rel2
    results_df$C_bypass[i] <- res$results$P$conc$C_bypass
    results_df$C_seep_recycle[i] <- res$results$P$conc$C_seep_recycle
    results_df$C_seep_discharge[i] <- res$results$P$conc$C_seep_discharge

    results_df$L_treated[i] <- res$results$P$loads$treated
    results_df$L_rel1[i] <- res$results$P$loads$rel1
    results_df$L_rel2[i] <- res$results$P$loads$rel2
    results_df$L_bypass[i] <- res$results$P$loads$bypass
    results_df$L_seep_recycle[i] <- res$results$P$loads$seep_recycle
    results_df$L_seep_discharge[i] <- res$results$P$loads$seep_discharge
    results_df$L_out[i] <-
      res$results$P$loads$treated +
      res$results$P$loads$rel1 +
      res$results$P$loads$rel2 +
      res$results$P$loads$bypass +
      res$results$P$loads$seep_discharge


    water_df$RainVol[i] <- res$budgets$water$RainVol
    water_df$EtVol[i] <- res$budgets$water$EtVol
    water_df$NetAtmo[i] <- res$budgets$water$NetAtmo
    water_df$WB_in[i] <- res$budgets$water$WB_in
    water_df$WB_out[i] <- res$budgets$water$WB_out
    water_df$WB_err[i] <- res$budgets$water$WB_err
    water_df$WB_rel[i] <- res$budgets$water$WB_rel

    mb <- res$budgets$mass
    if (!is.null(mb)) {
      mass_df$dP[i] <- mb$storage$dP
      mass_df$L_uptake[i] <- mb$mechanisms$L_uptake
      mass_df$L_recycle[i] <- mb$mechanisms$L_recycle
      mass_df$L_sed[i] <- mb$mechanisms$L_sed
      mass_df$L_direct[i] <- mb$mechanisms$L_direct
      mass_df$L_burial[i] <- mb$mechanims$L_burial
      mass_df$L_release[i] <- mb$mechanims$L_release
    }

    if (return_steps) steps_store[[i]] <- res$meta$steps
  }

  meta <- list(
    V_init = V_init,
    V_end = V,
    P_state_end = P_state,
    tanks = tanks,
    ppar = ppar,
    constants = constants,
    Nsteps = Nsteps,
    N_plant = N_plant,
    Qmethod = Qmethod,
    Pmethod = Pmethod,
    has_depth_constraint_series = has_depth_constraint_series,
    has_Qr0_series = has_Qr0_series
  )

  if (return_steps) meta$steps <- steps_store

  out <- list(
    results = results_df,
    budgets = list(water = water_df, mass = mass_df),
    meta = meta
  )

  class(out) <- c("dmsta_result", "list")
  out
}
