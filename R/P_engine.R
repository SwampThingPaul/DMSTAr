
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
dmsta_deriv_mass <- function(state, drivers,ppar,constants) {
  # state: list(M=kg, S=kg)
  # drivers: list with required fields described below

  M <- state$M
  S <- state$S

  A_tank <- drivers$A_tank
  A_cell <- drivers$A_cell

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
    if (ppar$Z_1[k] > 0 && zZ < ppar$Z_1[k]) Fz[k] <- zZ / ppar$Z_1[k]
    # Robust to NA
    # z1 <- ppar$Z_1[k]
    # if (is.finite(z1) && z1 > 0 && is.finite(zZ) && zZ < z1) {
    #   Fz[k] <- zZ / z1
    # } else {
    #   Fz[k] <- 1
    # }

    if (ppar$Chalf[k] > 0) {
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
  P_uptake <- P_recycle <- P_sed <- 0
  P_direct <- 0

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
    if (!is.finite(ppar$Ysigma) || ppar$Ysigma <= 0) {
      Ftrans <- 1
    } else {
      Ftrans <- 1 / (1 + exp(-(Y - ppar$Ytrans) / ppar$Ysigma))
    }

    k1A <- Ftrans * ppar$K1[1] * Fz[1] * fC[1] + (1 - Ftrans) * ppar$K1[2] * Fz[2] * fC[2]
    k2A <- Ftrans * ppar$K2[1]                 + (1 - Ftrans) * ppar$K2[2]
    k3A <- Ftrans * ppar$K3[1]                 + (1 - Ftrans) * ppar$K3[2]

    # reservoir depth penalty blend using Z_plant
    if (ppar$Z_3[1] > 0 && drivers$Z_plant > ppar$Z_2[1]) {
      if (drivers$Z_plant >= ppar$Z_3[1]) {
        fres <- 1
      } else {
        fres <- (drivers$Z_plant - ppar$Z_2[1]) / (ppar$Z_3[1] - ppar$Z_2[1])
      }
      k1A <- fres * ppar$K1[3] * Fz[3] * fC[3] + (1 - fres) * k1A
      k2A <- fres * ppar$K2[3]                 + (1 - fres) * k2A
      k3A <- fres * ppar$K3[3]                 + (1 - fres) * k3A
    }

    P_uptake  <- k1A * C * Y
    P_recycle <- k2A * Y * Y
    k3_used <- k3A
  }

  # SeasonalFactor = [cstar_2] multiplier on uptake
  if (is.finite(ppar$SeasonalFactor) && ppar$SeasonalFactor > 0) {
    xJ <- drivers$julian * 2 * pi / 365.25
    sF <- 1.03926 + ppar$SeasonalFactor * (0.02944 * cos(xJ) - 0.14851 * sin(xJ))
    P_uptake <- P_uptake * sF
  }

  # direct sedimentation (slot 1)
  P_direct <- ppar$K2Coef[1] * C * Fz[1]

  # biomass -> soil flux used as burial shortcut in VBA
  P_sed <- k3_used * Y

  # mass fluxes per area in VBA form
  M_outflow <- (drivers$Qo_tank / A_tank + drivers$Seepout / A_cell) * C
  M_inflow  <- drivers$Li_tank / A_tank +
    (constants$DryDepo + drivers$Rain * constants$C_rain) +
    (drivers$Seepin * constants$seepin_conc / A_cell)

  # derivatives (kg/day)
  dMdt <- (M_inflow + P_recycle - M_outflow - P_uptake - P_direct) * A_tank
  dSdt <- (P_uptake - P_recycle - P_sed) * A_tank

  # limiter (mratiO = 0.01) to prevent collapse below a fraction of Mo/So
  mratiO <- 0.01
  Mo <- drivers$Mo_fix
  So <- drivers$So_fix

  if (Mo > 0 && (Mo + dMdt * drivers$Ddt) <= mratiO * Mo) {
    dMdt2 <- -Mo * (1 - mratiO) / drivers$Ddt
    fF <- M_inflow + P_recycle - M_outflow - dMdt2 / A_tank
    if (fF > 0 && (P_uptake + P_direct) > 0) {
      fF <- fF / (P_uptake + P_direct)
      P_uptake <- P_uptake * fF
      P_direct <- P_direct * fF
    }
    dMdt <- dMdt2
    dSdt <- (P_uptake - P_recycle - P_sed) * A_tank
  }

  if (So > 0 && (So + dSdt * drivers$Ddt) <= mratiO * So) {
    dSdt <- -So * (1 - mratiO) / drivers$Ddt
  }

  list(
    dMdt = dMdt,
    dSdt = dSdt,
    flux = list(P_uptake = P_uptake, P_recycle = P_recycle, P_sed = P_sed, P_direct = P_direct),
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

  stage_call <- function(Ms, Ss, step_frac, Ddt_use) {
    ## internal function for dmsta_rk4_P_step function
    drv <- args_base$drivers
    drv$StepFrac <- step_frac
    drv$Ddt      <- Ddt_use
    drv$Mo_fix   <- Mo
    drv$So_fix   <- So
    dmsta_deriv_mass(state = list(M = Ms, S = Ss),
                     drivers = drv,
                     ppar = ppar,
                     constants = constants
    )
  }

  Mo <- M
  So <- S

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
    if (!is.finite(M_new)) M_new <- Mmin
    if (!is.finite(S_new)) S_new <- Smin
    if (M_new < Mmin) M_new <- Mmin
    if (S_new < Smin) S_new <- Smin
  }

  #  averaged flux diagnostics over step
  P_uptake_ts  <- (s1$flux$P_uptake  + 2*s2$flux$P_uptake  + 2*s3$flux$P_uptake  + s4$flux$P_uptake)  / 6
  P_recycle_ts <- (s1$flux$P_recycle + 2*s2$flux$P_recycle + 2*s3$flux$P_recycle + s4$flux$P_recycle) / 6
  P_sed_ts     <- (s1$flux$P_sed     + 2*s2$flux$P_sed     + 2*s3$flux$P_sed     + s4$flux$P_sed)     / 6
  P_direct_ts  <- (s1$flux$P_direct  + 2*s2$flux$P_direct  + 2*s3$flux$P_direct  + s4$flux$P_direct)  / 6

  #  averaged diag diagnostics over step (C, Y, z)
  C_ts <- (s1$diag$C + 2*s2$diag$C + 2*s3$diag$C + s4$diag$C) / 6
  Y_ts <- (s1$diag$Y + 2*s2$diag$Y + 2*s3$diag$Y + s4$diag$Y) / 6
  z_ts <- (s1$diag$z + 2*s2$diag$z + 2*s3$diag$z + s4$diag$z) / 6

  list(
    M_new = M_new,
    S_new = S_new,
    dMdt_ts = dMdt_ts,
    dSdt_ts = dSdt_ts,
    P_uptake_ts  = P_uptake_ts,
    P_recycle_ts = P_recycle_ts,
    P_sed_ts     = P_sed_ts,
    P_direct_ts  = P_direct_ts,
    C_ts = C_ts,
    Y_ts = Y_ts,
    z_ts = z_ts
  )
}

#' Run one day of coupled hydrology + phosphorus (internal)
#'
#' Integrates DMSTA hydrology over one day (via `dmsta_flow_day_steps()`)
#' and then integrates phosphorus mass balance across tanks using RK4 substeps.
#' Returns daily totals for flows and loads, updated tank states, and optional
#' water and mass budget diagnostics.
#'
#' @param V Numeric scalar. Starting cell volume (e.g., hm^3).
#' @param P_state Named list with `M` and `S` vectors (kg), each of
#'   length `tanks$Ntanks`.
#' @param tanks Tank geometry list as returned by `dmsta_build_tanks()`,
#'   including `Ntanks`, `A_Tank`, `F_Tank`, and `Fcum`.
#' @param inputs Named list of day forcings and controls. Must include at least
#'   `Date`, `Qi`, `Ci`, `Rain`, `Et`, and `Zcontrol`
#'   (plus optional releases and recycle fields used downstream).
#' @param params Hydrology parameter list passed to `dmsta_flow_day_steps()`.
#'   Must include `A_cell`.
#' @param ppar Phosphorus kinetic parameter list, typically from `build_P_kin_slots()`.
#' @param constants Constants list used by `dmsta_deriv_mass()`.
#' @param Nsteps Integer. Number of RK substeps per day.
#' @param Z_plant Numeric scalar. Rolling/representative plant depth (m) used for blending.
#'
#' @return A named list with components:
#' \describe{
#'   \item{results}{Daily totals and end-of-day state, including `P_state_end`.}
#'   \item{budgets}{Lists of `water` and `mass` budgets (mass may be `NULL`).}
#'   \item{meta}{Metadata and (by default) hydrology substep details for debugging/coupling.}
#' }
#'
#' @rdname internal_dmsta_phosphorus
#' @keywords internal
dmsta_flowP_day <- function(V, P_state, tanks, inputs, params,
                            ppar, constants,
                            Nsteps = 4L, Z_plant = 0) {

  A_cell <- params$A_cell
  isa_node <- dmsta_is_node(A_cell, params$IsaNode)

  if (isa_node) {
    # hydrology already returns node behavior (after patch above)
    hyd <- dmsta_flow_day_steps(V = 0, inputs = inputs, params = params, Nsteps = Nsteps)

    Qi_eff <- params$Qin_Frac * inputs$Qi
    Ci <- inputs$Ci
    RecycleQ <- if (is.null(inputs$RecycleQ)) 0 else inputs$RecycleQ
    RecycleM <- if (is.null(inputs$RecycleM)) 0 else inputs$RecycleM

    # node “outlet” concentration: inflow+recycle weighted (VBA effectively uses inflow conc)
    Qin_total <- Qi_eff - hyd$Bypass + RecycleQ
    Lin_total <- (Qi_eff - hyd$Bypass) * Ci + RecycleM
    Cout <- if (Qin_total > 0) Lin_total / Qin_total else 0

    # assign all non-bypass outflow to “treated” stream in node mode
    Q_treat <- hyd$Qout
    L_treat <- Q_treat * Cout
    Q_byp <- hyd$Bypass
    L_byp <- Q_byp * Ci

    seepC <- Cout
    if (!is.null(constants$seepout_conc_max) && constants$seepout_conc_max > 0) {
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
          treated = Q_treat, rel1 = 0, rel2 = 0, bypass = Q_byp,
          seep_recycle = Q_seep_rec, seep_discharge = Q_seep_dis
        ),
        loads = list(
          treated = L_treat, rel1 = 0, rel2 = 0, bypass = L_byp,
          seep_recycle = L_seep_rec, seep_discharge = L_seep_dis
        ),
        conc = list(
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

    water_budget <- list(
      RainVol = hyd$RainVol,
      EtVol   = hyd$EtVol,
      NetAtmo = hyd$NetAtmo,
      WB_in   = hyd$WB_in,
      WB_out  = hyd$WB_out,
      WB_err  = hyd$WB_err,
      WB_rel  = hyd$WB_rel
    )

    # node has no tanks, so storage is defined as zero and unchanged
    mass_budget <- list(
      storage = list(
        M_start = 0, S_start = 0, P_start = 0,
        M_end = 0,   S_end = 0,   P_end = 0,
        dM = 0, dS = 0, dP = 0
      ),
      inflow_tanks = list(
        Q_in_tanks = Qin_total,
        L_in_tanks = Lin_total,
        C_in_tanks = if (Qin_total > 0) Lin_total / Qin_total else NA_real_,
        Q_in_flow = max(0, Qi_eff - hyd$Bypass),
        L_in_flow = (max(0, Qi_eff - hyd$Bypass)) * Ci,
        C_in_flow = if (max(0, Qi_eff - hyd$Bypass) > 0) Ci else NA_real_,
        Q_in_recycle = RecycleQ,
        L_in_recycle = RecycleM
      ),
      inputs_external = list(
        # you can include these if you want; otherwise set to 0
        L_rain = (inputs$Rain * constants$C_rain) * 0,
        L_drydep = (constants$DryDepo) * 0,
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
        L_uptake = 0, L_recycle = 0, L_sed = 0, L_direct = 0
      ),
      closure = {
        Pin_total  <- Lin_total # watershed+recycle already rolled into Lin_total
        Pout_total <- (L_treat + L_byp + L_seep_dis + L_seep_rec)
        dP <- 0

        Perr_total <- dP - (Pin_total - Pout_total)
        Prel_total <- Perr_total / max(1e-12, max(Pin_total, Pout_total))

        list(
          Pin_total = Pin_total,
          Pout_total = Pout_total,
          Perr_total = Perr_total,
          Prel_total = Prel_total,
          Pin_external = Pin_total,
          Pout_external = (L_treat + L_byp + L_seep_dis),
          Perr_external = dP - (Pin_total - (L_treat + L_byp + L_seep_dis)),
          Prel_external = (dP - (Pin_total - (L_treat + L_byp + L_seep_dis))) /
            max(1e-12, max(Pin_total, (L_treat + L_byp + L_seep_dis)))
        )
      }
    )

    out <- list(
      results = results,
      budgets = list(
        water = water_budget,
        mass = mass_budget
      ),
      meta = list(Date = inputs$Date, IsaNode = TRUE)
    )
  }else{

    #  1) Hydrology with step series
    hyd <- dmsta_flow_day_steps(V, inputs, params, Nsteps = Nsteps)
    step_out <- hyd$steps

    # 1b) Daily-average water level diagnostics (VBA-style Vavg*Dt sum)
    V_cell_day <- sum(vapply(step_out, function(s) 0.5 * (s$Vo + s$V) * s$Dt, 0.0)) # hm3
    A_cell <- params$A_cell
    Z_avg <- if (A_cell > 0) V_cell_day / A_cell else NA_real_
    Z_end <- if (A_cell > 0) hyd$V_end / A_cell else NA_real_

    #  2) Ensure P state vectors match Ntanks
    Nt <- tanks$Ntanks
    if (length(P_state$M) != Nt || length(P_state$S) != Nt) {
      stop("P_state$M and P_state$S must have length equal to tanks$Ntanks.")
    }
    M <- P_state$M
    S <- P_state$S

    # P budget storage starts (optional)
    M_start <- sum(M);
    S_start <- sum(S);
    P_start <- M_start + S_start

    #  3) Daily accumulators (flows hm3, loads kg)
    Q_treat <- L_treat <- 0
    Q_r1    <- L_r1    <- 0
    Q_r2    <- L_r2    <- 0
    Q_byp   <- L_byp   <- 0

    # optional seepage bookkeeping (matches DMSTA seepage conc cap pattern)
    Q_seep_rec <- L_seep_rec <- 0
    Q_seep_dis <- L_seep_dis <- 0

    # P budget accumulators for budget
    # treated inflow into tanks (includes recycle)
    Q_in_tanks <- 0
    L_in_tanks <- 0

    # split inflow into external flow part vs internal recycle transfer
    Q_in_flow    <- 0
    L_in_flow    <- 0
    Q_in_recycle <- 0
    L_in_recycle <- 0

    # internal mechanisms
    L_uptake <- 0
    L_recycle <- 0
    L_sed <- 0
    L_direct <- 0

    # convenience
    jul <- julian_day(inputs$Date)
    A_cell <- params$A_cell

    #  4) Substep loop
    for (k in seq_along(step_out)) {
      hs <- step_out[[k]]
      Dt <- hs$Dt

      Vo <- hs$Vo
      V1 <- hs$V
      Vavg <- 0.5 * (Vo + V1)

      # hydrology rates for this substep
      qout_total <- hs$Qout
      byp_rate   <- hs$Bypass
      seepout    <- hs$SeepOut
      seepin     <- hs$SeepIn

      # treated/release split rates for this step
      q_treated <- hs$Q_treated
      q_rel1    <- hs$Q_rel1
      q_rel2    <- hs$Q_rel2

      # inflow concentration & (cell-scaled) inflow flow rate
      Ci <- inputs$Ci
      Qi_eff <- hs$Qi_eff

      # inflow to tank system after bypass + recycle (matches VBA for Tank=1)
      RecycleQ <- if (is.null(inputs$RecycleQ)) 0 else inputs$RecycleQ
      RecycleM <- if (is.null(inputs$RecycleM)) 0 else inputs$RecycleM

      Qin_total <- Qi_eff - byp_rate + RecycleQ
      if (!is.finite(Qin_total) || Qin_total < 0) Qin_total <- 0

      # bypass load uses inflow conc (VBA Mt(3)=Bypass*Ci)
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

      # routing between tanks
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

        # DMSTA tank outflow interpolation (Qo_tank)
        Qo_tank <- Qin_total * (1 - Fcum) + qout_total * Fcum
        if (!is.finite(Qo_tank) || Qo_tank < 0) Qo_tank <- 0

        # Build drivers for derivative/RK4
        # NOTE: dmsta_deriv_mass currently references 'ppar' and 'constants' from scope.
        #       Ensure those are visible (globals) OR wrap dmsta_deriv_mass in a closure/module.
        args_base <- list(drivers = list(
          A_tank = A_tk,
          A_cell = A_cell,
          V_tank_avg = V_tank_avg,
          Vdel = Vdel,
          Qo_tank = Qo_tank,
          Li_tank = Li_tank,
          Rain = inputs$Rain,
          Seepout = seepout,
          Seepin  = seepin,
          Z_plant = Z_plant,
          julian  = jul
        ))

        # Integrate M,S over this substep (RK4)
        Mo <- M[tk]
        So <- S[tk]
        rk <- dmsta_rk4_P_step(M = Mo, S = So, args_base = args_base, Dt = Dt,
                               ppar = ppar,constants = constants)

        M[tk] <- rk$M_new
        S[tk] <- rk$S_new

        # P budget: mechanism flux totals (optional)

        # rk$P_*_ts are mean rates per area in this Dt interval
        L_uptake  <- L_uptake  + rk$P_uptake_ts  * A_tk * Dt
        L_recycle <- L_recycle + rk$P_recycle_ts * A_tk * Dt
        L_sed     <- L_sed     + rk$P_sed_ts     * A_tk * Dt
        L_direct  <- L_direct  + rk$P_direct_ts  * A_tk * Dt


        # VBA routes using Cavg = Mavg / V_tank (not instantaneous C)
        Mavg <- 0.5 * (Mo + M[tk])
        Cavg <- Mavg / max(V_tank_avg, 1e-12)

        Lo_tank <- Qo_tank * Cavg
        Qo_prev <- Qo_tank
        Lo_prev <- Lo_tank

        # last tank = cell outlet concentration for this substep
        if (tk == Nt) {
          # treated/release loads (kg) for this step
        Q_treat <- Q_treat + q_treated * Dt
        L_treat <- L_treat + (q_treated * Cavg) * Dt

        Q_r1 <- Q_r1 + q_rel1 * Dt
        L_r1 <- L_r1 + (q_rel1 * Cavg) * Dt

        Q_r2 <- Q_r2 + q_rel2 * Dt
        L_r2 <- L_r2 + (q_rel2 * Cavg) * Dt

        # optional seepage recycle/discharge bookkeeping w/ conc cap
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
        L_seep_rec <- L_seep_rec + (Qsr * seepC) * Dt

        Q_seep_dis <- Q_seep_dis + Qsd * Dt
        L_seep_dis <- L_seep_dis + (Qsd * seepC) * Dt
      }
    } # tank loop
  } # substep loop

  # P budget: day-level atmospheric + seep inputs and storage end (optional)
  mass_budget <- NULL
  M_end <- sum(M); S_end <- sum(S); P_end <- M_end + S_end
  dM <- M_end - M_start
  dS <- S_end - S_start
  dP <- P_end - P_start

  # external (non-transfer) inputs to tanks are L_in_flow (treated inflow excluding recycle)
  # plus atmospheric + seep-in inputs (apply to the whole cell area)
  L_rain   <- (inputs$Rain * constants$C_rain) * A_cell
  L_drydep <- (constants$DryDepo) * A_cell
  L_seepin <- hyd$SeepIn * constants$seepin_conc

  # define outflows:
  L_out_external <- L_treat + L_r1 + L_r2 + L_byp + L_seep_dis
  L_out_transfer <- L_seep_rec  # internal transfer stream
  L_in_external  <- L_in_flow + L_rain + L_drydep + L_seepin
  L_in_transfer  <- L_in_recycle

  # closures:
  # total closure includes transfers
  Pin_total  <- L_in_external + L_in_transfer
  Pout_total <- L_out_external + L_out_transfer
  Perr_total <- dP - (Pin_total - Pout_total)
  Prel_total <- Perr_total / max(1e-12, max(Pin_total, Pout_total))

  # external-only closure (will not close if transfers exist and you don't track transit)
  Perr_external <- dP - (L_in_external - L_out_external)
  Prel_external <- Perr_external / max(1e-12, max(L_in_external, L_out_external))

  mass_budget <- list(
    storage = list(
      M_start = M_start, S_start = S_start, P_start = P_start,
      M_end   = M_end,   S_end   = S_end,   P_end   = P_end,
      dM = dM, dS = dS, dP = dP
    ),
    inflow_tanks = list(
      Q_in_tanks = Q_in_tanks,
      L_in_tanks = L_in_tanks,
      C_in_tanks = fw(L_in_tanks, Q_in_tanks),
      Q_in_flow  = Q_in_flow,
      L_in_flow  = L_in_flow,
      C_in_flow  = fw(L_in_flow, Q_in_flow),
      Q_in_recycle = Q_in_recycle,
      L_in_recycle = L_in_recycle
    ),
    inputs_external = list(
      L_rain = L_rain, L_drydep = L_drydep, L_seepin = L_seepin
    ),
    outputs_external = list(
      L_treated = L_treat, L_rel1 = L_r1, L_rel2 = L_r2,
      L_bypass = L_byp, L_seep_discharge = L_seep_dis
    ),
    transfers = list(
      L_seep_recycle_out = L_seep_rec,
      Q_seep_recycle_out = Q_seep_rec
    ),
    mechanisms = list(
      L_uptake = L_uptake, L_recycle = L_recycle,
      L_sed = L_sed, L_direct = L_direct
    ),
    closure = list(
      Pin_total = Pin_total, Pout_total = Pout_total,
      Perr_total = Perr_total, Prel_total = Prel_total,
      Pin_external = L_in_external, Pout_external = L_out_external,
      Perr_external = Perr_external, Prel_external = Prel_external
    )
  )

  #  5) Return combined hydro + P day outputs
  # WATER budget
  water_budget <- list(
    RainVol = hyd$RainVol,
    EtVol   = hyd$EtVol,
    NetAtmo = hyd$NetAtmo,
    WB_in   = hyd$WB_in,
    WB_out  = hyd$WB_out,
    WB_err  = hyd$WB_err,
    WB_rel  = hyd$WB_rel
  )

  # RESULTS
  results <- list(
    # state / levels
    V_end      = hyd$V_end,
    Z_end      = Z_end,
    Z_avg      = Z_avg,
    V_cell_day = V_cell_day,
    Qin        = hyd$Qin,

    # hydrology totals
    Qout      = hyd$Qout,
    Q_treated = hyd$Q_treated,
    Q_rel1    = hyd$Q_rel1,
    Q_rel2    = hyd$Q_rel2,
    SeepOut   = hyd$SeepOut,
    SeepIn    = hyd$SeepIn,
    Bypass    = hyd$Bypass,

    # P outputs (daily totals)
    P = list(
      flows = list(
        treated = Q_treat, rel1 = Q_r1, rel2 = Q_r2, bypass = Q_byp,
        seep_recycle = Q_seep_rec, seep_discharge = Q_seep_dis
      ),
      loads = list(
        treated = L_treat, rel1 = L_r1, rel2 = L_r2, bypass = L_byp,
        seep_recycle = L_seep_rec, seep_discharge = L_seep_dis
      ),
      conc = list(
        C_treated = fw(L_treat, Q_treat),
        C_rel1    = fw(L_r1, Q_r1),
        C_rel2    = fw(L_r2, Q_r2),
        C_out     = fw(L_treat + L_r1 + L_r2, Q_treat + Q_r1 + Q_r2),
        C_bypass  = fw(L_byp, Q_byp),
        C_seep_recycle   = fw(L_seep_rec, Q_seep_rec),
        C_seep_discharge = fw(L_seep_dis, Q_seep_dis)
      )
    ),

    # updated P state
    P_state_end = list(M = M, S = S)
  )

  # META
  meta <- list(
    Date   = inputs$Date,
    V_start = V,
    V_end   = hyd$V_end,
    A_cell  = A_cell,
    Nsteps  = Nsteps,
    Z_plant = Z_plant,
    steps   = hyd$steps,       # keep for debugging/coupling
    inputs_used = inputs,
    params_used = params
  )

  out <- list(
    results = results,
    budgets = list(
      water = water_budget,
      mass  = mass_budget
    ),
    meta = meta
  )
  }

  out
}

#' Run coupled DMSTA hydrology and phosphorus over a time series
#'
#' Runs a coupled DMSTA simulation over multiple days, integrating hydrology and
#' phosphorus mass balance through a series of conceptual tanks. Hydrology is
#' computed with sub-daily RK steps and phosphorus is integrated over the same
#' substeps for each tank.
#'
#' The function can build tank geometry, phosphorus kinetics, constants, and
#' initial states automatically from `params`, or accept precomputed objects.
#'
#' @param series A `data.frame` containing (at minimum) columns:
#'   \describe{
#'     \item{Date}{Date column coercible via `as.Date()`.}
#'     \item{Qi}{Daily inflow rate (volume/day).}
#'     \item{Ci}{Inflow concentration (ppb).}
#'     \item{Rain}{Rain rate (m/day).}
#'     \item{Et}{ET rate (m/day).}
#'     \item{Zcontrol}{Control depth (m).}
#'   }
#'   Optional columns include `Qr0`, `Qr1`, `Qr2` (release components),
#'   `RecycleQ` (recycle flow), and `RecycleM` (recycle mass load).
#'
#' @param params Named list of hydrology + phosphorus parameters. Must include at least:
#'   \describe{
#'     \item{A_cell}{Cell area.}
#'     \item{DutyCycle}{Duty-cycle multiplier used for kinetics building (if `ppar` not supplied).}
#'     \item{Zinit}{Initial depth (cm) used to derive `V_init` if not provided.}
#'     \item{C_init_ppb}{Initial water concentration (ppb) for state initialization.}
#'     \item{Y_init_mgm2}{Initial areal store (mg/m^2) for state initialization.}
#'   }
#'   Additional fields are used if `ppar` and/or `constants` are not supplied
#'   (e.g., deposition, seepage concentrations, and kinetic builder inputs).
#'
#' @param pparams Optional named list of phosphorus parameters to merge into `params`
#'   (via `modifyList()`); useful if P parameters are maintained separately.
#' @param ttankS Numeric scalar > 0. Effective number of tanks; may be fractional.
#' @param Nsteps Integer. RK substeps per day (default `4L`).
#' @param N_plant Integer. Window length (days) for rolling mean depth `Z_plant`.
#' @param ppar Optional precomputed kinetics parameter object from `build_P_kin_slots()`.
#'   If `NULL`, it is built internally for modules `STA`, `PSTA`, and `RES`.
#' @param constants Optional constants list used by `dmsta_deriv_mass()`.
#'   If `NULL`, defaults are constructed from `params`.
#' @param tanks Optional tank geometry from `dmsta_build_tanks()`. If `NULL`, built from
#'   `params$A_cell` and `ttankS`.
#' @param V_init Optional initial volume (same units as hydrology volume). If `NULL`,
#'   derived from `params$Zinit` and `params$A_cell`.
#' @param init_P_state Optional initial phosphorus state list `list(M, S)` with vectors
#'   of length `tanks$Ntanks`. If `NULL`, constructed from `params` via
#'   `dmsta_p_init_state()`.
#' @param return_steps Logical; if `TRUE`, store per-day hydrology substep details in
#'   `meta$steps`. This can be large.
#'
#' @importFrom utils modifyList
#' @details
#' Internally, control-depth neighbor values are obtained using `neighbors_zcontrol()`.
#' The function also uses `build_P_kin_slots()` and `validate_P_paramsK()` when
#' `ppar` is not supplied.
#'
#' @return An object of class `"dmsta_result"` (a list) with elements:
#' \describe{
#'   \item{results}{A data.frame of daily volumes, flows, concentrations, and loads.}
#'   \item{budgets}{A list with `water` and `mass` budget data.frames.}
#'   \item{meta}{Run metadata, including final states and optional per-day step details.}
#' }
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
#'  DutyCycle = 0.95, # duty [for P cycling]
#'  Cmax = 2000
#')
#'
#' pparams <- list(
#'  C1000 = 22,        # _c1000_1
#'  Cstar = 3,         # cstar_1
#'  Ks_per_yr = 16.8, # ks_1
#'  Z1 = 40,           # z_1
#'  Z2 = 100,          # z_2
#'  Z3 = 200,          # z_3
#'  Chalf = 300,       # _kc2
#'  K2Coef1 = 0,              # k2res
#'  Ytrans = 0,             # y00
#'  Ysigma = 0,             # ys00
#'  Czero = 0,                # CZero = 0
#'  C_rain = 10,              # rainfall P concentration (ug/L)
#'  DryDepo = 20,             # atmos (mg/m2-yr)
#'  SeasonalFactor = 0,       # cstar_2 ##OFF for base parity
#'  C1000_2 = NULL,           # _c1000_2
#'  ks_2 = 0,                 # ks_2
#'  zh_2 = 0,                 # zh_2
#'  k_depth_penalty = 1,      # m/yr
#'  seepage_c = 20,           # max seepage outflow concentration (ug/L)
#'  seepin_conc = 0,          # Inflow Seepage Conc (ug/L)
#'  C_init_ppb = 30,                # water column initial conc (ug/L)
#'  Y_init_mgm2 = 3387.67297548954, # biomass storage (mg/m2)
#'  n_tanks = 3,
#'  Nsteps = 4
#' )
#'
#' params <- modifyList(params,pparams)
#'
#' ##  choose model structure
#' ttankS  <- params$n_tanks    # tanks in series (can be fractional)
#' Nsteps  <- params$Nsteps     # RK substeps per day
#'
#' ##  build tank geometry once
#' tanks <- dmsta_build_tanks(params$A_cell, ttankS)
#'
#' ##  build kinetics once: 3 modules STA/PSTA/RES
#' P_MODEL_BUILDERS <- list(
#' STA  = build_STA,
#' PSTA = build_PSTA,
#' RES  = build_RES
#' )
#' ppar <- build_P_kin_slots(
#'   mods     = c("STA", "PSTA", "RES"),
#'   registry = NULL,
#'   pparams  = params,
#'   Dpy      = 365.25,
#'   DutyCycle = params$DutyCycle
#' )
#'  validate_P_paramsK(ppar)
#'
#' ##  constants used by dmsta_deriv_mass
#' constants <- list(
#'   Cmax = params$Cmax,
#'   C_rain = params$C_rain,
#'   DryDepo = params$DryDepo / 365.25,   # mg/m2-day (VBA divides by Dpy)
#'   seepin_conc = params$seepin_conc,
#'   seepout_conc_max = params$seepage_c, # cap seepage concentration
#'   fseep_recycle = 0,                   # set later if you add seepage recycle fractions
#'   fseep_out = 0
#' )
#'
#' ##  initial conditions
#' Z_init_m <- cm_to_m(params$Zinit)
#' V_init   <- params$A_cell * Z_init_m
#'
#' P_state0 <- dmsta_p_init_state(
#'   tanks,
#'   Z_init_m      = Z_init_m,
#'   C_init_ppb    = params$C_init_ppb,
#'   Y_init_mgm2   = params$Y_init_mgm2
#' )
#'
#' out <- dmsta_flowP_series(
#'   series = series,
#'   params = params,
#'   ttankS = ttankS,
#'   Nsteps = Nsteps,
#'   tanks  = tanks,
#'   ppar   = ppar,
#'   constants = constants,
#'   V_init = V_init,
#'   init_P_state = P_state0,
#'   return_steps  = FALSE
#'  )
#'  out$results
#'  }
#'
#' @export

dmsta_flowP_series <- function(
    series,               # data.frame with Date, Qi, Ci, Rain, Et, Zcontrol (+ optional Qr0/Qr1/Qr2/RecycleQ/RecycleM)
    params,               # hydrology + P params merged
    pparams = NULL,       # optional separate P params (if not already merged into params)
    ttankS = 3.0,         # tanks in series (can be fractional)
    Nsteps = 4L,          # RK substeps per day
    N_plant = 30L,        # rolling depth window (days) for Z_plant
    ppar = NULL,          # optional precomputed output of build_P_kin_slots(...)
    constants = NULL,     # optional constants list used by dmsta_deriv_mass
    tanks = NULL,         # optional pre-built tanks from dmsta_build_tanks
    V_init = NULL,        # optional initial volume hm3 (default from Zinit*A_cell)
    init_P_state = NULL,  # optional initial P state list(M,S); else uses params C_init_ppb/Y_init_mgm2
    return_steps = FALSE  # if TRUE, store hydrology substep list per day (big)
) {

  #  basic checks
  if (!is.data.frame(series)) stop("'series' must be a data.frame.")
  req <- c("Date", "Qi", "Ci", "Rain", "Et", "Zcontrol")
  miss <- setdiff(req, names(series))
  if (length(miss) > 0) stop("series is missing columns: ", paste(miss, collapse = ", "))

  n <- nrow(series)
  if (n < 1) stop("series has zero rows.")

  # Merge params + pparams if provided
  if (!is.null(pparams)) params <- modifyList(params, pparams)

  #  build tanks once
  if (is.null(tanks)) {
    tanks <- dmsta_build_tanks(params$A_cell, ttankS)
  }

  #  build kinetics ppar once (3 modules STA/PSTA/RES)
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

  #  constants (match what dmsta_deriv_mass expects)
  if (is.null(constants)) {
    constants <- list(
      Cmax = if (is.null(params$Cmax)) 2000 else params$Cmax,
      C_rain = if (is.null(params$C_rain)) 0 else params$C_rain,
      DryDepo = if (is.null(params$DryDepo)) 0 else (params$DryDepo / 365.25),  # mg/m2-day
      seepin_conc = if (!is.null(params$Seepin_Conc)) params$Seepin_Conc else 0,
      seepout_conc_max = if (!is.null(params$seepage_c)) params$seepage_c else 0,
      fseep_recycle = if (!is.null(params$fseep_recycle)) params$fseep_recycle else 0,
      fseep_out = if (!is.null(params$fseep_out)) params$fseep_out else 0
    )
  }

  #  initial volume
  if (is.null(V_init)) {
    Z_init_m <- (params$Zinit / 100)
    V_init <- params$A_cell * Z_init_m
  }
  V <- V_init

  #  initial P state
  if (is.null(init_P_state)) {
    Z_init_m <- (params$Zinit / 100)
    init_P_state <- dmsta_p_init_state(
      tanks,
      Z_init_m = Z_init_m,
      C_init_ppb = params$C_init_ppb,
      Y_init_mgm2 = params$Y_init_mgm2
    )
  }
  P_state <- init_P_state

  #  rolling depth history for Z_plant
  Z_hist <- numeric(n)
  # initialize with starting depth
  Z_hist[1] <- V / params$A_cell

  # RESULTS table (primary outputs)
  results_df <- data.frame(
    Date = as.Date(series$Date),

    # state/levels
    V_end = NA_real_, Z_end = NA_real_, Z_avg = NA_real_, V_cell_day = NA_real_,

    # hydrology
    Qin = NA_real_,
    Qout = NA_real_, Q_treated = NA_real_, Q_rel1 = NA_real_, Q_rel2 = NA_real_,
    SeepOut = NA_real_, SeepIn = NA_real_, Bypass = NA_real_,

    # P concentrations (ppb)
    Cin = series$Ci,
    C_out = NA_real_, C_treated = NA_real_, C_rel1 = NA_real_, C_rel2 = NA_real_, C_bypass = NA_real_,
    C_seep_recycle = NA_real_, C_seep_discharge = NA_real_,

    # P loads (kg)
    Lin = series$Qi * series$Ci,
    L_out = NA_real_, L_treated = NA_real_, L_rel1 = NA_real_, L_rel2 = NA_real_, L_bypass = NA_real_,
    L_seep_recycle = NA_real_, L_seep_discharge = NA_real_,

    stringsAsFactors = FALSE
  )

  # WATER budget table (separate)
  water_df <- data.frame(
    Date = as.Date(series$Date),
    RainVol = NA_real_, EtVol = NA_real_, NetAtmo = NA_real_,
    WB_in = NA_real_, WB_out = NA_real_, WB_err = NA_real_, WB_rel = NA_real_,
    stringsAsFactors = FALSE
  )

  # MASS budget table (separate; optional)
  # We'll store a summary row per day; optionally a full list column
  mass_df <- NULL
  mass_list <- NULL

  mass_df <- data.frame(
    Date = as.Date(series$Date),

    # storage change
    dP = NA_real_,

    # closures (total and external)
    Pin_total = NA_real_, Pout_total = NA_real_, Perr_total = NA_real_, Prel_total = NA_real_,
    Pin_external = NA_real_, Pout_external = NA_real_, Perr_external = NA_real_, Prel_external = NA_real_,

    # tank inflow totals
    Q_in_tanks = NA_real_, L_in_tanks = NA_real_, C_in_tanks = NA_real_,

    # external inputs (kg)
    L_rain = NA_real_, L_drydep = NA_real_, L_seepin = NA_real_,

    # external outputs (kg)
    L_treated = NA_real_, L_rel1 = NA_real_, L_rel2 = NA_real_, L_bypass = NA_real_, L_seep_discharge = NA_real_,

    # internal transfer (seep recycle)
    L_seep_recycle_out = NA_real_,

    # mechanisms
    L_uptake = NA_real_, L_recycle = NA_real_, L_sed = NA_real_, L_direct = NA_real_,

    stringsAsFactors = FALSE
  )

  if (return_steps) {
    steps_store <- vector("list", n)
  } else {
    steps_store <- NULL
  }

  #  main day loop
  for (i in seq_len(n)) {
    day_inputs <- as.list(series[i, ])

    # Zcontrol neighbors
    nz <- neighbors_zcontrol(i, series$Zcontrol)
    day_inputs$Zcontrol      <- nz$today
    day_inputs$Zcontrol_prev <- nz$prev
    day_inputs$Zcontrol_next <- nz$nxt
    day_inputs$has_depth_constraint <- is.finite(day_inputs$Zcontrol) && day_inputs$Zcontrol != 0

    # ensure optional fields exist
    if (is.null(day_inputs$RecycleQ)) day_inputs$RecycleQ <- 0
    if (is.null(day_inputs$RecycleM)) day_inputs$RecycleM <- 0
    if (is.null(day_inputs$Qr0)) day_inputs$Qr0 <- 0
    if (is.null(day_inputs$Qr1)) day_inputs$Qr1 <- 0
    if (is.null(day_inputs$Qr2)) day_inputs$Qr2 <- 0

    # compute Z_plant (rolling mean depth)
    if (i > 1) Z_hist[i] <- V / params$A_cell
    i0 <- max(1L, i - N_plant + 1L)
    Z_plant <- mean(Z_hist[i0:i])

    # run coupled day
    res <- dmsta_flowP_day(
      V = V,
      P_state = P_state,
      tanks = tanks,
      inputs = day_inputs,
      params = params,
      ppar = ppar,
      constants = constants,
      Nsteps = Nsteps,
      Z_plant = Z_plant
    )

    # advance state
    V <- res$results$V_end
    P_state <- res$results$P_state_end

    # store depth history for next day
    Z_hist[i] <- V / params$A_cell

    # fill RESULTS
    results_df$V_end[i] <- res$results$V_end
    results_df$Z_end[i] <- res$results$Z_end
    results_df$Z_avg[i] <- res$results$Z_avg
    results_df$V_cell_day[i] <- res$results$V_cell_day

    results_df$Qin[i]       <- res$results$Qin
    results_df$Qout[i]      <- res$results$Qout
    results_df$Q_treated[i] <- res$results$Q_treated
    results_df$Q_rel1[i]    <- res$results$Q_rel1
    results_df$Q_rel2[i]    <- res$results$Q_rel2
    results_df$SeepOut[i]   <- res$results$SeepOut
    results_df$SeepIn[i]    <- res$results$SeepIn
    results_df$Bypass[i]    <- res$results$Bypass

    # P outputs
    results_df$C_out[i]          <- res$results$P$conc$C_out
    results_df$C_treated[i]      <- res$results$P$conc$C_treated
    results_df$C_rel1[i]         <- res$results$P$conc$C_rel1
    results_df$C_rel2[i]         <- res$results$P$conc$C_rel2
    results_df$C_bypass[i]       <- res$results$P$conc$C_bypass
    results_df$C_seep_recycle[i] <- res$results$P$conc$C_seep_recycle
    results_df$C_seep_discharge[i] <- res$results$P$conc$C_seep_discharge

    results_df$L_treated[i]      <- res$results$P$loads$treated
    results_df$L_rel1[i]         <- res$results$P$loads$rel1
    results_df$L_rel2[i]         <- res$results$P$loads$rel2
    results_df$L_bypass[i]       <- res$results$P$loads$bypass
    results_df$L_seep_recycle[i] <- res$results$P$loads$seep_recycle
    results_df$L_seep_discharge[i] <- res$results$P$loads$seep_discharge
    results_df$L_out[i] <- res$results$P$loads$treated + res$results$P$loads$rel1 + res$results$P$loads$rel2

    # fill WATER budgets
    water_df$RainVol[i] <- res$budgets$water$RainVol
    water_df$EtVol[i]   <- res$budgets$water$EtVol
    water_df$NetAtmo[i] <- res$budgets$water$NetAtmo
    water_df$WB_in[i]   <- res$budgets$water$WB_in
    water_df$WB_out[i]  <- res$budgets$water$WB_out
    water_df$WB_err[i]  <- res$budgets$water$WB_err
    water_df$WB_rel[i]  <- res$budgets$water$WB_rel

    # fill MASS budgets
    mb <- res$budgets$mass

    # mb can be NULL if caller asked FALSE or something went wrong
    if (!is.null(mb)) {
      mass_df$dP[i] <- mb$storage$dP

      mass_df$Pin_total[i] <- mb$closure$Pin_total
      mass_df$Pout_total[i] <- mb$closure$Pout_total
      mass_df$Perr_total[i] <- mb$closure$Perr_total
      mass_df$Prel_total[i] <- mb$closure$Prel_total

      mass_df$Pin_external[i] <- mb$closure$Pin_external
      mass_df$Pout_external[i] <- mb$closure$Pout_external
      mass_df$Perr_external[i] <- mb$closure$Perr_external
      mass_df$Prel_external[i] <- mb$closure$Prel_external

      mass_df$Q_in_tanks[i] <- mb$inflow_tanks$Q_in_tanks
      mass_df$L_in_tanks[i] <- mb$inflow_tanks$L_in_tanks
      mass_df$C_in_tanks[i] <- mb$inflow_tanks$C_in_tanks

      mass_df$L_rain[i]   <- mb$inputs_external$L_rain
      mass_df$L_drydep[i] <- mb$inputs_external$L_drydep
      mass_df$L_seepin[i] <- mb$inputs_external$L_seepin

      mass_df$L_treated[i]        <- mb$outputs_external$L_treated
      mass_df$L_rel1[i]           <- mb$outputs_external$L_rel1
      mass_df$L_rel2[i]           <- mb$outputs_external$L_rel2
      mass_df$L_bypass[i]         <- mb$outputs_external$L_bypass
      mass_df$L_seep_discharge[i] <- mb$outputs_external$L_seep_discharge

      mass_df$L_seep_recycle_out[i] <- mb$transfers$L_seep_recycle_out

      mass_df$L_uptake[i]  <- mb$mechanisms$L_uptake
      mass_df$L_recycle[i] <- mb$mechanisms$L_recycle
      mass_df$L_sed[i]     <- mb$mechanisms$L_sed
      mass_df$L_direct[i]  <- mb$mechanisms$L_direct


    }

    if (return_steps) steps_store[[i]] <- res$steps
  }

  meta <- list(
    V_init = V_init,
    V_end = V,
    P_state_end = P_state,
    tanks = tanks,
    ppar = ppar,
    constants = constants,
    Nsteps = Nsteps,
    N_plant = N_plant
  )
  if (return_steps) meta$steps <- steps_store

  out <- list(
    results = results_df,
    budgets = list(
      water = water_df,
      mass  = mass_df
    ),
    meta = meta
  )

  class(out) <- c("dmsta_result", "list")
  out

}
