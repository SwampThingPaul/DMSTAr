#'
#' DMSTA hydrology engine (internal)
#'
#' Internal helpers implementing DMSTA hydrologic integration, including
#' within-day control depth interpolation, derivative evaluation, and RK4 stepping.
#'
#' @name internal_dmsta_hydro
#' @keywords internal
NULL

#'
#' Interpolate daily control depth within a day (internal)
#'
#' Internal helper to interpolate a daily control depth time series within a day
#' using one of three DMSTA interpolation conventions:
#' \itemize{
#'   \item `1`: inputs treated as daily means (no within-day variation)
#'   \item `2`: inputs treated as mid-day values (DMSTA default)
#'   \item `3`: inputs treated as end-of-day values
#' }
#'
#' Missing values are handled defensively:
#' \itemize{
#'   \item if today's `Zcontrol` is `NA`, it is treated as `0`
#'   \item if neighbor days are `NA`, they fall back to today's value
#' }
#' The interpolation fraction `delta` is clamped to \eqn{[0,1]}.
#'
#' @param Zcontrol1 Numeric scalar. Control depth for yesterday.
#' @param Zcontrol Numeric scalar. Control depth for today.
#' @param Zcontrol2 Numeric scalar. Control depth for tomorrow.
#' @param delta Numeric scalar in \eqn{[0,1]}. Fraction of day elapsed.
#' @param interp_option Integer. Interpolation convention (1, 2, or 3).
#'
#' @return Numeric scalar interpolated control depth for the within-day time point.
#'
#' @keywords internal
#' @rdname internal_dmsta_hydro

dmsta_interp_zcontrol <- function(
    Zcontrol1, # yesterday
    Zcontrol, # today
    Zcontrol2, # tomorrow
    delta, # delta in [0, 1]
    interp_option = 2L # 1=means, 2=mid-day (default), 3=end-of-day
) {
  # Zcontrol1 = yesterday, Zcontrol = today, Zcontrol2 = tomorrow
  # If today's control depth is missing, use 0
  if (is.na(Zcontrol))  Zcontrol  <- 0
  # If neighbors are missing, fall back to today's control depth
  if (is.na(Zcontrol1)) Zcontrol1 <- Zcontrol
  if (is.na(Zcontrol2)) Zcontrol2 <- Zcontrol

  # Be defensive about delta
  if (!is.finite(delta)) delta <- 0
  if (delta < 0) delta <- 0 else if (delta > 1) delta <- 1

  if (interp_option == 1L) {
    # Inputs treated as daily means: no within-day variation
    Zcontrol
  } else if (interp_option == 2L) {
    # Inputs treated as mid-day values (DMSTA Case 2)
    if (delta <= 0.5) Zcontrol1 + (Zcontrol  - Zcontrol1) * (delta + 0.5)
    else              Zcontrol  + (Zcontrol2 - Zcontrol ) * (delta - 0.5)
  } else if (interp_option == 3L) {
    # Inputs treated as end-of-day values (DMSTA Case 3)
    Zcontrol1 * (1 - delta) + Zcontrol * delta
  } else {
    # Fallback to DMSTA default (2)
    if (delta <= 0.5) Zcontrol1 + (Zcontrol  - Zcontrol1) * (delta + 0.5)
    else              Zcontrol  + (Zcontrol2 - Zcontrol ) * (delta - 0.5)
  }
}

#' Compute hydrologic fluxes and dV/dt for one derivative evaluation (internal)
#'
#' Internal helper corresponding to one derivative evaluation (DMSTA "DerivFlow")
#' used by the Runge-Kutta integrator. Given current volume and forcing inputs,
#' computes net fluxes, outflow, seepage, bypass, and the instantaneous rate of
#' change in volume (`dvdt`).
#'
#' Control depth is interpolated within the day via `dmsta_interp_zcontrol()`
#' using the RK stage timing defined by `step_index`, `nsteps`, and
#' `step_frac`.
#'
#' @param V Numeric scalar. Current volume (e.g., hm^3).
#' @param A_cell Numeric scalar. Cell area (e.g., km^2).
#' @param Qi Numeric scalar. Inflow rate (e.g., hm^3/day).
#' @param Rain Numeric scalar. Rainfall rate (e.g., m/day).
#' @param Et Numeric scalar. Evapotranspiration rate (e.g., m/day).
#' @param Zcontrol_t Numeric scalar. Today's control depth (m).
#' @param Zcontrol_t_minus Numeric scalar. Yesterday's control depth (m).
#' @param Zcontrol_t_plus Numeric scalar. Tomorrow's control depth (m).
#' @param params Named list of model parameters. Expected elements include
#'   (units in parentheses as used by the code):
#'   \describe{
#'     \item{Zmin}{minimum depth (cm; converted to m)}
#'     \item{Vmin}{optional minimum volume (if absent, computed as `Zmin*A_cell`)}
#'     \item{Q_a, Q_b}{weir/outflow coefficients}
#'     \item{Zweir}{weir crest (cm; converted to m)}
#'     \item{Q_zmin}{control depth offset (cm; converted to m)}
#'     \item{Qomax, Qimax}{optional max outflow / max inflow constraints}
#'     \item{Width}{weir width}
#'     \item{Bypass_elev}{bypass elevation (cm; converted to m)}
#'     \item{Seepout_Rate, Seepout_Elev}{seep-out parameters (rate; elev cm->m)}
#'     \item{Seepin_Rate, Seepin_Elev}{seep-in parameters (rate; elev cm->m)}
#'     \item{ShutdownET}{logical; reduce ET if minimum pool is violated}
#'     \item{force_Q_out}{logical; force outflow when no depth constraint}
#'     \item{Qr_0}{optional baseline release component}
#'     \item{interp_option}{optional integer; passed to `dmsta_interp_zcontrol`}
#'   }
#' @param step_index Integer. Current sub-step index in `1..nsteps`.
#' @param nsteps Integer. Total number of sub-steps per day.
#' @param step_frac Numeric. RK stage fraction (0, 0.5, 0.5, 1).
#' @param Ddt Numeric. Effective RK stage duration (days).
#' @param Qrelease Numeric scalar. Optional release rate (hm^3/day).
#' @param RecycleQ Numeric scalar. Optional recycle inflow (hm^3/day).
#' @param has_depth_constraint Logical. Whether a depth-control series is present.
#' @param ... Reserved for future extensions.
#'
#' @return A named list with instantaneous fluxes and diagnostics:
#' \describe{
#'   \item{dvdt}{Net `dV/dt`.}
#'   \item{Qo}{Outflow rate excluding `Qrelease`.}
#'   \item{Qout}{Total outflow rate including `Qrelease`.}
#'   \item{Etest}{Effective ET rate used after any shutdown logic.}
#'   \item{Seepout, Seepin}{Seepage rates.}
#'   \item{Bypass}{Bypass flow rate.}
#'   \item{Z}{Current depth estimate `V/A_cell`.}
#'   \item{Zcont}{Interpolated control depth used for this stage.}
#' }
#'
#' @keywords internal
#' @rdname internal_dmsta_hydro

dmsta_deriv_flow <- function(
    V,                   # current volume [hm3]
    A_cell,              # area [km2]
    Qi, Rain, Et,        # inflow [hm3/d], rain [m/d], ET [m/d]
    Zcontrol_t, Zcontrol_t_minus, Zcontrol_t_plus,  # control depths [m]
    params,             # list of parameters (see below)
    step_index, nsteps,  # current sub-step index (1..nsteps), total nsteps
    step_frac,          # fraction of step in RK4 staging: 0,0.5,0.5,1
    Ddt,                 # effective sub-step duration [day]; RK4 uses 0.5*Dt for stages 2-3
    Qrelease = 0,        # optional release [hm3/d]
    RecycleQ = 0,         # optional recycle inflow [hm3/d]
    has_depth_constraint = FALSE,  # <-- depth series present
    ...
) {

  ## parameters & defaults
  interp_option <- if (is.null(params$interp_option)) 2L else as.integer(params$interp_option)

  Zmin   <- params$Zmin / 100
  A      <- A_cell
  Vmin   <- if (is.null(params$Vmin)) Zmin * A else params$Vmin

  Q_a    <- params$Q_a
  Q_b    <- params$Q_b
  Zweir  <- params$Zweir / 100
  Q_zmin <- params$Q_zmin / 100
  Qomax  <- params$Qomax
  Qimax  <- params$Qimax
  Width  <- params$Width
  Bypass_elev <- params$Bypass_elev / 100

  Seepout_Rate <- params$Seepout_Rate
  Seepout_Elev <- params$Seepout_Elev / 100
  Seepin_Rate  <- params$Seepin_Rate
  Seepin_Elev  <- params$Seepin_Elev / 100

  ShutdownET  <- isTRUE(params$ShutdownET)
  force_Q_out <- isTRUE(params$force_Q_out)
  Qr_0        <- if (is.null(params$Qr_0)) 0 else params$Qr_0

  ##   derived state
  Z <- if (A > 0) V / A else 0

  delta <- ((step_index - 1) + step_frac) / nsteps
  Zcont <- dmsta_interp_zcontrol(
    Zcontrol_t_minus, Zcontrol_t, Zcontrol_t_plus,
    delta, interp_option
  )
  Zcont <- max(Zcont + Q_zmin, Zweir, Zmin)
  Vcontrol <- Zcont * A

  ##   seepage (always evaluated)
  Seepin  <- Seepin_Rate  * max(0, Seepin_Elev  - Z) * A
  Seepout <- Seepout_Rate * max(0, Z - Seepout_Elev) * A

  ##   atmosphere
  Etest <- Et
  AtmoS <- (Rain - Etest) * A

  ##   bypass logic (priority-based)
  if (Bypass_elev > 0 && Z > Bypass_elev) {
    Bypass <- Qi
  } else if (!is.null(Qimax) && Qimax > 0 && Qi > Qimax) {
    Bypass <- Qi - Qimax
  } else if (!is.null(Qomax) && Qomax < 0) {
    Bypass <- min(-Qomax, Qi)
  } else {
    Bypass <- 0
  }

  ##   net flux without overflow
  Qnet <- Qi + AtmoS - Seepout + Seepin + RecycleQ - Qrelease - Bypass
  Qo   <- 0

  ##   trial volume (no overflow)
  Vtrial <- V + Ddt * Qnet

  ##   minimum pool protection (VBA order)
  if (Vtrial < Vmin) {

    if (Seepout > 0) {
      Seepout <- max(Seepout + (Vtrial - Vmin) / Ddt, 0)
      Qnet    <- Qi + AtmoS - Seepout + Seepin + RecycleQ - Qrelease - Bypass
      Vtrial  <- V + Ddt * Qnet
    }

    if (Vtrial < Vmin && ShutdownET) {
      Etest <- max(Etest + (Vtrial - Vmin) / Ddt / A, 0)
      AtmoS <- (Rain - Etest) * A
      Qnet  <- Qi + AtmoS - Seepout + Seepin + RecycleQ - Qrelease - Bypass
      Vtrial <- V + Ddt * Qnet
    }

  } else if (force_Q_out && !has_depth_constraint) {

    Qo <- Qr_0

  } else if (!is.null(Q_a) && Q_a < 0) {

    Qo <- Qnet

  } else if (A > 0 && Z > Zcont) {

    if (!is.null(Q_a) && Q_a == 0) {
      Qo <- max((V - Vcontrol) / Ddt, 0)
    } else {
      Qo <- Q_a * Width * (Z - Zweir)^Q_b
    }

    if (!is.null(Qomax) && Qomax > 0 && Qo > Qomax) {
      Qo <- Qomax
    }
  }

  ##   final derivative (NO state clamping here)
  dvdt <- Qnet - Qo
  Qot  <- Qo + Qrelease

  ##   return fluxes & diagnostics
  list(
    dvdt    = dvdt,
    Qo      = Qo,
    Qout    = Qot,
    Etest   = Etest,
    Seepout = Seepout,
    Seepin  = Seepin,
    Bypass  = Bypass,
    Z       = Z,
    Zcont   = Zcont
  )
}

#' Advance volume by one RK4 sub-step (internal)
#'
#' Internal helper that performs one Runge-Kutta 4th order (RK4) sub-step for
#' hydrologic volume using `dmsta_deriv_flow()` for derivative evaluations.
#'
#' Returns the updated volume (with safety clamp to minimum volume) along with
#' stage-averaged fluxes for outflow, seepage, bypass, and ET.
#'
#' @param V Numeric scalar. Current volume at the start of the sub-step.
#' @param args_base Named list of arguments passed to `dmsta_deriv_flow()` via
#'   `do.call()`, typically containing `Qi`, `Rain`, `Et`,
#'   `Zcontrol_t`, neighbor control depths, `params`, `nsteps`,
#'   release components (`Qr0`, `Qr1`, `Qr2`), and flags.
#' @param step_index Integer. Sub-step index in `1..nsteps`.
#' @param Dt Numeric scalar. Sub-step duration in days (`1/nsteps`).
#'
#' @return A named list with elements:
#' \describe{
#'   \item{V_new}{Updated volume (clamped to minimum volume).}
#'   \item{Qout_ts}{Stage-averaged total outflow rate (including releases).}
#'   \item{SeepOut_ts, SeepIn_ts}{Stage-averaged seepage rates.}
#'   \item{Bypass_ts}{Stage-averaged bypass rate.}
#'   \item{Etest_ts}{Stage-averaged ET rate used.}
#'   \item{Q_treated_ts}{Portion of outflow assigned to the treated component.}
#'   \item{Q_rel1_ts, Q_rel2_ts}{Portions of outflow assigned to release components.}
#' }
#'
#' @keywords internal
#' @rdname internal_dmsta_hydro

dmsta_rk4_step <- function(V, args_base, step_index, Dt) {

  A_cell <- args_base$params$A_cell
  Zrel   <- if (is.null(args_base$params$Zrelease)) 0 else args_base$params$Zrelease
  Vmin   <- if (is.null(args_base$params$Vmin)) args_base$params$Zmin/100 * A_cell else args_base$params$Vmin

  Qr0 <- if (is.null(args_base$Qr0)) 0 else args_base$Qr0
  Qr1 <- if (is.null(args_base$Qr1)) 0 else args_base$Qr1
  Qr2 <- if (is.null(args_base$Qr2)) 0 else args_base$Qr2

  ##  RK stage 1
  Qrelease1 <- if (A_cell > 0 && V <= Zrel * A_cell) 0 else (Qr1 + Qr2)

  s1 <- do.call(dmsta_deriv_flow,
                c(list(V = V, step_frac = 0.0, Ddt = Dt,
                       step_index = step_index,
                       Qrelease = Qrelease1),
                  args_base))

  V2 <- V + s1$dvdt * 0.5 * Dt

  ##  RK stage 2
  Qrelease2 <- if (A_cell > 0 && V2 <= Zrel * A_cell) 0 else (Qr1 + Qr2)

  s2 <- do.call(dmsta_deriv_flow,
                c(list(V = V2, step_frac = 0.5, Ddt = Dt,
                       step_index = step_index,
                       Qrelease = Qrelease2),
                  args_base))

  V3 <- V + s2$dvdt * 0.5 * Dt

  ##  RK stage 3
  Qrelease3 <- if (A_cell > 0 && V3 <= Zrel * A_cell) 0 else (Qr1 + Qr2)

  s3 <- do.call(dmsta_deriv_flow,
                c(list(V = V3, step_frac = 0.5, Ddt = Dt,
                       step_index = step_index,
                       Qrelease = Qrelease3),
                  args_base))

  V4 <- V + s3$dvdt * Dt

  ##  RK stage 4
  Qrelease4 <- if (A_cell > 0 && V4 <= Zrel * A_cell) 0 else (Qr1 + Qr2)

  s4 <- do.call(dmsta_deriv_flow,
                c(list(V = V4, step_frac = 1.0, Ddt = Dt,
                       step_index = step_index,
                       Qrelease = Qrelease4),
                  args_base))

  ##  RK combine
  dv <- (s1$dvdt + 2*s2$dvdt + 2*s3$dvdt + s4$dvdt) / 6

  qout <- (s1$Qout + 2*s2$Qout + 2*s3$Qout + s4$Qout) / 6
  seepO <- (s1$Seepout + 2*s2$Seepout + 2*s3$Seepout + s4$Seepout) / 6
  seepI <- (s1$Seepin  + 2*s2$Seepin  + 2*s3$Seepin  + s4$Seepin ) / 6
  byps  <- (s1$Bypass  + 2*s2$Bypass  + 2*s3$Bypass  + s4$Bypass ) / 6
  Etest_m <- (s1$Etest + 2*s2$Etest + 2*s3$Etest + s4$Etest) / 6

  ##  state update + safety clamp (THIS is the right place)
  V_new <- max(V + dv * Dt, Vmin)

  ##  release splitting (use instantaneous end-of-step outflow)
  q_end <- s4$Qout
  Sspec <- Qr0 + Qr1 + Qr2

  if (q_end > Sspec && q_end > 0) {
    f1 <- Qr1 / q_end
    f2 <- Qr2 / q_end
  } else if (Sspec > 0) {
    f1 <- Qr1 / Sspec
    f2 <- Qr2 / Sspec
  } else {
    f1 <- 0; f2 <- 0
  }

  f0 <- max(0, 1 - f1 - f2)

  list(
    V_new = V_new,
    Qout_ts = qout,
    SeepOut_ts = seepO,
    SeepIn_ts = seepI,
    Bypass_ts = byps,
    Etest_ts = Etest_m,
    Q_treated_ts = qout * f0,
    Q_rel1_ts = qout * f1,
    Q_rel2_ts = qout * f2
  )
}

#' Integrate one day of DMSTA hydrology using RK4 (internal)
#'
#' Internal helper that integrates hydrologic volume and fluxes over one day
#' using `Nsteps` RK4 sub-steps. Returns daily totals/diagnostics and
#' (optionally) per-substep details in `meta$steps`.
#'
#' @param V Numeric scalar. Starting volume for the day (e.g., hm^3).
#' @param inputs Named list of day-specific forcings with elements:
#'   \describe{
#'     \item{Qi}{Inflow rate (hm^3/day).}
#'     \item{Rain}{Rainfall rate (m/day).}
#'     \item{Et}{ET rate (m/day).}
#'     \item{Zcontrol}{Today's control depth (m).}
#'     \item{Zcontrol_prev}{Yesterday's control depth (m); optional.}
#'     \item{Zcontrol_next}{Tomorrow's control depth (m); optional.}
#'     \item{RecycleQ}{Recycle inflow (hm^3/day); optional.}
#'     \item{Qr0, Qr1, Qr2}{Release components (hm^3/day); optional.}
#'     \item{has_depth_constraint}{Logical; optional.}
#'   }
#' @param params Named list of model parameters. Must include `A_cell` and
#'   `Zmin` at minimum; see `dmsta_deriv_flow()` for additional fields.
#' @param Nsteps Integer. Number of RK4 sub-steps per day.
#'
#' @return A structured list with components:
#' \describe{
#'   \item{results}{Named list of daily totals and end-of-day state (e.g., `V_end`, `Qout`, seepage, bypass).}
#'   \item{budgets}{Named lists for water budget diagnostics (mass budget is `NULL` for hydrology-only).}
#'   \item{meta}{Inputs used, parameters used, and per-substep records.}
#' }
#'
#' @keywords internal
#' @rdname internal_dmsta_hydro

dmsta_flow_day <- function(
    V,                      # starting volume [hm3]
    inputs,                  # list: Qi, Rain, Et, Zcontrol (today), Zcontrol_prev, Zcontrol_next,Qrelease (optional), RecycleQ (optional)
    params,                  # list of parameters (see below)
    Nsteps = 4             # steps per day (same meaning as DMSTA 'steps')
) {

  Dt <- 1 / Nsteps
  V_start <- V

  params2 <- params
  if (is.null(params2$interp_option)) params2$interp_option <- 2L
  if (is.null(params2$Qin_Frac)) params2$Qin_Frac <- 1.0

  if (is.null(inputs$Zcontrol_prev)) inputs$Zcontrol_prev <- inputs$Zcontrol
  if (is.null(inputs$Zcontrol_next)) inputs$Zcontrol_next <- inputs$Zcontrol
  if (is.null(inputs$RecycleQ)) inputs$RecycleQ <- 0

  Qi_eff <- params2$Qin_Frac * inputs$Qi

  args_base <- list(
    A_cell = params2$A_cell,
    Qi = Qi_eff,
    Rain = inputs$Rain,
    Et = inputs$Et,
    Zcontrol_t = inputs$Zcontrol,
    Zcontrol_t_minus = inputs$Zcontrol_prev,
    Zcontrol_t_plus = inputs$Zcontrol_next,
    params = params2,
    nsteps = Nsteps,
    RecycleQ = inputs$RecycleQ,
    Qr0 = if (is.null(inputs$Qr0)) 0 else inputs$Qr0,
    Qr1 = if (is.null(inputs$Qr1)) 0 else inputs$Qr1,
    Qr2 = if (is.null(inputs$Qr2)) 0 else inputs$Qr2,
    has_depth_constraint = isTRUE(inputs$has_depth_constraint)
  )

  # per-step records
  step_out <- vector("list", Nsteps)

  for (k in seq_len(Nsteps)) {
    Vo <- V
    step_res <- dmsta_rk4_step(V, args_base, step_index = k, Dt = Dt)
    V <- step_res$V_new

    step_out[[k]] <- list(
      step = k,
      Vo = Vo,
      V  = V,
      Dt = Dt,
      Qi_eff = Qi_eff,
      Qout = step_res$Qout_ts,
      Q_treated = step_res$Q_treated_ts,
      Q_rel1 = step_res$Q_rel1_ts,
      Q_rel2 = step_res$Q_rel2_ts,
      SeepOut = step_res$SeepOut_ts,
      SeepIn  = step_res$SeepIn_ts,
      Bypass  = step_res$Bypass_ts,
      Etest   = step_res$Etest_ts
    )
  }

  # daily aggregates (same as your dmsta_flow_day)
  Qout_day <- sum(vapply(step_out, function(x) x$Qout * x$Dt, 0.0))
  Q_treated_day <- sum(vapply(step_out, function(x) x$Q_treated * x$Dt, 0.0))
  Q_rel1_day    <- sum(vapply(step_out, function(x) x$Q_rel1    * x$Dt, 0.0))
  Q_rel2_day    <- sum(vapply(step_out, function(x) x$Q_rel2    * x$Dt, 0.0))

  SeepOut_day <- sum(vapply(step_out, function(x) x$SeepOut * x$Dt, 0.0))
  SeepIn_day  <- sum(vapply(step_out, function(x) x$SeepIn  * x$Dt, 0.0))
  Bypass_day  <- sum(vapply(step_out, function(x) x$Bypass  * x$Dt, 0.0))

  EtVol_day   <- sum(vapply(step_out, function(x) x$Etest * params2$A_cell * x$Dt, 0.0))
  RainVol_day <- inputs$Rain * params2$A_cell
  NetAtmo_day <- RainVol_day - EtVol_day

  RecycleVol_day <- inputs$RecycleQ*1.0

  # Water budget (DMSTA-style)
  WB_in  <- Qi_eff + RainVol_day + SeepIn_day + RecycleVol_day
  WB_out <- Qout_day + SeepOut_day + EtVol_day + Bypass_day
  WB_err <- (V - V_start) - (WB_in - WB_out)
  WB_rel <- WB_err / max(1e-12, max(WB_in, WB_out))

  # structured outputs
  results <- list(
    V_end = V,
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
    NetAtmo = NetAtmo_day
  )

  water_budget <- list(
    # main diagnostics
    WB_in = WB_in,
    WB_out = WB_out,
    WB_err = WB_err,
    WB_rel = WB_rel,

    # atmospheric terms (often treated as budget components)
    RainVol = RainVol_day,
    EtVol   = EtVol_day,
    NetAtmo = NetAtmo_day,

    # optional component breakdowns (very useful for debugging)
    in_components  = list(
      Qi_eff   = Qi_eff,
      RainVol  = RainVol_day,
      SeepIn   = SeepIn_day,
      Recycle  = RecycleVol_day
    ),
    out_components = list(
      Qout     = Qout_day,
      SeepOut  = SeepOut_day,
      EtVol    = EtVol_day,
      Bypass   = Bypass_day
    ),
    dV = V - V_start
  )
  mass_budget <- NULL  # hydrology-only

  meta <- list(
    V_start = V_start,
    V_end   = V,
    Nsteps  = Nsteps,
    Dt      = Dt,
    Qi_eff  = Qi_eff,
    params_used = params2,
    # useful for debugging / parity checks
    inputs_used = inputs,
    steps = step_out
  )

  list(
    results = results,
    budgets = list(
      water = water_budget,
      mass  = mass_budget
    ),
    meta = meta
  )
}

#' Integrate one day of DMSTA hydrology and return sub-step means (internal)
#'
#' Similar to `dmsta_flow_day()` but returns the per-substep mean rates
#' explicitly in the `steps` element (useful for coupling to constituent
#' transport/biogeochemistry that needs sub-daily hydrology).
#'
#' @param V Numeric scalar. Starting volume for the day as hm3
#' @param inputs Named list of forcings; see `dmsta_flow_day()`.
#' @param params Named list of parameters; see `dmsta_flow_day()`.
#' @param Nsteps Integer. Number of sub-steps per day.
#'
#' @return A named list including:
#' \describe{
#'   \item{V_end}{End-of-day volume.}
#'   \item{steps}{List of length `Nsteps` with per-substep mean rates and states.}
#'   \item{Qout, SeepOut, SeepIn, Bypass, RainVol, EtVol, NetAtmo}{Daily totals.}
#'   \item{WB_in, WB_out, WB_err, WB_rel}{Water balance diagnostics.}
#' }
#'
#' @keywords internal
#' @rdname internal_dmsta_hydro

dmsta_flow_day_steps <- function(
    V,          # starting volume [hm3]
    inputs,     # list: Qi, Rain, Et, Zcontrol (today), Zcontrol_prev, Zcontrol_next, optional: Qr0, Qr1, Qr2, RecycleQ
    params,     # parameter list
    Nsteps = 4L # steps per day (DMSTA 'steps')
) {

  Dt <- 1 / Nsteps
  V_start <- V

  params2 <- params
  if (is.null(params2$interp_option)) params2$interp_option <- 2L
  if (is.null(params2$Qin_Frac)) params2$Qin_Frac <- 1.0

  # fill optional input fields
  if (is.null(inputs$Zcontrol_prev)) inputs$Zcontrol_prev <- inputs$Zcontrol
  if (is.null(inputs$Zcontrol_next)) inputs$Zcontrol_next <- inputs$Zcontrol
  if (is.null(inputs$RecycleQ)) inputs$RecycleQ <- 0

  # scale inflow by cell inflow fraction
  Qi_eff <- params2$Qin_Frac * inputs$Qi

  # base args passed into dmsta_rk4_step -> dmsta_deriv_flow
  args_base <- list(
    A_cell = params2$A_cell,
    Qi = Qi_eff,
    Rain = inputs$Rain,
    Et = inputs$Et,
    Zcontrol_t = inputs$Zcontrol,
    Zcontrol_t_minus = inputs$Zcontrol_prev,
    Zcontrol_t_plus = inputs$Zcontrol_next,
    params = params2,
    nsteps = Nsteps,
    RecycleQ = inputs$RecycleQ,
    Qr0 = if (is.null(inputs$Qr0)) 0 else inputs$Qr0,
    Qr1 = if (is.null(inputs$Qr1)) 0 else inputs$Qr1,
    Qr2 = if (is.null(inputs$Qr2)) 0 else inputs$Qr2,
    has_depth_constraint = isTRUE(inputs$has_depth_constraint)
  )

  # store per-substep results
  step_out <- vector("list", Nsteps)

  for (k in seq_len(Nsteps)) {
    Vo <- V
    step_res <- dmsta_rk4_step(V, args_base, step_index = k, Dt = Dt)
    V <- step_res$V_new

    step_out[[k]] <- list(
      step = k,
      Vo = Vo,
      V  = V,
      Dt = Dt,
      Qi_eff = Qi_eff,

      # mean rates for this substep (hm3/day or m/day for Etest)
      Qout = step_res$Qout_ts,
      Q_treated = step_res$Q_treated_ts,
      Q_rel1 = step_res$Q_rel1_ts,
      Q_rel2 = step_res$Q_rel2_ts,
      SeepOut = step_res$SeepOut_ts,
      SeepIn  = step_res$SeepIn_ts,
      Bypass  = step_res$Bypass_ts,
      Etest   = step_res$Etest_ts
    )
  }

  # daily aggregates (integrate substep mean rates over the day)
  Qout_day      <- sum(vapply(step_out, function(x) x$Qout      * x$Dt, 0.0))
  Q_treated_day <- sum(vapply(step_out, function(x) x$Q_treated * x$Dt, 0.0))
  Q_rel1_day    <- sum(vapply(step_out, function(x) x$Q_rel1    * x$Dt, 0.0))
  Q_rel2_day    <- sum(vapply(step_out, function(x) x$Q_rel2    * x$Dt, 0.0))

  SeepOut_day <- sum(vapply(step_out, function(x) x$SeepOut * x$Dt, 0.0))
  SeepIn_day  <- sum(vapply(step_out, function(x) x$SeepIn  * x$Dt, 0.0))
  Bypass_day  <- sum(vapply(step_out, function(x) x$Bypass  * x$Dt, 0.0))

  # ET volume is Etest(m/day) * A_cell(km2) -> hm3/day; integrate over day
  EtVol_day <- sum(vapply(step_out, function(x) x$Etest * params2$A_cell * x$Dt, 0.0))

  # Rain volume is constant in DMSTA for the day: Rain(m/day) * A_cell(km2) = hm3/day; for 1 day => hm3
  RainVol_day <- inputs$Rain * params2$A_cell
  NetAtmo_day <- RainVol_day - EtVol_day

  # recycle volume (hm3) for the day
  RecycleVol_day <- inputs$RecycleQ * 1.0

  # water balance diagnostics
  WB_in  <- Qi_eff + RainVol_day + SeepIn_day + RecycleVol_day
  WB_out <- Qout_day + SeepOut_day + EtVol_day + Bypass_day
  WB_err <- (V - V_start) - (WB_in - WB_out)
  WB_rel <- WB_err / max(1e-12, max(WB_in, WB_out))

  list(
    # end-of-day state
    V_end = V,

    # **key addition** for P coupling
    steps = step_out,

    # daily totals
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

    # diagnostics
    WB_in = WB_in,
    WB_out = WB_out,
    WB_err = WB_err,
    WB_rel = WB_rel
  )
}




#' Run a DMSTA hydrology simulation over a time series
#'
#' Integrates DMSTA hydrology (volume and fluxes) over multiple days using
#' `dmsta_flow_day()` and returns daily totals and diagnostics.
#'
#' The input `series` must contain daily forcings (inflow, rain, ET) and a
#' daily control depth. Neighbor-day control depths are derived using
#' `neighbors_zcontrol()` (must be available in the package namespace).
#'
#' @param V_init Numeric scalar. Initial volume at the start of the series (e.g., hm^3).
#' @param series A `data.frame` with at least the following columns:
#'   \describe{
#'     \item{Date}{Date column coercible via `as.Date()`.}
#'     \item{Qi}{Daily inflow rate (hm^3/day).}
#'     \item{Rain}{Daily rainfall rate (m/day).}
#'     \item{Et}{Daily ET rate (m/day).}
#'     \item{Zcontrol}{Daily control depth (m).}
#'   }
#' @param params Named list of model parameters. Must include `A_cell`
#'   (area) and `Zmin` (cm, converted to m) at minimum; see
#'   `dmsta_deriv_flow()` for additional parameters.
#' @param Nsteps Integer. Number of RK4 sub-steps per day (default `4L`).
#'
#' @details
#' The function computes daily results and a water budget diagnostic including
#' `WB_in`, `WB_out`, and closure error `WB_err`.
#'
#' **Note:** In the current implementation, the final assignment
#' `results_df$Qin <- series$Qin` expects a column named `Qin`, but the
#' required input column is `Qi`. Consider changing this to
#' `results_df$Qin <- series$Qi` (or require `Qin` instead).
#'
#' @return An object of class `"dmsta_hydro_result"` (a list) with elements:
#' \describe{
#'   \item{results}{A data.frame of daily totals and end-of-day state variables.}
#'   \item{budgets}{A list containing a water budget data.frame and `mass = NULL`.}
#'   \item{meta}{Metadata including `V_init`, `V_end`, `Nsteps`, and `A_cell`.}
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
#'  RecycleQ = 0
#')
#'
#'  V_init <- (cm_to_m(params$Zinit) * params$A_cell)
#'  out <- dmsta_flow_series(V_init, series, params, Nsteps = 4)
#'  out$results
#'  }
#'
#' @export
dmsta_flow_series <- function(V_init, series, params, Nsteps = 4L){

  n <- nrow(series)
  if (n < 1) stop("series has zero rows.")
  if (!("Date" %in% names(series))) stop("series must include a 'Date' column.")
  if (!all(c("Qi","Rain","Et","Zcontrol") %in% names(series))) {
    stop("series must include columns: Qi, Rain, Et, Zcontrol")
  }

  # --- results data.frame (primary outputs) ---
  results_df <- data.frame(
    Date = as.Date(series$Date),
    Qin = series$Qi,
    Qout = NA_real_,
    Q_treated = NA_real_,
    Q_rel1 = NA_real_,
    Q_rel2 = NA_real_,
    SeepOut = NA_real_,
    SeepIn = NA_real_,
    Bypass = NA_real_,
    V_end = NA_real_,
    Z_end = NA_real_,
    RainVol = NA_real_,
    EtVol = NA_real_,
    NetAtmo = NA_real_,
    stringsAsFactors = FALSE
  )

  # water budgets data.frame (separate)
  water_df <- data.frame(
    Date = as.Date(series$Date),
    WB_in = NA_real_,
    WB_out = NA_real_,
    WB_err = NA_real_,
    WB_rel = NA_real_,
    stringsAsFactors = FALSE
  )


  ##  persistent hydrologic state (Module1 V)
  V <- V_init

  for (i in seq_len(n)) {

    in_i <- as.list(series[i, ])

    ##  control-depth neighbors (Module1 semantics)
    nz <- neighbors_zcontrol(i, series$Zcontrol)
    in_i$Zcontrol      <- nz$today
    in_i$Zcontrol_prev <- nz$prev_day
    in_i$Zcontrol_next <- nz$nxt

    ##  day-specific depth constraint (HydroIndex logic)
    in_i$has_depth_constraint <-
      is.finite(in_i$Zcontrol) && in_i$Zcontrol != 0

    ##  run one DMSTA hydrology day
    day_res <- dmsta_flow_day(
      V      = V,
      inputs = in_i,
      params = params,
      Nsteps = Nsteps
    )

    ##  advance state
    V <- day_res$results$V_end

    # Store results
    results_df$Qin[i]      <- day_res$results$Qin
    results_df$Qout[i]      <- day_res$results$Qout
    results_df$Q_treated[i] <- day_res$results$Q_treated
    results_df$Q_rel1[i]    <- day_res$results$Q_rel1
    results_df$Q_rel2[i]    <- day_res$results$Q_rel2
    results_df$SeepOut[i]   <- day_res$results$SeepOut
    results_df$SeepIn[i]    <- day_res$results$SeepIn
    results_df$Bypass[i]    <- day_res$results$Bypass
    results_df$V_end[i]     <- V
    results_df$Z_end[i]     <- if (params$A_cell > 0) V / params$A_cell else NA_real_

    # Store water budget diagnostics
    water_df$RainVol[i] <- day_res$budgets$water$RainVol
    water_df$EtVol[i]   <- day_res$budgets$water$EtVol
    water_df$NetAtmo[i] <- day_res$budgets$water$NetAtmo
    water_df$WB_in[i]   <- day_res$budgets$water$WB_in
    water_df$WB_out[i]  <- day_res$budgets$water$WB_out
    water_df$WB_err[i]  <- day_res$budgets$water$WB_err
    water_df$WB_rel[i]  <- day_res$budgets$water$WB_rel
  }

  # results_df$Qin <- series$Qin
  # Structured output
  out <- list(
    results = results_df,
    budgets = list(
      water = water_df,
      mass  = NULL
    ),
    meta = list(
      V_init = V_init,
      V_end  = V,
      Nsteps = Nsteps,
      A_cell = params$A_cell
    )
  )
  class(out) <- c("dmsta_hydro_result", "list")
  out
}
