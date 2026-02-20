#' LPWEM internal helper functions
#'
#' Internal building blocks used by the LPWEM (Lake Phosphorus Wetland
#' Ecosystem Model) implementation for transport and kinetic updates across
#' a tanks-in-series wetland representation.
#'
#' These functions are not intended to be called directly by end users.
#'
#' @name lpwem-internal
#' @keywords internal
NULL

#' Convert an areal flux to a mass rate (kg/day)
#'
#' Converts an areal flux expressed in \eqn{mg m^{-2} d^{-1}} to a mass rate
#' in \eqn{kg d^{-1}} by multiplying by area in \eqn{km^2}.
#'
#' @details
#' The conversion is implicit because:
#' \deqn{(mg/m^2/day) * (km^2) = (mg/m^2/day) * (10^6 m^2) = 10^6 mg/day = 1 kg/day.}
#'
#'
#' @param rate_mgm2d Numeric areal rate in \eqn{mg m^{-2} d^{-1}}.
#' @param A_km2 Numeric area in \eqn{km^2}.
#'
#' @return Numeric mass rate in \eqn{kg d^{-1}}.
#'
#' @rdname lpwem-internal
#' @keywords internal
lpwem_areal_to_kgd <- function(rate_mgm2d, A_km2) {
  rate_mgm2d <- safe_num_vec(rate_mgm2d, 0)
  A_km2 <- safe_num_vec(A_km2, 0)
  # mg/m2-day * km2 == kg/day numerically
  rate_mgm2d * A_km2
}

#' Build an equal-area tanks-in-series structure
#'
#' Constructs an equal-area tanks-in-series (TIS) representation for a wetland
#' cell, where each tank has area \eqn{A_i = A_{cell}/N}.
#'
#' @param A_cell_km2 Total wetland surface area in \eqn{km^2}.
#' @param N Integer number of tanks. Must be \eqn{\ge 1}.
#'
#' @return A list with:
#' \itemize{
#'   \item `N`: integer number of tanks
#'   \item `A_i`: numeric vector of length `N` of tank areas (\eqn{km^2})
#' }
#'
#' @examples
#' lpwem_build_tis(12, N = 6)
#'
#' @rdname lpwem-internal
#' @keywords internal
lpwem_build_tis <- function(A_cell_km2, N = 6L) {
  N <- as.integer(N)
  if (N < 1) stop("N must be >= 1")
  A_i <- rep(A_cell_km2 / N, N)
  list(N = N, A_i = A_i)
}


#' Compute LPWEM return-flux profile across tanks
#'
#' Computes the LPWEM return flux \eqn{r^*_i} (\eqn{mg m^{-2} d^{-1}}) across a
#' tanks-in-series chain, linearly interpolated from the first tank value
#' `r1_mgm2d` to the last tank value `rN_mgm2d`.
#'
#' @param N Integer number of tanks.
#' @param r1_mgm2d Return flux at the first tank (\eqn{mg m^{-2} d^{-1}}).
#' @param rN_mgm2d Return flux at the last tank (\eqn{mg m^{-2} d^{-1}}).
#'
#' @return Numeric vector of length `N` with the return flux for each tank.
#'
#' @examples
#' lpwem_rstar_vec(N = 4, r1_mgm2d = 10, rN_mgm2d = 2)
#'
#' @rdname lpwem-internal
#' @keywords internal
lpwem_rstar_vec <- function(N, r1_mgm2d, rN_mgm2d) {
  if (N == 1) return(rep(r1_mgm2d, 1))
  frac <- (0:(N-1)) / (N-1)
  r1_mgm2d + frac * (rN_mgm2d - r1_mgm2d)
}

#' Compute LPWEM uptake flux across tanks
#'
#' Computes the LPWEM areal uptake flux \eqn{J_i} (\eqn{mg m^{-2} d^{-1}}) from
#' tank concentrations and hydrologic forcing using:
#' \deqn{J = k_0 \left(1 + \alpha \frac{Q}{W}\right) TP^2}
#' where \eqn{Q} is in \eqn{m^3 d^{-1}}, \eqn{W} is in \eqn{m}, and \eqn{TP} is
#' in \eqn{ppb} (treated as \eqn{mg m^{-3}}). Unit consistency is carried by
#' \eqn{k_0}.
#'
#' @param TP_ppb Numeric vector of tank TP concentrations (ppb).
#' @param Q_hm3d Numeric vector of tank flows (\eqn{hm^3 d^{-1}}), typically the
#'   outflow from each tank.
#' @param W_m Numeric wetland width parameter (\eqn{m}).
#' @param k0 Numeric uptake coefficient \eqn{k_0}.
#' @param alpha Numeric coefficient \eqn{\alpha} applied to \eqn{Q/W}.
#'
#' @return Numeric vector of uptake flux \eqn{J_i} in \eqn{mg m^{-2} d^{-1}}.
#'
#' @rdname lpwem-internal
#' @keywords internal
lpwem_J_vec <- function(TP_ppb, Q_hm3d, W_m, k0, alpha) {
  TP_ppb  <- safe_num_vec(TP_ppb, 0)
  Q_hm3d  <- safe_num_vec(Q_hm3d, 0)

  W_m   <- safe_num(W_m, 1)
  k0    <- safe_num(k0, 0)
  alpha <- safe_num(alpha, 0)

  Q_m3d <- Q_hm3d * 1e6
  vel_term <- 1 + alpha * (Q_m3d / max(1e-12, W_m))

  k0 * vel_term * (pmax(TP_ppb, 0)^2)
}


#' Interpolate flow along tanks-in-series
#'
#' Constructs the sequence \eqn{Q_i} for \eqn{i = 0..N} (length \eqn{N+1}) by
#' linear interpolation between wetland inflow \eqn{Q_0} and wetland outflow
#' \eqn{Q_N}.
#'
#' @param Q0_hm3d Wetland inflow \eqn{Q_0} in \eqn{hm^3 d^{-1}}.
#' @param QN_hm3d Wetland outflow \eqn{Q_N} in \eqn{hm^3 d^{-1}}.
#' @param N Integer number of tanks.
#'
#' @return Numeric vector of length `N + 1` giving `Q_0, ..., Q_N`
#' in \eqn{hm^3 d^{-1}}.
#'
#' @examples
#' lpwem_Q_along_tanks(Q0_hm3d = 2, QN_hm3d = 1, N = 4)
#'
#' @rdname lpwem-internal
#' @keywords internal
lpwem_Q_along_tanks <- function(Q0_hm3d, QN_hm3d, N) {
  i <- 0:N
  Q0_hm3d + (QN_hm3d - Q0_hm3d) * (i / N)
}


#' One explicit Euler micro-step for LPWEM TP dynamics
#'
#' Advances the tank TP concentration state one micro time step using an explicit
#' Euler update of the mass-balance formulation across a tanks-in-series chain.
#'
#' @details
#' The implementation updates per-tank mass \eqn{M_i} and converts back to
#' concentration. Water depth is approximated as:
#' \deqn{h = V/A}
#' with \eqn{V} in \eqn{hm^3} and \eqn{A} in \eqn{km^2}, yielding \eqn{h} in
#' meters because \eqn{1\,m * 1\,km^2 = 1\,hm^3}.
#'
#' Flow is interpolated along tanks from \eqn{Q_0} to \eqn{Q_N}, and the return
#' flux \eqn{r^*} is interpolated linearly from the first to last tank values.
#'
#' @param TP Numeric vector of length `tis$N` of tank TP concentrations (ppb).
#' @param Vavg_hm3 Step-average wetland volume (\eqn{hm^3}).
#' @param Q0_hm3d Wetland inflow rate (\eqn{hm^3 d^{-1}}).
#' @param QN_hm3d Wetland outflow rate (\eqn{hm^3 d^{-1}}).
#' @param tis Tanks-in-series structure as returned by [lpwem_build_tis()].
#' @param pars List of LPWEM parameters; must include `Cin_ppb`, `W_m`,
#'   `k0`, `alpha`, `r1_mgm2d`, `rN_mgm2d`.
#' @param dt_day Micro-step duration in days.
#'
#' @return A list with:
#' \itemize{
#'   \item `TP_new`: updated TP vector (ppb)
#'   \item `diag`: diagnostics including `h_m`, `Q_along`,
#'     `J_mgm2d`, `r_mgm2d`, and per-tank mass-balance terms
#' }
#'
#' @rdname lpwem-internal
#' @keywords internal
lpwem_tp_euler_step <- function(TP, Vavg_hm3, Q0_hm3d, QN_hm3d, tis, pars, dt_day) {
  N <- tis$N
  A_i <- tis$A_i
  A_cell <- sum(A_i)

  h_m <- if (A_cell > 0) Vavg_hm3 / A_cell else 0
  Qvec <- lpwem_Q_along_tanks(Q0_hm3d, QN_hm3d, N)

  rvec <- lpwem_rstar_vec(N, pars$r1_mgm2d, pars$rN_mgm2d)

  Jvec <- lpwem_J_vec(TP_ppb = TP, Q_hm3d = Qvec[2:(N+1)],
                      W_m = pars$W_m, k0 = pars$k0, alpha = pars$alpha)

  V_i <- h_m * A_i
  M <- TP * V_i

  TP_in <- c(pars$Cin_ppb, TP[1:(N-1)])

  Lin  <- Qvec[1:N] * TP_in
  Lout <- Qvec[2:(N+1)] * TP

  Uptake_kgd <- lpwem_areal_to_kgd(Jvec, A_i)
  Return_kgd <- lpwem_areal_to_kgd(rvec, A_i)

  dMdt <- Lin - Lout - Uptake_kgd + Return_kgd

  M_new <- pmax(0, M + dMdt * dt_day)
  TP_new <- ifelse(V_i > 0, M_new / V_i, 0)

  list(
    TP_new = TP_new,
    diag = list(
      h_m = h_m,
      Q_along = Qvec,
      J_mgm2d = Jvec,
      r_mgm2d = rvec,
      Uptake_kgd = Uptake_kgd,
      Return_kgd = Return_kgd,
      Lin_kgd = Lin,
      Lout_kgd = Lout
    )
  )
}

#' Run one model day with DMSTA hydrology and LPWEM TP kinetics
#'
#' Runs a single day consisting of DMSTA hydrology substeps coupled with LPWEM
#' tank-chain phosphorus dynamics. Hydrology is advanced using
#' `DMSTAr:::dmsta_flow_day_steps()`, and TP is updated in each hydrology
#' substep using multiple LPWEM micro-steps (explicit Euler).
#'
#' @param V Starting wetland volume at the beginning of the day (\eqn{hm^3}).
#' @param TP_init Optional numeric vector of initial TP concentrations (ppb) in
#'   each tank (length `pars_lpwem$N`). If `NULL`, all tanks start at 20.
#' @param inputs List of day inputs. Must include at least:
#' \itemize{
#'   \item `Date` (date-like)
#'   \item `Qi` inflow (\eqn{hm^3 d^{-1}})
#'   \item `Ci` inflow TP concentration (ppb)
#'   \item `Rain` rainfall (units consistent with DMSTA hydrology)
#'   \item `Et` evapotranspiration (units consistent with DMSTA hydrology)
#'   \item `Zcontrol` stage/depth control (units consistent with DMSTA hydrology)
#' }
#' Additional optional hydrology fields (e.g., recycle flows) may be included.
#' @param params_hydro List of DMSTA hydrology parameters (e.g., `A_cell`).
#' @param pars_lpwem List of LPWEM parameters. Must include:
#' \itemize{
#'   \item `N` number of tanks
#'   \item `W_m`, `k0`, `alpha`
#'   \item `r1_mgm2d`, `rN_mgm2d`
#'   \item `Np_per_day` micro-steps per day (typical value ~20)
#' }
#' @param Nsteps_hydro Integer number of hydrology substeps per day.
#'
#' @return A list with:
#' \itemize{
#'   \item `results`: end-of-day state and summary fluxes (volume, TP outflow,
#'     flow-weighted outflow concentration, mass exports, uptake/return totals)
#'   \item `budgets`: water and mass budget summaries
#'   \item `meta`: run metadata including diagnostics per hydrology step
#' }
#'
#'
#' @references
#' Juston, J. M., & Kadlec, R. H. (2019).
#' Data-driven modeling of phosphorus (P) dynamics in low-P stormwater wetlands.
#' *Environmental Modelling & Software*, 118, 226–240.
#' https://doi.org/10.1016/j.envsoft.2019.05.002
#'
#'
#' @export
lpwem_day <- function(
    V,
    TP_init = NULL,
    inputs,
    params_hydro,
    pars_lpwem,
    Nsteps_hydro = 4L
) {
  # 1) Hydrology with your DMSTA function (returns per-substep step list)
  hyd <- DMSTAr:::dmsta_flow_day_steps(V,
                                       inputs,
                                       params_hydro,
                                       Nsteps = Nsteps_hydro)

  steps <- hyd$steps
  A_cell <- params_hydro$A_cell

  # 2) LPWEM TIS structure
  tis <- lpwem_build_tis(A_cell_km2 = A_cell, N = pars_lpwem$N)

  # 3) Init TP state
  if (is.null(TP_init)) TP <- rep(20, tis$N) else {
    if (length(TP_init) != tis$N) stop("TP_init length must equal N tanks.")
    TP <- TP_init
  }

  # Daily accumulators
  Qout_day <- 0; Lout_day <- 0
  Uptake_day <- 0; Return_day <- 0

  # 4) Step loop (each hydrology step is Dt days; LPWEM uses ~20 steps/day in the paper)
  Np_per_day <- as.integer(safe_num(pars_lpwem$Np_per_day, 20L))

  step_diag <- vector("list", length(steps))

  for (k in seq_along(steps)) {
    st <- steps[[k]]
    Dt <- st$Dt

    # Step-average volume (hm3) from DMSTA staging
    Vavg <- 0.5 * (st$Vo + st$V)

    # Effective inflow that actually enters the wetland (exclude bypass)
    Q0 <- max(0, st$Qi_eff - safe_num(st$Bypass, 0))

    # Outflow from wetland for this step (hm3/day)
    QN <- max(0, st$Qout)

    # Microsteps within this hydro substep
    n_micro <- max(1L, as.integer(round(Np_per_day * Dt)))
    dt_micro <- Dt / n_micro

    # Cin fixed over the day (LPWEM uses interpolated daily Cin time series)
    pars_lpwem$Cin_ppb <- inputs$Ci

    for (m in seq_len(n_micro)) {
      res <- lpwem_tp_euler_step(
        TP = TP,
        Vavg_hm3 = Vavg,
        Q0_hm3d = Q0,
        QN_hm3d = QN,
        tis = tis,
        pars = pars_lpwem,
        dt_day = dt_micro
      )
      TP <- res$TP_new

      # Accumulate vertical flux totals (kg) over the day
      Uptake_day <- Uptake_day + sum(res$diag$Uptake_kgd) * dt_micro
      Return_day <- Return_day + sum(res$diag$Return_kgd) * dt_micro
    }

    # Outflow load for this hydrology step based on last tank TP (flow-weighted)
    TP_out_step <- TP[tis$N]
    Qout_day <- Qout_day + QN * Dt
    Lout_day <- Lout_day + (QN * TP_out_step) * Dt

    step_diag[[k]] <- list(
      step = k, TP_out = TP_out_step, TP_vec = TP,
      Q0 = Q0, QN = QN, Vavg = Vavg, Dt = Dt
    )
  }

  Qin <- st$Qi_eff
  TP_in_day <- inputs$Ci

  TP_out_day <- if (Qout_day > 0) Lout_day / Qout_day else NA_real_

  list(
    results = list(
      Qin = Qin,
      TP_in = TP_in_day,
      V_end = hyd$V_end,
      TP_out = TP_out_day,
      TP_profile_end = TP,
      Qout = Qout_day,
      L_out = Lout_day,
      Uptake = Uptake_day,
      Return = Return_day
    ),
    budgets = list(
      water = list(
        RainVol = hyd$RainVol,
        EtVol = hyd$EtVol,
        NetAtmo = hyd$NetAtmo,
        WB_in = hyd$WB_in,
        WB_out = hyd$WB_out,
        WB_err = hyd$WB_err,
        WB_rel = hyd$WB_rel
      ),
      mass = list(
        L_out = Lout_day,
        Uptake = Uptake_day,
        Return = Return_day
      )
    ),
    meta = list(
      Date = inputs$Date,
      Nsteps_hydro = Nsteps_hydro,
      Np_per_day = Np_per_day,
      step_diag = step_diag,
      params_lpwem_used = pars_lpwem
    )
  )
}


#' Run a multi-day LPWEM simulation over a time series
#'
#' Wraps [lpwem_day()] over a daily time series, carrying forward
#' end-of-day volume and the tank TP profile as initial conditions for the next
#' day.
#'
#' @param series Data frame containing daily drivers. Must include columns:
#' `Date`, `Qi`, `Ci`, `Rain`, `Et`, `Zcontrol`.
#' @param params_hydro List of DMSTA hydrology parameters.
#' @param pars_lpwem List of LPWEM parameters (see [lpwem_day()]).
#' @param V_init Optional starting volume (\eqn{hm^3}). If `NULL`, defaults
#'   to `params_hydro$A_cell * (params_hydro$Zinit/100)`.
#' @param TP_init Optional initial TP tank profile (ppb). If `NULL`, defaults
#'   to 20 in all tanks.
#' @param Nsteps_hydro Integer number of hydrology substeps per day passed to
#'   [lpwem_day()].
#'
#' @return A list with:
#' \itemize{
#'   \item `results`: data frame of daily outputs including end-of-day volume,
#'     outflow TP, outflow, load, uptake, and return totals
#'   \item `meta`: run metadata including final state and parameters used
#' }
#'
#'
#' @references
#' Juston, J. M., & Kadlec, R. H. (2019).
#' Data-driven modeling of phosphorus (P) dynamics in low-P stormwater wetlands.
#' *Environmental Modelling & Software*, 118, 226–240.
#' https://doi.org/10.1016/j.envsoft.2019.05.002
#'
#' @examples
#' \dontrun{
#' # series is an input data.frame similar to
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
#' params <- list(
#'   A_cell = 9.29,          # cell area; km2
#'   Zmin   = 2,             # minimum depth; cm
#'   Vmin   = 0,             # minimum volume; hm3
#'   Q_a = 2.852,            # discharge coef
#'   Q_b = 4.0,              # discharge exponent
#'   Zweir = 0,              # depth offset for outflow computation; cm
#'   Q_zmin = 35.97,         # minimum depth of discharge; cm
#'   Qomax = 0.0,            # maximum discharge  hm3/day
#'   Qimax = 0,              # maximum inflow  (hm3/day)
#'   Width = 2.00,           # cell widthl km
#'   Bypass_elev = 0,        # mean depth at which bypass begins, cm
#'   Seepout_Rate = 0.00789, # outflow seepage rate per unit head; cm/d/cm
#'   Seepout_Elev = 0.0,     # elevation controlling outflow seepage rate; cm
#'   Seepin_Rate  = 0.0,     # seepage inflow rate; ; cm/d/cm
#'   Seepin_Elev  = 0.0,     # elevation controlling inflow seepage rate; cm
#'   ShutdownET = TRUE,
#'   force_Q_out = FALSE,
#'   wrap_interp = TRUE,
#'   Zinit = 40,             # initial water column depth; cm
#'   Qin_Frac = 0.15,        # fraction of basin flows going into this cell
#'   Zrelease = 0,           #  minimum depth for releases; cm
#'   RecycleQ = 0
#' )
#'
#' pars_lpwem <- list(
#'   N = 6L,
#'   W_m = params$Width*1000,
#'   k0 = 0.0052,
#'   alpha = 0.0091,
#'   r1_mgm2d = 7.4,
#'   rN_mgm2d = 1.6,
#'   Np_per_day = 20L
#' )
#'
#' res <- lpwem_series(series, params_hydro = params, pars_lpwem = pars_lpwem)


#'
#' }
#'
#' @export
lpwem_series <- function(
    series,             # data.frame with Date, Qi, Ci, Rain, Et, Zcontrol
    params_hydro,       # DMSTA hydro params
    pars_lpwem,         # LPWEM params list
    V_init = NULL,      # initial V (hm3), default Zinit*A
    TP_init = NULL,     # initial TP profile (ppb), default 20 in all tanks
    Nsteps_hydro = 4L
) {
  req <- c("Date","Qi","Ci","Rain","Et","Zcontrol")
  miss <- setdiff(req, names(series))
  if (length(miss) > 0) stop("series missing columns: ", paste(miss, collapse=", "))
  n <- nrow(series)
  if (n < 1) stop("series has zero rows.")

  if (is.null(V_init)) {
    Z0_m <- (params_hydro$Zinit / 100)
    V_init <- params_hydro$A_cell * Z0_m
  }

  tis <- lpwem_build_tis(params_hydro$A_cell, pars_lpwem$N)
  if (is.null(TP_init)) TP_state <- rep(20, tis$N) else TP_state <- TP_init

  V <- V_init

  out <- data.frame(
    Date = as.Date(series$Date),
    Qin = NA_real_,
    TP_in = NA_real_,
    V_end = NA_real_,
    TP_out = NA_real_,
    Qout = NA_real_,
    L_out = NA_real_,
    Uptake = NA_real_,
    Return = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n)) {
    in_i <- as.list(series[i, ])

    # Neighbor Zcontrol handling (your helper)
    nz <- DMSTAr:::neighbors_zcontrol(i, series$Zcontrol)
    in_i$Zcontrol <- nz$today
    in_i$Zcontrol_prev <- nz$prev_day
    in_i$Zcontrol_next <- nz$nxt
    in_i$has_depth_constraint <- is.finite(in_i$Zcontrol) && in_i$Zcontrol != 0

    # Ensure optional fields exist for hydrology
    if (is.null(in_i$RecycleQ)) in_i$RecycleQ <- 0
    if (is.null(in_i$Qr0)) in_i$Qr0 <- 0
    if (is.null(in_i$Qr1)) in_i$Qr1 <- 0
    if (is.null(in_i$Qr2)) in_i$Qr2 <- 0

    res <- lpwem_day(
      V = V,
      TP_init = TP_state,
      inputs = in_i,
      params_hydro = params_hydro,
      pars_lpwem = pars_lpwem,
      Nsteps_hydro = Nsteps_hydro
    )

    V <- res$results$V_end
    TP_state <- res$results$TP_profile_end

    out$Qin[i]    <- res$results$Qin
    out$TP_in[i]     <- res$results$TP_in
    out$V_end[i]  <- V
    out$Qout[i]   <- res$results$Qout
    out$L_out[i]  <- res$results$L_out
    out$TP_out[i] <- res$results$TP_out
    out$Uptake[i] <- res$results$Uptake
    out$Return[i] <- res$results$Return
  }

  list(
    results = out,
    meta = list(
      V_init = V_init,
      V_end = V,
      TP_profile_end = TP_state,
      params_hydro = params_hydro,
      pars_lpwem = pars_lpwem,
      Nsteps_hydro = Nsteps_hydro
    )
  )
}
