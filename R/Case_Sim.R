#' DMSTA network/case helpers (internal)
#'
#' Internal helpers used for DMSTA multi-cell / network (case) simulations,
#' including per-cell initialization, validation, per-cell kinetics preparation,
#' daily-mean volume calculations, and output column normalization.
#'
#' These functions are used by [dmsta_flowP_case()] and related
#' network orchestration logic.
#'
#' @name internal_dmsta_case
#' @keywords internal
NULL

#' Daily-average cell volume from RK substeps (internal)
#'
#' Computes the DMSTA-style daily average volume as the time integral of the
#' substep mean volume:
#' \deqn{V_{day} = \sum_k \left(\frac{V_{o,k} + V_{k}}{2}\right)\Delta t_k}
#' where each substep record provides `Vo`, `V`, and `Dt`.
#'
#' @param step_list A list of substep records. Each element should be a list
#'   containing numeric fields `Vo` (start volume), `V` (end volume),
#'   and `Dt` (substep duration in days).
#'
#' @return Numeric scalar. Daily-average volume integral (same volume units as
#'   `Vo` and `V`).
#'
#' @rdname internal_dmsta_case
#' @keywords internal
dmsta_daily_avg_volume <- function(step_list) {
  # VBA: V_CellDay = sum(Vavg * Dt) over substeps, where Vavg = (Vo + V)/2
  sum(vapply(step_list, function(s) 0.5 * (s$Vo + s$V) * s$Dt, 0.0))
}

#' Prepare per-cell kinetics and constants for a network simulation (internal)
#'
#' Populates each cell definition with:
#' \itemize{
#'   \item `ppar`: per-cell phosphorus kinetic parameters (typically from
#'     `build_P_kin_slots(mods = c("STA","PSTA","RES"), ...)`),
#'   \item `constants`: a per-cell constants list used by the phosphorus
#'     mass-balance derivative (e.g., concentration caps, deposition terms,
#'     seepage bookkeeping fractions).
#' }
#'
#' The function returns the updated `cells` object.
#'
#' @param cells A list of cell definitions. Each element must contain
#'   `$params` with the raw parameters needed by the kinetics builder
#'   (and optional constants inputs).
#'
#' @details
#' This function uses the phosphorus model registry exposed as
#' `P_MODEL_BUILDERS` and assumes the `STA`, `PSTA`, and `RES`
#' builders are available in the registry.
#'
#' @return The updated `cells` list with `$ppar` and `$constants`
#'   added to each element.
#'
#' @seealso [dmsta_validate_cells()] for strict validation of the
#'   resulting cell definitions.
#'
#' @rdname internal_dmsta_case
#' @keywords internal
dmsta_prepare_cells_modules <- function(cells) {
  for (ic in seq_along(cells)) {
    p <- cells[[ic]]$params

    # build ppar (3 kinetic modules)
    cells[[ic]]$ppar <- build_P_kin_slots(
      mods     = c("STA", "PSTA", "RES"),
      registry = NULL,
      pparams  = p,
      Dpy      = 365.25,
      DutyCycle = p$DutyCycle
    )
    validate_P_paramsK(cells[[ic]]$ppar)

    # build constants (what your dmsta_deriv_mass expects)
    cells[[ic]]$constants <- list(
      Cmax = if (is.null(p$Cmax)) 2000 else p$Cmax,
      C_rain = if (is.null(p$C_rain)) 0 else p$C_rain,
      DryDepo = if (is.null(p$DryDepo)) 0 else (p$DryDepo / 365.25),  # mg/m2-day
      seepin_conc = if (!is.null(p$seepin_conc)) p$seepin_conc else 0,
      seepout_conc_max = if (!is.null(p$seepage_c)) p$seepage_c else 0,
      fseep_recycle = if (!is.null(p$fseep_recycle)) p$fseep_recycle else 0,
      fseep_out = if (!is.null(p$fseep_out)) p$fseep_out else 0
    )
  }
  cells
}

#' Validate and normalize network cell definitions (internal)
#'
#' Checks that `cells` is a non-empty list of cell definitions and that each
#' cell contains the minimum required fields for a coupled hydrology + phosphorus
#' simulation.
#'
#' The validator enforces:
#' \itemize{
#'   \item required structural fields: `params`, `ttankS`, `ppar`, `constants`
#'   \item required `ppar` fields used by the phosphorus derivative
#'   \item required `constants` fields used by the phosphorus derivative
#'   \item DMSTA convention for recycle indices (0 means self)
#'   \item at most one splitter cell labeled `"SPLITTER"` (case-insensitive)
#' }
#'
#' The function also fills defaults for some routing fields if missing.
#'
#' @param cells A list of cell definitions (as created by [dmsta_make_cell()]
#'   or assembled manually).
#'
#' @return The normalized `cells` list (with defaults filled in). Invisibly
#'   returns `cells` but typically used as `cells <- dmsta_validate_cells(cells)`.
#'
#' @export
dmsta_validate_cells <- function(cells) {
  if (!is.list(cells) || length(cells) < 1) stop("'cells' must be a non-empty list.")
  n <- length(cells)

  for (i in seq_len(n)) {
    cd <- cells[[i]]

    # check if a node
    # cd$IsaNode <- dmsta_is_node(cd$A_cell, NULL)
    # cells[[i]]$params <- cd

    if (is.null(cd$params)) stop("cells[[", i, "]]$params is missing.")
    if (is.null(cd$ttankS)) stop("cells[[", i, "]]$ttankS is missing.")

    # REQUIRE ppar and constants (package-safe design)
    if (is.null(cd$ppar)) stop("cells[[", i, "]]$ppar is missing. Build per-cell kinetics with build_P_kin_slots().")
    if (is.null(cd$constants)) stop("cells[[", i, "]]$constants is missing. Build per-cell constants list (Cmax, C_rain, DryDepo, etc.).")

    # minimal sanity checks for ppar fields used in DerivMass
    need_ppar <- c("K1","K2","K3","Z_1","Z_2","Z_3","Chalf","K2Coef","PModel","SeasonalFactor","CZero")
    miss_ppar <- setdiff(need_ppar, names(cd$ppar))
    if (length(miss_ppar) > 0) {
      stop("cells[[", i, "]]$ppar missing fields: ", paste(miss_ppar, collapse = ", "))
    }

    # minimal sanity checks for constants fields used in DerivMass
    need_const <- c("Cmax","C_rain","DryDepo","seepin_conc","seepout_conc_max","fseep_recycle","fseep_out")
    miss_const <- setdiff(need_const, names(cd$constants))
    if (length(miss_const) > 0) {
      stop("cells[[", i, "]]$constants missing fields: ", paste(miss_const, collapse = ", "))
    }

    # defaults for routing fields (these are not part of kinetics)
    if (is.null(cd$DownCell)) cd$DownCell <- 0L
    if (is.null(cd$Qin_Frac)) cd$Qin_Frac <- 0

    # VBA: recycleto=0 means self
    if (is.null(cd$RecycleIndex) || is.na(cd$RecycleIndex) || cd$RecycleIndex == 0L) {
      cd$RecycleIndex <- i
    }

    cells[[i]] <- cd
  }

  # single splitter supported like DMSTA VBA (one splitter per case)
  spl <- which(vapply(cells, function(x) identical(toupper(x$label), "SPLITTER"), logical(1)))
  if (length(spl) > 1) stop("Only one SPLITTER cell is supported (matches DMSTA).")
  if (length(spl) == 1 && is.null(cells[[spl]]$SplitterFrac)) {
    stop("Splitter cell must provide SplitterFrac.")
  }

  cells
}

#' Create a DMSTA network cell definition
#'
#' Constructs a cell definition used by [dmsta_flowP_case()] network/case
#' simulations. The function attaches:
#' \itemize{
#'   \item per-cell tank specification `ttankS`,
#'   \item routing fields (`DownCell`, `Qin_Frac`, `RecycleIndex`),
#'   \item optional splitter specification (`SplitterFrac`),
#'   \item per-cell kinetics parameters (`ppar`) using built-in modules
#'     (`STA`, `PSTA`, `RES`),
#'   \item per-cell constants list used by the phosphorus derivative.
#' }
#'
#' @param label Character scalar. Cell label (use `"SPLITTER"` for the
#'   splitter cell if your case uses one).
#' @param params Named list of per-cell parameters (must include at least
#'   `A_cell`, `DutyCycle`, `Zinit`, `C_init_ppb`, `Y_init_mgm2`,
#'   plus fields required by the selected kinetics builders).
#' @param ttankS Numeric scalar > 0. Effective number of tanks (may be fractional).
#' @param DownCell Integer index of the downstream cell (0 means terminal/out of system).
#' @param Qin_Frac Numeric scalar. Fraction of watershed inflow assigned to this cell
#'   in network routing (case runner may override internally for “already allocated” inflows).
#' @param RecycleIndex Integer index of the cell receiving recycled seepage flow/mass.
#'   `0` or `NA` indicates default handling (often self).
#' @param SplitterFrac Optional. Splitter routing fractions (only used if
#'   `label` is `"SPLITTER"`). May be a numeric vector (length = number of cells)
#'   or a named vector/list mapping downstream indices to fractions.
#'
#' @return A named list representing a cell definition, suitable for inclusion
#'   in the `cells` argument of [dmsta_flowP_case()].
#'
#' @seealso [dmsta_flowP_case()] for running a multi-cell simulation.
#'
#' @export
dmsta_make_cell <- function(label, params, ttankS, DownCell = 0L, Qin_Frac = 0,
                            RecycleIndex = NULL, SplitterFrac = NULL) {
  if (is.null(RecycleIndex)) RecycleIndex <- NA_integer_

  if (is.null(params$Zrelease) || !is.finite(params$Zrelease)) {
    params$Zrelease <- 0
  }

  # build per-cell ppar (3 kinetic modules)
  ppar <- build_P_kin_slots(
    mods     = c("STA","PSTA","RES"),
    registry = NULL,
    pparams  = params,
    Dpy      = 365.25,
    DutyCycle = params$DutyCycle
  )
  validate_P_paramsK(ppar)

  # build per-cell constants (DerivMass needs mg/m2-day for DryDepo like VBA does /Dpy)
  constants <- list(
    Cmax = if (is.null(params$Cmax)) 2000 else params$Cmax,
    C_rain = if (is.null(params$C_rain)) 0 else params$C_rain,
    DryDepo = if (is.null(params$DryDepo)) 0 else (params$DryDepo / 365.25),
    seepin_conc = if (!is.null(params$seepin_conc)) params$seepin_conc else 0,
    seepout_conc_max = if (!is.null(params$seepage_c)) params$seepage_c else 0,
    fseep_recycle = if (!is.null(params$fseep_recycle)) params$fseep_recycle else 0,
    fseep_out = if (!is.null(params$fseep_out)) params$fseep_out else 0
  )

  list(
    label = label,
    params = params,
    ttankS = ttankS,
    DownCell = as.integer(DownCell),
    Qin_Frac = Qin_Frac,
    RecycleIndex = if (is.null(RecycleIndex) || is.na(RecycleIndex) || RecycleIndex == 0L) NA_integer_ else as.integer(RecycleIndex),
    SplitterFrac = SplitterFrac,
    ppar = ppar,
    constants = constants
  )
}

#' Initialize per-cell state for a network/case run
#'
#' Builds tank geometry and initializes hydrologic and phosphorus state for each
#' cell using the cell parameters:
#' \itemize{
#'   \item `V[ic]` is initialized as `A_cell * (Zinit/100)` (depth cm -> m),
#'   \item tank geometry is built using `dmsta_build_tanks(A_cell, ttankS)`,
#'   \item phosphorus states `M` and `S` are initialized with
#'     `dmsta_p_init_state()` using `C_init_ppb` and `Y_init_mgm2`.
#' }
#'
#' @param cells A validated list of cell definitions. Each cell must contain
#'   `$params` (with `A_cell`, `Zinit`, `C_init_ppb`, `Y_init_mgm2`)
#'   and `$ttankS`.
#'
#' @return A list with elements:
#' \describe{
#'   \item{V}{Numeric vector of initial volumes, length = number of cells.}
#'   \item{tanks}{List of per-cell tank geometry objects.}
#'   \item{Pstate}{List of per-cell phosphorus state objects (each with `M` and `S` vectors).}
#' }
#'
#' @export
dmsta_init_case_state <- function(cells) {
  ncell <- length(cells)
  tanks <- vector("list", ncell)
  V     <- numeric(ncell)
  Pstate<- vector("list", ncell)

  for (ic in seq_len(ncell)) {
    ## if its a node, should still work
    p <- cells[[ic]]$params
    Z0 <- p$Zinit / 100
    V[ic] <- p$A_cell * Z0

    tanks[[ic]] <- dmsta_build_tanks(p$A_cell, cells[[ic]]$ttankS)

    Pstate[[ic]] <- dmsta_p_init_state(
      tanks[[ic]],
      Z_init_m    = Z0,
      C_init_ppb  = p$C_init_ppb,
      Y_init_mgm2 = p$Y_init_mgm2
    )
  }

  list(V = V, tanks = tanks, Pstate = Pstate)
}

#' Normalize case/cell output columns to a common schema (internal)
#'
#' Utility to standardize output column names and structure between per-cell and
#' case-level result data frames. The function:
#' \itemize{
#'   \item maps known aliases (e.g., `RainVol_total` -> `RainVol`),
#'   \item builds a consistent set of inflow/outflow component columns,
#'   \item optionally keeps seep recycle (Q17/L17/C17),
#'   \item optionally carries through P-budget fields (e.g., columns prefixed with `P_`),
#'   \item optionally preserves specified list-columns (e.g., `mass_budget`).
#' }
#'
#' @param df A data.frame containing DMSTA-style result columns.
#' @param keep_Q17 Logical; if `TRUE`, keep seep recycle fields (Q17/L17/C17) if present.
#' @param keep_P Logical; if `TRUE`, carry through phosphorus budget fields.
#' @param P_prefix Regular expression used to match P-budget fields (default `"^P_"`).
#' @param keep_list_cols Character vector of list-column names to carry through if present.
#' @param keep_extra Logical; if `TRUE`, carry through any remaining fields not
#'   in the standardized output.
#'
#' @return A data.frame with standardized columns, plus optional carry-through fields.
#'
#' @rdname internal_dmsta_case
#' @keywords internal
dmsta_case_components <- function(
    df,
    keep_Q17 = TRUE,
    keep_P = TRUE,                 # keep P_* fields unchanged
    P_prefix = "^P_",              # regex to match P budget fields
    keep_list_cols = c("mass_budget"),# list-columns to carry through if present
    keep_extra = TRUE              # also keep any other unformatted fields (recommended)
) {
  ## internal function to help reformat output columns

  # base-R safe; works even if some columns were dropped
  has <- function(nm) nm %in% names(df)
  # get a single column or NA vector
  get <- function(nm) {
    if (has(nm)) df[[nm]] else rep(NA_real_, nrow(df))
  }

  # get first available among aliases (e.g., RainVol vs RainVol_total)
  get_first <- function(nms) {
    for (nm in nms) if (has(nm)) return(df[[nm]])
    rep(NA_real_, nrow(df))
  }
  # safe concentration helper

  ## map fields that differ between cell vs case
  # State / water levels
  V_end <- get_first(c("V_end", "V_total_end"))

  # Prefer patched (true) depths first; fall back to legacy totals if present
  Z_end <- get_first(c("Z_end", "Z_total_end"))
  Z_avg <- get_first(c("Z_avg", "Z_total_avg"))

  # DMSTA-style daily avg volume (midpoint-integrated); keep legacy alias too
  V_cell_day <- get_first(c("V_cell_day", "V_cell_day_total"))

  ## Z_* should always be present ... commenting out for now.
  ## Optional fallbacks if Z_* not present but volumes + A_cell are
  # A_cell <- get_first(c("A_cell", "A_cell_km2", "A_cell_total"))

  ## If Z_end missing but V_end and A_cell exist, compute Z_end = V_end / A_cell
  # if (all(!is.finite(Z_end)) && any(is.finite(V_end)) && any(is.finite(A_cell)) && all(A_cell > 0, na.rm = TRUE)) {
  #   Z_end <- V_end / A_cell
  # }

  ## If Z_avg missing but V_cell_day and A_cell exist, compute Z_avg = V_cell_day / A_cell
  # if (all(!is.finite(Z_avg)) && any(is.finite(V_cell_day)) && any(is.finite(A_cell)) && all(A_cell > 0, na.rm = TRUE)) {

  # Atmospheric totals
  RainVol <- get_first(c("RainVol", "RainVol_total"))
  EtVol   <- get_first(c("EtVol",   "EtVol_total"))
  NetAtmo <- get_first(c("NetAtmo", "NetAtmo_total"))

  # Water budget totals
  # WB_in  <- get_first(c("WB_in",  "WB_in_total"))
  # WB_out <- get_first(c("WB_out", "WB_out_total"))
  # WB_err <- get_first(c("WB_err", "WB_err_total"))
  # WB_rel <- get_first(c("WB_rel", "WB_rel_total"))

  # Inflow totals (cell has Qin/Lin/Cin; case has Qin/Lin/Cin too)
  Q_in_total <- get("Qin")
  L_in_total <- get("Lin")
  C_in_total <- get("Cin")

  # Basin/upstream decomposition may be cell-only
  Q_in_basin   <- get("Qin_basin")
  L_in_basin   <- get("Lin_basin")
  C_in_basin   <- get("Cin_basin")
  Q_in_upstream <- get("Qin_up")
  L_in_upstream <- get("Lin_up")
  C_in_upstream <- get("Cin_up")

  # Treated inflow term (Stream 7)
  Q_in_treated <- get("Q7")
  L_in_treated <- get("L7")
  C_in_treated <- get("C7")

  # Outflow terms (DMSTA-style)
  Q_out_bypass <- get("Q3")
  L_out_bypass <- get("L3")
  C_out_bypass <- get("C3")

  Q_out_treated <- get("Q13")
  L_out_treated <- get("L13")
  C_out_treated <- get("C13")

  Q_out_seep_discharge <- get("Q14")
  L_out_seep_discharge <- get("L14")
  C_out_seep_discharge <- get("C14")

  Q_out_release1 <- get("Q25")
  L_out_release1 <- get("L25")
  C_out_release1 <- get("C25")

  Q_out_release2 <- get("Q26")
  L_out_release2 <- get("L26")
  C_out_release2 <- get("C26")

  ## build formatted output
  out <- data.frame(
    Date = as.Date(df$Date),

    # state/level
    V_end      = V_end,
    Z_end      = Z_end,
    Z_avg      = Z_avg,
    Z_end_cm   = m_to_cm(Z_end),
    Z_avg_cm   = m_to_cm(Z_avg),
    V_cell_day = V_cell_day,

    # atmospheric
    RainVol = RainVol,
    EtVol   = EtVol,
    NetAtmo = NetAtmo,

    # water budget
    # WB_in  = WB_in,
    # WB_out = WB_out,
    # WB_err = WB_err,
    # WB_rel = WB_rel,

    # inflow totals and decomposition
    Q_in_total = Q_in_total,
    L_in_total = L_in_total,
    C_in_total = C_in_total,

    Q_in_basin = Q_in_basin,
    L_in_basin = L_in_basin,
    C_in_basin = C_in_basin,

    Q_in_upstream = Q_in_upstream,
    L_in_upstream = L_in_upstream,
    C_in_upstream = C_in_upstream,

    # treated inflow
    Q_in_treated = Q_in_treated,
    L_in_treated = L_in_treated,
    C_in_treated = C_in_treated,

    # outflows
    Q_out_bypass = Q_out_bypass,
    L_out_bypass = L_out_bypass,
    C_out_bypass = C_out_bypass,

    Q_out_treated = Q_out_treated,
    L_out_treated = L_out_treated,
    C_out_treated = C_out_treated,

    Q_out_seep_discharge = Q_out_seep_discharge,
    L_out_seep_discharge = L_out_seep_discharge,
    C_out_seep_discharge = C_out_seep_discharge,

    Q_out_release1 = Q_out_release1,
    L_out_release1 = L_out_release1,
    C_out_release1 = C_out_release1,

    Q_out_release2 = Q_out_release2,
    L_out_release2 = L_out_release2,
    C_out_release2 = C_out_release2,

    stringsAsFactors = FALSE
  )

  ## optional seep recycle term (Q17/L17/C17)
  if (keep_Q17 && has("Q17")) {
    out$Q_seep_recycle <- get("Q17")
    out$L_seep_recycle <- get("L17")
    out$C_seep_recycle <- get("C17")
  }

  ## rollups
  out$Q_out_total <- rowSums(
    out[, intersect(names(out), c("Q_out_bypass","Q_out_treated","Q_out_seep_discharge","Q_out_release1","Q_out_release2")),
        drop = FALSE],
    na.rm = TRUE
  )
  out$L_out_total <- rowSums(
    out[, intersect(names(out), c("L_out_bypass","L_out_treated","L_out_seep_discharge","L_out_release1","L_out_release2")),
        drop = FALSE],
    na.rm = TRUE
  )
  out$C_out_total <- fw(out$L_out_total, out$Q_out_total)

  ## NEW: carry-through P budget fields unchanged
  if (isTRUE(keep_P)) {
    p_cols <- grep(P_prefix, names(df), value = TRUE)
    p_cols <- setdiff(p_cols, names(out)) # avoid duplicates
    if (length(p_cols) > 0) out[p_cols] <- df[p_cols]

    # carry list-columns like mass_budget
    for (lc in keep_list_cols) {
      if (has(lc) && !(lc %in% names(out))) out[[lc]] <- df[[lc]]
    }
  }

  ## OPTIONAL: carry-through any other fields unchanged
  if (isTRUE(keep_extra)) {
    extra <- setdiff(names(df), names(out))
    if (length(extra) > 0) out[extra] <- df[extra]
  }

  out
}

#' Run a networked DMSTA hydrology–phosphorus simulation
#'
#' Simulates coupled hydrology and phosphorus dynamics for a network of
#' interconnected DMSTA cells over a daily time series. Each cell is
#' simulated sequentially within each day, with treated outflows routed
#' downstream according to network topology, splitter rules, and recycle
#' indices.
#'
#' This function manages network‑level state, routing, lagged recycle
#' bookkeeping, and optional convergence iteration, while delegating
#' per‑cell daily physics to `dmsta_flowP_day()`.
#'
#' @param series Data frame of daily watershed inputs. Must include
#'   `Date`, `Qi`, `Ci`, `Rain`, `Et`, and
#'   `Zcontrol`.
#' @param cells List of cell definitions created by
#'   `dmsta_make_cell()` and validated with
#'   `dmsta_validate_cells()`.
#' @param Nsteps Integer. Number of hydrology sub‑steps per day.
#' @param N_plant Integer. Window length (days) for rolling mean depth
#'   used in reservoir penalty blending.
#' @param Qmethod Character string specifying the hydrology integrator
#'   (`"RK4"`, `"Euler"`, `"RKF45"`, or `"custom"`).
#' @param Pmethod Character string specifying the phosphorus integrator
#'   (`"RK4"` or `"Euler"`).
#' @param integrator_fun Optional custom hydrology integrator function.
#' @param interp_option Control‑depth interpolation option (DMSTA semantics).
#' @param max_iter Integer. Maximum number of network convergence iterations.
#' @param conv_tol Numeric. Relative convergence tolerance on external
#'   phosphorus loads.
#' @param return_cell_series Logical. If `TRUE`, return per‑cell
#'   daily output series.
#' @param keep_Q17 Logical. If `TRUE`, retain seep‑recycle bookkeeping
#'   streams (Q17/L17/C17).
#' @param ... Additional arguments passed to daily hydrology integration.
#'
#' @return An object of class `"dmsta_network_result"` containing:
#' \describe{
#'   \item{results}{Case‑level and optional per‑cell daily output series.}
#'   \item{budgets}{Water and phosphorus mass budgets at case and cell level.}
#'   \item{meta}{Convergence diagnostics and model configuration metadata.}
#' }
#'
#' @details
#' Network routing follows DMSTA conventions. Treated outflows are routed
#' either to a downstream cell or out of system, while bypass, release,
#' and seepage discharge streams leave the system immediately. Lagged
#' seepage recycle is tracked explicitly as an internal transit reservoir.
#'
#' For strict DMSTA parity, use `Qmethod = "RK4"` and
#' `Pmethod = "RK4"` with no operational overrides.
#'
#'
#' @examples
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
#' # --- 1) Base hydrology params (shared structure) ---
#' hydro_base <- list(
#'   A_cell = 2.19, # km2
#'   # depths in cm (engine converts /100 to meters internally)
#'   Zmin   = 2,              # cm
#'   Zinit  = 40,             # cm
#'   Zweir  = 0,              # cm
#'   Q_zmin = 38,             # cm
#'   Zrelease = 0,            # cm
#'   Bypass_elev = 121.92,    # cm (≈ 1.2192 m)
#'   # hydraulics
#'   Q_a = 1.0,
#'   Q_b = 4.0,
#'   Width = 1.55,            # km
#'   Qomax = 0.0,             # hm3/day; 0 disables max cap in this implementation
#'   Qimax = 0.0,             # hm3/day; 0 disables inflow cap
#'   # seepage (rates in m/day per m head; elevations in cm)
#'   Seepout_Rate = 0.0,
#'   Seepout_Elev = 0.0,      # cm
#'   Seepin_Rate  = 0.0,
#'   Seepin_Elev  = 0.0,      # cm
#'   ShutdownET = TRUE,
#'   force_Q_out = FALSE,
#'   DutyCycle = 0.95,
#'   Cmax = 2000
#' )
#'
#' # 2) Base P params (shared structure)
#' P_base <- list(
#'   # STA module
#'   C1000 = 22,
#'   Cstar = 3,
#'   Ks_per_yr = 16.8,
#'   Z1 = 40,
#'   Z2 = 100,
#'   Z3 = 200,
#'   Chalf = 300,
#'   K2Coef1 = 0,
#'   SeasonalFactor = 0,    # keep 0 for base parity
#'   # PSTA (NEWS transition)
#'   Ytrans = 0,
#'   Ysigma = 0,
#'   Czero = 0,
#'   C1000_2 = NULL,
#'   ks_2 = 0,
#'   zh_2 = 0,
#'   # RES depth penalty
#'   k_depth_penalty = 1,
#'   # atmos + seepage water quality
#'   C_rain = 10,           # ppb (ug/L)
#'   DryDepo = 20,          # mg/m2-yr
#'   seepage_c = 20,        # ppb cap for seep outflow
#'   seepin_conc = 0,       # ppb
#'   # initial P state
#'   C_init_ppb = 30,
#'   Y_init_mgm2 = 1000
#'   )
#'
#'   #  3) Cell-specific params
#'   params_cell1 <- modifyList(hydro_base, modifyList(P_base, list(
#'   Qin_Frac = 0.22,
#'   Seepout_Rate = 0.00789,
#'   Ks_per_yr = 16.8,
#'   Y_init_mgm2 = 3387.67297548954
#'   )))
#'
#'   params_cell2 <- modifyList(hydro_base, modifyList(P_base, list(
#'   Qin_Frac = 0,
#'   Seepout_Rate = 0.00155,
#'   Ks_per_yr = 52.5,
#'   Y_init_mgm2 = 768.480041186681
#'   )))
#'
#'   # 4) Build cells
#'   cells <- list(
#'   dmsta_make_cell(
#'   label = "CELL1",
#'   params = params_cell1,
#'   ttankS = 3.0,
#'   DownCell = 2L,
#'   Qin_Frac = params_cell1$Qin_Frac,
#'   RecycleIndex = 1L  # self; can omit if your validator maps NA/0 -> self
#'   ),
#'   dmsta_make_cell(
#'   label = "CELL2",
#'   params = params_cell2,
#'   ttankS = 3.0,
#'   DownCell = 0L,
#'   Qin_Frac = params_cell2$Qin_Frac,
#'   RecycleIndex = 2L  # self
#'   )
#'   )
#'
#'   cells <- dmsta_validate_cells(cells)
#'
#'   # Run case
#'   out <- dmsta_flowP_case(
#'    series = series,
#'    cells  = cells,
#'    Nsteps = 4L,
#'    max_iter = 1L,
#'    return_cell_series = TRUE,
#'    keep_Q17 = TRUE,
#'    keep_extra = FALSE
#'   )
#'   head(out$results$case)
#'   head(out$results$cells[[1]])
#'   head(out$results$cells[[2]])
#'
#' @export
dmsta_flowP_case <- function(
  series,
  cells,
  Nsteps = 4L,
  N_plant = 30L,
  Qmethod = c("RK4", "Euler", "RKF45", "custom"),
  Pmethod = c("RK4", "Euler"),
  integrator_fun = NULL,
  interp_option = 2L,
  max_iter = 1L,
  conv_tol = 0.01,
  return_cell_series = TRUE,
  keep_Q17 = TRUE,
  ...
) {

  Qmethod <- match.arg(Qmethod)
  Pmethod <- match.arg(Pmethod)

  # Input checks
  req <- c("Date","Qi","Ci","Rain","Et","Zcontrol")
  miss <- setdiff(req, names(series))
  if (length(miss) > 0) stop("series missing: ", paste(miss, collapse=", "))

  cells <- dmsta_validate_cells(cells)
  ncell <- length(cells)
  nday  <- nrow(series)

  A_cells <- vapply(cells, function(cd) cd$params$A_cell, numeric(1))
  A_total <- sum(A_cells[is.finite(A_cells) & A_cells > 0], na.rm = TRUE)

  # Identify splitter (0 if none)
  spl_idx <- which(vapply(cells, function(x) identical(toupper(x$label), "SPLITTER"), logical(1)))
  spl_idx <- if (length(spl_idx) == 0) 0L else spl_idx[1]
  if (length(spl_idx) > 1) stop("Only one SPLITTER cell is supported (matches DMSTA).")
  if (spl_idx > 0L && is.null(cells[[spl_idx]]$SplitterFrac)) {
    stop("Splitter cell must provide SplitterFrac.")
  }

  # DMSTA convention: recycleto=0 means recycle to self
  for (i in seq_len(ncell)) {
    ri <- cells[[i]]$RecycleIndex
    if (is.null(ri) || is.na(ri) || ri == 0L) cells[[i]]$RecycleIndex <- i
  }

  # Output templates
  make_cell_df <- function() {
    data.frame(
      Date = as.Date(series$Date),
      # volume and water level
      V_end = NA_real_,
      Z_end = NA_real_, Z_avg = NA_real_,
      V_cell_day = NA_real_,
      # total inflow to cell before bypass
      Qin = NA_real_,  Cin = NA_real_,  Lin = NA_real_,
      # components of Qin
      Qin_basin = NA_real_, Cin_basin = NA_real_, Lin_basin = NA_real_,
      Qin_up    = NA_real_, Cin_up    = NA_real_, Lin_up    = NA_real_,
      # Term 7: stream in, flow into cell (treated inflow = Qin - bypass; exclude recycle)
      Q7 = NA_real_,  C7 = NA_real_,  L7 = NA_real_,
      # Term 3: bypass
      Q3 = NA_real_,  C3 = NA_real_,  L3 = NA_real_,
      # Term 13: treated outflow (routed downstream; terminal only to out-of-system)
      Q13 = NA_real_, C13 = NA_real_, L13 = NA_real_,
      # Term 14: seppage discharge out of system
      Q14 = NA_real_, C14 = NA_real_, L14 = NA_real_,
      # Term 25/26: release 1 and 2 out of system
      Q25 = NA_real_, C25 = NA_real_, L25 = NA_real_,
      Q26 = NA_real_, C26 = NA_real_, L26 = NA_real_,
      # Term 17: seepage recycling bookkeeping (lagged +1 day)
      Q17 = NA_real_, C17 = NA_real_, L17 = NA_real_,
      # offline values
      Qin_frac_base = NA_real_,
      Qin_frac_used = NA_real_,
      offline = NA,
      offline_index = NA_integer_,
      off_diff = NA_integer_,
      off_mod = NA_integer_,
      stringsAsFactors = FALSE
    )
  }
  make_cell_wb_df <- function() {
    data.frame(
      Date = as.Date(series$Date),
      RainVol = NA_real_, EtVol = NA_real_, NetAtmo = NA_real_,
      WB_in  = NA_real_, WB_out = NA_real_,
      WB_err = NA_real_, WB_rel = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  make_cell_mb_df <- function() {
    data.frame(
      Date = as.Date(series$Date),
      dP = NA_real_,
      Pin_total = NA_real_, Pout_total = NA_real_, Perr_total = NA_real_, Prel_total = NA_real_,
      Pin_external = NA_real_, Pout_external = NA_real_, Perr_external = NA_real_, Prel_external = NA_real_,
      Q_in_tanks = NA_real_, L_in_tanks = NA_real_, C_in_tanks = NA_real_,
      L_rain = NA_real_, L_drydep = NA_real_, L_seepin = NA_real_,
      L_treated = NA_real_, L_rel1 = NA_real_, L_rel2 = NA_real_,
      L_bypass = NA_real_, L_seep_discharge = NA_real_,
      L_seep_recycle_out = NA_real_,
      L_uptake = NA_real_, L_recycle = NA_real_, L_sed = NA_real_, L_direct = NA_real_,
      L_burial = NA_real_, L_release = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  make_case_df <- function() {
    data.frame(
      Date = as.Date(series$Date),
      Qin = series$Qi,
      Cin = NA_real_,
      Lin = series$Qi * series$Ci,
      Q7  = NA_real_, C7  = NA_real_, L7  = NA_real_,
      Q3  = NA_real_, C3  = NA_real_, L3  = NA_real_,
      Q13 = NA_real_, C13 = NA_real_, L13 = NA_real_,
      Q14 = NA_real_, C14 = NA_real_, L14 = NA_real_,
      Q25 = NA_real_, C25 = NA_real_, L25 = NA_real_,
      Q26 = NA_real_, C26 = NA_real_, L26 = NA_real_,
      Q17 = NA_real_, C17 = NA_real_, L17 = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  make_case_wb_df <- function() {
    data.frame(
      Date = as.Date(series$Date),
      RainVol_total = NA_real_,
      EtVol_total   = NA_real_,
      NetAtmo_total = NA_real_,
      WB_in_total   = NA_real_,
      WB_out_total  = NA_real_,
      WB_err_total  = NA_real_,
      WB_rel_total  = NA_real_,
      V_total_end   = NA_real_,
      dV_total      = NA_real_,
      Z_total_end   = NA_real_,
      Z_total_avg   = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  make_case_mb_df <- function() {
    data.frame(
      Date = as.Date(series$Date),
      P_in_watershed = NA_real_,
      P_in_rain_total = NA_real_,
      P_in_drydep_total = NA_real_,
      P_in_seepin_total = NA_real_,
      P_in_external = NA_real_,
      P_out_bypass = NA_real_,
      P_out_terminal_treated = NA_real_,
      P_out_releases = NA_real_,
      P_out_seep_discharge = NA_real_,
      P_out_external = NA_real_,
      P_transit_start = NA_real_,
      P_transit_end   = NA_real_,
      dP_transit      = NA_real_,
      P_cells_start = NA_real_,
      P_cells_end   = NA_real_,
      dP_cells      = NA_real_,
      P_err_case = NA_real_,
      P_rel_case = NA_real_,
      P_mech_uptake_total = NA_real_,
      P_mech_recycle_total = NA_real_,
      P_mech_sed_total = NA_real_,
      P_mech_direct_total = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  # helpers: force scalar numeric (DMSTA semantics: missing => 0)
  .scalar0 <- function(x, default = 0) {
    if (is.null(x) || length(x) == 0L) return(default)
    x <- suppressWarnings(as.numeric(x[1]))
    if (!is.finite(x)) default else x
  }
  # concentration helper: prefer explicit conc; else compute from L/Q
  .scalarC <- function(Cx, Lx, Qx, default = 0) {
    Cx0 <- .scalar0(Cx, default = NA_real_)
    if (is.finite(Cx0)) return(Cx0)
    Lx0 <- .scalar0(Lx, default = 0)
    Qx0 <- .scalar0(Qx, default = 0)
    if (Qx0 > 0) Lx0 / Qx0 else default
  }

  .scalarNA <- function(x) {
    if (is.null(x) || length(x) == 0L) return(NA_real_)
    x <- suppressWarnings(as.numeric(x[1]))
    if (!is.finite(x)) NA_real_ else x
  }


  cell_out <- cell_wb <- cell_mb <- NULL
  if (return_cell_series) {
    cell_out <- vector("list", ncell)
    cell_wb  <- vector("list", ncell)
    cell_mb  <- vector("list", ncell)
    for (ic in seq_len(ncell)) {
      cell_out[[ic]] <- make_cell_df()
      cell_wb[[ic]]  <- make_cell_wb_df()
      cell_mb[[ic]]  <- make_cell_mb_df()
    }
  }

  case_out <- make_case_df()
  case_wb  <- make_case_wb_df()
  case_mb  <- make_case_mb_df()

  case_out$Qin <- series$Qi
  case_out$Lin <- series$Qi * series$Ci
  case_out$Cin <- ifelse(case_out$Qin > 0, case_out$Lin / case_out$Qin, 0)

  for (nm in c("Q7","L7","Q3","L3","Q13","L13","Q14","L14","Q25","L25","Q26","L26","Q17","L17")) {
    case_out[[nm]] <- 0.0
  }

  conv_hist <- data.frame(iter=integer(0), totalL=numeric(0), conv_test=numeric(0))
  converged <- FALSE
  iterations_used <- 0L
  prev_totalL <- NA_real_

  ## Series‑level HydroIndex semantics
  # Per-cell HydroIndex semantics (structural flags)
  # Structural presence (HydroIndex-style) should be per cell.
  has_Qr0_cell <- vapply(cells, function(cd) {
    nm <- cd$params$QR0_name
    !is.null(nm) && is.character(nm) && nzchar(nm)
  }, logical(1))

  has_Qr1_cell <- vapply(cells, function(cd) {
    nm <- cd$params$QR1_name
    !is.null(nm) && is.character(nm) && nzchar(nm)
  }, logical(1))

  has_Qr2_cell <- vapply(cells, function(cd) {
    nm <- cd$params$QR2_name
    !is.null(nm) && is.character(nm) && nzchar(nm)
  }, logical(1))

  has_Zcon_cell <- vapply(cells, function(cd) {
    nm <- cd$params$Zcon_name
    !is.null(nm) && is.character(nm) && nzchar(nm)
  }, logical(1))

  # has_Qr0_series <- any(is.finite(series$Qr0) & series$Qr0 != 0)
  # has_Qr1_series <- any(is.finite(series$Qr1) & series$Qr1 != 0)
  # has_Qr2_series <- any(is.finite(series$Qr2) & series$Qr2 != 0)
  # has_depth_constraint_series <- any(is.finite(series$Zcontrol) & series$Zcontrol != 0)

  ## DMSTA HydroIndex(1)
  has_depth_constraint_series <- any(is.finite(series$Zcontrol) & series$Zcontrol != 0)

  # Effective depth-constraint flag per cell:
  # - if the cell structurally has a depth-series mapped (Zcon_name), honor it
  # - OR if the global series provides a nonzero Zcontrol (RS-style), honor it
  has_depth_constraint_cell <- has_Zcon_cell | has_depth_constraint_series


  # Offline schedules per cell (precomputed; DMSTA 2C2B) ---
  offline_by_cell <- vector("list", ncell)

  for (ic in seq_len(ncell)) {
    cd <- cells[[ic]]
    p  <- cd$params

    dmsta_version <- if (!is.null(p$dmsta_version)) as.character(p$dmsta_version) else "2E"
    base_frac <- if (!is.null(cd$Qin_Frac) && is.finite(as.numeric(cd$Qin_Frac))) as.numeric(cd$Qin_Frac) else 0.0

    # defaults
    off_start <- if (!is.null(p$offline_start) && !is.na(p$offline_start)) as.Date(p$offline_start) else as.Date("1965-03-15")
    off_freq  <- if (!is.null(p$offline_freq)  && is.finite(as.numeric(p$offline_freq)))  as.integer(p$offline_freq) else 3L
    off_dur   <- if (!is.null(p$offline_dur)   && is.finite(as.numeric(p$offline_dur)))   as.integer(p$offline_dur) else 45L

    fracs6 <- c(p$frac_1, p$frac_2, p$frac_3, p$frac_4, p$frac_5, p$frac_6)

    if (!identical(dmsta_version, "2C2B")) {
      # non-2C2B: always base fraction, no offline
      offline_by_cell[[ic]] <- data.frame(
        Date = as.Date(series$Date),
        Qin_frac_base = base_frac,
        Qin_frac_used = base_frac,
        offline = FALSE,
        offline_index = NA_integer_,
        off_ini = as.Date(NA),
        off_diff = NA_integer_,
        off_mod = NA_integer_,
        frac_selected = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      tmp <- lapply(seq_len(nday), function(d) {
        .dmsta_offline_qin_frac(
          date = series$Date[d],
          base_frac = base_frac,
          offline_trigger = isTRUE(p$offline_trigger),
          offline_start = off_start,
          offline_freq  = off_freq,
          offline_dur   = off_dur,
          offline_fracs = fracs6
        )
      })

      offline_by_cell[[ic]] <- data.frame(
        Date = as.Date(series$Date),
        Qin_frac_base = vapply(tmp, `[[`, numeric(1), "Qin_frac_base"),
        Qin_frac_used = vapply(tmp, `[[`, numeric(1), "Qin_frac_used"),
        offline       = vapply(tmp, `[[`, logical(1), "offline"),
        offline_index = vapply(tmp, `[[`, integer(1), "offline_index"),
        off_ini       = as.Date(vapply(tmp, function(x) as.character(x$off_ini), character(1))),
        off_diff      = vapply(tmp, `[[`, integer(1), "off_diff"),
        off_mod       = vapply(tmp, `[[`, integer(1), "off_mod"),
        frac_selected = vapply(tmp, `[[`, numeric(1), "frac_selected"),
        stringsAsFactors = FALSE
      )
    }
  }

  # Iteration loop
  for (iter in seq_len(max_iter)) {

    # Initialize states (base)
    st     <- dmsta_init_case_state(cells)
    V      <- st$V
    tanks  <- st$tanks
    Pstate <- st$Pstate

    # Initialize V and Pstate ONCE per iteration (DMSTA Initialize)
    for (ic in seq_len(ncell)) {
      p <- cells[[ic]]$params

      p$force_Q_out <- has_Qr0_cell[ic]# has_Qr0_series

      # if (isTRUE(has_depth_constraint_series)) {
      # if (isTRUE(has_Zcon_cell[ic])){
      if (isTRUE(has_depth_constraint_cell[ic])) {
        zc1 <- series$Zcontrol[1]   # meters
        if (!is.finite(zc1)) zc1 <- 0
        Z0_m <- max(zc1, p$Zmin / 100)
      } else {
        Z0_m <- max(p$Zinit / 100, p$Zmin / 100)
      }

      V[ic] <- p$A_cell * Z0_m
      Pstate[[ic]] <- dmsta_p_init_state(
        tanks[[ic]],
        Z_init_m    = Z0_m,
        C_init_ppb  = p$C_init_ppb,
        Y_init_mgm2 = p$Y_init_mgm2
      )
    }

    # Routed inflows
    Qi_cell <- matrix(0.0, nrow=nday, ncol=ncell)
    Mi_cell <- matrix(0.0, nrow=nday, ncol=ncell)
    Qi_basin_mat <- matrix(0.0, nrow=nday, ncol=ncell)
    Mi_basin_mat <- matrix(0.0, nrow=nday, ncol=ncell)

    # Lagged seep recycle arrays
    QRecycle <- matrix(0.0, nrow=nday+1L, ncol=ncell)
    MRecycle <- matrix(0.0, nrow=nday+1L, ncol=ncell)

    for (nm in c("Q7","L7","Q3","L3","Q13","L13","Q14","L14","Q25","L25","Q26","L26","Q17","L17")) {
      case_out[[nm]][] <- 0.0
    }
    case_wb[,!(names(case_wb)%in%c("Date"))] <- NA_real_
    case_mb[,!(names(case_mb)%in%c("Date"))] <- NA_real_

    # Rolling depth history per cell for Z_plant (seed from initialized V)
    Z_hist <- matrix(0.0, nrow=nday, ncol=ncell)
    for (ic in seq_len(ncell)) {
      # Z_hist[1, ic] <- V[ic] / cells[[ic]]$params$A_cell
      Aic <- cells[[ic]]$params$A_cell
      Z_hist[1, ic] <- if (is.finite(Aic) && Aic > 0) V[ic] / Aic else NA_real_

    }

    # Daily loop
    for (day in seq_len(nday)) {

      V_start_total <- sum(V)

      P_cells_start <- sum(vapply(Pstate, function(ps) sum(ps$M) + sum(ps$S), 0.0))
      P_transit_start <- sum(MRecycle[day, ])
      P_transit_end   <- sum(MRecycle[day + 1L, ])

      RainVol_total <- 0.0
      EtVol_total   <- 0.0
      NetAtmo_total <- 0.0
      WB_in_total   <- 0.0
      WB_out_total  <- 0.0
      V_cell_day_total <- 0.0

      P_in_rain_total   <- 0.0
      P_in_drydep_total <- 0.0
      P_in_seepin_total <- 0.0

      P_mech_uptake_total  <- 0.0
      P_mech_recycle_total <- 0.0
      P_mech_sed_total     <- 0.0
      P_mech_direct_total  <- 0.0

      # Add basin inflows (pre-splitter)
      for (ic in seq_len(ncell)) {
        cd <- cells[[ic]]
        basin_ok <- (spl_idx == 0L) || (ic <= spl_idx)
        Qin_frac_eff <- offline_by_cell[[ic]]$Qin_frac_used[day]

        Qin_basin <- if (basin_ok) series$Qi[day] * Qin_frac_eff else 0.0
        Lin_basin <- if (basin_ok) (series$Qi[day] * series$Ci[day]) * Qin_frac_eff else 0.0

        Qi_cell[day, ic] <- Qi_cell[day, ic] + Qin_basin
        Mi_cell[day, ic] <- Mi_cell[day, ic] + Lin_basin

        Qi_basin_mat[day, ic] <- Qi_basin_mat[day, ic] + Qin_basin
        Mi_basin_mat[day, ic] <- Mi_basin_mat[day, ic] + Lin_basin
      }

      # Simulate each cell
      for (ic in seq_len(ncell)) {
        cd <- cells[[ic]]
        p  <- cd$params

        ## Lagged recycle
        RecycleQ <- sum(QRecycle[day, cells[[ic]]$RecycleIndex], na.rm = TRUE)
        RecycleM <- sum(MRecycle[day, cells[[ic]]$RecycleIndex], na.rm = TRUE)

        nz <- dmsta_zneighbors(day, series$Zcontrol)

        # IMPORTANT: DO NOT reinitialize V/Pstate here

        Qin_total <- Qi_cell[day, ic]
        Lin_total <- Mi_cell[day, ic]
        Cin_total <- if (Qin_total > 0) Lin_total / Qin_total else 0.0

        Qin_basin <- Qi_basin_mat[day, ic]
        Lin_basin <- Mi_basin_mat[day, ic]
        Cin_basin <- if (Qin_basin > 0) Lin_basin / Qin_basin else 0.0

        Qin_up <- max(0.0, Qin_total - Qin_basin)
        Lin_up <- Lin_total - Lin_basin
        Cin_up <- if (Qin_up > 0) Lin_up / Qin_up else 0.0

        # Lagged seep recycle into this cell
        # RecycleQ <- 0.0
        # RecycleM <- 0.0
        # for (j in seq_len(ncell)) {
        #   if (!is.null(cells[[j]]$RecycleIndex) && cells[[j]]$RecycleIndex == ic) {
        #     RecycleQ <- RecycleQ + QRecycle[day, j]
        #     RecycleM <- RecycleM + MRecycle[day, j]
        #   }
        # }

        nz <- dmsta_zneighbors(day, series$Zcontrol)

        day_inputs <- list(
          Date = series$Date[day],
          Qi   = Qin_total,
          Ci   = Cin_total,
          Rain = series$Rain[day],
          Et   = series$Et[day],
          Zcontrol = nz$today,
          Zcontrol_prev = nz$prev_day,
          Zcontrol_next = nz$nxt,
          # has_depth_constraint = has_Zcon_cell[ic],
          has_depth_constraint = has_depth_constraint_cell[ic],
          Qr0 = if (has_Qr0_cell[ic]) series$Qr0[day] else 0,
          Qr1 = if (has_Qr1_cell[ic]) series$Qr1[day] else 0,
          Qr2 = if (has_Qr2_cell[ic]) series$Qr2[day] else 0,
          RecycleQ = RecycleQ,
          RecycleM = RecycleM
        )

        ## Rolling Z_plant
        if (day > 1) Z_hist[day, ic] <- V[ic] / p$A_cell
        i0 <- max(1L, day - N_plant + 1L)
        Z_plant <- mean(Z_hist[i0:day, ic], na.rm = TRUE)

        ## Prevent double Qin scaling
        p_run <- p
        p_run$Qin_Frac <- 1.0

        ## Updated call with integrator pass‑through
        res <- dmsta_flowP_day(
          V = V[ic],
          P_state = Pstate[[ic]],
          tanks = tanks[[ic]],
          inputs = day_inputs,
          params = p_run,
          ppar = cd$ppar,
          constants = cd$constants,
          Qmethod = Qmethod,
          Pmethod = Pmethod,
          Nsteps = Nsteps,
          Z_plant = Z_plant,
          integrator_fun = integrator_fun,
          interp_option = interp_option,
          ...
        )

        # Advance states
        V[ic] <- res$results$V_end
        Pstate[[ic]] <- res$results$P_state_end
        Z_hist[day, ic] <- V[ic] / p$A_cell

        RainVol_total <- RainVol_total + res$budgets$water$RainVol
        EtVol_total   <- EtVol_total   + res$budgets$water$EtVol
        NetAtmo_total <- NetAtmo_total + res$budgets$water$NetAtmo
        WB_in_total   <- WB_in_total   + res$budgets$water$WB_in
        WB_out_total  <- WB_out_total  + res$budgets$water$WB_out
        V_cell_day_total <- V_cell_day_total + res$results$V_cell_day

        # DMSTA-style terms
        # Term 3: bypass
        Q3  <- .scalar0(res$results$Bypass,0)
        L3  <- .scalar0(res$results$P$loads$bypass,0)
        C3  <- .scalarC(res$results$P$conc$C_bypass, L3,Q3,default = 0)

        # Stream 7: into cell (treated inflow excluding bypass; excludes recycle)
        # (Matches the concept of "Qt(1)-Qt(3)" for Tank 1 inflow excluding bypass.)
        Q7 <- max(0.0, Qin_total - Q3)
        L7 <- Q7 * Cin_total
        C7 <- fw(L7, Q7)

        # Term 13: treated outflow (ONLY routed downstream)
        Q13 <- .scalar0(res$results$Q_treated, 0)
        L13 <- .scalar0(res$results$P$loads$treated, 0)
        C13 <- .scalarC(res$results$P$conc$C_treated, L13, Q13, default = 0)

        # Term 25/26: releases (out-of-system)
        Q25 <- .scalar0(res$results$Q_rel1, 0)
        L25 <- .scalar0(res$results$P$loads$rel1, 0)
        C25 <- .scalarC(res$results$P$conc$C_rel1, L25, Q25, default = 0)

        Q26 <- .scalar0(res$results$Q_rel2, 0)
        L26 <- .scalar0(res$results$P$loads$rel2, 0)
        C26 <- .scalarC(res$results$P$conc$C_rel2, L26, Q26, default = 0)

        # Term 14: seepage discharge out-of-system
        Q14 <- .scalar0(res$results$P$flows$seep_discharge, 0)
        L14 <- .scalar0(res$results$P$loads$seep_discharge, 0)
        C14 <- .scalarC(res$results$P$conc$C_seep_discharge, L14, Q14, default = 0)

        # Term 17: seep recycle bookkeeping (lagged)
        Q17 <- .scalar0(res$results$P$flows$seep_recycle, 0)
        L17 <- .scalar0(res$results$P$loads$seep_recycle, 0)
        C17 <- .scalarC(res$results$P$conc$C_seep_recycle, L17, Q17, default = 0)

        pb <- res$budgets$mass
        if (is.null(pb) || is.null(pb$storage) || is.null(pb$closure)) {
          stop("dmsta_flowP_day() did not return a valid P_budget.")
        }

        # External inputs for case budget (rain + drydep + seepin)
        P_in_rain_total   <- P_in_rain_total   + pb$inputs_external$L_rain
        P_in_drydep_total <- P_in_drydep_total + pb$inputs_external$L_drydep
        P_in_seepin_total <- P_in_seepin_total + pb$inputs_external$L_seepin

        # Mechanism totals
        P_mech_uptake_total  <- P_mech_uptake_total  + pb$mechanisms$L_uptake
        P_mech_recycle_total <- P_mech_recycle_total + pb$mechanisms$L_recycle
        P_mech_sed_total     <- P_mech_sed_total     + pb$mechanisms$L_sed
        P_mech_direct_total  <- P_mech_direct_total  + pb$mechanisms$L_direct

        # Store per-cell series
        if (return_cell_series) {
          co   <- cell_out[[ic]]
          c_wb <- cell_wb[[ic]]
          c_mb <- cell_mb[[ic]]

          co$V_end[day]      <- res$results$V_end
          co$Z_end[day]      <- res$results$Z_end
          co$Z_avg[day]      <- res$results$Z_avg
          co$V_cell_day[day] <- res$results$V_cell_day

          co$Qin[day] <- Qin_total
          co$Lin[day] <- Lin_total
          co$Cin[day] <- Cin_total

          co$Qin_basin[day] <- Qin_basin
          co$Lin_basin[day] <- Lin_basin
          co$Cin_basin[day] <- Cin_basin

          co$Qin_up[day] <- Qin_up
          co$Lin_up[day] <- Lin_up
          co$Cin_up[day] <- Cin_up

          co$Q7[day] <- Q7;   co$L7[day] <- L7;   co$C7[day] <- C7
          co$Q3[day] <- Q3;   co$L3[day] <- L3;   co$C3[day] <- C3
          co$Q13[day] <- Q13; co$L13[day] <- L13; co$C13[day] <- C13
          co$Q14[day] <- Q14; co$L14[day] <- L14; co$C14[day] <- C14
          co$Q25[day] <- Q25; co$L25[day] <- L25; co$C25[day] <- C25
          co$Q26[day] <- Q26; co$L26[day] <- L26; co$C26[day] <- C26

          if (keep_Q17) {
            co$Q17[day] <- Q17; co$L17[day] <- L17; co$C17[day] <- C17
          }

          offrow <- offline_by_cell[[ic]][day, ]
          co$Qin_frac_base[day] <- offrow$Qin_frac_base
          co$Qin_frac_used[day] <- offrow$Qin_frac_used
          co$offline[day]       <- offrow$offline
          co$offline_index[day] <- offrow$offline_index
          co$off_diff[day]      <- offrow$off_diff
          co$off_mod[day]       <- offrow$off_mod

          c_wb$RainVol[day] <- res$budgets$water$RainVol
          c_wb$EtVol[day]   <- res$budgets$water$EtVol
          c_wb$NetAtmo[day] <- res$budgets$water$NetAtmo
          c_wb$WB_in[day]   <- res$budgets$water$WB_in
          c_wb$WB_out[day]  <- res$budgets$water$WB_out
          c_wb$WB_err[day]  <- res$budgets$water$WB_err
          c_wb$WB_rel[day]  <- res$budgets$water$WB_rel

          c_mb$dP[day] <- .scalarNA(pb$storage$dP)
          c_mb$Pin_external[day]  <- .scalarNA(pb$closure$Pin_external)
          c_mb$Pout_external[day] <- .scalarNA(pb$closure$Pout_external)
          c_mb$Perr_external[day] <- .scalarNA(pb$closure$Perr_external)
          c_mb$Prel_external[day] <- .scalarNA(pb$closure$Prel_external)
          c_mb$Pin_total[day]     <- .scalarNA(pb$closure$Pin_total)
          c_mb$Pout_total[day] <- .scalarNA(pb$closure$Pout_total)
          c_mb$Perr_total[day] <- .scalarNA(pb$closure$Perr_total)
          c_mb$Prel_total[day] <- .scalarNA(pb$closure$Prel_total)
          c_mb$Q_in_tanks[day] <- .scalarNA(pb$inflow_tanks$Q_in_tanks)
          c_mb$L_in_tanks[day] <- .scalarNA(pb$inflow_tanks$L_in_tanks)
          c_mb$C_in_tanks[day] <- .scalarNA(pb$inflow_tanks$C_in_tanks)
          c_mb$L_rain[day] <- .scalarNA(pb$inputs_external$L_rain)
          c_mb$L_drydep[day] <- .scalarNA(pb$inputs_external$L_drydep)
          c_mb$L_seepin[day] <- .scalarNA(pb$inputs_external$L_seepin)
          c_mb$L_treated[day] <- .scalarNA(res$results$P$loads$treated) # pb$outputs_external$L_treated
          c_mb$L_rel1[day] <- .scalarNA(res$results$P$loads$rel1) # pb$outputs_external$L_rel1
          c_mb$L_rel2[day] <- .scalarNA(res$results$P$loads$rel2) # pb$outputs_external$L_rel2
          c_mb$L_bypass[day] <- .scalarNA(res$results$P$loads$bypass) # pb$outputs_external$L_bypass
          c_mb$L_seep_discharge[day] <- .scalarNA(res$results$P$loads$seep_discharge) #pb$outputs_external$L_seep_discharge
          c_mb$L_seep_recycle_out[day] <- .scalarNA(res$results$P$loads$seep_recycle) # pb$transfers$L_seep_recycle_out
          c_mb$L_uptake[day] <- .scalarNA(pb$mechanisms$L_uptake)
          c_mb$L_recycle[day] <- .scalarNA(pb$mechanisms$L_recycle)
          c_mb$L_sed[day] <- .scalarNA(pb$mechanisms$L_sed)
          c_mb$L_direct[day] <- .scalarNA(pb$mechanisms$L_direct)
          c_mb$L_burial[day] <- .scalarNA(pb$mechanisms$L_burial)
          c_mb$L_release[day] <- .scalarNA(pb$mechanisms$L_release)

          cell_out[[ic]] <- co
          cell_wb[[ic]]  <- c_wb
          cell_mb[[ic]]  <- c_mb
        }

        # Case-level accumulation
        # Raw watershed inflow already in case_out$Qin/Lin/Cin
        # Accumulate "into cell" (Stream 7) using basin portion minus bypass (audit-style)
        # This mirrors the DMSTA "into cell" concept of external treated inflow (Qt(43)-Qt(3)),
        # but we clip at 0 for physical interpretability.
        case_out$Q7[day] <- case_out$Q7[day] + max(0.0, Qin_basin - Q3)
        case_out$L7[day] <- case_out$L7[day] + max(0.0, Qin_basin - Q3) * Cin_total

        case_out$Q3[day]  <- case_out$Q3[day]  + Q3
        case_out$L3[day]  <- case_out$L3[day]  + L3
        case_out$Q25[day] <- case_out$Q25[day] + Q25
        case_out$L25[day] <- case_out$L25[day] + L25
        case_out$Q26[day] <- case_out$Q26[day] + Q26
        case_out$L26[day] <- case_out$L26[day] + L26
        case_out$Q14[day] <- case_out$Q14[day] + Q14
        case_out$L14[day] <- case_out$L14[day] + L14
        if (keep_Q17) {
          case_out$Q17[day] <- case_out$Q17[day] + Q17
          case_out$L17[day] <- case_out$L17[day] + L17
        }

        if (spl_idx > 0L && ic == spl_idx) {
          frac <- cd$SplitterFrac
          if (is.null(frac)) stop("Splitter cell must provide SplitterFrac.")
          if (!is.null(names(frac))) {
            for (nm in names(frac)) {
              j <- as.integer(nm)
              f <- frac[[nm]]
              if (is.finite(f) && f != 0) {
                Qi_cell[day, j] <- Qi_cell[day, j] + f * Q13
                Mi_cell[day, j] <- Mi_cell[day, j] + f * L13
              }
            }
          } else {
            for (j in seq_len(ncell)) {
              f <- frac[j]
              if (is.finite(f) && f != 0) {
                Qi_cell[day, j] <- Qi_cell[day, j] + f * Q13
                Mi_cell[day, j] <- Mi_cell[day, j] + f * L13
              }
            }
          }
        } else {
          dcell <- cd$DownCell
          if (!is.null(dcell) && dcell > 0) {
            Qi_cell[day, dcell] <- Qi_cell[day, dcell] + Q13
            Mi_cell[day, dcell] <- Mi_cell[day, dcell] + L13
          } else {
            # terminal treated outflow leaves system
            case_out$Q13[day] <- case_out$Q13[day] + Q13
            case_out$L13[day] <- case_out$L13[day] + L13
          }
        }

        # Seep recycle arrays (lagged by 1 day)
        QRecycle[day + 1L, ic] <- Q17
        MRecycle[day + 1L, ic] <- L17


      } # end cell loop

      ## Finish case-level diagnostics for the day
      ## Concentrations for case terms
      case_out$C7[day]  <- fw(case_out$L7[day],  case_out$Q7[day])
      case_out$C3[day]  <- fw(case_out$L3[day],  case_out$Q3[day])
      case_out$C13[day] <- fw(case_out$L13[day], case_out$Q13[day])
      case_out$C14[day] <- fw(case_out$L14[day], case_out$Q14[day])
      case_out$C25[day] <- fw(case_out$L25[day], case_out$Q25[day])
      case_out$C26[day] <- fw(case_out$L26[day], case_out$Q26[day])
      if (keep_Q17) case_out$C17[day] <- fw(case_out$L17[day], case_out$Q17[day])

      ## Case-level hydrology totals
      V_end_total <- sum(V, na.rm = TRUE)
      dV_total <- V_end_total - V_start_total
      ## Case-level WB based on summing cell budgets (internal routed treated flows cancel)
      WB_err_total <- dV_total - (WB_in_total - WB_out_total)
      WB_rel_total <- WB_err_total / max(1e-12, max(WB_in_total, WB_out_total))
      # case-level water budget
      case_wb$RainVol_total[day] <- RainVol_total
      case_wb$EtVol_total[day]   <- EtVol_total
      case_wb$NetAtmo_total[day] <- NetAtmo_total
      case_wb$WB_in_total[day]   <- WB_in_total
      case_wb$WB_out_total[day]  <- WB_out_total
      case_wb$WB_err_total[day]  <- WB_err_total
      case_wb$WB_rel_total[day]  <- WB_rel_total
      case_wb$V_total_end[day]   <- V_end_total
      case_wb$dV_total[day]      <- dV_total
      case_wb$Z_total_end[day]   <- if (A_total > 0) V_end_total / A_total else NA_real_
      case_wb$Z_total_avg[day]   <- if (A_total > 0) V_cell_day_total / A_total else NA_real_

      ## Case-level external P budget
      # Update transit end now that MRecycle[day+1,] is fully assigned
      P_transit_end <- sum(MRecycle[day + 1L, ])
      dP_transit <- P_transit_end - P_transit_start
      # P in cells at end of day
      P_cells_end <- sum(vapply(Pstate, function(ps) sum(ps$M) + sum(ps$S), 0.0))
      dP_cells <- P_cells_end - P_cells_start
      # External P inputs (kg): watershed + rain + drydep + seepin
      P_in_watershed <- series$Qi[day] * series$Ci[day]
      P_in_external <- P_in_watershed + P_in_rain_total + P_in_drydep_total + P_in_seepin_total
      # External P outputs (kg): bypass + terminal treated + releases + seep discharge
      P_out_bypass <- case_out$L3[day]
      P_out_terminal_treated <- case_out$L13[day]
      P_out_releases <- case_out$L25[day] + case_out$L26[day]
      P_out_seep_discharge <- case_out$L14[day]
      P_out_external <- P_out_bypass + P_out_terminal_treated + P_out_releases + P_out_seep_discharge
      # Closure includes transit (internal lagged transfer reservoir)
      P_err_case <- (dP_cells + dP_transit) - (P_in_external - P_out_external)
      P_rel_case <- P_err_case / max(1e-12, max(P_in_external, P_out_external))
      # Store case P terms
      case_mb$P_in_watershed[day] <- P_in_watershed
      case_mb$P_in_rain_total[day] <- P_in_rain_total
      case_mb$P_in_drydep_total[day] <- P_in_drydep_total
      case_mb$P_in_seepin_total[day] <- P_in_seepin_total
      case_mb$P_in_external[day] <- P_in_external
      case_mb$P_out_bypass[day] <- P_out_bypass
      case_mb$P_out_terminal_treated[day] <- P_out_terminal_treated
      case_mb$P_out_releases[day] <- P_out_releases
      case_mb$P_out_seep_discharge[day] <- P_out_seep_discharge
      case_mb$P_out_external[day] <- P_out_external
      case_mb$P_transit_start[day] <- P_transit_start
      case_mb$P_transit_end[day] <- P_transit_end
      case_mb$dP_transit[day] <- dP_transit
      case_mb$P_cells_start[day] <- P_cells_start
      case_mb$P_cells_end[day] <- P_cells_end
      case_mb$dP_cells[day] <- dP_cells
      case_mb$P_err_case[day] <- P_err_case
      case_mb$P_rel_case[day] <- P_rel_case
      case_mb$P_mech_uptake_total[day] <- P_mech_uptake_total
      case_mb$P_mech_recycle_total[day] <- P_mech_recycle_total
      case_mb$P_mech_sed_total[day] <- P_mech_sed_total
      case_mb$P_mech_direct_total[day] <- P_mech_direct_total

    } # end day loop

    # Convergence metric
    totalL <- sum(case_out$L3 + case_out$L13 + case_out$L14 + case_out$L25 + case_out$L26, na.rm = TRUE)
    conv_test <- if (!is.na(prev_totalL) && prev_totalL > 0) abs(totalL - prev_totalL) / prev_totalL else NA_real_

    conv_hist <- rbind(conv_hist, data.frame(iter=iter, totalL=totalL, conv_test=conv_test))
    iterations_used <- iter

    if (is.finite(conv_test) && conv_test < conv_tol) {
      converged <- TRUE
      break
    }
    prev_totalL <- totalL
  } # end iteration loop

  # Drop Q17 if requested
  if (!keep_Q17) {
    drop_cols <- c("Q17","C17","L17")
    case_out <- case_out[, setdiff(names(case_out), drop_cols), drop=FALSE]
    if (return_cell_series) {
      cell_out <- lapply(cell_out, function(df) df[, setdiff(names(df), drop_cols), drop=FALSE])
    }
  }

  ## Build case-level results from raw pieces (case_out + budgets)
  stopifnot(all(as.Date(case_out$Date) == as.Date(case_wb$Date)))
  if (!is.null(case_mb)) stopifnot(all(as.Date(case_out$Date) == as.Date(case_mb$Date)))

  case_raw <- cbind(
    case_out,
    case_wb[, setdiff(names(case_wb), "Date"), drop = FALSE]
  )

  case_std <- dmsta_case_components(
    case_raw,
    keep_Q17 = keep_Q17,
    keep_P   = TRUE,
    ...
  )

  if (return_cell_series && !is.null(cell_out)) {
    cell_out <- lapply(cell_out, dmsta_case_components, keep_Q17 = keep_Q17, keep_P = TRUE, ...)
  }

  meta <- list(
    convergence = list(
      converged = converged,
      iterations_used = iterations_used,
      conv_tol = conv_tol,
      history = conv_hist
    ),
    Qmethod = Qmethod,
    Pmethod = Pmethod,
    Nsteps = Nsteps,
    N_plant = N_plant,
    max_iter = max_iter,
    conv_tol = conv_tol,
    cells = cells,
    offline_qin_frac = offline_by_cell
  )

  out <- list(
    results = list(
      case  = case_std,
      cells = cell_out
    ),
    budgets = list(
      water = list(
        case  = case_wb,
        cells = cell_wb
      ),
      mass = list(
        case  = case_mb,
        cells = cell_mb,
        transit = NULL
      )
    ),
    meta = meta
  )

  class(out) <- c("dmsta_network_result", "list")
  out
}
