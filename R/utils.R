
#' Calculate julian day from Date
#'
#' @param date
#'
#' @return Numeric vector of julian day
#'
#' @keywords internal
#' @noRd
julian_day <- function(date) {
  as.integer(format(as.Date(date), "%j"))
}

#' Flow weighted concentration calculator
#'
#' @param L TP load as kg/d
#' @param Q flow as hm3/day
#'
#' @return Numeric vector of concentration in ug/L or mg/m3
#'
#' @keywords internal
#' @noRd
fw <- function(L, Q){
  ## internal helper function
  ifelse(Q > 0, L / Q, NA_real_)
}



#' Convert Cubic Feet per Second to Cubic Hectometers per Day
#'
#' Converts volumetric flow from cubic feet per second (cfs)
#' to cubic hectometers per day (hm^3/day).
#'
#' @param x Numeric vector of flow values in cubic feet per second.
#' @return Numeric vector of flow values in cubic hectometers per day.
#' @export

cfs_to_hm3d <- function(x) {
  # m3_per_ft3  <- 0.3048^3
  # m3_per_day  <- m3_per_ft3 * 86400
  # hm3_per_day <- m3_per_day / 1e6
  # x * hm3_per_day
  x * 0.002448455 # consistent with VBA
}



#' centimeters to meters
#'
#' @param x numeric value
#'
#' @return converted numeric value
#' @keywords internal
#' @export
cm_to_m <- function(x){x/100}

#' meters to centimeters
#'
#' @param x numeric value
#'
#' @return converted numeric value
#' @keywords internal
#' @export
m_to_cm <- function(x){x*100}

#' inches to meters
#'
#' @param x numeric value
#'
#' @return converted numeric value
#' @keywords internal
#' @export
in_to_m <- function(x){x/39.37}


#' Safely Resolve a Numeric Scalar
#'
#' Returns a numeric value if it is non-NULL and finite.
#' Otherwise returns a specified default.
#'
#' This helper is used internally to guard against
#' NULL, NA, NaN, or infinite values in model calculations.
#'
#' @param x A candidate numeric value.
#' @param default A scalar value returned if `x` is NULL or not finite.
#'
#' @return A single numeric value.
#'
#' @keywords internal
#' @noRd
safe_num <- function(x, default = 0) {
  if (is.null(x) || !is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    return(default)
  }
  x
}

#' Safely Resolve a Numeric vector
#'
#' Returns a numeric value if it is non-NULL and finite.
#' Otherwise returns a specified default.
#'
#' This helper is used internally to guard against
#' NULL, NA, NaN, or infinite values in model calculations.
#'
#' @param x A candidate numeric value.
#' @param default A scalar value returned if `x` is NULL or not finite.
#'
#' @return A single numeric value.
#'
#' @keywords internal
#' @noRd
safe_num_vec <- function(x, default = 0){
  if (is.null(x)) return(default)
  x[!is.finite(x)] <- default
  x
}

#' Return First Non-NULL Finite Value
#'
#' Scans a list of candidate values and returns the first element
#' that is not `NULL` and is a finite numeric value.
#' If no such value is found, the supplied `default` is returned.
#'
#' This function is primarily used internally for resolving
#' parameter values across multiple model slots.
#'
#' @param values A list or vector of candidate values.
#' @param default A scalar value returned if no finite value is found.
#'
#' @return A single numeric value.
#'
#' @keywords internal
#' @noRd
first_non_null_finite <- function(values, default) {
  for (v in values) {
    if (!is.null(v) && is.numeric(v) && length(v) == 1L && is.finite(v)) {
      return(v)
    }
  }
  default
}

#' Retrieve Current and Neighboring Control Elevations
#'
#' Returns the control elevation at index `i` along with the
#' previous and next values in the vector `z`.
#'
#' For boundary indices (first or last element), the current value
#' is repeated in place of the missing neighbor to avoid
#' out-of-bounds indexing.
#'
#' This function is used internally to provide stable
#' neighbor-aware control elevations for hydrologic integration.
#'
#' @param i Integer index into `z`.
#' @param z Numeric vector of control elevations.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{today}{Value at index `i`.}
#'   \item{prev_day}{Value at index `i - 1`, or `today` if `i == 1`.}
#'   \item{nxt}{Value at index `i + 1`, or `today` if `i == length(z)`.}
#' }
#'
#' @keywords internal
#' @noRd
neighbors_zcontrol <- function(i, z) {
  zi  <- z[i]
  zim <- if (i > 1) z[i - 1] else zi
  zip <- if (i < length(z)) z[i + 1] else zi
  list(today = zi, prev_day = zim, nxt = zip)
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



#' Create a default DMSTAr parameter list
#'
#' Constructs a named list of model parameters used by DMSTAr. Most arguments
#' are scalar numeric values (or logical flags). Additional named parameters
#' may be supplied via `...`; these will be appended to the returned list
#' and will override any existing defaults with the same name.
#'
#' @details
#' \itemize{
#'   \item **Derived `IsaNode`:** If `IsaNode` is `NULL`,
#'   it is set to `(A_cell <= 0)`; otherwise the provided `IsaNode`
#'   value is used.
#'   \item **`MT` behavior:** When `MT = TRUE`, all numeric
#'   values in the returned parameter list are replaced with `0`,
#'   providing an empty list. Logical and non-numeric entries are left unchanged.
#'   \item **Extra parameters via `...`:** All arguments passed in
#'   `...` must be named. New names trigger a warning and are added to the
#'   returned list. Names that already exist in the defaults override the default
#'   value.
#' }
#'
#' @param MT Logical. If `TRUE`, replaces all numeric parameter values in
#'   the returned list with `0`. Defaults to `FALSE`.
#' @param Qin_Frac Numeric. Fraction of flow entering cell from basin.
#' @param A_cell Numeric. Cell effective treatment are (units: km2).
#' @param Width Numeric. mean width of flow path (units: km.
#' @param Ntanks Integer-ish numeric. Number of tanks in series.
#' @param Zrelease Numeric. Release elevation (units: cm).
#' @param Q_zmin Numeric. Zc; no outflow below this water level , added to control depth specified in input series (units: cm).
#' @param Zweir Numeric. Zw; fixed weir depth, use for reservoirs with outflow hydraulics controlled by outlet structure; = 0 for shallow systems when outflow is usually controlled by vegetation resistance (units: cm).
#' @param Q_b Numeric. Rating curve parameter b; q / w = a  (Z - Zw)  ^b  for Z >= Zc; 'typicallly ~  3 to 4 for marsh control; 1.5 for weir control (reservoirs & other deep cells).
#' @param Q_a Numeric. Rating curve parameter a; flow/width at water depth of 1 m (a); typically ~ 0.5 to 2.
#' @param Bypass_elev Numeric. depth at which bypass begins ( 0 = no limit ).
#' @param Qomax Numeric. Inflow capacity (triggers bypass)  ( 0 = no limit ) (units: hm3/d).
#' @param Qimax Numeric. Outflow capacity (0 = no limit)  (units: hm3/d).
#' @param Seepin_Rate Numeric. centimeters per day per centimeter of head, reflects transmissivity of soils (units: cm/d/cm).
#' @param Seepin_Elev Numeric. 'drives inflow seepage rate, depth relative to mean ground surface elev (units: cm).
#' @param seepin_conc Numeric. Concentration associated with seepage inflow (units: mg/m3).
#' @param Seepout_Rate Numeric. centimeters per day per centimeter of head, reflects transmissivity of soils (units: cm/d/cm).
#' @param Seepout_Elev Numeric. drives outflow seepage rate, relative to mean ground surface elev, can be < 0 (units: cm).
#' @param seepage_c Numeric. Seepage concentration term (units: mg/m3).
#' @param C_init_ppb Numeric. Initial concentration (units: mg/m3).
#' @param Y_init_mgm2 Numeric. Initial biomass P storage (units: mg/m^2).
#' @param Zinit Numeric. Initial water column depth relative to mean ground elevation (units: cm).
#' @param Cstar Numeric. Reference concentration parameter also defined as water column conc at storage = 0 mg/m2 at steady-state (units: mg/m3).
#' @param C1000 Numeric. water column conc at storage = 1000 mg/m2 at steady-state (units: mg/m3).
#' @param Chalf Numeric. Half-saturation / half-response concentration parameter or water column concentration at 1/2 maximum uptake (units: mg/m3).
#' @param Ks_per_yr Numeric. net settling rate in steady state in K/C* model (units: m/yr).
#' @param Z1 Numeric. Uptake rate decreases below this depth; flat between Z1 and Z2; =0 ignored (units: cm).
#' @param Z2 Numeric. uptake rate starts to decrease above this depth (30-day average); decreases linearly between Z2 & Z3; 0 = ignored; reflects damage to vegetation at high depths (units: cm)
#' @param Z3 Numeric. upper end of depth penalty range; K = 1 m/yr at depths above this value, regardless of calibration; 0 = ignored (units: cm).
#' @param K2Coef1 Numeric. Coefficient for secondary rate term (model-specific).
#' @param SeasonalFactor Numeric. Seasonal factor multiplier (model-specific).
#' @param Ytrans Numeric. Transform parameter for Y term (model-specific).
#' @param Ysigma Numeric. Sigma/spread parameter for Y term (model-specific).
#' @param Czero Numeric. Baseline concentration offset (model-specific).
#' @param ks_2 Numeric. Secondary coefficient (model-specific).
#' @param zh_2 Numeric. Secondary depth/elevation parameter (model-specific).
#' @param k_depth_penalty Numeric. Depth penalty multiplier/parameter.
#' @param C_rain Numeric. Rain concentration (model-specific units).
#' @param DryDepo Numeric. Dry deposition loading term (model-specific units).
#' @param ShutdownET Logical. If `TRUE`, ET is shut down under model-specific conditions.
#' @param force_Q_out Logical. If `TRUE`, forces outflow behavior (model-specific).
#' @param DutyCycle Numeric. Duty cycle (0–1) applied to relevant process(es).
#' @param Zmin Numeric. Minimum elevation (model-specific).
#' @param Cmax Numeric. Maximum concentration cap (model-specific).
#' @param IsaNode Logical or `NULL`. If `NULL`, derived as `(A_cell <= 0)`.
#' @param ... Additional named parameters to add to the returned list, or to
#'   override existing defaults by name. All must be named (e.g., `foo = 1`).
#'
#' @return A named list of DMSTAr parameters.
#'
#' @examples
#' # Get defaults
#' p <- dmstar_default_params()
#'
#' # Override a few defaults
#' p <- dmstar_default_params(Zmin = 5, DutyCycle = 0.9, Ntanks = 3)
#'
#' # Add a new parameter via ...
#' p <- dmstar_default_params(NewParam = 123)
#'
#' # Derive IsaNode automatically when IsaNode is NULL
#' p <- dmstar_default_params(A_cell = 0)   # IsaNode becomes TRUE
#' p <- dmstar_default_params(A_cell = 10)  # IsaNode becomes FALSE
#'
#' # MT = TRUE zeroes numeric parameters
#' p0 <- dmstar_default_params(MT = TRUE)
#'
#' @export

dmstar_default_params <- function(MT = FALSE,
                                  Qin_Frac = 0,
                                  A_cell = 0,
                                  Width = 0,
                                  Ntanks = 1,
                                  Zrelease = 0,
                                  Q_zmin = 0,
                                  Zweir = 0,
                                  Q_b = 0,
                                  Q_a = 0,
                                  Bypass_elev = 0,
                                  Qomax = 0,
                                  Qimax = 0,
                                  Seepin_Rate = 0,
                                  Seepin_Elev = 0,
                                  seepin_conc = 0,
                                  Seepout_Rate = 0,
                                  Seepout_Elev = 0,
                                  seepage_c = 20,
                                  C_init_ppb = 0,
                                  Y_init_mgm2 = 0,
                                  Zinit = 40,
                                  Cstar = 3,
                                  C1000 = 22,
                                  Chalf = 300,
                                  Ks_per_yr = 0,
                                  Z1 = 40,
                                  Z2 = 100,
                                  Z3 = 200,
                                  K2Coef1 = 0,
                                  SeasonalFactor = 0,
                                  Ytrans = 0,
                                  Ysigma = 0,
                                  Czero = 0,
                                  ks_2 = 0,
                                  zh_2 = 0,
                                  k_depth_penalty = 1,
                                  C_rain = 10,
                                  DryDepo = 20,
                                  ShutdownET = TRUE,
                                  force_Q_out = FALSE,
                                  DutyCycle = 0.95,
                                  Zmin = 2,
                                  Cmax = 2000,
                                  IsaNode = NULL,
                                  ...
){

  # Build parameter list
  param_temp <-
    list(A_cell = A_cell,
         Zmin = Zmin,
         Zinit = Zinit,
         Zweir = Zweir,
         Q_zmin = Q_zmin,
         Zrelease = Zrelease,
         Bypass_elev = Bypass_elev,
         Q_a = Q_a,
         Q_b = Q_b,
         Width = Width,
         Qomax = Qomax,
         Qimax = Qimax,
         Seepout_Rate = Seepout_Rate, Seepout_Elev = Seepout_Elev,
         Seepin_Rate = Seepin_Rate, Seepin_Elev = Seepin_Elev,
         ShutdownET = ShutdownET,
         force_Q_out = force_Q_out,
         DutyCycle = DutyCycle,
         Cmax = Cmax,
         C1000 = C1000,
         Cstar = Cstar,
         Ks_per_yr = Ks_per_yr,
         Z1 = Z1,
         Z2 = Z2,
         Z3 = Z3,
         Chalf = Chalf,
         K2Coef1 = K2Coef1,
         SeasonalFactor = SeasonalFactor,
         Ytrans = Ytrans,
         Ysigma = Ysigma,
         Czero = Czero,
         ks_2 = ks_2,
         zh_2 = zh_2,
         k_depth_penalty = k_depth_penalty,
         C_rain = C_rain,
         DryDepo = DryDepo,
         seepage_c = seepage_c,
         seepin_conc = seepin_conc,
         C_init_ppb = C_init_ppb,
         Y_init_mgm2 = Y_init_mgm2,
         Qin_Frac = Qin_Frac,
         Ntanks = Ntanks)


  # derive IsaNode if not supplied
  if (is.null(IsaNode)) {
    param_temp$IsaNode <- (A_cell <= 0)
  } else {
    param_temp$IsaNode <- IsaNode
  }

  if (isTRUE(MT)) {
    param_temp <- lapply(param_temp, function(x) {
      if (is.numeric(x)) 0 else x
    })
  }

  if (is.null(param_temp$Zrelease) || !is.finite(param_temp$Zrelease)) {
    param_temp$Zrelease <- 0
  }


  dots <- list(...)
  # enforce that extras are named (prevents accidental positional args)
  if (length(dots) && (is.null(names(dots)) || any(names(dots) == ""))) {
    stop("All arguments in ... must be named (e.g., new_param = 123).")
  }
  if (length(dots)) {
    new_names <- setdiff(names(dots), names(param_temp))
    if (length(new_names)) {
      warning("Adding new parameter(s): ", paste(new_names, collapse = ", "))
    }
  }
  # merge: overrides existing + adds new
  param_temp[names(dots)] <- dots

  param_temp
}



