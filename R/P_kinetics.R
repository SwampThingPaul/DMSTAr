#' Internal phosphorus kinetics helpers
#' @name internal_p_kinetics
#' @keywords internal
NULL


#' Compute DMSTAr kinetic coefficients (K1, K2, K3) and kinetic model flag
#'
#' Internal helper that derives kinetic coefficients from `C1000`,
#' `Cstar`, and `Ks`. Two formulations are supported:
#'
#' \itemize{
#'   \item **PModel 1 (standard)**: used when `C1000 > 0`, `Ks > 0`,
#'     and `Cstar > 0`.
#'   \item **PModel 2 (transformed/special)**: triggered when `Cstar < 0`,
#'     interpreted as `G = -Cstar`.
#' }
#'
#' If inputs do not define a valid kinetic model, the function returns
#' `K1 = K2 = K3 = 0` and `PModel = NA_integer_`.
#'
#' @param C1000 Numeric scalar. Concentration-like parameter (units consistent
#'   with the DMSTAr kinetic formulation).
#' @param Cstar Numeric scalar. Kinetic threshold parameter. If negative,
#'   the transformed kinetics (PModel 2) are used.
#' @param Ks Numeric scalar. Rate parameter on the time step of interest.
#'
#' @return A named `list` with elements:
#'   \describe{
#'     \item{K1}{Numeric scalar.}
#'     \item{K2}{Numeric scalar.}
#'     \item{K3}{Numeric scalar.}
#'     \item{PModel}{Integer scalar: `1L`, `2L`, or `NA_integer_`.}
#'   }
#'
#' @keywords internal
#' @rdname internal_p_kinetics

compute_DMSTA_kvals <- function(C1000, Cstar, Ks) { # internal function
  # if (Cstar >= C1000 && Ks > 0) {
  #   warning("Invalid parameter set: Cstar >= C1000")
  # }
  # default
  K1 <- K2 <- K3 <- 0
  PModel <- NA_integer_

  # PModel 2: special / transformed kinetics
  if (!is.na(Cstar) && Cstar < 0) {

    G <- -Cstar

    if (G > 0 && Ks > 0 && C1000 > 0) {
      K1 <- Ks / G
      K3 <- Ks * C1000 / 1000
      K2 <- (1 - G) / G * K3
    }

    PModel <- 2
    return(list(K1 = K1, K2 = K2, K3 = K3, PModel = PModel))
  }

  # PModel 1: standard Module 1
  if (C1000 > 0 && Ks > 0 && Cstar > 0) {

    # faithful translation of VB logic
    f <- max(1, C1000 - Cstar) / 1000

    K3 <- Ks * f
    K1 <- K3 / Cstar
    K2 <- (K3 * K1) / Ks

    PModel <- 1
  }

  list(K1 = K1, K2 = K2, K3 = K3, PModel = PModel)
}

#' Build phosphorus kinetics parameters for a registered model
#'
#' Builds a standardized parameter list for a selected phosphorus model type
#' (e.g., `"STA"`, `"PSTA"`, `"RES"`) using a registered
#' model-builder function, then derives kinetic coefficients `K1`, `K2`,
#' `K3` and the kinetic model indicator `PModel` via
#' `compute_DMSTA_kvals()`.
#'
#' Model builders are obtained from the package phosphorus model registry.
#'
#' @param mod_type Character scalar. Model identifier (e.g., `"STA"`).
#' @param Dpy Numeric scalar. Time steps per year (default `365.25`).
#' @param DutyCycle Numeric scalar in \eqn{[0, 1]}. Duty-cycle multiplier applied
#'   to time-scaled rate parameters.
#' @param pparams Named `list`. Raw parameters expected by the chosen model
#'   builder.
#' @param ... Additional arguments passed through to the underlying model builder.
#'
#' @details
#' The selected model builder should return a list containing at minimum
#' `C1000`, `Cstar`, and `Ks` (on the target time step).
#' Additional fields returned by the builder (e.g., `Z1`, `Z2`, `Chalf`)
#' are preserved. The function appends `K1`, `K2`, `K3`, and `PModel`.
#'
#' The result is assigned class `"P_kinetics"` (prepended to any existing classes).
#'
#' @return A named `list` of standardized parameters with elements
#' `K1`, `K2`, `K3`, and `PModel`. The result has class
#' `"P_kinetics"`.
#'
#' @examples
#' # Example assumes DMSTAr ships with a registered "STA" builder.
#' pparams <- list(
#'   C1000 = 1, Cstar = 0.2, Ks_per_yr = 0.5,
#'   Z1 = 10, Z2 = 30, Z3 = 60,
#'   K2Coef1 = 0.1, Chalf = 0.2, SeasonalFactor = 1,
#'   DutyCycle = 0.95
#' )
#' pars <- build_P_kinetics("STA", Dpy = 365.25, pparams = pparams)
#' pars$K1
#' pars$PModel
#'
#' @export
build_P_kinetics <- function(mod_type, Dpy = 365.25, DutyCycle = NULL, pparams, ...) {
  mod_type <- toupper(trimws(mod_type))

  builder <- .P_MODEL_BUILDERS_DEFAULT[[mod_type]]
  if (is.null(builder)) {
    stop("Unknown P model type: ", mod_type,
         ". Valid: ", paste(names(.P_MODEL_BUILDERS_DEFAULT), collapse = ", "),
         call. = FALSE)
  }

  DutyCycle <- pparams$DutyCycle

  pars <- builder(Dpy = Dpy, DutyCycle = DutyCycle,pparams = pparams, ...)

  K <- compute_DMSTA_kvals(
    C1000 = pars$C1000,
    Cstar = pars$Cstar,
    Ks    = pars$Ks
  )

  pars$K1 <- K$K1
  pars$K2 <- K$K2
  pars$K3 <- K$K3
  pars$PModel <- K$PModel

  class(pars) <- c("P_kinetics", class(pars))
  pars
}

#' Build multi-slot phosphorus kinetics parameter vectors
#'
#' Builds phosphorus kinetics for multiple model "slots" and assembles the
#' results into K-length vectors (where `K = length(mods)`).
#'
#' Each slot is constructed by calling [build_P_kinetics()] with the same
#' parameter list `pparams`. This is useful when a simulation uses multiple
#' phosphorus modules that differ by model type but share the same raw parameter set.
#'
#' @param mods Character vector of model identifiers (e.g., `c("STA","RES")`).
#' @param registry Optional model registry (named list of builder functions).
#'   If `NULL` (default), the package's default registry is used.
#' @param pparams Named `list` of raw parameters. The same `pparams` is
#'   passed to every slot.
#' @param Dpy Numeric scalar. Time steps per year (default `365.25`).
#' @param DutyCycle Numeric scalar in \eqn{[0,1]}. Duty cycle multiplier.
#' @param derive_PModel Logical; reserved for future use.
#' @param default_PModel Integer scalar; reserved for future use.
#' @param ... Additional arguments passed to [build_P_kinetics()].
#'
#' @return A named `list` containing:
#' \describe{
#'   \item{K1, K2, K3}{Numeric vectors of length `K`.}
#'   \item{Chalf, Z_1, Z_2, Z_3, K2Coef}{Numeric vectors of length `K`.}
#'   \item{Kslots}{Integer scalar giving `K`.}
#'   \item{SeasonalFactor, Ytrans, Ysigma, Czero, PModel}{Global scalars.}
#'   \item{mods}{Character vector of model identifiers.}
#' }
#'
#' @examples
#' \donttest{
#' pparams <- list(
#'   # shared / STA
#'   C1000 = 1000, Cstar = 100, Ks_per_yr = 1,
#'   Z1 = 10, Z2 = 30, Z3 = 60,
#'   K2Coef1 = 0.1, Chalf = 50, SeasonalFactor = 1,
#'   # PSTA
#'   Ytrans = 1, Ysigma = 1, C1000_2 = 1000, ks_2 = 1, zh_2 = 10,
#'   # RES
#'   k_depth_penalty = 0.5
#' )
#'
#' out <- build_P_kin_slots(
#'   mods = c("STA", "PSTA", "RES"),
#'   pparams = pparams,
#'   Dpy = 365.25
#' )
#' out$K1
#' out$Z_1
#' }
#'
#' @export
build_P_kin_slots <- function(
    mods,
    registry = NULL ,
    pparams,
    Dpy = 365.25,
    DutyCycle = NULL,
    derive_PModel = TRUE,
    default_PModel = 1L,
    ...
) {

  if (is.null(registry)) registry <- .P_MODEL_BUILDERS_DEFAULT

  mods <- toupper(trimws(mods))
  if (any(!mods %in% names(registry))) {
    stop("All 'mods' must exist as names in 'registry'.")
  }

  DutyCycle <- pparams$DutyCycle


  # 1) Build each slot with your existing single-slot function
  slot_pars <- lapply(mods, function(m) {
    build_P_kinetics(
      mod_type = m,
      Dpy = Dpy,
      DutyCycle = DutyCycle,
      pparams = pparams,
      ...
    )
  })

  # 2) Assemble K-length vectors
  K <- length(slot_pars)

  K1 <- sapply(slot_pars, function(p) p$K1) # sapply(slot_pars, function(p) safe_num(p$K1, 0))
  K2 <- sapply(slot_pars, function(p) p$K2) # sapply(slot_pars, function(p) safe_num(p$K2, 0))
  K3 <- sapply(slot_pars, function(p) p$K3) # sapply(slot_pars, function(p) safe_num(p$K3, 0))

  Chalf  <- sapply(slot_pars, function(p) p$Chalf)  # sapply(slot_pars, function(p) safe_num(p$Chalf, 0))
  Z_1    <- sapply(slot_pars, function(p) p$Z1)     # sapply(slot_pars, function(p) safe_num(p$Z1, 0))
  Z_2    <- sapply(slot_pars, function(p) p$Z2)     # sapply(slot_pars, function(p) safe_num(p$Z2, 0))
  Z_3    <- sapply(slot_pars, function(p) p$Z3)     # sapply(slot_pars, function(p) safe_num(p$Z3, 0))
  K2Coef <- sapply(slot_pars, function(p) safe_num(p$K2Coef, 0))

  # 3) Global scalars: take first non-null/finite across slots
  SeasonalFactor <- first_non_null_finite(lapply(slot_pars, function(p) p$SeasonalFactor), 0)
  Ytrans         <- first_non_null_finite(lapply(slot_pars, function(p) p$Ytrans),         NA_real_)
  Ysigma         <- first_non_null_finite(lapply(slot_pars, function(p) p$Ysigma),         NA_real_)
  CZero          <- first_non_null_finite(lapply(slot_pars, function(p) p$Czero),          0)
  PModel         <- slot_pars[[1]]$PModel

  list(
    # K-length vectors
    K1 = K1, K2 = K2, K3 = K3,
    Chalf = Chalf,
    Z_1 = Z_1, Z_2 = Z_2, Z_3 = Z_3,
    K2Coef = K2Coef,
    # Global scalars
    Kslots = K,
    SeasonalFactor = SeasonalFactor,
    Ytrans = Ytrans,
    Ysigma = Ysigma,
    CZero = CZero,
    PModel = PModel,
    mods = mods
  )
}

#' Validate alignment of K-length phosphorus parameter vectors
#'
#' Internal helper that checks whether the K-length parameter vectors in
#' `params` are aligned with `params$Kslots`.
#'
#' @param params A named `list` containing `Kslots` and K-length vectors:
#'   `K1`, `K2`, `K3`, `Chalf`, `Z_1`, `Z_2`,
#'   `Z_3`, and `K2Coef`.
#'
#' @return Invisibly returns `TRUE` if validation passes; otherwise throws
#'   an error.
#'
#' @export
validate_P_paramsK <- function(params) {
  K <- params$Kslots
  if (length(params$K1) != K ||
      length(params$K2) != K ||
      length(params$K3) != K ||
      length(params$Chalf) != K ||
      length(params$Z_1) != K ||
      length(params$Z_2) != K ||
      length(params$Z_3) != K ||
      length(params$K2Coef) != K) {
    stop("K-length parameter vectors are not aligned.")
  }
  invisible(TRUE)
}
