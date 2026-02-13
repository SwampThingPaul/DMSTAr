#' Build parameter lists for supported DMSTAr models
#'
#' These functions convert a raw parameter list (`pparams`) into a
#' standardized parameter list used by DMSTAr model implementations.
#'
#' Builders perform unit conversions and time-step scaling:
#' \itemize{
#'   \item rate parameters per year are converted to per-timestep using `Dpy`
#'   \item depth parameters in centimeters are converted to meters
#'   \item duty cycle is applied to rate parameters where applicable
#' }
#'
#' @param Dpy Number of time steps per year (e.g., 365 for daily).
#' @param DutyCycle Fraction of time active in each step (0–1).
#' @param pparams A named `list` of raw parameters for the selected model.
#' @param ... Reserved for future extensions; ignored.
#'
#' @return A named `list` of standardized model parameters.
#'
#' @name model_parameter_builders
NULL

#' Create an empty parameter list for DMSTAr model builders
#'
#' Returns a named list initialized with default values used by DMSTAr
#' parameter-builder functions (e.g., `build_STA()`, `build_PSTA()`,
#' `build_RES()`).
#'
#' This helper is primarily intended for internal use when constructing model
#' parameter objects with a consistent set of fields.
#'
#' @return A named `list` of parameter values. Numeric parameters are set
#'   to `0`; transformation fields are set to `NA`.
#'
#' @keywords internal
#' @noRd
empty_P_pars <- function() {
  list(
    C1000 = 0,
    Cstar = 0,
    Ks = 0,
    Z1 = 0,
    Z2 = 0,
    Z3 = 0,
    Kc = 0,
    Chalf = 0,
    SeasonalFactor = 0,
    Ytrans = NA,
    Ysigma = NA,
    Czero = 0,
    ks_2 = 0,
    zh_2 = 0
  )
}


#' @describeIn model_parameter_builders Build parameters for the STA model.
#' @export
build_STA <- function(Dpy, DutyCycle, pparams, ... ) {
  # STA
  c1000_1 <-  pparams$C1000
  cstar_1 <- pparams$Cstar
  ks_1 <- pparams$Ks_per_yr
  Z1 <- pparams$Z1
  Z2 <- pparams$Z2
  Z3 <- pparams$Z3
  k2res <- pparams$K2Coef1
  kc2 <- pparams$Chalf
  SeasonalFactor <- pparams$SeasonalFactor

  pars <- empty_P_pars()

  pars$C1000 <- c1000_1
  pars$Cstar <- cstar_1
  pars$Ks <- ks_1 / Dpy * DutyCycle
  pars$Z1 <- Z1 / 100 # cm -> m
  pars$Z2 <- Z2 / 100 # cm -> m
  pars$Z3 <- Z3 / 100 # cm -> m
  pars$K2Coef <- k2res / Dpy
  pars$Chalf <- kc2
  pars$SeasonalFactor <- SeasonalFactor

  pars
}

#' @describeIn model_parameter_builders Build parameters for the PSTA model.
#' @export
build_PSTA <- function(Dpy, DutyCycle, pparams, ... ) {
  # PSTA  part of NEWS
  y00 <- pparams$Ytrans
  ys00 <- pparams$Ysigma
  c1000_2 <- pparams$C1000_2
  cstar_2 <- pparams$SeasonalFactor
  ks_2 <- pparams$ks_2
  zh_2 <- pparams$zh_2

  pars <- empty_P_pars()

  pars$Ytrans <- y00
  pars$Ysigma <- ys00
  pars$Czero  <- 0

  if (y00 > 0) {
    pars$C1000 <- c1000_2
    pars$Cstar <- cstar_2
    pars$Ks <- ks_2 / Dpy * DutyCycle
    pars$Z1 <- zh_2 / 100 # cm -> m
  }

  pars
}

#' @describeIn model_parameter_builders Build parameters for the RES model.
#' @export
build_RES <- function(Dpy, DutyCycle, pparams, ...) {
  # RES
  k_depth_penalty <- pparams$k_depth_penalty
  ks_1 <- pparams$Ks_per_yr
  c1000_1 <- pparams$C1000
  cstar_1 <- pparams$Cstar
  Z1_1 <- pparams$Z1
  kc2 <- pparams$Chalf

  pars <- empty_P_pars()

  Ks_res <- k_depth_penalty / Dpy * DutyCycle
  Ks_sta <- ks_1 / Dpy * DutyCycle

  pars$C1000 <- c1000_1
  pars$Cstar <- cstar_1
  pars$Ks <- min(Ks_res, Ks_sta)
  pars$Z1 <- Z1_1 / 100 # cm -> m
  pars$Chalf <- kc2

  pars
}


#' Default phosphorus model builders (internal)
#'
#' Named list of built-in phosphorus model builder functions used as the default
#' registry for DMSTAr when no runtime registry is supplied.
#'
#' @format A named list of functions with elements `STA`, `PSTA`, and `RES`.
#' @keywords internal
#' @noRd
.P_MODEL_BUILDERS_DEFAULT <- list(
  STA  = build_STA,
  PSTA = build_PSTA,
  RES  = build_RES
)


## next - build a registry adder function for future model types.
