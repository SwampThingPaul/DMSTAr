#' DMSTA tank geometry and state initialization (internal)
#'
#' Internal helpers for constructing tank partitions and initializing
#' DMSTA phosphorus state vectors.
#'
#' @name internal_dmsta_tanks
#' @keywords internal
NULL

#' Build DMSTA tank partitioning for a cell
#'
#' Internal helper that partitions a single cell area into a series of
#' conceptual "tanks" used by DMSTA. The number of tanks is derived from
#' `ttankS`; if `ttankS` is fractional, the final tank receives
#' the fractional area share and all preceding tanks receive equal shares.
#'
#' The function also returns per-tank area fractions and cumulative fractions
#' (useful for mapping depth/area relationships). Optionally, the last cumulative
#' fraction can be "snapped" to exactly 1.0 for numerical stability.
#'
#' @param A_cell Numeric scalar > 0. Cell area (units consistent with the rest
#'   of the DMSTA implementation).
#' @param ttankS Numeric scalar > 0. Effective number of tanks. May be fractional.
#' @param snap_last Logical; if `TRUE`, force the last cumulative fraction
#'   to exactly `1.0`. Default is `TRUE`.
#'
#' @return A named `list` with elements:
#' \describe{
#'   \item{Ntanks}{Integer number of tanks.}
#'   \item{A_Tank}{Numeric vector of tank areas, length `Ntanks`.}
#'   \item{F_Tank}{Numeric vector of tank area fractions (`A_Tank / A_cell`).}
#'   \item{Fcum}{Numeric vector of cumulative area fractions (`cumsum(F_Tank)`).}
#' }
#'
#' @seealso [dmsta_p_init_state()] for initializing state vectors
#'   compatible with the returned `Ntanks`.
#'
#' @export
dmsta_build_tanks <- function(A_cell, ttankS,snap_last = TRUE) {
  # ttankS can be fractional; last tank gets the fraction, all others equal.
  if (A_cell <= 0) {
    return(list(Ntanks = 1L, A_Tank = 0, F_Tank = 1, Fcum = 1))
  }
  if (ttankS <1) ttankS <- 1

  Ntanks <- max(1L, floor(ttankS))
  frac <- ttankS - Ntanks
  if (frac > 0) Ntanks <- Ntanks + 1L

  A_Tank <- rep(A_cell / max(1.0, ttankS), Ntanks)
  if (frac > 0) A_Tank[Ntanks] <- A_cell * frac / ttankS

  F_Tank <- A_Tank / A_cell
  Fcum <- cumsum(F_Tank)
  # snap last cumulative fraction to exactly 1 for stability
  if(snap_last==TRUE){Fcum[Ntanks] <- 1.0}
  list(Ntanks = Ntanks, A_Tank = A_Tank, F_Tank = F_Tank, Fcum = Fcum)
}

#' Initialize DMSTA phosphorus state vectors for tanks
#'
#' Internal helper that initializes per-tank state vectors for the DMSTA
#' phosphorus module given tank geometry and initial conditions.
#'
#' For each tank `i`, the function computes:
#' \itemize{
#'   \item `M[i] = C_init_ppb * A_Tank[i] * Z_init_m`
#'   \item `S[i] = Y_init_mgm2 * A_Tank[i]`
#' }
#'
#' where `A_Tank[i]` is the area of tank `i`. Units are assumed to be
#' consistent with the DMSTA implementation (e.g., depth in meters).
#'
#' @param tanks A list as returned by [dmsta_build_tanks()] containing
#'   `Ntanks` and `A_Tank`.
#' @param Z_init_m Numeric scalar. Initial water column depth (meters).
#' @param C_init_ppb Numeric scalar. Initial concentration (ppb).
#' @param Y_init_mgm2 Numeric scalar. Initial areal mass/loading (mg/m^2).
#'
#' @return A named `list` with elements:
#' \describe{
#'   \item{M}{Numeric vector (length `tanks$Ntanks`) of initialized water-column masses.}
#'   \item{S}{Numeric vector (length `tanks$Ntanks`) of initialized sediment/areal stores.}
#' }
#'
#' @seealso [dmsta_build_tanks()] to generate tank geometry.
#'
#' @export
dmsta_p_init_state <- function(tanks, Z_init_m, C_init_ppb, Y_init_mgm2) {
  Nt <- tanks$Ntanks
  M <- numeric(Nt)
  S <- numeric(Nt)
  for (i in seq_len(Nt)) {
    A <- tanks$A_Tank[i]
    # VBA: M_tank = C_init * A_tank * Z  (ppb * hm3 = kg)
    M[i] <- C_init_ppb * A * Z_init_m
    # VBA: S_tank = Y_init * A_tank  (mg/m2 * km2 = kg)
    S[i] <- Y_init_mgm2 * A
  }
  list(M = M, S = S)
}
