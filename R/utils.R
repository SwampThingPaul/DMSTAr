
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
  m3_per_ft3  <- 0.3048^3
  m3_per_day  <- m3_per_ft3 * 86400
  hm3_per_day <- m3_per_day / 1e6
  x * hm3_per_day
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


