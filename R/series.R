#' Example DMSTA daily forcing time series
#'
#' A example time series of daily hydrology and phosphorus forcing inputs
#' suitable for demonstrating DMSTAr case and cell simulations.
#'
#' @format A data frame with daily records and the following columns:
#' \describe{
#'   \item{Date}{Date for the record. Stored as a character vector in `"YYYY-MM-DD"`
#'     format in the raw object shown; you may convert to `Date` with `as.Date()`.}
#'   \item{Flow}{Daily inflow/flow rate (cfs).}
#'   \item{Conc}{Daily inflow concentration (e.g., ug/L).}
#'   \item{Rainfall}{Daily rainfall rate (e.g., in/day).}
#'   \item{ET}{Daily evapotranspiration rate (e.g., in/day).}
#' }
#'
#' @details
#' This dataset is intended as a lightweight example for testing and
#' documentation. For use with functions that expect columns named
#' `Qi`, `Ci`, and `Rain`, you may need to rename columns
#' (e.g., `Flow -> Qi`, `Conc -> Ci`, `Rainfall -> Rain`).
#'
#' @source
#' Source workbook: `PROJECT_SFWMD_EC_01MAR2012_NET_EAA_STA1E.xls`
#'
#' sheet `Series_Input`
#'
#' @examples
#' data(series)
#'
#' # Convert Date column (if needed)
#' series$Date <- as.Date(series$Date)
#'
#' # If your simulation expects Qi/Ci/Rain column names:
#' sim_series <- within(series, {
#'   Qi   <- Flow
#'   Ci   <- Conc
#'   Rain <- Rainfall
#' })
#' sim_series$Flow <- NULL
#' sim_series$Conc <- NULL
#' sim_series$Rainfall <- NULL
#'
#' head(sim_series)
"series"
