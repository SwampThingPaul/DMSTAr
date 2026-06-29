#' Case-level network simulation utilities for DMSTAr
#'
#' A collection of helper functions that support **DMSTA-style CASE networks**
#' in \pkg{DMSTAr}. These utilities provide:
#'
#' \itemize{
#'   \item Parsing of DMSTA workbook "Networks" tables into routing definitions
#'   \item Topological ordering of CASEs for downstream execution
#'   \item Robust extraction of case- and cell-level output tables
#'   \item Stream-to-column mapping for routing discharge and load
#'   \item Vector lagging utilities for travel-time effects
#'   \item A high-level wrapper to run downstream CASE networks
#' }
#'
#' Together, these functions allow multiple DMSTA CASE simulations to be linked
#' into a strictly downstream network, emulating the behavior of the DMSTA VBA
#' workbook network framework while preserving DMSTAr's case- and cell-level APIs.
#'
#' @details
#' The network execution model follows standard DMSTA semantics:
#'
#' \enumerate{
#'   \item CASEs are executed in upstream-to-downstream order.
#'   \item Routed discharge and load from upstream CASE outputs are **added**
#'         to downstream CASE inflows.
#'   \item Concentrations are recomputed safely as \eqn{C = L / Q}.
#'   \item Routes may terminate at downstream CASEs or numeric outlet bins.
#' }
#'
#' These helpers are primarily used by `run_network_of_cases()`,
#' but many are also useful independently for diagnostics, testing, and
#' custom orchestration.
#'
#' @section Main user-facing functions:
#' \itemize{
#'   \item `build_routes_from_net_table()`
#'   \item `extract_df()`
#'   \item `run_network_of_cases()`
#' }
#'
#' @section Internal helpers:
#' \itemize{
#'   \item `extract_case_df`
#'   \item `stream_map_cols`
#'   \item `topo_order_cases`
#'   \item `lag_vec`
#' }
#'
#' @name dmsta_case_network
#' @keywords internal
NULL


#' Build routing table from a DMSTA-style network table
#'
#' Converts a DMSTA workbook "Networks" table into a normalized routing
#' `data.frame` suitable for downstream case-level network simulation.
#' Each non-blank destination in the network table becomes one route row.
#'
#' Destinations are interpreted as:
#' \itemize{
#'   \item A positive integer: route to an `"OUTLET"` bin.
#'   \item Otherwise: route to another `"CASE"` identified by string ID.
#'   \item Blank / whitespace: ignored (no route for that stream).
#' }
#'
#' The returned routes include `frac` (default 1) and `lag_days`
#' (default 0) to support future extensions such as flow splitting or travel-time.
#'
#' @param net_table A `data.frame` representing the DMSTA network table.
#'   Must contain `case_col` plus the stream destination columns
#'   `Bypass_to`, `Release1_to`, `Release2_to`, `Outflow_to`,
#'   and `Seepage_to`.
#' @param outlet_count Integer number of outlet bins. Numeric destinations must be
#'   within \eqn{1..outlet_count}.
#' @param case_col Character name of the column in `net_table` containing
#'   DMSTA CASE identifiers.
#'
#' @return A `data.frame` with one row per configured route, containing:
#' \describe{
#'   \item{from_case}{Upstream CASE ID (character).}
#'   \item{stream}{One of `"bypass"`, `"release1"`, `"release2"`,
#'                `"outflow"`, `"seepage"`.}
#'   \item{to_type}{`"CASE"` or `"OUTLET"`.}
#'   \item{to_id}{Downstream CASE ID (character) or outlet index (integer stored as character).}
#'   \item{frac}{Routing fraction (numeric), currently always 1.}
#'   \item{lag_days}{Integer travel-time lag in days, currently always 0.}
#' }
#'
#' @examples
#'
#' net <- data.frame(
#'   CaseName    = c("STA1_DW", "STA1W"),
#'   Bypass_to   = c("", ""),
#'   Release1_to = c("", ""),
#'   Release2_to = c("", ""),
#'   Outflow_to  = c("STA1W", "1"),
#'   Seepage_to  = c("", ""),
#'   stringsAsFactors = FALSE
#' )
#' build_routes_from_net_table(net, outlet_count = 1L)
#'
#' @rdname dmsta_case_network
#' @export
build_routes_from_net_table <- function(net_table,
                                        outlet_count = 5L,
                                        case_col = "CaseName") {
  req <- c(case_col, "Bypass_to", "Release1_to", "Release2_to", "Outflow_to", "Seepage_to")
  miss <- setdiff(req, names(net_table))
  if (length(miss)) stop("net_table missing columns: ", paste(miss, collapse = ", "))

  # DMSTA rule: blank means ignore stream.
  parse_dest <- function(x) {
    if (is.null(x) || length(x) == 0 || is.na(x) || trimws(as.character(x)) == "") return(NULL)
    x_chr <- trimws(as.character(x))
    # if numeric -> OUTLET bin
    if (grepl("^[0-9]+$", x_chr)) {
      k <- as.integer(x_chr)
      if (k < 1 || k > outlet_count) stop("OUTLET number out of range: ", k)
      return(list(to_type = "OUTLET", to_id = k))
    }
    # otherwise -> CASE
    list(to_type = "CASE", to_id = x_chr)
  }

  stream_cols <- c(
    bypass   = "Bypass_to",
    release1 = "Release1_to",
    release2 = "Release2_to",
    outflow  = "Outflow_to",
    seepage  = "Seepage_to"
  )

  out <- list()
  k <- 0L

  for (i in seq_len(nrow(net_table))) {
    from_case <- net_table[[case_col]][i]
    for (s in names(stream_cols)) {
      dest <- parse_dest(net_table[[stream_cols[[s]]]][i])
      if (!is.null(dest)) {
        k <- k + 1L
        out[[k]] <- data.frame(
          from_case = from_case,
          stream    = s,
          to_type   = dest$to_type,
          to_id     = dest$to_id,
          frac      = 1,
          lag_days  = 0L,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, out)
}

#' Extract a case- or cell-level daily output table from a DMSTAr result
#'
#' Flexible extractor for pulling either the case-level output (`level = "case"`)
#' or a specific cell-level output (`level = "cell"`) from a list returned by
#' `dmsta_flowP_case()`.
#'
#' @param res A list-like object, typically returned by `dmsta_flowP_case()`.
#' @param level Character, one of `"case"` or `"cell"`.
#' @param cell_index Integer index of the cell to extract when `level = "cell"`.
#'
#' @return A `data.frame` containing the requested daily time series.
#'
#' @details
#' Current DMSTAr documentation shows `dmsta_flowP_case()` returning a list
#' with `out$results$case` and `out$results$cells[[i]]` (or similarly named
#' cell outputs).
#' This function supports that nested style and also checks common legacy fields
#' such as `res$case` / `res$case_out`.
#'
#' @examples
#' \dontrun{
#' out <- dmsta_flowP_case(series, cells, return_cell_series = TRUE)
#' df_case <- extract_df(out, level = "case")
#' df_cell1 <- extract_df(out, level = "cell", cell_index = 1)
#' }
#' @rdname dmsta_case_network
#' @export
extract_df <- function(res,
                       level = c("case", "cell"),
                       cell_index = 1L) {
  level <- match.arg(level)

  if (!is.list(res)) {
    stop("`res` must be a list (dmsta_flowP_case result).")
  }

  # DMSTAr-style nesting: res$results is a list, containing $case and $cell
  if (!is.null(res$results) && is.list(res$results)) {

    if (level == "case") {
      if (!is.null(res$results$case) && is.data.frame(res$results$case)) {
        return(res$results$case)
      }
      # Some variants may use case_out
      if (!is.null(res$results$case_out) && is.data.frame(res$results$case_out)) {
        return(res$results$case_out)
      }
    }

    if (level == "cell") {
      if (!is.null(res$results$cell) && is.list(res$results$cell)) {
        cell_index <- as.integer(cell_index)
        if (cell_index < 1L || cell_index > length(res$results$cell)) {
          stop("cell_index out of range: 1..", length(res$results$cell))
        }
        if (is.data.frame(res$results$cell[[cell_index]])) {
          return(res$results$cell[[cell_index]])
        }
      }
    }
  }

  # Legacy fallbacks (keep these to be robust across versions)
  if (level == "case") {
    if (!is.null(res$case) && is.data.frame(res$case)) return(res$case)
    if (!is.null(res$case_out) && is.data.frame(res$case_out)) return(res$case_out)
  }

  # Last resort: find first data.frame at top level
  dfs <- vapply(res, is.data.frame, logical(1))
  if (any(dfs)) return(res[[which(dfs)[1]]])

  stop("Could not locate requested output data.frame in dmsta_flowP_case() result.")
}

#' Map a routed stream name to DMSTAr output column names
#'
#' Resolves a stream identifier (e.g., `"bypass"`, `"release1"`,
#' `"outflow"`) to the corresponding discharge (`Q_*`) and load
#' (`L_*`) columns in a DMSTAr case output table.
#'
#' @param stream Character stream name. Supported values:
#'   `"bypass"`, `"release1"`, `"release2"`, `"outflow"`, `"seepage"`.
#' @param outflow_def Character definition for `"outflow"` mapping:
#'   `"treated"` maps to treated outflow columns; `"total"` maps to total outflow columns.
#'
#' @return A named character vector of length 2 with names `Q` and `L`,
#'   giving the discharge and load column names.
#'
#' @details
#' DMSTAr case outputs include standardized `Q_*` and `L_*` fields for
#' bypass, releases, seepage discharge, and outflows.
#'
#' @keywords internal
#' @rdname dmsta_case_network
stream_map_cols <- function(stream, outflow_def = c("treated", "total")) {
  outflow_def <- match.arg(outflow_def)
  stream <- tolower(stream)

  if (stream == "outflow") {
    if (outflow_def == "treated") {
      return(c(Q="Q_out_treated", L="L_out_treated"))
    } else {
      return(c(Q="Q_out_total", L="L_out_total"))
    }
  }

  switch(
    stream,
    bypass   = c(Q="Q_out_bypass",         L="L_out_bypass"),
    release1 = c(Q="Q_out_release1",       L="L_out_release1"),
    release2 = c(Q="Q_out_release2",       L="L_out_release2"),
    seepage  = c(Q="Q_out_seep_discharge", L="L_out_seep_discharge"),
    stop("Unknown stream: ", stream)
  )
}

#' Compute a topological execution order for downstream-only case networks
#'
#' Given a routing table and a set of case names, computes a topological order
#' (upstream to downstream) under the assumption that the case network is a DAG
#' (i.e., contains no directed cycles).
#'
#' @param routes A `data.frame` of routes such as produced by
#'   `build_routes_from_net_table()`. Must include `from_case`,
#'   `to_type`, and `to_id`.
#' @param case_names Character vector of all CASE IDs participating in the network.
#'
#' @return Character vector of CASE IDs ordered such that upstream cases appear
#'   before any downstream cases they route to.
#'
#' @details
#' This uses Kahn's algorithm over CASE-to-CASE edges (routes with `to_type == "CASE"`).
#' If not all nodes can be ordered, the function stops, indicating a directed cycle.
#'
#' @importFrom stats setNames
#' @importFrom utils head tail
#'
#' @keywords internal
#' @rdname dmsta_case_network
topo_order_cases <- function(routes, case_names) {
  edges <- routes[routes$to_type == "CASE", c("from_case", "to_id"), drop = FALSE]
  names(edges) <- c("from", "to")
  edges <- unique(edges)

  incoming <- setNames(rep(0L, length(case_names)), case_names)
  adj <- setNames(vector("list", length(case_names)), case_names)

  for (i in seq_len(nrow(edges))) {
    a <- edges$from[i]; b <- edges$to[i]
    adj[[a]] <- unique(c(adj[[a]], b))
    incoming[[b]] <- incoming[[b]] + 1L
  }

  q <- names(incoming[incoming == 0L])
  out <- character(0)

  while (length(q)) {
    n <- q[1]; q <- q[-1]
    out <- c(out, n)
    for (m in adj[[n]]) {
      incoming[[m]] <- incoming[[m]] - 1L
      if (incoming[[m]] == 0L) q <- c(q, m)
    }
  }

  if (length(out) != length(case_names)) {
    stop("Network contains a cycle; strictly downstream requirement violated.")
  }
  out
}

#' Shift a numeric vector by an integer number of days
#'
#' Applies a simple discrete lag/lead to a numeric vector, padding with zeros.
#' Positive `lag_days` shifts values later in time (prepends zeros);
#' negative values shift earlier in time (appends zeros).
#'
#' @param x Numeric vector.
#' @param lag_days Integer number of days to shift. Positive shifts forward in time;
#'   negative shifts backward. Zero returns `x` unchanged.
#'
#' @return Numeric vector of the same length as `x`.
#'
#' @importFrom utils head tail
#'
#' @examples
#' \dontrun{
#' lag_vec(1:5, 2)   # 0 0 1 2 3
#' lag_vec(1:5, -2)  # 3 4 5 0 0
#' }
#' @keywords internal
#' @rdname dmsta_case_network
lag_vec <- function(x, lag_days) {
  lag_days <- as.integer(lag_days)
  if (lag_days == 0L) return(x)
  n <- length(x)
  if (lag_days > 0L) c(rep(0, lag_days), head(x, n - lag_days)) else {
    k <- abs(lag_days)
    c(tail(x, n - k), rep(0, k))
  }
}


#' Run a DMSTA-style downstream network simulation across CASEs
#'
#' Executes a set of CASE simulations in topological (upstream-to-downstream)
#' order and routes configured outflows/loads from upstream CASE outputs into
#' downstream CASE inflow forcing, emulating DMSTA workbook "Networks" behavior.
#'
#' For each CASE:
#' \enumerate{
#'   \item Base inflow discharge and concentration (`Qi`, `Ci`) are sanitized.
#'   \item Routed inflow discharge/load from upstream CASEs are added.
#'   \item A safe routed inflow concentration is computed (`Cin = Lin / Qin`).
#'   \item `dmsta_flowP_case()` is executed for that CASE.
#'   \item Specified stream outputs are routed to downstream CASEs or outlet bins.
#' }
#'
#' @param cases Named list of CASE definitions. Names must be DMSTA CASE IDs.
#'   Each element must contain at least:
#'   \itemize{
#'     \item `series_base`: a `data.frame` with columns `Date`, `Qi`, `Ci`,
#'           `Rain`, `Et`, `Zcontrol`
#'     \item `cells`: a list of cell definitions for `dmsta_flowP_case()`
#'   }
#' @param net_table Optional DMSTA-style network table. If provided and `routes` is `NULL`,
#'   routing is built with `build_routes_from_net_table()`.
#' @param routes Optional normalized routing table. If provided, `net_table` is ignored.
#' @param outlet_count Integer number of outlet bins (used when parsing numeric destinations).
#' @param verbose Logical; if `TRUE`, prints execution order and per-CASE summaries.
#' @param check_route Logical; if `TRUE`, prints per-CASE stream discharge sums for debugging.
#' @param ... Additional arguments passed through to `dmsta_flowP_case()`.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{order}{Character vector of CASE execution order.}
#'   \item{case_results}{Named list of raw `dmsta_flowP_case()` results for each CASE.}
#'   \item{routed_in}{Named list of routed inflow time series per CASE (`Q`, `L`).}
#'   \item{outlets}{List of outlet-bin time series (`Q`, `L`) of length `outlet_count`.}
#'   \item{outlet_summary}{`data.frame` summarizing total outlet discharge/load and FWC.}
#'   \item{ledger}{`data.frame` of per-route totals (`total_Q`, `total_L`) for auditing.}
#' }
#'
#' @details
#' This wrapper relies on the standardized case output columns produced by
#' `dmsta_flowP_case()` (e.g., `Q_out_bypass`, `Q_out_treated`,
#' `Q_out_total`, and corresponding `L_*` fields).
#' It also expects that all CASEs share the same `Date` sequence (same length and alignment).
#'
#' The outlet summary uses `fw()` to compute flow-weighted concentration (FWC) from total
#' load and discharge; ensure `fw()` is available in your namespace.
#'
#' @examples
#' \dontrun{
#' # 1) Build cases
#' cases <- list(
#'  STA1_DW = list(
#'   series_base = STA1_DW_input,   # data.frame with Date, Qi, Ci, Rain, Et, Zcontrol
#'   cells       = STA1_DW_cell     # list of dmsta_make_cell(...) objects
#' ),
#' STA1W = list(
#'   series_base = STA1W_input,
#'   cells       = STA1W_cell
#' )
#' )
#'
#' # 2) Build/parse routes
#' net <- data.frame(
#'   CaseName    = c("STA1_DW", "STA1W"),
#'   Bypass_to   = c("", ""),
#'   Release1_to = c("", ""),
#'   Release2_to = c("", ""),
#'   Outflow_to  = c("STA1W", "1"),
#'   Seepage_to  = c("", ""),
#'   stringsAsFactors = FALSE
#' )
#' routes <- build_routes_from_net_table(net_table, outlet_count = 1L)
#'
#' # 3) Run network
#' out <- run_network_of_cases(
#'   cases = cases,
#'   routes = routes,
#'   Nsteps = 4L,
#'   return_cell_series = TRUE
#' )
#'
#' out$outlet_summary
#' head(out$ledger)
#' }
#' @seealso `dmsta_flowP_case()`
#' @rdname dmsta_case_network
#' @export

run_network_of_cases <- function(
    cases,
    net_table = NULL,
    routes = NULL,
    outlet_count = 5L,
    verbose = TRUE,
    check_route = FALSE,
    Nsteps_case = NULL,
    ...
) {

  # Resolve case names
  case_names <- names(cases)
  if (is.null(case_names) || any(case_names == "")) {
    stop("cases must be a named list; names are DMSTA CASE IDs.")
  }

  # Build routing table if needed
  if (is.null(routes)) {
    if (is.null(net_table)) stop("Provide either net_table or routes.")
    routes <- build_routes_from_net_table(net_table, outlet_count = outlet_count)
  }

  # Validate route endpoints
  bad_from <- setdiff(unique(routes$from_case), case_names)
  if (length(bad_from)) stop("Unknown from_case in routes: ", paste(bad_from, collapse = ", "))

  bad_to <- setdiff(unique(routes$to_id[routes$to_type == "CASE"]), case_names)
  if (length(bad_to)) stop("Unknown downstream CASE in routes: ", paste(bad_to, collapse = ", "))

  # Determine execution order (topological)
  exec_order <- topo_order_cases(routes, case_names)
  if (verbose) {
    message("Case execution order: ", paste(exec_order, collapse = " -> "))
  }

  # Shared date sequence
  net_dates <- cases[[exec_order[1]]]$series_base$Date
  nday <- length(net_dates)

  # Routed inflows (Q, L) per case
  routed_in <- lapply(cases, function(cd) list(
    Q = numeric(nday),
    L = numeric(nday)
  ))

  # Outlet bins
  outlets <- lapply(seq_len(outlet_count), function(i) {
    list(Q = numeric(nday), L = numeric(nday))
  })
  names(outlets) <- paste0("Outlet", seq_len(outlet_count))

  # Output containers
  case_results <- vector("list", length(case_names))
  names(case_results) <- case_names

  ledger <- data.frame(
    from_case = character(0),
    stream = character(0),
    to_type = character(0),
    to_id = character(0),
    frac = numeric(0),
    lag_days = integer(0),
    total_Q = numeric(0),
    total_L = numeric(0),
    stringsAsFactors = FALSE
  )


  # NEW BLOCK: Extract integrator controls from
  dots <- list(...)

  # Defaults
  Qmethod <- "RK4"
  Pmethod <- "RK4"
  interp_option <- 2L
  integrator_fun <- NULL

  if (!is.null(dots$Qmethod)) {
    Qmethod <- dots$Qmethod
    dots$Qmethod <- NULL
  }
  if (!is.null(dots$Pmethod)) {
    Pmethod <- dots$Pmethod
    dots$Pmethod <- NULL
  }
  if (!is.null(dots$interp_option)) {
    interp_option <- as.integer(dots$interp_option)
    dots$interp_option <- NULL
  }
  if (!is.null(dots$integrator_fun)) {
    integrator_fun <- dots$integrator_fun
    dots$integrator_fun <- NULL
  }

  Qmethod <- match.arg(Qmethod, c("RK4", "Euler", "RKF45", "custom"))
  Pmethod <- match.arg(Pmethod, c("RK4", "Euler"))

  # Resolve Nsteps per case
  nsteps_default <- if (!is.null(dots$Nsteps)) {
    as.integer(dots$Nsteps)
  } else {
    4L
  }
  dots$Nsteps <- NULL

  resolve_nsteps <- function(case_id) {
    if (is.null(Nsteps_case)) return(nsteps_default)
    if (length(Nsteps_case) == 1L && is.null(names(Nsteps_case))) {
      return(as.integer(Nsteps_case))
    }
    if (!is.null(names(Nsteps_case)) && case_id %in% names(Nsteps_case)) {
      return(as.integer(Nsteps_case[[case_id]]))
    }
    nsteps_default
  }

  # Execute cases in downstream order
  for (cn in exec_order) {
    cd <- cases[[cn]]
    base <- cd$series_base

    needed <- c("Date", "Qi", "Ci", "Rain", "Et", "Zcontrol")
    if (!all(needed %in% names(base))) {
      stop("Case ", cn, " series_base must include: ", paste(needed, collapse = ", "))
    }
    if (length(base$Date) != nday) {
      stop("Case ", cn, " Date length differs from network calendar.")
    }

    # Combine base + routed inflows (DMSTA semantics)
    Qb <- base$Qi
    Cb <- base$Ci
    Qb[!is.finite(Qb)] <- 0
    Cb[!is.finite(Cb)] <- 0
    Lb <- Qb * Cb

    # sanitize routed inflows too (defensive)
    rq <- routed_in[[cn]]$Q
    rl <- routed_in[[cn]]$L
    rq[!is.finite(rq)] <- 0
    rl[!is.finite(rl)] <- 0

    # combine (DMSTA network semantics: add upstream to base)
    Qin <- Qb + rq
    Lin <- Lb + rl
    # compute Cin safely
    Cin <- ifelse(Qin > 0, Lin / Qin, 0)
    Cin[!is.finite(Cin)] <- 0

    # NOW build series_run using the sanitized values
    series_run <- base
    series_run$Qi <- Qin
    series_run$Ci <- Cin

    nsteps_i <- resolve_nsteps(cn)
    ## To check the process
    if (verbose) {
      message("Running CASE: ", cn, " (Nsteps=", nsteps_i, ")")
    }
    if (verbose) {
      message(sprintf(
        "CASE %s: baseQ sum=%.4f anyNA=%s | routedQ sum=%.4f anyNA=%s | routedL sum=%.4f anyNA=%s",
        cn,
        sum(base$Qi, na.rm=TRUE), anyNA(base$Qi),
        sum(routed_in[[cn]]$Q, na.rm=TRUE), anyNA(routed_in[[cn]]$Q),
        sum(routed_in[[cn]]$L, na.rm=TRUE), anyNA(routed_in[[cn]]$L)
      ))
    }

    # Run the case simulation
    res <- do.call(
      dmsta_flowP_case,
      c(
        list(
          series = series_run,
          cells = cd$cells,
          Nsteps = nsteps_i,
          Qmethod = Qmethod,
          Pmethod = Pmethod,
          integrator_fun = integrator_fun,
          interp_option = interp_option
        ),
        dots
      )
    )

    case_results[[cn]] <- res
    df <- extract_df(res)

    # Adds summary of routing printed in console
    if(check_route){
      message(cn," sums:",
              " bypass=", sum(df$Q_out_bypass, na.rm=TRUE),
              " treated=", sum(df$Q_out_treated, na.rm=TRUE),
              " rel1=", sum(df$Q_out_release1, na.rm=TRUE),
              " rel2=", sum(df$Q_out_release2, na.rm=TRUE),
              " seep=", sum(df$Q_out_seep_discharge, na.rm=TRUE),
              " total=", sum(df$Q_out_total, na.rm=TRUE))
    }

    # Route outputs downstream
    rsub <- routes[routes$from_case == cn, , drop = FALSE]
    if (nrow(rsub)) {
      for (ri in seq_len(nrow(rsub))) {
        r <- rsub[ri, ]
        cols <- stream_map_cols(r$stream, outflow_def = "treated")

        if (!(cols["Q"] %in% names(df) && cols["L"] %in% names(df))) {
          stop("Case ", cn, " output missing: ", cols["Q"], " / ", cols["L"])
        }

        Qs <- df[[cols["Q"]]] * r$frac
        Ls <- df[[cols["L"]]] * r$frac
        Qs[!is.finite(Qs)] <- 0
        Ls[!is.finite(Ls)] <- 0

        # optional lag
        if (!is.null(r$lag_days) && r$lag_days != 0L) {
          Qs <- lag_vec(Qs, r$lag_days)
          Ls <- lag_vec(Ls, r$lag_days)
        }

        if (r$to_type == "CASE") {
          to_case <- as.character(r$to_id)
          routed_in[[to_case]]$Q <- routed_in[[to_case]]$Q + Qs
          routed_in[[to_case]]$L <- routed_in[[to_case]]$L + Ls
        } else {
          k <- as.integer(r$to_id)
          outlets[[k]]$Q <- outlets[[k]]$Q + Qs
          outlets[[k]]$L <- outlets[[k]]$L + Ls
        }

        ledger <- rbind(
          ledger,
          data.frame(
            from_case = cn,
            stream = r$stream,
            to_type = r$to_type,
            to_id = as.character(r$to_id),
            frac = r$frac,
            lag_days = as.integer(r$lag_days),
            total_Q = sum(Qs, na.rm = TRUE),
            total_L = sum(Ls, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  # Outlet summary
  outlet_summary <- do.call(rbind, lapply(seq_len(outlet_count), function(k) {
    Q <- outlets[[k]]$Q
    L <- outlets[[k]]$L
    data.frame(
      outlet = paste0("Outlet ", k),
      total_Q = sum(Q, na.rm = TRUE),
      total_L = sum(L, na.rm = TRUE),
      FWC = fw(sum(L, na.rm = TRUE), sum(Q, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))

  out <- list(
    order = exec_order,
    case_results = case_results,
    routed_in = routed_in,
    outlets = outlets,
    outlet_summary = outlet_summary,
    ledger = ledger
  )

  attr(out, "Qmethod") <- Qmethod
  attr(out, "Pmethod") <- Pmethod
  attr(out, "interp_option") <- interp_option

  out
}


