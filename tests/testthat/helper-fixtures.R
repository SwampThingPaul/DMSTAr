
make_series_base <- function(n = 10) {
  data(series, package = "DMSTAr")  # example dataset exists in package [1](https://evergladesfoundation-my.sharepoint.com/personal/pjulian_evergladesfoundation_org/Documents/Documents/R/win-library/4.1/DMSTAr/html/series.html?web=1)
  s <- series[seq_len(n), ]

  # reshape per docs: Flow -> Qi, Conc -> Ci, Rainfall -> Rain, ET -> Et
  # keep Date; set Zcontrol to 0 by default
  out <- data.frame(
    Date     = as.Date(s$Date),
    Qi       = as.numeric(s$Flow),
    Ci       = as.numeric(s$Conc),
    Rain     = as.numeric(s$Rainfall),
    Et       = as.numeric(s$ET),
    Zcontrol = rep(0, n)
  )
  out
}

# Minimal parameter list: use dmstar_default_params and override what's needed.
make_min_params <- function(...) {
  p <- dmstar_default_params()  # shipped by package [7](https://evergladesfoundation-my.sharepoint.com/personal/pjulian_evergladesfoundation_org/Documents/Documents/R/win-library/4.1/DMSTAr/html/dmstar_default_params.html?web=1)
  modifyList(p, list(...))
}

# Create a "node" cell (A_cell <= 0 triggers IsaNode derivation per docs) [7](https://evergladesfoundation-my.sharepoint.com/personal/pjulian_evergladesfoundation_org/Documents/Documents/R/win-library/4.1/DMSTAr/html/dmstar_default_params.html?web=1)
make_node_cell <- function(label = "NODE") {
  params <- make_min_params(
    A_cell     = 0,     # node-like; avoids hydrology complexities
    DutyCycle  = 1,
    Zinit      = 0,
    C_init_ppb = 0,
    Y_init_mgm2= 0
  )

  cell <- dmsta_make_cell(
    label = label,
    params = params,
    ttankS = 1,
    DownCell = 0L,
    Qin_Frac = 1,
    RecycleIndex = 1L
  )
  dmsta_validate_cells(list(cell))  # normalize/validate [6](https://evergladesfoundation-my.sharepoint.com/personal/pjulian_evergladesfoundation_org/Documents/Documents/R/win-library/4.1/DMSTAr/html/dmsta_validate_cells.html?web=1)
}

make_case <- function(name, n = 10) {
  list(
    series_base = make_series_base(n),
    cells = make_node_cell(label = paste0(name, "_NODE"))
  )
}

make_cases_2 <- function(n = 10) {
  list(
    UP = make_case("UP", n),
    DN = make_case("DN", n)
  )
}

# Simple net_table fixture matching your parser requirements
make_net_table_simple <- function() {
  data.frame(
    CaseName     = c("UP", "DN"),
    Bypass_to    = c("", ""),
    Release1_to  = c("", ""),
    Release2_to  = c("", ""),
    Outflow_to   = c("DN", "1"),  # UP -> DN, DN -> Outlet 1
    Seepage_to   = c("", ""),
    stringsAsFactors = FALSE
  )
}
