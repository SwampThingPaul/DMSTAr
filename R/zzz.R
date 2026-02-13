# Private package environment
# .dmstar_env <- new.env(parent = emptyenv())

# .onLoad <- function(libname, pkgname) {
#
#   # Ensure registry exists
#   .dmstar_env$P_MODEL_BUILDERS <- new.env(parent = emptyenv())
#
#   # Register built-ins
#   register_P_model("STA",  build_STA)
#   register_P_model("PSTA", build_PSTA)
#   register_P_model("RES",  build_RES)
#
#   # --- Compatibility layer: P_MODEL_BUILDERS behaves like a list ---
#   ns <- asNamespace(pkgname)
#
#   # If active binding already exists (common in dev reload), do nothing.
#   if (exists("P_MODEL_BUILDERS", envir = ns, inherits = FALSE) &&
#       bindingIsActive("P_MODEL_BUILDERS", ns)) {
#     return(invisible())
#   }
#
#   # If a non-active binding exists, try to remove it; if locked, fail loudly
#   if (exists("P_MODEL_BUILDERS", envir = ns, inherits = FALSE)) {
#     if (bindingIsLocked("P_MODEL_BUILDERS", ns)) {
#       stop(
#         "Cannot install active binding 'P_MODEL_BUILDERS' because an existing locked ",
#         "object with that name is present in the namespace. ",
#         "Remove/rename the package-level object 'P_MODEL_BUILDERS' defined in R/*.R.",
#         call. = FALSE
#       )
#     }
#     rm("P_MODEL_BUILDERS", envir = ns)
#   }
#
#   makeActiveBinding(
#     sym = "P_MODEL_BUILDERS",
#     fun = function(value) {
#       if (!missing(value)) {
#         stop("`P_MODEL_BUILDERS` is read-only; use register_P_model() instead.", call. = FALSE)
#       }
#       reg <- .dmstar_env$P_MODEL_BUILDERS
#       if (!is.environment(reg)) {
#         # Defensive fallback (should not happen if .onLoad ran cleanly)
#         reg <- new.env(parent = emptyenv())
#         .dmstar_env$P_MODEL_BUILDERS <- reg
#       }
#       as.list(reg, all.names = TRUE)
#     },
#     env = ns
#   )
#
#   invisible()
# }
