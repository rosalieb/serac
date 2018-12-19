#' Lake Allos shortlived radionuclide data
#'
#' Example of data used for serac
#'
#' @docType data
#'
#' @usage data(grav)
#'
#' @format txt file
#'
#' @keywords datasets
#'
#' @references Wilhelm, B., et al. (2012) Quat. Res. 78, 1–12.
#' (\href{https://www.cambridge.org/core/journals/quaternary-research/article/1400-years-of-extreme-precipitation-patterns-over-the-mediterranean-french-alps-and-possible-forcing-mechanisms/D6DBE2AD889B14347B540CD4CAB0C6E2}{Cambridge Core})
#'
#' @source \href{https://www.sciencedirect.com/science/article/pii/S0033589412000294}{Science Direct}
#'
#' @examples
#' # Radionuclides data
#' data(serac_example_ALO09P12)
#'
#' # proxy data
#' data(serac_example_ALO09P12_proxy)
#'
#' View()
#' times <- attr(grav, "time")
#' phe <- grav$pheno
#' \donttest{iplotCurves(phe, times)}
"serac_example_ALO09P12"
