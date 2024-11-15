#' Lake Saint-André shortlived radionuclide
#'
#' @description Example of data used for serac
#'
#' @docType data
#'
#' @usage SAN
#'
#' @format rda file
#'
#' @keywords datasets
#'
#' @references Sabatier, P., Poulenard, J., Fanget, B., Reyss, J.-L., Develle, A.-L., Wilhelm, B., Ployon, E., Pignol, C., Naffrechoux, E., Dorioz, J.-M., Montuelle, B., Arnaud, F., 2014. Long-term relationships among pesticide applications, mobility, and soil erosion in a vineyard watershed. Proc. Natl. Acad. Sci. Unit. States Am. 111, 15647–15652. https://doi.org/10.1073/pnas.1411512111.
#' (\href{https://doi.org/10.1073/pnas.1411512111.})
#'
#' @source \href{https://doi.org/10.1073/pnas.1411512111.}{PNAS}
#'
#' @examples
#' # The file SAN must be in a folder Cores/SAN for this code to work (see Bruel and Sabatier 2020)
#' serac(name = "SAN", coring_yr = 2011, model = c("CFCS", "CIC", "CRS", "CRS_pw"),
#'       plot_Pb =TRUE, sedchange =c(165, 260), plot_Am =T,
#'       plot_Cs =TRUE, Cher =c(195, 205), Hemisphere =c("NH"), NWT =c(285, 295), FF =c(315, 325),
#'       plotpdf =TRUE, depth_forced_CRS = c(200, 290, 320), age_forced_CRS =c(1986, 1963, 1955),
#'       archive_metadata =TRUE)
"SAN"
