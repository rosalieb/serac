#' Lake Iseo shortlived radionuclides
#'
#' @description Example of data used for serac
#'
#' @docType data
#'
#' @usage Iseo
#'
#' @format rda file
#'
#' @keywords datasets
#'
#' @references Rapuc, W., Sabatier, P., Arnaud, F., Palumbo, A., Develle, A.-L., Reyss, J.-L., Augustin, L., Régnier, E., Piccin, A., Chapron, E., Dumoulin, J.-P., von Grafenstein, U., 2019. Holocene-long record of flood frequency in the Southern Alps (Lake Iseo, Italy) under human and climate forcing. Global Planet. Change 175, 160–172. https://doi.org/10.1016/j.gloplacha.2019.02.010
#'
#' @source \href{https://doi.org/10.1016/j.gloplacha.2019.02.010}{Rapuc et al (2019)}
#'
#' @examples
#' # The files Iseo and Iseo_varves must be in a folder Cores/Iseo for this code to work (see Bruel and Sabatier 2020)
#' serac(name = "Iseo", coring_yr = 2010, model = c("CFCS", "CIC", "CRS"),
#'       plot_Pb =TRUE, plot_Am =TRUE, plot_Cs =TRUE, Cher =c(70, 75),
#'       Hemisphere =c("NH"), NWT =c(130, 140), FF =c(164, 173), varves =TRUE,
#'       plotpdf =TRUE, stepout =5)
"Iseo"
