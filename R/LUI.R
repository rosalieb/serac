#' Lake Luitel shortlived radionuclides
#'
#' @description Example of data used for serac
#'
#' @docType data
#'
#' @usage LUI
#'
#' @format rda file
#'
#' @keywords datasets
#'
#' @references Guédron, S., Amouroux, D., Sabatier, P., Desplanque, C., Develle, A.-L., Barre, J., Feng, C., Guiter, F., Arnaud, F., Reyss, J.L., Charlet, L., 2016. A hundred year record of industrial and urban development in French Alps combining Hg accumulation rates and isotope composition in sediment archives from Lake Luitel. Chem. Geol. 431, 10–19. https://doi.org/10.1016/j.chemgeo.2016.03.016
#'
#' @source \href{https://doi.org/10.1016/j.chemgeo.2016.03.016}{Guédron et al (2016)}
#'
#' @examples
#' @examples
#' # The file LUI must be in a folder Cores/LUI for this code to work (see Bruel and Sabatier 2020)
#' serac(name ="LUI", coring_yr =2012, model =c("CFCS", "CRS", "CRS_pw"),
#'       mass_depth =TRUE, plot_Pb =TRUE, plot_Cs = TRUE, Cher = c(115, 125),
#'       Hemisphere = c("NH"), NWT =c(285, 295), FF =c(305, 315), plotpdf = TRUE,
#'       depth_forced_CRS = c(115, 285, 305), age_forced_CRS =c(1986, 1963, 1955))
"LUI"
