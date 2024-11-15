#' Pierre Blanche Lagoon radionuclides data
#'
#' @description Example of data used for serac
#'
#' @docType data
#'
#' @usage PB06
#'
#' @format rda file
#'
#' @keywords datasets
#'
#' @references Sabatier, P., Dezileau, L., Blanchemanche, P., Siani, G., Condomines, M., Bentaleb, I., Piquès, G., 2010b. Holocene variations of radiocarbon reservoir ages in a mediterranean lagoonal system. Radiocarbon 52, 91–102
#'
#' @source \href{https://www.sciencedirect.com/science/refhub/S0265-931X(20)30695-0/sref43}{Sabatier et al (2010)}
#'
#' @examples
#' # The files PB06 and PB06_proxy must be in a folder Cores/PB06 for this code to work (see Bruel and Sabatier 2020)
#' serac(name = "PB06", coring_yr = 2006, model = c("CFCS", "CRS", "CRS_pw"),
#'                                           plot_Pb = T, plot_Pb_inst_deposit = T, inst_deposit = c(315,
#'                                           350), SML = 30, plot_Cs = T, Cher = c(50, 60),
#'                                           Hemisphere = c("NH"), NWT = c(100, 120),
#'                                           suppdescriptor = T, descriptor_lab = c("Si/Al"),
#'                                           historic_d = c(315, 350), historic_a = c(1893),
#'                                           historic_n = c("1894 storm"), min_yr = 1870,
#'                                           dmax =c(350), plotpdf =TRUE, depth_forced_CRS =c
#'                                           (55, 105), age_forced_CRS =c(1986, 1963))
