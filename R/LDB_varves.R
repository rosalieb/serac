#' Lake Bourget varves data
#'
#' @description Example of data used for serac
#'
#' @docType data
#'
#' @usage LDB_varves
#'
#' @format rda file
#'
#' @keywords datasets
#'
#' @references
#'
#' @source \href{}{}
#'
#' @examples
#' # The files LDB and LDB_varves must be in a folder Cores/LDB for this code to work (see Bruel and Sabatier 2020)
#' serac(name = "LDB", coring_yr = 2004, model = c("CFCS"), plot_Pb =T,
#'       plot_Pb_inst_deposit = T, plot_Cs = T, plot_Am = T, Cher = c(75, 85),
#'       Hemisphere = c("NH"), NWT = c(172, 180),
#'       FF = c(220, 230), inst_deposit = c(197, 210),
#'       historic_d = c(197, 210), historic_a = c(1958),
#'       historic_n = c("earthquake 1958"), varves = T,
#'       plotpdf =T, stepout =1)
"LDB_varves"
