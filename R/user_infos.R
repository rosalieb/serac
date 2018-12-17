# serac -
#
# first time user function - enter some metadata
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

#' Information on the user
#'
#' Question the user on professional ID (name, ORCID, email, Institute). These metadata are then included in the general metadata file output every time a model is created.
#'
#' @export
#' @keywords age-depth modelling
#' @keywords visualisation
#' @examples user_infos()
#'

user_infos <- function()
{
  message("If you are using serac for the first time, this function will\ncreate a file with some metadata (all are optional!).")
  message("\nEnter an username. (recommended: first and last name)")
  username <- readline()
  message("\nWhat is your affiliation?")
  affiliation <- readline()
  message("\nORCID number:")
  ORCID <- readline()
  message("\nemail:")
  email <- readline()
  message("\nThanks for filling up this file; as long as you don't delete this file\nor run the function again, the information you entered will be added\nto the metadata of your projects.")
  metadata <- matrix(c("SECTIONS","USER INFORMATIONS",
                       "username",username,
                       "affiliation",affiliation,
                       "ORCID",ORCID,
                       "email",email),
                     ncol = 2, byrow = T)
  write.table(x = metadata, file = paste(getwd(),"/Cores/serac_metadata_",username,".txt",sep = ""),col.names = F, row.names = F)
}

