#' Information on the user
#'
#' @description Question the user on professional ID (name, ORCID, email, Institute), the first time serac is used.. These metadata are then included in the general metadata file output every time a model is created.
#'
#' @export
#' @keywords age-depth modelling
#' @keywords visualisation
#' @keywords reproducibility
#' @keywords data management
#' @examples user_infos()
#'

user_infos <- function()
{
  message("If you are using serac for the first time, this function will\ncreate a file with some metadata (all are optional!).")

  answer<-"1"

  # First, check whether a file already exists
  if(length(list.files(paste0(getwd(),"/Cores/"), pattern="serac_metadata*", full.names=TRUE))==1) {
    metadata <- read.delim(list.files(paste0(getwd(),"/Cores/"), pattern="serac_metadata*", full.names=TRUE), header=T, sep="")
    message("\nA file with metadata already exists in this folder, with the \ninformation below:\n")
    print(metadata)
    message("\nIf you want to replace it, type 1 in the console (then 'Enter'). \nPress any other key (then 'Enter') to exit the function.")
    answer <- readline()
  }

  if(answer=="1")  {
    message("\nEnter an username. (recommended: first and last name, e.g., Rosalie Bruel)")
    username <- readline()
    message("\nWhat is your affiliation?, e.g., University Savoie Mont-Blanc")
    affiliation <- readline()
    message("\nORCID number:")
    ORCID <- readline()
    message("\nemail:")
    email <- readline()
    metadata <- data.frame("SECTIONS"= c("username","affiliation","ORCID","email"),
                           "USER INFORMATIONS" = c(username,affiliation,ORCID,email))
    write.table(x = metadata, file = paste(getwd(),"/Cores/serac_metadata_",username,".txt",sep = ""),col.names = T, row.names = F)
    cat("\nYou're all set! Here are the information you entered\n ")
    cat("_______________________________________________________\n")
    print(metadata)
    cat("_______________________________________________________")
    message("\nThanks for filling up these information; as long as you don't delete this file\nor run the function again, the information you entered will be added\nto the metadata of your projects.")

  }

}

