#' Prepare data input file
#'
#' @description Function to prepare data input file. The input data file can easily be prepared in excel or other programs. The function below may help the user in checking the right input format. Function only for the general data file, 'varves' and 'proxy' input data have to be formated separately.
#'
#' @export
#' @param name Name of the core, given using quotes. Use preferably the published name of the core for traceability.
#' @keywords visualisation
#' @examples serac_input_formatting('PB06')
#'

serac_input_formatting <- function(name)
{
  answer=NULL;answer2=NULL
  message("This function will help you to format the input file for serac.")
  # read the data
  dt <- read.delim(file = paste(getwd(),"/Cores/",name,"/",name,".txt", sep=""))
  ncol_init <- ncol(dt)
  # write in the folder the raw data for archive
  if(length(list.files(paste(getwd(),"/Cores/",name,"/", sep=""), pattern=paste(name,"_raw.txt",sep=""), full.names=TRUE))==0) write.table(x = dt, file = paste(getwd(),"/Cores/",name,"/",name,"_raw.txt", sep=""),col.names = T, row.names = F, sep="\t")

  # Look at what the data look like
  message("This is what your raw input data file looks like (head and tail):")
  print(head(dt,4))
  cat("              [...........]\n")
  print(tail(dt,4))

  # multiplication factor if wrong unit used
  message("Note that serac uses mm as a default depth unit. \nIf your unit is different, what is the converter to get mm?\n e.g. 1 cm= 10 mm, enter 10 in the console\n Enter 1 if no conversion is needed.")
  multipli <- as.numeric(readline())

  # Mean depth
  message("\nWhich column has the average depth information? Enter the number \ncorresponding to the column name according to the following table.\nWrite 0 if you don't have this variable in your input file.")
  print(as.data.frame(matrix(c(colnames(dt)[1:ncol_init],1:ncol_init), ncol=2, byrow = F)))
  answer2 <- as.numeric(readline())
  if (answer2!=0) dt$Depth <- dt[,answer2]*multipli else dt$Depth <- rep(NA,nrow(dt))

  # Min depth
  message("\nWhich column has the minimum depth information (top section)? \nEnter the number corresponding to the column name according \nto the following table.\nWrite 0 if you don't have this variable in your input file.")
  print(as.data.frame(matrix(c(colnames(dt)[1:ncol_init],1:ncol_init), ncol=2, byrow = F)))
  answer2 <- as.numeric(readline())
  if (answer2!=0) dt$depth_min <- dt[,answer2]*multipli else dt$depth_min <- rep(NA,nrow(dt))

  # Max depth
  message("\nWhich column has the maximum depth information (bottom section)? \nEnter the number corresponding to the column name according \nto the following table.\nWrite 0 if you don't have this variable in your input file.")
  print(as.data.frame(matrix(c(colnames(dt)[1:ncol_init],1:ncol_init), ncol=2, byrow = F)))
  answer2 <- as.numeric(readline())
  if (answer2!=0)  dt$depth_max <- dt[,answer2]*multipli else dt$depth_max <- rep(NA,nrow(dt))

  # Thickness
  message("\nWhich column has the thickness information? Enter the number \ncorresponding to the column name according to the following table.\nWrite 0 if you don't have this variable in your input file.")
  print(as.data.frame(matrix(c(colnames(dt)[1:ncol_init],1:ncol_init), ncol=2, byrow = F)))
  answer2 <- as.numeric(readline())
  if (answer2!=0)  {
    dt$thickness <- dt[,answer2]*multipli
  } else {
    if (all(is.na(dt$depth_max))) stop("You need to enter either the min/max depth for each sample OR the thickness!")
    dt$thickness <- dt$depth_max-dt$depth_min
  }

  # Change min-max depth if NAs
  if (all(is.na(dt$depth_min))&all(is.na(dt$Depth))) cat(' Warning, you need to include the depth in your input file. \n You have the choice between average depth + thickness OR minimun and maximum depth of a sample.')
  if (all(is.na(dt$depth_max))&all(is.na(dt$Depth))) cat(' Warning, you need to include the depth in your input file. \n You have the choice between average depth + thickness OR minimun and maximum depth of a sample.')
  if (all(is.na(dt$thickness))&all(is.na(dt$Depth))) cat(' Warning, you need to include the depth in your input file. \n You have the choice between average depth + thickness OR minimun and maximum depth of a sample.')

  if (all(is.na(dt$depth_min))) dt$depth_min <- dt$Depth-dt$thickness/2
  if (all(is.na(dt$depth_max))) dt$depth_max <- dt$Depth+dt$thickness/2

  # Density
  message("\nDo you have the density information in your table and/or \ndo you wish to calculate it (write Y for yes or N for no \nin the console)?\nDensity is needed for inventory calculation and CRS model.")
  answer <- readline()
  if (answer=="y") answer="Y"
  if (answer=="n") answer="N"
  if (answer=="Y") {
    message("\nIs density already in your table? (Y/N)")
    answer <- readline()
    if (answer=="y") answer="Y"
    if (answer=="n") answer="N"
    if (answer=="Y") {
      message("\nWhich column has the density information? Enter the number \ncorresponding to the column name according to the following table.")
      print(as.data.frame(matrix(c(colnames(dt)[1:ncol_init],1:ncol_init), ncol=2, byrow = F)))
      answer2 <- as.numeric(readline())
      if (answer2!=0)  dt$density <- dt[,answer2]
    } else {
      message("\nWhich column has the dry sediment mass (g/cm3) information? \nEnter the number corresponding to the column name according \nto the following table.\nWrite 0 if you don't have this variable in your input file.")
      print(as.data.frame(matrix(c(colnames(dt)[1:ncol_init],1:ncol_init), ncol=2, byrow = F)))
      answer2 <- as.numeric(readline())
      if (answer2!=0)  dt$drymass <- dt[,answer2] else dt$drymass <- rep(NA,nrow(dt))
      if (all(is.na(dt$drymass))) {
        cat(" Warning - density was not calculated: need dry mass.")
        dt$density <- rep(NA,nrow(dt))
      } else {
        message("\nWhat is the diameter (in mm) of the core? e.g. 63")
        diameter <- as.numeric(readline())
        message("\nWhat portion of the core was prepared for dry mass? \ne.g. enter in the console 1/2 or 0.5 if half of the core.")
        subsetcore <- as.numeric(readline())
        dt$density <- dt$drymass/(pi*(diameter/2)^2*dt$thickness*subsetcore)
      }
    }
  } else {dt$density <- rep(NA,nrow(dt))}


  # 210 Pb ex
  message("\nWhich column has the 210Pb(excess) information? Enter the number \ncorresponding to the column name according to the following table.\nWrite 0 if you don't have this variable in your input file.")
  print(as.data.frame(matrix(c(colnames(dt)[1:ncol_init],1:ncol_init), ncol=2, byrow = F)))
  answer2 <- as.numeric(readline())
  if (answer2!=0)  dt$Pb210ex <- dt[,answer2] else dt$Pb210ex <- rep(NA,nrow(dt))

  # 210 Pb ex error
  message("\nWhich column has the 210Pb(excess) ERROR information? Enter the number \ncorresponding to the column name according to the following table.\nWrite 0 if you don't have this variable in your input file.")
  print(as.data.frame(matrix(c(colnames(dt)[1:ncol_init],1:ncol_init), ncol=2, byrow = F)))
  answer2 <- as.numeric(readline())
  if (answer2!=0)  dt$Pb210ex_er <- dt[,answer2] else dt$Pb210ex_er <- rep(NA,nrow(dt))

  # 137Cs
  message("\nWhich column has the 137Cs information? Enter the number \ncorresponding to the column name according to the following table.\nWrite 0 if you don't have this variable in your input file.")
  print(as.data.frame(matrix(c(colnames(dt)[1:ncol_init],1:ncol_init), ncol=2, byrow = F)))
  answer2 <- as.numeric(readline())
  if (answer2!=0)  dt$Cs137 <- dt[,answer2] else dt$Cs137 <- rep(NA,nrow(dt))

  # 137Cs error
  message("\nWhich column has the 137Cs ERROR information? Enter the number \ncorresponding to the column name according to the following table.\nWrite 0 if you don't have this variable in your input file.")
  print(as.data.frame(matrix(c(colnames(dt)[1:ncol_init],1:ncol_init), ncol=2, byrow = F)))
  answer2 <- as.numeric(readline())
  if (answer2!=0)  dt$Cs137_er <- dt[,answer2] else dt$Cs137_er <- rep(NA,nrow(dt))

  # 241Am
  message("\nWhich column has the 241Am information? Enter the number \ncorresponding to the column name according to the following table.\nWrite 0 if you don't have this variable in your input file.")
  print(as.data.frame(matrix(c(colnames(dt)[1:ncol_init],1:ncol_init), ncol=2, byrow = F)))
  answer2 <- as.numeric(readline())
  if (answer2!=0)  dt$Am241 <- dt[,answer2] else dt$Am241 <- rep(NA,nrow(dt))

  # 241Am error
  message("\nWhich column has the 241Am ERROR information? Enter the number \ncorresponding to the column name according to the following table.\nWrite 0 if you don't have this variable in your input file.")
  print(as.data.frame(matrix(c(colnames(dt)[1:ncol_init],1:ncol_init), ncol=2, byrow = F)))
  answer2 <- as.numeric(readline())
  if (answer2!=0)  dt$Am241_er <- dt[,answer2] else dt$Am241_er <- rep(NA,nrow(dt))

  # Final Output file
  dt2 <- data.frame("depth_min"= dt$depth_min,
                    "depth_max"= dt$depth_max,
                    "thickness"= dt$thickness,
                    "density"= dt$density,
                    "Pb210ex"= dt$Pb210ex,
                    "Pb210ex_er"= dt$Pb210ex_er,
                    "Cs137"= dt$Cs137,
                    "Cs137_er"= dt$Cs137_er,
                    "Am241"= dt$Am241,
                    "Am241_er"= dt$Am241_er)

  # Delete all NA columns and all NA rows
  dt2 <- dt2[,colSums(is.na(dt2))<nrow(dt2)]
  dt2 <- dt2[rowSums(is.na(dt2))<ncol(dt2),]

  # write final input file in the folder
  write.table(x = dt2, file = paste(getwd(),"/Cores/",name,"/",name,".txt", sep=""),col.names = T, row.names = F, sep="\t")
  # View result
  cat("_______________________________________________________")
  if(nrow(dt2)<=15) {
    message(paste("\nYou're all set! Below is the file you created. Now run serac. \nReminder, minimum function would be: serac('",name,"', coring_yr = ",format(Sys.Date(), "%Y"),") \nEdit coring year accordingly.",sep=""))
    print(dt2)
  } else {
    message(paste("\nYou're all set! Below is the file you created (first 15 rows). \nNow run serac.\nReminder, minimum function would be: serac('",name,"', coring_yr = ",format(Sys.Date(), "%Y"),") \nEdit coring year accordingly.",sep=""))
    print(dt2[1:15,])
    cat("              [...........]\n")
  }
  cat("_______________________________________________________")
}
