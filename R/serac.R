#' serac age-depth modelling function
#'
#' This is the main age-depth modelling function. The default values can be changed permanently within this file or temporarily when calling serac(). If there is any options you would like to see included in future version, please contact one of the authors
#'
#' @export
#' @param name Name of the core, given using quotes. Defaults to the core provided with serac. Use preferably the published name of the core for traceability.
#' @param coring_yr Year of coring.
#' @param model Select 1 to 3 item between c("CFCS", "CIC", "CRS"). If several models are selected, they will all be plotted together in the last window.
#' @param Cher If 137Cs measurement were done, where do you detect the Chernobyl peak? The argument is a vector of two depth given in millimeters giving the top and bottom threshold for the 1986 Chernobyl event. The user can run the model without giving any specification before making a decision. In such case, leave the argument empty. Note that the two depths needs to represent a sample, or more than a sample.
#' @param NWT If 137Cs measurement were done, where do you detect the Nuclear Weapon Test peak? The argument is a vector of two depth given in millimeters giving the top and bottom threshold for the 1960s Nuclear Weapon Test event. The user can run the model without giving any specification before making a decision. In such case, leave the argument empty. Note that the two depths needs to represent a sample, or more than a sample.
#' @param Hemisphere Chose between North Hemisphere "NH" and South Hemisphere "SH" depending on the location of your system. This argument is required if you chose to plot NWT.
#' @param FF If 137Cs measurement were done, where do you detect the First Fallout period? The argument is a vector of two depth given in millimeters giving the top and bottom threshold for the First Fallout period in 1955. The user can run the model without giving any specification before making a decision. In such case, leave the argument empty. Note that the two depths needs to represent a sample, or more than a sample.
#' @param inst_deposit Upper and lower depths (in mm) of sections of abrupt accumulation that inst_deposit c() should be excised, e.g., c(100, 120, 185, 195) for two sections of 10.0-12.0 cm and 18.5-19.5 cm depth
#' @param input_depth_mm Units for SML, radionuclides peaks (Chernobyl, Nuclear War Tests, and First Fallouts), instantaneous deposits, depths to ignore. Default is TRUE, and inputs must be in mm - unless "input_depth_mm=F". If turned to FALSE, input can be given in g/cm2.
#' @param ignore The depth (in mm - unless "input_depth_mm=F") of any sample that should be ignored from the age-depth model computation, e.g., c(55) will remove the measurement done at 5.5 cm. The data will be ploted by default in grey on the output graph (you can change this with the inst_depositcol argument)
#' @param plotpdf Logical argument to indicate whether you want the output graph to be saved to your folder.
#' @param preview Logical argument to indicate whether you want the output graph to be ploted. Default is TRUE, and the graph is ploted within your R session. It might be convenient to turn this argument to FALSE if errors keep coming telling you your R window is too small.
#' @param plotphoto Logical argument to indicate whether you want to plot the photo of the core along your age-model. If plotphoto=TRUE, you need to indicate the upper and lower limit of the photo in mm - unless "input_depth_mm=F" in following arguments.
#' @param minphoto Mandatory if plotphoto=TRUE. Lower limit of the core photo in mm - unless "input_depth_mm=F", e.g., minphoto=0 indicates that the photo starts at 0 mm. The photo will automatically be truncated acording to the minimum and maximum depth of the age model given in other arguments.
#' @param maxphoto Mandatory if plotphoto=TRUE. Upper limit of the core photo in mm - unless "input_depth_mm=F", e.g., maxphoto=320 indicates that the photo ends at 32 cm. The photo will automatically be truncated acording to the minimum and maximum depth of the age model given in other arguments.
#' @param Pbcol Vector of color to plot 210Pbex data. If length(Pbcol)>1, the different colors will be used to plot the different slopes in between change(s) in sedimentation rate. Example of color vector: Pbcol=c("black","midnightblue","darkgreen").
#' @param inst_depositcol A color to plot the data points within instantaneous deposit or ignored data. Example: inst_depositcol=grey(0.85).
#' @param modelcol Vector of color to plot different model if length(model)>1. If length(modelcol)>1, the different colors will be used to plot the different change in sedimentation rate. Example of color vector: modelcol=c("black","red","darkorange") to plot "CFCS", "CIC", "CRS" models in this order.
#' @param historic_d Vector with upper and lower depth of historical event(s), e.g., historic_d=c(120,130) will identify the event between 12 and 13 cm on the last window with the age model.
#' @param historic_a Vector of year of different historical events, e.g., historic_a=c(1895) will add the point 1895 on the last window with the age model. Historical events can be older than the dated section, in which case the depth is obtained from the model if historic_d is not specified. historic_a is a vector twice as short as historic_d, as each age correspond to an upper+lower limit in the vector 'historic_d'. If not all ages are known, put NA in the vector, e.g., historic_a=c(NA,1895)
#' @param historic_n Vector of names of different historical events, e.g., historic_n=c("1895 flood"). Optional. If you plot several events, and don't want to plot all the names, add a NA in the vector, e.g., historic_n=c(NA,"1895 flood") will understand that the first event doesn't have a name, but the second does.
#' @param historic_test Visualisation tool for known ages. This argument will plot a vertical line in the last window (the one with the age-depth model). Can be useful when the user know specific ages that may have resulted in changes in sedimentation rates. E.g., historic_test=c(1996).
#' @param suppdescriptor Up to two supplementary descriptor(s) to plot in an additional window. Logical argument. The decision on ploting more than one supplementary descriptor depends on the length of the vector descriptor_lab. An additional input file with these data should be included in the folder with the initial data.
#' @param descriptor_lab Label used on the axis, e.g., descriptor_lab=c("LOI", "Ca/Fe") if two supplementary descriptors are specified.
#' @param suppdescriptorcol Vector of color to plot different descriptor if length(descriptor_lab)>1. If length(descriptor_lab)>1, the different colors will be used to plot the different change in sedimentation rate. Example of color vector: suppdescriptorcol=c("black","purple").
#' @param plot_Am Logical argument indicating whether or not serac should plot 241Am.
#' @param plot_Cs Logical argument indicating whether or not serac should plot 137Cs.
#' @param plot_Pb Logical argument indicating whether or not serac should plot 210Pbex.
#' @param plot_Pb_inst_deposit Logical argument indicating whether or not serac should plot 210Pbex without instantaneous deposit. If TRUE, inst_deposit shouldn't be a null vector.
#' @param plot_CFCS_regression Whether to plot or not the linear regression. If the parameter is not specified, it will automatically turn to TRUE, but will also automatically turn to FALSE if instantaneous deposit are present but the argument 'plot_Pb_inst_deposit' is turned to FALSE. Linear regression won't match if there are some instantaneous deposit. In other words, in most cases, the user shouldn't need to modify this parameter.
#' @param varves Logical argument to indicate whether varve counting results should be ploted on the last window. An additional input file with these data should be included in the folder with the initial data.
#' @param dmin Maximum depth of age-depth model (useful if the user doesn't want to plot the lower region).
#' @param dmax Maximum depth of age-depth model (useful if the user doesn't want to plot the lower region). dmax cannot be in the middle of an instantaneous deposit. e.g. if there is an instantaneous deposit between 180 and 200 mm, dmax cannot be 190 mm, and will be converted to 200 mm automatically.
#' @param sedchange Up to two changes in sedimentation rate, e.g., sedchange=c(175,290) indicates two changes of sedimentation rate at 17.5 and 29.0 cm.
#' @param min_yr The minimum year limit for the age-depth model plot. The user can adjust this argument after a first computation of the model
#' @param SML Surface Mixed Layer: a depth in mm - unless "input_depth_mm=F" above which the sediment is considered to be mixed. E.g., SML=30 indicates that the first 3 cm are mixed sediment: the data point are ploted but not included in the Pb models.
#' @param stepout Depth resolution for the file out in mm - unless "input_depth_mm=F".
#' @param mycex Graphical parameter: a multiplication factor to increase (mycex>1) ou decrease (mycex<1) label sizes.
#' @param archive_metadata Logical argument. If TRUE, require fields regarding the measurements on the core. Allows missing information; just press 'ENTER' in your computer (leave an empty field).
#' @param mass_depth Logical argument. If TRUE, require density, and will plot the radionuclides against massic depth. Core photo and supplementary descriptor are not available under this option.
#' @keywords age-depth modelling
#' @keywords visualisation
#' @examples
#' # Lake Bourget
#' # serac(name="LDB",coring_yr=2004)
#' # serac(name="LDB",coring_yr=2004,model=c("CFCS"),plotphoto=TRUE,minphoto=c(0),maxphoto=c(370),plot_Pb=T,plot_Pb_inst_deposit=T,plot_Cs=T,plot_Am=T,Cher=c(75,85),Hemisphere=c("NH"),NWT=c(172,180),inst_deposit=c(197,210),historic_d=c(197,210),historic_a=c(1958),historic_n=c("earthquake 1958"),varves=T,plotpdf=T,preview=T,stepout=1)
#'
#' # Lake Iseo
#' # serac(name="Iseo",coring_yr=2010)
#' # serac(name="Iseo",coring_yr=2010,model=c("CFCS","CIC","CRS"),plotphoto=TRUE,minphoto=c(0),maxphoto=c(320),plot_Pb=T,plot_Am=T,plot_Cs=T,Cher=c(70,75),Hemisphere=c("NH"),NWT=c(130,140),FF=c(164,173),varves=TRUE,plotpdf=T,preview=T,stepout=5)
#'
#' # Lake Saint-Andre
#' # serac(name="SAN",coring_yr=2011,model=c("CFCS"),plotphoto=TRUE,minphoto=c(0),maxphoto=c(420),plot_Pb=T,plot_Am=T,plot_Cs=T,Cher=c(195,205),Hemisphere=c("NH"),NWT=c(275,295),FF=c(315,325),sedchange=c(165,260),plotpdf=TRUE)
#'
#' # Lake Allos
#' # serac(name="ALO09P12",coring_yr=2009,model=c("CFCS"),plotphoto=TRUE,minphoto=c(0),maxphoto=c(210),plot_Pb=T,plot_Am=T,plot_Cs=T,Cher=c(30,40),Hemisphere=c("NH"),NWT=c(51,61),sedchange=c(75.5),plot_Pb_inst_deposit=T,inst_deposit=c(20,28,100,107,135,142,158,186),suppdescriptor=TRUE,descriptor_lab=c("Ca/Fe"),historic_d=c(20,28,100,107,135,142,158,186),historic_a=c(1994,1920,1886,1868),historic_n=c("sept1 994 flood","1920 flood","1886 flood","1868 flood ?"), min_yr=c(1750),dmax=c(180), plotpdf=TRUE,preview=F)
#'
#' # Pierre-Blanche lagoon
#' # serac(name="PB06",coring_yr=2006,model=c("CFCS","CRS"),plotphoto=TRUE,minphoto=c(0),maxphoto=c(350),plot_Pb=T,plot_Cs=T,Cher=c(50,60),Hemisphere=c("NH"),NWT=c(100,120),suppdescriptor=T,descriptor_lab=c("Si/Al"),SML=30,inst_deposit=c(315,350),historic_d=c(315,350),historic_a=c(1893),historic_n=c("1894storm"),min_yr=1870,dmax=c(350),plotpdf=TRUE)
#'

serac <- function(name="", model=c("CFCS"),Cher=NA,NWT=NA,Hemisphere=NA,FF=NA,inst_deposit=c(0),
                  input_depth_mm = T,ignore=c(), mass_depth=FALSE,
                  plotpdf=FALSE,preview=TRUE,plotphoto=FALSE,minphoto=c(),maxphoto=c(),
                  Pbcol=c("black","midnightblue","darkgreen"),inst_depositcol=grey(0.85),
                  modelcol=c("black","red","darkorange"),
                  historic_d=NA,historic_a=NA,historic_n=NA,historic_test=NA,
                  suppdescriptor=FALSE,descriptor_lab=c(),suppdescriptorcol=c("black","purple"),
                  coring_yr=c(),plot_Am=FALSE,plot_Cs=FALSE,plot_Pb=TRUE,
                  plot_Pb_inst_deposit=FALSE,plot_CFCS_regression=c(),
                  varves=FALSE, dmin=c(),dmax=c(),sedchange=c(0),
                  min_yr=1880, SML=c(0), stepout=5, mycex=1,
                  archive_metadata=FALSE, save_code=TRUE)
.serac(name, model,Cher,NWT,Hemisphere,FF,inst_deposit,
       input_depth_mm,ignore,mass_depth,
       plotpdf,preview,plotphoto,minphoto,maxphoto,
       Pbcol,inst_depositcol,
       modelcol,
       historic_d,historic_a,historic_n,historic_test,
       suppdescriptor,descriptor_lab,suppdescriptorcol,
       coring_yr,plot_Am,plot_Cs,plot_Pb,
       plot_Pb_inst_deposit,plot_CFCS_regression,
       varves, dmin,dmax,sedchange,
       min_yr, SML,stepout, mycex,
       archive_metadata,save_code)

.serac <- function(name, model,Cher,NWT,Hemisphere,FF,inst_deposit,
                   input_depth_mm,ignore,mass_depth,
                   plotpdf,preview,plotphoto,minphoto,maxphoto,
                   Pbcol,inst_depositcol,
                   modelcol,
                   historic_d,historic_a,historic_n,historic_test,
                   suppdescriptor,descriptor_lab,suppdescriptorcol,
                   coring_yr,plot_Am,plot_Cs,plot_Pb,
                   plot_Pb_inst_deposit,plot_CFCS_regression,
                   varves, dmin,dmax,sedchange,
                   min_yr, SML,stepout, mycex,
                   archive_metadata,save_code) {


  # 0. INITIALIZE ####
  # Calculate how long the function took to run
  old_time <- Sys.time() # get start time

  # load packages
  pkgTest("Hmisc")
  pkgTest("jpeg")
  pkgTest("TeachingDemos")

  # create empty list where outputs will be saved
  out_list <- list()

  # Archive metadata
  if(archive_metadata) core_metadata(name=name)

  # warn and stop if abnormal settings are provided
  # coring year is one of the two mandatory argument
  if(is.null(coring_yr))     stop("\n Warning, please enter the 'coring_yr'.\n\n")

  # serac support 2 changes in sedimentation rates maximum
  if(length(sedchange)>2)    stop("\n Warning, serac only support two changes in sedimentation rate. Please check the manual.\n\n", call.=FALSE)

  # Chernobyl fallouts are dated at different years depending on the hemisphere
  if(!is.na(Cher) && (Hemisphere=="NH"|Hemisphere=="SH")==FALSE)  stop("\n Warning, please select the hemisphere where your system is located.\n 'NH' or 'SH' for Northern/Southern Hemisphere.\n\n")

  # if the argument plotphoto is true, then a photo with the exact same name and the extension .jpg must be provided in the folder
  # min and max depth must be provided so the image is automatically cropped.
  if(plotphoto==TRUE) {
    if(!file.exists(paste(getwd(),"/Cores/",name,"/",name,".jpg", sep=""))) stop("\n Warning, you asked to include the photo of the core but it was not found in the repository.\n Check the name (must be the same than input data name) and the extension (must be .jpg).\n\n")
    if(is.null(minphoto) | (is.null(maxphoto)))                             stop("\n Warning, you need to indicate upper (minphoto) and lower (maxphoto) depth_avg of the core (mm).\n\n")
  }

  # depth must be provided by pair (upper and lower depths)
  if(length(historic_d)>=1 && !all(is.na(historic_d))) {
    if (length(historic_a) != length(historic_d)/2) stop("\n Warning, length(historic_a) != length(historic_d)/2 \n Depths for historic events must be given by pair - upper and lower depth.\n Read the help section.\n\n")
  }

  # specific requirements to run CIC model.
  if(any(model=="CIC")) {
    if(!is.null(inst_deposit)&&max(inst_deposit)>0) cat("\n Warning, in most of the situations, CIC model should not be run if you assume the\n presence of instantaneous deposit. Be cautious while interpreting this output. \n\n")
    if(SML>0)                                       stop("\n Warning, CIC model should not be run if you assume the presence of a surface mixed layer. \n\n")
  }

  #### 1. READ DATA ----
  dt <- read.delim(file = paste(getwd(),"/Cores/",name,"/",name,".txt", sep=""))
  dt <- dt[,colSums(is.na(dt))<nrow(dt)]
  if(plotphoto) photo <- readJPEG(paste(getwd(),"/Cores/",name,"/",name,".jpg", sep=""))
  if(varves) varve <- read.delim(file = paste(getwd(),"/Cores/",name,"/",name,"_varves.txt", sep=""))
  if(suppdescriptor) dt_suppdescriptor <- read.delim(file = paste(getwd(),"/Cores/",name,"/",name,"_proxy.txt", sep=""))

  # 1.1. Flexibility in input columns format ####
  # I'm adding '[1]' at the end of each grep expressions, in case several column carry the name
  # This shouldn't ever cause an issue to the user, as long as they use a clean input data file
  #         (without overlap in column names - there should be only one column with the 'key'
  #         informations). Any other column can be added and won't be read by serac if they
  #         don't contain the keywords used below.
  # Depth columns
  if (length(intersect(grep("epth|EPTH",colnames(dt)),grep("top|bottom|min|max|mass",colnames(dt),invert=TRUE)))>=1){
    dt$depth_avg <- dt[,intersect(grep("epth|EPTH",colnames(dt)),grep("top|bottom|min|max|mass",colnames(dt),invert=TRUE))[1]]
  }
  if (length(grep("hickness|HICKNESS",colnames(dt))>=1)) {
    dt$thickness <- dt[,grep("hickness|HICKNESS",colnames(dt))[1]]
  }

  # For 210Pbex
  if(length(grep("Pb",colnames(dt)))>1) {
    dt$Pbex <- dt[,intersect(intersect(grep("Pb",colnames(dt)),grep("ex",colnames(dt))),grep("er",colnames(dt),invert=TRUE))[1]]
    dt$Pbex_er <- dt[,intersect(intersect(grep("Pb",colnames(dt)),grep("ex",colnames(dt))),grep("er",colnames(dt)))[1]]
  } else if(plot_Pb|plot_Pb_inst_deposit) message("\n Warning. We did not find the Lead column (+ error) in the input file.\n\n")
  # For 137Cs
  if(length(grep("Cs",colnames(dt)))>1) {
    dt$Cs <- dt[,intersect(grep("Cs",colnames(dt)),grep("er",colnames(dt),invert=TRUE))[1]]
    dt$Cs_er <- dt[,intersect(grep("Cs",colnames(dt)),grep("er",colnames(dt)))[1]]
  } else {
    dt$Cs <- rep(NA, nrow(dt))
    dt$Cs_er <- rep(NA, nrow(dt))
    if(plot_Cs) message("\n Warning. We did not find the Cesium column (+ error) in the input file.\n\n")
  }
  # For 241Am
  if(length(grep("Am",colnames(dt)))>1) {
    dt$Am <- dt[,intersect(grep("Am",colnames(dt)),grep("er",colnames(dt),invert=TRUE))[1]]
    dt$Am_er <- dt[,intersect(grep("Am",colnames(dt)),grep("er",colnames(dt)))[1]]
  } else {
    dt$Am <- rep(NA, nrow(dt))
    dt$Am_er <- rep(NA, nrow(dt))
    if(plot_Am) message("\n Warning. We did not find the Americium column (+ error) in the input file.\n\n")
  }

  # 1.2. Calculate thickness if missing, or conversely calculate upper and lower depth of samples ####
  if(is.null(dt$thickness) & is.null(dt$depth_top) & is.null(dt$depth_bottom) & is.null(dt$depth_min) & is.null(dt$depth_max)) stop("\n Warning, please indicate the thickness of each sample in mm or the top and \nbottom section of each sample so we can compute it for you.\n\n")

  # Change dt$depth_min/dt$depth_max by top and bottom - correction to international format
  if(length(dt$depth_min)==nrow(dt)) dt$depth_top <- dt$depth_min
  if(length(dt$depth_max)==nrow(dt)) dt$depth_bottom <- dt$depth_max

  if(length(dt$depth_avg)<nrow(dt)) dt$depth_avg <- (dt$depth_top+dt$depth_bottom)/2

  if(length(dt$thickness)<nrow(dt)) {
    dt$thickness <- rep(NA,nrow(dt))
    for (i in 1:nrow(dt)) dt$thickness[i] <- (abs(dt$depth_top[i]-dt$depth_bottom[i]))
  }

  if(is.null(dt$depth_top))         dt$depth_top <- dt$depth_avg-dt$thickness/2
  if(is.null(dt$depth_bottom))      dt$depth_bottom <- dt$depth_avg+dt$thickness/2

  # Same for varves file, if present
  if(varves) {
    varve$depth_avg <- varve[,grep("epth",colnames(varve))]
    varve$thickness <- varve[,grep("hickness",colnames(varve))]

    if(is.null(varve$thickness) & is.null(varve$depth_top) & is.null(varve$depth_bottom)) stop("\n Warning, please indicate in the 'varves' input data file the thickness \nof each sample in mm or the top and bottom section of each sample so \nwe can compute it for you.\n\n")

    if(is.null(varve$thickness)) {
      varve$thickness <- rep(NA,nrow(varve))
      for (i in 1:nrow(varve)) varve$thickness[i] <- (abs(varve$depth_top[i]-varve$depth_bottom[i]))
    }
  }

  # Additional warning not so related to this section
  # If density is missing, cannot calculate CRS
  if(any(model=="CRS")) if(is.null(dt$density)) stop("\n Warning, you need to include the density in g/cm2 for each sample to compute CRS model.\n\n")
  if(mass_depth) if(is.null(dt$density)) stop("\n Warning, you need to include the density in g/cm2 for each sample to calculate mass accumulation rate.\n\n")

  # Additional warning if density is not given continuously
  for (i in 2:nrow(dt)) {
    if(i==2) test1<-NULL
    test1 <- c(test1,dt$depth_top[i]-dt$depth_bottom[i-1])
  }
  if(!is.null(dt$density) && !all(test1==0)) {
    message(paste0("\n Warning, density is not given continuously for the whole core.\n Inventories, CRS model, and mass accumulation rate should be\n interpreted very carefully. Alternatively, enter the density\n for the whole core.\n Density is missing for ",length(test1[test1!=0])," layer(s) missing at: "))
    message(paste0("     - ",dt$depth_bottom[which(test1>0)],"-",dt$depth_top[which(test1>0)+1]," mm \n"))
  }

  # 1.3. Fill in missing data ####
  # 1.4 Complete core vector ####
  # Create the vector complete_core_depth when the measurements haven't been done for all the layers.
  # Necessary for inventory for instance.
  complete_core_temporary <- c(0,dt$depth_top, dt$depth_bottom,inst_deposit)
  # Remove the instantaneous deposit depth deeper than measured depth.
  # We do not want to extrapolate 210Pbex and 137Cs below the actual measurement.
  complete_core_temporary <- complete_core_temporary[complete_core_temporary<=max(dt$depth_bottom,na.rm=T)]
  complete_core_temporary <- unique(complete_core_temporary)
  complete_core_temporary <- complete_core_temporary[order(complete_core_temporary, decreasing = F)]
  complete_core_thickness <- NULL
  complete_core_depth <- NULL
  complete_core_depth_top <- NULL
  complete_core_depth_bottom <- NULL
  for (i in 2:length(complete_core_temporary)) {
    complete_core_thickness <- c(complete_core_thickness,complete_core_temporary[i]-complete_core_temporary[i-1])
    complete_core_depth <- c(complete_core_depth, complete_core_temporary[i]-complete_core_thickness[i-1]/2)
    complete_core_depth_top <- c(complete_core_depth_top, complete_core_temporary[i-1])
    complete_core_depth_bottom <- c(complete_core_depth_bottom, complete_core_temporary[i])
  }
  rm(complete_core_temporary)

  # Generated the complete 210Pbex and 137Cs profile (in case the sampling was not continuous)
  complete_core_Pbex <- approx(x= dt$depth_avg, dt$Pbex, xout= complete_core_depth, rule = 2)$y
  complete_core_Pbex_err <- approx(x= dt$depth_avg, dt$Pbex_er, xout= complete_core_depth, rule = 2)$y
  if(any(!is.na(dt$Cs))) complete_core_Cs <- approx(x= dt$depth_avg, dt$Cs, xout= complete_core_depth, rule = 2)$y
  if(any(!is.na(dt$Cs))) complete_core_Cs_err <- approx(x= dt$depth_avg, dt$Cs_er, xout= complete_core_depth, rule = 2)$y

  # Generate the complete density (in case the sampling was not continuous).
  # It is just a linear interpolation.
  if(length(grep("density",x = colnames(dt)))>=1) complete_core_density <- approx(x= dt$depth_avg, dt$density, xout= complete_core_depth, rule = 2)$y

  # 1.5. Which keep ####
  # When calculating the inventories, we don't want to take in account the depth included
  # in a instantaneous deposit, or the depth explicitely 'ignored' by the operator.
  # Here, I'm creating an index of data that will be ignored from the inventory calculation.
  myvec <- 1:length(complete_core_depth)
  whichkeep <- NULL
  if(exists("inst_deposit") && max(length(inst_deposit))>0 && length(inst_deposit) %% 2 == 0){
    # Not the best way to do it but I need this vector early
    myseq <- seq(1, length(inst_deposit), by = 2)
    for (i in seq(1, length(inst_deposit), by = 2)) {
      if (i == myseq[1]) whichkeep <- c(whichkeep, myvec[complete_core_depth<inst_deposit[i]])
      whichkeep <- c(whichkeep, myvec[complete_core_depth>inst_deposit[i-1] & complete_core_depth<inst_deposit[i]])
      if (i == rev(myseq)[1]) whichkeep <- c(whichkeep, myvec[complete_core_depth>inst_deposit[i+1]])
    }
    rm(myseq)
  } else whichkeep <- myvec
  whichkeep <- whichkeep[!whichkeep %in% which(complete_core_depth %in% ignore)]
  #whichkeep <- whichkeep[!whichkeep %in% which(complete_core_depth %in% dt$depth_avg[is.na(dt$Pbex)])]
  rm(myvec)
  whichkeep <- whichkeep[!is.na(whichkeep)]
  # whichkeep tells you which "complete_core_depth" are not in an instantaneous deposit
  # Only these depths will be included when computing the inventory

  # Create the complete core depth 2 (for CRS calculations)
  # Everything but ignore
  complete_core_depth_2 <- complete_core_depth
  complete_core_depth_2[!complete_core_depth_2 %in% complete_core_depth_2[whichkeep]] <- NA

  # 1.6. Set some parameters ####
  if(!is.na(NWT) && Hemisphere=="NH") NWT_a <- 1963
  if(!is.na(NWT) && Hemisphere=="SH") NWT_a <- 1965
  if(is.null(dmin)) dmin <- min(dt$depth_avg,na.rm=T)
  if(!is.null(dmax) && length(inst_deposit)>1 && dmax<max(inst_deposit)) dmax <- max(inst_deposit)
  if(is.null(dmax)) if(exists("historic_d")) dmax <- max(c(dt$depth_avg,historic_d),na.rm=T) else dmax <- max(c(dt$depth_avg),na.rm=T)
  if(is.null(sedchange)) sedchange <- 0
  if(is.null(inst_deposit)) inst_deposit <- 0
  if(is.null(SML)) SML=0
  # Two next lines to prevent plotting the linear regression on raw 210Pb if instantaneous deposit are present
  if(is.null(plot_CFCS_regression) & plot_Pb_inst_deposit) plot_CFCS_regression=TRUE
  if(is.null(plot_CFCS_regression) & plot_Pb & length(inst_deposit) %% 2 != 1 & min(inst_deposit)<= max(dt$depth_avg)) plot_CFCS_regression=FALSE else plot_CFCS_regression=TRUE
  myylim <- c(-mround(dmax,10),-mround(dmin,10)+10)
  dates <- NULL
  dates_depth_avg <- NULL
  err_dated_depth_avg <- matrix(nrow=2)
  err_dates_avg <- NULL
  mtext_Cs <- NULL # text legend in case mass_depth=T and there are some peaks detected
  mylegend <- NULL
  mypchlegend <- NULL
  myltylegend <- NULL
  mycollegend <- NULL

  # 1.7 Create the composite free depth_avg ####
  ### Create the composite free depth_avg - step 1
  if(!exists("ignore")) ignore <- NULL
  if(SML>0) ignore <- c(ignore,dt$depth_avg[dt$depth_avg<=SML])

  # For Iseo for example, some NA were added? Quick fix.
  ignore <- ignore[!is.na(ignore)]
  # It is very circular, but if creating that made ignore == logical(0),
  #    then we just want it to be back at NULL. Very quick fix.
  if(length(ignore)==0) ignore <- NULL

  dt$depth_avg_2 <- rep(NA,nrow(dt))
  for (i in 1:nrow(dt)) {
    if(length(ignore)>0|max(inst_deposit)>0) {
      if(any(ignore==dt$depth_avg[i])) {
        dt$depth_avg_2[i] <- NA
      } else {dt$depth_avg_2[i] <- dt$depth_avg[i]}
    } else {dt$depth_avg_2[i] <- dt$depth_avg[i]}
  }

  ### Create the composite free depth_avg - step 2: inst_deposit
  if (length(sedchange)==1 && sedchange == 0) sedchange_corr=max(dt$depth_avg,na.rm = T) else sedchange_corr=sedchange

  if(exists("inst_deposit")&length(inst_deposit) > 1)
  {
    if(length(inst_deposit) %% 2 == 1) stop("\n Warning, inst_deposits need both upper and lower depths. Please check the manual.", call.=FALSE)
    inst_deposit_present = TRUE # Argument inst_deposit_present = FALSE decided elsewhere if no inst_deposit
    inst_deposit <- matrix(sort(inst_deposit), ncol=2, byrow=TRUE)
    for(i in 1:nrow(inst_deposit)) {
      if(length(dt$depth_avg[dt$depth_avg >= min(inst_deposit[i,]) & dt$depth_avg <= max(inst_deposit[i,])])>0) ignore <- c(ignore,dt$depth_avg[dt$depth_avg >= min(inst_deposit[i,]) & dt$depth_avg <= max(inst_deposit[i,])])
    }
    if(length(ignore)>0) ignore <- ignore[!duplicated(ignore)]
    if(length(ignore)>0) ignore <- ignore[order(ignore)]
    for (i in 1:nrow(dt)) {
      if((length(ignore)>0&&max(ignore, na.rm=T)>0)|(max(inst_deposit)>0)) {
        if(length(ignore)>0&&any(ignore==dt$depth_avg[i])) {
          dt$depth_avg_2[i] <- NA
        } else {dt$depth_avg_2[i] <- dt$depth_avg[i]}
      } else {dt$depth_avg_2[i] <- dt$depth_avg[i]}
    }
    d <- dt$depth_avg_2[!is.na(dt$depth_avg_2)]
    if(!is.na("historic_d")) dmax_corr=max(c(dt$depth_avg,historic_d),na.rm=T) else dmax_corr=max(c(dt$depth_avg),na.rm=T)
    inst_deposit_corr <- inst_deposit
    complete_core_depth_corr <- complete_core_depth[!is.na(complete_core_depth_2)]
    if(inst_deposit_present) for(i in 1:nrow(inst_deposit))
    {
      d[d > min(inst_deposit_corr[i,])] <- d[d > min(inst_deposit_corr[i,])] - (max(inst_deposit_corr[i,]) - min(inst_deposit_corr[i,]))
      # The depth for CRS should also be corrected for instantaneous events
      complete_core_depth_corr[complete_core_depth_corr > min(inst_deposit_corr[i,])] <- complete_core_depth_corr[complete_core_depth_corr > min(inst_deposit_corr[i,])] - (max(inst_deposit_corr[i,]) - min(inst_deposit_corr[i,]))

      dmax_corr <- dmax_corr - (max(inst_deposit_corr[i,])-min(inst_deposit_corr[i,]))
      for (n in seq_along(sedchange_corr)) {
        if(sedchange_corr[n] > min(inst_deposit_corr[i,],na.rm=T)) sedchange_corr[n] <- sedchange_corr[n][sedchange_corr[n] > min(inst_deposit_corr[i,])] - (max(inst_deposit_corr[i,]) - min(inst_deposit_corr[i,]))
      }
      if((1+i)<=nrow(inst_deposit)) inst_deposit_corr[c(1+i):nrow(inst_deposit),] <- inst_deposit_corr[c(1+i):nrow(inst_deposit),] - (max(inst_deposit_corr[i,])-min(inst_deposit_corr[i,]))
    }
    dt <- dt[order(dt$depth_avg_2),]
    dt$d <- c(d,rep(NA,length(dt$depth_avg)-length(d)))
    dt <- dt[order(dt$depth_avg),]
  } else {
    inst_deposit_present = FALSE
    dt$d=dt$depth_avg_2
    complete_core_depth_corr <- complete_core_depth[!is.na(complete_core_depth_2)]
  }

  # By the end here, you should have 3 columns for depth_avg: 1 with original depth_avg, 1 with removed events + suspicious data, 1 with event free depth_avg

  # 1.8 Mass_depth ####
  {
    # 1.8.1 mass_depth - Create the composite mass depth corrected for event ####
    if(mass_depth) {
      # Step 1.1: mass thickness (epaisseur massique)
      dt$mass_depth_top    <- rep(NA, nrow(dt))
      dt$mass_depth_bottom <- dt$density/10 * (dt$depth_bottom - dt$depth_top)
      for(i in 2:nrow(dt)) {
        if(dt$depth_top[i]==dt$depth_bottom[i-1]) {
          dt$mass_depth_top[i] = dt$mass_depth_bottom[i-1]
        } else {
          dt$mass_depth_top[i] = complete_core_density[which(complete_core_depth_bottom==dt$depth_top[i])] *
            (complete_core_depth_bottom[which(complete_core_depth_bottom==dt$depth_top[i])] - complete_core_depth_top[which(complete_core_depth_bottom==dt$depth_top[i])])
        }
      }
      dt$mass_depth_top[1]=0
      if(dt$depth_top[1]!=0) message(paste0("\n Warning. Mass depth for your first sample (",dt$depth_top[1],"-",dt$depth_bottom[1], " mm) was set\n to 0 to allow further calculation, but you did not provide\n density for the surface layer (0-",dt$depth_top[1]," mm). Include density for\n the surface layer if you can, or interpret the results with care.\n\n"))

      dt$mass_depth_avg         <- rep(NA, nrow(dt))

      # Step 2 (optional): if inst deposit calculate corrected mass_depth
      if(inst_deposit_present) { # If inst_deposit, create specific columns
        # Start from same vectors
        dt$mass_depth_top_corr    <- dt$mass_depth_top
        dt$mass_depth_bottom_corr <- dt$mass_depth_bottom

        # Put NA for the ones with missing values
        dt$mass_depth_top_corr[is.na(dt$depth_avg_2)]    <- NA
        dt$mass_depth_bottom_corr[is.na(dt$depth_avg_2)] <- NA
        dt <- dt[c(order(dt$depth_avg)[!is.na(dt$depth_avg_2)],which(is.na(dt$depth_avg_2))),]

        # Finally for the corrected vector, compute the actual mass depth (integral starting from the surface)
        for(i in 2:nrow(dt[which(!is.na(dt$depth_avg_2)),])) {
          dt$mass_depth_top_corr[i]    <- dt$mass_depth_top_corr[i-1]    + dt$mass_depth_top_corr[i]
          dt$mass_depth_bottom_corr[i] <- dt$mass_depth_bottom_corr[i-1] + dt$mass_depth_bottom_corr[i]
        }
        dt$mass_depth_avg_corr         <- (dt$mass_depth_bottom_corr + dt$mass_depth_top_corr)/2
        dt <- dt[order(dt$depth_avg),]

      }

      # Step 3: calculate actual mass depth, integral starting from the surface
      for(i in 2:nrow(dt)) {
        dt$mass_depth_top[i]    <- dt$mass_depth_top[i-1]    + dt$mass_depth_top[i]
        dt$mass_depth_bottom[i] <- dt$mass_depth_bottom[i-1] + dt$mass_depth_bottom[i]
      }
      dt$mass_depth_avg         <- (dt$mass_depth_bottom + dt$mass_depth_top)/2

      # Step 4 (optional, only if step 2 was skipped)
      # If no inst deposit, we still need these vector for plotting
      # Since there's no need to correct, create the corrected as the exact same than the normal
      if(!inst_deposit_present) {
        dt$mass_depth_top_corr    <- dt$mass_depth_top
        dt$mass_depth_bottom_corr <- dt$mass_depth_bottom
        dt$mass_depth_avg_corr    <- dt$mass_depth_avg
      }

      # Lastly, ylim for mass depth plots
      myylim_md <- c(-ceiling(max(c(dt$mass_depth_bottom,dt$mass_depth_top), na.rm=T)),0)
    }
    # By the end here, you should have 6 columns for mass_depth:
    #    2 times 3 columns. The 3 columns are top, average, and bottom mass depth
    #    replicate because there's actual depth (for plot_Pb)
    #    and depth with inst_deposit (for plot_Pb_inst_deposit)

    # 1.8.2 mass_depth - Create an extra column 'which_scale' for depth, according to mass_depth==T/F ####
    if(!mass_depth) dt$which_scale <- dt$d else dt$which_scale <- dt$mass_depth_avg_corr
    # 1.8.3 mass_depth - Create an interpolated mass_depth vector ####
    # If mass_depth=T, we'll need to match depths in g/cm2 to depths in mm.
    # CFCS ages between two depths are easy to find (linear relationship)
    # For mass_depth, if the interval is too big, we can really lose a lot
    #      of info.
    # This temporary vector will just interpolate mass_depth between two
    #      calculated values.
    # It's an apporximation; the higher the resolution of input data,
    #      the better (less approximation are needed)
    step_out_md   <- 1
    md_interp     <- c(seq(min(c(dt$depth_top,dt$depth_bottom), na.rm=T), max(c(dt$depth_top,dt$depth_bottom), na.rm=T), by = step_out_md))
    md_interp     <- matrix(c(md_interp,rep(NA,length(md_interp)*3)), ncol = 4)
    md_interp[,2] <- approx(x= dt$depth_top, dt$mass_depth_top, xout= md_interp[,1], rule=2)$y
    md_interp[,3] <- approx(x= dt$depth_bottom, dt$mass_depth_top, xout= md_interp[,1], rule=2)$y
    md_interp[,4] <- approx(x= dt$depth_avg, dt$mass_depth_top, xout= md_interp[,1], rule=2)$y
    md_interp     <- as.data.frame(md_interp)
    colnames(md_interp) <- c("depth_mm", "md_top", "md_bott", "md_avg")
  }
  # 1.9. If scale was given in mass_depth, convert in mm so the rest of the script work####
  if(input_depth_mm==F) {
    # User gave the input in g/cm2. Converting it in mm, because the script was initially built that way.
    # message if conversion
    msg_conversion <- " (depth converted from g/cm2)\n     "
    # SML
    if(SML!=0) SML <- md_interp$depth_mm[which.min(abs(md_interp$md_avg - SML))]
    # Cher
    if(!is.null(Cher)&&!is.na(Cher)) for(i in seq_along(Cher)) Cher[i] <- md_interp$depth_mm[which.min(abs(md_interp$md_avg - Cher[i]))]
    # NWT
    if(!is.null(NWT)&&!is.na(NWT))  for(i in seq_along(NWT))  NWT[i] <- md_interp$depth_mm[which.min(abs(md_interp$md_avg - NWT[i]))]
    # FF
    if(!is.null(FF)&&!is.na(FF))  for(i in seq_along(FF))    FF[i] <- md_interp$depth_mm[which.min(abs(md_interp$md_avg - FF[i]))]
    # inst_deposit
    if(!is.null(inst_deposit)&&!is.na(inst_deposit))  for(i in 1:nrow(inst_deposit))    inst_deposit[i,1] <- md_interp$depth_mm[which.min(abs(md_interp$md_avg - inst_deposit[i,1]))]
    if(!is.null(inst_deposit)&&!is.na(inst_deposit))  for(i in 1:nrow(inst_deposit))    inst_deposit[i,2] <- md_interp$depth_mm[which.min(abs(md_interp$md_avg - inst_deposit[i,2]))]
    # inst_deposit_corr
    if(!is.null(inst_deposit_corr)&&!is.na(inst_deposit_corr))  for(i in 1:nrow(inst_deposit_corr))    inst_deposit_corr[i,1] <- md_interp$depth_mm[which.min(abs(md_interp$md_avg - inst_deposit_corr[i,1]))]
    if(!is.null(inst_deposit_corr)&&!is.na(inst_deposit_corr))  for(i in 1:nrow(inst_deposit_corr))    inst_deposit_corr[i,2] <- md_interp$depth_mm[which.min(abs(md_interp$md_avg - inst_deposit_corr[i,2]))]

    # ignore
    if(!is.null(ignore)&&!is.na(ignore))  for(i in seq_along(ignore))    ignore[i] <- md_interp$depth_mm[which.min(abs(md_interp$md_avg - ignore[i]))]
  } else msg_conversion<-NULL
  # 1.10. Prepare depth vectors for CFCS model ####
  # Step somehow related to the creation of the composite free depht_avg
  # Here, we are looking to get two vectors:
  #     - (a) One vector of the actual depths on the core we are trying to date (upper and lower limits of instantaneous deposit for instance)
  #     - (b) The corrected version of this 1st vector, with instantaneous deposit removed
  #     - (c) Depth vector that will be used to visualise CFCS model
  d_for_CFCS <- unique(c(inst_deposit,max(dt$depth_avg[!is.na(dt$d)])))
  if(SML!=0) d_for_CFCS <- c(d_for_CFCS,SML)
  if(mass_depth) {
    d_for_CFCS <- c(d_for_CFCS,dt$depth_bottom[!is.na(dt$mass_depth_bottom)])
    d_for_CFCS <- unique(d_for_CFCS)
  }
  d_for_CFCS <- c(0,d_for_CFCS,sedchange,dmax)
  d_for_CFCS <- d_for_CFCS[order(d_for_CFCS)]

  # (a) Final vector of depths that will be used for linear model (This is the 1/2 vector we're creating in this section)
  depth_avg_to_date <- d_for_CFCS

  # The next steps are for the corrected vector (b)
  # When there is an instantaneous deposit, the top depth is used.
  if(max(inst_deposit,na.rm = T)>0) {
    for(r in 1:nrow(inst_deposit)) {
      d_for_CFCS[d_for_CFCS > min(inst_deposit[r,]) & d_for_CFCS <= max(inst_deposit[r,])]  <- d_for_CFCS[d_for_CFCS == min(inst_deposit[r,])]
    }
  }

  # (c) d_for_CFCS will be the vector used for visualization
  # (b) depth_avg_to_date_corr will be the totally corrected vector, and is created in the next steps.
  depth_avg_to_date_corr <- d_for_CFCS[order(d_for_CFCS)]

  # If there are more than 1 change in sedimentation rate,we also need to
  #   substract the sediment layer to the corrected depth average to date
  #   vector 'depth_avg_to_date_corr'
  # Create a corrected table of inst_deposit that take in account instantaneous deposits
  inst_deposit_corr2 <- inst_deposit
  if(exists("inst_deposit")&&length(inst_deposit) > 1)
    for(i in 1:nrow(inst_deposit)) {
      depth_avg_to_date_corr[depth_avg_to_date_corr>min(inst_deposit_corr2[i,])] <-
        depth_avg_to_date_corr[depth_avg_to_date_corr>min(inst_deposit_corr2[i,])] - (max(inst_deposit_corr2[i,])-min(inst_deposit_corr2[i,]))
      if((1+i)<=nrow(inst_deposit))  {
        inst_deposit_corr2[c(1+i):nrow(inst_deposit),] <-
          inst_deposit_corr2[c(1+i):nrow(inst_deposit),] - (max(inst_deposit_corr[i,])-min(inst_deposit_corr[i,]))
      }
    }
  rm(inst_deposit_corr2)

  # 1.11 Create separate datasets for different sedimentation rates ####
  if(length(sedchange)==1 && sedchange==0) {dt_sed1=dt} else {
    if(length(sedchange)==1) {
      dt_sed1 <- dt[dt$depth_avg<sedchange,]
      dt_sed2 <- dt[dt$depth_avg>=sedchange,]
    } else {
      dt_sed1 <- dt[dt$depth_avg<=min(sedchange),]
      dt_sed2 <- dt[dt$depth_avg>=min(sedchange) & dt$depth_avg<=max(sedchange),]
      dt_sed3 <- dt[dt$depth_avg>=max(sedchange),]
    }
  }

  # 1.12. Save data to the output list ####
  out_list$data <- dt[-grep("which_scale",colnames(dt))]
  if(suppdescriptor) out_list$data_suppdescriptor <- dt_suppdescriptor
  if(varves) out_list$data_varves <- varve

  # 1.13. Save the code to the output file with the code history ####
  # save the model attempt in a file
  # Row with all parameters that will be incremented:
  this_code_history <- c(name,coring_yr,as.character(Sys.time()),paste(model, collapse = ", "),
                         paste(Cher, collapse = ", "), paste(NWT, collapse = ", "), Hemisphere, paste(FF, collapse = ", "),
                         paste(inst_deposit, collapse = ", "),
                         paste(ignore, collapse = ", "),
                         paste(historic_d, collapse = ", "),
                         paste(historic_a, collapse = ", "),
                         paste(historic_n, collapse = ", "),
                         suppdescriptor,
                         paste(descriptor_lab, collapse = ", "),
                         varves,
                         paste(sedchange, collapse = ", "),
                         SML)
  this_code_history[this_code_history==""]=NA
  this_code_history<- as.data.frame(matrix(this_code_history, nrow=1))
  colnames(this_code_history) <- c("name","coring_yr","date_computation", "model_tested",
                                   "Chernobyl","NWT","Hemisphere","FF",
                                   "inst_deposit","ignore_depths",
                                   "historic_depth","historic_age","historic_name",
                                   "suppdescriptor","descriptor_lab",
                                   "varves","sedchange","SML")
  # First, check whether a file already exists
  if(length(list.files(paste(getwd(),"/Cores/",name,"/", sep=""), pattern="serac_model_history*", full.names=TRUE))==1) {
    # read previous file
    code_history <- read.delim(list.files(paste(getwd(),"/Cores/",name,"/", sep=""), pattern="serac_model_history*", full.names=TRUE))
    # increment new code
    if(ncol(code_history) == ncol(this_code_history))
      code_history <- rbind(code_history,this_code_history) else
        code_history <- this_code_history
  } else {
    code_history <- this_code_history
  }
  colnames(code_history) <- colnames(this_code_history)
  write.table(x = code_history, file = paste0(getwd(),"/Cores/",name,"/serac_model_history_",name,".txt"),col.names = T, row.names = F, sep = "\t")
  #Check whether the code is a duplicate from a previous code (has this combination been tested before)
  # First reread the file so all have the same format (some columns are turned in logical argument)
  code_history <- read.delim(list.files(paste(getwd(),"/Cores/",name,"/", sep=""), pattern="serac_model_history*", full.names=TRUE))
  if(nrow(code_history)>1 && all(sapply(code_history, function(x) duplicated(x))[nrow(code_history),-3])==TRUE)
    cat(paste0("\n General message: It seems you already tried this code combination. \n A historic of parameters tested can be looked up in the file\n 'serac_model_history_", name,".txt' (in the core directory).\n ________________________\n \n\n"))

  #### 2. LEAD 210 MODEL -----
  if(length(grep("Pb",x = colnames(dt)))>1 & length(grep("density",x = colnames(dt)))>=1) {
    # Inventory = sum(activity layer z * dry sediment accumulated at layer z * thickness layer z)
    # The inventory should account only for the continuous deposition:
    # [whichkeep] allows to keep only the data for the depth that are not in an instantaneous deposit
    Inventory_CRS <- complete_core_Pbex[whichkeep]*complete_core_density[whichkeep]*complete_core_thickness[whichkeep]
    Inventory_CRS_error <- complete_core_Pbex_err[whichkeep]*complete_core_density[whichkeep]*complete_core_thickness[whichkeep]
    Inventory_CRS_error[is.na(Inventory_CRS_error)] <- 0
    # Inventory: sum from depth to the bottom
    for(i in 1:length(Inventory_CRS)) {
      Inventory_CRS[i] <- sum(Inventory_CRS[i:length(Inventory_CRS)])
      Inventory_CRS_error[i] <- sum(Inventory_CRS_error[i:length(Inventory_CRS_error)])
    }
  }

  if(length(model)>=1) {
    # Write the result of the model
    # Calculate sedimentation rate
    # Le 0.031 dans l'équation du taux de sédimentation est issus de la periode de demi-vie du 210Pb qui est de 22.3 ans (T) donc il faut calculer le lambda qui est la constante de désintégration lambda = ln(2)/T=0.031.
    # Puis cf l'équation de désintégration du Plomb : plomb excess(t) = plomb excess (0) * e^(lambda*t)
    # Par analogie la pente = -lambda/V, V étant le taux de sédimentation
    lambda = log(2)/22.3
    lambda_err = 0.00017

    if(any(model=="CFCS")) {
      # Linear model and V calculation (sedimentation rate)
      lm_sed1 <- lm(log(dt_sed1$Pbex[!is.na(dt_sed1$d)&dt_sed1$Pbex>0]) ~ dt_sed1$which_scale[!is.na(dt_sed1$d)&dt_sed1$Pbex>0])
      sr_sed1 <- lambda/lm_sed1$coefficients[2]
      sr_sed1_err = sr_sed1*((lambda_err/lambda)^2+(summary(lm_sed1)$coefficients[2,2]/lm_sed1$coefficients[2])^2)^(0.5)

      # Save output
      if(!mass_depth) {
        out_list$`CFCS sediment accumulation rate` <- data.frame("SAR_mm.yr-1"=as.numeric(sr_sed1),"error_mm.yr-1"=as.numeric(sr_sed1_err),"R2"=summary(lm_sed1)$r.squared)
        rownames(out_list$`CFCS sediment accumulation rate`) <- "sedchange1"
      } else {
        out_list$`CFCS mass accumulation rate` <- data.frame("MAR_g.mm-1.yr-1"=as.numeric(sr_sed1),"error_g.mm-1.yr-1"=as.numeric(sr_sed1_err),"R2"=summary(lm_sed1)$r.squared)
        rownames(out_list$`CFCS mass accumulation rate`) <- "sedchange1"
      }

      # Print sed rate and error
      if (max(sedchange)==0) {
        if(!mass_depth) {
          cat(paste("\n Sediment accumulation rate (CFCS model): SAR= ",abs(round(sr_sed1,3))," mm/yr, R2= ", round(summary(lm_sed1)$r.squared,4),"\n", sep=""))
          cat(paste("                          Error:     +/- ",abs(round(sr_sed1_err,3))," mm/yr\n", sep=""))
        } else {
          cat(paste("\n Mass accumulation rate (CFCS model): MAR= ",abs(round(sr_sed1,3))," g/mm/yr, R2= ", round(summary(lm_sed1)$r.squared,4),"\n", sep=""))
          cat(paste("                              Error:     +/- ",abs(round(sr_sed1_err,3))," g/mm/yr\n", sep=""))
        }
      }

      if (max(sedchange)>0) {
        if(length(sedchange)==1) {
          # Linear model and V calculation (sedimentation rate)
          lm_sed2 <- lm(log(dt_sed2$Pbex[!is.na(dt_sed2$d)&dt_sed2$Pbex>0]) ~ dt_sed2$which_scale[!is.na(dt_sed2$d)&dt_sed2$Pbex>0])
          sr_sed2 <- lambda/lm_sed2$coefficients[2]
          sr_sed2_err = sr_sed2*((lambda_err/lambda)^2+(summary(lm_sed2)$coefficients[2,2]/lm_sed2$coefficients[2])^2)^(0.5)

          # Print sed rate and error
          if(!mass_depth) {
            cat(paste("\n Sediment accumulation rate (CFCS model) ", SML,"-",sedchange[1],"mm: SAR= ",abs(round(sr_sed1,3))," mm/yr, R2= ", round(summary(lm_sed1)$r.squared,4),"\n", sep=""))
            cat(paste("                          Error:     +/- ",abs(round(sr_sed1_err,3))," mm/yr\n", sep=""))
            cat(paste("\n Sediment accumulation rate (CFCS model) ", sedchange[1],"mm-bottom",": SAR= ",abs(round(sr_sed2,3))," mm/yr, R2= ", round(summary(lm_sed2)$r.squared,4),"\n", sep=""))
            cat(paste("                          Error:     +/- ",abs(round(sr_sed2_err,3))," mm/yr\n", sep=""))
          } else {
            cat(paste("\n Mass accumulation rate (CFCS model) ", SML,"-",sedchange[1],"mm: MAR= ",abs(round(sr_sed1,3))," g/mm/yr, R2= ", round(summary(lm_sed1)$r.squared,4),"\n", sep=""))
            cat(paste("                              Error:     +/- ",abs(round(sr_sed1_err,3))," g/mm/yr\n", sep=""))
            cat(paste("\n Mass accumulation rate (CFCS model) ", sedchange[1],"mm-bottom",": MAR= ",abs(round(sr_sed2,3))," g/mm/yr, R2= ", round(summary(lm_sed2)$r.squared,4),"\n", sep=""))
            cat(paste("                              Error:     +/- ",abs(round(sr_sed2_err,3))," g/mm/yr\n", sep=""))
          }

          # Save output
          if(!mass_depth) { # default
            out_list$`CFCS sediment accumulation rate` <- rbind(out_list$`CFCS sediment accumulation rate`,
                                                                c(as.numeric(sr_sed2),as.numeric(sr_sed2_err),summary(lm_sed2)$r.squared))
            rownames(out_list$`CFCS sediment accumulation rate`) <- c("sedchange1","sedchange2")
          } else {
            out_list$`CFCS mass accumulation rate` <- rbind(out_list$`CFCS mass accumulation rate`,
                                                            c(as.numeric(sr_sed2),as.numeric(sr_sed2_err),summary(lm_sed2)$r.squared))
            rownames(out_list$`CFCS mass accumulation rate`) <- c("sedchange1","sedchange2")
          }
        }
        if(length(sedchange)==2) {
          ## 2nd change in sedimentation rate
          # Linear model and V calculation (sedimentation rate)
          lm_sed2 <- lm(log(dt_sed2$Pbex[!is.na(dt_sed2$d)&dt_sed2$Pbex>0]) ~ dt_sed2$which_scale[!is.na(dt_sed2$d)&dt_sed2$Pbex>0])
          sr_sed2 <- lambda/lm_sed2$coefficients[2]
          sr_sed2_err = sr_sed2*((lambda_err/lambda)^2+(summary(lm_sed2)$coefficients[2,2]/lm_sed2$coefficients[2])^2)^(0.5)

          ## 3rd change in sedimentation rate
          # Linear model and V calculation (sedimentation rate)
          lm_sed3 <- lm(log(dt_sed3$Pbex[!is.na(dt_sed3$d)&dt_sed3$Pbex>0]) ~ dt_sed3$which_scale[!is.na(dt_sed3$d)&dt_sed3$Pbex>0])
          sr_sed3 <- lambda/lm_sed3$coefficients[2]
          sr_sed3_err = sr_sed3*((lambda_err/lambda)^2+(summary(lm_sed3)$coefficients[2,2]/lm_sed3$coefficients[2])^2)^(0.5)

          # Print sed rate and error
          if (!mass_depth) {
            cat(paste("\n Sediment accumulation rate (CFCS model) ", SML,"-",sedchange[1],"mm: SAR= ",abs(round(sr_sed1,3))," mm/yr, R2= ", round(summary(lm_sed1)$r.squared,4),"\n", sep=""))
            cat(paste("                          Error:     +/- ",abs(round(sr_sed1_err,3))," mm/yr\n", sep=""))
            cat(paste("\n Sediment accumulation rate (CFCS model) ", sedchange[1],"-",sedchange[2],"mm: SAR= ",abs(round(sr_sed2,3))," mm/yr, R2= ", round(summary(lm_sed2)$r.squared,4),"\n", sep=""))
            cat(paste("                          Error:     +/- ",abs(round(sr_sed2_err,3))," mm/yr\n", sep=""))
            cat(paste("\n Sediment accumulation rate (CFCS model) ", sedchange[2],"mm-bottom",": SAR= ",abs(round(sr_sed3,3))," mm/yr, R2= ", round(summary(lm_sed3)$r.squared,4),"\n", sep=""))
            cat(paste("                          Error:     +/- ",abs(round(sr_sed3_err,3))," mm/yr\n", sep=""))
          } else {
            cat(paste("\n Mass accumulation rate (CFCS model) ", SML,"-",sedchange[1],"mm: MAR= ",abs(round(sr_sed1,3))," g/mm/yr, R2= ", round(summary(lm_sed1)$r.squared,4),"\n", sep=""))
            cat(paste("                              Error:     +/- ",abs(round(sr_sed1_err,3))," g/mm/yr\n", sep=""))
            cat(paste("\n Mass accumulation rate (CFCS model) ", sedchange[1],"-",sedchange[2],"mm: MAR= ",abs(round(sr_sed2,3))," g/mm/yr, R2= ", round(summary(lm_sed2)$r.squared,4),"\n", sep=""))
            cat(paste("                              Error:     +/- ",abs(round(sr_sed2_err,3))," g/mm/yr\n", sep=""))
            cat(paste("\n Mass accumulation rate (CFCS model) ", sedchange[2],"mm-bottom",": MAR= ",abs(round(sr_sed3,3))," g/mm/yr, R2= ", round(summary(lm_sed3)$r.squared,4),"\n", sep=""))
            cat(paste("                              Error:     +/- ",abs(round(sr_sed3_err,3))," g/mm/yr\n", sep=""))
          }

          # Save output
          if(!mass_depth) {
            out_list$`CFCS sediment accumulation rate` <- rbind(out_list$`CFCS sediment accumulation rate`,
                                                                c(as.numeric(sr_sed2),as.numeric(sr_sed2_err),summary(lm_sed2)$r.squared))
            out_list$`CFCS sediment accumulation rate` <- rbind(out_list$`CFCS sediment accumulation rate`,
                                                                c(as.numeric(sr_sed3),as.numeric(sr_sed3_err),summary(lm_sed3)$r.squared))
            rownames(out_list$`CFCS sediment accumulation rate`) <- c("sedchange1","sedchange2","sedchange3")
          } else {
            out_list$`CFCS mass accumulation rate` <- rbind(out_list$`CFCS mass accumulation rate`,
                                                            c(as.numeric(sr_sed2),as.numeric(sr_sed2_err),summary(lm_sed2)$r.squared))
            out_list$`CFCS mass accumulation rate` <- rbind(out_list$`CFCS mass accumulation rate`,
                                                            c(as.numeric(sr_sed3),as.numeric(sr_sed3_err),summary(lm_sed3)$r.squared))
            rownames(out_list$`CFCS mass accumulation rate`) <- c("sedchange1","sedchange2","sedchange3")
          }
        }
      }
    }

    if(any(model=="CIC")) {
      # Calculation of age to be substracted
      Tm_CIC <- (1/lambda)*log(dt$Pbex[1]/dt$Pbex[!is.na(dt$Pbex)])
      # calculation age error: delta(tx)= 1/lambda*[(lambda_err*t)^2+(delta(A0)/A0)^2+(delta(Ax)/Ax)^2]^(0.5)
      # with Ax: activity at depth x; A0: initial activity
      # Two steps, 1 and 2
      # 1) replace error NA per 0 (just for this calculation, in temporary vectors)
      # We are not selecting the Pbex data with NAs, the model will interpolate between this ages
      # Having NA data breaks the assumptions of CIC model, so it's not
      #   recommended to use it in the first place if NAs are present.
      Pbex    <- dt$Pbex[!is.na(dt$Pbex)]
      Pbex_er <- dt$Pbex_er[!is.na(dt$Pbex)&!is.na(dt$Pbex_er)]
      Pbex_er <- approx(x = dt$depth_avg[!is.na(dt$Pbex)&!is.na(dt$Pbex_er)], y = Pbex_er,
                        xout = dt$depth_avg[!is.na(dt$Pbex)])$y
      # 2) Actual error
      Tm_CIC_err <- (1/lambda)*((lambda_err*Tm_CIC)^2+(Pbex_er[1]/Pbex[1])^2+(Pbex_er/Pbex)^2)^(0.5)

      # calculation of Best Age and errors
      m_CIC <- coring_yr - Tm_CIC
      m_CIC_low <- m_CIC - Tm_CIC_err
      m_CIC_high <- m_CIC + Tm_CIC_err

      # Save output
      out_list$`CIC model` <- data.frame("m_CIC"      = m_CIC,
                                         "m_CIC_low"  = m_CIC_low,
                                         "m_CIC_high" = m_CIC_high)
    }

    if(any(model=="CRS")) {
      if(rev(dt$Pbex)[1] >= dt$Pbex[1]/16) cat("\n Warning, it seems that 210Pb_excess has not reached equilibrium. \n Make sure the conditions of application for CRS model are fulfilled.")

      Tm_CRS <- 1/lambda*log(Inventory_CRS[1]/Inventory_CRS)
      # calculation age error: delta(tx)=1/lambda*((0.00017*t)^2+(delta(I0)/I0)^2+(1-2*Ix/Io)*(delta(Ix)/Ix)^2)^(0.5)
      # with I0: iInventory, Ix= Inventory below depth x
      Tm_CRS_err <- 1/lambda*((lambda_err*Tm_CRS)^2+(Inventory_CRS_error[1]/Inventory_CRS[1])^2+(1-2*Inventory_CRS/Inventory_CRS[1])*(Inventory_CRS_error/Inventory_CRS)^2)^(0.5)

      # calculation of Best Age and errors
      m_CRS <- coring_yr - Tm_CRS
      m_CRS_low <- m_CRS - Tm_CRS_err
      m_CRS_high <- m_CRS + Tm_CRS_err

      # Save output
      out_list$`CRS model` <- data.frame("m_CRS"      = m_CRS,
                                         "m_CRS_low"  = m_CRS_low,
                                         "m_CRS_high" = m_CRS_high)
    }
  }

  #### 3. CESIUM -----
  if(plot_Cs) {

    #Chernobyl
    if (all(!is.na(Cher))) {
      peakCher <- -mean(Cher)
      dates <- c(dates,1986)
      dates_depth_avg <- c(dates_depth_avg,peakCher)
      err_dated_depth_avg <- cbind(err_dated_depth_avg,Cher)
      err_dates_avg <- c(err_dates_avg,NA)
      if(mass_depth) mtext_Cs <- ("Cher")
    }
    #NWT
    if (all(!is.na(NWT))) {
      peakNWT <- -mean(NWT)
      dates <- c(dates,NWT_a)
      dates_depth_avg <- c(dates_depth_avg,peakNWT)
      err_dated_depth_avg <- cbind(err_dated_depth_avg,NWT)
      if (Hemisphere =="NH") err_dates_avg <- c(err_dates_avg,NA)
      if (Hemisphere =="SH") err_dates_avg <- c(err_dates_avg,.5)
      if(mass_depth)         mtext_Cs <- c(mtext_Cs, "NWT")
    }
    #First radionuclides fallout
    if (all(!is.na(FF))) {
      peakFF <- -mean(FF)
      dates <- c(dates,1955)
      dates_depth_avg <- c(dates_depth_avg,peakFF)
      err_dated_depth_avg <- cbind(err_dated_depth_avg,FF)
      err_dates_avg <- c(err_dates_avg,NA)
      if(mass_depth) mtext_Cs <- c(mtext_Cs, "FF")
    }


  }

  #### 4. AGE DEPTH MODEL -----
  # 4. Prep If mass_depth==T, convert here every depth in cm into depth g.cm2 ####
  if(mass_depth) {
    sedchange_corr_allscales <- NULL
    for(i in seq_along(sedchange_corr)) sedchange_corr_allscales <- c(sedchange_corr_allscales,md_interp$md_avg[which.min(abs(md_interp$depth_mm - sedchange_corr[i]))])
    sedchange_allscales      <- NULL
    for(i in seq_along(sedchange)) sedchange_allscales <- c(sedchange_allscales,md_interp$md_avg[which.min(abs(md_interp$depth_mm - sedchange[i]))])
    Cher_allscales           <- NULL
    for(i in seq_along(Cher)) Cher_allscales <- c(Cher_allscales,md_interp$md_avg[which.min(abs(md_interp$depth_mm - Cher[i]))])
    NWT_allscales            <- NULL
    for(i in seq_along(NWT)) NWT_allscales <- c(NWT_allscales,md_interp$md_avg[which.min(abs(md_interp$depth_mm - NWT[i]))])
    FF_allscales             <- NULL
    for(i in seq_along(FF)) FF_allscales <- c(FF_allscales,md_interp$md_avg[which.min(abs(md_interp$depth_mm - FF[i]))])
    depth_avg_to_date_allscales <- NULL
    for(i in seq_along(depth_avg_to_date)) c(depth_avg_to_date_allscales,depth_avg_to_date_allscales <- c(depth_avg_to_date_allscales,md_interp$md_avg[which.min(abs(md_interp$depth_mm - depth_avg_to_date[i]))]))
    depth_avg_to_date_allscales[1] <- 0
    depth_avg_to_date_corr_allscales <- NULL
    for(i in seq_along(depth_avg_to_date_corr)) depth_avg_to_date_corr_allscales <- c(depth_avg_to_date_corr_allscales,md_interp$md_avg[which.min(abs(md_interp$depth_mm - depth_avg_to_date_corr[i]))])
    depth_avg_to_date_corr_allscales[1] <- 0
    if (all(!is.na(Cher))) {
      peakCher_allscales <- NULL
      for(i in seq_along(peakCher)) peakCher_allscales <- c(peakCher_allscales,-md_interp$md_avg[which.min(abs(md_interp$depth_mm + peakCher[i]))])
    }
    if (all(!is.na(NWT))) {
      peakNWT_allscales <- NULL
      for(i in seq_along(peakNWT)) peakNWT_allscales <- c(peakNWT_allscales,-md_interp$md_avg[which.min(abs(md_interp$depth_mm + peakNWT[i]))])
    }
    if (all(!is.na(FF))) {
      peakFF_allscales <- NULL
      for(i in seq_along(peakFF)) peakFF_allscales <- c(peakFF_allscales,-md_interp$md_avg[which.min(abs(md_interp$depth_mm + peakFF[i]))])
    }
  } else {
    sedchange_corr_allscales <- sedchange_corr
    sedchange_allscales      <- sedchange
    Cher_allscales           <- Cher
    NWT_allscales            <- NWT
    FF_allscales             <- FF
    depth_avg_to_date_allscales <- depth_avg_to_date
    depth_avg_to_date_corr_allscales <- depth_avg_to_date_corr
    if (all(!is.na(Cher))) peakCher_allscales       <- peakCher
    if (all(!is.na(NWT))) peakNWT_allscales        <- peakNWT
    if (all(!is.na(FF))) peakFF_allscales         <- peakFF
  }

  # 4. Actual code for depth age model ####
  for (i in which(dt$depth_avg>=SML)){
    if (is.na(dt$d[i])) {
      if(i==1) dt$d[i] <- dt$depth_avg[1] else dt$d[i] <- dt$d[i-1]
    }
  }

  if(any(model=="CFCS")) {
    if(max(sedchange)>0) {
      age_break <- coring_yr-sedchange_corr_allscales[1]/abs(sr_sed1)
      age_break_low <- age_break-sedchange_corr_allscales[1]*abs(sr_sed1_err)/abs(sr_sed1)^2
      age_break_high <- age_break+sedchange_corr_allscales[1]*abs(sr_sed1_err)/abs(sr_sed1)^2
      cat(paste(" Approximation of age at change(s) in sedimentation rate:\n"))
      if(length(sedchange)==2) {
        age_break2 <- age_break-(sedchange_corr_allscales[2]-sedchange_corr_allscales[1])/abs(sr_sed2)
        age_break2_low <- age_break2-((sedchange_corr_allscales[2]-sedchange_corr_allscales[1]))*abs(sr_sed2_err)/abs(sr_sed2)^2
        age_break2_high <- age_break2+((sedchange_corr_allscales[2]-sedchange_corr_allscales[1]))*abs(sr_sed2_err)/abs(sr_sed2)^2
        cat(paste("     Best Age (1st change): ",abs(round(age_break,0))," (uncertainty: ",abs(round(age_break_low,0)),"-",abs(round(age_break_high,0)),")\n",sep=""))
        cat(paste("     Best Age (2nd change): ",abs(round(age_break2,0))," (uncertainty: ",abs(round(age_break2_low,0)),"-",abs(round(age_break2_high,0)),")\n\n",sep=""))
      } else {
        cat(paste("     Best Age: ",abs(round(age_break,0))," (uncertainty: ",abs(round(age_break_low,0)),"-",abs(round(age_break_high,0)),")\n\n",sep=""))
      }
    }

    output_agemodel_CFCS <- matrix(rep(NA,length(depth_avg_to_date_allscales)*4), ncol=4)
    for(i in seq_along(depth_avg_to_date_allscales)){
      output_agemodel_CFCS[i,1] <- depth_avg_to_date[i]
      output_agemodel_CFCS[i,2] <- coring_yr-depth_avg_to_date_corr_allscales[i]/abs(sr_sed1)
      output_agemodel_CFCS[i,3] <- output_agemodel_CFCS[i,2]-depth_avg_to_date_corr_allscales[i]*abs(sr_sed1_err)/abs(sr_sed1)^2
      output_agemodel_CFCS[i,4] <- output_agemodel_CFCS[i,2]+depth_avg_to_date_corr_allscales[i]*abs(sr_sed1_err)/abs(sr_sed1)^2

      if(max(sedchange)>0 && depth_avg_to_date[i]>=sedchange[1]) {
        output_agemodel_CFCS[i,2] <- age_break-(depth_avg_to_date_corr_allscales[i]-sedchange_corr_allscales[1])/abs(sr_sed2)
        output_agemodel_CFCS[i,3] <- output_agemodel_CFCS[i,2]-(depth_avg_to_date_corr_allscales[i])*abs(sr_sed2_err)/abs(sr_sed2)^2
        output_agemodel_CFCS[i,4] <- output_agemodel_CFCS[i,2]+(depth_avg_to_date_corr_allscales[i])*abs(sr_sed2_err)/abs(sr_sed2)^2
      }
      if(length(sedchange)>1 && depth_avg_to_date[i]>=sedchange[2]) {
        output_agemodel_CFCS[i,2] <- age_break2-(depth_avg_to_date_corr_allscales[i]-sedchange_corr_allscales[2])/abs(sr_sed3)
        output_agemodel_CFCS[i,3] <- output_agemodel_CFCS[i,2]-(depth_avg_to_date_corr_allscales[i])*abs(sr_sed3_err)/abs(sr_sed3)^2
        output_agemodel_CFCS[i,4] <- output_agemodel_CFCS[i,2]+(depth_avg_to_date_corr_allscales[i])*abs(sr_sed3_err)/abs(sr_sed3)^2
      }
    }
    output_agemodel_CFCS <- as.data.frame(output_agemodel_CFCS)
    colnames(output_agemodel_CFCS) <- c("depth", "BestAD", "MinAD", "MaxAD")
    if(mass_depth) output_agemodel_CFCS$mass_depth <- depth_avg_to_date_allscales
    output_agemodel_CFCS <- output_agemodel_CFCS[!duplicated(output_agemodel_CFCS[,1]),]

    # if mass depth, we do not want to interpolate below the last depth with measurement
    if(mass_depth) output_agemodel_CFCS <- output_agemodel_CFCS[output_agemodel_CFCS$depth<=max(dt$depth_bottom[!is.na(dt$Pbex)],na.rm=T),]

    output_agemodel_CFCS_inter <- as.data.frame(seq(0,max(output_agemodel_CFCS$depth,na.rm = T),by=stepout))
    if (length(historic_d)>=1 && !all(is.na(historic_d)) && any(is.na(historic_a))) {
      whichNA <- which(is.na(historic_a))
      historic_d_dt <- matrix(historic_d, ncol = 2, byrow = T)
      myage_low <- approx(x= output_agemodel_CFCS$depth, output_agemodel_CFCS$MinAD, xout= historic_d_dt[is.na(historic_a),2])$y
      myage_high <- approx(x= output_agemodel_CFCS$depth, output_agemodel_CFCS$MaxAD, xout= historic_d_dt[is.na(historic_a),1])$y
      cat(paste("\n Age approximation of non-dated historical events from CFCS model:\n"))
      for (i in whichNA) {
        cat(paste("     The historical event at ",historic_d_dt[whichNA,1][i],"-", historic_d_dt[whichNA,2][i]," mm has an estimated range of: ",round(myage_low[i]),"-",round(myage_high[i]),".\n",sep=""))
      }
    }
    if(exists("inst_deposit")&&length(inst_deposit)>1) {
      cat(paste("\n Age approximation of instantaneous deposit(s) from CFCS model:\n"))
      for (i in 1:nrow(inst_deposit)) {
        if(length(round(output_agemodel_CFCS[which(output_agemodel_CFCS[,1]==inst_deposit[i,1]),3]))>0)
          cat(paste0("     The instantaneous deposit at ",inst_deposit[i,1],"-", inst_deposit[i,2]," mm",msg_conversion," has an estimated range of: ", round(output_agemodel_CFCS[which(output_agemodel_CFCS[,1]==inst_deposit[i,1]),3]),"-",round(output_agemodel_CFCS[which(output_agemodel_CFCS[,1]==inst_deposit[i,1]),4]),".\n")) else
            cat(paste0("    /!\\  Instantaneous deposit at ",inst_deposit[i,1],"-", inst_deposit[i,2]," mm: not enough data to calculate an age-range.  "))
        if(i==nrow(inst_deposit)) cat("\n")
      }
    }

    output_agemodel_CFCS_inter <- as.data.frame(output_agemodel_CFCS_inter)

    # Interpolate to get the age-depth model with the input stepout
    # We were extra-cautious and first interpolated to a 0.1 mm resolution to be sure we wouldn't miss a change in sedimentation rate.
    temporary <- approx(x= output_agemodel_CFCS$depth, output_agemodel_CFCS$BestAD, xout= seq(0,max(output_agemodel_CFCS$depth,na.rm = T),.1))
    output_agemodel_CFCS_inter <- cbind(output_agemodel_CFCS_inter,approx(x= temporary$x, temporary$y, xout= seq(0,max(output_agemodel_CFCS$depth,na.rm = T),stepout))$y)
    temporary <- approx(x= output_agemodel_CFCS$depth, output_agemodel_CFCS$MinAD, xout= seq(0,max(output_agemodel_CFCS$depth,na.rm = T),.1))
    output_agemodel_CFCS_inter <- cbind(output_agemodel_CFCS_inter,approx(x= temporary$x, temporary$y, xout= seq(0,max(output_agemodel_CFCS$depth,na.rm = T),stepout))$y)
    temporary <- approx(x= output_agemodel_CFCS$depth, output_agemodel_CFCS$MaxAD, xout= seq(0,max(output_agemodel_CFCS$depth,na.rm = T),.1))
    output_agemodel_CFCS_inter <- cbind(output_agemodel_CFCS_inter,approx(x= temporary$x, temporary$y, xout= seq(0,max(output_agemodel_CFCS$depth,na.rm = T),stepout))$y)

    colnames(output_agemodel_CFCS_inter) <- c("depth_avg", "BestAD", "MinAD", "MaxAD")
    write.table(x = output_agemodel_CFCS[order(output_agemodel_CFCS$depth, decreasing = F),], file = paste(getwd(),"/Cores/",name,"/",name,"_CFCS.txt",sep = ""),col.names = T, row.names = F)
    write.table(x = output_agemodel_CFCS_inter[order(output_agemodel_CFCS_inter$depth_avg, decreasing = F),], file = paste(getwd(),"/Cores/",name,"/",name,"_CFCS_interpolation.txt",sep = ""),col.names = T, row.names = F)

    # Save output in the list
    out_list$'CFCS age-depth model interpolated' <- output_agemodel_CFCS_inter

    # Parameters for legend
    mylegend <- c(mylegend, "CFCS")
    mypchlegend <- c(mypchlegend,NA)
    myltylegend <- c(myltylegend,1)
    mycollegend <- c(mycollegend,modelcol[1])
  }

  if(any(model=="CIC")) {
    output_agemodel_CIC <- as.data.frame(matrix(c(c(0,dt$depth_avg[!is.na(dt$Pbex)]),c(coring_yr,m_CIC),c(coring_yr,m_CIC_low),c(coring_yr,m_CIC_high)), byrow = F, ncol=4))
    colnames(output_agemodel_CIC) <- c("depth_avg", "BestAD_CIC", "MinAD_CIC", "MaxAD_CIC")
    output_agemodel_CIC_inter <- as.data.frame(seq(0,max(output_agemodel_CIC$depth_avg,na.rm = T),stepout))
    output_agemodel_CIC_inter <- cbind(output_agemodel_CIC_inter,approx(x= output_agemodel_CIC$depth_avg, output_agemodel_CIC$BestAD_CIC, xout= seq(0,max(output_agemodel_CIC$depth_avg,na.rm = T),stepout))$y)
    output_agemodel_CIC_inter <- cbind(output_agemodel_CIC_inter,approx(x= output_agemodel_CIC$depth_avg, output_agemodel_CIC$MinAD_CIC, xout= seq(0,max(output_agemodel_CIC$depth_avg,na.rm = T),stepout))$y)
    output_agemodel_CIC_inter <- cbind(output_agemodel_CIC_inter,approx(x= output_agemodel_CIC$depth_avg, output_agemodel_CIC$MaxAD_CIC, xout= seq(0,max(output_agemodel_CIC$depth_avg,na.rm = T),stepout))$y)
    colnames(output_agemodel_CIC_inter) <- c("depth_avg", "BestAD_CIC", "MinAD_CIC", "MaxAD_CIC")
    write.table(x = output_agemodel_CIC[order(output_agemodel_CIC$depth_avg, decreasing = F),], file = paste(getwd(),"/Cores/",name,"/",name,"_CIC.txt",sep = ""),col.names = T, row.names = F)
    write.table(x = output_agemodel_CIC_inter[order(output_agemodel_CIC_inter$depth_avg, decreasing = F),], file = paste(getwd(),"/Cores/",name,"/",name,"_CIC_interpolation.txt",sep = ""),col.names = T, row.names = F)

    # Save output in the list
    out_list$'CIC age-depth model interpolated' <- output_agemodel_CIC_inter

    # Parameters for legend
    mylegend <- c(mylegend, "CIC")
    mypchlegend <- c(mypchlegend,NA)
    myltylegend <- c(myltylegend,1)
    mycollegend <- c(mycollegend,modelcol[2])
  }

  if(any(model=="CRS")) {
    if(exists("inst_deposit")&&length(inst_deposit)>1) {
      # Create the depth for CRS model plotting
      new_y_CRS <- c(complete_core_depth[!is.na(complete_core_depth_2)], as.vector(inst_deposit)[inst_deposit<=max(dt$depth_avg,na.rm=T)])
      new_y_CRS_corr <- c(complete_core_depth_corr, inst_deposit_corr[,1][inst_deposit[,1]<=max(dt$depth_avg,na.rm=T)], inst_deposit_corr[,1][inst_deposit[,2]<=max(dt$depth_avg,na.rm=T)])

      # Sort in correct order
      new_y_CRS <- new_y_CRS[order(new_y_CRS)]
      new_y_CRS_corr <- new_y_CRS_corr[order(new_y_CRS_corr)]

      # Create the new ages
      new_x_CRS <- approx(x = complete_core_depth_corr,
                          y = m_CRS,
                          xout = new_y_CRS_corr, rule = 2)$y
      new_x_CRS_low <- approx(x = complete_core_depth_corr,
                              y = m_CRS_low,
                              xout = new_y_CRS_corr, rule = 2)$y
      new_x_CRS_high <- approx(x = complete_core_depth_corr,
                               y = m_CRS_high,
                               xout = new_y_CRS_corr, rule = 2)$y
    }

    output_agemodel_CRS <- as.data.frame(matrix(c(0,complete_core_depth_top[order(complete_core_depth_top, decreasing = F)][whichkeep],c(coring_yr,m_CRS),c(coring_yr,m_CRS_low),c(coring_yr,m_CRS_high)), byrow = F, ncol=4))
    colnames(output_agemodel_CRS) <- c("depth", "BestAD_CRS", "MinAD_CRS", "MaxAD_CRS")
    output_agemodel_CRS_inter <- as.data.frame(seq(0,max(output_agemodel_CRS$depth,na.rm = T),stepout))
    output_agemodel_CRS_inter <- cbind(output_agemodel_CRS_inter,approx(x= output_agemodel_CRS$depth, output_agemodel_CRS$BestAD_CRS, xout= seq(0,max(output_agemodel_CRS$depth,na.rm = T),stepout))$y)
    output_agemodel_CRS_inter <- cbind(output_agemodel_CRS_inter,approx(x= output_agemodel_CRS$depth, output_agemodel_CRS$MinAD_CRS, xout= seq(0,max(output_agemodel_CRS$depth,na.rm = T),stepout))$y)
    output_agemodel_CRS_inter <- cbind(output_agemodel_CRS_inter,approx(x= output_agemodel_CRS$depth, output_agemodel_CRS$MaxAD_CRS, xout= seq(0,max(output_agemodel_CRS$depth,na.rm = T),stepout))$y)
    colnames(output_agemodel_CRS_inter) <- c("depth", "BestAD_CRS", "MinAD_CRS", "MaxAD_CRS")
    write.table(x = output_agemodel_CRS[order(output_agemodel_CRS$depth, decreasing = F),], file = paste(getwd(),"/Cores/",name,"/",name,"_CRS.txt",sep = ""),col.names = T, row.names = F)
    write.table(x = output_agemodel_CRS_inter[order(output_agemodel_CRS_inter$depth, decreasing = F),], file = paste(getwd(),"/Cores/",name,"/",name,"_CRS_interpolation.txt",sep = ""),col.names = T, row.names = F)

    # Save output in the list
    out_list$'CRS age-depth model interpolated' <- output_agemodel_CRS_inter

    # Parameters for legend
    mylegend <- c(mylegend, "CRS")
    mypchlegend <- c(mypchlegend,NA)
    myltylegend <- c(myltylegend,1)
    mycollegend <- c(mycollegend,modelcol[3])
  }

  if(varves) {
    # Parameters for legend
    mylegend <- c(mylegend, "varves")
    mypchlegend <- c(mypchlegend,4)
    myltylegend <- c(myltylegend,NA)
    mycollegend <- c(mycollegend,"black")
  }


  if(exists("Cher") | exists("NWT") | exists("FF")) {
    if(mass_depth) {
      err_dated_depth_avg_allscales <- NULL
      for(i in seq_along(err_dated_depth_avg)) err_dated_depth_avg_allscales <- c(err_dated_depth_avg_allscales,md_interp$md_avg[which.min(abs(md_interp$depth_mm - err_dated_depth_avg[i]))])
      err_dated_depth_avg_allscales <- matrix(err_dated_depth_avg_allscales[!is.na(err_dated_depth_avg_allscales)],nrow=2,byrow = F)
    }
    err_dated_depth_avg <- matrix(-err_dated_depth_avg[!is.na(err_dated_depth_avg)],nrow=2,byrow = F)
  }

  # Various parameter to print (e.g. Inventories)
  # Inventory Lead
  if(length(grep("Pb",x = colnames(dt)))>1 & length(grep("density",x = colnames(dt)))>=1) {
    # We multiply the value by 10 because we ask for the depth in mm, and the density in g/cm3
    cat(paste(" Inventory (Lead): ",round(Inventory_CRS[1],3), " Bq/m2 (range: ",round(Inventory_CRS[1]-Inventory_CRS_error[1]), "-",round(Inventory_CRS[1]+Inventory_CRS_error[1])," Bq/m2)\n", sep=""))

    # Save output in list
    out_list$`Inventories` <- rbind(out_list$`Inventories`,
                                    data.frame("Inventory" = Inventory_CRS[1],
                                               "min" = Inventory_CRS[1]-Inventory_CRS_error[1],
                                               "max" = Inventory_CRS[1]+Inventory_CRS_error[1]))
    rownames(out_list$`Inventories`)[nrow(out_list$`Inventories`)] <- "Lead   Bq.m-2"
  }

  # Inventory Cesium
  if(length(grep("Cs",x = colnames(dt)))>1 & length(grep("density",x = colnames(dt)))>=1) {
    # Inventory = sum(activity layer z * dry sediment accumuated at layer z * thickness layer z)
    # The inventory should account only for the continuous deposition:
    # [whichkeep] allows to keep only the data for the depth that are not in an instantaneous deposit
    Inventory_Cesium <- complete_core_Cs[whichkeep]*complete_core_density[whichkeep]*complete_core_thickness[whichkeep]
    Inventory_Cesium_low <- (complete_core_Cs-complete_core_Cs_err)[whichkeep]*complete_core_density[whichkeep]*complete_core_thickness[whichkeep]
    Inventory_Cesium_low[is.na(Inventory_Cesium_low)] <- Inventory_Cesium[is.na(Inventory_Cesium_low)]
    Inventory_Cesium_high <- (complete_core_Cs+complete_core_Cs_err)[whichkeep]*complete_core_density[whichkeep]*complete_core_thickness[whichkeep]
    Inventory_Cesium_high[is.na(Inventory_Cesium_high)] <- Inventory_Cesium[is.na(Inventory_Cesium_high)]
    # We multiply the value by 10 because we ask for the depth in mm, and the density in g/cm3
    cat(paste(" Inventory (Cesium): ",round(sum(Inventory_Cesium,na.rm=T),3), " Bq/m2 (range: ",round(sum(Inventory_Cesium_low,na.rm=T)), "-",round(sum(Inventory_Cesium_high,na.rm=T))," Bq/m2)\n", sep=""))

    # Save output in list
    out_list$`Inventories` <- rbind(out_list$`Inventories`,
                                    data.frame("Inventory" = sum(Inventory_Cesium,na.rm=T),
                                               "min" = sum(Inventory_Cesium_low,na.rm=T),
                                               "max" = sum(Inventory_Cesium_high,na.rm=T)))
    rownames(out_list$`Inventories`)[nrow(out_list$`Inventories`)] <- "Cesium Bq.m-2"
  }

  #### 5. OUTPUT FILE METADATA -----
  # Read supp metadata if the file already exists
  if(length(list.files(paste(getwd(),"/Cores/",name,"/", sep=""), pattern="serac_metadata_suppmetadata*", full.names=TRUE))==1) suppmetadata <- read.delim(list.files(paste(getwd(),"/Cores/",name,"/", sep=""), pattern="serac_metadata_suppmetadata*", full.names=TRUE), header=T, sep="")

  if(length(list.files(paste(getwd(),"/Cores", sep=""), pattern="serac_metadata*", full.names=TRUE))==1) {
    mmetadata <- read.delim(list.files(paste(getwd(),"/Cores", sep=""), pattern="serac_metadata*", full.names=TRUE), header=T, sep="")
  } else {
    mmetadata <- matrix(rep(NA,30), ncol=3)
  }

  metadata <- matrix(c("GENERAL_INFORMATIONS","",
                       "Core_code",name,
                       "Sampled",coring_yr,
                       "","",
                       "AGE_MODEL_COMPUTATION", "",
                       "Computation_date", paste(Sys.Date()),
                       "User",paste(mmetadata[1,2]),
                       "Computer_LogName",Sys.getenv("LOGNAME"),
                       "Affiliation",paste(mmetadata[2,2]),
                       "ORCID",paste(mmetadata[3,2]),
                       "email",paste(mmetadata[4,2]),
                       "","",
                       "AGE_MODEL_PARAMETERS", "",
                       "Method_selected",paste(model, collapse = " "),
                       "Varves",paste(varves),
                       "Change_sed_rate", paste(sedchange, collapse = ","),
                       "Chernobyl", paste(Cher, collapse = "-"),
                       "Nuclear_War_Test", paste(paste(NWT, collapse = "-"), " (",Hemisphere,")", sep=""),
                       "First_Fallout", paste(FF, collapse = "-"),
                       "Surface_Mixed_Layer",paste(SML, collapse = ""),
                       "Instantaneous_deposit_up_and_low_limits", paste(inst_deposit, collapse = ", "),
                       "Ignore_up_and_low_limits",paste(inst_deposit, collapse = ", "),
                       "Supplementary_descriptor",paste(descriptor_lab, collapse = ", "),
                       "Historic_depth_up_and_low_limits",paste(historic_d, collapse = ", "),
                       "Historic_age",paste(historic_a, collapse = ", "),
                       "Historic_name",paste(historic_n, collapse = ", "),
                       "Step_out",paste(stepout, collapse = ""),
                       "","",
                       "AGE_MODEL_OUTPUT", ""),
                     ncol = 2, byrow = T)

  if(exists("suppmetadata")) {
    metadata <- rbind(c("INFORMATIONS_REGARDING_MEASUREMENTS",""),
                      as.matrix(suppmetadata),
                      c("",""),
                      metadata)
  }

  # Add in output the results of the sedimentation rate
  if(any(model=="CFCS")) {
    if(!mass_depth) {
      if (max(sedchange)>0) {
        metadata <- rbind(metadata,
                          c(paste("Sediment accumulation rate (CFCS model) ", SML,"-",sedchange[1],"mm",sep=""),paste("SAR= ",abs(round(sr_sed1,3)),"mm/yr, R2= ", round(summary(lm_sed1)$r.squared,4),", Error +/- ",abs(round(sr_sed1_err,3)),"mm/yr", sep="")))
      } else {
        metadata <- rbind(metadata,
                          c("Sediment accumulation rate (CFCS model)",paste("SAR= ",abs(round(sr_sed1,3)),"mm/yr, R2= ", round(summary(lm_sed1)$r.squared,4),", Error +/- ",abs(round(sr_sed1_err,3)),"mm/yr", sep="")))
      }
      if (max(sedchange)>0) {
        if(length(sedchange)==1) {
          metadata <- rbind(metadata,
                            c(paste("Sediment accumulation rate (CFCS model) ", sedchange[1],"mm-bottom", sep=""),paste("SAR= ",abs(round(sr_sed2,3)),"mm/yr, R2= ", round(summary(lm_sed2)$r.squared,4), ", Error +/- ",abs(round(sr_sed1_err,3)),"mm/yr",sep="")))
        }
        if(length(sedchange)==2) {
          metadata <- rbind(metadata,
                            c(paste("Sediment accumulation rate (CFCS model) ", sedchange[1],"-",sedchange[2],"mm", sep=""),paste("SAR= ",abs(round(sr_sed2,3)),"mm/yr, R2= ", round(summary(lm_sed2)$r.squared,4), ", Error +/- ",abs(round(sr_sed1_err,3)),"mm/yr",sep="")),
                            c(paste("Sediment accumulation rate (CFCS model) ", sedchange[2],"mm-bottom", sep=""),paste("SAR= ",abs(round(sr_sed3,3)),"mm/yr, R2= ", round(summary(lm_sed3)$r.squared,4), ", Error +/- ",abs(round(sr_sed1_err,3)),"mm/yr",sep="")))
        }
      }
    } else {
      if (max(sedchange)>0) {
        metadata <- rbind(metadata,
                          c(paste("Mass accumulation rate (CFCS model) ", SML,"-",sedchange[1],"mm",sep=""),paste("MAR= ",abs(round(sr_sed1,3)),"g/cm2/yr, R2= ", round(summary(lm_sed1)$r.squared,4),", Error +/- ",abs(round(sr_sed1_err,3)),"g/cm2/yr", sep="")))
      } else {
        metadata <- rbind(metadata,
                          c("Mass accumulation rate (CFCS model)",paste("MAR= ",abs(round(sr_sed1,3)),"g/cm2/yr, R2= ", round(summary(lm_sed1)$r.squared,4),", Error +/- ",abs(round(sr_sed1_err,3)),"g/cm2/yr", sep="")))
      }
      if (max(sedchange)>0) {
        if(length(sedchange)==1) {
          metadata <- rbind(metadata,
                            c(paste("Mass accumulation rate (CFCS model) ", sedchange[1],"mm-bottom", sep=""),paste("MAR= ",abs(round(sr_sed2,3)),"g/cm2/yr, R2= ", round(summary(lm_sed2)$r.squared,4), ", Error +/- ",abs(round(sr_sed1_err,3)),"g/cm2/yr",sep="")))
        }
        if(length(sedchange)==2) {
          metadata <- rbind(metadata,
                            c(paste("Mass accumulation rate (CFCS model) ", sedchange[1],"-",sedchange[2],"mm", sep=""),paste("MAR= ",abs(round(sr_sed2,3)),"g/cm2/yr, R2= ", round(summary(lm_sed2)$r.squared,4), ", Error +/- ",abs(round(sr_sed1_err,3)),"g/cm2/yr",sep="")),
                            c(paste("Mass accumulation rate (CFCS model) ", sedchange[2],"mm-bottom", sep=""),paste("MAR= ",abs(round(sr_sed3,3)),"g/cm2/yr, R2= ", round(summary(lm_sed3)$r.squared,4), ", Error +/- ",abs(round(sr_sed1_err,3)),"g/cm2/yr",sep="")))
        }
      }
    }
  }

  # Add in the output data file the ages estimations for changes in sed rate
  if(max(sedchange)>0) {
    if(length(sedchange)==2) {
      metadata <- rbind(metadata,
                        c("Best Age (1st change)",paste(abs(round(age_break,0))," (uncertainty: ",abs(round(age_break_low,0)),"-",abs(round(age_break_high,0)),")",sep="")))
      metadata <- rbind(metadata,
                        c("Best Age (2nd change)",paste(abs(round(age_break2,0))," (uncertainty: ",abs(round(age_break2_low,0)),"-",abs(round(age_break2_high,0)),")",sep="")))
    } else {
      metadata <- rbind(metadata,
                        c("Best Age",paste(abs(round(age_break,0))," (uncertainty: ",abs(round(age_break_low,0)),"-",abs(round(age_break_high,0)),")",sep="")))
    }
  }
  # Add in output the inventory of Lead if CRS hypothesis was selected
  if(length(grep("Pb",x = colnames(dt)))>1 & length(grep("density",x = colnames(dt)))>=1) {
    metadata <- rbind(metadata,
                      c("Inventory (Lead)",paste(round(Inventory_CRS[1],3), " Bq/m2 (range: ",round(Inventory_CRS[1]-Inventory_CRS_error[1]), "-",round(Inventory_CRS[1]+Inventory_CRS_error[1])," Bq/m2)\n", sep="")))
  }

  # Add in output the inventory of Cesium if Cs and density available
  if(length(grep("Cs",x = colnames(dt)))>1 & length(grep("density",x = colnames(dt)))>=1) {
    metadata <- rbind(metadata,
                      c("Inventory (Cesium)",paste(round(sum(Inventory_Cesium,na.rm=T),3), " Bq/m2 (range: ",round(sum(Inventory_Cesium_low,na.rm=T)), "-",round(sum(Inventory_Cesium_high,na.rm=T))," Bq/m2)\n", sep="")))
  }

  # Add in the output data file the ages estimations for instantaneous deposit
  if(length(inst_deposit)>1) {
    for (i in 1:nrow(inst_deposit)) {
      metadata <- rbind(metadata,
                        c(paste("Age instantenous deposit ",i," (",inst_deposit[i,1],"-",inst_deposit[i,2],"mm)",sep=""),
                          paste("Estimated range (from CFCS model): ", round(output_agemodel_CFCS[which(output_agemodel_CFCS[,1]==inst_deposit[i,1]),3]),"-",round(output_agemodel_CFCS[which(output_agemodel_CFCS[,1]==inst_deposit[i,1]),4]),sep="")))
    }
  }

  if(save_code) {
    # Save the code that was used
    # First save the history to folder
    savehistory(file = "myhistory.Rhistory")

    # The code may extend on several lines (true for RStudio users at least)
    # This loop find the beginning of the function, serac
    whichline=NULL
    for (i in 1:30) {
      if (length(grep(pattern = "serac", rev(readLines(con = "myhistory.Rhistory"))[i]))>0) whichline <- c(whichline,i)
    }
    if(!is.null(whichline))
    {
      whichline <- min(whichline, na.rm=T)
      mycode=NULL
      for (i in whichline:1) {
        mycode <- paste(mycode, rev(readLines(con = "myhistory.Rhistory"))[i], sep="")
      }
    } else {
      mycode <- "serac could not read the code you used to produce your model"
    }

    out_list$Rcode <- mycode

    # Remove the history
    if (file.exists("myhistory.Rhistory")) file.remove("myhistory.Rhistory")

    metadata <- rbind(metadata,
                      c("",""),
                      c("code",mycode)) # Add line to metadata
  }

  # Write final output
  write.table(x = metadata, file = paste(getwd(),"/Cores/",name,"/",name,"_Metadata_",Sys.Date(),".txt",sep = ""),col.names = F, row.names = F, sep = "\t")


  #### 6. FINAL PLOT -----
  if(preview|plotpdf) {
    # First, to avoid any potential issue or past parameters from the user, create and delete a plot
    {plot(0,0, axes=F, xlab="",ylab="",pch=NA); dev.off(); dev.new()}
    # size for output plot
    cex_1=.8*mycex
    cex_2=1.1*mycex
    cex_4=1*mycex #Writing within plot of sedimentation rate

    # Set plot panels
    mylayout <- NULL # mylayout vector will set the width of the different windows within the plot

    if(!mass_depth) { # Layout if depth in mm is used
      if(plotphoto) mylayout <- c(mylayout,.2)
      if(suppdescriptor) mylayout <- c(mylayout,1)
      if(plot_Pb) mylayout <- c(mylayout,1)
      if(plot_Pb_inst_deposit) mylayout <- c(mylayout,1.3)
      if(plot_Cs) mylayout <- c(mylayout,1.3)
      mylayout <- c(mylayout,1.6)
    } else {          # Layout if mass_depth is used
      if(plot_Pb) mylayout <- c(mylayout,1)
      if(plot_Pb_inst_deposit) mylayout <- c(mylayout,1.3)
      if(plot_Cs) mylayout <- c(mylayout,1.3)
      if(plotphoto) mylayout <- c(mylayout,.6,.6) # more space than for the mm situation because we need place here for the second scale in g.cm2
      if(suppdescriptor&plotphoto) mylayout <- c(mylayout,.9)
      if(suppdescriptor&!plotphoto) mylayout <- c(mylayout,1.8)
      mylayout <- c(mylayout,2.3)
    }

    mylayout[1] <- mylayout[1]+.3 #Add margin to the right windows to include the depth_avg scale
    nwindows <- length(mylayout)

    # If inst_deposit present, we do not always want to plot all of the inst_deposit
    # This length_id vector identifies which instantaneous deposit have corresponding density | radionuclides data
    if(length(inst_deposit)>1) length_id <- length(inst_deposit[inst_deposit[,1]<=max(dt$depth_top[!is.na(dt$density)],na.rm=T),1]) else length_id = 0

    # Save plot to an object using a null PDF device
    pdf(NULL)
    dev.control(displaylist="enable")

    # Layout
    layout(matrix(c(1:nwindows),1,nwindows), widths = mylayout)
    # 6.1. Add core photo ####
    if(plotphoto & !mass_depth) {
      par(mar=c(4.1,3.3,4.1,0))
      plot(c(0,1),myylim, xlab="",ylab="", axes=F, type="n",ylim=myylim)
      axis(2, at = seq(min(myylim),0,by=10), NA, cex.axis=cex_2, lwd=.3)
      axis(2, at = -(pretty(seq(dmin,dmax,5))), labels=pretty(seq(dmin,dmax,5)), cex.axis=cex_2)
      mtext(text = "Depth (mm)", side = 2, line=2.2, cex=cex_1)

      par(xpd=TRUE)
      if(inst_deposit_present) for (i in 1:nrow(inst_deposit)) rect(xleft = .5, ybottom = -inst_deposit[i,2], xright = 3, ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
      if(inst_deposit_present) rect(xleft = -2,ybottom = -dmax*1.2,xright = 3,ytop = -dmax,col = "white",border = "white", density = 1)
      if(SML>0) rect(xleft = .5, ybottom = -SML, xright = 3, ytop = 0, col=grey(0.97), border=NA)
      par(xpd=FALSE)

      rasterImage(photo,xleft = 0,xright = 1,ytop = -minphoto, ybottom = -maxphoto)
    }

    # 6.2 Descriptor ####
    if(suppdescriptor & !mass_depth) {
      if(!exists("suppdescriptorcol")) suppdescriptorcol=c("black","purple")
      dt_suppdescriptor <- dt_suppdescriptor[dt_suppdescriptor$Depth<=max(abs(myylim)),]
      if(plotphoto) {
        par(mar=c(4.1,1.1,4.1,1.1))
        plot(dt_suppdescriptor[,2],-dt_suppdescriptor[,1], xlab="",ylab="", axes=F, type="n",ylim=myylim)
        myxlim_min=min(dt_suppdescriptor[,2],na.rm=T)-1*(max(dt_suppdescriptor[,2],na.rm=T)-min(dt_suppdescriptor[,2],na.rm=T))
        myxlim_max=max(dt_suppdescriptor[,2],na.rm=T)+1*(max(dt_suppdescriptor[,2],na.rm=T)-min(dt_suppdescriptor[,2],na.rm=T))

        par(xpd=TRUE)
        if(inst_deposit_present) for (i in 1:nrow(inst_deposit)) rect(xleft = myxlim_min, ybottom = -inst_deposit[i,2], xright = myxlim_max, ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
        if(SML>0) rect(xleft = myxlim_min, ybottom = -SML, xright = myxlim_max, ytop = 0, col=grey(0.97), border=NA)
        par(xpd=FALSE)

        points(dt_suppdescriptor[,2],-dt_suppdescriptor[,1], pch=16, cex=.8, col=suppdescriptorcol[1])
        lines(dt_suppdescriptor[,2],-dt_suppdescriptor[,1], col=suppdescriptorcol[1])
      } else {
        par(mar=c(4.1,4.1,4.1,1.1))
        plot(dt_suppdescriptor[,2],-dt_suppdescriptor[,1], xlab="",ylab="", axes=F, type="n",ylim=myylim)
        myxlim_min=min(dt_suppdescriptor[,2],na.rm=T)-.5*(max(dt_suppdescriptor[,2],na.rm=T)-min(dt_suppdescriptor[,2],na.rm=T))
        myxlim_max=max(dt_suppdescriptor[,2],na.rm=T)+.5*(max(dt_suppdescriptor[,2],na.rm=T)-min(dt_suppdescriptor[,2],na.rm=T))

        if(inst_deposit_present) for (i in 1:nrow(inst_deposit)) rect(xleft = myxlim_min, ybottom = -inst_deposit[i,2], xright = max(dt_suppdescriptor[,2],na.rm=T), ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
        if(SML>0) {
          rect(xleft = myxlim_min, ybottom = -SML, xright = max(dt_suppdescriptor[,2],na.rm=T), ytop = 0, col=grey(0.97), border=NA)
          abline(h=-SML, lwd=.6, col="darkgrey")
        }
        par(xpd=TRUE)
        if(inst_deposit_present) for (i in 1:nrow(inst_deposit)) rect(xleft = max(dt_suppdescriptor[,2],na.rm=T), ybottom = -inst_deposit[i,2], xright = myxlim_max, ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
        if(SML>0) rect(xleft = myxlim_min, ybottom = -SML, xright = myxlim_max, ytop = 0, col=grey(0.97), border=NA)
        par(xpd=FALSE)

        points(dt_suppdescriptor[,2],-dt_suppdescriptor[,1], pch=16, cex=.8, col=suppdescriptorcol[1])
        lines(dt_suppdescriptor[,2],-dt_suppdescriptor[,1], col=suppdescriptorcol[1])
        #add y axis if first window to be plotted
        axis(2, at = seq(min(myylim),0,by=10), NA, cex.axis=cex_2, lwd=.5)
        axis(2, at = -(pretty(seq(dmin,dmax,5))), labels=pretty(seq(dmin,dmax,5)), cex.axis=cex_2)
        mtext(text = "Depth (mm)", side = 2, line=2.2, cex=cex_1)
      }
      axis(3, cex.axis=cex_2)
      mtext(text = descriptor_lab[1], side = 3, line=2.2, cex=cex_1)

      if(length(descriptor_lab)>1) {
        points(dt_suppdescriptor[,3],-dt_suppdescriptor[,1], pch=1, cex=.8, col=suppdescriptorcol[2])
        lines(dt_suppdescriptor[,3],-dt_suppdescriptor[,1], col=suppdescriptorcol[2])
        points(dt_suppdescriptor[,3],-dt_suppdescriptor[,1], pch=20, cex=.95,col="white")
        axis(1)
        mtext(text = descriptor_lab[2], side = 1, line=2.2, cex=cex_1)
        legend("bottomright", legend = descriptor_lab, bty="n", pch=c(16,1), col=suppdescriptorcol, cex=mycex)
      }
    }

    # 6.3.a Plot 210Pb ####
    # 6.3.a.1 Plot 210Pb in mm ####
    if(plot_Pb) {
      if(!mass_depth) { # default
        if(plotphoto || suppdescriptor) {
          par(mar=c(4.1,1.1,4.1,1.1))
          plot(dt$Pbex,-dt$depth_avg,xlab="",ylab="", axes="F", type="n",xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim)
          myxlim_min=min(log(dt$Pbex),na.rm=T)-.5*(max(log(dt$Pbex),na.rm=T)-min(log(dt$Pbex),na.rm=T))
          myxlim_max=max(log(dt$Pbex),na.rm=T)+.5*(max(log(dt$Pbex),na.rm=T)-min(log(dt$Pbex),na.rm=T))

          par(xpd=TRUE)
          if(inst_deposit_present)  for (i in 1:nrow(inst_deposit)) rect(xleft = log(.1), ybottom = -inst_deposit[i,2], xright = log(15000), ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
          if(SML>0) rect(xleft = log(.1), ybottom = -SML, xright = log(15000), ytop = 0, col=grey(0.97), border=NA)
          par(xpd=FALSE)

        } else {
          par(mar=c(4.1,4.1,4.1,1.1))
          plot(dt$Pbex,-dt$depth_avg,xlab="",ylab="", axes="F", type="n",xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim)
          myxlim_min=min(log(dt$Pbex),na.rm=T)-.5*(max(log(dt$Pbex),na.rm=T)-min(log(dt$Pbex),na.rm=T))
          myxlim_max=max(log(dt$Pbex),na.rm=T)+.5*(max(log(dt$Pbex),na.rm=T)-min(log(dt$Pbex),na.rm=T))

          if(inst_deposit_present)  for (i in 1:nrow(inst_deposit)) rect(xleft = log(.1), ybottom = -inst_deposit[i,2], xright = log(max(log(dt$Pbex),na.rm=T)), ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
          if(SML>0) rect(xleft = log(.1), ybottom = -SML, xright = log(max(log(dt$Pbex),na.rm=T)), ytop = 0, col=grey(0.97), border=NA)
          par(xpd=T)
          if(inst_deposit_present)  for (i in 1:nrow(inst_deposit)) rect(xleft = log(15000), ybottom = -inst_deposit[i,2], xright = log(max(log(dt$Pbex),na.rm=T)), ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
          if(SML>0) rect(xleft = log(15000), ybottom = -SML, xright = log(max(log(dt$Pbex),na.rm=T)), ytop = 0, col=grey(0.97), border=NA)
          par(xpd=F)

          axis(2, at = seq(min(myylim),0,by=10), NA, cex.axis=cex_2, lwd=.5)
          axis(2, at = -(pretty(seq(dmin,dmax,5))), labels=pretty(seq(dmin,dmax,5)), cex.axis=cex_2)
          mtext(text = "Depth (mm)", side = 2, line=2.2, cex=cex_1)
        }

        par(new=T)
        with (
          data=dt_sed1[!is.na(dt_sed1$depth_avg_2),]
          , expr = errbar(log(Pbex),-depth_avg,c(-depth_avg+thickness/2),c(-depth_avg-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim, col=Pbcol[1], errbar.col = Pbcol[1], cex=.8)
        )
        for (i in which(dt_sed1$depth_avg_2>0 & !is.na(dt_sed1$Pbex_er))) {
          lines(c(log(dt_sed1$Pbex[i]+dt_sed1$Pbex_er[i]),log(dt_sed1$Pbex[i]-dt_sed1$Pbex_er[i])),
                rep(-dt_sed1$depth_avg[i],2), type="o", pch="|", cex=.5, col=Pbcol[1])
        }

        if (max(sedchange)>0) {
          par(new=T)
          with (
            data=dt_sed2[!is.na(dt_sed2$depth_avg_2),]
            , expr = errbar(log(Pbex),-depth_avg,c(-depth_avg+thickness/2),c(-depth_avg-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim, col=Pbcol[2], errbar.col = Pbcol[2])
          )
          for (i in which(dt_sed2$depth_avg>0 & !is.na(dt_sed2$Pbex_er))) {
            lines(c(log(dt_sed2$Pbex[i]+dt_sed2$Pbex_er[i]),log(dt_sed2$Pbex[i]-dt_sed2$Pbex_er[i])),
                  rep(-dt_sed2$depth_avg[i],2), type="o", pch="|", cex=.5, col=Pbcol[2])
          }
          if (length(sedchange)==2) {
            par(new=T)
            with (
              data=dt_sed3[!is.na(dt_sed3$depth_avg_2),]
              , expr = errbar(log(Pbex),-depth_avg,c(-depth_avg+thickness/2),c(-depth_avg-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim, col=Pbcol[3], errbar.col = Pbcol[3])
            )
            for (i in which(dt_sed3$depth_avg>0 & !is.na(dt_sed3$Pbex_er))) {
              lines(c(log(dt_sed3$Pbex[i]+dt_sed3$Pbex_er[i]),log(dt_sed3$Pbex[i]-dt_sed3$Pbex_er[i])),
                    rep(-dt_sed3$depth_avg[i],2), type="o", pch="|", cex=.5, col=Pbcol[3])
            }
          }
        }

        # Add 'ignore' values
        par(new=T)
        with (data=dt[is.na(dt$depth_avg_2),]
              , expr = errbar(log(Pbex),-depth_avg,c(-depth_avg+thickness/2),c(-depth_avg-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim, col=grey(.65), errbar.col = grey(.65), cex=.8)
        )
        for (i in which(is.na(dt$depth_avg_2))) {
          lines(c(log(dt$Pbex[i]+dt$Pbex_er[i]),log(dt$Pbex[i]-dt$Pbex_er[i])),
                rep(-dt$depth_avg[i],2), type="o", pch="|", cex=.5, col=grey(.65))
        }
      } else { # if plot against massic depth
        # 6.3.a.2 Plot 210Pb in g/cm2 ####
        par(mar=c(4.1,4.1,4.1,1.1))
        plot(dt$Pbex,-dt$mass_depth_avg,xlab="",ylab="", axes="F", type="n",xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim_md)
        myxlim_min=min(log(dt$Pbex),na.rm=T)-.5*(max(log(dt$Pbex),na.rm=T)-min(log(dt$Pbex),na.rm=T))
        myxlim_max=max(log(dt$Pbex),na.rm=T)+.5*(max(log(dt$Pbex),na.rm=T)-min(log(dt$Pbex),na.rm=T))

        if(inst_deposit_present&length_id>0)  for (i in 1:length_id) rect(xleft = log(.1), ybottom = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,2]))], xright = log(max(log(dt$Pbex),na.rm=T)), ytop = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,1]))],col=inst_depositcol, border=inst_depositcol, lwd=.4)
        if(SML>0) rect(xleft = log(.1), ybottom = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - SML))], xright = log(max(log(dt$Pbex),na.rm=T)), ytop = 0, col=grey(0.97), border=NA)
        par(xpd=T)
        if(inst_deposit_present&length_id>0)  for (i in 1:length_id) rect(xleft = log(15000), ybottom = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,2]))], xright = log(max(log(dt$Pbex),na.rm=T)), ytop = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,1]))],col=inst_depositcol, border=inst_depositcol, lwd=.4)
        if(SML>0) rect(xleft = log(15000), ybottom = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - SML))], xright = log(max(log(dt$Pbex),na.rm=T)), ytop = 0, col=grey(0.97), border=NA)
        par(xpd=F)

        axis(2, at = pretty(seq(myylim_md[1], myylim_md[2], length.out = 20),n=40), labels = NA, cex.axis=cex_2, lwd=.5)
        axis(2, at = pretty(seq(myylim_md[1], myylim_md[2], length.out = 5)), labels=-(pretty(seq(myylim_md[1], myylim_md[2], length.out = 5))), cex.axis=cex_2)
        mtext(text = bquote("Mass depth (g.cm"*~""^-2*")"), side = 2, line=2.2, cex=cex_1)


        par(new=T)
        with (
          data=dt_sed1[!is.na(dt_sed1$depth_avg_2),]
          , expr = errbar(log(Pbex),-mass_depth_avg,-mass_depth_bottom,-mass_depth_top, pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim_md, col=Pbcol[1], errbar.col = Pbcol[1], cex=.8)
        )
        for (i in which(dt_sed1$depth_avg_2>0 & !is.na(dt_sed1$Pbex_er))) {
          lines(c(log(dt_sed1$Pbex[i]+dt_sed1$Pbex_er[i]),log(dt_sed1$Pbex[i]-dt_sed1$Pbex_er[i])),
                rep(-dt_sed1$mass_depth_avg[i],2), type="o", pch="|", cex=.5, col=Pbcol[1])
        }

        if (max(sedchange)>0) {
          par(new=T)
          with (
            data=dt_sed2[!is.na(dt_sed2$depth_avg_2),]
            , expr = errbar(log(Pbex),-mass_depth_avg,-mass_depth_top,-mass_depth_bottom, pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim_md, col=Pbcol[2], errbar.col = Pbcol[2])
          )
          for (i in which(dt_sed2$depth_avg>0 & !is.na(dt_sed2$Pbex_er))) {
            lines(c(log(dt_sed2$Pbex[i]+dt_sed2$Pbex_er[i]),log(dt_sed2$Pbex[i]-dt_sed2$Pbex_er[i])),
                  rep(-dt_sed2$mass_depth_avg[i],2), type="o", pch="|", cex=.5, col=Pbcol[2])
          }
          if (length(sedchange)==2) {
            par(new=T)
            with (
              data=dt_sed3[!is.na(dt_sed3$depth_avg_2),]
              , expr = errbar(log(Pbex),-mass_depth_avg,-mass_depth_top,-mass_depth_bottom, pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim_md, col=Pbcol[3], errbar.col = Pbcol[3])
            )
            for (i in which(dt_sed3$depth_avg>0 & !is.na(dt_sed3$Pbex_er))) {
              lines(c(log(dt_sed3$Pbex[i]+dt_sed3$Pbex_er[i]),log(dt_sed3$Pbex[i]-dt_sed3$Pbex_er[i])),
                    rep(-dt_sed3$mass_depth_avg[i],2), type="o", pch="|", cex=.5, col=Pbcol[3])
            }
          }
        }

        # Add 'ignore' values
        par(new=T)
        with (data=dt[is.na(dt$depth_avg_2),]
              , expr = errbar(log(Pbex),-mass_depth_avg,-mass_depth_bottom,-mass_depth_top, pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim_md, col=grey(.65), errbar.col = grey(.65), cex=.8)
        )
        for (i in which(is.na(dt$depth_avg_2))) {
          lines(c(log(dt$Pbex[i]+dt$Pbex_er[i]),log(dt$Pbex[i]-dt$Pbex_er[i])),
                rep(-dt$mass_depth_avg[i],2), type="o", pch="|", cex=.5, col=grey(.65))
        }
      }

      # create the flexible axis ticks for 210Pb
      power <- 1
      for (i in 0:3) power <- c(power, 2:10*10^c(i))
      label <- power
      label[grep("1",power,invert = TRUE)] <- NA
      axis(3, at=log(power[power<exp(myxlim_max)]), labels=label[power<exp(myxlim_max)])
      # axis label
      mtext(text = bquote(~""^210*Pb[ex]*" (mBq/g)"), side = 3, line=2.2, cex=cex_1)

    }

    # 6.3.b Plot 210Pb without inst_deposits ####
    # 6.3.b.1 Plot 210Pb without inst_deposits in mm ####
    if(plot_Pb_inst_deposit) {
      if (inst_deposit_present) {
        if(!mass_depth) { # default, not plotting against mass depth
          if(plotphoto || suppdescriptor || plot_Pb) {
            par(mar=c(4.1,1.1,4.1,1.1))
            with (
              data=dt_sed1
              , expr = errbar(log(Pbex),-d,c(-d+thickness/2),c(-d-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim, col=Pbcol[1], errbar.col = Pbcol[1], cex=.8)
            )

            par(xpd=T)
            if(SML>0) rect(xleft = log(.1), ybottom = -SML, xright = log(18000), ytop = 0, col=grey(0.97), border=NA)
            if(inst_deposit_present) {
              for (i in 1:nrow(inst_deposit)) rect(xleft = log(.1), ybottom = -inst_deposit[i,2], xright = log(.8), ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
              for (i in 1:nrow(inst_deposit)) {
                pol_x <- c(log(.8),log(2),log(max(dt$Pbex,na.rm=T)),log(max(dt$Pbex,na.rm=T))+log(2)-log(.8),log(max(dt$Pbex,na.rm=T))+log(2)-log(.8),log(max(dt$Pbex,na.rm=T)),log(2),log(.8))
                pol_y <- c(-inst_deposit[i,1],-inst_deposit_corr[i,1],-inst_deposit_corr[i,1],-inst_deposit[i,1],-inst_deposit[i,2],-inst_deposit_corr[i,1],-inst_deposit_corr[i,1],-inst_deposit[i,2])
                polygon(x=pol_x, y = pol_y, col=inst_depositcol, border=NA)
                lines(c(log(2),log(max(dt$Pbex,na.rm=T))), c(-inst_deposit_corr[i,1],-inst_deposit_corr[i,1]),col=inst_depositcol, lwd=.5)
              }
              for (i in 1:nrow(inst_deposit)) rect(xleft = log(max(dt$Pbex,na.rm=T))+log(2)-log(.8), ybottom = -inst_deposit[i,2], xright = log(18000), ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
              points(log(dt_sed1$Pbex),-dt_sed1$d, pch=16, cex=.8)
            }
            par(xpd=F)


          } else {
            par(mar=c(4.1,4.1,4.1,1.1))
            with (
              data=dt_sed1
              , expr = errbar(log(Pbex),-d,c(-d+thickness/2),c(-d-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim, col=Pbcol[1], errbar.col = Pbcol[1], cex=.8)
            )

            if(inst_deposit_present) {
              if(SML>0) rect(xleft = log(.1), ybottom = -SML, xright = log(18000), ytop = 0, col=grey(0.97), border=NA)
              par(xpd=T)
              for (i in 1:nrow(inst_deposit)) {
                pol_x <- c(log(.5),log(2),log(max(dt$Pbex,na.rm=T)),log(max(dt$Pbex,na.rm=T))+log(2)-log(.5),log(max(dt$Pbex,na.rm=T))+log(2)-log(.5),log(max(dt$Pbex,na.rm=T)),log(2),log(.5))
                pol_y <- c(-inst_deposit_corr[i,1],-inst_deposit_corr[i,1],-inst_deposit_corr[i,1],-inst_deposit[i,1],-inst_deposit[i,2],-inst_deposit_corr[i,1],-inst_deposit_corr[i,1],-inst_deposit_corr[i,1])
                polygon(x=pol_x, y = pol_y, col=inst_depositcol, border=NA)
              }
              for (i in 1:nrow(inst_deposit)) rect(xleft = log(max(dt$Pbex,na.rm=T))+log(2)-log(.5), ybottom = -inst_deposit[i,2], xright = log(18000), ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
              par(xpd=F)
              points(log(dt_sed1$Pbex),-dt_sed1$d, pch=16, cex=.8)
            }

          }
        } else {
          # 6.3.b.2 Plot 210Pb without inst_deposits in g/cm2 ####
          if(plot_Pb) {
            par(mar=c(4.1,1.1,4.1,1.1))
            with (
              data=dt_sed1
              , expr = errbar(log(Pbex),-mass_depth_avg_corr,-mass_depth_bottom_corr,-mass_depth_top_corr, pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim_md, col=Pbcol[1], errbar.col = Pbcol[1], cex=.8)
            )

            par(xpd=T)
            if(!is.na(SML) & SML>0) rect(xleft = log(.1), ybottom = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - SML))] , xright = log(18000), ytop = 0, col=grey(0.97), border=NA)
            if(inst_deposit_present&length_id>0) {
              for (i in 1:length_id) rect(xleft = log(.1), ybottom = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,2]))], xright = log(.8), ytop = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,1]))],col=inst_depositcol, border=inst_depositcol, lwd=.4)
              for (i in 1:length_id) {
                pol_x <- c(log(.8),log(2),log(max(dt$Pbex,na.rm=T)),log(max(dt$Pbex,na.rm=T))+log(2)-log(.8),
                           log(max(dt$Pbex,na.rm=T))+log(2)-log(.8),log(max(dt$Pbex,na.rm=T)),log(2),log(.8))
                pol_y <- c(-md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,1]))],
                           -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit_corr[i,1]))],
                           -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit_corr[i,1]))],
                           -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,1]))],
                           -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,2]))],
                           -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit_corr[i,1]))],
                           -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit_corr[i,1]))],
                           -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,2]))])
                polygon(x=pol_x, y = pol_y, col=inst_depositcol, border=NA)
                lines(c(log(2),log(max(dt$Pbex,na.rm=T))), c(-md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit_corr[i,1]))],-md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit_corr[i,1]))]),col=inst_depositcol, lwd=.5)
              }
              for (i in 1:nrow(inst_deposit)) rect(xleft = log(max(dt$Pbex,na.rm=T))+log(2)-log(.8),
                                                   ybottom = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,2]))], xright = log(18000),
                                                   ytop = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,1]))],col=inst_depositcol, border=inst_depositcol, lwd=.4)
              points(log(dt_sed1$Pbex),-dt_sed1$mass_depth_avg_corr, pch=16, cex=.8)
            }
            par(xpd=F)


          } else {
            par(mar=c(4.1,4.1,4.1,1.1))
            with (
              data=dt_sed1
              , expr = errbar(log(Pbex),-mass_depth_avg_corr,-mass_depth_bottom_corr,-mass_depth_top_corr, pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim_md, col=Pbcol[1], errbar.col = Pbcol[1], cex=.8)
            )

            if(inst_deposit_present&length_id>0) {
              if(SML>0) rect(xleft = log(.1), ybottom = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - SML))] , xright = log(18000), ytop = 0, col=grey(0.97), border=NA)
              par(xpd=T)
              for (i in 1:length_id) {
                pol_x <- c(log(.5),log(2),log(max(dt$Pbex,na.rm=T)),log(max(dt$Pbex,na.rm=T))+log(2)-log(.5),log(max(dt$Pbex,na.rm=T))+log(2)-log(.5),log(max(dt$Pbex,na.rm=T)),log(2),log(.5))
                pol_y <- c(-md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,1]))],
                           -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit_corr[i,1]))],
                           -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit_corr[i,1]))],
                           -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,1]))],
                           -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,2]))],
                           -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit_corr[i,1]))],
                           -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit_corr[i,1]))],
                           -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,2]))])
                polygon(x=pol_x, y = pol_y, col=inst_depositcol, border=NA)
              }
              for (i in 1:length_id) rect(xleft = log(max(dt$Pbex,na.rm=T))+log(2)-log(.5), -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,2]))], xright = log(18000), ytop = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,1]))],col=inst_depositcol, border=inst_depositcol, lwd=.4)
              par(xpd=F)
              points(log(dt_sed1$Pbex),-dt_sed1$mass_depth_avg_corr, pch=16, cex=.8)
            }

            axis(2, at = pretty(seq(myylim_md[1], myylim_md[2], length.out = 20),n=40), labels = NA, cex.axis=cex_2, lwd=.5)
            axis(2, at = pretty(seq(myylim_md[1], myylim_md[2], length.out = 5)), labels=-(pretty(seq(myylim_md[1], myylim_md[2], length.out = 5))), cex.axis=cex_2)
            mtext(text = bquote("Mass depth (g.cm"*~""^-2*")"), side = 2, line=2.2, cex=cex_1)

          }
        }

        if(!mass_depth) which_scale = dt_sed1$d else which_scale = dt_sed1$mass_depth_avg_corr
        for (i in which(dt_sed1$d>0)) {
          lines(c(log(dt_sed1$Pbex[i]+dt_sed1$Pbex_er[i]),log(dt_sed1$Pbex[i]-dt_sed1$Pbex_er[i])),
                rep(-which_scale[i],2), type="o", pch="|", cex=.5, col=Pbcol[1])
        }

        # This was the base graph for the points until the 1st sedimentation rate
        # Now doing the same if changes in sed rate were identified.
        if (max(sedchange)>0) {
          par(new=T)
          if(!mass_depth) {
            with (
              data=dt_sed2
              , expr = errbar(log(Pbex),-d,c(-d+thickness/2),c(-d-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim, col=Pbcol[2], errbar.col = Pbcol[2])
            )
          } else {
            with (
              data=dt_sed2
              , expr = errbar(log(Pbex),-mass_depth_avg_corr,-mass_depth_bottom_corr,-mass_depth_top_corr, pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim_md, col=Pbcol[2], errbar.col = Pbcol[2])
            )
          }

          if(!mass_depth) which_scale = dt_sed2$d else which_scale = dt_sed2$mass_depth_avg_corr
          for (i in which(dt_sed2$d>0)) {
            lines(c(log(dt_sed2$Pbex[i]+dt_sed2$Pbex_er[i]),log(dt_sed2$Pbex[i]-dt_sed2$Pbex_er[i])),
                  rep(-which_scale[i],2), type="o", pch="|", cex=.5, col=Pbcol[2])
          }

          # If there's a second change in sedimentation
          if(length(sedchange)==2) {
            par(new=T)
            if(!mass_depth) {
              with (
                data=dt_sed3
                , expr = errbar(log(Pbex),-d,c(-d+thickness/2),c(-d-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim_md, col=Pbcol[3], errbar.col = Pbcol[3])
              )
            } else {
              with (
                data=dt_sed3
                , expr = errbar(log(Pbex),-mass_depth_avg_corr,-mass_depth_bottom_corr,-mass_depth_top_corr, pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim_md, col=Pbcol[3], errbar.col = Pbcol[3])
              )
            }

            if(!mass_depth) which_scale = dt_sed3$d else which_scale = dt_sed3$mass_depth_avg_corr
            for (i in which(dt_sed3$d>0)) {
              lines(c(log(dt_sed3$Pbex[i]+dt_sed3$Pbex_er[i]),log(dt_sed3$Pbex[i]-dt_sed3$Pbex_er[i])),
                    rep(-which_scale[i],2), type="o", pch="|", cex=.5, col=Pbcol[3])
            }
          }
        }
        # create the flexible axis ticks
        power <- 1
        for (i in 0:3) power <- c(power, 2:10*10^c(i))
        label <- power
        label[grep("1",power,invert = TRUE)] <- NA
        axis(3, at=log(power[power<exp(myxlim_max)]), labels=label[power<exp(myxlim_max)])
        # Axis label
        mtext(text = bquote(~""^210*Pb[ex]*" (mBq/g)"), side = 3, line=2.2, cex=cex_1)

        #Add depth label if last plot
        if(!plot_Cs&mass_depth) {
          #par(xpd=T)
          axis(4, at = pretty(seq(myylim_md[1], myylim_md[2], length.out = 20),n=40), labels = NA, cex.axis=cex_2, lwd=.5,line=1.1)
          axis(4, at = pretty(seq(myylim_md[1], myylim_md[2], length.out = 5)), labels=-(pretty(seq(myylim_md[1], myylim_md[2], length.out = 5))), cex.axis=cex_2,line=1.1)
          mtext(text = bquote("Mass depth (g.cm"*~""^-2*")"), side = 4, line=3.3, cex=cex_1)
          #par(xpd=F)
        }
      }
    }

    # 6.3.c Plot 210Pb model on top of 210Pb measurements ####
    if(any(model=="CFCS")) {
      # Warning, you shouldn't plot the linear regression on the point without instantaneous deposit while you mentioned there were some
      if(exists("inst_deposit")&length(inst_deposit)>1&min(inst_deposit)<max(dt$depth_avg)&plot_Pb_inst_deposit==F&plot_CFCS_regression==F) {
        message("\n Warning. It is not possible to visualise the linear regression calculated from\n the CFCS model on the 210Pbex curve without instantaneous deposits, on the\n initial depth. Add the argument plot_Pb_inst_deposit=TRUE to visualise the\n regression line on 210Pbex profile corrected from instananeous events.\n Keep in mind the regression line won't necesseraly match the points.\n\n")
      }
      # If you decide to do it anyway, this message will display
      if(exists("inst_deposit")&length(inst_deposit)>1&min(inst_deposit)<max(dt$depth_avg)&plot_Pb_inst_deposit==F&plot_CFCS_regression==T) {
        message("\n Warning. It seems you requested to visualise the linear regression calculated\n from the CFCS model on the 210Pbex curve without instantaneous deposits,\n while you specified there were some.\n Please bear in mind the linear regression does not correspond to the points. \n Turn plot_Pb_inst_deposit to TRUE or plot_CFCS_regression to FALSE.\n\n")
      }
      if(plot_Pb & plot_CFCS_regression | plot_Pb_inst_deposit & plot_CFCS_regression) {

        # 210Pb linear model needs to stop in case deepest depth have been set to be 'ignored'
        # (bug observed on 2018-10-26 by PS on CA08 - answer in email RB to PS on 2018-11-13)
        # Also adding more condition to still plot the linear model when a element is set to ignore but there are still values to be consider afterwards
        # When plotting the linear model, we'll use either of these two following lines (after toplm definition)
        # One other special case regarding the next lines (before if()): when the surface samples are ignored, the regression line must not extend to the surface
        # We then define the value 'top linear model' to decrease the number of conditons in the code...
        if (length(ignore)>0) toplm <- max(c(SML,min(dt$depth_avg[!dt$depth_avg %in% ignore], na.rm = T)), na.rm = T) else toplm <- SML
        if (mass_depth) toplm <- md_interp$md_avg[which.min(abs(md_interp$depth_mm - toplm))]

        if (is.null(ignore) || !is.null(is.null(ignore)) && which.min(abs(dt$depth_avg-max(ignore, na.rm=T)))<nrow(dt) && any(!is.na(dt$depth_avg_2[which.min(abs(dt$depth_avg-max(ignore, na.rm=T))):nrow(dt)]))) {
          lines(c(-toplm,-sedchange_corr_allscales[1])~ c(lm_sed1$coefficients[1]+toplm*lm_sed1$coefficients[2],lm_sed1$coefficients[1]+sedchange_corr_allscales[1]*lm_sed1$coefficients[2]), col=Pbcol[1], lwd=2)
        }
        if (length(ignore)>0 && sedchange_corr[1] <= max(ignore, na.rm=T)) {
          if(!mass_depth) lines(c(-toplm,-max(dt$d[!dt$depth_avg %in% ignore & dt$d<sedchange_corr[1]],na.rm=T))~ c(lm_sed1$coefficients[1]+toplm*lm_sed1$coefficients[2],lm_sed1$coefficients[1]+max(dt$d[!dt$depth_avg %in% ignore & dt$d<sedchange_corr[1]],na.rm=T)*lm_sed1$coefficients[2]), col=Pbcol[1], lwd=2)
          if(mass_depth)  lines(c(-toplm,-max(dt$mass_depth_avg_corr[!dt$depth_avg %in% ignore & dt$d<sedchange_corr[1]],na.rm=T))~ c(lm_sed1$coefficients[1]+toplm*lm_sed1$coefficients[2],lm_sed1$coefficients[1]+max(dt$mass_depth_avg_corr[!dt$depth_avg %in% ignore & dt$d<sedchange_corr[1]],na.rm=T)*lm_sed1$coefficients[2]), col=Pbcol[1], lwd=2)
        }

        if(!mass_depth) {
          d_legend <- mean(c(min(dt_sed1$d,na.rm = T),max(dt_sed1$d,na.rm = T)))*.8
          shadowtext(x = 0,y = -d_legend-.06*dmax,labels = bquote(r^2 ~ "=" ~ .(round(summary(lm_sed1)$r.squared,4))), pos = 4, col=Pbcol[1], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4)
          shadowtext(x = 0,y = -d_legend,labels = bquote(SAR ~ "=" ~ .(abs(round(sr_sed1,3))) ~ mm.yr^-1), pos = 4, col=Pbcol[1], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4)
        } else {
          d_legend <- mean(c(min(dt_sed1$mass_depth_avg_corr,na.rm = T),max(dt_sed1$mass_depth_avg_corr,na.rm = T)))*.8
          shadowtext(x = 0,y = -d_legend+.06*myylim_md[1],labels = bquote(r^2 ~ "=" ~ .(round(summary(lm_sed1)$r.squared,4))), pos = 4, col=Pbcol[1], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4*.9)
          shadowtext(x = 0,y = -d_legend,labels = bquote(MAR ~ "=" ~ .(abs(round(sr_sed1,3))) ~ g.mm^-1*.yr^-1), pos = 4, col=Pbcol[1], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4*.9)
        }

        if (max(sedchange)>0) {
          if(length(sedchange)==1) {
            if(!mass_depth) {
              if (is.null(ignore) || !is.null(is.null(ignore)) && max(dt_sed2$depth_avg, na.rm = T) > max(ignore, na.rm=T)) lines(c(-min(dt_sed2$d, na.rm = T),-max(dt_sed2$d, na.rm = T))~ c(lm_sed2$coefficients[1]+min(dt_sed2$d, na.rm = T)*lm_sed2$coefficients[2],lm_sed2$coefficients[1]+max(dt_sed2$d, na.rm = T)*lm_sed2$coefficients[2]), lwd=2, col=Pbcol[2])
              if (length(ignore)>0 && max(dt_sed2$depth_avg, na.rm = T) <= max(ignore, na.rm=T)) lines(c(-sedchange_corr[1],-max(dt$d[!dt$depth_avg %in% ignore],na.rm=T))~ c(lm_sed2$coefficients[1]+min(dt_sed2$d, na.rm = T)*lm_sed2$coefficients[2],lm_sed2$coefficients[1]+max(dt$d[!dt$depth_avg %in% ignore],na.rm=T)*lm_sed2$coefficients[2]), lwd=2, col=Pbcol[2])
              d_legend <- mean(c(min(dt_sed2$d,na.rm = T),max(dt_sed2$d,na.rm = T)))*.8
              shadowtext(x = 0,y = -d_legend-.06*dmax,labels = bquote(r^2 ~ "=" ~ .(round(summary(lm_sed2)$r.squared,4))), pos = 4, col=Pbcol[2], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4)
              shadowtext(x = 0,y = -d_legend,labels = bquote(SAR ~ "=" ~ .(abs(round(sr_sed2,3))) ~ mm.yr^-1), pos = 4, col=Pbcol[2], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4)
            } else {
              if (is.null(ignore) || !is.null(is.null(ignore)) && max(dt_sed2$depth_avg, na.rm = T) > max(ignore, na.rm=T)) lines(c(-min(dt_sed2$mass_depth_avg_corr, na.rm = T),-max(dt_sed2$mass_depth_avg_corr, na.rm = T))~ c(lm_sed2$coefficients[1]+min(dt_sed2$mass_depth_avg_corr, na.rm = T)*lm_sed2$coefficients[2],lm_sed2$coefficients[1]+max(dt_sed2$mass_depth_avg_corr, na.rm = T)*lm_sed2$coefficients[2]), lwd=2, col=Pbcol[2])
              if (length(ignore)>0 && max(dt_sed2$mass_depth_avg_corr, na.rm = T) <= max(ignore, na.rm=T)) lines(c(-min(dt_sed2$mass_depth_avg_corr, na.rm = T),-max(dt$mass_depth_avg_corr[!dt$depth_avg %in% ignore],na.rm=T))~ c(lm_sed2$coefficients[1]+min(dt_sed2$mass_depth_avg_corr, na.rm = T)*lm_sed2$coefficients[2],lm_sed2$coefficients[1]+max(dt$mass_depth_avg_corr[!dt$depth_avg %in% ignore],na.rm=T)*lm_sed2$coefficients[2]), lwd=2, col=Pbcol[2])
              d_legend <- mean(c(min(dt_sed2$mass_depth_avg_corr,na.rm = T),max(dt_sed2$mass_depth_avg_corr,na.rm = T)))*.8
              shadowtext(x = 0,y = -d_legend+.06*myylim_md[1],labels = bquote(r^2 ~ "=" ~ .(round(summary(lm_sed2)$r.squared,4))), pos = 4, col=Pbcol[2], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4*.9)
              shadowtext(x = 0,y = -d_legend,labels = bquote(MAR ~ "=" ~ .(abs(round(sr_sed2,3))) ~ g.mm^-1*.yr^-1), pos = 4, col=Pbcol[2], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4*.9)
            }
          }
          if(length(sedchange)==2) {
            if(!mass_depth) {
              lines(c(-min(dt_sed2$d, na.rm = T),-max(dt_sed2$d, na.rm = T))~ c(lm_sed2$coefficients[1]+min(dt_sed2$d, na.rm = T)*lm_sed2$coefficients[2],lm_sed2$coefficients[1]+max(dt_sed2$d, na.rm = T)*lm_sed2$coefficients[2]), lwd=2, col=Pbcol[2])
              d_legend <- mean(c(min(dt_sed2$d,na.rm = T),max(dt_sed2$d,na.rm = T)))*.8
              shadowtext(x = 0,y = -d_legend-.06*dmax,labels = bquote(r^2 ~ "=" ~ .(round(summary(lm_sed2)$r.squared,4))), pos = 4, col=Pbcol[2], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1)
              shadowtext(x = 0,y = -d_legend,labels = bquote(SAR ~ "=" ~ .(abs(round(sr_sed2,3))) ~ mm.yr^-1), pos = 4, col=Pbcol[2], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1)

              if (is.null(ignore) || !is.null(is.null(ignore)) && max(dt_sed3$depth_avg, na.rm = T) > max(ignore, na.rm=T))  lines(c(-min(dt_sed3$d, na.rm = T),-max(dt_sed3$d, na.rm = T))~ c(lm_sed3$coefficients[1]+sedchange_corr[2]*lm_sed3$coefficients[2],lm_sed3$coefficients[1]+max(dt_sed3$d, na.rm = T)*lm_sed3$coefficients[2]), lwd=2, col=Pbcol[3])
              if (length(ignore)>0 && max(dt_sed3$depth_avg, na.rm = T) <= max(ignore, na.rm=T)) lines(c(-min(dt_sed3$d, na.rm = T),-max(dt$d[!dt$depth_avg %in% ignore],na.rm=T))~ c(lm_sed3$coefficients[1]+max(dt_sed3$d, na.rm = T)*lm_sed3$coefficients[2],lm_sed3$coefficients[1]+max(dt$d[!dt$depth_avg %in% ignore],na.rm=T)*lm_sed3$coefficients[2]), lwd=2, col=Pbcol[2])
              d_legend <- mean(c(min(dt_sed3$d,na.rm = T),max(dt_sed3$d,na.rm = T)))*.8
              shadowtext(x = 0,y = -d_legend-.06*dmax,labels = bquote(r^2 ~ "=" ~ .(round(summary(lm_sed3)$r.squared,4))), pos = 4, col=Pbcol[3], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4)
              shadowtext(x = 0,y = -d_legend,labels = bquote(SAR ~ "=" ~ .(abs(round(sr_sed3,3))) ~ mm.yr^-1), pos = 4, col=Pbcol[3], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4)
            } else {
              lines(c(-min(dt_sed2$mass_depth_avg_corr, na.rm = T),-max(dt_sed2$mass_depth_avg_corr, na.rm = T))~ c(lm_sed2$coefficients[1]+min(dt_sed2$mass_depth_avg_corr, na.rm = T)*lm_sed2$coefficients[2],lm_sed2$coefficients[1]+max(dt_sed2$mass_depth_avg_corr, na.rm = T)*lm_sed2$coefficients[2]), lwd=2, col=Pbcol[2])
              d_legend <- mean(c(min(dt_sed2$mass_depth_avg_corr,na.rm = T),max(dt_sed2$mass_depth_avg_corr,na.rm = T)))*.8
              shadowtext(x = 0,y = -d_legend+.06*myylim_md[1],labels = bquote(r^2 ~ "=" ~ .(round(summary(lm_sed2)$r.squared,4))), pos = 4, col=Pbcol[2], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4*.9)
              shadowtext(x = 0,y = -d_legend,labels = bquote(MAR ~ "=" ~ .(abs(round(sr_sed2,3))) ~ g.mm^-1*.yr^-1), pos = 4, col=Pbcol[2], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4*.9)

              if (is.null(ignore) || !is.null(is.null(ignore)) && max(dt_sed3$depth_avg, na.rm = T) > max(ignore, na.rm=T))  lines(c(-min(dt_sed3$mass_depth_avg_corr, na.rm = T),-max(dt_sed3$mass_depth_avg_corr, na.rm = T))~ c(lm_sed3$coefficients[1]+min(dt_sed3$mass_depth_avg_corr, na.rm = T)*lm_sed3$coefficients[2],lm_sed3$coefficients[1]+max(dt_sed3$mass_depth_avg_corr, na.rm = T)*lm_sed3$coefficients[2]), lwd=2, col=Pbcol[3])
              if (length(ignore)>0 && max(dt_sed3$depth_avg, na.rm = T) <= max(ignore, na.rm=T)) lines(c(-min(dt_sed3$mass_depth_avg_corr, na.rm = T),-max(dt$mass_depth_avg_corr[!dt$depth_avg %in% ignore],na.rm=T))~ c(lm_sed3$coefficients[1]+min(dt_sed3$mass_depth_avg_corr, na.rm = T)*lm_sed3$coefficients[2],lm_sed3$coefficients[1]+max(dt$mass_depth_avg_corr[!dt$depth_avg %in% ignore],na.rm=T)*lm_sed3$coefficients[2]), lwd=2, col=Pbcol[2])
              d_legend <- mean(c(min(dt_sed3$mass_depth_avg_corr,na.rm = T),max(dt_sed3$mass_depth_avg_corr,na.rm = T)))*.8
              shadowtext(x = 0,y = -d_legend+.06*myylim_md[1],labels = bquote(r^2 ~ "=" ~ .(round(summary(lm_sed3)$r.squared,4))), pos = 4, col=Pbcol[3], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4*.9)
              shadowtext(x = 0,y = -d_legend,labels = bquote(MAR ~ "=" ~ .(abs(round(sr_sed3,3))) ~ g.mm^-1*.yr^-1), pos = 4, col=Pbcol[3], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4*.9)
            }
          }
        }
      }
    }


    # 6.4. 137Cs ####
    if(plot_Cs) {
      if(plotphoto || suppdescriptor || plot_Pb || plot_Pb_inst_deposit) par(mar=c(4.1,1.1,4.1,1.1)) else par(mar=c(4.1,4.1,4.1,1.1))

      if(!mass_depth) myxlim_max <- max(dt$Cs, na.rm=T)*1.2+max(dt$Cs_er, na.rm = T) else myxlim_max <- max(dt$Cs, na.rm=T)*1.4+max(dt$Cs_er, na.rm = T)
      myxlim_min <- min(dt$Cs, na.rm=T)-max(dt$Cs_er, na.rm = T)

      if(!mass_depth) {
        with (
          data=dt[dt$depth_avg<SML,]
          , expr = errbar(Cs,-depth_avg,c(-depth_avg+thickness/2),c(-depth_avg-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,ylim=myylim, xlim=c(myxlim_min,myxlim_max),col=grey(.65), errbar.col = grey(.65), cex=.8)
        )
        par(xpd=TRUE)
        if(inst_deposit_present) for (i in 1:nrow(inst_deposit)) rect(xleft = -2000, ybottom = -inst_deposit[i,2], xright = max(dt$Cs,na.rm=T)*1.7+2000, ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
        if(SML>0) rect(xleft = -2000, ybottom = -SML, xright = max(dt$Cs,na.rm=T)*1.5, ytop = 0, col=grey(0.97), border=NA)
        par(xpd=FALSE)
      } else {
        with (
          data=dt[dt$depth_avg<SML,]
          , expr = errbar(Cs,-mass_depth_avg,-mass_depth_top,-mass_depth_bottom, pch=16, cap=.01, xlab="",ylab="", axes=F,ylim=myylim_md, xlim=c(myxlim_min,myxlim_max),col=grey(.65), errbar.col = grey(.65), cex=.8)
        )
        par(xpd=TRUE)
        if(inst_deposit_present&length_id>0) for (i in 1:length_id) rect(xleft = -2000, ybottom = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,2]))], xright = max(dt$Cs,na.rm=T)*1.7+2000, ytop = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - inst_deposit[i,1]))],col=inst_depositcol, border=inst_depositcol, lwd=.4)
        if(SML>0) rect(xleft = -2000, ybottom = -md_interp$md_avg[which.min(abs(md_interp$depth_mm - SML))], xright = max(dt$Cs,na.rm=T)*1.5, ytop = 0, col=grey(0.97), border=NA)
        par(xpd=FALSE)
      }

      par(new=T)
      if(!mass_depth) {
        with (
          data=dt[dt$depth_avg<SML,]
          , expr = errbar(Cs,-depth_avg,c(-depth_avg+thickness/2),c(-depth_avg-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,ylim=myylim, xlim=c(myxlim_min,myxlim_max),col=grey(.65), errbar.col = grey(.65), cex=.8)
        )
      } else {
        with (
          data=dt[dt$depth_avg<SML,]
          , expr = errbar(Cs,-mass_depth_avg,-mass_depth_top,-mass_depth_bottom, pch=16, cap=.01, xlab="",ylab="", axes=F,ylim=myylim_md, xlim=c(myxlim_min,myxlim_max),col=grey(.65), errbar.col = grey(.65), cex=.8)
        )
      }

      if(!mass_depth) which_scale = dt$depth_avg else which_scale = dt$mass_depth_avg
      for (i in which(dt$Cs>=0 & !is.na(dt$Cs_er) & dt$depth_avg<SML)) {
        lines(c(dt$Cs[i]+dt$Cs_er[i],dt$Cs[i]-dt$Cs_er[i]),
              rep(-which_scale[i],2), type="o", pch="|", cex=.5, col=grey(.65))
      }
      lines(dt$Cs[which(dt$depth_avg<SML+2)],-which_scale[which(dt$depth_avg<SML+2)], col=grey(.65), lwd=.5)

      par(new=T)
      if(!mass_depth) {
        with (
          data=dt[dt$depth_avg>=SML,]
          , expr = errbar(Cs,-depth_avg,c(-depth_avg+thickness/2),c(-depth_avg-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,ylim=myylim, xlim=c(myxlim_min,myxlim_max), col=Pbcol[1], errbar.col = Pbcol[1], cex=.8)
        )
      } else {
        with (
          data=dt[dt$depth_avg>=SML,]
          , expr = errbar(Cs,-mass_depth_avg,-mass_depth_top,-mass_depth_bottom, pch=16, cap=.01, xlab="",ylab="", axes=F,ylim=myylim_md, xlim=c(myxlim_min,myxlim_max), col=Pbcol[1], errbar.col = Pbcol[1], cex=.8)
        )
      }
      for (i in which(dt$Cs>0 & !is.na(dt$Cs_er) & dt$depth_avg>=SML)) {
        lines(c(dt$Cs[i]+dt$Cs_er[i],dt$Cs[i]-dt$Cs_er[i]),
              rep(-which_scale[i],2), type="o", pch="|", cex=.5, col=Pbcol[1])
      }
      lines(dt$Cs[which(dt$depth_avg>=SML&!is.na(dt$Cs))],-which_scale[which(dt$depth_avg>=SML&!is.na(dt$Cs))], lwd=.5)
      lines(dt$Cs[which(dt$depth_avg>=SML)],-which_scale[which(dt$depth_avg>=SML)])

      axis(3,  cex.axis=cex_2)
      mtext(text = bquote(~""^137*"Cs (mBq/g)"), side = 3, line=2.2, cex=cex_1)

      #Add depth label if last plot
      if(mass_depth) {
        #par(xpd=T)
        axis(4, at = pretty(seq(myylim_md[1], myylim_md[2], length.out = 20),n=40), labels = NA, cex.axis=cex_2, lwd=.5,line=1.1)
        axis(4, at = pretty(seq(myylim_md[1], myylim_md[2], length.out = 5)), labels=-(pretty(seq(myylim_md[1], myylim_md[2], length.out = 5))), cex.axis=cex_2, line=1.1)
        mtext(text = bquote("Mass depth (g.cm"*~""^-2*")"), side = 4, line=3.3, cex=cex_1)
        #par(xpd=F)
      }

      # Add text
      par(xpd=TRUE)
      # determine which depth is used according to mass_depth==T/F
      if(mass_depth)  which_scale=dt$mass_depth_avg else which_scale=dt$depth_avg
      #Chernobyl
      if (exists("Cher")&&!is.na(Cher)) {
        lines(rep(max(dt$Cs[which_scale>min(Cher_allscales-.01*(max(dt$which_scale,na.rm=T))) & which_scale<max(Cher_allscales+.01*(max(dt$which_scale,na.rm=T)))],na.rm = T)*1.1,2),c(-Cher_allscales[1],-Cher_allscales[2]), lwd=1.5)
        if(!mass_depth) {
          shadowtext(max(dt$Cs[which_scale>min(Cher_allscales-.01*(max(dt$which_scale,na.rm=T))) & which_scale<max(Cher_allscales+.01*(max(dt$which_scale,na.rm=T)))],na.rm = T)+0.1*max(dt$Cs,na.rm=T),-(min(Cher_allscales)),
                     labels = c("C 1986"), pos = 3,col="black", bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=mycex)
          lines(c(max(dt$Cs[which_scale>min(Cher_allscales-.01*(max(dt$which_scale,na.rm=T))) & which_scale<max(Cher_allscales+.01*(max(dt$which_scale,na.rm=T)))],na.rm = T)*1.1, max(dt$Cs,na.rm = T)*2),rep(peakCher_allscales,2), lty=2)
        } else {
          shadowtext(max(dt$Cs[which_scale>min(Cher_allscales-.01*(max(dt$which_scale,na.rm=T))) & which_scale<max(Cher_allscales+.01*(max(dt$which_scale,na.rm=T)))],na.rm = T)+0.1*max(dt$Cs,na.rm=T),peakCher_allscales,
                     labels = c("C 1986"), pos = 4,col="black", bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=mycex)
        }
      }
      #NWT
      if (exists("NWT")&&!is.na(NWT)) {
        lines(rep(max(dt$Cs[which_scale>min(NWT_allscales-.01*(max(dt$which_scale,na.rm=T))) & which_scale<max(NWT_allscales+.01*(max(dt$which_scale,na.rm=T)))],na.rm = T)*1.2,2),c(-NWT_allscales[1],-NWT_allscales[2]), lwd=1.5)
        if(!mass_depth) {
          if (Hemisphere == "NH") shadowtext(max(dt$Cs,na.rm = T)+0.1*max(dt$Cs,na.rm=T),-(min(NWT_allscales)),
                                             labels = "NWT 1963", pos = 3, col="black",bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=mycex)
          if (Hemisphere == "SH") shadowtext(max(dt$Cs,na.rm = T)+0.1*max(dt$Cs,na.rm=T),-(min(NWT_allscales)),
                                             labels = "NWT 1964/1965", pos = 3, col="black",bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=mycex)
          lines(c(max(dt$Cs[which_scale>min(NWT_allscales-.01*(max(dt$which_scale,na.rm=T))) & which_scale<max(NWT_allscales+.01*(max(dt$which_scale,na.rm=T)))],na.rm = T)*1.2, max(dt$Cs,na.rm = T)*2),rep(peakNWT_allscales,2), lty=2)
        } else {
          if (Hemisphere == "NH") shadowtext(max(dt$Cs,na.rm = T)+0.1*max(dt$Cs,na.rm=T),peakNWT_allscales,
                                             labels = "NWT 1963", pos = 3, col="black",bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=mycex)
          if (Hemisphere == "SH") shadowtext(max(dt$Cs,na.rm = T)+0.1*max(dt$Cs,na.rm=T),peakNWT_allscales,
                                             labels = "NWT 1964/1965", pos = 3, col="black",bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=mycex)
        }
      }
      #First radionuclides fallout
      if (exists("FF")&&!is.na(FF)) {
        lines(rep(max(dt$Cs[which_scale>min(FF_allscales-.01*(max(dt$which_scale,na.rm=T))) & which_scale<max(FF_allscales+.01*(max(dt$which_scale,na.rm=T)))],na.rm = T)*1.2,2),c(-FF_allscales[1],-FF_allscales[2]), lwd=1.5)
        if(!mass_depth) {
          shadowtext(max(dt$Cs,na.rm = T)+0.1*max(dt$Cs,na.rm=T),-(max(FF_allscales)),
                     labels = c("FF 1955"), pos = 1, col="black",bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=mycex)
          lines(c(max(dt$Cs[which_scale>min(FF_allscales-.01*(max(dt$which_scale,na.rm=T))) & which_scale<max(FF_allscales+.01*(max(dt$which_scale,na.rm=T)))],na.rm = T)*1.2, max(dt$Cs,na.rm = T)*2),rep(peakFF_allscales,2), lty=2)
        } else {
          shadowtext(max(dt$Cs,na.rm = T)+0.1*max(dt$Cs,na.rm=T),peakFF_allscales,
                     labels = c("FF 1955"), pos = 1, col="black",bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=mycex)
        }
      }
      par(xpd=FALSE)

      # 6.5. 241Am ####
      if (plot_Am) {
        legend("bottomright", legend = c("Cesium", "Americium"), bty="n", pch=c(16,1), cex=mycex)
        par(new=T, mar=c(4.1,1.1,4.1,6.1))
        myxlim_max <- max(dt$Am, na.rm=T)*1.2+max(dt$Am_er, na.rm=T)
        myxlim_min <- min(dt$Am, na.rm=T)-max(dt$Am_er, na.rm=T)
        if(!mass_depth) {
          with (
            data=dt[which(!is.na(dt$Am)&dt$Am>0&dt$depth_avg>SML),]
            , expr = errbar(Am,-depth_avg,c(-depth_avg+thickness/2),c(-depth_avg-thickness/2), pch=1, cap=.01, xlab="",ylab="", axes=F,ylim=myylim, xlim=c(myxlim_min,myxlim_max),col=Pbcol[1], errbar.col = Pbcol[1],cex=.8)
          )
        } else {
          with (
            data=dt[which(!is.na(dt$Am)&dt$Am>0&dt$depth_avg>SML),]
            , expr = errbar(Am,-mass_depth_avg,-mass_depth_top,-mass_depth_bottom, pch=1, cap=.01, xlab="",ylab="", axes=F,ylim=myylim_md, xlim=c(myxlim_min,myxlim_max),col=Pbcol[1], errbar.col = Pbcol[1],cex=.8)
          )
        }
        axis(1, cex.axis=cex_2)
        if(mass_depth)  which_scale=dt$mass_depth_avg else which_scale=dt$depth_avg
        mtext(text = bquote(~""^241*"Am (mBq/g)"), side = 1, line=2.4, cex=cex_1)
        for (i in which(dt$Am>0 & !is.na(dt$Am_er) & dt$depth_avg>SML)) {
          lines(c(dt$Am[i]+dt$Am_er[i],dt$Am[i]-dt$Am_er[i]),
                rep(-which_scale[i],2), type="o", pch="|", cex=.5, col=Pbcol[1])
        }
        points(dt$Am[which(dt$Am>0&dt$depth_avg>SML)],-which_scale[which(dt$Am>0&dt$depth_avg>SML)], pch=20, col="white")
        points(dt$Am[which(dt$Am>0&dt$depth_avg>SML)],-which_scale[which(dt$Am>0&dt$depth_avg>SML)])
      }
    }


    # 6.6. if(mass_depth) Add core photo ####
    if(plotphoto & mass_depth) {
      par(mar=c(4.1,1,4.1,1))
      plot(c(0,1),myylim, xlab="",ylab="", axes=F, type="n",ylim=myylim)
      par(mar=c(4.1,2.1,4.1,0))
      plot(c(0,1),myylim, xlab="",ylab="", axes=F, type="n",ylim=myylim)
      axis(2, at = seq(min(myylim),0,by=10), NA, cex.axis=cex_2, lwd=.3)
      axis(2, at = -(pretty(seq(dmin,dmax,5))), labels=pretty(seq(dmin,dmax,5)), cex.axis=cex_2)
      mtext(text = "Depth (mm)", side = 2, line=2.2, cex=cex_1)

      par(xpd=TRUE)
      if(inst_deposit_present) for (i in 1:nrow(inst_deposit)) rect(xleft = .5, ybottom = -inst_deposit[i,2], xright = 3, ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
      if(inst_deposit_present) rect(xleft = -2,ybottom = -dmax*1.2,xright = 3,ytop = -dmax,col = "white",border = "white", density = 1)
      if(SML>0) rect(xleft = .5, ybottom = -SML, xright = 3, ytop = 0, col=grey(0.97), border=NA)
      par(xpd=FALSE)

      rasterImage(photo,xleft = 0,xright = 1,ytop = -minphoto, ybottom = -maxphoto)
    }

    # 6.7 if(mass_depth) Descriptor ####
    if(suppdescriptor & mass_depth) {
      if(!exists("suppdescriptorcol")) suppdescriptorcol=c("black","purple")
      dt_suppdescriptor <- dt_suppdescriptor[dt_suppdescriptor$Depth<=max(abs(myylim)),]
      if(plotphoto) {
        par(mar=c(4.1,0.3,4.1,0.3))
        plot(dt_suppdescriptor[,2],-dt_suppdescriptor[,1], xlab="",ylab="", axes=F, type="n",ylim=myylim)
        myxlim_min=min(dt_suppdescriptor[,2],na.rm=T)-2*(max(dt_suppdescriptor[,2],na.rm=T)-min(dt_suppdescriptor[,2],na.rm=T))
        myxlim_max=max(dt_suppdescriptor[,2],na.rm=T)+2*(max(dt_suppdescriptor[,2],na.rm=T)-min(dt_suppdescriptor[,2],na.rm=T))

        par(xpd=TRUE)
        if(inst_deposit_present) for (i in 1:nrow(inst_deposit)) rect(xleft = myxlim_min, ybottom = -inst_deposit[i,2], xright = myxlim_max, ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
        if(SML>0) rect(xleft = myxlim_min, ybottom = -SML, xright = myxlim_max, ytop = 0, col=grey(0.97), border=NA)
        par(xpd=FALSE)

        points(dt_suppdescriptor[,2],-dt_suppdescriptor[,1], pch=16, cex=.8, col=suppdescriptorcol[1])
        lines(dt_suppdescriptor[,2],-dt_suppdescriptor[,1], col=suppdescriptorcol[1])
      } else {
        par(mar=c(4.1,8.1,4.1,0.3))
        plot(dt_suppdescriptor[,2],-dt_suppdescriptor[,1], xlab="",ylab="", axes=F, type="n",ylim=myylim)
        myxlim_min=min(dt_suppdescriptor[,2],na.rm=T)-.5*(max(dt_suppdescriptor[,2],na.rm=T)-min(dt_suppdescriptor[,2],na.rm=T))
        myxlim_max=max(dt_suppdescriptor[,2],na.rm=T)+.5*(max(dt_suppdescriptor[,2],na.rm=T)-min(dt_suppdescriptor[,2],na.rm=T))

        if(inst_deposit_present) for (i in 1:nrow(inst_deposit)) rect(xleft = myxlim_min, ybottom = -inst_deposit[i,2], xright = max(dt_suppdescriptor[,2],na.rm=T), ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
        if(SML>0) {
          rect(xleft = myxlim_min, ybottom = -SML, xright = max(dt_suppdescriptor[,2],na.rm=T), ytop = 0, col=grey(0.97), border=NA)
          abline(h=-SML, lwd=.6, col="darkgrey")
        }
        par(xpd=TRUE)
        if(inst_deposit_present) for (i in 1:nrow(inst_deposit)) rect(xleft = max(dt_suppdescriptor[,2],na.rm=T), ybottom = -inst_deposit[i,2], xright = myxlim_max, ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
        if(SML>0) rect(xleft = myxlim_min, ybottom = -SML, xright = myxlim_max, ytop = 0, col=grey(0.97), border=NA)
        par(xpd=FALSE)

        points(dt_suppdescriptor[,2],-dt_suppdescriptor[,1], pch=16, cex=.8, col=suppdescriptorcol[1])
        lines(dt_suppdescriptor[,2],-dt_suppdescriptor[,1], col=suppdescriptorcol[1])
        #add y axis if first window to be plotted
        axis(2, at = seq(min(myylim),0,by=10), NA, cex.axis=cex_2, lwd=.5)
        axis(2, at = -(pretty(seq(dmin,dmax,5))), labels=pretty(seq(dmin,dmax,5)), cex.axis=cex_2)
        mtext(text = "Depth (mm)", side = 2, line=2.2, cex=cex_1)
      }
      axis(3, cex.axis=cex_2)
      mtext(text = descriptor_lab[1], side = 3, line=2.2, cex=cex_1)

      if(length(descriptor_lab)>1) {
        points(dt_suppdescriptor[,3],-dt_suppdescriptor[,1], pch=1, cex=.8, col=suppdescriptorcol[2])
        lines(dt_suppdescriptor[,3],-dt_suppdescriptor[,1], col=suppdescriptorcol[2])
        points(dt_suppdescriptor[,3],-dt_suppdescriptor[,1], pch=20, cex=.95,col="white")
        axis(1)
        mtext(text = descriptor_lab[2], side = 1, line=2.2, cex=cex_1)
        legend("bottomright", legend = descriptor_lab, bty="n", pch=c(16,1), col=suppdescriptorcol, cex=mycex)
      }
    }


    # 6.8.a plot Age Model ####
    if(!mass_depth) par(mar=c(4.1,1.1,4.1,4.1))
    if(mass_depth&&plotphoto|suppdescriptor) par(mar=c(4.1,0.1,4.1,4.1)) else par(mar=c(4.1,4.1,4.1,4.1))
    plot(c(-min_yr,-mround(coring_yr,10)),c(-dmin, -dmax), xlab="",ylab="", axes=F, type="n",ylim=myylim)

    # Plot the 'historic_test' argument i.e. know dates we want to add
    if (!is.na(historic_test)){
      abline(v = -historic_test, col = adjustcolor("grey", alpha.f = .3), lwd=2)
      shadowtext(-historic_test,rep(-dmin, length(historic_test)),
                 labels = as.character(historic_test), col="black",bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1,cex = .9*mycex)
    }


    if(inst_deposit_present)  {
      for (i in 1:nrow(inst_deposit)) rect(xleft = -coring_yr, ybottom = -inst_deposit[i,2], xright = -min_yr+20, ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
    }

    par(xpd=T)
    if(inst_deposit_present)  {
      if(!mass_depth) for (i in 1:nrow(inst_deposit)) rect(xleft = -2100, ybottom = -inst_deposit[i,2], xright = -coring_yr, ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
      if(mass_depth&&plotphoto|suppdescriptor) for (i in 1:nrow(inst_deposit)) rect(xleft = -2100, ybottom = -inst_deposit[i,2], xright = -coring_yr, ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
    }
    par(xpd=F)


    axis(4, at = seq(min(myylim),0,by=10), NA, cex.axis=cex_2, lwd=.3)
    axis(4, at = -(pretty(seq(dmin,dmax,5))), labels=pretty(seq(dmin,dmax,5)), cex.axis=cex_2)
    mtext(text = "Depth (mm)", side = 4, line=2.2, cex=cex_1)
    axis(3, at = seq(-mround(coring_yr,10),-min_yr+20,20), labels = seq(mround(coring_yr,10),min_yr-20,-20),cex.axis=cex_2)
    mtext(text = "Year (C.E.)", side = 3, line=2.2, cex=cex_1)

    if(any(model=="CFCS")) {
      which_CFCS <- output_agemodel_CFCS$depth<=max(dt$depth_avg[!is.na(dt$d)])
      if(!mass_depth) {
        lines(-output_agemodel_CFCS$BestAD,-output_agemodel_CFCS$depth, col=modelcol[1],lty=2,lwd=.5)
        pol_x <- c(-output_agemodel_CFCS$MinAD, rev(-output_agemodel_CFCS$MaxAD))
        pol_y <- c(-output_agemodel_CFCS$depth, rev(-output_agemodel_CFCS$depth))
      } else {
        pol_x <- c(-output_agemodel_CFCS$MinAD[which_CFCS], rev(-output_agemodel_CFCS$MaxAD[which_CFCS]))
        pol_y <- c(-output_agemodel_CFCS$depth[which_CFCS], rev(-output_agemodel_CFCS$depth[which_CFCS]))
      }
      polygon(x=pol_x, y = pol_y, col=adjustcolor(modelcol[1], alpha.f=0.2), border=NA)
      lines(-output_agemodel_CFCS$BestAD[output_agemodel_CFCS$depth>=SML&which_CFCS],-output_agemodel_CFCS$depth[output_agemodel_CFCS$depth>=SML&which_CFCS], col=modelcol[1])
    }

    if(any(model=="CIC")) {
      # Creating the logical vector which_CIC to plot only the CIC model point with data.
      which_CIC <- output_agemodel_CIC$depth_avg[-1] %in% dt$depth_avg[!is.na(dt$depth_avg_2)] & !is.na(m_CIC_low)& !is.na(m_CIC_high)
      pol_x <- c(-c(coring_yr,m_CIC_low[which_CIC]), c(rev(-m_CIC_high[which_CIC]), -coring_yr))
      pol_y <- c(-output_agemodel_CIC$depth_avg[which_CIC], rev(-output_agemodel_CIC$depth_avg[which_CIC]))
      polygon(x=pol_x, y = pol_y, col=adjustcolor(modelcol[2], alpha.f=0.2), border=NA)
      which_CIC <- output_agemodel_CIC$depth_avg[-1] %in% dt$depth_avg[!is.na(dt$depth_avg_2)] &!is.na(m_CIC)
      lines(-c(coring_yr,m_CIC[which_CIC]),-output_agemodel_CIC$depth_avg[which_CIC], col=modelcol[2])
    }

    if(any(model=="CRS")) {
      if(exists("inst_deposit")&&length(inst_deposit)>1) {
        lines(-new_x_CRS,-new_y_CRS, col=modelcol[3],lty=2,lwd=.5)
        pol_x <- c(-new_x_CRS_low, rev(-new_x_CRS_high))
        pol_y <- c(-new_y_CRS, rev(-new_y_CRS))
        polygon(x=pol_x, y = pol_y, col=adjustcolor(modelcol[3], alpha.f=0.2), border=NA)
        lines(-new_x_CRS[new_y_CRS>=SML&new_y_CRS<=max(dt$depth_avg)],-new_y_CRS[new_y_CRS>=SML&new_y_CRS<=max(dt$depth_avg)], col=modelcol[3])
      } else {
        depth_CRS_plot = -complete_core_depth[whichkeep]
        lines(-m_CRS,depth_CRS_plot, col=modelcol[3],lty=2,lwd=.5)
        pol_x <- c(-m_CRS_low, rev(-m_CRS_high))
        pol_y <- c(depth_CRS_plot, rev(depth_CRS_plot))
        polygon(x=pol_x, y = pol_y, col=adjustcolor(modelcol[3], alpha.f=0.2), border=NA)
        lines(-m_CRS[-depth_CRS_plot>=SML&-depth_CRS_plot<=max(dt$depth_avg)],depth_CRS_plot[-depth_CRS_plot>=SML&-depth_CRS_plot<=max(dt$depth_avg)], col=modelcol[3])
      }

    }

    if(varves) {
      points(-varve$Age,-varve$depth_avg, pch=4)
    }

    if(exists("Cher") | exists("NWT") | exists("FF")) {
      err_dated_depth_avg <- matrix(err_dated_depth_avg[!is.na(err_dated_depth_avg)],nrow=2,byrow = F)
      err_dated_depth_avg <- -abs(err_dated_depth_avg)

      ## plot the age and depth error
      if(!is.null(dates)) {
        for (i in 1:length(dates)) {
          lines(rep(-dates[i],2), c(err_dated_depth_avg[1,i], err_dated_depth_avg[2,i]), type="o", pch="_", col="black")
        }
        for (i in 1:length(dates)) {
          if(!is.na(err_dates_avg[i])) lines(c(-dates[i]+err_dates_avg[i],-dates[i]-err_dates_avg[i]), rep(dates_depth_avg[i],2), col="black")
          if(!is.na(err_dates_avg[i])) points(c(-dates[i]+err_dates_avg[i],-dates[i]-err_dates_avg[i]), pch= "|", rep(dates_depth_avg[i],2), col="black")
        }
        points(-dates,dates_depth_avg, pch=16, cex=.8)

        if(!mass_depth) {
          par(xpd=T)
          for (i in 1:length(dates)) {
            lines(c(-2100,-dates[i]), rep(dates_depth_avg[i],2), lty=2)
          }
          par(xpd=F)
        } else {
          text(-dates,dates_depth_avg,labels = seq_along(mtext_Cs), pch=16, cex=.8, offset = .3,pos = 4)
          mylegend <- c(paste(seq_along(mtext_Cs),mtext_Cs,sep=": ", collapse = "; "),
                        mylegend)
          mypchlegend <- c(16, mypchlegend)
          myltylegend <- c(NA, myltylegend)
          mycollegend <- c( 1, mycollegend)
        }
      }
    }

    # Plot legend
    if(length(model)>1|varves|(mass_depth&any(!is.na(Cher),!is.na(NWT),!is.na(FF)))) {
      legend("bottomleft", legend = mylegend, pch = mypchlegend,lty = myltylegend, col = mycollegend, bty='n', cex=mycex)
    }

    # 6.8.b Plot historic event on age model ####
    if(length(historic_d)>=1 && !all(is.na(historic_d))) {
      historic_d_dt <- matrix(abs(historic_d), ncol = 2, byrow = T)
      for (i in seq_along(historic_a)) {
        if (!is.na(historic_a[i])) {
          lines(x = c(rep(-historic_a[i],2)), y = c(20,-(historic_d_dt[i,1])),lty=3)
          par(xpd=T)
          if (!is.na(historic_n[i])) {
            shadowtext(-(min_yr+(coring_yr-min_yr)*.15),-mean(historic_d_dt[i,], na.rm=T),
                       labels = as.character(historic_n[i]), col="black",bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1,cex = 1*mycex)
          }
          par(xpd=F)
        }
      }
    }

    # Save plot as a pseudo-object
    out_list$plot <- recordPlot()
    invisible(dev.off())

    # Save the plot to pdf
    if(plotpdf) {
      pdf(paste(getwd(),"/Cores/",name,"/",name,".pdf",sep=""),width = 1.2+1.8*nwindows, height = 5,family = "Helvetica")
      replayPlot(out_list$plot)
      dev.off()
    }

    # View the plot if preview==TRUE
    if(preview) {
      grid::grid.newpage()
      replayPlot(out_list$plot)
    }


  }

  # print elapsed time
  new_time <- Sys.time() - old_time # calculate difference
  # print in nice format
  cat("\n\n ________________________\n")
  cat(paste0("\n The calculation took ", round(new_time,4), " ",units(new_time),".\n\n"))

  # Return outlist
  return(invisible(out_list))

}

