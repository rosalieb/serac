#' serac age-depth modelling function
#'
#' This is the main age-depth modelling function. The default values can be changed permanently within this file or temporarily when calling serac(). If there is any options you would like to see included in future version, please contact one of the authors
#'
#' @export
#' @param name Name of the core, given using quotes. Defaults to the core provided with serac. Use preferably the published name of the core for traceability.
#' @param coring_yr Year of coring.
#' @param model Select 1 to 3 item between c("CFCS", "CIC", "CRS"). If several models are selected, they will all be plotted together in the last window.
#' @param Cher If 137Cs measurement were done, where do you detect the Chernobyl peak? The argument is a vector of two depth given in millimeters giving the top and bottom threshold for the 1986 Chernobyl event. The user can run the model without giving any specification before making a decision. In such case, leave the argument empty.
#' @param NWT If 137Cs measurement were done, where do you detect the Nuclear Weapon Test peak? The argument is a vector of two depth given in millimeters giving the top and bottom threshold for the 1960s Nuclear Weapon Test event. The user can run the model without giving any specification before making a decision. In such case, leave the argument empty.
#' @param Hemisphere Chose between North Hemisphere "NH" and South Hemisphere "SH" depending on the location of your system. This argument is required if you chose to plot NWT.
#' @param FF If 137Cs measurement were done, where do you detect the First Fallout period? The argument is a vector of two depth given in millimeters giving the top and bottom threshold for the First Fallout period in 1955. The user can run the model without giving any specification before making a decision. In such case, leave the argument empty.
#' @param inst_deposit Upper and lower depths (in mm) of sections of abrupt accumulation that inst_deposit c() should be excised, e.g., c(100, 120, 185, 195) for two sections of 10.0-12.0 cm and 18.5-19.5 cm depth
#' @param ignore The depth (in mm) of any sample that should be ignored from the age-depth model computation, e.g., c(55) will remove the measurement done at 5.5 cm. The data will be ploted by default in grey on the output graph (you can change this with the inst_depositcol argument)
#' @param plotpdf Logical argument to indicate whether you want the output graph to be saved to your folder.
#' @param preview Logical argument to indicate whether you want the output graph to be ploted. Default is TRUE, and the graph is ploted within your R session. It might be convenient to turn this argument to FALSE if errors keep coming telling you your R window is too small.
#' @param plotphoto Logical argument to indicate whether you want to plot the photo of the core along your age-model. If plotphoto=TRUE, you need to indicate the upper and lower limit of the photo in mm in following arguments.
#' @param minphoto Mandatory if plotphoto=TRUE. Lower limit of the core photo in mm, e.g., minphoto=0 indicates that the photo starts at 0 mm. The photo will automatically be truncated acording to the minimum and maximum depth of the age model given in other arguments.
#' @param maxphoto Mandatory if plotphoto=TRUE. Upper limit of the core photo in mm, e.g., maxphoto=320 indicates that the photo ends at 32 cm. The photo will automatically be truncated acording to the minimum and maximum depth of the age model given in other arguments.
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
#' @param SML Surface Mixed Layer: a depth in mm above which the sediment is considered to be mixed. E.g., SML=30 indicates that the first 3 cm are mixed sediment: the data point are ploted but not included in the Pb models.
#' @param stepout Depth resolution for the file out in mm.
#' @param mycex Graphical parameter: a multiplication factor to increase (mycex>1) ou decrease (mycex<1) label sizes.
#' @param archive_metadata Logical argument. If TRUE, require fields regarding the measurements on the core. Allows missing information; just press 'ENTER' in your computer (leave an empty field).
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

serac <- function(name="", model=c("CFCS"),Cher=c(),NWT=c(),Hemisphere=c(),FF=c(),inst_deposit=c(0),
                  ignore=c(),plotpdf=FALSE,preview=TRUE,plotphoto=FALSE,minphoto=c(),maxphoto=c(),
                  Pbcol=c("black","midnightblue","darkgreen"),inst_depositcol=grey(0.85),
                  modelcol=c("black","red","darkorange"),
                  historic_d=c(),historic_a=c(),historic_n=c(),historic_test=c(),
                  suppdescriptor=FALSE,descriptor_lab=c(),suppdescriptorcol=c("black","purple"),
                  coring_yr=c(),plot_Am=FALSE,plot_Cs=FALSE,plot_Pb=TRUE,
                  plot_Pb_inst_deposit=FALSE,plot_CFCS_regression=c(),
                  varves=FALSE, dmin=c(),dmax=c(),sedchange=c(0),
                  min_yr=1880, SML=c(0), stepout=5, mycex=1,
                  archive_metadata=FALSE)
  .serac(name, model,Cher,NWT,Hemisphere,FF,inst_deposit,
         ignore,plotpdf,preview,plotphoto,minphoto,maxphoto,
         Pbcol,inst_depositcol,
         modelcol,
         historic_d,historic_a,historic_n,historic_test,
         suppdescriptor,descriptor_lab,suppdescriptorcol,
         coring_yr,plot_Am,plot_Cs,plot_Pb,
         plot_Pb_inst_deposit,plot_CFCS_regression,
         varves, dmin,dmax,sedchange,
         min_yr, SML,stepout, mycex,
         archive_metadata)

.serac <- function(name, model,Cher,NWT,Hemisphere,FF,inst_deposit,
                   ignore,plotpdf,preview,plotphoto,minphoto,maxphoto,
                   Pbcol,inst_depositcol,
                   modelcol,
                   historic_d,historic_a,historic_n,historic_test,
                   suppdescriptor,descriptor_lab,suppdescriptorcol,
                   coring_yr,plot_Am,plot_Cs,plot_Pb,
                   plot_Pb_inst_deposit,plot_CFCS_regression,
                   varves, dmin,dmax,sedchange,
                   min_yr, SML,stepout, mycex,
                   archive_metadata) {

  # Calculate how long the function took to run
  old_time <- Sys.time() # get start time

  # load packages
  pkgTest("Hmisc")
  pkgTest("jpeg")
  pkgTest("TeachingDemos")

  # Archive metadata
  if(archive_metadata) core_metadata(name=name)

  # warn and stop if abnormal settings are provided
  # coring year is one of the two mandatory argument
  if(is.null(coring_yr))     stop("\n Warning, please enter the 'coring_yr'.\n\n")

  # serac support 2 changes in sedimentation rates maximum
  if(length(sedchange)>2)    stop("\n Warning, serac only support two changes in sedimentation rate. Please check the manual.\n\n", call.=FALSE)

  # Chernobyl fallouts are dated at different years depending on the hemisphere
  if(!is.null(Cher) && (Hemisphere=="NH"|Hemisphere=="SH")==FALSE)  stop("\n Warning, please select the hemisphere where your system is located.\n\n")

  # if the argument plotphoto is true, then a photo with the exact same name and the extension .jpg must be provided in the folder
  # min and max depth must be provided so the image is automatically cropped.
  if(plotphoto==TRUE) {
    if(!file.exists(paste(getwd(),"/Cores/",name,"/",name,".jpg", sep=""))) stop("\n Warning, you asked to include the photo of the core but it was not found in the repository.\n Check the name and the extension (must be .jpg).\n\n")
    if(is.null(minphoto) | (is.null(maxphoto)))                             stop("\n Warning, you need to indicate upper (minphoto) and lower (maxphoto) depth_avg of the core (mm).\n\n")
  }

  # depth must be provided by pair (upper and lower depths)
  if(length(historic_d)>=1) {
    if (length(historic_a) != length(historic_d)/2) stop("\n Warning, length(historic_a) != length(historic_d)/2 \n Read the help section.\n\n")
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
  if (length(intersect(grep("epth|EPTH",colnames(dt)),grep("top|bottom|min|max",colnames(dt),invert=TRUE)))>=1){
    dt$depth_avg <- dt[,intersect(grep("epth|EPTH",colnames(dt)),grep("top|bottom|min|max",colnames(dt),invert=TRUE))[1]]
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

  # 1.3. Fill in missing data ####
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

  # 1.4. Which keep ####
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
  rm(myvec)
  whichkeep <- whichkeep[!is.na(whichkeep)]
  # whichkeep tells you which "complete_core_depth" are not in an instantaneous deposit
  # Only these depths will be included when computing the inventory

  # Create the complete core depth 2 (for CRS calculations)
  # Everything but ignore
  complete_core_depth_2 <- complete_core_depth
  complete_core_depth_2[!complete_core_depth_2 %in% complete_core_depth_2[whichkeep]] <- NA

  # 1.5. Set some parameters ####
  if(!is.null(NWT) && Hemisphere=="NH") NWT_a <- 1963
  if(!is.null(NWT) && Hemisphere=="SH") NWT_a <- 1965
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
  mylegend <- NULL
  mypchlegend <- NULL
  myltylegend <- NULL
  mycollegend <- NULL

  # 1.6. Create the composite free depth_avg ####
  # Create the composite free depth_avg - step 1
  if(!exists("ignore")) ignore <- NULL
  if(SML>0) ignore <- c(ignore,dt$depth_avg[dt$depth_avg<=SML])

  dt$depth_avg_2 <- rep(NA,nrow(dt))
  for (i in 1:nrow(dt)) {
    if(!is.null(ignore)|max(inst_deposit)>0) {
      if(any(ignore==dt$depth_avg[i])) {
        dt$depth_avg_2[i] <- NA
      } else {dt$depth_avg_2[i] <- dt$depth_avg[i]}
    } else {dt$depth_avg_2[i] <- dt$depth_avg[i]}
  }

  # Create the composite free depth_avg - step 2: inst_deposit
  if (length(sedchange)==1 && sedchange == 0) sedchange_corr=max(dt$depth_avg,na.rm = T) else sedchange_corr=sedchange

  if(exists("inst_deposit")&&length(inst_deposit) > 1)
  {
    if(length(inst_deposit) %% 2 == 1) stop("\n Warning, inst_deposits need both upper and lower depth_avgs. Please check the manual.", call.=FALSE)
    inst_deposit_present = TRUE # Argument inst_deposit_present = FALSE decided elsewhere if no inst_deposit
    inst_deposit <- matrix(sort(inst_deposit), ncol=2, byrow=TRUE)
    for(i in 1:nrow(inst_deposit)) {
      if(length(dt$depth_avg[dt$depth_avg >= min(inst_deposit[i,]) & dt$depth_avg <= max(inst_deposit[i,])])>0) ignore <- c(ignore,dt$depth_avg[dt$depth_avg >= min(inst_deposit[i,]) & dt$depth_avg <= max(inst_deposit[i,])])
    }
    if(!is.null(ignore)) ignore <- ignore[!duplicated(ignore)]
    if(!is.null(ignore)) ignore <- ignore[order(ignore)]
    for (i in 1:nrow(dt)) {
      if(!is.null(ignore)&&max(ignore)>0|max(inst_deposit)>0) {
        if(!is.null(ignore)&&any(ignore==dt$depth_avg[i])) {
          dt$depth_avg_2[i] <- NA
        } else {dt$depth_avg_2[i] <- dt$depth_avg[i]}
      } else {dt$depth_avg_2[i] <- dt$depth_avg[i]}
    }
    d <- dt$depth_avg_2[!is.na(dt$depth_avg_2)]
    if(exists("historic_d")) dmax_corr=max(c(dt$depth_avg,historic_d),na.rm=T) else dmax_corr=max(c(dt$depth_avg),na.rm=T)
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

  # 1.7. Create separate datasets for different sedimentation rates) ####
  # PREPARATION FOR CFCS model
  # Here, we are looking to get three vectors:
  #     - One vector of the actual depths on the core we are trying to date (upper and lower limits of instantaneous deposit for instance)
  #     - The corrected version of this 1st vector, with instantaneous deposit removed
  #     - Depths that will be used to build visualise CFCS model
  d_for_CFCS <- unique(c(inst_deposit,max(dt$depth_avg[!is.na(dt$d)])))
  if(SML!=0) d_for_CFCS <- c(d_for_CFCS,SML)
  d_for_CFCS <- d_for_CFCS[order(d_for_CFCS)]

  # Final vector of depths that will be used for linear model (This is the 1/2 vector we're creating in this section)
  depth_avg_to_date <- c(0,d_for_CFCS,sedchange,dmax)
  depth_avg_to_date <- depth_avg_to_date[order(depth_avg_to_date)]

  # When there is an instantaneous deposit, the bottom depth is used.
  for (i in 2:length(d_for_CFCS)) {
    if(max(inst_deposit,na.rm = T)>0) {
      if(any(d_for_CFCS[i]==inst_deposit[,2])) {d_for_CFCS[i] <- d_for_CFCS[i-1]}  #else {d_for_CFCS[i] <- d_for_CFCS[i-1]+d_temp1[i]-d_temp1[i-1]}
    } else {d_for_CFCS[i] <- d_for_CFCS[i]}
  }
  depth_avg_to_date_corr <- c(0,d_for_CFCS,sedchange,dmax)
  depth_avg_to_date_corr <- depth_avg_to_date_corr[order(depth_avg_to_date_corr)]

  # Create a corrected vector that take in account instantaneous deposits, in the same way that before
  inst_deposit_corr2 <- inst_deposit
  if(exists("inst_deposit")&&length(inst_deposit) > 1) for(i in 1:nrow(inst_deposit))
  {
    depth_avg_to_date_corr[depth_avg_to_date_corr > min(inst_deposit_corr2[i,])]  <- depth_avg_to_date_corr[depth_avg_to_date_corr > min(inst_deposit_corr2[i,])] - (max(inst_deposit_corr2[i,]) - min(inst_deposit_corr2[i,]))
    if((1+i)<=nrow(inst_deposit)) inst_deposit_corr2[c(1+i):nrow(inst_deposit),] <- inst_deposit_corr2[c(1+i):nrow(inst_deposit),] - (max(inst_deposit_corr[i,])-min(inst_deposit_corr[i,]))
  }
  rm(inst_deposit_corr2)

  # Create sub-dataset if different sedimentation rates
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
      lm_sed1 <- lm(log(dt_sed1$Pbex[!is.na(dt_sed1$d)&dt_sed1$Pbex>0]) ~ dt_sed1$d[!is.na(dt_sed1$d)&dt_sed1$Pbex>0])
      sr_sed1 <- lambda/lm_sed1$coefficients[2]
      sr_sed1_err = sr_sed1*((lambda_err/lambda)^2+(summary(lm_sed1)$coefficients[2,2]/lm_sed1$coefficients[2])^2)^(0.5)

      # Print sed rate and error
      if (max(sedchange)==0) {
        cat(paste("\n Sedimentation rate (CFCS model): V= ",abs(round(sr_sed1,3)),"mm/yr, R2= ", round(summary(lm_sed1)$r.squared,4),"\n", sep=""))
        cat(paste("                          Error:     +/- ",abs(round(sr_sed1_err,3)),"mm/yr\n", sep=""))
      }

      if (max(sedchange)>0) {
        if(length(sedchange)==1) {
          # Linear model and V calculation (sedimentation rate)
          lm_sed2 <- lm(log(dt_sed2$Pbex[!is.na(dt_sed2$d)&dt_sed2$Pbex>0]) ~ dt_sed2$d[!is.na(dt_sed2$d)&dt_sed2$Pbex>0])
          sr_sed2 <- lambda/lm_sed2$coefficients[2]
          sr_sed2_err = sr_sed2*((lambda_err/lambda)^2+(summary(lm_sed2)$coefficients[2,2]/lm_sed2$coefficients[2])^2)^(0.5)

          # Print sed rate and error
          cat(paste("\n Sedimentation rate (CFCS model) ", SML,"-",sedchange[1],"mm: V= ",abs(round(sr_sed1,3)),"mm/yr, R2= ", round(summary(lm_sed1)$r.squared,4),"\n", sep=""))
          cat(paste("                          Error:     +/- ",abs(round(sr_sed1_err,3)),"mm/yr\n", sep=""))
          cat(paste("\n Sedimentation rate (CFCS model) ", sedchange[1],"mm-bottom",": V= ",abs(round(sr_sed2,3)),"mm/yr, R2= ", round(summary(lm_sed2)$r.squared,4),"\n", sep=""))
          cat(paste("                          Error:     +/- ",abs(round(sr_sed2_err,3)),"mm/yr\n", sep=""))
        }
        if(length(sedchange)==2) {
          ## 2nd change in sedimentation rate
          # Linear model and V calculation (sedimentation rate)
          lm_sed2 <- lm(log(dt_sed2$Pbex[!is.na(dt_sed2$d)&dt_sed2$Pbex>0]) ~ dt_sed2$d[!is.na(dt_sed2$d)&dt_sed2$Pbex>0])
          sr_sed2 <- lambda/lm_sed2$coefficients[2]
          sr_sed2_err = sr_sed2*((lambda_err/lambda)^2+(summary(lm_sed2)$coefficients[2,2]/lm_sed2$coefficients[2])^2)^(0.5)

          ## 3rd change in sedimentation rate
          # Linear model and V calculation (sedimentation rate)
          lm_sed3 <- lm(log(dt_sed3$Pbex[!is.na(dt_sed3$d)&dt_sed3$Pbex>0]) ~ dt_sed3$d[!is.na(dt_sed3$d)&dt_sed3$Pbex>0])
          sr_sed3 <- lambda/lm_sed3$coefficients[2]
          sr_sed3_err = sr_sed3*((lambda_err/lambda)^2+(summary(lm_sed3)$coefficients[2,2]/lm_sed3$coefficients[2])^2)^(0.5)

          # Print sed rate and error
          cat(paste("\n Sedimentation rate (CFCS model) ", SML,"-",sedchange[1],"mm: V= ",abs(round(sr_sed1,3)),"mm/yr, R2= ", round(summary(lm_sed1)$r.squared,4),"\n", sep=""))
          cat(paste("                          Error:     +/- ",abs(round(sr_sed1_err,3)),"mm/yr\n", sep=""))
          cat(paste("\n Sedimentation rate (CFCS model) ", sedchange[1],"-",sedchange[2],"mm: V= ",abs(round(sr_sed2,3)),"mm/yr, R2= ", round(summary(lm_sed2)$r.squared,4),"\n", sep=""))
          cat(paste("                          Error:     +/- ",abs(round(sr_sed2_err,3)),"mm/yr\n", sep=""))
          cat(paste("\n Sedimentation rate (CFCS model) ", sedchange[2],"mm-bottom",": V= ",abs(round(sr_sed3,3)),"mm/yr, R2= ", round(summary(lm_sed3)$r.squared,4),"\n", sep=""))
          cat(paste("                          Error:     +/- ",abs(round(sr_sed3_err,3)),"mm/yr\n", sep=""))

        }
      }
    }

    if(any(model=="CIC")) {
      # Calculation of age to be substracted
      Tm_CIC <- (1/lambda)*log(dt$Pbex[1]/dt$Pbex)
      # calculation age error: delta(tx)= 1/lambda*[(lambda_err*t)^2+(delta(A0)/A0)^2+(delta(Ax)/Ax)^2]^(0.5)
      # with Ax: activity at depth x; A0: initial activity
      # Two steps, 1 and 2
      # 1) replace any NA per 0 (just for this calculation, in temporary vectors)
      Pbex <- dt$Pbex
      Pbex[is.na(Pbex)] <- 0
      Pbex_er <- dt$Pbex_er
      Pbex_er[is.na(Pbex_er)] <- 0
      # 2) Actual error
      Tm_CIC_err <- (1/lambda)*((lambda_err*Tm_CIC)^2+(Pbex_er[1]/Pbex[1])^2+(Pbex_er/Pbex)^2)^(0.5)

      # calculation of Best Age and errors
      m_CIC <- coring_yr - Tm_CIC
      m_CIC_low <- m_CIC - Tm_CIC_err
      m_CIC_high <- m_CIC + Tm_CIC_err
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
    }
  }

  #### 3. CESIUM -----
  if(plot_Cs) {

    #Chernobyl
    if (exists("Cher")&&!is.null(Cher)) {
      peakCher <- -mean(Cher)
      dates <- c(dates,1986)
      dates_depth_avg <- c(dates_depth_avg,peakCher)
      err_dated_depth_avg <- cbind(err_dated_depth_avg,Cher)
    }
    #NWT
    if (exists("NWT")&&!is.null(NWT)) {
      peakNWT <- -mean(NWT)
      dates <- c(dates,NWT_a)
      dates_depth_avg <- c(dates_depth_avg,peakNWT)
      err_dated_depth_avg <- cbind(err_dated_depth_avg,NWT)
    }
    #First radionuclides fallout
    if (exists("FF")&&!is.null(FF)) {
      peakFF <- -mean(FF)
      dates <- c(dates,1955)
      dates_depth_avg <- c(dates_depth_avg,peakFF)
      err_dated_depth_avg <- cbind(err_dated_depth_avg,FF)
    }


  }

  #### 4. AGE DEPTH MODEL -----
  for (i in which(dt$depth_avg>=SML)){
    if (is.na(dt$d[i])) {
      if(i==1) dt$d[i] <- dt$depth_avg[1] else dt$d[i] <- dt$d[i-1]
    }
  }

  if(any(model=="CFCS")) {
    if(max(sedchange)>0) {
      age_break <- coring_yr-sedchange_corr[1]/abs(sr_sed1)
      age_break_low <- age_break-sedchange_corr[1]*abs(sr_sed1_err)/abs(sr_sed1)^2
      age_break_high <- age_break+sedchange_corr[1]*abs(sr_sed1_err)/abs(sr_sed1)^2
      cat(paste(" Approximation of age at change(s) in sedimentation rate:\n"))
      if(length(sedchange)==2) {
        age_break2 <- age_break-(sedchange_corr[2]-sedchange_corr[1])/abs(sr_sed2)
        age_break2_low <- age_break2-((sedchange_corr[2]-sedchange_corr[1]))*abs(sr_sed2_err)/abs(sr_sed2)^2
        age_break2_high <- age_break2+((sedchange_corr[2]-sedchange_corr[1]))*abs(sr_sed2_err)/abs(sr_sed2)^2
        cat(paste("     Best Age (1st change): ",abs(round(age_break,0))," (incertitude: ",abs(round(age_break_low,0)),"-",abs(round(age_break_high,0)),")\n",sep=""))
        cat(paste("     Best Age (2nd change): ",abs(round(age_break2,0))," (incertitude: ",abs(round(age_break2_low,0)),"-",abs(round(age_break2_high,0)),")\n\n",sep=""))
      } else {
        cat(paste("     Best Age: ",abs(round(age_break,0))," (incertitude: ",abs(round(age_break_low,0)),"-",abs(round(age_break_high,0)),")\n\n",sep=""))
      }
    }

    output_agemodel_CFCS <- matrix(rep(NA,length(depth_avg_to_date)*4), ncol=4)
    for(i in seq_along(depth_avg_to_date)){
      output_agemodel_CFCS[i,1] <- depth_avg_to_date[i]
      output_agemodel_CFCS[i,2] <- coring_yr-depth_avg_to_date_corr[i]/abs(sr_sed1)
      output_agemodel_CFCS[i,3] <- output_agemodel_CFCS[i,2]-depth_avg_to_date_corr[i]*abs(sr_sed1_err)/abs(sr_sed1)^2
      output_agemodel_CFCS[i,4] <- output_agemodel_CFCS[i,2]+depth_avg_to_date_corr[i]*abs(sr_sed1_err)/abs(sr_sed1)^2

      if(max(sedchange)>0 && depth_avg_to_date[i]>sedchange[1]) {
        output_agemodel_CFCS[i,2] <- age_break-(depth_avg_to_date_corr[i]-sedchange_corr[1])/abs(sr_sed2)
        output_agemodel_CFCS[i,3] <- output_agemodel_CFCS[i,2]-(depth_avg_to_date_corr[i])*abs(sr_sed2_err)/abs(sr_sed2)^2
        output_agemodel_CFCS[i,4] <- output_agemodel_CFCS[i,2]+(depth_avg_to_date_corr[i])*abs(sr_sed2_err)/abs(sr_sed2)^2
      }
      if(length(sedchange)>1 && depth_avg_to_date[i]>sedchange[2]) {
        output_agemodel_CFCS[i,2] <- age_break2-(depth_avg_to_date_corr[i]-sedchange_corr[2])/abs(sr_sed3)
        output_agemodel_CFCS[i,3] <- output_agemodel_CFCS[i,2]-(depth_avg_to_date_corr[i])*abs(sr_sed3_err)/abs(sr_sed3)^2
        output_agemodel_CFCS[i,4] <- output_agemodel_CFCS[i,2]+(depth_avg_to_date_corr[i])*abs(sr_sed3_err)/abs(sr_sed3)^2
      }
    }
    output_agemodel_CFCS <- as.data.frame(output_agemodel_CFCS)
    colnames(output_agemodel_CFCS) <- c("depth", "BestAD", "MinAD", "MaxAD")
    output_agemodel_CFCS <- output_agemodel_CFCS[!duplicated(output_agemodel_CFCS[,1]),]
    output_agemodel_CFCS_inter <- as.data.frame(seq(0,max(output_agemodel_CFCS$depth,na.rm = T),by=stepout))
    if (length(historic_d)>=1 && any(is.na(historic_a))) {
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
        cat(paste("     The instantaneous deposit at ",inst_deposit[i,1],"-", inst_deposit[i,2]," mm has an estimated range of: ", round(output_agemodel_CFCS[which(output_agemodel_CFCS[,1]==inst_deposit[i,1]),3]),"-",round(output_agemodel_CFCS[which(output_agemodel_CFCS[,1]==inst_deposit[i,1]),4]),".\n\n",sep=""))
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

    # Parameters for legend
    mylegend <- c(mylegend, "CFCS")
    mypchlegend <- c(mypchlegend,NA)
    myltylegend <- c(myltylegend,1)
    mycollegend <- c(mycollegend,modelcol[1])
  }

  if(any(model=="CIC")) {
    output_agemodel_CIC <- as.data.frame(matrix(c(c(0,dt$depth_avg),c(coring_yr,m_CIC),c(coring_yr,m_CIC_low),c(coring_yr,m_CIC_high)), byrow = F, ncol=4))
    colnames(output_agemodel_CIC) <- c("depth_avg", "BestAD_CIC", "MinAD_CIC", "MaxAD_CIC")
    output_agemodel_CIC_inter <- as.data.frame(seq(0,max(output_agemodel_CIC$depth_avg,na.rm = T),stepout))
    output_agemodel_CIC_inter <- cbind(output_agemodel_CIC_inter,approx(x= output_agemodel_CIC$depth_avg, output_agemodel_CIC$BestAD_CIC, xout= seq(0,max(output_agemodel_CIC$depth_avg,na.rm = T),stepout))$y)
    output_agemodel_CIC_inter <- cbind(output_agemodel_CIC_inter,approx(x= output_agemodel_CIC$depth_avg, output_agemodel_CIC$MinAD_CIC, xout= seq(0,max(output_agemodel_CIC$depth_avg,na.rm = T),stepout))$y)
    output_agemodel_CIC_inter <- cbind(output_agemodel_CIC_inter,approx(x= output_agemodel_CIC$depth_avg, output_agemodel_CIC$MaxAD_CIC, xout= seq(0,max(output_agemodel_CIC$depth_avg,na.rm = T),stepout))$y)
    colnames(output_agemodel_CIC_inter) <- c("depth_avg", "BestAD_CIC", "MinAD_CIC", "MaxAD_CIC")
    write.table(x = output_agemodel_CIC[order(output_agemodel_CIC$depth_avg, decreasing = F),], file = paste(getwd(),"/Cores/",name,"/",name,"_CIC.txt",sep = ""),col.names = T, row.names = F)
    write.table(x = output_agemodel_CIC_inter[order(output_agemodel_CIC_inter$depth_avg, decreasing = F),], file = paste(getwd(),"/Cores/",name,"/",name,"_CIC_interpolation.txt",sep = ""),col.names = T, row.names = F)

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
    err_dated_depth_avg <- matrix(-err_dated_depth_avg[!is.na(err_dated_depth_avg)],nrow=2,byrow = F)
  }


  # Various parameter to print (e.g. Inventories)
  # Inventory Lead
  if(length(grep("Pb",x = colnames(dt)))>1 & length(grep("density",x = colnames(dt)))>=1) {
    # We multiply the value by 10 because we ask for the depth in mm, and the density in g/cm3
    cat(paste(" Inventory (Lead): ",round(Inventory_CRS[1],3), " Bq/m2 (range: ",round(Inventory_CRS[1]-Inventory_CRS_error[1]), "-",round(Inventory_CRS[1]+Inventory_CRS_error[1])," Bq/m2)\n", sep=""))
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
    if (max(sedchange)>0) {
      metadata <- rbind(metadata,
                        c(paste("Sedimentation rate (CFCS model) ", SML,"-",sedchange[1],"mm",sep=""),paste("V= ",abs(round(sr_sed1,3)),"mm/yr, R2= ", round(summary(lm_sed1)$r.squared,4),", Error +/- ",abs(round(sr_sed1_err,3)),"mm/yr", sep="")))
    } else {
      metadata <- rbind(metadata,
                        c("Sedimentation rate (CFCS model)",paste("V= ",abs(round(sr_sed1,3)),"mm/yr, R2= ", round(summary(lm_sed1)$r.squared,4),", Error +/- ",abs(round(sr_sed1_err,3)),"mm/yr", sep="")))
    }
    if (max(sedchange)>0) {
      if(length(sedchange)==1) {
        metadata <- rbind(metadata,
                          c(paste("Sedimentation rate (CFCS model) ", sedchange[1],"mm-bottom", sep=""),paste("V= ",abs(round(sr_sed2,3)),"mm/yr, R2= ", round(summary(lm_sed2)$r.squared,4), ", Error +/- ",abs(round(sr_sed1_err,3)),"mm/yr",sep="")))
      }
      if(length(sedchange)==2) {
        metadata <- rbind(metadata,
                          c(paste("Sedimentation rate (CFCS model) ", sedchange[1],"-",sedchange[2],"mm", sep=""),paste("V= ",abs(round(sr_sed2,3)),"mm/yr, R2= ", round(summary(lm_sed2)$r.squared,4), ", Error +/- ",abs(round(sr_sed1_err,3)),"mm/yr",sep="")),
                          c(paste("Sedimentation rate (CFCS model) ", sedchange[2],"mm-bottom", sep=""),paste("V= ",abs(round(sr_sed3,3)),"mm/yr, R2= ", round(summary(lm_sed3)$r.squared,4), ", Error +/- ",abs(round(sr_sed1_err,3)),"mm/yr",sep="")))
      }
    }
  }

  # Add in the output data file the ages estimations for changes in sed rate
  if(max(sedchange)>0) {
    if(length(sedchange)==2) {
      metadata <- rbind(metadata,
                        c("Best Age (1st change)",paste(abs(round(age_break,0))," (incertitude: ",abs(round(age_break_low,0)),"-",abs(round(age_break_high,0)),")",sep="")))
      metadata <- rbind(metadata,
                        c("Best Age (2nd change)",paste(abs(round(age_break2,0))," (incertitude: ",abs(round(age_break2_low,0)),"-",abs(round(age_break2_high,0)),")",sep="")))
    } else {
      metadata <- rbind(metadata,
                        c("Best Age",paste(abs(round(age_break,0))," (incertitude: ",abs(round(age_break_low,0)),"-",abs(round(age_break_high,0)),")",sep="")))
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

  # Add the code that was used
  # First save the history to folder
  savehistory(file = "myhistory.Rhistory")

  # The code may extend on several lines (true for RStudio users at least)
  # This loop find the beginning of the function, serac
  whichline=NULL
  for (i in 1:30) {
    if (length(grep(pattern = "serac", rev(readLines(con = "myhistory.Rhistory"))[i]))>0) whichline <- c(whichline,i)
  }
  whichline <- min(whichline, na.rm=T)
  mycode=NULL
  for (i in whichline:1) {
    mycode <- paste(mycode, rev(readLines(con = "myhistory.Rhistory"))[i], sep="")
  }
  # Remove the history
  if (file.exists("myhistory.Rhistory")) file.remove("myhistory.Rhistory")

  metadata <- rbind(metadata,
                    c("",""),
                    c("code",mycode)) # Add line to metadata

  # Write final output
  write.table(x = metadata, file = paste(getwd(),"/Cores/",name,"/",name,"_Metadata_",Sys.Date(),".txt",sep = ""),col.names = F, row.names = F, sep = "\t")


  #### 6. FINAL PLOT -----
  run <- NULL
  if(preview) run <- c(run,1)
  if(plotpdf) run <- c(run,2)
  if(!is.null(run)) {
    for(whichrun in c(run)) {
      if (whichrun==1) {plot(0,0, axes=F, xlab="",ylab="",pch=NA); dev.off(); dev.new()}
      # size for output plot
      cex_1=.8*mycex
      cex_2=1.1*mycex
      cex_4=1.1*mycex #Writing within plot of sedimentation rate

      # set plot window
      mylayout <- NULL # mylayout vector will set the width of the different windows within the plot
      if(plotphoto) mylayout <- c(mylayout,.2)
      if(suppdescriptor) mylayout <- c(mylayout,1)
      if(plot_Pb) mylayout <- c(mylayout,1)
      if(plot_Pb_inst_deposit) mylayout <- c(mylayout,1.3)
      if(plot_Cs) mylayout <- c(mylayout,1.3)
      mylayout <- c(mylayout,1.6)

      mylayout[1] <- mylayout[1]+.3 #Add margin to the right windows to include the depth_avg scale
      nwindows <- length(mylayout)

      if(whichrun==1) {plot.new()}
      if(whichrun==2) pdf(paste(getwd(),"/Cores/",name,"/",name,".pdf",sep=""),width = 1.2+1.8*nwindows, height = 5,family = "Helvetica")

      layout(matrix(c(1:nwindows),1,nwindows), widths = mylayout)
      # 6.1. Add core photo ####
      if(plotphoto) {
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
      if(suppdescriptor) {
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
      if(plot_Pb) {
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

        # create the flexible axis ticks
        power <- 1
        for (i in 0:3) power <- c(power, 2:10*10^c(i))
        label <- power
        label[grep("1",power,invert = TRUE)] <- NA
        axis(3, at=log(power[power<exp(myxlim_max)]), labels=label[power<exp(myxlim_max)])
        # axis label
        mtext(text = bquote(~""^210*Pb[ex]*" (mBq/g)"), side = 3, line=2.2, cex=cex_1)

        # Add 'ignore' values
        par(new=T)
        with (data=dt[is.na(dt$depth_avg_2),]
              , expr = errbar(log(Pbex),-depth_avg,c(-depth_avg+thickness/2),c(-depth_avg-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim, col=grey(.65), errbar.col = grey(.65), cex=.8)
        )
        for (i in which(is.na(dt$depth_avg_2))) {
          lines(c(log(dt$Pbex[i]+dt$Pbex_er[i]),log(dt$Pbex[i]-dt$Pbex_er[i])),
                rep(-dt$depth_avg[i],2), type="o", pch="|", cex=.5, col=grey(.65))
        }
      }

      # 6.3.b Plot 210Pb without inst_deposits ####
      if(plot_Pb_inst_deposit) {
        if (inst_deposit_present) {
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


          for (i in 1:length(which(dt_sed1$d>0))) {
            lines(c(log(dt_sed1$Pbex[which(dt_sed1$d>0)][i]+dt_sed1$Pbex_er[which(dt_sed1$d>0)][i]),log(dt_sed1$Pbex[which(dt_sed1$d>0)][i]-dt_sed1$Pbex_er[which(dt_sed1$d>0)][i])),
                  rep(-dt_sed1$d[which(dt_sed1$d>0)][i],2), type="o", pch="|", cex=.5, col=Pbcol[1])
          }

          if (max(sedchange)>0) {
            par(new=T)
            with (
              data=dt_sed2
              , expr = errbar(log(Pbex),-d,c(-d+thickness/2),c(-d-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim, col=Pbcol[2], errbar.col = Pbcol[2])
            )
            for (i in which(dt_sed2$d>0)) {
              lines(c(log(dt_sed2$Pbex[i]+dt_sed2$Pbex_er[i]),log(dt_sed2$Pbex[i]-dt_sed2$Pbex_er[i])),
                    rep(-dt_sed2$d[i],2), type="o", pch="|", cex=.5, col=Pbcol[2])
            }
            if(length(sedchange)==2) {
              par(new=T)
              with (
                data=dt_sed3
                , expr = errbar(log(Pbex),-d,c(-d+thickness/2),c(-d-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,xlim=c(log(1),log(mround(max(dt$Pbex,na.rm=T),1000))),ylim=myylim, col=Pbcol[3], errbar.col = Pbcol[3])
              )
              for (i in which(dt_sed3$d>0)) {
                lines(c(log(dt_sed3$Pbex[i]+dt_sed3$Pbex_er[i]),log(dt_sed3$Pbex[i]-dt_sed3$Pbex_er[i])),
                      rep(-dt_sed3$d[i],2), type="o", pch="|", cex=.5, col=Pbcol[3])
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
        }
      }


      # 6.3.c Plot 210Pb model on top of 210Pb measurements ####
      if(any(model=="CFCS")) {
        # Warning, you shouldn't plot the linear regression on the point without instantaneous deposit while you mentioned there were some
        if(whichrun==1&exists("inst_deposit")&length(inst_deposit)>1&min(inst_deposit)<max(dt$depth_avg)&plot_Pb_inst_deposit==F&plot_CFCS_regression==F) {
          message("\n Warning. It is not possible to visualise the linear regression calculated from\n the CFCS model on the 210Pbex curve without instantaneous deposits, on the\n initial depth. Add the argument plot_Pb_inst_deposit=TRUE to visualise the\n regression line on 210Pbex profile corrected from instananeous events.\n Keep in mind the regression line won't necesseraly match the points.\n\n")
        }
        # If you decide to do it anyway, this message will display
        if(whichrun==1&exists("inst_deposit")&length(inst_deposit)>1&min(inst_deposit)<max(dt$depth_avg)&plot_Pb_inst_deposit==F&plot_CFCS_regression==T) {
          message("\n Warning. It seems you requested to visualise the linear regression calculated\n from the CFCS model on the 210Pbex curve without instantaneous deposits,\n while you specified there were some.\n Please bear in mind the linear regression does not correspond to the points. \n Turn plot_Pb_inst_deposit to TRUE or plot_CFCS_regression to FALSE.\n\n")
        }
        if(plot_Pb & plot_CFCS_regression | plot_Pb_inst_deposit & plot_CFCS_regression) {

          # 210Pb linear model needs to stop in case deepest depth have been set to be 'ignored'
          # (bug observed on 2018-10-26 by PS on CA08 - answer in email RB to PS on 2018-11-13)
          # When plotting the linear model, we'll use either of these two following lines (after toplm definition)
          # One other special case regarding the next lines (before if()): when the surface samples are ignored, the regression line must not extend to the surface
          # We then define the value 'top linear model' to decrease the number of conditons in the code...
          if (!is.null(ignore)) toplm <- max(c(SML,min(dt$depth_avg[!dt$depth_avg %in% ignore], na.rm = T)), na.rm = T) else toplm <- SML
          if (is.null(ignore) || !is.null(is.null(ignore)) && sedchange_corr[1] > max(ignore)) lines(c(-toplm,-sedchange_corr[1])~ c(lm_sed1$coefficients[1]+toplm*lm_sed1$coefficients[2],lm_sed1$coefficients[1]+sedchange_corr[1]*lm_sed1$coefficients[2]), col=Pbcol[1], lwd=2)
          if (!is.null(ignore) && sedchange_corr[1] <= max(ignore)) lines(c(-toplm,-max(dt$d[!dt$depth_avg %in% ignore & dt$d<sedchange_corr[1]],na.rm=T))~ c(lm_sed1$coefficients[1]+toplm*lm_sed1$coefficients[2],lm_sed1$coefficients[1]+max(dt$d[!dt$depth_avg %in% ignore & dt$d<sedchange_corr[1]],na.rm=T)*lm_sed1$coefficients[2]), col=Pbcol[1], lwd=2)

          d_legend <- mean(c(min(dt_sed1$d,na.rm = T),max(dt_sed1$d,na.rm = T)))*.8
          shadowtext(x = 0,y = -d_legend-.06*dmax,labels = bquote(r^2 ~ "=" ~ .(round(summary(lm_sed1)$r.squared,4))), pos = 4, col=Pbcol[1], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4)
          shadowtext(x = 0,y = -d_legend,labels = bquote(V ~ "=" ~ .(abs(round(sr_sed1,3))) ~ mm.yr^-1), pos = 4, col=Pbcol[1], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4)

          if (max(sedchange)>0) {
            if(length(sedchange)==1) {
              if (is.null(ignore) || !is.null(is.null(ignore)) && max(dt_sed2$depth_avg, na.rm = T) > max(ignore)) lines(c(-sedchange_corr[1],-max(dt_sed2$d, na.rm = T))~ c(lm_sed2$coefficients[1]+sedchange_corr[1]*lm_sed2$coefficients[2],lm_sed2$coefficients[1]+max(dt_sed2$d, na.rm = T)*lm_sed2$coefficients[2]), lwd=2, col=Pbcol[2])
              if (!is.null(ignore) && max(dt_sed2$depth_avg, na.rm = T) <= max(ignore)) lines(c(-sedchange_corr[1],-max(dt$d[!dt$depth_avg %in% ignore],na.rm=T))~ c(lm_sed2$coefficients[1]+sedchange_corr[1]*lm_sed2$coefficients[2],lm_sed2$coefficients[1]+max(dt$d[!dt$depth_avg %in% ignore],na.rm=T)*lm_sed2$coefficients[2]), lwd=2, col=Pbcol[2])
              d_legend <- mean(c(min(dt_sed2$d,na.rm = T),max(dt_sed2$d,na.rm = T)))*.8
              shadowtext(x = 0,y = -d_legend-.06*dmax,labels = bquote(r^2 ~ "=" ~ .(round(summary(lm_sed2)$r.squared,4))), pos = 4, col=Pbcol[2], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4)
              shadowtext(x = 0,y = -d_legend,labels = bquote(V ~ "=" ~ .(abs(round(sr_sed2,3))) ~ mm.yr^-1), pos = 4, col=Pbcol[2], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=cex_4)
            }
            if(length(sedchange)==2) {
              lines(c(-sedchange_corr[1],-max(dt_sed2$d, na.rm = T))~ c(lm_sed2$coefficients[1]+sedchange_corr[1]*lm_sed2$coefficients[2],lm_sed2$coefficients[1]+max(dt_sed2$d, na.rm = T)*lm_sed2$coefficients[2]), lwd=2, col=Pbcol[2])
              d_legend <- mean(c(min(dt_sed2$d,na.rm = T),max(dt_sed2$d,na.rm = T)))*.8
              shadowtext(x = 0,y = -d_legend-.06*dmax,labels = bquote(r^2 ~ "=" ~ .(round(summary(lm_sed2)$r.squared,4))), pos = 4, col=Pbcol[2], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1)
              shadowtext(x = 0,y = -d_legend,labels = bquote(V ~ "=" ~ .(abs(round(sr_sed2,3))) ~ mm.yr^-1), pos = 4, col=Pbcol[2], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1)

              if (is.null(ignore) || !is.null(is.null(ignore)) && max(dt_sed3$depth_avg, na.rm = T) > max(ignore))  lines(c(-sedchange_corr[2],-max(dt_sed3$d, na.rm = T))~ c(lm_sed3$coefficients[1]+sedchange_corr[2]*lm_sed3$coefficients[2],lm_sed3$coefficients[1]+max(dt_sed3$d, na.rm = T)*lm_sed3$coefficients[2]), lwd=2, col=Pbcol[3])
              if (!is.null(ignore) && max(dt_sed3$depth_avg, na.rm = T) <= max(ignore)) lines(c(-sedchange_corr[2],-max(dt$d[!dt$depth_avg %in% ignore],na.rm=T))~ c(lm_sed3$coefficients[1]+sedchange_corr[2]*lm_sed3$coefficients[2],lm_sed3$coefficients[1]+max(dt$d[!dt$depth_avg %in% ignore],na.rm=T)*lm_sed3$coefficients[2]), lwd=2, col=Pbcol[2])
              d_legend <- mean(c(min(dt_sed3$d,na.rm = T),max(dt_sed3$d,na.rm = T)))*.8
              shadowtext(x = 0,y = -d_legend-.06*dmax,labels = bquote(r^2 ~ "=" ~ .(round(summary(lm_sed3)$r.squared,4))), pos = 4, col=Pbcol[3], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1)
              shadowtext(x = 0,y = -d_legend,labels = bquote(V ~ "=" ~ .(abs(round(sr_sed3,3))) ~ mm.yr^-1), pos = 4, col=Pbcol[3], bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1)
            }
          }
        }
      }


      # 6.4. 137Cs ####
      if(plot_Cs) {
        if(plotphoto || suppdescriptor || plot_Pb || plot_Pb_inst_deposit) par(mar=c(4.1,1.1,4.1,1.1)) else par(mar=c(4.1,4.1,4.1,1.1))
        myxlim_max <- max(dt$Cs, na.rm=T)*1.2+max(dt$Cs_er, na.rm = T)
        myxlim_min <- min(dt$Cs, na.rm=T)-max(dt$Cs_er, na.rm = T)
        with (
          data=dt[dt$depth_avg<SML,]
          , expr = errbar(Cs,-depth_avg,c(-depth_avg+thickness/2),c(-depth_avg-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,ylim=myylim, xlim=c(myxlim_min,myxlim_max),col=grey(.65), errbar.col = grey(.65), cex=.8)
        )

        par(xpd=TRUE)
        if(inst_deposit_present) for (i in 1:nrow(inst_deposit)) rect(xleft = -2000, ybottom = -inst_deposit[i,2], xright = max(dt$Cs,na.rm=T)*1.5, ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
        if(SML>0) rect(xleft = -2000, ybottom = -SML, xright = max(dt$Cs,na.rm=T)*1.5, ytop = 0, col=grey(0.97), border=NA)
        par(xpd=FALSE)

        par(new=T)
        with (
          data=dt[dt$depth_avg<SML,]
          , expr = errbar(Cs,-depth_avg,c(-depth_avg+thickness/2),c(-depth_avg-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,ylim=myylim, xlim=c(myxlim_min,myxlim_max),col=grey(.65), errbar.col = grey(.65), cex=.8)
        )


        for (i in which(dt$Cs>=0 & !is.na(dt$Cs_er) & dt$depth_avg<SML)) {
          lines(c(dt$Cs[i]+dt$Cs_er[i],dt$Cs[i]-dt$Cs_er[i]),
                rep(-dt$depth_avg[i],2), type="o", pch="|", cex=.5, col=grey(.65))
        }
        lines(dt$Cs[which(dt$depth_avg<SML+2)],-dt$depth_avg[which(dt$depth_avg<SML+2)], col=grey(.65), lwd=.5)

        par(new=T)
        with (
          data=dt[dt$depth_avg>=SML,]
          , expr = errbar(Cs,-depth_avg,c(-depth_avg+thickness/2),c(-depth_avg-thickness/2), pch=16, cap=.01, xlab="",ylab="", axes=F,ylim=myylim, xlim=c(myxlim_min,myxlim_max), col=Pbcol[1], errbar.col = Pbcol[1], cex=.8)
        )
        for (i in which(dt$Cs>0 & !is.na(dt$Cs_er) & dt$depth_avg>=SML)) {
          lines(c(dt$Cs[i]+dt$Cs_er[i],dt$Cs[i]-dt$Cs_er[i]),
                rep(-dt$depth_avg[i],2), type="o", pch="|", cex=.5, col=Pbcol[1])
        }
        lines(dt$Cs[which(dt$depth_avg>=SML)],-dt$depth_avg[which(dt$depth_avg>=SML)])

        axis(3,  cex.axis=cex_2)
        mtext(text = bquote(~""^137*"Cs (mBq/g)"), side = 3, line=2.2, cex=cex_1)

        # Add text
        par(xpd=TRUE)
        #Chernobyl
        if (exists("Cher")&&!is.null(Cher)) {
          lines(rep(max(dt$Cs[dt$depth_avg>min(Cher-10) & dt$depth_avg<max(Cher+10)],na.rm = T)*1.1,2),c(-Cher[1],-Cher[2]), lwd=2)
          shadowtext(max(dt$Cs[dt$depth_avg>min(Cher-10) & dt$depth_avg<max(Cher+10)],na.rm = T)+0.1*max(dt$Cs,na.rm=T),-(min(Cher)),
                     labels = c("C 1986"), pos = 3,col="black", bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=mycex)
          lines(c(max(dt$Cs[dt$depth_avg>min(Cher-10) & dt$depth_avg<max(Cher+10)],na.rm = T)*1.1, max(dt$Cs,na.rm = T)*2),rep(peakCher,2), lty=2)
        }
        #NWT
        if (exists("NWT")&&!is.null(NWT)) {
          lines(rep(max(dt$Cs[dt$depth_avg>min(NWT-10) & dt$depth_avg<max(NWT+10)],na.rm = T)*1.1,2),c(-NWT[1],-NWT[2]), lwd=2)
          shadowtext(max(dt$Cs,na.rm = T)+0.1*max(dt$Cs,na.rm=T),-(min(NWT)),
                     labels = paste("NWT",NWT_a,sep=" "), pos = 3, col="black",bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=mycex)
          lines(c(max(dt$Cs[dt$depth_avg>min(NWT-10) & dt$depth_avg<max(NWT+10)],na.rm = T)*1.1, max(dt$Cs,na.rm = T)*2),rep(peakNWT,2), lty=2)
        }
        #First radionuclides fallout
        if (exists("FF")&&!is.null(FF)) {
          lines(rep(max(dt$Cs[dt$depth_avg>min(FF-10) & dt$depth_avg<max(FF+10)],na.rm = T)*1.1,2),c(-FF[1],-FF[2]), lwd=2)
          shadowtext(max(dt$Cs,na.rm = T)+0.1*max(dt$Cs,na.rm=T),-(max(FF)),
                     labels = c("FF 1955"), pos = 1, col="black",bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, cex=mycex)
          lines(c(max(dt$Cs[dt$depth_avg>min(FF-10) & dt$depth_avg<max(FF+10)],na.rm = T)*1.1, max(dt$Cs,na.rm = T)*2),rep(peakFF,2), lty=2)
        }
        par(xpd=FALSE)

        # 6.5. 241Am ####
        if (plot_Am) {
          legend("bottomright", legend = c("Cesium", "Americium"), bty="n", pch=c(16,1), cex=mycex)
          par(new=T, mar=c(4.1,1.1,4.1,6.1))
          myxlim_max <- max(dt$Am, na.rm=T)*1.2+max(dt$Am_er, na.rm=T)
          myxlim_min <- min(dt$Am, na.rm=T)-max(dt$Am_er, na.rm=T)
          with (
            data=dt[which(!is.na(dt$Am)&dt$Am>0&dt$depth_avg>SML),]
            , expr = errbar(Am,-depth_avg,c(-depth_avg+thickness/2),c(-depth_avg-thickness/2), pch=1, cap=.01, xlab="",ylab="", axes=F,ylim=myylim, xlim=c(myxlim_min,myxlim_max),col=Pbcol[1], errbar.col = Pbcol[1],cex=.8)
          )
          axis(1, cex.axis=cex_2)
          mtext(text = bquote(~""^241*"Am (mBq/g)"), side = 1, line=2.4, cex=cex_1)
          for (i in which(dt$Am>0 & !is.na(dt$Am_er) & dt$depth_avg>SML)) {
            lines(c(dt$Am[i]+dt$Am_er[i],dt$Am[i]-dt$Am_er[i]),
                  rep(-dt$depth_avg[i],2), type="o", pch="|", cex=.5, col=Pbcol[1])
          }
          points(dt$Am[which(dt$Am>0&dt$depth_avg>SML)],-dt$depth_avg[which(dt$Am>0&dt$depth_avg>SML)], pch=20, col="white")
          points(dt$Am[which(dt$Am>0&dt$depth_avg>SML)],-dt$depth_avg[which(dt$Am>0&dt$depth_avg>SML)])

        }
      }

      # 6.5.a plot Age Model ####
      par(mar=c(4.1,1.1,4.1,4.1))
      plot(c(-min_yr,-mround(coring_yr,10)),c(-dmin, -dmax), xlab="",ylab="", axes=F, type="n",ylim=myylim)

      # Plot the 'historic_test' argument i.e. know dates we want to add
      if (!is.null(historic_test)){
        abline(v = -historic_test, col = adjustcolor("grey", alpha.f = .3), lwd=2)
        shadowtext(-historic_test,rep(-dmin, length(historic_test)),
                   labels = as.character(historic_test), col="black",bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1,cex = .9*mycex)
      }


      if(inst_deposit_present)  {
        for (i in 1:nrow(inst_deposit)) rect(xleft = -coring_yr, ybottom = -inst_deposit[i,2], xright = -min_yr+20, ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
      }

      par(xpd=T)
      if(inst_deposit_present)  {
        for (i in 1:nrow(inst_deposit)) rect(xleft = -2100, ybottom = -inst_deposit[i,2], xright = -coring_yr, ytop = -inst_deposit[i,1],col=inst_depositcol, border=inst_depositcol, lwd=.4)
      }
      par(xpd=F)


      axis(4, at = seq(min(myylim),0,by=10), NA, cex.axis=cex_2, lwd=.3)
      axis(4, at = -(pretty(seq(dmin,dmax,5))), labels=pretty(seq(dmin,dmax,5)), cex.axis=cex_2)
      mtext(text = "Depth (mm)", side = 4, line=2.2, cex=cex_1)
      axis(3, at = seq(-mround(coring_yr,10),-min_yr+20,20), labels = seq(mround(coring_yr,10),min_yr-20,-20),cex.axis=cex_2)
      mtext(text = "Year (C.E.)", side = 3, line=2.2, cex=cex_1)

      if(any(model=="CFCS")) {
        lines(-output_agemodel_CFCS$BestAD,-output_agemodel_CFCS$depth, col=modelcol[1],lty=2,lwd=.5)
        pol_x <- c(-output_agemodel_CFCS$MinAD, rev(-output_agemodel_CFCS$MaxAD))
        pol_y <- c(-output_agemodel_CFCS$depth, rev(-output_agemodel_CFCS$depth))
        polygon(x=pol_x, y = pol_y, col=adjustcolor(modelcol[1], alpha.f=0.2), border=NA)
        lines(-output_agemodel_CFCS$BestAD[output_agemodel_CFCS$depth>=SML&output_agemodel_CFCS$depth<=max(dt$depth_avg[!is.na(dt$d)])],-output_agemodel_CFCS$depth[output_agemodel_CFCS$depth>=SML&output_agemodel_CFCS$depth<=max(dt$depth_avg[!is.na(dt$d)])], col=modelcol[1])
      }

      if(any(model=="CIC")) {
        pol_x <- c(-m_CIC_low[!is.na(dt$depth_avg_2)], rev(-m_CIC_high[!is.na(dt$depth_avg_2)]))
        pol_y <- c(-dt$depth_avg[!is.na(dt$depth_avg_2)], rev(-dt$depth_avg[!is.na(dt$depth_avg_2)]))
        polygon(x=pol_x, y = pol_y, col=adjustcolor(modelcol[2], alpha.f=0.2), border=NA)
        lines(-m_CIC[!is.na(dt$depth_avg_2)],-dt$depth_avg[!is.na(dt$depth_avg_2)], col=modelcol[2])
      }

      if(any(model=="CRS")) {
        if(exists("inst_deposit")&&length(inst_deposit)>1) {
          lines(-new_x_CRS,-new_y_CRS, col=modelcol[3],lty=2,lwd=.5)
          pol_x <- c(-new_x_CRS_low, rev(-new_x_CRS_high))
          pol_y <- c(-new_y_CRS, rev(-new_y_CRS))
          polygon(x=pol_x, y = pol_y, col=adjustcolor(modelcol[3], alpha.f=0.2), border=NA)
          lines(-new_x_CRS[new_y_CRS>=SML&new_y_CRS<=max(dt$depth_avg)],-new_y_CRS[new_y_CRS>=SML&new_y_CRS<=max(dt$depth_avg)], col=modelcol[3])
        } else {
          lines(-m_CRS,-complete_core_depth[whichkeep], col=modelcol[3],lty=2,lwd=.5)
          pol_x <- c(-m_CRS_low, rev(-m_CRS_high))
          pol_y <- c(-complete_core_depth[whichkeep], rev(-complete_core_depth[whichkeep]))
          polygon(x=pol_x, y = pol_y, col=adjustcolor(modelcol[3], alpha.f=0.2), border=NA)
          lines(-m_CRS[complete_core_depth[whichkeep]>=SML&complete_core_depth[whichkeep]<=max(dt$depth_avg)],-complete_core_depth[whichkeep][complete_core_depth[whichkeep]>=SML&complete_core_depth[whichkeep]<=max(dt$depth_avg)], col=modelcol[3])
        }

      }

      if(varves) {
        points(-varve$Age,-varve$depth_avg, pch=4)
      }

      if(length(model)>1|varves) {
        legend("bottomleft", legend = mylegend, pch = mypchlegend,lty = myltylegend, col = mycollegend, bty='n', cex=mycex)
      }


      if(exists("Cher") | exists("NWT") | exists("FF")) {
        err_dated_depth_avg <- matrix(err_dated_depth_avg[!is.na(err_dated_depth_avg)],nrow=2,byrow = F)
        err_dated_depth_avg <- -abs(err_dated_depth_avg)

        ##
        if(!is.null(dates)) {
          for (i in 1:length(dates)) {
            lines(rep(-dates[i],2), c(err_dated_depth_avg[1,i], err_dated_depth_avg[2,i]), type="o", pch="_", col="black")
          }
          points(-dates,dates_depth_avg, pch=16, cex=.8)

          par(xpd=T)
          for (i in 1:length(dates)) {
            lines(c(-2100,-dates[i]), rep(dates_depth_avg[i],2), lty=2)
          }
          par(xpd=F)
        }
      }

      # 6.5.b Plot historic event on age model ####
      if(length(historic_d)>=1) {
        historic_d_dt <- matrix(abs(historic_d), ncol = 2, byrow = T)
        for (i in seq_along(historic_a)) {
          if (!is.na(historic_a[i])) {
            lines(x = c(rep(-historic_a[i],2)), y = c(20,-(historic_d_dt[i,1])),lty=3)
            par(xpd=T)
            if (!is.null(historic_n[i]) && !is.na(historic_n[i])) {
              shadowtext(-(min_yr+(coring_yr-min_yr)*.15),-mean(historic_d_dt[i,], na.rm=T),
                         labels = as.character(historic_n[i]), col="black",bg = "white", theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1,cex = 1*mycex)
            }
            par(xpd=F)
          }
        }
      }


      if(whichrun==2) dev.off()
    }
  }

  # print elapsed time
  new_time <- Sys.time() - old_time # calculate difference
  # print in nice format
  cat("\n\n ________________________\n")
  cat(paste("\n The calculation took ", round(new_time,3), " seconds.", sep=""))
}

