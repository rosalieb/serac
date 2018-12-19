\name{serac}
\alias{serac}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function for age-depth modelling
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
serac(name = "", model = c("CFCS"), Cher = c(), NWT = c(), Hemisphere = c(), FF = c(), inst_deposit = c(0), ignore = c(), plotpdf = FALSE, preview = TRUE, plotphoto = FALSE, minphoto = c(), maxphoto = c(), Pbcol = c("black", "midnightblue", "darkgreen"), inst_depositcol = grey(0.85), modelcol = c("black", "red", "darkorange"), historic_d = c(), historic_a = c(), historic_n = c(), historic_test = c(), suppdescriptor = FALSE, descriptor_lab = c(), suppdescriptorcol = c("black", "purple"), coring_yr = c(), plot_Am = FALSE, plot_Cs = FALSE, plot_Pb = TRUE, plot_Pb_inst_deposit = FALSE, plot_CFCS_regression = c(), varves = FALSE, dmin = c(), dmax = c(), sedchange = c(0), min_yr = 1880, SML = c(0), stepout = 5, mycex = 1, archive_metadata = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
Name of the core, given using quotes. Defaults to the core provided with serac. Use preferably the published name of the core for traceability.
}
  \item{coring_yr}{
Year of coring.
}
  \item{model}{
Select 1 to 3 item between c("CFCS", "CIC", "CRS"). If several models are selected, they will all be plotted together in the last window.
}
  \item{Cher}{
If 137Cs measurement were done, where do you detect the Chernobyl peak? The argument is a vector of two depth given in millimeters giving the top and bottom threshold for the 1986 Chernobyl event. The user can run the model without giving any specification before making a decision. In such case, leave the argument empty.
}
  \item{NWT}{
If 137Cs measurement were done, where do you detect the Nuclear Weapon Test peak? The argument is a vector of two depth given in millimeters giving the top and bottom threshold for the 1960s Nuclear Weapon Test event. The user can run the model without giving any specification before making a decision. In such case, leave the argument empty.
}
  \item{Hemisphere}{
Chose between North Hemisphere "NH" and South Hemisphere "SH" depending on the location of your system. This argument is required if you chose to plot NWT.
}
  \item{FF}{
If 137Cs measurement were done, where do you detect the First Fallout period? The argument is a vector of two depth given in millimeters giving the top and bottom threshold for the First Fallout period in 1955. The user can run the model without giving any specification before making a decision. In such case, leave the argument empty.
}
  \item{inst_deposit}{
Upper and lower depths (in mm) of sections of abrupt accumulation that inst_deposit c() should be excised, e.g., c(100, 120, 185, 195) for two sections of 10.0-12.0 cm and 18.5-19.5 cm depth
}
  \item{ignore}{
The depth (in mm) of any sample that should be ignored from the age-depth model computation, e.g., c(55) will remove the measurement done at 5.5 cm. The data will be ploted by default in grey on the output graph (you can change this with the inst_depositcol argument)
}
  \item{plotpdf}{
Logical argument to indicate whether you want the output graph to be saved to your folder.
}
  \item{preview}{
Logical argument to indicate whether you want the output graph to be ploted. Default is TRUE, and the graph is ploted within your R session. It might be convenient to turn this argument to FALSE if errors keep coming telling you your R window is too small.
}
  \item{plotphoto}{
Logical argument to indicate whether you want to plot the photo of the core along your age-model. If plotphoto=TRUE, you need to indicate the upper and lower limit of the photo in mm in following arguments.
}
  \item{minphoto}{
Mandatory if plotphoto=TRUE. Lower limit of the core photo in mm, e.g., minphoto=0 indicates that the photo starts at 0 mm. The photo will automatically be truncated acording to the minimum and maximum depth of the age model given in other arguments.
}
  \item{maxphoto}{
Mandatory if plotphoto=TRUE. Upper limit of the core photo in mm, e.g., maxphoto=320 indicates that the photo ends at 32 cm. The photo will automatically be truncated acording to the minimum and maximum depth of the age model given in other arguments.
}
  \item{Pbcol}{
Vector of color to plot 210Pbex data. If length(Pbcol)>1, the different colors will be used to plot the different slopes in between change(s) in sedimentation rate. Example of color vector: Pbcol=c("black","midnightblue","darkgreen").
}
  \item{inst_depositcol}{
A color to plot the data points within instantaneous deposit or ignored data. Example: inst_depositcol=grey(0.85).
}
  \item{modelcol}{
Vector of color to plot different model if length(model)>1. If length(modelcol)>1, the different colors will be used to plot the different change in sedimentation rate. Example of color vector: modelcol=c("black","red","darkorange") to plot "CFCS", "CIC", "CRS" models in this order.
}
  \item{historic_d}{
Vector with upper and lower depth of historical event(s), e.g., historic_d=c(120,130) will identify the event between 12 and 13 cm on the last window with the age model.
}
  \item{historic_a}{
Vector of year of different historical events, e.g., historic_a=c(1895) will add the point 1895 on the last window with the age model. Historical events can be older than the dated section, in which case the depth is obtained from the model if historic_d is not specified. historic_a is a vector twice as short as historic_d, as each age correspond to an upper+lower limit in the vector 'historic_d'. If not all ages are known, put NA in the vector, e.g., historic_a=c(NA,1895)
}
  \item{historic_n}{
Vector of names of different historical events, e.g., historic_n=c("1895 flood"). Optional. If you plot several events, and don't want to plot all the names, add a NA in the vector, e.g., historic_n=c(NA,"1895 flood") will understand that the first event doesn't have a name, but the second does.
}
  \item{historic_test}{
Visualisation tool for known ages. This argument will plot a vertical line in the last window (the one with the age-depth model). Can be useful when the user know specific ages that may have resulted in changes in sedimentation rates. E.g., historic_test=c(1996).
}
  \item{suppdescriptor}{
Up to two supplementary descriptor(s) to plot in an additional window. Logical argument. The decision on ploting more than one supplementary descriptor depends on the length of the vector descriptor_lab. An additional input file with these data should be included in the folder with the initial data.
}
  \item{descriptor_lab}{
Label used on the axis, e.g., descriptor_lab=c("LOI", "Ca/Fe") if two supplementary descriptors are specified.
}
  \item{suppdescriptorcol}{
Vector of color to plot different descriptor if length(descriptor_lab)>1. If length(descriptor_lab)>1, the different colors will be used to plot the different change in sedimentation rate. Example of color vector: suppdescriptorcol=c("black","purple").
}
  \item{plot_Am}{
Logical argument indicating whether or not serac should plot 241Am.
}
  \item{plot_Cs}{
Logical argument indicating whether or not serac should plot 137Cs.
}
  \item{plot_Pb}{
Logical argument indicating whether or not serac should plot 210Pbex.
}
  \item{plot_Pb_inst_deposit}{
Logical argument indicating whether or not serac should plot 210Pbex without instantaneous deposit. If TRUE, inst_deposit shouldn't be a null vector.
}
  \item{plot_CFCS_regression}{
Whether to plot or not the linear regression. If the parameter is not specified, it will automatically turn to TRUE, but will also automatically turn to FALSE if instantaneous deposit are present but the argument 'plot_Pb_inst_deposit' is turned to FALSE. Linear regression won't match if there are some instantaneous deposit. In other words, in most cases, the user shouldn't need to modify this parameter.
}
  \item{varves}{
Logical argument to indicate whether varve counting results should be ploted on the last window. An additional input file with these data should be included in the folder with the initial data.
}
  \item{dmin}{
Minimum depth of age-depth model (useful if the user doesn't want to plot the upper region).
}
  \item{dmax}{
Maximum depth of age-depth model (useful if the user doesn't want to plot the lower region). dmax cannot be in the middle of an instantaneous deposit. e.g. if there is an instantaneous deposit between 180 and 200 mm, dmax cannot be 190 mm, and will be converted to 200 mm automatically.
}
  \item{dmax}{
Maximum depth of age-depth model (useful if the user doesn't want to plot the lower region).
}
  \item{sedchange}{
Up to two changes in sedimentation rate, e.g., sedchange=c(175,290) indicates two changes of sedimentation rate at 17.5 and 29.0 cm.
}
  \item{min_yr}{
The minimum year limit for the age-depth model plot. The user can adjust this argument after a first computation of the model
}
  \item{SML}{
Surface Mixed Layer: a depth in mm above which the sediment is considered to be mixed. E.g., SML=30 indicates that the first 3 cm are mixed sediment: the data point are ploted but not included in the Pb models.
}
  \item{stepout}{
Depth resolution for the file out in mm.
}
  \item{mycex}{
Graphical parameter: a multiplication factor to increase (mycex>1) ou decrease (mycex<1) label sizes.
}
  \item{archive_metadata}{
Logical argument. If TRUE, require fields regarding the measurements on the core. Allows missing information; just press 'ENTER' in your computer (leave an empty field).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Rosalie Bruel, Pierre Sabatier
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## First time user - setting up the folders
## Line below will create the 'Cores' folder if it doesn't exist already
dir.create(file.path(getwd(), 'Cores'), showWarnings = FALSE)

# Within the 'Cores' folder, user need to create one folder per core
dir.create(file.path(paste(getwd(),'/Cores',sep=""), 'serac_example_ALO09P12'), showWarnings = FALSE)

## We are writing in this subfolder an example dataset
## Format your own data following this template, or get help with formatting using
#    the serac_input_formatting() function
?serac_example_ALO09P12
write.table(x = serac_example_ALO09P12, file = paste(getwd(),'/Cores/serac_example_ALO09P12/serac_example_ALO09P12.txt',sep=""),col.names = T, row.names = F)
## Including proxy data for this core too
write.table(x = serac_example_ALO09P12_proxy, file = paste(getwd(),'/Cores/serac_example_ALO09P12/serac_example_ALO09P12_proxy.txt',sep=""),col.names = T, row.names = F)

## Now run serac
## At the minimum, your function will be:
serac(name="ALO09P12",coring_yr=2009)

## Several sedimentation hypotheses were tested, and this is what the author chose
#    as the best model
serac(name="ALO09P12",coring_yr=2009,model=c("CFCS"),plotphoto=FALSE,minphoto=c(0),maxphoto=c(210),plot_Pb=T,plot_Am=T,plot_Cs=T,Cher=c(30,40),Hemisphere=c("NH"),NWT=c(51,61),sedchange=c(75.5),plot_Pb_inst_deposit=T,inst_deposit=c(20,28,100,107,135,142,158,186),suppdescriptor=TRUE,descriptor_lab=c("Ca/Fe"),historic_d=c(20,28,100,107,135,142,158,186),historic_a=c(1994,1920,1886,1868),historic_n=c("sept1 994 flood","1920 flood","1886 flood","1868 flood ?"), min_yr=c(1750),dmax=c(180), plotpdf=TRUE,preview=F)
## If you want to see the preview, change preview=T in the code
## Before that, make sure to extend your window in RStudio (large plot!)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ age-depth modelling }% use one of  RShowDoc("KEYWORDS")
\keyword{ visualisation }% __ONLY ONE__ keyword per line
