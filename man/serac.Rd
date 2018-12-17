\name{serac}
\alias{serac}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
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
%%     ~~Describe \code{name} here~~
}
  \item{model}{
%%     ~~Describe \code{model} here~~
}
  \item{Cher}{
%%     ~~Describe \code{Cher} here~~
}
  \item{NWT}{
%%     ~~Describe \code{NWT} here~~
}
  \item{Hemisphere}{
%%     ~~Describe \code{Hemisphere} here~~
}
  \item{FF}{
%%     ~~Describe \code{FF} here~~
}
  \item{inst_deposit}{
%%     ~~Describe \code{inst_deposit} here~~
}
  \item{ignore}{
%%     ~~Describe \code{ignore} here~~
}
  \item{plotpdf}{
%%     ~~Describe \code{plotpdf} here~~
}
  \item{preview}{
%%     ~~Describe \code{preview} here~~
}
  \item{plotphoto}{
%%     ~~Describe \code{plotphoto} here~~
}
  \item{minphoto}{
%%     ~~Describe \code{minphoto} here~~
}
  \item{maxphoto}{
%%     ~~Describe \code{maxphoto} here~~
}
  \item{Pbcol}{
%%     ~~Describe \code{Pbcol} here~~
}
  \item{inst_depositcol}{
%%     ~~Describe \code{inst_depositcol} here~~
}
  \item{modelcol}{
%%     ~~Describe \code{modelcol} here~~
}
  \item{historic_d}{
%%     ~~Describe \code{historic_d} here~~
}
  \item{historic_a}{
%%     ~~Describe \code{historic_a} here~~
}
  \item{historic_n}{
%%     ~~Describe \code{historic_n} here~~
}
  \item{historic_test}{
%%     ~~Describe \code{historic_test} here~~
}
  \item{suppdescriptor}{
%%     ~~Describe \code{suppdescriptor} here~~
}
  \item{descriptor_lab}{
%%     ~~Describe \code{descriptor_lab} here~~
}
  \item{suppdescriptorcol}{
%%     ~~Describe \code{suppdescriptorcol} here~~
}
  \item{coring_yr}{
%%     ~~Describe \code{coring_yr} here~~
}
  \item{plot_Am}{
%%     ~~Describe \code{plot_Am} here~~
}
  \item{plot_Cs}{
%%     ~~Describe \code{plot_Cs} here~~
}
  \item{plot_Pb}{
%%     ~~Describe \code{plot_Pb} here~~
}
  \item{plot_Pb_inst_deposit}{
%%     ~~Describe \code{plot_Pb_inst_deposit} here~~
}
  \item{plot_CFCS_regression}{
%%     ~~Describe \code{plot_CFCS_regression} here~~
}
  \item{varves}{
%%     ~~Describe \code{varves} here~~
}
  \item{dmin}{
%%     ~~Describe \code{dmin} here~~
}
  \item{dmax}{
%%     ~~Describe \code{dmax} here~~
}
  \item{sedchange}{
%%     ~~Describe \code{sedchange} here~~
}
  \item{min_yr}{
%%     ~~Describe \code{min_yr} here~~
}
  \item{SML}{
%%     ~~Describe \code{SML} here~~
}
  \item{stepout}{
%%     ~~Describe \code{stepout} here~~
}
  \item{mycex}{
%%     ~~Describe \code{mycex} here~~
}
  \item{archive_metadata}{
%%     ~~Describe \code{archive_metadata} here~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (name = "", model = c("CFCS"), Cher = c(), NWT = c(), 
    Hemisphere = c(), FF = c(), inst_deposit = c(0), ignore = c(), 
    plotpdf = FALSE, preview = TRUE, plotphoto = FALSE, minphoto = c(), 
    maxphoto = c(), Pbcol = c("black", "midnightblue", "darkgreen"), 
    inst_depositcol = grey(0.85), modelcol = c("black", "red", 
        "darkorange"), historic_d = c(), historic_a = c(), historic_n = c(), 
    historic_test = c(), suppdescriptor = FALSE, descriptor_lab = c(), 
    suppdescriptorcol = c("black", "purple"), coring_yr = c(), 
    plot_Am = FALSE, plot_Cs = FALSE, plot_Pb = TRUE, plot_Pb_inst_deposit = FALSE, 
    plot_CFCS_regression = c(), varves = FALSE, dmin = c(), dmax = c(), 
    sedchange = c(0), min_yr = 1880, SML = c(0), stepout = 5, 
    mycex = 1, archive_metadata = FALSE) 
.serac(name, model, Cher, NWT, Hemisphere, FF, inst_deposit, 
    ignore, plotpdf, preview, plotphoto, minphoto, maxphoto, 
    Pbcol, inst_depositcol, modelcol, historic_d, historic_a, 
    historic_n, historic_test, suppdescriptor, descriptor_lab, 
    suppdescriptorcol, coring_yr, plot_Am, plot_Cs, plot_Pb, 
    plot_Pb_inst_deposit, plot_CFCS_regression, varves, dmin, 
    dmax, sedchange, min_yr, SML, stepout, mycex, archive_metadata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
