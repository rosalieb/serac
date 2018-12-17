\name{help_serac}
\alias{help_serac}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
help_serac()
}
%- maybe also 'usage' for other objects documented here.
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
function () 
{
    cat("Hi there, welcome to serac.\n1) GET STARTED\n   In the directory (folder) your working in, you should have created a folder named 'Cores'.\n   Within this folder, you will be able to store all your different projects.\n   Pick a name for your project - we recommend using the name of the sediment sequence.\n   Within this folder, you need to have at least one .txt file, named the exact same way\n   than your sub-folder.\n   Here are some requirements for this .txt file:\n          - saved using tabulations as delimitators;\n          - points for decimal delimitators;\n          - all depths must be in mm;\n          - there are two format possible to enter depths:\n                (1) one column depth + one column thickness\n                (2) upper and lower limits of your sample, \n                    ID with 'depth_top' OR 'depth_min' and 'depth_bottom' OR 'depth_max';\n          - 'density' must be in g/cm3;\n          - your different radioelement measurements, named 'Pb210', 'Pb210ex', 'Cs137',\n             'Am241';\n          - the error for each measurements, named 'Pb210_er', 'Pb210ex_er','Cs137_er', \n            'Am241_er'.\n2) HAVE A FIRST LOOK AT YOUR DATA\n   You can have a first look at your data just by typing into your console:\n           serac(name='mydata',coring_yr=2010)\n   Of course, 'mydata' must match the folder name you just created, and the coring\n   year might not be 2010.\n3) TEST HYPOTHESES\n   Using the different argument, you can test different hypotheses, and chose what\n   to plot on your output figure.\n4) OTHER FILES FORMAT\n   Core picture extension should be .jpg (not .jpeg, not .png, not .pdf).\n   If you're adding a proxy or varves counting, the depth column must contain the\n   average depth of the sample (in mm) and then your data in other columns.\n   Again, you must save as .txt (tabulations separators).\n5) COMMON ERRORS\n   If you get errors, first make sure your data are in the correct format\n   If the message is related to the size of the window being too small, try extending\n   your plot window (you can do that at least in R studio).\n   Try setting preview to FALSE too.\n   It is a silly problem we could not fix - solutions depends on the platform you use.\n   We decided we didn't want serac to be utilisable on only one platform.\n6) PACKAGES\n   serac uses packages, they will be downloaded automatically.\n   Make sure to be connected to internet the first time you run the code so the\n   packages can be added to your library.\n7) CONTACT US\n   If the problem persists and it feels like you tried everything, we're happy to help.\n   Send us your file in the correct format (as described above) + the line of code you\n   used + the error message you got.")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
