\name{serac_map}
\alias{serac_map}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate map of sites in Cores folder
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
serac_map(which_lakes = c(), output_name = c())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{which_lakes}{
Name of the sites to be plotted. Must match folders names. Default: all the folder in the Cores folder.
}
  \item{output_name}{
Extension to be added to the output name (can be left blank).
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
Rosalie Bruel
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
function (which_lakes = c(), output_name = c())
.serac_map(which_lakes, output_name)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
