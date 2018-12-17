\name{serac_map}
\alias{serac_map}
\title{Map sites location.}
\usage{
serac_map(which_lakes,output_name)
}
\description{
Map the location of the sites, if core_metadata have been correctly informed (WGS84 format).
}
\arguments{
\item{which_lakes}{Name of the sites to be plotted. Must match folders names. Default: all the folder in the Cores folder.}
\item{output_name}{Extension to be added to the output name.}
}
\examples{
serac_map()
serac_map(which_lakes = c('PB06', 'Iseo'), output_name = 'subset map')
}
