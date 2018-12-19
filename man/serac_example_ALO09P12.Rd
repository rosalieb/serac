\name{serac_example_ALO09P12}
\alias{serac_example_ALO09P12}
\docType{data}
\title{
Lake Allos shortlived radionuclide and proxy data
}
\description{
Example of data input needed for serac.
}
\usage{serac_example_ALO09P12}
\format{
  A data frame with 24 observations on the following 9 variables.
  \describe{
    \item{\code{depth_min}}{a numeric vector of the upper depth of the sample}
    \item{\code{depth_max}}{a numeric vector of the lower depth of the sample}
    \item{\code{density}}{a numeric vector of the density calculated for the sample}
    \item{\code{X210Pbex}}{a numeric vector of excess 210Pb}
    \item{\code{X210Pbex_er}}{a numeric vector of error excess 210Pb}
    \item{\code{X137Cs}}{a numeric vector}
    \item{\code{X137Cs_er}}{a numeric vector}
    \item{\code{X241Am}}{a numeric vector}
    \item{\code{X241Am_er}}{a numeric vector}
  }
}
\details{
Lake Allos (2230 m a.s.l., 44° 14' N, 6° 42' 35" E) is a 1-km-long by 700-m-wide, high-altitude lake in the French Alps. It is located ap- proximately 100 km from the Mediterranean Sea and situated in a 5-km2-north-facing cirque, surrounded by sandstone mountains that rise to 2672 m a.s.l.
}
\source{
\href{https://www.sciencedirect.com/science/article/pii/S0033589412000294}{Science Direct}
}
\references{
Wilhelm, B., et al. (2012) Quat. Res. 78, 1–12.
}
\examples{
serac_example_ALO09P12
str(serac_example_ALO09P12)
}
\keyword{datasets}
