\name{mround}
\alias{mround}
\title{Round values}
\usage{
mround(x, base)
}
\description{
Question the user on professional ID (name, ORCID, email, Institute). These metadata are then included in the general metadata file output every time a model is created.
}
\arguments{
\item{x}{The value to round.}
\item{base}{What base value should be used for rounding.}
}
\examples{
mround(142, 10)
mround(-43, 5)
}
