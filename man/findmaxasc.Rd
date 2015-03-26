\name{findmaxasc}
\alias{findmaxasc}
\title{Find Local Maxima on a Map of Class 'asc'}
\description{
  \code{findmaxasc} finds the local maxima on a map of class
  \code{asc}.
}
\usage{
findmaxasc(asc)
}
\arguments{
  \item{asc}{a map of class \code{asc}}
}
\details{
  This function may be useful, among other things, to identify the local
  modes of the utilization distribution of an animal (e.g. estimated
  using \code{kernelUD}). 
}
\value{
  a data frame with two columns containing the x and y coordinates of
  the local maxima.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{ \code{\link{as.asc}}, \code{\link{kernelUD}} }
\examples{

## example, to find the mode of the utilization distribution (UD) of an
## animal:

## load the data
data(puechabon)
loc <- puechabon$locs[, c("X", "Y")]
id <- puechabon$locs[, "Name"]

## Estimation of UD for the four animals
(ud <- kernelUD(loc, id))
image(ud)


## Now consider the UD of Chou
map <- ud$Chou$UD
image(map)

## Find the local maxima:
maxim <- findmaxasc(map)
points(maxim, col="red", pch=16, cex=1.5)



}
\keyword{spatial}
