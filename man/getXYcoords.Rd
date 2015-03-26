\name{getXYcoords}
\alias{getXYcoords}
\title{Computes the X and Y Coordinates of the Pixels of a Raster Map}
\description{
  \code{getXYcoords} computes the geographical coordinates of the rows
  and columns of pixels of a raster map.
}
\usage{
getXYcoords(w)
}
\arguments{
  \item{w}{an object of class \code{asc}, \code{kasc}, \code{sahrlocs},
    or \code{mapattr}.}
}
\value{
  Returns a list with two components:
  \item{x}{the x coordinates of the columns of pixels of the map}
  \item{y}{the y coordinates of the rows of pixels of the map}
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{ \code{\link{asc}}, \code{\link{kasc}},
  \code{\link{as.sahrlocs}}, \code{\link{storemapattr}} for additionnal
  informations on objects of class \code{area}, \code{kasc},
  \code{sahrlocs}, \code{mapattr} respectively}
\examples{
data(puechabon)
(elev <- getkasc(puechabon$kasc, "Elevation"))
(coords <- getXYcoords(elev))
nrow(elev) == length(coords$x)
ncol(elev) == length(coords$y)
}
\keyword{programming}


