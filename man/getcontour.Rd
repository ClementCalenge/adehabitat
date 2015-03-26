\name{getcontour}
\alias{getcontour}
\title{Computes the Contour Polygon of a Raster Object}
\description{
  \code{getcontour} computes the contour polygon of a raster object
  of class \code{asc}.  When the object is made of several parts, the
  function returns one polygon per part.
}
\usage{
getcontour(x)
}

\arguments{
  \item{x}{an object of class \code{asc}}
}
\value{
  Returns an object of class \code{area}.
}
\section{Warning }{
  Holes in the polygons are not taken into account by the function.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{ \code{\link{asc}} for additionnal information on
  objects of class \code{asc}, \code{\link{area}} for
  additional information on objects of class \code{area}, and
  \code{\link{area2dxf}} to export the results toward a GIS.}

\examples{

data(puechabon)
kasc <- puechabon$kasc

#########################################
## Example with one object:
## Gets the first map of the "kasc" object
## Map of the elevation
elev <- getkasc(kasc, "Elevation")
image(elev)
 
## Get the contour polygon
cpol <- getcontour(elev)

## Draw the polygon
i <- cpol[,2:3]
polygon(i, col = "green", lwd = 2)

#########################################
## Example with two objects:
## home ranges of wild boar

hr <- getsahrlocs(puechabon$sahr, "hr")
u <- getkasc(hr, "Jean")
image(u)

## Get the contour polygons
p <- getcontour(u)
plot(p, lwd = 2)

}
\keyword{programming}
\keyword{hplot}
