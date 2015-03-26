\name{hr.rast}
\alias{hr.rast}
\title{Rasterisation of Objects of Class 'area'}
\description{
Converts an object of class \code{area} (used in many functions of the
package \code{ade4}) to an object of class \code{kasc} (rasterisation).
}
\usage{
hr.rast(mcp, w, border=c("include", "exclude"))
}

\arguments{
  \item{mcp}{an object of class \code{area}}
  \item{w}{a raster map of class \code{kasc} or of class \code{asc}}
  \item{border}{a character string indicating what happens when the
    center of the pixel is located exactly on the limit of the polygon
    (\code{"include"} indicates that the pixel is considered to be
    inside the polygon). }
}
\value{
  Returns an object of class \code{kasc}.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{  \code{\link{kasc}} for additional information on
  objects of class \code{kasc}, \code{\link{area}} for further
  information on the class \code{area} }
\examples{
data(puechabon)
kasc <- puechabon$kasc
locs <-  puechabon$locs

## Computes the home range of the animals
cp <- mcp(locs[,4:5], locs[,1])
area.plot(cp)

## Converts the home range to raster
cprast <- hr.rast(cp, kasc)
image(cprast)
}
\keyword{spatial}
