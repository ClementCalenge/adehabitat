\name{area2asc}
\alias{area2asc}
\title{Converts a Polygon to Raster}
\description{
  \code{area2asc} converts a polygon map to a raster map of class
  \code{asc}. It is an alias for \code{mcp.rast}.
}
\usage{
  area2asc(poly, w, border=c("include", "exclude"))
}
\arguments{
  \item{poly}{a data frame with 2 columns giving the coordinates
              of a polygon object}
  \item{w}{an object of class \code{kasc}, or of class \code{asc}}
  \item{border}{a character string indicating what happens when the
    center of the pixel is located exactly on the limit of the polygon
    (\code{"include"} indicates that the pixel is considered to be
    inside the polygon). }
}
\details{
  The raster map is needed to pass the format for the output raster
  object to the function.
}
\value{
  Returns an object of class \code{asc}.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{ \code{\link{hr.rast}}}
\examples{
data(puechabon)
toto <- puechabon$kasc
loc <- puechabon$locs

## gets the coordinates of the relocations for the wild boar #1
wb1 <- loc[loc$Name == "Chou",]
wb1 <- cbind(wb1$X, wb1$Y)
nbpol <- chull(wb1)
xycoord <- wb1[nbpol,]

## rasterization of wb1
tutu <- area2asc(xycoord, toto)
image(tutu)

polygon(xycoord, lwd = 2)
}
\keyword{spatial}
