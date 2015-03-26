\name{count.points.id}
\alias{count.points.id}
\alias{count.points}
\title{Number of Points in Each Pixel of a Raster Map}
\description{
  \code{count.points} counts the number of points in
  each pixel of a raster map of class \code{kasc} or \code{asc}.\cr
  \code{count.points.id} counts the number of points in
  each pixel of a raster map of class \code{kasc} or \code{asc}, for
  different sets of points (e.g. the relocations of several animals
  monitored using radio-tracking) 
}
\usage{
count.points(xy, w)
count.points.id(xy, id, w)
}
\arguments{
  \item{xy}{a data frame with 2 columns containing the x and y
    coordinates of the points.}
  \item{id}{a factor giving, for each point, the membership of a point
    to a set.}
  \item{w}{an object of class \code{asc}, \code{kasc} or \code{mapattr}.}
}
\value{
  \code{count.points} returns an object of class \code{asc}
  containing the number of points in each cell of the raster
  map.\cr
  \code{count.points.id} returns an object of class \code{kasc},
  with one column per level of the factor \code{id},
  containing the number of points numbered in each cell of the raster
  map.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link{kasc}} for additionnal information on objects
  of class \code{kasc}, and \code{\link{storemapattr}} for further
  information on objects of class \code{mapattr}.}

\examples{
data(puechabon)
kasc <- puechabon$kasc
locs <-  puechabon$locs

## Counts the number of relocations of each wild boar
## per pixel of the raster map
(nlocrast <- count.points.id(locs[,4:5], locs[,1], kasc))
image(nlocrast)

## Counts the number of all relocations
## per pixel of the raster map
(nlocrast <- count.points(locs[,c("X","Y")], kasc))
image(nlocrast)

}
\keyword{spatial}
