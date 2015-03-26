\name{ascgen}
\alias{ascgen}
\title{Creation of Raster Maps}
\description{
  \code{ascgen} creates an object of class \code{asc} using
  a set of points.
}
\usage{
ascgen(xy = NULL, cellsize = NULL, nrcol = 10, count = TRUE)
}

\arguments{
  \item{xy}{a data frame with two columns containing the x and y
    coordinates of the points}
  \item{cellsize}{the \code{cellsize} attribute
    of the object of class \code{asc} to be created}
  \item{nrcol}{the size of the square raster map to
    be created (number of rows and columns)}
  \item{count}{logical.  If \code{TRUE}, the object of class \code{asc}
    contains the number of points in each cell.  If \code{FALSE}, all
    the cells are set to zero} 
}
\value{
  Returns an object of class \code{asc}.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{\code{\link{asc}} for additional information on objects
  of class \code{asc}.}

\examples{
data(puechabon)
lo <- puechabon$locs[,c("X","Y")]
plot(lo, asp = 1, pch = 16)

## lo contains the relocations of wild boars
rast <- ascgen(lo, cellsize = 100)
image(rast)

## Alternatively, one can specify the size of the square raster map
rast <- ascgen(lo, nrcol = 10)
rast
image(rast)

## can be used for further analyses
## (e.g. correspondence analyses)
locs <- puechabon$locs[, c("Name", "X", "Y")]
o <- count.points.id(locs[,2:3], locs[,1], rast)
image(o)
}
\keyword{spatial}
