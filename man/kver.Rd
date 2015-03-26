\name{kver}
\alias{kver}
\alias{kver.rast}
\alias{kver2shapefile}
\alias{plot.kver}
\title{Handling of Objects of Class kver}
\description{
  Objects of class \code{kver} are created by \code{getverticeshr} and
  \code{getverticesclusthr}: they contain the vertices of the home range
  estimated using either nearest neighbour clustering or kernel
  estimation.  An object of class \code{kver} is list of data frames of
  class \code{area}.\cr
  
  \code{plot.kver} is used to display the home-range contours defined in
  an object of clss \code{kver}.\cr

  \code{kver.rast} is used to rasterize the home ranges.\cr

  \code{kver2shapefile} is used to convert the home ranges into a
  shapefile object (for exportation to a GIS).
}
\usage{
kver.rast(kv, asc)
kver2shapefile(kv, which = names(kv))
\method{plot}{kver}(x, which = names(x), colpol = rainbow(length(which)),
          colborder = rep("black", length(which)), lwd = 2,
          add = FALSE, \dots)
}
\arguments{
  \item{kv,x}{an object of class \code{kver}}
  \item{asc}{a matrix of class \code{asc}}
  \item{which}{a vector of character indicating the polygons to be
    plotted}
  \item{colpol}{a vector of the color for filling the polygon.  The default,
    \code{NA}, is to leave polygons unfilled }
  \item{colborder}{a vector of the color to draw the border (black, by
    default).  Use \code{border = NA} to omit borders} 
  \item{lwd}{the border width, a positive number}
  \item{add}{logical. If \code{TRUE}, the polygons are added to a
    previous plot}
  \item{\dots}{additional arguments to be passed to the function
    \code{plot.area}}
}
\value{
  \code{kver.rast} returns an object of class \code{asc}.\cr
  \code{kver2shapefile} returns a shapefile object (can be exported
  using the function \code{write.shapefile} of the package
  \code{shapefiles}.\cr
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\section{Warning}{
  \code{kver2shapefile} requires the packages \code{shapefiles}.  
}
\seealso{\code{\link{getverticesclusthr}} and
  \code{\link{getverticeshr}} for functions creating this class of
  objects. \code{\link{area}} for more information on objects of
  class \code{area}.  \code{\link{asc}} for more information on
  objects of class \code{asc}.
  \code{\link[shapefiles]{write.shapefile}} for exportation of shape
  files.
}
\examples{
data(puechabon)
lo<-puechabon$locs[,c("X","Y")]

## Home Range Estimation
res <- clusthr(lo, puechabon$locs$Name)

## gets the vertices
vec <- getverticesclusthr(res)
plot(vec)
image(kver.rast(vec, getkasc(puechabon$kasc,1)))

}
\keyword{spatial}
