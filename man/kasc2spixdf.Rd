\name{kasc2spixdf}
\alias{kasc2spixdf}
\alias{asc2spixdf}
\alias{spixdf2kasc}
\alias{area2spol}
\alias{spol2area}
\alias{attpol2area}
\alias{traj2spdf}
\alias{traj2sldf}
\alias{ltraj2spdf}
\alias{ltraj2sldf}
\alias{kver2spol}
\title{ Conversion of maps from/to the package "sp" }
\description{
  These functions convert maps of classes available in adehabitat toward
  classes available in the package \code{sp} and conversely.\cr
  \cr
  \code{kasc2spixdf} converts an object of class \code{kasc} into an
  object of class \code{SpatialPixelsDataFrame}.\cr
  \cr
  \code{asc2spixdf} converts an object of class \code{asc} into an
  object of class \code{SpatialGridDataFrame}.\cr
  \cr
  \code{spixdf2kasc} converts an object of class
  \code{SpatialPixelsDataFrame} or \code{SpatialGridDataFrame} into an
  object of class \code{asc} or \code{kasc}.\cr 
  \cr
  \code{area2spol} converts an object of class \code{area} into an
  object of class \code{SpatialPolygons}.\cr
  \cr
  \code{spol2area} converts an object of class \code{SpatialPolygons} or
  \code{SpatialPolygonsDataFrame} into an object of class
  \code{area}.\cr
  \cr
  \code{attpol2area} gets the data attribute of an object of class
  \code{SpatialPolygonsDataFrame} and stores it into a data frame.\cr
  \cr
  \code{traj2spdf} converts an object of class \code{traj} into an
  object of class \code{SpatialPointsDataFrame}.\cr
  \cr
  \code{traj2sldf} converts an object of class \code{traj} into an
  object of class \code{SpatialLinesDataFrame}.\cr
  \cr
  \code{ltraj2spdf} converts an object of class \code{ltraj} into an
  object of class \code{SpatialPointsDataFrame}.\cr
  \cr
  \code{ltraj2sldf} converts an object of class \code{ltraj} into an
  object of class \code{SpatialLinesDataFrame}.\cr
  \cr
  \code{kver2spol} converts an object of class \code{kver} into an
  object of class \code{SpatialPolygons}.
}
\usage{
kasc2spixdf(ka)
asc2spixdf(a)
spixdf2kasc(sg)
area2spol(ar)
spol2area(sr)
attpol2area(srdf)
traj2spdf(tr)
traj2sldf(tr, byid = FALSE)
ltraj2spdf(ltr)
ltraj2sldf(ltr, byid = FALSE)
kver2spol(kv)
}
\arguments{
  \item{ka}{an object of class \code{kasc}.}
  \item{a}{an object of class \code{asc}.}
  \item{sg}{an object of class \code{SpatialPixelsDataFrame} or
    \code{SpatialGridDataFrame}.}
  \item{ar}{an object of class \code{area}.}
  \item{sr}{an object of class \code{SpatialPolygons} or
    \code{SpatialPolygonsDataFrame}.}
  \item{srdf}{an object of class \code{SpatialPolygonsDataFrame}.}
  \item{tr}{an object of class \code{traj}.}
  \item{ltr}{an object of class \code{ltraj}.}
  \item{kv}{an object of class \code{kver}.}
  \item{byid}{logical.  If \code{TRUE}, one objects of class
    \code{Lines} correspond to one animal.  if \code{FALSE},
    one object of class \code{Lines} correspond to one burst.}
}
\details{
  We describe here more in detail the functions \code{spol2area} and
  \code{attpol2area}.  Objects of class \code{area} do not deal with
  holes in the polygons, whereas the objects of class
  \code{SpatialPolygons} do.  Therefore, when holes are present in the
  \code{SpatialPolygons} object passed as argument, the function
  \code{spol2area} ignore them and returns only the external contour of
  the polygon (though a warning is returned).
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr}, code
  of \code{spixdf2kasc} kindly provided by Roger Bivand
  \email{Roger.Bivand@nhh.no} }
\seealso{ \code{\link{asc}} for information on objects of
  class \code{asc}, \code{\link{kasc}} for info on objects of
  class \code{kasc}, \code{\link{area}} for info on objects of
  class \code{area}, \code{\link{traj}} for objects of class
  \code{traj}, \code{\link{ltraj}} for objects of class
  \code{ltraj}.}
\examples{
\dontrun{
if (require(sp)) {

#########################################
##
## Conversion kasc -> SpatialPixelsDataFrame
##

data(puechabon)
toto <- kasc2spixdf(puechabon$kasc)
image(toto)
summary(toto)

#### and conversely
toto <- spixdf2kasc(toto)
image(toto)
hist(toto)

data(meuse.grid)
m <- SpatialPixelsDataFrame(points = meuse.grid[c("x", "y")],
                           data = meuse.grid)
i <- spixdf2kasc(m)
image(i)


### conversion asc -> SpatialPixelsDataFrame
cuicui <- asc2spixdf(getkasc(toto,1))
image(cuicui)


#########################################
##
## Conversion kver -> SpatialPolygons
##

hr <- getverticeshr(kernelUD(puechabon$locs[,c("X","Y")], puechabon$locs$Name,
                    grid=100))

plot(hr)
class(hr)
spo <- kver2spol(hr)
plot(spo)
class(spo)

#########################################
##
## Conversion area -> SpatialPolygons
##

data(elec88)
ar <- as.area(elec88$area)
plot(ar)
toto <- area2spol(ar)
plot(toto)


#########################################
##
## Conversion SpatialPolygons -> area
##

## First create an object of class "SpatialRingsDataFrame"
tutu <- SpatialPolygonsDataFrame(toto, elec88$tab)

## and then conversion:
coincoin <- spol2area(tutu)
plot(coincoin)

## gets the attributes
haha <- attpol2area(tutu)
area.plot(coincoin, values = haha$Waechter)


#########################################
##
## Conversion ltraj -> SpatialPointsDataFrame
##

data(puechcirc)
plot(puechcirc)

toto <- ltraj2spdf(puechcirc)
plot(toto)



#########################################
##
## Conversion ltraj -> SpatialLinesDataFrame
##

toto <- ltraj2sldf(puechcirc)
plot(toto)


}
}
}
\keyword{hplot}
