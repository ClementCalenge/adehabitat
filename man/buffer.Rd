\name{buffer}
\alias{buffer}
\alias{buffer.ani}
\alias{buffer.line}
\alias{buffer.ltraj}
\title{Compute Buffers}
\description{
  \code{buffer} compute buffers around a set of
  locations. \cr
  \code{buffer.ani} is to be used when the
  points can be grouped into several categories  (e.g. the relocations
  of several animals monitored using radio-tracking; the function
  \code{buffer} is then applied to each animal).\cr
  \code{buffer.line} compute buffers around a line.
  \code{buffer.ltraj} compute buffers around one or several animals
  trajectories.
}
\usage{
buffer(pts, x, dist)
buffer.ani(pts, fac, x, dist)
buffer.line(xy, x, dist)
buffer.ltraj(ltraj, x, dist, perani=FALSE)
}

\arguments{
  \item{pts}{a data frame with two columns (x and y coordinates of the
    points)}
  \item{x}{either an object of class \code{asc} or \code{kasc} with the same
    attributes as those desired for the output, or an object of class
    \code{mapattr} (see \code{storemapattr})}
  \item{dist}{a value of distance}
  \item{fac}{a factor defining the categories of the points}
  \item{xy}{a data frame containing the coordinates of the vertices of
    the lines}
  \item{ltraj}{an object of class \code{ltraj}}
  \item{perani}{logical.  If \code{FALSE}, one buffer is computed for
    each burst of relocation.  If \code{TRUE}, the buffers are computed
    per levels of the attribute \code{id}.}
}
\value{
  \code{buffer} and \code{buffer.line} return an object of class
  \code{asc}, with 1 for pixels 
  located within a specified distance of given points, and \code{NA}
  otherwise.\cr 
  \code{buffer.ani} returns a data frame of class \code{kasc},
  with each column corresponding to one level of the factor
  \code{fac}.\cr
  \code{buffer.ltraj} returns a data frame of class \code{kasc},
  with each column corresponding to one level of the attribute \code{id}
  of the object of class \code{ltraj} passed as argument.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link{kasc}} for additionnal information on
  objects of class \code{kasc}, \code{\link{asc}} for further
  information on objects of class \code{asc}, \code{\link{storemapattr}}
  for further information on objects of class \code{mapattr}.
}
\examples{
data(puechabon)

# locs is the data frame containing the
# relocations of wild boars monitored by radio-tracking 
locs <- puechabon$locs

# sa is the "kasc" object of maps of the study area
sa <- puechabon$kasc

# Buffer of 500 m around all relocations
bu <- buffer(locs[,4:5], sa, 500)
image(bu)

# displays all the pixels of the study area within 500 m
# of a relocation of each monitored wild boar
buani <- buffer.ani(locs[,4:5], locs[,1], sa, 500)
image(buani)


## buffer around a trajectory
data(puechcirc)
image(buffer.ltraj(puechcirc,sa, 100))

## buffer around a line
gc <- getcontour(getkasc(sa,1))
out <- buffer.line(gc[,2:3], sa, 300)
image(out)
lines(gc[,2:3], lwd=2)



}
\keyword{spatial}
