\name{squirreloc}
\alias{squirreloc}
\docType{data}
\title{Radio-tracking of squirrels}
\description{
  This data set contains the trajectories of 15 radio-monitored
  squirrels, as well as the vector maps of habitat composition.
}
\usage{data(squirreloc)}
\format{
  This data set is a list of two objects:
  \itemize{
    \item{locs}{is a data frame containing the relocations of 15
      squirrels (with three columns: \code{id, x, y}}
    \item{map}{is an object of class \code{area} containing the habitat
      composition of the area. The habitat types and colour coding are
      stored in the attribute "info" of this object}
  }
}
\details{
  The dataset \code{squirreloc} comes from the Ranges VI software.  It
  has been used to illustrate the compositional analysis (see
  \code{?compana}) and the eigenanalysis of selection ratios (see
  \code{?eisera}).  See also the dataset \code{squirrel}.
}
\source{
  Kenward, R.E., South, A.B. and Walls, S.S. (2003). Ranges6 v1.2 : For
  the analysis of tracking and location data. Online manual. Anatrack
  Ltd. Wareham, UK. ISBN 0-9546327-0-2.
}
\examples{
data(squirreloc)


## habitat:
are <- squirreloc$map
co <- attr(are, "info")
plot(are, colp = co[,2])

legend(-10, 210, unique(co[,1]),
       unique(co[,2]), bg="white")



## relocations
locs <- squirreloc$locs
li <- split(locs[,2:3], locs[,1])
opar <- par(mfrow=n2mfrow(length(li)), mar=c(0,0,2,0))
lapply(1:length(li), function(i) {
plot(are, colp = co[,2], main=names(li)[i], axes=FALSE)
points(li[[i]], pch=16, cex=1.5)
box()
})
par(opar)

}
\keyword{datasets}
