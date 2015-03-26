\name{as.area}
\alias{as.area}
\alias{area}
\title{Objects of Class "area"}
\description{
  Objects of class \code{area} are used to store the information on
  the shape of vectorised objects.\cr
  An object of class \code{area} is a data frame with three variables.
  The first variable is a factor defining the polygons.\cr
  The second and third variables are the xy coordinates of the
  polygon vertices in the order where they are found.  This kind
  of objects are current in the package \code{ade4}, though this package
  does not define \code{area} as a special class.
}
\usage{
as.area(x)
}
\arguments{
  \item{x}{a data frame with three columns}
}
\value{
  Returns an object of class \code{area}
}
\seealso{ \code{\link[ade4]{area.plot}} (package ade4) for other
  operations on objects of class \code{area}, \code{\link{plot.area}} to
  display objects of this class, and \code{\link{area2dxf}} for
  exportation of the objects \code{area} toward a GIS. }
\examples{
data(elec88)
ar <- as.area(elec88$area)
plot(ar)
}
\keyword{spatial}
