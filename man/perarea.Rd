\name{perarea}
\alias{perarea}
\alias{ararea}
\title{Compute Areas and Perimeters of Objects of Class "area" }
\description{
  \code{perarea} computes the perimeters of polygons in objects of class
  \code{area}.\cr
  \code{ararea} computes the areas of polygons in objects of class
  \code{area}.\cr
}
\usage{
perarea(x)
ararea(x)
}
\arguments{
  \item{x}{object of class \code{area} }
}
\value{
  a vector.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{area}}, \code{\link{plot.area}}}
\examples{

data(puechabon)
locs <-  puechabon$locs

cp <- mcp(locs[,4:5], locs[,1])

perarea(cp)

ararea(cp)

}
\keyword{hplot}

