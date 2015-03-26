\name{distfacmap}
\alias{distfacmap}
\title{Compute distances to the different levels of a factor map }
\description{
  This function computes map of distances to patches belonging to the
  different levels of a map of class \code{asc} and of type
  \code{factor}.
}
\usage{
distfacmap(x)
}
\arguments{
  \item{x}{an object of class \code{asc} and of type \code{factor} }
}
\value{
  An object of class \code{kasc}.
}
\author{  Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{ \code{\link{asc}} for further info on objects of class
  \code{asc}. To compute distance maps from points, lines or polygons,
  see \code{\link[spatstat]{distmap}} in the package \code{spatstat}.
}
\examples{
\dontrun{
data(puechabon)
asp <- getkasc(puechabon$kasc, "Aspect")
image(asp)
sor <- distfacmap(asp)
image(sor)
}
}
\keyword{spatial}
