\name{join.asc}
\alias{join.asc}
\alias{join.kasc}
\title{Finds the Value of Mapped Variables at some Specified Locations
  (Spatial Join)}
\description{
  \code{join.asc} finds the value of a mapped variable at some specified
  locations.\cr
  \code{join.kasc} is the same function as \code{join.asc}, 
  with several maps.
}
\usage{
 join.asc(pts, x) 
 join.kasc(pts, w) 
}

\arguments{
  \item{x}{an object of class \code{asc}}
  \item{w}{an object of class \code{kasc}}
  \item{pts}{a data frame with two columns containing the x and y
    coordinates of the points to be placed on the map}
}


\value{
  \code{join.asc} returns a vector with length equals to the number of
  points in \code{pts}.\cr
  \code{join.kasc} returns a data frame with a number of columns equals 
  to the number of variables in the object of class \code{kasc}, and
  with each row corresponding to the rows of \code{pts}.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link{kasc}} and \code{\link{asc}}}

\examples{
data(puechabon)
x <- puechabon$kasc

## for each relocation, finds the values of the variables in x
toto <- join.kasc(pts = cbind(puechabon$locs$X, puechabon$locs$Y), x)
toto[1:4,]
}
\keyword{programming}
\keyword{hplot}
