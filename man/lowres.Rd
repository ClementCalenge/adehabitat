\name{lowres}
\alias{lowres}
\alias{lowres.asc}
\alias{lowres.kasc}
\title{Reducing the Resolution of a Map}
\description{
  \code{lowres} is a generic function, having methods for the classes
  \code{asc} and \code{kasc}. It is used to reduce the resolution of
  the maps.
}
\usage{
lowres(x, np = 2, ...)
}
\arguments{
  \item{x}{an object of class \code{asc} or \code{kasc}}
  \item{np}{a number giving the number of pixels to merge together
    (see below)}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  The function merges together squares of \code{np * np} pixels.  For
  maps of type \code{"numeric"} (see
  \code{help(asc)}), the function averages the value of the variable.
  For maps of type \code{"factor"}, the function gives the most
  frequent level in the square of \code{np * np} pixels.  When several
  levels are equally represented in the square of \code{np * np} pixels,
  the function randomly samples one of these levels.
}
\value{
  Returns an object of class \code{asc} or \code{kasc}.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr}}
\seealso{\code{\link{asc}} and \code{\link{kasc}} for further
  information on objects of class \code{asc} and \code{kasc}
  respectively.}
\examples{
data(puechabon)
kasc <- puechabon$kasc

## The initial image
image(kasc)

## The transformed image
m <- lowres(kasc, np = 4)
image(m)

}
\keyword{spatial}
