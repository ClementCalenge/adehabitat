\name{asc2im}
\alias{asc2im}
\alias{im2asc}
\title{Conversion of Maps of Class 'asc' and 'im' (Package spatstat)}
\description{
  These functions convert maps of class \code{asc} to objects of
  class \code{im} (package spatstats) and conversely.
}
\usage{
asc2im(x)
im2asc(x)
}
\arguments{
  \item{x}{an object of class \code{asc} or \code{im}}
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link{asc}} for additionnal information on
  objects of class \code{asc}, and \code{\link[spatstat]{im}} for
  additionnal information on objects of class \code{im}}
\examples{
if (require(spatstat)) {
  
#######################
### Conversion asc -> im

  data(puechabon)
  el <- getkasc(puechabon$kasc, "Elevation")
  image(el, main = "An object of class \"asc\"")
  elim <- asc2im(el)
  image(elim, main = "An object of class \"im\"")

#######################
### Conversion im -> asc
  u <- matrix(rnorm(10000), 100, 100)
  haha <- im(u)
  image(haha, main = "class im")
  hihi <- im2asc(haha)
  image(hihi, main = "class asc")

}
}
\keyword{spatial}

