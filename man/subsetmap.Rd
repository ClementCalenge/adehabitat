\name{subsetmap}
\alias{subsetmap}
\alias{subsetmap.asc}
\alias{subsetmap.kasc}
\title{Storing a Part of a Map}
\description{
  \code{subsetmap} is a generic function.  It has methods for the classes
  \code{asc} and \code{kasc}.  It is used to store a part of any given
  map into an other object.
}
\usage{
subsetmap(x, xlim = NULL, ylim = NULL, ...)
}
\arguments{
  \item{x}{an object of class \code{asc} or \code{kasc}}
  \item{xlim}{numerical vector of length 2.  The x limits of the rectangle
    including the new map}
  \item{ylim}{numerical vector of length 2.  The y limits of the rectangle
    including the new map}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  If \code{xlim} or \code{ylim} are not provided, the function asks the
  user to click on the map to delimit the lower left corner and the
  higher right corner of the new map (see Examples).
}
\value{
  Returns an object of class \code{asc} or \code{kasc}
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr},
  improvements by Jon Olav Vik}
\seealso{\code{\link{asc}}, \code{\link{kasc}}}
\examples{

data(puechabon)
kasc <- puechabon$kasc
el <- getkasc(kasc, "Elevation")

## limits of the new map:
xl <- c(701561, 704017)
yl <- c(3160560, 3162343)

## computation of the new map:
su <- subsetmap(el, xlim = xl, ylim = yl)

## Display
opar <- par(mar = c(0,0,0,0))
layout(matrix(c(1,1,1,1,1,1,1,1,2), byrow = TRUE, ncol = 3))
image(el, axes = FALSE)
polygon(c(xl[1], xl[2], xl[2], xl[1]),
        c(yl[1], yl[1], yl[2], yl[2]))
image(su, axes = FALSE)
box()

par(opar)
par(mfrow = c(1,1))

### Gets this part for the whole kasc object
m <- subsetmap(kasc, xlim = xl, ylim = yl)
image(m)


\dontrun{
 ## Interactive example 
 su <- subsetmap(kasc)

 image(su)
}

}
\keyword{spatial}
