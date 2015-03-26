\name{morphology}
\alias{morphology}
\title{Morphology: Erosion or Dilatation of Features on a Raster Map}
\description{
  \code{morphology} performs morphological operations on images
  of class \code{asc}.
}
\usage{
morphology(x, operation = c("erode", "dilate"), nt = 5)
}
\arguments{
  \item{x}{a matrix of class \code{asc}}
  \item{operation}{a character string indicating the operation to be
    processed: either \code{"erode"} or \code{"dilate"}}
  \item{nt}{the number of times that the operation should be processed}
}
\value{
Returns a matrix of class \code{asc}, containing 1 when the pixel belong
to one feature of the image and \code{NA} otherwise (see examples).  
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{\code{\link{asc}} for further information on objects of
  class \code{asc}.} 
\examples{
data(puechabon)
a <- getkasc(puechabon$kasc,"Elevation")

## dilatation
toto1 <- morphology(a, operation = "dilate", nt = 1)
toto2 <- morphology(a, operation = "dilate", nt = 2)
toto3 <- morphology(a, operation = "dilate", nt = 3)
toto5 <- morphology(a, operation = "dilate", nt = 5)
colo <- grey((1:5)/6)
image(toto5, col = colo[1])
image(toto3, add = TRUE, col = colo[2])
image(toto2, add = TRUE, col = colo[3])
image(toto1, add = TRUE, col = colo[4])
image(a, add = TRUE)

## erosion
colo <- grey((1:20)/21)
image(a, col = 1)
for (i in 1:19) {
  toto <- morphology(a, operation = "erode", nt = i)
  image(toto, add = TRUE, col = colo[i])
}


}
\keyword{spatial}

