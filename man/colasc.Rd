\name{colasc}
\alias{colasc}
\title{Creates a Vector of Colors for a Raster Map of Type 'factor' }
\description{
  \code{colasc} creates a vector of colors for a raster map of class
  \code{asc} and of type \code{"factor"}. 
}
\usage{
colasc(x, \dots)
}

\arguments{
  \item{x}{an object of class \code{asc}. }
  \item{\dots}{arguments named as the levels of the factor, with
    character values equal 
    to the colors for these levels (see examples)}

}
\value{
  Returns a character vector.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link{asc}} }

\examples{
data(puechabon)

## gets the aspect in an asc object
asp <- getkasc(puechabon$kasc, "Aspect")

## creates the vector of colors
cl <- colasc(asp, NorthEast = "blue", SouthEast = "red", 
             SouthWest = "orange", NorthWest = "green")

## graphical display
image(asp, clfac = cl, main = "Aspect", xlab = "Lambert X",
      ylab = "Lambert Y")
legend(706500, 3162000, legend=levels(asp), fill = cl, cex = 0.7)
 
}
\keyword{programming}
