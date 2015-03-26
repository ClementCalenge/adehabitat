\name{domain}
\alias{domain}
\title{Estimation of the Potential Distribution of a Species }
\description{
\code{domain} uses the DOMAIN algorithm to estimate the potential
distribution of a species based on a list of species occurrences and on
maps of the area.
}
\usage{
domain(kasc, pts, type = c("value", "potential"), thresh = 0.95)
}
\arguments{
  \item{kasc}{an object of class \code{kasc}}
  \item{pts}{a data frame giving the x and y coordinates of the species
    occurrences.}
  \item{type}{a character string. The \code{"value"} of the suitability
    may be returned 
    or the \code{"potential"} area of distribution}
  \item{thresh}{if \code{value = "potential"}, a threshold value should be
    supplied for the suitability (by default 0.95)}
}
\details{
  This function implements the DOMAIN algorithm described in Carpenter
  et al. (1993).
}
\value{
  Returns a matrix of class \code{asc}.
}
\section{Warning }{
  \code{domain} is restricted to maps
  containing only numerical variables (i.e. no factors).
}
\references{
  Carpenter, G., Gillison, A.N. and Winter, J. (1993)
  DOMAIN: a flexible modelling
  procedure for mapping potential distributions of plants and
  animals. \emph{Biodiversity and conservation}, \bold{2}, 667--680.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{ \code{\link{kasc}} for additionnal information on objects
  of class \code{kasc}, \code{\link{asc}} for additionnal information on
  matrices of class \code{asc}.  }

\examples{

## Preparation of the data
data(puechabon)
kasc <- puechabon$kasc
kasc$Aspect <- NULL
pts <- puechabon$locs[puechabon$locs$Name == "Brock", 4:5]

## View of the data
elevation <- getkasc(kasc, "Elevation")
image(elevation)
points(pts, col = "red", pch = 16)

## Estimation of habitat suitability map
hsm <- domain(kasc, pts)
image(hsm, col = grey((1:256)/256))
contour(hsm, add = TRUE)
## Lighter areas are the most preferred areas

## Potential distribution
hsm <- domain(kasc, pts, type = "potential")
image(elevation, main = "Habitat suitability map")
image(hsm, add = TRUE, col = "orange")
points(pts, col = "red", pch = 16)

}
\keyword{spatial}
\keyword{multivariate}
