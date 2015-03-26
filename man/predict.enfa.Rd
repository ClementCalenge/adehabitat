\name{predict.enfa}
\alias{predict.enfa}
\title{Habitat Suitability Maps Built from the ENFA}
\description{
  \code{predict.enfa} computes habitat suitability maps using the
  Ecological-Niche Factor Analysis and the Mahalanobis distances
  method.
}
\usage{
\method{predict}{enfa}(object, index, attr, nf, \dots)
}
\arguments{
  \item{object}{an object of class \code{enfa}}
  \item{index}{an integer vector giving the position of the rows of
    \code{tab} in the  initial object of class \code{kasc}.}
  \item{attr}{an object of class \code{kasc} or \code{mapattr}.} 
  \item{nf}{the number of axes of specialization kept for the
    predictions.  By default, all axes kept in \code{object} are
    used}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
The predictions are based on the position of the niche defined by the
ENFA within the multidimensional space of environmental variables.  The
ENFA produces row coordinates for each pixel, which are used with the
function \code{mahalanobis}.  For each pixel, this function computes the
Mahalanobis distances from the barycentre of the niche.\cr

Actually, the function \code{predict.enfa} is identical to the function
\code{mahasuhab}, except that the habitat suitability map is computed
using the axes of the ENFA, instead of the raw data.\cr

Note that the MADIFA allows a more consistent factorial decomposition of
the Mahalanobis distances.
}
\value{
Returns a raster map of class \code{kasc}.
}
\references{
  Clark, J.D., Dunn, J.E. and Smith, K.G. (1993) A multivariate model of
  female black bear habitat use for a geographic information
  system. \emph{Journal of Wildlife Management}, \bold{57}, 519--526.

  Hirzel, A.H., Hausser, J., Chessel, D. & Perrin, N. (2002)
  Ecological-niche factor analysis: How to compute habitat-suitability
  maps without absence data? \emph{Ecology}, \bold{83}, 2027--2036.  
}
\author{Mathieu Basille \email{basille@ase-research.org}}
\seealso{
\code{\link[stats]{mahalanobis}} for information on the computation of
Mahalanobis distances.  \code{\link{mahasuhab}} for more details on the
computation of habitat suitability maps using the Mahalanobis distances.
\code{\link{madifa}} for a more consistent factorial decomposition of
the Mahalanobis distances
}
\examples{
\dontrun{
data(lynxjura)

map <- lynxjura$map

## We keep only "wild" indices.
tmp <- lynxjura$loc[,4] != "D"
locs <- lynxjura$locs[tmp, c("X","Y")]
dataenfa1 <- data2enfa(map, locs)

(enfa1 <- enfa(dudi.pca(dataenfa1$tab, scannf=FALSE),
               dataenfa1$pr, scannf = FALSE))

## Compute the prediction
pred <- predict(enfa1, dataenfa1$index, dataenfa1$attr)
image(pred)
contour(pred, col="green", add=T)
points(locs, col = "red", pch = 16)
## Lighter areas are the most preferred areas
}
}
\keyword{multivariate}
