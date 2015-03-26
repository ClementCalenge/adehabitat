\name{scatter.enfa}
\alias{scatter.enfa}
\title{Scatter Plot of the Results of the ENFA}
\description{
Performs the scatter diagrams of objects of class \code{enfa}.
}
\usage{
\method{scatter}{enfa}(x, xax = 1, yax = 2, pts = FALSE, nc = TRUE, 
             percent = 95, clabel = 1, side = c("top", "bottom", "none"),
             Adensity, Udensity, Aangle, Uangle, Aborder, Uborder, 
             Acol, Ucol, Alty, Ulty, Abg, Ubg, Ainch, Uinch, \dots) 

}
\arguments{
  \item{x}{an object of class \code{enfa}}
  \item{xax}{the column number for the x-axis}
  \item{yax}{the column number for the y-axis}
  \item{pts}{logical. Whether the points should be drawn.  If
    \code{FALSE}, minimum convex polygons are displayed}
  \item{nc}{whether or not the niche center should be displayed}
  \item{percent}{100 minus the proportion of outliers to be excluded from the
    computation of the minimum convex polygons}
  \item{clabel}{a character size for the columns}
  \item{side}{if \code{"top"}, the legend of the kept axis is upside, if
    \code{"bottom"} it is downside, if \code{"none"} no legend}
  \item{Adensity}{the density of shading lines, in lines per inch, for the 
    available pixels polygon. See \code{\link{polygon}} for more details}
  \item{Udensity}{the density of shading lines, in lines per inch, for the 
    used pixels polygon. See \code{\link{polygon}} for more details}
  \item{Aangle}{the slope of shading lines, given as an angle in degrees
    (counter-clockwise), for the available pixels polygon}
  \item{Uangle}{the slope of shading lines, given as an angle in degrees
    (counter-clockwise), for the used pixels polygon}
  \item{Aborder}{the color for drawing the border of the available pixels 
    polygon. See \code{\link{polygon}} for more details}
  \item{Uborder}{the color for drawing the border of the used pixels polygon. 
    See \code{\link{polygon}} for more details}
  \item{Acol}{the color for filling the available pixels polygon. 
    if \code{pts == FALSE}, the color for the points corresponding to available
    pixels}
  \item{Ucol}{the color for filling the used pixels polygon. 
    if \code{pts == FALSE}, the color for the points corresponding to used 
    pixels}
  \item{Alty}{the line type for the available pixels polygon, as in \code{par}}
  \item{Ulty}{the line type for the used pixels polygon, as in \code{par}}
  \item{Abg}{if \code{pts == TRUE}, background color for open plot symbols of 
    available pixels}
  \item{Ubg}{if \code{pts == TRUE}, background color for open plot symbols of 
    used pixels}
  \item{Ainch}{if \code{pts == TRUE}, heigth in inches of the available pixels} 
  \item{Uinch}{if \code{pts == TRUE}, heigth in inches of the largest used pixels}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
\code{scatter.enfa} displays a factorial map of pixels, as well as the
projection of the vectors of the canonical basis multiplied by a
constant of rescaling.
The kept axes for the plot are specified in a corner.
}
\references{
  Basille, M., Calenge, C., Marboutin, E., Andersen, R. & Gaillard,
  J.M. (2008) Assessing habitat selection using multivariate statistics:
  Some refinements of the ecological-niche factor
  analysis. \emph{Ecological Modelling}, \bold{211}, 233--240.

}
\author{Mathieu Basille \email{basille@ase-research.org}}
\seealso{\code{\link{enfa}}, \code{\link[ade4]{scatter}}}
\examples{
data(lynxjura)

map <- lynxjura$map

## We keep only "wild" indices.
tmp <- lynxjura$locs[,4]!="D"
locs <- lynxjura$locs[tmp, c("X","Y")]

## We perform a square root transformation
## of the variable to normalize it
map[,4] <- sqrt(map[,4])

## We perform the ENFA
tmp <- data2enfa(map, locs[tmp, c("X","Y")])
(enfa1 <- enfa(dudi.pca(tmp$tab, scannf=FALSE),
               tmp$pr, scannf = FALSE))
scatter(enfa1)

}
\keyword{multivariate}
\keyword{hplot}
