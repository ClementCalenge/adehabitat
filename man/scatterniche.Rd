\name{scatterniche}
\alias{scatterniche}
\title{Display the Niche in the Ecological Space }
\description{
  \code{scatterniche} displays the niche in the Ecological space
  (multidimensional space defined by habitat variables).
}
\usage{
scatterniche(x, pr, xax = 1, yax = 2, pts = FALSE,
             percent = 95, clabel = 1,
             side = c("top", "bottom", "none"),
             Adensity, Udensity, Aangle, Uangle, Aborder,
             Uborder, Acol, Ucol, Alty, Ulty, Abg,
             Ubg, Ainch, Uinch, \dots)
}
\arguments{
  \item{x}{a data frame giving the value of environmental variables
    (columns) in resource units (rows, e.g. the pixels of a raster map)}
  \item{pr}{a vector giving the utilisation weight for each resource
    unit}
  \item{xax}{the column number for the x-axis}
  \item{yax}{the column number for the y-axis}
  \item{pts}{logical. Whether the points should be drawn.  If
    \code{FALSE}, minimum convex polygons are displayed}
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
  \item{Aborder}{the color to draw the border of the available pixels 
    polygon. See \code{\link{polygon}} for more details}
  \item{Uborder}{the color to draw the border of the used pixels polygon. 
    See \code{\link{polygon}} for more details}
  \item{Acol}{the color for filling the available pixels polygon. 
    if \code{pts==FALSE}, the color for the points corresponding to
    available pixels}
  \item{Ucol}{the color for filling the used pixels polygon. 
    if \code{pts==FALSE}, the color for the points corresponding to used 
    pixels}
  \item{Alty}{the line type for the available pixels polygon, as in
    \code{par}.}
  \item{Ulty}{the line type for the used pixels polygon, as in \code{par}.}
  \item{Abg}{if \code{pts==TRUE}, background color for open plot symbols of 
    available pixels}
  \item{Ubg}{if \code{pts==TRUE}, background color for open plot symbols of 
    used pixels}
  \item{Ainch}{if \code{pts==TRUE}, heigth in inches of the available pixels} 
  \item{Uinch}{if \code{pts==TRUE}, heigth in inches of the largest used
    pixels}
  \item{\dots}{further arguments passed to or from other methods}
}
\author{Mathieu Basille \email{basille@ase-research.org}\cr
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\examples{
data(chamois)
cp <- count.points(chamois$locs, chamois$map)
chamois$map
li <- kasc2df(chamois$map)
cpi <- c(cp)[li$index]
## we focus on the distance to ecotone and on the slope,
## after centring and scaling (with the help of a PCA)
scatterniche(dudi.pca(li$tab[,2:3], scannf=FALSE)$tab, cpi)
scatterniche(dudi.pca(li$tab[,2:3], scannf=FALSE)$tab, cpi, pts=TRUE)

}
\keyword{hplot}
