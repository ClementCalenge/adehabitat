\name{madifa}
\alias{madifa}
\alias{print.madifa}
\alias{hist.madifa}
\alias{scatter.madifa}
\alias{predict.madifa}
\alias{plot.madifa}
\alias{s.madifa}
\title{The MADIFA: a Factorial Decomposition of the Mahalanobis Distances }
\description{
  The MADIFA allows a factorial decomposition of the Mahalanobis
  distances.  This method is presented here in the framework of
  niche-environment studies.\cr
  \code{predict.madifa} allows the computation of the Mahalanobis
  Distances based on a restricted number of factorial axes.\cr
  All other functions allow various graphical displays of the results of
  the MADIFA.
}
\usage{
madifa(dudi, pr, scannf = TRUE, nf = 2)
\method{print}{madifa}(x, \dots)
\method{scatter}{madifa}(x, xax = 1, yax = 2, pts = FALSE, percent = 95,
               clabel = 1, side = c("top", "bottom", "none"),
               Adensity, Udensity, Aangle, Uangle, Aborder,
               Uborder, Acol, Ucol, Alty,
               Ulty, Abg, Ubg, Ainch, Uinch, \dots)
\method{hist}{madifa}(x, scores = TRUE, type = c("h", "l"), adjust = 1, Acol,
            Ucol, Aborder, Uborder, Alwd = 1, Ulwd = 1, \dots)
\method{predict}{madifa}(object, index, attr, nf, \dots)
s.madifa(x, xax = 1, yax = 2, cgrid = 1, clab = 1, \dots)
\method{plot}{madifa}(x, index, attr, xax = 1, yax = 2, cont = FALSE, \dots)
}
\arguments{
  \item{dudi}{a duality diagram, an object of class \code{dudi}}
  \item{pr}{a vector giving the utilization weights associated to each
    unit}
  \item{scannf}{logical.  Whether the eigenvalues barplot should be
    displayed}
  \item{nf}{an integer indicating the number of kept
    factorial axes}
  \item{x,object}{an object of class \code{madifa}}
  \item{xax}{the column number for the x-axis}
  \item{yax}{the column number for the y-axis}
  \item{pts}{logical. Whether the points should be drawn.  If
    \code{FALSE}, minimum convex polygons are displayed}
  \item{percent}{100 minus the proportion of outliers to be excluded
    from the computation of the minimum convex polygons}
  \item{clabel}{a character size for the columns}
  \item{side}{if \code{"top"}, the legend of the kept axis is upside, if
    \code{"bottom"} it is downside, if \code{"none"} no legend}
  \item{Adensity}{the density of shading lines, in lines per inch, for
    the available pixels polygon. See \code{\link{polygon}} for more
    details}
  \item{Udensity}{the density of shading lines, in lines per inch, for
    the used pixels polygon. See \code{\link{polygon}} for more details}
  \item{Aangle}{the slope of shading lines, given as an angle in degrees
    (counter-clockwise), for the available pixels polygon}
  \item{Uangle}{the slope of shading lines, given as an angle in degrees
    (counter-clockwise), for the used pixels polygon}
  \item{Aborder}{the color for drawing the border of the available pixels 
    polygon (or of the bars of the histogram) . See
    \code{\link{polygon}} for more details}
  \item{Uborder}{the color for drawing the border of the used pixels
    polygon (or of the bars of the histogram).  See
    \code{\link{polygon}} for more details}
  \item{Acol}{the color for filling the available pixels polygon. 
    if \code{pts == FALSE}, the color for the points corresponding to
    available pixels}
  \item{Ucol}{the color for filling the used pixels polygon. 
    if \code{pts == FALSE}, the color for the points corresponding to used
    pixels}
  \item{Alty}{the line type for the available pixels polygon, as in
    \code{par}}
  \item{Ulty}{the line type for the used pixels polygon, as in
    \code{par}}
  \item{Abg}{if \code{pts == TRUE}, background color for open plot symbols
    of available pixels}
  \item{Ubg}{if \code{pts == TRUE}, background color for open plot symbols
    of  used pixels}
  \item{Ainch}{if \code{pts == TRUE}, heigth in inches of the available
    pixels}
  \item{Uinch}{if \code{pts == TRUE}, heigth in inches of the largest used
    pixels}
    \item{scores}{logical.  If \code{TRUE}, the histograms display
    the row scores of the MADIFA.  If \code{FALSE}, they display the
    niche on the environmental variables (in this case, this is
    equivalent to \code{histniche})}
  \item{type}{what type of plot should be drawn.  Possible types are:\cr
    * \code{"h"} for histograms,\cr
    * \code{"l"} for kernel density estimates (see \code{?density}).\cr
    By default, \code{type = "h"} is used.  If \code{type = "l"} is
    used, the position of the mean of each distribution is indicated by
    dotted lines}
  \item{adjust}{if \code{type = "l"}, a parameter used to control the
    bandwidth of the density estimates (see \code{?density})}
  \item{Alwd}{if \code{type = "l"}, the line width of the kernel density 
    estimates of the available pixels}
  \item{Ulwd}{if \code{type = "l"}, the line width of the kernel density 
    estimates of the used pixels}
  \item{cgrid}{a character size, parameter used with par("cex")*
    \code{cgrid} to indicate the mesh of the grid}
  \item{clab}{if not NULL, a character size for the labels, used with
    \code{par("cex")*clab}}
  \item{index}{an integer vector giving the position of the rows of
    \code{tab} in the  initial object of class \code{kasc}}
  \item{attr}{an object of class \code{kasc} or \code{mapattr}}
  \item{cont}{logical.  Whether contour lines should be added to the maps}
  \item{\dots}{additional arguments to be passed to the functions
    \code{print}, \code{scatter}, and \code{plot}}
}
\details{
  The Mahalanobis distances are often used in the context of
  niche-environment studies (Clark et al. 1993, see the function 
  \code{mahasuhab}).  Each environmental variable defines a dimension in
  a multidimensionnal space, namely the ecological space.  The
  Mahalanobis distance between any resource unit in this space (e.g. the
  point defined by the values of environmental variables in a pixel of a
  raster map) and the centroid of the niche (the distribution of used
  resource units) can be used to give a value of eccentricity to this
  point.\cr
  
  For a given distribution of available resource units, for which a
  measure of Mahalanobis distances is desired, the MADIFA (MAhalanobis
  DIstances Factor Analysis) partitions the ecological space into a set
  of axes, so that the first axes maximises the average proportion of
  their squared Mahalanobis distances.  Note that the sum of the squared
  scores of any resource unit on all the axes of the analysis is equal
  to the squared Mahalanobis distances for this resource unit.  Thus,
  the MADIFA partitions the Mahalanobis distances into several axes of
  biological meaning (see examples).  \code{predict.madifa} allows to
  compute approximate Mahalanobis distances from the axes of the
  MADIFA.\cr

  \code{plot.madifa} returns a graphical summary of the analysis: it
  returns graphs of (i) the eigenvalues of the analysis (each eigenvalue
  measures the average Mahalanobis distance explained by each factorial
  axis); (ii) the scores of the habitat variables (i.e. the coefficients
  associated to each environmental variable in the linear combination
  defining the axes) - note that as the ecological space is distorted to
  "sphericize" the niche, the factorial axes are no longer orthogonals,
  and the scores of the variables are distributed within an ellipsoid
  instead of an hypersphere of radius equal to one in classical PCA. The
  limits of this ellipsoid is displayed on this graph, to see the amount
  of distortion done by the analysis (further research needs yet to be
  done on this graph); (iii) The projection of the available and used
  points on the factorial plane of the MADIFA; (iv) The map of the
  Mahalanobis distances computed from the original environmental
  variables; (v) the map of the approximated Mahalanobis distances
  computed from the two axes displayed in this plot; the correlations
  between the original environmental variables and the factorial axes;
  (v) the map of the first factorial axis and (vi) the map of the second
  factorial axis.\cr

  \code{hist.madifa} returns a graph of the niche and the available
  resource units on the factorial axes of the analysis.\cr
}
\value{
  \code{madifa} returns a list of class \code{madifa} containing the
  following components: 
  \item{call}{original call.}
  \item{tab}{a data frame with n rows and p columns.}
  \item{pr}{a vector of length n containing the number of points in each
    pixel of the map.}
  \item{nf}{the number of kept factorial axes.}
  \item{eig}{a vector with all the eigenvalues of the analysis.}
  \item{lw}{row weights, a vector with n components.}
  \item{li}{row coordinates, data frame with n rows and nf columns.}
  \item{l1}{row normed coordinates, data frame with n rows and nf columns.}
  \item{cw}{column weights, a vector with p components.}
  \item{co}{column coordinates, data frame with p rows and nf columns.}
  \item{mahasu}{a vector of length n containing the squared Mahalanobis
    distances for the n units.}
  \item{cor}{the correlation between the MADIFA axes and the original
    variable}
  \code{predict.madifa} returns a matrix of class \code{asc}. 
}
\references{
  Clark, J.D., Dunn, J.E. and Smith, K.G. (1993) A multivariate model of
  female black bear habitat use for a geographic information
  system. \emph{Journal of Wildlife Management}, \bold{57}, 519--526.\cr

  Calenge, C., Darmon, G., Basille, M., Loison, A. and Jullien
  J.M. (2008) The factorial decomposition of the Mahalanobis distances
  in habitat selection studies. \emph{Ecology}, \bold{89}, 555--566.\cr
}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\seealso{\code{\link{mahasuhab}} for a detailed description of the
  Mahalanobis distances, \code{\link{enfa}} and \code{\link{gnesfa}} for
  closely related methods.}
\examples{
\dontrun{

data(bauges)

map <- bauges$kasc
locs <- bauges$loc

## We prepare the data for the MADIFA
(datamad1 <- data2enfa(map, locs))

## We then perform the PCA before the MADIFA
pc <- dudi.pca(datamad1$tab, scannf = FALSE)
(mad <- madifa(pc, datamad1$pr, nf=7,
               scannf = FALSE))

#########################################
##                                     ##
## Graphical exploration of the MADIFA ##
##                                     ##
#########################################

hist(mad)

plot(mad, datamad1$index, datamad1$attr)

## this plot represents:
##  - the eigenvalues diagram
##  - the scores of the columns on the axes
##  - a graph of the niche in the available space
##  - a map of the Mahalanobis distances computed
##    using all environmental variables
##  - a map of the Mahalanobis distances computed
##    using the two factorial axes used in the
##    previous graphs
##  - the correlation between habitat variables
##    and factorial axes
##  - the geographical maps of the two
##    factorial axes

## predict with just the first axis
plot(sqrt(predict(mad, datamad1$index, datamad1$attr, 1)))



#########################################
##                                     ##
## Mathematical properties of  MADIFA  ##
##                                     ##
#########################################

## mad$li is equal to mad$l1, up to a constant (mad$l1 is normed)
plot(mad$li[,1],mad$l1[,1])

## This constant is the square root of the corresponding eigenvalue:
## the variance of mad$l1 is equal to the eigenvalue
apply(mad$l1,2,function(x) sum(x^2))/nrow(mad$li)

## the variance of mad$l1 weighted by pr is equal to 1
apply(mad$l1,2,function(x) sum(mad$pr*x^2)/sum(mad$pr))

## Therefore, the eigenvalues are equal to the average of Mahalanobis
## distance for the available resource units on each axis
mean(mahalanobis(matrix(mad$l1[,1], ncol=1), 0, 1))
mad$eig[1]

## Computation of the Mahalanobis distances
ma1 <- c(mahasuhab(map, locs))[datamad1$index]

## The sum of squared scores for a given Resource unit is equal to the
## Mahalanobis distances
ma2 <- apply(mad$l1,1, function(x) sum(x^2))
plot(ma1,ma2)

}

}
\keyword{multivariate}
