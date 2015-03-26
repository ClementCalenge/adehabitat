\name{kselect}
\alias{kselect}
\alias{print.kselect}
\alias{hist.kselect}
\alias{kplot.kselect}
\alias{plot.kselect}
\title{ K-Select Analysis: a Method to Analyse the Habitat Selection by
  Animals  }
\description{
Performs a multivariate analysis of ecological data (K-select analysis).
}
\usage{
kselect(dudi, factor, weight, scannf = TRUE, nf = 2, ewa = FALSE)
\method{print}{kselect}(x, \dots)
\method{kplot}{kselect}(object, xax = 1, yax = 2, csub = 2, possub = c("topleft",
              "bottomleft", "bottomright", "topright"),
              addval = TRUE, cpoint = 1, csize = 1, clegend = 2, \dots)
\method{hist}{kselect}(x, xax = 1, mar=c(0.1,0.1,0.1,0.1),
             ncell=TRUE, csub=2,
             possub=c("bottomleft", "topleft",
                      "bottomright", "topright"),
             ncla=15, ...)
\method{plot}{kselect}(x, xax = 1, yax = 2, \dots)
}
\arguments{
  \item{dudi}{an object of class \code{dudi}}
  \item{factor}{a factor with the same length as \code{nrow(dudi$tab)}}
  \item{weight}{a numeric vector of integer values giving the weight
    associated to the rows of \code{dudi$tab}}
  \item{scannf}{logical. Whether the eigenvalues bar plot
    should be displayed}
  \item{nf}{if \code{scannf = FALSE}, an integer indicating the number
    of kept axes}
  \item{ewa}{logical. If \code{TRUE}, uniform weights are given
    to all animals in the analysis. If \code{FALSE}, animal weights are
    given by the proportion of relocations of each animal (i.e. an
    animal with 10 relocations has a weight 10 times lower than an
    animal with 100 relocations)}
  \item{x, object}{an object of class \code{kselect}}
  \item{xax}{the column number for the x-axis}
  \item{yax}{the column number for the y-axis}
  \item{addval}{logical.  If \code{TRUE}, the frequency of the
    relocations per animal is displayed (see examples)}
  \item{cpoint}{the size of the points (if 0, the points where no
    relocations are found are not displayed)}
  \item{mar}{the margin parameter (see \code{help(par)})}
  \item{ncell}{logical. If \code{TRUE}, the histogram shows the
    distribution of the cells of 
    the raster map where at least one relocation is found. If
    \code{FALSE}, the histogram shows the distribution of the
    relocations}
  \item{csub}{the character size for the legend, used with
    \code{par("cex")*csub}}
  \item{csize}{the size coefficient for the points}
  \item{clegend}{the character size for the legend used by
    \code{par("cex")*clegend}}
  \item{possub}{a character string indicating the sub-title position
    \code{("topleft", "topright", "bottomleft", "bottomright")}}
  \item{ncla}{the number of classes of the histograms}
  \item{\dots}{additional arguments to be passed to the generic function
  \code{histniche}, \code{print} or, in the case of \code{plot.kselect},
  \code{s.distri}}
}
\details{
  The K-select analysis is intended for hindcasting studies of habitat
  selection by animals using radio-tracking data.  Each habitat variable
  defines one dimension in the ecological space.  For each animal, the
  difference between the vector of average available habitat conditions
  and the vector of average used conditions defines the marginality
  vector. Its size is proportional to the importance of habitat
  selection, and its direction indicates which variables are
  selected. By performing a non-centered principal component analysis of
  the table containing the coordinates of the marginality vectors of
  each animal (row) on the habitat variables (column), the K-select
  analysis returns a linear combination of habitat variables for which
  the average marginality is greatest. It is a synthesis of variables
  which contributes the most to the habitat selection. As with principal
  component analysis, the biological significance of the factorial axes is
  deduced from the loading of variables.\cr
  
  \code{plot.kselect} returns a summary of the analysis: it displays (i)
  a graph of the correlations between the principal axes of the PCA of
  the objects of class \code{dudi} passed as argument and the factorial
  axes of the K-select analysis; (ii) a graph giving the scores of the
  habitat variables on the factorial axes of the K-select analysis;
  (iii) the barplot of the eigenvalues of the analysis (each eigenvalue
  measure the mean marginality explained by the axis; (iv) the
  projection of the non-recentred marginality vectors on the factorial
  plane (the origin of the arrow indicates the average available habitat
  conditions, and the end of the arrow indicates the average used
  conditions); (v) the projection of the resource units available to
  each animal on the first factorial plane and (vi) the coordinates of
  the recentred marginality vectors (i.e. recentred so that they have a
  common origin) on the first factorial plane.\cr

  \code{kplot.kselect} returns one graph per animal showing the
  projections of the available resource units on the factorial plane, as
  well as their use by the animal.  \code{hist.kselect} does the same
  thing, but on one dimension instead of two.
}
\value{
  \code{kselect} returns a list of the class \code{kselect} and
\code{dudi} (see \code{\link[ade4]{dudi}}). 
}
\references{
  Calenge, C., Dufour, A.B. and Maillard, D. (2005) K-select analysis: a
  new method to analyse habitat selection in radio-tracking
  studies. \emph{Ecological modelling}, \bold{186}, 143--153.
}
\author{Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }

\seealso{\code{\link{sahrlocs2kselect}} for
  conversion of objects class \code{sahrlocs} to objects suitable for a
  K-select analysis, \code{\link[ade4]{s.distri}}, and
  \code{\link[ade4]{dudi}} for class \code{dudi}.}

\examples{

\dontrun{
## Loads the data
data(puechabon)
sahr <- puechabon$sahr

## prepares the data for the kselect analysis
x <- sahrlocs2kselect(sahr)
tab <- x$tab

## Example of analysis with two variables: the slope and the elevation.
## Have a look at the use and availability of the two variables
## for the 4 animals
tab <- tab[,((names(tab) == "Slope")|(names(tab) == "Elevation"))]
tab <- scale(tab)
tmp <- split.data.frame(tab, x$factor)
wg <- split(x$weight, x$factor)
opar <- par(mfrow = n2mfrow(nlevels(x$factor)))
for (i in names(tmp))
  s.distri(scale(tmp[[i]]), wg[[i]])
par(opar)

## We call a new graphic window
x11()
## A K-select analysis
acp <- dudi.pca(tab, scannf = FALSE, nf = 2)
kn <- kselect(acp, x$factor, x$weight,
 scannf = FALSE, nf = 2)

# use of the generic function scatter
scatter(kn)

# Displays the first factorial plane
kplot(kn)
kplot(kn, cellipse = 0, cpoint = 0)
kplot(kn, addval = FALSE, cstar = 0)

# this factorial plane can be compared with
# the other graph to see the rotation proposed by
# the analysis
graphics.off()

# Displays the first factorial axis
hist(kn)

# Displays the second factorial axis
hist(kn, xax = 2)

# Summary of the analysis
plot(kn)
}
}
\keyword{multivariate}
\keyword{hplot}
