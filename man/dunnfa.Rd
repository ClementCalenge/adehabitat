\name{dunnfa}
\alias{dunnfa}
\alias{print.dunnfa}
\title{Factorial Analysis of the Specialization in Habitat Selection
  Studies.  Unpublished Work of James Dunn (University of Arkansas)}
\description{
  \code{dunnfa} performs a factorial decomposition of the Mahalanobis
  distances in habitat selection studies (see details). 
}
\usage{
dunnfa(dudi, pr, scannf = TRUE, nf = 2)
\method{print}{dunnfa}(x, \dots)
}
\arguments{
  \item{dudi}{an object of class \code{pca}}
  \item{pr}{a vector giving the utilization weights associated to
    each unit}
  \item{scannf}{logical.  Whether the eigenvalues barplot should be
    displayed}
  \item{nf}{an integer indicating the number of kept axes}
  \item{x}{an object of class \code{dunnfa}}
  \item{\dots}{additional arguments to be passed to the function
    \code{print}}
}
\details{
  This analysis is in essence very similar to the MADIFA (see
  \code{?madifa}).   The Mahalanobis distances are often used in the
  context of niche-environment studies (Clark et al. 1993, see the
  function  \code{mahasuhab}).  Each resource unit takes a value on a
  set of environmental variables.  Each environmental variable defines a
  dimension in a multidimensionnal space, namely the ecological space.
  A set of points (resource units) describes what is available to the
  species.  For each point, a "utilization weight" measures the
  intensity of use of the point by the species.  The set of points for
  which the utilization weight is greater than zero defines the
  "niche". The Mahalanobis distance between any resource unit in this
  space (e.g. the point defined by the values of environmental variables
  in a pixel of a raster map) and the centroid of the niche (the
  distribution of used resource units) can be used to give a value of
  eccentricity to this point.\cr
  
  For a given distribution of available resource units, for which a
  measure of Mahalanobis distances is desired, the MADIFA (MAhalanobis
  DIstances Factor Analysis) partitions the ecological space into a set
  of axes, so that the first axes maximises the average proportion of
  their squared Mahalanobis distances.  James Dunn (formerly University
  of Arkansas) proposed the analysis programmed in the function
  \code{dunnfa}, as an alternative to the MADIFA (unpublished
  results).  This analysis is closely related to both the ENFA
  (Ecological niche factor analysis, Hirzel et al. 2002) and the
  MADIFA.\cr

  The analysis proposed by James Dunn searches, in the multidimensional
  space defined by environmental variables, synthesis variables which
  maximise the ratio (variance of the scores of available resource
  units) / (variance of the scores of used resource units).  This ratio
  is sometimes called "specialization" in the ecological literature
  (Hirzel et al. 2002).  It is therefore very similar to the ENFA (which
  also maximises the specialization), except that the factorial axes
  returned by this analysis are not required to be *orthogonal to the
  marginality axis*.\cr

  James Dunn demonstrated that this analysis also partitions the
  Mahalanobis distances into uncorrelated axes, which makes it similar
  to the MADIFA (the difference is that the MADIFA maximises the mean
  squared Mahalanobis distances on the first axes, whereas the DUNNFA
  maximises the specialization on the first axes).  Therefore, as for
  the MADIFA, the DUNNFA can be used to build reduced rank habitat
  suitability map.\cr

  Note that although this analysis could theoretically be used with all
  kinds of variables, it it currently implemented only for numeric
  variables.
}
\value{
  \code{dunnfa} returns a list of class \code{dunnfa} containing the
  following components: 
  \item{call}{original call.}
  \item{tab}{a data frame with n rows and p columns (original data frame
    centered by column for the uniform weighting).}
  \item{pr}{a vector of length n containing the number of points in each
    pixel of the map.}
  \item{nf}{the number of kept axes.}
  \item{eig}{a vector with all the eigenvalues of the analysis.}
  \item{liA}{row coordinates (centering on the centroid of the cloud of
  available points), data frame with n rows and nf columns.}
  \item{liU}{row coordinates (centering on the centroid of the cloud of
  available points), data frame with p rows and nf columns.}
  \item{mahasu}{a vector of length n containing the reduced-rank squared
  Mahalanobis distances for the n units.}
\item{co}{column (environmental variables) coordinates, a data frame
  with p rows and nf columns}
\item{cor}{the correlation between the DUNNFA axes and the original
  variable}
}
\references{
  Clark, J.D., Dunn, J.E. and Smith, K.G. (1993) A multivariate model of
  female black bear habitat use for a geographic information
  system. \emph{Journal of Wildlife Management}, \bold{57}, 519--526.\cr

  Hirzel, A.H., Hausser, J., Chessel, D. & Perrin, N. (2002)
  Ecological-niche factor analysis: How to compute habitat-suitability
  maps without absence data? \emph{Ecology}, \bold{83}, 2027--2036.  
  
  Calenge, C., Darmon, G., Basille, M., Loison, A. and Jullien
  J.M. (2008) The factorial decomposition of the Mahalanobis distances
  in habitat selection studies. \emph{Ecology}, \bold{89}, 555--566.\cr

}
\author{
  Clement Calenge \email{clement.calenge@oncfs.gouv.fr}
}
\note{
  This analysis was developed by James Dunn during an e-mail discussion
  on the MADIFA, and is still unpublished work.  Implemented in
  adehabitat with his autorization.
}
\seealso{ \code{\link{madifa}}, \code{\link{enfa}} and
  \code{\link{gnesfa}} for related methods.  \code{\link{mahasuhab}} for
  details about the Mahalanobis distances.}
\examples{

\dontrun{
data(bauges)

map <- bauges$kasc
locs <- bauges$loc

## We prepare the data for the analysis
(datadun1 <- data2enfa(map, locs))

## We then perform the PCA before the analysis
pc <- dudi.pca(datadun1$tab, scannf = FALSE)
(dun <- dunnfa(pc, datadun1$pr, nf=2,
               scannf = FALSE))

## We should keep only one axis:
barplot(dun$eig)


## The correlation of the variables with the first two axes:
s.arrow(dun$cor)

## A factorial map of the niche (centering on the available points)
scatterniche(dun$liA, dun$pr, pts=TRUE)

## a map of the reduced rank Maalanobis distances
## (here, with one axis)
dun2 <- dunnfa(pc, datadun1$pr, nf=1,
              scannf = FALSE)
kas <- getkasc(df2kasc(data.frame(dun2$mahasu,dun2$mahasu),
               datadun1$index, map), 1)
image(kas)


## Compute the specialization on the row scores of
## the analysis:
apply(dun$liA, 2, function(x) {
   varav <- sum((x - mean(x))^2) / length(x)
   meanus <- sum(dun$pr*x)/sum(dun$pr)
   varus <- sum(dun$pr * (x - meanus )^2)/sum(dun$pr)
   return(varav/varus)
})
## The eigenvalues:
dun$eig

}


}
\keyword{multivariate}
