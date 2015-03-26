\name{eisera}
\alias{eisera}
\alias{print.esr}
\alias{scatter.esr}
\title{Eigenanalysis of Selection Ratios}
\description{
  Performs an eigenanalysis of selection ratios.
}
\usage{
eisera(used, available, scannf = TRUE, nf = 2)
\method{print}{esr}(x, \dots)
\method{scatter}{esr}(x, xax = 1, yax = 2, 
            csub = 1, possub = "bottomleft", \dots)
}
\arguments{
  \item{used}{a data frame containing the *number* of relocations of
    each animal (rows) in each habitat type (columns) }
  \item{available}{a data frame containing the *proportion* of
    availability of each habitat type (columns) to each animal (rows)}
  \item{scannf}{logical.  Whether the eigenvalues bar plot should be
    displayed}
  \item{nf}{if \code{scannf = FALSE}, an integer indicating the number
    of kept axes}
  \item{x}{an object of class \code{esr}}
  \item{xax}{the column number for the x-axis}
  \item{yax}{the column number for the y-axis}
  \item{csub}{a character size for the legend, used with 
              \code{par("cex")*csub}}
  \item{possub}{a string of characters indicating the sub-title 
                position ("topleft", "topright", "bottomleft", 
                "bottomright")}
  \item{\dots}{further arguments passed to or from other methods}
}
\details{
  The eigenanalysis of selection ratios has been developped to explore
  habitat selection by animals monitored using radio-tracking, when
  habitat is defined by several categories (e.g. several vegetation
  types, see Calenge and Dufour 2006).\cr

  This analysis can be used for both designs II (same availability for
  all animals, e.g. selection of the home range within the study area)
  and designs III (different availability, e.g. selection of the sites
  within the home range).  In the latter case, when some available
  proportions are equal to zero, the selection ratios are replaced by
  their expectation under random habitat use, following the
  recommendations of Calenge and Dufour (2006).
}
\value{
  A list of class \code{esr} and \code{dudi} containing also:
  \item{available}{available proportions}
  \item{used}{number of relocations}
  \item{wij}{selection ratios}
}
\references{
  Calenge, C. and Dufour, A.B. (2006) Eigenanalysis of selection ratios
  from animal radio-tracking data. \emph{Ecology}. \bold{87},
  2349--2355.
}
\author{ Clement Calenge \email{clement.calenge@oncfs.gouv.fr} }
\seealso{ \code{\link{wi}} for further information about the
  selection ratios, \code{\link{compana}} for compositional analysis.}
\examples{

###########################################################
###########################################################
###
###  Example given in Calenge and Dufour 2006 (design II)


data(squirrel)

## computation of the number of relocations in each habitat type
## from the data given by Aebischer et al. (1993).
## squirrel$locs give the percentage of relocations in each habitat
## type, and Aebischer et al. (1993) indicate that there are 30
## relocations per animal.
## We therefore compute the number of relocations in each habitat type
## using:
us <- round(30 * squirrel$locs / 100)

## Habitat availability 
av <- squirrel$studyarea

## Eigenanalysis of selection ratios
ii <- eisera(us, av, scannf = FALSE)

scatter(ii, grid = FALSE, clab = 0.7)

## The following graph may help the interpretation
## (see Calenge and Dufour 2006)
data(squirreloc)
locs <- squirreloc$locs
are <- squirreloc$map
co <- attr(are, "info")

li <- split(locs[,2:3], locs[,1])
opar <- par(mfrow=n2mfrow(length(li)), mar=c(0,0,2,0))
lapply(1:length(li), function(i) {
plot(are, colp = co[,2], main=names(li)[i], axes=FALSE)
points(li[[i]], pch=16, cex=1.5)
box()
})
plot(0,0, axes=FALSE, ty="n", xlim=c(-1,1), asp=1)
legend(-0.8,0.8, unique(co[,1]), fill=unique(co[,2]))
par(opar)


###########################################################
###########################################################
###
###  Example of design III

iii <- eisera(us, squirrel$mcp, scannf = FALSE)
scatter(iii, grid = FALSE, clab = 0.7)

}
\keyword{multivariate}
\keyword{hplot}
